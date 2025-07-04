#---------------------------------------------------------
#
#		 perl modules for some geographic functions and extraction
# 	 of BERNESE CRD and COV files
#
#			&gdatum(): 													extracting datum information from BERNESE DATUM. definitioin file
#																						hardcoded for IGS05 at the moment
#     &rdcov(FIL): 												Reads BERNESE covariance file and returnes covariance matrix
#     &xyzell(datum,xyzcoor): 						XYZ -> llh using the datum designed for using 
#     																			&datum() as input
#     &ellecc(llh,vector,optional)
#			&eccell(llh,vector, swich ):			converting a vector in xyz coord. system to local enu system in m
#     &xyz2enumatr(refarray,swich):				rotation matrix to convert xyz coord. system to local enu system or vice versa
#			&redmatr(in, refind,covmat):
#		  &greatcircd(coor1,coor2):						distance along great circle
#			&M:
#			&rotpole(plate,xyz,swich):
#			&fitPeriod(xdata,ydata,weight)
#
#
#    10.11.07 BGO.
#
#
#
#
#---------------------------------------------------------

use strict;
use PDL;
use PDL::Core ':Internal';
#use PDL::LinearAlgebra;
#use Math::Cephes qw(atan2);
use PDL::NiceSlice;
use PDL::Fit::Linfit;
use Geo::Constants qw(PI);
use PDL::Options;
###------------------------------------------------------------------------------------
#special functions for extracting data from bernese files 

#---------------------------------------------------------
# 
#    "&gatum()"
# Extracts Datum information from bernese DATUM. file
# returns major axes ($ae), minor axes ($be), scale ($sc), shift (@dx), rotation (@rx)  
#
#
#
#---------------------------------------------------------
sub gdatum {
	 #defining variables
	 my @datum;
	 my (@dx,@rx) = [0,0,0];
	 my ($ae, $invf, $sc, $be) = undef;

	 my $DATUM="IGS05"; #defining the datum
	 #my $DATUM="PZ - 90";
	 open my $FIL, 'DATUM.' or die "Could not open file DATUM.: $!"; #the bernese file for Datum def.
	 while (<$FIL>) {
			if (/^$DATUM/) { #Found the DATUM
				 chomp( @datum[0,1,2] = ($_,<$FIL>));
         ($ae, $dx[0], $rx[0], $invf, $dx[1], $rx[1], $sc, $dx[2],$rx[2]) 
				                = join (" ", @datum) =~ m/(-?\d+\.\d+D\+?\-?\d+|-?\d+\.\d+)/g;
				 $sc =~ s/D/e/;
				 $be=$ae*(1-1/$invf); #minor axis of ellipsoid ($a is major axis)
				 $sc=1+$sc; #scale of ellipsoid
				 @rx = map $_/206264.806, @rx;
				 #print $ae," ", $be," ", $invf, " ", $sc, " ","@dx"," ","@rx", "\n";
				 last;
			}
	 }
	 close $FIL;
	 die "Datum ($DATUM) not found" unless $ae;
	 #returning the values
	 return [$ae, $be, $sc, \@dx, \@rx ];
}

#--------------------------------------------------------------
#
# The module reads a covariance file in bernese format (expects and 
# open file) and prints out the covariance matrics in upper triangular
# form.
# 			&rdcov(FIL) name of covariance file
# 			VERY PORLY WRITTEN NEEDS TO BE REVISED
#
#-------------------------------------------------------------

sub rdcov {
	 my $cov = pdl;
	 my (@tmp,@cova);
	 my @name1;
	 my %nameind;
	 my $line;
	 my ($head,$rmsunit, $obs, $unknowns) = undef; 
	 my $COV =  ${ shift() };
	 # The header-----------------------------------------------
	 for (my $i =1;$i<=6;$i++) {
			$head = <$COV>;
	 }
	 $head = <$COV>;
	 # extracting parameters
	 ($rmsunit, $obs, $unknowns) = $head =~ m/(\d+\.\d+|\d+)/g;
##	 print  $rmsunit, " ", $obs, " " ,$unknowns,"\n";
	 for (my $i =1;$i<=4;$i++) {
			$head = <$COV>;
	 }
	 #-----------------------------------------------------------
	 my $i=0;
	 while (<$COV>) { #looping over the matrix
			 unless (/^\s+$/) { #If the line has any content
					(@tmp) = $_ =~m/^(\w+\s+(?:\w+)?)\s+(\w)\s+(\w+\s+(?:\w+)?)\s+(\w)\s+(-?\d\.\d+D(?:\-?|\+?)\d+)/;
					#after every empty line a new colum in the matrix starts
					if ($i == 0) { 
						 unless (@name1) {
								push @name1 , $tmp[0];
								$nameind{substr($name1[0],0,4)} = 0;
						 }	else  {
								unless ($name1[$#name1] eq $tmp[0] ) {
									 push @name1, $tmp[0];
								   $nameind{substr ($name1[$#name1],0,4)} = $#name1;
								}
					   }
					}
					$tmp[4] =~s/D/e/;
					push @cova, $tmp[4]; 
					$i++
			 } else { #empty lines
					$i = 0;
			 }

	 } 
   my $dimen = (sqrt(1+8*($#cova+1))-1)/2;
	 $cov = zeroes(($dimen),($dimen));
	 my $k = 0;
	 for(my $i=0;$i<$dimen;$i++) {
			for (my $j=0;$j<=$i;$j++) {
				 $cov(($i),($j)) .= double($cova[$k]);
				 $cov(($j),($i)) .= double($cova[$k]) if ($i!=$j);
				 $k++;
			}
	 }
   $cov = $cov*$rmsunit**2; 
	 return ( $rmsunit,$obs,$unknowns, \%nameind, \@name1, $cov );
}

##
#####----------------------------------------------------------------------------------
#



#----------------------------------------------------------------------------------
#
#
#  x,y,z => lat,long,height
#  usage: &xyzell(datum,xyzcoor)
#  			 datum: reference to datum definition 
#  			 				expects &datum
#  			 xyzcoor: reference to a xyz coordinates 
#  			 returns reference to a list 
#  			 (LATITUDE, LONGITUDE, HEIGHT(over ellipsoid)) 
#
#
#
#----------------------------------------------------------------------------------

sub xyzell {
	 my ($dxell, $drell); # (Shift of the ellipsoid, Rotation of ellipsoid)
	 my ($aell,$bell,$scell); # Ellipsoid (major axis,minor axis, scale)
	 my ($RLAM, $PHI, $H); # return elements (Longditude, lattitude, heigth over ellipsoid )

	 ($aell, $bell, $scell, $dxell, $drell) = @{ shift() } ; #reading in the datum definition
	 my $xp = double (shift()); #cartesian coordinates
	 my $piddle = shift();
	 $xp = (pdl($xp) - pdl($dxell))/$scell; #Shifting and scaling
	 my ($sa, $ca, $sb, $cb, $sc, $cc ) = map {sin($_),cos($_)} @$drell; #preparing the rotation
	 
	 #The rotation matrix
	 my $rot = pdl( [ 
			[	 $cb*$cc, $ca*$sc+$sa*$sb*$cc, $sa*$sc-$ca*$sb*$cc ],
			[ -$cb*$sc, $ca*$cc-$sa*$sb*$sc, $sa*$cc+$ca*$sb*$sc ],
			[  		 $sb, 					 -$sa*$cb,             $ca*$cb ],					
	 ] );
	 $xp =$xp x $rot; #rotating

	 #parameters for the calculations
	 my $E2 = ($aell*$aell-$bell*$bell)/($aell*$aell);
	 my $S= sqrt inner($xp(0:1,0),$xp(0:1,0));
	 
	 #
	 #CHECK IF  X AND Y ARE BOTH ZERO, COMPUTE LONGITUDE
	 #
	 if ($S != 0.0) {
	 		$RLAM = atan2( $xp(1),$xp(0) );
	 } else {
	 		$RLAM = 0.0;
	 }

	 #
	 #COMPUTE APPROXIMATE HEIGHT AND LATITUDE
	 #
	 $H=sqrt(inner($xp,$xp))-$aell;
	 if ($H < -1e6) {
			warn "The antenna is more than 1000 Km below the earths surface\n",
			     "The height is ", $H->at(0), " m.";
	 } 
	 if ($H > 20e6) {
			warn "A priori coordinates higher tha 20000 Km above the earths surfac\n",
					 "The height is $H m";
	 }
	 if ($aell+$H != 0.0) {
			$PHI = atan2($xp(2), $S*(1-$E2*$aell/($aell+$H)));
	 } else {
			$PHI = 0.0;
	 }

	 #
	 # NORTH OR SOUTH POLE?
	 #
	 if ($S == 0.0) {
			$H = abs($xp(2))-$bell;
	 }

	 #
	 # ITERATE TO COMPUTE HEIGHT AND LATITUDE
	 # 
	 my ($PHIP,$HP) = null;
	 do { 
			my $H0 = $H;
	 		my $PHI0 = $PHI;
	 		my $N=$aell/sqrt(1-$E2*sin($PHI)**2);
   		$HP=$H;
   		$PHIP=$PHI;
   		$H=$S/cos($PHI)-$N;
   		if ($N+$H != 0) {
				$PHI = atan2($xp(2),$S*(1-$E2*$N/($N+$H)));
	 		} else {
				$PHI = 0;
	 		}
	 } while (abs($PHIP-$PHI) > 1e-11 || abs($HP-$H) > 1e-5);
	 
 	 if ($piddle){
			return pdl($PHI->at(0,0),$RLAM->at(0,0),$H->at(0,0));
	 } else {
	 		return [ ($PHI->at(0,0),$RLAM->at(0,0),$H->at(0,0)) ];
	 }

}

#-------------------------------------------------------------------
#

#-------------------------------------------------------------------
#
#	neu displacement -> xyz displacement
#
#	usage: &ellecc(llh,vector,optional)	
#					llh: reference to coordinates in ll or llh
# 			 vector: reference to vector array in neu coor
# 			 optional: returns a piddle if optional returns 
# 			 true on logiacal test
#
# 			 returns a re	ference to a list of displacements
# 			 vectors in xyz 
# 			 
# 			 or a vector piddle xyz in case of optional ret. true
#
#
#-------------------------------------------------------------------

sub ellecc {
	 my $xstell =  shift() ;
	 my $locecc = pdl(  map {[ $_ ]} @{shift()}  );
	 my $piddle = shift();
	 #calculating transformation matrix
	 my $rmat = &xyz2enumatr( [@{$xstell}[0,1]],1 ); 
my $xstecc =  transpose($rmat) x $locecc;
   if ($piddle){
			return $xstecc;
	 } else {
			return  [ list $xstecc ];
	 }
}

#-------------------------------------------------------------------
#
# xyz displacement -> neu displacement
#
# usage: &eccell(llh,vector, optional )
# 			 llh: reference to coordinates in ll or llh
# 			 vector: reference to vector array in xyz coor
# 			 optional: returns a piddle if optional returns 
# 			 true on logiacal test
#
# 			 returns a re	ference to a list of displacements
# 			 vectors in neu
# 			 
# 			 or a vector piddle neu in case of optional ret. true
#
#
#-------------------------------------------------------------------

sub eccell {
	 my $xstell =  shift() ;
	 my $xstecc = pdl(  map {[ $_ ]} @{shift()}  );
	 my $piddle = shift();
	 #calculating transformation matrix
	 my $rmat = &xyz2enumatr( [@{$xstell}[0,1]],1 ); 
	 my $locecc =  $rmat x $xstecc;
   if ($piddle){
			return $locecc;
	 } else {
			return  [ list $locecc ];
	 }
}

#------------------------------------------------------
# calculate tranformation matrix to convert  xyz to 
# local coordinate system enu or vice versa
# usage: &xyz2enumatr(refarray,swich) 
# 				refarray: refference to (latitude,longitude)
# 				swich: desites format of the return value
# 							if it returns true a 3x3 piddle is returned
# 							(needs PDL)
# 							if it returns false a reference to a 9 element
# 							list is returned
#
# 
#
#
#------------------------------------------------------

sub xyz2enumatr {
   my $phi = shift();
	 my $piddle = shift();
	 my ($sphi, $cphi, $slmb, $clmb) = map {sin, cos} @{$phi}[0,1];
 	 my $rmat = pdl 
			 				[ -$sphi*$clmb, -$sphi*$slmb, $cphi ],
			 				[       -$slmb,        $clmb,     0 ],
			 				[  $cphi*$clmb,  $cphi*$slmb, $sphi ],
			 ]);
   if ($piddle){
			return $rmat;
	 } else {
			return  [ list $rmat ];
	 }
}

#----------------------------------------------------------------
#    Redusing block matrix of arbitrary size to to 
#
#
#
#
#
#
#----------------------------------------------------------------

sub redmatr {
	 my $stat = pdl @{shift()};
	 my $cov = shift();
	 my $B = pdl;
	 unless ($stat(0) eq $stat(1)) {
      $B = zeroes(6,6);
			$B(0:2,0:2) .= $cov( (3*$stat(0)):(3*$stat(0))+2,(3*$stat(0)):(3*$stat(0))+2);
   		$B(3:5,0:2) .= $cov( (3*$stat(1)):(3*$stat(1))+2,(3*$stat(0)):(3*$stat(0))+2);
	 		$B(0:2,3:5) .= $cov( (3*$stat(0)):(3*$stat(0))+2,(3*$stat(1)):(3*$stat(1))+2);
	 		$B(3:5,3:5) .= $cov( (3*$stat(1)):(3*$stat(1))+2,(3*$stat(1)):(3*$stat(1))+2);
			return &M x $B x transpose(&M);
	 } else {
			$B = zeroes(3,3);
			$B(0:2,0:2) .= $cov( (3*$stat(0)):(3*$stat(0))+2,(3*$stat(0)):(3*$stat(0))+2);
			return $B;
	 }
}

#-----------------------------------------------------------------
#
#		Calculate the  distance along great circle between two points
#		on Earths surface
#		usage:  &greatcircd(coor1,coor2)
#								where coor1-2 are references to two points on earths
#								surface given in long,lat radiance
#
#
#
#
#------------------------------------------------------------------

sub greatcircd {
	 my $PI = 3.141592653589793;
	 my $earthr=6371000;
#	 my $trad = pdl $PI/180;

	 my ($long1, $lat1) = map {$_} @{shift()};
	 my ($long2,$lat2) = map {$_} @{shift()};
	 my ($Dlong,$a,$b);

	 $Dlong = $long2 - $long1;
	 $a = sqrt( ( cos($lat2)*sin($Dlong) )**2 + ( cos($lat1)*sin($lat2) - sin($lat1)*cos($lat2)*cos($Dlong) )**2 );
	 $b = sin($lat1)*sin($lat2) + cos($lat1)*cos($lat2)*cos($Dlong);
	 return   $earthr*atan2($a,$b) ;

}

#-----------------------------------------------------------
#
#
#
#
#
#
#----------------------------------------------------------




#------------------------------------
#
# Matrix of convieniense
#
#
#
#-----------------------------------

sub M {

	 my $M = pdl ([
					[ 1, 0, 0, -1, 0, 0 ],
					[ 0, 1, 0, 0, -1, 0 ],
					[ 0, 0, 1, 0, 0, -1 ],
			]);
}


#-------------------------------------------------------
#
# Calculating the nnr-nuvela velocity of a given coordinates in xyz
# usage: &rotpole($plate,$xyz,$wxwywz,option)
# 								$plate -> four word short name for the given plate
# 								$xyz -> coordinates in meters
# 								option -> true for piddle
# 													false for refr to list
#
#
#
#------------------------------------------------------

sub rotpole {
   my $plate = shift();
	 my $xyz =  double (shift());
	 my $wxwywz = double (shift()) unless $plate;
	 my $piddle = shift();

	 my $index;
	 my $nuvela = double 0.9562; #constant for changing NNR-NUVEL1 -> NNR-NUVEL1A
	 my $NNRNUVEL;
	 my $vel;

if ($plate) {
	 my @plates = (
			"PCFC",
			"AFRC",
			"ANTA",
			"ARAB",
			"AUST",
			"CARB",
			"COCO",
			"EURA",
			"NAZC",
			"NOAM",
			"SOAM",
			"JUFU",
			"INDI",
			"PHIL"
			);


	 #     DFA: NNR-NUVEL1
	 #     OMX(I),     OMY(I),    OMZ(I) #DEGREE/MY
	 $NNRNUVEL = double ([
	 		[ -0.0907,     0.2902,   -0.5976 ],
	 		[  0.0532,    -0.1856,    0.2348 ],
	 		[ -0.0494,    -0.1018,    0.2218 ],
	 		[  0.4003,    -0.0311,    0.4049 ],
	 		[  0.4695,     0.3072,    0.3762 ],
	 		[ -0.0109,    -0.2027,    0.0945 ],
	 		[ -0.6249,    -1.2944,    0.6544 ],
	 		[ -0.0590,    -0.1434,    0.1887 ],
	 		[ -0.0921,    -0.5138,    0.5756 ],
	 		[  0.0152,    -0.2155,   -0.0094 ],
	 		[ -0.0624,    -0.0906,   -0.0523 ],
	 		[  0.2995,     0.4805,   -0.2936 ],
	 		[  0.3995,     0.0026,    0.4066 ],
	 		[  0.5913,    -0.4412,   -0.5976 ],
	 ])*1.7453292e-08; # changing to #RAD/YR
	
		($index) =  grep (  $plate eq $plates[$_] ,  0..$#plates) or die "can't find plate $plate: $!";
  	$vel =  crossp ($NNRNUVEL(:,$index)*$nuvela,$xyz);
  } else {
		 #$NNRNUVEL = double ([
		 #	$wxwyw,
		 #]);
  	$vel =  crossp ($wxwywz,$xyz);
	}


	 if ($piddle){
			return $vel
	 } else {
			return [ list $vel ]
	 }
}


#--------------------------------------------------------------
#
#
#-------------------------------------------------------------

sub fitPeriod {
	  my $opthash = ref($_[-1]) eq "HASH" ? pop(@_) : {} ; 
    my %opt = parse( { Weights=>ones(1) }, $opthash ) ;
		my ($x, $y, $fpar) = @_;

		my $wt = $opt{Weights};
		my $ser_leng=$x->dim(0);
		my $fitFuncs;
    die "you need to specify a function:" unless (defined $fpar);
		if ($fpar == 4) {
			 $fitFuncs = cat ones($ser_leng),$x, cos($x*2*PI), sin($x*2*PI);
		} elsif ($fpar == 2) {
			 $fitFuncs = cat ones($ser_leng),$x;
		} else {
			die "there is no function defined"; 
		}
		return  linfit1d  $x, $y, $fitFuncs , {Weights=>$wt};
 }



#Returning true
1;
