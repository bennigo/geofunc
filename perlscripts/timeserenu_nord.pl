#!/usr/bin/perl -w
#
#
#
#
#-------------------------------------------------------
#
use strict;
use PDL;
use PDL::NiceSlice;
#use PDL::Core ':Internal';
use PDL::Slatec qw( matinv polyfit polfit );
use PDL::LinearAlgebra qw(mnorm);
use PDL::MatrixOps qw( inv );
use PDL::Fit::Polynomial;
use File::Basename qw(dirname basename fileparse);
use Date::Pcalc qw(Days_in_Year);
use Geo::Functions qw(rad_deg deg_rad);
use Geo::Constants qw(PI);

require "geo.pm";

#defining an initialising parameters
my $scale = '18/8';
#my $outdir = "kara/psfiles/itrf05/";
#my $outdir = "nord/psfiles/";
my $outdir = "/home/bgo/verkefni/gps_timar/disp/nord/psfiles/";

my ($paxis,$saxis, $area);
my ($year,$month,$day, $doy,$long,$lat,$height,$station,$N,$sigmamed,$ressq,$oversigma,$sc,$unctemp,$ind);
my ($filbase, $path,$suff,$psfileb,$psend,$psfile ); 
my ($start, $stop,$pdata,$coeff, $coeffuns, $label,$title,$nuvvel);
my (@start,@stop,@minmax,@minmaxt,@tmp,@int,@coeff,@coeffuns,@VDVA,@VDVB,@data,@xyz,@MEAN,@RMS);
my ($mean,$prms,$median,$min,$max,$adev,$rms,$days);

# scaling to mm or whatever
my $scaling = 1000;
my @comp = qw/East North Up/;
my @xlabel = ("","","Year");
my @pos = ("0.8 3.6", "5.6 0.6" );
my @sc = (5,4,3);


foreach  my $ARG (@ARGV) {
	 print $ARG,"\n";

	($filbase, $path,$suff ) = &fileparse($ARG,qw{enu.disp});
  $title = $filbase; #Station name
	#$psfileb = "test"; # psfile prefix 
	#$psend = "\.ps";
	#$psfile = ${psfileb}.${psend}; # psfile for gmt output
	$psfile = "${outdir}${filbase}enu\.ps"; # psfile for gmt output ####CHECK enu  and enu.disp make more general########3
	
	
	#changing variables through loop @tmp @minmax
 	
	#initiating plot
	system ("echo '0 0' | pstext  -JX${scale} -R-1/1.5/-1/3 -Y29.1 -P -K << EOF > ${psfile}\n0.3 -0.8 16 0 4 TC ${title}\nEOF\n");
	
	
	#reading in the data
	@data = rcols $ARG;
	
	#getting station information
	open FIL, $ARG;
	while (<FIL>) {
		 if (/^#{2}\s/) {
				chomp ((undef,undef,$xyz[0],$xyz[1],$xyz[2]) = split );
		 }
		 if (/^#{3}/) {
		 		chomp ((undef, $station, $long, $lat, $height ) = split);
				last;
		 }
		 
	}
	close FIL;

	
	my $interv = 0.01918;
	$start = &min($data[0]);
	$stop = &max($data[0]);
	$int[0] = which( $data[0] < 2007.0);
	$int[1] = setops( which( 2008.0 < $data[0]), 'AND', which( 2010.0 > $data[0]) ); 
	$int[2] = which( $data[0] > 2009.0);
	print $#int,"\n";
	#Fixing the Eurasian plate----------------------
	#Altamimi et al 2007 Absolute rotation pole fo Eurasia in ITRF2005 (56.330 deg,-95.979 deg,0.261 deg/my)
	#Altamimi et al 2007 Absolute rotation pole fo N. Amerika in ITRF2005 (-4.291 deg,-87.385 deg,0.192 deg/my)
	my @eurpole = map &rad_deg($_),(56.330,-95.979,0.261/1e6);
#	my @nAmerpole = map &rad_deg($_),(-4.291,-87.385,0.192/1e6);
	@eurpole = @{&ellecc([$eurpole[0],$eurpole[1]],[0,0,$eurpole[2]])};
#	@nAmerpole = @{&ellecc([$nAmerpole[0],$nAmerpole[1]],[0,0,$nAmerpole[2]])};
#	print join(" ",@eurpole),"\n";
#	print matinv(&xyz2enumatr(\@eurpole,1)) x ;
	$nuvvel = &eccell(&xyzell(&gdatum,\@xyz),&rotpole("",\@xyz,\@eurpole),1);
	print $nuvvel*1000;
#	my $namvel = &eccell(&xyzell(&gdatum,\@xyz),&rotpole("",\@xyz,\@nAmerpole),1);
#	print "-----------------\n";
	#$nuvvel = &eccell(&xyzell(&gdatum,\@xyz),&rotpole("EURA",\@xyz,),1);
	# my $vel=(pdl($namvel-$nuvvel)*$scaling)->squeeze ;
	#print $vel,"\n"; 
	#print mnorm($vel),"\n";
	#my $angle =  -(atan2($vel(0),$vel(1)))->sclr + &PI/2;
	#print  360+deg_rad($angle),"\n";
	#print "-----------------\n";
	$data[1] .= $data[1] - $nuvvel(0,1)*($data[0]-$start);
	$data[3] .= $data[3] - $nuvvel(0,0)*($data[0]-$start);
	#-------------------------------------------------------------------	

	for (my $i=0;$i<=2;$i++) {

		 # subtracting the average of first few days"
    If plate is passed
        Calculate the nnr-nuvela velocity of a given plate in geocentric coordinates.
    If plate is not passed
        Calculate the velocity in a given point in geocentric coordinates xyz with for a given rotation pole wxwywz
    """
    if plate != None:
        pass
    else:
        if xyz != None and wxwywz != None:
 
		 $days=0.002739726*7;
		 $ind = which($data[0] < $data[0]->at(0)+$days);
		 ($mean,$prms,$median,$min,$max,$adev,$rms) = 
		 stats($data[2*$i+1]->index($ind),1/($data[2*$i+2]->index($ind))**2);
		 $data[2*$i+1]=$data[2*$i+1]-$mean;

		 #the scale of axis
		 $paxis = "pa1.0f0.1:\"$xlabel[$i]\":/a20f5:\"$comp[$i]\"::.\"\":WSen"; #labels for axis
		 #		$saxis = "sa1of0.1"; #labels for axis

		 #data range
		 @minmaxt = ($data[0]->at(0)-0.2,$data[0]->at(-1)+0.2);
		 @minmax = (min($data[2*$i+1])*1000-5,max($data[2*$i+1])*1000+10);
		 $area = "$minmaxt[0]\/$minmaxt[1]\/$minmax[0]\/$minmax[1]"; #range
		 #----------------------------------------------------------------------------------------------	

		 #scaling of unsertainties----------------------------
		 #Fitting 
		 for (my $j =0 ; $j <= $#int ; $j++) { #fitting each year separetly
				#$unctemp = $data[2*$i+2];
				if ($j == 0) {
					 $data[2*$i+2] = $data[2*$i+2] * $sc[$i];
				}
				($coeff[$j], $coeffuns[$j]) = &smplsq($data[0]->index($int[$j]),
					 $data[2*$i+1]->index($int[$j]),
					 $data[2*$i+2]->index($int[$j]));
				($coeff[$j], $coeffuns[$j]) = ($coeff[$j]*${scaling},
					 $coeffuns[$j]*${scaling}); #changing to mm
				push @VDVB, ($coeff[$j],$coeffuns[$j]);

				#scaling of unsertainties----------------------------
				#if ($j == 0) {
				#($N) = dims($data[2*$i+1]->index(pdl[1..50]));
				#$sigmamed = sumover($data[2*$i+2]->index(pdl[1..50])*${scaling})/$N;
				#$ressq = ($data[2*$i+1]->index(pdl[1..50])*${scaling}-($data[0]->index(pdl[1..50])*($coeff[$j]->slice(1))+$coeff[$j]->slice(0)) )**2;
				#$oversigma = 1/(${scaling}*$data[2*$i+2]->index(pdl[1..50]));
#					print $sigmamed, " ", $ressq," ",$oversigma,"\n";

				#$sc = sqrt($N*sumover($ressq*$oversigma)/(($N-1)*sumover($oversigma)))/sqrt($sigmamed);
				#$unctemp = $data[2*$i+2]; 
				#$data[2*$i+2] = $unctemp * $sc;
				#print $sc,"\n";
#					print sqrt($ressq),"\n";
#					print $data[2*$i+2]*${scaling},"\n";
#					}

				#fiting again-----------------------------------------------------------------
#					($coeff[$j], $coeffuns[$j]) = &smplsq($data[0]->index($int[$j]),
#																					$data[2*$i+1]->index($int[$j]),
#																					$data[2*$i+2]->index($int[$j]));
#				  ($coeff[$j], $coeffuns[$j]) = ($coeff[$j]*${scaling},
#																		 $coeffuns[$j]*${scaling});
		 }	

		 #fitting everything

     my $cut = intersect(which($data[0] > 2007.1), which( $data[0] < 2007.7));
		 my $use = setops( which( 2010.3 > $data[0]), 'XOR' , $cut);
		 #my $use = which( 2010.3 > $data[0]);
		 ($coeff, $coeffuns) = &smplsq($data[0]->index($use),
				$data[2*$i+1]->index($use),
				$data[2*$i+2]->index($use));
		 ($coeff, $coeffuns) = ($coeff*${scaling},
				$coeffuns*${scaling}); #changing to mm
		 push @VDVA, ($coeff,$coeffuns);

		 #saving the average values of the last campaign
		 $ind = which($data[0] > 2008.0);
		 ($mean,$prms,$median,$min,$max,$adev,$rms) = 
		 stats($data[2*$i+1]->index($ind),1/($data[2*$i+2]->index($ind))**2);
		 push @MEAN, ($mean*${scaling},$rms*${scaling});

		 # preparing for ploting
		 @tmp = ([list $data[0]], [list $data[2*$i+1]], [list $data[2*$i+2]]);
		 $pdata="";
		 for (my $i=0;$i<=$#{$tmp[0]};$i++) { #input for psxy
				$pdata  .= sprintf "%s %s %s\n",  ${$tmp[0]}[$i], 
				${$tmp[1]}[$i]*${scaling}, 
				${$tmp[2]}[$i]*${scaling};
		 }

		 #ploting time series for each component
		 system ("psxy -JX$scale -R$area -B$paxis -Sc0.05 -W0.1/0 -G60 -Ey0.25c/2 -Y-9.01 -P -O -K << END >> ${psfile}\n${pdata}\nEND\n") == 0 or die "system failed: $?";

		 #plot the zero line
		 system("psxy -JX${scale} -R${area} -W2/0 -P -O -K << END >> ${psfile}\n$minmaxt[0] 0\n$minmaxt[1] 0\nEND\n");

		 #For plotting the fit	
		 for (my $j=0;$j<= $#int ;$j++) {
				@start = ($start-0.1, $coeff[$j]->at(0,0)+($start-0.1)*$coeff[$j]->at(1,0) );
				@stop = ($stop+0.1, $coeff[$j]->at(0,0)+($stop+0.1)*$coeff[$j]->at(1,0) );
				$label = sprintf( "%2.1f\261%2.1f mm/yr",$coeff[$j]->at(1,0),$coeffuns[$j]->at(1,0));

				#plot the fitt lines
				#if ($j == 0 ) {
				#					system("psxy -JX${scale} -R${area} -W2/0 -P -O -K << END >> ${psfile}\n$start[0] $start[1]\n$stop[0] $stop[1]\nEND\n");
				#put the fitt parameters on labels
				#			 system("pstext -JX${scale} -R0/8/0/4 -P -O -K << END >> ${psfile}\n$pos[$j] 13 0 0 5  $label\nEND\n") == 0 or die "system failed: $?";
				#					}
		 }
	}
	#finish the plot
	system "echo '0 0' | psxy -JX -R -Ss0.1 -G255 -P -O >> ${psfile}";

	my $file1 = "NORD0509enu.vel";
	#open for appending
	open FIL1, '>>' , $file1;
	# the fit for the hole period
	printf FIL1 "%2.6f %2.6f\t%2.4f %2.4f\t%2.4f %2.4f\t%2.4f %2.4f\t%s\n",
	$long,$lat,
	$VDVA[0]->at(1,0),$VDVA[1]->at(1,0),
	$VDVA[2]->at(1,0),$VDVA[3]->at(1,0),
	$VDVA[4]->at(1,0),$VDVA[5]->at(1,0),
	$station;
	#open for appending
	close FIL1;


	my $file2 = "NORD0506enu.vel";
	#open for appending
	open FIL2, '>>' , $file2;
	#printing the first part of the period
	 printf FIL2 "%2.6f %2.6f\t%2.4f %2.4f\t%2.4f %2.4f\t%2.4f %2.4f\t%s\n",
																						 							$long,$lat,
																 $VDVB[0]->at(1,0),$VDVB[1]->at(1,0),
																 $VDVB[6]->at(1,0),$VDVB[7]->at(1,0),
																 $VDVB[12]->at(1,0),$VDVB[13]->at(1,0),
																														$station;
	close FIL2;

	my $file3 = "NORD0809enu.vel";
	#open for appending
	open FIL3, '>>' , $file3;
	#printing the last part of the period
	printf FIL3 "%2.6f %2.6f\t%2.4f %2.4f\t%2.4f %2.4f\t%2.4f %2.4f\t%s\n",
																													$long,$lat,
																 $VDVB[2]->at(1,0),$VDVB[3]->at(1,0),
																 $VDVB[8]->at(1,0),$VDVB[9]->at(1,0),
															 $VDVB[14]->at(1,0),$VDVB[15]->at(1,0),
																														$station;

	my $file4 = "NORD0910enu.vel";
	#open for appending
	open FIL4, '>>' , $file4;
	#printing the last part of the period
	printf FIL4 "%2.6f %2.6f\t%2.4f %2.4f\t%2.4f %2.4f\t%2.4f %2.4f\t%s\n",
																													$long,$lat,
																 $VDVB[4]->at(1,0),$VDVB[5]->at(1,0),
																 $VDVB[10]->at(1,0),$VDVB[11]->at(1,0),
															 $VDVB[16]->at(1,0),$VDVB[17]->at(1,0),
																														$station;


	$year = 2007;
	$doy = 243;
	#	$day = $year+$doy/Days_in_Year($year,12); #calculating the difference at day $day
	$day = $stop; #calculating the difference at day $day
	my $Dday = ($day - $start)/2; 
	
	
	print join(" ",@MEAN),"\n";
	print $MEAN[4]-($VDVB[8]->(0)+$VDVB[8]->(1)*$day),"\n";

	$file4 = "NORDenu.dsp";
	#open for appending
	open FIL4, '>>' , $file4;
	#printing the last part of the period
						  printf FIL4 "%2.6f %2.6f\t%2.4f %2.4f\t%2.4f %2.4f\t%2.4f %2.4f\t%s\n",
																																		 $long,$lat,
																															$MEAN[0],$MEAN[1],
																															$MEAN[2],$MEAN[3],
	 														list($MEAN[4]-($VDVB[8]->(0)+$VDVB[8]->(1)*$day)),
										 list( sqrt(($MEAN[5] )**2 +  ( $VDVB[9]->(1)*$Dday )**2 )),
				 																															 $station;

#						  printf FIL4 "%2.6f %2.6f\t%2.4f %2.4f\t%2.4f %2.4f\t%2.4f %2.4f\t%s\n",
#																																		 $long,$lat,
#		list(($VDVB[2]->(0)+$VDVB[2]->(1)*$day)-($VDVB[0]->(0)+$VDVB[0]->(1)*$day)),
#					  list( sqrt( ($VDVB[3]->(1)*$Dday)**2 + ( $VDVB[1]->(1)*$Dday)**2 )),
#	 list( ($VDVB[6]->(0)+$VDVB[6]->(1)*$day)-($VDVB[4]->(0)+$VDVB[4]->(1)*$day)),
#						list(sqrt( ( $VDVB[7]->(1)*$Dday)**2 + ( $VDVB[5]->(1)*$Dday)**2 )),
#	list(($VDVB[10]->(0)+$VDVB[10]->(1)*$day)-($VDVB[8]->(0)+$VDVB[8]->(1)*$day)),
#				 list( sqrt(($VDVB[11]->(1)*$Dday )**2 +  ( $VDVB[9]->(1)*$Dday )**2 )),
#				 																															 $station;
	close FIL4;
	@MEAN = ();
  @VDVB = ();
	@VDVA = ();
}







# subroutines ------------------------------------------------------------


#----------------------------------------------------------------
#
# Least squere fit
#
#
#
#------------------------------------------------------------------

sub smplsq { 
	 my $M = zeroes(2,2);
	 my $V = zeroes(1,2);
	 my $X = shift();
	 my $Y = shift();
	 my $chi = shift();

	 $M(0,0) .= sumover($chi**-2); 
	 $M(0,1) .= sumover($X*$chi**-2); 
	 $M(1,0) .= $M(0,1);
	 $M(1,1) .= sumover($X**2*$chi**-2);
	 $V(0,0) .= sumover($Y*$chi**-2);
	 $V(0,1) .= sumover($Y*$X*$chi**-2);
	 if ($M(0,0)*$M(1,1)-$M(0,1)*$M(1,0)) {
			return ( transpose(matinv($M) x $V), sqrt(matinv($M)->diagonal(0,1)) );
	 } else {
			return pdl(0 , 0), pdl(0,0);
	 }
}
