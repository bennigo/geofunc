#!/usr/bin/perl -w
#-----------------------------------------
#
#
#
#-------------------------------------------------------
#
use strict;
use PDL;
use PDL::NiceSlice;
#use PDL::Core ':Internal';
use PDL::Slatec qw( matinv);
use File::Basename qw(dirname basename fileparse);
#use Date::Pcalc qw(Days_in_Year);
use Geo::Functions qw(rad_deg);
use PDL::Graphics::PLplot;
use PDL::Fit::Linfit;
use PDL::Fit::LM;
use Geo::Constants qw(PI);
#use PDL::Primitive;
use PDL::LinearAlgebra qw(mnorm  msolve);
use PDL::Fit::Polynomial;
use Time::Elapse;
use PDL::Slatec qw(polyfit polyvalue); 
use Math::Random qw/random_permuted_index/;
use Carp;


require "geo.pm";

my @data;
my ($paxis,$saxis, $area);
my ($mean,$prms,$median,$min,$max,$adev,$adev1,$rms,$omega,$halsminmax);
my ($ind,$ser_leng,$min_max);
my ($station, $long, $lat, $height); 
my ($pl,$sig,$start_t,$stop_t,$days,$sc,$scale);
my ($filbase, $path,$suff,$psfile);
my (@xyz,@stations,@minmax,@minmaxt,@coeff,@start,@stop);
my ($yfit, $coeffs,$secular,$use,$rem,$ind1,$halsvel,$res_ind,$residual);
my ($eventsta,$i0,$i1,$time0,$time1,$fitFuncs,$yfit0, $coeffs0,$yfit1, $coeffs1,$tmpdata,$head,$covar,$iters); #for the fit
my ($pdata,@tmp,@head);
my ($ndeg, $r, $ierr, $a,$v,$Funcs);
my (@secular,@extra);

# scaling to mm or whatever
my $scaling = 1000;
my @comp = qw/East North Up/;
my @xlabel = ("","","Year");
my @pos = ("0.8 3.6", "5.6 0.6" );
my @sc = (5,4,3);

$scale = '18/8';

my $swarm=pdl(2008.224043716,2008.25136612);
#$eventsta = 2007.4;
$eventsta = 2007.33137;
my $vikd_sec=pdl(0.00175545495752973, 0.00410243260798486, 0.0161129116270589);
my $vikd_var=pdl(0.00763232329656145, 0.00863380587599824, 0.00811667295410776);
#my $alfd_sec=pdl(0.00179518917855505,0.00187397582058953,0.01815677384357);
my $alfd_sec=pdl(0.00179518917855505,0.00160397582058953,0.02530677384357);
my $alfd_var=pdl(0.00540273123446652, 0.0068435970071225, 0.0061381003857843);

my $amp=pdl ([1.5587921,0.25702068,6.1270153])/1000;
my $phi=pdl ([1.4367988,1.048543,0.24641968]);

my $amp_bru=pdl ([2.0819406374,1.084016699,8.90617785])/1000;
my $phi_bru=pdl ([1.2870044,1.4126034,0.26338735]);

my $amp_karv=pdl ([2.6330469,1.7921579,5.1403398])/1000;
my $phi_karv=pdl ([1.6474249,1.1399135,0.32451402]);

#my $intervals = ;
foreach  my $ARG (@ARGV) {
#	 print $ARG,"\n";
		
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
#----------------------------------------------
($filbase, $path,$suff ) = &fileparse($ARG,qw{enu.disp});
@data = rcols $ARG;
push @stations,  $filbase;
	$psfile = "${filbase}\.ps"; # psfile for gmt output ####CHECK enu  and enu.disp make more general########3

	#initiating plot
	system ("echo '0 0' | pstext  -JX${scale} -R-1/1.5/-1/3 -Y29.1 -P -K << EOF > ${psfile}\n0.3 -0.8 16 0 4 TC ${filbase}\nEOF\n");

$omega=2*PI;
$scaling=1000;

$ser_leng=$data[0]->dim(0); #number of ephocs
$start_t = &min($data[0]);
$stop_t = &max($data[0]);

$fitFuncs = cat ones($ser_leng),$data[0],cos($data[0]*$omega),sin($data[0]*$omega);

$days=0.002739726*7;
$ind = which($data[0] < $data[0]->at(0)+$days);


#Fixing the Eurasian plate----------------------
#Altamimi et al 2007 Absolute rotation pole for Eurasia in ITRF2005 (56.330 deg,-95.979 deg,0.261 deg/my)
my @eurpole = map &rad_deg($_),(56.330,-95.979,0.261/1e6);
@eurpole = @{&ellecc([$eurpole[0],$eurpole[1]],[0,0,$eurpole[2]])};
my $nuvvel = &eccell(&xyzell(&gdatum,\@xyz),&rotpole("",\@xyz,\@eurpole),1);
$data[1] .= $data[1] - $nuvvel(0,1)*($data[0]-$start_t);
$data[3] .= $data[3] - $nuvvel(0,0)*($data[0]-$start_t);
#-------------------------------------------------------------------	


for(my $i=0;$i <= 2;$i++) {
   
   # Mean of the first few days as a refference 
   $sig = 1/($data[2*$i+2]**2); #forming 1/chi**2
	 $sig = $sig/sumover($sig);
   #Subtacting the mean of the first ephcs #----------------------------------------------------------
   ($mean,$prms,$median,$min,$max,$adev,$rms) = stats($data[2*$i+1]->index($ind),$sig->index($ind));
   $data[2*$i+1]=($data[2*$i+1]-$mean);

			#the scale of axis
			$paxis = "pa2.0f0.5:\"$xlabel[$i]\":/a10f1:\"$comp[$i]\"::.\"\":WSen"; #labels for axis
	 

			#$data[$i]*=$scaling; #changing to mm
			#$data[$i+1]*=$scaling; #changing to mm
#-----removing annual signal---------------------------------------------------------------------------------------------
#	 $data[2*$i+1]*=$scale; #changing to mm
#	 $data[2*$i+2]*=$scale; #changing to mm
			#	$data[2*$i+1]=$data[2*$i+1]-$amp($i)*cos($omega*$data[0]+$phi($i));
	 
	 $ind1 = setops( which( 2004.0 < $data[0]), 'AND', which( 2010.0 > $data[0]) ); 
	 print $ind1;
#   $Funcs = cat ones($ind1->dims(0)),$data[0]->index($ind1),cos($data[0]*$omega+$phi($i))->index($ind1);
   $Funcs = cat ones($ind1->dims(0)),$data[0]->index($ind1);
	 ($yfit, $coeffs) = linfit1d $data[0]->index($ind1), $data[2*$i+1]->index($ind1), $Funcs, { Weights=>$sig->index($ind1) };
	 #$Funcs = cat ones($ser_leng),$data[0],cos($data[0]*$omega+$phi($i));
	 $Funcs = cat ones($ser_leng),$data[0];
	 $data[2*$i+1] .= ( $data[2*$i+1]-$coeffs(1) x $Funcs(:,1) )-> clump(-1) ; # subtracting anual sicle
	 # $data[2*$i+1] .= ( $data[2*$i+1]-$amp($i)*cos($omega*$data[0]+$phi($i)) ) -> clump(-1) if ($station eq "SAUD");
	 #$data[2*$i+1] .= ( $data[2*$i+1]-$amp_bru($i)*cos($omega*$data[0]+$phi_bru($i)) ) -> clump(-1) if ($station eq "BRUJ");
	 #$data[2*$i+1] .= ( $data[2*$i+1]-$amp_karv($i)*cos($omega*$data[0]+$phi_karv($i)) ) -> clump(-1) if ($station eq "KARV");
#-----------------------------------------------------
			
	 
	
			#	if ( dims($secular) ) {
			# $secular = $secular->glue(1,$coeffs);
			#} else {
			# $secular = $coeffs;
			#}

			#----- fitting ---------------------------------------------
			
			my @color = ("red","green","yellow");
      my $tmt0=$data[0] - $eventsta;
      $tmt0->where($tmt0 < 0).=0;
			$use = which($data[0] < 2010.00);
			#my $cut = intersect(which($data[0] > 2006.0), which( $data[0] < 2007.2)) if ($filbase eq "BRUJ");
			#$use = setops( which( 2008.35 > $data[0]), 'XOR' , $cut) if  ($filbase eq "BRUJ");
			$Funcs = cat ones($use->dims(0)),$data[0]->index($use),$tmt0->index($use);
			unless ($station eq "UPP3") {
				 ($yfit, $coeffs) = linfit1d $data[0]->index($use), $data[2*$i+1]->index($use), $Funcs, { Weights=>$sig->index($use) };
				 ($yfit,$coeffs,$covar,$iters) = lmfit $data[0]->index($use), $data[2*$i+1]->index($use), $sig->index($use), \&twolinefit, $coeffs, {Maxiter => 900, Eps => 1e-9} ;
			   $v = bootstr($data[0]->index($use), $data[2*$i+1]->index($use), $sig->index($use),$coeffs,1000);
				 #print join( ", ", bootstrap($data[2*$i+1]->index($use),1000)),"\n";
			   ($mean,$prms,$median,$min,$max,$adev,$rms) = stats($v(0,:));
				 push @secular, ($coeffs->at(1,0),$adev->at(0));
			   ($mean,$prms,$median,$min,$max,$adev,$rms) = stats($v(1,:));
				 push @extra, ($coeffs->at(2,0),$adev->at(0));
			}
			$Funcs = cat ones($data[0]->dim(0)),$data[0];
			$data[2*$i+1] .= ($data[2*$i+1]-$coeffs(0:1) x $Funcs) -> clump(-1);
      print join(" ",$coeffs),"\n";
		  
			
 			# removing outliers ---------------------
#			$res_ind=which(abs($data[2*$i+1]-$coeffs(2)x$tmt0)<0.005) and which($data[0] < 2008.3) unless ($i == 2);
#			$res_ind=which(abs( ( $data[2*$i+1]-$coeffs(2)x$tmt0)<0.005 and $data[0] < 2008.3 )) unless ($i == 2);
#      $res_ind=which(abs($data[2*$i+1]-$coeffs(2)x$tmt0)<0.02) if ($i == 2);
#			$data[0] = $data[0]->index($res_ind);
#			$data[1] = $data[1]->index($res_ind);
#      $data[2] = $data[2]->index($res_ind);     
#			$data[3] = $data[3]->index($res_ind);
#			$data[4] = $data[4]->index($res_ind); 
#			$data[5] = $data[5]->index($res_ind);
#			$data[6] = $data[6]->index($res_ind);
#			$tmt0 = $tmt0->index($res_ind);
#		  $ser_leng=$data[0]->dim(0); #number of ephocs

			# data points on the fit line
			($use,$rem)=which_both($eventsta > $data[0]);
			$yfit=$coeffs(2) x $tmt0;
			$yfit=$yfit->index($use)->glue(0,$coeffs(2) x pdl(0))->glue(0,$yfit->index($rem));
			my $x = $data[0]->index($use)->glue(0,pdl $eventsta)->glue(0,$data[0]->index($rem)); 

			


	 		#data range
			@minmaxt = ($data[0]->at(0)-0.2,$data[0]->at(-1)+0.2) if ($i==0);
			@minmax = (min($data[2*$i+1])*1000-5,max($data[2*$i+1])*1000+10);
			$area = "$minmaxt[0]\/$minmaxt[1]\/$minmax[0]\/$minmax[1]"; #range
			#----------------------------------------------------------------------------------------------	

      #-----Ploting data----------------------------------------------------------------------------------------------------------
			# preparing for ploting
			@tmp = ([list $data[0]], [list $data[2*$i+1]], [list $sc[$i]*$data[2*$i+2]]);
			$pdata="";
			for (my $i=0;$i<=$#{$tmp[0]};$i++) { #input for psxy
				 $pdata  .= sprintf "%s %s %s\n",  ${$tmp[0]}[$i], 
				 												${$tmp[1]}[$i]*${scaling}, 
																${$tmp[2]}[$i]*${scaling};
			}
				 #ploting time series for each component with error bars
				 system ("psxy -JX$scale -R$area -B$paxis -Sc0.05 -W0.1/0 -G60 -Ey0.25c/2 -Y-9.01 -P -O -K << END >> ${psfile}\n${pdata}\nEND\n") == 0 or die "system failed: $?";
				 #ploting time series for each component without error bars
			   #plot the zero line
				 system("psxy -JX${scale} -R${area} -W2/0 -P -O -K << END >> ${psfile}\n$minmaxt[0] 0\n$minmaxt[1] 0\nEND\n");
#----------------------------------------------------------------------------------------------------------------------


			
			# preparing for ploting fit lines
			@tmp = ([list $x], [list $yfit]);
			$pdata="";
			for (my $i=0;$i<=$#{$tmp[0]};$i++) { #input for psxy
				 $pdata  .= sprintf "%s %s\n",  ${$tmp[0]}[$i], 
				 												${$tmp[1]}[$i]*${scaling}, 
			}
#		  system("psxy -JX${scale} -R${area} -W2/0   -P -O -K << END >> ${psfile}\n$minmaxt[0] 0\n$minmaxt[1] 0\nEND\n");
			#system ("psxy -JX$scale -R$area -B$paxis  -W2/$color[0]   -P -O -K << END >> ${psfile}\n${pdata}\nEND\n") == 0 or die "system failed: $?";

			#@tmp = ([list $swarm], [list $coeffs(2) x pdl($swarm-$eventsta) ]);
			@tmp = ([ ( $swarm->at(0),$swarm->at(0) )],[ (-1, 1)]);
			$pdata="";
			for (my $i=0;$i<=$#{$tmp[0]};$i++) { #input for psxy
				 $pdata  .= sprintf "%s %s\n",  ${$tmp[0]}[$i], 
				 												${$tmp[1]}[$i]*${scaling}, 
			}
		  print $pdata,"\n";
			#system ("psxy -JX$scale -R$area -B$paxis -Sx0.2 -W$color[1] -P -O -K << END >> ${psfile}\n${pdata}\nEND\n") == 0 or die "system failed: $?";
			#system ("psxy -JX$scale -R$area -B$paxis -W9/$color[1] -P -O -K << END >> ${psfile}\n${pdata}\nEND\n") == 0 or die "system failed: $?";

			#@tmp = ([list $eventsta], [list $coeffs(2) x pdl(0) ]);
			@tmp = ([list $eventsta], [list $coeffs(2) x pdl(0) ]);
			$pdata="";
			$pdata  .= sprintf "%s %s\n",${$tmp[0]}[0], 
				 												${$tmp[1]}[0]*${scaling};
		  print $pdata,"\n";
		  system ("psxy -JX$scale -R$area -B$paxis -Sh0.3 -W9/$color[1]  -P -O -K << END >> ${psfile}\n${pdata}\nEND\n") == 0 or die "system failed: $?";
			

			# special for SAUD and KARV -------------------------------------------------------------------
			if ( ( ($station eq "KARV") or ($station eq "SAUD") ) and ($i != 2 ) ) {
				 @tmp = ([list $halsminmax], [list $halsvel(0:1) x pdl(ones(2),$halsminmax) ]);
				 $pdata="";
				 for (my $i=0;$i<=$#{$tmp[0]};$i++) { #input for psxy
						$pdata  .= sprintf "%s %s\n",  ${$tmp[0]}[$i], 
				 												           ${$tmp[1]}[$i]*${scaling}, 
				 }
				 system ("psxy -JX$scale -R$area -B$paxis -W9/$color[2] -P -O -K << END >> ${psfile}\n${pdata}\nEND\n") == 0 or die "system failed: $?";
			}
 			#---------------------------------------------------------
			#
}
#finish the plot
system "echo '0 0' | psxy -JX -R -Ss0.1 -G255 -P -O >> ${psfile}";

	open FIL, $ARG;
	@head = <FIL>;
	$head = join "", @head[0..7] ;
	close FIL;

	#print $data[0]->dim(0),"\n"; #number of ephocs
	my $datafile = "${station}_filt\.dat";
	wcols "%.4f\t%7.4f %7.4f\t%7.4f %7.4f\t%7.4f %7.4f", @data, $datafile, { HEADER => $head };
  
	#print $head;
#  print join (" ",@secular)," ",$station,"\n";
	my $file = "upp_secular.vel";
	open FIL, '>>' , $file;
	 	 printf   "%2.6f %2.6f\t%2.4f %2.4f\t%2.4f %2.4f\t%2.4f %2.4f\t%s\n",
		 																						           $long,$lat,
																										$secular[0],$secular[1],
																										$secular[2],$secular[3],
																										$secular[4],$secular[5],
																														 $station;
	close FIL;
#	print join (" ",@extra)," ",$station,"\n";
	$file = "upp_extra.vel";
	open FIL, '>>' , $file;
	 	 printf FIL  "%2.6f %2.6f\t%2.4f %2.4f\t%2.4f %2.4f\t%2.4f %2.4f\t%s\n",
		 																						           $long,$lat,
																										$extra[0],$extra[1],
																										$extra[2],$extra[3],
																										$extra[4],$extra[5],
																														 $station;
	close FIL;

	@secular = ();
	@extra = ();

}


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
	 return ( transpose(matinv($M) x $V), sqrt(matinv($M)->diagonal(0,1)) );
}

sub  twolinefit {
	 my ($x,$par,$ym,$dyda) = @_;
	 my $eventsta =2007.4;
	 my ($b,$vs,$ve) = map { $par->slice("($_)") } (0..2);
   my $tmt0;

	 $tmt0=$x - $eventsta;
   $tmt0->where($tmt0 < 0).=0;

	 $ym .= $b + $vs*$x + $ve*$tmt0;
	 my (@dy) = map {$dyda -> slice(",($_)") } (0..2);

	 $dy[0] .= 1;
	 $dy[1] .= $x;
	 $dy[2] .= $tmt0;


}
sub  threelinefit {
	 my ($x,$par,$ym,$dyda) = @_;
	 my $eventsta1 =2007.4;
	 my $eventsta2 = 2008.25136612;
	 my ($b,$vs,$ve,$vds) = map { $par->slice("($_)") } (0..3);
   my $tmt0;
   my $tmt1;

	 $tmt0=$x - $eventsta1;
   $tmt0->where($tmt0 < 0).=0;
   $tmt0->where($tmt0 > $eventsta2).= $eventsta2 - $eventsta1;
	 $tmt1=$x - $eventsta2;
   $tmt1->where($tmt1 < 0).=0;

	 $ym .= $b + $vs*$x + $ve*$tmt0 + $vds*$tmt1;
	 my (@dy) = map {$dyda -> slice(",($_)") } (0..3);

	 $dy[0] .= 1;
	 $dy[1] .= $x;
	 $dy[2] .= $tmt0;
	 $dy[3] .= $tmt1;

}
sub  linefit {
	 my ($x,$par,$ym,$dyda) = @_;
	 my ($b,$vs) = map { $par->slice("($_)") } (0..1);
	 $ym .= $b + $vs*$x;
	 my (@dy) = map {$dyda -> slice(",($_)") } (0..1);

	 $dy[0] .= 1;
	 $dy[1] .= $x;


}

sub bootstrap {
	my ($x, $rounds) = @_;

	my @N = dims($x);
	print join( ", ",@N),"\n";
	croak "Can only use a single-dimensional xset." unless @N == 1;

	my $davg = avg($x);

	$rounds ||= 1000;
	my $b = int(0.1 * $N[0]+0.5);
	my $mask = zeroes($N[0]);

	my $d = zeroes($rounds);

	Time::Elapse->lapse(my $now = "processing");
	for my $i(1..$rounds)
	{
		my @index = random_permuted_index($N[0]);
		$mask .= 0;
		$mask->index(pdl(@index[0..$b-1])) .= 1;
		my $r = $x->where($mask);
		$d->set($i-1, avg($r) - $davg);
	}
	print "Time Wasted: $now\n";
	my $u = $davg - pct($d, 0.025);
	my $l = $davg - pct($d, 0.975);

	return (($davg-$l) + ($u-$davg))/2 unless wantarray;
	return ($l, $davg, $u);
}

sub bootstr {
	 my ($x ,$y ,$wt ,$m ,$iter) = @_;
	 my ($ix);
	 my $n_image=$y->dim(0);
	 my $v=zeros(2,$iter);
	 my $fitdim = $m->getdim(0);
	
	 Time::Elapse->lapse(my $now = "processing");
	for (my $i=0;$i< $iter;$i++) {
		   $ix=ceil(random($n_image,1)*$n_image-1);
			 #(undef, $m) = fitpoly1d $x($ix1),$y($ix1),2;
			 unless ($fitdim == 2) {
	 		 (undef,$m,undef,undef) = lmfit $x($ix), $y($ix), $wt($ix), \&twolinefit, $m, {Maxiter => 900, Eps => 1e-9} ;
			 #(undef, $m) = linfit1d $x1($ix), $y($ix), $G, { Weights => $wt($ix) };
		 	  $v(:,$i).=$m(1:2);
		 } else {
	 		 (undef,$m,undef,undef) = lmfit $x($ix), $y($ix), $wt($ix), \&linefit, $m, {Maxiter => 900, Eps => 1e-9} ;
		 	  $v(:,$i).=$m(1);
		 }
		 }
		 print "Time Wasted: $now\n";
	return $v;
}
