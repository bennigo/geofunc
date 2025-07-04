#!/usr/bin/perl -w
#--------------------------------------------------------
#
#
#
#
#
#
#---------------------------------------------------------

use strict;
use Date::Calc qw(Moving_Window);

require "displ.pm";

my @flist = ();
my @stlist = ();
my $enu = 1;
my $dir = "../NORD/";

my $LIST;
my $listfile = "lists/nord-camp.list";

#reading in list of stations
open $LIST, $listfile or die "can't open file $listfile: $!\n";
while (<$LIST>) {
	 next if /#/;
	 chomp;
	 push @stlist, $_;
}
@stlist = qw/OLAF MASK/;
#@stlist = qw/DYNG/;

my $files = [glob $dir."F1_*.CRD"];
$files = [ sort( { Moving_Window( substr($a,-10,2) ) <=> Moving_Window( substr($b,-10,2)) } @{$files}) ];
my $outfile;

#my $reffile = "F1_052410.CRD";
my ($station,$refst) = (undef,"DYNG");
my $ref = 0;
if (defined $refst) {
	 print STDERR "Reference station $refst used\n","-"x46,"\n";
} else {
	 print STDERR "No reference station selected\n","-"x46,"\n";
}

foreach (@stlist) { #looping through stations
		
		$station = $_;
		$refst = $_ unless( $ref == 1) && (defined $refst );

		if ($refst eq $station) { #making outfile date +%j -d20070401

			 $outfile = "nord/${station}enu\.disp";
		} else {
			 $outfile = "nord/${station}-${refst}enu\.disp";
		}
		#next if (-s $outfile); #next if the file exists
		
		@ARGV = @{$files};
		while (<>) {
			 push @flist, $ARGV if grep( {/$station/ and /\ \ \ \ [AW]/} ($_));
		}

		unless (@flist) {
			 print STDERR $station," is not available\n";
		   print STDERR "-"x46,"\n";
			 next;
		}
    
		unless ($station eq $refst) {
			 @ARGV = @flist;
			 @flist = ();
		   while (<>) {
					push @flist, $ARGV if grep( {/$refst/ and /\ \ \ \ [AW]/} ($_));
			 }

			 unless (@flist) {
					print STDERR "Reference station ",$refst," is not available\n";
		   		print STDERR "-"x46,"\n";
			 		next;
			 }

		}
    


		print STDERR "Reference epoch for $station in file ",$flist[0],"\n" if  @flist;
		print STDERR "Displacement file writen to ",$outfile,"\n";
		print STDERR "-"x46,"\n";
			
			open STDOUT, ">$outfile";
			&displ($flist[0],$refst,$station,\@flist,$enu);
} continue {
	 @flist = ();
}
