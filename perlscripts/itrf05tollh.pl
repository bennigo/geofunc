#!/usr/bin/perl -w
# converts geocentric coordinates to geographic coordinats
# Assumes bernese coordinate file i.e. ITR05_R.CRD
#
# usage ./itrf052llh.pl "coordinate file(s)"
#
# 	02.10.07:  BGO
#
#--------------------------------------------

use strict;
use Geo::Functions qw( deg_rad dms_deg  ); 
require "geo.pm";

my @xstat;
my ($xyzell,$station,$stat,$flag);
while (<>) {
	 next if $. <= 6; # WARNIG can only handle one file
	 next if $_ =~ /^\s*$/;
	 chomp; # (@xstat = split);
	 ( undef,undef,$station,$xstat[0], $xstat[1], $xstat[2], undef) = 
	 split /(^\s*\d+)\s+([\d\D]{14})\s+(-?\d+\.\d{4})\s+(-?\d+\.\d{4})\s+(-?\d+\.\d{4})\s+(\w)?$/, $_; 
	 ($stat,undef) = split / /,$station;
															 # 			print $stat, " ",$xstat[0]," ", $xstat[1], " ", $xstat[2],"\n";
	 		$xyzell =  &xyzell(&gdatum,\@xstat);
				print deg_rad(@{$xyzell}[1])," ",deg_rad(@{$xyzell}[0])," ",@{$xyzell}[2]," ",$stat,"\n";
} #continue {
#close ARGV if eof;
#}
