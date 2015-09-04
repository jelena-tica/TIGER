#!/usr/bin/perl

use warnings;
use strict;

my $file_regions = $ARGV[0] or die "Can't open!\n";
my $file_blat = $ARGV[1]  or die "Can't open!\n";
my $output = $ARGV[2] or die "Can't open!\n";
#my $stat = $ARGV[3]  or die "Can't open!\n";

open REGIONS, $file_regions;
open BLAT, $file_blat;
open OUTPUT, ">".$output;
#open STATOUT, $stat;

my @blat = <BLAT>;

while (<REGIONS>) {
	my $regline = $_;
	chomp $regline;
	#print $regline;
	my @regparts = split(/\t/,$regline);	
	my $count=0;
	
	foreach my $blatline (@blat) {
		chomp $blatline;
	#print $blatline; exit;
		my @blatparts = split(/\t/,$blatline);
	#print $blatparts[4]."\t".$regparts[0]."\t".$blatparts[5]."\t".($regparts[1]-100)."\t".($regparts[2]+100)."\n"; exit;
		if ($blatparts[4] eq $regparts[0]) {
			if ($blatparts[5]>=($regparts[1]-100)) {
				if ($blatparts[5]<=($regparts[2]+100)) {
	
			print OUTPUT $regparts[0]."\t".$regparts[1]."\t".$regparts[2]."\t".$blatparts[0]."\t".$blatparts[1]."\t".$blatparts[2]."\t".$blatparts[3]."\t".$blatparts[4]."\t".$blatparts[5]."\t".$blatparts[6]."\t".$blatparts[7]."\t".$blatparts[8]."\t".$blatparts[9]."\t".$blatparts[10]."\t".$blatparts[11]."\t".$blatparts[12]."\t".$blatparts[13]."\t".$blatparts[14]."\t".$blatparts[15]."\t".$blatparts[16]."\t".$blatparts[17]."\n";
			$count++;
		}
		}
		}
		else {next;}
		
		}
		#print STATOUT $regparts[0]."\t".$regparts[1]."\t".$regparts[2]."\t".$count."\n";
	}


close REGIONS;
close BLAT;
close OUTPUT;
#close STATOUT;
