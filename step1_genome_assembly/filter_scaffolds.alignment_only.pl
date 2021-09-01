#!/usr/bin/perl -w
use strict;

die "perl $0 [BLAST input] [blob_taxonomy.blob_out.blobDB.table.txt] [output file]\n Only consider alignment\n" unless @ARGV == 3;

my %hash;
my %blob;
open (LS, "$ARGV[1]") or die "$ARGV[1] $!\n";
while(<LS>){
	next if(/^#/);
	chomp;
	my $all = $_;
	my @line = split;
	if($line[4] < 70){
		$hash{$line[0]} = 1;
		$blob{$line[0]} = $all;
	}
}
close LS;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OT1, ">$ARGV[2].Blast") or die "$ARGV[2].tmp1 $!\n";
open (OT11, ">$ARGV[2].Blast.short") or die "$ARGV[2].tmp1 $!\n";
open (OT2, ">$ARGV[2].depth_noBlast") or die "$ARGV[2].tmp2 $!\n";
open (OT3, ">$ARGV[2].blast_depth") or die "$ARGV[2].Pem $!\n";
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split;
	my $pem = 0;
	if(($line[4] >= 95)&&($line[5] >= 1000)){
		print OT1 "$all\n";

		if($hash{$line[2]}){
			print OT3 "$all\n";
			$pem = 1;
			$hash{$line[2]} = 2;
		#}else{
		#	print OT1 "$all\n";
		}
	}elsif(($line[0] < 1000)&&($line[4] >= 95)){
		print OT11 "$all\n";
	}
}
close IN;

foreach my $k (sort keys %hash){
	if($hash{$k} == 1){
		print OT2 "$blob{$k}\n";
	}
}

print "$ARGV[2].Blast\n$ARGV[2].Blast.short\n$ARGV[2].depth_noBlast\n$ARGV[2].blast_depth\n";
