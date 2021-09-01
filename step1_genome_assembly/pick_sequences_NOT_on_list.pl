#!/usr/bin/perl -w
use strict;

die "perl $0 [fasta] [name list] [output file]\n Cha aphid version\n" unless @ARGV == 3;

my %hash;
open (LS, "$ARGV[1]") or die "$ARGV[1] $!\n";
while(<LS>){
	chomp;
	s/\>//;
	my @line = split;
	$hash{$line[0]} = 1;
}
close LS;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[2]") or die "$ARGV[2] $!\n";
$/=">";
<IN>;
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split /\n+/;
	my $id = shift @line;
	my $seq = join "", @line;
#	my @names = split /\_/, $line[0];
#	if($hash{$names}){
	my @id = split /\s+/, $id;

	if($hash{$id[0]}){

	}else{
		$/="\n";
#		print OUT ">$id\n$seq\n";
		print OUT ">$all";
		$/=">";
	}

}
close IN;
