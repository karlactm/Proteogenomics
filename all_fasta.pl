#!/usr/bin/env perl
#use strict;

my $in = $ARGV[0];
my $out = $ARGV[1];
opendir(my $dir, $in) || die "Can't opendir $in: $!";
my @fasta = grep { /\.faa$/ } readdir ($dir);
my $path_fasta = undef;

foreach (@fasta) {
    $path_fasta .=  $in.$_." ";
}

my $arg = $out." ".$in." ".$path_fasta;
print $arg."\n\n";
system("./round_rand.sh $arg");
