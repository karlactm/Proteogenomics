#!/usr/bin/env perl
use strict;

my $in = $ARGV[0];
my $out = $ARGV[1];
opendir(my $dh, $in) || die "Can't opendir $in: $!";
my @lineages = grep { /[^\.]/ } readdir($dh);
closedir $dh;
my @rand = (5, 10, 20, 30, 65);
for (my $i = 0; $i < scalar (@lineages); $i++) {
    my $path = $in."/".$lineages [$i]."/";
    for (my $j = 0; $j < scalar (@rand); $j++) {
        my @faa = ();
        opendir(my $ln, $path) || die "Can't opendir $path: $!";
        my @fasta = grep { /\.faa$/ } readdir ($ln);
        foreach (@fasta) {
            if ($_ =~ /GCA.*/) {
                push (@faa, $path.$_);
            }
        } 
        my $arg = undef;
        $arg .= $out."/".$lineages [$i]." ";
        $arg .= $in."/".$lineages [$i]." ";
        for (my $k = 0; $k < $turn; $k++) {
             my $random = int rand ($#faa);
             $arg .= $faa[$random]." ";
             splice @faa, $random, 1;
        }
        system("./round_rand.sh $arg");
        $arg = undef;       
   }      
}
