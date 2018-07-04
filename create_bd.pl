#!/usr/bin/env perl
use strict;

my $dir = $ARGV[0];
print $dir."\n";
opendir(my $dh, $dir) || die "Can't opendir $dir: $!";
my @lineages = grep { /[^\.]/ } readdir($dh);
print join ("\n", @lineages)."\n";
closedir $dh;
my @rand = (5, 10, 20, 30, 65);
for (my $i = 0; $i < scalar (@lineages); $i++) {
    for (my $j = 0; $j < scalar (@rand); $j++) {
        my $turn = $rand [$j];
        my $arg = $dir."/".$lineages [$i]."/".$turn."/output_find_homologous/"; 
        system ("mkdir $arg/output_peptide_trip/");
        system("./peptide_trip.pl " . $arg);
        $arg = undef;       
   }        
}
