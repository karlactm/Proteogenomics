#!/usr/bin/env perl
use strict;
use FASTAParser;

my @parsers = ();
my $output_ref_path = $ARGV[0];
my $output_parser_path = $ARGV[1];
my $ref_fasta = FASTAParser->new ($ARGV[2]);
if (!$ref_fasta->parse) {
    print "-> Insira um arquivo no formato FASTA!";
    exit 1; 
}
my $ref_path = $ref_fasta->getFilePath;
my $ref_fileName = $ref_fasta->getFileName;

for (my $i = 3; $i <= $#ARGV ; $i++) {
    my $parser = FASTAParser->new($ARGV[$i]);
    if (!$parser->parse) {
        print "-> Insira um arquivo no formato FASTA!";
        exit 1; 
    }
    push @parsers, $parser;
}
my $y = 0;
print "->> ".scalar(@parsers)."\n";
for (my $i = 0; $i < scalar(@parsers); $i++) {
    my $parser_path = $parsers[$i]->getFilePath;
    my $parser_fileName = $parsers[$i]->getFileName;
    my $parser_output = $output_parser_path."/".$y."_PARSER_GCA_".$parser_fileName."_GCA_".$ref_fileName."_output.txt";
    my $ref_output = $output_ref_path."/".$y."_REF_GCA_".$ref_fileName."_GCA_".$parser_fileName."_output.txt";
    `makeblastdb -in $ref_path -dbtype prot -title GCA_$ref_fileName`;
    `makeblastdb -in $parser_path -dbtype prot -title GCA_$parser_fileName`;   
    `blastp -query $ref_path -db $parser_path -num_alignments 3 -num_threads 8 -out $ref_output -outfmt '7 qseqid sseqid pident ppos length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp'`;
    `blastp -query $parser_path -db $ref_path -num_alignments 3 -num_threads 8 -out $parser_output -outfmt '7 qseqid sseqid pident ppos length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp'`;
   $y++;
}  

