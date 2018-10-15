#!/usr/bin/env perl
use strict;
use FASTAParser;
use BLASTParser;

my %fasta = ();
my $number = $ARGV[0];
my $fasta_subpath = $ARGV[1];
my $blast_subpath = $ARGV[2];
my $homologous_subpath = $ARGV[3];
my $directory = $fasta_subpath."/";
opendir(directory, "$directory");
my @list = readdir(directory);
closedir(directory);
foreach (@list) {	
    if($_=~/(.faa)$/){
        my $parser = FASTAParser->new($directory.$_);
        if (!$parser->parse) {
            print "-> Insert FASTA file!";
            exit 1; 
        }
        my $proteins = $parser->getProteins;

        for(my $i = 0; $i < scalar (@$proteins); $i++) {
            my $format_sequence = &formatFasta(@$proteins[$i]->getSequence);
            $fasta{@$proteins[$i]->getId} = ">".@$proteins[$i]->getHeader."\n".$format_sequence;     
        }              
   }	
}

my $ref_proteins = undef;
my $parser_proteins = undef;
my %proteins_not_found_ref = ();
my %proteins_found = ();
my $count = 0;

for (my $i = 4; $i <= $#ARGV ; $i++) {
    my %proteins_not_found_parser = ();
    #print "-> $i $ARGV[$i]\n";
    $ref_proteins = BLASTParser->new($ARGV[$i]);
    if (!$ref_proteins->parse) {
        print "-> Insert blast file!";
        exit 1; 
    }
    my ($ref, $not_found_ref) = &align ($ref_proteins);
    my %ref_id = %$ref;
    %proteins_not_found_ref = %$not_found_ref;
    my @file_name = split (/\//,$ARGV[$i]);
    @file_name = split (/\_/,$file_name[$#file_name]);
    my $parser_path = $blast_subpath."/PARSER_".$number."/".$file_name[0]."_PARSER_GCA_".$file_name[5]."_GCA_".$file_name[3]."_output.txt";

    $parser_proteins = BLASTParser->new($parser_path);
    if (!$parser_proteins->parse) {
        print "-> Insert blast file!";
        exit 1; 
    }  
    my ($parser, $not_found_parser) = &align ($parser_proteins);  
    my %parser_id = %$parser;
    %proteins_not_found_parser = %$not_found_parser;
    foreach my $key (keys %ref_id) {
        if (exists $parser_id {$key}) {
            foreach my $value (@{$parser_id{$key}}) {
                   if ($ref_id {$key} eq $value) {
                        if (exists $proteins_not_found_parser {$parser_id{$key}}) {
                            delete $proteins_not_found_parser {$parser_id{$key}};
                        }
                        push @{$proteins_found{$key}}, $value;
                    } else {
                        $count++;
                        $proteins_not_found_parser {$value} = $file_name[5];
                    }
             }
             if (scalar (@{$parser_id{$key}}) == $count) {
                 $proteins_not_found_ref {$key} = $file_name[3];
             }
         } else {
               $proteins_not_found_ref {$key} = $file_name[3]; 
         }
         $count = 0;
     }
                                                               
    foreach my $key (keys %parser_id) {
        if (!(exists $ref_id {$key})) {
            $proteins_not_found_ref {$key} = $file_name[3]; 
         }
         foreach my $value (@{$parser_id{$key}}) {
            if (exists $proteins_found{$key}) {
                 if(!(grep {/^$value$/} @{$proteins_found {$key}})) {
                     $proteins_not_found_parser {$value} = $file_name[5];
                 }
            } else {
                $proteins_not_found_parser {$value} = $file_name[5];
            }
          }
     }
  
     foreach my $key (keys %proteins_found) {
         if (exists $proteins_not_found_ref{$key}) {
             delete $proteins_not_found_ref{$key};
         }
     }
     print "Proteins_not_found_ref = ".scalar(keys %proteins_not_found_ref)."\n";
     &writeGCA (%proteins_not_found_parser);
}  
&writeInputs (%proteins_found);
&writeGCA (%proteins_not_found_ref);

sub align {
    my ($proteins) = $_[0];
    my $proteins_fields = $proteins->getProteins;  
    my %proteins_not_found = ();
    my $query_id = @$proteins_fields[0]->getQueryId;
    my $subject_id = @$proteins_fields[0]->getSubjectId;
    my $bit_score = @$proteins_fields[0]->getBitScore;
    my $subject_cov = @$proteins_fields[0]->getSubjectCov;
    my $hsp_cov = @$proteins_fields[0]->getHspCov;
    my $index = 0;
    my %id = ();
    my @file_name = split(/\_/,$proteins->getFileName);
    my $found = 0;
    for (my $k = 1; $k < scalar (@$proteins_fields); $k++) {   
        if ($query_id eq @$proteins_fields[$k]->getQueryId) {
            if ($bit_score <  @$proteins_fields[$k]->getBitScore) {
                    $bit_score = @$proteins_fields[$k]->getBitScore;
                    $index = $k;
            }        
        } else {
            if (@$proteins_fields[$index]->getIdentity > 50 && @$proteins_fields[$index]->getPositives > 70 
                && @$proteins_fields[$index]->getSubjectCov > 70) {
                 if (exists $proteins_not_found{@$proteins_fields[$index]->getQueryId}) {
                     delete $proteins_not_found{@$proteins_fields[$index]->getQueryId};
                 }                                      
                 if ($file_name[1] eq "REF") {
                     $id {@$proteins_fields[$index]->getQueryId} = @$proteins_fields[$index]->getSubjectId;   
                 } elsif ($file_name[1] eq "PARSER") {
                     push @{$id {@$proteins_fields[$index]->getSubjectId}}, @$proteins_fields[$index]->getQueryId;
                 }            
            } else {                                                                      
                $proteins_not_found {@$proteins_fields[$index]->getQueryId} = $file_name[3];           
            }
            $query_id = @$proteins_fields[$k]->getQueryId;
            $subject_id = @$proteins_fields[$k]->getSubjectId;
            $bit_score = @$proteins_fields[$k]->getBitScore;
            $index = $k;
        }   
    }     
    if (@$proteins_fields[$index]->getIdentity > 50 && @$proteins_fields[$index]->getPositives > 70 
       && @$proteins_fields[$index]->getSubjectCov > 70) {
       if (exists $proteins_not_found{@$proteins_fields[$index]->getQueryId}) {
           delete $proteins_not_found{@$proteins_fields[$index]->getQueryId};
       } 
       if ($file_name[1] eq "REF") {
           $id {@$proteins_fields[$index]->getQueryId} = @$proteins_fields[$index]->getSubjectId;   
       } elsif ($file_name[1] eq "PARSER") {
           push @{$id {@$proteins_fields[$index]->getSubjectId}}, @$proteins_fields[$index]->getQueryId;
       }    
  } else {                        
      $proteins_not_found {@$proteins_fields[$index]->getQueryId} = $file_name[3];
  }
   return (\%id, \%proteins_not_found);
}

sub formatFasta {
    my $sequence = $_[0];
    my $format_sequence = undef;
    for (my $i = 0; $i <= length ($sequence); $i+=80) {
        $format_sequence .= substr ($sequence,$i,80)."\n";
    }
    return $format_sequence;
}  

sub writeInputs {
    my (%proteins_found) = @_;
    my $k = -1;
    my $path = $homologous_subpath."/round_".$number; 
    
    foreach my $key (keys %proteins_found) {
         $k++;         
         my $output = undef;
         open($output, ">", $path."/input_".$k.".faa");
         print $output $fasta{$key};   
         foreach my $value (@{$proteins_found{$key}}) {
             print $output $fasta{$value};
         }
    }
}                                                                                                 
sub writeGCA {
    my (%proteins_not_found) = @_;
    my %proteins = ();
    foreach my $key (keys %proteins_not_found) {     
        #print $fasta{$key}."\n";
        push @{$proteins{$proteins_not_found{$key}}}, $fasta{$key};
    }
    foreach my $key (keys %proteins) {
        my $output = undef;
        open ($output, ">", $homologous_subpath."/round_".$number."/output/GCA_".$key."_protein.faa");
        foreach my $value (@{$proteins{$key}}) {
            print $output $value;
        }
    }
}
