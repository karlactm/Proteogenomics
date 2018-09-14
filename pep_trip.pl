#!/usr/bin/env perl
use strict;
use FASTAParser;

my $dir = $ARGV[0];
my $dir_db = $dir."output_pep_trip/database.txt";
my $dir_map = $dir."output_pep_trip/map.txt";
my $dir_id = $dir."output_pep_trip/identification.txt";
my $dir_class = $dir."output_pep_trip/classification.txt";
open(my $db, ">", $dir_db);
open(my $map, ">", $dir_map);
open(my $id, ">", $dir_id);
open(my $class, ">", $dir_class);
print $map "Log Id\tFasta Id\n";

my $value = 0;
$value = sprintf("%05d", $value);
opendir(my $dh, $dir) || die "Can't opendir $dir: $!";
my @round = grep { /^round/ } readdir($dh);
closedir $dh;

for (my $j = 0; $j < scalar (@round); $j++) {
        my $dir_input = $dir.$round[$j]."/";
        opendir(my $input, $dir_input) || die "Can't opendir $dir_input: $!";
        print $dir_input."\n";
        my @inputs = ();
        foreach (readdir ($input)) {
            if ($_ =~/^GCA/) {
                my $gca_protein = $dir_input.$_;
                &gca ($gca_protein);
            } elsif  ($_ =~/^input/) {
                push (@inputs, $dir_input.$_);
            }
        }
        &fasta (@inputs);
        #last;
}

#Processa as proteínas únicas de cada linhagem
sub gca {
    my $gca_protein = $_[0];
    print "> ".$gca_protein."\n";   
    my $parser = FASTAParser->new($gca_protein);
    if (!$parser->parse) {
        print "-> Insira um arquivo no formato FASTA!";
        exit 1; 
    }
    my $proteins = $parser->getProteins;
    for (my $i = 0; $i < scalar (@$proteins); $i++) {
        my $description = @$proteins[$i]->getDescription;
        my $sample = @$proteins[$i]->getSample;
        my $id_protein = @$proteins[$i]->getId;
        my $sequence = @$proteins[$i]->getSequence;
        my $split_protein = &splitProtein ($sequence);
        #@pep = &tripticPeptide($sequence);
        $value++;
        #print "GCA >KT$value \n";
        if ($gca_protein eq "/home/kctmachado/bacteria/dados_blastp/1_random/mycobacterium_tuberculosis/65/output_find_homologous/round_10/GCA_000195955.2_protein.faa") {
            print "KT$value\n";
 
        }
        print $class "KT$value\tReference\n";
        print $id "# Accession Number KT$value\n";
        print $id "# Single entry: $sample [".$id_protein."]\n";
        print $db ">KT$value ".$description."\n".$split_protein;
        print $map "KT$value\t".$id_protein."\tReference\n"; 
    }
}

#Processas as proteínas homólogas entre as linhagens
sub fasta {
    my @inputs = @_;
    my @parsers = ();
    for (my $i = 0; $i < scalar (@inputs) ; $i++) {
        #print "->  $i - $inputs[$i]\n";
        my $parser = FASTAParser->new($inputs[$i]);
        if (!$parser->parse) {
            print "-> Insira um arquivo no formato FASTA!";
            exit 1; 
        }
        push @parsers, $parser;
    }
                                    
    for (my $k = 0; $k < scalar (@parsers); $k++) {
        my %peptides;
        my %classification;
        my %sample;
        my %mutation;
        my %amino;
        my @ref_peptide = ();
        my $ref_sequence = undef;
        my @parser_peptide = ();
        my $parser_sequence = undef;
        my $sequence_size = 0;
        my $index_ref = 0;
        my $index_parser = 0;
        my $proteins = $parsers[$k]->getProteins;
        my @samples_identity = ();
        my %homologous = ();
        #Escolhe a maior proteína do arquivo como referência
        #print $k." ".scalar (@$proteins)."\n";
        for (my $w = 0; $w < scalar (@$proteins); $w++) {
            #print "$k >> ".@$proteins[$w]->getHeader." ".@$proteins[$w]->getSequenceSize." ".$sequence_size."\n";       
            if (@$proteins[$w]->getSequenceSize > $sequence_size) {
                $sequence_size = @$proteins[$w]->getSequenceSize;
                $index_ref = $w;
                $ref_sequence = @$proteins[$index_ref]->getSequence;
                @ref_peptide = &tripticPeptide ($ref_sequence);
                #print join ("\t", @ref_peptide)."\n";
             } elsif (@$proteins[$w]->getSequenceSize == $sequence_size) {
                 if ($ref_sequence =~ /([X|U])/) {
             	     $index_ref = $w;
                     $ref_sequence = @$proteins[$index_ref]->getSequence;
                     @ref_peptide = &tripticPeptide ($ref_sequence);
             	 }
            }
        }
        $homologous {@$proteins[$index_ref]->getId} = scalar (@$proteins);

        #Percorre todas as proteínas
        for (my $w = 0; $w < scalar (@$proteins); $w++) { 
            if ($w == $index_ref) {
                next;
            }
            $parser_sequence = @$proteins[$w]->getSequence;
            @parser_peptide = &tripticPeptide ($parser_sequence);
            my $identity = 1;
            my $last_index = undef;
            #Percorre todos os peptídeos da proteína parse
            for (my $y = 0; $y < scalar (@parser_peptide); $y++) {
                #Percorre todos os pepídeos da proteína referência            
                OUTER: for (my $i = 0; $i < scalar(@ref_peptide); $i++) {
                     #Verifica se os peptídeo são iguais
                     if ($parser_peptide[$y] eq $ref_peptide[$i]) {
                         #Classifica o peptídeo como livre de mutação
                         $peptides {$parser_peptide[$y]} = 1;                    
                         last OUTER;
                     } else {
                         if ($y == 0 && $i == 0) {
                             if ($parser_peptide[$y] =~ /^M/) {
                                 my $peptide = substr ($parser_peptide[$y], 1, length ($parser_peptide[$y])); 
                                 $peptides {$parser_peptide[$y]} = 0;
                                 $classification {$parser_peptide[$y]} = "TSS";
                                 $sample {$parser_peptide[$y]} .= @$proteins[$w]->getSample." [".@$proteins[$w]->getId."]|"; 
                                 $amino {$parser_peptide[$y]} = "-";
                                 $mutation {$parser_peptide[$y]} = "-";
                                 $index_parser = $w;
                                 $identity = 0;
                                 if (length ($peptide) > 6) {
                                     $peptides {$peptide} = 0;
                                     $classification {$peptide} = "TSS";
                                     $sample {$peptide} .= @$proteins[$w]->getSample." [".@$proteins[$w]->getId."]|"; 
                                     $mutation {$peptide} = "-";
                                     $amino {$peptide} = "-";
                                }   
                             } #else {
                               #  $identity = 0;
                             #}
                             #next; 
                         }                                                                                     
                         #Classifica o peptídeo como mutado
                         $peptides {$parser_peptide[$y]} = 0;
                         #Quebra os peptídeos em aminoácidos
                         my @array_ref = split (//, $ref_peptide[$i]);
                         my @array_parser = split (//, $parser_peptide[$y]);  

                         if (scalar (@array_ref) == scalar (@array_parser)) {  
                             my @index_aminoacids = split(",", &comparePeptide(\@array_ref, \@array_parser));
                             if ($index_aminoacids [$#index_aminoacids] eq " ") {
                                 pop (@index_aminoacids);
                             } 
                             if (scalar (@index_aminoacids) <= 0.2 * length ($ref_peptide[$i])) {                             	
                                foreach (@index_aminoacids) {                                                   
                                     if (&mutationInvalid ($array_ref[$_],$array_parser[$_])) {
                                          #last OUTER;
                                          #$peptides {$parser_peptide[$y]} = 1; 
                                          $mutation {$parser_peptide[$y]} = $ref_peptide[$i];
                                          next;
                                      } else {
                                         #$peptides {$parser_peptide[$y]} = 0;
                                         $identity = 0;
                                         my $aminoacids .= $array_ref[$_]."->".$array_parser[$_]; 
                                         $amino {$parser_peptide[$y]} = $aminoacids;
                                         $mutation {$parser_peptide[$y]} = $ref_peptide[$i];
                                         $sample {$parser_peptide[$y]} .= @$proteins[$w]->getSample." [".@$proteins[$w]->getId."]|";
                                         #print "->".$sample {$parser_peptide[$y]}."\n";
                                         $classification {$parser_peptide[$y]} = "SAP";
                                         $index_parser = $w;
                                         last OUTER;
                                     }    
                                }
                            }
                        } elsif (scalar (@array_ref) > scalar (@array_parser))  {
                              my @index_aminoacids = split(",", &comparePeptide(\@array_ref, \@array_parser));
                              if ($index_aminoacids [$#index_aminoacids] eq "-1") {
                                  pop (@index_aminoacids);
                                  #print "=>> ".join ("", @array_parser)." y = $y i = $i\n";
                                  if ((scalar (@index_aminoacids)) <= (int (0.2 * length ($parser_peptide[$y])))) {
                                     foreach (@index_aminoacids) {                                                   
                                         if (&mutationInvalid ($array_ref[$_],$array_parser[$_])) {
                                            next; 
                                         } else {                                 	
			                                 $identity = 0;
                                             my $aminoacids .= $array_ref[$_]."->".$array_parser[$_];
                                             $amino {$parser_peptide[$y]} = $aminoacids;
                                             $mutation {$parser_peptide[$y]} = $ref_peptide[$i];
                                             $sample {$parser_peptide[$y]} .= @$proteins[$w]->getSample." [".@$proteins[$w]->getId."]|";
                                             #print "->".$sample {$parser_peptide[$y]}."\n";
                                             $classification {$parser_peptide[$y]} = "SAP R|K";
                                             $index_parser = $w;
                                             $last_index = ($index_aminoacids[$#index_aminoacids] + 1);
                                             last OUTER;
                                        }     
                                    }

                                  } else { #Processa a segunda parte do peptídeo que foi quebrado
                                      if (defined $last_index) {
                                          my @div_ref = ();
                                          my $new_peptide = undef;
                                          for (my $x = $last_index; $x <= $#array_ref; $x++) {
                                              push (@div_ref, $array_ref[$x]);
                                              $new_peptide .= $array_ref[$x];
                                          }
                                          if ($ref_peptide [$i] eq "HFYSQAVEEXNHAMMLVQHLLDR") {
                                               print $ref_peptide [$i]." ".$parser_peptide [$y]." ".$new_peptide."\n";
                                          }
                                          if (scalar (@div_ref) == scalar (@array_parser)) {
                                              my @index_aminoacids = split(",", &comparePeptide(\@div_ref, \@array_parser));
                                              if ((scalar (@index_aminoacids)) <= (int (0.2 * length ($parser_peptide[$y])))) {
                                              	  #foreach (@index_aminoacids) {  
                                                  #    if (&mutationInvalid ($array_ref[$_],$array_parser[$_])) {
                                                  #        next;
                                                  #    } else { 
                                                          $identity = 0;
                                                          my $aminoacids .= $array_ref[$_]."->".$array_parser[$_];
                                                          $amino {$parser_peptide[$y]} = $aminoacids;
                                                          $mutation {$parser_peptide[$y]} = $ref_peptide[$i];
                                                          $sample {$parser_peptide[$y]} .= @$proteins[$w]->getSample." [".@$proteins[$w]->getId."]|";
                                                          #print "->".$sample {$parser_peptide[$y]}."\n";
                                                          $classification {$parser_peptide[$y]} = "SAP R|K";
                                                          $index_parser = $w;
                                                          $last_index = undef;
                                                          last OUTER;
                                                      #}
                                                   #}
                                             }
                                        }
                                    } 
                                }  
                              } #else {
                                  #Inserção e deleção
                              #}                            
                        } elsif (scalar (@array_ref) < scalar (@array_parser))  {
                              my @index_aminoacids = split(",", &comparePeptide(\@array_ref, \@array_parser));
                              if ($index_aminoacids [$#index_aminoacids] eq "-1") {
                                  pop (@index_aminoacids);
                                  if (scalar (@index_aminoacids) <= 0.2 * length ($ref_peptide[$i])) {
                                       foreach (@index_aminoacids) {                                                   
                                          if (&mutationInvalid ($array_ref[$_],$array_parser[$_])) {
                                              next;
                                          } else {
                                              $identity = 0;
                                              my $aminoacids .= $array_ref[$_]."->".$array_parser[$_];
                                              $amino {$parser_peptide[$y]} = $aminoacids;
                                              $mutation {$parser_peptide[$y]} = $ref_peptide[$i];
                                              $sample {$parser_peptide[$y]} .= @$proteins[$w]->getSample." [".@$proteins[$w]->getId."]|";
                                              #print "->".$sample {$parser_peptide[$y]}."\n";
                                              $classification {$parser_peptide[$y]} = "SAP R|K";
                                              $index_parser = $w;
                                              last OUTER;
                                        }
                                      }
                                  }
                              } #else { 
                                  #Inserção e Deleção
                              #}
                      
                     }
                  } #Fim do else
               } #Fim do for que percorre ref_peptide
               
                 
                if (($peptides {$parser_peptide[$y]} == 0) && (($mutation {$parser_peptide[$y]} eq undef) || ($mutation {$parser_peptide[$y]} eq "NOT FOUND"))) {
                      #print "y = $y i =  ".$parser_peptide[$y]."\t".$peptides{$parser_peptide[$y]}."\n";
                      $identity = 0;
                      $mutation {$parser_peptide[$y]} = "NOT FOUND";
                      $amino {$parser_peptide[$y]} = "-";
                      $sample {$parser_peptide[$y]} .= @$proteins[$w]->getSample." [".@$proteins[$w]->getId."]|";
                      $classification {$parser_peptide[$y]} = "NC";
                      $peptides {$parser_peptide[$y]} = 0;
                }
           } #Fim do for que percorre parser_peptide
           if ($identity) {
               my $sample_id = @$proteins[$w]->getSample." [".@$proteins[$w]->getId."]";
               push (@samples_identity, $sample_id); 
           } 
           #last;
        } #Fim do for que precorre proteins

        my $artificial_peptide = undef;
        my $c_terminal = undef;
        foreach (keys %peptides) {
        	my $boolean = $peptides {$_};
            if (!$boolean) {
        	    if ($_ =~ /([^RK]+(R|K))/gci) {
        		    $artificial_peptide = $_;
        	    } else {
                    $c_terminal = $_;
        	    }
        	}
        } 
        if (defined ($c_terminal)){
            $artificial_peptide .= $c_terminal; 
        }  
        $value++;	
        #print ">KT$value ".$artificial_peptide."\n";
        my $split_protein = &splitProtein (@$proteins[$index_ref]->getSequence);
        print $db ">KT$value ".@$proteins[$index_ref]->getDescription."\n".$split_protein; 
        #print $k." >KT$value ".@$proteins[$index_ref]->getDescription."\n".$split_protein."\n";
        my $s = undef;
        for (my $h = 0; $h < scalar (@samples_identity); $h++){
            $s .= $samples_identity[$h]."|";
        }
        if (defined $artificial_peptide) {
            print $class "KT$value\tReference\t";
            print $id "# Accession Number KT$value\n";
            print $id "# Reference entry: ".@$proteins[$index_ref]->getSample." [".@$proteins[$index_ref]->getId."]\n";
            print $id "# Number of strains in homologue file = ".$homologous {@$proteins[$index_ref]->getId}."\n";
            print $id "# Number of strains with 100% identity = ".scalar (@samples_identity)."\n";
            print $map "KT$value\t".@$proteins[$index_ref]->getId."\tReference\n";
            if (defined $s) {
                #$s = &splitSamples($s);
                print $id "# Fields: entries with 100% identity\n$s\n";
            }
            $value++;
            print $class "KT$value\tMutated\n";
            if (length ($artificial_peptide) > 70) {
                $artificial_peptide = &splitProtein ($artificial_peptide);
            } else {
            	$artificial_peptide .= "\n";
            }
            #print "PEP DIF >KT$value\n";
            print $db ">KT$value ".@$proteins[$index_parser]->getDescription."\n".$artificial_peptide;
            #print $class "KT$value\tMutated\n";
            print $id "# Accession Number KT$value\n";
            print $id "# Artificial entry\n";
            print $id "# Fields: modification, peptide, peptide mutated, aminoacid mutation, strain\n";
            print $map "KT$value\t".@$proteins[$index_ref]->getId."\tMutated\n";
            foreach (keys %classification) {
                #print "-".$sample{$_}."\n";
                #my $split_sample = &splitSamples ($sample{$_});
                print $id $classification{$_}."\t".$mutation {$_}."\t".$_."\t".$amino {$_}."\t".$sample {$_}."\n";  
            }
        } else {
            if (defined $s) {
                #$s = &splitSamples ($s);
                print $class "KT$value\tReference\n";
                print $id "# Accession Number KT$value\n";
                print $id "# Reference entry\n";
                print $id "# Number of strains in homologue file = ".$homologous {@$proteins[$index_ref]->getId}."\n";
                print $id "# Number of strains with 100% identity = ".scalar(@samples_identity)."\n";
                print $id "# Fields: Entries with 100% identity\n".$s."\n";;
                print $map "KT$value\t".@$proteins[$index_ref]->getId."\tReference\n";
            } else {
              print $class "KT$value\tReference\n";
              print $id "# Accession Number KT$value\n";
              print $id "# Reference entry\n";
              print $id "# Number of strains in homologue file = ".$homologous {@$proteins[$index_ref]->getId}."\n";
              print $id "# Number of strains with 100% identity = ".scalar (@samples_identity)."\n";
              print $map "KT$value\t".@$proteins[$index_ref]->getId."\tReference\n";
            }
        } 

    }
}

sub tripticPeptide {
    my $sequence = $_[0];
    my @peptide = ();
    $sequence =~ s/\\|\R//g;	
    while ($sequence =~ /([^RK]+(R|K|$))/gci) {
        if (length ($1) >= 7 && length ($1) <= 35) { 
            #print "$1\n";
            push (@peptide, $1); 
        }
    }
    return @peptide;
}

sub mutationInvalid {
   my ($aminoacid_1, $aminoacid_2) = @_;
   if ($aminoacid_1 =~ /(I|L)/ && $aminoacid_2 =~/(I|L)/) {
     return 1;
   } elsif ($aminoacid_1 =~ /(X|U)/ || $aminoacid_2 =~ /(X|U)/ ) {
     return 1;
   }
}


sub splitProtein {
    my ($protein) = $_[0];
    my $split_protein = "";
    for (my $i = 0; $i < length ($protein); $i += 80) {
        $split_protein .= substr ($protein, $i, 80)."\n";
     }
     return $split_protein;
}

sub comparePeptide {
    my ($peptide_1, $peptide_2) = @_;
    my $index_aminoacids = undef;
    for (my $i = 0; $i < scalar (@$peptide_1); $i++) {
        if (@$peptide_1[$i] ne @$peptide_2[$i]) {
            $index_aminoacids .= $i.",";
            if (@$peptide_1[$i] =~ /(R|K)/ || @$peptide_2[$i] =~/(R|K)/) {
                #print "-> ".@$peptide_1[$i]." ".@$peptide_2[$i]."\n";
                #Se mutação de gerar um novo peptídeo tríptico setar $index_aminoacids para -1
                $index_aminoacids .= "-1"; 
                return $index_aminoacids;
            }               
        }
    }
    return $index_aminoacids;
}

sub splitSamples {
    my @samples = split (//,$_[0]);
    my $counter = 0;
    my $size = 0;
    my $split_sample = undef;
    for (my $k = 0; $k < scalar (@samples); $k++) {
        $counter++;
        $split_sample .= $samples[$k];
        if ($samples[$k] eq "|") {
            if ($counter > 70) {
                #$size = (length ($split_sample) - 1);
                #$split_sample = substr ($split_sample , 0, $size);
                #print $split_sample."\n";
                $split_sample .= "\n";
                $counter = 0;
             } else {
                 if ($k == (scalar (@samples) - 1)) {
                      $split_sample .= "\n";
                  }
            }
        }
        #$counter++;
    }
    #$size = (length ($split_sample) - 1);
    #$split_sample = substr ($split_sample, 0, $size);
    return $split_sample;
}

print "pronto\n";
