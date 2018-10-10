# Proteogenomics

This computational strategy builds proteins sequence databases customized,through processing and analysis of protein sequence data from several strains of the same bacterial species.

For the construction of databases, the approach performs the alignment of protein sequences of bacteria strains. Then, identifies and compares homologous and uniquely annotated proteins in all strains. And finally, reports those sequences in a non-redundant manner, which means, sequences extensively repeated among annotations are reported only once in order to keep the size database under control. Databases also report sequence variations, whether they result from genetic variations or annotation divergences.

## Get Started

Click [here](https://github.com/karlactm/Proteogenomics.git) to download or clone the source code.

### Dependencies

Download and install [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi). 

### Next steps

Download / copy FASTA files of interest into input folder 

Execute rand.pl

```
rand.pl /input_folder/ /output_folder/
```

The following output file will be created into output folder. The output files are the result of the bidirectional alignment performed in rand.pl.

Run create_database.pl to create database

```
create_database.pl /output_folder/
```







