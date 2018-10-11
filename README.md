# Proteogenomics

The bioinformatic tool builds proteins sequence databases customized,through processing and analysis of protein sequence data from several strains of the same bacterial species.

For the construction of databases, the tool performs the alignment of protein sequences of bacteria strains. Then, identifies and compares homologous and uniquely annotated proteins in all strains. And finally, reports those sequences in a non-redundant manner, which means, sequences extensively repeated among annotations are reported only once in order to keep the size database under control. Databases also report sequence variations, whether they result from genetic variations or annotation divergences.

## Script design and availability

The tool was designed in PERL and is present as two modules: rand.pl provides the sequence alignment and creates the outputs with unique entries and homologues; create_bd.pl process the homologues output, and create the final database and the log file. In house BLAST installation is required.

![diagrama](https://user-images.githubusercontent.com/8170234/46799804-26eb2580-cd2c-11e8-95bb-4a2da2b573d9.png)

Click [here](https://github.com/karlactm/Proteogenomics.git) to download or clone the scripts.

### Dependencies

Download and install [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi). 

### Get Started

Download / copy FASTA files of interest into input folder 

Execute rand.pl

```
rand.pl /input_folder/ /output_folder/
```

The following output file will be created into output folder. The output files are the result of the bidirectional alignment performed in rand.pl.

Run create_database.pl to create database

```
create_db.pl /output_folder/
```







