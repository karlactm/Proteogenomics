# Proteogenomics

The bioinformatic tool builds proteins sequence databases customized,through processing and analysis of protein sequence data from several strains of the same bacterial species.

The tool defines unique annotated proteins as well as protein homologues across strains, adding all unique sequence information in the final database on a non-redundant manner.

## Script design and availability

The tool was designed in PERL and is present as two modules: rand.pl provides the sequence alignment and creates the outputs with unique entries and homologues; create_bd.pl process the homologues output, and create the final database and the log file. In house [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) installation is required.

![diagrama](https://user-images.githubusercontent.com/8170234/46799804-26eb2580-cd2c-11e8-95bb-4a2da2b573d9.png)

Click [here](https://github.com/karlactm/Proteogenomics.git) to download or clone the scripts.

### Get Started

Download / copy FASTA files of interest into input folder 

Execute rand.pl

```
rand.pl /input_folder/ /output_folder/
```

The following homologues output will be created into output folder. 

Run create_database.pl to create database

```
create_db.pl /output_folder/
```







