# Proteogenomics

The bioinformatic tool builds proteins sequence databases customized, through processing and analysis of protein sequence data from several strains of the same bacterial species.

The tool defines unique annotated proteins as well as protein homologues across strains, adding all unique sequence information in the final database on a non-redundant manner.

## Script design and download

The tool was designed in PERL and is present as two modules: 
1. **all_fasta.pl** provides the sequence alignment and creates outputs with unique entries and homologues using all fasta files in folder. **rand.pl** provides the sequence alignment and creates the outputs with unique entries and homologues using random fasta inputs, to reproduce data analysis from Machado et al., bioRx 2018.
2. **pep_trip.pl** process the homologues output, and creates a database with all fasta files in folder. **create_bd.pl** process the homologues output, and create different sized databases. 

In house [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) installation is required.

![diagrama](https://user-images.githubusercontent.com/8170234/46799804-26eb2580-cd2c-11e8-95bb-4a2da2b573d9.png)

Click [here](https://github.com/karlactm/Proteogenomics.git) to download or clone the scripts.

### Get Started

Download / copy FASTA files of interest into input folder 

Execute **all_fasta.pl** to create a database with all fasta files in folder. Use **rand.pl** to create different sized databases using random fasta inputs.

```
all_fasta.pl /input_folder/ /output_folder/
```

The following homologues output will be created into output folder. 

Execute **pep_trip.pl** to create a database with all fasta files in folder. Use **create_database.pl** to create different sized databases. 

```
pep_trip.pl /output_folder/output_find_homologous/
```







