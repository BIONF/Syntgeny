# Syntgeny
Tool to detect if the gene order is conserved over a set of taxa

# Installation
Syntegny requires the targeted ortholog search tool fDOG. For installation and setup of fDOG, you can just look on the [fDOG GitHub page](https://github.com/BIONF/fDOG).

# Usage
Syntgeny requires the same data structure as described in the fDOG manual. 

    python3 syntgeny.py --seed 'maker-scaffold8-snap-gene-67.113' --jobName <Name> --refSpec <reference Species> --gff .<gff_folder> --neighbours <number neighbours upstream and downstream> --out <folder> --searchpath <searchTaxa_dir> --corepath <coreTaxa_dir> --annopath <annotation_dir>

Example Output

    Extracting neighbours of seed gene maker-scaffold8-snap-gene-67.113
    Starting fDOG for seed gene and identified neighbours
    fdogs.run --seqFolder ./tmp/seeds/ --jobName test --refspec DIPRO@73405@230629 --outpath . --corepath ../data/coreTaxa_dir --searchpath ../data/searchTaxa_dir --annopath ../data/annotation_dir
    100%|████████████████████████████████████████████████████████████████████████████████████| 5/5 [04:09<00:00, 49.91s/it]  
    fDOG finished. Parsing fDOG output and extracting ortholog locations
    Creating output table  
    Run time -265.59861493110657`sec

Parameters:

    usage: syntgeny.py [-h] [--version] --seed SEED --refSpec REFSPEC --neighbours NEIGHBOURS --gff GFF --jobName JOBNAME [--fasta FASTA] [--searchTaxa SEARCHTAXA] [--out OUT]
                   [--searchpath SEARCHPATH] [--corepath COREPATH] [--annopath ANNOPATH]

    You are running Syntgeny version 0.0.1.

    options:
      -h, --help            show this help message and exit
      --version             show program's version number and exit

    Required arguments:
      --seed SEED           Seed gene name
      --refSpec REFSPEC     Reference taxon
      --neighbours NEIGHBOURS
                        Number of neighbouring genes up and downstream that should be checked
      --gff GFF             Path to gff file of reference species
      --jobName JOBNAME     Job name

    Optional arguments:
      --fasta FASTA         Path to protein fasta file
      --searchTaxa SEARCHTAXA
                        File containing search species line by line
      --out OUT             Output folder
      --searchpath SEARCHPATH
                        Path for the search taxa directory
      --corepath COREPATH   Path for the core taxa directory
      --annopath ANNOPATH   Path for the pre-calculated feature annotion directory

