# KEGGChem
## Summary
KEGGChem is a Python program for the batch retrieval of KEGG compound details using the KEGG API. The primary usage is to retrieve KEGG compounds using an inputted list of KEGG orthologue, module, or reaction codes, as well as


## Dependencies
KEGGChem requires the following Python libraries:
+ collections
+ requests
+ PubChemPy


## How to use:
KEGGChem is accessed by running the Python script using two required arguments, mode and input file. 
`python KEGGChem.py --mode <mode> --file <input_file>`
Optional arguments include setting the output directory with `--outdir`, downloading compound structural data `--structure`, or changing the downloads directory `--download`.
 

Help can be accessed py calling `python KEGGChem.py --help` or `python KEGGChem.py -h`. 

```
usage: KEGGChem [-h] -m {ko,module,compound,reaction,mdata} -f  [-o] [-d] [-s] [--pubchem] [--sdf] [--homolog] [--homolog-organism] [--quiet | --verbose]

KEGGChem is a simple web scraper of the KEGG API, that extracts information from the KEGG compound database using either KEGG orthologues or modules as input.

options:
  -h, --help            show this help message and exit
  -m {ko,module,compound,reaction,mdata}, --mode {ko,module,compound,reaction,mdata}
                        The mode to run KEGGChem in. Select either 'ko' or 'module' to search for compounds using either ko or module codes as input. Use 'compound' or 'reaction' mode to download related data from listed data entries.
  -f , --file           The input file
  -o , --out            The output directory to be used. If this argument is not used, then a default directory will be created using the input filename + '_output'.
  -d , --download       Download directory. If none is given, a new one will be made if not present already. Default = KEGGChem_downloads
  -s, --structure       Download structural data for compounds, such as weight and mass, and .mol and .kcf files. Default = False
  --pubchem             Search PubChem database using SID codes retrieved from compounds. Default = False
  --sdf                 Download sdf file from PubChem. Requires --pubchem argument. Default = False
  --homolog             Find genes on the KEGG database with the same KO codes. Requires --mode=ko. Default = False
  --homolog-organism    Restrict the genes to be downloaded with the --homolog command. This can be done with one of the following options: -Input a KEGG organism code (e.g., hsa, ggo, ptr). -Input a file containing a list of KEGG
                        organism codes. -Input a search keyword, which corresponds to a premade list of organism codes for all KEGG organisms in a given clade. E.g., animals, mammals, bacteria.
  --quiet               Run program quietly and reduce data printed to screen. Default = False
  --verbose             Print extra program details to screen. Default = False

Thank you for using KEGGChem. More details can be found at https://github.com/Lachiemckbioinfo/KEGGChem. This program is neither endorsed by or affiliated with KEGG. Please cite the relevant KEGG literature when using this program.
```
## Command Arguments

## Mode
KEGGChem has several different modes of function, which determine how it behaves. The modes currently available include: ko, module, compound, and reaction. The modes ko and module are used to retrieve KEGG compound codes, while the modes compound and reaction retrieve additional metadata for those codes.

### ko
KO mode retrieves KEGG compound codes using an input of KEGG orthologue (KO) codes. This is done by retrieving the KEGG reaction codes from the KO API page, then retrieving the KEGG compound codes from each linked module and reactionpage.

KO results will be outputted into the supplied output directory containing two subdirectories, titled "Results" and "Summaries". The results directory contains the subdirectory "Individual results", which stores results for individual orthologues, as well as a number of text files. The files "compounds.txt", "glycans.txt", and "reactions.txt" contain all compounds, glycans and reactions retrieved by the program through all methods. Additional results files include results such as "compounds_through_reactions.txt" and "glycans_from_modules", which give resultsing compounds retrieved from KO entries either throught their linked module or reaction entries. Finally, the ko results directory also contains the files "all_module_results.txt" and "all_reaction_results.txt" which contain all the respective compounds retrieved for each method for each KEGG orthologue. The "Summaries" directory will contain a run summary, as well as summaries of the number and frequency of returns for compounds, glycans, and reactions.

### module
Module mode retrieves compound codes directly from the KEGG compound pages. Note that the list of compound on module pages may be shorter than the list you will get from retrieving the compound codes from each reaction page associated with said module.

Running KEGGChem on module mode will output the results

The module mode will output the results into a single summary file and a subdirectory titled "Results". The summary file will summaries the number of modules given and the number of unique compounds retrieved, as well as the total number of times each compound was retrieved. The "Results" subdirectory will contain files for each module showing the compounds retrieved in the format of <module>: <compound>, <compound>, <compound>, and so on. There will also be a file containing the results for all modules together in a single file called all_results.txt

### compound
Compound mode retrieves data related to KEGG compounds from the KEGG database. This includes the molecular formula, molecular weight, and exact mass. This can be further expanded using the ```--structure``` argument, which downloads .mol and .kcf files as well. This mode also be expanded with PubChemPy integratration using the ```--pubchem``` argument, which searches for SID and CID codes and associated data, as well as the ```--sdf``` command, which downloads sdf files from PubChem.

### reaction
Reaction mode retrieves data related to KEGG reaction entries, including definitions, molecular equations, listed KO codes, reaction pathways and reaction classes.

## Other Arguments
 
### Input
The input file command (`-f` or `--file`) is a required option. All three modes use a raw text files as input, in the form of a list of KEGG orthologue, compound, or module codes, with one code per line. There is no file extension requirement.

### Output directory
The output directory, or out option (`-o` or `--out`) is an optional argument that sets the output directory. When run, KEGGChem will output the results in the given directory, and will create the directory if necessary. If no outdir is specified, KEGGChem will create a directory using the name of the input file plus "_output". 

### Structure
 Compound .mol and .kcf files can be downloaded from the KEGG database using the structure argument (`-s` or `--structure`). These are saved in the KEGG_downloads directory, and local versions will be saved in the results directory. This argument requires KEGGChem to be run on module mode.
### Pubchem
Using the ```--pubchem``` arguement enables KEGGChem to retrieve data from Pubchem using PubChemPy, including SIDs (substance IDs), CIDs (compound IDs), and SMILES. This argument is required in order to use the ```--sdf``` argument. The SID search will search each compound for SID codes, and each SID for CID codes. CIDs are then searched for molecular formulas, and both canonical and isomeric SMILES. These results are saved in the "pubchem_substances.txt" and "pubchem_compounds.txt" files, respectively.

### SDF
SDF files can be downloaded from PubChem using the ```--sdf``` argument. This requires the ```--pubchem``` argument, which requires KEGGChem to be on compound mode. SDF files are stored in the SDF subdirectory of the Results directory.

### Homolog
KEGG genes can be downloaded from KEGG using the homolog argument (`--homolog`) when using KO mode. These are saved in the Genes subdirectory of the Results directory. This can be refined using the **homolog-organism** argument.

### Homolog-organism
Gene homolog searches of the KEGG database can be further refined using the homolog-organism argument (`--homolog-organism`). There are three ways to use this argument:
1. Input an organism code from the KEGG database. E.g., hsa, ggo, ptr. The input command would be `--homolog-organism hsa`.
2. Input a file containing a list of organism codes from KEGG. For example: `--homolog-organism organism_list.txt`
3. Input a keyword corresponding to a major category of organisms, as described by the KEGG organism list. These correspond to major categories of organisms, as well as important subgroups, including prokaryotes, eukaryotes, mammals, birds, bacteria, insects. These are not case sensitive. A list of keywords is below.

**Organism keywords** (slashes indicate alternate spellings)
+ Prokaryotes/Prokaryota
+ Eukaryotes/Eukaryota
+ Animals
+ Plants
+ Bacteria
+ Fungi
+ Protists
+ Archaea
+ Cyanobacteria/Cyanobacteriota
+ Vertebrates
+ Mammals
+ Birds
+ Fish/Fishes
+ Reptiles
+ Amphibians
+ Insects
+ Molluscs/Mollusks
+ Crustaceans
+ Chelicerates
+ Green_algae
+ Red_algae
+ Cnidaria
+ Eudicots
+ Monocots
+ Tunicates
+ Ameobozoa
+ Alveolates
+ Acidobacteriota
+ Escherichia

## Legal
KEGGChem is published under a MIT license (details in license.txt). 
Neither the author nor KEGGChem is not associated with, nor endorsed by the KEGG team in any way. Please abide by the copyright and licensing terms of the KEGG program and associated programs as described by their products. 



## Author details
Lachlan McKinnie

ORCiD: [0000-0002-4996-5941](https://orcid.org/0000-0002-4996-5941)

This program was written as part of a PhD project at the University of the Sunshine Coast, which was supervised by Dr Min Zhao, Professor Scott Cummins, and Dr Sankar Subramanian.


