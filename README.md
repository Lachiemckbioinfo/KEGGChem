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
usage: KEGGChem [-h] -m {ko,module,compound,reaction,mdata} -f  [-o] [-d] [-s] [--pubchem] [--sdf] [--quiet | --verbose]

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
  --quiet               Run program quietly and reduce data printed to screen. Default = False
  --verbose             Print extra program details to screen. Default = False

Thank you for using KEGGChem. More details can be found at https://github.com/Lachiemckbioinfo/KEGGChem. This program is neither endorsed by or affiliated with KEGG. Please cite the relevant KEGG literature when using this program.
```
## Command Arguments

## Mode
Mode is a positional argument with three options: ko, module, and compound. This is given by inputting the mode right after calling KEGGChem:
`python KEGGChem.py <mode>`

### ko
KO mode retrieves KEGG compound codes using an input of KEGG orthologue (KO) codes. This is done by retrieving the KEGG reaction codes from the KO API page, then retrieving the KEGG compound codes from each reaction page. This mode will also function if module codes are input instead of KO codes, using the same function.

KO results will be outputted into the supplied output directory, and will contain a run_summary.txt file, as well as two subdirectorries, titled "Results" and "Summaries". The Results directory will contain the raw results for each individual KEGG orthologue as well as a file containing all the results, while the Summaries directory will contain summary files showing the counts of different reactions and compounds.

### module
Module mode retrieves compound codes directly from the KEGG compound pages. Note that the list of compound on module pages may be shorter than the list you will get from retrieving the compound codes from each reaction page associated with said module. If you wish to retrieve the compound via the reaction pages associated with each compound given, use the ko mode, but provide a list of modules instead of KEGG orthologues.

The module mode will output the results into a single summary file and a subdirectory titled "Results". The summary file will summaries the number of modules given and the number of unique compounds retrieved, as well as the total number of times each compound was retrieved. The "Results" subdirectory will contain files for each module showing the compounds retrieved in the format of <module>: <compound>, <compound>, <compound>, and so on. There will also be a file containing the results for all modules together in a single file called all_results.txt

### compound
Compound mode is a convenience script function that retrieves the formula for each given KEGG compound. 

Results are outputted in a .txt file in the format of KEGG compound:compound formula. Compounds without formulas will have their formula result be NULL.

## Other Arguments
 
### Input
The input file command (`-f` or `--file`) is a required option. All three modes use a raw text files as input, in the form of a list of KEGG orthologue, compound, or module codes, with one code per line. There is no file extension requirement.

### Output directory
The output directory, or out option (`-o` or `--out`) is an optional argument that sets the output directory. When run, KEGGChem will output the results in the given directory, and will create the directory if necessary. If no outdir is specified, KEGGChem will create a directory using the name of the input file plus "_output". 

### Structure
 Compound .mol and .kcf files can be downloaded from the KEGG database using the structure argument (`-s` or `--structure`). These are saved in the KEGG_downloads directory, and local versions will be saved in the results directory. This command cannot be used with the compound or reaction modes. 


## Legal
The author of these scripts is not associated with, nor endorsed by the KEGG team in any way.



## Author details
Lachlan McKinnie

ORCiD: [0000-0002-4996-5941](https://orcid.org/0000-0002-4996-5941)


