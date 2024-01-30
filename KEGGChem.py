#!/usr/bin/env python3
#version 0.9.2

import requests
import re
import sys
import os
from collections import Counter
import argparse
import statistics as stat
import datetime
import pubchempy as pcp
import csv
import signal
#import logging


github = "https://github.com/Lachiemckbioinfo/KEGGChem"
citation = "XXXXX"
#KEGGcitation = (f"KEGG citation goes here")
appname = 'KEGGChem'
version = 'KEGGChem 0.9.2'



#Initiate argparse
parser = argparse.ArgumentParser(
    prog = f'{appname}',
    description = f'{appname} is a simple web scraper of the KEGG API, that extracts information from the KEGG compound database using either KEGG orthologues or modules as input.',
    epilog = f'Thank you for using {appname}.\nMore details can be found at {github}.\nThis program is neither endorsed by or affiliated with KEGG. Please cite the relevant KEGG literature when using this program.'
)

parser.add_argument("-m", "--mode", 
                    choices=["ko", "module", "compound", "reaction", "mdata"],
                    required = True,
                    #help = "The mode to run in. Select either 'ko', 'module', to search for compounds using the relevant data. Use 'compound', or 'reaction' to download data for items from those databases. Use 'download'"
                    help = f"""
                    The mode to run {appname} in. Select either 'ko' or 'module' to search for compounds using either ko or module codes as input. 
                    Use 'compound' or 'reaction' mode to download related data from listed data entries.
                    """)
                    

parser.add_argument("-i", "--input",
                    required = True,
                    help = "The search input. Either select a file, or input a searchterm or KEGG entry code.",
                    #type = argparse.FileType('r'),
                    #Blank metavar argument to make help output read a bit better
                    metavar='',
                    dest = 'infile')

parser.add_argument("-o", "--out",
                    required = False,
                    help = """
                    The output directory to be used. If this argument is not used, then a default directory will be created using
                    the input filename + '_output'.
                    """,
                    dest = 'outdir',
                    metavar = '')

parser.add_argument("-d", "--download",
                    required = False,
                    help = f"Download directory. If none is given, a new one will be made if not present already. Default = {appname}_downloads",
                    metavar = '')

parser.add_argument("-s", "--structure",
                    required = False,
                    help = "Download structural data for compounds, such as weight and mass, and .mol and .kcf files. Default = False",
                    action = "store_true")

#PubChem search argument using PubChemyPy
parser.add_argument("--pubchem",
                    required = False,
                    help = "Search PubChem database using SID codes retrieved from compounds. Default = False",
                    action = "store_true")

parser.add_argument("--sdf",
                    required = False,
                    help = "Download sdf file from PubChem. Requires --pubchem argument. Default = False",
                    action = "store_true")


#Retrieve genes
parser.add_argument("--homolog",
                    required = False,
                    help = "Find genes on the KEGG database with the same KO codes. Requires --mode=ko. Default = False",
                    action = "store_true")


parser.add_argument("--homolog-organism",
                    required = False,
                    help = """
                    Restrict the genes to be downloaded with the --homolog command. This can be done with one of the following options:
                    -Input a KEGG organism code (e.g., hsa, ggo, ptr).
                    -Input a file containing a list of KEGG organism codes.
                    -Input a search keyword, which corresponds to a premade list of organism codes for all KEGG organisms in a given clade. E.g., animals, mammals, bacteria.
                    """,
                    #type = argparse.FileType('r'),
                    metavar = '',
                    dest = "homologorg")

parser.add_argument("-w", "--overwrite",
                    required = False,
                    help = "Download and overwrite stored files. Default = False.",
                    action = "store_true")





#Mutually exlclusive arguments --quite and --verbose, for printing less or more output, respectively
group = parser.add_mutually_exclusive_group()

group.add_argument("--quiet",
                    required = False,
                    help = "Run program quietly and reduce data printed to screen. Default = False",
                    action = "store_true")

group.add_argument("--verbose",
                   required = False,
                   help = "Print extra program details to screen. Default = False",
                   action = "store_true")



#Parse argparse parser parse
args = parser.parse_args()
mode = args.mode
structure = args.structure
infile = args.infile
quiet = args.quiet
verbose = args.verbose
pubchemarg = args.pubchem
sdfarg = args.sdf
homolog = args.homolog
homologorg = args.homologorg
overwrite = args.overwrite
download = args.download


linebreak = f"\n{'-'*100}\n"


#Set start and finish times with datetime
starttime = datetime.datetime.now()
timestamp = starttime.strftime("%b-%d %H:%M:%S")
filedatestamp = starttime.strftime("%b-%d")
#Function to set endtime and calculate runtime
def timetaken():
    endtime = datetime.datetime.now()
    timestamp_end = endtime.strftime("%b-%d %H:%M:%S")
    runtime = endtime - starttime
    return runtime, timestamp_end


   



#Set output directory. If none given, use input file as output directory name
def set_outdir():
    #Function to check if outdir of same name exists. If so, add a date, and if still exists, add a number (starting at one, goes up with repeats)
    def check_outdir(outdir):
        if os.path.exists(outdir) == True:
            #Set filename of default name + date
            date_outdir = f"{outdir}_{filedatestamp}"
            if os.path.exists(date_outdir) == False:
                outdir = date_outdir
            else:
                outnum = 1
                while os.path.exists(f"{date_outdir}_{outnum}") == True:
                    outnum += 1
                outdir = f"{date_outdir}_{outnum}"
        return outdir

    #If args.outdir is a value, use that value as outdir
    #If directory does not exist, proceed as normal
    if args.outdir is not None:
        if os.path.exists(args.outdir) == False:
            outdir = str(args.outdir)
        #If outdir does exist, create another directory inside it with the same name
        else:
            outdir_internal = os.path.join(str(args.outdir), str(args.outdir))
            outdir = check_outdir(outdir_internal)
            #raise SystemError("Error: Result directory already exists")
    else:
        outdir = str(infile + "_output")
        outdir = check_outdir(outdir)
    return outdir
outdir = set_outdir()
            

        

#Set downloads directory
if args.download is not None:
    dir_download = str(args.download)
else:
    dir_download = "keggchem_downloads"
dirs = ["KEGG_entries/compounds", "KEGG_entries/orthologues", "KEGG_entries/reactions", "KEGG_entries/genes", "KEGG_lists",
        "KEGG_entries/modules", "Structure_files/Downloaded", "Structure_files/nullfiles", "SDFfiles", "KEGG_links/genes"]
dir_download = os.path.abspath(dir_download)
if verbose == True:
        print(f"{linebreak}Checking if download subdirectories exist{linebreak}")
#Create each subdirectory in download directory (ok is exists already)
for item in dirs:
    path = os.path.join(dir_download, item)
    os.makedirs(path, exist_ok = True)

#Set separator (aka ignore everything after separator to avoid false positives)
if mode == "ko" or mode == "module":
    separator = "GENES"
 

input = []
input_invalid = {}

#Extract input codes and appent to input list
def openfile(x):
    if os.path.exists(infile) == True:
        with open(infile, "r") as file:
            line_number = 0
            while (line := file.readline().strip()):
                line = line.upper()
                line_number += 1
                if "KO:" in line:
                    line = line.replace("KO:", "")
                #If search == [LETTER] or RC + 5 digits, proceed
                if mode != "download":
                    if re.search(rf"{x}\d{{5}}\b", line):
                    #if re.search(r"\b([MKTCGRNHD]|RC)\d{5}\b", line):
                        input.append(line)
                    else:
                        input_invalid[line_number] = line
                        if  verbose == True:
                            print(f"Invalid search term: {line}\n")
                else:
                    if re.search(r"\b([MKTCGRNHD]|RC)\d{5}\b", line):
                        input.append(line)
                    else:
                        input_invalid[line_number] = line
                        if verbose == True:
                            print(f"Invalid search term {line}")
        input_mode = "file"
    else:
        #Search KEGG list 
        input_dict = {}
        def request_input(mode):
            if mode == "mdata":
                modeterm = "module"
            else:
                modeterm = mode
            url = f"https://rest.kegg.jp/list/{modeterm}"
            input_list_file = os.path.join(dir_download, "KEGG_lists", f"{modeterm}.txt")
            return url, input_list_file
        #Return url and input_list_file
        url, input_list_file = request_input(mode)
        #Process input
        if os.path.exists(input_list_file) == False or overwrite == True:
            req = requests.get(url).text
            #Write request to file
            with open(input_list_file, "w") as filehandle:
                filehandle.write(req)
            inputdata = req.splitlines()
        else:
            with open (input_list_file, "r") as filehandle:
                inputdata = filehandle.readlines()
        #Process ko
        for item in inputdata:
            data = item.split("\t")
            input_code = data[0]
            input_searchterm = data[1].lower()
            input_dict[input_code] = input_searchterm
        #Search ko_dict for KO number or search term
        for input_code, input_searchterm in input_dict.items():
            infile_list = infile.split(",")
            #Lowercase infile_list
            infile_list = [item.lower().strip() for item in infile_list]
            if input_code.lower() in infile_list:
                input.append(input_code)
                if verbose == True:
                    print(f"Selecting {mode} entry {input_code}")
            else:
                for infile_term in infile_list:
                    if infile_term.lower() in input_searchterm:
                        input.append(input_code)
                
        if len(input) == 0:
            raise SystemExit(f"No {mode} results or files were found for the searchterm {infile}")
        else:
            if quiet == False:
                print(f"Error: Searching {mode} list for {infile} returned {len(input)} KEGG entries")
        input_mode = "search"
    return input_mode


#----------------------------------------Process input----------------------------------------#
#openfile()
if mode == "ko":
    input_mode = openfile("K")
elif mode == "compound":
    input_mode = openfile("C")
elif mode == "reaction":
    input_mode = openfile("R")
elif mode == "module":
    input_mode = openfile("M")
elif mode == "mdata":
    input_mode = openfile("M")
elif mode == "homolog":
    input_mode = openfile("K")
input_unique = sorted([*set(input)])
input_total = len(input)
input_unique_total =len(input_unique)
total_reactions_list = []
total_modules_list = []
total_compounds_list = []
total_glycans_list = []
dict_pathways = {}
orgdict = {}
homolog_orglist = []
#----------------------------------------Process orglist for homolog arguments----------------------------------------#
def get_orglist():
    orgfile = os.path.join(dir_download, "KEGG_lists", "organisms.txt")
    if os.path.exists(orgfile) == False or overwrite == True:
        #Retrieve the KEGG organisms list from API and save to downloads folder
        url = "https://rest.kegg.jp/list/organism"
        if verbose == True:
            print(f"Downloading KEGG organisms list from {url}")
        req = requests.get(url).text
        with open(orgfile, "w") as handle:
            handle.write(req)
            if verbose == True:
                print(f"Wrote KEGG organisms list to {orgfile}")
            orgdata = req.splitlines()
    else:
        #Retrieve the KEGG organisms list from the downloads folder
        with open(orgfile, "r") as req:
            if verbose == True:
                print(f"Reading KEGG organisms list from {orgfile}")
            orgdata = req.readlines()
    
    #Process organism list data
    for item in orgdata:
        data = item.split("\t")
        orgcode = data[1]
        orgstring = f"{data[0]};{data[2]};{data[3]}"
        orgdict[orgcode] = orgstring



def open_homologfile(homologfile):
    with open(homologfile, "r") as file:
            while (line := file.readline().strip()):
                homolog_orglist.append(line)
                
if homolog == True:
    if homologorg is not None:
        get_orglist()
        if os.path.isfile(homologorg) == True:
            homologfile = homologorg
            open_homologfile(homologfile)
        elif os.path.isfile(homologorg) == False:
            homologfile = None
            homologorg_list = homologorg.split(",")
            #Strip whitespace
            homologorg_list = [homologitem.strip() for homologitem in homologorg_list]

            #Search descriptions of organism codes and append organism codes to homolog_orglist
            for homologorgitem in homologorg_list:
                for orgcode, orgstring in orgdict.items():
                    if homologorgitem.lower() == orgcode:
                        homolog_orglist.append(orgcode)
                    elif homologorgitem.lower() in orgstring.lower():
                        homolog_orglist.append(orgcode)
            #Raise error and exit if search term found no results
            if len(homolog_orglist) == 0:
                raise SystemExit(f"Error: No organism codes or descriptions matched the keyword {homologorg}")
        #No directories. Only files or keywords
        elif os.path.isdir(homologorg):
            raise SystemExit(f"Error: --homolog-organism command unable to process {homologorg} as it is a directory. Please select a file or enter a keyword search.")
                
            



#Create output directory structure
def create_outdirs(dir_list):
    for item in dir_list:
        path = os.path.join(outdir, item)
        os.makedirs(path, exist_ok = True)
    if verbose == True:
        print(f"Create output directories {dir_list} in {outdir}")
    

#Pathway summary writer
#Write KEGG pathway summary
def write_pathway_summary(filename):
    pathway_writer = csv.writer(filename, delimiter='\t')
    pathway_writer.writerow(["Pathway", "Name", "Count"])
    pathway_data = []
    for key, value in dict_pathways.items():
        pathway_data.append([key, value[0], value[1]])
    pathway_writer.writerows(pathway_data)



#----------------------------------------Print Header----------------------------------------#
if quiet == False:
    print(f"{linebreak}KEGGChem{linebreak}")
    print("Arguments given:")
    print(f"Mode: {mode}, Input mode: '{input_mode}', Input: {infile}, Output directory: {outdir}, Download directory: {dir_download}",
        f"Quiet: {quiet}, Verbose: {verbose}, Structure: {structure}",
        f"Pubchem: {pubchemarg}, SDF: {sdfarg}, Homolog: {homolog}, Homolog organism: {homologorg}",
        sep="\n")
    if homologorg is not None:
        if homologfile is None:
            print(f"Filtering homolog genes using the keyword {homologorg}")
        else:
            print(f"Filtering homolog genes by entries in file {homologfile}")
    


#----------------------------------------Specific mode functions----------------------------------------#

#--------------------KEGG KO/Module to compound--------------------
# Function to extract reaction codes from KEGG orthologue pages
def get_reaction_codes(KO):
    ko_file = os.path.join(dir_download, "KEGG_entries/orthologues", KO)
    #Check if file exists in download folder
    if os.path.exists(ko_file) == False or overwrite == True:
        url = f"https://rest.kegg.jp/get/{KO}"
        req = requests.get(url).text
        #Strip everything in req from "GENES" onwards
        r = req.split(separator, 1)[0]
        with open(ko_file, "w") as file:
            file.write(req)
    else:
        with open(ko_file, "r") as file:
            r = file.read().split(separator, 1)[0]
        if verbose == True:
            print(f"Extracted {KO} from file\n")

    # Extract reaction codes using regular expressions
    reaction_codes = [*set(re.findall(r"R\d{5}", r))]
    module_codes = [*set(re.findall(r"M\d{5}", r))]
    return reaction_codes, module_codes

#Function to retrieve genes from the list function of the KEGG API using KO codes
def retrieve_genes_from_ko(KO):
    url = f"https://rest.kegg.jp/link/genes/{KO}"
    linktext = os.path.join(dir_download, "KEGG_links/genes", KO)
    try:
        if os.path.exists(linktext) == False or overwrite == True:
            req = requests.get(url).text
            with open(linktext, "w") as filehandle:
                filehandle.write(req)
                if quiet == False:
                    print(f"Writing KO and gene links to file {linktext}")
        else:
            with open(linktext, "r") as filehandle:
                req = filehandle.read()
                if quiet == False:
                    print(f"Reading KO and gene links from file {linktext}")
    except:
        print("Error: File handling in recording KO-Gene links.")
        
    req = req.replace(f"ko:{KO}\t", "")
    gene_dict = {}
    #Set up counter and length for showing progress 
    length = len(req.splitlines())
    genecount = 1
    
    #Function to retrieve genes. Functioned so that it can be modified for homologfile command
    def retrieve_gene(genename):
        geneurl = f"https://rest.kegg.jp/get/{genename}"
        genefile = os.path.join(dir_download, "KEGG_entries/genes", genename)

        #Collect gene data
        def collect_gene_data(gene):
            AASEQ_pattern = re.compile(r'^AASEQ\s+(\d+)', re.MULTILINE)
            AASEQ_match = AASEQ_pattern.search(gene)
            NTSEQ_pattern = re.compile(r'^NTSEQ\s+(\d+)', re.MULTILINE)
            NTSEQ_match = NTSEQ_pattern.search(gene)
            NTSEQ_length = int(NTSEQ_match.group(1))
            AASEQ_length = int(AASEQ_match.group(1))
            return NTSEQ_length, AASEQ_length
        
        if os.path.exists(genefile) == False:
            gene = requests.get(geneurl).text
            NTSEQ_length, AASEQ_length = collect_gene_data(gene)
            with open(genefile, "w") as handle:
                if verbose == True:
                    print(f"Writing gene {genename} to file {genefile}")
                handle.write(gene)
                gene_dict[genename] = [gene, NTSEQ_length, AASEQ_length] 
        else:
            with open(genefile, "r") as handle:
                gene = handle.read()
                NTSEQ_length, AASEQ_length = collect_gene_data(gene)
                if verbose == True:
                    print(f"Retrieving gene file {genename} from download directory")
                gene_dict[genename] = [gene, NTSEQ_length, AASEQ_length]
        if quiet == False:
            print(f"Retrieved gene {genename} ({genecount}/{length})")
    
    #Parse lines using the retrieve_gene function
    genenames = []
    for line in req.splitlines():
        genename = line
        if homologorg is None:
            retrieve_gene(genename)
            genenames.append(genename)
        else: 
            #Find organism code
            org = line.split(":", 1)[0]
            if org in homolog_orglist:
                retrieve_gene(genename)
                genenames.append(line)
            else:
                if verbose == True:
                    print(f"{genename} ignored as it was not in the organism list")
        #retrieve_gene(genename)
        genecount += 1
    return gene_dict, genenames


#Homolog extraction
def extract_homolog(homolog_dict):
    if quiet == False:
        print(f"{linebreak}Extracting homolog data from orthologues{linebreak}")
    count_current = 1
    ko_genecount = {}
    for orthologue_id in input_unique:
        if quiet == False:
            print(f"{linebreak}Retrieving genes for {orthologue_id} ({count_current}/{input_unique_total})")
        
        gene_out = os.path.join(outdir, "Results/Genes", orthologue_id)
        
        
        try:
            gene_dict, genenames = retrieve_genes_from_ko(orthologue_id)
            ko_genecount[orthologue_id] = genenames
            homolog_dict[orthologue_id] = gene_dict
            if len(gene_dict) > 0:
                #Make output folder for orthologue
                os.makedirs(gene_out)
                #Write gene results to results folder for orthologue
                if verbose == True:
                    print(f"\nWriting genes to results folder {gene_out}\n")
                for genename, generesult in gene_dict.items():
                    #genefile = genename + ".txt"
                    with open(os.path.join(gene_out, genename + ".txt"), "w") as handle:
                        handle.write(generesult[0])
            else:
                with open(os.path.join(outdir, "Results/Genes", "KOs_without_genes.txt"), "a") as filehandle:
                    filehandle.write(f"{orthologue_id}\n")
            
        except:
            print(f"Error: Failed to retrieve genes for {orthologue_id}\n")
        
        if quiet == False:
            print(f"\nRetrieved genes for {orthologue_id} ({count_current}/{input_unique_total})\n")
        count_current += 1
        #End homolog extraction step
    return ko_genecount



# Function to retrieve compound/glycan codes associated with a reaction or module code
def get_compound_codes(query_code):
    #Ensure that queries are saved to/loaded from the correct directory
    if re.search(rf"R\d{{5}}\b", query_code):
        query_file = os.path.join(dir_download, "KEGG_entries/reactions", query_code)
    elif re.search(rf"M\d{{5}}\b", query_code):
        query_file = os.path.join(dir_download, "KEGG_entries/modules", query_code)

    if os.path.exists(query_file) == False or overwrite == True:
        url = f"http://rest.kegg.jp/get/{query_code}"
        req = requests.get(url).text
        #Strip everything in req from "GENES" onwards
        if mode == "ko":
            r = req.split(separator, 1)[0]
        else:
            if "COMPOUND" in req:
                r = req.split("COMPOUND", 1)[1]
            else:
                r = req
        with open(query_file, "w") as file:
            file.write(req)
        if verbose == True:
            print(f"Wrote {query_code} to file\n")
    else:
        with open(query_file, "r") as file:
            if mode == "ko":
                r = file.read().split(separator, 1)[0]
            else:
                req = file.read()
                if "COMPOUND" in req:
                    r = req.split("COMPOUND", 1)[1]
                else:
                    r = req
        if verbose == True:
            print(f"Extracted {query_code} from file\n")
    #Extract compound codes from page
    compound_codes = [*set(re.findall(r"\b[CG]\d{5}\b", r))]
    return compound_codes


#Function to retrieve compound data from compound codes
def get_compound_data(compound_code):
    compound_file = os.path.join(dir_download, "KEGG_entries/compounds", compound_code)
    if os.path.exists(compound_file) == False  or overwrite == True:
        url = f"http://rest.kegg.jp/get/{compound_code}"
        req = requests.get(url).text   
        #r = req.split(separator, 1)[0]
        with open(compound_file, "w") as file:
            file.write(req)
        if verbose == True:
            print(f"Wrote {compound_code} to file\n")
    else:
        with open(compound_file, "r") as file:
            req = file.read()#.split(separator, 1)[0]
        if verbose == True:
            print(f"Extracting {compound_code} from file")
    #Retrieve compound formula
    formula_string = re.findall(r"FORMULA     .*", req)
    if len(formula_string) > 0:
        compound_formula = formula_string[0].replace("FORMULA     ", "")
    else:
        compound_formula = "NULL"

    #Retrieve exact mass
    exact_mass_string = re.findall(r"EXACT_MASS  .*", req)
    if len (exact_mass_string) > 0:
        exact_mass = float(exact_mass_string[0].replace("EXACT_MASS  ", ""))
    else:
        exact_mass = "NULL"

    #Retrieve mol. weight
    mol_weight_string = re.findall(r"MOL_WEIGHT  .*", req)
    if len(mol_weight_string) > 0:
        mol_weight = float(mol_weight_string[0].replace("MOL_WEIGHT  ", ""))
    else:
        mol_weight = "NULL"
    
    #Retrieve pathway data
    map_string = re.findall(r"map\d{5}  .*", req)
    map_codes = []
    if len(map_string) > 0:
        #Split pathway map names and append to dictionary
        for map in map_string:
            split_map = map.split("  ")
            map_code = split_map[0]
            map_name = split_map[1]
            if map_code not in dict_pathways:
                dict_pathways[map_code] = [map_name, 1]
            else:
                dict_pathways[map_code][1] += 1
            map_codes.append(map_code)
        #Append map codes only to list and return
        pathway_map = map_codes
    else:
        pathway_map = "NULL"

    

    #Retrieve database data (ChEBI, PubChem)
    def find_database_data(database):
        database_string = re.findall(rf"{database}: .*", req)
        #print(database_string)
        if len(database_string) > 0:
            #Remove leading header (eg. ChEBI:)
            database = database_string[0].replace(f"{database}: ", "")
            database = database.split(" ")
            #database = [i for i in database]
            #Split codes into list, remove leading and trailing spaces
            #database = database.strip()
            return database
        else:
            return "NULL"            
    ChEBI = find_database_data("ChEBI")
    PubChem = find_database_data("PubChem")
    


    return compound_formula, exact_mass, mol_weight, ChEBI, PubChem, pathway_map

def get_glycan_data(glycan_code):
    compound_file = os.path.join(dir_download, "KEGG_entries/compounds", glycan_code)
    if os.path.exists(compound_file) == False or overwrite == True:
        url = f"http://rest.kegg.jp/get/{glycan_code}"
        req = requests.get(url).text   
        #r = req.split(separator, 1)[0]
        with open(compound_file, "w") as file:
            file.write(req)
        if verbose == True:
            print(f"Wrote {glycan_code} to file\n")
    else:
        with open(compound_file, "r") as file:
            req = file.read()#.split(separator, 1)[0]
        if verbose == True:
            print(f"Extracting {glycan_code} from file\n")
    #Retrieve glycan composition
    composition_string = re.findall(r"COMPOSITION .*", req)
    if len(composition_string) > 0:
        composition = composition_string[0].replace("COMPOSITION ", "")
    else:
        composition = "NULL"

    #Retrieve exact mass
    exact_mass_string = re.findall(r"EXACT_MASS  .*", req)
    if len (exact_mass_string) > 0:
        exact_mass = float(exact_mass_string[0].replace("EXACT_MASS  ", ""))
    else:
        exact_mass = "NULL"


#Function to retrieve reaction data from reaction codes
def get_reaction_data(reaction_code):
    reaction_file = os.path.join(dir_download, "KEGG_entries/reactions", reaction_code)
    if os.path.exists(reaction_file) == False or overwrite == True:
        url = f"http://rest.kegg.jp/get/{reaction_code}"
        req = requests.get(url).text
        with open(reaction_file, "w") as file:
            file.write(req)
        if verbose == True:
            print(f"Wrote {reaction_code} to file\n")
    else:
        with open (reaction_file, "r") as file:
            req = file.read()
        if verbose == True:
            print(f"Extracted {reaction_code} from file\n")
    #Check if data present
    if len(req) > 0:
        reaction_result = True
    else:
        reaction_result = False
    #Extract name
    name_string = re.findall(r"NAME        .*", req)
    if len(name_string) > 0:
        r_name = name_string[0].replace("NAME        ", "")
    else:
        r_name = "NULL"
    
    #Extract equation
    equation_string = re.findall(r"EQUATION    .*", req)
    if len(equation_string) > 0:
        equation = equation_string[0].replace("EQUATION    ", "")
    else:
        equation = "NULL"
    #Extract definition
    definition_string = re.findall(r"DEFINITION  .*", req)
    if len(definition_string) > 0:
        definition = definition_string[0].replace("DEFINITION  ", "")
    else:
        definition = "NULL"
    rclass = [*set(re.findall(r"RC\d{5} .*", req))]
    ko_codes = [*set(re.findall(r"K\d{5}", req))]
    pathways = [*set(re.findall(r"rn\d{5}  .*", req))]
    return reaction_result, r_name, equation, definition, rclass, ko_codes, pathways


#---------------------Structure download functions----------------------#
def download_structure_files(compound_code, mode):
    filename = f"{compound_code}.{mode}"
    filedir = os.path.join(dir_download, "Structure_files/Downloaded", filename)
    nulldir = os.path.join(dir_download, "Structure_files/nullfiles", filename)
    if mode == "mol":
        urlcode = "-f+m+"
    elif mode == "kcf":
        urlcode = "-f+k+"

    if (os.path.exists(filedir) == False and os.path.exists(nulldir) == False) or overwrite == True:
        url = f"https://www.kegg.jp/entry/{urlcode}{compound_code}"
        req = requests.get(url).text
        if len(req) > 0:
            with open(filedir, "w") as file:
                file.write(req)
            mol_status = "Saved"
            mol_output = req
            if verbose == True:
                print(f"Saved {compound_code} to {filedir}\n")
        else:
            with open(nulldir, "w"):
                pass
            mol_status = "NULL"
            mol_output = "NULL"
            if verbose == True:
                print(f"Saved null result for {compound_code} to {nulldir}")
    elif os.path.exists(filedir) == True:
        with open(filedir, "r") as file:
            req = file.read()
        mol_output = req
        mol_status = "Loaded"
        if verbose == True:
            print(f"Read {compound_code} from {filedir}\n")
    elif os.path.exists(nulldir) == True:
        mol_status = "NULL"
        mol_output = "NULL"
        if verbose == True:
            print(f"Read null result for {compound_code} from {nulldir}")
    return mol_status, mol_output, filename

#Download structure files
def structure_download(compound):
    path = os.path.join(outdir, "Results", "Structure_files")
    os.makedirs(path, exist_ok=True)
    molstatus, moldata, molfile = download_structure_files(compound, "mol")
    if molstatus != "NULL":
        with open(os.path.join(path, molfile), "w") as file:
            file.write(moldata)

    kcfdata, kcfdata, kcffile = download_structure_files(compound, "kcf")
    if molstatus != "NULL":
        with open(os.path.join(path, kcffile), "w") as file:
            file.write(kcfdata)

def structure_file_download(totallist):
    #Download structure files
    if structure == True:
        structurecount = 1
        structuretotal = len(totallist)
        if quiet == False:
            print(f"{linebreak}Downloading compound structure data for {structuretotal} structures{linebreak}")       
        for compound in totallist:
            if quiet == False:
                print(f"Retrieving compound structure data for: {compound} ({structurecount}/{structuretotal})")
            structure_download(compound) 
            structurecount += 1
        if quiet == False:
            print(f"{linebreak}Finished downloading compound structure data{linebreak}")

#----------------------------------------Function - Retrieve homologs from KO----------------------------------------#
def ko_to_homolog():
    #Create outdir structure
    dir_list = ["Results/Genes"]
    create_outdirs(dir_list)



#----------------------------------------Function - Run KO to compound----------------------------------------#
def ko_to_compound():
    #Create outdir structure
    dir_list = ["Results/Individual_results", "Summaries"]
    if homolog == True:
        dir_list.append("Results/Genes")
    create_outdirs(dir_list)

    count_current = 1
    dict_reaction_codes = {}
    dict_module_codes = {}
    module_compounds_list = []
    reaction_compounds_list = []
    module_glycan_list = []
    reaction_glycan_list = []
    homolog_dict = {}

    if quiet == False:
        print(f"{linebreak}Analysing {input_unique_total} unique KEGG orthologues from ",\
            f"{input_mode} {infile} ({input_total} given){linebreak}", sep="")
    
    #Iterate over all orthologues and extract reaction codes and associated compound codes
    for orthologue_id in input_unique:
        #Print output header
        if quiet == False:
            print(f"Extracting reaction and compound codes from {orthologue_id}")
        reaction_codes, module_codes = get_reaction_codes(orthologue_id)
        if verbose == True:
            print(f"Extracted reaction and module codes from {orthologue_id}")
        outfile = orthologue_id + ".txt"
        
        with open(os.path.join(outdir, "Results/Individual_results", outfile), "a") as out_all,\
        open(os.path.join(outdir, "Results", "all_reaction_results.txt"), "a") as all_reaction_results,\
            open(os.path.join(outdir, "Results", "all_module_results.txt"), "a") as all_module_results:
            out_all.write(f"{orthologue_id}\n")
            all_reaction_results.write(f"{orthologue_id}\n")
            all_module_results.write(f"{orthologue_id}\n")
            #Iterate through reaction/module codes and retrieve compound codes before writing to file
            def extract_code(codetype, listname, dictionary, filename, totalcompounds, totalglycans):
                for itemcode in codetype:
                    #If reaction alrady stored, use previously searched results
                    if itemcode in dictionary:
                        listname.append(itemcode)
                        compound_codes = dictionary[itemcode]
                    else:
                        compound_codes = get_compound_codes(itemcode)
                        listname.append(itemcode)
                        #Add to dictionary of reactions so as to not repeat reaction page search
                        dictionary[itemcode] = compound_codes
                    #Add compound codes to set
                    for compound_code in compound_codes:
                        if compound_code[0] == "C":
                            #Append compound to total list for both types
                            total_compounds_list.append(compound_code)
                            totalcompounds.append(compound_code)
                            #Append compound to total list for individual mode
                            #totalcompounds.append(compound_code)
                        elif compound_code[0] == "G":
                            total_glycans_list.append(compound_code)
                            totalglycans.append(compound_code)
                    if verbose == True:
                        print(f"Extracted compound codes from {itemcode}")
                    out_all.write(f"{itemcode}:{','.join(compound_codes)}\n")
                    #all_reaction_results.write(f"{itemcode}:{','.join(compound_codes)}\n")
                    filename.write(f"{itemcode}:{','.join(compound_codes)}\n")
                #return compound_codes
            extract_code(reaction_codes, total_reactions_list, dict_reaction_codes, all_reaction_results, reaction_compounds_list, reaction_glycan_list)
            extract_code(module_codes, total_compounds_list, dict_module_codes, all_module_results, module_compounds_list, module_glycan_list)

        #Print current count in order to show progress
        if quiet == False:
            print(f"\nWrote {orthologue_id} to file ({count_current}/{input_unique_total}){linebreak}")
        count_current += 1



        
    #Set unique lists of total compounds/glycans, as well as those sourced through modules or reactions
    compounds_unique = sorted([*set(total_compounds_list)])
    glycans_unique = sorted([*set(total_glycans_list)])
    reactions_unique = sorted([*set(total_reactions_list)])
    reaction_compounds_unique = sorted([*set(reaction_compounds_list)])
    module_compounds_unique = sorted([*set(module_compounds_list)])
    reaction_glycans_unique = sorted([*set(reaction_glycan_list)])
    module_glycans_unique = sorted([*set(module_glycan_list)])

    #Write unique compounds to list
    with open(os.path.join(outdir, "Results", "compounds.txt"), "w") as result,\
        open(os.path.join(outdir, "Results", "glycans.txt"), "w") as glycans,\
            open(os.path.join(outdir,"Results","reactions.txt"), "w") as reactions,\
                open(os.path.join(outdir, "Results", "compounds_from_reactions.txt"), "w") as reaction_compounds,\
                    open(os.path.join(outdir, "Results", "compounds_from_modules.txt"), "w") as module_compounds,\
                        open(os.path.join(outdir, "Results", "glycans_from_reactions.txt"), "w") as reaction_glycans,\
                            open(os.path.join(outdir, "Results", "glycans_from_modules.txt"), "w") as module_glycans:
        for compound in compounds_unique:
            result.write(f"{compound}\n")
        for glycan in glycans_unique:
            glycans.write(f"{glycan}\n")
        for reaction in reactions_unique:
            reactions.write(f"{reaction}\n")
        for compound in reaction_compounds_unique:
            reaction_compounds.write(f"{compound}\n")
        for compound in module_compounds_unique:
            module_compounds.write(f"{compound}\n")
        for glycan in reaction_glycans_unique:
            reaction_glycans.write(f"{glycan}\n")
        for glycan in module_glycans_unique:
            module_glycans.write(f"{glycan}\n")

    total_compounds_dict = dict(Counter(total_compounds_list))
    compounds_count = len(total_compounds_dict)
    total_reactions_dict = dict(Counter(total_reactions_list))
    reactions_count = len(total_reactions_dict)
    total_glycans_dict = dict(Counter(total_glycans_list))
    glycans_count = len(total_glycans_dict)

    

    #Maths! (Intersection of module/reaction lists)
    def intersection(lst1, lst2):
        lst3 = [value for value in lst1 if value in lst2]
        lst3_count = len(lst3)
        lst1_unique = len(lst1) - lst3_count
        lst2_unique = len(lst2) - lst3_count
        return lst3_count, lst1_unique, lst2_unique
    
    compound_common_count, compound_module_count, compound_reaction_count = intersection(module_compounds_unique, reaction_compounds_unique)
    glycan_common_count, glycan_module_count, glycan_reaction_count = intersection(module_glycans_unique, reaction_glycans_unique)

    

    #Download structure files
    structure_file_download(compounds_unique)

    #Write out summary text
    with open(os.path.join(outdir, "Summaries", "compound_summary.txt"), "a") as sum_compounds,\
    open(os.path.join(outdir, "Summaries", "reaction_summary.txt"), "a") as sum_reactions,\
        open(os.path.join(outdir, "Summaries", "glycan_summary.txt"), "w") as sum_glycans:
        

        #Write reaction summary results to file
        sum_reactions.write(f"{linebreak}Reactions summary:{linebreak}")
        sum_reactions.write(f"Total unique reactions: {reactions_count}\n\n")
        sum_reactions.write("Reaction: count\n")
        for key, value in total_reactions_dict.items():
            sum_reactions.write(f"{key}: {value}\n")

        #Function to write compound/glycan summary results to file
        def write_summaries_to_file(itemtype, filename, total_count, module_total, reaction_total, 
                                    count_module, count_reaction, count_common, dict_name):
            filename.write(f"{linebreak}{itemtype}s summary:{linebreak}")
            filename.write(f"Total unique {itemtype}s: {total_count}\n\n")
            filename.write(f"Total unique {itemtype}s inferred through modules: {len(module_total)}\n")
            filename.write(f"Total unique {itemtype}s inferred through reactions: {len(reaction_total)}\n")
            filename.write(f"Unique {itemtype}s inferred from modules only: {count_module}\n")
            filename.write(f"Unique {itemtype}s inferred from reactions only: {count_reaction}\n")
            filename.write(f"{itemtype}s inferred through both modules and reactions: {count_common}\n")
            filename.write(f"{linebreak}Number of occurences for each unique {itemtype}{linebreak}")
            filename.write(f"{itemtype}: count\n")
            for key, value in dict_name.items():
                filename.write(f"{key}:{value}\n")
        
        #Write compound/glycan summary results to filegene_dict
        write_summaries_to_file("Compound", sum_compounds, compounds_count, module_compounds_unique,
                                reaction_compounds_unique, compound_module_count, compound_reaction_count,
                                compound_common_count, total_compounds_dict)
        write_summaries_to_file("Glycan", sum_glycans, glycans_count, module_glycans_unique, 
                                reaction_glycans_unique, glycan_module_count, glycan_reaction_count, 
                                glycan_common_count, total_glycans_dict)

    #Append summarise results dictionaries into list
    ko_results = [total_reactions_dict, total_compounds_dict, total_glycans_dict]
    #Extract homologs for each gene
    if homolog == True:
        ko_genelist = extract_homolog(homolog_dict)
        
        #Do some stats on the genes returned
        ko_genes = []
        for geneitem in ko_genelist.values():
            for gene in geneitem:
                ko_genes.append(gene)
        total_genes_dict = dict(Counter(ko_genes))
        #Append ko_genelist to ko_results so that it can be displayed in the log
        ko_results.append(total_genes_dict)

        with open(os.path.join(outdir, "Summaries/homolog_summary.txt"), "w") as handle:
            handle.write(f"{linebreak}Homolog gene summary{linebreak}")
            handle.write(f"KO\tGenes returned{linebreak}")
            for key, value in ko_genelist.items():
                handle.write(f"{key}\t{len(value)}\n")
            handle.write(f"{linebreak}Gene\tCount\n")
            for key, value in total_genes_dict.items():
                handle.write(f"{key}\t{value}\n")


    #Return results list for log summary
    return ko_results



#----------------------------------------Function - Run module to compound----------------------------------------#
def module_to_compound():
    count_current = 1
    #Create outdir structure
    dir_list = ["Results/Individual_results", "Summaries"]
    create_outdirs(dir_list)
    
    if quiet == False:
        print(f"{linebreak}Analysing {input_unique_total} unique KEGG modules from ",\
            f"{input_mode} {infile} ({input_total} given){linebreak}", sep="")
    
    for module in input_unique:
        print(f"Extracting compounds from KEGG module {module}\n")
        compound_codes = get_compound_codes(module)
        outfile = f"{module}.txt"

        with open(os.path.join(outdir, "Results/Individual_results", outfile),"a") as result,\
                open(os.path.join(outdir, "Results", "all_results.txt"), "a") as summary:
                    #Write file header
                    for compound_code in compound_codes:
                        if compound_code[0] == "C":
                            total_compounds_list.append(compound_code)
                        elif compound_code[0] == "G":
                            total_glycans_list.append(compound_code)
                        
                    result.write(f"{module}: {','.join(compound_codes)}\n")
                    summary.write(f"{module}: {','.join(compound_codes)}\n")
                    if quiet == False:
                        print(f"Extracted compound codes from {module} ({count_current}/{input_unique_total}){linebreak}")
                    count_current += 1

    compounds_unique = sorted([*set(total_compounds_list)])
    glycans_unique = sorted([*set(total_glycans_list)])
    #Write unique compounds to list
    with open(os.path.join(outdir, "Results", "compounds.txt"), "w") as result,\
        open(os.path.join(outdir, "Results", "glycans.txt"), "w") as glycans:
        for compound in compounds_unique:
            result.write(f"{compound}\n")
        for glycan in glycans_unique:
            glycans.write(f"{glycan}\n")

    total_compounds_dict = dict(Counter(total_compounds_list))
    compounds_count = len(total_compounds_dict)
    total_glycans_dict = dict(Counter(total_glycans_list))
    glycans_count = len(total_glycans_dict)
    
    #Download structure files
    structure_file_download(compounds_unique)
    if structure == True:
        structurecount = 1
        structuretotal = len(compounds_unique)
        if quiet == False:
            print(f"{linebreak}Downloading compound structure data for {structuretotal} structures{linebreak}")       
        for compound in compounds_unique:
            if quiet == False:
                print(f"Retrieving compound structure data for: {compound} ({structurecount}/{structuretotal})")
            structure_download(compound) 
            structurecount += 1
            if quiet == False:
                print(f"{linebreak}Finished downloading compound structure data{linebreak}")
    
    #Write out summary text
    with open(os.path.join(outdir, "Summaries", "compound_summary.txt"), "w") as sum_compound,\
            open(os.path.join(outdir, "Summaries", "glycan_summary.txt"), "w") as sum_glycan:

        #Write compound summary results to file
        sum_compound.write(f"{linebreak}Compounds summary:{linebreak}")
        sum_compound.write(f"Total unique compounds: {compounds_count}\n\n")
        sum_compound.write("Compound\tcount\n")
        for key, value in total_compounds_dict.items():
            sum_compound.write(f"{key}\t{value}\n")
        #Write gycan summary results to file
        sum_glycan.write(f"{linebreak}Glycan summary:{linebreak}")
        sum_glycan.write(f"Total unique glycans: {glycans_count}\n\n")
        sum_glycan.write("Glycan\tcount\n")
        for key, value in total_glycans_dict.items():
            sum_glycan.write(f"{key}\t{value}\n")
    
    module_results = [total_compounds_dict, total_glycans_dict]
    return module_results
    


#----------------------------------------Function - Run compound data retrieve----------------------------------------
def compound_data():
    count_current = 1
    #Create outdir
    dir_list = ["Results", "Summaries"]
    if sdfarg == True:
        dir_list.append("Results/SDF")
    create_outdirs(dir_list)

    compound_formula_list = []
    compound_formula_null = []
    lst_exact_mass = []
    lst_mol_weight = []
    lst_chebi = []
    lst_pubchemsid = []
    lst_pathways = []

    print(f"{linebreak}Retrieving data for {input_unique_total} unique KEGG compounds from ",\
                    f"{input_mode} {infile} ({input_total} given){linebreak}", sep="")

    with open(os.path.join(outdir, "Results", "Compound_data.tsv"), "a") as out:
        out_csv = []
        out_write = csv.writer(out, delimiter='\t')
        #out.write("Compound;formula;exact mass;molecular weight;ChEBI ID;PubChem SID;pathways\n")
        out_header = ["compound", "formula", "exact mass", "molecular weight", "ChEBI ID", "PubChem SID", "KEGG pathways"]
        #formulas.write("Compound:formula\n")
        #weight.write("compound;exact_mass;mol_weight\n")
        for compound in input_unique:
            compound_formula,exact_mass,mol_weight, ChEBI, PubChem, pathway_map = get_compound_data(compound)
            #Add weight to lists for maths later
            if exact_mass != "NULL":
                lst_exact_mass.append(exact_mass)
            if mol_weight != "NULL":
                lst_mol_weight.append(mol_weight)
            #Add pathways to list
            if pathway_map != 'NULL':
                for item in pathway_map:
                    lst_pathways.append(item)
                pathways = ",".join(pathway_map)
            #Write formula result
            if compound_formula != 'NULL':
                #formulas.write(f"{compound}:{compound_formula}\n")
                compound_formula_list.append(compound)
                if quiet == False:
                    print(f"Data retrieved for {compound}: {compound_formula} ({count_current}/{input_unique_total})\n")
            else:
                #formulas.write(f"{compound}:NULL\n")
                compound_formula_null.append(compound)
                if quiet == False:
                    print(f"Data retrieved for {compound}: NULL ({count_current}/{input_unique_total})\n")
            #Write weight results
            #weight.write(f"{compound};{exact_mass};{mol_weight}\n")
            if ChEBI != "NULL":
                chebi = ';'.join(ChEBI)
                #chebi = ChEBI
                for item in ChEBI:
                    lst_chebi.append(item)
            else:
                chebi = 'NULL'
            if PubChem != 'NULL':
                pubchem = ';'.join(PubChem)
                #pubchem = PubChem
                for item in PubChem:
                    lst_pubchemsid.append(item)
            else:
                pubchem = "NULL"
            out_data = [compound, compound_formula, exact_mass, mol_weight, chebi, pubchem, pathways]
            out_csv.append(out_data)
            #out.write(f"{compound};{compound_formula};{exact_mass};{mol_weight};{chebi};{pubchem};{pathways}\n")

            count_current += 1
        #Write to file
        if verbose == True:
            print(f"{linebreak}Writing compound results to file{linebreak}")
        out_write.writerow(out_header)
        out_write.writerows(out_csv)
        print(f"{linebreak}Finished retrieving compound data{linebreak}")

    #Retrieve PubChem compound data using SID codes via PubChemPy
    #Create unique list of pubchemsid codes
    lst_pubchemsid_unique = sorted([*set(lst_pubchemsid)])
    if pubchemarg == True:
        lst_CID = []
        #dict_CID = {}
        dict_SID = {}
        
        #Establish counter for sid codes for progress reports
        sidcount = 1
        sidcount_total = len(lst_pubchemsid_unique)
        #Open PubChem results file and search PubChem Substance database
        with open(os.path.join(outdir, "Results", "PubChem_substances.txt"), "a") as SIDout:
            print(f"{linebreak}Retrieving PubChem SID data from PubChem Substance database{linebreak}")
            if verbose == True:
                print(f"Opened PubChem substance results file at {SIDout}")
            SIDout.write("SID\tSID source ID\tCID\n")
            for sid in lst_pubchemsid_unique:
                try:
                    if sid in dict_SID:
                        sid_source_id = dict_SID[sid][0]
                        sid_CID = dict_SID[sid][1]
                    else:
                        sid_data = pcp.Substance.from_sid(sid)
                        #sid_synonyms = sid_data.synonyms
                        sid_source_id = sid_data.source_id
                        sid_CID = sid_data.standardized_cid
                        dict_SID[sid] = [sid_source_id, sid_CID]
                    if sid_CID is not None:
                        lst_CID.append(sid_CID)
                    if quiet == False:
                        print(f"Retrieved substance data for SID: {sid} from PubChem Substance database ({sidcount}/{sidcount_total})")
                    SIDout.write(f"{sid}\t{sid_source_id}\t{sid_CID}\n")
                    sidcount += 1
                except:
                    print(f"Invalid SID code searched: {sid}")
        with open(os.path.join(outdir, "Results", "PubChem_compounds.txt"), "a") as CIDout:
            print(f"{linebreak}Retrieving PubChem CID data from PubChem Compound database{linebreak}")
            if verbose == True:
                print(f"Opened PubChem Compound results file at {CIDout}")
            CIDout.write("CID\tMolecular formula\tCanonical SMILES\tIsomeric SMILES\n")
            lst_CID_unique = sorted([*set(lst_CID)])
            cidcount = 1
            cidcounttotal = len(lst_CID_unique)
            for cid in lst_CID_unique:
                try:
                    if sdfarg == True:
                        with open(os.path.join(outdir, "Results/SDF", f"{cid}.txt"), "a") as sdfout:
                            sdf_data = pcp.get_sdf(cid)
                            sdfout.write(sdf_data)
                    cid_data = pcp.Compound.from_cid(cid)
                    cid_formula = cid_data.molecular_formula
                    cid_canon_SMILES = cid_data.canonical_smiles
                    cid_isomeric_SMILES = cid_data.isomeric_smiles
                    CIDout.write(f"{cid}\t{cid_formula}\t{cid_canon_SMILES}\t{cid_isomeric_SMILES}\n")
                    if quiet == False:
                        print(f"Retrieved compound data for CID: {cid} from PubChem Compound database ({cidcount}/{cidcounttotal})")
                    cidcount += 1
                except:
                    print(f"Invalid CID code searched: {cid}")
        print(f"{linebreak}Completed retrieving SID and CID data from PubChem{linebreak}")

    #Download structure files
    structure_file_download(input_unique)

    #Write summary
    with open(os.path.join(outdir, "Summaries", "pathway_summary.tsv"), "a") as path_summary:
        write_pathway_summary(path_summary)
    compound_results = [compound_formula_list, compound_formula_null, lst_exact_mass, lst_mol_weight]
    return compound_results
        
#----------------------------------------Function - Run reaction data retrieve----------------------------------------
def reaction_data():
    count_current = 1
    #Create outdir
    dir_list = ["Results", "Summaries"]
    create_outdirs(dir_list)

    reaction_names_list = []
    reaction_names_null = []
    reaction_result_true = []
    reaction_result_false = []
    lst_rclass = []
    lst_ko = []
    lst_pathways = []
    lst_pathway_codes = []

    print(f"{linebreak}Retrieving data for {input_unique_total} unique KEGG reactions from ",\
                    f"{input_mode} {infile} ({input_total} given){linebreak}", sep="")

    with open(os.path.join(outdir, "Results", "Equations.txt"), "w") as f_equation,\
        open(os.path.join(outdir, "Results", "Definitions.txt"), "w") as f_definition,\
            open(os.path.join(outdir, "Results", "Rclass.txt"), "w") as f_rclass,\
                open(os.path.join(outdir, "Results", "KO_codes.txt"), "w") as f_ko,\
                    open(os.path.join(outdir, "Results", "Pathways.txt"), "w") as f_pathway:
        f_equation.write("Reaction\tequation\n")
        f_definition.write("Reaction\tdefinition\n")
        f_rclass.write("Reaction\treaction classes\n")
        f_ko.write("Reaction\tKO codes\n")
        f_pathway.write("Reaction\tPathways\n")
        
        for reaction in input_unique:
            #Write equation results
            reaction_result, r_name, equation, definition, rclasses, ko_codes, pathways = get_reaction_data(reaction)
            #Check if any result found
            if reaction_result == True:
                reaction_result_true.append(reaction)
            else:
                reaction_result_false.append(reaction)
            
            #Reaction name
            if len(r_name) > 0:
                reaction_names_list.append(reaction)
            else:
                reaction_names_null.append(reaction)
            #Equation
            if len(equation) > 0:
                f_equation.write(f"{reaction}\t{equation}\n")
            else:
                f_equation.write(f"{reaction}\tNULL\n")
            #Write definition results
            if len(definition) > 0:
                f_definition.write(f"{reaction}\t{definition}\n")
            else:
                f_definition.write(f"{reaction}\tNULL\n")
            #Write RCLASS results
            if len(rclasses) > 0:
                for rclass in rclasses:
                    lst_rclass.append(rclass)
                f_rclass.write(f"{reaction}\t{', '.join(rclasses)}\n")
            else:
                f_rclass.write(f"{reaction}\tNULL\n")
            #Write KO codes
            if len(ko_codes) > 0:
                for ko_code in ko_codes:
                    lst_ko.append(ko_code)
                f_ko.write(f"{reaction}\t{','.join(ko_codes)}\n")
            else:
                f_ko.write(f"{reaction}\tNULL\n")
            #Write pathways
            if len(pathways) > 0:
                pathway_codes = [re.sub("  .*", "", item) for item in pathways]
                for pathway_code in pathway_codes:
                    lst_pathway_codes.append(pathway_code)
                f_pathway.write(f"{reaction}\t{','.join(pathway_codes)}\n")
                for pathway in pathways:
                    lst_pathways.append(pathway)
            else:
                f_pathway.write(f"{reaction}\tNULL\n")
                
            if quiet == False:
                print(f"Data retrieved for {reaction}: ({count_current}/{input_unique_total})\n")
            count_current += 1
        print(f"{linebreak}Finished retrieving reaction data{linebreak}")

    #Download structure files
    if structure == True:
        print("Structure retrieval is not run with KEGG reaction retrieve")
    reaction_results_list = [reaction_result_true, reaction_result_false, reaction_names_list, reaction_names_null]
    return reaction_results_list

        
#------------------------------------------------------------Run module data retrieve---------------------------------------------------------------#
#Function for retrieving module data
def get_module_data(module_code):
    module_file = os.path.join(dir_download, "KEGG_entries/modules", module_code)
    if os.path.exists(module_file) == False or overwrite == True:
        url = f"http://rest.kegg.jp/get/{module_code}"
        req = requests.get(url).text
        with open(module_file, "w") as handle:
            handle.write(req)
    else:
        with open(module_file, "r") as file:
            req = file.read()
    if len(req) > 0:
        entry = re.search(r"ENTRY       .*", req)
        moduletype = entry.group(0).replace("ENTRY       ", "")
        moduletype = moduletype.replace("            ",":")
        name = re.search(r"NAME        .*", req)
        modulename = name.group(0).replace("NAME        ","")
        print(f"{moduletype}:{modulename}")
    else:
        print(f"{module_code}: No data found")

def module_data():
    #count_current = 1
    #Create outdir
    #dir_list = ["Results", "Summaries"]
    #create_outdirs(dir_list)



    print(f"{linebreak}Retrieving data for {input_unique_total} unique KEGG modules from ",\
                    f"{infile} ({input_total} given){linebreak}", sep="")
    for module in input_unique:
        get_module_data(module)
        #print(f"{module}: {moduletype}: {name}")
        #count_current += 1
        #print(f"{linebreak}Finished retrieving module data {module}{linebreak}")

    #Download structure files
    if structure == True:
        print("Structure retrieval is not run with KEGG module retrieve")


 #----------------------------------------Function to write log----------------------------------------#
def write_log():
    #Calculate runtime
    runtime, timestamp_end = timetaken()

    keywords = {}
    keywords["arguments"] = args
    keywords["infile"] = infile
    keywords['mode'] = mode
    keywords["overwrite"] = overwrite
    keywords["verbose"] = verbose
    keywords["filepath"] = os.path.abspath(outdir)
    keywords["download"] = download
    keywords["timestamp_start"] = timestamp
    keywords["runtime"] = runtime
    keywords["timestamp_end"] = timestamp_end
    keywords["pubchemarg"] = pubchemarg
    keywords["sdfarg"] = sdfarg
    keywords["homolog"] = homolog
    keywords["homologorg"] = homologorg
    
    if input_mode == "file":
        input_term = os.path.abspath(infile)
    elif input_mode == "search":
        input_term = (infile)
    with open(os.path.join(outdir, "run_summary.txt"), "w") as outlog:
        outlog.write(f"{linebreak}KEGGChem log{linebreak}")
        outlog.write(f"Arguments: {sys.argv}\n")
        outlog.write(f"Start time: {timestamp}\nEnd time: {timestamp_end}\nRun time: {runtime}\n")
        outlog.write(f"Input mode: {input_mode}\nInput: {input_term}\nMode: {mode}\nOutput directory: {os.path.abspath(outdir)}\nDownload directory: {os.path.abspath(dir_download)}\n")
        outlog.write(f"Verbose: {verbose}\nQuiet: {quiet}\nOverwrite: {overwrite}\n\n")
        outlog.write(f"Total number of valid entries: {input_total}\nTotal number of unique entries: {input_unique_total}\nInvalid entries: {len(input_invalid)}\n")

        #If invalid entries found, write details to log
        if len(input_invalid) > 0:
            outlog.write(f"{linebreak}Invalid entries{linebreak}")
            for key, value in input_invalid.items():
                outlog.write(f"Line {key}: {value}\n")

        #Write compound mode cummary
        if mode == "compound":
            
            outlog.write(f"{linebreak}Compound results summary{linebreak}")
            outlog.write(f"Compound codes with formulas retrieved: {len(compound_results[0])}\n")
            outlog.write(f"Compound codes without formulas: {len(compound_results[1])}\n\n")

            #Maths!
            outlog.write(f"\nCompound weights summary\n")
            outlog.write(f"Compounds with exact mass: {len(compound_results[2])}\n")
            if len(compound_results[2]) > 0:
                outlog.write(f"Mean compound exact mass: {stat.mean(compound_results[2])}\n")
                outlog.write(f"Median compound exact mass: {stat.median(compound_results[2])}\n")
            else:
                if quiet != False:
                    outlog.write("No compounds detected with exact mass listed\n")
            outlog.write(f"\nCompounds with mol. weight: {len(compound_results[3])}\n")
            if len(compound_results[3]) > 0:
                outlog.write(f"Mean compound mol. weight: {stat.mean(compound_results[3])}\n")
                outlog.write(f"Median compound mol.weight: {stat.median(compound_results[3])}\n")
            else:
                if quiet != False:
                    outlog.write("No compounds detected with molecular weight listed\n")
        #Write KO mode summary
        elif mode == "ko":
            outlog.write(f"{linebreak}KO results summary{linebreak}")
            outlog.write(f"Total unique reactions retrieved: {len(ko_results[0])}\n")
            outlog.write(f"Total unique compounds retrieved: {len(ko_results[1])}\n")
            outlog.write(f"Total unique glycans retrieved: {len(ko_results[2])}\n")
            if homolog == True:
                outlog.write(f"Total unique genes retrieved: {len(ko_results[3])}")
                #outlog.write(f"{linebreak}KO Gene retrieval{linebreak}")
                #outlog.write(f"KO: Genes retrieved\n")
                #for key, value in ko_results[3].items():
                #    outlog.write(f"{key}: {len(value)}\n")
        #Write module mode summary
        elif mode == "module":
            outlog.write((f"{linebreak}Module results summary{linebreak}"))
            outlog.write(f"Total unique compounds retrieved: {len(module_results[0])}\nTotal unique glycans retrieved: {len(module_results[1])}")
        #Write reaction mode summary
        elif mode == "reaction":
            outlog.write((f"{linebreak}Reaction results summary{linebreak}"))
            outlog.write(f"Reaction codes with data retrieved: {len(reaction_results_list[0])}\nReaction codes without data retrieved: {len(reaction_results_list[1])}\n")
        elif mode == "mdata":
            pass


    
#------------------------------------------------------------Run program---------------------------------------------------------------#

if mode == "compound":
    compound_results = compound_data()
elif mode == "ko":
    ko_results = ko_to_compound()
elif mode == "module":
    module_results = module_to_compound()
elif mode == "reaction":
    reaction_results_list = reaction_data()
elif mode == 'mdata':
    module_data()
if mode != "mdata":
    write_log()

print(f"{linebreak}Thank you for using {appname}. More details regarding {appname} can be found at {github}.",
      f"{appname} is neither endorsed by, nor associated with KEGG. Please cite the relevant KEGG literature:",
      "Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000). https://doi.org/10.1093/nar/28.1.27",
      "Kanehisa, M; Toward understanding the origin and evolution of cellular organisms. Protein Sci. 28, 1947-1951 (2019). https://doi.org/10.1002/pro.3715",
      "Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M. and Ishiguro-Watanabe, M.; KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Res. 51, D587-D592 (2023). https://doi.org/10.1093/nar/gkac963",
      sep="\n")
if pubchemarg == True:
    print(f"\nThis program can retrieve a subset of data from the PubChem database using PUG-REST via the PubChemPy library. Please cite the relevanet PubChem literature if publishing.",
        f"Kim S, Chen J, Cheng T, et al. PubChem 2023 update. Nucleic Acids Res. 2023;51(D1):D1373-D1380. doi:10.1093/nar/gkac956", 
        f"Kim S, Thiessen PA, Cheng T, Yu B, Bolton EE. An update on PUG-REST: RESTful interface for programmatic access to PubChem. Nucleic Acids Res. 2018 July 2; 46(W1):W563-570. [PubMed PMID: 29718389] doi: 10.1093/nar/gky294.",
        "PubChemPy can be found at https://pubchempy.readthedocs.io/en/latest", 
        sep="\n")
print(f"{linebreak}")

#Written by Lachlan McKinnie


