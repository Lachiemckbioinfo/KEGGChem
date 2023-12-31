#!/usr/bin/env python3
#version 0.9

import requests
import re
#import sys
import os
from collections import Counter
import argparse
import statistics as stat
import random
import datetime
import pubchempy as pcp
import csv

github = "https://github.com/Lachiemckbioinfo/KEGGChem"
citation = "XXXXX"
#KEGGcitation = (f"KEGG citation goes here")
appname = 'KEGGChem'
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
                    

parser.add_argument("-f", "--file",
                    required = True,
                    help = "The input file",
                    type = argparse.FileType('r'),
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
        outdir = str(infile.name + "_output")
        outdir = check_outdir(outdir)
    return outdir
outdir = set_outdir()

            
            

        

#Set downloads directory
if args.download is not None:
    dir_download = str(args.download)
else:
    dir_download = "keggchem_downloads"
dirs = ["KEGG_entries/compounds", "KEGG_entries/orthologues", "KEGG_entries/reactions",
        "KEGG_entries/modules", "Structure_files/Downloaded", "Structure_files/nullfiles", "SDFfiles"]
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

#Extract input codes and appent to input list
def openfile(x):
    with infile as file:
        while (line := file.readline().strip()):
            line = line.upper()
            if "KO:" in line:
                line = line.replace("KO:", "")
            #If search == [LETTER] or RC + 5 digits, proceed
            if mode != "download":
                if re.search(rf"{x}\d{{5}}\b", line):
                #if re.search(r"\b([MKTCGRNHD]|RC)\d{5}\b", line):
                    input.append(line)
                else:
                    if quiet != False:
                        print(f"Invalid search term: {line}\n")
            else:
                if re.search(r"\b([MKTCGRNHD]|RC)\d{5}\b", line):
                    input.append(line)
                else:
                    if quiet != False:
                        print(f"Invalid search term {line}")


#Create output directory structure
def create_outdirs(dir_list):
    for item in dir_list:
        path = os.path.join(outdir, item)
        os.makedirs(path, exist_ok = True)
    if verbose == True:
        print(f"Create output directories {dir_list} in {outdir}")

#----------------------------------------Process input----------------------------------------#
#openfile()
if mode == "ko":
    openfile("K")
elif mode == "compound":
    openfile("C")
elif mode == "reaction":
    openfile("R")
elif mode == "module":
    openfile("M")
input_unique = sorted([*set(input)])
input_total = len(input)
input_unique_total =len(input_unique)
total_reactions_list = []
total_modules_list = []
total_compounds_list = []
total_glycans_list = []
dict_pathways = {}


#Pathway summary writer
#Write KEGG pathway summary
def write_pathway_summary(filename):
    pathway_writer = csv.writer(filename, delimiter='\t')
    pathway_writer.writerow(["Pathway", "Name", "Count"])
    pathway_data = []
    for key, value in dict_pathways.items():
        pathway_data.append([key, value[0], value[1]])
    pathway_writer.writerows(pathway_data)


#----------------------------------------Specific mode functions----------------------------------------#

#--------------------KEGG KO/Module to compound--------------------
# Function to extract reaction codes from KEGG orthologue pages
def get_reaction_codes(KO):
    ko_file = os.path.join(dir_download, "KEGG_entries/orthologues", KO)
    #Check if file exists in download folder
    if os.path.exists(ko_file) == False:
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




# Function to retrieve compound/glycan codes associated with a reaction or module code
def get_compound_codes(reaction_code):
    reaction_file = os.path.join(dir_download, "KEGG_entries/reactions", reaction_code)
    if os.path.exists(reaction_file) == False:
        url = f"http://rest.kegg.jp/get/{reaction_code}"
        req = requests.get(url).text
        #Strip everything in req from "GENES" onwards
        if mode == "ko":
            r = req.split(separator, 1)[0]
        else:
            if "COMPOUND" in req:
                r = req.split("COMPOUND", 1)[1]
            else:
                r = req
        with open(reaction_file, "w") as file:
            file.write(req)
        if verbose == True:
            print(f"Wrote {reaction_code} to file\n")
    else:
        with open(reaction_file, "r") as file:
            if mode == "ko":
                r = file.read().split(separator, 1)[0]
            else:
                req = file.read()
                if "COMPOUND" in req:
                    r = req.split("COMPOUND", 1)[1]
                else:
                    r = req
        if verbose == True:
            print(f"Extracted {reaction_code} from file\n")
    #Extract compound codes from page
    compound_codes = [*set(re.findall(r"\b[CG]\d{5}\b", r))]
    return compound_codes


#Function to retrieve compound data from compound codes
def get_compound_data(compound_code):
    compound_file = os.path.join(dir_download, "KEGG_entries/compounds", compound_code)
    if os.path.exists(compound_file) == False:
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
    if os.path.exists(compound_file) == False:
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
    if os.path.exists(reaction_file) == False:
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
    return equation, definition, rclass, ko_codes, pathways


#---------------------Structure download functions----------------------#
def download_structure_files(compound_code, mode):
    filename = f"{compound_code}.{mode}"
    filedir = os.path.join(dir_download, "Structure_files/Downloaded", filename)
    nulldir = os.path.join(dir_download, "Structure_files/nullfiles", filename)
    if mode == "mol":
        urlcode = "-f+m+"
    elif mode == "kcf":
        urlcode = "-f+k+"

    if os.path.exists(filedir) == False and os.path.exists(nulldir) == False:
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
#----------------------------------------Function - Run KO to compound----------------------------------------#
def ko_to_compound():
    #Create outdir structure
    dir_list = ["Results/Individual_results", "Summaries"]
    create_outdirs(dir_list)
    count_current = 1
    dict_reaction_codes = {}
    dict_module_codes = {}
    module_compounds_list = []
    reaction_compounds_list = []
    module_glycan_list = []
    reaction_glycan_list = []

    print(f"{linebreak}Analysing {input_unique_total} unique KEGG orthologues from ",\
        f"{infile.name} ({input_total} given){linebreak}", sep="")
    
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
    with open(os.path.join(outdir, "run_summary.txt"), "a") as summary,\
    open(os.path.join(outdir, "Summaries", "compound_summary.txt"), "a") as sum_compounds,\
    open(os.path.join(outdir, "Summaries", "reaction_summary.txt"), "a") as sum_reactions,\
        open(os.path.join(outdir, "Summaries", "glycan_summary.txt"), "w") as sum_glycans:
        
        summary.write(f"{linebreak}Overall summary:{linebreak}")
        filepath = os.path.abspath(outdir)
        summary.write(f"Input file: {infile.name}\n")
        summary.write(f"Results saved to {filepath}\n")
        runtime, timestamp_end = timetaken()
        summary.write(f"Start time: {timestamp}\n")
        summary.write(f"End time: {timestamp_end}\n")
        summary.write(f"Runtime: {runtime}\n")
        summary.write(f"Total number of orthologues: {input_total}\n")
        summary.write(f"Total number of unique orthologues: {input_unique_total}\n")
        summary.write(f"Total unique reactions: {reactions_count}\n")
        summary.write(f"Total unique compounds: {compounds_count}\n")
        summary.write(f"Total unique glycans count: {glycans_count}\n")
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
        
        #Write compound/glycan summary results to file
        write_summaries_to_file("Compound", sum_compounds, compounds_count, module_compounds_unique,
                                reaction_compounds_unique, compound_module_count, compound_reaction_count,
                                compound_common_count, total_compounds_dict)
        write_summaries_to_file("Glycan", sum_glycans, glycans_count, module_glycans_unique, 
                                reaction_glycans_unique, glycan_module_count, glycan_reaction_count, 
                                glycan_common_count, total_glycans_dict)


    



#----------------------------------------Function - Run module to compound----------------------------------------#
def module_to_compound():
    count_current = 1
    #Create outdir structure
    dir_list = ["Results", "Summaries"]
    create_outdirs(dir_list)

    print(f"{linebreak}Analysing {input_unique_total} unique KEGG modules from ",\
        f"{infile.name} ({input_total} given){linebreak}", sep="")
    
    for module in input_unique:
        print(f"Extracting compounds from KEGG module {module}\n")
        compound_codes = get_compound_codes(module)
        outfile = f"{module}.txt"

        with open(os.path.join(outdir, "Results", outfile),"a") as result,\
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
    with open(os.path.join(outdir, "run_summary.txt"), "w") as summary,\
        open(os.path.join(outdir, "Summaries", "compound_summary.txt"), "w") as sum_compound,\
            open(os.path.join(outdir, "Summaries", "glycan_summary.txt"), "w") as sum_glycan:
        summary.write(f"{linebreak}Overall summary:{linebreak}")
        filepath = os.path.abspath(outdir)
        summary.write(f"Input file: {infile.name}\n")
        summary.write(f"Results saved to {filepath}\n")
        runtime, timestamp_end = timetaken()
        summary.write(f"Start time: {timestamp}\n")
        summary.write(f"End time: {timestamp_end}\n")
        summary.write(f"Runtime: {runtime}\n")
        summary.write(f"Total number of modules {input_total}\n")
        summary.write(f"Total number of unique modules: {input_unique_total}\n")
        summary.write(f"Total unique compounds: {compounds_count}\n")
        summary.write(f"Total number of unique glycans {glycans_count}\n")
        #Write compound summary results to file
        sum_compound.write(f"{linebreak}Compounds summary:{linebreak}")
        sum_compound.write(f"Total unique compounds: {compounds_count}\n\n")
        sum_compound.write("Compound: count\n")
        for key, value in total_compounds_dict.items():
            sum_compound.write(f"{key}: {value}\n")
        #Write gycan summary results to file
        sum_glycan.write(f"{linebreak}Glycan summary:{linebreak}")
        sum_glycan.write(f"Total unique glycans: {glycans_count}\n\n")
        sum_glycan.write("Glycan: count\n")
        for key, value in total_glycans_dict.items():
            sum_glycan.write(f"{key}: {value}\n")
    
    


#----------------------------------------Function - Run compound data retrieve----------------------------------------
def compound_data():
    count_current = 1
    #Create outdir
    dir_list = ["Results", "Summaries"]
    if sdfarg == True:
        dir_list.append("Results/SDF")
    create_outdirs(dir_list)

    compound_names_list = []
    compound_names_null = []
    lst_exact_mass = []
    lst_mol_weight = []
    lst_chebi = []
    lst_pubchemsid = []
    lst_pathways = []

    print(f"{linebreak}Retrieving data for {input_unique_total} unique KEGG compounds from ",\
                    f"{infile.name} ({input_total} given){linebreak}", sep="")

    with open(os.path.join(outdir, "Results", "Formulas.txt"), "a") as formulas,\
    open(os.path.join(outdir, "Results", "Weights.txt"), "a") as weight,\
        open(os.path.join(outdir, "Results", "Compound_data.tsv"), "a") as out:
        out_csv = []
        out_write = csv.writer(out, delimiter='\t')
        #out.write("Compound;formula;exact mass;molecular weight;ChEBI ID;PubChem SID;pathways\n")
        out_header = ["compound", "formula", "exact mass", "molecular weight", "ChEBI ID", "PubChem SID", "KEGG pathways"]
        formulas.write("Compound:formula\n")
        weight.write("compound;exact_mass;mol_weight\n")
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
                formulas.write(f"{compound}:{compound_formula}\n")
                compound_names_list.append(compound)
                if quiet == False:
                    print(f"Data retrieved for {compound}: {compound_formula} ({count_current}/{input_unique_total})\n")
            else:
                formulas.write(f"{compound}:NULL\n")
                compound_names_null.append(compound)
                if quiet == False:
                    print(f"Data retrieved for {compound}: NULL ({count_current}/{input_unique_total})\n")
            #Write weight results
            weight.write(f"{compound};{exact_mass};{mol_weight}\n")
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
            SIDout.write("SID;SID source ID;CID\n")
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
                    SIDout.write(f"{sid};{sid_source_id};{sid_CID}\n")
                    sidcount += 1
                except:
                    print(f"Invalid SID code searched: {sid}")
        with open(os.path.join(outdir, "Results", "PubChem_compounds.txt"), "a") as CIDout:
            print(f"{linebreak}Retrieving PubChem CID data from PubChem Compound database{linebreak}")
            if verbose == True:
                print(f"Opened PubChem Compound results file at {CIDout}")
            CIDout.write("CID;Molecular formula;Canonical SMILES;Isomeric SMILES\n")
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
                    CIDout.write(f"{cid};{cid_formula};{cid_canon_SMILES};{cid_isomeric_SMILES}\n")
                    if quiet == False:
                        print(f"Retrieved compound data for CID: {cid} from PubChem Compound database ({cidcount}/{cidcounttotal})")
                    cidcount += 1
                except:
                    print(f"Invalid CID code searched: {cid}")
        print(f"{linebreak}Completed retrieving SID and CID data from PubChem{linebreak}")

    #Download structure files
    structure_file_download(input_unique)

    #Write summary
    with open(os.path.join(outdir, "Summaries", "run_summary.txt"), "a") as summary,\
        open(os.path.join(outdir, "Summaries", "Pathway_summary.tsv"), "a") as path_summary:
        if verbose == True:
            print(f"\nWriting summaries to {os.path.join(outdir, 'Summaries')}\n")
        summary.write(f"{linebreak}KEGG compound data retrieval summary:{linebreak}")
        summary.write("\nGeneral information\n")
        filepath = os.path.abspath(outdir)
        summary.write(f"Input file: {infile.name}\n")
        summary.write(f"Results saved to {filepath}\n")
        runtime, timestamp_end = timetaken()
        summary.write(f"Start time: {timestamp}\n")
        summary.write(f"End time: {timestamp_end}\n")
        summary.write(f"Runtime: {runtime}\n")
        summary.write(f"Compound codes given: {input_total}\n")
        summary.write(f"Unique compound codes given: {input_unique_total}\n")
        summary.write(f"Output directiory: {outdir}\n")
        summary.write(f"\nCompound formula summary\n")
        summary.write(f"Compound codes with formulas retrieved: {len(compound_names_list)}\n")
        summary.write(f"Compound codes without formulas: {len(compound_names_null)}\n\n")

        #Maths!
        summary.write(f"\nCompound weights summary\n")
        summary.write(f"Compounds with exact mass: {len(lst_exact_mass)}\n")
        if len(lst_exact_mass) > 0:
            summary.write(f"Mean compound exact mass: {stat.mean(lst_exact_mass)}\n")
            summary.write(f"Median compound exact mass: {stat.median(lst_exact_mass)}\n")
        else:
            if quiet != False:
                print("No compounds detected with exact mass listed")
        summary.write(f"Compounds with mol. weight: {len(lst_mol_weight)}\n")
        if len(lst_mol_weight) > 0:
            summary.write(f"Mean compound mol. weight: {stat.mean(lst_mol_weight)}\n")
            summary.write(f"Median compound mol.weight: {stat.median(lst_mol_weight)}\n")
        else:
            if quiet != False:
                print("No compounds detected with molecular weight listed")
        write_pathway_summary(path_summary)
        
        
#----------------------------------------Function - Run reaction data retrieve----------------------------------------
def reaction_data():
    count_current = 1
    #Create outdir
    dir_list = ["Results", "Summaries"]
    create_outdirs(dir_list)

    reaction_names_list = []
    reaction_names_null = []
    lst_rclass = []
    lst_ko = []
    lst_pathways = []
    lst_pathway_codes = []

    print(f"{linebreak}Retrieving data for {input_unique_total} unique KEGG reactions from ",\
                    f"{infile.name} ({input_total} given){linebreak}", sep="")

    with open(os.path.join(outdir, "Results", "Equations.txt"), "w") as f_equation,\
        open(os.path.join(outdir, "Results", "Definitions.txt"), "w") as f_definition,\
            open(os.path.join(outdir, "Results", "Rclass.txt"), "w") as f_rclass,\
                open(os.path.join(outdir, "Results", "KO_codes.txt"), "w") as f_ko,\
                    open(os.path.join(outdir, "Results", "Pathways.txt"), "w") as f_pathway:
        f_equation.write("Reaction:equation\n")
        f_definition.write("Reaction:definition\n")
        f_rclass.write("Reaction:reaction classes\n")
        f_ko.write("Reaction:KO codes\n")
        f_pathway.write("Reaction:Pathways\n")
        
        for reaction in input_unique:
            #Write equation results
            equation, definition, rclasses, ko_codes, pathways = get_reaction_data(reaction)
            if len(equation) > 0:
                f_equation.write(f"{reaction}:{equation}\n")
            else:
                f_equation.write(f"{reaction}:NULL\n")
            #Write definition results
            if len(definition) > 0:
                f_definition.write(f"{reaction}:{definition}\n")
            else:
                f_definition.write(f"{reaction}:NULL\n")
            #Write RCLASS results
            if len(rclasses) > 0:
                for rclass in rclasses:
                    lst_rclass.append(rclass)
                f_rclass.write(f"{reaction}:{', '.join(rclasses)}\n")
            else:
                f_rclass.write(f"{reaction}:NULL\n")
            #Write KO codes
            if len(ko_codes) > 0:
                for ko_code in ko_codes:
                    lst_ko.append(ko_code)
                f_ko.write(f"{reaction}:{','.join(ko_codes)}\n")
            else:
                f_ko.write(f"{reaction}:NULL\n")
            #Write pathways
            if len(pathways) > 0:
                pathway_codes = [re.sub("  .*", "", item) for item in pathways]
                for pathway_code in pathway_codes:
                    lst_pathway_codes.append(pathway_code)
                f_pathway.write(f"{reaction}:{','.join(pathway_codes)}\n")
                for pathway in pathways:
                    lst_pathways.append(pathway)
            else:
                f_pathway.write(f"{reaction}:NULL\n")
                
            if quiet == False:
                print(f"Data retrieved for {reaction}: ({count_current}/{input_unique_total})\n")
            count_current += 1
        print(f"{linebreak}Finished retrieving reaction data{linebreak}")

    #Download structure files
    if structure == True:
        print("Structure retrieval is not run with KEGG reaction retrieve")

    #Write summary
    with open(os.path.join(outdir, "Summaries", "run_summary.txt"), "a") as summary:
        summary.write(f"{linebreak}KEGG reaction data retrieval summary:{linebreak}")
        summary.write("\nGeneral information\n")
        filepath = os.path.abspath(outdir)
        summary.write(f"Input file: {infile.name}\n")
        summary.write(f"Results saved to {filepath}\n")
        summary.write(f"File input: {infile.name}\n")
        runtime, timestamp_end = timetaken()
        summary.write(f"Start time: {timestamp}\n")
        summary.write(f"End time: {timestamp_end}\n")
        summary.write(f"Runtime: {runtime}\n")
        summary.write(f"Reaction codes given: {input_total}\n")
        summary.write(f"Unique reaction codes given: {input_unique_total}\n")
        summary.write(f"Output directory: {outdir}\n")
        summary.write(f"\nReaction formula summary\n")
        summary.write(f"Reaction codes with data retrieved: {len(reaction_names_list)}\n")
        summary.write(f"Reaction codes without data retrieved: {len(reaction_names_null)}\n\n")

        
#------------------------------------------------------------Run module data retrieve---------------------------------------------------------------#
#Function for retrieving module data
def get_module_data(module_code):
    module_file = os.path.join(dir_download, "KEGG_entries/modules", module_code)
    if os.path.exists(module_file) == False:
        url = f"http://rest.kegg.jp/get/{module_code}"
        req = requests.get(url).text
    else:
        with open(module_file, "r") as file:
            req = file.read()
    entry = re.search(r"ENTRY       .*", req)
    moduletype = entry.group(0).replace("ENTRY       ", "")
    moduletype = moduletype.replace("            ",":")
    name = re.search(r"NAME        .*", req)
    modulename = name.group(0).replace("NAME        ","")
    print(f"{moduletype}:{modulename}")

def module_data():
    #count_current = 1
    #Create outdir
    #dir_list = ["Results", "Summaries"]
    #create_outdirs(dir_list)



    print(f"{linebreak}Retrieving data for {input_unique_total} unique KEGG modules from ",\
                    f"{infile.name} ({input_total} given){linebreak}", sep="")
    for module in input_unique:
        get_module_data(module)
        #print(f"{module}: {moduletype}: {name}")
        #count_current += 1
        #print(f"{linebreak}Finished retrieving module data {module}{linebreak}")

    #Download structure files
    if structure == True:
        print("Structure retrieval is not run with KEGG module retrieve")


 

    
#------------------------------------------------------------Run program---------------------------------------------------------------#

if mode == "compound":
    compound_data()
elif mode == "ko":
    ko_to_compound()
elif mode == "module":
    module_to_compound()
elif mode == "reaction":
    reaction_data()
elif mode == 'mdata':
    module_data()

print(f"{linebreak}Thank you for using {appname}. More details regarding {appname} can be found at {github}.",
      f"{appname} is neither endorsed by, nor associated with KEGG. Please cite the relevant KEGG literature:",
      "Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000). https://doi.org/10.1093/nar/28.1.27",
      "Kanehisa, M; Toward understanding the origin and evolution of cellular organisms. Protein Sci. 28, 1947-1951 (2019). https://doi.org/10.1002/pro.3715",
      f"Kanehisa, M., Furumichi, M., Sato, Y., Kawashima, M. and Ishiguro-Watanabe, M.; KEGG for taxonomy-based analysis of pathways and genomes. Nucleic Acids Res. 51, D587-D592 (2023). https://doi.org/10.1093/nar/gkac963",
      sep="\n")
if pubchemarg == True:
    print(f"\nThis program can retrieve a subset of data from the PubChem database using PUG-REST via the PubChemPy library. Please cite the relevanet PubChem literature if publishing.",
        f"Kim S, Chen J, Cheng T, et al. PubChem 2023 update. Nucleic Acids Res. 2023;51(D1):D1373-D1380. doi:10.1093/nar/gkac956", 
        f"Kim S, Thiessen PA, Cheng T, Yu B, Bolton EE. An update on PUG-REST: RESTful interface for programmatic access to PubChem. Nucleic Acids Res. 2018 July 2; 46(W1):W563-570. [PubMed PMID: 29718389] doi: 10.1093/nar/gky294.",
        "PubChemPy can be found at https://pubchempy.readthedocs.io/en/latest", 
        sep="\n")
print(f"{linebreak}")

#Written by Lachlan McKinnie


