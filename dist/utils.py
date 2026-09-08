"""
MetaMiner - Streamlined GUI Tool for Retrieving, Normalizing and Exploring Metadata
Copyright (C) [2025] [Patel Jaykumar Kiritkumar, 
Molecular Imaging and Molecular Diagnostics Lab, 
Indian Institute of Technology Delhi, India]

Licensed under the GNU Affero General Public License v3.0 (AGPL-3.0)
with additional restrictions: **Non-Commercial Use Only**.

For full terms, see the LICENSE and LICENSE-EXCEPTIONS.md files.
For commercial licensing inquiries, contact: prekijpatel2.0@gmail.com
"""

import json
from tkinter import messagebox
import subprocess
import socket
import pandas as pd
import pycountry as pc
from datetime import datetime
from geopy.geocoders import Nominatim
import re
import logging
import time
import unicodedata as uni
import numpy as np
from rapidfuzz import fuzz
import os
from dash import dcc, html, Dash, Input, Output, ctx
import dash_bootstrap_components as dbc
from dash_graphs import Choropleth_map, Assembly_level_bar, Annotation_bar, Submission_year_line, Sequencing_technologies_scatter, Coverage_bar, ANI_scatter, Total_genes_hist, CDSs_hist, Non_coding_hist, Pseudogenes_hist, Isolation_source_treemap, N50L50_scatter
import dash_daq as daq
import os
import webbrowser
import sys


logging.basicConfig(filename="log_file.log", level=logging.DEBUG, format='%(asctime)s:%(levelname)s - %(message)s')

def fix_gb_rf_annotation(df:pd.DataFrame):
    annotation_df =  df.copy()
    annotation_df[['prefix','suffix']] = annotation_df['Accession_ID'].str.split('_', expand=True)
    annotation_df = annotation_df.sort_values(by='prefix', ascending=False)
    annotation_df = annotation_df.drop_duplicates(subset='suffix', keep='first')

    condition = (
    annotation_df['Accession_ID'].str.startswith('GCA') &
    annotation_df['Paired_assembly_accession'].str.startswith('GCF') &
    (annotation_df['Annotation_from'] == 'NCBI RefSeq')
    )

    annotation_df.loc[condition, 'Annotation_from'] = 'GenBank'

    return annotation_df

def is_internet_there(host="8.8.8.8", port=53, timeout=3):
    """
    Check for an active internet connection by trying to connect to a public DNS server.
    """

    logging.info("Checking for active internet connection.")
    
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        
        logging.debug("Active internet connection found.")
        return True

    except socket.error:
        logging.error("No active internet connection found.")
        messagebox.showerror('No active internet connection!', 'Please check your internet connection.')
        return False

def pass_to_cmd(command):
    try:
        logging.debug(f"Executing command in pass_to_cmd function: {command}")
        process = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        return True, process.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred while executing command: {e}")
        self.end_download_logging = True
        return False, e.stderr

def load_json(filepath:str):
    """
    Loads JSON data from a file. It takes any of the two-types of multi-jsonl files(provided by NCBI) as input and convert them to uniform output. 
    """
    
    try:
        # Attempt to read the second-type of multi-JSON file
        with open(filepath, 'r', encoding='utf-8') as json_file:
            json_content = json_file.read()
            logging.debug(f"Successfully read the file: {filepath} with utf-8.")
    
    except UnicodeDecodeError as e:
        with open(filepath, 'r', encoding='latin1') as json_file:
            json_content = json_file.read()
            logging.debug(f"Successfully read the file: {filepath} with latin1 encoding.")
            
    except UnicodeDecodeError as e:
        with open(filepath, 'r', encoding='utf-16') as json_file:
            json_content = json_file.read()
            logging.debug(f"Successfully read the file: {filepath} with utf-16 encoding.")

    except UnicodeDecodeError as e:
        messagebox.showerror('Unicode decode error!', f'The JSON file is not in expected format. {e}')
        logging.error(f"The JSON file is not in expected format: {e}")
        return None

    except IOError as IOe:
        messagebox.showerror('Unable to read the file!', IOe)
        logging.error(f"Error reading the JSON file: {IOe}")
        return None
        
    except Exception as e:
        messagebox.showerror('An Unexpected error occurred!', e)
        logging.error(f"An unexpected error occurred while reading the JSON file. Make sure the file is json file.: {e}")
        return None

    if 'total_count' in json_content: # loading the read `json_content` of the second-type `multi-json` directly into `data` variable
        
        try:
            data = json.loads(json_content)
            logging.debug("Second-type multi-JSON file detected.")
            logging.info("Successfully loaded the multi-JSON file.")
            logging.debug(f"Data type of 'data' variable: {type(data)}")
            logging.debug(f"Total records in the multi-JSON file: {data['total_count']}")
            return data
            
        except json.JSONDecodeError as e:
            messagebox.showerror('JSON decode error!', f'The JSON file is not in expected format. {e}')
            logging.error(f"The JSON file is not in expected format: {e}")
            return None
    
    else: # loading the read `json_content` of the first-type `multi-json` (aka json-lines) into `data` variable with some addition of keys, as it requires to add some extra keys to match with the output of second-type `multi-json`.
        
        try:
            # Attempt to read the first-type of multi-JSON file
            with open(filepath, 'r') as json_file:
                json_content = json_file.readlines()
                logging.debug("First-type multi-JSON file detected.")
                logging.debug(f"Successfully read the JSON lines.")
        
        except UnicodeDecodeError as e:
            with open(filepath, 'r', encoding='utf-8') as json_file:
                json_content = json_file.readlines()
                logging.debug("First-type multi-JSON file detected.")
                logging.debug(f"Successfully read the JSON lines with utf-8 encoding.")
        
        except UnicodeDecodeError as e:
            with open(filepath, 'r', encoding='latin1') as json_file:
                json_content = json_file.readlines()
                logging.debug("First-type multi-JSON file detected.")
                logging.debug(f"Successfully read the JSON lines with latin1 encoding.")
            
        except UnicodeDecodeError as e:
            with open(filepath, 'r', encoding='utf-16') as json_file:
                json_content = json_file.readlines()
                logging.debug("First-type multi-JSON file detected.")
                logging.debug(f"Successfully read the JSON lines with utf-16 encoding.")

        except IOError as IOe:
            messagebox.showerror('Unable to read the file!', IOe)
            logging.error(f"Error reading the JSON file: {IOe}")
            return None
        
        except Exception as e:
            messagebox.showerror('An Unexpected error occurred!', e)
            logging.error(f"An unexpected error occurred while reading the JSON file: {e}")
            return None
 
        try:
            reports = []
            for line in json_content:
                line = line.strip()
                # print(line)
                if line: 
                    report = json.loads(line)
                    # print(report)
                    reports.append(report)
                    # print(reports)

            logging.debug("Successfully read the JSON lines into 'reports' variable.")
            logging.debug(f"The total lines loaded in the `reports` variable: {len(reports)}")
            

            data = {
            'reports': reports,
            'total_count': len(json_content) 
            }

            logging.info("Successfully loaded the multi-JSON file.")
            logging.debug(f"Data type of 'data' variable: {type(data)}")
            logging.debug(f"Total records in the multi-JSON file: {data['total_count']}. This number is from load_json() function.")
            
            return data
        
        except json.JSONDecodeError as e:
            messagebox.showerror('JSON decode error!', f'The JSON file is not in expected format. {e}')
            logging.error(f"The JSON file is not in expected format: {e}")
            return None

def changing_keys_to_camelcase(key:str):
    """
    Changes the snakecase keys to camelcase keys. For example, 'assembly_info' to 'assemblyInfo'.
    """
    try:
        parts = key.split('_')
        camelcased_key = parts[0] + ''.join(word.capitalize() for word in parts[1:])
        return camelcased_key
    except:
        return key

def get_value(basic_string, nested_key_chain, when_key_not_found=''):
    """
    Gets values from loaded multi-JSON (dictionary) gracefully and helps in avoiding errors,
    especially when data is not uniform for each genome. 

    - basic_string: all the data within the multi-JSON file is in `['reports'][i]` section. Where `i` is ith genome entry in multi-JSON file. 
    - nested_key_chain: List of nested entries one after another from which we need information. Helps in avoiding errors when any intermediate key is not found in JSON file (dictionary).
    - when_key_not_found: Default option is `''` aka - empty string. When entry for specific key is not found in JSON file, this is returned.  
    
    """
    
    value = basic_string
    
    try:      
        for key in nested_key_chain:
            try:
                if isinstance(key, int):
                    if key < len(value):
                        value = value[key]
                    else:
                        return when_key_not_found
                else:
                    value = value[key]
                
            except (KeyError, IndexError):
                try:
                    value = value[changing_keys_to_camelcase(key)]
                except (KeyError, IndexError):
                    return when_key_not_found
        return value
            
    except (KeyError, TypeError, IndexError):
        return when_key_not_found

def record_numbers(data:dict):
    """
    Extracts the value for total number of genomes in the loaded multi-JSON file.
    """
    
    if isinstance(data, dict) and 'total_count' in data:
        total_count = data['total_count']
        logging.debug(f"Total records in the multi-JSON file: {total_count}. This should match one with the load_json() function.")
        return total_count
    else:
        messagebox.showerror('Error!', "The loaded file may not be multi-JSON file.")
        return 0

def strip_gca_gcf_id(id):
    try:
        if 'GCA' in id:
            return id.replace('GCA_', '')
        else:
            return id.replace('GCF_', '')
    except: # in cases where the given `id` is not string.
        return ''

def biosample_attributes_from_json(data:dict, total_count:int):
    """
    Extracts various attributes of Biosample from loaded Multi-JSON file with the help of `get_value` function.
    """
    if data:
        logging.info("Extracting Biosample attributes from JSON file.")

    BioSampleData = []

    if total_count > 0:
        for i in range(total_count):
            biosample_metadata = {
                'Accession_ID': get_value(data['reports'][i], ['accession']),
                'Organism': get_value(data['reports'][i], ['organism', 'organism_name']),
                'strain': get_value(data['reports'][i], ['organism', 'infraspecific_names', 'strain']),
                'Common_name': get_value(data['reports'][i], ['organism', 'common_name']),
                'Taxonomy_ID': get_value(data['reports'][i], ['organism', 'tax_id']),
                'Sequencing_Technology' : get_value(data['reports'][i], ['assembly_info', 'sequencing_tech']),
                'Owner/Submitter': get_value(data['reports'][i], ['assembly_info', 'submitter']),
                'Bioproject_ID': get_value(data['reports'][i], ['assembly_info', 'bioproject_accession']),
                'Bioproject_parentID':get_value(data['reports'][i], ['assembly_info', 'bioproject_lineage', 0, 'bioprojects', 0, 'parent_accessions', 0]),
                'Bioproject_title': get_value(data['reports'][i], ['assembly_info', 'bioproject_lineage', 0, 'bioprojects', 0, 'title']),
                'Biosample_ID': get_value(data['reports'][i], ['assembly_info', 'biosample', 'accession']),
                'Biosample_title': get_value(data['reports'][i], ['assembly_info', 'biosample', 'description', 'title']),

                # not keeping the following things as they are not present in all the genomes and will create a lot of empty columns in the final dataframe. Plus, if any of them are there, they usually are also there in the 'attributes' section of biosample. So, we will extract them from there anyways.

                # # since the first one would return empty string (falsy value in python), we will try to get the value from the second one.
                # 'Biosample_breed': get_value(data['reports'][i], ['assembly_info', 'biosample', 'breed']) or get_value(data['reports'][i], ['organism', 'infraspecific_names', 'breed']),
                # 'Biosample_cultivar': get_value(data['reports'][i], ['assembly_info', 'biosample', 'cultivar']) or get_value(data['reports'][i], ['organism', 'infraspecific_names', 'cultivar']),
                # 'Biosample_ecotype': get_value(data['reports'][i], ['assembly_info', 'biosample', 'ecotype']) or get_value(data['reports'][i], ['organism', 'infraspecific_names', 'ecotype']),
                # 'Biosample_isolate': get_value(data['reports'][i], ['assembly_info', 'biosample', 'isolate']) or get_value(data['reports'][i], ['organism', 'infraspecific_names', 'isolate']),

                # 'Biosample_sex': get_value(data['reports'][i], ['assembly_info', 'biosample', 'sex']),
                # 'Biosample_developmetal_stage': get_value(data['reports'][i], ['assembly_info', 'biosample', 'dev_stage']),
                # 'Biosample_tissue': get_value(data['reports'][i], ['assembly_info', 'biosample', 'tissue']),
                # 'Biosample_biomaterial_provider': get_value(data['reports'][i], ['assembly_info', 'biosample', 'biomaterial_provider']),
                # 'Collected_by': get_value(data['reports'][i], ['assembly_info', 'biosample', 'collected_by']),
                # 'Collection_date': get_value(data['reports'][i], ['assembly_info', 'biosample', 'collection_date']),
                # 'Host': get_value(data['reports'][i], ['assembly_info', 'biosample', 'host']),
                # 'Geo_loc_name': get_value(data['reports'][i], ['assembly_info', 'biosample', 'geo_loc_name']),

                'Model': get_value(data['reports'][i], ['assembly_info', 'biosample', 'models', 0]),
                'IFSAC_category': get_value(data['reports'][i], ['assembly_info', 'biosample', 'ifsac_category']),
                'Owner': get_value(data['reports'][i], ['assembly_info', 'biosample', 'owner', 'name']),
                'Owner_contact': get_value(data['reports'][i], ['assembly_info', 'biosample', 'owner', 'contacts', 0]),
                'Submission_date': get_value(data['reports'][i], ['assembly_info', 'biosample', 'submission_date']),
                'Last_update_date': get_value(data['reports'][i], ['assembly_info', 'biosample', 'last_update_date']),
                'Release_date': get_value(data['reports'][i], ['assembly_info', 'biosample', 'release_date']),
            }
            
            # since the attributes are really random and vary a lot from genome to genome, we will try to take all the possible attributes from all the assemblies   
            biosample_attributes = get_value(data['reports'][i], ['assembly_info', 'biosample', 'attributes'])
            # print(biosample_attributes)
            for attribute in biosample_attributes:
                # print(attribute)
                if 'name' in attribute and 'value' in attribute:
                    biosample_metadata.update({attribute['name']: attribute['value']})
                else:
                    logging.warning(f"Attribute missing 'name' or 'value' key: {attribute} in genome: {biosample_metadata['Accession_ID']}")


            BioSampleData.append(biosample_metadata)
            # print(biosample_metadata)

        if len(BioSampleData) == total_count:
            logging.debug(f"Biosample attributes extracted successfully. len(BioSampleData) = {len(BioSampleData)}")

    return BioSampleData

def assembly_attributes_from_json(data:dict, total_count:int):
    """
    Extracts various assembly attributes from loaded Multi-JSON file with the help of `get_value` function.
    """
    if data:
        logging.info("Extracting Assembly attributes from JSON file.")

    AssemblyData = []

    if total_count > 0:
        for i in range(total_count):
            Assembly_metadata = {
                'Accession_ID': get_value(data['reports'][i], ['accession']),
                'Organism': get_value(data['reports'][i], ['organism', 'organism_name']),
                'strain': get_value(data['reports'][i], ['organism', 'infraspecific_names', 'strain']),
                'Taxonomy_ID': get_value(data['reports'][i], ['organism', 'tax_id']),
                'Sequencing_Technology' : get_value(data['reports'][i], ['assembly_info', 'sequencing_tech']),
                'Owner/Submitter': get_value(data['reports'][i], ['assembly_info', 'submitter']),

                'Current_accession': get_value(data['reports'][i], ['current_accession']),
                'Paired_assembly_accession': get_value(data['reports'][i], ['paired_accession']),
                'Assembly_release_date': get_value(data['reports'][i], ['assembly_info', 'release_date']),

                'Assembly_name': get_value(data['reports'][i], ['assembly_info', 'assembly_name']),
                'Assembly_level': get_value(data['reports'][i], ['assembly_info', 'assembly_level']),
                'Assembly_method': get_value(data['reports'][i], ['assembly_info', 'assembly_method']),
                'Assembly_status': get_value(data['reports'][i], ['assembly_info', 'assembly_status']),
                'Assembly_type': get_value(data['reports'][i], ['assembly_info', 'assembly_type']),
                
                'Total_chromosomes': get_value(data['reports'][i], ['assembly_stats', 'total_number_of_chromosomes']),
                'Total_length': get_value(data['reports'][i], ['assembly_stats', 'total_sequence_length']),
                'Contig_count': get_value(data['reports'][i], ['assembly_stats', 'number_of_contigs']),
                'Contig_L50': get_value(data['reports'][i], ['assembly_stats', 'contig_l50']),
                'Contig_N50': get_value(data['reports'][i], ['assembly_stats', 'contig_n50']),
                'Coverage_Depth': get_value(data['reports'][i], ['assembly_stats', 'genome_coverage']),
                'Scaffold_count': get_value(data['reports'][i], ['assembly_stats', 'number_of_scaffolds']),
                'Scaffold_L50': get_value(data['reports'][i], ['assembly_stats', 'scaffold_l50']),
                'Scaffold_N50': get_value(data['reports'][i], ['assembly_stats', 'scaffold_n50']),
                'Total_ungapped_length': get_value(data['reports'][i], ['assembly_stats', 'total_ungapped_length']),

                'Total_component_sequences': get_value(data['reports'][i], ['assembly_stats', 'number_of_component_sequences']),
                'GC_Percent': get_value(data['reports'][i], ['assembly_stats', 'gc_percent']),
                'GC_Count': get_value(data['reports'][i], ['assembly_stats', 'gc_count']),

                'Number_of_organelles': get_value(data['reports'][i], ['assembly_info', 'number_of_organelles']),
                'Organelle1' : get_value(data['reports'][i], ['organelle_info', 0, 'description']),
                'Organelle1_Sequence_Length': get_value(data['reports'][i], ['organelle_info', 0, 'total_seq_length']),
                'Organelle2' : get_value(data['reports'][i], ['organelle_info', 1, 'description']),
                'Organelle2_Sequence_Length': get_value(data['reports'][i], ['organelle_info', 1, 'total_seq_length']),

                # final attributes of assembly
                'Reference_genome?': 'Yes' if get_value(data['reports'][i], ['assembly_info', 'refseq_category']) == 'reference genome' else 'No',
                'Assmbly_atypical?': 'No' if get_value(data['reports'][i], ['assembly_info', 'atypical']) == '' else 'Yes',
                'Assembly_atypical_warning': 'NA' if get_value(data['reports'][i], ['assembly_info', 'atypical']) == '' else get_value(data['reports'][i], ['assembly_info', 'atypical', 'warnings', 0]),
                'Assembly_from_metagenome?': 'Yes' if 'metagenome' in get_value(data['reports'][i], ['assembly_info', 'genome_notes']) else 'No',
                'Linked_assembly': 'yes' if get_value(data['reports'][i], ['assembly_info', 'linked_assemblies', 0, 'linked_assembly']) else 'no',
                'Linked_assembly_accession_ID': get_value(data['reports'][i], ['assembly_info', 'linked_assemblies', 0, 'linked_assembly']),
                'Linked_assembly_type': get_value(data['reports'][i], ['assembly_info', 'linked_assemblies', 0, 'assembly_type']),
                'Source_Database': get_value(data['reports'][i], ['assembly_info', 'source_database']),
                'Bioproject_ID': get_value(data['reports'][i], ['assembly_info', 'bioproject_accession']),
                'Biosample_ID': get_value(data['reports'][i], ['assembly_info', 'biosample', 'accession']),
            }
            
            AssemblyData.append(Assembly_metadata)

            #print(AssemblyData)
        #print(AssemblyData)
        if len(AssemblyData) == total_count:
            logging.debug(f"Assembly attributes extracted successfully. len(AssemblyData) = {len(AssemblyData)}")

    return AssemblyData

def annotation_attributes_from_json(data:dict, total_count:int):
    """
    Extracts various annotation attributes from loaded Multi-JSON file with the help of `get_value` function.
    """
    if data:
        logging.info("Extracting Annotation attributes from JSON file.")

    AnnotationData = []

    if total_count > 0:
        for i in range(total_count):
            Annotation_metadata = {
                'Accession_ID': get_value(data['reports'][i], ['accession']),
                'Organism': get_value(data['reports'][i], ['organism', 'organism_name']),
                'strain': get_value(data['reports'][i], ['organism', 'infraspecific_names', 'strain']),
                'Taxonomy_ID': get_value(data['reports'][i], ['organism', 'tax_id']),
                'Sequencing_Technology': get_value(data['reports'][i], ['assembly_info', 'sequencing_tech']),
                'Owner/Submitter': get_value(data['reports'][i], ['assembly_info', 'submitter']),
                
                # based on an observation, the following keys are not present in all the genomes. So, we will try to get them from the 'annotation_info' section of the genome.    
                'Annotation_available?': 'No' if get_value(data['reports'][i], ['annotation_info']) == '' else 'Yes',
                
                'Annotation_from': (
                    'NCBI RefSeq' if (get_value(data['reports'][i], ['annotation_info']) != '' and 
                                    get_value(data['reports'][i], ['annotation_info', 'provider']) == 'NCBI RefSeq') or 
                                    (get_value(data['reports'][i], ['annotation_info']) != '' and 
                                    strip_gca_gcf_id(get_value(data['reports'][i], ['accession'])) == strip_gca_gcf_id(get_value(data['reports'][i], ['paired_accession'])))
                    else 'GenBank' if (get_value(data['reports'][i], ['annotation_info']) != '' and 
                                    get_value(data['reports'][i], ['annotation_info', 'provider']) == 'NCBI' and 
                                    strip_gca_gcf_id(get_value(data['reports'][i], ['accession'])) != strip_gca_gcf_id(get_value(data['reports'][i], ['paired_accession'])))
                    else get_value(data['reports'][i], ['annotation_info', 'provider'])
                ),

                'Annotation_method': get_value(data['reports'][i], ['annotation_info','method']),
                'Annotation_name': get_value(data['reports'][i], ['annotation_info', 'name']),
                'Annotation_pipeline': get_value(data['reports'][i], ['annotation_info', 'pipeline']),
                'Annotation_provider': get_value(data['reports'][i], ['annotation_info', 'provider']),
                'Annotation_version': get_value(data['reports'][i], ['annotation_info', 'release_version']),
                'Annotation_ReleaseDate': get_value(data['reports'][i], ['annotation_info', 'release_date']),
                'Annotation_status': get_value(data['reports'][i], ['annotation_info', 'status']),
                'Annotation_software_version': get_value(data['reports'][i], ['annotation_info', 'software_version']),                
                'Total_genes': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'total']),
                'Protein-coding_genes': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'protein_coding']),
                'Non-coding_genes': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'non_coding']),
                'Pseudogenes': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'pseudogene']),
                'other': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'other']),
        }

            if get_value(data['reports'][i], ['annotation_info', 'busco']) != '':
                busco_data = {
                    'Complete_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'complete']),
                    'Single-copy_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'single_copy']),
                    'Duplicated_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'duplicated']),
                    'Fragmented_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'fragmented']),
                    'Missing_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'missing']),
                    'Total_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'total_count']),
                    'BUSCO_lineage': get_value(data['reports'][i], ['annotation_info', 'busco', 'busco_lineage']),
                    'BUSCO_version': get_value(data['reports'][i], ['annotation_info', 'busco', 'busco_ver'])
                }
                Annotation_metadata.update(busco_data)

            if get_value(data['reports'][i], ['average_nucleotide_identity']) != '':
                ani_data = {
                    'ANI_best_match_score': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'ani']),
                    'ANI_best_matched_assembly': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'assembly']),
                    'ANI_best_matched_assembly\'s_coverage': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'assembly_coverage']),
                    'ANI_best_match_category': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'category']),
                    'ANI_best_matched_Organism': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'organism_name']),
                    'ANI_best_match_type_assembly\'s_coverage': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'type_assembly_coverage']),

                    'submitted_ANI_score': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'ani']),
                    'submitted_ANI\'s_matched_assembly': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'assembly']),
                    'submitted_ANI\'s_matched_assembly\'s_coverage': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'assembly_coverage']),
                    'submitted_ANI_category': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'category']),
                    'submitted_ANI\'s_matched_Organism': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'organism_name']),
                    'submitted_ANI\'s_type_assembly\'s_coverage': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'type_assembly_coverage']),

                    'ANI_match_status': get_value(data['reports'][i], ['average_nucleotide_identity', 'match_status']),
                    'ANI_comment': get_value(data['reports'][i], ['average_nucleotide_identity', 'comment']),
                    'Taxonomy_check_status':get_value(data['reports'][i], ['average_nucleotide_identity', 'taxonomy_check_status']),
                }

                Annotation_metadata.update(ani_data)
            
            if get_value(data['reports'][i], ['checkm_info']) != '':
                checkm_data = {
                    'CheckM_completeness': get_value(data['reports'][i], ['checkm_info', 'completeness']),
                    'CheckM_contamination': get_value(data['reports'][i], ['checkm_info', 'contamination']),
                    'CheckM_completeness_percentile': get_value(data['reports'][i], ['checkm_info', 'completeness_percentile']),
                    'CheckM_markerset': get_value(data['reports'][i], ['checkm_info', 'checkm_marker_set']),
                    'CheckM_markerset_rank/level': get_value(data['reports'][i], ['checkm_info', 'checkm_marker_set_rank']),
                    'CheckM_version': get_value(data['reports'][i], ['checkm_info', 'checkm_version']),
                }

                Annotation_metadata.update(checkm_data)

            AnnotationData.append(Annotation_metadata)
        
        if len(AnnotationData) == total_count:
            logging.debug(f"Annotation attributes extracted successfully. len(AnnotationData) = {len(AnnotationData)}")
    
    return AnnotationData

def get_allmetadata(data:dict, total_count:int, saving_file_path:str=None): 
    """
    Extracts all the metadata from the genome and will append in single file.
    """

    if saving_file_path is None:
        saving_file_path = os.getcwd()
        logging.debug(f"Since the `saving file path` is `None`, the current working directory is choosen. The WD is {saving_file_path}")
    else:
        logging.debug(f"the current saving path: {saving_file_path}")

    
    if data:
        logging.info("Extracting all metadata from JSON file.")

    allmetadata = []

    if total_count > 0:
        for i in range(total_count):
            all_metadata = {
                'Accession_ID': get_value(data['reports'][i], ['accession']),
                'Organism': get_value(data['reports'][i], ['organism', 'organism_name']),
                'strain': get_value(data['reports'][i], ['organism', 'infraspecific_names', 'strain']),
                'Common_name': get_value(data['reports'][i], ['organism', 'common_name']),
                'Taxonomy_ID': get_value(data['reports'][i], ['organism', 'tax_id']),
                'Sequencing_Technology' : get_value(data['reports'][i], ['assembly_info', 'sequencing_tech']),
                'Owner/Submitter': get_value(data['reports'][i], ['assembly_info', 'submitter']),
                'Bioproject_ID': get_value(data['reports'][i], ['assembly_info', 'bioproject_accession']),
                'Bioproject_parentID':get_value(data['reports'][i], ['assembly_info', 'bioproject_lineage', 0, 'bioprojects', 0, 'parent_accessions', 0]),
                'Bioproject_title': get_value(data['reports'][i], ['assembly_info', 'bioproject_lineage', 0, 'bioprojects', 0, 'title']),
                'Biosample_ID': get_value(data['reports'][i], ['assembly_info', 'biosample', 'accession']),
                'Biosample_title': get_value(data['reports'][i], ['assembly_info', 'biosample', 'description', 'title']),

                'Current_accession': get_value(data['reports'][i], ['current_accession']),
                'Paired_assembly_accession': get_value(data['reports'][i], ['paired_accession']),
                'Assembly_release_date': get_value(data['reports'][i], ['assembly_info', 'release_date']),

                'Assembly_name': get_value(data['reports'][i], ['assembly_info', 'assembly_name']),
                'Assembly_level': get_value(data['reports'][i], ['assembly_info', 'assembly_level']),
                'Assembly_method': get_value(data['reports'][i], ['assembly_info', 'assembly_method']),
                'Assembly_status': get_value(data['reports'][i], ['assembly_info', 'assembly_status']),
                'Assembly_type': get_value(data['reports'][i], ['assembly_info', 'assembly_type']),

                'Annotation_available?': 'No' if get_value(data['reports'][i], ['annotation_info']) == '' else 'Yes',
                
                'Annotation_from': (
                    'NCBI RefSeq' if (get_value(data['reports'][i], ['annotation_info']) != '' and 
                                    get_value(data['reports'][i], ['annotation_info', 'provider']) == 'NCBI RefSeq') or 
                                    (get_value(data['reports'][i], ['annotation_info']) != '' and 
                                    strip_gca_gcf_id(get_value(data['reports'][i], ['accession'])) == strip_gca_gcf_id(get_value(data['reports'][i], ['paired_accession'])))
                    else 'GenBank' if (get_value(data['reports'][i], ['annotation_info']) != '' and 
                                    get_value(data['reports'][i], ['annotation_info', 'provider']) == 'NCBI' and 
                                    strip_gca_gcf_id(get_value(data['reports'][i], ['accession'])) != strip_gca_gcf_id(get_value(data['reports'][i], ['paired_accession'])))
                    else get_value(data['reports'][i], ['annotation_info', 'provider'])
                ),

                'Annotation_method': get_value(data['reports'][i], ['annotation_info','method']),
                'Annotation_name': get_value(data['reports'][i], ['annotation_info', 'name']),
                'Annotation_pipeline': get_value(data['reports'][i], ['annotation_info', 'pipeline']),
                'Annotation_provider': get_value(data['reports'][i], ['annotation_info', 'provider']),
                'Annotation_version': get_value(data['reports'][i], ['annotation_info', 'release_version']),
                'Annotation_ReleaseDate': get_value(data['reports'][i], ['annotation_info', 'release_date']),
                'Annotation_status': get_value(data['reports'][i], ['annotation_info', 'status']),
                'Annotation_software_version': get_value(data['reports'][i], ['annotation_info', 'software_version']),                
                'Total_genes': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'total']),
                'Protein-coding_genes': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'protein_coding']),
                'Non-coding_genes': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'non_coding']),
                'Pseudogenes': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'pseudogene']),
                'other': get_value(data['reports'][i], ['annotation_info', 'stats', 'gene_counts', 'other']),
                
                'Total_chromosomes': get_value(data['reports'][i], ['assembly_stats', 'total_number_of_chromosomes']),
                'Total_length': get_value(data['reports'][i], ['assembly_stats', 'total_sequence_length']),
                'Contig_count': get_value(data['reports'][i], ['assembly_stats', 'number_of_contigs']),
                'Contig_L50': get_value(data['reports'][i], ['assembly_stats', 'contig_l50']),
                'Contig_N50': get_value(data['reports'][i], ['assembly_stats', 'contig_n50']),
                'Coverage_Depth': get_value(data['reports'][i], ['assembly_stats', 'genome_coverage']),
                'Scaffold_count': get_value(data['reports'][i], ['assembly_stats', 'number_of_scaffolds']),
                'Scaffold_L50': get_value(data['reports'][i], ['assembly_stats', 'scaffold_l50']),
                'Scaffold_N50': get_value(data['reports'][i], ['assembly_stats', 'scaffold_n50']),
                'Total_ungapped_length': get_value(data['reports'][i], ['assembly_stats', 'total_ungapped_length']),

                'Total_component_sequences': get_value(data['reports'][i], ['assembly_stats', 'number_of_component_sequences']),
                'GC_Percent': get_value(data['reports'][i], ['assembly_stats', 'gc_percent']),
                'GC_Count': get_value(data['reports'][i], ['assembly_stats', 'gc_count']),

                'Number_of_organelles': get_value(data['reports'][i], ['assembly_info', 'number_of_organelles']),
                'Organelle1' : get_value(data['reports'][i], ['organelle_info', 0, 'description']),
                'Organelle1_Sequence_Length': get_value(data['reports'][i], ['organelle_info', 0, 'total_seq_length']),
                'Organelle2' : get_value(data['reports'][i], ['organelle_info', 1, 'description']),
                'Organelle2_Sequence_Length': get_value(data['reports'][i], ['organelle_info', 1, 'total_seq_length']),

                # final attributes of assembly
                'Reference_genome?': 'Yes' if get_value(data['reports'][i], ['assembly_info', 'refseq_category']) == 'reference genome' else 'No',
                'Assmbly_atypical?': 'No' if get_value(data['reports'][i], ['assembly_info', 'atypical']) == '' else 'Yes',
                'Assembly_atypical_warning': 'NA' if get_value(data['reports'][i], ['assembly_info', 'atypical']) == '' else get_value(data['reports'][i], ['assembly_info', 'atypical', 'warnings', 0]),
                'Assembly_from_metagenome?': 'Yes' if 'metagenome' in get_value(data['reports'][i], ['assembly_info', 'genome_notes']) else 'No',
                'Linked_assembly': 'yes' if get_value(data['reports'][i], ['assembly_info', 'linked_assemblies', 0, 'linked_assembly']) else 'no',
                'Linked_assembly_accession_ID': get_value(data['reports'][i], ['assembly_info', 'linked_assemblies', 0, 'linked_assembly']),
                'Linked_assembly_type': get_value(data['reports'][i], ['assembly_info', 'linked_assemblies', 0, 'assembly_type']),
                'Source_Database': get_value(data['reports'][i], ['assembly_info', 'source_database']),
                'Bioproject_ID': get_value(data['reports'][i], ['assembly_info', 'bioproject_accession']),
                'Biosample_ID': get_value(data['reports'][i], ['assembly_info', 'biosample', 'accession']),

                'Model': get_value(data['reports'][i], ['assembly_info', 'biosample', 'models', 0]),
                'IFSAC_category': get_value(data['reports'][i], ['assembly_info', 'biosample', 'ifsac_category']),
                'Owner': get_value(data['reports'][i], ['assembly_info', 'biosample', 'owner', 'name']),
                'Owner_contact': get_value(data['reports'][i], ['assembly_info', 'biosample', 'owner', 'contacts', 0]),
                'Submission_date': get_value(data['reports'][i], ['assembly_info', 'biosample', 'submission_date']),
                'Last_update_date': get_value(data['reports'][i], ['assembly_info', 'biosample', 'last_update_date']),
                'Release_date': get_value(data['reports'][i], ['assembly_info', 'biosample', 'release_date']),
            }
            
            if get_value(data['reports'][i], ['annotation_info', 'busco']) != '':
                busco_data = {
                    'Complete_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'complete']),
                    'Single-copy_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'single_copy']),
                    'Duplicated_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'duplicated']),
                    'Fragmented_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'fragmented']),
                    'Missing_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'missing']),
                    'Total_buscos': get_value(data['reports'][i], ['annotation_info', 'busco', 'total_count']),
                    'BUSCO_lineage': get_value(data['reports'][i], ['annotation_info', 'busco', 'busco_lineage']),
                    'BUSCO_version': get_value(data['reports'][i], ['annotation_info', 'busco', 'busco_ver'])
                }
                all_metadata.update(busco_data)

            if get_value(data['reports'][i], ['average_nucleotide_identity']) != '':
                ani_data = {
                    'ANI_best_match_score': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'ani']),
                    'ANI_best_matched_assembly': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'assembly']),
                    'ANI_best_matched_assembly\'s_coverage': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'assembly_coverage']),
                    'ANI_best_match_category': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'category']),
                    'ANI_best_matched_Organism': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'organism_name']),
                    'ANI_best_match_type_assembly\'s_coverage': get_value(data['reports'][i], ['average_nucleotide_identity', 'best_ani_match', 'type_assembly_coverage']),

                    'submitted_ANI_score': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'ani']),
                    'submitted_ANI\'s_matched_assembly': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'assembly']),
                    'submitted_ANI\'s_matched_assembly\'s_coverage': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'assembly_coverage']),
                    'submitted_ANI_category': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'category']),
                    'submitted_ANI\'s_matched_Organism': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'organism_name']),
                    'submitted_ANI\'s_type_assembly\'s_coverage': get_value(data['reports'][i], ['average_nucleotide_identity', 'submitted_ani_match', 'type_assembly_coverage']),

                    'ANI_match_status': get_value(data['reports'][i], ['average_nucleotide_identity', 'match_status']),
                    'ANI_comment': get_value(data['reports'][i], ['average_nucleotide_identity', 'comment']),
                    'Taxonomy_check_status':get_value(data['reports'][i], ['average_nucleotide_identity', 'taxonomy_check_status']),
                }

                all_metadata.update(ani_data)
            
            if get_value(data['reports'][i], ['checkm_info']) != '':
                checkm_data = {
                    'CheckM_completeness': get_value(data['reports'][i], ['checkm_info', 'completeness']),
                    'CheckM_contamination': get_value(data['reports'][i], ['checkm_info', 'contamination']),
                    'CheckM_completeness_percentile': get_value(data['reports'][i], ['checkm_info', 'completeness_percentile']),
                    'CheckM_markerset': get_value(data['reports'][i], ['checkm_info', 'checkm_marker_set']),
                    'CheckM_markerset_rank/level': get_value(data['reports'][i], ['checkm_info', 'checkm_marker_set_rank']),
                    'CheckM_version': get_value(data['reports'][i], ['checkm_info', 'checkm_version']),
                }

                all_metadata.update(checkm_data)
            
            # since the attributes are really random and vary a lot from genome to genome, we will try to take all the possible attributes from all the assemblies   
            biosample_attributes = get_value(data['reports'][i], ['assembly_info', 'biosample', 'attributes'])
            # print(biosample_attributes)
            for attribute in biosample_attributes:
                # print(attribute)
                if 'name' in attribute and 'value' in attribute:
                    all_metadata.update({attribute['name']: attribute['value']})
                else:
                    logging.debug(f"Attribute missing 'name' or 'value' key: {attribute} in genome: {all_metadata['Accession_ID']}")


            allmetadata.append(all_metadata)

        if len(allmetadata) == total_count:
            logging.info(f"All metadata extracted successfully.")

    if allmetadata:
        all_metadata_df = pd.DataFrame(allmetadata)
        logging.debug(f"All metadata converted to dataframe.")
        logging.debug(f"shape of the all metadata dataframe: {all_metadata_df.shape}")
        
        os.makedirs("tmp", exist_ok=True)
        
        logging.info(f"saving the raw data for your future use.")

        date_time = datetime.now()
        date_time = str(date_time).replace(" ", "_")
        date_time = date_time.replace(":", "-")
        
        
        saving_file_name = "All_raw_metadata_" + date_time + ".tsv"
        saving_file_as = str(os.path.join(saving_file_path, saving_file_name)) 

        # df.to_csv(saving_file_as, index=False, sep='\t')

        all_metadata_df.to_csv(saving_file_as, sep='\t',index=False)
        logging.info(f"all raw metadata saved to {saving_file_as}.")

        logging.info(f"Removing the redundant data from the dataframe.")
        all_metadata_df = fix_gb_rf_annotation(all_metadata_df)
        logging.info(f"Redundant data removed from the dataframe.")
        logging.debug(f"Fixed the 'Annotation_from' column in the all metadata dataframe. Now the non-redundant data in dataframe is {all_metadata_df.shape}")
        logging.debug(f"Added `Annotation_category` column in the all metadata dataframe from `Annotation_from`.")

        logging.info(f"We have metadata for {all_metadata_df.shape[0]} genomes.")
        return all_metadata_df
    else:
        logging.info("Check the input file. It must be JSON from NCBI datasets.")

def normalizing_lon_lat(longitude_latitude_togethor:str):
    """
    This functions is to separate the longitude and latitude from the same column into two separate ones. The expected format for this is '40.4230 N 98.7372 W.' 
    """
    try:
        splitted_parts = longitude_latitude_togethor.split(' ')

        if len(splitted_parts) < 4:
            return (None, None)

        latitude_value = float(splitted_parts[0].strip())
        latitude_direction = splitted_parts[1].strip().upper()
        longitude_value = float(splitted_parts[2].strip())
        longitude_direction = splitted_parts[3].strip().upper()

        if latitude_direction == 'S':
            latitude_value = -latitude_value
        if longitude_direction == 'W':
            longitude_value = -longitude_value

        return (latitude_value, longitude_value)
    
    except:
        return (None, None)

from datetime import datetime
def extract_year(value):

    if pd.isna(value):
        return pd.NA

    value = str(value).strip()

    if value == "":
        return pd.NA

    # Missing values
    if value.lower() in {
        "na", "n/a", "nan", "none",
        "unknown", "missing", "not collected",
        "not available"
    }:
        return pd.NA

    current_year = datetime.now().year
    current_yy = current_year % 100

    # ----------------------------
    # Four-digit year
    # ----------------------------
    if re.fullmatch(r"\d{4}", value):
        year = int(value)
        return year if 1800 <= year <= current_year else pd.NA

    # ----------------------------
    # YYYY-MM or YYYY-MM-DD
    # YYYY/MM or YYYY/MM/DD
    # ----------------------------
    m = re.match(r"^(\d{4})[-/]", value)

    if m:
        year = int(m.group(1))
        return year if 1800 <= year <= current_year else pd.NA

    # ----------------------------
    # YYYYMMDD
    # ----------------------------
    if re.fullmatch(r"\d{8}", value):
        year = int(value[:4])
        return year if 1800 <= year <= current_year else pd.NA

    # ----------------------------
    # MM-YYYY / MM/YYYY
    # Example: 06-2001 -> 2001
    # ----------------------------
    m = re.fullmatch(r"\d{1,2}[-/](\d{4})", value)

    if m:
        year = int(m.group(1))
        return year if 1800 <= year <= current_year else pd.NA

    # ----------------------------
    # MM-DD-YYYY / MM/DD/YYYY
    # Example: 06-01-2001 -> 2001
    # ----------------------------
    m = re.fullmatch(r"\d{1,2}[-/]\d{1,2}[-/](\d{4})", value)

    if m:
        year = int(m.group(1))
        return year if 1800 <= year <= current_year else pd.NA

    # ----------------------------
    # MM-YY / MM/YY
    # Example: 06-79 -> 1979
    # ----------------------------
    m = re.fullmatch(r"\d{1,2}[-/](\d{2})", value)

    if m:
        yy = int(m.group(1))

        if yy <= current_yy:
            year = 2000 + yy
        else:
            year = 1900 + yy

        return year

    # ----------------------------
    # Two-digit year alone
    # Example: 79 -> 1979
    # ----------------------------
    if re.fullmatch(r"\d{2}", value):
        yy = int(value)

        if yy <= current_yy:
            year = 2000 + yy
        else:
            year = 1900 + yy

        return year

    # ----------------------------
    # General parser
    # ----------------------------
    try:
        dt = parser.parse(value, fuzzy=True)

        if 1800 <= dt.year <= current_year:
            return dt.year

    except Exception:
        pass

    # IMPORTANT:
    # Return pd.NA rather than None if nothing worked
    return pd.NA

def add_normalized_country_state_columns(df:pd.DataFrame, country_state_map:dict):
    """
    This function is to normalize the `geo_loc` data into country names, state names and other related parameters from local database - country_state_map.
    """

    logging.debug(f"Adding normalized country-state columns to the dataframe.")
    logging.debug(f"Length of the dataframe: {len(df)}")
    logging.debug(f"Length of the country_state_map: {len(country_state_map)}")
    try:
        df['country_name'] = df['geo_loc_name'].map(
            lambda x: country_state_map.get(x, (None, None, None, None, None))[0] if x != None else None)
        df['country_common_name'] = df['geo_loc_name'].map(
            lambda x: country_state_map.get(x, (None, None, None, None, None))[1] if x != None else None)
        df['country_three_lettered_name'] = df['geo_loc_name'].map(
            lambda x: country_state_map.get(x, (None, None, None, None, None))[2] if x != None else None)
        df['state_name'] = df['geo_loc_name'].map(
            lambda x: country_state_map.get(x, (None, None, None, None, None))[3] if x != None else None)
        df['state_code'] = df['geo_loc_name'].map(
        lambda x: country_state_map.get(x, (None, None, None, None, None))[4] if x != None else None)
        logging.info(f"Nomalized country-state data columns are added to the dataframe.")
    except Exception as e:
        df['country_name'] = None
        df['country_common_name'] = None
        df['country_three_lettered_name'] = None
        df['state_name'] = None
        df['state_code'] = None
        logging.error(f"Error in adding the normalized country-state columns to the dataframe: {e}")  

def categorizing_annotation_for_graph(value:str):
    """
    This function is to just categorize annotation data for better filtering for plot.
    """
    if value == "GenBank":
        return "GenBank"
    elif value == "NCBI RefSeq":
        return "NCBI RefSeq"
    elif value is None or pd.isna(value):
        return "No Annotation"
    else:
        return "Others"

def extract_year_from_dandt_stamp(data, column_name:str, new_column_name:str):

    data[column_name] = pd.to_datetime(data[column_name], errors='coerce')

    data[new_column_name] = data[column_name].dt.year

    return data

def get_country_name(fuzzy_country_names:str):
    """
    This function is to normalize the `geo_loc` data into country names and if possible, states names. 
    """
    logging.debug(f"Normalizing the country names from the column: {fuzzy_country_names}")

    if isinstance(fuzzy_country_names, str):

        if not fuzzy_country_names:
            return (None, None, None, None, None)

        fuzzy_country_names = fuzzy_country_names.strip().title()
        
        if fuzzy_country_names.lower() in ["na", "n/a", "not collected", "not applicable", "missing", "unknown", "not provided", "none", " ", "nan", "not determined", "not collected", "null"]: # if the country name is not provided or is missing
            logging.debug(f"Country name is not provided or is missing.")
            return (None, None, None, None, None)
        
        splitted_parts = [part.strip().replace('.', '').replace(')', '').replace('(', '').title() for part in fuzzy_country_names.split(':')]
        logging.debug(f"Splitting the country names into parts: {splitted_parts}")
        
        if len(splitted_parts) == 1: # if the country name is not separated by a colon
            logging.debug(f"Length of the splitted parts is 1.")
            try: # opening the try block for looking for the country name in the splitted_parts[0]
                country_fuzzy = pc.countries.search_fuzzy(splitted_parts[0])
                country_name = country_fuzzy[0].name
                logging.debug(f"Country name is: {country_name}")
                country_three_lettered_name = country_fuzzy[0].alpha_3
                logging.debug(f"Country three-lettered name is: {country_three_lettered_name}")
                try:
                    country_common_name = country_fuzzy[0].common_name
                    logging.debug(f"Country common name is: {country_common_name}")
                except:
                    country_common_name = country_name
                    logging.debug(f"Country common name is same as country name.")
                try: # opening the try block for looking for the state name in order to check weather the identified country name is derived from the state's name or country's name
                    state_fuzzy = pc.subdivisions.search_fuzzy(splitted_parts[0])
                    state_name = state_fuzzy[0].name
                    logging.debug(f"State name is: {state_name}")
                    if state_fuzzy[0].country_code == country_fuzzy[0].alpha_2: # if the state name is in the same country
                        
                        if len(state_name) > len(splitted_parts[0]):
                            logging.debug(f"State belongs to the country. However, the query is a country name and due to fuzzy search is giving some random state search. returning just the country name.")
                            return (country_name, country_common_name, country_three_lettered_name, None, None)
                        else:
                            logging.debug(f"State name is in the same country. Returning the country and state names.")
                            state_code = state_fuzzy[0].code
                            return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                    else: # if the state name is not in the same country return just the country name as the splitted_parts[0] is a country's name only
                        logging.debug(f"State doesn't belong to the country. Due to the fuzzy search, it has selected some state. Returning just the country name.")
                        return (country_name, country_common_name, country_three_lettered_name, None, None)
                except: # closing the try block for looking for the state name in the splitted_parts[0]
                    logging.debug(f"State name couldn't be found. Returning just the country name.")
                    return (country_name, country_common_name, country_three_lettered_name, None, None)   
            except: # closing the try block for looking for the country name in the splitted_parts[0]
                try: # trying to get the country name from the geolocator
                    logging.debug(f"Trying to get the country name from the geolocator for entire fuzzy_country_name.")
                    metageolocator = Nominatim(user_agent="metaminerGUI")

                    
                    address = metageolocator.geocode(fuzzy_country_names, addressdetails=True, language='en', exactly_one=True)
                    time.sleep(1)
                    logging.debug(f"Address is: {address}")

                    if address:
                        formatted_address = address.raw.get('address')
                        ad_country_name = formatted_address.get('country')
                        
                        country_name = pc.countries.search_fuzzy(ad_country_name)[0].name
                        logging.debug(f"Country name is: {country_name}")

                        country_three_lettered_name = pc.countries.search_fuzzy(country_name)[0].alpha_3
                        logging.debug(f"Country three-lettered name is: {country_three_lettered_name}")
                        try:
                            country_common_name = pc.countries.search_fuzzy(country_name)[0].common_name
                            logging.debug(f"Country common name is: {country_common_name}")
                        except:
                            country_common_name = country_name
                            logging.debug(f"Country common name is same as country name.")
                        
                        ad_state_name = formatted_address.get('state')
                        state_name = pc.subdivisions.search_fuzzy(ad_state_name)[0].name
                        

                        if not state_name:
                            try:
                                state_name = formatted_address.get('county')
                                logging.debug(f"State name is: {state_name}")
                                try:
                                    state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                    logging.debug(f"State code is: {state_code}")
                                except:
                                    state_code = None
                                    logging.debug(f"State code couldn't be found.")
                            except:
                                state_name = formatted_adress.get('island')
                                logging.debug(f"State name is: {state_name}")
                                try:
                                    state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                    logging.debug(f"State code is: {state_code}")
                                except:
                                    state_code = None
                                    logging.debug(f"State code couldn't be found.")
                        
                        logging.debug(f"Country name is: {country_name} and state name is: {state_name}")
                        state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                    else:
                        return (None, None, None, None, None)

                except: # if the country name couldn't be found using the geolocator
                    logging.debug(f"Country couldn't be found using the geolocator. returning 'Unclassified' for both country and state.")
                    return (None, None, None, None, None)

        elif len(splitted_parts) == 2: # if the country name and state name is separated by a colon
            logging.debug(f"Length of the splitted parts is 2.")
            try: # opening the try block for looking for the country name in fuzzy_country_names
                try: # opening the try block for looking for the country name in the splitted_parts[0]
                    country_fuzzy = pc.countries.search_fuzzy(splitted_parts[0])
                    country_name = country_fuzzy[0].name
                    logging.debug(f"Country name is: {country_name}")
                    country_three_lettered_name = country_fuzzy[0].alpha_3
                    logging.debug(f"Country three-lettered name is: {country_three_lettered_name}")
                    try:
                        country_common_name = country_fuzzy[0].common_name
                        logging.debug(f"Country common name is: {country_common_name}")
                    except:
                        country_common_name = country_name
                        logging.debug(f"Country common name is same as country name.")
                
                    if country_name: # if country name is identified in the splitted_parts[0]
                        
                        try:
                            logging.debug("testing if the splitted_parts[0] from which the country names is extracted is a country name or state's name.")
                            state_fuzzy = pc.subdivisions.search_fuzzy(splitted_parts[0])
                            state_name = state_fuzzy[0].name

                            if state_name: # weather the idenfied country name is derived from the state's name or country's name
                                
                                if state_fuzzy[0].country_code == country_fuzzy[0].alpha_2: # if the state name is in the same country
                                    
                                    if len(state_name) > len(splitted_parts[0]):
                                        logging.debug(f"State belongs to the country. However, the query is a country name and due to fuzzy search is giving some random state search.")
                                        
                                    else:
                                        state_code = state_fuzzy[0].code
                                        logging.debug(f"State code is: {state_code}")
                                        logging.debug(f"State name is in the same country. Returning the country and state names.")
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                            
                                else: # if the state name is not in the same country return just the country name as the splitted_parts[0] is a country's name only
                                    logging.debug(f"State doesn't belong to the country. Due to the fuzzy search, it has selected some state.")
                                    logging.debug(f"splitted_parts[0] is a country's name only!")
                                
                            else:
                                logging.debug("splitted_parts[0] is a country's name only!")
                        
                        except: # closing the try block for looking for the state name in the splitted_parts[0]
                            logging.debug("splitted_parts[0] is a country's name only!")

                        splitting_parts_into_more = [part.strip() for part in splitted_parts[1].split(',')]
                        logging.debug(f"Splitting the {splitted_parts[1]} with ',': {splitting_parts_into_more}")

                        try: # opening the try block for looking for the state name in the splitted_parts[1]
                        
                            if len(splitting_parts_into_more) == 1: # if some codes are there in the splitted_parts[1] for state names
                                logging.debug(f"Length of the splitting_parts_into_more is 1.")
                            
                                if len(splitting_parts_into_more[0]) > 3: # if it's not code and some state name
                                    logging.debug(f"Length of the splitting_parts_into_more[0] is greater than 3.")
                                    state_fuzzy = pc.subdivisions.search_fuzzy(splitting_parts_into_more[0])
                                    state_name = state_fuzzy[0].name
                                    logging.debug(f"State name is: {state_name}")
                                    if state_fuzzy[0].country_code == country_fuzzy[0].alpha_2:
                                        logging.debug(f"State name is in the same country. Returning the country and state names.")
                                        returned = True
                                        state_code = state_fuzzy[0].code
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                    
                                    if not returned:
                                        logging.debug("No state name was able to identified using the pycountry. Finding the state name using the geolocator.")
                                        try: # trying to get the state name from the geolocator from splitted_parts[1]
                                            logging.debug(f"Trying to get the state name from the geolocator for {splitted_parts[1]}.") 
                                            metageolocator = Nominatim(user_agent="metaminerGUI")
                                    
                                            address = metageolocator.geocode(fuzzy_country_names, addressdetails=True, language='en', exactly_one=True)
                                            logging.debug(f"Address is: {address.raw.get('address')}")

                                            time.sleep(1)

                                            if address:
                                                formatted_address = address.raw.get('address')
                                                ad_state_name = formatted_address.get('state')
                                                state_name = pc.subdivisions.search_fuzzy(ad_state_name)[0].name
                                                if not state_name:
                                                    logging.debug("State name couldn't be found. Trying to get the county name.")
                                                    try:
                                                        state_name = formatted_address.get('county')
                                                        logging.debug(f"State name is: {state_name}")
                                                        try:
                                                            state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                            logging.debug(f"State code is: {state_code}")
                                                        except:
                                                            state_code = None
                                                            logging.debug(f"State code couldn't be found.")
                                                    except:
                                                        state_name = formatted_address.get('island')
                                                        logging.debug(f"State name is: {state_name}")
                                                        try:
                                                            state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                            logging.debug(f"State code is: {state_code}")
                                                        except:
                                                            state_code = None
                                                            logging.debug(f"State code couldn't be found.")

                                                if country_fuzzy[0].alpha_2 == formatted_address.get('country_code').upper():
                                                    logging.debug(f"State name is: {state_name}")
                                                    returned = True
                                                    state_code = pc.subdivisons.search_fuzzy(state_name)[0].code
                                                    return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                                else:
                                                    logging.debug(f"State doesn't belong to the country. The splitted_parts[1] may contain something other than the state name.")
                                                    returned = True
                                                    return (country_name, country_common_name, country_three_lettered_name, None, None)
                                            else:
                                                return (country_name, country_common_name, country_three_lettered_name, None, None)

                                        except: # if the state name couldn't be found using the geolocator
                                            logging.debug(f"State couldn't be found using the geolocator. Returning just the country name.")
                                            return (country_name, country_common_name, country_three_lettered_name, None, None)
                                        
                                elif len(splitting_parts_into_more[0]) <= 2: # if it's a two-letter code for the state
                                    logging.debug(f"Length of the splitting_parts_into_more[0] is less than or equal to 2.")
                                    two_letter_country_code = country_fuzzy[0].alpha_2
                                    state_code = two_letter_country_code + "-" + splitting_parts_into_more[0].strip().upper()
                                    logging.debug(f"State code is: {state_code}")
                                    state_name = pc.subdivisions.get(code = state_code).name
                                    if state_name:
                                        logging.debug(f"State name is: {state_name}")
                                        returned = True
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                elif len(splitting_parts_into_more[0]) == 3: # if it's a three-letter code for the state
                                    logging.debug(f"Length of the splitting_parts_into_more[0] is 3.")
                                    two_letter_country_code = country_fuzzy[0].alpha_2
                                    state_code = two_letter_country_code + "-" + splitting_parts_into_more[0:2].strip().upper()
                                    logging.debug(f"State code is: {state_code}")
                                    state_name = pc.subdivisions.get(code = state_code).name
                                    if state_name:
                                        logging.debug(f"State name is: {state_name}")
                                        returned = True
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                            else: # if states name is in more than few parts
                                logging.debug(f"Length of the splitting_parts_into_more is greater than 1.")
                                logging.debug(f"Looping through the {splitting_parts_into_more}.")
                                for i in range(len(splitting_parts_into_more)):
                                    logging.debug(f"part: {splitting_parts_into_more[i]}")
                                    try:
                                        if len(splitting_parts_into_more[i]) <= 2: # if it's a two-letter code for the state
                                            logging.debug(f"Length of the splitting_parts_into_more[i] is less than or equal to 2.")
                                            two_letter_country_code = country_fuzzy[0].alpha_2
                                            state_code = two_letter_country_code + "-" + splitting_parts_into_more[i].strip().upper()
                                            logging.debug(f"State code is: {state_code}")
                                            state_name = pc.subdivisions.get(code = state_code).name
                                            if state_name:
                                                logging.debug(f"State name is: {state_name}")
                                                returned = True
                                                return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                        elif len(splitting_parts_into_more[i]) == 3: # if it's a three-letter code for the state
                                            logging.debug(f"Length of the splitting_parts_into_more[i] is 3.")
                                            two_letter_country_code = country_fuzzy[0].alpha_2
                                            state_code = two_letter_country_code + "-" + splitting_parts_into_more[i][0:2].strip().upper()
                                            logging.debug(f"State code is: {state_code}")
                                            state_name = pc.subdivisions.get(code = state_code).name
                                            if state_name:
                                                logging.debug(f"State name is: {state_name}")
                                                returned = True
                                                return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                        else: 
                                            logging.debug(f"Length of the {splitting_parts_into_more[i]} is greater than 3.")                                       
                                            state_fuzzy = pc.subdivisions.search_fuzzy(splitting_parts_into_more[i])
                                            state_name = state_fuzzy[0].name
                                            logging.debug(f"State name is: {state_name}")
                                            if state_fuzzy[0].country_code == country_fuzzy[0].alpha_2:
                                                logging.debug(f"State name is in the same country. Returning the country and state names.")
                                                returned = True
                                                state_code = state_fuzzy[0].code
                                                logging.debug(f"State code is: {state_code}")
                                                return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                            else:
                                                logging.debug(f"Identified state doesn't belong to the same country.")
                                    except:
                                        logging.debug(f"State name couldn't be found. Looping through the next part.")
                                        continue

                                if not returned:
                                    logging.debug("No state name was able to identified using the pycountry and looping in splitted_parts[1]")
                                    try: # trying to get the state name from the geolocator from splitted_parts[1]
                                        logging.debug(f"Trying to get the state name from the geolocator for splitted_parts[1].") 
                                        metageolocator = Nominatim(user_agent="metaminerGUI")
                                
                                        address = metageolocator.geocode(fuzzy_country_names, addressdetails=True, language='en', exactly_one=True)
                                        logging.debug(f"Address is: {address.raw.get('address')}")

                                        time.sleep(1)

                                        if address:
                                            formatted_address = address.raw.get('address')
                                            ad_state_name = formatted_address.get('state')
                                            state_name = pc.subdivisions.search_fuzzy(ad_state_name)[0].name
                    
                                            if not state_name:
                                                logging.debug("State name couldn't be found. Trying to get the county name.")
                                                try:
                                                    state_name = formatted_address.get('county')
                                                    logging.debug(f"State name is: {state_name}")
                                                    try:
                                                        state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                        logging.debug(f"State code is: {state_code}")
                                                    except:
                                                        state_code = None
                                                        logging.debug(f"State code couldn't be found.")
                                                except:
                                                    state_name = formatted_address.get('island')
                                                    logging.debug(f"State name is: {state_name}")
                                                    try:
                                                        state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                        logging.debug(f"State code is: {state_code}")
                                                    except:
                                                        state_code = None
                                                        logging.debug(f"State code couldn't be found.") 

                                            if country_fuzzy[0].alpha_2 == formatted_address.get('country_code').upper():
                                                logging.debug(f"State name is: {state_name}")
                                                returned = True
                                                state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                logging.debug(f"State code is: {state_code}")
                                                return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                            else:
                                                logging.debug(f"State doesn't belong to the country. The splitted_parts[1] may contain something other than the state name.")
                                                returned = True
                                                return (country_name, country_common_name, country_three_lettered_name, None, None)
                                        else:
                                            return (country_name, country_common_name, country_three_lettered_name, None, None)

                                    except: # if the state name couldn't be found using the geolocator
                                        logging.debug(f"State couldn't be found using the geolocator. Returning just the country name.")
                                        return (country_name, country_common_name, country_three_lettered_name, None, None)

                        except: # closing the try block for looking for the state name in the splitted_parts[1] with geolocator
                        
                            try: # trying to get the state name from the geolocator from splitted_parts[1]
                                logging.debug(f"Trying to get the state name from the geolocator for splitted_parts[1].") 
                                metageolocator = Nominatim(user_agent="metaminerGUI")
                                
                                address = metageolocator.geocode(fuzzy_country_names, addressdetails=True, language='en', exactly_one=True)
                                logging.debug(f"Address is: {address.raw.get('address')}")

                                time.sleep(1)

                                if address:
                                    formatted_address = address.raw.get('address')
                                    ad_state_name = formatted_address.get('state')
                                    state_name = pc.subdivisions.search_fuzzy(ad_state_name)[0].name

                                    if not state_name:
                                        logging.debug("State name couldn't be found. Trying to get the county name.")
                                        try:
                                            state_name = formatted_address.get('county')
                                            logging.debug(f"State name is: {state_name}")
                                            try:
                                                state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                logging.debug(f"State code is: {state_code}")
                                            except:
                                                state_code = None
                                                logging.debug(f"State code couldn't be found.")
                                        except:
                                            state_name = formatted_address.get('island')
                                            logging.debug(f"State name is: {state_name}")
                                            try:
                                                state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                logging.debug(f"State code is: {state_code}")
                                            except:
                                                state_code = None
                                                logging.debug(f"State code couldn't be found.")

                                    if country_fuzzy[0].alpha_2 == formatted_address.get('country_code').upper():
                                        logging.debug(f"State name is: {state_name}")
                                        state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                        logging.debug(f"State code is: {state_code}")
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                    else:
                                        logging.debug(f"State doesn't belong to the country. The splitted_parts[1] may contain something other than the state name.")
                                        return (country_name, country_common_name, country_three_lettered_name, None, None)
                                else:
                                    return (country_name, country_common_name, country_three_lettered_name, None, None)

                            except: # if the state name couldn't be found using the geolocator
                                logging.debug(f"State couldn't be found using the geolocator. Returning just the country name.")
                                return (country_name, country_common_name, country_three_lettered_name, None, None)
                
                except: # closing the try block via looking for the country name in the splitted_parts[1]
                    logging.debug("Trying to find country name in splitted_parts[1].")
                    country_fuzzy = pc.countries.search_fuzzy(splitted_parts[1])
                    country_name = country_fuzzy[0].name
                    logging.debug(f"Country name is: {country_name}")
                    country_three_lettered_name = country_fuzzy[0].alpha_3
                    logging.debug(f"Country three-lettered name is: {country_three_lettered_name}")
                    try:
                        country_common_name = country_fuzzy[0].common_name
                        logging.debug(f"Country common name is: {country_common_name}")
                    except:
                        country_common_name = country_name
                        logging.debug(f"Country common name is same as country name.")

                    if country_name: # if country name is identified in the splitted_parts[1]
                        logging.debug(f"country name is found in {splitted_parts[1]}.")
                        
                        splitting_parts_into_more = [part.strip() for part in splitted_parts[1].split(',')]
                        logging.debug(f"Splitting the first part into more parts: {splitting_parts_into_more}")
                        
                        try: # opening the try block for looking for the state name in the splitted_parts[0]
                            
                            if len(splitting_parts_into_more) == 1: # if some codes are there in the splitted_parts[0] for state names
                                logging.debug(f"Length of the splitting_parts_into_more is 1.")
                                
                                if len(splitting_parts_into_more[0]) > 3: # if it's not code and some state name
                                    logging.debug(f"Length of the splitting_parts_into_more[0] is greater than 3.")
                                    state_fuzzy = pc.subdivisions.search_fuzzy(splitting_parts_into_more[0])
                                    state_name = state_fuzzy[0].name
                                    logging.debug(f"State name is: {state_name}")
                                    if state_fuzzy[0].country_code == country_fuzzy[0].alpha_2:
                                        logging.debug(f"State name is in the same country. Returning the country and state names.")
                                        returned = True
                                        state_code = state_fuzzy[0].code
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                elif len(splitting_parts_into_more[0]) <= 2: # if it's a two-letter code for the state
                                    logging.debug(f"Length of the splitting_parts_into_more[0] is less than or equal to 2.")
                                    two_letter_country_code = country_fuzzy[0].alpha_2
                                    state_code = two_letter_country_code + "-" + splitting_parts_into_more[0].strip().upper()
                                    logging.debug(f"State code is: {state_code}")
                                    state_name = pc.subdivisions.get(code = state_code).name
                                    if state_name:
                                        logging.debug(f"State name is: {state_name}")
                                        returned = True
                                        state_code = state_fuzzy[0].code
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                elif len(splitting_parts_into_more[0]) == 3: # if it's a three-letter code for the state
                                    logging.debug(f"Length of the splitting_parts_into_more[0] is 3.")
                                    two_letter_country_code = country_fuzzy[0].alpha_2
                                    state_code = two_letter_country_code + "-" + splitting_parts_into_more[0:2].strip().upper()
                                    logging.debug(f"State code is: {state_code}")
                                    state_name = pc.subdivisions.get(code = state_code).name
                                    if state_name:
                                        logging.debug(f"State name is: {state_name}")
                                        returned = True
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                        
                            
                            else: # if states name is in more than few parts
                                logging.debug(f"Length of the splitting_parts_into_more[0] is greater than 1.")
                                logging.debug(f"Looping through the {splitting_parts_into_more}.")
                                for i in range(len(splitting_parts_into_more)): # looping through the splitted_parts[0] to find the state name
                                    logging.debug(f"part: {splitting_parts_into_more[i]}")
                                    try:
                                        if len(splitting_parts_into_more[i]) <= 2: # if it's a two-letter code for the state
                                            logging.debug(f"Length of the splitting_parts_into_more[i] is less than or equal to 2.")
                                            two_letter_country_code = country_fuzzy[0].alpha_2
                                            state_code = two_letter_country_code + "-" + splitting_parts_into_more[i].strip().upper()
                                            logging.debug(f"State code is: {state_code}")
                                            state_name = pc.subdivisions.get(code = state_code).name
                                            if state_name:
                                                logging.debug(f"State name is: {state_name}")
                                                return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                        elif len(splitting_parts_into_more[i]) == 3: # if it's a three-letter code for the state
                                            logging.debug(f"Length of the splitting_parts_into_more[i] is 3.")
                                            two_letter_country_code = country_fuzzy[0].alpha_2
                                            state_code = two_letter_country_code + "-" + splitting_parts_into_more[i][0:2].strip().upper()
                                            logging.debug(f"State code is: {state_code}")
                                            state_name = pc.subdivisions.get(code = state_code).name
                                            if state_name:
                                                logging.debug(f"State name is: {state_name}")
                                                return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                        else: 
                                            logging.debug(f"Length of the splitting_parts_into_more[i] is greater than 3.")                                       
                                            state_fuzzy = pc.subdivisions.search_fuzzy(splitting_parts_into_more[i])
                                            state_name = state_fuzzy[0].name
                                            logging.debug(f"State name is: {state_name}")
                                            if state_fuzzy[0].country_code == country_fuzzy[0].alpha_2:
                                                logging.debug(f"State name is in the same country. Returning the country and state names.")
                                                state_code = state_fuzzy[0].code
                                                return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                            else:
                                                logging.debug(f"Identified state doesn't belong to the same country.")
                                    except:
                                        logging.debug(f"State name couldn't be found. Looping through the next part.")
                                        continue
                                
                                if not state_name:
                                    logging.debug("No state name was able to identified using the pycountry and looping in splitted_parts[0]")
                                    try: # trying to get the state name from the geolocator from splitted_parts[1]
                                        logging.debug(f"Trying to get the state name from the geolocator for splitted_parts[0].") 
                                        metageolocator = Nominatim(user_agent="metaminerGUI")
                                
                                        address = metageolocator.geocode(fuzzy_country_names, addressdetails=True, language='en', exactly_one=True)
                                        logging.debug(f"Address is: {address.raw.get('address')}")

                                        time.sleep(1)

                                        if address:
                                            formatted_address = address.raw.get('address')
                                            ad_state_name = formatted_address.get('state')
                                            state_name = pc.subdivisions.search_fuzzy(ad_state_name)[0].name
                                            if not state_name:
                                                logging.debug("State name couldn't be found. Trying to get the county name.")
                                                try:
                                                    state_name = formatted_address.get('county')
                                                    logging.debug(f"State name is: {state_name}")
                                                    try:
                                                        state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                        logging.debug(f"State code is: {state_code}")
                                                    except:
                                                        state_code = None
                                                        logging.debug(f"State code couldn't be found.")
                                                except:
                                                    state_name = formatted_address.get('island')
                                                    logging.debug(f"State name is: {state_name}")
                                                    try:
                                                        state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                        logging.debug(f"State code is: {state_code}")
                                                    except:
                                                        state_code = None
                                                        logging.debug(f"State code couldn't be found.")
                                            
                                            if country_fuzzy[0].alpha_2 == formatted_address.get('country_code').upper():
                                                logging.debug(f"State name is: {state_name}")
                                                state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                logging.debug(f"State code is: {state_code}")
                                                return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                                
                                            else:
                                                logging.debug(f"State doesn't belong to the country. The splitted_parts[1] may contain something other than the state name.")
                                                return (country_name, country_common_name, country_three_lettered_name, None, None)
                                        else:
                                            return (country_name, country_common_name, country_three_lettered_name, None, None)

                                    except: # if the state name couldn't be found using the geolocator
                                        logging.debug(f"State couldn't be found using the geolocator. Returning just the country name.")
                                        return (country_name, country_common_name, country_three_lettered_name, None, None)

                        except: # closing the try block for looking for the state name in the splitted_parts[0] with geolocator

                            try: # trying to get the state name from the geolocator from splitted_parts[0]
                                logging.debug(f"Trying to get the country name from the geolocator for splitted_parts[0].") 
                                metageolocator = Nominatim(user_agent="metaminerGUI")
                                
                                address = metageolocator.geocode(fuzzy_country_names, addressdetails=True, language='en', exactly_one=True)
                                logging.debug(f"Address is: {address.raw.get('address')}")

                                time.sleep(1)
                                
                                if address:
                                    formatted_address = address.raw.get('address')
                                    ad_state_name = formatted_address.get('state')
                                    state_name = pc.subdivisions.search_fuzzy(ad_state_name)[0].name
                                    if not state_name:
                                        logging.debug("State name couldn't be found. Trying to get the county name.")
                                        try:
                                            state_name = formatted_address.get('county')
                                            logging.debug(f"State name is: {state_name}")
                                            try:
                                                state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                logging.debug(f"State code is: {state_code}")
                                            except:
                                                state_code = None
                                                logging.debug(f"State code couldn't be found.")
                                        except:
                                            state_name = formatted_address.get('island')
                                            logging.debug(f"State name is: {state_name}")
                                            try:
                                                state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                                logging.debug(f"State code is: {state_code}")
                                            except:
                                                state_code = None
                                                logging.debug(f"State code couldn't be found.")
                                    if country_fuzzy[0].alpha_2 == formatted_address.get('country_code').upper():
                                        logging.debug(f"State name is: {state_name}")
                                        state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                        logging.debug(f"State code is: {state_code}")
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                    else:
                                        return (country_name, country_common_name, country_three_lettered_name, None, None)
                                else:
                                    return (country_name, country_common_name, country_three_lettered_name, None, None)

                            except: # in case the state name couldn't be found using the geolocator
                                logging.debug(f"State couldn't be found using the geolocator. Returning just the country name.")
                                return (country_name, country_common_name, country_three_lettered_name, None, None)

            except: # closing the try block as the country name couldn't be found in the splitted_parts[0] and splitted_parts[1] and trying to get the country name from the geolocator

                try: # trying to get the country name from the geolocator
                    metageolocator = Nominatim(user_agent="metaminerGUI")
                    logging.debug(f"Trying to get the country name from the geolocator for entire fuzzy_country_name. Primary try block failed when len(splitted_parts) ==2.")
                    address = metageolocator.geocode(fuzzy_country_names, addressdetails=True, language='en', exactly_one=True)
                    logging.debug(f"Address is: {address.raw.get('address')}")
                    time.sleep(1)
                    if address:
                        formatted_address = address.raw.get('address')
                        ad_country_name = formatted_address.get('country')
                        
                        country_name = pc.countries.search_fuzzy(ad_country_name)[0].name
                        logging.debug(f"Country name is: {country_name}")
                        
                        country_three_lettered_name = pc.countries.search_fuzzy(country_name)[0].alpha_3
                        logging.debug(f"Country three-lettered name is: {country_three_lettered_name}")

                        try:
                            country_common_name = pc.countries.search_fuzzy(country_name)[0].common_name
                            logging.debug(f"Country common name is: {country_common_name}")
                        except:
                            country_common_name = country_name
                            logging.debug(f"Country common name is same as country name.")

                        ad_state_name = formatted_address.get('state')
                        state_name = pc.subdivisions.search_fuzzy(ad_state_name)[0].name

                        if not state_name:
                            logging.debug("State name couldn't be found. Trying to get the county name.")
                            try:
                                state_name = formatted_address.get('county')
                                logging.debug(f"State name is: {state_name}")
                                try:
                                    state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                    logging.debug(f"State code is: {state_code}")
                                except:
                                    state_code = None
                                    logging.debug(f"State code couldn't be found.")
                            except:
                                state_name = formatted_address.get('island')
                                logging.debug(f"State name is: {state_name}")
                                try:
                                    state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                    logging.debug(f"State code is: {state_code}")
                                except:
                                    state_code = None
                                    logging.debug(f"State code couldn't be found.")

                        state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                    else:
                        return (None, None, None, None, None)

                except: # if the country name couldn't be found using the geolocator
                    logging.debug(f"Country couldn't be found using the geolocator. returning 'Unclassified' for both country and state.")
                    return (None, None, None, None, None)
        
        elif len(splitted_parts) > 2: # if the country name and state name is separated by a colon and there are more than 2 parts
            logging.debug(f"Length of the splitted parts is greater than 2.")
            splitting_more = [part.split(',') for part in splitted_parts]
            all_flattened = [part.strip() for sublist in splitting_more for part in sublist]
            logging.debug(f"Flattening the splitted parts: {all_flattened}")
            try: # opening the try block for looking for the country name in the all_flattened
        
                country_name = None # initializing the country name

                for flattened_part in all_flattened: # looping through the all_flattened to find the country name
                    try:
                        if country_name:
                            country_three_lettered_name = country_fuzzy[0].alpha_3
                            logging.debug(f"Country three-lettered name is: {country_three_lettered_name}")
                            try:
                                country_common_name = country_fuzzy[0].common_name
                                logging.debug(f"the common name for country: {country_common_name}")
                            except:
                                country_common_name = country_name
                                logging.debug(f"Country common name is same as country name.")
                            state_fuzzy = pc.subdivisions.search_fuzzy(flattened_part)
                            state_name = state_fuzzy[0].name
                            logging.debug(f"State name is: {state_name}")
                            if state_fuzzy[0].country_code == country_fuzzy[0].alpha_2: # if the state name is in the same country
                                if len(state_name) > len(flattened_part):
                                    logging.debug(f"State belongs to the country. However, the query is a country name and due to fuzzy search is giving some random state search.")
                                else:
                                    logging.debug(f"State name is in the same country. Returning the country and state names.")
                                    state_code = state_fuzzy[0].code
                                    logging.debug(f"State code is: {state_code}")
                                    return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                            else: # if the state name is not in the same country return just the country name as the splitted_parts[0] is a country's name only
                                logging.debug(f"State doesn't belong to the country. Due to the fuzzy search, it has selected some state.")
                        
                        else: 
                            country_fuzzy = pc.countries.search_fuzzy(flattened_part)
                            country_name = country_fuzzy[0].name
                            logging.debug(f"Country name is: {country_name}")

                            country_three_lettered_name = country_fuzzy[0].alpha_3
                            logging.debug(f"Country three-lettered name is: {country_three_lettered_name}")
                            
                            try:
                                country_common_name = country_fuzzy[0].common_name
                                logging.debug(f"the common name for country: {country_common_name}")
                            except:
                                country_common_name = country_name
                                logging.debug(f"Country common name is same as country name.")
                            
                            
                            if country_name:
                                state_fuzzy = pc.subdivisions.search_fuzzy(flattened_part)
                                state_name = state_fuzzy[0].name
                                logging.debug(f"State name is: {state_name}")
                                if state_fuzzy[0].country_code == country_fuzzy[0].alpha_2: # if the state name is in the same country
                                    if len(state_name) > len(flattened_part):
                                        logging.debug(f"State belongs to the country. However, the query is a country name and due to fuzzy search is giving some random state search.")               
                                    else:
                                        logging.debug(f"State name is in the same country. Returning the country and state names.")
                                        state_code = state_fuzzy[0].code
                                        logging.debug(f"State code is: {state_code}")
                                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                                else: # if the state name is not in the same country return just the country name as the splitted_parts[0] is a country's name only
                                    logging.debug(f"State doesn't belong to the country. Due to the fuzzy search, it has selected some state.")
                                
                    except: # closing the try block for looking for the state name in the splitted_parts[0]
                        logging.debug(f"State name couldn't be found.")
                        continue
                        
            except: # closing the try block for looking for the country name in the splitted_parts[0]
                try: # trying to get the country name from the geolocator
                    logging.debug(f"Trying to get the country name from the geolocator.")
                    metageolocator = Nominatim(user_agent="metaminerGUI")

                    
                    address = metageolocator.geocode(fuzzy_country_names, addressdetails=True, language='en', exactly_one=True)
                    time.sleep(1)
                    logging.debug(f"Address is: {address}")

                    if address:
                        formatted_address = address.raw.get('address')
                        ad_country_name = formatted_address.get('country')
                        
                        country_name = pc.countries.search_fuzzy(ad_country_name)[0].name
                        logging.debug(f"Country name is: {country_name}")

                        country_three_lettered_name = pc.countries.search_fuzzy(country_name)[0].alpha_3
                        logging.debug(f"Country three-lettered name is: {country_three_lettered_name}")

                        try:
                            country_common_name = pc.countries.search_fuzzy(country_name)[0].common_name
                            logging.debug(f"Country common name is: {country_common_name}")
                        except:
                            country_common_name = country_name
                            logging.debug(f"Country common name is same as country name.")

                        
                        ad_state_name = formatted_address.get('state')
                        state_name = pc.subdivisions.search_fuzzy(ad_state_name)[0].name
                        if not state_name:
                            logging.debug("State name couldn't be found. Trying to get the county name.")
                            try:
                                state_name = formatted_address.get('county')
                                logging.debug(f"State name is: {state_name}")
                                try:
                                    state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                    logging.debug(f"State code is: {state_code}")
                                except:
                                    state_code = None
                                    logging.debug(f"State code couldn't be found.")
                            except:
                                state_name = formatted_address.get('island')
                                logging.debug(f"State name is: {state_name}")
                                try:
                                    state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                                    logging.debug(f"State code is: {state_code}")
                                except:
                                    state_code = None
                                    logging.debug(f"State code couldn't be found.")
                        
                        state_code = pc.subdivisions.search_fuzzy(state_name)[0].code
                        logging.debug(f"State code is: {state_code}")
                        return (country_name, country_common_name, country_three_lettered_name, state_name, state_code)
                    else:
                        return (None, None, None, None, None)

                except: # if the country name couldn't be found using the geolocator
                    logging.debug(f"Country couldn't be found using the geolocator. returning 'Unclassified' for both country and state.")
                    return (None, None, None, None, None)
        
        else:
            return (None, None, None, None, None)
    else:
        return (None, None, None, None, None)

def remove_accents(input):
    nfkd_form = uni.normalize('NFKD', input)
    return u"".join([c for c in nfkd_form if not uni.combining(c)])

def dict_update_from_new_loc(df:pd.DataFrame, dict:dict):

    logging.info(f"Updating the local database for normalization of geographical data.")
    length_of_records = len(df['geo_loc_name'].unique())
    # logging.debug(f"Length of df records to be normlized: {length_of_records}")
    counter = 1
    # logging.debug(f"counter set to 1.")
    for fuzzy_country_names in df['geo_loc_name'].unique():
        # logging.debug(f"Fuzzy country name: {fuzzy_country_names}")
        logging.info(f"processing: {counter}/{length_of_records}")
        if fuzzy_country_names is not None:
            # logging.debug(f"fuzzy country name is not None")
            if fuzzy_country_names in dict:
                logging.debug(f"Country name is already in the database.")
                counter += 1
                continue
            else:
                logging.debug(f"Country data is not in database. Looking for the country it can belong to!")
                dict[fuzzy_country_names] = get_country_name(fuzzy_country_names)
                counter += 1
        else:
            counter += 1
            continue    

def closing_function(country_state_map:dict, isolation_source_df:pd.DataFrame):
    try:
        updated_country_state_database = pd.DataFrame.from_dict(country_state_map, orient='index', columns=['country_name', 'country_common_name', 'country_three_lettered_name', 'state_name', 'state_code'])
        logging.debug(f"The updated country_state_map has been updated and transformed into dataframe.")

        updated_country_state_database.to_csv(os.path.abspath('./data/geo_loc_name_transforming_database.tsv'), sep='\t', index=True)
        isolation_source_df.to_csv(os.path.abspath('./data/categorized_isolation_sources.tsv'), sep='\t', index=False)
        logging.debug(f"Database updated and saved.")
    except Exception as e:
        logging.error(f"An error occurred while updating the database: {e}")
    
def categorize_sequencing_technology(given_value):

    illumina_platforms = ('illumina', 'miseq', 'hiseq', 'nextseq', 'novaseq', 'truseq', 'ilumina', 'slx')
    pacbio_platforms = ('pacbio', 'sequel', 'sequel-ii', 'rs-ii', 'hifi', 'hi-fi', 'rs', 'smrt', 'pacific', 'rsii')
    oxford_nanopore_platforms = ('nanopore', 'oxford', 'minion', 'gridion', 'promethion', 'flongle', 'ont')
    bgi_platforms = ('bgi', 'mgiseq', 'bgiseq', 'dnbseq', 'mgi', 'dnb-seq', 'dnbseq', 'dnb')
    thermo_fisher_platforms = ('thermo', 'ion', 'torrent', 'ion proton', 'pgm', 'studio', 'iontorrent')
    sanger_platforms = ('sanger', 'abi', 'capillary')
    roche_454_platforms = ('454', 'roche', 'gs', 'flx', 'titanium', 'junior')
    solid_platforms = ('solid', 'abi solid')
    genolab_platforms = ('genolab', 'genolab M', 'genemind', 'genemind genolab', 'genemind genolab M')

    try:
        matched_platforms = set()

        if given_value == None:
            return 'Unknown'
        elif given_value == '454' or given_value == 454:
            matched_platforms.add('Roche 454')
        
        if isinstance(given_value, str):
            part = given_value.lower().replace('(', ' ').replace(')', ' ')

            if '454' in part:
                matched_platforms.add('Roche 454')

            part = re.sub(r'\d+', '', part)
            parts = re.split(r'[-+:;/_,\s]|\band\b', part)
            parts = [part.strip() for part in parts]
            parts = [str(part) for part in parts]
            for part in parts:
                if part != '':
                    if part in illumina_platforms:
                        matched_platforms.add('Illumina')
                    if part in pacbio_platforms:
                        matched_platforms.add('PacBio')
                    if part in oxford_nanopore_platforms:
                        matched_platforms.add('Oxford Nanopore')
                    if part in bgi_platforms:
                        matched_platforms.add('MGI')
                    if part in thermo_fisher_platforms:
                        matched_platforms.add('Ion Torrent')
                    if part in sanger_platforms:
                        matched_platforms.add('Sanger')
                    if part in roche_454_platforms:
                        matched_platforms.add('Roche 454')
                    if part in solid_platforms:
                        matched_platforms.add('SOLiD')
                    if part in genolab_platforms:
                        matched_platforms.add('GenoLab')
            
            if len(matched_platforms) == 1:
                return matched_platforms.pop()
            elif len(matched_platforms) > 1:
                return ' and '.join(sorted(matched_platforms))
            else:
                return 'Unknown'
        else:
            return 'Unknown'
    except:
        return 'Unknown'

def find_the_host(fuzzy_str:str):
    logging.debug(f"Started! finding the host for string `{fuzzy_str}`")

    host_categories = {
        "Unknown" : ["na", "missing", "not provided", "none", "nan", "unknown", "no data", "missing", "none", "not applicable", "not recorded", "not available", "not found", "unspecified", "not specified", "not collected", "not determined",],
        
        "Hospital-associated" : ['human', 'homo sapiens', 'homo sapiens sapiens', 'man', 'woman', 'male', 'female', 'adult', 'child', 'infant', 'newborn', 'person', 'individual', 'human being', 'person of interest', 'human host', 'human subject', 'patient', 'patient sample', 'volunteer', 'test subject', 'human model', 'human sample', 'homo sapiens individual', 'homo sapiens sample', 'homo sapiens subject', 'human volunteer', 'medical personal', 'peramedical personal', 'Bed', 'hospital bed', 'adjustable bed', 'electric hospital bed', 'manual hospital bed', 'patient bed', 'Care Cart', 'nurse care cart', 'medical cart', 'hospital trolley', 'medication cart', 'equipment trolley', 'Wheelchair', 'hospital wheelchair', 'manual wheelchair', 'electric wheelchair', 'patient wheelchair', 'mobility chair', 'Hospital Stretcher', 'medical stretcher', 'emergency stretcher', 'patient stretcher', 'transport stretcher', 'Examination Table', "doctor's table", 'patient exam table', 'adjustable examination table', 'clinic exam table', 'Incubator', 'intensive care unit incubator', 'premature baby incubator', 'newborn care incubator', 'neonatal incubator', 'Gurney', 'hospital gurney', 'emergency gurney', 'transport gurney', 'patient transport gurney', 'soap dispenser', 'hand soap dispenser', 'automatic soap dispenser', 'disinfectant dispenser', 'disinfectant station', 'sanitizer dispenser', 'disinfectant bottle', 'sanitizing station', 'alcohol dispenser', 'disinfectant', 'hospital disinfectant', 'antiseptic', 'sanitizing agent', 'sterilizing solution', 'hand sanitizer', 'sanitizing wipes', 'cleaning wipes', 'antibacterial wipes', 'hospital wipes', 'surface sanitizing wipes', 'sterilizing station', 'disinfection station', 'medical sterilizing area', 'sanitization station', 'gloves', 'medical gloves', 'disposable gloves', 'latex gloves', 'non-latex gloves', 'sterile gloves', 'face mask', 'medical face mask', 'surgical mask', 'disposable face mask', 'N95 mask', 'respirator mask', 'ecg device', 'electrocardiogram device', 'ecg machine', 'electrocardiograph', 'heart monitor', 'cardiac monitor', 'pulse oximeter', 'blood oxygen monitor', 'oxygen saturation monitor', 'pulse ox', 'oximeter probe', 'pulse oximetry device', 'blood pressure cuff', 'sphygmomanometer', 'manual blood pressure cuff', 'automatic blood pressure cuff', 'bp cuff', 'patient monitor', 'heart rate monitor', 'syringe', 'syringe pump', 'syringe driver', 'insulin syringe', 'syringe for injections', 'sterile syringe', 'stethoscope', 'manual stethoscope', 'digital stethoscope', 'electronic stethoscope', 'doctor stethoscope', 'thermometer', 'digital thermometer', 'oral thermometer', 'infrared thermometer', 'temperature probe', 'thermometer gun', 'rectal thermometer', 'glucometer', 'blood glucose meter', 'diabetes meter', 'glucose monitor', 'diabetic glucose test', 'blood sugar meter', 'defibrillator', 'automated external defibrillator', 'aed', 'manual defibrillator', 'external defibrillator', 'cardiac defibrillator', 'AED pads', 'defibrillator electrodes', 'cardiac defibrillator pads', 'defibrillator patches', 'external defibrillator pads', 'endotracheal tube', 'et tube', 'intubation tube', 'artificial airway', 'tracheal tube', 'endotracheal intubation', 'ventilator', 'mechanical ventilator', 'respiratory ventilator', 'breathing machine', 'oxygen concentrator', 'ventilator machine', 'infusion pump', 'iv infusion pump', 'drug infusion pump', 'medication infusion pump', 'volume-controlled infusion pump', 'nebulizer', 'nebulizer machine', 'compressor nebulizer','nebuliser', 'breathing treatment device', 'inhalation therapy device', 'respiratory nebulizer', 'anesthesia machine', 'anaesthesia workstation', 'anesthesia delivery system', 'gas anesthesia machine', 'intubation anesthesia machine', 'surgical instruments', 'scalpel', 'surgical scissors', 'surgical forceps', 'needle holder', 'surgical clamps', 'laser therapy device', 'medical laser', 'laser therapy machine', 'diode laser', 'surgical laser device', 'laser wound healing machine', 'sterile field', 'sterile drape', 'sterile surgical field', 'surgical sterile area', 'sterilized surgical area', 'electrosurgical unit', 'electrocautery unit', 'diathermy machine', 'electrosurgical pencil', 'electrocoagulation unit', 'electrosurgical system', 'x-ray machine', 'x-ray equipment', 'x-ray imaging device', 'ultrasound machine', 'ultrasonography', 'portable ultrasound', 'diagnostic ultrasound', 'ultrasound scanner', 'ultrasound imaging system', 'ct scanner', 'computed tomography', 'ct scan', 'CAT scan', 'CT imaging system', '3D CT scanner', 'mri machine', 'magnetic resonance imaging', 'mri scanner', 'MRI equipment', 'MRI system', 'magnetic resonance scanner', 'endoscope', 'medical endoscope', 'fiber-optic endoscope', 'gastroscope', 'bronchoscope', 'colonoscope', 'mammography machine', 'breast x-ray machine', 'mammogram machine', 'digital mammography', 'mammography system', 'pet scanner', 'positron emission tomography scanner', 'PET imaging machine', 'PET scanner', 'nuclear medicine scanner', 'intravenous pole', 'iv pole', 'iv stand', 'intravenous drip stand', 'drip pole', 'iv support stand', 'intravenous infusion pole', 'iv infusion set', 'intravenous infusion set', 'IV drip set', 'infusion tubing', 'IV fluid administration set', 'infusion line', 'catheter', 'urinary catheter', 'foley catheter', 'iv catheter', 'peripheral catheter', 'central line catheter', 'blood bag', 'blood transfusion bag', 'blood collection bag', 'iv blood bag', 'donated blood bag', 'transfusion bag', 'autoclave', 'steam autoclave', 'pressure cooker sterilizer', 'uv sterilizer', 'ultraviolet sterilizer', 'uv sanitizing machine', 'disinfection UV device', 'UV-C sterilizer', 'ultraviolet light sterilizer', 'refrigerator', 'medical refrigerator', 'pharmaceutical refrigerator', 'vaccine refrigerator', 'cold storage', 'medicine fridge', 'hospital surface', 'hospital floor', 'hospital counter', 'bedside table', 'hospital furniture', 'infant incubator', 'neonatal incubator', 'baby incubator', 'newborn incubator', 'medical cooler', 'medicine cooler', 'cooling box', 'drug cooler', 'temperature-controlled medical cooler', 'oxygen concentrator', 'portable oxygen concentrator', 'oxygen generation unit', 'medical oxygen machine', 'dialysis machine', 'hemodialysis machine', 'kidney dialysis machine', 'dialysis equipment', 'renal dialysis unit', 'cryotherapy machine', 'cryosurgery unit', 'cryotherapy device', 'cold therapy machine', 'cryotherapy equipment', 'trauma kit', 'emergency trauma kit', 'first aid kit', 'trauma response kit', 'medical emergency kit', 'splint', 'splinting device', 'leg splint', 'arm splint', 'trauma splint', 'splint bandage', 'tourniquet', 'medical tourniquet', 'blood pressure tourniquet', 'hemostatic tourniquet', 'trauma tourniquet', 'ambulance stretcher', 'ambulance transport stretcher', 'field stretcher', 'mobile stretcher', 'emergency stretcher', 'nosocomial', 'nosocomial infections', 'bacterial infection', 'hospital infection', 'hospital', 'phc', 'primary health centre', 'primary hospital', 'secondary health centre', 'surgery', 'tertiary hospital', 'tertiary care', 'icu', 'intesive care units', 'emergency ward', 'ward', 'hospital ward', 'general ward', 'hospital environnment', 'sputum', 'hospital envrionment', 'clinical', 'clinical isolate', 'clinical material', 'Hospital ventilator', 'bedsheet', 'bed bar', 'bed side rail', 'bed rails', 'bed controllers', 'crash cart', 'emergency cart', 'Incubator swab', 'Incubator sample', 'alcohol', 'hospital alcohol', 'alcohol foam', 'cleaning equipment', 'hospital cleaning tools', 'Mopper', 'Electric mop', 'Swiffer', 'mop handle', 'Cleaning trolly', 'Cleaning cart', 'Cleaning bucket', 'Cleaning mop', 'Cleaning mop head', 'Cleaning mop handle', 'Cleaning mop bucket', 'ECG clip', 'Electrocardiogram clip', 'ecg monitor', 'endobronchial tube', 'bronchial tube', 'bronchial intubation', 'endobronchial intubation', 'bronchial intubation tube', 'endobronchial intubation tube', 'respiratory assistance device', 'respiratory support device', 'respiratory therapy device', 'respiratory machine', 'respiratory support machine', 'respiratory therapy machine', 'rad', 'hospital ventilator', 'ventilator shelf', 'tracheostomy tube', 'trach tube', 'tracheostomy tube holder', 'tracheostomy tube cuff', 'tracheostomy tube connector', 'trachy', 'Infusion stand', 'arterial catheter', 'cvc', 'central venous catheter', 'Foley Catheter', 'catheter container', 'catheter bag', 'nasogastric tube', 'NG tube', 'nasal tube', 'nasal catheter', 'crbsi', 'catheter associated bloodstream infection', 'Central Venous Access Devices', 'CVADs', 'Peripherally Inserted Central Catheters', 'PICC', 'Central Venous Catheters', 'CVCs', 'Tunneled Central Venous Catheters', 'shaldon catheter', 'refrigator handle', 'phone', 'dial device', 'landline phone', 'intrahospital phone line', 'Computers', 'Computer screens', 'TV', 'Television', 'Remote control', 'Keyboard', 'Computer mouse', 'Computer keyboard', 'Computer monitor', 'Computer CPU', 'Computer tower', 'Computer screen', 'Computer display', 'Switch button', 'hospital material', 'Intesive care units', 'icus', 'nicu', 'hospital drawer', 'door', 'sliding door', 'window', 'window door', 'door handle', 'tap', 'water tap', 'sink', 'sink inside', 'outside sink', 'sink drain', 'Washroom sink', 'hospital icu', 'hospital intensive care unit', 'Dialysate', 'hospital drain', 'sink countertop', 'hospital washbasin', 'washbasin', 'basin', 'hospital basin', 'hospital basin water', 'basin drainage water', 'hospital surface surveillance', 'crash trolley','hospital surface swab', 'hospital surface sample', 'hospital surface culture', 'hospital surface isolate', 'hospital surface cleaning', 'hospital surface disinfection', 'hospital surface sterilization', 'hospital surface wipe', 'hospital surface wipe sample', 'hospital surface wipe culture', 'surveillance', 'surface surveillance', 'surface swab', 'surface sample', 'surface culture', 'surface isolate', 'surface cleaning', 'surface disinfection', 'surface sterilization', 'surface wipe', 'surface wipe sample', 'surface wipe culture', "carious dentine", "sicu nursing station", "human body fluids", "human blood culture"],
        
        "Animal-associated" : ["poultry flies", "blood sausage", "duck incubator", "rodent feces", "cattle feces", "Chicken fecal", "cattle fecal", "Poultry flies", "Poultry shed", "Livestock feces", "pig feces", "wildlife feces", "dog feces", "Animal feces", "Capra", "lowland paca", "agouti lowland paca", "milvus migrans", "milvus lineatus lineatus", "pandion haliaetus", "black-eared kite", 'pork burger', 'Bactrian Camel', 'Non-Domesticated Animal', 'Eurycantha calcarata', 'Sus Scrofa Domesticus', 'Larus dominicanus', 'Oriental Cockroach', 'ox', 'african lion', 'patty', 'Hydrophilus piceus', 'Common Goldfish', 'E. testudinaria', 'Beluga Dolphin', 'Tribolium castaneum larva', 'Hog', 'Impala', 'venison steak', 'Bovine', 'pork hot dog', 'jellyfish tentacle', 'wet market', 'Eel meat', 'cnidaria', 'insect larva', 'Caribou Reindeer', 'Lygaeus equestris', 'turtle meat', 'chroicocephalus novaehollandiae', 'Vicugna Pacos', 'salmon', 'Water Animal', 'silk moth larva', 'Dromaius Novaehollandiae', 'White Whale', 'phascolarctos cinereus', 'mutton hot dog', 'Picidae Family', 'Ovis Aries', 'Cygnus Atratus', 'Boar meat', 'insect body part', 'fish', 'sumatran orangutan', 'white rhino', 'Albatross', 'Murine', 'rattus', 'Accipitridae family', 'Dor Beetle', 'Cow', 'Large Milkweed Bug', 'Skate meat', 'Pythonidae', "Grévy's Zebra", 'Wild Animals', 'Dairy Cow', 'Lucanus cervus larva', 'clam', 'lung', 'Phascolarctos Cinereus', 'veal burger', 'Vermilion Flycatcher', 'insect tarsal pulvilli', 'Phyllium bioculatum', 'giant squid', 'Orangutan', 'Green-Winged Teal', 'lamb pastrami', 'insect tarsal aroliae', 'Domestic Rat', 'Salamander meat', 'Working Ox', 'puma', 'Antelope meat', 'Green-winged Teal', 'lamb nugget', 'pacific salmon', 'crab meat', 'veal sausage roll', 'Lizard meat', 'Prawn', 'L. equestris', 'Periplaneta americana nymph', 'Snowy Owl', 'Halyomorpha halys', 'Corvus corax', 'Corvus frugilegus', 'Blackfish', 'insect body swab', 'mouse', 'C. campestris', 'Marine Creature', 'Dacelo novaeguineae', 'asteroidea', 'Scarabaeus sacer', 'raw chicken meat', 'turkey burger', 'Chimpanzee', 'lamb pepperoni', 'ass', 'Galleria mellonella larva', 'rhinoceros beetle grub', 'Meat/Organ', 'Steller Sea Lion', 'European Rhinoceros Beetle', 'Megaceryle alcyon', 'Setophaga Petechia', 'flying bird', 'Turtle meat', 'chicken salami', 'Green Lacewing', 'cockatoo', 'oreamnos americanus', 'Narwhal', 'Sialia currucoides', 'European Stag Beetle', 'Granary Weevil', 'insect egg', 'trash panda', 'G. mellonella', 'delphinapterus leucas', 'Morus Bassanus', 'insect thorax', 'Fauna', 'B. orientalis', 'bug', 'Rodent', 'Bear', 'decomposed larva', 'Sloth', 'Laridae family', 'mutton leg', 'sea otter', 'Great Ape', 'sperm whale', 'Ovine', 'Wolf', 'Slug meat', 'Serinus canaria', 'domestic yak', 'Snow Goose', 'P. interpunctella', 'spheniscus humboldti', 'insect patella', 'Llama Alpaca', 'Burro', 'Timber Wolf', 'Newt meat', 'pongo', 'A. cephalotes', 'earwig nymph', 'leafcutter ant larva', 'Crow', 'D. parallelipipedus', 'Sunbird', 'Rice Weevil', 'skunk', 'insect abdomen', 'Miscellaneous', 'Ephestia kuehniella', 'starfish arm', 'House Cat', 'Platalea Species', 'Guinea Pig', 'Indian Elephant', 'Humboldt penguin', 'rabbit', 'Sitophilus granarius larva', 'virginia opossum', 'Folivora', 'Crocodile meat', 'Pyrrhocoris apterus', 'Fawn', 'Silverback Gorilla', 'Bubo Virginianus', 'asiatic lion', 'Chilean Flamingo', 'starfish', 'Lion', 'saltwater fish', 'Golden Eagle', 'arctic fox', 'H. halys', 'Lobster', 'Brachyura', 'Cetacea', 'Arctic Wolf', 'Insect larva/Nymph', 'Tribolium castaneum', 'African Lion', 'Golden Pheasant', 'endangered species', 'Kangaroo', 'Conehead Mantis', 'mustela', 'Carassius Auratus', 'pink shrimp', 'pony', 'raccoon', 'Bronze Carabid', 'sea lion', 'Mouse', 'Donkey', 'leporidae', 'bison roast', 'mareca strepera', 'Hagfish meat', 'raw pork sausage', 'insect body isolate', 'Green Peach Aphid', 'Falco Tinnunculus', 'X. violacea', 'Bald-headed Ibis', 'Cuttlefish meat', 'impala', 'Poult', 'domestic cat', 'Tettigonia viridissima', 'Mountain Goat', 'Fruit Bat', 'chicken ham', 'pork bacon', 'Cyanocitta cristata', 'Southern Green Stink Bug', 'macropodidae', 'vampire bat', 'Canine', 'G. italicum', 'ladybug larva', "Boyd's forest dragon", 'chicken nugget', 'yak', 'Rat', 'Cuttlefish', 'insect tibia', 'Phoenicopterus Chilensis', 'A. nerii', 'fowl', 'Draught Ox', 'Bee', 'Ox', 'Phoenicopterus chilensis', 'Gromphadorhina portentosa', 'Cardinalis cardinalis', 'Culex pipiens', 'leopard', 'Woolly Apple Aphid', 'Setophaga petechia', 'Camel', 'box jellyfish', 'N. glauca', 'two-spotted cricket nymph', 'beetle larva', 'Bunny', 'Tomicus piniperda', 'domestic dog', 'grey heron', 'C. septempunctata', 'common goldfish', 'mole', 'Gannet', 'red fox', 'wildcat', 'green turtle', 'Mussel meat', 'cougar', 'European Spruce Bark Beetle', 'wolverine', 'King Of The Jungle', 'Anopheles gambiae larva', 'C. punctata', 'Female Bovine', 'mutton burger', 'insecta', 'Cattle', 'livestock', 'silkworm', 'manatee', 'Gavia Pacifica', 'coyote', 'Bombus terrestris', 'Macropodidae', 'arthropods', 'mountain goat', 'Harlequin Ladybird', 'Hermit Crab', 'Common Octopus', 'mallard', 'spinus spinus', 'Coccinella septempunctata larva', 'Amazonian Manatee', 'squirrel', 'Ram', 'Falco Peregrinus', 'Crowned Crane', 'Harmonia axyridis', 'turkey pastrami', 'chimp', 'bottlenose dolphin', 'Rose Chafer', 'Snake', 'asian elephant', 'hippopotamus', 'marine creature', 'Octopus', 'Cathartes aura', 'Caribou meat', 'chicken pastrami', 'farm animal', 'chicken wings', 'Unicorn Whale', 'cavia porcellus', 'E. tiaratum', 'Ostrich', 'lobster tail', 'Chroicocephalus ridibundus', 'Culex pipiens larva', 'cavy', 'Pinnipedia', 'Skink meat', 'american badger', 'Stork', 'Beef Cattle', 'camel', 'timber wolf', 'Aphis nerii', 'white grub', 'Grizzly Bear', 'orca', 'Armadillo', 'Phalacrocorax', 'cynomys', 'Casuarius casuarius', 'Fish', 'Calosoma sycophanta', 'Pan Troglodytes', 'red flour beetle larva', 'jackrabbit', 'Scallop meat', 'Giant Cave Cockroach', 'annelid', 'beef bologna', 'kangaroo', 'Domestic Yak', 'veal bologna', 'hen', 'Lucilia sericata larva', 'stag beetle larva', 'Nectariniidae family', 'spiny anteater', 'beef bacon', 'zebra', 'venomous snake', 'Hierodula membranacea', 'S. granarius', 'turkey pepperoni', 'camelus', 'Canada Goose', 'Mule Horse', 'Nine-Banded Armadillo', 'squid tube', 'green lacewing larva', 'Yoke Ox', 'T. viridissima', 'White-Tailed Deer', 'crocodile tail', 'lumbricina', 'Bison Bison', 'pork pastrami', 'bengal tiger', 'goose leg', 'Squirrel meat', 'lizard', 'Rangifer Tarandus', 'mutton meatball', 'Balearica Pavonina', 'Coccinella septempunctata', 'Equus Zebra', 'Gray Wolf', 'spotted skunk', 'Malaria Mosquito', 'great ape', 'Domestic Rabbit', 'Yellow Fever Mosquito', 'Dasypodidae', 'orcinus orca', 'Caprimulgidae family', 'Red-Crowned Crane', 'Colt', 'Corythucha ciliata', 'Puffin', 'Squid meat', 'stallion', 'stag', 'three-toed sloth', 'olive baboon', 'Calliphora vomitoria larva', 'Domestic Cat', 'panda', 'acinonyx jubatus', 'canis lupus familiaris', 'insect tarsal pad', 'Peruphasma schultei', 'hyena', 'Water Buffalo', 'Moon Jellyfish', 'C. pipiens', 'P. dominula', 'panthera tigris', 'Bubo scandiacus', 'bunny', 'Laridae Family', 'jaguar', 'Southern Armadillo', 'vulnerable species', 'insect leg', 'Insect', 'capra aegagrus hircus', 'Mahi-mahi', 'canis lupus', 'Bison meat', 'Others', 'Great Egret', 'Hare meat', 'kitten', 'Jellyfish', 'Tame Animals', 'European Hare', 'paper wasp larva', 'Reindeer meat', 'H. membranacea', 'earthworm', 'Woodcock', 'larva', 'felis catus', 'Kookaburra', 'Antelope', 'Farm Animals', 'H. axyridis', 'Threatened', 'american jaguar', 'sea urchin roe', 'Anser Caerulescens', 'Untamed Animal', 'Sun Beetle', 'mare', 'helarctos malayanus', 'Gavia immer', 'Chlamydosaurus kingii', 'F. auricularia', 'Wool Carder Bee', 'Monkey', 'Boar', 'Fowl', 'Peruvian Guinea Pig', 'Goose', 'Canis Lupus', 'beef ham', 'Ocean Fauna', 'Eremophila Alpestris', 'equus mulus', 'Agrilus planipennis', 'lobster', 'bugs', 'lamb chop', 'alpaca', 'Geronticus Eremita', 'marsupial', 'Corvus Brachyrhynchos', 'P. apterus', 'ailuropoda melanoleuca', 'Hamster', 'Iguana meat', 'Mountain Bluebird', 'insect trochanter', 'Eremophila alpestris', 'veal sausage', 'uncooked meat', 'crested gecko', 'Black Garden Ant', 'antelope', 'Cathartes Aura', 'gulo gulo', 'turkey sausage roll', 'bee', 'insect tarsal claw', 'Green Turtle', 'armadillo', 'Marsupial', 'cow', 'lamb salami', 'Eurygaster testudinaria', 'Ruffed Grouse', 'Atta cephalotes larva', 'Rock Goat', 'hammerhead shark', 'Clam meat', 'arthropod', 'European Mantis', 'Falco naumanni', 'koala', 'wolf', 'chick', 'veal meatball', 'fawn', 'Tettigonia viridissima nymph', 'Polistes dominula', 'Barnacle Goose', 'pork meatball', 'Gryllus bimaculatus nymph', 'Koala', 'Heifer', 'House Pet', 'European Starling', 'Companion', 'Loon', 'American Cockroach', 'Porcine', 'Sagittarius Serpentarius', 'Budgerigar', 'Siberian Tiger', 'insect proboscis', 'Black-Headed Gull', 'Pets', 'Equus Mulus', 'greater wax moth larva', 'I. typographus', 'A. planipennis', 'Buffalo Bison', 'ursus americanus', 'bubalus bubalis', 'bluebottle fly maggot', 'Vampire Bat', 'Southern Cassowary', 'rhinoceros', 'Nautilus meat', 'cetacea', 'scallop meat', 'Alpaca', 'Leptinotarsa decemlineata', 'tuscan sausage', 'C. nemoralis', 'blowfly maggot', 'Fur Seal', 'chicken hot dog', 'Bombus terrestris larva', 'violet carpenter bee larva', 'Acheta domesticus', 'Oryctes nasicornis larva', 'cricetinae', 'Mole', 'turkey jerky', 'Common Green Bottle Fly', 'chicken bologna', 'chiroptera', 'turkey ham', 'Feral Animal', 'Red-Tailed Hawk', 'mammal', 'Apis Mellifera', 'otter', 'Charadrius melodus', 'Toad meat', 'Seal', 'Diomedeidae family', 'Wandering Albatross', 'mephitidae', 'kitty', 'Vespa crabro', 'Chrysoperla carnea larva', 'Loxodonta Africana', 'Moose', 'H. piceus', 'S. oryzae', 'psittaciformes', 'giant anteater', 'poultry meat', 'Red Flour Beetle', 'Ursus Arctos', 'Megaceryle Alcyon', 'Myzus persicae', 'Bald-Headed Ibis', 'Scallop', 'pork sausage roll', 'black jaguar', 'Scolopax rusticola', 'malaria mosquito larva', 'Rodentia', 'caridea', 'insect chitin', 'insectoid', 'Sea Star', 'Leptinotarsa decemlineata larva', 'Syrian Hamster', 'Asian Longhorned Beetle', 'canis aureus', 'Cat', 'West Indian Manatee', 'brachyura', 'Haliaeetus leucocephalus', 'desert hare', 'Harmonia axyridis larva', 'tiger', 'ornithorhynchus anatinus', 'Forest Caterpillar Hunter', 'american lobster', 'Orcinus Orca', 'gorilla gorilla', 'Migratory bird', 'Whale', 'mosquito larva', 'Blatta orientalis', 'Polecat', 'Drosophila melanogaster larva', 'D. melanogaster', 'buffalo', 'lamb ham', 'L. niger', 'Yellow-billed Hornbill', "grévy's zebra", 'pantry moth larva', 'Sheep meat', 'chicken pepperoni', 'turkey', 'Armadillo meat', 'carassius auratus', 'teuthida', 'mussel meat', 'Geronticus eremita', 'Eudocimus ruber', 'insect femur', 'soft-shelled clam', 'Violet Carpenter Bee', 'Alces Alces', 'guinea pig', 'Sternidae Family', 'Thorny Devil Stick Insect', 'V. crabro', 'Dorcus parallelipipedus', 'indochinese tiger', 'dog', 'Bald Eagle', 'Pavo Cristatus', 'otariidae', 'eastern gray kangaroo', 'Farm Stock', 'Lesser Stag Beetle', 'Red-winged Blackbird', 'Crayfish meat', 'uncooked beef', 'Panda', 'canine', 'insect body organ', 'erethizontidae', 'Rock Lobster', 'Prawn meat', 'sloth', 'Oyster meat', 'chicken burger', 'primates', 'E. kuehniella', 'rock goat', 'Mediterranean flour moth larva', 'Yellow Warbler', 'tiger cat', 'Lucanus cervus', 'Bombyx mori larva', 'eastern chipmunk', 'felis silvestris catus', 'meles meles', 'beef sausage', 'B. terrestris', 'pet', 'Elephas Maximus', 'Corvus Corax', 'snake', 'Elk meat', 'X. dispar', 'kangaroo steak', 'orangutan', 'Raccoon meat', 'Myrmecophagidae', 'Red-breasted Nuthatch', 'tachyglossidae', 'Red Kangaroo', 'Pisces', 'ewe', 'buck', 'beef salami', 'Nightjar', 'Ips typographus', 'chick box', 'chicken breast', 'Lesser Water Boatman', 'Tockus flavirostris', 'L. cervus', 'sterile boot kit', 'shark', 'Falco peregrinus', 'Phoenicopterus roseus', 'turkey nugget', 'Chrysolophus pictus', 'african cheetah', 'harlequin ladybird larva', 'Dytiscus marginalis', 'Giraffa Camelopardalis', 'Cyanocitta Cristata', 'alces alces', 'Alfalfa Leafcutter Bee', 'equus asinus', 'Great Green Bush-Cricket', 'White-tailed Kite', 'A. manicatum', 'Dacelo Novaeguineae', 'Caridea', 'cricket', 'Elanus leucurus', 'Red Wolf', 'amphibian', 'Starfish', 'Echinoderm', 'pig', 'G. bimaculatus', 'chipmunk', 'Struthio camelus', 'african buffalo', 'dasypodidae', 'A. glabripennis', 'undercooked meat', 'lamb bacon', 'grizzly bear', 'tobacco hornworm', 'columbidae', 'Ardea Herodias', 'Madagascar Hissing Cockroach', 'bivalvia', 'stoat', 'Gallus gallus', 'raw meat', 'duck breast', 'Black-and-Red Bug', 'Chameleon meat', 'silverback gorilla', 'Horse', 'Coturnix Coturnix', 'Salmo Salar', 'Eriosoma lanigerum', 'Salami', 'Bivalvia', 'L. decemlineata', 'Red-crowned Crane', 'humpback whale', 'american bison', 'pork chop', 'fruit fly larva', 'Black Carpenter Ant', 'Squirrel Monkey', 'Western Honeybee', 'Ewe', 'Emerald Ash Borer', 'Camelus', 'red kangaroo', 'Belted Kingfisher', 'Domestic Chicken', 'Flamingo', 'Land Iguanas', 'Common Fruit Fly', 'Piranga ludoviciana', 'Chroicocephalus Ridibundus', 'chimpanzee', 'gray wolf', 'bornean orangutan', 'Muskrat meat', 'turtle', 'Manatee', 'Swallow', 'Livestock', 'swine liver', 'Oncopeltus fasciatus', 'gazelle', 'Ara macao', 'Morus bassanus', 'Water-Dwelling Animal', 'Dendrocygna viduata', 'Waxwing', 'Harpia Harpyja', 'Elanoides forficatus', 'monodon monoceros', 'Nezara viridula', 'beef pepperoni', 'Lepidoptera larva', 'C. auratus', 'lamb bologna', 'macaw', 'Sturnella magna', 'talpidae', 'turkey salami', 'mutton ham', 'beluga whale', 'Monitor meat', 'sciuridae', 'insect nymph', 'Puppy', 'castor canadensis', 'Caiman', 'Spiny Lobster', 'Blatta orientalis nymph', 'Sea Cow', 'decomposed grub', 'tamias', 'horse', 'Sardines', 'Turtle', 'pork ham', 'Monodon Monoceros', 'Empusa pennata', 'Manduca sexta larva', 'gerbil', 'mule', 'Stag', 'Swine', 'european elk', 'Aedes aegypti', 'insect body fluid', 'Plodia interpunctella', 'Freshwater Shrimp', 'brush wolf', 'Rock Pigeon', 'Eastern Gray Kangaroo', 'Duck', 'Miscellaneous Birds', 'Struthio Camelus', 'Indian meal moth larva', 'Ramphastos toco', 'Chinook Salmon', 'wings', 'narwhal', 'panthera pardus', 'house cricket nymph', 'Ape', 'D. marginalis', 'insect exoskeleton', 'Camponotus pennsylvanicus', 'Stallion', 'T. castaneum', 'monkey', 'rat snake', 'fox', 'Colorado potato beetle larva', 'beef hot dog', 'correlophus ciliatus', 'elk chop', 'glutton', 'Cervidae', 'Peregrine Falcon', 'Dairy Cattle', 'hyaenidae', 'rose chafer', 'california sea lion', 'moon jellyfish', 'Wild', 'polar bear', 'Elanoides Forficatus', 'insect body parts', 'Kitten', 'Snowshoe Hare', 'C. sycophanta', 'Rabbit', 'King Crab', 'Venison meat', 'Toco Toucan', 'P. prasina', 'Jack', 'Musca domestica larva', 'Pig meat', 'rangifer tarandus', 'Chrysoperla carnea', 'Caprimulgidae Family', 'Giant Asian Mantis', 'Acheta domesticus nymph', 'Quail', 'golden jackal', 'Heron', 'Musca domestica', 'Opossum meat', 'turkey meatball', 'Xylocopa violacea larva', 'Insects', 'Pyrrhalta luteola', 'Apis mellifera larva', 'Teuthida', 'Domestic Dog', 'Three-Toed Sloth', 'Camelid', 'Retriever', 'Poodle', 'cuttlefish steak', 'Central bearded dragon', 'Zebra', 'Killer Whale', 'quail breast', 'Beluga Whale', 'Kingfisher', 'Fratercula Arctica', 'clouded leopard', 'Sheep', 'king crab', 'puma concolor', 'Buck', 'Serinus Canaria', 'Domesticated Animals', 'Beef Cow', 'Canis Lupus Familiaris', 'Ciconia Family', 'sus scrofa domesticus', 'Conservation-Dependent', 'ground beef', 'Tyto alba', 'Ciconia family', 'lophosaurus boydii', 'Pyrocephalus rubinus', 'Tortoise Bug', 'pork salami', 'ape', 'Wild Animal', 'Bald Ibis', 'insect pupa', 'rattle snake', 'Xylocopa violacea', 'Madagascar hissing cockroach nymph', 'okapi', 'Ocean Creature', 'Buteo jamaicensis', 'dwarf hamster', 'insect maxilla', 'Bos Grunniens', 'Colorado Potato Beetle', 'Wildlife', 'Hammerhead Shark', 'Sausage', 'Anas crecca', 'Anteater', 'Aquatic Animals', 'Bowerbird', 'Tenebrio molitor', 'veal pepperoni', 'moth caterpillar', 'mustelidae', 'fruit bat', 'striped skunk', 'maggot', 'Primates', 'Great Diving Beetle', 'Indian Meal Moth', 'pig production settings', 'Bos Taurus', 'pork', 'M. rotundata', 'Anser caerulescens', 'Giant Panda', 'water animal', 'Lamb meat', 'hen farm', 'freshwater fish', 'Chiroptera', 'Mountain Gorilla', 'Anopheles gambiae', 'Cuculidae Family', 'myrmecophagidae', 'Sea Unicorn', 'shrimp', 'Scolopacidae family', 'giraffa camelopardalis', 'Frilled-neck lizard', 'bee brood', 'Pyrocephalus Rubinus', 'iguana tail', 'female lion', 'beef', 'ovis aries', 'Donkey-Horse Hybrid', 'Leporidae', 'lamb hot dog', 'Asiatic Lion', 'insect palpus', 'Ardeidae Family', 'Pachnoda marginata', 'Mus musculus', 'Asteroidea', 'white stork', 'frog legs', 'Poultry meat', 'Selachimorpha', 'Melopsittacus Undulatus', 'insect ovipositor', 'Pet Rabbit', 'L. sericata', 'Equus Ferus Caballus', 'Mantle Squid', 'Ardea herodias', 'Hare', 'sheep', 'deer', 'European Hornet', 'didelphimorphia', 'green bottle fly maggot', 'hippo', 'Diomedeidae Family', 'okapia johnstoni', 'avian', 'seal', 'insect tarsal arolia', 'european wildcat', 'Sitta canadensis', 'beef pastrami', 'Extatosoma tiaratum', 'worm', 'Dwarf Hamster', 'felis silvestris', 'killer whale', 'turkey sausage', 'Gerris lacustris', 'hard-shelled clam', 'Dromedary', 'Nutria meat', 'Pet', 'milk snake', 'black garden ant larva', 'rock pigeon', 'Fratercula arctica', 'Family Pets', 'Sus Scrofa', 'hedgehog', 'Bluebottle Fly', 'Brown Marmorated Stink Bug', 'antilopinae', 'Orca', 'insect body sample', 'Shrew', 'Cnidaria', 'poultry production setting', 'Alligator meat', 'Beaver meat', 'insect tarsal empodium', 'P. viburni', 'Porcupine meat', 'whale', 'Four-toed hedgehog', 'american black bear', 'Jackrabbit', 'platypus', 'Wren', 'insect body tissue', 'Bull', 'insectum', 'trichechus', 'Tamandua', "Man O' War", 'ground turkey', 'mutton pastrami', 'Branta leucopsis', 'goat', 'prawn', 'Dog', 'megaptera novaeangliae', 'Bucerotidae family', 'Northern Flicker', 'crab', 'American cockroach nymph', 'American Bison', 'Golden Hamster', 'Common Raven', 'Haliaeetus Leucocephalus', 'Trichechus', 'elephas maximus', 'Rook', 'mutton bacon', 'Black-headed Gull', 'Carabus auratus', 'Felis Catus', 'european mole', 'Plains Zebra', 'reindeer', 'Masai Giraffe', 'Elm Leaf Beetle', 'cervidae', 'Doe', 'Branta Canadensis', 'Blue Whale', 'house sparrow', 'Rooster', 'elk', 'insect antenna', 'Anas platyrhynchos', 'Squid', 'equus ferus caballus', 'Red-tailed Hawk', 'butterfly caterpillar', 'Barnacle meat', 'Turdus migratorius', 'Marine Life', 'mutton pepperoni', 'Lasius niger larva', 'Agelaius Phoeniceus', 'wapiti', 'Ardeidae family', 'mutton bologna', 'Whistling Duck', 'Golden Ground Beetle', 'Palomena prasina', 'non-venomous snake', 'Greater Wax Moth', 'insectae', 'Domestic Ferret', 'Bucerotidae Family', 'Viburnum Leaf Beetle', 'Cylindera germanica', 'bovine', 'mutton salami', 'B. giganteus', 'Flying Mammal', 'Eastern Meadowlark', 'Tiger', 'octopus', 'Clam', 'folivora', 'parrot', 'Pink Shrimp', 'Accipitridae Family', 'Alcedinidae Family', 'Deer', 'Forficula auricularia nymph', 'white whale', 'silver gull', 'Lucilia sericata', 'ground chicken', 'cormorant', 'Common Backswimmer', 'syrian hamster', 'Plodia interpunctella larva', 'House Cricket', 'Frilled Lizard', 'Asian Buffalo', 'chicken production setting', 'Anoplophora glabripennis', 'chicken bacon', "Rothschild's Giraffe", 'Geotrupes stercorarius', 'domestic alpaca', 'insect body specimen', 'mutton patty', 'Columba livia', 'Gryllus bimaculatus', 'Great White Shark', 'beef patty', 'veal hot dog', 'vicugna pacos', 'european badger', 'beauty rat snake', 'Deer meat', 'Ardea alba', 'ram', 'Feline', 'Goat meat', 'Turdus Migratorius', 'Great Blue Heron', 'S. sacer', 'Canary', 'doe', 'tenderized squid roll', 'Pogona vitticeps', 'bison bison', 'common raccoon', 'Axolotl meat', 'Melopsittacus undulatus', 'swine', 'Sycamore Lace Bug', 'yoke ox', 'blue crab', 'honeybee larva', 'Eurasian siskin', 'possum', 'choanal', 'Corixa punctata', 'Mutton Sheep', 'animal', 'house cat', 'wildlife', 'aphid lion', 'black-tailed prairie dog', 'Colaptes auratus', 'insect wing', 'Giant Squid', 'european hedgehog', 'A. gambiae', 'lamb burger', 'giant panda', 'giant pacific octopus', 'canis latrans', 'Oryctes nasicornis', 'Anthidium manicatum', 'European Elk', 'Kelp gull', 'mink', 'Giant Spiny Stick Insect', 'Crab', 'Giant Pacific Octopus', 'Ass', 'Crab meat', 'anteater', 'Polistes dominula larva', 'Grus Japonensis', 'striped hyena', 'insect body', 'selachimorpha', 'undercooked beef', 'Burro Donkey', 'Grus japonensis', 'boot swab kit', 'Lesser Kestrel', 'Bamboo Bear', 'bat', 'Sitophilus oryzae', 'desert rat', 'steed', 'Domestic Alpaca', 'Freshwater Fish', 'Endangered Wildlife', 'cattle', 'pan troglodytes', 'Northern Cardinal', 'Aedes aegypti larva', 'African Buffalo', 'Atlantic Salmon', 'Agelaius phoeniceus', 'Vulnerable', 'elephant', 'Pony', 'Cattle Stock', 'Earth Mole', 'gecko', 'female bovine', 'forest giraffe', 'conolophus subcristatus', 'Two-Toed Sloth', 'house fly maggot', 'chicken sausage roll', 'Drosophila melanogaster', 'Western Gray Kangaroo', 'lioness', 'Secretary Bird', 'atelerix albiventris', 'spiny lobster', 'mutton sausage roll', 'goose', 'donkey', 'threatened species', 'norway rat', 'Balearica pavonina', 'waxworm', 'African Elephant', 'Parrot', 'insect', 'grub', 'Aquatic Animal', 'feline', 'European Mole', 'ferret', 'Octopus meat', 'Mantis religiosa', 'Columba Livia', 'Meleagris gallopavo', 'insects', 'Great Silver Water Beetle', 'yellow baboon', 'insect body culture', 'prairie dog', 'California Sea Lion', 'Alces Alces Americanus', 'Water Creatures', 'Bottlenose Dolphin', 'Woodpecker', 'Companion Animal', 'Lamb', 'cat', 'Cetonia aurata', 'Gazelle', 'Harpy Eagle', 'Domestic Goat', 'M. domestica', 'oyster meat', 'German Tiger Beetle', 'Aquila Chrysaetos', 'A. aegypti', 'beef burger', 'Mussel', 'giraffe', 'Horned Lark', 'Whale Shark', 'Shark', 'Turkey Vulture', 'whistling duck', 'Farming Animals', 'Fish meat', 'Ray meat', 'beef nugget', 'pork bologna', 'Leafcutter Ant', 'capuchin monkey', 'Platalea species', 'Panthera Leo', 'insect isolate', 'Oreamnos Americanus', 'insect coxa', 'Koala Bear', 'feral animal', 'cervus canadensis', 'chicken', 'Husbandry Animal', 'lampropeltis triangulum', 'Moose meat', 'Tern', 'Sagittarius serpentarius', 'Gorilla', 'mustela putorius furo', 'Cuculidae family', 'Gadwall', 'lepus californicus', 'Trochilidae family', 'Casuarius Casuarius', 'Periplaneta americana', 'European Paper Wasp', 'gallus gallus domesticus', 'Pigeon', 'black bear', 'mealworm', 'Antilopinae', 'Gerbil', 'Falco tinnunculus', 'rhinocerotidae', 'Nature', 'cattle stock', 'Chrysolophus Pictus', 'burro', 'Elephant', 'Cavia Porcellus', 'Talpidae', 'Rabbit meat', 'Mule', 'Steed', 'insect claw', 'Hornbill', 'Elk (In Europe)', 'Giraffe', 'octopus tentacle', 'penguin', 'mutton sausage', 'house lizard', 'Tyto Alba', 'ciconia ciconia', 'squirrel monkey', 'insect cuticle', 'koala bear', 'water buffalo', 'Golden Fish', 'rat', 'Apis mellifera', 'Mealworm Beetle', 'domestic rabbit', 'lamb sausage roll', 'Snail meat', 'octopoda', 'caribou', 'hog', 'insect head', 'Cephalopod meat', 'companion animal', 'Seagull', 'organ', 'hippopotamus amphibius', 'north american porcupine', 'C. ciliata', 'American Robin', 'ostrich fillet', 'T. piniperda', 'Ephestia kuehniella larva', 'broiler chicken', 'Octopoda', 'Virginia opossum', 'Humpback Whale', 'Diomedea exulans', 'Sperm Whale', 'papio', 'Gallus Gallus Domesticus', 'Reticulated Giraffe', 'trachea', 'aquatic animal', 'Hummingbird', 'veal cutlet', 'bos grunniens', 'bos taurus', 'Otariidae', 'guinea fowl', 'Merops apiaster', 'insect labium', 'mealworm larva', 'Scolopacidae Family', 'intestine', 'chicken meatball', 'turkey hot dog', 'Sitophilus granarius', 'Ptilonorhynchidae Family', 'atlantic salmon', 'House Pets', 'Tuscan sausage', 'Dendrocygna Viduata', 'baboon', 'E. lanigerum', 'Gromphadorhina portentosa nymph', 'sea cow', 'donkey-horse hybrid', 'field mouse', 'house mouse', 'pig production setting', 'M. persicae', 'T. molitor', 'badger', 'beef meatball', 'M. melolontha', 'House Fly', 'fish meat', 'poultry', 'Companion Animals', 'Hog meat', 'Buffalo meat', 'ground pork', 'procyon lotor', 'Hen', 'Frog meat', 'Scarlet Ibis', 'Sitophilus oryzae larva', 'pork pepperoni', 'domestic animal', 'Bengal Tiger', 'Equine', 'equine', 'European hornet larva', 'north american beaver', 'Mountain Zebra', 'bison', 'O. nasicornis', 'domestic pet', 'Nephropidae', 'Soricidae', 'Polar Bear', 'Goat', 'beef steak', 'Oleander Aphid', 'Pine Shoot Beetle', 'caterpillar', 'Sloth Bear', 'insect tarsal arolium', 'domestic ferret', 'caiman', 'Indian Peafowl', 'C. pennsylvanicus', 'rabbit leg', 'Bison', 'bear', 'great white shark', 'Firebug', 'squid', 'seven-spotted ladybird larva', 'Working Animals', 'Brown Bear', 'insect integument', 'unicorn whale', 'Tenebrio molitor larva', 'veal pastrami', 'Sialia Currucoides', 'G. stercorarius', 'C. germanica', 'Western Tanager', 'Ursidae', 'weasel', 'Peacock', 'mus musculus', 'vertebrate', 'P. marginata', 'Green Tiger Beetle', 'dromedary', 'nine-banded armadillo', 'Family Pet', 'farm stock', 'Chicken', 'oriental cockroach nymph', 'beaver', 'Domestic Pet', 'Carabus nemoralis', 'Troglodytidae family', 'N. viridula', 'Himalayan Yak', 'Capra Aegagrus Hircus', 'panthera onca', 'Poultry', 'Melolontha melolontha larva', 'loxodonta africana', 'snake meat', 'mountain lion', 'Leaf Insect', 'Sow', 'pisces', 'meriones unguiculatus', 'pigeon', 'cockchafer grub', 'pork nugget', 'Xyleborus dispar', 'Spoonbill', 'porcupine', 'Mare', 'sun bear', 'tabby', 'boot swab', 'Dorcus parallelipipedus larva', 'Goldfish', 'delphinidae', 'C. vomitoria', 'turkey bologna', 'sus scorfa', 'Troglodytidae Family', 'Harbor Seal', 'Emu', 'goldfish', 'moose', 'veal bacon', 'Bubo virginianus', 'snow leopard', 'ursidae', 'Titmouse', 'Filly', 'guinea hen', 'hornworm', 'turkey drumstick', 'Ardea Alba', 'elephant seal', 'Domestic Animals', 'P. bioculatum', 'Red-Winged Blackbird', 'Shrimp meat', 'salmo salar', 'Pyrrhalta viburni', 'Feedlot Water Bowl', 'Mice', 'Syncerus Caffer', 'Anthus species', 'Common Pond Skater', 'rodent', 'Notonecta glauca', 'hamster', 'didelphis virginiana', 'side-striped jackal', 'Pig', 'Alcedinidae family', 'ursus arctos horribilis', 'erinaceidae', 'ardea cinerea', 'Aquatic', 'Eurasian Kestrel', 'Blaberus giganteus', 'Cardinalis Cardinalis', 'Picidae family', 'Coleoptera larva', 'Elephant Seal', 'chicken sausage', 'Endangered', 'Lobster meat', 'Reindeer', 'liver', 'river otter', 'Bat', 'Bubo Scandiacus', 'fish meal', 'malayan sun bear', 'Sacred Scarab', 'prairie wolf', 'wild animal', 'wet market samples', 'Gorilla Gorilla', 'House Sparrow', 'P. americana', 'yellow fever mosquito larva', 'insect tarsus', 'Camponotus pennsylvanicus larva', 'veal ham', 'beef sausage roll', 'Octopus Vulgaris', 'Mustela Putorius Furo', 'Common House Mosquito', 'meat', 'jellyfish', 'cheetah', 'Dolphin', 'insect adult', 'Giant Anteater', 'Coryphaena hippurus', 'cuttlefish', 'Anas Crecca', 'Macaw', 'black rat', 'C. aurata', 'annelida', 'G. lacustris', 'Sturnus vulgaris', 'Malayan Tiger', 'stork', 'Barn Owl', 'Ptilonorhynchidae family', 'P. luteola', 'Mediterranean Flour Moth', 'Seven-Spotted Ladybird', 'Crustacean meat', 'lamb meatball', 'gorilla', 'alligator mississippiensis', 'Gavia Immer', 'turkey bacon', 'green sea turtle', 'Yak', 'clam meat', 'Desert Camel', 'brown bear', 'Mustelid', 'feedlot water bowl', 'tame animal', 'Unknown', 'breaching whale', 'great green bush-cricket nymph', 'G. portentosa', 'Pacific Salmon', 'Capuchin Monkey', 'Calliphora vomitoria', 'Falco Naumanni', 'avian incubator', 'avian environmental incubator', 'farm bed', 'farm maternity bed', 'Cavy', 'Caribou', 'Cicindela campestris', 'Gecko meat', 'Delphinapterus Leucas', 'two-toed sloth', 'Buffalo', 'tree squirrel', 'pinnipedia', 'raw beef', 'annelids', 'Atta cephalotes', 'rice weevil larva', 'pork patty', 'Vulture', 'Panthera Tigris', 'choanal swab', 'P. schultei', 'Squid roll', 'Blue Crab', 'Two-Spotted Cricket', 'carpenter ant larva', 'Springbok', 'Bonasa umbellus', 'Great Horned Owl', 'insect mandible', 'wiggler', 'jackal', 'Ferret', 'Bear meat', 'Bubalus Bubalis', 'Bombycilla species', 'harbor seal', 'Chimp', 'Cricetinae', 'tabby cat', 'European Bee-eater', 'Passer domesticus', 'veal patty', 'Graphosoma italicum', 'Aquatic Fish', 'Baeolophus family', 'Blue Jay', 'Harpia harpyja', 'Coturnix coturnix', 'granary weevil larva', 'alligator leg', 'lamb sausage', 'lion', 'lesser stag beetle larva', 'chicken patty', 'Pipit', 'Branta canadensis', 'Black Bear', 'Sturnus Vulgaris', 'nephropidae', 'Indochinese Tiger', 'Melolontha melolontha', 'echidna', 'poultry production settings', 'M. religiosa', 'Shrimp', 'Megachile rotundata', 'syncerus caffer', 'insectoids', 'boar ribs', 'Salmon', 'dolphin', 'Ciconia ciconia', 'ground squirrel', 'Asian Elephant', 'spotted hyena', 'Buteo Jamaicensis', 'Hirundinidae family', 'Cuckoo', 'Black Beauty Stick Insect', 'Sandpiper', 'Colaptes Auratus', 'Soft-Shelled Clam', 'Agricultural Animal', 'Lasius niger', 'equus zebra', 'Lepus', 'Forficula auricularia', 'American Lobster', 'turkey patty', 'Lamprey meat', 'panthera leo', 'Hard-Shelled Clam', 'bactrian camel', 'Swallowtail Kite', 'Oryx', 'blue whale', 'Marine Wildlife', 'C. carnea', 'Corvus brachyrhynchos', 'Piping Plover', 'European Shot Hole Borer', 'Cavy Pig', 'E. calcarata', 'shrimp meat', 'reticulated giraffe', 'Macaque', 'Mollusk meat', 'city pigeon', 'Sternidae family', 'Equus Asinus', 'Buff-Tailed Bumblebee', 'Black Swan', 'Saltwater Fish', 'E. pennata', 'nymph', 'A. mellifera', 'veal salami', 'Branta Leucopsis', 'black rhino', 'Pacific Loon', 'plains zebra', 'Green Shield Bug', 'buff-tailed bumblebee larva', 'Box Jellyfish', 'Sea Lion', 'bird', 'Ara Macao', 'Pavo cristatus', 'Ailuropoda Melanoleuca', 'elaphe taeniura', 'Common Earwig', 'Italian Striped Bug', 'Vespa crabro larva', 'lamb patty', 'Cygnus atratus', 'Green Sea Turtle', 'Wild Goat', 'Snake meat', 'west indian manatee', 'Common Cockchafer', 'Aquila chrysaetos', 'insectums', 'pork sausage', 'Blind Mole', 'sea star', 'Shark meat', 'Galleria mellonella', 'A. domesticus', 'duck-billed platypus', 'small chicks', 'vulpes vulpes', "Fox", "Vulpes Vulpes", "Red Fox", "Arctic Fox", "Fennec Fox", "Swift Fox", "Kit Fox", "Gray Fox", "Corsac Fox", "Bengal Fox", "Blanford's Fox", 'Gavia pacifica', 'Bonobo', 'Mutton meat', 'Delphinidae', 'broiler', 'O. fasciatus', 'horse feces', 'black-tufted marmoset', 'callithrix penicillata', "columbiformes", "columba livia", "rock pigeon", "columba palumbus", "common wood pigeon", "columba oenas", "stock dove", "columba livia domestica", "domestic pigeon", "columba palumbus", "common wood pigeon", "columba oenas", "stock dove", "anas", "duck", "Anatidae", "Anatidae Family", "Anas Platyrhynchos", "Mallard Duck", "Anas Crecca", "Green-Winged Teal", "Anas Penelope", "Eurasian Wigeon", "Anas Acuta", "Northern Pintail", "Anas Clypeata", "Northern Shoveler", "Anas Discors", "Blue-Winged Teal", "Anas Americana", "American Black Duck", "organic drumett", "marine ascidian", "beef enrichment", "cricket powder", "cricket", "oncorhynchus mykiss", "Rainbow trout", "Salmo trutta", "Brown trout", "Salvelinus fontinalis", "Brook trout", "Oncorhynchus mykiss", "Rainbow trout", "Salmo salar", "Atlantic salmon", "Salvelinus namaycush", "Lake trout", "Micropterus salmoides", "Largemouth bass", "Micropterus dolomieu", "Smallmouth bass", "salvelinus namaycush", "lake trout", "salmo trutta", "brown trout", "salmo salar", "atlantic salmon", "oncorhynchus mykiss", "rainbow trout", "carassius auratus", "goldfish", "coryphaena hippurus", "mahi-mahi", "black bass", "oncorhynchus mykiss", "rainbow trout", "marmota himalayana", "himalayan marmot", "marmot", "cooked turkey stuffed with mushroms and chestnuts (pork in the meat)", "scarlet snapper", "lutjanus campechanus", "red snapper", "lutjanus campechanus", "lutjanidae", "snapper", "lutjanidae family", "lutjanus campechanus", "red drum", "sciaenops ocellatus", "redfish", "sciaenidae", "sciaenidae family", "sciaenops ocellatus", "frozen meat minced", "bream", "sparidae", "sparidae family", "sardine", "clupeidae", "clupeidae family", "sardina pilchardus", "pilchard", "sardina pilchardus", "pilchardus", "sardina pilchardus", "pilchardus", "poultry flock base material and dust", "Poultry barn", "periwinkle snail", "myodes glareolus", "bank vole", "clethrionomys glareolus", "clethrionomys", "clethrionomys sp", "clethrionomys sp.", "Field vole", "Vole", "drosophila", "melipona lateralis", "stingless bee", "melipona", "stingless bees", "hymenoptera", "hymenopteran", "hymenopterans", "hymenopteroid", "hymenopteroids", "hymenopterous", "hymenoptera insects", "hymenoptera insect", "melipona seminigra", "scaptotrigona polysticta", "scaptotrigona", "apidae", "apidae bees", "apidae bee", "hymenoptera insects", "benjoi", "bijui", "frieseomelitta varia", "frieseomelitta", "abelha marmelada-amarela", "yellow marmalade bee", "melipona interrupta", "jandaira", "perch", "frozen tuna saku aaa grade", "tuna", "frozen tuna", "frozen tuna fish", "frozen tuna fillet", "tuna fish", "longaniza", "spanish sausage", "embutido", "eidolon helvum", "africal fruit bat", "Golden lion tamarin", "Callitrichidae", "mico", "marmosets", "rectal swab of live bovine cattle", "bos indicus", "zebu", "animal non-clinical, bovine feces", "wombat", "vombatidae", "ovine", "seafood", "duck eggs", "carcass cleaning wipe", "boneless pork picnic", "dairy equipment filters", "meat processing equipment", "poultry processing equipment", "tilapia", "oreochromis sp.", "sphynx cat", "phocarctos hookeri", "oncorhynchus mykiss", "oncorhynchus kisutch", "diseased pig", "oreochromis niloticus", "nile tilapia", "cleaned poultry house", "skink", "equine hospital environmental", "insect", "planococcus ficus", "vine mealybug", "apodemus flavicollis", "yellow-necked field mouse", "apodemus sylvaticus", "wood mouse", "apodemus agrarius", "striped field mouse", "apodemus chevrieri", "chevrier's field mouse", "apodemus speciosus", "japanese wood mouse", "apodemus uralensis", "ural field mouse", "apodemus peninsulae", "korean field mouse", "apodemus agrarius", "striped field mouse", "apodemus chevrieri", "chevrier's field mouse", "apodemus speciosus", "japanese wood mouse", "apodemus uralensis", "ural field mouse", "apodemus peninsulae", "korean field mouse", "apodemus sylvaticus", "wood mouse", "apodemus flavicollis", "yellow-necked field mouse", "labeo rohita", "rodent deer mouse", "deer mouse", "canine endotracheal tube site", "feline endotracheal tube site", "Pied oystercatche", "Haematopus longirostris", "oystercatcher", "sooty oystercather", "haemotopus sp.", "grass carp", "Ctenopharyngodon idell", "black carp", "silver carp", "ctemopharyngodon sp.", "minced chicken chilled", "chicken thigh chilled", "chicken drumstick chilled", "lagria villosa", "chicharron with chili", "cygnus olor", "mute swans", "eagle", "broiler fillet", "lagria grenieri", "lagria atripes", "ecnolagria sp.", "lagria rufipennis", "lagria okinawana", "lagria", "lagria sp.", "mirounga leonina", "leptonychotes weddellii", "weddell seal", "eptesicus fuscus", "big brown bat", "lasionycteris noctivagans", "silver-haired bat", "myotis ciliolabrum", "western small footed bat", "corynorhinus townsendii", "Townsend's big-eared bat", "myotis volans", "long-legged myotis", "chrysemys picta", "painted turtle", "magellanic penguin", "penguin", "ruminant", "milking machine", "stegodyphus dumicola", "Nest Material", "nest component", "nest", "pig farm environment","Silkworms", "pantherophis guttatus", "centroscymnus coelolepis", "dogfish", "bathysaurus ferox", "deepsea lizard fish", "agalychnis callidryas", "craugastor fitzingeri", "ictalurus punctatus", "channel catfish", "apodemus uralensis", "ural field mouse", "mus spicilegus", "steppe mouse", "pontoporia blainvillei", "la plata dolphin", "hermetia illucens", "black soldier fly", "surface of larvae", "surface veterinary clinic", "veterinary clinic", "slaughter house", "meat", "pork", "poultry meat", "laboratory mouse", "animal cecal content", "a feces sample of chicken origin", "Black bone chicken", "Ayam Cemani", "boneless skinless chicken thigh", "breaded chicken breast nugget", "breast from whole chicken", "chicken bedding sample", "chicken broccoli and cheese", "chicken carcass rinse", "chicken carcasses showing gastrointestinal lesions", "chicken from supermarket", "chicken from wild market", "chicken frozen form", "chicken giblets", "chicken gizzards", "chicken gizzards and liver", "chicken gizzards and necks", "chicken gizzards and feet", "chicken gizzards and wings", "chicken gizzards and thighs", "chicken gizzards and drumsticks", "chicken gizzards and backs", "chicken gizzards and livers", "chicken gizzards and hearts", "chicken gizzards and kidneys", "chicken gizzards and lungs", "chicken gizzards and spleens", "chicken gizzards and pancreases", "chicken gizzards and intestines", "chicken house fly trap", "chicken house fly trap inside", "chicken house fly trap outside", "chicken house front inside", "chicken house front outside", "chicken house litter", "chicken house litter back", "chicken house litter front", "chicken house manure pile", "chicken house middle", "chicken house swab", "chicken liver", "chicken liver with fowl typhoid", "chicken neck skin", 'chicken slaughterhouse plucked skin', 'chicken slaughterhouse saw', 'chicken slaughterhouse tool', "chicken tender raw", "chicken tetrazzini", "chicken thigh and drumstick", "chicken thigh skin on", "chicken whole broiler carcass", "chicken whole cut", "chicken whole-cut in the lab", "cooked chicken thigh", "curry raw chicken skewer", "drag swab poultry chicken", "dried vegetable chicken soup", "emulsified chicken", "chicken fluff", "chicken fluffy", "feces chicken house back", "feces chicken house front", "feces chicken house middle", "finished chicken", "finished chicken fat", "food pre-packaged", "fresh whole chicken carcass pluck shop", "freshly slaughtered native chicken", "freshly slaughtered silkie chicken", "frozen chicken carcass", "frozen chicken carcass", "frozen chicken part supermarket", "frozen chicken shawrma", "frozen chicken tender", "frozen chicken thigh and leg", "frozen raw chicken grind pet food", "frozen raw stuffed chicken", "frozen raw stuffed chicken", "frozen raw stuffed chicken product", "frozen raw stuffed chicken product", "frozen whole chicken carcasses", "frozen whole chicken carcasses", "mayonnaise chicken risotto", "nonintact chicken", "non-organic retail chicken", "packaged chicken", "packaged chicken broiler leg", "peanut chicken soup", "popcorn chicken", "post chill chicken carcass", "pre-chill chicken carcass", "product raw intact chicken", "product chicken", "product raw intact chicken", "raw stuffed breaded chicken", "raw stuffed chicken products", "ready-to-cook chicken", "ready-to-cook chicken cutlets", "refrigerated raw chicken thigh and leg", "retail chicken carcasse", "retail chicken gizzard", "retail chicken gizzards", "retail chicken quarter leg", "retail chicken skin", "rice stew with chicken", "rinse from chicken carcass", "roasted chicken pea powder", "roasted chicken salad", "salted chicken breast from poultry slaughterhouse", "seasoned chicken fillet", "sick organs of a dead chicken", "stewed chicken", "stuffed chicken products", "uncooked chicken tikka", "wild chicken liver avian", "myulchi jeot", "anchovy", "anchovy jeot", "jeotgal", "jeot", "anchovy fish", "anchovy snack", "animal liver", "animal origin", "animal stool", "animal-snake", "animal-turtle", "Apple snail", "Pomacea canaliculata", "Pomacea bridgesii", "Golden Apple Snail", "Pomacea diffusa", "Pomacea paludosa", "Florida Apple Snail", "Ampullariidae", "Ampullariidae Family", "avian", "avian fluff", "Frog", "Toad", "Anura", "Rana", "Bufo", "Lithobates", "Hyla", "Pseudacris", "Pelophylax", "Acris", "Anaxyrus", "Ranitomeya", "Adenomera", "Aparasphenodon", "Aparasphenodon brunoi", "Aparasphenodon arapuca", "Aparasphenodon sp.", "Aparasphenodon sp. nov.", "bufo bufo", "common toad", "bufo viridis", "green toad", "bufo calamita", "natterjack toad", "bufo japonicus", "japanese toad", "bufo boreas", "western toad", "bufo cognatus", "great plains toad", "rana pipiens", "northern leopard frog", "rana catesbeiana", "bullfrog", "rana sylvatica", "wood frog", "rana clamitans", "green frog", "rana temporaria", "common frog", "rana aurora", "red-legged frog", "rana muscosa", "mountain yellow-legged frog", "rana draytonii", "bullfrog", "rana boylii", "foothill yellow-legged frog", "rana pipiens", "northern leopard frog", "rana catesbeiana", "bullfrog", "rana sylvatica", "wood frog", "rana clamitans", "green frog", "rana temporaria", "common frog", "rana aurora", "red-legged frog", "rana muscosa", "mountain yellow-legged frog", "rana draytonii", "bullfrog", "american bullfrog", "rana boylii", "foothill yellow-legged frog", "baby toad from irrigation pond", "baby toad", "baby frog", "frog tadpole", "frog tadpoles", "frog spawn", "frog eggs", "froglet", "froglets", "toad tadpole", "toad tadpoles", "toad spawn", "toad eggs", "toadlet", "toadlets", "beef carcass", "beef carcass trim", "beef for fajita", "beef trimmings", "boiled pork with mustard greens", "boneless beef", "boneless beef", "boneless pork", "boneless pork butt", "boneless pork shoulder", "bovine booties", "bovine feces", "bovine insect fly composite", "bovine post anaerobic digestion", "bovine stool", "broadhead fish", "broiler carcass", "cage surface", "calf feces", "canned fish", "canned meat", "canned seafood", "frozen seafood", "frozen meat", "frozen poultry", "frozen fish", "frozen shrimp", "frozen crab", "frozen lobster", "frozen scallops", "frozen clams", "frozen mussels", "frozen squid", "frozen octopus", "caprine lung", "carcass of slaughterhouse", "carcass rinse", "carcass swab", "carcass chiller", "calf dung", "calf excreta", "calf manure", "calf droppings", "calf stool", "bovine feces", "bovine dung", "bovine excreta", "bovine manure", "bovine droppings", "bovine stool", "bovine fecal matter", "bovine fecal sample", "bovine fecal swab", "bovine fecal content", "bovine fecal sample", "bovine fecal swab", "bovine fecal content", "bovine fecal matter", "bovine feces sample", "bovine feces swab", "bovine feces content", "bovine feces matter", "bovine feces sample", "bovine feces swab", "bovine feces content", "bovine feces matter", "bovine slurry", "bovine manure slurry", "bovine fecal slurry", "bovine fecal matter slurry", "cattle slurry", "cattle manure slurry", "cattle fecal slurry", "cattle fecal matter slurry", "cecal feces of duck slaughterhouse", "Celebes crested macaque", "Hylobatidae", "Gibbons", "Hominidae", "Homininae", "Hominini", "Pongo", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "celebes ape",],

        "Environment-associated" : ["Marine slduge", "lemon bar", "cinnamon halva bar", "brown meal", "terminalia chebula", "tamarindus indica", "barn air", "quercus mongolica", "buahinia variegata", "morella rubra", "terminalia ivorensis", "zea mays", "coal bed", "crutose lichen", "lichen", "gaillardia", "Blanket flowers", "saccharina japonica", "Kombu", "ulva", "Sea lettuce", "cos lettuce", "romaine lettuce", "Corals", "pocillopora grandis", "Antler corals", "plant epidermis", "pesticide contaminated farm soil","Pesticide contaminated soil", "Cream products","cream", "red oak lettuce","tango lettuce","shredded lettuce","red leaf lettuce","loose leaf lettuce","organic red romaine lettuce","white mushroom","mushroom","cultivated mushroom","Fried snack food","Fried snack pack","Fried snakc","Food","urtica dioica", "Soil sample", "Soil isolate", 'environmental', 'environment', 'grass containing manure','grass and manure','grass with manure','grass without manure', 'grass no manure','ground grass', 'background grass', 'lawn', 'grass', 'water', 'freshwater', 'tap water', 'surface water', 'stormwater', 'drinking water', 'mineral water', 'desalinated water', 'filtered water', 'sea water', 'saltwater', 'ocean water', 'marine water', 'coastal water', 'saltwater lagoon', 'brine water', 'groundwater', 'subsurface water', 'aquifer water', 'well water', 'spring water', 'artesian water', 'borehole water', 'wastewater', 'sewage', 'municipal wastewater', 'graywater', 'blackwater', 'septic tank water', 'contaminated water', 'recycled wastewater', 'waterwater treatment', 'wastewater treatment effluent', "wwtp", "chlorinated wwtp", "non clorinated wwtp", 'treated effluent', 'river water', 'stream water', 'creek water', 'brook water', 'riverside water', 'tributary water', 'river runoff', 'river', 'stream', 'tributary', 'lake water', 'lake sample', 'freshwater lake', 'saltwater lake', 'lake runoff', 'lake basin water', 'rainwater', 'precipitation', 'catchment rainwater', 'weathered rainwater', 'rainfall', 'downpour water', 'pond water', 'shallow water pond', 'wetland water', 'swamp pond water', 'still water pond', 'pond effluent', 'industrial effluent', 'chemical effluent', 'polluted industrial water', 'industrial discharge', 'chemical plant effluent', 'effluent', 'soil', 'generic soil', 'earth', 'soil particles', 'dirt', 'sand', 'loam', 'rocky soil', 'non-specific soil', 'soil debris', 'topsoil', 'subsoil', 'soil matter', 'earthy matter', 'natural soil', 'soil mixture', 'agricultural soil', 'farmland soil', 'farming soil', 'cultivated soil', 'crop soil', 'planting soil', 'agricultural land soil', 'crop-field soil', 'farm soil', 'tillable soil', 'field soil', 'farmland', 'soil for agriculture', 'irrigated soil', 'cultivated land soil', 'forest soil', 'woodland soil', 'rainforest soil', 'temperate forest soil', 'tropical forest soil', 'deciduous forest soil', 'coniferous forest soil', 'pine forest soil', 'wooded area soil', 'mountain forest soil', 'tropical rainforest soil', 'boreal forest soil', 'forest floor soil', 'old growth forest soil', 'forest soil type', 'natural forest soil', 'desert soil', 'arid soil', 'sandy soil', 'dry soil', 'dusty soil', 'desert sand', 'desert terrain soil', 'saline soil', 'semi-arid soil', 'desertified soil', 'cactus soil', 'shrubland soil', 'bare soil', 'dune soil', 'desert ecosystem soil', 'low organic soil','urban','rural','from city','from village', 'urban soil', 'city soil', 'town soil', 'urban environment soil', 'paved area soil', 'built-up soil', 'urban landscape soil', 'street soil', 'roadside soil', 'park soil', 'urban garden soil', 'city agricultural soil', 'residential area soil', 'commercial area soil', 'soil in urban areas', 'soil in suburbs', 'industrial soil', 'factory soil', 'polluted soil', 'contaminated soil', 'toxic soil', 'mining soil', 'industrial waste soil', 'manufacturing soil', 'industrial plant soil', 'chemical soil', 'soil in industrial zones', 'pollution affected soil', 'urban industrial soil', 'waste disposal area soil', 'mine tailings', 'excavated soil', 'processed soil', 'mineral soil', 'heavy metal soil', 'mining waste soil', 'mine waste soil', 'extracted soil', 'mining dump soil', 'waste rock soil', 'tailings pile', 'ore residue soil', 'residual mining soil', 'volcanic soil', 'lava soil', 'basaltic soil', 'tuff soil', 'ash soil', 'pumice soil', 'volcanic ash soil', 'volcanic eruption soil', 'volcanic rock soil', 'fertile volcanic soil', 'tephra soil', 'hawaiian volcanic soil', 'andesite soil', 'volcanic island soil', 'compost', 'organic compost', 'soil enrichment', 'soil amendment', 'decomposed organic matter', 'humus', 'compost soil', 'mulch', 'organic matter', 'plant-based compost', 'garden compost', 'organic fertilizer', 'green waste compost', 'kitchen waste compost', 'composted material', 'composted soil', 'soil conditioner', 'clay soil', 'heavy soil', 'sticky soil', 'loamy clay soil', 'clay-rich soil', 'clay loam soil', 'red clay soil', 'yellow clay soil', 'black clay soil', 'clay texture soil', 'plastic soil', 'clayey soil', 'dense clay soil', 'moist clay soil', 'clay-based soil', 'stiff soil', 'permafrost', 'frozen soil', 'cold soil', 'tundra soil', 'frozen earth', 'glacial soil', 'arctic soil', 'polar soil', 'deep frozen soil', 'subsurface ice soil', 'permanent ice soil', 'arctic tundra soil', 'permafrost layers', 'ice-rich soil', 'permafrost soil formation', 'peat soil', 'moorland soil', 'sphagnum soil', 'bog soil', 'swamp soil', 'peat moss soil', 'organic soil', 'hydric peat soil', 'peaty soil', 'high organic matter soil', 'wetland peat soil', 'tropical peat soil', 'humic soil', 'turf soil', 'acidic peat soil', 'sedimentary soil', 'sandstone soil', 'limestone soil', 'shale soil', 'claystone soil', 'fossiliferous soil', 'gravelly soil', 'riverbed soil', 'alluvial soil', 'soil from sedimentary layers', 'depositional soil', 'soil from river deposits', 'sediment-based soil', 'deltaic soil', 'sedimentary basin soil', 'coastal soil', 'beach soil', 'sandy beach soil', 'saltwater soil', 'marine soil', 'estuarine soil', 'shoreline soil', 'tidal soil', 'salt marsh soil', 'mangrove soil', 'coastal wetland soil', 'seashore soil', 'rocky coast soil', 'seaside soil', 'coastal ecosystem soil', 'salty soil', 'beach sand', 'coastal sand', 'seaside sand', 'marine sand', 'tidal sand', 'sandy beach', 'shore sand', 'sea sand', 'coastal sediment', 'white sand', 'golden sand', 'coral sand', 'sand dune', 'beach area sand', 'landfill soil', 'waste soil', 'garbage soil', 'landfill waste', 'polluted landfill soil', 'landfill waste soil', 'municipal waste soil', 'waste disposal soil', 'trash heap soil', 'dump soil', 'hazardous waste soil', 'recycled landfill soil', 'garbage pit soil', 'landfill contaminated soil', 'soil in landfill sites', 'decomposed waste soil', 'rhizosphere soil', 'root nodule soil', 'rhizobium soil', 'leguminous soil', 'nodular soil', 'nitrogen fixing soil', 'rhizobial soil', 'soil with nitrogen-fixing bacteria', 'rhizobium inoculated soil', 'root symbiosis soil', 'root zone soil', 'legume-rhizobium soil', 'air', 'indoor air', 'home air', 'office air', 'building air', 'indoor air quality', 'room air', 'living room air', 'bedroom air', 'indoor environment air', 'indoor air samples', 'indoor air composition', 'air in homes', 'indoor climate', 'indoor air quality assessment', 'heated indoor air', 'outdoor air', 'city air', 'outdoor air quality', 'urban air', 'country air', 'park air', 'outdoor environment air', 'open environment air', 'outdoor air samples', 'hvac system', 'heating, ventilation, and air conditioning', 'hvac air', 'air system', 'ducted air system', 'ventilation system', 'central air', 'hvac system air', 'building air system', 'air circulation system', 'hvac duct air', 'air filtered system', 'commercial hvac system', 'residential hvac system', 'ventilation ducts', 'ducted air', 'air ducts', 'airflow ducts', 'hvac ducts', 'ductwork', 'ducted ventilation', 'air duct system', 'hvac duct system', 'ventilation duct air', 'ductwork contamination', 'duct air samples', 'residential duct system', 'commercial duct system', 'aerosols', 'bioaerosols', 'biological aerosols', 'microbial aerosols', 'fungal aerosols', 'bacterial aerosols', 'viruses in air', 'spores in air', 'microorganisms in air', 'bio-contaminants in air', 'indoor bioaerosols', 'outdoor bioaerosols', 'airborne pathogens', 'airborne microorganisms', 'factory emissions', 'industrial emissions', 'industrial pollution', 'factory air', 'manufacturing emissions', 'emissions from factories', 'factory exhaust', 'industrial air pollutants', 'airborne factory waste', 'pollutant particles from factories', 'chemical emissions', 'factory gaseous emissions', 'smoke particles', 'smoke', 'smoke in air', 'particulate smoke', 'pollution from smoke', 'fire smoke', 'wildfire smoke', 'smoke exposure', 'airborne smoke particles', 'soot in air', 'carcinogenic smoke', 'aerosolized smoke', 'burning smoke', 'pollen', 'pollen particles', 'flower pollen', 'tree pollen', 'grass pollen', 'pollen grains', 'airborne pollen', 'pollination particles', 'pollen in air', 'pollen exposure', 'pollen in atmosphere', 'allergen pollen', 'seasonal pollen', 'pollen clouds', 'greenhouse air', 'controlled greenhouse air', 'greenhouse air quality', 'horticultural air', 'indoor greenhouse air', 'greenhouse environment', 'greenhouse ventilation air', 'hydroponic air', 'agricultural greenhouse air', 'greenhouse gas concentration', 'atmosphere in greenhouses', 'high humidity greenhouse air', 'lab air', 'laboratory air', 'research lab air', 'laboratory environment', 'cleanroom air', 'air in labs', 'chemical lab air', 'biological lab air', 'controlled lab air', 'indoor laboratory air', 'air quality in labs', 'sterile lab air', 'research environment air', 'pollution-free lab air', 'smog', 'pollution smog', 'industrial smog', 'smoggy air', 'fossil fuel emissions', 'particulate smog', 'smog from cars', 'urban smog', 'heavy smog', 'fog and smog', 'smog in cities', 'haze and smog', 'chemical smog', 'photochemical smog', 'airborne pollutants', 'pollutants in air', 'air pollution', 'toxic air particles', 'airborne toxic materials', 'environmental pollutants', 'pollution particles', 'pollutants from industry', 'car emissions in air', 'airborne chemicals', 'hazardous pollutants', 'particulate matter in air', 'airborne contaminants', 'pollution aerosols', 'chemical contaminants in air', 'dust particles', 'dust', 'dust in air', 'dust samples', 'particulate matter', 'pm2.5', 'pm10', 'indoor dust', 'outdoor dust', 'fine dust', 'coarse dust', 'respirable dust', 'dust accumulation', 'dust mites', 'dust storm air', 'root', 'root sample', 'edible root', 'taproot', 'fibrous root', 'tuber', 'starchy root', 'underground stem', 'bulb', 'rhizome', 'root system', 'root cutting', 'carrot', 'potato', 'sweet potato', 'beetroot', 'radish', 'turnip', 'ginger', 'garlic', 'onion', 'yam', 'cassava', 'parsnip', 'daikon', 'cassava root', 'garlic bulb', 'radish root', 'turmeric', 'daikon root', 'sugar beet', 'jerusalem artichoke', 'yam root', 'salsify root', 'chicory root', 'artichoke root', 'burdock root', 'dandelion root', 'sunchoke', 'lotus root', 'sweet yam', 'beet', 'cucumber root', 'mallow root', 'arrowroot', 'black salsify', 'chufa tuber', 'water chestnut', 'burdock', 'konjac root', 'water yam', 'kohlrabi root', 'salsify', 'crosne', 'black garlic', 'jerusalem artichoke tuber', 'water lily root', 'leaf', 'leaf sample', 'green leaf', 'plant leaf', 'tender leaf', 'edible leaf', 'vegetable leaf', 'spinach', 'lettuce', 'kale', 'cabbage', 'broccoli', 'cauliflower', 'mint leaf', 'herb leaf', 'tropical leaf', 'herbaceous leaf', 'leafy greens', 'basil', 'parsley', 'oregano', 'chard', 'mustard greens', 'collard greens', 'dandelion greens', 'moringa leaf', 'brussels sprouts', 'bean leaves', 'tomato leaf', 'alfalfa', 'watercress', 'artichoke leaf', 'radicchio', 'romaine lettuce', 'endive', 'escarole', 'mache', 'sorrel', 'arugula', 'chamomile', 'lemon balm', 'thyme leaf', 'bitter melon leaf', 'ginger leaf', 'paprika leaf', 'sage', 'oregano leaf', 'curry leaf', 'chili leaf', 'apple leaf', 'pear leaf', 'peach leaf', 'almond leaf', 'lemon leaf', 'citrus leaf', 'green tea leaf', 'aloe vera leaf', 'spearmint', 'peppermint', 'bay leaf', 'coriander leaf', 'lime leaf', 'grape leaf', 'pineapple leaf', 'sorghum leaf', 'sunflower leaf', 'tobacco leaf', 'cassava leaf', 'passionfruit leaf', 'bamboo leaf', 'palm leaf', 'mango leaf', 'chili pepper leaf', 'seed', 'seed sample', 'plant seed', 'tree seed', 'flower seed', 'corn seed', 'bean seed', 'wheat seed', 'rice seed', 'soybean seed', 'cotton seed', 'pea seed', 'sunflower seed', 'pumpkin seed', 'melon seed', 'chili seed', 'tomato seed', 'cucumber seed', 'apple seed', 'grape seed', 'pepper seed', 'seedling', 'seed pod', 'radish seed', 'alfalfa seed', 'mustard seed', 'carrot seed', 'sorghum seed', 'barley seed', 'canola seed', 'watermelon seed', 'papaya seed', 'avocado seed', 'tomato seeds', 'bell pepper seed', 'squash seed', 'cabbage seed', 'broccoli seed', 'aubergine seed', 'sweet corn seed', 'green bean seed', 'beetroot seed', 'cantaloupe seed', 'kiwifruit seed', 'fennel seed', 'coriander seed', 'dandelion seed', 'asparagus seed', 'okra seed', 'wheatgrass seed', 'taro seed', 'carambola seed', 'pomegranate seed', 'mango seed', 'passionfruit seed', 'lemongrass seed', 'apple seedling', 'pear seed', 'banana seed', 'peach seed', 'pineapple seed', 'sweet pepper seed', 'citrus seed', 'apricot seed', 'plum seed', 'olive seed', 'cherry seed', 'stem', 'stem sample', 'plant stem', 'tree stem', 'woody stem', 'herbaceous stem', 'flower stem', 'shrub stem', 'cactus stem', 'corn stalk', 'sugarcane stem', 'cucumber vine', 'tomato stem', 'potato stem', 'pumpkin vine', 'watermelon vine', 'grapevine', 'cotton stem', 'sorghum stem', 'coffee plant stem', 'cotton plant stem', 'chili pepper plant stem', 'tobacco plant stem', 'banana plant stem', 'kale stalk', 'brussels sprouts stalk', 'cabbage stalk', 'carrot stalk', 'ginseng root stalk', 'yucca stem', 'dragon fruit stem', 'cabbage stem', 'okra stem', 'chili stem', 'sunflower stalk', 'sorghum cane', 'bamboo stalk', 'zucchini stalk', 'eggplant stem', 'peanut plant stem', 'artichoke stem', 'artichoke heart stalk', 'sugar beet stem', 'jute stem', 'mango tree stem', 'coffee bean stem', 'pineapple stem', 'palm stem', 'soybean stem', 'banana stalk', 'avocado tree stem', 'pomegranate tree stem', 'cassava stem', 'cassowary stem', 'hemp stem', 'lemon tree stem', 'palm leaf stem', 'walnut stem', 'fig tree stem', 'fruit', 'fruit sample', 'edible fruit', 'apple', 'banana', 'grape', 'orange', 'pear', 'cherry', 'melon', 'strawberry', 'peach', 'nectarine', 'pineapple', 'watermelon', 'kiwi', 'mango', 'blueberry', 'avocado', 'papaya', 'lemon', 'lime', 'plum', 'grapefruit', 'pomegranate', 'apricot', 'fig', 'coconut', 'lychee', 'passionfruit', 'olive', 'blackberry', 'gooseberry', 'dragon fruit', 'starfruit', 'tangerine', 'persimmon', 'mulberry', 'black currant', 'elderberry', 'soursop', 'jackfruit', 'mandarin', 'durian', 'longan', 'rambutan', 'custard apple', 'cantaloupe', 'kumquat', 'pawpaw', 'tamarind', 'barbados cherry', 'cherry plum', 'green apple', 'saskatoon berry', 'red currant', 'cranberry', 'blackcurrant', 'goji berry', 'sea buckthorn', 'chayote', 'pluot', 'clementine', 'lingonberry', 'juneberry', 'longan fruit', 'white mulberry', 'yellow watermelon', 'fuzzy melon', 'quince', 'marula', 'bitter orange', 'mango chutney fruit', 'flowering plant', 'flowering species', 'blooming plant', 'flower-bearing plant', 'annual flowering plant', 'perennial flowering plant', 'sunflower', 'tulip', 'rose', 'orchid', 'lily', 'daisy', 'carnation', 'wildflower', 'lilac', 'hibiscus', 'geranium', 'chrysanthemum', 'lavender', 'violet', 'begonia', 'calendula', 'marigold', 'zinnia', 'morning glory', 'bougainvillea', 'jasmine', 'petunia', 'poppy', 'plumeria', 'freesia', 'pansy', 'lotus flower', 'camellia', 'daffodil', 'azalea', 'orchidaceae', 'magnolia', 'lotus', 'water lily', 'fuchsia', 'heliotrope', 'crocus', 'wild orchid', 'peony', 'dandelion', 'iris', 'angelonia', 'snapdragon', 'bluebell', 'cherry blossom', 'hibiscus flower', 'saffron crocus', 'snapdragon flower', 'petunia flower', 'morning glory vine', 'magnolia flower', 'fuchsia plant', 'camellia plant', 'honeysuckle', 'red hot poker', 'hydrangea', 'alstroemeria', 'gladiolus', 'periwinkle', 'maranta', 'helichrysum', 'mock orange', 'indigofera', 'algae', 'green algae', 'blue-green algae', 'red algae', 'brown algae', 'marine algae', 'algal bloom', 'algal sample', 'pond algae', 'freshwater algae', 'seaweed', 'kelp', 'sargassum', 'spirulina', 'chlorella', 'diatom', 'phyto plankton', 'euglena', 'porphyra', 'gracilaria', 'nori', 'agar', 'fucus', 'ulva', 'carrageenan', 'zostera', 'ascophyllum', 'fucus vesiculosus', 'sargassum muticum', 'ulva lactuca', 'spirogyra', 'chara', 'volvox', 'red tide algae', 'tetraselmis', 'euglena gracilis', 'coccolithophores', 'chlorophyta', 'sargassum fusiforme', 'sea lettuce', 'kelp forest', 'edible seaweed', 'algal mat', 'sea moss', 'rockweed', 'bladderwrack', 'nori algae', 'wakame', 'kombu', 'dulse', 'irish moss', 'agar-agar', 'green algal bloom', 'blue-green phytoplankton', 'bark', 'tree bark', 'woody bark', 'tree trunk bark', 'plant bark', 'cork bark', 'outer bark', 'bark sample', 'birch bark', 'oak bark', 'pine bark', 'maple bark', 'cedar bark', 'willow bark', 'eucalyptus bark', 'cinnamon bark', 'cherry bark', 'bamboo bark', 'rubber tree bark', 'mango bark', 'teak bark', 'mahogany bark', 'alder bark', 'redwood bark', 'poplar bark', 'chestnut bark', 'beech bark', 'birchwood', 'carob bark', 'juniper bark', 'linden bark', 'larch bark', 'sequoia bark', 'buckthorn bark', 'sweet birch bark', 'black cherry bark', 'bamboo sheath', 'hickory bark', 'sandalwood bark', 'camphor bark', 'neem bark', 'sassafras bark', 'lemon bark', 'palo santo bark', 'mesquite bark', 'prickly ash bark', 'sandpaper tree', 'waste', 'garbage', 'trash', 'rubbish', 'refuse', 'scrap', 'discarded materials', 'general waste', 'unsorted waste', 'unclassified waste', 'uncategorized refuse', 'municipal solid waste', 'household waste', 'residential waste', 'domestic waste', 'city trash', 'urban waste', 'city garbage', 'neighborhood waste', 'street trash', 'public waste', 'community waste', 'municipal refuse', 'city refuse', 'town waste', 'household plastic waste', 'plastic packaging', 'plastic bottles', 'plastic bags', 'industrial waste', 'factory waste', 'manufacturing waste', 'chemical waste', 'polluted industrial waste', 'toxic industrial waste', 'mining waste', 'industrial byproducts', 'production waste', 'waste from factories', 'metal scrap waste', 'chemical byproducts', 'industrial residues', 'plastic waste from factories', 'plastic manufacturing waste', 'textile waste', 'fabric waste', 'clothing waste', 'discarded textiles', 'worn-out clothes', 'clothing scraps', 'textile byproducts', 'old clothing', 'damaged textiles', 'fashion waste', 'used fabric', 'textile disposal', 'scrapped textiles', 'waste fabric', 'second-hand clothing waste', 'hazardous waste', 'toxic waste', 'dangerous waste', 'flammable waste', 'radioactive waste', 'poisonous waste', 'contaminated waste', 'corrosive waste', 'biohazard waste', 'explosive waste', 'pollutant waste', 'hazardous material waste', 'dangerous goods waste', 'dangerous chemicals waste', 'food waste', 'organic waste', 'kitchen waste', 'leftover food', 'spoiled food', 'food scraps', 'compostable waste', 'discarded food', 'uneaten food', 'food debris', 'waste from food processing', 'food leftovers', 'perishable waste', 'unwanted food', 'rotting food', 'construction waste', 'demolition waste', 'building debris', 'construction debris', 'demolition rubble', 'building materials waste', 'concrete waste', 'metal scrap', 'renovation waste', 'site waste', 'construction leftovers', 'demolition materials', 'building rubble', 'construction site waste', 'medical waste', 'hospital waste', 'biomedical waste', 'clinical waste', 'infectious waste', 'sharps', 'used syringes', 'pharmaceutical waste', 'expired medications', 'discarded medical instruments', 'toxic healthcare waste', 'patient waste', 'used bandages', 'used needles', 'surgical waste', 'agricultural waste', 'farm waste', 'crop residue', 'animal manure', 'plant debris', 'harvest waste', 'crop byproducts', 'agricultural runoff', 'pesticide containers', 'fertilizer packaging', 'agriculture-related waste', 'agricultural chemicals waste', 'agro-waste', 'farm animal waste', 'waste from farming', 'wastewater sludge', 'sewage sludge', 'sewage waste', 'treated sludge', 'wastewater treatment byproducts', 'sludge from sewage plants', 'effluent sludge', 'sewage treatment sludge', 'municipal sludge', 'organic sludge', 'industrial sludge', 'sludge from water treatment', 'biosolids', 'wastewater solids', 'treatment plant waste', 'wood waste', 'wooden debris', 'scrap wood', 'wood scraps', 'wooden pallets', 'old timber', 'discarded wood products', 'sawdust', 'wood chips', 'wooden packaging waste', 'scrapped wood', 'unused timber', 'waste wood planks', 'wooden material waste', 'epilobium', 'manure field', 'field manure', 'manure on field', 'freshwater lake sample', 'plant', 'tankmilk', 'dairymilk', 'stored milk','soybean hull', 'cotton seed hull', 'hull pellets','blood meal','feed pellet','raw feed','chicken feed', 'farm feed', 'animal feed', 'poultry feed', 'livestock feed', 'cattle feed', 'pig feed', 'horse feed', 'goat feed', 'sheep feed', 'rabbit feed', 'fish feed', 'aquatic feed', 'pet feed', 'dog feed', 'cat feed', 'bird feed', 'wildlife feed', 'zoo feed', 'aquarium feed', 'insect feed', 'rodent feed',"Pistachios", "Almonds", "Cashews", "Peanuts", "Walnuts", "Hazelnuts", "Macadamia Nuts", "Pecans", "Brazil Nuts", "Pine Nuts", "Chestnuts", "Pistachio Nuts", "Almond Nuts", "Cashew Nuts", "Peanut Nuts", "Walnut Nuts", "Hazelnut Nuts", "Macadamia Nuts", "Pecan Nuts", "peppercorn", "cinnamon", "cardamom", "coriander", "basil", "ocimum basilicum", "cilantro",'fishmeal','poultry feather meal', 'bonemeal', "Galia muskmelon", "muskmelon", "Galia melon", "melon", "Sow manure", "milk powder", "dried milk powder","macadamia nut", "ready to eat food","pre-made food","readymade food","Pasta", 'meringue', "breading mixture", "breading mix", "Sprout", "moong sprout", "moong bean sprout", "pulse sprout", "butterscotch pudding", "pudding", "chocolate pudding", "Meal", "Dinner", "Breakfast", "yam ring", "taro basket", "taro", "yam", "anise spice", "star anise", "whey protein concentrate", "Pate", "Bread", "dried noodle", "romaine", "kale", "scarlet kale", "chard", "red chard", "green chard", "neem powder", "dessert", "marine sediment", "Cyanobacterial mats", "bolivina costata", "environmental", "coconut meat", 'tamazula river', 'spring roll', 'hot underground spring', 'hot water spring', 'spring', 'mountain spring', 'lichen', 'cladonia cristatella', 'artisanal cheese', 'organic buckwheat noodle', "organic buckwheat flour", "organic buckwheat", "organic buckwheat grain", "organic buckwheat groats", "organic buckwheat flakes", "organic buckwheat cereal", "organic buckwheat pasta", "organic buckwheat noodles", "organic buckwheat bread", "organic buckwheat pancake mix", "organic buckwheat granola", "organic buckwheat muesli", "organic buckwheat porridge", "organic buckwheat flour mix", "organic buckwheat flour blend", "vegetarian instant noodle", "vegetarian noodle", "vegetarian ramen", "vegetarian udon", "vegetarian soba", "vegetarian pasta", "vegetarian spaghetti", "vegetarian macaroni", "vegetarian fusilli", "vegetarian penne", "vegetarian rotini", "vegetarian farfalle", "vegetarian orecchiette", "vegetarian gnocchi", "vegetarian ravioli", "vegetarian tortellini", "vegetarian lasagna", "vegetarian fettuccine", "vegetarian linguine", "sugar cane", "frozen sugar cane", "multigrain hakka noodle", "multigrain noodle", "multigrain ramen", "multigrain udon", "multigrain soba", "multigrain pasta", "multigrain spaghetti", "multigrain macaroni", "multigrain fusilli", "multigrain penne", "multigrain rotini", "multigrain farfalle", "multigrain orecchiette", "multigrain gnocchi", "multigrain ravioli", "multigrain tortellini", "multigrain lasagna", "multigrain fettuccine", "multigrain linguine", "nor-mai-dong; sour bamboo", "sour bamboo", "pumiao rice noodle", "aegiceras corniculatum", "mangrove", "mangrove plant", "mangrove roots", "marine sponge", "cheese core stage", "cheese rind stage", "desert plant free sand dune", "desert plant free sand dune", "desert rock", "vachellia gummifera", "vachellia tortilis subsp. raddiana", "rhizosphere of crop plant from saline desert area", "waste lagoon stage", "elephant ginger", "air from international space station", "ared foot platform of international space station", "atlantic intertidal shore", "campsite water hose", "chocolate sandwich cookies with peanut cream fill", "chocolate sandwich", "cooked peanuts from wild market", "diseased rice seeds and leaf sheaths", "international space station", "international space station surface", "pmm port 1 of international space station", "rhizosphere soil of crown gall infected cherry rootstock-colt", "rhizosphere soil", "crown gall", "crown gall infection", "rice", "rice seed", "root knot nematode egg mass", "meloidogyne sp.", "tomato", "rice", "phaseolus vulgaris", "black rice plant", "skim milk clt", "surface of international space station", "water sample of lake-2 sach pass", "western north pacific station s1", "zoo composting manure", "populus x jackii", "cat food raw chicken", "daily vitamin for cat and kitten", "finished chubs of kitten grind", "kitten grind", "pharmaceutical production line", "raw cat food chicken", "raw cat food duck", "raw cat food tuna", "raw cat food", "fish meal", "fish meal feed", "fish meal for animal feed", "fish masala", "freshwater ditch", "freshwater from boreal shield lake", "freshwater sediment", "freshwater stream sediment", "freshwater stream system, activated sludge", "freshwater: river", "sediment of a fish dumping yard at balugaon near chilika lake", "sediment of a fish dumping yard", "spice mixture for lahori fish", "spices and herbs", "spices and salt-ground black pepper", "water area where the ocean and a freshwater spring meet", "gold and copper mine", "industry eqipments", "commercial building tap", "chloraminated drinking water", "commercial tattoo ink", "cuatro cienegas basin", "frozen animal feed adult mouse", "frozen animal feed baby mouse", "frozen animal feed feeder mouse", "frozen animal feed young mouse", "frozen feeder mouse", "frozen fuzzie feeder mouse", "frozen pet food hopper mouse", "frozen pinkie feeder mouse", "heterotrigona itama honey", "Honey", "phytophthora alni", "kitchen ware", "kitchen countertop", "kitchen sink taps", "kitchen unused sink", "kitchen unused sink drain", "kitchen sink", "kitchen cutting board", "kitchen utensils", "kitchen appliances", "kitchen equipment", "kitchen tools", "kitchen surfaces", "kitchen area", "kitchen environment", "kitchen space", "kitchen setting", "kitchen worktop", "kitchen work surface", "kitchen work area", "kitchen work environment", "kitchen work setting", "kitchen work space", "Morus", "Mulberry", "Mulberry fruit", "morus alba", "marine hydrothermal chimney wall", "marine hydrothermal sulfide sediment", "marine hydrothermal sulfide sediment", "orinoco river basin", "pintassilgo field from potiguar basin", "raw gallnut honey", "residence cold water tap", "sediment associated with methane hydrate from krishna godavari basin", "sufu", "sufu sample", "fermented sufu", "fermented tofu", "surface of lake", "lake", "surface seawater sample", "surface snow", "surface soil", "residential aged care facility", "ventilago sp.", "water stream", "sponge", "ground hot chili", "hot chili", "ground chili", "a fermented food zha-chili", "fermented food", "medicago", "amylase concentrate", "Industrial product", "artificial lake water", "bakery environment", "burger", "Chili", "coastal surface sea water", "coldroom for food storage", "complementary food of infant", "Compost", "dehydrate small chili", "medicinal plant", "endophytic isolate", "endophyte", "physalis ixocarpa", "endosphere", "Feminine hygiene products", "Sanitary pads", "Tampons", "Menstrual cups", "Whisper", "Stayfree","Make-up removers", "Cleansing oil", "Micellar water", "Perfumes", "Eau de parfum", "Body spray", "Fragrance mist", "Bath & Body Works","Deodorants", "Colognes", "Aftershave cologne", "Body wash", "Shower gel", "Shower cream", "Soaps", "Bar soap", "Toilet soap","Bath salts", "Soaking salt","Body scrubs", "Exfoliating scrub", "Sugar scrub","Loofahs", "Bath sponge", "Body puff", "Toothbrushes", "Brush", "Electric toothbrush", "Toothpaste", "Tooth gel", "Fluoride paste", "Mouthwash", "Oral rinse", "Dental floss", "Flossing thread", "Floss picks","Teeth whitening kits", "Whitening strips", "Charcoal toothpaste", "Shampoos", "Hair shampoo", "Head wash", "Conditioners", "Hair conditioner", "Silkening cream","Hair oils", "Coconut oil", "Argan oil", "Deep conditioning mask", "Hair styling gels", "Hair wax", "Face wash", "Facial cleanser", "Moisturizers", "Face cream", "Lotion", "Sunscreen", "Sunblock", "SPF lotion", "Exfoliators", "Face scrub", "Scrubbing gel", "Facial masks", "Sheet masks", "Face packs", "Foundation sets", "Base makeup", "Lipstick combos", "Lip colors", "Eye shadow palettes", "Eyeshadow", "Blush & contour kits", "Cheek palette", "Blusher", "Contour powder", "Make-up brushes", "Brush sets", "Blending brushes", "Make-up sponges", "Beauty blender", "Makeup puff", "Hair brushes", "Hairbrush", "Combs", "Hair comb", "Nail clippers", "Nail cutter", "Manicure tools", "Pedicure kit", "Tweezers", "Pluckers", "Shaving kits", "Shave set", "Shaving set", "Razors", "Shavers", "Electric razors", "Gillette", "Wilkinson Sword", "Beard trimmers", "Hair trimmers", "Glacial ice", "glacial ice water", "glacial water", "glacial surface ice", "glacial", "ground medium hot molido chili", "ground red chili", "hot spring environment", "Coal bed", "coal seam", "coal deposit", "carboniferous bed", "peat layer", "bituminous coal zone", "Deep sea", "abyss", "oceanic trench", "hadal zone", "mariana trench", "deep ocean", "benthic zone", "Hydrothermal vent", "black smoker", "white smoker", "seafloor vent", "deep-sea vent", "geothermal vent", "marine hydrothermal chimney wall", "marine hydrothermal sulfide sediment", "marine hydrothermal sulfide sediment", "Hot spring", "thermal spring", "geyser", "boiling spring", "fumarole", "geothermal spring", "hot spring environment", "Acid mine drainage", "AMD", "acidic runoff", "mine water", "sulfuric drainage", "contaminated mine water", "Desert", "arid zone", "drought region", "sand dunes", "semi-arid land", "barren land", "Hypersaline lake", "salt lake", "brine pool", "saline lagoon", "soda lake", "alkaline lake", "Antarctic", "polar region", "south pole", "ice shelf", "glacial zone", "permafrost area", "Deep subsurface", "deep biosphere", "geothermal subsurface", "rock pore system", "subterranean habitat", "High radiation zone", "radioactive site", "chernobyl exclusion zone", "gamma-ray exposed area", "nuclear waste site", "Permafrost", "frozen soil", "tundra", "ice-bound ground", "periglacial zone", "cryotic soil", "Alkaline soil", "soda soil", "high pH soil", "basic soil", "calcareous soil", "salt-affected soil", "Volcanic crater", "lava crater", "caldera", "active volcano", "magma chamber", "pyroclastic zone", "Black shale", "organic-rich shale", "kerogen rock", "carbonaceous shale", "petroleum source rock", "Petroleum reservoir", "oil field", "crude oil deposit", "hydrocarbon reservoir", "tar sand", "bitumen field", "oil sands", "tar sands", "Soda lake", "alkaline lake", "high pH water body", "sodium carbonate lake", "brine lake", "Extreme high altitude", "mountain peak", "oxygen-deprived zone", "Himalayan plateau", "andes summit", "Deep cave", "karst cave", "limestone cave", "subterranean cavern", "underground chamber", "Household", "Home", "Residential Area", "Domestic Environment", "Living Space", "Private Residence", "Family Home", "Apartment", "Condominium", "Townhouse", "Single-Family Home", "Duplex", "Bungalow", "Cottage", "Villa", "Mansion", "Ranch House", "Mobile Home", "Tiny House", "Care Facility", "Nursing Home", "Assisted Living Facility", "Long-Term Care Facility", "Senior Living Community", "Retirement Home", "Skilled Nursing Facility", "Memory Care Facility", "Rehabilitation Center", "Group Home", "Residential Care Facility", "Personal Care Home", "Elderly Care Facility", "Adult Foster Care", "Continuing Care Retirement Community", "Senior Care Center", "Alzheimer's Care Facility", "Palliative Care Facility", "Hospice Care Facility", "Respite Care Facility", "Transitional Care Facility", "Subacute Care Facility", "Independent Living Facility", "Senior Housing Community", "Assisted Living Residence", "Nursing Care Facility", "Long-Term Care Home", "Elderly Residential Facility", "Senior Assisted Living Facility", "Skilled Nursing Home", "Memory Support Facility", "Rehabilitation Hospital", "Group Living Facility", "Residential Rehabilitation Facility", "Personal Care Residence", "Elderly Group Home", "Adult Day Care Facility", "Continuing Care Facility", "Senior Living Residence", "Alzheimer's Care Home", "Palliative Care Home", "Hospice Residence", "Respite Care Home", "Transitional Care Home", "Subacute Care Home", "Independent Living Residence", "Senior Housing Facility", "Assisted Living Community", "residential aged care facility", "nereocystis luetkeana", "subtidal kelp", "kimchi system", "rapeseed", "winter rape", "Medicago lupulina", "maize stem", "marine sample from the pacific ocean", "marine sponge", "alectra sessiliflora", "medicinal plant", "pothos", "heliconia", "anthurium", "mud from rhine river", "nut raisin blend trail mixture", "oil field", "oil refineries", "rhaphiolepis umbellata", "pyrus", "prunus avium", "broussonetia kazinoki", "castanea crenata", "prunus yedoensis", "daphniphyllum teijsmannii", "aesculus hippocastanum", "dendropanax trifidus", "eriobotrya japonica", "morella rubra", "fraxinus excelsior", "nerium oleander", "medicago", "medicago falcata", "solanum tuberosum", "plant nodule", "plant root nodules", "raisin", "Residential Areas", "rice field", "rice paddy", "rind of surface ripened cheese", "ripe olives", "root surface", "root tuber", "glycine max", "wenyujin", "sesame field", "soil inside the cave", "soft rock", "soil surface", "soil subsurface", "standing rain water", "protoceratium reticulatum", "dinoflagellate", "polyalae radix", "surface mud", "surface sediment", "surface soil of a manganese mine", "surface sterilized root", "oryza sativa", "callitris preissii", "native pine tree", "jasmin rice", "water for irrigation", "aegiceras corniculatum", "black mangrove", "river mangrove", "jatropha curcas", "Sedum sp." , "ripariosida hermaphrodita", "tannery", "tannery bath containers", "tannery lime containers", "water pipe", "water spray", "water statue", "water tank", "zha-chili", "maple", "medicago plant", "commercial cultured foods", "commercial dietary supplements", "cultured food", "humus", "hummus", "humus sample", "hummus sample", "goat cheese", "long-term care facility", "river", "river water", "traditional korean fermented vegetable", "rain sample", "human sludge", "Soil Crust", "Biological Soil Crust", "Cryptobiotic Soil Crust", "Soil Surface Crust", "Soil Microbial Crust", "Soil Biological Layer", "Soil Surface Layer", "Soil Microbial Layer", "Soil Surface Biofilm", "Soil Surface Microbiome", "Soil Surface Community", "Incubated in Space", "Space Incubation", "Microgravity Incubation", "Spacecraft Incubation", "Extraterrestrial Incubation", "Zero Gravity Incubation", "Space Environment Incubation", "Space Habitat Incubation", "Spacecraft Culture", "Space", "extraterrestrial culture", "zero gravity culture", "outer space culture", "mars surface culture", "ISS culture", "degreasser", "kitchen floor", "kitchen wall", "kitchen ceiling", "kitchen cabinet", "kitchen drawer", "kitchen shelf", "kitchen countertop surface", "kitchen sink surface", "kitchen cutting board surface", "kitchen utensil surface", "kitchen appliance surface", "kitchen equipment surface", "kitchen tool surface", "kitchen area surface", "kitchen environment surface", "kitchen space surface", "kitchen setting surface", "kitchen worktop surface", "kitchen work area surface", "kitchen work environment surface", "kitchen work setting surface", "kitchen work space surface", "degresser sample", "urban sediment", "urban soil", "urban soil sample", "urban soil core", "urban soil profile", "urban soil horizon", "urban soil layer", "urban soil column", "urban soil block", "urban soil aggregate", "urban soil clod", "urban soil crumb", "urban soil particle", "urban soil grain", "urban soil fragment", "urban soil chunk", "urban soil piece", "sediment samples", "soil from a culture of crop", "chicken from pet food", "chicken jerky dog treat", "chicken jerky pet treat", "chicken quail blend for dog food", "cocktail snacks chicken flavor cashew shaped biscuits", "curry spice mixture chicken ginger masala", "raw dog food chicken organic vegetable", "turkey chicken formula cat kitten food", "raw dehydrated chicken supreme dog food", "pet food with turkey chicken whitefish", "pet treat chicken jerky", "fat milk", "all purpose flour", "all purpose wheat flour", "spring roll", "Fried snack food", "Fried snack pack", "Fried snakc", "snack", "snack pack", "snack size", "snack platter", "snack tray", "snack box", "snack bag", "snack pouch", "snack cup", "snack bowl", "snack plate", "snack stick", "corn chips snack", "chip coconut chips", "chip", "chips", "cocktail snacks chicken flavor cashew shaped biscuits", "nachos", "trail mix", "snacks", "snack food", "snack mix", "snack bar", "almond flatbread snack badam lachha", "almond powder", "unused razor head", "unused razor", "unused razor blade", "unused razor blades", "unused razor head sample", "unused razor head isolate", "sewage sludge", "sludge from wastewater treatment plant", "sludge from wastewater treatment", "sludge from sewage treatment plant", "sludge from sewage treatment", "anaerobic activated sludge of molasses wastewater in a continuous stirred-tank reactor", "anaerobic activated sludge of molasses wastewater", "anaerobic activated sludge", "anaerobic sludge", "anaerobic sludge of molasses wastewater in a continuous stirred-tank reactor", "anaerobic sludge of molasses wastewater", "anaerobic sludge of molasses wastewater", "anaerobic sludge of molasses wastewater in a continuous stirred-tank reactor", "anaerobic sludge of molasses wastewater", "anaheim pepper", "jalapeno pepper", "habanero pepper", "serrano pepper", "banana pepper", "bell pepper", "poblano pepper", "cayenne pepper", "chipotle pepper", "green chili", "red chili", "yellow chili", "orange chili", "purple chili", "white chili", "black chili", "brown chili", "chili pod", "chili seed", "chili plant", "chili bush", "chili tree", "chili flower", "apple cider", "cider", "cider vinegar", "apple cider vinegar", "apple cider sauce", "apple cider dressing", "apple cider marinade", "apple cider dip", "apple cider glaze", "apple cider syrup", "apple cider jelly", "apple cider jam", "apple cider chutney", "apple cider relish", "apple cider salsa", "apple cider sauce mix", "apple cider spice mix", "apple cider spice blend", "apple cider spice rub", "apple cider spice seasoning", "apple cider spice paste", "apple cider spice marinade", "arthicoke", "artichoke", "artichoke heart", "artichoke leaf", "artichoke stem", "artichoke bud", "artichoke flower", "artichoke root", "artichoke tuber", "artichoke stalk", "artichoke petal", "artichoke bract", "artichoke scale", "artichoke leaf scale", "artichoke leaf bract", "artichoke leaf petal", "artichoke leaf bud", "artisanal dairy product", "asiago cheese sauce", "asiago cheese dip", "asiago cheese dressing", "asiago cheese marinade", "asiago cheese spread", "asiago cheese paste", "asiago cheese mix", "asiago cheese blend", "asiago cheese rub", "asiago cheese seasoning", "asiago cheese spice mix", "asiago cheese spice blend", "asiago cheese spice rub", "asiago cheese spice seasoning", "asiago cheese spice paste", "asiago cheese spice marinade", "baking flour", "baking wheat flour", "bread flour", "bread wheat flour", "cake flour", "cake wheat flour", "pastry flour", "pastry wheat flour", "self rising flour", "self rising wheat flour", "whole wheat flour", "whole grain wheat flour", "whole grain flour", "whole grain all purpose flour", "whole grain bread flour", "whole grain cake flour", "whole grain pastry flour", "whole grain self rising flour", "barbeque seasoning", "barbeque seasoning", "barbeque spice mix", "barbeque rub", "bay water spring sample", "beef strip pet treat", "beetroot sprout", "broccoli sprout", "radish sprout", "sunflower sprout", "pea sprout", "lentil sprout", "chickpea sprout", "mung bean sprout", "soybean sprout", "alfalfa sprout", "clover sprout", "mustard sprout", "cabbage sprout", "cauliflower sprout", "kale sprout", "brussels sprouts", "broccoli rabe sprout", "broccolini sprout", "rapini sprout", "arugula sprout", "watercress sprout", "radish seedling", "bhel spicy snack mixture", "biogas plant", "anaerobic digestion plant", "anaerobic digestion facility", "biogas facility", "biogas production plant", "anaerobic digestion system", "biogas generation plant", "biogas production facility", "anaerobic digestion unit", "biogas unit", "anaerobic digestion reactor", "biogas reactor", "anaerobic digester", "biogas digester", "anaerobic sludge of molasses wastewater in a continuous stirred-tank reactor", "Anaerobic Digestion", "Biogas Plant", "Biogas Reactor", "Biogas Fermentation", "Biogas Generation", "Biogas System", "Biogas Facility", "Anaerobic Digestion Sludge", "Anaerobic Digestion Residue", "Anaerobic Digestion Byproducts", "Anaerobic Digestion Waste", "Biogas Treatment Plant", "Biogas Production Facility", "Anaerobic Digestion Process", "Anaerobic Digestion Technology", "Biogas-related", "Biogas", "biomethane", "green gas", "biomethane gas", "biogas fuel", "Biogas Byproducts", "Biogas Byproduct Residue", "Biogas Byproduct Waste", "Biogas Byproduct Slurry", "Biogas Byproduct Effluent", "Biogas Effluent", "Biogas Slurry", "Biogas Residue", "Biogas Waste", "Biogas Digestate", "Biogas Residuals", "Biogas Byproduct Treatment", "Biogas Byproduct Management", "Biogas Byproduct Utilization", "Biogas Byproduct Recycling", "Biogas Byproduct Disposal", "Biogas Byproduct Processing", "Biogas Byproduct Handling", "Biogas Byproduct Storage", "Biogas Byproduct Separation", "Biogas Byproduct Extraction", "Biogas Byproduct Analysis", "Biogas Byproduct Characterization", "Biogas Digestate", "peanut kernel", "nut kernel", "almond kernel", "cashew kernel", "walnut kernel", "pecan kernel", "hazelnut kernel", "macadamia kernel", "pistachio kernel", "brazil nut kernel", "pine nut kernel", "chestnut kernel", "nut butter", "almond butter", "peanut butter", "cashew butter", "walnut butter", "pecan butter", "hazelnut butter", "macadamia butter", "pistachio butter", "brazil nut butter", "pine nut butter", "chestnut butter", "blended whey", "blueberry fermentation", "Boiled rice", "cooked rice", "Cooked Grains", "cacao-beans", "caciotta cheese", "grass", "grassland", "grass sample", "grass isolate", "grass mucilage", "grass exudate", "grass exudates", "grass exudate sample", "grass exudate isolate", "grass material", "grasses", "grasses of the family Poaceae", "grasses of the family Cyperaceae", "canadian high arctic grass", "arctic grass", "arctic grass sample", "arctic grass isolate", "arctic grass mucilage", "arctic grass exudate", "arctic grass exudates", "arctic grass exudate sample", "arctic grass exudate isolate", "arctic grass material", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Cyperaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "cactus", "cactus sample", "cactus isolate", "cactus mucilage", "cactus exudate", "cactus exudates", "cactus exudate sample", "cactus exudate isolate", "cactus material", "cacti", "canned food", "canned vegetables", "canned fruits", "canned beans", "canned soup", "canned tomatoes", "canned corn", "canned peas", "canned carrots", "canned mushrooms", "canned olives", "canned pickles", "canned fruit cocktail", "canned chili", "canned pasta sauce", "canned onions", "canned garlic", "canned artichokes", "canned asparagus", "canned beets", "canned pumpkin", "canned sweet potatoes", "canned green beans", "canned mixed vegetables", "canned fruit salad", "canned fruit cocktail syrup", "canned fruit juice", "canned fruit puree", "canned fruit compote", "canned fruit preserves", "canned fruit spread", "frozen food", "frozen vegetables", "frozen fruits", "frozen meals", "frozen dinners", "frozen pizza", "frozen desserts", "frozen snacks", "frozen appetizers", "frozen breakfast foods", "frozen entrees", "frozen side dishes", "canola meal feed bulk", "canola meal feed", "canola meal animal feed", "canola meal for animal feed", "canola meal for livestock feed", "canola meal for poultry feed", "canola meal for fish feed", "canola meal for pet food", "canola meal for dog food", "canola meal for cat food", "canola meal for wildlife feed", "canola meal for zoo feed", "canola meal for aquarium feed", "canola meal for insect feed", "canola meal for rodent feed", "canola meal for animal food", "canola meal for livestock food", "canola meal for poultry food", "canola meal for fish food", "canola meal for pet food", "capsicum powder", "caribe pepper", "cayenne pepper", "chili powder", "cashew piece", "walnut piece", "pecan piece", "hazelnut piece", "macadamia piece", "pistachio piece", "brazil nut piece", "pine nut piece", "chestnut piece", "almond piece", "peanut piece", "nut piece", "nut blend", "nut mix", "nut assortment", "nut medley", "nut combination", "nut variety", "nut selection", "nut collection", "nut platter", "nut tray", "nut bowl", "nut basket", "nut jar", "nut bag", "nut pouch", "raw compost", "asterionella formosa", "Skeletonema costatum", "thalassiosira pseudonana", "thalassiosira oceanica", "thalassiosira weissflogii", "thalassiosira rotula", "thalassiosira gravida", "thalassiosira nordenskioeldii", "thalassiosira punctigera", "thalassiosira decipiens", "thalassiosira delicatula", "thalassiosira symmetrica", "thalassiosira subtilis", "thalassiosira partheneia", "thalassiosira lacustris", "thalassiosira minima", "Haptophytes", "Phaeocystis globosa", "Phaeocystis pouchetii", "Phaeocystis antarctica", "Phaeocystis globosa", "Phaeocystis pouchetii", "Phaeocystis antarctica", "Chrysophytes", "Synura", "Dinobryon", "Ochromonas", "Mallomonas", "Chrysosphaerella", "Chrysococcus", "Chrysodidymus", "Chrysonebula", "Chrysotilos", "Chrysocapsa", "Chrysocystis", "Chrysocystis fragilis", "Chrysocystis minor", "Chrysocystis major", "Chrysocystis parva", "pavlova", "pavlova lutheri", "pavlova salina", "pavlova sp.", "pavlova sp. lutheri", "pavlova sp. salina", "pavlova sp. lutheri", "pavlova sp. salina", "pavlova sp. lutheri", "pavlova sp. salina", "Parks", "Public Parks", "City Parks", "Urban Parks", "Community Parks", "Recreational Parks", "National Parks", "State Parks", "Botanical Gardens", "Playgrounds", "Green Spaces", "Parkland", "Public Transportation", "Buses", "Trains", "Subways", "Light Rail", "Streetcars", "Ferries", "Public Transit Systems", "Mass Transit", "Public Transit", "Public Transportation Systems", "Public Transit Networks", "Public Transportation Services", "Public Transit Routes", "Public Transportation Infrastructure", "Public Transit Vehicles", "Public Transportation Facilities", "Public Transit Stations", "Public Transportation Hubs", "Public Transit Stops", "Public Transportation Schedules", "Public Transit Maps", "Public Transportation Accessibility", "Shopping Centers", "Malls", "Retail Complexes", "Commercial Centers", "Shopping Malls", "Retail Centers", "Shopping Plazas", "Commercial Complexes", "Shopping Districts", "Retail Parks", "Shopping Arcades", "Commercial Districts", "Shopping Streets", "Retail Streets", "Shopping Avenues", "Commercial Streets", "Shopping Boulevards", "Retail Boulevards", "Streets", "Roads", "Avenues", "Boulevards", "Lanes", "Highways", "Public Roads", "City Streets", "Urban Roads", "Residential Streets", "Commercial Streets", "Rural Roads", "Suburban Streets", "Street Corners", "Street Intersections", "Street Crossings", "Street Medians", "Street Sidewalks", "Public Buildings", "Government Buildings", "Libraries", "Museums", "Community Centers", "City Halls", "Courthouses", "Public Schools", "Public Hospitals", "Public Libraries", "Public Museums", "Public Community Centers", "Public Government Buildings", "Public City Halls", "Public Courthouses", "Public Schools", "Public Hospitals", "eucalyptus", "ipomoea pes-caprae", "pigeon pea", "Marine", "chickpea", ],
        
        "Laboratory-based" : ["lab strain", "lab mutant", "laboratory strain","culture", "pure culture","lab derived", "lab constructed", "lab synthesized", "synthetic sample", "cultured sample", "experimental sample", "isolated sample", "reconstructed sample", "engineered sample", "laboratory sample", "laboratory isolate", "in vitro sample", "lab-grown", "synthetically produced", "lab-based isolate", "bioengineered sample", "cultivated isolate", "artificial sample", "extracted sample", "processed sample", "isolated culture", "cell-cultured sample", "lab-generated", "in vivo-derived sample", "constructed sample", "experimental isolate", "recombinant sample", "biological isolate", "isolated culture specimen", "synthetic isolate", "chemically synthesized sample", "genetically engineered sample", "recombinant isolate", "broth", "lb media", "lysogeny media", "agar plate", "solid culture", "liquid culture", "lysogeny broth", "LB media", "human cell culture", "lab collection", "atcc", "mtcc", "laboratory derivative", "laboratory mutant", "laboratory", "laboratory mutant", "laboratory passage", "laboratory strain", "lab-derived", "plate contaminant", "lab contaminant", "proficiency panel isolate", "wild type strain", "wildtype obtained from arisolatebank from cdc", "culture:lab stock", "culture:lab strain", "culture:lab isolate", "culture:lab mutant", "culture:lab derived", "culture:lab constructed", "culture:lab synthesized", "culture:synthetic sample", "culture:cultured sample", "culture:experimental sample", "culture:isolated sample", "culture:reconstructed sample", "culture:engineered sample", "culture:laboratory sample", "culture:laboratory isolate", "culture:in vitro sample", "culture:lab-grown", "culture:synthetically produced", "culture:lab-based isolate", "culture:bioengineered sample", "culture:cultivated isolate", "culture:artificial sample", "culture:extracted sample", "culture:processed sample", "culture:isolated culture", "culture:cell-cultured sample", "culture:lab-generated", "culture:in vivo-derived sample", "culture:constructed sample", "culture:experimental isolate", "culture:recombinant sample", "culture:biological isolate", "culture:isolated culture specimen", "culture:synthetic isolate", "culture:chemically synthesized sample", "culture:genetically engineered sample", "culture:recombinant isolate"],
    }

    best_match_score = 0
    best_host_match = "Unknown"

    try:
        if isinstance(fuzzy_str, str):
            logging.debug(f"Entered the str loop in find_host function with string: `{fuzzy_str}`.")
            normalized_str = fuzzy_str.strip("(").strip(")").strip().lower()

            for category, terms in host_categories.items():
                for term in terms:
                    similarity_score = fuzz.ratio(normalized_str, term.lower())  
                    if similarity_score > best_match_score:
                        best_match_score = similarity_score
                        best_host_match = category
                        best_matched_term = term
            
            if best_match_score >= 70:
                logging.debug(f"Host identified for string: `{fuzzy_str}`. The identified host is `{best_host_match}` with a score of {best_match_score}. The best matched term was `{best_matched_term}`.")
                return best_host_match, best_match_score
            
            logging.debug(f"The non-splitted form of text didn't gave the match. Trying with the splitted form of text.")

            part_digit_free = re.sub(r'\d+', '', fuzzy_str)
            part = re.sub(r'[\[\]{}()]', '', part_digit_free)
            part = re.sub(r'"', '', part)
            part = re.sub(r'\'', '', part)
            parts = re.split(r'[-+:;/_,\s]|\band\b', part)
            logging.debug(f"Splitting the string: `{fuzzy_str}` into parts: {parts}")

            searches = []
            for part_i in range(0, len(parts), 2):
                searches.append(" ".join(parts[part_i:part_i+2]))


            for search_word in searches:
                normalized_str = search_word.strip().lower()
                for category, terms in host_categories.items():
                    for term in terms:
                        similarity_score = fuzz.ratio(normalized_str, term.lower())  
                        if similarity_score > best_match_score:
                            best_match_score = similarity_score
                            best_host_match = category
                            best_matched_term = term

            if best_match_score >= 80:
                logging.debug(f"Host identified for string: `{fuzzy_str}`. The identified host is `{best_host_match}` with a score of {best_match_score}. The best matched term was `{best_matched_term}`.")
                return best_host_match, best_match_score
            
            else:
                logging.debug(f"The string: `{fuzzy_str}` didn't match any entry from the database. Highest match was with `{best_host_match}` with {best_match_score} score. Returning `Unknown`.")
                return "Unknown", 0

        else:
            logging.debug(f"Given string: `{fuzzy_str}` is not a string, returning `Unknown`.")
            return "Unknown", 0
    
    except Exception as e:
        logging.error(f"An unexpected error occurred while finding host for string `{fuzzy_str}`: {e}")
        return "Unknown", 0

def find_human_disease(fuzzy_str):
    logging.debug(f"Started! finding the human disease for string `{fuzzy_str}`")
    host_disease_terms = {
        "Unknown" : {
            'Unknown' : ["Clinical material"]
            },
        "Asymptomatic/Healthy" : {
            'Asymptomatic' : ["Asymptomatic", "non-symtomatic", "Healthy", "no-disease", "carrier", "healthy", "healthy carrier", "healthy individual", "healthy person", "healthy subject", "healthy volunteer", "non-symptomatic individual", "non-symptomatic person", "non-symptomatic subject", "non-symptomatic volunteer", "asymptomatic carrier", "asymptomatic individual", "asymptomatic person", "asymptomatic subject", "asymptomatic volunteer", "asymptomatic carriage", "asymptomatic state", "asymptomatic condition", "asymptomatic status", "asymptomatic infection", "asymptomatic disease", "asymptomatic presentation", "asymptomatic phase", "asymptomatic period", "asymptomatic stage", "asymptomatic form", "asymptomatic manifestation", "asymptomatic response", "asymptomatic carrier state", "asymptomatic carrier condition", "healthy volounteers", "healthy subject", "healthy person", "healthy adult", "healthy male", "healthy female", "healthy human", "asymptomatic infection", "asymptomatic healthy carrier"],
            },

        "Bacterial Infections (BI)": {
            "UTIs": ["uti", "urinary tract", "urinary tract infection", "urinary infection", "urinary tract infections", "urinary System Infection", "complicated urinary tract infection", "cystitis", "chronic cystitis", "urethritis", "acute cystitis", "human uti", "community-acquired urinary", "bacturia", "right ureter", "left ureter", "ureter", "anterior urethritis", "recurrent uti", "complicated UTI", "uncomplicated UTI", "uncomplicated urinary tract infection", "uncomplicated urinary tract infections", "uncomplicated cystitis",],
            "Bacteremia": ["infectious shock", "septicopyemia", "blood isolate", "bloodstream Infection", "bacteremia", "blood stream infection", "septicemia",  "BSI", "meningococcemia", "bacteremic pneumonia", "fatal septicemia"],
            "Sepsis/Septic Shock" : ["sepsis", "septic shock", "sepsis and the multiorgan failure", "sepsis, septic shock", "severe sepsis", "septicemia", "severe sepsis with septic shock", "sepsis-related organ failure", "puerperal sepsis", "meningococcemia", "bacterial sepsis", "fatal sepsis", "puerperal sepsis", "bacterial sepsis", "fatal sepsis", "sepsis patient", "urosepsis", "uro-sepsis", "uro sepsis", "neonatal sepsis", "septic infection"],
            "Respiratory Infections/Illnesses": ["respiratory","pulmonary contusion", "pulmonary", "pulmonary sample", "pulmonary Biopsy", "pleural disease", "tuberculosis", "Mycobacterium", "bronchus", "Tuberculosis patient", "TB", "TB patient", "patient with tuberculosis", "patient with TB", "Mycobacterium tuberculosis", "Nasopharyngeal infection","Nasal infection","nasal passage","bronchial infection","pneumonia", "hospital-acquired pneumonia", "respiratory infection", "pulmonary tuberculosis", "Pulmonary TB", "PTB", "Extra pulmonary TB", "Extra pulmonary tuberculosis", "BCGitis", "respiratory tract infections", "bronchopneumonia", "chronic respiratory failure", "pneumonia, sepsis", "invasive klebsiella pneumonia syndrome", "severe pneumonia", "polysegmental community acquired pneumonia", "severe acute pneumonia", "Respiratory tract colonization", "Pneumocystis pneumonia", "Bacteremic Pneumonia", 'obstructive lung disease', 'Chronic obstructive pulmonary disease', 'copd', 'expectoration', 'bronchiolitis', "upper respiratory tract infection", "upper respiratory infection", "sinusitis", "pharyngitis", "tonsillitis", "laryngitis", "Stuffed nose", "nasopharyngitis", "rhinitis", "acute rhinosinusitis", "acute sinusitis", "acute pharyngitis", "acute tonsillitis", "acute laryngitis", "acute nasopharyngitis", "type i respiratory failure", "type i respiratory failure", "aspiratory pneumonia", "respiratory", "acute respiratory failure", "brochial infection", "pulmonary infection", "Community-acquired pneumonia", "cap", "bilateral pneumonia", "disseminated tb", "Tuberculosis", "acute respiratory insufficiency", "chronic respiratory insufficiency", "cystic fibrosis", "non-cystic fibrosis bronchiectasis", "bronchiectasis", "bronchiolitis obliterans", "bronchiolitis obliterans organizing pneumonia", "bronchiolitis", "bronchitis", "bronchopneumonia", "cepacia syndrome", "cepacia", "empyema", "penumonia with empyema", "respiratory distress syndrome", "delayed hypersensitivity pneumonitis", "pulmonary emboli", "tuberculous pleurisy", "chest infection", "CF", "cap patient", "community-acquired pneumonia", "community-acquired pneumonia patient", "community-acquired pneumonia with bacteremia", "community-acquired pneumonia with sepsis", "community-acquired pneumonia with septic shock", "community-acquired pneumonia with respiratory failure", "community-acquired pneumonia with pleural effusion", "community-acquired pneumonia with empyema", "community-acquired pneumonia with lung abscess", "community-acquired pneumonia with parapneumonic effusion", "community-acquired pneumonia with necrotizing pneumonia", "copd patient", "inhalation syndrome", "idiopathic pulmonary fibrosis", "pneumoconiosis", "dyspnea", "copd with acute exacerbation", "pneumothorax", "hemopneumothorax", "upper right atelectasis", "atelenctasis", "emphysema", "empyema thoracis", "emphysema mediastinum", "emptysis", "silicosis", "lung infection", "abscess peritonsillar", "abscess retropharyngeal", "peritonsillar abscess", "retropharyngeal abscess", "lung abscess", "lung infection", "pulmonary abscess", "pulmonary infection", "lung disease", "pulmonary disease", "nasal abscess", "chronic purulent sinusitis", "primary bacterial intercostal pyomyositis", "intercoastal pyomyisitis", "bronchiectasis", "adenoiditis", "purulent chronic rhinosinusitis", "non-cystic fibrosis bronchiectasis", "cystic fibrosis infant bronchiectasis", "cystic fibrosis adult severe", "diffuse panbronchiolitis", "hemoptysis", "chronic rhinosinusitis", "pleuritis", "patient with bronchiectasis", "lobar pneumonia", "pulmonary edema", "pleuropneumoniae", "respiratory airways", "respiratory culture", "respiratory isolate", "respiratory processed", "respiratory specimen", "respiratory tract", "respiratory tract of human", "respiratory tract sample", "pulmonary nodules", "CRSsNP", "CRSwNP", "asphyxia neonatalis gravis", "URT", "subglottis stenosis", "resp culture"],
            "Soft Tissue Infections/Colonization": ["wound", "skin and soft tissue", "sacrum abscess", "decubitus sacrum", "wound infection", "soft tissue infection", "wound swab", "skin infection", "wound infection", "necrotizing fasciitis",  "pressure ulcer", "cellulitis", "chronic wound infection", "wound fluid", "wound drain fluid", "subphrenic abscess", "carbuncle", "ulcer", "punctate", "incisions", "tracheostoma", "pressure ulcer", "decubitus", "Sacrum ulcer", "sacral ulcer", "pressure sore", "soft tissue", "hip wound", "leg wound", "feet wound", "foot wound", "cutaneous abscess", "necrotizing infection", "non-healing wound", "infected ulcer", "cutaneous lesion", "skin abscess", "dermal infection", "fasciitis", "skin boil", "superficial wound infection", "deep wound infection", "infected pressure ulcer", "venous ulcer", "traumatic wound",  "bacterial skin infection", "open wound", "infected dermal tissue", "infected skin graft", "infected bite wound", "skin breakdown", "infected amputation site", "purulent wound", "ulcerative lesion", "gangrenous wound", "left hand wound", "right hand wound", "hand wound", "soft tissue", "Pitted keratolysis", "Erysipelas", "spongiotic dermatitis", "Nail bed", "superficial abscess", "dermal abscess", "breast abscess", "rectal abscess", "renal abscess", "cervical abscess", "groin abscess", "inguinal abscess", "perianal abscess", "perineal abscess", "abdominal abscess", "pelvic abscess", "subcutaneous abscess", "subcutaneous tissue infection", "subcutaneous tissue abscess", "hand abscess", "foot abscess", "leg abscess", "arm abscess", "thigh abscess", "buttock abscess", "toe infecttion", "finger infection", "hand infection", "leg infection", "arm infection", "thigh infection", "buttock infection", "toe abscess", "finger abscess", "ulcer infection", "hemorrhagic ulcera", "acne vulgaris", "skin disease", "dermatitis", "infected nail well", "paronychia", "impetigo", "scarlet fever", "scarlatiniform rash", "scarlatina", "actinomycosis", "acne pustules", "Necrotizing Soft Tissue Infection", "NSTI", "paronychia", "drain inserting site", "drain insertion site infection", "drain site infection", "drain site abscess", "drain site", "toxic epidermal necrolysis", "face swab", "facial acne", "facial cutaneous", "facial skin", "facial/lip", "finger", "hand", "leg", "arm", "thigh", "buttock", "pododermatitis", "retro auricular cellulitis", "skin lesion", "furoncle", "furunculosis", "gangrenous lesion", "gangrenous ulcer", "gangrenous wound", "gangrene", "gangrene wound", "gangrene ulcer", "gangrene lesion", "subcutaneous infection", "folliculitis", "Pemphigus vulgaris", "popliteal_fossa", "popliteal fossa", "purpuric lesion", "retroauricular crease", "right shin blister", "flank", "acne", "eczema herpeticum", "folliculitis", "skin and soft tissue infection", "gluteal skin infection", "pyoderma", "pyoderma gangrenosum", "pyoderma infections", "toe web infection", "toe web", "toe web space infection", "ulcus cruris"],
            "CNS infections": ["meningoencephalitis", "encephalitis", "bacterial meningitis", "meningitis", "brain abscess", "cerebral abscess", "brain infection", "melioidosis encephalomyelitis", "intracranial infection", "meningococcal meningitis", "streptococcal meningitis", "neuroinfection", "intracranial infection", "meningococcal disease", "central nervous system infection", "cns infection", "cerebellum abscess", "epidural abscess", "epidural Biopsy", "meningococcal disease", "invasive meningococcal disease", "meningococcal", "meningococcal infection", "invasive meningococcal infection", "invasive meningococcal disease", "post-traumatic meningitis"],
            "Bone Infection/Inflammation" : ['Bone isolate','orthopedic infections', 'Osteitis', 'osteomyelitis', 'Joint', 'Hip joint', 'pelvic joints', 'metatarsals', 'sarcum', "osteitis", "emphysematous osteomyelitis", "Septic arthritis", "joint infection", "cynovial fluid infection", "arthritis", "mastoid", "mastoiditis", "spinal osteomyelitis",  "knee infection", "medulary cavity", "sacroiliac abscess", "mandibular osteomyelitis", "osteoarticular", "articular cavity infection", "destructive osteomyelitis", "septic shoulder arthrtitis", "right ischial Biopsy", "right ischial tissue", "right hip sonicate", "septic arthritis", "spondylodiscitis", "glenoid membrane", "glenoid cavity", "glenoid fossa", "glenoid labrum", "glenoid cavity infection", "glenoid cavity disease", "glenoid cavity disorder", "glenoid cavity abscess", "glenoid cavity inflammation", "humeral membrane", "humeral cavity", "humeral fossa", "humeral labrum", "humeral cavity infection", "humeral cavity disease", "humeral cavity disorder", "humeral cavity abscess", "humeral cavity inflammation", "shoulder collar membrane", "shoulder collar cavity", "shoulder collar fossa", "shoulder collar labrum", "shoulder collar cavity infection", "shoulder collar cavity disease", "shoulder collar cavity disorder", "shoulder collar cavity abscess", "shoulder collar cavity inflammation", "shoulder joint membrane", "shoulder joint cavity", "shoulder joint fossa", "shoulder joint labrum", "shoulder joint cavity infection", "shoulder joint cavity disease", "shoulder joint cavity disorder", "shoulder joint cavity abscess", "shoulder joint cavity inflammation", "mastoiditis"],
            "Unsorted Bacterial infections": ["infection", "bacterial infection", "halscysta", "other isolate", "other infection", "throglossal duct cyst", "diseased", "infectious disease", "infections diseases", "bacterial infectious disease", "aneurysm with ventriculostomy", "pneumococcal infection", "fungal-bacterial co-infection", "acinetobacter infections", "mrsa infection", "clostridium difficile", "Colonization", "Bacterial colonization", "Microbial colonization", "Opportunistic infections", 'salmonellosis','Typhoid', 'typhoid fever', 'typhoid septicemia', 'typhoid infection', 'typhoid enteric fever', 'typhoid fever with perforation', 'typhoid fever with complications', 'typhoid fever with septicemia', 'enteric fever', 'connective tissue infection', 'CRAB infection', 'MDR infection', 'XDR infection', 'MDR', 'CRAB', 'XDR', 'AMR', 'Antibiotic resistant infection', 'Antimicrobial resistant infection', 'chorioamnionitis', 'Periprostheric', 'Periprostheric inflammation', 'periprosthetic biopsy', 'periprosthetic infection', 'prosthetic infection', 'Hydrothorax', 'Melioidosis', 'spondylodiscitis', 'mediastinitis', "streptococcal infectious disease", "streptococcal infection", "invasive neonatal disease", "methicillin-resistant staphylococcus aureus carriage", "drain infection", "streptococcus suis infection", "methicillin resistant staph aureus infection", "MRSA infection", "endometritisis", "escherichia coli infections", "ecoli infection", "escherichia coli disease", "ecoli disease", "escherichia coli infectious disease", "ecoli infectious disease", "gram-negative bacterial infections", "recurrent salmonellosis", "staphylococcal infection", "periodontal disease", "gas infection", "group a streptococcus infection", "Menstrual toxic shock syndrome", "sterile fluid infection", "infected cyst", "infected pilonidal cyst", "dissaminated infection", "pseudomonas aeruginosa infection", "pseudomonas infection", "klebsiella pneumoniae infection", "klebsiella infection", "staphylococcal infection", "methicillin-resistant staphylococcus aureus colonization", "infection with burkholderia thailandensis", "invasive pneumococcal disease", "systemic salmonella infection", "anthrax", "anthrax disease", "glanders", "Rheumatic fever",  "otomastoidite", "parotid infection", "parotitis", "parotiditis", "parotid gland infection", "parotid gland disease", "parotid duct infection", "parotid duct disease", "parotid duct disorder", "invasive infection", "invasive disease", "invasive infections", "salmonella foodborne diseases", "Leptospirosis", "Weil disease", "cervico-facial actinomycosis", "chorioamnionitis"],
            "Ear Infection" : ['otitis media', 'otitis externa', 'otitis', 'ear infection', 'ear disorder', 'ear disease', 'acute otitis', 'otorrhea', 'tympanitis', 'eardrum infection', 'ear canal infection', 'relapse otitis', 'otomastoidite'],
            "Oral Infection/Colonization": ['gingivitis', 'cavity', 'carious lesion', 'cheek abscess', 'endodontic infection', 'root canal infection', 'enamel caries', 'oral infection', 'oral disease', 'oral disorder', 'oral cavity infection', 'oral cavity disease', 'oral cavity disorder', 'oral mucosa infection', 'oral mucosa disease', 'oral mucosa disorder', 'gingival infection', 'gingival disease', 'gingival disorder', 'periodontal disease', 'periodontal infection', 'periodontal disease', 'periodontitis', 'dental abscess', 'Refractory endodontic lesions', 'supragingival dental plaque', 'supragingival plaque', 'subgingival plaque biofilm', 'sub-gingival plaque oral cavity', 'supragingival plaque biofilm of a molar'],
            "Eye Infection" : ["eye infection", "conjuctivitis", "eye inflammation", "eye disorder", "eye disease", "entophthalmia", "endophthalmitis", "Conjunctival Swab", "Conjunctival Sample","Conjunctival Discharge", "Cornea", "Corneal Sample", "Corneal Discharge", "Corneal Secretion", "Corneal Fluid", "Corneal Culture", "Eye isolate", "eye abscess", "eye secretion", "eye culture", "eye aspirate", "eye drainage", "eye swab isolate", "eye swab sample", "eye aspirate sample", "eye aspirate fluid", "eye aspirate specimen", "eye aspirate isolate", "eye aspirate culture", "keratitis", "microbial keratitis", "injected sclera", "contact lens related keratitis", "microbial keratitis", "conjunctiva infection", "ocular hordeolum", "stye", "stye infectioon", "scleritis", "scleral infection", "ocular infections", "ophtalmia"],
            "Ventilator-associated Pneumonia": [ "ventilator-associated pneumonia", "ventilator-related pneumonia", "ventilator-associated pneumoniae", "probable vap", "VAP"],
            "Nosocomial Infections": ["hospital-acquired infection",  "nosocomial pneumonia", "infections, nosocomial", "nosocomial infection", "healthcare-associated infection", "infections, nosocomial", "healthcare-associated infection", "nosocomial", "hospital acquired infection", "hospital-acquired infection after orthopedic intervention", "hospital-acquired infection after surgery", "hospital acquired pneumonia", "nosocomial infection", "nocomomial pathogen"],
            "STDs" : ['sexually trasmitted diseases', 'stds', 'gonorrhea', 'syphilis', 'chlamydia', 'genital herpes', 'pubic lice', 'salpingitis, stss', "asymptomatic gonorrhea", "asymptomatic chlamydia", "asymptomatic syphilis", "asymptomatic genital herpes", "asymptomatic trichomoniasis", "asymptomatic sexually transmitted infection", "neisseria gonorrhoeae positive", "std eval", "sexually transmitted infection", "neisseria infection", "uncomplicated gonorrhoea", "VDG", "venerial disease gonorrhea", "venerial disease", ],
            "ETEC Infection": ["enterotoxigenic e. coli", "enterotoxigenic escherichia coli", "enterotoxigenic e. coli infection", "enterotoxigenic escherichia coli infection", "etec infection", "enterotoxigenic e. coli disease", "enterotoxigenic escherichia coli disease", "etec infection"],
            "ExPEC Infection": ["extraintestinal pathogenic e. coli", "extraintestinal pathogenic escherichia coli", "expec infection", "extraintestinal pathogenic e. coli infection", "extraintestinal pathogenic escherichia coli infection", "expec disease", "extraintestinal pathogenic e. coli disease", "extraintestinal pathogenic escherichia coli disease", "extraintestinal e coli", "extraintestinal e. coli", ],
            "STEC Infection" : ["shiga toxin-producing e. coli", "shiga toxin-producing escherichia coli", "stec infection", "shiga toxin-producing e. coli infection", "shiga toxin-producing escherichia coli infection", "stec disease", "shiga toxin-producing e. coli disease", "shiga toxin-producing escherichia coli disease", "enterohemorrhagic e. coli", "enterohemorrhagic escherichia coli", "verotoxin-producing e.coli", "enterohemorrhagic escherichia coli", "stec_foodborne_diseases",  "stec diarrhoeal disease"],
            "EAEC Infection" : ["enteroaggregative e. coli", "enteroaggregative escherichia coli", "eaec infection", "enteroaggregative e. coli infection", "enteroaggregative escherichia coli infection", "eaec disease", "enteroaggregative e. coli disease", "enteroaggregative escherichia coli disease"],
            "Mastitis" : ["mastitis", "mastitis infection",],
            "UPEC Infection" : ["uropathogenic e. coli", "uropathogenic escherichia coli", "upec infection", "uropathogenic e. coli infection", "uropathogenic escherichia coli infection", "upec disease", "uropathogenic e. coli disease", "uropathogenic escherichia coli disease"],
            "Catheter-associated Infections": ["catheter-related infection", "catheter-associated infection", "catheter-associated urinary tract infection", "catheter-associated bloodstream infection", "catheter-associated pneumonia", "catheter-associated sepsis", "catheter-related bloodstream infection", "catheter-related urinary tract infection", "catheter-related pneumonia", "catheter-related sepsis", "central line associated bloodstream infection", "catheter-related bloodstream infection", "central venous line infection", "central line-associated bacterial bloodstream infection", "central line infection", "bacteremia associated with central line infection", "catheter site", "umbilical vein catheterization", "bactermia associated with dialysis line infection"],
            "Reproductive System Infections/Illnesses" : ["genital infection", "genital disease", "genital disorder", "genital tract infection", "genital tract disease", "genital tract disorder", "genital ulcer", "genital warts",  "genital candidiasis", "vulvovaginal candidiasis", "vulvovaginal infection", "vulvovaginal disease", "vulvovaginal disorder", "vulvar infection", "vulvar disease", "vulvar disorder", "salpingitis and oophoritis", "Oophoritis", "salpingitis", "bartholin gland infection", "bartholinitis", "vulvitis", "vaginitis", "cervicitis", "endometritis", "tubo-ovarian abscess", "gynaecological infection", "placental infection", "chronical prostatitis", "prostatitis", "prostate infection", "prostate disease", "prostate disorder", "endometrium", "endometrial infection", "endometrial disease", "endometrial disorder", "endometrial infection", "endometrial disease", "endometrial disorder", "cervical infection", "cervical disease", "cervical disorder", "cervical infection", "cervical disease", "cervical disorder", "vaginal infection", "vaginal disease", "vaginal disorder", "cervicitis", "chronic prostatitis", "right hydrocele", "hydrocele", "left hydrocele", "adenomyosis", "endometritis", "orchiepididymitis", "orchitis", "epididymitis", "prostatitis", "vulvovaginits", "fibroid uterus"],
            "Lymph Node Infections" : ["lymphadenitis", "lymphadenopathy", "lymph node infection", "lymph node disease", "lymph node disorder", "lymphatic infection", "lymphatic disease", "lymphatic disorder", "lymphatic filariasis", "lymphatic filariasis infection", "lymphatic filariasis disease", "lymphatic filariasis disorder", "Lymph node abscess", "lymph node abscess infection", "lymphatic filariasis abscess", "adenitis", "bubo", "adenoflegmon"],
            "Implant-associated Infections": ["implant infection", "implant-associated infection", "prosthetic joint infection", "prosthetic device infection", "orthopedic implant infection", "orthopedic device infection", "surgical site infection", "surgical wound infection", "surgical implant infection", "surgical device infection", "surgical mesh infection", "surgical prosthesis infection", "surgical implant-associated infection", "surgical device-associated infection", "surgical mesh-associated infection", "surgical prosthesis-associated infection", "wound humeral stem explant", "wound explant humeral stem", "Infected hemiarthroplasty", "orthopedic device related infection", "ODRI", "hip prothesis biopsy", "human infected knee prosthesis", "vascular prosthesis infection", "vascular prosthesis-associated infection", "vascular device infection", "periprosthetic liquid", "periprosthetic material", "vascular prosthesis", "vascular prosthesis infection", "prosthetic joint infection_hip", "prothetic joint infection", "prosthetic joint infection_knee", "prosthetic joint infection_ankle", "prosthetic joint infection_shoulder", "prosthetic joint infection_elbow", "prosthetic joint infection_wrist", "prosthetic joint infection_finger", "prosthetic joint infection_toe", "prosthetic joint infection_foot", "prosthetic joint infection_hand", "prosthetic joint infection_spine", "prosthetic joint infection_sternum", "prosthetic joint infection_rib", "prosthetic joint infection_clavicle", "prosthetic joint infection_scapula", "prosthetic joint infection_pelvis", "chronic prosthetic infection", "prosthesis-related infections"],
            "Bursitis" : ["bursitis", "infection of bursa", "septic bursitis", "bursa abscess", "septic arthritis/bursitis",], 
        },

        "Parasitic Infections" : {
            "Pinworm Infections" :["Enterobiasis", "pinworm infections", "enterobius vermicularis", "enterobius infection", "enterobiasis vermicularis", "enterobiasis infection", "pinworm disease", "pinworm infection", "pinworm disease", "enterobius disease"],
            "Babesia Infection" : ["Babesiosis", "babesia", "babesia infection", "babesia microti", "babesia divergens", "babesia bovis", "babesia bigemina", "babesia caballi", "babesia equi", "babesia gibsoni", "babesia canis", "babesia vogeli", "babesia duncani", "babesia venatorum", "babesia odocoilei", "babesia conradae", "babesia rodhaini", "babesia motasi", "babesia bennetti", "babesia felis", "babesia felis"]
        },

        "Viral Infections (VIs)": {
            "COVID-19": ["Corona patient","covid patient","corona","covid","covid-19", "covid-19 disease", "covid-19 virus identified", "covid-19, virus identified", "covid-19 pneumonia", "covid-19 infection", "covid-19 with pneumonia", "severe covid-19", "covid-19, virus identified", "long covid", "post-acute sequelae of covid-19", "coronavirus infection", "sars-cov-2 interstitial pneumonia", ],
            "Influenza and other Respiratory Viral Infections": ["h1n1", "influenza", "avian influenza h7n9 virus infection", "respiratory viral infection", "influenza infection", "seasonal flu", "influenza a virus infection", "influenza b virus infection", "swine flu", "bird flu", "flu-like illness", "h7n9", "h3n2", "h5n1", "h5n6", "h7n7", "h9n2", "influenza a virus", "influenza b virus", "influenza c virus", "influenza d virus", "respiratory syncytial virus infection", "parainfluenza virus infection", "human metapneumovirus infection", "avian influenza", "swine influenza", "seasonal influenza", "influenza-like illness", "viral pneumonia", "viral respiratory infection", "acute respiratory infection", "acute febrile respiratory illness", "acute febrile respiratory syndrome", "acute febrile respiratory disease", "acute febrile respiratory infection", "acute febrile respiratory syndrome", "acute febrile respiratory disease", "h1n1 pneumoniae"],
            "Encephalitis and Meningitis": ["viral encephalitis", "herpes simplex virus encephalitis", "enteroviral encephalitis", "meningitis caused by viruses", "japanese encephalitis", "west nile virus encephalitis", "arboviral encephalitis", "viral meningitis", "aseptic meningitis"], 
            "Hepatitis": ["hepatitis a", "hepatitis b", "hepatitis c", "viral hepatitis", "acute hepatitis a", "chronic hepatitis b", "hepatitis e virus infection", "hepatitis d", "non-a non-b hepatitis", "viral liver inflammation", "severe hepatitis", "severe hepatitis", "hcv-related liver cirrhosis"], 
            "Parvovirus" : ["parvovirus b19", "parvovirus infection", "human parvovirus", "parvovirus b19 infection", "fifth disease", "erythema infectiosum", "slapped cheek syndrome", "parvovirus b19 disease", "parvovirus infection"],
            "Human Immunodeficiency Virus (HIV)": ["hiv", "hiv infection", "aids", "human immunodeficiency virus", "hiv disease", "advanced hiv", "asymptomatic hiv infection"], 
            "Herpesvirus Infections": ["herpes simplex", "herpes simplex virus", "cold sores", "genital herpes", "human herpesvirus", "herpes zoster", "shingles", "varicella-zoster virus infection", "chickenpox", "epstein-barr virus infection", "infectious mononucleosis", "cytomegalovirus infection"], 
            "Measles and Rubella": ["measles", "rubeola", "rubella", "german measles", "measles virus infection", "rubella virus infection", "congenital rubella syndrome"],
            "Mumps": ["mumps", "mumps virus infection", "epidemic parotitis", "parotid gland swelling"], 
            "Arboviral Infections": ["dengue", "dengue fever", "zika virus", "chikungunya", "yellow fever", "west nile virus", "arboviral fever", "arboviral infection"], 
            "Rabies": ["rabies", "rabies virus infection", "hydrophobia", "lyssavirus infection"], 
            "Poliomyelitis": ["poliomyelitis", "polio", "acute poliomyelitis", "poliovirus infection"], 
            "Ebola and Other Hemorrhagic Fevers": ["ebola", "ebola virus disease", "marburg virus", "lassa fever", "hemorrhagic fever", "crimean-congo hemorrhagic fever", "viral hemorrhagic fever"], 
            "Uncategorized Viral Infections": ["viral fever", "viral infection", "unspecified viral infection", "chronic viral infection", "systemic viral infection", "viral syndrome", "fever of viral origin", "unknown viral disease", "nonspecific viral infection", "viral infection with rash","viral exanthema", "subclinical viral infection", "low-grade viral infection",  "virus-associated fever", "viral rash", "transient viral illness", "viral flare-up", "viremia", "mild viral illness", "unspecified viral condition", "febrile illness of viral nature", "self-limiting viral infection", "latent viral infection", ]
        },
        "Fungal Infections": {
            "Cryptococcal and Fungal Meningitis": ["cryptococcal meningitis", "candida meningitis", "aspergillosis", "C. neoformans", "meningitis caused by fungi", "fungal brain infection", "cryptococcus infection"],
            "Fungemia": ["fungemia", "fungal sepsis", "fungal bloodstream infection", "invasive candidiasis", "candidemia", "systemic fungal infection", "fungal septicemia"],
            "Opportunistic Fungal Infections": ["opportunistic fungal infection", "Fungal Infections", "candidiasis", "aspergillosis", "mucormycosis", "pneumocystis pneumonia", "fungal infections in immunocompromised patients", "HIV-associated fungal infections", "pulmonary aspergillosis", "pulmonary candidiasis"],
            "Superficial Fungal Infections": ["tinea corporis", "tinea pedis", "tinea capitis", "ringworm", "athlete's foot", "jock itch", "cutaneous mycoses", "superficial mycosis"],
            "Systemic Fungal Infections": ["histoplasmosis", "blastomycosis", "coccidioidomycosis", "paracoccidioidomycosis", "fungal sepsis", "systemic candidiasis", "invasive aspergillosis", " invasive fungal infections", "invasive fungal disease"],
            "Dermatophyte Infections": ["dermatophytes", "epidermophyton infection", "microsporum infection", "trichophyton infection"],
            "Yeast Infections": ["candidiasis", "thrush", "vaginal yeast infection", "oral candidiasis", "esophageal candidiasis", "cutaneous candidiasis"]
        },

        "Majorly operated": {
            "Transplant Recepient": ["Kidney Transplant", "Renal Transplant", "kidney replacement", "Liver Transplant", "Hepatic Transplant", "liver replacement", "Heart Transplant", "Cardiac Transplant", "heart replacement", "Lung Transplant", "Pulmonary Transplant", "lung replacement", "Pancreas Transplant", "pancreatic replacement", "Bone Marrow Transplant", "Hematopoietic Stem Cell Transplant", "HSCT", "Corneal Transplant", "Keratoplasty", "corneal graft", "Intestinal Transplant", "Small Bowel Transplant", "intestinal graft", "kidney transplantation recipient", "transplanted organ and Biopsy"],
            "Transplant Donor" : ["Kidney Donor", "Renal Donor", "liver donor", "Liver Donor", "Heart Donor", "Cardiac Donor", "lung donor", "Lung Donor", "Pancreas Donor", "Bone Marrow Donor", "Hematopoietic Stem Cell Donor", "HSCT Donor", "Corneal Donor", "Intestinal Donor", "donor organ", ],
            "Amputations" : ["Amputation", "Limb Removal", "Disarticulation", "Severing", "Stump Formation", "Limb Loss", "Excision", "Resection", "Medical Amputation", "Traumatic Amputation", "Surgical Amputation", "Dismemberment", "Extremity Removal", "Limb Resection", "Limb Severance", "Amputee Condition", "Limb Sacrifice", "amputation wound", "amputation wound infection", ], 
            "Cardiovascular Surgeries": ["Coronary Artery Bypass Grafting", "CABG", "heart bypass surgery", "Heart Valve Replacement", "valvular heart surgery", "Aneurysm Repair", "aortic aneurysm surgery", "Pacemaker Implantation", "cardiac pacemaker surgery", "abdominal aortic dissection", "cardiovascular surgery", "aortic dissection surgery", "cardiac surgery", "vascular surgery"],  
            "Neurological Surgeries": ["Brain Tumor Removal", "neurosurgical tumor resection", "Spinal Fusion", "spinal stabilization surgery", "Deep Brain Stimulation", "DBS", "neurological stimulator implant", "Craniotomy", "brain surgery", "Neurosurgery"],  
            "Orthopedic Surgeries": ["Total Knee Replacement", "TKR", "knee arthroplasty", "Total Hip Replacement", "THR", "hip arthroplasty", "Spinal Fusion", "vertebral fusion", "Arthroscopy", "minimally invasive joint surgery", "Shoulder arthroplasty", "arthroplasty", ],  
            "Gastrointestinal Surgeries": ["Gastric Bypass Surgery", "bariatric surgery", "Colectomy", "colon removal surgery", "Appendectomy", "appendix removal", "abdominal surgery", "ileostomy", "colonostomy", "colonostomy", "gastrostomy", "Percutaneous Endoscopic Gastrostomy", "PEG", "gastrostomy tube placement", "cholecystectomy", "bariatric surgery", "gastrointestinal surgery", "gastrointestinal bypass surgery", "gastrointestinal resection", "gastrointestinal reconstruction", "gastrointestinal diversion surgery", "gastrointestinal bypass", "postoperative of sigmoid colon torsion", "postoperative of sigmoid colon torsion with obstruction", "gastrostomy border", "visceral and digestive surgery"],
            "Biliary Tract Surgeries" : ["Cholecystectomy", "gallbladder removal", "hepatobiliary surgery", "postcholecystectomy syndrome", ],   
            "Oncological Surgeries/Treatments": ["Mastectomy", "breast cancer surgery", "Lobectomy", "lung cancer surgery", "Prostatectomy", "prostate cancer surgery", "Hysterectomy", "uterus removal", "carcinoma resection", "radical hysterectomy for cervial carcinoma", "postoperative breast cancer", "postoperative gastric cancer", "post excision carotid body tumor", "non cf tumor irradiated", "tumor irradiated", "postoperative hilar cholangiocarcinoma"],  
            "Urological Surgeries": ["Nephrectomy", "kidney removal surgery", "Cystectomy", "bladder removal", "Prostatectomy", "prostate removal", "urology surgery", "urostomy", "pyeloplasty", "renal surgery", "surgery for rt vuj calculi"],
            "Gynecological Surgeries": ["Hysterectomy", "uterus removal surgery", "Oophorectomy", "ovary removal surgery", "Salpingectomy", "fallopian tube removal surgery", "panhysterectomy", "total hysterectomy", "myomectomy", "fibroid removal surgery", "bleeding after endometrial polyp surgery", "endometrial polyp removal", "endometrial polyp surgery", "postoperative endometrial polyp excision", "cesarian"],
            "Ocular Surgeries": ["Cataract Surgery", "lens replacement surgery", "Corneal Transplant", "keratoplasty", "Retinal Detachment Repair", "retinal surgery", "Glaucoma Surgery", "intraocular pressure surgery"],
            "Dental Surgeries": ["Tooth Extraction", "dental surgery", "Wisdom Tooth Removal", "third molar extraction", "Root Canal Treatment", "endodontic therapy", "Dental Implant Surgery", "dental implant placement"],
            "Plastic and Reconstructive Surgeries": ["Rhinoplasty", "nose surgery", "Breast Augmentation", "breast enhancement surgery", "Liposuction", "liposculpture", "Facelift Surgery", "rhytidectomy", "Abdominoplasty", "tummy tuck surgery"],
            "Trauma Surgeries": ["Trauma Surgery", "trauma management", "Emergency Surgery", "trauma care", "Acute Care Surgery", "trauma and emergency surgery", "Trauma Resuscitation", "trauma stabilization", "Trauma Reconstruction", "trauma repair", "Trauma Decompression", "trauma decompression surgery", "post-traumatic surgery"],
            "Thoracic Surgeries": ["Thoracotomy", "chest surgery", "Lung Resection", "lung removal surgery", "Mediastinal Surgery", "mediastinal tumor surgery", "Pleural Surgery", "pleural drainage surgery"],
            "Miscellaneous Surgery Sites": ["Surgical Drain Placement", "drainage surgery", "Fistula Repair", "fistula closure surgery", "Abscess Drainage", "abscess incision and drainage", "Wound Debridement", "wound cleaning surgery", "Surgical Biopsy", "biopsy surgery", "Skin Grafting", "skin graft surgery", "Tissue Expansion Surgery", "tissue expander placement", "surgical wound", "surgical incision", "surgical site", "surgical procedure", "surgical intervention", "surgical management", "surgical treatment", "surgical approach", "infected incision", "post-surgical wound infection", "surgical injury", "surgical wound", "postoperative wound infection", "surgical site infection", "surgical puncture", "surgical wound infection", "postoperative infection", "surgical injuries", "operative incision infection", "wound surgical", "fasciotomy", "general practice surgery", "general surgery", "intraoperative injury", "surgical site infection", "intrasurgical", "post-op", "postoperative maxillary cyst", "secretion surgical", "Miscellaneous Surgery Sites", "surg aspir", "surgical bandage", "surgical puncture fluid", "surgical specimen from rectum", "several surgeries"],
        },

        "Trauma and Injury": {
            "Head Injury and Traumatic Brain Injury": ["head injury", "concussion", "traumatic brain injury", "brain trauma", "craniofacial injury", "intracranial hemorrhage", "cerebral contusion", "cerebral hemorrhage", "brainstem hemorrhage", "hemorrhagic shock", "severe open traumatic brain injury", "craniofacial injury", "closed head injury", "subdural hematoma", "epidural hematoma", "skull fracture", "diffuse axonal injury", "Cerebral hemorrahge", "Cerebral infarction", "cerebral trauma", "hematencephalon", "cerebral bleeding", "cerebral hematoma", "post-traumatic hydrocephalus", "basilar artery occlusion", "hypoxic ischemic encephalopathy", "Cerebral artery occlusion", "traumatic intracranial hematoma", "traumatic subdural hemorrhage", "closed craniocerebral injury"],
            "Fractures and Bone Injuries": ["fracture", "fracture healing", "bone fracture", "fracture of thoracic vertebra", "fracture of occipital bone", "fracture of the cervical spine", "pelvic fracture", "compound fracture", "displaced fracture", "fracture of femur", "spinal fracture", "broken bone", "stress fracture", "hairline fracture", "rib fracture", "fracture of radius", "fracture of tibia", "fracture of ulna", "greenstick fracture", "osteoarthritis", "intertrochantric fracture", "pelvic fracture", "femoral shaft fracture", "intertrochanteric injury of femur", "femoral neck fracture", "postoperative femoral trochanter fracture"],
            "Burns and Soft Tissue Injuries": ["burn injury", "severe burns", "burnt patient", "electrical burns", "chemical burns", "third-degree burns", "second-degree burns", "thermal burns", "burn shock", "scalding injury", "skin burns", "superficial burns", "soft tissue injury", "abrasion", "laceration", "contusion", "infected burn", "post-burn infection", "burnt patient", "burn injury", "electrical burn", "burn disease", "third-degree scald", "second-degree scald", "thermal scald", "chemical scald", "electrical scald", "chemical burn", "electrical burn", "thermal burn", "superficial burn", "partial-thickness burn", "full-thickness burn", "deep partial-thickness burn", "superficial partial-thickness burn", "deep second-degree burn", "superficial second-degree burn", "burned", "burn"],
            "Other Injuries": ["Trauma", "multiple trauma", "polytrauma", "concomitant injury", "traumatic shock", "crush injury", "penetrating injury", "blunt trauma", "trauma-related hemorrhage", "vascular trauma", "soft tissue trauma", "gunshot wound", "stab wound", "Carbon monoxide poisoning", "multiple injuries", "hematoma", "unwitnessed fall", "finger injury", "compound trauma", "blast", "bomb blast", "multiple injuries", "external injury", "multiple site damage", "gun shots", "gunshot injury", "stab injury", "stab wound", "gunshot wound", "penetrating trauma", "blunt force trauma", "crush syndrome", "crush injury", "blast injury", "blast trauma", "traumatic amputation", "traumatic dislocation", "traumatic rupture", "electrical injury", "hematoma", "injured elbow", "traumatic discharge"],
            "Sports Injuries": ["sports injury", "ACL tear", "ligament injury", "rotator cuff injury", "tendon tear", "meniscal tear", "sprained ankle", "hamstring injury", "shoulder dislocation", "knee dislocation", "stress injury", "sports-related fracture", "Achilles tendon rupture"],
            "Spinal and Neck Injuries": ["spinal cord injury", "cervical injury", "neck trauma", "vertebral fracture", "herniated disc", "spinal shock", "thoracic injury", "lumbar injury", "whiplash injury", "cervical strain", "compression fracture", "spinal instability"],
            "Abdominal and Thoracic Injuries": ["abdominal trauma", "thoracic trauma", "chest injury", "blunt abdominal injury", "penetrating abdominal injury", "diaphragmatic rupture", "liver laceration", "splenic injury", "thoracic contusion", "rib injury", "pulmonary contusion", "flail chest"],
            "Animal Bites and Stings": ["animal bite", "dog bite", "cat bite", "snake bite", "insect sting", "spider bite", "scorpion sting", "animal-related injury", "wildlife encounter injury", "venomous bite", "rabid animal bite", "animal attack", "insect bite", "bee sting", "wasp sting", "ant bite", "tick bite", "mosquito bite", "fire ant sting", "brown recluse spider bite", "black widow spider bite", "scorpion sting", "jellyfish sting", "stingray sting", "centipede bite", "horsefly bite", "deer tick bite", "chigger bite", "bite wound", "sting wound", "animal-related injury", "animal attack", "animal encounter injury", "animal-related trauma", "animal-inflicted injury", "animal bite injury", "animal bite trauma", "animal bite wound", "animal bite syndrome", "animal bite infection", "animal bite abscess", "animal bite cellulitis", "animal bite laceration", "animal bite puncture wound", "animal bite scratch", "animal bite dermatitis", "animal bite allergy", "animal bite hypersensitivity", "animal bite anaphylaxis", "animal bite immunization", "animal bite prophylaxis", "animal bite treatment", "animal bite management", "animal bite care", "animal bite first aid", "animal bite prevention", "animal bite education", "animal bite awareness", "bite wound tissue", "abscess following dog bite", "fatal septicemia after dog bite"],
        },

        "Abdominal Disorders": {
            "Pancreatic Disorders": ["acute pancreatitis", "pancreatitis", "chronic pancreatitis", "pancreatic necrosis", "pancreatic pseudocyst", "pancreatic abscess", "hereditary pancreatitis", "necrotizing pancreatitis", "exocrine pancreatic insufficiency", "autoimmune pancreatitis", "severe pancreatitis", "pancreatitis colonization", "idiopathic acute pancreatitis", "necrotic pancreatitis", "acute necrotizing pancreatitis", "acute hemorrhagic pancreatitis", "acute edematous pancreatitis", "acute interstitial pancreatitis", "necrotic pancreatic", "pancreatic fluid"],
            "Liver Disorders": ["liver cirrhosis", "hepatitis", "liver benign tumor", "liver metastases", "hepatic failure", "non-alcoholic fatty liver disease", "cirrhosis due to alcohol", "alcoholic hepatitis", "viral hepatitis", "fatty liver", "autoimmune hepatitis", "drug-induced liver injury", "primary biliary cholangitis", "primary sclerosing cholangitis", "hepatic steatosis", "liver fibrosis", "liver abscess", "hepatocellular carcinoma", "liver fluid", "liver insufficiency", "liver lesions", "liver function impairment", "pyogenic liver abscess", "hepatic insufficiency", "post-hepatitis cirrhosis", "hepato-splenomegally", "enlarged liver", "hyperbilirubinemia", "jaundice", "liver lesions", "hepatic abscess", "liver abscess", "cirrhosis", "liver diseases", "autoimmune liver disease", "perihepatic fluid", "abnormal liver function", "hepatic polycyst"],
            "Gastrointestinal Disorders": ["intra-abdominal infection", "intraabdominal", "intra abdominal abscess", "bloody diarrheal patient", "abdominal distension", "abdominal infection","Watery diarrhea","Abdominal crmaps", "complete obstruction", "Stomachache","nsaid induced ulcer","salmonella gastroenteritis","cholera","gastritis", "peptic ulcer", "stomach ulcer", "gastric cancer", "Helicobacter pylori infection", "dyspepsia", "gastroparesis", "Zollinger-Ellison syndrome", "stomach polyps", "gastric bleeding", "acute gastritis", "chronic gastritis", "colitis", "diverticular disease", "anal fissure", "rectal prolapse", "constipation", "diarrhea", "intestinal polyps", "rectal bleeding", "megacolon", "ischemic colitis", "perianal abscess", "irritable bowel syndrome", "IBS", "Crohn's disease", "ulcerative colitis", "inflammatory bowel disease", "diverticulitis", "diverticulosis", "intestinal obstruction", "intestinal ischemia", "Celiac disease", "short bowel syndrome", "intestinal perforation", "small intestinal bacterial overgrowth", "malabsorption syndrome", "gastroesophageal reflux disease", "GERD", "esophagitis", "Barrett's esophagus", "esophageal cancer", "achalasia", "esophageal stricture", "esophageal varices", "esophageal perforation", "hiatal hernia", "esophageal spasm", "gastroenteritis", "infectious diarrhea", "food poisoning", "Clostridium difficile infection", "norovirus infection", "bacterial enteritis", "viral gastroenteritis", "parasitic infections", "amebiasis", "giardiasis", "functional dyspepsia", "functional abdominal pain", "functional constipation", "functional diarrhea", "chronic idiopathic constipation", "abdominal bloating", "rumination syndrome", "functional nausea and vomiting", "gastrointestinal bleeding", "hematemesis", "melena", "ascites", "abdominal pain", "intestinal perforation", "mesenteric ischemia", "intestinal fistula", "gastrointestinal stromal tumors", "intestinal malrotation", "peritonitis", "hirschsprung disease", "subcutaneous intestinal obstruction", "gut", "bloody diarrhea", "intraabdominal infection", "appendicitis", "perforated appendicitis", "necrotizing enterocolitis", "intestinal obstruction", "e. coli gastroenteritis", "dysentery", "colitis", "ulcerative colitis", "hemorrhagic colitis", "peritonitis disease", "peritonitis", "hematochezia", "non-ulcer dyspepsia", "dyspepsia", "perisplenic abscess", "acute diffuse peritonitis", "intersphinteric abscess", "acute gastrojejunal ulcer with hemorrhage", "chronic or unspecified peptic ulcer", "chronic peptide ulcer", "chronic gastric ulcer", "chronic duodenal ulcer", "chronic gastric ulcer with hemorrhage", "chronic gastric ulcer with obstruction", "chronic gastric ulcer with perforation", "chronic gastric ulcer with hemorrhage and perforation", "chronic gastric ulcer with hemorrhage and obstruction", "chronic gastric ulcer with obstruction and perforation", "chronic gastric ulcer with hemorrhage, obstruction, and perforation", "chronic duodenal ulcer with hemorrhage", "chronic duodenal ulcer with obstruction", "chronic duodenal ulcer with perforation", "shigellosis", "upper gastrointestinal hemorrhage", "congenital pouch colon", "type 4 congenital pouch colon", "enterocolitis", "diarrhea from food poisoning", "rectal-space occupying lesion", "colonic-space occupying lesion", "abdominal cavity infection", "celiac disease", "esophagus aspirate", "esophagus brush", "primary ciliary dyskinesia", "intrabdominal sepsis", "abdominal sepsis", "g-tube drainage", "travellers diarhhea", "infant with community-acquired diarrhea", "infant with diarrhea", "intraperitoneal", "small bowl obstruction", "gastrointestinal system disease", "enteritis", "phlegmonous appendicitis", "acute peritonitis", "diffuse peritonitis", "diffuse fibrinopurulent peritonitis", "acute gangrenous-perforated appendicitis", "alimentary tract hemorrhage", "upper gastrointestinal hemorrhage", "food poisoning fatality", "IAI", "acute abdomen", "acute abdominal pain", "acute abdominal syndrome", "acute abdomen syndrome", "acute abdomen disease", "acute abdomen infection", "acute abdomen disorder", "acute abdomen complication", "acute abdomen event", "acute abdomen episode", "acute abdomen case", "acute abdomen situation", "acute abdomen problem", "acute abdomen ailment", "acute abdomen illness", "pouchitis", "ileal pouch inifection"],
            "Gallbladder and Biliary Tract Disorders": ["chronic cholecystitis", "acute suppurative cholangitis", "ASC", "gallstones", "cholecystitis", "cholangitis", "biliary colic", "bile duct obstruction", "gallbladder polyps", "primary biliary cholangitis", "choledocholithiasis", "biliary dyskinesia", "gallbladder cancer", "biliary stricture", "bile infections", "choledochocyst", "acute cholangitis", "choledocholithiasis", "gall bladder infection", "bile duct infection", "gallstones with acute cholecystitis", "calculus of bile duct", "abscess gallbladder", "bile duct obstruction", "bile duct stricture", "biliary obstruction", "biliary atresia", "gallbladder disease", "gallbladder inflammation", "gallbladder sludge", "gallbladder empyema", "gallbladder perforation", "gallbladder wall thickening", "cholecystitis with gallstones", "cholecystitis without gallstones", "cholecystitis with bile duct stones", "cholecystitis with pancreatitis", "biloma", "gall bladder infection", "acute acalculous cholecystitis", "acalculous cholecystitis", "cholelithiasis", "calculus of gallbladder", "obstructive jaundice", "biliary tract infection", "alimentary tract hemorrhage"],
            "Gut Microbiome" : ['human gut',"gut microbiome", "ed gut","gut flora", "gut bacteria", "gut microflora", "intestinal microbiota", "intestinal flora", "intestinal bacteria", "intestinal microflora", "intestinal microbiome", "gut microorganisms", "intestinal microorganisms", "gut microbial community", "intestinal microbial community", "gut microbial diversity", "intestinal microbial diversity", "gut microbial composition", "intestinal microbial composition"],  
        },

        "Cardiac Diseases": {
            "Cardiac Disorders": ["Heart Attack", "Myocardial Infarction", "MI", "Acute Coronary Syndrome", "ACS", "Coronary Artery Disease","Coronary Disease", "CAD", "Ischemic Heart Disease", "IHD", "Atherosclerosis", "Heart Failure", "Congestive Heart Failure", "CHF", "Cardiac Insufficiency", "Arrhythmia", "Irregular Heartbeat", "Dysrhythmia", "Heart Rhythm Disorder", "Atrial Fibrillation", "AFib", "AF", "Irregular Atrial Rhythm", "Bradycardia", "Slow Heart Rate", "Sinus Bradycardia", "AV Block", "Tachycardia", "Fast Heart Rate", "Supraventricular Tachycardia", "SVT", "Ventricular Tachycardia", "VT", "Cardiomyopathy", "Enlarged Heart", "Dilated Cardiomyopathy", "Hypertrophic Cardiomyopathy", "HCM", "Restrictive Cardiomyopathy", "Congenital Heart Defects", "CHD", "Congenital Heart Disease", "Birth Heart Defect", "Heart Valve Disease", "Mitral Valve Disease", "Aortic Stenosis", "Valve Regurgitation", "Valve Prolapse", "Pericarditis", "Inflammation of Heart Sac", "Pericardial Disease", "Viral Pericarditis" ,"Endocarditis", "Infective Endocarditis", "Bacterial Endocarditis", "aortal abscess", "Heart Infection", "Hypertensive Heart Disease", "High Blood Pressure Heart Disease", "HHD", "Rheumatic Heart Disease", "RHD", "Rheumatic Fever Heart Damage", "Pulmonary Hypertension", "PH", "High Lung Blood Pressure", "Pulmonary Arterial Hypertension", "PAH", "Sudden Cardiac Arrest", "SCA", "Cardiac Arrest", "Sudden Heart Failure", "Angina", "Chest Pain", "Stable Angina", "Unstable Angina", "Prinzmetal's Angina", "Myocarditis", "Heart Muscle Inflammation", "Viral Myocarditis", "coronary artery atherosclerosis", "hydropericardium", "Hydropericardium syndrome", "Angara disease", "Litchi heart disease", "congestive heart failure", "cardiac arrhythmia", "cardiac dysrhythmia", "cardiac ischemia", "cardiac tamponade", "cardiogenic shock", "cardiomyopathy", "coronary artery disease", "coronary artery spasm", "coronary artery thrombosis", "coronary artery vasospasm", "coronary heart disease", "coronary insufficiency", "coronary occlusion", "coronary vasospasm", "heart attack", "heart failure", "heart murmur", "heart valve disease", "myocardial infarction", "portal hypertensive gastroenteropathy", "pulmonary hypertension", "sudden cardiac death", "ventricular fibrillation", "ventricular tachycardia", "phg", "ventricular septal defect", "mitral insufficiency", "valvular heart disease", "heart valve abscess", "infective endocarditis", "aortic valve stenosis", "subacute case of bacterial endocarditis", "cardiac respiratory arrest", "acute decompensated heart failure", "cardiovascular disease", "cardiac ischemia"],
            "Vascular Disorders": ["aortic dissecting aneurysm", "Aortic Aneurysm", "Aortic Dilatation", "aortic dissection", "Aneurysm", "AAA", "Thoracic Aneurysm", "Peripheral Artery Disease", "PAD", "Peripheral Vascular Disease", "PVD", "Claudication", "Deep Vein Thrombosis", "DVT", "Blood Clot in Leg", "Venous Thrombosis","Pulmonary Embolism", "PE", "Lung Blood Clot", "Venous Thromboembolism", "VTE","Varicose Veins", "Swollen Veins", "Venous Insufficiency", "Spider Veins", "Carotid Artery Disease", "Carotid Stenosis", "Stroke Risk", "Blocked Neck Artery", "Pseudoaneurysm", "vasculitis", "systemic vasculitis", "levamisole vasculitis", "thrombotic microangiopathy", "embolism", "thrombus", "mesenteric artery thrombosis", "thrombophlebitis", "arteriosclerosis", "abdominal aortic aneurysm", "femoral haematoma"],
            "Circulatory Conditions": ["Shock", "Cardiogenic Shock", "Hypovolemic Shock", "Circulatory Collapse", "Hypertension", "High Blood Pressure", "HTN", "Essential Hypertension","Hypotension", "Low Blood Pressure", "Postural Hypotension", "Orthostatic Hypotension", "Raynaud's Disease", "Raynaud's Phenomenon", "Cold-Induced Vasospasm"]
        },


        "Neurological Disorders": {
            "Cerebrovascular Diseases": ["cerebrovascular accident", "subarachnoid hemorrhage", "spontaneous subarachnoid hemorrhage","intracranial aneurysm", "ischemic stroke", "hemorrhagic stroke", "transient ischemic attack", "cerebral infarction", "brain ischemia", "lacunar stroke", "cerebral venous thrombosis", "carotid artery stenosis", "brain aneurysm rupture", "vascular malformation", "stroke-related paralysis", "basilar ganglion hemorrhage", "non traumatic intracranial hamorrahge", "hypertensive cerebral hemorrhage", "Anoxic Brain Damage", "encephalorrhagia", "brain abscess", "brain hemorrhage", "acute toxic encephalopathy", "brain bleed", "brain stroke"],
            "Neurodegenerative Disorders": ["Alzheimer's disease", "Parkinson's disease", "vascular dementia", "frontotemporal dementia", "Lewy body dementia", "amyotrophic lateral sclerosis", "ALS", "motor neuron disease", "Huntington's disease", "Creutzfeldt-Jakob disease", "prion disease", "progressive supranuclear palsy", "spinal muscular atrophy", "multiple system atrophy", "dementia", "cerebral infarction with hemorrhage"],
            "Epilepsy and Seizure Disorders": ["epilepsy", "seizures", "absence seizures", "tonic-clonic seizures", "partial seizures", "febrile seizures", "status epilepticus", "photosensitive epilepsy", "focal seizures", "temporal lobe epilepsy", "myoclonic epilepsy", "juvenile myoclonic epilepsy", "epileptic encephalopathy"],
            "Headache Disorders": ["migraine", "cluster headache", "tension-type headache", "chronic migraine", "rebound headache", "trigeminal neuralgia", "sinus headache", "thunderclap headache", "post-traumatic headache", "hemiplegic migraine", "new daily persistent headache"],
            "Movement Disorders": ["Parkinsonism", "essential tremor", "dystonia", "chorea", "athetosis", "tardive dyskinesia", "Tourette syndrome", "restless legs syndrome", "myoclonus", "spasticity", "rigidity", "cervical dystonia", "motor dysfunction", "dyskinesia", "abnormalities of gait and mobility", "tetraplegia"],
            "Demyelinating Diseases": ["multiple sclerosis", "MS", "neuromyelitis optica", "acute disseminated encephalomyelitis", "transverse myelitis", "chronic inflammatory demyelinating polyneuropathy", "Guillain-Barré syndrome", "demyelinating neuropathy", "optic neuritis"],
            "Peripheral Nervous System Disorders": ["peripheral neuropathy", "carpal tunnel syndrome", "Guillain-Barré syndrome", "sciatica", "brachial plexus injury", "radiculopathy", "entrapment neuropathy", "nerve compression", "polyneuropathy", "diabetic neuropathy", "charcot arthropathy"],
            "Neuromuscular Disorders": ["myasthenia gravis", "muscular dystrophy", "polymyositis", "dermatomyositis", "inclusion body myositis", "spinal muscular atrophy", "Lambert-Eaton myasthenic syndrome", "congenital myopathy"],
            "Disorders of Consciousness": ["coma", "persistent vegetative state", "minimally conscious state", "brain death", "disorders of consciousness", "locked-in syndrome"],
            # "Infectious Neurological Disorders": ["meningitis", "encephalitis", "brain abscess", "neurocysticercosis", "tuberculous meningitis", "viral encephalitis", "bacterial meningitis", "progressive multifocal leukoencephalopathy", "spinal cord infection"],
            "Spinal Cord Disorders": ["spinal cord injury", "spinal stenosis", "herniated disc", "cauda equina syndrome", "spinal myelopathy", "spinal cord tumor", "syringomyelia", "spinal cord compression", "paraplegia", "quadriplegia"],
            "Uncategorized Neurological Disorders": ["encephalopathy", "hydrocephalus", "normal pressure hydrocephalus", "cerebral palsy", "brain tumor", "cranial nerve disorders", "Arnold-Chiari malformation", "posterior fossa syndrome", "intracranial hypertension", "neuralgia", "neurologic alteration", "cerebral hernia", "subdural edema", "neurological disorder", "intracranial space-occupying lesion", "diencephalic syndrome", "sellar lesions", "Botulism", "infant botulism", "spastic quadriplegic cerebral palsy", "schizophrenia", "schizophrenic disorder", "schizophreniform disorder", "schizoaffective disorder", "schizotypal personality disorder", "schizoid personality disorder", "delusional disorder", "brief psychotic disorder", "substance-induced psychotic disorder", "psychotic disorder due to another medical condition", "psychotic disorder not otherwise specified"],
            "Spinal Disorders": ["spinal stenosis", "spinal disc herniation", "spondylolisthesis", "spinal deformity", "spinal cord compression", "sciatica", "cervical radiculopathy", "lumbar radiculopathy", "spinal instability", "spinal tumors", "spinal fractures", "spinal infections", "spinal arthritis", "spinal degeneration", "kyphoscolosis", "scoliosis", "spinal deformities", "spinal cord injury", "spinal stenosis", "spinal disc disease", "spinal fusion", "vertebral destruction"],
        },

        "Cancer": {
            "Solid Tumors": ["tumor","solid tumor","cancer","lung cancer", "esophageal cancer", "gastric cancer", "prostate cancer", "liver cancer", "colorectal cancer", "stomach cancer", "pancreatic cancer", "breast cancer", "bladder cancer", "ovarian cancer", "kidney cancer", "renal cell carcinoma", "melanoma", "thyroid cancer", "head and neck cancer", "sarcoma", "brain tumor", "glioblastoma", "endometrial cancer", "cervical cancer", "testicular cancer", "soft tissue sarcoma", "neuroblastoma", "cholangiocarcinoma", "nasopharyngeal carcinoma", "adrenocortical carcinoma", "gastrointestinal stromal tumor", "mesothelioma", "retinoblastoma", "carcinoid tumor", "Merkel cell carcinoma", "malignant phyllodes tumor", "angiosarcoma", "clear cell carcinoma", "desmoid tumor", "anaplastic thyroid cancer", "mucinous carcinoma", "vulvar cancer", "peritoneal cancer", "Wilms tumor", "rhabdomyosarcoma", "Ewing's sarcoma", "medulloblastoma", "osteosarcoma", "germ cell tumors", "hepatoblastoma", "bone tumor", "urinary system tumors", "Esophageal cancer", "rhabdoid tumor", "small cell lung carcinoma", "slcc", "adenocarcinoma", "gastric cancers", "colorectal cancer", "colon cancer", "metastases", "metastasised cancer", "metastasised tumor", "primary hepatic carcinoma", "primary carcinoma of the liver", "ampulla of vater cancer", "ampullary cancer", "bladder tumor", "urogenital tumour", "abdomninal cancer", "gallbladder cancer", "bile duct cancer", "intra-abdominal malignancy", "esophagus cancer", "rectal cancer", "malignant neoplasm", "gastric carcinoma", "pancreatic carcinoma", "ovarian carcinoma", "intestinal cancer", "meningiomas", "renal tumor", "gallbladder carcinoma", "hepatocellular carcinoma", "adenocarcinoma of esophagus", "lung malignant tumor", "sclc", "small cell lung carcinoma", "non small cell lung cancer", "nsclc", "laryngeal neoplasm", "laryngeal cancer", "cancer oesophagus", "gallbladder carcinoma", "gastric carcinoma", "oesophageal cancer", "oesophageal carcinoma", "oesophageal neoplasm", "oesophageal malignancy", "oesophageal adenocarcinoma", "oesophageal squamous cell carcinoma", "oesophageal squamous carcinoma", "oesophageal squamous cell neoplasm", "oesophageal squamous cell malignancy", "oesophageal squamous cell tumor", "oesophageal squamous cell tumor", "oesophageal squamous cell tumor", "oesophageal squamous cell tumor", "pleura tumor", "malignant tumor of colon", "buccal mucosal cancer", "squamous cell carcinoma", "duodenal neoplasms", "bladder tumor", "human abdominal tumor", "malignant tumor of the tonsil", "gallbladder carcinoma", "malignant tumor of bladder", "advanced stage adenocarcinoma of lung", "adenocarcinoma of lung", "squamous cell carcinoma - skin", "oral squamous cell carcinoma patient", "saliva of oral squamous cell carcinoma patient", "squamous cell carcinoma of the skin", "squamous cell carcinoma of the head and neck", "squamous cell carcinoma of the oral cavity", "squamous cell carcinoma of the larynx", "squamous cell carcinoma of the pharynx", "squamous cell carcinoma of the esophagus", "squamous cell carcinoma of the cervix", "squamous cell carcinoma of the vulva", "carcinoma of descending colon", "tumor Biopsy of colon cancer", "malignant tumor of kidney", "benign prostatic hyperplasia", "benign prostratic hypertrophy", "vestibular schwannoma", "brain tumor"],
            "Hematologic Cancers": ["plasma cell leukemia", "leukemia", "acute myeloid leukemia", "AML", "acute lymphoblastic leukemia", "chronic myeloid leukemia", "chronic lymphocytic leukemia","lymphoma", "Hodgkin's lymphoma", "non-Hodgkin lymphoma", "multiple myeloma", "myelodysplastic syndrome", "MSD", "plasma cell neoplasm", "hairy cell leukemia", "T-cell lymphoma", "Burkitt's lymphoma", "lymphoblastic lymphoma", "Waldenström's macroglobulinemia", "large B-cell lymphoma", "acute monocytic leukemia", "leukemic", "hematological malignancy", "follicular lymphoma grade ii", "follicular lymphoma grade iii", "chronic lymphocytic leukaemia", "chronic lymphocytic leukaemia with 17p deletion", "chronic lymphocytic leukaemia with 11q deletion", "chronic lymphocytic leukaemia with 13q deletion", "chronic lymphocytic leukaemia with unmutated igvh gene", ]
        },

        "Metabolic and Systemic Diseases": {
            "Kidney Disorders": ["chronic kidney disease", "Kidney stones", "acute kidney injury", "renal insufficiency", "kidney failure", "renal tumor", "end-stage renal disease", "polycystic kidney disease", "dialysis", "glomerulonephritis", "nephrotic syndrome", "acute tubular necrosis", "hydronephrosis", "nephrolithiasis", "renal artery stenosis", "renal vein thrombosis", "kidney transplant", "lupus nephritis", "IgA nephropathy", "rta", "renal tubular acidosis", "Uremia", "proteinuria", "pauci-immune crescentic glomerulonephritis", "picgn", "glomerulonephritis", "pyelonephritis", "acute pyelonephritis", "pyelitis", "polycystic kidney disease", "acute tubulointerstitial nephritis", "renal failure", "Urolithiasis", "kidney stone", "nephrolithiasis", "renal anemia", "chronic renal failure", "capd", "peritoneal dialysis", "hemodialysis", "kidney transplant rejection", "chronic kidney disease stage 1", "chronic kidney disease stage 2", "chronic kidney disease stage 3", "chronic kidney disease stage 4", "chronic kidney disease stage 5", "kidney disease", "double renal cyst", "vesicoureteral reflux", "nephropathy", "Acute renal failure", "huge kidney stones", "the huge stone left kidney", "the huge stone right kidney", "end stage renal failure", "urolithiasis", "nephrolithiasis", "renal colic", "renal stone", "kidney stone disease", "kidney stone disease", "nephrolithiasis", "urolithiasis", "renal calculi", "urinary calculi", "urinary stones", "urinary tract stones", "urinary tract calculi", "kidney infection", "kidney", "kidney isolate", "nephretic drain fluid", "end stage kidney disease", "serum creatinine elevation", "acute post-streptococcal glomerulonephritis", "PSGN", "acute bilateral pyelonephritis", "acute obstructive pyelonephritis"],
            "Endocrine Disorders": ["diabetes", "diabetes mellitus", "type 1 diabetes", "type 2 diabetes", "gestational diabetes", "diabetic foot", "thyroid disorder", "hypothyroidism", "hyperthyroidism", "Cushing's syndrome", "Addison's disease", "Graves' disease", "Hashimoto's thyroiditis", "pituitary adenoma", "acromegaly", "hyperparathyroidism", "hypoparathyroidism", "pheochromocytoma", "adrenal insufficiency", "diabetic ketoacidosis", "hyperaldosteronism", "paraneoplastic syndrome", "diabetes ketosis", "pituitary macroadenoma", "diabetic hyperosmolar coma", "Hyperosmolar Hyperglycemic State", "HHS", "Goiter", "substernal goiter", "diabetic nephropathy", "diabetic retinopathy", "diabetic neuropathy", "diabetic foot ulcer", "diabetic ketoacidosis", "diabetes insipidus", "insulinoma", "thyroiditis", "thyroid nodule", "thyroid storm", "thyrotoxicosis", "hyperparathyroidism", "hypoparathyroidism", "primary hyperparathyroidism", "secondary hyperparathyroidism", "tertiary hyperparathyroidism", "diabetic food ulcer", "diabetic foot infection", "diabetic foot abscess", "DFU"],
            "Hematologic Disorders": ["anemia", "iron deficiency anemia", "megaloblastic anemia", "folate deficiency anemia", "sickle cell anemia", "aplastic anemia", "thalassemia", "polycythemia", "polycythemia vera", "thrombocytopenia", "immune thrombocytopenic purpura", "leukopenia", "neutropenia", "hemolytic anemia", "hemophilia", "disseminated intravascular coagulation", "myelofibrosis", "antiphospholipid syndrome", "extremely severe aplastic anemia", "Hemolytic Anemia", "Thrombocytopenic Purpura", "Thrombotic Thrombocytopenic Purpura", "TTP", "Hemolytic Uremic Syndrome", "HUS", "Sickle Cell Disease", "SCD", "Sickle Cell Crisis", "Sickle Cell Anemia", "Hematologic diseases", "febrile neutropenia"],
            "Metabolic Disorders": ["metabolic syndrome", "obesity", "hyperlipidemia", "hypercholesterolemia", "hypertriglyceridemia", "phenylketonuria", "Gaucher disease", "Fabry disease", "lysosomal storage disorders", "mitochondrial diseases", "Wilson disease", "hemochromatosis", "galactosemia", "Maple Syrup Urine Disease", "ketoacidosis"],
            "Nutritional Disorders": ["malnutrition", "kwashiorkor", "marasmus", "vitamin B12 deficiency", "beriberi", "scurvy", "rickets", "pellagra", "iron overload", "zinc deficiency", "vitamin D deficiency", "vitamin A deficiency", "iodine deficiency disorders", "vitamin d dependency rickets",],
        },

        'Others': {
            "Congenital Disorders" : ["pierre robin syndrome", "gastroschisis", "multiple congenital malformations", "entophthalmia", "multiple congenital malformations", "megalourethra"],
            "Pregnancy-related Disorders" : ["gestational diabetes", "preeclampsia", "eclampsia", "hyperemesis gravidarum", "placenta previa", "placental abruption", "preterm labor", "postpartum hemorrhage", "gestational hypertension", "intrauterine growth restriction", "oligohydramnios", "polyhydramnios", "fetal distress", "chorioamnionitis", "amniotic fluid embolism"],
            "Poisoned/Intoxicated" : ["lead poisoning", "carbon monoxide poisoning", "mercury poisoning", "arsenic poisoning", "alcohol poisoning", "drug overdose", "toxic shock syndrome", "heavy metal toxicity", "chemical exposure", "pesticide poisoning", "organophosphate poisoning", "acetaminophen overdose", "salicylate toxicity", "copper sulphate poisoning", "anxiolytic poisoning", "organophosphorus pesticide poisoning"],
            "Anastomosis" : ["anastomosis"],
        },  


        'Immonological Diseases' : {
            "Autoimmune Diseases": ['Rheumatoid Arthritis', 'RA', 'rheumatoid arthritis', 'Systemic Lupus Erythematosus', 'SLE', 'lupus', 'Type 1 Diabetes', 'juvenile diabetes', 'insulin-dependent diabetes', 'Multiple Sclerosis', 'MS', 'sclerosis', 'Psoriasis', 'skin psoriasis', "Hashimoto's Thyroiditis", "hashimoto's disease", "Graves' Disease", 'hyperthyroidism', 'Celiac Disease', 'celiac sprue', 'gluten intolerance', 'Myasthenia Gravis', 'MG', "Sjogren's Syndrome", "Sjogren's disease", 'Pernicious Anemia', 'Myositis', 'autoimmune liver disease', 'bullous pemphigoid', "autoimmune disease"],
            "Allergic Disorders": ['Asthma', 'asthmic', 'asthmatic', 'asthma', 'Anaphylaxis', 'anaphylactic shock', 'Hay Fever', 'Allergic Rhinitis', 'seasonal allergies', 'Eczema', 'Atopic Dermatitis', 'eczema dermatitis', 'Food Allergies', 'food intolerance', 'Drug Allergies', 'medication allergy', 'Urticaria', 'Hives', 'skin rash', 'Angioedema', 'swelling', 'Contact Dermatitis', 'skin contact allergy', 'Allergic Conjunctivitis', 'eye allergy', 'Allergic Bronchopulmonary Aspergillosis', 'ABPA', 'allergic aspergillosis',],
            "Immunodeficiency Disorders": ['Severe Combined Immunodeficiency', 'SCID', 'bubble boy disease', 'Common Variable Immunodeficiency', 'CVID', 'common variable immunodeficiency disease', 'HIV/AIDS', 'HIV', 'AIDS', 'human immunodeficiency virus', 'acquired immunodeficiency syndrome', 'IgA Deficiency', 'IgA deficiency syndrome', "Bruton's Agammaglobulinemia", "Bruton's disease", 'Hyper IgE Syndrome', "Job's Syndrome", 'hyper-IgE syndrome', 'Wiskott-Aldrich Syndrome', 'WAS', 'DiGeorge Syndrome', '22q11.2 deletion syndrome', 'Chronic Granulomatous Disease', 'CGD', 'granulomatous disease', 'Leukocyte Adhesion Deficiency', 'LAD', 'Chediak-Higashi Syndrome', 'CHS', 'Ataxia-Telangiectasia', 'AT', 'ataxia-telangiectasia syndrome', 'Autoimmune Lymphoproliferative Syndrome', 'ALPS', 'lymphoproliferative syndrome', 'Severe Neutropenia', 'neutropenic', 'Severe Eosinophilia', 'eosinophilic', 'patient immunosuppressed', 'immunosuppressed patient', ],
            "Autoinflammatory Diseases": ['Gout', 'uric acid arthritis', "Still's Disease", "adult-onset still's disease", 'Familial Mediterranean Fever', 'FMF', 'Mediterranean fever', 'Cryopyrin-Associated Periodic Syndromes', 'CAPS', 'familial cold autoinflammatory syndrome', "Behçet's Disease", "Behçet's syndrome", 'Sarcoidosis', 'granulomatosis', 'Granulomatosis with Polyangiitis', 'GPA', 'Wegener granulomatosis', 'Eosinophilic Granulomatosis with Polyangiitis', 'EGPA', 'Churg-Strauss syndrome', 'Polyarteritis Nodosa', 'PAN', 'Kawasaki Disease', 'mucocutaneous lymph node syndrome', 'Takayasu Arteritis', 'Takayasu disease', 'Giant Cell Arteritis', 'temporal arteritis', 'Polymyalgia Rheumatica', 'PMR', "Adult-Onset Still's Disease"],
            "Transplant Rejection": ['Acute Transplant Rejection', 'acute organ rejection', 'Chronic Transplant Rejection', 'chronic organ rejection', 'Graft-versus-Host Disease', 'GVHD', 'graft-host disease', 'acute cellular rejection'],
            "Hypersensitivity Reactions": ['Serum Sickness', 'immune complex disease', 'Contact Dermatitis', 'skin contact allergy', 'Hypersensitivity Pneumonitis', "bird fancier's lung", "Drug-Induced Hypersensitivity Syndrome", 'DIHS', 'drug rash with eosinophilia and systemic symptoms', ],
            "Cancer Immunodeficiencies": ['Leukemia', 'blood cancer', 'Lymphoma', 'lymphatic cancer', 'Multiple Myeloma', 'plasma cell myeloma'],
            
        }

    }

    best_match_score = 0
    best_category = "Unknown"
    best_disease = "Unknown"
    

    try:
        if isinstance(fuzzy_str, str):
            logging.debug(f"Entered the str loop in find_human_disease function with string: `{fuzzy_str}`.")
            # normalizing the string
            part_digit_free = re.sub(r'\d+', '', fuzzy_str)
            part = re.sub(r'[\[\]{}()]', '', part_digit_free)
            part = re.sub(r'"', '', part)
            part = re.sub(r'\'', '', part)
            parts = [part]
            logging.debug(f"Normalized string: {part}")

            for normalized_str in parts:
                normalized_str = normalized_str.lower()
                logging.debug(f"Normalized string: {normalized_str}")
                for categories_of_diseases, subcategories_of_diseases in host_disease_terms.items():
                    for disease, fuzzy_diseases_names in subcategories_of_diseases.items():
                        for disease_name in fuzzy_diseases_names:
                            similarity_score = fuzz.ratio(normalized_str.lower(), disease_name.lower())           
                            if similarity_score > best_match_score:
                                best_match_score = similarity_score
                                best_category = categories_of_diseases
                                best_disease = disease
                                best_matched_string = disease_name
                                normalized_str_matched = normalized_str
                            else:
                                continue
                
            if best_match_score >= 80:
                logging.debug(f"Human Disease identified for string: `{fuzzy_str}`. The identified disease is `{best_disease}` in category `{best_category}` with a score of {best_match_score}. The best matched string was `{best_matched_string}` with `{normalized_str_matched}`.")
                return best_disease, best_match_score
            
            logging.debug(f"No match was found with unsplitted Form. Trying to split the string.")

            parts = re.split(r'[-+:;/_,]|\band\b', part)
            logging.debug(f"Split string: {parts}")
            parts = [part.strip() for part in parts]
            parts = [str(part) for part in parts]
            logging.debug(f"Normalized string: {parts}")

            for normalized_str in parts:
                normalized_str = normalized_str.lower()
                logging.debug(f"Normalized string: {normalized_str}")
                for categories_of_diseases, subcategories_of_diseases in host_disease_terms.items():
                    for disease, fuzzy_diseases_names in subcategories_of_diseases.items():
                        for disease_name in fuzzy_diseases_names:
                            similarity_score = fuzz.ratio(normalized_str.lower(), disease_name.lower())           
                            if similarity_score > best_match_score:
                                best_match_score = similarity_score
                                best_category = categories_of_diseases
                                best_disease = disease
                                best_matched_string = disease_name
                                normalized_str_matched = normalized_str
                            else:
                                continue
                
            if best_match_score >= 80:
                logging.debug(f"Human Disease identified for string: `{fuzzy_str}`. The identified disease is `{best_disease}` in category `{best_category}` with a score of {best_match_score}. The best matched string was `{best_matched_string}` with `{normalized_str_matched}`.")
                return best_disease, best_match_score

            else:
                logging.debug(f"The string: `{fuzzy_str}` didn't match any entry from the database. Highest match was with `{best_category}` with {best_match_score} score. Returning `Unknown`.")
                return "Unknown", 0
        else:   

            logging.debug(f"Given string: `{fuzzy_str}` is not a string, returning `Unknown`.")
            return "Unknown", 0
    
    except Exception as e:
        logging.debug(f"An unexpected error occurred while finding disease for string `{fuzzy_str}`: {e}")
        return "Unknown", 0

def find_human_sample_source(fuzzy_str:str):
    
    human_body_samples = {
        "Blood": ["Human Blood","Blood", "Peripheral Blood", "Whole Blood", "Blood Sample", "Venous Blood", "Capillary Blood", "blood isolate", "Arterial Blood", "Clotted Blood", "Anticoagulated Blood", "Blood Culture", "Blood Component", "Red Blood Cells", "White Blood Cells", "Blood Cells", "aneurysm sac", "neonatal blood", "human neonatal blood", "blood specimen", "hemoculture", "ductus arantii", "ductus venosus", "emocolture", "venous", "human blood culture"],
        "Plasma": ["Human Plasma", "Plasma", "Blood Plasma", "Plasma Sample", "Serum-Free Plasma", "Plasma Fraction",],
        "Serum": ["Human Serum", "Serum", "Blood Serum", "Serum Sample", "Post-Clotting Plasma", "Serum Fluid", "Serum Separation", "Serum Extraction", "Serum Culture"],
        "Urine": [ "Urine specimen", "Human Urine", "Urine", "Mid-Stream Urine", "Mid Stream Urine", "Urinary Sample", "Urine Sample", "Pee", "Urinary Tract Fluid", "Urinary Secretion", "Urinary Excretion", "Urine Culture", "urine clean catch", "clean voided urine", "midstream clean catch", "nephrostomy urine", "patient urine"],
        "Bladder" : ["bladder wall", "bladder", "bladder wall mucus", "urinary bladder", "suprapubic swab", "suprapubic aspirate", "suprapubic aspirate sample", "suprapubic aspirate fluid", "suprapubic aspirate isolate", "suprapubic aspirate culture", "bladder aspirate", "bladder aspirate sample", "bladder aspirate fluid", "bladder aspirate isolate", "bladder aspirate culture", "urinary bladder aspirate"],
        "CSF": ["Human CSF", "CSF", "Cerebrospinal Fluid", "Spinal Fluid", "Lumbar Puncture Fluid", "Cerebral Spinal Fluid", "Brain Fluid", "Spinal Cord Fluid", "Cerebrospinal Sample", "human cerebrospinal fluid", "Ventricular liquor", "vent liquor", "liquor", "cerebrospinal", "cerebro-spinal fluid", "fluid cerebrospinal", "lumbar fluid", "lumbar puncture"],
        "Sputum": ["expectoration","Human Sputum", "Sputum", "Phlegm", "Expectorated Sputum", "Respiratory Sputum", "Bronchial Secretion", "Saliva With Mucus", "Cough Sample", "Sputum Culture", "Mucoid sputum", "Patient Sputum", "sputamentum", "induced sputum", "sputum specimen", "sputum isolate", "mucopurulent sputum", "expectorate", "sputa", "sputa-trachel"],
        "Pus": ["Brain Abscess","Skin Abscess","Abscess","Human wound pus", "Human pus", "Wound pus", "Pus", "Putulent Discharge", "Infected Discharge", "Wound Drainage", "Pus Culture", "Abscess Fluid", "Pus Sample", "Suppurative Discharge", "carbuncle exudate", "Pus aspirate", "abscess aspirate", "purulent material", "superficial abscess", "swab from purulent wound", "purulent wound", "purulent exudate", "purulent discharge", "purulent secretion", "purulent fluid", "purulent drainage", "purulent aspirate", "purulent sample", "pus aspirate", "pus specimen", "renal abscess", "breast abscess", "purulent secretion", "ear pus", "cerebral abscess", "cerebral pus", "pus from brain", "abscess biopsy", "Abscess Swab", "carbancule", "carbuncle pus", "carbuncle abscess", "carbuncle drain", "carbuncle fluid", "carbuncle culture", "carbuncle sample", "groin abscess",  "perirectal abscess", "perianal abscess", "perianal abscess swab", "leg abscess", "rectal abscess", "hand abscess", "skin abscess", "nares abscess", "apostem", "acne pustules", "furuncle", "boils", "furoncle", "furunculosis", "human furuncles", "human purulent material", "human pus aspirate", "human pus specimen", "human pus sample", "human pus culture", "human pus isolate", "human pus swab", "human pus", "human abscess", "puds swab", "purulent abscess", "purulent specimen", "pustule", "pyoderma", "scrotal abscess", "scrotal abs swab", "gluteus furuncle", "gluteal abscess",  "gluteal skin abscess", "gluteal skin furuncle", "leg furncles", "leg abscess", "gluteal boil", "gluteal skin boil", "leg boil", "soft Biopsy abscess", "soft tissue abscess", "suppuration", "tissue_abscess"],
        "Blister" : ["Blister", "Blister fluid", "right shin blister", "shin blister"],
        "Bleb": ["Bleb", "Bleb Fluid", "Bleb Sample", "Bleb Culture", "Subconjunctival Bleb", "Bleb Aspirate", "Bleb Exudate", "Bleb Discharge", "Bleb Secretion", "Bleb Drainage"],
        "Stool": ["human faeces", "Patinet fecal","Patient stool", "faecal","Human Stool", "Stool", "Feces", "Fecal Sample", "Bowel Movement", "Poop", "Fecal Matter", "Stool Specimen", "Bowel Fluid", "Meconium", "Newborn Stool", "Baby Stool", "First Stool", "Meconial Sample", "Baby's First Poop", "Neonatal Stool", "excrement", "excreta", "Fecal material", "Fecal swab", "Feces swab", "Fecal isolate", "fecal swab isolate", "fecal sample", "stool sample", "human stool sample", "human fecal sample", "adult fecal", "diarrheal sample", "human faecal", "human faecal sample", "human-fecal", ],
        "Pleural Fluid": ["Pleural Fluid", "Pleural Effusion", "Pleural Exudate", "Pleural Aspirate", "Fluid In Pleura", "Thoracic Fluid", "Pleural Cavity Fluid", "pleural puncture fluid", "Hydrothorax", "pleural fluid specimen", "Lung Lavage", "Lung Fluid Sample", "Lung Aspirate", "fluid, pleural", "pleu./f", "thoracentesis fld", "plueral fluid isolate", "thoracentesis fluid peritoneal fluid", "pleural effusion fluid", "pleural effusion aspirate", "peritoneal fluid thoracentesis fluid", "thorax drainage", "pleural puncture", "PLF"],
        "Peritoneal Fluid": ["Peritoneal Fluid", "Ascitic Fluid", "Peritoneal Effusion", "Peritoneal Sample", "Fluid In Peritoneum", "Abdominal Cavity Fluid", "seroperitoneum", "Peritoneal Aspirate", "Ascitic Fluid", "Peritoneal Fluid",  "Ascites", "Subhepatic Fluid", "Subhepatic Sample", "Subhepatic Drainage", "Subhepatic Aspirate", "Subhepatic Fluid Collection", "pelvic fluid", "seroperitoneum", "peritoneal liquor", "peritoneal aspirate", "peritoneal exudate", "peritoneal effusion", "peritoneal lavage", "peritoneal wash", "peritoneal drainage", "peritoneal culture", "abdominal ascites fluid", "Peritoneal abscess", "peritoneal", "abdominal pig tail", "pigtail drain", "abdominal dropsy", "fluid - peritoneal", "fluid, peritoneal", "human peritoneal fluid", "human ascitic fluid", "hydroperitoneum", "intraperitoneal", "intraperitoneal fluid", "intraperitoneal drainage", "peritoneal drainage fluid", "peritoneal fluid ascitic", "peritoneal hematoma", "seroperitoneum", "fluid - peritoneal"],
        "Synovial Fluid": ["Synovial Fluid", "Joint Fluid", "Knee Fluid", "Synovial Aspirate", "Synovial Sample", "Synovial Joint Fluid", "Synovial Exudate", "joint aspirate", "Joint", "collar-bone fluid", "Hip joint", "effusion of articular cavity", "elbow fluid", "fluid knee", "fluid wrist", "fluid_synovial", "hip joint fluid", "hip fluid", "hip joint aspirate", "knee aspirate", "bursitis", "knee fluid", "ankle fluid", "elbow fluid", "shoulder fluid", "fluid, synovial", "fluid, joint", "joint aspirate", "knee joint aspiration", "knee joint fluid", "knee synovial fluid", "knee joint fld", "joint fld", "right elbow aspirate", "left elbow aspirate", "right knee aspirate", "left knee aspirate", "right ankle aspirate", "left ankle aspirate", "right shoulder aspirate", "left shoulder aspirate", "right wrist aspirate", "left wrist aspirate", "right hip aspirate", "left hip aspirate", "right joint aspirate", "left joint aspirate", "joint fluid sample", "joint fluid isolate", "joint fluid specimen", "joint fluid culture", "purulent synovial fluid", "r elbow synovial fluid", "l elbow synovial fluid", "r knee synovial fluid", "l knee synovial fluid", "r ankle synovial fluid", "l ankle synovial fluid", "r shoulder synovial fluid", "l shoulder synovial fluid", "r wrist synovial fluid", "l wrist synovial fluid", "r hip synovial fluid", "l hip synovial fluid", "human joint aspirate", "human joint fluid", "human knee aspirate", "human knee fluid", "human elbow aspirate", "human elbow fluid", "human shoulder aspirate", "human shoulder fluid", "right knee fluid", "left knee fluid", "right knee joint", "right knee joint sonicate", "right knee synovial fluid", "left knee joint", "Left knee joint sonicate", "left knee synovial fluid", "right leg aspirate", "right leg fluid", "left leg aspirate", "left leg fluid", "arthritis aspirate", "hip fluid", "hip joint"],
        "Amniotic Fluid": ["Amniotic Fluid", "Amniotic Sac Fluid", "Amniotic Sac Sample", "Fluid From Amniotic Sac", "Amniotic Discharge", "Fetal Fluid", ],
        "Hair": ["Human Hair", "Hair", "Scalp Hair", "Body Hair", "Hair Sample", "Follicular Sample", "Hairs", "Head Hair", "Human Hair", "Hair Shaft", "Hair follicle", "Hair Root", "Hair Bulb", "Hair Culture", "hair sample", "hair isolate", "hair specimen"],
        "Skin": ["Human Skin", "Skin", "Epidermal Sample", "Skin Scraping", "Cutaneous Tissue", "Skin Cells", "Epidermis", "Dermis", "Skin Biopsy", "Foreskin", "human inguinal skin", "leg skin"],
        "Nail": ["Human Nail", "Nail", "Nail Clipping", "Nail Sample", "Toenail", "Fingernail", "Nail Scrapings", "Nail Biopsy", "Finger Nail", "Toe Nail", "feet nail", "toe nail"],
        "Saliva": ["Human Saliva", "Saliva", "Salivary Fluid", "Oral Fluid", "Spit", "Salivary Sample", "Salivary Secretion", "Mouth Fluid", "Buccal Fluid", "unstimulated saliva"],
        "Muscle": ["Human Muscle", "Muscle", "Muscle Tissue", "Muscle Biopsy", "Skeletal Muscle", "Muscular Tissue", "Cardiac Muscle", "Smooth Muscle", "Muscle Sample", "Muscle Fiber", "hamstring", "hypothenar palm"],
        "Tears": ["Human Tears", "Tears", "Ocular Fluid", "Tear Sample", "Eye Secretion", "Tear Secretion", "Lacrimal Fluid", "Tear Drop", "Eye Fluid"],
        "Vitreous Humor": ["Vitreous Humor", "Vitreous Fluid", "Vitreous Sample", "Ocular Vitreous", "Eye Vitreous", "Vitreous Body", "Vitreous Gel", "Vitreous Humor Sample", "eye gel", "Vitreous", "vitreum"],
        "Breast Milk": ["Human Breast Milk","Breast Milk", "Human Milk", "Lactating Fluid", "Nursing Fluid", "Colostrum", "Breastfeeding Fluid", "Milk From Breast", "Mammary Fluid", "milk", "milk sample"],
        "Sweat": ["Sweat", "Perspiration", "Sweat Sample", "Excreted Fluid", "Skin Secretion", "Sweat Glands Secretion", "Perspiratory Fluid", "Sweat Droplets"],
        "Placental": ["Placental Tissue", "Placenta", "Placental Sample", "Placenta Biopsy", "Placental Cells", "placenta Biopsy", "placenta fetal side", "placenta swab"],
        "Umbilical Cord": ["Umbilical Cord", "Umbilical Tissue", "Placental Cord", "Cord Sample", "Umbilical Fluid", "Umbilical Cord Tissue", "Umbilical Cord Biopsy", "Umbilical Cord Cells", "Umblilical cord secretion", "Umbilical Cord Blood", "Cord Blood", "Blood From Umbilical Cord", "Newborn Cord Blood", "Cord Swab", "Umbilical", "umbilicus swab", "swab umbilical", "umbilical cord secretions", "umbilical secretions", "umbilical swab", "umbilicus", "umbilicus swab", "umbilical stump", "umbilical stump tissue", "umbilical stump biopsy"],
        "Semen": ["Male Semen", "Human Semen","Semen", "Seminal Fluid", "Ejaculate", "Sperm Sample", "Spermatozoa", "Male Ejaculate", "Semen Specimen", "Semen Culture"],
        "Gastric": ["Gastric Lavage", "Stomach Wash", "Gastric Fluid", "Gastric Sample", "Gastric Aspiration", "Gastric Irrigate", "Gastric Secretion", "Gastric Fluid", "Stomach Fluid", "Gastric Sample", "Stomach Aspirate", "Gastric Secretion", "Gastric Juice", "gastric aspirate", "stomach content", "vomit", "vomitus", "stomach mucosa", "Gastric Mucosa", "Gastric Biopsy", "gastric", "gastric contents", "human stomach", "Stomach", "stomach biopsy", "stomach content"],
        "Spleen" : ["Spleen", "Spleen sample"],
        "Lung" : ["Lung", "Lung Sample", "Lung Tissue", "Pulmonary Sample", "lung tissue", "lung biopsy", "cystic fibrosis lung", "lung aspirate", "lung lower lobe Biopsy", "lung swab", "lung washing", "lungs", "right lung cavity swab", "right lung swab", "left lung cavity swab", "left lung swab"],
        "Bronchoalveolar Lavage": ["Bronchoalveolar Lavage", "Bal", "Minibal", "Mini Bronchoalveolar Lavage", "Bronchial Lavage", "Alveolar Lavage", 'balf', 'bronchoalveolar lavage fluid', 'bronchial', 'bronchial aspirate', "broncheal washing", "br wash", "bronchoaspirate"],
        "Tracheal": ["Tracheal Aspirate", "Tracheal Secretion", "Respiratory Aspirate", "Tracheal Fluid", "Endotracheal Aspirate", "ET tube aspirate", "endotracheal aspirate", "ta", "etta", "trach asp", "tracheal asp", "transtracheal aspirate sample", "tracheal swab", "treachea swab", "transtracheal aspirate sample", "tracheal secretions", "tracheal wash", "tracheal aspirate", "tracheal aspirate sample", "tracheal aspirate fluid", "tracheal aspirate specimen", "tracheal aspirate isolate", "tracheal aspirate culture", "tracheal intubation", "tracheobronchial aspirate", "endotrach", "endotracheal", "human trachea", "orotracheal aspirate", "orotracheal secretion", "resp endotrach", "secretion from the trachea branch", "trach. aspirate", "trachaeal secretion", "trachael aspirate", "trachea", "tracheal", "tracheobronchial", "tracheal discharge", "tracheal fluid", "tracheal aspirate cultures", "tracheal aspirate isolate", "tracheal aspirates", "tracheal aspiration", "tracheal aspiration fluid", "tracheal secretion", "tracheal specimen", "tracheal suction", "tracheal swab", "tracheal swab", "trachealaspirate", "trans tracheal aspirate", "transtracheal aspirate", "transtracheal aspirate fluid", "traqueal aspirate", "traqueal secretion", "aspirate-tracheal", "swab trachea"],
        "Cartilage": ["Cartilage", "Articular Cartilage", "Cartilage Tissue", "Joint Cartilage", "Cartilage Sample", "Hyaline Cartilage", "Fibrocartilage", "pinna cartilage"],
        "Fat Tissue": ["Fat Tissue", "Adipose Tissue", "Fat Sample", "Lipid Tissue", "Adipose Sample", "Fat Biopsy", "Adipocyte Tissue", "Fat"],
        "Swab": ["Swab", "Surveillance Swab", "Cotton Swab", "Sample Swab", "Swabbing", "Medical Swab", "Clinical Swab", "Swab Sample", "swallow swab", "urogenital", "urogenital tract", "axilla/groin", "axilla/groin swab", "axilla/groin/perianal Swab", "axilla & groin", "axilla + groin", "colostomy bag swab", "colostomy swab", "drainage swab", "discharge swab", "groin and nares swab", "groin and rectal swab", "intra-operative swab", "nares/axilla/groin", "nasal and perineal swab", "nasal perirectal swab", "neck/axilla swab", "nose & groin swab", "nose and groin swab", "nose and rectal swab", "nose/rectal", "nose, throat inguinal swab", "r. ischium swab", "rectal & groin swab", "mucus", "stoma swab", "swab urostoma", "vaginal and perianal swab", "Biopsy swab", "deep Biopsy by swab", "drain swab", "nasal and perineal swab", "nasal perirectal swab"],
        "Nasal": ["Nasal Swab", "Nose Swab", "Nasal Sample", "Nasal Fester cultures", "Nasal Fluid Sample", "Nasal Secretion", "Mucus", "Nasal Discharge", "Runny Nose Fluid", "Nasal Fluid", "Nasal Drainage", "Nasal Culture", "Nose isolate", "Nasal isolate", "nasal passage", "human nose", "nares", "human nares", "Nasal", "nasal", "nasal cavity", "nose", "anterior nares", "ethmoid sinus", "Nasal Lavage", "nasal specimen", "Ethmoid", "human nasal swab", "human nose", "maxillary sinus", "middle meatus", "middle nasal concha", "middle turbinate", "nares swab", "nares -skin", "nostril swab", "nostril", "paranasal sinus", "patient nose", "r_nares", "l_nares"],
        "Throat": ["Throat Swab", "Oropharyngeal Swab", "Oropharyngeal", "Pharyngeal Swab", "Oral Pharynx Swab", "Pharyngeal Culture", "Throat Culture", "Nasopharyngeal Swab", "Nasal Pharynx Swab", "Nasopharyngeal Sample", "nose throat swab", "Throat", "Nasopharngeal", "Oropharyngeal", "pharynx", "Pharyngeal", "nasopharynx", "Throat", "endotracheal throat swab", "pharynx swab", "pharyngeal exudade", "nasopharyngeal aspiration", "nasopharyngeal secretion", "nasopharyngeal secretions", "deep throat", "human pharyngeal swab", "human pharynx", "human throat swab", "nasopharyngeal", "nasopharyngeal aspirate", "nasopharynx", "nasotrach aspirate", "orofaryngeal", "oropharynx", "phanrynx swab", "pharyngeal mucosa", "pharyngeal secretion", "pharyngeal tonsil", "pharynx mucosa", "pharynx swab", "human, np swab", "NP swab", "NP", "OP", "OP swab", "OP aspirate", "OP aspirate sample", "OP aspirate fluid", "OP aspirate specimen", "NP aspirate", "NP aspirate sample", "NP aspirate fluid", "NP aspirate specimen", "NP aspirate isolate", "NP aspirate culture", "OP aspirate isolate", "OP aspirate culture", "OP swab isolate", "OP swab sample", "tonsil swab", "tonsiles swab", "tonsils"],
        "Wound Swab": ["Wound","Wound Swab", "Wound Sample", "Infected Wound Swab", "Wound Culture", "Wound Exudate", "wound discharge", "Wound Drainage", "punctate", "incision", "Ulcer", "Wound isolate", "wound abscess", "wound secretion", "hip wound", "leg wound", "surgical injury", "surgical wound", "feet wound", "foot wound", "ulcer isolate", "hand wound", "left hand wound", "right hand wound", "skin wound", "sternal wound", "sanies", "post surgical secretion", "surgical secretion","sore", "intercostal tube site", "fester", "festering wound", "slough from ulcer", "tracheostoma", "tracheostome secretion", "tracheostomy", "tracheostomy site", "tracheostomy site swab", "tracheostomy swab", "trach site", "trach swab", "drain inserting site", "drain inserting site swab", "drain site swab", "gangrenous lesion", "gangrenous ulcer", "gangrenous wound", "gangrene", "gangrene wound", "gangrene ulcer", "gangrene lesion", "gastrostomy site swab", "g-tube site", "hickmans site swab", "human lesion", "human ulcer", "human wound", "human sore", "human fester", "human festering wound", "human gangrenous lesion", "human gangrenous ulcer", "human gangrenous wound", "human gangrene", "human gangrene ulcer", "human gangrene lesion", "human sore", "ileostomy site", "impetigo lesion", "lesion", "infusaport access site", "insect bite", "j-tube site", "knee lesion swab", "lesion sample", "lesion swab", "mycetoma lesion", "mycetoma on foot", "nephrostomy swab", "p.e.g site swab", "peg site swab", "diabetic foot ulcer", "pej site", "pej site swab", "pej swab", "postoperative transudate", "purulent infected wound", "sanies", "wound sanies", "diabetic food ulcer", "diabetic foot infection", "diabetic foot abscess", "DFU", "gluteal lesion", "gluteal skin lesion", "leg lesion", "skin ulcer", "skin wound", "toe Biopsy ulcer", "ulcus cruris", "drain would fluid", "gastrostomy secretion"],
        "Ear Swab": ["Ear Swab", "Aural Swab", "Ear Canal Swab", "Ear Sample", "Ear Discharge", "Auricular Swab", "middle ear effusion", "Otic fluid", "auditory canal", "Sample from auditory canal", "middle ear", "otic fluid", "ear isolate", "ear abscess", "ear secretion", "ear culture", "ear aspirate", "ear drainage", "ear swab isolate", "ear swab sample", "ear aspirate sample", "ear aspirate fluid", "ear aspirate specimen", "ear aspirate isolate", "ear aspirate culture", "ear", "ear canal", "ear canal secretions", "middle ear effusion fluid", "middle ear effusion", "swab from ear", "swab of ear", "ear secretions", "external auditory canal", "external ear", "ear canal swab", "external ear swab", "human otorrhea", "inner/middle ear", "left ear", "right ear", "middle ear fluid", "middle ear swab", "otic swab", "aural swab", "aural aspirate", "aural aspirate sample", "aural aspirate fluid", "aural aspirate isolate", "aural aspirate culture", "aural secretions", "aural swab sample", "aural swab isolate", "human ear", "human ear canal", "human external ear", "human external auditory canal", "secretion of auditory canal", "swab ear", "swab ear left", "swab ear right", "swab left eyelid"],
        "Eye Swab": ["Eye Swab", "Ocular Swab", "Conjunctival Swab", "Corneal Swab", "Eyelid Swab", "Eye Secretion", "Conjunctival Sample", "Tear Swab", "Ocular Secretion", "Eye Discharge",  "Pinkeye", "Pink Eye", "Conjunctivitis", "Conjunctival Swab", "Conjunctival Sample","Conjunctival Discharge", "Cornea", "Corneal Sample", "Corneal Discharge", "Corneal Secretion", "Corneal Fluid", "Corneal Culture", "Eye isolate", "eye abscess", "eye secretion", "eye culture", "eye aspirate", "eye drainage", "eye swab isolate", "eye swab sample", "eye aspirate sample", "eye aspirate fluid", "eye aspirate specimen", "eye aspirate isolate", "eye aspirate culture", "eye secreta", "eye", "clinical sample from eye infection", "ocular", "conjunctival secretion", "discharge of the right eye", "general eye", "intraocular contents", "left eye", "right eye", "ocular", "ocular cornea", "ocular hordeolum", "cornea secretion"],
        "Rectal Swab": ["Rectal Swab", "Anal Swab", "Rectum Swab", "Rectal Sample", "Anal Sample", "Rectal Discharge", "rectal swab isolate", "rectal", "rectal surveillance swabs", "rectal screen swab", "rectum", "rectal infection", "rectal culture", "rectal aspirate", "rectal drainage", "rectal swab sample", "anal swab sample", "anorectal", "anal smear sample", "anorectal swab", "human rectal sample", "human rectal swab", "newborn, rectal swab", "routine rectal swab"],
        "Perianal Swab": ["Perianal Swab", "Anal Region Swab", "Perianal Sample", "Perianal Discharge", "Anal Margin Swab", "perirectal", "perirectal swab", "perianal", "anal margin", "swab perianal"],
        "Laryngeal" : ["laryngeal suction tube", "laryngeal suction liquid", "layrnx", "laryngeal aspirate", "laryngeal aspirate sample", "laryngeal aspirate fluid", "laryngeal aspirate isolate", "laryngeal aspirate culture", "laryngeal secretions", "larynx swab", "larynx aspirate"],
        "Cervical/Vaginal Swab": ["Vaginal Swab", "Vaginal Sample", "Cervical Swab", "Vaginal Discharge Swab", "Vaginal Fluid Sample", "Cervical Mucus Swab", "Cervical Mucus", "Vaginal Mucus", "Cervical Fluid", "Cervical Discharge", "Cervical Secretion", "cervical abscess", "Endocervical Swab", "labia", "labial", "genital exudate", "cervical", "endocervix", "endocervical swab", "cervix swab", "cervix", "cervical", "ectocervical mucosa", "ectocervical", "ectocervical swab", "cervical swab isolate", "cervical swab sample", "cervical aspirate sample", "cervical aspirate fluid", "cervical aspirate specimen", "cervical aspirate isolate", "cervical aspirate culture", "vaginal swab isolate", "vaginal swab sample", "vaginal aspirate sample", "vaginal aspirate fluid", "vaginal aspirate specimen", "vaginal aspirate isolate", "vaginal aspirate culture", "endocervical", "endocervix", "female genital", "genital female", "genital vagina swab", "high vaginal swab", "human vagina", "human vaginal discharge", "human vaginal swab", "introitus", "leucorrhea", "leucorrhea sample", "leucorrhea swab", "lochia", "lochia discharge", "lochia swab", "maternal birth canal", "uterine cervix", "vagina", "vaginal", "vaginal discharge", "vaginal discharge", "vaginal environment", "vaginal fluid", "vaginal secretions", "vaginal specimen"],
        "Vulval Swab": ["genital vulva", "vulva swab", "vulval secretions", "vulval discharge"],
        "Genital Swab": ["Genital Swab", "Genital Sample", "Genital Discharge", "Genital secretion", "genital", "genital tract", "swab of genitalia"],
        "Urethral Swab": ["Urethral Swab", "Urethra Swab", "Urethral Sample", "Urinary Swab", "Urethral Fluid Sample", "Urethral Discharge", "mucous urethral discharge", "urethral exudate", "urethra", "urethral", "fluid urethra", "urethral secretion"],
        "Penile Swab": ["Penile Swab", "penile discharge", "male genital", "male genital tract", "genital male", "penis swab", "Penile secretion", "penis", "penile lesion"],
        "Scrotal/Testicular": ["scrotal swab", "scrotum swab", "scrotum", "scrotal", "testicle", "testis"],
        "Axillary Swab": ["Axillary Swab", "Armpit Swab", "Underarm Swab", "Axilla Swab", "Axillary Fluid Sample", "Axilla", "Axilla samples", "Axilla isolates", "armpit", "axillary smear"],
        "Oral Swab": ["Oral Swab", "Mouth Swab", "Buccal Swab", "Oral Cavity Swab", "dental biofilm", "oral cavity", "Salivary Swab", "Oral Fluid Sample", "Gingival Swab", "Gum Swab", "mouth", "Oral Mucosa Swab", "Gums Swab", "Buccal Swab", "Cheek Swab", "Mouth Swab", "Buccal Cavity Swab", "Oral Mucosa Sample", "Buccal mucosa", "coronal", "gingival", "oral mucosa", "oral cavity swab", "Oral", "periodontal pocket", "subgingival dental plaque", "subgingival plaque", "subgingival", "subgiingival swab", "dental calculus", "tartar", "gingiva", "gingival crevice", "gingival margin", "gingival sulcus", "gingivitis sample", "human oral cavity", "human orapharyngeals", "cheek", "l_cheek", "r_cheek", "buccal mucosa", "buccal mucosa swab", "buccal mucosa sample", "buccal mucosa aspirate", "buccal mucosa drainage", "buccal mucosa secretion", "buccal mucosa isolate", "buccal mucosa culture", "subgiingival plaque", "subgiingival plaque swab", "oral", "oral mucosal discharge", "oral cavity swab", "oral cavity- mouth", "oral mucosal discharge", "oral swab sample", "peridontal pocket", "periodontal lesion", "dental caries", "caries swab", "plaque swab", "r_cheek", "l_cheek", "smooth surface caries of tooth", "swab of buccal mucosa", "swab oral", "tongue", "tongue coating", "tongue dorsum", "tongue oral cavity", "tongue scraping", "human tooth surface", "whole mouth swab"],
        "Root Canal": ["Root Canal", "Root Canal Sample", "Root Canal Discharge", "Root Canal Secretion", "Root Canal Fluid", "Root Canal Aspirate", "Root Canal Culture", "Root Canal isolate", "root canal aspirate", "root canal aspirate sample", "root canal aspirate fluid", "root canal aspirate specimen", "root canal aspirate isolate", "root canal aspirate culture", "root caries"],
        "Skin Swab": ["Skin Swab", "Dermal Swab", "Epidermal Swab", "Cutaneous Swab", "Skin sample", "Skin isolate", "Skin", "skin lesion", "skin culture", "skin scrape", "skin biopsy", "skin tissue", "skin aspirate", "skin drainage", "skin swab isolate", "skin swab sample", "shin swab", "skin exudate", "alar crease", "alar crease swab", "ankle swab", "antecubital fossa swab", "back", "back swab", "finger swab", "toe swab", "foot swab", "forearm swab", "forehead swab", "foreskin swab", "gluteal crease", "gluteal sulcus", "gluteal fold", "great toe swab", "hand swab", "skin surface", "human skin surface", "healthy toe", "hip swab", "knee swab", "l. elbow swab", "l_index", "index", "r_index", "l_index finger", "l_palm", "l_thumb", "l_thumb finger", "l_thumb finger swab", "l_thumb swab", "l_thumb finger swab", "l_wrist", "l_wrist swab", "l_wrist swab", "r_index finger", "r_index finger swab", "r_palm", "r_palm swab", "r_thumb", "r_thumb finger", "r_thumb finger swab", "r_thumb swab", "r_thumb finger swab", "left arm swab", "arm swab", "wrist swab", "left foot swab", "foot swab", "feet swab", "right foot swab", "right feed swab", "left feed swab", "left heel swab", "right heel swab", "heel swab", "left hip swab", "right hip swab", "hip swab", "left knee swab", "right knee swab", "knee swab", "left leg swab", "right leg swab", "leg swab", "left palm swab", "right palm swab", "palm swab", "left wrist swab", "right wrist swab", "wrist swab", "left elbow swab", "right elbow swab", "elbow swab", "left shoulder swab", "right shoulder swab", "shoulder swab", "left thumb finger swab", "right thumb finger swab", "limb swab", "lip swab", "male behind the ears", "behind the ears", "middle finger swab", "ring finger swab", "index swab", "index finger swab", "neck swab", "neck skin swab", "plantar swab", "r_index", "l_index", "r_palm", "r_thumb", "l_palm", "l_thumb", "rash", "rash swab", "retroauricular crease", "right elbow swab", "right index finger swab", "right palm swab", "right thumb finger swab", "right wrist swab", "right heel swab", "left elbow swab", "left index finger swab", "left palm swab", "left thumb finger swab", "left wrist swab", "left heel swab", "right forehead swab", "right hand palm swab", "left forehead swab", "left hand palm swab", "facial acne", "facial cutaneous", "flank swab", "forehead swab", "frontal scalp", "frontal scalp swab", "human skin swab", "shin swab", "skin exudate", "skin swab sample", "surface of healthy skin", "surface of skin", "toe web", "toe web space swab", "toe web-space skin", "toes swab", "toeweb", "upper arm swab", "upper front left leg swab", "upper lip swab", "volar forearm swab", "glabell"],
        "Mammilla Swab": ["mammilla swab", "nipple swab", "nipple aspirate", "nipple aspirate sample", "nipple aspirate fluid", "nipple aspirate isolate", "nipple aspirate culture", "nipple secretions", "nipple discharge", "human mammilla"],
        "Mediastinal": ["Mediastinal Fluid", "Mediastinal Sample", "Mediastinal Aspirate", "Mediastinal Drainage", "Mediastinal Secretion", "mediastinum aspirate", "mediastinum fluid", "mediastinum secretion", "mediastinal swab"],
        "Perineal Swab" : ["Perineal Swab", "Perineum Swab", "Perineal Sample", "Perineal Discharge", "Perineal Fluid Sample", "Perineal Secretion", "perineal abscess", "perineal", "Perineum", "human perineal skin"],
        "Intestinal": ["Intestinal Fluid", "Intestinal Secretion", "Gastrointestinal Fluid", "Intestinal Aspirate", "Intestinal Drainage", "small colon", "small intestine", "appendix", "appendix fluid", "Intestine", "large intestine", "small intestine", "duodenal mucosa", "duodenal aspirate", "caecal effluent", "distal colon", "duodenal mucosa", "duodenum", "ileum", "ileal", "ileal aspirate", "ileal effluent", "ileal fluid", "ileal secretion", "ileal sample", "ileal aspirate", "ileal drainage", "ileal secretion", "ileal sample", "ileostomy effluent", "ileostomy fluid", "ileostomy aspirate", "ileostomy drainage", "ileostomy secretion", "ileostomy sample", "fetal intestine", "human duodenal aspirate", "intestine from healthy child", "intestine of a healthy child", "intestine of adult", "intestine of an adult", "intestine of human", "intestine of patients with crohn's disease", "intestines", "intestines of patients", "jejunal aspirate", "large colon", "small colon", "small bowel", "large bowel", "large intestine", "small intestine", "large intestinal", "small intestinal", "lumen descending colon", "lumen ileum", "lumen transverse colon", "mucosal surface ascending colon", "mucosal surface descending colon", "mucosal surface ileum", "mucosal surface transverse colon", "Sigmoid colon", "Sigmoid colon aspirate", "Sigmoid colon effluent", "Sigmoid colon fluid", "Sigmoid colon secretion", "Sigmoid colon sample", "Sigmoid colon aspirate", "Sigmoid colon drainage", "Sigmoid colon secretion", "Sigmoid colon sample", "Sigmoid Fluid", "Sigmoid aspirate", "Sigmoid secretion", "Sigmoid sample", "Sigmoid aspirate", "Sigmoid drainage", "Sigmoid secretion", "small bowel", "small bowel brush", "terminal ileum", "transcending colon aspirate", "transcending colon brush", "transverse colon", "transverse colon lumen"],
        "Abdominal" : ["Abdominal Cavity Fluid", "Abdominal Cavity Sample", "Abdominal Cavity Aspirate", "Fluid From Abdominal Cavity", "Abdominal Cavity Drainage", "abdominal cavity fluid", "abdominal cavity aspirate", "abdominal cavity sample", "abdominal cavity drainage", "Abdominal Fluid","Abdomen Fluid", "Fluid From Abdomen", "Abdominal Effusion", "abdominal drainage", "abdominal", "abdominal cavity", "abdominal clot", "abdominal drain", "abdominal drain fluid", "abdominal drainage fluid", "punctate from the abdomen", "abd", "abdominal secretion", "abdominal cavity", "drainage abdomen", "fluid; abdomen", "intraabdominal", "intra-abdominal drainage", "abdomen", "abdominal", "abdomen swab", "abdominal swab"],
        "Inguinal Swab": ["Inguinal Swab", "Groin Swab", "Inguinal Region Swab", "Groin Area Swab", "Inguinal Discharge", "groin skin", "inguinal", "groin isolate", "groin", "groin sinus swab", "human groin", "inguinal folds", "inguinal liquid", "inguinal skin", "inguinal region", "inguinal area", "inguinal skin swab", "inguinal skin sample", "human inguinal skin", "human groin swab", "human groin sample", ],
        "Rectovaginal Swab": ["Rectovaginal Swab", "Vaginal-Rectal Swab", "Vaginal-Rectal Sample", "Vaginal-Rectal Culture", "Vaginal-Rectal Isolate","Vaginal-Rectal Swabbing", "vaginal rectal screen", "rectovaginal", "vaginal-rectal swab", "vaginal-rectal", "vaginal/rectal", "anal/vaginal swab"],
        "Lymph": ["Human Lymph", "Lymph", "Lymphatic Fluid", "Lymph Sample", "Lymph Node Aspirate", "Lymph Node Fluid", "Lymphatic Drainage"],
        "Lymph Node" : ["lymph node", "lymph node sample", "lymph node biopsy", "lymph node culture", "lymph node swab", "lymph node tissue", "lymph node fluid", "lymph node aspirate", "lymph node lavage", "lymph node wash", "mesentric lymph node", "submandibular lymph node", "axillary lymph node", "inguinal lymph node", "cervical lymph node", "popliteal lymph node", "bronchial lymph node", "mediastinal lymph node", "thoracic lymph node", "abdominal lymph node", "retroperitoneal lymph node", "pelvic lymph node", "Lymph gland", "mesenteric lymph node biopsies", "Adenoid", "Adenoid tissue", "Adenoid aspirate", "Adenoid culture", "Adenoid biopsy", "Adenoid swab", "Adenoid fluid", "Adenoid sample", "nodi lymphatici mesenterici", "right neck ln aspirate", "left neck ln aspirate", "neck ln aspirate", "neck lymph node aspirate"],
        "Bile": ["Human Bile", "Bile", "Gallbladder Fluid", "Bile Sample", "Biliary Fluid", "Gallbladder Aspirate", "Bile Fluid", "Bile Culture", "gallstone", "fluid gall bladder", "ptbd fluid", "percutaneous transhepatic biliary drainage", "percutaneous transhepatic cholangial drainage", "percutaneous transhepatic cholangial", "PTCD fluid", "PTCD", "PTBD", "PTC", "PTC fluid", "fluid (body)_bile"],
        "Gallbladder" : ["gallbladder", "gall bladder", "gall", "gall bladder drainage", "haematoma gall bladder"],
        "Biopsy": ["Biopsy", "toe Biopsy", " Biopsy", "Biopsy Sample", "soft tissue", "Organ Biopsy", "Liver Biopsy", "superficial Biopsy", "Intestinal Biopsy", "Gut Biopsy", "Bowel Biopsy", "Intestinal Tissue Sample", "Intestinal Tract Biopsy", "toe tissue", "bronchoscopic brushing", "periprosthetic biopsy", "prosthetic biopsy", "Abscess tissue", "sterile tissue", "tissue", "tissues", "liquid tissue", "synovial Biopsy", "Nail bed", "Nail bed tissue", "Nail bed biopsy", "heart valve", "trycuspid valve", "aortic valve", "mitral valve", "valve biopsy", "carious dentine", "carious dentine biopsy", "carious dentine sample", "corneal scrapings", "colon epithelium", "ankle wound Biopsy", "wound biopsy", "debridement", "esophagus brush", "fascia", "fine needle aspirate", "finger flexor tendon sheath", "foot Biopsy", "forearm Biopsy", "glenoid membrane", "gluteal_Biopsy", "groin sinus", "Biopsy, rectum", "heel Biopsy", "hip tissue", "humeral membrane", "glenoid membrane", "shoulder collar membrane", "debridement tissue", "knee Biopsy", "tibia debrid", "oral lichen planus biopsy", "lichen planus biopsy", "pleural Biopsy", "prostate biopsy", "proximal phalynx Biopsy", "distal phalynx biopsy", "pulmonary Biopsy", "right ischial Biopsy", "right ischial tissue", "right hip sinus tract Biopsy", "shin Biopsy", "sphacelus", "transthoracic needle aspiration", "bicep sheath of a left shoulder", "liver puncture"],
        "Graft" : ["graft", "graft Biopsy", "graft sample", "graft site", "graft swab", "devitalized Biopsy", "diaphragm Biopsy", "mitral valve vegetation"],
        "Intraocular Fluid": ["Intraocular Fluid", "Eye Fluid", "Ocular Fluid", "Fluid From Eye", "Eye Aspirate"],
        "Catheter": ["Catheter Fluid", "Catheter Fluid Collection", "Catheter Fluid Sample", "Catheter Fluid Isolate", "Catheter Fluid Specimen", "Catheter liquid", "cvp", "central line", "cvc", "catheter", "catheter tip", "catheter fluid", "catheter drainage", "catheter aspirate", "central venous catheter", "gastric tube", "central vein catheter", "peripheral venous catheter", "peripheral intravenous catheter", "peripheral venous line", "peripheral intravenous line", "peripheral venous access device", "peripheral intravenous access device", "peripheral venous cannula", "peripheral intravenous cannula", "peripheral venous catheterization", "peripheral intravenous catheterization", "peripheral venous access procedure", "peripheral intravenous access procedure", "peripheral venous cannulation", "peripheral intravenous cannulation", "peripheral venous access technique", "peripheral intravenous access technique", "deep-venous catheter", "Nasointestinal Tube", "Nasoenteric tube", "catheter indwelling", "catheter liquid", "protected specimen brushes", "urethra catheter", "gastric intubation", "nasogastric tube", "nasointestinal tube", "nasojejunal tube", "tracheostomy tube", "intubation tube", "cholecystostomy tube fluid", "art line", "choli tube drain", "drain tube tip", "drainge tube tip", "nasointestinal tube", "nasojejunal tube", "nephrostomy tube", "pej", "percutaneous endoscopic gastrostomy tube", "permacath", "suction tip", "transcatheter bc", "uvc tip", "umbilical vein catheter", "umbilical vein catheter tip", "umbilical venous catheter", "umbilical venous catheter tip", "pig-tail ureteric stent", "foley"],
        "Drainage Fluid" : ["Drainage Fluid", "Drainage Sample", "Drainage Aspirate", "Drainage Culture", "Drainage Fluid Collection", "Drainage Fluid Sample", "Drainage Fluid Isolate", "Drain", "IR Drain", "liver puncture fluid", "drainage aspirate", "drainage fluid", "drainage sample", "drainage culture", "drainage isolate", "drainage aspirate sample", "drain fluid", "drain liquid", "drainage fluid", "drainage liquid", "percutaneous drainage", "choli tube drainage fluid", "choli tube drainage", "drain", "drain fluid culture", "drainage", "drainage of fluid", "drainage or drn", "drainage secretion", "drainage jackson pratt", "drainage jackson pratt fluid", "drained fluid", "drained material", "fluid from drain", "gastrostomy tube drainage", "g-tube drainage", "Jackson-Pratt drain", "jp drain", "aspirate, jp drain", "liver drainage", "nephretic drain fluid", "liquid drain", "peg tube drainage", "pelvic drainage", "pelvic drain fluid", "pelvic drain", "pelvic drainage", "percutaneous drainage", "percutaneous drain fluid", "percutaneous drain", "percutaneous drainage fluid", "right chest tube fluid", "right flank drainage", "tubal drain fluid", "tubal drain", "tubal drainage", "tubal drainage fluid", "tubal drainage aspirate", "tubal drainage sample", "tubal drainage isolate", "tubal drainage culture", ],
        "Decubitus Swab" : ["decubitus swab", "pressure ulcer", "decubitus", "Sacrum ulcer", "sacral ulcer", "pressure sore", "Sacral abcess", "Sacrum abcess", "sacrum swab",  "coccygeal lesions", "coccyx lesion", "coccyx abscess", "coccyx wound", "decubitus pressure ulcers", "decubitus ulcer", "ischium wound", "bedsore", "sacral lesion", "sacral sore", "sacral wound", "sacral ulcer", "sacral abscess", "ischial ulcer", "ischial sore", "ischial wound", "ischial abscess", "ischium ulcer", "ischium sore", "ischium wound", "ischium abscess"],
        "Uncategorized Fluid": ["body fluids", "asp fluid", "fluid", "fluids", "fluid sample", "fluid isolate", "fluid culture", "fluid specimen", "fluid aspirate", "effusion", "fluid aspirate", "fluid collection", "fluid drainage", "fluid wash", "fluid swab", "fluid lavage", "fluid biopsy", "fluid puncture", "secretion", "aspirate", "bodily fluid", "body fluid", "lavage fluid", "secreta", "diversion fluid", "Exudate", "punctate", "puncture fluid", "secretion sample", "biological fluid", "breast fluid", "deep fluid", "enema fluid", "esophagus aspirate", "extravasate fluid", "fasciotomy fluid", "fluid from conjunctiva", "gluteal fluid", "hip fluid", "hip secretion", "human fluid", "human body fluid", "human effusion", "human exudate", "human secretion", "humeral canal fluid", "infective secretion", "knee secretion", "lavage", "lavage fluid", "leg secretion", "liquid", "liver fluid isolate", "miscellaneous body fluid", "osteomyelitis aspirate", "pancreatic secretions", "other fluid", "other body fluid", "other bodily fluid", "other discharge", "other secretion", "other exudate", "other aspirate", "pancreatic fluid", "paravertebral fluid", "patient hydatid fluid", "hydatid fluid", "pelvic collection aspirate", "pelvic secretion", "pelvic fl", "perihepatic fluid", "periprosthetic liquid", "puncture solution", "puncturesolution", "puncturinee fluid", "secretion surgical", "punctuate liquid", "diversion fluid", "endosecret", "human, fluids", "shunt fluid", "surg aspir", "surgical puncture fluid", "traumatic discharge", "urinogenital fluid", "withdrawing fluid", "aspirate_other", "aspiration", "biological secretion", "bodyfluid", "debridement fluid", "hip secretion", "sterile body fluid"],
        "Cyst" : ["Cyst Fluid", "Cyst Sample", "Cyst Aspirate", "Cyst Culture", "thyroglossal cyst", "halscysta", "cystic fluid", "cystic aspirate", "cystic sample", "cystic culture", "cystic drainage", "cystic wash", "cystic swab", "cystic biopsy", "cystic puncture", "cystic tap", "liver cyst", "ovarian cyst", "postoperative maxillary cyst", "maxillary cyst", "right inner gluteal pilonidal", "gluteal cyst", "gluteal skin cyst", "Cyst Swab", "paraprostatic cyst"],
        "Uterine": ["Uterine Fluid", "Uterine Sample", "Endometrial Sample", "Uterine Swab", "Uterine Aspirate", "Uterine Lavage", "Uterine Culture", "Uterine Biopsy", "Endometrial Fluid", "Endometrial Swab", "Endometrial Aspirate", "Endometrial Lavage", "Endometrial Culture", "Endometrial Biopsy", "Uterine Wash", "Uterine Exudate", "Uterine Discharge", "uterine rinsing", "uterus", "endometrium", "pyometra", "equine infectious endometritis", "swab from uterus", "uterus discharge", "uterine aspirate", "uterine aspirate sample", "uterine aspirate fluid", "uterine aspirate isolate", "uterine aspirate culture", "uterine swab isolate", "uterine swab sample", "uterine swab aspirate", "uterine swab drainage", "uterine swab secretion", "uterine swab isolate", "uterine swab culture", "endometrial aspirate sample", "endometrial aspirate fluid", "endometrial aspirate isolate", "endometrial aspirate culture", "endometrial swab isolate", "endometrial swab sample", "endometrial swab aspirate", "endometrial swab drainage", "endometrial swab secretion"],
        "IUDs": ["IUD", "IUD swab", "intrauterine device", "intrauterine device swab", "intrauterine device sample", "IUD wash"],
        # "Coelomic Cavity" : ["coelomic cavity", "coelomic fluid", "coelomic sample", "coelomic aspirate", "coelomic lavage", "coelomic wash", "coelomic swab", "coelomic culture", "coelomic biopsy", "coelomic puncture", "coelomic tap"],
        # "Unknown": ["General Sample", "Unspecified Sample", "Unclassified Sample", "Other Sample", "Miscellaneous Sample", "Non-Specific Sample", "Missing", "Unknown", "Na", "Missing", "Not Provided", "None", "Nan", "No Data", "Not Applicable", "Not Recorded", "Not Available","Not Found", "Unspecified", "Not Specified", "Not Collected", "Not Determined", "Clinical Sample", "Pathological Sample"],
        "Metagenomic" : ['metagenomic', 'metagenomic sample', 'microbiome', 'microbiome sample', 'microbial sample', 'microbial community', 'microbial community sample', 'microbial metagenome', 'microbial metagenome sample', 'microbial community', 'microbiome', 'microbiome metagenome sample', 'microbiome community', 'microbiome community sample', 'microbiome community metagenome', 'microbiome', 'oral microbiome', 'resident flora'],
        "Cecal" : ['Cecal', "caecum", "ceca", "cecal content", "caecal effluent", "caecal swabs"],
        "Bone" : ['sacrum','Bone isolate', 'bone', 'bone sample', 'bone marrow', 'bone marrow sample', 'bone biopsy', 'bone tissue', 'bone tissue sample', 'bone marrow biopsy', 'bone marrow tissue', 'bone marrow tissue sample', 'bone marrow biopsy', "Patient Bone Marrow","Human Bone Marrow", "Bone Marrow", "Bone Marrow Aspirate", "Bone Marrow Biopsy", "Hematopoietic Tissue", "Marrow Sample", "Hematopoietic Fluid", "Bone Marrow Culture", "marrow", "strenum", "fibia", "carpals", "tubula", "bone fragment", "metatarsal", "tarsal", "sarcum", "bone aspiration", "bone marrow aspirate", "biopsy of vertebra", "vertebral biopsy", "calcaneous", "distal fibula", "fibula", "distal phalanx", "proximal phalanx", "distal femur", "femur", "proximal femur", "distal radius", "radius", "proximal radius", "distal ulna", "ulna", "proximal ulna", "distal tibia", "tibia", "proximal tibia", "distal fibula", "fibula", "proximal fibula", "ischium", "left ischium", "right ischium", "ischial tuberosity", "left lateral malleolus", "right lateral malleolus", "lateral malleolus", "left medial malleolus", "right medial malleolus", "medial malleolus", "left metacarpal", "right metacarpal", "metacarpal", "left metatarsal", "right metatarsal", "metatarsal", "left phalanx", "right phalanx", "phalanx", "left proximal phalanx", "right proximal phalanx", "proximal phalanx", "manubrium", "manubrium sterni", "sternum", "sternal", "sternal biopsy", "olecranon", "patella", "patellar", "patellar biopsy", "pelvis", "pelvic bone", "pubis", "pubic bone", "sacrum", "sacral biopsy", "scapula", "scapular", "scapular biopsy", "spine", "spinal biopsy", "sternal biopsy", "tarsal bone", "tibia", "tibial", "right scapula", "left scapula", "scapula", "right femur", "left femur", "femur", "right humerus", "left humerus", "humerus", "right radius", "left radius", "radius", "right ulna", "left ulna", "ulna", "right tibia", "left tibia", "tibia", "right fibula", "left fibula", "fibula", "right ischium", "left ischium", "ischium", "rt ischium", "lt ischium", "sacral", "sacral isolate", "sacral tissue", "sacro/coccyx", "shoulder", "talus", "tarsal", "tibia dist", "tibia distal", "tibia midshaft", "tibia prox", "tibie", "spinous process"],
        "Implant" : ["screws", "implant", "infected implant", "bone screw", "femoral pin", "femoral plate", "implant site", "implant swab", "surgical implant",  "bar", "rod", "prosthetic", "artificial joint", "artificial limb", "artificial device", "artificial implant", "artificial prosthesis", "artificial replacement", "artificial organ", "artificial body part", "artificial body implant", "artificial body prosthesis", "artificial body replacement", "prosthetic implant hip", "prosthetic implant knee", "prosthetic implant joint", "prosthetic implant limb", "prosthetic implant device", "prosthetic implant organ", "prosthetic implant body part", "prosthetic implant body prosthesis", "prosthetic implant body replacement", "prosthetic implant body part replacement", "prosthetic implant body part prosthesis", "rod spine", "screw", "screw spine", "glenoid explant", "explant", "humeral explant", "humeral head explant", "implant sonicate", "knee prothesis", "vascular prosthesis", "prosthesis", "prosthesis liquid", "prosthetic heart valve", "prosthetic joint right knee", "prosthetic material", "protected distal sampling", "sonication prosthesis", "total knee hardware", "vascular prosthesis", "acetabular cup", "back, ortho screw", "joint prosthesis", "sonicate fluid"],
        "Pericardial": ["Pericardial Fluid", "Pericardial Effusion", "Heart Fluid", "Pericardial Sample", "Fluid In Pericardium", "Pericardial Aspirate", "Pericardial Drainage", "pericardial effusion", "pericardial aspirate", "pericardium", "pericardial cavity", "fluid, pericardial", "pericardial"],
        "Dialysis Fluid": ["Dialysis Fluid", "Dialysis Sample", "dialysate", "peritoneal dialysis fluid"],
        "Autopsy" : ["autopsy", "autopsy sample", "post morterm", "post mortem sample", "Post Morterm Swab"],
        "Fistula" : ["Fistula","Fistula Sample", "Fistula Swab", "Fistula Culture", "Fistula tract swab", "fistula tract", "fistula tract swab", "fistula tract culture", "fistula tract sample", "fistula tract aspirate", "fistula tract biopsy", "fistula tract fluid", "fistula tract discharge", "fistula tract drainage", "fistula tract exudate", "fistula tract pus", "fistula tract secretion", "oroantral fistula", "oroantral fistula swab", "anal fistula", "perianal fistula", "fistula fluid"],
        "Stones" : ["Stones","Calculi", "Urinary Stones", "Kidney Stones", "Gallstones", "Bladder Stones", "Stones Sample", "Stones", "biliary stones", "bile stone", "kidney matrix stone", "calculus", "renal calculi", "cystolith", "bladder cystolith", "ureteral stone", "urinary calculus", "urinary calculi", "Urolith"],
        "Bursa Sample": ["Bursa fluid", "bursa", "bursa content", "right hip bursa Biopsy", "left hip bursa biopsy", "hip bursa biopsy", "bursa biopsy"],
        "Tooth" : ["tooth", "teeth", "milk tooth",],
        "Pacemaker" : ["pacemaker", "pacemaker lead", "pacemaker wire", "pacemaker site", "pacemaker swab", "pacemaker culture", "pacemaker aspirate", "pacemaker sample", "pacemaker fluid", "pacemaker drainage", "pacemaker exudate", ],
        "Prostatic": ["Prostate", "prostatic fluid", "prosate fluid",],
        "Fetus": ["Fetus", "Fetal Sample", "Fetal Fluid", "Fetal Aspirate", "Fetal Culture", "Fetal Drainage", "fetal tissue", "fetal blood", "fetal swab", "fetal biopsy", "fetal aspirate", "fetal fluid sample", "fetal fluid culture", "fetal fluid drainage", "fetal fluid aspirate", "fetal fluid biopsy"],
        "Scalp Swab": ["Scalp Swab", "Scalp Sample", "Scalp Culture", "Scalp Aspirate", "Scalp Fluid", "Scalp Drainage", "scalp skin", ],

    }

    logging.debug(f"Starting the search for sample source for human samples.")
    
    best_match_score = 0
    best_sample = "Unknown"

    try:
        if isinstance(fuzzy_str, str):
            logging.debug(f"Entered the str loop in find_human_sample function with string: `{fuzzy_str}`.")

            part_digit_free = re.sub(r'\d+', '', fuzzy_str)
            part = re.sub(r'[\[\]{}()]', '', part_digit_free)
            part = re.sub(r'"', '', part)
            part = re.sub(r'\'', '', part)
            logging.debug(f"Normalized string: {part}")

            
            best_match_score = 0
            best_sample = "Unknown"
            best_matched_string = ""
            normalized_str_matched = ""

            normalized_str = part.lower()  
            logging.debug(f"Checking full string: {normalized_str}")

            for sample, sample_fuzzy_words in human_body_samples.items():
                for sample_fuzzy_word in sample_fuzzy_words:
                    similarity_score = fuzz.ratio(normalized_str, sample_fuzzy_word.lower())
                    if similarity_score > best_match_score:
                        best_match_score = similarity_score
                        best_sample = sample
                        best_matched_string = sample_fuzzy_word
                        normalized_str_matched = normalized_str

            
            if best_match_score >= 85:
                logging.debug(f"Matched without splitting: {best_sample} with {best_match_score} score. The best matched string was `{best_matched_string}` with `{normalized_str_matched}`.")
                return best_sample, best_match_score

            
            logging.debug("No match found in unsplit form. Trying with splitting strings.")

            parts = re.split(r'[-+:;/_,\s]|\band\b', part)
            parts = [part.strip() for part in parts if part.strip()] 
            logging.debug(f"Split string parts: {parts}")

            for normalized_str in parts:
                normalized_str = normalized_str.lower()
                logging.debug(f"Checking part: {normalized_str}")
                for sample, sample_fuzzy_words in human_body_samples.items():
                    for sample_fuzzy_word in sample_fuzzy_words:
                        similarity_score = fuzz.ratio(normalized_str, sample_fuzzy_word.lower())
                        if similarity_score > best_match_score:
                            best_match_score = similarity_score
                            best_sample = sample
                            best_matched_string = sample_fuzzy_word
                            normalized_str_matched = normalized_str

            if best_match_score >= 85:
                logging.debug(f"Matched after splitting: {best_sample} with {best_match_score} score.")
                return best_sample, best_match_score
            else:
                logging.debug(f"No match found. Returning `Unknown`.")
                return "Unknown", 0

        else:
            logging.debug(f"Given input is not a string. Returning `Unknown`.")
            return "Unknown", 0

    except Exception as e:
        logging.error(f"An error occurred while finding human sample: {e}")
        return "Unknown", 0

def find_disease_from_sample(identified_sample: str):
    if identified_sample in ['Blood', 'Plasma', 'Serum']:
        best_disease = 'Bacteremia'
    elif identified_sample in ['Urine', 'Urethral Swab', 'Bladder']:
        best_disease = 'UTIs'
    elif identified_sample in ['Sputum', 'Pleural Fluid', 'Bronchoalveolar Lavage', 'Tracheal', 'Nasal', 'Throat', 'Laryngeal']:
        best_disease = 'Respiratory Infections/Illnesses'
    elif identified_sample in ['Pus', 'Wound Swab', 'Perianal Swab', 'Skin Swab', 'Axillary Swab', 'Inguinal Swab', 'Skin', 'Decubitus Swab', 'Perineal Swab', 'Mammilla Swab']:
        best_disease = 'Soft Tissue Infections/Colonization'
    elif identified_sample in ['Bile']:
        best_disease = 'Gallbladder and Biliary Tract Disorders'
    elif identified_sample in ['Stool', 'Peritoneal Fluid', 'Rectal Swab', 'Intestinal', 'Gastric', 'Abdominal']:
        best_disease = 'Gastrointestinal Disorders'
    elif identified_sample in ['CSF']:
        best_disease = 'CNS infections'
    elif identified_sample in ['Bone', 'Synovial Fluid', 'Cartilage']:
        best_disease = 'Bone Infection/Inflammation'
    elif identified_sample in ['Eye Swab', 'Tears', 'Vitreous Humor',]:
        best_disease = 'Eye Infection'
    elif identified_sample in ['Ear Swab']:
        best_disease = 'Ear Infection'
    elif identified_sample in ['Oral Swab', 'Saliva', 'Root Canal', 'Tooth']:
        best_disease = 'Oral Infection/Colonization'
    elif identified_sample in ['Cervical/Vaginal Swab', 'IUDs', 'Genital Swab', 'Penile Swab', 'Prostatic', 'Semen', 'Scrotal/Testicular', 'Vulval Swab']:
        best_disease = 'Reproductive System Infections/Illnesses'
    elif identified_sample in ['Bursa Sample']:
        best_disease = 'Bursitis'
    elif identified_sample in ['Lymph Node', 'Lymph']:
        best_disease = 'Lymph Node Infections'
    elif identified_sample in ['Pericardial']:
        best_disease = 'Cardiac Disorders'
    elif identified_sample in ['Implant']:
        best_disease = 'Implant-associated Infections'
    else:
        best_disease = 'Unsorted Bacterial infections'
    
    return best_disease

def find_hospital_equipment_category_and_name(fuzzy_str):
    logging.debug(f"Started! finding the hospital equiment for string `{fuzzy_str}`")
    hospital_equipments = {
        "Unknown": {
            "Unknown": ["environment", "env", "environmental swab", "environmental sample", "environmental isolate",],
        },
        "Beds and Support Equipment": {
            "Bed": ["Bed", "hospital bed", "adjustable bed", "electric hospital bed", "manual hospital bed", "patient bed", "bedsheet", "bed bar", "bed side rail", "bed rails", "bed controllers", "bed side switch", "bedside button", "sickbed swab"],
            "Bedside Table": ["Bedside Table", "hospital bedside table", "patient bedside table", "medical bedside table", "bed side table", "bed table", "overbed table", "over bed table", "hospital overbed table", "patient overbed table",],
            "Care Cart": ["Care Cart", "nurse care cart", "medical cart", "hospital trolley", "medication cart", "equipment trolley", "crash cart", "emergency cart", "crash trolley", "dressing trolley"],
            "Wheelchair": ["Wheelchair", "hospital wheelchair", "manual wheelchair", "electric wheelchair", "patient wheelchair", "mobility chair"],
            "Hospital Stretcher": ["Hospital Stretcher", "medical stretcher", "emergency stretcher", "patient stretcher", "transport stretcher"],
            "Examination Table": ["Examination Table", "doctor's table", "patient exam table", "adjustable examination table", "clinic exam table",  "exam table"],
            "Treatment Furniture": ["Treatment Furniture", "hospital treatment table", "medical treatment table", "patient treatment table", "treatment chair", "treatment couch", "treatment bed", "treatment table"],
            "Walker" : ["walker", "hospital walker"],
            "Incubator": ["Incubator", "intensive care unit incubator", "premature baby incubator", "newborn care incubator", "neonatal incubator", "Incubator swab", "Incubator sample", "infant incubator"],
            "Gurney": ["Gurney", "hospital gurney", "emergency gurney", "transport gurney", "patient transport gurney"]
        },
        "Infection Control and Hygiene": {
            "Soap Dispenser": ["soap dispenser", "hand soap dispenser", "automatic soap dispenser", "liquid soap", "soap station", "soap bottle", "liquid soap dispenser",],
            "Disinfectant Dispenser": ["disinfectant dispenser", "disinfectant station", "sanitizer dispenser", "disinfectant bottle", "sanitizing station", "alcohol dispenser",],
            "Disinfectant": ["disinfectant", "hospital disinfectant", "antiseptic", "sanitizing agent", "sterilizing solution", "hand sanitizer", "alcohol", "hospital alcohol", "alcohol foam"],
            "Sanitizing Wipes": ["sanitizing wipes", "cleaning wipes", "antibacterial wipes", "hospital wipes", "surface sanitizing wipes"],
            "Sterilizing Station": ["sterilizing station", "disinfection station", "medical sterilizing area", "sanitization station"],
            "Gloves": ["gloves", "medical gloves", "disposable gloves", "latex gloves", "non-latex gloves", "sterile gloves"],
            "Face Mask/Respirator": ["face mask", "medical face mask", "surgical mask", "disposable face mask", "N95 mask", "respirator mask", "respirator"],
            "Cleaning equipment": ["cleaning equipment", "hospital cleaning tools", "Mopper", "Electric mop", "Swiffer", "mop handle", "Cleaning trolly", "Cleaning trolly", "Cleaning cart", "Cleaning bucket", "Cleaning mop", "Cleaning mop head", "Cleaning mop handle", "Cleaning mop bucket", "Brush", "mop handle"],
        },
        "Monitoring and Diagnostic Equipment": {
            "ECG device": ["ecg device", "electrocardiogram device", "ecg machine", "electrocardiograph", "heart monitor", "cardiac monitor", "ECG clip", "Electrocardiogram clip", "ecg monitor"],
            "Vital Monitor": ["vital monitor device", "vital monitoring device", "vital monitor button", "vitals monitor screen", "vital monitor display", "vital monitor button", "vital monitor knob"],
            "Pulse oximeter": ["pulse oximeter", "blood oxygen monitor", "oxygen saturation monitor", "pulse ox", "oximeter probe", "pulse oximetry device"],
            "Blood Pressure cuff": ["blood pressure cuff", "sphygmomanometer", "manual blood pressure cuff", "automatic blood pressure cuff", "bp cuff"],
            "Patient Monitor": ["patient monitor", "heart rate monitor"],
            "Syringe": ["syringe", "syringe pump", "syringe driver", "insulin syringe", "syringe for injections", "sterile syringe"],
            "Canula" : ["canula", "canula sample", "canuala cap", "luer lock plug", "canula wings", 'cannula', 'secrete from cannula'],
            "Stethoscope": ["stethoscope", "manual stethoscope", "digital stethoscope", "electronic stethoscope", "doctor stethoscope"],
            "Thermometer": ["thermometer", "digital thermometer", "oral thermometer", "infrared thermometer", "temperature probe", "thermometer gun", "rectal thermometer"],
            "Glucometer": ["glucometer", "blood glucose meter", "diabetes meter", "glucose monitor", "diabetic glucose test", "blood sugar meter"],
            "Diagnostic Kit": ["diagnostic kit", "medical diagnostic kit", "rapid diagnostic test", "point-of-care diagnostic kit", "diagnostic test kit", "diagnostic assay kit"],
        },
        "Surgical and Therapeutic Equipment": {
            "Defibrillator": ["defibrillator", "automated external defibrillator", "aed", "manual defibrillator", "external defibrillator", "cardiac defibrillator", "AED pads", "defibrillator electrodes", "cardiac defibrillator pads", "defibrillator patches", "external defibrillator pads"],
            "Endotracheal Tube": ["endotracheal tube", "et tube", "intubation tube", "artificial airway", "tracheal tube", "endotracheal intubation"],
            "Endobronchial Tube": ["endobronchial tube", "bronchial tube", "bronchial intubation", "endobronchial intubation", "bronchial intubation tube", "endobronchial intubation tube"],
            "Ventilator": ["ventilator", "mechanical ventilator", "respiratory ventilator", "breathing machine", "oxygen concentrator", "ventilator machine", "hospital ventilator", "ventilator shelf"],
            "Nebulizer": ["nebulizer", "nebulizer machine", "compressor nebulizer", "breathing treatment device", "inhalation therapy device", "respiratory nebulizer", "Nebuliser", "nebuliser machine", "compressor nebuliser", "respiratory nebuliser"],
            "Anesthesia Machine": ["anesthesia machine", "anaesthesia workstation", "anesthesia delivery system", "gas anesthesia machine", "intubation anesthesia machine"],
            "Surgical Instruments": ["surgical instruments", "scalpel", "surgical scissors", "surgical forceps", "needle holder", "surgical clamps"],
            "Laser Therapy Device": ["laser therapy device", "medical laser", "laser therapy machine", "diode laser", "surgical laser device", "laser wound healing machine"],
            "Sterile Field": ["sterile field", "sterile drape", "sterile surgical field", "surgical sterile area", "sterilized surgical area", "sterile site"],
            "Electrosurgical Unit": ["electrosurgical unit", "electrocautery unit", "diathermy machine", "electrosurgical pencil", "electrocoagulation unit", "electrosurgical system"],
            "Tracheostomy Tube": ["tracheostomy tube", "trach tube", "tracheal tube", "tracheostomy tube holder", "tracheostomy tube cuff", "tracheostomy tube connector", "trachy", "tracheoflex"],
            "Driveline" : ["driveline", "ventricular assist device driveline", "left ventricular assist device driveline", "right ventricular assist device driveline", "VAD driveline", "LVAD driveline", "RVAD driveline"],
            "Ventricular Assist Device" : ["ventricular assist device", "left ventricular assist device", "right ventricular assist device", "VAD", "LVAD", "RVAD"],
            "Support Equipments": ["right pelvis external fixator", "pelvis external fixator"],
            "Suture": ["Suture", "stiches", "surgical suture", "suture material"],
        },
        "Medical Imaging Equipment": {
            "X-ray machine": ["x-ray machine", "x-ray equipment", "x-ray imaging device"],
            "Ultrasound machine": ["ultrasound machine", "ultrasonography", "portable ultrasound", "diagnostic ultrasound", "ultrasound scanner", "ultrasound imaging system", "ultrasound gel"],
            "CT scanner": ["ct scanner", "computed tomography", "ct scan", "CAT scan", "CT imaging system", "3D CT scanner"],
            "MRI machine": ["mri machine", "magnetic resonance imaging", "mri scanner", "MRI", "MRI system", "magnetic resonance scanner"],
            "Endoscope": ["endoscope", "medical endoscope", "fiber-optic endoscope", "gastroscope", "bronchoscope", "colonoscope", "duodenoscope", "ureteroscope"],
            "Mammography machine": ["mammography machine", "breast x-ray machine", "mammogram machine", "digital mammography", "mammography system"],
            "Pet scanner": ["pet scanner", "positron emission tomography scanner", "PET imaging machine", "PET scanner", "nuclear medicine scanner"]
        },
        "Infusion and IV Equipment": {
            "Intravenous pole": ["intravenous pole", "iv pole", "iv stand", "intravenous drip stand", "drip pole", "iv support stand", "intravenous infusion pole", "Infusion stand", "drip stand"],
            "IV Infusion set": ["iv infusion set", "intravenous infusion set", "IV drip set", "infusion tubing", "IV fluid administration set", "infusion line"],
            "Catheter": ["catheter", "urinary catheter", "arterial catheter", "foley catheter", "iv catheter", "peripheral catheter", "central line catheter", "cvc" ,"central venous catheter", "Foley Catheter", "catheter container", "catheter bag", "nasogastric tube", "NG tube", "nasal tube", "nasal catheter", "crbsi", "catheter associated bloodstream infection", "Central Venous Access Devices", "CVADs", "Peripherally Inserted Central Catheters", "PICC", "Central Venous Catheters", "CVCs", "Tunneled Central Venous Catheters", "shaldon catheter", "swab from catheter", "catheter isolate", "catheter swab", "gastric tube", "cardiovascular puncture (cvp) tip", "central line", "central venous line", "choli tube drain", "folyes tip", "hick brown lumen", "hickman red lumen", 'IV line', "line tip", "nasointestinal tube", "nasojejunal tube", "nephrostomy tube", "pej", "percutaneous endoscopic gastrostomy tube", "permacath", "siphon", "suction tip", "transcatheter bc", "tube tip", "uvc tip", "umbilical vein catheter", "umbilical vein catheter tip", "umbilical venous catheter", "umbilical venous catheter tip", ],
            "Infusion Pump": ["infusion pump", "iv infusion pump", "drug infusion pump", "medication infusion pump", "volume-controlled infusion pump"],
            "Urine Bag": ["urine bag", "urinary collection bag", "urine collection bag", "urine drainage bag", "urine container", ],
            "Blood bag": ["blood bag", "blood transfusion bag", "blood collection bag", "iv blood bag", "donated blood bag", "transfusion bag"],
            "Perfusion Fluid": ["perfusion fluid", "perfusion solution", "perfusion liquid", "perfusion medium", "perfusion buffer", "perfusion culture", "perfusion sample",],
            "Platlet Concentrate": ["platlet concentrate", "platelet-rich plasma", "platelet transfusion", "platelet preparation", "platelet product", "platelet suspension", "platelet unit"],
            "RBC Bag" : ["red cell bag", "red blood cell bag", "RBC bag"],
        },
        "Sterilization and Disinfection Equipment": {
            "Autoclave": ["autoclave", "steam autoclave", "pressure cooker sterilizer", "autoclave water"],
            "UV sterilizer": ["uv sterilizer", "ultraviolet sterilizer", "uv sanitizing machine", "disinfection UV device", "UV-C sterilizer", "ultraviolet light sterilizer"], 
            "Instrument Washers": ["dental instruments washing machine", ],
            
        },
        "Miscellaneous": {
            "Refrigerator": ["refrigerator", "medical refrigerator", "pharmaceutical refrigerator", "vaccine refrigerator", "cold storage", "medicine fridge", "refrigator handle"],
            "Phone" : ["phone", "dial device", "landline phone", "intrahospital phone line", "call button", "room telephone", "hospital phone", "hospital telephone", "hospital intercom", "hospital intercom system", "hospital intercom device", "hospital intercom unit", "hospital intercom phone", "hospital intercom line", "hospital intercom connection", "hospital intercom network"],
            "Other Instruments": ["Computers", "Computer screens", "TV", "Television", "Remote control", "Keyboard", "Mouse", "Computer mouse", "Computer keyboard", "Computer monitor", "Computer CPU", "Computer tower", "Computer screen", "Computer display", "device infections", "medical equipment", "medical device", "portable suction machine"],
            "Hospital surface": ["Switch button","hospital material","hospital","hospital surface", "hospital floor", "hospital counter", "hospital furniture", "hospital envrionment", "Intesive care units", "icus", "nicu", "hospital drawer", "door", "sliding door", "window", "window door", "door handle", "tap", "water tap", "sink", "sink inside", "outside sink", "sink drain", "Washroom sink", "hospital icu", "hospital intensive care unit", "hospital sink drainage", "sink drainage", "hospital drainage water", "Hospital basin", "basin inside", "outside basin", "washroom basin", "basin drainage", 'sink countertop', 'hospital washbasin', 'washbasin', 'basin', 'hospital surveillance', 'hospital surface surveillance', 'hospital surface swab', 'hospital surface sample', 'hospital surface culture', 'hospital surface isolate', 'hospital surface cleaning', 'hospital surface disinfection', 'hospital surface sterilization', 'hospital surface wipe', 'hospital surface wipe sample', 'hospital surface wipe culture', 'surface surveillance', 'surface swab', 'surface sample', 'surface culture', 'surface isolate', 'surface cleaning', 'surface disinfection', 'surface sterilization', 'surface wipe', 'surface wipe sample', 'surface wipe culture', 'inner wall of basin nurse station', "floor trap bathroom room swab", "bathroom swab", "bathroom room swab", "patient room", "rm drain", "rm sink", "rm floor trap", "rm floor trap swab", "rm floor trap sample", "rm floor trap culture", "rm floor trap isolate", "rm floor trap cleaning", "rm floor trap disinfection", "rm floor trap sterilization", "rm floor trap wipe", "rm floor trap wipe sample", "rm floor trap wipe culture", "hospital room drain", "hospital room sink", "hospital room floor trap", "toilet", "hospital toilet", "toilet handle", "hospital toilet handle", "toilet seat", "hospital toilet seat", "toilet bowl", "hospital toilet bowl", "toilet paper holder", "hospital toilet paper holder", "toilet paper dispenser", "hospital toilet paper dispenser", "toilet tank", "hospital toilet tank", "toilet floor", "hospital toilet floor", "toilet wall", "hospital toilet wall", "toilet door", "hospital toilet door", "Hospital chair"],            
            "Infant incubator": ["infant incubator", "neonatal incubator", "baby incubator", "newborn incubator"],
            "Medical cooler": ["medical cooler", "medicine cooler", "cooling box", "drug cooler", "temperature-controlled medical cooler"],
            "Oxygen concentrator": ["oxygen concentrator", "portable oxygen concentrator", "oxygen generation unit", "medical oxygen machine"],
            "Vaccum Extrcator": ["vaccum extractor", "vacuum extraction device", "vacuum delivery system", "vacuum-assisted delivery device", "vacuum-assisted birth device"],
            "Hospital Food" : ["Hospital eggs", "hospital food", "hospital meal", "hospital catering", "hospital kitchen", "hospital cafeteria", "hospital dining", ],
            "Hospital Water": ["hospital water", "hospital drinking water", "hospital tap water", "hospital water supply", "hospital water source", "hospital water system", "hospital water treatment", "hospital water filtration", ],
            "Contaminated Medicine": ["hospital medicine", "hospital drug", "hospital medication", "artificial tears", "medicine bioflor", "medicine capsule", "clade", "geriatric medicine", "orine"],
            "Miscellaneous Fluids" : ["hospital fluid", "hospital liquid", "hospital solution", "sterile fluid", "sterile liquid", "sterile solution", "hospital sterile fluid", "hospital sterile liquid", "hospital sterile solution", "sterile water", "sterile saline", "sterile solution", "sterile saline solution", "sterile water for injection", "sterile saline for injection", "Transplant transfer fluid", "transport medium", "preservation fluid", "preservation solution", "transportation fluid", "transportation solution", "transportation medium", "transportation liquid", "transportation water", "transportation saline", "transportation buffer", "transportation gel", "transportation agent", "washing liquid"],
            "Medical Waste": ["Medical Waste", "Hospital Waste", "Biomedical Waste", "Clinical Waste", "Infectious Waste", "Sharps", "Used Syringes", "Pharmaceutical Waste", "Expired Medications", "Discarded Medical Instruments", "Toxic Healthcare Waste", "Patient Waste", "Used Bandages", "Used Needles", "Surgical Waste"],
            "Hospital Sewage" : ["Hospital sewage water", "hosptial waste water", "hospital wastewater"],
            "Comtaminated Supplements" : ["probiotic supplement", "dietary supplement", "nutritional supplement", "vitamin supplement", "herbal supplement", "nutraceutical", "probiotic capsule", "probiotic powder", "probiotic drink", "probiotic yogurt", "probiotic beverage"],
            "Intensive Care Unit" : ["Intensive Care Unit", "ICU", "intensive care unit", "critical care unit", "intensive therapy unit", "critical care area", "intensive care area", "hospital icu", "hospital icu ward", "hospital icu room", "medical icu", "surgical icu", "neuro icu", "cardiac icu", "pediatric icu", "adult icu", "icu room", "icu ward", "icu area", "micu"],
            "Neonatal Intensive Care Unit" : ["Neonatal Intensive Care Unit", "NICU", "neonatal intensive care unit", "neonatal therapy unit", "neonatal critical care unit", "neonatal care area"],
            "Live Vaccines": ["vivotif", "Varivax", "RotaTeq", "Rotarix", "Sabin vaccine", "YF-VAX", "Stamaril", "FluMist", "BCG Vaccine", "Dengvaxia", "Qdenga", "Zostavax"],
            "OT Room" : ["OT Room", "Operation Theatre", "Surgical Theatre", "Surgical Suite", "Operating Room", "Surgical Room", "Surgical Unit", "Operating Theatre"],
            "Surgical Ward" : ["Surgical Ward", "Surgical Unit", "Surgical Floor", "Surgical Department", "Surgical Care Unit", "Surgical Area", "Surgical Section"],
            "Newborn Nursery" : ["Newborn Nursery", "Neonatal Nursery", "Infant Nursery", "Baby Nursery", "Newborn Care Unit", "Neonatal Care Unit", "Infant Care Unit", "Newborn Unit"],
        },
        "Specialized Medical Equipment": {
            "Dialysis Machine": ["dialysis machine", "hemodialysis machine", "kidney dialysis machine", "dialysis equipment", "renal dialysis unit", "Dialysate", "capd", "peritoneal dialysis", "peritoneal dialysis machine", "peritoneal dialysis unit", "peritoneal dialysis equipment", "Continuous Ambulatory Peritoneal Dialysis", "CAPD", "Automated Peritoneal Dialysis", "APD", "dialysis liquid", "dialysis soution", "haemodialysis unit monfalcone"],
            "cryotherapy machine": ["cryotherapy machine", "cryosurgery unit", "cryotherapy device", "cold therapy machine", "cryotherapy equipment"],
            "Urostomy Device": ["urostomy device", "urostomy bag", "urostomy pouch", "urostomy collection bag", "urostomy drainage bag", "urostomy container"],
        },
        "Emergency and Trauma Equipment": {
            "Trauma kit": ["trauma kit", "emergency trauma kit", "first aid kit", "trauma response kit", "medical emergency kit"],
            "Splint": ["splint", "splinting device", "leg splint", "arm splint", "trauma splint", "splint bandage"],
            "Tourniquet": ["tourniquet", "medical tourniquet", "blood pressure tourniquet", "hemostatic tourniquet", "trauma tourniquet"],
            "Ambulance stretcher": ["ambulance stretcher", "ambulance transport stretcher", "field stretcher", "mobile stretcher", "emergency stretcher"],
            "Resuscitation bag": ["resuscitation bag", "bag-valve-mask", "ambu bag", "manual resuscitator", "self-expanding resuscitator",],
            "Emergency Room" : ["Emergency Room", "trauma room", "trauma bay", "emergency department", "emergency unit", "emergency care area", "emergency treatment room"],
            "Transport Bag" : ["Transport Bag", "Medical Transport Bag", "Emergency Transport Bag", "Patient Transport Bag", "Medical Emergency Bag", "First Aid Transport Bag", "Trauma Transport Bag", "Portable Medical Bag"],
        }
    }

    best_match_score = 0
    best_category = "Unknown"
    best_equipment = "Unknown"

    try:
        if isinstance(fuzzy_str, str):
            logging.debug(f"Entered the str loop in find_hospital_equipment function with string: `{fuzzy_str}`.")
            
            part_digit_free = re.sub(r'\d+', '', fuzzy_str)
            part = re.sub(r'[\[\]{}()]', '', part_digit_free)
            part = re.sub(r'"', '', part)
            part = re.sub(r'\'', '', part)
            logging.debug(f"Normalized string: {part}")
            space_count = part.count(' ')
            if space_count > 2:
                parts = re.split(r'[-+:;/_,\s]|\band\b', part)
                parts = [part.strip() for part in parts]
            else:
                parts = [part]
            logging.debug(f"Split string: {parts}")

            for normalized_str in parts:
                normalized_str = normalized_str.lower()

                for equipment_category, equipment_dictionaries in hospital_equipments.items():
                    for equipment_name, fuzzy_equipment_names in equipment_dictionaries.items():
                        for synonym in fuzzy_equipment_names:
                            similarity_score = fuzz.ratio(normalized_str.lower(), synonym.lower())
                            if similarity_score > best_match_score:
                                best_match_score = similarity_score
                                best_category = equipment_category
                                best_equipment = equipment_name
                                best_matched_string = synonym
                                normalized_str_matched = normalized_str
                            else:
                                continue

            if best_match_score >= 80:
                logging.debug(f"Hospital source identified for string: `{fuzzy_str}`. The identified isolation source is `{best_equipment}` in category `{best_category}` with a score of {best_match_score}. The best matched string was `{best_matched_string}` with `{normalized_str_matched}`.")
                return best_category, best_equipment, best_match_score
            
            logging.debug(f"Could't identify the source with the non-splitted form. Trying with splitted form.")
            parts = re.split(r'[-+:;/_,\s]|\band\b', part)
            parts = [part.strip() for part in parts]

            for normalized_str in parts:
                normalized_str = normalized_str.lower()

                for equipment_category, equipment_dictionaries in hospital_equipments.items():
                    for equipment_name, fuzzy_equipment_names in equipment_dictionaries.items():
                        for synonym in fuzzy_equipment_names:
                            similarity_score = fuzz.ratio(normalized_str.lower(), synonym.lower())
                            if similarity_score > best_match_score:
                                best_match_score = similarity_score
                                best_category = equipment_category
                                best_equipment = equipment_name
                                best_matched_string = synonym
                                normalized_str_matched = normalized_str
                            else:
                                continue

            if best_match_score >= 80:
                logging.debug(f"Hospital source identified for string: `{fuzzy_str}`. The identified isolation source is `{best_equipment}` in category `{best_category}` with a score of {best_match_score}. The best matched string was `{best_matched_string}` with `{normalized_str_matched}`.")
                return best_category, best_equipment, best_match_score

            else:
                logging.debug(f"The string: `{fuzzy_str}` didn't match any entry from the database. Highest match was with `{best_category}` with {best_match_score} score. Returning `Unknown`.")
                return "Unknown", "Unknown", 0
        else:   
            logging.debug(f"Given string: `{fuzzy_str}` is not a string, returning `Unknown`.")
            return "Unknown", "Unknown", 0
    
    except Exception as e:
        logging.error(f"An unexpected error occurred while finding hospital equipment for string `{fuzzy_str}`: {e}")
        return "Unknown", "Unknown", 0

def find_animal_category(fuzzy_str: str):

    animals_categorized = {
        "Unknown": {
            "Unknown": ["bird", "animal", "flying bird", "mammal", "amphibian", "vertebrate",],
        },
        "Insects":{
            "Miscellaneous": ["Cetonia aurata", "Rose Chafer", "C. aurata", "Tenebrio molitor", "Mealworm Beetle", "T. molitor","Gryllus bimaculatus", "Two-Spotted Cricket", "G. bimaculatus","Acheta domesticus", "House Cricket", "A. domesticus","Blatta orientalis", "Oriental Cockroach", "B. orientalis","Periplaneta americana", "American Cockroach", "P. americana","Lucilia sericata", "Common Green Bottle Fly", "L. sericata","Calliphora vomitoria", "Bluebottle Fly", "C. vomitoria","Musca domestica", "House Fly", "M. domestica","Apis mellifera", "Western Honeybee", "A. mellifera","Bombus terrestris", "Buff-Tailed Bumblebee", "B. terrestris","Vespa crabro", "European Hornet", "V. crabro","Lasius niger", "Black Garden Ant", "L. niger","Camponotus pennsylvanicus", "Black Carpenter Ant", "C. pennsylvanicus","Anopheles gambiae", "Malaria Mosquito", "A. gambiae","Culex pipiens", "Common House Mosquito", "C. pipiens","Aedes aegypti", "Yellow Fever Mosquito", "A. aegypti","Tettigonia viridissima", "Great Green Bush-Cricket", "T. viridissima","Drosophila melanogaster", "Common Fruit Fly", "D. melanogaster","Pyrrhocoris apterus", "Firebug", "P. apterus","Leptinotarsa decemlineata", "Colorado Potato Beetle", "L. decemlineata","Harmonia axyridis", "Harlequin Ladybird", "H. axyridis","Coccinella septempunctata", "Seven-Spotted Ladybird", "C. septempunctata","Melolontha melolontha", "Common Cockchafer", "M. melolontha","Chrysoperla carnea", "Green Lacewing", "C. carnea","Forficula auricularia", "Common Earwig", "F. auricularia","Galleria mellonella", "Greater Wax Moth", "G. mellonella","Plodia interpunctella", "Indian Meal Moth", "P. interpunctella","Ephestia kuehniella", "Mediterranean Flour Moth", "E. kuehniella","Tribolium castaneum", "Red Flour Beetle", "T. castaneum","Sitophilus oryzae", "Rice Weevil", "S. oryzae","Sitophilus granarius", "Granary Weevil", "S. granarius","Atta cephalotes", "Leafcutter Ant", "A. cephalotes","Polistes dominula", "European Paper Wasp", "P. dominula","Nezara viridula", "Southern Green Stink Bug", "N. viridula","Graphosoma italicum", "Italian Striped Bug", "G. italicum","Pyrrhocoris apterus", "Firebug", "P. apterus","Eurygaster testudinaria", "Tortoise Bug", "E. testudinaria","Palomena prasina", "Green Shield Bug", "P. prasina","Oncopeltus fasciatus", "Large Milkweed Bug", "O. fasciatus","Lygaeus equestris", "Black-and-Red Bug", "L. equestris","Oryctes nasicornis", "European Rhinoceros Beetle", "O. nasicornis","Lucanus cervus", "European Stag Beetle", "L. cervus","Dorcus parallelipipedus", "Lesser Stag Beetle", "D. parallelipipedus","Geotrupes stercorarius", "Dor Beetle", "G. stercorarius","Scarabaeus sacer", "Sacred Scarab", "S. sacer","Gerris lacustris", "Common Pond Skater", "G. lacustris","Notonecta glauca", "Common Backswimmer", "N. glauca","Corixa punctata", "Lesser Water Boatman", "C. punctata","Dytiscus marginalis", "Great Diving Beetle", "D. marginalis","Hydrophilus piceus", "Great Silver Water Beetle", "H. piceus","Cicindela campestris", "Green Tiger Beetle", "C. campestris","Carabus nemoralis", "Bronze Carabid", "C. nemoralis","Carabus auratus", "Golden Ground Beetle", "C. auratus","Calosoma sycophanta", "Forest Caterpillar Hunter", "C. sycophanta","Cylindera germanica", "German Tiger Beetle", "C. germanica","Mantis religiosa", "European Mantis", "M. religiosa","Empusa pennata", "Conehead Mantis", "E. pennata","Phyllium bioculatum", "Leaf Insect", "P. bioculatum","Peruphasma schultei", "Black Beauty Stick Insect", "P. schultei","Extatosoma tiaratum", "Giant Spiny Stick Insect", "E. tiaratum","Eurycantha calcarata", "Thorny Devil Stick Insect", "E. calcarata","Hierodula membranacea", "Giant Asian Mantis", "H. membranacea","Blaberus giganteus", "Giant Cave Cockroach", "B. giganteus","Gromphadorhina portentosa", "Madagascar Hissing Cockroach", "G. portentosa","Pachnoda marginata", "Sun Beetle", "P. marginata","Eriosoma lanigerum", "Woolly Apple Aphid", "E. lanigerum","Aphis nerii", "Oleander Aphid", "A. nerii","Myzus persicae", "Green Peach Aphid", "M. persicae","Corythucha ciliata", "Sycamore Lace Bug", "C. ciliata","Halyomorpha halys", "Brown Marmorated Stink Bug", "H. halys","Pyrrhalta viburni", "Viburnum Leaf Beetle", "P. viburni","Pyrrhalta luteola", "Elm Leaf Beetle", "P. luteola","Xyleborus dispar", "European Shot Hole Borer", "X. dispar","Ips typographus", "European Spruce Bark Beetle", "I. typographus","Agrilus planipennis", "Emerald Ash Borer", "A. planipennis","Anoplophora glabripennis", "Asian Longhorned Beetle", "A. glabripennis","Tomicus piniperda", "Pine Shoot Beetle", "T. piniperda","Megachile rotundata", "Alfalfa Leafcutter Bee", "M. rotundata","Anthidium manicatum", "Wool Carder Bee", "A. manicatum","Xylocopa violacea", "Violet Carpenter Bee", "X. violacea", "cricket powder", "cricket", 'insect', 'insects', 'bug', 'bugs', 'arthropod', 'arthropods', 'insecta', 'insectae', 'insectum', 'insectums', 'insectoid', 'insectoids', 'cricket', 'buffalo fly', 'drosophila', 'fruit fly', 'housefly', 'house fly', 'house-fly', 'house-fly', "melipona lateralis", "stingless bee", "melipona", "stingless bees", "hymenoptera", "hymenopteran", "hymenopterans", "hymenopteroid", "hymenopteroids", "hymenopterous", "hymenoptera insects", "hymenoptera insect", "melipona seminigra", "scaptotrigona polysticta", "scaptotrigona", "apidae", "apidae bees", "apidae bee", "hymenoptera insects", "benjoi", "bijui", "frieseomelitta varia", "frieseomelitta", "abelha marmelada-amarela", "yellow marmalade bee", "melipona interrupta", "jandaira", "planococcus ficus", "vine mealybug", "lagria villosa", "bettle", "lagria grenieri", "lagria atripes", "ecnolagria sp.", "lagria rufipennis", "lagria okinawana", "lagria", "lagria sp.", "stegodyphus dumicola", "african social spider", "Silkworms", "hermetia illucens", "black soldier fly", "isolate from fly samples", "Dryocosmus kuriphilus", "chestnut gall wasp", "D. kuriphilus", "Dryocosmus", "kuriphilus", "gall wasp", "gall wasps", "gall insect", "gall insects", "gall-inducing insect", "gall-inducing insects", "gall-forming insect", "gall-forming insects", "gall-maker insect", "gall-maker insects", "gall-inducer insect", "gall-inducer insects", "cynipidae", "cynipid wasps", "cynipid wasp", "cynipids", "cynipidae family", "nezara viridula", "southern green stink bug", "N. viridula", "graphosoma italicum", "italian striped bug", "G. italicum", "pyrrhocoris apterus", "firebug", "P. apterus", "eurygaster testudinaria", "tortoise bug", "E. testudinaria", "palomena prasina", "green shield bug", "P. prasina", "oncopeltus fasciatus", "large milkweed bug", "O. fasciatus", "lygaeus equestris", "black-and-red bug", "L. equestris", "orius albidipennis", "whitefly", "O. albidipennis", "trialeurodes vaporariorum", "greenhouse whitefly", "T. vaporariorum", "bemisia tabaci", "silverleaf whitefly", "B. tabaci", "aphididae", "aphids", "aphid", "aphid family", "corythucha ciliata", "sycamore lace bug", "C. ciliata", "halyomorpha halys", "brown marmorated stink bug", "H. halys", "pyrrhalta viburni", "viburnum leaf beetle", "P. viburni", "pyrrhalta luteola", "elm leaf beetle", "P. luteola", "bombyx mori", "silkworm", "B. mori", "bombyx", "mori", "silkworms", "hermetia illucens", "black soldier fly", "flesh fly"],    
        },
        "Livestock": {   
            "Bovine": ["Bovine", "Cow", "Cattle", "Bos Taurus", "Ox", "Beef Cattle", "Dairy Cattle", "Bull", "Heifer", "Dairy Cow", "Beef Cow", "Female Bovine", "Bison", "American Bison", "Bison Bison", "Buffalo", "Buffalo Bison", "Water Buffalo", "Bubalus Bubalis", "African Buffalo", "Syncerus Caffer", "Asian Buffalo", "Yoke Ox", "Working Ox", "Draught Ox", "beef", "buff", "buffalo meat", "domestic cattle", "bos indicus", "zebu", "bovine booties", "bos taurus indicus", "steers", "bubaline", "bos taurus indicus", "bos taurus jersey", "jercey cow", "jersey buffalo"],
            "Equine": ["Equine", "Horse", "Pony", "Equus Ferus Caballus", "Steed", "Mare", "Stallion", "Colt", "Equus Ferus Caballus", "Filly", "Donkey", "Ass", "Equus Asinus", "Burro", "Jack", "Burro Donkey", "Mule", "Donkey-Horse Hybrid", "Equus Mulus", "Mule Horse", "horse meat", "horse leg", "hoof", "equine meat"],
            "Poultry": ["Poult","guinea hen","guinea fowl","Gallus gallus","Goose","hen farm","hen","turkey","Meleagris gallopavo","Poultry", "Fowl", "Chicken", "Gallus Gallus Domesticus", "Domestic Chicken", "Hen", "Rooster", "Duck", "broiler", "chick", "small chicks", "chick box", "chicken meat", "chicken eggs", "duck meat","duck eggs", "Black bone chicken", "Ayam Cemani", "chicken fluff", "chicken fluffy"],
            "Sheep": ["Sheep", "Ovis Aries", "Ram", "Ewe", "Lamb", "Mutton Sheep", "Ovine", "sheep meat", "sheep bone", "lamb"],
            "Goat": ["Goat", "Capra Aegagrus Hircus", "Buck", "Doe", "Domestic Goat", "goat meat", "Capra", "Capra Aegagrus", "chevon", "capra hircus"],
            "Porcine": ["Porcine","Pig", "Swine", "Hog", "Sus Scrofa Domesticus", "Sus Scrofa", "Boar", "Sow", "pork", "pig meat", "pork", "pork sausage", "ham", "longaniza", "spanish sausage", "embutido", "sus scrofa cristatus"],
            "Alpaca": ["Alpaca", "Vicugna Pacos", "Domestic Alpaca", "Alpaca Camelid"],
            "Camelid": ["Camelid", "Camelidae", "South American Camelid", "Andean Camelid", "New World Camelid"],
            "Llama": ["Llama", "Lama Glama", "Domestic Llama", "Llama Camelid"],
            "Deer": ["Deer", "Cervidae", "Stag", "Doe", "Fawn", "White-Tailed Deer", "mule deer", "Odocoileus hemionus"],
            "Moose": ["Moose", "Alces Alces", "Elk (In Europe)", "European Elk", "Alces Alces Americanus"],
            "Caribou": ["Caribou", "Rangifer Tarandus", "Reindeer", "Caribou Reindeer"],
            "Mountain Goat": ["Mountain Goat", "Oreamnos Americanus", "Rock Goat", "Wild Goat"],
            "Camel": ["Camel", "Camelus", "Dromedary", "Bactrian Camel", "Desert Camel", "Camelus dromedarius", "Camelus bactrianus"],
            "Yak": ["Yak", "Bos Grunniens", "Domestic Yak", "Himalayan Yak"],
            "Miscellaneous Birds": ["Peacock", "Pavo Cristatus", "Canary", "Serinus Canaria", "Ostrich", "Struthio Camelus", "Quail", "Coturnix Coturnix", "Pigeon", "Columba Livia", "Indian Peafowl", "Pavo Cristatus", "Rock Pigeon", "Columba Livia", "Budgerigar", "Melopsittacus Undulatus", "domestic bird"],
            "Others": ["Livestock", "Farm Animals", "Tame Animals", "Domestic Animals", "Working Animals", "Agricultural Animal", "Husbandry Animal", "Livestock", "Farm Stock", "Cattle Stock", "Farming Animals", "Domesticated Animals", "Feedlot Water Bowl"],
            "Livestock Breeding Environment/Equipments" : ["Livestock Breeding Environment", "Livestock Farm", "Livestock Breeding Farm", "Livestock Production Setting", "Livestock Farming Environment", "Livestock Breeding Facility", "Livestock Farming Facility", "Livestock Production Facility", "Livestock Farming Setting", "Livestock Breeding Setting", "Boot Swab", "Poultry farm", "pig farm", "Barn", "Avian incubator", "Avian environmental incubator", "farm maternity bed", "duck incubator", "buffalo production setting", "beef production setting", "cattle production setting", "poultry production setting", "chicken production setting", "pig production setting", "egg farm cage shed", "carcass cleaning wipe", "dairy equipment filters", "meat processing equipment", "milking machine", "colostrum feed bucket", "boot cover", "boot cover swab", "boot cover sample", "boot cover culture", "boot cover isolate", "boot cover cleaning", "boot cover disinfection", "boot cover sterilization", "boot cover wipe", "boot cover wipe sample", "boot cover wipe culture", "farm animal", "farm animal environment", "farm animal equipment", "farm animal setting", "farm animal facility", "farm animal breeding environment", "farm animal breeding facility", "bovine calf shed", "bovine calf barn", "bovine calf pen", "bovine calf house", "bovine calf rearing facility", "bovine calf rearing environment", "bovine calf rearing setting", "bovine calf rearing farm", "bovine calf rearing unit", "bovine calf rearing area", "bovine calf rearing equipment", "bovine calf rearing system", "bovine calf rearing shed", "bovine calf rearing barn", "bovine calf rearing pen", "bovine farm", "cattle farm", "cattle breeding farm", "cattle production farm", "cattle rearing farm", "cattle housing farm", "cattle raising farm", "cattle feeding farm", "cattle fattening farm", "cattle finishing farm", "cattle grazing farm", "cattle pasture farm", "cattle feedlot farm", "cattle feedlot environment", "cattle feedlot facility", "cattle feedlot setting", "cattle feedlot equipment", "arsenic contaminated cattle tick dip sites", "cattle tick dip site soil", "cattle tick dip site water", "poultry house dust", "poultry house dust sample", "poultry house dust swab", "poultry house environment", "poultry house environment sample", "poultry house environmental swab", "poultry house litter sample", "poultry house litter swab", "poultry house manure sample", "poultry house manure swab", "poultry house floor swab", "poultry house wall swab", "poultry house ceiling swab", "poultry house equipment swab", "poultry house tool swab", "poultry house surface swab", "poultry house bedding swab", "poultry house air sample", "poultry house ventilation system swab", "poultry house heating system swab", "poultry house cooling system swab", "poultry house lighting system swab", "tick dip site", "tick dip site soil", "tick dip site water", "veterinary clinic environment", "veterinary clinic environment sample", "veterinary clinic environmental swab", "veterinary hospital environment", "veterinary hospital environment sample", "veterinary hospital environmental swab", "animal hospital environment", "animal hospital environment sample", "animal hospital environmental swab", "boot swab sample", "boot swab environmental sample", "boot swab environmental swab", "boot swab farm environment", "boot swab farm environmental sample", "boot swab farm environmental swab", "cow paddock straw bedding bovine", "paddock straw bedding bovine", "cow paddock straw bedding", "paddock straw bedding", "livestock barn environment", "livestock barn environment sample", "livestock barn environmental swab", "livestock pen environment", "livestock pen environment sample", "livestock pen environmental swab", "livestock shed environment", "livestock shed environment sample", "livestock shed environmental swab", "livestock facility environment", "livestock facility environment sample", "livestock facility environmental swab", "livestock equipment swab", "livestock tool swab", "livestock surface swab", "livestock bedding swab", "dairy farm", "dairy farm soil", "dairy farm field soil", ""],
            
        },
        
        "Wild": {
            "Others": ["lowland paca","agouti lowland paca","alligator mississippiensis", "caiman","Central bearded dragon", "Pogona vitticeps", "lophosaurus boydii", "Boyd's forest dragon", "elaphe taeniura", "correlophus ciliatus", "crested gecko","Land Iguanas","conolophus subcristatus","Wild Animal", "Wildlife", "Feral Animal", "Untamed Animal", "Non-Domesticated Animal", "Endangered", "Threatened", "Vulnerable", "Conservation-Dependent", "Wildlife", "Wild Animals", "Nature", "Endangered Wildlife", "Fauna", "Pythonidae", "Soricidae", "Shrew", "Apis Mellifera", "Bee", "Frilled Lizard", "Frilled-neck lizard", "Chlamydosaurus kingii", "Four-toed hedgehog", "atelerix albiventris", "Four-toed hedgehog", "atelerix albiventris", "lizard", "house lizard","lumbricina", "earthworm", "annelida", "annelid", "annelids", "gecko", "reptile", "callithrix penicillata", "black-tufted marmoset", "tapirus terrestris", "brazilian tapir", "tapir", "tapirus", "tapirus sp", "dassie", "rock hyrax", "Suricata suricatta", "Meercat", "neogale vison", "equus asinus africanus", "african wild ass", "canis latrans", "coyote", "pogona", "skink", "lampropholis guichenoti", "agalychnis callidryas", "craugastor fitzingeri"],
            "Deer" : ["axis porcinus", "axis porcinus", "chital deer", "axis axis", "sambar deer", "cervus unicolor", "cervus elaphus", "red deer", "cervus canadensis", "elk (in North America)", "wapiti", "odocoileus virginianus", "white-tailed deer", "odocoileus hemionus", "hog deer", "Indian hog deer"],
            "Antelope": ["Antelope", "blesbok", "oryx", "common eland", "Roan antelope", "Antilopinae", "Gazelle", "Impala", "Springbok", "Oryx", "kudu", "tragelaphus oryx", "oryx gazella", "Gemsbok", "Gazella subgutturosa", "Gazella sp.", "Persian gazelles", "gazella cuvieri", "blesbok", "blesbuck", "Damaliscus pygargus phillipsi", "Damaliscus pygargus", "Damaliscus lunatus jimela", "topi antelope"],
            "Procyonidae": ["Procyonidae", "Raccoon", "Procyon Lotor", "Ringtail", "Bassariscus astutus", "Cacomistle", "Bassariscus sumichrasti", "Kinkajou", "Potos flavus"],
            "Mustelidae": ["Mastelidae", "Mustelid", "Weasel", "Ferret", "american mink", "mink", "Badger", "Otter", "Mink", "Wolverine", "Martes", "Lutra", "Lutreolina", "Neovison", "Mustela", "Gulo", "Martes", "Meles", "Lutra", "Lutreolina", "Neovison", "Mustela", "European pine marten", "Martes martes", "American marten", "Martes americana", "European badger", "Meles meles", "American badger", "Taxidea taxus", "European otter", "Lutra lutra", "North American river otter", "Lontra canadensis", "European polecat", "Mustela putorius", "North American river otter", "Lontra canadensis"],
            "Snake": ["Snake", "Serpentes", "Python", "Boa", "Viper", "Cobra", "Anaconda", "Rattlesnake", "Garter Snake", "King Cobra", "Coral Snake", "Copperhead", "Cottonmouth", "python curtus", "pantherophis guttatus", "beauty rat snake", "rat snake", "rattle snake", "venomous snake", "non-venomous snake", "lampropeltis triangulum", "milk snake", "animal-snake"],
            "Elephant": ["Elephant", "Loxodonta Africana", "African Elephant", "Asian Elephant", "Elephas Maximus", "Indian Elephant",],
            "Rodent":["Beaver", "microtus agrestis", "vole", "wild vole", "hydrochoerus hydrochaeris", "capybara", "marmota himalayana", "himalayan marmot", "marmot", "myodes glareolus", "bank vole", "clethrionomys glareolus", "clethrionomys", "clethrionomys sp", "clethrionomys sp.", "Field vole", "Vole", "apodemus flavicollis", "yellow-necked field mouse", "apodemus sylvaticus", "wood mouse", "apodemus agrarius", "striped field mouse", "apodemus chevrieri", "chevrier's field mouse", "apodemus speciosus", "japanese wood mouse", "apodemus uralensis", "ural field mouse", "apodemus peninsulae", "korean field mouse", "apodemus agrarius", "striped field mouse", "apodemus chevrieri", "chevrier's field mouse", "apodemus speciosus", "japanese wood mouse", "apodemus uralensis", "ural field mouse", "apodemus peninsulae", "korean field mouse", "apodemus sylvaticus", "wood mouse", "apodemus flavicollis", "yellow-necked field mouse", "rodent deer mouse", "deer mouse", "apodemus uralensis", "ural field mouse", "mus spicilegus", "steppe mouse", "rattus rattus", "black rat", "ship rat", "house rat", "roof rat"],

            "Frog/Toad": ["Frog", "Toad", "Anura", "Rana", "Bufo", "Lithobates", "Hyla", "Pseudacris", "Pelophylax", "Acris", "Anaxyrus", "Ranitomeya", "Adenomera", "Aparasphenodon", "Aparasphenodon brunoi", "Aparasphenodon arapuca", "Aparasphenodon sp.", "Aparasphenodon sp. nov.", "bufo bufo", "common toad", "bufo viridis", "green toad", "bufo calamita", "natterjack toad", "bufo japonicus", "japanese toad", "bufo boreas", "western toad", "bufo cognatus", "great plains toad", "rana pipiens", "northern leopard frog", "rana catesbeiana", "bullfrog", "rana sylvatica", "wood frog", "rana clamitans", "green frog", "rana temporaria", "common frog", "rana aurora", "red-legged frog", "rana muscosa", "mountain yellow-legged frog", "rana draytonii", "bullfrog", "rana boylii", "foothill yellow-legged frog", "rana pipiens", "northern leopard frog", "rana catesbeiana", "bullfrog", "rana sylvatica", "wood frog", "rana clamitans", "green frog", "rana temporaria", "common frog", "rana aurora", "red-legged frog", "rana muscosa", "mountain yellow-legged frog", "rana draytonii", "bullfrog", "american bullfrog", "rana boylii", "foothill yellow-legged frog", "baby toad from irrigation pond", "baby toad", "baby frog", "frog tadpole", "frog tadpoles", "frog spawn", "frog eggs", "froglet", "froglets", "toad tadpole", "toad tadpoles", "toad spawn", "toad eggs", "toadlet", "toadlets"],
            "Felidae": ["felidae"],
            "Wallabies" : ["Tammar wallaby", "macropus eugenii", "parma wallaby", "wallaby", "agile wallaby", "swamp wallaby"],

            "Marsupial": ["Marsupial", "Marsupialia", "Macropodidae", "Wallaby", "Wombat", "Tasmanian Devil", "Thylacine", "Dasyuromorphia", "Didelphimorphia", "Peramelemorphia", "Diprotodontia", "Petauridae", "Phalangeridae", "Macropodinae", "opossum", "didelphidae", "dasyuridae", "peramelidae", "petauridae", "phalangeridae", "macropodinae", "wombat", "vombatidae", "didelphis marsupialis", "common opossum", "didelphis virginiana", "virginia opossum", "opossum sp.", "opossum sp", "opossum", "wallaby", "wallabies", "swamp wallaby", "agile wallaby", "tammar wallaby", "macropus eugenii", "parma wallaby", "Virginia opossum",],
            
            "Murine": ["Mice", "Rat", "Murine", "Mus musculus"],
            "Terrestrial Mollusks": ["Terrestrial Mollusks", "Gastropoda", "Pulmonata", "Stylommatophora", "Eupulmonata", "Heterobranchia", "Basommatophora", "Eupulmonata", "Stylommatophora", "achatina achatina", "Giant African Snail", "Achatina fulica", "Giant African Land Snail", "Helix aspersa", "Common Garden Snail", "Cornu aspersum", "Helix pomatia", "Roman Snail", "Eobania vermiculata", "Mediterranean Land Snail", "helix aspersa muller", "Apple snail", "Pomacea canaliculata", "Pomacea bridgesii", "Golden Apple Snail", "Pomacea diffusa", "Pomacea paludosa", "Florida Apple Snail", "Ampullariidae", "Ampullariidae Family"],
            "Tiger": ["Tiger", "Panthera Tigris", "Bengal Tiger", "Indochinese Tiger", "Siberian Tiger", "Malayan Tiger", "panthera tigris altaica"],
            "Jackal": ["Jackal", "Canis Aureus", "Golden Jackal", "Black-Backed Jackal", "Side-Striped Jackal"],
            "Lion": ["Lion", "Panthera Leo", "African Lion", "Asiatic Lion", "King Of The Jungle"],
            "Jaguar": ["Jaguar", "Panthera Onca", "Panthera Onca Palustris"],
            "American Lion": ["Panthera Leo Atrox", "American Lion", "Pleistocene Lion"],
            "Leopard" : ["Leopard", "Panthera Pardus", "Black Leopard", "panthera pardus japonensis", ],
            "Snow Leopard" : ["Panthera uncia", "Snow Leopard"],
            "Bear": ["Bear", "Ursidae", "Grizzly Bear", "Black Bear", "Polar Bear", "Brown Bear", "Ursus Arctos", "ursus", "ursus sp"],
            "Wolf": ["Wolf", "Canis Lupus Dingo", "Gray Wolf", "Timber Wolf", "Arctic Wolf", "Red Wolf"],
            "Kangaroo": ["Kangaroo", "Macropodidae", "Red Kangaroo", "Eastern Gray Kangaroo", "Western Gray Kangaroo", ""],
            "Rat Kangaroo" : ["dipodomys sp."],
            "Koala": ["Koala", "Phascolarctos Cinereus", "Koala Bear", ],
            "Giraffe": ["Giraffe", "Giraffa Camelopardalis", "Reticulated Giraffe", "Masai Giraffe", "Rothschild's Giraffe"],
            "Panda": ["Panda", "Ailuropoda Melanoleuca", "Giant Panda", "Bamboo Bear"],
            "Zebra": ["Zebra", "Equus Zebra", "Plains Zebra", "Grévy's Zebra", "Mountain Zebra"],
            "Bat": ["Bat", "Chiroptera", "Fruit Bat", "Vampire Bat", "Flying Mammal", "eidolon helvum", "africal fruit bat", "eptesicus fuscus", "big brown bat", "lasionycteris noctivagans", "silver-haired bat", "myotis ciliolabrum", "western small footed bat", "corynorhinus townsendii", "Townsend's big-eared bat", "myotis volans", "long-legged myotis", ],
            "Fox": ["Fox", "Vulpes Vulpes", "Red Fox", "Arctic Fox", "Fennec Fox", "Swift Fox", "Kit Fox", "Gray Fox", "Corsac Fox", "Bengal Fox", "Blanford's Fox"],
            "Mole": ["Mole", "Talpidae", "European Mole", "Blind Mole", "Earth Mole"],
            "Anteater": ["Anteater", "Myrmecophagidae", "Giant Anteater", "Tamandua"],
            "Armadillo": ["Armadillo", "Dasypodidae", "Nine-Banded Armadillo", "Southern Armadillo"],
            "Sloth": ["Sloth", "Folivora", "Two-Toed Sloth", "Three-Toed Sloth", "Sloth Bear"],
            "Monkey": ["Monkey", "Primates", "Squirrel Monkey", "Capuchin Monkey", "Macaque", "Golden lion tamarin", "Callitrichidae", "mico", "marmosets", "lion-tailed macaque", "Macaca silenus", ],
            "Ape": ["Ape", "Great Ape", "Chimpanzee", "Gorilla", "Orangutan", "Celebes crested macaque", "Hylobatidae", "Gibbons", "Hominidae", "Homininae", "Hominini", "Pongo", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "Orangutan", "Pongo Pygmaeus", "Pongo Abelii", "Pongo Tapanuliensis", "celebes ape",],
            "Gorilla": ["Gorilla", "Gorilla Gorilla", "Silverback Gorilla", "Mountain Gorilla"],
            "Chimpanzee": ["Chimpanzee", "Pan Troglodytes", "Chimp", "Bonobo", "Ape"],
            "Hare" : ["Hare", "Lepus", "Jackrabbit", "Snowshoe Hare", "European Hare"],

            "Miscellaneous Birds": ["black-eared kite", "milvus migrans", "milvus lineatus lineatus","Osprey","pandion haliaetus","Phalacrocorax","cormorant", "house sparrow","Passer domesticus","penguin","Humboldt penguin","spheniscus humboldti","silver gull","chroicocephalus novaehollandiae","American Robin", "Turdus Migratorius", "Bald Eagle", "Haliaeetus Leucocephalus", "Barn Owl", "Tyto Alba", "Blue Jay", "Cyanocitta Cristata", "Canada Goose", "Branta Canadensis", "Common Raven", "Corvus Corax", "Crow", "Corvus Brachyrhynchos", "European Starling", "Sturnus Vulgaris", "Golden Eagle", "Aquila Chrysaetos", "Great Horned Owl", "Bubo Virginianus", "Northern Cardinal", "Cardinalis Cardinalis", "Peregrine Falcon", "Falco Peregrinus", "Red-Tailed Hawk", "Buteo Jamaicensis", "Snowy Owl", "Bubo Scandiacus", "Southern Cassowary", "Casuarius Casuarius", "Turkey Vulture", "Cathartes Aura", "Woodpecker", "Picidae Family", "Yellow Warbler", "Setophaga Petechia", "Kingfisher", "Alcedinidae Family", "Golden Eagle", "Aquila Chrysaetos", "Vulture", "Accipitridae Family", "Kookaburra", "Dacelo Novaeguineae", "Cuckoo", "Cuculidae Family", "Wren", "Troglodytidae Family", "Gannet", "Morus Bassanus", "Bowerbird", "Ptilonorhynchidae Family", "Puffin", "Fratercula Arctica", "Secretary Bird", "Sagittarius Serpentarius", "Tern", "Sternidae Family", "Nightjar", "Caprimulgidae Family", "Hornbill", "Bucerotidae Family", "Crowned Crane", "Balearica Pavonina", "Red-Crowned Crane", "Grus Japonensis", "Northern Flicker", "Colaptes Auratus", "Vermilion Flycatcher", "Pyrocephalus Rubinus", "Harpy Eagle", "Harpia Harpyja", "Bald-Headed Ibis", "Geronticus Eremita", "Eurasian Kestrel", "Falco Tinnunculus", "Barnacle Goose", "Branta Leucopsis", "Horned Lark", "Eremophila Alpestris", "Red-Winged Blackbird", "Agelaius Phoeniceus", "Mountain Bluebird", "Sialia Currucoides", "Lesser Kestrel", "Falco Naumanni", "Golden Pheasant", "Chrysolophus Pictus", "Emu", "Dromaius Novaehollandiae", "Larus dominicanus", "Kelp gull", "spinus spinus", "Eurasian siskin", "mareca strepera", "Gadwall", "white stork", "ciconia ciconia", "Migratory bird", "bubulcus ibis", "cattle egret", "columbiformes", "columba livia", "rock pigeon", "columba palumbus", "common wood pigeon", "columba oenas", "stock dove", "columba livia domestica", "domestic pigeon", "columba palumbus", "common wood pigeon", "columba oenas", "stock dove", "anas", "duck", "Anatidae", "Anatidae Family", "Anas Platyrhynchos", "Mallard Duck", "Anas Crecca", "Green-Winged Teal", "Anas Penelope", "Eurasian Wigeon", "Anas Acuta", "Northern Pintail", "Anas Clypeata", "Northern Shoveler", "Anas Discors", "Blue-Winged Teal", "Anas Americana", "American Black Duck", "phasianus colchicus", "house finch", "carpodacus mexicanus", "common redpoll", "carduelis flammea", "redpoll", "carduelis flammea", "common redpoll", "carduelis flammea", "redpoll", "carduelis flammea", "common redpoll", "carduelis flammea", "redpoll", "carduelis flammea", "common redpoll", "carduelis flammea", "columbidae", "american goldfinch", "spinus tristis", "goldfinch", "carduelis carduelis", "european goldfinch", "carduelis carduelis", "goldfinch", "carduelis carduelis", "european goldfinch", "carduelis carduelis", "goldfinch", "carduelis carduelis", "european goldfinch", "carduelis carduelis", "Pied oystercatche", "Haematopus longirostris", "oystercatcher", "sooty oystercather", "haemotopus sp.", "cygnus olor", "mute swans", "eagle", "magellanic penguin", "penguin", "magnificent frigatebird", "Demoiselle crane", "Grus Virgo", "Ruddy Turnstone", "Arenaria Interpres", "Pigeon Guillemot", "Cepphus Columba", "Black Skimmer", "Rynchops Niger", "Common Murre", "Uria Aegis", "Razorbill", "Alca Torda", "Atlantic Puffin", "Fratercula Arctica", "Tufted Puffin", "Fratercula Cirrhata", "anthropoides virgo", "Magellanic penguin", "spheniscus magellanicus", "anas poecilorhyncha", "Indian Runner Duck", "anas platyrhynchos domesticus", "domestic duck", "anas platyrhynchos", "mallard duck", "anas crecca", "green-winged teal", "anas penelope", "eurasian wigeon", "anas acuta", "northern pintail", "anas clypeata", "northern shoveler", "anas discors", "blue-winged teal", "anas americana", "american black duck"],
        },

        "Aquatic": {
            "Others": ["Aquatic Animal", "Marine Creature", "Water Animal", "Ocean Creature", "Water-Dwelling Animal", "Aquatic Animals", "Marine Life", "Water Creatures", "Ocean Fauna", "Marine Wildlife",],
            "Turtle" : ["Turtle", "Chelonia", "Testudines", "Green Sea Turtle", "Green Turtle", "Chrysemys Picta", "Painted Turtle", "Chelonia Mydas", "Loggerhead Sea Turtle", "Caretta Caretta", "Leatherback Sea Turtle", "Dermochelys Coriacea", "Hawksbill Sea Turtle", "Eretmochelys Imbricata"],
            "Fish": ["Fried fish","Fish", "Pisces", "Freshwater Fish", "Saltwater Fish", "Aquatic Fish", "Sardines", "Goldfish", "Carassius Auratus", "Common Goldfish", "Golden Fish", "Mahi-mahi", "Coryphaena hippurus", "Bass", "Black bass", "oncorhynchus mykiss", "Rainbow trout", "Salmo trutta", "Brown trout", "Salvelinus fontinalis", "Brook trout", "Oncorhynchus mykiss", "Rainbow trout", "Salmo salar", "Atlantic salmon", "Salvelinus namaycush", "Lake trout", "Micropterus salmoides", "Largemouth bass", "Micropterus dolomieu", "Smallmouth bass", "salvelinus namaycush", "lake trout", "salmo trutta", "brown trout", "salmo salar", "atlantic salmon", "oncorhynchus mykiss", "rainbow trout", "som fak", "carassius auratus", "goldfish", "coryphaena hippurus", "mahi-mahi", "black bass", "oncorhynchus mykiss", "rainbow trout", "scarlet snapper", "lutjanus campechanus", "red snapper", "lutjanus campechanus", "lutjanidae", "snapper", "lutjanidae family", "lutjanus campechanus", "red drum", "sciaenops ocellatus", "redfish", "sciaenidae", "sciaenidae family", "sciaenops ocellatus", "bream", "sparidae", "sparidae family", "sardine", "clupeidae", "clupeidae family", "sardina pilchardus", "pilchard", "sardina pilchardus", "pilchardus", "sardina pilchardus", "pilchardus", "frozen striped bass", "bass", "striped bass", "mullet","ray finned fish", "ayer fish", "bacha fish", "baby mullet fish", "baila fish", "kajoli fish", "barb fish", "barramundi fillet fish", "basa fish", "black drum fish", "drum fish", "black pomfret fish", "blolovi smoked fish", "blue mackerel fish", "boal fish", "climbing perch fish", "climbing perch", "cobia fish", "cobias fish", "cod fish", "common carp fish", "corydoras fish", "corydoras catfish", "crucian carp fish", "dace fish", "dorado fish", "drum fish", "eel fish", "eel catfish", "eel fish fillet", "eel fillet fish", "eel smoked fish", "eel smoked fillet fish", "eel smoked whole fish", "eel whole fish", "frozen eel fillet fish", "frozen eel whole fish", "frozen eel smoked fillet fish", "marine fish", "fresh water fish", "salt water fish", "freshwater fish", "saltwater fish", "cooked mackerel fish", "etroplus suratensis", "atlantic saury", "scomberesox saurus", "mahi mahi", "coryphaena hippurus", "dolphin fish", "mahi-mahi", "coryphaena hippurus", "dolphin fish", "mahi-mahi", "coryphaena hippurus", "mugil cephalus", "Flathead grey mullet", "mugil cephalus", "grey mullet", "mugil cephalus", "flathead grey mullet", "mugil cephalus", "grey mullet", "mugil cephalus", "flathead grey mullet", "mugil cephalus", "grey mullet", "mugil cephalus", "perch", "perch fillet", "tuna", "frozen tuna", "frozen tuna fish", "frozen tuna fillet", "tuna fish", "nile tilapia", "oreochromis niloticus", "tilapia", "oreochromis", "oreochromis niloticus", "tilapia", "oreochromis", "oreochromis niloticus", "tilapia", "oreochromis", "oreochromis niloticus", "tilapia", "oreochromis", "non ictalurus fish", "oncorhynchus tshawytscha", "labeo rohita", "grass carp", "Ctenopharyngodon idell", "black carp", "silver carp", "ctemopharyngodon sp.", "finless porpoise", "centroscymnus coelolepis", "dogfish", "bathysaurus ferox", "deepsea lizard fish", "ictalurus punctatus", "channel catfish", "Anchovy", "anchovy fish", "broadhead fish", "Heteropneustes fossilis", "stingging catfish", "lotia fish", "lotiidae", "bombay duck", "bummalo", "Harpadon nehereus", "Rainbow Trout", "Oncorhynchus mykiss", "Salmonidae", "Salmo gairdneri", "Salmo irideus", "Parasalmo mykiss", "Steelhead", "Kamloops trout", "Freshwater / Anadromous fish", "Tilapia", "Oreochromis sp.", "Cichlidae", "Tilapia sp.", "Nile tilapia", "Mango fish", "St. Peter's fish", "Freshwater fish", "Nile Tilapia", "Oreochromis niloticus", "Tilapia nilotica", "Chromis niloticus", "Bolti", "Climbing Perch", "Anabas testudineus", "Anabantidae", "Anabas scandens", "Perca testudinea", "Walking fish", "Climbing fish", "Bathysaurus ferox", "Bathysauridae", "Macristium chavesi", "Deepsea lizardfish", "Marine / Deep-sea fish", "Goldfish", "Carassius auratus", "Cyprinidae", "Carassius auratus auratus", "Cyprinus auratus", "Common goldfish", "Fancy goldfish", "Carassius auratus auratus", "Centroscymnus coelolepis", "Somniosidae", "Scymnodon coelolepis", "Portuguese shark", "Cirrhinus mrigala", "Cirrhina mrigala", "Cyprinus mrigala", "Mrigal", "White carp", "Clarias magur", "Clariidae", "Clarias batrachus", "Indian walking catfish", "Magur", "Cyprinus carpio var. specularis", "Cyprinus specularis", "Mirror carp", "Scale carp", "Ictalurus punctatus", "Ictaluridae", "Silurus punctatus", "Pimelodus punctatus", "Channel cat", "Spotted catfish", "Labeo catla", "Cyprinus catla", "Gibelion catla", "Catla", "Major Indian carp", "Labeo rohita", "Cyprinus rohita", "Rohita buchananit", "Rohu", "Rui", "Tapra", "Lepidocephalichthys guntea", "Cobitidae", "Cobitis guntea", "Lepidocephalus guntea", "Guntea loach", "Peppered loach", "Ctenopharyngodon idella", "Leuciscus idella", "Ctenopharyngodon idellus", "Grass Carp", "White amur", "Ophichthys cuchia", "Ophichthidae", "Amphipnous cuchia", "Monopterus cuchia", "Cuchia", "Swamp eel", "Sperata seenghala", "Bagridae", "Mystus seenghala", "Aorichthys seenghala", "Seenghala catfish", "Giant river catfish", "Seriola dumerili", "Carangidae", "Seriola purpurascens", "Caranx dumerili", "Greater Amberjack", "Medregal", "Marine fish", "Seriola lalandi", "Seriola aureovittata", "Seriola dorsalis", "Yellowtail Amberjack", "Yellowtail kingfish", "Hiramasa", "California yellowtail", "Acanthopagrus latus", "Sparidae", "Sparus latus", "Mylio latus", "Yellowfin Seabream", "Yellowfin bream", "Silver bream", "Mugil cephalus", "Mugilidae", "Mugil braziliensis", "Mugil lineatus", "Mullet", "Flathead grey mullet", "Striped mullet", "Black mullet", "Tetraodon sp.", "Lagocephalus sp.", "Tetraodontidae", "Fugu sp.", "Spheroides sp.", "Puffer Fish", "Blowfish", "Fugu", "Globefish", "Toadfish", "Perca fluviatilis", "Perca flavescens", "Percidae", "Perch", "European perch", "River perch", "Yellow perch", "Micropterus salmoides", "Dicentrarchus labrax", "Centrarchidae", "Moronidae", "Huro salmoides", "Morone labrax", "Bass", "Largemouth bass", "Striped bass", "Sea bass", "European bass", "Grystes salmoides", "Black Bass", "Bigmouth bass", "Green bass", "dhani koral", "channa punctata", "spotted snakehead", "Channa striata", "striped snakehead", "Channa marulius", "bullseye snakehead", "Channa micropeltes", "giant snakehead", "Channa argus", "northern snakehead", "Channa asiatica", "small snakehead", "Channa bleheri", "Rainbow Snakehead", "Channa gachua", "dwarf snakehead", "ompok bimaculatus", "giant snakehead", "sheatfish", "butter catfish", "puffer skin", "Salmon", "Salmo Salar", "Atlantic Salmon", "Pacific Salmon", "Chinook Salmon", "Oncorhynchus tshawytscha", "Coho Salmon", "Oncorhynchus kisutch", "Pink Salmon", "Oncorhynchus gorbuscha", "Sockeye Salmon", "Oncorhynchus nerka", "Chum Salmon", "Oncorhynchus keta", "Masu Salmon", "Oncorhynchus masou"],
            "Marine": ["Marine Animal", "Marine Life", "Aquatic Animal", "Oceanic Animal", "Saltwater Animal", "Marine Mammal", "Marine Reptile", "Marine Invertebrate", "Marine Fish", "Marine Amphibian", "marine ascidian", "marine invertebrate", "marine vertebrate", "marine reptile", "marine mammal", "marine fish", "marine amphibian", ],
            "Aquaculture": ["Aquaculture", "Aquafarming", "Fish Farming", "Aquatic Farming", "Aquatic Cultivation", "Aquatic Husbandry", "Aquatic Agriculture", "Aquatic Production", "Aquaculture System", "Aquaculture Facility", "Aquaculture Environment", "Aquaculture Equipment", "Culture pond", "Aquaculture Pond", "Aquaculture Tank", "Aquaculture System", "Aquaculture Facility", "Aquaculture Environment", "Aquaculture Equipment", "Aquaculture Production Setting", "Aquaculture Breeding Environment", "Aquaculture Breeding Farm", "Aquaculture Production Facility", "Aquaculture Farming Facility", "Aquaculture Breeding Setting", "Fish Farm", "Fish Hatchery", "Fish Rearing Facility", "Fish Cultivation Environment", "Fish Breeding Environment", "sand sediment of culture pond"],
            "Shark": ["Shark", "Selachimorpha", "Great White Shark", "Hammerhead Shark", "Whale Shark"],
            "Whale": ["Whale", "Cetacea", "Blue Whale", "Sperm Whale", "Humpback Whale", "Killer Whale"],
            "Porpoise": ["Porpoise", "Phocoenidae", "Harbor Porpoise", "Dall's Porpoise", "Vaquita"],
            "Dolphin": ["Dolphin", "Delphinidae", "Bottlenose Dolphin", "Orca", "Beluga Dolphin", "stenella attenuata", "tursiops truncatus", "delphinus delphis", "delphinidae", "delphinidae family", "delphinus", "tursiops", "stenella", "pontoporia blainvillei", "la plata dolphin"],
            "Seal": ["Seal", "Pinnipedia", "Harbor Seal", "Elephant Seal", "Fur Seal", "mirounga leonina", "leptonychotes weddellii", "weddell seal", "arctocephalus gazella", "antarctic fur seal", "phoca vitulina", "harbor seal", "phoca largha", "spotted seal", "phoca hispida", "ringed seal"],
            # "Seafood": ["Seafood", "Edible Marine Life", "Marine Food", "Aquatic Food", "Oceanic Food",],
            "Sea Lion": ["Sea Lion", "Otariidae", "California Sea Lion", "Steller Sea Lion", "phocarctos hookeri", ],
            "Jellyfish": ["Jellyfish", "Cnidaria", "Box Jellyfish", "Moon Jellyfish", "Man O' War"],
            "Octopus": ["Octopus", "Octopoda", "Giant Pacific Octopus", "Common Octopus", "Octopus Vulgaris"],
            "Crab": ["Crab", "Brachyura", "Blue Crab", "King Crab", "Hermit Crab"],
            "Lobster": ["Lobster", "Nephropidae", "American Lobster", "Spiny Lobster", "Rock Lobster"],
            "Shrimp": ["Shrimp", "Caridea","Pink Shrimp", "Freshwater Shrimp", "Prawn"],
            "Shellfish": ["Clam", "Bivalvia", "Soft-Shelled Clam", "Hard-Shelled Clam", "Mussel", "Scallop", "Oyster"],
            "Sea Cucumber" : ["Sea Cucumber", "Holothuroidea", "Sea Apple", "Sea Slug", "Holothurian"],
            "Sea Snail" : ["Sea snail", "sea slug", 'batillaria multiformis'],
            "Squid": ["Squid", "Teuthida", "Giant Squid", "Cuttlefish", "Mantle Squid"],
            "Starfish": ["Starfish", "Asteroidea", "Sea Star", "Echinoderm"],
            "Narwhal": ["Narwhal", "Monodon Monoceros", "Unicorn Whale", "Sea Unicorn"],
            "Orca": ["Orca", "Orcinus Orca", "Killer Whale", "Blackfish"],
            "Beluga Whale": ["Beluga Whale", "Delphinapterus Leucas", "White Whale"],
            "Walrus": ["Walrus", "Odobenus Rosmarus", "Pinniped", "Odobenidae"],
            "Manatee": ["Manatee", "Trichechus", "Sea Cow", "West Indian Manatee", "Amazonian Manatee"],
            "Miscellaneous Birds": ["Canada Goose", "Branta Canadensis", "Seagull", "Laridae Family", "Albatross", "Diomedeidae Family", "Heron", "Ardeidae Family", "Stork", "Ciconia Family", "Loon", "Gavia Immer", "Sandpiper", "Scolopacidae Family", "Gannet", "Morus Bassanus", "Great Blue Heron", "Ardea Herodias", "Black Swan", "Cygnus Atratus", "Snow Goose", "Anser Caerulescens", "Green-Winged Teal", "Anas Crecca", "Red-Crowned Crane", "Grus Japonensis", "Black-Headed Gull", "Chroicocephalus Ridibundus", "Chilean Flamingo", "Phoenicopterus Chilensis", "Great Egret", "Ardea Alba", "Spoonbill", "Platalea Species", "Belted Kingfisher", "Megaceryle Alcyon", "Swallowtail Kite", "Elanoides Forficatus", "Pacific Loon", "Gavia Pacifica", "Dendrocygna Viduata", "Whistling Duck", "Parrot", "ardea cinerea", "grey heron"],
        },
        "Companion": {
            "Others": ["Pet", "Domestic Pet", "Companion Animal", "House Pet", "Family Pet", "Pets", "Companion Animals", "House Pets", "Family Pets"],
            "Feline": ["Feline", "Cat", "Domestic Cat", "Felis Catus", "House Cat", "Kitten", "Sphynx cat", "felis silvestris catus", "domestic cat", "house cat", "kitten", "kitty", "tabby", "tabby cat", "tiger cat", "felis catus"],
            "Canine": ["Canine", "Dog", "Domestic Dog", "Canis Lupus Familiaris", "Puppy", "Poodle", "Retriever", ],
            "Rodent": ["Rodent", "Mouse", "Rat", "Gerbil", "Hamster", "Rodentia", "Gerbil", "Domestic Rat", "Laboratory mouse"],
            "Rabbit": ["Rabbit", "Leporidae", "Bunny", "Domestic Rabbit", "Hare", "Pet Rabbit"],
            "Guinea Pig": ["Guinea Pig", "Cavia Porcellus", "Cavy", "Cavy Pig", "Peruvian Guinea Pig"],
            "Hamster": ["Hamster", "Cricetinae", "Syrian Hamster", "Dwarf Hamster", "Golden Hamster"],
            "Ferret": ["Ferret", "Mustela Putorius Furo", "Domestic Ferret", "Polecat", "Mustelid"],
            "Miscellaneous Birds": ["Macaw", "Ara Macao", "Canary", "Serinus Canaria", "Budgerigar", "Melopsittacus Undulatus"]
        }
    }


    best_match_score = 0
    best_category = "Unknown"
    best_animal = "Unknown"
    
    try:
        if isinstance(fuzzy_str, str):
            logging.debug(f"Entered the str loop in find_animal_category function with string: `{fuzzy_str}`.")
            part_digit_free = re.sub(r'\d+', '', fuzzy_str)
            part = re.sub(r'[\[\]{}()]', '', part_digit_free)
            part = re.sub(r'"', '', part)
            part = re.sub(r'\'', '', part)            
            
            logging.debug(f"Normalized string: {part}")
            space_count = part.count(' ')
            if space_count > 2:
                parts = re.split(r'[-+:;/_,\s]|\band\b', part)
                parts = [part.strip() for part in parts]
            else:
                parts = [part]
            logging.debug(f"Split string: {parts}")

            for normalized_str in parts:
                normalized_str = normalized_str.lower()
                logging.debug(f"Normalized string: {normalized_str}")
                for categories_of_animals, subcategories_of_animal in animals_categorized.items():
                    for animal_type, fuzzy_animal_names in subcategories_of_animal.items():
                        for fuzzy_animal_name in fuzzy_animal_names:
                            similarity_score = fuzz.ratio(normalized_str.lower(), fuzzy_animal_name.lower())           
                            if similarity_score > best_match_score:
                                best_match_score = similarity_score
                                best_category = categories_of_animals
                                best_animal = animal_type
                                best_matched_string = fuzzy_animal_name
                                normalized_str_matched = normalized_str
                            else:
                                continue

            if best_match_score >= 80:
                logging.debug(f"Human Disease identified for string: `{fuzzy_str}`. The identified animal is `{best_animal}` in category `{best_category}` with a score of {best_match_score}. Matched string: {best_matched_string} with {normalized_str_matched}.")
                return best_category, best_animal, best_match_score
            
            
            logging.debug(f"Couldn't identify animal with non-splitted form. Trying with splitted string.")
            parts = re.split(r'[-+:;/_,\s]|\band\b', part)
            parts = [part.strip() for part in parts]

            logging.debug(f"Split string: {parts}")

            for normalized_str in parts:
                normalized_str = normalized_str.lower()
                logging.debug(f"Normalized string: {normalized_str}")
                for categories_of_animals, subcategories_of_animal in animals_categorized.items():
                    for animal_type, fuzzy_animal_names in subcategories_of_animal.items():
                        for fuzzy_animal_name in fuzzy_animal_names:
                            similarity_score = fuzz.ratio(normalized_str.lower(), fuzzy_animal_name.lower())           
                            if similarity_score > best_match_score:
                                best_match_score = similarity_score
                                best_category = categories_of_animals
                                best_animal = animal_type
                                best_matched_string = fuzzy_animal_name
                                normalized_str_matched = normalized_str
                            else:
                                continue

            if best_match_score >= 80:
                logging.debug(f"Human Disease identified for string: `{fuzzy_str}`. The identified animal is `{best_animal}` in category `{best_category}` with a score of {best_match_score}. Matched string: {best_matched_string} with {normalized_str_matched}.")
                return best_category, best_animal, best_match_score

            else:
                logging.debug(f"The string: `{fuzzy_str}` didn't match any entry from the database. Highest match was with `{best_animal}` with {best_match_score} score. Returning `Unknown`.")
                return "Unknown", "Unknown", 0
        else:   

            logging.debug(f"Given string: `{fuzzy_str}` is not a string, returning `Unknown`.")
            return "Unknown", "Unknown", 0
    
    except Exception as e:
        logging.error(f"An unexpected error occurred while finding animal category from string `{fuzzy_str}`: {e}")
        return "Unknown", "Unknown", 0

def find_animal_sample(fuzzy_str: str):
    logging.debug(f"Starting to look for animal sample source for {fuzzy_str}")

    animal_body_samples = {
        "Blood": ["Serum", "Plasma", "Whole Blood", "Blood Sample", "Venous Blood", "blood whole", "Capillary Blood", "Blood Clot", "Hematologic Sample", "Blood Serum", "Blood Culture", "Blood Smear", "Whole Blood (Aquatic)", "Blood", "blood stream", "heart blood", "blood from heart", "blood from vein", "blood from artery", "blood from capillary", "blood from clot", "blood from hematologic sample", "blood from serum", "blood from culture", "blood from smear", ],
        "Hemolymph": ["Hemolymph", "Insect Hemolymph", "Arthropod Hemolymph", "Hemolymph Sample", "Hemolymph Fluid", "Hemolymph Extract", "Hemolymph Collection", "Hemolymph Specimen", "Hemolymph Sample from Insect", "Hemolymph Sample from Arthropod"],
        "Feces": ["Cloaca","broiler feces", "cloacal swab", "fecal composite (dairy cow farm)", "fecal composite", "fecal grab", "vent swab", "Stool", "polar bear stool", "Dung", "Poop", "Droppings", "Manure", "fecal swab horse", "cattle stool", "deer rectum", "doe feces", "cattle feces", "Excrement", "excreta", "Faeces", "Rectal Sample", "Intestinal Content", "Fecal Swab", "Fecal Matter", "Fish Feces", "Feces", "Shrimp Feces", "feces sample", "pellet", "stool pellet", "chimpanze feces", "guano", "dairy calf feces", "feces from pen", "feces sample", "stool sample", "Rectal Swab", "camel feces", "camel dung", "camel excreta", "camel manure", "camel droppings", "camel stool", "cow feces", "cow dung", "cow excreta", "cow manure", "cow droppings", "cow stool", "goat feces", "goat dung", "goat excreta", "goat manure", "goat droppings", "goat stool", "sheep feces", "sheep dung", "sheep excreta", "sheep manure", "sheep droppings", "sheep stool", "buffalo feces", "buffalo dung", "buffalo excreta", "buffalo manure", "buffalo droppings", "buffalo stool", "horse feces", "horse dung", "horse excreta", "horse manure", "horse droppings", "horse stool", "donkey feces", "donkey dung", "donkey excreta", "donkey manure", "donkey droppings", "donkey stool", "yak feces", "yak dung", "yak excreta", "yak manure", "yak droppings", "yak stool", "reindeer feces", "reindeer dung", "reindeer excreta", "reindeer manure", "reindeer droppings", "reindeer stool", "moose feces", "moose dung", "moose excreta", "moose manure", "moose droppings", "moose stool", "elk feces", "elk dung", "elk excreta", "elk manure", "elk droppings", "elk stool", "deer feces", "deer dung", "deer excreta", "deer manure", "deer droppings", "deer stool", "anal swab", "anus", "stray dog feces", "white-tailed deer feces", "elephant manure", "equine feces", "horse dung", "horse manure", "horse feces", "horse stool", "rectum", "feline cat feces", "feces goat", "feces goat siberian ibex", "goat feces", "feces from goat", "goat feces", "goat dung", "goat excreta", "goat manure", "goat droppings", "goat stool", "sheep feces", "sheep dung", "sheep excreta", "sheep manure", "sheep droppings", "sheep stool", "buffalo feces", "buffalo dung", "buffalo excreta", "buffalo manure", "buffalo droppings", "buffalo stool", "animal fecal content", "a feces sample of chicken origin", "feces from chicken", "feces from pen", "feces from pig", "feces from sheep", "feces from goat", "feces from horse", "feces from cow", "feces from donkey", "feces from yak", "feces from reindeer", "feces from moose", "feces from elk", "feces from deer", "feces from camel", "feces from dog", "feces from cat", "feces from rabbit", "feces from guinea pig", "feces from hamster", "feces from ferret", "feces from chinchilla", "feces from gerbil", "feces from mouse", "feces from rat", "feces from hedgehog", "feces from ferret", "feces from rabbit", "feces from guinea pig", "feces from hamster", "feces from chinchilla", "feces from gerbil", "feces from mouse", "feces from rat", "feces from hedgehog", "a feces sample of rabbit origin", "a feces sample of guinea pig origin", "a feces sample of hamster origin", "a feces sample of ferret origin", "a feces sample of chinchilla origin", "a feces sample of gerbil origin", "a feces sample of mouse origin", "a feces sample of rat origin", "a feces sample of hedgehog origin", "fecal swab from pen", "animal stool", "animal feces", "animal excreta", "animal manure", "animal droppings", "animal dung", "animal waste", "animal excrement", "animal fecal matter", "animal fecal sample", "animal fecal swab", "animal fecal content", "Bovine feces", "bovine stool", "calf feces", "calf dung", "calf excreta", "calf manure", "calf droppings", "calf stool", "bovine feces", "bovine dung", "bovine excreta", "bovine manure", "bovine droppings", "bovine stool", "bovine fecal matter", "bovine fecal sample", "bovine fecal swab", "bovine fecal content", "bovine fecal sample", "bovine fecal swab", "bovine fecal content", "bovine fecal matter", "bovine feces sample", "bovine feces swab", "bovine feces content", "bovine feces matter", "bovine feces sample", "bovine feces swab", "bovine feces content", "bovine feces matter", "bovine slurry", "bovine manure slurry", "bovine fecal slurry", "bovine fecal matter slurry", "cattle slurry", "cattle manure slurry", "cattle fecal slurry", "cattle fecal matter slurry", "defecation", "defecate", "defecating", "defecation sample", "defecation swab", "defecation content", "defecation matter", "dung", "litter", "wet droppings"],
        "Blowhole": ["Blowhole", "Blowhole Sample", "Blowhole Swab", "Blowhole Culture", "Blowhole Secretion", "Dolphin Blowhole", "Whale Blowhole", "Blowhole Fluid", "Blowhole Mucus", "Blowhole Secretion Sample", "Blowhole Swab Sample", "Blowhole Culture Sample", "Blowhole Fluid Sample", "Blowhole Mucus Sample", "blowhole swab from dolphin", "blowhole swab from whale"],
        "Urine": ["Urine","Pee", "Urinary Sample", "Urine Sample", "Urinary Tract Sample", "Urinary Secretion", "Urine Culture", "Urine", "urine swab", "Canine Urine", "Feline Urine", "Urinary Bladder Sample", "Urinary Excretion", "Urinary Output", "Urinary Analysis", "Urine Test", "Urine Specimen", "Urine Collection", "Urine Dipstick", "Urine Strip Test", "Urine Microscopy", "Urine Culture and Sensitivity", "Urinary Tract Infection Sample", "urine sample from pen", "bladder", "mid-stream urine", "feline cat urine", "urine cystocentesis", "urine culture", "urine swab", ],
        "Saliva": ["saliva","Spit", "Oral Fluid", "Salivary Sample", "Saliva Culture", "Salivary Secretion", "Saliva", "human saliva"],
        "Skin": ["skin","Dermal Sample", "Epidermis", "Skin Scrapings", "Skin Biopsy", "Cutaneous Sample", "Skin", "Skin Lesion Sample", "Hair Follicle Sample", "Fish Skin", "Marine Mammal Skin", "hide", "bovine hide", "skin lesion", "pastern skin", "pastern", "cat skin", "chicken neck skin"],
        "Anal Gland": ["anal gland","Anal Sac", "Anal Gland Sample", "Anal Gland Fluid", "Anal Gland Secretion", "Anal Gland Culture", "Anal Gland Swab", "anal sac aspirate", "anal bag", "anal gland discharge swab", "anal gland discharge", "anal gland sinus", ""],
        "Adrenal Gland" : ["adrenal gland","Adrenal Sample", "Adrenal Gland Fluid", "Adrenal Gland Secretion", "Adrenal Gland Culture", "Adrenal Gland Swab", "adrenal gland aspirate", "adrenal gland biopsy", "adrenal gland fluid", "adrenal gland secretion"],
        "Stones" : ["Stones","Calculi", "Urinary Stones", "Kidney Stones", "Gallstones", "Bladder Stones", "Stones Sample", "Stones Analysis", "Stones Culture", "Stones Swab", "Stones", "biliary stones", "bile stone", "plus stones"],
        "Accessory Gland": ["accessory gland", ],
        "Ancient" : ["ancient sample", "ancient specimen", "ancient material", "ancient remains", "ancient tissue", "ancient bone", "ancient hair", "ancient skin", "ancient blood", "ancient feces", "ancient urine", "ancient saliva", "ancient mucus", "ancient bile", "ancient pus", "ancient semen", "ancient eggs", "Ancient", "mummified"],
        "Bile": ["bile","Gall", "Gall Bladder", "Biliary Sample", "Bile Sample", "Bile Culture", "Bile Fluid", "Bile", "bile duct swab", "biliary drain", "biliary drainage", "biloma drainage", "biloma", "bovine bile", "gall bladder bile", "gall bladder fluid", ],
        "Mucus": ["mucus","Mucous", "Mucosa", "Mucosal Sample", "Mucous Membrane", "Mucous Secretion", "Mucous Culture", "Mucous"],
        "Hair": ["hair","Fur", "Pelage", "Fur Sample", "Hair Follicle", "Hair Tuft", "Hair Swab", "Fur Biopsy", "Hair Root", "Hair Shaft", "dorsal fur"],
        "Sputum": ["Sputum", "Cough Sample", "Respiratory Secretion", "Lung Sample", "Expectorated Sputum", "Sputum Swab", "Sputum Culture", "Sputum Sample"],
        "Milk": ["Milk", "Udder Secretion", "Lactation Sample", "Milk Sample", "Mammary Secretion", "Colostrum", "Milk Culture", "Milk Analysis", "bovine milk", "dairy cow feces", "colostrum", "Pooled milk", "milk of cow", "cow milk", "milk of goat", "goat milk", "milk of sheep", "sheep milk", "milk of buffalo", "buffalo milk", "milk of camel", "camel milk", "raw milk", "raw cow milk", "raw goat milk", "raw sheep milk", "raw buffalo milk", "raw camel milk", "pasteurized milk", "pasteurized cow milk", "pasteurized goat milk", "pasteurized sheep milk", "pasteurized buffalo milk", "pasteurized camel milk", "milk of horse", "horse milk", "milk of donkey", "donkey milk", "milk of yak", "yak milk", "milk of reindeer", "reindeer milk", "milk of moose", "moose milk", "milk of deer", "deer milk", "milk of elk", "elk milk", "Cow milk fresh", "cow milk evaporated", "cow milk condensed", "cow milk skimmed", "cow milk whole", "cow milk low fat", "cow milk non fat", "cow milk full cream", "cow milk reduced fat", "cow milk lactose free", "cow milk organic", "cow milk raw", "cow milk pasteurized", "cow milk homogenized", "cow milk fortified", "Mastitis", "Mastitis Milk", "Mastitis Cow", "mastitis-milk", "dromedary milk", "milk swab", "goat milk", "high altitude milk"], 
        "Semen": ["semen","Ejaculate", "Sperm Sample", "Spermatozoa", "Semen Analysis", "Sperm Count", "Sperm Extraction", "Fish Semen", "Semen", "Sperm", "Semen", "Ejaculate", "Ejaculate Sample", "Ejaculate Isolate", "Ejaculate Specimen", "Ejaculate Fluid",],
        "Eggs": ["egg", "Yolk","Ova", "Egg Sample", "Clutch", "Reproductive Eggs", "Ovarian Sample", "Egg Fluid", "Egg Yolk", "Aquatic Eggs", "Fish Eggs", "Marine Eggs", "Egg shells", "roes", "roe", "fish roe", "caviar", "fish eggs", "fish roe", "caviar", "paulty eggs", ],
        "Pus": ["pus","abscess", "abscess drain", "abscess pus", "abscess swab", "Putulent Discharge", "Infected Discharge", "Wound Drainage", "Exudate", "Pus Culture", "Abscess Fluid", "Pus Sample", "pus swab", "Suppurative Discharge", "carbancule", "carbuncle pus", "carbuncle abscess", "carbuncle drain", "carbuncle fluid", "carbuncle culture", "carbuncle sample", "foot abscess", "foot abscess swab", "purulent material", "pustule", "pustule pus", "pustule drain", "pustule fluid", "pustule sample", "pus swab", "pus sample", "pus culture", "cheek abscess", "chest abscess", "frontal sinus abscess swab", "incision abscess", "jaw abscess swab", "leg abscess", "lung abscess swab", "lymph node abscess", "neck abscess", "parotid abscess", "perineal abscess", "peritoneal abscess", "pleural abscess", "pneumonic abscess", "retropharyngeal abscess", "sacral abscess", "submandibular abscess", "subcutaneous abscess", "surgical site abscess", "thoracic abscess", "ventral neck abscess", "wound abscess", "purulent discharge of mammary gland", "purulent fluid", "purulent material", "purulent discharge", "pus swabs", "intrauterine purulent content", "liver abscess", "liver abscess swab", "liver abscess fluid", "liver abscess culture", "liver abscess sample", "liver abscess drain", "liver abscess pus", "liver abscess aspirate", "liver abscess biopsy"],
        "Horn": ["horn","Antler", "Tusk", "Horn Sample", "Antler Sample", "Tusk Sample", "Horn Biopsy", "Antler Biopsy", "Tusk Biopsy", "Horn Swab", "Antler Swab", "Tusk Swab"],
        "Feather": ["Feather","Plumage", "Bird Feathers", "Feather Sample", "Avian Feathers", "Feather Swab", "Feather Biopsy", "Aquatic Bird Feathers"],
        "Tissues": ["Multiple tissue", "tissue", "Biopsy", "Muscle Sample", "Tissue Biopsy", "Fetal Tissue", "Tissue sample", "elbow tissue", "necropsy", "canine necroscopy", "tissue pool", "pooled tissues", "bovine multiple tissue", "muscle", "liver tissue", "tissue composite", "canine fetal tissue", "fetal tissue", "fetal tissue sample", "fetal tissue biopsy", "autopsy", "autopsy samples", "biopsy punch", "punch biopsy", "bladder biopsy", "bladder tissue", "callus tissue", "canine muscle", "cervical mass", "foot pad mass tissue", "forelimb amputation tissue", "necrotizing tissue", "open fracture biopsy", "subcutaneous tissues", "surgical site tissue", "tissue composite", "atrioventricular valve", "biopsy coronary mass", "equine necropsy", "fracture site tissue", "mass", "feline necropsy kidney", "feline tissue", "heart tissue", "incision tissue", "wound tissue", "wound tissue pool", "wound tissue sample", "wound tissue biopsy", "chicken dirty pool tissue", "chicken pool tissue", "bovine mixed tissue"],
        "Placental": ["Placental Tissue", "Placenta", "Placental Sample", "Placenta Biopsy", "Placental Cells", "placenta Biopsy", "placenta fetal side", "placenta swab"],
        "Umbilical Cord": ["Umbilical Cord", "Umbilical Tissue", "Placental Cord", "Cord Sample", "Umbilical Fluid", "Umbilical Cord Tissue", "Umbilical Cord Biopsy", "Umbilical Cord Cells", "Umblilical cord secretion", "Umbilical Cord Blood", "Cord Blood", "Blood From Umbilical Cord", "Newborn Cord Blood", "Cord Swab", "Umbilical", "umbilicus swab", "swab umbilical", "umbilical cord secretions", "umbilical secretions", "umbilical swab", "umbilicus", "umbilicus swab", "umbilical stump", "umbilical stump tissue", "umbilical stump biopsy", "omphalitis swab", "omphalitis", "omphalitis biopsy", "omphalitis tissue", "omphalitis sample", "omphalitis culture", "omphalitis fluid", "omphalitis secretion"],
        "Body Fluid": [],
        "Coelomic Cavity" : ["coelomic cavity", "coelomic fluid", "coelomic sample", "coelomic aspirate", "coelomic lavage", "coelomic wash", "coelomic swab", "coelomic culture", "coelomic biopsy", "coelomic puncture", "coelomic tap"],
        "Fetus": ["Fetus", "aborted fetus", "fetal tissue", "fetal tissue sample", "fetal tissue biopsy", "aborted donkey", "aborted horse", "aborted cow", "aborted goat", "aborted sheep", "aborted buffalo", "aborted camel", "aborted yak", "aborted reindeer", "aborted moose", "aborted elk", "aborted deer", "aborted dog", "aborted cat", "aborted rabbit", "aborted guinea pig", "aborted hamster", "aborted ferret", "aborted chinchilla", "aborted gerbil", "aborted mouse", "aborted rat", "aborted hedgehog"],
        "Fistula" : ["Fistula","Fistula Sample", "Fistula Swab", "Fistula Culture", "Fistula tract swab", "oroantral fisutla of a horse", "fistula tract", "fistula tract swab", "fistula tract culture", "fistula tract sample", "fistula tract aspirate", "fistula tract biopsy", "fistula tract fluid", "fistula tract discharge", "fistula tract drainage", "fistula tract exudate", "fistula tract pus", "fistula tract secretion", "oroantral fistula", "oroantral fistula swab", ],
        "Vaginal": ["Vaginal", "Vaginal Swab", "Vaginal Culture", "Vaginal Secretion", "Vaginal Sample", "Vaginal Fluid", "vagina", "vagina", "vaginal", "vaginal mass fluid", "vaginal swab", "vaginal culture", "vaginal secretion", "vaginal sample", "vaginal fluid"],
        "Hair Follicles": ["Hair follicles","Roots", "Hair Samples", "Dermal Follicle", "Hair Bulb", "Hair Root", "Hair Shaft"],
        "Swab": ["swab","drag swab", "equine withers", "withers", "suture swab","pin tract", "Pin tract swab", "Cotton Swab", "Culture Swab", "Oral Swab", "Skin Swab", "Swabbing Sample", "Sterile Boot Kit", "dental plaque of cows", "Choanal Swab", "Choanal", "Choana", "Choanal sample", "Choanal isolate", "swab", "swab isolate", "oral cavity", "brisket swab", "rump swab", "perineum", "perineal", "perineal swab", "udder", "udder skin", "udder swab", "udder skin swab", "udder area swab", "hip joint swab", "mammary gland", "recto-anal mucosal swab", "nipple swab", "nipple skin swab", "teat", "teat skin swab", "Teat swab", "cow udder", "bovine udder", "wounds", "teat apex", "toungue swab", "bovine toungue", "toungue", "incision swab", "incision", "incision swab", "incision dehiscence", "incision site", "incisional infection", "incision site swab", "armpit", "axillary swab", "axilla", "axillary", "axillary swab", "bladder swab plus pseudomembrane", "bladder mucosa swab", "bladder mucosa", "bladder swab", "bladder wall swab", "bladder wall", "bulla swab", "bullae", "bullae swan", "canine foot", "foot", "leg", "foot swab", "leg swab", "foot pad", "foot pad swab", "foot pad skin", "foot pad skin swab", "foot pad area swab", "foot pad area", "foot pad area swab", "foot pad area skin", "foot pad area skin swab", "foot pad area skin swab", "foot pad area skin swab", "canine mouth", "cheek", "cheek swab", "cheek pouch", "cheek pouch swab", "canine cheek", "canine cheek swab", "canine cheek pouch", "canine cheek pouch swab", "canine mouth swab", "wound swab", "wound", "wound exudate swab", "wound hock", "digit", "amputation discharge", "discharge from amputation", "discharge from amputation site", "distal right forelimb", "elbow swab", "elbow joint", "elbow joint swab", "elbow surgical site", "elbow surgical site swab", "elbow wound swab", "epididymis swab", "face", "facial", "facial swab", "facial area swab", "facial area", "foot wound", "hoof", "hoof wound", "hoof swab", "hoof area", "hoof area swab", "hoof area skin", "genital", "genital swab", "genital lesion swab", "groin", "groin swab", "inguinal", "inguinal swab", "inguinal area swab", "inguinal area", "hock", "hock swab", "interdigital furunculosis", "joint swab", "left stifle swab", "stifle swab", "abdominal incision", "lip", "mass caudal swab", "multisite swab", "neck", "neck swab", "neck area swab", "neuter site swab", "open wound", "open wound swab", "oral cavity swab", "palate", "palate swab", "Oral cavity", "oral mucosa", "oral mucosa swab", "pad swab", "Paw", "paw pad swab", "paw pad", "paw pad lesion", "paw pad lesion swab", "perianal", "pododermatitis", "pyoderma", "pyoderma swab", "Scrotal swab", "scrotum", "tail base", "tail wound", "tail base swab", "tail wound swab", "ulcer", "ulcer swab", "ulcerated mass", "ulcerated dermal mass", "urethra", "urethral", "urethral swab", "urethral area swab", "urethral area", "vulva", "vulvar", "vulvar swab", "vulvar area swab", "vulvar area", "wound site swab", "wound site", "Dental plaque", "dental plaque swab", "dental plaque sample", "dental plaque culture", "dental plaque isolate",  "aspirate swab", "aspirate", "chronic wound", "clitoris", "clitoral", "clitoral swab", "clitoral area swab", "clitoral area", "draining wound swab", "eczema", "joint swab", "left hindlimb fetlock wound", "left hock swab", "lesion swab", "lesion swab from chest", "mammary", "mass swab", "mass swab from left forelimb", "mass swab from right forelim", "Neck", "nodule swab", "omphalitis - swab", "omphalitis swab", "omphalitis swab from umbilicus", "perineal area swab", "perineal area", "perineal region swab", "perineal region", "perineum swab", "perineum area swab", "perineum area", "rectal mucosa swab", "rectal mucosa", "rectal mucosa swab from rectum", "rectum swab", "rectum area swab", "rectum area", "surgical site swab", "surgical site culture", "surgical site sample", "sub mandibular swab", "subcutaneous wound", "submandibular", "udder discharge", "urachus swab", "withers", "withers swab", "bite", "bite wound", "bite wound swab", "bite wound discharge", "deep ulceration", "e-tube site", "exudate swab", "esophagostomy tube", "e tube", "feline swab", "mass fluid", "Oral", "surgical infection site", "surgical site", "surgical site swab", "surgical site culture", "surgical site sample", "surgical site aspirate", "surgical site fluid", "fluid from surgical incision", "surgical site discharge", "surgical site drainage", "surgical site exudate", "surgical site pus", "surgical site secretion", "surgical site drainage swab", "surgical site exudate swab", "ulcerated maxilla", "footpad infection", "chicken leg", "crust swab", "catheter insertion site", "tplo site", "tplo swab", ],
        "Implant" : ["screws", "implant", "infected implant", "bone screw", "femoral pin", "femoral plate", "implant site", "implant swab", "surgical implant", "tta implant", "tplo implant", "tplo screw swab", "tplo screw", "screw surgery", "explant hardware", "tplo swab plus screw"],
        "Abdominal": ["abdominal fluid","Abdominal Fluid swab", "Abdominal Sample", "Abdominal Swab", "Abdominal Aspirate", "Abdominal Culture", "Abdominal Biopsy", "Abdominal Wash", "Abdominal Drainage", "Abdominal Exudate", "abdominal cavity", "abdominal abcess", "abdominal effusion", "abdominal tap", "peritoneum", "peritoneal fluid", "abdominal fluid","peritoneal aspirate", "peritoneal lavage", "peritoneal swab", "peritoneal culture", "peritoneal biopsy", "peritoneal puncture", "peritoneal tap", "abdominal tissue", "abdominal incision", "abdominal incision", "abdominal biopsy", "abdominal muscle", "abdominal muscle tissue", "abdominal muscle biopsy", "abdominal muscle wall", "gut", "omentum", "abdomen", "caudoventral abdomen", "digestive", "gut contents", "Gut Tissue", "fish gut", "foregut sample", "foregut", "hindgut", "hindgut sample", "human abdominal fluid", "ascites", "ascitic fluid", "ascitic fluid swab", "ascitic fluid culture", "ascitic fluid sample", "ascitic fluid aspirate", "ascitic fluid biopsy", "ascitic fluid exudate", "ascitic fluid drainage", "ascitic fluid discharge", "punctate abdominal cavity", "gizzard"],
        "Tear Fluid": ["tear", "Tears", "Lacrimal Fluid", "Tear Sample",],
        "Metagenomic/Cecal" : ['metagenomic', 'metagenomic sample', 'microbiome', 'microbiome sample', 'microbial sample', 'microbial community', 'microbial community sample', 'microbial metagenome', 'microbial metagenome sample', 'microbial community', 'microbiome', 'microbiome metagenome sample', 'microbiome community', 'microbiome community sample', 'microbiome community metagenome', 'microbiome', 'community', 'Cecal', "caecum", "equine fecal microbiome", "ceca", "cecal content", "animal cecal content",  "cecal","cecal sample","ceca", "cecal feces of duck slaughterhouse", "broiler cecum", "cecum",],
        "Draining Tract": ["Draining tract", "draining tract sample", "draining tract swab", "draining tract aspirate", "draining tract culture", "draining tract biopsy", "draining tract fluid", "draining tract exudate", "draining tract discharge", "draining tract pus", "drain tube", "Draining Tract"],
        "Bone": ["Bone","Bone Marrow", "Bone Sample", "Skeletal Tissue", "Bone Biopsy", "Bone Fragment", "Bone Swab", "Osseous Sample", "Fish Bone", "spinal cord", "carpus", "femoral pin", "tibia", "tarsus", "tarsal", "tarsal bone", "tarsal joint", "feline bone biopsy", "glenoid biopsy", "humeral head biopsy", "tympanic bulla", "vertebral bone",],
        "Cartilage" : ["cartilage","Cartilage Sample", "Cartilage Biopsy","Cartilage Fluid", "cartilaginous tissue", "cartilaginous", "chondroid", "chondroid tissue", "chondroid sample", "chondroid biopsy", "chondroid fluid", "cartilage swab", "cartilage aspirate", "cartilage culture", "cartilage biopsy", "cartilage puncture", "cartilage tap"],
        "Respiratory Samples": ["tracheal isolate","Tracheal","Trachea", "equine lower respiratory tract", "paranasal sinuses","infraorbital sinuses", "lung", "lung swab", "bovine lung", "Lower Respiratory tract","upper respiratory tract","Lung Lavage", "Airway Sample", "Bronchoalveolar Lavage", "Respiratory Tract Sample", "Tracheal Sample", "Bronchial Secretion", "Bronchial wash","nose swab", "nasopharyngeal swab", "nostril sample", "Throat Swab", "Nasal Swab",  "goat nasal swab", "pleural fluid", "pleural aspirate", "pleural lavage", "pleural swab", "pleural culture", "pleural biopsy", "pleural puncture", "pleural tap", "lung tissue", "tracheal wash", "tracheal aspirate", "tracheal lavage", "tracheal swab", "tracheal culture", "tracheal biopsy", "tracheal puncture", "tracheal tap", "bronchial aspirate", "bronchial lavage", "bronchial swab", "bronchial culture", "bronchial biopsy", "bronchial fluid", "bronchial lavage", "bronchi", "bronchus", "bronchial puncture", "airway", "anterior nares", "bal", "bronchoalveolar lavage", "Minibal", "Mini Bronchoalveolar Lavage", "Lung Lavage", "Lung Fluid Sample", "Bronchial Lavage", "Alveolar Lavage", 'balf', 'bronchoalveolar lavage fluid', 'bronchial', 'bronchial aspirate', "broncheal washing", "br wash", "canine necropsy lung", "canine lung", "lung necroscopy", "chest", "throat", "oropharynx", "air sac", "horse lung", "lung fragment", "lung tissue", "lung tissue sample", "maxillary sinus of a horse", "nasal cavity", "nasal discharge", "nasal drainage swab", "nasal exudate", "nasal flush", "nasal", "nasal mucosa", "nasal pharyngeal", "nasal pharyngeal wash", "nasal swab", "nasopharygneal lavage", "nasopharyngeal", "nasopharyngeal swab", "nose swab", "nose", "nose swab", "nostril", "pharyngeal", "pharyngeal swab", "pharyngeal wash sample", "pharynx swab", "pleura","pleural effusion", "respiratory", "respiratory sample", "retropharyngeal swab", "sinus", "sinus aspirate", "sinus cyst swab", "sinus exudate", "sinus fluid", "sinus swab", "thoracic aspirate", "thoracic fluid", "thoracic mass", "thoracic cavity", "thoracocentesis", "trach wash", "tracheal aspirate", "tracheal wash", "tracheal wash sample", "tracheobronchial aspirate", "tracheobronchial lavage", "tracheobronchial swab", "tracheobronchial culture", "tracheobronchial biopsy", "tracheobronchial puncture", "tracheobronchial tap", "ventral nasal cavity", "transedoscopic tracheal wash", "transtracheal wash fluid", "transtracheal aspirate", "transtracheal lavage", "transtracheal swab", "transtracheal culture", "transtracheal biopsy", "transtracheal puncture", "transtracheal tap", "ventral nasal cavity swab", "ventral nasal cavity aspirate", "ventral nasal cavity lavage", "ventral nasal cavity culture", "ventral nasal cavity biopsy", "ventral nasal cavity puncture", "ventral nasal cavity tap", "transtracheal wash", "TTW", "endotracheal tube", "tracheal tube", "ET tube", "feline cat lung", "frontal sinus", "mustelidae mink lung", "nasal flush", "pharynx", "dog tonsil", "nasal fossa", "sinus", "sinus tract", "sinus tract swab" , "caprine lung", "endotracheal tube site",  "lung fragment", "lung tissue", "Thoracic Fluid"],
        "Cerebrospinal Fluid": ["cerebrospinal fluid","Csf", "Spinal Fluid", "Brain Fluid", "Lumbar Puncture", "Cerebral Fluid", "Spinal Tap Sample", ],
        "Cyst" : ["Cyst Fluid", "Cyst Sample", "Cyst Aspirate", "Cyst Culture", "thyroglossal cyst", "halscysta", "cystic fluid", "cystic aspirate", "cystic sample", "cystic culture", "cystic drainage", "cystic wash", "cystic swab", "cystic biopsy", "cystic puncture", "cystic tap", "liver cyst", "ovarian cyst", "postoperative maxillary cyst", "maxillary cyst", "right inner gluteal pilonidal", "gluteal cyst", "gluteal skin cyst", "Cyst Swab", "paraprostatic cyst"],
        "Nail/Nail-like":["Claw", "Nail Sample", "Talon", "Nail Clipping", "Nail Biopsy", "Nail Swab", "Nail Clipping", "hoof", "nailbed", "nail bed", "claw", "claw sample", "claw biopsy"],
        "Intestinal": ["Gut Fluid", "Intestinal Secretion", "Intestinal Sample", "Gut Lavage", "Intestinal Aspirate", "Appendix fluid", "small intestine", "large intestine", "intestinal fluid", "colon", "jejunum", "bovine necropsy intestine", "ileum", "ruminant gastrointestinal tract", "canine colon", "canine intestine", "canine intestinal fluid", "canine intestinal content", "canine intestinal aspirate", "canine intestinal lavage", "canine intestinal swab", "canine intestinal culture", "equine intestine", "equine necropsy colon", "equine necropsy intestine", "fish intestine", "intestine tissue", "goat intestine", "large intestine goat", "goat large intestiine", "cattle intestinal", "intestine cattle", "bovine intestinal", "intestinal content", "intestinal aspirate", "intestinal lavage", "intestinal swab", "intestinal culture", "intestinal biopsy", "intestinal puncture", "intestinal tap", "intestinal tissue", "intestinal fluid", "intestinal exudate", "intestinal discharge", "intestinal drainage", "intestinal secretion", "Colon Biopsy", "colon tissue", "intestine"],
        "Scab" : ["Scab", "Scab swab"],
        "Uterine": ["Uterine Fluid", "Uterine Sample", "Endometrial Sample", "Uterine Swab", "Uterine Aspirate", "Uterine Lavage", "Uterine Culture", "Uterine Biopsy", "Endometrial Fluid", "Endometrial Swab", "Endometrial Aspirate", "Endometrial Lavage", "Endometrial Culture", "Endometrial Biopsy", "Uterine Wash", "Uterine Exudate", "Uterine Discharge", "uterine rinsing", "uterus", "endometrium", "pyometra", "equine infectious endometritis", "swab from uterus"],
        "Stomach": ["Stomach fluid","Gastric Juice", "Stomach Content", "Rumen Fluid", "Gastric Fluid", "Rumen Aspirate", "abomasum", "abomasum content", "abomasum fluid", "Rumen", "goat rumen", "Rumen liqour", "Rumen fluid", "Rumen aspirate", "Rumen content", "Rumen puncture", "Rumen tap", "vomit","Regurgitation", "Emesis", "Stomach Contents", "Vomitus", "Gastric Sample", "Vomit Sample", "Stomach" , "pyloric mass"],
        "Synovial Fluid": ["Joint","Joint Fluid","Synovial Fluid", "Joint Fluid", "Knee Fluid", "Synovial Aspirate", "Synovial Sample", "Synovial Joint Fluid", "Synovial Exudate", "stifle fluid", "joint aspirate", "joint swab", "joint culture", "joint puncture", "joint tap", "synovial fluid aspirate", "synovial fluid swab", "synovial fluid culture", "synovial fluid biopsy", "synovial fluid tap", "femorotibial joint", "fluid left elbow", "fluid on mandible", "hock joint", "hock joint fluid", "joint cavity", "femorotibial joint", "mandible", "swab from tap of coffin joint", "coffin joint fluid", "coffin joint aspirate", "coffin joint swab", "coffin joint culture", "coffin joint biopsy", "coffin joint puncture", "coffin joint tap"],        
        "Ear": ["Earwax","Cerumen", "Aural Secretion", "Ear Discharge", "Ear Swab", "Aural Fluid", "Ear Sample", "Earwax Culture", "Earwax Analysis", "Earwax Swab", "Earwax Biopsy", "Earwax Fluid", "Ear biopsy", "Auricular", "acute otitis externa","middle ear fluid", "Otitis", "ear mass swab", "ear sample", "ear exudate", "feline ear", "mass on the auricle", "myringotomy fluid", "otopyorrhea", "Ear", "auricular swab", "auricular culture", "auricular biopsy", "auricular puncture", "auricular tap", "ear swab", "ear culture", "ear biopsy", "ear puncture", "ear tap", "Ear piece"],
        "Ocular": ["ocular swab","Eye Swab", "Conjunctival Sample", "Tear Swab", "Ocular Secretion", "Eye Discharge",  "Pinkeye", "Pink Eye", "Conjunctivitis", "Conjunctival Swab", "Conjunctival Sample","Conjunctival Discharge", "ocular discharge", "cornea", "corneal swab", "corneal ulcer", "corneal ulcer swab", "ocular", "conjunctiva", "corneal swab", "eye", "eye swab", "eye discharge", "eye exudate",  "Ocular Secretion", "Ocular Fluid", "Eye Discharge", "retrobulbar",],
        "Feathers": ["Feather","Plumage", "Bird Feathers", "Feather Sample", "Avian Feathers", "Feather Swab", "Feather Biopsy", "Aquatic Bird Feathers"],
        "Whole Animal": ["Entire Animal", "Whole Body Sample", "Full Animal Specimen", "Complete Animal", "Entire Organism", "Whole Carcass", "Whole Body Collection", "Live Animal Sample", "Whole Aquatic Animal", "Whole Fish", "Dried fish", ],
        "Lymph Node" : ["lymph node", "lymph node sample", "lymph node biopsy", "lymph node necropsy","lymph node culture", "lymph node swab", "lymph node tissue", "lymph node fluid", "lymph node aspirate", "lymph node lavage", "lymph node wash", "mesentric lymph node", "submandibular lymph node", "axillary lymph node", "inguinal lymph node", "cervical lymph node", "popliteal lymph node", "bronchial lymph node", "mediastinal lymph node", "thoracic lymph node", "abdominal lymph node", "retroperitoneal lymph node", "pelvic lymph node", "mesenteric lymph node biopsies", "tracheobronchial lymph node", "laryngopharyngeal lymph node", "laryngopharyngeal lymph node", "pectoral lymph nodes", "mandibular lymph node", "retropharyngeal lymph node", "mesentery ln", "pre-scapular ln", "subiliac lymph node", "mesenteric lymph node", "subiliac lymph nodes", "bovine subiliac lymph nodes", "bovine lymph node", "intestine lymph node", "prefemoral lymph node", "liver lymph node", "calf lymph node", "lymph node draining tract", "ln", "ln fluid", "submandibular lymph node abscess", "submandibular lymph node aspirate"],
        "Meat/Organ": ["tuscan sausage", "beef from market", "animal organ", "rump swab", "raw pork sausage","Sausage", "goat meat", "som fak", "the corpse of a deer", "Salami","Poultry meat", "Squid roll","tenderized squid roll", "ground beef patties", "raw ground beef", "heart", "roasted beef", "beef trim", "pancreas", "bob veal", "testicle", "spleen", "kidney", "bovine minced meat", "minced meat", "spiced and fried sole fish", "spiced and fried fish", "undercooked meat", "bovine adipose trim", "ground beef", "uncooked meat","undercooked beef","uncooked beef","raw beef","raw chicken meat","raw meat","beef", "fish meal","poultry meat","meat", 'organ', "fish meat", "chicken breast", "wings", "chicken wings","liver","beef steak", "pork chop", "lamb chop", "veal cutlet", "mutton leg", "turkey drumstick", "duck breast", "goose leg",  "gooose meat", "quail breast", "rabbit leg", "venison steak", "elk chop", "bison roast", "boar ribs", "ostrich fillet", "kangaroo steak", "crocodile tail", "alligator leg", "iguana tail", "turtle meat", "snake meat", "frog legs", "crab meat", "lobster tail", "shrimp meat", "clam meat", "oyster meat", "mussel meat", "scallop meat", "octopus tentacle", "squid tube", "cuttlefish steak", "starfish arm", "sea urchin roe", "jellyfish tentacle", "ground beef", "ground chicken", "ground pork", "ground turkey", "trachea", "lung and intestine", "lung live pool", "lung intestinal pool", "lung and liver", "turkey jerky", "patty", "turkey patty", "beef patty", "chicken patty", "pork patty", "lamb patty", "veal patty", "mutton patty", "turkey burger", "beef burger", "chicken burger", "pork burger", "lamb burger", "veal burger", "mutton burger", "turkey sausage", "beef sausage", "chicken sausage", "pork sausage", "lamb sausage", "veal sausage", "mutton sausage", "turkey bacon", "beef bacon", "chicken bacon", "pork bacon", "lamb bacon", "veal bacon", "mutton bacon", "turkey ham", "beef ham", "chicken ham", "pork ham", "lamb ham", "veal ham", "mutton ham", "turkey salami", "beef salami", "chicken salami", "pork salami", "lamb salami", "veal salami", "mutton salami", "turkey pepperoni", "beef pepperoni", "chicken pepperoni", "pork pepperoni", "lamb pepperoni", "veal pepperoni", "mutton pepperoni", "turkey pastrami", "beef pastrami", "chicken pastrami", "pork pastrami", "lamb pastrami", "veal pastrami", "mutton pastrami", "turkey bologna", "beef bologna", "chicken bologna", "pork bologna", "lamb bologna", "veal bologna", "mutton bologna", "turkey hot dog", "beef hot dog", "chicken hot dog", "pork hot dog", "lamb hot dog", "veal hot dog", "hot dog", "mutton hot dog", "turkey sausage roll", "beef sausage roll", "chicken sausage roll", "pork sausage roll", "lamb sausage roll", "veal sausage roll", "mutton sausage roll", "turkey meatball", "beef meatball", "chicken meatball", "pork meatball", "lamb meatball", "veal meatball", "mutton meatball", "turkey nugget", "beef nugget", "chicken nugget", "pork nugget", "lamb nugget", "Buffalo meat", "Bison meat", "Venison meat", "Elk meat", "Moose meat", "Reindeer meat", "Caribou meat", "Boar meat", "Hog meat", "Pig meat", "Rabbit meat", "Hare meat", "Venison meat", "Elk meat", "Moose meat", "Bison meat", "Buffalo meat", "Caribou meat", "Reindeer meat", "Deer meat", "Antelope meat", "Goat meat", "Sheep meat", "Lamb meat", "Mutton meat", "Boar meat", "Hog meat", "Pig meat", "Bear meat", "Raccoon meat", "Opossum meat", "Squirrel meat", "Rabbit meat", "Hare meat", "Beaver meat", "Muskrat meat", "Nutria meat", "Porcupine meat", "Armadillo meat", "Alligator meat", "Crocodile meat", "Turtle meat", "Snake meat", "Frog meat", "Toad meat", "Lizard meat", "Iguana meat", "Monitor meat", "Gecko meat", "Skink meat", "Chameleon meat", "Salamander meat", "Newt meat", "Axolotl meat", "Fish meat", "Shark meat", "Ray meat", "Skate meat", "Eel meat", "Lamprey meat", "Hagfish meat", "Mollusk meat", "Clam meat", "Oyster meat", "Mussel meat", "Scallop meat", "Snail meat", "Slug meat", "Squid meat", "Octopus meat", "Cuttlefish meat", "Nautilus meat", "Cephalopod meat", "Crustacean meat", "Lobster meat", "Crab meat", "Shrimp meat", "Prawn meat", "Crayfish meat", "Barnacle meat", "blood sausage", "carcass", "bovine pre-evisceration carcass", "post-intervention beef carcass", "goat carcass", "bovine liver", "bovine kidney", "ready-to-cook beef cutlets", "strip loin", "bovine strip loin", "bovine ribeye", "bovine rib eye", "bovine rib-eye", "bovine rib eye", "beef cubed steak", "organic drumett", "beef enrichment", "gill", "fin", "fish skin", "fish fillet", "fresh crab meat", "crab meat", "Cooked crab meat", "pre-cooked crab meat", "retail deer", "spleen of pantropical spotted dolphin", "liver of pantropical spotted dolphin", "kidney of pantropical spotted dolphin", "heart of pantropical spotted dolphin", "lung of pantropical spotted dolphin", "frozen striped bass", "bass", "striped bass", "liver lung intestine", "fish sauce mash", "anchovy chrunch fish", "anchovy dried fish", "anchovy dried fish", "fillet fish", "fish fillet", "fish meat", "Fermented Freshwater fish", "fermented fish", "salted fermented fish", "salted fish", "fish roe", "fish egg", "caviar", "fish balls", "fish paste", "fish sauce", "fish oil", "fish liver oil", "fish meal", "fish powder", "fish protein concentrate", "fish protein isolate", "fish protein hydrolysate", "fish protein digest", "fish protein extract", "fish protein concentrate powder", "fish protein isolate powder", "fish protein hydrolysate powder", "fish protein digest powder", "fish protein extract powder", "Frozen fish", "frozen anchovy fish",  "frozen baby eel fih fish",  "frozen baby gold fish",  "frozen baila fish",  "frozen bangamary fish",  "frozen bangladeshi fish",  "frozen barramundi fish",  "frozen bonito fish",  "frozen broadhead fish",  "frozen caesio fish",  "frozen clarias fish",  "frozen clarias fish",  "frozen conger pike sifh fish",  "frozen corvina fish",  "frozen croaker fish",  "frozen cutlass fish",  "frozen dace fish",  "frozen devil fish",  "frozen doby fish",  "frozen dried salted cobia fish",  "frozen eel fish",  "frozen farmed headless vacuum packed bullhead fish",  "frozen feather back fish meat",  "frozen featherback fish",  "frozen featherback fish meat",  "frozen fish",  "frozen fish",  "frozen fish bone",  "frozen fish pabda",  "frozen fusilier fish",  "frozen goby fish",  "frozen golden threadfin bream fish",  "frozen grouper fish",  "frozen hilsa fish",  "frozen hilsa fish",  "frozen ladyfish fish",  "frozen long tailed anvoy fish",  "frozen lotia fish",  "frozen mackeral fish",  "frozen mackerel fish",  "frozen mackerel fish",  "frozen mahi fish",  "frozen mandarin fish",  "frozen moilla fish",  "frozen mullet fish",  "frozen pangasius fish",  "frozen pangasius fish maw",  "frozen pomfret fish",  "frozen pompano fish",  "frozen raw escolar steak fish",  "frozen raw grouper fillet fish",  "frozen raw whole baby sand goby fish",  "frozen red drum whole fish",  "frozen red nilotica whole round tilapia fish",  "frozen river tinfoil barb fish",  "frozen riverbarb fish",  "frozen rohu fish",  "frozen rohu fish",  "frozen seabass fish",  "frozen sharptooth fish",  "frozen sheet fish",  "frozen sillago fish",  "frozen sillago fish",  "frozen silver barb gutted fish",  "frozen silver fish",  "frozen silver fish",  "frozen smelt fish",  "frozen snakehead fish",  "frozen tongue fish",  "frozen walking fish",  "frozen white fish",  "frozen whole indo pacific king mackerel fish",  "frozen whole layang scad fish",  "frozen whole mackerel vacuum packed wild fish",  "frozen whole shari puti fish",  "frozen whole spiny goby fish",  "frozen whole trumpet emperor fish",  "frozen whole waliking yellow fish",  "frozen wild caught scad fish",  "frozen yellow cleaned walking fish",  "frozen yellow headless walking fish",  "frozen yellow tail fusilier fish", "frozen tuna saku aaa grade", "longaniza", "spanish sausage", "embutido", "boneless pork picnic", "minced chicken chilled", "chicken thigh chilled", "chicken drumstick chilled", "beef chicken soy composite", "boneless skinless chicken thigh", "breaded chicken breast nugget", "breast from whole chicken", "chicken broccoli and cheese", "chicken carcass rinse", "chicken carcasses showing gastrointestinal lesions", "chicken from supermarket", "chicken from wild market", "chicken frozen form", "chicken giblets", "chicken gizzards", "chicken gizzards and hearts", "chicken gizzards and liver", "chicken gizzards and necks", "chicken gizzards and feet", "chicken gizzards and wings", "chicken gizzards and thighs", "chicken gizzards and drumsticks", "chicken gizzards and backs", "chicken gizzards and livers", "chicken gizzards and hearts", "chicken gizzards and kidneys", "chicken gizzards and lungs", "chicken gizzards and spleens", "chicken gizzards and pancreases", "chicken gizzards and intestines", "chicken liver", "chicken liver with fowl typhoid", "chicken tender raw", "chicken tetrazzini", "chicken thigh and drumstick", "chicken thigh skin on", "chicken whole broiler carcass", "chicken whole cut", "chicken whole-cut in the lab", "cooked chicken thigh", "curry raw chicken skewer", "dried vegetable chicken soup", "emulsified chicken", "finished chicken", "finished chicken fat", "food pre-packaged", "fresh whole chicken carcass pluck shop", "freshly slaughtered native chicken", "freshly slaughtered silkie chicken", "frozen chicken carcass", "frozen chicken carcass", "frozen chicken part supermarket", "frozen chicken shawrma", "frozen chicken tender", "frozen chicken thigh and leg", "frozen raw chicken grind pet food", "frozen raw stuffed chicken", "frozen raw stuffed chicken", "frozen raw stuffed chicken product", "frozen raw stuffed chicken product", "frozen whole chicken carcasses", "frozen whole chicken carcasses", "mayonnaise chicken risotto", "nonintact chicken", "non-organic retail chicken", "packaged chicken", "packaged chicken broiler leg", "peanut chicken soup", "popcorn chicken", "post chill chicken carcass", "pre-chill chicken carcass", "product raw intact chicken", "product chicken", "product raw intact chicken", "raw stuffed breaded chicken", "raw stuffed chicken products", "ready-to-cook chicken", "ready-to-cook chicken cutlets", "refrigerated raw chicken thigh and leg", "retail chicken carcasse", "retail chicken gizzard", "retail chicken gizzards", "retail chicken quarter leg", "retail chicken skin", "rice stew with chicken", "rinse from chicken carcass", "roasted chicken pea powder", "roasted chicken salad", "salted chicken breast from poultry slaughterhouse", "seasoned chicken fillet", "sick organs of a dead chicken", "stewed chicken", "stuffed chicken products", "uncooked chicken tikka", "wild chicken liver avian", "myulchi jeot", "anchovy", "anchovy jeot", "jeotgal", "jeot", "anchovy snack", "animal liver", "animal origin", "beef carcass", "beef carcass trim", "beef for fajita", "beef trimmings", "boiled pork with mustard greens", "boneless beef", "boneless beef", "boneless pork", "boneless pork butt", "boneless pork shoulder", "broiler carcass", "canned fish", "canned meat", "canned seafood", "frozen seafood", "frozen meat", "frozen poultry", "frozen fish", "frozen shrimp", "frozen crab", "frozen lobster", "frozen scallops", "frozen clams", "frozen mussels", "frozen squid", "frozen octopus", "carcass of slaughterhouse", "carcass rinse", "carcass swab", "carcass chiller", "pork curry", "pork loin", "pork loin chop", "pork loin roast", "pork rib chop", "pork rib roast", "pork shoulder chop", "pork shoulder roast", "pork tenderloin", "pork tenderloin roast", "pork belly", "pork belly roast", "pork butt roast", "pork picnic roast", "pork ham", "pork ham roast", "pork leg", "pork leg roast", "pork leg chop", "pork shank", "pork shank roast", "lamb loin chop", "lamb loin roast", "lamb rib chop", "lamb rib roast", "lamb shoulder chop", "lamb shoulder roast", "lamb leg chop", "lamb leg roast", "lamb shank", "lamb shank roast", "veal loin chop", "veal loin roast", "veal rib chop", "veal rib roast", "veal shoulder chop", "veal shoulder roast", "veal leg chop", "veal leg roast", "veal shank", "veal shank roast", "chevon", "chevon chop", "chevon roast", "goat loin chop", "goat loin roast", "goat rib chop", "goat rib roast", "goat shoulder chop", "goat shoulder roast", "goat leg chop", "goat leg roast", "goat shank", "goat shank roast", "mutton loin chop", "mutton loin roast", "mutton rib chop", "mutton rib roast", "mutton shoulder chop", "mutton shoulder roast", "mutton leg chop", "mutton leg roast", "mutton shank", "mutton shank roast", "Bovine toungue", "Bovine heart", "Bovine liver", "Bovine kidney", "Bovine spleen", "Bovine tripe", "Bovine brain", "Porcine tongue", "Porcine heart", "Porcine liver", "Porcine kidney", "Porcine spleen", "Porcine tripe", "Porcine brain", "Caprine tongue", "Caprine heart", "Caprine liver", "Caprine kidney", "Caprine spleen", "Caprine tripe", "Caprine brain", "Ovine tongue", "Ovine heart", "Ovine liver", "Ovine kidney", "Ovine spleen", "Ovine tripe", "Ovine brain", "brisket swab", "carcass swab", "carcass rinse", "meat swab", "meat rinse", "organ swab", "organ rinse", "beef from market", "pork from market", "lamb from market", "goat from market", "veal from market", "mutton from market", "frozen kangaroo meat", "frozen kangaroo loin", "frozen kangaroo tail", "frozen kangaroo steak", "frozen kangaroo fillet", "chicharron with chili", "broiler fillet", "chilled raw goat", "frozen raw goat"],
        "Insect Body" : ["insect body", "insect body part", "insect body parts", "insect body sample", "insect body specimen", "insect body fluid", "insect body tissue", "insect body organ", "insect body fluid", "insect isolate", "insect body isolate", "insect body culture", "insect body swab", "insect wing", "insect leg", "insect head", "insect thorax", "insect abdomen", "insect antenna", "insect proboscis", "insect ovipositor", "insect egg", "insect larva", "insect pupa", "insect nymph", "insect adult", "insect exoskeleton", "insect integument", "insect cuticle", "insect chitin", "insect mandible", "insect maxilla", "insect labium", "insect proboscis", "insect palpus", "insect tarsus", "insect tibia", "insect femur", "insect coxa", "insect trochanter", "insect patella", "insect claw", "insect tarsal claw", "insect tarsal pad", "insect tarsal pulvilli", "insect tarsal empodium", "insect tarsal arolium", "insect tarsal arolia", "insect tarsal aroliae", ],
        "Fish Market" : ["fish market", "fish market sample", "fish market swab", "retail fish market", "fish market environment", "fish market environment sample", "fish market environmental swab", "fish market environmental sample", "whole fish sold as", "fish sold as", ], 
        "Incubator" : ["duck incubator", "avian incubator", "avian environmental incubator"],
        "Insect larva/Nymph" : ["decomposed grub","decomposed larva","insect larva", "larva", "grub", "maggot", "caterpillar", "worm", "nymph","Tenebrio molitor larva", "mealworm", "mealworm larva","Galleria mellonella larva", "waxworm", "greater wax moth larva","Lucilia sericata larva", "green bottle fly maggot", "maggot","Calliphora vomitoria larva", "bluebottle fly maggot", "blowfly maggot","Musca domestica larva", "house fly maggot","Drosophila melanogaster larva", "fruit fly larva","Manduca sexta larva", "tobacco hornworm", "hornworm","Bombyx mori larva", "silkworm", "silk moth larva","Lepidoptera larva", "butterfly caterpillar", "moth caterpillar","Culex pipiens larva", "mosquito larva", "wiggler","Aedes aegypti larva", "yellow fever mosquito larva","Anopheles gambiae larva", "malaria mosquito larva","Coleoptera larva", "beetle larva","Leptinotarsa decemlineata larva", "Colorado potato beetle larva","Harmonia axyridis larva", "harlequin ladybird larva", "ladybug larva","Coccinella septempunctata larva", "seven-spotted ladybird larva","Melolontha melolontha larva", "cockchafer grub", "white grub","Oryctes nasicornis larva", "rhinoceros beetle grub","Lucanus cervus larva", "stag beetle larva","Dorcus parallelipipedus larva", "lesser stag beetle larva","Tribolium castaneum larva", "red flour beetle larva","Sitophilus oryzae larva", "rice weevil larva","Sitophilus granarius larva", "granary weevil larva","Ephestia kuehniella larva", "Mediterranean flour moth larva","Plodia interpunctella larva", "Indian meal moth larva", "pantry moth larva","Chrysoperla carnea larva", "green lacewing larva", "aphid lion","Forficula auricularia nymph", "earwig nymph","Lasius niger larva", "black garden ant larva","Camponotus pennsylvanicus larva", "carpenter ant larva","Atta cephalotes larva", "leafcutter ant larva","Vespa crabro larva", "European hornet larva","Apis mellifera larva", "honeybee larva", "bee brood","Bombus terrestris larva", "buff-tailed bumblebee larva","Polistes dominula larva", "paper wasp larva","Xylocopa violacea larva", "violet carpenter bee larva","Gryllus bimaculatus nymph", "two-spotted cricket nymph","Acheta domesticus nymph", "house cricket nymph","Tettigonia viridissima nymph", "great green bush-cricket nymph","Blatta orientalis nymph", "oriental cockroach nymph","Periplaneta americana nymph", "American cockroach nymph","Gromphadorhina portentosa nymph", "Madagascar hissing cockroach nymph", "surface of larvae"],
        "Slaughterhouse": ['Slaughterhouse', 'abattoir', 'chicken slaughterhouse plucked skin', 'chicken slaughterhouse saw', 'chicken slaughterhouse tool', 'chicken-abbatoir-sporadic', "bleeding room floor of slaughterhouse", "chicken slaughterhouse", "chicken slaughterhouse floor", "chicken slaughterhouse equipment", "chicken slaughterhouse environment", "chicken slaughterhouse environment sample", "chicken slaughterhouse environmental swab", "chicken slaughterhouse environmental sample", "chicken slaughterhouse equipment swab", "chicken slaughterhouse equipment sample", "chicken slaughterhouse equipment culture", "chicken slaughterhouse equipment analysis", "chicken slaughterhouse equipment examination", "chicken slaughterhouse equipment testing", "chicken slaughterhouse equipment aspirate", "chicken slaughterhouse equipment biopsy", "chicken slaughterhouse equipment wash", "chicken slaughterhouse equipment drainage", "chilling room floor of slaughterhouse", "chilling room wall of slaughterhouse", "environment_(farm/slaughterhouse)", "slaughterhouse environment", "slaughterhouse environment sample", "slaughterhouse environmental swab", "slaughterhouse environmental sample", "slaughterhouse equipment swab", "slaughterhouse equipment sample", "slaughterhouse equipment culture", "slaughterhouse equipment analysis", "slaughterhouse equipment examination", "slaughterhouse equipment testing", "slaughterhouse equipment aspirate", "slaughterhouse equipment biopsy", "slaughterhouse equipment wash", "slaughterhouse equipment drainage"],
        "Insects as Food" : ["Cricket powder", "mealworm powder", "grasshopper powder", "ant powder", "beetle powder", "maggot powder", "fly larvae powder", "insect protein powder", "insect meal", "insect protein", "insect flour", "insect-based food", "edible insects", "cricket flour", "mealworm flour", "grasshopper flour", "ant flour", "beetle flour", "maggot flour", "fly larvae flour", "insect-based protein powder", "insect-based protein meal", "insect-based protein flour", "insect-based protein supplement", "insect-based protein source", "insect-based protein ingredient", "insect-based protein product", "insect-based protein food", "insect-based protein feed"],
        "Breeding/Hospital Environment Sample": ["milk filter", "farm", "pig farm", "environmental", "feed lot", "slurry tank",  "slurry sample", "slurry solids", "muck heap effluent", "heifer shed", "calf shed", "calf pen", "calf barn", "calf house", "calf hutch", "calf rearing unit", "calf rearing facility", "calf rearing environment", "calf rearing area", "calf rearing system", "calf rearing operation", "calf rearing farm", "calf rearing site", "calf rearing location", "calf rearing establishment", "calf rearing business", "calf rearing enterprise", "milk filter", "Boot Swab", "Boot Swab Kit", "poultry flock base material and dust", "colic surgical floor", "ncsu equine educational unit", "swiffer", "swiffer equine colic recovery stall under mat", "swiffer equine colic surgical suite", "swiffer equine iso stall floor", "swiffer equine isolation nurses station", "swiffer iso 5 stall","veterinary hospital", "egg farm cage shed", "carcass cleaning wipe", "Dairy equipment filters", "meat processing equipment", "equine hospital environmental", "feedlot catch basin", "clinical", "hospital", "milking machine", "colostrum feed bucket", "chicken bedding sample", "chicken house back inside", "chicken house fly trap", "chicken house fly trap inside", "chicken house fly trap outside", "chicken house front inside", "chicken house front outside", "chicken house litter", "chicken house litter back", "chicken house litter front", "chicken house manure pile", "chicken house middle", "chicken house swab", "boot covers", "bovine calf shed", "bovine calf pen", "bovine calf barn", "bovine calf house", "bovine calf hutch", "bovine calf rearing unit", "bovine calf rearing facility", "bovine calf rearing environment", "bovine calf rearing area", "bovine calf rearing system", "bovine calf rearing operation", "bovine calf rearing farm", "bovine calf rearing site", "bovine calf rearing location", "bovine calf rearing establishment", "bovine calf rearing business", "bovine calf rearing enterprise", "bovine environmental milk filter", "bovine insect fly composite", "cage surface", "arsenic contaminated cattle tick dip sites", "poultry house dust", "poultry house dust sample", "poultry house dust swab", "poultry house environment", "poultry house environment sample", "poultry house environmental swab", "poultry house litter sample", "poultry house litter swab", "poultry house manure sample", "poultry house manure swab", "poultry house floor swab", "poultry house wall swab", "poultry house ceiling swab", "poultry house equipment swab", "poultry house tool swab", "poultry house surface swab", "poultry house bedding swab", "poultry house air sample", "poultry house ventilation system swab", "poultry house heating system swab", "poultry house cooling system swab", "poultry house lighting system swab", "tick dip site", "veterinary clinic environment", "veterinary clinic environment sample", "veterinary clinic environmental swab", "veterinary hospital environment", "veterinary hospital environment sample", "veterinary hospital environmental swab", "animal hospital environment", "animal hospital environment sample", "animal hospital environmental swab", "boot swab sample", "boot swab environmental sample", "boot swab environmental swab", "boot swab farm environment", "boot swab farm environmental sample", "boot swab farm environmental swab", "bovine farm", "cattle farm", "cattle breeding farm", "cattle production farm", "cattle rearing farm", "cattle housing farm", "cattle raising farm", "cattle feeding farm", "cattle fattening farm", "cattle finishing farm", "cattle grazing farm", "cattle pasture farm", "cattle feedlot farm", "cattle feedlot environment", "cattle feedlot facility", "cattle feedlot setting", "cattle feedlot equipment", "cattle feedlot tool", "cattle feedlot surface", "cattle feedlot bedding", "cattle feedlot air", "cattle feedlot ventilation system", "cattle feedlot heating system", "cattle feedlot cooling system", "cattle feedlot lighting system", "Rumen pill", "seed lot", "Bathroom", "Nylon in Surgery", "Surgical Nylon", "Surgical Nylon Suture", "Surgical Nylon Thread", "Surgical Nylon Material", "Suture"],
        "Breeding Environment Soil" : ["soil", "farm soil", "barn soil", "pen soil", "pasture soil", "grazing soil", "feedlot soil", "manure soil", "compost soil", "bedding soil", "soil from animal housing", "soil from animal pen", "soil from animal barn", "soil from animal pasture", "soil from animal grazing area", "soil from animal feedlot", "soil from animal manure pile", "soil from animal compost pile", "soil from animal bedding area", "cattle tick dip site soil", "tick dip site soil", "dairy farm soil", "dairy farm field soil", "dairy farm pasture soil", "dairy farm grazing soil", "dairy farm feedlot soil", "dairy farm manure soil", "dairy farm compost soil", "dairy farm bedding soil", "swine farm soil", "swine farm field soil", "swine farm pasture soil", "swine farm grazing soil", "swine farm feedlot soil", "swine farm manure soil", "swine farm compost soil", "swine farm bedding soil", "poultry farm soil", "poultry farm field soil", "poultry farm pasture soil", "poultry farm grazing soil", "poultry farm feedlot soil", "poultry farm manure soil", "poultry farm compost soil", "poultry farm bedding soil", "organic dairy farm soil", "organic dairy farm field soil", "organic dairy farm pasture soil", "organic dairy farm grazing soil", "organic dairy farm feedlot soil", "organic dairy farm manure soil", "organic dairy farm compost soil", "organic dairy farm bedding soil", "organic swine farm soil", "organic swine farm field soil", "organic swine farm pasture soil", "organic swine farm grazing soil", "organic swine farm feedlot soil", "organic swine farm manure soil", "organic swine farm compost soil", "organic swine farm bedding soil", "organic poultry farm soil", "organic poultry farm field soil", "organic poultry farm pasture soil", "organic poultry farm grazing soil", "organic poultry farm feedlot soil", "organic poultry farm manure soil", "organic poultry farm compost soil", "organic poultry farm bedding soil", "soil adjecent to bovine farm", "soil adjacent to bovine farm", "soil near bovine farm", "soil surrounding bovine farm", "soil around bovine farm", "soil beside bovine farm", "soil close to bovine farm", "soil in proximity to bovine farm", "soil at bovine farm", "soil on bovine farm", "soil within bovine farm", "soil next to bovine farm", "soil adjecent to cattle farm", "soil adjacent to cattle farm", "soil near cattle farm", "soil surrounding cattle farm", "soil around cattle farm", "soil beside cattle farm", "soil close to cattle farm", "soil in proximity to cattle farm", "soil at cattle farm", "soil on cattle farm", "soil within cattle farm", "soil next to cattle farm", "soil adjecent to swine farm", "soil adjacent to swine farm", "soil near swine farm", "soil surrounding swine farm", "soil around swine farm", "soil beside swine farm", "soil close to swine farm", "soil in proximity to swine farm", "soil at swine farm", "soil on swine farm", "soil within swine farm", "soil next to swine farm", "soil adjecent to poultry farm", "soil adjacent to poultry farm", "soil near poultry farm", "soil surrounding poultry farm", "soil around poultry farm", "soil beside poultry farm", "soil close to poultry farm", "soil in proximity to poultry farm", "soil at poultry farm", "soil on poultry farm", "soil within poultry farm", "soil next to poultry farm", "soil contaminated with killed cattle", "soil in kangaroo area of zoo"],
        "Breeding Environment Water" : ["water", "farm water", "barn water", "pen water", "pasture water", "grazing water", "feedlot water", "manure water", "compost water", "bedding water", "water from animal housing", "water from animal pen", "water from animal barn", "water from animal pasture", "water from animal grazing area", "water from animal feedlot", "water from animal manure pile", "water from animal compost pile", "water from animal bedding area", "cattle tick dip site water", "tick dip site water", "dairy farm water", "dairy farm field water", "dairy farm pasture water", "dairy farm grazing water", "dairy farm feedlot water", "trough water", "water from tanks in pens", "trough water sample", "water from tanks in pens sample", "dairy farm manure water", "dairy farm compost water", "dairy farm bedding water", "swine farm water", "swine farm field water", "swine farm pasture water", "swine farm grazing water", "swine farm feedlot water", "swine farm manure water", "swine farm compost water", "swine farm bedding water", "poultry farm water", "poultry farm field water", "poultry farm pasture water", "poultry farm grazing water", "poultry farm feedlot water", "poultry farm manure water", "poultry farm compost water", "poultry farm bedding water", "organic dairy farm water", "organic dairy farm field water", "organic dairy farm pasture water", "organic dairy farm grazing water", "organic dairy farm feedlot water", "organic dairy farm manure water", "organic dairy farm compost water", "organic dairy farm bedding water", "organic swine farm water", "organic swine farm field water", "organic swine farm pasture water", "organic swine farm grazing water", "organic swine farm feedlot water", "organic swine farm manure water", "organic swine farm compost water", "organic swine farm bedding water", "organic poultry farm water", "organic poultry farm field water", "organic poultry farm pasture water", "organic poultry farm grazing water", "organic poultry farm feedlot water", "organic poultry farm manure water", "organic poultry farm compost water", "organic poultry farm bedding water", "water adjecent to bovine farm", "water adjacent to bovine farm", "water near bovine farm", "water surrounding bovine farm", "water around bovine farm", "water beside bovine farm", "water close to bovine farm", "water in proximity to bovine farm", "water at bovine farm", "water on bovine farm", "water within bovine farm", "water next to bovine farm", "maturation tank water", "spawning tank water"],
        "Catheter": ["Catheter Fluid", "Catheter Fluid Collection", "Catheter Fluid Sample", "Catheter Fluid Isolate", "Catheter Fluid Specimen", "Catheter liquid", "cvp", "central line", "cvc", "catheter", "catheter tip", "catheter fluid", "catheter drainage", "catheter aspirate", "central venous catheter", "gastric tube", "central vein catheter", "peripheral venous catheter", "peripheral intravenous catheter", "peripheral venous line", "peripheral intravenous line", "peripheral venous access device", "peripheral intravenous access device", "peripheral venous cannula", "peripheral intravenous cannula", "peripheral venous catheterization", "peripheral intravenous catheterization", "peripheral venous access procedure", "peripheral intravenous access procedure", "peripheral venous cannulation", "peripheral intravenous cannulation", "peripheral venous access technique", "peripheral intravenous access technique", "indwelling suction", "jugular catheter", "ureteral stent", ],
        "Uncategorized Fluid" : ["Fluid", "Fluid Sample", "Fluid Isolate", "Fluid Specimen", "Fluid Collection", "Fluid Culture", "Fluid Analysis", "Fluid Examination", "Fluid Testing", "Fluid Swab", "Fluid Aspirate", "Fluid Biopsy", "Fluid Wash", "Fluid Drainage", "Fluid Discharge", "Fluid Exudate", "Bodily fluid", "shoulder aspirate", "swelling aspirate", "Body Fluid","Body Fluid Sample", "Body Fluid Culture", "Body Fluid Analysis", "Body Fluid Aspirate", "Body Fluid Swab", "Body Fluid", "antebrachium fine needle aspiration", "aspirate", "cervical fluid", "Penial discharge", "preputial discharge", "preputial discharge swab", "preputial fluid", "preputial fluid swab", "seroma aspirate", "seroma", "surgical drain", "renal pelvic fluid", "body fluid from infection", "aspirate swab",],
        "Penile/Prostatic" : ["Penial Discharge", "Penis", "Penile Discharge", "preputial discharge", "preputial discharge swab", "preputial fluid", "preputial fluid swab", "prostatic fluid", "prostatic wash", "prostatoc wash fluid", "prepuce", "foreskin", "prepuce swab", "prepuce area swab", "prepuce area", "preputial discharge", "preputial discharge swab", "preputial swab", "preputial wash", "Prostatic cyst", "Prostatic Abscess"],
        "Guttural Pouch" : ["gluttural pouch", "guteral pouch", "gutt pouch", "gutteral pouch contents", "gutteral pouch fluid", "gutteral pouch fluid", "gutteral pouch aspirate", "guttoral pouch lavage", "guttural pouch flush", "guttural pouch lavage fluid", "guttural pouch material", "right guttural pouch"],
        "Pericardial": ["Pericardial Fluid", "Pericardial Effusion", "Heart Fluid", "Pericardial Sample", "Fluid In Pericardium", "Pericardial Aspirate", "Pericardial Drainage", "pericardial effusion", "pericardial aspirate", "pericardium",],
        "Aquacultured" : ["Aquacultured", "Aquaculture", "Aquaculture Sample", "Aquaculture Isolate", "Aquaculture Specimen", "Aquaculture Fluid", "Aquaculture farm", "aquaculture sourced", "aquaculture facility", "Aquaculture water source", "Aquaculture Water", "Fish farm", "water of a fish culture pond", "fish culture pond"],
        "Nest" : ["Nest Material", "nest component", "nest", "ant nest", "wasp nest", ],
    }

    best_match_score = 0
    best_sample = "Unknown"

    try:
        if isinstance(fuzzy_str, str):
            logging.debug(f"Entered the str loop in find_animal_sample function with string: `{fuzzy_str}`.")
            
            part_digit_free = re.sub(r'\d+', '', fuzzy_str)
            part = re.sub(r'[\[\]{}()]', '', part_digit_free)            
            part = re.sub(r'"', '', part)
            part = re.sub(r'\'', '', part)
            logging.debug(f"Normalized string: {part}")
            space_count = part.count(' ')
            if space_count > 2:
                parts = re.split(r'[-+:;/_,\s]|\band\b', part)
                parts = [part.strip() for part in parts]
            else:
                parts = [part]
            logging.debug(f"Split string: {parts}")

            for normalized_str in parts:
                normalized_str = normalized_str.lower()
                logging.debug(f"Normalized string: {normalized_str}")
            
                for sample, sample_fuzzy_words in animal_body_samples.items():
                    for sample_fuzzy_word in sample_fuzzy_words:
                        similarity_score = fuzz.ratio(normalized_str.lower(), sample_fuzzy_word.lower())
                        if similarity_score > best_match_score:
                            best_match_score = similarity_score
                            best_sample = sample
                            best_matched_string = sample_fuzzy_word
                            normalized_str_matched = normalized_str
                        else:
                            continue

            if best_match_score >= 80:
                logging.debug(f"Animal sample source identified: {best_sample} with {best_match_score} score. Matched string: {best_matched_string} with {normalized_str_matched}.")
                return best_sample, best_match_score
            else:
                logging.debug(f"The string: `{fuzzy_str}` didn't match any entry from the database. Highest match was with `{best_sample}` with {best_match_score} score. Returning `Unknown`.")
                return "Unknown", 0
        else:   
            logging.debug(f"Given string: `{fuzzy_str}` is not a string, returning `Unknown`.")
            return "Unknown", 0
    
    except Exception as e:
        logging.error(f"An unexpected error occurred while finding animal sample for string `{fuzzy_str}`: {e}")
        return "Unknown", 0

def find_environmental_source(fuzzy_str):
    logging.debug(f"Starting to find the environmental source from str {fuzzy_str}")
    
    environmental_sphere = {
        "Unknown": {
            "Unknown": ['Environment', 'env', 'environmental', "environmental isolate", "environmental sample on unit"]
        },
        "Food": {
            "Baby Food" : ["Formula Feed", "Formula milk", "Baby formula"],
            "Dairy": [
                "Cream products","cream","Milk", "Tankmilk", "Dairymilk", "Stored Milk", "Freshmilk", "Skimmedmilk", "Goatmilk", "Cheese", "Butter", "Yogurt", "Cream", "Ice Cream","Milk Powder", "Evaporated Milk", "Condensed Milk", "Cottage Cheese", "Paneer", "Sour Cream", "Kefir", "Whey", "Fermented Milk", "Lactose-Free Milk", "Dairy-Free Milk", "Dairy Farm", "Dairy Product", "milk powder", "dried milk powder", "buttermilk powder", "home-made dairy foods", "home made dairy", "skim milk clt", "Dairy", "rind of surface ripened cheese", "goat cheese", "whey culture", "fat milk", "milk fat", "milk protein", "milk protein concentrate", "milk protein isolate", "milk permeate", "milk solids", "milk solids not fat", "non-fat dry milk", "non-fat milk solids", "skim milk powder", "skimmed milk powder", "whole milk powder", "whole milk solids", "whole milk solids not fat", "artisanal dairy product", "blended whey", "caciotta cheese", "caciotta", "caciotta cheese rind", "caciotta cheese surface", "caciotta cheese culture", "caciotta cheese sample", "caciotta cheese isolate", "ricotta cheese", "ricotta", "ricotta cheese rind", "ricotta cheese surface", "ricotta cheese culture", "ricotta cheese sample", "ricotta cheese isolate", "ricotta cheese product", "ricotta cheese product sample", "ricotta cheese product isolate", "pecorino cheese", "pecorino", "pecorino cheese rind", "pecorino cheese surface", "pecorino cheese culture", "pecorino cheese sample", "pecorino cheese isolate", "pecorino cheese product", "pecorino cheese product sample", "pecorino cheese product isolate", "blue cheese", "blue cheese rind", "blue cheese surface", "blue cheese culture", "blue cheese sample", "blue cheese isolate", "blue cheese product", "blue cheese product sample", "blue cheese product isolate", "brie cheese", "brie", "brie cheese rind", "brie cheese surface", "brie cheese culture", "brie cheese sample", "brie cheese isolate", "brie cheese product", "brie cheese product sample", "brie cheese product isolate", "parmesan cheese", "parmesan", "parmesan cheese rind", "parmesan cheese surface", "parmesan cheese culture", "parmesan cheese sample", "parmesan cheese isolate", "parmesan cheese product", "parmesan cheese product sample", "parmesan cheese product isolate", "cheddar cheese", "cheddar", "cheddar cheese rind", "cheddar cheese surface", "cheddar cheese culture", "cheddar cheese sample", "cheddar cheese isolate", "cheddar cheese product", "cheddar cheese product sample", "cheddar cheese product isolate", "gouda cheese", "gouda", "gouda cheese rind", "gouda cheese surface", "gouda cheese culture", "gouda cheese sample", "gouda cheese isolate", "gouda cheese product", "gouda cheese product sample", "gouda cheese product isolate", "feta cheese", "feta", "feta cheese rind", "feta cheese surface", "feta cheese culture", "feta cheese sample", "feta cheese isolate", "feta cheese product", "feta cheese product sample", "feta cheese product isolate", "mozzarella cheese", "mozzarella", "mozzarella cheese rind", "mozzarella cheese surface", "mozzarella cheese culture", "mozzarella cheese sample", "mozzarella cheese isolate", "mozzarella cheese product", "mozzarella cheese product sample", "mozzarella cheese product isolate", "ricotta salata cheese", "ricotta salata", "ricotta salata cheese rind", "ricotta salata cheese surface", "ricotta salata cheese culture", "ricotta salata cheese sample", "ricotta salata cheese isolate", "ricotta salata cheese product", "ricotta salata cheese product sample", "ricotta salata cheese product isolate", "cheese whey", "chinese fermented milk", "fermented milk", "fermented milk product", "fermented milk sample", "fermented milk isolate", "fermented milk culture", "fermented milk product sample", "fermented milk product isolate",  
                ],

            "Contaminated Dairy": [
                "Sour milk", "Spoiled milk", "Curdled milk", "Expired milk", "Moldy cheese", "Rotten butter", "Fermented butter", "Spoiled yogurt", "Foul cream", "Off cheese", "Staleice cream", "Sour kefir", "Expired paneer", "Rotten eggs", "Moldy cream", "Sour whey", "Curdled butter",
                ],
            
            "Contaminated Fruits": [
                "Rottenfruit", "Overripebananas", "Soggyapples", "Foulberries", "Rottenavocados", "Moldygrapes"
                ],
            
            "Contaminated Vegetables": [
                "cos lettuce", "romaine lettuce", "red oak lettuce","tango lettuce","shredded lettuce","red leaf lettuce","loose leaf lettuce","organic red romaine lettuce","lettuce","Soggy vegetables", "Wilted lettuce", "Spoiled potatoes", "Rotted onions", "Edible Leaf", "Vegetable Leaf", "Spinach", "Lettuce", "Kale", "Cabbage", "Broccoli", "Cauliflower", "Foul garlic", "Moldycarrots", "Spoiledspinach", "Wiltedtomatoes", "Spoiledzucchini", "Moldypeppers", "Soggymushrooms", "rotten cucumbers", "Celery stalk", "Celery leaf", "romaine", "kale", "scarlet kale", "chard", "red chard", "green chard", "Edible Leaf", "yam ring", "taro basket", "taro", "yam", "rotting potato", "arthicoke", "artichoke", "artichoke heart", "artichoke leaf", "artichoke stem", "artichoke bud", "artichoke flower", "artichoke root", "artichoke tuber", "artichoke stalk", "artichoke petal", "artichoke bract", "artichoke scale", "artichoke leaf scale", "artichoke leaf bract", "artichoke leaf petal", "artichoke leaf bud", "celery", "celery microgreens", 
                ],
            
            "Contaminated Nuts": [
                "macadamia nut","Pistachios", "Almonds", "Cashews", "Peanuts", "Walnuts", "Hazelnuts", "Macadamia Nuts", "Pecans", "Brazil Nuts", "Pine Nuts", "Chestnuts", "Pistachio Nuts", "Almond Nuts", "Cashew Nuts", "Peanut Nuts", "Walnut Nuts", "Hazelnut Nuts", "Macadamia Nuts", "Pecan Nuts", "nut raisin blend trail mixture", "raisin", "almond powder", "peanut kernel", "nut kernel", "almond kernel", "cashew kernel", "walnut kernel", "pecan kernel", "hazelnut kernel", "macadamia kernel", "pistachio kernel", "brazil nut kernel", "pine nut kernel", "chestnut kernel", "nut butter", "almond butter", "peanut butter", "cashew butter", "walnut butter", "pecan butter", "hazelnut butter", "macadamia butter", "pistachio butter", "brazil nut butter", "pine nut butter", "chestnut butter", "cashew piece", "walnut piece", "pecan piece", "hazelnut piece", "macadamia piece", "pistachio piece", "brazil nut piece", "pine nut piece", "chestnut piece", "almond piece", "peanut piece", "nut piece", "nut blend", "nut mix", "nut assortment", "nut medley", "nut combination", "nut variety", "nut selection", "nut collection", "nut platter", "nut tray", "nut bowl", "nut basket", "nut jar", "nut bag", "nut pouch"
                ],
            
            "Contaminated Herbs/Spices": [
                "sage","dried sage","dehydrated garlic","dried garlic","garlic powder","garlic","peppercorn", "cinnamon", "cardamom", "coriander", "basil", "ocimum basilicum", "cilantro", "anise spice", "star anise", "neem powder", "mint leaf", "peppermint", "spearmint", "oregano", "thyme", "rosemary", "sage", "parsley", "dill", "chive", "chervil", "tarragon", "marjoram", "bay leaf", "bay laurel", "cumin", "cumin seed", "cumin powder", "turmeric", "turmeric powder", "ginger", "ginger powder", "ginger root", "garlic", "garlic powder", "garlic clove", "onion", "onion powder", "onion flake", "green cardamom", "rotten onion bulb", "chili powder", "red chili powder", "milk masala spice mixture powder", "red paprika powder", "mint leaf powder", "mint powder", "pepper powder", "black pepper powder", "white pepper powder", "cayenne pepper powder", "chili flakes", "chili paste", "chili sauce", "chili oil", "chili powder mix", "chili seasoning", "chili blend", "chili rub", "chili marinade", "chili dressing", "chili dip", "elephant ginger", "fish masala", "spice mixture for lahori fish", "spices and herbs", "spices and salt-ground black pepper", "green pepper", "black cumin", "ground hot chili", "hot chili", "ground chili", "chili", "dehydrate small chili", "small dried chili", "ground medium hot molido chili", "ground red chili", "curry spice mixture chicken ginger masala", "anaheim pepper", "jalapeno pepper", "habanero pepper", "serrano pepper", "banana pepper", "bell pepper", "poblano pepper", "cayenne pepper", "chipotle pepper", "green chili", "red chili", "yellow chili", "orange chili", "purple chili", "white chili", "black chili", "brown chili", "chili pod", "chili seed", "chili plant", "chili bush", "chili tree", "chili flower", "barbeque seasoning", "barbeque spice mix", "barbeque rub", "Basil", "Parsley", "Oregano", "caribe pepper", "cayenne pepper", "chili powder", "chili powder super hot", "chili thai hot pepper", 
                ],


            "Animal Feed" : ['bonemeal', 'fishmeal','poultry feather meal','soybean hull', 'cotton seed hull', 'hull pellets','blood meal','feed pellet','raw feed','chicken feed', 'farm feed', 'animal feed', 'poultry feed', 'livestock feed', 'cattle feed', 'pig feed', 'horse feed', 'goat feed', 'sheep feed', 'rabbit feed', 'fish feed', 'aquatic feed', 'pet feed', 'dog feed', 'cat feed', 'bird feed', 'wildlife feed', 'zoo feed', 'aquarium feed', 'insect feed', 'rodent feed', 'dog food', 'cat food', 'pet food', 'animal food', 'livestock food', 'poultry food', 'fish food', 'wildlife food', 'zoo food', 'aquarium food', 'insect food', 'rodent food', 'dog treat', 'cat treat', 'pet treat', 'animal treat', 'livestock treat', 'poultry treat', 'fish treat', 'aquarium treat', 'insect treat', 'rodent treat', 'silage', 'rawhide moccasin dog treat chews', 'rawhide dog chews', 'rawhide dog treat', 'rawhide dog bone', 'rawhide dog bone chews', 'rawhide dog bone treats', 'dog chews', 'animal feed dried spinach powder', "cat food raw chicken", "daily vitamin for cat and kitten", "finished chubs of kitten grind", "kitten grind", "raw cat food chicken", "raw cat food duck", "raw cat food tuna", "raw cat food", "fish meal", "fish meal feed", "fish meal for animal feed", "sunflower meal pellets", "frog food", "infant food", "wellness petite treat lamb apple and cinnamon", "wellness petite treat", "wellness petite", "frozen animal feed adult mouse", "frozen animal feed baby mouse", "frozen animal feed feeder mouse", "frozen animal feed young mouse", "frozen feeder mouse", "frozen fuzzie feeder mouse", "frozen pet food hopper mouse", "frozen pinkie feeder mouse", "chicken from pet food", "chicken jerky dog treat", "chicken jerky pet treat", "chicken quail blend for dog food", "pet food with turkey chicken whitefish", "pet treat chicken jerky", "raw dog food chicken organic vegetable", "turkey chicken formula cat kitten food", "raw dehydrated chicken supreme dog food", "beef strip pet treat", "canola meal feed bulk", "canola meal feed", "canola meal animal feed", "canola meal for animal feed", "canola meal for livestock feed", "canola meal for poultry feed", "canola meal for fish feed", "canola meal for pet food", "canola meal for dog food", "canola meal for cat food", "canola meal for wildlife feed", "canola meal for zoo feed", "canola meal for aquarium feed", "canola meal for insect feed", "canola meal for rodent feed", "canola meal for animal food", "canola meal for livestock food", "canola meal for poultry food", "canola meal for fish food", "canola meal for pet food", "water buffalo dried meat dog chews"],

            "Forage" : [
                "Forage", "forages", 
                ],

            "Snacks": [
                "spring roll", "Fried snack food", "Fried snack pack", "Fried snakc", "snack", "snack pack", "snack size", "snack platter", "snack tray", "snack box", "snack bag", "snack pouch", "snack cup", "snack bowl", "snack plate", "snack stick", "corn chips snack", "chip coconut chips", "chip", "chips", "cocktail snacks chicken flavor cashew shaped biscuits", "nachos", "trail mix", "snacks", "snack food", "snack mix", "snack bar", "almond flatbread snack badam lachha", "bhel spicy snack mixture", "corn flakes snack", "corn flakes snack", "corn chips snack", "corn chips", "potato chips", "tortilla chips", "pretzels", "popcorn", "rice cakes snack", "rice crackers snack", "rice crackers", "snack crackers", "snack cookies", "snack biscuits", "snack bars", "granola bars", "energy bars", "protein bars", "fruit snacks", "fruit leather snacks", "fruit roll-ups snacks", "fruit roll-ups", "fruit snacks mix", "fruit snacks variety pack", "fruit snacks assortment pack"
            ],
            
            "Soy-Based & Fermented Soy": [
                "soybean", "hydrolyzed vegetable protein", "HVP", "Meju", "fermented soy product", "fermented soy food", "fermented soy product sample", "fermented soy food sample", "fermented soy isolate", "fermented soy specimen", "fermented sufu", "fermented tofu", "fermented bean curd", "fermented soy", "chungkukjang", "chungkukjang powder fermentation","chungkukjang powder culture", "chungkukjang powder starter culture", "korean fermented soybean food", "tofu", "spoiled tofu", "sufu", "fermented tofu", "spoiled tofu",
            ],

            "Fermented/Cultured Foods": [
                "Kombucha", "kimchi", "traditional fermented food in korea", "home-made water kefir", "sufu", "chungkukjang", "commercial cultured foods", "cultured food", "beer contaminant", "korean traditional fermented food", "pickled food", "blueberry fermentation", "traditional korean fermented vegetable", "fermented vegetable", "fermented food", "fermented food sample", "idli batter", "dosa batter", "fermented food culture", "fermented food isolate", "fermented food specimen", "fermented food swab",
            ],

            "Mushrooms": [
                "white mushroom", "mushroom", "cultivated mushroom", "shiitake", "enoki", "portobello", "oyster mushroom"
            ],

            "Ready-to-Eat / Processed Foods": [
                "ready to eat food", "pre-made food", "readymade food", "Processed Food", "Preserved Food", "cooked peanuts from wild market", "burger", "hot dog", "sausage", "processed meat", "processed food sample", "processed food isolate", "processed food specimen", "processed food culture", "processed food swab", "processed food aspirate", "processed food biopsy", "processed food wash", "processed food drainage", "processed food discharge", "processed food exudate", "ready to eat meal", "ready to eat dish", "ready to eat product", "ready to eat item", "ready to eat cuisine"
            ],

            "Pasta & Noodles": [
                "Pasta", "dried noodle", "vegetarian instant noodle", "vegetarian noodle", "vegetarian ramen", "vegetarian udon", "vegetarian soba", "vegetarian pasta", "vegetarian spaghetti", "vegetarian macaroni", "vegetarian fusilli", "vegetarian penne", "vegetarian rotini", "vegetarian farfalle", "vegetarian orecchiette", "vegetarian gnocchi", "vegetarian ravioli", "vegetarian tortellini", "vegetarian lasagna", "vegetarian fettuccine", "vegetarian linguine", "multigrain hakka noodle", "multigrain noodle", "multigrain ramen", "multigrain udon", "multigrain soba", "multigrain pasta", "multigrain spaghetti", "multigrain macaroni", "multigrain fusilli", "multigrain penne", "multigrain rotini", "multigrain farfalle", "multigrain orecchiette", "multigrain gnocchi", "multigrain ravioli", "multigrain tortellini", "multigrain lasagna", "multigrain fettuccine", "multigrain linguine", "pumiao rice noodle","instant noodle", "ramen", "glass noodle", "rice vermicelli", "ramen noodle", "udon noodle", "soba noodle", "egg noodle", "rice noodle", "wheat noodle","whole wheat noodle", "spaghetti noodle", "fettuccine noodle", "linguine noodle", "macaroni noodle", "fusilli noodle", "penne noodle", "rotini noodle", "farfalle noodle", "orecchiette noodle", "gnocchi noodle", "ravioli noodle", "tortellini noodle", "lasagna noodle", "pasta noodle", "noodle soup", "instant ramen", "instant udon", "instant soba", "instant pasta", "instant spaghetti", "instant macaroni", "instant fusilli", "instant penne", "instant rotini", "instant farfalle", "instant orecchiette", "instant gnocchi", "instant ravioli", "instant tortellini", "penne", "fusilli", "rotini", "farfalle", "orecchiette", "gnocchi", "ravioli", "tortellini", "lasagna", "fettuccine", "linguine"
            ],

            "Grains & Flours": [
                "organic buckwheat", "organic buckwheat grain", "organic buckwheat flakes", "organic buckwheat groats", "organic buckwheat cereal", "organic buckwheat noodles", "organic buckwheat bread", "organic buckwheat pancake mix", "organic buckwheat granola", "organic buckwheat muesli", "organic buckwheat porridge", "organic buckwheat flour", "organic buckwheat flour mix", "organic buckwheat flour blend", "Fonio", "quinoa", "millet", "sorghum", "barley", "oat flour", "rice flour", "amaranth flour", "spelt flour", "rye flour", "buckwheat flour", "corn flour", "potato flour", "tapioca flour", "arrowroot flour", "almond flour", "coconut flour", "chickpea flour", "besan", "gram flour", "lentil flour", "pea flour", "soybean flour", "wheat flour", "whole wheat flour", "all-purpose flour", "self-rising flour", "bread flour", "cake flour", "pastry flour", "gluten-free flour", "gluten-free all-purpose flour", "gluten-free bread flour", "gluten-free cake flour", "gluten-free pastry flour", "gluten-free oat flour", "gluten-free rice flour", "gluten-free almond flour", "gluten-free coconut flour", "gluten-free chickpea flour", "gluten-free lentil flour", "gluten-free pea flour", "gluten-free soybean flour", "gluten-free buckwheat flour", "gluten-free corn flour", "gluten-free potato flour", "gluten-free tapioca flour", "gluten-free arrowroot flour", "gluten-free spelt flour", "gluten-free rye flour", "gluten-free amaranth flour", "gluten-free quinoa flour", "gluten-free millet flour", "gluten-free sorghum flour", "gluten-free barley flour", "gluten-free fonio flour", "gluten-free teff flour", "gluten-free cassava flour", "gluten-free plantain flour", "gluten-free sweet potato flour", "gluten-free yam flour", "gluten-free taro flour", "gluten-free pumpkin flour", "gluten-free squash flour", "gluten-free zucchini flour", "gluten-free carrot flour", "gluten-free beetroot flour", "gluten-free parsnip flour", "gluten-free turnip flour", "taro flour", "yam flour", "sweet potato flour", "plantain flour", "cassava flour", "arrowroot starch", "potato starch", "tapioca starch", "cornstarch", "rice starch", "wheat starch", "soybean starch", "pea starch", "lentil starch", "chickpea starch", "amaranth starch", "quinoa starch", "millet starch", "sorghum starch", "barley starch", "oat starch", "buckwheat starch", "spelt starch", "rye starch", "fonio starch", "teff starch", "cassava starch", "plantain starch", "sweet potato starch", "taro starch", "pumpkin starch", "squash starch", "zucchini starch", "carrot starch", "beetroot starch", "parsnip starch", "turnip starch", "potato flour", "tapioca flour", "arrowroot flour", "almond flour", "coconut flour", "chickpea flour", "lentil flour", "pea flour", "soybean flour", "wheat flour", "whole wheat flour", "all-purpose flour", "self-rising flour", "bread flour", "cake flour", "pastry flour", "fonio", "wheat", "rice", "barley", "oat", "rye", "sorghum", "millet", "quinoa", "amaranth", "buckwheat", "spelt", "teff", "wild rice", "cereal grains", "rice seed", "jasmine rice", "oryza sativa", "Wheat Seed", "Rice Seed", "all purpose flour", "all purpose wheat flour", "baking flour", "baking wheat flour", "bread flour", "bread wheat flour", "cake flour", "cake wheat flour", "pastry flour", "pastry wheat flour", "self rising flour", "self rising wheat flour", "whole wheat flour", "whole grain wheat flour", "whole grain flour", "whole grain all purpose flour", "whole grain bread flour", "whole grain cake flour", "whole grain pastry flour", "whole grain self rising flour", "bulk flour", "cereals", "cereal grains", "cereal grain flour", "cereal grain seed", "cereal grain powder", "cereal grain starch", "cereal grain meal", "cereal grain mix", "cereal grain blend", "cereal grain product", "cereal grain food", "cereal grain feed", "cereal grain supplement", "cereal grain extract", "cereal grain concentrate", "cereal grain isolate", "cereal grain protein", "cereal grain fiber", "cereal grain starches"
            ],

            "Desserts & Sweets": [
                "Cake", "Chocolate cake", "vanilla cake", "lemon bar", "cinnamon halva bar", "dessert", "diabetic freindly dessert", "brown meal", "strawberry pie", "butterscotch pudding", "pudding", "chocolate pudding", "ice cream", "kulfi", "sorbet", "gelato", "frozen yogurt", "sundae", "milkshake", "smoothie", "fruit salad", "fruit cocktail", "fruit parfait", "fruit tart", "fruit cake", "fruit pudding", "fruit mousse", "fruit sorbet", "fruit gelato", "fruit ice cream", "fruit popsicle", "fruit freeze", "fruit slushie", "fruit smoothie bowl", "gulab jamun", "rasgulla", "jalebi", "barfi", "kheer", "halwa", "peda", "modak", "ladoo", "besan ladoo", "coconut ladoo", "motichoor ladoo", "boondi ladoo", "rava ladoo", "dry fruit ladoo", "chocolate ladoo", "milk cake", "sooji halwa", "gajar halwa", "moong dal halwa", "besan halwa", "suji halwa", "carrot halwa", "pistachio halwa", "almond halwa", "cashew halwa", "walnut halwa", "peanut halwa", "coconut halwa", "chocolate halwa", "fruit halwa", "saffron halwa", "rose halwa", "cardamom halwa", "ginger halwa", "turmeric halwa", "cinnamon halwa", "clove halwa", "nutmeg halwa", "mace halwa", "pepper halwa", "vanilla halwa", "orange halwa", "lemon halwa", "lime halwa", "bucklava", "knafeh", "maamoul", "qatayef", "basbousa", "konafa", "muhallebi", "kadayif", "saffron ice cream", "rose ice cream", "pistachio ice cream", "almond ice cream", "cashew ice cream", "walnut ice cream", "peanut ice cream", "coconut ice cream", "chocolate ice cream", "fruit ice cream", "saffron kulfi", "rose kulfi", "pistachio kulfi", "almond kulfi", "cashew kulfi", "walnut kulfi", "peanut kulfi", "coconut kulfi", "chocolate kulfi", "fruit kulfi", "saffron sorbet", "rose sorbet", "pistachio sorbet", "almond sorbet", "cashew sorbet", "walnut sorbet", "peanut sorbet", "coconut sorbet", "chocolate sorbet", "fruit sorbet", "saffron gelato", "rose gelato", "pistachio gelato", "almond gelato", "cashew gelato", "walnut gelato", "peanut gelato", "coconut gelato", "chocolate gelato", "fruit gelato", "chilled dessert"
            ],

            "Dips, Spreads & Sauces": [
                "peanut butter", "tahini", "hummus", "humus", "humus sample", "hummus sample", "chutney", "chutney powder", "peanut sauce", "mayonnaise", "guacamole", "salsa", "tzatziki", "pesto", "sriracha", "hot sauce", "barbecue sauce", "ketchup", "mustard", "soy sauce", "teriyaki sauce", "fish sauce", "vinegar", "balsamic vinegar", "apple cider vinegar", "red wine vinegar", "white wine vinegar", "rice vinegar", "distilled vinegar", "malt vinegar", "coconut aminos", "liquid aminos", "tamari sauce", "hoisin sauce", "sweet and sour sauce", "plum sauce", "duck sauce", "chili garlic sauce", "chili oil", "chili paste", "chili powder", "chili flakes", "chili seasoning", "chili rub", "chili marinade", "chili dressing", "chili dip", "szechuan sauce", "black bean sauce", "oyster sauce", "curry paste", "red curry paste", "green curry paste", "yellow curry paste", "massaman curry paste", "panang curry paste", "vindaloo curry paste", "korma curry paste", "chutney", "chutney powder", "chutney sauce", "chutney dip", "chutney dressing", "chutney marinade", "chutney spread", "chutney paste", "chutney mix", "chutney blend", "chutney rub", "chutney seasoning", "chutney spice mix", "chutney spice blend", "chutney spice rub", "chutney spice seasoning", "chutney spice paste", "chutney spice marinade", "coriander chutney", "mint chutney", "mango chutney", "tomato chutney", "onion chutney", "coconut chutney", "peanut chutney", "ginger chutney", "garlic chutney", "green chili chutney", "red chili chutney", "sweet chutney", "spicy chutney", "sour chutney", "tamarind chutney", "date chutney", "apple cider", "cider", "cider vinegar", "apple cider vinegar", "apple cider sauce", "apple cider dressing", "apple cider marinade", "apple cider dip", "apple cider glaze", "apple cider syrup", "apple cider jelly", "apple cider jam", "apple cider chutney", "apple cider relish", "apple cider salsa", "apple cider sauce mix", "apple cider spice mix", "apple cider spice blend", "apple cider spice rub", "apple cider spice seasoning", "apple cider spice paste", "apple cider spice marinade", "asiago cheese sauce", "asiago cheese dip", "asiago cheese dressing", "asiago cheese marinade", "asiago cheese spread", "asiago cheese paste", "asiago cheese mix", "asiago cheese blend", "asiago cheese rub", "asiago cheese seasoning", "asiago cheese spice mix", "asiago cheese spice blend", "asiago cheese spice rub", "asiago cheese spice seasoning", "asiago cheese spice paste", "asiago cheese spice marinade"
            ],

            "Powdered-products/supplements": [
                "whey protein concentrate", "pea powder", "pea protein", "chia flax broccoli seed smoothie mixture", "melon powder", "coconut milk powder", "corn silk powder", "commercial dietary supplements", "dried egg powder", "ceremonial kava powder"
            ],

            "Oil & Fats": [
                "oil", "rancid oil", "butter", "ghee", "olive oil", "coconut oil", "vegetable oil", "canola oil", "sunflower oil", "soybean oil", "peanut oil", "sesame oil", "corn oil", "grapeseed oil", "avocado oil", "palm oil", "safflower oil", "flaxseed oil", "hempseed oil", "walnut oil", "almond oil", "hazelnut oil", "macadamia nut oil", "pumpkin seed oil", "chia seed oil", "coconut cream", "coconut milk",
            ],

            "Bakery & Bread": [
                "Bread", "breading mixture", "breading mix", "organic buckwheat bread", "pancakes", "muffin", "croissant", "bagel", "rolls", "brioche", "ciabatta", "focaccia", "sourdough", "flatbread", "naan", "pita", "tortilla", "wraps", "lavash", "chapati", "paratha", "roti", "pancake mix", "waffle mix", "crepe mix", "pizza dough", "bread dough", "pastry dough", "cookie dough", "cake batter", "pastry", "tarts", "pies", "cookies", "brownies", "blondies", "scones", "biscuits", "crackers", "pretzels", "breadsticks", "baguette", "cake fondent" ,"brown bread", "multi-grain bread", "whole wheat bread", "rye bread", "sourdough bread", "gluten-free bread", "white bread", "sweet bread", "savory bread", "herb bread", "cheese bread", "garlic bread", "onion bread", "olive bread", "sun-dried tomato bread", "pumpkin bread", "zucchini bread", "banana bread", "apple cinnamon bread", "blueberry bread", "cranberry orange bread", "chocolate chip banana bread", "cake batter", "eggless cake batter", "eggless cake mix", "eggless pancake mix", "eggless waffle mix", "eggless crepe mix", "eggless pizza dough", "eggless bread dough", "eggless pastry dough", "eggless cookie dough", "eggless brownie batter", "eggless blondie batter", "eggless muffin batter", "eggless scone batter", "eggless biscuit batter", "eggless cracker dough", "eggless pretzel dough", "eggless breadstick dough",
            ],

            "Egg-related Foods" : [
                "omlette", "omelet", "scrambled egg", "boiled egg", "poached egg", "fried egg", "deviled egg", "egg salad", "egg sandwich", "egg wrap", "egg roll", "egg curry", 
                 "egg wash", 
            ],

            "Fruit-derived Products": [
                "coconut meat", "raw gallnut honey", "gallnut honey", "heterotrigona itama honey", "Honey", "ripe olives", "table olives", "brine from table olives", "olive brine"
            ],

            "Miscellaneous": [
                "produce", "raw food", "Food additives", "food and non contact surface", "kimchi system", "complementary food of infant", "zha-chili", "Meal", "Dinner", "Breakfast", "brown meal", "spoiled food", "poisoned food", "cheese ravioli" 
            ],

            "Chocolate-related": ["cacao-beans", "cacao nibs", "cacao powder", "cacao butter", "cacao mass", "cacao liquor", "chocolate liquor", "chocolate powder", "chocolate syrup", "chocolate sauce", "chocolate spread", "chocolate chips", "chocolate chunks", "chocolate bars", "chocolate truffles", "chocolate bonbons", "chocolate-covered nuts", "chocolate-covered fruits", "chocolate-covered pretzels"], 

            "Cooked Grains": ["Boiled rice", "cooked rice", "Cooked Grains"],

            "Canned Foods": ["canned food", "canned vegetables", "canned fruits", "canned beans", "canned soup", "canned tomatoes", "canned corn", "canned peas", "canned carrots", "canned mushrooms", "canned olives", "canned pickles", "canned fruit cocktail", "canned chili", "canned pasta sauce", "canned onions", "canned garlic", "canned artichokes", "canned asparagus", "canned beets", "canned pumpkin", "canned sweet potatoes", "canned green beans", "canned mixed vegetables", "canned fruit salad", "canned fruit cocktail syrup", "canned fruit juice", "canned fruit puree", "canned fruit compote", "canned fruit preserves", "canned fruit spread"],
            
            "Frozen Foods": ["frozen food", "frozen vegetables", "frozen fruits", "frozen meals", "frozen dinners", "frozen pizza", "frozen desserts", "frozen snacks", "frozen appetizers", "frozen breakfast foods", "frozen entrees", "frozen side dishes", ],

            "Sprouts" : ["sprouts", "sprouted grains", "sprouted seeds", "sprouted legumes", "sprouted beans", "sprouted lentils", "sprouted chickpeas", "sprouted peas", "sprouted nuts", "sprouted quinoa", "sprouted buckwheat", "sprouted wheat", "sprouted barley", "sprouted millet", "sprouted amaranth", "sprouted rice", "sprouted oats", "sprouted corn", "sprouted soybeans", "alfalfa sprout", "Sprout", "moong sprout", "moong bean sprout", "pulse sprout", "beetroot sprout", "broccoli sprout", "radish sprout", "sunflower sprout", "pea sprout", "lentil sprout", "chickpea sprout", "mung bean sprout", "soybean sprout", "alfalfa sprout", "clover sprout", "mustard sprout", "cabbage sprout", "cauliflower sprout", "kale sprout", "brussels sprouts", "broccoli rabe sprout", "broccolini sprout", "rapini sprout", "arugula sprout", "watercress sprout", "radish seedling",],

            "Cooking Place" : ["kitchen table", "dining table", "kitchen knife", "cookware", "kitchen ware", "kitchen countertop", "kitchen sink taps", "kitchen unused sink", "kitchen unused sink drain", "kitchen sink", "kitchen cutting board", "kitchen utensils", "kitchen appliances", "kitchen equipment", "kitchen tools", "kitchen surfaces", "kitchen area", "kitchen environment", "kitchen space", "kitchen setting", "kitchen worktop", "kitchen work surface", "kitchen work area", "kitchen work environment", "kitchen work setting", "kitchen work space", "bakery Environment", "degreasser", "kitchen floor", "kitchen wall", "kitchen ceiling", "kitchen cabinet", "kitchen drawer", "kitchen shelf", "kitchen countertop surface", "kitchen sink surface", "kitchen cutting board surface", "kitchen utensil surface", "kitchen appliance surface", "kitchen equipment surface", "kitchen tool surface", "kitchen area surface", "kitchen environment surface", "kitchen space surface", "kitchen setting surface", "kitchen worktop surface", "kitchen work area surface", "kitchen work environment surface", "kitchen work setting surface", "kitchen work space surface", "degresser sample"], 
            "Restaurants" : ["restaurant", "food court", "fast food restaurant", "dining hall", "cafeteria", "buffet restaurant", "takeout restaurant", "food truck", "food stall", "food stand", "food kiosk", "food bar", "food market", "food festival", "food fair", "coldroom for food storage"],

            "Animal Feed Processing Plant" : ["animal feed processing plant", "animal feed production facility", "animal feed manufacturing plant", "animal feed processing site", "animal feed production site", "animal feed manufacturing facility", "animal feed processing facility", "animal feed production plant", "animal feed manufacturing site"],

            "Food Processing Plant" : ["Food Processing Plant"],

            "Imitation Meat" : ["imitation meat", "plant-based meat", "vegan meat", "vegetarian meat", "meat substitute", "meat alternative", "meat analog", "meatless meat", "fake meat", "mock meat", "simulated meat", "vegetable protein meat", "soy-based meat", "pea-based meat", "mushroom-based meat", "seitan-based meat", "jackfruit-based meat", "tempeh-based meat", "tofu-based meat", "mycoprotein-based meat"],

            

        },

        "Water": {
            "Others": ['spring',"aquarium","aquarium water", "Water", "Freshwater", "Surface Water", "Stormwater", "Drinking Water", "Mineral Water", "Desalinated Water", "Filtered Water", "Purified Water", "Processed Water", "Table Water", "Potable Water", "agricultureal water", "freshwater ditch", "ditch water", "ditch", "flood water", "water for irrigation", "water pipe", "water spray", "water statue"],
            "Tap Water": ["tap water", "municipal water", "city water", "household water", "domestic water", "urban water", "public water supply", "drinking tap water", "tap water sample", "tap water source", "tap water supply", "tap water quality", "tap water treatment", "residence cold water tap"],
            "Seawater": ["seawater", "Sea Water", "Saltwater", "Ocean Water", "Marine Water", "Coastal Water", "Saltwater Lagoon", "Brine Water", "Deep Sea Water", "Tidal Water", "Estuarine Water", "Marine Brine", "Ocean Brine", "Open Sea Water", "surface sea water", "costal surface sea water", "marine sample from the pacific ocean", "marine", "marinewater"],
            "Groundwater": ['mountain spring', "Groundwater", "Subsurface Water", "Aquifer Water", "Well Water", "Spring Water", "Artesian Water", "Borehole Water", "Ground Water Table", "Subterranean Water", "Underground Water", "Confined Aquifer Water", "Unconfined Aquifer Water", "Natural Spring Water"],
            "Wastewater": ["Wastewater", "Sewage", "Municipal Wastewater", "Graywater", "Blackwater", "Septic Tank Water", "Contaminated Water", "Recycled Wastewater", "Wastewater Treatment", "Treated Effluent", "WWTP", "Chlorinated WWTP", "Non Chlorinated WWTP", "Industrial Wastewater", "Sanitary Effluent", "Effluent Water", "Bio-Treated Effluent", "Processed Wastewater", "Sewage Effluent", "Bioreactor Effluent", "Chemical Wastewater", "Dechlorinated Effluent", "Sewagewater", "Water from Sewage", "Raw Sewage", "sewage influent", "domestic sewage"],
            "Riverwater": ["River Water", "Stream Water", "Creek Water", "Brook Water", "Riverside Water", "Tributary Water", "River Runoff", "River", "Stream", "Tributary", "Floodplain Water", "River Basin Water", "Freshwater River", "Streambed Water", "River Current Water", "Delta Water", "Mountain Stream Water", "freshwater stream sediment", "waterfall", "river sample", "water stream"],
            "Lakewater": ["Lake Water", "Lake Sample", "Freshwater Lake", "Saltwater Lake", "Lake Runoff", "Lake Basin Water", "Freshwater Lake Sample", "Reservoir Water", "Dam Water", "Inland Lake Water", "Natural Lake Water", "Artificial Lake Water", "Lake Surface Water", "Calm Lake Water", "surface of lake", "lake", "artificial lake water",],
            "Rainwater": ["Rainwater", "Precipitation", "Catchment Rainwater", "Weathered Rainwater", "Rainfall", "Downpour Water", "Condensed Rainwater", "Drizzle Water", "Storm Rainwater", "Rainstorm Water", "Rain Droplets", "Collected Rainwater", "standing rain water", "rain sample"],
            "Pondwater": ["Pond Water", "Shallow Water Pond", "Swamp Pond Water", "Still Water Pond", "Pond Effluent", "Manmade Pond Water", "Natural Pond Water", "Garden Pond Water", "Freshwater Pond", "Stagnant Pond Water", "Eutrophic Pond Water", "Aquatic Pond Water", "Fish Pond Water"],
            "Industrial Effluent": ["Industrial Effluent", "Chemical Effluent", "Polluted Industrial Water", "Industrial Discharge", "Chemical Plant Effluent", "environmental effluent", "Effluent", "Factory Wastewater", "Manufacturing Effluent", "Heavy Metal Effluent", "Industrial Byproduct Water", "Thermal Power Plant Effluent", "Refinery Effluent", "Oil and Gas Effluent", "Industrial Sludge Water"],
            "Metagenomic": ["Marine Metagenome", "Freshwater metagenome"],
            "Agricultural Water": ["water from the field", "water from farm", ],
            "Biofilm" : ["Biofilm", "Biofilm Sample", "Biofilm Isolate", "Biofilm Specimen", "Biofilm Culture", "Biofilm Swab", "Biofilm Slime", "Biofilm Matrix", "Biofilm Layer", "Biofilm Community", "Marine Biofilm"],
            "Polluted Water": ["Polluted Water", "Contaminated Water", "Toxic Water", "Hazardous Water", "Chemical Contaminated Water", "Heavy Metal Contaminated Water", "Industrial Contaminated Water", "Water Pollution", "Water Contamination", "Environmental Contaminated Water", "Water Toxicity", "Water Hazardous Waste", "burned water", "oil contaminated water"],
            
        },

        "Soil" : {
            "Others" : ["Farm water", "soil microorganisms","Farm flooded water", "Soil sample","Soil", "Generic Soil", "Earth", "Ground", "Soil Particles", "Dirt", "Sand", "Loam", "Rocky Soil", "Non-Specific Soil", "Soil Debris", "Topsoil", "Subsoil", "Soil Matter", "Earthy Matter", "Natural Soil", "Soil Mixture", "background grass", "not manured soil", "cave soil", "black soil", "maybe soil", "lake soil", "red latosol", 'Grass Without Manure', 'Grass No Manure', 'alkaline soil', "archived soil", "arid soil", "bare soil", "dry soil", "hypersaline soil", "island soil", "non-saline soil", "red soil", "saline soil", "saline-alkali soil", "sandy soil", "wood soil", "surface soil", "subsurface soil", "soil sample", "soil core", "soil profile", "soil horizon", "soil layer", "soil column", "soil block", "soil aggregate", "soil clod", "soil crumb", "soil particle", "soil grain", "soil fragment", "soil chunk", "soil piece", "soil inside the cave", "soft rock", "soil surface", "soil subsurface", "soils"],
            "Moraine Soil" : ["Moraine Soil", "Glacial Till Soil", "Glacial Drift Soil", "Glacially Deposited Soil", "Till Soil", "Glacially Affected Soil", "Glacially Formed Soil", "Glacially Derived Soil", "Glacially Transported Soil", "Glacially Modified Soil", "Glacially Influenced Soil"],
            "Alpine Soil" : ["Alpine Soil", "Mountain Soil", "High-Altitude Soil", "Subalpine Soil", "Alpine Tundra Soil", "Alpine Meadow Soil", "Alpine Grassland Soil", "Alpine Ecosystem Soil", "Alpine Region Soil", "Alpine Habitat Soil", "Alpine Biome Soil"],
            "Paddy Soil": ["Paddy Soil", "Rice Paddy Soil", "Wetland Rice Soil", "Flooded Rice Soil", "Irrigated Rice Soil", "Paddy Field Soil", "Rice Cultivation Soil", "Wet Rice Soil", "Paddy Field Sediment", "Rice Growing Soil", "Paddy Land Soil", "Rice Cultivation Area Soil", "Paddy Waterlogged Soil", "rice field"],
            "Contaminated Soil": ["Contaminated Soil", "Polluted Soil", "Toxic Soil", "Hazardous Soil", "Chemical Contaminated Soil", "Heavy Metal Contaminated Soil", "Industrial Contaminated Soil", "Soil Pollution", "Soil Contamination", "Environmental Contaminated Soil", "Soil Toxicity", "Soil Hazardous Waste", "burned soil", "oil contaminated soil"],
            "Agricultural Soil": ["pesticide contaminated farm soil", "Sugarcane soil", "sugar cane soil", "Pesticide contaminated soil","Field","Agricultural Soil", "Farmland Soil", "Farming Soil", "Cultivated Soil", "Crop Soil", "Planting Soil", "Agricultural Land Soil", "Crop-Field Soil", "Farm Soil", "Tillable Soil", "Field Soil", "Farmland", "Soil For Agriculture", "Irrigated Soil", "Cultivated Land Soil", "commercial crop production field", "chestnut soil", "hazelnut soil", "almond soil", "walnut soil", "pecan soil", "macadamia soil", "pistachio soil", "cacao soil", "coffee soil", "tea soil", "fruit tree soil", "vegetable garden soil", "herb garden soil", "flower garden soil", "ornamental plant soil", "landscaping soil", "lettuce soil", "beet soil", "Citrus soil", "garden soil", "orchard", "mulberry soil", "sesame field", "soil sample from farm", "soil, beetroot field"],
            "Forest Soil": ["Forest","Forest Soil", "Woodland Soil", "Rainforest Soil", "Temperate Forest Soil", "Tropical Forest Soil", "Deciduous Forest Soil", "Coniferous Forest Soil", "Pine Forest Soil", "Wooded Area Soil", "Mountain Forest Soil", "Tropical Rainforest Soil", "Boreal Forest Soil", "Forest Floor Soil", "Old Growth Forest Soil", "Forest Soil Type", "Natural Forest Soil", "grassland soil", "jungle soil", "rain forest", "alpine forest soil"],
            "Desert Soil": ["Desert Soil", "Arid Soil", "Sandy Soil", "Dry Soil", "Dusty Soil", "desert rock", "Desert Sand", "Desert Terrain Soil", "Saline Soil", "Semi-Arid Soil", "Desertified Soil", "Cactus Soil", "Shrubland Soil", "Bare Soil", "Dune Soil", "Desert Ecosystem Soil", "Low Organic Soil", "cold desert"],
            "Urban Soil/Fields": ["Urban Soil", "City Soil", "Town Soil", "Urban Environment Soil", "Paved Area Soil", "Built-Up Soil", "Urban Landscape Soil", "Street Soil", "Roadside Soil", "Park Soil", "Urban Garden Soil", "City Agricultural Soil", "Residential Area Soil", "Commercial Area Soil", "Soil In Urban Areas", "Soil In Suburbs", 'Ground Grass', "Grass Ground", "urban sediment", "urban soil", "urban soil sample", "urban soil core", "urban soil profile", "urban soil horizon", "urban soil layer", "urban soil column", "urban soil block", "urban soil aggregate", "urban soil clod", "urban soil crumb", "urban soil particle", "urban soil grain", "urban soil fragment", "urban soil chunk", "urban soil piece"],
            "Manured Lands": ['Manure Field', 'Field Manure', 'Manure On Field', 'Grass Containing Manure', 'Grass And Manure', 'Grass With Manure',  'manure soil'],
            "Industrial Soil": ["Industrial Soil", "Factory Soil", "Polluted Soil", "Toxic Soil", "Mining Soil", "Industrial Waste Soil", "Manufacturing Soil", "Industrial Plant Soil", "Chemical Soil", "Soil In Industrial Zones", "Pollution Affected Soil", "Urban Industrial Soil", "Industrial Waste Disposal Area Soil", "Mine Tailings", "Excavated Soil", "Processed Soil", "Mineral Soil", "Heavy Metal Soil", "Mining Waste Soil", "Mine Waste Soil", "Extracted Soil", "Mining Dump Soil", "Waste Rock Soil", "Tailings Pile", "Ore Residue Soil", "Residual Mining Soil", "Mine soil", "Mining Site Soil", "Mining Area Soil", "Mining Region Soil", "Mining District Soil", "Mining Zone Soil", "Mining Field Soil", "Mining Land Soil", "Mining Territory Soil", "Mining Ground Soil", "Mining Landscape Soil", "Mining Environment Soil", "coal mine soil", "mangenese mine soil"],
            "Volcanic Soil": ["Volcanic Soil", "Lava Soil", "Basaltic Soil", "Tuff Soil", "Ash Soil", "Pumice Soil", "Volcanic Ash Soil", "Volcanic Eruption Soil", "Volcanic Rock Soil", "Fertile Volcanic Soil", "Tephra Soil", "Hawaiian Volcanic Soil", "Andesite Soil", "Volcanic Island Soil"],
            "Compost": ["Compost", "Organic Compost", "Soil Enrichment", "Soil Amendment", "Decomposed Organic Matter", "Humus", "Compost Soil", "Mulch", "Organic Matter", "Plant-Based Compost", "Garden Compost", "Organic Fertilizer", "Green Waste Compost", "Kitchen Waste Compost", "Composted Material", "Composted Soil", "Soil Conditioner", "composted wastes", "composting sample", "raw compost", "vermicompost", "vermicompost", "vermi-compost", "vermi compost", "vermicomposted soil", "vermicomposted material", "vermicomposted organic matter", "vermicomposted plant material", "vermicomposted kitchen waste", "vermicomposted garden waste", "vermicomposted green waste", "vermicomposted food waste"],
            "Clay Soil": ["Clay Soil", "Heavy Soil", "Sticky Soil", "Loamy Clay Soil", "Clay-Rich Soil", "Clay Loam Soil", "Red Clay Soil", "Yellow Clay Soil", "Black Clay Soil", "Clay Texture Soil", "Plastic Soil", "Clayey Soil", "Dense Clay Soil", "Moist Clay Soil", "Clay-Based Soil", "Stiff Soil"],
            "Permafrost": ["Permafrost", "Frozen Soil", "Cold Soil", "Tundra Soil", "Frozen Earth", "Glacial Soil", "Arctic Soil", "Polar Soil", "Deep Frozen Soil", "Subsurface Ice Soil", "Permanent Ice Soil", "Arctic Tundra Soil", "Permafrost Layers", "Ice-Rich Soil", "Permafrost Soil Formation"],
            "Peat Soil": ["Peat Soil", "Moorland Soil", "Sphagnum Soil", "Bog Soil", "Swamp Soil", "Peat Moss Soil", "Organic Soil", "Hydric Peat Soil", "Peaty Soil", "High Organic Matter Soil", "Wetland Peat Soil", "Tropical Peat Soil", "Humic Soil", "Turf Soil", "Acidic Peat Soil"],
            "Sedimentary Soil": ["Sedimentary Soil", "Sandstone Soil", "Shale Soil", "Claystone Soil", "Fossiliferous Soil", "Gravelly Soil", "Riverbed Soil", "Alluvial Soil", "Soil From Sedimentary Layers", "Depositional Soil", "Soil From River Deposits", "Sediment-Based Soil", "Deltaic Soil", "Sedimentary Basin Soil", "Sediment", "Sediment Sample", "River Sediment", "Lake Sediment", "Marine Sediment", "Aquatic Sediment", "Water Sediment", "Bottom Sediment", "Sediment Layer", "Sediment Core", "Sediment Deposit", "Pond Sediment", "marine sediment", "freshwater sediment", "freshwater stream sediment", "sediment associated", "sediment associated with basin", "sediment", "river mud", "mud from river", "surface mud", "sediment samples", "lagoon sediment", "sediment from 275 cmbsf of arabian sea"],
            "Coastal Soil": ["Coastal Soil", "Beach Soil", "Sandy Beach Soil", "Saltwater Soil", "Marine Soil", "Estuarine Soil", "Shoreline Soil", "Tidal Soil", "Salt Marsh Soil", "Mangrove Soil", "Coastal Wetland Soil", "Seashore Soil", "Rocky Coast Soil", "Seaside Soil", "Coastal Ecosystem Soil", "Salty Soil", "Beach Sand", "Coastal Sand", "Seaside Sand", "Marine Sand", "Tidal Sand", "Sandy Beach", "Shore Sand", "Sea Sand", "Coastal Sediment", "White Sand", "Golden Sand", "Coral Sand", "Beach Area Sand", "atlantic intertidal shore"],
            "Radioactive Soil": ["uranium ore deposits soil", "uranium mine soil", "uranium mining site soil", "uranium-contaminated soil", "radioactive soil", "radioactively contaminated soil", "radioactive waste site soil", "radioactive contamination site soil", "radioactive material site soil", "radioactive substance site soil", "radioactive element site soil", "radioactive isotope site soil", "radioactive nuclide site soil", "radioactive particle site soil", "radioactive dust site soil", "radioactive debris site soil"],
            "Salt-based Soil": ["Salt Based Soil", "Saline Soil", "Sodic Soil", "Salinized Soil", "Salt-Affected Soil", "Salinity-Impacted Soil", "Soil Salinity", "High Salinity Soil", "Soil With High Salt Content", "Salty Soil", "Soil With Salt Accumulation", "Soil With High Sodium Content", "Beach"],
            "Landfill Soil": ["Landfill Soil", "Waste Soil", "Garbage Soil", "Landfill Waste", "Polluted Landfill Soil", "Landfill Waste Soil", "Municipal Waste Soil", "Waste Disposal Soil", "Trash Heap Soil", "Dump Soil", "Hazardous Waste Soil", "Recycled Landfill Soil", "Garbage Pit Soil", "Landfill Contaminated Soil", "Soil In Landfill Sites", "Decomposed Waste Soil", "landfill site for plastic",],
            "Rhizosphere Soil": ["Rhizosphere harvested","Rhizosphere Soil", "Root Nodule Soil", "Rhizobium Soil", "Leguminous Soil", "Nodular Soil", "Nitrogen Fixing Soil", "Rhizobial Soil", "Soil With Nitrogen-Fixing Bacteria", "Rhizobium Inoculated Soil", "Root Symbiosis Soil", "Root Zone Soil", "Legume-Rhizobium Soil", "Rhizosphere",],
            "Soil Crust" : ["Soil Crust", "Biological Soil Crust", "Cryptobiotic Soil Crust", "Soil Surface Crust", "Soil Microbial Crust", "Soil Biological Layer", "Soil Surface Layer", "Soil Microbial Layer", "Soil Surface Biofilm", "Soil Surface Microbiome", "Soil Surface Community"],
            "Soil Microbiome": ["Soil Microbiome", "Soil Microbial Community", "Soil Microbial Diversity", "Soil Microbial Ecology", "Soil Microbial Analysis", "Soil Microbial Structure", "Soil Microbial Function", "Soil Microbial Interactions", "Soil Microbial Dynamics", "Soil Microbial Composition", "Forest soil community", "forest biome", "forest microbiome"],
            "Metagenomic": ["Soil Metagenome", "Soil Microbiome", "Soil Microbial Community", "Soil Microbial Diversity", "Soil Microbial Ecology", "Soil Microbial Analysis", "Soil Microbial Structure", "Soil Microbial Function", "Soil Microbial Interactions", "Soil Microbial Dynamics", "Soil Microbial Composition", "Forest soil community", "forest biome", "forest microbiome"],
        },

        "Natural Habitats" : {
            "Marsh" : ["Marsh", "Wetland", "Swamp", "Bogs", "Wetland Ecosystem", "Wetland Habitat", "Wetland Biome", "Wetland Area", "Wetland Environment", "Wetland Flora and Fauna", "Wetland Biodiversity"],
            "Basin" : ["Basin", "Basin Ecosystem", "Basin Habitat", "Basin Biome", "Basin Area", "Basin Environment", "Basin Flora and Fauna", "Basin Biodiversity", "river basin", "lake basin", "orinoco river basin", "amazon river basin", ],
            "Salt Pan": ["Salt Pan", "Salt Flat", "Salt Lake", "Salt Desert", "natural saltern", "Salt marsh"],
            "Salt Mine" : ["Salt Mine", "Salt Mine Environment", "Salt Mine Habitat", "Salt Mine Biome", "Salt Mine Area", "Salt Mine Flora and Fauna", "Salt Mine Biodiversity"],
            "Snowy": ["Snowy", "Snowy Environment", "Snowy Habitat", "Snowy Biome", "Snowy Area", "Snowy Flora and Fauna", "Snowy Biodiversity", "snowy region", "snowy landscape", "snowy mountain", "snowy forest", "snowy tundra", "snowy biome", "snowy ecosystem", "snowy habitat", "snowy area", "snowy environment", "surface snow"],
            "Glacial" : ["Glacial ice", "glacial ice water", "glacial water", "glacial surface ice"],
            "Corals": ["Corals", "pocillopora grandis", "Antler corals", "surface of a marine sponge", "Coral", "Coral Reef", "Coral Sample", "Coral Polyp", "Coral Skeleton", "Coral Fragment", "Coral Tissue", "Coral Algae", "Coral Microbiome", "Coral Larvae", "Coral Bleaching", "Coral Ecosystem", "Coral Community", "Coral Species", "Coral Genotype", "Coral Morphology", "Coral Diversity", "Coral Health", "marine sponge", "sponge", "demosponge", "hexactinellid", "glass sponge", "marine sponge", "acropora", "coral reef", "coral sample", "coral polyp", "coral skeleton", "coral fragment", "coral tissue", "coral algae", "coral microbiome", "coral larvae", "coral bleaching", "coral ecosystem", "coral community", "coral species", "coral genotype", "coral morphology", "coral diversity", "coral health"],
        },

        "Space" : {
            "International Space Station": ["International Space Station", "ISS Air", "Space Station Air", "Microgravity Air", "Orbital Air", "Spacecraft Air", "Extraterrestrial Air", "Zero Gravity Air", "Space Environment Air", "Space Habitat Air", "Space", "extraterrestrial environment", "zero gravity", "outer space", "mars surface", "ISS habitat"],
            "Incubated in Space": ["Incubated in Space", "Space Incubation", "Microgravity Incubation", "Spacecraft Incubation", "Extraterrestrial Incubation", "Zero Gravity Incubation", "Space Environment Incubation", "Space Habitat Incubation", "Spacecraft Culture", "Space", "extraterrestrial culture", "zero gravity culture", "outer space culture", "mars surface culture", "ISS culture"],
            "Spacecraft": ["Spacecraft", "Space Shuttle", "Orbital Module", "Lunar Module", "Mars Rover", "Space Capsule", "Satellite", "Interplanetary Spacecraft", "Extraterrestrial Vehicle", "Space Exploration Vehicle", "Space Probe"],
            "Extraterrestrial": ["Extraterrestrial", "Alien Environment", "Extraterrestrial Atmosphere", "Extraterrestrial Surface", "Extraterrestrial Soil", "Extraterrestrial Habitat", "Extraterrestrial Ecosystem", "Extraterrestrial Life", "Extraterrestrial Microbial Life", "Extraterrestrial Organisms", "Extraterrestrial Samples"],
            "Mars": ["Mars", "Mars Surface", "Mars Soil", "Mars Atmosphere", "Mars Environment", "Mars Habitat", "Mars Ecosystem", "Mars Microbial Life", "Mars Samples", "Martian Soil", "Martian Atmosphere", "Martian Environment", "Martian Habitat", "Martian Ecosystem"],
            "Lunar": ["Lunar", "Moon Surface", "Lunar Soil", "Lunar Atmosphere", "Lunar Environment", "Lunar Habitat", "Lunar Ecosystem", "Lunar Microbial Life", "Lunar Samples", "Moon Soil", "Moon Atmosphere", "Moon Environment", "Moon Habitat", "Moon Ecosystem"],
            "Asteroid": ["Asteroid", "Asteroid Surface", "Asteroid Soil", "Asteroid Atmosphere", "Asteroid Environment", "Asteroid Habitat", "Asteroid Ecosystem", "Asteroid Microbial Life", "Asteroid Samples", "Space Rock"],
            "Comet": ["Comet", "Comet Surface", "Comet Soil", "Comet Atmosphere", "Comet Environment", "Comet Habitat", "Comet Ecosystem", "Comet Microbial Life", "Comet Samples"],
            "Planetary": ["Planetary", "Planetary Surface", "Planetary Soil", "Planetary Atmosphere", "Planetary Environment", "Planetary Habitat", "Planetary Ecosystem", "Planetary Microbial Life", "Planetary Samples"],
        },

        "Air": {
            "Others": ["Air"],
            "Indoor Air": ["Indoor Air", "Home Air", "Office Air", "Building Air", "Indoor Air Quality", "Room Air", "Living Room Air", "Bedroom Air", "Indoor Environment Air", "Indoor Air Samples", "Indoor Air Composition", "Air In Homes", "Indoor Climate", "Indoor Air Quality Assessment", "Heated Indoor Air", "barn air",],
            "Outdoor Air": ["Outdoor Air", "City Air", "Outdoor Air Quality", "Urban Air", "Country Air", "Park Air", "Outdoor Environment Air", "Open Environment Air", "Outdoor Air Samples"],
            "HVAC System": ["Hvac System", "Heating, Ventilation, And Air Conditioning", "Hvac Air", "Air System", "Ducted Air System", "Ventilation System", "Central Air", "Hvac System Air", "Building Air System", "Air Circulation System", "Hvac Duct Air", "Air Filtered System", "Commercial Hvac System", "Residential Hvac System", "Ventilation Ducts", "Ducted Air", "Air Ducts", "Airflow Ducts", "Hvac Ducts", "Ductwork", "Ducted Ventilation", "Air Duct System", "Hvac Duct System", "Ventilation Duct Air", "Ductwork Contamination", "Duct Air Samples", "Residential Duct System", "Commercial Duct System", "air conditioner"],
            "Bioaerosols": ["Aerosols", "Bioaerosols", "Biological Aerosols", "Microbial Aerosols", "Fungal Aerosols", "Bacterial Aerosols", "Viruses In Air", "Spores In Air", "Microorganisms In Air", "Bio-Contaminants In Air", "Indoor Bioaerosols", "Outdoor Bioaerosols", "Airborne Pathogens", "Airborne Microorganisms"],
            "Factory Emissions": ["Factory Emissions", "Industrial Emissions", "Industrial Pollution", "Factory Air", "Manufacturing Emissions", "Emissions From Factories", "Factory Exhaust", "Industrial Air Pollutants", "Airborne Factory Waste", "Pollutant Particles From Factories", "Chemical Emissions", "Factory Gaseous Emissions"],
            "Smoke Particles": ["Smoke Particles", "Smoke", "Smoke In Air", "Particulate Smoke", "Pollution From Smoke", "Fire Smoke", "Wildfire Smoke", "Smoke Exposure", "Airborne Smoke Particles", "Soot In Air", "Carcinogenic Smoke", "Aerosolized Smoke", "Burning Smoke"],
            "Pollen": ["Pollen", "Pollen Particles", "Flower Pollen", "Tree Pollen", "Grass Pollen", "Pollen Grains", "Airborne Pollen", "Pollination Particles", "Pollen In Air", "Pollen Exposure", "Pollen In Atmosphere", "Allergen Pollen", "Seasonal Pollen", "Pollen Clouds"],
            "Greenhouse Air": ["Greenhouse Air", "Controlled Greenhouse Air", "Greenhouse Air Quality", "Horticultural Air", "Indoor Greenhouse Air", "Greenhouse Environment", "Greenhouse Ventilation Air", "Hydroponic Air", "Agricultural Greenhouse Air", "Greenhouse Gas Concentration", "Atmosphere In Greenhouses", "High Humidity Greenhouse Air"],
            "Lab Air": ["Lab Air", "Laboratory Air", "Research Lab Air", "Laboratory Environment", "Cleanroom Air", "Air In Labs", "Chemical Lab Air", "Biological Lab Air", "Controlled Lab Air", "Indoor Laboratory Air", "Air Quality In Labs", "Sterile Lab Air", "Research Environment Air", "Pollution-Free Lab Air"],
            "Smog": ["Smog", "Pollution Smog", "Industrial Smog", "Smoggy Air", "Fossil Fuel Emissions", "Particulate Smog", "Smog From Cars", "Urban Smog", "Heavy Smog", "Fog And Smog", "Smog In Cities", "Haze And Smog", "Chemical Smog", "Photochemical Smog"],
            "Dust And Pollutants": ["Airborne Pollutants", "Pollutants In Air", "Air Pollution", "Toxic Air Particles", "Airborne Toxic Materials", "Environmental Pollutants", "Pollution Particles", "Pollutants From Industry", "Car Emissions In Air", "Airborne Chemicals", "Hazardous Pollutants", "Particulate Matter In Air", "Airborne Contaminants", "Pollution Aerosols", "Chemical Contaminants In Air", "Dust Particles", "Dust", "Dust In Air", "Dust Samples", "Particulate Matter", "Pm2.5", "Pm10", "Indoor Dust", "Outdoor Dust", "Fine Dust", "Coarse Dust", "Respirable Dust", "Dust Accumulation", "Dust Mites", "Dust Storm Air"],
            "High Altitude Air": ["High Altitude Air", "Mountain Air", "High Elevation Air", "Thin Air", "Alpine Air", "High Mountain Air", "High Altitude Atmosphere", "High Altitude Environment", "High Altitude Samples", "Air At High Elevation"],
        },

        "Plant": {
            "Phyllosphere": ["phyllosphere","phyllosphere sample","phyllosphere isolate","phyllosphere mucilage", "vine", "phyllosphere exudate","phyllosphere exudates","phyllosphere exudate sample","phyllosphere exudate isolate", "phyllosphere material", "aboveground tissues", "mulberry wine", "wine", "petiole"],

            "Nectar" : ["nectar", "nectaries", "plant nectar", "nectary"],
            
            "Root": ["Aerial root sample","aerial root mucilage","Arial Root","root isolate","plant toot","Root", "Root Sample", "Edible Root", "Taproot", "Fibrous Root", "Tuber", "Starchy Root", "Underground Stem", "Bulb", "Rhizome", "Root System", "Root Cutting", "Carrot", "Potato", "Sweet Potato", "Beetroot", "Radish", "Turnip", "Ginger", "Garlic", "Onion", "Yam", "Cassava", "Parsnip", "Daikon", "Cassava Root", "Garlic Bulb", "Radish Root", "Turmeric", "Daikon Root", "Sugar Beet", "Jerusalem Artichoke", "Yam Root", "Salsify Root", "Chicory Root", "Artichoke Root", "Burdock Root", "Dandelion Root", "Sunchoke", "Lotus Root", "Sweet Yam", "Beet", "Cucumber Root", "Mallow Root", "Arrowroot", "Black Salsify", "Chufa Tuber", "Water Chestnut", "Burdock", "Konjac Root", "Water Yam", "Kohlrabi Root", "Salsify", "Crosne", "Black Garlic", "Jerusalem Artichoke Tuber", "Water Lily Root", "mangrove roots", "mangrove root", "mangrove root sample", "mangrove root isolate", "mangrove root mucilage", "mangrove root exudate", "mangrove root exudates", "mangrove root exudate sample", "mangrove root exudate isolate", "wing root", "wingstem root", "wingstem root sample", "wingstem root isolate", "wingstem root mucilage", "wingstem root exudate", "root tissues of peanut", "peanut roots", "peanut root nodule", "roots of healthy hydroponic lamb's lettuce plants", "corn roots", "roots of a desert plant", "endophytic isolate", "endophyte", "endosphere", "plant nodule", "plant root nodules", "root surface", "root tuber", "surface sterilized root", "Potato scab lesion", "chickpea root nodules", "soybean root nodules", "legume root nodules", "legume nodule", "root nodule", "root nodule sample", "root nodule isolate", "root nodule bacteria", "root nodule fungi", "root nodule microorganisms", "root nodule microbes", "root nodule"],
            
            "Root Pathogen": ["Root Pathogen", "Root Disease", "Root Rot", "Root Infection", "Root Fungus", "Root Bacteria", "Root Nematodes", "Root Parasites", "Root Pathogenic Fungi", "Root Pathogenic Bacteria", 
            "Root Pathogenic Nematodes", "Root Pathogen Samples", "meloidogyne sp.", "meloidogyne incognita", "meloidogyne javanica", "meloidogyne hapla", "meloidogyne arenaria", "root-knot nematodes", "root-knot nematode", "root-knot nematode sample", "root knot nematode egg mass", "Phytophthora alni", "root bacterial wilt"], 
            
            "Leaf": ["plant leaf","Leaf", "Leaf Sample", "pear tree leaf", "strawberry leaves", "Green Leaf", "Plant Leaf", "Tender Leaf", "banana leaf", "Tropical Leaf", "Herbaceous Leaf", "Leafy Greens", "Chard", "Mustard Greens", "Collard Greens", "Dandelion Greens", "Moringa Leaf", "Brussels Sprouts", "Bean Leaves", "Tomato Leaf", "Alfalfa", "Watercress", "Artichoke Leaf", "Radicchio", "Romaine Lettuce", "Endive", "Escarole", "Mache", "Sorrel", "Arugula", "Chamomile", "Lemon Balm", "Thyme Leaf", "Bitter Melon Leaf", "Ginger Leaf", "Paprika Leaf", "Sage", "Oregano Leaf", "Curry Leaf", "Chili Leaf", "Apple Leaf", "Pear Leaf", "Peach Leaf", "Almond Leaf", "Lemon Leaf", "Citrus Leaf", "Green Tea Leaf", "Aloe Vera Leaf", "Spearmint", "Peppermint", "Bay Leaf", "Coriander Leaf", "Lime Leaf", "Grape Leaf", "Pineapple Leaf", "Sorghum Leaf", "Sunflower Leaf", "Tobacco Leaf", "Cassava Leaf", "Passionfruit Leaf", "Bamboo Leaf", "Palm Leaf", "Mango Leaf", "Chili Pepper Leaf", "Okra Leaf", "Cucumber Leaf", "Squash Leaf", "Pumpkin Leaf", "Zucchini Leaf", "Cabbage Leaf", "Broccoli Leaf", "Cauliflower Leaf", "Kale Leaf", "Spinach Leaf", "Chard Leaf", "Beet Greens", "the aging flue-cured tobacco leaves", "the flue-cured tobacco leaves", "oak leaf", "oak leaves", "oak leaf sample", "oak leaf isolate", "diseased rice seeds and leaf sheaths", "ventilago sp. leaf", "ventilago sp. leaf sample", "ventilago sp. leaf isolate", "leaves", "leaf lesion"],
            
            "Seed/Seedling": ["Seedling", "Tomato seedling", "Seed", "Seed Sample", "Plant Seed", "Tree Seed", "Flower Seed", "Corn Seed", "Bean Seed", "Soybean Seed", "Cotton Seed", "Pea Seed", "Sunflower Seed", "Pumpkin Seed", "Melon Seed", "Chili Seed", "Tomato Seed", "Cucumber Seed", "Apple Seed", "Grape Seed", "Pepper Seed", "Seedling", "Seed Pod", "Radish Seed", "Alfalfa Seed", "Mustard Seed", "Carrot Seed", "Sorghum Seed", "Barley Seed", "Canola Seed", "Watermelon Seed", "Papaya Seed", "Avocado Seed", "Tomato Seeds", "Bell Pepper Seed", "Squash Seed", "Cabbage Seed", "Broccoli Seed", "Aubergine Seed", "Sweet Corn Seed", "Green Bean Seed", "Beetroot Seed", "Cantaloupe Seed", "Kiwifruit Seed", "Fennel Seed", "Coriander Seed", "Dandelion Seed", "Asparagus Seed", "Okra Seed", "Wheatgrass Seed", "Taro Seed", "Carambola Seed", "Pomegranate Seed", "Mango Seed", "Passionfruit Seed", "Lemongrass Seed", "Apple Seedling", "Pear Seed", "Banana Seed", "Peach Seed", "Pineapple Seed", "Sweet Pepper Seed", "Citrus Seed", "Apricot Seed", "Plum Seed", "Olive Seed", "Cherry Seed", "Fig Seed", "Mulberry Seed", "Persimmon Seed", "Quince Seed", "Lychee Seed", "Longan Seed", "Durian Seed", "Jackfruit Seed", "Starfruit Seed", "Sapodilla Seed", "Guava Seed", "Papaya Seeds", "alfalfa seed", "endosperm", "seed coat", "rapeseed",],
            
            "Stem": ["Stem", "Stem Sample", "Plant Stem", "Tree Stem", "Woody Stem", "Herbaceous Stem", "Flower Stem", "Shrub Stem", "Cactus Stem", "Corn Stalk", "Sugarcane Stem", "Cucumber Vine", "Tomato Stem", "Potato Stem", "Pumpkin Vine", "Watermelon Vine", "Grapevine", "Cotton Stem", "Sorghum Stem", "Coffee Plant Stem", "Cotton Plant Stem", "Chili Pepper Plant Stem", "Tobacco Plant Stem", "Banana Plant Stem", "Kale Stalk", "Brussels Sprouts Stalk", "Cabbage Stalk", "Carrot Stalk", "Ginseng Root Stalk", "Yucca Stem", "Dragon Fruit Stem", "Cabbage Stem", "Okra Stem", "Chili Stem", "Sunflower Stalk", "Sorghum Cane", "Bamboo Stalk", "Zucchini Stalk", "Eggplant Stem", "Peanut Plant Stem", "Artichoke Stem", "Artichoke Heart Stalk", "Sugar Beet Stem", "Jute Stem", "Mango Tree Stem", "Coffee Bean Stem", "Pineapple Stem", "Palm Stem", "Soybean Stem", "Banana Stalk", "Avocado Tree Stem", "Pomegranate Tree Stem", "Cassava Stem", "Cassowary Stem", "Hemp Stem", "Lemon Tree Stem", "Palm Leaf Stem", "Walnut Stem", "Fig Tree Stem", "Cherry Tree Stem", "Peach Tree Stem", "Plum Tree Stem", "Olive Tree Stem", "Mulberry Tree Stem", "Persimmon Tree Stem", "Quince Tree Stem", "Lychee Tree Stem", "Longan Tree Stem", "Durian Tree Stem", "Jackfruit Tree Stem", "Starfruit Tree Stem", "Sapodilla Tree Stem", "Guava Tree Stem", "stem of peanut", "hazelnut twig", "symptomatic shoots/leaves", "mulberry stem segments", "maize stem"],
            
            "Fruit": ["Tomato", "ear", "zea mays ear","Fruit", "Fruit Sample", "Edible Fruit", "Apple", "Banana", "Grape", "Orange", "Pear", "Cherry", "Melon", "Strawberry", "Peach", "Nectarine", "Pineapple", "Watermelon", "Kiwi", "Mango", "Blueberry", "Avocado", "Papaya", "Lemon", "citrus lemon", "Lime", "Plum", "Grapefruit", "Pomegranate", "Apricot", "Fig", "Coconut", "Lychee", "Passionfruit", "Olive", "Blackberry", "Gooseberry", "Dragon Fruit", "Starfruit", "Tangerine", "Persimmon", "Mulberry", "Black Currant", "Elderberry", "Soursop", "Jackfruit", "Mandarin", "Durian", "Longan", "Rambutan", "Custard Apple", "Cantaloupe", "Kumquat", "Pawpaw", "Tamarind", "Barbados Cherry", "Cherry Plum", "Green Apple", "Saskatoon Berry", "Red Currant", "Cranberry", "Blackcurrant", "Goji Berry", "Sea Buckthorn", "Chayote", "Pluot", "Clementine", "Lingonberry", "Juneberry", "Longan Fruit", "White Mulberry", "Yellow Watermelon", "Fuzzy Melon", "Quince", "Marula", "Bitter Orange", "Mango Chutney Fruit", "Galia muskmelon", "muskmelon", "Galia melon", "melon", "waxberry", "Chinese bayberry", "yangmei", "yumberry", "Japanese bayberry", "bayberry", "Chinese waxberry", "Japanese waxberry", "Japanese yumberry", "Japanese yangmei", "Japanese yumberry fruit", "Japanese waxberry fruit", "yangmei fruit", "yumberry fruit", "Morus", "Mulberry", "Mulberry fruit", "morus alba", "cherry tomato"],
            
            "Flowering Plant/Flower": ["gaillardia","Blanket flowers","Flowering Plant", "strawberry pollen", "Flowering Species", "Blooming Plant", "Flower-Bearing Plant", "Annual Flowering Plant", "Perennial Flowering Plant", "Sunflower", "Tulip", "Rose", "Orchid", "Lily", "Daisy", "Carnation", "Wildflower", "Lilac", "Hibiscus", "Geranium", "Chrysanthemum", "Lavender", "Violet", "Begonia", "Calendula", "Marigold", "Zinnia", "Morning Glory", "Bougainvillea", "Jasmine", "Petunia", "Poppy", "Plumeria", "Freesia", "Pansy", "Lotus Flower", "Camellia", "Daffodil", "Azalea", "Orchidaceae", "Magnolia", "Lotus", "Water Lily", "Fuchsia", "Heliotrope", "Crocus", "Orchid", "Peony", "Dandelion", "Iris", "Angelonia", "Snapdragon", "Bluebell", "Cherry Blossom", "Hibiscus Flower", "Saffron Crocus", "Snapdragon Flower", "Petunia Flower", "Morning Glory Vine", "Magnolia Flower", "Fuchsia Plant", "Camellia Plant", "Honeysuckle", "Red Hot Poker", "Hydrangea", "Alstroemeria", "Gladiolus", "Periwinkle", "Maranta", "Helichrysum", "Mock Orange", "Indigofera", "Willowherbs", "Epilobium", "floweing plant", "suaeda maritima", "Amaranthaceae", "Medicago lupulina"],
            
            "Algae": ["Rock weed", "rockweed", "seagrapes", "fucus vesiculosus", "saccharina japonica", "Kombu", "ulva", "Sea lettuce", "Cyanobacterial mats","Algae", "Green Algae", "Blue-Green Algae", "Red Algae", "Brown Algae", "Marine Algae", "Algal Bloom", "Algal Sample", "Pond Algae", "Freshwater Algae", "Seaweed", "Kelp", "Sargassum", "Spirulina", "Chlorella", "Diatom", "Phyto Plankton", "Euglena", "Porphyra", "Gracilaria", "Nori", "Agar", "Fucus", "Ulva", "Carrageenan", "Zostera", "Ascophyllum", "Fucus Vesiculosus", "Sargassum Muticum", "Ulva Lactuca", "Spirogyra", "Chara", "Volvox", "Red Tide Algae", "Tetraselmis", "Euglena Gracilis", "Coccolithophores", "Chlorophyta", "Sargassum Fusiforme", "Sea Lettuce", "Kelp Forest", "Edible Seaweed", "Algal Mat", "Sea Moss", "Rockweed", "Bladderwrack", "Nori Algae", "Wakame", "Kombu", "Dulse", "Irish Moss", "Agar-Agar", "Green Algal Bloom", "Blue-Green Phytoplankton", "sargassum polycystum", " marine Phaeophyta", "brown algae", "nereocystis luetkeana", "subtidal kelp", "protoceratium reticulatum", "dinoflagellate", "asterionella formosa", "Skeletonema costatum", "thalassiosira pseudonana", "thalassiosira oceanica", "thalassiosira weissflogii", "thalassiosira rotula", "thalassiosira gravida", "thalassiosira nordenskioeldii", "thalassiosira punctigera", "thalassiosira decipiens", "thalassiosira delicatula", "thalassiosira symmetrica", "thalassiosira subtilis", "thalassiosira partheneia", "thalassiosira lacustris", "thalassiosira minima", "Haptophytes", "Phaeocystis globosa", "Phaeocystis pouchetii", "Phaeocystis antarctica", "Phaeocystis globosa", "Phaeocystis pouchetii", "Phaeocystis antarctica", "Chrysophytes", "Synura", "Dinobryon", "Ochromonas", "Mallomonas", "Chrysosphaerella", "Chrysococcus", "Chrysodidymus", "Chrysonebula", "Chrysotilos", "Chrysocapsa", "Chrysocystis", "Chrysocystis fragilis", "Chrysocystis minor", "Chrysocystis major", "Chrysocystis parva", "pavlova", "pavlova lutheri", "pavlova salina", "pavlova sp.", "pavlova sp. lutheri", "pavlova sp. salina", "pavlova sp. lutheri", "pavlova sp. salina", "pavlova sp. lutheri", "pavlova sp. salina", "jania sp."],
            
            "Epidermis/Bark": ["plant epidermis","Barberry Bark","Barberry","Bark", "Tree Bark", "Woody Bark", "Tree Trunk Bark", "Plant Bark", "Cork Bark", "Outer Bark", "Bark Sample", "Birch Bark", "Oak Bark", "Pine Bark", "Maple Bark", "Cedar Bark", "Willow Bark", "Eucalyptus Bark", "Cinnamon Bark", "Cherry Bark", "Bamboo Bark", "Rubber Tree Bark", "Mango Bark", "Teak Bark", "Mahogany Bark", "Alder Bark", "Redwood Bark", "Poplar Bark", "Chestnut Bark", "Beech Bark", "Birchwood", "Carob Bark", "Juniper Bark", "Linden Bark", "Larch Bark", "Sequoia Bark", "Buckthorn Bark", "Sweet Birch Bark", "Black Cherry Bark", "Bamboo Sheath", "Hickory Bark", "Sandalwood Bark", "Camphor Bark", "Neem Bark", "Sassafras Bark", "Lemon Bark", "Palo Santo Bark", "Mesquite Bark", "Prickly Ash Bark", "Sandpaper Tree"],
            
            "Desert Plant" : ["carex pumila", "desert plant", "desert plant sample", "desert plant isolate", "desert plant mucilage", "desert plant exudate", "desert plant exudates", "desert plant exudate sample", "desert plant exudate isolate", "desert plant material", "plant from desert", "plant of moving sand dunes"],
            
            "Lichen" : ["cladonia cristatella", "lichen"], 
            
            "Sap": ["birch sap", "sap", "tree sap"],
            
            "Slime Mold" : ["slime mold"],
            
            "Mangrove" : ["avicennia alba", "mangrove", "aegiceras corniculatum", "black mangrove", "river mangrove"], 
            
            "Gall" : ["tumor", "crown gall", "cane gall", "crown gall tissues", "gall", "stem gall tissue", "stem gall", "tumor tissue"],
            
            "Others" : ["Plant material", "plants", "plant sample", "plant isolate", "plant exudate", "plant exudates", "plant exudate sample", "plant exudate isolate", "plant tissues", "populus x jackii", "plantation", "powder plant", "alectra sessiliflora", "pothos", "heliconia", "anthurium", "rhaphiolepis umbellata", "pyrus", "prunus avium", "broussonetia kazinoki", "castanea crenata", "prunus yedoensis", "daphniphyllum teijsmannii", "aesculus hippocastanum", "dendropanax trifidus", "eriobotrya japonica", "morella rubra", "fraxinus excelsior", "nerium oleander", "medicago", "medicago plant", "medicago falcata", "solanum tuberosum", "glycine max", "wenyujin", "polyalae radix", "oryza sativa", "callitris preissii", "native pine tree", "jasmin rice", "jatropha curcas", "Sedum sp.", "maple", "eucalyptus", ],
            
            "Grass/Grassland": ["grass", "grassland", "grass sample", "grass isolate", "grass mucilage", "grass exudate", "grass exudates", "grass exudate sample", "grass exudate isolate", "grass material", "grasses", "grasses of the family Poaceae", "grasses of the family Cyperaceae", "canadian high arctic grass", "arctic grass", "arctic grass sample", "arctic grass isolate", "arctic grass mucilage", "arctic grass exudate", "arctic grass exudates", "arctic grass exudate sample", "arctic grass exudate isolate", "arctic grass material", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Cyperaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae", "grasses of the family Cyperaceae", "grasses of the family Poaceae", "grasses of the family Juncaceae", "grasses of the family Restionaceae"],
            
            "Cactus": ["cactus", "cactus sample", "cactus isolate", "cactus mucilage", "cactus exudate", "cactus exudates", "cactus exudate sample", "cactus exudate isolate", "cactus material", "cacti"],

            "Medicinal Plants": ["medicinal plant", "medicinal plant sample", "medicinal plant isolate", "medicinal plant mucilage", "medicinal plant exudate", "medicinal plant exudates", "medicinal plant exudate sample", "medicinal plant exudate isolate", "medicinal plant material", "humulus lupulus"],

            "Herbal Medicine": ["herbal medicine", "herbal medicine sample", "herbal medicine isolate", "chinese herbal medicine", "asian herbal medicine", "traditional herbal medicine", "herbal remedy", "herbal extract", "herbal supplement", "herbal formulation", "herbal product", "herbal preparation", "herbal tea", "herbal tincture", "herbal infusion", "herbal decoction", "herbal poultice", "herbal ointment", "herbal balm", "herbal salve", "herbal syrup", "herbal capsule", "herbal tablet", "indian herbal medicine", "western herbal medicine", "herbal medicine formulation", "herbal medicine preparation", ],
        },

        "Solidwaste": {
            "Others": ["Marine sludge","Waste", "Garbage", "Trash", "Rubbish", "Refuse", "Scrap", "Discarded Materials", "General Waste", "Unsorted Waste", "Unclassified Waste", "Uncategorized Refuse", ],
            "Treatment Plant" : ["solid waste treatment site", "solid waste treatment plant", "solid waste treatment facility", ],
            "Municipal Solid Waste": ["Municipal Solid Waste", "Household Waste", "Residential Waste", "Domestic Waste", "City Trash", "Urban Waste", "City Garbage", "Neighborhood Waste", "Street Trash", "Public Waste", "Community Waste", "Municipal Refuse", "City Refuse", "Town Waste", "Household Plastic Waste", "Plastic Packaging", "Plastic Bottles", "Plastic Bags", "plastic"],
            "Industrial And Hazardous Waste": ["Industrial Waste", "Factory Waste", "Manufacturing Waste", "Chemical Waste", "Polluted Industrial Waste", "Toxic Industrial Waste", "Mining Waste", "Industrial Byproducts", "Production Waste", "Waste From Factories", "Metal Scrap Waste", "Chemical Byproducts", "Industrial Residues", "Plastic Waste From Factories", "Plastic Manufacturing Waste", "Textile Waste", "Fabric Waste", "Clothing Waste", "Discarded Textiles", "Worn-Out Clothes", "Clothing Scraps", "Textile Byproducts", "Old Clothing", "Damaged Textiles", "Fashion Waste", "Used Fabric", "Textile Disposal", "Scrapped Textiles", "Waste Fabric", "Second-Hand Clothing Waste", "Hazardous Waste", "Toxic Waste", "Dangerous Waste", "Flammable Waste", "Radioactive Waste", "Poisonous Waste", "Contaminated Waste", "Corrosive Waste", "Biohazard Waste", "Explosive Waste", "Pollutant Waste", "Hazardous Material Waste", "Dangerous Goods Waste", "Dangerous Chemicals Waste", "tannery sludge"],
            "Food Waste": ["Food Waste", "Organic Waste", "Kitchen Waste", "Leftover Food", "Spoiled Food", "Food Scraps", "Compostable Waste", "Discarded Food", "Uneaten Food", "Food Debris", "Waste From Food Processing", "Food Leftovers", "Perishable Waste", "Unwanted Food", "Rotting Food"],
            "Construction And Demolition Waste": ["Construction Waste", "Demolition Waste", "Building Debris", "Construction Debris", "Demolition Rubble", "Building Materials Waste", "Concrete Waste", "Metal Scrap", "Renovation Waste", "Site Waste", "Construction Leftovers", "Demolition Materials", "Building Rubble", "Construction Site Waste"],
            "Agricultural Waste": ["Agricultural Waste", "Farm Waste", "Crop Residue", "Animal Manure", "Plant Debris", "Harvest Waste", "Crop Byproducts", "Agricultural Runoff", "Pesticide Containers", "Fertilizer Packaging", "Agriculture-Related Waste", "Agricultural Chemicals Waste", "Agro-Waste", "Farm Animal Waste", "Waste From Farming","Sow manure",],
            "Hazardous Waste": [],
            "Wood Waste": ["Wood Waste", "Wooden Debris", "Scrap Wood", "Wood Scraps", "Wooden Pallets", "Old Timber", "Discarded Wood Products", "Sawdust", "Wood Chips", "Wooden Packaging Waste", "Scrapped Wood", "Unused Timber", "Waste Wood Planks", "Wooden Material Waste"],
            "Dumpyard": ["Dump Yard", "Landfill", "Garbage Dump", "Waste Disposal Site", "Refuse Dump", "Rubbish Dump", "Trash Dump", "Waste Management Facility", "Landfill Site", "Garbage Disposal Site", "Refuse Disposal Site", "Rubbish Disposal Site", "Trash Disposal Site"],
        },

        "Biogas/Wastewater-related": {
            "Biogas Plant": ["biogas plant", "anaerobic digestion plant", "anaerobic digestion facility", "biogas facility", "biogas production plant", "anaerobic digestion system", "biogas generation plant", "biogas production facility", "anaerobic digestion unit", "biogas unit", "anaerobic digestion reactor", "biogas reactor", "anaerobic digester", "biogas digester", "anaerobic sludge of molasses wastewater in a continuous stirred-tank reactor", "Anaerobic Digestion", "Biogas Plant", "Biogas Reactor", "Biogas Fermentation", "Biogas Generation", "Biogas System", "Biogas Facility", "Anaerobic Digestion Sludge", "Anaerobic Digestion Residue", "Anaerobic Digestion Byproducts", "Anaerobic Digestion Waste", "Biogas Treatment Plant", "Biogas Production Facility", "Anaerobic Digestion Process", "Anaerobic Digestion Technology",],
            "Bioreactors" : ["Bioreactors", "Anaerobic Bioreactor", "Anaerobic Digestion Bioreactor", "Biogas Bioreactor", "Anaerobic Reactor", "Biogas Reactor", "Anaerobic Digestion Reactor", "Biogas Production Reactor", "Anaerobic Digestion System", "Biogas System", "Anaerobic Digestion Unit", "Biogas Unit", "Anaerobic Digester", "Biogas Digester", "Anaerobic Treatment Reactor", "Biogas Treatment Reactor", "Anaerobic Fermentation Reactor", "Biogas Fermentation Reactor", "CSTRs", "UASB Reactors", "Anaerobic Fluidized Bed Reactor", "Anaerobic Sequencing Batch Reactor", "Anaerobic Packed Bed Reactor", "Anaerobic Hybrid Reactor", "Anaerobic Membrane Bioreactor", "Anaerobic Fixed Film Reactor", "Anaerobic Contact Reactor", "Anaerobic Expanded Bed Reactor", "Anaerobic Sludge Blanket Reactor", "Anaerobic Continuous Stirred-Tank Reactor",],
            "Wastewater Sludge": ["Wastewater Sludge", "Sewage Sludge", "Sewage Waste", "Treated Sludge", "Wastewater Treatment Byproducts", "Sludge From Sewage Plants", "Effluent Sludge", "Sewage Treatment Sludge", "Municipal Sludge", "Organic Sludge", "Industrial Sludge", "Sludge From Water Treatment", "Biosolids", "Wastewater Solids", "Treatment Plant Waste", "isolated from sludge collected at the drain outlet", "Waste at lagoon stage", "activated sludge", "human sludge", "sewage sludge", "sludge from wastewater treatment plant", "sludge from wastewater treatment", "sludge from sewage treatment plant", "sludge from sewage treatment", "anaerobic activated sludge of molasses wastewater in a continuous stirred-tank reactor", "anaerobic activated sludge of molasses wastewater", "anaerobic activated sludge", "anaerobic sludge", "anaerobic sludge of molasses wastewater", "anaerobic sludge of molasses wastewater",  "anaerobic sludge of molasses wastewater"],
            "Reactor Sludge": ["Reactor Sludge", "Bioreactor Sludge", "Anaerobic Reactor Sludge", "Biogas Reactor Sludge", "Anaerobic Digestion Reactor Sludge", "Biogas Production Reactor Sludge", "Anaerobic Digestion System Sludge", "Biogas System Sludge", "Anaerobic Digestion Unit Sludge", "Biogas Unit Sludge", "Anaerobic Digester Sludge", "Biogas Digester Sludge", "Anaerobic Treatment Reactor Sludge", "Biogas Treatment Reactor Sludge", "Anaerobic Fermentation Reactor Sludge", "Biogas Fermentation Reactor Sludge", ],
            "Biogas" : ["Biogas", "biomethane", "green gas", "biomethane gas", "biogas fuel", ],
            "Biogas " : ["Biogas Byproducts", "Biogas Byproduct Residue", "Biogas Byproduct Waste", "Biogas Byproduct Slurry", "Biogas Byproduct Effluent", "Biogas Effluent", "Biogas Slurry", "Biogas Residue", "Biogas Waste", "Biogas Digestate", "Biogas Residuals", "Biogas Byproduct Treatment", "Biogas Byproduct Management", "Biogas Byproduct Utilization", "Biogas Byproduct Recycling", "Biogas Byproduct Disposal", "Biogas Byproduct Processing", "Biogas Byproduct Handling", "Biogas Byproduct Storage", "Biogas Byproduct Separation", "Biogas Byproduct Extraction", "Biogas Byproduct Analysis", "Biogas Byproduct Characterization", "Biogas Digestate", ],
            "Wastewater Treatment Plant": ["wastewater treatment plant", "sewage treatment plant", "wastewater treatment facility", "sewage treatment facility", "wastewater treatment site", "sewage treatment site", "wastewater treatment system", "sewage treatment system", "wastewater treatment unit", "sewage treatment unit", "wastewater treatment process", "sewage treatment process", "wastewater treatment technology", "sewage treatment technology"],
        },


        "Industrial":{
            "Unknown" : ["isolate from industry", "industry isolate", "industry sample", "isolate from industry sample", "industry sample isolate", "isolate from industry isolate", "industry isolate sample", "isolate from industry isolate sample", "industry isolate isolate", "isolate from industry isolate isolate", "factory"],
            "Equipments" : ["Equipments", "Industrial Equipment", "Machinery", "Industrial Machinery", "Manufacturing Equipment", "Production Equipment", "Processing Equipment", "Heavy Machinery", "Industrial Tools", "Mechanical Equipment", "Industrial Devices", "Manufacturing Machinery", "almond hulling equipment", "almond hulling machine", "almond hulling machines", "almond hulling process", "almond hulling processes", "almond hulling system", "almond hulling systems", "almond hulling technology", "almond hulling technologies", "almond hulling unit", "almond hulling units", "almond hulling line", "almond hulling lines", "pharmaceutical production line", "Industry equipment"],
            "Mines" : ["gold and copper mine", "coal mine", ],
            "Industrial Products": ["amylase concentrate", "cellulase concentrate", "galactooligosaccharide used in manufacture"],
            "Oil Fiels" : ["oil field", "oil refineries", ],
            "Tanneries": ["tannery", "tannery bath containers", "tannery lime containers", ],
            "Cement" : ["cement", "cement plant", "cement factory", "cement production site", "cement manufacturing facility", "cement processing plant", "cement production plant", "cement manufacturing site", "cement processing facility", "cement production facility", "cement manufacturing plant", "cement processing site"],
            "Raw Material" : ["raw material", "industrial raw material", "manufacturing raw material", "production raw material", "processing raw material", "industrial feedstock", "manufacturing feedstock", "production feedstock", "processing feedstock", "industrial input", "manufacturing input", "production input", "processing input"],
            "Effluent Treatment Plant" : ["effluent treatment plant", "effluent treatment facility", "effluent treatment site", "industrial effluent treatment plant", "industrial effluent treatment facility", "industrial effluent treatment site", "wastewater effluent treatment plant", "wastewater effluent treatment facility", "wastewater effluent treatment site"],
            
        },

        "Rocks" : {
            "Sedimentary Rocks": ["Sedimentary Rock", "Sandstone", "Shale", "Limestone", "Conglomerate", "Breccia", "Siltstone", "Mudstone", "Chalk", "Gypsum", "Halite", "Dolomite", "Coal", "Oil Shale", "Phosphorite", "Chert"],
            "Igneous Rocks": ["Igneous Rock", "Granite", "Basalt", "Andesite", "Diorite", "Gabbro", "Rhyolite", "Obsidian", "Pumice", "Scoria", "Volcanic Rock"],
            "Metamorphic Rocks": ["Metamorphic Rock", "Marble", "Slate", "Gneiss", "Schist", "Quartzite", "Phyllite", "Amphibolite", "Eclogite", "Hornfels"],
            "Minerals": ["Mineral", "Quartz", "Feldspar", "Mica", "Calcite", "Gypsum", "Halite", "Dolomite", "Pyrite", "Hematite", "Magnetite", "Sphalerite", "Galena", "Bauxite", "Fluorite", "Apatite"],
            "Fossils": ["Fossil", "Petrified Wood", "Amber", "Coprolite", "Gastrolith", "Trace Fossil", "Body Fossil", "Microfossil", "Fossilized Shell", "Fossilized Bone", "Fossilized Plant Material"],
            "Volcanic Rocks": ["Volcanic Rock", "Basalt", "Andesite", "Diorite", "Gabbro", "Rhyolite", "Obsidian", "Pumice", "Scoria"],
            "Other Rocks": ["Rock", "Stone", "Pebble", "Boulder", "Gravel", "Cobble", "Rock Fragment", "Rock Sample", "Rock Material", "Rock Specimen", "Rock Outcrop", "Rock Formation", "Rock Layer", "Rock Stratum", "Rock Bed", "Rock Deposit"],},
            
        "Extreme Environments" : {
            "Coal bed": ["Coal bed", "coal seam", "coal deposit", "carboniferous bed", "peat layer", "bituminous coal zone"],
            "Deep sea": ["Deep sea", "abyss", "oceanic trench", "hadal zone", "mariana trench", "deep ocean", "benthic zone"],
            "Hydrothermal vent": ["Hydrothermal vent", "black smoker", "white smoker", "seafloor vent", "deep-sea vent", "geothermal vent", "marine hydrothermal chimney wall", "marine hydrothermal sulfide sediment", "marine hydrothermal sulfide sediment",],
            "Hot spring": ["Hot spring", "thermal spring", "geyser", "boiling spring", "fumarole", "geothermal spring", 'hot underground spring', 'hot water spring'],
            "Acid mine drainage": ["Acid mine drainage", "AMD", "acidic runoff", "mine water", "sulfuric drainage", "contaminated mine water"],
            "Desert": ["Desert", "arid zone", "drought region", "sand dunes", "semi-arid land", "barren land"],
            "Antarctic": ["Antarctic", "polar region", "south pole", "ice shelf", "glacial zone", "permafrost area"],
            "Deep subsurface": ["Deep subsurface", "deep biosphere", "geothermal subsurface", "rock pore system", "subterranean habitat"],
            "High radiation zone": ["High radiation zone", "radioactive site", "chernobyl exclusion zone", "gamma-ray exposed area", "nuclear waste site"],
            "Permafrost": ["Permafrost", "frozen soil", "tundra", "ice-bound ground", "periglacial zone", "cryotic soil"],
            "Alkaline soil": ["Alkaline soil", "soda soil", "high pH soil", "basic soil", "calcareous soil", "salt-affected soil"],
            "Volcanic crater": ["Volcanic crater", "lava crater", "caldera", "active volcano", "magma chamber", "pyroclastic zone"],
            "Volcanic ash": ["Volcanic ash", "tephra", "volcanic dust", "ash fall", "pyroclastic material", "volcanic eruption residue", "Volcanic soil"],
            "Black shale": ["Black shale", "organic-rich shale", "kerogen rock", "carbonaceous shale", "petroleum source rock"],
            "Petroleum reservoir": ["Petroleum reservoir", "oil field", "crude oil deposit", "hydrocarbon reservoir", "tar sand", "bitumen field", "oil sands", "tar sands"],
            "Soda lake": ["Soda lake", "alkaline lake", "high pH water body", "sodium carbonate lake", "brine lake", "Hypersaline lake", "salt lake", "brine pool", "saline lagoon", "soda lake", "alkaline lake"],
            "Extreme high altitude samples": ["Extreme high altitude", "mountain peak", "oxygen-deprived zone", "Himalayan plateau", "andes summit"],
            "Deep cave": ["Deep cave", "karst cave", "limestone cave", "subterranean cavern", "underground chamber"],
            "Endolithic" : ["Endolithic", "Endolithic Microorganisms", "Endolithic Community", "Endolithic Diversity", "Endolithic Ecology", "Endolithic Analysis", "Endolithic Structure", "Endolithic Function", "Endolithic Interactions", "Endolithic Dynamics", "Endolithic Composition"],
        },
        "Residential Areas":{
            "Household": ["Household", "Home", "Residential Area", "Domestic Environment", "Living Space", "Private Residence", "Family Home", "Apartment", "Condominium", "Townhouse", "Single-Family Home", "Duplex", "Bungalow", "Cottage", "Villa", "Mansion", "Ranch House", "Mobile Home", "Tiny House", "residential environment"],
            "Care Facility": ["Care Facility", "Nursing Home", "Assisted Living Facility", "Long-Term Care Facility", "Senior Living Community", "Retirement Home", "Skilled Nursing Facility", "Memory Care Facility", "Rehabilitation Center", "Group Home", "Residential Care Facility", "Personal Care Home", "Elderly Care Facility", "Adult Foster Care", "Continuing Care Retirement Community", "Senior Care Center", "Alzheimer's Care Facility", "Palliative Care Facility", "Hospice Care Facility", "Respite Care Facility", "Transitional Care Facility", "Subacute Care Facility", "Independent Living Facility", "Senior Housing Community", "Assisted Living Residence", "Nursing Care Facility", "Long-Term Care Home", "Elderly Residential Facility", "Senior Assisted Living Facility", "Skilled Nursing Home", "Memory Support Facility", "Rehabilitation Hospital", "Group Living Facility", "Residential Rehabilitation Facility", "Personal Care Residence", "Elderly Group Home", "Adult Day Care Facility", "Continuing Care Facility", "Senior Living Residence", "Alzheimer's Care Home", "Palliative Care Home", "Hospice Residence", "Respite Care Home", "Transitional Care Home", "Subacute Care Home", "Independent Living Residence", "Senior Housing Facility", "Assisted Living Community", "residential aged care facility"],            
        },

        "Personal Care Products": {
            "Grooming Accessories": ["Hair brushes", "Hairbrush", "Combs", "Hair comb", "Nail clippers", "Nail cutter", "Manicure tools", "Pedicure kit", "Tweezers", "Pluckers", "Shaving kits", "Shave set", "Shaving set", "Razors", "Shavers", "Electric razors", "Gillette", "Wilkinson Sword", "Beard trimmers", "Hair trimmers", "unused razor head", "unused razor", "unused razor blade", "unused razor blades", "unused razor head sample", "unused razor head isolate"],
            "Make-up Kits": ["Foundation sets", "Base makeup", "Lipstick combos", "Lip colors", "Eye shadow palettes", "Eyeshadow", "Blush & contour kits", "Cheek palette", "Blusher", "Contour powder", "Make-up brushes", "Brush sets", "Blending brushes", "Make-up sponges", "Beauty blender", "Makeup puff"],
            "Skincare Products": ["Face wash", "Facial cleanser", "Moisturizers", "Face cream", "Lotion", "Sunscreen", "Sunblock", "SPF lotion", "Exfoliators", "Face scrub", "Scrubbing gel", "Facial masks", "Sheet masks", "Face packs",],
            "Hair Care Products": ["Shampoos", "Hair shampoo", "Head wash", "Conditioners", "Hair conditioner", "Silkening cream","Hair oils", "Coconut oil", "Argan oil", "Deep conditioning mask", "Hair styling gels", "Hair wax", ],
            "Oral Care": ["Toothbrushes", "Brush", "Electric toothbrush", "Toothpaste", "Tooth gel", "Fluoride paste", "Mouthwash", "Oral rinse", "Dental floss", "Flossing thread", "Floss picks","Teeth whitening kits", "Whitening strips", "Charcoal toothpaste"],
            "Bath & Body": ["Body wash", "Shower gel", "Shower cream", "Soaps", "Bar soap", "Toilet soap","Bath salts", "Soaking salt","Body scrubs", "Exfoliating scrub", "Sugar scrub","Loofahs", "Bath sponge", "Body puff",],
            "Fragrances": ["Perfumes", "Eau de parfum", "Body spray", "Fragrance mist", "Bath & Body Works","Deodorants", "Colognes", "Aftershave cologne", ],
            "Women's Care": ["Feminine hygiene products", "Sanitary pads", "Tampons", "Menstrual cups", "Whisper", "Stayfree","Make-up removers", "Cleansing oil", "Micellar water", ],
            "Baby Care" : ["Baby shampoo", "Baby wash", "Baby lotion", "Baby oil", "Diaper rash cream", "Baby wipes", "Baby powder", "Baby sunscreen", "Baby soap", "Baby moisturizer", "Baby wipes", ],
        },

        "Personal Iterms": {
            "Mobile Phones": ["Mobile Phones", "Smartphones", "Cell Phones", "Handsets", "Mobile Devices", "Mobile Communication Devices", "Mobile Telephones", "Mobile Phone Accessories", "Mobile Phone Cases", "Mobile Phone Chargers", "Mobile Phone Screens", "Mobile Phone Batteries"],
            "Laptops": ["Laptops", "Notebooks", "Ultrabooks", "Gaming Laptops", "Business Laptops", "Convertible Laptops", "2-in-1 Laptops", "Laptop Accessories", "Laptop Cases", "Laptop Chargers", "Laptop Screens", "Laptop Batteries"],
            "Tablets": ["Tablets", "iPads", "Android Tablets", "Windows Tablets", "Tablet Accessories", "Tablet Cases", "Tablet Chargers", "Tablet Screens", "Tablet Batteries"],

        },

        "Public Spaces" : {
            "Parks": ["Parks", "Public Parks", "City Parks", "Urban Parks", "Community Parks", "Recreational Parks", "National Parks", "State Parks", "Botanical Gardens", "Playgrounds", "Green Spaces", "Parkland"],
            "Public Transportation": ["Public Transportation", "Buses", "Trains", "Subways", "Light Rail", "Streetcars", "Ferries", "Public Transit Systems", "Mass Transit", "Public Transit", "Public Transportation Systems", "Public Transit Networks", "Public Transportation Services", "Public Transit Routes", "Public Transportation Infrastructure", "Public Transit Vehicles", "Public Transportation Facilities", "Public Transit Stations", "Public Transportation Hubs", "Public Transit Stops", "Public Transportation Schedules", "Public Transit Maps", "Public Transportation Accessibility"],
            "Shopping Centers": ["Shopping Centers", "Malls", "Retail Complexes", "Commercial Centers", "Shopping Malls", "Retail Centers", "Shopping Plazas", "Commercial Complexes", "Shopping Districts", "Retail Parks", "Shopping Arcades", "Commercial Districts", "Shopping Streets", "Retail Streets", "Shopping Avenues", "Commercial Streets", "Shopping Boulevards", "Retail Boulevards"],
            "Streets": ["Streets", "Roads", "Avenues", "Boulevards", "Lanes", "Highways", "Public Roads", "City Streets", "Urban Roads", "Residential Streets", "Commercial Streets", "Rural Roads", "Suburban Streets", "Street Corners", "Street Intersections", "Street Crossings", "Street Medians", "Street Sidewalks"],
            "Public Buildings": ["Public Buildings", "Government Buildings", "Libraries", "Museums", "Community Centers", "City Halls", "Courthouses", "Public Schools", "Public Hospitals", "Public Libraries", "Public Museums", "Public Community Centers", "Public Government Buildings", "Public City Halls", "Public Courthouses", "Public Schools", "Public Hospitals"], 
            "Beach Utilities": ["Beaches", "Sandy Beaches", "Coastal Beaches", "Seaside Beaches", "Ocean Beaches", "Beachfronts", "Beach Resorts", "Beach Parks", "Beach Areas", "Beach Zones", "Beachfront Properties", "Beachfront Locations", "Beachfront Communities", "Beachfront Destinations", "Beachfront Attractions", "Beachfront Activities", "Beach showers"],
        },

        "Biofouling": {
            "Ship Hulls": ["Ship Hulls", "Boat Hulls", "Vessel Hulls", "Marine Vessel Hulls", "Ship Bottoms", "Boat Bottoms", "Vessel Bottoms", "Marine Vessel Bottoms", "Ship Keels", "Boat Keels", "Vessel Keels", "Marine Vessel Keels"],
            "Marine Structures": ["Marine Structures", "Offshore Platforms", "Underwater Pipelines", "Subsea Cables", "Marine Docks", "Underwater Structures", "Offshore Wind Turbines", "Marine Buoys", "Underwater Foundations", "Offshore Oil Rigs", "Marine Bridges", "Underwater Piers"],
            "Aquaculture Equipment": ["Aquaculture Equipment", "Fish Farm Equipment", "Shellfish Farm Equipment", "Aquaculture Structures", "Fish Farm Structures", "Shellfish Farm Structures", "Aquaculture Cages", "Fish Farm Cages", "Shellfish Farm Cages", "Aquaculture Nets", "Fish Farm Nets", "Shellfish Farm Nets"],
            "Water Intake Systems": ["Water Intake Systems", "Cooling Water Systems", "Desalination Plant Intake Systems", "Power Plant Cooling Water Systems", "Industrial Water Intake Systems", "Municipal Water Intake Systems", "Water Treatment Plant Intake Systems", "Hydroelectric Power Plant Intake Systems", "Water Distribution System Intake Systems"],
            "Marine Sensors": ["Marine Sensors", "Underwater Sensors", "Oceanographic Sensors", "Marine Monitoring Sensors", "Underwater Monitoring Sensors", "Oceanographic Monitoring Sensors", "Marine Environmental Sensors", "Underwater Environmental Sensors", "Oceanographic Environmental Sensors"],
            
        }
    }

    best_match_score = 0
    best_domain = "Unknown"
    best_source = "Unknown"

    try:
        if isinstance(fuzzy_str, str):
            logging.debug(f"Entered the str loop in find_environmental_source function with string: `{fuzzy_str}`.")
            
            part_digit_free = re.sub(r'\d+', '', fuzzy_str)
            part = re.sub(r'[\[\]{}()]', '', part_digit_free)
            part = re.sub(r'"', '', part)
            part = re.sub(r'\'', '', part)
            logging.debug(f"Normalized string: {part}")
            space_count = part.count(' ')
            if space_count > 3:
                parts = re.split(r'[-+:;/_,\s]|\band\b', part)
                parts = [part.strip() for part in parts]
                splitted = True
            else:
                parts = [part]
                splitted = False
            
            logging.debug(f"Split string: {parts}")

            for normalized_str in parts:
                normalized_str = normalized_str.lower()
                logging.debug(f"Normalized string: {normalized_str}")

                for environmental_domain, sources in environmental_sphere.items():
                    for source, fuzzy_source_names in sources.items():
                        for fuzzy_source_name in fuzzy_source_names:
                            similarity_score = fuzz.ratio(normalized_str.lower(), fuzzy_source_name.lower())
                            if similarity_score > best_match_score:
                                best_match_score = similarity_score
                                best_domain = environmental_domain
                                best_source = source
                                best_matched_str = fuzzy_source_name
                                normalized_str_matched = normalized_str
                            else:
                                continue

            if best_match_score >= 73:
                logging.debug(f"Environmental source identified for string: `{fuzzy_str}`. The identified environment is `{best_source}` in category `{best_domain}` with a score of {best_match_score}. The best matched string is `{best_matched_str}` with normalized string `{normalized_str_matched}`.")
                return best_domain, best_source, best_match_score

            if not splitted:
                logging.debug(f"Couldn't find any match with environmental source with non-splitted string. Trying with the splitted string.")
                parts = re.split(r'[-+:;/_,\s]|\band\b', part)
                parts = [part.strip() for part in parts]

                logging.debug(f"Split string: {parts}")

                for normalized_str in parts:
                    normalized_str = normalized_str.lower()
                    logging.debug(f"Normalized string: {normalized_str}")

                    for environmental_domain, sources in environmental_sphere.items():
                        for source, fuzzy_source_names in sources.items():
                            for fuzzy_source_name in fuzzy_source_names:
                                similarity_score = fuzz.ratio(normalized_str.lower(), fuzzy_source_name.lower())
                                if similarity_score > best_match_score:
                                    best_match_score = similarity_score
                                    best_domain = environmental_domain
                                    best_source = source
                                    best_matched_str = fuzzy_source_name
                                    normalized_str_matched = normalized_str
                                else:
                                    continue

            if best_match_score >= 73:
                logging.debug(f"Environmental source identified for string: `{fuzzy_str}`. The identified environment is `{best_source}` in category `{best_domain}` with a score of {best_match_score}. The best matched string is `{best_matched_str}` with normalized string `{normalized_str_matched}`.")
                return best_domain, best_source, best_match_score


            else:
                logging.debug(f"The string: `{fuzzy_str}` didn't match any entry from the database. Highest match was with `{best_source}` with {best_match_score} score. Returning `Unknown`.")
                return "Unknown", "Unknown", 0
        
        else:   
            logging.debug(f"Given string: `{fuzzy_str}` is not a string, returning `Unknown`.")
            return "Unknown", "Unknown", 0
    
    except Exception as e:
        logging.error(f"An unexpected error occurred while finding environmental source of string `{fuzzy_str}`: {e}")
        return "Unknown", "Unknown", 0

def get_normalized_isolation_source(combined_columns_list:list):

    logging.debug(f"Startiing the normalization of {combined_columns_list}")
    try:
        host_is, host_match_score = "Unknown", 0

        if isinstance(combined_columns_list, list):

            try:
                for item in combined_columns_list:
                    logging.debug(f"Trying to identify host from {combined_columns_list} ")
                    host_is, host_match_score = find_the_host(fuzzy_str=item)
                    if host_is != "Unknown":
                        break
                
                logging.debug(f"Found the host: {host_is} with {host_match_score} score. Moving forward with further categorization.")

                if combined_columns_list in ([None, None, 'urine'], [None, None, 'Urine']):
                    logging.debug(f"Urine is given without any other entries. Saving it from identifying as animal source as it matches with `Murine`")
                    host_is, host_match_score = "Unknown", 0
                    logging.debug(f"Done with urine. Returning `Unknown`.")

                
                if host_is == "Hospital-associated":
                    try:
                        logging.debug(f"Trying to identify pediatric samples!")
                        pediatric_keywords = {'pediatric', 'paediatric', 'child', 'children', 'kids', 'kid','infant', 'infants', 'neonate', 'neonates', 'neonatal', 'newborn', 'newborns', 'baby', 'babies', 'toddler', 'toddlers', 'adolescent', 'adolescents', 'preterm', 'premature', 'preemie', 'preemies', 'nicu', 'neonatal icu', 'neonatal intensive care unit', 'pediatric icu', 'pediatric intensive care unit', 'pediatric ward', 'pediatric unit', 'pediatric hospital', 'children hospital', 'child hospital', 'pediatric clinic', 'pediatric department', 'pediatric center', 'pediatric care', 'pediatric patient', 'pediatric sample', 'pediatric case', 'pediatric study', 'pediatric research', 'pediatric cohort', 'pediatric population', 'pediatric group', 'Omphalitis', "meconium", 'ductus venosus', 
                        'asphyxia neonatorum', 'hydrops fetalis', 'hydrops', 'neonatalis', 'neonatology'}

                        for item in combined_columns_list:
                            
                            if item is None:
                                continue
                            
                            pediatric_sample = any(keyword in item.lower() for keyword in pediatric_keywords)
                            
                            if pediatric_sample:
                                pediatric_sample = True
                                break
                            else:
                                pediatric_sample = False
                    
                    except Exception as e:
                        logging.error(f"An unexpected error occured while looking for pediatric samples in the given list {combined_columns_list}: {e}")
                        pediatric_sample = False
                        
                    try:
                        for item in combined_columns_list:
                            human_disease, human_disease_score = find_human_disease(fuzzy_str=item)
                            hospital_equipment_category, hospital_equipment, equipment_score = find_hospital_equipment_category_and_name(fuzzy_str=item)
                            if human_disease_score >= 80:
                                break
                            elif equipment_score >= 80:
                                break
                            else:
                                continue
                        
                        if equipment_score >= 80 and human_disease_score < 80:
                            logging.debug(f"Equipment score is {equipment_score} and human disease score is {human_disease_score}.")
                            # just to be on safe side, trying to identify the human disease from the sampple, in cases where we want to avoid samples to be identified as equipments
                            for item in combined_columns_list:

                                human_sample, human_sample_score = find_human_sample_source(fuzzy_str=item)
                                logging.debug(f"Equiment score is {equipment_score} and human sample score is {human_sample_score}.")
                                if int(human_sample_score) == int(equipment_score):
                                    logging.debug(f"The human disease score is equal to equipment score.")
                                    return host_is, "Hospital Environment Sample", hospital_equipment_category, hospital_equipment

                                elif int(human_sample_score) != int(equipment_score):
                                    logging.debug(f"The human sample score is not equal to equipment score. equiment score is {equipment_score} and {human_sample_score}.")
                                    try:
                                        if human_sample_score >=80 and human_sample != "Unknown" and int(human_sample_score) > int(equipment_score):
                                            logging.debug(f"Trying to find the disease from the sample.")
                                            human_disease = find_disease_from_sample(human_sample)
                                            if human_disease != "Unsorted Bacterial infections":
                                                if pediatric_sample:
                                                    logging.debug(f"""
                                                    Returning `Host` : {host_is}, 
                                                    `Source Category` : 'Pediatric Clinical Sample',
                                                    `Identified host disease` : {human_disease},
                                                    `sample` : {human_sample}
                                                    """)
                                                    return host_is, "Pediatric Clinical Sample", human_disease, human_sample
                                                else:
                                                    logging.debug(f"""
                                                    Returning `Host` : {host_is}, 
                                                    `Source Category` : 'Human Clinical Sample',
                                                    `Identified host disease` : {human_disease},
                                                    `sample` : {human_sample}
                                                    """)
                                                    return host_is, "Human Clinical Sample", human_disease, human_sample
                                            else:
                                                logging.debug(f"""
                                                        Returning `Host` : {host_is}, 
                                                        `Source Category` : 'Human Clinical Sample',
                                                        `Identified host disease` : Unknown,
                                                        `sample` : {human_sample}
                                                        """)
                                                if pediatric_sample:
                                                    return host_is, "Pediatric Clinical Sample", "Unknown", human_sample
                                                else:
                                                    return host_is, "Human Clinical Sample", "Unknown", human_sample
                                        else:
                                            continue
                                    except:
                                        logging.error(f"An unexpected error occured while looking for human disease/hospital source for given list {combined_columns_list}")        

                            # for item in combined_columns_list:
                            #     hospital_equipment_category, hospital_equipment, equipment_score = find_hospital_equipment_category_and_name(fuzzy_str=item)
                            #     if equipment_score >= 80:
                            if int(human_sample_score) < int(equipment_score) or human_sample == "Unknown":
                                logging.debug(f"equipment score is {equipment_score}.")
                                logging.debug(f"""
                                        Returning `Host` : {host_is}, `Source Category` : Hospital Envrironment Sample,
                                        `Source category 2` : {hospital_equipment_category}, `Identified host disease` : {hospital_equipment},
                                        """)
                                return host_is, "Hospital Environment Sample", hospital_equipment_category, hospital_equipment
                                # else:
                                #     continue

                            # if equipment_score < 80:
                            #     logging.debug(f"""
                            #                 Returning `Host` : {host_is}, `Source Category` : Hospital Envrironment Sample,
                            #                 `Source category 2` : {hospital_equipment_category}, `Identified host disease` : {hospital_equipment},
                            #                 """)
                            #     return host_is, "Hospital Environment Sample", "Unknown", "Unknown"
                        
                        else:
                            logging.debug(f"The human disease score is higher than equipment score.")
                            if human_disease_score >= 80:
                                logging.debug(f" human disease score is {human_disease_score}.")
                                for item in combined_columns_list:
                                    human_sample, human_sample_score = find_human_sample_source(fuzzy_str=item)
                                    if human_sample_score >=80:
                                        if human_disease == "Unsorted Bacterial infections":
                                            logging.debug(f"Trying to find the disease from the sample.")
                                            human_disease = find_disease_from_sample(human_sample)

                                            if human_disease != "Unsorted Bacterial infections":
                                                logging.debug(f"""
                                                Returning `Host` : {host_is},
                                                `Source Category` : Human Clinical Sample,
                                                `Identified host disease` : {human_disease},
                                                `sample` : {human_sample}
                                                """)
                                                if pediatric_sample:
                                                    return host_is, "Pediatric Clinical Sample", human_disease, human_sample
                                                else:
                                                    return host_is, "Human Clinical Sample", human_disease, human_sample,
                                            else:
                                                logging.debug(f"""
                                                Returning `Host` : {host_is},
                                                `Source Category` : Human Clinical Sample,
                                                `Identified host disease` : {human_disease},
                                                `sample` : {human_sample}
                                                """)
                                                if pediatric_sample:
                                                    return host_is, "Pediatric Clinical Sample", "Unsorted Bacterial infections", human_sample
                                                else:
                                                    return host_is, "Human Clinical Sample", "Unsorted Bacterial infections", human_sample
                                        else:
                                            if pediatric_sample:
                                                return host_is, "Pediatric Clinical Sample", human_disease, human_sample
                                            else:
                                                return host_is, "Human Clinical Sample", human_disease, human_sample
                                    else:
                                        continue
                                
                                if human_sample_score < 80:
                                    logging.debug(f"Couldn't identify the sample based on the database for human clinical sample.")
                                    logging.debug(f"""
                                                Returning `Host` : {host_is}, 
                                                `Source Category` : 'Human Clinical Sample',
                                                `Identified host disease` : {human_disease},
                                                `sample` : "Unknown"
                                                """)
                                    if pediatric_sample:
                                        return host_is, "Pediatric Clinical Sample", human_disease, human_sample
                                    else:
                                        return host_is, "Human Clinical Sample", human_disease, human_sample

                            elif human_disease_score < 80:
                                logging.debug(f" human disease score is {human_disease_score}. and still trying to find the sample and may be disease from the sample.")

                                for item in combined_columns_list:
                                    human_sample, human_sample_score = find_human_sample_source(fuzzy_str=item)

                                    if human_sample_score >=80 and human_sample != "Unknown":
                                        logging.debug(f"Trying to find the disease from the sample.")
                                        human_disease = find_disease_from_sample(human_sample)
                                        if human_disease != "Unsorted Bacterial infections":
                                            logging.debug(f"""
                                                Returning `Host` : {host_is}, 
                                                `Source Category` : 'Human Clinical Sample',
                                                `Identified host disease` : {human_disease},
                                                `sample` : {human_sample}
                                                """)
                                            if pediatric_sample:
                                                return host_is, "Pediatric Clinical Sample", human_disease, human_sample
                                            else:
                                                return host_is, "Human Clinical Sample", human_disease, human_sample
                                        else:
                                            logging.debug(f"""
                                                    Returning `Host` : {host_is}, 
                                                    `Source Category` : 'Human Clinical Sample',
                                                    `Identified host disease` : Unknown,
                                                    `sample` : {human_sample}
                                                    """)
                                            if pediatric_sample:
                                                return host_is, "Pediatric Clinical Sample", "Unknown", human_sample
                                            else:
                                                return host_is, "Human Clinical Sample", "Unknown", human_sample
                                    else:
                                        continue
                                
                                logging.debug(f"testing if the sample is from hospital environment.")

                                # for item in combined_columns_list:
                                #     hospital_equipment_category, hospital_equipment, equipment_score = find_hospital_equipment_category_and_name(fuzzy_str=item)
                                #     if equipment_score >= 80:
                                #         logging.debug(f"equipment score is {equipment_score}.")
                                #         logging.debug(f"""
                                #                 Returning `Host` : {host_is}, `Source Category` : Hospital Envrironment Sample,
                                #                 `Source category 2` : {hospital_equipment_category}, `Identified host disease` : {hospital_equipment},
                                #                 """)
                                #         return host_is, "Hospital Environment Sample", hospital_equipment_category, hospital_equipment
                                #     else:
                                #         continue
                                
                                # if equipment_score < 80:
                                if human_sample_score < 80 or human_sample == "Unknown":
                                    if pediatric_sample:
                                        return host_is, "Pediatric Clinical Sample", "Unknown", "Unknown"
                                    else:
                                        return host_is, "Human Clinical Sample", "Unknown", "Unknown"
                                
                    except Exception as e:
                        logging.debug(f"An unexpected error occured while looking for human disease/hospital source for given list {combined_columns_list}. The error is {e}")
                        if pediatric_sample:
                            return host_is, "Pediatric Clinical Sample", "Unknown", "Unknown"
                        else:
                            return host_is, "Human Clinical Sample", "Unknown", "Unknown"
                
                elif host_is == "Animal-associated":
                    logging.debug(f"Host is animal. Now moving forward...")
                    try:
                        for item in combined_columns_list:
                            animal_category, animal, animal_score = find_animal_category(fuzzy_str=item)
                            if animal_score >= 80:
                                break
                            else:
                                continue

                        if animal_score >= 80:
                            logging.debug(f"animal score is {animal_score}.")
                            for item in combined_columns_list:
                                animal_sample, animal_sample_score = find_animal_sample(fuzzy_str=item)
                                if animal_sample_score >= 80:
                                    logging.debug(f"animal sample score is {animal_sample_score}.")
                                    logging.debug(f"""
                                        Returning `Host`: {host_is}, `animal category`: {animal_category},
                                        `animal` : {animal}, `animal sample` : {animal_sample}
                                            """)
                                    return host_is, animal_category, animal, animal_sample
                                else:
                                    continue 
                            
                            if animal_sample_score < 80:
                                logging.debug(f"animal sample score is {animal_sample_score}.")
                                logging.debug(f"""
                                        Returning `Host`: {host_is}, `animal category`: {animal_category},
                                        `animal` : {animal}, `animal sample` : {animal_sample}
                                            """)
                                return host_is, animal_category, animal, animal_sample
                                
                        elif animal_score < 80:
                            logging.debug(f"animal score is {animal_score}.")
                            
                            for item in combined_columns_list:
                                animal_sample, animal_sample_score = find_animal_sample(fuzzy_str=item)
                                if animal_sample_score >= 80:
                                    logging.debug(f"animal sample score is {animal_sample_score}.")
                                    logging.debug(f"""
                                        Returning `Host`: {host_is}, `animal category`: {animal_category},
                                        `animal` : {animal}, `animal sample` : {animal_sample}
                                            """)
                                    return host_is, animal_category, animal, animal_sample
                                else:
                                    continue 
                            
                            if animal_sample_score < 80:
                                logging.debug(f"animal sample score is {animal_sample_score}.")
                                logging.debug(f"""
                                        Returning `Host`: {host_is}, `animal category`: {animal_category},
                                        `animal` : {animal}, `animal sample` : {animal_sample}
                                            """)
                                return host_is, animal_category, animal, animal_sample
                            
                    except Exception as e:
                        logging.debug(f"An unexpected error occured while looking for animal source for given list {combined_columns_list}. The error is {e}")
                        return host_is, "Unknown", "Unknown", "Unknown"

                elif host_is == "Environment-associated":
                    logging.debug("Host is environment. Now moving forward...")
                    try:
                        for item in combined_columns_list:
                            environment_domain, environment, environment_score = find_environmental_source(fuzzy_str=item)
                            if environment_score >= 80:
                                logging.debug(f"environment score is {environment_score}.")
                                logging.debug(f"""
                                        Returning `Host`: {host_is}, `environment domain`: {environment_domain},
                                        `environment` : {environment}
                                            """)
                                return host_is, environment_domain, environment, environment
                            else:
                                continue
                        
                        if environment_score < 80:
                            logging.debug(f"environment score is {environment_score}.")
                            logging.debug(f"""
                                        Returning `Host`: {host_is}, `environment domain`: {environment_domain},
                                        `environment` : {environment}
                                            """)
                            return host_is, environment_domain, environment, environment
                    except Exception as e:
                        logging.debug(f"An unexpected error occured while looking for environmental source for given list {combined_columns_list}. The error is {e}")
                        return host_is, "Unknown", "Unknown", "Unknown"
                
                elif host_is == "Laboratory-based":
                    logging.debug(f"Host is Laboratory-based. Done!")
                    return host_is, "Unknown", "Unknown", "Unknown"

                elif host_is == "Unknown":
                    logging.debug(f"Algorithm is not able to categorize isolation source given.")
                    return host_is, "Unknown", "Unknown", "Unknown"
                else:
                    return host_is, "Unknown", "Unknown", "Unknown"

            except Exception as e:
                logging.debug(f"An unexpected error occured while going through list for identifying host. The error is {e}.")
                return host_is, "Unknown", "Unknown", "Unknown"
        else:
            return host_is, "Unknown", "Unknown", "Unknown"
    
    except Exception as e:
        logging.debug(f"An unexpected error occured while normalizing the isolation source. The error is {e}.")
        return "Unknown", "Unknown", "Unknown", "Unknown"

def dict_update_from_new_isolation_sources(df:pd.DataFrame, dictionary: dict):

    logging.info(f"Updating the local database to categorize the isolation sources.")
    iso_df_just = df[['host', 'host_disease', 'isolation_source']].drop_duplicates()
    length_of_records = len(iso_df_just)


    counter = 1    
    # for fuzzy_lists_of_isolation_sources in enumerate(tqdm(iso_df_just.values.tolist(), desc="Updating Isolation Source dictionary:")):
    #     fuzzy_tuple = tuple(fuzzy_lists_of_isolation_sources[1])
    #     print(fuzzy_tuple)

    for fuzzy_lists_of_isolation_sources in iso_df_just.values.tolist():
        fuzzy_tuple = tuple(fuzzy_lists_of_isolation_sources)
        # print(fuzzy_tuple)

        logging.info(f"processing: {counter}/{length_of_records}")

      
        if fuzzy_tuple not in dictionary:
            result = get_normalized_isolation_source(fuzzy_lists_of_isolation_sources)
            logging.debug(f"The result from get_normalized isolation source is {result}")
        
            if result is None or len(result) != 4:
                logging.debug(f"The result is either none or the length of result is more than eight.")
                identified_host, source_category, source, sample = [None] * 4
            else:
                identified_host, source_category, source, sample = result
                logging.debug(f"The {result} mapped to specific variables.")

            dictionary[fuzzy_tuple] = {
                    'host': fuzzy_tuple[0],
                    'host_disease': fuzzy_tuple[1], 
                    'isolation_source': fuzzy_tuple[2],
                    'identified_host': identified_host,
                    'source_category' : source_category,
                    'source': source,
                    'sample': sample,
                }
            logging.debug(f"Added new entry for: {fuzzy_tuple}")
            counter += 1
        
        else:
            logging.debug(f"Skipping existing entry: {fuzzy_tuple}")
            counter += 1
    
    values = list(dictionary.values())
    fuzzy_dict_df = pd.DataFrame(values)

    logging.debug(f"Fuzzy_dict for isolation source has been converted to dataframe.")
    logging.debug(f"{fuzzy_dict_df.columns}")
    logging.debug(f"All set! returning the updated dictionary.")
    return fuzzy_dict_df

def map_isolation_sources(df:pd.DataFrame, map_df:pd.DataFrame):
    logging.debug(f"length of dataframe {len(df)} and map_df {len(map_df)}")
    # print(map_df)
    merged_df = df.merge(map_df, on=['host', 'host_disease', 'isolation_source'], how='left')
    logging.debug(f"length of dataframe {len(merged_df)}.")

    return merged_df

def all_normalization_operations(df:pd.DataFrame, saving_file_path:str = None):

    # print("Here")
    logging.debug(f"The length of dataframe for normalization is {len(df)}")

    if saving_file_path is None:
        saving_file_path = os.getcwd()
        logging.debug(f"Since the `saving file path` is `None`, the current working directory is choosen. The WD is {saving_file_path}")
    else:
        logging.debug(f"the current saving path: {saving_file_path}")
    
    
    logging.info(f"Starting normalizing/categorizing various data for visulization and exploration.")
    try:
        logging.info(f"Starting with geographic locations.")
        logging.info(f"Loading the local database of geographical locations for normalization.")
        
        try:
            # print("Here_try1")
            df_country_state_data = pd.read_csv(os.path.abspath('./data/geo_loc_name_transforming_database.tsv'), sep='\t', low_memory=False)
        except UnicodeEncodeError as e:

            logging.error(f"An error occurred while loading the local database of country names and state names for normalization. The error is: {e}")
            try:
                # print("Here_try2")
                logging.info(f"Trying to load the file with encoding: ISO-8859-1")
                df_country_state_data = pd.read_csv(os.path.abspath('./data/geo_loc_name_transforming_database.tsv'), sep='\t', low_memory=False, encoding='ISO-8859-1')
            except Exception as e:
                # print("Here_try3")
                logging.info(f"Loading the local database of country names and state names for normalization from the backup file.")
                df_country_state_data = pd.read_csv(os.path.abspath('./data/geo_loc_name_transforming_database_backup.tsv'), sep='\t', low_memory=False)
            except Exception as e:
                # print("here_try4")
                logging.info(f"Exiting the normalization process.")
                messagebox.showerror('Error!', f"An error occurred while loading the local database of country names and state names for normalization. The error is: {e}")  
                raise e

        country_state_map = {row[0]: (row[1], row[2], row[3], row[4], row[5]) for row in df_country_state_data.itertuples(index=False)}
        logging.debug(f"Local database of country names and state names has been loaded. The length of the database now is {len(country_state_map)}")
        country_state_map = {key: value for key, value in country_state_map.items() if value != (None, None, None, None, None)}
        logging.debug(f"Local database has been stripped of the None values. The length of the database now is {len(country_state_map)}")
        dict_update_from_new_loc(df, country_state_map)
        logging.debug(f"Size of the local database after the update is: {len(country_state_map)}")

        logging.debug(f"Adding the normalized country and state columns to the dataframe for plotting.")
        add_normalized_country_state_columns(df, country_state_map)
        # print("here try7")

        # lowering all the data for uniformity in isolation source
        try:
            logging.debug(f"Lowering the 'model' data for uniformity in isolation source.")
            df['Model'] = df['Model'].apply(lambda x: str(x).lower() if pd.notnull(x) else None)
        except:
            logging.debug(f"Model column not found in the dataframe. Adding empty column for model.")
            df['Model'] = None

        try:
            logging.debug(f"Lowering the 'host' data for uniformity in isolation source.")
            df['host'] = df['host'].apply(lambda x: str(x).lower() if pd.notnull(x) else None)
        except:
            logging.debug(f"host column not found in the dataframe. Adding empty column for host.")
            df['host'] = None
        
        try:
            logging.debug(f"Lowering the 'host_disease' data for uniformity in isolation source.")
            df['host_disease'] = df['host_disease'].apply(lambda x: str(x).lower() if pd.notnull(x) else None)
        except:
            logging.debug(f"host_disease column not found in the dataframe. Adding empty column for host_disease.")
            df['host_disease'] = None
        
        try:
            logging.debug(f"Lowering the 'isolation_source' data for uniformity in isolation source.")
            df['isolation_source'] = df['isolation_source'].apply(lambda x: str(x).lower() if pd.notnull(x) else None)
        except:
            logging.debug(f"isolation_source column not found in the dataframe. Adding empty column for isolation_source.")
            df['isolation_source'] = None


        try:
            logging.info(f"Loading the local database of isolation sources for categorizing the isolation sources of isolates.")
            try:
                logging.debug(f"Trying to load the file with default encoding: utf-8")
                df_isolation_sources_categorized = pd.read_csv(os.path.abspath('./data/categorized_isolation_sources.tsv'), sep='\t', low_memory=False)
                logging.debug(f"Length of the local isolation source categorization database is {len(df_isolation_sources_categorized)}")            
            except UnicodeDecodeError as e:
                logging.error(f"UTF-8 decode failed: {e}")
                logging.info("Retrying with encoding: ISO-8859-1")
                df_isolation_sources_categorized = pd.read_csv(os.path.abspath('./data/categorized_isolation_sources.tsv'), sep='\t', low_memory=False, encoding='ISO-8859-1')
            except FileNotFoundError as e:
                logging.warning("Primary file not found. Trying backup file...")
                logging.info(f"Loading the local database of isolation sources for categorizing the isolation sources of isolates from the backup file.")
                df_isolation_sources_categorized = pd.read_csv(os.path.abspath('./data/categorized_isolation_sources_backup.tsv'), sep='\t', low_memory=False)
            except Exception as e:
                logging.error(f"All attempts failed while loading local database: {e}")
                logging.info(f"Exiting the isolation source categorization process.")
                messagebox.showerror('Error!', f"An error occurred while loading the local database of isolation sources for categorizing the isolation sources of isolates. The error is: {e}")
                raise e
        except Exception as e:
            logging.error(f"An error occurred while loading the local database of isolation sources for categorizing the isolation sources of isolates. The error is: {e}")
            df_isolation_sources_categorized = None
            # messagebox.showerror('Error!', f"An error occurred while loading the local database of isolation sources for categorizing the isolation sources of isolates. The error is: {e}")
            # raise e
            
        if df_isolation_sources_categorized is not None and not df_isolation_sources_categorized.empty:
            logging.info(f"Local database of isolation sources for categorizing the isolation sources of isolates has been loaded successfully.")
            try:
                logging.debug(f"Making `None` in place of `NaN` for the isolation source categorization database.")
                df_isolation_sources_categorized = df_isolation_sources_categorized.where(pd.notnull(df_isolation_sources_categorized),None)
                logging.debug(f"Isolation source categorization database has been made to include `None` in place of `NaN`.")
                
                fuzzy_to_isolation_source_map = {tuple(row[['host', 'host_disease', 'isolation_source']].values): row.to_dict() for _, row in df_isolation_sources_categorized.iterrows()}
                logging.debug(f"Local database of categories of isolation source has been loaded into a functional map. The length of the database now is {len(fuzzy_to_isolation_source_map)}.")
            except Exception as e:
                logging.error(f"An error occurred while creating the fuzzy to isolation source map. The error is: {e}")
                messagebox.showerror('Error!', f"An error occurred while creating the fuzzy to isolation source map. The error is: {e}")
                raise e
        else:
            logging.error(f"Local database for isolation source categorization was not found or is empty.")
            messagebox.showerror('Error!', "Local database for isolation source categorization was not found or is empty.")
            raise FileNotFoundError("Local database for isolation source categorization was not found or is empty.")

        # else:               
        #     logging.debug(f"local database for isolation source categorization was not found. Starting a new dictionary")
        #     fuzzy_to_isolation_source_map = {}

        updated_df_isolation_sources_categorized = dict_update_from_new_isolation_sources(df, fuzzy_to_isolation_source_map)
        logging.debug(f"Size of the local database for isolation source categorization after update is: {len(updated_df_isolation_sources_categorized)}")

        logging.debug(f"Now, categorizing the isolation sources in df via mapping!")
        df = map_isolation_sources(df, map_df=updated_df_isolation_sources_categorized)

        logging.info(f"Now, categorizing the annotations.")
        df['Annotation_category'] = df['Annotation_from'].apply(categorizing_annotation_for_graph)
        logging.info(f"Annotation categories have been added to the dataframe.")
        
        logging.info(f"Now, categorizing the sequencing technologies used for sequencing.")
        df['Categorized_sequencing_technologies'] = df['Sequencing_Technology'].apply(lambda x: categorize_sequencing_technology(str(x)) if x != None else 'Unknown')
        logging.info(f"Sequencing technologies are categorized and added to the dataframe.")

        df = extract_year_from_dandt_stamp(df, 'Submission_date', 'Submission_year')
        logging.info(f"Assembly submission date and time stamp has been converted to years and added to dataframe.")

        try:
            logging.debug(f"Removing the `X` from the coverage column for plotting.")
            df['Coverage_Depth'] = df['Coverage_Depth'].str.replace('x', '')
            df['Coverage_Depth'] = pd.to_numeric(df['Coverage_Depth'], errors='coerce')
            logging.info(f"Sequencing Coverage data has been cleaned and binned for plotting.")
        except:
            pass

        closing_function(country_state_map, updated_df_isolation_sources_categorized)
        logging.debug(f"Local database of country names and state names has been updated and saved.")

        date_time = datetime.now()
        date_time = str(date_time).replace(" ", "_")
        date_time = date_time.replace(":", "-")
        
        
        saving_file_name = "Allmetadata_" + date_time + ".tsv"
        saving_file_as = str(os.path.join(saving_file_path, saving_file_name)) 

        df.to_csv(saving_file_as, index=False, sep='\t')
        
        logging.info(f"Metadata file saved successfully at {os.path.abspath(saving_file_as)}")
        
        return df
    
    except Exception as e:
        logging.error(f"An error occurred while normalizing the data. {e}")
        # messagebox.showerror('Error!', e)
        return None

def meta_mined_app(df: pd.DataFrame, saving_file_path:str = None):
    logging.info(f"Dash started!")

    if saving_file_path is None:
        saving_file_path = os.getcwd()
        logging.debug(f"Since the `saving file path` is `None`, the current working directory is choosen. The WD is {saving_file_path}")
    else:
        logging.debug(f"the current saving path: {saving_file_path}")

    # for country drop-down
    df_for_cdo = df[pd.notnull(df['country_common_name'])]

    country_dropdown_options = df_for_cdo['country_common_name'].dropna().astype(str).unique()

    country_dropdown_options.sort()
    country_dropdown_options = np.insert(country_dropdown_options, 0, 'World')

    # for sequencing technology drop-down
    sequencing_technology_options = df['Categorized_sequencing_technologies'].dropna().astype(str).unique()
    sequencing_technology_options.sort()
    sequencing_technology_options = np.insert(sequencing_technology_options, 0, 'All')

    # rename the column if source_y is present
    if 'source_y' in df.columns:
        logging.debug(f"`source_y` column is present in the dataframe. Renaming it to `source`")
        df.rename(columns={'source_y': 'source'}, inplace=True)
    elif 'source' in df.columns:
        logging.debug(f"`source` column only is there in the dataframe.")

    if 'sample_y' in df.columns:
        logging.debug(f"`sample_y` column is present in the dataframe. Renaming it to `sample`")
        df.rename(columns={'sample_y': 'sample'}, inplace=True)
    elif 'sample' in df.columns:
        logging.debug(f"`sample` column only is there in the dataframe.")

    # for sliders making sure the columns are numeric
    df['Coverage_Depth'] = pd.to_numeric(df['Coverage_Depth'], errors='coerce') # making sure this column is numeric

    logging.debug(f"Coverage Depth column has been converted to numeric.")


    df['Total_genes'] = pd.to_numeric(df['Total_genes'], errors='coerce')
    logging.debug(f"Total genes column has been converted to numeric.")

    df['Protein-coding_genes'] = pd.to_numeric(df['Protein-coding_genes'], errors='coerce')
    logging.debug(f"Protein coding genes column has been converted to numeric.")

    df['Non-coding_genes'] = pd.to_numeric(df['Non-coding_genes'], errors='coerce')
    logging.debug(f"Non-coding genes column has been converted to numeric.")

    df['Protein-coding_genes'] = pd.to_numeric(df['Protein-coding_genes'], errors='coerce')
    logging.debug(f"Protein coding genes column has been converted to numeric.")

    df['Pseudogenes'] = pd.to_numeric(df['Pseudogenes'], errors='coerce')
    logging.debug(f"Pseudogenes column has been converted to numeric.")
    
    df['ANI_best_match_score'] = pd.to_numeric(df['ANI_best_match_score'], errors='coerce')
    logging.debug(f"ANI best match score column has been converted to numeric.")

    df['ANI_best_matched_assembly\'s_coverage'] = pd.to_numeric(df['ANI_best_matched_assembly\'s_coverage'], errors='coerce')
    logging.debug(f"ANI best matched assembly's coverage column has been converted to numeric.")

    df['Contig_N50'] = pd.to_numeric(df['Contig_N50'], errors='coerce')
    logging.debug(f"Contig N50 column has been converted to numeric.")

    df['Contig_L50'] = pd.to_numeric(df['Contig_L50'], errors='coerce')
    logging.debug(f"Contig L50 column has been converted to numeric.")

    flask_logger = logging.getLogger('werkzeug')
    flask_logger.setLevel(logging.ERROR)
    
    # print("Initializing app...")
    meta_mined = Dash(__name__, external_stylesheets=[dbc.themes.LUX])

    meta_mined.layout = dcc.Loading(
        id = 'loading-fullpage',
        type='circle',
        fullscreen=True,
        overlay_style={"visibility": "visible"},
        delay_hide=800,
        delay_show=800,
        children=[
            html.Div(
                children=[
                    html.H1(
                        'Meta-mined!',
                        style={'textAlign': 'center', 'fontSize': '2.5em', 'color': '#333', 'margin': '20px 0'}
                    ),
                    html.Hr(style={'border': '1px solid #ccc', 'margin': '10px 0'}),
                    
                    # LED Display for genome count
                    html.Div(
                        style={'display': 'flex', 'justifyContent': 'center', 'marginBottom': '20px'},
                        children=[
                            daq.LEDDisplay(
                                id='genome-count',
                                value=len(df),
                                color="#5e5eff",
                                backgroundColor="#f9f9f9",
                                size=30,
                            ),
                        ],
                    ),
                    
                    # First row: Radio buttons for atypical and suppressed assembly options
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 1, 'padding': '10px', 'backgroundColor': '#f9f9f9', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    dcc.RadioItems(
                                        id='atypical-radio',
                                        options=[
                                            {'label': 'Include Atypical Assemblies', 'value': 'all'},
                                            {'label': 'Exclude Atypical Assemblies', 'value': 'no_atypical'},
                                            {'label': 'Only Atypical Assemblies', 'value': 'only_atypical'},
                                        ],
                                        value='all',
                                        labelStyle={'display': 'block', 'margin': '5px 0'},
                                    ),
                                    html.Div(
                                        id='atypical-radio-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555'},
                                    ),
                                ],
                            ),
                            html.Div(
                                style={
                                    'flex': 1, 'padding': '10px', 'backgroundColor': '#f9f9f9', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    dcc.RadioItems(
                                        id='suppressed-radio',
                                        options=[
                                            {'label': 'Include Suppressed Assemblies', 'value': 'all'},
                                            {'label': 'Exclude Suppressed Assemblies', 'value': 'no_suppressed'},
                                            {'label': 'Only Suppressed Assemblies', 'value': 'only_suppressed'},
                                        ],
                                        value='all',
                                        labelStyle={'display': 'block', 'margin': '5px 0'},
                                    ),
                                    html.Div(
                                        id='suppressed-radio-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555'},
                                    ),
                                ],
                            ),
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '10px', 'backgroundColor': '#f9f9f9', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Assembly Submission Year', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='submission-year-slider',
                                        min=1980,
                                        max=2025,
                                        step=1,
                                        marks={i: str(i) for i in range(1980, 2026, 5)},
                                        value=[1980, 2025],
                                    ),
                                    html.Div(
                                        id='submission-year-slider-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                ],
                            ),
                        ],
                    ),
                    
                    # Second row: Choropleth map for isolate's geographical location
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 3, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Geographical Distribution of Isolates', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Dropdown(
                                        id='country-dropdown',
                                        options=country_dropdown_options,
                                        value='World',
                                        placeholder='Select a country',
                                        style={'marginBottom': '5px', 'justifyContent': 'center'}
                                    ),
                                    html.Div(
                                        style={'display': 'flex', 'justifyContent': 'center'},
                                        children=[
                                            dcc.Graph(
                                                id='choropleth-map',
                                                figure=Choropleth_map(df=df, selected_country='World'),
                                            )
                                        ]
                                    ),
                                ],
                            ),
                        ],
                    ),

                    # Third row: Graph for assembly level and annotations
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Assembly Level:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Checklist(
                                        id='assembly-level-checklist',
                                        options=[
                                            {'label': 'Complete Genome', 'value': 'Complete Genome'},
                                            {'label': 'Chromosome', 'value': 'Chromosome'},
                                            {'label': 'Scaffold', 'value': 'Scaffold'},
                                            {'label': 'Contig', 'value': 'Contig'},
                                        ],
                                        value=['Complete Genome', 'Chromosome', 'Scaffold', 'Contig'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                    ),
                                    html.Div(
                                        style={'display': 'flex', 'justifyContent': 'center'},
                                        children=[
                                            dcc.Graph(
                                                id='assembly-level-graph',
                                                figure=Assembly_level_bar(
                                                    df=df,
                                                    selected_assembly_levels=['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']
                                                ),
                                            ),
                                        ]
                                    ),
                                ],
                            ),
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Annotation From:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Checklist(
                                        id='annotation-checklist',
                                        options=[
                                            {'label': 'GenBank', 'value': 'GenBank'},
                                            {'label': 'NCBI RefSeq', 'value': 'NCBI RefSeq'},
                                            {'label': 'Others', 'value': 'Others'},
                                            {'label': 'No Annotation', 'value': 'No Annotation'},
                                        ],
                                        value=['GenBank', 'NCBI RefSeq', 'Others', 'No Annotation'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                    ),
                                    html.Div(
                                        style={'display': 'flex', 'justifyContent': 'center'},
                                        children=[
                                            dcc.Graph(
                                                id='annotation-graph',
                                                figure=Annotation_bar(
                                                    df=df,
                                                    show_annotations_from=['GenBank', 'NCBI RefSeq', 'Others', 'No Annotation']
                                                ),
                                            ),
                                        ],
                                    ),
                                ],
                            ),

                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Submissions over the years:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Graph(
                                        id='submission-year-scatter',
                                        figure=Submission_year_line(df=df),
                                    ),
                                ],
                            ),
                        ],
                    ),
                    
                    # Fourth row: Sequencing technologies and coverage depth
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Sequencing Technologies', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Dropdown(
                                        id='sequencing-technology-dropdown',
                                        options=[{'label': i, 'value': i} for i in sequencing_technology_options],
                                        multi=True,
                                        value=['All'],
                                        placeholder='Select a sequencing technology',
                                        style={'marginBottom': '5px'}
                                    ),
                                    dcc.Graph(
                                        id='sequencing-technologies-line',
                                        figure=Sequencing_technologies_scatter(df=df, selected_sequencing_technologies=['All']),
                                    ),
                                ],
                            ),
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Coverage Depth', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='coverage-slider',
                                        min=0,
                                        max=5000,
                                        step=25,
                                        marks={i: str(i) for i in range(0, 5001, 500)},
                                        value=[0, 5000],
                                    ),
                                    dcc.Checklist(
                                        id='coverage-checklist',
                                        options=[
                                            {'label': 'Include coverage depth of >5000', 'value': '> 5000'},
                                        ],
                                        value=['> 5000'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='coverage-slider-output1',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='coverage-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no coverage data', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='coverage-slider-output2',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    html.Div(
                                        style={'display': 'flex', 'justifyContent': 'center'},
                                        children=[
                                            dcc.Graph(
                                                id='coverage-depth-bar',
                                                figure=Coverage_bar(df=df),
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),

                    # Fifth row: Average Nucleotide Identity and Contig L50 - N50
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Average Nucleotide Identity', style={'marginBottom': '10px', 'color': '#333'}),
                                    html.H6('% Identity', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='ani-identity-slider',
                                        min=0,
                                        max=100,
                                        step=1,
                                        marks={i: str(i) for i in range(0, 101, 10)},
                                        value=[0, 100],
                                    ),
                                    dcc.Checklist(
                                        id='ani-identity-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no ANI % Identity data', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='ani-identity-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    html.H6('% Coverage', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='ani-coverage-slider',
                                        min=0,
                                        max=100,
                                        step=1,
                                        marks={i: str(i) for i in range(0, 101, 10)},
                                        value=[0, 100],
                                    ),
                                    dcc.Checklist(
                                        id='ani-coverage-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no ANI %coverage data', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='ani-coverage-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    html.Div(
                                        id='ani-sliders-outputs',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='ani-scatter',
                                        figure=ANI_scatter(df=df),
                                    ),
                                ],
                            ),

                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Contig L50 - N50:', style={'marginBottom': '10px', 'color': '#333'}),
                                    html.H6('Contig N50', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='n50-slider',
                                        min=int(df['Contig_N50'].min()),
                                        max=int(df['Contig_N50'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Contig_N50'].max()), 1000000)},
                                        value=[0, int(df['Contig_N50'].max())],
                                    ),
                                    html.Br(),
                                    html.H6('Contig L50', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='l50-slider',
                                        min=int(df['Contig_L50'].min()),
                                        max=int(df['Contig_L50'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Contig_L50'].max()), int(df['Contig_L50'].max()/10))},
                                        value=[0, int(df['Contig_L50'].max())],
                                    ),
                                    html.Div(
                                        id='n50l50-slider-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    html.Br(),
                                    dcc.Graph(
                                        id='n50l50-scatter',
                                        figure=N50L50_scatter(df=df),
                                    ),
                                ],
                            ),
                        ],
                    ),

                    # Sixth row: histograms for total genes, CDSs, non-coding genes, and pseudogenes
                    html.Div(
                        style={
                            'display': 'grid',
                            'gridTemplateColumns': 'repeat(4, 1fr)',  
                            'gap': '20px',  
                            'marginBottom': '20px',
                            # # 'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                            # # 'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                            'padding': '20px',
                            'backgroundColor': '#f8f9fa',
                            'borderRadius': '10px'
                        },
                        children=[
                            
                            html.H5('Annotations:', style={'gridColumn': 'span 4', 'marginBottom': '10px', 'color': '#333', 'marginLeft':'20px'}),
                            # html.Br(),

                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    dcc.RangeSlider(
                                        id = 'total-gene-slider',
                                        min = int(df['Total_genes'].min()),
                                        max = int(df['Total_genes'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Total_genes'].max()), 2000)},
                                        value=[0, int(df['Total_genes'].max())],
                                    ),
                                    html.Div(
                                        id='total-gene-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='total-gene-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no total gene count', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='total-gene-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='total-gene-hist',
                                        figure=Total_genes_hist(df=df),
                                    ),
                                ]
                            ),

                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    dcc.RangeSlider(
                                        id='cds-slider',
                                        min=int(df['Protein-coding_genes'].min()),
                                        max=int(df['Protein-coding_genes'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Protein-coding_genes'].max()), 2000)},
                                        value=[0, int(df['Protein-coding_genes'].max())],
                                    ),
                                    html.Div(
                                        id='cds-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='cds-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no CDSs count', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='cds-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='cds-hist',
                                        figure=CDSs_hist(df=df),
                                    ),
                                ]
                            ),

                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    dcc.RangeSlider(
                                        id='non-coding-slider',
                                        min=int(df['Non-coding_genes'].min()),
                                        max=int(df['Non-coding_genes'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Non-coding_genes'].max()), 50)},
                                        value=[0, int(df['Non-coding_genes'].max())],
                                    ),
                                    html.Div(
                                        id='non-coding-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='non-coding-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no non-coding gene count', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='non-coding-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='non-coding-hist',
                                        figure=Non_coding_hist(df=df),
                                    ),
                                ]
                            ),

                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    dcc.RangeSlider(
                                        id='pseudogene-slider',
                                        min=int(df['Pseudogenes'].min()),
                                        max=int(df['Pseudogenes'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Pseudogenes'].max()), 500)},
                                        value=[0, int(df['Pseudogenes'].max())],
                                    ),
                                    html.Div(
                                        id='pseudogene-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='pseudogene-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no Pseudogene count', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='pseudogene-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='pseudogene-hist',
                                        figure=Pseudogenes_hist(df=df),
                                    ),
                                ]
                            ),
                        ]
                    ),


                    html.Div(
                        style={
                            'display': 'block',
                            'gap': '20px',
                            'marginBottom': '20px',
                            'padding': '20px',
                        },
                        children=[
                            html.H5('Filter BioProject:', style={'marginBottom': '10px', 'color': '#333',}),
                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    html.H6('Search BioProject by keywords:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Input(
                                        id='bioproject-input',
                                        type='text',
                                        placeholder='Search Keywords separated by commas...',
                                        debounce=True,
                                        style={
                                            'width': '700px',
                                            'height': '40px',         
                                            'fontSize': '16px',       
                                            'padding': '5px 10px',    
                                            'borderRadius': '5px',
                                            'marginBottom': '15px',     
                                        }
                                    ),
                                    html.H6('Selected BioProjects:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Dropdown(
                                        id='bioproject-dropdown',
                                        options=[{'label': d, 'value': d} for d in df['Bioproject_title'].unique()],
                                        placeholder='Filtered BioProjects will appear here...',
                                        searchable=True,
                                        multi=True,
                                    ),
                                    html.Div(
                                        id='bioproject-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                ],
                            ),
                        ],
                    ),

                    html.Div(
                        style={
                            'display': 'block',
                            'gap': '20px',
                            'marginBottom': '20px',
                            'padding': '20px',
                        },
                        children=[
                            html.H5('Filter BioSample:', style={'marginBottom': '10px', 'color': '#333', }),
                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    html.H6('Search BioSample by keywords:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Input(
                                        id='biosample-input',
                                        type='text',
                                        placeholder='Search Keywords separated by commas...',
                                        debounce=True,
                                        style={
                                            'width': '700px',
                                            'height': '40px',         
                                            'fontSize': '16px',       
                                            'padding': '5px 10px',    
                                            'borderRadius': '5px',
                                            'marginBottom': '15px'    
                                        }
                                    ),
                                    html.H6('Selected BioSamples:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Dropdown(
                                        id='biosample-dropdown',
                                        options=[{'label': d, 'value': d} for d in df['Biosample_title'].unique()],
                                        placeholder='Filtered BioSamples will appear here...',
                                        searchable=True,
                                        multi=True,
                                    ),
                                    html.Div(
                                        id='biosample-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                ],
                            ),
                        ],
                    ),

                    # Eighth row: Isolation source treemap and dropdowns
                    html.Div(
                        style={
                            'display': 'flex',
                            'gap': '20px',
                            'marginBottom': '20px',
                            'padding': '20px',
                            'backgroundColor': '#f8f9fa',
                            'borderRadius': '10px'
                        },
                        children=[
                            # Left Panel - Dropdowns in Two Sections
                            html.Div(
                                style={
                                    'flex': 2,
                                    'padding': '20px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '10px',
                                    'boxShadow': '0 4px 8px rgba(0, 0, 0, 0.1)',
                                    'display': 'flex',
                                    'flexDirection': 'column',
                                    'gap': '15px'
                                },
                                children=[
                                    html.H5('Isolation Source', style={'marginBottom': '10px', 'color': '#333', 'fontWeight': 'bold'}),

                                    # All Dropdowns in One Column First
                                    html.Div(
                                        style={'display': 'flex', 'flexDirection': 'column', 'gap': '15px'},
                                        children=[
                                            html.Div(
                                                children=[
                                                    html.H6('Select Host', style={'marginBottom': '8px', 'color': '#444'}),
                                                    dcc.Dropdown(
                                                        id='identified-host-dropdown',
                                                        options=[
                                                            {'label': 'Hospital-associated', 'value': 'Hospital-associated'},
                                                            {'label': 'Animal-associated', 'value': 'Animal-associated'},
                                                            {'label': 'Environment-associated', 'value': 'Environment-associated'},
                                                            {'label': 'Laboratory-based', 'value': 'Laboratory-based'},
                                                            {'label': 'Unknown', 'value': 'Unknown'}
                                                        ],
                                                        multi=True,
                                                        value=['Hospital-associated', 'Animal-associated', 'Environment-associated', 'Laboratory-based','Unknown'],
                                                        placeholder="Select Host",
                                                    )
                                                ]
                                            ),
                                            html.Div(
                                                children=[
                                                    html.H6('Select Source/Disease Category', style={'marginBottom': '8px', 'color': '#444'}),
                                                    dcc.Dropdown(
                                                        id='source-category-dropdown',
                                                        multi=True,
                                                        placeholder="Select Source/Disease Category",
                                                    )
                                                ]
                                            ),
                                            html.Div(
                                                children=[
                                                    html.H6('Select Source/Disease', style={'marginBottom': '8px', 'color': '#444'}),
                                                    dcc.Dropdown(
                                                        id='source-dropdown',
                                                        multi=True,
                                                        placeholder="Select Source/Disease",
                                                    )
                                                ]
                                            ),
                                            html.Div(
                                                children=[
                                                    html.H6('Select Sample', style={'marginBottom': '8px', 'color': '#444'}),
                                                    dcc.Dropdown(
                                                        id='sample-dropdown',
                                                        multi=True,
                                                        placeholder="Select Sample",
                                                    )
                                                ]
                                            ),
                                        ]
                                    ),
                                ],
                            ),
                            # Right Panel - Treemap Output on Top, Graph Below
                            html.Div(
                                style={
                                    'flex': 1,
                                    'display': 'flex',
                                    'flexDirection': 'column',
                                    'alignItems': 'center',
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '10px',
                                    'boxShadow': '0 4px 8px rgba(0, 0, 0, 0.1)',
                                    'gap': '15px'
                                },
                                children=[
                                    # Treemap Output (Text or Additional Info)
                                    html.Div(
                                        id='treemap-output',
                                        style={
                                            'textAlign': 'center',
                                            'fontSize': '0.9em',
                                            'color': '#555',
                                            'padding': '10px',
                                            'backgroundColor': '#f1f1f1',
                                            'borderRadius': '8px',
                                            'width': '100%'
                                        },
                                    ),

                                    # Treemap Graph
                                    dcc.Graph(
                                        id='isolation-source-treemap',
                                        figure=Isolation_source_treemap(df=df),
                                        style={'width': '100%', 'height': '100%'}
                                    ),
                                ]
                            ),
                        ]
                    ),

        
                    html.Div(
                        children=[
                            html.Div(
                                children=[
                                    html.Button(
                                        "Save data",
                                        id='save-button',
                                        n_clicks=0,
                                        style={
                                            'backgroundColor': '#007bff',  
                                            'color': 'white', 
                                            'border': 'none', 
                                            'borderRadius': '5px', 
                                            'padding': '10px 20px', 
                                            'fontSize': '1em', 
                                            'cursor': 'pointer',
                                            'boxShadow': '0px 4px 6px rgba(0, 0, 0, 0.1)'
                                        }
                                    ),
                                ],
                                style={
                                    'display': 'flex',  
                                    'justifyContent': 'center',  
                                    'alignItems': 'center',  
                                    'height': '100px'  
                                }
                            ),
                            html.Div(
                                id='save-button-output',
                                style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'}
                            ),
                        ]
                    ),
                ],
            ),
        ],
    )
    # Define callback
    @meta_mined.callback(
        [Output('choropleth-map', 'figure'),
        Output('assembly-level-graph', 'figure'),
        Output('annotation-graph', 'figure'),
        Output('submission-year-slider-output', 'children'),
        Output('atypical-radio-output', 'children'),
        Output('suppressed-radio-output', 'children'),
        Output('genome-count', 'value'),
        Output('sequencing-technologies-line', 'figure'),
        Output('submission-year-scatter', 'figure'),
        Output('coverage-depth-bar', 'figure'),
        Output('coverage-slider-output1', 'children'),
        Output('coverage-slider-output2', 'children'),
        Output('ani-scatter', 'figure'),
        Output('ani-identity-null-output', 'children'),
        Output('ani-sliders-outputs', 'children'),
        Output('ani-coverage-null-output', 'children'),
        Output('n50l50-slider-output', 'children'),
        Output('n50l50-scatter', 'figure'), 
        [Output('total-gene-output', 'children'),
         Output('total-gene-hist', 'figure')],
        Output('total-gene-null-output', 'children'),
        [Output('cds-output', 'children'),
         Output('cds-hist', 'figure')],
        Output('cds-null-output', 'children'),
        [Output('non-coding-output', 'children'),
         Output('non-coding-hist', 'figure')],
        Output('non-coding-null-output', 'children'),
        [Output('pseudogene-output', 'children'),
         Output('pseudogene-hist', 'figure')],
        Output('pseudogene-null-output', 'children'),
        [Output('bioproject-dropdown', 'options'),
        Output('bioproject-dropdown', 'value')],
        Output('bioproject-output', 'children'),
        [Output('biosample-dropdown', 'options'),
        Output('biosample-dropdown', 'value')],
        Output('biosample-output', 'children'),
        Output('source-category-dropdown', 'options'),
        Output('source-dropdown', 'options'),
        Output('sample-dropdown', 'options'),
        [Output('treemap-output', 'children'),
        Output('isolation-source-treemap', 'figure')],
        Output('save-button-output', 'children'),],
        [Input('country-dropdown', 'value'),
        Input('assembly-level-checklist', 'value'),
        Input('annotation-checklist', 'value'),
        Input('submission-year-slider', 'value'),
        Input('atypical-radio', 'value'),
        Input('suppressed-radio', 'value'),
        Input('genome-count', 'id'),
        Input('sequencing-technology-dropdown', 'value'),
        Input('submission-year-scatter', 'id'),
        Input('coverage-slider', 'value'),
        Input('coverage-checklist', 'value'),
        Input('coverage-include-null-checklist', 'value'),
        Input('ani-identity-slider', 'value'),
        Input('ani-identity-include-null-checklist', 'value'),
        Input('ani-coverage-slider', 'value'),
        Input('ani-coverage-include-null-checklist', 'value'),
        Input('n50-slider', 'value'),
        Input('l50-slider', 'value'),
        Input('total-gene-slider', 'value'),
        Input('total-gene-include-null-checklist', 'value'),
        Input('cds-slider', 'value'),
        Input('cds-include-null-checklist', 'value'),
        Input('non-coding-slider', 'value'),
        Input('non-coding-include-null-checklist', 'value'),
        Input('pseudogene-slider', 'value'),
        Input('pseudogene-include-null-checklist', 'value'),
        Input('bioproject-input', 'value'),
        Input('bioproject-dropdown', 'value'),
        Input('biosample-input', 'value'),
        Input('biosample-dropdown', 'value'),
        Input('identified-host-dropdown', 'value'),
        Input('source-category-dropdown', 'value'),
        Input('source-dropdown', 'value'),
        Input('sample-dropdown', 'value'),
        Input('save-button', 'n_clicks')]
    )

    # print("we are at update_dash")
    def update_dash(selected_country, selected_assembly_level, show_annotations_from, year_range, atypical_radio, suppressed_radio, genome_count_id, selected_sequencing_technologies, submission_year_line_id, coverage_range, coverage_checklist, coverage_include_null_checklist, ani_identity_range, ani_identity_include_null_checklist, ani_coverage_range, ani_coverage_include_null_checklist, n50_range, l50_range,total_genes_range, total_genes_include_null_checklist, cds_range, cds_include_null_checklist, non_coding_range, non_coding_include_null_checklist, pseudogene_range, pseudogene_include_null_checklist, bioproject_input, bioproject_dropdown_values, biosample_input, biosample_dropdown_values, identified_host, source_category, source, sample, n_save_button_clicks):
        modified_df = df.copy()

        logging.debug(f"Shape of the modified_df before filtering: {modified_df.shape}")

        for col in modified_df.columns:
            if modified_df[col].apply(lambda x: isinstance(x, dict)).any():
                logging.debug(f"Converting {col} data to string values as dictionaries may cause problem in filtering/removing duplicates.")
                modified_df[col] = modified_df[col].apply(str)


        if selected_country == 'World':
            modified_df = modified_df
            logging.debug(f"Seleted country is World. Shape of the modified_df after filtering by country: {modified_df.shape}")
        else:
            modified_df = modified_df[modified_df['country_common_name'] == selected_country]
            logging.debug(f"Selected country is {selected_country}. Shape of the modified_df after filtering by country: {modified_df.shape}")

        modified_df = modified_df[modified_df['Assembly_level'].isin(selected_assembly_level)]
        logging.debug(f"Selected assembly levels: {selected_assembly_level}")
        logging.debug(f"Shape of the modified_df after filtering by assembly level: {modified_df.shape}")

        modified_df = modified_df[modified_df['Annotation_category'].isin(show_annotations_from)]
        logging.debug(f"Selected annotation categories: {show_annotations_from}")
        logging.debug(f"Shape of the modified_df after filtering by annotation category: {modified_df.shape}")

        
        modified_df = modified_df[(modified_df['Submission_year'] >= year_range[0]) & (modified_df['Submission_year'] <= year_range[1])]


        if atypical_radio == 'no_atypical':
            earlier_length = modified_df.shape[0]
            modified_df = modified_df[modified_df['Assmbly_atypical?'] == 'No']
            atypical_output = f"Number of atypical assemblies excluded: {earlier_length - modified_df.shape[0]}"

        elif atypical_radio == 'only_atypical':
            earlier_length = modified_df.shape[0]
            modified_df = modified_df[modified_df['Assmbly_atypical?'] == 'Yes']
            atypical_output = f"Number of atypical assemblies: {modified_df.shape[0]}"
        else:
            modified_df = modified_df
            atypical_output = None

        if suppressed_radio == 'no_suppressed':
            all_length = modified_df.shape[0]
            modified_df = modified_df[modified_df['Assembly_status'] == 'current']
            suppressed_output = f"Number of suppressed assemblies excluded: {all_length - modified_df.shape[0]}"
        elif suppressed_radio == 'only_suppressed':
            all_length = modified_df.shape[0]
            modified_df = modified_df[modified_df['Assembly_status'] == 'suppressed']
            suppressed_output = f"Number of suppressed assemblies: {modified_df.shape[0]}"
        else:
            modified_df = modified_df
            suppressed_output = None

        if 'All' in selected_sequencing_technologies:
            modified_df = modified_df
            logging.debug(f"Selected sequencing technology: {selected_sequencing_technologies}")
            logging.debug(f"Shape of the modified_df after filtering by sequencing technology: {modified_df.shape}")
        else:
            modified_df = modified_df[modified_df['Categorized_sequencing_technologies'].isin(selected_sequencing_technologies)]
            logging.debug(f"Selected sequencing technology: {selected_sequencing_technologies}")
            logging.debug(f"Shape of the modified_df after filtering by sequencing technology: {modified_df.shape}")
        
        logging.debug(f"Shape of the modified_df after sorting sequencing technologies: {modified_df.shape}")
        coverage_null_text=""
        if coverage_range[0] == 0 & coverage_range[1] == 5000:            
            if '> 5000' in coverage_checklist and 'all' in coverage_include_null_checklist:
                modified_df = modified_df  
            elif '> 5000' in coverage_checklist and 'all' not in coverage_include_null_checklist:
                coverage_null_text = f"Genomes with no coverage data are excluded."
                coverage_df_above_5000 = modified_df[modified_df['Coverage_Depth'] > 5000]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_df_above_5000, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' not in coverage_checklist and 'all' in coverage_include_null_checklist:
                coverage_null_df = modified_df[pd.isnull(modified_df['Coverage_Depth'])]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_null_df, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' not in coverage_checklist and 'all' not in coverage_include_null_checklist:
                coverage_null_text = f"Genomes with no coverage data are excluded."
                modified_df = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
        else:
            if '> 5000' in coverage_checklist and 'all' in coverage_include_null_checklist:
                coverage_null_df = modified_df[pd.isnull(modified_df['Coverage_Depth'])]
                coverage_df_above_5000 = modified_df[modified_df['Coverage_Depth'] > 5000]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_null_df, coverage_df_above_5000, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' in coverage_checklist and 'all' not in coverage_include_null_checklist:
                coverage_null_text = f"Genomes with no coverage data are excluded."
                coverage_df_above_5000 = modified_df[modified_df['Coverage_Depth'] > 5000]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_df_above_5000, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' not in coverage_checklist and 'all' in coverage_include_null_checklist:
                coverage_null_df = modified_df[pd.isnull(modified_df['Coverage_Depth'])]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_null_df, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' not in coverage_checklist and 'all' not in coverage_include_null_checklist:
                coverage_null_text = f"Genomes with no coverage data are excluded."
                modified_df = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]

        for col in modified_df.columns:
            if modified_df[col].apply(lambda x: isinstance(x, dict)).any():
                logging.debug(f"Column '{col}' contains dictionaries.")

        logging.debug(f"Shape of the modified_df after coverage filters: {modified_df.shape}")
        ani_identity_null_text=""
        if ani_identity_range[0] == 0 & ani_identity_range[1] == 100:
            if 'all' in ani_identity_include_null_checklist:
                modified_df = modified_df
            else:
                ani_identity_null_text = f"{len(modified_df[pd.isnull(modified_df['ANI_best_match_score'])])} Genomes with no % ANI identity data are excluded."
                modified_df = modified_df[pd.notnull(modified_df['ANI_best_match_score'])]
        else:
            if 'all' in ani_identity_include_null_checklist:
                ANI_identity_null_df = modified_df[pd.isnull(modified_df['ANI_best_match_score'])]
                selected_ANI_identity_df = modified_df[(modified_df['ANI_best_match_score'] >= ani_identity_range[0]) & (modified_df['ANI_best_match_score'] <= ani_identity_range[1])]
                modified_df = pd.concat([selected_ANI_identity_df, ANI_identity_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                ani_identity_null_text = f"{len(modified_df[pd.isnull(modified_df['ANI_best_match_score'])])} Genomes with no % ANI identity data are excluded."
                modified_df = modified_df[(modified_df['ANI_best_match_score'] >= ani_identity_range[0]) & (modified_df['ANI_best_match_score'] <= ani_identity_range[1])]
        
        ani_coverage_null_text = ""
        if ani_coverage_range[0] == 0 & ani_coverage_range[1] == 100:
            if 'all' in ani_coverage_include_null_checklist:
                modified_df = modified_df                
            else:
                ani_coverage_null_count = len(modified_df[pd.isnull(modified_df['ANI_best_matched_assembly\'s_coverage'])])
                ani_coverage_null_text = f"{ani_coverage_null_count} genomes with no % coverage data are excluded."
                modified_df = modified_df[pd.notnull(modified_df['ANI_best_matched_assembly\'s_coverage'])]
        else:
            if 'all' in ani_coverage_include_null_checklist:
                ANI_coverage_null_df = modified_df[pd.isnull(modified_df['ANI_best_matched_assembly\'s_coverage'])]
                selected_ANI_coverage_df = modified_df[(modified_df['ANI_best_matched_assembly\'s_coverage'] >= ani_coverage_range[0]) & (modified_df['ANI_best_matched_assembly\'s_coverage'] <= ani_coverage_range[1])]
                modified_df = pd.concat([selected_ANI_coverage_df, ANI_coverage_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                ani_coverage_null_count = len(modified_df[pd.isnull(modified_df['ANI_best_matched_assembly\'s_coverage'])])
                ani_coverage_null_text = f"{ani_coverage_null_count} genomes with no % coverage data are excluded."
                modified_df = modified_df[(modified_df['ANI_best_matched_assembly\'s_coverage'] >= ani_coverage_range[0]) & (modified_df['ANI_best_matched_assembly\'s_coverage'] <= ani_coverage_range[1])]

        if n50_range[0] == 0 & n50_range[1] == int(df['Contig_N50'].max()):
            modified_df = modified_df
        else:
            modified_df = modified_df[(modified_df['Contig_N50'] >= n50_range[0]) & (modified_df['Contig_N50'] <= n50_range[1])]
        
        if l50_range[0] == 0 & l50_range[1] == int(df['Contig_L50'].max()):
            modified_df = modified_df
        else:
            modified_df = modified_df[(modified_df['Contig_L50'] >= l50_range[0]) & (modified_df['Contig_L50'] <= l50_range[1])]

        logging.debug(f"Shape of the modified_df after ani filters: {modified_df.shape}")
        total_gene_null_text = ""
        if total_genes_range[0] == 0 & total_genes_range[1] == int(df['Total_genes'].max()):
            if 'all' in total_genes_include_null_checklist:
                modified_df = modified_df
            else:
                total_gene_null_text = f"{len(modified_df[pd.isnull(modified_df['Total_genes'])])} genomes excluded that didn't have total gene count."
                modified_df = modified_df[pd.notnull(modified_df['Total_genes'])]
        else:
            if 'all' in total_genes_include_null_checklist:
                total_gene_null_df = modified_df[pd.isnull(modified_df['Total_genes'])]
                selected_total_gene_df = modified_df[(modified_df['Total_genes'] >= total_genes_range[0]) & (modified_df['Total_genes'] <= total_genes_range[1])]
                modified_df = pd.concat([selected_total_gene_df, total_gene_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                total_gene_null_text = f"{len(modified_df[pd.isnull(modified_df['Total_genes'])])} genomes excluded that didn't have total gene count."
                modified_df = modified_df[(modified_df['Total_genes'] >= total_genes_range[0]) & (modified_df['Total_genes'] <= total_genes_range[1])]

        cds_null_text = ""
        if cds_range[0] == 0 & cds_range[1] == int(df['Protein-coding_genes'].max()):
            if 'all' in cds_include_null_checklist:
                modified_df = modified_df
            else:
                cds_null_text = f"{len(modified_df[pd.isnull(modified_df['Protein-coding_genes'])])} genomes excluded that didn't have CDSs count."
                modified_df = modified_df[pd.notnull(modified_df['Protein-coding_genes'])]
        else:
            if 'all' in cds_include_null_checklist:
                cds_null_df = modified_df[pd.isnull(modified_df['Protein-coding_genes'])]
                selected_cds_df = modified_df[(modified_df['Protein-coding_genes'] >= cds_range[0]) & (modified_df['Protein-coding_genes'] <= cds_range[1])]
                modified_df = pd.concat([selected_cds_df, cds_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                cds_null_text = f"{len(modified_df[pd.isnull(modified_df['Protein-coding_genes'])])} genomes excluded that didn't have CDSs count."
                modified_df = modified_df[(modified_df['Protein-coding_genes'] >= cds_range[0]) & (modified_df['Protein-coding_genes'] <= cds_range[1])]

        non_coding_null_text = ""
        if non_coding_range[0] == 0 & non_coding_range[1] == int(df['Non-coding_genes'].max()):
            if 'all' in non_coding_include_null_checklist:
                modified_df = modified_df
            else:
                non_coding_null_text = f"{len(modified_df[pd.isnull(modified_df['Non-coding_genes'])])} genomes excluded that didn't have non-coding gene count."
                modified_df = modified_df[pd.notnull(modified_df['Non-coding_genes'])]
        else:
            if 'all' in non_coding_include_null_checklist:
                non_coding_null_df = modified_df[pd.isnull(modified_df['Non-coding_genes'])]
                selected_non_coding_df = modified_df[(modified_df['Non-coding_genes'] >= non_coding_range[0]) & (modified_df['Non-coding_genes'] <= non_coding_range[1])]
                modified_df = pd.concat([selected_non_coding_df, non_coding_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                non_coding_null_text = f"{len(modified_df[pd.isnull(modified_df['Non-coding_genes'])])} genomes excluded that didn't have non-coding gene count."
                modified_df = modified_df[(modified_df['Non-coding_genes'] >= non_coding_range[0]) & (modified_df['Non-coding_genes'] <= non_coding_range[1])]
        
        pseudogene_null_text = ""
        if pseudogene_range[0] == 0 & pseudogene_range[1] == int(df['Pseudogenes'].max()):
            if 'all' in pseudogene_include_null_checklist:
                modified_df = modified_df
            else:
                pseudogene_null_text = f"{len(modified_df[pd.isnull(modified_df['Pseudogenes'])])} genomes excluded that didn't have pseudogene count"
                modified_df = modified_df[pd.notnull(modified_df['Pseudogenes'])]
        else:
            if 'all' in pseudogene_include_null_checklist:
                pseudogene_null_df = modified_df[pd.isnull(modified_df['Pseudogenes'])]
                selected_pseudogene_df = modified_df[(modified_df['Pseudogenes'] >= pseudogene_range[0]) & (modified_df['Pseudogenes'] <= pseudogene_range[1])]
                modified_df = pd.concat([selected_pseudogene_df, pseudogene_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                pseudogene_null_text = f"{len(modified_df[pd.isnull(modified_df['Pseudogenes'])])} genomes excluded that didn't have pseudogene count"
                modified_df = modified_df[(modified_df['Pseudogenes'] >= pseudogene_range[0]) & (modified_df['Pseudogenes'] <= pseudogene_range[1])]
        
        logging.debug(f"Shape of the modified_df after gene count filters: {modified_df.shape}")

        # filtering the dataframe based on keywords in bioproject and biosample names  
        
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None

        logging.debug(f"Triggered ID: {triggered_id}")
        
        bioproject_options = [{'label': i, 'value': i} for i in modified_df['Bioproject_title'].unique()]

        if bioproject_input:

            bioproject_list = modified_df['Bioproject_title'].tolist()

            keywords = [word.strip().lower() for word in bioproject_input.split(',') if word.strip()]
            
            logging.debug(f"Input for bioproject: {bioproject_input}")

            list_based_on_bioproject_input = [title for title in bioproject_list if any(kw in title.lower() for kw in keywords)]

        else:
            bioproject_dropdown = []
            list_based_on_bioproject_input = []
            logging.debug(f"No filter was used for Bioproject")

        if triggered_id == 'bioproject-input':
            bioproject_dropdown = list(set(bioproject_dropdown_values or []) | set(list_based_on_bioproject_input))
            logging.debug("User typed in input box — adding matches to current dropdown selection")
        elif triggered_id == 'bioproject-dropdown':
            bioproject_dropdown = list(set(bioproject_dropdown_values)) or []
            logging.debug("User changed dropdown — using dropdown values only")
        else:
            if bioproject_dropdown_values:
                # bioproject_dropdown = list(set(bioproject_dropdown_values)) or []
                bioproject_dropdown = [title for title in modified_df['Bioproject_title'].unique() if title in bioproject_dropdown_values]
                logging.debug("Fallback: using whatever is available in dropdown")

        if bioproject_dropdown != []:
            modified_df = modified_df[modified_df['Bioproject_title'].isin(bioproject_dropdown)]
            logging.debug(f"Filtered {len(modified_df)} rows using: {triggered_id}")
        else:
            modified_df = modified_df
            logging.debug("No filtering applied for bioproject")
        

        biosample_options = [{'label': i, 'value': i} for i in modified_df['Biosample_title'].unique()]

        if biosample_input:
            biosample_list = modified_df['Biosample_title'].tolist()

            keywords_2 = [word.strip().lower() for word in biosample_input.split(',') if word.strip()]
            
            logging.debug(f"Input for biosample: {biosample_input}")

            list_based_on_biosample_input = [title for title in biosample_list if any(kw in title.lower() for kw in keywords_2)]
        
        else:
            biosample_dropdown = []
            list_based_on_biosample_input = []
            logging.debug(f"No filter was used for Biosample")
        
        if triggered_id == 'biosample-input':
            biosample_dropdown = list(set(biosample_dropdown_values or []) | set(list_based_on_biosample_input))
            logging.debug("User typed in input box — adding matches to current dropdown selection")
        elif triggered_id == 'biosample-dropdown':
            biosample_dropdown = list(set(biosample_dropdown_values)) or []
            logging.debug("User changed dropdown — using dropdown values only")
        else:
            if biosample_dropdown_values:
                # biosample_dropdown = list(set(biosample_dropdown_values)) or []
                biosample_dropdown = [title for title in modified_df['Biosample_title'].unique() if title in biosample_dropdown_values]
                logging.debug("Fallback: using whatever is available in dropdown")
        
        if biosample_dropdown != []:
            modified_df = modified_df[modified_df['Biosample_title'].isin(biosample_dropdown)]
            logging.debug(f"Filtered {len(modified_df)} rows using: {triggered_id}")

        else:
            modified_df = modified_df
            logging.debug("No filtering applied for biosample")

        modified_df = modified_df[modified_df['identified_host'].isin(identified_host)]
        logging.debug(f"Selected hosts: {identified_host}")
        logging.debug(f"Shape of the modified_df after filtering by host: {modified_df.shape}")

        source_category_options = [{'label': i, 'value': i} for i in modified_df['source_category'].unique()]

        if source_category:
            modified_df = modified_df[modified_df['source_category'].isin(source_category)]
            logging.debug(f"Selected source categories: {source_category}")
            logging.debug(f"Shape of the modified_df after filtering by source category: {modified_df.shape}")
            source_options = [{'label': i, 'value': i} for i in modified_df['source'].unique()]
        else:
            modified_df = modified_df
            source_options = []
        
        if source:
            modified_df = modified_df[modified_df['source'].isin(source)]
            logging.debug(f"Selected sources: {source}")
            logging.debug(f"Shape of the modified_df after filtering by source: {modified_df.shape}")

            sample_options = [{'label': i, 'value': i} for i in modified_df['sample'].unique()]
        else:
            modified_df = modified_df
            sample_options = []

        if sample:
            modified_df = modified_df[modified_df['sample'].isin(sample)]
            logging.debug(f"Selected samples: {sample}")
            logging.debug(f"Shape of the modified_df after filtering by sample: {modified_df.shape}")
        else:
            modified_df = modified_df

        genome_count = modified_df.shape[0]

        logging.debug(f"The genome count is {genome_count}")

        location_map = Choropleth_map(df=modified_df, selected_country=selected_country)
        
        assembly_level_graph = Assembly_level_bar(df=modified_df, selected_assembly_levels=selected_assembly_level)
        
        annotation_graph = Annotation_bar(df=modified_df, show_annotations_from=show_annotations_from)

        submission_year_text = f"Assembly Submission Year Range: {year_range[0]} - {year_range[1]}"

        sequencing_technologies_over_the_years = Sequencing_technologies_scatter(df=modified_df, selected_sequencing_technologies=selected_sequencing_technologies)

        submissions_over_the_years = Submission_year_line(df=modified_df)

        coverage_depth_bar = Coverage_bar(df=modified_df)

        coverage_depth_text = f"Coverage Depth range: {modified_df['Coverage_Depth'].min()}X to {modified_df['Coverage_Depth'].max()}X"
        
        ani_score_scatter = ANI_scatter(df=modified_df)
        
        ani_score_text = f"Selected % Identity: {ani_identity_range[0]}  to {ani_identity_range[1]}, Selected % Coverage: {ani_coverage_range[0]} to {ani_coverage_range[1]}" 

        n50l50_output_text= f"Contig N50: {n50_range[0]} - {n50_range[1]}, Contig L50: {l50_range[0]} - {l50_range[1]}"

        n50l50_scatter = N50L50_scatter(df=modified_df)
        
        total_gene_dist = Total_genes_hist(df=modified_df)
        total_gene_text = f"Selected 'Total gene' range: {total_genes_range[0]} to {total_genes_range[1]}"

        cds_dist = CDSs_hist(df=modified_df)
        cds_text = f"Selected 'CDSs' range: {cds_range[0]} to {cds_range[1]}"
        
        non_coding_dist = Non_coding_hist(df=modified_df)
        non_coding_text = f"Selected 'Non-coding gene' range: {non_coding_range[0]} to {non_coding_range[1]}"

        pseudogene_dist = Pseudogenes_hist(df=modified_df)
        pseudogene_text = f"Selected 'Pseudogene' range: {pseudogene_range[0]} to {pseudogene_range[1]}"

        bioproject_output = f"{len(modified_df['Bioproject_title'].unique())} BioProjects are selected. They contain a total of {len(modified_df)}assemblies."
        
        biosample_output = f"Selected BioSamples contain a total of {len(modified_df)} assemblies."
        
        treemap_df = modified_df.copy()

        treemap_df = treemap_df[treemap_df['identified_host'] != 'Unknown']

        treemap_text = f"Isolation source is Unknown for {len(modified_df[modified_df['identified_host'] == 'Unknown'])} genomes among all selected genomes."
        isolation_treemap = Isolation_source_treemap(df=treemap_df)

        if n_save_button_clicks > 0:
            logging.debug("Button Pressed!")

            try:
                date_time = datetime.now()
                date_time = str(date_time).replace(" ", "_")
                date_time = date_time.replace(":", "-")
                # print(date_time)
                saving_file_name = "filtered_data_" + date_time + ".tsv"
                # print(saving_file_name)
                
                saving_file_as = str(os.path.join(saving_file_path, saving_file_name)) 
                # print(saving_file_as)
                modified_df.to_csv(saving_file_as, index=False, sep='\t')
                logging.info(f"File saved successfully at {os.path.abspath(saving_file_as)}")
                save_message = f"File saved successfully at {os.path.abspath(saving_file_as)}"
            except Exception as e:
                logging.info(f"Error saving file: {str(e)}")
                save_message =  f"Error saving file: {str(e)}"
            
        else:
            save_message = ""



        return location_map, assembly_level_graph, annotation_graph, submission_year_text, atypical_output, suppressed_output, genome_count, sequencing_technologies_over_the_years, submissions_over_the_years, coverage_depth_bar, coverage_depth_text, coverage_null_text, ani_score_scatter, ani_identity_null_text, ani_score_text, ani_coverage_null_text, n50l50_output_text, n50l50_scatter, [total_gene_text, total_gene_dist], total_gene_null_text, [cds_text, cds_dist], cds_null_text, [non_coding_text, non_coding_dist], non_coding_null_text, [pseudogene_text, pseudogene_dist], pseudogene_null_text, [bioproject_options, bioproject_dropdown], bioproject_output, [biosample_options, biosample_dropdown], biosample_output, source_category_options, source_options, sample_options, [treemap_text, isolation_treemap], save_message

    
        # webbrowser.open("http://localhost:8050")
    
    
    # meta_mined.run_server(debug=True, port=8050)

    webbrowser.open("http://localhost:8050")
    logging.debug("Opening the web browser...")
    
    # meta_mined.run_server(debug=True, port=8050)
    logging.debug("Running the server...")
    meta_mined.run(debug=True, port=8050, use_reloader=False, threaded=True) 

# jay
