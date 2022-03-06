#! /usr/bin/env python3

import sys
import kirjava
import atomium
from tqdm import tqdm
from common import structure_site_to_sample, sequence_site_to_sample

import pandas as pd

from utilities import *

API_URL = "https://api.zincbind.net/"

FAMILY_SITES_QUERY = """query familySites($family: String) {
    zincsites(family: $family, pdb__resolution__lt: 2) { edges { node { 
        id pdb { id assembly } chainInteractions { edges { node { sequence } } }
        residues(primary: true) { edges { node {
            id chainSignature atomiumId atoms(name: "CA") {
                edges { node { x y z } }
            }
        } } }
    } } }
}"""

# Get options
families = parse_data_args(sys.argv, clustering=False)

df_positives = pd.DataFrame()
df_negatives = pd.DataFrame()

# Produce dataset for all families
for family in families:
    print(f"Fetching {family} data...")

    # Get binding sites
    sites = fetch_data(API_URL, FAMILY_SITES_QUERY, {"family": family})
    
    # Create dataset
    positives, negatives, failures = [], [], []
    with tqdm(total=len(sites) * 2) as pbar:
        
        # Positives
        for site in sites:
            try:

                atomium_site = get_atomium_site(site)
                sample_struct = structure_site_to_sample(atomium_site)

                sequence = site["chainInteractions"][0]["sequence"]
                sample_seq = sequence_site_to_sample(sequence)

                sample = {}
                for key, value in sample_struct.items():
                    sample[key] = value

                for key, value in sample_seq.items():
                    sample[key] = value

                
                if (sample["ca_max"] <= 30): 
                #     # positives.append(sample)

                    sample["pdb"] = site["id"]
                    sample["positive"] = 1
                     
                    df_positives = df_positives.append(sample, ignore_index=True)
                
            except Exception: 
                failures.append(site["id"])
            pbar.update()

        df_positives = df_positives.dropna()
        print(df_positives)
        
        # Negatives
        all_codes = get_all_pdb_codes()
        while df_negatives.shape[1] != df_positives.shape[1]:
            code = random.choice(all_codes)
            try:
                atomium_site = get_random_atomium_site(code, family)
                if site:
                    sample_struct = structure_site_to_sample(atomium_site)

                    sequence = site["chainInteractions"][0]["sequence"]
                    sample_seq = sequence_site_to_sample(sequence)

                    sample = {}
                    for key, value in sample_struct.items():
                        sample[key] = value

                    for key, value in sample_seq.items():
                        sample[key] = value   

                    if sample["ca_max"] <= 30:
                        # negatives.append(sample)
                        sample["pdb"] = site["id"]
                        sample["positive"] = 0
                        df_negatives = df_negatives.append(sample, ignore_index=True)
                        pbar.update()
            except Exception: pass

    # Report any sites which could not be turned into a sample
    if len(failures): print("Could not process:", ", ".join(failures))
    
    # Create CSV file for family
    save_csv(
        df_positives, df_negatives, family, None,
        os.path.join("data", "csv", "structure")
    )