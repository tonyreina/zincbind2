#! /usr/bin/env python3

import sys
sys.path.append("../zincbindpredict")
import kirjava
import random
from tqdm import tqdm
import subprocess
import re
from common import sequence_site_to_sample
from utilities import *

API_URL = "https://api.zincbind.net/"

FAMILY_CHAINS_QUERY = """query familySites($family: String) {
    zincsites(family: $family) { edges { node { 
        id pdb { id assembly } chainInteractions { edges { node { sequence } } }
    } } }
}"""

# Get options
families, clustering = parse_data_args(sys.argv)

# Produce dataset for all families
for family in families:
    print(f"Fetching {family} data...")

    # Get sequence per binding site
    sites = fetch_data(API_URL, FAMILY_CHAINS_QUERY, {"family": family})
    pdbs = [site["pdb"] for site in sites if\
        len(site["chainInteractions"]) == 1]
    sequences = [site["chainInteractions"][0]["sequence"] for site in sites if\
        len(site["chainInteractions"]) == 1]
    sequences = [s for s in sequences if sequence_contains_family(s, family)]

    # Create dataset for each cluster
    for similarity in clustering:
        if similarity:
            print("similarity")
            clustered_sequences = cluster_sequences(sequences, similarity)
        else:
            print("no similarity")
            clustered_sequences = sequences[:]

        # Create dataset
        positives, negatives = [], []
        print(f"{family} (clustering={similarity})")
        with tqdm(total=len(clustered_sequences) * 2) as pbar:

            # Positives
            for idx, sequence in enumerate(clustered_sequences):
                sample = sequence_site_to_sample(sequence)

                sample["pdb"] = pdbs[idx]["id"]

                positives.append(sample)
                pbar.update()

            # Negatives
            all_sequences = get_all_uniprot_sequences()
            while len(negatives) != len(positives):
                sequence = random.choice(all_sequences).lower()
                site = get_random_sequence_site(sequence, family)
                if site:
                    sample = sequence_site_to_sample(site)
                    sample["pdb"] = site["pdb"]["id"]

                    negatives.append(sample)
                    pbar.update()
        
        # Create CSV file for family and similarity
        save_csv(
            positives, negatives, family, similarity,
            os.path.join("data", "csv", "sequence")
        )