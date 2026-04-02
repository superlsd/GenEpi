import numpy as np
import networkx as nx
import pickle as pk
import pandas as pd
from Bio import Entrez
import time
import re
from collections import Counter

Entrez.email = "yourmail@gmail.com"
Entrez.tool = "drug_project"
# Entrez.api_key = "YOUR_KEY"  # recommended for large runs (10 req/sec vs 1 req/sec)

# -----------------------------
# Constants
# -----------------------------
GENERIC_SINGLE_WORDS = {
    "simple", "other", "unspecified", "related", "associated",
    "mixed", "acute", "chronic", "secondary", "primary", "pain",
    "mental", "general", "multiple", "various", "certain", "due"
}

# -----------------------------
# Core helpers
# -----------------------------
def contains_boolean_words(term):
    return bool(re.search(r'\b(and|or|not)\b', term, flags=re.IGNORECASE))

def split_on_boolean_words(term):
    """Split a phrase on AND/OR/NOT, returning only the text chunks."""
    parts = re.split(r'\b(and|or|not)\b', term, flags=re.IGNORECASE)
    chunks = [p.strip() for p in parts
              if not re.match(r'^(and|or|not)$', p.strip(), re.IGNORECASE)]
    return [c for c in chunks if c]

def is_meaningful_fragment(fragment):
    """
    Multi-word fragments are always kept.
    Single-word fragments are kept only if clinically specific
    (i.e. not in the generic stoplist).
    """
    words = fragment.split()
    if len(words) > 1:
        return True
    return words[0].lower() not in GENERIC_SINGLE_WORDS

# -----------------------------
# Query builders
# -----------------------------
def build_disease_clause(disease_name):
    if not contains_boolean_words(disease_name):
        # Unquoted, no field tag → automatic term mapping, mirrors web interface
        return f'{disease_name}'

    # Boolean words present: split into fragments and OR them explicitly.
    # This prevents PubMed from parsing the raw 'and'/'or' words in the disease
    # string as uncontrolled Boolean operators with wrong precedence.
    # Each fragment is still unquoted → automatic term mapping applies per fragment.
    parts = split_on_boolean_words(disease_name)
    meaningful_parts = [p for p in parts if is_meaningful_fragment(p)]
    if not meaningful_parts:
        meaningful_parts = parts

    clauses = ' OR '.join(f'{p}' for p in meaningful_parts)
    return f'({clauses})'
def build_drug_clause(drug_name):
    # Unquoted, no field tag → automatic term mapping
    return f'{drug_name}'

def build_pubmed_query(drug_name, disease_name):
    drug_q    = build_drug_clause(drug_name)
    disease_q = build_disease_clause(disease_name)
    return f'({drug_q} AND {disease_q})'

# -----------------------------
# PubMed fetch
# -----------------------------
def get_pubmed_count(query):
    handle = Entrez.esearch(db="pubmed", term=query, retmode="xml", retmax=0)
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])

# -----------------------------
# Spot-check helper
# -----------------------------
def spot_check(drug_name, disease_name):
    """Print the query and its count for manual verification."""
    query = build_pubmed_query(drug_name, disease_name)
    count = get_pubmed_count(query)
    print(f"Drug   : {drug_name}")
    print(f"Disease: {disease_name}")
    print(f"Query  : {query}")
    print(f"Count  : {count}")
    print()
    return count

# -----------------------------
# Load pairs
# -----------------------------
pairs = set()

with open("disease_pairs_for_pubmed.txt", "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        drug, disease = parts
        pairs.add((drug, disease))

print(f"Loaded {len(pairs)} unique drug-disease pairs")


# -----------------------------
# Main loop
# -----------------------------
drugbank_lit_associations_chem_dict = {}

print("=" * 60)
print("MAIN RUN")
print("=" * 60)

for idx, (drug, disease) in enumerate(pairs):
    try:
        query = build_pubmed_query(drug, disease)
        count = get_pubmed_count(query)
        drugbank_lit_associations_chem_dict[drug, disease] = count
        time.sleep(0.34)
    except Exception as e:
        print(f"ERROR [{drug} | {disease}]: {e}")
        drugbank_lit_associations_chem_dict[drug, disease] = None

    if (idx + 1) % 100 == 0:
        print(f"  Processed {idx + 1}/{len(pairs)} pairs...")

# -----------------------------
# Save
# -----------------------------
with open('drugbank_lit_associations_chem_dict_claude.pickle', 'wb') as handle:
    pk.dump(drugbank_lit_associations_chem_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

print("\nSaved to drugbank_lit_associations_chem_dict_claude.pickle")