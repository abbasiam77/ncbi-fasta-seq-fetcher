#!/usr/bin/env python3
"""
Fetch protein FASTA sequences from NCBI Entrez for a gene (or protein) of interest in one or more species.

Usage:
  python fetch_ncbi_fasta.py <gene> <species1> [<species2> ...] <number> <outfile_base> [--api-key API_KEY] [--outdir OUTDIR] [--all-fields] [--longest]

Example:
  python fetch_ncbi_fasta.py GLI3 "Homo sapiens" "Mus musculus" 10 GLI3_Seqs --longest
"""

import requests
import os
import argparse
import sys

def parse_fasta(fasta):
    """Yield (header, sequence) tuples from a multi-FASTA string."""
    header, seq = None, []
    for line in fasta.splitlines():
        if line.startswith(">"):
            if header:
                yield (header, ''.join(seq))
            header = line
            seq = []
        else:
            seq.append(line.strip())
    if header:
        yield (header, ''.join(seq))

def fetch_proteins(gene, species, n, outfile_base, outdir=".", api_key=None, all_fields=False, longest=False):
    # Build the search term
    if all_fields:
        term = f"{gene} AND {species}[Organism]"
    else:
        term = f"{gene}[Gene] AND {species}[Organism]"
    # Search for protein IDs
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "protein",
        "term": term,
        "retmax": n,
        "retmode": "xml"
    }
    if api_key:
        params["api_key"] = api_key
    headers = {'User-Agent': 'NCBI-Seq-Fetcher/1.0 (your_email@example.com)'}
    print(f"Searching for up to {n} protein sequences of {gene} in {species}...")
    try:
        r = requests.get(base_url, params=params, headers=headers, timeout=20)
        r.raise_for_status()
    except Exception as e:
        print(f"Network or server error during ESearch: {e}")
        return
    ids = []
    for line in r.text.splitlines():
        if "<Id>" in line:
            ids.append(line.strip().replace("<Id>", "").replace("</Id>", ""))
    if not ids:
        print(f"No protein IDs found for {gene} in {species}.")
        return

    print(f"Found {len(ids)} IDs. Fetching sequences...")

    # Fetch the FASTA sequences
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fetch_params = {
        "db": "protein",
        "id": ",".join(ids),
        "rettype": "fasta",
        "retmode": "text"
    }
    if api_key:
        fetch_params["api_key"] = api_key
    try:
        fasta = requests.get(fetch_url, params=fetch_params, headers=headers, timeout=60).text
    except Exception as e:
        print(f"Error during sequence retrieval: {e}")
        return

    # If --longest is set, keep only the longest sequence
    if longest:
        seqs = list(parse_fasta(fasta))
        if not seqs:
            print("No sequences fetched.")
            return
        # Find the longest sequence
        header, seq = max(seqs, key=lambda x: len(x[1]))
        fasta = f"{header}\n"
        # Split long sequence into lines of 60 chars (standard FASTA format)
        fasta += '\n'.join([seq[i:i+60] for i in range(0, len(seq), 60)])

    # Build unique output filename
    safe_species = species.replace(" ", "_")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, f"{outfile_base}_{safe_species}.fasta")
    suffix = 1
    while os.path.exists(outfile):
        outfile = os.path.join(outdir, f"{outfile_base}_{safe_species}_{suffix}.fasta")
        suffix += 1

    # Write to output file
    with open(outfile, "w") as out:
        out.write(fasta)
    if longest:
        print(f"Downloaded longest protein FASTA sequence to {outfile}")
    else:
        print(f"Downloaded {len(ids)} protein FASTA sequences to {outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch protein FASTA sequences from NCBI Entrez for a gene or protein of interest in one or more species.")
    parser.add_argument("gene", help="Gene or protein name (query)")
    parser.add_argument("species", nargs="+", help="Species name(s), quoted if contains spaces")
    parser.add_argument("number", type=int, help="Max number of protein isoforms to download per species")
    parser.add_argument("outfile_base", help="Base name for output file(s)")
    parser.add_argument("--api-key", help="NCBI API key (optional, speeds up large queries)")
    parser.add_argument("--outdir", default=".", help="Output directory (default: current)")
    parser.add_argument("--all-fields", action="store_true", help="Search all fields (not just gene name, for viral proteins etc.)")
    parser.add_argument("--longest", action="store_true", help="Download only the single longest protein isoform (wild type)")
    args = parser.parse_args()

    gene = args.gene
    n = args.number
    outfile_base = args.outfile_base
    if len(args.species) > 2:
        species_list = args.species[:-2]
    else:
        species_list = args.species
    for species in species_list:
        fetch_proteins(gene, species, n, outfile_base, args.outdir, args.api_key, args.all_fields, args.longest)
