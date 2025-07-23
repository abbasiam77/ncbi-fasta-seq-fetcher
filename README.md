NCBI FASTA Protein Sequence Fetcher
===================================

A Python script to download protein FASTA sequences from NCBI Entrez by gene or protein name and species.
Supports standard gene/protein queries, multiple species, custom output folders, the option to fetch only the longest (canonical) isoform, and an NCBI API key for fast/bulk usage.

Author
------
Dr. Amir Ali Abbasi
National Center for Bioinformatics
Quaid-i-Azam University, Islamabad, Pakistan
Email: abbasiam@qau.edu.pk

Requirements
------------
- Python 3.x
- requests Python package

Install via:
    pip install requests
Or, to install all dependencies from requirements.txt:
    pip install -r requirements.txt

Usage
-----
python3 fetch_ncbi_fasta.py <gene> <species1> [<species2> ...] <number> <outfile_base> [--api-key API_KEY] [--outdir OUTDIR] [--all-fields] [--longest]

- Use double quotes for any species name with spaces.
- Output files will be auto-named and saved in the specified folder (default: current directory).

Examples
--------
1. Download up to 10 GLI3 protein isoforms for human and mouse:
    python3 fetch_ncbi_fasta.py GLI3 "Homo sapiens" "Mus musculus" 10 GLI3_Seqs --api-key YOUR_API_KEY --outdir results

2. Download only the single longest GLI3 protein isoform per species:
    python3 fetch_ncbi_fasta.py GLI3 "Homo sapiens" "Mus musculus" 10 GLI3_Seqs --api-key YOUR_API_KEY --outdir results --longest

3. For viral proteins or when gene name isn’t unique, add --all-fields:
    python3 fetch_ncbi_fasta.py nsp1 "Severe acute respiratory syndrome coronavirus 2" 5 NSP1_Seqs --all-fields

Options
-------
--api-key:   NCBI API key for higher request limits (Get yours here: https://www.ncbi.nlm.nih.gov/account/settings/)
--outdir:    Folder for saving results (created if not present)
--all-fields: Search all NCBI fields (useful for viral proteins or non-standard gene names)
--longest:   Download only the single longest (full-length) protein isoform per species

NCBI API Key
------------
For higher throughput and reliability, create your free NCBI API key (https://www.ncbi.nlm.nih.gov/account/settings/) and provide it with --api-key YOUR_API_KEY.

To install all requirements:
    pip install -r requirements.txt

Contact
-------
For questions, feedback, or contributions, please contact:

Dr. Amir Ali Abbasi
National Center for Bioinformatics
Quaid-i-Azam University, Islamabad, Pakistan
Email: abbasiam@qau.edu.pk
