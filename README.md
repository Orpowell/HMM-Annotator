# HMM-Annotator

HMM-annotator annotates DNA sequences using HMM profiles.

## Installation

HMM-annotator has only been tested using Python 3.10.14. Therefore, we recommend installing HMM-anntotar in an environment using conda or venv.

### Download HMM-Annotator

    git clone https://github.com/orpowell/HMM-annotator.git
    cd HMM-annotator

### Installion with venv

    python3 -m venv venv
    source activate venv/bin/activate
    pip install .

### Installation with conda (recommended)

    conda env create -n hmm-annotator python=3.10.14
    conda activate hmm-annotator
    pip install .


# Usage

HMM-Annotator can be run using the following command:

    hmm-annotator \
    --pfam profile.hmm \
    --sequences sequence.fasta \
    --cores 4 \
    --output annotations.txt \
    --bed annotations.bed 

## Parameters

|Parameter|Description|
|---|---|
|--pfam| Path to the HMM profile(s).|
|--sequences| Path to fasta file of sequences (genome or transcriptome).|
|--window| Window size used for splitting large sequences such as chromosomes (default and maximum: 100000 bp) |
|--overlap| Overlap between windows (Default: 100 bp) |
|--output| Path to output txt file of HMM annotations |
|--bed| Path to output txt file of HMM annotations. (optional) |
|--cores| Number of Cores. Note! This should a 25% of the total cores you wish to allocate to the tool.|

# Acknowledgements

Many thanks to Burkhard Steuernagel for his suggestions on how to implement the annotation method!

# Citing HMM-Annotator

TBA

# Contact

If you wish to discuss HMM-annotator, raise issues or collaborate please contact me [here](mailto:mail@oliverpowell.com)
