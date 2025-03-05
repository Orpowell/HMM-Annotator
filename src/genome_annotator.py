import multiprocessing
import shutil
import os
import glob
from functools import partial
from chunk_annotator import GenomeChunk
from chopper import genome_chopper


import pandas as pd

def analyse_chunk(PFAM: str, tmp: str, chunk_path: str) -> None:
    try:
        chunk = GenomeChunk(chunk_path=chunk_path, tmp=tmp)
        chunk.run(PFAM)
    except ValueError:
        print(f"WARNING! Could not process {chunk_path} - likely Alphabet error...")

def merge_annotations(dir: str, output: str):
    data_files = glob.glob(rf"{dir}/*.bed")

    with open(output, "w") as outfile:
        for fname in data_files:
            with open(fname) as infile:
                outfile.write(infile.read())
    
    df = pd.read_csv(output, sep="\t", header=None)
    df.sort_values(by=[0,1], inplace=True)
    print(df.head(10))
    df[3] = [f"{i}_bf_{n}" for n,i in enumerate(df[0])]
    print(df.head(10))
    df.to_csv(output, header=False, index=False, sep="\t")

def annotate_genome(
    pfam: str, genome: str, window_size: int, overlap: int, cores: int
) -> None:
    
    tmpdir = f"{os.getcwd()}/tmp_hmm_annotator"

    os.mkdir(tmpdir)

    genome_chopper(genome=genome, window_size=window_size, overlap=overlap, tmp=tmpdir)

    if cores > multiprocessing.cpu_count():
        print(
            f"Too many cores specified {cores}! {multiprocessing.cpu_count()} cores will be used"
        )
        cores = multiprocessing.cpu_count()

    chunk_pool = glob.glob(rf"{tmpdir}/*.fasta")

    pool = multiprocessing.get_context("spawn").Pool(processes=cores)
    pool.map(partial(analyse_chunk, pfam, tmpdir), chunk_pool)

    merge_annotations(dir=tmpdir, output="annotations.bed")

    try:
        shutil.rmtree(tmpdir)
    except OSError as e:
        print(f"Error: {e.filename} - {e.strerror}.")