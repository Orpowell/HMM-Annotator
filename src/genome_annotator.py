import multiprocessing
import shutil
import os
import glob
from functools import partial
from chunk_annotator import GenomeChunk
from chopper import genome_chopper
import pandas as pd
import logging
import sys

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO)

def analyse_chunk(PFAM: str, tmp: str, chunk_path: str) -> None:
    try:
        chunk = GenomeChunk(chunk_path=chunk_path, tmp=tmp)
        chunk.run(PFAM)
    except ValueError:
        print(f"WARNING! Could not process {chunk_path} - likely Alphabet error...")

def merge_annotations(dir: str, output: str, bed: str|None):
    data_files = glob.glob(rf"{dir}/*.bed")

    with open(output, "w") as outfile:
        for fname in data_files:
            with open(fname) as infile:
                outfile.write(infile.read())
    
    df = pd.read_csv(output, sep="\t", header=None)
    df.sort_values(by=[0,1], inplace=True)
    df.to_csv(output, header=False, index=False, sep="\t", columns=[0,3,5,1,2])

    if bed is not None:
        df.reset_index(inplace=True, drop=True)
        df[3] = [f"{row[3]}_{index+1}" for index,row in df.iterrows()]
        df.to_csv(bed, header=False, index=False, sep="\t")

def annotate_genome(
    pfam: str, genome: str, window_size: int, overlap: int, cores: int, output: str, bed: str|None
) -> None:
    
    tmpdir = f"{os.getcwd()}/tmp_hmm_annotator"

    logging.info(f"Generating temporary directory... ({tmpdir})")

    os.makedirs(tmpdir, exist_ok=True)

    logging.info(f"Splitting sequences into chunks... ({window_size} bp with {overlap} bp overlaps)")

    genome_chopper(genome=genome, window_size=window_size, overlap=overlap, tmp=tmpdir)

    if cores > multiprocessing.cpu_count():

        available_cores = multiprocessing.cpu_count()

        requested_cores = round(available_cores / 4) - 1

        logging.info(
            f"More cores specified ({cores}) than available! A total of {requested_cores*4} cores will be used instead. To avoid this message, specify {requested_cores} for the --cores parameter nextime)"
        )
        
        if requested_cores < 1:
            cores = 1

        cores = requested_cores


    chunk_pool = glob.glob(rf"{tmpdir}/*.fasta")

    logging.info("Annotating chunks...")

    pool = multiprocessing.get_context("spawn").Pool(processes=cores)
    pool.map(partial(analyse_chunk, pfam, tmpdir), chunk_pool)

    logging.info("Merging annotations...")

    merge_annotations(dir=tmpdir, output=output, bed=bed)

    logging.info(f"Removing temporary directory... ({tmpdir})")

    try:
        shutil.rmtree(tmpdir)
    except OSError as e:
        logging.error(f"Error: {e.filename} - {e.strerror}.")
    
    logging.info("All processes finsihed...")