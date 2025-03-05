import os
from Bio import SeqIO
import pyhmmer
import pandas as pd
import numpy as np

pd.options.mode.copy_on_write = True


class GenomeChunk:
    def __init__(self, chunk_path: str, tmp: str) -> None:
        record = next(SeqIO.parse(chunk_path, "fasta"))

        chunk_info = record.id.split("-")

        self.chunk = record.id
        self.chromosome: str = chunk_info[0]
        self.start: int = int(chunk_info[1])
        self.end: int = int(chunk_info[2])
        self.forward_sequence = record.seq
        self.reverse_sequence = record.seq.reverse_complement()
        self.tmp = tmp
        self.forward = f"{self.tmp}/{self.chunk}_forward.fasta"
        self.reverse = f"{self.tmp}/{self.chunk}_reverse.fasta"

    def translate_ORFs(self):
        with open(self.forward, "w+") as out:
            out.write(">f-0\n")
            out.write(f"{self.forward_sequence.translate()}\n")
            out.write(">f-1\n")
            out.write(f"{self.forward_sequence[1:].translate()}\n")
            out.write(">f-2\n")
            out.write(f"{self.forward_sequence[2:].translate()}\n")

        with open(self.reverse, "w+") as out:
            out.write(">r-0\n")
            out.write(f"{self.reverse_sequence.translate()}\n")
            out.write(">r-1\n")
            out.write(f"{self.reverse_sequence[1:].translate()}\n")
            out.write(">r-2\n")
            out.write(f"{self.reverse_sequence[2:].translate()}\n")

    def run_pyhmmer_hmmsearch(self, input, pfam) -> list[pyhmmer.plan7.Hit]:
        pressed = [f"{pfam}.h3p", f"{pfam}.h3m", f"{pfam}.h3i", f"{pfam}.h3f"]

        if not all(map(os.path.isfile, pressed)):
            pfam_load = pyhmmer.plan7.HMMFile(pfam)
            pyhmmer.hmmer.hmmpress(pfam_load, pfam)

        with pyhmmer.easel.SequenceFile(input, digital=True) as seq_file:
            sequences = list(seq_file)

        with pyhmmer.plan7.HMMFile(pfam) as hmm_file:
            targets = hmm_file.optimized_profiles()

            results = list(
                pyhmmer.hmmsearch(sequences=sequences, queries=targets, cpus=4)
            )

        return results

    def generate_annotation_info(
        self, results: list[pyhmmer.plan7.Hit], strand: bool
    ) -> pd.DataFrame:
        data_array = []

        # Forward strand processing
        if strand:
            for hits in results:
                for hit in hits.reported:
                    for domain in hit.domains.reported:
                        start = domain.env_from
                        end = domain.env_to
                        name = domain.alignment.hmm_name.decode()
                        score = domain.i_evalue

                        chromStart = (int(start) * 3) + self.start - 2
                        chromEnd = (int(end) * 3) + self.start
                        strand = "+"

                        row = [
                            self.chromosome,
                            self.chromosome,
                            name,
                            score,
                            chromStart,
                            chromEnd,
                            strand,
                        ]
                        data_array.append(row)

        # Reverse strand processing
        else:
            for hits in results:
                for hit in hits.reported:
                    for domain in hit.domains.reported:
                        start = domain.env_from
                        end = domain.env_to
                        name = domain.alignment.hmm_name.decode()
                        score = domain.i_evalue

                        chromStart = self.end - (int(end) * 3) - 2
                        chromEnd = self.end - (int(start) * 3)
                        strand = "-"

                        row = [
                            self.chromosome,
                            self.chromosome,
                            name,
                            score,
                            chromStart,
                            chromEnd,
                            strand,
                        ]
                        data_array.append(row)

        return pd.DataFrame(data_array)

    def clean_up(self):
        os.remove(self.forward)
        os.remove(self.reverse)

    def run(self, PFAM):
        self.translate_ORFs()

        # Forward
        forward_hits = self.run_pyhmmer_hmmsearch(self.forward, PFAM)
        forward_annotations = self.generate_annotation_info(forward_hits, True)

        # Reverse
        reverse_hits = self.run_pyhmmer_hmmsearch(self.reverse, PFAM)
        reverse_annotations = self.generate_annotation_info(reverse_hits, False)

        # Merge forward and reverse strand annotations
        df = pd.concat([forward_annotations, reverse_annotations], axis=0)

        if len(df) == 0:
            print(f"No annotations in : {self.chunk}")

        else:
            df.drop_duplicates(subset=[0, 2, 4, 5, 6], inplace=True, keep="first")
            df[9] = np.where(df[6] == "+", "255,0,0", "0,0,225")

            filtered_df = df[df[3] < 1e-3]
            filtered_df.sort_values(by=4, inplace=True)
            filtered_df.to_csv(
                f"{self.tmp}/{self.chunk}.bed",
                sep="\t",
                columns=[0, 4, 5, 2, 3, 6, 4, 5, 9],
                index=False,
                header=False,
            )

        self.clean_up()
