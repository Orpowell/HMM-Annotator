from Bio import SeqIO

def genome_chopper(genome: str, window_size: int, overlap: int, tmp: str) -> None:
    pusher: int = window_size - overlap

    records = SeqIO.parse(open(genome), "fasta")

    for record in records:
        chromosome: str = record.id
        sequence_length: int = len(record.seq)

        i: int = 0
        while i < sequence_length:
            start: int = i
            end: int = i + window_size

            if end > sequence_length:
                end = sequence_length

            chunk_name: str = f"{chromosome}-{start}-{end}"

            with open(f"{tmp}/{chunk_name}.fasta", "w+") as file:
                file.write(f">{chunk_name}\n")
                file.write(f"{str(record.seq[start:end])}\n")

            i += pusher