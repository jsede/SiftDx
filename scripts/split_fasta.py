import sys
import os

def split_fasta(input_fasta, chunk_size):
    if not os.path.exists(input_fasta):
        print(f"Error: Input file {input_fasta} not found.")
        sys.exit(1)

    base = os.path.splitext(os.path.basename(input_fasta))[0]
    out_dir = os.path.join(os.path.dirname(input_fasta))

    chunk_index = 1
    seq_count = 0
    out_handle = None

    with open(input_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Start new chunk if needed
                if seq_count % chunk_size == 0:
                    if out_handle:
                        out_handle.close()
                    chunk_name = os.path.join(out_dir, f"{base}_chunk_{chunk_index}.fa")
                    out_handle = open(chunk_name, "w")
                    chunk_index += 1
                seq_count += 1
            out_handle.write(line)

    if out_handle:
        out_handle.close()

    print(f"Created {chunk_index-1} chunks in: {out_dir}")

if __name__ == "__main__":
    split_fasta(sys.argv[1], 20000)