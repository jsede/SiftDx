import sys
import os

MB = 1024 * 1024


def decide_chunks(file_size):
    """Return number of splits based on size."""
    if file_size < 10 * MB:
        return 1
    elif file_size <= 30 * MB:
        return 2
    else:
        return 3


def split_fastq(input_fastq):
    if not os.path.exists(input_fastq):
        print(f"Error: {input_fastq} not found")
        sys.exit(1)

    file_size = os.path.getsize(input_fastq)
    n_chunks = decide_chunks(file_size)

    print(f"{input_fastq} → {file_size} bytes → {n_chunks} chunk(s)")

    # No split needed
    if n_chunks == 1:
        print("File < 10 MB. No splitting required.")
        return

    # Count reads
    total_reads = sum(1 for _ in open(input_fastq)) // 4
    reads_per_chunk = total_reads // n_chunks

    base = os.path.splitext(os.path.basename(input_fastq))[0]
    out_dir = os.path.dirname(input_fastq)

    handles = []
    for i in range(n_chunks):
        out_file = os.path.join(out_dir, f"{base}_part_{i+1}.fastq")
        handles.append(open(out_file, "w"))

    current_chunk = 0
    read_counter = 0

    with open(input_fastq, "r") as f:
        while True:
            read = [f.readline() for _ in range(4)]
            if not read[0]:
                break

            for line in read:
                handles[current_chunk].write(line)

            read_counter += 1

            if read_counter >= reads_per_chunk and current_chunk < n_chunks - 1:
                current_chunk += 1
                read_counter = 0

    for h in handles:
        h.close()

    print("Splitting complete.")


if __name__ == "__main__":
    split_fastq(sys.argv[1])