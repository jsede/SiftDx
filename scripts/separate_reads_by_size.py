import os
import sys
import gzip

def open_file(file_path, mode='r'):
    """
    Opens a file, handling both plain text and gzip-compressed files.
    """
    if file_path.endswith('.gz') or is_gzipped(file_path):
        return gzip.open(file_path, mode + 't')  # Text mode for gzip
    else:
        return open(file_path, mode)

def is_gzipped(file_path):
    """
    Checks if a file is gzip-compressed by reading its magic number.
    """
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def separate_reads_by_size(input_reads, longer_reads, shorter_reads, min_len, max_len):
    """
    Splits a FASTQ file into two based on read length:
    - Reads with specified length go to longer_reads
    - Reads with specified length go to shorter_reads
    """

    with open_file(input_reads, 'r') as infile, \
         open(longer_reads, 'w') as out_long, \
         open(shorter_reads, 'w') as out_short:

        while True:
            header = infile.readline()
            if not header:
                break  # EOF

            seq = infile.readline()
            plus = infile.readline()
            qual = infile.readline()

            read_length = len(seq.strip())

            # Discard reads smaller than min_len
            if read_length < min_len:
                continue

            # Long reads (>= max)
            elif read_length >= max_len:
                out_long.write(header + seq + plus + qual)

            # Between min_len and max_len
            else:
                out_short.write(header + seq + plus + qual)

if __name__ == "__main__":
    separate_reads_by_size(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))