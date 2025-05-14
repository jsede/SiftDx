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

def separate_reads_by_size(input_reads, longer_reads, shorter_reads):
    """
    Splits a FASTQ file into two based on read length:
    - Reads with length >= 100 bp go to longer_reads
    - Reads with length < 100 bp go to shorter_reads
    """
    with open_file(input_reads, 'r') as infile, \
         open(longer_reads, 'w') as out_ge_100, \
         open(shorter_reads, 'w') as out_lt_100:

        while True:
            header = infile.readline()
            if not header:
                break  # EOF
            seq = infile.readline()
            plus = infile.readline()
            qual = infile.readline()

            read_length = len(seq.strip())

            if read_length >= 100:
                out_ge_100.write(header + seq + plus + qual)
            else:
                out_lt_100.write(header + seq + plus + qual)

if __name__ == "__main__":
    separate_reads_by_size(sys.argv[1], sys.argv[2], sys.argv[3])