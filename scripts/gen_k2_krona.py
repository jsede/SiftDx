import os
import sys
import pandas as pd
import cleaners

def k2_krona(fasta_file, covstats_file, kraken_file):
    # get outdir
    dirpath = os.path.dirname(os.path.abspath(fasta_file))

    # prep the initial dataframe
    headers, lengths = cleaners.convert_fasta_to_list(fasta_file)
    fasta_df = pd.DataFrame({"Fasta_Headers": headers, "Seq_Length": lengths})
    covstats_df = cleaners.extract_name_numreads(covstats_file)
    fasta_df = fasta_df.merge(
        covstats_df, left_on="Fasta_Headers", right_on="#rname", how='left'
    )
    fasta_df['numreads'] = fasta_df['numreads'].fillna("1").astype(int)

    # prep the kraken dataframe
    kraken_df = pd.read_csv(
        kraken_file,
        sep="\t",
        header=None,
        names=["ID", "Kraken_Species"],
        usecols=[1, 2],
    )
    kraken_df["Kraken_taxid"] = kraken_df["Kraken_Species"].apply(cleaners.get_taxid).astype(int)
    kraken_df["Kraken_Species"] = kraken_df["Kraken_Species"].apply(
        cleaners.remove_brackets
    )
    desired_cols = ["ID", "Kraken_taxid", "Kraken_Species"]
    kraken_df = kraken_df.reindex(columns=desired_cols)

    # merge dataframes
    merged_df = fasta_df.merge(
        kraken_df, left_on="Fasta_Headers", right_on="ID", how="outer"
    )

    # collapse by taxid
    species_df = kraken_df[["Kraken_taxid", "Kraken_Species"]].drop_duplicates() # get unique species
    collapsed_df = merged_df.groupby(["Kraken_taxid"], as_index=False).agg(
        {"numreads": "sum"}
    )
    collapsed_df = collapsed_df.merge(
        species_df, left_on="Kraken_taxid", right_on="Kraken_taxid", how="outer"
    )

    desired_cols = ["Kraken_taxid", "Kraken_Species", "numreads"]
    collapsed_df["Kraken_Species"] = collapsed_df["Kraken_Species"].str.strip() # remove the trailing space in species name
    collapsed_df = collapsed_df.reindex(columns=desired_cols)
    report = dirpath + "/k2_krona_input.tsv"
    collapsed_df.to_csv(report, sep="\t", index=None, header=None)


if __name__ == "__main__":
    k2_krona(sys.argv[1], sys.argv[2], sys.argv[3])