import argparse
import pandas as pd
import numpy as np
import os

def __has_required_cols(d):
        return {"qseqid", "sseqid", "pident", "qcovs"}.issubset(set(d.columns))

# função para estimar cobertura (fração 0..1)
def __estimate_coverage_fraction(row):
    if pd.notna(row.get("qcovs")):
        return float(row["qcovs"]) / 100.0   
    return np.nan

def _getMatrix(df, bed):

    # converte colunas numéricas se existirem
    for col in ["pident", "qcovs"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df["pident_frac"] = pd.to_numeric(df["pident"], errors="coerce") / 100.0
    df["coverage_frac"] = df.apply(__estimate_coverage_fraction, axis=1)

    # lista de IDs (união de queries e subjects), ordenada
    ids = [f"{c}:{s}-{e}" for c, s, e in  bed[["chrom", "chromStart", "chromEnd"]].itertuples(index=False, name=None)]

    # inicializa matriz com NaN
    matrix = pd.DataFrame(0.0, index=ids, columns=ids, dtype=float)

    # mapa id -> posição
    id_to_pos = {id_: pos for pos, id_ in enumerate(ids)}

    for _, row in df.iterrows():
        a = str(row["qseqid"])
        b = str(row["sseqid"])
        pid = row.get("pident_frac", 0.0)
        cov = row.get("coverage_frac", 0.0)
        i = id_to_pos[a]
        j = id_to_pos[b]
        if i <= j:
            r_upper, c_upper = ids[i], ids[j]
            r_lower, c_lower = ids[j], ids[i]
        else:
            r_upper, c_upper = ids[j], ids[i]
            r_lower, c_lower = ids[i], ids[j]
        
        if pd.notna(pid):
            prev = matrix.at[r_upper, c_upper]
            if pd.isna(prev) or pid > prev:
                matrix.at[r_upper, c_upper] = float(pid)
        
        if pd.notna(cov):
            prev = matrix.at[r_lower, c_lower]
            if pd.isna(prev) or cov > prev:
                matrix.at[r_lower, c_lower] = float(cov)

    # diagonal = 1.0
    for idv in ids:
        matrix.at[idv, idv] = 1.0

    return matrix
    
def main():
    parser = argparse.ArgumentParser(description='A simple program that will create a Square matrix with coverage and identity from centromeres blastn results.')
    parser.add_argument('--input', help='File to process.')
    parser.add_argument('--outname', help='The basename of the file after processing.', default="results")
    parser.add_argument('--sorted_bed_file', help='a sorted bed file with chrom, chromStart, and chromEnd column names')

    args = parser.parse_args()

    sorted_bed_file = pd.read_csv(args.sorted_bed_file, sep="\t")
    assert set(["chrom", "chromStart", "chromEnd"]).issubset(sorted_bed_file.columns), "header must contain 'chrom', 'chromStart', and 'chromEnd' column names"

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Not found: {args.input}")

    df = pd.read_csv(args.input, sep="\t")
    
    if not __has_required_cols(df):
        raise ValueError("Missing columns: qseqid, sseqid, pident, qcovs")
    
    matrix = _getMatrix(df=df, bed=sorted_bed_file)
    matrix.to_csv(f"{args.outname}.tsv", sep="\t", float_format="%.3f")

    print(f"Matrix {args.outname}.tsv was saved !")

if __name__ == "__main__":
    main()
