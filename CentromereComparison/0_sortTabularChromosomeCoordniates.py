import argparse
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='A simple program that will sort a tabular file with Chromosome coordinates')
    parser.add_argument('--input', help='File to process.')
    parser.add_argument('--output', help='The name of the file after processing.')
    parser.add_argument('--ascending', action='store_true', help='ascending sort ?')
    
    args = parser.parse_args()
    
    filename_without_ext_os = os.path.splitext(os.path.basename(args.input))[0]
    ext_os = os.path.splitext(os.path.basename(args.input))[1]
    df_input = pd.read_csv(args.input, sep="\t")
    assert set(["chrom", "chromStart", "chromEnd"]).issubset(df_input.columns), "header must contain 'chrom', 'chromStart', and 'chromEnd' column names"

    df_input["LENGTH"] = df_input.chromEnd - df_input.chromStart

    df_input.sort_values("LENGTH", ascending=args.ascending, inplace=True)
    df_input = df_input.drop('LENGTH', axis=1)
    output = args.output
    
    if not output:
        output = filename_without_ext_os+"_sorted"+ext_os

    df_input.to_csv(output, sep="\t", index=False)
    
    print(f"The file {output} was saved!")

if __name__ == "__main__":
    main()
