import argparse
import warnings
import os

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    import pandas as pd

except ValueError as e:
    print(e)


def global_identity(df):
    # ordenar intervalos por início
    hsps = []
    for _, row in df.iterrows():
        start, end = sorted([row.qstart, row.qend])
        hsps.append((start, end, row.pident, row.length))
    hsps.sort(key=lambda x: x[0])
    
    # unir intervalos e calcular identidade
    total_matches = 0
    total_bases = 0
    covered = []  # intervalos já mesclados
    
    for start, end, pident, length in hsps:
        matches = pident/100 * length
        
        # se o intervalo não sobrepõe com o último, adiciona direto
        if not covered or start > covered[-1][1]:
            covered.append([start, end, matches, length])
        else:
            # sobreposição → junta com o último
            last = covered[-1]
            overlap = min(end, last[1]) - start + 1
            
            # adicionar só a parte nova do intervalo
            new_len = length - overlap if overlap > 0 else length
            new_matches = (pident/100) * new_len
            
            last[1] = max(end, last[1])           # estende o intervalo
            last[2] += new_matches                # acumula matches sem duplicar
            last[3] += new_len                    # acumula length sem duplicar
    
    # soma os intervalos mesclados
    for _, _, matches, bases in covered:
        total_matches += matches
        total_bases += bases
    
    return total_matches / total_bases if total_bases > 0 else 0

def _save2Fasta(fasta_dict, output):

    out = open(output, "w")

    for k, v in fasta_dict.items():
        out.write(f">{k}\n{v}\n")
    out.close()

def _parseFastaFile(fasta):
    
    final = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        if ID not in final:
            final[ID] = SEQ
        else:
            warnings.warn(f"The ID {ID} is duplicated. Considering the first occurrence.")

    return final

def _parseBlastOutput(blastout6, output):
    
    cols = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "qcovhsp", "qcovs", "evalue", "bitscore"]
    df = pd.read_csv(blastout6, sep="\t", names=cols)
    
    df["qcov"] = (df["length"] / df["qlen"]) * 100   # Cobertura da query
    df["scov"] = (df["length"] / df["slen"]) * 100
    
    print(df.head())

    df.to_csv(f"{output}.tsv", sep="\t", index=False)

    VALUES = []
    qseqids = list(set(df.qseqid))
    sseqids = list(set(df.sseqid))
    for qid in qseqids:
        for sid in sseqids:
            mask = (df["qseqid"] == qid) & (df["sseqid"] == sid)
            if True in set(mask):
                pident = global_identity(df.loc[mask]) * 100
                VALUES.append({"qseqid": qid, "sseqid": sid, "pident": pident, "qcovs": df.loc[mask]["qcovs"].max()})
                #print(VALUES)
    df2 = pd.DataFrame(VALUES)
    df2.to_csv(f"{output}_filtered.tsv", sep="\t", index=False)

    
def _createDB(fasta: str = None):

    try:
        command = f"makeblastdb -in {fasta} -dbtype nucl -out blastDB/blastdb"
        os.system(command)

    except ValueError as e:
        print(e)


def _runBlastn(fasta, output, blastdb, cpu):

    try:
        os.system(f"blastn -query {fasta} -db {blastdb} -max_hsps 500 -out {output} -outfmt \'6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovhsp qcovs evalue bitscore\'")

    except ValueError as e:
        print(e)

def main():

    parser = argparse.ArgumentParser(description='A simple program that will sort a tabular file with Chromosome coordinates')
    parser.add_argument('--input', help='File to process.', type=str, default=None)
    parser.add_argument('--outname', help='The name of the file after processing without extension. Default is results', type=str, default='results')
    parser.add_argument("--cpu", help="Optional - number of threads to be used in each step default=1", type=int, default=1)
    
    args = parser.parse_args()
    
    filename_without_ext_os = os.path.splitext(os.path.basename(args.input))[0]
    ext_os = os.path.splitext(os.path.basename(args.input))[1]
    
    fasta_dict = _parseFastaFile(args.input)
    _save2Fasta(fasta_dict=fasta_dict, output="temp.fa")

    _createDB(fasta="temp.fa")

    _runBlastn(fasta="temp.fa", output=args.outname+".tsv", blastdb="blastDB/blastdb", cpu=args.cpu)
    
    _parseBlastOutput(blastout6=args.outname+".tsv", output=args.outname+"_coverage")
        
    print(f"The file {args.outname}_coverage.tsv was saved!")

if __name__ == "__main__":
    main()
