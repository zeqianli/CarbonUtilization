import click
import pandas as pd
import gzip, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

@click.command()
@click.argument('input',nargs=1)
@click.argument('output',nargs=1)
def main(input,output):
    rast_df=pd.read_csv(input,sep="\t")
    rast_df=rast_df[rast_df['aa_sequence'].notna()]

    seqs=[]
    for _, row in rast_df.iterrows():
        seq=SeqRecord(Seq(row['aa_sequence']),id=row['location'],description=row['function'])
        seqs.append(seq)

    with gzip.open(output,'wt') as f:
        SeqIO.write(seqs,f,'fasta')

if __name__ == '__main__':
    main()