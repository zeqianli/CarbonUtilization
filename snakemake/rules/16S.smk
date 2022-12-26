from Bio import SeqIO

rule barrnap:
    input: 
        #fna=os.path.join(config['DIR_fna'], '{sample}.fna')
        assembly=os.path.join(config['DIR_assembly'], '{sample}.fasta')
    output:
        sixteens=os.path.join(config['DIR_16S'],'{sample}.16S.fna'),
    conda: 'barrnap'
    threads: config['threads']
    shell:
        """
        barrnap --kingdom bac --threads {threads} --outseq {output.sixteens} {input.assembly}
        """

# samples=pd.read_csv("output/samples.txt",header=None).values.flatten()

# rule concatenate_16s:
#     # Pick the first 16S sequence for each sample
#     input:
#         sixteens=expand(os.path.join(config['DIR_16S'],'{sample}.16S.fna'),sample=samples)
#     output:
#         all_16s=config['FILE_16S']
#     params:
#         samples=samples
#     run:
#         records=[]
#         for ff_16S, sample in zip(input.sixteens, params.samples):
#             for record in SeqIO.parse(ff_16S, "fasta"):
#                 record.id = sample
#                 records.append(record)
#                 break
#         SeqIO.write(records, output.all_16s, "fasta")


# rule msa:
#     input:
#         sixteenS=config['FILE_16S']
#     output:
#         msa=config['FILE_16S_msa']
#     threads: config['threads']
#     conda: 'clustalo'
#     shell:
#         """
#         clustalo -i {input.sixteenS} -o {output.msa} --threads={threads}
#         """

# rule fasttree:
#     input:
#         msa=config['FILE_16S_msa']
#     output:
#         tree=config['FILE_16S_tree']
#     threads: config['threads']
#     conda: 'fasttree'
#     shell:
#         """
#         fasttree -nt -gtr -gamma -out {output.tree} {input.msa}
#         """