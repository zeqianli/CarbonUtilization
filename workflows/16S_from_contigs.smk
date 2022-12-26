from Bio import SeqIO

configfile: "/home/zeqianli/Kuehn/Carbon/workflows/16S.yml"

samples=[ff.replace('.fsa','') for ff in os.listdir(config['DIR_assembly']) if ff.endswith(".fsa")]

rule all:
    input:
        expand(os.path.join(config['DIR_16S'],"{sample}.16S.fna"), sample=samples),
        config['FILE_16S'],
        config['FILE_16S_msa'],
        config['FILE_16S_tree']

# rule prodigal:
#     input:
#         assembly=os.path.join(config['DIR_assembly'],'{sample}.fsa')
#     output:
#         fna=os.path.join(config['DIR_prodigal'],'{sample}.fna'),
#         gbk=os.path.join(config['DIR_prodigal'],'{sample}.gbk')
#     threads: 1
#     shell:
#         "prodigal -i {input.assembly} -o {output.gbk} -d {output.fna}"
        

rule barrnap:
    input: 
        #fna=os.path.join(config['DIR_prodigal'], '{sample}.fna')
        fna=os.path.join(config['DIR_assembly'], '{sample}.fsa')
    output:
        sixteens=os.path.join(config['DIR_16S'],'{sample}.16S.fna'),
    conda: 'barrnap'
    threads: config['threads']
    shell:
        """
        barrnap --kingdom bac --threads {threads} --outseq {output.sixteens} {input.fna}
        """

rule concatenate_16s:
    # Pick the first 16S sequence for each sample
    input:
        sixteens=expand(os.path.join(config['DIR_16S'],'{sample}.16S.fna'),sample=samples)
    output:
        all_16s=config['FILE_16S']
    params:
        samples=samples
    run: # TODO: check that the first one is actually 16S. 
        records=[]
        for ff_16S, sample in zip(input.sixteens, params.samples):
            for record in SeqIO.parse(ff_16S, "fasta"):
                record.id = sample
                records.append(record)
                break
        SeqIO.write(records, output.all_16s, "fasta")

rule sina:
    input:
        sixteenS=config['FILE_16S']
    output:
        msa=config['FILE_16S_msa']
    params:
        db=config['DB_SILVA']
    threads: config['threads']
    conda: 'sixteenS'
    shell:
        """
        sina -i {input.sixteenS} -o {output.msa} -r {params.db} --search
        """

rule fasttree:
    input:
        msa=config['FILE_16S_msa']
    output:
        tree=config['FILE_16S_tree']
    threads: config['threads']
    conda: 'sixteenS'
    shell:
        """
        fasttree -nt -gtr -gamma -out {output.tree} {input.msa}
        """