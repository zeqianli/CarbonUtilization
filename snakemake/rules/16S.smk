from Bio import SeqIO

rule barrnap:
    input: 
        #fna=os.path.join(config['DIR_fna'], '{sample}.fna')
        assembly=os.path.join(config['DIR_assembly'], '{sample}.fasta')
    output:
        sixteens=os.path.join(config['DIR_16S'],'{sample}.16S.fna'),
    conda: 'sixteenS'
    threads: config['threads']
    shell:
        """
        barrnap --kingdom bac --threads {threads} --outseq {output.sixteens} {input.assembly}
        """