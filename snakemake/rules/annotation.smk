# rule bakta: 
#     input:
#         contig=os.path.join(config['DIR_assembly'],'{sample}.fasta')
#     output:
#         bakta=directory(os.path.join(config['DIR_bakta'],'{sample}')),
#         faa=os.path.join(config['DIR_bakta'],'{sample}/{sample}.faa'),
#         fna=os.path.join(config['DIR_bakta'],'{sample}/{sample}.fna')  
#     params:
#         db=config['PATH_bakta_db'],
#         sample='{sample}'
#     conda: 'bakta'
#     threads: config['threads']
#     shell:
#         """
#         bakta --db {params.db} --output {output.bakta} --threads {threads} --prefix {params.sample} {input.contig}
#         """

# rule copy_faa:
#     input:
#         faa=os.path.join(config['DIR_bakta'],'{sample}/{sample}.faa')
#     output:
#         faa=os.path.join(config['DIR_faa'],'{sample}.faa')
#     shell:
#         """
#         cp {input.faa} {output.faa}
#         """

# rule copy_fna:
#     input:
#         fna=os.path.join(config['DIR_bakta'],'{sample}/{sample}.fna')
#     output:
#         fna=os.path.join(config['DIR_fna'],'{sample}.fna')
#     shell:
#         """
#         cp {input.fna} {output.fna}
#         """


# rule unzip_faa:
#     input: 
#         faa=os.path.join(config["DIR_faa"],'{sample}.faa.gz')
#     output:
#         faa_unzipped=os.path.join(config["DIR_faa"],'{sample}.faa')
#     shell:
#         "gunzip -c {input.faa} > {output.faa_unzipped}"

rule prodigal:
    input:
        contig=os.path.join(config['DIR_assembly'],'{sample}.fasta')
    output:
        # bakta=directory(os.path.join(config['DIR_bakta'],'{sample}')),
        faa=os.path.join(config['DIR_faa'],'{sample}.faa'),
        gff=os.path.join(config['DIR_faa'],'{sample}.gff'),
        fna=os.path.join(config['DIR_fna'],'{sample}.fna')  
    params:
        #db=config['PATH_bakta_db'],
        sample='{sample}'
    conda: 'prodigal'
    threads: 1
    shell:
        """
        prodigal -i {input.contig} -o {output.gff} -f gff -a {output.faa} -d {output.fna}
        """

rule kofamscan:
    input:
        faa_unzipped=os.path.join(config["DIR_faa"],'{sample}.faa')
    output:
        #model=os.path.join(config['DIR_CARVEME'],'{sample}.xml'),
        kofamscan=os.path.join(config['DIR_KOFAMSCAN'],'{sample}.ko')
    params:
        profiles=config['PATH_KOFAMSCAN_PROFILES'],
        ko_list=config['PATH_KOFAMSCAN_KO_LIST'],
        temp_dir=config['DIR_TEMP']
    conda: 'kofamscan'
    threads: config['threads']
    shell:
        """
        exec_annotation -o {output.kofamscan} -p {params.profiles} -k {params.ko_list} -f mapper --cpu={threads} --tmp-dir={params.temp_dir} {input.faa_unzipped}
        """
