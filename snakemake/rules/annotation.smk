rule prodigal:
    input:
        contig=os.path.join(config['DIR_assembly'],'{sample}.fasta')
    output:
        faa=os.path.join(config['DIR_faa'],'{sample}.faa'),
        gff=os.path.join(config['DIR_faa'],'{sample}.gff'),
        fna=os.path.join(config['DIR_fna'],'{sample}.fna')  
    params:
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
