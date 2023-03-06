rule carveme:
    input:
        faa=os.path.join(config["DIR_faa"],'{sample}.faa')
    output:
        model=os.path.join(config['DIR_CARVEME'],'{sample}.xml'),
        model_gapfill=os.path.join(config['DIR_CARVEME'],'{sample}_gapfill.xml')
    params:
        media='M9',
        PATH_CPLEX=config['PATH_CPLEX']
    conda: 'fba'
    threads:
        config['threads']
    shell:
        """
        export PYTHONPATH="$PYTHONPATH:{params.PATH_CPLEX}"
        carve {input.faa} --output {output.model_gapfill} -g {params.media} -i "{params.media}"
        carve {input.faa} --output {output.model} 
        """
