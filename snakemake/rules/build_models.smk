rule carveme:
    input:
        faa=os.path.join(config["DIR_faa"],'{sample}.faa')
    output:
        model=os.path.join(config['DIR_CARVEME'],'{sample}.xml'),
        model_gapfill=os.path.join(config['DIR_CARVEME'],'{sample}_gapfill.xml')
    params:
        media='M9',
    conda: 'fba'
    threads:
        config['threads']
    shell:
        """
        export PYTHONPATH=/home/zeqianli/scratch-midway3/bin/cplex/cplex/python/3.9/x86-64_linux 
        carve {input.faa} --output {output.model_gapfill} -g {params.media} -i "{params.media}"
        carve {input.faa} --output {output.model} 
        """
        # I wish I don't need to add the first line but UChicago cluster is a piece of fucking garbage and this is the only way I can make it work. 