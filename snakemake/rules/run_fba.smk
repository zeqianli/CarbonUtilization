rule run_fba:
    input:
        model=os.path.join(config['DIR_CARVEME'],'{sample}.xml'),
        model_gapfill=os.path.join(config['DIR_CARVEME'],'{sample}_gapfill.xml')
    output:
        growth=os.path.join(config['DIR_FBA_PREDICTION'],'{sample}.csv'),
        growth_force_uptake=os.path.join(config['DIR_FBA_PREDICTION'],'{sample}_force_uptake.csv'),
        growth_gapfill=os.path.join(config['DIR_FBA_PREDICTION'],'{sample}_gapfill.csv'),
        growth_gapfill_force_uptake=os.path.join(config['DIR_FBA_PREDICTION'],'{sample}_gapfill_force_uptake.csv')
    params:
        initial_medium=config['FILE_INITIAL_MEDIUM']
    conda: 'fba'
    shell:
        """ 
        python scripts/run_fba.py --output {output.growth_gapfill} --initial_medium {params.initial_medium} {input.model_gapfill}
        python scripts/run_fba.py --output {output.growth_gapfill_force_uptake} --force_uptake --initial_medium {params.initial_medium} {input.model_gapfill}
        python scripts/run_fba.py --output {output.growth} --initial_medium {params.initial_medium} {input.model}
        python scripts/run_fba.py --output {output.growth_force_uptake} --force_uptake --initial_medium {params.initial_medium} {input.model}
        """
        
