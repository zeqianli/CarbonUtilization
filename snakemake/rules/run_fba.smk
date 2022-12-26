import re 

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
        threshold=config['fba_growth_threshold'],
        initial_medium=config['FILE_INITIAL_MEDIUM']
    conda: 'fba'
    shell:
        """ 
        python scripts/run_fba.py --output {output.growth_gapfill} --initial_medium {params.initial_medium} {input.model_gapfill}
        python scripts/run_fba.py --output {output.growth_gapfill_force_uptake} --force_uptake --initial_medium {params.initial_medium} {input.model_gapfill}
        python scripts/run_fba.py --output {output.growth} --initial_medium {params.initial_medium} {input.model}
        python scripts/run_fba.py --output {output.growth_force_uptake} --force_uptake --initial_medium {params.initial_medium} {input.model}
        """
        

# samples=pd.read_csv('output/samples.txt',header=None).values.flatten()

# rule concatenate_results:
#     input:
#         growth_gapfill=expand(os.path.join(config['DIR_FBA_PREDICTION'],'{sample}_gapfill.csv'),sample=samples),
#         growth_gapfill_force_uptake=expand(os.path.join(config['DIR_FBA_PREDICTION'],'{sample}_gapfill_force_uptake.csv'),sample=samples)
#     output:
#         matrix_gapfill=os.path.join(config['DIR_FBA_PREDICTION'],'growth_matrix_gapfill.csv'),
#         matrix_gapfill_force_uptake=os.path.join(config['DIR_FBA_PREDICTION'],'growth_matrix_gapfill_force_uptake.csv')
#     run:
#         def concatenate(ffs,output):
#             growths=[pd.read_csv(ff,index_col=0,names=[re.findall('/([^/]*)_[\w_]+.csv',ff)[0]]) for ff in ffs]
#             growth_matrix=pd.concat(growths,axis=1)
#             growth_matrix.to_csv(output)
#         print(input.growth_gapfill,output.matrix_gapfill)
#         concatenate(input.growth_gapfill,output.matrix_gapfill)
#         concatenate(input.growth_gapfill_force_uptake,output.matrix_gapfill_force_uptake)

# TODO: generate a binary matrix of growth