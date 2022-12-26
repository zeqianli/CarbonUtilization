rule carveme:
    input:
        faa=os.path.join(config["DIR_faa"],'{sample}.faa')
    output:
        model=os.path.join(config['DIR_CARVEME'],'{sample}.xml'),
        model_gapfill=os.path.join(config['DIR_CARVEME'],'{sample}_gapfill.xml')
    params:
        media='M9',
        #media_file=config['FILE_CARVEME_MEDIA']
    conda: 'fba'
    threads:
        config['threads']
    shell:
        """
        export PYTHONPATH=/home/zeqianli/project/bin/cplex/cplex/python/3.9/x86-64_linux 
        carve {input.faa} --output {output.model_gapfill} -g {params.media} -i "{params.media}"
        carve {input.faa} --output {output.model} 
        """
        # I wish I don't need to add the first line but UChicago cluster is a piece of fucking garbage and this is the only way I can make it work. 
        
        # carve {input.faa} --output {output.model_gapfill} -g "{params.media}" -i "{params.media}" --mediadb {params.media_file} 


rule memote:
    input:
        # model=os.path.join(config['DIR_CARVEME'],'{sample}.xml'), 
        model_gapfill=os.path.join(config['DIR_CARVEME'],'{sample}_gapfill.xml')
    output:
        #report=os.path.join(config['DIR_MEMOTE'],'{sample}.html'),
        report_gapfill=os.path.join(config['DIR_MEMOTE'],'{sample}_gapfill.html'),
        #report_diff=os.path.join(config['DIR_MEMOTE'],'{sample}_diff.html')
    conda: 'fba'
    threads:
        config['threads']
    shell:
        """
        memote report snapshot --filename {output.report_gapfill} {input.model_gapfill} 
        """
        # memote report snapshot --filename {output.report} {input.model} 
        # memote report diff --filename {output.report_diff} {input.model} {input.model_gapfill}
