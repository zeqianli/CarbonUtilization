import click
import pandas as pd, numpy as np
import cobra
import yaml

# load bigg metabolite database
def load_bigg():
    with open('config/config.yml') as f:
        config=yaml.load(f,Loader=yaml.FullLoader)
    bigg = pd.read_csv(config['FILE_BIGG_METABOLITES'],sep = '\t')
    bigg=bigg.drop_duplicates(subset=['bigg_id']) # only one duplicate bigg_id "(GlcNAc)7 (Man)3 (Asn)1"
    bigg=bigg.set_index('bigg_id',drop=False,verify_integrity=True)
    return bigg

bigg=load_bigg()


def load_fba_model(model):
    return cobra.io.read_sbml_model(model)

def find_c_condidates(model):
    # Carbon source candiates: all bigg metabolites (exchange, external; EX_[...]_e) that are not already in the model medium. 
    bigg_e=('EX_'+bigg[bigg['bigg_id'].str.endswith('_e')]['bigg_id']).values
    c_candidates=np.setdiff1d(bigg_e,list(model.medium.keys())) # model.medium should be same for all models. This is initialized by the -i tag in carveme? 
    return c_candidates

def set_force_uptake(model, metabolite_id):
    # Copied from Ahmed's code. 
    metabolite_name = bigg.at[metabolite_id+'_e','name']
    ## Adding the external metabolite to the model
    metabolite_cobra = cobra.Metabolite(metabolite_id + '_e', name = metabolite_name, compartment= 'C_e')
    try:
        x = model.metabolites.get_by_id(metabolite_id + '_e')
    except KeyError:
        model.add_metabolites(metabolite_cobra)

    ## Adding periplasmic metabolite
    metabolite_p_cobra = cobra.Metabolite(metabolite_id + '_p', name = metabolite_name, compartment= 'C_p')
    try:
        x = model.metabolites.get_by_id(metabolite_id + '_p')
    except KeyError:
        model.add_metabolites(metabolite_p_cobra)  

    ## Adding cytoplasmic metabolite
    metabolite_c_cobra = cobra.Metabolite(metabolite_id + '_c', name = metabolite_name, compartment= 'C_c')
    try:
        x = model.metabolites.get_by_id(metabolite_id + '_c')
    except KeyError:
        model.add_metabolites(metabolite_c_cobra)   
        
    ## Adding the exchange reaction
    try:
        x = model.reactions.get_by_id("EX_" + metabolite_id + "_e")
    except KeyError: 
        reaction_EX = cobra.Reaction("EX_" + metabolite_id + "_e")
        reaction_EX.name = metabolite_name + ' exchange'
        reaction_EX.lower_bound = -10
        reaction_EX.upper_bound = 1000
        reaction_EX.add_metabolites({metabolite_cobra : -1})
        model.add_reactions([reaction_EX])

    ## Adding the diffusion reaction
    try:
        x = model.reactions.get_by_id(metabolite_id.upper() + "tex")
    except KeyError:     
        reaction_diff = cobra.Reaction(metabolite_id.upper() + "tex")
        reaction_diff.name = metabolite_name + ' diffusion'
        reaction_diff.lower_bound = -1000
        reaction_diff.upper_bound = 1000
        reaction_diff.add_metabolites({metabolite_cobra: -1, metabolite_p_cobra: 1})
        model.add_reactions([reaction_diff])
    
    ## Adding the ATP Uptake reaction
    reaction_uptake = cobra.Reaction(metabolite_id + "_UP")
    reaction_uptake.name = metabolite_name + ' uptake'
    reaction_uptake.lower_bound = -10
    reaction_uptake.upper_bound = 1000
    reaction_uptake.add_metabolites({metabolite_p_cobra: -1, model.metabolites.get_by_id('atp_c') :-1, metabolite_c_cobra: 1,  model.metabolites.get_by_id('adp_c'): 1, model.metabolites.get_by_id('pi_c'): 1})
    model.add_reactions([reaction_uptake])

def add_external_metabolite(model, metabolite_id):
    metabolite_name = bigg.at[metabolite_id+'_e','name']
    metabolite_cobra = cobra.Metabolite(metabolite_id, name = metabolite_name, compartment= 'C_e')
    try:
        x = model.metabolites.get_by_id(metabolite_id+'_e')
    except KeyError:
        model.add_metabolites(metabolite_cobra)

    try:
        x = x = model.reactions.get_by_id("EX_" + metabolite_id + "_e")
    except KeyError:
        reaction_EX = cobra.Reaction("EX_" + metabolite_id+'_e')
        reaction_EX.name = metabolite_name + ' Exchange'
        reaction_EX.lower_bound = 0
        reaction_EX.upper_bound = 1000
        reaction_EX.add_metabolites({metabolite_cobra : -1})
        model.add_reactions([reaction_EX])

def initialize_medium(model, metabolite_ids, default_flux=10):
    for metabolite_id in metabolite_ids:
        add_external_metabolite(model, metabolite_id)
    
    model.medium={'EX_'+metabolite_id+'_e':default_flux for metabolite_id in metabolite_ids}

def run_fba(model, c_candidates,force_uptake=False):
    growths={c:np.nan for c in c_candidates}

    with model: # negative control
        growth=model.slim_optimize()
        growths['cfree']=growth

    for c in c_candidates:
        try:
            with model: # Ahmed's code is wrong here. He put with outside the loop. That forces uptakes of all bigg metabolites. The results might be same, since the lower bound is changed to zero at the end, but it might confuses the optimization.
                add_external_metabolite(model, c[3:-2])
                model.reactions.get_by_id(c).lower_bound=-10 # Turn on the candidate exchange (uptake) # TODO: is this the same thing as setting medium? 
                if force_uptake:
                    set_force_uptake(model, c[3:-2])
                growth=model.slim_optimize()
                growths[c]=growth
                model.reactions.get_by_id(c).lower_bound=0 # Turn off the candidate exchange (uptake). Although this shouldn't be necessary. 
        except KeyError:
            print('Error with',c)
            growths[c]=np.nan
    growths=pd.Series(growths)
    # growths=(growths>=growth_threshold).astype(int)
    return growths


@click.command()
@click.option("--force_uptake",is_flag=True)
@click.option("--output",required=True)
@click.option("--initial_medium",required=True)
@click.argument('model',nargs=1)
def main(force_uptake,output,initial_medium,model):
    base_medium=pd.read_csv(initial_medium,sep='\t',header=None, names=['media','description','bigg_id','metabolite_description'])['bigg_id'].values
    model=load_fba_model(model)
    initialize_medium(model, base_medium,default_flux=10)
    c_condidates=find_c_condidates(model)
    growths=run_fba(model,c_condidates,force_uptake=force_uptake)
    growths.to_csv(output,header=False)

if __name__=="__main__":
    main()
