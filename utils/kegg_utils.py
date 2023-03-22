import numpy as np, pandas as pd
import os, re, requests
from tqdm import tqdm

class Kegg:
    @staticmethod
    def download_csv(query, labels=None, get_col=None,max_queries=10):
        if isinstance(query,(list,tuple,np.ndarray,pd.Series)):
            # query[0] is the rest api commmand, and query[1:] are the arguments
            cmd=query[0]
            out=[]
            for i, qs in tqdm(enumerate(np.split(query[1:],np.arange(len(query[1:]))[::max_queries])[1:])):
                # print(f"{cmd}/{'+'.join(qs)}")
                try:
                    res=Kegg.download_csv(f"{cmd}/{'+'.join(qs)}",labels=labels,get_col=get_col)
                    out.append(res)
                except pd.errors.EmptyDataError:
                    pass
            if out:    
                if get_col is None:
                    out=pd.concat([o for o in out if len(o)>0],axis=0,ignore_index=True)
                else:
                    out=np.concatenate(out)
            return out
        else:
            try:
                out=pd.read_csv(f"http://rest.kegg.jp/{query}",sep="\t",header=None)
                if get_col:
                    return out.iloc[:,get_col].values
                elif labels:
                    out.columns=labels
                return out
            except pd.errors.EmptyDataError:
                return []

    @staticmethod
    def download_genes_under_kos(kos,max_queries=10):
        if max_queries:
            out=[]
            for i,ko in enumerate(np.split(kos,np.arange(len(kos))[::max_queries])[1:]):
                res=Kegg.download_csv("link/genes/"+('+'.join(ko)),labels=['ko','genes'])
                out.append(res)
                print(f"Progress: {(i+1)*max_queries}/{len(kos)}")
            out=pd.concat(out,axis=0,ignore_index=True)
            return out
        else:
            return Kegg.download_csv("link/genes/"+('+'.join(kos)),labels=['ko','genes'])

        
    @staticmethod
    def get_entry(ids, entry, p=None,max_query=10,labels=None):
        success=[]
        fail=[]
        if p is not None:
            p=entry.upper()+r'[\s]+(.*)\n'
        for i,genes in enumerate(np.split(ids,np.arange(len(ids))[::max_query])[1:]):
            try:
                res=requests.get("http://rest.kegg.jp/get/%s" % ("+".join(genes)))
                if res:
                    texts=res.text.split(r"///")[:-1]
                    for gene_id, text in zip(genes,texts):
                        name=re.findall(p,text)
                        if len(name)!=1:
                            fail.append(gene_id)
                        else:
                            success.append([gene_id,name[0]])
                else:
                    fail.extend(genes)
                print("Successful: %d/%d" % (len(success),(i+1)*max_query))
            except Exception as e:
                print(e)
                break
        success=pd.DataFrame(success)
        if labels:
            success.columns=labels

        return success, fail
    
    def download_fasta(ko,max_queries=10,fa_type='aaseq',ff_out=None):
        
        if ff_out is None:
            ff_out=f'{ko}_{fa_type}.fa'

        if os.path.exists(ff_out):
            raise FileExistsError

        print(f'Downloading genes under {ko}')
        genes_all=Kegg.download_csv(f'link/genes/{ko}',get_col=1)

        for i,genes in enumerate(np.split(genes_all,np.arange(len(genes_all))[::max_queries])[1:]):
            #print(genes)
            res=requests.get("http://rest.kegg.jp/get/%s/%s" % ("+".join(genes), fa_type))
            if res:
                with open(ff_out,'a') as f:
                    f.write(res.text)
            else:
                print("Query failed: %s" % genes)
            print(f'Progress: {max_queries*(i+1)}/{len(genes_all)}')

        # check completeness
        from genomics_utils import IO
        print("Total genes: %d; downloaded: %d" % (len(genes_all),IO.read_fasta(ff_out,df=True).shape[0]))