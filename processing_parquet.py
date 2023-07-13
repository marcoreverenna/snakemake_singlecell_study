from src import snparray_to_vcf_mod as snplib
from seqseek import Chromosome, BUILD37
import os
import pandas as pd

list_families = ['PGD036']#, 'PGD043']
dir = 'data/raw/'

if __name__=='__main__':
    for family in list_families:
        vcf_ref= snplib.import_vcf_reference()
        vcf_ref['REFALT_DBSNP']=vcf_ref['REF']+vcf_ref['ALT_1']
        vcf_ref=vcf_ref.rename(columns={'ID':'Name'})[['REFALT_DBSNP','Name']]
        to_refalt=snplib.load_manifest('data/raw/HumanCytoSNP-12v2-1_hg19.bpm')
        df = pd.read_parquet(f'data/raw/{family}_stage2_nofilter_hg19.df.parquet',engine='pyarrow')
        df=df.set_index('Name').join(to_refalt.set_index('Name'))
        df.rename(columns={'variant/REFALT':'REFALT_MANIFEST'},inplace=True)
        df=df.join(vcf_ref.set_index('Name'))
        df['gtype_vcf']=df.apply(lambda x: snplib.map_ab_to_bin(x['REFALT_MANIFEST'],x['REFALT_DBSNP'],x['gtype'],x['Chr'],x['Position'],BUILD37),axis=1)
        df['gtype_reconstructed']=df.apply(lambda x: snplib.map_bin_to_ab(x['REFALT_MANIFEST'],x['REFALT_DBSNP'],x['gtype_vcf'],x['Chr'],x['Position'],BUILD37),axis=1)
        parq = df.reset_index().copy()
        parq.to_parquet(os.path.join(dir, f'{family}_stage2_nofilter_hg19_processed.parquet'), engine='pyarrow')
