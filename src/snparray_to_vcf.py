#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for conversion of VCF-coded genotypes to Illumina A/B and vice versa

__author__ = Ivan Vogel, Lishan Cai
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Ivan Vogel
__email__ = ivogel@sund.ku.dk
__status__ = Dev
"""

import pandas as pd
import allel
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype
from seqseek import Chromosome, BUILD37
import sys


complement={
    'CT': 'AG'
}


complement_nucl={
    'A': 'T',
    'G': 'C',
    'C': 'G',
    'T': 'A'
}

#def map_bin_to_ab_2(refalt,gt,chr,pos,assembly = BUILD37):
def map_bin_to_ab(refalt_manifest, refalt_dbsnp, gt_bin, chr, pos,assembly):
    """ Function mapping genotype call in nucleotides to Illumina AB genotypes

    Args:
        refalt (str): reference and alternate allele (concatenated, i.e. AT, GC etc.)
        gt (str): genotype call in nucleotide format (i.e. AT)
        chr (str): _description_ ()
        pos (int): _description_
        assembly (_type_, optional): version of genome assembly. Defaults to BUILD37.

    Returns:
        str: genotype call in Illumina format (AA,BB,AB or NC)
    """
    def revert(r):
        if r=='AA': return 'BB'
        elif r=='BB': return 'AA'
        elif r=='AB': return 'AB'
        else: return 'NC'
        
    def complement(r):
        if r=='A': return 'T'
        elif r=='T': return 'A'
        elif r=='G': return 'C'
        elif r=='C': return 'G'
        else: return '-'
        
    retval='NC'
    # translate to refalt nucleotides
    if chr=='XY' or gt_bin[0]==-1 or gt_bin[1]==-1 or pd.isna(refalt_dbsnp): retval='NC'
    
    else:
        gt_nucl=refalt_dbsnp[gt_bin[0]] + refalt_dbsnp[gt_bin[1]]
        if gt_nucl[0] not in refalt_manifest or gt_nucl[1] not in refalt_manifest:
            old_val=gt_nucl
            gt_nucl=complement(old_val[0])+ complement(old_val[1])
            
        if gt_nucl == 'AT'  or gt_nucl == 'TA' or gt_nucl == 'AC' or gt_nucl == 'CA' or gt_nucl == 'AG' or gt_nucl == 'GA' or gt_nucl == 'TC' or gt_nucl == 'CT' or gt_nucl == 'TG' or gt_nucl == 'GT' or gt_nucl == 'CG' or gt_nucl == 'GC' or gt_nucl[0] != gt_nucl[1]: retval='AB'
        elif refalt_manifest == 'AT' or refalt_manifest =='TA' or refalt_manifest == 'CG' or refalt_manifest == 'GC':
            # these calls are all homozygous, so gt_nucl[0]==gt_nucl[1]
            #
            assert(gt_nucl[0]==gt_nucl[1])
            
            planking1 = Chromosome(chr,assembly).sequence(pos-100, pos-1)[::-1]
            planking2 = Chromosome(chr,assembly).sequence(pos, pos+100)
            strand = walking(list(map(sum_1,planking1,planking2)))
            if strand == 'Top' and 'A' in gt_nucl: retval='AA'
            elif strand == 'Bot' and 'A' in gt_nucl: retval='BB'
            elif strand == 'Top' and 'T' in gt_nucl: retval='BB'
            elif strand == 'Bot' and 'T' in gt_nucl: retval='AA'
            elif strand == 'Top' and 'C' in gt_nucl: retval='AA'
            elif strand == 'Bot' and 'C' in gt_nucl: retval='BB'
            elif strand == 'Top' and 'G' in gt_nucl: retval='BB'
            elif strand == 'Bot' and 'G' in gt_nucl: retval='AA'	
            else: retval='NC'
            
        elif gt_nucl == 'AA' or gt_nucl == 'TT': retval='AA'
        elif gt_nucl == 'CC' or gt_nucl == 'GG': retval='BB'            
        else: retval='NC'
        
    return retval 

def walking(seq):
    """ Sequence walking to determine TOP/BOT strand
    The method iterates list of alleles around the SNP of interest to determine the whether the SNP is on TOP or BOT strand
    Documentation here: https://www.illumina.com/documents/products/technotes/technote_topbot.pdf

    Args:
        seq (list): list of nucl[n-k] + nucl[n+k] where n is the locus of interest and k is offset

    Returns:
        str: returns Top or Bot
    """
    for i in seq:
        if i == 'AG' or i == 'AC' or i == 'TC' or i == 'TG': return'Top'
        elif i == 'GA' or i == 'GT' or i == 'CA' or i == 'CT': return'Bot'
        elif i =='multi': return 'Top'
        else: continue
    return seq


def sum_1(a,b):
    """Concatenate two nucleotides

    Args:
        a (str): first nucleotide
        b (str): second nucleotide

    Returns:
        str: concatenation of two nucleotides
    """
    return a+b 


def load_manifest(manifest_path='../snakemake_workflow_snp_array_data/data/HumanCytoSNP-12v2-1_hg19.bpm'):
    """ Generate dictionary of variants from Illumina manifest file

    Args:
        manifest_path (str, optional): _description_. Defaults to '/data/gqc954/TMP/CLUSTER_MANIFEST/HumanCytoSNP-12v2-1_hg19.bpm'.

    Returns:
        dic: dictionary where key is SNP ID and value is reference and variant allele
    """
    
    try:
        manifest=BeadPoolManifest(manifest_path)
    except:
        raise Exception('load_manifest: Error loading manifest')

    #df_manifest=pd.DataFrame.from_dict({'variant/REFALT':manifest.snps,'Name':manifest.names,'Strand': manifest.})
    df_manifest = pd.DataFrame.from_dict({'variant/REFALT': manifest.snps,
                                          'Name': manifest.names,
                                          'source_strands': manifest.source_strands,
                                          'ref_strands': manifest.ref_strands
                                          })
    
    df_manifest['variant/REFALT'] = df_manifest['variant/REFALT'].str.replace(r'\[|\]|\/', '')
    # according to BeadPoolManifest module from IlluminaBeadArray, https://github.com/Illumina/BeadArrayFiles
    #array_to_swap=df_manifest.loc[df_manifest.source_strands==2,'variant/REFALT']
    #df_manifest.loc[df_manifest.source_strands==2,'variant/REFALT']=array_to_swap.apply(lambda x: x[1]+x[0])
    return df_manifest
    #return dict(zip(df_manifest['Name'],df_manifest['variant/REFALT']))


def map_ab_to_bin(refalt_manifest, refalt_dbsnp, gt_illumina, chrom, pos, assembly):
    """ Map Illumina AB call to vcf format

    Args:
        ref (str): reference allele (nucleotide)
        alt (str): alternative allele (nucleotide)
        gt_illumina (str): illumina call (AA|BB|AB|NC)
        chrom (str): chromosome
        pos (int): position of the SNP

    Returns:
        list: call where 0 is reference allele, 1 is alternative allel, i.e. [0,0]
    """        
    # this is default for het calls
    vcf_gtypes = {'major':(0,0),
                  'minor':(1,1),
                  'het':(0,1),
                  'nc':(-1,-1)}
    
    retval = vcf_gtypes['het']
    retnucl = '--'
    
    if chrom == 'XY' or gt_illumina == 'NC' or pd.isna(refalt_dbsnp):
        retnucl = '--'
        retval = vcf_gtypes['nc']
    elif gt_illumina != 'AB':
        transl = {nucl: refalt_dbsnp.rfind(nucl) for nucl in ['A','T','G','C']}

        # is homozygous and amibiguous 
        if refalt_manifest == 'AT' or refalt_manifest == 'TA' or refalt_manifest == 'CG' or refalt_manifest == 'GC':
            # Implemented according to https://www.illumina.com/documents/products/technotes/technote_topbot.pdf
            planking1 = Chromosome(chrom,assembly).sequence(pos-21, pos-1)[::-1]
            planking2 = Chromosome(chrom,assembly).sequence(pos, pos+20)
            strand = walking(list(map(sum_1, planking1, planking2)))
            #assert(strand.lower()==source_strand.lower())
            
            # C or A
            if strand == 'Top' and gt_illumina == 'AA': 
                retnucl = 'CC' if 'C' in refalt_manifest else 'AA'
                #retval = [transl['C'],transl['C']] if transl['C']>= 0 else [transl['A'],transl['A']]
            # C or A
            elif strand == 'Bot' and gt_illumina =='BB':
                retnucl = 'CC' if 'C' in refalt_manifest else 'AA'
                #retval = [transl['C'],transl['C']] if transl['C']>= 0 else [transl['A'],transl['A']]
            # G or T
            elif strand == 'Top' and gt_illumina == 'BB':
                retnucl ='GG' if 'G' in refalt_manifest else 'TT'
                #retval = [transl['G'],transl['G']] if transl['G']>= 0 else [transl['T'],transl['T']]
            # G or T
            elif strand == 'Bot' and gt_illumina == 'AA':
                retnucl = 'GG' if 'G' in refalt_manifest else 'TT'
                #retval = [transl['G'],transl['G']] if transl['G']>= 0 else [transl['T'],transl['T']]
                #retval=vcf_gtypes['minor']              #'T' in gt: return 'AA'
            #elif strand == 'Top' and gt_illumina=='AA':retval=vcf_gtypes['major']              #'C' in gt: return 'AA'
            #elif strand == 'Bot' and gt_illumina=='BB':retval=vcf_gtypes['major']              #'C' in gt: return 'BB'
            #elif strand == 'Top' and gt_illumina=='BB':retval=vcf_gtypes['minor']              #'G' in gt: return 'BB'
            #elif strand == 'Bot' and gt_illumina=='AA':retval=vcf_gtypes['minor']              #'G' in gt: return 'AA'	
                
        elif (gt_illumina == 'AA'): #and ilmn_strand=='TOP') or (gt_illumina=='BB' and ilmn_strand=='BOT'): #retval=vcf_gtypes['major']
            retnucl = refalt_manifest[0] + refalt_manifest[0]
        elif (gt_illumina == 'BB'): #and ilmn_strand=='TOP') or (gt_illumina=='AA' and ilmn_strand=='BOT'): #retval=vcf_gtypes['minor']
            retnucl = refalt_manifest[1] + refalt_manifest[1]
        else: #retval=vcf_gtypes['nc']
            retnucl = '--'
        # translate retnucl
        if '-' not in retnucl and transl[retnucl[0]] >= 0:
                retval = transl[retnucl[0]], transl[retnucl[1]]
        elif '-' not in retnucl:
            # complement
            retval = transl[complement_nucl[retnucl[0]]], transl[complement_nucl[retnucl[0]]]
        else:
            retval = vcf_gtypes['nc']
    
    else:
        retval = vcf_gtypes['het']
    return retval

    
def get_strand_from_refalt(chrom, pos, refalt, assembly = BUILD37):
    if chrom == 'XY':
        return 'N'
    #transl={nucl: refalt.rfind(nucl) for nucl in ['A','T','G','C']}
    # Implemented according to https://www.illumina.com/documents/products/technotes/technote_topbot.pdf
    planking1 = Chromosome(chrom, assembly).sequence(pos-21, pos-1)[::-1]
    planking2 = Chromosome(chrom, assembly).sequence(pos, pos+20)
    strand = walking(list(map(sum_1, planking1, planking2)))
    if isinstance(strand, list):
        return 'N'
    else:
        return strand
    
    

def map_to_nucleotide(gt, ref, alt):
    """ This is an intermediate procedure for converting binary vcf to Illumina calls 
    The function converts binary vcf notation to nucleotide using ref and alt allele as values
    In case it is a heterozygous call, returns directly AB and in case the call cannot be resolved, it returns NC

    Args:
        gt (_type_): call in binary notation in an array, i.e. [0,0]
        ref (_type_): reference nucleotide 
        alt (_type_): alternative nucleotide

    Returns:
        _type_: _description_
    """
    gt_new = []
    if gt[0] != gt[1] and gt[0] >= 0 and gt[1] >= 0: return 'AB'
    for i in range(len(gt)):
        # 0 means it is reference allele
        if gt[i] == 0:
            gt_new.append(ref)
        # less than zero means there is no allele   
        elif gt[i] < 0:
            gt_new.append('')
        # it is an alternative allele, starting with 1 (0 is reference)
        # therefore, in case there is multiple alt alleles, index them with i-1
        elif gt[i] != 0:
            alt = list(alt)
            try:
                gt_new.append(alt[gt[i]-1])
            except:
                return 'NC'
    return ''.join(gt_new)


def load_params():
    """ Routine that loads the params from the command line
    Currently, there is only input parquet file as input parameters
    TODO: use argparse if we want to extend the number of parameters
    """
    params = dict()
    if len(sys.argv) == 2:
        params['parquet']
    elif len(sys.argv) == 1:
        params['parquet'] = 'PGD036_stage2_nofilter_hg19.df.parquet'
    else:
        raise ValueError('Command line error - wrong number of arguments')


def import_vcf_reference():
    """_summary_

    Returns:
       pandas.Dataframe: vcf file in dataframe format
    """
    return allel.vcf_to_dataframe('../snakemake_workflow_snp_array_data/data/GM12878_gDNA_3.vcf.gz')


def extract_dbsnp_refalt(x):
    pass
    

if __name__=='__main__':    
    # load command line parameters
    ## params=load_params()
    
    # load csv manifest
    #manifest_csv=pd.read_csv('metadata/HumanCytoSNP-12v2-1_L.csv',skiprows=7).set_index('Name')
    #manifest_csv[['REF','ALT']]=manifest_csv['TopGenomicSeq'].str.extract('\[([ATGC])/([ATGC])\]',expand=True)
    #manifest_csv['REFALT_MANIFEST']=manifest_csv['REF']+manifest_csv['ALT']
    
    # import an example of vcf file to load the ref and alt allele as used in dbSNP database
    vcf_ref=import_vcf_reference()
    vcf_ref['REFALT_DBSNP']=vcf_ref['REF']+vcf_ref['ALT_1']
    vcf_ref=vcf_ref.rename(columns={'ID':'Name'})[['REFALT_DBSNP','Name']]
        
    # load manifest
    to_refalt = load_manifest('../snakemake_workflow_snp_array_data/data/HumanCytoSNP-12v2-1_hg19.bpm')    
    
    # load parquet file
    # use params['parquet'] if loading from command line argument
    df = pd.read_parquet('data/PGD036_stage2_nofilter_hg19.df.parquet',engine='pyarrow')

    # map SNP ID to ref and alt alleles
    #df['refalt']=df.Name.map(to_refalt)
    
    
    df=df.set_index('Name').join(to_refalt.set_index('Name'))
    df.rename(columns={'variant/REFALT':'REFALT_MANIFEST'},inplace=True)
    df=df.join(vcf_ref.set_index('Name'))
    
    
    # convert to binary vcf format
    df['gtype_vcf']=df.apply(lambda x: map_ab_to_bin(x['REFALT_MANIFEST'],x['REFALT_DBSNP'],x['gtype'],x['Chr'],x['Position'],BUILD37),axis=1)

   
    # # check if the conversion was correct by converting back to nucleotides ...
    #df['gtype_nucleotide']=df.apply(lambda x: map_to_nucleotide(x['gtype_vcf'],x['refalt'][0],x['refalt'][1]),axis=1)
    
    # ... and then to Illumina calls
    df['gtype_reconstructed']=df.apply(lambda x: map_bin_to_ab(x['REFALT_MANIFEST'],x['REFALT_DBSNP'],x['gtype_vcf'],x['Chr'],x['Position'],BUILD37),axis=1)
    
    
    # check concordance
    print(df[['gtype','gtype_reconstructed']].value_counts())
    print(df.gtype_vcf.value_counts())
