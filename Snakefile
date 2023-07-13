import os
import pandas as pd
from src import snparray_to_vcf as snplib
from seqseek import Chromosome, BUILD37

#configfile: "config.yaml"

CHROMOSOME_STR = [str(chrom) for chrom in list(range(1,23))]

rule all:
    input:
        'data/processed/PGD036_nomissing_sorted.parquet',
        expand('data/processed/PGD036_validated_{validation}.parquet', validation = ['correct', 'incorrect']),
        expand('data/processed/PGD036_chr{chrom}_correct.parquet', chrom=CHROMOSOME_STR),
        expand('data/processed/PGD036_chr{chrom}_incorrect.parquet', chrom=CHROMOSOME_STR),
        expand('data/processed/PGD036_chr{chrom}_correct.vcf', chrom=CHROMOSOME_STR),
        expand('data/processed/PGD036_chr{chrom}_correct_imputed.vcf.gz', chrom = CHROMOSOME_STR),
        expand('data/processed/PGD036_chr{chrom}_correct_imputed_flt.vcf.gz', chrom=CHROMOSOME_STR),
        expand('data/processed/PGD036_chr{chrom}_incorrect_imputed_flt.parquet', chrom=CHROMOSOME_STR),
        expand('data/processed/PGD036_chr{chrom}_incorrect_imputed_flt_gtimp.parquet', chrom=CHROMOSOME_STR),
        expand('data/processed/PGD036_chr{chrom}_correct_merged_gtimp.parquet', chrom=CHROMOSOME_STR),
        'data/final/PGD036_all_chroms_merged_imputed.parquet',
        'data/final/PGD036_stage2_hg19_imputation_filtered.parquet'

rule cleaning:
    """
    This rule performs a preliminary cleaning of the parquet file

    - Step 1 = removing NC
    - Step 2 = considering only autosomes
    - Step 3 = sorting values
    - Step 4 = resetting the index
    """
    input:
        parquet = 'data/raw/PGD036_stage2_nofilter_hg19_processed.parquet'
    output:
        parquet_sorted = 'data/processed/PGD036_nomissing_sorted.parquet'
    run:
        parquet = pd.read_parquet(input.parquet, engine='pyarrow')
        parquet = parquet[~parquet[['gtype_reconstructed', 'mother_gtype', 'father_gtype']].apply(lambda x: x.str.contains('NC')).any(axis=1)]
        parquet = parquet[parquet['Chr'].isin(CHROMOSOME_STR)]
        parquet['Chr'] = parquet['Chr'].astype(int)
        parquet.sort_values(by='Chr', ascending=True, inplace=True)
        parquet.reset_index(inplace=True, drop=True)
        parquet.to_parquet(output.parquet_sorted, engine='pyarrow')
        return parquet


rule check_genotypes:
    """
    This rule performs an embryo check by considering parental information
    """
    input:
        #'data/processed/PGD036_nomissing_sorted.parquet'
        parquet = rules.cleaning.output
    output:
        expand('data/processed/PGD036_validated_{validation}.parquet', validation = ['correct', 'incorrect'])
    run:
        parquet = pd.read_parquet(input.parquet, engine='pyarrow')
        parquet['validation'] = ['correct' if \
            (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'AA' and row['gtype_reconstructed'] == 'AA') or \
            (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AA') or \
            (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'AA' and row['father_gtype'] == 'BB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AA' and row['gtype_reconstructed'] == 'AA') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AA' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AA') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'BB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'BB' and row['gtype_reconstructed'] == 'BB') or \
            (row['mother_gtype'] == 'AB' and row['father_gtype'] == 'BB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'AA' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'AB') or \
            (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'AB' and row['gtype_reconstructed'] == 'BB') or \
            (row['mother_gtype'] == 'BB' and row['father_gtype'] == 'BB' and row['gtype_reconstructed'] == 'BB') \
            else 'incorrect' for _, row in parquet.iterrows()]
        parquet_correct = parquet[parquet['validation'] == 'correct']
        parquet_incorrect = parquet[parquet['validation'] == 'incorrect']
        parquet_correct.to_parquet(output[0], engine='pyarrow')
        parquet_incorrect.to_parquet(output[1], engine='pyarrow')
        return parquet_correct, parquet_incorrect


rule filt_autosomes_correct:
    """
    This rule creates a file for each chromosome containing the correct embryonic information.
    """
    input:
        correct_version = 'data/processed/PGD036_validated_correct.parquet'
    output:
        expand('data/processed/PGD036_chr{chrom}_correct.parquet', chrom=CHROMOSOME_STR)
    run:
        for chrom in list(range(1,23)):
            parquet = pd.read_parquet(input.correct_version, engine='pyarrow')
            parquet_chr = parquet[parquet['Chr'] == int(chrom)]
            # remember to specify the index
            parquet_chr.to_parquet(output[chrom -1], engine='pyarrow')


rule filt_autosomes_incorrect:
    """
    This rule creates a file for each chromosome containing the incorrect embryonic information.
    """
    input:
        correct_version = 'data/processed/PGD036_validated_incorrect.parquet'
    output:
        expand('data/processed/PGD036_chr{chrom}_incorrect.parquet', chrom=CHROMOSOME_STR)
    run:
        for chrom in list(range(1,23)):
            parquet = pd.read_parquet(input.correct_version, engine='pyarrow')
            parquet_chr = parquet[parquet['Chr'] == int(chrom)]
            # remember to specify the index
            parquet_chr.to_parquet(output[chrom -1], engine='pyarrow')


rule vcf_conversion:
    """
    This rule defines all the variables necessary to create a VCF (variant calling format).
    """
    input:
        expand('data/processed/PGD036_chr{chrom}_correct.parquet', chrom=CHROMOSOME_STR)
    output:
        expand('data/processed/PGD036_chr{chrom}_correct.vcf', chrom=CHROMOSOME_STR)
    run:
        for chrom in list(range(1,23)):
            # define the paths
            path_parquet_file = input[chrom-1]
            path_vcf_file = output[chrom-1]
            # make the conversion in the same format
            parq_chr = pd.read_parquet(path_parquet_file, engine='pyarrow')
            parq_chr.rename(columns={'Chr': '#CHROM', 'Position': 'POS', 'Name': 'ID'}, inplace=True)
            parq_chr['QUAL'] = '.'
            parq_chr['FILTER'] = 'PASS'
            parq_chr['INFO'] = '.'
            parq_chr['FORMAT'] = 'GT'
            parq_chr['REF'] = parq_chr['REFALT_DBSNP'].apply(lambda x: x[0])
            parq_chr['ALT'] = parq_chr['REFALT_DBSNP'].apply(lambda x: x[1])
            parq_chr['SAMPLE'] = parq_chr['gtype_vcf'].apply(lambda x: f"{x[0]}/{x[1]}")
            vcf_sort = parq_chr.loc[:, ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']]
            # final adjustments
            vcf_sort.sort_values('POS', inplace=True)
            vcf_sort.reset_index(drop=True, inplace=True)

            vcf_sort.to_csv(path_vcf_file, sep='\t', index=False)

rule beagle_gt:
    """
    This rule take into account 3 different inputs to run the beagle algorithm.

    Below the files required to perform the imputation:
    - File1 = VCF file (sample file)
    - File2 = map file, containing all the positions referring to HG19/GRCh37
    - File3 = reference file, containing the SNPs from the population of individuals

    Please consider that in this preliminary version all the SNPs from Reference panel has been imputed.
    In this way we could have even the possibility to increment the density of our data.
    By the way, it is possible setting a parameter and specificy all the positions to exclude (--exclude-markers)
    """
    input:
        gtfile = 'data/processed/PGD036_chr{chrom}_correct.vcf',
        mapfile = 'map/plink.chr{chrom}.GRCh37.map',
        reffile = 'reference/chr{chrom}.1kg.phase3.v5a.vcf.gz'
    output:
        outfile = 'data/processed/PGD036_chr{chrom}_correct_imputed.vcf.gz'
    params:
        # this step is required as beagle command automatically doesn't need the extension but the output file is compressed in gz
        # in this specific case it is useful add params parameter
        outname = lambda wildcards, output: output.outfile.replace('.vcf.gz', '')
    shell:
        """
        echo 'Imputation of PGD036_chr{wildcards.chrom}_correct.vcf started...'
        java -jar beagle.22Jul22.46e.jar gt={input.gtfile} ref={input.reffile} out={params.outname} map={input.mapfile}
        """

rule filtering_post_imputation:
    """
    This rule filters the positions of incorrect embryos on imputed VCF, at the same time filter the positions on parquet file which beagle was not able to impute.
    Beagle is not always able to impute all the SNPs assayed on SNP-array data.

    It does the following:
    - Step 1 = dropping the duplicates
    - Step 2 = creating a list of SNP positions on bulk file and filter them on VCF
    - Step 3 = creating a list of SNP positions on VCF filtered and filter parquet file
    """
    input:
        vcf_files = expand('data/processed/PGD036_chr{chrom}_correct_imputed.vcf.gz', chrom=CHROMOSOME_STR),
        parquet_files = expand('data/processed/PGD036_chr{chrom}_incorrect.parquet', chrom=CHROMOSOME_STR)
    output:
        vcf_filtered = expand('data/processed/PGD036_chr{chrom}_correct_imputed_flt.vcf.gz', chrom=CHROMOSOME_STR),
        parquet_filtered = expand('data/processed/PGD036_chr{chrom}_incorrect_imputed_flt.parquet', chrom=CHROMOSOME_STR)
    run:
        for chrom in range(1, 23):
            path_parquet_file = input.parquet_files[chrom - 1]
            path_vcf_file = input.vcf_files[chrom - 1]

            parq_chr = pd.read_parquet(path_parquet_file, engine='pyarrow')
            vcf_chr = pd.read_csv(path_vcf_file, header=8, sep='\t', dtype='object')
            # Step 1
            vcf_chr.POS.drop_duplicates(inplace=True)
            # Step 2
            list_inc = parq_chr.Position.astype(str).tolist()
            vcf_chr_flt = vcf_chr[vcf_chr['POS'].isin(list_inc)]
            # Step 3
            list_imp_flt = vcf_chr_flt.POS.astype(int).tolist()
            parq_chr_flt = parq_chr[parq_chr['Position'].isin(list_imp_flt)]

            vcf_chr_flt.to_csv(output.vcf_filtered[chrom - 1], sep='\t', index=False)
            parq_chr_flt.to_parquet(output.parquet_filtered[chrom - 1], engine='pyarrow')


rule conversion_gt_imp:
    """
    This rule converts the binary format 1/0 into genotype format AB considering parquet and VCF.
    """
    input:
        parquet_files = expand('data/processed/PGD036_chr{chrom}_incorrect_imputed_flt.parquet', chrom=CHROMOSOME_STR),
        vcf_files = expand('data/processed/PGD036_chr{chrom}_correct_imputed_flt.vcf.gz', chrom=CHROMOSOME_STR)
    output:
        parquet_final = expand('data/processed/PGD036_chr{chrom}_incorrect_imputed_flt_gtimp.parquet', chrom=CHROMOSOME_STR)
    params:
        valid_values = [(0, 1), (1, 0), (0, 0), (1, 1)]
    run:
        for chrom in list(range(1,23)):
            path_parquet_file = input.parquet_files[chrom - 1]
            path_vcf_file = input.vcf_files[chrom - 1]
            prq_chr = pd.read_parquet(path_parquet_file, engine='pyarrow')
            vcf_chr = pd.read_csv(path_vcf_file, header=0, sep='\t', dtype='object')

            # sort all the values considering the positions
            prq_df_sort = prq_chr.sort_values('Position', ascending=True).reset_index(drop=True)
            # copy the variable in parquet file from vcf
            prq_df_sort['SAMPLE'] = vcf_chr['SAMPLE']
            # take the new variable (SAMPLE) and split
            prq_df_sort.iloc[:,-1] = prq_df_sort.iloc[:,-1].apply(lambda x: x.split(':')[0])
            # add the brackets
            prq_df_sort['SAMPLE'] = prq_df_sort['SAMPLE'].apply(lambda x: [int(i) for i in x.split('|')])
            # convert the chromosome type
            prq_df_sort['Chr'] = prq_df_sort['Chr'].astype(str)
            # it is important to filter considering the parameters
            prq_df_sort['SAMPLE'] = prq_df_sort['SAMPLE'].apply(tuple)
            """
            #prq_df_sort['gtype_vcf'] = prq_df_sort['gtype_vcf'].apply(tuple)
            #prq_df_sort['SAMPLE'] = prq_df_sort['SAMPLE'].apply(tuple)
            """
            # considering the params object
            prq_df_sort = prq_df_sort[prq_df_sort['SAMPLE'].apply(lambda x: x in params.valid_values)]
            # converted considering unpashed even if beagle already has phased everything to do the imputation
            prq_df_sort['SAMPLE'] = prq_df_sort['SAMPLE'].apply(lambda x: (0, 1) if x == (1, 0) else x)
            prq_df_sort['sample_adapted'] = [([x[0], x[1]]) for x in prq_df_sort['SAMPLE']]
            # reconvert the chromosome tuple in string
            prq_df_sort['Chr'] = prq_df_sort['Chr'].astype(str)
            # converting in gtype imputed
            prq_df_sort['gtype_imputed'] = prq_df_sort.apply(
                                                            lambda x: snplib.map_bin_to_ab(
                                                                refalt_manifest = x['REFALT_MANIFEST'],
                                                                refalt_dbsnp = x['REFALT_DBSNP'],
                                                                gt_bin = x['sample_adapted'],
                                                                chr = x['Chr'],
                                                                pos = x['Position'],
                                                                assembly = BUILD37),
                                                                axis=1
                                                            )
            prq_df_sort.to_parquet(output.parquet_final[chrom - 1], engine='pyarrow')

# at this point we should have corrected the incorrect information, we can merge it with the correct one
# we need to filter the initial dataset with the same positions as the final ones
# then take exactly the initial one and run HMM, compare the error rate between case 1 and case 2

rule concat_parquets:
    """
    This rule concatenates the correct parquet file obtained by checking the parental information and the imputed information filtered.
    """
    input:
        parquet_files_gtype_imp = expand('data/processed/PGD036_chr{chrom}_incorrect_imputed_flt_gtimp.parquet', chrom=CHROMOSOME_STR),
        parquet_files_correct = expand('data/processed/PGD036_chr{chrom}_correct.parquet', chrom=CHROMOSOME_STR)
    output:
        parquet_merged = expand('data/processed/PGD036_chr{chrom}_correct_merged_gtimp.parquet', chrom=CHROMOSOME_STR)
    run:
        for chrom in list(range(1,23)):
            path_parquet_correct = input.parquet_files_correct[chrom-1]
            path_parquet_gtype_imp = input.parquet_files_gtype_imp[chrom-1]
            prq_chr_cor = pd.read_parquet(path_parquet_correct, engine='pyarrow')
            prq_chr_imp = pd.read_parquet(path_parquet_gtype_imp, engine='pyarrow')

            prq_chr_imp.drop(columns=['validation', 'SAMPLE', 'sample_adapted', 'gtype_imputed'], axis=1, inplace=True)
            prq_chr_cor.drop(columns=['validation'], axis=1, inplace=True)

            prq_chr_cor['Chr'] = prq_chr_cor['Chr'].astype(str)

            prq_chr_cor.reset_index(drop=True, inplace=True)
            prq_chr_imp.reset_index(drop=True, inplace=True)

            merged_df = pd.concat([prq_chr_cor, prq_chr_imp]).reset_index(drop=True)
            merged_df.to_parquet(output.parquet_merged[chrom - 1], engine='pyarrow')


rule concat_merged_parquets:
    """
    This rule concatenates all the parquet files containing both correct and imputed embryos genotypes.
    """
    input:
        parquet_merged = expand('data/processed/PGD036_chr{chrom}_correct_merged_gtimp.parquet', chrom=CHROMOSOME_STR)
    output:
        final_merged = 'data/final/PGD036_all_chroms_merged_imputed.parquet'
    run:
        chromosome_list = [str(chrom) for chrom in range(1, 23)]
        dfs = []

        for chrom in list(range(1,23)):
            parquet_file = input.parquet_merged[chrom-1]
            df = pd.read_parquet(parquet_file, engine='pyarrow')
            dfs.append(df)

        merged_df = pd.concat(dfs, ignore_index=True)
        merged_df.to_parquet(output.final_merged, engine='pyarrow')


rule same_positions_original:
    """
    This rule takes the same positions of the original parquet file.
    """
    input:
        parquet_merged_imputed = 'data/final/PGD036_all_chroms_merged_imputed.parquet',
        parquet_original = 'data/raw/PGD036_stage2_nofilter_hg19_processed.parquet'
    output:
        parquet_original_filtered = 'data/final/PGD036_stage2_hg19_imputation_filtered.parquet'
    run:
        merged_df = pd.read_parquet(input.parquet_merged_imputed, engine='pyarrow')
        original_df = pd.read_parquet(input.parquet_original, engine='pyarrow')
        filtered_df = original_df[original_df['Position'].isin(merged_df['Position'])]
        filtered_df.to_parquet(output.parquet_original_filtered, engine='pyarrow')


        #'data/PGD036_all_chroms_merged_imputed.parquet'
        #'data/PGD036_stage2_hg19_imputation_filtered.parquet'
