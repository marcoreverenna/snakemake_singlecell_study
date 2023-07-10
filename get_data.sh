#!/bin/bash

# Create directories
mkdir reference map

# Download files
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr1.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr2.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr3.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr4.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr5.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr6.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr7.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr8.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr9.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr10.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr11.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr12.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr13.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr14.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr15.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr16.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr17.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr18.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference/https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr19.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr20.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr21.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P reference https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr22.1kg.phase3.v5a.vcf.gz
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P map https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar

# Move downloaded data to the proper directories
mv chr* reference
mv plink.GRCh37.map.zip map

unzip map/plink.GRCh37.map.zip
mv *.map map

rm README.txt
rm plink.README.txt
