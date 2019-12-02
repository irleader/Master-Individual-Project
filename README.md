# Master-Individual-Project
This whole set of programs is able to run DeepVariant together with Hap.py and do a set of analysis with/without presence of truth VCF.

**Requirements:**

Unix-like operating system, Python 2.7, Docker

**Methods:**




1. generate artificial genomes:

    •	cut_fa.py: extract specific chromosomes from reference FASTA files.

    •	generate_chr20_SNP_fa.py: Generate artificial genome with SNP only for chromosome 20. For experiment 2 only.

    •	generate_artificial_SNP_vcf.py: Generate ground truth VCF for chr20 with reference and artificial genome. For experiment 2 only.

    •	generate_artificial_SNP_INDEL_vcf_fa.py: Generate artificial genome with SNP and INDEL for chr1-22 together with ground truth VCF for chr20. For experiments 3 and 4.



2. simulate Illumina and PacBio aligned reads:

    •	Aligned_Illumina_reads.sh: for simulating aligned Illumina reads from an artificial genome

    •	Aligned_PacBio_reads.sh: for simulating aligned PacBio reads from an artificial genome



3. run different bash scripts for different needs with DeepVariant and Hap.py:

    •	run_wgs_real_docker.sh: run DeepVariant and Hap.py with Illumina reads from human genome HG002

    •	run_pacbio_real_docker.sh: run DeepVariant and Hap.py with PacBio reads from human genome HG002

    •	run_illumina_simulate_docker.sh: run DeepVariant and Hap.py with Illumina reads from artificial genome

    •	run_pacbio_simulate_docker.sh: run DeepVariant and Hap.py with PacBio reads from artificial genome



4. with all Query VCFs ready, you can perform any of these analysis:

    •	Chr20_Truth_csv_real.py: Clean and Extract chromosome 20 variant calls only from Illumina, PacBio and Truth VCF. Store them into CSV for further process. For experiment 1 only.

    •	Chr20_Truth_csv_simulate.py: Clean and Extract chromosome 20 variant calls only from Illumina, PacBio and Truth VCF. Store them into CSV for further process. For experiments 2 and 3.

    •	Data_cleaning_real.py: Clean and extract all important info from Hap.py VCF file, remove duplicated variant calls. Convert Hap.py VCF files into csv files for further process. For experiment 1 only.

    •	Data_cleaning_simulate.py: Clean and extract all important info from Hap.py VCF file, remove duplicated variant calls. Convert Hap.py VCF files into csv files for further process. For experiments 2 and 3.

    •	Regions_real.py: Read CSV variant calls from Illumina, PacBio, and Truth. Divide all three files by intersections into 7 regions. For experiment 1 only.

    •	Regions_experiment2.py: Read CSV variant calls from Illumina, PacBio, and Truth. Divide all three files by intersections into 7 regions. Draw distribution graphs for variant information. Generate CSV for region C and E SNP. For experiment 2 only.

    •	Regions_experiment3.py: Read CSV variant calls from Illumina, PacBio, and Truth. Divide all three files by intersections into 7 regions. Draw distribution graphs for variant information. Generate CSV for region C and E SNP. For experiment 3 only.

    •	Repeats_real.py: Check proportion of repeats in each region, for experiment 1 only.

    •	Repeats_experiment2.py: Check proportion of repeats in each region, for experiment 2 only.

    •	Repeats_experiment3.py: Check proportion of repeats in each region, for experiment 3 only.

    •	SDs_real.py:  Check proportion of segmental duplications in each region, for experiment 1 only.

    •	SDs_experiment2.py: Check proportion of segmental duplications in each region, for experiment 2 only.

    •	SDs_experiment3.py: Check proportion of segmental duplications in each region, for experiment 3 only.


5. After all analysis, classifiers can differentiate True Positives and False Positives:

    •	Classify_C_E.py
