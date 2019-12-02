cd /mnt/volume/jiajia/bin\
./pbsim --data-type CCS --depth 30 --prefix p30.SNP_INDEL_chr20 --model_qc model_qc_ccs /mnt/volume/jiajia/data/ref/SNP_INDEL_chr20_hs37d5.fa\
cat p30.SNP_INDEL_chr20* > p30.SNP_INDEL_chr20.fastq\
bwa mem -t 10 -x pacbio -R '@RG\\tID:88888888\\tPL:PACBIO\\tSM:HG666\\tPM:SEQUEL' /mnt/volume/jiajia/data/ref/chr20_hs37d5.fa p30.SNP_INDEL_chr20.fastq -o p30.SNP_INDEL_chr20.sam\
samtools view -@11 -bS -o p30.SNP_INDEL_chr20.bam p30.SNP_INDEL_chr20.sam\
samtools sort -@11 -o p30.SNP_INDEL_chr20.sorted.bam p30.SNP_INDEL_chr20.bam\
samtools index -b -@11 p30.SNP_INDEL_chr20.sorted.bam\
cp p30.SNP_INDEL_chr20.sorted.bam* /mnt/volume/jiajia/pacbio_simulate/input/data\
cd /mnt/volume/jiajia/pacbio_simulate/input/data\
