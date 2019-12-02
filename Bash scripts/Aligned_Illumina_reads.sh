
cd /mnt/volume/jiajia/bin\
./art_illumina -ss HS25 -sam -i /mnt/volume/jiajia/data/ref/SNP_INDEL_chr20_hs37d5.fa -l 150 -f 50 -o i.SNP_INDEL_chr20\
rm i.SNP_INDEL_chr20.sam\
bwa mem -t 10 -R '@RG\\tID:NIST_150BP\\tSM:HG666' /mnt/volume/jiajia/data/ref/chr20_hs37d5.fa i.SNP_INDEL_chr20.fq > i.SNP_INDEL_chr20.sam\
samtools view -@11 -bS -o i.SNP_INDEL_chr20.bam i.SNP_INDEL_chr20.sam\
samtools sort -@11 -o i.SNP_INDEL_chr20.sorted.bam i.SNP_INDEL_chr20.bam\
samtools index -b -@11 i.SNP_INDEL_chr20.sorted.bam\
mv i.SNP_INDEL_chr20.sorted.bam* /mnt/volume/jiajia/illumina_simulate/input/data/.\
cd /mnt/volume/jiajia/illumina_simulate/input/data
