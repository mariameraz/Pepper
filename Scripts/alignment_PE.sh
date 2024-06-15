#PBS -N alig.235215
#PBS -q default
#PBS -l nodes=1:ppn=8,walltime=999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load hisat2/2.0.4 samtools/1.9


INDEX_PATH=/LUSTRE/usuario/mtorres/pepper/idx

for i in $(ls out* | grep _1 | sed 's/_1.*//g')
do

hisat2 -p 8 --dta -q $INDEX_PATH/pepper.idx -1 ${i}_1.fastq -2 ${i}_2.fastq -S ${i}.sam

done


for i in $(ls out* | sed 's/.sam//g')
do

#sam to bam file
samtools view -bS ${i}.sam > ${i}.bam

#sorting bam files
samtools sort ${i}.bam -o ${i}.sorted.bam

done
