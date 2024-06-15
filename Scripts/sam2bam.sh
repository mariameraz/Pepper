#PBS -N bam.505972
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load samtools/1.9

for i in $(ls out* | sed 's/.sam//g')
do

#sam to bam file
samtools view -bS ${i}.sam > ${i}.bam

#sorting bam files
samtools sort ${i}.bam -o ${i}.sorted.bam

done
