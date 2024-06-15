#PBS -N fasta.235215
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=99999:999:00,mem=36gb,vmem=36gb
#PBS -V

cd $PBS_O_WORKDIR

module load sratoolkit/2.8.2

#Single end
#for i in $(more ../SraRunInfo.csv | sed 's/,/\t/g' | cut -f 1,16 | grep -i single | cut -f 1)
#do
#fastq-dump ../$i
#done

#Paired end
#for i in $(more ../SraRunInfo.csv | sed 's/,/\t/g' | cut -f 1,16 | grep -i paired | cut -f 1)
#do 
#fastq-dump --split-files ../$i
#done


fastq-dump --split-files SRR*
