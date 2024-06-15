#PBS -N fastp.223222
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=99999:999:00,mem=36gb,vmem=36gb
#PBS -V

cd $PBS_O_WORKDIR

module load fastp/0.20

#SE
for i in $(ls | grep _1 | sed 's/_1\.fastq//g')
do
fastp -i ${i}_1.fastq -o out.${i}_1.fastq -q 30
done

