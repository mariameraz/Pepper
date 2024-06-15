#PBS -N fastp.235215
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=99999:999:00,mem=36gb,vmem=36gb
#PBS -V

cd $PBS_O_WORKDIR

module load fastp/0.20

for i in $(ls | grep _1 | sed 's/_1\.fastq//g')
do
fastp -i ${i}_1.fastq -I ${i}_2.fastq -o out.${i}_1.fastq -O out.${i}_2.fastq -q 30 -h ${i}.html -j ${i}.json 
done
