#PBS -N fastp.223222
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=99999:999:00,mem=36gb,vmem=36gb
#PBS -V

cd $PBS_O_WORKDIR

module load fastp/0.20

#SE
for i in $(more Sra.RunInfo.csv | sed 's/,/\t/g' | cut -f 2,17 | grep SINGLE | cut -f 1 | sed 's/^.//g' | sed 's/.$//')
do
fastp -i ${i}_1.fastq -o out.${i}_1.fastq -q 30
done

#PE
for i in $(more Sra.RunInfo.csv | sed 's/,/\t/g' | cut -f 2,17 | grep PAIRED | cut -f 1 | sed 's/^.//g' | sed 's/.$//')
do
fastp -i ${i}_1.fastq -I ${i}_2.fastq -o out.${i}_1.fastq -O out.${i}_2.fastq -q 30
done
