#PBS -N alig.223222.se
#PBS -q default
#PBS -l nodes=1:ppn=8,walltime=999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load hisat2/2.0.4


INDEX_PATH=/LUSTRE/usuario/mtorres/pepper/idx

for i in $(ls out* | grep _1 | sed 's/_1.*//g')
do

hisat2 -p 8 --dta -q $INDEX_PATH/pepper.idx -U ${i}_1.fastq -S ${i}.sam

done
