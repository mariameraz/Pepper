#PBS -N fasta.235215
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=99999:999:00,mem=36gb,vmem=36gb
#PBS -V

cd $PBS_O_WORKDIR

module load sratoolkit/2.8.2

fastq-dump --split-files SRR*
