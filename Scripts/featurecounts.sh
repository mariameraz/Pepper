#PBS -N fc.235215
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=999:999:00,mem=30gb,vmem=30gb
#PBS -V

cd $PBS_O_WORKDIR

module load subread/1.5.1

PATH=/LUSTRE/usuario/mtorres/pepper
FC=/data/software/subread-1.5.1-source/bin

$FC/featureCounts -T 4 -t exon -g gene_id -a $PATH/GCA_000512255.2_ASM51225v2_genomic.gtf -o counts.235215.txt *.bam
