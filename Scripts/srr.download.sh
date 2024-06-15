#PBS -N data.235215
#PBS -q default
#PBS -l nodes=1:ppn=1,walltime=99999:999:00,mem=8gb,vmem=8gb
#PBS -V

cd $PBS_O_WORKDIR

#module load R/3.6.0

#R CMD BATCH download.R

for i in $(more SraRunInfo.csv | sed 's/,/\t/g' | cut -f 1 | grep SRR)
do
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/$i/$i
done
