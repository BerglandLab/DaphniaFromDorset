module load gcc

git clone https://github.com/lh3/minimap2
cd minimap2 && make


wget http://server7.wfleabase.org/genome/Daphnia_species_genomes/dplx20pubasm/daphplx16asm.fa.gz

wget \
-O /scratch/aob2x/daphnia_pafs/daphplx_gasm16ml.fa.gz \
http://arthropods.eugenes.org/EvidentialGene/daphnia/daphnia_pulex/daphnia_pulex_genes2017/chromosomes/daphplx16chr/daphplx_gasm16ml.fa.gz

gunzip /scratch/aob2x/daphnia_pafs/daphplx_gasm16ml.fa.gz

mkdir /scratch/aob2x/daphnia_pafs/

~/minimap2/minimap2 \
-x asm20 \
/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.interleaved.fa \
/scratch/aob2x/daphnia_pafs/daphplx_gasm16ml.fa > /scratch/aob2x/daphnia_pafs/D84A.daphplx_gasm16ml.paf

cut -f1-12 /scratch/aob2x/daphnia_pafs/D84A.daphplx_gasm16ml.paf > /scratch/aob2x/daphnia_pafs/D84A.daphplx_gasm16ml.taglesspaf

cut -f1-12 /scratch/aob2x/daphnia_pafs/D84A.daphplx16asm.paf > /scratch/aob2x/daphnia_pafs/D84A.daphplx16asm.tagless.paf
