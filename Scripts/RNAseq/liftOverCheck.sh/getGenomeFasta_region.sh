~/seqtk/seqtk \
subseq \
/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.interleaved.fa \
~/region.D84a.bed > ~/region.D84a.fasta


cp ~/region.D84a.bed ~/chr2217.D84a.bed


~/seqtk/seqtk \
subseq \
/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.interleaved.fa \
~/chr2217.D84a.bed > ~/chr2217.D84a.fasta

grep -A1 ">Scaffold_2217_HRSCAF_2652" totalHiCwithallbestgapclosed.fa > ~/chr2217.D84a.fasta

scp aob2x@rivanna.hpc.virginia.edu:~/region.D84a.fasta ~/.
