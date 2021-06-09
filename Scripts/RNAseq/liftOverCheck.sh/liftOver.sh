gene=Daphnia00786

chr=$( grep   ${gene} /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff | grep -E "mRNA" | cut -f1 )
start=$( grep ${gene} /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff | grep -E "mRNA" | cut -f4 )

echo ${start}


gene=Daphnia00789
stop=$( grep  ${gene} /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff | grep -E "mRNA" | cut -f5 | head -n1)

echo ${stop}

echo ${chr},${start},${stop} | tr ',' '\t' > ~/region.D84a.bed

cat ~/region.D84a.bed

~/liftOver \
~/region.D84a.bed \
liftover/reference_to_pa42/out.liftOver \
~/region.PA42.bed \
~/region.unmapped.bed

cat ~/region.PA42.bed: FLTH02000094.1  126368  149371
cat ~/region.unmapped.bed
