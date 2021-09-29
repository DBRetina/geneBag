# Best Mutual Gene Bag Match
conda create --yes -n geneMatcher python=3.7
conda activate geneMatcher
pip install kSpider2 kProcessor


## Prepare the FASTA transcriptome files  where the header has the format ">transcriptID|geneID|optionalMetadata"
bash prep_refSeqTrans.sh ## output is REFSEQ_CHRLatest.fa
bash prep_refSeqGAnn.sh  ## output is REFSEQ_GAnn.fa
bash prep_gencode.sh     ## output is GENCODE_CHR38.fa

# Compare Refseq transcripts to GENCODE annotation
KSIZE=21
REF1=REFSEQ_CHRLatest
REF2=GENCODE_CHR38
bash geneMatcher.sh ${REF1} ${REF2} ${KSIZE}

# Compare Refseq genome annotation to GENCODE annotation
KSIZE=21
REF1=REFSEQ_GAnn
REF2=GENCODE_CHR38
bash geneMatcher.sh ${REF1} ${REF2} ${KSIZE}


KSIZE=8
REF1=REFSEQ_GAnn.miRNA
REF2=GENCODE_CHR38.miRNA
bash geneMatcher.sh ${REF1} ${REF2} ${KSIZE}

KSIZE=8
REF1=REFSEQ_GAnn.immSeg
REF2=GENCODE_CHR38.immSeg
bash geneMatcher.sh ${REF1} ${REF2} ${KSIZE}

KSIZE=21
REF1=REFSEQ_GAnn.simple
REF2=GENCODE_CHR38.simple
bash geneMatcher.sh ${REF1} ${REF2} ${KSIZE}



