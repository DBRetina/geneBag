## Convert the refSeq genome annotation (from https://www.ncbi.nlm.nih.gov/genome/guide/human/) to reference transcriptome by gffread
## Annotate by the current NCBI Entrez gene IDs

conda install -c bioconda gffread
conda install -c bioconda samtools

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
gunzip GRCh38_latest_genomic.gff.gz
samtools faidx GRCh38_latest_genomic.fna
gffread -w GRCh38_latest_gffread.fa -g GRCh38_latest_genomic.fna GRCh38_latest_genomic.gff

echo "gff_id Entrez Symbol hgnc refseq gffAnn gbkey" | tr ' ' '\t' > GRCh38_latest_genomic.map
cat GRCh38_latest_genomic.gff | awk 'BEGIN{FS="\t";OFS=";"}!/#/{print $3,$9}' | awk 'BEGIN{FS=";";OFS="\t"}{ \
  sub("ID=","",$2); \
  delete vars; \
  for(i = 1; i <= NF; ++i){ \
     if($i ~ "^Dbxref="){sub("Dbxref=","",$i);split($i,A,","); \
                         for(a in A){n= index(A[a],":"); if(n){x= substr(A[a], n+1); vars[substr(A[a], 1, n-1)]=substr(A[a], n+1, length(x))}} \
                         id = vars["GeneID"]; acc = vars["Genbank"]; hgnc = vars["HGNC"];} \
     if($i ~ "^gbkey="){sub("gbkey=","",$i); gbkey = $i;} \
     if($i ~ "^gene="){sub("gene=","",$i); gene = $i;}} \
  if(vars["GeneID"])print $2,"Entrez:"id,(gene==""?"-":gene),(hgnc==""?"-":hgnc),(acc==""?"-":acc),$1,(gbkey==""?"-":gbkey)}' >> GRCh38_latest_genomic.map

awk 'BEGIN{S="|";}FNR==NR{a[">"$1]=">"$1 S $2 S $3 S $4 S $5 S $6 S $7 S;next;} \
     {if($0 ~ "^>"){if(a[$1]){print a[$1];safe=1}else safe=0;} \
     else if(safe==1) print $0;}' GRCh38_latest_genomic.map GRCh38_latest_gffread.fa > REFSEQ_GAnn.fa ## Sequences for 59677 Enterz gene IDs (4312 ID are still missing from the Homo_sapiens.gene_info file. These missing IDs do not have sequencing records e.g. disease associted loci)

cat REFSEQ_GAnn.fa | seqkit grep -r -p "\|MIR.*\|primary_transcript\|" -p "\|MIR.*\|miRNA\|" > REFSEQ_GAnn.miRNA.fa
cat REFSEQ_GAnn.fa | seqkit grep -r -p "\|[CDJV]_gene_segment\|" > REFSEQ_GAnn.immSeg.fa
cat REFSEQ_GAnn.fa | seqkit grep -v -r -p "\|MIR.*\|primary_transcript\|" -p "\|MIR.*\|miRNA\|" -p "\|[CDJV]_gene_segment\|" > REFSEQ_GAnn.simple.fa


