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
cat GRCh38_latest_genomic.gff | awk 'BEGIN{FS="\t";OFS=";"}!/#/{print $3,$9}' | awk 'BEGIN{FS=";";OFS="\t"}{for(i = 1; i <= NF; ++i){if($i ~ "^Dbxref="){sub("ID=","",$2);sub("Dbxref=","",$i);split($i,A,",");delete vars; for(a in A){n= index(A[a],":"); if(n){x= substr(A[a], n+1); vars[substr(A[a], 1, n-1)] = substr(A[a], n+1, length(x))}} id = vars["GeneID"]; acc = vars["Genbank"]; hgnc = vars["HGNC"]; if(vars["GeneID"])print $1,$2,"Entrez:"id,hgnc,acc;}}}'  > GRCh38_latest_genomic.map
awk 'BEGIN{S="|";}FNR==NR{a[">"$2]=">"$2 S $3 S $1 S $4 S $5 S;next;} \
     {if($0 ~ "^>"){if(a[$1]){print a[$1];safe=1}else safe=0;} \
     else if(safe==1) print $0;}' GRCh38_latest_genomic.map GRCh38_latest_gffread.fa > REFSEQ_GAnn.fa

