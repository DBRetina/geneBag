## Obtain the refSeq transcriptome (from https://www.ncbi.nlm.nih.gov/genome/guide/human/)
## Annotate by the current NCBI Entrez gene IDs
## Exclude non-current and non-chromosomal sequences 

REF=REFSEQ_CHRLatest
REF_URL=https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz ## matching: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_rna.fna.gz
wget ${REF_URL} -O ${REF}.temp.fa.gz
gunzip ${REF}.temp.fa.gz

wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz
gunzip GRCh38_latest_genomic.gff.gz
echo "refseq GeneID HGNC" | tr ' ' '\t' > refseq.map
cat GRCh38_latest_genomic.gff | awk -F"\t" '!/#/{if($9 ~ "^ID=rna-*")print $9}' | awk -F";" '{for(i = 1; i <= NF; ++i){if($i ~ "^Dbxref=")print $i}}' | grep "Genbank:" | sed 's/^Dbxref=//' | awk 'BEGIN{FS=",";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, ":"); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["GeneID"]; acc = vars["Genbank"]; hgnc = vars["HGNC"]; print acc,"Entrez:"id,hgnc; }' | sort | uniq >> refseq.map
## merge with Homo_sapiens.gene_info to select chr genes [1:22],MT,X,Y,X|Y
wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
gunzip Homo_sapiens.gene_info.gz
awk 'BEGIN{FS=OFS="\t";a["GeneID"]="Chromosome";}FNR==NR{a["Entrez:"$2]=$7;next;}{print $0,a[$2]}' Homo_sapiens.gene_info refseq.map > refseq_loc.map
awk 'BEGIN{FS=OFS="\t";for(i=1; i<=22; ++i)a[i]=1;a["X"]=a["Y"]=a["X|Y"]=a["MT"]=a["Chromosome"]=1;}{if(a[$4])print $0}' refseq_loc.map > refseq_chr.map
awk 'BEGIN{S="|";}FNR==NR{a[">"$1]=">"$1 S $2 S $3 S;next;} \
     {if($0 ~ "^>"){if(a[$1]){id=$1;sub($1 FS,"",$0);print a[id] $0 S;safe=1}else safe=0;} \
     else if(safe==1) print $0;}' refseq_chr.map ${REF}.temp.fa > ${REF}.fa  ## 103 (out of 164352) sequences were excluded. the remaining 164249 transcripts belong to 38469 GeneIDs from the 63981 IDs in Homo_sapiens.gene_info

