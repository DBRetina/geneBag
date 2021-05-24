# Best Mutual Gene Bag Match

KSIZE=25

VER38_GRCh38_URL=http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
VER28_GRCh37_URL=http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.transcripts.fa.gz

VER38_GRCh38=VER38_GRCh38.fa
VER28_GRCh37=VER28_GRCh37.fa

wget ${VER38_GRCh38_URL} -O ${VER38_GRCh38}.gz
wget ${VER28_GRCh37_URL} -O ${VER28_GRCh37}.gz

gunzip *gz

#1 Transform first annotation -------------------------------------------------------------
sed -i 's/^>/>VER38_GRCh38|/' ${VER38_GRCh38}

grep ">" ${VER38_GRCh38} | cut -c2- |  awk -F'|' '{print $0"\tver38_GRCh38."$3}' > ${VER38_GRCh38}.names

#2 Transform second annotation
sed -i 's/^>/>VER28_GRCh37|/' ${VER28_GRCh37}

grep ">" ${VER28_GRCh37} | cut -c2- |  awk -F'|' '{print $0"\tver28_GRCh37."$3}' > ${VER28_GRCh37}.names

#3 Indexing & Pairwise generation -------------------------------------------------------------

cat ${VER38_GRCh38} ${VER28_GRCh37} > merged.fa
cat ${VER38_GRCh38}.names ${VER28_GRCh37}.names > merged.fa.names

python index.py merged.fa ${KSIZE}

kSpider2 pairwise -i idx_merged

#4 Annotation & Filteration -------------------------------------------------------------

# Let us pool the output files into one final output 
# 1. merge the names map and nodes sizes in one output file
kPro_index="idx_merged"
paste <(tail -n+2 ${kPro_index}.namesMap |cut -d" " -f1)  <(tail -n+2 ${kPro_index}.namesMap |cut -d" " -f2-) > ${kPro_index}.namesMap.tmp
echo "node_id node_name size" | tr ' ' '\t' > ${kPro_index}_nodes_size.tsv
awk 'BEGIN{FS=OFS="\t";}FNR==NR{a[$2]=$3;next;}{if(a[$1]!="")print $0,a[$1]}' ${kPro_index}_kSpider_seqToKmersNo.tsv ${kPro_index}.namesMap.tmp >> ${kPro_index}_nodes_size.tsv
rm ${kPro_index}.namesMap.tmp*

# Let us have a look
head -n5 ${kPro_index}_nodes_size.tsv

# 2. Annotation of the numerically coded association file
# a) Add items names and no of associated items 
# b) calc jaccard distance and containment ratio for each pair
# c) Additionally, we can filter out those with minimal similarities 
md=20   ## minimum jaccard distance (as a percentage) to keep the record 
mc=80  ## minimum containment ratio (as a percentage) to keep the record 
echo ":START_ID|START_name|START_size|shared_count:int|jDist:float|smPerc:float|END_name|END_size|:END_ID" > ${kPro_index}_relations.csv
awk -v md=$md -v mc=$mc 'BEGIN{FS="\t";S="|";}FNR==NR{a[$1]=$3;b[$1]=$2S$3;next;}{
   g1=a[$2]; g2=a[$3]; min=g1;min=(min < g2 ? min : g2); 
   jDist=$4*100/(g1+g2-$4); smPerc=$4*100/min; 
   if(jDist>md || smPerc>mc)
     printf("%s%s%s%s%s%s%.1f%s%.1f%s%s%s%s\n", $2,S,b[$2],S,$4,S,jDist,S,smPerc,S,b[$3],S,$3)}' \
   ${kPro_index}_nodes_size.tsv <(tail -n+2 ${kPro_index}_kSpider_pairwise.tsv) >> ${kPro_index}_relations.csv

# Let us have a look
wc -l ${kPro_index}_relations.csv 
echo "----------------------"
head -n5 ${kPro_index}_relations.csv

#5 Automatic Detection of Best Mutual GeneBag Match -------------------------------------------------------------



#6  Best Match GeneBag query
