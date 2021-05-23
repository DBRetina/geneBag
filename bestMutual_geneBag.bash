# Best Mutual Gene Bag Match

KSIZE=25

VER38_GRCh38_URL=http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
VER28_GRCh37_URL=http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.transcripts.fa.gz

VER38_GRCh38=VER38_GRCh38.fa
VER28_GRCh37=VER28_GRCh37.fa

wget ${VER38_GRCh38_URL} -O ${VER38_GRCh38}.gz
wget ${VER28_GRCh37_URL} -O ${VER28_GRCh37}.gz

gunzip *gz

#1 Transform first annotation
sed -i 's/^>/>VER38_GRCh38|/' ${VER38_GRCh38}

grep ">" ${VER38_GRCh38} | cut -c2- |  awk -F'|' '{print $0"\tver38_GRCh38."$3}' > ${VER38_GRCh38}.names

#2 Transform second annotation
sed -i 's/^>/>VER28_GRCh37|/' ${VER28_GRCh37}

grep ">" ${VER28_GRCh37} | cut -c2- |  awk -F'|' '{print $0"\tver28_GRCh37."$3}' > ${VER28_GRCh37}.names

#3 Indexing & Pairwise generation

cat ${VER38_GRCh38} ${VER28_GRCh37} > merged.fa
cat ${VER38_GRCh38}.names ${VER28_GRCh37}.names > merged.fa.names

python index.py merged.fa ${KSIZE}


#4 Annotation & Filteration



#5 Automatic Detection of Best Mutual GeneBag Match



#6  Best Match GeneBag query
