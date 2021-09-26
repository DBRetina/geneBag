REF=GENCODE_CHR38
REF_URL=http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
wget ${REF_URL} -O ${REF}.temp.fa.gz
gunzip ${REF}.temp.fa.gz

safe=1
cat ${REF}.temp.fa | while read line;do
if [[ "$line" == ">"* ]];then if [[ "$line" == *"_PAR_Y"* ]];then safe=0;else safe=1;fi;fi
if [ $safe -eq 1 ];then echo $line;fi
done > ${REF}.fa

