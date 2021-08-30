#### 1. Ambiguity of human gene symbols within the same database:
#### Multiple genes in the same database might have the same gene symbol. 
#### The problem is much worse when we consider the alternative gene alias. 
#### Moreover, shuffling of symbols between genes add another hard-to-fix source of errors for meta-analysis studies. 
#### Some gene catalogs report the previous gene symbols but without clear versioning scheme, it is hard to resolve ambiguities. 
##########################

#### HGNC 
echo "Explore the gene symbols ambiguity in HGNC"
echo "=========================================="
## Download HGNC dataset (link in the "Statistics & download files" page)
if [ ! -f hgnc_complete_set.txt ];then 
  wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt;fi
if [ ! -f withdrawn.txt ];then
  wget http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/withdrawn.txt;fi

echo "Basic check of HGNC Ids:"
echo "-----------------------------"
app_IDs=$(tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS="\t";}{print $1}' | sort | uniq | wc -l)
app_tot=$(tail -n+2 hgnc_complete_set.txt | wc -l)
if (( app_tot != app_IDs));then echo "WARNING: The are $app_tot IDs in the approved HGNC but only $app_IDs are uniq. The approved HGNC has duplicate IDs!!";fi
nonApp=$(tail -n+2 hgnc_complete_set.txt | awk -F"\t" '{if($6!="Approved")print}' | wc -l)
if (( nonApp > 0));then echo "WARNING: There are $nonApp non-approved records in the approved dataset";fi

wd_IDs=$(tail -n+2 withdrawn.txt | awk 'BEGIN{FS="\t";}{print $1}' | sort | uniq | wc -l)
wd_tot=$(tail -n+2 withdrawn.txt | wc -l)
if (( wd_tot != wd_IDs));then echo "WARNING: The are $wd_tot IDs in the withdrawn HGNC dataset but only $wd_IDs are uniq. The withdrawn HGNC has duplicate IDs!!";fi

echo "Basic check is done"
echo "-------------------"
## Generate a map of gene IDs to symbols
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,"status",$2}' > hgnc.ID_to_Current
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,"Approved","<"$2">"}' >> hgnc.ID_to_Current 

## Generate a map of gene IDs to symbols and each one of the alias
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$9}' > hgnc.ID_to_EachAlias
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($9!="")print $1,$2,$9}' | sed 's/"//g' | awk 'BEGIN{FS="\t";OFS="\n";}{split($3,a,"|");for(i in a)print $1"\t"$2"\t<"a[i]">";}' >> hgnc.ID_to_EachAlias  

## Generate a map of gene IDs to symbols and each one of the previous symbols
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$11}' > hgnc.ID_to_EachPrev
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($11!="")print $1,$2,$11}' | sed 's/"//g' | awk 'BEGIN{FS="\t";OFS="\n";}{split($3,a,"|");for(i in a)print $1"\t"$2"\t<"a[i]">";}' >> hgnc.ID_to_EachPrev  

## Generate a map for genes IDs withdrawn without approved replacement
## Known limitation: We are considering genes replaced by withdrawn genes to be discontinued. This is a possiblity that this new withdrawn gene is also replaced by new approved gene 
##                   but I do not see any example in the current DB (if a gene is replaced, the new genes (column 4) are either "Approved" or "Entry Withdrawn" but not "Merged/Split") 
head -n1 withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$3}' > hgnc.ID_to_discontinued 
cat withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{if($2=="Entry Withdrawn")print $1,$2,"<"$3">"}' >> hgnc.ID_to_discontinued 
cat withdrawn.txt | grep -v "Approved" | awk 'BEGIN{FS=OFS="\t";}{if($2=="Merged/Split")print $1,$2,"<"$3">"}'  >> hgnc.ID_to_discontinued  

## Generate a map for genes IDs withdrawn but replaced by new approved IDs
head -n1 withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$3}' > hgnc.ID_to_replaced 
cat withdrawn.txt | grep "Approved" | awk 'BEGIN{FS=OFS="\t";}{if($2=="Merged/Split")print $1,$2,"<"$3">"}' >> hgnc.ID_to_replaced 


## Identify genes with ambiguous alias or previous symbol (the alias or previous symbol is ambiguous if it matches another alias, previous or current gene symbol)
# create list of all gene symbols
tail -n+2 hgnc_complete_set.txt | awk -F"\t" '{print "<"$2">"}' | sort > hgnc.Symbols  
# create list of all alias symbols
tail -n+2 hgnc.ID_to_EachAlias | awk -F "\t" '{print $3}' | sort > hgnc.Alias 
# create list of all previous symbols
tail -n+2 hgnc.ID_to_EachPrev | awk -F "\t" '{print $3}' | sort > hgnc.Prev
# create list of withdrawn symbols without approved replacement 
tail -n+2 hgnc.ID_to_discontinued  | awk -F "\t" '{print $3}' | sort > hgnc.discontinued 
# create list of withdrawn symbols with approved replacement 
tail -n+2 hgnc.ID_to_replaced | awk -F "\t" '{print $3}' | sort > hgnc.replaced

## stats
echo "HGNC approved symbols     = " $(cat hgnc.Symbols | wc -l)        ## 42698
echo "HGNC alias symbols        = " $(cat hgnc.Alias | wc -l)          ## 42347
echo "HGNC previous symbols     = " $(cat hgnc.Prev | wc -l)           ## 15180
echo "HGNC discontinued symbols = " $(cat hgnc.discontinued | wc -l)   ## 1826
echo "HGNC replaced symbols     = " $(cat hgnc.replaced | wc -l)       ## 3314

## Gene ambiguity venn diagram
cat hgnc.Symbols | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.01.Current_symbols_matching_other_current_symbols
cat hgnc.Alias | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.02.Alias_symbols_matching_other_alias_symbols
comm -12 <(cat hgnc.Alias | uniq) <(cat hgnc.Symbols | uniq) > HGNC.03.Alias_symbols_matching_current_symbols
cat hgnc.Prev | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.04.Previous_symbols_matching_other_previous_symbols
comm -12 <(cat hgnc.Prev | uniq) <(cat hgnc.Symbols | uniq) > HGNC.05.Previous_symbols_matching_current_symbols
comm -12 <(cat hgnc.Prev | uniq) <(cat hgnc.Alias | uniq)  > HGNC.06.Previous_symbols_matching_alias_symbols
cat hgnc.discontinued | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.07.Discontinued_symbols_matching_other_discontinued_symbols
comm -12 <(cat hgnc.discontinued | uniq) <(cat hgnc.Symbols | uniq) > HGNC.08.Discontinued_symbols_matching_current_symbols
comm -12 <(cat hgnc.discontinued | uniq) <(cat hgnc.Alias | uniq) > HGNC.09.Discontinued_symbols_matching_alias_symbols
comm -12 <(cat hgnc.discontinued | uniq) <(cat hgnc.Prev | uniq) > HGNC.10.Discontinued_symbols_matching_previous_symbols
cat hgnc.replaced | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNC.11.Replaced_symbols_matching_other_replaced_symbols
comm -12 <(cat hgnc.replaced | uniq) <(cat hgnc.Symbols | uniq) > HGNC.12.Replaced_symbols_matching_current_symbols
comm -12 <(cat hgnc.replaced | uniq) <(cat hgnc.Alias | uniq) > HGNC.13.Replaced_symbols_matching_alias_symbols
comm -12 <(cat hgnc.replaced | uniq) <(cat hgnc.Prev | uniq) > HGNC.14.Replaced_symbols_matching_previous_symbols
comm -12 <(cat hgnc.replaced | uniq) <(cat hgnc.discontinued | uniq) > HGNC.15.Replaced_symbols_matching_discontinued_symbols

wc -l HGNC.*_matching_*_symbols 

cat hgnc.{Symbols,Alias,Prev,discontinued,replaced} | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -nr  > HGNC.ambiguous_freq
cat HGNC.ambiguous_freq | awk '{print $2}' | grep -Fwf - <(cat hgnc.ID_to_{Current,EachAlias,EachPrev,discontinued,replaced}) | sort -t$'\t' -k3,3 >  hgnc.ambiguous.temp 
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$6,$9,$11}' > hgnc.complete_withdrawn.temp
tail -n+2 withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$3,$2}' >> hgnc.complete_withdrawn.temp
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print "<Ambiguous_Symbol>",$1,$2,$6,$9,$11}' > HGNC.ambiguous
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{print $3,a[$1]}' hgnc.complete_withdrawn.temp hgnc.ambiguous.temp >> HGNC.ambiguous

echo "HGNC has "$(cat HGNC.ambiguous_freq | wc -l)" ambigious symbols causing "$(tail -n+2 HGNC.ambiguous | wc -l)" ambigious records."
echo "Here are the most 10 ambiguous symbols and how many time do they show up among all gene symbols:"
head HGNC.ambiguous_freq
echo "-------------------------"

echo "Track replaced HGNC ID"
echo "----------------------"
cat withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{if($2=="Merged/Split"){split($4,a,", ");for(i in a){split(a[i],b,"|");if(b[3]=="Approved")print b[1],b[2],"<"$3">";}}}' | sort | uniq > wdHGNC.ID_to_EachPrev ## the terminal sort|uniq is to overcome a bug in the current HGNC report where the HGNC:35188 shows up twice in the same record of the withdrawn HGNC:21128 ID
echo "No. of withdrawn symbols replaced by new approved symbols = "$(cat wdHGNC.ID_to_EachPrev | wc -l)

comm -12 <(sort hgnc.ID_to_EachAlias) <(sort wdHGNC.ID_to_EachPrev) > wdHGNC.ID_to_EachPrev_AsAlias
echo "No. of withdrawn symbols used as alias symbols = "$(cat wdHGNC.ID_to_EachPrev_AsAlias | wc -l)

comm -13 <(sort hgnc.ID_to_EachAlias) <(sort wdHGNC.ID_to_EachPrev) > wdHGNC.ID_to_EachPrev2
comm -12 <(sort hgnc.ID_to_EachPrev) <(sort wdHGNC.ID_to_EachPrev2) > wdHGNC.ID_to_EachPrev_AsPrev
echo "No. of withdrawn symbols used as previous symbols = "$(cat wdHGNC.ID_to_EachPrev_AsPrev | wc -l)

comm -13 <(sort hgnc.ID_to_EachPrev) <(sort wdHGNC.ID_to_EachPrev2) > wdHGNC.ID_to_EachPrev.missing
echo "No. of missing withdrawn symbols = "$(cat wdHGNC.ID_to_EachPrev.missing | wc -l)
echo "----------------------"
##################################################################################################################
#### NCBI genes
echo "Explore the gene symbols ambiguity in Entrez gene DB"
echo "===================================================="
## Download
if [ ! -f Homo_sapiens.gene_info ];then
  wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
  gunzip Homo_sapiens.gene_info.gz
fi
if [ ! -f human_gene_history_track ];then
  wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz
  gunzip gene_history.gz
  head -n1 gene_history > human_gene_history
  grep -w ^9606 gene_history >> human_gene_history
  head -n1 human_gene_history | awk 'BEGIN{FS=OFS="\t";}{print $0,"status"}' > human_gene_history_track
  awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$2;next;}!/^#/{if($2=="-"){print $0;}
                                               else{new_id=$2;while(1){ \
                                                 if(a[new_id]=="-"){print $0,"<discontinued>";break;} \
                                                 else if(!a[new_id]){print $0,"<replaced>";break;} \
                                                 else new_id=a[new_id];}}}' human_gene_history human_gene_history >> human_gene_history_track
fi

echo "Basic check of Entrez Ids:"
echo "--------------------------"
ncbi_IDs=$(tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{print $2}' | sort | uniq | wc -l) ## 61622
ncbi_tot=$(tail -n+2 Homo_sapiens.gene_info | wc -l)
if (( ncbi_tot != ncbi_IDs));then echo "WARNING: The are $ncbi_tot IDs for NCBI genes but only $ncbi_IDs are uniq. There are NCBI genes with duplicate IDs!!";fi

ncbiHis_IDs=$(tail -n+2 human_gene_history | awk 'BEGIN{FS="\t";}{print $3}' | sort | uniq | wc -l)
ncbiHis_tot=$(tail -n+2 human_gene_history | wc -l)
if (( ncbiHis_tot != ncbiHis_IDs));then echo "WARNING: The are $ncbiHis_tot IDs for discontinued NCBI genes but only $ncbiHis_IDs are uniq. There are discontinued NCBI genes with duplicate IDs!!";fi

echo "Basic check is done"
echo "-------------------"

## Generate a map of gene IDs to symbols
head -n1 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,"status",$3}' > entrez.ID_to_Current
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,"Official","<"$3">"}' >> entrez.ID_to_Current 

## Generate a map of gene IDs to symbols and each one of the alias
head -n1  Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,$5}' > entrez.ID_to_EachAlias
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{if($5!="-")print $2,$3,$5}' | awk 'BEGIN{FS="\t";OFS="\n";}{split($3,a,"|");for(i in a)print $1"\t"$2"\t<"a[i]">";}' >> entrez.ID_to_EachAlias 

## Generate a map for genes IDs withdrawn without replacement 
head -n1 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{print $3,$6,$4}' > entrez.ID_to_discontinued 
cat human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($2=="-" || $6=="<discontinued>")print $3,"discontinued","<"$4">"}' >> entrez.ID_to_discontinued 

## Generate a map for genes IDs withdrawn but replaced by new IDs
head -n1 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{print $3,$6,$4}' > entrez.ID_to_replaced
cat human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($6=="<replaced>")print $3,"replaced","<"$4">"}' >> entrez.ID_to_replaced


## Identify genes with ambiguous alias (the alias is ambiguous if it matches another alias, previous or current gene symbol)
# create list of all gene symbols
tail -n+2 Homo_sapiens.gene_info | awk -F"\t" '{print "<"$3">"}' | sort > entrez.Symbols 
# create list of all alias symbols
tail -n+2 entrez.ID_to_EachAlias | awk -F "\t" '{print $3}' | sort > entrez.Alias 
# create list of withdrawn symbols without new replacement 
tail -n+2 entrez.ID_to_discontinued  | awk -F "\t" '{print $3}' | sort > entrez.discontinued 
# create list of withdrawn symbols with new replacement 
tail -n+2 entrez.ID_to_replaced | awk -F "\t" '{print $3}' | sort > entrez.replaced

## stats
echo "Entrez official symbols     = " $(cat entrez.Symbols | wc -l)        ## 63881
echo "Entrez alias symbols        = " $(cat entrez.Alias | wc -l)          ## 71213
echo "Entrez discontinued symbols = " $(cat entrez.discontinued | wc -l)   ## 142208
echo "Entrez replaced symbols     = " $(cat entrez.replaced | wc -l)       ## 21743

## Gene ambiguity venn diagram
cat entrez.Symbols | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.01.Current_symbols_matching_other_current_symbols
cat entrez.Alias | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.02.Alias_symbols_matching_other_alias_symbols
comm -12 <(cat entrez.Alias | uniq) <(cat entrez.Symbols | uniq) > Entrez.03.Alias_symbols_matching_current_symbols
cat entrez.discontinued | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.07.Discontinued_symbols_matching_other_discontinued_symbols
comm -12 <(cat entrez.discontinued | uniq) <(cat entrez.Symbols | uniq) > Entrez.08.Discontinued_symbols_matching_current_symbols
comm -12 <(cat entrez.discontinued | uniq) <(cat entrez.Alias | uniq) > Entrez.09.Discontinued_symbols_matching_alias_symbols
cat entrez.replaced | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.11.Replaced_symbols_matching_other_replaced_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.Symbols | uniq) > Entrez.12.Replaced_symbols_matching_current_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.Alias | uniq) > Entrez.13.Replaced_symbols_matching_alias_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.discontinued | uniq) > Entrez.15.Replaced_symbols_matching_discontinued_symbols

wc -l Entrez.*_matching_*_symbols

cat entrez.{Symbols,Alias,discontinued,replaced} | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -nr  > Entrez.ambiguous_freq
cat Entrez.ambiguous_freq | awk '{print $2}' | grep -Fwf - <(cat entrez.ID_to_{Current,EachAlias,discontinued,replaced}) | sort -t$'\t' -k3,3 >  entrez.ambiguous.temp 
head -n1 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,"status",$5}' > entrez.complete_withdrawn.temp
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,"Official",$5}' >> entrez.complete_withdrawn.temp
tail -n+2 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($2=="-")print $3,$4,"discontinued";else if($6=="<discontinued>")print $3,$4,"replaced then discontinued";else print $3,$4,"replaced";}' >> entrez.complete_withdrawn.temp
head -n1 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print "<Ambiguous_Symbol>",$2,$3,"status",$5}' > Entrez.ambiguous
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{print $3,a[$1]}' entrez.complete_withdrawn.temp entrez.ambiguous.temp >> Entrez.ambiguous

echo "Entrez Gene DB has "$(cat Entrez.ambiguous_freq | wc -l)" ambigious symbols causing "$(tail -n+2 Entrez.ambiguous | wc -l)" ambigious records."
echo "Here are the most 10 ambiguous symbols and how many time do they show up among all gene symbols:"
head Entrez.ambiguous_freq
echo "-------------------------"


echo "Track replaced Entrez ID"
echo "-------------------------"
tail -n+2 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($2!="-" && $6!="<discontinued>")print $2,"<"$4">";}' > wdEntrez.ID_to_EachPrev
echo "No. of withdrawn symbols replaced by new approved symbols = "$(cat wdEntrez.ID_to_EachPrev | wc -l)

comm -12 <(cat entrez.ID_to_EachAlias | awk 'BEGIN{FS=OFS="\t";}{print $1,$3}' | sort) <(sort wdEntrez.ID_to_EachPrev) > wdEntrez.ID_to_EachPrev_AsAlias
echo "No. of withdrawn symbols used as alias symbols = "$(cat wdEntrez.ID_to_EachPrev_AsAlias | wc -l)

comm -13 <(cat entrez.ID_to_EachAlias | awk 'BEGIN{FS=OFS="\t";}{print $1,$3}' | sort) <(sort wdEntrez.ID_to_EachPrev) > wdEntrez.ID_to_EachPrev.missing
echo "No. of missing withdrawn symbols = "$(cat wdEntrez.ID_to_EachPrev.missing | wc -l)

grep -v "<LOC.*>" wdEntrez.ID_to_EachPrev.missing > wdEntrez.ID_to_EachPrev.missing.knownIDs
echo "No. of missing withdrawn symbols after exclusion of LOC* IDs = "$(cat wdEntrez.ID_to_EachPrev.missing.knownIDs | wc -l)
echo "-------------------------"

##################################################################################################################
#### Gencode_human
vcur=38
vold=7; dateOld="Dec, 2010"; 
echo "Explore the gene symbols ambiguity in Gencode $vcur"
echo "=========================================="
## Download
if [ ! -f gencode.v${vcur}.annotation.gtf ];then
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${vcur}/gencode.v${vcur}.annotation.gtf.gz
  gunzip gencode.v${vcur}.annotation.gtf.gz
fi

## Generate a map of gene IDs to symbols and list of cross references to other DB
echo "gene_id gene_type gene_name hgnc_id havana_gene" | tr ' ' '\t' > gencode.ann
cat gencode.v${vcur}.annotation.gtf | awk -F"\t" '!/#/{if($3=="gene")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); key = substr($i, 1, n - 1); val = substr($i, n + 1, length(x));if(vars[key]=="")vars[key] = val;else vars[key] = vars[key]","val;} } id = vars["gene_id"]; ann = vars["gene_type"]; name = vars["gene_name"]; hgnc = vars["hgnc_id"]; hav = vars["havana_gene"]; print id,ann,name,hgnc,hav; }' >> gencode.ann

gencode_tot=$(tail -n+2 gencode.ann | wc -l)
gencode_IDs=$(tail -n+2 gencode.ann | awk 'BEGIN{FS="\t";}{print $1}' | sort | uniq | wc -l)
gencode_sym=$(tail -n+2 gencode.ann | awk 'BEGIN{FS="\t";}{print $3}' | sort | uniq | wc -l)
if (( gencode_tot != gencode_IDs));then echo "WARNING: The are $gencode_tot IDs in gencode annotation but only $gencode_IDs are uniq. There are gencode genes with duplicate IDs";fi
PAR_Y=$(grep "_PAR_Y" gencode.ann | wc -l)
echo "If we considered the unversioned IDs, rhere are $PAR_Y duplicate IDs with the PAR_Y suffix for Y chromosome version of genes"


gencode_status=""
> gencode.ambiguous_report
if (( gencode_sym != gencode_IDs));then gencode_status="WARNING: There are duplicate gene symbols in the current gencode annotation.";
  tail -n+2 gencode.ann | awk -F"\t" '{print $3}' | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -nr > gencode.ambiguous_freq
  gencode_dedup_ids=$(cat gencode.ambiguous_freq | awk '{a+=$1}END{print a}')
  gencode_dedup_sym=$(wc -l gencode.ambiguous_freq)
  echo "There are $gencode_dedup_sym symbols assigned to $gencode_dedup_ids genes" > gencode.ambiguous_report
  echo "Here are the most ambiguous symbols in the current gencode annotation:" >> gencode.ambiguous_report
  head gencode.ambiguous_freq >> gencode.ambiguous_report
else gencode_status="There are no duplicate gene symbols in the current gencode annotation.";fi
echo "$gencode_status The current gencode annotation have $gencode_IDs IDs and corresponding $gencode_sym symbols."
cat gencode.ambiguous_report
echo "-------------------------"
echo Note: Gencode annotation files do not have gene information about gene alias or previous gene symbols.
echo "-------------------------"

echo "Compare Gencode v${vcur} versus v${vold}"
echo "=========================================="
if [ ! -f gencode.v${vold}.annotation.gtf ];then
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${vold}/gencode.v${vold}.annotation.gtf.gz
  gunzip gencode.v${vold}.annotation.gtf.gz
fi

echo "gene_id gene_type gene_name hgnc_id havana_gene" | tr ' ' '\t' > gencode.oldAnn
cat gencode.v${vold}.annotation.gtf | awk -F"\t" '!/#/{if($3=="gene")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); key = substr($i, 1, n - 1); val = substr($i, n + 1, length(x));if(vars[key]=="")vars[key] = val;else vars[key] = vars[key]","val;} } id = vars["gene_id"]; ann = vars["gene_type"]; name = vars["gene_name"]; hgnc = vars["hgnc_id"]; hav = vars["havana_gene"]; print id,ann,name,hgnc,hav; }' >> gencode.oldAnn

genOld_tot=$(tail -n+2 gencode.oldAnn | wc -l)
echo "While Gencode ${vcur} has $gencode_tot genes, Gencode ${vold} had $genOld_tot genes"

head -n1 gencode.ann | awk 'BEGIN{FS=OFS="\t"}{$1=$1 FS "version";print}' > gencode.ann2
tail -n+2 gencode.ann | sed 's/\./|/' | tr "|" "\t" >> gencode.ann2
tail -n+2 gencode.ann2 | awk -F"\t" '{print $1}' | sort | uniq > gen_${vcur}.IDs

head -n1 gencode.oldAnn | awk 'BEGIN{FS=OFS="\t"}{$1=$1 FS "version";print}' > gencode.oldAnn2
tail -n+2 gencode.oldAnn | sed 's/\./|/' | tr "|" "\t" >> gencode.oldAnn2
tail -n+2 gencode.oldAnn2 | awk -F"\t" '{print $1}' | sort | uniq > gen_${vold}.IDs

echo "After exclusion of duplicate Y chromosome IDs, we have 3 categories of IDs"
echo "shared IDs in both :" $(comm -12 gen_${vcur}.IDs gen_${vold}.IDs | wc -l) ## 43531 shared IDs
echo "novel IDs in ${vcur}:" $(comm -23 gen_${vcur}.IDs gen_${vold}.IDs | wc -l) ## 17074 novel IDs
echo "discontinued IDs in ${vold}:" $(comm -13 gen_${vcur}.IDs gen_${vold}.IDs | wc -l) ## 7551 discontinued IDs


echo "gene_id version.${vold} symbol.${vold} version.${vcur} symbol.${vcur}" | tr " " "\t" > gencode.v${vold}.vs.v${vcur}
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2 FS $4;next;}{if(a[$1]!="")print $1,$2,$4,a[$1]}' <(tail -n+2 gencode.ann2) <(tail -n+2 gencode.oldAnn2) >> gencode.v${vold}.vs.v${vcur}
updated_sym=$(tail -n+2 gencode.${vold}.vs.${vcur} | awk -F"\t" '{if($3!=$5)a+=1}END{print a}') ## 22402
echo "Among the gene with shared IDs, $updated_sym has updated gene symbols"

##################################################################################################################
##################################################################################################################
##################################################################################################################
#### 2. Different catalogs are built on different gene annotations: 
#### There is a significant non-overlapping between different gene catalogs 
#### Evenmore there is inconsistent cross-referencing between different catalogs

#### HGNC
## Explore
# Explore types of genes
echo "Approved HGNC genes can be calssified as:"
tail -n+2 hgnc_complete_set.txt | awk -F"\t" '{if($6=="Approved")print $4 " - " $5}' | sort | uniq -c 

## Generate a map of approved HGNC gene IDs to dbXrefs
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$19,$20,$9,$11}' > hgnc_approved.map ## hgnc_id symbol  entrez_id       ensembl_gene_id alias_symbol    prev_symbol
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($6=="Approved")print $1,$2,$19,$20,$9,$11}' >> hgnc_approved.map
 
## Stats
# tail -n+2 hgnc_approved.map | awk -F"\t" '{a+=1;if($3!="")b+=1;if($4!="")c+=1;if($5!="")d+=1;if($6!="")e+=1;}END{print "hgnc_id=",a,"entrez_id=",b,"ensembl_id=",c,"alias_symbol=",d,"prev_symbol=",e}'  ## hgnc_id= 42181 entrez_id= 42110 ensembl_id= 39199 alias_symbol= 21740 prev_symbol= 12049
tail -n+2 hgnc_approved.map | awk 'BEGIN{FS="\t";OFS="\n";}{a+=1;if($3!="")b+=1;if($4!="")c+=1;if($5!="")d+=1;if($6!="")e+=1;}END{print "There are "a" approved HGNC ids with "d" aliases and "e" previous symbols", "They are cross referenced with:", "Entrez IDs = "b, "Ensembl IDs = "c}' 
tail -n+2 hgnc_approved.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "HGNC genes with entrez and ensembl ids =",a}'  ## hgnc_id with entrez and ensembl ids = 39578
#tail -n+2 hgnc_approved.map | awk -F"\t" '{a[$2]=1;if($3!="")b[$3]=1;if($4!="")c[$4]=1;}END{print "uniq_symb=",length(a),"uniq_entrez=",length(b),"uniq_ens=",length(c)}' ## uniq_symb= 42181 uniq_entrez= 42110 uniq_ens= 39196



#### NCBI (Entrez gene ids from Homo_sapiens.gene_info which does not have the discontinued genes)
## Explore
## col 10 = type_of_gene, 11 = Symbol_from_nomenclature_authority, 13 = Nomenclature_status
# Explore types of genes:
echo "Current NCBI genes can be calssified as:"
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{print $10}' | sort | uniq -c
echo "NOTE: NCBI IDs labeled as [biological-region] are actually no genes but usually refer to genomic regions with know biological function e.g. trascription regulatory regions"
echo "NOTE: NCBI IDs labeled as [unknown] are usually refering to loci but the exact gene coordinates are not known yet"

# genes with non authorized names
#tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($3!=$11)print;}' > non_authority_naming #19553
# Classification of genes based on naming source (either "O" or "-". The "O" means nomenclature_authority i.e. HGN)
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{print $13}' | sort | uniq -c ## 19516 - && 42106 O   
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($13=="O")print;}' > Nomenclature_status_O #42106
grep "HGNC:" Nomenclature_status_O | wc -l     ## 42106
grep "Ensembl:" Nomenclature_status_O | wc -l  ## 33165
#tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($13=="O")print $10}' | sort | uniq -c
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($13=="O" && $3!=$11)print}' | wc -l # 37
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($13!="O")print;}' > Nomenclature_status_NA #42106
grep "HGNC:" Nomenclature_status_NA | wc -l     ## 0
grep "Ensembl:" Nomenclature_status_NA | wc -l  ## 1898

## Generate a map of gene IDs to dbXrefs
head -n1 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,"HGNC","Ensembl","MIM",$5}' > Homo_sapiens.gene_info.map ## GeneID  Symbol  HGNC    Ensembl MIM     Synonyms
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{if($6!="-")print $2,$3,$6,$5}' | awk 'BEGIN{FS=OFS="\t"}{ delete vars; split($3,a,"|");for(i in a) { n = index(a[i], ":"); if(n) { x = substr(a[i], n + 1); key = substr(a[i], 1, n - 1); val = substr(a[i], n + 1, length(x)); if(vars[key]=="")vars[key] = val;else vars[key] = vars[key]","val; } } MIM = vars["MIM"]; HGNC = vars["HGNC"]; Ensembl = vars["Ensembl"]; print $1,$2,HGNC,Ensembl,MIM,$4; }' >> Homo_sapiens.gene_info.map
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{if($6=="-")print $2,$3,"","","",$5}' >> Homo_sapiens.gene_info.map ## $6 = dbXrefs

## Stats
#tail -n+2 Homo_sapiens.gene_info.map | awk -F"\t" '{a+=1;if($3!="")b+=1;if($4!="")c+=1;if($5!="")d+=1;if($6!="-")e+=1;}END{print "entrez_id=",a,"hgnc_id=",b,"ensembl_id=",c,"mim_id=",d,"alias_symbol=",e}'  ## entrez_id= 61622 hgnc_id= 42106 ensembl_id= 35063 mim_id= 17648 alias_symbol= 26915
tail -n+2 Homo_sapiens.gene_info.map | awk 'BEGIN{FS="\t";OFS="\n";}{a+=1;if($3!="")b+=1;if($4!="")c+=1;if($5!="")d+=1;if($6!="-")e+=1;}END{print "There are "a" current Entrez ids with "e" aliases", "They are cross referenced with:", "HGNC IDs = "b, "Ensembl IDs = "c, "MIM IDs = "d}'
tail -n+2 Homo_sapiens.gene_info.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "Entrez genes with hgnc and ensembl IDs =",a}'  ## entrez_id with hgnc and ensembl ids= 33165
#tail -n+2 Homo_sapiens.gene_info.map | awk -F"\t" '{a[$2]=1;if($3!="")b[$3]=1;if($4!="")c[$4]=1;}END{print "uniq_symb=",length(a),"uniq_hgnc=",length(b),"uniq_ens=",length(c)}' ## uniq_symb= 61563 uniq_hgnc= 42106 uniq_ens= 34966



#### Gencode_human
if [ ! -f gencode.v${vcur}.metadata.EntrezGene ];then
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${vcur}/gencode.v${vcur}.metadata.EntrezGene.gz
  gunzip gencode.v${vcur}.metadata.EntrezGene.gz
fi

## Explore
# Explore types of genes:
echo "Current Gencode genes can be calssified as:"
tail -n+2 gencode.ann | awk 'BEGIN{FS="\t";}{print $2}' | sort | uniq -c

## Generate a map of approved Gencide gene IDs to dbXrefs
echo "gene_id transcript_id gene_name hgnc_id havana_gene" | tr ' ' '\t' > gencode.trans_ann
cat gencode.v${vcur}.annotation.gtf | awk -F"\t" '!/#/{if($3=="transcript")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; trans = vars["transcript_id"]; name = vars["gene_name"]; hgnc = vars["hgnc_id"]; hav = vars["havana_gene"]; print id,trans,name,hgnc,hav; }' >> gencode.trans_ann

head -n1 gencode.trans_ann | awk 'BEGIN{FS=OFS="\t"}{print $0,"EntrezGene"}' > gencode.trans_ann.map
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{ print $0, a[$2]}' gencode.v${vcur}.metadata.EntrezGene <(tail -n+2 gencode.trans_ann) >> gencode.trans_ann.map

cat gencode.trans_ann.map | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$4,$6}' | uniq > gencode.gene_ann.map

## Stats
#tail -n+2 gencode.gene_ann.map | awk -F"\t" '{a+=1;if($3!="")b+=1;if($4!="")c+=1;}END{print "ensembl_id=",a,"hgnc_id=",b,"entrez_id=",c}'  ## ensembl_id= 60656 hgnc_id= 38596 entrez_id= 25600
tail -n+2 gencode.gene_ann.map | awk 'BEGIN{FS="\t";OFS="\n";}{a+=1;if($3!="")b+=1;if($4!="")c+=1;}END{print "There are "a" current Ensembl ids", "They are cross referenced with:", "HGNC IDs = "b, "Entrez IDs = "c}'  
tail -n+2 gencode.gene_ann.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "Gencode gene with hgnc and entrez  ids =",a}'  ## gencode_id with hgnc and entrez  ids= 24528
#tail -n+2 gencode.v35.map | awk -F"\t" '{a[$2]=1;if($3!="")b[$3]=1;if($4!="")c[$4]=1;}END{print "uniq_symb=",length(a),"uniq_hgnc=",length(b),"uniq_entrez=",length(c)}' ## uniq_symb= 59609 uniq_hgnc= 38543 uniq_entrez= 25532

##################################################################################################################
##################################################################################################################
##################################################################################################################
#### 3. discrepancies between databases

## NCBI and HGNC
## merge HGNC ids and symbols into the NCBI DB to compare
head -n1 Homo_sapiens.gene_info.map | awk 'BEGIN{FS=OFS="\t"}{print $0,"hgnc-hgnc_id","hgnc-symbol"}' > NCBI.map.hgncExt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$1 FS $2;next}{print $0,a[$1]}' <(tail -n+2 hgnc_approved.map) <(tail -n+2 Homo_sapiens.gene_info.map) >> NCBI.map.hgncExt

## Table 1
echo "- NCBI (Entrez) genes which have HGNC ids and symbols in NCBI different from those in HGNC database:"
echo "Entrez_id Entrez_Symbol ncbi-hgnc_id hgnc-hgnc_id hgnc-Symbol" | tr ' ' '\t'
tail -n+2 NCBI.map.hgncExt | awk 'BEGIN{FS=OFS="\t"}{if($3!=$7)print $1,$2,($3==""?"-":$3),($7==""?"-":$7),($8==""?"-":$8)}' ## (7 records)
echo "-------------------------"

## Table 2
echo "- NCBI (Entrez) genes with the same HGNC ids in NCBI and HGNC database but with different symbols: "
echo "Entrez_id Entrez_Symbol ncbi-hgnc_id hgnc-hgnc_id hgnc-Symbol" | tr ' ' '\t'
tail -n+2 NCBI.map.hgncExt | awk 'BEGIN{FS=OFS="\t"}{if($3!="" && $3==$7 && $2!=$8)print $1,$2,$3,$7,($8==""?"-":$8)}' | head ## (37 mitochondrial genes)
echo "-------------------------"

## Table 3
echo "- Current approved HGNC genes with discontinued Enterz ID"
head -n1 hgnc_approved.map | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,"discontinued_entrez_id","discontinued_entrez_symbol"}'
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$4;next}{if(a[$3])print $1,$2,$3,a[$3]}' <(tail -n+2 human_gene_history) <(tail -n+2 hgnc_approved.map) 
echo "-------------------------"


## Gencode and HGNC
## Table 1
echo "- In HGNC database, multiple HGNC IDs map to the same Ensembl ID"
head -n1 hgnc_approved.map
tail -n+2 hgnc_approved.map | awk -F"\t" '{if($4!="")print $4}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' | grep -Fwf - hgnc_approved.map ## 
echo "-------------------------"

## Table 2
echo "- In current Gencode annotation, multiple Ensembl IDs map to the same HGNC ID "
head -n1 gencode.gene_ann.map
tail -n+2 gencode.gene_ann.map | grep -v "_PAR_Y" | awk -F"\t" '{if($3!="")print $3}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' | grep -Fwf - gencode.gene_ann.map | sort -k3,3 ## there are 17 pairs of gencode IDs where each pair maps to one  HGNC IDs. Each pair has the gene symbol as well. Similarly, each pair has the same Entrez ID but - unexpectedly - one member in 12 (out of the 17) pairs is missing the Entrez ID.
echo "-------------------------"

## merge HGNC ids and symbols into Gencode to compare
cat gencode.gene_ann.map | awk 'BEGIN{FS=OFS="\t"}{split($1,a,".");print a[1],$2,$3,$4}' > gencode.noGeneVer.map
head -n1 gencode.noGeneVer.map | awk 'BEGIN{FS=OFS="\t"}{print $0,"hgnc-hgnc_id","hgnc-symbol","hgnc-hgnc_id2","hgnc-symbol2"}' > gencode.map.hgncExt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if(a[$4]=="")a[$4]=$1 FS $2;else a[$4]=a[$4] FS $1 FS $2;next}{print $0,a[$1]}' <(tail -n+2 hgnc_approved.map) <(tail -n+2 gencode.noGeneVer.map) >> gencode.map.hgncExt

## Table 3
echo "- Gencode genes which have HGNC ids and symbols in Gencode annotation different from those in HGNC database:"
echo "Gencode_id Gencode_Symbol Gencode-hgnc_id Gencode-Entrez_id hgnc-hgnc_id hgnc-Symbol" | tr ' ' '\t'
tail -n+2 gencode.map.hgncExt | awk -F"\t" '{if(3!="" && $5!="" && $3!=$5 && $3!=$7)print}' 
echo "-------------------------"

## Table 4
echo "- Gencode genes with the same HGNC ids in Gencode and HGNC database but with different symbols:"
echo "Gencode_id Gencode_Symbol Gencode-hgnc_id Gencode-Entrez_id hgnc-hgnc_id hgnc-Symbol" | tr ' ' '\t'
tail -n+2 gencode.map.hgncExt | awk -F"\t" '{if($3!="" && (($3==$5 && $2!=$6) || ($3==$7 && $2!=$8)))print}' 
echo "-------------------------"

## Table 5
echo "- Current approved HGNC genes have Gencode in the HGNC database but Gencode failed to report this asscoiation"
echo "HGNC_id HGNC_Symbol HGNC-Entrez_id Gencode_id Gencode_Symbol" | tr ' ' '\t'
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if(!$3)a[$1]=$2;next}{if(a[$4])print $1,$2,$3,$4,a[$4]}' <(tail -n+2 gencode.noGeneVer.map) <(tail -n+2 hgnc_approved.map)
echo "-------------------------"


## Gencode and NCBI
## Table 1
echo "One Entrez id might have several corresponding Gencode IDs in both NCBI and Gencode annotations"
echo "a) This is in NCBI "
head -n1 Homo_sapiens.gene_info.map
grep "ENSG.*ENSG" Homo_sapiens.gene_info.map 

## Table 2 
echo "b) This is in Gencode"
tail -n+2 gencode.gene_ann.map | grep -v "_PAR_Y" | awk -F"\t" '{if($4!="")print $4}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' > dup.Enterz_ids
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$4]==1)print $1,$2,$4,$3}' dup.Enterz_ids gencode.gene_ann.map | sort -k3,3nr  ## There are 45 sets (44 pairs and an additional set of 3) of Gencode IDs where each set maps to one Entrez IDs. In most sets, the genes have different symbols and HGNC IDs

##################################################################################################################
##################################################################################################################
##################################################################################################################
#### 4. Ambiguity of gene symbols across multiple species

# Download HGNC dataset from the "Custom downloads" page
# Select these columns: HGNC ID Approved symbol Status  Mouse genome database ID        Mouse genome database ID(supplied by MGI)       Rat genome database ID(supplied by RGD)
# Submit and save output webpage as custom_ortho.txt


## a) analysis for the column of Mouse genome database ID(supplied by MGI)
st=$(head -n1 custom_ortho.txt | awk -v col="Status" -F"\t" '{for(i=1; i<=NF; i++){if($i==col)print i}}');
MGD=$(head -n1 custom_ortho.txt | awk -v col="Mouse genome database ID" -F"\t" '{for(i=1; i<=NF; i++){if($i==col)print i}}');
MGI=$(head -n1 custom_ortho.txt | awk -v col="Mouse genome database ID(supplied by MGI)" -F"\t" '{for(i=1; i<=NF; i++){if($i==col)print i}}');
tail -n+2 custom_ortho.txt | awk -v st=$st -v MGD=$MGD -v MGI=$MGI 'BEGIN{FS=OFS="\t"}{if($st=="Approved" && ($MGD!="" || $MGI!=""))print $1,$MGD,$MGI}' > HGNC-to-Mouse.hgnc-custom

tot=$(tail -n+2 custom_ortho.txt | awk -v st=$st -F"\t" '{if($st=="Approved")print $1}' | wc -l)
hum=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | wc -l) ## 18206
mse=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | sort | uniq | wc -l) ## 17976
echo "Among the $tot approved gene IDs in HGNC on $(date '+%b %d, %Y'), $hum genes in HGNC have $mse corresponding genes in the Mouse genome database (MGI) as supplied by MGI."

echo "The relationship is not always one-to-one." 
m_to_one=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)a+=$1}END{print a}') ## 413 
#cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $0}' | wc -l ## 183
echo "There are $m_to_one HGNC genes where each 2 or more genes share the same MGI ID. For example"
topX=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | sort | uniq -c | sort -k1,1nr | head -n1 | awk '{print $2}')
grep "$topX" HGNC-to-Mouse.hgnc-custom | awk -v top=$topX -F"\t" 'BEGIN{FS="\t";output=top"\t"}{output=output" "$1", "}END{print output}'

one_to_m=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | grep "," | wc -l) ## 416
echo "Also, there are another $one_to_m HGNC genes that could be assigned to several MGI IDs. for example"
topY=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $1,$3}' | awk -F"," '{print NF}' | sort -nr | head -n1) ## 15
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $1,$3}' | awk -v top=$topY -F"," '{if(NF==top) print $0}' 


echo "In addition to the Mouse genome database IDs supplied by MGI, HGNC has its own curated Mouse IDs."
hum2=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="")print $2}' | wc -l) ## 17692
mse2=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="")print $2}' | sort | uniq | wc -l) ## 17561
echo "There are $hum2 approved gene IDs in HGNC assigned $mse2 curated Mouse IDs"
iden=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="" && $2==$3)print $0}' | wc -l) ## 17515 
echo "Among these $hum2 relationships, $iden are identical to those provided by MGI "
curated_only=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="" && $3=="")print $0}' | wc -l) ## 3
supplied_only=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2=="" && $3!="")print $0}' | wc -l) ## 515
diff=$(cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="" && $3!="" && $2!=$3)print $0}' | wc -l) ## 174 
echo "But there are $diff IDs different, $curated_only genes with only curated IDs, and $supplied_only genes with only IDs supplied by MGI"
