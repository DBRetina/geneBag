#### 1. Ambiguity of human gene symbols within the same database:
#### Multiple genes in the same database might have the same gene symbol. 
#### The problem is much worse when we consider the alternative gene alias. 
#### Moreover, shuffling of symbols between genes add another hard-to-fix source of errors for meta-analysis studies. 
#### Some gene catalogs report the previous gene symbols but without clear versioning scheme, it is hard to resolve ambiguities. 
##########################

#### HGNC 
echo "## Explore the gene symbols ambiguity in HGNC"
echo "============================================="
#### Download HGNC dataset (link in the "Statistics & download files" page)
if [ ! -f hgnc_complete_set.txt ];then 
  wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt;fi
if [ ! -f withdrawn.txt ];then
  wget http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/withdrawn.txt;fi


#### Generate maps of gene IDs to symbols
## Generate a map of current gene IDs to official symbols
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2}' > hgnc.ID_to_Current
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,"<"$2">"}' >> hgnc.ID_to_Current 

## Generate a map of current gene IDs to aliases 
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$9}' > hgnc.ID_to_EachAlias
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($9!="")print $1,$9}' | sed 's/"//g' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> hgnc.ID_to_EachAlias  

## Generate a map of current gene IDs to previous symbols
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$11}' > hgnc.ID_to_EachPrev
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($11!="")print $1,$11}' | sed 's/"//g' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> hgnc.ID_to_EachPrev  

## Generate a map for genes IDs withdrawn without approved replacement to their symbols
## Known limitation: We are considering genes replaced by withdrawn genes to be discontinued. There is a possiblity that this new withdrawn gene is also replaced by new approved gene 
##                   but I do not see any example in the current DB (if a gene is replaced, the new genes (column 4) are either "Approved" or "Entry Withdrawn" but not "Merged/Split") 
head -n1 withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$3}' > hgnc.ID_to_discontinued
cat withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{if($2=="Entry Withdrawn")print $1,"<"$3">"}' >> hgnc.ID_to_discontinued;
cat withdrawn.txt | grep -v "Approved" | awk 'BEGIN{FS=OFS="\t";}{if($2=="Merged/Split")print $1,"<"$3">"}'  >> hgnc.ID_to_discontinued  

## Generate a map for genes IDs withdrawn but replaced by new approved IDs to their symbols
head -n1 withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$3}' > hgnc.ID_to_replaced 
cat withdrawn.txt | grep "Approved" | awk 'BEGIN{FS=OFS="\t";}{if($2=="Merged/Split")print $1,"<"$3">"}' >> hgnc.ID_to_replaced 


#### Generate lists of gene symbols
## create list of all current symbols
tail -n+2 hgnc.ID_to_Current | awk -F"\t" '{print $2}' | sort > hgnc.Symbols  
## create list of all alias symbols
tail -n+2 hgnc.ID_to_EachAlias | awk -F "\t" '{print $2}' | sort > hgnc.Alias 
## create list of all previous symbols
tail -n+2 hgnc.ID_to_EachPrev | awk -F "\t" '{print $2}' | sort > hgnc.Prev
## create list of withdrawn symbols without approved replacement 
tail -n+2 hgnc.ID_to_discontinued  | awk -F "\t" '{print $2}' | sort > hgnc.discontinued 
## create list of withdrawn symbols with approved replacement 
tail -n+2 hgnc.ID_to_replaced | awk -F "\t" '{print $2}' | sort > hgnc.replaced


#### Basic check
echo "1. Basic check of HGNC Ids:"
echo "---------------------------"
app_IDs=$(tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS="\t";}{print $1}' | sort | uniq | wc -l)
app_tot=$(tail -n+2 hgnc_complete_set.txt | wc -l)
if (( app_tot != app_IDs));then echo "WARNING: The are $app_tot IDs in the approved HGNC but only $app_IDs are uniq. The master HGNC dataset has duplicate IDs!!";
else echo "OK! The master HGNC dataset has no duplicate IDs";fi
nonApp=$(tail -n+2 hgnc_complete_set.txt | awk -F"\t" '{if($6!="Approved")print}' | wc -l)
if (( nonApp > 0));then echo "WARNING: There are $nonApp non-approved records in the master HGNC dataset";
else echo "OK! All records are approved in the master HGNC dataset";fi

wd_IDs=$(tail -n+2 withdrawn.txt | awk 'BEGIN{FS="\t";}{print $1}' | sort | uniq | wc -l)
wd_tot=$(tail -n+2 withdrawn.txt | wc -l)
if (( wd_tot != wd_IDs));then echo "WARNING: The are $wd_tot IDs in the withdrawn HGNC dataset but only $wd_IDs are uniq. The withdrawn HGNC dataset has duplicate IDs!!";
else echo "OK! The withdrawn HGNC dataset has no duplicate IDs";fi
echo "Basic check is done"
echo "-------------------"


#### Basic statistics
echo "2. Basic statistics"
echo "-------------------"
echo "HGNC approved symbols     = " $(cat hgnc.Symbols | wc -l)        ## 42698
echo "HGNC alias symbols        = " $(cat hgnc.Alias | wc -l)          ## 42347
echo "HGNC previous symbols     = " $(cat hgnc.Prev | wc -l)           ## 15180
echo "HGNC discontinued symbols = " $(cat hgnc.discontinued | wc -l)   ## 1826
echo "HGNC replaced symbols     = " $(cat hgnc.replaced | wc -l)       ## 3314
echo "-------------------"


#### check for repeated symbols within the same gene record 
echo "3. check for repeated symbols within the same gene record:"
echo "----------------------------------------------------------"
cat hgnc.ID_to_Current hgnc.ID_to_EachAlias | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNCsame.01.Alias_symbols_matching_current_symbols
cat hgnc.ID_to_Current hgnc.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNCsame.02.Previous_symbols_matching_current_symbols
cat hgnc.ID_to_EachAlias hgnc.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > HGNCsame.03.Previous_symbols_matching_alias_symbols

wc -l HGNCsame.*_matching_*_symbols 
echo "Examples for these records (if any):"
cat HGNCsame.*_matching_*_symbols | head
echo "-------------------"

echo "4. Tracking of replaced HGNC IDs"
echo "--------------------------------"
cat withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{if($2=="Merged/Split"){split($4,a,", ");for(i in a){split(a[i],b,"|");if(b[3]=="Approved")print b[1],"<"$3">";}}}' | sort | uniq > wdHGNC.ID_to_EachPrev ## the terminal sort|uniq is to overcome a bug in the current HGNC report where the HGNC:35188 shows up twice in the same record of the withdrawn HGNC:21128 ID
echo "No. of withdrawn UIDs replaced by new approved UIDs = "$(cat wdHGNC.ID_to_EachPrev | wc -l)

comm -12 <(sort hgnc.ID_to_Current) <(sort wdHGNC.ID_to_EachPrev) > wdHGNC.ID_to_EachPrev_AsCurrent
echo "No. of withdrawn symbols used as current official symbols = "$(cat wdHGNC.ID_to_EachPrev_AsCurrent | wc -l)

comm -12 <(sort hgnc.ID_to_EachAlias) <(sort wdHGNC.ID_to_EachPrev) > wdHGNC.ID_to_EachPrev_AsAlias
echo "No. of withdrawn symbols used as current alias symbols = "$(cat wdHGNC.ID_to_EachPrev_AsAlias | wc -l)

comm -12 <(sort hgnc.ID_to_EachPrev) <(sort wdHGNC.ID_to_EachPrev) > wdHGNC.ID_to_EachPrev_AsPrev
echo "No. of withdrawn symbols reported in previous symbols of current genes = "$(cat wdHGNC.ID_to_EachPrev_AsPrev | wc -l)

comm -23 <(sort wdHGNC.ID_to_EachPrev) <(cat wdHGNC.ID_to_EachPrev_AsCurrent wdHGNC.ID_to_EachPrev_AsAlias wdHGNC.ID_to_EachPrev_AsPrev | sort | uniq) > PrevSymbols.missing_From_Current
echo "No. of missing withdrawn symbols that do not show up in the new approved gene records = "$(cat PrevSymbols.missing_From_Current | wc -l)

comm -23 <(sort hgnc.ID_to_EachPrev) <(sort wdHGNC.ID_to_EachPrev) > PrevSymbols.missing_From_Withdrawn
echo "No. of previous symbols of current genes that do not show up in the withdrawn symbols = "$(cat PrevSymbols.missing_From_Withdrawn | wc -l)
echo "----------------------"


#### Gene ambiguity
## Identify genes with ambiguous alias or previous symbol (the alias or previous symbol is ambiguous if it matches another alias, previous or current gene symbol)
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

echo "5. Gene ambiguity venn diagram"
echo "------------------------------"
wc -l HGNC.*_matching_*_symbols 

cat hgnc.{Symbols,Alias,Prev,discontinued,replaced} | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -nr  > HGNC.ambiguous_freq.txt
cat HGNC.ambiguous_freq.txt | awk -F"<" '{print "<"$2}' | grep -Fwf - <(cat hgnc.ID_to_{Current,EachAlias,EachPrev,discontinued,replaced}) | sort -t$'\t' -k2,2 >  hgnc.ambiguous.temp 
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$6,$9,$11,"new_symbol_ifReplaced"}' > hgnc.complete_withdrawn.temp
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$6,$9,$11,"-"}' >> hgnc.complete_withdrawn.temp
tail -n+2 withdrawn.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$3,$2,"-","-",$4}' >> hgnc.complete_withdrawn.temp
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print "<Ambiguous_Symbol>",$1,$2,$6,$9,$11,"new_symbol_ifReplaced"}' > HGNC.ambiguous.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{print $2,a[$1]}' hgnc.complete_withdrawn.temp hgnc.ambiguous.temp | sort | uniq >> HGNC.ambiguous.tab

echo " "
echo "HGNC has "$(cat HGNC.ambiguous_freq.txt | wc -l)" ambigious symbols causing "$(tail -n+2 HGNC.ambiguous.tab | wc -l)" ambigious records."
echo "Here are the most 10 ambiguous symbols and how many time do they show up among all gene symbols:"
head HGNC.ambiguous_freq.txt
echo "-------------------------"

##################################################################################################################
#### NCBI genes
echo "## Explore the gene symbols ambiguity in Entrez gene DB"
echo "======================================================="
#### Download
if [ ! -f Homo_sapiens.gene_info ];then
  wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
  gunzip Homo_sapiens.gene_info.gz
fi
if [ ! -f gene_history ];then
  wget ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz
  gunzip gene_history.gz
fi

## Update Master files (Part 1):
head -n1 gene_history > human_gene_history 
grep -w ^9606 gene_history | sed 's/~withdrawn//' >> human_gene_history ## The sed command is used to remove "~withdrawn" appended - most probably by mistake - to 17 Discontinued_Symbol 
head -n1 human_gene_history | awk 'BEGIN{FS=OFS="\t";}{print $0,"status","current_GeneID"}' > human_gene_history_track
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$2;next;}!/^#/{if($2=="-"){print $0,"<discontinued>","-";}
                                               else{new_id=$2;while(1){ \
                                                 if(a[new_id]=="-"){print $0,"<discontinued>","-";break;} \
                                                 else if(!a[new_id]){print $0,"<replaced>",new_id;break;} \
                                                 else new_id=a[new_id];}}}' human_gene_history human_gene_history >> human_gene_history_track

#### Generate maps of gene IDs to symbols
## Generate a map of current gene IDs to official symbols
head -n1 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$3}' > entrez.ID_to_Current
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,"<"$3">"}' >> entrez.ID_to_Current 

## Generate a map of current gene IDs to aliases
head -n1  Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$5}' > entrez.ID_to_EachAlias
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{if($5!="-")print $2,$5}' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> entrez.ID_to_EachAlias 

## Generate a map of current gene IDs to previous symbols
cat human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($6=="<replaced>")print $7,"<"$4">";}' | sort | uniq > wdEntrez.ID_to_EachPrev
echo "GeneID" "prev_symbol" | tr ' ' '\t' > entrez.ID_to_EachPrev
comm -23 <(sort wdEntrez.ID_to_EachPrev) <(cat entrez.ID_to_Current entrez.ID_to_EachAlias | sort | uniq) >> entrez.ID_to_EachPrev 

## Generate a map for genes IDs withdrawn without approved replacement to their symbols
head -n1 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{print $3,$4}' > entrez.ID_to_discontinued 
cat human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($6=="<discontinued>")print $3,"<"$4">"}' >> entrez.ID_to_discontinued 

## Generate a map for genes IDs withdrawn but replaced by new approved IDs to their symbols
head -n1 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{print $3,$4}' > entrez.ID_to_replaced
cat human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{if($6=="<replaced>")print $3,"<"$4">"}' >> entrez.ID_to_replaced


#### Generate lists of gene symbols
## create list of all current symbols
tail -n+2 entrez.ID_to_Current | awk -F"\t" '{print $2}' | sort > entrez.Symbols 
## create list of all alias symbols
tail -n+2 entrez.ID_to_EachAlias | awk -F "\t" '{print $2}' | sort > entrez.Alias 
## create list of all previous symbols
tail -n+2 entrez.ID_to_EachPrev | awk -F "\t" '{print $2}' | sort > entrez.Prev
## create list of withdrawn symbols without new replacement 
tail -n+2 entrez.ID_to_discontinued  | awk -F "\t" '{print $2}' | sort > entrez.discontinued 
## create list of withdrawn symbols with new replacement 
tail -n+2 entrez.ID_to_replaced | awk -F "\t" '{print $2}' | sort > entrez.replaced


## Update Master files (Part 2):
tail -n+2 entrez.ID_to_EachPrev | awk 'BEGIN{FS="[\t<>]";OFS="\t";S="|"}{if(!a[$1])a[$1]=$3;else a[$1]=a[$1]S$3;}END{for(i in a)print i,a[i]}' > entrez.ID_to_EachPrev_aggSyn
awk 'BEGIN{FS=OFS="\t";a["GeneID"]="PrevSymbol";}FNR==NR{a[$1]=$2;next;}{if(a[$2])print $0,a[$2];else print $0,"-";}'  entrez.ID_to_EachPrev_aggSyn Homo_sapiens.gene_info > Homo_sapiens.gene_info.wPrev


echo "1. Basic check of Entrez Ids:"
echo "-----------------------------"
ncbi_IDs=$(tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{print $2}' | sort | uniq | wc -l) ## 61622
ncbi_tot=$(tail -n+2 Homo_sapiens.gene_info | wc -l)
if ((ncbi_tot != ncbi_IDs));then echo "WARNING: The are $ncbi_tot IDs for NCBI genes but only $ncbi_IDs are uniq. The master NCBI dataset has duplicate IDs!!";
else echo "OK! The master NCBI dataset has no duplicate IDs";fi

ncbiHis_IDs=$(tail -n+2 human_gene_history | awk 'BEGIN{FS="\t";}{print $3}' | sort | uniq | wc -l)
ncbiHis_tot=$(tail -n+2 human_gene_history | wc -l)
if ((ncbiHis_tot != ncbiHis_IDs));then echo "WARNING: The are $ncbiHis_tot IDs for discontinued NCBI genes but only $ncbiHis_IDs are uniq. There are discontinued NCBI dataset has duplicate IDs!!";
else echo "OK! The withdrawn NBCI dataset has no duplicate IDs";fi
echo "Basic check is done"
echo "-------------------"


#### Basic statistics
echo "2. Basic statistics"
echo "-------------------"
echo "Entrez official symbols     = " $(cat entrez.Symbols | wc -l)        ## 63881
echo "Entrez alias symbols        = " $(cat entrez.Alias | wc -l)          ## 71213
echo "Entrez previous symbols     = " $(cat entrez.Prev | wc -l)           ## 15723
echo "Entrez discontinued symbols = " $(cat entrez.discontinued | wc -l)   ## 142208
echo "Entrez replaced symbols     = " $(cat entrez.replaced | wc -l)       ## 21743
echo "-------------------"


## check for repeated symbols in the same gene record
echo "3. check for repeated symbols within the same gene record:"
echo "----------------------------------------------------------"
cat entrez.ID_to_Current entrez.ID_to_EachAlias | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrezsame.01.Alias_symbols_matching_current_symbols
cat entrez.ID_to_Current entrez.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrezsame.02.Previous_symbols_matching_current_symbols
cat entrez.ID_to_EachAlias entrez.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrezsame.03.Previous_symbols_matching_alias_symbols

wc -l Entrezsame.*_matching_*_symbols
echo "Examples for these records (if any):"
cat Entrezsame.*_matching_*_symbols | head
echo "-------------------"

echo "4. Tracking of replaced Entrez IDs"
echo "----------------------------------"
echo "No. of withdrawn UIDs replaced by new approved UIDs = "$(cat wdEntrez.ID_to_EachPrev | wc -l)

comm -12 <(sort entrez.ID_to_Current) <(sort wdEntrez.ID_to_EachPrev) > wdEntrez.ID_to_EachPrev_AsCurrent
echo "No. of withdrawn symbols used as current official symbols = "$(cat wdEntrez.ID_to_EachPrev_AsCurrent | wc -l)

comm -12 <(sort entrez.ID_to_EachAlias) <(sort wdEntrez.ID_to_EachPrev) > wdEntrez.ID_to_EachPrev_AsAlias
echo "No. of withdrawn symbols used as current alias symbols = "$(cat wdEntrez.ID_to_EachPrev_AsAlias | wc -l)

comm -23 <(sort wdEntrez.ID_to_EachPrev) <(cat wdEntrez.ID_to_EachPrev_AsCurrent wdEntrez.ID_to_EachPrev_AsAlias | sort | uniq) > PrevSymbols.missing_From_Current
echo "No. of missing withdrawn symbols that do not show up in the new approved gene records = "$(cat PrevSymbols.missing_From_Current | wc -l)

grep -v "<LOC.*>" PrevSymbols.missing_From_Current | grep -v "<FLJ.*>" | grep -v "<KIAA.*>" | grep -v "<MGC.*>" | grep -v "<DKFZ[Pp].*>" | grep -vi "<c.*orf.*>" | grep -v "<DJ.*\..*>" | grep -v "<HSA.*>" | grep -v "<PRO.*>" | grep -v "<RP.*-.*\..*>" | grep -v "<AC.*\..*>" | grep -vi "<none>" | grep -vi "<null>" > PrevSymbols.missing_From_Current.knownIDs
echo "No. of missing withdrawn symbols after exclusion of LOC* IDs  = "$(cat PrevSymbols.missing_From_Current.knownIDs | wc -l)
echo "-------------------------"



#### Gene ambiguity
## Identify genes with ambiguous alias (the alias is ambiguous if it matches another alias, previous or current gene symbol)
## Gene ambiguity venn diagram
cat entrez.Symbols | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.01.Current_symbols_matching_other_current_symbols
cat entrez.Alias | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.02.Alias_symbols_matching_other_alias_symbols
comm -12 <(cat entrez.Alias | uniq) <(cat entrez.Symbols | uniq) > Entrez.03.Alias_symbols_matching_current_symbols
cat entrez.Prev | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.04.Previous_symbols_matching_other_previous_symbols
comm -12 <(cat entrez.Prev | uniq) <(cat entrez.Symbols | uniq) > Entrez.05.Previous_symbols_matching_current_symbols
comm -12 <(cat entrez.Prev | uniq) <(cat entrez.Alias | uniq)  > Entrez.06.Previous_symbols_matching_alias_symbols
cat entrez.discontinued | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.07.Discontinued_symbols_matching_other_discontinued_symbols
comm -12 <(cat entrez.discontinued | uniq) <(cat entrez.Symbols | uniq) > Entrez.08.Discontinued_symbols_matching_current_symbols
comm -12 <(cat entrez.discontinued | uniq) <(cat entrez.Alias | uniq) > Entrez.09.Discontinued_symbols_matching_alias_symbols
comm -12 <(cat entrez.discontinued | uniq) <(cat entrez.Prev | uniq) > Entrez.10.Discontinued_symbols_matching_previous_symbols
cat entrez.replaced | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Entrez.11.Replaced_symbols_matching_other_replaced_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.Symbols | uniq) > Entrez.12.Replaced_symbols_matching_current_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.Alias | uniq) > Entrez.13.Replaced_symbols_matching_alias_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.Prev | uniq) > Entrez.14.Replaced_symbols_matching_previous_symbols
comm -12 <(cat entrez.replaced | uniq) <(cat entrez.discontinued | uniq) > Entrez.15.Replaced_symbols_matching_discontinued_symbols

echo "5. Gene ambiguity venn diagram"
echo "------------------------------"
wc -l Entrez.*_matching_*_symbols

cat entrez.{Symbols,Alias,Prev,discontinued,replaced} | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -nr  > Entrez.ambiguous_freq.txt
cat Entrez.ambiguous_freq.txt | awk -F"<" '{print "<"$2}' | grep -Fwf - <(cat entrez.ID_to_{Current,EachAlias,EachPrev,discontinued,replaced}) | sort -t$'\t' -k2,2 >  entrez.ambiguous.temp 
head -n1 Homo_sapiens.gene_info.wPrev | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,"status",$5,$17,"new_symbol_ifReplaced"}' > entrez.complete_withdrawn.temp
tail -n+2 Homo_sapiens.gene_info.wPrev | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,"Official",$5,$17,"-"}' >> entrez.complete_withdrawn.temp
tail -n+2 human_gene_history_track | awk 'BEGIN{FS=OFS="\t";}{print $3,$4,$6,$7;}' > withdrawn.temp
awk 'BEGIN{FS=OFS="\t";a["-"]="-";}FNR==NR{a[$2]=$3;next;}{print $1,$2,$3,"-","-",a[$4]}' Homo_sapiens.gene_info withdrawn.temp >> entrez.complete_withdrawn.temp
head -n1 Homo_sapiens.gene_info.wPrev | awk 'BEGIN{FS=OFS="\t";}{print "<Ambiguous_Symbol>",$2,$3,"status",$5,$17,"new_symbol_ifReplaced"}' > Entrez.ambiguous.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{print $2,a[$1]}' entrez.complete_withdrawn.temp entrez.ambiguous.temp | sort |  uniq | sort -k1,1 >> Entrez.ambiguous.tab

echo " "
echo "Entrez Gene DB has "$(cat Entrez.ambiguous_freq.txt | wc -l)" ambigious symbols causing "$(tail -n+2 Entrez.ambiguous.tab | wc -l)" ambigious records."
echo "Here are the most 10 ambiguous symbols and how many time do they show up among all gene symbols:"
head Entrez.ambiguous_freq.txt
echo "-------------------------"
##################################################################################################################
#### GENCODE genes
echo "Explore the gene symbols ambiguity in GENCODE/Ensembl gene DB"
echo "============================================================="
## Download
wget -O ens_current.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" completionStamp = "1" >
            
    <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "external_gene_name" />
        <Attribute name = "external_synonym" />
    </Dataset>
</Query>'
completionStamp=$(tail -n1 ens_current.txt)
if [ "$completionStamp" == "[success]" ]; then
  echo "Ensembl biomart query was done successfuly";
  head -n -1 ens_current.txt > ens_current.txt.temp
  mv ens_current.txt.temp ens_current.txt
else echo "Ensembl biomart query is incomplete. Repeat the download step";fi
echo "GeneID" "ensSymbol" "Aliases" | tr ' ' '\t' > ens_current_aggSyn.txt
cat ens_current.txt | awk 'BEGIN{FS=OFS="\t";S="|"}{if(!a[$1FS$2])a[$1FS$2]=$3;else a[$1FS$2]=a[$1FS$2]S$3;}END{for(i in a)print i,a[i]}' >> ens_current_aggSyn.txt

## Find the version of most recent Gencode annotation
curGTF=$(wget -q -O- http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/ | \
         grep -o "gencode.v[0-9]\+.chr_patch_hapl_scaff.annotation.gtf.gz" | head -n1)
prefix="gencode.v";
suffix=".chr_patch_hapl_scaff.annotation.gtf.gz";
vcur=${curGTF#$prefix};
vcur=${vcur%$suffix};
#vLast=$(($vcur - 1))
echo $vcur #$vLast

# Download all assembly report to make assembly map to differentiate 1ry assembly from patches and alternative loci
mkdir -p assemblies_reports && cd assemblies_reports
wget -q -O- https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/ | \
         grep -o "GCA_000001405.[0-9]\+_GRCh3[78].p[0-9]\+" | sort | uniq > assemblies.lst
cat assemblies.lst | while read asm;do
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/${asm}/${asm}_assembly_report.txt
done

for sym in M X Y $(seq 1 22);do echo "chr"$sym"|Primary Assembly" | tr '|' '\t';done > assembly_map
for rep in *_assembly_report.txt;do cat $rep | awk 'BEGIN{FS=OFS="\t";}!/#/{print $5,$8}';done | sort | uniq >> assembly_map
cd ../
  
# Dowanload all previous Gencode GTFs
mkdir -p gencode_gtf && cd gencode_gtf
genFTP="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"
for i in $(seq 20 ${vcur});do echo $i;
  wget -O gencode.v${i}.ALL.GRCh38.gtf.gz $genFTP/release_${i}/gencode.v${i}.chr_patch_hapl_scaff.annotation.gtf.gz 2>/dev/null;done
for i in $(seq 16 19);do echo $i;
  wget -O gencode.v${i}.ALL.GRCh37.gtf.gz $genFTP/release_${i}/gencode.v${i}.chr_patch_hapl_scaff.annotation.gtf.gz 2>/dev/null;done
for i in $(seq 5 15);do echo $i;
  wget -O gencode.v${i}.CHR.GRCh37.gtf.gz $genFTP/release_${i}/gencode.v${i}.annotation.gtf.gz 2>/dev/null;done
wget -O gencode.v4.CHR.GRCh37.gtf.gz $genFTP/release_4/gencode_v4.annotation.GRCh37.gtf.gz 2>/dev/null;
wget -O gencode.v3d.CHR.GRCh37.gtf.gz $genFTP/release_3d/gencode.v3d.gtf.gz 2>/dev/null;
wget -O gencode.v3c.CHR.GRCh37.gtf.gz $genFTP/release_3c/gencode.v3c.annotation.GRCh37.gtf.gz 2>/dev/null;
wget -O gencode.v3b.CHR.GRCh37.gtf.gz $genFTP/release_3b/gencode.v3b.annotation.GRCh37.gtf.gz 2>/dev/null;
wget -O gencode.v2a.CHR.GRCh37.gtf.gz $genFTP/release_2/gencode_data.rel2a.gtf.gz 2>/dev/null;
wget -O gencode.v2.CHR.GRCh37.gtf.gz $genFTP/release_2/gencode_data.rel2.gtf.gz 2>/dev/null;
wget -O gencode.v1.CHR.GRCh37.gtf.gz $genFTP/release_1/gencode_data.rel1.v2.gtf.gz 2>/dev/null;

r=2
for i in 2 2a 3b 3c 3d $(seq 4 ${vcur});do
  gtf=(gencode.v${i}.*.gtf.gz)
  if [ -f "${gtf}" ];then
    echo $gtf
    gunzip ${gtf}
    output=${gtf%.gtf.gz}
    cat ${gtf%.gz} | awk -F"\t" '!/#/{if($3=="gene")print $1":"$4"-"$5";"$1";"$9}' | grep -v "_PAR_Y" | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" -v ann_version=${i} -v rank=${r} 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; name = vars["gene_name"]; type = vars["gene_type"]; hgnc = vars["hgnc_id"]; sub(/\..*/,"",id); print id,name,type,hgnc,$1,$2,ann_version,rank; }' | grep -v "^ENSGR" > $output.genes
    echo "GeneID" "Symbol" "gene_type" "HGNC" "Location" "Assembly_type" "gencode_version" "rank" | tr ' ' '\t' > $output.genes.ann
    awk 'BEGIN{FS=OFS="\t";}FNR==NR{a[$1]=$2;next;}{$6=a[$6];print $0}' ../assemblies_reports/assembly_map $output.genes >> $output.genes.ann
  fi
  ((r=r+1));
done


# special processing for release 1
# note: There are 4262 gene IDs in this release without gene symbols
gtf=gencode.v1.CHR.GRCh37.gtf.gz
gunzip ${gtf}
> v1.start; > v1.end;
cat ${gtf%.gz} | awk -F"\t" '!/#/{if($3=="exon")print $1";"$4";"$5";"$9}' | grep -v "_PAR_Y" | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" -v ann_version=1 -v rank=1 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; name = vars["gene_name"]; type = vars["gene_type"]; hgnc = vars["hgnc_id"]; sub(/\..*/,"",id); print id,name,type,hgnc,$1,$2,ann_version,rank >> "v1.start"; print id,$3 >> "v1.end"; }'
cat v1.start | sort -t$'\t' -k1,1 -k6,6  | sort -t$'\t' -k1,1 -u > v1.start.uq
cat v1.end | sort -t$'\t' -k1,1 -k2,2nr  | sort -t$'\t' -k1,1 -u > v1.end.uq
output=${gtf%.gtf.gz}
paste v1.start.uq v1.end.uq | awk 'BEGIN{FS=OFS="\t";}!/#/{print $1,$2,$3,$4,$5":"$6"-"$10,$5,$7,$8;}' | grep -v "^ENSGR" > $output.genes
echo "GeneID" "Symbol" "gene_type" "HGNC" "Location" "Assembly_type" "gencode_version" "rank" | tr ' ' '\t' > $output.genes.ann
awk 'BEGIN{FS=OFS="\t";}FNR==NR{a[$1]=$2;next;}{$6=a[$6];print $0}' ../assemblies_reports/assembly_map $output.genes >> $output.genes.ann
rm v1.end v1.start v1.end.uq v1.start.uq


cur_Ann=(gencode.v${vcur}.*.genes.ann)
head -n1 $cur_Ann > ../gencode.gene.track
for ann in *.genes.ann;do tail -n+2 $ann;done | sort -t$'\t' -k1,1 -k8,8nr  >> ../gencode.gene.track
cd ../

## generate list of previous symbols for each gene ID (check to exclude entries from v1 with no gene symbols)
tail -n+2 gencode.gene.track | sort -t$'\t' -k1,2 -u | sort -t$'\t' -k1,1 -k8,8nr | awk 'BEGIN{FS=OFS="\t"}{if($2)print $1,$2}' > gencode.gene.uniqIDs
echo "GeneID" "ensSymbol" "Aliases" "PrevSymbols" | tr ' ' '\t' > ens_current_aggSyn_aggPrev.txt
awk 'BEGIN{FS=OFS="\t";S="|"}FNR==NR{if(!a[$1])a[$1]=$2;else a[$1]=a[$1]S$2;next;}{prev=a[$1];n=index(prev,"|");if(n==0)prev="";print $0,substr(prev, n+1);}' gencode.gene.uniqIDs <(tail -n+2 ens_current_aggSyn.txt) >> ens_current_aggSyn_aggPrev.txt

## update Ensembl gene symbols from current gencode annotation && add all annotation fields from the GTF
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;b[$1]=$3 FS $4 FS $5 FS $6 FS $7;next;}{if(a[$1])print $1,$2,a[$1],$3,$4,b[$1];else print $1,$2,"Not_in_Gencode",$3;}' gencode_gtf/$cur_Ann ens_current_aggSyn_aggPrev.txt > ens_current_aggSyn_aggPrev_genAnn.txt

## Add Entrez IDs
if [ ! -f gencode.v${vcur}.metadata.EntrezGene ];then
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${vcur}/gencode.v${vcur}.metadata.EntrezGene.gz
  gunzip gencode.v${vcur}.metadata.EntrezGene.gz
fi

cur_gtf=(gencode_gtf/gencode.v${vcur}.*.gtf)
cat $cur_gtf | awk -F"\t" '!/#/{if($3=="transcript")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; trans = vars["transcript_id"]; print id,trans; }' > gencode.trans_to_gene.map
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{ print $0, a[$2]}' gencode.v${vcur}.metadata.EntrezGene gencode.trans_to_gene.map > gencode.trans_to_entrez.map
echo "GeneID EntrezGene" | tr ' ' '\t' > gencode.gene_to_entrez.map
cat gencode.trans_to_entrez.map | awk 'BEGIN{FS="[.\t]";OFS="\t"}{print $1,$5}' | sort | uniq >> gencode.gene_to_entrez.map
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next;}{$7=a[$1] FS $7;print $0;}' gencode.gene_to_entrez.map ens_current_aggSyn_aggPrev_genAnn.txt > ens_current_aggSyn_aggPrev_genAnn_dbXrefs.txt
gencode_master="ens_current_aggSyn_aggPrev_genAnn_dbXrefs.txt"

## generate a list of symbols for discontinued gene IDs
head -n1 gencode.gene.track > gencode.gene.discontinued
tail -n+2 gencode.gene.track | grep -v ^OTTHUMG | grep -v -Fwf <(tail -n+2 ens_current_aggSyn.txt | cut -f1) | sort -t$'\t' -k1,2 -u >> gencode.gene.discontinued


#### Generate maps of gene IDs to symbols
## Generate a map of current gene IDs to official symbols
echo "GeneID" "Symbol" | tr ' ' '\t' > gencodePrim.ID_to_Current
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $9=="Primary Assembly")print $1,"<"$3">"}' >> gencodePrim.ID_to_Current ## 60665

## Generate a map of current gene IDs to aliases 
head -n1 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print $1,$4}' > gencodePrim.ID_to_EachAlias
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $4!="" && $10=="Primary Assembly")print $1,$4}' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> gencodePrim.ID_to_EachAlias

## Generate a map of current gene IDs to previous symbols
head -n1 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print $1,$5}' > gencodePrim.ID_to_EachPrev
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $5!="" && $10=="Primary Assembly")print $1,$5}' | awk 'BEGIN{FS="\t";OFS="\n";}{split($2,a,"|");for(i in a)print $1"\t<"a[i]">";}' >> gencodePrim.ID_to_EachPrev

## Generate a map for genes IDs withdrawn to their symbols
echo "GeneID" "Symbol" | tr ' ' '\t' > gencodePrim.ID_to_discontinued
tail -n+2 gencode.gene.discontinued | awk 'BEGIN{FS=OFS="\t";}{if($2!="" && $6=="Primary Assembly")print $1,"<"$2">"}' >> gencodePrim.ID_to_discontinued


#### Generate lists of gene symbols
## create list of all current symbols
tail -n+2 gencodePrim.ID_to_Current | awk -F"\t" '{print $2}' | sort > gencodePrim.Symbols
## create list of all alias symbols
tail -n+2 gencodePrim.ID_to_EachAlias | awk -F "\t" '{print $2}' | sort > gencodePrim.Alias
## create list of all previous symbols
tail -n+2 gencodePrim.ID_to_EachPrev | awk -F "\t" '{print $2}' | sort > gencodePrim.Prev
## create list of withdrawn symbols
tail -n+2 gencodePrim.ID_to_discontinued  | awk -F "\t" '{print $2}' | sort > gencodePrim.discontinued


#### Basic check
echo "1. Basic check of Ensembl and gencode annotations:"
echo "--------------------------------------------------"
## Problem exploration: Difference between Ensembl and gencode annotations
# 1. Ensembl doesn't have symbols while Gencode does
cat ens_current_aggSyn_aggPrev_genAnn.txt | awk 'BEGIN{FS=OFS="\t"}{if($2!=$3 && $3!="Not_in_Gencode")print}' > ens_current_aggSyn_aggPrev_genAnn_unmatchedSym.txt ##  21606
echo "No of IDs where Ensembl doesn't have symbols while Gencode does     = " $(tail -n+2 ens_current_aggSyn_aggPrev_genAnn_unmatchedSym.txt | wc -l)       
# 2. Ensembl IDs not found in Gencode annotation:
# Example "ENSG00000274081" belongs to ALT_REF_LOCI. The Ensembl website states that its 1st transcript is part of the Gencode basic annotation
head -n1 ens_current_aggSyn_aggPrev_genAnn.txt | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$4}' > ens_current_aggSyn_aggPrev_genAnn_NotinGencode.txt
tail -n+2 ens_current_aggSyn_aggPrev_genAnn.txt | awk 'BEGIN{FS=OFS="\t"}{if($3=="Not_in_Gencode")print $1,$2,$4}' >> ens_current_aggSyn_aggPrev_genAnn_NotinGencode.txt ## 123
echo "No of Ensembl IDs not found in Gencode annotation                   = " $(tail -n+2 ens_current_aggSyn_aggPrev_genAnn_NotinGencode.txt | wc -l)       
# 3. Gencode IDs not found in Ensembl annotation
tail -n+2 gencode_gtf/$cur_Ann | awk -F"\t" '{print $1}' | sort > id.gen  ## 67005
tail -n+2 ens_current_aggSyn_aggPrev.txt | awk -F"\t" '{print $1}' | sort > id.ens     ## 67128
comm -23 id.gen id.ens > id.gen.sp  ## 0  ##i.e. all gencode IDs present in ensembl
missing_ids=$(cat id.gen.sp | wc -l)
if [ $missing_ids -gt 0 ];then echo "These GENCODE IDS are missing from Ensembl annotation"; cat id.gen.sp;else echo "No GENCODE IDS are missing from Ensembl annotation";fi
echo "Basic check is done"
echo "-------------------"


#### Basic statistics
echo "2. Basic statistics"
echo "-------------------"
echo "Gencode official symbols        = " $(cat gencodePrim.Symbols | wc -l)        ## 60664
echo "Gencode(Ensembl) alias symbols  = " $(cat gencodePrim.Alias | wc -l)          ## 54944
echo "Gencode previous symbols        = " $(cat gencodePrim.Prev | wc -l)           ## 80774
echo "Gencode discontinued symbols    = " $(cat gencodePrim.discontinued | wc -l)   ## 38409
echo "-------------------"


#### check for repeated symbols within the same gene record 
echo "3. check for repeated symbols within the same gene record:"
echo "----------------------------------------------------------"
cat gencodePrim.ID_to_Current gencodePrim.ID_to_EachAlias | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencodesame.01.Alias_symbols_matching_current_symbols
cat gencodePrim.ID_to_Current gencodePrim.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencodesame.02.Previous_symbols_matching_current_symbols
cat gencodePrim.ID_to_EachAlias gencodePrim.ID_to_EachPrev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencodesame.03.Previous_symbols_matching_alias_symbols

wc -l Gencodesame.*_matching_*_symbols
echo "Examples for these records (if any):"
cat Gencodesame.*_matching_*_symbols | head
echo "-------------------"

echo "4. Tracking of replaced HGNC IDs"
echo "--------------------------------"
echo "TO BE DONE"
echo "----------------------"


#### Gene ambiguity
## Identify genes with ambiguous alias (the alias is ambiguous if it matches another alias, previous or current gene symbol)
## Gene ambiguity venn diagram
cat gencodePrim.Symbols | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencode.01.Current_symbols_matching_other_current_symbols
cat gencodePrim.Alias | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencode.02.Alias_symbols_matching_other_alias_symbols
comm -12 <(cat gencodePrim.Alias | uniq) <(cat gencodePrim.Symbols | uniq) > Gencode.03.Alias_symbols_matching_current_symbols
cat gencodePrim.Prev | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencode.04.Previous_symbols_matching_other_previous_symbols
comm -12 <(cat gencodePrim.Prev | uniq) <(cat gencodePrim.Symbols | uniq) > Gencode.05.Previous_symbols_matching_current_symbols
comm -12 <(cat gencodePrim.Prev | uniq) <(cat gencodePrim.Alias | uniq)  > Gencode.06.Previous_symbols_matching_alias_symbols
cat gencodePrim.discontinued | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Gencode.07.Discontinued_symbols_matching_other_discontinued_symbols
comm -12 <(cat gencodePrim.discontinued | uniq) <(cat gencodePrim.Symbols | uniq) > Gencode.08.Discontinued_symbols_matching_current_symbols
comm -12 <(cat gencodePrim.discontinued | uniq) <(cat gencodePrim.Alias | uniq) > Gencode.09.Discontinued_symbols_matching_alias_symbols

echo "5. Gene ambiguity venn diagram"
echo "------------------------------"
wc -l Gencode.*_matching_*_symbols

cat gencodePrim.{Symbols,Alias,Prev,discontinued} | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -nr  > Gencode.ambiguous_freq.txt
cat Gencode.ambiguous_freq.txt | awk '{print $2}' | grep -Fwf - <(cat gencodePrim.ID_to_{Current,EachAlias,EachPrev,discontinued}) | sort -t$'\t' -k2,2 >  gencodePrim.ambiguous.temp
cat $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print $1,$3,"status",$4,$5,"new_symbol_ifReplaced"}' > gencodePrim.complete_withdrawn.temp
cat $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $9=="Primary Assembly") print $1,$3,"Current",$4,$5,"-"}' >> gencodePrim.complete_withdrawn.temp
tail -n+2 gencode.gene.discontinued | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,"discontinued","-","-","-";}' >> gencodePrim.complete_withdrawn.temp
head -n1 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print "<Ambiguous_Symbol>",$1,$3,"status",$4,$5,"new_symbol_ifReplaced"}' > Gencode.ambiguous.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{print $2,a[$1]}' gencodePrim.complete_withdrawn.temp gencodePrim.ambiguous.temp | sort | uniq >> Gencode.ambiguous.tab

echo " "
echo "GENCODE has "$(cat Gencode.ambiguous_freq.txt | wc -l)" ambigious symbols causing "$(tail -n+2 Gencode.ambiguous.tab | wc -l)" ambigious records."
echo "Here are the most 10 ambiguous symbols and how many time do they show up among all gene symbols:"
head Gencode.ambiguous_freq.txt
echo "-------------------------"

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
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$19,$20,$9,$11,$4 " - " $5}' > hgnc_approved.map ## hgnc_id symbol  entrez_id   ensembl_gene_id    alias_symbol    prev_symbol    locus_group-locus_type(i.e. gene_type)
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($6=="Approved")print $1,$2,$19,$20,$9,$11,$4 " - " $5}' >> hgnc_approved.map
 
## Stats
tail -n+2 hgnc_approved.map | awk 'BEGIN{FS="\t";OFS="\n";}{a+=1;if($3!="")b+=1;if($4!="")c+=1;}END{print "There are "a" approved HGNC ids.", "They are cross referenced with:", "Entrez IDs = "b, "Ensembl IDs = "c}' 
tail -n+2 hgnc_approved.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "HGNC genes with entrez and ensembl ids =",a}'  ## hgnc_id with entrez and ensembl ids = 39578



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
head -n1 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,"HGNC","Ensembl",$5,$10}' > Homo_sapiens.gene_info.map ## GeneID  Symbol  HGNC    Ensembl     Synonyms    type_of_gene
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{if($6!="-")print $2,$3,$6,$5,$10}' | awk 'BEGIN{FS=OFS="\t"}{ delete vars; split($3,a,"|");for(i in a) { n = index(a[i], ":"); if(n) { x = substr(a[i], n + 1); key = substr(a[i], 1, n - 1); val = substr(a[i], n + 1, length(x)); if(vars[key]=="")vars[key] = val;else vars[key] = vars[key]","val; } } HGNC = vars["HGNC"]; Ensembl = vars["Ensembl"]; print $1,$2,HGNC,Ensembl,$4,$5; }' >> Homo_sapiens.gene_info.map
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{if($6=="-")print $2,$3,"","",$5,$10}' >> Homo_sapiens.gene_info.map ## $6 = dbXrefs

## Stats
tail -n+2 Homo_sapiens.gene_info.map | awk 'BEGIN{FS="\t";OFS="\n";}{a+=1;if($3!="")b+=1;if($4!="")c+=1;}END{print "There are "a" current Entrez ids.", "They are cross referenced with:", "HGNC IDs = "b, "and Ensembl IDs = "c}'
tail -n+2 Homo_sapiens.gene_info.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "Entrez genes with hgnc and ensembl IDs =",a}'  ## entrez_id with hgnc and ensembl ids= 33165



#### Gencode_human
## Explore
# Explore types of genes:
echo "Current Gencode genes can be calssified as:"
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $10=="Primary Assembly")print $6}' | sort | uniq -c

## Generate a map of Gencode gene IDs on the Primary Assembly to dbXrefs
head -n1 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{print $1,$3,$7,$8,$4,$5,$6}' > gencode_primary.map ## GeneID Symbol  EntrezGene   HGNC    Aliases    PrevSymbols    gene_type
tail -n+2 $gencode_master | awk 'BEGIN{FS=OFS="\t";}{if($3!="Not_in_Gencode" && $10=="Primary Assembly")print $1,$3,$7,$8,$4,$5,$6}' >> gencode_primary.map

## Stats
tail -n+2 gencode_primary.map | awk 'BEGIN{FS="\t";OFS="\n";}{a+=1;if($3!="")b+=1;if($4!="")c+=1;}END{print "There are "a" current Ensembl ids", "They are cross referenced with:", "HGNC IDs = "b, "Entrez IDs = "c}'  
tail -n+2 gencode_primary.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "Gencode gene with hgnc and entrez  ids =",a}'  ## gencode_id with hgnc and entrez  ids= 24528
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
