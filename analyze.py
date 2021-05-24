"""
tr1 | tr1(jac, cont) | tr2(jac, cont)
geneX | geneY(80, 100) | geneX(100, 80), geneY(80, 100), .... 
geneY | geneX(80, 100) | geneX(100, 80), geneY(80, 100), ....
"""

tr1_tr1 = dict()
tr1_tr2 = dict()

tr2_tr2 = dict()
tr2_tr1 = dict()

tr1_prefix = "ver28_GRCh37."
tr2_prefix = "ver38_GRCh38."

def keepTheGreatest_local(current_gene, new_gene):
    # print(f"checking {current_gene} VS. {new_gene}")
    newGeneName, newGeneJacard, newGeneCont = new_gene
    currentGeneName, currentGeneJacard, currentGeneCont = current_gene
    newGeneScore = newGeneJacard + newGeneCont
    currentGeneScore = currentGeneJacard + currentGeneCont
    return new_gene if newGeneScore > currentGeneScore else current_gene

def keepTheGreatest_across(current_list, new_gene):
    newGeneName, newGeneJacard, newGeneCont = new_gene
    newGeneScore = newGeneJacard + newGeneCont
    maxScoreFound = float()
    if not len(current_list):
        return []

    for currentGene in current_list:
        currentGeneName, currentGeneJacard, currentGeneCont = currentGene
        currentGeneScore = currentGeneJacard + currentGeneCont
        maxScoreFound = max(currentGeneScore, maxScoreFound)
    
    if newGeneScore == maxScoreFound:
        current_list.append(new_gene)
    elif newGeneScore > maxScoreFound:
        current_list = [new_gene]
    
    return current_list
        
        

    

file1 = "sample.csv"
file2 = "idx_merged_relations.csv"

with open(file2) as READER:
    next(READER)
    for line in READER:
        line = line.split('|')
        containment, jaccard = float(line[4]), float(line[5])
        gene1, gene2 = line[1], line[6]
        gene1_full = (gene1, jaccard, containment)
        gene2_full = (gene2, jaccard, containment)


        # gene1 & gene2 belongs to tr1
        if tr1_prefix in gene1 and tr1_prefix in gene2:
            if gene1 not in tr1_tr1:
                tr1_tr1[gene1] = gene2_full
            else:
                tr1_tr1[gene1] = keepTheGreatest_local(tr1_tr1[gene1], gene2_full)
                continue
        
        # gene1 & gene2 belongs to tr2
        elif tr2_prefix in gene1 and tr2_prefix in gene2:
            if gene1 not in tr2_tr2:
                tr2_tr2[gene1] = gene2_full
            else:
                tr2_tr2[gene1] = keepTheGreatest_local(tr2_tr2[gene1], gene2_full)
        
        
        # gene1 ∈ tr1, and gene2 ∈ tr2
        elif tr1_prefix in gene1 and tr2_prefix in gene2:
            if gene1 not in tr1_tr2:
                tr1_tr2[gene1] = [gene2_full]
            else:
                tr1_tr2[gene1] = keepTheGreatest_across(tr1_tr2[gene1], gene2_full)

        # gene1 ∈ tr2, and gene2 ∈ tr1
        elif tr2_prefix in gene1 and tr1_prefix in gene2:
            print("OKKKKK")
            if gene1 not in tr2_tr1:
                tr2_tr1[gene1] = [gene2_full]
            else:
                tr2_tr1[gene1] = keepTheGreatest_across(tr2_tr1[gene1], gene2_full)

def getResult_tr1_tr2(gene):
    if gene in tr1_tr2:
        return tr1_tr2[gene]
    else:
        return "-"

def getResult_tr2_tr1(gene):
    if gene in tr2_tr1:
        return tr2_tr1[gene]
    else:
        return "-"

_del = '\t'

with open("results_tr1.tsv", 'w') as R:
    R.write(f"tr1{_del}tr1_local{_del}tr2_across\n")
    for tr1_gene1, tr1_gene2 in tr1_tr1.items():
        R.write(f"{tr1_gene1}{_del}{tr1_gene2}{_del}{getResult_tr1_tr2(tr1_gene1)}\n")

with open("results_tr2.tsv", 'w') as R:
    R.write(f"tr2{_del}tr2_local{_del}tr1_across\n")
    for tr2_gene1, tr2_gene2 in tr2_tr2.items():
        R.write(f"{tr2_gene1}{_del}{tr2_gene2}{_del}{getResult_tr2_tr1(tr1_gene1)}\n")


# print("\nTR1_TR1\n")
# print(list(tr1_tr1.keys())[0], tr1_tr1[list(tr1_tr1.keys())[0]])
# print("\nTR2_TR2\n")
# print(list(tr2_tr2.keys())[0], tr2_tr2[list(tr2_tr2.keys())[0]])
# print("\nTR1_TR2\n")
# print(list(tr1_tr2.keys())[0], tr1_tr2[list(tr1_tr2.keys())[0]])
# print("\nTR2_TR1\n")
# print(list(tr2_tr1.keys())[0], tr2_tr1[list(tr2_tr1.keys())[0]])