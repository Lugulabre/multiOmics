setwd("/Users/MAEL/Documents/M2_BI/Genomique/omiques_floobits/droit_joly/ParisDiderot_202010")

require(MASS)
library(mixOmics)

# Préparation des données
mirna = read.csv("mirna.csv")
row.names(mirna) = mirna[,1]
mirna = mirna[,-1]
mrna = read.csv("mrna.csv")
row.names(mrna) = mrna[,1]
mrna = mrna[,-1]
protein = read.csv("protein.csv")
row.names(protein) = protein[,1]
protein = protein[,-1]
sample_group = read.csv("sample_group.csv")

# Fonction de coefficient de variation
coeff_var = function(vec_val){
  return(sd(vec_val)/mean(vec_val))
}

truehist(apply(mirna, 2, coeff_var), xlab = "Coefficient variation mirna")
truehist(apply(mrna, 2, coeff_var), xlab = "Coefficient variation mrna")
truehist(apply(protein, 2, coeff_var), xlab = "Coefficient variation protein")

# Filtrage des données
new_mirna = mirna[,(which(abs(apply(mirna, 2, coeff_var)) > 0.15))]
new_mrna = mrna[,(which(abs(apply(mrna, 2, coeff_var)) > 0.15))]
new_protein = protein[,(which(abs(apply(protein, 2, coeff_var)) > 0.15))]

# Ttrouver le gène avec la plus grand coefficient de variation
which(max(apply(new_mrna, 2, coeff_var)) == apply(new_mrna, 2, coeff_var))
# Gène avec le plus grand coefficient de variance : PLCD4
# Phospholipase C enzymes play a critical role in many cellular 
# processes by hydrolyzing phosphatidylinositol 4,5-bisphosphate 
# into two intracellular second messengers, inositol 1,4,5-trisphosphate
# and diacylglycerol. Expression of this gene may be a marker for cancer.
# Chromosome 2
# Longueur séquence 30 749 nucléotides
# Longueur acides aminés de la protéine Q9BRC7 762

# PCA
pca.mrna = pca(new_mrna, ncomp = 30, center = TRUE, scale = TRUE)
pca.mrna$cum.var
plotVar(pca.mrna, var.names = FALSE)
vec_col = as.numeric(as.factor(sample_group$Y))
plotIndiv(pca.mrna, col = vec_col+1)
selectVar(pca.mrna, comp = 1)$value

# SPCA
spca.mrna = spca(new_mrna, ncomp = 3, center = TRUE, scale = TRUE,
                    keepX = c(10, 5, 15))
plotVar(spca.mrna)
plotIndiv(spca.mrna, col = vec_col+1)

# Gènes sélectionnés
selectVar(spca.mrna, comp = 1)$value
selectVar(spca.mrna, comp = 2)$value


# PLS-DA
plsda.mrna = plsda(new_mrna, sample_group$Y, ncomp = 2)
plotIndiv(plsda.mrna)
plotVar(plsda.mrna)
splsda.mrna = splsda(new_mrna, sample_group$Y, ncomp = 2,
                     keepX = c(10, 10))
plotIndiv(splsda.mrna)
plotVar(splsda.mrna)
selectVar(splsda.mrna, comp = 1)$value
selectVar(splsda.mrna, comp = 2)$value


#Block PLS-DA
list.of.all = list(mrna = new_mrna, prot = new_protein, mirna = new_mirna)
block.splsda.all = block.splsda(list.of.all,
                                Y = sample_group$Y, ncomp = 5)
plotIndiv(block.splsda.all)
plotVar(block.splsda.all)
perf.splsda = perf(block.splsda.all)

#Block PLS-DA avec keepX
list.keepX = list(mrna = c(15,10), prot = c(15,10), mirna = c(15,10))
block.splsda.keepX = block.splsda(list.of.all,
                                Y = sample_group$Y, ncomp = 2,
                                keepX = list.keepX)
plotIndiv(block.splsda.keepX)
plotVar(block.splsda.keepX)


# Circosplot
circosPlot(block.splsda.keepX, cutoff = 0.5)

# dev.off()
