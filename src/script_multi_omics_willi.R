#setwd("/Users/MAEL/Documents/M2_BI/Genomique/omiques_floobits/droit_joly/ParisDiderot_202010")

require(MASS)
library(mixOmics)

                                        # partie I
# Préparation des données
mirna = read.csv("../ParisDiderot_202010/mirna.csv")
row.names(mirna) = mirna[,1]
mirna = mirna[,-1]
mrna = read.csv("../ParisDiderot_202010/mrna.csv")
row.names(mrna) = mrna[,1]
mrna = mrna[,-1]
protein = read.csv("../ParisDiderot_202010/protein.csv")
row.names(protein) = protein[,1]
protein = protein[,-1]
sample_group = read.csv("../ParisDiderot_202010/sample_group.csv")

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

## Partie II
## PCA
pca.mrna = pca(new_mrna, ncomp = 30, center = TRUE, scale = TRUE)
pca.mrna$cum.var
plotVar(pca.mrna, var.names = FALSE)
vec_col =  factor(sample_group$Y)
plotIndiv(pca.mrna, col.per.group = sample_group$Y)



## SPCA
spca.results = spca(new_mrna, ncomp = 3, center = TRUE, scale = TRUE,
                    keepX = c(10, 5, 15))

plotVar(spca.results)

## Gènes sélectionnés

selectVar(spca.results, comp = 1)$value
selectVar(spca.results, comp = 2)$value

## 2 projection latent structures (PLS)
pls_result = pls(X = new_mrna,Y = new_protein,ncomp = 3)

pls_result$loadings
pls_result$variates
pls_result$names
pls_result$loadings$X[,1]

par(mfrow = c(2,3))
############################################################

plot(pls_result$variates$X[,2],pls_result$variates$X[,3],
     col =sample_group$Y ,
     ylab = "composante 3",xlab = "composante 2",
     main = "comparaisons des composante 1 et 2 avec les groupe\n d'echantillons afficher en couleur.")
legend(x = 11, y = 8, names(table(sample_group$Y)), col = sample_group$X)
col = as.numeric(sample_group$Y)

plotArrow(pls_result,comp = c(1,3), X.label =  "comp1", Y.label = "comp3",col = as.numeric(sample_group$Y), title = "Arrow plot de composante 1 et 3")
l  = list()
l$new_mrna
spls_resulte = spls(X = new_mrna ,Y = new_protein ,ncomp =3,keepX =c(10,5,1) ,keepY =c(9,5,1) )
spls_resulte
cim(spls_resulte)
network(spls_resulte,cutoff = 0.65)
###########################################################
## PLS-DA
pls.mrna = plsda(new_mrna, sample_group$Y, ncomp = 2)
plotIndiv(pls.mrna)
plotVar(pls.mrna)

pls.protein = plsda(new_protein, sample_group$Y, ncomp = 3)
plotIndiv(pls.protein, comp = c(2,3), X.label = "comp 2",
          Y.label = "comp 3")
plotArrow(pls.protein, comp = c(1,3), X.label = "comp 1",
          Y.label = "comp 3")
