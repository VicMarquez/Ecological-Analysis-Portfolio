rm(list=ls())
objects()

## GENERAR ARBOL FILOGENÉTICO CON U.PHYLOMAKER

#### INSTALACION DE U.PHYLOMAKER #########
devtools::install_github("jinyizju/U.PhyloMaker")

# load the package
library (ape)
library("U.PhyloMaker")

library(phytools)
library (maps)
library(V.PhyloMaker)

### AQUI VEO DONDE DEBEN IR LOS ARCHIVOS: POR DEFAULT "C:/Users/Ramiro/Documents"
getwd()

# input the sample species list, the megatree and genus-family relationship ﬁles ###
sp.list <- read.csv("forrajeras.csv")
megatree <- read.tree("plant_megatree.tre")
gen.list <- read.csv("plant_genus_list.csv") ### ASEGURARSE QUE LA GEN.LIST TENGA DOS COLUMNAS (SPECIES Y FAMILIA)



# generate a phylogeny for the sample species list
result <- phylo.maker(sp.list, megatree, gen.list, nodes.type=1, scenario=3)
write.tree(result$phylo, "output_tree.tre")
write.csv(result$sp.list, "output_splist.csv")

### Con phytools puedo calcular señal filogenética (K de Bloomberg)

remotes::install_github("liamrevell/phytools")

library(phytools)

## load data from Garland et al. (1992)

data(mammal.tree)
data(mammal.data)

## extract characters of interest
ln.bodyMass<-log(setNames(mammal.data$bodyMass,
                          rownames(mammal.data)))
ln.homeRange<-log(setNames(mammal.data$homeRange,
                           rownames(mammal.data)))
## compute phylogenetic signal K
K.bodyMass<-phylosig(mammal.tree,ln.bodyMass,
                     test=TRUE)
print(K.bodyMass)
plot(K.bodyMass)

###PARA DIBUJAR EL �RBOL FILOGENETICO ####

library(ape)
getwd()

tree <- read.tree("forrajeras_tree.newick")
plot(tree)

### Customize the tree visualization: Add titles, labels, or modify the style.
plot(tree, main = "Phylogenetic Tree", edge.width = 1, cex = 0.5)

## Export the tree: Save it as an image or PDF:
pdf("phylogenetic_tree.pdf")
plot(tree)
dev.off()

#### Using the "ggtree Package"
### Install and load the necessary packages:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
library(ggtree)
library(ape) # Also needed to read the tree
library(ggplot2)

tree <- read.tree("forrajeras_tree.newick")

## Visualize the tree:
ggtree(tree) + 
  theme_tree() + 
  ggtitle("Phylogenetic Tree")

## Add annotations and customizations: For example, add node labels or customize branch colors:
ggtree(tree, layout = "circular") + 
  geom_tiplab(size = 3) + 
  ggtitle("Circular Phylogenetic Tree")

## 2. Join Family Data to the Tree
## Use left_join to add the family information to the tree object.
library(ggtree)
library(dplyr)

# Read your tree
# tree <- read.tree("your_tree.newick")

# Add family data to the tree
tree_data <- as_tibble(tree) %>%
  left_join(family_data, by = c("label" = "tip"))

## 3. Visualize with Colors by Family
## Use the aes(color = family) aesthetic in geom_tiplab():
ggtree(tree, layout = "circular") +
  geom_tiplab(aes(color = family), size = 1) +
  ggtitle("Circular Phylogenetic Tree") +
  theme(legend.position = "right") # Optional: Show legend

## 5. If Families Are Already Encoded
##If your tree object already contains family information in its tip labels, you can directly extract and map it:
# Extract family from tip labels (example: "species_Family")
tree$data <- tree$data %>%
  mutate(family = sub(".*_", "", label)) # Extract family info after "_"

# Plot
ggtree(tree, layout = "circular") +
  geom_tiplab(aes(color = family), size = 1) +
  ggtitle("Circular Phylogenetic Tree")

