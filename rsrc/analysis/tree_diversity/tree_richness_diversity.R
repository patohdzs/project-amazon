
library(raster)
library(sp)

# Data Source: https://figshare.com/articles/dataset/_b_Mapping_density_diversity_and_species-richness_of_the_Amazon_tree_flora_b_/22786166

# Tree species-richness (species/ha) Plot
tree_richness <- raster("data/raw/tree_diversity/TreeRichness_ha.asc")

tree_richness[tree_richness == -9999] <- NA
png("plots/tree_diversity/tree_richness_plot.png", width = 800, height = 600)

plot(tree_richness, col = terrain.colors(100), main = "Tree Species-Richness (species/ha) Amazon Rainforest", axes = TRUE)

dev.off()

# Tree alpha-diversity (Fisher's alpha) Plot
# S=a*ln(1+n/a) where S is number of taxa, n is number of individuals and a is the Fisher's alpha
tree_richness <- raster("data/raw/tree_diversity/TreeDiversity.asc")

tree_richness[tree_richness == -9999] <- NA
png("plots/tree_diversity/tree_diversity_plot.png", width = 800, height = 600)

plot(tree_richness, col = terrain.colors(100), main = "Tree alpha-diversity (Fisher's alpha) in Amazon Rainforest", axes = TRUE)

dev.off()
