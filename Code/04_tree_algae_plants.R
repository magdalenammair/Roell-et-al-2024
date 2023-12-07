### Create phylogenetic tree from only plant and algal ott IDs and information on genomic resources and include tip labels ###

# load packages
library(ggtree)
library(ggplot2)
library(rotl)

# import data---------

# phylogenetic tree species list with correct ott numbers:
resolved_species = readRDS("Data/resolved_species_for_trees.rds")
# trait list with adjusted species names (synonyms removed)
trait_data = readRDS("Data/spp_trait_data.rds")


# create tree - algae and plants-------
plant_trait = trait_data[trait_data$tax_category %in% c("algae", "plant", "procaryote"),]

#change sci_name to all small letters
plant_trait$sci_name = tolower(plant_trait$sci_name)

# subset ott list
resolved_plants = resolved_species[resolved_species$search_string %in% plant_trait$sci_name,]
plant_ott = resolved_plants$ott_id
plants_in_tree = is_in_tree(plant_ott)
resolved_plants[!plants_in_tree,] # these three species are not in the tree, mrca will be used in tree

# create tree
planttree = tol_induced_subtree(ott_ids = plant_ott)
length(planttree$tip.label) # 102

# merge name info from resolved species list in order of tip labels ----
plants_merged = data.frame(
  tip.label = planttree$tip.label
)
plants_merged$tip.ott = unlist(lapply(plants_merged$tip.label, function(x) strsplit(x, "_ott")[[1]][2]))
#there are NA values

plants_merged[is.na(plants_merged$tip.ott),]
# tip.ott           tip.label search_string unique_name approximate_match is_synonym flags number_matches
# 100     <NA> mrcaott5016ott17974          <NA>        <NA>                NA         NA  <NA>             NA
# tol_node_info(5016) #> Microcystis wesenbergii
# tol_node_info(17974) #> Microcystis aeroginosa > change tip.ott to ott from resolved_plants
plants_merged$tip.ott[plants_merged$tip.label == "mrcaott5016ott17974"] = resolved_plants$ott_id[grepl("icrocystis",resolved_plants$unique_name)]

# 100    <NA>    mrcaott89ott1043          <NA>        <NA>                NA         NA  <NA>             NA
#tol_node_info(89) #Rhizobium sp.
#tol_node_info(1043) #Rhizobium sp.
plants_merged$tip.ott[plants_merged$tip.label == "mrcaott89ott1043"] = resolved_plants$ott_id[grepl("izobium",resolved_plants$unique_name)]

# 101    <NA> mrcaott1976ott20943          <NA>        <NA>                NA         NA  <NA>             NA
#tol_node_info(1976) # Dolichospermum sigmoideum
#tol_node_info(20943) # Dolichospermum flosaquae
plants_merged$tip.ott[plants_merged$tip.label == "mrcaott1976ott20943"] = resolved_plants$ott_id[grepl("olichospermum",resolved_plants$unique_name)]

# merge name info 
plants_merged = merge(plants_merged, as.data.frame(resolved_plants), by.x = "tip.ott", by.y = "ott_id", all.y = FALSE, all.x = TRUE,sort = FALSE) 

##Which species are missing in tree? - 3 tips are missing now
resolved_plants[!(resolved_plants$unique_name %in% plants_merged$unique_name),]

# search_string     unique_name approximate_match  ott_id is_synonym                    flags number_matches
# 32    brassica rapa   Brassica rapa             FALSE  833632      FALSE                                       1
# Brassica rapa subsp. pekinensis is in the tree > change label to Brassica rapa
plants_merged$unique_name[grepl("Brassica rapa",plants_merged$unique_name)] = "Brassica rapa"
# 246           vigna           Vigna             FALSE  560323      FALSE                                       1
##>> Vigna is not in the tree

# merge trait info - in order of tip labels-----

# using search_string:
# merge trait data based on sci_nem to search_string in merged:
merged = merge(plants_merged, plant_trait, by.x = "search_string", by.y = "sci_name", sort = FALSE, all.y = FALSE, all.x = TRUE,)
dim(merged) #

# check order of species:
head(merged$unique_name, 20)
head(planttree$tip.label, 20)
# ok worked.

#check lengths:
nrow(merged)
length(planttree$tip.label)
## worked

### plot tree ----


# change tip label names to unique_name
planttree$tip.label <- merged$unique_name


# genome_info into data frame and add row names
genome_info = data.frame(Available = merged$genome_species, 'High quality' = merged$good_quality)
#there are NA values in unique_name for last three tips >use tiplabels as rownames for all others:
rownames(genome_info) <- planttree$tip.label



## try out tree fan tree with additional group information
circ = ggtree(planttree, layout = "circular")
circ
# Add genome information
(p1 = gheatmap(circ, genome_info, offset = -0.5, width = .2,
               colnames_angle = 0, colnames_offset_y = c(1,1.1), 
               colnames_offset_x = c(-1,1.2), colnames_position = "top",
               font.size = 3)) #, hjust = c(1,0)
(p1 = p1 + theme(legend.position = 'none'))
(p1 = p1 + scale_fill_manual(values = c("grey90","darkgreen")))
(p1 + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 4, fontface = 3, colour = "black"))

# Add branch colors for habitat
grp = list(freshwater = planttree$tip.label[merged$habitat == "freshwater"],
           saltwater = planttree$tip.label[merged$habitat == "saltwater"],
           terrestrial = planttree$tip.label[merged$habitat == "terrestrial"])
(p3 = groupOTU(p1, grp, "Habitat") + aes(color = Habitat)) 
(p3 = p3 + scale_color_manual(values = c("lightskyblue2", "steelblue3", "burlywood4")))
(p4 = p3 + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 4, fontface = 3, colour = "black"))

png("Plots/planttree.png", res = 1000, units = "cm", height = 20, width = 22)
p4
dev.off()

### END ###



