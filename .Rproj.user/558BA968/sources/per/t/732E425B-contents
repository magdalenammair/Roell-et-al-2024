### Create phylogenetic tree from only invertebrate ott IDs and information on genomic resources and include tip labels ###

# load packages
library(rotl)
library(ggtree)
library(ggplot2)

# import data---------

# phylogenetic tree species list with correct ott numbers:
resolved_species = readRDS("Data/resolved_species_for_trees.rds")
# trait list with adjusted species names (synonyms removed)
trait_data = readRDS("Data/spp_trait_data.rds")

# create tree - algae and plants-------
invert_trait = trait_data[trait_data$tax_category == "invertebrate",]

#change sci_name to all small letters
invert_trait$sci_name = tolower(invert_trait$sci_name)

# subset ott list
resolved_invert = resolved_species[resolved_species$search_string %in% invert_trait$sci_name,]
invert_ott = resolved_invert$ott_id
invert_in_tree = is_in_tree(invert_ott)
resolved_invert[!invert_in_tree,] # all species are in the tree
dim(resolved_invert)

# create tree
inverttree = tol_induced_subtree(ott_ids = invert_ott)
length(inverttree$tip.label) # 88 tips left, 1 tip is missing!

# merge name info from resolved species list in order of tip labels ----
invert_merged = data.frame(
  tip.label = inverttree$tip.label
)
(invert_merged$tip.ott = unlist(lapply(invert_merged$tip.label, function(x) strsplit(x, "_ott")[[1]][2])))
#ok. no NA values

# merge name info 
invert_merged = merge(invert_merged, as.data.frame(resolved_invert), by.x = "tip.ott", by.y = "ott_id", all.y = FALSE, all.x = TRUE,sort = FALSE) 

##Which species is missing in tree? - 1 tip is missing:
resolved_invert[!(resolved_invert$unique_name %in% invert_merged$unique_name),]

# search_string unique_name approximate_match ott_id is_synonym flags number_matches
# 48 chironomus riparius  Chironomus             FALSE 269685      FALSE                    2
## >> try whether tree uses mrca


# merge trait info - in order of tip labels-----

# now, merge information from dat using search_string:
# merge trait data based on sci_nem to search_string in merged:
merged = merge(invert_merged, invert_trait, by.x = "search_string", by.y = "sci_name", sort = FALSE, all.y = FALSE, all.x = TRUE,)
dim(merged) #

# check order of species:
head(merged$unique_name, 20)
head(inverttree$tip.label, 20)
# great, worked.

#check lengths:
nrow(merged)
length(inverttree$tip.label)
## worked

### plot tree ----


# change tip label names to unique_name
inverttree$tip.label <- merged$unique_name


# genome_info into data frame and add row names
genome_info = data.frame(Available = merged$genome_species,
                         'High quality' = merged$good_quality)
#there are NA values in unique_name for last three tips >use tiplabels as rownames for all others:
rownames(genome_info) <- inverttree$tip.label
#colnames(genome_info) = c("Genus", "Species")

## try out tree fan tree with additional group information
circ = ggtree(inverttree, layout = "circular")
circ
# Add genome information
(p1 = gheatmap(circ, genome_info, offset = -0.5, width = .2,
               colnames_angle = 0, colnames_offset_y = c(1,1.1), 
               colnames_offset_x = c(-1, 1.2), colnames_position = "top",
               font.size = 3)) #, hjust = c(1,0)
(p1 = p1 + theme(legend.position = 'none'))
(p1 = p1 + scale_fill_manual(values = c("grey90","darkslateblue")))

(p2 = p1 + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 4, fontface = 3, colour = "black"))

# Add branch colors for habitat
grp = list(freshwater = inverttree$tip.label[merged$habitat == "freshwater"],
           saltwater = inverttree$tip.label[merged$habitat == "saltwater"],
           terrestrial = inverttree$tip.label[merged$habitat == "terrestrial"])
(p3 = groupOTU(p1, grp, "Habitat") + aes(color = Habitat)) 
(p3 = p3 + scale_color_manual(values = c("lightskyblue2", "steelblue3", "burlywood4")))
(p4 = p3 + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 4, fontface = 3, colour = "black"))

png("Plots/inverttree.png", res = 1000, units = "cm", height = 20, width = 22)
p4
dev.off()

### END ###