### Create phylogenetic tree from only vertebrate ott IDs and information on genomic resources and include tip labels ###

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
verteb_trait = trait_data[trait_data$tax_category == "vertebrates",]

#change sci_name to all small letters
verteb_trait$sci_name = tolower(verteb_trait$sci_name)

# subset ott list
resolved_verteb = resolved_species[resolved_species$search_string %in% verteb_trait$sci_name,]
verteb_ott = resolved_verteb$ott_id
verteb_in_tree = is_in_tree(verteb_ott)
resolved_verteb[!verteb_in_tree,] # all species are in the tree
dim(resolved_verteb)

# create tree
vertebtree = tol_induced_subtree(ott_ids = verteb_ott)
length(vertebtree$tip.label) # all 56 tips present, no tip lost

# merge name info from resolved species list in order of tip labels ----
verteb_merged = data.frame(
  tip.label = vertebtree$tip.label
)
verteb_merged$tip.ott = unlist(lapply(verteb_merged$tip.label, function(x) strsplit(x, "_ott")[[1]][2]))
#ok. no NA values

# merge name info 
verteb_merged = merge(verteb_merged, as.data.frame(resolved_verteb), by.x = "tip.ott", by.y = "ott_id", all.y = FALSE, all.x = TRUE,sort = FALSE) 
# invert_merged$tip.ott
# tail(plants_merged) 
# str(plants_merged)

# merge trait info - in order of tip labels-----

# now, merge information from dat using search_string:
# merge trait data based on sci_nem to search_string in merged:
merged = merge(verteb_merged, verteb_trait, by.x = "search_string", by.y = "sci_name", sort = FALSE, all.y = FALSE, all.x = TRUE,)
dim(merged) #

# check order of species:
head(merged$unique_name, 20)
head(vertebtree$tip.label, 20)
# great, worked.

#check lengths:
nrow(merged)
length(vertebtree$tip.label)
## worked

### plot tree ----


# change tip label names to unique_name
vertebtree$tip.label <- merged$unique_name


# genome_info into data frame and add row names
genome_info = data.frame(Available = merged$genome_species,
                         'High quality' = merged$good_quality)
#there are NA values in unique_name for last three tips >use tiplabels as rownames for all others:
rownames(genome_info) <- vertebtree$tip.label
# colnames(genome_info) = c("Genus", "Species")



## try out tree fan tree with additional group information
circ = ggtree(vertebtree, layout = "circular")
circ + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 3, fontface = 3, colour = "black")

# flip to show birds next to frogs and toads:
# check node numbers
circ + geom_text(aes(label=node), hjust=-.3) + geom_tiplab()
# Flip
circ = flip(circ, 59,71)
circ

# Add genome information
(p1 = gheatmap(circ, genome_info, offset = -0.4, width = .2,
               colnames_angle = 0, colnames_offset_y = c(1,1.1), 
               colnames_offset_x = c(-1, 1.2), colnames_position = "top",
               font.size = 3)) #, hjust = c(1,0)
(p1 = p1 + theme(legend.position = 'none'))
(p1 = p1 + scale_fill_manual(values = c("grey90","darkorange3")))

(p1 = p1 + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 3, fontface = 3, colour = "black"))

# Add branch colors for habitat
grp = list(freshwater = vertebtree$tip.label[merged$habitat == "freshwater"],
           saltwater = vertebtree$tip.label[merged$habitat == "saltwater"],
           terrestrial = vertebtree$tip.label[merged$habitat == "terrestrial"])
(p2 = groupOTU(p1, grp, "Habitat") + aes(color = Habitat)) 
p2 = p2 + scale_color_manual(values = c("lightskyblue2", "steelblue3", "burlywood4"))
(p2 = p2 + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 3, fontface = 3, colour = "black"))

png("Plots/vertebtree.png", res = 1000, units = "cm", height = 20, width = 22)
p2
dev.off()


### END ###