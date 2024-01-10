### Create phylogenetic tree from ott IDs and information on genomic resources ###

# Aim: Create synthetic tree based on ott numbers in resolved_species. Merge species names and genome information to tips of the synthetic tree
# Steps: 
# 1. create synthetic tree and check tip labels
# 2. merge species names from resolved_species to labelled tree based on ott IDs
# 3. merge genome information from trait_data to labelled tree based on species names

# load packages
library(rotl)
library(ggtree)
library(ggplot2)

# import data---------

# phylogenetic tree species list with correct ott numbers:
resolved_species = readRDS("Data/resolved_species_for_trees.rds")
# trait list with adjusted species names (synonyms removed)
trait_data = readRDS("Data/spp_trait_data.rds")

#change sci_name to all small letters
trait_data$sci_name = tolower(trait_data$sci_name)

# Exclude samples with genus names only (no species names)
# number of genus-only entries
length(resolved_species$unique_name[!grepl(" ", resolved_species$unique_name)]) # 13
# remove genus-only entries from dataset
resolved_species = resolved_species[grepl(" ", resolved_species$unique_name),]
nrow(resolved_species)

# 1. Create tree and check tip labels ----------
# check otts in tree again
ott = resolved_species$ott_id
length(ott) #236
in_tree = is_in_tree(ott) 
resolved_species[!in_tree,] # these 2 species are not found in tree. Their most recent common ancestor (mrca) will be used for plotting.

# create tree
fulltree = tol_induced_subtree(ott_ids = ott)
length(fulltree$tip.label) # 234 tips left, 2 tips are missing

# The tip labels are now created by rotl----
# To change this add name info from resolved species list in order of tip labels 
full_merged = data.frame(
  tip.label = fulltree$tip.label
)
(full_merged$tip.ott = unlist(lapply(full_merged$tip.label, function(x) strsplit(x, "_ott")[[1]][2])))
#there are 2 NA values

# check NAs - get info for the two merged ott numbers
full_merged[is.na(full_merged$tip.ott),]
# tip.label tip.ott
# 230 mrcaott1976ott20943    <NA>
tol_node_info(1976) # Dolichospermum sigmoideum
tol_node_info(20943) # Dolichospermum flosaquae > change tip.ott to ott from resolved_species. This does not distort the structure of the tree.
full_merged$tip.ott[full_merged$tip.label == "mrcaott1976ott20943"] = resolved_species$ott_id[grepl("olichospermum",resolved_species$unique_name)]

# 231 mrcaott5016ott17974    <NA>
tol_node_info(5016) #> Microcystis wesenbergii
tol_node_info(17974) #> Microcystis aeroginosa > change tip.ott to ott from resolved_species. This does not distort the structure of the tree.
full_merged$tip.ott[full_merged$tip.label == "mrcaott5016ott17974"] = resolved_species$ott_id[grepl("icrocystis",resolved_species$unique_name)]


# Add species names from resolved_species to dataset according to otts
full_merged = merge(full_merged, as.data.frame(resolved_species), by.x = "tip.ott", by.y = "ott_id", all.y = FALSE, all.x = TRUE,sort = FALSE) 
tail(full_merged) 

# Check dataset
str(full_merged)
# 2 species are now missing in the tree

## Check species names - Which 2 tips are missing?
resolved_species[!(resolved_species$unique_name %in% full_merged$unique_name),]
# brassica rapa
# chironomus riparius

full_merged[grepl("Brassica", full_merged$unique_name),]
resolved_species[grepl("Brassica", resolved_species$unique_name),]

# Brassica rapa is not in the tree, but:
# Brassica rapa subsp. pekinensis is in the tree 
# We want to display B. rapa instead of B rapa subsp. pekinensis
# > change unique_name to Brassica rapa. This is used as tip labels in the plot
full_merged$unique_name[grepl("Brassica rapa",full_merged$unique_name)] = "Brassica rapa"
# > change search_string to brassica rapa. This is used for merging genome info
full_merged$search_string[grepl("Brassica rapa",full_merged$unique_name)] = "brassica rapa"

resolved_species[grepl("Chironomus", resolved_species$unique_name),]
#Chironomus riparius is not in the tree, but Ch. dilutus is. Set tip label to "Chironomus spp.". There are 3 Chironomus species of which only Ch. dilutus is in the treer right now 
# > change unique_name to Chironomus spp. This is used as tip labels in the plot
full_merged$unique_name[grepl("Chironomus dilutus",full_merged$unique_name)] = "Chironomus spp."
# > change search_string to Chironomus riparius. This is used for merging genome info
full_merged$search_string[grepl("Chironomus",full_merged$unique_name)] = "chironomus riparius"


# merge trait info - in order of tip labels-----

# now, merge information from dat using search_string:
# merge trait data based on sci_nem to search_string in merged:
merged = merge(full_merged, trait_data, by.x = "search_string", by.y = "sci_name", sort = FALSE, all.y = FALSE, all.x = TRUE,)
dim(merged) #

# check order of species:
head(merged$unique_name, 20)
head(fulltree$tip.label, 20)
# ok, worked.

#check lengths:
nrow(merged)
length(fulltree$tip.label)
## worked


# change tip label names to unique_name
fulltree$tip.label <- merged$unique_name

### plot tree ----

# genome_info into separate data frame and add row names
genome_info = data.frame(merged$genome_species, merged$good_quality, merged$tax_category)
#there are NA values in unique_name for last three tips >use tiplabels as rownames for all others:
rownames(genome_info) <- fulltree$tip.label
colnames(genome_info) = c("Available", "Quality", "tax_category")

genome_info$Quality = as.character(genome_info$Quality)
genome_info$Quality[genome_info$Quality == "yes" & genome_info$tax_category == "plant"] = "plant2"
genome_info$Quality[genome_info$Quality == "yes" & genome_info$tax_category == "algae"] = "plant2"
genome_info$Quality[genome_info$Quality == "yes" & genome_info$tax_category == "procaryote"] = "plant2"
genome_info$Quality[genome_info$Quality == "yes" & genome_info$tax_category == "invertebrate"] = "invertebrate2"
genome_info$Quality[genome_info$Quality == "yes" & genome_info$tax_category == "vertebrates"] = "vertebrates2"
genome_info$Quality = as.factor(genome_info$Quality)

genome_info$Available = as.character(genome_info$Available)
genome_info$Available[genome_info$Available == "yes" & genome_info$tax_category == "plant"] = "plant"
genome_info$Available[genome_info$Available == "yes" & genome_info$tax_category == "algae"] = "plant"
genome_info$Available[genome_info$Available == "yes" & genome_info$tax_category == "procaryote"] = "plant"
genome_info$Available[genome_info$Available == "yes" & genome_info$tax_category == "invertebrate"] = "invertebrate"
genome_info$Available[genome_info$Available == "yes" & genome_info$tax_category == "vertebrates"] = "vertebrates"
genome_info$Available = as.factor(genome_info$Available)
  
# remove column with taxon groups. Necessary for plotting
genome_info = genome_info[,-3]

#check sum available genomes and sum high quality again:
# table(genome_info$Available)
# table(genome_info$Quality) #55 with high quality; 2 are missing: Chironomus riparius is not in tree, B. rapa subsp. 

# save genome info table for summary statistics in 03_data_table.R
saveRDS(genome_info, "Data/genome_info.RDS")

## plot fan tree with additional group information
circ = ggtree(fulltree, layout = "circular")
circ

# Add genome information
(p1 = gheatmap(circ, genome_info, offset = -0.5, width = .2,
               colnames_angle = 0, colnames_offset_y = 1, 
               colnames_offset_x = c(-0.8,1), colnames_position = "top",
               font.size = 2.5))#, hjust = c(1,0)))
(p1 = p1 + theme(legend.position = 'none'))
(p1 = p1 + scale_fill_manual(values = c("darkslateblue",adjustcolor("darkslateblue", alpha.f = 0.6), 
                                        "grey85","darkgreen", adjustcolor("darkgreen", alpha.f = 0.6),
                                        "darkorange3",adjustcolor("darkorange3", alpha.f = 0.6))))

# check labels
p1 + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 4, fontface = 3, colour = "black")

# Flip node to bring all Brassica species next to each other:
# check node numbers:
p1 + geom_text(aes(label=node), hjust=-.3, size = 1.5)
# Flip
p1 = flip(p1,20,259)
# check labels again
p1 + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 4, fontface = 3, colour = "black")

# Add branch colors for habitat
grp = list(freshwater = fulltree$tip.label[merged$habitat == "freshwater"],
           saltwater = fulltree$tip.label[merged$habitat == "saltwater"],
           terrestrial = fulltree$tip.label[merged$habitat == "terrestrial"])
(p3 = groupOTU(p1, grp, "Habitat") + aes(color = Habitat)) 
(p3 = p3 + scale_color_manual(values = c("lightskyblue2", "steelblue3", "burlywood4")))
(p4 = p3 + geom_tiplab(size = 2, align=FALSE, linesize=.5, offset = 4, fontface = 3, colour = "black"))

png("Plots/fulltree.png",res = 500, units = "cm", height = 15, width = 15 )
p3
dev.off()

###END####
