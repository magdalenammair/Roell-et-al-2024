## Data preparation ##

## Aim: Get ott IDs for all species names extracted from the guidance documents. Identify synonyms and check accurate identifiaction ##

#Load packages ------
library(rotl)
library(dplyr)

#Import data-----
dat <- read.csv("Data/2022_guidline_surrogate_sp_genomes.csv", header = T, stringsAsFactors = TRUE)
str(dat)
length(unique(dat$sci_name))
#257 unique species names


# 0. Get initial ott IDs from ToL ------
#remove microbial community and protozoan community:
dat <- dat[!(dat$sci_name %in% c("Microbial community", "Protozoan community")),]

# Create species list:
species = as.character(dat$sci_name)

# get ott IDs for species
resolved.species <-  tnrs_match_names(names = species)
resolved.species # gives ott_id (Open Tree Taxonomy ID) and whether the species name is a synonym

# 1. check for species not found: NA in unique_names ------
resolved.species[is.na(resolved.species$unique_name),]
# none

# 2. check for duplicate otts: ------
nrow(resolved.species) == length(unique(resolved.species$ott_id))
#>> FALSE -> there are duplicates

# get duplicates and adjust according to open taxonomic databases
dups <-resolved.species[duplicated(resolved.species$ott_id) |duplicated(resolved.species$ott_id, fromLast=TRUE),]
dups 

# check species status of duplicate otts in external taxonomy databases: WoRMS, CoL

# search_string          unique_name approximate_match ott_id is_synonym flags number_matches
# 90   farfantepenaeus aztecus      Penaeus aztecus             FALSE 736322       TRUE                    1 >> synonym to Penaeus aztecus (WoRMS, https://www.marinespecies.org/index.php), remove
# 91  farfantepenaeus duorarum     Penaeus duorarum             FALSE 163894       TRUE                    1 >> synonym to Penaeus duorarum (WoRMS), remove
# 123     lemna aequinoctialis Lemna aequinoctialis             FALSE    290      FALSE                    1 >> correct name, keep
# 127       lemna paucicostata Lemna aequinoctialis             FALSE    290       TRUE                    1 >> synonym of L. aequinoctialis, remove
# 139    litopenaeus setiferus    Penaeus setiferus             FALSE 163881       TRUE                    1 >> synonym to Penaeus setiferus (WoRMS), remove
# 190          penaeus aztecus      Penaeus aztecus             FALSE 736322      FALSE                    1 >> accepted name on World Register of Marine Species (WoRMS), keep 
# 191         penaeus duorarum     Penaeus duorarum             FALSE 163894      FALSE                    1 >> correct name (WoRMS), keep
# 192        penaeus setiferus    Penaeus setiferus             FALSE 163881      FALSE                    1 >> correct name (WoRMS), keep


# remove synonymous species from dataset
notsel = dat$sci_name %in% c("Farfantepenaeus aztecus", "Farfantepenaeus duorarum", "Litopenaeus setiferus", 
                             "Lemna paucicostata")
dat <- dat[!notsel,]
remove(notsel)
# drop unused factor levels:
dat <- droplevels(dat)

str(dat) # 251 species names, but 253 rows

# 3. check for duplicates in sci_name: ----
dat[duplicated(dat$sci_name) | duplicated(dat$sci_name, fromLast = TRUE),]
# remove second mention of species:
notsel = duplicated(dat$sci_name, fromLast = TRUE)
dat =dat[!notsel,]
str(dat) # prefect, 251 species, 251 rows

# get ott numbers for the species and check, if duplicates are now absent:
species = as.character(dat$sci_name)
resolved_species <-  tnrs_match_names(names = species)
resolved_species 

# check for duplicate otts:
nrow(resolved_species) == length(unique(resolved_species$ott_id))
#>> TRUE - no more duplicates

# save dat for merging traits to tree later:
saveRDS(dat, "Data/spp_trait_data.rds")

# save resolved_species list for manual checks: TP checks all OTT numbers again manually
saveRDS(resolved_species, "Data/resolved_species_list_for_manual_check.rds")

# 4. Adjust ott according to manual checks by TP----
# 
# Suggestions for changes based on manual ott checks are in: resolved_species_list_manual_correction_TP.csv
# unique_name is later used as tip labels, not for finding positions in the tree. they can be changed
# search_string is later used to merge the grouping information. This should not be changed
corrections = read.csv("Data/resolved_species_list_manual_correction_TP.csv", na.strings = "")
resolved_species = readRDS("Data/resolved_species_list_for_manual_check.rds")

species_correction_list = as.character(corrections$search_string[!is.na(corrections$new.ott)])
species_correction_list 

lines.to.correct = as.numeric(row.names(corrections[!is.na(corrections$new.ott),]))

for(i in species_correction_list){
  resolved_species$ott_id[resolved_species$search_string == i] = corrections$new.ott[corrections$search_string == i]
}

dups <-resolved_species[duplicated(resolved_species$ott_id) |duplicated(resolved_species$ott_id, fromLast=TRUE),]
dups
# There are two Chironomus species set to genus ott now. Remove Chironomus yoshimatsui later:
# do this later, because this seems to distort sorting in the inspect function ## 

# change Spirodela polyrhiza name to lemna major in unique_name,because the species is commonly known as L. major:
resolved_species$unique_name[grepl("Spirodela",resolved_species$unique_name)] = "Lemna major"

# Check unique_names again
resolved_species$unique_name

# 5. check species names for which is_synonym is TRUE: check on CoL, WoRMs -----
# adjust ott number where needed, adjust unique_name where appropriate
resolved_species[resolved_species$is_synonym == TRUE,]

# search_string               unique_name approximate_match  ott_id is_synonym                    flags number_matches
# 2              acartia clausi           Acartia clausii             FALSE 1057154       TRUE           sibling_higher              1 > correct species name is Acartia clausi (CoL); change to A. clausi
resolved_species$unique_name[resolved_species$unique_name == "Acartia clausii"] = "Acartia clausi"
# 5             agrostis tenuis           Agrostis canina             FALSE  550221       TRUE                                       3 > Agrostis tenuis is synonym of either A. canina, A. capillaris or A. idahoensis (CoL) > keep A. canina
# 14            amnicola limosa          Lyogyrus limosus             FALSE 1044239       TRUE                                       1 > both Amnicola limosa and Lyogorus limosus are synonym of A. limosus (CoL); change name to A. limosus
resolved_species$unique_name[resolved_species$unique_name == "Lyogyrus limosus"] = "Amnicola limosus"
# 16        anabaena flos-aquae Dolichospermum flos-aquae             FALSE  368766       TRUE           sibling_higher              2 > Anabaena flos-aquae is basionym of Dolichospermum flosaquae; D. flos-aquae is synonym of D. flosaquae (gbif); change to D. flosaquae
resolved_species$unique_name[resolved_species$unique_name == "Dolichospermum flos-aquae"] = "Dolichospermum flosaquae"
# 42      catostomus commersoni    Catostomus commersonii             FALSE  955101       TRUE                                       1 > C. commersonii is correct name (CoL); keep
# 60          crassostrea gigas           Magallana gigas             FALSE  987409       TRUE                                       1 > Crassostrea gigas is synonym of Magellana gigas (CoL); keep M. gigas
# 78      echinochloa crusgalli    Echinochloa crus-galli             FALSE    2587       TRUE                                       1 > correct species name is E. crus-galli (CoL); keep
# 115       hypoaspis aculeifer     Gaeolaelaps aculeifer             FALSE  158257       TRUE                                       1 > Hypoaspis aculeifer is synonym of Gaeolaelaps aculeifer (gbif); keep G. aculeifer
# 141       lychnis flos-cuculi        Silene flos-cuculi             FALSE  569458       TRUE                                       1 > Lychnis flos-cuculi is synonym of Silene flos-cuculi (CoL); keep S. flos-cuculi
# 156      navicula pelliculosa   Fistulifera pelliculosa             FALSE  525510       TRUE           sibling_higher              1 > both are synonyms of Synedra minutissima var. pelliculosa (WoRMS), but seems to be unresolved; keep ott and change name to species name Synedra minutissima
resolved_species$unique_name[resolved_species$unique_name == "Fistulifera pelliculosa"] = "Synedra minutissima"
# 161 oncorhynchus tschawytscha  Oncorhynchus tshawytscha             FALSE  730762       TRUE                                       1 > O. tschawytscha is synonym to O. tshawytscha (CoL); keep O. tshawytscha
# 164    orthonychiurus folsomi        Onychiurus folsomi             FALSE 5011511       TRUE                                       1 > not the same species; change ott to genus Orthonychiurus spp. and change name to Orthonychiurus folsomi (accepted species, CoL)
orthonychiurus = tnrs_match_names("orthonychiurus")
is_in_tree(orthonychiurus$ott_id) #ott is in tree
resolved_species$ott_id[resolved_species$unique_name == "Onychiurus folsomi"] = orthonychiurus$ott_id
resolved_species$unique_name[resolved_species$unique_name == "Onychiurus folsomi"] = "Orthonychiurus folsomi"
# 171  palaemonetes intermedius      Palaemon mundusnovus             FALSE  382634       TRUE           sibling_higher              1 > Palaemonetes intermedius is synonym of Palaemon mundusnovus (CoL); keep Palaemon mundusnovus
# 172        palaemonetes pugio            Palaemon pugio             FALSE   41512       TRUE           sibling_higher              1 > Palaemonetes pugio is synonym of Palaemon pugio (CoL): keep Palaemon pugio
# 173     palaemonetes vulgaris         Palaemon vulgaris             FALSE  648346       TRUE           sibling_higher              1 > Palaemonetes vulgaris is synonym of Palaemon vulgaris (CoL): keep Palaemon vulgaris
# 193            physella acuta               Physa acuta             FALSE   55904       TRUE           sibling_higher              1 > Physa acuta is synonym of Physella acuta (CoL): change to Physella acuta
resolved_species$unique_name[resolved_species$unique_name == "Physa acuta"] = "Physella acuta"
# 199          poephila guttata       Taeniopygia guttata             FALSE  708327       TRUE                                       1 > Poephila guttata is synonym of T. guttata (gbif): keep T. guttata
# 200     polygonum convolvulus      Fallopia convolvulus             FALSE  362478       TRUE                                       1 > Polygonum convolvulus is synonym to Fallopia convolvulus (CoL): keep Fallopia convolvulus
# 201   polygonum lapathifolium   Persicaria lapathifolia             FALSE 1006197       TRUE                                       1 > Polygonum lapathifolium is synonym to Persicaria lapathifolia subsp. lapathifolia (CoL); keep Persicaria lapathifolia
# 202  polygonum pennsylvanicum   Persicaria pensylvanica             FALSE  701572       TRUE                                       1 > Polygonum pennsylvanicum is synonym to Persicaria pensylvanica (CoL): keep Persicaria pensylvanica
# 225         sesbania exaltata         Sesbania herbacea             FALSE   72192       TRUE           sibling_higher              1 > Sesbania exaltata is synonym to S. herbacea (CoL): keep S. herbacea
# 234       stachys officinalis      Betonica officinalis             FALSE  103745       TRUE incertae_sedis_inherited              1 > Stachys officinalis is synonym to  Betonica officinalis (CoL): keep  Betonica officinalis
# 238    trichogramma cacoeciae    Trichogramma cacaeciae             FALSE 1003031       TRUE                                       1 > names match
# 247    xanthium pensylvanicum          Xanthium pungens             FALSE 7051609       TRUE                                       2 > Xanthium pensylvanicum is synonym to X. orientale (CoL, gbif): change ott and name to X. orientale
xanthium = tnrs_match_names("xanthium orientale")
is_in_tree(xanthium$ott_id) #is in tree
resolved_species$ott_id[resolved_species$unique_name == "Xanthium pungens"] = xanthium$ott_id
resolved_species$unique_name[resolved_species$unique_name == "Xanthium pungens"] = "Xanthium orientale"

# 6. Check whether species/genera are in correct domain -----
# check again species with information in brackets in unique_name:
resolved_species[grepl(" (",resolved_species$unique_name, fixed = TRUE),]

# 
# search_string                                            unique_name approximate_match     score  ott_id is_synonym          flags number_matches
# 20       apis mellifera                    Apis mellifera (in domain Bacteria)             FALSE 1.0000000 5900566      FALSE                             2 > should be in domain Eukaryota
inspect(resolved_species, taxon_name = "apis mellifera")
resolved_species$ott_id[resolved_species$search_string == "apis mellifera"] = 461645
# 28    bombus terrestris                 Bombus terrestris (in domain Bacteria)             FALSE 1.0000000 5901391      FALSE                             2 > should be in domain Eukaryota
inspect(resolved_species, taxon_name = "bombus terrestris")
resolved_species$ott_id[resolved_species$search_string == "bombus terrestris"] = 161197
# 48  chironomus riparius               Chironomus riparius (in domain Bacteria)             FALSE 1.0000000  269685      FALSE                             2 > should be in domain Eukaryota
inspect(resolved_species, taxon_name = "chironomus riparius")
is_in_tree(209337) #not in tree. Leave genus ott 
# 62      cucumis sativus                   Cucumis sativus (in domain Bacteria)             FALSE 1.0000000 5893675      FALSE                             2 > should be in domain Eukaryota
inspect(resolved_species, taxon_name = "cucumis sativus")
resolved_species$ott_id[resolved_species$search_string == "cucumis sativus"] = 1006246
# 75  digitalis purpurea  Digitalis purpurea (species in kingdom Archaeplastida)              TRUE 0.9444444 6072863      FALSE                             4 > ok
# 88          ephemerella          Ephemerella (genus in kingdom Archaeplastida)             FALSE 1.0000000  662492      FALSE         barren              2 > should be in Ephemeroptera
inspect(resolved_species, taxon_name = "ephemerella") #should be the second entry
resolved_species$ott_id[resolved_species$search_string == "ephemerella"] = 944045
# 96         gadus morhua                      Gadus morhua (in domain Bacteria)             FALSE 1.0000000 5905500      FALSE                             2 > should be in domain Eukaryota
inspect(resolved_species, taxon_name = "gadus morhua") #should be the second entry
resolved_species$ott_id[resolved_species$search_string == "gadus morhua"] = 114170
# 160 oncorhynchus mykiss               Oncorhynchus mykiss (in domain Bacteria)             FALSE 1.0000000 5256670      FALSE                             2 > should be in domain Eukaryota
inspect(resolved_species, taxon_name = "oncorhynchus mykiss") #should be the second entry
resolved_species$ott_id[resolved_species$search_string == "oncorhynchus mykiss"] = 165368
# 168        pachygrapsus               Pachygrapsus (subgenus in genus Grapsus)             FALSE 1.0000000 7117774      FALSE         barren              2 >ok
# 215           rhizobium               Rhizobium (genus in family Rhizobiaceae)             FALSE 1.0000000  263867      FALSE                             2 > ok
# 226       sida spinosa        Sida spinosa (species in kingdom Archaeplastida)              TRUE 0.9166667  991909      FALSE sibling_higher              2 > ok


# Remove information in brackets
resolved_species$unique_name = unlist(lapply(1:nrow(resolved_species), function(k) strsplit(resolved_species$unique_name, " (", fixed = TRUE)[[k]][1]))

# 7. check whether all ott numbers are in synthetic tree:----------
sel = is_in_tree(resolved_species$ott_id)
resolved_species[!sel,]

# Find solutions for species not in tree:
# > resolved_species[!sel,]
# search_string              unique_name approximate_match  ott_id is_synonym                    flags number_matches
# 16     anabaena flos-aquae Dolichospermum flosaquae             FALSE  368766       TRUE           sibling_higher              2 > neglect in tree- tree finds common ancestor
#check genus name:
(anabaena = tnrs_match_names("anabaena"))
is_in_tree(anabaena$ott_id) #also not in tree
# check family:
(nostocaceae = tnrs_match_names("Nostocaceae"))
is_in_tree(nostocaceae$ott_id) #also not in tree 
#check other species in same genus:
(anabaenacyl = tnrs_match_names("Anabaena cylindrica"))
is_in_tree(anabaenacyl$ott_id) #also not in tree

# 36                cambarus                 Cambarus             FALSE  876591      FALSE                                       1 > change to family ott
#check family
(cambaridae = tnrs_match_names("cambaridae"))
is_in_tree(cambaridae$ott_id) #is in tree > 
# change to family cambaridae
resolved_species$ott_id[resolved_species$search_string == "cambarus"] = cambaridae$ott_id

# 66   cyprinodon variegatus    Cyprinodon variegatus             FALSE   82335      FALSE                                       1 > change to genus ott
(cyprinodon = tnrs_match_names("cyprinodon"))
is_in_tree(cyprinodon$ott_id) # is in tree
resolved_species$ott_id[resolved_species$search_string == "cyprinodon variegatus"] = cyprinodon$ott_id

# 150 microcystis aeruginosa   Microcystis aeruginosa             FALSE  269069      FALSE                                       1 > neglect in tree- tree finds common ancestor
# check genus
(microcystis = tnrs_match_names("microcystis"))
is_in_tree(microcystis$ott_id) # not in tree
#check family
(chroococcaceae = tnrs_match_names("chroococcaceae"))
is_in_tree(chroococcaceae$ott_id) #not in tree

# 152       musca autumnalis         Musca autumnalis             FALSE 4364297      FALSE           incertae_sedis              1 > change to genus ott
# check genus
(musca = tnrs_match_names("musca"))
is_in_tree(musca$ott_id) # is in tree
resolved_species$ott_id[resolved_species$search_string == "musca autumnalis"] = musca$ott_id

# 168           pachygrapsus             Pachygrapsus             FALSE 7117774      FALSE                   barren              2 > change to family ott
# check family
(grapsidae = tnrs_match_names("grapsidae"))
is_in_tree(grapsidae$ott_id) # is in tree
resolved_species$ott_id[resolved_species$search_string == "pachygrapsus"] = grapsidae$ott_id

# 205            procambarus              Procambarus             FALSE   79859      FALSE                                       1 > change to family ott
# check family
(cambaridae = tnrs_match_names("cambaridae"))
is_in_tree(cambaridae$ott_id) # is in tree
resolved_species$ott_id[resolved_species$search_string == "procambarus"] = cambaridae$ott_id

# 209                   rana                     Rana             FALSE  364550      FALSE                                       1 > change to family ott
# check family
(ranidae = tnrs_match_names("ranidae"))
is_in_tree(ranidae$ott_id) # is in tree
resolved_species$ott_id[resolved_species$search_string == "rana"] = ranidae$ott_id

# 215              rhizobium                Rhizobium             FALSE  263867      FALSE                                       2 > neglect in tree - tree finds common ancestor
# check family
(rhizobiaceae = tnrs_match_names("rhizobiaceae"))
is_in_tree(rhizobiaceae$ott_id) # not in tree

# 246                  vigna                    Vigna             FALSE  904024      FALSE                                       1 > change to family ott
# check family
(fabaceae = tnrs_match_names("fabaceae"))
is_in_tree(fabaceae$ott_id) # is in tree
resolved_species$ott_id[resolved_species$search_string == "vigna"] = fabaceae$ott_id

# check again for duplicate otts:
dups <-resolved_species[duplicated(resolved_species$ott_id) |duplicated(resolved_species$ott_id, fromLast=TRUE),]
dups

# 7. Remove duplicates
# remove second Chironomus spp.:
resolved_species = resolved_species[!grepl("yoshimatsui",resolved_species$unique_name),] 
# remove second Procamburus spp. within family Cambaridae:
resolved_species = resolved_species[!grepl("Proc",resolved_species$unique_name),] 

# 8. Final check:
resolved_species$unique
intree = is_in_tree(resolved_species$ott_id)
sum(!intree) # only 3 items left. Most recent common ancestor nodes are used for them 
resolved_species[!intree,] #Dolichospermum, Microcystis, Rhizobium
# final duplicates check
(dups <-resolved_species[duplicated(resolved_species$ott_id) |duplicated(resolved_species$ott_id, fromLast=TRUE),])

# 8. Save resolved species data -------
saveRDS(resolved_species, "Data/resolved_species_for_trees.rds")

### END ##### 







