### Get number of species with available and good genomic resources ###

# import data
trait_data = readRDS("Data/spp_trait_data.rds")
str(trait_data)

# remove genus level entries
# number of genus-only entries
length(trait_data$sci_name[!grepl(" ", trait_data$sci_name)]) # 14
# remove genus-only entries from dataset
trait_data = trait_data[grepl(" ", trait_data$sci_name),]
nrow(trait_data) #237

#genome availability table:
t = as.data.frame.matrix(table(trait_data$genome_species, trait_data$tax_category))
t$plantalgae = t$plant + t$algae
t = t[,c("plantalgae", "invertebrate", "vertebrates")]
t

sums = sapply(t, sum)

t.prop = t
t.prop$plantalgae.prop = t.prop$plantalgae/sums[1]
t.prop$invertebrate.prop = t.prop$invertebrate/sums[2]
t.prop$vertebrates.prop = t.prop$vertebrates/sums[3]
t.prop
###(t.prop = t(t.prop))
#(t.prop[,4:6] = round(t.prop[,4:6]*100, digits = 1))

# genome quality table:
t2 = as.data.frame.matrix(table(trait_data$good_quality, trait_data$tax_category))
t2$plantalgae = t2$plant + t2$algae
t2 = t2[,c("plantalgae", "invertebrate", "vertebrates")]
t2
sums2 = sapply(t2, sum)

t.prop2 = t2
t.prop2$plantalgae.prop = t.prop2$plantalgae/sums[1]
t.prop2$invertebrate.prop = t.prop2$invertebrate/sums[2]
t.prop2$vertebrates.prop = t.prop2$vertebrates/sums[3]
t.prop2
#(t.prop2 = round(t.prop2*100, digits = 1))
##(t.prop2 = t(t.prop2))

# final table

tab = data.frame(
  Taxon = rep(NA, 3),
  Availability = rep(NA, 3),
  Quality = rep(NA, 3)
)

tab$Taxon = c("Plants and algae", "Invertebrates", "Vertebrates")
tab$Availability[1] = paste0(round(t.prop$plantalgae.prop[2]*100, digits = 1), 
                              "% (", 
                              t.prop$plantalgae[2], 
                              "/",
                              sum(t.prop$plantalgae), 
                             ")") 
tab$Availability[2] = paste0(round(t.prop$invertebrate.prop[2]*100, digits = 1), 
                             "% (", 
                             t.prop$invertebrate[2], 
                             "/",
                             sum(t.prop$invertebrate), 
                             ")") 
tab$Availability[3] = paste0(round(t.prop$vertebrates.prop[2]*100, digits = 1), 
                             "% (", 
                             t.prop$vertebrates[2], 
                             "/",
                             sum(t.prop$vertebrates), 
                             ")") 
tab$Quality[1] = paste0(round(t.prop2$plantalgae.prop[2]*100, digits = 1), 
                             "% (", 
                             t.prop2$plantalgae[2], 
                             "/",
                             sum(t.prop2$plantalgae), 
                             ")") 
tab$Quality[2] = paste0(round(t.prop2$invertebrate.prop[2]*100, digits = 1), 
                             "% (", 
                             t.prop2$invertebrate[2], 
                             "/",
                             sum(t.prop2$invertebrate), 
                             ")") 
tab$Quality[3] = paste0(round(t.prop2$vertebrates.prop[2]*100, digits = 1), 
                             "% (", 
                             t.prop2$vertebrates[2], 
                             "/",
                             sum(t.prop2$vertebrates), 
                             ")") 
tab


tab2 = data.frame(
  Taxon = rep(NA, 3),
  Availability = rep(NA, 3),
  Quality = rep(NA, 3)
)

tab2$Taxon = c(paste0("Plants and algae (", sum(t.prop$plantalgae), ")"), 
               paste0("Invertebrates (", sum(t.prop$invertebrate), ")"),
               paste0("Vertebrates (", sum(t.prop$vertebrates), ")"))
tab2$Availability[1] = paste0(round(t.prop$plantalgae.prop[2]*100, digits = 1), 
                             "% (", 
                             t.prop$plantalgae[2], 
                             ")") 
tab2$Availability[2] = paste0(round(t.prop$invertebrate.prop[2]*100, digits = 1), 
                             "% (", 
                             t.prop$invertebrate[2], 
                             ")") 
tab2$Availability[3] = paste0(round(t.prop$vertebrates.prop[2]*100, digits = 1), 
                             "% (", 
                             t.prop$vertebrates[2], 
                             ")") 
tab2$Quality[1] = paste0(round(t.prop2$plantalgae.prop[2]*100, digits = 1), 
                        "% (", 
                        t.prop2$plantalgae[2], 
                        ")") 
tab2$Quality[2] = paste0(round(t.prop2$invertebrate.prop[2]*100, digits = 1), 
                        "% (", 
                        t.prop2$invertebrate[2], 
                        ")") 
tab2$Quality[3] = paste0(round(t.prop2$vertebrates.prop[2]*100, digits = 1), 
                        "% (", 
                        t.prop2$vertebrates[2], 
                        ")") 
tab2

saveRDS(list(tab1 = tab, tab2 = tab2), "Data/coverage_table.RDS")

###END###

#add-on - manual pie charts

op = par(mfrow = c(3,3), mar = c(0,0,0,0), oma = c(0,0,2,0), cex = 1)

plot(c(0,0), type = "n", bty = "")
text(0, 0.5, "Plants & Algae", pos = 4)
pie(c(37,101-37), col = c(adjustcolor("darkgreen", alpha.f = 0.6),"grey85"), labels = "", clockwise = TRUE, border = "grey10")
mtext("Availability", side = 3, line = -0.5, cex = 1)
pie(c(18, 101-18), col = c(adjustcolor("darkgreen", alpha.f = 0.6),"grey85"), labels = "", clockwise = TRUE, border = "grey10")
mtext("Quality", side = 3, line = -0.5, cex = 1)


plot(c(0,0), type = "n", bty = "")
text(0,0.5, "Invertebrates", pos = 4)
pie(c(25, 82-25), col = c(adjustcolor("darkslateblue", alpha.f = 0.6),"grey85"), labels = "", clockwise = TRUE, border = "grey10")
pie(c(14, 82-14), col = c(adjustcolor("darkslateblue", alpha.f = 0.6),"grey85"), labels = "", clockwise = TRUE, border = "grey10")

plot(c(0,0), type = "n", bty = "")
text(0,0.5, "Vertebrates", pos = 4)
pie(c(38, 54-38), col = c(adjustcolor("darkorange3", alpha.f = 0.6),"grey85"), labels = "", clockwise = TRUE, border = "grey10")
pie(c(24, 54-24), col = c(adjustcolor("darkorange3", alpha.f = 0.6),"grey85"), labels = "", clockwise = TRUE, border = "grey10")

par(op)    

### END ###
