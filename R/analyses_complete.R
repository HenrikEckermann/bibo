## ----import data, echo=FALSE, message=FALSE, warning=FALSE---------------
library(microbiome)
library(ggfortify)
library(cluster)
library(MASS)
library(plot3D)
library(HITChipDB)
library(readxl)
library(tidyverse)
library(here)

# L1 abundances 
L1_nolog <- as.matrix(read.delim(here("data/hitchip/L1-rpa.tab"), row.names=1))
L1 <- log10(L1_nolog)

# L1 abundances 
L2_nolog <- as.matrix(read.delim(here("data/hitchip/L2-rpa.tab"), row.names=1))
L2 <- log10(L2_nolog)

dim(L2)
L2 

as.data.frame(L2_nolog) %>%
  rownames_to_column("species") %>% gather(key, value, -species) %>%
  mutate(value = scale(value)[,1]) %>%
  ggplot(aes(value)) +
    geom_density()
# relative abundances
L2_rel <- relative.abundance(L2_nolog)
L1_rel <- relative.abundance(L1_nolog)

L2_rel %>% as.data.frame() %>%
  rownames_to_column("genus") %>%
  gather(sample, value, -genus) %>%
  spread(genus, value) %>%
  qplot(Bifidobacterium, data = ., geom = "density")
# import metadata 
my.metadata <- data.frame(read_excel(here("data/meta_variables/my.metadata.xlsx")))
rownames(my.metadata) <- my.metadata$Var.2
sum(L2_rel["Bifidobacterium", ])
L2_rel[, 1] %>% sum()

#
my.metadata$correlation_L4 <- as.numeric(my.metadata$correlation_L4)
my.metadata$correlation_L3 <- as.numeric(my.metadata$correlation_L3)
my.metadata$correlation_L2 <- as.numeric(my.metadata$correlation_L2)
my.metadata[,36:41] <- sapply(my.metadata[,36:41], as.numeric)

my.metadata$groupcode <- as.factor(my.metadata$groupcode)
my.metadata$subject <- as.factor(my.metadata$subject)
my.metadata$childcarecenter <- as.factor(my.metadata$childcarecenter)



my.metadata <- arrange(my.metadata,subject,sampleID)
rownames(my.metadata) <- my.metadata$Var.2
my.metadata$feeding_type <- ""


all.samples <- intersect(rownames(my.metadata),colnames(L2))
my.metadata <- my.metadata[all.samples,]

myd <- L2[,all.samples]
duplicates <- my.metadata$sampleID[duplicated(my.metadata$sampleID)]
duplicates
resulting_samples <- as.character(rownames(subset(my.metadata, sampleID %in% duplicates)))
resulting_samples

metadata <- my.metadata[my.metadata$Begingroup == 1,]
myd <- L2[,rownames(metadata)]

#choose yes or no antibiotics
metadata <- metadata[metadata$antibio != 1, ]

pre_KDV <- metadata[metadata$groupcode == "KDV1",][,c("Var.2","subject")]
post_KDV <- metadata[metadata$groupcode == "KDV5",][,c("Var.2","subject")]

matching_subjects <- Reduce(intersect, list(pre_KDV$subject,post_KDV$subject))

matching_subjects
selected_samples <- as.character(subset(my.metadata, subject %in% matching_subjects)[,"Var.2"])

pre_KDV <- as.character(subset(my.metadata[resulting_samples,], groupcode == "KDV1")[,c("Var.2")])
post_KDV <- as.character(subset(my.metadata[resulting_samples,], groupcode == "KDV5")[,c("Var.2")])

pre_KDV_CC <- as.character(subset(my.metadata[resulting_samples,], childcarecenter == "CC")[,c("Var.2")])
post_KDV_noCC <- as.character(subset(my.metadata[resulting_samples,], childcarecenter == "noCC")[,c("Var.2")])


resulting_samples <- as.character(subset(my.metadata[selected_samples,], groupcode == "KDV1" | groupcode == "KDV5")[,c("Var.2")])

### FILTER TAXA BY LOG ABUNDANCE
# Check dimensionality 



plot(density(L2),main = "distribution of probesignal in all samples") +
  abline(v = 1.74, col ="red") +
  abline(v = quantile(L2,0.8), col ="green")


plot(density(L2[,resulting_samples]),main = "distribution of probesignal in selected samples") +
  abline(v = 1.74, col ="red") +
  abline(v = quantile(L2[,resulting_samples],0.8), col ="green")



plot(density(L2_rel[,resulting_samples]),main = "distribution of relative abundance in selected samples") +
  abline(v = 1.74, col ="red") +
  abline(v = quantile(L2_rel[,resulting_samples],0.8), col ="green")

plot(density(L2_rel[,resulting_samples]),main = "distribution of relative abundance in selected samples, zoomed in",xlim=c(-0.005,0.01)) +
  abline(v = 1.74, col ="red") +
  abline(v = quantile(L2_rel[,resulting_samples],0.8), col ="green")

loghitchipF <- t(L2_rel[,resulting_samples])
# filter bacteria to get those in which at least 20% of the samples have abundance higher than 2.5 log
loghitchipF.mean <- apply((loghitchipF > 0.005), 2, mean)
# here you get boolean values: 0 if log < 2.4 or 1 if log > 2.5
# then with the mean you actually get the percentage of taxa with abundance > 2.5
# and add these value for each of the taxa (this case 130)
loghitchipF1.7 <- cbind(t(loghitchipF), loghitchipF.mean)
# then you set the threshold 
loghitchipF1.7 <- loghitchipF1.7[which(loghitchipF1.7[,ncol(loghitchipF1.7)] > 0.2),]
# and get rid of that column 
abundant_taxa <- rownames(loghitchipF1.7[,-ncol(loghitchipF1.7)])

loghitchipF <- t(L2[,resulting_samples])
# filter bacteria to get those in which at least 20% of the samples have abundance higher than 2.5 log
loghitchipF.mean <- apply((loghitchipF > 1.8), 2, mean)
# here you get boolean values: 0 if log < 2.4 or 1 if log > 2.5
# then with the mean you actually get the percentage of taxa with abundance > 2.5
# and add these value for each of the taxa (this case 130)
loghitchipF1.7 <- cbind(t(loghitchipF), loghitchipF.mean)
# then you set the threshold 
loghitchipF1.7 <- loghitchipF1.7[which(loghitchipF1.7[,ncol(loghitchipF1.7)] > 0.2),]
# and get rid of that column 
L2 <- loghitchipF1.7[,-ncol(loghitchipF1.7)]
## ---- echo=FALSE, warning=FALSE,message=FALSE----------------------------
L2
library(dplyr)
# change filepath to the metadata_long.csv
file <- here("data/meta_variables/metadata_long.csv")
meta_long <- read.csv(file)

# convert to formats to work with
meta_long <- 
  dplyr::select(meta_long, -X) %>%
    dplyr::select(subject, week, bf, expressed_bf, formula, solid_food, type_solid_food, pet) %>%
    arrange(subject) %>%
    mutate(
      bf = as.numeric(bf),
      expressed_bf = as.numeric(expressed_bf),
      formula = as.numeric(formula),
      solid_food = as.numeric(solid_food),
      type_solid_food = ifelse(
        as.character(type_solid_food) == "", NA, 
        as.character(type_solid_food)
      )
    )

str(meta_long)     

feeding_mode <- meta_long

feeding_mode_df <- tbl_df(feeding_mode)
subjects <- group_by(feeding_mode_df, subject)


#trying to salvage infants that have weird things with age
#my.metadata[is.na(my.metadata$Age_plus4weeks),] <- my.metadata$Begin_age_weeks[is.na(my.metadata$Age_plus4weeks)] + 4

#determine which class they're in at post childcare
input <- as.data.frame(my.metadata[my.metadata$groupcode== "KDV5",c("subject","Age_plus4weeks")])
#remove subject 451 has no age
input <- input[complete.cases(input), ]
input$Age_plus4weeks <- round(input$Age_plus4weeks)



my.data <-  subjects[subjects$subject %in% as.vector(input$subject),]

mylist <- list()
for (z in 1:nrow(input)){
  x = as.character(input[z,c("subject")])
  d = input[z,c("Age_plus4weeks")]
  m=my.data[my.data$subject == x,]
  m=m[as.numeric(paste(1:d)),]
  mylist[[z]] = m
}

feedingmode_mysamples <- do.call(rbind, mylist)


df <- as.data.frame(summarise(feedingmode_mysamples,age_weeks= n(),bf= sum(bf,na.rm = T),formula=sum(formula,na.rm = T),solid_food=sum(solid_food,na.rm = T),type_solid_food=n_distinct(type_solid_food,na.rm = TRUE)))
df <- df[complete.cases(df), ]
rownames(df) <- df$subject

meta <- my.metadata[my.metadata$groupcode == "KDV5",]

all.data1 <- merge(x = df,y = meta,by = "subject", all.x = T )


########BEWARE I CHANGED the solid food column. removed dots and a comma and also 0
all.data1$feeding_type[all.data1$bf > 0 & all.data1$formula == 0 & all.data1$solid_food==0 & all.data1$type_solid_food == 0] <- "BF" 
all.data1$feeding_type[all.data1$bf > 0 & all.data1$formula > 0 & all.data1$solid_food==0 & all.data1$type_solid_food == 0] <- "MF"
all.data1$feeding_type[all.data1$bf == 0 & all.data1$formula > 0 & all.data1$solid_food==0 & all.data1$type_solid_food == 0] <- "IF"
all.data1$feeding_type[all.data1$solid_food>0 | !(all.data1$type_solid_food == 0)] <- "MF+SF"
all.data1$feeding_type[all.data1$bf == 0 & all.data1$formula == 0 & all.data1$solid_food==0] <- "no_data"
rownames(all.data1) <- all.data1$Var.2

post <- all.data1

#determine which class they're in at pre childcare
input <- as.data.frame(my.metadata[my.metadata$groupcode== "KDV1",c("subject","Begin_age_weeks")])
#remove subject 451 has no age
input <- input[complete.cases(input), ]
input$Begin_age_weeks <- round(input$Begin_age_weeks)


my.data <-  subjects[subjects$subject %in% as.vector(input$subject),]

mylist <- list()
for (z in 1:nrow(input)){
  x = as.character(input[z,c("subject")])
  d = input[z,c("Begin_age_weeks")]
  m=my.data[my.data$subject == x,]
  m=m[as.numeric(paste(1:d)),]
  mylist[[z]] = m
}

feedingmode_mysamples <- do.call(rbind, mylist)


df <- as.data.frame(summarise(feedingmode_mysamples,age_weeks= n(),bf= sum(bf,na.rm = T),formula=sum(formula,na.rm = T),solid_food=sum(solid_food,na.rm = T),type_solid_food=n_distinct(type_solid_food,na.rm = TRUE)))
df <- df[complete.cases(df), ]
rownames(df) <- df$subject

meta <- my.metadata[my.metadata$groupcode == "KDV1",]

all.data1 <- merge(x = df,y = meta,by = "subject",all.x = T)


########BEWARE I CHANGED the solid food column. removed dots and a comma and also 0
all.data1$feeding_type[all.data1$bf > 0 & all.data1$formula == 0 & all.data1$solid_food==0 & all.data1$type_solid_food == 0] <- "BF" 
all.data1$feeding_type[all.data1$bf > 0 & all.data1$formula > 0 & all.data1$solid_food==0 & all.data1$type_solid_food == 0] <- "MF"
all.data1$feeding_type[all.data1$bf == 0 & all.data1$formula > 0 & all.data1$solid_food==0 & all.data1$type_solid_food == 0] <- "IF"
all.data1$feeding_type[all.data1$solid_food>0 | !(all.data1$type_solid_food == 0) ] <- "MF+SF"
all.data1$feeding_type[all.data1$bf == 0 & all.data1$formula == 0 & all.data1$solid_food==0] <- "no_data"
rownames(all.data1) <- all.data1$Var.2

pre <- all.data1

####################create table with feedin mode up to that time
all.data <- rbind(pre,post)

##############add real age
all.data[all.data$groupcode =="KDV5",54] <- all.data[all.data$groupcode =="KDV5",]$ExactageCCplus28
all.data[all.data$groupcode =="KDV1",54] <- all.data[all.data$groupcode =="KDV1",]$ExactageCCmin2d

colnames(all.data)[54] <- "real_age"
all.data$real_age <- as.numeric(all.data$real_age)
################

## ---- echo=FALSE, fig.height=7, fig.width=12, message=FALSE, warning=FALSE----
library(plyr)
detach("package:plyr", unload=TRUE)


dat1 <- all.data[resulting_samples,]
dat2 <- L2[,resulting_samples]

dat <- cbind(dat1,t(dat2))

dat <- dat[complete.cases(dat$childcarecenter),]
dat <- dat[dat$feeding_type != "no_data",]

library("gdata")
dat <- drop.levels(dat)

dat$age_weeks <- as.numeric(dat$age_weeks)
dat$time <- as.factor(dat$time)


diversity_input <- dat
rownames(diversity_input) <- diversity_input$Var.2 
diversity_input$Var.2 <- NULL

theme_set(theme_bw())
library(ggplot2)
library(ggsignif)
p <- ggplot(diversity_input, aes(x = feeding_type, y = diversity) )
p <- p + geom_boxplot(aes(color = feeding_type))
p <- p + geom_point(aes(color=feeding_type),position = position_jitterdodge())
p <- p + facet_wrap("time")
comparisons <- combn(levels(as.factor(diversity_input$feeding_type)),2,simplify = F)
p <- p + geom_signif(comparisons = comparisons,map_signif_level=F,step_increase=0.05,test = "wilcox.test",tip_length = 0.01,textsize = 3,vjust = 0.4)
p

p <- ggplot(diversity_input, aes(x = feeding_type, y = diversity) )
p <- p + geom_boxplot(aes(color = feeding_type))
p <- p + geom_point(aes(color=feeding_type),position = position_jitterdodge())
p <- p + facet_grid(time~childcarecenter)
comparisons <- combn(levels(as.factor(diversity_input$feeding_type)),2,simplify = F)
p <- p + geom_signif(comparisons = comparisons,map_signif_level=F,step_increase=0.05,test = "wilcox.test",tip_length = 0.01,textsize = 3,vjust = 0.4)
p

p <- ggplot(diversity_input, aes(x = feeding_type, y = diversity.invsimpson) )
p <- p + geom_boxplot(aes(color = feeding_type))
p <- p + geom_point(aes(color=feeding_type),position = position_jitterdodge())
p <- p + facet_wrap("time")
comparisons <- combn(levels(as.factor(diversity_input$feeding_type)),2,simplify = F)
p <- p + geom_signif(comparisons = comparisons,map_signif_level=F,step_increase=0.05,test = "wilcox.test",tip_length = 0.01,textsize = 3,vjust = 0.4)
p

p <- ggplot(diversity_input, aes(x = feeding_type, y = diversity.invsimpson) )
p <- p + geom_boxplot(aes(color = feeding_type))
p <- p + geom_point(aes(color=feeding_type),position = position_jitterdodge())
p <- p + facet_grid(time~childcarecenter)
comparisons <- combn(levels(as.factor(diversity_input$feeding_type)),2,simplify = F)
p <- p + geom_signif(comparisons = comparisons,map_signif_level=F,step_increase=0.05,test = "wilcox.test",tip_length = 0.01,textsize = 3,vjust = 0.4)
p

p <- ggplot(diversity_input, aes(x = feeding_type, y = richness) )
p <- p + geom_boxplot(aes(color = feeding_type))
p <- p + geom_point(aes(color=feeding_type),position = position_jitterdodge())
p <- p + facet_wrap("time")
comparisons <- combn(levels(as.factor(diversity_input$feeding_type)),2,simplify = F)
p <- p + geom_signif(comparisons = comparisons,map_signif_level=F,step_increase=0.05,test = "wilcox.test",tip_length = 0.01,textsize = 3,vjust = 0.4)
p

p <- ggplot(diversity_input, aes(x = feeding_type, y = richness) )
p <- p + geom_boxplot(aes(color = feeding_type))
p <- p + geom_point(aes(color=feeding_type),position = position_jitterdodge())
p <- p + facet_grid(time~childcarecenter)
comparisons <- combn(levels(as.factor(diversity_input$feeding_type)),2,simplify = F)
p <- p + geom_signif(comparisons = comparisons,map_signif_level=F,step_increase=0.05,test = "wilcox.test",tip_length = 0.01,textsize = 3,vjust = 0.4)
p

## ---- echo=FALSE, fig.height=7, fig.width=12, message=FALSE, warning=FALSE----
###########################################
L2_rel_phyloseq_cut <- L2_rel[abundant_taxa,]
L2_rel_phyloseq_all <- L2_rel

library(readr)
tax_table <- read_delim(here("data/hitchip/tax.table.txt"),"\t", escape_double = FALSE, trim_ws = TRUE)
tax_table <- data.frame(tax_table)
rownames(tax_table) <- tax_table$Genus
L2_rel_phyloseq_all

data(dietswap) 
dietswap@sam_data <- sample_data(all.data)
dietswap@otu_table <- otu_table(L2_rel_phyloseq_cut[,rownames(all.data)], taxa_are_rows = TRUE)
dietswap@tax_table <- tax_table(as.matrix(tax_table))
dietswap@otu_table
ps1_cut <- dietswap

data(dietswap) 
dietswap@sam_data <- sample_data(all.data)
dietswap@otu_table <- otu_table(L2_rel_phyloseq_all[,rownames(all.data)], taxa_are_rows = TRUE)
dietswap@tax_table <- tax_table(as.matrix(tax_table))

ps1_all <- dietswap
############################################
theme_set(theme_bw())

set.seed(666)

ordu.wt.uni_all = ordinate(ps1_all, "NMDS", "bray")

z <- plot_ordination(ps1_all, ordu.wt.uni_all)
# Oh no, the table wasn't ordered
library("data.table")
newtab = data.table(z$data)
setorder(newtab,subject,time)
z$data <- newtab
p <- ggplot(data = z$data,aes(x = NMDS1,y = NMDS2))
#p <- p + scale_fill_brewer(type="qual", palette="Set3")
#p <- p + scale_colour_brewer(type="qual", palette="Set3")
p <- p + ggtitle("NMDS bray, all taxa")
p <- p + geom_point(size = 2,alpha=0.3)
#p <- p + stat_ellipse(type = "norm", linetype = 2, level = 0.95)
#p <- p + geom_text(aes(label=subject,colour=feeding_type),size=3,hjust=-0.5, vjust=0)
#p <- p + geom_segment(aes(x=x_coord_start,y=y_coord_start, col=breastfeeding, xend = x_coord_end, yend = y_coord_end),size =0.5, arrow=arrow(length=unit(0.1,"cm")))
#p <- p + geom_text(aes(label=groupcode,colour=groupcode),size=3,hjust=-0.5, vjust=0)
p <- p + geom_path(aes(group=subject, col=feeding_type),alpha=0.3,size =1,arrow=arrow(length=unit(0.35,"cm")),lineend = "butt",linemitre = 10)
p <- p +facet_grid(feeding_type~childcarecenter)
print(p)



z <- plot_ordination(ps1_all, ordu.wt.uni_all)
# Oh no, the table wasn't ordered
library("data.table")
newtab = data.table(z$data)
setorder(newtab,subject,time)
z$data <- newtab
p <- ggplot(data = z$data,aes(x = NMDS1,y = NMDS2))
#p <- p + scale_fill_brewer(type="qual", palette="Set3")
#p <- p + scale_colour_brewer(type="qual", palette="Set3")
p <- p + ggtitle("NMDS bray, all taxa")
p <- p + geom_point(size = 2,alpha=0.3)
#p <- p + stat_ellipse(type = "norm", linetype = 2, level = 0.95)
#p <- p + geom_text(aes(label=subject,colour=feeding_type),size=3,hjust=-0.5, vjust=0)
#p <- p + geom_segment(aes(x=x_coord_start,y=y_coord_start, col=breastfeeding, xend = x_coord_end, yend = y_coord_end),size =0.5, arrow=arrow(length=unit(0.1,"cm")))
#p <- p + geom_text(aes(label=groupcode,colour=groupcode),size=3,hjust=-0.5, vjust=0)
p <- p + geom_path(aes(group=subject, col=feeding_type),alpha=0.3,size =1,arrow=arrow(length=unit(0.35,"cm")),lineend = "butt",linemitre = 10)
p1 <- p + facet_grid("feeding_type")

ordu.wt.uni_cut = ordinate(ps1_cut, "NMDS", "bray")
z <- plot_ordination(ps1_all, ordu.wt.uni_cut)
# Oh no, the table wasn't ordered
library("data.table")
newtab = data.table(z$data)
setorder(newtab,subject,time)
z$data <- newtab
p <- ggplot(data = z$data,aes(x = NMDS1,y = NMDS2))
#p <- p + scale_fill_brewer(type="qual", palette="Set3")
#p <- p + scale_colour_brewer(type="qual", palette="Set3")
p <- p + ggtitle("NMDS bray, only abundant taxa")
p <- p + geom_point(size = 2,alpha=0.3)
#p <- p + stat_ellipse(type = "norm", linetype = 2, level = 0.95)
#p <- p + geom_text(aes(label=subject,colour=feeding_type),size=3,hjust=-0.5, vjust=0)
#p <- p + geom_segment(aes(x=x_coord_start,y=y_coord_start, col=breastfeeding, xend = x_coord_end, yend = y_coord_end),size =0.5, arrow=arrow(length=unit(0.1,"cm")))
#p <- p + geom_text(aes(label=groupcode,colour=groupcode),size=3,hjust=-0.5, vjust=0)
p <- p + geom_path(aes(group=subject, col=feeding_type),alpha=0.3,size =1,arrow=arrow(length=unit(0.35,"cm")),lineend = "butt",linemitre = 10)
p2 <- p + facet_grid("feeding_type")


p <- ggplot(data = z$data,aes(x = NMDS1,y = NMDS2))
#p <- p + scale_fill_brewer(type="qual", palette="Set3")
#p <- p + scale_colour_brewer(type="qual", palette="Set3")
p <- p + ggtitle("NMDS bray, just playing density plot")
p <- p + stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200) +
  scale_fill_continuous(low = "white", high = "red")
p <- p + geom_point(size = 2,alpha=0.3)
#p <- p + stat_ellipse(type = "norm", linetype = 2, level = 0.95)
#p <- p + geom_text(aes(label=subject,colour=feeding_type),size=3,hjust=-0.5, vjust=0)
#p <- p + geom_segment(aes(x=x_coord_start,y=y_coord_start, col=breastfeeding, xend = x_coord_end, yend = y_coord_end),size =0.5, arrow=arrow(length=unit(0.1,"cm")))
#p <- p + geom_text(aes(label=groupcode,colour=groupcode),size=3,hjust=-0.5, vjust=0)
p <- p + geom_path(aes(group=subject, col=feeding_type),alpha=0.3,size =1,arrow=arrow(length=unit(0.35,"cm")),lineend = "butt",linemitre = 10)
p <- p + facet_grid("feeding_type")
print(p)


#with arrows
p <- plot_ordination(ps1_all, ordu.wt.uni_all,type = "taxa",label="Genus", color="Order")
p3 <- p + ggtitle("NMDS bray, all taxa")


p <- plot_ordination(ps1_cut, ordu.wt.uni_cut,type = "taxa",label="Genus", color="Order")
p4 <- p + ggtitle("NMDS bray, only abundant taxa")


library("cowplot")
plot_grid(p1 + theme(legend.position="none"), p2,p3 + theme(legend.position="none"),p4, labels = c("A", "B", "C", "D"), ncol = 2,nrow = 2 )


library(vegan)
loadings <- scores(ordu.wt.uni_all,display=c("species"))
top.coef <- loadings[,1][rev(order(abs(loadings[,1])))][1:30]
names(top.coef) <- gsub("Clostridium_sensu_stricto","Clostridium_(sensu_stricto)", names(top.coef))
top.coef <- cbind(top.coef,as.data.frame(ps1_all@tax_table)[names(top.coef),])


library(RColorBrewer)
cols <- brewer.pal(5,"Set3")
p <- ggplot(top.coef, aes(x=reorder(Genus, top.coef),y=top.coef))
p <- p + scale_fill_brewer(palette="Set1")
p <- p + geom_bar(stat = "identity",aes(fill=Phylum),colour="black")
p <- p + labs(y = "NMDS1 score",x = NULL)
p5 <- p + coord_flip()



library(vegan)
loadings <- scores(ordu.wt.uni_cut,display=c("species"))
top.coef <- loadings[,1][rev(order(abs(loadings[,1])))][1:nrow(loadings)]
names(top.coef) <- gsub("Clostridium_sensu_stricto","Clostridium_(sensu_stricto)", names(top.coef))
top.coef <- cbind(top.coef,as.data.frame(ps1_cut@tax_table)[names(top.coef),])

## ---- echo=FALSE, fig.height=3.5, fig.width=14, message=FALSE, warning=FALSE----
library(RColorBrewer)
cols <- brewer.pal(5,"Set3")
p <- ggplot(top.coef, aes(x=reorder(Genus, top.coef),y=top.coef))
p <- p + scale_fill_brewer(palette="Set1")
p <- p + geom_bar(stat = "identity",aes(fill=Order),colour="black")
p <- p + labs(y = "NMDS1 score",x = NULL)
p6 <- p + coord_flip()

library("cowplot")
plot_grid(p5 , p6,  
          labels = c("A", "B"),
          ncol = 2 )

## ---- echo=FALSE, fig.height=8, fig.width=12, message=FALSE, warning=FALSE----
pre_KDV <- all.data[all.data$groupcode == "KDV1",][,c("Var.2","subject")]
post_KDV <- all.data[all.data$groupcode == "KDV5",][,c("Var.2","subject")]

matching_subjects <- Reduce(intersect, list(pre_KDV$subject,post_KDV$subject))

selected_samples <- as.character(subset(all.data, subject %in% matching_subjects)[,"Var.2"])

resulting_samples2 <- all.data[selected_samples,]


subject_list <- split(as.character(resulting_samples2$Var.2),as.character(resulting_samples2$subject), drop= TRUE)

myd <- L1_rel[,rownames(resulting_samples2)]

mylist <- list()
for(i in 1:length(subject_list)){
  x=subject_list[[names(subject_list)[i]]]
  y=as.character(x)
  temp_matrix=myd[,y]
  mylist[[i]] <- temp_matrix
}

names(mylist) <- names(subject_list)

values <- list()
for(i in 1:length(subject_list)){
  x=names(mylist)[i]
  pre_pro=mylist[[x]][,2] - mylist[[x]][,1]
  values[[i]] <- pre_pro
}

names(values) <- names(subject_list)

values <- as.data.frame(values)

values.matrix <- as.matrix(values)

#values.matrix <- 10^(values.matrix)
matrix_pre_pro = values.matrix

input <- matrix_pre_pro

#rownames(input) <- df$taxa
#input_ordered <- input[order(rownames(input)),]
#######################################################
library(gplots)
#heatmap.2(input, main="",col=bluered(100),
#          symbreaks=TRUE,trace="none",margin=c(5,13),Rowv=T,
#          Colv=T,cexRow=0.8,cexCol=0.8,keysize = 1,density.info="none",key.title="NA")

#heatmap.2(input, main="",col=bluered(100),
#          symbreaks=TRUE,trace="none",margin=c(5,13),Rowv=T,scale = "row",
#          Colv=T,cexRow=0.8,cexCol=0.8,keysize = 1,density.info="none",key.title="NA")

############################
library(RColorBrewer)
col_groups <- as.factor(resulting_samples2$feeding_type)[which(resulting_samples2$groupcode == "KDV5")]

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- resulting_samples2$subject[which(resulting_samples2$groupcode == "KDV5")]

mat_colors <- list(group = brewer.pal(length(unique(col_groups)), "Set1"))
names(mat_colors$group) <- unique(col_groups)


#min(scale(myd[abundant_taxa,])
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 29)
breaks = c(seq(-15,-7,length=10),seq(-4,4,length=10),seq(5,15,length=6))

library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

d <- stats::dist(t(input[c("Actinobacteria","Bacilli","Bacteroidetes","Proteobacteria","Clostridium_cluster_IV","Clostridium_cluster_XIVa"),]),method = 'manhattan' )
hpws <- hclust(d, method = "mcquitty")
hpws_ordered <- sort_hclust(hpws)

library(pheatmap)
pheatmap(input[,],color = my_palette, clustering_distance_cols = 'minkowski',
         clustering_distance_rows = 'maximum',cutree_cols = 3, annotation_col= mat_col,annotation_colors = mat_colors, main = "try whether clustering works")

## ---- echo=FALSE, fig.height=10, fig.width=10, message=FALSE, warning=FALSE----
#r <- c("Actinobacteria","Bacilli","Bacteroidetes","Proteobacteria","Clostridium_cluster_IV","Clostridium_cluster_XIVa")
r <- c("Actinobacteria","Bacilli","Proteobacteria","Bacteroidetes")


myd <- L1_rel[,rownames(resulting_samples2)]


# calculate which are the most abundant states in the subgroups and visualise with heatmap

library(arules)


par(mfrow=c(3,3))
mylist <- list()
for (i in 1:length(r)){
  taxa = r[i]
  plot(hist(input[c(taxa),],breaks = 20), main = paste(taxa,"equal frequency")) +
    abline(v = discretize(input[c(taxa),], method = "interval", breaks = 3, 
                          onlycuts = TRUE,include.lowest = T), col = "red")
  mylist[[taxa]] <- as.data.frame(discretize(input[c(taxa),], "interval",breaks=3,labels=c(paste("low",taxa, sep="."),paste("no",taxa, sep="."),paste("high",taxa, sep="."))))
}



bimod <- do.call(cbind,mylist) 
colnames(bimod) <- r
rownames(bimod) <- colnames(input[])

#combine the data of all grouping into 1 barcode
mylist <- c()
for (i in 1:nrow(bimod)){
  x=paste(bimod[c(i),], collapse= "")
  mylist[i]=x
}
names(mylist)<- rownames(bimod)


Tipping.elements <- as.data.frame(table(mylist))
Tipping.elements <- Tipping.elements[order(Tipping.elements$Freq, decreasing=TRUE),]
Tipping.elements <- cbind((Tipping.elements$Freq/sum(Tipping.elements$Freq)*100),Tipping.elements)
colnames(Tipping.elements) <- c("frequency %","barcode","Freq")
Tipping.elements$barcode <- as.character(Tipping.elements$barcode)
print(Tipping.elements)



mylist2 <- list()
for (i in 1:length(r)){
  taxa = r[i]
  minima <- discretize(input[c(taxa),], "frequency",breaks=3,onlycuts = TRUE)[2]
  mylist2[[taxa]] <- as.matrix(input[paste(taxa),])
  #in the original function
  #mylist2[[taxa]] <- as.matrix(input[paste(taxa),]-minima)
}

bimod.scores <- do.call(cbind,mylist2) 
colnames(bimod.scores) <- r

for.plot <- as.data.frame(cbind(mylist,bimod.scores))

ordering <- list()
for (i in 1:nrow(Tipping.elements)){
  x=subset(for.plot, mylist == Tipping.elements$barcode[i])
  ordering[[i]]=x
}
names(ordering)<- Tipping.elements$barcode

library(plyr)
df <- rbind.fill(ordering)
df <- cbind(rownames(bimod),df)

#turn the factors into numeric
dat1 <- data.frame(lapply(df[,3:ncol(df)], function(x) as.numeric(as.character(x))))
df <- cbind(df[,1:2],dat1)

## ---- echo=FALSE, fig.height=2.5, fig.width=10, message=FALSE, warning=FALSE----
library(ggplot2)
library(reshape)

input.sampleID <- melt(df, id.vars='mylist', measure.vars=c("rownames(bimod)"))
input_bc <- melt(df, id.vars='mylist', measure.vars=c(colnames(df[,3:ncol(df)])))

input_bc <- cbind(input_bc,input.sampleID[,c("value")])
colnames(input_bc) <- c("barcode","variable","value","Subject")

#make the plot
p <- ggplot(input_bc, aes(x = Subject,y = variable,fill=value))
#p <- p + scale_fill_gradientn(colors = c("red","white","blue"),values=c(1,0.55,0.5,0.45,0),limits = c(-1, 1))
p <- p + scale_fill_gradient2(low="blue",high= "red" ,guide = "colorbar")
p <- p + geom_tile()
p <- p + ggtitle("Change in bacterial groups")
p <- p + theme(axis.text.x=element_text(angle = -90, hjust = 0), axis.title.y = element_blank())
p


#---------------------------------------------------

barcodes <- df
colnames(barcodes)[1:2] <- c("subject","barcode")
barcodes$subject <- as.character(gsub(x = barcodes$subject,pattern = "X", replacement = ""))

barcode.metadata <- merge(barcodes,all.data[rownames(resulting_samples2),], by="subject",  all.x=TRUE)
barcode.metadata$barcode <- as.factor(barcode.metadata$barcode)


xx <- barcode.metadata[barcode.metadata$childcarecenter == "CC" & barcode.metadata$groupcode == "KDV5",c(colnames(barcodes))]
xx <- xx[complete.cases(xx),]


input.sampleID <- melt(xx, id.vars='barcode', measure.vars=c("subject"))
input_bc1 <- melt(xx, id.vars='barcode', measure.vars=c(colnames(df[,c(3:ncol(xx))])))
input_bc1 <- cbind(input_bc1,input.sampleID[,c("value")])
colnames(input_bc1) <- c("barcode","variable","value","Subject")
input_bc1$Subject <- as.factor(input_bc1$Subject)
input_bc1$value <- as.numeric(input_bc1$value)


#make the plot
p <- ggplot(input_bc1, aes(x = Subject,y = variable,fill=value))
#p <- p + scale_fill_gradientn(colors = c("red","white","blue"),values=c(1,0.55,0.5,0.45,0),limits = c(-1, 1))
p <- p + scale_fill_gradient2(low="blue",high= "red" ,guide = "colorbar")
p <- p + geom_tile()
p <- p + ggtitle("Change in bacterial groups childcare")
p <- p + theme(axis.text.x=element_text(angle = -90, hjust = 0), axis.title.y = element_blank())
p


#---------------------------------------------------------------------------

xx <- barcode.metadata[barcode.metadata$childcarecenter == "noCC" & barcode.metadata$groupcode == "KDV5",c(colnames(barcodes))]
xx <- xx[complete.cases(xx),]


input.sampleID <- melt(xx, id.vars='barcode', measure.vars=c("subject"))
input_bc1 <- melt(xx, id.vars='barcode', measure.vars=c(colnames(df[,c(3:ncol(xx))])))
input_bc1 <- cbind(input_bc1,input.sampleID[,c("value")])
colnames(input_bc1) <- c("barcode","variable","value","Subject")
input_bc1$Subject <- as.factor(input_bc1$Subject)
input_bc1$value <- as.numeric(input_bc1$value)


#make the plot
p <- ggplot(input_bc1, aes(x = Subject,y = variable,fill=value))
#p <- p + scale_fill_gradientn(colors = c("red","white","blue"),values=c(1,0.55,0.5,0.45,0),limits = c(-1, 1))
p <- p + scale_fill_gradient2(low="blue",high= "red" ,guide = "colorbar")
p <- p + geom_tile()
p <- p + ggtitle("Change in bacterial groups no childcare")
p <- p + theme(axis.text.x=element_text(angle = -90, hjust = 0), axis.title.y = element_blank())
p

#########
xx <- barcode.metadata[barcode.metadata$feeding_type == "BF" & barcode.metadata$groupcode == "KDV5",c(colnames(barcodes))]
xx <- xx[complete.cases(xx),]


input.sampleID <- melt(xx, id.vars='barcode', measure.vars=c("subject"))
input_bc1 <- melt(xx, id.vars='barcode', measure.vars=c(colnames(df[,c(3:ncol(xx))])))
input_bc1 <- cbind(input_bc1,input.sampleID[,c("value")])
colnames(input_bc1) <- c("barcode","variable","value","Subject")
input_bc1$Subject <- as.factor(input_bc1$Subject)
input_bc1$value <- as.numeric(input_bc1$value)


#make the plot
p <- ggplot(input_bc1, aes(x = Subject,y = variable,fill=value))
#p <- p + scale_fill_gradientn(colors = c("red","white","blue"),values=c(1,0.55,0.5,0.45,0),limits = c(-1, 1))
p <- p + scale_fill_gradient2(low="blue",high= "red" ,guide = "colorbar")
p <- p + geom_tile()
p <- p + ggtitle("Change in bacterial groups BF")
p <- p + theme(axis.text.x=element_text(angle = -90, hjust = 0), axis.title.y = element_blank())
p


#########
xx <- barcode.metadata[barcode.metadata$feeding_type == "IF" & barcode.metadata$groupcode == "KDV5",c(colnames(barcodes))]
xx <- xx[complete.cases(xx),]


input.sampleID <- melt(xx, id.vars='barcode', measure.vars=c("subject"))
input_bc1 <- melt(xx, id.vars='barcode', measure.vars=c(colnames(df[,c(3:ncol(xx))])))
input_bc1 <- cbind(input_bc1,input.sampleID[,c("value")])
colnames(input_bc1) <- c("barcode","variable","value","Subject")
input_bc1$Subject <- as.factor(input_bc1$Subject)
input_bc1$value <- as.numeric(input_bc1$value)


#make the plot
p <- ggplot(input_bc1, aes(x = Subject,y = variable,fill=value))
#p <- p + scale_fill_gradientn(colors = c("red","white","blue"),values=c(1,0.55,0.5,0.45,0),limits = c(-1, 1))
p <- p + scale_fill_gradient2(low="blue",high= "red" ,guide = "colorbar")
p <- p + geom_tile()
p <- p + ggtitle("Change in bacterial groups IF")
p <- p + theme(axis.text.x=element_text(angle = -90, hjust = 0), axis.title.y = element_blank())
p


#########
xx <- barcode.metadata[barcode.metadata$feeding_type == "MF" & barcode.metadata$groupcode == "KDV5",c(colnames(barcodes))]
xx <- xx[complete.cases(xx),]


input.sampleID <- melt(xx, id.vars='barcode', measure.vars=c("subject"))
input_bc1 <- melt(xx, id.vars='barcode', measure.vars=c(colnames(df[,c(3:ncol(xx))])))
input_bc1 <- cbind(input_bc1,input.sampleID[,c("value")])
colnames(input_bc1) <- c("barcode","variable","value","Subject")
input_bc1$Subject <- as.factor(input_bc1$Subject)
input_bc1$value <- as.numeric(input_bc1$value)


#make the plot
p <- ggplot(input_bc1, aes(x = Subject,y = variable,fill=value))
#p <- p + scale_fill_gradientn(colors = c("red","white","blue"),values=c(1,0.55,0.5,0.45,0),limits = c(-1, 1))
p <- p + scale_fill_gradient2(low="blue",high= "red" ,guide = "colorbar")
p <- p + geom_tile()
p <- p + ggtitle("Change in bacterial groups MF")
p <- p + theme(axis.text.x=element_text(angle = -90, hjust = 0), axis.title.y = element_blank())
p

## ---- echo=FALSE, fig.height=6, fig.width=8, message=FALSE, warning=FALSE----
library(plyr)
detach("package:plyr", unload=TRUE)

dat1 <- all.data[resulting_samples,]
dat2 <- L2[,resulting_samples]

dat <- cbind(dat1,t(dat2))

dat <- dat[complete.cases(dat$childcarecenter),]
dat <- dat[dat$feeding_type != "no_data" & dat$feeding_type != "MF+SF",]


library("gdata")
dat <- drop.levels(dat)

dat$age_weeks <- as.numeric(dat$age_weeks)
dat$time <- as.factor(dat$time)

library(dplyr)
dat <- 
  dat %>% dplyr::select(
    rownames(dat2), 
    subject,
    Var.2,
    age_weeks, 
    feeding_type, 
    time, 
    childcarecenter
  ) %>% na.omit() %>%
  mutate(
    time = droplevels(time),
    age_weeks.c = age_weeks - mean(age_weeks)
  )

group_by(dat, feeding_type) %>% summarise(n = n())

# how many oberservation per subgroup left?
group_by(dat, childcarecenter, time, feeding_type) %>%
  summarise(mean = mean(L2[names(mylist[1])]), n = n())

model_input <- dat
rownames(model_input) <- model_input$Var.2 
model_input$Var.2 <- NULL


library(vegan)
childcarecenter <- model_input$childcarecenter
feeding_type <- as.factor(model_input$feeding_type)
time <- model_input$time
real_age <- as.factor(model_input$real_age)
response <- model_input[,1:nrow(L2)]


#mod <- prc(model_input[,1:130], real_age , time)
mod <- prc(response, feeding_type , time)
anova.cca(mod, strata = time, first=TRUE)
#anova.cca(mod,by="terms", strata = time, first=TRUE)
RsquareAdj(mod)
sum_mod <- summary(mod) #PRC
p <- plot(mod, species = TRUE, select = abs(sum_mod$sp) > 0.5 , scaling = 2, axis = 1, type = "l", 
     ylim, lty = 1:5, col = 1:6, cex = 0.8,main= "Euclidean based PRC") 
print(p)
################

#############################
##############################
x <- deparse(substitute(feeding_type))
z <- deparse(substitute(time))

##############################
fla <- as.formula(paste("~", x, "+", z))
mf <- model.frame(fla, response, na.action = na.pass)
fla.zx <- as.formula(paste("~", z, ":", x))
fla.z <- as.formula(paste("~", z))
X = model.matrix(fla.zx, mf)[, -c(seq_len(nlevels(time) + 1))]
Z = model.matrix(fla.z, mf)[, -1]

d <- as.matrix(vegdist(response, method = "bray"))
dbRDA <- capscale(d ~ X + Condition(Z), data = model_input[(nrow(L2)+1):ncol(model_input)], comm = response)

dbRDA$terminfo$xlev = list(levels(time), levels(feeding_type))
names(dbRDA$terminfo$xlev) = c(paste(z), paste(x))
dbRDA$call <- match.call()
class(dbRDA) <- c("prc", class(dbRDA))
sum_dbRDA <- summary(mod)

anova.cca(mod, strata = time, first=TRUE)
#anova.cca(mod,by="terms", strata = time, first=TRUE)
RsquareAdj(mod)
p <- plot(dbRDA, species = TRUE, select = abs(sum_dbRDA$sp) > 0.5 , scaling = 2, axis = 1, type = "l", 
     ylim = c(-0.25,1), lty = 1:5, col = 1:6, cex = 0.8, xlab = "Age (weeks)",main= "Bray Curtis based PRC") 
print(p)
response
## ---- echo=FALSE, fig.height=8, fig.width=10, message=FALSE, warning=FALSE----
RDA <- rda(response ~ childcarecenter * feeding_type * time + Condition(time))
anova.cca(RDA,by = "terms", strata = time)
RDA

library(ggvegan)
autoplot(RDA, arrows = TRUE, geom = c("point", "text"),
         layers = c("species", "sites", "biplot", "centroids"),
         legend.position = "right")


d <- as.matrix(vegdist(response, method = "bray"))
dbRDA <- capscale(d ~ childcarecenter * feeding_type * time + Condition(time))
anova.cca(dbRDA,by = "terms", strata = time)
dbRDA

#doesn't work
#library(ggvegan)
#autoplot(dbRDA, arrows = TRUE, geom = c("point", "text"),
#         layers = c("species", "sites", "biplot", "centroids"),
#         legend.position = "right")

## ----create models and pictures L1, echo=FALSE, fig.height=5, fig.width=10, message=FALSE, warning=FALSE----
library(plyr)
detach("package:plyr", unload=TRUE)

dat1 <- all.data[resulting_samples,]
dat2 <- L1[,resulting_samples]

dat <- cbind(dat1,t(dat2))

dat <- dat[complete.cases(dat$childcarecenter),]
dat <- dat[dat$feeding_type != "no_data" & dat$feeding_type != "MF+SF",]


library("gdata")
dat <- drop.levels(dat)

dat$age_weeks <- as.numeric(dat$age_weeks)
dat$time <- as.factor(dat$time)

library(dplyr)
dat <- 
  dat %>% dplyr::select(
    rownames(dat2), 
    subject,
    Var.2,
    age_weeks, 
    feeding_type, 
    time, 
    childcarecenter
  ) %>% na.omit() %>%
  mutate(
    time = droplevels(time),
    age_weeks.c = age_weeks - mean(age_weeks)
  )

group_by(dat, feeding_type) %>% summarise(n = n())

# how many oberservation per subgroup left?
group_by(dat, childcarecenter, time, feeding_type) %>%
  summarise(mean = mean(L2[names(mylist[1])]), n = n())

model_input <- dat
rownames(model_input) <- model_input$Var.2 
model_input$Var.2 <- NULL


model_input %>% group_by(time, childcarecenter) %>% summarise(n = n())

model_input %>% group_by(feeding_type, childcarecenter, time, subject) %>% summarise(n = n()) %>% arrange(subject)

library(afex) 
source(here("R/reporting.R"))


mylist_childcarecenter <- list()
mylist_feeding_type <- list()
mylist_time <- list()
mylist_age_weeks_c <- list()
mylist_childcarecenter.time <- list()
mylist_childcarecenter.feeding_type <- list()
mylist_feeding_type.time <- list()
mylist_childcarecenter.feeding_type.time <- list()
mylist_diagnostic.plots <- list()

mylist <- list()

model_input %>% dim()


for (i in rownames(dat2)){
  dillie <- colnames(model_input[(nrow(dat2)+1):ncol(model_input)])[c(5,3,4)]
  form<- reformulate(paste(paste(dillie, collapse="*"),"+ age_weeks.c + (1|subject)"), i)
  fit<- mixed(form, data = model_input, method = 'KR', type = 2,check_contrasts = TRUE)
  
  if(fit$anova_table[1,4] <=0.05) {mylist_childcarecenter[[i]] <- fit$anova_table[1,4,drop=FALSE]}
  if(fit$anova_table[2,4] <=0.05) {mylist_feeding_type[[i]] <- fit$anova_table[2,4,drop=FALSE]}
  if(fit$anova_table[3,4] <=0.05) {mylist_time[[i]] <- fit$anova_table[3,4,drop=FALSE]}
  if(fit$anova_table[4,4] <=0.05) {mylist_age_weeks_c[[i]] <- fit$anova_table[4,4,drop=FALSE]}
  if(fit$anova_table[5,4] <=0.05) {mylist_childcarecenter.feeding_type[[i]] <- fit$anova_table[5,4,drop=FALSE]}
  if(fit$anova_table[6,4] <=0.05) {mylist_childcarecenter.time[[i]] <- fit$anova_table[6,4,drop=FALSE]}
  if(fit$anova_table[7,4] <=0.05) {mylist_feeding_type.time[[i]] <- fit$anova_table[7,4,drop=FALSE]}
  if(fit$anova_table[8,4] <=0.05) {mylist_childcarecenter.feeding_type.time[[i]] <- fit$anova_table[8,4,drop=FALSE]}

  if(sum(fit$anova_table$`Pr(>F)` <=0.05) >0 ) {
    diag_df <- 
      tibble( 
        sresid = resid(fit$full_model, scale=T), 
        fitted = fitted(fit$full_model),
        observed = as.numeric(unlist(fit$data[c(paste(i))]))
      )
  # fitted versus standardized residuals
  p1 <- ggplot(diag_df, aes(fitted, sresid)) + geom_point() + geom_smooth(method = 'loess') + ggtitle(paste(i))
  # fitted versus observed
  p2 <- ggplot(diag_df, aes(fitted, observed)) + geom_point() + geom_smooth(method = 'loess') + ggtitle(paste(i))
  # normality of residuals density
  p3 <- ggplot(diag_df, aes(sresid)) +  geom_density() + ggtitle(paste(i))
  # normality of residuals qqplot
  p4 <- gg_qq(diag_df$sresid)
  

  mylist_diagnostic.plots[[i]] <- 
    cowplot::plot_grid(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2)
  
  taxa <- L1_rel[i,dat$Var.2]
  df <- cbind(taxa,dat)
  
  dat.all <- melt(df, id.vars='time', measure.vars=c('taxa'))
  dat.all2 <- melt(df, id.vars='subject', measure.vars=c("taxa"))
  dat.all3 <- melt(df, id.vars='feeding_type', measure.vars=c("taxa"))
  dat.all4 <- melt(df, id.vars='childcarecenter', measure.vars=c("taxa"))
  
  #get rid of all the NA's
  dat.all <- cbind(dat.all,dat.all2[,1],dat.all3[,1],dat.all4[,1])
  dat.all <- dat.all[rowSums(is.na(dat.all)) < 1, ]
  
  dat.all[ order(dat.all[,4], dat.all[,1],decreasing = TRUE), ]
  
  colnames(dat.all) <- c("time","variable","value","subject","feeding_type","childcarecenter")
  dat.all$time <- as.factor(dat.all$time)
  dat.all$subject <- as.factor(dat.all$subject)
  
  
  #make the plot for "group"
  p <- ggplot(dat.all, aes(time, value, fill=feeding_type)) 
  p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
  p <- p + geom_boxplot(alpha=.2, width=0.4)
  p <- p + geom_text(aes(label=subject,colour=subject),size=2,hjust=-0.1, vjust=0) + ylab("Relative abundance")
  p <- p + geom_path(aes(group=subject, col=subject))
  #p <- p + geom_smooth(aes(group=childcarecenter, fill = childcarecenter, alpha=0.5),colour="black", method="loess")
  p <- p + theme_bw()
  p <- p + facet_grid("feeding_type")
  p <- p + ggtitle(paste(i))
  #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
  p1 <- p + theme(legend.position = "none")
  

  p <- ggplot(dat.all, aes(time, value, fill=childcarecenter)) 
  p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
  p <- p + geom_boxplot(alpha=.2, width=0.4)
  p <- p + geom_text(aes(label=subject,colour=subject),size=2,hjust=-0.1, vjust=0) + ylab("Relative abundance")
  p <- p + geom_path(aes(group=subject, col=subject))
  #p <- p + geom_smooth(aes(group=childcarecenter, fill = childcarecenter, alpha=0.5),colour="black", method="loess")
  p <- p + theme_bw()
  p <- p + facet_grid("childcarecenter")
  p <- p + ggtitle(paste(i))
  #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
  p2 <- p + theme(legend.position = "none")
  
  
  print(plot_grid(p1, p2,labels = c("A", "B"), ncol = 
  2))

  
  mylist[[i]] <- anova(fit)
  }
}  

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
l <- list(mylist_childcarecenter,mylist_feeding_type,mylist_time,mylist_age_weeks_c,mylist_childcarecenter.feeding_type,mylist_childcarecenter.time,mylist_feeding_type.time,mylist_childcarecenter.feeding_type.time)
names(l) <- c("mylist_childcarecenter","mylist_feeding_type","mylist_time","mylist_age_weeks_c","mylist_childcarecenter.feeding_type","mylist_childcarecenter.time","mylist_feeding_type.time","mylist_childcarecenter.feeding_type.time")



keys <- unique(unlist(lapply(l, names)))
zz <- setNames(do.call(mapply, c(FUN=c, lapply(l, `[`, keys))), keys)

list_df <- function(list) {
  df <- as.data.frame(list)
  names(df) <- list_names(list)
  return (df)
}

list_names <- function(list) {
  
  recursor <- function(list, names) {
    if (is.list(list)) {
      new_names <- paste(names, names(list), sep = ".")
      out <- unlist(mapply(list, new_names, FUN = recursor))
    } else {
      out <- names
    }
    return(out)
  }
  
  new_names <- unlist(mapply(list, names(list), FUN = recursor))
  return(new_names)
}

df <- as.data.frame(t(list_df(zz)))

df$taxa <- rownames(df)


library(dplyr)
library(tidyr)

df <- df %>%
  separate(taxa, c("taxa", "variable"), ".mylist_")
df$variable <- gsub(pattern = '.{7}$',replacement = '',x = df$variable,)
colnames(df)[1] <- "p.value"
print(df)

## ---- eval=FALSE, message=FALSE, warning=FALSE, include=FALSE------------
## myd <- L1_rel[names(mylist),dat$Var.2]
## 
## all.graphs <- list()
## for (i in rownames(myd)) {
##   taxa <- myd[i,dat$Var.2]
##   df <- cbind(taxa,dat)
## 
##   dat.all <- melt(df, id.vars='time', measure.vars=c('taxa'))
##   dat.all2 <- melt(df, id.vars='subject', measure.vars=c("taxa"))
##   dat.all3 <- melt(df, id.vars='feeding_type', measure.vars=c("taxa"))
##   dat.all4 <- melt(df, id.vars='childcarecenter', measure.vars=c("taxa"))
## 
##   #get rid of all the NA's
##   dat.all <- cbind(dat.all,dat.all2[,1],dat.all3[,1],dat.all4[,1])
##   dat.all <- dat.all[rowSums(is.na(dat.all)) < 1, ]
## 
##   dat.all[ order(dat.all[,4], dat.all[,1],decreasing = TRUE), ]
## 
##   colnames(dat.all) <- c("time","variable","value","subject","feeding_type","childcarecenter")
##   dat.all$time <- as.factor(dat.all$time)
##   dat.all$subject <- as.factor(dat.all$subject)
## 
## 
##   #make the plot for "group"
##   p <- ggplot(dat.all, aes(time, value, fill=feeding_type))
##   p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
##   p <- p + geom_boxplot(alpha=.2, width=0.4)
##   p <- p + geom_text(aes(label=subject,colour=subject),size=2,hjust=-0.1, vjust=0) + ylab("Relative abundance")
##   p <- p + geom_path(aes(group=subject, col=subject))
##   #p <- p + geom_smooth(aes(group=childcarecenter, fill = childcarecenter, alpha=0.5),colour="black", method="loess")
##   p <- p + theme_bw()
##   p <- p + facet_grid("feeding_type")
##   p <- p + ggtitle(paste(i))
##   #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
##   p <- p + theme(legend.position = "none")
##   p1 <- p
## 
##   all.graphs[[i]] <- p1
## 
## }
## 
## all.graphs

## ----create models and pictures L2, echo=FALSE, fig.height=8, fig.width=10, message=FALSE, warning=FALSE----
library(plyr)
detach("package:plyr", unload=TRUE)

dat1 <- all.data[resulting_samples,]
dat2 <- L2[,resulting_samples]

dat <- cbind(dat1,t(dat2))

dat <- dat[complete.cases(dat$childcarecenter),]
dat <- dat[dat$feeding_type != "no_data" & dat$feeding_type != "MF+SF",]


library("gdata")
dat <- drop.levels(dat)

dat$age_weeks <- as.numeric(dat$age_weeks)
dat$time <- as.factor(dat$time)

library(dplyr)
dat <- 
  dat %>% dplyr::select(
    rownames(dat2), 
    subject,
    Var.2,
    age_weeks, 
    feeding_type, 
    time, 
    childcarecenter
  ) %>% na.omit() %>%
  mutate(
    time = droplevels(time),
    age_weeks.c = age_weeks - mean(age_weeks)
  )

group_by(dat, feeding_type) %>% summarise(n = n())

# how many oberservation per subgroup left?
group_by(dat, childcarecenter, time, feeding_type) %>%
  summarise(mean = mean(L2[names(mylist[1])]), n = n())

model_input <- dat
rownames(model_input) <- model_input$Var.2 
model_input$Var.2 <- NULL

library(afex)
source(here("R/reporting.R"))

mylist_childcarecenter <- list()
mylist_feeding_type <- list()
mylist_time <- list()
mylist_childcarecenter.time <- list()
mylist_childcarecenter.feeding_type <- list()
mylist_feeding_type.time <- list()
mylist_childcarecenter.feeding_type.time <- list()
mylist_diagnostic.plots <- list()
mylist <- list()
all.models <- list()
for (i in rownames(dat2)){
  dillie <- colnames(model_input[(nrow(dat2)+1):ncol(model_input)])[c(5,3,4)]
  form<- reformulate(paste(paste(dillie, collapse="*"),"+ age_weeks.c + (1|subject)"), i)
  fit<- mixed(form, data = model_input, method = 'KR', type = 2,check_contrasts = TRUE)
  # store all models in a list
  all.models[[i]] <- fit
  
  if(fit$anova_table[1,4] <=0.05) {mylist_childcarecenter[[i]] <- fit$anova_table[1,4,drop=FALSE]}
  if(fit$anova_table[2,4] <=0.05) {mylist_feeding_type[[i]] <- fit$anova_table[2,4,drop=FALSE]}
  if(fit$anova_table[3,4] <=0.05) {mylist_time[[i]] <- fit$anova_table[3,4,drop=FALSE]}
  if(fit$anova_table[4,4] <=0.05) {mylist_age_weeks_c[[i]] <- fit$anova_table[4,4,drop=FALSE]}
  if(fit$anova_table[5,4] <=0.05) {mylist_childcarecenter.feeding_type[[i]] <- fit$anova_table[5,4,drop=FALSE]}
  if(fit$anova_table[6,4] <=0.05) {mylist_childcarecenter.time[[i]] <- fit$anova_table[6,4,drop=FALSE]}
  if(fit$anova_table[7,4] <=0.05) {mylist_feeding_type.time[[i]] <- fit$anova_table[7,4,drop=FALSE]}
  if(fit$anova_table[8,4] <=0.05) {mylist_childcarecenter.feeding_type.time[[i]] <- fit$anova_table[8,4,drop=FALSE]}

  if(sum(fit$anova_table$`Pr(>F)` <=0.05) >0 ) {
    diag_df <- 
      tibble( 
        sresid = resid(fit$full_model, scale=T), 
        fitted = fitted(fit$full_model),
        observed = as.numeric(unlist(fit$data[c(paste(i))]))
      )
  # fitted versus standardized residuals
  p1 <- ggplot(diag_df, aes(fitted, sresid)) + geom_point() + geom_smooth(method = 'loess') + ggtitle(glue("Fitted vs Residuals {i}"))
  # fitted versus observed
  p2 <- ggplot(diag_df, aes(fitted, observed)) + geom_point() + geom_smooth(method = 'loess') + ggtitle(glue("Fitted vs Observed {i}"))
  # normality of residuals
  p3 <- ggplot(diag_df, aes(sresid)) +  geom_density() + ggtitle(glue("Density Residuals {i}"))
  # normality of residuals 
  p4 <- gg_qq(diag_df$sresid) + ggtitle(glue("QQplot {i}"))
  

  mylist_diagnostic.plots[[i]] <- 
    cowplot::plot_grid(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), ncol = 2)
  
  taxa <- L2_rel[i,dat$Var.2]
  df <- cbind(taxa,dat)
  dat.all <- melt(df, id.vars='time', measure.vars=c('taxa'))
  dat.all2 <- melt(df, id.vars='subject', measure.vars=c("taxa"))
  dat.all3 <- melt(df, id.vars='feeding_type', measure.vars=c("taxa"))
  dat.all4 <- melt(df, id.vars='childcarecenter', measure.vars=c("taxa"))
  
  #get rid of all the NA's
  dat.all <- cbind(dat.all,dat.all2[,1],dat.all3[,1],dat.all4[,1])
  dat.all <- dat.all[rowSums(is.na(dat.all)) < 1, ]
  
  dat.all[ order(dat.all[,4], dat.all[,1],decreasing = TRUE), ]
  
  colnames(dat.all) <- c("time","variable","value","subject","feeding_type","childcarecenter")
  dat.all$time <- as.factor(dat.all$time)
  dat.all$subject <- as.factor(dat.all$subject)
  
  
  #make the plot for "group"
  p <- ggplot(dat.all, aes(time, value, fill=feeding_type)) 
  p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
  p <- p + geom_boxplot(alpha=.2, width=0.4)
  p <- p + geom_text(aes(label=subject,colour=subject),size=2,hjust=-0.1, vjust=0) + ylab("Relative abundance")
  p <- p + geom_path(aes(group=subject, col=subject))
  #p <- p + geom_smooth(aes(group=childcarecenter, fill = childcarecenter, alpha=0.5),colour="black", method="loess")
  p <- p + theme_bw()
  p <- p + facet_grid("feeding_type")
  p <- p + ggtitle(paste(i))
  #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
  p1 <- p + theme(legend.position = "none")
  

  p <- ggplot(dat.all, aes(time, value, fill=childcarecenter)) 
  p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
  p <- p + geom_boxplot(alpha=.2, width=0.4)
  p <- p + geom_text(aes(label=subject,colour=subject),size=2,hjust=-0.1, vjust=0) + ylab("Relative abundance")
  p <- p + geom_path(aes(group=subject, col=subject))
  #p <- p + geom_smooth(aes(group=childcarecenter, fill = childcarecenter, alpha=0.5),colour="black", method="loess")
  p <- p + theme_bw()
  p <- p + facet_grid("childcarecenter")
  p <- p + ggtitle(paste(i))
  #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
  p2 <- p + theme(legend.position = "none")
  
  
  print(plot_grid(p1, p2,labels = c("A", "B"), ncol = 
  2))

  
  mylist[[i]] <- anova(fit)
  }
} 

## ------------------------------------------------------------------------
mylist_diagnostic.plots
summary(all.models$Bifidobacterium$full_model)
all.models$

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
l <- list(mylist_childcarecenter,mylist_feeding_type,mylist_time,mylist_age_weeks_c,mylist_childcarecenter.feeding_type,mylist_childcarecenter.time,mylist_feeding_type.time,mylist_childcarecenter.feeding_type.time)
names(l) <- c("mylist_childcarecenter","mylist_feeding_type","mylist_time","mylist_age_weeks_c","mylist_childcarecenter.feeding_type","mylist_childcarecenter.time","mylist_feeding_type.time","mylist_childcarecenter.feeding_type.time")



keys <- unique(unlist(lapply(l, names)))
zz <- setNames(do.call(mapply, c(FUN=c, lapply(l, `[`, keys))), keys)

list_df <- function(list) {
  df <- as.data.frame(list)
  names(df) <- list_names(list)
  return (df)
}

list_names <- function(list) {
  
  recursor <- function(list, names) {
    if (is.list(list)) {
      new_names <- paste(names, names(list), sep = ".")
      out <- unlist(mapply(list, new_names, FUN = recursor))
    } else {
      out <- names
    }
    return(out)
  }
  
  new_names <- unlist(mapply(list, names(list), FUN = recursor))
  return(new_names)
}

df <- as.data.frame(t(list_df(zz)))

df$taxa <- rownames(df)


library(dplyr)
library(tidyr)

df <- df %>%
  separate(taxa, c("taxa", "variable"), ".mylist_")
df$variable <- gsub(pattern = '.{7}$',replacement = '',x = df$variable,)
colnames(df)[1] <- "p.value"
print(df)

## ---- eval=FALSE, fig.height=8, fig.width=10, message=FALSE, warning=FALSE, include=FALSE----
## 
## myd <- L2_rel[names(mylist),dat$Var.2]
## 
## all.graphs <- list()
## for (i in rownames(myd)) {
##   taxa <- myd[i,dat$Var.2]
##   df <- cbind(taxa,dat)
## 
##   dat.all <- melt(df, id.vars='time', measure.vars=c('taxa'))
##   dat.all2 <- melt(df, id.vars='subject', measure.vars=c("taxa"))
##   dat.all3 <- melt(df, id.vars='feeding_type', measure.vars=c("taxa"))
##   dat.all4 <- melt(df, id.vars='childcarecenter', measure.vars=c("taxa"))
## 
##   #get rid of all the NA's
##   dat.all <- cbind(dat.all,dat.all2[,1],dat.all3[,1],dat.all4[,1])
##   dat.all <- dat.all[rowSums(is.na(dat.all)) < 1, ]
## 
##   dat.all[ order(dat.all[,4], dat.all[,1],decreasing = TRUE), ]
## 
##   colnames(dat.all) <- c("time","variable","value","subject","feeding_type","childcarecenter")
##   dat.all$time <- as.factor(dat.all$time)
##   dat.all$subject <- as.factor(dat.all$subject)
## 
## 
##   #make the plot for "group"
##   p <- ggplot(dat.all, aes(time, value, fill=feeding_type))
##   p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
##   p <- p + geom_boxplot(alpha=.2, width=0.4)
##   p <- p + geom_text(aes(label=subject,colour=subject),size=2,hjust=-0.1, vjust=0) + ylab("Relative abundance")
##   p <- p + geom_path(aes(group=subject, col=subject))
##   #p <- p + geom_smooth(aes(group=childcarecenter, fill = childcarecenter, alpha=0.5),colour="black", method="loess")
##   p <- p + theme_bw()
##   p <- p + facet_grid(feeding_type~childcarecenter)
##   p <- p + ggtitle(paste(i))
##   #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
##   p <- p + theme(legend.position = "none")
##   p1 <- p
## 
##   all.graphs[[i]] <- p1
## 
## }
## 
## all.graphs

## ---- eval=FALSE, message=FALSE, warning=FALSE, include=FALSE------------
## 
## library(reshape2)
## #myd <- L2_rel
## 
## all.graphs <- list()
## for (i in names(mylist)) {
## #make the plot for "group"
## p <- ggplot(dat.all, aes(real_age, value, fill=feeding_type))
## p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
## #p <- p + geom_boxplot(alpha=.2, width=0.4)
## p <- p + geom_text(aes(label=subject,colour=feeding_type),size=2,hjust=-0.1, vjust=0) + ylab("Abundance (log10)")
## p <- p + geom_path(aes(group=subject, col=feeding_type))
## #p <- p + geom_smooth(aes(group=feeding_type, colour = feeding_type), method="loess")
## p <- p + theme_bw()
## p <- p + facet_grid("feeding_type")
## p <- p + ggtitle(paste(i))
## #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
## p <- p + theme(legend.position = "none")
## p1 <- p
## 
## all.graphs[[i]] <- p1
## 
## }
## 
## 
## plot_grid(plotlist=all.graphs[1:4], nrow=2)
## plot_grid(plotlist=all.graphs[5:8], nrow=2)
## plot_grid(plotlist=all.graphs[9:12], nrow=2)
## plot_grid(plotlist=all.graphs[13:15], nrow=2)
## plot_grid(plotlist=all.graphs[16:19], nrow=2)
## plot_grid(plotlist=all.graphs[20:23], nrow=2)

## ---- eval=FALSE, message=FALSE, warning=FALSE, include=FALSE------------
## library(reshape2)
## #myd <- L2_rel
## 
## all.graphs <- list()
## for (i in names(mylist)) {
## #make the plot for "group"
## p <- ggplot(dat.all, aes(real_age, value, fill=feeding_type))
## p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
## #p <- p + geom_boxplot(alpha=.2, width=0.4)
## p <- p + geom_text(aes(label=subject,colour=feeding_type),size=2,hjust=-0.1, vjust=0) + ylab("Abundance (log10)")
## p <- p + geom_path(aes(group=subject, col=feeding_type))
## #p <- p + geom_smooth(aes(group=feeding_type, colour = feeding_type), method="loess")
## p <- p + theme_bw()
## p <- p + facet_grid("feeding_type")
## p <- p + ggtitle(paste(i))
## #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
## p <- p + theme(legend.position = "none")
## p1 <- p
## 
## all.graphs[[i]] <- p1
## 
## }
## 
## 
## plot_grid(plotlist=all.graphs[1:4], nrow=2)
## plot_grid(plotlist=all.graphs[5:8], nrow=2)
## plot_grid(plotlist=all.graphs[9:12], nrow=2)
## plot_grid(plotlist=all.graphs[13:15], nrow=2)
## plot_grid(plotlist=all.graphs[16:19], nrow=2)
## plot_grid(plotlist=all.graphs[20:23], nrow=2)

## ---- fig.height=8, fig.width=10, message=FALSE, warning=FALSE-----------
library(microbiome)

xx <- dat[,(nrow(L2)+1):ncol(dat)]
rownames(xx) <-(xx$Var.2)

data(dietswap)
data(dietswap)
dietswap@sam_data <- sample_data(xx)
dietswap@otu_table <- otu_table(L2_rel_phyloseq_all[,dat$Var.2], taxa_are_rows = TRUE)
dietswap@tax_table <- tax_table(as.matrix(tax_table))

ps1_dat <- dietswap
############################################
theme_set(theme_bw())

set.seed(666)

ordu.bray_dat = ordinate(ps1_dat, "RDA", "bray",formula = ~ childcarecenter * feeding_type * time + Condition(time) )

z <- plot_ordination(ps1_dat, ordu.bray_dat)
# Oh no, the table wasn't ordered
library("data.table")
newtab = data.table(z$data)
setorder(newtab,subject,time)
z$data <- newtab
p <- ggplot(data = z$data,aes(x = RDA1,y = RDA2))
#p <- p + scale_fill_brewer(type="qual", palette="Set3")
#p <- p + scale_colour_brewer(type="qual", palette="Set3")
p <- p + ggtitle("NMDS bray, only abundant taxa")
p <- p + geom_point(size = 2,alpha=0.3)
#p <- p + stat_ellipse(type = "norm", linetype = 2, level = 0.95)
#p <- p + geom_text(aes(label=subject,colour=feeding_type),size=3,hjust=-0.5, vjust=0)
#p <- p + geom_segment(aes(x=x_coord_start,y=y_coord_start, col=breastfeeding, xend = x_coord_end, yend = y_coord_end),size =0.5, arrow=arrow(length=unit(0.1,"cm")))
#p <- p + geom_text(aes(label=groupcode,colour=groupcode),size=3,hjust=-0.5, vjust=0)
p <- p + geom_path(aes(group=subject, col=feeding_type),alpha=0.3,size =1,arrow=arrow(length=unit(0.35,"cm")),lineend = "butt",linemitre = 10)
#p2 <- p + facet_grid("feeding_type")
print(p)


z <- plot_ordination(ps1_dat, ordu.bray_dat)
# Oh no, the table wasn't ordered
library("data.table")
newtab = data.table(z$data)
setorder(newtab,subject,time)
z$data <- newtab
p <- ggplot(data = z$data,aes(x = RDA1,y = RDA2))
p <- p + stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200) +
  scale_fill_continuous(low = "white", high = "red")
#p <- p + scale_fill_brewer(type="qual", palette="Set3")
#p <- p + scale_colour_brewer(type="qual", palette="Set3")
p <- p + ggtitle("RDA bray , ~ childcarecenter * feeding_type * time + Condition(time) ")
p <- p + geom_point(size = 2,alpha=0.3)
#p <- p + stat_ellipse(type = "norm", linetype = 2, level = 0.95)
#p <- p + geom_text(aes(label=subject,colour=feeding_type),size=3,hjust=-0.5, vjust=0)
#p <- p + geom_segment(aes(x=x_coord_start,y=y_coord_start, col=breastfeeding, xend = x_coord_end, yend = y_coord_end),size =0.5, arrow=arrow(length=unit(0.1,"cm")))
#p <- p + geom_text(aes(label=groupcode,colour=groupcode),size=3,hjust=-0.5, vjust=0)
p <- p + geom_path(aes(group=subject, col=feeding_type),alpha=0.3,size =1,arrow=arrow(length=unit(0.35,"cm")),lineend = "butt",linemitre = 10)
#p2 <- p + facet_grid("feeding_type")
print(p)

p <- ggplot(data = z$data,aes(x = RDA1,y = RDA2))
#p <- p + scale_fill_brewer(type="qual", palette="Set3")
#p <- p + scale_colour_brewer(type="qual", palette="Set3")
p <- p + ggtitle("RDA bray , ~ childcarecenter * feeding_type * time + Condition(time),  density plot")
p <- p + stat_density2d(aes(fill = ..density..^0.5), geom = "tile", contour = FALSE, n = 200) +
  scale_fill_continuous(low = "white", high = "red")
p <- p + geom_point(size = 2,alpha=0.3)
#p <- p + stat_ellipse(type = "norm", linetype = 2, level = 0.95)
#p <- p + geom_text(aes(label=subject,colour=feeding_type),size=3,hjust=-0.5, vjust=0)
#p <- p + geom_segment(aes(x=x_coord_start,y=y_coord_start, col=breastfeeding, xend = x_coord_end, yend = y_coord_end),size =0.5, arrow=arrow(length=unit(0.1,"cm")))
#p <- p + geom_text(aes(label=groupcode,colour=groupcode),size=3,hjust=-0.5, vjust=0)
p <- p + geom_path(aes(group=subject, col=feeding_type),alpha=0.3,size =1,arrow=arrow(length=unit(0.35,"cm")),lineend = "butt",linemitre = 10)
p <- p + facet_grid("feeding_type")
print(p)


#with arrows
p <- plot_ordination(ps1_dat, ordu.bray_dat,type = "taxa",label="Genus", color="Order")
# Oh no, the table wasn't ordered
#p <- p + scale_fill_brewer(type="qual", palette="Set3")
#p <- p + scale_colour_brewer(type="qual", palette="Set3")
p3 <- p + ggtitle("RDA bray , ~ childcarecenter * feeding_type * time + Condition(time) ")
#p <- p + geom_point(size = 2)
#p <- p + stat_ellipse(type = "norm", linetype = 2, level = 0.95)
#p <- p + geom_text(aes(label=subject,colour=feeding_type),size=3,hjust=-0.5, vjust=0)
#p <- p + geom_segment(aes(x=x_coord_start,y=y_coord_start, col=breastfeeding, xend = x_coord_end, yend = y_coord_end),size =0.5, arrow=arrow(length=unit(0.1,"cm")))
#p <- p + geom_text(aes(label=groupcode,colour=groupcode),size=3,hjust=-0.5, vjust=0)
#p <- p + geom_path(aes(group=subject, col=feeding_type),size =1,arrow=arrow(length=unit(0.35,"cm")),lineend = "butt",linemitre = 10)
#p <- p +facet_grid("feeding_type")
print(p3)

## ---- echo=FALSE, fig.height=12, fig.width=12, message=FALSE, warning=FALSE----
metadata <- all.data
metadata <- metadata[metadata$antibio != 1,]
#metadata <- metadata[metadata$time <76 & metadata$time > 3 | metadata$time == 105 ,]

bp.data <- metadata[,c("real_age","subject","feeding_type","childcarecenter")]
bp.data1 <- bp.data[complete.cases(bp.data),]


library(reshape2)
myd <- L1_rel

all.graphs <- list()
for (i in rownames(myd)) {
  taxa <- myd[i,rownames(bp.data)]
  df <- cbind(taxa,bp.data)

  dat.all <- melt(df, id.vars='real_age', measure.vars=c('taxa'))
  dat.all2 <- melt(df, id.vars='subject', measure.vars=c("taxa"))
  dat.all3 <- melt(df, id.vars='feeding_type', measure.vars=c("taxa"))
  dat.all4 <- melt(df, id.vars='childcarecenter', measure.vars=c("taxa"))

  #get rid of all the NA's
  dat.all <- cbind(dat.all,dat.all2[,1],dat.all3[,1],dat.all4[,1])
  dat.all <- dat.all[rowSums(is.na(dat.all)) < 1, ]

  dat.all[ order(dat.all[,4], dat.all[,1],decreasing = TRUE), ]

  colnames(dat.all) <- c("real_age","variable","value","subject","feeding_type","childcarecenter")
  dat.all$real_age <- as.factor(dat.all$real_age)

  #make the plot for "group"
  p <- ggplot(dat.all, aes(real_age, value, fill=feeding_type)) 
  p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
  #p <- p + geom_boxplot(alpha=.2, width=0.4)
  p <- p + geom_text(aes(label=subject,colour=feeding_type),size=2,hjust=-0.1, vjust=0) + ylab("Abundance (log10)")
  p <- p + geom_path(aes(group=subject, col=feeding_type))
  #p <- p + geom_smooth(aes(group=feeding_type, colour = feeding_type), method="loess")
  p <- p + theme_bw()
  p <- p + facet_grid("feeding_type")
  p <- p + ggtitle(paste(i))
  #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
  p <- p + theme(legend.position = "none")
  p1 <- p                                                       

  all.graphs[[i]] <- p1

}

#plot_grid(plotlist=all.graphs[1:(length(all.graphs)/3)], nrow=3)
#plot_grid(plotlist=all.graphs[(length(all.graphs)/3 +1):length(all.graphs)], nrow=3)

plot_grid(plotlist=all.graphs[1:9], nrow=3)
plot_grid(plotlist=all.graphs[10:18], nrow=3)
plot_grid(plotlist=all.graphs[19:22], nrow=2)

## ---- echo=FALSE, fig.height=8, fig.width=12, message=FALSE, warning=FALSE----
metadata <- all.data
metadata <- metadata[metadata$antibio != 1,]
#metadata <- metadata[metadata$time <76 & metadata$time > 3 | metadata$time == 105 ,]

bp.data <- metadata[,c("real_age","subject","feeding_type","childcarecenter")]
bp.data1 <- bp.data[complete.cases(bp.data),]


library(reshape2)
myd <- L2_rel[abundant_taxa,]

all.graphs <- list()
for (i in rownames(myd)) {
  taxa <- myd[i,rownames(bp.data)]
  df <- cbind(taxa,bp.data)

  dat.all <- melt(df, id.vars='real_age', measure.vars=c('taxa'))
  dat.all2 <- melt(df, id.vars='subject', measure.vars=c("taxa"))
  dat.all3 <- melt(df, id.vars='feeding_type', measure.vars=c("taxa"))
  dat.all4 <- melt(df, id.vars='childcarecenter', measure.vars=c("taxa"))

  #get rid of all the NA's
  dat.all <- cbind(dat.all,dat.all2[,1],dat.all3[,1],dat.all4[,1])
  dat.all <- dat.all[rowSums(is.na(dat.all)) < 1, ]

  dat.all[ order(dat.all[,4], dat.all[,1],decreasing = TRUE), ]

  colnames(dat.all) <- c("real_age","variable","value","subject","feeding_type","childcarecenter")
  dat.all$real_age <- as.factor(dat.all$real_age)

  #make the plot for "group"
  p <- ggplot(dat.all, aes(real_age, value, fill=feeding_type)) 
  p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
  #p <- p + geom_boxplot(alpha=.2, width=0.4)
  p <- p + geom_text(aes(label=subject,colour=feeding_type),size=2,hjust=-0.1, vjust=0) + ylab("Abundance (log10)")
  p <- p + geom_path(aes(group=subject, col=feeding_type))
  #p <- p + geom_smooth(aes(group=feeding_type, colour = feeding_type), method="loess")
  p <- p + theme_bw()
  p <- p + facet_grid("feeding_type")
  p <- p + ggtitle(paste(i))
  #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
  p <- p + theme(legend.position = "none")
  p1 <- p                                                       

  all.graphs[[i]] <- p1

}


plot_grid(plotlist=all.graphs[1:9], nrow=3)
plot_grid(plotlist=all.graphs[10:13], nrow=2)

plot_grid(plotlist=all.graphs[c("Bifidobacterium","Streptococcus_bovis_et_rel","Streptococcus_intermedius_et_rel","Streptococcus_mitis_et_rel","Lactobacillus_plantarum_et_rel","Enterococcus","Escherichia_coli_et_rel")], nrow=3)

plot_grid(plotlist=all.graphs[setdiff(names(all.graphs),c("Bifidobacterium","Streptococcus_bovis_et_rel","Streptococcus_intermedius_et_rel","Streptococcus_mitis_et_rel","Lactobacillus_plantarum_et_rel","Enterococcus","Escherichia_coli_et_rel"))], nrow=3)


#plot_grid(plotlist=all.graphs[19:27], nrow=3)
#plot_grid(plotlist=all.graphs[28:36], nrow=3)
#plot_grid(plotlist=all.graphs[37:45], nrow=3)
#plot_grid(plotlist=all.graphs[46:55], nrow=3)

## ---- echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
used_samples <- as.character(rownames(subset(my.metadata, !(sampleID %in% duplicates))))

metadata <- my.metadata[used_samples,]
metadata <- metadata[metadata$antibio != 1,]
metadata <- metadata[metadata$time <76 & metadata$time > 3 | metadata$time == 105 ,]


bp.data1 <- metadata[,c("time","subject","feeding_type","childcarecenter")]
bp.data2 <- bp.data1[complete.cases(bp.data1),]


myd <- L1_rel[,rownames(bp.data2)]

all.graphs <- list()
for (i in rownames(myd)) {
  taxa <- myd[i,rownames(bp.data2)]
  df <- cbind(taxa,bp.data2)
  
  dat.all <- melt(df, id.vars='time', measure.vars=c('taxa'))
  dat.all2 <- melt(df, id.vars='subject', measure.vars=c("taxa"))
  dat.all3 <- melt(df, id.vars='feeding_type', measure.vars=c("taxa"))
  dat.all4 <- melt(df, id.vars='childcarecenter', measure.vars=c("taxa"))
  
  #get rid of all the NA's
  dat.all <- cbind(dat.all,dat.all2[,1],dat.all3[,1],dat.all4[,1])
  dat.all <- dat.all[rowSums(is.na(dat.all)) < 1, ]
  
  dat.all[ order(dat.all[,4], dat.all[,1],decreasing = TRUE), ]
  
  colnames(dat.all) <- c("time","variable","value","subject","feeding_type","childcarecenter")
  dat.all$time <- as.factor(dat.all$time)
  
  #make the plot for "group"
  p <- ggplot(dat.all, aes(time, value, fill=feeding_type)) 
  p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
  p <- p + geom_boxplot(alpha=.2, width=0.4)
  p <- p + geom_text(aes(label=subject,colour=subject),size=2,hjust=-0.1, vjust=0) + ylab("Abundance (log10)")
  p <- p + geom_path(aes(group=subject, col=subject))
  #p <- p + geom_smooth(aes(group=breastfeeding, colour = breastfeeding), method="loess")
  p <- p + theme_bw()
  p <- p + facet_grid("childcarecenter")
  p <- p + ggtitle(paste(i))
  #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
  p <- p + theme(legend.position = "none")
  p1 <- p                                                       
  
  all.graphs[[i]] <- p1
  
}

all.graphs$Actinobacteria
all.graphs$Bacilli
all.graphs$Bacteroidetes
all.graphs$Proteobacteria

## ---- echo=FALSE, fig.height=10, fig.width=10, message=FALSE, warning=FALSE----
used_samples <- as.character(rownames(subset(my.metadata, !(sampleID %in% duplicates))))

metadata <- my.metadata[used_samples,]
metadata <- metadata[metadata$antibio != 1,]
metadata <- metadata[metadata$time <76 & metadata$time > 3 | metadata$time == 105 ,]


bp.data1 <- metadata[,c("time","subject","feeding_type","childcarecenter")]
bp.data2 <- bp.data1[complete.cases(bp.data1),]


myd <- L2_rel[,rownames(bp.data2)]

all.graphs <- list()
for (i in rownames(myd)) {
  taxa <- myd[i,rownames(bp.data2)]
  df <- cbind(taxa,bp.data2)
  
  dat.all <- melt(df, id.vars='time', measure.vars=c('taxa'))
  dat.all2 <- melt(df, id.vars='subject', measure.vars=c("taxa"))
  dat.all3 <- melt(df, id.vars='feeding_type', measure.vars=c("taxa"))
  dat.all4 <- melt(df, id.vars='childcarecenter', measure.vars=c("taxa"))
  
  #get rid of all the NA's
  dat.all <- cbind(dat.all,dat.all2[,1],dat.all3[,1],dat.all4[,1])
  dat.all <- dat.all[rowSums(is.na(dat.all)) < 1, ]
  
  dat.all[ order(dat.all[,4], dat.all[,1],decreasing = TRUE), ]
  
  colnames(dat.all) <- c("time","variable","value","subject","feeding_type","childcarecenter")
  dat.all$time <- as.factor(dat.all$time)
  
  #make the plot for "group"
  p <- ggplot(dat.all, aes(time, value, fill=feeding_type)) 
  p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
  p <- p + geom_boxplot(alpha=.2, width=0.4)
  p <- p + geom_text(aes(label=subject,colour=subject),size=2,hjust=-0.1, vjust=0) + ylab("Abundance (log10)")
  p <- p + geom_path(aes(group=subject, col=subject))
  #p <- p + geom_smooth(aes(group=breastfeeding, colour = breastfeeding), method="loess")
  p <- p + theme_bw()
  p <- p + facet_grid("childcarecenter")
  p <- p + ggtitle(paste(i))
  #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
  p <- p + theme(legend.position = "none")
  p1 <- p                                                       
  
  all.graphs[[i]] <- p1
  
}

all.graphs$Bifidobacterium
all.graphs$Streptococcus_bovis_et_rel
all.graphs$Streptococcus_mitis_et_rel
all.graphs$Enterococcus
all.graphs$Staphylococcus
all.graphs$Enterococcus

## ---- echo=FALSE, fig.height=8, fig.width=8, message=FALSE, warning=FALSE----
used_samples <- as.character(rownames(subset(my.metadata, !(sampleID %in% duplicates))))

metadata <- my.metadata[used_samples,]
metadata <- metadata[metadata$antibio != 1,]
metadata <- metadata[metadata$time <76 & metadata$time > 3 | metadata$time == 105 ,]


bp.data1 <- metadata[,c("time","subject","feeding_type","childcarecenter")]
bp.data2 <- bp.data1[complete.cases(bp.data1),]


myd <- L2_rel[,rownames(bp.data2)]

all.graphs <- list()
for (i in rownames(myd)) {
  taxa <- myd[i,rownames(bp.data2)]
  df <- cbind(taxa,bp.data2)
  
  dat.all <- melt(df, id.vars='time', measure.vars=c('taxa'))
  dat.all2 <- melt(df, id.vars='subject', measure.vars=c("taxa"))
  dat.all3 <- melt(df, id.vars='feeding_type', measure.vars=c("taxa"))
  dat.all4 <- melt(df, id.vars='childcarecenter', measure.vars=c("taxa"))
  
  #get rid of all the NA's
  dat.all <- cbind(dat.all,dat.all2[,1],dat.all3[,1],dat.all4[,1])
  dat.all <- dat.all[rowSums(is.na(dat.all)) < 1, ]
  
  dat.all[ order(dat.all[,4], dat.all[,1],decreasing = TRUE), ]
  
  colnames(dat.all) <- c("time","variable","value","subject","feeding_type","childcarecenter")
  dat.all$time <- as.factor(dat.all$time)
  
  #make the plot for "group"
  p <- ggplot(dat.all, aes(time, value, fill=feeding_type)) 
  p <- p + geom_point(position = position_dodge(width=1),alpha=0.5)
  p <- p + geom_boxplot(alpha=.2, width=0.4)
  p <- p + geom_text(aes(label=subject,colour=subject),size=2,hjust=-0.1, vjust=0) + ylab("Abundance (log10)")
  p <- p + geom_path(aes(group=subject, col=subject))
  p <- p + geom_smooth(aes(group=childcarecenter, fill = childcarecenter, alpha=0.5),colour="black", method="loess")
  p <- p + theme_bw()
  p <- p + facet_grid(feeding_type~childcarecenter)
  p <- p + ggtitle(paste(i))
  #p <- p + scale_x_discrete(limits=c("0","2","7","14","28","75","79","84","91","105","2193"))
  p <- p + theme(legend.position = "none")
  p1 <- p                                                       
  
  all.graphs[[i]] <- p1
  
}

all.graphs$Bifidobacterium
all.graphs$Streptococcus_bovis_et_rel
all.graphs$Streptococcus_mitis_et_rel
all.graphs$Enterococcus
all.graphs$Staphylococcus
all.graphs$Enterococcus

## ------------------------------------------------------------------------
# # uncomment to store as r script
# knitr::purl("analyses_complete.Rmd", here("R/analyses_complete.R"))

all.models$Bifidobacterium

model_input %>% group_by(feeding_type) %>% summarise(n = n())
model_input[1, 1:112] %>% sum()
L2_rel
