## In utero exposure to caffeine and acetaminophen, the gut microbiome, and neurodevelopmental outcomes, a prospective birth cohort study: Cross-sectional analysis notebook
## Author: Yike Shen (ys3419@cumc.columbia.edu)

setwd(“/Users/yikeshen/CaffeineMicrobiome”)

#Housekeeping
``` r
rm(list=ls())
#intall and load packages
require()
require("readxl")
require(plyr)
require(dplyr)
require(devtools)
require(phyloseq)
require(readr)
require(tidyverse)
require(tidyr)
require(Maaslin2)
require(scales)
require(microbiome)
require(ape)
require(vegan)
require(RColorBrewer)
require(lmtest)
R.version.string
```

    ## [1] "R version 4.1.2 (2021-11-01)"

#Metadata cleanup

``` r
#Due to patient privacy, we can not share any input files, please contact Dr. Yike Shen (ys3419@cumc.columbia.edu) to connect with GESTE cohort to obtain input data files.

RAW <- read.csv("Caffine_data_raw.csv")

RAW <- RAW %>% as.data.frame()
Caffineprocessing<- apply(RAW, 2, function(x) gsub("^$|^ $", NA, x)) %>% as.data.frame()

Caffineprocessing$meco_caffeine[Caffineprocessing$meco_caffeine == "."] <- NA

Caffineprocessing <- Caffineprocessing %>% as.data.frame() %>% 
  dplyr::rename(SampleID=record_id) %>%
  dplyr::rename(Caffeine=meco_caffeine) %>% 
  dplyr::rename(Acetaminophen=meco_acetaminophene) %>%
  dplyr::rename(Maternal_recruitment_age=momagerecrutmcalcul) %>% 
  dplyr::rename(Maternal_BMI_prepreg=bmipregros_cal) %>% 
  dplyr::rename(Family_Income_Y6=revenufami) %>%
  dplyr::rename(Parity=primipar) %>% 
  dplyr::rename(Gestational_age=ga_wks_bb) %>% 
  dplyr::rename(Delivery_mode=acc_cs) %>% 
  dplyr::rename(Sex=sexe_bb) %>% 
  dplyr::rename(Child_birthweight=masse_bb) %>% 
  dplyr::rename(Breastfed=maternelbf) %>% 
  dplyr::rename(Block_Design=bls_pond) %>% 
  dplyr::rename(Coding=code_pond) %>% 
  dplyr::rename(Digit_span=sct_pond) %>% 
  dplyr::rename(Information=connaiss_p) %>% 
  dplyr::rename(Vocabulary=vocab_pond) %>% 
  dplyr::rename(QTAC=qtac) 
  

Caffineprocessing <- Caffineprocessing %>% 
  dplyr::select(SampleID,Caffeine,Acetaminophen,Maternal_recruitment_age,Maternal_BMI_prepreg,Family_Income_Y6,Parity,Gestational_age,Delivery_mode,Delivery_mode,Sex,Child_birthweight,Breastfed,Block_Design,Coding,Digit_span,Information,Vocabulary,QTAC)

Caffineprocessing$Parity[Caffineprocessing$Parity == " 1"] <- 'Yes'
Caffineprocessing$Parity[Caffineprocessing$Parity == " 0"] <- 'No'

Caffineprocessing$Delivery_mode[Caffineprocessing$Delivery_mode == " 1"] <- 'Vaginal'
Caffineprocessing$Delivery_mode[Caffineprocessing$Delivery_mode == " 2"] <- 'C_section'

Caffineprocessing$Sex[Caffineprocessing$Sex == " 0"] <- 'Female'
Caffineprocessing$Sex[Caffineprocessing$Sex == " 1"] <- 'Male'

Caffineprocessing$Breastfed[Caffineprocessing$Breastfed == " 1"] <- 'Ever'
Caffineprocessing$Breastfed[Caffineprocessing$Breastfed == " 0"] <- 'Never'
Caffineprocessing$Breastfed[Caffineprocessing$Breastfed == "99"] <- NA

Caffineprocessing$Maternal_BMI_prepreg[Caffineprocessing$Maternal_BMI_prepreg == " 0.00000"] <- NA
Caffineprocessing$Maternal_BMI_prepreg[Caffineprocessing$Maternal_BMI_prepreg == "     Inf"] <- NA


Caffineprocessing$Acetaminophen[Caffineprocessing$Acetaminophen == "<LOD"] <-'ND'
Caffineprocessing$Acetaminophen[Caffineprocessing$Acetaminophen == "LOD<"] <-'ND'
Caffineprocessing$Acetaminophen[Caffineprocessing$Acetaminophen == "<LOQ"] <-'ND'

#Participant X,Y,Z NA imputation with population mean
#removed in md due to patient information protection

SampleIDFULL <- Caffineprocessing$SampleID
rownames(Caffineprocessing) <- SampleIDFULL

Caffineprocessing[,c(2,4:6,8,11,13:18)] <- lapply(Caffineprocessing[,c(2,4:6,8,11,13:18)], function(x) as.numeric(as.character(x)))

WISCSUM <- Caffineprocessing[,13:17]
WISCSUM <- rowSums(WISCSUM) %>% as.data.frame()
colnames(WISCSUM) <- "WISC_sum"

Caffineprocessing <- cbind(Caffineprocessing,WISCSUM)
Caffineprocessing <- Caffineprocessing[,-1]

#Table1:Summary Statistics-All
#table1::table1(~., data = Caffineprocessing)

Caffineprocessing_compeletecace <- Caffineprocessing[complete.cases(Caffineprocessing[ ,2]),]
Caffineprocessing_compeletecaffeine <- Caffineprocessing[complete.cases(Caffineprocessing[ ,1]),]

#Table1-intersect:Caffeine-microbiome

PathwayClean <- read.csv("Pathway_clean.csv")
PathwayClean <- PathwayClean %>% dplyr::rename(SampleID=X)

Y6Con_raw <- read.csv("Y6_df_cleaned.csv")
Y6Con_raw <- Y6Con_raw[,-1]
Y6Con_process <- Y6Con_raw


ID_mecocaffiene <- rownames(Caffineprocessing_compeletecaffeine) %>% as.data.frame()
colnames(ID_mecocaffiene) <- 'SampleID'

#Rawdata has empty space for all character string, SampleID name fix
ID_mecocaffiene$SampleID <- as.numeric(ID_mecocaffiene$SampleID)
ID_mecocaffiene$SampleID <- as.character(ID_mecocaffiene$SampleID)
Caffineprocessing_compeletecaffeine <- cbind(ID_mecocaffiene,Caffineprocessing_compeletecaffeine)


PathwayClean$SampleID <- as.character(PathwayClean$SampleID)
Y6Con_process$SampleID <- as.character(Y6Con_process$SampleID)
microbiome_85 <- match_df(Y6Con_process,PathwayClean,on ="SampleID")

ID_365 <- rownames(Caffineprocessing)
ID_365 <- ID_365 %>% as.numeric() %>% as.character() %>% as.data.frame()
colnames(ID_365) <- "SampleID"
Caffineprocessing <- cbind(ID_365,Caffineprocessing)

CaffeineMatch <- match_df(microbiome_85,Caffineprocessing, on ="SampleID")
CaffeineMatch_var <- match_df(Caffineprocessing, microbiome_85, on ="SampleID")

CaffeineMatch <- cbind(CaffeineMatch_var,CaffeineMatch)
CaffeineMatch <- CaffeineMatch[,-c(1:3)]
CaffeineMatch <- CaffeineMatch[,-17]

#table1::table1(~., data = CaffeineMatch)


CaffeineMatchModel <- CaffeineMatch %>% select(CaffeineY6,Breastfed,Delivery_mode,Family_Income_Y6,Sex)
AceMatchMatchModel <- CaffeineMatch %>% select(AcetaminophenY6,Breastfed,Delivery_mode,Family_Income_Y6,Sex)
```

#Metaphlan table cleanup

``` r
SPECIES <- lapply(list.files(pattern = "*merged_abundance_table.tsv"), function(i){read_tsv(i)}) 

IDS <- lapply(SPECIES,function(x) substring(names(x),1,regexpr("_dehost_metaphlan_bugs_list",names(x))-1))
IDS <- unlist(IDS)[-1]

SPECIES <- lapply(SPECIES, function(x) 
  separate(data = x, col = "#clade_name", into = c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep = "\\|"))

SPECIES <- lapply(SPECIES, subset)
SPECIES <- SPECIES %>% reduce(full_join, by = "Domain")
colnames(SPECIES)[8:ncol(SPECIES)] <- IDS
SPECIES <- SPECIES[,-8] %>% as.data.frame()

SPECIES_Species <- SPECIES[complete.cases(SPECIES[ ,7]),]
SPECIES_Species_relab <- SPECIES_Species
dim(SPECIES_Species[duplicated(SPECIES_Species$Species),])[1]

Species_phyla <- SPECIES[!rowSums(is.na.data.frame(SPECIES))<5,]
Species_phyla <- Species_phyla[complete.cases(Species_phyla[ ,2]),]
Species_phyla <- Species_phyla[,-c(3:7)]
```

#Caffeine and acetaminophen concentration normalization

``` r
#log transfer caffeine
#turn acetaminophen ND as refernce group

CaffeineMatchModel$CaffeineY6 <- log2(CaffeineMatchModel$CaffeineY6)
AceMatchMatchModel$AcetaminophenY6 <- as.factor(AceMatchMatchModel$AcetaminophenY6)
AceMatchMatchModel <- within(AceMatchMatchModel, AcetaminophenY6 <- relevel(AcetaminophenY6, ref = "ND"))
```

#Turn to phyloseq object-caffeine

``` r
#Note: Name it OTU, but it's not OTU! It's species relative abundance
OTU <- SPECIES_Species_relab
ROWNAMES_Species <- OTU$Species
OTU <- OTU[,-c(1:7)]

row.names(OTU) <- NULL
row.names(OTU) <- ROWNAMES_Species
OTU <- OTU %>% as.matrix()

colnames(OTU) <- paste("P",colnames(OTU),sep="_")

Taxonomy <- SPECIES_Species_relab %>% 
  dplyr::select("Domain","Phylum","Class","Order","Family","Genus","Species")

TAX <- Taxonomy %>% as.data.frame() %>% 
  mutate(Domain=gsub("d__","",Domain)) %>% 
  mutate(Phylum=gsub("p__","",Phylum)) %>% 
  mutate(Class=gsub("c__","",Class)) %>% 
  mutate(Order=gsub("o__","",Order)) %>% 
  mutate(Family=gsub("f__","",Family)) %>% 
  mutate(Genus=gsub("g__","",Genus)) %>% 
  mutate(Species=gsub("s__","",Species)) %>% 
  as.matrix() %>% as.data.frame()

TAX <- TAX %>% as.matrix()
row.names(TAX) <- ROWNAMES_Species

map <- CaffeineMatchModel
map <- cbind(rownames(CaffeineMatchModel),map) 
map$`rownames(CaffeineMatchModel)` <- as.numeric(map$`rownames(CaffeineMatchModel)`)
List_49 <- map$`rownames(CaffeineMatchModel)`
row.names(map) <- List_49
row.names(map) <- paste("P",row.names(map),sep="_")
map <- map[,-1]

OTU=otu_table(OTU,taxa_are_rows = TRUE)
rownames(TAX) <- row.names(OTU)
TAX=tax_table(TAX)
map <- map %>% as.data.frame()
MAP <- sample_data(map)

physeq_caffeine = merge_phyloseq(OTU, TAX, MAP)
Tree=rtree(ntaxa(physeq_caffeine),rooted=TRUE,tip.label = taxa_names(physeq_caffeine))

physeq_caffeine = phyloseq(OTU, TAX, MAP,Tree)
```

#Turn to phyloseq object-acetaminophen

``` r
#Note: Name it OTU, but it's not OTU! It's species relative abundance
OTU <- SPECIES_Species_relab
ROWNAMES_Species <- OTU$Species
OTU <- OTU[,-c(1:7)]

row.names(OTU) <- NULL
row.names(OTU) <- ROWNAMES_Species
OTU <- OTU %>% as.matrix()

colnames(OTU) <- paste("P",colnames(OTU),sep="_")

Taxonomy <- SPECIES_Species_relab %>% 
  dplyr::select("Domain","Phylum","Class","Order","Family","Genus","Species")

TAX <- Taxonomy %>% as.data.frame() %>% 
  mutate(Domain=gsub("d__","",Domain)) %>% 
  mutate(Phylum=gsub("p__","",Phylum)) %>% 
  mutate(Class=gsub("c__","",Class)) %>% 
  mutate(Order=gsub("o__","",Order)) %>% 
  mutate(Family=gsub("f__","",Family)) %>% 
  mutate(Genus=gsub("g__","",Genus)) %>% 
  mutate(Species=gsub("s__","",Species)) %>% 
  as.matrix() %>% as.data.frame()

TAX <- TAX %>% as.matrix()
row.names(TAX) <- ROWNAMES_Species

map_ace <- AceMatchMatchModel
map_ace <- cbind(rownames(AceMatchMatchModel),map_ace) 
map_ace$`rownames(AceMatchMatchModel)` <- as.numeric(map_ace$`rownames(AceMatchMatchModel)`)
List_46 <- map_ace$`rownames(AceMatchMatchModel)`
row.names(map_ace) <- List_46
row.names(map_ace) <- paste("P",row.names(map_ace),sep="_")
map_ace <- map_ace[,-1]

OTU=otu_table(OTU,taxa_are_rows = TRUE)
rownames(TAX) <- row.names(OTU)
TAX=tax_table(TAX)
map_ace <- map_ace %>% as.data.frame()
MAP <- sample_data(map_ace)

physeq_ace = merge_phyloseq(OTU, TAX, MAP)
Tree=rtree(ntaxa(physeq_ace),rooted=TRUE,tip.label = taxa_names(physeq_ace))

physeq_ace = phyloseq(OTU, TAX, MAP,Tree)
```

#Alpha diversity

``` r
AlphaEvenness <- evenness(physeq_caffeine, 'pielou')
AlphaRichness <- alpha(physeq_caffeine,'shannon')
ALPHA_variables <- cbind(map,AlphaRichness,AlphaEvenness)

#alpha diversity distribution check
ALPHA_variables$diversity_shannon <- scale(ALPHA_variables$diversity_shannon)
Shannondensity <- ggplot(ALPHA_variables, aes(x=diversity_shannon)) + 
    geom_histogram(aes(y=..density..), sbinwidth=.5,
                   colour="black", fill="white") +
      geom_density(alpha=.2, fill="#FF6666")+ggtitle("Alpha Diversity Density Plot - Shannon")+
  theme(legend.position="top",axis.text=element_text(size=16),
        axis.title=element_text(size=16),plot.title = element_text(size=22, hjust=0.5),text = element_text(size = 16),axis.text.x = element_text(angle = 0))+labs(x="Shannon Diversity", y="Count")

ALPHA_variables$pielou <- scale(ALPHA_variables$pielou)
Pieloudensity <- ggplot(ALPHA_variables, aes(x=pielou)) + 
    geom_histogram(aes(y=..density..), sbinwidth=.5,
                   colour="black", fill="white") +
      geom_density(alpha=.2, fill="#FF6666")+ggtitle("Alpha Diversity Density Plot - Pielou")+
  theme(legend.position="top",axis.text=element_text(size=16),
        axis.title=element_text(size=16),plot.title = element_text(size=22, hjust=0.5),text = element_text(size = 16),axis.text.x = element_text(angle = 0))+labs(x="Pielou Evenness", y="Count")

#alpha diversity association calculation
summary(shanCaff <- lm(diversity_shannon ~ CaffeineY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex, data=ALPHA_variables))
```

    ## 
    ## Call:
    ## lm(formula = diversity_shannon ~ CaffeineY6 + Breastfed + Delivery_mode + 
    ##     Family_Income_Y6 + Sex, data = ALPHA_variables)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8394 -0.5381  0.2155  0.6292  1.6812 
    ## 
    ## Coefficients:
    ##                        Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)           7.948e-01  4.147e-01   1.916   0.0589 .
    ## CaffeineY6           -5.660e-02  4.501e-02  -1.258   0.2123  
    ## BreastfedNever        2.725e-01  2.935e-01   0.929   0.3560  
    ## Delivery_modeVaginal -2.572e-01  2.806e-01  -0.916   0.3623  
    ## Family_Income_Y6     -2.296e-06  2.417e-06  -0.950   0.3449  
    ## SexMale              -3.614e-01  2.193e-01  -1.648   0.1034  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9868 on 79 degrees of freedom
    ## Multiple R-squared:  0.08421,    Adjusted R-squared:  0.02625 
    ## F-statistic: 1.453 on 5 and 79 DF,  p-value: 0.2148

``` r
par(mfrow=c(2,2))
plot(shanCaff,which=1:4,main="regression diagnose_Caffeine~Shannon")

summary(piCaff <- lm(pielou ~ CaffeineY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex, data=ALPHA_variables))
```

    ## 
    ## Call:
    ## lm(formula = pielou ~ CaffeineY6 + Breastfed + Delivery_mode + 
    ##     Family_Income_Y6 + Sex, data = ALPHA_variables)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.94649 -0.45790  0.06161  0.71364  1.87023 
    ## 
    ## Coefficients:
    ##                        Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)           6.290e-01  4.144e-01   1.518    0.133
    ## CaffeineY6           -4.868e-02  4.497e-02  -1.082    0.282
    ## BreastfedNever        4.193e-01  2.932e-01   1.430    0.157
    ## Delivery_modeVaginal -1.511e-01  2.804e-01  -0.539    0.591
    ## Family_Income_Y6     -2.207e-06  2.415e-06  -0.914    0.363
    ## SexMale              -3.306e-01  2.191e-01  -1.509    0.135
    ## 
    ## Residual standard error: 0.9859 on 79 degrees of freedom
    ## Multiple R-squared:  0.08586,    Adjusted R-squared:  0.02801 
    ## F-statistic: 1.484 on 5 and 79 DF,  p-value: 0.2045

``` r
par(mfrow=c(2,2))
plot(piCaff,which=1:4,main="regression diagnose_Caffeine~Pielou")

AlphaEvenness_ace <- evenness(physeq_ace, 'pielou')
AlphaRichness_ace <- alpha(physeq_ace,'shannon')

ALPHA_variables_ace <- cbind(map_ace,AlphaRichness_ace,AlphaEvenness_ace)

summary(shanAce <- lm(diversity_shannon ~ AcetaminophenY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex, data=ALPHA_variables_ace))

par(mfrow=c(2,2))
plot(shanCaff,which=1:4,main="regression diagnose_Acetaminophen~Shannon")

summary(piAce <- lm(pielou ~ AcetaminophenY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex, data=ALPHA_variables_ace))

par(mfrow=c(2,2))
plot(piCaff,which=1:4,main="regression diagnose_Acetaminophen~Pielou")

AlphaFig <- data.frame(Method = factor(),
                       Exposure = factor(),
                       beta = numeric(),
                       se = numeric(),
                       ci_low=numeric(),
                       ci_high=numeric(),
                       pval=numeric(),
                       FDR_qval=numeric(),
                       stringsAsFactors=FALSE)

AlphaFig[1:4,] <- NA
AlphaFig$Method <- rep(c("Shannon","Pielou"))
AlphaFig$Method <- factor(AlphaFig$Method, levels = c("Shannon", "Pielou"))

AlphaFig$Exposure <- c("Caffeine","Caffeine","Acetaminophen","Acetaminophen")

AlphaFig$beta[1] <- shanCaff$coefficients["CaffeineY6"]
AlphaFig$se[1] <- coef(summary(shanCaff))["CaffeineY6", "Std. Error"] 
AlphaFig$pval[1] <- coef(summary(shanCaff))["CaffeineY6", "Pr(>|t|)"]
AlphaFig$beta[2] <- piCaff$coefficients[shanCaff$coefficients["CaffeineY6"]]
AlphaFig$se[2] <- coef(summary(piCaff))["CaffeineY6", "Std. Error"]
AlphaFig$pval[2] <- coef(summary(piCaff))["CaffeineY6", "Pr(>|t|)"]
AlphaFig$beta[3] <- shanAce$coefficients["AcetaminophenY6D"]
AlphaFig$se[3] <- coef(summary(shanAce))["AcetaminophenY6D", "Std. Error"]
AlphaFig$pval[3] <- coef(summary(shanAce))["AcetaminophenY6D", "Pr(>|t|)"]
AlphaFig$beta[4] <- piAce$coefficients["AcetaminophenY6D"]
AlphaFig$se[4] <- coef(summary(piAce))["AcetaminophenY6D", "Std. Error"]
AlphaFig$pval[4] <- coef(summary(piAce))["AcetaminophenY6D", "Pr(>|t|)"]

AlphaFig$ci_low <- AlphaFig$beta-1.96*AlphaFig$se
AlphaFig$ci_high <- AlphaFig$beta+1.96*AlphaFig$se

AlphaFig$FDR_qval <- p.adjust(AlphaFig$pval, method = "fdr")

AlphaFigure <- ggplot(AlphaFig, aes(Method, beta), group=Exposure) +
                geom_point(size = 3,  position = position_dodge(width = 1)) +
                geom_errorbar(data = AlphaFig, aes(ymin = beta - 1.96*se,
                                  ymax = beta + 1.96*se, linetype = NULL),
                              width = 0.3, size =0.6,
                              position=position_dodge(width=1),
                              group = AlphaFig$Exposure) +
                geom_hline(yintercept = 0, size = 1) +
                scale_x_discrete(name = "Alpha Diversity Indices", labels = c("Shannon","Pielou")) +
                scale_y_continuous(name = "Adjusted Differences in Alpha Diversity \nAmong Exposed Compared to Unexposed") +
                theme_bw() +
                theme(panel.grid.major =element_blank(), 
                      panel.grid.minor =element_blank(),
                      axis.text = element_text(size = 18),
                      axis.title = element_text(size = 18)) +
                facet_wrap( ~ Exposure,scales = "free") +
                theme(strip.background = element_rect(fill = "white", color = "black", size = 1),strip.text = element_text(size = 18))+
  theme(panel.spacing = unit(2, "cm"))
```

#Beta diversity

``` r
adonisTABLE <- data.frame(Exposure = factor(),
                      Method = character(),
                      Adjust_covariates = numeric(),
                      FDR_qval=numeric(),
                      Adjust_covariates_r2 = numeric())
adonisTABLE[1:8,] <- NA
adonisTABLE$Exposure <- factor(adonisTABLE$Exposure, levels = c("CaffeineY6", "AcetaminophenY6"))
adonisTABLE$Exposure[1:4] <- "CaffeineY6"
adonisTABLE$Exposure[5:8] <- "AcetaminophenY6"

adonisTABLE$Method <- rep(c("Weighted UniFrac", "Bray-Curtis", "Unweighted UniFrac","Jaccard"))
adonisTABLE$Method <- factor(adonisTABLE$Method, levels = c("Weighted UniFrac","Bray-Curtis","Unweighted UniFrac", "Jaccard"))

#Weighted Unifraq
yy_weighted_caff <- UniFrac(physeq_caffeine, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
yy_weighted_ace <- UniFrac(physeq_ace, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)

set.seed(4)
weighted_caff <- adonis(formula=yy_weighted_caff~CaffeineY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex,data=CaffeineMatchModel, permutations = 999, contr.unordered = "contr.sum")
adonisTABLE[1,3] <- weighted_caff$aov.tab["CaffeineY6","Pr(>F)"]
adonisTABLE$Adjust_covariates_r2[1] <- weighted_caff$aov.tab["CaffeineY6","R2"]

set.seed(4)
weighted_ace <- adonis(formula=yy_weighted_ace~AcetaminophenY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex,data=AceMatchMatchModel, permutations = 999, contr.unordered = "contr.sum")
adonisTABLE[5,3] <- weighted_ace$aov.tab["AcetaminophenY6","Pr(>F)"]
adonisTABLE$Adjust_covariates_r2[5] <- weighted_ace$aov.tab["AcetaminophenY6","R2"]

#Unweighted unifraq
yy_unweighted_caff <- UniFrac(physeq_caffeine, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
yy_unweighted_ace <- UniFrac(physeq_ace, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)

set.seed(4)
unweighted_caff <- adonis(formula=yy_unweighted_caff~CaffeineY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex, data=CaffeineMatchModel, permutations = 999, contr.unordered = "contr.sum")
adonisTABLE[3,3] <- unweighted_caff$aov.tab["CaffeineY6","Pr(>F)"]
adonisTABLE$Adjust_covariates_r2[3] <- unweighted_caff$aov.tab["CaffeineY6","R2"]

set.seed(4)
unweighted_ace <- adonis(formula=yy_unweighted_ace~AcetaminophenY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex, data=AceMatchMatchModel, permutations = 999, contr.unordered = "contr.sum")
adonisTABLE[7,3] <- unweighted_ace$aov.tab["AcetaminophenY6","Pr(>F)"]
adonisTABLE$Adjust_covariates_r2[7] <- unweighted_ace$aov.tab["AcetaminophenY6","R2"]

#Bray-Curtis Distance
YY_bray_caff <- distance(physeq_caffeine, "bray")
YY_bray_ace <- distance(physeq_ace, "bray")

set.seed(4)
bray_caff <- adonis(formula=YY_bray_caff~CaffeineY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex,data=CaffeineMatchModel, permutations = 999, contr.unordered = "contr.sum")
adonisTABLE[2,3] <- bray_caff$aov.tab["CaffeineY6","Pr(>F)"]
adonisTABLE$Adjust_covariates_r2[2] <- bray_caff$aov.tab["CaffeineY6","R2"]

set.seed(4)
bray_ace <- adonis(formula=YY_bray_ace~AcetaminophenY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex,data=AceMatchMatchModel, permutations = 999, contr.unordered = "contr.sum")
adonisTABLE[6,3] <- bray_ace$aov.tab["AcetaminophenY6","Pr(>F)"]
adonisTABLE$Adjust_covariates_r2[6] <- bray_ace$aov.tab["AcetaminophenY6","R2"]

#jaccard distance
yy_jaccard_caff <- distance(physeq_caffeine,"jaccard",binary = TRUE)
yy_jaccard_ace <- distance(physeq_ace,"jaccard",binary = TRUE)

set.seed(4)
jaccard_caff <- adonis(formula=yy_jaccard_caff~CaffeineY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex, data=CaffeineMatchModel, permutations = 999, contr.unordered = "contr.sum")
adonisTABLE[4,3] <- jaccard_caff$aov.tab["CaffeineY6","Pr(>F)"]
adonisTABLE$Adjust_covariates_r2[4] <- jaccard_caff$aov.tab["CaffeineY6","R2"]

set.seed(4)
jaccard_ace <- adonis(formula=yy_jaccard_ace~AcetaminophenY6+Breastfed+Delivery_mode+Family_Income_Y6+Sex, data=AceMatchMatchModel, permutations = 999, contr.unordered = "contr.sum")
adonisTABLE[8,3] <- jaccard_ace$aov.tab["AcetaminophenY6","Pr(>F)"]
adonisTABLE$Adjust_covariates_r2[8] <- jaccard_ace$aov.tab["AcetaminophenY6","R2"]


adonisTABLE$FDR_qval <- p.adjust(adonisTABLE$Adjust_covariates, method = "fdr")
adonisTABLE[,3:4] <- adonisTABLE[,3:4] %>% round(3)
#write.csv(BetasigTable,"betaccassociations.csv")
```

#Species association

``` r
SPECIES_Species <- SPECIES_Species[,-c(1:6)]
rownames(SPECIES_Species) <- NULL
SPECIESLIST <- SPECIES_Species$Species
rownames(SPECIES_Species) <- SPECIESLIST
SPECIES_Species <- SPECIES_Species[,-1]
SPECIES_Species_2 <- t(SPECIES_Species)
Samplelist_85 <- rownames(SPECIES_Species_2) %>% as.data.frame()
colnames(Samplelist_85) <- 'SampleID'
SPECIES_Species_2 <- cbind(Samplelist_85,SPECIES_Species_2)
rownames(SPECIES_Species_2) <- NULL

#Caffeine species association
CaffeineMatch_species <- match_df(SPECIES_Species_2,Samplelist_85, on ="SampleID")
CaffeineMatch_species <- CaffeineMatch_species[,colSums(CaffeineMatch_species == 0) <= nrow(CaffeineMatch_species)*0.9]

input_metadata_species <- CaffeineMatchModel %>% as.data.frame()
input_data_species <- CaffeineMatch_species %>% as.data.frame()
rownames(input_data_species) <- rownames(input_metadata_species)
input_data_species <- input_data_species[,-1]

fit_data3 = Maaslin2(
    input_data = input_data_species,
    input_metadata = input_metadata_species,
    output = "masslin_species_CaffeineY6",
    fixed_effects = c("CaffeineY6","Breastfed","Delivery_mode","Family_Income_Y6","Sex"))

FDRq_Caffeine_species <- read.table(file="./masslin_species_CaffeineY6/all_results.tsv",sep = '\t', header = TRUE)
FDRq_Caffeine_species <- FDRq_Caffeine_species[which(FDRq_Caffeine_species$metadata=='CaffeineY6', ),]
correctedq_Caffeine_species <- p.adjust(FDRq_Caffeine_species$pval, method = "fdr")
FDRq_Caffeine_species <- cbind(FDRq_Caffeine_species,correctedq_Caffeine_species)

#acetaminophen species association
AceMatch_species <- match_df(SPECIES_Species_2,Samplelist_85, on ="SampleID")
AceMatch_species <- AceMatch_species[,colSums(AceMatch_species == 0) <= nrow(AceMatch_species)*0.9]

input_metadata_species_ace <- AceMatchMatchModel %>% as.data.frame()
input_data_species_ace <- AceMatch_species %>% as.data.frame()
rownames(input_data_species_ace) <- rownames(input_metadata_species_ace)
input_data_species_ace <- input_data_species_ace[,-1]

fit_data4 = Maaslin2(
    input_data = input_data_species_ace,
    input_metadata = input_metadata_species_ace,
    output = "masslin_species_AceY6",
    fixed_effects = c("AcetaminophenY6","Breastfed","Delivery_mode","Family_Income_Y6","Sex"))

FDRq_Ace_species <- read.table(file="./masslin_species_AceY6/all_results.tsv",sep = '\t', header = TRUE)
FDRq_Ace_species <- FDRq_Ace_species[which(FDRq_Ace_species$metadata=='AcetaminophenY6', ),]
correctedq_Ace_species <- p.adjust(FDRq_Ace_species$pval, method = "fdr")
FDRq_Ace_species <- cbind(FDRq_Ace_species,correctedq_Ace_species)
```

#Phylum association

``` r
Species_phyla <- Species_phyla[,-1]
rownames(Species_phyla) <- NULL
SPECIES_phylaLIST <- Species_phyla$Phylum
rownames(Species_phyla) <- SPECIES_phylaLIST
Species_phyla <- Species_phyla[,-1]
SPECIES_phyla_2 <- t(Species_phyla)
Samplelist_85 <- rownames(SPECIES_phyla_2) %>% as.data.frame()
colnames(Samplelist_85) <- 'SampleID'
SPECIES_phyla_2 <- cbind(Samplelist_85,SPECIES_phyla_2)
rownames(SPECIES_phyla_2) <- NULL

######Caffeine phylum association
CaffeineMatch_phy <- SPECIES_phyla_2
CaffeineMatch_phy <- CaffeineMatch_phy[,colSums(CaffeineMatch_phy == 0) <= nrow(CaffeineMatch_phy)*0.9]

input_metadata_phy <- CaffeineMatchModel %>% as.data.frame()
input_data_phy <- CaffeineMatch_phy %>% as.data.frame()
rownames(input_data_phy) <- rownames(input_metadata_phy)
input_data_phy <- input_data_phy[,-1]

fit_data5 = Maaslin2(
    input_data = input_data_phy,
    input_metadata = input_metadata_phy,
    output = "masslin_phy_CaffeineY6",
    fixed_effects = c("CaffeineY6","Breastfed","Delivery_mode","Family_Income_Y6","Sex"))

FDRq_Caffeine_phy <- read.table(file="./masslin_phy_CaffeineY6/all_results.tsv",sep = '\t', header = TRUE)
FDRq_Caffeine_phy <- FDRq_Caffeine_phy[which(FDRq_Caffeine_phy$metadata=='CaffeineY6', ),]
correctedq_Caffeine_phy <- p.adjust(FDRq_Caffeine_phy$pval, method = "fdr")
FDRq_Caffeine_phy <- cbind(FDRq_Caffeine_phy,correctedq_Caffeine_phy)
#no significant

#acetaminophen phylum association
AceMatch_phy <- SPECIES_phyla_2
AceMatch_phy <- AceMatch_phy[,colSums(AceMatch_phy == 0) <= nrow(AceMatch_phy)*0.9]

input_metadata_phy_ace <- AceMatchMatchModel %>% as.data.frame()
input_data_phy_ace <- AceMatch_phy %>% as.data.frame()
rownames(input_data_phy_ace) <- rownames(input_metadata_phy_ace)
input_data_phy_ace <- input_data_phy_ace[,-1]

fit_data6 = Maaslin2(
    input_data = input_data_phy_ace,
    input_metadata = input_metadata_phy_ace,
    output = "masslin_phy_AceY6",
    fixed_effects = c("AcetaminophenY6","Breastfed","Delivery_mode","Family_Income_Y6","Sex"))

FDRq_Ace_phy <- read.table(file="./masslin_phy_AceY6/all_results.tsv",sep = '\t', header = TRUE)
FDRq_Ace_phy <- FDRq_Ace_phy[which(FDRq_Ace_phy$metadata=='AcetaminophenY6', ),]
correctedq_Ace_phy <- p.adjust(FDRq_Ace_phy$pval, method = "fdr")
FDRq_Ace_phy <- cbind(FDRq_Ace_phy,correctedq_Ace_phy)
```

#Pathway association

``` r
#caffeine pathway
CaffeineMatch_pathway <- PathwayClean

CaffeineMatch_pathway <- CaffeineMatch_pathway[,colSums(CaffeineMatch_pathway == 0) <= nrow(CaffeineMatch_pathway)*0.9]

input_metadata_path <- CaffeineMatchModel %>% as.data.frame()
input_data_path <- CaffeineMatch_pathway %>% as.data.frame()
rownames(input_data_path) <- rownames(input_metadata_path)
input_data_path <- input_data_path[,-1]

fit_data1 = Maaslin2(
    input_data = input_data_path,
    input_metadata = input_metadata_path,
    output = "masslin_pathway_CaffeineY6",
    fixed_effects = c("CaffeineY6","Breastfed","Delivery_mode","Family_Income_Y6","Sex"))

FDRq_Caffeine <- read.table(file="./masslin_pathway_CaffeineY6/all_results.tsv",sep = '\t', header = TRUE)
FDRq_Caffeine <- FDRq_Caffeine[which(FDRq_Caffeine$metadata=='CaffeineY6', ),]
correctedq_Caffeine <- p.adjust(FDRq_Caffeine$pval, method = "fdr")
FDRq_Caffeine <- cbind(FDRq_Caffeine,correctedq_Caffeine)
#3 sig
#PWY.5345 :superpathway of L-methionine biosynthesis (by sulfhydrylation)
#SO4ASSIM.PWY :sulfate reduction I (assimilatory)
#SULFATE.CYS.PWY :superpathway of sulfate assimilation and cysteine biosynthesis 
#grep "PWY.5345" humann_pathway_relab.tsv|column -t -s $'\t' |less -S


#acetaminophen pathway
AceMatch_pathway <- PathwayClean
AceMatch_pathway <- AceMatch_pathway[,colSums(AceMatch_pathway == 0) <= nrow(AceMatch_pathway)*0.9]

input_metadata_path_ace <- AceMatchMatchModel %>% as.data.frame()
input_data_path_ace <- AceMatch_pathway %>% as.data.frame()
rownames(input_data_path_ace) <- rownames(input_metadata_path_ace)
input_data_path_ace <- input_data_path_ace[,-1]

fit_data2 = Maaslin2(
    input_data = input_data_path_ace,
    input_metadata = input_metadata_path_ace,
    output = "masslin_pathway_AceY6",
    fixed_effects = c("AcetaminophenY6","Breastfed","Delivery_mode","Family_Income_Y6","Sex"))

FDRq_Ace <- read.table(file="./masslin_pathway_AceY6/all_results.tsv",sep = '\t', header = TRUE)
FDRq_Ace <- FDRq_Ace[which(FDRq_Ace$metadata=='AcetaminophenY6', ),]
correctedq_Ace <- p.adjust(FDRq_Ace$pval, method = "fdr")
FDRq_Ace <- cbind(FDRq_Ace,correctedq_Ace)

#write.csv(FDRq_Ace,"FDRq_Ace_pathY6.csv")
#write.csv(FDRq_Caffeine,"FDRq_Caff_pathY6.csv")

Pathwaymean <- CaffeineMatch_pathway[,-1]
Pathwaymean <- colMeans(Pathwaymean) %>% as.matrix()
Pathwaymean <- Pathwaymean %>% as.data.frame()
colnames(Pathwaymean) <- 'AverageRelativeAbundance'
```

#pathway volcano plot

``` r
Valcanoraw <- FDRq_Caffeine
volcanoraw_ace <- FDRq_Ace

Valcanoraw['q_value'] <- NA
Valcanoraw$q_value[Valcanoraw$correctedq_Caffeine< 0.1] <- "q value<0.1"
Valcanoraw$q_value[Valcanoraw$correctedq_Caffeine> 0.1] <- "q value>0.1"
Valcanoraw['Label'] <- NA

Valcanoraw$Label[Valcanoraw$pval> 0.0001 & Valcanoraw$pval<0.000108] <- "L-methionine biosynthesis"
Valcanoraw$Label[Valcanoraw$pval> 0.00014 & Valcanoraw$pval<0.00015] <- "Sulfate reduction I (assimilatory)"
Valcanoraw$Label[Valcanoraw$pval> 0.00011 & Valcanoraw$pval<0.00012] <-"Sulfate assimilation and cysteine biosynthesis"

Valcanoraw_RA <- Valcanoraw
rownames(Valcanoraw_RA) <- Valcanoraw_RA$feature
Valcanoraw_RA_ <- merge(Valcanoraw_RA, Pathwaymean,by = 'row.names', all = TRUE)
Valcanoraw_RA_ <- Valcanoraw_RA_[-c(360,361),]

library("ggrepel")
mycolors <- c("blue","grey")

Caffpath_plot <- ggplot(data=Valcanoraw_RA_, aes(x=coef, y=-log10(pval), col=q_value, label = Label)) + 
  theme_bw()+geom_point()+ 
  geom_hline(yintercept=-log10(0.1), col="blue",linetype='longdash')+
  scale_x_continuous(name = "Difference in Pathway Relative Abundance \nper Doubling of Caffeine (ng/g)")+
  scale_y_continuous(name = "-log10(p-value)")+ 
  geom_point(aes(size = AverageRelativeAbundance), color=8)+
  scale_size_continuous(range = c(1.5, 10))+
  #geom_text_repel(size = 5,hjust=0.95,nudge_y = 0.1)+
  geom_text_repel(size = 5)+
  #ggrepel::geom_text_repel(aes(size = 1,hjust=0.95,nudge_y = 3.8,label = Label))+
  scale_color_manual(breaks = c("q value<0.1","q value>0.1"),
                     values=c("blue","grey"))+
  ggtitle("Caffeine Pathway Potential")+theme(legend.position="none")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16, hjust=0.5),
        text = element_text(size = 16))
Caffpath_plot

volcanoraw_ace_RA <- volcanoraw_ace
rownames(volcanoraw_ace_RA) <- volcanoraw_ace$feature
volcanoraw_ace_RA <- merge(volcanoraw_ace_RA, Pathwaymean,by = 'row.names', all = TRUE)

mycolors <- c("blue","grey")
volcanoraw_ace$q_value[volcanoraw_ace$correctedq_Ace< 0.1] <- "q value<0.1"
volcanoraw_ace['Label'] <- NA
volcanoraw_ace_RA <- volcanoraw_ace_RA[-c(360,361),]

Acepath_plot <- ggplot(data=volcanoraw_ace_RA, aes(x=coef, y=-log10(pval))) + 
  theme_bw()+geom_point()+ 
  geom_hline(yintercept=-log10(0.1), col="blue",linetype='longdash')+
  scale_x_continuous(name = "Difference in Pathway Relative Abundance \nComparing Exposed to Unexposed")+scale_y_continuous(name = "-log10(p-value)")+ 
  geom_point(aes(size = AverageRelativeAbundance), color=8)+
  scale_size_continuous(range = c(1.5, 10))+
  scale_color_manual(breaks = c("q value<0.1","q value>0.1"),
                     values=c("blue","grey"))+
  ggtitle("Acetaminophen Pathway Potential")+theme(legend.position="none")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        plot.title = element_text(size=16, hjust=0.5),
        text = element_text(size = 16))
Acepath_plot
```
