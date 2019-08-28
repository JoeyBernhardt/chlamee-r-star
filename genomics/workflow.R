
## Download the vcf files from Dryad; name their folder "vcfs"

library(tidyverse)
library(fs)
library(beyonce)
library(car)
library(VariantAnnotation)

########################################################################
## Define convenience functions for reading the variant calling files ##
########################################################################

get_variant_table <- function(vcf_object)
{
    genome_matrix <- geno(vcf_object)$GT
    read_depth <- geno(vcf_object)$DP
    var_df <- data.frame(ID=rownames(genome_matrix), Variant=genome_matrix,
                        Depth=read_depth,
                        row.names=NULL) %>%
        separate(ID, into=c("Chromosome", "Loc"), sep=":") %>%
        separate(Loc, into=c("Locus", "Variants"), sep="_") # %>%
    names(var_df)[4] <- "Variant_type"
    names(var_df)[5] <- "Depth"
    var_df
}

get_no_alternative_table <- function(vcf_object)
{
    alternate <- unlist(geno(vcf_object)$AO[,1])
    var_df <- data.frame(ID=names(alternate), Alternate=alternate,
                        row.names=NULL) %>%
        separate(ID, into=c("Chromosome", "Loc"), sep=":") %>%
        separate(Loc, into=c("Locus", "Variants"), sep="_") # %>%
    names(var_df)[4] <- "Alternate"
    var_df
}

get_vcf_table <- function(vcf_object)
{
    var_df1 <- get_variant_table(vcf_object)
    var_df2 <- get_no_alternative_table(vcf_object)
    merge(var_df1, var_df2, by=c("Chromosome", "Locus", "Variants")) %>%
        mutate(Proportion=Alternate/Depth*100)
}

##################################################
## Prepare the mutation table.                  ##
## Takes a while to run; EXECUTE WITH CAUTION!! ##
##################################################

dir_ls("vcfs") %>%
    keep(~grepl("filt", .)) %>% 
    keep(~grepl("vcf", .)) %>% 
    map_dfr(~get_vcf_table(readVcf(.)),
            .id = "File") %>% 
    separate(File, c("junk1", "Name"), sep="filt_") %>%
    separate(Name, c("junk2", "Ancestor", "Treatment", "junk3"), sep="_") %>% 
    dplyr::select(-junk1, -junk2, -junk3) %>%
    unite(Strain, Ancestor, Treatment, sep="_", remove=FALSE) %>% 
    write_csv("long_variant_table.csv")

############################################
## Computationally intense part ends here ##
############################################


##################################################################################
## Prepare the comparisons of SNP variability between Ancestors and Descendants ##
##################################################################################

long_variant_table <-
    read_csv("long_variant_table.csv") %>% 
    mutate(Variant_type = ifelse(Proportion > 70, 1, 0)) %>% 
    filter(Depth > 20,
           Depth < 250)

mut_table <- 
    long_variant_table %>% 
    unite(Loc, Chromosome, Locus, sep="_") %>%
    dplyr::select(Loc, Variant_type, Strain) %>%
    spread(Strain, Variant_type) %>%
    mutate(rsum = rowSums(.[-1], na.rm = TRUE)) %>%
    filter(rsum != 37) %>%
    dplyr::select(-rsum) %>% 
    gather(Str, Val, -Loc) %>%
    separate(Str, c("Strain", "Treatment"), "_")


########################################################################################
## Prepare pairwise comparisons of each ancestor to all of its respective descendants ##
########################################################################################

mut_count <- function(strain)
{
    exp_subset <- 
        filter(mut_table,
               Strain == strain,
               !is.na(Val))
    left_join(
        filter(exp_subset, Treatment == "A"),
        filter(exp_subset, Treatment != "A"),
        by="Loc") %>%
        filter(Val.x != Val.y) %>%
        group_by(Treatment.y) %>%
        summarise(n=n()) %>%
        mutate(Strain = strain) %>%
        dplyr::rename(Treatment = Treatment.y)
}

cbbPalette <- c("#6b6b6b",
                "#f9a729",
                "#97cfd0",
                "#00a2b3",
                "#f1788d",
                "#cf3e53",
                "#b9ca5d")

map_dfr(
    c("CC1691", "Anc2", "Anc3", "Anc4", "Anc5"),
    mut_count) %>%
    spread(Treatment, n)

## 111100715 bases in Chlamy genome
## 285 generations
map_dfr(
    c("CC1691", "Anc2", "Anc3", "Anc4", "Anc5"),
    mut_count) %>%
    mutate(mu = n / 285 / 111100715) %>% 
    ggplot(aes(x=Treatment, y=mu)) +
    geom_boxplot() +
    geom_point(aes(color=Strain)) +
    scale_color_manual(values = beyonce_palette(type = "discrete", 18), name = "")
ggsave("treatment_effect_on_mu.pdf", last_plot(), useDingbats=FALSE)

map_dfr(
    c("CC1691", "Anc2", "Anc3", "Anc4", "Anc5"),
    mut_count) %>%
    mutate(mu = n / 285 / 111100715) %>% 
    ggplot(aes(x=Strain, y=mu)) +
    geom_boxplot() +
    geom_point(aes(color=Treatment)) +
    scale_color_manual(values = cbbPalette)
ggsave("ancestor_effect_on_mu.pdf", last_plot(), useDingbats=FALSE)

mut_data <- 
    map_dfr(
        c("CC1691", "Anc2", "Anc3", "Anc4", "Anc5"),
        mut_count) %>%
    mutate(mu = n / 285 / 111100715)

## Prepare the ANOVAs
tr_aov <- aov(mut_data$mu ~ mut_data$Treatment)
str_aov <- aov(mut_data$mu ~ mut_data$Strain)

## Checking the variance homogeneity
leveneTest(mu ~ Treatment, data=mut_data)
leveneTest(mu ~ Strain, data=mut_data)
## OK!

## Checking the normality of the factor levels
tr_aov %>% residuals %>% shapiro.test
str_aov %>% residuals %>% shapiro.test
## OK!

## Not significant; p = 0.788
tr_aov %>% summary

## Highly significant; p = 5.3e-8
str_aov %>% summary

################################################################
## Prepare pairwise comparisons of ancestors Anc2-5 to CC1690 ##
################################################################

comp_count <- function(strain1, strain2)
{
    cc <- filter(mut_table, Strain == strain1, Treatment == "A")
    anc <- filter(mut_table, Strain == strain2, Treatment == "A")
    left_join(cc, anc, by="Loc") %>%
        filter(complete.cases(.),
               Val.x != Val.y) %>%
        dim %>%
        .[1]
}

snp_table <- 
    c("CC1690", "Anc2", "Anc3", "Anc4", "Anc5") %>% 
    combn(2) %>%
    t %>%
    data.frame %>%
    mutate(Diff = map2(X1, X2, ~comp_count(.x, .y)))

################################################
## Check the genomic organization of the SNPs ##
################################################

complete_mutations <-
    filter(long_variant_table,
           Depth > 20,
           Depth < 250,
           Treatment != "A") %>% 
    mutate(Mutation = ifelse(Proportion > 70, 1, 0)) %>%
    unite(ID, Chromosome, Locus, sep=":") %>%
    unite(ID, ID, Variants, sep="_") %>%
    dplyr::select(Ancestor, Treatment, ID, Mutation) %>%
    unite(Strain, Ancestor, Treatment, sep="-") %>%
    spread(Strain, Mutation, fill=NA) %>%
    filter(complete.cases(.)) %>%
    gather(Strain, Mutation, -ID) %>%
    group_by(ID) %>%
    mutate(Mut_sums = sum(Mutation)) %>%
    filter(Mut_sums != 32, Mut_sums != 0) %>%
    dplyr::select(-Mut_sums)

complete_mutations %>%
    spread(ID, Mutation) %>%
    dplyr::select(-Strain) %>%
    chisq.test

