make_test_phyloseq <- function() {
  counts <- matrix(c(10,0,5,3,4,8,0,2,1,4,7,6,2,1,3,9,8,5,2,1,5,6,4,2), nrow = 4, byrow = FALSE,
    dimnames = list(paste0("ASV", 1:4), paste0("S", 1:6)))
  tax <- matrix(c("Bacteria","Firmicutes","GenusA","SpeciesA","Bacteria","Firmicutes","GenusB","SpeciesB","Bacteria","Bacteroidota","GenusC","SpeciesC","Bacteria","Actinobacteriota","GenusD","SpeciesD"), nrow = 4, byrow = TRUE,
    dimnames = list(paste0("ASV", 1:4), c("Kingdom", "Phylum", "Genus", "Species")))
  meta <- data.frame(group = factor(rep(c("Control", "Treatment"), each = 3)), subject = factor(paste0("P", 1:6)), age = c(20,22,24,21,23,25), row.names = paste0("S", 1:6))
  phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE), phyloseq::tax_table(tax), phyloseq::sample_data(meta))
}
