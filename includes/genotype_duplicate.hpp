#include "genotype_metrics.hpp"

std::vector<Genotype> GenomesDuplication(std::vector<Genotype> genomes);
std::vector<Genotype> GeneDuplication(Genotype& genotype);
void GenomesJiggleDuplication(std::vector<Genotype>& genomes, std::vector<Genotype>& jiggle_genomes, std::vector<Genotype>& duplicates);
