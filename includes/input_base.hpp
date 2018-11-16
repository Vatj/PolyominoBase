#include "genotype_metrics.hpp"

// The purpose of this extension is to create an interactive app which will take a genome as input and will spit out the information

namespace simulation_params
{
  extern uint16_t n_genes, colours, metric_colours;
  extern uint32_t n_jiggle;
}

namespace io_params
{
  extern std::string set_metric_file, genome_metric_file, neighbour_file;
}


std::pair<Genotype_Metrics, Genome_to_Set> single_genome_to_metric(Genotype genome, PhenotypeTable* pt);

void multiple_genomes_to_metric(std::vector<Genotype> genomes, PhenotypeTable* pt);

void genome_to_pID_distribution(Genotype genome, PhenotypeTable* pt);
