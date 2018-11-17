#include "input_base.hpp"


std::pair<Genotype_Metrics, Genome_to_Set> single_genome_to_metric(Genotype genome, PhenotypeTable* pt)
{
  Genotype genotype = genome;

  std::map<Phenotype_ID, uint8_t> pID_counter = GetPIDCounter(genotype, pt);

  std::vector<Phenotype_ID> ref_pIDs;
  for(auto pID: pID_counter)
    ref_pIDs.emplace_back(pID.first);

  double neutral_weight = ((double) NeutralSize(genotype, 1, simulation_params::metric_colours - 1));

  Genotype_Metrics genome_metric(simulation_params::n_genes, simulation_params::metric_colours);
  genome_metric.pID_counter.insert(std::begin(pID_counter), std::end(pID_counter));
  genome_metric.set_reference(genotype, ref_pIDs, neutral_weight);

  Genome_to_Set neighbourhood;

  for(Genotype neighbour : genotype_neighbourhood(genotype))
  {
     std::vector<Phenotype_ID> neighbour_pIDs = GetSetPIDs(neighbour, pt);
     neighbourhood[neighbour] = neighbour_pIDs;
     genome_metric.analyse_pIDs(neighbour_pIDs);
   }

   // std::pair<Genotype_Metrics, Genome_to_Set> data = std::make_pair(genome_metric, neighbourhood);
   //
   // std::ofstream genome_metric_out(io_params::genome_metric_file, std::ios_base::app);
   // genome_metric.save_to_file(genome_metric_out);
   // genome_metric_out.close();
   //
   // return data;
   return std::make_pair(genome_metric, neighbourhood);
}

void multiple_genomes_to_metric(std::vector<Genotype> genomes, PhenotypeTable* pt)
{
  std::vector<std::pair<Genotype_Metrics, Genome_to_Set>> data_genomes;

  for (auto genome: genomes)
    // single_genome_to_metric(genome, pt);
    data_genomes.emplace_back(single_genome_to_metric(genome, pt));

  //
  // for (std::vector<Genotype>::iterator iter = std::begin(genomes); iter != std::end(genomes); iter++)
  //   data_genomes.emplace_back(single_genome_to_metric(*iter, pt));

  std::cout << "Printing to files : \n";
  std::cout << io_params::genome_metric_file << std::endl;
  std::cout << io_params::neighbour_file << std::endl;

  std::ofstream set_metric_out(io_params::set_metric_file);
  std::ofstream genome_metric_out(io_params::genome_metric_file);
  header_metric_files(set_metric_out, genome_metric_out);
  std::ofstream neighbour_out(io_params::neighbour_file);

  for(std::vector<std::pair<Genotype_Metrics, Genome_to_Set>>::iterator iter = std::begin(data_genomes);
  iter != std::end(data_genomes); iter++)
  {
    (iter->first).save_to_file(genome_metric_out);
    PrintNeighbourhood(neighbour_out, iter->second);
  }
}

void genome_to_pID_distribution(Genotype genome, PhenotypeTable* pt)
{
  double neutral_weight = 0;

  Genotype genotype = genome;
  std::vector <Phenotype_ID> ref_pIDs = GetSetPIDs(genotype, pt);

  std::cout << "Metric distribution of {";
  for(auto pID: ref_pIDs)
    std::cout << "(" <<+ pID.first << ", " <<+ pID.second << "), ";
  std::cout << "} over " <<+ simulation_params::n_jiggle << " jiggles! \n";

  // Create new files and open them for writing (erase previous data)
  std::ofstream set_metric_out(io_params::set_metric_file);
  std::ofstream genome_metric_out(io_params::genome_metric_file);
  header_metric_files(set_metric_out, genome_metric_out);

  Set_Metrics set_metrics(simulation_params::n_genes, simulation_params::metric_colours);
  set_metrics.ref_pIDs = ref_pIDs;

  neutral_weight = ((double) NeutralSize(genotype, 1, simulation_params::metric_colours - 1)) / simulation_params::n_jiggle; // This weight will be counted n_jiggle time when added to form the total neutral weight
  set_metrics.originals.emplace_back(genotype);

  #pragma omp parallel for schedule(dynamic) firstprivate(genotype, neutral_weight)
  for(uint32_t nth_jiggle=0; nth_jiggle<simulation_params::n_jiggle; ++nth_jiggle)
  {
    Clean_Genome(genotype, false);
    JiggleGenotype(genotype);

    Genotype_Metrics genome_metric(simulation_params::n_genes, simulation_params::metric_colours);
    genome_metric.set_reference(genotype, ref_pIDs, neutral_weight);

    std::map<Phenotype_ID, uint8_t> pID_counter = GetPIDCounter(genotype, pt);
    genome_metric.pID_counter.insert(std::begin(pID_counter), std::end(pID_counter));

    std::vector<Phenotype_ID> pIDs;
    for(auto pID: pID_counter)
      pIDs.emplace_back(pID.first);

    if(pIDs != ref_pIDs)
      set_metrics.misclassified[genotype] = pIDs;

    for(Genotype neighbour : genotype_neighbourhood(genotype))
    {
      std::vector<Phenotype_ID> neighbour_pIDs = GetSetPIDs(neighbour, pt);
      genome_metric.analyse_pIDs(neighbour_pIDs);
    }
    #pragma omp critical
    {
      set_metrics.add_genotype_metrics(genome_metric);
    }

    // metrics.emplace_back(set_metrics);
    set_metrics.save_to_file(set_metric_out, genome_metric_out);
  }
}
