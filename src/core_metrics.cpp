#include "core_metrics.hpp"

// Constructor of the Genotype_Metrics structure
Genotype_Metrics::Genotype_Metrics(uint16_t n_genes, uint16_t colours):
n_genes(n_genes), colours(colours)
{
  number_of_neighbours = (colours - 1) * n_genes * 4.;
}

void Genotype_Metrics::set_reference(Genotype& genotype, std::vector<Phenotype_ID> pIDs, double neutral)
{
  ref_genotype = genotype;
  ref_pIDs = pIDs;
  original = genotype;
  Clean_Genome(original, false);
  max_size = 0;

  for(auto pID: ref_pIDs)
    max_size = std::max(max_size, pID.first);

  neutral_weight = neutral;

  // for (auto pID: pIDs)
  //   shapes.emplace_back(Shape_Metrics(pID));
}

void Genotype_Metrics::analyse_pIDs(std::vector <Phenotype_ID>& pIDs)
{
  if (std::find(std::begin(pIDs), std::end(pIDs), unbound_pID) != std::end(pIDs))
  {
    unbound += 1.;
    return;  // Unbound is an exclusive label
  }

  if (pIDs == ref_pIDs)
    strict_robustness += 1.;

  if (std::find(std::begin(pIDs), std::end(pIDs), rare_pID) != std::end(pIDs))
  {
    rare += 1.;
    pIDs.erase(std::find(std::begin(pIDs), std::end(pIDs), rare_pID), std::end(pIDs));
  }

  std::vector <Phenotype_ID> intersection, union_set;

  std::set_intersection(std::begin(pIDs), std::end(pIDs), std::begin(ref_pIDs), std::end(ref_pIDs), std::back_inserter(intersection));
  std::set_union(std::begin(pIDs), std::end(pIDs), std::begin(ref_pIDs), std::end(ref_pIDs), std::back_inserter(union_set));

  if(intersection.size() > 0)
    // intersection_robustness += (double) intersection.size() / (double) ref_pIDs.size();
    intersection_robustness += 1;
  // else
  //   intersection_robustness += 0;

  if(union_set.size() > ref_pIDs.size())
  {
    union_evolvability += 1;
    if(intersection.size() > 0)
      robust_evolvability += 1;
  }

  // for (auto shape: shapes)
  //   shape.robust_pID(pIDs);

  uint8_t new_max = 0;

  for (auto pID: pIDs)
  {
    diversity.insert(pID);
    new_max = std::max(new_max, pID.first);
  }

  if(new_max > max_size)
    complex_evolvability += 1.;
}

void Genotype_Metrics::save_to_file(std::ofstream& fout)
{
  fout << "(";
  for (auto face: ref_genotype)
    fout <<+ face << ",";
  fout.seekp((long) fout.tellp() - 1);
  fout << ") ";

  fout << "(";
  for (auto face: original)
    fout <<+ face << ",";
  fout.seekp((long) fout.tellp() - 1);
  fout << ") ";

  for (auto pID: diversity)
    if (pID.first > max_size)
      complex_diversity++;

  fout <<+ strict_robustness / number_of_neighbours << " ";
  fout <<+ intersection_robustness / number_of_neighbours << " ";
  fout <<+ union_evolvability / number_of_neighbours << " ";
  fout <<+ complex_evolvability / number_of_neighbours << " ";
  fout <<+ robust_evolvability / number_of_neighbours << " ";
  fout <<+ rare / number_of_neighbours << " ";
  fout <<+ unbound / number_of_neighbours << " ";
  fout <<+ complex_diversity << " ";
  fout <<+ diversity.size() << " ";
  fout <<+ neutral_weight << " ";

  fout << "(";
  for(auto paired: pID_counter)
    fout <<+ paired.second << ",";
  fout.seekp((long) fout.tellp() - 1);
  fout << ") ";

  fout << "{";
  for (auto paired: pID_counter)
    fout <<+ "(" <<+ paired.first.first << "," <<+ paired.first.second << "),";
  fout.seekp((long) fout.tellp() - 1);
  fout << "}\n";

  // fout << "{";
  // for (auto pID: ref_pIDs)
  //   fout <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
  // fout.seekp((long) fout.tellp() - 1);
  // fout << "}\n";
}

void Genotype_Metrics::clear()
{
  strict_robustness = 0, intersection_robustness = 0, union_evolvability = 0;
  robust_evolvability = 0, complex_evolvability = 0, complex_diversity = 0;
  rare = 0, unbound = 0;

  // shapes.clear();
  diversity.clear();
}

Shape_Metrics::Shape_Metrics(Phenotype_ID pID): pID(pID), robustness(0)
{}

void Shape_Metrics::robust_pID(std::vector <Phenotype_ID> pIDs)
{
  auto presence = std::find(std::begin(pIDs), std::end(pIDs), pID);

  if (presence != std::end(pIDs))
    robustness++;
}


// Constructor of the Set_Metrics structure
Set_Metrics::Set_Metrics(uint16_t n_genes, uint16_t colours):
n_genes(n_genes), colours(colours), analysed(0)
{
  diversity_tracker.emplace_back(0);
}

// Member functions of the Set_Metrics structure

void Set_Metrics::add_genotype_metrics(Genotype_Metrics& genome_metric)
{
  double number_of_neighbours = (colours - 1) * n_genes * 4.;

  analysed++;

  strict_robustnesses.emplace_back(genome_metric.strict_robustness / number_of_neighbours);
  intersection_robustnesses.emplace_back(genome_metric.intersection_robustness / number_of_neighbours);
  union_evolvabilities.emplace_back(genome_metric.union_evolvability / number_of_neighbours);
  robust_evolvabilities.emplace_back(genome_metric.complex_evolvability / number_of_neighbours);
  complex_evolvabilities.emplace_back(genome_metric.robust_evolvability / number_of_neighbours);
  complex_diversities.emplace_back(genome_metric.complex_diversity);
  rares.emplace_back(genome_metric.rare / number_of_neighbours);
  unbounds.emplace_back(genome_metric.unbound / number_of_neighbours);

  neutral_weightings.emplace_back(genome_metric.neutral_weight);

  for (auto pID: genome_metric.diversity)
    diversity.insert(pID);
  diversity_tracker.emplace_back(diversity.size());

  genome_metrics.emplace_back(genome_metric);
}

void Set_Metrics::save_to_file(std::ofstream& set_out, std::ofstream& genome_out)
{
  for(auto metric: genome_metrics)
    metric.save_to_file(genome_out);

  double total_neutral_size = std::accumulate(neutral_weightings.begin(), neutral_weightings.end(), uint64_t(0));

  double average_strict_robustness = std::inner_product(std::begin(strict_robustnesses), std::end(strict_robustnesses), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_intersection_robustness = std::inner_product(std::begin(intersection_robustnesses), std::end(intersection_robustnesses), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_union_evolvability = std::inner_product(std::begin(union_evolvabilities), std::end(union_evolvabilities), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_complex_evolvability = std::inner_product(std::begin(complex_evolvabilities), std::end(complex_evolvabilities), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_robust_evolvability = std::inner_product(std::begin(robust_evolvabilities), std::end(robust_evolvabilities), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_complex_diversity = std::inner_product(std::begin(complex_diversities), std::end(complex_diversities), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_rare = std::inner_product(std::begin(rares), std::end(rares), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_unbound = std::inner_product(std::begin(unbounds), std::end(unbounds), std::begin(neutral_weightings), 0) / total_neutral_size;

  set_out <<+ average_strict_robustness << " " <<+ average_intersection_robustness << " ";
  set_out <<+ average_union_evolvability << " ";
  set_out <<+ average_robust_evolvability << " ";
  set_out <<+ average_complex_evolvability << " ";
  set_out <<+ average_complex_diversity << " ";
  set_out <<+ average_rare << " " <<+ average_unbound << " ";
  set_out <<+ analysed << " " <<+ misclassified.size() << " ";
  set_out <<+ total_neutral_size << " " << diversity.size() << " ";

  set_out << "(";
  for(auto value: diversity_tracker)
    set_out <<+ value << ",";
  set_out.seekp((long) set_out.tellp() - 1);
  set_out << ") ";

  set_out << "{";
  for (auto original: originals)
  {
    set_out << "(";
    for (auto face: original)
      set_out <<+ face << ",";
    set_out.seekp((long) set_out.tellp() - 1);
    set_out << "),";
  }
  set_out.seekp((long) set_out.tellp() - 1);
  set_out << "} ";

  set_out << "[";
  if(misclassified.size()) {
    for (auto mishap: misclassified)
    {
      set_out << "(";
      for (auto face: mishap.first)
        set_out <<+ face << ",";
      set_out.seekp((long) set_out.tellp() - 1);
      set_out << "),{";
      for (auto pID: mishap.second)
        set_out <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
      set_out.seekp((long) set_out.tellp() - 1);
      set_out << "},";
    }
    set_out.seekp((long) set_out.tellp() - 1);
    set_out << "] ";
  } else {
    set_out << "] ";
  }

  set_out << "{";
  for (auto pID: ref_pIDs)
    set_out <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
  set_out.seekp((long) set_out.tellp() - 1);
  set_out << "}\n";
}

void Set_Metrics::clear()
{
  strict_robustnesses.clear(), union_evolvabilities.clear(), rares.clear();
  robust_evolvabilities.clear(), complex_evolvabilities.clear();
  intersection_robustnesses.clear(), neutral_weightings.clear();
  complex_diversities.clear(), unbounds.clear();

  genome_metrics.clear();
  diversity.clear();
}
