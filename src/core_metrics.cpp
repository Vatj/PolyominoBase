#include "core_metrics.hpp"

// Constructor of the Genotype_Metrics structure
Genotype_Metrics::Genotype_Metrics(uint8_t n_genes, uint8_t colours):
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

  neutral_weight = neutral;

  // for (auto pID: pIDs)
  //   shapes.emplace_back(Shape_Metrics(pID));
}

void Genotype_Metrics::analyse_pIDs(std::vector <Phenotype_ID>& pIDs)
{
  if (pIDs == ref_pIDs)
    strict_robustness+= 1.;

  std::vector <Phenotype_ID> intersection, union_set;

  std::set_intersection(std::begin(pIDs), std::end(pIDs), std::begin(ref_pIDs), std::end(ref_pIDs), std::back_inserter(intersection));
  std::set_union(std::begin(pIDs), std::end(pIDs), std::begin(ref_pIDs), std::end(ref_pIDs), std::back_inserter(union_set));

  if(ref_pIDs.size() > 0)
    intersection_robustness += (double) intersection.size() / (double) ref_pIDs.size();
  else
    intersection_robustness = 0;

  union_evolvability += (double) (union_set.size() - ref_pIDs.size());

  if (std::find(std::begin(pIDs), std::end(pIDs), rare_pID) != std::end(pIDs))
    rare += 1.;

  if (std::find(std::begin(pIDs), std::end(pIDs), loop_pID) != std::end(pIDs))
    loop += 1.;

  // for (auto shape: shapes)
  //   shape.robust_pID(pIDs);

  for (auto pID: pIDs)
    diversity.insert(pID);
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

  fout <<+ strict_robustness / number_of_neighbours << " ";
  fout <<+ intersection_robustness / number_of_neighbours << " ";
  fout <<+ union_evolvability / number_of_neighbours << " ";
  fout <<+ (union_evolvability - rare - loop) / number_of_neighbours << " ";
  fout <<+ rare / number_of_neighbours << " ";
  fout <<+ loop / number_of_neighbours << " ";
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
  strict_robustness = 0, union_evolvability = 0, rare = 0;
  intersection_robustness = 0, loop = 0;

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
Set_Metrics::Set_Metrics(uint8_t n_genes, uint8_t colours):
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
  rares.emplace_back(genome_metric.rare / number_of_neighbours);
  loops.emplace_back(genome_metric.loop / number_of_neighbours);

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
  double average_rare = std::inner_product(std::begin(rares), std::end(rares), std::begin(neutral_weightings), 0) / total_neutral_size;
  double average_loop = std::inner_product(std::begin(loops), std::end(loops), std::begin(neutral_weightings), 0) / total_neutral_size;

  set_out <<+ average_strict_robustness << " " <<+ average_intersection_robustness << " ";
  set_out <<+ average_union_evolvability << " ";
  set_out <<+ average_union_evolvability - average_loop - average_rare << " ";
  set_out <<+ average_rare << " " <<+ average_loop << " ";
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
    set_out << "}" << std::endl;
  }
  set_out.seekp((long) set_out.tellp() - 1);
  set_out << "] ";

  set_out << "{";
  for (auto pID: ref_pIDs)
    set_out <<+ "(" <<+ pID.first << "," <<+ pID.second << "),";
  set_out.seekp((long) set_out.tellp() - 1);
  set_out << "}" << std::endl;
}

void Set_Metrics::clear()
{
  strict_robustnesses.clear(), union_evolvabilities.clear(), rares.clear();
  intersection_robustnesses.clear(), neutral_weightings.clear(), loops.clear();

  genome_metrics.clear();
  diversity.clear();
}
