#include <stochastic_model.hpp>

namespace Stochastic
{

  std::vector<Phenotype_ID> AssemblePlasticGenotype(Genotype genotype, PhenotypeTable* pt) {

    std::vector<Phenotype_ID> plastic_phenotypes;
    uint16_t rare_phenotypes=0;
    std::map<Phenotype_ID,uint8_t> ID_counter;

    for(uint8_t seed=0;seed<genotype.size()/4;++seed) {
      std::vector<Phenotype> raw_phenotypes;raw_phenotypes.reserve(simulation_params::phenotype_builds);
      std::vector<Phenotype_ID> Phenotype_IDs;Phenotype_IDs.reserve(simulation_params::phenotype_builds);

      for(uint8_t kth=0;kth<simulation_params::phenotype_builds;++kth) {
	       std::vector<int8_t> placed_tiles=Stochastic_Polyomino_Builder(genotype, seed);
	        if(placed_tiles.size()>0)
	         raw_phenotypes.emplace_back(Generate_Spatial_Occupancy(placed_tiles));
	        else {
            plastic_phenotypes.clear();
            plastic_phenotypes.emplace_back(std::make_pair(255,0));

            return plastic_phenotypes;
          }
      }

      #pragma omp critical(phenotype_lookup)
      {
	       for(Phenotype& phen : raw_phenotypes)
          Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));

	       pt->RelabelPhenotypes(Phenotype_IDs);
	       ID_counter=pt->PhenotypeFrequencies(Phenotype_IDs,rare_phenotypes);
      }

      for(auto kv : ID_counter) {
	       if(kv.second >=static_cast<uint16_t>(simulation_params::phenotype_builds*simulation_params::UND_threshold)) {
	          if(std::find(plastic_phenotypes.begin(),plastic_phenotypes.end(),kv.first)==plastic_phenotypes.end())
	           plastic_phenotypes.emplace_back(kv.first);
	       }
	       else
	        rare_phenotypes++;
      }
    }
    if(rare_phenotypes)
      plastic_phenotypes.emplace_back(std::make_pair(0,0));

    return plastic_phenotypes;
  }

  std::map<Phenotype_ID, uint8_t> AssemblePlasticGenotypeFrequency(Genotype genotype, PhenotypeTable* pt) {

    std::map<Phenotype_ID,uint8_t> pID_counter_plastic;
    uint16_t rare_phenotypes=0;
    std::map<Phenotype_ID,uint8_t> ID_counter;

    for(uint8_t seed=0;seed<genotype.size()/4;++seed) {
      std::vector<Phenotype> raw_phenotypes;raw_phenotypes.reserve(simulation_params::phenotype_builds);
      std::vector<Phenotype_ID> Phenotype_IDs;Phenotype_IDs.reserve(simulation_params::phenotype_builds);

      for(uint8_t kth=0;kth<simulation_params::phenotype_builds;++kth) {
	       std::vector<int8_t> placed_tiles=Stochastic_Polyomino_Builder(genotype, seed);
	        if(placed_tiles.size()>0)
	         raw_phenotypes.emplace_back(Generate_Spatial_Occupancy(placed_tiles));
	        else {
            pID_counter_plastic.clear();
            pID_counter_plastic[std::make_pair(255,0)] = 0;

            return pID_counter_plastic;
          }
      }

      #pragma omp critical(phenotype_lookup)
      {
	       for(Phenotype& phen : raw_phenotypes)
          Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));

	       pt->RelabelPhenotypes(Phenotype_IDs);
	       ID_counter=pt->PhenotypeFrequencies(Phenotype_IDs, rare_phenotypes);
      }

      for(auto kv : ID_counter) {
	       if(kv.second >=static_cast<uint16_t>(simulation_params::phenotype_builds*simulation_params::UND_threshold)) {
           if(pID_counter_plastic.count(kv.first))
            pID_counter_plastic[kv.first] += kv.second;
          else
            pID_counter_plastic[kv.first] = kv.second;
	       }
	       else
	        rare_phenotypes++;
      }
    }

    if(rare_phenotypes)
      pID_counter_plastic[std::make_pair(0,0)] = rare_phenotypes;

    return pID_counter_plastic;
  }

  bool FastNoPIDs2(Genotype genotype, PhenotypeTable* pt) {
    /* Return true if unbound or rare, interrupting the process */
    std::map<Phenotype_ID,uint8_t> pID_counter_plastic;
    std::map<Phenotype_ID,uint8_t> ID_counter;
    uint16_t rare_phenotypes;
    // Clear all known shapes from previous iteration
    (pt->known_phenotypes).clear();

    // Looping the building using each gene as a seed
    for(uint8_t seed=0;seed<genotype.size()/4;++seed) {
      std::vector<Phenotype> raw_phenotypes;raw_phenotypes.reserve(simulation_params::phenotype_builds);
      std::vector<Phenotype_ID> Phenotype_IDs;Phenotype_IDs.reserve(simulation_params::phenotype_builds);

      // Build everything first, before matching with pIDs
      for(uint8_t kth=0;kth<simulation_params::phenotype_builds;++kth) {
	       std::vector<int8_t> placed_tiles=Stochastic_Polyomino_Builder(genotype, seed);
	        if(placed_tiles.size()>0)
	         raw_phenotypes.emplace_back(Generate_Spatial_Occupancy(placed_tiles));
	        else
            // If the genotype is unbound it is rejected
	         return true;
      }

      // Match shapes with pID and compute frequencies
      #pragma omp critical(phenotype_lookup)
      {
	       for(Phenotype& phen : raw_phenotypes)
          Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));

         pt->RelabelPhenotypes(Phenotype_IDs);
	       ID_counter=pt->PhenotypeFrequencies(Phenotype_IDs, rare_phenotypes);
      }
      if(rare_phenotypes)
        return true;
      }
    return false;
  }

  bool FastNoPIDs(Genotype genotype, PhenotypeTable* pt) {

    std::vector<Phenotype_ID> plastic_phenotypes;
    uint16_t rare_phenotypes=0;
    std::map<Phenotype_ID,uint8_t> ID_counter;
    (pt->known_phenotypes).clear();

    for(uint8_t seed=0;seed<genotype.size()/4;++seed) {
      std::vector<Phenotype> raw_phenotypes;raw_phenotypes.reserve(simulation_params::phenotype_builds);
      std::vector<Phenotype_ID> Phenotype_IDs;Phenotype_IDs.reserve(simulation_params::phenotype_builds);

      for(uint8_t kth=0;kth<simulation_params::phenotype_builds;++kth) {
	       std::vector<int8_t> placed_tiles=Stochastic_Polyomino_Builder(genotype, seed);
	        if(placed_tiles.size()>0)
	         raw_phenotypes.emplace_back(Generate_Spatial_Occupancy(placed_tiles));
	        else
	         return true;
      }

      #pragma omp critical(phenotype_lookup)
      {
	       for(Phenotype& phen : raw_phenotypes)
          Phenotype_IDs.emplace_back(pt->GetPhenotypeID(phen));

	       pt->RelabelPhenotypes(Phenotype_IDs);
	       ID_counter=pt->PhenotypeFrequencies(Phenotype_IDs, rare_phenotypes);
      }

      for(auto kv : ID_counter) {
	       if(kv.second >=static_cast<uint16_t>(simulation_params::phenotype_builds*simulation_params::UND_threshold)) {
	          if(std::find(plastic_phenotypes.begin(),plastic_phenotypes.end(),kv.first)==plastic_phenotypes.end())
	           plastic_phenotypes.emplace_back(kv.first);
	       }
	       else
	        return true;
      }
    }
    if(rare_phenotypes)
      return true;

    return false;
  }

  std::vector<int8_t> Stochastic_Polyomino_Builder(const Genotype& genome, uint8_t initial_Tile) {
    const uint8_t THRESHOLD_SIZE=(genome.size()*genome.size())/4;
    std::vector<int8_t> Placed_Tiles{0,0,static_cast<int8_t>(initial_Tile),0},next_binds; //DEFINED AS (X,Y,Tile Type Number, Tile Rotation[in CW rotation])
    for(uint8_t face=0; face<4; ++face)
      if(genome[initial_Tile*4+face]!=0)
        Interacting_Adjacency(next_binds,genome[initial_Tile*4+face],face,0,0);

    while(!next_binds.empty()) {
      for(int Noptions=(next_binds.size()/4)-1; Noptions >=0; --Noptions) {
        bool overlapped=false;
        for(std::vector<int8_t>::iterator occupied_iter = Placed_Tiles.begin(); occupied_iter!=Placed_Tiles.end(); occupied_iter+=4) {
          if(next_binds[Noptions*4]==*occupied_iter && next_binds[Noptions*4+1]==*(occupied_iter+1)) {
            if(simulation_params::STERIC_FORBIDDEN && next_binds[Noptions*4+2]!=genome[*(occupied_iter+2)*4+(next_binds[Noptions*4+3]-*(occupied_iter+3)+4)%4])
              return {};
            next_binds.erase(next_binds.end()-4,next_binds.end());
            overlapped=true;
            break;
          }
        }
        if(overlapped)
          continue;

        uint8_t conjugate_count=std::count(genome.begin(),genome.end(),next_binds[Noptions*4+2]);
        std::uniform_int_distribution<uint8_t> Random_Count(0,conjugate_count-1);
        int nth_conjugate=Random_Count(simulation_params::RNG_Engine);
        auto current_conjugate=std::find(genome.begin(),genome.end(),next_binds[Noptions*4+2]);
        for(uint8_t conj_cnt=1;conj_cnt<=nth_conjugate;++conj_cnt)
          current_conjugate=std::find(current_conjugate+1,genome.end(),next_binds[Noptions*4+2]);
        uint8_t new_Tile=(current_conjugate-genome.begin())/4,new_Face=(current_conjugate-genome.begin())%4,rotation= (next_binds[Noptions*4+3]-new_Face+4)%4;

        Placed_Tiles.insert(Placed_Tiles.end(),{next_binds[Noptions*4],next_binds[Noptions*4+1],static_cast<int8_t>(new_Tile),static_cast<int8_t>(rotation)});
        if(Placed_Tiles.size()/4 > THRESHOLD_SIZE)
          return {};//Placed_Tiles;
        uint8_t placed_X=next_binds[Noptions*4],placed_Y=next_binds[Noptions*4+1];
        next_binds.erase(next_binds.begin()+Noptions*4,next_binds.begin()+Noptions*4+4);
        for(uint8_t face=1;face<4;++face) {
          uint8_t temp_Face=genome[new_Tile*4+(new_Face+face)%4];
          if(temp_Face!=0)
            Interacting_Adjacency(next_binds,temp_Face,(new_Face+face+rotation)%4,placed_X,placed_Y);
        }
        break;
      }
    }
    return Placed_Tiles;
  }

  void Interacting_Adjacency(std::vector<int8_t>& next_binds, uint8_t interacting_Face, uint8_t face_index, int8_t X, int8_t Y) {
    int8_t X_OFFSET=0,Y_OFFSET=0;
    switch(face_index) {
    case 0:
      X_OFFSET=0;Y_OFFSET=1;
      break;
    case 1:
      X_OFFSET=1;Y_OFFSET=0;
      break;
    case 2:
      X_OFFSET=0;Y_OFFSET=-1;
      break;
    case 3:
      X_OFFSET=-1;Y_OFFSET=0;
      break;
    }
    uint8_t conjugate_Face=Interaction_Matrix(interacting_Face);
    next_binds.insert(next_binds.begin()+std::uniform_int_distribution<size_t>{0,next_binds.size()/4}(simulation_params::RNG_Engine)*4,{static_cast<int8_t>(X+X_OFFSET),static_cast<int8_t>(Y+Y_OFFSET),static_cast<int8_t>(conjugate_Face),static_cast<int8_t>((face_index+2)%4)});
  }


  Phenotype Generate_Spatial_Occupancy(std::vector<int8_t>& Placed_Tiles_Check) {
    std::vector<int8_t> X_Locs_Check, Y_Locs_Check;
    std::vector<uint8_t>  Tile_Type_Check,Tile_Orientation_Check;
    for(std::vector<int8_t>::iterator check_iter = Placed_Tiles_Check.begin();check_iter!=Placed_Tiles_Check.end();check_iter+=4) {
      X_Locs_Check.emplace_back(*check_iter);
      Y_Locs_Check.emplace_back(*(check_iter+1));
      Tile_Type_Check.emplace_back(*(check_iter+2));
      Tile_Orientation_Check.emplace_back(*(check_iter+3));
    }
    std::vector<int8_t>::iterator lx,rx,ty,by;
    std::tie(lx,rx)=std::minmax_element(X_Locs_Check.begin(),X_Locs_Check.end());
    std::tie(by,ty)=std::minmax_element(Y_Locs_Check.begin(),Y_Locs_Check.end());
    uint8_t dx=*rx-*lx+1;
    uint8_t dy=*ty-*by+1;
    std::vector<uint8_t> Spatial_Occupancy_Check(dx*dy);

    for(uint8_t tileIndex=0;tileIndex<X_Locs_Check.size();++tileIndex)
    	Spatial_Occupancy_Check[(*ty-Y_Locs_Check[tileIndex])*dx + (X_Locs_Check[tileIndex]-*lx)]=1+Tile_Type_Check[tileIndex]*(DETERMINISM_LEVEL==2+4*(DETERMINISM_LEVEL==3))+Tile_Orientation_Check[tileIndex]*(DETERMINISM_LEVEL==4);
    return Phenotype{dx,dy,Spatial_Occupancy_Check};
  }

}
