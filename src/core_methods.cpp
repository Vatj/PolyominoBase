#include "core_methods.hpp"

uint8_t Interaction_Matrix(uint8_t input_face) {
  return input_face>0 ?  (1-input_face%2)*(input_face-1)+(input_face%2)*(input_face+1) : input_face;
}

void Clean_Genome(std::vector<uint8_t>& genome,int secondNonInteracting=0,bool Remove_Duplicates=true) {
  //defaults set to -1 and true//
  //Removes any non-bonding interfaces and reduces additional non-interacting faces to 0s
  //Note, no need to check for negatives, as either present and able to self-interact, or absent and no need to act
  if(secondNonInteracting!=0) 
    std::replace(genome.begin(),genome.end(),secondNonInteracting,0);
  
  for(int t=1;t<=*std::max_element(genome.begin(),genome.end());t+=2) { 
    if(std::count(genome.begin(),genome.end(),t)==0) //genotype doens't contain this face
        std::replace(genome.begin(),genome.end(),t+1,0);
    else //genotype does contain this face
      if(std::find(genome.begin(),genome.end(),t+1)==genome.end()) //genotype doesn't contain conjugate
        std::replace(genome.begin(),genome.end(),t,0);
  }
  Minimize_Tile_Set(genome);
  if(Remove_Duplicates) 
    DuplicateGenes(genome);
}

void Minimize_Tile_Set(std::vector<uint8_t>& genome) {
  std::vector<uint8_t> Minimal_Genome;
  for(std::vector<uint8_t>::iterator genome_iter=genome.begin();genome_iter!=genome.end();genome_iter+=4) {
    Minimal_Genome.assign(genome_iter,genome_iter+4);
    if(std::count(Minimal_Genome.begin(),Minimal_Genome.end(),0)==4)
      continue;
    for(uint8_t rotation=1;rotation<=3;++rotation) {
      std::rotate(genome_iter,genome_iter+1,genome_iter+4);
      for(uint8_t index=0;index<4;++index) {
        if(Minimal_Genome[index]<*(genome_iter+index))
          break;
        if(Minimal_Genome[index]==*(genome_iter+index))
          continue;
        else {
          Minimal_Genome.assign(genome_iter,genome_iter+4);
          break;
        }
      }
    }
    std::vector<uint8_t>::iterator Minimal_beginning=Minimal_Genome.begin();
    for(uint8_t swapping=0;swapping<4;++swapping) {
      *(genome_iter+swapping)=*(Minimal_beginning+swapping);
    }
  }
}

std::map<uint8_t,uint8_t> DuplicateGenes(std::vector<uint8_t>& genome) {
  std::map<uint8_t,uint8_t> dups;
  for(int check_index=genome.size()/4-1;check_index>0;--check_index)
    for(int compare_index=0;compare_index<check_index;++compare_index)
      if(std::equal(genome.begin()+check_index*4,genome.begin()+check_index*4+4,genome.begin()+compare_index*4)) {
        genome.erase(genome.begin()+check_index*4,genome.begin()+check_index*4+4);
        dups[check_index]=compare_index;
        break;
      }
  return dups;
}

bool Disjointed_Check(std::vector<uint8_t>& genome) {
  std::vector<uint8_t> Unvisited(genome.size()/4-1);
  std::map<uint8_t, std::vector<uint8_t> > Connected_Components;
  Connected_Components[0].insert(Connected_Components[0].end(),genome.begin(),genome.begin()+4);
  std::iota(Unvisited.begin(),Unvisited.end(),1);
  Search_Next_Tile(genome,Unvisited,Connected_Components[0],0);
  while(!Unvisited.empty()) {
    uint8_t CC_Size=Connected_Components.size();
    uint8_t New_Launch_Point=Unvisited[0];
    Connected_Components[CC_Size].insert(Connected_Components[CC_Size].end(),genome.begin()+ New_Launch_Point*4,genome.begin()+New_Launch_Point*4+4);
    Unvisited.erase(Unvisited.begin());
    Search_Next_Tile(genome,Unvisited,Connected_Components[CC_Size], New_Launch_Point);
  }
  
  if(Connected_Components.size()>1) {
    genome.clear();
    for(auto& KV_pair : Connected_Components) {
      genome.insert(genome.end(),KV_pair.second.begin(),KV_pair.second.end());
      genome.emplace_back(-1);
    }
    return true;
  }
  return false;  
}

void Search_Next_Tile(std::vector<uint8_t>& genome, std::vector<uint8_t>& Unvisited, std::vector<uint8_t>& Connected_Components,uint8_t tile) {
  for(uint8_t face=0;face<4;++face) {
    if(genome[tile*4+face]==0)
      continue;
    uint8_t conjugate_Face=Interaction_Matrix(genome[tile*4+face]);
    std::vector<uint8_t>::iterator conjugate_iter=std::find(genome.begin(),genome.end(),conjugate_Face);
    while(conjugate_iter!=genome.end()) {
      uint8_t corresponding_Tile=(conjugate_iter-genome.begin())/4;
      if(std::find(Unvisited.begin(),Unvisited.end(),corresponding_Tile)!=Unvisited.end()) { //new tile
        Connected_Components.insert(Connected_Components.end(),genome.begin()+corresponding_Tile*4,genome.begin()+corresponding_Tile*4+4);
        Unvisited.erase(std::find(Unvisited.begin(),Unvisited.end(),corresponding_Tile));
        Search_Next_Tile(genome,Unvisited,Connected_Components,corresponding_Tile);
      } 
      conjugate_iter=std::find(conjugate_iter+1,genome.end(),conjugate_Face);
    }
  }
}
