#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

struct ec_info {
  std::vector<bool> pseudoalignment;
  uint32_t count = 0;
  uint16_t last_val = 0;
};

class Pseudoalignment {
public:
  std::vector<bool> alignment;
  std::vector<uint32_t> read_ids;
  uint16_t last_val = 0;
  Pseudoalignment(std::vector<bool> alignment);
};

struct PseudoalignmentHash {
  size_t operator()(const Pseudoalignment &obj) const {
    return std::hash<std::vector<bool>>()(obj.alignment);
  }
};

struct PseudoalignmentEqual_to {
  bool operator()(const Pseudoalignment &obj1, const Pseudoalignment &obj2) const {
    return obj1.alignment == obj2.alignment;
  }
};

void ReadFullAlignments(const std::string &path, std::unordered_map<std::vector<bool>, ec_info> &ecs, uint32_t n_refs) {
  std::ifstream file(path);
  if (file.is_open()) {
    std::string line;
    while (getline(file, line)) {
      std::vector<bool> alignment(n_refs, false);
      std::string part;
      std::stringstream partition(line);
      getline(partition, part, ' ');
      uint32_t ec_id = std::stoi(part);
      uint16_t cluster_id;
      bool any_aligned = false;
      while (getline(partition, part, ' ')) {
	cluster_id = std::stoi(part);
	alignment[cluster_id] = true;
	any_aligned = true;
      }
      if (any_aligned) {
	ecs[alignment].count += 1;
	ecs[alignment].last_val = std::max(cluster_id, ecs[alignment].last_val);
      }
    }
  }
}

void ReadIntersectionAlignments(const std::string &path, std::unordered_map<uint32_t, ec_info> &ecs, uint32_t n_refs) {
  std::ifstream file(path);
  if (file.is_open()) {
    std::string line;
    while (getline(file, line)) {
      std::vector<bool> alignment(n_refs, false);
      std::string part;
      std::stringstream partition(line);
      getline(partition, part, ' ');
      uint32_t ec_id = std::stoi(part);
      uint16_t cluster_id;
      bool any_aligned = false;
      while (getline(partition, part, ' ')) {
	cluster_id = std::stoi(part);
	alignment[cluster_id] = true;
	any_aligned = true;
      }
      if (any_aligned) {
	if (ecs.find(ec_id) == ecs.end()) {
	  ecs[ec_id].pseudoalignment = alignment;
	  ecs[ec_id].count = 1;
	  ecs[ec_id].last_val = std::max(cluster_id, ecs[ec_id].last_val);
	} else {
	  any_aligned = false;
	  for (size_t i = 0; i < n_refs; ++i) {
	    ecs[ec_id].pseudoalignment[i] = ecs[ec_id].pseudoalignment[i] && alignment[i];
	    ecs[ec_id].last_val = (ecs[ec_id].pseudoalignment[i] ? i : ecs[ec_id].last_val);
	    any_aligned = ecs[ec_id].pseudoalignment[i] || any_aligned;
	  }
	  if (!any_aligned) {
	    ecs.erase(ec_id);
	  }
	}
      }
    }
  }
}

void ReadUnionAlignments(const std::string &path, std::unordered_map<uint32_t, ec_info> &ecs, uint32_t n_refs) {
  std::ifstream file(path);
  if (file.is_open()) {
    std::string line;
    while (getline(file, line)) {
      std::vector<bool> alignment(n_refs, false);
      std::string part;
      std::stringstream partition(line);
      getline(partition, part, ' ');
      uint32_t ec_id = std::stoi(part);
      uint16_t cluster_id;
      bool any_aligned = false;
      while (getline(partition, part, ' ')) {
	cluster_id = std::stoi(part);
	alignment[cluster_id] = true;
	any_aligned = true;
      }
      if (any_aligned) {
	if (ecs.find(ec_id) == ecs.end()) {
	  ecs[ec_id].pseudoalignment = alignment;
	  ecs[ec_id].count = 1;
	  ecs[ec_id].last_val = std::max(cluster_id, ecs[ec_id].last_val);
	} else {
	  any_aligned = false;
	  for (size_t i = 0; i < n_refs; ++i) {
	    ecs[ec_id].pseudoalignment[i] = ecs[ec_id].pseudoalignment[i] || alignment[i];
	    ecs[ec_id].last_val = (ecs[ec_id].pseudoalignment[i] ? i : ecs[ec_id].last_val);
	    any_aligned = ecs[ec_id].pseudoalignment[i] || any_aligned;
	  }
	  if (!any_aligned) {
	    ecs.erase(ec_id);
	  }
	}
      }
    }
  }
}

std::unordered_map<std::vector<bool>, ec_info> CompressAlignments(const std::unordered_map<uint32_t, ec_info> &ecs) {
  std::unordered_map<std::vector<bool>, ec_info> compressed_ecs;
  for (auto kv : ecs) {
    if (compressed_ecs.find(kv.second.pseudoalignment) == compressed_ecs.end()) {
      compressed_ecs[kv.second.pseudoalignment].count = 1;
    } else {
      compressed_ecs[kv.second.pseudoalignment].count += 1;
    }
    compressed_ecs[kv.second.pseudoalignment].last_val = kv.second.last_val;
  }
  return compressed_ecs;
}

// void ReadFullAlignments2(const std::string &path, std::unordered_map<Pseudoalignment, uint32_t> &ecs) {
//   Reader reader(path);
//   if (reader.is_ready()) {
//     while (reader.next_line()) {
//       reader.partition_line();
//       reader.next_part(' ');
//       std::vector<bool> alignment(3815, false);
//       uint32_t ec_id = std::stoi(reader.current_part());
//       uint16_t cluster_id;
//       bool any_aligned = false;
//       while (reader.next_part(' ')) {
// 	cluster_id = std::stoi(reader.current_part());
// 	alignment[cluster_id] = true;
// 	any_aligned = true;
//       }
//       if (any_aligned) {
// 	Pseudoalignment pseudo(alignment);

// 	ecs[alignment].count += 1;
// 	ecs[alignment].last_val = std::max(cluster_id, ecs[alignment].last_val);
//       }
//     }
//   }
// }

void WriteFullAlignments(const std::string &path, const std::unordered_map<std::vector<bool>, ec_info> &ecs) {
  std::ofstream ec_file(path + ".ec");
  std::ofstream tsv_file(path + ".tsv");
  if (ec_file.is_open() && tsv_file.is_open()) {
    uint32_t ec_id = 0;
    for (auto ec : ecs) {
      ec_file << ec_id << '\t';
      tsv_file << ec_id << '\t';
      tsv_file << ec.second.count << '\n';
      for (size_t i = 0; i < ec.first.size(); ++i) {
	if (ec.first.at(i)) {
	  ec_file << i << (i == ec.second.last_val ? '\n' : ',');
	}
      }
      ++ec_id;
    }
    ec_file.close();
    tsv_file.close();
  }
}

// void ReadFullAlignments(const std::string &path, std::unordered_set<Pseudoalignment> &ecs) {
//   std::ifstream file(path);
//   if (file.is_open()) {
//     std::string line;
//     while (getline(file, line)) {
//       std::vector<bool> alignment(3815, false);
//       std::string part;
//       std::stringstream partition(line);
//       getline(partition, part, ' ');
//       short unsigned ec_id = std::stoi(part);
//       short unsigned cluster_id;
//       bool any_aligned = false;
//       while (getline(partition, part, ' ')) {
// 	cluster_id = std::stoi(part);
// 	alignment[cluster_id] = true;
// 	any_aligned = true;
//       }
//       if (any_aligned) {
// 	Pseudoalignment pseudo;
// 	ecs.insert(pseudo);
// 	//	ecs[alignment].alignment = alignment;
// 	//	ecs[alignment].count += 1;
// 	//	ecs[alignment].last_val = cluster_id;
//       }
//     }
//   }
// }

// void WriteFullAlignments(const std::string &path, const std::unordered_set<Pseudoalignment> &ecs) {
//   std::ofstream ec_file(path + ".ec");
//   std::ofstream tsv_file(path + ".tsv");
//   if (ec_file.is_open() && tsv_file.is_open()) {
//     short unsigned ec_id = 0;
//     for (auto ec : ecs) {
//       ec_file << ec_id << '\t';
//       tsv_file << ec_id << '\t';
//       tsv_file << ec.count << '\n';
//       for (size_t i = 0; i < ec.alignment.size(); ++i) {
// 	if (ec.alignment.at(i)) {
// 	  ec_file << i << (i == ec.last_val ? '\n' : ',');
// 	}
//       }
//       ++ec_id;
//     }
//     ec_file.close();
//     tsv_file.close();
//   }
// }
 
int main(int argc, char* argv[]) {
  std::string infile_1 = argv[1];
  std::string infile_2 = argv[2];
  std::string outfile = argv[3];
  uint32_t n_refs = std::stoi(argv[4]);
  std::unordered_map<std::vector<bool>, ec_info> ecs;
  // std::unordered_map<uint32_t, ec_info> ecs;
  // std::unordered_set<Pseudoalignment, PseudoalignmentHash, PseudoalignmentEqual_to> ec_pseudos;

  ReadFullAlignments(infile_1, ecs, n_refs);
  ReadFullAlignments(infile_2, ecs, n_refs);
  // ReadIntersectionAlignments(infile_1, ecs, n_refs);
  // ReadIntersectionAlignments(infile_2, ecs, n_refs);
  // const std::unordered_map<std::vector<bool>, ec_info> &cmp_ecs = CompressAlignments(ecs);
  WriteFullAlignments(outfile, ecs);
  std::cout << "Finished, exiting..." << std::endl;
  return 0;
}
