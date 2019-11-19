#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <exception>

#include "cxxargs/include/cxxargs.hpp"

#include "file.hpp"

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
  File::In file(path);
  if (file.stream().good()) {
    std::string line;
    while (getline(file.stream(), line)) {
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
  } else {
    throw std::runtime_error("File " + path + " is not readable.");
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

void WriteFullAlignments(const std::string &path, const std::unordered_map<std::vector<bool>, ec_info> &ecs) {
  File::Out ec_file(path + ".ec");
  File::Out tsv_file(path + ".tsv");
  if (ec_file.stream().good() && tsv_file.stream().good()) {
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
  } else {
    throw std::runtime_error("Pseudoalignment files in " + path + " are not readable.");
  }
}

void parse_args(int argc, char* argv[], cxxargs::Arguments &args) {
  args.add_short_argument<std::string>('1', "Pseudoalignments for strand 1.");
  args.add_short_argument<std::string>('2', "Pseudoalignments for strand 2.");
  args.add_short_argument<std::string>('o', "Output files prefix.");
  args.add_long_argument<uint32_t>("n-refs", "Number of reference sequences in the pseudoalignment.");
  args.parse(argc, argv);
}

int main(int argc, char* argv[]) {
  cxxargs::Arguments args("telescope-unpaired", "Usage:...");
  try {
    std::cerr << "Parsing arguments" << '\n';
    parse_args(argc, argv, args);
  } catch (std::exception &e) {
    std::cerr << "Parsing arguments failed:"
	      << '\t' << e.what()
	      << "\nexiting" << std::endl;
    return 0;
  }
  std::unordered_map<std::vector<bool>, ec_info> ecs;

  ReadFullAlignments(args.value<std::string>('1'), ecs, args.value<uint32_t>("n-refs"));
  ReadFullAlignments(args.value<std::string>('2'), ecs, args.value<uint32_t>("n-refs"));
  WriteFullAlignments(args.value<std::string>('o'), ecs);
  std::cerr << "Finished, exiting..." << std::endl;
  return 0;
}
