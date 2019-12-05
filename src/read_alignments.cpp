#include "read_alignments.hpp"

#include <utility>
#include <algorithm>
#include <sstream>
#include <exception>

void read_unpaired(const uint16_t cluster_id, const uint32_t ec_id, std::unordered_map<uint32_t, ec_info>* ecs) {
  ecs->at(ec_id).last_val = std::max(cluster_id, ecs->at(ec_id).last_val);
}

void read_paired(const Mode &mode, const std::vector<bool> &alignment, const uint32_t ec_id, const uint32_t n_refs, bool* any_aligned, std::unordered_map<uint32_t, ec_info>* ecs) {
  *any_aligned = false;
  for (size_t i = 0; i < n_refs; ++i) {
    if (mode == m_intersection) {
      ecs->at(ec_id).pseudoalignment[i] = ecs->at(ec_id).pseudoalignment[i] && alignment[i];
    } else if (mode == m_union) {
      ecs->at(ec_id).pseudoalignment[i] = ecs->at(ec_id).pseudoalignment[i] || alignment[i];
    }		      
    ecs->at(ec_id).last_val = (ecs->at(ec_id).pseudoalignment[i] ? i : ecs->at(ec_id).last_val);
    *any_aligned = ecs->at(ec_id).pseudoalignment[i] || *any_aligned;
  }
  if (!(*any_aligned)) {
    ecs->erase(ec_id);
  }
}

void read_alignments(const Mode &mode, const uint32_t n_refs, std::istream *stream, std::unordered_map<uint32_t, ec_info>* ecs) {
  std::string line;
  while (getline(*stream, line)) {
    std::vector<bool> alignment(n_refs, false);
    std::string part;
    std::stringstream partition(line);
    getline(partition, part, ' ');
    uint32_t ec_id = std::stoul(part);
    uint16_t cluster_id;
    bool any_aligned = false;
    while (getline(partition, part, ' ')) {
      cluster_id = std::stoul(part);
      alignment[cluster_id] = true;
      any_aligned = true;
    }
    if (any_aligned) {
      if (ecs->find(ec_id) == ecs->end()) {
	ec_info info;
	info.pseudoalignment = alignment;
	info.count = 1;
	info.last_val = cluster_id;
	ecs->insert(std::make_pair(ec_id, info));
      } else {
	ecs->at(ec_id).count += 1;
	switch (mode) {
	 case m_unpaired :
	  read_unpaired(cluster_id, ec_id, ecs);
	  break;
	 default :
	  read_paired(mode, alignment, ec_id, n_refs, &any_aligned, ecs);
	  break;
	}
      }
    }
  }
}

void CompressAlignments(const std::unordered_map<uint32_t, ec_info> &ecs, std::unordered_map<std::vector<bool>, ec_info>* compressed_ecs) {
  for (auto kv : ecs) {
    if (compressed_ecs->find(kv.second.pseudoalignment) == compressed_ecs->end()) {
      ec_info info;
      info.pseudoalignment = kv.second.pseudoalignment;
      info.count = 1;
      info.last_val = kv.second.last_val;
      compressed_ecs->insert(std::make_pair(kv.second.pseudoalignment, info));
    } else {
      compressed_ecs->at(kv.second.pseudoalignment).count += 1;
    }
    compressed_ecs->at(kv.second.pseudoalignment).last_val = kv.second.last_val;
  }
}

void ReadToRef(const std::unordered_map<uint32_t, ec_info> &ecs, std::unordered_map<uint32_t, std::vector<uint16_t>>* read_to_ref) {
  for (auto ec : ecs) {
    if (read_to_ref->find(ec.first) == read_to_ref->end()) {
      read_to_ref->insert(std::make_pair(ec.first, std::vector<uint16_t>()));
    }
    for (size_t i = 0; i < ec.second.pseudoalignment.size(); ++i) {
      if (ec.second.pseudoalignment.at(i)) {
	read_to_ref->at(ec.first).push_back(i);
      }
    }
  }
}

KAlignment ReadAlignments(const Mode &mode, const uint32_t n_refs, std::istream* strand_1, std::istream* strand_2) {
  std::unordered_map<uint32_t, ec_info> ecs_by_id;
  read_alignments(mode, n_refs, strand_1, &ecs_by_id);
  read_alignments(mode, n_refs, strand_2, &ecs_by_id);

  KAlignment alignments;
  ReadToRef(ecs_by_id, &alignments.read_to_ref);
  CompressAlignments(ecs_by_id, &alignments.ecs);

  return alignments;
}
