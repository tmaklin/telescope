// telescope: convert between Themisto and kallisto pseudoalignments
// Copyright (C) 2019 Tommi Mäklin (tommi@maklin.fi)
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

#include "read_alignments.hpp"

#include <utility>
#include <algorithm>
#include <sstream>
#include <exception>

void insert_read(const std::vector<bool> &alignment, const uint16_t cluster_id, const uint32_t read_id, std::unordered_map<uint32_t, ec_info>* reads) {
  ec_info info;
  info.pseudoalignment = alignment;
  info.count = 1;
  info.last_val = cluster_id;
  reads->insert(std::make_pair(read_id, info));
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

void read_alignments(const Mode &mode, const uint32_t n_refs, std::istream *stream, std::unordered_map<uint32_t, ec_info>* reads, uint32_t *max_read_id) {
  std::string line;
  while (getline(*stream, line)) {
    std::vector<bool> alignment(n_refs, false);
    std::string part;
    std::stringstream partition(line);
    getline(partition, part, ' ');
    uint32_t read_id = std::stoul(part);
    *max_read_id = (read_id > *max_read_id ? read_id : *max_read_id);
    uint16_t cluster_id;
    bool any_aligned = false;
    while (getline(partition, part, ' ')) {
      cluster_id = std::stoul(part);
      alignment[cluster_id] = true;
      any_aligned = true;
    }
    if (any_aligned) {
      if (reads->find(read_id) == reads->end()) {
	insert_read(alignment, cluster_id, read_id, reads);
      } else {
	switch (mode) {
	 case m_unpaired :
	  insert_read(alignment, cluster_id, *max_read_id + 1 + read_id, reads);
	  break;
	 default :
	  read_paired(mode, alignment, read_id, n_refs, &any_aligned, reads);
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
    read_to_ref->insert(std::make_pair(ec.first, std::vector<uint16_t>()));
    for (size_t i = 0; i < ec.second.pseudoalignment.size(); ++i) {
      if (ec.second.pseudoalignment.at(i)) {
	read_to_ref->at(ec.first).push_back(i);
      }
    }
  }
}

// KAlignment ReadAlignments(const Mode &mode, const uint32_t n_refs, std::istream* strand_1, std::istream* strand_2) {
KAlignment ReadAlignments(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*>* strands) {
  std::unordered_map<uint32_t, ec_info> ecs_by_id;
  uint32_t max_read_id = 0;
  for (size_t i = 0; i < strands->size(); ++i) {
    read_alignments(mode, n_refs, strands->at(i), &ecs_by_id, &max_read_id);
  }

  KAlignment alignments;
  ReadToRef(ecs_by_id, &alignments.read_to_ref);
  CompressAlignments(ecs_by_id, &alignments.ecs);

  // Update the metadata
  alignments.n_targets = n_refs;
  alignments.n_processed = (mode == m_unpaired ? 2*max_read_id + 1 : max_read_id + 1);
  alignments.n_pseudoaligned = ecs_by_id.size();
  alignments.n_unique = alignments.ecs.size();

  return alignments;
}
