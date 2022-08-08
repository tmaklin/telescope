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

#include "telescope.hpp"

#include <string>
#include <sstream>
#include <unordered_map>

#include "bm.h"

namespace telescope {
uint32_t ReadAlignments(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, bm::bvector<> *all_ecs) {
  // Returns the number of reads processed
  uint32_t n_reads = 0;
  uint8_t n_streams = streams.size();
  std::vector<std::string> lines(n_streams);

  std::vector<bm::bvector<>> allbits(n_streams);
  std::vector<bm::bvector<>::bulk_insert_iterator> its;
  for (uint8_t i = 0; i < n_streams; ++i) {
    its.emplace_back(bm::bvector<>::bulk_insert_iterator(allbits[i]));
  }

  while (std::getline(*streams[0], lines[0])) {
    for (uint8_t i = 1; i < n_streams; ++i) {
      std::getline(*streams[i], lines[i]);
    }

    for (uint8_t i = 0; i < n_streams; ++i) {
      std::string part;
      std::stringstream partition(lines[i]);
      getline(partition, part, ' ');
      while (getline(partition, part, ' ')) {
	  its[i] = n_reads*n_refs + std::stoul(part);
      }
    }
    ++n_reads;
  }

  (*all_ecs) = bm::bvector<>(allbits[0].size(), bm::BM_GAP);
  if (mode == m_intersection) {
    all_ecs->set();
  }
  for (uint8_t i = 0; i < n_streams; ++i) {
    if (mode == m_intersection) {
      all_ecs->bit_and(allbits[i], bm::bvector<>::opt_compress);
    } else {
      all_ecs->bit_or(allbits[i], bm::bvector<>::opt_compress);
    }
  }

  return n_reads;
}

uint32_t ReadGroupedAlignments(const Mode &mode, const std::vector<uint16_t> &group_indicators, const uint32_t n_refs, const uint16_t n_groups, std::vector<std::istream*> &streams, std::vector<std::vector<uint16_t>> *ec_group_counts, std::vector<uint32_t> *ec_counts) {
  // Returns the number of reads processed
  uint32_t n_reads = 0;
  uint8_t n_streams = streams.size();
  std::vector<std::string> lines(n_streams);

  std::unordered_map<bm::bvector<>, uint32_t> ec_to_pos;

  uint32_t ec_pos;
  while (std::getline(*streams[0], lines[0])) {
    for (uint8_t i = 1; i < n_streams; ++i) {
      std::getline(*streams[i], lines[i]);
    }

    bm::bvector<> current_ec(n_refs);
    if (mode == m_intersection)
      current_ec.set_range(0, n_refs, true);

    std::vector<bm::bvector<>> proposed_ecs(n_streams, bm::bvector<>(n_refs));
    for (uint8_t i = 0; i < n_streams; ++i) {
      std::string part;
      std::stringstream partition(lines[i]);
      getline(partition, part, ' ');
      while (getline(partition, part, ' ')) {
	proposed_ecs[i][std::stoul(part)] = true;
      }
    }
    for (uint8_t i = 0; i < n_streams; ++i) {
      if (mode == m_intersection)
	current_ec &= proposed_ecs[i];
      else
	current_ec |= proposed_ecs[i];
    }
    if (current_ec.any()) {
      if (ec_to_pos.find(current_ec) == ec_to_pos.end()) {
        ec_group_counts->emplace_back(std::vector<uint16_t>(n_groups, 0));
	ec_counts->push_back(0);
	ec_to_pos.insert(std::make_pair(current_ec, ec_pos));
	++ec_pos;
	std::vector<uint16_t> *current_read = &ec_group_counts->back();
	for (uint32_t j = 0; j  < n_refs; ++j) {
	  (*current_read)[group_indicators[j]] += current_ec[j];
	}
      }
      (*ec_counts)[ec_to_pos[current_ec]] += 1;
    }
    ++n_reads;
  }
  return n_reads;
}

uint32_t ReadAndAssign(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, std::vector<bm::bvector<>> *ec_configs, std::vector<uint32_t> *aligned_reads) {
  // Returns the number of reads processed
  uint32_t n_reads = 0;
  uint8_t n_streams = streams.size();
  std::vector<std::string> lines(n_streams);

  while (std::getline(*streams[0], lines[0])) {
    for (uint8_t i = 1; i < n_streams; ++i) {
      std::getline(*streams[i], lines[i]);
    }
    bm::bvector<> current_ec(n_refs);
    if (mode == m_intersection)
      current_ec.set_range(0, n_refs, true);

    uint32_t read_id;
    std::vector<bm::bvector<>> proposed_ecs(n_streams, bm::bvector<>(n_refs));
    for (uint8_t i = 0; i < n_streams; ++i) {
      std::string part;
      std::stringstream partition(lines[i]);
      getline(partition, part, ' ');
      read_id = std::stoul(part);
      while (getline(partition, part, ' ')) {
	proposed_ecs[i][std::stoul(part)] = true;
      }
    }
    for (uint8_t i = 0; i < n_streams; ++i) {
      if (mode == m_intersection)
	current_ec &= proposed_ecs[i];
      else
	current_ec |= proposed_ecs[i];
    }
    if (current_ec.any()) {
      ec_configs->emplace_back(current_ec);
      aligned_reads->emplace_back(read_id);
    }
    ++n_reads;
  }
  return n_reads;
}

std::vector<bm::bvector<>> CompressAlignment(const bm::bvector<> &all_ecs, const size_t n_refs, const size_t num_alns, std::vector<uint32_t> *ec_counts) {
  // Compress the alignment into equivalence classes
  std::vector<bm::bvector<>> compressed_ec_configs;
  std::unordered_map<bm::bvector<>, uint32_t> ec_to_pos;
  uint32_t ec_pos = 0;
  std::unordered_map<bm::bvector<>, uint32_t>::iterator it;
  for (unsigned i = 0; i < num_alns; ++i) {
    bm::bvector<> current_ec(n_refs);
    bool any_aligned = false;
    for (size_t j = 0; j < n_refs; ++j) {
      current_ec[j] = all_ecs[i*n_refs + j];
      any_aligned |= current_ec[j];
    }
    if (any_aligned) {
      it = ec_to_pos.find(current_ec);
      if (it == ec_to_pos.end()) {
	compressed_ec_configs.emplace_back(bm::bvector<>(current_ec, bm::finalization::READONLY));
	ec_counts->emplace_back(0);
	it = ec_to_pos.insert(std::make_pair(current_ec, ec_pos)).first; // return iterator to inserted element
	++ec_pos;
      }
      (*ec_counts)[it->second] += 1;
    }
  }
  return compressed_ec_configs;
}

std::vector<bm::bvector<>> CompressAndAssign(const std::vector<uint32_t> &read_ids, const std::vector<bm::bvector<>> &ec_configs, std::vector<uint32_t> *ec_counts, std::vector<std::vector<uint32_t>> *aln_reads) {
  // Compress the alignment into equivalence classes and store the read assignments to equivalence classes
  std::vector<bm::bvector<>> compressed_ec_configs;
  uint32_t num_alns = ec_configs.size();
  std::unordered_map<bm::bvector<>, uint32_t> ec_to_pos;
  uint32_t ec_pos = 0;
  for (unsigned i = 0; i < num_alns; ++i) {
    if (ec_to_pos.find(ec_configs[i]) == ec_to_pos.end()) {
      compressed_ec_configs.push_back(ec_configs[i]);
      ec_counts->push_back(0);
      ec_to_pos.insert(std::make_pair(ec_configs[i], ec_pos));
      aln_reads->emplace_back(std::vector<uint32_t>());
      ++ec_pos;
    }
    uint32_t current_pos = ec_to_pos[ec_configs[i]];
    (*ec_counts)[current_pos] += 1;
    (*aln_reads)[current_pos].emplace_back(read_ids[i]);
  }
  return compressed_ec_configs;
}

namespace read {
void Themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, CompressedAlignment *aln) {
  // Read in only the ec_configs
  bm::bvector<> all_ecs;
  aln->n_processed = ReadAlignments(mode, n_refs, streams, &all_ecs);
  aln->ec_configs = CompressAlignment(all_ecs, n_refs, aln->n_processed, &aln->ec_counts);
}

void ThemistoGrouped(const Mode &mode, const std::vector<uint16_t> &group_indicators, const uint32_t n_refs, const uint16_t n_groups, std::vector<std::istream*> &streams, GroupedAlignment *aln) {
  // Read in group counts
  aln->n_processed = ReadGroupedAlignments(mode, group_indicators, n_refs, n_groups, streams, &aln->ec_group_counts, &aln->ec_counts);
}

void ThemistoAlignedReads(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, ThemistoAlignment *taln) {
  // Read in the ec_configs and assign reads to equivalence classes
  std::vector<uint32_t> aligned_reads_ids;
  taln->n_processed = ReadAndAssign(mode, n_refs, streams, &taln->ec_configs, &aligned_reads_ids);
  taln->ec_configs = CompressAndAssign(aligned_reads_ids, taln->ec_configs, &taln->ec_counts, &taln->aligned_reads);
}

void ThemistoToKallisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, KallistoAlignment *aln) {
  // Read in the ec_configs and fill the ec_ids vector
    //  aln->n_processed = ReadAlignments(mode, n_refs, streams, &aln->ec_configs);
    //aln->ec_configs = CompressAlignment(aln->ec_configs, &aln->ec_counts);
  for (uint32_t i = 0; i < aln->size(); ++i) {
    aln->ec_ids.emplace_back(i);
  }
}
}
}
