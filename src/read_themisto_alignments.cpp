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

#include "read_themisto_alignments.hpp"

#include <string>
#include <sstream>
#include <unordered_map>

#include "common.hpp"

namespace telescope {
uint32_t ReadAlignments(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, std::vector<std::vector<bool>> *ec_configs) {
  // Returns the number of reads processed
  uint32_t n_reads = 0;
  uint8_t n_streams = streams.size();
  std::vector<std::string> lines(n_streams);

  while (std::getline(*streams[0], lines[0])) {
    for (uint8_t i = 1; i < n_streams; ++i) {
      std::getline(*streams[i], lines[i]);
    }
    ec_configs->emplace_back(std::vector<bool>(n_refs, mode == m_intersection));
    std::vector<bool> *current_ec = &ec_configs->back();

    std::vector<std::vector<bool>> proposed_ecs(n_streams, std::vector<bool>(n_refs, false));
    for (uint8_t i = 0; i < n_streams; ++i) {
      std::string part;
      std::stringstream partition(lines[i]);
      getline(partition, part, ' ');
      while (getline(partition, part, ' ')) {
	proposed_ecs[i][std::stoul(part)] = true;
      }
    }
    for (uint8_t i = 0; i < n_streams; ++i) {
      for (uint32_t j = 0; j < n_refs; ++j) {
	(*current_ec)[j] = (mode == m_intersection ? ((*current_ec)[j] && proposed_ecs[i][j]) : ((*current_ec)[j] || proposed_ecs[i][j]));
      }
    }
    bool any_aligned = false;
    for (uint32_t j = 0; j < n_refs; ++j) {
      any_aligned = (*current_ec)[j] || any_aligned;
    }
    if (!any_aligned) {
      ec_configs->pop_back();
    }
    ++n_reads;
  }
  return n_reads;
}

uint32_t ReadGroupedAlignments(const Mode &mode, const std::vector<uint16_t> &group_indicators, const uint32_t n_refs, const uint16_t n_groups, std::vector<std::istream*> &streams, std::vector<std::vector<uint16_t>> *ec_group_counts, std::vector<uint32_t> *ec_counts) {
  // Returns the number of reads processed
  uint32_t n_reads = 0;
  uint8_t n_streams = streams.size();
  std::vector<std::string> lines(n_streams);

  std::unordered_map<std::vector<bool>, uint32_t> ec_to_pos;

  uint32_t ec_pos;
  while (std::getline(*streams[0], lines[0])) {
    for (uint8_t i = 1; i < n_streams; ++i) {
      std::getline(*streams[i], lines[i]);
    }

    std::vector<bool> current_ec(n_refs, mode == m_intersection);
    std::vector<std::vector<bool>> proposed_ecs(n_streams, std::vector<bool>(n_refs, false));
    for (uint8_t i = 0; i < n_streams; ++i) {
      std::string part;
      std::stringstream partition(lines[i]);
      getline(partition, part, ' ');
      while (getline(partition, part, ' ')) {
	proposed_ecs[i][std::stoul(part)] = true;
      }
    }
    for (uint8_t i = 0; i < n_streams; ++i) {
      for (uint32_t j = 0; j < n_refs; ++j) {
	current_ec[j] = (mode == m_intersection ? (current_ec[j] && proposed_ecs[i][j]) : (current_ec[j] || proposed_ecs[i][j]));
      }
    }
    bool any_aligned = false;
    for (uint32_t j = 0; j < n_refs; ++j) {
      any_aligned = current_ec[j] || any_aligned;
    }
    if (any_aligned) {
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

uint32_t ReadAndAssign(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, std::vector<std::vector<bool>> *ec_configs, std::vector<uint32_t> *aligned_reads) {
  // Returns the number of reads processed
  uint32_t n_reads = 0;
  uint8_t n_streams = streams.size();
  std::vector<std::string> lines(n_streams);

  while (std::getline(*streams[0], lines[0])) {
    for (uint8_t i = 1; i < n_streams; ++i) {
      std::getline(*streams[i], lines[i]);
    }
    ec_configs->emplace_back(std::vector<bool>(n_refs, mode == m_intersection));
    std::vector<bool> *current_ec = &ec_configs->back();

    uint32_t read_id;
    std::vector<std::vector<bool>> proposed_ecs(n_streams, std::vector<bool>(n_refs, false));
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
      for (uint32_t j = 0; j < n_refs; ++j) {
	(*current_ec)[j] = (mode == m_intersection ? ((*current_ec)[j] && proposed_ecs[i][j]) : ((*current_ec)[j] || proposed_ecs[i][j]));
      }
    }
    bool any_aligned = false;
    for (uint32_t j = 0; j < n_refs; ++j) {
      any_aligned = (*current_ec)[j] || any_aligned;
    }
    if (!any_aligned) {
      ec_configs->pop_back();
    } else {
      aligned_reads->emplace_back(read_id);
    }
    ++n_reads;
  }
  return n_reads;
}

std::vector<std::vector<bool>> CompressAlignment(const std::vector<std::vector<bool>> &ec_configs, std::vector<uint32_t> *ec_counts) {
  // Compress the alignment into equivalence classes
  std::vector<std::vector<bool>> compressed_ec_configs;
  uint32_t num_alns = ec_configs.size();
  std::unordered_map<std::vector<bool>, uint32_t> ec_to_pos;
  uint32_t ec_pos = 0;
  for (unsigned i = 0; i < num_alns; ++i) {
    if (ec_to_pos.find(ec_configs[i]) == ec_to_pos.end()) {
      compressed_ec_configs.push_back(ec_configs[i]);
      ec_counts->push_back(0);
      ec_to_pos.insert(std::make_pair(ec_configs[i], ec_pos));
      ++ec_pos;
    }
    (*ec_counts)[ec_to_pos[ec_configs[i]]] += 1;
  }
  return compressed_ec_configs;
}

std::vector<std::vector<bool>> CompressAndAssign(const std::vector<uint32_t> &read_ids, const std::vector<std::vector<bool>> &ec_configs, std::vector<uint32_t> *ec_counts, std::vector<std::vector<uint32_t>> *aln_reads) {
  // Compress the alignment into equivalence classes and store the read assignments to equivalence classes
  std::vector<std::vector<bool>> compressed_ec_configs;
  uint32_t num_alns = ec_configs.size();
  std::unordered_map<std::vector<bool>, uint32_t> ec_to_pos;
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
  aln->n_processed = ReadAlignments(mode, n_refs, streams, &aln->ec_configs);
  aln->ec_configs = CompressAlignment(aln->ec_configs, &aln->ec_counts);
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
  aln->n_processed = ReadAlignments(mode, n_refs, streams, &aln->ec_configs);
  aln->ec_configs = CompressAlignment(aln->ec_configs, &aln->ec_counts);
  for (uint32_t i = 0; i < aln->size(); ++i) {
    aln->ec_ids.emplace_back(i);
  }
}
}
}
