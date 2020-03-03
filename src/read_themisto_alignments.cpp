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

void ReadAlignment(const uint32_t n_refs, std::istream &stream, std::vector<std::vector<bool>> *ec_configs) {
  std::string line;
  while (getline(stream, line)) {
    ec_configs->emplace_back(std::vector<bool>(n_refs, false));
    std::string part;
    std::stringstream partition(line);
    getline(partition, part, ' ');
    while (getline(partition, part, ' ')) {
      uint16_t cluster_id = std::stoul(part);
      ec_configs->back()[cluster_id] = true;
    }
  }
}

void ReadIntersection(const uint32_t n_refs, std::vector<std::istream*> &streams, std::vector<std::vector<bool>> *ec_configs) {
  ReadAlignment(n_refs, *streams[0], ec_configs);
  uint32_t num_ecs = ec_configs->size();
  for (uint8_t i = 1; i < streams.size(); ++i) {
    std::vector<std::vector<bool>> pair_alns;
    ReadAlignment(n_refs, *streams[i], &pair_alns);
    for (uint32_t j = 0; j < num_ecs; ++j) {
      for (uint32_t k = 0; k < n_refs; ++k) {
	(*ec_configs)[j][k] = (*ec_configs)[j][k] && pair_alns[j][k];
      }
    }
  }
}

void ReadUnion(const uint32_t n_refs, std::vector<std::istream*> &streams, std::vector<std::vector<bool>> *ec_configs) {
  ReadAlignment(n_refs, *streams[0], ec_configs);
  uint32_t num_ecs = ec_configs->size();
  for (uint8_t i = 1; i < streams.size(); ++i) {
    std::vector<std::vector<bool>> pair_alns;
    ReadAlignment(n_refs, *streams[i], &pair_alns);
    for (uint32_t j = 0; j < num_ecs; ++j) {
      for (uint32_t k = 0; k < n_refs; ++k) {
	(*ec_configs)[j][k] = (*ec_configs)[j][k] || pair_alns[j][k];
      }
    }
  }
}

void CompressAlignment(std::vector<std::vector<bool>> &ec_configs, std::vector<uint32_t> *ec_counts, std::vector<std::vector<bool>> *compressed_ec_configs) {
  uint32_t num_alns = ec_configs.size();
  std::unordered_map<std::vector<bool>, uint32_t> ec_to_pos;
  uint32_t ec_pos = 0;
  for (unsigned i = 0; i < num_alns; ++i) {
    if (ec_to_pos.find(ec_configs[i]) == ec_to_pos.end()) {
      compressed_ec_configs->push_back(ec_configs[i]);
      ec_counts->push_back(1);
      ec_to_pos.insert(std::make_pair(ec_configs[i], ec_pos));
      ++ec_pos;
    } else {
      (*ec_counts)[ec_to_pos[ec_configs[i]]] += 1;
    }
  }
}

void ReadThemistoFiles(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, std::vector<uint32_t> *ec_counts, std::vector<std::vector<bool>> *compressed_ec_configs) {
  std::vector<std::vector<bool>> ec_configs;
  if (mode == m_intersection) {
    ReadIntersection(n_refs, streams, &ec_configs);
  } else if (mode == m_union) {
    ReadUnion(n_refs, streams, &ec_configs);
  }
  CompressAlignment(ec_configs, ec_counts, compressed_ec_configs);
}
