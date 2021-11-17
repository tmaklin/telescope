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

#include "read_kallisto_alignments.hpp"

#include <vector>
#include <string>
#include <sstream>

#include "common.hpp"

namespace telescope {
void ReadEquivalenceClasses(const std::vector<uint32_t> &ec_ids, const uint32_t n_refs, std::istream &stream, std::vector<std::vector<bool>> *ec_configs) {
  std::string line;
  
  uint32_t current_ec_pos = 0;
  while(getline(stream, line)) {
    std::string part;
    std::stringstream partition(line);
    getline(partition, part, '\t');
    uint32_t ec_id = std::stoul(part);
    if (ec_id == ec_ids[current_ec_pos]) {
      ec_configs->emplace_back(std::vector<bool>(n_refs, false));
      getline(partition, part, '\t');
      std::string aln;
      std::stringstream alns(part);
      while(getline(alns, aln, ',')) {
  	ec_configs->back()[std::stoul(aln)] = true;
      }
      ++current_ec_pos;
    }
  }
}

void ReadAlignmentCounts(std::istream &stream, std::vector<uint32_t> *ec_ids, std::vector<uint32_t> *ec_counts) {
  std::string line;
  while(getline(stream, line)) {
    std::string part;
    std::stringstream partition(line);
    getline(partition, part, '\t');
    uint32_t ec_id = std::stoul(part);
    getline(partition, part, '\t');
    uint32_t ec_count = std::stoul(part);
    if (ec_count > 0) {
      ec_ids->push_back(ec_id);
      ec_counts->push_back(ec_count);
    }
  }
}

namespace read {
void Kallisto(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, CompressedAlignment *aln) {
  std::vector<uint32_t> ec_ids;
  ReadAlignmentCounts(tsv_file, &ec_ids, &aln->ec_counts);
  ReadEquivalenceClasses(ec_ids, n_refs, ec_file, &aln->ec_configs);
}

void KallistoEcIds(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, KallistoAlignment *kaln) {
  ReadAlignmentCounts(tsv_file, &kaln->ec_ids, &kaln->ec_counts);
  ReadEquivalenceClasses(kaln->ec_ids, n_refs, ec_file, &kaln->ec_configs);
}
}
}
