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

#include <vector>
#include <string>
#include <sstream>

#include "bm64.h"

namespace telescope {
void ReadEquivalenceClasses(const std::vector<uint32_t> &ec_ids, const uint32_t n_refs, std::istream &stream, bm::bvector<> *ec_configs) {
  std::string line;
  ec_configs->set_new_blocks_strat(bm::BM_GAP);
  bm::bvector<>::bulk_insert_iterator it(*ec_configs);

  uint32_t current_ec_pos = 0;
  while(getline(stream, line)) {
    std::string part;
    std::stringstream partition(line);
    getline(partition, part, '\t');
    uint32_t ec_id = std::stoul(part);
    if (ec_id == ec_ids[current_ec_pos]) {
      getline(partition, part, '\t');
      std::string aln;
      std::stringstream alns(part);
      while(getline(alns, aln, ',')) {
	it = current_ec_pos*n_refs + std::stoul(aln);
      }
      ++current_ec_pos;
    }
  }
  if (ec_configs->size() != current_ec_pos*n_refs) {
    ec_configs->resize(current_ec_pos*n_refs); // add trailing zeros
  }
  ec_configs->optimize();
  ec_configs->freeze();
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
  ReadEquivalenceClasses(ec_ids, n_refs, ec_file, aln->get());
}

void KallistoEcIds(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, KallistoAlignment *kaln) {
  ReadAlignmentCounts(tsv_file, &kaln->ec_ids, &kaln->ec_counts);
  ReadEquivalenceClasses(kaln->ec_ids, n_refs, ec_file, kaln->get());
}
}
}
