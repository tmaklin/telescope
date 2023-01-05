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
void ReadEquivalenceClasses(std::istream &stream, KallistoAlignment *aln) {
  std::string line;
  bm::bvector<>::bulk_insert_iterator it = aln->get_iterator();

  uint32_t current_ec_pos = 0;
  while(getline(stream, line)) {
    std::string part;
    std::stringstream partition(line);
    getline(partition, part, '\t');
    uint32_t ec_id = std::stoul(part);
    if (ec_id == aln->get_ec_id(current_ec_pos)) {
      getline(partition, part, '\t');
      std::string line;
      std::stringstream alns(part);
      while(getline(alns, line, ',')) {
	it = current_ec_pos*aln->n_targets() + std::stoul(line);
      }
      ++current_ec_pos;
    }
  }
}

void ReadAlignmentCounts(std::istream &stream, KallistoAlignment *kaln) {
  std::string line;
  while(getline(stream, line)) {
    std::string part;
    std::stringstream partition(line);
    getline(partition, part, '\t');
    uint32_t ec_id = std::stoul(part);
    getline(partition, part, '\t');
    uint32_t ec_count = std::stoul(part);
    kaln->insert(ec_id, ec_count);
  }
}

namespace read {
void Kallisto(std::istream &ec_file, std::istream &tsv_file, KallistoAlignment *aln) {
  ReadAlignmentCounts(tsv_file, aln);
  ReadEquivalenceClasses(ec_file, aln);
}
}
}
