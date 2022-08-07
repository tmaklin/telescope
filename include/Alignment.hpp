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

#ifndef TELESCOPE_ALIGNMENT_HPP
#define TELESCOPE_ALIGNMENT_HPP

#include <vector>
#include <cstddef>

#include "bm.h"

template<> struct std::hash<bm::bvector<>> {
    std::size_t operator()(bm::bvector<> const& vec) const noexcept {
	std::vector<bool> cnv(vec.size());
	for (size_t i = 0; i < vec.size(); ++i) {
	    cnv[i] = vec[i];
	}
	return std::hash<std::vector<bool>>{}(cnv);
    }
};

namespace telescope {
struct Alignment {
  std::vector<uint32_t> ec_counts;

  uint32_t n_processed;

  uint32_t size() const { return ec_counts.size(); }
};

struct CompressedAlignment : public Alignment{
    std::vector<bm::bvector<>> ec_configs;

  uint32_t n_targets() const { return this->ec_configs.at(0).size(); }
};

struct GroupedAlignment : public Alignment{
  std::vector<std::vector<uint16_t>> ec_group_counts;
};

struct ThemistoAlignment : public CompressedAlignment{
  std::vector<std::vector<uint32_t>> aligned_reads;
};
}

#endif
