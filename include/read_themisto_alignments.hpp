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

#ifndef TELESCOPE_READ_THEMISTO_ALIGNMENTS_HPP
#define TELESCOPE_READ_THEMISTO_ALIGNMENTS_HPP

#include <cstddef>
#include <vector>
#include <fstream>

#include "Alignment.hpp"
#include "KallistoAlignment.hpp"

namespace telescope {
namespace read {
// Read equivalence classes
ThemistoAlignment Themisto(const bm::set_operation &merge_op, const size_t n_refs, std::vector<std::istream*> &streams);

// Read the plain alignment (n_reads x n_targets)
ThemistoAlignment ThemistoPlain(const bm::set_operation &merge_op, const size_t n_refs, std::vector<std::istream*> &streams);

// Group the targets by some `group_indicators` and read corresponding equivalence classes
GroupedAlignment ThemistoGrouped(const bm::set_operation &merge_op, const size_t n_refs, const size_t n_groups, const std::vector<uint32_t> &group_indicators, std::vector<std::istream*> &streams);

// Read Themisto format alignment and convert it to Kallisto format
KallistoAlignment ThemistoToKallisto(const bm::set_operation &merge_op, const size_t n_refs, std::vector<std::istream*> &streams);

}
}

#endif
