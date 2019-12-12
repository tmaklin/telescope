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

#ifndef TELESCOPE_WRITE_ALIGNMENTS_HPP
#define TELESCOPE_WRITE_ALIGNMENTS_HPP

#include <unordered_map>
#include <vector>
#include <fstream>

#include "common.hpp"

void WriteAlignments(const std::unordered_map<std::vector<bool>, ec_info> &ecs, std::ostream* ec_file, std::ostream* tsv_file);
void WriteReadToRef(const std::unordered_map<uint32_t, std::vector<uint16_t>> &read_to_ref, std::ostream* out);
void WriteRunInfo(const KAlignment &alignment, std::ostream* out, uint8_t indent_len = 4);

#endif
