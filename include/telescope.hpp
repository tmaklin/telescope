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

#ifndef TELESCOPE_TELESCOPE_HPP
#define TELESCOPE_TELESCOPE_HPP

#include "write_alignments.hpp"
#include "common.hpp"

namespace telescope {
namespace read {
// Themisto input
void Themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, CompressedAlignment *aln);
void ThemistoGrouped(const Mode &mode, const std::vector<uint16_t> &group_indicators, const uint32_t n_refs, const uint16_t n_groups, std::vector<std::istream*> &streams, GroupedAlignment *aln);
void ThemistoAlignedReads(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, ThemistoAlignment *aln);
void ThemistoToKallisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, KallistoAlignment *aln);

// Kallisto input
// Read pseudoalignments and counts
void Kallisto(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, CompressedAlignment *aln);
// Also read the equivalence class IDs
void KallistoEcIds(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, KallistoAlignment *aln);
}
}

#endif
