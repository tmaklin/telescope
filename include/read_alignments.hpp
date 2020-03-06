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

#ifndef TELESCOPE_READ_ALIGNMENTS_HPP
#define TELESCOPE_READ_ALIGNMENTS_HPP

#include <fstream>

#include "common.hpp"
#include "read_kallisto_alignments.hpp"
#include "read_themisto_alignments.hpp"

inline void ReadThemisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &strands, CompressedAlignment* aln) { ReadThemistoFiles(mode, n_refs, strands, aln); }
inline void ReadThemisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, ThemistoAlignment *aln) { ReadThemistoFiles(mode, n_refs, streams, aln); }
inline void ReadKallisto(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, CompressedAlignment *aln) { ReadKallistoFiles(n_refs, ec_file, tsv_file, aln); }
inline void ReadKallisto(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, KallistoAlignment *aln) { ReadKallistoFiles(n_refs, ec_file, tsv_file, aln); }

#endif
