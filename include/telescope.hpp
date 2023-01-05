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

#include <fstream>
#include <cstddef>

#include "read_themisto_alignments.hpp"
#include "Alignment.hpp"
#include "KallistoAlignment.hpp"

namespace telescope {
namespace read {
  // Check `read_themisto_alignments.hpp`
}

namespace write {
void ThemistoToKallisto(const ThemistoAlignment &aln, std::ostream* ec_file, std::ostream* tsv_file);
void ThemistoReadAssignments(const ThemistoAlignment &aln, std::ostream* out);
void KallistoInfoFile(const KallistoRunInfo &run_info, const uint8_t indent_len, std::ostream *out);
}
}

#endif
