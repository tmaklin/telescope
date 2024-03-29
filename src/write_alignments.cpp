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

#include <string>

namespace telescope {
namespace write {
void ThemistoToKallisto(const ThemistoAlignment &aln, std::ostream* ec_file, std::ostream* tsv_file) {
  // telescope::write::ThemistoToKallisto
  //
  // Writes the alignment contained in `aln` into files matching the Kallisto format.
  //
  // Input:
  //   `aln`: The pseudoalignment to write.
  //   `ec_file`: Pointer to the file that will store the equivalence class configurations.
  //   `tsv_file`: Pointer to the file that will contain the observation counts of each equivalence class.
  //
  uint32_t ec_id = 0;
  for (uint32_t i = 0; i < aln.n_ecs(); ++i) {
    std::string aligneds("");
    for (uint32_t j = 0; j < aln.n_targets(); ++j) {
      if (aln(i, j)) {
	aligneds += std::to_string(j);
	aligneds += ',';
      }
    }
    aligneds.pop_back();
    *ec_file << ec_id << '\t' << aligneds << '\n';
    *tsv_file << ec_id << '\t' << aln.reads_in_ec(i) << '\n';
    ++ec_id;
  }
  ec_file->flush();
  tsv_file->flush();
}

void ThemistoReadAssignments(const ThemistoAlignment &aln, std::ostream* out) {
  // telescope::write::ThemistoReadAssignments
  //
  // Writes the alignment of each read against the reference sequences.
  //
  // Input:
  //   `aln`: The pseudoalignment to write.
  //   `out`: Pointer to the output file stream.
  //
  for (uint32_t i = 0; i < aln.n_ecs(); ++i) {
    for (uint32_t j = 0; j < aln.reads_assigned_to_ec(i).size(); ++j) {
      *out << aln.reads_assigned_to_ec(i)[j] << ' ';
      std::string aligned_to("");
      for (uint32_t k = 0; k < aln.n_targets(); ++k) {
	if (aln(i, k)) {
	  aligned_to += std::to_string(k);
	  aligned_to += ' ';
	}
      }
      aligned_to.pop_back();
      *out << aligned_to << '\n';
    }
  }
  out->flush();
}

void KallistoInfoFile(const KallistoRunInfo &run_info, const uint8_t indent_len, std::ostream *out) {
  // telescope::write::KallistoInfoFile
  //
  // Writes the Kallisto run_info.json file from the `run_info` object.
  //
  // Input:
  //   `run_info`: The KallistoRunInfo object to write.
  //   `indent_len`: Indent length in the written json file.
  //   `out`: Pointer to a file stream to write into.
  //
  std::string indent;
  for (uint8_t i = 0; i < indent_len; ++i) {
    indent += " ";
  }
  *out << "{" << '\n';
  *out << indent << "\"n_targets\": " << run_info.n_targets << ',' << '\n';
  *out << indent << "\"n_bootstraps\": " << run_info.n_bootstraps << ',' << '\n';
  *out << indent << "\"n_processed\": " << run_info.n_processed << ',' << '\n';
  *out << indent << "\"n_pseudoaligned\": " << run_info.n_pseudoaligned << ',' << '\n';
  *out << indent << "\"n_unique\": " << run_info.n_unique << ',' << '\n';
  *out << indent << "\"p_pseudoaligned\": " << run_info.p_pseudoaligned << ',' << '\n';
  *out << indent << "\"p_unique\": " << run_info.p_unique << ',' << '\n';
  *out << indent << "\"kallisto_version\": \"" << run_info.kallisto_version << '\"' << ',' << '\n';
  *out << indent << "\"index_version\": " << run_info.index_version << ',' << '\n';
  *out << indent << "\"start_time\": \"" << run_info.start_time << '\"' << ',' << '\n';
  *out << indent << "\"call\": \"" << run_info.call << '\"' << '\n';
  *out << "}" << '\n';
  out->flush();
}
} // namespace write
} // namespace telescope
