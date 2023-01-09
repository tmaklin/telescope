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
// Functions for reading pseudoalignment files into telescope::Alignment objects.

// telescope::read::Themisto
//
// Read in a Themisto pseudoalignment and collapse it into
// equivalence classes. Reads that align to exactly same reference
// sequences are assigned to the same equivalence class.
//
// Input:
//   `merge_op`: bm::set_OR for union or bm::set_AND for intersection of multiple alignmnet files
//   `n_refs`: number of pseudoalignment targets (reference
//                sequences). It's not possible to infer this from the plaintext Themisto
//                file format so has to be provided separately. If the file is in the
//                compact format will check that the numbers match.
//   `streams`: vector of pointers to the istreams opened on the pseudoalignment files.
// Output:
//   `aln`: The pseudoalignment as a telescope::ThemistoAlignment object.
ThemistoAlignment Themisto(const bm::set_operation &merge_op, const size_t n_refs, std::vector<std::istream*> &streams);

// telescope::read::ThemistoPlain
//
// Read in a Themisto pseudoalignment in the plain format
// i. e. without collapsing it into equivalence classes.
//
// Input:
//   `merge_op`: bm::set_OR for union or bm::set_AND for intersection of multiple alignmnet files
//   `n_refs`: number of pseudoalignment targets (reference
//                sequences). It's not possible to infer this from the plaintext Themisto
//                file format so has to be provided separately. If the file is in the
//                compact format will check that the numbers match.
//   `streams`: vector of pointers to the istreams opened on the pseudoalignment files.
// Output:
//   `aln`: The pseudoalignment as a telescope::ThemistoAlignment object.
ThemistoAlignment ThemistoPlain(const bm::set_operation &merge_op, const size_t n_refs, std::vector<std::istream*> &streams);

// telescope::read::ThemistoGrouped
//
// Read in a Themisto pseudoalignment and collapse it into
// equivalence classes. Reads that align to the same number of
// reference sequences in each reference group (defined by
// `group_indicators`) are assigned to the same equivalence class.
//
// Input:
//   `merge_op`: bm::set_OR for union or bm::set_AND for intersection of multiple alignmnet files
//   `n_refs`: number of pseudoalignment targets (reference
//                sequences). It's not possible to infer this from the plaintext Themisto
//                file format so has to be provided separately. If the file is in the
//                compact format will check that the numbers match.
//   `group_indicators`: Vector assigning each reference sequence to a reference group. The group
//                       of the n:th sequence is the value at the (n - 1):th position in the vector.
//   `streams`: vector of pointers to the istreams opened on the pseudoalignment files.
// Output:
//   `aln`: The pseudoalignment as a telescope::GroupedAlignment object.
GroupedAlignment ThemistoGrouped(const bm::set_operation &merge_op, const size_t n_refs, const std::vector<uint32_t> &group_indicators, std::vector<std::istream*> &streams);

// telescope::read::ThemistoToKallisto
//
// Read in a Themisto pseudoalignment and convert it into a Kallisto pseudoalignment.
//
// Input:
//   `merge_op`: bm::set_OR for union or bm::set_AND for intersection of multiple alignmnet files
//   `n_refs`: number of pseudoalignment targets (reference
//                sequences). It's not possible to infer this from the plaintext Themisto
//                file format so has to be provided separately. If the file is in the
//                compact format will check that the numbers match.
//   `streams`: vector of pointers to the istreams opened on the pseudoalignment files.
// Output:
//   `aln`: The pseudoalignment as a telescope::KallistoAlignment object.
KallistoAlignment ThemistoToKallisto(const bm::set_operation &merge_op, const size_t n_refs, std::vector<std::istream*> &streams);

}
}

#endif
