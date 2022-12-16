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
//
#include "read_themisto_alignments.hpp"

#include <string>
#include <sstream>
#include <unordered_map>
#include <functional>

#include "bm64.h"

#include "telescope.hpp"

namespace telescope {
void ReadAlignmentFile(std::istream *stream, bm::bvector<> *ec_configs, Alignment *alignment) {
  // telescope::ReadAlignmentFile
  //
  // Reads in a istream pointing to a single themisto pseudoalignment
  // as a contiguous bit vector.
  //
  // Input:
  //   `stream`: pointer to an istream opened on the pseudoalignment file.
  //   `alignment`: pointer to the object that will contain the results.
  //
  bm::bvector<>::bulk_insert_iterator it(*ec_configs);
  std::string line;
  while (std::getline(*stream, line)) {
    if (alignment->parse_from_buffered()) {
      alignment->parse(line, stream, ec_configs);
    } else {
      alignment->parse(line, &it);
    }
  }

  //alignment->add_trailing_zeros(alignment->n_reads(), alignment->n_targets());
  ec_configs->optimize();
}

void ReadPairedAlignments(const Mode &mode, std::vector<std::istream*> &streams, bm::bvector<> *ec_configs, Alignment *alignment) {
  // telescope::ReadPairedAlignments
  //
  // Reads in one or more pseudoalignment files from themisto for paired reads.
  //
  // Input:
  //   `mode`: intersect, union, or unpair the paired alignments (one of m_intersection, m_union, m_unpaired).
  //   `streams`: pointers to istreams opened on the pseudoalignment files.
  //   `alignment`: pointer to an object that will contain the results.
  //
  uint8_t n_streams = streams.size(); // Typically 1 (unpaired reads) or 2 (paired reads).

  for (uint8_t i = 0; i < n_streams; ++i) {
    if (i == 0) {
      // Read the first alignments in-place to the output variable.
      ReadAlignmentFile(streams[i], ec_configs, alignment);
    } else {
      // Read subsequent alignments into a new bitvector.
      // Size is now known since the paired files should have the same numbers of reads.
      CompressedAlignment pair_alignment(alignment->n_targets(), alignment->n_reads()); // Initialize with known final size
      if (alignment->parse_from_buffered()) {
	pair_alignment.set_parse_from_buffered();
      }
      ReadAlignmentFile(streams[i], pair_alignment.get(), &pair_alignment);

      if (alignment->n_reads() != pair_alignment.n_reads()) {
	// Themisto's output from paired-end reads should contain the same amount of reads.
	throw std::runtime_error("Pseudoalignment files have different numbers of pseudoalignments.");
      }
      if (mode == m_intersection) {
	// m_intersection: both reads in a pair should align to be considered a match.
	(*ec_configs) &= *pair_alignment.get();
      } else {
	// m_union or m_unpaired: count alignments regardless of pair's status.
	(*ec_configs) |= *pair_alignment.get();
      }
    }
  }
}

void CompressAlignment(bm::bvector<> &ec_configs, Alignment *full_alignment) {
  // telescope::CompressAlignment
  //
  // Compresses the full pseudoalignment data into equivalence
  // classes, meaning unique pseudoalignment patterns and the numbers
  // of times they were observed.
  //
  // Input:
  //   `full_alignment`: pointer to the alignment object (after running ReadPairedAlignment or ReadAlignmentFile).
  //
  bm::bvector<> compressed_ec_configs;
  compressed_ec_configs.set_new_blocks_strat(bm::BM_GAP); // Store data in compressed format.
  bm::bvector<>::bulk_insert_iterator bv_it(compressed_ec_configs);

  // Need to hash the alignment patterns to count the times they appear.
  std::unordered_map<std::vector<bool>, uint32_t> ec_to_pos;

  size_t ec_id = 0;
  for (size_t i = 0; i < full_alignment->n_reads(); ++i) {
    // Check if the current read aligned against any reference and
    // discard the read if it didn't.
    if (ec_configs.any_range(i*full_alignment->n_targets(), i*full_alignment->n_targets() + full_alignment->n_targets() - 1)) {
      // Copy the current alignment into a std::vector<bool> for hashing.
      //
      // TODO: implement the std::vector<bool> hash function
      // for bm::bvector<> and use bm::copy_range?
      //
      std::vector<bool> current_ec(full_alignment->n_targets(), false);
      for (size_t j = 0; j < full_alignment->n_targets(); ++j) {
	current_ec[j] = ec_configs[i*full_alignment->n_targets() + j];
      }

      // Insert the current equivalence class to the hash map or
      // increment its observation count by 1 if it already exists.
      full_alignment->insert(current_ec, i, &ec_id, &ec_to_pos, &bv_it);
    }
  }
  bv_it.flush(); // Insert everything
  ec_configs.swap(compressed_ec_configs);
}

namespace read {
void Themisto(const Mode &mode, std::vector<std::istream*> &streams, CompressedAlignment *aln) {
  // Read in only the ec_configs
  ReadPairedAlignments(mode, streams, aln->get(), aln);
  CompressAlignment(*aln->get(), aln);

  aln->add_trailing_zeros(aln->compressed_size(), aln->n_targets());
  aln->make_read_only();
}

void ThemistoGrouped(const Mode &mode, std::vector<std::istream*> &streams, GroupedAlignment *aln) {
  // Read in group counts
  bm::bvector<> ec_configs;
  bm::sparse_vector<uint16_t, bm::bvector<>> sparse_counts;
  aln->sparse_group_counts = &sparse_counts;
  ReadPairedAlignments(mode, streams, &ec_configs, aln);
  CompressAlignment(ec_configs, aln);

  aln->build_group_counts();
  // aln->make_read_only();
}

void ThemistoAlignedReads(const Mode &mode, std::vector<std::istream*> &streams, ThemistoAlignment *taln) {
  // Read in the ec_configs and which reads are assigned to which equivalence classes
  ReadPairedAlignments(mode, streams, taln->get(), taln);
  taln->fill_read_ids();
  CompressAlignment(*taln->get(), taln);

  taln->add_trailing_zeros(taln->compressed_size(), taln->n_targets());
  taln->make_read_only();
}

void ThemistoToKallisto(const Mode &mode, std::vector<std::istream*> &streams, KallistoAlignment *aln) {
  // Read in the ec_configs and fill the ec_ids vector
  Themisto(mode, streams, aln);

  aln->ec_ids = std::vector<uint32_t>(aln->compressed_size(), 0);
  for (uint32_t i = 0; i < aln->compressed_size(); ++i) {
    aln->ec_ids[i] = i;
  }
}
}
}
