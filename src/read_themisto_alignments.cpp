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
#include <sstream>
#include <unordered_map>
#include <functional>

#include "bm64.h"

namespace telescope {

void ParseLine(const std::string &line, const size_t &n_refs, CompressedAlignment *alignment, bm::bvector<>::bulk_insert_iterator *it) {
  // telescope::ParseLine
  //
  // Parses a line in the pseudoalignment file.
  //
  std::string part;
  std::stringstream partition(line);
  // Skip read id (first column)
  std::getline(partition, part, ' ');
  while (std::getline(partition, part, ' ')) {
    *it = alignment->n_processed*n_refs + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
  }
  ++alignment->n_processed; // assumes --sort-output was used when running `themisto pseudoalign`
}

void ParseLine2(const std::string &line, const size_t &n_refs, CompressedAlignment *alignment, bm::bvector<>::bulk_insert_iterator *it) {
  // telescope::ParseLine
  //
  // Parses a line in the pseudoalignment file and includes the read ids.
  //
  std::string part;
  std::stringstream partition(line);
  std::getline(partition, part, ' ');
  // Store the read id
  static_cast<ThemistoAlignment*>(alignment)->read_ids.emplace_back(std::stoul(part));
  while (std::getline(partition, part, ' ')) {
    *it = alignment->n_processed*n_refs + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
  }
  ++alignment->n_processed; // assumes --sort-output was used when running `themisto pseudoalign`
}

template<typename Function>
void ReadAlignmentFile(const size_t &n_refs, const Function &Parser, std::istream *stream, CompressedAlignment *_alignment) {
  // telescope::ReadAlignmentFile
  //
  // Reads in a istream pointing to a single themisto pseudoalignment
  // as a contiguous bit vector.
  //
  // Input:
  //   `n_refs`: number of reference sequences in the used themisto index.
  //   `stream`: pointer to an istream opened on the pseudoalignment file.
  //   `alignment`: pointer to a bitvector that will contain the results (can be uninitialized).
  //
  // Output:
  //   `read_id`: total number of reads (aligned and unaligned) in the pseudoalignment file.
  //
  bm::bvector<> *alignment = _alignment->get();
  bm::bvector<>::bulk_insert_iterator it(*alignment);
  std::string line;
  _alignment->n_processed = 0;
  while (std::getline(*stream, line)) {
    Parser(line, n_refs, _alignment, &it);
  }

  if (alignment->size() != _alignment->n_processed*n_refs) {
    // Add trailing zeros if last alignment did not hit the last reference sequence.
    alignment->resize(_alignment->n_processed*n_refs);
  }
  alignment->optimize(); // Conserve memory.
}

template<typename Function>
void ReadPairedAlignments(const Mode &mode, const uint32_t &n_refs, Function &Parser, std::vector<std::istream*> &streams, CompressedAlignment *alignment) {
  // telescope::ReadPairedAlignments
  //
  // Reads in one or more pseudoalignment files from themisto for paired reads.
  //
  // Input:
  //   `n_refs`: number of reference sequences in the used themisto index.
  //   `streams`: pointers to istreams opened on the pseudoalignment files.
  //   `alignment`: pointer to a bitvector that will contain the results (can be uninitialized).
  //
  // Output:
  //   `n_reads`: total number of paired reads (aligned and unaligned) in the pseudoalignment files.
  //
  uint8_t n_streams = streams.size(); // Typically 1 (unpaired reads) or 2 (paired reads).
  alignment->n_processed = 0;
  alignment->get()->set_new_blocks_strat(bm::BM_GAP);

  for (uint8_t i = 0; i < n_streams; ++i) {
    if (i == 0) {
      // Read the first alignments in-place to the output variable.
      ReadAlignmentFile(n_refs, Parser, streams[i], alignment);
    } else {
      // Read subsequent alignments into a new bitvector.
      // Size is now known since the paired files should have the same numbers of reads.
      CompressedAlignment pair_alignment;
      *pair_alignment.get() = bm::bvector<>(alignment->n_processed*n_refs, bm::BM_GAP);
      ReadAlignmentFile(n_refs, Parser, streams[i], &pair_alignment);

      if (alignment->n_processed != pair_alignment.n_processed) {
	// Themisto's output from paired-end reads should contain the same amount of reads.
	throw std::runtime_error("Pseudoalignment files have different numbers of pseudoalignments.");
      }
      pair_alignment.get()->freeze(); // Make the alignment read-only.

      if (mode == m_intersection) {
	// m_intersection: both reads in a pair should align to be considered a match.
	*alignment->get() &= *pair_alignment.get();
      } else {
	// m_union or m_unpaired: count alignments regardless of pair's status.
	*alignment->get() |= *pair_alignment.get();
      }
    }
  }
  // Conserve memory and make alignments read-only.
  alignment->get()->optimize();
  alignment->get()->freeze();
}

void InsertPseudoalignment(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it, CompressedAlignment *alignment) {
  // Check if the pattern has been observed
  std::unordered_map<std::vector<bool>, uint32_t>::iterator it = ec_to_pos->find(current_ec);
  if (it == ec_to_pos->end()) {
    // Add new patterns to compressed_ec_configs.
    for (size_t j = 0; j < alignment->n_refs; ++j) {
      if (current_ec[j]) {
	*bv_it = (*ec_id)*alignment->n_refs + j;
      }
    }
    // Add a new counter for the new pattern
    alignment->ec_counts.emplace_back(0);
    // Insert the new pattern into the hashmap
    it = ec_to_pos->insert(std::make_pair(current_ec, *ec_id)).first; // return iterator to inserted element
    ++(*ec_id);
  }
  alignment->ec_counts[it->second] += 1; // Increment number of times the pattern was observed
}

void InsertPseudoalignment2(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it, CompressedAlignment *alignment) {
  // Check if the pattern has been observed
  std::unordered_map<std::vector<bool>, uint32_t>::iterator it = ec_to_pos->find(current_ec);
  if (it == ec_to_pos->end()) {
    // Add new patterns to compressed_ec_configs.
    for (size_t j = 0; j < alignment->n_refs; ++j) {
      if (current_ec[j]) {
	*bv_it = (*ec_id)*alignment->n_refs + j;
      }
    }
    // Add a new counter for the new pattern
    alignment->ec_counts.emplace_back(0);
    // Insert the new pattern into the hashmap
    it = ec_to_pos->insert(std::make_pair(current_ec, *ec_id)).first; // return iterator to inserted element
    static_cast<ThemistoAlignment*>(alignment)->aligned_reads.emplace_back(std::vector<uint32_t>());
    ++(*ec_id);
  }
  alignment->ec_counts[it->second] += 1; // Increment number of times the pattern was observed
  static_cast<ThemistoAlignment*>(alignment)->aligned_reads[it->second].emplace_back(static_cast<ThemistoAlignment*>(alignment)->read_ids[i]);
}

void InsertPseudoalignment3(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it, CompressedAlignment *alignment) {

  // Check if the pattern has been observed
  std::unordered_map<std::vector<bool>, uint32_t>::iterator it = ec_to_pos->find(current_ec);
  if (it == ec_to_pos->end()) {
    alignment->ec_counts.emplace_back(0);
    it = ec_to_pos->insert(std::make_pair(current_ec, *ec_id)).first;
    ++(*ec_id);

    static_cast<GroupedAlignment*>(alignment)->ec_group_counts.emplace_back(std::vector<uint16_t>(static_cast<GroupedAlignment*>(alignment)->n_groups, 0));
    std::vector<uint16_t> *current_read = &static_cast<GroupedAlignment*>(alignment)->ec_group_counts.back();
    for (uint32_t j = 0; j  < alignment->n_refs; ++j) {
      (*current_read)[static_cast<GroupedAlignment*>(alignment)->group_indicators[j]] += current_ec[j];
    }
  }
  alignment->ec_counts[it->second] += 1;
}

template<typename Function>
void CompressAlignment(const Function &Inserter, CompressedAlignment *full_alignment) {
  // telescope::CompressAlignment
  //
  // Compresses the full pseudoalignment data into equivalence
  // classes, meaning unique pseudoalignment patterns and the numbers
  // of times they were observed.
  //
  // Input:
  //   `full_alignment`: the pseudoalignment as a bitvector (output from ReadPairedAlignment or ReadAlignmentFile).
  //   `n_refs`: number of reference sequences in the themisto index.
  //   `num_alns`: number of reads in the pseudoalignment file(s).
  //   `ec_counts`: pointer to a vector for storing the equivalence class observation counts.
  //
  // Output:
  //   `compressed_ec_configs`: unique pseudoalignment patterns that were observed
  //                            at least once. Guaranteed to be in the same order
  //                            as `ec_counts`.
  //
  bm::bvector<> compressed_ec_configs;
  compressed_ec_configs.set_new_blocks_strat(bm::BM_GAP); // Store data in compressed format.
  bm::bvector<>::bulk_insert_iterator bv_it(compressed_ec_configs);

  // Need to hash the alignment patterns to count the times they appear.
  std::unordered_map<std::vector<bool>, uint32_t> ec_to_pos;

  size_t ec_id = 0;
  for (size_t i = 0; i < full_alignment->n_processed; ++i) {
    // Check if the current read aligned against any reference and
    // discard the read if it didn't.
    if (full_alignment->get()->any_range(i*full_alignment->n_refs, i*full_alignment->n_refs + full_alignment->n_refs - 1)) {
      // Copy the current alignment into a std::vector<bool> for hashing.
      //
      // TODO: implement the std::vector<bool> hash function
      // for bm::bvector<> and use bm::copy_range?
      //
      std::vector<bool> current_ec(full_alignment->n_refs, false);
      for (size_t j = 0; j < full_alignment->n_refs; ++j) {
	current_ec[j] = (*full_alignment->get())[i*full_alignment->n_refs + j];
      }

      // Insert the current equivalence class to the hash map or
      // increment its observation count by 1 if it already exists.
      Inserter(current_ec, i, &ec_id, &ec_to_pos, &bv_it, full_alignment);
    }
  }
  full_alignment->get()->swap(compressed_ec_configs);
}

void OptimizeStorage(const size_t &n_ecs, const size_t &n_refs, bm::bvector<> *bvec) {
  if (bvec->size() != n_ecs*n_refs) {
    bvec->resize(n_ecs*n_refs); // add trailing zeros
  }
  bvec->optimize();
  bvec->freeze();
}

namespace read {
void Themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, CompressedAlignment *aln) {
  // Read in only the ec_configs
  aln->n_refs = n_refs;
  ReadPairedAlignments(mode, n_refs, ParseLine, streams, aln);
  CompressAlignment(InsertPseudoalignment, aln);

  OptimizeStorage(aln->ec_counts.size(), n_refs, aln->get());
}

void ThemistoGrouped(const Mode &mode, const std::vector<uint16_t> &group_indicators, const uint32_t n_refs, const uint16_t n_groups, std::vector<std::istream*> &streams, GroupedAlignment *aln) {
  // Read in group counts
  aln->n_refs = n_refs;
  aln->n_groups = n_groups;
  aln->group_indicators = group_indicators;
  ReadPairedAlignments(mode, n_refs, ParseLine, streams, aln);
  CompressAlignment(InsertPseudoalignment3, aln);
}

void ThemistoAlignedReads(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, ThemistoAlignment *taln) {
  // Read in the ec_configs and which reads are assigned to which equivalence classes
  taln->n_refs = n_refs;
  ReadPairedAlignments(mode, n_refs, ParseLine2, streams, taln);
  CompressAlignment(InsertPseudoalignment2, taln);

  OptimizeStorage(taln->ec_counts.size(), n_refs, taln->get());
}

void ThemistoToKallisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, KallistoAlignment *aln) {
  // Read in the ec_configs and fill the ec_ids vector
  Themisto(mode, n_refs, streams, aln);

  aln->ec_ids = std::vector<uint32_t>(aln->ec_counts.size(), 0);
  for (uint32_t i = 0; i < aln->ec_counts.size(); ++i) {
    aln->ec_ids[i] = i;
  }
}
}
}
