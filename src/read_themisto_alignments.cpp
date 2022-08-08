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

#include "bm64.h"

namespace telescope {

size_t ReadAlignmentFile(const size_t &n_refs, std::istream *stream, bm::bvector<> *alignment) {
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
  bm::bvector<>::bulk_insert_iterator it(*alignment);
  std::string line;
  size_t n_reads = 0;
  while (std::getline(*stream, line)) {
    std::string part;
    std::stringstream partition(line);
    std::getline(partition, part, ' ');
    while (std::getline(partition, part, ' ')) {
      it = n_reads*n_refs + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
    }
    ++n_reads; // assumes --sort-output was used when running `themisto pseudoalign`
  }

  if (alignment->size() != n_reads*n_refs) {
    // Add trailing zeros if last alignment did not hit the last reference sequence.
    alignment->resize(n_reads*n_refs);
  }
  alignment->optimize(); // Conserve memory.
  return n_reads;
}

uint32_t ReadPairedAlignments(const Mode &mode, const uint32_t &n_refs, std::vector<std::istream*> &streams, bm::bvector<> *alignment) {
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
  size_t n_reads = 0;

  for (uint8_t i = 0; i < n_streams; ++i) {
    if (i == 0) {
      // Read the first alignments in-place to the output variable.
      n_reads = ReadAlignmentFile(n_refs, streams[i], alignment);
    } else {
      // Read subsequent alignments into a new bitvector.
      // Size is now known since the paired files should have the same numbers of reads.
      bm::bvector<> pair_alignment(n_reads*n_refs, bm::BM_GAP);
      size_t n_reads_new = ReadAlignmentFile(n_refs, streams[i], &pair_alignment);

      if (n_reads != n_reads_new) {
	// Themisto's output from paired-end reads should contain the same amount of reads.
	throw std::runtime_error("Pseudoalignment files have different numbers of pseudoalignments.");
      }
      pair_alignment.freeze(); // Make the alignment read-only.

      if (mode == m_intersection) {
	// m_intersection: both reads in a pair should align to be considered a match.
	*alignment &= pair_alignment;
      } else {
	// m_union or m_unpaired: count alignments regardless of pair's status.
	*alignment |= pair_alignment;
      }
    }
  }
  // Conserve memory and make alignments read-only.
  alignment->optimize();
  alignment->freeze();

  return n_reads;
}

uint32_t ReadGroupedAlignments(const Mode &mode, const std::vector<uint16_t> &group_indicators, const uint32_t n_refs, const uint16_t n_groups, std::vector<std::istream*> &streams, std::vector<std::vector<uint16_t>> *ec_group_counts, std::vector<uint32_t> *ec_counts) {
  // Returns the number of reads processed
  uint32_t n_reads = 0;
  uint8_t n_streams = streams.size();
  std::vector<std::string> lines(n_streams);

  std::unordered_map<bm::bvector<>, uint32_t> ec_to_pos;

  uint32_t ec_pos;
  while (std::getline(*streams[0], lines[0])) {
    for (uint8_t i = 1; i < n_streams; ++i) {
      std::getline(*streams[i], lines[i]);
    }

    bm::bvector<> current_ec(n_refs);
    if (mode == m_intersection)
      current_ec.set_range(0, n_refs, true);

    std::vector<bm::bvector<>> proposed_ecs(n_streams, bm::bvector<>(n_refs));
    for (uint8_t i = 0; i < n_streams; ++i) {
      std::string part;
      std::stringstream partition(lines[i]);
      getline(partition, part, ' ');
      while (getline(partition, part, ' ')) {
	proposed_ecs[i][std::stoul(part)] = true;
      }
    }
    for (uint8_t i = 0; i < n_streams; ++i) {
      if (mode == m_intersection)
	current_ec &= proposed_ecs[i];
      else
	current_ec |= proposed_ecs[i];
    }
    if (current_ec.any()) {
      if (ec_to_pos.find(current_ec) == ec_to_pos.end()) {
        ec_group_counts->emplace_back(std::vector<uint16_t>(n_groups, 0));
	ec_counts->push_back(0);
	ec_to_pos.insert(std::make_pair(current_ec, ec_pos));
	++ec_pos;
	std::vector<uint16_t> *current_read = &ec_group_counts->back();
	for (uint32_t j = 0; j  < n_refs; ++j) {
	  (*current_read)[group_indicators[j]] += current_ec[j];
	}
      }
      (*ec_counts)[ec_to_pos[current_ec]] += 1;
    }
    ++n_reads;
  }
  return n_reads;
}

uint32_t ReadAndAssign(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, std::vector<bm::bvector<>> *ec_configs, std::vector<uint32_t> *aligned_reads) {
  // Returns the number of reads processed
  uint32_t n_reads = 0;
  uint8_t n_streams = streams.size();
  std::vector<std::string> lines(n_streams);

  while (std::getline(*streams[0], lines[0])) {
    for (uint8_t i = 1; i < n_streams; ++i) {
      std::getline(*streams[i], lines[i]);
    }
    bm::bvector<> current_ec(n_refs);
    if (mode == m_intersection)
      current_ec.set_range(0, n_refs, true);

    uint32_t read_id;
    std::vector<bm::bvector<>> proposed_ecs(n_streams, bm::bvector<>(n_refs));
    for (uint8_t i = 0; i < n_streams; ++i) {
      std::string part;
      std::stringstream partition(lines[i]);
      getline(partition, part, ' ');
      read_id = std::stoul(part);
      while (getline(partition, part, ' ')) {
	proposed_ecs[i][std::stoul(part)] = true;
      }
    }
    for (uint8_t i = 0; i < n_streams; ++i) {
      if (mode == m_intersection)
	current_ec &= proposed_ecs[i];
      else
	current_ec |= proposed_ecs[i];
    }
    if (current_ec.any()) {
      ec_configs->emplace_back(current_ec);
      aligned_reads->emplace_back(read_id);
    }
    ++n_reads;
  }
  return n_reads;
}

bm::bvector<> CompressAlignment(const bm::bvector<> &full_alignment, const size_t &n_refs, const size_t &num_alns, std::vector<uint32_t> *ec_counts) {
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
  std::unordered_map<std::vector<bool>, uint32_t>::iterator it;

  size_t ec_id = 0;
  for (size_t i = 0; i < num_alns; ++i) {
    // Check if the current read aligned against any reference and
    // discard the read if it didn't.
    if (full_alignment.any_range(i*n_refs, i*n_refs + n_refs - 1)) {
      // Copy the current alignment into a std::vector<bool> for hashing.
      //
      // TODO: implement the std::vector<bool> hash function
      // for bm::bvector<> and use bm::copy_range?
      //
      std::vector<bool> current_ec(n_refs, false);
      for (size_t j = 0; j < n_refs; ++j) {
	current_ec[j] = full_alignment[i*n_refs + j];
      }

      // Check if the pattern has been observed
      it = ec_to_pos.find(current_ec);
      if (it == ec_to_pos.end()) {
	// Add new patterns to compressed_ec_configs.
	for (size_t j = 0; j < n_refs; ++j) {
	  if (current_ec[j]) {
	    bv_it = ec_id*n_refs + j;
	  }
	}
	// Add a new counter for the new pattern
	ec_counts->emplace_back(0);
	// Insert the new pattern into the hashmap
	it = ec_to_pos.insert(std::make_pair(current_ec, ec_id)).first; // return iterator to inserted element
	++ec_id;
      }
      (*ec_counts)[it->second] += 1; // Increment number of times the pattern was observed
    }
  }
  return compressed_ec_configs;
}

std::vector<bm::bvector<>> CompressAndAssign(const std::vector<uint32_t> &read_ids, const std::vector<bm::bvector<>> &ec_configs, std::vector<uint32_t> *ec_counts, std::vector<std::vector<uint32_t>> *aln_reads) {
  // Compress the alignment into equivalence classes and store the read assignments to equivalence classes
  std::vector<bm::bvector<>> compressed_ec_configs;
  uint32_t num_alns = ec_configs.size();
  std::unordered_map<bm::bvector<>, uint32_t> ec_to_pos;
  uint32_t ec_pos = 0;
  for (unsigned i = 0; i < num_alns; ++i) {
    if (ec_to_pos.find(ec_configs[i]) == ec_to_pos.end()) {
      compressed_ec_configs.push_back(ec_configs[i]);
      ec_counts->push_back(0);
      ec_to_pos.insert(std::make_pair(ec_configs[i], ec_pos));
      aln_reads->emplace_back(std::vector<uint32_t>());
      ++ec_pos;
    }
    uint32_t current_pos = ec_to_pos[ec_configs[i]];
    (*ec_counts)[current_pos] += 1;
    (*aln_reads)[current_pos].emplace_back(read_ids[i]);
  }
  return compressed_ec_configs;
}

namespace read {
void Themisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, CompressedAlignment *aln) {
  // Read in only the ec_configs
  bm::bvector<> alignment;
  alignment.set_new_blocks_strat(bm::BM_GAP);
  aln->n_processed = ReadPairedAlignments(mode, n_refs, streams, &alignment);
  *aln->get() = CompressAlignment(alignment, n_refs, aln->n_processed, &aln->ec_counts);

  if (aln->get()->size() != aln->ec_counts.size()*n_refs) {
    aln->get()->resize(aln->ec_counts.size()*n_refs); // add trailing zeros
  }
  aln->get()->optimize();
  aln->get()->freeze();

  aln->n_refs = n_refs;
}

void ThemistoGrouped(const Mode &mode, const std::vector<uint16_t> &group_indicators, const uint32_t n_refs, const uint16_t n_groups, std::vector<std::istream*> &streams, GroupedAlignment *aln) {
  // Read in group counts
  aln->n_processed = ReadGroupedAlignments(mode, group_indicators, n_refs, n_groups, streams, &aln->ec_group_counts, &aln->ec_counts);
}

void ThemistoAlignedReads(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, ThemistoAlignment *taln) {
  // Read in the ec_configs and assign reads to equivalence classes
  std::vector<uint32_t> aligned_reads_ids;
  taln->n_processed = ReadAndAssign(mode, n_refs, streams, &taln->ec_configs, &aligned_reads_ids);
  taln->ec_configs = CompressAndAssign(aligned_reads_ids, taln->ec_configs, &taln->ec_counts, &taln->aligned_reads);
}

void ThemistoToKallisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, KallistoAlignment *aln) {
  // Read in the ec_configs and fill the ec_ids vector
    //  aln->n_processed = ReadAlignments(mode, n_refs, streams, &aln->ec_configs);
    //aln->ec_configs = CompressAlignment(aln->ec_configs, &aln->ec_counts);
  for (uint32_t i = 0; i < aln->size(); ++i) {
    aln->ec_ids.emplace_back(i);
  }
}
}
}
