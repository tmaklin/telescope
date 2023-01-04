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
#include <memory>

#include "bm64.h"
#include "unpack.hpp"

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
      size_t next_buffer_size = std::stoul(line);
      alignment_writer::DeserializeBuffer(next_buffer_size, stream, ec_configs);
    } else {
      std::string part;
      std::stringstream partition(line);

      std::getline(partition, part, ' ');
      size_t read_id = std::stoul(part);

      while (std::getline(partition, part, ' ')) {
	it = read_id*alignment->n_targets() + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
      }
      alignment->add_read(read_id);
    }
  }
  it.flush();
}

void Merge(const bm::set_operation &op, const bm::bvector<> &tmp, const size_t n_refs, const size_t read_id, bm::bvector<> *ec_configs) {
  for (size_t j = 0; j < n_refs; ++j) {
    size_t pos = read_id*n_refs + j;
    if (op == bm::set_AND) {
      // Intersection mode
      (*ec_configs)[pos] = (*ec_configs)[pos] && tmp[pos];
    } else {
      // Union or unpaired mode
      (*ec_configs)[pos] =(*ec_configs)[pos] || tmp[pos];
    }
  }
}

void MergePairedAlignment(const bm::set_operation &op, const size_t n_refs, const bool parse_from_buffered, std::istream *stream, bm::bvector<> *ec_configs) {
  std::string line;
  while (std::getline(*stream, line)) {
    if (parse_from_buffered) {
      // buffereiden väliin jääviä linjaamattomia ei nyt intersektoida
      size_t next_buffer_size = std::stoul(line);
      bm::bvector<> tmp;
      alignment_writer::DeserializeBuffer(next_buffer_size, stream, &tmp);

      bm::bvector<>::size_type first;
      bm::bvector<>::size_type last;
      tmp.find(first);
      tmp.find_reverse(last);

      size_t read_id = std::floor(first/n_refs);
      for (size_t i = first; i <= last; ++i) {
	size_t next_read_id = std::floor((i + 1)/n_refs);
	if (next_read_id != read_id) {
	  Merge(op, tmp, n_refs, read_id, ec_configs);
	}
	read_id = next_read_id;
      }
    } else {
      std::string part;
      std::stringstream partition(line);

      std::getline(partition, part, ' ');
      size_t read_id = std::stoul(part);

      bm::bvector<> tmp;
      bm::bvector<>::bulk_insert_iterator it(tmp);
      while (std::getline(partition, part, ' ')) {
	it = read_id*n_refs + std::stoul(part);
      }
      it.flush();
      Merge(op, tmp, n_refs, read_id, ec_configs);
    }
  }
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

  size_t n_reads;
  size_t n_refs;
  if (alignment->parse_from_buffered()) {
    alignment_writer::ReadHeader(streams[0], &n_reads, &n_refs);
    ec_configs->resize(n_refs*n_reads);
    alignment->set_n_reads(n_reads);
  }

  for (uint8_t i = 0; i < n_streams; ++i) {
    if (i == 0) {
      // Read the first alignments in-place to the output variable.
      ReadAlignmentFile(streams[i], ec_configs, alignment);
    } else {
      if (alignment->parse_from_buffered()) {
       	alignment_writer::ReadHeader(streams[i], &n_reads, &n_refs); // TODO check that the numbers are equal?
	if (alignment->n_reads() != n_reads || alignment->n_targets() != n_refs) {
	  // Themisto's output from paired-end reads should contain the same amount of reads.
	  throw std::runtime_error("Pseudoalignment files have different numbers of pseudoalignments.");
	}

      }
      MergePairedAlignment((mode == m_intersection ? bm::set_AND : bm::set_OR), alignment->n_targets(), alignment->parse_from_buffered(), streams[i], ec_configs);
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
void Themisto(const Mode &mode, std::vector<std::istream*> &streams, ThemistoAlignment *aln) {
  // Read in only the ec_configs
  ReadPairedAlignments(mode, streams, aln->get(), aln);
  aln->fill_read_ids();
  CompressAlignment(*aln->get(), aln);

  aln->add_trailing_zeros(aln->compressed_size(), aln->n_targets());
  aln->make_read_only();
}

void ThemistoPlain(const Mode &mode, std::vector<std::istream*> &streams, ThemistoAlignment *aln) {
  // Read in the plain alignment without compacting to equivalence classes
  ReadPairedAlignments(mode, streams, aln->get(), aln);
  aln->add_trailing_zeros(aln->n_reads(), aln->n_targets());
  aln->make_read_only();
}

void ThemistoGrouped(const Mode &mode, std::vector<std::istream*> &streams, GroupedAlignment *aln) {
  // Read in group counts
  bm::bvector<> ec_configs;
  bm::sparse_vector<uint16_t, bm::bvector<>> sparse_counts;
  aln->sparse_group_counts = &sparse_counts;
  ReadPairedAlignments(mode, streams, &ec_configs, aln);
  aln->fill_read_ids();
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
