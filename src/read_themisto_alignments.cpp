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

#include "bm64.h"
#include "unpack.hpp"

#include "telescope.hpp"

namespace telescope {
void ReadCompactAlignment(std::istream *stream, bm::bvector<> *ec_configs) {
  std::string line;
  while (std::getline(*stream, line)) {
    size_t next_buffer_size = std::stoul(line);
    alignment_writer::DeserializeBuffer(next_buffer_size, stream, ec_configs);
  }
}

void ReadPlaintextLine(const size_t n_targets, std::string &line, bm::bvector<>::bulk_insert_iterator &it) {
  std::string part;
  std::stringstream partition(line);

  // First column is read id
  std::getline(partition, part, ' ');
  size_t read_id = std::stoul(part);

  while (std::getline(partition, part, ' ')) {
    *it = read_id*n_targets + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
  }
}

size_t ReadPlaintextAlignment(const size_t n_targets, std::string &line, std::istream *stream, bm::bvector<> *ec_configs) {
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
  ReadPlaintextLine(n_targets, line, it);
  size_t n_reads = 1;

  while (std::getline(*stream, line)) {
    ReadPlaintextLine(n_targets, line, it);
    ++n_reads;
  }
  return n_reads;
}

size_t ReadAlignmentFile(const size_t n_targets, std::istream *stream, bm::bvector<> *ec_configs) {
  // Wrapper for determining which file format is used
  std::string line;
  std::getline(*stream, line);
  size_t n_reads;
  if (line.find(',') != std::string::npos) {
    // First line contains a ','; stream could be in the compact format.
    size_t n_refs;
    alignment_writer::ReadHeader(line, &n_reads, &n_refs);
    ec_configs->resize(n_reads*n_refs);
    ReadCompactAlignment(stream, ec_configs);
  } else {
    // Stream could be in the plaintext format.
    n_reads = ReadPlaintextAlignment(n_targets, line, stream, ec_configs);
  }
  return n_reads;
}

size_t ReadPairedAlignments(const Mode &mode, const size_t n_targets, std::vector<std::istream*> &streams, bm::bvector<> *ec_configs) {
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

  for (uint8_t i = 0; i < n_streams; ++i) {
    if (i == 0) {
      // Read the first alignments in-place to the output variable.
      n_reads = ReadAlignmentFile(n_targets, streams[i], ec_configs);
    } else {
      // Read subsequent alignments into a ThemistoAlignment object.
      bm::bvector<> new_configs(n_reads*n_targets, bm::BM_GAP);
      size_t n_processed;
      n_processed = ReadAlignmentFile(n_targets, streams[i], &new_configs);

      if (n_processed != n_reads) {
	// Themisto's output from paired-end reads should contain the same amount of reads.
	throw std::runtime_error("Pseudoalignment files have different numbers of pseudoalignments.");
      }
      if (mode == m_intersection) {
	// m_intersection: both reads in a pair should align to be considered a match.
	(*ec_configs) &= new_configs;
      } else {
	// m_union or m_unpaired: count alignments regardless of pair's status.
	(*ec_configs) |= new_configs;
      }
    }
  }
  return n_reads;
}

namespace read {
ThemistoAlignment Themisto(const Mode &mode, const size_t n_refs, std::vector<std::istream*> &streams) {
  // Read in only the ec_configs
  bm::bvector<> ec_configs(bm::BM_GAP);
  size_t n_reads = ReadPairedAlignments(mode, n_refs, streams, &ec_configs);
  ThemistoAlignment aln(n_refs, n_reads, ec_configs);
  aln.collapse();
  return aln;
}

ThemistoAlignment ThemistoPlain(const Mode &mode, const size_t n_refs, std::vector<std::istream*> &streams) {
  // Read in the plain alignment without compacting to equivalence classes
  bm::bvector<> ec_configs(bm::BM_GAP);
  size_t n_reads = ReadPairedAlignments(mode, n_refs, streams, &ec_configs);
  ThemistoAlignment aln(n_refs, n_reads, ec_configs);
  return aln;
}

GroupedAlignment ThemistoGrouped(const Mode &mode, const size_t n_refs, const size_t n_groups, const std::vector<uint32_t> &group_indicators, std::vector<std::istream*> &streams) {
  // Read in group counts
  bm::bvector<> ec_configs(bm::BM_GAP);
  size_t n_reads = ReadPairedAlignments(mode, n_refs, streams, &ec_configs);
  GroupedAlignment aln(n_refs, n_groups, n_reads, group_indicators);
  aln.collapse(ec_configs);

  return aln;
}

KallistoAlignment ThemistoToKallisto(const Mode &mode, const size_t n_refs, std::vector<std::istream*> &streams) {
  // Read in the ec_configs and fill the ec_ids vector
  // Read in only the ec_configs
  bm::bvector<> ec_configs(bm::BM_GAP);
  size_t n_reads = ReadPairedAlignments(mode, n_refs, streams, &ec_configs);
  KallistoAlignment aln(n_refs, n_reads, ec_configs);
  aln.collapse();

  aln.ec_ids = std::vector<uint32_t>(aln.n_ecs(), 0);
  for (uint32_t i = 0; i < aln.n_ecs(); ++i) {
    aln.ec_ids[i] = i;
  }
  return aln;
}
}
}
