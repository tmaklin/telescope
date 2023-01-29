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

#include "unpack.hpp"

#include "telescope.hpp"

namespace telescope {
void ReadCompactAlignment(std::istream *stream, bm::bvector<> *ec_configs) {
  // telescope::ReadCompactAlignment
  //
  // Reads an alignment file that has been compacted with
  // alignment-writer (https://github.com/tmaklin/alignment-writer)
  // into `*ec_configs`.
  //
  // Input:
  //   `stream`: pointer to an istream opened on the pseudoalignment file.
  //     NOTE:   Use alignment_writer::ReadHeader before calling this function!
  //   `ec_configs`: pointer to the output variable that will contain the alignment.
  //
  std::string line;
  while (std::getline(*stream, line)) {
    // Deserialize each chunk in the file by ORing into ec_configs
    size_t next_buffer_size = std::stoul(line);
    alignment_writer::DeserializeBuffer(next_buffer_size, stream, ec_configs);
  }
}

void ReadPlaintextLine(const size_t n_targets, std::string &line, bm::bvector<>::bulk_insert_iterator &it) {
  // telescope::ReadPlaintextLine
  //
  // Reads a line in a plaintext alignment file from Themisto
  // (https://github.com/algbio/themisto) into the bm::bvector<> `it`
  // inserts to.
  //
  // Input:
  //   `n_targets`: number of pseudoalignment targets (reference
  //                sequences). It's not possible to infer this from the Themisto
  //                file format so has to be provided separately.
  //   `line`: the line from the alignment file to read in.
  //   `it`: insert iterator to the bm::bvector<> variable for storing the alignment.
  //
  std::string part;
  std::stringstream partition(line);

  // First column is read id (0-based indexing).
  std::getline(partition, part, ' ');
  size_t read_id = std::stoul(part);

  // Next columns contain the target sequence id (0-based indexing).
  while (std::getline(partition, part, ' ')) {
    *it = read_id*n_targets + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
  }
}

size_t ReadPlaintextAlignment(const size_t n_targets, std::string &line, std::istream *stream, bm::bvector<> *ec_configs) {
  // telescope::ReadPlaintextAlignment
  //
  // Reads a plaintext alignment file from Themisto
  // (https://github.com/algbio/themisto) into `*ec_configs` and
  // return the number of reads in the file (both unaligned and
  // aligned).
  //
  // Input:
  //   `n_targets`: number of pseudoalignment targets (reference
  //                sequences). It's not possible to infer this from the Themisto
  //                file format so has to be provided separately.
  //   `line`: read each line into this variable.
  //     NOTE: the contents of the *first* line should already be stored in this variable.
  //   `stream`: pointer to an istream opened on the pseudoalignment file.
  //   `ec_configs`: pointer to the output variable that will contain the alignment.
  // Output:
  //   `n_reads`: total number of reads in the pseudoalignment (unaligned + aligned).
  //
  bm::bvector<>::bulk_insert_iterator it(*ec_configs); // Bulk insert iterator buffers the insertions

  // Contents of the first line is already stored in `line`
  ReadPlaintextLine(n_targets, line, it);
  size_t n_reads = 1;

  while (std::getline(*stream, line)) {
    // Insert each line into the alignment
    ReadPlaintextLine(n_targets, line, it);
    ++n_reads;
  }
  return n_reads;
}

size_t ReadAlignmentFile(const size_t n_targets, std::istream *stream, bm::bvector<> *ec_configs) {
  // telescope::ReadAlignmentFile
  //
  // Wrapper for determining which file format (alignment-writer or
  // plaintext) is used. Returns the number of aligned + unaligned reads in the
  // pseudoalignment.
  //
  // Input:
  //   `n_targets`: number of pseudoalignment targets (reference
  //                sequences). It's not possible to infer this from the plaintext Themisto
  //                file format so has to be provided separately. If the file is in the
  //                compact format will check that the numbers match.
  //   `stream`: pointer to an istream opened on the pseudoalignment file.
  //   `ec_configs`: pointer to the output variable that will contain the alignment.
  // Output:
  //   `n_reads`: total number of reads in the pseudoalignment (unaligned + aligned).
  //
  std::string line;
  std::getline(*stream, line); // Read the first line to check the format
  size_t n_reads;
  if (line.find(',') != std::string::npos) {
    // First line contains a ','; stream could be in the compact format.
    size_t n_refs;
    alignment_writer::ReadHeader(line, &n_reads, &n_refs);
    if (n_refs > n_targets) {
      throw std::runtime_error("Pseudoalignment file has more target sequences than expected.");
    } else if (n_targets < n_refs) {
      throw std::runtime_error("Pseudoalignment file has less target sequences than expected.");
    }
    // Size is given on the header line.
    ec_configs->resize(n_reads*n_refs);
    ReadCompactAlignment(stream, ec_configs);
  } else {
    // Stream could be in the plaintext format.
    // Size is unknown.
    n_reads = ReadPlaintextAlignment(n_targets, line, stream, ec_configs);
  }
  return n_reads;
}

size_t ReadPairedAlignments(const bm::set_operation &merge_op, const size_t n_targets, std::vector<std::istream*> &streams, bm::bvector<> *ec_configs) {
  // telescope::ReadPairedAlignments
  //
  // Reads one or more pseudoalignment files from Themisto for
  // paired reads into `ec_configs`. Can be in plaintext or alignment-writer
  // format. Returns the number of reads (unaligned + aligned) in the
  // alignment.
  //
  // Input:
  //   `merge_op`: bm::set_OR for union or bm::set_AND for intersection of multiple alignmnet files
  //   `n_targets`: number of pseudoalignment targets (reference
  //                sequences). It's not possible to infer this from the plaintext Themisto
  //                file format so has to be provided separately. If the file is in the
  //                compact format will check that the numbers match.
  //   `streams`: vector of pointers to the istreams opened on the pseudoalignment files.
  //   `ec_configs`: pointer to the output variable that will contain the alignment.
  // Output:
  //   `n_reads`: total number of reads in the pseudoalignment (unaligned + aligned).
  //
  uint8_t n_streams = streams.size(); // Typically 1 (unpaired reads) or 2 (paired reads).
  size_t n_reads;

  for (uint8_t i = 0; i < n_streams; ++i) {
    if (i == 0) {
      // Read the first alignments in-place to the output variable.
      n_reads = ReadAlignmentFile(n_targets, streams[i], ec_configs);
    } else {
      // Initialize a temporary object for storing the alignments.
      bm::bvector<> new_configs(n_reads*n_targets, bm::BM_GAP);
      size_t n_processed;
      n_processed = ReadAlignmentFile(n_targets, streams[i], &new_configs);

      // Themisto's output from paired-end reads should contain the same amount of reads.
      if (n_processed != n_reads) {
	throw std::runtime_error("Pseudoalignment files have different numbers of pseudoalignments.");
      }

      // Merge the temporary into the the output `ec_configs`.
      if (merge_op == bm::set_AND) {
	// m_intersection: both reads in a pair should align to be considered a match.
	(*ec_configs) &= new_configs;
      } else if (merge_op == bm::set_OR){
	// m_union: count alignments regardless of pair's status.
	(*ec_configs) |= new_configs;
      } else {
	throw std::runtime_error("Unknown paired alignment merge mode.");
      }
    }
  }
  return n_reads;
}

namespace read {
ThemistoAlignment Themisto(const bm::set_operation &merge_op, const size_t n_refs, std::vector<std::istream*> &streams) {
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
  //
  bm::bvector<> ec_configs(bm::BM_GAP);
  size_t n_reads = ReadPairedAlignments(merge_op, n_refs, streams, &ec_configs);
  ThemistoAlignment aln(n_refs, n_reads, ec_configs);
  aln.collapse();
  return aln;
}

ThemistoAlignment ThemistoPlain(const bm::set_operation &merge_op, const size_t n_refs, std::vector<std::istream*> &streams) {
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
  //
  bm::bvector<> ec_configs(bm::BM_GAP);
  size_t n_reads = ReadPairedAlignments(merge_op, n_refs, streams, &ec_configs);
  ThemistoAlignment aln(n_refs, n_reads, ec_configs);
  return aln;
}

KallistoAlignment ThemistoToKallisto(const bm::set_operation &merge_op, const size_t n_refs, std::vector<std::istream*> &streams) {
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
  //
  bm::bvector<> ec_configs(bm::BM_GAP);
  size_t n_reads = ReadPairedAlignments(merge_op, n_refs, streams, &ec_configs);
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
