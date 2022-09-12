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

#ifndef TELESCOPE_ALIGNMENT_HPP
#define TELESCOPE_ALIGNMENT_HPP

#include <vector>
#include <cstddef>
#include <unordered_map>
#include <sstream>
#include <cmath>

#include "bm64.h"
#include "bmserial.h"
#include "bmsparsevec.h"

namespace telescope {
class Alignment {
protected:
  uint32_t n_processed;
  size_t n_refs;
  std::vector<uint32_t> ec_counts;
  bool parsing_from_buffered = false; // TODO use a function pointer to the correct parse function?

public:
  size_t compressed_size() const { return ec_counts.size(); }
  uint32_t size() const { return this->compressed_size(); } // Backwards compatibility.
  uint32_t n_targets() const { return this->n_refs; }
  size_t n_reads() const { return this->n_processed; }
  size_t reads_in_ec(const size_t &ec_id) const { return this->ec_counts[ec_id]; }

  // Parse plaintext alignments (themisto format)
  virtual void parse(const std::string &line, bm::bvector<>::bulk_insert_iterator *it) =0;
  // Parse alignments compressed with alignment_writer::BufferedPack.
  // See https://github.com/tmaklin/alignmen-writer for details.
  virtual void parse(const std::string &buffer_size_line, std::istream *in, bm::bvector<> *out) =0;
  virtual void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it) =0;

  void clear_counts() {
    this->ec_counts.clear();
    this->ec_counts.shrink_to_fit();
  }

  void add_counts(const size_t &count) { this->ec_counts.emplace_back(count); }

  void set_parse_from_buffered() { this->parsing_from_buffered = true; }
  bool parse_from_buffered() const { return this->parsing_from_buffered; }

  std::vector<uint32_t>::iterator ec_counts_begin() { return this->ec_counts.begin(); }
  std::vector<uint32_t>::iterator ec_counts_end() { return this->ec_counts.end(); }
};

class CompressedAlignment : public Alignment{
protected:
  bm::bvector<> ec_configs;

public:
  CompressedAlignment() {
    this->n_processed = 0;
    this->ec_configs.set_new_blocks_strat(bm::BM_GAP);
  };
  CompressedAlignment(const size_t &_n_refs) {
    this->n_refs = _n_refs;
    this->n_processed = 0;
    this->ec_configs.set_new_blocks_strat(bm::BM_GAP);
  }
  CompressedAlignment(const size_t &_n_refs, const size_t &n_to_process) {
    // Constructor with known final size for ec_configs
    this->n_refs = _n_refs;
    this->n_processed = 0;
    this->ec_configs = bm::bvector<>(_n_refs*n_to_process, bm::BM_GAP);
  }

  bool operator()(const size_t row, const size_t col) const { return this->ec_configs[row*this->n_refs + col]; }
  bm::bvector<>* get() { return &this->ec_configs; }
  const bm::bvector<>& get() const { return this->ec_configs; }

  void add_trailing_zeros(const size_t &rows, const size_t &cols) {
    if (this->ec_configs.size() != rows*cols) {
      this->ec_configs.resize(rows*cols); // add trailing zeros
    }
  }

  void make_read_only() {
    this->ec_configs.optimize();
    this->ec_configs.freeze();
  }

  void optimize() { this->ec_configs.optimize(); }

  void clear_configs() { this->ec_configs.clear(true); } // free memory

  void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it) override {
    // Check if the pattern has been observed
    std::unordered_map<std::vector<bool>, uint32_t>::iterator it = ec_to_pos->find(current_ec);
    if (it == ec_to_pos->end()) {
      // Add new patterns to compressed_ec_configs.
      for (size_t j = 0; j < this->n_refs; ++j) {
	if (current_ec[j]) {
	  *bv_it = (*ec_id)*this->n_refs + j;
	}
      }
      // Add a new counter for the new pattern
      this->ec_counts.emplace_back(0);
      // Insert the new pattern into the hashmap
      it = ec_to_pos->insert(std::make_pair(current_ec, *ec_id)).first; // return iterator to inserted element
      ++(*ec_id);
    }
    this->ec_counts[it->second] += 1; // Increment number of times the pattern was observed
  }

  void parse(const std::string &line, bm::bvector<>::bulk_insert_iterator *it) override {
    // telescope::ParseLine
    //
    // Parses a line in the pseudoalignment file.
    //
    std::string part;
    std::stringstream partition(line);
    // Skip read id (first column)
    std::getline(partition, part, ' ');
    while (std::getline(partition, part, ' ')) {
      *it = this->n_processed*this->n_refs + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
    }
    ++this->n_processed; // assumes --sort-output was used when running `themisto pseudoalign`
  }

  void parse(const std::string &buffer_size_line, std::istream *in, bm::bvector<> *out) override {
    // telescope::ParseLine
    //
    // Parses a line in the pseudoalignment file.
    //
    size_t next_buffer_size = std::stoul(buffer_size_line);

    // Allocate space for the block
    char* cbuf = new char[next_buffer_size + 1];

    // Read the next block into buf
    in->read(cbuf, next_buffer_size);
    unsigned char* buf = reinterpret_cast<unsigned char*>(const_cast<char*>(cbuf));

    // Deserialize block (OR with old data in bits)
    bm::deserialize(*out, buf);

    bm::bvector<>::size_type last;
    out->find_reverse(last);
    size_t last_in_batch = std::ceil(last/n_refs) + 1;
    this->n_processed = (this->n_processed > last_in_batch ? this->n_processed : last_in_batch);
    delete[] cbuf;
  }

  void merge_pair(const Mode &mode, const CompressedAlignment &pair) {
    if (mode == m_intersection) {
      // m_intersection: both reads in a pair should align to be considered a match.
      this->ec_configs &= pair.get();
    } else {
      // m_union or m_unpaired: count alignments regardless of pair's status.
      this->ec_configs |= pair.get();
    }
  }
};

struct GroupedAlignment : public Alignment {
private:
  uint16_t n_groups;
  std::vector<uint32_t> group_indicators;

public:
  bm::sparse_vector<uint16_t, bm::bvector<>> *sparse_group_counts;
  std::vector<uint16_t> ec_group_counts;

  GroupedAlignment() {
    this->n_processed = 0;
  }

  GroupedAlignment(const size_t _n_refs, const size_t _n_groups, const std::vector<uint32_t> _group_indicators) {
    this->n_refs = _n_refs;
    this->n_groups = _n_groups;
    this->group_indicators = _group_indicators;
    this->n_processed = 0;
  }

  std::vector<size_t> ec_ids;

  void build_group_counts() {
    this->ec_group_counts.resize(this->ec_ids.size()*this->n_groups);
#pragma omp parallel for schedule(static)
    for (size_t j = 0; j < this->n_groups; ++j) {
      for (size_t i = 0; i < this->ec_ids.size(); ++i) {
	this->ec_group_counts[j*this->ec_ids.size() + i] = (*this->sparse_group_counts)[this->ec_ids[i]*this->n_groups + j];
      }
    }
  }

  void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it) override {

    // Check if the pattern has been observed
    std::unordered_map<std::vector<bool>, uint32_t>::iterator it = ec_to_pos->find(current_ec);
    if (it == ec_to_pos->end()) {
      this->ec_counts.emplace_back(0);
      this->ec_ids.emplace_back(*ec_id);
      it = ec_to_pos->insert(std::make_pair(current_ec, *ec_id)).first;

      for (size_t i = 0; i < this->n_groups; ++i) {
	this->ec_group_counts.emplace_back(0);
      }

      size_t read_start = (*ec_id)*this->n_groups;
      for (uint32_t j = 0; j  < this->n_refs; ++j) {
	if (current_ec[j]) {
	  this->sparse_group_counts->inc(read_start + this->group_indicators[j]);
	}
      }
      ++(*ec_id);
    }
    this->ec_counts[it->second] += 1;
  }

  const std::vector<uint16_t>& get_group_counts() const { return this->ec_group_counts; }
  uint16_t get_group_count(const size_t row, const size_t col) {
    return this->ec_group_counts[row*this->ec_ids.size() + col];
  }

  void reset_group_indicators(const std::vector<uint32_t> &new_group_indicators) {
    this->group_indicators.assign(new_group_indicators.begin(), new_group_indicators.end());
  }

  void resize(const size_t new_n_refs, const size_t new_n_groups) {
    this->n_refs = new_n_refs;
    this->n_groups = new_n_groups;
  }

  const std::vector<uint16_t>& get_ec_group_counts() const {
    return this->ec_group_counts;
  }

  void free_counts() {
    this->ec_group_counts.clear();
    this->ec_group_counts.shrink_to_fit();
  }

  void parse(const std::string &line, bm::bvector<>::bulk_insert_iterator *it) override {
    // telescope::ParseLine
    //
    // Parses a line in the pseudoalignment file.
    //
    std::string part;
    std::stringstream partition(line);
    // Skip read id (first column)
    std::getline(partition, part, ' ');
    while (std::getline(partition, part, ' ')) {
      *it = this->n_processed*this->n_refs + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
    }
    ++this->n_processed; // assumes --sort-output was used when running `themisto pseudoalign`
  }

  void parse(const std::string &buffer_size_line, std::istream *in, bm::bvector<> *out) override {
    // telescope::ParseLine
    //
    // Parses a line in the pseudoalignment file.
    //
    size_t next_buffer_size = std::stoul(buffer_size_line);

    // Allocate space for the block
    char* cbuf = new char[next_buffer_size + 1];

    // Read the next block into buf
    in->read(cbuf, next_buffer_size);
    unsigned char* buf = reinterpret_cast<unsigned char*>(const_cast<char*>(cbuf));

    // Deserialize block (OR with old data in bits)
    bm::deserialize(*out, buf);

    bm::bvector<>::size_type last;
    out->find_reverse(last);
    size_t last_in_batch = std::ceil(last/n_refs) + 1;
    this->n_processed = (this->n_processed > last_in_batch ? this->n_processed : last_in_batch);
    delete[] cbuf;
  }
};

struct ThemistoAlignment : public CompressedAlignment{
private:
  std::vector<uint32_t> read_ids;
  std::vector<std::vector<uint32_t>> aligned_reads;

public:
  using CompressedAlignment::CompressedAlignment;

  const std::vector<uint32_t>& reads_assigned_to_ec(const size_t &ec_id) const { return this->aligned_reads[ec_id]; }
  void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it) override {
    // Check if the pattern has been observed
    std::unordered_map<std::vector<bool>, uint32_t>::iterator it = ec_to_pos->find(current_ec);
    if (it == ec_to_pos->end()) {
      // Add new patterns to compressed_ec_configs.
      for (size_t j = 0; j < this->n_refs; ++j) {
	if (current_ec[j]) {
	  *bv_it = (*ec_id)*this->n_refs + j;
	}
      }
      // Add a new counter for the new pattern
      this->ec_counts.emplace_back(0);
      // Insert the new pattern into the hashmap
      it = ec_to_pos->insert(std::make_pair(current_ec, *ec_id)).first; // return iterator to inserted element
      this->aligned_reads.emplace_back(std::vector<uint32_t>());
      ++(*ec_id);
    }
    this->ec_counts[it->second] += 1; // Increment number of times the pattern was observed
    this->aligned_reads[it->second].emplace_back(this->read_ids[i]);
  }

  void parse(const std::string &line, bm::bvector<>::bulk_insert_iterator *it) override {
    // telescope::ParseLine
    //
    // Parses a line in the pseudoalignment file and includes the read ids.
    //
    std::string part;
    std::stringstream partition(line);
    std::getline(partition, part, ' ');

    while (std::getline(partition, part, ' ')) {
      *it = this->n_processed*this->n_refs + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
    }
    ++this->n_processed; // assumes --sort-output was used when running `themisto pseudoalign`
  }

  void parse(const std::string &buffer_size_line, std::istream *in, bm::bvector<> *out) override {
    // telescope::ParseLine
    //
    // Parses a line in the pseudoalignment file.
    //
    size_t next_buffer_size = std::stoul(buffer_size_line);

    // Allocate space for the block
    char* cbuf = new char[next_buffer_size + 1];

    // Read the next block into buf
    in->read(cbuf, next_buffer_size);
    unsigned char* buf = reinterpret_cast<unsigned char*>(const_cast<char*>(cbuf));

    // Deserialize block (OR with old data in bits)
    bm::deserialize(*out, buf);

    bm::bvector<>::size_type first;
    bm::bvector<>::size_type last;
    out->find(first);
    out->find_reverse(last);

    size_t last_in_batch = std::ceil(last/n_refs) + 1;

    this->n_processed = (this->n_processed > last_in_batch ? this->n_processed : last_in_batch);
    delete[] cbuf;
  }
    void fill_read_ids() {
	this->read_ids = std::vector<uint32_t>(this->n_processed, 0);
	for (size_t i = 0; i < this->n_processed; ++i) {
	    this->read_ids[i] = i;
	}
    }
};
}

#endif
