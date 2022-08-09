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

#include "bm64.h"

namespace telescope {
class Alignment {
protected:
  uint32_t n_processed;
  size_t n_refs;
  std::vector<uint32_t> ec_counts;

public:
  size_t compressed_size() const { return ec_counts.size(); }
  uint32_t size() const { return this->compressed_size(); } // Backwards compatibility.
  uint32_t n_targets() const { return this->n_refs; }
  size_t n_reads() const { return this->n_processed; }
  size_t reads_in_ec(const size_t &ec_id) const { return this->ec_counts[ec_id]; }

  virtual void parse(const std::string &line, bm::bvector<>::bulk_insert_iterator *it) =0;

  void clear_counts() {
    this->ec_counts.clear();
    this->ec_counts.shrink_to_fit();
  }

  void add_counts(const size_t &count) { this->ec_counts.emplace_back(count); }

  std::vector<uint32_t>::iterator ec_counts_begin() { return this->ec_counts.begin(); }
  std::vector<uint32_t>::iterator ec_counts_end() { return this->ec_counts.end(); }
};

class CompressedAlignment : public Alignment{
private:
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

  void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it) {
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

struct GroupedAlignment : public CompressedAlignment {
private:
  uint16_t n_groups;
  std::vector<uint16_t> group_indicators;
  std::vector<std::vector<uint16_t>> ec_group_counts;

public:
  void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it) {

    // Check if the pattern has been observed
    std::unordered_map<std::vector<bool>, uint32_t>::iterator it = ec_to_pos->find(current_ec);
    if (it == ec_to_pos->end()) {
      this->ec_counts.emplace_back(0);
      it = ec_to_pos->insert(std::make_pair(current_ec, *ec_id)).first;
      ++(*ec_id);

      this->ec_group_counts.emplace_back(std::vector<uint16_t>(this->n_groups, 0));
      std::vector<uint16_t> *current_read = &this->ec_group_counts.back();
      for (uint32_t j = 0; j  < this->n_refs; ++j) {
	(*current_read)[this->group_indicators[j]] += current_ec[j];
      }
    }
    this->ec_counts[it->second] += 1;
  }

};

struct ThemistoAlignment : public CompressedAlignment{
private:
  std::vector<uint32_t> read_ids;
  std::vector<std::vector<uint32_t>> aligned_reads;

public:
  using CompressedAlignment::CompressedAlignment;

  const std::vector<uint32_t>& reads_assigned_to_ec(const size_t &ec_id) const { return this->aligned_reads[ec_id]; }
  void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it) {
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
    // Store the read id
    this->read_ids.emplace_back(std::stoul(part));
    while (std::getline(partition, part, ' ')) {
      *it = this->n_processed*this->n_refs + std::stoul(part); // set bit `n_reads*n_refs + std::stoul(part)` as true
    }
    ++this->n_processed; // assumes --sort-output was used when running `themisto pseudoalign`
  }

};
}

#endif
