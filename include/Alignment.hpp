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

#include "bm64.h"

namespace telescope {
struct Alignment {
  std::vector<uint32_t> ec_counts;

  uint32_t n_processed;
  size_t n_refs;

  uint32_t size() const { return ec_counts.size(); }
  uint32_t n_targets() const { return this->n_refs; }
};

class CompressedAlignment : public Alignment{
private:
  bm::bvector<> ec_configs;

public:
  bool operator()(const size_t row, const size_t col) const { return this->ec_configs[row*this->n_refs + col]; }
  bm::bvector<>* get() { return &this->ec_configs; }

  void add_trailing_zeros() {
    if (this->ec_configs.size() != this->ec_counts.size()*this->n_refs) {
      this->ec_configs.resize(this->ec_counts.size()*this->n_refs); // add trailing zeros
    }
  }

  void optimize_storage() {
    this->ec_configs.optimize();
    this->ec_configs.freeze();
  }

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
};

struct GroupedAlignment : public CompressedAlignment {
  uint16_t n_groups;
  std::vector<uint16_t> group_indicators;
  std::vector<std::vector<uint16_t>> ec_group_counts;

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
  std::vector<uint32_t> read_ids;
  std::vector<std::vector<uint32_t>> aligned_reads;

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

};
}

#endif
