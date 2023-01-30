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

#include <cstddef>
#include <vector>
#include <unordered_map>

#include "bm64.h"
#include "bmsparsevec.h"

namespace telescope {
class Alignment {
private:
  // Insert a pseudoalignment into the equivalence class format (varies by alignment type, implement in children).
  // Used by the public collapse() method to create the equivalence classes.
  virtual void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it) =0;

protected:
  // Number of reads in the alignment
  uint32_t n_processed;

  // Number of alignment targets
  size_t n_refs;

  // Number of times an alignment corresponding to each equivalence class was observed
  std::vector<uint32_t> ec_counts;

  // IDs of reads that are assigned to each equivalence class
  std::vector<std::vector<uint32_t>> aligned_reads;

public:
  // Collapse the argument alignment into equivalence classes and their observation counts.
  // Assumes that the internal variables `n_refs` and `n_processed` are the same as in the argument.
  // The logic for creating the equivalence classes must be implemented in the insert() method in
  // each realization of the base class.
  void collapse(bm::bvector<> &ec_configs) {
    bm::bvector<> compressed_ec_configs;
    compressed_ec_configs.set_new_blocks_strat(bm::BM_GAP); // Store data in compressed format.
    bm::bvector<>::bulk_insert_iterator bv_it(compressed_ec_configs);

    // Need to hash the alignment patterns to count the times they appear.
    std::unordered_map<std::vector<bool>, uint32_t> ec_to_pos;

    size_t ec_id = 0;
    for (size_t i = 0; i < this->n_reads(); ++i) {
      // Check if the current read aligned against any reference and
      // discard the read if it didn't.
      if (ec_configs.any_range(i*this->n_refs, i*this->n_refs + this->n_refs - 1)) {
	// Copy the current alignment into a std::vector<bool> for hashing.
	//
	// TODO: implement the std::vector<bool> hash function
	// for bm::bvector<> and use bm::copy_range?
	//
	std::vector<bool> current_ec(this->n_refs, false);
	for (size_t j = 0; j < this->n_refs; ++j) {
	  current_ec[j] = ec_configs[i*this->n_refs + j];
	}

	// Insert the current equivalence class to the hash map or
	// increment its observation count by 1 if it already exists.
	this->insert(current_ec, i, &ec_id, &ec_to_pos, &bv_it);
      }
    }
    bv_it.flush(); // Insert everything

    ec_configs.swap(compressed_ec_configs);
    ec_configs.optimize();
    ec_configs.freeze();
  }

  // Check if `row` aligned against `col`.
  virtual size_t operator()(const size_t row, const size_t col) const =0;

  // Get the total number of equivalence classes in the alignment
  size_t n_ecs() const { return ec_counts.size(); }

  // Get the dimensions of the alignment
  uint32_t n_targets() const { return this->n_refs; }
  size_t n_reads() const { return this->n_processed; }

  // Get number times an equivalence class was observed
  size_t reads_in_ec(const size_t &ec_id) const { return this->ec_counts[ec_id]; }

  // Get the IDs of reads assigned to an equivalence class
  const std::vector<uint32_t>& reads_assigned_to_ec(const size_t &ec_id) const { return this->aligned_reads[ec_id]; }

  // Get all aligned reads
  const std::vector<std::vector<uint32_t>>& get_aligned_reads() const { return this->aligned_reads; }
};

class ThemistoAlignment : public Alignment{
private:
  // Store the pseudoalignment as a n_reads (rows) x n_refs (columns) matrix
  bm::bvector<> ec_configs;

  // Implement insert() from the base class
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
    this->aligned_reads[it->second].emplace_back(i);
  }

public:
  ThemistoAlignment() = default;

  ThemistoAlignment(const size_t &_n_refs, bm::bvector<> &ec_configs) {
    this->n_refs = _n_refs;
    this->n_processed = 0;
    this->ec_configs = std::move(ec_configs);
  }

  ThemistoAlignment(const size_t &_n_refs, const size_t &_n_reads, bm::bvector<> &ec_configs) {
    // Constructor with known final size for ec_configs
    this->n_refs = _n_refs;
    this->n_processed = _n_reads;
    this->ec_configs = std::move(ec_configs);
  }

  // Check if ec_id `row` aligned against group `col`.
  size_t operator()(const size_t row, const size_t col) const override { return this->ec_configs[row*this->n_refs + col]; }

  // Collapse the stored pseudoalignment into equivalence classes and their observation counts.
  void collapse() { Alignment::collapse(this->ec_configs); }

  // Get the ec_configs
  const bm::bvector<> &get_configs() const { return this->ec_configs; }
};

template <typename T>
struct GroupedAlignment : public Alignment {
private:
  // Total number of reference groups
  uint16_t n_groups;

  // Vector reference sequence at <position> to the group at <value>
  std::vector<uint32_t> group_indicators;

  // Number of sequences in each group that reads belonging to an
  // equivalence class aligned against.
  bm::sparse_vector<T, bm::bvector<>> sparse_group_counts;

  // Implement insert() from the base class
  void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator*) override {
    // Check if the pattern has been observed
    std::unordered_map<std::vector<bool>, uint32_t>::iterator it = ec_to_pos->find(current_ec);
    if (it == ec_to_pos->end()) {
      this->ec_counts.emplace_back(0);
      it = ec_to_pos->insert(std::make_pair(current_ec, *ec_id)).first;

      size_t read_start = (*ec_id)*this->n_groups;
      for (uint32_t j = 0; j  < this->n_refs; ++j) {
	if (current_ec[j]) {
	  this->sparse_group_counts.inc(read_start + this->group_indicators[j]);
	}
      }
      this->aligned_reads.emplace_back(std::vector<uint32_t>());
      ++(*ec_id);
    }
    this->ec_counts[it->second] += 1;
    this->aligned_reads[it->second].emplace_back(i);
  }

public:
  // Default constructor
  GroupedAlignment() {
    this->sparse_group_counts = bm::sparse_vector<T, bm::bvector<>>();
  }

  GroupedAlignment(const size_t _n_refs, const size_t _n_groups, const std::vector<uint32_t> _group_indicators) {
    this->n_refs = _n_refs;
    this->n_groups = _n_groups;
    this->group_indicators = _group_indicators;
    this->n_processed = 0;
    this->sparse_group_counts = bm::sparse_vector<T, bm::bvector<>>();
  }

  GroupedAlignment(const size_t _n_refs, const size_t _n_groups, const size_t _n_reads, const std::vector<uint32_t> _group_indicators) {
    this->n_refs = _n_refs;
    this->n_groups = _n_groups;
    this->group_indicators = _group_indicators;
    this->n_processed = _n_reads;
    this->sparse_group_counts = bm::sparse_vector<T, bm::bvector<>>();
  }

  // Get the number of sequences in group_id that the ec_id aligned against.
  T get_group_count(const size_t group_id, const size_t ec_id) const {
    size_t pos = ec_id*this->n_groups + group_id;
    return this->sparse_group_counts[pos];
  }

  size_t operator()(const size_t row, const size_t col) const override { return this->get_group_count(row, col); }

};
}

#endif
