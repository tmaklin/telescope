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
#include <memory>

#include "bm64.h"
#include "bmserial.h"
#include "bmsparsevec.h"
#include "unpack.hpp"

namespace telescope {
class Alignment {
protected:
  uint32_t n_processed;
  size_t n_refs;
  std::vector<uint32_t> ec_counts;
  bool parsing_from_buffered = false; // TODO use a function pointer to the correct parse function?

  std::vector<uint32_t> read_ids;
  std::vector<std::vector<uint32_t>> aligned_reads;

public:
  size_t compressed_size() const { return ec_counts.size(); }
  uint32_t size() const { return this->compressed_size(); } // Backwards compatibility.
  uint32_t n_targets() const { return this->n_refs; }
  size_t n_reads() const { return this->n_processed; }
  size_t reads_in_ec(const size_t &ec_id) const { return this->ec_counts[ec_id]; }
  const std::vector<uint32_t>& reads_assigned_to_ec(const size_t &ec_id) const { return this->aligned_reads[ec_id]; }

  virtual void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator *bv_it) =0;

  void add_counts(const size_t &count) { this->ec_counts.emplace_back(count); }

  void set_parse_from_buffered(const bool _parse_from_buffered) { this->parsing_from_buffered = _parse_from_buffered; }
  bool parse_from_buffered() const { return this->parsing_from_buffered; }

  void set_n_reads(const size_t _n_reads) { this->n_processed = _n_reads; }

  void add_read(const size_t _read_id) { ++this->n_processed; }

};

class ThemistoAlignment : public Alignment{
protected:
  std::unique_ptr<bm::bvector<>> ec_configs; // TODO implement copy constructor

public:
  ThemistoAlignment() {
    this->n_processed = 0;
    this->ec_configs->set_new_blocks_strat(bm::BM_GAP);
    this->ec_configs.reset(new bm::bvector<>());
  };
  ThemistoAlignment(const size_t &_n_refs, bm::bvector<> *ec_configs) {
    this->n_refs = _n_refs;
    this->n_processed = 0;
    this->ec_configs.reset(ec_configs);
  }
  ThemistoAlignment(const size_t &_n_refs, const size_t &n_to_process, bm::bvector<> *ec_configs) {
    // Constructor with known final size for ec_configs
    this->n_refs = _n_refs;
    this->n_processed = 0;
    this->ec_configs.reset(ec_configs);
  }

  bool operator()(const size_t row, const size_t col) const { return (*this->ec_configs)[row*this->n_refs + col]; }

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
    for (size_t j = 0; j < this->n_groups; ++j) {
      for (size_t i = 0; i < this->ec_ids.size(); ++i) {
	this->ec_group_counts[j*this->ec_ids.size() + i] = (*this->sparse_group_counts)[this->ec_ids[i]*this->n_groups + j];
      }
    }
  }

  void insert(const std::vector<bool> &current_ec, const size_t &i, size_t *ec_id, std::unordered_map<std::vector<bool>, uint32_t> *ec_to_pos, bm::bvector<>::bulk_insert_iterator*) override {

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
      this->aligned_reads.emplace_back(std::vector<uint32_t>());
      ++(*ec_id);
    }
    this->ec_counts[it->second] += 1;
    this->aligned_reads[it->second].emplace_back(i);
  }

  uint16_t get_group_count(const size_t row, const size_t col) {
    return this->ec_group_counts[row*this->ec_ids.size() + col];
  }
};
}

#endif
