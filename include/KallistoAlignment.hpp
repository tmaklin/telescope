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

#ifndef TELESCOPE_KALLISTO_RUN_INFO_HPP
#define TELESCOPE_KALLISTO_RUN_INFO_HPP

#include <cstddef>
#include <string>
#include <chrono>

#include "Alignment.hpp"

namespace telescope {
struct KallistoRunInfo {
  KallistoRunInfo() = default;
  KallistoRunInfo(uint32_t n_targets, uint32_t n_processed, uint32_t n_pseudoaligned) :
    n_targets(n_targets), n_processed(n_processed), n_pseudoaligned(n_pseudoaligned), p_pseudoaligned(((double)n_pseudoaligned/n_processed)*100) {};
  KallistoRunInfo(const ThemistoAlignment &aln) {
    n_targets = aln.n_targets();
    n_processed = aln.n_reads();
    n_pseudoaligned = 0;
    n_unique = 0;
    for (uint32_t i = 0; i < aln.n_ecs(); ++i) {
      n_unique += (aln.reads_in_ec(i) == 1);
      n_pseudoaligned += aln.reads_in_ec(i);
    }
    p_unique = (double)n_unique/n_processed;
    p_pseudoaligned = (double)n_pseudoaligned/n_processed;
  }

  uint32_t n_targets;
  uint32_t n_bootstraps = 0;
  uint32_t n_processed;
  uint32_t n_pseudoaligned;
  uint32_t n_unique;
  double p_pseudoaligned;
  double p_unique;
  std::string kallisto_version = "0.45.0";
  std::string index_version = "0";
  std::time_t start_time;
  std::string call;
};

struct KallistoAlignment : public ThemistoAlignment{
  using ThemistoAlignment::ThemistoAlignment;
  std::vector<uint32_t> ec_ids;
  KallistoRunInfo run_info;

  void fill_info() {
    run_info = KallistoRunInfo(this->n_targets(), this->n_reads(), this->n_ecs());
  }

  void insert(const size_t &ec_id, const size_t &ec_count) {
    if (ec_count > 0) {
      this->ec_ids.emplace_back(ec_id);
      this->ec_counts.emplace_back(ec_count);
    }
  }

  size_t get_ec_id(const size_t &ec_pos) const { return this->ec_ids[ec_pos]; }
  bm::bvector<>::bulk_insert_iterator get_iterator() { return bm::bvector<>::bulk_insert_iterator(this->ec_configs); }

};
}

#endif
