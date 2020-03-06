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

#ifndef TELESCOPE_COMMON_HPP
#define TELESCOPE_COMMON_HPP

#include <exception>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <algorithm>

enum Mode { m_unpaired, m_union, m_intersection };
inline Mode get_mode(const std::string &mode_str) {
    if (mode_str == "unpaired") return m_unpaired;
    if (mode_str == "union") return m_union;
    if (mode_str == "intersection") return m_intersection;
    throw std::runtime_error("Unrecognized paired-end mode.");
}

namespace cxxargs {
  inline std::istream& operator>> (std::istream &in, Mode &t) {
    std::string in_val;
    in >> in_val;
    t = get_mode(in_val);
    return in;
  }
}

struct CompressedAlignment {
  std::vector<std::vector<bool>> ec_configs;
  std::vector<uint32_t> ec_counts;

  uint32_t n_processed;

  const uint32_t size() const { return ec_configs.size(); }
  const uint32_t n_targets() const { return ec_configs.at(0).size(); }
};

struct KallistoRunInfo {
  KallistoRunInfo() = default;
  KallistoRunInfo(uint32_t n_targets, uint32_t n_processed, uint32_t n_pseudoaligned) :
    n_targets(n_targets), n_processed(n_processed), n_pseudoaligned(n_pseudoaligned), p_pseudoaligned(((double)n_pseudoaligned/n_processed)*100) {};
  KallistoRunInfo(const CompressedAlignment &aln) {
    n_targets = aln.n_targets();
    n_processed = aln.n_processed;
    n_pseudoaligned = 0;
    n_unique = 0;
#pragma omp parallel for schedule(static) reduction(+:n_pseudoaligned) reduction(+:n_unique)
    for (uint32_t i = 0; i < aln.size(); ++i) {
      n_unique += (aln.ec_counts[i] == 1);
      n_pseudoaligned += aln.ec_counts[i];
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


struct KallistoAlignment : public CompressedAlignment{
  std::vector<uint32_t> ec_ids;
  KallistoRunInfo run_info;

  void fill_info() {
    run_info = KallistoRunInfo(n_targets(), n_processed, size());
  }
};

struct ThemistoAlignment : public CompressedAlignment{
  std::vector<std::vector<uint32_t>> aligned_reads;
};

#endif
