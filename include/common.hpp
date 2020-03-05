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

struct ec_info {
  std::vector<bool> pseudoalignment;
  uint32_t count = 0;
  uint16_t last_val = 0;
};

struct CompressedAlignment {
  std::vector<std::vector<bool>> ec_configs;
  std::vector<uint32_t> ec_counts;
};

class KallistoAlignment {
private:
  CompressedAlignment aln;
  std::vector<uint32_t> ec_ids;
public:
  void clear_ids() { ec_ids.clear(); }
  inline uint32_t size() { return this->aln.ec_configs.size(); }
  inline uint32_t n_targets() { return this->aln.ec_configs.at(0).size(); }

  inline const std::vector<uint32_t>& get_ec_ids() const { return this->ec_ids; }
  inline const std::vector<uint32_t>& get_ec_counts() const { return this->aln.ec_counts; }
  inline const std::vector<std::vector<bool>>& get_ec_configs() const { return this->aln.ec_configs; }
  inline CompressedAlignment* access_aln() { return &this->aln; }
  inline std::vector<uint32_t>* access_ec_ids() { return &this->ec_ids; }
  inline std::vector<std::vector<bool>>* access_ec_configs() { return &this->aln.ec_configs; }
  inline std::vector<uint32_t>* access_ec_counts() { return &this->aln.ec_counts; }
};

struct KAlignment {
  // Kallisto-style alignments
  const uint32_t n_bootstraps = 0;
  const double p_pseudoaligned = 0.0;
  const double p_unique = 0.0;
  const std::string kallisto_version = "0.45.0";
  const std::string index_version = "0";
  const std::time_t start_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  uint32_t n_targets;
  uint32_t n_processed;
  uint32_t n_pseudoaligned;
  uint32_t n_unique;
  std::string call;

  std::unordered_map<std::vector<bool>, ec_info> ecs;
  std::unordered_map<uint32_t, std::vector<uint16_t>> read_to_ref;
};

#endif
