#ifndef TELESCOPE_COMMON_HPP
#define TELESCOPE_COMMON_HPP

enum Mode { m_unpaired, m_union, m_intersection };
inline Mode get_mode(const std::string &mode_str) {
    if (mode_str == "unpaired") return m_unpaired;
    if (mode_str == "union") return m_union;
    if (mode_str == "intersection") return m_intersection;
    throw std::runtime_error("Unrecognized paired-end mode.");
}

struct ec_info {
  std::vector<bool> pseudoalignment;
  uint32_t count = 0;
  uint16_t last_val = 0;
};

struct KAlignment {
  // Kallisto-style alignments
  std::unordered_map<std::vector<bool>, ec_info> ecs;
  std::unordered_map<uint32_t, std::vector<uint16_t>> read_to_ref;
};

#endif
