#ifndef TELESCOPE_STRUCTS_HPP
#define TELESCOPE_STRUCTS_HPP

struct ec_info {
  std::vector<bool> pseudoalignment;
  uint32_t count = 0;
  uint16_t last_val = 0;
};

#endif
