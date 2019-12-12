#ifndef TELESCOPE_WRITE_ALIGNMENTS_HPP
#define TELESCOPE_WRITE_ALIGNMENTS_HPP

#include <unordered_map>
#include <vector>
#include <fstream>

#include "common.hpp"

void WriteAlignments(const std::unordered_map<std::vector<bool>, ec_info> &ecs, std::ostream* ec_file, std::ostream* tsv_file);
void WriteReadToRef(const std::unordered_map<uint32_t, std::vector<uint16_t>> &read_to_ref, std::ostream* out);
void WriteRunInfo(const KAlignment &alignment, std::ostream* out, uint8_t indent_len = 4);

#endif
