#ifndef TELESCOPE_READ_THEMISTO_ALIGNMENTS_HPP
#define TELESCOPE_READ_THEMISTO_ALIGNMENTS_HPP

#include <vector>
#include <fstream>

#include "common.hpp"

void ReadThemistoFiles(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, std::vector<uint32_t> *ec_counts, std::vector<std::vector<bool>> *compressed_ec_configs);

#endif
