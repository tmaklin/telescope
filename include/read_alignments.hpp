#ifndef TELESCOPE_READ_ALIGNMENTS_HPP
#define TELESCOPE_READ_ALIGNMENTS_HPP

#include <unordered_map>
#include <vector>
#include <fstream>

#include "common.hpp"

KAlignment ReadAlignments(const Mode &mode, const uint32_t n_refs, std::istream* strand_1, std::istream* strand_2);

#endif
