#ifndef TELESCOPE_READ_ALIGNMENTS_HPP
#define TELESCOPE_READ_ALIGNMENTS_HPP

#include <fstream>

#include "common.hpp"

KAlignment ReadAlignments(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*>* strands);

#endif
