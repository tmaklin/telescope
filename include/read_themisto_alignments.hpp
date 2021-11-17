#ifndef TELESCOPE_READ_THEMISTO_ALIGNMENTS_HPP
#define TELESCOPE_READ_THEMISTO_ALIGNMENTS_HPP

#include <cstddef>
#include <vector>
#include <fstream>

#include "common.hpp"

namespace telescope {
void ReadThemistoFiles(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, CompressedAlignment *aln);
void ReadThemistoFiles(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, KallistoAlignment *aln);
void ReadThemistoFiles(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, ThemistoAlignment *aln);
}

#endif
