#ifndef TELESCOPE_READ_THEMISTO_ALIGNMENTS_HPP
#define TELESCOPE_READ_THEMISTO_ALIGNMENTS_HPP

#include <vector>
#include <fstream>

#include "common.hpp"

void ReadThemistoFiles(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, CompressedAlignment *aln);
void ThemistoToKallisto(const Mode &mode, const uint32_t n_refs, std::vector<std::istream*> &streams, KallistoAlignment *aln);

#endif
