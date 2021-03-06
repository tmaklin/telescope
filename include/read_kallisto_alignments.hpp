#ifndef TELESCOPE_READ_KALLISTO_ALIGNMENTS_HPP
#define TELESCOPE_READ_KALLISTO_ALIGNMENTS_HPP

#include <fstream>
#include <vector>

#include "common.hpp"

// Read pseudoalignments and counts
void ReadKallistoFiles(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, CompressedAlignment *aln);
// Also read the equivalence class IDs
void ReadKallistoFiles(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, KallistoAlignment *aln);

#endif
