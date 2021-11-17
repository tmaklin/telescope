#ifndef TELESCOPE_READ_KALLISTO_ALIGNMENTS_HPP
#define TELESCOPE_READ_KALLISTO_ALIGNMENTS_HPP

#include <cstddef>
#include <fstream>

#include "common.hpp"

namespace telescope {
namespace read {
// Read pseudoalignments and counts
void Kallisto(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, CompressedAlignment *aln);
// Also read the equivalence class IDs
void KallistoEcIds(const uint32_t n_refs, std::istream &ec_file, std::istream &tsv_file, KallistoAlignment *aln);
}
}

#endif
