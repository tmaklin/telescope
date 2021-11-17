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
void ReadThemistoFiles(const Mode &mode, const std::vector<uint16_t> &group_indicators, const uint32_t n_refs, const uint16_t n_groups, std::vector<std::istream*> &streams, GroupedAlignment *aln);

struct VectorHasher {
  // See https://stackoverflow.com/questions/10405030/c-unordered-map-fail-when-used-with-a-vector-as-key
  template <typename T>
  int operator()(const std::vector<T> &V) const {
    int hash = V.size();
    for(auto &i : V) {
      hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};
}

#endif
