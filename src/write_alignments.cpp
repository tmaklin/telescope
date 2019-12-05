#include "write_alignments.hpp"

#include <exception>

void WriteAlignments(const std::unordered_map<std::vector<bool>, ec_info> &ecs, std::ostream* ec_file, std::ostream* tsv_file) {
  if (!ec_file->fail() && !tsv_file->fail()) {
    uint32_t ec_id = 0;
    for (auto ec : ecs) {
      *ec_file << ec_id << '\t';
      *tsv_file << ec_id << '\t';
      *tsv_file << ec.second.count << '\n';
      for (size_t i = 0; i < ec.first.size(); ++i) {
	if (ec.first.at(i)) {
	  *ec_file << i << (i == ec.second.last_val ? '\n' : ',');
	}
      }
      ++ec_id;
    }
    ec_file->flush();
    tsv_file->flush();
  } else {
    throw std::runtime_error("Output files are not writable.");
  }
}

void WriteReadToRef(const std::unordered_map<uint32_t, std::vector<uint16_t>> &read_to_ref, std::ostream* out) {
  if (!out->fail()) {
    for (auto ec : read_to_ref) {
      *out << ec.first << ' ';
      size_t n_reads = ec.second.size();
      for (size_t i = 0; i < n_reads; ++i) {
	*out << i << (i == n_reads ? '\n' : ' ');
      }
    }
    out->flush();
  } else {
    throw std::runtime_error("Output directory is not writable.");
  }
}
