#include "write_alignments.hpp"

#include <exception>
#include <iomanip>
#include <chrono>

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
	*out << ec.second[i] << (i == n_reads - 1 ? '\n' : ' ');
      }
    }
    out->flush();
  } else {
    throw std::runtime_error("Output directory is not writable.");
  }
}

void WriteRunInfo(const KAlignment &alignment, std::ostream* out, uint8_t indent_len) {
  std::string indent;
  for (uint8_t i = 0; i < indent_len; ++i) {
    indent += " ";
  }
  if (!out->fail()) {
    *out << "{" << '\n';
    *out << indent << "\"n_targets\": " << alignment.n_targets << ',' << '\n';
    *out << indent << "\"n_bootstraps\": " << alignment.n_bootstraps << ',' << '\n';
    *out << indent << "\"n_processed\": " << alignment.n_processed << ',' << '\n';
    *out << indent << "\"n_pseudoaligned\": " << alignment.n_pseudoaligned << ',' << '\n';
    *out << indent << "\"n_unique\": " << alignment.n_unique << ',' << '\n';
    *out << indent << "\"p_pseudoaligned\": " << std::fixed << std::setprecision(1) << alignment.p_pseudoaligned << ',' << '\n';
    *out << indent << "\"p_unique\": " << alignment.p_unique << ',' << '\n';
    *out << indent << "\"kallisto_version\": \"" << alignment.kallisto_version << '\"' << ',' << '\n';
    *out << indent << "\"index_version\": " << alignment.index_version << ',' << '\n';
    std::string time = std::asctime(std::localtime(&alignment.start_time));
    time.pop_back();
    *out << indent << "\"start_time\": \"" << time << '\"' << ',' << '\n';
    *out << indent << "\"call\": \"" << alignment.call << '\"' << '\n';
    *out << "}" << '\n';
    out->flush();
  } else {
    throw std::runtime_error("Output directory is not writable.");
  }
}
