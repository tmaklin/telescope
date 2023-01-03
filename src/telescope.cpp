// telescope: convert between Themisto and kallisto pseudoalignments
// Copyright (C) 2019 Tommi Mäklin (tommi@maklin.fi)
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

#include "telescope.hpp"

#include <string>
#include <algorithm>
#include <vector>
#include <exception>
#include <cstddef>
#include <chrono>

#include "cxxargs.hpp"
#include "cxxio.hpp"

#include "version.h"
#include "log.hpp"

namespace telescope {
bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

void parse_args(int argc, char* argv[], cxxargs::Arguments &args, cxxio::Out &log) {
  args.add_short_argument<std::vector<std::string>>('r', "Themisto pseudoalignment(s)");
  args.add_short_argument<std::string>('o', "Output file directory.");
  args.add_long_argument<uint32_t>("n-refs", "Number of reference sequences in the pseudoalignment.");
  args.add_long_argument<std::string>("format", "Output format (kallisto or themisto, default: kallisto", "kallisto");
  args.add_long_argument<telescope::Mode>("mode", "How to merge paired-end alignments (one of unpaired, union, intersection; default: unpaired)", telescope::m_unpaired);
  args.add_long_argument<bool>("read-compact", "Read alignments that have been compressed with alignment-writer (default: false).", false);
  args.add_long_argument<bool>("write-compact", "Write themisto format alignments in alignment-writer compressed format (default: true).", true);
  args.add_long_argument<bool>("cin", "Read the last alignment file from cin (default: false).", false);
  args.add_long_argument<bool>("silent", "Suppress status messages (default: false)", false);
  args.add_long_argument<bool>("help", "Print the help message.", false);
  if (CmdOptionPresent(argv, argv+argc, "--help")) {
    log << "\n" + args.help() << '\n' << '\n';
    log.flush();
  }
  args.parse(argc, argv);
  if (args.value<std::string>("format") != "kallisto" || args.value<std::string>("format") != "themisto") {
    throw std::invalid_argument("--format must be one of kallisto or themisto.");
  }
}
}

int main(int argc, char* argv[]) {
  telescope::Log log(std::cerr, !telescope::CmdOptionPresent(argv, argv+argc, "--silent"));
  cxxargs::Arguments args("telescope-" + std::string(TELESCOPE_BUILD_VERSION), "Usage: telescope -r <strand_1>,<strand_2> -o <output prefix> --n-refs <number of pseudoalignment targets>");
  log << args.get_program_name() + '\n';
  try {
    log << "Parsing arguments\n";
    parse_args(argc, argv, args, log);

    // Check that the input directories  exist and are accessible
    cxxio::directory_exists(args.value<std::string>('o'));
  } catch (std::exception &e) {
    log.verbose = true;
    log << "Parsing arguments failed:\n"
	<< std::string("\t") + std::string(e.what()) + "\n"
	<< "\trun telescope with the --help option for usage instructions.\n";
    log.flush();
    return 1;
  }

  log << "Reading Themisto alignments\n";
  std::vector<cxxio::In> infiles(args.value<std::vector<std::string>>('r').size());
  std::vector<std::istream*> infile_ptrs(infiles.size());
  for (size_t i = 0; i < args.value<std::vector<std::string>>('r').size(); ++i) {
    infiles.at(i).open(args.value<std::vector<std::string>>('r').at(i));
    infile_ptrs.at(i) = &infiles.at(i).stream();
  }
  if (args.value<bool>("cin")) {
    infile_ptrs.push_back(&std::cin);
  }

  uint32_t n_refs = args.value<uint32_t>("n-refs");
  telescope::ThemistoAlignment alignments(n_refs);

  if (args.value<bool>("read-compact")) {
    alignments.set_parse_from_buffered();
  }
  telescope::read::ThemistoAlignedReads(args.value<telescope::Mode>("mode"), infile_ptrs, &alignments);

  if (args.value<std::string>("format") == "kallisto") {
    log << "Writing Kallisto format alignments\n";
    telescope::KallistoRunInfo run_info(alignments);
    run_info.call = "";
    run_info.start_time = std::chrono::system_clock::to_time_t(log.start_time);
    for (int i = 0; i < argc; ++i) {
      run_info.call += argv[i];
      run_info.call += (i == argc - 1 ? "" : " ");
    }

    log << "Writing converted alignment\n";
    cxxio::Out ec_file(args.value<std::string>('o') + "/pseudoalignments.ec");
    cxxio::Out tsv_file(args.value<std::string>('o') + "/pseudoalignments.tsv");
    telescope::write::ThemistoToKallisto(alignments, &ec_file.stream(), &tsv_file.stream());

    log << "Writing read assignments to equivalence classes\n";
    cxxio::Out read_to_ref_file(args.value<std::string>('o') + "/read-to-ref.txt");
    telescope::write::ThemistoReadAssignments(alignments, &read_to_ref_file.stream());

    cxxio::Out run_info_file(args.value<std::string>('o') + "/run_info.json");
    telescope::write::KallistoInfoFile(run_info, 4, &run_info_file.stream());
  } else {
    log << "Writing Themisto format alignment\n";
    cxxio::Out alignment_file(args.value<std::string>('o') + ".aln");
    if (args.value<bool>("write-compact")) {
      alignments.write_alignment(&alignment_file.stream());
    } else {
      throw std::runtime_error("Writing plaintext Themisto alignments is currently unsupported, use alignment-writer to decompress the files.");
    }
  }

  log << "Done\n";
  log.flush();

  return 0;
}
