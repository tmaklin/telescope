#include "telescope.hpp"

#include <string>

#include "cxxargs/include/cxxargs.hpp"
#include "cxxio/file.hpp"
#include "cxxio/log.hpp"

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

void parse_args(int argc, char* argv[], cxxargs::Arguments &args, File::Out &log) {
  args.add_short_argument<std::vector<std::string>>('r', "Themisto pseudoalignment(s)");
  args.add_short_argument<std::string>('o', "Output file directory.");
  args.add_long_argument<uint32_t>("n-refs", "Number of reference sequences in the pseudoalignment.");
  args.add_long_argument<Mode>("mode", "How to merge paired-end alignments (one of unpaired, union, intersection; default: unpaired)", m_unpaired);
  args.add_long_argument<bool>("silent", "Suppress status messages (default: false)", false);
  args.add_long_argument<bool>("help", "Print the help message.", false);
  if (CmdOptionPresent(argv, argv+argc, "--help")) {
    log << "\n" + args.help() << '\n' << '\n';
    log.flush();
  }
  args.parse(argc, argv);
}

int main(int argc, char* argv[]) {
  Log log(std::cerr, !CmdOptionPresent(argv, argv+argc, "--silent"));
  cxxargs::Arguments args("telescope", "Usage: telescope -r <strand_1>,<strand_2> -o <output prefix> --mode <merge mode> --n-refs <number of references>");
  log << args.get_program_name() + '\n';
  try {
    log << "Parsing arguments\n";
    parse_args(argc, argv, args, log);
  } catch (std::exception &e) {
    log.verbose = true;
    log << "Parsing arguments failed:\n"
	<< std::string("\t") + std::string(e.what()) + "\n"
	<< "\trun telescope with the --help option for usage instructions.\n";
    log.flush();
    return 1;
  }
  log << "Reading Themisto alignments\n";
  std::vector<File::In> infiles(args.value<std::vector<std::string>>('r').size());
  std::vector<std::istream*> infile_ptrs(infiles.size());
  for (size_t i = 0; i < args.value<std::vector<std::string>>('r').size(); ++i) {
    infiles.at(i).open(args.value<std::vector<std::string>>('r').at(i));
    infile_ptrs.at(i) = &infiles.at(i).stream();
  }
  KAlignment alignments = ReadAlignments(args.value<Mode>("mode"), args.value<uint32_t>("n-refs"), &infile_ptrs);
  alignments.call = "";
  for (size_t i = 0; i < argc; ++i) {
    alignments.call += argv[i];
    alignments.call += (i == argc - 1 ? "" : " ");
  }

  log << "Writing converted alignments\n";
  File::Out ec_file(args.value<std::string>('o') + "/pseudoalignments.ec");
  File::Out tsv_file(args.value<std::string>('o') + "/pseudoalignments.tsv");
  WriteAlignments(alignments.ecs, &ec_file.stream(), &tsv_file.stream());

  log << "Writing read assignments to equivalence classes\n";
  File::Out read_to_ref_file(args.value<std::string>('o') + "/read-to-ref.txt");
  WriteReadToRef(alignments.read_to_ref, &read_to_ref_file.stream());

  File::Out run_info_file(args.value<std::string>('o') + "/run_info.json");
  WriteRunInfo(alignments, &run_info_file.stream());

  log << "Done\n";
  log.flush();

  return 0;
}
