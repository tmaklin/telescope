#include "telescope.hpp"

#include <string>

#include "cxxargs/include/cxxargs.hpp"

#include "file.hpp"

bool CmdOptionPresent(char **begin, char **end, const std::string &option) {
  return (std::find(begin, end, option) != end);
}

void parse_args(int argc, char* argv[], cxxargs::Arguments &args, File::Out &log) {
  args.add_short_argument<std::string>('1', "Themisto pseudoalignments for strand 1.");
  args.add_short_argument<std::string>('2', "Themisto pseudoalignments for strand 2.");
  args.add_short_argument<std::string>('o', "Output files prefix.");
  args.add_long_argument<uint32_t>("n-refs", "Number of reference sequences in the pseudoalignment.");
  args.add_long_argument<Mode>("mode", "How to merge paired-end alignments (one of unpaired, union, intersection; default: unpaired)", m_unpaired);
  args.add_long_argument<bool>("help", "Print the help message.", false);
  if (CmdOptionPresent(argv, argv+argc, "--help")) {
    log << args.help() << '\n' << '\n';
    log.flush();
  }
  args.parse(argc, argv);
}

int main(int argc, char* argv[]) {
  File::Out log(std::cerr);
  cxxargs::Arguments args("telescope", "Usage: telescope -1 <strand_1> -2 <strand_2> -o <output prefix> --n-refs <number of references>");
  try {
    log << "Parsing arguments" << '\n';
    parse_args(argc, argv, args, log);
  } catch (std::exception &e) {
    log << "Parsing arguments failed:"
	<< ' ' << e.what()
	<< "\nexiting" << '\n';
    log.flush();
    return 1;
  }
  log << "Reading Themisto alignments" << '\n';
  File::In strand_1(args.value<std::string>('1')); 
  File::In strand_2(args.value<std::string>('2'));
  const KAlignment &alignments = ReadAlignments(args.value<Mode>("mode"), args.value<uint32_t>("n-refs"), &strand_1.stream(), &strand_2.stream());

  log << "Writing converted alignments" << '\n';
  File::Out ec_file(args.value<std::string>('o') + ".ec");
  File::Out tsv_file(args.value<std::string>('o') + ".tsv");
  WriteAlignments(alignments.ecs, &ec_file.stream(), &tsv_file.stream());

  log << "Writing read assignments to equivalence classes" << '\n';
  File::Out read_to_ref_file(args.value<std::string>('o') + "_read-to-ref.txt");
  WriteReadToRef(alignments.read_to_ref, &read_to_ref_file.stream());

  log << "Done" << '\n';
  log.flush();

  return 0;
}
