#include <cstdlib>
#include <fstream>

#include "cxxopts.hpp"
#include "edlib.h"
#include "fmt/core.h"
#include "sniff/io.h"

int main(int argc, char** argv) {
  try {
    auto options = cxxopts::Options(
        "sniff", "add edit distance to reverse complement pairs");
    /* clang-format off */
    options.add_options("general")
      ("h,help", "print help")
      ("t,threads", "number of threads to use",
        cxxopts::value<std::uint32_t>()->default_value("1"));
    options.add_options("input")
      ("reads", "input fasta/fastq reads", cxxopts::value<std::string>())
      ("pairs",
       "csv file containing reverse complements;"
       "each read should appear in no more than one pair",
        cxxopts::value<std::string>());
    /* clang-format on */

    options.positional_help("<reads> <pairs>");
    options.parse_positional({"reads", "pairs"});
    options.show_positional_help();

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      fmt::print(stderr, "{}\n", options.help());
      return EXIT_SUCCESS;
    }

  } catch (std::exception const& e) {
    fmt::print("{}\n", e.what());
  }

  return EXIT_SUCCESS;
}
