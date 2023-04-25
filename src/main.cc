#include <cstdint>
#include <cstdlib>
#include <iostream>

#include "biosoup/nucleic_acid.hpp"
#include "cxxopts.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects = 0;

int main(int argc, char** argv) {
  auto options =
      cxxopts::Options("sniff", "pair up potential reverse complement reads");
  /* clang-format off */
  options.add_options()
    ("h,help", "display help")
    ("t,threads", "number of threads available",
      cxxopts::value<std::uint32_t>()->default_value("1"))
    ("input", "input fasta/fastq file", cxxopts::value<std::string>());
  /* clang-format on */

  options.parse_positional({"input-file"});
  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return EXIT_SUCCESS;
  }

  auto const n_threads = result["threads"].as<std::uint32_t>();

  return EXIT_SUCCESS;
}
