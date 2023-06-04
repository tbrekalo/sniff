#include <cstdlib>

#include "cxxopts.hpp"
#include "edlib.h"
#include "fmt/core.h"
#include "sniff/io.h"

int main(int argc, char** argv) {
  try {
    auto options = cxxopts::Options(
        "sniff", "add edit distance to reverse complement pairs");
    /* clang-format off */
    options.add_options("input")
      ("reads", "input fasta/fastq reads", cxxopts::value<std::string>())
      ("pairs", "csv file containing reverse complements",
       cxxopts::value<std::string>());
    /* clang-format on */
  } catch (std::exception const& e) {
    fmt::print("{}\n", e.what());
  }

  return EXIT_SUCCESS;
}
