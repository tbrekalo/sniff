#include <cstdint>
#include <cstdlib>
#include <iostream>

// 3rd party dependencies
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "version.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects = 0;

// sniff
#include "algo.h"
#include "io.h"

int main(int argc, char** argv) {
  try {
    auto options =
        cxxopts::Options("sniff", "pair up potential reverse complement reads");
    /* clang-format off */
    options.add_options("general")
      ("h,help", "print help")
      ("v,version", "print version")
      ("t,threads", "number of threads to use",
        cxxopts::value<std::uint32_t>()->default_value("1"));
    options.add_options("heuristic")
      ("p,percent",
       "maximum allowed difference in length as % of shorter read's length",
       cxxopts::value<double>()->default_value("0.05"))
      ("l,sample-length", 
       "maximum sample length from beginning/end of  sequence",
       cxxopts::value<std::uint32_t>()->default_value("5000"))
      ("e,max-edit-distance", "maximum allowed edit distance between samples",
       cxxopts::value<std::uint32_t>()->default_value("100"));
    options.add_options("input")
      ("input", "input fasta/fastq file", cxxopts::value<std::string>());
    /* clang-format on */

    options.parse_positional({"input"});
    auto result = options.parse(argc, argv);

    auto early_quit = false;
    if (result.count("version")) {
      /* clang-format off */
      std::cout << sniff_VERSION_MAJOR << '.'
                << sniff_VERSION_MINOR << '.'
                << sniff_VERSION_PATCH << std::endl;
      /* clang-format on */
      early_quit = true;
    }

    if (result.count("help")) {
      std::cout << options.help() << std::endl;
      early_quit = true;
    }

    if (early_quit) {
      return EXIT_SUCCESS;
    }

    auto const n_threads = result["threads"].as<std::uint32_t>();
    auto const reads_path =
        std::filesystem::path(result["input"].as<std::string>());

    auto timer = biosoup::Timer();
    timer.Start();

    auto pairs = sniff::FindReverseComplementPairs(
        sniff::Config{.p = result["percent"].as<double>(),
                      .length = result["sample_length"].as<std::uint32_t>(),
                      .max_edit_distance =
                          result["max_edit_distance"].as<std::uint32_t>()},
        sniff::LoadReads(reads_path));

    fmt::print(stderr, "[sniff::main]({:12.3f})\n", timer.Stop());

  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
  }

  return EXIT_SUCCESS;
}
