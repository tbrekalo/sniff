#include <sys/resource.h>

#include <cstdlib>

// 3rd party dependencies
#include "biosoup/nucleic_acid.hpp"
#include "biosoup/timer.hpp"
#include "cxxopts.hpp"
#include "fmt/core.h"
#include "tbb/task_arena.h"
#include "version.h"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects = 0;

// sniff
#include "sniff/algo.h"
#include "sniff/io.h"

static auto GetPeakMemoryUsageKB() -> std::uint32_t {
  struct rusage rusage_info;
  getrusage(RUSAGE_SELF, &rusage_info);

  return rusage_info.ru_maxrss;
}

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
      ("a,alpha",
       "shorter read length as percentage of longer read lenght in pair",
        cxxopts::value<double>()->default_value("0.10"))
      ("b,beta", "minimum required coverage on each read",
        cxxopts::value<double>()->default_value("0.90"));
    options.add_options("mapping")
      ("k,kmer-length", "kmer length used in mapping",
        cxxopts::value<std::uint32_t>()->default_value("15"))
      ("w,window-length", "window length used in mapping",
        cxxopts::value<std::uint32_t>()->default_value("5"))
      ("f,frequent", "filter f most frequent kmers",
        cxxopts::value<double>()->default_value("0.0002"));
    options.add_options("input")
      ("input", "input fasta/fastq file", cxxopts::value<std::string>());
    /* clang-format on */

    options.positional_help("<reads>");
    options.parse_positional({"input"});
    options.show_positional_help();
    auto result = options.parse(argc, argv);

    auto early_quit = false;
    if (result.count("version")) {
      fmt::print(stderr, "{}.{}.{}\n", sniff_VERSION_MAJOR, sniff_VERSION_MINOR,
                 sniff_VERSION_PATCH);
      early_quit = true;
    }

    if (result.count("help")) {
      fmt::print(stderr, "{}\n", options.help());
      early_quit = true;
    }

    if (early_quit) {
      return EXIT_SUCCESS;
    }

    auto const n_threads = result["threads"].as<std::uint32_t>();
    auto const reads_path =
        std::filesystem::path(result["input"].as<std::string>());

    auto task_arena = tbb::task_arena(n_threads);
    auto timer = biosoup::Timer();
    timer.Start();

    task_arena.execute([&] {
      auto const cfg = sniff::Config{
          .alpha_p = result["alpha"].as<double>(),
          .beta_p = result["beta"].as<double>(),
          .filter_freq = result["frequent"].as<double>(),
          .kmer_len = result["kmer-length"].as<std::uint32_t>(),
          .window_len = result["window-length"].as<std::uint32_t>()
          };

      /* clang-format off */
        fmt::print(stderr,
          "[sniff]\n"
          "\tthreads: {}\n"
          "\talpha: {:1.2f}; beta: {:1.2f}\n"
          "\tfilter-freq: {}; k: {}; w: {};\n",
          n_threads,
          cfg.alpha_p, cfg.beta_p,
          cfg.filter_freq, cfg.kmer_len, cfg.window_len);
      /* clang-format on */

      auto pairs =
          sniff::FindReverseComplementPairs(cfg, sniff::LoadReads(reads_path));
      for (auto const& [lhs, rhs] : pairs) {
        fmt::print("{},{}\n", lhs, rhs);
      }
    });

    fmt::print(stderr, "[sniff::main]({:12.3f}) peak rss {:0.3f} GB\n",
               timer.Stop(), static_cast<double>(GetPeakMemoryUsageKB()) / 1e6);
  } catch (std::exception const& e) {
    fmt::print(stderr, "{}\n", e.what());
  }

  return EXIT_SUCCESS;
}
