// MIT Licenced, Copyright (c) 2023-2025 CERN
//
// Code to run and time the jet finding of against various
// HepMC3 input files

// Original version of this code Philippe Gras, IRFU
// Modified by Graeme A Stewart, CERN

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>   // needed for io
#include <iostream> // needed for io
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <vector>

// Program Options Parser Library (https://github.com/badaix/popl)
#include "popl.hpp"

#include "JetReconstruction.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/ReaderAscii.h"

#include "cjetreconstruction-utils.hh"

#ifdef JETRECONSTRUCTION_COMPILER_PACKAGECOMPILER
extern "C" {
#include "julia_init.h"
}
#endif

using namespace std;
using namespace popl;
using Time = std::chrono::high_resolution_clock;
using us = std::chrono::microseconds;

// sorted into decreasing transverse momentum as in fastjet::sorted_by_pt
void sorted_by_pt(jetreconstruction_JetsResult &jets) {
  std::sort(jets.data, jets.data + jets.length,
            [](const auto &left, const auto &right) {
              return left._pt2 > right._pt2;
            });
}

jetreconstruction_ClusterSequence
run_clustering(std::vector<jetreconstruction_PseudoJet> input_particles,
               jetreconstruction_RecoStrategy strategy,
               jetreconstruction_JetAlgorithm algorithm, double R, double p) {

  auto clust_seq = jetreconstruction_ClusterSequence{};
  auto retv = jetreconstruction_jet_reconstruct(
      input_particles.data(), input_particles.size(), algorithm, R, strategy,
      &clust_seq);
  assert(retv == jetreconstruction_StatusCode::JETRECONSTRUCTION_STATUSCODE_OK);
  return clust_seq;
}

int main(int argc, char *argv[]) {
#ifdef JETRECONSTRUCTION_COMPILER_PACKAGECOMPILER
  init_julia(0, nullptr);
#endif

  // Default values
  int maxevents = -1;
  int skip_events = 0;
  int trials = 1;
  string mystrategy = "Best";
  double power = -1.0;
  string alg = "";
  double R = 0.4;
  string dump_file = "";

  OptionParser opts("Allowed options");
  auto help_option = opts.add<Switch>("h", "help", "produce help message");
  auto max_events_option = opts.add<Value<int>>("m", "maxevents", "Maximum events in file to process (-1 = all events)", maxevents, &maxevents);
  auto skip_events_option = opts.add<Value<int>>("", "skipevents", "Number of events to skip over (0 = none)", skip_events, &skip_events);
  auto trials_option = opts.add<Value<int>>("n", "trials", "Number of repeated trials", trials, &trials);
  auto strategy_option = opts.add<Value<string>>("s", "strategy", "Valid values are 'Best' (default), 'N2Plain', 'N2Tiled'", mystrategy, &mystrategy);
  auto power_option = opts.add<Value<double>>("p", "power", "Algorithm p value: -1=antikt, 0=cambridge_aachen, 1=inclusive kt; otherwise generalised Kt", power, &power);
  auto alg_option = opts.add<Value<string>>("A", "algorithm", "Algorithm: AntiKt CA Kt GenKt EEKt Durham (overrides power)", alg, &alg);
  auto radius_option = opts.add<Value<double>>("R", "radius", "Algorithm R parameter", R, &R);
  auto ptmin_option = opts.add<Value<double>>("", "ptmin", "pt cut for inclusive jets");
  auto dijmax_option = opts.add<Value<double>>("", "dijmax", "dijmax value for exclusive jets");
  auto njets_option = opts.add<Value<int>>("", "njets", "njets value for exclusive jets");
  auto dump_option = opts.add<Value<string>>("d", "dump", "Filename to dump jets to");
  auto debug_clusterseq_option = opts.add<Switch>("c", "debug-clusterseq", "Dump cluster sequence jet and history content");


  opts.parse(argc, argv);

  if (help_option->count() == 1) {
    cout << argv[0] << " [options] HEPMC3_INPUT_FILE" << endl;
    cout << endl;
	  cout << opts << "\n";
    cout << "Note the only one of ptmin, dijmax or njets can be specified!\n" << endl;
    exit(EXIT_SUCCESS);
  }

  const auto extra_args = opts.non_option_args();
  std::string input_file{};
  if (extra_args.size() == 1) {
    input_file = extra_args[0];
  } else if (extra_args.size() == 0) {
    std::cerr << "No <HepMC3_input_file> argument after options" << std::endl;
  } else {
    std::cerr << "Only one <HepMC3_input_file> supported" << std::endl;
  }

  // Check we only have 1 option for final jet selection
  auto sum = int(njets_option->is_set()) + int(dijmax_option->is_set()) + int(ptmin_option->is_set());
  if (sum != 1) {
    cerr << "One, and only one, of ptmin, dijmax or njets needs to be specified (currently " <<
      sum << ")" << endl;
    exit(EXIT_FAILURE);
  }

  // read in input events
  //----------------------------------------------------------
  auto events = read_input_events(input_file.c_str(), maxevents);

  // Set strategy
  auto strategy =
      jetreconstruction_RecoStrategy::JETRECONSTRUCTION_RECOSTRATEGY_BEST;
  if (mystrategy == string("N2Plain")) {
    strategy =
        jetreconstruction_RecoStrategy::JETRECONSTRUCTION_RECOSTRATEGY_N2PLAIN;
  } else if (mystrategy == string("N2Tiled")) {
    strategy =
        jetreconstruction_RecoStrategy::JETRECONSTRUCTION_RECOSTRATEGY_N2TILTED;
  }

  auto algorithm =
      jetreconstruction_JetAlgorithm::JETRECONSTRUCTION_JETALGORITHM_ANTIKT;
  if (alg != "") {
    if (alg == "AntiKt") {
      algorithm =
          jetreconstruction_JetAlgorithm::JETRECONSTRUCTION_JETALGORITHM_ANTIKT;
      power = -1.0;
    } else if (alg == "CA") {
      algorithm =
          jetreconstruction_JetAlgorithm::JETRECONSTRUCTION_JETALGORITHM_CA;
      power = 0.0;
    } else if (alg == "Kt") {
      algorithm =
          jetreconstruction_JetAlgorithm::JETRECONSTRUCTION_JETALGORITHM_KT;
      power = 1.0;
      //} else if (alg == "GenKt") {
      //  algorithm = fastjet::genkt_algorithm;
      //} else if (alg == "Durham") {
      //  algorithm = fastjet::ee_kt_algorithm;
      //  power = 1.0;
      //} else if (alg == "EEKt") {
      //  algorithm = fastjet::ee_genkt_algorithm;
    } else {
      std::cout << "Unknown algorithm type: " << alg << std::endl;
      exit(1);
    }
  } else {
    if (power == -1.0) {
      algorithm =
          jetreconstruction_JetAlgorithm::JETRECONSTRUCTION_JETALGORITHM_ANTIKT;
    } else if (power == 0.0) {
      algorithm =
          jetreconstruction_JetAlgorithm::JETRECONSTRUCTION_JETALGORITHM_CA;
    } else if (power == 1.0) {
      algorithm =
          jetreconstruction_JetAlgorithm::JETRECONSTRUCTION_JETALGORITHM_KT;
      //} else {
      //  algorithm = fastjet::genkt_algorithm;
    }
  }

  std::cout << "Strategy: " << mystrategy << "; Power: " << power << "; Algorithm " << algorithm << std::endl;

  auto dump_fh = stdout;
  if (dump_option->is_set()) {
    if (dump_option->value() != "-") {
      dump_fh = fopen(dump_option->value().c_str(), "w");
    }
  }

  double time_total = 0.0;
  double time_total2 = 0.0;
  double sigma = 0.0;
  double time_lowest = 1.0e20;
  for (long trial = 0; trial < trials; ++trial) {
    std::cout << "Trial " << trial << " ";
    auto start_t = std::chrono::steady_clock::now();
    for (size_t ievt = skip_events_option->value(); ievt < events.size(); ++ievt) {
      auto cluster_sequence =
          run_clustering(events[ievt], strategy, algorithm, R, power);
      auto final_jets = jetreconstruction_JetsResult{nullptr,0};
      if (ptmin_option->is_set()) {
        auto retv = jetreconstruction_inclusive_jets(
            &cluster_sequence, ptmin_option->value(), &final_jets);
        assert(retv ==
               jetreconstruction_StatusCode::JETRECONSTRUCTION_STATUSCODE_OK);
      } else if (dijmax_option->is_set()) {
        auto retv = jetreconstruction_exclusive_jets_dcut(
            &cluster_sequence, dijmax_option->value(), &final_jets);
        assert(retv ==
               jetreconstruction_StatusCode::JETRECONSTRUCTION_STATUSCODE_OK);
      } else if (njets_option->is_set()) {
        auto retv = jetreconstruction_exclusive_jets_njets(
            &cluster_sequence, njets_option->value(), &final_jets);
        assert(retv ==
               jetreconstruction_StatusCode::JETRECONSTRUCTION_STATUSCODE_OK);
      }
      sorted_by_pt(final_jets);

      if (dump_option->is_set() && trial == 0) {
        fprintf(dump_fh, "Jets in processed event %zu\n", ievt + 1);

        // print out the details for each jet
        for (unsigned int i = 0; i < final_jets.length; i++) {
          const auto &jet = final_jets.data[i];
          auto perp = std::sqrt(jet.px * jet.px + jet.py * jet.py);
          fprintf(dump_fh, "%5u %15.10f %15.10f %15.10f\n", i, jet._rap,
                  jet._phi, perp);
        }

        // Dump the cluster sequence history content as well?
         if (debug_clusterseq_option->is_set()) {
        //  dump_clusterseq(cluster_sequence);
        }
        jetreconstruction_ClusterSequence_free_members(&cluster_sequence);
        jetreconstruction_JetsResult_free_members(&final_jets);
      }
    }
    auto stop_t = std::chrono::steady_clock::now();
    auto elapsed = stop_t - start_t;
    auto us_elapsed = double(chrono::duration_cast<chrono::microseconds>(elapsed).count());
    std::cout << us_elapsed << " us" << endl;
    time_total += us_elapsed;
    time_total2 += us_elapsed*us_elapsed;
    if (us_elapsed < time_lowest) time_lowest = us_elapsed;
  }
  time_total /= trials;
  time_total2 /= trials;
  if (trials > 1) {
    sigma = std::sqrt(double(trials)/(trials-1) * (time_total2 - time_total*time_total));
  } else {
    sigma = 0.0;
  }
  double mean_per_event = time_total / events.size();
  double sigma_per_event = sigma / events.size();
  time_lowest /= events.size();
  std::cout << "Processed " << events.size() << " events, " << trials << " times" << endl;
  std::cout << "Total time " << time_total << " us" << endl;
  std::cout << "Time per event " << mean_per_event << " +- " << sigma_per_event << " us" << endl;
  std::cout << "Lowest time per event " << time_lowest << " us" << endl;

#ifdef JETRECONSTRUCTION_COMPILER_PACKAGECOMPILER
  shutdown_julia(0);
#endif
  return 0;
}
