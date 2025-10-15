// fastjet-finder.cc
// MIT Licenced, Copyright (c) 2023-2024 CERN
//
// Code to run and time the jet finding of against various
// HepMC3 input files

// Original version of this code Philippe Gras, IRFU
// Modified by Graeme A Stewart, CERN

#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <string>
#include <vector>
#include <chrono>

#include <unistd.h>
#include <stdlib.h>

// Program Options Parser Library (https://github.com/badaix/popl)
#include "popl.hpp"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/ReaderAscii.h"

#include "fastjet-utils.hh"

using namespace std;
using namespace popl;
using Time = std::chrono::high_resolution_clock;
using us = std::chrono::microseconds;

fastjet::ClusterSequence run_fastjet_clustering(std::vector<fastjet::PseudoJet> input_particles,
  fastjet::Strategy strategy, fastjet::JetAlgorithm algorithm, fastjet::RecombinationScheme recombine_scheme, 
  double R, double p) {

  fastjet::JetDefinition jet_definition;
  if (algorithm == fastjet::genkt_algorithm || algorithm == fastjet::ee_genkt_algorithm) {
    jet_definition = fastjet::JetDefinition(algorithm, R, p, recombine_scheme, strategy);
  } else if (algorithm == fastjet::ee_kt_algorithm) {
    jet_definition = fastjet::JetDefinition(algorithm, recombine_scheme, strategy);
  } else {
    jet_definition = fastjet::JetDefinition(algorithm, R, recombine_scheme, strategy);
  }

  // run the jet clustering with the above jet definition
  fastjet::ClusterSequence clust_seq(input_particles, jet_definition);

  return clust_seq;
}

void dump_clusterseq(fastjet::ClusterSequence clust_seq) {
  // Print out the contents of the cluster sequence, for debug purposes
  // N.B. Indexes counted from 1 (to match Julia)
  // Jets
  auto jets = clust_seq.jets();
  auto ijets = 1;
  for (auto jet: jets) {
    std::cout << ijets << ": px=" << jet.px() << " py=" << jet.py() << " pz=" << jet.pz() << " E=" << jet.E() << std::endl;
    ijets++;
  } 
  // History
  auto history = clust_seq.history();
  auto ihistory = 1;
  for (auto he: history) {
    std::cout << ihistory << ": " <<
      he.parent1+1 << " " << he.parent2+1 << " " << he.child+1 << " " << 
      he.dij << " " << he.max_dij_so_far << std::endl;
    ihistory++;
  }
}

int main(int argc, char* argv[]) {
  // Default values
  int maxevents = -1;
  int skip_events = 0;
  int trials = 1;
  string mystrategy = "Best";
  double power = -1.0;
  string alg = "";
  string recombine = "";
  double R = 0.4;
  string dump_file = "";

  OptionParser opts("Allowed options");
  auto help_option = opts.add<Switch>("h", "help", "produce help message");
  auto max_events_option = opts.add<Value<int>>("n", "maxevents", "Maximum events in file to process (-1 = all events)", maxevents, &maxevents);
  auto skip_events_option = opts.add<Value<int>>("", "skipevents", "Number of events to skip over (0 = none)", skip_events, &skip_events);
  auto trials_option = opts.add<Value<int>>("m", "nsamples", "Number of repeated trials", trials, &trials);
  auto strategy_option = opts.add<Value<string>>("s", "strategy", "Valid values are 'Best' (default), 'N2Plain', 'N2Tiled'", mystrategy, &mystrategy);
  auto alg_option = opts.add<Value<string>>("A", "algorithm", "Algorithm: AntiKt CA Kt GenKt EEKt Durham (overrides power)", alg, &alg);
  auto power_option = opts.add<Value<double>>("p", "power", "Algorithm p value, only used/needed for GenKt and EEKt", power, &power);
  auto radius_option = opts.add<Value<double>>("R", "radius", "Algorithm R parameter", R, &R);
  auto recombine_option = opts.add<Value<string>>("", "recombine", "Recombination scheme for jet merging", recombine, &recombine);
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
  fastjet::Strategy strategy = fastjet::Best;
  if (mystrategy == string("N2Plain")) {
    strategy = fastjet::N2Plain;
  } else if (mystrategy == string("N2Tiled")) {
    strategy = fastjet::N2Tiled;
  }

  auto algorithm = fastjet::antikt_algorithm;
  if (alg != "") {
    if (alg == "AntiKt") {
      algorithm = fastjet::antikt_algorithm;
      power = -1.0;
    } else if (alg == "CA") {
      algorithm = fastjet::cambridge_aachen_algorithm;
      power = 0.0;
    } else if (alg == "Kt") {
      algorithm = fastjet::kt_algorithm;
      power = 1.0;
    } else if (alg == "GenKt") {
      algorithm = fastjet::genkt_algorithm;
    } else if (alg == "Durham") {
      algorithm = fastjet::ee_kt_algorithm;
      power = 1.0;
    } else if (alg == "EEKt") {
      algorithm = fastjet::ee_genkt_algorithm;
    } else {
      std::cout << "Unknown algorithm type: " << alg << std::endl;
      exit(1);
    }
  }

  auto recombine_scheme = fastjet::RecombinationScheme::E_scheme;
  std::cout << recombine << std::endl;
  if (recombine == "pt_scheme") {
    recombine_scheme = fastjet::RecombinationScheme::pt_scheme;
  } else if (recombine == "pt2_scheme") {
    recombine_scheme = fastjet::RecombinationScheme::pt2_scheme;
  }

  std::cout << "Strategy: " << mystrategy << "; Power: " << power << "; Algorithm " << algorithm << 
    "; Recombine " << recombine_scheme << std::endl;

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
      auto cluster_sequence = run_fastjet_clustering(events[ievt], strategy, algorithm, recombine_scheme, R, power);

      vector<fastjet::PseudoJet> final_jets;
      if (ptmin_option->is_set()) {
        final_jets = cluster_sequence.inclusive_jets(ptmin_option->value());
      } else if (dijmax_option->is_set()) {
        final_jets = cluster_sequence.exclusive_jets(dijmax_option->value());
      } else if (njets_option->is_set()) {
        final_jets = cluster_sequence.exclusive_jets(njets_option->value());
      }

      if (dump_option->is_set() && trial==0) {
        // sort by pt so files can be compared
        final_jets = fastjet::sorted_by_pt(final_jets);
        fprintf(dump_fh, "Jets in processed event %zu\n", ievt+1);
    
        // print out the details for each jet
        for (unsigned int i = 0; i < final_jets.size(); i++) {
          fprintf(dump_fh, "%5u %15.10f %15.10f %15.10f\n",
          i, final_jets[i].rap(), final_jets[i].phi(),
          final_jets[i].perp());
        }

        // Dump the cluster sequence history content as well?
        if (debug_clusterseq_option->is_set()) {
          dump_clusterseq(cluster_sequence);
        }
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

  return 0;
}
