// MIT Licenced, Copyright (c) 2023-2025 CERN
//
// Utilities common to all fastjet codes used for benchmarking and validation

// Original version of this code Philippe Gras, IRFU
// Modified by Graeme A Stewart, CERN
#include "cjetreconstruction-utils.hh"
#include "JetReconstruction.h"
#include <cassert>
#include <iostream> // needed for io
#include <string>
#include <vector>

#include <unistd.h>

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/ReaderAscii.h"


std::vector<std::vector<jetreconstruction_PseudoJet>> read_input_events(const char* fname, long maxevents) {
  // Read input events from a HepMC3 file, return the events in a vector
  // Each event is a vector of initial particles
  
  HepMC3::ReaderAscii input_file (fname);

  int events_parsed = 0;

  std::vector<std::vector<jetreconstruction_PseudoJet>> events;

  while(!input_file.failed()) {
    if (maxevents >= 0 && events_parsed >= maxevents) break;

    std::vector<jetreconstruction_PseudoJet> input_particles;

    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);

    // Read event from input file
    input_file.read_event(evt);

    // If reading failed - exit loop
    if (input_file.failed()) break;

    ++events_parsed;
    input_particles.clear();
    input_particles.reserve(evt.particles().size());
    for(auto p: evt.particles()){
      if(p->status() == 1){
        auto jet = jetreconstruction_PseudoJet{};
        auto retv = jetreconstruction_PseudoJet_init(
            &jet, p->momentum().px(), p->momentum().py(), p->momentum().pz(),
            p->momentum().e());
        assert(retv ==
               jetreconstruction_StatusCode::JETRECONSTRUCTION_STATUSCODE_OK);
        input_particles.push_back(std::move(jet));
      }
    }
    events.push_back(input_particles);
  }

  std::cout << "Read " << events_parsed << " events from " << fname << std::endl;
  return events;
}
