#ifndef C_JETRECONSTRUCTION_UTILS_H_
#define C_JETRECONSTRUCTION_UTILS_H_
#include "JetReconstruction.h"
#include <vector>

std::vector<std::vector<jetreconstruction_PseudoJet>> read_input_events(const char* fname, long maxevents = -1);

#endif // C_JETRECONSTRUCTION_UTILS_H_
