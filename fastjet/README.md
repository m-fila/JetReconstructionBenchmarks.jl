# FastJet Test Applications

This code compiles a few small FastJet applications that are used for
benchmarking and validation of alternative implementations of clustering
algorithms.

The code requires the fastjet libraries (<https://fastjet.fr/>) as well as those
from HepMC3 (<https://gitlab.cern.ch/hepmc/HepMC3>).

Depending on your system setup you make need to add the path to the HepMC3 and
Fastjet CMake setup files to `CMAKE_PREFIX_PATH`. It may also be necessary to
set `siscone_DIR`.

## Compilation

Configure and compile using CMake in the standard way, e.g.,

```sh
cmake -S . -B build
cmake --build build
```

## Applications

### `fastjet-finder`

`fastjet-finder` is the main application and will run fastjet with the standard
suite of pp algorithms, optionally outputting a list of inclusive or exclusive
jets. This is the standard executable used to benchmark fastjet.

```sh
./fastjet-finder [options] HEPMC3_INPUT_FILE

Allowed options:
  -h, --help                  produce help message
  -n, --maxevents arg (=-1)   Maximum events in file to process (-1 = all events)
  -m, --nsamples arg (=1)     Number of repeated trials
  -s, --strategy arg (=Best)  Valid values are 'Best' (default), 'N2Plain', 'N2Tiled'
  -A, --algorithm arg         Valid values are 'AntiKt' 'CA' 'Kt' 'GenKt' 'EEKt' 'Durham'
  -p, --power arg             Algorithm p value, only needed for 'GenKt' and 'EEKt'
  -R, --radius arg (=0.4)     Algorithm R parameter
  --ptmin arg                 pt cut for inclusive jets
  --dijmax arg                dijmax value for exclusive jets
  --njets arg                 njets value for exclusive jets
  -d, --dump arg              Filename to dump jets to
  -c, --debug-clusterseq      Dump cluster sequence history content

Note that only one of ptmin, dijmax or njets can be specified!
```

### `fastjet2json.jl`

`fastjet2json.jl` script converts the text output from the fastjet applications
into a JSON format, more suitable for integration tests.

```sh
julia --project src/fastjet2json.jl fastjet_input_file json_output_file
```
