# JetReconstruction.jl C-bindings Test Applications

Benchmark application using the JetReconstruction.jl statically compiled to a shared library.

## Compilation

Configure and compile using CMake in the standard way, e.g.,

```sh
cmake -S . -B build
# or if libjectreconstruction is installed in non-standard location
# JetReconstruction_DIR=<JetReconstruction-installation-location> cmake -S . -B build
cmake --build build
```

> [!NOTE]
> Make sure to find the same version of Julia that was used to compile the JetReconstruction package.
> Custom search path can be set with `Julia_ROOT`, for instance:
>
> ```sh
> Julia_ROOT=${HOME}/.julia/juliaup/julia-1.12/ cmake -S . -B build
> ```
>
> Also make sure that for Julia 1.12 at least GLIBCXX_3.4.30 (gcc 12.1.0) is used.
