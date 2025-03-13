# JetReconstruction.jl C-bindings Test Applications

Benchmark application using the JetReconstruction.jl C-bindings.

## Compilation

Configure and compile using CMake in the standard way, e.g.,

```sh
cmake -S . -B build
# or if libjectreconstion is installed in non-standard location
# JetReconstruction_DIR=<JetReconstruction-installation-location> cmake -S . -B build
cmake --build build
```

> [!NOTE]
> Make sure to find the same version of Julia that was used to compile the JetReconstruction package.
> Custom search path can be set with `Julia_ROOT`, for instance:
>
> ```sh
> JULIA_ROOT=${HOME}/.julia/juliaup/julia-nightly/ cmake -S . -B build
> ```
