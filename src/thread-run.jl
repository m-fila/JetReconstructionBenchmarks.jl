#! /usr/bin/env julia
"""
Run Julia jet reconstruction over input events, running
ovar all threads used.

This script is not intended to be run directly, but 
controlled by the thread-scan.{jl,sh} script.
"""

using ArgParse
using Logging
using JSON
using CSV
using DataFrames
using EnumX
using CodecZlib

using LorentzVectorHEP
using JetReconstruction

# Backends for the jet reconstruction
@enumx T=Backend Backends Julia FastJet
const AllBackends = [String(Symbol(x)) for x in instances(Backends.Backend)]

# Parsing for EnumX types
function ArgParse.parse_item(opt::DataType, x::AbstractString)
    s = tryparse(opt, x)
    if s === nothing
        throw(ErrorException("Invalid value for enum: $(x)"))
    end
    s
end

function julia_jet_process_threads(events::Vector{Vector{T}};
                                    ptmin::Float64 = 5.0,
                                    distance::Float64 = 0.4,
                                    p::Union{Real, Nothing} = nothing,
                                    algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt,
                                    strategy::RecoStrategy.Strategy,
                                    nsamples::Integer = 1,
                                    repeats::Int = 1) where T <: JetReconstruction.FourMomentum
    @info "Will process $(size(events)[1]) events"

    # Set consistent algorithm and power
    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
                                                                       algorithm = algorithm)

    n_events = length(events)

    # Warmup code if we are doing a multi-sample timing run
    if nsamples > 1
        @info "Doing initial warm-up run"
        for event in events
            _ = inclusive_jets(jet_reconstruct(event, R = distance, p = p,
                                               strategy = strategy); ptmin = ptmin)
        end
    end

    # Threading
    nthreads = Threads.nthreads()
    if nthreads > 1
        @info "Will use $nthreads threads"
    end

    # Now setup timers and run the loop
    cummulative_time = 0.0
    cummulative_time2 = 0.0
    lowest_time = typemax(Float64)

    for irun in 1:nsamples
        t_start = time_ns()
        Threads.@threads for event_counter ∈ 1:n_events * repeats
            event_idx = mod1(event_counter, n_events)
            inclusive_jets(jet_reconstruct(events[event_idx], R = distance, p = p,
                                           strategy = strategy), ptmin = ptmin)
        end
        t_stop = time_ns()
        dt_μs = convert(Float64, t_stop - t_start) * 1.e-3
        if nsamples > 1
            @info "$(irun)/$(nsamples) $(dt_μs)"
        end
        cummulative_time += dt_μs
        cummulative_time2 += dt_μs^2
        lowest_time = dt_μs < lowest_time ? dt_μs : lowest_time
    end

    mean = cummulative_time / nsamples
    cummulative_time2 /= nsamples
    if nsamples > 1
        sigma = sqrt(nsamples / (nsamples - 1) * (cummulative_time2 - mean^2))
    else
        sigma = 0.0
    end
    # Event rate in Hz (lowest time is in μs)
    event_rate = (n_events * repeats / lowest_time) * 1_000_000

    # Average time per event
    mean /= n_events * repeats
    sigma /= n_events * repeats
    lowest_time /= n_events * repeats
    # Why also record the lowest time? 
    # 
    # The argument is that on a "busy" machine, the run time of an application is
    # always TrueRunTime+Overheads, where Overheads is a nuisance parameter that
    # adds jitter, depending on the other things the machine is doing. Therefore
    # the minimum value is (a) more stable and (b) reflects better the intrinsic
    # code performance.
    nthreads, lowest_time, event_rate
end

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
        "--ptmin"
        help = "Minimum p_t for final jets (GeV)"
        arg_type = Float64
        default = 5.0

        "--distance", "-R"
        help = "Distance parameter for jet merging"
        arg_type = Float64
        default = 0.4

        "--algorithm", "-A"
        help = """Algorithm to use for jet reconstruction: $(join(JetReconstruction.AllJetRecoAlgorithms, ", "))"""
        arg_type = JetAlgorithm.Algorithm
        default = JetAlgorithm.AntiKt

        "--power", "-p"
        help = "Power value for jet reconstruction"
        arg_type = Float64

        "--strategy", "-S"
        help = """Strategy for the algorithm, valid values: $(join(JetReconstruction.AllJetRecoStrategies, ", "))"""
        arg_type = RecoStrategy.Strategy
        default = RecoStrategy.Best

        "--nsamples", "-m"
        help = "Number of measurement points to acquire."
        arg_type = Int
        default = 16

        "--repeats"
        help = "Run over whole event sample this number of times"
        arg_type = Int
        default = 1

        "--backend"
        help = """Backend to use for the jet reconstruction: $(join(AllBackends, ", "))"""
        arg_type = Backends.Backend
        default = Backends.Julia

        "--info"
        help = "Print info level log messages"
        action = :store_true

        "--debug"
        help = "Print debug level log messages"
        action = :store_true

        "file"
        help = "HepMC3 event file in to process"
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)
    if args[:debug]
        logger = ConsoleLogger(stdout, Logging.Debug)
    elseif args[:info]
        logger = ConsoleLogger(stdout, Logging.Info)
    else
        logger = ConsoleLogger(stdout, Logging.Warn)
    end
    global_logger(logger)

    # Try to read events into the correct type!
    if JetReconstruction.is_ee(args[:algorithm])
        JetType = EEJet
    else
        JetType = PseudoJet
    end
    events::Vector{Vector{JetType}} = read_final_state_particles(args[:file], JetType)
    if isnothing(args[:algorithm]) && isnothing(args[:power])
        @warn "Neither algorithm nor power specified, defaulting to AntiKt"
        args[:algorithm] = JetAlgorithm.AntiKt
    end

    nthreads, event_rate, time_per_event = julia_jet_process_threads(events, ptmin = args[:ptmin],
                                                distance = args[:distance],
                                                algorithm = args[:algorithm],
                                                p = args[:power],
                                                strategy = args[:strategy],
                                                nsamples = args[:nsamples], repeats = args[:repeats])

    println("$(nthreads),$(event_rate),$(time_per_event)")
end

main()
