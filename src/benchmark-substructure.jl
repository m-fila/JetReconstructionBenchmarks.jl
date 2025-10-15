using ArgParse
using JetReconstruction
using CSV
using DataFrames
# using Infiltrator

function makeJets(data::Vector{PseudoJet})
    finaljets = jet_reconstruct(data; p=0, R = 0.4)
    
    finaljets
end

"""
Read a file, reconstruct and select exclusive jets
"""
function get_exclusive_jet_selection(filename, ptmin, R)
    exclusive_jets = Vector{PseudoJet}[]
    cluster_seqs = ClusterSequence[]
    
    events = read_final_state_particles(filename)
    for event in events
        cs = jet_reconstruct(event; p=0, R = R)
        push!(exclusive_jets, inclusive_jets(cs, PseudoJet; ptmin=ptmin))
        push!(cluster_seqs, cs)
    end
    (cluster_seqs, exclusive_jets)
end

"""
Read all files in turn, constructing a dictionary of inclusive jets
"""
function process_input_files(filenames, ptmin, R)
    final_dict = Dict{String, Tuple{Vector{ClusterSequence}, Vector{Vector{PseudoJet}}}}()
    for file in filenames
        final_dict[file] = get_exclusive_jet_selection(file, ptmin, R)
    end
    final_dict
end

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
    "--filter"
    help = "Set filtering parameters"
    arg_type = JetFilter
    default = JetFilter(0.3, 3)
    
    "--trim"
    help = "Set trimming parameters"
    arg_type = JetTrim
    default = JetTrim(0.3, 0.3, JetAlgorithm.CA)
    
    "--massdrop"
    help = "Set mass drop parameters"
    arg_type = MassDropTagger
    default = MassDropTagger(0.67, 0.09)
    
    "--softdrop"
    help = "Set soft drop parameters"
    arg_type = SoftDropTagger
    default = SoftDropTagger(0.1, 2.0)
    
    "--ptmin"
    help = "pT cut for exclusive jet selection"
    arg_type = Float64
    default = 2.0
    
    "--distance", "-R"
    help = "Distance parameter for jet merging"
    arg_type = Float64
    default = 0.4
    
    "--nsamples", "-n"
    help = "Number of measurement points to acquire"
    arg_type = Int
    default = 16
    
    "--results"
    help = "Write results in CSV format to this file"
    
    "files"
    help = "Event files to process"
    required = true
    nargs = '+'
    
end
return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)
    
    if endswith(args[:files][1], ".csv")
        if length(args[:files]) > 1
            println("When CSV input is given, only one file can be used")
            exit(1)
        end
        hepmc3_files_df = CSV.read(args[:files][1], DataFrame)
        input_file_path = dirname(args[:files][1])
        input_files = hepmc3_files_df[:, :File]
        input_files_full = map(fname -> joinpath(input_file_path, fname), input_files)
        hepmc3_files_df[:, :File_path] = input_files_full
    else
        input_files_full = args[:files]
        hepmc3_files_df = DataFrame("File_path" => input_files_full)
        hepmc3_files_df[:, :File] = map(basename, hepmc3_files_df[:, :File_path])
    end
    
    excl_jets_dict = process_input_files(input_files_full, args[:ptmin], args[:distance])
    
    total_ejets_per_sample = Int[]
    for k in input_files_full
        jet_count = 0
        for event_jets in excl_jets_dict[k][2]
            # println(event_jets)
            jet_count += length(event_jets)
        end
        push!(total_ejets_per_sample, jet_count)
        @info "File $k, got $jet_count exclusive jets"
    end
    hepmc3_files_df[:, :ExJets] = total_ejets_per_sample
    
    mass_drop_sample = Float64[]
    soft_drop_sample = Float64[]
    filter_sample = Float64[]
    trim_sample = Float64[]

    for (f, (event_cseqs, event_jets)) in excl_jets_dict
        mass_drop_times = Float64[]
        soft_drop_times = Float64[]
        filtering = Float64[]
        trimming = Float64[]
        # @infiltrate
        n_ex_jets_v = hepmc3_files_df[hepmc3_files_df[!, :File_path] .== f, :ExJets]
        @assert length(n_ex_jets_v) == 1
        n_ex_jets = n_ex_jets_v[1]

        for _ in 1:args[:nsamples]
            md = sd = fi = tr = 0.0
            for (cs, ejets) in zip(event_cseqs, event_jets)
                fi += @elapsed [jet_filtering(jet, cs, args[:filter]) for jet in ejets]
                tr += @elapsed [jet_trimming(jet, cs, args[:trim]) for jet in ejets]
                md += @elapsed [mass_drop(jet, cs, args[:massdrop]) for jet in ejets]
                sd += @elapsed [soft_drop(jet, cs, args[:softdrop]) for jet in ejets]
            end
            push!(mass_drop_times, md * 10^6 / n_ex_jets)
            push!(soft_drop_times, sd * 10^6 / n_ex_jets)
            push!(filtering, fi * 10^6 / n_ex_jets)
            push!(trimming, tr * 10^6 / n_ex_jets)
        end
        push!(mass_drop_sample, minimum(mass_drop_times))
        push!(soft_drop_sample, minimum(soft_drop_times))
        push!(filter_sample, minimum(filtering))
        push!(trim_sample, minimum(trimming))
    end
    hepmc3_files_df[!, :filter_time] = filter_sample
    hepmc3_files_df[!, :trimming_time] = trim_sample
    hepmc3_files_df[!, :massdrop_time] = mass_drop_sample
    hepmc3_files_df[!, :softdrop_time] = soft_drop_sample
    
    if !isnothing(args[:results])
        CSV.write(args[:results], hepmc3_files_df)
    end
    
    # Now calculate the global average for the 4 substructure functions
    # @infiltrate
    gavg_filter = sum(hepmc3_files_df[:, :ExJets] .* hepmc3_files_df[:, :filter_time]) / sum(hepmc3_files_df[:, :ExJets])
    gavg_trim = sum(hepmc3_files_df[:, :ExJets] .* hepmc3_files_df[:, :trimming_time]) / sum(hepmc3_files_df[:, :ExJets])
    gavg_softd = sum(hepmc3_files_df[:, :ExJets] .* hepmc3_files_df[:, :softdrop_time]) / sum(hepmc3_files_df[:, :ExJets])
    gavg_massd = sum(hepmc3_files_df[:, :ExJets] .* hepmc3_files_df[:, :massdrop_time]) / sum(hepmc3_files_df[:, :ExJets])

    # Construct a new dataframe with the summary values
    summary_df = DataFrame("Method" => ["Filter", "Trim", "Massdrop", "Softdrop"],
                            "Time" => [gavg_filter, gavg_trim, gavg_massd, gavg_softd])
    println(summary_df)
    if !isnothing(args[:results])
        CSV.write("summary-" * args[:results], summary_df)
    end
end

main()
