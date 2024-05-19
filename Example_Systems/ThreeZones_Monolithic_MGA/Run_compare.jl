case = dirname(@__FILE__);

import Pkg
using DataFrames
using CSV
Pkg.activate("../../../GenX-Benders")
include("../../../GenX-Benders/src/GenX-Benders.jl")
# using Distributed
# addprocs(5)
# @everywhere begin
#     import Pkg
#     # Pkg.activate("/home/fp0820/GenX-Benders")
#     Pkg.activate("C:/Users/pecci/Code/GenX-Benders/")
# end

# println("Number of procs: ", nprocs())
# println("Number of workers: ", nworkers())

# # @everywhere include("/home/fp0820/GenX-Benders/src/GenX-Benders.jl")
# @everywhere include("C:/Users/pecci/Code/GenX-Benders/src/GenX-Benders.jl")
function run_cost_optimal_monolithic(case)
    genx_settings = get_settings_path(case, "genx_settings.yml") #Settings YAML file path
    mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters
    settings_path = get_settings_path(case)

    ### Cluster time series inputs if necessary and if specified by the user
    TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])

    if mysetup["TimeDomainReduction"] == 1
        prevent_doubled_timedomainreduction(case)
        println("Clustering Time Series Data (Grouped)...")
        if !time_domain_reduced_files_exist(TDRpath)
            if mysetup["NumScenarios"]<=1
                cluster_inputs(case, settings_path, mysetup)
            else
                unpack_scenarios(case,mysetup) 
                for s in 1:mysetup["NumScenarios"]
                    scenario_path = joinpath(joinpath(case,"Scenarios"),"Scenario_$s");
                    cluster_inputs(scenario_path, settings_path, mysetup);
                end
                pack_scenarios(case,mysetup)
            end
        else
            println("Time Series Data Already Clustered.")
        end
    end

    ### Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)

    #### Running a case

    ### Load inputs
    println("Loading Inputs")
    myinputs = load_inputs(mysetup, case)

    println("Generating the Optimization Model")
    time_elapsed = @elapsed EP = generate_model(mysetup, myinputs, OPTIMIZER)
    println("Time elapsed for model building is")
    println(time_elapsed)

    println("Solving Model")
    EP, solve_time = solve_model(EP, mysetup)
    myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

    # Run MGA if the MGA flag is set to 1 else only save the least cost solution
    println("Writing Output")
    outputs_path = get_default_output_folder(case)
    elapsed_time = @elapsed write_outputs(EP, outputs_path, mysetup, myinputs)
    println("Time elapsed for writing is")
    println(elapsed_time)
    lc_cost = value(EP[:eObj])
    return EP, outputs_path, mysetup, myinputs, lc_cost
end

function compare_mga_runner()
    EP, outputs_path, mysetup, myinputs, budget = run_cost_optimal_monolithic(case)

    @time capacities, myRandArray, myTechTypes, mga_summary,costs = compare_benders_mga(case, budget)
    println(myRandArray)
    println(capacities', myRandArray, myTechTypes, mga_summary)
    @time GenX.compare_mga(EP, case, mysetup, myinputs, outputs_path, myRandArray)

    outputs_df = DataFrame(capacities'.*10^3,:auto)
    outputs_df[!,:cost] = costs*10^6
    outpath = "../../../GenX-Benders/Example_Systems/ThreeZones_Benders_MGA/Outputs/Summary.csv"
    CSV.write(outpath,outputs_df)
end

function test_operations()
    EP, outputs_path, mysetup, myinputs, budget = run_cost_optimal_monolithic(case)
    test_ops(EP,outputs_path,mysetup,myinputs)
    
end

test_operations()