##### Remember to modify the settings file to:
##### Fix threads for Gurobi Master to 1
##### Fix threads for Gurobi subproblem to 1

case = dirname(@__FILE__);

import Pkg

Pkg.activate("../../../GenX-Benders")
include("../../src/GenX-Benders.jl")

genx_settings = get_settings_path(case, "genx_settings.yml") # Settings YAML file path
mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

using Distributed
addprocs(5)
@everywhere begin
    import Pkg
    Pkg.activate("../../../GenX-Benders")
end
println("Number of procs: ", nprocs())
println("Number of workers: ", nworkers())

@everywhere include("../../src/GenX-Benders.jl")

@time myRandArray, myTechTypes, mga_summary,Results_df= run_benders_mga(case)
