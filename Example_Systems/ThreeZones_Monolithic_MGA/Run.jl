case = dirname(@__FILE__);

import Pkg
println("Activating Package")

Pkg.activate("../../../GenX-Benders")
include("../../../GenX-Benders/src/GenX-Benders.jl")

println("Success")

@time myRandArray, myTechTypes, mga_summary,Results_df= run_benders_mga(case)