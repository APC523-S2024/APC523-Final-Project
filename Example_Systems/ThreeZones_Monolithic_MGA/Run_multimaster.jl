case = dirname(@__FILE__);

import Pkg


Pkg.activate("../../../GenX-Benders/")
include("../../../GenX-Benders/src/GenX-Benders.jl")

using Distributed,ClusterManagers
n_master = 2
n_threads_master = 8
n_sub = 5
n_threads_sub = 1
addprocs(n_sub;exeflags=["-t $n_threads_sub"])
addprocs(n_master;exeflags=["-t $n_threads_master"])

@everywhere begin
    import Pkg
    Pkg.activate("../../../GenX-Benders/")
end

println("Number of procs: ", nprocs())
println("Number of workers: ", nworkers())
for i in workers()
    id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
    println(id, " " , pid, " ", host, " ")
end
@everywhere include("../../../GenX-Benders/src/GenX-Benders.jl")

@time myRandArray, myTechTypes, mga_summary,Results_df= benders_mga_multi_master(case)