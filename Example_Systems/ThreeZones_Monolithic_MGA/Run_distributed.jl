case = dirname(@__FILE__);

import Pkg


ENV["GRB_LICENSE_FILE"] = "/scratch/gpfs/ryzhu/gurobi.lic"
ENV["SLURM_CPUS_PER_TASK"] = 5
ENV["JULIA_NUM_THREADS"] = 5
Pkg.activate("../../../GenX-Benders")
include("../../../GenX-Benders/src/GenX-Benders.jl")

using Distributed,ClusterManagers
#ntasks = parse(Int, ENV["SLURM_NTASKS"]);
cpus_per_task = parse(Int, ENV["SLURM_CPUS_PER_TASK"]);
addprocs(cpus_per_task, exeflags=["--threads=5", "--project=$(Base.active_project())"])

@everywhere begin
    import Pkg
    Pkg.activate("../../../GenX-Benders/")
end

println("Number of procs: ", nprocs())
println("Number of workers: ", nworkers())
for i in workers()
    id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
    println(id, " " , pid, " ", host)
end
@everywhere include("../../../GenX-Benders/src/GenX-Benders.jl")

@time run_genx_case!(case)