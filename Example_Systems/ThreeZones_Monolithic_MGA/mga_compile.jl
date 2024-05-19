using CSV
using DataFrames

file_loc() = "C:\\Users\\mike_\\Documents\\ZeroLab\\GenX-Benders\\Example_Systems\\ThreeZones_Benders_MGA"

function compile()
    max = "MGAResults_max"
    min = "MGAResults_min"
    capacities = Array{Float64,2}(undef,(0,11))
    costs = []
    for folder in readdir(file_loc(), join = false)
        if folder == max || folder == min
            for case in readdir(joinpath(file_loc(),folder), join=true)
                for file in readdir(case,join=false)
                    name = joinpath(case,file)
                    if file == "capacity.csv"
                        df = DataFrame(CSV.File(name))
                        capacities = [capacities;df.EndCap']
                    elseif file == "costs.csv"
                        df = DataFrame(CSV.File(name))
                        costs = [costs;df.Total[1]]
                    end
                end
            end
        end
    end
    df_summary = DataFrame(capacities,:auto)
    df_summary[!,:cost] = costs
    outpath = "C:\\Users\\mike_\\Documents\\ZeroLab\\GenX-Benders\\Example_Systems\\ThreeZones_Benders_MGA\\Outputs\\Summary_mono.csv"
    CSV.write(outpath,df_summary)
    return 
end

compile()