
# function initialize_dist_helpers(setup::Dict,inputs_decomp::Dict,master_vars::Vector{String},master_cons::Vector{String})
function initialize_dist_helpers()
    ##### Initialize a distributed arrays of JuMP models
	## Start pre-solve timer
	subproblem_generation_time = time()
    
    num_procs = ENV["SLURM_CPUS_PER_TASK"]
    helpers_all = distribute([Dict() for i in 1:num_procs]);
	workers_all = workers()

    @sync for i in 1:num_procs
		p = workers_all[i]
        @async @spawnat p begin
            W_local = localindices(helpers_all)[1];
            inputs_local = [inputs_decomp[k] for k in W_local];
			SUBPROB_OPTIMIZER =  configure_benders_subprob_solver(setup["Solver"], setup["settings_path"]);
            init_local_helper!(setup,inputs_local,localpart(helpers_all),master_vars,master_cons,SUBPROB_OPTIMIZER);
        end
    end

	# master_vars_sub = Dict();
	# for i in eachindex(helpers_all)
	# 	w = helpers_all[i]["SubPeriod"];
	# 	master_vars_sub[w] = helpers_all[i]["master_vars_sub"];
	# end

	p_id = workers_all[1:num_sub];
    np_id = length(p_id);

    master_vars_sub = [Dict() for k in 1:np_id];

    @sync for k in 1:np_id
              @async master_vars_sub[k]= @fetchfrom p_id[k] get_local_master_vars(localpart(helpers_all))
    end

	master_vars_sub = merge(master_vars_sub...);

    ## Record pre-solver time
	subproblem_generation_time = time() - subproblem_generation_time
	println("Distributed operational subproblems generation took $subproblem_generation_time seconds")

    return helpers_all,master_vars_sub

end

function init_local_helper!(setup::Dict,inputs_local::Vector{Dict{Any,Any}},helper_local::Vector{Dict{Any,Any}},master_vars::Vector{String},master_cons::Vector{String},OPTIMIZER::MOI.OptimizerWithAttributes)

    nW = length(inputs_local)

    for i=1:nW
		EP, master_vars_sub = generate_subproblem(setup,inputs_local[i],master_vars,master_cons,OPTIMIZER);
        helper_local[i]["Model"] = EP;
        helper_local[i]["master_vars_sub"] = master_vars_sub
        helper_local[i]["SubPeriod"] = inputs_local[i]["SubPeriod"];
    end
end


function solve_dist_helpers(EP_helpers::DArray{Dict{Any, Any}, 1, Vector{Dict{Any, Any}}},master_sol::NamedTuple,inputs::Dict)

    p_id = workers();
    np_id = length(p_id);

    sub_results = [Dict() for k in 1:np_id];

    @sync for k in 1:np_id
              @async sub_results[k]= @fetchfrom p_id[k] solve_local_helper(localpart(EP_helpers),master_sol,inputs); ### This is equivalent to fetch(@spawnat p .....)
    end

	sub_results = merge(sub_results...);

    return sub_results
end
 
function solve_local_helper(helper_local::Vector{Dict{Any,Any}},master_sol::NamedTuple,inputs::Dict)

    local_sol=Dict();
    for m in helper_local
        EP = m["Model"];
        master_vars_sub = m["master_vars_sub"]
        w = m["SubPeriod"];
		local_sol[w] = solve_subproblem(EP,master_sol,master_vars_sub,inputs);
    end
    return local_sol
end