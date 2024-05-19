@doc raw"""
	mga(EP::Model, path::AbstractString, setup::Dict, inputs::Dict, outpath::AbstractString)

We have implemented an updated Modeling to Generate Alternatives (MGA) Algorithm proposed by [Berntsen and Trutnevyte (2017)](https://www.sciencedirect.com/science/article/pii/S0360544217304097) to generate a set of feasible, near cost-optimal technology portfolios. This algorithm was developed by [Brill Jr, E. D., 1979](https://pubsonline.informs.org/doi/abs/10.1287/mnsc.25.5.413) and introduced to energy system planning by [DeCarolia, J. F., 2011](https://www.sciencedirect.com/science/article/pii/S0140988310000721).

To create the MGA formulation, we replace the cost-minimizing objective function of GenX with a new objective function that creates multiple generation portfolios by zone. We further add a new budget constraint based on the optimal objective function value $f^*$ of the least-cost model and the user-specified value of slack $\delta$. After adding the slack constraint, the resulting MGA formulation is given as:

```math
\begin{aligned}
	\text{max/min} \quad
	&\sum_{z \in \mathcal{Z}}\sum_{r \in \mathcal{R}} \beta_{z,r}^{k}P_{z,r}\\
	\text{s.t.} \quad
	&P_{zr} = \sum_{y \in \mathcal{G}}\sum_{t \in \mathcal{T}} \omega_{t} \Theta_{y,t,z,r}  \\
	& f \leq f^* + \delta \\
	&Ax = b
\end{aligned}
```

where, $\beta_{zr}$ is a random objective fucntion coefficient betwen $[0,100]$ for MGA iteration $k$. $\Theta_{y,t,z,r}$ is a generation of technology $y$ in zone $z$ in time period $t$ that belongs to a resource type $r$. We aggregate $\Theta_{y,t,z,r}$ into a new variable $P_{z,r}$ that represents total generation from technology type $r$ in a zone $z$. In the second constraint above, $\delta$ denote the increase in budget from the least-cost solution and $f$ represents the expression for the total system cost. The constraint $Ax = b$ represents all other constraints in the power system model. We then solve the formulation with minimization and maximization objective function to explore near optimal solution space.
"""
function mga(EP::Model, path::AbstractString, setup::Dict, inputs::Dict, outpath::AbstractString)

    if setup["ModelingToGenerateAlternatives"]==1
        # Start MGA Algorithm
	    println("MGA Module")

	    # Objective function value of the least cost problem
	    Least_System_Cost = objective_value(EP)

	    # Read sets
	    dfGen = inputs["dfGen"]
	    T = inputs["T"]     # Number of time steps (hours)
	    Z = inputs["Z"]     # Number of zonests
	    G = inputs["G"]

	    # Create a set of unique technology types
	    TechTypes = unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])

	    # Read slack parameter representing desired increase in budget from the least cost solution
	    slack = setup["ModelingtoGenerateAlternativeSlack"]

	    ### Variables ###

	    #@variable(EP, vSumvP[TechTypes = 1:length(TechTypes), z = 1:Z] >= 0) # Variable denoting total generation from eligible technology of a given type

	    ### End Variables ###

	    ### Constraints ###

	    # Constraint to set budget for MGA iterations
	    @constraint(EP, budget, EP[:eObj] <= Least_System_Cost * (1 + slack) )

        # Constraint to compute total generation in each zone from a given Technology Type
	    #@constraint(EP,cGeneration[tt = 1:length(TechTypes), z = 1:Z], vSumvP[tt,z] == sum(EP[:vP][y,t] * inputs["omega"][t]
	    #for y in dfGen[(dfGen[!,:Resource_Type] .== TechTypes[tt]) .& (dfGen[!,:Zone] .== z), :R_ID], t in 1:T))

	    ### End Constraints ###
		#@expression(EP,eTotalCapByType[type in TechTypes,z=1:Z],sum(EP[:eTotalCap][y] for y in dfGen[(dfGen[!,:Resource_Type] .== type)  .& (dfGen[!,:Zone] .== z), :R_ID]))
		#println(EP[:eTotalCapByType])
	    ### Create Results Directory for MGA iterations
        outpath_max = joinpath(path, "MGAResults_max")
	    if !(isdir(outpath_max))
	    	mkdir(outpath_max)
	    end
        outpath_min = joinpath(path, "MGAResults_min")
	    if !(isdir(outpath_min))
	    	mkdir(outpath_min)
	    end

	    ### Begin MGA iterations for maximization and minimization objective ###
	    mga_start_time = time()

	    print("Starting the first MGA iteration")
		pRand = rand(length(unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])),length(unique(dfGen[!,:Zone])), setup["ModelingToGenerateAlternativeIterations"])

		#rand_vectors = Dict((TechTypes[i],z,k)=> pRand[i,z,k] for i in eachindex(TechTypes), z in 1:inputs["Z"], k in 1:setup["ModelingToGenerateAlternativeIterations"]);
		
		EP_copy = Vector{Model}(undef,Threads.nthreads())
		settings_path = joinpath(path, "Settings")
		OPTIMIZER = configure_solver("gurobi", settings_path)
		for i in 1:length(EP_copy)
		    EP_copy[i] = copy(EP)
		    set_optimizer(EP_copy[i], OPTIMIZER)
		end
		
		Threads.@threads for i in 1:setup["ModelingToGenerateAlternativeIterations"]
	    	# Create random coefficients for the generators that we want to include in the MGA run for the given budget
	    	#pRand = rand(length(unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])),length(unique(dfGen[!,:Zone])))

	    	### Maximization objective
			println(pRand)
	
	    	@objective(EP_copy[i], Max, sum(pRand[tt,z,i] * EP_copy[i][:vSumvCap][tt,z] for tt in 1:length(TechTypes), z in 1:Z ))

	    	# Solve Model Iteration
	    	status = optimize!(EP_copy[i])

            # Create path for saving MGA iterations
	    	mgaoutpath_max = joinpath(outpath_max, string("MGA", "_", slack,"_", i))

	    	# Write results
	    	write_outputs(EP_copy[i], mgaoutpath_max, setup, inputs)

	    	### Minimization objective 
	    	@objective(EP_copy[i], Min, sum(pRand[tt,z,i] * EP_copy[i][:vSumvCap][tt,z] for tt in 1:length(TechTypes), z in 1:Z ))

	    	# Solve Model Iteration
	    	status = optimize!(EP_copy[i])

            # Create path for saving MGA iterations
	    	mgaoutpath_min = joinpath(outpath_min, string("MGA", "_", slack,"_", i))

	    	# Write results
	    	write_outputs(EP_copy[i], mgaoutpath_min, setup, inputs)

	    end

	    total_time = time() - mga_start_time
	    ### End MGA Iterations ###
	end

end

function mono_mga(EP::Model, path::AbstractString, setup::Dict, inputs::Dict, outpath::AbstractString)
		
	include(joinpath(path, "../../src/additional_tools/mono_mga_helpers.jl"))
    if setup["ModelingToGenerateAlternatives"]==1
        # Start MGA Algorithm
	    println("MGA Module")

	    # Objective function value of the least cost problem
	    Least_System_Cost = objective_value(EP)

	    # Read sets
	    dfGen = inputs["dfGen"]
	    T = inputs["T"]     # Number of time steps (hours)
	    Z = inputs["Z"]     # Number of zonests
	    G = inputs["G"]

	    # Create a set of unique technology types
	    TechTypes = unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])

	    # Read slack parameter representing desired increase in budget from the least cost solution
	    slack = setup["ModelingtoGenerateAlternativeSlack"]

	    ### Variables ###

	    #@variable(EP, vSumvP[TechTypes = 1:length(TechTypes), z = 1:Z] >= 0) # Variable denoting total generation from eligible technology of a given type

	    ### End Variables ###

	    ### Constraints ###

	    # Constraint to set budget for MGA iterations
	    @constraint(EP, budget, EP[:eObj] <= Least_System_Cost * (1 + slack) )

        # Constraint to compute total generation in each zone from a given Technology Type
	    #@constraint(EP,cGeneration[tt = 1:length(TechTypes), z = 1:Z], vSumvP[tt,z] == sum(EP[:vP][y,t] * inputs["omega"][t]
	    #for y in dfGen[(dfGen[!,:Resource_Type] .== TechTypes[tt]) .& (dfGen[!,:Zone] .== z), :R_ID], t in 1:T))

	    ### End Constraints ###
		#@expression(EP,eTotalCapByType[type in TechTypes,z=1:Z],sum(EP[:eTotalCap][y] for y in dfGen[(dfGen[!,:Resource_Type] .== type)  .& (dfGen[!,:Zone] .== z), :R_ID]))
		#println(EP[:eTotalCapByType])
	    ### Create Results Directory for MGA iterations
        outpath_max = joinpath(path, "MGAResults_max")
	    if !(isdir(outpath_max))
	    	mkdir(outpath_max)
	    end
        outpath_min = joinpath(path, "MGAResults_min")
	    if !(isdir(outpath_min))
	    	mkdir(outpath_min)
	    end

	    ### Begin MGA iterations for maximization and minimization objective ###
	    mga_start_time = time()

	    print("Starting the first MGA iteration")
		pRand = rand(length(unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])),length(unique(dfGen[!,:Zone])), setup["ModelingToGenerateAlternativeIterations"])

		#rand_vectors = Dict((TechTypes[i],z,k)=> pRand[i,z,k] for i in eachindex(TechTypes), z in 1:inputs["Z"], k in 1:setup["ModelingToGenerateAlternativeIterations"]);
		
		EP_copy = Vector{Model}(undef,setup["ModelingToGenerateAlternativeIterations"])

		settings_path = joinpath(path, "Settings")
		OPTIMIZER = configure_solver("gurobi", settings_path)
		for i in 1:length(EP_copy)
		    EP_copy[i] = copy(EP)
		    set_optimizer(EP_copy[i], OPTIMIZER)
		end

		println("breakpoint")
		
		Threads.@threads for i in 1:setup["ModelingToGenerateAlternativeIterations"]
	    	# Create random coefficients for the generators that we want to include in the MGA run for the given budget
	    	#pRand = rand(length(unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])),length(unique(dfGen[!,:Zone])))
			tid = Threads.threadid()
			print("Thread ID: $tid")
	    	### Maximization objective
			# println(pRand)
	
	    	@objective(EP_copy[i], Max, sum(pRand[tt,z,i] * EP_copy[i][:vSumvCap][tt,z] for tt in 1:length(TechTypes), z in 1:Z ))

	    	# Solve Model Iteration
	    	status = optimize!(EP_copy[i])

            # Create path for saving MGA iterations
	    	mgaoutpath_max = joinpath(outpath_max, string("MGA", "_", slack,"_", i))

	    	# Write results
	    	write_outputs(EP_copy[i], mgaoutpath_max, setup, inputs)

	    	### Minimization objective 
	    	@objective(EP_copy[i], Min, sum(pRand[tt,z,i] * EP_copy[i][:vSumvCap][tt,z] for tt in 1:length(TechTypes), z in 1:Z ))

	    	# Solve Model Iteration
	    	status = optimize!(EP_copy[i])

            # Create path for saving MGA iterations
	    	mgaoutpath_min = joinpath(outpath_min, string("MGA", "_", slack,"_", i))

	    	# Write results
	    	write_outputs(EP_copy[i], mgaoutpath_min, setup, inputs)

	    end

	    total_time = time() - mga_start_time
	    ### End MGA Iterations ###
	end

end

function compare_mga(EP::Model, path::AbstractString, setup::Dict, inputs::Dict, outpath::AbstractString, pRand::AbstractArray)

    if setup["ModelingToGenerateAlternatives"]==1
        # Start MGA Algorithm
	    println("MGA Module")

	    # Objective function value of the least cost problem
	    Least_System_Cost = objective_value(EP)

	    # Read sets
	    dfGen = inputs["dfGen"]
	    T = inputs["T"]     # Number of time steps (hours)
	    Z = inputs["Z"]     # Number of zonests
	    G = inputs["G"]

	    # Create a set of unique technology types
	    TechTypes = unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])

	    # Read slack parameter representing desired increase in budget from the least cost solution
	    slack = setup["ModelingtoGenerateAlternativeSlack"]

	    ### Variables ###

	    #@variable(EP, vSumvP[TechTypes = 1:length(TechTypes), z = 1:Z] >= 0) # Variable denoting total generation from eligible technology of a given type

	    ### End Variables ###

	    ### Constraints ###

	    # Constraint to set budget for MGA iterations
	    @constraint(EP, budget, EP[:eObj] <= Least_System_Cost * (1 + slack) )

        # Constraint to compute total generation in each zone from a given Technology Type
	    #@constraint(EP,cGeneration[tt = 1:length(TechTypes), z = 1:Z], vSumvP[tt,z] == sum(EP[:vP][y,t] * inputs["omega"][t]
	    #for y in dfGen[(dfGen[!,:Resource_Type] .== TechTypes[tt]) .& (dfGen[!,:Zone] .== z), :R_ID], t in 1:T))

	    ### End Constraints ###
		@expression(EP,eTotalCapByType[type in TechTypes,z=1:Z],sum(EP[:eTotalCap][y] for y in dfGen[(dfGen[!,:Resource_Type] .== type)  .& (dfGen[!,:Zone] .== z), :R_ID]))
	    ### Create Results Directory for MGA iterations
        outpath_max = joinpath(path, "MGAResults_max")
	    if !(isdir(outpath_max))
	    	mkdir(outpath_max)
	    end
        outpath_min = joinpath(path, "MGAResults_min")
	    if !(isdir(outpath_min))
	    	mkdir(outpath_min)
	    end

	    ### Begin MGA iterations for maximization and minimization objective ###
	    mga_start_time = time()

	    print("Starting the first MGA iteration")
		"""
		pRand = Array{Float64,3}(undef,(3,3,3))
		pRand[:,:,1] = [0.9941433520374436 0.9722530809556771 0.005613401863890477; 0.5384442571726544 0.5416959163836843 0.10273178350840828; 0.596469687081784 0.17242019223228944 0.4670661081818075]
		pRand[:,:,2] = [0.6208070401227073 0.705847807622118 0.5037772809444974; 0.43198353057476746 0.11994611553124401 0.5701523148533167; 0.341542218512473 0.6519646218223503 0.39696854376708157]
		pRand[:,:,3] = [0.6632099672363165 0.48500172265450514 0.854403449992875; 0.4002080492200133 0.00018749765849035427 0.6935243148345092; 0.48195575476135954 0.1811854849371659 0.4251400756432784]
		"""
		rand_vectors = Dict((TechTypes[i],z,k)=> pRand[i,z,k] for i in eachindex(TechTypes), z in 1:inputs["Z"], k in 1:setup["ModelingToGenerateAlternativeIterations"]);
		for i in 1:setup["ModelingToGenerateAlternativeIterations"]

	    	# Create random coefficients for the generators that we want to include in the MGA run for the given budget
	    	#pRand = rand(length(unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])),length(unique(dfGen[!,:Zone])))

	    	### Maximization objective
			println(pRand)
			println(EP[:eTotalCapByType])
	    	@objective(EP, Max, sum(rand_vectors[type,z,i] * EP[:eTotalCapByType][type,z] for type in TechTypes, z in 1:Z ))

	    	# Solve Model Iteration
	    	status = optimize!(EP)

            # Create path for saving MGA iterations
	    	mgaoutpath_max = joinpath(outpath_max, string("MGA", "_", slack,"_", i))

	    	# Write results
	    	write_outputs(EP, mgaoutpath_max, setup, inputs)

	    	### Minimization objective 
	    	@objective(EP, Min, sum(rand_vectors[type,z,i] * EP[:eTotalCapByType][type,z] for type in TechTypes, z in 1:Z ))

	    	# Solve Model Iteration
	    	status = optimize!(EP)

            # Create path for saving MGA iterations
	    	mgaoutpath_min = joinpath(outpath_min, string("MGA", "_", slack,"_", i))

	    	# Write results
	    	write_outputs(EP, mgaoutpath_min, setup, inputs)

	    end

	    total_time = time() - mga_start_time
	    ### End MGA Iterations ###
	end

end

function hsj_obj(new_point::AbstractArray, indic::AbstractArray)
    nrow, ncol = size(new_point)
    for j in 1:ncol, k in 1:nrow
		if new_point[k, j] >= 0.01
			indic[k,j] += 1
		end
    end
    return indic
end

function test_ops(EP::Model,path::AbstractString, setup::Dict, inputs::Dict)
    Least_System_Cost = objective_value(EP)
    MASTER_OPTIMIZER =  optimizer_with_attributes(Gurobi.Optimizer,MOI.Silent() => true)
    
    EP_new = generate_model(setup, inputs, MASTER_OPTIMIZER)
    
    # Read sets
    dfGen = inputs["dfGen"]
    T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zonests
    G = inputs["G"]
    
    
	EP_m, m_vars, m_cons = generate_master_problem(setup,inputs,MASTER_OPTIMIZER)
	
	m_vars_mga = intersect(name.(all_variables(EP_new)),m_vars)

    # Create a set of unique technology types
    TechTypes = unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])

    # Read slack parameter representing desired increase in budget from the least cost solution
    slack = setup["ModelingtoGenerateAlternativeSlack"]

	# Setup storage
	Indic = fill(0, (length(TechTypes), Z))
	point_list = fill(0.0, (length(TechTypes), Z, setup["ModelingToGenerateAlternativeIterations"]+1))
	point_list[:,:,1] = value.(EP[:vSumvCap])
	
	
    # Constraint to set budget for MGA iterations
    @constraint(EP_new, budget, EP_new[:eObj] <= Least_System_Cost * (1 + slack) )


    ### End Constraints ###

    ### Create Results Directory for MGA iterations
    outpath_hsj = joinpath(path, "MGAResults")
    if !(isdir(outpath_hsj))
    	mkdir(outpath_hsj)
    end
    """
    println("In Min all Base")
    Indic = hsj_obj(point_list[:,:,1],Indic)
    @objective(EP_new, Min, sum(EP_new[:vSumvCap][tt,z] for tt in 1:length(TechTypes), z in 1:Z ))#
    # Solve Model Iteration
	status = optimize!(EP_new)
	mgaoutpath_hsj = joinpath(outpath_hsj, string("MGA", "_", slack,"_", "MinAllBase"))
	println("Solved HSJ Base")
	# Write results
	write_outputs(EP_new, mgaoutpath_hsj, setup, inputs)
	vals = Vector{Float64}(undef,0)
	println("Fixing")
	for y in m_vars_mga
	    vy = variable_by_name(EP_new,y)
	    push!(vals,value(vy))
	end
	for i in 1:length(m_vars_mga)
	    vy = variable_by_name(EP_new,m_vars_mga[i])
	    print(vy)
	    fix(vy, vals[i]; force=true)
	 end
	println("Solving min all Opt")
	@objective(EP_new, Min, EP_new[:eObj])
	optimize!(EP_new)
	mgaoutpath_hsj = joinpath(outpath_hsj, string("MGA", "_", slack,"_", "MinAllOpt"))

	# Write results
	write_outputs(EP_new, mgaoutpath_hsj, setup, inputs)
	println("Unfixing")
	for i in 1:length(m_vars_mga)
	    vy = variable_by_name(EP_new,m_vars_mga[i])
	    print(vy)
	    unfix(vy)
	 end
	"""
	pRand = rand(length(unique(dfGen[dfGen[!, :MGA] .== 1, :Resource_Type])),length(unique(dfGen[!,:Zone])))
	@expression(EP_new, e_vSumvP[tt = 1:length(TechTypes), z = 1:Z], sum(EP_new[:vP][y,t] * inputs["omega"][t]
	    for y in dfGen[(dfGen[!,:Resource_Type] .== TechTypes[tt]) .& (dfGen[!,:Zone] .== z), :R_ID], t in 1:T))
	println("Solving production max")
	@objective(EP_new, Max, sum(pRand[tt,z] * EP_new[:e_vSumvP][tt,z] for tt in 1:length(TechTypes), z in 1:Z ))
    # Solve Model Iteration
	status = optimize!(EP_new)
	mgaoutpath_hsj = joinpath(outpath_hsj, string("MGA", "_", slack,"_", "BaseProd"))
	# Write results
	write_outputs(EP_new, mgaoutpath_hsj, setup, inputs)
	vals = Vector{Float64}(undef,0)
	"""
	Generate model
	least cost optimization set up
	
	Array structure
	    Rows - MGA iterations
	    Columns - Capacity variables or resource ids
	Make sure you match the capacity variables in this to the actual model - it has to be right match    
	
	
	loop over rows
	    fix capacity variables to capacities from mga iteration in question
	    optimize!
	    write outputs
	    unfix capacity variables
    end loop
    
    
	
	"""
	println("Fixing prod cap")
	for y in m_vars_mga
	    vy = variable_by_name(EP_new,y)
	    push!(vals,value(vy))
	end
	for i in 1:length(m_vars_mga)
	    vy = variable_by_name(EP_new,m_vars_mga[i])
	    println(vy)
	    println(vals[i])
	    fix(vy, vals[i]; force=true)
	 end
	println("Solving prod Opt")
	@objective(EP_new, Min, EP_new[:eObj])
	optimize!(EP_new)
	mgaoutpath_hsj = joinpath(outpath_hsj, string("MGA", "_", slack,"_", "OptimalProd"))
	for i in 1:length(m_vars_mga)
	    vy = variable_by_name(EP_new,m_vars_mga[i])
	    println(vy)
	    println(vals[i])
	 end
	write_outputs(EP_new, mgaoutpath_hsj, setup, inputs)
	
end