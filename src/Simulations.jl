###################################################################################################
# SIMULATION FUNCTIONS
###################################################################################################
"""
    `gen_x_mc(N,α,β)`

    Function that generates a methylation vector from an Ising model with `[N1,...,NK]`, and parameters
    `[α1,...,αK]` and `β`.

    # Examples
    ```julia-repl
    julia> using Random; Random.seed!(1234);
    julia> CpelTdm.gen_x_mc([5],[0.0,0.0])
    5-element Array{Int64,1}:
     1
     1
     1
     -1
     1
    ```
"""
function gen_x_mc(n::Vector{Int64},θ::Vector{Float64})::Vector{Int64}

    # Ge parameters
    α = θ[1:(end-1)]
    β = θ[end]

    # Find subregion label for each CpG site
    subid = vcat([i*ones(Int64,n[i]) for i in 1:length(n)]...)

    # Sample first CpG site
    x = Vector{Int64}()
    p = comp_lkhd(pushfirst!(zeros(Int64,sum(n)-1),1),n,α,β)
    push!(x,sample([-1,1],Distributions.Weights([1-p,p])))
    sum(n)>1 || return x

    # Sequentially sample from Markov chain
    @inbounds for i=2:sum(n)

        # Add boundary function g() if necessary
        expaux = exp(α[subid[i]]+β*x[i-1])
        if i < sum(n)
            ap1p = 2.0*β + α[subid[i]]
            ap1m = -2.0*β + α[subid[i]]
            n_miss = [count(x->x==sr,subid[(i+1):sum(n)]) for sr in unique(subid[(i+1):sum(n)])]
            gp = comp_g(n_miss,α[unique(subid[(i+1):sum(n)])],β,ap1p,α[end])
            gm = comp_g(n_miss,α[unique(subid[(i+1):sum(n)])],β,ap1m,α[end])
            p = gp*expaux/(gp*expaux+gm/expaux)
        else
            p = expaux/(expaux+1.0/expaux)
        end

        # Add i-th CpG site
        push!(x,sample([-1,1],Distributions.Weights([1-p,p])))

    end

    # Return methylation vector
    return x

end # end gen_x_mc
"""
    `unmatched_sim(m,n,θ1,θ2)`

    Function that simulates hypothesis testing in a (balanced) unmatched group comparison.

    # Examples
    ```julia-repl
    julia> m=5; n=[5]; θ1=[0.5,0.5,0.25]; θ2=[0.5,0.5,0.25];
    julia> tmml,tnme,tpdm = CpelTdm.unmatched_sim(m,n,θ1,θ2)
    ```
"""
function unmatched_sim(m::Int64,n::Vector{Int64},θ1::Vector{Float64},θ2::Vector{Float64})::NTuple{3,Vector{NTuple{2,Float64}}}

    # Generate 1k data-sets
    tmml = Vector{NTuple{2,Float64}}()
    tnme = Vector{NTuple{2,Float64}}()
    tpdm = Vector{NTuple{2,Float64}}()
    while length(tmml)<1000
        
        ## Generate data
        
        # g1
        g1_samples = []
        for i=1:m
            push!(g1_samples,[gen_x_mc(n,θ1) for j=1:rand(10:30)])
        end

        # g2
        g2_samples = []
        for i=1:m
            push!(g2_samples,[gen_x_mc(n,θ2) for j=1:rand(10:30)])
        end

        ## Estimate vectors

        # g1
        g1_θs = Vector{Vector{Float64}}()
        for xobs in g1_samples
            push!(g1_θs,est_theta_sa(n,xobs))
        end

        # g2
        g2_θs = Vector{Vector{Float64}}()
        for xobs in g2_samples
            push!(g2_θs,est_theta_sa(n,xobs))
        end

        ## Compute unmatched pvals
        tmml_data,tnme_data,tpdm_data = unmat_tests(n,g1_θs,g2_θs)

        # Append
        push!(tmml,tmml_data)
        push!(tnme,tnme_data)
        push!(tpdm,tpdm_data)

    end
    
    # Return p-values
    return tmml,tnme,tpdm

end
"""
    `matched_sim(m,n,θ1,θ2)`

    Function that simulates hypothesis testing in a matched group comparison.

    # Examples
    ```julia-repl
    julia> m=5; n=[5]; θ1=[0.5,0.5,0.25]; θ2=[0.5,0.5,0.25];
    julia> tmml,tnme,tpdm = CpelTdm.matched_sim(m,n,θ1,θ2)
    ```
"""
function matched_sim(m::Int64,n::Vector{Int64},θ1s::Vector{Vector{Float64}},θ2s::Vector{Vector{Float64}})::NTuple{3,Vector{NTuple{2,Float64}}}

    # Generate 1k data-sets
    tmml = Vector{NTuple{2,Float64}}()
    tnme = Vector{NTuple{2,Float64}}()
    tpdm = Vector{NTuple{2,Float64}}()
    while length(tmml)<1000
        
        ## Generate data

        # g1
        g1_samples = []
        for i=1:m
            push!(g1_samples,[gen_x_mc(n,θ1s[i]) for j=1:rand(10:30)])
        end

        # g2
        g2_samples = []
        for i=1:m
            push!(g2_samples,[gen_x_mc(n,θ2s[i]) for j=1:rand(10:30)])
        end

        ## Estimate vectors

        # g1
        g1_θs = Vector{Vector{Float64}}()
        for xobs in g1_samples
            push!(g1_θs,est_theta_sa(n,xobs))
        end

        # g2
        g2_θs = Vector{Vector{Float64}}()
        for xobs in g2_samples
            push!(g2_θs,est_theta_sa(n,xobs))
        end

        ## Compute unmatched pvals
        tmml_data,tnme_data,tpdm_data = mat_tests(n,g1_θs,g2_θs)

        # Append
        push!(tmml,tmml_data)
        push!(tnme,tnme_data)
        push!(tpdm,tpdm_data)

    end
    
    # Return p-values
    return tmml,tnme,tpdm

end