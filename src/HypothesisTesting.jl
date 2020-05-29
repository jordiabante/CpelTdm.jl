###################################################################################################
# COMMON
###################################################################################################
"""
    `get_all_∇logZs(n,θs)`

    Function that returns gradient vectors for all parameter vectors θ passed.

    # Examples
    ```julia-repl
    julia> n=[5]; θs=[[0.0,0.0],[0.0,0.0]];
    julia> CpelRrbs.get_all_∇logZs(n,θs)
    [[0.0,0.0],[0.0,0.0]]
    ```
"""
function get_all_∇logZs(n::Vector{Int64},θs::Vector{Vector{Float64}})::Vector{Vector{Float64}}

    # Generate all gradients
    ∇logZs = []
    @inbounds for s=1:length(θs)
        push!(∇logZs,get_∇logZ(n,θs[s]))
    end

    # Return vector of gradients
    return ∇logZs

end
###################################################################################################
# UNMATCHED TESTS
###################################################################################################
"""
    `comp_unmat_stat_mml(n,∇logZ1s,∇logZ2s)`

    Function that computes MML difference between groups (unmatched case).

    # Examples
    ```julia-repl
    julia> n=[5]; θ1s=[[0.0,0.0],[0.0,0.0]]; θ2s=[[1.0,1.0],[1.0,1.0]];
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> CpelRrbs.comp_unmat_stat_mml(n,∇logZ1s,∇logZ2s)
    -0.48867352
    ```
"""
function comp_unmat_stat_mml(n::Vector{Int64},∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}})::Float64
    
    # Compute mean for group 1
    mml1 = 0.0
    @inbounds for ∇logZ in ∇logZ1s
        mml1 += comp_mml(n,∇logZ)
    end
    mml1 /= length(∇logZ1s)

    # Compute mean for group 2
    mml2 = 0.0
    @inbounds for ∇logZ in ∇logZ2s
        mml2 += comp_mml(n,∇logZ)
    end
    mml2 /= length(∇logZ2s)

    # Return Tmml
    return mml1-mml2

end
"""
    `comp_unmat_perm_stat_mml(n,∇logZ1s,∇logZ2s,perm_ids)`

    Function that produces a permutation statistic for Tmml in unmatched case.

    # Examples
    ```julia-repl
    julia> n=[10]; θ1s=fill([0.5,0.5],5); θ2s=fill([-0.5,0.5],5); perm_ids=[1,3,5,7,9]; 
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> CpelRrbs.comp_unmat_perm_stat_mml(n,∇logZ1s,∇logZ2s,perm_ids)
    0.15564650799999996
    ```
"""
function comp_unmat_perm_stat_mml(n::Vector{Int64},∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}},perm_ids::Vector{Int64})::Float64

    # Get vectors for each group
    ∇logZ2sp = vcat(∇logZ1s,∇logZ2s)
    ∇logZ1sp = ∇logZ2sp[perm_ids]
    deleteat!(∇logZ2sp,perm_ids)

    # Return
    return comp_unmat_stat_mml(n,∇logZ1sp,∇logZ2sp)
    
end
"""
    `comp_unmat_stat_nme(n,θ1s,θ2s,∇logZ1s,∇logZ2s)`

    Function that computes NME difference between groups (unmatched case).

    # Examples
    ```julia-repl
    julia> n=[5]; θ1s=[[0.0,0.0],[0.0,0.0]]; θ2s=[[1.0,1.0],[1.0,1.0]];
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> CpelRrbs.comp_unmat_stat_nme(n,θ1s,θ2s,∇logZ1s,∇logZ2s)
    0.91986634
    ```
"""
function comp_unmat_stat_nme(n::Vector{Int64},θ1s::Vector{Vector{Float64}},θ2s::Vector{Vector{Float64}},
    ∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}})::Float64
    
    # Compute mean for group 1
    nme1 = 0.0
    @inbounds for s=1:length(θ1s)
        nme1 += comp_nme(n,θ1s[s],∇logZ1s[s])
    end
    nme1 /= length(θ1s)

    # Compute mean for group 2
    nme2 = 0.0
    @inbounds for s=1:length(θ2s)
        nme2 += comp_nme(n,θ2s[s],∇logZ2s[s])
    end
    nme2 /= length(θ2s)

    # Return difference in nme
    return nme1-nme2

end
"""
    `comp_unmat_perm_stat_nme(n,θ1s,θ2s,perm_ids)`

    Function that produces a permutation statistic for Tnme in unmatched case.

    # Examples
    ```julia-repl
    julia> n=[10]; θ1s=fill([0.5,0.5],5); θ2s=fill([-0.5,0.5],5); perm_ids=[1,3,5,7,9];
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> CpelRrbs.comp_unmat_perm_stat_nme(n,θ1s,θ2s,∇logZ1s,∇logZ2s,perm_ids)
    0.0
    ```
"""
function comp_unmat_perm_stat_nme(n::Vector{Int64},θ1s::Vector{Vector{Float64}},θ2s::Vector{Vector{Float64}},
    ∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}},perm_ids::Vector{Int64})::Float64

    # Get vectors for each group
    θ2sp = vcat(θ1s,θ2s)
    θ1sp = θ2sp[perm_ids]
    deleteat!(θ2sp,perm_ids)

    # Get vectors for each group
    ∇logZ2sp = vcat(∇logZ1s,∇logZ2s)
    ∇logZ1sp = ∇logZ2sp[perm_ids]
    deleteat!(∇logZ2sp,perm_ids)    

    # Return
    return comp_unmat_stat_nme(n,θ1sp,θ2sp,∇logZ1sp,∇logZ2sp)
    
end
"""
    `comp_unmat_stat_cmd(n,θ1s,θ2s,∇logZ1s,∇logZ2s)`

    Function that computes test statistic Tgsd in unmatched case.

    # Examples
    ```julia-repl
    julia> n=[5]; θ1s=[[-1.0,1.0],[-1.0,1.0]]; θ2s=[[1.0,1.0],[1.0,1.0]];
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> CpelRrbs.comp_unmat_stat_cmd(n,θ1s,θ2s,∇logZ1s,∇logZ2s)
    0.7888652058295635
    ```
"""
function comp_unmat_stat_cmd(n::Vector{Int64},θ1s::Vector{Vector{Float64}},θ2s::Vector{Vector{Float64}},
    ∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}})::Float64
    
    ## Compute stat
    cmds = []
    @inbounds for s1=1:length(θ1s)
        @inbounds for s2=1:length(θ2s)
            push!(cmds,comp_cmd(n,θ1s[s1],θ2s[s2],∇logZ1s[s1],∇logZ2s[s2]))
        end
    end

    # Return sum of GJSD
    return sum(cmds)/(length(θ1s)*length(θ2s))

end
"""
    `comp_unmat_perm_stat_cmd(n,θ1s,θ2s,∇logZ1s,∇logZ2s,perm_ids)`

    Function that produces a permutation statistic for Tcmd in unmatched case.

    # Examples
    ```julia-repl
    julia> n=[10]; θ1s=fill([0.5,0.5],5); θ2s=fill([-0.5,0.5],5); perm_ids=[1,3,5,7,9];
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> CpelRrbs.comp_unmat_perm_stat_cmd(n,θ1s,θ2s,∇logZ1s,∇logZ2s,perm_ids)
    0.1691594717040136
    ```
"""
function comp_unmat_perm_stat_cmd(n::Vector{Int64},θ1s::Vector{Vector{Float64}},θ2s::Vector{Vector{Float64}},
    ∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}},perm_ids::Vector{Int64})::Float64

    # Get vectors for each group
    θ2sp = vcat(θ1s,θ2s)
    θ1sp = θ2sp[perm_ids]
    deleteat!(θ2sp,perm_ids)

    # Get vectors for each group
    ∇logZ2sp = vcat(∇logZ1s,∇logZ2s)
    ∇logZ1sp = ∇logZ2sp[perm_ids]
    deleteat!(∇logZ2sp,perm_ids)  

    # Return
    return comp_unmat_stat_cmd(n,θ1sp,θ2sp,∇logZ1s,∇logZ2s)
    
end
"""
    `unmat_tests(n,θ1s,θ2s)`

    Function that performs hypothesis testing in unmatched samples group comparison.

    # Examples
    ```julia-repl
    julia> using Random; Random.seed!(1234);
    julia> n=[10]; θ1s=fill([0.5,0.5],5); θ2s=fill([-0.5,0.5],5);
    julia> tmml_test,tnme_test,tcmd_test = CpelRrbs.unmat_tests(n,θ1s,θ2s)
    ((0.7782325400000001, 0.007936507936507936), (0.0, 1.0), (0.3253066763538721, 0.003968253968253968))
    julia> n=[10]; θ1s=fill([0.5,0.5],8); θ2s=fill([-0.5,0.5],8);
    julia> tmml_test,tnme_test,tcmd_test = CpelRrbs.unmat_tests(n,θ1s,θ2s)
    ((0.77823254, 0.0020833333333333333), (0.0, 1.0), (0.325306676353872, 0.0010416666666666667))
    julia> n=[10]; θ1s=fill([0.0,0.0],8); θ2s=fill([0.0,2.5],8);
    julia> tmml_test,tnme_test,tcmd_test = CpelRrbs.unmat_tests(n,θ1s,θ2s)
    ((0.0, 1.0), (0.84782978, 0.0010309278350515464), (0.43786215630140907, 0.0010309278350515464))
    ```
"""
function unmat_tests(n::Vector{Int64},θ1s::Vector{Vector{Float64}},θ2s::Vector{Vector{Float64}};Lmax::Int64=1000)::NTuple{3,NTuple{2,Float64}}

    # Compute gradients
    ∇logZ1s = get_all_∇logZs(n,θ1s)
    ∇logZ2s = get_all_∇logZs(n,θ2s)

    # Compute observed stats
    tmml_obs = comp_unmat_stat_mml(n,∇logZ1s,∇logZ2s)
    tnme_obs = comp_unmat_stat_nme(n,θ1s,θ2s,∇logZ1s,∇logZ2s)
    tcmd_obs = comp_unmat_stat_cmd(n,θ1s,θ2s,∇logZ1s,∇logZ2s)

    # Compute number of possible randomizations
    L = binomial(length(θ1s)+length(θ2s),length(θ1s))
    exact = L<=Lmax

    # Create iteratable object with all combinations
    comb_iter = combinations(1:(length(θ1s)+length(θ2s)),length(θ1s))

    # Get group label combinations to use
    comb_iter_used = []
    if exact
        # Use all group assignments
        comb_iter_used = comb_iter
    else
        # Use Lmax group assignments
        ind_subset = rand(1:L,Lmax)
        @inbounds for (ind,comb) in enumerate(comb_iter)
            (ind in ind_subset) && push!(comb_iter_used,comb)
        end
    end

    # Use method for random permutation
    tmml_perms = map(x->comp_unmat_perm_stat_mml(n,∇logZ1s,∇logZ2s,x),comb_iter_used)
    tnme_perms = map(x->comp_unmat_perm_stat_nme(n,θ1s,θ2s,∇logZ1s,∇logZ2s,x),comb_iter_used)
    tcmd_perms = map(x->comp_unmat_perm_stat_cmd(n,θ1s,θ2s,∇logZ1s,∇logZ2s,x),comb_iter_used)

    # Compute p-values two-sided test
    if exact
        # Compute exact p-value
        tmml_pval = sum(abs.(tmml_perms).>=abs(tmml_obs))/length(tmml_perms)
        tnme_pval = sum(abs.(tnme_perms).>=abs(tnme_obs))/length(tnme_perms)
        tcmd_pval = sum(tcmd_perms.>=tcmd_obs)/length(tcmd_perms)
    else
        # Compute p-value in MC setting
        tmml_pval = (1.0+sum(abs.(tmml_perms).>=abs(tmml_obs)))/(1.0+length(tmml_perms))
        tnme_pval = (1.0+sum(abs.(tnme_perms).>=abs(tnme_obs)))/(1.0+length(tnme_perms))
        tcmd_pval = (1.0+sum(tcmd_perms.>=tcmd_obs))/(1.0+length(tcmd_perms))
    end

    # Return stat-pval pairs
    return (tmml_obs,tmml_pval),(tnme_obs,tnme_pval),(tcmd_obs,tcmd_pval)
    
end
###################################################################################################
# MATCHED TESTS
###################################################################################################
"""
    `comp_mat_diff_mml(n,∇logZ1s,∇logZ2s)`

    Function that computes MML differences between pairs (matched case).

    # Examples
    ```julia-repl
    julia> n=[5]; θ1s=[[0.0,0.0],[0.0,0.0]]; θ2s=[[1.0,1.0],[1.0,1.0]];
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> CpelRrbs.comp_mat_diff_mml(n,∇logZ1s,∇logZ2s)
    2-element Array{Float64,1}:
     -0.48867352
     -0.48867352
    ```
"""
function comp_mat_diff_mml(n::Vector{Int64},∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}})::Vector{Float64}
    
    # Compute mean differences
    diffs = []
    @inbounds for s=1:length(∇logZ1s)
        push!(diffs,comp_mml(n,∇logZ1s[s])-comp_mml(n,∇logZ2s[s]))
    end

    # Return vector of differences
    return diffs

end
"""
    `comp_mat_diff_nme(n,θ1s,θ2s,∇logZ1s,∇logZ2s)`

    Function that computes NME differences between pairs (matched case).

    # Examples
    ```julia-repl
    julia> n=[5]; θ1s=[[0.0,0.0],[0.0,0.0]]; θ2s=[[1.0,1.0],[1.0,1.0]];
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> CpelRrbs.comp_mat_diff_nme(n,θ1s,θ2s,∇logZ1s,∇logZ2s)
    2-element Array{Float64,1}:
     0.91986634
     0.91986634
    ```
"""
function comp_mat_diff_nme(n::Vector{Int64},θ1s::Vector{Vector{Float64}},θ2s::Vector{Vector{Float64}},
    ∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}})::Vector{Float64}
    
    # Compute entropy differences
    diffs = []
    @inbounds for s=1:length(θ1s)
        push!(diffs,comp_nme(n,θ1s[s],∇logZ1s[s])-comp_nme(n,θ2s[s],∇logZ2s[s]))
    end

    # Return vector of differences
    return diffs

end
"""
    `comp_mat_j_stat(diffs,j)`

    Function that computes permutation statistic for j-th sign assigment given vector of differences
    between matched pairs (matched case).

    # Examples
    ```julia-repl
    julia> n=[10]; θ1s=fill([0.5,0.5],5); θ2s=fill([-0.5,0.5],5);
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> mml_diffs = CpelRrbs.comp_mat_diff_mml(n,∇logZ1s,∇logZ2s)
    julia> CpelRrbs.comp_mat_j_stat(mml_diffs,0)
    -0.7782325400000001
    ```
"""
function comp_mat_j_stat(diffs::Vector{Float64},j::Int64)::Float64
    
    # Init diff in permutation
    diffs_perm = 0.0
    
    # Change sign if pertinent
    @inbounds for i=1:length(diffs)
        diffs_perm += Bool(j & 1) ? diffs[i] : -diffs[i]
        j >>= 1
    end

    # Return permutation statistic
    return diffs_perm/length(diffs)

end
"""
    `comp_mat_perm_stats(n,θ1s,θ2s)`

    Function that computes permutation statistics given a vector of differences between pairs (matched case).

    # Examples
    ```julia-repl
    julia> n=[10]; θ1s=fill([0.5,0.5],2); θ2s=fill([-0.5,0.5],2);
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> mml_diffs = CpelRrbs.comp_mat_diff_mml(n,∇logZ1s,∇logZ2s);
    julia> js = collect(0:2^length(θ1s)-1);
    julia> CpelRrbs.comp_mat_perm_stats(mml_diffs,js)
    4-element Array{Float64,1}:
     -0.7782325400000001
     0.0               
     0.0               
     0.7782325400000001
    ```
"""
function comp_mat_perm_stats(diffs::Vector{Float64},js::Vector{Int64})::Vector{Float64}
    
    # Return all possible signed sums
    return [comp_mat_j_stat(diffs,j) for j in js]

end
"""
    `comp_mat_stat_cmd(n,θ1s,θ2s)`

    Function that computes GJSD between pairs (matched case).

    # Examples
    ```julia-repl
    julia> n=[5]; θ1s=[[-1.0,1.0],[-1.0,1.0]]; θ2s=[[1.0,1.0],[1.0,1.0]];
    julia> ∇logZ1s = [CpelRrbs.get_∇logZ(n,θ1) for θ1 in θ1s];
    julia> ∇logZ2s = [CpelRrbs.get_∇logZ(n,θ2) for θ2 in θ2s];
    julia> CpelRrbs.comp_mat_stat_cmd(n,θ1s,θ2s,∇logZ1s,∇logZ2s)
     0.7888652058295635
    ```
"""
function comp_mat_stat_cmd(n::Vector{Int64},θ1s::Vector{Vector{Float64}},θ2s::Vector{Vector{Float64}},
    ∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}})::Float64
    
    # Compute stat
    cmds = 0.0
    @inbounds for s=1:length(θ1s)
        cmds += comp_cmd(n,θ1s[s],θ2s[s],∇logZ1s[s],∇logZ2s[s])
    end

    # Return mean GJSD
    return cmds/length(θ1s)

end
"""
    `mat_tests(n,θ1s,θ2s)`

    Function that performs hypothesis testing in matched samples group comparison. Note that the matched test
    is only performed for Tmml and Tnme.

    # Examples
    ```julia-repl
    julia> using Random; Random.seed!(1234);
    julia> n=[10]; θ1s=fill([0.5,0.5],5); θ2s=fill([-0.5,0.5],5);
    julia> tmml_test,tnme_test,tcmd_test = CpelRrbs.mat_tests(n,θ1s,θ2s)
    ((0.7782325400000001, 0.0625), (0.0, 1.0), (0.3253066763538722, 0.03125))
    julia> n=[10]; θ1s=fill([0.5,0.5],8); θ2s=fill([-0.5,0.5],8);
    julia> tmml_test,tnme_test,tcmd_test = CpelRrbs.mat_tests(n,θ1s,θ2s)
    ((0.7782325400000002, 0.00390625), (0.0, 1.0), (0.3253066763538723, 1.0))
    ```
"""
function mat_tests(n::Vector{Int64},θ1s::Vector{Vector{Float64}},θ2s::Vector{Vector{Float64}};Lmax::Int64=1000)::NTuple{3,NTuple{2,Float64}}

    # Compute gradients
    ∇logZ1s = get_all_∇logZs(n,θ1s)
    ∇logZ2s = get_all_∇logZs(n,θ2s)

    # Compute observed stats
    mml_diffs = comp_mat_diff_mml(n,∇logZ1s,∇logZ2s)
    nme_diffs = comp_mat_diff_nme(n,θ1s,θ2s,∇logZ1s,∇logZ2s)
    tcmd_obs = comp_mat_stat_cmd(n,θ1s,θ2s,∇logZ1s,∇logZ2s)

    # Get group label combinations to use
    exact = 2^length(θ1s)<=Lmax
    js = exact ? collect(0:2^length(θ1s)-1) : rand(0:2^length(θ1s)-1,Lmax)

    # Compute permutation stats
    tmml_perms = comp_mat_perm_stats(mml_diffs,js)
    tnme_perms = comp_mat_perm_stats(nme_diffs,js)

    # Compute observed stats
    tmml_obs = sum(mml_diffs)/length(θ1s)
    tnme_obs = sum(nme_diffs)/length(θ1s)
    
    # Compute p-values
    if exact
        # Compute exact p-value
        tmml_pval = sum(abs.(tmml_perms).>=abs(tmml_obs))/length(tmml_perms)
        tnme_pval = sum(abs.(tnme_perms).>=abs(tnme_obs))/length(tnme_perms)
    else
        # Compute p-value in MC setting
        tmml_pval = (1.0+sum(abs.(tmml_perms).>=abs(tmml_obs)))/(1.0+length(tmml_perms))
        tnme_pval = (1.0+sum(abs.(tnme_perms).>=abs(tnme_obs)))/(1.0+length(tnme_perms))
    end

    # Return (stat,pval) pairs
    return (tmml_obs,tmml_pval),(tnme_obs,tnme_pval),(tcmd_obs,NaN)
    
end