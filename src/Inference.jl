###################################################################################################
# FUNCTIONS
###################################################################################################
"""
    `create_Ux([N1,...,NK],[α1,...,αK],β)`

    Function that creates a function to compute potential energy for a region with [N1,...,NK],
    parameters [α1,...,αK], and correlation β.

    # Examples
    ```julia-repl
    julia> Ux_fun = CpelTdm.create_Ux([2,2],[1.0,1.0],1.0);
    julia> Ux_fun([1,1,1,1])
    -7.0
    ```
"""
function create_Ux(n::Vector{Int64},a::Vector{Float64},b::Float64)

    # Define function given n, α, and β
    function Ux(x::Vector{Int64})::Float64
        u = a[1] * sum(x[1:n[1]]) + b * sum(x[1:(end-1)] .* x[2:end])
        @inbounds for i in 2:length(n)
            u += a[i] * sum(x[(sum(n[1:(i-1)])+1):(sum(n[1:i]))])
        end
        return -u
    end

    # Return function
    return Ux

end # end create_Ux
"""
    `get_W(N,α,β)`

    Function that returns W^{N-1} computed efficiently via sum of rank-1 matrices
    assuming α and β.

    # Examples
    ```julia-repl
    julia> CpelTdm.get_W(4,0.0,0.0)
    2×2 Array{Float64,2}:
    4.0  4.0
    4.0  4.0
    ```
"""
function get_W(n::Int64,a::Float64,b::Float64)::Array{Float64,2}

    # if n=1 then return identity
    n==1 && return [1.0 0.0;0.0 1.0]

    # Compute relevant quantities
    exp_a = exp(a)
    exp_b = exp(b)
    cosh_a = 0.5*(exp_a+1.0/exp_a)
    sinh_a = exp_a-cosh_a

    # Eigenvalues
    aux1 = exp_b * cosh_a
    aux2 = sqrt(1.0 + exp_b^4*sinh_a^2)
    lambda1N = (aux1-1.0/exp_b*aux2)^(n-1)
    lambda2N = (aux1+1.0/exp_b*aux2)^(n-1)

    # Eigenvectors
    aux1 = -exp_b^2 * sinh_a
    e1 = [aux1-aux2; 1.0]
    e1 /= sqrt(e1'*e1)
    e2 = [aux1+aux2; 1.0]
    e2 /= sqrt(e2'*e2)

    # Return W^{N-1}
    return e1*e1'*lambda1N + e2*e2'*lambda2N

end
"""
    `get_V([α1,α2],β)`

    Function that returns V matrix assuming [α1,α2] and β.

    # Examples
    ```julia-repl
    julia> CpelTdm.get_V([1.0,1.0],1.0)
    2×2 Array{Float64,2}:
    1.0       0.367879
    0.367879  7.38906
    ```
"""
function get_V(a::Vector{Float64},b::Float64)::Array{Float64,2}

    # Compute auxiliary vars
    exp_b = exp(b)
    exp_a_p = exp(0.5*sum(a))
    exp_a_m = exp(0.5*(a[2]-a[1]))

    # Return V
    return [exp_b/exp_a_p exp_a_m/exp_b;1.0/(exp_b*exp_a_m) exp_b*exp_a_p]

end
"""
    `get_u(α)`

    Function that boundary vector for Z computation.

    # Examples
    ```julia-repl
    julia> CpelTdm.get_u(0.0)
    2-element Vector{Float64}:
     1.0
     1.0
    ```
"""
function get_u(a::Float64)::Vector{Float64}

    # Compute auxiliary
    exp_aux = exp(a/2.0)

    # Return u
    return [1.0/exp_aux; exp_aux]

end
"""
    `comp_Z([N1,...,NK],[α1,...,αK],β)`

    Compute partition function of a model with [N1,...,NK] CpG cites and with
    parameters [α1,...,αK] and β.

    # Examples
    ```julia-repl
    julia> CpelTdm.comp_Z([1,1,1],[1.0,1.0,1.0],1.0)
    155.37102759254836
    ```
"""
function comp_Z(n::Vector{Int64},a::Vector{Float64},b::Float64)::Float64

    # Boundary conditions.
    y = get_u(a[1])'*get_W(n[1],a[1],b)
    if length(n)>1
        y *= prod([get_V(a[(i-1):i],b)*get_W(n[i],a[i],b) for i in 2:length(n)])
    end
    y *= get_u(a[end])

    # Return Z
    return max(1e-100,y)

end
"""
    `get_∇logZ([N1,...,NK],θhat)`

    Numerically computes the gradient of the log partition function of a model with [N1,...,NK] 
    CpG cites and estimated parameter vector θhat. This is equivalent to computing the expected 
    value of sufficient statistis (SS).

    # Examples
    ```julia-repl
    julia> CpelTdm.get_∇logZ([1,1,1],[1.0,-1.0,1.0,0.0])
    4-element Array{Float64,1}:
     0.7615941559702503
     -0.7615941559702503
     0.7615941559335818
     -1.1600513167692226
    ```
"""
function get_∇logZ(n::Vector{Int64},θhat::Vector{Float64})::Vector{Float64}

    # Define function
    function f(θ::Vector{Float64})
        log(comp_Z(n,θ[1:(end-1)],θ[end]))
    end

    # Return ∇logZ(θ)
    return Calculus.gradient(f,θhat)

end
"""
    `check_boundary(θhat)`

    Function that returns a bool indicating whether model with parameter estimate vector θhat is on the
    boundary of the parameter space in any of its components.

    # Examples
    ```julia-repl
    julia> CpelTdm.check_boundary([1.0,1.0,1.0,1.0])
    false
    julia> CpelTdm.check_boundary([1.0,5.0,1.0,1.0])
    true
    ```
"""
function check_boundary(θhat::Vector{Float64})::Bool

    # Return true θhat on boundary.
    return any(isapprox.(abs.(θhat),ETA_MAX_ABS;atol=5e-2)) || any(abs.(θhat).>ETA_MAX_ABS)

end
"""
    `comp_g([R1,...,RK],[α1,...,αK],β,αp1,αp2)`

    Compute scaling factor in a model with [R1,...,RK] unobserved CpG cites from each block, with
    parameters [α1,...,αK] and β. The transfer-matrix method is used to obtain a computationally
    efficient expression without the need of recursive expressions. αp1 and αp2 are determined by the
    respective kind of boundary.

    # Examples
    ```julia-repl
    julia> CpelTdm.comp_g([1],[1.0],1.0,3.0,3.0)
    20.135323991555527
    ```
"""
function comp_g(r::Vector{Int64},a::Vector{Float64},b::Float64,ap1::Float64,ap2::Float64)::Float64

    # Return one if r=0
    r==0 && return 1.0

    # Boundary conditions
    y = get_u(ap1)'*get_W(r[1],a[1],b)
    if length(r)>1
        y *= prod([get_V(a[(i-1):i],b)*get_W(r[i],a[i],b) for i in 2:length(r)])
    end
    y *= get_u(ap2)

    # Return scaling factor
    return max(1e-100,y)

end
"""
    `comp_lkhd(X,[N1,...,NK],[α1,...,αK],β)`

    Compute likelihood of a partial/full observation X that can be missing values anywhere (given by
    0's).

    # Examples
    ```julia-repl
    julia> CpelTdm.comp_lkhd([1,1,0,1,1],[5],[1.0],1.0)
    0.953646032691218
    ```
"""
function comp_lkhd(x::Vector{Int64},n::Vector{Int64},a::Vector{Float64},b::Float64)::Float64

    # Avoid the rest if N=1
    length(x)==1 && return 0.5 * exp(x[1]*a[1]) / cosh(a[1])

    # Get partition function and energy function
    Z = comp_Z(n,a,b)
    Ux = create_Ux(n,a,b)

    # Find changes to/from 0 in x vector
    zerost = findall(isequal(-1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1
    zeroend = findall(isequal(1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1

    # Determine whether it starts/finishes with 0
    x[1]==0 && pushfirst!(zerost,1)
    x[end]==0 && push!(zeroend,length(x)+1)

    # Find subregion label for each CpG site
    subid = vcat([i*ones(Int64,n[i]) for i in 1:length(n)]...)

    # Get overall scaling factor as product of individual factors
    sf = 1.0
    @inbounds for i in 1:length(zerost)
        # Find boundary conditions (other X's or boundary of X)
        ap1 = zerost[i]==1 ? a[1] : 2.0*x[zerost[i]-1]*b+a[subid[zerost[i]]]
        ap2 = zeroend[i]==sum(n)+1 ? a[end] : 2.0*x[zeroend[i]]*b+a[subid[zeroend[i]-1]]

        # Figure out b (block IDs of missing) and r (α indices)
        b_id = subid[zerost[i]:(zeroend[i]-1)]
        n_miss = [count(x->x==n,b_id) for n in unique(b_id)]

        # Call scaling factor function
        sf *= comp_g(n_miss,a[unique(b_id)],b,ap1,ap2)

    end

    # Return energy function properly scaled.
    return exp(-Ux(x)) * sf / Z

end
"""
    `create_Llkhd([N1,...,NK],XOBS)`

    Create function to compute the minus log-likelihood function for a region with N CpG sites given
    the M partial observations XOBS.

    # Examples
    ```julia-repl
    julia> using Random; Random.seed!(1234); n=[4];
    julia> xobs=[CpelTdm.gen_x_mc(n,[0.0,0.0]) for _ in 1:20];
    julia> LogLike=CpelTdm.create_Llkhd(n,xobs);
    julia> LogLike([0.0,0.0,0.0])
    55.45177444479562
    ```
"""
function create_Llkhd(n::Vector{Int64},xobs::Array{Vector{Int64},1})

    # Define minus log-likelihood function
    function Llkhd_fun(theta::Vector{Float64})::Float64
        # Get parameters
        aux = 0.0
        a = theta[1:(end-1)]
        b = theta[end]

        # Get energy function and partition function
        Ux = create_Ux(n,a,b)
        logZ = log(comp_Z(n,a,b))

        # Initialize variables used in for loop
        ap1 = ap2 = 0.0                 # Boundary values left and right
        b_id = n_miss = []              # Block IDs of missing CpG; number missing per ID
        zerost = zeroend = [] # Indices of zeros
        subid = vcat([i*ones(Int64,n[i]) for i in 1:length(n)]...)    # ID's of CpGs

        # Contribution of each observation
        for x in xobs

            # Find changes to/from 0 in x vector
            zerost = findall(isequal(-1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1
            zeroend = findall(isequal(1),abs.(x[2:end]) - abs.(x[1:(end-1)])) .+ 1

            # Determine whether it starts/finishes with 0
            x[1]==0 && pushfirst!(zerost,1)
            x[end]==0 && push!(zeroend,length(x)+1)

            # Get scaling factors due to all missing data
            @inbounds for i in 1:length(zerost)
                # Find boundary conditions (other X's or boundary of X)
                ap1 = zerost[i]==1 ? a[1] : 2.0*x[zerost[i]-1]*b+a[subid[zerost[i]]]
                ap2 = zeroend[i]==sum(n)+1 ? a[end] : 2.0*x[zeroend[i]]*b+a[subid[zeroend[i]-1]]

                # Figure out b (block IDs of missing) and r (α indices)
                b_id = subid[zerost[i]:(zeroend[i]-1)]
                n_miss = [count(x->x==n,b_id) for n in unique(b_id)]

                # Call scaling factor function
                aux += log(comp_g(n_miss,a[unique(b_id)],b,ap1,ap2))
            end

            # Add log of it to Llkhd
            aux += -Ux(x)
        end

        # Return MINUS log-likelihood. Add as many logZ as samples we have
        -aux + length(xobs) * logZ
    end

    # Return function
    return Llkhd_fun

end
"""
    `est_alpha(XOBS)`

    Estimate parameter α in N=1 case (first entry in returned vector).

    # Examples
    ```julia-repl
    julia> using Random; Random.seed!(1234); n=[1];
    julia> xobs=[CpelTdm.gen_x_mc(n,[0.0,0.0]) for _ in 1:20];
    julia> CpelTdm.est_alpha(xobs)
    2-element Array{Float64,1}:
     0.10033534773107566
     0.0
    ```
"""
function est_alpha(xobs::Array{Vector{Int64},1})::Vector{Float64}
    # Derivation of estimate:
        # log(p)+log(2) = α-log(cosh(α)); log(1-p)+log(2) = -α-log(cosh(α))
        # α = (log(p)-log(1-p))/2

    # Proportion of X=1
    phat = length(findall(x->x==[1],xobs)) / length(xobs)
    a = 0.5*(log(phat)-log(1.0-phat))

    # Return estimate
    return [min(max(-ETA_MAX_ABS,a),ETA_MAX_ABS),0.0]

end
"""
    `est_theta_sa([N1,...,NK],XOBS)`

    Estimate parameter vector [α1,...,αK,β] using simulated annealing.

    # Examples
    ```julia-repl
    julia> using Random; Random.seed!(1234); n=[4];
    julia> xobs=[CpelTdm.gen_x_mc(n,[0.0,0.0]) for _ in 1:100];
    julia> CpelTdm.est_theta_sa(n,xobs)
    2-element Array{Float64,1}:
     -0.08515675262874925 
     -0.040190205815463766
    ```
"""
function est_theta_sa(n::Vector{Int64},xobs::Array{Vector{Int64},1})::Vector{Float64}

    # If N=1, then estimate α
    sum(n)==1 && return est_alpha(xobs)

    # Define boundaries and initialization
    L = create_Llkhd(n,xobs)
    init = zeros(Float64,length(n)+1)
    lower = -ETA_MAX_ABS * ones(Float64,length(n)+1)
    upper = ETA_MAX_ABS * ones(Float64,length(n)+1)
    opts = Optim.Options(iterations=10^4,show_trace=false,store_trace=false)

    # Boxed Simulated Annealing (SAMI)
    optim = Optim.optimize(L,lower,upper,init,SAMIN(rt=1e-4;f_tol=1e-3,verbosity=0),opts)

    # Return estimate
    return optim.minimizer

end
###################################################################################################
# OUTPUT QUANTITIES
###################################################################################################
"""
    `comp_mml([N1,...,NK],∇logZ)`

    Function that computes mean methylation level (MML) over the entire region of interest.

    # Examples
    ```julia-repl
    julia> n=[1,1,1]; θ=[1.0,-1.0,1.0,0.0];
    julia> ∇logZ = CpelTdm.get_∇logZ(n,θ);
    julia> CpelTdm.comp_mml(n,∇logZ)
    0.62693236
    ```
"""
function comp_mml(n::Vector{Int64},∇logZ::Vector{Float64})::Float64

    # Return
    return abs(round(0.5*(1.0+1.0/sum(n)*sum(∇logZ[1:length(n)]));digits=8))

end
"""
    `comp_nme([N1,...,NK],θ,∇logZ)`

    Function that computes normalized methylation entropy (NME) using the ∇logZ information.

    # Examples
    ```julia-repl
    julia> n=[5]; α=[0.0]; β=0.0; θ=vcat([α,β]...);
    julia> ∇logZ=CpelTdm.get_∇logZ(n,θ);
    julia> CpelTdm.comp_nme(n,θ,∇logZ)
    1.0
    ```
"""
function comp_nme(n::Vector{Int64},θ::Vector{Float64},∇logZ::Vector{Float64})::Float64

    # Return
    return abs(round(1.0/(sum(n)*LOG2)*(log(comp_Z(n,θ[1:(end-1)],θ[end]))-θ'*∇logZ);digits=8))

end
"""
    `comp_cmd([N1,...,NK],θ1,θ2,∇logZ1s,∇logZ2s)`

    Function that computes the normalized geometric Jensen-Shannon divergence between two Ising models 
    parametrized by θ1 and θ2, respectively.

    # Examples
    ```julia-repl
    julia> n=[20]; θ1=θ2=[1.0,1.0];
    julia> ∇logZ1s = CpelTdm.get_∇logZ(n,θ1);
    julia> ∇logZ2s = CpelTdm.get_∇logZ(n,θ2);
    julia> CpelTdm.comp_cmd(n,θ1,θ2,∇logZ1s,∇logZ2s)
    0.0
    julia> n=[20]; θ1=[-5.0,5.0]; θ2=[5.0,5.0];
    julia> ∇logZ1s = CpelTdm.get_∇logZ(n,θ1);
    julia> ∇logZ2s = CpelTdm.get_∇logZ(n,θ2);
    julia> CpelTdm.comp_cmd(n,θ1,θ2,∇logZ1s,∇logZ2s)
    1.0
    ```
"""
function comp_cmd(n::Vector{Int64},θ1::Vector{Float64},θ2::Vector{Float64},∇logZ1::Vector{Float64},∇logZ2::Vector{Float64})::Float64

    # Parameter geometric mixture of Ising
    θγ = 0.5.*(θ1+θ2)

    # Get partition functions
    logZ1 = log(comp_Z(n,θ1[1:(end-1)],θ1[end]))
    logZ2 = log(comp_Z(n,θ2[1:(end-1)],θ2[end]))
    logZγ = log(comp_Z(n,θγ[1:(end-1)],θγ[end]))

    # Compute numerator of 1-CMD
    cmd = logZ1 + logZ2 - (θ1'*∇logZ1 + θ2'*∇logZ2)

    # Normalize 1-CMD
    cmd /= 2.0*logZγ - θγ'* (∇logZ1 + ∇logZ2)
    
    # Return GJSD
    return 1.0 - cmd

end
"""
    `comp_cmd_unnorm([N1,...,NK],θ1,θ2,∇logZ1s,∇logZ2s)`

    Function that computes the unnormalized geometric Jensen-Shannon divergence between two Ising models 
    parametrized by θ1 and θ2, respectively.

    # Examples
    ```julia-repl
    julia> n=[20]; θ1=θ2=[1.0,1.0];
    julia> ∇logZ1s = CpelTdm.get_∇logZ(n,θ1);
    julia> ∇logZ2s = CpelTdm.get_∇logZ(n,θ2);
    julia> CpelTdm.comp_cmd_unnorm(n,θ1,θ2,∇logZ1s,∇logZ2s)
    0.0
    julia> n=[20]; θ1=[-5.0,5.0]; θ2=[5.0,5.0];
    julia> ∇logZ1s = CpelTdm.get_∇logZ(n,θ1);
    julia> ∇logZ2s = CpelTdm.get_∇logZ(n,θ2);
    julia> CpelTdm.comp_cmd_unnorm(n,θ1,θ2,∇logZ1s,∇logZ2s)
    1.3880197058418844
    ```
"""
function comp_cmd_unnorm(n::Vector{Int64},θ1::Vector{Float64},θ2::Vector{Float64},∇logZ1::Vector{Float64},∇logZ2::Vector{Float64})::Float64

    # Parameter geometric mixture of Ising
    θγ = 0.5.*(θ1+θ2)

    # Get partition functions
    logZ1 = log(comp_Z(n,θ1[1:(end-1)],θ1[end]))
    logZ2 = log(comp_Z(n,θ2[1:(end-1)],θ2[end]))
    logZγ = log(comp_Z(n,θγ[1:(end-1)],θγ[end]))

    # Return unnormalized GJSD
    return 2.0*logZγ - θγ'* (∇logZ1 + ∇logZ2) - (logZ1 + logZ2 - (θ1'*∇logZ1 + θ2'*∇logZ2))

end