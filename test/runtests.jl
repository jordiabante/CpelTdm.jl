using CpelTdm, Test

@testset "Bioinformatics" begin
    # Aligning strand
    @test CpelTdm.get_align_strand(true,UInt16(99),UInt16(147)) == "OT"
    @test CpelTdm.get_align_strand(true,UInt16(83),UInt16(163)) == "OB"
    @test CpelTdm.get_align_strand(true,UInt16(147),UInt16(99)) == "CTOT"
    @test CpelTdm.get_align_strand(true,UInt16(163),UInt16(83)) == "CTOB"
    # Mean coverage
    xobs=[[1,-1] for i=1:10]; append!(xobs,[[1,0] for i=1:10]);
    @test CpelTdm.mean_cov(xobs)==15.0
    @test CpelTdm.get_obs_per_cpg(xobs)==[20,10]
end

@testset "Inference" begin
    # Energy function
    Ux_fun = CpelTdm.create_Ux([2,2],[1.0,1.0],1.0);
    @test Ux_fun([1,1,1,1]) == -7.0
    # Matrices
    @test isapprox(CpelTdm.get_W(4,0.0,0.0),[4.0 4.0; 4.0 4.0];atol=2)
    @test isapprox(CpelTdm.get_V([0.0,0.0],0.0),[1.0 1.0; 1.0 1.0];atol=2)
    @test CpelTdm.get_u(0.0) == [1.0,1.0]
    @test CpelTdm.comp_Z([4],[0.0],0.0) ≈ 16
    @test CpelTdm.comp_g([1],[0.0],0.0,0.0,0.0) == 2.0
    @test isapprox(CpelTdm.comp_lkhd([1,1,0,1,1],[5],[1.0],1.0),0.9536;atol=2)
    # MML
    n=[1,1,1]; θ=[1.0,-1.0,1.0,0.0]; ∇logZ = CpelTdm.get_∇logZ(n,θ);
    @test isapprox(CpelTdm.comp_mml(n,∇logZ),0.62693236)
    # NME
    n=[5]; α=[0.0]; β=0.0; θ=vcat([α,β]...); ∇logZ=CpelTdm.get_∇logZ(n,θ);
    @test isapprox(CpelTdm.comp_nme(n,θ,∇logZ),1.0)
    # CMD
    n=[20]; θ1=[-5.0,5.0]; θ2=[5.0,5.0];
    ∇logZ1s = CpelTdm.get_∇logZ(n,θ1); ∇logZ2s = CpelTdm.get_∇logZ(n,θ2);
    @test isapprox(CpelTdm.comp_cmd(n,θ1,θ2,∇logZ1s,∇logZ2s),1.0;atol=2)
end

@testset "Hypothesis Testing" begin
    # Unmatched test
    n=[10]; θ1s=fill([0.0,0.0],5); θ2s=fill([0.0,2.5],5);
    tmml_test,tnme_test,tpdm_test = CpelTdm.unmat_tests(n,θ1s,θ2s)
    @test tmml_test == (0.0, 1.0)
    @test tnme_test == (0.84782978, 0.007936507936507936)
    @test tpdm_test == (0.43786215630140846, 0.007936507936507936)
    # Matched test
    n=[10]; θ1s=fill([0.0,0.0],5); θ2s=fill([0.0,2.5],5);
    tmml_test,tnme_test,tpdm_test = CpelTdm.mat_tests(n,θ1s,θ2s)
    @test tmml_test == (0.0, 1.0)
    @test tnme_test == (0.84782978, 0.0625)
    @test tpdm_test[1] == 0.4378621563014087
end