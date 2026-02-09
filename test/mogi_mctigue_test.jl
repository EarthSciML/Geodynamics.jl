@testsnippet MogiMcTigueSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: mtkcompile
    using DynamicQuantities
    using Geodynamics

    # Table 1 parameters from Taylor et al. (2021)
    const G_val = 4.0e9    # Pa (4.0 GPa)
    const ν_val = 0.25     # dimensionless
    const ΔP_val = 2.0e6   # Pa (2.0 MPa)

    # Helper to create a NonlinearProblem using the non-deprecated API
    function make_prob(compiled, params)
        NonlinearProblem(compiled, Dict(params))
    end
end

@testitem "Mogi Model Structure" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    sys=MogiModel()

    # Should have 3 equations (R, Ur, Uz)
    @test length(equations(sys)) == 3

    # Should have 3 unknowns
    @test length(unknowns(sys)) == 3

    # Check unknown names
    unames=Set(string.(Symbol.(unknowns(sys))))
    @test "R(t)" in unames
    @test "Ur(t)" in unames
    @test "Uz(t)" in unames

    # Should compile to zero dynamic unknowns (all algebraic)
    compiled=mtkcompile(sys)
    @test length(unknowns(compiled)) == 0
end

@testitem "McTigue Model Structure" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    sys=McTigueModel()

    # Should have 3 equations (R, Ur, Uz)
    @test length(equations(sys)) == 3

    # Should have 3 unknowns
    @test length(unknowns(sys)) == 3

    # Check unknown names
    unames=Set(string.(Symbol.(unknowns(sys))))
    @test "R(t)" in unames
    @test "Ur(t)" in unames
    @test "Uz(t)" in unames

    # Should compile to zero dynamic unknowns (all algebraic)
    compiled=mtkcompile(sys)
    @test length(unknowns(compiled)) == 0
end

@testitem "Mogi Unit Consistency" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    sys=MogiModel()

    # Check that all unknowns have units of meters
    for v in unknowns(sys)
        @test dimension(ModelingToolkit.get_unit(v)) == dimension(u"m")
    end
end

@testitem "McTigue Unit Consistency" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    sys=McTigueModel()

    # Check that all unknowns have units of meters
    for v in unknowns(sys)
        @test dimension(ModelingToolkit.get_unit(v)) == dimension(u"m")
    end
end

@testitem "Mogi Equation Correctness" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    # Test with Table 1 parameters: G=4 GPa, ν=0.25, ΔP=2 MPa
    # Source: a=1 km, d=4 km; evaluate at x=5 km

    a_val=1000.0   # m
    d_val=4000.0   # m
    x_val=5000.0   # m

    # Expected values (hand-computed from Eqs. 1-2)
    R_expected=sqrt(x_val^2+d_val^2)
    Ur_expected=a_val^3*ΔP_val*(1-ν_val)*x_val/(G_val*R_expected^3)
    Uz_expected=a_val^3*ΔP_val*(1-ν_val)*d_val/(G_val*R_expected^3)

    # Evaluate the compiled model
    sys=MogiModel()
    compiled=mtkcompile(sys)

    prob=make_prob(compiled,
        [
            compiled.G=>G_val,
            compiled.ν=>ν_val,
            compiled.ΔP=>ΔP_val,
            compiled.a=>a_val,
            compiled.d=>d_val,
            compiled.x=>x_val
        ])
    sol=solve(prob)

    @test sol[compiled.R] ≈ R_expected rtol=1e-10
    @test sol[compiled.Ur] ≈ Ur_expected rtol=1e-10
    @test sol[compiled.Uz] ≈ Uz_expected rtol=1e-10
end

@testitem "McTigue Equation Correctness" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    # Test with Table 1 parameters: G=4 GPa, ν=0.25, ΔP=2 MPa
    # Source: a=1 km, d=4 km; evaluate at x=5 km

    a_val=1000.0   # m
    d_val=4000.0   # m
    x_val=5000.0   # m

    # Expected values (hand-computed from Eqs. 3-4)
    R_expected=sqrt(x_val^2+d_val^2)
    ε=a_val/d_val
    correction=1+ε^3*(
        (1+ν_val)/(2*(-7+5ν_val))+
    15*d_val^2*(-2+ν_val)/(4*R_expected^2*(-7+5ν_val))
    )
    Ur_expected=a_val^3*ΔP_val*(1-ν_val)*x_val/(G_val*R_expected^3)*correction
    Uz_expected=a_val^3*ΔP_val*(1-ν_val)*d_val/(G_val*R_expected^3)*correction

    # Evaluate the compiled model
    sys=McTigueModel()
    compiled=mtkcompile(sys)

    prob=make_prob(compiled,
        [
            compiled.G=>G_val,
            compiled.ν=>ν_val,
            compiled.ΔP=>ΔP_val,
            compiled.a=>a_val,
            compiled.d=>d_val,
            compiled.x=>x_val
        ])
    sol=solve(prob)

    @test sol[compiled.R] ≈ R_expected rtol=1e-10
    @test sol[compiled.Ur] ≈ Ur_expected rtol=1e-10
    @test sol[compiled.Uz] ≈ Uz_expected rtol=1e-10
end

@testitem "McTigue Reduces to Mogi for Small ε" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    # When ε = a/d → 0, the McTigue correction → 1 and results should match Mogi
    # Use a very small source radius relative to depth

    a_val=10.0      # m (very small)
    d_val=10000.0   # m (10 km) → ε = 0.001
    x_val=5000.0    # m

    # Compile both models
    mogi_sys=MogiModel()
    mogi_compiled=mtkcompile(mogi_sys)
    mctigue_sys=McTigueModel()
    mctigue_compiled=mtkcompile(mctigue_sys)

    params_mogi=[
        mogi_compiled.G=>G_val, mogi_compiled.ν=>ν_val,
        mogi_compiled.ΔP=>ΔP_val, mogi_compiled.a=>a_val,
        mogi_compiled.d=>d_val, mogi_compiled.x=>x_val
    ]
    params_mctigue=[
        mctigue_compiled.G=>G_val, mctigue_compiled.ν=>ν_val,
        mctigue_compiled.ΔP=>ΔP_val, mctigue_compiled.a=>a_val,
        mctigue_compiled.d=>d_val, mctigue_compiled.x=>x_val
    ]

    sol_mogi=solve(make_prob(mogi_compiled, params_mogi))
    sol_mctigue=solve(make_prob(mctigue_compiled, params_mctigue))

    # For ε = 0.001, McTigue correction should be negligible
    @test sol_mctigue[mctigue_compiled.Ur] ≈ sol_mogi[mogi_compiled.Ur] rtol=1e-6
    @test sol_mctigue[mctigue_compiled.Uz] ≈ sol_mogi[mogi_compiled.Uz] rtol=1e-6
end

@testitem "Displacement Symmetry and Signs" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    # For a positive overpressure (inflation):
    # - Vertical displacement Uz should be positive (uplift)
    # - Horizontal displacement Ur should be positive (radially outward) for positive x

    a_val=1000.0
    d_val=4000.0
    x_val=5000.0

    sys=MogiModel()
    compiled=mtkcompile(sys)

    params=[
        compiled.G=>G_val, compiled.ν=>ν_val,
        compiled.ΔP=>ΔP_val, compiled.a=>a_val,
        compiled.d=>d_val, compiled.x=>x_val
    ]
    sol=solve(make_prob(compiled, params))

    # Positive overpressure → positive (outward) radial displacement
    @test sol[compiled.Ur] > 0

    # Positive overpressure → positive (upward) vertical displacement
    @test sol[compiled.Uz] > 0

    # Uz at x=0 (directly above source) should be maximum vertical displacement
    params_x0=[
        compiled.G=>G_val, compiled.ν=>ν_val,
        compiled.ΔP=>ΔP_val, compiled.a=>a_val,
        compiled.d=>d_val, compiled.x=>0.001 # near zero (avoid x=0 singularity)
    ]
    sol_x0=solve(make_prob(compiled, params_x0))
    @test sol_x0[compiled.Uz] > sol[compiled.Uz]
end

@testitem "Displacement Scaling Properties" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    # Mogi displacements scale linearly with ΔP and as a^3
    a_val=1000.0
    d_val=4000.0
    x_val=5000.0

    sys=MogiModel()
    compiled=mtkcompile(sys)

    # Base case
    params1=[
        compiled.G=>G_val, compiled.ν=>ν_val,
        compiled.ΔP=>ΔP_val, compiled.a=>a_val,
        compiled.d=>d_val, compiled.x=>x_val
    ]
    sol1=solve(make_prob(compiled, params1))

    # Double overpressure → double displacement
    params2=[
        compiled.G=>G_val, compiled.ν=>ν_val,
        compiled.ΔP=>2*ΔP_val, compiled.a=>a_val,
        compiled.d=>d_val, compiled.x=>x_val
    ]
    sol2=solve(make_prob(compiled, params2))
    @test sol2[compiled.Ur] ≈ 2 * sol1[compiled.Ur] rtol=1e-10
    @test sol2[compiled.Uz] ≈ 2 * sol1[compiled.Uz] rtol=1e-10

    # Double radius → 8× displacement (a^3 scaling)
    params3=[
        compiled.G=>G_val, compiled.ν=>ν_val,
        compiled.ΔP=>ΔP_val, compiled.a=>2*a_val,
        compiled.d=>d_val, compiled.x=>x_val
    ]
    sol3=solve(make_prob(compiled, params3))
    @test sol3[compiled.Ur] ≈ 8 * sol1[compiled.Ur] rtol=1e-10
    @test sol3[compiled.Uz] ≈ 8 * sol1[compiled.Uz] rtol=1e-10
end

@testitem "Ur/Uz Ratio Equals x/d" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    # For the Mogi model, Ur/Uz = x/d (from Eqs. 1-2)
    a_val=1000.0
    d_val=4000.0

    sys=MogiModel()
    compiled=mtkcompile(sys)

    for x_val in [1000.0, 3000.0, 5000.0, 10000.0, 20000.0]
        params=[
            compiled.G=>G_val, compiled.ν=>ν_val,
            compiled.ΔP=>ΔP_val, compiled.a=>a_val,
            compiled.d=>d_val, compiled.x=>x_val
        ]
        sol=solve(make_prob(compiled, params))
        @test sol[compiled.Ur] / sol[compiled.Uz] ≈ x_val / d_val rtol=1e-10
    end
end

@testitem "Far-Field Decay" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    # Displacements should decrease with increasing distance
    # For large x, Ur ~ x/R^3 ~ 1/x^2 and Uz ~ d/R^3 ~ 1/x^3 approximately

    a_val=1000.0
    d_val=4000.0

    sys=MogiModel()
    compiled=mtkcompile(sys)

    distances=[5000.0, 10000.0, 15000.0, 20000.0, 25000.0]
    Ur_vals=Float64[]
    Uz_vals=Float64[]

    for x_val in distances
        params=[
            compiled.G=>G_val, compiled.ν=>ν_val,
            compiled.ΔP=>ΔP_val, compiled.a=>a_val,
            compiled.d=>d_val, compiled.x=>x_val
        ]
        sol=solve(make_prob(compiled, params))
        push!(Ur_vals, sol[compiled.Ur])
        push!(Uz_vals, sol[compiled.Uz])
    end

    # Both should be monotonically decreasing
    for i in 1:(length(distances) - 1)
        @test Ur_vals[i] > Ur_vals[i + 1]
        @test Uz_vals[i] > Uz_vals[i + 1]
    end
end

@testitem "McTigue Always Greater Than Mogi for ε > 0" setup=[MogiMcTigueSetup] tags=[:mogi_mctigue] begin
    # For ν = 0.25, the McTigue correction is > 1, so McTigue displacements
    # should be slightly larger than Mogi displacements

    d_val=4000.0
    x_val=5000.0

    mogi_sys=MogiModel()
    mogi_compiled=mtkcompile(mogi_sys)
    mctigue_sys=McTigueModel()
    mctigue_compiled=mtkcompile(mctigue_sys)

    for a_val in [500.0, 1000.0, 1500.0, 2000.0]
        params_mogi=[
            mogi_compiled.G=>G_val, mogi_compiled.ν=>ν_val,
            mogi_compiled.ΔP=>ΔP_val, mogi_compiled.a=>a_val,
            mogi_compiled.d=>d_val, mogi_compiled.x=>x_val
        ]
        params_mctigue=[
            mctigue_compiled.G=>G_val, mctigue_compiled.ν=>ν_val,
            mctigue_compiled.ΔP=>ΔP_val, mctigue_compiled.a=>a_val,
            mctigue_compiled.d=>d_val, mctigue_compiled.x=>x_val
        ]

        sol_mogi=solve(make_prob(mogi_compiled, params_mogi))
        sol_mctigue=solve(make_prob(mctigue_compiled, params_mctigue))

        # McTigue correction is positive for ν=0.25, so magnitudes should be larger
        @test abs(sol_mctigue[mctigue_compiled.Ur]) > abs(sol_mogi[mogi_compiled.Ur])
        @test abs(sol_mctigue[mctigue_compiled.Uz]) > abs(sol_mogi[mogi_compiled.Uz])
    end
end
