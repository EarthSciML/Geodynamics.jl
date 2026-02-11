@testsnippet MogiMcTigueSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: mtkcompile
    using DynamicQuantities
    using NonlinearSolve
    using Geodynamics

    # Table 1 parameters from Taylor et al. (2021)
    G_val = 4.0e9    # Pa (4.0 GPa)
    ν_val = 0.25     # dimensionless
    ΔP_val = 2.0e6   # Pa (2.0 MPa)

    # Helper to create a NonlinearProblem and solve it
    function make_prob(compiled, params)
        NonlinearProblem(compiled, Dict(params))
    end
end

@testitem "Mogi Model Structure" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    sys = MogiModel()

    # Should have 3 equations (R, Ur, Uz)
    @test length(equations(sys)) == 3

    # Should have 3 unknowns
    @test length(unknowns(sys)) == 3

    # Check unknown names
    unames = Set(string.(Symbol.(unknowns(sys))))
    @test "R(t)" in unames
    @test "Ur(t)" in unames
    @test "Uz(t)" in unames

    # Should compile to zero dynamic unknowns (all algebraic)
    compiled = mtkcompile(sys)
    @test length(unknowns(compiled)) == 0
end

@testitem "McTigue Model Structure" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    sys = McTigueModel()

    # Should have 3 equations (R, Ur, Uz)
    @test length(equations(sys)) == 3

    # Should have 3 unknowns
    @test length(unknowns(sys)) == 3

    # Check unknown names
    unames = Set(string.(Symbol.(unknowns(sys))))
    @test "R(t)" in unames
    @test "Ur(t)" in unames
    @test "Uz(t)" in unames

    # Should compile to zero dynamic unknowns (all algebraic)
    compiled = mtkcompile(sys)
    @test length(unknowns(compiled)) == 0
end

@testitem "Mogi Unit Consistency" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    sys = MogiModel()

    # Check that all unknowns have units of meters
    for v in unknowns(sys)
        @test dimension(ModelingToolkit.get_unit(v)) == dimension(u"m")
    end
end

@testitem "McTigue Unit Consistency" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    sys = McTigueModel()

    # Check that all unknowns have units of meters
    for v in unknowns(sys)
        @test dimension(ModelingToolkit.get_unit(v)) == dimension(u"m")
    end
end

@testitem "Mogi Equation Correctness" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # Test with Table 1 parameters: G=4 GPa, ν=0.25, ΔP=2 MPa
    # Source: a=1 km, d=4 km; evaluate at x=5 km
    # Reference values independently computed from Eqs. 1-2:
    #   R = sqrt(5000^2 + 4000^2) = 6403.12423743...
    #   Ur = 1e9 * 2e6 * 0.75 * 5000 / (4e9 * R^3) = 7.13834...e-3 m
    #   Uz = 1e9 * 2e6 * 0.75 * 4000 / (4e9 * R^3) = 5.71068...e-3 m

    a_val = 1000.0   # m
    d_val = 4000.0   # m
    x_val = 5000.0   # m

    R_ref = 6403.1242374328485  # m, precomputed
    Ur_ref = 0.00714209276929601  # m, precomputed
    Uz_ref = 0.0057136742154368075  # m, precomputed

    # Evaluate the compiled model
    sys = MogiModel()
    compiled = mtkcompile(sys)

    prob = make_prob(
        compiled,
        [
            compiled.G => G_val,
            compiled.ν => ν_val,
            compiled.ΔP => ΔP_val,
            compiled.a => a_val,
            compiled.d => d_val,
            compiled.x => x_val,
        ]
    )
    sol = solve(prob)

    @test sol[compiled.R] ≈ R_ref rtol = 1.0e-10
    @test sol[compiled.Ur] ≈ Ur_ref rtol = 1.0e-10
    @test sol[compiled.Uz] ≈ Uz_ref rtol = 1.0e-10
end

@testitem "McTigue Equation Correctness" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # Test with Table 1 parameters: G=4 GPa, ν=0.25, ΔP=2 MPa
    # Source: a=1 km, d=4 km; evaluate at x=5 km
    # Reference values independently computed from Eqs. 3-4:
    #   R = 6403.124..., correction = 1.005261...,
    #   Ur = 0.007180..., Uz = 0.005744...

    a_val = 1000.0   # m
    d_val = 4000.0   # m
    x_val = 5000.0   # m

    R_ref = 6403.1242374328485  # m, precomputed
    Ur_ref = 0.0071796659144155  # m, precomputed
    Uz_ref = 0.005743732731532399  # m, precomputed

    # Evaluate the compiled model
    sys = McTigueModel()
    compiled = mtkcompile(sys)

    prob = make_prob(
        compiled,
        [
            compiled.G => G_val,
            compiled.ν => ν_val,
            compiled.ΔP => ΔP_val,
            compiled.a => a_val,
            compiled.d => d_val,
            compiled.x => x_val,
        ]
    )
    sol = solve(prob)

    @test sol[compiled.R] ≈ R_ref rtol = 1.0e-10
    @test sol[compiled.Ur] ≈ Ur_ref rtol = 1.0e-10
    @test sol[compiled.Uz] ≈ Uz_ref rtol = 1.0e-10
end

@testitem "McTigue Reduces to Mogi for Small ε" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # When ε = a/d → 0, the McTigue correction → 1 and results should match Mogi
    # Use a very small source radius relative to depth

    a_val = 10.0      # m (very small)
    d_val = 10000.0   # m (10 km) → ε = 0.001
    x_val = 5000.0    # m

    # Compile both models
    mogi_sys = MogiModel()
    mogi_compiled = mtkcompile(mogi_sys)
    mctigue_sys = McTigueModel()
    mctigue_compiled = mtkcompile(mctigue_sys)

    params_mogi = [
        mogi_compiled.G => G_val, mogi_compiled.ν => ν_val,
        mogi_compiled.ΔP => ΔP_val, mogi_compiled.a => a_val,
        mogi_compiled.d => d_val, mogi_compiled.x => x_val,
    ]
    params_mctigue = [
        mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
        mctigue_compiled.ΔP => ΔP_val, mctigue_compiled.a => a_val,
        mctigue_compiled.d => d_val, mctigue_compiled.x => x_val,
    ]

    sol_mogi = solve(make_prob(mogi_compiled, params_mogi))
    sol_mctigue = solve(make_prob(mctigue_compiled, params_mctigue))

    # For ε = 0.001, McTigue correction should be negligible
    @test sol_mctigue[mctigue_compiled.Ur] ≈ sol_mogi[mogi_compiled.Ur] rtol = 1.0e-6
    @test sol_mctigue[mctigue_compiled.Uz] ≈ sol_mogi[mogi_compiled.Uz] rtol = 1.0e-6
end

@testitem "Displacement Symmetry and Signs" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # For a positive overpressure (inflation):
    # - Vertical displacement Uz should be positive (uplift)
    # - Horizontal displacement Ur should be positive (radially outward) for positive x

    a_val = 1000.0
    d_val = 4000.0
    x_val = 5000.0

    sys = MogiModel()
    compiled = mtkcompile(sys)

    params = [
        compiled.G => G_val, compiled.ν => ν_val,
        compiled.ΔP => ΔP_val, compiled.a => a_val,
        compiled.d => d_val, compiled.x => x_val,
    ]
    sol = solve(make_prob(compiled, params))

    # Positive overpressure → positive (outward) radial displacement
    @test sol[compiled.Ur] > 0

    # Positive overpressure → positive (upward) vertical displacement
    @test sol[compiled.Uz] > 0

    # Uz at x=0 (directly above source) should be maximum vertical displacement
    params_x0 = [
        compiled.G => G_val, compiled.ν => ν_val,
        compiled.ΔP => ΔP_val, compiled.a => a_val,
        compiled.d => d_val, compiled.x => 0.001, # near zero (avoid x=0 singularity)
    ]
    sol_x0 = solve(make_prob(compiled, params_x0))
    @test sol_x0[compiled.Uz] > sol[compiled.Uz]
end

@testitem "Displacement Scaling Properties" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # Mogi displacements scale linearly with ΔP and as a^3
    a_val = 1000.0
    d_val = 4000.0
    x_val = 5000.0

    sys = MogiModel()
    compiled = mtkcompile(sys)

    # Base case
    params1 = [
        compiled.G => G_val, compiled.ν => ν_val,
        compiled.ΔP => ΔP_val, compiled.a => a_val,
        compiled.d => d_val, compiled.x => x_val,
    ]
    sol1 = solve(make_prob(compiled, params1))

    # Double overpressure → double displacement
    params2 = [
        compiled.G => G_val, compiled.ν => ν_val,
        compiled.ΔP => 2 * ΔP_val, compiled.a => a_val,
        compiled.d => d_val, compiled.x => x_val,
    ]
    sol2 = solve(make_prob(compiled, params2))
    @test sol2[compiled.Ur] ≈ 2 * sol1[compiled.Ur] rtol = 1.0e-10
    @test sol2[compiled.Uz] ≈ 2 * sol1[compiled.Uz] rtol = 1.0e-10

    # Double radius → 8× displacement (a^3 scaling)
    params3 = [
        compiled.G => G_val, compiled.ν => ν_val,
        compiled.ΔP => ΔP_val, compiled.a => 2 * a_val,
        compiled.d => d_val, compiled.x => x_val,
    ]
    sol3 = solve(make_prob(compiled, params3))
    @test sol3[compiled.Ur] ≈ 8 * sol1[compiled.Ur] rtol = 1.0e-10
    @test sol3[compiled.Uz] ≈ 8 * sol1[compiled.Uz] rtol = 1.0e-10
end

@testitem "Ur/Uz Ratio Equals x/d" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # For the Mogi model, Ur/Uz = x/d (from Eqs. 1-2)
    a_val = 1000.0
    d_val = 4000.0

    sys = MogiModel()
    compiled = mtkcompile(sys)

    for x_val in [1000.0, 3000.0, 5000.0, 10000.0, 20000.0]
        params = [
            compiled.G => G_val, compiled.ν => ν_val,
            compiled.ΔP => ΔP_val, compiled.a => a_val,
            compiled.d => d_val, compiled.x => x_val,
        ]
        sol = solve(make_prob(compiled, params))
        @test sol[compiled.Ur] / sol[compiled.Uz] ≈ x_val / d_val rtol = 1.0e-10
    end
end

@testitem "Far-Field Decay" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # Displacements should decrease with increasing distance
    # For large x, Ur ~ x/R^3 ~ 1/x^2 and Uz ~ d/R^3 ~ 1/x^3 approximately

    a_val = 1000.0
    d_val = 4000.0

    sys = MogiModel()
    compiled = mtkcompile(sys)

    distances = [5000.0, 10000.0, 15000.0, 20000.0, 25000.0]
    Ur_vals = Float64[]
    Uz_vals = Float64[]

    for x_val in distances
        params = [
            compiled.G => G_val, compiled.ν => ν_val,
            compiled.ΔP => ΔP_val, compiled.a => a_val,
            compiled.d => d_val, compiled.x => x_val,
        ]
        sol = solve(make_prob(compiled, params))
        push!(Ur_vals, sol[compiled.Ur])
        push!(Uz_vals, sol[compiled.Uz])
    end

    # Both should be monotonically decreasing
    for i in 1:(length(distances) - 1)
        @test Ur_vals[i] > Ur_vals[i + 1]
        @test Uz_vals[i] > Uz_vals[i + 1]
    end
end

@testitem "Kīlauea Deflation Event (Table 3)" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # Test with Table 3 parameters: Kīlauea June 2007 deflation event
    # ΔP = -3.73 MPa (negative = deflation), G = 4.0 GPa, ν = 0.25
    ΔP_kil = -3.73e6  # Pa
    a_val = 800.0      # m (within 0.40-1.10 km range)
    d_val = 2000.0     # m (within 1.00-3.00 km range)
    x_val = 3000.0     # m

    # Mogi model with deflation
    sys = MogiModel()
    compiled = mtkcompile(sys)

    params = [
        compiled.G => G_val, compiled.ν => ν_val,
        compiled.ΔP => ΔP_kil, compiled.a => a_val,
        compiled.d => d_val, compiled.x => x_val,
    ]
    sol = solve(make_prob(compiled, params))

    # Negative overpressure (deflation) → negative displacements
    @test sol[compiled.Ur] < 0  # inward radial displacement
    @test sol[compiled.Uz] < 0  # subsidence

    # Verify against precomputed reference values from Eqs. 1-2 with Table 3 params
    @test sol[compiled.Ur] ≈ -0.022918505338191925 rtol = 1.0e-10
    @test sol[compiled.Uz] ≈ -0.015279003558794616 rtol = 1.0e-10

    # McTigue with same parameters
    mctigue_sys = McTigueModel()
    mctigue_compiled = mtkcompile(mctigue_sys)

    params_mc = [
        mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
        mctigue_compiled.ΔP => ΔP_kil, mctigue_compiled.a => a_val,
        mctigue_compiled.d => d_val, mctigue_compiled.x => x_val,
    ]
    sol_mc = solve(make_prob(mctigue_compiled, params_mc))

    # McTigue also negative for deflation
    @test sol_mc[mctigue_compiled.Ur] < 0
    @test sol_mc[mctigue_compiled.Uz] < 0

    # McTigue magnitudes should be larger than Mogi for ε = 0.4
    @test abs(sol_mc[mctigue_compiled.Ur]) > abs(sol[compiled.Ur])
    @test abs(sol_mc[mctigue_compiled.Uz]) > abs(sol[compiled.Uz])
end

@testitem "Mogi Ur Maximum at x_crit" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # The Mogi model predicts Ur is maximized at x_crit = d/sqrt(2)
    # This follows from setting dUr/dx = 0 in Eq. 1
    a_val = 1000.0
    d_val = 4000.0

    sys = MogiModel()
    compiled = mtkcompile(sys)

    x_crit = d_val / sqrt(2)

    # Evaluate at x_crit and nearby points
    params_crit = [
        compiled.G => G_val, compiled.ν => ν_val,
        compiled.ΔP => ΔP_val, compiled.a => a_val,
        compiled.d => d_val, compiled.x => x_crit,
    ]
    sol_crit = solve(make_prob(compiled, params_crit))
    Ur_max = sol_crit[compiled.Ur]

    # Check that Ur_max matches precomputed reference value
    @test Ur_max ≈ 0.009021097956087904 rtol = 1.0e-10

    # Verify Ur at x_crit is greater than at nearby points
    for x_offset in [x_crit - 500.0, x_crit + 500.0]
        params_off = [
            compiled.G => G_val, compiled.ν => ν_val,
            compiled.ΔP => ΔP_val, compiled.a => a_val,
            compiled.d => d_val, compiled.x => x_offset,
        ]
        sol_off = solve(make_prob(compiled, params_off))
        @test sol_off[compiled.Ur] < Ur_max
    end
end

@testitem "Mogi Uz_max at x=0" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # At x≈0, Uz = a^3 * ΔP * (1-ν) / (G * d^2)
    # Precomputed reference value for a=1km, d=4km, Table 1 params
    a_val = 1000.0
    d_val = 4000.0

    sys = MogiModel()
    compiled = mtkcompile(sys)

    params = [
        compiled.G => G_val, compiled.ν => ν_val,
        compiled.ΔP => ΔP_val, compiled.a => a_val,
        compiled.d => d_val, compiled.x => 0.001,  # near x=0
    ]
    sol = solve(make_prob(compiled, params))

    # Uz_max = a^3 * ΔP * (1-ν) / (G * d^2) = 0.0234375 m
    @test sol[compiled.Uz] ≈ 0.0234375 rtol = 1.0e-6
end

@testitem "McTigue Correction Magnitude" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # Verify the McTigue correction factor magnitude for specific ε values
    # For ν=0.25, ε=0.25 (a=1km, d=4km), the correction should be ≈1.0053
    a_val = 1000.0
    d_val = 4000.0
    x_val = 5000.0

    mogi_sys = MogiModel()
    mogi_compiled = mtkcompile(mogi_sys)
    mctigue_sys = McTigueModel()
    mctigue_compiled = mtkcompile(mctigue_sys)

    params_mogi = [
        mogi_compiled.G => G_val, mogi_compiled.ν => ν_val,
        mogi_compiled.ΔP => ΔP_val, mogi_compiled.a => a_val,
        mogi_compiled.d => d_val, mogi_compiled.x => x_val,
    ]
    params_mctigue = [
        mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
        mctigue_compiled.ΔP => ΔP_val, mctigue_compiled.a => a_val,
        mctigue_compiled.d => d_val, mctigue_compiled.x => x_val,
    ]

    sol_mogi = solve(make_prob(mogi_compiled, params_mogi))
    sol_mctigue = solve(make_prob(mctigue_compiled, params_mctigue))

    # McTigue/Mogi ratio should equal the correction factor ≈ 1.005261
    correction = sol_mctigue[mctigue_compiled.Ur] / sol_mogi[mogi_compiled.Ur]
    @test correction ≈ 1.0052608032873807 rtol = 1.0e-10

    # Same correction applies to Uz
    correction_uz = sol_mctigue[mctigue_compiled.Uz] / sol_mogi[mogi_compiled.Uz]
    @test correction_uz ≈ correction rtol = 1.0e-10
end

@testitem "McTigue Always Greater Than Mogi for ε > 0" setup = [MogiMcTigueSetup] tags = [:mogi_mctigue] begin
    # For ν = 0.25, the McTigue correction is > 1, so McTigue displacements
    # should be slightly larger than Mogi displacements

    d_val = 4000.0
    x_val = 5000.0

    mogi_sys = MogiModel()
    mogi_compiled = mtkcompile(mogi_sys)
    mctigue_sys = McTigueModel()
    mctigue_compiled = mtkcompile(mctigue_sys)

    for a_val in [500.0, 1000.0, 1500.0, 2000.0]
        params_mogi = [
            mogi_compiled.G => G_val, mogi_compiled.ν => ν_val,
            mogi_compiled.ΔP => ΔP_val, mogi_compiled.a => a_val,
            mogi_compiled.d => d_val, mogi_compiled.x => x_val,
        ]
        params_mctigue = [
            mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
            mctigue_compiled.ΔP => ΔP_val, mctigue_compiled.a => a_val,
            mctigue_compiled.d => d_val, mctigue_compiled.x => x_val,
        ]

        sol_mogi = solve(make_prob(mogi_compiled, params_mogi))
        sol_mctigue = solve(make_prob(mctigue_compiled, params_mctigue))

        # McTigue correction is positive for ν=0.25, so magnitudes should be larger
        @test abs(sol_mctigue[mctigue_compiled.Ur]) > abs(sol_mogi[mogi_compiled.Ur])
        @test abs(sol_mctigue[mctigue_compiled.Uz]) > abs(sol_mogi[mogi_compiled.Uz])
    end
end
