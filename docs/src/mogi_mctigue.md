# Mogi & McTigue Volcanic Deformation Models

## Overview

The Mogi (1958) and McTigue (1987) models are analytical solutions for computing
surface displacements caused by a pressurized spherical magma chamber embedded
in an elastic, homogeneous half-space with a flat free surface. These models are
widely used in volcano geodesy to interpret ground deformation observations and
to infer magma reservoir properties (location, size, and overpressure).

The **Mogi model** approximates the magma chamber as a point source (small sphere
approximation, ``\varepsilon = a/d \ll 1``), where ``a`` is the source radius
and ``d`` is the source depth. The **McTigue model** extends Mogi by incorporating
higher-order correction terms that account for the finite dimensions of the
spherical source, providing improved accuracy for larger ``\varepsilon`` values.

Taylor et al. (2021) demonstrated that:
- The Mogi model agrees with Finite Element models within 5% when ``\varepsilon < 0.37``
- The McTigue model agrees within 5% when ``\varepsilon < 0.59``

**Reference**: Taylor, N.C., Johnson, J.H., Herd, R.A. (2021). "Making the most
of the Mogi model: Size matters." *Journal of Volcanology and Geothermal Research*,
419, 107380. [https://doi.org/10.1016/j.jvolgeores.2021.107380](https://doi.org/10.1016/j.jvolgeores.2021.107380)

```@docs
MogiModel
McTigueModel
```

## Implementation

The models compute horizontal-radial (``U_r``) and vertical (``U_z``) surface
displacements as functions of horizontal distance ``x`` from the source.

### State Variables

```@example mogi_mctigue
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using Geodynamics

sys = MogiModel()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example mogi_mctigue
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [let u = ModelingToolkit.get_unit(p); u === nothing ? "dimensionless" : dimension(u) end for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

**Mogi Model** (Equations 1-2 from Taylor et al., 2021):

```@example mogi_mctigue
equations(sys)
```

**McTigue Model** (Equations 3-4 from Taylor et al., 2021):

```@example mogi_mctigue
mctigue_sys = McTigueModel()
equations(mctigue_sys)
```

## Analysis

### Displacement Profiles (cf. Figure 1 geometry)

Surface displacement profiles for both models, using the reference parameters from
Table 1 of Taylor et al. (2021): ``G = 4.0`` GPa, ``\nu = 0.25``, ``\Delta P = 2.0`` MPa.

```@example mogi_mctigue
using CairoMakie
using ModelingToolkit: mtkcompile

# Compile both models
mogi_compiled = mtkcompile(MogiModel())
mctigue_compiled = mtkcompile(McTigueModel())

# Parameters from Table 1
G_val = 4.0e9    # Pa
ν_val = 0.25
ΔP_val = 2.0e6   # Pa
d_val = 4000.0    # m (4 km depth)

# Evaluate displacement profiles for different source radii
x_range = range(200, 25000, length=200)

fig = Figure(size=(900, 600))
ax_ur = Axis(fig[1, 1], xlabel="Distance x (km)", ylabel="Ur (mm)",
    title="Horizontal-Radial Displacement")
ax_uz = Axis(fig[2, 1], xlabel="Distance x (km)", ylabel="Uz (mm)",
    title="Vertical Displacement")

colors_mogi = [:blue, :red, :green]
colors_mctigue = [:cyan, :orange, :lightgreen]
radii = [500.0, 1000.0, 2000.0]
labels_a = ["a = 0.5 km", "a = 1.0 km", "a = 2.0 km"]

for (i, a_val) in enumerate(radii)
    Ur_mogi = Float64[]
    Uz_mogi = Float64[]
    Ur_mctigue = Float64[]
    Uz_mctigue = Float64[]

    for x_val in x_range
        # Mogi
        prob_m = NonlinearProblem(mogi_compiled, Dict([
            mogi_compiled.G => G_val, mogi_compiled.ν => ν_val,
            mogi_compiled.ΔP => ΔP_val, mogi_compiled.a => a_val,
            mogi_compiled.d => d_val, mogi_compiled.x => x_val,
        ]))
        sol_m = solve(prob_m)
        push!(Ur_mogi, sol_m[mogi_compiled.Ur] * 1000)  # convert to mm
        push!(Uz_mogi, sol_m[mogi_compiled.Uz] * 1000)

        # McTigue
        prob_mc = NonlinearProblem(mctigue_compiled, Dict([
            mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
            mctigue_compiled.ΔP => ΔP_val, mctigue_compiled.a => a_val,
            mctigue_compiled.d => d_val, mctigue_compiled.x => x_val,
        ]))
        sol_mc = solve(prob_mc)
        push!(Ur_mctigue, sol_mc[mctigue_compiled.Ur] * 1000)
        push!(Uz_mctigue, sol_mc[mctigue_compiled.Uz] * 1000)
    end

    x_km = collect(x_range) ./ 1000
    lines!(ax_ur, x_km, Ur_mogi, color=colors_mogi[i], label="Mogi " * labels_a[i])
    lines!(ax_ur, x_km, Ur_mctigue, color=colors_mctigue[i], linestyle=:dash,
        label="McTigue " * labels_a[i])
    lines!(ax_uz, x_km, Uz_mogi, color=colors_mogi[i], label="Mogi " * labels_a[i])
    lines!(ax_uz, x_km, Uz_mctigue, color=colors_mctigue[i], linestyle=:dash,
        label="McTigue " * labels_a[i])
end

axislegend(ax_ur, position=:rt)
axislegend(ax_uz, position=:rt)
fig
```

### McTigue Correction Magnitude vs. ε (cf. Figure 2)

The percentage difference between Mogi and McTigue predictions increases with
``\varepsilon = a/d``, illustrating the growing importance of finite-body corrections.

```@example mogi_mctigue
fig2 = Figure(size=(700, 500))
ax = Axis(fig2[1, 1],
    xlabel="ε = a/d",
    ylabel="(McTigue - Mogi) / Mogi × 100 (%)",
    title="McTigue vs Mogi: Relative Difference at Ur_max and Uz_max")

depths = [1000.0, 4000.0, 7000.0, 10000.0, 13000.0]
depth_colors = [:purple, :blue, :green, :orange, :gold]
depth_labels = ["d = 1 km", "d = 4 km", "d = 7 km", "d = 10 km", "d = 13 km"]

ε_range = range(0.05, 0.70, length=50)

for (j, d_val) in enumerate(depths)
    pct_diff_ur = Float64[]
    pct_diff_uz = Float64[]

    for ε in ε_range
        a_val = ε * d_val

        # Evaluate at x = x_crit for Mogi (where Ur is maximum): x_crit = d/sqrt(2)
        x_crit = d_val / sqrt(2)

        # Mogi at x_crit
        prob_m = NonlinearProblem(mogi_compiled, Dict([
            mogi_compiled.G => G_val, mogi_compiled.ν => ν_val,
            mogi_compiled.ΔP => ΔP_val, mogi_compiled.a => a_val,
            mogi_compiled.d => d_val, mogi_compiled.x => x_crit,
        ]))
        sol_m = solve(prob_m)

        # McTigue at x_crit
        prob_mc = NonlinearProblem(mctigue_compiled, Dict([
            mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
            mctigue_compiled.ΔP => ΔP_val, mctigue_compiled.a => a_val,
            mctigue_compiled.d => d_val, mctigue_compiled.x => x_crit,
        ]))
        sol_mc = solve(prob_mc)

        push!(pct_diff_ur,
            (sol_mc[mctigue_compiled.Ur] - sol_m[mogi_compiled.Ur]) /
            sol_m[mogi_compiled.Ur] * 100)
        push!(pct_diff_uz,
            (sol_mc[mctigue_compiled.Uz] - sol_m[mogi_compiled.Uz]) /
            sol_m[mogi_compiled.Uz] * 100)
    end

    lines!(ax, collect(ε_range), pct_diff_ur, color=depth_colors[j],
        label=depth_labels[j] * " (Ur)")
    lines!(ax, collect(ε_range), pct_diff_uz, color=depth_colors[j],
        linestyle=:dash, label=depth_labels[j] * " (Uz)")
end

axislegend(ax, position=:lt)
fig2
```

### Displacement as a Function of Depth (cf. Figure 4)

The effect of source depth on median differences between Mogi and McTigue models.

```@example mogi_mctigue
fig3 = Figure(size=(800, 400))
ax_ur3 = Axis(fig3[1, 1], xlabel="ε", ylabel="ΔUr (%)",
    title="Mogi-McTigue Ur Difference")
ax_uz3 = Axis(fig3[1, 2], xlabel="ε", ylabel="ΔUz (%)",
    title="Mogi-McTigue Uz Difference")

for (j, d_val) in enumerate(depths)
    diff_ur = Float64[]
    diff_uz = Float64[]

    for ε in ε_range
        a_val = ε * d_val
        x_vals = range(200, 25000, length=50)
        ur_diffs = Float64[]
        uz_diffs = Float64[]

        for x_val in x_vals
            prob_m = NonlinearProblem(mogi_compiled, Dict([
                mogi_compiled.G => G_val, mogi_compiled.ν => ν_val,
                mogi_compiled.ΔP => ΔP_val, mogi_compiled.a => a_val,
                mogi_compiled.d => d_val, mogi_compiled.x => x_val,
            ]))
            sol_m = solve(prob_m)

            prob_mc = NonlinearProblem(mctigue_compiled, Dict([
                mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
                mctigue_compiled.ΔP => ΔP_val, mctigue_compiled.a => a_val,
                mctigue_compiled.d => d_val, mctigue_compiled.x => x_val,
            ]))
            sol_mc = solve(prob_mc)

            ur_m = sol_m[mogi_compiled.Ur]
            ur_mc = sol_mc[mctigue_compiled.Ur]
            uz_m = sol_m[mogi_compiled.Uz]
            uz_mc = sol_mc[mctigue_compiled.Uz]

            if abs(ur_m) > 1e-15
                push!(ur_diffs, abs(ur_mc - ur_m) / abs(ur_m) * 100)
            end
            if abs(uz_m) > 1e-15
                push!(uz_diffs, abs(uz_mc - uz_m) / abs(uz_m) * 100)
            end
        end

        push!(diff_ur, isempty(ur_diffs) ? 0.0 : median(ur_diffs))
        push!(diff_uz, isempty(uz_diffs) ? 0.0 : median(uz_diffs))
    end

    lines!(ax_ur3, collect(ε_range), diff_ur, color=depth_colors[j],
        label=depth_labels[j])
    lines!(ax_uz3, collect(ε_range), diff_uz, color=depth_colors[j],
        label=depth_labels[j])
end

axislegend(ax_ur3, position=:lt)
axislegend(ax_uz3, position=:lt)
fig3
```
