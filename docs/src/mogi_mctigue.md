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
param_units = String[]
for p in params
    u = ModelingToolkit.get_unit(p)
    push!(param_units, u === nothing ? "dimensionless" : string(dimension(u)))
end
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => param_units,
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

### Reference Parameters (Table 1, Taylor et al., 2021)

The following parameters, based on Kilauea Volcano, Hawaii, are used in the
analysis below:

| Symbol          | Definition                     | Units | Value        |
|:--------------- |:------------------------------ |:----- |:------------ |
| ``G``           | Shear modulus                  | GPa   | 4.0          |
| ``\nu``         | Poisson's ratio                | --    | 0.25         |
| ``\Delta P``    | Overpressure                   | MPa   | 2.0          |
| ``a``           | Source radius                  | km    | variable     |
| ``d``           | Source depth                   | km    | variable     |
| ``\varepsilon`` | ``a/d``                        | --    | 0.05 -- 0.70 |
| ``U_r``         | Horizontal-radial displacement | m     | result       |
| ``U_z``         | Vertical displacement          | m     | result       |

## Analysis

### Displacement Profiles (cf. Figure 1 geometry)

Surface displacement profiles for both models, using the reference parameters from
Table 1 of Taylor et al. (2021): ``G = 4.0`` GPa, ``\nu = 0.25``, ``\Delta P = 2.0`` MPa.

```@example mogi_mctigue
using CairoMakie
using ModelingToolkit: mtkcompile
using NonlinearSolve

# Compile both models
mogi_compiled = mtkcompile(MogiModel())
mctigue_compiled = mtkcompile(McTigueModel())

# Parameters from Table 1
G_val = 4.0e9    # Pa
ν_val = 0.25
ΔP_val = 2.0e6   # Pa
d_val = 4000.0    # m (4 km depth)

# Evaluate displacement profiles for different source radii
x_range = range(200, 25000, length = 200)

fig = Figure(size = (900, 600))
ax_ur = Axis(fig[1, 1], xlabel = "Distance x (km)", ylabel = "Ur (mm)",
    title = "Horizontal-Radial Displacement")
ax_uz = Axis(fig[2, 1], xlabel = "Distance x (km)", ylabel = "Uz (mm)",
    title = "Vertical Displacement")

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
        prob_m = NonlinearProblem(mogi_compiled,
            Dict([
                mogi_compiled.G => G_val, mogi_compiled.ν => ν_val,
                mogi_compiled.ΔP => ΔP_val, mogi_compiled.a => a_val,
                mogi_compiled.d => d_val, mogi_compiled.x => x_val]))
        sol_m = solve(prob_m)
        push!(Ur_mogi, sol_m[mogi_compiled.Ur] * 1000)  # convert to mm
        push!(Uz_mogi, sol_m[mogi_compiled.Uz] * 1000)

        # McTigue
        prob_mc = NonlinearProblem(mctigue_compiled,
            Dict([
                mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
                mctigue_compiled.ΔP => ΔP_val, mctigue_compiled.a => a_val,
                mctigue_compiled.d => d_val, mctigue_compiled.x => x_val]))
        sol_mc = solve(prob_mc)
        push!(Ur_mctigue, sol_mc[mctigue_compiled.Ur] * 1000)
        push!(Uz_mctigue, sol_mc[mctigue_compiled.Uz] * 1000)
    end

    x_km = collect(x_range) ./ 1000
    lines!(ax_ur, x_km, Ur_mogi, color = colors_mogi[i], label = "Mogi " * labels_a[i])
    lines!(ax_ur, x_km, Ur_mctigue, color = colors_mctigue[i], linestyle = :dash,
        label = "McTigue " * labels_a[i])
    lines!(ax_uz, x_km, Uz_mogi, color = colors_mogi[i], label = "Mogi " * labels_a[i])
    lines!(ax_uz, x_km, Uz_mctigue, color = colors_mctigue[i], linestyle = :dash,
        label = "McTigue " * labels_a[i])
end

axislegend(ax_ur, position = :rt)
axislegend(ax_uz, position = :rt)
fig
```

### McTigue Correction Magnitude vs. epsilon (cf. Figure 2)

The percentage difference between Mogi and McTigue predictions increases with
``\varepsilon = a/d``, illustrating the growing importance of finite-body corrections.

```@example mogi_mctigue
fig2 = Figure(size = (700, 500))
ax = Axis(fig2[1, 1],
    xlabel = "epsilon = a/d",
    ylabel = "(McTigue - Mogi) / Mogi x 100 (%)",
    title = "McTigue vs Mogi: Relative Difference at Ur_max and Uz_max")

depths = [1000.0, 4000.0, 7000.0, 10000.0, 13000.0]
depth_colors = [:purple, :blue, :green, :orange, :gold]
depth_labels = ["d = 1 km", "d = 4 km", "d = 7 km", "d = 10 km", "d = 13 km"]

eps_range = range(0.05, 0.70, length = 50)

for (j, d_val) in enumerate(depths)
    pct_diff_ur = Float64[]
    pct_diff_uz = Float64[]

    for eps_val in eps_range
        a_val = eps_val * d_val

        # Evaluate at x = x_crit for Mogi (where Ur is maximum): x_crit = d/sqrt(2)
        x_crit = d_val / sqrt(2)

        # Mogi at x_crit
        prob_m = NonlinearProblem(mogi_compiled,
            Dict([
                mogi_compiled.G => G_val, mogi_compiled.ν => ν_val,
                mogi_compiled.ΔP => ΔP_val, mogi_compiled.a => a_val,
                mogi_compiled.d => d_val, mogi_compiled.x => x_crit]))
        sol_m = solve(prob_m)

        # McTigue at x_crit
        prob_mc = NonlinearProblem(mctigue_compiled,
            Dict([
                mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
                mctigue_compiled.ΔP => ΔP_val, mctigue_compiled.a => a_val,
                mctigue_compiled.d => d_val, mctigue_compiled.x => x_crit]))
        sol_mc = solve(prob_mc)

        push!(pct_diff_ur,
            (sol_mc[mctigue_compiled.Ur] - sol_m[mogi_compiled.Ur]) /
            sol_m[mogi_compiled.Ur] * 100)
        push!(pct_diff_uz,
            (sol_mc[mctigue_compiled.Uz] - sol_m[mogi_compiled.Uz]) /
            sol_m[mogi_compiled.Uz] * 100)
    end

    lines!(ax, collect(eps_range), pct_diff_ur, color = depth_colors[j],
        label = depth_labels[j] * " (Ur)")
    lines!(ax, collect(eps_range), pct_diff_uz, color = depth_colors[j],
        linestyle = :dash, label = depth_labels[j] * " (Uz)")
end

axislegend(ax, position = :lt)
fig2
```

### Maximum Displacement vs. Source Depth

The maximum vertical displacement ``U_{z,max}`` (at ``x=0``) decreases rapidly with
increasing source depth for both models.

```@example mogi_mctigue
fig3 = Figure(size = (700, 500))
ax3 = Axis(fig3[1, 1], xlabel = "Source depth d (km)", ylabel = "Uz_max (mm)",
    title = "Maximum Vertical Displacement vs. Depth")

d_range = range(1000, 13000, length = 50)
radii = [500.0, 1000.0, 2000.0]
colors_m = [:blue, :red, :green]
colors_mc = [:cyan, :orange, :lightgreen]
labels_a = ["a = 0.5 km", "a = 1.0 km", "a = 2.0 km"]

for (i, a_val) in enumerate(radii)
    Uz_max_mogi = Float64[]
    Uz_max_mctigue = Float64[]

    for d_val in d_range
        x_near = 1.0  # near x=0

        prob_m = NonlinearProblem(mogi_compiled,
            Dict([
                mogi_compiled.G => G_val, mogi_compiled.ν => ν_val,
                mogi_compiled.ΔP => ΔP_val, mogi_compiled.a => a_val,
                mogi_compiled.d => d_val, mogi_compiled.x => x_near]))
        sol_m = solve(prob_m)
        push!(Uz_max_mogi, sol_m[mogi_compiled.Uz] * 1000)

        prob_mc = NonlinearProblem(mctigue_compiled,
            Dict([
                mctigue_compiled.G => G_val, mctigue_compiled.ν => ν_val,
                mctigue_compiled.ΔP => ΔP_val, mctigue_compiled.a => a_val,
                mctigue_compiled.d => d_val, mctigue_compiled.x => x_near]))
        sol_mc = solve(prob_mc)
        push!(Uz_max_mctigue, sol_mc[mctigue_compiled.Uz] * 1000)
    end

    d_km = collect(d_range) ./ 1000
    lines!(ax3, d_km, Uz_max_mogi, color = colors_m[i],
        label = "Mogi " * labels_a[i])
    lines!(ax3, d_km, Uz_max_mctigue, color = colors_mc[i], linestyle = :dash,
        label = "McTigue " * labels_a[i])
end

axislegend(ax3, position = :rt)
fig3
```
