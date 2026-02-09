"""
    Mogi and McTigue Volcanic Deformation Models

ModelingToolkit.jl implementation of the Mogi (1958) and McTigue (1987) analytical
models for volcanic surface deformation, as described in:

Taylor, N.C., Johnson, J.H., Herd, R.A. (2021). "Making the most of the Mogi model:
Size matters." Journal of Volcanology and Geothermal Research, 419, 107380.

These models compute surface displacements (horizontal-radial and vertical) due to
a pressurized spherical magma chamber embedded in an elastic, homogeneous half-space
with a flat free surface.

The Mogi model approximates the source as a point (small sphere, ε = a/d ≪ 1).
The McTigue model extends Mogi by incorporating higher-order corrections for the
finite radius of a spherical source.

**Equations implemented:**
- Equations 1-2: Mogi forward model (horizontal-radial and vertical displacement)
- Equations 3-4: McTigue forward model (with finite-sphere corrections)
"""

export MogiModel, McTigueModel

"""
$(TYPEDSIGNATURES)

Mogi (1958) point-source model for volcanic surface deformation.

Computes horizontal-radial (`Ur`) and vertical (`Uz`) surface displacements
due to a pressurized spherical cavity in an elastic half-space, using the
small-sphere (point-source) approximation (ε = a/d ≪ 1).

**Reference**: Mogi, K., 1958. Relations between the eruptions of various
volcanoes and the deformations of the ground surfaces around them. Bull.
Earthq. Res. Inst. 36, 99–134. As presented in Taylor et al. (2021),
Eqs. 1–2.

**Parameters:**

  - `G`: Shear modulus (Pa)
  - `ν`: Poisson's ratio (dimensionless)
  - `ΔP`: Source overpressure (Pa)
  - `a`: Source radius (m)
  - `d`: Source depth (m)
  - `x`: Horizontal distance from source (m)
"""
@component function MogiModel(; name = :MogiModel)
    @constants begin
        one_m = 1.0, [description = "Unit length for non-dimensionalization", unit = u"m"]
    end

    @parameters begin
        G, [description = "Shear modulus", unit = u"Pa"]
        ν, [description = "Poisson's ratio (dimensionless)"]
        ΔP, [description = "Source overpressure", unit = u"Pa"]
        a, [description = "Source radius", unit = u"m"]
        d, [description = "Source depth", unit = u"m"]
        x, [description = "Horizontal distance from source", unit = u"m"]
    end

    @variables begin
        R(t), [description = "Distance from source to surface point", unit = u"m"]
        Ur(t), [description = "Horizontal-radial surface displacement", unit = u"m"]
        Uz(t), [description = "Vertical surface displacement", unit = u"m"]
    end

    eqs = [
        # Auxiliary: Distance from source center to surface point
        R ~ (x^2 / one_m^2 + d^2 / one_m^2)^0.5 * one_m,

        # Eq. 1 — Mogi horizontal-radial displacement
        Ur ~ a^3 * ΔP * (1 - ν) * x / (G * R^3),

        # Eq. 2 — Mogi vertical displacement
        Uz ~ a^3 * ΔP * (1 - ν) * d / (G * R^3)
    ]

    return System(eqs, t; name)
end

"""
$(TYPEDSIGNATURES)

McTigue (1987) finite-sphere model for volcanic surface deformation.

Extends the Mogi model by incorporating higher-order correction terms that
account for the finite radius of a spherical deformation source. The correction
terms depend on ε = a/d and provide improved accuracy for larger ε values.

The McTigue model agrees with FE models within 5% for ε < 0.59, compared to
ε < 0.37 for the Mogi model (Taylor et al., 2021).

**Reference**: McTigue, D.F., 1987. Elastic stress and deformation near a finite
spherical magma body: resolution of the point source paradox. J. Geophys. Res.
92 (B12), 12931. As presented in Taylor et al. (2021), Eqs. 3–4.

**Parameters:**

  - `G`: Shear modulus (Pa)
  - `ν`: Poisson's ratio (dimensionless)
  - `ΔP`: Source overpressure (Pa)
  - `a`: Source radius (m)
  - `d`: Source depth (m)
  - `x`: Horizontal distance from source (m)
"""
@component function McTigueModel(; name = :McTigueModel)
    @constants begin
        one_m = 1.0, [description = "Unit length for non-dimensionalization", unit = u"m"]
    end

    @parameters begin
        G, [description = "Shear modulus", unit = u"Pa"]
        ν, [description = "Poisson's ratio (dimensionless)"]
        ΔP, [description = "Source overpressure", unit = u"Pa"]
        a, [description = "Source radius", unit = u"m"]
        d, [description = "Source depth", unit = u"m"]
        x, [description = "Horizontal distance from source", unit = u"m"]
    end

    @variables begin
        R(t), [description = "Distance from source to surface point", unit = u"m"]
        Ur(t), [description = "Horizontal-radial surface displacement", unit = u"m"]
        Uz(t), [description = "Vertical surface displacement", unit = u"m"]
    end

    # McTigue correction factor (Eqs. 3-4):
    #   correction = 1 + (a/d)^3 * ((1+ν)/(2(-7+5ν)) + (15d^2(-2+ν))/(4R^2(-7+5ν)))
    eqs = [
        # Auxiliary: Distance from source center to surface point
        R ~ (x^2 / one_m^2 + d^2 / one_m^2)^0.5 * one_m,

        # Eq. 3 — McTigue horizontal-radial displacement
        Ur ~
        a^3 * ΔP * (1 - ν) * x / (G * R^3) * (
            1 +
            (a / d)^3 * (
            (1 + ν) / (2 * (-7 + 5ν)) +
            15 * d^2 * (-2 + ν) / (4 * R^2 * (-7 + 5ν))
        )
        ),

        # Eq. 4 — McTigue vertical displacement
        Uz ~
        a^3 * ΔP * (1 - ν) * d / (G * R^3) * (
            1 +
            (a / d)^3 * (
            (1 + ν) / (2 * (-7 + 5ν)) +
            15 * d^2 * (-2 + ν) / (4 * R^2 * (-7 + 5ν))
        )
        )
    ]

    return System(eqs, t; name)
end
