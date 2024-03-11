# Microfacet Models for Refraction through Rough Surfaces

## Microfacet Theory

### A BSDF describes how light scatters from a surface. It's defined as the ratio of scattered radiance in a direction `o` caused per unit irradiance incident from direction `i` denoted as `fs(i, o, n)`. It's dependent on the local surface normal `n`.

### BSDF is the sum of BRDF, `fr`, and BTDF, `ft`. So a BSDF can correctly handle directions directions on either side of the surface.

### A detailed mcrosurface is replaced by a simplified macrosurface with a modified scattering function (BSDF) that matches the aggregate directional scattering of the microsurface.

### The microsurface can be described by two statistical measures, a microfacet distribution function D and a shadowing-masking function G and a microsurface BSDF `fsm`

### Microfacet Distribution Function `D`

#### The microfacet normal distribution `D(m)` describes the statistical distribution of surface normals `m` over the microsurface. Given an solid angal `dwm` centered on `m`, and a macrosurface area `dA`, `D(m)dwmdA` is the total area of the portion of the corresponding microsurface whose normals lie within that specified solid angle.

#### A plausible microfacet distribtion should obey at leas the following properties:

#### Microfacet density is positive
$$
D(m) > 0
$$

#### Total microsurface area is at least as large as the corresponding macrosurface's area
$$
1 \leq \int D(m) \, d\omega_m
$$

#### The (signed) projected area of the microsurface is the same as the projected area of the macrosurface for any direction `v`
$$
(v \cdot n) = \int D(m) (v \cdot m) d\omega_m
$$