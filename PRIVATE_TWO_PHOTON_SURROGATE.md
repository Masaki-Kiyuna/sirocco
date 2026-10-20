# Private H I two-photon line surrogate

Status: incomplete diagnostic model; not validated production physics.

The private collapsed-hydrogen atomic-data variants containing
`h20_2ph_lines.dat` represent the H I 2s -> 1s two-photon decay with an
effective rate of `8.226 / 4 = 2.0565 s^-1`.  The factor 1/4 is the statistical
2s fraction of the collapsed n=2 level.

This is only a line surrogate.  The macro-atom emits one packet carrying the
full transition energy at the Ly-alpha frequency.  The code bypasses the
Sobolev Ly-alpha escape probability for this tagged transition, but it does
not generate two photons or sample the physical two-photon continuum.
Therefore it does not correctly model:

- the two-photon spectral shape;
- frequency-dependent propagation and reabsorption, including absorption in
  the Balmer continuum;
- the division between true Ly-alpha and two-photon emergent luminosity; or
- feedback of the redistributed continuum on ionization and thermal balance.

The runtime switch is:

```text
@Hydrogen.two_photon_emission(on,off) off
```

`off` sets the tagged surrogate transition rate to zero while leaving true
Ly-alpha and all other transitions unchanged.  Use `off` for production runs
until a continuum emission channel has been implemented and validated with
rate, energy-conservation, transport, and spectrum tests.
