# HADDOCK3 2026.7.0 seed compatibility correction

The supplied `PROTOCOL_LOCK.md` places `seed = 42` at the top level. A
setup-only parse with the pinned HADDOCK3 2026.7.0 rejects that key because
this version has no general `seed` parameter. Its CNS modules expose the
version-specific parameter `iniseed` instead.

To preserve the intended reproducibility constraint in an executable config,
`guided_ppil4_haddock.cfg` pins `iniseed = 42` independently in `topoaa`,
`rigidbody`, `flexref`, and `emref`. No sampling, scoring, restraint, topology,
or input setting was otherwise changed. This compatibility correction is
included in the protocol digest and must be byte-identical on both machines.
