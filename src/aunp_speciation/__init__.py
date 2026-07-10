"""AuNP speciation: monomer/dimer/trimer optical model for UV-Vis broadening.

Layers:
  dielectric  - gold permittivity (Etchegoin analytic; swap for J&C in production)
  mie         - single-sphere Mie cross sections + dipole polarizability
  clusters    - coupled-dipole dimer/trimer spectra with orientation averaging
  equilibrium - monomer<->dimer<->trimer association vs temperature (isosbestic)
  spectra     - polydispersity + population-weighted ensemble spectra
"""
from . import dielectric, mie, clusters, equilibrium, spectra  # noqa: F401

__all__ = ["dielectric", "mie", "clusters", "equilibrium", "spectra"]
