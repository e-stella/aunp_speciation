"""Layer 1: monomer <-> dimer <-> trimer association equilibrium.

Schematic thermodynamic model. Monomer number concentration c1 (in units of a
reference concentration c0), with

    c2 = K2 * c1**2        (dimer)
    c3 = K3 * c1**3        (trimer)

and total-gold conservation (in monomer equivalents):

    C_tot = c1 + 2*c2 + 3*c3.

Temperature enters through van 't Hoff:  K(T) = exp(-(dH - T dS)/(R T)),
so association constants change with T and the mix of species shifts. When the
total spectrum is a population-weighted sum of species spectra that share
crossing wavelengths, this naturally produces an isosbestic point.

Concentrations are dimensionless (relative to c0); absolute values are absorbed
into K. This is a demonstration model, not fitted physical constants.
"""

from __future__ import annotations
import numpy as np
from scipy.optimize import brentq

R_KJ = 8.314462618e-3  # kJ / (mol K)


def association_constants(T, dH2=-40.0, dS2=-0.10, dH3=-70.0, dS3=-0.20):
    """van 't Hoff K2(T), K3(T). dH in kJ/mol, dS in kJ/(mol K).

    Defaults: exothermic association (dH<0) with entropy penalty (dS<0), i.e.
    higher T favors monomers -> aggregation decreases with heating (typical of
    enthalpy-driven association). Tune to your system.
    """
    K2 = np.exp(-(dH2 - T * dS2) / (R_KJ * T))
    K3 = np.exp(-(dH3 - T * dS3) / (R_KJ * T))
    return K2, K3


def solve_populations(C_tot, K2, K3):
    """Return number concentrations (c1, c2, c3) given total and constants."""
    def mass_balance(c1):
        return c1 + 2 * K2 * c1**2 + 3 * K3 * c1**3 - C_tot

    # c1 is between 0 and C_tot
    c1 = brentq(mass_balance, 0.0, C_tot, xtol=1e-14, rtol=1e-12)
    c2 = K2 * c1**2
    c3 = K3 * c1**3
    return c1, c2, c3


def gold_fractions(C_tot, K2, K3):
    """Fraction of total gold residing in monomer / dimer / trimer."""
    c1, c2, c3 = solve_populations(C_tot, K2, K3)
    f1 = c1 / C_tot
    f2 = 2 * c2 / C_tot
    f3 = 3 * c3 / C_tot
    return dict(monomer=f1, dimer=f2, trimer=f3, number=(c1, c2, c3))
