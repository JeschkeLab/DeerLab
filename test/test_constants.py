import numpy as np
import deerlab.constants as con


def test_planck():
    h = 6.62607015e-34       # Planck constant, J/Hz (CODATA 2018)
    assert con.h==h
    
def test_reducedplanck():
    hbar = con.h/2/con.π
    assert con.hbar==hbar
    
def test_ge():
    ge = 2.00231930436256    # free-electron g factor (CODATA 2018 value)
    assert con.ge==ge

def test_avogadro():
    Nav = 6.02214076e23      # Avogadro constant, mol^-1
    assert con.Nav==Nav

def test_bohrmagneton():
    μB = 9.2740100783e-24    # Bohr magneton, J/T (CODATA 2018 value)
    assert con.μB==μB
    
def test_magneticonstant():
    μ0 = 1.25663706212e-6    # magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
    assert con.μ0==μ0
