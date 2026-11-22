"""Physical constants and orbital naming helpers."""

# CODATA 2018/2022 (identical to the digits used here)
HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903
FWHM_PER_SIGMA = 2.3548200450309493  # 2*sqrt(2 ln 2)

L_LETTER = "spdfghi"

# CP2K spherical order m = -l..l (see kpoint_mo_dump.F / sgf_symbol)
_M_NAMES = {
    0: {0: "s"},
    1: {-1: "py", 0: "pz", 1: "px"},
    2: {-2: "dxy", -1: "dyz", 0: "dz2", 1: "dxz", 2: "dx2-y2"},
    3: {-3: "fy(3x2-y2)", -2: "fxyz", -1: "fyz2", 0: "fz3", 1: "fxz2", 2: "fz(x2-y2)", 3: "fx(x2-3y2)"},
}

# LOBSTER style orbital names (without principal quantum number)
_M_LOBSTER = {
    0: {0: "s"},
    1: {-1: "p_y", 0: "p_z", 1: "p_x"},
    2: {-2: "d_xy", -1: "d_yz", 0: "d_z^2", 1: "d_xz", 2: "d_x^2-y^2"},
    3: {-3: "f_y(3x^2-y^2)", -2: "f_xyz", -1: "f_yz^2", 0: "f_z^3", 1: "f_xz^2", 2: "f_z(x^2-y^2)", 3: "f_x(x^2-3y^2)"},
}

# CP2K raw PDOS component labels -> friendly names
CP2K_COMPONENT_ALIASES = {
    "s": "s",
    "py": "py", "pz": "pz", "px": "px",
    "d-2": "dxy", "d-1": "dyz", "d0": "dz2", "d+1": "dxz", "d+2": "dx2-y2",
    "f-3": "fy(3x2-y2)", "f-2": "fxyz", "f-1": "fyz2", "f0": "fz3", "f+1": "fxz2", "f+2": "fz(x2-y2)", "f+3": "fx(x2-3y2)",
}


def m_name(l, m):
    try:
        return _M_NAMES[l][m]
    except KeyError:
        return f"{L_LETTER[l] if l < len(L_LETTER) else 'l' + str(l)}{m:+d}"


def m_name_lobster(l, m):
    try:
        return _M_LOBSTER[l][m]
    except KeyError:
        return f"{L_LETTER[l]}_{m:+d}"


def l_name(l):
    return L_LETTER[l] if l < len(L_LETTER) else f"l{l}"


ELEMENTS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
]


def nominal_principal_n(z, l):
    """Valence principal quantum number used only for LOBSTER-style labels (e.g. '5d')."""
    period_edges = [2, 10, 18, 36, 54, 86, 118]
    period = next(i + 1 for i, e in enumerate(period_edges) if z <= e) if z > 0 else 1
    if l == 0:
        return period
    if l == 1:
        return max(2, period)
    if l == 2:
        return max(3, period - 1)
    if l == 3:
        return max(4, period - 2)
    return l + 1
