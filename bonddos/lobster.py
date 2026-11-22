"""Writers for LOBSTER-style text files (readable by pymatgen / LobsterPy).

The numbers are CP2K Gaussian-AO (or Löwdin) quantities, not LOBSTER projections; the first line of
every file says so. Energies are relative to the chosen zero (E_F), like LOBSTER.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .constants import ELEMENTS, m_name_lobster, nominal_principal_n

NOTE = "written by CP2K-BondDOS from CP2K Gaussian-AO data (not a LOBSTER projection)"


def _atom_label(st, i):
    return f"{st.symbol(i)}{i}"


def write_coxxcar(path, kind, energies, efermi_abs, bonds, curves, cumul, spins, structure):
    """kind: COHP|COOP|COBI. curves[bi][s] / cumul[bi][s] arrays on `energies`.

    Column layout (LOBSTER): E, avg, iavg, (X_i, IX_i)... for spin 1, then the same for spin 2.
    """
    nb = len(bonds)
    lines = [f"{kind}CAR.lobster {NOTE}",
             f"{nb + 1:7d}{len(spins):8d}{len(energies):8d}{energies[0]:13.5f}{energies[-1]:13.5f}{efermi_abs:13.5f}",
             "Average"]
    for i, b in enumerate(bonds, 1):
        lines.append(f"No.{i}:{_atom_label(structure, b.a)}->{_atom_label(structure, b.b)}({b.length:.5f})")
    cols = [energies]
    for s in spins:
        avg = np.mean([curves[i][s] for i in range(nb)], axis=0) if nb else np.zeros_like(energies)
        iavg = np.mean([cumul[i][s] for i in range(nb)], axis=0) if nb else np.zeros_like(energies)
        cols += [avg, iavg]
        for i in range(nb):
            cols += [curves[i][s], cumul[i][s]]
    data = np.column_stack(cols)
    body = "\n".join(" ".join(f"{x:.5f}" for x in row) for row in data)
    Path(path).write_text("\n".join(lines) + "\n" + body + "\n", encoding="utf-8")


def write_icoxxlist(path, kind, bonds, integrals, spins, structure, nspin_total):
    """integrals[i][s] per spin channel (restricted: total incl. spin degeneracy)."""
    head = {"COHP": "ICOHP", "COOP": "ICOOP", "COBI": "ICOBI"}[kind]
    lines = []
    for s in spins:
        sp = s + 1
        lines.append(f"  {kind}#  atomMU  atomNU     distance      translation      {head} (eF) for spin  {sp}")
        for i, b in enumerate(bonds, 1):
            t = b.r_raw
            lines.append(f"{i:7d}{_atom_label(structure, b.a):>8s}{_atom_label(structure, b.b):>8s}"
                         f"{b.length:13.5f}{t[0]:5d}{t[1]:5d}{t[2]:5d}{integrals[i - 1][s]:13.5f}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_charge(path, structure, mulliken, loewdin):
    lines = ["", "#      Atom    Mulliken   Loewdin", "-" * 40]
    for i in range(1, len(structure.atoms) + 1):
        lw = loewdin[i - 1] if loewdin is not None else float("nan")
        lines.append(f"{i:6d}{structure.symbol(i):>6s}{mulliken[i - 1]:12.2f}{lw:12.2f}")
    lines += ["-" * 40, "", ""]
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def _lobster_orbital_columns(system, atom):
    """Orbital columns for one atom: shells of the same (l, m) are summed (LOBSTER-like minimal basis)."""
    d = system.d
    a = d.atoms[atom - 1]
    idx = system.atom_aos(atom)
    cols = {}
    order = []
    for mu in idx:
        l, m = int(d.ao_l[mu]), int(d.ao_m[mu])
        n = nominal_principal_n(a.z_nuc, l)
        name = f"{n}{m_name_lobster(l, m)}"
        if name not in cols:
            cols[name] = []
            order.append(name)
        cols[name].append(mu)
    return [(n, np.array(cols[n])) for n in order]


def write_grosspop(path, system, pop_m, pop_l):
    """GROSSPOP-like table: per orbital Mulliken and Löwdin gross populations."""
    lines = [f"# GROSSPOP {NOTE}",
             ("#    Atom  Orbital   Mulliken GP(up)  Mulliken GP(down)  Loewdin GP(up)  Loewdin GP(down)" if system.ns == 2
              else "#    Atom  Orbital   Mulliken GP   Loewdin GP")]
    for atom in range(1, len(system.d.atoms) + 1):
        el = system.d.atoms[atom - 1].element
        first = True
        tm = np.zeros(system.ns)
        tl = np.zeros(system.ns)
        for name, idx in _lobster_orbital_columns(system, atom):
            vm = pop_m.ao_pop[:, idx].sum(axis=1)
            vl = pop_l.ao_pop[:, idx].sum(axis=1) if pop_l is not None else np.full(system.ns, np.nan)
            tm += vm
            tl += vl
            pre = f"{atom:5d}{el:>6s}" if first else " " * 11
            first = False
            lines.append(pre + f"  {name:<12s}" + "".join(f"{x:12.4f}" for x in vm) + "".join(f"{x:12.4f}" for x in vl))
        lines.append(" " * 11 + f"  {'total':<12s}" + "".join(f"{x:12.4f}" for x in tm) + "".join(f"{x:12.4f}" for x in tl))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_doscar(path, system, grid, tdos, pdos_cols, spins):
    """DOSCAR.lobster-like file. tdos[s] curve; pdos_cols[atom] = list of (name, {s: curve})."""
    nat = len(system.d.atoms)
    de = grid[1] - grid[0] if grid.size > 1 else 1.0
    head = f"{grid[-1]:15.8f}{grid[0]:15.8f}{grid.size:6d}{0.0:15.8f}{1.0:15.8f}"
    lines = [f"{nat:4d}{nat:4d}{1:4d}{0:4d}", "  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.1000000E-15",
             "  1.0000000000000000E-004", "  CAR ", f" DOSCAR.lobster {NOTE}", head]
    integ = {s: np.cumsum(tdos[s]) * de for s in spins}
    for i in range(grid.size):
        vals = [grid[i]] + [tdos[s][i] for s in spins] + [integ[s][i] for s in spins]
        lines.append(" ".join(f"{x:.8f}" for x in vals))
    for atom in range(1, nat + 1):
        cols = pdos_cols[atom]
        z = system.d.atoms[atom - 1].z_nuc
        lines.append(head + f"; Z= {z}; " + " ".join(n for n, _ in cols))
        for i in range(grid.size):
            vals = [grid[i]]
            for _, cur in cols:
                for s in spins:
                    vals.append(cur[s][i])
            lines.append(" ".join(f"{x:.8f}" for x in vals))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_poscar(path, structure, title="CP2K structure (raw coordinates, consistent with ICOXXLIST translations)"):
    lines = [title, "1.0"]
    for v in structure.cell:
        lines.append("".join(f"{x:16.10f}" for x in v))
    syms = [a.element for a in structure.atoms]
    groups = []
    for s in syms:
        if groups and groups[-1][0] == s:
            groups[-1][1] += 1
        else:
            groups.append([s, 1])
    lines.append(" ".join(g[0] for g in groups))
    lines.append(" ".join(str(g[1]) for g in groups))
    lines.append("Direct")
    for f in structure.frac_raw:
        lines.append("".join(f"{x:16.10f}" for x in f))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def element_z(sym):
    try:
        return ELEMENTS.index(sym)
    except ValueError:
        return 0
