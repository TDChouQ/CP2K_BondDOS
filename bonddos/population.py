"""Mulliken / Löwdin populations, charges and atom/orbital-resolved PDOS from MO_KP data."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .constants import l_name, m_name
from .kspace import KSystem
from .spectra import broaden


@dataclass(eq=False)
class PopulationResult:
    method: str
    # per AO gross population per spin (nspin, nao) and total
    ao_pop: np.ndarray
    atom_pop: np.ndarray            # (nspin, natom)
    charges: np.ndarray             # (natom,) Z_eff - total population
    moments: np.ndarray | None      # (natom,) up - down for UKS
    labels: list = field(default_factory=list)


def state_weights(system: KSystem, k, s, method="mulliken"):
    """(nao, nmo) projection weights of each state on each AO (columns sum to 1)."""
    C = system.C(k, s)
    if method == "lowdin":
        sh, _, _, _ = system.lowdin(k)
        Q = sh @ C
        return np.abs(Q) ** 2
    SC = system.S(k) @ C
    return np.real(C.conj() * SC)


def populations(system: KSystem, method="mulliken", progress=None) -> PopulationResult:
    d = system.d
    nao = system.nao
    ao_pop = np.zeros((system.ns, nao))
    for k in range(system.nk):
        if progress:
            progress(k / system.nk, f"布居分析 k {k + 1}/{system.nk}")
        for s in range(system.ns):
            W = state_weights(system, k, s, method)
            ao_pop[s] += system.wts[k] * (W @ system.occ(k, s))
    natom = len(d.atoms)
    atom_pop = np.zeros((system.ns, natom))
    for i, a in enumerate(d.atoms):
        atom_pop[:, i] = ao_pop[:, a.first_ao - 1:a.last_ao].sum(axis=1)
    zeff = np.array([a.z_eff for a in d.atoms])
    charges = zeff - atom_pop.sum(axis=0)
    moments = atom_pop[0] - atom_pop[1] if system.ns == 2 else None
    return PopulationResult(method, ao_pop, atom_pop, charges, moments, [a.label for a in d.atoms])


def gross_populations(system: KSystem, pop: PopulationResult):
    """Rows (atom, element, shell, l, orbital, pop_up, pop_down/None, total)."""
    d = system.d
    rows = []
    for a in d.atoms:
        for mu in range(a.first_ao - 1, a.last_ao):
            l, m, sh = int(d.ao_l[mu]), int(d.ao_m[mu]), int(d.ao_shell[mu])
            up = float(pop.ao_pop[0, mu])
            dn = float(pop.ao_pop[1, mu]) if system.ns == 2 else None
            rows.append((a.index, a.element, sh, l, m_name(l, m), up, dn, up + (dn or 0.0)))
    return rows


@dataclass
class PDOSGroup:
    name: str
    atoms: list                # 1-based
    resolution: str = "total"  # total | l | lm | shell


def group_ao_channels(system: KSystem, group: PDOSGroup):
    """[(channel_name, ao_indices)] for a group."""
    d = system.d
    idx = np.concatenate([system.atom_aos(a) for a in group.atoms]) if group.atoms else np.array([], int)
    res = group.resolution
    if res == "total":
        return [("total", idx)]
    chans = {}
    order = []
    for mu in idx:
        l, m, sh = int(d.ao_l[mu]), int(d.ao_m[mu]), int(d.ao_shell[mu])
        if res == "l":
            key = l_name(l)
        elif res == "lm":
            key = m_name(l, m)
        else:
            key = f"{l_name(l)}({sh})"
        if key not in chans:
            chans[key] = []
            order.append(key)
        chans[key].append(mu)
    out = [("total", idx)]
    out += [(k, np.array(chans[k], int)) for k in order]
    return out


def kpoint_pdos(system: KSystem, groups: list[PDOSGroup], eref, grid, fwhm, shape="gaussian",
                method="mulliken", include_tdos=True, progress=None, occupied=False):
    """Atom/orbital-resolved PDOS curves.

    Returns {spin: {series_name: curve}} with states/eV per spin channel (restricted: both spins included).
    """
    plan = []
    for g in groups:
        for ch, idx in group_ao_channels(system, g):
            plan.append((f"{g.name}" if ch == "total" else f"{g.name}·{ch}", idx))
    out = {}
    for s in range(system.ns):
        energies = []
        weights = []
        for k in range(system.nk):
            if progress:
                progress((s * system.nk + k) / (system.ns * system.nk), f"PDOS k {k + 1}/{system.nk}")
            W = state_weights(system, k, s, method) if plan else None
            e = system.eps_ev(k, s) - eref
            fac = system.g * system.wts[k] * (system.occ(k, s) / system.g if occupied else 1.0)
            cols = []
            if include_tdos:
                cols.append(np.ones(system.nmo) * fac)
            for _, idx in plan:
                cols.append(W[idx].sum(axis=0) * fac)
            energies.append(e)
            weights.append(np.column_stack(cols) if cols else np.zeros((system.nmo, 0)))
        E = np.concatenate(energies)
        Wt = np.vstack(weights)
        curves = broaden(E, Wt, grid, fwhm, shape)
        names = (["TDOS"] if include_tdos else []) + [p[0] for p in plan]
        out[s] = {n: curves[:, i] for i, n in enumerate(names)}
    return out


def stick_weights(system: KSystem, group: PDOSGroup, method="mulliken", channel="total"):
    """Discrete (energy_abs_eV, weight) arrays for descriptor statistics of one channel."""
    idx = dict(group_ao_channels(system, group)).get(channel)
    if idx is None:
        raise KeyError(channel)
    es, ws = [], []
    for s in range(system.ns):
        for k in range(system.nk):
            W = state_weights(system, k, s, method)
            es.append(system.eps_ev(k, s))
            ws.append(system.g * system.wts[k] * W[idx].sum(axis=0))
    return np.concatenate(es), np.concatenate(ws)
