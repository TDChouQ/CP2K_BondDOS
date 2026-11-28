"""CP2K input fragments that CP2K BondDOS needs (verified with CP2K 2026.2).

``python main.py template [--mesh 3 3 3] [--uks] [--smear] [--added-mos 50] [--no-dos]`` prints them.
Comments are on their own lines on purpose: CP2K takes the rest of a FILENAME line literally
(an inline "! comment" there ends up in the file name).  The same text is kept in ``examples/cp2k_input/bonddos_print.inc`` (PRINT_BLOCK only).
"""
from __future__ import annotations

# goes inside  &FORCE_EVAL / &DFT / &PRINT   (file can be used with:  @INCLUDE 'bonddos_print.inc')
PRINT_BLOCK = """\
      ! ---- CP2K BondDOS: k-point MO coefficients, eigenvalues, occupations and S(k)  ->  *.mokp
      &MO_KP ON
        ! -> cohp_mo-1_0.mokp
        FILENAME ./cohp_mo
        ! REQUIRED: the default GTO_BASIS does not write S(k)
        AO_EXPORT_TYPE OVERLAP_MATRIX
        ! default 9 also works; 12-14 gives tighter consistency checks
        NDIGITS 12
        UNIT ANGSTROM
      &END MO_KP
      ! ---- CP2K BondDOS: Kohn-Sham matrix H(k), one text file per spin and k-point  ->  *.csr
      &KS_CSR_WRITE ON
        ! -> cohp_ks-KS_SPIN_<s>_K_<k>-1_0.csr
        FILENAME ./cohp_ks
        ! H(k), not H(R)          (default)
        REAL_SPACE F
        ! text files              (default)
        BINARY F
        ! T is also read (detected automatically)
        UPPER_TRIANGULAR F
      &END KS_CSR_WRITE
"""

DOS_BLOCK = """\
      ! ---- optional: CP2K's own DOS / PDOS per KIND (first page of BondDOS)
      &DOS
        ! Hartree (= 0.01 eV)
        DELTA_E 0.0003674932
        NLUMO {nlumo}
        &CURVE
          ENERGY_UNIT EV
          ENERGY_ZERO FERMI
          &BROADEN
            TYPE GAUSSIAN
            ! FWHM
            WIDTH [eV] 0.10
          &END BROADEN
        &END CURVE
        &PDOS
          ! T: px/py/pz, dxy ... resolved
          COMPONENTS F
        &END PDOS
      &END DOS
"""


def cp2k_template(mesh=(3, 3, 3), added_mos=50, uks=False, smear=False, dos=True) -> str:
    """Commented CP2K input fragment (everything that has to be added to a normal single-point input)."""
    gamma_only = tuple(mesh) == (1, 1, 1)
    even = any(n % 2 == 0 for n in mesh)
    o = ["! ====================================================================",
         "!  CP2K BondDOS: what to add to a CP2K 2026.2+ single-point input",
         "!  (RUN_TYPE ENERGY on the final geometry; three places, all inside &DFT)",
         "! ====================================================================",
         "  &DFT"]
    if uks:
        o += ["    ! spin polarised: BondDOS gives up/down/total and moments", "    UKS T",
              "    ! MULTIPLICITY ... and MAGNETIZATION in the &KIND sections as usual"]
    o += ["",
          "    ! (1) k-points. REQUIRED even for a Gamma-only supercell (SCHEME GAMMA):",
          "    !     without &KPOINTS no .mokp file is written at all.",
          "    &KPOINTS"]
    if gamma_only:
        o += ["      ! only PAIR mode (all images summed) is meaningful", "      SCHEME GAMMA"]
    else:
        o += [f"      SCHEME MONKHORST-PACK {mesh[0]} {mesh[1]} {mesh[2]}"]
        if even:
            o += ["      ! even meshes: keeps the lattice symmetry (fcc, hexagonal ...)", "      GAMMA_CENTERED T"]
    o += ["      ! default; do NOT use point-group reduction", "      SYMMETRY F",
          "      ! FULL_GRID T                   ! optional: BondDOS unfolds the default time-reversal reduced mesh",
          "      ! WAVEFUNCTIONS COMPLEX         ! default",
          "    &END KPOINTS",
          "",
          "    ! (2) SCF: diagonalisation (OT cannot be used with k-points) and enough empty bands",
          "    &SCF",
          f"      ! MO_KP exports only these bands: they must cover the antibonding region you want to plot (-1 = all)", f"      ADDED_MOS {added_mos:<6d}",
          "      &DIAGONALIZATION",
          "      &END DIAGONALIZATION",
          "      &MIXING",
          "        METHOD BROYDEN_MIXING",
          "        ALPHA 0.2",
          "      &END MIXING"]
    if smear:
        o += ["      ! metals / small-gap magnetic oxides", "      &SMEAR ON",
              "        METHOD FERMI_DIRAC",
              "        ELECTRONIC_TEMPERATURE [K] 300",
              "      &END SMEAR"]
    o += ["    &END SCF",
          "",
          "    ! (3) output for BondDOS: put these sections into the existing &DFT / &PRINT",
          "    &PRINT",
          PRINT_BLOCK.rstrip("\n")]
    if dos:
        o += [DOS_BLOCK.format(nlumo=added_mos).rstrip("\n")]
    o += ["    &END PRINT",
          "  &END DFT",
          "",
          "! After the run give BondDOS:  cohp_mo-1_0.mokp  +  all cohp_ks-KS_SPIN_*_K_*-1_0.csr",
          "!   (same folder; optional: *.pdos / *.dos for the DOS page)"]
    return "\n".join(o) + "\n"
