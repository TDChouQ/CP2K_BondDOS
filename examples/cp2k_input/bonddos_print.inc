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
