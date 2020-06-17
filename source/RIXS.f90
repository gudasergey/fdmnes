! FDMNES subroutine

! Set of subroutines for RIXS

!***********************************************************************

! nlm_f = nlmomax: number of harmonics in the outer sphere (= number of states at the current energy)

subroutine main_RIXS(Ampl_abs,Coef_g,Core_resolved,dV0bdcF,E_cut,E_cut_imp,E_Fermi,E_cut_man,Eclie,Eneg,Energ_rixs, &
        Energ_s,Epsii,Epsii_ref,Epsii_ref_man,Eseuil,Estart,File_name_rixs,FDM_comp_m,Full_potential,Gamma_hole, &
        Gamma_hole_man,Hubb_a,Hubb_d,iabsorig, &
        icheck,igreq,iprabs_nonexc,isymeq,ip0,ip_max,is_g,jseuil,lmoins1,lplus1,lseuil,m_g,m_hubb, &
        Moment_conservation,multi_run, &
        Multipole,n_atom_proto,n_atom_uc,n_multi_run,n_oo,n_q_rixs,n_rel,n_theta_in,n_two_theta,natomsym,nbseuil, &
        nenerg_in_rixs,nenerg_s,neqm,ngamh,ninit1,ninit,ninitr,nlm_pot,nlm_probe,nlm_f,nomfich, &
        nr,nrm,ns_dipmag,nseuil,nspin,nspino,nspinp,numat,Orthmati,Posn,Powder, &
        psii,q_rixs,r,Relativiste,Renorm,RIXS_fast,Rmtg,Rmtsd,Rot_atom_abs,Rot_int,Spin_channel,Spinorbite, &
        Step_loss,Theta_in,Two_theta,V_intmax,V_hubb,Vcato,Vhbdc,VxcbdcF,Vxcato,Ylm_comp)

  use declarations
  implicit none

  integer, parameter:: n_Ec = 2
  integer, parameter:: ndim2 = 1
  integer, parameter:: nenerg_tddft = 0
  integer, parameter:: nlmamax = 0

  real(kind=db), parameter:: quatre_mc2 = 4 * 2 / alfa_sf**2  ! in Rydberg
  
  integer:: iabsorig, i_q, i_Spin_channel, i_theta, i_theta_1, i_theta_2, ia, ich, icheck, ie, &
     ie_fast, ie_loss, ie_not_fast, ie_oc, ie_s, ii, ip0, ip_max, initr, iprabs_nonexc, iseuil, isp, isym, &
     js, jseuil, Length, lmax_pot, lmax, lseuil, m_hubb, multi_run, n_atom_proto, n_atom_uc, &
     n_isc, n_multi_run, n_oo, n_q_dim, n_q_rixs, n_rel, n_Spin_channel, n_theta_in, n_two_theta, n_theta, natomsym, nbseuil, &
     ne, ne_loss, nenerg_in, nenerg_in_rixs, nenerg_oc, nenerg_unoc, nenerg_s, neqm, ngamh, ninit, ninit1, ninitr, nlm_f, &
     nlm_p_fp, nlm_pot, nlm_probe, npl, nr, nrm, &
     ns_dipmag, nseuil, nspin, nspino, nspinp, numat

  integer, dimension(natomsym):: isymeq
  integer, dimension(ninit,2):: m_g
  integer, dimension(ninit):: is_g
  integer, dimension(0:n_atom_proto,neqm):: igreq

  character(len=9):: nomsym
  character(len=Length_name):: File_name_rixs, nomfich

  complex(kind=db):: CFac, Lorentzian
  complex(kind=db), dimension(nenerg_s,nlm_probe,nspinp,nlm_f,nspinp):: Ampl_abs
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull, Taull_Spin_channel
  complex(kind=db), dimension(3,3,2,ninitr):: secmd
  complex(kind=db), dimension(3,3,ninitr):: secmm
  complex(kind=db), dimension(3,3,3,2,ninitr):: secdq
  complex(kind=db), dimension(3,3,3,3,2,ninitr):: secdo
  complex(kind=db), dimension(3,3,3,3,ninitr):: secqq
  complex(kind=db), dimension(3,3,n_rel,ninitr):: secdd
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr):: secoo
  complex(kind=db), dimension(3,3):: RIXS_M1M1_tens_xtal
  complex(kind=db), dimension(3,3,2):: RIXS_E1M1_tens_xtal
  complex(kind=db), dimension(3,3,n_rel):: RIXS_E1E1_tens_xtal
  complex(kind=db), dimension(3,3,3,2):: RIXS_E1E2_tens_xtal
  complex(kind=db), dimension(3,3,3,3):: RIXS_E2E2_tens_xtal, RIXS_E3E3_tens_xtal
  complex(kind=db), dimension(3,3,3,3,2):: RIXS_E1E3_tens_xtal

  complex(kind=db), dimension(:), allocatable:: Amplitude
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: RIXS_E1M1_tens, RIXS_M1M1_tens
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: RIXS_ampl, RIXS_E1E1_tens, RIXS_E1E2_tens, rof, rof_e
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: RIXS_E1E3_tens, RIXS_E2E2_tens, RIXS_E3E3_tens

  complex(kind=db), dimension(:,:,:,:), allocatable:: secmm_xtal
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: secdd_xtal, secmd_xtal
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: secdq_xtal, secqq_xtal, secoo_xtal
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: secdo_xtal

  logical:: Circular, Core_resolved, Dip_rel, E_cut_man, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, &
            Epsii_ref_man, Extract_ten, Eneg, FDM_comp_m, Final_optic, Final_tddft, Full_potential, Gamma_hole_man, &
            Hubb_a, Hubb_d, k_not_possible, lmoins1, lplus1, M1M1, Magn_sens, Moment_conservation, &
            Monocrystal, Only_Intensity, &
            Powder, Relativiste, Renorm, RIXS, RIXS_fast, Run_fast, Spin_channel, Spinorbite, Tddft, Time_reversal, Write_bav, &
            Ylm_comp 

  Logical, save:: Epsii_ref_RIXS_man

  logical, dimension(10):: Multipole
  logical, dimension(:), allocatable:: Probed_L

  logical, dimension(natomsym):: Atom_done
  
  real(kind=db):: conv_mbarn_nelec, cst, Delta_E_inf, Delta_E_sup, Deltar, E_cut, E_cut_imp, E_cut_o, E_cut_rixs, E_Fermi, E_i, &
            E_max, E_s, Ecinetic_i, Ecinetic_s, Eclie, Eimag, Energ_unoc, Energ_unoc_1, Ephoton, Epsii_ref, Estart, &
            Gamma_hole_o, Gamma_max, k_elec_i, k_elec_s, Lorentzian_i, Lorentzian_r, Mod_Q, Omega_i, Omega_s, Pas, pds, &
            Rmtg, Rmtsd, Rot_sample, Step_loss, Tab_width, Theta, Two_theta_L, V_intmax, V0, Vhbdc
  real(kind=db), dimension(3):: Pol_i_s, Pol_i_p, Pol_s_s, Pol_s_p, Q_vec, Vec_i, Vec_s
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitr):: Epsii
  real(kind=db), dimension(ninit,2):: coef_g
  real(kind=db), dimension(nspin):: dV0bdcF, V0bdc, VxcbdcF
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(n_theta_in):: Theta_in
  real(kind=db), dimension(n_two_theta):: Two_theta
  real(kind=db), dimension(nenerg_in_rixs):: Energ_rixs
  real(kind=db), dimension(3,3):: Mat_orth, matopsym , Orthmat, Orthmati, Rot_atom, Rot_atom_abs, Rot_int, &
                                  Trans_Lab_Dir_Q, Trans_rec_dir, Trans_rec_lab
  real(kind=db), dimension(3,n_atom_uc):: posn
  real(kind=db), dimension(4,n_q_rixs):: q_rixs
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nlm_pot):: Vcato
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato, Vxcato

  real(kind=db), dimension(:), allocatable:: E_inf, E_loss, E_sup, Energ_in, Energ_oc, Gamma, Shift
  real(kind=db), dimension(:,:), allocatable:: Angle, Conv_nelec
  real(kind=db), dimension(:,:,:,:), allocatable:: Mod_k

  if( icheck > 0 ) write(3,110)
  
! When no momentum transfer direction given, calculation is done for powder 
  Monocrystal = .not. Powder

  Eimag = 0._db
  Extract_ten = .false.
  Final_optic = .false.
  Final_tddft = .false.
  RIXS = .true.
  Tddft = .false.
  
 ! With Run_fast, Tau_loss does not depend on incoming energy, and consequently the tensors
  Run_fast = Powder .or. RIXS_fast .or. .not. Moment_conservation 
  if( multi_run == 1 ) Epsii_ref_RIXS_man = Epsii_ref_man

  ich = max(icheck-1,1)
  Magn_sens = nspin == 2 .or. Spinorbite
  Circular = Magn_sens .and. Monocrystal

  Only_Intensity = .false.

! output file name
  File_name_rixs = nomfich
  Length = len_trim(File_name_rixs)

  File_name_rixs(Length+1:Length+5) = '_rixs'

  if( n_multi_run > 1 ) then
    Length = Length + 6
    File_name_rixs(Length:Length) = '_'
    call ad_number(iabsorig,File_name_rixs,Length_name)
  endif
  Length = len_trim(File_name_rixs)
  File_name_rixs(Length+1:Length+4) = '.txt'
  
  if( Only_Intensity ) return
  
  Dip_rel = n_rel > 1
  
  if( Circular ) then
    npl = 8
  elseif( Powder ) then
    npl = 6
  else
    npl = 4
  endif
  allocate( Amplitude(npl) )

  if( Hubb_a .or. Full_potential ) then
    nlm_p_fp = nlm_probe
  else
    nlm_p_fp = 1
  endif
  allocate( rof_e(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,nenerg_s)  )  
  allocate( rof(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,n_Ec)  )  

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  if( Spin_channel ) then
    n_Spin_channel = 5
  else
    n_Spin_channel = 1
  endif

  allocate( secdd_xtal(3,3,n_rel,ninitr,n_spin_channel) )
  allocate( secmd_xtal(3,3,2,ninitr,n_spin_channel) )
  allocate( secmm_xtal(3,3,ninitr,n_spin_channel) )
  allocate( secdq_xtal(3,3,3,2,ninitr,n_spin_channel) )
  allocate( secdo_xtal(3,3,3,3,2,ninitr,n_spin_channel) )
  allocate( secqq_xtal(3,3,3,3,ninitr,n_spin_channel) )
  allocate( secoo_xtal(3,n_oo,3,n_oo,ninitr,n_spin_channel) )
  
  n_isc = ninitr * n_Spin_channel

  allocate( Shift(ninitr) )
  call Cal_shift(Epsii,Epsii_ref,Epsii_ref_RIXS_man,icheck,ninit,ninit1,ninitr,Shift)

  lmax = nint( sqrt( real( nlm_probe ) ) ) - 1
  lmax_pot = nint( sqrt( real( nlm_pot ) ) ) - 1

  allocate( Probed_L(0:Lmax) )
  Probed_L(:) = .false.
  if( E1M1 .or. M1M1 .or. ( ( E2E2 .or. E1E2 ) .and. Lseuil == 2 ) ) Probed_L(Lseuil) = .true.
  if( E1E1 .or. E1M1 .or. E1E2 .or. E1E3 .or. E3E3 ) Probed_L(Lseuil+1) = .true.
  if( ( E1E1 .or. E1E3 .or. E3E3 ) .and. ( Lseuil == 1 .or. Lseuil == 3 ) ) Probed_L(Lseuil-1) = .true. 
  if( E1E2 .or. E2E2 ) Probed_L(Lseuil+2) = .true.
  if( ( E1E2 .or. E2E2 ) .and. Lseuil == 2 ) Probed_L(Lseuil-2) = .true.
  if( E1E3 .or. E3E3 ) Probed_L(Lseuil+3) = .true.

  if( E_cut_man ) then
    E_cut_rixs = E_cut_imp
    E_cut_o = E_cut_imp
  else
    E_cut_rixs = E_cut
    E_cut_o = - 1001._db
  endif

  if( Gamma_hole_man ) then
    Gamma_hole_o = Gamma_hole(1) * Rydb
  else
    Gamma_hole_o = -1001._db
  endif
  Gamma_max = -1001._db
  Deltar = 0._db
  
  do ie = 1,nenerg_s
    if( Energ_s(ie) > E_cut_RIXS + eps10 ) exit
  end do
  
  if( ie <= nenerg_s ) then
    Energ_unoc_1 = Energ_s(ie)
  else
    Energ_unoc_1 = Energ_s(nenerg_s)
  endif 

! E_cut rixs is not let equal to a grid point energy. 
  if( ie > 1 ) then
    if( abs( Energ_s(ie-1) - E_cut_RIXS ) < eps10 ) E_cut_RIXS = ( Energ_s(ie-1) + Energ_s(ie) ) / 2 
  endif
  
  nenerg_oc = ie - 1
  nenerg_unoc = nenerg_s - nenerg_oc 

  if( nenerg_in_rixs > 0 ) then
    nenerg_in = nenerg_in_rixs
  else 
    Estart = max( Energ_s(1) + eps10, Estart )  
    if( Estart < E_cut_RIXS .and. ie > 1 ) then
      Pas =  Energ_s(ie) - Energ_s(ie-1)
      ne = nint( ( Energ_s(ie) - Estart ) / pas )
    else
      ne = 0
    endif
    nenerg_in = nenerg_unoc + ne
  endif

  if( nenerg_in_rixs > 0 ) then
    E_max = Energ_rixs(nenerg_in_rixs) - eps10
  else
    E_max = Energ_s(nenerg_s) - eps10  
  endif

  if( Step_loss < 0._db ) Step_loss = Energ_s(2) - Energ_s(1) 

  do ie = nenerg_s, 1, -1
    if( Energ_s(ie) > E_max ) cycle
    ne_loss = nint( ( Energ_s(ie) - Energ_s(1) ) / Step_loss )
    exit
  end do

  write(6,'(a30,i5)') ' Number of incoming energies =', nenerg_in
  write(6,'(a43,i5)') ' Number of energies below the Fermi level =', nenerg_oc
  
  allocate( E_inf(nenerg_oc) )
  allocate( E_sup(nenerg_oc) )
  allocate( Energ_oc(nenerg_oc) )
  allocate( Energ_in(nenerg_in) )
  allocate( Gamma(ninitr) )
  allocate( E_loss(ne_loss) )
  allocate( Conv_nelec(nenerg_in,ninitr) )

  if( E1E1 ) then
    allocate( RIXS_E1E1_tens(3,3,n_rel,n_isc,nenerg_in*ne_loss,natomsym) )
    RIXS_E1E1_tens(:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E1E2 ) then
    allocate( RIXS_E1E2_tens(3,3,3,n_isc,nenerg_in*ne_loss,natomsym) )
    RIXS_E1E2_tens(:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E1E3 ) then
    allocate( RIXS_E1E3_tens(3,n_oo,3,n_oo,n_isc,nenerg_in*ne_loss,natomsym) )
    RIXS_E1E3_tens(:,:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E1M1 ) then
    allocate( RIXS_E1M1_tens(3,3,n_isc,nenerg_in*ne_loss,natomsym) )
    RIXS_E1M1_tens(:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E2E2 ) then
    allocate( RIXS_E2E2_tens(3,3,3,3,n_isc,nenerg_in*ne_loss,natomsym) )
    RIXS_E2E2_tens(:,:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E3E3 ) then
    allocate( RIXS_E3E3_tens(3,n_oo,3,n_oo,n_isc,nenerg_in*ne_loss,natomsym) )
    RIXS_E1E3_tens(:,:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( M1M1 ) then
    allocate( RIXS_M1M1_tens(3,3,n_isc,nenerg_in*ne_loss,natomsym) )
    RIXS_M1M1_tens(:,:,:,:,:) = (0._db, 0._db)
  endif

  n_theta = max( n_two_theta, 1 ) *  max( n_theta_in, 1 ) 
  n_q_dim = max( n_q_rixs, 1 )
  allocate( Angle(2,n_theta) )
  allocate( Mod_k(nenerg_in,ne_loss,n_theta,n_q_dim) )
  
  allocate( RIXS_ampl(ne_loss,nenerg_in,npl,n_theta,n_q_dim,n_isc) )
  RIXS_ampl(:,:,:,:,:,:) = (0._db, 0._db) 

  call Energy_Grid_RIXS(E_loss,Energ_in,Energ_oc,Energ_rixs,Energ_s,icheck,ne_loss,nenerg_in,nenerg_in_rixs,nenerg_oc, &
                        nenerg_s,Step_loss)

  do initr = 1,ninitr       ! ----------> Loop over edges or core states
    if( Core_resolved ) then
      if( initr <= ninit1 ) then
        iseuil = 1
      else
        iseuil = min(2, nbseuil)
      endif
    else
       iseuil = min(initr, nbseuil)
    endif
    
! Value if core state width
    if( Gamma_hole_man ) then
      if( ngamh == 1 ) then
        Gamma(initr) = Gamma_hole(1)
      elseif( ngamh == ninitr ) then
        Gamma(initr) = Gamma_hole(initr)
      elseif( initr <= ninit1 ) then
        Gamma(initr) = Gamma_hole(1)
      else
        Gamma(initr) = Gamma_hole(2)
      endif
    else
      js = jseuil + iseuil - 1 
      Gamma_hole(1) = Tab_width(Eseuil(iseuil),js,nbseuil,numat)
      Gamma(initr) = Gamma_hole(1)
    endif
    Gamma(initr) = 0.5_db * Gamma(initr)  ! It is in fact Gamma / 2 in the formula  

! Conversion factor in number of electron
    do ie = 1,nenerg_in
      Ephoton = Eseuil(iseuil) + Energ_in(ie) + Shift(initr)
! For very low energy edges
      Ephoton = max(0.001_db/Rydb, Ephoton)
! alfa_sf = e*e/(2*epsilon0*h*c) is the fine structure constant
      cst = quatre_pi * pi * alfa_sf * Ephoton
! To get result in Mbarn (10E-18 cm2)
      cst = 100 * bohr**2 * cst
      Conv_nelec(ie,initr) = cst * conv_mbarn_nelec(Ephoton) / pi
    end do

  end do

! One takes the exchange correlation potential independant from energy
  do isp = 1,nspin
    Vrato(1:nr,:,isp) = Vcato(1:nr,:) + Vxcato(1:nr,:,isp)
    V0bdc(isp) = Vhbdc + VxcbdcF(isp) + dV0bdcF(isp)
  end do

! Calculation of the radial matrices
  call cal_rof(E_Fermi,Eclie,Energ_s,Eseuil,Full_potential,Hubb_a,Hubb_d,icheck,ip0,ip_max,lmax,lmax_pot, &
                   m_hubb,n_Ec,nbseuil,nenerg_s,ninit1,nlm_pot,nlm_p_fp,nlm_probe,nr,nrm,nspin,nspino,nspinp,numat,psii,r, &
                   Relativiste,Renorm,Rmtg,Rmtsd,rof_e,Spinorbite,V_hubb,V_intmax,V0bdc,Vrato,Ylm_comp)
  
  if( nenerg_oc == 1 ) then
    E_sup(1) = Energ_oc(1) + 0.05_db / Rydb
    E_inf(1) = Energ_oc(1) - 0.05_db / Rydb
  else
    E_inf(1) = 1.5_db * Energ_oc(1) - 0.5_db * Energ_oc(2) 
    do ie_oc = 2,nenerg_oc
      E_inf(ie_oc) = ( Energ_oc(ie_oc-1) + Energ_oc(ie_oc) ) / 2
    end do
    do ie_oc = 1,nenerg_oc-1
      E_sup(ie_oc) = ( Energ_oc(ie_oc) + Energ_oc(ie_oc+1) ) / 2
    end do
    E_sup(nenerg_oc) = E_cut_Rixs 
  endif

! Rotation to go from the local atom basis to the tensors R1 basis
  Rot_atom = matmul( Rot_int, transpose(Rot_atom_abs) )

! Transformation to go from crystal basis to orthonormal basis
  call invermat(Orthmati,Orthmat)
  Mat_orth = matmul( Rot_int, Orthmat )

  call Reciprocal_Vector(icheck,Mat_orth,Trans_rec_dir,Trans_rec_lab)

  do i_q = 1,n_q_dim

    if( n_q_rixs > 0 ) then
      Q_vec(:) = Q_rixs(1:3,i_q)
      Rot_sample = Q_rixs(4,i_q)
      Q_vec = matmul( Trans_rec_dir, Q_vec )
      Q_vec(:) = Q_vec(:) / sqrt( sum( Q_vec(:)**2 ) ) 
    else   ! default 
      Q_vec(1:2) = 0._db; Q_vec(3) = 1._db
      Rot_sample = 0._db
    endif

    call Cal_Mat_rot(icheck,Mat_orth,Q_vec,Rot_sample,Trans_Lab_Dir_Q,Trans_rec_dir,Trans_rec_lab)      

    i_theta = 0
    do i_theta_1 = 1,max( n_theta_in, 1 )
    do i_theta_2 = 1,max( n_two_theta, 1 )
      i_theta = i_theta + 1
      
      if( n_two_theta == 0 .and. n_theta_in == 0) then
        Theta = pi / 4._db
        Two_theta_L = 2 * Theta   
      elseif( n_theta_in == 0 ) then
        Two_theta_L = Two_theta(i_theta_2)
        Theta = Two_theta_L / 2   
      elseif( n_two_theta == 0 ) then
        Theta = Theta_in(i_theta_1)   
        Two_theta_L = 2 * Theta   
      else
        Theta = Theta_in(i_theta_1)
        Two_theta_L = Two_theta(i_theta_2)   
      endif
      Angle(1,i_theta) = Theta
      Angle(2,i_theta) = Two_theta_L
      
      call Cal_Pol(icheck,Pol_i_s,Pol_i_p,Pol_s_s,Pol_s_p,Theta,Trans_Lab_Dir_Q,Two_theta_L,Vec_i,Vec_s)

      do ie_not_fast = 1,nenerg_in  ! loop on incoming energy

        if( Run_fast .and. ie_not_fast > 1 ) exit

        if( .not. Run_fast ) then
          write(6,'(a27,f8.3,a3)') ' Incoming energy - E_edge =', Energ_in(ie_not_fast) * Rydb, ' eV' 
          if( icheck > 0 ) then
            write(3,'(/A)') ' ------------------------------------------------------------------------------------------ '
            write(3,'(/a27,f8.3,a3)') ' Incoming energy - E_edge =', Energ_in(ie_not_fast) * Rydb, ' eV'
          endif 
          E_i = Eseuil(nbseuil) + Energ_in(ie_not_fast)
        endif
      
        do ie_oc = nenerg_oc,1,-1 ! Loop over occupied states

          Loop_loss: do ie_loss = 1,ne_loss
          
            Energ_unoc = Energ_oc(ie_oc) + E_loss(ie_loss)

            if( Energ_unoc < E_cut_RIXS - eps10 ) cycle
            if( Energ_unoc > Energ_s(nenerg_s) + eps10 ) cycle

            if( Run_fast ) E_i = Eseuil(nbseuil) + Energ_unoc

            if( icheck > 1 ) then
              if( Run_fast ) then
                write(3,120) ie_oc, ie_loss, E_i*rydb
              else
                write(3,130) ie_not_fast, ie_oc, ie_loss, E_i*rydb
              endif
            endif 
            write(3,140) Energ_oc(ie_oc)*rydb, Energ_unoc*rydb, E_loss(ie_loss)*rydb
              
            Write_bav = icheck > 1 .or. ( icheck == 1 .and. ie_not_fast == 1 .and. ie_oc == 1 .and. ie_loss == 1 )
                        
            E_s = E_i - E_loss(ie_loss)
            Mod_Q = 0.5_db * alfa_sf * sqrt( E_s**2 + E_i**2 - 2 * E_s * E_i * cos( pi - Two_theta_L ) )
            if( Run_fast ) then
              do ie = 1,nenerg_in
                if( abs( Energ_in(ie) - Energ_unoc ) < eps10 ) then
                  Mod_k(ie,ie_loss,i_theta,i_q) = Mod_Q
                  if( abs( Energ_unoc - Energ_unoc_1 ) < eps10 ) Mod_k(1:ie-1,ie_loss,i_theta,i_q) = Mod_Q
                  if( ie_oc == nenerg_oc ) Mod_k(ie+1:nenerg_in,ie_loss,i_theta,i_q) = Mod_Q 
                  if( ie == nenerg_in .and. ie_oc == nenerg_oc ) Mod_k(1:nenerg_in,ie_loss,i_theta,i_q) = Mod_Q
                elseif( Energ_in(ie) > Energ_unoc ) then
                  exit
                endif
              end do   
            else
              Mod_k(ie_not_fast,ie_loss,i_theta,i_q) = Mod_Q
            endif

            V0 = sum( V0bdc(:) ) / nspin
            Ecinetic_i = Energ_unoc + E_Fermi - V0
            Ecinetic_s = Energ_oc(ie_oc) + E_Fermi - V0
            if( .not. Eneg ) Ecinetic_i = max( Ecinetic_i, Eclie ) 
            if( .not. Eneg ) Ecinetic_s = max( Ecinetic_s, Eclie ) 
            k_elec_i = sqrt( Ecinetic_i )
            k_elec_s = sqrt( Ecinetic_s )

 ! index 1: low energy, occupied              
            rof(:,:,:,:,:,1) = rof_e(:,:,:,:,:,ie_oc)
                               
 ! index 2: high energy, not occupied 
            do ie_s = 1,nenerg_s
              if( Energ_s(ie_s) > Energ_unoc - eps10 ) exit
            end do

 ! When energy step not uniform, the energy point can be between two points of the energy grid.             
            if( abs( Energ_s(ie_s) - Energ_unoc ) > eps10 ) then
              ie_s = max( ie_s - 1, 1 )
              pds = ( Energ_s(ie_s+1) - Energ_unoc ) / ( Energ_s(ie_s+1) - Energ_s(ie_s) )  
              rof(:,:,:,:,:,2) = pds * rof_e(:,:,:,:,:,ie_s) + ( 1 - pds ) * rof_e(:,:,:,:,:,ie_s+1)
            else
              rof(:,:,:,:,:,2) = rof_e(:,:,:,:,:,ie_s)
            endif  

            if( E1E1 ) secdd_xtal(:,:,:,:,:) = ( 0._db, 0._db )
            if( E1E2 ) secdq_xtal(:,:,:,:,:,:) = ( 0._db, 0._db )
            if( E1M1 ) secmd_xtal(:,:,:,:,:) = ( 0._db, 0._db )
            if( E2E2 ) secqq_xtal(:,:,:,:,:,:) = ( 0._db, 0._db )
            if( E1E3 ) secdo_xtal(:,:,:,:,:,:,:) = ( 0._db, 0._db )
            if( E3E3 ) secoo_xtal(:,:,:,:,:,:) = ( 0._db, 0._db )
            if( M1M1 ) secmm_xtal(:,:,:,:) = ( 0._db, 0._db )

            Atom_done(:) = .false.

            do ia = 1,natomsym
            
              isym = abs( isymeq(ia) )
              call opsym(isym,matopsym)
              Time_reversal = isymeq(ia) < 0 

              if( Atom_done(ia) ) cycle 

              if( icheck > 1 ) then
                if( isym > 0 ) then
                  write(3,150) ia, adjustl( nomsym( isym ) )
                else
                  write(3,160) ia, adjustl( nomsym( isym ) )
                endif
              endif
              
              call Cal_Tau_loss(Ampl_abs,E_loss(ie_loss),Energ_unoc,Energ_s,icheck,isymeq(ia),k_elec_i,k_elec_s, &
                        k_not_possible,Lmax,matopsym,Mod_Q,Moment_conservation,Monocrystal,nenerg_s,ndim2,nlm_f, &
                        nlm_probe,ns_dipmag,nspin,nspino,nspinp,Probed_L,Q_vec,Rot_atom, &
                        Spinorbite,Taull_Spin_channel,Time_reversal,Ylm_comp)

              if( k_not_possible ) cycle Loop_loss

              do i_Spin_channel = 1,n_Spin_channel       ! ----------> Loop over spin transition channels
      
                select case(i_Spin_channel)
                  case(1)
                    Taull(:,:,:,:,:,:,:) = Taull_Spin_channel(:,:,:,:,:,:,:) 
                  case(2)
                    Taull(:,:,:,:,:,:,:) = (0._db,0._db)
                    Taull(:,:,:,:,1,1,:) = Taull_Spin_channel(:,:,:,:,1,1,:)
                  case(3)
                    Taull(:,:,:,:,:,:,:) = (0._db,0._db)
                    Taull(:,:,:,:,1,2,:) = Taull_Spin_channel(:,:,:,:,1,2,:)
                  case(4)
                    Taull(:,:,:,:,:,:,:) = (0._db,0._db)
                    Taull(:,:,:,:,2,1,:) = Taull_Spin_channel(:,:,:,:,2,1,:)
                  case(5)
                    Taull(:,:,:,:,:,:,:) = (0._db,0._db)
                    Taull(:,:,:,:,2,2,:) = Taull_Spin_channel(:,:,:,:,2,2,:)
                end select    

                call Tensor_rixs(coef_g,Core_resolved,FDM_comp_m,Final_optic,Final_tddft,Full_potential, &
                  i_Spin_channel,ich,ip_max,ip0,is_g,lmax,lmoins1,lplus1,lseuil,m_g,Multipole, &
                  n_Ec,n_oo,n_rel,ns_dipmag,ndim2,ninit1,ninit,ninitr,nlm_probe,nlm_p_fp,nspino,nspinp,numat,RIXS,rof, &
                  Rot_atom,secdd,secdo,secdq,secmd,secmm,secoo, &
                  secqq,Spinorbite,Taull,Time_reversal,Write_bav,Ylm_comp)
                
                call Tensor_xtal(Atom_done,i_spin_channel,ia,isymeq,Multipole,n_oo,n_rel, &
                       n_spin_channel,natomsym,ninitr,secdd,secdd_xtal,secdq,secdq_xtal,secdo,secdo_xtal, &
                       secmd,secmd_xtal,secmm,secmm_xtal,secoo,secoo_xtal,secqq,secqq_xtal,Write_bav)

              end do   ! i_spin_channel
              
            end do  ! end of loop on ia
              
            do initr = 1,ninitr       ! ----------> Loop over edges or core states
         
              do ie_fast = 1,nenerg_in
              
                if( .not. Run_fast .and. ie_fast > 1 ) exit
 
                if( Run_fast ) then
                  ie = ie_fast
                else
                  ie = ie_not_fast
                endif

! The factor CFac = Delta_E / (  Energ_in(ie) - Energ_unoc + Shift(initr) + img * Gamma(initr) )
! is replaced by the corresponding integral over energy for this energy step
                Delta_E_inf = E_loss(ie_loss) + E_inf(ie_oc) - Energ_in(ie) - Shift(initr)
                Delta_E_sup = E_loss(ie_loss) + E_sup(ie_oc) - Energ_in(ie) - Shift(initr)

                if( Gamma(initr) > eps10 ) then                 
                  Lorentzian_i = atan( Delta_E_inf / Gamma(initr) ) - atan( Delta_E_sup / Gamma(initr) )
                  Lorentzian_r = 0.5_db * log( ( Gamma(initr)**2 + Delta_E_inf**2 ) / ( Gamma(initr)**2 + Delta_E_sup**2 ) )
                else
                  if( abs( Energ_in(ie) + Shift(initr) - Energ_unoc ) < eps10 ) then
                    if( abs( Delta_E_inf + Delta_E_sup ) < eps10 ) then
                      Lorentzian_r = 0._db
                    elseif( abs( Delta_E_inf ) < eps10 ) then
                      Lorentzian_r = - log( 4._db )   ! limitation on divergence
                    else
                      Lorentzian_r = log( - Delta_E_inf / Delta_E_sup )
                    endif
                  else
                    Lorentzian_r = log( Delta_E_inf / Delta_E_sup )
                  endif
                  if( Energ_in(ie) + Shift(initr) > Energ_unoc - eps10 ) then 
                    Lorentzian_i = pi
                  else
                    Lorentzian_i = 0._db
                  endif
                endif  

                Lorentzian = cmplx( Lorentzian_r, Lorentzian_i, db )                               

                CFac = Lorentzian * Conv_nelec(ie,initr)
                 
                Omega_i = E_i / quatre_mc2 
                Omega_s = E_s / quatre_mc2

                do i_spin_channel = 1,n_spin_channel

                  if( E1E1 ) RIXS_E1E1_tens_xtal(:,:,:) = CFac * Secdd_xtal(:,:,:,initr,i_spin_channel)
                  if( E1E2 ) RIXS_E1E2_tens_xtal(:,:,:,:) = CFac * Secdq_xtal(:,:,:,:,initr,i_spin_channel)
                  if( E1M1 ) RIXS_E1M1_tens_xtal(:,:,:) = CFac * Secmd_xtal(:,:,:,initr,i_spin_channel)
                  if( E2E2 ) RIXS_E2E2_tens_xtal(:,:,:,:) = CFac * Secqq_xtal(:,:,:,:,initr,i_spin_channel)
                  if( E1E3 ) RIXS_E1E3_tens_xtal(:,:,:,:,:) = CFac * Secdo_xtal(:,:,:,:,:,initr,i_spin_channel)
                  if( E3E3 ) RIXS_E3E3_tens_xtal(:,:,:,:) = CFac * Secoo_xtal(:,:,:,:,initr,i_spin_channel)
                  if( M1M1 ) RIXS_M1M1_tens_xtal(:,:) = CFac * Secmm_xtal(:,:,initr,i_spin_channel)

! Calculation of RIXS amplitude
                  call Cal_RIXS_ampl(Amplitude,Dip_rel,E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,M1M1, &
                         Monocrystal,n_rel,npl,Omega_i,Omega_s,pol_i_s,pol_i_p,pol_s_s,pol_s_p, &
                         RIXS_E1E1_tens_xtal,RIXS_E1E2_tens_xtal,RIXS_E1E3_tens_xtal,&
                         RIXS_E2E2_tens_xtal,RIXS_E3E3_tens_xtal,RIXS_E1M1_tens_xtal,RIXS_M1M1_tens_xtal, &
                         Theta,Vec_i,Vec_s)
 
                  ii = initr + ( i_Spin_channel - 1 ) * ninitr 
                  RIXS_ampl(ie_loss,ie,:,i_theta,i_q,ii) = RIXS_ampl(ie_loss,ie,:,i_theta,i_q,ii) + Amplitude(:)   

                end do

              end do   ! end of loop on ie_fast

            end do  ! end of loop over initr      
 
          end do Loop_loss ! end of loop over ie_loss

          write(6,'(a13,f8.3)') ' E_oc - E_F =', Energ_oc(ie_oc) * Rydb 

        end do  ! end of loop over ie_oc
      
      end do ! end of loop over incoming energy ie_fast
              
    end do ! end of loop over i_theta_2
    end do ! end of loop over i_theta_1
  end do ! end of loop over i_q 

! Writing
  call Write_rixs(Angle,Deltar,E_cut_o,E_loss,Energ_in,Epsii,Eseuil(nbseuil),File_name_rixs,Gamma_hole_o,Gamma_max, &
                  npl,jseuil,Mod_k,Monocrystal,n_q_dim,n_Spin_channel,n_theta,ne_loss,nenerg_in,ninit1,ninitr, &
                  nseuil,numat,RIXS_ampl)
  
  deallocate( Amplitude, Angle, Conv_nelec, E_inf, E_sup, Energ_in, Energ_oc, E_loss, Gamma, Mod_k )
  deallocate( Probed_L, RIXS_ampl, rof, rof_e, Shift )
  deallocate( secdd_xtal, secmd_xtal, secmm_xtal, secdq_xtal, secdo_xtal, secqq_xtal, secoo_xtal )
  
  if( E1E1 ) deallocate( RIXS_E1E1_tens )
  if( E1E2 ) deallocate( RIXS_E1E2_tens )
  if( E1E3 ) deallocate( RIXS_E1E3_tens )
  if( E1M1 ) deallocate( RIXS_E1M1_tens )
  if( E2E2 ) deallocate( RIXS_E2E2_tens )
  if( E3E3 ) deallocate( RIXS_E3E3_tens )
  if( M1M1 ) deallocate( RIXS_M1M1_tens )

  return
  110 format(/150('-'),//' RIXS step')
  120 format(/120('-'),//' Run fast, ie_oc =',i4,', ie_loss =',i4,', E_in approximated by',f9.3,' eV')
  130 format(/120('-'),//' ie_not_fast =',i4,', ie_oc =',i4,', ie_loss =',i4,', E_in =',f9.3,' eV')
  140 format('    E_oc =',f7.3,', E_unoc =',f7.3,', E_loss =',f7.3,' eV')
  150 format(/2x,100('-'),//2x,'ia =',i3,', Sym = ',A) 
  160 format(/2x,100('-'),//2x,'ia =',i3,', Sym = ',A,' x Time reversal') 
end

!***********************************************************************

! Energy grid for RIXS

subroutine Energy_Grid_RIXS(E_loss,Energ_in,Energ_oc,Energ_rixs,Energ_s,icheck,ne_loss,nenerg_in,nenerg_in_rixs,nenerg_oc, &
                            nenerg_s,Step_loss)

  use declarations
  implicit none

  integer,intent(in):: icheck, ne_loss, nenerg_in, nenerg_in_rixs, nenerg_oc, nenerg_s
  integer:: ie

  real(kind=db):: Step_loss

  real(kind=db), dimension(nenerg_s), intent(in):: Energ_s
  real(kind=db), dimension(nenerg_in_rixs) :: Energ_rixs

  real(kind=db), dimension(nenerg_in), intent(out):: Energ_in
  real(kind=db), dimension(nenerg_oc), intent(out):: Energ_oc
  real(kind=db), dimension(ne_loss), intent(out):: E_loss

  do ie = 1,nenerg_oc
    Energ_oc(ie) = Energ_s(ie)
  end do
  if( nenerg_in_rixs > 0 ) then
    Energ_in(:) = Energ_rixs(:)
  else
    do ie = 1,nenerg_in
      Energ_in(ie) = Energ_s(ie + nenerg_s - nenerg_in)
    end do
  endif

  do ie = 1,ne_loss
    E_loss(ie) = ie * Step_loss 
  end do

  if( icheck > 1 ) then
    write(3,100)
    write(3,'(/A)') '  Incoming energy grid:'
    write(3,110) Energ_in(:)*rydb
    write(3,'(/A)') '  Energy loss grid:'
    write(3,110) E_loss(:)*rydb
    write(3,'(/A)') '  Occupied state energy grid:'
    write(3,110) Energ_oc(:)*rydb
  end if

  return
  100 format(/' ---- Energy_Grid_RIXS -',100('-'))
  110 format(10f13.7)
end

!***************************************************************

! Calculation of the shifts

subroutine Cal_shift(Epsii,Epsii_ref,Epsii_ref_RIXS_man,icheck,ninit,ninit1,ninitr,Shift)

  use declarations
  implicit none

  integer:: icheck, ninit, ninit1, ninitr
  
  logical:: Epsii_ref_RIXS_man
  
  real(kind=db):: Epsii_moy, Epsii_ref
  real(kind=db), dimension(ninitr):: Epsii, Shift
   
  if( .not. Epsii_ref_RIXS_man ) then
    select case(ninitr)
      case(1)
        Epsii_moy = Epsii(1)
      case(2)
        if( ninit1 /= ninit ) then
          Epsii_moy = Epsii(2)
        else
          Epsii_moy = sum( Epsii(1:ninitr) ) / ninitr
        endif
      case(4,6,10)
 ! in these case ninit = ninitr
        if( ninit1 /= ninit ) then
          Epsii_moy = sum( Epsii(ninit1+1:ninitr) ) / ( ninitr - ninit1 )
        else
          Epsii_moy = sum( Epsii(1:ninitr) ) / ninitr
        endif
    end select
    Epsii_ref = Epsii_moy
    Epsii_ref_RIXS_man = .true.
  endif
  Shift(:) = Epsii(:) - Epsii_ref 

  if( icheck > 0 ) then
    write(3,110) 
    write(3,'(/A)') ' Shift(initl = 1,ninit)'
    write(3,'(10f9.5)') Shift(:) * rydb
  endif  
    
  return
110 format(/'---------- RIXS shift ',80('-'))
end

!*****************************************************************************************************************

! Calculation of the radial matrix rof

subroutine Cal_rof(E_Fermi,Eclie,Energ_s,Eseuil,Full_potential,Hubb_a,Hubb_d,icheck,ip0,ip_max,lmax,lmax_pot, &
                   m_hubb,n_Ec,nbseuil,nenerg_s,ninit1,nlm_pot,nlm_p_fp,nlm_probe,nr,nrm,nspin,nspino,nspinp,numat,psii,r, &
                   Relativiste,Renorm,Rmtg,Rmtsd,rof_e,Spinorbite,V_hubb,V_intmax,V0bdc,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, ie, ip, ip0, ip_max, lmax, lmax_pot, m_hubb, n_Ec, nbseuil, nenerg_s, ninit1, nlm_pot, nlm_p_fp, &
            nlm_probe, nr, nrm, nspin, nspino, nspinp, numat

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,nenerg_s):: rof_e
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: rof

  logical:: Final_tddft, Full_potential, Hubb_a, Hubb_d, NRIXS, Relativiste, Renorm, Spinorbite, Ylm_comp
       
  real(kind=db):: E_Fermi, Eclie, Eimag, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nspin,n_Ec):: Ecinetic
  real(kind=db), dimension(n_Ec):: Enervide
  real(kind=db), dimension(nspin):: V0bdc
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato
  real(kind=db), dimension(nr,ip0:ip_max):: r_or_bess

  Eimag = 0._db
  Final_tddft = .false.
  NRIXS = .false.

  if( icheck > 1 ) write(3,110)

  allocate( rof(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,n_Ec) )

  do ip = ip0,ip_max
    select case(ip)
      case(0)
        r_or_bess(1:nr,ip) = 1._db
      case(1)
        r_or_bess(1:nr,ip) = r(1:nr)
      case default
        r_or_bess(1:nr,ip) = r(1:nr)**ip
    end select
  end do
 
  do ie = 1,nenerg_s

    Enervide(1) = Energ_s(ie) + E_Fermi
    Ecinetic(:,1) = max( Enervide(1) - V0bdc(:), Eclie )

    rof(:,:,:,:,:,:) = (0._db, 0._db)
 
    call radial_RIXS(E_Fermi,Ecinetic,Eimag,Energ_s(ie),Enervide(1),Eseuil,Final_tddft,Full_potential,Hubb_a,Hubb_d,icheck, &
         1,ip_max,ip0,lmax,lmax_pot,m_hubb,nbseuil,ninit1,n_Ec,nlm_pot,nlm_probe,nlm_p_fp,nr,NRIXS,nrm,nspin,nspino, &
         nspinp,numat,psii,r,r_or_bess,Relativiste,Renorm,Rmtg,Rmtsd,rof,Spinorbite,V_hubb,V_intmax,V0bdc,Vrato, &
         Ylm_comp)

    rof_e(:,:,:,:,:,ie) = rof(:,:,:,:,:,1)

  end do
    
  deallocate( rof )

  return
110 format(/'---------- Cal_rof in RIXS ',80('-'))
end

!**********************************************************************

! Calculation of radial integrals 

subroutine radial_RIXS(E_Fermi,Ecinetic,Eimag,Energ,Enervide,Eseuil,Final_tddft,Full_potential,Hubb_a,Hubb_d,icheck, &
         initlv,ip_max,ip0,Lmax,lmax_pot,m_hubb,nbseuil,ninit1,ninitv,nlm_pot,nlma,nlma2,nr,NRIXS,nrm,nspin,nspino, &
         nspinp,numat,psii,r,r_or_bess,Relativiste,Renorm,Rmtg,Rmtsd,rof,Spinorbite,V_hubb,V_intmax,V0bd,Vrato, &
         Ylm_comp)

  use declarations
  implicit none

  integer:: i, ich, icheck, initl, initlv, ip, ip_max, ip0, iseuil, isp, L, l_hubbard, lfin, lm, lmax, &
    lmax_pot, lmp, lp, m, m_hubb, nlm1, nlm2, mp, nbseuil, ninit1, ninitv, nlm, nlm_pot, nlma, nlma2, nr, nrm, &
    nrmtsd, nrmtg, nspin, nspino, nspinp, numat

  character(len=104):: mot

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlma,nlma2,nspinp,nspino,ip0:ip_max,ninitv):: rof
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau

  logical:: Ecomp, Final_tddft, Full_potential, Hubb_a, &
    Hubb_d, Hubb_m, NRIXS, Radial_comp, Relativiste, Renorm, Spinorbite, Ylm_comp

  real(kind=db):: E_Fermi, Eimag, Energ, Enervide, Ephoton, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato

  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nbseuil):: Eseuil, Vecond
  real(kind=db), dimension(nr,nspin):: g0, gm, gpi, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,ip0:ip_max):: r_or_bess
  real(kind=db), dimension(nspinp*(Lmax+1)**2):: E_kin, E_V

  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur

  E_kin = 0._db; E_V(:) = 0._db
  
  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

  if( nbseuil == 2 ) then
    if( Final_tddft .and. nbseuil /= ninitv ) then
      ninit1 = ( ninitv - 2 ) / 2
      if( initlv <= ninit1 ) then
        iseuil = 1
      else
        iseuil = 2
      endif
    else ! not used in DFT (only for TDDFT)
      iseuil = initlv
    endif
  else
    iseuil = 1
  endif

  ich = icheck - 1
  
  if( icheck > 1 ) then
    if( nspin == 1 ) then
      write(3,110) Energ*rydb, Ecinetic(:)*rydb, ( V0bd(:) - E_Fermi )*rydb, real(konde(:))
    else
      write(3,120) Energ*rydb, Ecinetic(:)*rydb, ( V0bd(:) - E_Fermi )*rydb, real(konde(:))
    endif
  endif

  if( Full_potential ) then
    nlm = ( lmax + 1 )**2
  else
    nlm = 1
  endif
  allocate( V(nr,nlm,nlm,nspin) )

  do i = 1,nbseuil
    if( Final_tddft ) then
      if( i /= nbseuil ) cycle
      Ephoton = Energ + Eseuil(nbseuil)
    else
      Ephoton = Energ + Eseuil(i)
    endif
    Ephoton = max( 0.001_db, Ephoton ) ! Optic case
! Multiplicative term for quadrupolar transitions
! In S.I. vecond = k = E*alfa_sf*4*pi*epsilon0 / (e*e)
! In a.u. and Rydberg : k = 0.5_db * alfa_sf * E
    Vecond(i) = 0.5_db * alfa_sf * Ephoton
  end do

  call mod_V(ich,lmax,lmax_pot,nlm,nlm_pot,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)

  gpi(:,:) = 1 / gp(:,:)

  if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10 .or. Ecinetic(nspin) < eps10 ) then
    Ecomp = .true.
  else
    Ecomp = .false.
  endif
!  Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

  if( Full_potential ) then
    lfin = 0
  else
    lfin = lmax
  endif

  do L = 0,lfin

    if( Hubb_a .and. L == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif
    Radial_comp = Ecomp .or. ( Hubb_m .and. Ylm_comp )

    if( Full_potential ) then
      nlm1 = nlm    ! = ( lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*L + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*L + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,ich,konde, &
         L,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

!    Call Kinetic_and_total_energy(E_elecZ,E_Fermi,E_kin2,E_tot,E_V2,Energ,Full_potential,Hubb_m,icheck,L,Lmax,nlm,nlm1, &
!                        nlm2,nr,nspin,nspino,nspinp,r,Rmtsd,Spinorbite,ur,V,Vc,numat)

! Radial integral for the golden rule
    call radial_matrix(Final_tddft,initlv,ip_max,ip0,iseuil,L,nlm1,nlm2,nbseuil, &
           ninitv,nlma,nlma2,nr,NRIXS,nrm,nrmtsd,nspino,nspinp,psii,r,r_or_bess,Radial_comp,Rmtsd,rof,ui,ur,Vecond)

    deallocate( Tau, ui, ur )

  end do   ! end of loop over L

  if( icheck > 1 ) then
    mot = ' '
    if( Ecomp .and. nspino == 2 ) then
      mot(11:17) = 'up sol1'
      mot(37:43) = 'up sol2'
      mot(63:69) = 'dn sol1'
      mot(89:95) = 'dn sol2'
    elseif( nspino == 2 ) then
      mot(5:50) = 'up sol1      up sol2      dn sol1      dn sol2'
    elseif( nspinp == 2 ) then
      mot = '      up           dn'
    else
      mot = ' '
    endif

    do i = 1,nbseuil
      if( Final_tddft ) then
        if( i /= nbseuil ) cycle
        initl = initlv
      else
        initl = i
        if( nbseuil > 1 ) write(3,155) i
      endif

      do ip = ip0,ip_max
        select case(ip)
          case(0)
            Write(3,'(/A)') '    Radial monopole matrix'
          case(1)
            Write(3,'(/A)') '    Radial dipole matrix'
          case(2)
            Write(3,'(/A)') '    Radial quadrupole matrix'
          case(3)
            Write(3,'(/A)') '    Radial octupole matrix'
        end select
        if( Full_potential .or. Hubb_a ) then
          write(3,170) '  L  m lp mp', mot
        else
          write(3,175) '  L  m', mot
        endif
        lm = 0
        do L = 0,lmax
          do m = -L,L
            lm = lm + 1
            if( .not. ( Hubb_a .or. Full_potential .or. Spinorbite ) .and. m /= 0 ) cycle
            lmp = 0
            do lp = 0,lmax
              if( .not. Full_potential .and. L /= lp ) cycle
              do mp = -lp,lp
                if( .not. ( Hubb_a .or. Full_potential .or. m == mp ) ) cycle
                lmp = lp**2 + lp + 1 + mp
                lmp = min(nlma2,lmp)
                if( Full_potential .or. Hubb_a ) then
                  if( Ecomp) then
                    write(3,180) L, m, lp, mp, ( rof(lm,lmp,isp,:,ip,initl), isp = 1,nspinp)
                  else
                    write(3,180) L, m, lp, mp, ( real(rof(lm,lmp,isp,:,ip,initl)), isp = 1,nspinp)
                  endif
                else
                  if( Ecomp) then
                    write(3,185) L, m, ( rof(lm,lmp,isp,:,ip,initl), isp = 1,nspinp )
                  else
                    write(3,185) L, m, ( real(rof(lm,lmp,isp,:,ip,initl)), isp = 1,nspinp )
                  endif
                endif
              end do
            end do
          end do
        end do
      end do

    end do

  endif

  deallocate( V )

  return
  110 format(/' Energ =',f8.3,', Ecinetic =',f8.3,', V0bd =',f8.3,' eV, konde =',f8.5)
  120 format(/' Energ =',f8.3,', Ecinetic =',2f8.3,', V0bd =',2f8.3,' eV, konde =',2f8.5)
  155 format(/' iseuil =',i2)
  170 format(3x,a12,a104)
  175 format(3x,a6,a104)
  180 format(3x,4i3,1p,8e13.5)
  185 format(3x,2i3,1p,8e13.5)
end

!*****************************************************************************************************************

! Not used

! E_kin = E_KS - E_coulomb - E_Vxc
! E_T = E_kin + E_elec + E_Z
! E_elec = 0.5*( E_coulomb - E_Z )
 
subroutine Kinetic_and_total_energy(E_elecZ,E_Fermi,E_kin,E_tot,E_V,Energ,Full_potential,Hubb_m,icheck,L_in,Lmax,nlm,nlm1, &
                          nlm2,nr,nspin,nspino,nspinp,r,Rmtsd,Spinorbite,ur,V,Vc,Z)

  use declarations
  implicit none

  character(len=3), dimension(2):: Spin
  
  integer:: i, icheck, is, iso, isp, L, L_in, L2, Lm, Lmax, m, mv, n, nlm, nlm1, nlm2, np, nr, nspin, nspino, &
            nspinp, Z
  integer, dimension(nspinp):: ispp, mm
  
  logical:: Full_potential, Hubb_m, Spinorbite
  
  real(kind=db):: E_c, E_eZ, E_Fermi, Energ, f_integr3, psi_psi, Rmtsd 
  real(kind=db), dimension(nr):: fct1, fct2, r, r2, Vc, VcZ
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ur
  real(kind=db), dimension(nspinp*(Lmax+1)**2):: E_elecZ, E_kin, E_Tot, E_V
  real(kind=db), dimension(nr,nlm,nlm,nspin):: V

  data Spin/ ' up', ' dn'/

  r2(:) = r(:)**2
  
  VcZ(:) = Vc(:) - 2 * Z / r(:)
  
  do L = 0,Lmax
    if( .not. Full_potential .and. L /= L_in ) cycle

    Lm = nspinp * L**2
              
    L2 = L**2 + L

    do m = -L-1,L    

      if( nlm1 == 1 .and. .not. Spinorbite .and. m /= 0 ) cycle
      if( .not. Spinorbite .and. m == -L-1 ) cycle

      do iso = 1,nspino  

        if( Spinorbite .and. ( ( iso == 1 .and. m == -L-1 ) .or. (  iso == 2 .and. m == L ) ) ) cycle

        if( Spinorbite ) then
          Lm = Lm + 1
          E_c = 0._db
          E_eZ = 0._db
          psi_psi = 0._db
        endif
        
        do isp = 1,nspinp
          is = min( isp, nspin)
          
          if( Spinorbite ) then
            mv = m + isp - 1
            if( mv < - L .or. mv > L ) cycle
            n = L + 1 + mv
          elseif( Hubb_m ) then
            n = L + 1 + m
            Lm = Lm + 1
          else
            Lm = Lm + 1
            n = 1
          endif
          if( nlm2 == nlm1 ) then
            np = n
          else
            np = 1
          endif
          
          if( .not. Spinorbite ) then
            E_c = 0._db
            E_eZ = 0._db
            psi_psi = 0._db
          endif
          
          fct1(:) = ( ur(:,n,np,isp,iso) * r(:) )**2
          psi_psi = psi_psi + quatre_pi * f_integr3(r,fct1,1,nr,Rmtsd)

! E_c
          fct2(1:nr) = fct1(1:nr) * V(1:nr,1,1,is)
          E_c = E_c + quatre_pi * f_integr3(r,fct2,1,nr,Rmtsd) 

! Elec energy + E_Z
          fct2(1:nr) = fct1(1:nr) * VcZ(1:nr)
          E_eZ = E_eZ + 0.5_db * quatre_pi * f_integr3(r,fct2,1,nr,Rmtsd) 

          if( .not. Spinorbite ) then
            E_V(Lm) = E_c / psi_psi - E_Fermi
            E_elecZ(Lm) = E_eZ / psi_psi - E_Fermi
            E_kin(Lm) = Energ - E_V(Lm)
            E_tot(Lm) = E_kin(Lm) + E_elecZ(Lm)
            if( .not. ( Hubb_m .or. Full_potential ) ) then
              do n = 1,2*L
                E_V(Lm+n*nspinp) = E_V(Lm)
                E_elecZ(Lm+n*nspinp) = E_elecZ(Lm)
                E_kin(Lm+n*nspinp) = E_kin(Lm)
                E_tot(Lm+n*nspinp) = E_tot(Lm)
              end do
            endif
          endif

        end do

        if( Spinorbite ) then
          E_V(Lm) = E_c / psi_psi - E_Fermi
          E_elecZ(Lm) = E_eZ / psi_psi - E_Fermi
          E_kin(Lm) = Energ - E_V(Lm)
          E_tot(Lm) = E_kin(Lm) + E_elecZ(Lm)
        endif

      end do
    end do
  end do            

  if( icheck > 1 .and. ( Full_potential .or. L_in == Lmax )  ) then

    write(3,'(/3x,A)') '  L     Kinetic energy       E_Tot           E_V          E_elecZ           k' 
    Lm = 0
    do L = 0,Lmax
      do m = -L-1,L
        if( .not. Spinorbite .and. m == -L-1 ) cycle
        do iso = 1,nspino  
          if( Spinorbite ) then
            if( ( iso == 1 .and. m == -L-1 ) .or. (  iso == 2 .and. m == L ) ) cycle
            Lm = Lm + 1
            n = 0
            do isp = 1,nspinp
              mv = m + isp - 1
              if( mv < - L .or. mv > L ) cycle
              n = n + 1
              mm(n) = mv
              ispp(n) = isp
            end do
            if( n == 1 ) then
              write(3,120) L, ( mm(i), spin(ispp(i)), i = 1,n ), iso, E_kin(Lm)*rydb, E_Tot(Lm)*rydb, E_V(Lm)*rydb, &
                                                                                      E_elecZ(Lm)*rydb, sqrt(E_kin(Lm)) 
            else
              write(3,130) L, ( mm(i), spin(ispp(i)), i = 1,n ), iso, E_kin(Lm)*rydb, E_Tot(Lm)*rydb, E_V(Lm)*rydb, &
                                                                                      E_elecZ(Lm)*rydb, sqrt(E_kin(Lm)) 
            endif
          else
            do isp = 1,nspinp
              Lm = Lm + 1
              if( .not. Spinorbite .and. m /= 0 ) cycle
              if( nspinp == 1 ) then
                write(3,140) L, '   ', E_kin(Lm)*rydb, E_Tot(Lm)*rydb, E_V(Lm)*rydb, E_elecZ(Lm)*rydb, sqrt(E_kin(Lm))
              else
                write(3,140) L, spin(isp), E_kin(Lm)*rydb, E_Tot(Lm)*rydb, E_V(Lm)*rydb, E_elecZ(Lm)*rydb, sqrt(E_kin(Lm))
              endif 
            end do
          endif
        end do
      end do
    end do 
  endif 

  return
  120 format(3x,i3,1x,i3,a3,6x,i3,2x,7f15.5)
  130 format(3x,i3,1x,2(i3,a3),i3,2x,7f15.5)
  140 format(3x,i3,a3,7f15.5)
end

!***************************************************************

! Calculation of the reciprocal vector in the lab frame

subroutine Reciprocal_Vector(icheck,Mat_orth,Trans_rec_dir,Trans_rec_lab)

  use declarations
  implicit none

  integer:: i, icheck
  
  real(kind=db):: Cos_a, Cos_b, Cos_t, Phi, Sin_t, Theta, Vol, Wx_mod, Wy_mod, Wz_mod 
  real(kind=db), dimension(3):: Vec, Vx, Vy, Vz, Wx, Wy, Wx_lab, Wy_lab, Wz, Wz_lab
  real(kind=db), dimension(3,3):: Mat_orth, Trans_rec_dir, Trans_rec_lab 
  
! Diffraction vector in the internal orthogonal basis
  Vx(:) = Mat_orth(:,1)
  Vy(:) = Mat_orth(:,2)
  Vz(:) = Mat_orth(:,3)

! Wx, Wy, Wz : reciprocal cell basis
  call prodvec(Wx,vy,vz)

  vol = sum( Wx(:) * vx(:) )
  Wx(:) = Wx(:) / vol
  call prodvec(Wy,vz,vx)
  Wy(:) = Wy(:) / vol
  call prodvec(Wz,vx,vy)
  Wz(:) = Wz(:) / vol
  Trans_rec_dir(1,:) = Wx(:)
  Trans_rec_dir(2,:) = Wy(:)
  Trans_rec_dir(3,:) = Wz(:)

! Calculation of the matrix to transform the diffraction vector in the Lab basis, at origin, that is before any rotation:
! x: vertical
! y: horizontal, along the incoming beam
! z: horizontal, perpendicular to the incoming beam

  Wx_mod = sqrt( dot_product( Wx, Wx ) )
  Wy_mod = sqrt( dot_product( Wy, Wy ) )
  Wz_mod = sqrt( dot_product( Wz, Wz ) )

  call prodvec(Vec,Wx,Wy)
  Vec(:) = Vec(:) / sqrt( dot_product( Vec, Vec ) )
        
! Basis Vx vertical, Wy in the (incoming, Vx) plane
  Wx_lab(1) = Wx_mod; Wx_lab(2) = 0._db; Wx_lab(3) = 0._db;
  Cos_t = dot_product( Wx, Wy ) / ( Wx_mod * Wy_mod ) 
  Theta = acos( Cos_t )
  Sin_t = sin( Theta )
  Wy_lab(1) = Cos_t * Wy_mod; Wy_lab(2) = sin_t * Wy_mod; Wy_lab(3) = 0._db;

  call prodvec(Vz,Wx,Wy)
  Vz(:) = Vz(:) / sqrt( dot_product( Vz, Vz ) )
  call prodvec(Vy,Vz,Wx)
  Vy(:) = Vy(:) / sqrt( dot_product( Vy, Vy ) )
  Cos_t = dot_product( Wz, Vz ) / Wz_mod
  Theta = acos( Cos_t )
  Sin_t = sin( Theta ) 

  if( abs( Cos_t - 1._db ) < eps10 ) then
    Wz_lab(1) = 0._db; Wz_lab(2) = 0._db; Wz_lab(3) = Wz_mod
  elseif(  abs( Cos_t + 1._db ) < eps10 ) then
    Wz_lab(1) = 0._db; Wz_lab(2) = 0._db; Wz_lab(3) = -Wz_mod
  else
    Cos_a = dot_product( Wz, Vx ) / Wz_mod 
    Cos_b = dot_product( Wz, Vy ) / Wz_mod
    if( abs( Cos_a ) < eps10 ) then
      Wz_lab(1) = 0._db; Wz_lab(2) = Sin_t * Wz_mod; Wz_lab(3) = Cos_t * Wz_mod
    else
      Phi = atan( Cos_b / Cos_a )
      if( Cos_a < 0._db ) Phi = Phi + pi  
      Wz_lab(1) = Sin_t * cos( Phi ) * Wz_mod ; Wz_lab(2) = Sin_t * sin( Phi ) * Wz_mod; Wz_lab(3) = Cos_t * Wz_mod
    endif
  endif

  Trans_rec_lab(:,1) = Wx_lab(:)
  Trans_rec_lab(:,2) = Wy_lab(:)
  Trans_rec_lab(:,3) = Wz_lab(:)
  
  if( icheck > 0 ) then
    write(3,110)
    write(3,'(/A)') ' Reciprocal vector in lab frame:'
    write(3,'(A)') '      A         B        C'
    do i = 1,3
      write(3,'(3f10.5)') Trans_rec_lab(i,:)
    end do
  endif 

  return
110 format(/'---------- Basis in RIXS ',80('-'))
end

!***************************************************************

! Calculation of transormation matrix to have polarizations and wave vectors in crystal basis

subroutine Cal_Mat_rot(icheck,Mat_orth,Q_vec,Rot_sample,Trans_Lab_Dir_Q,Trans_rec_dir,Trans_rec_lab)

  use declarations
  implicit none

  integer:: i, icheck
  
  real(kind=db):: Cos_p, Cos_t, Phi, Q_mod, Rot_sample, Sin_p, Sin_t, Theta 
  real(kind=db), dimension(3):: Q_lab, Q_dir, Q_vec, Vec
  real(kind=db), dimension(3,3):: Mat, Mat_orth, Rot, Rot_1, Rot_2, Trans_Lab_Dir_Q, Trans_rec_dir, Trans_rec_lab
  
  Q_lab(:) = Matmul( Trans_rec_lab, Q_vec )
  Q_mod = sqrt( sum( Q_lab(:)**2 ) )
  if( abs(Q_mod) > eps10 ) Q_lab(:) = Q_lab(:) / Q_mod
    
! Crystal rotation to be with Q along z lab
  if( abs( Q_lab(2) ) > eps10 ) then
    if( abs( Q_lab(3) ) < eps10 .and. Q_lab(2) > eps10 ) then
      Phi = pi / 2
    elseif( abs( Q_lab(3) ) < eps10 .and. Q_lab(2) < - eps10 ) then
      Phi = - pi / 2
    else
      Phi = atan( Q_lab(2) / Q_lab(3) )
      if( Q_lab(2) < 0._db ) Phi = Phi + pi
      Phi = pi / 2 - Phi
    endif
    Cos_p = cos( Phi )
    Sin_p = sin( Phi )
    Rot_1(1,1) = 1._db; Rot_1(1,2) = 0._db; Rot_1(1,3) = 0._db  
    Rot_1(2,1) = 0._db; Rot_1(2,2) = Cos_p; Rot_1(2,3) = - Sin_p  
    Rot_1(3,1) = 0._db; Rot_1(3,2) = Sin_p; Rot_1(3,3) = Cos_p  
  else
    Rot_1(1,1) = 1._db; Rot_1(1,2) = 0._db; Rot_1(1,3) = 0._db  
    Rot_1(2,1) = 0._db; Rot_1(2,2) = 1._db; Rot_1(2,3) = 0._db  
    Rot_1(3,1) = 0._db; Rot_1(3,2) = 0._db; Rot_1(3,3) = 1._db
  endif
  
  Theta = - pi/2 + acos( Q_lab(1) )   ! rotation is negative
  Sin_t = sin( Theta )
  Cos_t = cos( Theta )
  Rot_2(1,1) = Cos_t; Rot_2(1,2) = 0._db; Rot_2(1,3) = Sin_t  
  Rot_2(2,1) = 0._db; Rot_2(2,2) = 1._db; Rot_2(2,3) = 0._db  
  Rot_2(3,1) = - Sin_t; Rot_2(3,2) = 0._db; Rot_2(3,3) = Cos_t
  
  Rot = Matmul( Rot_2, Rot_1 ) 

! Origin of crystal rotation around Q is such that a is perpendicular to the incoming beam
  Q_dir = Matmul( Trans_rec_dir, Q_vec )

  Mat = Matmul( Rot, Mat_orth )
  
  if( abs( Q_dir(2) ) > eps10 .or. abs( Q_dir(3) ) > eps10 ) then
    Vec(1) = 1._db; Vec(2) = 0._db;; Vec(3) = 0._db; 
  else
    Vec(1) = 0._db; Vec(2) = 0._db;; Vec(3) = 1._db; 
  endif  
  Vec = Matmul( Mat, Vec )
  Vec(:) = Vec(:) / sqrt( sum( Vec(:)**2 ) )

  if( abs( Vec(2) ) > eps10 ) then
  
    if( abs( Vec(1) ) > eps10 ) then
      Phi = atan( Vec(2) / Vec(1) )
      if( Vec(2) < 0._db ) Phi = Phi + pi
    elseif( Vec(2) > eps10 ) then
      Phi = - pi / 2
    else
      Phi = pi / 2
    endif
  else
    Phi = 0._db
  endif
  
  Phi = Phi + Rot_sample 

  if( abs( Phi ) > eps10 ) then
    Cos_p = cos( Phi )
    Sin_p = sin( Phi )
    Rot_1(1,1) = Cos_p; Rot_1(1,2) = - Sin_p; Rot_1(1,3) = 0._db  
    Rot_1(2,1) = Sin_p; Rot_1(2,2) = Cos_p;   Rot_1(2,3) = 0._db   
    Rot_1(3,1) = 0._db; Rot_1(3,2) = 0._db;   Rot_1(3,3) = 1._db   

    Rot = Matmul( Rot_1, Rot ) 
  endif
  
  Trans_Lab_Dir_Q = Transpose( Rot )  ! it is to transform vec and pol in lab frame to crystal frame 

  if( icheck > 0 ) then
    write(3,110) Q_vec(:), Rot_sample / radian
    write(3,'(/A)') ' Tr Lab frame - Othonormal unit cell basis'
    do i = 1,3
      write(3,'(3f10.5)') Trans_Lab_Dir_Q(i,:)
    end do 
  endif
  
  return
  110 format(/' Q =',3f10.5,', Azimuthal rotation sample = ',f8.3,' degrees')
end

!***************************************************************

! Calculation of the Polarization and wave vector

subroutine Cal_Pol(icheck,Pol_i_s,Pol_i_p,Pol_s_s,Pol_s_p,Theta,Trans_Lab_Dir_Q,Two_theta_L,Vec_i,Vec_s)

  use declarations
  implicit none

  integer:: i, icheck
  
  real(kind=db):: Cos_2t, Cos_t, Sin_2t, Sin_t, Theta, Two_theta_L 
  real(kind=db), dimension(3):: Pol_i_s, Pol_i_p, Pol_s_s, Pol_s_p, Vec_i, Vec_s
  real(kind=db), dimension(3,3):: Rot, Trans_Lab_Dir_Q
  
  Cos_2t = cos( pi - Two_theta_L )
  Sin_2t = sin( pi - Two_theta_L )

  Cos_2t = cos( Two_theta_L )
  Sin_2t = sin( Two_theta_L )

! Wave vector and polarization in lab frame
  Vec_i(1) = 0._db;   Vec_i(2) = 1._db;  Vec_i(3) = 0._db  
  Vec_s(1) = 0._db;   Vec_s(2) = Cos_2t; Vec_s(3) = Sin_2t
  Pol_i_s(1) = 0._db; Pol_i_s(2) = 0._db; Pol_i_s(3) = 1._db   
  Pol_i_p(1) = 1._db; Pol_i_p(2) = 0._db; Pol_i_p(3) = 0._db   
  Pol_s_s(1) = 0._db; Pol_s_s(2) = - Sin_2t; Pol_s_s(3) = Cos_2t   
  Pol_s_p(1) = 1._db; Pol_s_p(2) = 0._db; Pol_s_p(3) = 0._db   

  Cos_t = cos( pi/2 - Theta )
  Sin_t = sin( pi/2 - Theta )  

  Cos_t = cos( Theta )
  Sin_t = sin( Theta )  

  Rot(1,1) = 1._db; Rot(1,2) = 0._db; Rot(1,3) = 0._db
  Rot(2,1) = 0._db; Rot(2,2) = Cos_t; Rot(2,3) = - Sin_t
  Rot(3,1) = 0._db; Rot(3,2) = Sin_t; Rot(3,3) = Cos_t
  Rot = Transpose( Rot )  ! it is to transform vec and pol in lab frame to crystal frame 
  
  Rot = Matmul( Trans_Lab_Dir_Q, Rot )
  
  Pol_i_s = Matmul( Rot, Pol_i_s )
  Pol_i_s(:) = Pol_i_s(:) / sqrt( sum( Pol_i_s(:)**2 ) ) 
  Pol_i_p = Matmul( Rot, Pol_i_p )
  Pol_i_p(:) = Pol_i_p(:) / sqrt( sum( Pol_i_p(:)**2 ) ) 
  Pol_s_s = Matmul( Rot, Pol_s_s )
  Pol_s_s(:) = Pol_s_s(:) / sqrt( sum( Pol_s_s(:)**2 ) ) 
  Pol_s_p = Matmul( Rot, Pol_s_p )
  Pol_s_p(:) = Pol_s_p(:) / sqrt( sum( Pol_s_p(:)**2 ) ) 
  Vec_i = Matmul( Rot, Vec_i )
  Vec_i(:) = Vec_i(:) / sqrt( sum( Vec_i(:)**2 ) ) 
  Vec_s = Matmul( Rot, Vec_s )
  Vec_s(:) = Vec_s(:) / sqrt( sum( Vec_s(:)**2 ) ) 

  if( icheck > 0 ) then
    write(3,110) Theta/radian, Two_theta_L/radian 
    write(3,'(/A)') '   Pol_i_sig  Pol_i_pi   Pol_s_sig  Pol_s_pi   WaveVec_i  WaveVec_s    (in crytal basis)'
    do i = 1,3
      write(3,'(6f11.5)') Pol_i_s(i), Pol_i_p(i), Pol_s_s(i), Pol_s_p(i), Vec_i(i), Vec_s(i) 
    end do   
  endif
  return
  110 format(/' Theta_in =',f10.5,', 2Theta =',f10.5,' degrees')
end

!***************************************************************

! Calculation of the pseudo inelastic scattering amplitude

subroutine Cal_Tau_loss(Ampl_abs,E_loss,Energ_unoc,Energ_s,icheck,isym,k_elec_i,k_elec_s, &
                        k_not_possible,Lmax,matopsym,Mod_Q,Moment_conservation,Monocrystal,nenerg_s,ndim2,nlm_f, &
                        nlm_probe,ns_dipmag,nspin,nspino,nspinp,Probed_L,Q_vec,Rot_atom, &
                        Spinorbite,Taull,Time_reversal,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, ie_unoc, ie_oc, iso1, iso2, isp, isp1, isp2, isp_f, ispp1, ispp2, iss1, iss2, isym, je, &
            je_unoc, L, L_f, L1, L2, Lm, Lm_f, Lm0, Lm1, Lm1_0, Lm2, Lm2_0, Lmax, Lmax_f, Lmax_u, &
            Lmp, Lms1, Lms2, m, m_F, m1, m2, mp, mp1, mp2, ndim2, ne, nenerg_s, nlm_f, nlm_probe, nlm_u, &
            ns_dipmag, nspin, nspino, nspinp

  complex(kind=db):: Cfac
  complex(kind=db), dimension(nenerg_s,nlm_probe,nspinp,nlm_f,nspinp):: Ampl_abs
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull
  complex(kind=db), dimension(:,:), allocatable:: Tau
  complex(kind=db), dimension(:,:,:), allocatable:: Dlmm
  complex(kind=db), dimension(:,:,:,:), allocatable:: Ampl_oc
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: Ampl_unoc

  logical:: First, k_not_possible, Moment_conservation, Monocrystal, Rotation, Spinorbite, Time_reversal, &
            Title, Ylm_comp
  logical, dimension(0:Lmax):: Probed_L
  logical, dimension(0:Lmax,0:Lmax):: k_usable
  
  real(kind=db):: alfa_i, alfa_s, cos_i, cos_s, d_angle_Q, d_solid, E_loss, Energ_oc, Energ_unoc, &
                  k_elec_i, k_elec_s, Mod_Q, sin_i, sin_s
  real(kind=db), dimension(2):: p
  real(kind=db), dimension(3):: Q_vec, Q_vec_i, V, W
  real(kind=db), dimension(3,3):: Mat_one, Mat_rot, Matopsym, Matopsym_i, Rot_atom
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(:), allocatable:: YlmYlm

  Mat_one(1,1) = 1._db; Mat_one(1,2) = 0._db; Mat_one(1,3) = 0._db
  Mat_one(2,1) = 0._db; Mat_one(2,2) = 1._db; Mat_one(2,3) = 0._db
  Mat_one(3,1) = 0._db; Mat_one(3,2) = 0._db; Mat_one(3,3) = 1._db

  Mat_rot = Mat_one
  Rotation = .false.
     
! Aperture (resolution) on Q angle, taken as 1 degree
  d_angle_Q = 1._db * pi / 180 
  
  k_not_possible = .false.

  Lmax_f = nint( sqrt( real( nlm_f ) ) ) - 1

  if( icheck > 1 ) then
    write(3,'(/5x,A/)') '--- Cal_tau_loss ----------'
    write(3,110) k_elec_s, k_elec_i, Mod_Q
  endif

! One uses the occupied and unoccupied states at respectively the incoming energy (without edge energy)
! and the same minus the energy loss to calculate the angles of the wave vector versus Q  
  Energ_oc = Energ_unoc - E_loss

  do ie_oc = 1,nenerg_s
    if( Energ_s(ie_oc) > Energ_oc - eps10 ) exit
  end do
  
  do ie_unoc = 1,nenerg_s
    if( Energ_s(ie_unoc) > Energ_unoc - eps10 ) exit
  end do
    
  if( abs( Energ_s(ie_unoc) - Energ_unoc ) > eps10 ) then
    ne = 2
    ie_unoc = max( ie_unoc - 1, 1 )
    p(1) = ( Energ_s(ie_unoc+1) - Energ_unoc ) / ( Energ_s(ie_unoc+1) - Energ_s(ie_unoc) )  
    p(2) = 1._db - p(1)  
  else
    ne = 1
    p(1) = 1._db
  endif  

! Second index is incoming 'i' and unoc
  do L1 = 0,Lmax
    do L2 = 0,Lmax
      k_usable(L1,L2) = Probed_L(L1) .and. Probed_L(L2)
    end do
  end do     

  if( Moment_conservation ) then

    First = .true.
  
    if( k_elec_s + k_elec_i < Mod_Q ) then
      if( icheck > 0 ) then
        if( First ) then
          Write(3,*)
          First = .false.
        endif 
        Write(3,120) Energ_unoc*rydb, Energ_oc*rydb, k_elec_s, '+', k_elec_i, '<', Mod_Q
      endif
      k_not_possible = .true.
    endif
    
    if( abs( k_elec_i - k_elec_s ) > Mod_Q ) then
      if( icheck > 0 ) then
        if( First ) then
          Write(3,*)
          First = .false.
        endif 
        Write(3,130) Energ_unoc*rydb, Energ_oc*rydb, k_elec_s, '-', k_elec_i, '>', Mod_Q
      endif
      k_not_possible = .true.
    endif
    
    if( k_not_possible ) then
      if( icheck > 0 ) Write(3,'(5x,A)') 'No possible channel at this energy loss'    
      return
    endif

    do L = Lmax,0,-1
      if( Probed_L(L) ) then
        Lmax_u = L
        exit 
      endif
    end do
    nlm_u = ( Lmax_u + 1 )**2

    cos_s = ( Mod_Q**2 + k_elec_s**2 - k_elec_i**2 ) / ( 2 * k_elec_s * Mod_Q )      
    cos_i = ( Mod_Q - k_elec_s * cos_s ) / k_elec_i
    alfa_s = acos( cos_s )
    if( alfa_s < 0._db ) alfa_s = alfa_s + pi 
    alfa_i = acos( cos_i )
    if( alfa_i < 0._db ) alfa_i = alfa_i + pi 
    sin_s = sin( alfa_s )
    sin_i = sin( alfa_i )
    if( Monocrystal ) then
      d_solid = 2 * pi * d_angle_Q * ( sin_i + sin_s ) / 2 
    else
! Proportion of the 4*pi full solid angle
      d_solid = 0.5_db * d_angle_Q * ( sin_i + sin_s ) / 2
    endif

  else

    Lmax_u = Lmax 
    nlm_u = ( Lmax_u + 1 )**2
    d_solid = 0.5_db * d_angle_Q * 2 / pi   ! 2 / pi is the average of ( sin_i + sin_s ) / 2 
    
  endif
  
  d_solid = d_solid * quatre_pi**2   ! 4*pi factor for relation exp(ikr) versus Bessel

  allocate( Ampl_unoc(ne,nlm_u,nspinp,nlm_f,nspinp) )
  allocate( Ampl_oc(nlm_u,nspinp,nlm_f,nspinp) )
  Ampl_unoc(:,:,:,:,:) = ( 0._db, 0._db )
  Ampl_oc(:,:,:,:) = ( 0._db, 0._db )
   
  if( Moment_conservation .and. Monocrystal ) then

! YlmYlm is the product of both harmonics for phi = 0    
    allocate( YlmYlm(nlm_f) )

    call Cylm_phi0(icheck,Lmax_f,Cos_i,Cos_s,YlmYlm,nlm_f)  

! Rotation of Q_vec to get it in the cluster basis 
    Q_vec_i(:) = matmul( Rot_atom, Q_vec(:) )
    if( abs( Q_vec_i(3) - 1._db ) > eps10 ) then
      V(1) = 0._db; W(2) = 0._db; W(3) = 1._db
    else
      V(1) = 1._db; W(2) = 0._db; W(3) = 0._db
    endif
    call prodvec(W,Q_vec_i,V)
    W(:) = W(:) / sqrt( sum( W(:)**2 ) ) 
    call prodvec(W,Q_vec_i,V)

! Rotation of the atomic basis
    Mat_rot(:,1) = V(:)  
    Mat_rot(:,2) = W(:) 
    Mat_rot(:,3) = Q_vec_i(:)              

    if( abs(isym) /= 1 ) then
      call invermat( Matopsym, Matopsym_i )
      Mat_rot = Matmul( Matopsym_i, Mat_rot )
    endif
    
    Rotation = sum( abs( Mat_rot(:,:) - Mat_one(:,:) ) ) > eps10

  endif
  
  if( Rotation ) then
      
    allocate( Dlmm(-Lmax_u:Lmax_u,-Lmax_u:Lmax_u,0:Lmax_u) )
    call Dlmm_COOP(Dlmm,icheck,Lmax_u,Mat_rot,Ylm_comp)
    Lm = 0
    do L = 0,Lmax_u
      Lm0 = L**2 + L + 1
      do m = -L,L
        Lm = Lm + 1
        do mp = -L,L
          Lmp = Lm0 + mp
          do je = 1,ne
            je_unoc = ie_unoc + je - 1
            Ampl_unoc(je,Lm,:,:,:) = Ampl_unoc(je,Lm,:,:,:) + Dlmm(m,mp,L) * Ampl_abs(je_unoc,Lmp,:,:,:)
          end do
          Ampl_oc(Lm,:,:,:) = Ampl_oc(Lm,:,:,:) + Dlmm(m,mp,L) * Ampl_abs(ie_oc,Lmp,:,:,:)
        end do
        do je = 1,ne
          Ampl_unoc(je,Lm,:,:,:) = p(je) * Ampl_unoc(je,Lm,:,:,:)
        end do 
      end do
    end do
  
    if( icheck > 2 ) then
      write(3,140)
      if( nspinp == 1 ) then
        write(3,150) (( L, m, m = -L,L ), L = 0,Lmax_f )
      else
        write(3,160) ((( L, m, isp, m = -L,L ), L = 0,Lmax_f ), isp = 1,nspinp )
      endif
      do isp = 1,nspinp
        Lm = 0
        do L = 0,Lmax_u
          do m = -L,L
            Lm = Lm + 1
            do je = 1,ne
              write(3,170) L, m, isp, 'un', (Ampl_unoc(je,Lm,isp,1:nlm_f,isp_f), isp_f = 1,nspinp )
            end do
            write(3,170) L, m, isp, 'oc', (Ampl_oc(Lm,isp,1:nlm_f,isp_f), isp_f = 1,nspinp )
          end do
        end do
      end do
    endif

  else
  
    do Lm = 1,nlm_u
      do je = 1,ne
        je_unoc = ie_unoc + je - 1
        Ampl_unoc(je,Lm,:,:,:) = p(je) *  Ampl_abs(je_unoc,Lm,:,:,:)
      end do
      Ampl_oc(Lm,:,:,:) = Ampl_abs(ie_oc,Lm,:,:,:)
    end do
 
  endif
  
  Taull(:,:,:,:,:,:,:) = ( 0._db, 0._db )

  do isp1 = 1,nspinp

    if( Time_reversal ) then
      ispp1 = 3 - isp1
    else
      ispp1 = isp1
    endif
          
    do L1 = 0,Lmax_u
      Lm1_0 = L1**2 + L1 + 1
      do m1 = -L1,L1
        Lm1 = Lm1_0 + m1
        do iso1 = 1,nspino
          if( Spinorbite ) then
            m = m1 + ispp1 - iso1  ! real m
            if( m > L1 .or. m < - L1 ) cycle
            if( Time_reversal ) then
              Lms1 = 2 * ( Lm1_0 - m - 1 ) + 3 - iso1
            else
              Lms1 = 2 * ( Lm1_0 + m - 1 ) + iso1
            endif
          else
            Lms1 = Lm1
          endif

          if( Spinorbite ) then
            iss1 = iso1
          else
            iss1 = min( isp1, nspin )
          endif

          do isp2 = 1,nspinp
      
            if( Time_reversal ) then
              ispp2 = 3 - isp2
            else
              ispp2 = isp2
            endif
            
            do L2 = 0,Lmax_u

              if( .not. k_usable(L1,L2) ) cycle

              Lm2_0 = L2**2 + L2 + 1
              do m2 = -L2,L2
                Lm2 = Lm2_0 + m2
                do iso2 = 1,nspino
                  if( Spinorbite ) then
                    m = m2 + ispp2 - iso2  ! real m
                    if( m > L2 .or. m < - L2 ) cycle
                    if( Time_reversal ) then
                      Lms2 = 2 * ( Lm2_0 - m - 1 ) + 3 - iso2
                    else
                      Lms2 = 2 * ( Lm2_0 + m - 1 ) + iso2
                    endif
                  else
                    Lms2 = Lm2
                  endif

                  if( Spinorbite  ) then
                    iss2 = iso2
                  else
                    iss2 = min( isp2, nspin)
                  endif
      
                  Title = .true.
      
                  do je = 1,ne
                    je_unoc = ie_unoc + je - 1
      
                    do isp_f = 1,nspinp
                      if( .not. Spinorbite .and. isp_f /= isp1 ) cycle
                      if( .not. Spinorbite .and. isp_f /= isp2 ) cycle

                      Lm_f = 0      
                      do L_f = 0,Lmax_f
                        do m_f = -L_f,L_f
                          Lm_f = Lm_f + 1
      
                          if( Monocrystal .and. Moment_conservation ) then
                            Cfac = d_solid * YlmYlm(Lm_f) * Ampl_oc(Lm1,iss1,Lm_f,isp_f) &
                                                          * Conjg( Ampl_unoc(je,Lm2,iss2,Lm_f,isp_f) )
                          else
                            Cfac = d_solid * Ampl_oc(Lm1,iss1,Lm_f,isp_f) * Conjg( Ampl_unoc(je,Lm2,iss2,Lm_f,isp_f) )
                          endif      
                          
                          Taull(Lms1,Lms2,1,1,ispp1,ispp2,1) = Taull(Lms1,Lms2,1,1,ispp1,ispp2,1) + Cfac
                          
                          if( icheck > 2 .and. abs( Cfac ) > eps10 ) then
                            if( Spinorbite) then
                              if( Title ) write(3,210)
                              write(3,220) Lms1, Lm1, L1, m1, iso1, ispp1, Lms2, Lm2, L2, m2, iso2, ispp2, &
                                       L_f, m_f, isp_f, &
                                       Taull(Lms1,Lms2,1,1,ispp1,ispp2,1), &
                                       Cfac, Ampl_oc(Lm1,iss1,Lm_f,isp_f), &                  
                                       Conjg( Ampl_unoc(je,Lm2,iss2,Lm_f,isp_f) ), YlmYlm(Lm_f)
                            else
                              if( Title ) write(3,230)
                              write(3,240) Lm1, L1, m1, ispp1, Lm2, L2, m2, ispp2, L_f, m_f, isp_f, &
                                       Taull(Lms1,Lms2,1,1,ispp1,ispp2,1), &
                                       Cfac, Ampl_oc(Lm1,iss1,Lm_f,isp_f), &                  
                                       Conjg( Ampl_unoc(je,Lm2,iss2,Lm_f,isp_f) )
                            endif
                            Title = .false.
                          endif
      
                        end do
                      end do
                        
                    end do
                  end do
      
                end do
              end do
            end do
          end do  ! end of loop over m2
        end do ! end of loop over L2 
      end do
    end do  ! end of loop over m1  
  end do ! end of loop over L1  

! One comes back to the local cluster basis
  if( Rotation ) then

    allocate( Tau(-Lmax_u:Lmax_u,-Lmax_u:Lmax_u) )

    do L1 = 0,Lmax_u
      Lm1_0 = L1**2 + L1 + 1
      do L2 = 0,Lmax_u
        if( .not. k_usable(L1,L2) ) cycle
        Lm2_0 = L2**2 + L2 + 1
    
        do isp1 = 1,nspinp
          do isp2 = 1,nspinp
            if( .not. Spinorbite .and. isp1 /= isp2 ) cycle
            do iso1 = 1,nspino
              do iso2 = 1,nspino

                do m1 = -L1,L1
                  if( Spinorbite ) then
                    Lms1 = 2 * ( Lm1_0 + m1 - 1 ) + iso1
                  else
                    Lms1 = Lm1_0 + m1
                  endif
                  do m2 = -L2,L2
                    if( Spinorbite ) then
                      Lms2 = 2 * ( Lm2_0 + m2 - 1 ) + iso2
                    else
                      Lms2 = Lm2_0 + m2
                    endif
                    Tau(m1,m2) = Taull(Lms1,Lms2,1,1,isp1,isp2,1)
                    Taull(Lms1,Lms2,1,1,isp1,isp2,1) = ( 0._db, 0._db ) 
                  end do
                end do

          ! Inverse rotation
                do m1 = -L1,L1
                  if( Spinorbite ) then
                    Lms1 = 2 * ( Lm1_0 + m1 - 1 ) + iso1
                  else
                    Lms1 = Lm1_0 + m1
                  endif
                  do m2 = -L2,L2
                    if( Spinorbite ) then
                      Lms2 = 2 * ( Lm2_0 + m2 - 1 ) + iso2
                    else
                      Lms2 = Lm2_0 + m2
                    endif
                    do mp1 = -L1,L1
                      do mp2 = -L2,L2
                        Taull(Lms1,Lms2,1,1,isp1,isp2,1) = Taull(Lms1,Lms2,1,1,isp1,isp2,1) + conjg( Dlmm(mp1,m1,L1) ) &
                                                                                      * Tau(m1,m2) * Dlmm(mp2,m2,L2)
                      end do
                    end do
                  end do  
                end do 
 
              end do
            end do
          end do  
        end do 
      end do  
    end do 

    deallocate( Tau )   
  endif

  if( nspinp == 1 ) Taull(:,:,:,:,2,2,:) = Taull(:,:,:,:,1,1,:) 
    
  if( icheck > 1 ) then
    if( icheck > 2 ) write(3,*)
    write(3,250) p(1:ne)
    
    if( nspinp == 1 ) then
      write(3,260) (( L, m, m = -L,L ), L = 0,Lmax_u )
    elseif( Spinorbite ) then
      write(3,270) (((( L, m, iso1, isp, iso1 = 1,nspino), m = -L,L ), L = 0,Lmax_u ), isp = 1,nspinp )
    else
      write(3,280) ((( L, m, isp, m = -L,L ), L = 0,Lmax_u ), isp = 1,nspinp )
    endif
    
    do isp1 = 1,nspinp
      lm1 = 0
      do L = 0,Lmax_u
        do m = -L, L
          do iso1 = 1,nspino
            lm1 = lm1 + 1
            if( nspinp == 1 ) then
              write(3,290) L, m, (( Taull(lm1,lm2,1,1,isp1,isp2,1), lm2 = 1,nlm_u ), isp2 = 1, nspinp )
            elseif( Spinorbite ) then
              write(3,300) L, m, iso1, isp1, ((( Taull(lm1,2*(lm2-1)+iso2,1,1,isp1,isp2,1), iso2 = 1,nspino), &
                                                 lm2 = 1,nlm_u ), isp2 = 1, nspinp )
            else
              write(3,310) L, m, isp1, (( Taull(lm1,lm2,1,1,isp1,isp2,1), lm2 = 1,nlm_u ), isp2 = 1, nspinp )
            endif
          end do 
        end do
      end do
    end do
  endif
  
  if( Monocrystal .and. Moment_conservation ) deallocate( YlmYlm )
  if( Rotation ) deallocate( Dlmm )
  deallocate( Ampl_oc, Ampl_unoc )
        
  return
  110 format(5x,'k_elec_s =',f10.5,', k_elec_s =',f10.5,', Q =',f10.5)
  120 format(5x,'E_unoc =',f10.5,', E_occ =',f10.5,' eV, k_e_s =',f10.5,' ',a1,' k_e_i =',f10.5,' ',a1,' Q =',f10.5)
  130 format(5x,'E_unoc =',f10.5,', E_occ =',f10.5,' eV, abs( k_e_s =',f10.5,' ',a1,' k_e_i =',f10.5,') ',a1,' Q =',f10.5)
  140 format(/' Atomic amplitude after rotation')
  150 format('  L  m sp',24200(14x,2i3,7x))
  160 format('  L  m sp',24200(12x,3i3,6x))
  170 format(3i3,1x,a2,1p,24200(1x,2e13.5))
  210 format(/' Lms1 Lm1  L1  m1 io1 is1 Lms2 Lm2  L2  m2 io2 is2   Lf  mf isf',12x,'Taull',21x,'dTaull',20x,'Ampl_occ', &
                                                                                    16x,'conj( Ampl_unoc )') 
  220 format(2(1x,6i4),1x,3i4,1p,8(1x,2e13.5)) 
  230 format(/'  Lm1  L1  m1 is1  Lm2  L2  m2 is2   Lf  mf isf',12x,'Taull',21x,'dTaull',20x,'Ampl_occ', &
                                                                                    16x,'conj( Ampl_unoc )') 
  240 format(2(1x,4i4),1x,3i4,1p,8(1x,2e13.5)) 
  250 format(/' Tau inelatic, p =',2f7.4)
  260 format('  L  m',100(12x,2i2,11x))  
  270 format('  L  m io is',200(10x,4i2,9x))  
  280 format('  L  m is',100(11x,3i2,10x))  
  290 format(2i3,1p,100(1x,2e13.5))
  300 format(4i3,1p,200(1x,2e13.5))
  310 format(3i3,1p,100(1x,2e13.5))
  
end

!***********************************************************************

! Calculation of complex Ylm in the case phi = 0
! It is then real
! Only m >= 0 are calculated
! The index of Ylm is lm = L(L+1)/2 + 1 + m

subroutine Cylm_phi0(icheck,Lmax,Cos_i,Cos_s,YlmYlm,nlm)

  use declarations
  implicit none

  integer:: icheck, L, Lm, Lm0, Lm1, Lm1_0, Lm2, Lmax, m, nlm
  
  real(kind=db):: Cos_i, Cos_s, Cot_i, Cot_s, f, g, Sin_i, Sin_s
  real(kind=db), dimension( (Lmax+1)*(Lmax+2)/2 ):: Ylm_i, Ylm_s
  real(kind=db), dimension(nlm):: YlmYlm

  Sin_i = sqrt( 1._db - Cos_i**2 )
  if( abs( Sin_i ) > eps10 ) Cot_i = Cos_i / Sin_i
  Sin_s = sqrt( 1._db - Cos_s**2 )
  if( abs( Sin_s ) > eps10 ) Cot_s = Cos_s / Sin_s
    
! Calcul de Y(0,0) :
  f = 1._db / sqrt( quatre_pi )
  Ylm_i(1) = f
  Ylm_s(1) = f

! Calcul de Y(1,0) :
  if( Lmax > 0 ) then
    f = sqrt( 3 / quatre_pi )
    Ylm_i(2) = sqrt( 3 / quatre_pi ) * Cos_i
    Ylm_s(2) = sqrt( 3 / quatre_pi ) * Cos_s
  endif

! Calcul des Y(L,L) :
  do L = 1,lmax
    Lm = ( (L+1) * (L+2) ) / 2
    Lm1 = ( L * ( L+1) ) / 2
    f = - sqrt( 1._db + 0.5_db/L )
    Ylm_i(Lm) = f * Sin_i * Ylm_i(Lm1)
    Ylm_s(Lm) = f * Sin_s * Ylm_s(Lm1)
  end do

! Calcul des Y(L,L-1) :
  do L = 2,lmax
    lm = ( ( L + 1 ) * ( L + 2 ) ) / 2 - 1
    lm1 = ( L**2 + L ) / 2 - 1
    f = - sqrt( (2*L + 1._db) / (2*L - 2) )
    Ylm_i(Lm) = f * Sin_i * Ylm_i(Lm1)
    Ylm_s(Lm) = f * Sin_s * Ylm_s(Lm1)
  end do

  do L = 2,lmax
    Lm0 = ( L**2 + L ) / 2
    do m = L-2,1,-1
      Lm = Lm0 + m + 1
      f = - 2 * (m + 1._db) / sqrt( L*(L+1._db) - m*(m+1._db) )
      g = - sqrt( ( L*(L+1._db) - (m+1._db)*(m+2._db) ) / ( L*(L+1._db) - m*(m+1._db) ) )
      if( abs( Sin_i ) > eps10 ) then
        Ylm_i(Lm) = Cot_i * f * Ylm_i(Lm+1) + g * Ylm_i(Lm+2)
      else
        Ylm_i(Lm) = 0._db
      endif
      if( abs( Sin_s ) > eps10 ) then
        Ylm_s(Lm) = Cot_s * f * Ylm_s(Lm+1) + g * Ylm_s(Lm+2)
      else
        Ylm_s(Lm) = 0._db
      endif
    end do
  end do

! Calcul des Y(L,0) :
  do L = 2,lmax
    Lm = ( L * ( L + 1 ) ) / 2 + 1
    Lm1 = ( ( L - 1 ) * L ) / 2 + 1
    Lm2 = ( ( L - 2 ) * ( L - 1 ) ) / 2 + 1
    f = sqrt( (2*L + 1._db) / (2*L - 1._db) ) * (2*L - 1._db) / L
    g = - sqrt( (2*L + 1._db) / (2*L - 3._db) ) * (L - 1._db) / L
    Ylm_i(Lm) = f * Cos_i * Ylm_i(Lm1) + g * Ylm_i(Lm2)
    Ylm_s(Lm) = f * Cos_s * Ylm_s(Lm1) + g * Ylm_s(Lm2)
  end do

  do L = 0,Lmax
    Lm0 = L**2 + L + 1
    Lm1_0 = ( L**2 + L ) / 2 + 1
    do m = 0,L
      Lm = Lm0 + m
      Lm1 = Lm1_0 + m
      YlmYlm(Lm) = Ylm_i(Lm1) * Ylm_s(Lm1) 
      if( m > 0 ) YlmYlm(Lm0-m) = YlmYlm(Lm) 
    end do
  end do
  
  if( icheck > 2 ) then
    write(3,'(/A)') ' YlmYlm'
    write(3,'(a8,1x,i2,5x,1000(4x,i4,5x))') '  L  m =', ( m, m = 0,Lmax )
    do L = 0,Lmax
      Lm0 = L**2 + L + 1
      write(3,110) L, YlmYlm(Lm0:Lm0+L)
    end do  
  endif
  
  return
  110 format(i3,1p,1000e13.5)
end

!***********************************************************************

! Calculation of the cartesian tensors

! < g_1 | o*( l_s m_s, k_s, j_s, l_s, irang ) | f_1 > < f_2 | o*( l_i, m_i, k_i, j_i, l_i, jrang  ) | g_2 >  

subroutine Tensor_rixs(coef_g,Core_resolved,FDM_comp_m,Final_optic,Final_tddft,Full_potential, &
                i_Spin_channel,icheck,ip_max,ip0,is_g,lmax,lmoins1,lplus1,Lseuil,m_g,Multipole, &
                n_Ec,n_oo,n_rel,ns_dipmag,ndim2,ninit1,ninit,ninitr,nlm_probe,nlm_p_fp,nspino,nspinp,numat,RIXS,rof, &
                Rot_atom,secdd,secdo,secdq,secmd,secmm,secoo, &
                secqq,Spinorbite,Taull,Time_reversal,Write_bav,Ylm_comp)

  use declarations
  implicit none

  integer:: h_s, h_i, i, i_Spin_channel, icheck, index_cross, initr, ip_max, ip0, ipr, irang, irang1, &
    isp1, isp2, isp3, isp4, ispfg, j, j_i, j_s, jh_i, jh_s, jrang, k, k_i, k_s, &
    l_i, l_s, lm_i, lm_s, lmax, lmomax, lomax, lseuil, m_i, m_s, n_Ec, n_oo, n_rel, ndim2, ninit1, ninit, &
    ninitr, nh_i, nh_s, nj_i, nj_s, nlm_probe, nlm_p_fp, nrang, ns_dipmag, nspino, nspinp, numat

  parameter( lomax = 3, lmomax = ( lomax + 1 )**2 )

  complex(kind=db), dimension(ninitr):: Ten, Tens
  complex(kind=db), dimension(3,3) :: Mat2
  complex(kind=db), dimension(3,3,3) :: Mat3
  complex(kind=db), dimension(3,3,3,3) :: Mat4
  complex(kind=db), dimension(3,3,3,3,3,3) :: Mat6
  complex(kind=db), dimension(3,3,ninitr):: secmm
  complex(kind=db), dimension(3,3,2,ninitr):: secmd
  complex(kind=db), dimension(3,3,3,2,ninitr):: secdq
  complex(kind=db), dimension(3,3,3,3,ninitr):: secqq
  complex(kind=db), dimension(3,3,3,3,2,ninitr):: secdo
  complex(kind=db), dimension(3,3,n_rel,ninitr):: secdd
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr):: secoo
  complex(kind=db), dimension(lmomax,lmomax,n_rel,ninitr):: Tens_lm
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull
  complex(kind=db), dimension(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,n_Ec):: rof

  integer, dimension(ninit,2):: m_g
  integer, dimension(ninit):: is_g

  logical:: Core_resolved, Dip_rel, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, FDM_comp_m, Final_optic, &
    Final_tddft, Full_potential, lmoins1, lplus1, M1M1, RIXS, Spinorbite, Time_reversal, Ylm_comp

  logical, dimension(10):: Multipole
  logical:: Write_bav

  real(kind=db):: c, c0, c1, c12, c120, c3, c5, c8, clm_s, clm_i
  real(kind=db), dimension(0:lomax,3,3,3,lmomax):: clm
  real(kind=db), dimension(3,3):: Rot_atom
  real(kind=db), dimension(ninit,2):: coef_g

  if( icheck > 1 .or. Write_bav ) write(3,'(/5x,A)') '--- Tensor_rixs -----------'

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  Dip_rel = n_rel > 1

  if( E1E3 .or. E3E3 ) then
    nrang = 3
  elseif( E1E2 .or. E2E2 ) then
    nrang = 2
  else
    nrang = 1
  endif
  if( E1M1 .or. M1M1 ) then
    irang1 = 0
  else
    irang1 = 1
  endif

  if( nrang > lomax ) then
    call write_error
    do ipr = 6,9,3
      write(ipr,'(/A)') ' nrang > lomax in coabs.f !'
    end do
    stop
  endif

! Calcul des composantes de la transf base-cartesienne - base-spherique:

  do irang = irang1,nrang

    clm(irang,:,:,:,:) = 0._db

    select case(irang)

      case(0)
! Magnetic dipole
        c = sqrt( 4 * pi )
! Dans ce cas, correspond a l=0,m=0 mais sert a definir x, y ou z.
        clm(irang,1,:,:,4) = c
        clm(irang,2,:,:,2) = c
        clm(irang,3,:,:,3) = c

      case(1)
! Dipole
        c = sqrt( 4 * pi / 3 )
        clm(irang,1,:,:,4) = c
        clm(irang,2,:,:,2) = c
        clm(irang,3,:,:,3) = c

      case(2)
! Quadrupole
        c0 = sqrt( 4 * pi ) / 3
        c = sqrt( 4 * pi / 15 )
        c3 = c / sqrt( 3._db )

        clm(irang,1,1,:,1) = c0
        clm(irang,1,1,:,7) = - c3;   clm(irang,1,1,:,9) = c
        clm(irang,1,2,:,5) = c
        clm(irang,1,3,:,db) = c
        clm(irang,2,2,:,1) = c0
        clm(irang,2,2,:,7) = - c3;   clm(irang,2,2,:,9) = - c
        clm(irang,2,3,:,6) = c
        clm(irang,3,3,:,1) = c0
        clm(irang,3,3,:,7) = 2 * c3

        do i = 1,3
          do j = i+1,3
            clm(irang,j,i,:,:) = clm(irang,i,j,:,:)
          end do
        end do

      case(3)
! Octupole
        c1 = sqrt( 4 * pi / 75 )

        c = sqrt( 4 * pi / 35 )
        c5 = c / sqrt( 5._db )
        c8 = c / sqrt( 2._db )
        c12 = c / sqrt( 3._db )
        c120 = c / sqrt( 30._db )

        clm(irang,1,1,1,4) = 3 * c1
        clm(irang,1,1,1,14) = - 3 * c120; clm(irang,1,1,1,16) = c8
        clm(irang,2,2,2,2) = 3 * c1
        clm(irang,2,2,2,12) = - 3 * c120; clm(irang,2,2,2,10) = - c8
        clm(irang,3,3,3,3) = 3 * c1
        clm(irang,3,3,3,13) = 2 * c5
        clm(irang,1,1,2,2) = c1
        clm(irang,1,1,2,10) = c8;        clm(irang,1,1,2,12) = - c120
        clm(irang,1,2,2,4) = c1
        clm(irang,1,2,2,16) = - c8;      clm(irang,1,2,2,14) = - c120
        clm(irang,1,1,3,3) = c1
        clm(irang,1,1,3,13) = - c5;      clm(irang,1,1,3,15) = c12
        clm(irang,2,2,3,3) = c1
        clm(irang,2,2,3,13) = - c5;      clm(irang,2,2,3,15) = - c12
        clm(irang,2,3,3,2) = c1
        clm(irang,2,3,3,12) = 4 * c120
        clm(irang,1,3,3,4) = c1
        clm(irang,1,3,3,14) = 4 * c120
        clm(irang,1,2,3,11) = c12

        do i = 1,3
          do j = i,3
            do k = j,3
              if( i == j .and. i == k ) cycle
              clm(irang,i,k,j,:) = clm(irang,i,j,k,:)
              clm(irang,j,i,k,:) = clm(irang,i,j,k,:)
              clm(irang,k,i,j,:) = clm(irang,i,j,k,:)
              clm(irang,j,k,i,:) = clm(irang,i,j,k,:)
              clm(irang,k,j,i,:) = clm(irang,i,j,k,:)
            end do
          end do
        end do

    end select

  end do

  secmm(:,:,:) = (0._db,0._db)
  secmd(:,:,:,:) = (0._db,0._db)
  secdq(:,:,:,:,:) = (0._db,0._db)
  secqq(:,:,:,:,:) = (0._db,0._db)
  secdo(:,:,:,:,:,:) = (0._db,0._db)
  secdd(:,:,:,:) = (0._db,0._db)

  if( E3E3 ) secoo(:,:,:,:,:) = (0._db,0._db)

  if( icheck > 1 .and. .not. Final_optic ) write(3,120) lseuil

! Loops over tensor rank
  do irang = irang1,nrang
    do jrang = irang1,nrang

      Tens_lm(:,:,:,:) = (0._db,0._db)
      if( .not. E1E1 .and. irang == 1 .and. jrang == 1 ) cycle
      if( .not. E1E2 .and. ( ( irang == 1 .and. jrang == 2 ) .or. ( irang == 2 .and. jrang == 1 ) ) ) cycle
      if( .not. E2E2 .and. irang == 2 .and. jrang == 2 ) cycle
      if( .not. E1E3 .and. ( ( irang == 1 .and. jrang == 3 ) .or. ( irang == 3 .and. jrang == 1 ) ) ) cycle
      if( .not. M1M1 .and. irang == 0 .and. jrang == 0 ) cycle
      if( .not. E1M1 .and. ( ( irang == 0 .and. jrang == 1 ) .or. ( irang == 1 .and. jrang == 0 ) ) ) cycle
      if( .not. E3E3 .and. irang == 3 .and. jrang == 3 ) cycle
      if( ( irang == 0 .and. jrang > 1 ) .or. ( irang > 1 .and. jrang == 0 ) ) cycle
      if( ( jrang == 3 .and. irang == 2 ) .or. ( jrang == 2 .and. irang == 3 ) ) cycle

      if( icheck > 1 ) write(3,130) irang, jrang

      if( irang < jrang ) then
        index_cross = 1
      else
        index_cross = 2
      endif
      
! Loops over tensor indices
      do k_s = 1,3

        if( irang > 1 ) then
          nj_s = 3
        else
          nj_s = 1
        endif

        do j_s = 1,nj_s

          if( irang == 3 ) then
            nh_s = 3
          else
            nh_s = 1
          endif

          do h_s = 1,nh_s

            jh_s = 3 * ( j_s - 1 ) + h_s

            do k_i = 1,3

              if( jrang > 1 ) then
                nj_i = 3
              else
                nj_i = 1
              endif

              do j_i = 1,nj_i

                if( jrang == 3 ) then
                  nh_i = 3
                else
                  nh_i = 1
                endif

                do h_i = 1,nh_i

                  jh_i = 3 * ( j_i - 1 ) + h_i

                  if( irang == 1 .and. jrang == 1 ) then
                    if( sum( abs( secdd(k_s,k_i,1,:) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 1 .and. jrang == 2 ) then
                    if( sum( abs( secdq(k_s,k_i,j_i,index_cross,:) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 2 .and. jrang == 2 ) then
                    if( sum( abs( secqq(k_s,j_s,k_i,j_i,:) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 1 .and. jrang == 3 ) then
                    if( sum( abs( secdo(k_s,k_i,j_i,h_i,index_cross,:) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 3 .and. jrang == 3 ) then
                    if( sum( abs( secoo(k_s,jh_s,k_i,jh_i,:) ) ) > 1.e-15_db ) cycle
                  endif

! To take into account the relativistic transition channel (just for E1E1).
                  ispfg = 0
                  boucle_dip_pos: do isp1 = 1,2
                  do isp2 = 1,2
                  do isp3 = 1,2
                  do isp4 = 1,2
                    if( lseuil == 0 .and. isp1 /= isp4 ) cycle  ! at K edge, core states are mono-spin. 
                    ispfg = ispfg + 1
                    if( ispfg > 1 .and. ( ( irang /= 1 .or. jrang /= 1 ) .or. .not. dip_rel ) ) exit boucle_dip_pos 

                    Tens(:) = (0._db,0._db)

                    if( icheck > 1 ) write(3,140) k_s, j_s, h_s, k_i, j_i, h_i, ispfg

! Loops over spherical components of the tensors
                  lm_s = 0
                  do l_s = 0,max(irang,1)
                    do m_s = -l_s,l_s
                      lm_s = lm_s + 1

                      clm_s = clm(irang,k_s,j_s,h_s,lm_s)
                      if( abs(clm_s) < eps10 ) cycle

                      lm_i = 0
                      do l_i = 0,max(jrang,1)
                        do m_i = -l_i,l_i
                          lm_i = lm_i + 1

                          clm_i = clm(jrang,k_i,j_i,h_i,lm_i)
                          if( abs(clm_i) < eps10 ) cycle

! Calcul de la composante du tenseur
                          if( sum( abs( Tens_lm(lm_s,lm_i,ispfg,:) ) ) < 1.e-15_db ) then

                            if( icheck > 1 ) write(3,150) l_s, m_s, l_i, m_i, clm_s, clm_i
                            
                            call tens_ab_rixs(coef_g,Core_resolved,Dip_rel,FDM_comp_m,Final_tddft,Full_potential, &
                                icheck,ip_max, &
                                ip0,irang,is_g,isp1,isp2,isp3,isp4,jrang,l_s,m_s,l_i,m_g,m_i,lmax,lmoins1,lplus1,lseuil, &
                                ns_dipmag,ndim2,ninit1,ninit,n_Ec,ninitr,nlm_probe,nlm_p_fp,nspinp,nspino,RIXS,rof, &
                                Spinorbite,Taull,Ten,Time_reversal,Ylm_comp)

                           ! For hydrogen there is only 1 core state
                            if( numat == 1) Ten(:) = Ten(:) / 2

                            Tens_lm(lm_s,lm_i,ispfg,:) = Ten(:)

                          endif

                          Tens(:) = Tens(:) + clm_s * clm_i * Tens_lm(lm_s,lm_i,ispfg,:)
                        end do
                      end do
                    end do
                  end do ! fin boucle l_s

! Remplissage de la valeur calculee dans tous les elements du tenseur equivalent par symetrie.

! M1-M1 (dipole magnetique - dipole magnetique)
                  if( irang == 0 .and. jrang == 0 ) secmm(k_s,k_i,:) = Tens(:)

! M1-E1 (dipole magnetique - dipole electrique)
                  if( ( irang == 0 .and. jrang == 1 ) .or. ( irang == 1 .and. jrang == 0 ) ) &
                    secmd(k_s,k_i,index_cross,:) = Tens(:)

! E1-E1 (Dipole-dipole)
                  if( irang == 1 .and. jrang == 1 ) secdd(k_s,k_i,ispfg,:) = Tens(:)

! Dipole-quadrupole
                  if( ( irang == 1 .and. jrang == 2 ) .or. ( irang == 2 .and. jrang == 1 ) ) then
                    secdq(k_s,k_i,j_i,index_cross,:) = Tens(:)
                    if( index_cross == 1 ) then 
                      secdq(k_s,j_i,k_i,index_cross,:) = Tens(:)
                    else
                      secdq(j_i,k_i,k_s,index_cross,:) = Tens(:)
                    endif
                  endif

! Dipole-octupole
                  if( ( irang == 1 .and. jrang == 3 ) .or. ( irang == 3 .and. jrang == 1 ) ) then
                    secdo(k_s,k_i,j_i,h_i,index_cross,:) = Tens(:)
                    secdo(k_s,k_i,h_i,j_i,index_cross,:) = Tens(:)
                    if( index_cross == 1 ) then
                      secdo(k_s,j_i,k_i,h_i,index_cross,:) = Tens(:)
                      secdo(k_s,h_i,k_i,j_i,index_cross,:) = Tens(:)
                      secdo(k_s,j_i,h_i,k_i,index_cross,:) = Tens(:)
                      secdo(k_s,h_i,j_i,k_i,index_cross,:) = Tens(:)
                    else
                      secdo(j_i,k_i,k_s,h_i,index_cross,:) = Tens(:)
                      secdo(h_i,k_i,k_s,j_i,index_cross,:) = Tens(:)
                      secdo(j_i,k_i,h_i,k_s,index_cross,:) = Tens(:)
                      secdo(h_i,k_i,j_i,k_s,index_cross,:) = Tens(:)
                    endif
                  endif

! Quadrupole-quadrupole
                  if( irang == 2 .and. jrang == 2 ) then
                    secqq(k_s,j_s,k_i,j_i,:) = Tens(:)
                    secqq(j_s,k_s,k_i,j_i,:) = Tens(:)
                    secqq(k_s,j_s,j_i,k_i,:) = Tens(:)
                    secqq(j_s,k_s,j_i,k_i,:) = Tens(:)
                  endif
                  
! Octupole-octupole
                  if( irang == 3 .and. jrang == 3 ) secoo(k_s,jh_s,k_i,jh_i,:) = Tens(:)

                  end do
                  end do
                  end do 
                  end do boucle_dip_pos

                end do
              end do
            end do
          end do
        end do
      end do

      if( icheck < 2 ) cycle

      if( .not. E1E2 .and. irang == 1 .and. jrang == 2 ) cycle
      if( jrang == 3 .and. irang == 2 ) cycle
      if( irang == 0 .and. jrang > 1 ) cycle
      write(3,160)
      lm_s = 0
      do l_s = 0,max(irang,1)
        do m_s = -l_s,l_s
         lm_s = lm_s + 1
         lm_i = 0
         do l_i = 0,max(jrang,1)
            do m_i = -l_i,l_i
              lm_i = lm_i + 1
              do ispfg = 1,n_rel
                if( ( irang /= 1 .or. jrang /= 1 ) .and. ispfg > 1 ) exit
                if( sum( abs(Tens_lm(lm_s,lm_i,ispfg,:)) ) > eps15 ) write(3,170) irang, jrang, l_s, m_s, l_i, m_i, &
                    ispfg, ( Tens_lm(lm_s,lm_i,ispfg,initr), initr = 1,ninitr )
              end do
            end do
          end do
        end do
      end do

    end do
  end do

! Rotation to get the tensors in R1 basis

  do initr = 1,ninitr

    if( E1E1 ) then
      do ispfg = 1,n_rel
        mat2(:,:) = secdd(:,:,ispfg,initr)
        call rot_tensor_2( mat2, Rot_atom )
        secdd(:,:,ispfg,initr) = mat2(:,:)
      end do
    endif

    if( E1E2 ) then
      do index_cross = 1,2
        mat3(:,:,:) = secdq(:,:,:,index_cross,initr)
        call rot_tensor_3( mat3, Rot_atom )
        secdq(:,:,:,index_cross,initr) = mat3(:,:,:)
      end do
    endif

    if( E2E2 ) then
      mat4(:,:,:,:) = secqq(:,:,:,:,initr)
      call rot_tensor_4( mat4, Rot_atom )
      secqq(:,:,:,:,initr) = mat4(:,:,:,:)
    endif

    if( E1E3 ) then
      do index_cross = 1,2
        mat4(:,:,:,:) = secdo(:,:,:,:,index_cross,initr)
        call rot_tensor_4( mat4, Rot_atom )
        secdo(:,:,:,:,index_cross,initr) = mat4(:,:,:,:)
      end do
    endif

    if( E1M1 ) then
      do index_cross = 1,2
        mat2(:,:) = secmd(:,:,initr,index_cross)
        call rot_tensor_2( mat2, Rot_atom )
        secmd(:,:,index_cross,initr) = mat2(:,:)
      end do
    endif

    if( M1M1 ) then
      mat2(:,:) = secmm(:,:,initr)
      call rot_tensor_2( mat2, Rot_atom )
      secmm(:,:,initr) = mat2(:,:)
    endif

    if( E3E3 ) then
      jh_s = 0
      do j_s = 1,3
        do h_s = 1,3
          jh_s = jh_s + 1
          jh_i = 0
          do j_i = 1,3
            do h_i = 1,3
              jh_i = jh_i + 1
              Mat6(:,j_s,h_s,:,j_i,h_i) = secoo(:,jh_s,:,jh_i,initr)
            end do
          end do
        end do
      end do
      call rot_tensor_6( Mat6, Rot_atom )
      jh_s = 0
      do j_s = 1,3
        do h_s = 1,3
          jh_s = jh_s + 1
          jh_i = 0
          do j_i = 1,3
            do h_i = 1,3
              jh_i = jh_i + 1
              secoo(:,jh_s,:,jh_i,initr) = Mat6(:,j_s,h_s,:,j_i,h_i)
            end do
          end do
        end do
      end do
    endif

  end do

  if( Write_bav ) Call Write_tensor_RIXS(i_Spin_channel,E1E1,E2E2,n_rel,ninitr,secdd,secqq)
                                                         
  return
  110 format(/'  --- Tensor_rixs -----------')
  120 format(/' lseuil =',i2)
  130 format(/' -- irang =',i2,', jrang =',i2,' --')
  140 format(/' k_s, j_s, h_s =',3i3,',  k_i, j_i, h_i =',3i3,',  ispfg =',i2)
  150 format(/' l_s, m_s =',2i3,',  l_i, m_i =',2i3,', clm_s, clm_i =',1p, 2e13.5)
  160 format(/' Tensor by harmonics (basis R4) :'/, ' ir jr  l_s m_s  l_i m_i ispfg     Tens(Ylm,i=1,ninitr)')
  170 format(2i3,i4,i3,2x,i4,i3,i5,1p,20e15.7)
  180 format('i',i2,i3,i4,i3,2x,i4,i3,1p,20e15.7)

end

!***********************************************************************

subroutine Tens_ab_rixs(coef_g,Core_resolved,Dip_rel,FDM_comp_m,Final_tddft,Full_potential,icheck,ip_max, &
                              ip0,irang,is_g,isp1,isp2,isp3,isp4,jrang,l_s,m_s,l_i,m_g,m_i,lmax,lmoins1,lplus1,lseuil,ns_dipmag, &
                              ndim2,ninit1,ninit,ninitv,ninitr,nlm_probe,nlm_p_fp,nspinp,nspino,RIXS,rof, &
                              Spinorbite,Taull,Ten,Time_reversal,Ylm_comp)

  use declarations
  implicit none

  integer, intent(in):: icheck, ip_max, ip0, irang, jrang, l_i, l_s, m_i, m_s, lseuil, ndim2, ninit1, &
      ninit, ninitv, ninitr, nlm_probe, nlm_p_fp, ns_dipmag, nspinp, nspino
  integer, dimension(ninit,2), intent(in):: m_g
  integer, dimension(ninit), intent(in):: is_g

  integer:: i_g_1, i_g_2, initl1, initl2, initr, is_dipmag, is_r1, is_r2, iseuil1, iseuil2, iso1, iso2, &
    isp1, isp2, isp3, isp4, &
    ispf1, ispf2, ispinf1, ispinf2, isping1, isping2, ispp_f1, ispp_f2, l_f1, l_f2, li, lm01, lm02, lm_f1, lm_f2, lmax, lmp01, &
    lmp02, lmp_f1, lmp_f2, lms_f1, lms_f2, lp_f1, lp_f2, m_f1, m_f2, mi1, mi2, mp_f1, mp_f2, mv1, mv2

  complex(kind=db):: Cg, rof_1, rof_2, Tau_rad

  complex(kind=db):: Gaunt_i, Gaunt_s, Gaunt_xrc, Gauntmag 
  complex(kind=db), dimension(ninitr):: Ten
  complex(kind=db), dimension(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,ninitv):: rof
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull

  logical:: Core_resolved, Dip_rel, FDM_comp_m, Final_tddft, Full_potential, lmoins1, lplus1, RIXS, &
            Spinorbite, Time_reversal, Titre, Ylm_comp

  real(kind=db):: Ci_1, Ci_2, Ci2, J_initl1, J_initl2, Jz1, Jz2
  real(kind=db), dimension(ninit,2):: coef_g
  
  li = lseuil

  Ten(:) = (0._db,0._db)

! Loop over core states < g_1 ...  > < ... >
  do i_g_1 = 1,ninit

    if( Final_tddft ) then
      initl1 = i_g_1
    else
      initl1 = 1
    endif

    J_initl1 = li + 0.5_db * is_g(i_g_1)
    if( i_g_1 <= ninit1 ) then
      Jz1 = - J_initl1 + i_g_1 - 1
    else
      Jz1 = - J_initl1 + i_g_1 - ninit1 - 1
    endif

    if( i_g_1 <= ninit1 ) then
      iseuil1 = 1
    else
      iseuil1 = 2
    endif

    if( Final_tddft ) then
      initr = 1
    elseif( Core_resolved ) then
      initr = i_g_1
    else
      initr = iseuil1
    endif

    if( ninitv == ninit .and. Final_tddft ) then
      is_r1 = i_g_1
    elseif( RIXS ) then
      is_r1 = 2
    else
      is_r1 = iseuil1
    endif

! Loop over core state spin < g_1 ...  > < ... >
    do isping1 = 1,2 
    
      if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp1 /= isping1 ) cycle

      mi1 = m_g(i_g_1,isping1)
      Ci_1 = Coef_g(i_g_1,isping1)

      if( abs( Ci_1 ) < eps6 ) cycle

      do i_g_2 = 1,ninit

        if( .not. Final_tddft .and. i_g_1 /= i_g_2 ) cycle

        if( Final_tddft ) then
          initl2 = i_g_2
        else
          initl2 = 1
        endif

        J_initl2 = li + 0.5_db * is_g(i_g_2)
        if( i_g_2 <= ninit1 ) then
          Jz2 = - J_initl2 + i_g_2 - 1
        else
          Jz2 = - J_initl2 + i_g_2 - ninit1 - 1
        endif

        if( i_g_2 <= ninit1 ) then
          iseuil2 = 1
        else
          iseuil2 = 2
        endif

        if( ninitv == ninit .and. Final_tddft ) then
          is_r2 = i_g_2
        elseif( RIXS ) then
          is_r2 = 1
        else
          is_r2 = iseuil2
        endif

! Loop over core state spin < ...  > < ... g_2 >
        do isping2 = 1,2  
        
          if( .not. ( RIXS .or. Final_tddft .or. li /= 0 ) .and. isping1 /= isping2 ) cycle
          if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp4 /= isping2 ) cycle

          mi2 = m_g(i_g_2,isping2)
          Ci_2 = Coef_g(i_g_2,isping2)

          if( abs( Ci_2 ) < eps6 ) cycle

          Ci2 = Ci_1 * Ci_2

         if( icheck > 1 ) then
           write(3,120)
           write(3,130) i_g_1, isping1, Jz1, mi1, Ci_1, i_g_2, isping2, Jz2, mi2, Ci_2
         endif
         Titre = .true.

! Loop over spherical harmonics of final states < ... f_1 > < ... >
          do l_f1 = 0,lmax

            if( lplus1 .and. l_f1 < li + 1 ) cycle
            if( lmoins1 .and. l_f1 > li ) cycle
            if( l_f1 > li + irang .or. l_f1 < li - irang .or. mod(l_f1,2) /= mod(li+irang,2) ) cycle

            lm01 = l_f1**2 + l_f1 + 1
            do m_f1 = -l_f1,l_f1
              if( Time_reversal .and. Spinorbite ) then
                lm_f1 = lm01 - m_f1    ! used only for the radial integral rof which is calculated for the not time_reversal
              else
                lm_f1 = lm01 + m_f1
              endif
              if( lm_f1 > nlm_probe ) cycle

! Loop over spin of final states < ... f_1 > < ... >
              do ispinf1 = 1,2

                if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp2 /= ispinf1 ) cycle

                ispf1 = min( ispinf1, nspinp )

! Il n'y a que pour l_s dipole magnetique qu'on peut avoir du spin-flip (ou pour l_s terme dipolaire relativiste)
                if( ispinf1 /= isping1 .and. ( irang > 1 .or. ( irang == 1 .and. &
                                                                           ( .not. Dip_rel .or. jrang /= 1 ) ) ) ) cycle

                if( irang == 0 ) then
                  Gaunt_s = Gauntmag(l_f1,m_f1,ispinf1,m_s,li,mi1,isping1,Ylm_comp,.true.)
                else
                  Gaunt_s = Gaunt_xrc(l_f1,m_f1,l_s,m_s,li,mi1,Ylm_comp)
                endif
                if( abs(Gaunt_s) < eps10 ) cycle

! Loop over spin of final states < ... > < f_2 ... >
                do l_f2 = 0,lmax
                  if( lplus1 .and. l_f2 < li + 1 ) cycle
                  if( lmoins1 .and. l_f2 > li ) cycle
                  if( l_f2 > li + jrang .or. l_f2 < li - jrang .or. mod(l_f2,2) /= mod(li+jrang,2) ) cycle

                  lm02 = l_f2**2 + l_f2 + 1
                  do m_f2 = -l_f2,l_f2
                    lm_f2 = lm02 + m_f2
                    if( Time_reversal .and. Spinorbite ) then
                      lm_f2 = lm02 - m_f2    ! used only for the radial integral rof which is calculated for the not time_reversal
                    else
                      lm_f2 = lm02 + m_f2
                    endif
                    if( lm_f2 > nlm_probe ) cycle

                    do ispinf2 = 1,2  ! spin etat final en sortie
                      ispf2 = min( ispinf2, nspinp )

                      if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp3 /= ispinf2 ) cycle

                      if( ispinf2 /= isping2 .and. ( jrang > 1 .or. ( jrang == 1 &
                                                                             .and. ( .not. Dip_rel .or. irang /= 1 ) ) ) ) cycle 

                      if( .not. ( RIXS .or. Final_tddft .or. Spinorbite ) .and. ispinf2 /= ispinf1 ) cycle

                      if( jrang == 0 ) then
                        Gaunt_i = Gauntmag(l_f2,m_f2,ispinf2,m_i,li,mi2,isping2,Ylm_comp,.true.)
                      else
                        Gaunt_i = Gaunt_xrc(l_f2,m_f2,l_i,m_i,li,mi2,Ylm_comp)
                      endif
                      if( abs(Gaunt_i) < eps10 ) cycle

                      Cg = Ci2 * conjg( Gaunt_s ) * Gaunt_i

                      if( Final_tddft .and. ispinf1 /= isping1 ) then
                        is_dipmag = 2
                      else
                        is_dipmag = 1
                      endif

! Loop over harmonics of non spherical final states < ... f_1 > < ... >
          do lp_f1 = 0,lmax
            if( .not. Full_potential .and. lp_f1 /= l_f1 ) cycle
            lmp01 = lp_f1**2 + lp_f1 + 1

            do mp_f1 = -lp_f1,lp_f1
              if( nlm_p_fp == 1 .and. mp_f1 /= m_f1 ) cycle
              if( Time_reversal .and. Spinorbite ) then
                lmp_f1 = min( lmp01 - mp_f1, nlm_p_fp )
              else
                lmp_f1 = min( lmp01 + mp_f1, nlm_p_fp )
              endif 

! Loop over spin-orbit solution index < ... f_1 > < ... >
              do ispp_f1 = 1,nspinp
                if( .not. Spinorbite .and. ispp_f1 /= ispf1 ) cycle
                iso1 = min( ispp_f1, nspino )

!                if( Spinorbite .and. nlm_p_fp == 1 ) then
!                  mv1 = mp_f1 - ispinf1 + iso1
!                  if( mv1 > l_f1 .or. mv1 < -l_f1 ) cycle
!                else
                  mv1 = mp_f1
!                endif
                lms_f1 = ( lmp01 + mv1 - 1 ) * nspino + iso1

! Loop over harmonics of non spherical final states < ... > < f_2 ... >
                do lp_f2 = 0,lmax
                  if( .not. Full_potential .and. lp_f2 /= l_f2 ) cycle
                  lmp02 = lp_f2**2 + lp_f2 + 1

                  do mp_f2 = -lp_f2,lp_f2
                    if( nlm_p_fp == 1 .and. mp_f2 /= m_f2 ) cycle
                    if( Time_reversal .and. Spinorbite ) then
                      lmp_f2 = min( lmp02 - mp_f2, nlm_p_fp ) 
                    else
                      lmp_f2 = min( lmp02 + mp_f2, nlm_p_fp ) 
                    endif

! Loop over spin-orbit solution index < ... > < f_2 ... >
                    do ispp_f2 = 1,nspinp
                      if( .not. Spinorbite .and. ispp_f2 /= ispf2 ) cycle
                      iso2 = min( ispp_f2, nspino )

!                      if( Spinorbite .and. nlm_p_fp == 1 ) then
!                        mv2 = mp_f2 - ispinf2 + iso2
!                        if( mv2 > l_f2 .or. mv2 < -l_f2 ) cycle
!                      else
                        mv2 = mp_f2
!                      endif
                      lms_f2 = ( lmp02 + mv2 - 1 ) * nspino + iso2

                      if( Time_reversal .and. Spinorbite ) then
                        rof_1 = rof(lm_f1,lmp_f1,3-ispf1,iso1,irang,is_r1)
                        rof_2 = rof(lm_f2,lmp_f2,3-ispf2,iso2,jrang,is_r2)
                      elseif( Time_reversal ) then
                        rof_1 = rof(lm_f1,lmp_f1,3-ispf1,iso1,irang,is_r1)
                        rof_2 = rof(lm_f2,lmp_f2,3-ispf2,iso2,jrang,is_r2)
                      else
                        rof_1 = rof(lm_f1,lmp_f1,ispf1,iso1,irang,is_r1)
                        rof_2 = rof(lm_f2,lmp_f2,ispf2,iso2,jrang,is_r2)
                      endif
 
                      if( FDM_comp_m ) then
                        Tau_rad = rof_1 * rof_2 * Taull(lms_f1,lms_f2,initl1,initl2,ispinf1,ispinf2,is_dipmag)
                      else
                        Tau_rad = conjg( rof_1 ) * rof_2 * Taull(lms_f1,lms_f2,initl1,initl2,ispinf1,ispinf2,is_dipmag)
                      endif
                      
                      Ten(initr) = Ten(initr) + Cg * Tau_rad

                      if( icheck > 1 .and. abs( Tau_rad ) > eps15 ) then
                        if( Titre ) then
                          Titre = .false.
                          if( nlm_p_fp == 1 ) then
                            write(3,131)
                          else
                            write(3,132)
                          endif
                        endif
                        if( nlm_p_fp == 1 ) then
                          write(3,140) l_f1, m_f1, ispinf1, iso1, l_f2, m_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, Gaunt_i, Gaunt_s, &
                            rof_1, rof_2, Taull(lms_f1,lms_f2,initl1,initl2,ispp_f1,ispp_f2,is_dipmag)
                        else
                          write(3,145) l_f1, m_f1, lp_f1, mp_f1, ispinf1, iso1, l_f2, m_f2, lp_f2, mp_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, Gaunt_i, Gaunt_s, &
                            rof_1, rof_2, Taull(lms_f1,lms_f2,initl1,initl2,ispp_f1,ispp_f2,is_dipmag)
                        endif 
                      endif

                    end do ! end of loop over spin-orbit solution f_i ( incoming )
                  end do ! end of loop over mp_f2 non-spherical
                end do ! end of loop over lp_f2 non-spherical

              end do ! end of loop over spin-orbit solution f_s ( outgoing )
            end do ! end of loop over mp_f1 non-spherical
          end do ! end of loop over lp_f1 non-spherical


                    end do ! end of loop over ispinf2, spin of final state ( incoming )
                  end do !  end of loop over  m_f2 final state
                end do !  end of loop over  l_f2 final state

              end do  !  end of loop over spin of final state ( outgoing )
            end do  !  end of loop over  m_f1 final state
          end do !  end of loop over  l_f1 final state

        end do    !  end of loop over spin of core states ( incoming )
      end do    !  end of loop over core states ( incoming )
    end do    !  end of loop over spin of core states ( outgoing )
  end do   !  end of loop over core states  ( outgoing )

! One multiplies by "i" in order the real part of the tensor is the absorption.
  Ten(:) = img * Ten(:)

  return
  120 format(/' ini1 isg1 Jz1 mi1  Ci1  ini2 isg2 Jz2 mi2  Ci2')
  130 format(2(2i4,f6.1,i3,f7.3))
  131 format('  l1  m1 is1 io1  l2  m2 is2 io2',15x,'Ten',27x,'dTen',28x, &
                'Gaunt_i',24x,'Gaunt_s',25x,'rof_s',26x,'rof_i',26x,'Tau_s',26x,'Tau_i')
  132 format('  l1  m1 lp1 mp1 is1 io1  l2  m2 lp2 mp2 is2 io2',15x,'Ten',27x,'dTen',28x, &
                'Gaunt_i',24x,'Gaunt_s',25x,'rof_s',26x,'rof_i',26x,'Tau_s',26x,'Tau_i')
  140 format(8i4,1p,16(1x,2e15.7))
  145 format(12i4,1p,16(1x,2e15.7))
end

!***************************************************************

! Writing of the tensor

subroutine Write_tensor_rixs(i_Spin_channel,E1E1,E2E2,n_rel,ninitr,secdd,secqq)

  use declarations
  implicit none

  integer:: i, i_Spin_channel, initr, j, k, n_rel, ninitr

  complex(kind=db), dimension(3,3,3,3,ninitr):: secqq
  complex(kind=db), dimension(3,3,n_rel,ninitr):: secdd
  
  logical:: E1E1, E2E2
  
  select case(i_Spin_channel)
    case(2)
      write(3,'(/5x,A)') 'Spin channel up-up'
    case(3)
      write(3,'(/5x,A)') 'Spin channel up-dn'
    case(4)
      write(3,'(/5x,A)') 'Spin channel dn-up'
    case(5)
      write(3,'(/5x,A)') 'Spin channel dn-udn'
  end select
  
  if( E1E1 ) then
    if( ninitr > 1 ) then
      write(3,'(/A,10(25x,a7,i2,49x))') '   Tensor E1E1', ( 'initl =',initr, initr = 1,ninitr )
    else
      write(3,'(/A)') '   Tensor E1E1'
    endif
    do i = 1,3
      write(3,'(1p,10(2x,3(1x,2e13.5)))') ( secdd(i,:,1,initr), initr  = 1,ninitr )
    end do
  endif
  if( E2E2 ) then
    write(3,'(/A)') ' Tensor E2E2'
    do k = 1,3
      write(3,120) k, k, k
      do i = 1,3
        write(3,'(1p,10(2x,3(3x,3(1x,2e13.5))))') (( secqq(i,:,j,k,initr), j = 1,3 ), initr = 1,ninitr )
      end do
    end do
  endif

  return
  120 format(40x,'(i,j,1,',i1,')',75x,'(i,j,1,',i1,')',75x,'(i,j,1,',i1,')')
end

!**************************************************************************************************************************************

! RIXS matrix elements for the equivalent atoms

 subroutine Tensor_xtal(Atom_done,i_spin_channel,ia,isymeq,Multipole,n_oo,n_rel, &
                       n_spin_channel,natomsym,ninitr,secdd,secdd_xtal,secdq,secdq_xtal,secdo,secdo_xtal, &
                       secmd,secmd_xtal,secmm,secmm_xtal,secoo,secoo_xtal,secqq,secqq_xtal,Write_bav)

  use declarations
  implicit none

  integer:: i_spin_channel, ia, ib, index_cross, initr, ispfg, isym, je, jhe, jhs, js, he, hs, &
            n_oo, n_rel, n_spin_channel, n_translation, natomsym, ninitr

  integer, dimension(natomsym):: isymeq
  
  complex(kind=db):: Fac
  complex(kind=db), dimension(3,3):: Tens2
  complex(kind=db), dimension(3,3,3):: Tens3
  complex(kind=db), dimension(3,3,3,3):: Tens4
  complex(kind=db), dimension(3,3,3,3,3,3):: Tens6
  complex(kind=db), dimension(3,3,n_rel,ninitr,n_spin_channel):: secdd_xtal
  complex(kind=db), dimension(3,3,ninitr,n_spin_channel):: secmm_xtal
  complex(kind=db), dimension(3,3,2,ninitr,n_spin_channel):: secmd_xtal
  complex(kind=db), dimension(3,3,3,2,ninitr,n_spin_channel):: secdq_xtal
  complex(kind=db), dimension(3,3,3,3,ninitr,n_spin_channel):: secqq_xtal
  complex(kind=db), dimension(3,3,3,3,2,ninitr,n_spin_channel):: secdo_xtal
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,n_spin_channel):: secoo_xtal
  complex(kind=db), dimension(3,3,n_rel,ninitr):: secdd
  complex(kind=db), dimension(3,3,ninitr):: secmm
  complex(kind=db), dimension(3,3,2,ninitr):: secmd
  complex(kind=db), dimension(3,3,3,2,ninitr):: secdq
  complex(kind=db), dimension(3,3,3,3,ninitr):: secqq
  complex(kind=db), dimension(3,3,3,3,2,ninitr):: secdo
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr):: secoo

  complex(kind=db), dimension(:,:,:,:), allocatable:: secdd_x
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: secqq_x

  logical:: All_done, Dip_rel, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, M1M1, Write_bav

  logical, dimension(10):: Multipole
  logical, dimension(natomsym):: Atom_done

  real(kind=db), dimension(3,3):: matopsym

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  Dip_rel = n_rel > 1
  
  if( i_spin_channel == n_spin_channel ) Atom_done(ia) = .true.
  n_translation = 1
  do ib = ia+1,natomsym
    if( Atom_done(ib) ) cycle
    if( isymeq(ia) == isymeq(ib) ) then 
      n_translation = n_translation + 1
      if( i_spin_channel == n_spin_channel ) Atom_done(ib) = .true.
    endif
  end do
      
  isym = abs( isymeq(ia) )
  call opsym(isym,matopsym)

  do initr = 1,ninitr       ! ----------> Loop over edges or core states

! Time reversal already done in Cal_Tau_loss (with inverse rotation on k_elec direction) and Tensor_rixs for rof

    if( E1E1 ) then
      do ispfg = 1,n_rel
        Tens2(:,:) = secdd(:,:,ispfg,initr)
        if( isym /= 1 ) call rot_tensor_2( Tens2, matopsym )
        secdd_xtal(:,:,ispfg,initr,i_Spin_channel) = secdd_xtal(:,:,ispfg,initr,i_Spin_channel) + n_translation * Tens2(:,:)
      end do
    endif

    if( E1E2 ) then
      do index_cross = 1,2
   ! the img of (i/2).kr is set now.
        if( index_cross == 1 ) then
          Fac = img
        else
          Fac = - img
        endif
        Tens3(:,:,:) = secdq(:,:,:,index_cross,initr)
        if( isym /= 1 ) call rot_tensor_3( Tens3, matopsym )
   ! contrary to xanes or DAFS, the E1E2 tensor are multiplied by img in cal_tau_loss
        secdq_xtal(:,:,:,index_cross,initr,i_spin_channel) = secdq_xtal(:,:,:,index_cross,initr,i_spin_channel) &
                                                           + Fac * n_translation * Tens3(:,:,:)
      end do

    endif

    if( E2E2 ) then
      Tens4(:,:,:,:) = secqq(:,:,:,:,initr)
      if( isym /= 1 ) call rot_tensor_4( Tens4, matopsym )
      secqq_xtal(:,:,:,:,initr,i_spin_channel) = secqq_xtal(:,:,:,:,initr,i_spin_channel) + n_translation * Tens4(:,:,:,:)
    endif

    if( E1E3 ) then
      do index_cross = 1,2
        Tens4(:,:,:,:) = secdo(:,:,:,:,index_cross,initr)
        if( isym /= 1 ) call rot_tensor_4( Tens4, matopsym )
        secdo_xtal(:,:,:,:,index_cross,initr,i_spin_channel) = secdo_xtal(:,:,:,:,index_cross,initr,i_spin_channel) &
                                                             + n_translation * Tens4(:,:,:,:)
      end do
    endif

    if( E3E3 ) then
      jhe = 0
      do je = 1,3
        do he = 1,3
          jhe = jhe + 1
          jhs = 0
          do js = 1,3
            do hs = 1,3
              jhs = jhs + 1
              Tens6(:,je,he,:,js,hs) = secoo(:,jhe,:,jhs,initr)
            end do
          end do
        end do
      end do
      if( isym /= 1 ) call rot_tensor_6( Tens6, Matopsym )
      jhe = 0
      do je = 1,3
        do he = 1,3
          jhe = jhe + 1
          jhs = 0
          do js = 1,3
            do hs = 1,3
              jhs = jhs + 1
              secoo_xtal(:,jhe,:,jhs,initr,i_spin_channel) = secoo_xtal(:,jhe,:,jhs,initr,i_spin_channel) &
                                                           + n_translation * Tens6(:,je,he,:,js,hs)
            end do
          end do
        end do
      end do
    endif

    if( E1M1 ) then
      do index_cross = 1,2
        Tens2(:,:) = - secmd(:,:,index_cross,initr)   ! note the change of sign
        if( isym /= 1 ) call rot_tensor_2( Tens2, matopsym )
   ! Plane and inverse symmetries change the sign
        if( ( isym >= 25 .and. isym <= 48 ) .or. ( isym >= 53 .and. isym <= 57 ) .or. isym == 59  .or. isym == 61 &
                        .or. isym == 63 ) Tens2(:,:) = - Tens2(:,:)
   
        secmd_xtal(:,:,index_cross,initr,i_spin_channel) = secmd_xtal(:,:,index_cross,initr,i_spin_channel) &
                                                         + n_translation * Tens2(:,:)
      end do
    endif

    if( M1M1 ) then
      Tens2(:,:) = secmm(:,:,initr)
      if( isym /= 1 ) call rot_tensor_2( Tens2, matopsym )
      secmm_xtal(:,:,initr,i_spin_channel) = secmm_xtal(:,:,initr,i_spin_channel) + n_translation * Tens2(:,:)
    endif

  end do ! end of loop over initr
 
  All_done = .true.
  do ib = 2,natomsym
    if( Atom_done(ib) ) cycle
    All_done = .false.
    exit
  end do
  
  if( Write_bav .and. All_done ) then
    write(3,'(/5x,A)') '--- Crystal tensor --------'
    allocate( secdd_x(3,3,n_rel,ninitr) )
    allocate( secqq_x(3,3,3,3,ninitr) )
    if( E1E1 ) secdd_x(:,:,:,:) = secdd_xtal(:,:,:,:,i_spin_channel)  
    if( E2E2 ) secqq_x(:,:,:,:,:) = secqq_xtal(:,:,:,:,:,i_spin_channel)  
    Call Write_tensor_RIXS(i_Spin_channel,E1E1,E2E2,n_rel,ninitr,secdd_x,secqq_x)
    deallocate( secdd_x, secqq_x )
  endif
                                                         
  return
end

!*****************************************************************************************************************

subroutine Cal_RIXS_ampl(Amplitude,Dip_rel,E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,M1M1, &
                         Monocrystal,n_rel,npl,Omega_i,Omega_s,pol_i_s,pol_i_p,pol_s_s,pol_s_p, &
                         RIXS_E1E1_tens_xtal,RIXS_E1E2_tens_xtal,RIXS_E1E3_tens_xtal,&
                         RIXS_E2E2_tens_xtal,RIXS_E3E3_tens_xtal,RIXS_E1M1_tens_xtal,RIXS_M1M1_tens_xtal, &
                         Theta,Vec_i,Vec_s)

  use declarations
  implicit none

  integer:: he, hs, ipl, isp1, isp2, isp3, isp4, ispfg, j1, je, jhe, jhs, js, ke, ks, n_rel, npl  
   
  complex(kind=db):: matpole, matpols
  complex(kind=db), dimension(npl):: Amplitude
  complex(kind=db), dimension(3):: pol_i, pol_s, u_i, u_s
  complex(kind=db), dimension(3,3):: RIXS_M1M1_tens_xtal
  complex(kind=db), dimension(3,3,2):: RIXS_E1M1_tens_xtal
  complex(kind=db), dimension(3,3,n_rel):: RIXS_E1E1_tens_xtal
  complex(kind=db), dimension(3,3,3,2):: RIXS_E1E2_tens_xtal
  complex(kind=db), dimension(3,3,3,3):: RIXS_E2E2_tens_xtal, RIXS_E3E3_tens_xtal
  complex(kind=db), dimension(3,3,3,3,2):: RIXS_E1E3_tens_xtal
  
  logical:: Dip_rel, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, M1M1, Monocrystal
  
  real(kind=db):: Fac, Omega_i, Omega_s, sin2_2t, Theta
  real(kind=db), dimension(3):: pol_i_s, pol_i_p, pol_s_s, pol_s_p, Vec_i, Vec_s

  Amplitude(:) = (0._db, 0._db )
    
  if( .not. Monocrystal ) then
 
    if( E1E1 ) then
      sin2_2t = sin( 2 * Theta )**2 
      Fac = sqrt( ( 1._db / 5._db ) * ( 1 - ( 2._db / 3._db ) * sin2_2t ) + 1._db / 15._db )  
      do ipl = 1,3
        Amplitude(ipl) = Amplitude(ipl) + Fac * RIXS_E1E1_tens_xtal(ipl,ipl,1)
      end do
      Fac = sqrt( ( 1._db / 30._db ) * ( 2 + 3 * sin2_2t ) )  
      Amplitude(4) = Amplitude(4) + Fac * ( RIXS_E1E1_tens_xtal(1,2,1) + RIXS_E1E1_tens_xtal(2,1,1) )
      Amplitude(5) = Amplitude(5) + Fac * ( RIXS_E1E1_tens_xtal(1,3,1) + RIXS_E1E1_tens_xtal(3,1,1) )
      Amplitude(6) = Amplitude(6) + Fac * ( RIXS_E1E1_tens_xtal(2,3,1) + RIXS_E1E1_tens_xtal(3,2,1) )
! I do the hypothesis that there is no spin crossing contribution. I am not sure.
      if( n_rel /= 1 ) then
        Fac = sqrt( ( 1._db / 5._db ) * ( 1 - ( 2._db / 3._db ) * sin2_2t ) + 1._db / 15._db )  
        do ipl = 1,3
          Amplitude(ipl) = Amplitude(ipl) + Fac * RIXS_E1E1_tens_xtal(ipl,ipl,n_rel)
        end do
        Fac = sqrt( ( 1._db / 30._db ) * ( 2 + 3 * sin2_2t ) )  
        Amplitude(4) = Amplitude(4) + Fac * ( RIXS_E1E1_tens_xtal(1,2,n_rel) + RIXS_E1E1_tens_xtal(2,1,n_rel) )
        Amplitude(5) = Amplitude(5) + Fac * ( RIXS_E1E1_tens_xtal(1,3,n_rel) + RIXS_E1E1_tens_xtal(3,1,n_rel) )
        Amplitude(6) = Amplitude(6) + Fac * ( RIXS_E1E1_tens_xtal(2,3,n_rel) + RIXS_E1E1_tens_xtal(3,2,n_rel) )
      endif
    endif

! not checked, I just take the spherical contribution
    if( E2E2 ) then
      Fac = 1 / 45._db
      Fac = Fac / ( 2._db * sqrt( 5._db) ) 
      Amplitude(:) = Amplitude(:) &
                   + Fac * ( 6 * ( RIXS_E2E2_tens_xtal(1,3,1,3) + RIXS_E2E2_tens_xtal(2,3,2,3) + RIXS_E2E2_tens_xtal(1,2,1,2) ) &
                           + 2 * ( RIXS_E2E2_tens_xtal(1,1,1,1) + RIXS_E2E2_tens_xtal(2,2,2,2) + RIXS_E2E2_tens_xtal(3,3,3,3) ) &
                               - ( RIXS_E2E2_tens_xtal(1,1,2,2) + RIXS_E2E2_tens_xtal(1,1,3,3) + RIXS_E2E2_tens_xtal(3,3,2,2) &
                                 + RIXS_E2E2_tens_xtal(2,2,1,1) + RIXS_E2E2_tens_xtal(3,3,1,1) + RIXS_E2E2_tens_xtal(2,2,3,3) ) )
    endif

    if( M1M1 ) then
      sin2_2t = sin( 2 * Theta )**2 
      Fac = sqrt( ( 1._db / 5._db ) * ( 1 - ( 2._db / 3._db ) * sin2_2t ) + 1._db / 15._db )  
      do ipl = 1,3
        Amplitude(ipl) = Amplitude(ipl) + Fac * RIXS_M1M1_tens_xtal(ipl,ipl)
      end do
      Fac = sqrt( ( 1._db / 30._db ) * ( 2 + 3 * sin2_2t ) )  
      Amplitude(4) = Amplitude(4) + Fac * ( RIXS_M1M1_tens_xtal(1,2) + RIXS_M1M1_tens_xtal(2,1) )
      Amplitude(5) = Amplitude(5) + Fac * ( RIXS_M1M1_tens_xtal(1,3) + RIXS_M1M1_tens_xtal(3,1) )
      Amplitude(6) = Amplitude(6) + Fac * ( RIXS_M1M1_tens_xtal(2,3) + RIXS_M1M1_tens_xtal(3,2) )
    endif

    if( E3E3 ) then  ! not checked : taken like E2E2...
      Fac = 1 / 45._db
      Fac = Fac / ( 2._db * sqrt( 5._db) ) 
      Amplitude(:) = Amplitude(:) &
                   + Fac * ( 6 * ( RIXS_E3E3_tens_xtal(1,3,1,3) + RIXS_E3E3_tens_xtal(2,3,2,3) + RIXS_E3E3_tens_xtal(1,2,1,2) ) &
                           + 2 * ( RIXS_E3E3_tens_xtal(1,1,1,1) + RIXS_E3E3_tens_xtal(2,2,2,2) + RIXS_E3E3_tens_xtal(3,3,3,3) ) &
                               - ( RIXS_E3E3_tens_xtal(1,1,2,2) + RIXS_E3E3_tens_xtal(1,1,3,3) + RIXS_E3E3_tens_xtal(3,3,2,2) &
                                 + RIXS_E3E3_tens_xtal(2,2,1,1) + RIXS_E3E3_tens_xtal(3,3,1,1) + RIXS_E3E3_tens_xtal(2,2,3,3) ) )
    endif

  else    ! Q_rsolved

    do ipl = 1,npl
 
      select case(ipl)
        case(1)
          pol_i(:) = cmplx( pol_i_s(:), 0._db, db )
          pol_s(:) = cmplx( pol_s_s(:), 0._db, db )
        case(2)
          pol_i(:) = cmplx( pol_i_s(:), 0._db, db )
          pol_s(:) = cmplx( pol_s_p(:), 0._db, db )
        case(3)
          pol_i(:) = cmplx( pol_i_p(:), 0._db, db )
          pol_s(:) = cmplx( pol_s_s(:), 0._db, db )
        case(4)
          pol_i(:) = cmplx( pol_i_p(:), 0._db, db )
          pol_s(:) = cmplx( pol_s_p(:), 0._db, db )
        case(5)
          pol_i(:) = ( pol_i_s(:) + img * pol_i_p(:) ) / sqrt( 2._db ) 
          pol_s(:) = cmplx( pol_s_s(:), 0._db, db )
        case(6)
          pol_i(:) = ( pol_i_s(:) - img * pol_i_p(:) ) / sqrt( 2._db ) 
          pol_s(:) = cmplx( pol_s_s(:), 0._db, db )
        case(7)
          pol_i(:) = ( pol_i_s(:) + img * pol_i_p(:) ) / sqrt( 2._db ) 
          pol_s(:) = cmplx( pol_s_p(:), 0._db, db )
        case(8)
          pol_i(:) = ( pol_i_s(:) - img * pol_i_p(:) ) / sqrt( 2._db ) 
          pol_s(:) = cmplx( pol_s_p(:), 0._db, db )
      end select
    
      if( M1M1 .or. E1M1 ) then
        u_i(1) = Vec_i(2) * Pol_i(3) - Vec_i(3) * Pol_i(2)
        u_i(2) = Vec_i(3) * Pol_i(1) - Vec_i(1) * Pol_i(3)
        u_i(3) = Vec_i(1) * Pol_i(2) - Vec_i(2) * Pol_i(1)
        u_s(1) = Vec_s(2) * Pol_s(3) - Vec_s(3) * Pol_s(2)
        u_s(2) = Vec_s(3) * Pol_s(1) - Vec_s(1) * Pol_s(3)
        u_s(3) = Vec_s(1) * Pol_s(2) - Vec_s(2) * Pol_s(1)
      endif
    
      if( E1E1 ) then
    
        if( .not. Dip_rel ) then
    
          do ke = 1,3
            Amplitude(ipl) = Amplitude(ipl) + pol_i(ke) * sum( conjg(pol_s(:)) *  RIXS_E1E1_tens_xtal(:,ke,1) )
          end do
    
        else
    
          ispfg = 0
          do isp1 = 1,2
            do isp2 = 1,2
              do isp3 = 1,2
                do isp4 = 1,2
                  if( n_rel == 8 .and. isp1 /= isp4 ) cycle 
                  ispfg = ispfg + 1
      
                  do ks = 1,3
                    if( isp1 == isp2 ) then
                      matpols = Pol_s(ks)
                    else
                      matpols = (0._db, 0._db)
                    endif
                    if( isp1 == 1 .and. isp2 == 1 ) then
                      if( ks == 1 ) then
                        matpols = matpols - img * Omega_s * Pol_s(2)
                      elseif( ks == 2 ) then
                        matpols = matpols + img * Omega_s * Pol_s(1)
                      endif
                    elseif( isp1 == 1 .and. isp2 == 2 ) then
                      if( ks == 1 ) then
                        matpols = matpols + Omega_s * Pol_s(3)
                      elseif( ks == 2 ) then
                        matpols = matpols - img * Omega_s * Pol_s(3)
                      else
                        matpols = matpols + img * Omega_s * ( Pol_s(2) + img * Pol_s(1) )
                      endif
                    elseif( isp1 == 2 .and. isp2 == 1 ) then
                      if( ks == 1 ) then
                        matpols = matpols - Omega_s * Pol_s(3)
                      elseif( ks == 2 ) then
                        matpols = matpols - img * Omega_s * Pol_s(3)
                      else
                        matpols = matpols + img * Omega_s * ( Pol_s(2) - img * Pol_s(1) )
                      endif
                    elseif( isp1 == 2 .and. isp2 == 2 ) then
                      if( ks == 1 ) then
                        matpols = matpols + img * Omega_s * Pol_s(2)
                      elseif( ks == 2 ) then
                        matpols = matpols - img * Omega_s * Pol_s(1)
                      endif
                    endif                      
      
                    do ke = 1,3
                      if( isp3 == isp4 ) then
                        matpole = Pol_i(ke)
                      else
                        matpole = (0._db, 0._db)
                      endif
                      if( isp3 == 1 .and. isp4 == 1 ) then
                        if( ke == 1 ) then
                          matpole = matpole - img * Omega_i * Pol_i(2)
                        elseif( ke == 2 ) then
                          matpole = matpole + img * Omega_i * Pol_i(1)
                        endif
                      elseif( isp3 == 1 .and. isp4 == 2 ) then
                        if( ke == 1 ) then
                          matpole = matpole + Omega_i * Pol_i(3)
                        elseif( ke == 2 ) then
                          matpole = matpole - img * Omega_i * Pol_i(3)
                        else
                          matpole = matpole + img * Omega_i * ( Pol_i(2) + img * Pol_i(1) )
                        endif
                      elseif( isp3 == 2 .and. isp4 == 1 ) then
                        if( ke == 1 ) then
                          matpole = matpole - Omega_i * Pol_i(3)
                        elseif( ke == 2 ) then
                          matpole = matpole - img * Omega_i * Pol_i(3)
                        else
                          matpole = matpole + img * Omega_i * ( Pol_i(2) - img * Pol_i(1) )
                        endif
                      elseif( isp3 == 2 .and. isp4 == 2 ) then
                        if( ke == 1 ) then
                          matpole = matpole + img * Omega_i * Pol_i(2)
                        elseif( ke == 2 ) then
                          matpole = matpole - img * Omega_i * Pol_i(1)
                        endif
                      endif
      
                      Amplitude(ipl) = Amplitude(ipl) + conjg( matpols ) * matpole * RIXS_E1E1_tens_xtal(ks,ke,ispfg)
                    
                    end do
                  end do
      
                end do
              end do
            end do
          end do
    
        endif
    
      endif
    
      if( E1E2 ) then
        do ke = 1,3
          do ks = 1,3   ! img was included in the RIXS_E1E2_tens_xtal calculation 
            Amplitude(ipl) = Amplitude(ipl) + conjg( Pol_s(ks) ) * Pol_i(ke) &
                      * sum( Vec_i(:) * RIXS_E1E2_tens_xtal(ke,ks,:,1) + Vec_s(:) * RIXS_E1E2_tens_xtal(ke,ks,:,2) )
          end do
        end do
      endif
    
      if( E2E2 ) then
        do ke = 1,3
          do je = 1,3
            do ks = 1,3
              Amplitude(ipl) = Amplitude(ipl) &
                             + conjg( Pol_s(ks) ) * Pol_i(ke) * Vec_i(je) * sum( Vec_s(:) * RIXS_E2E2_tens_xtal(ks,:,ke,je) )
            end do
          end do
        end do
      endif
    
      if( E1E3 ) then
        do ke = 1,3
          do ks = 1,3
            do j1 = 1,3
              Amplitude(ipl) = Amplitude(ipl) + conjg( Pol_s(ks) ) * Pol_i(ke) &
                        * ( Vec_i(j1) * sum( Vec_i(:) * RIXS_E1E3_tens_xtal(ks,ke,j1,:,1) ) &
                          + Vec_s(j1) * sum( Vec_s(:) * RIXS_E1E3_tens_xtal(ks,ke,j1,:,2) ) )
            end do
          end do
        end do
      endif
      
      if( E3E3 ) then
        do ke = 1,3
          do je = 1,3
            do he = 1,3
              jhe = 3 * ( je - 1 ) + he
              do ks = 1,3
                do js = 1,3
                  do hs = 1,3
                    jhs = 3 * ( js - 1 ) + hs
                    Amplitude(ipl) = Amplitude(ipl) &
                                   + conjg( Pol_s(ks) ) * Vec_s(js) * Vec_s(hs) * Pol_i(ke) * Vec_i(je) * Vec_i(he) &
                                   * RIXS_E3E3_tens_xtal(ks,jhs,ke,jhe)
                  end do
                end do
              end do
            end do
          end do
        end do
      endif
    
      if( M1M1 ) then
        do ke = 1,3
          Amplitude(ipl) = Amplitude(ipl) + u_i(ke) * sum( conjg(u_s(:)) * RIXS_M1M1_tens_xtal(:,ke) )
        end do
      endif
    
      if( E1M1 ) then
        do ke = 1,3
          Amplitude(ipl) = Amplitude(ipl) + conjg( u_s(ke) ) * sum( Pol_i(:) * RIXS_E1M1_tens_xtal(ke,:,1) ) &
                                          + conjg( Pol_s(ke) ) * sum( u_i(:) * RIXS_E1M1_tens_xtal(ke,:,2) )  
        end do
      endif
              
    end do  ! end of loop over polarization ipl

  endif
  
  return
end

!*****************************************************************************************************************

subroutine Write_rixs(Angle,Deltar,E_cut,E_loss,Energ_in,Epsii,Eseuil,File_name_rixs,Gamma_hole,Gamma_max, &
                  npl,jseuil,Mod_k,Monocrystal,n_q_dim,n_Spin_channel,n_theta,ne_loss,nenerg_in,ninit1,ninitr, &
                  nseuil,numat,RIXS_ampl)

  use declarations
  implicit none
 
  character(len=Length_name):: File_name_rixs
  
  integer:: i_q, i_Monocrystal, i_t, ie, ie_loss, ii, ipl, jseuil, n_isc, n_q_dim, &
            n_Spin_channel, n_theta, ne_loss, nenerg_in, ninit1, ninitr, nseuil, numat, npl  
   
  complex(kind=db), dimension(ne_loss,nenerg_in,npl,n_theta,n_q_dim,n_Spin_channel*ninitr):: RIXS_ampl

  logical:: Monocrystal
   
  real(kind=db):: Deltar, E_cut, E_in, E_out, Eseuil, Fac, Gamma_hole, Gamma_max 
  real(kind=db), dimension(ninitr):: Epsii
  real(kind=db), dimension(nenerg_in):: Energ_in
  real(kind=db), dimension(ne_loss):: E_loss
  real(kind=db), dimension(2,n_theta):: Angle
  real(kind=db), dimension(nenerg_in,ne_loss,n_theta,n_q_dim):: Mod_k

  do ie = 1,nenerg_in
    E_in = Eseuil + Energ_in(ie)
    do ie_loss = 1,ne_loss
      E_out = E_in - E_loss(ie_loss)
      Fac = sqrt( E_out / E_in )
      RIXS_ampl(ie_loss,ie,:,:,:,:) = Fac * RIXS_ampl(ie_loss,ie,:,:,:,:) 
    end do
  end do
  
  n_isc = n_Spin_channel * ninitr
    
  if( Monocrystal ) then
    i_Monocrystal = 1
  else
    i_Monocrystal = 0
  endif

  open(2, file = File_name_rixs )
  write(2,110) Eseuil*rydb, numat, nseuil, jseuil, i_Monocrystal, nenerg_in, ne_loss, npl, n_theta, n_q_dim, ninit1, ninitr, &
               n_Spin_channel, Epsii(:)*rydb
  write(2,115) Gamma_hole, Gamma_max, E_cut, Deltar 
  write(2,'(A)') ' Theta_in, 2.Theta'
  write(2,120) ( Angle(:,i_t) / radian, i_t = 1,n_theta )
  if( Monocrystal ) then
    write(2,'(A)') ' k'
    write(2,125) ((( Mod_k(ie,:,i_t,i_q), i_t = 1,n_theta), i_q = 1,n_q_dim ), ie = 1,nenerg_in )
  endif
  write(2,'(A)') ' Energy'
  write(2,120) Energ_in(:)*rydb

  do ie_loss = 1,ne_loss
    do ipl = 1,npl
      do i_t = 1,n_theta 
        do i_q = 1,n_q_dim
          do ii = 1,n_isc
            if( abs( E_loss(ie_loss)*Rydb ) < 9.999995_db ) then
              write(2,165) E_loss(ie_loss)*Rydb, RIXS_ampl(ie_loss,:,ipl,i_t,i_q,ii)
            elseif( abs( E_loss(ie_loss)*Rydb ) < 99.99995_db ) then
              write(2,175) E_loss(ie_loss)*Rydb, RIXS_ampl(ie_loss,:,ipl,i_t,i_q,ii)
            elseif( abs( E_loss(ie_loss)*Rydb ) < 999.9995_db ) then
              write(2,185) E_loss(ie_loss)*Rydb, RIXS_ampl(ie_loss,:,ipl,i_t,i_q,ii)
            else
              write(2,195) E_loss(ie_loss)*Rydb, RIXS_ampl(ie_loss,:,ipl,i_t,i_q,ii)
            endif 
          end do
        end do
      end do
    end do
  end do
  
  Close(2)

  return
  110 format(f10.3,i5,3i3,2i6,i3,i6,4i3,1p,10e15.7)
  115 format(10x,4f12.5,' = Gamma_hole, Gamma_max, E_cut, Deltar')
  120 format(10x,100000f15.5)
  125 format(10x,1p,100000e13.5)
  165 format(f10.5,1p,100000(1x,2e13.5))
  175 format(f10.4,1p,100000(1x,2e13.5))
  185 format(f10.3,1p,100000(1x,2e13.5))
  195 format(f10.2,1p,100000(1x,2e13.5))
end

!*****************************************************************************************************************

subroutine Write_Int_rixs(File_name,n_File,RIXS_core)

  use declarations
  implicit none
 
  integer:: i, i_File, i_q, i_Monocrystal, i_t, ie, ie_loss, i_Spin_channel, ii, initr, ipl, ipr, istat, je_loss, &
            jseuil, k, L, Length, &
            n_File, n_i, n_isc, n_Spin_channel, n_Spin_channel_out, n_theta, n_q_dim, ne_loss, nenerg_in, ninit1, ninitr, npl, &
            npl_out, nseuil, numat  
 
  character(len=3):: mot3
  character(len=8):: dat
  character(len=10):: mot10, tim
  character(len=13):: Word
  character(len=19):: com_date
  character(len=24):: com_time
  character(len=3), dimension(7):: suf_spin
  character(len=13), dimension(2:7):: Channel_name
  character(len=13), dimension(0:10):: init_name
  character(len=13), dimension(8):: pol_name
  character(len=13), dimension(8):: pol_suf
  character(len=Length_name):: File_name_int_sum, File_name_o, File_name_out, mot
  character(len=Length_name), dimension(n_File):: File_name
  character(len=13), dimension(:), allocatable:: E_string, E_string_full
  character(len=Length_name), dimension(:,:), allocatable:: File_name_int
    
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: RIXS_ampl
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: RIXS_ampl_tot

  logical:: E_cut_man, Gamma_hole_man, Gamma_max_man, Powder, RIXS_core
  
  real(kind=db):: E, Delta_E, Deltar, E_cut, Eseuil, Gamma_hole, Gamma_max, Range_E_in, Range_E_Loss, x, y 
  real(kind=db), dimension(:), allocatable:: Ampl_i, Ampl_r, E_loss, Energ_in
  real(kind=db), dimension(:,:), allocatable:: Angle
  real(kind=db), dimension(:,:,:,:), allocatable:: Mod_k, Sum_E_in 
  real(kind=db), dimension(:,:,:,:,:), allocatable:: RIXS_int
  real(kind=db), dimension(:,:,:,:,:,:), allocatable:: Sum_E_loss

  data suf_spin/ '   ', '_uu', '_ud', '_du', '_dd', '_yy', '_xx' /
  data channel_name/ '    up-up    ', '    up-down  ', '   down-up   ', '  down-down  ','  same-spin  ', 'cross-channel' /
  data init_name/ '       Total ', '      init_1 ', '      init_2 ', '      init_3 ', '      init_4 ', '      init_5 ', &
                  '      init_6 ','      init_7 ', '      init_8 ', '      init_9 ', '     init_10 '/
  data pol_name/ '  sigma-sigma', '    sigma-pi ', '    pi-sigma ', '     pi-pi   ','  right-sigma', '   left-sigma', &
                 '    right-pi ', '    left-pi  ' /
  data pol_suf/ '_ss', '_sp', '_ps', '_pp','_rs', '_ls', '_rp', '_lp' /

! Check of input file names 
  do i_File = 1,n_File
    open(2, file = File_name(i_File), status='old', iostat=istat)
    if( istat /= 0 ) then
      mot = File_name(i_File)
      Length = len_trim(mot)
      if( mot(Length-3:Length) /= '.txt' ) then
        mot(Length+1:Length+4) = '.txt'
        Close(2) 
        open(2, file = mot, status='old', iostat=istat)
        if( istat == 0 ) File_name(i_File) = mot
      endif  
      if( istat /= 0 ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,110) File_name(i_File)
        end do
        stop
      endif
    endif
    read(2,*) Eseuil, numat, nseuil, jseuil, i_Monocrystal, nenerg_in, ne_loss, npl, n_theta, n_q_dim, ninit1, ninitr, &
              n_Spin_channel
    Powder = i_Monocrystal == 0
    Close(2)
  end do
  
  if( Powder ) then
    npl_out = 1
  else
    npl_out = npl
  endif

  if( RIXS_core .and. ninitr > 1 ) then
    n_i = ninitr
  else
    n_i = 0
  endif

  n_isc = ninitr * n_Spin_channel
   
  allocate( Energ_in(nenerg_in) )
  if( nenerg_in > 1 ) then
    allocate( E_string(nenerg_in+1) )
    allocate( E_string_full((nenerg_in+1)*npl_out) )
  else
    allocate( E_string(nenerg_in) )
    allocate( E_string_full(nenerg_in*npl_out) )
  endif
  allocate( E_loss(ne_loss) )
  allocate( Ampl_i(nenerg_in) )
  allocate( Ampl_r(nenerg_in) )
  allocate( RIXS_ampl(nenerg_in,npl,n_theta,n_q_dim,n_isc) )
  allocate( Angle(2,n_theta) )
  if( .not. Powder ) allocate( Mod_k(nenerg_in,ne_loss,n_theta,n_q_dim) )

  if( n_Spin_channel == 1 ) then
    n_Spin_channel_out = n_Spin_channel 
  else
    n_Spin_channel_out = n_Spin_channel + 2 
  endif
  allocate( RIXS_ampl_tot(nenerg_in,npl,n_theta,n_q_dim,0:n_i,n_Spin_channel_out) )
  allocate( File_name_int(n_q_dim,n_Spin_channel_out) )
  
  open(2, file = File_name(1), status='old')
  read(2,*)  
  read(2,*) Gamma_hole, Gamma_max, E_cut, Deltar
  read(2,*) 
  read(2,*) ( Angle(:,i_t), i_t = 1,n_theta )
  if( .not. Powder ) then
    read(2,*) 
    read(2,*) ((( Mod_k(ie,:,i_t,i_q), i_t = 1,n_theta), i_q = 1,n_q_dim ), ie = 1,nenerg_in )
  endif
  read(2,*) 
  read(2,*) Energ_in(:)
  do ie_loss = 1,ne_loss
    do ipl = 1,npl
      do i_t = 1,n_theta 
        do i_q = 1,n_q_dim
          do ii = 1,n_isc
            read(2,*) E_loss(ie_loss)
          end do
        end do
      end do
    end do  
  end do
  Close(2)
  
  Gamma_hole_man = Gamma_hole > -1000._db
  Gamma_max_man = Gamma_max > -1000._db
  E_cut_man = E_cut > -1000._db

! File_name out
  call date_and_time( date = dat, time = tim ) 
  com_date = ', date = ' // dat(7:8) // ' ' // dat(5:6) // ' ' // dat(1:4)
  com_time = ', time = ' // tim(1:2)  // ' h ' // tim(3:4) // ' mn ' // tim(5:6) // ' s'
  
  File_name_o = File_name(1)
  Length = len_trim(File_name_o)
  if( File_name_o(Length-5:Length-3) == '_1') then
    File_name_o(Length-5:Length) = '      '
  else
    File_name_o(Length-3:Length) = '    '
  endif

  File_name_int_sum = File_name_o
   
  do i_Spin_channel = 1,n_Spin_channel_out
       
    do i_q = 1,n_q_dim

      File_name_out = File_name_o
      Length = len_trim(File_name_out)

      if( n_q_dim > 1 ) then
        File_name_out(Length+1:Length+3) = '_iq'
        call ad_number(i_q,File_name_out,Length_name)
        Length = len_trim(File_name_out)
      endif

      if( i_Spin_channel > 1 ) then
        File_name_out(Length+1:Length+3) = suf_spin( i_Spin_channel )
        Length = len_trim(File_name_out)
      endif

      do k = 1,4
      
        select case(k)
          case(1)
           if( .not. E_cut_man ) cycle
           y = E_cut
          case(2)
           if( .not. Gamma_max_man ) cycle
           y = Gamma_max
          case(3)
           if( .not. Gamma_hole_man ) cycle
           y = Gamma_hole
          case(4)
           if( abs( Deltar ) < eps10 ) cycle
           y = Deltar          
       end select
      
        x = y / 10
        do i = 0,4
          x = 10 * x
          if( abs( x - nint(x) ) < eps10 ) exit
        end do
      
        select case(i)
          case(0)
            write(mot10,'(i10)') nint( y )
          case(1)
            write(mot10,'(f10.1)') y
          case(2)
            write(mot10,'(f10.2)') y
          case(3)
            write(mot10,'(f10.3)') y
          case default
            write(mot10,'(f10.4)') y
        end select
      
        mot10 = adjustl( mot10 )   
        L = len_trim( mot10 )
        if( i > 0 ) mot10(L-i:L-i) = 'p'
     
        select case(k)
          case(1)          
            mot3 = '_EF'
          case(2)
            mot3 = '_Gm'
          case(3)
            mot3 = '_GH'
          case(4)
            mot3 = '_Ga'
        end select
      
        Length = len_trim( File_name_out )
        File_name_out(Length+1:Length+3) = mot3
        Length = Length + 3
        File_name_out(Length+1:Length+L) = mot10(1:L)

        if( i_spin_channel > 1 .or. i_q > 1 ) cycle
        
        Length = len_trim( File_name_int_sum )
        File_name_int_sum(Length+1:Length+3) = mot3
        Length = Length + 3
        File_name_int_sum(Length+1:Length+L) = mot10(1:L)
              
      end do

      Length = len_trim(File_name_out)
          
      File_name_out(Length+1:Length+8) = '_int.txt'
      File_name_int(i_q,i_spin_channel) = File_name_out 
      
    end do
  end do
  
  Length = len_trim(File_name_int_sum)
  File_name_int_sum(Length+1:Length+17) = '_int_sum_loss.txt'

  i = 0
  do ipl = 1,npl_out
    do ie = 1,nenerg_in+1
      if( nenerg_in == 1 .and. ie == 2 ) cycle
      i = i + 1
      Word = ' '
      if( ie == nenerg_in+1 ) then 
        Word = 'Av_on_Ei'
      else
        E = Energ_in(ie)
        if( abs( nint( 10*E ) - 10*E ) < eps10 ) then
          write(Word,'(f13.1)') E
        elseif( abs( nint( 100*E ) - 100*E ) < eps10 ) then
          write(Word,'(f13.2)') E
        elseif( abs( nint( 1000*E ) - 1000*E ) < eps10 ) then
          write(Word,'(f13.3)') E
        elseif( abs( nint( 1000*E ) - 1000*E ) < eps10 ) then
          write(Word,'(f13.4)') E
        endif
      endif

      Word = adjustr(Word)
      E_string(ie) = Word

      if( npl_out > 1 ) then
        Word = adjustl(Word)
        Length = len_trim(Word)
        if( Length+3 < 13 ) Word(Length+1:Length+3) = pol_suf(ipl)
      endif
      
      Word = adjustr(Word)
      E_string_full(i) = Word
    end do
  end do

  allocate( RIXS_int(nenerg_in,npl_out,n_theta,n_q_dim,0:n_i) )
  allocate( Sum_E_in(npl_out,n_theta,n_q_dim,0:n_i) )
  allocate( Sum_E_loss(nenerg_in,npl_out,n_theta,n_q_dim,0:n_i,n_spin_channel_out) )
  Sum_E_loss(:,:,:,:,:,:) = 0._db
  if( ne_loss > 1 ) &
    Range_E_Loss = E_loss(ne_loss) - E_loss(1) + ( E_loss(2) - E_loss(1) + E_loss(ne_loss) - E_loss(ne_loss-1) ) / 2
  if( nenerg_in > 1 ) &
    Range_E_in = Energ_in(nenerg_in) - Energ_in(1) + ( Energ_in(2) - Energ_in(1) + Energ_in(nenerg_in) - Energ_in(nenerg_in-1) ) / 2

  do ie_loss = 1,ne_loss

    RIXS_ampl_tot(:,:,:,:,:,:) = ( 0._db, 0._db )

    do i_File = 1,n_File
      open(2, file = File_name(i_File), status='old')
      do i = 1,6
        read(2,*)
      end do
      if( .not. Powder ) then
        read(2,*)
        read(2,*)
      endif
      do je_loss = 1,(ie_loss - 1) * npl * n_theta * n_q_dim * n_isc
        read(2,*) 
      end do
      do ipl = 1,npl
        do i_t = 1,n_theta 
          do i_q = 1,n_q_dim
            do ii = 1,n_isc
              read(2,*) E_loss(ie_loss), ( Ampl_r(ie), Ampl_i(ie), ie = 1,nenerg_in )
              RIXS_ampl(:,ipl,i_t,i_q,ii) = cmplx( Ampl_r(:), Ampl_i(:), db )
            end do
          end do
        end do
      end do  

      do i_Spin_channel = 1,n_Spin_channel
        do initr = 1,ninitr
          ii = initr + ( i_Spin_channel - 1 ) * ninitr

          if( n_i > 0 ) then
            RIXS_ampl_tot(:,:,:,:,initr,i_Spin_channel) = RIXS_ampl_tot(:,:,:,:,initr,i_Spin_channel) + RIXS_ampl(:,:,:,:,ii)
            if( i_Spin_channel == 2 .or. i_Spin_channel == 5 ) then
              RIXS_ampl_tot(:,:,:,:,initr,6) = RIXS_ampl_tot(:,:,:,:,initr,6) + RIXS_ampl(:,:,:,:,ii)
            elseif( i_Spin_channel == 3 .or. i_Spin_channel == 4 ) then
              RIXS_ampl_tot(:,:,:,:,initr,7) = RIXS_ampl_tot(:,:,:,:,initr,7) + RIXS_ampl(:,:,:,:,ii) 
            endif
          endif

          RIXS_ampl_tot(:,:,:,:,0,i_Spin_channel) = RIXS_ampl_tot(:,:,:,:,0,i_Spin_channel) + RIXS_ampl(:,:,:,:,ii) 
          if( i_Spin_channel == 2 .or. i_Spin_channel == 5 ) then
            RIXS_ampl_tot(:,:,:,:,0,6) = RIXS_ampl_tot(:,:,:,:,0,6) + RIXS_ampl(:,:,:,:,ii)
          elseif( i_Spin_channel == 3 .or. i_Spin_channel == 4 ) then
            RIXS_ampl_tot(:,:,:,:,0,7) = RIXS_ampl_tot(:,:,:,:,0,7) + RIXS_ampl(:,:,:,:,ii) 
          endif

        end do
      end do

      close(2)
    end do ! end of loop on files

    do i_Spin_channel = 1,n_Spin_channel_out
   
      if( Powder ) then
! summation on the tensors components (fake polarization directions)
        RIXS_int(:,:,:,:,:) = 0._db
        do ipl = 1,npl
          RIXS_int(:,1,:,:,:) = RIXS_int(:,1,:,:,:) + real( RIXS_ampl_tot(:,ipl,:,:,:,i_Spin_channel), db )**2 &
                                                       + aimag( RIXS_ampl_tot(:,ipl,:,:,:,i_Spin_channel) )**2  
        end do
      else
        RIXS_int(:,:,:,:,:) = real( RIXS_ampl_tot(:,:,:,:,:,i_Spin_channel), db )**2 &
                            + aimag( RIXS_ampl_tot(:,:,:,:,:,i_Spin_channel) )**2 
      endif

      if( nenerg_in > 1 ) then
        Sum_E_in(:,:,:,:) = 0._db
        do ie = 1,nenerg_in
          if( ie == 1 ) then
            Delta_E = Energ_in(2) - Energ_in(1) 
          elseif( ie == nenerg_in ) then
            Delta_E = Energ_in(nenerg_in) - Energ_in(nenerg_in-1) 
          else
            Delta_E = ( Energ_in(ie + 1) - Energ_in(ie - 1) ) / 2 
          endif 
          Sum_E_in(:,:,:,:) = Sum_E_in(:,:,:,:) + RIXS_int(ie,:,:,:,:) * Delta_E / Range_E_in
        end do
      endif

      if( ne_loss > 1 ) then
        if( ie_loss == 1 ) then
          Delta_E = E_loss(2) - E_loss(1) 
        elseif( ie_loss == ne_loss ) then
          Delta_E = E_loss(ne_loss) - E_loss(ne_loss-1) 
        else
          Delta_E = ( E_loss(ie_loss + 1) - E_loss(ie_loss - 1) ) / 2 
        endif 
        Sum_E_loss(:,:,:,:,:,i_Spin_channel) = Sum_E_loss(:,:,:,:,:,i_Spin_channel) + RIXS_int(:,:,:,:,:) * Delta_E / Range_E_Loss 
      endif
      
! Writing
    
      do i_q = 1,n_q_dim

        if( ie_loss == 1 ) then    

          open(2, file = File_name_int(i_q,i_spin_channel) )
        
          write(2,'(a15,a19,a24)') ' RIXS intensity', com_date, com_time
          write(2,120) nenerg_in, ne_loss, npl_out, n_theta, n_q_dim, n_i, n_Spin_channel
    
          if( i_Spin_channel > 1 ) write(2,130) '  Channel: ', Channel_name(i_spin_Channel) 
          if( n_q_dim > 1 ) write(2,130) '  q index =', i_q 

          if( n_i > 0 ) write(2,130) ' init:    ', (((( init_name(initr), ie = 1,nenerg_in+1 ), ipl = 1,npl_out ), &
                                                                                             i_t = 1,n_theta ), initr = 0, n_i ) 
          if( npl_out > 1 ) write(2,130) ' Polariza:', (((( pol_name(ipl), ie = 1,nenerg_in+1 ), ipl = 1,npl_out ), &
                                                                                             i_t = 1,n_theta ), initr = 0, n_i ) 
          write(2,140) ' Theta_in:',(((( Angle(1,i_t), ie = 1,nenerg_in+1), ipl = 1,npl_out ), i_t = 1,n_theta ), initr = 0, n_i )
          write(2,140) '  2.Theta:',(((( Angle(2,i_t), ie = 1,nenerg_in+1), ipl = 1,npl_out ), i_t = 1,n_theta ), initr = 0, n_i )
          write(2,130) '   E_loss ', ( E_string(:), ii = 1, npl_out*n_theta*( 1 + n_i ) )
          if( npl_out > 1 ) write(2,130) '   E_loss ', ( E_string_full(:), ii = 1, ( 1 + n_i )*n_theta )
         
        else
        
          open(2, file = File_name_int(i_q,i_spin_channel), position='append' )
        
        endif

        E = - E_loss(ie_loss) 
        if( nenerg_in > 1 ) then        
          if( abs( E ) < 9.999995_db ) then
            write(2,160) E, ((( RIXS_int(:,ipl,i_t,i_q,initr), Sum_E_in(ipl,i_t,i_q,initr), ipl = 1,npl_out ), &
                                                                                        i_t = 1,n_theta ), initr = 0, n_i ) 
          elseif( abs( E ) < 99.99995_db ) then
            write(2,170) E, ((( RIXS_int(:,ipl,i_t,i_q,initr), Sum_E_in(ipl,i_t,i_q,initr), ipl = 1,npl_out ), &
                                                                                        i_t = 1,n_theta ), initr = 0, n_i )
          elseif( abs( E ) < 999.9995_db ) then
            write(2,180) E, ((( RIXS_int(:,ipl,i_t,i_q,initr), Sum_E_in(ipl,i_t,i_q,initr), ipl = 1,npl_out ), &
                                                                                        i_t = 1,n_theta ), initr = 0, n_i )
          else
            write(2,190) E, ((( RIXS_int(:,ipl,i_t,i_q,initr), Sum_E_in(ipl,i_t,i_q,initr), ipl = 1,npl_out ), &
                                                                                        i_t = 1,n_theta ), initr = 0, n_i ) 
          endif
        else
          if( abs( E ) < 9.999995_db ) then
            write(2,160) E, ((( RIXS_int(:,ipl,i_t,i_q,initr), ipl = 1,npl_out ), i_t = 1,n_theta ), initr = 0, n_i ) 
          elseif( abs( E ) < 99.99995_db ) then
            write(2,170) E, ((( RIXS_int(:,ipl,i_t,i_q,initr), ipl = 1,npl_out ), i_t = 1,n_theta ), initr = 0, n_i )
          elseif( abs( E ) < 999.9995_db ) then
            write(2,180) E, ((( RIXS_int(:,ipl,i_t,i_q,initr), ipl = 1,npl_out ), i_t = 1,n_theta ), initr = 0, n_i )
          else
            write(2,190) E, ((( RIXS_int(:,ipl,i_t,i_q,initr), ipl = 1,npl_out ), i_t = 1,n_theta ), initr = 0, n_i ) 
          endif
        endif 

        Close(2)

      end do ! end of loop on i_q
    end do
    
  end do ! end of loop on ie_loss

! Writing of the average over energy loss
  if( ne_loss > 1 ) then
    open(2, file = File_name_int_sum)
    write(2,'(a38,a19,a24)') ' RIXS, average along loss energy range', com_date, com_time
    if( npl_out == 1 .and. n_theta == 1 .and. n_i == 0 .and. n_q_dim == 1 .and. n_spin_channel_out == 1 ) then
      write(2,'(A)') '  Energ_in Average_RIXS'
    elseif( npl_out == 1 .and. n_theta == 1 .and. n_q_dim == 1 .and. n_spin_channel_out == 1 ) then
      write(2,'(a10,10(a12,i1))') ' Energ_in', ('   Av_RIXS_i',initr, initr = 0,n_i )
    elseif( n_theta == 1 .and. n_i == 0 .and. n_q_dim == 1 .and. n_spin_channel_out == 1 ) then
      write(2,'(a10,10(a10,a3))') ' Energ_in', ('  Av_RIXS',pol_suf(ipl), ipl = 1,npl_out )
    elseif( npl_out == 1 .and. n_i == 0 .and. n_q_dim == 1 .and. n_spin_channel_out == 1 ) then
      write(2,'(a10,10(a12,i1))') ' Energ_in', ('   Av_RIXS_t',i_t, i_t = 1,n_theta )
    elseif( npl_out == 1 .and. n_theta == 1 .and. n_i == 0 .and. n_spin_channel_out == 1 ) then
      write(2,'(a10,10(a12,i1))') ' Energ_in', ('   Av_RIXS_q',i_q, i_q = 1,n_q_dim )
    elseif( npl_out == 1 .and. n_theta == 1 .and. n_i == 0 .and. n_q_dim == 1 ) then
      write(2,'(a10,10(a10,a3))') ' Energ_in', ('   Av_RIXS', suf_spin(i_spin_channel), i_spin_channel = 1,n_spin_channel_out )
    else
      write(2,'(a10,10(a4,a3,3(a1,i1)))') ' Energ_in', (((('  AR',pol_suf(ipl),'_t',i_t,'i',i_t,'q',i_q, ipl = 1,npl_out ), &
                                                           i_t = 1,n_theta ), initr = 0, n_i ),  i_q = 1,n_q_dim )
    endif     
    do ie = 1,nenerg_in
      if( abs( Energ_in(ie) ) < 9.999995_db ) then
        write(2,160) Energ_in(ie), ((((( Sum_E_loss(ie,ipl,i_t,i_q,initr,i_spin_channel), ipl = 1,npl_out ), i_t = 1,n_theta ), &
                                                  initr = 0, n_i ), i_q = 1,n_q_dim ), i_spin_channel = 1,n_spin_channel_out )
      elseif( abs( Energ_in(ie) ) < 99.99995_db ) then
        write(2,170) Energ_in(ie), ((((( Sum_E_loss(ie,ipl,i_t,i_q,initr,i_spin_channel), ipl = 1,npl_out ), i_t = 1,n_theta ), &
                                                  initr = 0, n_i ), i_q = 1,n_q_dim ), i_spin_channel = 1,n_spin_channel_out )
      elseif( abs( Energ_in(ie) ) < 999.9995_db ) then
        write(2,180) Energ_in(ie), ((((( Sum_E_loss(ie,ipl,i_t,i_q,initr,i_spin_channel), ipl = 1,npl_out ), i_t = 1,n_theta ), &
                                                  initr = 0, n_i ), i_q = 1,n_q_dim ), i_spin_channel = 1,n_spin_channel_out )
      else
        write(2,190) Energ_in(ie), ((((( Sum_E_loss(ie,ipl,i_t,i_q,initr,i_spin_channel), ipl = 1,npl_out ), i_t = 1,n_theta ), &
                                                  initr = 0, n_i ), i_q = 1,n_q_dim ), i_spin_channel = 1,n_spin_channel_out )
      endif
    end do
    Close(2)
  endif

! Writing of Mod(Q)  
  if( .not. Powder ) then
  
    do i_q = 1,n_q_dim

      File_name_out = File_name_o
      Length = len_trim(File_name_o)
  
      if( n_q_dim > 1 ) then
        File_name_out(Length+1:Length+3) = '_iq'
        call ad_number(i_q,File_name_out,Length_name)
        Length = len_trim(File_name_out)
      endif

      File_name_out = File_name_o 
      Length = len_trim(File_name_out)

      File_name_out(Length+1:Length+9) = '_modQ.txt'
  
      open(2, file = File_name_out)

      write(2,'(a7,a19,a24)') ' Mod(Q)', com_date, com_time
      write(2,200) nenerg_in, ne_loss, n_theta, n_q_dim
      write(2,140) ' Theta_in:', (( Angle(1,i_t), ie = 1,nenerg_in), i_t = 1,n_theta )
      write(2,140) '  2.Theta:', (( Angle(2,i_t), ie = 1,nenerg_in), i_t = 1,n_theta )
      write(2,140) '   E_loss ', ( Energ_in(:), ii = 1, n_theta ) 
  
      do ie_loss = 1,ne_loss
        E = - E_loss(ie_loss) 
        if( abs( E ) < 9.999995_db ) then
          write(2,160) E, ( Mod_k(:,ie_loss,i_t,i_q), i_t = 1,n_theta ) 
        elseif( abs( E ) < 99.99995_db ) then
          write(2,170) E, ( Mod_k(:,ie_loss,i_t,i_q), i_t = 1,n_theta )
        elseif( abs( E ) < 999.9995_db ) then
          write(2,180) E, ( Mod_k(:,ie_loss,i_t,i_q), i_t = 1,n_theta )
        else
          write(2,190) E, ( Mod_k(:,ie_loss,i_t,i_q), i_t = 1,n_theta ) 
        endif 
      end do
  
      Close(2)

    end do
    
  endif
  
  deallocate( Ampl_i, Ampl_r, Angle, E_loss, E_string, E_string_full, Energ_in, File_name_int, RIXS_ampl, RIXS_ampl_tot, &
              RIXS_int, Sum_E_in, Sum_E_loss )
  if( .not. Powder ) deallocate( Mod_k )
  
  return
  110 format(//' Error opening the file:',//3x,A,//)  
  120 format(7i5,' / n_energy_in, n_energy_loss, npl, n_theta, n_q_dir, n_core-state_print, n_Spin_channel')
  130 format(a10,1000000a13)
  140 format(a10,1000000f13.5)
  150 format(a10,1p,1000000e13.5)
  160 format(f10.5,1p,100000e13.5)
  170 format(f10.4,1p,100000e13.5)
  180 format(f10.3,1p,100000e13.5)
  190 format(f10.2,1p,100000e13.5)
  200 format(4i5,' / n_energy_in, n_energy_loss, n_theta, n_q_dir')
end

