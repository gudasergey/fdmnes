! FDMNES subroutine

! Set of subroutines for RIXS

!***********************************************************************

subroutine main_RIXS(Ampl_abs,Classic_irreg,Coef_g,Core_resolved,dV0bdcF,E_cut,E_cut_imp,E_Fermi,E_cut_man,Eclie,Eneg,Energ_s, &
        Epsii,Epsii_ref,Epsii_ref_man,Eseuil,Estart,File_name_rixs,FDM_comp_m,Full_potential,Gamma_hole, &
        Gamma_hole_man,Hubb_a,Hubb_d,iabsorig, &
        icheck,igreq,iprabs_nonexc,isymeq,ip0,ip_max,is_g,jseuil,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb,msymdd, &
        msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi, &
        Multipole,n_atom_proto,n_atom_uc,n_multi_run,n_oo,n_q_rixs,n_rel,natomsym,nbseuil, &
        nenerg_s,neqm,ngamh,ninit1,ninitl,ninitr,nlm_pot,nlm_probe,nlmomax,nomfich, &
        nr,nrm,ns_dipmag,nseuil,nspin,nspino,nspinp,numat,Orthmati,Posn, &
        psii,q_rixs,r,Relativiste,Renorm,Rmtg,Rmtsd,Rot_atom_abs,Rot_int,Spin_channel,Spinorbite,V_intmax, &
        V_hubb,Vcato,Vhbdc,VxcbdcF,Vxcato,Ylm_comp)

  use declarations
  implicit none

  integer, parameter:: mpinodes = 1
  integer, parameter:: mpirank = 0
  integer, parameter:: mpirank0 = 0
  integer, parameter:: n_Ec = 2
  integer, parameter:: n_V = 1
  integer, parameter:: ndim2 = 1
  integer, parameter:: nenerg_tddft = 0
  integer, parameter:: nlmamax = 0

  real(kind=db), parameter:: quatre_mc2 = 4 * 2 / alfa_sf**2  ! in Rydberg
  
  integer:: iabsorig, i_q, i_Spin_channel, i_theta, ia, ich, icheck, ie, ie_computer, ie_cond, ie_loss, ii, ip0, &
     ip_max, initr, iprabs_nonexc, iseuil, isp, &
     je, js, jseuil, ke, lmax_pot, lmax_probe, lseuil, m_hubb, &
     n_atom_proto, n_atom_uc, n_isc, n_multi_run, n_oo, n_q_dim, n_q_rixs, n_rel, n_Spin_channel, n_theta, natomsym, nbseuil, &
     ne, ne_loss, nenerg, nenerg_cond, nenerg_s, neqm, ngamh, ninit1, ninitl, ninitr, ninitlv, nlm_p_fp, nlm_pot, &
     nlm_probe, nlmomax, npl, nr, nrm, &
     ns_dipmag, nseuil, nspin, nspino, nspinp, numat

  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3):: msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi
  integer, dimension(natomsym):: isymeq
  integer, dimension(ninitl,2):: m_g
  integer, dimension(ninitl):: is_g
  integer, dimension(0:n_atom_proto,neqm):: igreq

  character(len=Length_name):: File_name_rixs, nomfich

  complex(kind=db):: CFac, Ph, Ph_m
  complex(kind=db), dimension(nenerg_s,nlm_probe,nspinp,nlmomax,nspinp):: Ampl_abs
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull, Taull_Spin_channel
  complex(kind=db), dimension(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil):: rof0   ! empty
  complex(kind=db), dimension(3,3,ninitr,0:mpinodes-1):: secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitr,0:mpinodes-1):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitr,0:mpinodes-1):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,3,n_rel,ninitr,0:mpinodes-1):: secdd, secdd_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:mpinodes-1):: secoo, secoo_m
  complex(kind=db), dimension(3,3,n_rel,ninitr,0:natomsym):: secddia, secddia_m
  complex(kind=db), dimension(3,3,ninitr,0:natomsym):: secmdia, secmdia_m, secmmia, secmmia_m
  complex(kind=db), dimension(3,3,3,ninitr,0:natomsym):: secdqia, secdqia_m
  complex(kind=db), dimension(3,3,3,3,ninitr,0:natomsym):: secdoia, secdoia_m, secqqia, secqqia_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:natomsym):: secooia, secooia_m
  complex(kind=db), dimension(3,3):: RIXS_E1M1_tens_xtal, RIXS_E1M1_tens_xtal_m, RIXS_M1M1_tens_xtal
  complex(kind=db), dimension(3,3,n_rel):: RIXS_E1E1_tens_xtal
  complex(kind=db), dimension(3,3,3):: RIXS_E1E2_tens_xtal, RIXS_E1E2_tens_xtal_m
  complex(kind=db), dimension(3,3,3,3):: RIXS_E1E3_tens_xtal, RIXS_E2E2_tens_xtal, RIXS_E3E3_tens_xtal

  complex(kind=db), dimension(:), allocatable:: Amplitude
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: RIXS_E1M1_tens, RIXS_E1M1_tens_m, RIXS_M1M1_tens
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: RIXS_ampl, RIXS_E1E1_tens, RIXS_E1E2_tens, RIXS_E1E2_tens_m
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: RIXS_E1E3_tens, RIXS_E2E2_tens, RIXS_E3E3_tens

  logical:: Circular, Classic_irreg, Core_resolved, Dip_rel, E_cut_man, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Epsii_ref_man, &
            Extract_ten, Eneg, FDM_comp_m, Final_optic, Final_tddft, Full_potential, Gamma_hole_man, &
            Green, Green_int, Hubb_a, Hubb_d, lmoins1, lplus1, M1M1, Magn_sens, NRIXS, Q_resolved, Relativiste, Renorm, RIXS, &
            Solsing, Solsing_only, Spin_channel, Spinorbite, Tddft, Ylm_comp 

  logical, dimension(10):: Multipole

  real(kind=db):: conv_mbarn_nelec, Cos_2t, cst, dTheta, E_cut, E_cut_imp, E_cut_rixs, E_Fermi, E_i, &
                  E_s, Eclie, Eimag, Emin, Ephoton, Epsii_ref, Estart, Mod_Q, Omega_i, Omega_s, Pas, Rmtg, Rmtsd, &
                  Tab_width, Theta, V_intmax, Vhbdc
  real(kind=db), dimension(3):: P, Pol_i_s, Pol_i_p, Pol_s_s, Pol_s_p, Vec_i, Vec_s, Wx_lab, Wy_lab, Wz_lab
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitr):: Epsii
  real(kind=db), dimension(ninitl,2):: coef_g
  real(kind=db), dimension(nspin):: dV0bdcF, VxcbdcF
  real(kind=db), dimension(nspin,n_V):: V0bdc
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(n_Ec):: Enervide
  real(kind=db), dimension(3,n_atom_uc):: posn
  real(kind=db), dimension(3,n_q_rixs):: q_rixs
  real(kind=db), dimension(3,3):: Mat_rot, Orthmatt, Orthmati, Rot_atom_abs, Rot_int
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nspin,n_Ec):: Ecinetic
  real(kind=db), dimension(nr,nlm_pot):: Vcato
  real(kind=db), dimension(nr,nlm_pot,nspin,n_V):: Vrato
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vxcato

  real(kind=db), dimension(:), allocatable:: Delta_E, E_loss, Energ, Energ_cond, Gamma, Shift

! When no momentum transfer direction given, calculation is done for powder 
  Q_resolved = n_q_rixs /= 0

  Eimag = 0._db
  Extract_ten = .false.
  Final_optic = .false.
  Final_tddft = .false.
  Green = .false.
  Green_int = .false.
  ie_computer = 0
  ninitlv = n_Ec
  NRIXS = .false.
  RIXS = .true.
  Solsing = .false.
  Solsing_only = .false.
  Tddft = .false.

  ich = max(icheck-1,1)
  Magn_sens = nspin == 2 .or. Spinorbite
  Circular = Magn_sens
  
  Dip_rel = n_rel > 1
  
  if( .not. Q_resolved ) then
    npl = 1
  elseif( Circular ) then
    npl = 8   ! number of polarization
  else
    npl = 4
  endif
  allocate( Amplitude(npl) )

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  if( Spin_channel ) then
    n_Spin_channel = 5
  else
    n_Spin_channel = 1
  endif
  
  n_isc = ninitr * n_Spin_channel

  allocate( Shift(ninitr) )
  call Cal_shift(Epsii,Epsii_ref,Epsii_ref_man,icheck,ninit1,ninitr,Shift)

  lmax_probe = nint( sqrt( real( nlm_probe, db ) ) ) - 1
  lmax_pot = nint( sqrt( real( nlm_pot, db ) ) ) - 1

  if( E_cut_man ) then
    E_cut_rixs = E_cut_imp
  else
    E_cut_rixs = E_cut
  endif

  if( nenerg_s > 1 ) then      
    Emin = Energ_s(1) - ( Energ_s(2) - Energ_s(1) ) / 2
  else
    Emin = Energ_s(1) - eps10
  endif
  do ie = 1,nenerg_s
    if( Energ_s(ie) > E_cut_RIXS + eps10 ) exit
  end do
  nenerg_cond = nenerg_s - ie + 1 
  ne_loss = nenerg_s - 1
  
  if( Estart < E_cut_RIXS .and. ie > 1 ) then
    Pas =  Energ_s(ie) - Energ_s(ie-1)
    ne = nint( ( Energ_s(ie) - Estart ) / pas )
  else
    ne = 0
  endif
  nenerg = nenerg_cond + ne
  
  allocate( Energ_cond(nenerg_cond) )
  allocate( Energ(nenerg) )
  allocate( Gamma(nenerg) )
  allocate( Delta_E(nenerg_cond) )
  allocate( E_loss(ne_loss) )

  if( E1E1 ) then
    allocate( RIXS_E1E1_tens(3,3,n_rel,n_isc,nenerg*ne_loss,natomsym) )
    RIXS_E1E1_tens(:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E1E2 ) then
    allocate( RIXS_E1E2_tens(3,3,3,n_isc,nenerg*ne_loss,natomsym) )
    RIXS_E1E2_tens(:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E1E2 .and. Magn_sens ) then
    allocate( RIXS_E1E2_tens_m(3,3,3,n_isc,nenerg*ne_loss,natomsym) )
    RIXS_E1E2_tens_m(:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E1E3 ) then
    allocate( RIXS_E1E3_tens(3,n_oo,3,n_oo,n_isc,nenerg*ne_loss,natomsym) )
    RIXS_E1E3_tens(:,:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E1M1 ) then
    allocate( RIXS_E1M1_tens(3,3,n_isc,nenerg*ne_loss,natomsym) )
    RIXS_E1M1_tens(:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E1M1 .and. Magn_sens ) then
    allocate( RIXS_E1M1_tens_m(3,3,n_isc,nenerg*ne_loss,natomsym) )
    RIXS_E1M1_tens_m(:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E2E2 ) then
    allocate( RIXS_E2E2_tens(3,3,3,3,n_isc,nenerg*ne_loss,natomsym) )
    RIXS_E2E2_tens(:,:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( E3E3 ) then
    allocate( RIXS_E3E3_tens(3,n_oo,3,n_oo,n_isc,nenerg*ne_loss,natomsym) )
    RIXS_E1E3_tens(:,:,:,:,:,:,:) = (0._db, 0._db)
  endif
  if( M1M1 ) then
    allocate( RIXS_M1M1_tens(3,3,n_isc,nenerg*ne_loss,natomsym) )
    RIXS_M1M1_tens(:,:,:,:,:) = (0._db, 0._db)
  endif

  if( Q_resolved ) then
    n_theta = 1
    dtheta = pi / n_theta 
  else
    n_theta = 1
  endif
  n_q_dim = max( n_q_rixs, 1 )
  
  allocate( RIXS_ampl(ne_loss,nenerg,npl,n_theta,n_q_dim,n_isc) )
  RIXS_ampl(:,:,:,:,:,:) = (0._db, 0._db) 

  call Energy_Grid_RIXS(E_loss,Energ,Energ_cond,Energ_s,icheck,ne_loss,nenerg,nenerg_cond,nenerg_s)

! One takes the exchange correlation potential independant from energy
  do isp = 1,nspinp
    Vrato(1:nr,:,isp,1) = Vcato(1:nr,:) + Vxcato(1:nr,:,isp)
    V0bdc(isp,1) = Vhbdc + VxcbdcF(isp) + dV0bdcF(isp)
  end do

  if( Hubb_a .or. Full_potential ) then
    nlm_p_fp = nlm_probe
  else
    nlm_p_fp = 1
  endif

  do ie_cond = 1,nenerg_cond
    if( nenerg_cond == 1 ) then
      Delta_E(ie_cond) = 0.1_db / Rydb
    elseif( ie_cond == nenerg_cond ) then
      Delta_E(ie_cond) = ( Energ_cond(ie_cond) - Energ_cond(ie_cond-1) ) / 2
    elseif( ie_cond == 1 ) then
      Delta_E(ie_cond) = ( Energ_cond(2) + Energ_cond(1) ) / 2 - E_cut_Rixs 
    else
      Delta_E(ie_cond) = ( Energ_cond(ie_cond+1) - Energ_cond(ie_cond-1) ) / 2
    endif
  end do

  if( Q_resolved ) then
    call invermat(Orthmatt,Orthmati)
    call Reciprocal_Vector(icheck,Orthmatt,Wx_lab,Wy_lab,Wz_lab)
  endif

  do i_q = 1,n_q_dim

    if( Q_resolved ) call Cal_Mat_rot(i_q,icheck,Mat_rot,n_q_rixs,Orthmati,Wx_lab,Wy_lab,Wz_lab)      

    do i_theta = 1,n_theta

      if( Q_resolved ) then
        Theta = ( i_theta - 1 ) * dtheta
        call Cal_Pol(icheck,Mat_rot,Pol_i_s,Pol_i_p,Pol_s_s,Pol_s_p,Theta,Vec_i,Vec_s)
      endif     
      
      do ie_cond = 1,nenerg_cond ! Loop over electron above Fermi level energy

        do ie_loss = 1,ne_loss
          ke = ie_loss + ( ie - 1 ) * ne_loss

          if( Energ_cond(ie_cond) - E_loss(ie_loss) < Emin ) exit
          if( Energ_cond(ie_cond) - E_loss(ie_loss) > E_cut_RIXS + eps10 ) cycle
         
          call Cal_Tau_loss(Ampl_abs,E_loss(ie_loss),Energ_cond(ie_cond),Energ_s,icheck,ndim2,nenerg_s,nlm_probe,nlmomax, &
                        ns_dipmag,nspin,nspino,nspinp,Spinorbite,Taull_Spin_channel)

          Enervide(1) = Energ_cond(ie_cond) + E_Fermi
          Enervide(2) = Enervide(1) - E_loss(ie_loss)
          do je = 1,n_Ec 
            Ecinetic(:,je) = Enervide(je) - V0bdc(:,1)
            if( .not. Eneg ) Ecinetic(:,je) = max( Ecinetic(:,je), Eclie )
          end do

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
       
            call Tenseur_car(Classic_irreg,coef_g,Core_resolved,Ecinetic, &
                    Eimag,Energ_cond(ie_cond),Enervide,Eseuil,FDM_comp_m,Final_optic,Final_tddft,Full_potential,Green,Green_int, &
                    Hubb_a,Hubb_d,ich,ie_cond,ip_max,ip0,is_g,lmax_probe,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                    mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi, &
                    Multipole, &
                    n_Ec,n_oo,n_rel,n_V,nbseuil,ns_dipmag,ndim2,nenerg_tddft,ninit1,ninitl,ninitr,ninitlv,nlm_pot,nlm_probe, &
                    nlm_p_fp,nlmamax,nr,nrm,nspin,nspino,nspinp,numat,psii,r,Relativiste,Renorm,RIXS,Rmtg,Rmtsd,rof0,rot_atom_abs, &
                    Rot_int,secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m, &
                    secqq,secqq_m,Solsing,Solsing_only,Spinorbite,Taull,Tddft,V_hubb,V_intmax,V0bdc,Vrato,Ylm_comp)
       
            if( icheck > 1 ) Call write_tensor_RIXS(i_Spin_channel,ie_cond,ie_loss,E_loss(ie_loss),E1E1,E2E2,Energ_cond(ie_cond), &
                                                    mpinodes,n_rel,ninitr,secdd,secqq)
                                                           
            call Tensor_equiv(Green_int,ie_computer, &
                isymeq,mpinodes,Multipole,n_oo,n_rel,natomsym,ninitr,nspin, &
                secdd,secdd_m,secddia,secddia_m,secdq,secdq_m,secdqia,secdqia_m,secdo,secdo_m,secdoia,secdoia_m, &
                secmd,secmd_m,secmdia,secmdia_m,secmm,secmm_m,secmmia,secmmia_m,secoo,secoo_m,secooia,secooia_m, &
                secqq,secqq_m,secqqia,secqqia_m,Spinorbite)

            do initr = 1,ninitr       ! ----------> Loop over edges or core states

              ii = initr + ( i_Spin_channel - 1 ) * ninitr
 
              if( Core_resolved ) then
                if( initr <= ninit1 ) then
                  iseuil = 1
                else
                  iseuil = min(2, nbseuil)
                endif
              else
                 iseuil = min(initr, nbseuil)
              endif
      
              if( Gamma_hole_man ) then
                if( ngamh == 1 ) then
                  Gamma(:) = Gamma_hole(1)
                elseif( ngamh == ninitr ) then
                  Gamma(:) = Gamma_hole(initr)
                elseif( initr <= ninit1 ) then
                  Gamma(:) = Gamma_hole(1)
                else
                  Gamma(:) = Gamma_hole(2)
                endif
              else
                js = jseuil + iseuil - 1 
                Gamma_hole(1) = Tab_width(Eseuil(iseuil),js,nbseuil,numat)
                Gamma(:) = Gamma_hole(1)
              endif
          
              do ie = 1,nenerg  ! loop on incoming energy

                ke = ie_loss + ( ie - 1 ) * ne_loss
                CFac = Delta_E(ie_cond) / (  Energ(ie) + Shift(initr) - Energ_cond(ie_cond) + img * Gamma(ie) )

                Ephoton = Eseuil(iseuil) + Energ(ie) + Shift(initr)
! For very low energy edges:
                Ephoton = max(0.001_db/Rydb, Ephoton)

! alfa_sf = e*e/(2*epsilon0*h*c) is the fine structure constant
                cst = quatre_pi * pi * alfa_sf * Ephoton
! To get result in Mbarn (10E-18 cm2)
                cst = 100 * bohr**2 * cst
    
                CFac = CFac * cst * conv_mbarn_nelec(Ephoton) / pi
                 
                if( Q_resolved ) then
                  E_i = Eseuil(iseuil) + Energ(ie) + Shift(initr)
                  E_s = E_i - E_loss(ie_loss)
                  Mod_Q = 0.5_db * alfa_sf * sqrt( E_s**2 + E_i**2 - 2 * E_s * E_i * cos_2t )
                  Omega_i = E_i / quatre_mc2 
                  Omega_s = E_s / quatre_mc2
                endif 
 
                if( E1E1 ) RIXS_E1E1_tens_xtal(:,:,:) = ( 0._db, 0._db )
                if( E1E2 ) RIXS_E1E2_tens_xtal(:,:,:) = ( 0._db, 0._db )
                if( E1E2 .and. Magn_sens ) RIXS_E1E2_tens_xtal_m(:,:,:) = ( 0._db, 0._db )
                if( E1M1 ) RIXS_E1M1_tens_xtal(:,:) = ( 0._db, 0._db )
                if( E1M1 .and. Magn_sens ) RIXS_E1M1_tens_xtal_m(:,:) = ( 0._db, 0._db )
                if( E2E2 ) RIXS_E2E2_tens_xtal(:,:,:,:) = ( 0._db, 0._db )
                if( E1E3 ) RIXS_E1E3_tens_xtal(:,:,:,:) = ( 0._db, 0._db )
                if( E3E3 ) RIXS_E3E3_tens_xtal(:,:,:,:) = ( 0._db, 0._db )
                if( M1M1 ) RIXS_M1M1_tens_xtal(:,:) = ( 0._db, 0._db )

! Calculation of crystal tensor
                do ia = 1,natomsym
                  if( Q_resolved ) then
                    P(:) = Posn(:,igreq(iprabs_nonexc,ia))
                    Ph = exp( img * Mod_Q * sum( Q_rixs(:,i_q) * P(:) ) )
                    Ph_m = img * Ph
                  else
                    Ph = ( 1._db, 0._db )
                  endif               
                  Ph_m = img * Ph

                  if( E1E1 ) RIXS_E1E1_tens_xtal(:,:,:) = RIXS_E1E1_tens_xtal(:,:,:) + Ph * CFac * Secddia(:,:,:,initr,ia)
                  if( E1E2 ) RIXS_E1E2_tens_xtal(:,:,:) = RIXS_E1E2_tens_xtal(:,:,:) &
                                                        + Ph * CFac * real( Secdqia(:,:,:,initr,ia), db )
                  if( E1E2 .and. Magn_sens ) RIXS_E1E2_tens_xtal_m(:,:,:) = RIXS_E1E2_tens_xtal_m(:,:,:) &
                                                                          + Ph_m * CFac * aimag( Secdqia(:,:,:,initr,ia) )
                  if( E1M1 ) RIXS_E1M1_tens_xtal(:,:) = RIXS_E1M1_tens_xtal(:,:) + Ph_m * CFac * aimag( Secmdia(:,:,initr,ia) )
                  if( E1M1 .and. Magn_sens ) RIXS_E1M1_tens_xtal_m(:,:) = RIXS_E1M1_tens_xtal_m(:,:) &
                                                                        + Ph * CFac * real( Secmdia(:,:,initr,ia), db )
                  if( E2E2 ) RIXS_E2E2_tens_xtal(:,:,:,:) = RIXS_E2E2_tens_xtal(:,:,:,:) + Ph * CFac * Secqqia(:,:,:,:,initr,ia)
                  if( E1E3 ) RIXS_E1E3_tens_xtal(:,:,:,:) = RIXS_E1E3_tens_xtal(:,:,:,:) + Ph * CFac * Secdoia(:,:,:,:,initr,ia)
                  if( E3E3 ) RIXS_E3E3_tens_xtal(:,:,:,:) = RIXS_E3E3_tens_xtal(:,:,:,:) + Ph * CFac * Secooia(:,:,:,:,initr,ia)
                  if( M1M1 ) RIXS_M1M1_tens_xtal(:,:) = RIXS_M1M1_tens_xtal(:,:) + Ph * CFac * Secmmia(:,:,initr,ia)
                end do

! Calculation of RIXS amplitude
                call Cal_RIXS_ampl(Amplitude,Dip_rel,E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,M1M1,Magn_sens, &
                         n_rel,npl,Omega_i,Omega_s,pol_i_s,pol_i_p,pol_s_s,pol_s_p,Q_resolved, &
                         RIXS_E1E1_tens_xtal,RIXS_E1E2_tens_xtal,RIXS_E1E2_tens_xtal_m,RIXS_E1E3_tens_xtal,&
                         RIXS_E2E2_tens_xtal,RIXS_E3E3_tens_xtal,RIXS_E1M1_tens_xtal,RIXS_E1M1_tens_xtal_m,RIXS_M1M1_tens_xtal, &
                         Vec_i,Vec_s)
                RIXS_ampl(ie_loss,ie,:,i_theta,i_q,ii) = RIXS_ampl(ie_loss,ie,:,i_theta,i_q,ii) + Amplitude(:)   
                
              end do ! end of loop over incoming energy ie
              
            end do  ! end of loop over initr      
 
          end do ! end of loop over i_Spin_channel

        end do ! end of loop over ie_loss

        if( E1E1 ) write(6,'(f10.5,1p,1x,2e13.5)') Energ_cond(ie_cond) * Rydb 

      end do  ! end of loop over ie_cond
      
    end do ! end of loop over i_theta
  end do ! end of loop ovzer i_q 

! Writing
  call Write_rixs(E_loss,Energ,Epsii,Eseuil(nbseuil),iabsorig,File_name_rixs,jseuil,n_multi_run,n_q_dim,n_Spin_channel,n_theta, &
                      ne_loss,nenerg,ninit1,ninitr,nomfich,npl,nseuil,numat,RIXS_ampl)
  
  deallocate( Amplitude, Delta_E, Energ, Energ_cond, E_loss, Gamma, Shift )
  deallocate( RIXS_ampl )
  
  if( E1E1 ) deallocate( RIXS_E1E1_tens )
  if( E1E2 ) deallocate( RIXS_E1E2_tens )
  if( E1E2 .and. Magn_sens ) deallocate( RIXS_E1E2_tens_m )
  if( E1E3 ) deallocate( RIXS_E1E3_tens )
  if( E1M1 ) deallocate( RIXS_E1M1_tens )
  if( E1M1 .and. Magn_sens ) deallocate( RIXS_E1M1_tens_m )
  if( E2E2 ) deallocate( RIXS_E2E2_tens )
  if( E3E3 ) deallocate( RIXS_E3E3_tens )
  if( M1M1 ) deallocate( RIXS_M1M1_tens )
  
  return
end

!***********************************************************************

! Energy grid for RIXS

subroutine Energy_Grid_RIXS(E_loss,Energ,Energ_cond,Energ_s,icheck,ne_loss,nenerg,nenerg_cond,nenerg_s)

  use declarations
  implicit none

  integer,intent(in):: icheck, ne_loss, nenerg, nenerg_cond, nenerg_s
  integer:: ie, je

  real(kind=db), dimension(nenerg_s), intent(in):: Energ_s

  real(kind=db), dimension(nenerg), intent(out):: Energ
  real(kind=db), dimension(nenerg_cond), intent(out):: Energ_cond
  real(kind=db), dimension(ne_loss), intent(out):: E_loss

  do ie = 1,nenerg_cond
    Energ_cond(ie) = Energ_s(ie + nenerg_s - nenerg_cond)
  end do
  do ie = 1,nenerg
    Energ(ie) = Energ_s(ie + nenerg_s - nenerg)
  end do

  je = 0
  do ie = nenerg_s - nenerg, 1, -1
    je = je + 1
    E_loss(je) = Energ(1) - Energ_s(ie)
  end do
  do ie = 2,nenerg
    je = je + 1
    E_loss(je) = Energ(ie) - Energ_s(1)
  end do

  if( icheck > 2 ) then
    write(3,100)
    write(3,'(/A)') '  The incoming energy grid for the RIXS calculation:'
    write(3,120) Energ(:)*rydb
    write(3,'(/A)') '  The energy loss grid for the RIXS calculation:'
    write(3,120) E_loss(:)*rydb
  end if

  return
  100 format(/' ---- Energy_Grid_RIXS -',100('-'))
  110 format()
  120 format(5f13.7)
end

!***************************************************************

! Calculation of the shifts

subroutine Cal_shift(Epsii,Epsii_ref,Epsii_ref_man,icheck,ninit1,ninitr,Shift)

  use declarations
  implicit none

  integer:: icheck, ninit1, ninitr
  
  logical:: Epsii_ref_man
  
  real(kind=db):: Epsii_moy, Epsii_ref
  real(kind=db), dimension(ninitr):: Epsii, Shift
   
  if( .not. Epsii_ref_man ) then
    select case(ninitr)
      case(1)
        Epsii_moy = Epsii(1)
      case(2)
        if( ninit1 == 2 ) then
          Epsii_moy = sum( Epsii(1:ninit1) ) / ninit1
        else
          Epsii_moy = Epsii(1)
        endif
      case(4,6,10)
         if( ninit1 /= ninitr ) then
           Epsii_moy = sum( Epsii(ninit1+1:ninitr) ) / ( ninitr - ninit1 )
         else
           Epsii_moy = sum( Epsii(1:ninitr) ) / ninitr
         endif
    end select
    Epsii_ref = Epsii_moy
    Epsii_ref_man = .true.
  endif
  Shift(:) = Epsii(:) - Epsii_ref 

  if( icheck > 1 ) then
    write(3,'(/A)') ' Shift(initl = 1,ninitl)'
    write(3,'(10f9.5)') Shift(:) * rydb
  endif  
    
  return
end

!***************************************************************

! Calculation of the reciprocal vector in the lab frame

subroutine Reciprocal_Vector(icheck,Orthmatt,Wx_lab,Wy_lab,Wz_lab)

  use declarations
  implicit none

  integer:: i, icheck
  
  real(kind=db):: Cos_a, Cos_b, Cos_t, Phi, Sin_t, Theta, Vol, Wx_mod, Wy_mod, Wz_mod 
  real(kind=db), dimension(3):: Vec, Vx, Vy, Vz, Wx, Wy, Wx_lab, Wy_lab, Wz, Wz_lab
  real(kind=db), dimension(3,3):: Orthmatt
  
! Diffraction vector in the internal orthogonal basis
  Vx(:) = orthmatt(:,1)
  Vy(:) = orthmatt(:,2)
  Vz(:) = orthmatt(:,3)

! wx, wy, wz : reciprocal cell basis
  call prodvec(Wx,vy,vz)

  vol = sum( Wx(:) * vx(:) )
  Wx(:) = Wx(:) / vol
  call prodvec(Wy,vz,vx)
  Wy(:) = Wy(:) / vol
  call prodvec(Wz,vx,vy)
  Wz(:) = Wz(:) / vol

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
    Wz_lab(1) = 0._db; Wz_lab(2) = 0._db; Wz_lab(3) = 1._db
  elseif(  abs( Cos_t + 1._db ) < eps10 ) then
    Wz_lab(1) = 0._db; Wz_lab(2) = 0._db; Wz_lab(3) = -1._db
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

  if( icheck > 1 ) then
    write(3,'(/A)') ' Reciprocal vector in lab frame:'
    write(3,'(A)') '      A         B        C'
    do i = 1,3
      write(3,'(3f10.5)') Wx_lab(i), Wy_lab(i), Wz_lab(i)
    end do
  endif 

  return
end

!***************************************************************

! Calculation of the inelastic scattering amplitude

subroutine Cal_Mat_rot(i_q,icheck,Mat_rot,n_q_rixs,Orthmati,Wx_lab,Wy_lab,Wz_lab)

  use declarations
  implicit none

  integer:: i, i_q, icheck, n_q_rixs
  
  real(kind=db):: Cos_p, Cos_t, Phi, Q_mod, Sin_p, Sin_t, Theta 
  real(kind=db), dimension(3):: Q_lab, Wx_lab, Wy_lab, Wz_lab
  real(kind=db), dimension(3,3):: Mat_rot, Orthmati, Rot, Rot_1, Rot_2
  real(kind=db), dimension(3,n_q_rixs):: Q_rixs
  
  Q_lab(:) = Q_rixs(1,i_q) * Wx_lab(:) + Q_rixs(2,i_q) * Wy_lab(:) + Q_rixs(3,i_q) * Wz_lab(:)
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
  
  Cos_t = Q_lab(1)
  Theta = - pi/2 + acos( Q_lab(1) )   ! rotation is negative
  Sin_t = sin( Theta )
  Rot_2(1,1) = Cos_t; Rot_2(1,2) = 0._db; Rot_2(1,3) = Sin_t  
  Rot_2(2,1) = 0._db; Rot_2(2,2) = 1._db; Rot_2(2,3) = 0._db  
  Rot_2(3,1) = - Sin_t; Rot_2(3,2) = 0._db; Rot_2(3,3) = Cos_t
  
  Rot = Matmul( Rot_2, Rot_1 ) 
  Rot = Transpose( Rot )  ! it is to transform vec and pol in lab frame to crystal frame 

  Mat_rot = Matmul( Orthmati, Rot )

  if( icheck > 1 ) then
    write(3,110) i_q, Q_rixs(:,i_q)
    write(3,'(/A)') ' Rotation matrix'
    do i = 1,3
      write(3,'(3f10.5)') Mat_rot(i,:)
    end do 
  endif
  
  return
  110 format(' Q direction number,',i2,' Q =',3f10.5)
end

!***************************************************************

! Calculation of the inelastic scattering amplitude

subroutine Cal_Pol(icheck,Mat_rot,Pol_i_s,Pol_i_p,Pol_s_s,Pol_s_p,Theta,Vec_i,Vec_s)

  use declarations
  implicit none

  integer:: i, icheck
  
  real(kind=db):: Cos_2t, Cos_t, Sin_t, Theta 
  real(kind=db), dimension(3):: Pol_i_s, Pol_i_p, Pol_s_s, Pol_s_p, Vec_i, Vec_s
  real(kind=db), dimension(3,3):: Mat_rot, Rot
  
  Cos_t = cos( Theta )
  Sin_t = sin( Theta )  
  Cos_2t = cos( 2 * Theta )

  Rot(1,1) = 1._db; Rot(1,2) = 0._db; Rot(1,3) = 0._db
  Rot(2,1) = 0._db; Rot(2,2) = Cos_t; Rot(2,3) = - Sin_t
  Rot(3,1) = 0._db; Rot(3,2) = Sin_t; Rot(3,3) = Cos_t
  Rot = Transpose( Rot )  ! it is to transform vec and pol in lab frame to crystal frame 
  
  Rot = Matmul( Mat_rot, Rot )
  
  Vec_i(1) = 0._db;   Vec_i(2) = 1._db;  Vec_i(3) = 0._db  
  Vec_s(1) = 0._db;   Vec_s(2) = Cos_2t; Vec_s(3) = Sin( 2 * Theta )
  Pol_i_s(1) = 0._db; Pol_i_s(2) = 0._db; Pol_i_s(3) = -1._db   
  Pol_i_p(1) = 1._db; Pol_i_p(2) = 0._db; Pol_i_p(3) = 0._db   
  Pol_s_s(1) = Sin( 2 * Theta ); Pol_s_s(2) = 0._db; Pol_s_s(3) = - Cos_2t   
  Pol_s_p(1) = 1._db; Pol_s_p(2) = 0._db; Pol_s_p(3) = 0._db   

  Vec_i = Matmul( Rot, Vec_i )
  Vec_s = Matmul( Rot, Vec_s )
  Pol_i_s = Matmul( Rot, Pol_i_s )
  Pol_i_p = Matmul( Rot, Pol_i_p )
  Pol_s_s = Matmul( Rot, Pol_s_s )
  Pol_s_p = Matmul( Rot, Pol_s_p )

  if( icheck > 1 ) then
    write(3,'(/A)') '   Pol_i_s   Pol_i_p   Pol_s_s   Pol_s_p    Vec_i     Vec_s'
    do i = 1,3
      write(3,'(6f10.5)'), Pol_i_s(i), Pol_i_p(i), Pol_s_s(i), Pol_s_p(i), Vec_i(i), Vec_s 
    end do   
  endif
  return
  110 format(' Q direction number,',i2,' Q =',3f10.5)
end

!***************************************************************

! Calculation of the inelastic scattering amplitude

subroutine Cal_Tau_loss(Ampl_abs,E_loss,Energ,Energ_s,icheck,ndim2,nenerg_s,nlm_probe,nlmomax,ns_dipmag,nspin,nspino, &
                         nspinp,Spinorbite,Taull)

  use declarations
  implicit none

  integer:: icheck, ie_cond, ie_val, iso1, iso2, isp, isp1, isp2, ispf1, ispf2, iss1, iss2, je, je_loss, L, lm1, &
            lmf1, lms1, lm2, Lmax, lmf2, lms2, m, &
            ndim2, ne, nenerg_s, nlm_probe, nlmomax, ns_dipmag, nspin, nspino, nspinp

  complex(kind=db), dimension(nenerg_s,nlm_probe,nspinp,nlmomax,nspinp):: Ampl_abs
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull

  logical:: Spinorbite
  
  real(kind=db):: E_loss, E_val, Energ
  real(kind=db), dimension(2):: p
  real(kind=db), dimension(nenerg_s):: Energ_s
  
  E_val = Energ - E_loss

  do ie_cond = 1,nenerg_s
    if( Energ_s(ie_cond) > Energ - eps10 ) exit
  end do
    
  do ie_val = 1,nenerg_s
    if( Energ_s(ie_val) > E_val - eps10 ) exit
  end do
  
  if( abs( Energ_s(ie_val) - E_val ) > eps10 ) then
    ne = 2
    ie_val = max( ie_val - 1, 1 )
    p(1) = ( Energ_s(ie_val+1) - E_val ) / ( Energ_s(ie_val+1) - Energ_s(ie_val) )  
    p(2) = 1._db - p(1)  
  else
    ne = 1
    p(1) = 1._db
  endif  

  Taull(:,:,:,:,:,:,:) = ( 0._db, 0._db )
  
  lms1 = 0
  do lm1 = 1,nlm_probe
    do iso1 = 1,nspino
      lms1 = lms1 + 1
      lms2 = 0
      do lm2 = 1,nlm_probe
        do iso2 = 1,nspino
          lms2 = lms2 + 1
          do isp1 = 1,2
            do isp2 = 1,2
              if( Spinorbite ) then
                iss1 = iso1; iss2 = iso2
              else
!                if( isp1 /= isp2 ) cycle
                iss1 = min(isp1,nspin)
                iss2 = min(isp2,nspin)
              endif

              do je = 1,ne
                je_loss = ie_val + je - 1
                do ispf1 = 1,nspinp
                  if( .not. Spinorbite .and. ispf1 /= isp1 ) cycle
                  do ispf2 = 1,nspinp
                    if( .not. Spinorbite .and. ispf2 /= isp2 ) cycle
                    do lmf1 = 1,nlmomax
                      do lmf2 = 1,nlmomax
   ! multiplication by "- img" in order it looks like an inelatic multiple scattering amplitude
                        Taull(lms1,lms2,1,1,isp1,isp2,1) = Taull(lms1,lms2,1,1,isp1,isp2,1) &
                                                         - img *  p(je) * Ampl_abs(je_loss,lm1,iss1,lmf1,ispf1) &
                                                                        * Conjg( Ampl_abs(ie_cond,lm2,iss2,lmf2,ispf2) )
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
  end do  
  
  if( icheck > 1 ) then
    write(3,110) Energ*Rydb, E_val*Rydb, E_loss*Rydb, ie_cond, ie_val, p(1:ne)
    Lmax = nint( sqrt( real( nlm_probe, db ) ) ) - 1
    
    if( nspinp == 1 ) then
      write(3,120) (( L, m, m = -L,L ), L = 0,lmax )
    else
      write(3,130) ((( L, m, isp, m = -L,L ), L = 0,lmax ), isp = 1,nspinp )
    endif
    
    do isp1 = 1,nspinp
      lm1 = 0
      do L = 0,Lmax
        do m = -L, L
          lm1 = lm1 + 1
          if( nspinp == 1 ) then
            write(3,140) L, m, (( Taull(lm1,lm2,1,1,isp1,isp2,1), lm2 = 1, nlm_probe ), isp2 = 1, nspinp )
          else
            write(3,150) L, m, isp1, (( Taull(lm1,lm2,1,1,isp1,isp2,1), lm2 = 1, nlm_probe ), isp2 = 1, nspinp )
          endif 
        end do
      end do
    end do
  endif
  110 format(/' Tau inelatic, E_cond =',f10.5,' eV, E_val =',f10.5,' eV,   E_loss =',f10.5,' eV,   ie_cond, ie_val =',2i5, &
              ',   p =',2f7.4)
  120 format(6x,100(12x,2i2,11x))  
  130 format(9x,100(11x,3i2,10x))  
  140 format(2i3,1p,100(1x,2e13.5))
  150 format(3i3,1p,100(1x,2e13.5))
  return
end

!***************************************************************

! Writing of the tensor

subroutine write_tensor_rixs(i_Spin_channel,ie_cond,ie_loss,E_loss,E1E1,E2E2,Energ_cond, &
                             mpinodes,n_rel,ninitr,secdd,secqq)

  use declarations
  implicit none

  integer:: i, i_Spin_channel, ie_cond, ie_loss, initr, j, k, n_rel, ninitr, mpinodes

  character(len=19):: mot19

  complex(kind=db), dimension(3,3,3,3,ninitr,0:mpinodes-1):: secqq
  complex(kind=db), dimension(3,3,n_rel,ninitr,0:mpinodes-1):: secdd
  
  logical:: E1E1, E2E2
  
  real(kind=db):: E_loss, Energ_cond

  select case(i_Spin_channel)
    case(1)
      mot19 = ' '
    case(2)
      mot19 = ' Spin channel up-up'
    case(3)
      mot19 = ' Spin channel up-dn'
    case(4)
      mot19 = ' Spin channel dn-up'
    case(5)
      mot19 = ' Spin channel dn-dn'
  end select
  write(3,110) ie_cond, ie_loss, Energ_cond*Rydb, E_loss*Rydb, mot19
  
  if( E1E1 ) then
    if( ninitr > 1 ) then
      write(3,'(/A,10(25x,a7,i2,49x))') '   Tensor E1E1', ( 'initl =',initr, initr = 1,ninitr )
    else
      write(3,'(/A)') '   Tensor E1E1'
    endif
    do i = 1,3
      write(3,'(1p,10(2x,3(1x,2e13.5)))') ( secdd(i,:,1,initr,0), initr  = 1,ninitr )
    end do
  endif
  if( E2E2 ) then
    write(3,'(/A)') ' Tensor E2E2'
    do k = 1,3
      write(3,120) k, k, k
      do i = 1,3
        write(3,'(1p,10(2x,3(3x,3(1x,2e13.5))))') (( secqq(i,:,j,k,initr,0), j = 1,3 ), initr = 1,ninitr )
      end do
    end do
  endif

  return
  110 format(/' ie_cond, ie_loss =',2i4,', E_cond =',f10.5,', E_loss =',f10.5,', ',a19) 
  120 format(40x,'(i,j,1,',i1,')',75x,'(i,j,1,',i1,')',75x,'(i,j,1,',i1,')')
end
!**************************************************************************************************************************************

! RIXS matrix elements for the equivalent atoms

subroutine Tensor_equiv(Green_int,ie_computer, &
            isymeq,mpinodes,Multipole,n_oo,n_rel,natomsym,ninitr,nspin, &
            secdd_a,secdd_m_a,secddia,secddia_m,secdq_a,secdq_m_a,secdqia,secdqia_m,secdo_a,secdo_m_a,secdoia,secdoia_m, &
            secmd_a,secmd_m_a,secmdia,secmdia_m,secmm_a,secmm_m_a,secmmia,secmmia_m,secoo_a,secoo_m_a,secooia,secooia_m, &
            secqq_a,secqq_m_a,secqqia,secqqia_m,Spinorbite)

  use declarations
  implicit none

  integer:: ia, ie_computer, initr, ispfg, isym, je, jhe, jhs, js, he, hs, mpinodes, n_oo, n_rel, natomsym, ninitr,  nspin

  integer, dimension(natomsym):: isymeq

  complex(kind=db), dimension(3,3):: sec_dd, sec_md, sec_mm
  complex(kind=db), dimension(3,3,3):: sec_dq
  complex(kind=db), dimension(3,3,3,3):: sec_do, sec_qq
  complex(kind=db), dimension(3,3,3,3,3,3):: Mat6
  complex(kind=db), dimension(3,3,n_rel,ninitr,0:natomsym):: secddia, secddia_m
  complex(kind=db), dimension(3,3,ninitr,0:natomsym):: secmdia, secmdia_m, secmmia, secmmia_m
  complex(kind=db), dimension(3,3,3,ninitr,0:natomsym):: secdqia, secdqia_m
  complex(kind=db), dimension(3,3,3,3,ninitr,0:natomsym):: secdoia, secdoia_m, secqqia, secqqia_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:natomsym):: secooia, secooia_m
  complex(kind=db), dimension(3,3,n_rel,ninitr,0:mpinodes-1):: secdd_a, secdd_m_a
  complex(kind=db), dimension(3,3,ninitr,0:mpinodes-1):: secmd_a, secmd_m_a, secmm_a, secmm_m_a
  complex(kind=db), dimension(3,3,3,ninitr,0:mpinodes-1):: secdq_a, secdq_m_a
  complex(kind=db), dimension(3,3,3,3,ninitr,0:mpinodes-1):: secdo_a, secdo_m_a, secqq_a, secqq_m_a
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:mpinodes-1):: secoo_a, secoo_m_a

  logical:: Dip_rel, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Green_int, Green_int_mag, M1M1, Magn_sens, Spinorbite

  logical, dimension(10):: Multipole

  real(kind=db), dimension(3,3):: matopsym

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  Dip_rel = n_rel > 1
  
  if( nspin == 2 .or. Spinorbite ) then
    Magn_sens = .true.
  else
    Magn_sens = .false.
  endif

  if( Green_int .and. Magn_sens ) then
    Green_int_mag = .true.
  else
    Green_int_mag = .false.
  endif
  
  if( E1E1 ) secddia(:,:,:,:,:) = (0._db, 0._db)
  if( E1E1 .and. Green_int_mag ) secddia_m(:,:,:,:,:) = (0._db, 0._db)
  if( E1M1 ) secmdia(:,:,:,:) = (0._db, 0._db)
  if( E1M1 ) secmdia_m(:,:,:,:) = (0._db, 0._db)
  if( M1M1 ) secmmia(:,:,:,:) = (0._db, 0._db)
  if( M1M1 .and. Green_int_mag ) secmmia_m(:,:,:,:) = (0._db, 0._db)
  if( E1E2 ) secdqia(:,:,:,:,:) = (0._db, 0._db)
  if( E1E2 ) secdqia_m(:,:,:,:,:) = (0._db, 0._db)
  if( E2E2 ) secqqia(:,:,:,:,:,:) = (0._db, 0._db)
  if( E2E2 .and. Green_int_mag ) secqqia_m(:,:,:,:,:,:) = (0._db, 0._db)
  if( E1E3 ) secdoia(:,:,:,:,:,:) = (0._db, 0._db)
  if( E1E3 ) secdoia_m(:,:,:,:,:,:) = (0._db, 0._db)
  if( E3E3 ) secooia(:,:,:,:,:,:) = (0._db, 0._db)
  if( E3E3 .and. Green_int_mag ) secooia_m(:,:,:,:,:,:) = (0._db, 0._db)

  do initr = 1,ninitr       ! ----------> Loop over edges or core states

    do ia = 1,natomsym

      isym = abs( isymeq(ia) )
      call opsym(isym,matopsym)

      if( E1E1 ) then
        do ispfg = 1,n_rel
          sec_dd(:,:) = secdd_a(:,:,ispfg,initr,ie_computer)
          if( isym /= 1 ) call rot_tensor_2( sec_dd, matopsym )
          if( isymeq(ia) < 0 .and. .not. Green_int ) sec_dd(:,:) = conjg( sec_dd(:,:) )
          secddia(:,:,ispfg,initr,ia) = sec_dd(:,:)
          if( Green_int_mag ) then
            sec_dd(:,:) = secdd_m_a(:,:,ispfg,initr,ie_computer)
            if( isym /= 1 ) call rot_tensor_2( sec_dd, matopsym )
            if( isymeq(ia) < 0 ) sec_dd(:,:) = - sec_dd(:,:)
            secddia_m(:,:,ispfg,initr,ia) = sec_dd(:,:)
          endif
        end do
      endif

      if( E1E2 ) then
        sec_dq(:,:,:) = secdq_a(:,:,:,initr,ie_computer)
        if( isym /= 1 ) call rot_tensor_3( sec_dq, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_dq(:,:,:) = conjg( sec_dq(:,:,:) )
        secdqia(:,:,:,initr,ia) = sec_dq(:,:,:)
        if( Green_int_mag ) then
          sec_dq(:,:,:) = secdq_m_a(:,:,:,initr,ie_computer)
          if( isym /= 1 ) call rot_tensor_3( sec_dq, matopsym )
          if( isymeq(ia) < 0 ) sec_dq(:,:,:) = - sec_dq(:,:,:)
          secdqia_m(:,:,:,initr,ia) = sec_dq(:,:,:)
        else
          secdqia_m(:,:,:,initr,ia) = (0._db,0._db)
        endif
      endif

      if( E2E2 ) then
        sec_qq(:,:,:,:) = secqq_a(:,:,:,:,initr,ie_computer)
        if( isym /= 1 ) call rot_tensor_4( sec_qq, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_qq(:,:,:,:) = conjg( sec_qq(:,:,:,:) )
        secqqia(:,:,:,:,initr,ia) = sec_qq(:,:,:,:)
        if( Green_int_mag ) then
          sec_qq(:,:,:,:) = secqq_m_a(:,:,:,:,initr,ie_computer)
          if( isym /= 1 ) call rot_tensor_4( sec_qq, matopsym )
          if( isymeq(ia) < 0 ) sec_qq(:,:,:,:) = - sec_qq(:,:,:,:)
          secqqia_m(:,:,:,:,initr,ia) = sec_qq(:,:,:,:)
        endif
      endif

      if( E1E3 ) then
        sec_do(:,:,:,:) = secdo_a(:,:,:,:,initr,ie_computer)
        if( isym /= 1 ) call rot_tensor_4( sec_do, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_do(:,:,:,:) = conjg( sec_do(:,:,:,:) )
        secdoia(:,:,:,:,initr,ia) = sec_do(:,:,:,:)
        if( Green_int_mag ) then
          sec_do(:,:,:,:) = secdo_m_a(:,:,:,:,initr,ie_computer)
          if( isym /= 1 ) call rot_tensor_4( sec_do, matopsym )
          if( isymeq(ia) < 0 ) sec_do(:,:,:,:) = - sec_do(:,:,:,:)
          secdoia_m(:,:,:,:,initr,ia) = sec_do(:,:,:,:)
        endif
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
                Mat6(:,je,he,:,js,hs) = secoo_a(:,jhe,:,jhs,initr,ie_computer)
              end do
            end do
          end do
        end do
        if( isym /= 1 ) call rot_tensor_6( Mat6, Matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) Mat6 = conjg( Mat6 )
        jhe = 0
        do je = 1,3
          do he = 1,3
            jhe = jhe + 1
            jhs = 0
            do js = 1,3
              do hs = 1,3
                jhs = jhs + 1
                secooia(:,jhe,:,jhs,initr,ia) = Mat6(:,je,he,:,js,hs)
              end do
            end do
          end do
        end do
        if( Green_int_mag ) then
          jhe = 0
          do je = 1,3
            do he = 1,3
              jhe = jhe + 1
              jhs = 0
              do js = 1,3
                do hs = 1,3
                  jhs = jhs + 1
                  Mat6(:,je,he,:,js,hs) = secoo_m_a(:,jhe,:,jhs,initr,ie_computer)
                end do
              end do
            end do
          end do
          if( isym /= 1 ) call rot_tensor_6( Mat6, Matopsym )
          if( isymeq(ia) < 0 .and. .not. Green_int ) Mat6 = conjg( Mat6 )
          jhe = 0
          do je = 1,3
            do he = 1,3
              jhe = jhe + 1
              jhs = 0
              do js = 1,3
                do hs = 1,3
                  jhs = jhs + 1
                  secooia_m(:,jhe,:,jhs,initr,ia) = Mat6(:,je,he,:,js,hs)
                end do
              end do
            end do
          end do
        endif
      endif

      if( E1M1 ) then
        sec_md(:,:) = - secmd_a(:,:,initr,ie_computer)   ! note the change of sign
        if( isym /= 1 ) call rot_tensor_2( sec_md, matopsym )
   ! Plane and inverse symmetries change the sign
        if( ( isym >= 25 .and. isym <= 48 ) .or. ( isym >= 53 .and. isym <= 57 ) .or. isym == 59  .or. isym == 61 &
                          .or. isym == 63 ) sec_md(:,:) = - sec_md(:,:)
   
   ! Magnetism in on the real part
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_md(:,:) = - conjg( sec_md(:,:) )
        secmdia(:,:,initr,ia) = sec_md(:,:)
        if( Green_int_mag ) then
          sec_md(:,:) = - secmd_m_a(:,:,initr,ie_computer) ! note the change of sign
          if( isym /= 1 ) call rot_tensor_2( sec_md, matopsym )
          if( ( isym >= 25 .and. isym <= 48 ) .or. ( isym >= 53 .and. isym <= 57 ) .or. isym == 59  .or. isym == 61 &
                          .or. isym == 63 ) sec_md(:,:) = - sec_md(:,:)
          if( isymeq(ia) < 0 ) sec_md(:,:) = - sec_md(:,:)
          secmdia_m(:,:,initr,ia) = sec_md(:,:)
        else
          secmdia_m(:,:,initr,ia) = (0._db,0._db)
        endif
      endif

      if( M1M1 ) then
        sec_mm(:,:) = secmm_a(:,:,initr,ie_computer)
        if( isym /= 1 ) call rot_tensor_2( sec_mm, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_mm(:,:) = conjg( sec_mm(:,:) )
        secmmia(:,:,initr,ia) = sec_mm(:,:)
        if( Green_int_mag ) then
          sec_mm(:,:) = secmm_m_a(:,:,initr,ie_computer)
          if( isym /= 1 ) call rot_tensor_2( sec_mm, matopsym )
          if( isymeq(ia) < 0 ) sec_mm(:,:) = - sec_mm(:,:)
          secmmia_m(:,:,initr,ia) = sec_mm(:,:)
        endif
      endif

    end do ! end of loop over ia

  end do ! end of loop over initr

  return

end

!*****************************************************************************************************************

subroutine Cal_RIXS_ampl(Amplitude,Dip_rel,E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,M1M1,Magn_sens, &
                         n_rel,npl,Omega_i,Omega_s,pol_i_s,pol_i_p,pol_s_s,pol_s_p,Q_resolved, &
                         RIXS_E1E1_tens_xtal,RIXS_E1E2_tens_xtal,RIXS_E1E2_tens_xtal_m,RIXS_E1E3_tens_xtal,&
                         RIXS_E2E2_tens_xtal,RIXS_E3E3_tens_xtal,RIXS_E1M1_tens_xtal,RIXS_E1M1_tens_xtal_m,RIXS_M1M1_tens_xtal, &
                         Vec_i,Vec_s)

  use declarations
  implicit none

  integer:: he, hs, ipl, isp1, isp2, isp3, isp4, ispfg, j1, je, jhe, jhs, js, ke, ks, n_rel, npl  
   
  complex(kind=db):: matpole, matpols
  complex(kind=db), dimension(npl):: Amplitude
  complex(kind=db), dimension(3):: pol_i, pol_s, u_i, u_s
  complex(kind=db), dimension(3,3):: RIXS_E1M1_tens_xtal, RIXS_E1M1_tens_xtal_m, RIXS_M1M1_tens_xtal
  complex(kind=db), dimension(3,3,n_rel):: RIXS_E1E1_tens_xtal
  complex(kind=db), dimension(3,3,3):: RIXS_E1E2_tens_xtal, RIXS_E1E2_tens_xtal_m
  complex(kind=db), dimension(3,3,3,3):: RIXS_E1E3_tens_xtal, RIXS_E2E2_tens_xtal, RIXS_E3E3_tens_xtal
  
  logical:: Dip_rel, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, M1M1, Magn_sens, Q_resolved
  
  real(kind=db):: Fac, Omega_i, Omega_s
  real(kind=db), dimension(3):: pol_i_s, pol_i_p, pol_s_s, pol_s_p, Vec_i, Vec_s

  Amplitude(:) = (0._db, 0._db )
    
  if( .not. Q_resolved ) then
 
    if( E1E1 ) then
      Fac = 1 / 3._db
      Amplitude(1) = Amplitude(1) + Fac * ( RIXS_E1E1_tens_xtal(1,1,1) + RIXS_E1E1_tens_xtal(2,2,1) + RIXS_E1E1_tens_xtal(3,3,1) )
! I do the hypothesis than in that there is no spin crossing contribution. I am not sure.
      if( n_rel /= 1 ) Amplitude(1) = Amplitude(1) + Fac * ( RIXS_E1E1_tens_xtal(1,1,n_rel) + RIXS_E1E1_tens_xtal(2,2,n_rel) &
                                                                                            + RIXS_E1E1_tens_xtal(3,3,n_rel) )
    endif

    if( E2E2 ) then
      Fac = 1 / 45._db
      Fac = Fac / ( 2._db * sqrt( 5._db) ) 
      Amplitude(1) = Amplitude(1) &
                   + Fac * ( 6 * ( RIXS_E2E2_tens_xtal(1,3,1,3) + RIXS_E2E2_tens_xtal(2,3,2,3) + RIXS_E2E2_tens_xtal(1,2,1,2) ) &
                           + 2 * ( RIXS_E2E2_tens_xtal(1,1,1,1) + RIXS_E2E2_tens_xtal(2,2,2,2) + RIXS_E2E2_tens_xtal(3,3,3,3) ) &
                               - ( RIXS_E2E2_tens_xtal(1,1,2,2) + RIXS_E2E2_tens_xtal(1,1,3,3) + RIXS_E2E2_tens_xtal(3,3,2,2) &
                                 + RIXS_E2E2_tens_xtal(2,2,1,1) + RIXS_E2E2_tens_xtal(3,3,1,1) + RIXS_E2E2_tens_xtal(2,2,3,3) ) )
    endif

    if( M1M1 ) then
      Fac = 1 / 3._db
      Amplitude(1) = Amplitude(1) + Fac * ( RIXS_M1M1_tens_xtal(1,1) + RIXS_M1M1_tens_xtal(2,2) + RIXS_M1M1_tens_xtal(3,3) )
    endif

    if( E3E3 ) then  ! not checked : taken like E2E2...
      Fac = 1 / 45._db
      Fac = Fac / ( 2._db * sqrt( 5._db) ) 
      Amplitude(1) = Amplitude(1) &
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
          pol_i(:) = ( pol_i_s(:) + img * pol_i_p(:) ) / sqrt( 2._db ) 
          pol_s(:) = cmplx( pol_s_p(:), 0._db, db )
        case(7)
          pol_i(:) = ( pol_i_s(:) - img * pol_i_p(:) ) / sqrt( 2._db ) 
          pol_s(:) = cmplx( pol_s_s(:), 0._db, db )
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
          do ks = 1,3   ! Here one gets back img
            Amplitude(ipl) = Amplitude(ipl) - img * conjg( Pol_s(ks) ) * Pol_i(ke) &
                      * sum( Vec_i(:) * RIXS_E1E2_tens_xtal(ks,ke,:) - Vec_s(:) * RIXS_E1E2_tens_xtal(ke,ks,:) )
            if( Magn_sens ) Amplitude(ipl) = Amplitude(ipl) - img *  conjg( Pol_s(ks) ) * Pol_i(ke) &
                      * sum( Vec_i(:) * RIXS_E1E2_tens_xtal_m(ks,ke,:) + Vec_s(:) * RIXS_E1E2_tens_xtal_m(ke,ks,:) )
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
                        * ( Vec_i(j1) * sum( Vec_i(:) * RIXS_E1E3_tens_xtal(ks,ke,j1,:) ) &
                          + Vec_s(j1) * sum( Vec_s(:) * conjg( RIXS_E1E3_tens_xtal(ke,ks,j1,:) ) ) )
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
  ! on met le signe moins devant uae car la partie imaginaire contient le terme non-magnetique, equivaut au complexe conjugue
  ! sauf si secmdia est calcule pour une reflexion dafs ( la partie reelle est nulle sauf dafs )
          Amplitude(ipl) = Amplitude(ipl) + conjg( u_s(ke) ) * sum( Pol_i(:) * RIXS_E1M1_tens_xtal(:,ke) ) &
                                          - u_i(ke) * sum( conjg(Pol_s(:)) * RIXS_E1M1_tens_xtal(:,ke) )  
  ! on met le signe plus devant uae car la partie reelle contient le terme magnetique, (le complexe conjugue est donc inutile)
          if( Magn_sens ) Amplitude(ipl) = Amplitude(ipl) + conjg( u_s(ke) ) * sum( Pol_i(:) * RIXS_E1M1_tens_xtal_m(:,ke) ) &
                                                          + u_i(ke) * sum( conjg( Pol_s(:) ) * RIXS_E1M1_tens_xtal_m(:,ke) )
        end do
      endif
              
    end do  ! end of loop over polarization ipl

  endif
  
  return
end

!*****************************************************************************************************************

subroutine Write_rixs(E_loss,Energ,Epsii,Eseuil,iabsorig,File_name_rixs,jseuil,n_multi_run,n_q_dim,n_Spin_channel,n_theta, &
                      ne_loss,nenerg,ninit1,ninitr,nomfich,npl,nseuil,numat,RIXS_ampl)

  use declarations
  implicit none
 
  character(len=Length_name):: File_name_rixs, nomfich
  
  integer:: i_q, i_t, iabsorig, ie_loss, ii, ipl, jseuil, Length, n_isc, n_multi_run, n_q_dim, n_Spin_channel, n_theta, ne_loss, &
            nenerg, ninit1, ninitr, nseuil, numat, npl  
   
  complex(kind=db), dimension(ne_loss,nenerg,npl,n_theta,n_q_dim,n_Spin_channel*ninitr):: RIXS_ampl
 
  real(kind=db):: Eseuil 
  real(kind=db), dimension(ninitr):: Epsii
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(ne_loss):: E_loss

  n_isc = n_Spin_channel * ninitr
    
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

  open(2, file = File_name_rixs )
  write(2,110) Eseuil*rydb, numat, nseuil, jseuil, nenerg, ne_loss, npl, n_theta, n_q_dim, ninit1, ninitr, n_Spin_channel, &
               Epsii(:)*rydb
  write(2,120) ( Energ(:)*rydb, ii = 1,ninitr * ninitr )

  do ie_loss = 1,ne_loss
    if( abs( E_loss(ie_loss)*Rydb ) < 9.999995_db ) then
      write(2,165) E_loss(ie_loss)*Rydb, (((( RIXS_ampl(ie_loss,:,ipl,i_t,i_q,ii), ipl = 1,npl ), i_t = 1,n_theta ), &
                                                                                i_q = 1,n_q_dim ), ii = 1,n_isc ) 
    elseif( abs( E_loss(ie_loss)*Rydb ) < 99.99995_db ) then
      write(2,175) E_loss(ie_loss)*Rydb, (((( RIXS_ampl(ie_loss,:,ipl,i_t,i_q,ii), ipl = 1,npl ), i_t = 1,n_theta ), &
                                                                                i_q = 1,n_q_dim ), ii = 1,n_isc ) 
    elseif( abs( E_loss(ie_loss)*Rydb ) < 999.9995_db ) then
      write(2,185) E_loss(ie_loss)*Rydb, (((( RIXS_ampl(ie_loss,:,ipl,i_t,i_q,ii), ipl = 1,npl ), i_t = 1,n_theta ), &
                                                                                i_q = 1,n_q_dim ), ii = 1,n_isc )
    else
      write(2,195) E_loss(ie_loss)*Rydb, (((( RIXS_ampl(ie_loss,:,ipl,i_t,i_q,ii), ipl = 1,npl ), i_t = 1,n_theta ), &
                                                                                i_q = 1,n_q_dim ), ii = 1,n_isc )
    endif 
  end do
  
  Close(2)

  return
  110 format(f10.3,i5,2i3,2i6,i3,i6,4i3,1p,10e15.7)
  120 format(10x,100000(8x,f15.5,8x))
  165 format(f10.5,1p,100000(1x,2e15.7))
  175 format(f10.4,1p,100000(1x,2e15.7))
  185 format(f10.3,1p,100000(1x,2e15.7))
  195 format(f10.2,1p,100000(1x,2e15.7))
end

!*****************************************************************************************************************

subroutine Write_Int_rixs(File_name,n_File,RIXS_core)

  use declarations
  implicit none
 
  integer:: i_File, i_q, i_t, ie, ie_loss, i_Spin_channel, ii, initr, ipl, ipr, istat, jseuil, Length, n_File, n_i, &
            n_isc, n_Spin_channel, n_theta, n_q_dim, ne_loss, nenerg, ninit1, ninitr, npl, nseuil, numat  
 
  character(len=Length_word):: Word  
  character(len=Length_name):: File_name_out, mot
  character(len=Length_name), dimension(n_File):: File_name
  character(len=15), dimension(:), allocatable:: E_name
  
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: RIXS_ampl
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: RIXS_ampl_tot

  logical:: RIXS_core
  
  real(kind=db):: E, Eseuil 
  real(kind=db), dimension(:), allocatable:: E_loss, Energ
  real(kind=db), dimension(:,:,:,:,:), allocatable:: Ampl_i, Ampl_r

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
    read(2,*) Eseuil, numat, nseuil, jseuil, nenerg, ne_loss, npl, n_theta, n_q_dim, ninit1, ninitr, n_Spin_channel
    Close(2)
  end do

  if( RIXS_core .and. ninitr > 1 ) then
    n_i = ninitr
  else
    n_i = 0
  endif

  n_isc = ninitr * n_Spin_channel
   
  allocate( Energ(nenerg) )
  allocate( E_loss(ne_loss) )
  allocate( Ampl_i(nenerg,npl,n_theta,n_q_dim,n_isc) )
  allocate( Ampl_r(nenerg,npl,n_theta,n_q_dim,n_isc) )
  allocate( RIXS_ampl(ne_loss,nenerg,npl,n_theta,n_q_dim,n_isc) )
  allocate( RIXS_ampl_tot(ne_loss,nenerg,npl,n_theta,n_q_dim,0:n_i,n_Spin_channel) )
  allocate( E_name(nenerg*(n_i+1)*n_Spin_channel) )
  
  do i_File = 1,n_File
    open(2, file = File_name(i_File), status='old')
    read(2,*) Eseuil, numat, nseuil, jseuil, nenerg, ne_loss, npl, n_theta, n_q_dim, ninit1, ninitr
    read(2,*) Energ(:)
    Energ(:) = Energ(:) / Rydb
    Close(2)
  end do
    
  RIXS_ampl_tot(:,:,:,:,:,:,:) = ( 0._db, 0._db )
    
  do i_File = 1,n_File
    open(2, file = File_name(i_File), status='old')
    read(2,*)
    read(2,*)
    do ie_loss = 1,ne_loss
      read(2,*) E_loss(ie_loss), ((((( Ampl_r(ie,ipl,i_t,i_q,ii), Ampl_i(ie,ipl,i_t,i_q,ii), ie = 1,nenerg ), &
                                ipl = 1,npl ), i_t = 1,n_theta ), i_q = 1,n_q_dim ), ii = 1,n_isc)
      RIXS_ampl(ie_loss,:,:,:,:,:) = cmplx( Ampl_r(:,:,:,:,:), Ampl_i(:,:,:,:,:), db )  
    end do
    E_loss(:) = E_loss(:) / Rydb

    do ie_loss = 1,ne_loss
      
      do i_Spin_channel = 1,n_Spin_channel
        do initr = 1,ninitr
          ii = initr + ( i_Spin_channel - 1 ) * ninitr

          if( n_i > 0 ) RIXS_ampl_tot(ie_loss,:,:,:,:,initr,i_Spin_channel) = RIXS_ampl_tot(ie_loss,:,:,:,:,initr,i_Spin_channel) &
                                                                            + RIXS_ampl(ie_loss,:,:,:,:,ii)
          RIXS_ampl_tot(ie_loss,:,:,:,:,0,i_Spin_channel) = RIXS_ampl_tot(ie_loss,:,:,:,:,0,i_Spin_channel) &
                                                          + RIXS_ampl(ie_loss,:,:,:,:,ii) 

        end do    
      end do
    end do
    close(2)
  end do

! Writing

  ii = 0
  do i_Spin_channel = 1,n_Spin_channel
    do initr = 0,n_i
      do ie = 1,nenerg
        ii = ii + 1
        E_name(ii) = ' '
        E = Energ(ie) *rydb
        Word = ' '
        if( abs( nint( 10*E ) - 10*E ) < eps10 ) then
          if( Length_word == 15 ) then
            write(Word,'(f15.1)') E
          elseif( Length_word == 17 ) then
            write(Word,'(f17.1)') E
          endif
        elseif( abs( nint( 100*E ) - 100*E ) < eps10 ) then
          if( Length_word == 15 ) then
            write(Word,'(f15.2)') E
          elseif( Length_word == 17 ) then
            write(Word,'(f17.2)') E
          endif
        elseif( abs( nint( 1000*E ) - 1000*E ) < eps10 ) then
           if( Length_word == 15 ) then
            write(Word,'(f15.3)') E
          elseif( Length_word == 17 ) then
            write(Word,'(f17.3)') E
          endif
        elseif( abs( nint( 1000*E ) - 1000*E ) < eps10 ) then
          if( Length_word == 15 ) then
            write(Word,'(f15.4)') E
          elseif( Length_word == 17 ) then
            write(Word,'(f17.4)') E
          endif
        endif
        Word = adjustl(Word)
        
        if( initr > 0 ) then
          Length = len_trim(Word)
          Word(Length+1:Length+2) = '_i'
          call ad_number(initr,Word,length_word)
        endif
        if( i_Spin_channel > 1 ) then
          Length = len_trim(Word)
          select case(i_Spin_channel)
            case(2)
              Word(Length+1:Length+3) = '_uu'
            case(3)
              Word(Length+1:Length+3) = '_ud'
            case(4)
              Word(Length+1:Length+3) = '_du'
            case(5)
              Word(Length+1:Length+3) = '_dd'
          end select
        endif
        call center_word(Word,Length_word)
        E_name(ii) = Word
      end do
    end do
  end do

  do i_Spin_channel = 1,n_Spin_channel
   
    File_name_out = File_name(1)
    Length = len_trim(File_name_out)
    if( File_name_out(Length-5:Length-3) == '_1') then
      Length = Length-5
    else
      Length = Length-3
    endif
    
    select case(i_Spin_channel)
      case(1)
        File_name_out(Length:Length+7) = '_int.txt'
      case(2)
        File_name_out(Length:Length+10) = '_int_uu.txt'
      case(3)
        File_name_out(Length:Length+10) = '_int_ud.txt'
      case(4)
       File_name_out(Length:Length+10) = '_int_du.txt'
      case(5)
       File_name_out(Length:Length+10) = '_int_dd.txt'
    end select
  
    open(2, file = File_name_out)
    write(2,120) nenerg, ne_loss, npl, n_theta, n_q_dim, n_i, n_Spin_channel
    write(2,130) ( E_name(ii), ii = (i_Spin_channel-1) * nenerg * ( 1 + n_i ) + 1, i_Spin_channel * nenerg * ( 1 + n_i ) )
 
    do ie_loss = 1,ne_loss
      if( abs( E_loss(ie_loss)*Rydb ) < 9.999995_db ) then
        write(2,165) E_loss(ie_loss)*Rydb, (((( abs(RIXS_ampl_tot(ie_loss,:,ipl,i_t,i_q,initr,i_Spin_channel) )**2, ipl = 1,npl ), &
                                                 i_t = 1,n_theta ), i_q = 1,n_q_dim ), initr = 0, n_i ) 
      elseif( abs( E_loss(ie_loss)*Rydb ) < 99.99995_db ) then
        write(2,175) E_loss(ie_loss)*Rydb, (((( abs(RIXS_ampl_tot(ie_loss,:,ipl,i_t,i_q,initr,i_Spin_channel) )**2, ipl = 1,npl ), &
                                                 i_t = 1,n_theta ), i_q = 1,n_q_dim ), initr = 0, n_i )
      elseif( abs( E_loss(ie_loss)*Rydb ) < 999.9995_db ) then
        write(2,185) E_loss(ie_loss)*Rydb, (((( abs(RIXS_ampl_tot(ie_loss,:,ipl,i_t,i_q,initr,i_Spin_channel) )**2, ipl = 1,npl ), &
                                                 i_t = 1,n_theta ), i_q = 1,n_q_dim ), initr = 0, n_i )
      else
        write(2,195) E_loss(ie_loss)*Rydb, (((( abs(RIXS_ampl_tot(ie_loss,:,ipl,i_t,i_q,initr,i_Spin_channel) )**2, ipl = 1,npl ), &
                                                 i_t = 1,n_theta ), i_q = 1,n_q_dim ), initr = 0, n_i ) 
      endif 
    end do
    Close(2)
  end do
  
  deallocate( Ampl_i, Ampl_r, E_loss, E_name, Energ, RIXS_ampl, RIXS_ampl_tot )
  
  return
  110 format(//' Error opening the file:',//3x,A,//)  
  120 format(7i5,' / ne, ne_loss, npl, n_theta, n_q_dir, n_i, n_Spin_channel')
  130 format('  E_loss ',1000000a15 )
  165 format(f10.5,1p,100000e15.7)
  175 format(f10.4,1p,100000e15.7)
  185 format(f10.3,1p,100000e15.7)
  195 format(f10.2,1p,100000e15.7)
end

