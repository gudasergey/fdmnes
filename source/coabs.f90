! FDMNES subroutine

! Calculate the absorption cross sections and the DAFS amplitudes

subroutine Write_coabs(Abs_U_iso,Allsite,angxyz,axyz,Bragg_abs,Bulk_step,Cartesian_tensor, &
            Core_resolved,Dafs,Dafs_bio,E_cut,Energ,Energphot,Extract_ten,Epsii,Eseuil,Final_tddft,First_E, &
            f_avantseuil,Full_self_abs,Green_int,hkl_dafs,i_range,iabsorig,icheck,ie,ie_computer,igr_bulk_z_abs, &
            Int_tens,isigpi,isymeq,jseuil,Length_abs,Length_rel,Length_rel_abs,ltypcal,Matper,Moyenne,mpinodee, &
            multi_0,Multipole,n_abs_rgh,n_bulk_sup,n_multi_run, &
            n_bulk_z,n_bulk_z_abs,n_bulk_z_max_abs,n_bulk_zc_abs,n_max,n_oo,n_rel,n_tens_max,natomsym,nbseuil, &
            ncolm,ncolr,ncolt,nenerg,ninit1,ninitlr,nomabs,nomfich,nomfich_cal_conv,nomfich_s,nphi_dafs,npldafs, &
            nphim,nplr,nplrm,nseuil,nspin,numat_abs,nxanout,pdp,phdf0t,phdt,pol,poldafse,poldafss, &
            sec_atom,secdd_a,secdd_m_a,secdq_a,secdq_m_a,secdo_a,secdo_m_a, &
            secmd_a,secmd_m_a,secmm_a,secmm_m_a,secoo_a,secoo_m_a,secqq_a,secqq_m_a,Self_abs,Spherical_signal, &
            Spherical_tensor,Spinorbite_p,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_ref,Taux_eq,V0muf,Vecdafse,Vecdafss, &
            Vec,Volume_maille,Xan_atom)

  use declarations
  implicit none

  real(kind=db), parameter:: quatre_mc2 = 4 * 2 / alfa_sf**2  ! in Rydberg
  
  integer:: he, hs, i, i_bulk_z, i_range, ia, iabsorig, ib, ic1, ic2, icheck, icn1, icn2, ie, &
    ie_computer, ii, ind_mu, initlr, ip, ipl, ipldafs, iseuil, isp1, isp2, isp3, isp4, ispfg, &
    isym, ixandafs, j, j1, je, jhe, jhs, jpl, js, jseuil, ke, ks, l, ll, long, mpinodee, multi_0, n_abs_rgh, &
    n_bulk_sup, n_bulk_z, n_bulk_z_abs, n_bulk_z_max_abs, n_dim, n_mat_cal, n_max, n_multi_run, n_oo, &
    n_rel, n_tens, n_tens_max, n1, n2, na, &
    nab, natomsym, nb, nbseuil, nc, nccm, ncolm, ncolr, ncolt, nenerg, ninit1, ninitlr, nl, np, npldafs, &
    nphim, nplt, nplr, nplrm, npps, nseuil, nspin, numat_abs, nxanout, nw

  integer, dimension(natomsym):: isymeq
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(2,npldafs):: isigpi
  integer, dimension(n_bulk_z_abs):: n_bulk_zc_abs
  integer, dimension(n_bulk_z_max_abs,n_bulk_z_abs):: igr_bulk_z_abs

  character(len=length_word):: nomab
  character(len=length_word), dimension(ncolm):: nomabs
  character(len=length_word), dimension(ncolm*ninitlr):: title
  character(len=5), dimension(nplrm):: ltypcal
  character(len=Length_name):: nomfich, nomfich_s, nomfichdafst, nomficht
  character(len=Length_name), dimension(n_multi_run+n_bulk_sup+n_abs_rgh):: nomfich_cal_conv
  character(len=2310):: mot

  complex(kind=db):: cf, f_avantseuil, matpole, matpols, ph, ph_m, sec
  complex(kind=db), dimension(3):: plae, plas, uae, uas
  complex(kind=db), dimension(8*ninitlr):: compnum
  complex(kind=db), dimension(3,3):: sec_dd, sec_md, sec_mm
  complex(kind=db), dimension(3,3,3):: sec_dq
  complex(kind=db), dimension(3,3,3,3):: sec_do, sec_qq
  complex(kind=db), dimension(3,3,3,3,3,3):: Mat6
  complex(kind=db), dimension(3,3,n_rel,ninitlr,0:natomsym):: secddia, secddia_m
  complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secmdia, secmdia_m, secmmia, secmmia_m
  complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia, secdqia_m
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secdoia, secdoia_m, secqqia, secqqia_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:natomsym):: secooia, secooia_m
  complex(kind=db), dimension(3,3,n_rel,ninitlr,0:mpinodee-1):: secdd_a, secdd_m_a
  complex(kind=db), dimension(3,3,ninitlr,0:mpinodee-1):: secmd_a, secmd_m_a, secmm_a, secmm_m_a
  complex(kind=db), dimension(3,3,3,ninitlr,0:mpinodee-1):: secdq_a, secdq_m_a
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:mpinodee-1):: secdo_a, secdo_m_a, secqq_a, secqq_m_a
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:mpinodee-1):: secoo_a, secoo_m_a
  complex(kind=db), dimension(npldafs):: phdtem, phdf0t1, phdt1
  complex(kind=db), dimension(3,nplrm):: pol
  complex(kind=db), dimension(npldafs,nphim):: phdf0t
  complex(kind=db), dimension(npldafs,nphim,n_max):: phdt
  complex(kind=db), dimension(npldafs,n_max):: Sum_Bragg_bulk_abs_f0
  complex(kind=db), dimension(npldafs,n_bulk_z):: Sum_Bragg_nonabs_f
  complex(kind=db), dimension(natomsym,npldafs):: Bragg_abs
  complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss

  complex(kind=db), dimension(:,:,:,:), allocatable :: Ampldafs, Ampldafs_dd, Ampldafs_do, Ampldafs_dq, Ampldafs_md, &
    Ampldafs_mm, Ampldafs_oo, Ampldafs_qq
  complex(kind=db), dimension(:,:,:,:,:), allocatable :: mu, mu_dd, mu_do, mu_dq, mu_md, mu_mm, mu_oo, mu_qq

  logical:: Allsite, Bulk_roughness_case, Bulk_step, Cartesian_tensor, Cor_abs, Core_resolved, Diag_spin, Dip_rel, E1E1, E1E2, &
    E1E3, E1M1, E2E2, E3E3, E_vec, Dafs, Dafs_bio, Final_tddft, Energphot, Extract_ten, First_E, &
    Full_self_abs, Green_int, Green_int_mag, idafs, M1M1, Magn_sens, Matper, Moyenne, &
    mu_cal, Self_abs, Spherical_signal, Spherical_tensor, Spinorbite_p, Tens_comp, Tensor_eval, Xan_atom

  logical, dimension(10):: Multipole

  real(kind=db):: Abs_U_iso, ang, c_micro, c_milli, cst, conv_mbarn_nelec, dang, E_cut, &
    eph2, Ephoton, Ephseuil, natomsym_f, Omega, ptrans, Surface_ref, V0muf, Volume_maille

  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(ninitlr):: ct_nelec, Epsii
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(3):: angxyz, axyz, voae, voas
  real(kind=db), dimension(3,3):: matopsym
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(ninitlr):: sec_atom
  real(kind=db), dimension(3,nplrm):: vec
  real(kind=db), dimension(nplrm,2):: pdp
  real(kind=db), dimension(3,npldafs):: hkl_dafs
  real(kind=db), dimension(npldafs):: Length_Abs
  real(kind=db), dimension(n_bulk_z):: Length_rel
  real(kind=db), dimension(n_bulk_z_abs):: Length_rel_abs
  real(kind=db), dimension(ncolm*ninitlr):: tens
  real(kind=db), dimension(natomsym):: Taux_eq
  real(kind=db), dimension(n_tens_max*ninitlr,0:natomsym):: Int_tens
  real(kind=db), dimension(ncolr,ninitlr,0:natomsym):: Secabs, Secabs_dd, Secabs_dq, Secabs_do, Secabs_md, Secabs_mm, Secabs_oo, &
                                                       Secabs_qq
  real(kind=db), dimension(3,npldafs,nphim):: vecdafse, vecdafss
  real(kind=db), dimension(:,:), allocatable:: Abs_sum, hkl_dafs_fake

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  Dip_rel = n_rel > 1
  
  if( E1E1 ) Secabs_dd(:,:,:) = ( 0._db, 0._db )
  if( E2E2 ) Secabs_qq(:,:,:) = ( 0._db, 0._db )
  if( E1E3 ) Secabs_do(:,:,:) = ( 0._db, 0._db )
  if( E3E3 ) Secabs_oo(:,:,:) = ( 0._db, 0._db )
  if( M1M1 ) Secabs_mm(:,:,:) = ( 0._db, 0._db )
  if( E1M1 ) Secabs_md(:,:,:) = ( 0._db, 0._db )
  if( E1E2 ) Secabs_dq(:,:,:) = ( 0._db, 0._db )

  if( nspin == 2 .or. Spinorbite_p ) then
    Magn_sens = .true.
  else
    Magn_sens = .false.
  endif

  if( Green_int .and. Magn_sens ) then
    Green_int_mag = .true.
  else
    Green_int_mag = .false.
  endif

  Cor_abs = Full_self_abs .or. Self_abs

  Bulk_roughness_case = n_bulk_z_abs > 0 .and. .not. Bulk_step

  if( icheck > 0 .and. Bulk_roughness_case ) then
    write(3,110)
  elseif( icheck > 0 ) then
    write(3,115)
  endif
 
  natomsym_f = sum( Taux_eq(:) )
    
  if( Allsite ) then
    na = natomsym
    nb = natomsym
  else
    na = 0
    nb = 1
  endif
  if( n_bulk_z_abs > 1 .and. .not. Bulk_roughness_case ) then
    nab = natomsym
  else
    nab = na
  endif

  n_mat_cal = 0
  do ipl = 1,nplr
    if( ltypcal(ipl) == 'ss   ' .or. ltypcal(ipl) == 'sp   ' .or. ltypcal(ipl) == 'ps   ' &
                                .or. ltypcal(ipl) == 'pp   ' ) n_mat_cal = n_mat_cal + 1
  end do
 
  if( Dafs ) then

    allocate( Ampldafs(npldafs,nphim,ninitlr,0:natomsym) )
    if( E1E1 ) allocate( Ampldafs_dd(npldafs,nphim,ninitlr,0:natomsym) )
    if( E1E2 ) allocate( Ampldafs_dq(npldafs,nphim,ninitlr,0:natomsym) )
    if( E2E2 ) allocate( Ampldafs_qq(npldafs,nphim,ninitlr,0:natomsym) )
    if( E1E3 ) allocate( Ampldafs_do(npldafs,nphim,ninitlr,0:natomsym) )
    if( E3E3 ) allocate( Ampldafs_oo(npldafs,nphim,ninitlr,0:natomsym) )
    if( M1M1 ) allocate( Ampldafs_mm(npldafs,nphim,ninitlr,0:natomsym) )
    if( E1M1 ) allocate( Ampldafs_md(npldafs,nphim,ninitlr,0:natomsym) )
    if( E1E1 ) Ampldafs_dd(:,:,:,:) = ( 0._db, 0._db ) 
    if( E1E2 ) Ampldafs_dq(:,:,:,:) = ( 0._db, 0._db )
    if( E2E2 ) Ampldafs_qq(:,:,:,:) = ( 0._db, 0._db )
    if( E1E3 ) Ampldafs_do(:,:,:,:) = ( 0._db, 0._db )
    if( E3E3 ) Ampldafs_oo(:,:,:,:) = ( 0._db, 0._db )
    if( M1M1 ) Ampldafs_mm(:,:,:,:) = ( 0._db, 0._db )
    if( E1M1 ) Ampldafs_md(:,:,:,:) = ( 0._db, 0._db )
    if( Cor_abs ) then
      allocate( mu(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E1E1 ) allocate( mu_dd(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E1E2 ) allocate( mu_dq(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E2E2 ) allocate( mu_qq(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E1E3 ) allocate( mu_do(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E3E3 ) allocate( mu_oo(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( M1M1 ) allocate( mu_mm(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E1M1 ) allocate( mu_md(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E1E1 ) mu_dd(:,:,:,:,:) = ( 0._db, 0._db )
      if( E1E2 ) mu_dq(:,:,:,:,:) = ( 0._db, 0._db )
      if( E2E2 ) mu_qq(:,:,:,:,:) = ( 0._db, 0._db )
      if( E1E3 ) mu_do(:,:,:,:,:) = ( 0._db, 0._db )
      if( E3E3 ) mu_oo(:,:,:,:,:) = ( 0._db, 0._db )
      if( M1M1 ) mu_mm(:,:,:,:,:) = ( 0._db, 0._db )
      if( E1M1 ) mu_md(:,:,:,:,:) = ( 0._db, 0._db )
    endif

  endif

  do initlr = 1,ninitlr       ! ----------> Loop over edges or core states

    if( Core_resolved .and. .not. Final_tddft ) then
      if( initlr <= ninit1 ) then
        iseuil = 1
      else
        iseuil = min(2, nbseuil)
      endif
    elseif( Final_tddft ) then
      iseuil = min(2, nbseuil)
    else
      iseuil = min(initlr, nbseuil)
    endif

    Ephseuil = Energ(ie)
    Ephoton = Ephseuil + Eseuil(iseuil)
! For very low energy edges:
    Ephoton = max(0.001_db/Rydb, Ephoton)
    if( Energphot ) Ephseuil = Ephoton

    ct_nelec(initlr) = conv_mbarn_nelec(Ephoton)

    eph2 = 0.5_db * Ephoton**2
! To get tensors and cross sections in Megabarn
    if( Extract_ten ) then
      cst = 1._db
    else
! ptrans = S02 fixed at 1.
      ptrans = 1._db
! alfa_sf = e*e/(2*epsilon0*h*c) is the fine structure constant
      cst = quatre_pi * pi * ptrans * alfa_sf * Ephoton
! To get result in Mbarn (10E-18 cm2)
      cst = 100 * bohr**2 * cst
    endif
    
    do ia = 1,natomsym

      isym = abs( isymeq(ia) )
      call opsym(isym,matopsym)

      if( E1E1 ) then
        do ispfg = 1,n_rel
          sec_dd(:,:) = secdd_a(:,:,ispfg,initlr,ie_computer) * cst
          if( isym /= 1 ) call rot_tensor_2( sec_dd, matopsym )
          if( isymeq(ia) < 0 .and. .not. Green_int ) sec_dd(:,:) = conjg( sec_dd(:,:) )
          secddia(:,:,ispfg,initlr,ia) = sec_dd(:,:)
          if( Green_int_mag ) then
            sec_dd(:,:) = secdd_m_a(:,:,ispfg,initlr,ie_computer) * cst
            if( isym /= 1 ) call rot_tensor_2( sec_dd, matopsym )
            if( isymeq(ia) < 0 ) sec_dd(:,:) = - sec_dd(:,:)
            secddia_m(:,:,ispfg,initlr,ia) = sec_dd(:,:)
          endif
        end do
      endif

      if( E1E2 ) then
        sec_dq(:,:,:) = secdq_a(:,:,:,initlr,ie_computer) * cst
        if( isym /= 1 ) call rot_tensor_3( sec_dq, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_dq(:,:,:) = conjg( sec_dq(:,:,:) )
        secdqia(:,:,:,initlr,ia) = sec_dq(:,:,:)
        if( Green_int_mag ) then
          sec_dq(:,:,:) = secdq_m_a(:,:,:,initlr,ie_computer) * cst
          if( isym /= 1 ) call rot_tensor_3( sec_dq, matopsym )
          if( isymeq(ia) < 0 ) sec_dq(:,:,:) = - sec_dq(:,:,:)
          secdqia_m(:,:,:,initlr,ia) = sec_dq(:,:,:)
        else
          secdqia_m(:,:,:,initlr,ia) = (0._db,0._db)
        endif
      endif

      if( E2E2 ) then
        sec_qq(:,:,:,:) = secqq_a(:,:,:,:,initlr,ie_computer) * cst
        if( isym /= 1 ) call rot_tensor_4( sec_qq, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_qq(:,:,:,:) = conjg( sec_qq(:,:,:,:) )
        secqqia(:,:,:,:,initlr,ia) = sec_qq(:,:,:,:)
        if( Green_int_mag ) then
          sec_qq(:,:,:,:) = secqq_m_a(:,:,:,:,initlr,ie_computer) * cst
          if( isym /= 1 ) call rot_tensor_4( sec_qq, matopsym )
          if( isymeq(ia) < 0 ) sec_qq(:,:,:,:) = - sec_qq(:,:,:,:)
          secqqia_m(:,:,:,:,initlr,ia) = sec_qq(:,:,:,:)
        endif
      endif

      if( E1E3 ) then
        sec_do(:,:,:,:) = secdo_a(:,:,:,:,initlr,ie_computer) * cst
        if( isym /= 1 ) call rot_tensor_4( sec_do, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_do(:,:,:,:) = conjg( sec_do(:,:,:,:) )
        secdoia(:,:,:,:,initlr,ia) = sec_do(:,:,:,:)
        if( Green_int_mag ) then
          sec_do(:,:,:,:) = secdo_m_a(:,:,:,:,initlr,ie_computer) * cst
          if( isym /= 1 ) call rot_tensor_4( sec_do, matopsym )
          if( isymeq(ia) < 0 ) sec_do(:,:,:,:) = - sec_do(:,:,:,:)
          secdoia_m(:,:,:,:,initlr,ia) = sec_do(:,:,:,:)
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
                Mat6(:,je,he,:,js,hs) = secoo_a(:,jhe,:,jhs,initlr,ie_computer) * cst
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
                secooia(:,jhe,:,jhs,initlr,ia) = Mat6(:,je,he,:,js,hs)
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
                  Mat6(:,je,he,:,js,hs) = secoo_m_a(:,jhe,:,jhs,initlr,ie_computer) * cst
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
                  secooia_m(:,jhe,:,jhs,initlr,ia) = Mat6(:,je,he,:,js,hs)
                end do
              end do
            end do
          end do
        endif
      endif

      if( E1M1 ) then
        sec_md(:,:) = - secmd_a(:,:,initlr,ie_computer) * cst   ! note the change of sign
        if( isym /= 1 ) call rot_tensor_2( sec_md, matopsym )
   ! Plane and inverse symmetries change the sign
        if( ( isym >= 25 .and. isym <= 48 ) .or. ( isym >= 53 .and. isym <= 57 ) .or. isym == 59  .or. isym == 61 &
                          .or. isym == 63 ) sec_md(:,:) = - sec_md(:,:)
   
   ! Magnetism in on the real part
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_md(:,:) = - conjg( sec_md(:,:) )
        secmdia(:,:,initlr,ia) = sec_md(:,:)
        if( Green_int_mag ) then
          sec_md(:,:) = - secmd_m_a(:,:,initlr,ie_computer) * cst ! note the change of sign
          if( isym /= 1 ) call rot_tensor_2( sec_md, matopsym )
          if( ( isym >= 25 .and. isym <= 48 ) .or. ( isym >= 53 .and. isym <= 57 ) .or. isym == 59  .or. isym == 61 &
                          .or. isym == 63 ) sec_md(:,:) = - sec_md(:,:)
          if( isymeq(ia) < 0 ) sec_md(:,:) = - sec_md(:,:)
          secmdia_m(:,:,initlr,ia) = sec_md(:,:)
        else
          secmdia_m(:,:,initlr,ia) = (0._db,0._db)
        endif
      endif

      if( M1M1 ) then
        sec_mm(:,:) = secmm_a(:,:,initlr,ie_computer) * cst
        if( isym /= 1 ) call rot_tensor_2( sec_mm, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) sec_mm(:,:) = conjg( sec_mm(:,:) )
        secmmia(:,:,initlr,ia) = sec_mm(:,:)
        if( Green_int_mag ) then
          sec_mm(:,:) = secmm_m_a(:,:,initlr,ie_computer) * cst
          if( isym /= 1 ) call rot_tensor_2( sec_mm, matopsym )
          if( isymeq(ia) < 0 ) sec_mm(:,:) = - sec_mm(:,:)
          secmmia_m(:,:,initlr,ia) = sec_mm(:,:)
        endif
      endif

    end do ! end of loop over ia

  end do ! end of loop over initlr

  if( Dafs ) then
    phdt1(:) = phdt(:,1,1)
    if( Bulk_step ) then
      phdf0t1(:) = Sum_Bragg_bulk_abs_f0(:,1)
    else
      phdf0t1(:) = phdf0t(:,1)
    endif
  endif

  if( Spherical_tensor ) call Spherical_tensor_cal(Abs_U_iso,Bragg_abs,ct_nelec,Core_resolved,Volume_maille,E_cut,E1E1,E1E2, &
            E1M1,E2E2,Energ,Ephseuil,Epsii,Eseuil(nbseuil),First_E,i_range,icheck,ie,Int_tens,jseuil,Magn_sens,Moyenne,n_rel, &
            n_tens_max,natomsym,nenerg,ninit1,ninitlr,nomfich_s,nphim,npldafs,nplr,nplrm,nplt,nseuil,numat_abs,pdp,phdf0t1, &
            phdt1,pol,Poldafse,Poldafss,secddia,secdqia,secmdia,secqqia,Spherical_signal,Surface_ref,Taux_eq,V0muf,Vec,Vecdafse, &
            Vecdafss)

  E_vec = E1E2 .or. E2E2 .or. E1E3 .or. E3E3 .or. E1M1 .or. M1M1

  jpl = 0

  if( Cor_abs ) then
    nplt = nplr + 3*npldafs
  else
    nplt = nplr + npldafs
  endif
  
  do ixandafs = 1,2 !  = 1 absorption, = 2 DAFS

    do ipl = 1,nplt ! Loop over polarizations

      mu_cal = .false.
      idafs = .false.

      if( ipl > nplr ) then
        if( Cor_abs ) then
          ipldafs = ( ipl - nplr + 2 ) / 3
          ind_mu = mod(ipl - nplr + 2, 3)
          if( ind_mu == 0 ) then
            idafs = .true.
          else
            mu_cal = .true.
          endif
        else
          ipldafs = ipl - nplr
          idafs = .true.
        endif
      else
        if( ixandafs == 2 ) cycle
        idafs = .false.
        jpl = jpl + 1
        ipldafs = 0
        if( ipl == ncolr - 3 * n_mat_cal / 2 .and. Xan_atom ) jpl = jpl + 1
        if( ipl == ncolr - 3 * n_mat_cal / 2 .and. Moyenne ) jpl = jpl + 1
      endif

      if( idafs .and. ixandafs == 1 ) cycle
      if( .not. idafs .and. ixandafs == 2 ) cycle

      Tens_comp = Magn_sens .or. Green_int .or. idafs

      Tensor_eval = .true.
      if( idafs .and. ipldafs > 1 ) then
        if(     ( abs( hkl_dafs(1,ipldafs) - hkl_dafs(1,ipldafs-1) ) < eps10 ) &
          .and. ( abs( hkl_dafs(2,ipldafs) - hkl_dafs(2,ipldafs-1) ) < eps10 ) &
          .and. ( abs( hkl_dafs(3,ipldafs) - hkl_dafs(3,ipldafs-1) ) < eps10 ) ) Tensor_eval = .false.
      endif
      if( ipl > 1 .and. .not. idafs ) Tensor_eval = .false.

      if( Tensor_eval ) then

        secddia(:,:,:,:,0) = (0._db,0._db)
        secdqia(:,:,:,:,0) = (0._db,0._db)
        secdqia_m(:,:,:,:,0) = (0._db,0._db)
        secqqia(:,:,:,:,:,0) = (0._db,0._db)
        secdoia(:,:,:,:,:,0) = (0._db,0._db)
        if( E3E3 ) secooia(:,:,:,:,:,0) = (0._db,0._db)
        secmdia(:,:,:,0) = (0._db,0._db)
        secmdia_m(:,:,:,0) = (0._db,0._db)
        secmmia(:,:,:,0) = (0._db,0._db)
        if( Green_int_mag ) then
          secddia_m(:,:,:,:,0) = (0._db,0._db)
          secqqia_m(:,:,:,:,:,0) = (0._db,0._db)
          secdoia_m(:,:,:,:,:,0) = (0._db,0._db)
          if( E3E3 ) secooia_m(:,:,:,:,:,0) = (0._db,0._db)
          secmmia_m(:,:,:,0) = (0._db,0._db)
        endif

        do ia = 1,natomsym

          if( idafs ) then
! exp(iQr) is conjugated. One gets back the complex conjugate in convolution.
! Bragg_abs already contains Taux_eq.
            ph = conjg( Bragg_abs(ia,ipldafs) )
          else
            ph = (1._db, 0._db) * Taux_eq(ia)
          endif
          ph_m = img * ph

          if( E1E1 ) secddia(:,:,:,:,0) = secddia(:,:,:,:,0) + ph * secddia(:,:,:,:,ia)
          
          if( E1E2 ) then
            if( Green_int ) then
              secdqia(:,:,:,:,0) = secdqia(:,:,:,:,0) + ph * secdqia(:,:,:,:,ia)
            else
              secdqia(:,:,:,:,0) = secdqia(:,:,:,:,0) + ph * real( secdqia(:,:,:,:,ia), db)
              if( Magn_sens ) secdqia_m(:,:,:,:,0) = secdqia_m(:,:,:,:,0) + ph_m * aimag( secdqia(:,:,:,:,ia) )
            endif
          endif

          if( E2E2 ) secqqia(:,:,:,:,:,0) = secqqia(:,:,:,:,:,0) + ph * secqqia(:,:,:,:,:,ia)

          if( E1E3 ) secdoia(:,:,:,:,:,0) = secdoia(:,:,:,:,:,0) + ph * secdoia(:,:,:,:,:,ia)

          if( E3E3 ) secooia(:,:,:,:,:,0) = secooia(:,:,:,:,:,0) + ph * secooia(:,:,:,:,:,ia)

          if( E1M1 ) then
          ! Magnetism is on real part
            if( Green_int ) then
              secmdia(:,:,:,0) = secmdia(:,:,:,0) + ph * secmdia(:,:,:,ia)
            else
              secmdia(:,:,:,0) = secmdia(:,:,:,0) + ph_m * aimag( secmdia(:,:,:,ia) )
              if( Magn_sens ) secmdia_m(:,:,:,0) = secmdia_m(:,:,:,0) + ph * real( secmdia(:,:,:,ia), db )
            endif
          endif

          if( M1M1 ) secmmia(:,:,:,0) = secmmia(:,:,:,0) + ph * secmmia(:,:,:,ia)

          if( Green_int_mag ) then
            if( E1E1 ) secddia_m(:,:,:,:,0) = secddia_m(:,:,:,:,0) + ph * secddia_m(:,:,:,:,ia)
            if( E1E2 ) secdqia_m(:,:,:,:,0) = secdqia_m(:,:,:,:,0) + ph * secdqia_m(:,:,:,:,ia)
            if( E2E2 ) secqqia_m(:,:,:,:,:,0) = secqqia_m(:,:,:,:,:,0) + ph * secqqia_m(:,:,:,:,:,ia)
            if( E1E3 ) secdoia_m(:,:,:,:,:,0) = secdoia_m(:,:,:,:,:,0) + ph * secdoia_m(:,:,:,:,:,ia)
            if( E3E3 ) secooia_m(:,:,:,:,:,0 ) = secooia_m(:,:,:,:,:,0) + ph * secooia_m(:,:,:,:,:,ia)
            if( E1M1 ) secmdia_m(:,:,:,0) = secmdia_m(:,:,:,0) + ph * secmdia_m(:,:,:,ia)
            if( M1M1 ) secmmia_m(:,:,:,0) = secmmia_m(:,:,:,0) + ph * secmmia_m(:,:,:,ia)
          endif

        end do

        if( Cartesian_tensor ) then
          do ib = 0,nb
            if( natomsym == 1 .and. ib > 0 ) cycle
            if( ib == 0 ) then
              ia = 1
            elseif( ib == 1 ) then
              ia = 0
            else
              ia = ib
            endif
            if( ia /= 0 .and. ( ipl > 1 .and. .not. idafs ) ) cycle
            call Write_cartesian_tensor(Abs_U_iso,Volume_maille,E_cut,E1E2,E2E2,Ephseuil,Epsii,Eseuil(nbseuil),First_E,i_range, &
               ia,ipldafs,jseuil,M1M1,Magn_sens,n_rel,natomsym,ninit1,ninitlr,nomfich_s,nseuil,numat_abs,secddia,secdqia, &
               secdqia_m,secqqia,secmdia,Tens_comp,V0muf,Core_resolved)
          end do
        endif

      endif

      if( .not. Bulk_roughness_case .and. &
          (     ( ie == 1 .and. ( ( icheck > 0 .and. ipl == 1 ) .or. ( icheck > 1 .and. idafs ) )  ) &
           .or. (                 ( icheck > 1 .and. ipl == 1 ) .or. ( icheck > 2 .and. idafs )    )     ) ) &
        call write_ten_bav(Core_resolved, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Final_tddft, Green_int, Green_int_mag, icheck, &
                         idafs, ipl, ipldafs, jseuil, M1M1, Magn_sens, n_oo, n_rel, natomsym, nb, nbseuil, ninit1, ninitlr, &
                         nseuil, secddia, secddia_m, secdoia, secdoia_m, secdqia, secdqia_m, secmdia, secmdia_m, secmmia, &
                         secmmia_m, secooia, secooia_m, secqqia, secqqia_m, Tens_comp)

      if( idafs .or. mu_cal ) then
        np = nphi_dafs(ipldafs)
      else
        np = 1
      endif

      do ip = 1,np

        if( .not. ( idafs .or. mu_cal ) ) then

          if( ltypcal(ipl) == 'sp   ' ) then
            plae(:) = pol(:,ipl-1)
            plas(:) = pol(:,ipl)
          elseif( ltypcal(ipl) == 'ps   ' ) then
            plae(:) = pol(:,ipl+1)
            plas(:) = pol(:,ipl)
          else
            plae(:) = pol(:,ipl)
            plas(:) = pol(:,ipl)
          endif
          if( E_vec ) voae(:) = vec(:,ipl)
          if( E_vec ) voas(:) = vec(:,ipl)

        elseif( idafs ) then

          plae(:) = poldafse(:,ipldafs,ip)
          plas(:) = poldafss(:,ipldafs,ip)
          if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)
          if( E_vec ) voas(:) = vecdafss(:,ipldafs,ip)

        else ! calculation of mu

          if( Full_self_abs ) then

            select case(mod(ipldafs,4))
              case(1,0)
                if( ind_mu == 1 ) then   ! incoming
                  plae(:) = poldafse(:,ipldafs,ip)
                  plas(:) = poldafse(:,ipldafs,ip)
                  if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)
                  if( E_vec ) voas(:) = vecdafse(:,ipldafs,ip)
                else
                  plae(:) = poldafss(:,ipldafs,ip)
                  plas(:) = poldafss(:,ipldafs,ip)
                  if( E_vec ) voae(:) = vecdafss(:,ipldafs,ip)
                  if( E_vec ) voas(:) = vecdafss(:,ipldafs,ip)
                endif
              case(2)
                if( ind_mu == 1 ) then   ! incoming
                  plae(:) = poldafse(:,ipldafs,ip)
                  plas(:) = poldafse(:,ipldafs+1,ip)
                  if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)
                  if( E_vec ) voas(:) = vecdafse(:,ipldafs+1,ip)
                else
                  plae(:) = poldafss(:,ipldafs+1,ip)
                  plas(:) = poldafss(:,ipldafs+2,ip)
                  if( E_vec ) voae(:) = vecdafss(:,ipldafs+1,ip)
                  if( E_vec ) voas(:) = vecdafss(:,ipldafs+2,ip)
                endif
              case(3)
                if( ind_mu == 1 ) then   ! incoming
                  plae(:) = poldafse(:,ipldafs,ip)
                  plas(:) = poldafse(:,ipldafs-1,ip)
                  if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)
                  if( E_vec ) voas(:) = vecdafse(:,ipldafs-1,ip)
                else
                  plae(:) = poldafss(:,ipldafs+1,ip)
                  plas(:) = poldafss(:,ipldafs,ip)
                  if( E_vec ) voae(:) = vecdafss(:,ipldafs+1,ip)
                  if( E_vec ) voas(:) = vecdafss(:,ipldafs,ip)
                endif
            end select

          else  ! Self_abs

            if( ind_mu == 1 ) then
              plae(:) = poldafse(:,ipldafs,ip)
              plas(:) = poldafse(:,ipldafs,ip)
              if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)
              if( E_vec ) voas(:) = vecdafse(:,ipldafs,ip)
            else
              plae(:) = poldafss(:,ipldafs,ip)
              plas(:) = poldafss(:,ipldafs,ip)
              if( E_vec ) voae(:) = vecdafss(:,ipldafs,ip)
              if( E_vec ) voas(:) = vecdafss(:,ipldafs,ip)
            endif

          endif

        endif

        if( idafs .or. mu_cal ) then
! Dans convolution on reprend le complexe conjugue de l'ensemble polarisation x diffusion
          plae(:) = conjg( plae(:) )
          plas(:) = conjg( plas(:) )
        endif

        if( M1M1 .or. E1M1 ) then
          uae(1) = voae(2) * plae(3) - voae(3) * plae(2)
          uae(2) = voae(3) * plae(1) - voae(1) * plae(3)
          uae(3) = voae(1) * plae(2) - voae(2) * plae(1)
          uas(1) = voas(2) * plas(3) - voas(3) * plas(2)
          uas(2) = voas(3) * plas(1) - voas(1) * plas(3)
          uas(3) = voas(1) * plas(2) - voas(2) * plas(1)
        endif

        do ia = 0,nab
          do initlr = 1,ninitlr

            if( Core_resolved .and. .not. Final_tddft ) then
              if( initlr <= ninit1 ) then
                iseuil = 1
              else
                iseuil = min(2, nbseuil)
              endif
            elseif( Final_tddft ) then
              iseuil = min(2, nbseuil)
            else
              iseuil = initlr
            endif

            Ephoton = Energ(ie) + Eseuil(iseuil)
! For very low energy edges (to avoid negative photon energy) :
            Ephoton = max( 0.001_db/Rydb, Ephoton)
            Omega = Ephoton / quatre_mc2 

            if( E1E1 ) then
              sec = (0._db,0._db)
 
              if( .not. Dip_rel ) then
                do ke = 1,3
                  sec = sec + plae(ke) * sum( conjg(plas(:)) * secddia(:,ke,1,initlr,ia))
                end do
                if( Green_int_mag ) then
                  do ke = 1,3
                    sec = sec + plae(ke) * sum( conjg(plas(:)) * secddia_m(:,ke,1,initlr,ia) )
                  end do
                endif
                
              else

                Diag_spin = .false.
                
                ispfg = 0
                do isp1 = 1,2
                  do isp2 = 1,2
                    do isp3 = 1,2
                      do isp4 = 1,2
                        if( n_rel == 8 .and. isp1 /= isp4 ) cycle 
                        ispfg = ispfg + 1

                        if( isp2 /= isp3 .and. Diag_spin ) cycle

                        do ks = 1,3
                          if( isp1 == isp2 ) then
                            matpols = plas(ks)
                          else
                            matpols = (0._db, 0._db)
                          endif
                          if( isp1 == 1 .and. isp2 == 1 ) then
                            if( ks == 1 ) then
                              matpols = matpols - img * Omega * plas(2)
                            elseif( ks == 2 ) then
                              matpols = matpols + img * Omega * plas(1)
                            endif
                          elseif( isp1 == 1 .and. isp2 == 2 ) then
                            if( ks == 1 ) then
                              matpols = matpols + Omega * plas(3)
                            elseif( ks == 2 ) then
                              matpols = matpols - img * Omega * plas(3)
                            else
                              matpols = matpols + img * Omega * ( plas(2) + img * plas(1) )
                            endif
                          elseif( isp1 == 2 .and. isp2 == 1 ) then
                            if( ks == 1 ) then
                              matpols = matpols - Omega * plas(3)
                            elseif( ks == 2 ) then
                              matpols = matpols - img * Omega * plas(3)
                            else
                              matpols = matpols + img * Omega * ( plas(2) - img * plas(1) )
                            endif
                          elseif( isp1 == 2 .and. isp2 == 2 ) then
                            if( ks == 1 ) then
                              matpols = matpols + img * Omega * plas(2)
                            elseif( ks == 2 ) then
                              matpols = matpols - img * Omega * plas(1)
                            endif
                          endif                      

                          do ke = 1,3
                            if( isp3 == isp4 ) then
                              matpole = plae(ke)
                            else
                              matpole = (0._db, 0._db)
                            endif
                            if( isp3 == 1 .and. isp4 == 1 ) then
                              if( ke == 1 ) then
                                matpole = matpole - img * Omega * plae(2)
                              elseif( ke == 2 ) then
                                matpole = matpole + img * Omega * plae(1)
                              endif
                            elseif( isp3 == 1 .and. isp4 == 2 ) then
                              if( ke == 1 ) then
                                matpole = matpole + Omega * plae(3)
                              elseif( ke == 2 ) then
                                matpole = matpole - img * Omega * plae(3)
                              else
                                matpole = matpole + img * Omega * ( plae(2) + img * plae(1) )
                              endif
                            elseif( isp3 == 2 .and. isp4 == 1 ) then
                              if( ke == 1 ) then
                                matpole = matpole - Omega * plae(3)
                              elseif( ke == 2 ) then
                                matpole = matpole - img * Omega * plae(3)
                              else
                                matpole = matpole + img * Omega * ( plae(2) - img * plae(1) )
                              endif
                            elseif( isp3 == 2 .and. isp4 == 2 ) then
                              if( ke == 1 ) then
                                matpole = matpole + img * Omega * plae(2)
                              elseif( ke == 2 ) then
                                matpole = matpole - img * Omega * plae(1)
                              endif
                            endif

                            sec = sec + conjg( matpols ) * matpole * secddia(ks,ke,ispfg,initlr,ia)

                            if( Green_int_mag ) sec = sec + conjg( matpols ) * matpole * secddia_m(ks,ke,ispfg,initlr,ia)
                          
                          end do
                        end do

                      end do
                    end do
                  end do
                end do

              endif
              
! A factor pi is missing, it has aleady been taken into account in routine Tens_ab.
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                Ampldafs_dd(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                mu_dd(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                Secabs_dd(jpl,initlr,ia) = real( sec, db )
                if( ltypcal(ipl) == 'sp   ' .or. ltypcal(ipl) == 'ps   ' ) Secabs_dd(jpl+1,initlr,ia) = aimag( sec )
              endif
            endif

            if( E1E2 ) then
              sec = (0._db,0._db)
              do ke = 1,3
                do ks = 1,3
                  sec = sec + conjg( plas(ks) ) * plae(ke) &
                            * sum( voae(:) * secdqia(ks,ke,:,initlr,ia) - voas(:) * secdqia(ke,ks,:,initlr,ia) )
                  if( Magn_sens ) sec = sec + conjg( plas(ks) ) * plae(ke) &
                            * sum( voae(:) * secdqia_m(ks,ke,:,initlr,ia) + voas(:) * secdqia_m(ke,ks,:,initlr,ia) )
                end do
              end do

              sec = img * sec      ! Here one gets back img
              if( idafs ) sec = - sec

              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                Ampldafs_dq(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                mu_dq(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                Secabs_dq(jpl,initlr,ia) = real( sec, db )
                if( ltypcal(ipl) == 'sp   ' .or. ltypcal(ipl) == 'ps   ' ) Secabs_dq(jpl+1,initlr,ia) = aimag( sec )
              endif
            endif

            if( E2E2 ) then
              sec = (0._db,0._db)
              do ke = 1,3
                do je = 1,3
                  do ks = 1,3
                    sec = sec + conjg( plas(ks) ) * plae(ke) * voae(je) * sum( voas(:) * secqqia(ks,:,ke,je,initlr,ia) )
                    if( Green_int_mag ) &
                      sec = sec + conjg( plas(ks) ) * plae(ke) * voae(je) * sum( voas(:) * secqqia_m(ks,:,ke,je,initlr,ia) )
                  end do
                end do
              end do
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                Ampldafs_qq(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                mu_qq(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                Secabs_qq(jpl,initlr,ia) = real( sec,db )
                if( ltypcal(ipl) == 'sp   ' .or. ltypcal(ipl) == 'ps   ' ) Secabs_qq(jpl+1,initlr,ia) = aimag( sec )
              endif
            endif

            if( E1E3 ) then
              sec = (0._db,0._db)
              do ke = 1,3
                do ks = 1,3
                  do j1 = 1,3
                    sec = sec + conjg( plas(ks) ) * plae(ke) * ( voae(j1) * sum( voae(:) * secdoia(ks,ke,j1,:,initlr,ia) ) &
                        + voas(j1) * sum( voas(:) * conjg( secdoia(ke,ks,j1,:,initlr,ia) ) ) )
                    if( Green_int_mag ) &
                      sec = sec + conjg( plas(ks) ) * plae(ke) * ( voae(j1) * sum( voae(:) &
                          * secdoia_m(ks,ke,j1,:,initlr,ia) ) + voas(j1) * sum( voas(:) &
                          * conjg( secdoia_m(ke,ks,j1,:,initlr,ia) )))
                  end do
                end do
              end do
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                Ampldafs_do(ipldafs,ip,initlr,ia)=sec
              elseif( mu_cal ) then
                mu_do(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                Secabs_do(jpl,initlr,ia) = real( sec, db )
                if( ltypcal(ipl) == 'sp   ' .or. ltypcal(ipl) == 'ps   ' ) Secabs_do(jpl+1,initlr,ia) = aimag( sec )
              endif
            endif

            if( E3E3 ) then
              sec = (0._db,0._db)
              do ke = 1,3
                do je = 1,3
                  do he = 1,3
                    jhe = 3 * ( je - 1 ) + he
                    do ks = 1,3
                      do js = 1,3
                        do hs = 1,3
                          jhs = 3 * ( js - 1 ) + hs
                          sec = sec +conjg( plas(ks) ) * voas(js) * voas(hs) * plae(ke) * voae(je) * voae(he) &
                            * secooia(ks,jhs,ke,jhe,initlr,ia)
                          if( Green_int_mag ) &
                            sec = sec + conjg( plas(ks) ) * voas(js) * voas(hs) * plae(ke) * voae(je) * voae(he) &
                                      * secooia_m(ks,jhs,ke,jhe,initlr,ia)
                        end do
                      end do
                    end do
                  end do
                end do
              end do
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                Ampldafs_oo(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                mu_oo(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                Secabs_oo(jpl,initlr,ia) = real( sec, db )
                if( ltypcal(ipl) == 'sp   ' .or. ltypcal(ipl) == 'ps   ' ) Secabs_oo(jpl+1,initlr,ia) = aimag( sec )
              endif
            endif

            if( M1M1 ) then
              sec = ( 0._db, 0._db )
              
              do ke = 1,3
                sec = sec + uae(ke) * sum( conjg(uas(:)) * secmmia(:,ke,initlr,ia) )
               if( Green_int_mag ) &
                 sec = sec + uae(ke) * sum( conjg(uas(:)) * secmmia_m(:,ke,initlr,ia) )
              end do
! A factor pi is missing, it has aleady been taken into account in routine Tens_ab.
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                Ampldafs_mm(ipldafs,ip,initlr,ia)=sec
              elseif( mu_cal ) then
                mu_mm(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                Secabs_mm(jpl,initlr,ia) = real( sec, db )
                if( ltypcal(ipl) == 'sp   ' .or. ltypcal(ipl) == 'ps   ' ) Secabs_mm(jpl+1,initlr,ia) = aimag( sec )
              endif
            endif

            if( E1M1 ) then
              sec = (0._db,0._db)
              do ke = 1,3
                if( ia == 0 ) then
  ! on met le signe moins devant uae car la partie imaginaire contient le terme non-magnetique, equivaut au complexe conjugue
  ! sauf si secmdia est calcule pour une reflexion dafs ( la partie reelle est nulle sauf dafs )
                  sec = sec + conjg( uas(ke) ) * sum( plae(:) * secmdia(:,ke,initlr,ia) ) &
                            - uae(ke) * sum( conjg(plas(:)) * secmdia(:,ke,initlr,ia) )  
  ! on met le signe plus devant uae car la partie reelle contient le terme magnetique, (le complexe conjugue est donc inutile)
                  if( Magn_sens ) sec = sec + conjg( uas(ke) ) * sum( plae(:) * secmdia_m(:,ke,initlr,ia) ) &
                                                + uae(ke) * sum( conjg( plas(:) ) * secmdia_m(:,ke,initlr,ia) )
                else
                  sec = sec + conjg( uas(ke) ) * sum( plae(:) * secmdia(:,ke,initlr,ia) ) &
                          + uae(ke) * sum( conjg(plas(:)) * conjg( secmdia(:,ke,initlr,ia) ) )  
               endif 
              end do

! A factor pi is missing, it has aleady been taken into account in routine Tens_ab.
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                Ampldafs_md(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                mu_md(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                Secabs_md(jpl,initlr,ia) = real( sec, db )
                if( ltypcal(ipl) == 'sp   ' .or. ltypcal(ipl) == 'ps   ' ) Secabs_md(jpl+1,initlr,ia) = aimag( sec )
              endif
            endif

          end do   ! End of loop over core states (initlr)
        end do   ! End of loop over atoms

        if( ipl > nplr ) cycle
        if( ltypcal(ipl) == 'sp   ' .or. ltypcal(ipl) == 'ps   ') jpl = jpl + 1

      end do

    end do   ! End of loop over polarizations
  end do

  if( Moyenne ) then
    i = ncolr - 3 * n_mat_cal / 2
    if( Xan_atom ) i = i - 1
    ipl = 0
    do j = 1,i
      ipl = ipl + 1
      if( E1E1 ) Secabs_dd(i,:,:) = Secabs_dd(i,:,:) + pdp(ipl,1) * Secabs_dd(j,:,:)
      if( E2E2 ) Secabs_qq(i,:,:) = Secabs_qq(i,:,:) + pdp(ipl,2) * Secabs_qq(j,:,:)
      if( E1E3 ) Secabs_do(i,:,:) = Secabs_do(i,:,:) + pdp(ipl,1) * Secabs_do(j,:,:)
      if( E3E3 ) Secabs_oo(i,:,:) = Secabs_oo(i,:,:) + pdp(ipl,1) * Secabs_oo(j,:,:)
      if( M1M1 ) Secabs_mm(i,:,:) = Secabs_mm(i,:,:) + pdp(ipl,1) * Secabs_mm(j,:,:)
      if( E1M1 ) Secabs_md(i,:,:) = Secabs_md(i,:,:) + pdp(ipl,1) * Secabs_md(j,:,:)
      if( ipl == nplr ) exit
    end do
    if( E1E2 ) Secabs_dq(i,:,:) = (0._db,0._db)
  endif

  do ipl = 1,nplr
    if( ltypcal(ipl)(1:4) /= 'left' ) cycle
    allocate( Abs_sum(ninitlr,0:natomsym) )
    if( E1E1 ) then
      Abs_sum(:,:) = 0.5_db * ( Secabs_dd(ipl,:,:) + Secabs_dd(ipl+1,:,:) )
      Secabs_dd(ipl+1,:,:) = Secabs_dd(ipl,:,:) - Secabs_dd(ipl+1,:,:)
      Secabs_dd(ipl,:,:) = Abs_sum(:,:)
    endif
    if( E1E2 ) then
      Abs_sum(:,:) = 0.5_db * ( Secabs_dq(ipl,:,:) + Secabs_dq(ipl+1,:,:) )
      Secabs_dq(ipl+1,:,:) = Secabs_dq(ipl,:,:) - Secabs_dq(ipl+1,:,:)
      Secabs_dq(ipl,:,:) = Abs_sum(:,:)
    endif
    if( E2E2 ) then
      Abs_sum(:,:) = 0.5_db * ( Secabs_qq(ipl,:,:) + Secabs_qq(ipl+1,:,:) )
      Secabs_qq(ipl+1,:,:) = Secabs_qq(ipl,:,:) - Secabs_qq(ipl+1,:,:)
      Secabs_qq(ipl,:,:) = Abs_sum(:,:)
    endif
    if( E1E3 ) then
      Abs_sum(:,:) = 0.5_db * ( Secabs_do(ipl,:,:) + Secabs_do(ipl+1,:,:) )
      Secabs_do(ipl+1,:,:) = Secabs_do(ipl,:,:) - Secabs_do(ipl+1,:,:)
      Secabs_do(ipl,:,:) = Abs_sum(:,:)
    endif
    if( E3E3 ) then
      Abs_sum(:,:) = 0.5_db * ( Secabs_oo(ipl,:,:) + Secabs_oo(ipl+1,:,:) )
      Secabs_oo(ipl+1,:,:) = Secabs_oo(ipl,:,:) - Secabs_oo(ipl+1,:,:)
      Secabs_oo(ipl,:,:) = Abs_sum(:,:)
    endif
    if( M1M1 ) then
      Abs_sum(:,:) = 0.5_db * ( Secabs_mm(ipl,:,:) + Secabs_mm(ipl+1,:,:) )
      Secabs_mm(ipl+1,:,:) = Secabs_mm(ipl,:,:) - Secabs_mm(ipl+1,:,:)
      Secabs_mm(ipl,:,:) = Abs_sum(:,:)
    endif
    if( E1M1 ) then
      Abs_sum(:,:) = 0.5_db * ( Secabs_md(ipl,:,:) + Secabs_md(ipl+1,:,:) )
      Secabs_md(ipl+1,:,:) = Secabs_md(ipl,:,:) - Secabs_md(ipl+1,:,:)
      Secabs_md(ipl,:,:) = Abs_sum(:,:)
    endif
    deallocate( Abs_sum )
  end do

  if( xan_atom ) then
    Secabs_dd(ncolr-3*n_mat_cal/2,:,0) = sec_atom(:) * natomsym_f * cst
    do ia = 1,nab
      Secabs_dd(ncolr-3*n_mat_cal/2,:,ia) = sec_atom(:) * Taux_eq(ia) * cst
    end do
  endif

  Secabs(:,:,:) = 0._db
  if( E1E1 ) Secabs(:,:,:) = Secabs(:,:,:) + Secabs_dd(:,:,:)
  if( E1E2 ) Secabs(:,:,:) = Secabs(:,:,:) + Secabs_dq(:,:,:)
  if( E2E2 ) Secabs(:,:,:) = Secabs(:,:,:) + Secabs_qq(:,:,:)
  if( E1E3 ) Secabs(:,:,:) = Secabs(:,:,:) + Secabs_do(:,:,:)
  if( E3E3 ) Secabs(:,:,:) = Secabs(:,:,:) + Secabs_oo(:,:,:)
  if( M1M1 ) Secabs(:,:,:) = Secabs(:,:,:) + Secabs_mm(:,:,:)
  if( E1M1 ) Secabs(:,:,:) = Secabs(:,:,:) + Secabs_md(:,:,:)

  if( Dafs ) then
    Ampldafs(:,:,:,:) = (0._db,0._db)
    if( E1E1 ) Ampldafs(:,:,:,:) = Ampldafs(:,:,:,:) + Ampldafs_dd(:,:,:,:)
    if( E1E2 ) Ampldafs(:,:,:,:) = Ampldafs(:,:,:,:) + Ampldafs_dq(:,:,:,:)
    if( E2E2 ) Ampldafs(:,:,:,:) = Ampldafs(:,:,:,:) + Ampldafs_qq(:,:,:,:)
    if( E1E3 ) Ampldafs(:,:,:,:) = Ampldafs(:,:,:,:) + Ampldafs_do(:,:,:,:)
    if( E3E3 ) Ampldafs(:,:,:,:) = Ampldafs(:,:,:,:) + Ampldafs_oo(:,:,:,:)
    if( M1M1 ) Ampldafs(:,:,:,:) = Ampldafs(:,:,:,:) + Ampldafs_mm(:,:,:,:)
    if( E1M1 ) Ampldafs(:,:,:,:) = Ampldafs(:,:,:,:) + Ampldafs_md(:,:,:,:)
  endif

  if( Cor_abs ) then
    mu(:,:,:,:,:) = (0._db,0._db)
    if( E1E1 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mu_dd(:,:,:,:,:)
    if( E1E2 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mu_dq(:,:,:,:,:)
    if( E2E2 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mu_qq(:,:,:,:,:)
    if( E1E3 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mu_do(:,:,:,:,:)
    if( E3E3 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mu_oo(:,:,:,:,:)
    if( M1M1 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mu_mm(:,:,:,:,:)
    if( E1M1 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mu_md(:,:,:,:,:)
  endif

! Conversion in number of electrons
  if( dafs ) then
    do initlr = 1,ninitlr
      cst = ct_nelec(initlr) / pi 
      Ampldafs(:,:,initlr,:) = cst * Ampldafs(:,:,initlr,:)
      if( E1E1 ) Ampldafs_dd(:,:,initlr,:) = cst * Ampldafs_dd(:,:,initlr,:)
      if( E1E2 ) Ampldafs_dq(:,:,initlr,:) = cst * Ampldafs_dq(:,:,initlr,:)
      if( E2E2 ) Ampldafs_qq(:,:,initlr,:) = cst * Ampldafs_qq(:,:,initlr,:)
      if( M1M1 ) Ampldafs_mm(:,:,initlr,:) = cst * Ampldafs_mm(:,:,initlr,:)
      if( E1M1 ) Ampldafs_md(:,:,initlr,:) = cst * Ampldafs_md(:,:,initlr,:)
    end do
  endif

! Conversion in micrometers
  if( Cor_abs ) then
    c_micro = 100 / ( Volume_maille * bohr**3 )
    mu(:,:,:,:,:) = c_micro * mu(:,:,:,:,:)
    if( E1E1 ) mu_dd(:,:,:,:,:) = c_micro * mu_dd(:,:,:,:,:)
    if( E1E2 ) mu_dq(:,:,:,:,:) = c_micro * mu_dq(:,:,:,:,:)
    if( E2E2 ) mu_qq(:,:,:,:,:) = c_micro * mu_qq(:,:,:,:,:)
    if( E1E3 ) mu_do(:,:,:,:,:) = c_micro * mu_do(:,:,:,:,:)
    if( E3E3 ) mu_oo(:,:,:,:,:) = c_micro * mu_oo(:,:,:,:,:)
    if( M1M1 ) mu_mm(:,:,:,:,:) = c_micro * mu_mm(:,:,:,:,:)
    if( E1M1 ) mu_md(:,:,:,:,:) = c_micro * mu_md(:,:,:,:,:)
  endif

  if( n_mat_cal > 0 .and. Matper .and. .not. nseuil == 0 ) then
    c_micro = 100 / ( Volume_maille * bohr**3 )
    i = ncolr - 3 * n_mat_cal / 2 + 1
    Secabs(i:ncolr,:,:) = c_micro * Secabs(i:ncolr,:,:)
    if( E1E1 ) Secabs_dd(i:ncolr,:,:) = c_micro * Secabs_dd(i:ncolr,:,:)
    if( E1E2 ) Secabs_dq(i:ncolr,:,:) = c_micro * Secabs_dq(i:ncolr,:,:)
    if( E2E2 ) Secabs_qq(i:ncolr,:,:) = c_micro * Secabs_qq(i:ncolr,:,:)
    if( E1E3 ) Secabs_do(i:ncolr,:,:) = c_micro * Secabs_do(i:ncolr,:,:)
    if( E3E3 ) Secabs_oo(i:ncolr,:,:) = c_micro * Secabs_oo(i:ncolr,:,:)
    if( M1M1 ) Secabs_mm(i:ncolr,:,:) = c_micro * Secabs_mm(i:ncolr,:,:)
    if( E1M1 ) Secabs_md(i:ncolr,:,:) = c_micro * Secabs_md(i:ncolr,:,:)
  endif

 ! Optic case: output in millimeter^-1, possible only for periodical system
  if( nseuil == 0 .and. Matper ) then
    c_milli = 100000 / ( Volume_maille * bohr**3 )
    Secabs(:,:,:) = c_milli * Secabs(:,:,:)
    if( E1E1 ) Secabs_dd(:,:,:) = c_milli * Secabs_dd(:,:,:)
    if( E1E2 ) Secabs_dq(:,:,:) = c_milli * Secabs_dq(:,:,:)
    if( E2E2 ) Secabs_qq(:,:,:) = c_milli * Secabs_qq(:,:,:)
    if( E1E3 ) Secabs_do(:,:,:) = c_milli * Secabs_do(:,:,:)
    if( E3E3 ) Secabs_oo(:,:,:) = c_milli * Secabs_oo(:,:,:)
    if( M1M1 ) Secabs_mm(:,:,:) = c_milli * Secabs_mm(:,:,:)
    if( E1M1 ) Secabs_md(:,:,:) = c_milli * Secabs_md(:,:,:)
  endif

! Writing ----------------------------------------------------------------------------------------------------

  if( icheck > 0 ) then
    if( Bulk_roughness_case ) write(3,'(/A)') ' Bulk roughness case'
    if( Dip_rel ) write(3,'(/1x,a7,1p,e12.5)') 'Omega =', Omega
    do ia = 0,nab
      if( ia == 0 ) then
        if( .not. ( nseuil == 0 .and. Matper ) ) then
          write(3,120) ct_nelec(:)
          if( n_mat_cal > 0 .and. Matper ) write(3,130) c_micro 
        else
          write(3,140) c_milli
       endif
        write(3,150)
      else
        write(3,290) ia
      endif
      do ipl = 1,ncolt
        nomab = nomabs(ipl)
        call center_word(nomab,Length_word)
        nomabs(ipl) = nomab
      end do
      nccm = 360
      nl = 1 + ( ncolr - 1 ) / nccm

      do initlr = 1,ninitlr

        if( ninitlr > 1 ) then
          if( nseuil /= 0 ) then
            write(3,295) initlr
          elseif( initlr > 1 ) then
            write(3,297) initlr - 2
          endif
        endif

        do i = 1,nl
          ic1 = 1 + ( i - 1 ) * nccm
          ic2 = min( i * nccm, ncolr )
          write(3,300) nomabs(ic1:ic2)
          write(3,310) Ephseuil*rydb, Secabs(ic1:ic2,initlr,ia)
          if( E1E1 .and. ( E1E2 .or. E2E2 .or. E1E3 .or. M1M1 .or. E3E3 .or. E1M1) ) write(3,320) Secabs_dd(ic1:ic2,initlr,ia)
          if( E1E2 .and. ( E1E1 .or. E2E2 .or. E1E3 .or. M1M1 .or. E3E3 .or. E1M1) ) write(3,330) Secabs_dq(ic1:ic2,initlr,ia)
          if( E2E2 .and. ( E1E1 .or. E2E2 .or. E1E3 .or. M1M1 .or. E3E3 .or. E1M1) ) write(3,340) Secabs_qq(ic1:ic2,initlr,ia)
          if( E1E3 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. M1M1 .or. E3E3 .or. E1M1) ) write(3,350) Secabs_do(ic1:ic2,initlr,ia)
          if( E3E3 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. M1M1 .or. E1E3 .or. E1M1) ) write(3,351) Secabs_oo(ic1:ic2,initlr,ia)
          if( M1M1 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. E1E3 .or. E3E3 .or. E1M1) ) write(3,352) Secabs_mm(ic1:ic2,initlr,ia)
          if( E1M1 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. E1E3 .or. E3E3 .or. M1M1) ) write(3,354) Secabs_md(ic1:ic2,initlr,ia)
        end do
        if( dafs ) then
          if( self_abs ) then
            nc = 4
          elseif( Full_self_abs ) then
            nc = 6
          else
            nc = 2
          endif
          nl = 1 + ( nc * npldafs - 1 ) / nccm
          do i = 1,nl
            icn1 = 1 + ( i - 1 ) * nccm
            icn2 = min( i * nccm, nc * npldafs )
            ic1 = 1 + ( i - 1 ) * (nccm/nc)
            ic2 = min( i * (nccm/nc), npldafs )
            write(3,360) nomabs(ncolr+icn1:ncolr+icn2)
            if( self_abs ) then
              write(3,370) (Ampldafs(j,1,initlr,ia), Real(mu(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E1E1 .and. ( E1E2 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,320) (Ampldafs_dd(j,1,initlr,ia), Real(mu_dd(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E1E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,330) (Ampldafs_dq(j,1,initlr,ia), Real(mu_dq(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E2E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,340) (Ampldafs_qq(j,1,initlr,ia), Real(mu_qq(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E1E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) (Ampldafs_do(j,1,initlr,ia), Real(mu_do(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E3E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E1E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) (Ampldafs_oo(j,1,initlr,ia), Real(mu_do(j,1,:,initlr,ia)), j = ic1,ic2)
              if( M1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. E1M1 ) ) &
                write(3,352) (Ampldafs_mm(j,1,initlr,ia), Real(mu_mm(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. M1M1 ) ) &
                write(3,354) (Ampldafs_md(j,1,initlr,ia), Real(mu_md(j,1,:,initlr,ia)), j = ic1,ic2)
            elseif(Full_self_abs ) then
              write(3,370) (Ampldafs(j,1,initlr,ia), mu(j,1,:,initlr,ia), j = ic1,ic2)
              if( E1E1 .and. ( E1E2 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,320) (Ampldafs_dd(j,1,initlr,ia), mu_dd(j,1,:,initlr,ia), j = ic1,ic2)
              if( E1E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,330) (Ampldafs_dq(j,1,initlr,ia), mu_dq(j,1,:,initlr,ia), j = ic1,ic2)
              if( E2E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,340) (Ampldafs_qq(j,1,initlr,ia), mu_qq(j,1,:,initlr,ia), j = ic1,ic2)
              if( E1E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) (Ampldafs_do(j,1,initlr,ia), mu_do(j,1,:,initlr,ia), j = ic1,ic2)
              if( E3E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E1E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) (Ampldafs_oo(j,1,initlr,ia), mu_do(j,1,:,initlr,ia), j = ic1,ic2)
              if( M1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. E1M1 ) ) &
                write(3,352) (Ampldafs_mm(j,1,initlr,ia), mu_mm(j,1,:,initlr,ia), j = ic1,ic2)
              if( E1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. M1M1 ) ) &
                write(3,354) (Ampldafs_md(j,1,initlr,ia), mu_md(j,1,:,initlr,ia), j = ic1,ic2)
            else
              write(3,370) Ampldafs(ic1:ic2,1,initlr,ia)
              if( E1E1 .and. ( E1E2 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,320) Ampldafs_dd(ic1:ic2,1,initlr,ia)
              if( E1E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,330) Ampldafs_dq(ic1:ic2,1,initlr,ia)
              if( E2E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,340) Ampldafs_qq(ic1:ic2,1,initlr,ia)
              if( E1E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) Ampldafs_do(ic1:ic2,1,initlr,ia)
              if( E3E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E1E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) Ampldafs_oo(ic1:ic2,1,initlr,ia)
              if( M1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. E1M1 ) ) &
                write(3,352) Ampldafs_mm(ic1:ic2,1,initlr,ia)
              if( E1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. M1M1 ) ) &
                write(3,354) Ampldafs_md(ic1:ic2,1,initlr,ia)
            endif
          end do
        endif
      end do
    end do
  endif

  if( ie == 1 ) write(6,392) nenerg
  n_dim = ncolm * ninitlr

  do ia = 0,na

    i_bulk_z = 0
    do
      i_bulk_z = i_bulk_z + 1
    
    nomficht = nomfich
    nomfichdafst = nomfich

    if( ia > 0 ) then
      long = len_trim(nomficht)
      nomficht(long+1:long+5) = '_atom'
      call ad_number(ia,nomficht,Length_name)
      nomfichdafst(long+1:long+5) = '_atom'
      call ad_number(ia,nomfichdafst,Length_name)
    endif
    long = len_trim(nomficht)

    if( Final_tddft .and. .not. Extract_ten ) then
      nomficht(long+1:long+6) = '_tddft'
      nomfichdafst(long+1:long+11) = '_tddft_scan'
    else
      nomfichdafst(long+1:long+5) = '_scan'
    end if

    if( n_multi_run > 1 ) then
      l = len_trim(nomficht)
      nomficht(l+1:l+1) = '_'
      call ad_number(iabsorig,nomficht,Length_name)
      l = len_trim(nomfichdafst)
      nomfichdafst(l+1:l+1) = '_'
      call ad_number(iabsorig,nomfichdafst,Length_name)
    endif
    if( n_bulk_z_abs > 1 .and. .not. Bulk_roughness_case ) then
      l = len_trim(nomficht)
      nomficht(l+1:l+1) = '_'
      call ad_number(i_bulk_z,nomficht,Length_name)
      l = len_trim(nomfichdafst)
      nomfichdafst(l+1:l+1) = '_'
      call ad_number(i_bulk_z,nomfichdafst,Length_name)
    elseif( Bulk_roughness_case ) then ! bulk roughness case
      l = len_trim(nomficht)
      nomficht(l+1:l+2) = '_0'
    endif
    l = len_trim(nomficht)
    nomficht(l+1:l+4) = '.txt'
    l = len_trim(nomfichdafst)
    nomfichdafst(l+1:l+4) = '.txt'

    if( ie == 1 .and. ia == 0 ) nomfich_cal_conv(multi_0+i_bulk_z) = nomficht

    n_tens = 0

    do initlr = 1,ninitlr

      n1 = ( initlr - 1 ) * ( ncolt-nxanout+1 ) + 1
      n2 = initlr * ( ncolt-nxanout+1 )
      title(n1:n2) = nomabs(nxanout:ncolt)

      if( ninitlr > 1 .and. .not. (initlr == 1 .and. nseuil == 0 ) ) then
        do i = n1, n2
          nomab = adjustl(title(i))
          ll = len_trim( nomab )
          if( ll > length_word - 3 ) cycle
          if( nomab(ll:ll) /= '>' .and. n2-n1+1 /= ncolr ) cycle
          ll = ll + 1
          nomab(ll:ll) = '_'
          if( .not. Core_resolved ) then
            ll = ll + 1
            nomab(ll:ll) = achar(nseuil+74)
            ll = ll + 1
            nomab(ll:ll) = achar(jseuil+initlr+47)
          else
            if( nseuil == 0 ) then
              call ad_number(initlr-2,nomab,length_word)
            else
              call ad_number(initlr,nomab,length_word)
            endif
          endif
          call center_word(nomab,Length_word)
          title(i) = nomab
        end do
      endif

      ipl = n_tens + ncolr - nxanout + 1
      if( ia == 0 .and. n_bulk_z_abs > 1 .and. .not. Bulk_roughness_case ) then
        Tens(n_tens+1:ipl) = 0._db
        do ib = 1,natomsym
          do i = 1,n_bulk_zc_abs(i_bulk_z)
            if( igr_bulk_z_abs(i,i_bulk_z) == ib ) Tens(n_tens+1:ipl) = Tens(n_tens+1:ipl) + Secabs(nxanout:ncolr,initlr,ib)  
          end do
        end do
      else
        Tens(n_tens+1:ipl) = Secabs(nxanout:ncolr,initlr,ia)
      endif
      do ipldafs = 1,npldafs
        if( ia == 0 .and. ( n_bulk_z_abs <= 1 .or. Bulk_roughness_case ) ) then
          cf = Ampldafs(ipldafs,1,initlr,ia)
        elseif( ia == 0 ) then
          cf = ( 0._db, 0._db )
          do ib = 1,natomsym
            do i = 1,n_bulk_zc_abs(i_bulk_z)
              if( igr_bulk_z_abs(i,i_bulk_z) == ib ) cf = cf + conjg( Bragg_abs(ib,ipldafs) ) * Ampldafs(ipldafs,1,initlr,ib)  
            end do
          end do
        else
! Le exp(iQr) est converti. On recupere le complexe conjugue dans convolution.
          cf = conjg( Bragg_abs(ia,ipldafs) ) * Ampldafs(ipldafs,1,initlr,ia)
        endif
        ipl = ipl + 1
        Tens(ipl) = real( cf,db )
        ipl = ipl + 1
        Tens(ipl) = aimag( cf )
        if( Full_self_abs ) then
          do i = 1,2
            ipl = ipl + 1
            Tens(ipl) = real( mu(ipldafs,1,i,initlr,ia), db )
            ipl = ipl + 1
            Tens(ipl) = aimag( mu(ipldafs,1,i,initlr,ia) )
          end do
        elseif( self_abs ) then
          do i = 1,2
            ipl = ipl + 1
            Tens(ipl) = real( mu(ipldafs,1,i,initlr,ia), db )
          end do
        endif
      end do

      n_tens = ipl

    end do

    if( ia == 0 ) then
      phdtem(:) = phdt(:,1,i_bulk_z)
    else
      phdtem(:) = Bragg_abs(ia,:)
    endif

! Writing in the output file
    if( Full_self_abs .or. Self_abs ) then
      call write_out(Abs_U_iso,angxyz,axyz,f_avantseuil,E_cut,Ephseuil,Epsii,Eseuil(nbseuil),First_E,Green_int,hkl_dafs, &
            i_range,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title,npldafs,npldafs,3,npldafs,nseuil,numat_abs, &
            phdtem,phdf0t1,tens,v0muf,Core_resolved,natomsym_f,Surface_ref,Volume_maille)
    elseif( Bulk_step .and. ie == 1 ) then
      phdf0t1(:) = Sum_Bragg_bulk_abs_f0(:,i_bulk_z)
      if( i_bulk_z == 1 ) then
        npps = (n_bulk_z + 1) * npldafs 
      else
        npps = npldafs
      endif
      allocate( hkl_dafs_fake(3,npps) )
      if( npldafs > 0 ) then
        hkl_dafs_fake(1,1:npldafs) = Length_abs(1:npldafs) 
        hkl_dafs_fake(2,1:npldafs) = Length_rel_abs(i_bulk_z) * Length_abs(1:npldafs)   ! fake value used in write_out, signature of truncation 
        hkl_dafs_fake(2,1) = hkl_dafs_fake(2,1) + 10000._db   ! fake value used in write_out, signature of truncation 
        hkl_dafs_fake(3,1:npldafs) = hkl_dafs(3,1:npldafs)
        if( i_bulk_z == 1 ) then
          do i = 1,n_bulk_z
            hkl_dafs_fake(1,i*npldafs+1:i*npldafs+npldafs) = Length_rel(i) * Length_abs(1:npldafs) 
            hkl_dafs_fake(2,i*npldafs+1:i*npldafs+npldafs) = real( Sum_Bragg_nonabs_f(1:npldafs,i), db ) 
            hkl_dafs_fake(3,i*npldafs+1:i*npldafs+npldafs) = aimag( Sum_Bragg_nonabs_f(1:npldafs,i) ) 
          end do
        endif
      endif 
      call write_out(Abs_U_iso,rdum,rdum,f_avantseuil,E_cut,Ephseuil,Epsii,Eseuil(nbseuil),First_E,Green_int,hkl_dafs_fake, &
            i_range,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title,npldafs,npldafs,0,npps,nseuil,numat_abs, &
            phdtem,phdf0t1,tens,v0muf,Core_resolved,natomsym_f,Surface_ref,Volume_maille)
      deallocate( hkl_dafs_fake )
    else
      call write_out(Abs_U_iso,rdum,rdum,f_avantseuil,E_cut,Ephseuil,Epsii,Eseuil(nbseuil),First_E,Green_int,rdum, &
            i_range,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title,npldafs,npldafs,0,0,nseuil,numat_abs, &
            phdtem,phdf0t1,tens,v0muf,Core_resolved,natomsym_f,Surface_ref,Volume_maille)

    endif

! Writing on the screen
    if( ia == 0 ) call write_out(Abs_U_iso,rdum,rdum,f_avantseuil,E_cut,Ephseuil,Epsii,Eseuil(nbseuil),First_E,Green_int,rdum, &
            i_range,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title,npldafs,npldafs,0,0,nseuil,-1, &
            phdtem,phdf0t1,tens,v0muf,Core_resolved,natomsym_f,Surface_ref,Volume_maille)

    if( Dafs .and. nphim > 1 ) then
      if( ie == 1 ) then
        open(7, file = nomfichdafst)
        do ipl = 1,npldafs
          write(7,400) nphi_dafs(ipl)
        end do
      else
        open(7, file = nomfichdafst, position='append')
      endif
      mot = ' '
      i = 6
      do initlr = 1,ninitlr
        do ii = 1,3
          select case(ii)
            case(1)
              nomab = 'Amplitude_'
            case(2)
              if( .not. ( Full_self_abs .or. Self_abs ) ) cycle
              if( Full_self_abs ) then
                nomab = 'mu_in_'
              else
                nomab = 'mu_io_'
              endif
            case(3)
              if( .not. Full_self_abs ) cycle
              nomab = 'mu_ou_'
          end select
          if( .not. Core_resolved ) then
            nomab(11:11) = achar(nseuil+74)
            nomab(12:12) = achar(jseuil+initlr+47)
          else
            call ad_number(initlr,nomab,length_word)
          endif
          l = len_trim(nomab)
          mot(i:i+l-1) = nomab(1:l)
          i = i+26
        end do
      end do
      i = i - 9
      mot(i+1:i+26) = '   Non-resonant Amplitude '
      i = i + 26
      mot(i+1:i+26) = ' e_s.e_i * Somme_exp(iQR) '
      i = i + 26

      l = len_trim(mot)
      write(7,405) Ephseuil*rydb, mot(1:l)

      do ipl = 1,npldafs

        if( Dafs_bio ) then
          write(7,407) nint( hkl_dafs(:,ipl) )
        else
          write(7,410) nint( hkl_dafs(:,ipl) ), isigpi(:,ipl)
        endif
        dang = 360._db / nphi_dafs(ipl)

        do ip = 1,nphi_dafs(ipl)
          ang = ( ip - 1 ) * dang
          if( ia == 0 ) then
            cf = (1._db,0._db)
          else
            cf = conjg( Bragg_abs(ia,ipl) )
          endif
          nw = 0
          do initlr = 1,ninitlr
            nw = nw + 1
            compnum(nw) = cf * Ampldafs(ipl,ip,initlr,ia)
            if( Full_self_abs ) then
              do ind_mu = 1,2
                nw = nw + 1
                compnum(nw) = mu(ipl,ip,ind_mu,initlr,ia)
              end do
            elseif( Self_abs ) then
              nw = nw + 1
              compnum(nw) = Cmplx( Real(mu(ipl,ip,1,initlr,ia), db), Real(mu(ipl,ip,2,initlr,ia),db),db)
            endif
          end do
          nw = nw + 1
          compnum(nw) = phdf0t(ipl,ip)
          nw = nw + 1
          compnum(nw) = phdt(ipl,ip,i_bulk_z)
          if( Dafs_bio ) then
            write(7,420) compnum(1:nw)
          else
            write(7,430) ang, compnum(1:nw)
          endif
        end do
      end do
      close(7)
    endif

      if( i_bulk_z >= n_bulk_z_abs .or. Bulk_roughness_case ) exit
    end do
    
  end do  ! end of loop over atoms

  if( dafs ) then
    deallocate( Ampldafs )
    if( E1E1 ) deallocate( Ampldafs_dd )
    if( E1E2 ) deallocate( Ampldafs_dq )
    if( E2E2 ) deallocate( Ampldafs_qq )
    if( E1E3 ) deallocate( Ampldafs_do )
    if( E3E3 ) deallocate( Ampldafs_oo )
    if( M1M1 ) deallocate( Ampldafs_mm )
    if( E1M1 ) deallocate( Ampldafs_md )
    if( Cor_abs ) then
      deallocate( mu )
      if( E1E1 ) deallocate( mu_dd )
      if( E1E2 ) deallocate( mu_dq )
      if( E2E2 ) deallocate( mu_qq )
      if( E1E3 ) deallocate( mu_do )
      if( E3E3 ) deallocate( mu_oo )
      if( M1M1 ) deallocate( mu_mm )
      if( E1M1 ) deallocate( mu_md )
    endif
  endif

  return
  110 format(/' ---- Coab_rgh -----',100('-'))
  115 format(/' ---- Coabs --------',100('-'))
  120 format(/' Conversion factor (numb. of electron/Mbarn) =',10f10.5)
  130 format(' Conversion factor (microm^(-1)/Mbarn)       =',10f10.5)
  140 format(/' Output in mm^(-1), conversion factor =',f13.5,' mm^(-1)/Mbarn' )
  150 format(/'   Total signal')
  290 format(/'   Signal atom',i3)
  295 format(/'   Core state or edge',i3)
  297 format(/'   Contribution for valence state l =',i3)
  300 format(/4x,'Energy',3200a13)
  310 format(f10.3,1p,3200e13.5)
  320 format(6x,'E1E1',1p,3200e13.5)
  330 format(6x,'E1E2',1p,3200e13.5)
  340 format(6x,'E2E2',1p,3200e13.5)
  350 format(6x,'E1E3',1p,3200e13.5)
  351 format(6x,'E3E3',1p,3200e13.5)
  352 format(6x,'M1M1',1p,3200e13.5)
  354 format(6x,'E1M1',1p,3200e13.5)
  360 format(/10x,1p,3200a13)
  370 format('  Ampldafs',1p,3200e13.5)
  392 format(' Number of Energies =',i5)
  400 format(i5,4x,' = Number of angles')
  405 format(f10.3,A)
  407 format(3i5,' = (h,k,l)')
  410 format(' (h,k,l) = ',3i3,', sigpi =',2i3)
  420 format(7x,1p,3200e13.5)
  430 format(f7.1,1p,3200e13.5)
end

!***********************************************************************

! Writing of the cartesian tensors in the bav file

subroutine Write_ten_bav(Core_resolved, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Final_tddft, Green_int, Green_int_mag, icheck, &
                         idafs, ipl, ipldafs, jseuil, M1M1, Magn_sens, n_oo, n_rel, natomsym, nb, nbseuil, ninit1, ninitlr, &
                         nseuil, secddia, secddia_m, secdoia, secdoia_m, secdqia, secdqia_m, secmdia, secmdia_m, secmmia, &
                         secmmia_m, secooia, secooia_m, secqqia, secqqia_m, Tens_comp) 

  use declarations
  implicit none

  integer:: he, hs, ia, ib, icheck, initlr, ipl, ipldafs, iseuil, isp1, isp2, isp3, isp4, ispfg, j1, je, jhe, jhs, js, jseuil, &
            ke, ks, n_oo, n_rel, natomsym, nb, nbseuil, ninit1, ninitlr, nseuil

  complex(kind=db), dimension(3,3,n_rel,ninitlr,0:natomsym):: secddia, secddia_m
  complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secmdia, secmdia_m, secmmia, secmmia_m
  complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia, secdqia_m
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secdoia, secdoia_m, secqqia, secqqia_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:natomsym):: secooia, secooia_m

  logical:: Core_resolved, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Final_tddft, Green_int, Green_int_mag, idafs, M1M1, Magn_sens, &
            Tens_comp
  
  do initlr = 1,ninitlr

    if( Core_resolved .and. .not. Final_tddft ) then
      if( initlr <= ninit1 ) then
        iseuil = 1
      else
        iseuil = min(2, nbseuil)
      endif
    elseif( Final_tddft ) then
      iseuil = min(2, nbseuil)
    else
      iseuil = initlr
    endif

    if( nseuil == 0 ) then  ! optic
      if( Core_resolved ) then
        if( initlr == 1 ) then
          write(3,123) ninitlr
        else
          write(3,125) initlr-2, ninitlr
        endif
      else
        write(3,120) 'Opt'
      endif
    elseif( Final_tddft ) then
      if( nbseuil == 2 ) then
        write(3,120) achar(nseuil+74) // achar(jseuil+iseuil+46) // achar(jseuil+iseuil+47)
      else
        write(3,120) achar(nseuil+74) // achar(jseuil+iseuil+46)
      endif
    else
      if( Core_resolved ) then
        write(3,130) initlr, ninitlr
      elseif( ninitlr > 1 ) then
        write(3,135) initlr, ninitlr
      endif
    endif

    do ib = 0,nb
      if( natomsym == 1 .and. ib > 0 ) cycle
      if( ib == 0 ) then
        ia = 1
      elseif( ib == 1 ) then
        ia = 0
      else
        ia = ib
      endif
      if( ia /= 0 .and. ipl > 1 ) cycle
      if( ipl > 1 .and. .not. idafs ) cycle
      if( idafs .and. icheck < 2 ) cycle
      if( ipl > 1 ) write(3,140) ipldafs
      if( E1E1 ) then
        if( ia == 1 ) then
          if( Green_int_mag ) then
            write(3,141)
          elseif( Green_int ) then
            write(3,142)
          else
            write(3,143)
          endif
        elseif( ia == 0 ) then
          if( Green_int_mag ) then
            write(3,144)
          elseif( Green_int ) then
            write(3,145)
          else
            write(3,146)
          endif
        else
          if( Green_int_mag ) then
            write(3,147) ia
          elseif( Green_int ) then
            write(3,148) ia
          else
            write(3,149) ia
          endif
        endif
        ispfg = 0
        loop_dip_spinpos: do isp1 = 1,2
        do isp2 = 1,2
        do isp3 = 1,2
        do isp4 = 1,2
          if( n_rel == 8 .and. isp1 /= isp4 ) cycle
          ispfg = ispfg + 1
          if( ispfg > n_rel ) exit loop_dip_spinpos
          if( ispfg > 1 ) write(3,'(a26,4i2)') 'isp1, isp2, isp3, isp4 =', isp1, isp2, isp3, isp4 
          do ke = 1,3
            if( Green_int_mag ) then
              write(3,150) secddia(ke,:,ispfg,initlr,ia), secddia_m(ke,:,ispfg,initlr,ia)
            elseif( Tens_comp ) then
              write(3,150) secddia(ke,:,ispfg,initlr,ia)
            else
              write(3,150) real( secddia(ke,:,ispfg,initlr,ia) )
            endif
          end do
        end do
        end do
        end do
        end do loop_dip_spinpos
      endif
      if( E1E2 ) then
        do ke = 1,3
          if( ia == 1 ) then
            if( Green_int_mag ) then
              write(3,160) ke
            elseif( Green_int ) then
              write(3,161) ke
            else
              write(3,162) ke
            endif
          elseif( ia == 0 ) then
            if( Green_int_mag ) then
              write(3,163) ke
            elseif( Green_int ) then
              write(3,164) ke
            elseif( Magn_sens ) then
              write(3,165) ke
            else
              write(3,166) ke
            endif
          else
            if( Green_int_mag ) then
              write(3,167) ia, ke
            elseif( Green_int ) then
              write(3,168) ia, ke
            else
              write(3,169) ia, ke
            endif
          endif
          do ks = 1,3
            if( ( Magn_sens .and. ia == 0 ) .or. Green_int_mag ) then
              write(3,150) secdqia(ke,ks,:,initlr,ia), secdqia_m(ke,ks,:,initlr,ia)
            elseif( Tens_comp ) then
              write(3,150) secdqia(ke,ks,:,initlr,ia)
            else
              write(3,150) real( secdqia(ke,ks,:,initlr,ia), db )
            endif
          end do
        end do
      endif
      if( E2E2 ) then
        do js = 1,3
          do ks = 1,3
            if( ia == 1 ) then
              if( Green_int_mag ) then
                write(3,171) ks, js
              elseif( Green_int ) then
                write(3,172) ks, js
              else
                write(3,173) ks, js
              endif
            elseif( ia == 0 ) then
              if( Green_int_mag ) then
                write(3,174) ks, js
              elseif( Green_int ) then
                write(3,175) ks, js
              else
                write(3,176) ks, js
              endif
            else
              if( Green_int_mag ) then
                write(3,177) ia, ks, js
              elseif( Green_int ) then
                write(3,178) ia, ks, js
              else
                write(3,179) ia, ks, js
              endif
            endif
            do ke = 1,3
              if( Green_int_mag ) then
                write(3,150) secqqia(ke,1:3,ks,js,initlr,ia), secqqia_m(ke,1:3,ks,js,initlr,ia)
              elseif( Tens_comp ) then
                write(3,150) secqqia(ke,1:3,ks,js,initlr,ia)
              else
                write(3,150) real( secqqia(ke,1:3,ks,js,initlr,ia), db )
              endif
            end do
          end do
        end do
      endif
      if( E1E3 ) then
        do ke = 1,3
          do ks = 1,3
            if( ia == 1 ) then
              if( Green_int_mag ) then
                write(3,181) ke, ks
              elseif( Green_int ) then
                write(3,182) ke, ks
              else
                write(3,183) ke, ks
              endif
            elseif( ia == 0 ) then
              if( Green_int_mag ) then
                write(3,184) ke, ks
              elseif( Green_int ) then
                write(3,185) ke, ks
              else
                write(3,186) ke, ks
              endif
            else
              if( Green_int_mag ) then
                write(3,187) ia, ks, js
              elseif( Green_int ) then
                write(3,188) ia, ke, ks
              else
                write(3,189) ia, ke, ks
              endif
            endif
            do j1 = 1,3
              if( Green_int_mag ) then
                write(3,150) secdoia(ke,ks,j1,:,initlr,ia), secdoia_m(ke,ks,j1,:,initlr,ia)
              elseif( Tens_comp ) then
                write(3,150) secdoia(ke,ks,j1,:,initlr,ia)
              else
                write(3,150) real( secdoia(ke,ks,j1,:,initlr,ia) )
              endif
            end do
          end do
        end do
      endif
      if( E1M1 ) then
        if( ia == 1 ) then
          if( Green_int_mag ) then
            write(3,190)
          elseif( Green_int ) then
            write(3,191)
          else
            write(3,192)
          endif
        elseif( ia == 0 ) then
          if( Green_int_mag ) then
            write(3,193)
          elseif( Green_int ) then
            write(3,194)
          elseif( Magn_sens ) then
            write(3,195)
          else
            write(3,196)
          endif
        else
          if( Green_int_mag ) then
            write(3,197) ia
          elseif( Green_int ) then
            write(3,198) ia
          else
            write(3,199) ia
          endif
        endif
        do ke = 1,3
          if( ( Magn_sens .and. ia == 0 ) .or. Green_int_mag ) then
            write(3,150) secmdia(ke,:,initlr,ia), secmdia_m(ke,:,initlr,ia)
          elseif( Tens_comp ) then
            write(3,150) secmdia(ke,:,initlr,ia)
          else
            write(3,150) aimag( secmdia(ke,:,initlr,ia) )
          endif
        end do
      endif
      if( M1M1 ) then
        if( ia == 1 ) then
          if( Green_int_mag ) then
            write(3,201)
          elseif( Green_int ) then
            write(3,202)
          else
            write(3,203)
          endif
        elseif( ia == 0 ) then
          if( Green_int_mag ) then
            write(3,204)
          elseif( Green_int ) then
            write(3,205)
          else
            write(3,206)
          endif
        else
          if( Green_int_mag ) then
            write(3,207) ia
          elseif( Green_int ) then
            write(3,208) ia
          else
            write(3,209) ia
          endif
        endif
        do ke = 1,3
          if( Green_int_mag ) then
            write(3,150) secmmia(ke,:,initlr,ia), secmmia_m(ke,:,initlr,ia)
          elseif( Tens_comp ) then
            write(3,150) secmmia(ke,:,initlr,ia)
          else
            write(3,150) real( secmmia(ke,:,initlr,ia) )
          endif
        end do
      endif

      if( E3E3 ) then
        do hs = 1,3
          do js = 1,3
            jhs = 3 * ( js - 1 ) + hs
            if( ia == 1 ) then
              if( Green_int_mag ) then
                write(3,211) (ks, js, hs, ks = 1,3)
              elseif( Green_int ) then
                write(3,212) (ks, js, hs, ks = 1,3)
              else
                write(3,213) (ks, js, hs, ks = 1,3)
              endif
            elseif( ia == 0 ) then
              if( Green_int_mag ) then
                write(3,214) (ks, js, hs, ks = 1,3)
              elseif( Green_int ) then
                write(3,215) (ks, js, hs, ks = 1,3)
              else
                write(3,216) (ks, js, hs, ks = 1,3)
              endif
            else
              if( Green_int_mag ) then
                write(3,217) ia, (ks, js, hs, ks = 1,3)
              elseif( Green_int ) then
                write(3,218) ia, (ks, js, hs, ks = 1,3)
              else
                write(3,219) ia, (ks, js, hs, ks = 1,3)
              endif
            endif
            do he = 1,3
              do je = 1,3
                jhe = 3 * ( je - 1 ) + he
                if( Green_int_mag ) then
                  write(3,150) ( secooia(1:3,jhe,ks,jhs,initlr,ia), secooia_m(1:3,jhe,ks,jhs,initlr,ia), ks = 1,3 )
                elseif( tens_comp ) then
                  write(3,150) ( secooia(1:3,jhe,ks,jhs,initlr,ia), ks = 1,3 )
                else
                  write(3,150) ( real( secooia(1:3,jhe,ks,jhs,initlr,ia) ), ks = 1,3 )
                endif
              end do
            end do
          end do
        end do
      endif

    end do

  end do

  return
  120 format(/' Threshold ',a3)
  123 format(/'   Valence state, all on',i3)
  125 format(/'   Valence state l =',i3,' on',i3)
  130 format(/'   Initial state =',i3,' on',i3)
  135 format(/'            Edge =',i3,' on',i3)
  140 format(/8x,' RXS polarization number',i3)
  141 format(/' Tensor_dd(ke,ks), Prototypical atom,     Green integral, Non magnetic part',82x,'Magnetic part')
  142 format(/' Tensor_dd(ke,ks), Prototypical atom,     Green integral')
  143 format(/' Tensor_dd(ke,ks), Prototypical atom')
  144 format(/' Crystal Tensor_dd(ke,ks), Green integral, Non magnetic part',95x,'Magnetic part')
  145 format(/' Crystal Tensor_dd(ke,ks), Green integral')
  146 format(/' Crystal Tensor_dd(ke,ks)')
  147 format(/' Atom ',i3,' Tensor_dd(ke,ks), Green integral, Non magnetic part',80x,'Magnetic part')
  148 format(/' Atom ',i3,' Tensor_dd(ke,ks), Green integral')
  149 format(/' Atom ',i3,' Tensor_dd(ke,ks)')
  150 format(1p,12(6e18.10,2x))
  160 format(/' Tensor_dq(',i1,',ks,j2), Prototypical atom,   Green integral, Non magnetic part',80x,'Magnetic part')
  161 format(/' Tensor_dq(',i1,',ks,j2), Prototypical atom,   Green integral')
  162 format(/' Tensor_dq(',i1,',ks,j2), Prototypical atom')
  163 format(/' Crystal Tensor_dq(',i1,',ks,j2), Green integral, Non magnetic part',89x,'Magnetic part')
  164 format(/' Crystal Tensor_dq(',i1,',ks,j2), Green integral')
  165 format(/' Crystal Tensor_dq(',i1,',ks,j2),                 Non magnetic part',89x,'Magnetic part')
  166 format(/' Crystal Tensor_dq(',i1,',ks,j2)')
  167 format(/' Atom ',i3,' Tensor_dq(ke,ks), Green integral, Non magnetic part',80x,'Magnetic part')
  168 format(/' Atom ',i3,' Tensor_dq(',i1,',ks,j2), Green integral')
  169 format(/' Atom ',i3,' Tensor_dq(',i1,',ks,j2)')
  171 format(/' Tensor_qq(ke,je,',i1,',',i1,'), Prototypical atom, Green integral, Non magnetic part',78x,'Magnetic part')
  172 format(/' Tensor_qq(ke,je,',i1,',',i1,'), Prototypical atom, Green integral')
  173 format(/' Tensor_qq(ke,je,',i1,',',i1,'), Prototypical atom')
  174 format(/' Crystal Tensor_qq(ke,je,',i1,',',i1,'), Green integral, Non magnetic part',87x,'Magnetic part')
  175 format(/' Crystal Tensor_qq(ke,je,',i1,',',i1,'), Green integral')
  176 format(/' Crystal Tensor_qq(ke,je,',i1,',',i1,')')
  177 format(/' Atom ',i3,' Tensor_qq(ke,je,',i1,',',i1,'), Green integral, Non magnetic part',80x,'Magnetic part')
  178 format(/' Atom ',i3,' Tensor_qq(ke,je,',i1,',',i1,'), Green integral')
  179 format(/' Atom ',i3,' Tensor_qq(ke,je,',i1,',',i1,')')
  181 format(/' Tensor_do(',i1,',',i1,',j1,j2), Prototypical atom, Green integral, Non magnetic part',78x,'Magnetic part')
  182 format(/' Tensor_do(',i1,',',i1,',j1,j2), Prototypical atom, Green integral')
  183 format(/' Tensor_do(',i1,',',i1,',j1,j2), Prototypical atom')
  184 format(/' Crystal Tensor_do(',i1,',',i1,',j1,j2), Green integral, Non magnetic part',87x,'Magnetic part')
  185 format(/' Crystal Tensor_do(',i1,',',i1,',j1,j2), Green integral')
  186 format(/' Crystal Tensor_do(',i1,',',i1,',j1,j2)')
  187 format(/' Atom ',i3,' Tensor_do(',i1,',',i1,',j1,j2), Green integral, Non magnetic part',80x,'Magnetic part')
  188 format(/' Atom ',i3,' Tensor_do(',i1,',',i1,',j1,j2), Green integral')
  189 format(/' Atom ',i3,' Tensor_do(',i1,',',i1,',j1,j2)')
  190 format(/' Tensor_md(ke,ks), Prototypical atom,     Green integral, Non magnetic part',82x,'Magnetic part')
  191 format(/' Tensor_md(ke,ks), Prototypical atom,     Green integral')
  192 format(/' Tensor_md(ke,ks), Prototypical atom')
  193 format(/' Crystal Tensor_md(ke,ks), Green integral, Non magnetic part',89x,'Magnetic part')
  194 format(/' Crystal Tensor_md(ke,ks), Green integral')
  195 format(/' Crystal Tensor_md(ke,ks),                 Non magnetic part',89x,'Magnetic part')
  196 format(/' Crystal Tensor_md(ke,ks)')
  197 format(/' Atom ',i3,' Tensor_md(ke,ks), Green integral, Non magnetic part',80x,'Magnetic part')
  198 format(/' Atom ',i3,' Tensor_md(ke,ks), Green integral')
  199 format(/' Atom ',i3,' Tensor_md(ke,ks)')
  201 format(/' Tensor_mm(ke,ks), Prototypical atom, Green integral, Non magnetic part',82x,'Magnetic part')
  202 format(/' Tensor_mm(ke,ks), Prototypical atom, Green integral')
  203 format(/' Tensor_mm(ke,ks), Prototypical atom')
  204 format(/' Crystal Tensor_mm(ke,ks), Green integral, Non magnetic part',91x,'Magnetic part')
  205 format(/' Crystal Tensor_mm(ke,ks), Green integral')
  206 format(/' Crystal Tensor_mm(ke,ks)')
  207 format(/' Atom ',i3,' Tensor_mm(ke,ks), Green integral, Non magnetic part',80x,'Magnetic part')
  208 format(/' Atom ',i3,' Tensor_mm(ke,ks), Green integral')
  209 format(/' Atom ',i3,' Tensor_mm(ke,ks)')
  211 format(/3(' Tensor_oo(ke,je,he',3(',',i1),')',36x), ' Prototypical atom,', &
          ' Green integral, Non magnetic part',78x,'Magnetic part')
  212 format(/3(' Tensor_oo(ke,je,he',3(',',i1),')',36x), ' Prototypical atom, Green integral')
  213 format(/3(' Tensor_oo(ke,je,he',3(',',i1),')',36x), ' Prototypical atom')
  214 format(/3(' Crystal Tensor_oo(ke,je,he',3(',',i1),')',36x), ' Green integral, Non magnetic part',87x,'Magnetic part')
  215 format(/3(' Crystal Tensor_oo(ke,je,he',3(',',i1),')',36x), ' Green integral')
  216 format(/3(' Crystal Tensor_oo(ke,je,he',3(',',i1),')',36x))
  217 format(/' Atom ',i3,3(' Tensor_oo(ke,je,he',3(',',i1),')',36x), &
          ' Green integral, Non magnetic part',80x,'Magnetic part')
  218 format(/' Atom ',i3,3(' Tensor_oo(ke,je,he',3(',',i1),')',36x), ' Green integral')
  219 format(/' Atom ',i3,3(' Tensor_oo(ke,je,he',3(',',i1),')',36x))
end
      
!******************************************************************************************************************************

subroutine Write_nrixs(Abs_U_iso,All_nrixs,Allsite,axyz,Core_resolved,E_cut,Energ,Energphot,Extract_ten, &
                  Epsii,Eseuil,f_avantseuil,Final_tddft,First_E,Green_int,i_range, &
                  iabsorig,icheck,ie,ie_computer,isymeq,jseuil,l0_nrixs,lmax_nrixs,mpinodee,n_multi_run, &
                  natomsym,nbseuil,nenerg,ninit1,ninitlr,nomfich,nomfich_cal_convt,nq_nrixs,nseuil,nspinp, &
                  numat_abs,Orthmatt,q_nrixs,Rot_atom_abs,Rot_int,S_nrixs,S_nrixs_m,Spinorbite,Taux_eq,V0muf,Volume_maille)
                  
  use declarations
  implicit none

  integer:: i, i_range, ia, iabsorig, icheck, ie, ie_computer, initlr, iq, isym, jseuil, l, l_i, l_s, l0_nrixs, &
    lm, lm_i, lm_s, lmax_nrixs, long, m_i, m_s, mpinodee, n, n_dim, n_multi_run, n_tens, natomsym, nbseuil, nenerg, ninit1, &
    ninitlr, nlmc, npldafs, nq_nrixs, nseuil, nspinp, numat_abs

  integer, dimension(natomsym):: isymeq

  character(len=length_word):: nomab
  character(len=Length_name) nomfich, nomfich_cal_convt, nomficht
  character(len=length_word), dimension(:), allocatable:: title

  complex(kind=db):: cfac, f_avantseuil, Ylm_q_i, Ylm_q_s
  complex(kind=db), dimension(0):: cdum
  complex(kind=db), dimension(nq_nrixs,(1+lmax_nrixs)**2,(1+lmax_nrixs)**2,ninitlr,0:mpinodee-1):: S_nrixs, S_nrixs_m
  complex(kind=db), dimension(:), allocatable:: Ylm_q
  
  logical:: All_nrixs, Allsite, Core_resolved, Final_tddft, First_E, Energphot, Extract_ten, Green_int, Green_int_mag, &
    Magn_sens, Monocrystal, Spinorbite

  real(kind=db):: Abs_U_iso, Surface_ref, E_cut, Ephoton, Ephseuil, natomsym_f, q, V0muf, Volume_maille

  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(3):: axyz, q_vec, q_vec_a
  real(kind=db), dimension(3,3):: Matopsym, Orthmatt, Rot_atom_abs, Rot_int
  real(kind=db), dimension(ninitlr) :: Epsii
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(natomsym):: Taux_eq
  real(kind=db), dimension(nq_nrixs,ninitlr,0:natomsym):: S_nrixs_T
  real(kind=db), dimension(nq_nrixs,l0_nrixs:lmax_nrixs,ninitlr,0:natomsym):: S_nrixs_T_l
  real(kind=db), dimension(4,nq_nrixs):: q_nrixs
  real(kind=db), dimension(:), allocatable:: Tens

  if( icheck > 1 ) write(3,110)
  
  npldafs = 0
  Surface_ref = 0._db
  natomsym_f = sum( Taux_eq(:) )

  Ephseuil = Energ(ie)
  Ephoton = Ephseuil + Eseuil(nbseuil)
! Pour les seuils de tres basse Energie:
  Ephoton = max(0.001_db/Rydb, Ephoton)
  if( Energphot ) Ephseuil = Ephoton

  if( ( jseuil > 1 .and. nspinp == 2 ) .or. Spinorbite ) then
    Magn_sens = .true.
  else
    Magn_sens = .false.
  endif

  if( Green_int .and. Magn_sens ) then
    Green_int_mag = .true.
  else
    Green_int_mag = .false.
  endif

  S_nrixs_T(:,:,:) = 0._db
  S_nrixs_T_l(:,:,:,:) = 0._db

  do iq = 1,nq_nrixs

    q_vec(:) = q_nrixs(1:3,iq)
    Monocrystal = sum( abs(q_vec(:)) ) > eps10

    if( icheck > 1 ) write(3,120) iq, q_nrixs(4,iq) / bohr

    if( Monocrystal ) then    
 
      if( icheck > 1 ) write(3,130) iq, q_nrixs(4,iq), q_vec(:) 
    
! Conversion to internal orthonormalized basis
      q_vec(:) = q_vec(:) / axyz(:)
      q_vec = matmul( Orthmatt, q_vec )
      q_vec(:) = q_vec(:) / sqrt( sum( q_vec(:)**2 ) )

      if( icheck > 1 ) write(3,140) q_vec(:) 

      nlmc = ( ( lmax_nrixs + 1 ) * ( lmax_nrixs + 2 ) ) / 2
      allocate( Ylm_q(nlmc) )
    
    else

      if( icheck > 1 ) write(3,140) iq, q_nrixs(4,iq)
        
    endif

    do ia = 1,natomsym

      if( Monocrystal ) then
        isym = abs( isymeq(ia) )
        call opsym(isym,Matopsym)
        q_vec_a = Matmul( Transpose(Matopsym), q_vec )
! Rotation to the absorbing atom basis
        q_vec_a = Matmul( Rot_atom_abs, Matmul( Transpose( Rot_int ), q_vec_a ) )

        call cYlm(lmax_nrixs,q_vec_a,1._db,Ylm_q,nlmc)
      endif

      lm_s = 0   
      do l_s = l0_nrixs,lmax_nrixs
        do m_s = -l_s,l_s
          lm_s = lm_s + 1

          if( Monocrystal ) then    
            lm = l_s * ( l_s + 1 ) / 2 + 1 + abs(m_s)
            if( m_s >= 0 ) then
              Ylm_q_s = Ylm_q(lm)
            else
              Ylm_q_s = (-1)**m_s * conjg( Ylm_q(lm) )
            endif
          endif
        
          lm_i = 0
          do l_i = l0_nrixs,lmax_nrixs
            do m_i = -l_i,l_i
              lm_i = lm_i + 1
    
              if( ( l_s /= l_i .or. m_s /= m_i ) .and. .not. Monocrystal ) cycle

              if( Monocrystal ) then    
                lm = l_i * ( l_i + 1 ) / 2 + 1 + abs(m_i)
                if( m_i >= 0 ) then
                  Ylm_q_i = Ylm_q(lm)
                else
                  Ylm_q_i = (-1)**m_i * conjg( Ylm_q(lm) )
                endif
              endif

              if( Monocrystal ) then
                cfac = Taux_eq(ia) * conjg( Ylm_q_s ) * Ylm_q_i 
              else
                cfac = cmplx( Taux_eq(ia) / quatre_pi, 0._db, db )  ! 4 pi is for the average  
              endif
              
              S_nrixs_T(iq,:,ia) = S_nrixs_T(iq,:,ia) + real( cfac * S_nrixs(iq,lm_s,lm_i,:,ie_computer), db ) 
              if( l_i == l_s ) &
                S_nrixs_T_l(iq,l_i,:,ia) = S_nrixs_T_l(iq,l_i,:,ia) + real( cfac * S_nrixs(iq,lm_s,lm_i,:,ie_computer), db )  
            
              if( Green_int_mag ) then
                if( isymeq(ia) < 0 ) then
                  S_nrixs_T(iq,:,ia) = S_nrixs_T(iq,:,ia) - real( cfac * S_nrixs_m(iq,lm_s,lm_i,:,ie_computer), db )
                  if( l_i == l_s ) &
                    S_nrixs_T_l(iq,l_i,:,ia) = S_nrixs_T_l(iq,l_i,:,ia) - real( cfac * S_nrixs_m(iq,lm_s,lm_i,:,ie_computer), db )
                else
                  S_nrixs_T(iq,:,ia) = S_nrixs_T(iq,:,ia) + real( cfac * S_nrixs_m(iq,lm_s,lm_i,:,ie_computer), db )
                  if( l_i == l_s ) &
                    S_nrixs_T_l(iq,l_i,:,ia) = S_nrixs_T_l(iq,l_i,:,ia) + real( cfac * S_nrixs_m(iq,lm_s,lm_i,:,ie_computer), db )
                endif
              endif
              
            end do
          end do
        end do
      end do
      
      S_nrixs_T(iq,:,0) =  S_nrixs_T(iq,:,0) + S_nrixs_T(iq,:,ia)
      S_nrixs_T_l(iq,:,:,0) =  S_nrixs_T_l(iq,:,:,0) + S_nrixs_T_l(iq,:,:,ia)
    
    end do

    if( Monocrystal ) deallocate( Ylm_q )
            
  end do

  nomficht = nomfich
  long = len_trim(nomficht)
  nomficht(long+1:long+6) = '_nrixs'

  if( All_nrixs ) then
    n_dim = nq_nrixs * ninitlr * ( lmax_nrixs - l0_nrixs + 2 )
  else
    n_dim = nq_nrixs * ninitlr
  endif
  n_tens = n_dim
  allocate( Title( n_dim ) )
  allocate( Tens( n_dim ) )
    
  do ia = 0,natomsym
    if( .not. Allsite .and. ia == 1 ) exit

    nomficht = nomfich
    long = len_trim(nomficht)
    nomficht(long+1:long+6) = '_nrixs'

    if( ia > 0 ) then
      long = len_trim(nomficht)
      nomficht(long+1:long+5) = '_atom'
      call ad_number(ia,nomficht,Length_name)
    endif
    long = len_trim(nomficht)

    if( Final_tddft .and. .not. Extract_ten ) nomficht(long+1:long+6) = '_tddft'

    if( n_multi_run > 1 ) then
      long = len_trim(nomficht)
      nomficht(long+1:long+1) = '_'
      call ad_number(iabsorig,nomficht,Length_name)
    endif

    long = len_trim(nomficht)
    nomficht(long+1:long+4) = '.txt'

    if( ie == 1 .and. ia == 0 ) nomfich_cal_convt = nomficht

    i = 0
    do initlr = 1,ninitlr

      do l = l0_nrixs - 1, lmax_nrixs

        if( ( .not. All_nrixs ) .and. l > l0_nrixs - 1 ) exit 

        do iq = 1, nq_nrixs 
          i = i + 1
          nomab = 'q='
          q = q_nrixs(4,iq) / bohr
          n = int( q + eps10 )   
          call ad_number(n,nomab,length_word)
          long = len_trim( nomab )
          long = long + 1
          nomab(long:long) = '.'
          if( q < 0.001_db ) then
            n = nint( 10000 * q )
            call ad_number(0,nomab,length_word)
            call ad_number(0,nomab,length_word)
            call ad_number(0,nomab,length_word)
          elseif( q < 0.01_db ) then
            n = nint( 1000 * q )
            call ad_number(0,nomab,length_word)
            call ad_number(0,nomab,length_word)
          elseif( q < 0.1_db ) then
            n = nint( 100 * q )
            call ad_number(0,nomab,length_word)
          elseif( n < 10 ) then
            n = int( 1000 * ( q - n ) )
          elseif( n < 100 ) then
            n = int( 100 * ( q - n ) )
          elseif( n < 1000 ) then
            n = int( 10 * ( q - n ) )
          endif 
          call ad_number(n,nomab,length_word)
          if( l > l0_nrixs - 1 ) then
            long = len_trim( nomab )
            long = long + 1
            nomab(long:long+1) = '_l'
            call ad_number(l,nomab,length_word)
          endif
          if( ninitlr > 1 .and. .not. (initlr == 1 .and. nseuil == 0 ) ) then
            long = len_trim( nomab )
            if( long > length_word - 3 ) cycle
            long = long + 1
            nomab(long:long) = '_'
            if( .not. Core_resolved ) then
              long = long + 1
              nomab(long:long) = achar(nseuil+74)
              long = long + 1
              nomab(long:long) = achar(jseuil+initlr+47)
            else
              if( nseuil == 0 ) then
                call ad_number(initlr-2,nomab,length_word)
              else
                call ad_number(initlr,nomab,length_word)
              endif
            endif
          endif
          call center_word(nomab,Length_word)
          title(i) = nomab

          if( l == l0_nrixs - 1 ) then
            Tens(i) = S_nrixs_T(iq,initlr,ia)
          else
            Tens(i) = S_nrixs_T_l(iq,l,initlr,ia)
          endif

        end do
      end do
    end do

! Writing in the output file
    call write_out(Abs_U_iso,rdum,rdum,f_avantseuil,E_cut,Ephseuil,Epsii,Eseuil(nbseuil),First_E,Green_int,rdum, &
            i_range,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,Title,npldafs,npldafs,0,0,nseuil,numat_abs, &
            cdum,cdum,Tens,V0muf,Core_resolved,natomsym_f,Surface_ref,Volume_maille)

! Writing on the screen
    if( ia == 0 ) call write_out(Abs_U_iso,rdum,rdum,f_avantseuil,E_cut,Ephseuil,Epsii,Eseuil(nbseuil),First_E,Green_int,rdum, &
            i_range,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,Title,npldafs,npldafs,0,0,nseuil,-1, &
            cdum,cdum,Tens,V0muf,Core_resolved,natomsym_f,Surface_ref,Volume_maille)

  end do
  
  deallocate( Tens, Title )
   
  return
  110 format(/' ---- Write_nrixs --------',100('-'))
  120 format(/'  iq =',i3,', q =',f6.3)
  130 format(/'  iq =',i3,', q =',f6.3,', q_vec =',3f10.5)
  140 format(/' Normalized q_vec in the absorbing atom local basis =',3f10.5)
end

!***********************************************************************

! Conversion factor Megabarn --> nbr. of electron 

function Conv_Mbarn_nelec(Ephoton)

  use declarations
  implicit none
  
  real(kind=db):: Conv_Mbarn_nelec, Ephoton 

! alfa_sf = e*e/(2*epsilon0*h*c) is the fine structure constant
! Factor 100 * bohr**2 is to get in megaBarn (10E-18 cm2)
  Conv_Mbarn_nelec = 0.5_db * Ephoton / ( quatre_pi * alfa_sf * 100 * bohr**2  )

  return
end

!***********************************************************************

subroutine Write_cartesian_tensor(Abs_U_iso,Volume_maille,E_cut,E1E2,E2E2,Ephseuil,Epsii,Eseuil,First_E,i_range, &
                 ia,ipldafs,jseuil,M1M1,magn_sens,n_rel,natomsym,ninit1,ninitlr,nomfich_s,nseuil,numat_abs,secddia,secdqia, &
                 secdqia_m,secqqia,secmdia,tens_comp,v0muf,Core_resolved)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( n_dim=10*(168+12) )

  character(len=Length_name) nomficht, nomfich_s
  character(len=Length_word) mot
  character(len=Length_word), dimension(n_dim):: nomtens

  integer:: i_range, n_rel
  
  complex(kind=db):: zero_c
  complex(kind=db), dimension(1):: cdum
  complex(kind=db), dimension(3,3,n_rel,ninitlr,0:natomsym):: secddia
  complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secmdia
  complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia, secdqia_m
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secqqia

  logical:: Core_resolved, E1E2, E2E2, First_E, M1M1, magn_sens, tens_comp

  real(kind=db):: Surface_ref
  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(ninitlr):: Epsii
  real(kind=db), dimension(n_dim):: Tens

  nomficht = nomfich_s
  long = len_trim(nomficht)
  if( ia == 0 ) then
    nomficht(long+1:long+9) = '_car_xtal'
  else
    nomficht(long+1:long+9) = '_car_atom'
    call ad_number(ia,nomficht,Length_name)
  endif
  if( ipldafs > 0 ) then
    long = len_trim(nomficht)
    nomficht(long+1:long+4) = '_rxs'
    call ad_number(ipldafs,nomficht,Length_name)
  endif
  long = len_trim(nomficht)
  nomficht(long+1:long+4) = '.txt'

  index = 0
  Surface_ref = 0._db

  do initlr = 1,ninitlr

    mot = ' '
    mot(1:2) = 'D_'
    do i = 1,3
      mot(3:3) = achar(i+119)
      do j = i,3
        index = index + 1
        mot(4:4) = achar(j+119)
        if( tens_comp ) mot(5:6) = '_r'
        Tens(index) = real( secddia(i,j,1,initlr,ia),db )
        nomtens(index) = mot
        if( tens_comp ) then
          index = index + 1
          mot(6:6) = 'i'
          nomtens(index) = mot
          Tens(index) = aimag( secddia(i,j,1,initlr,ia) )
        endif
      end do
    end do
    if( E1E2 ) then
      mot(1:2) = 'I_'
      do i = 1,3
        mot(3:3) = achar(i+119)
        do j = 1,3
          mot(4:4) = achar(j+119)
          do k = j,3
            mot(5:5) = achar(k+119)
            index = index + 1
            if( tens_comp ) mot(6:7) = '_r'
            nomtens(index) = mot
            Tens(index) = real( secdqia(i,j,k,initlr,ia),db )
            if( tens_comp ) then
              index = index + 1
              mot(7:7) = 'i'
              nomtens(index) = mot
              Tens(index) = aimag( secdqia(i,j,k,initlr,ia) )
            endif
            if( ia == 0 .and. magn_sens ) then
              index = index+1
              mot(7:11) = 'r_mag'
              nomtens(index) = mot
              Tens(index) = real( secdqia_m(i,j,k,initlr,ia),db )
              index = index+1
              mot(7:7) = 'i'
              nomtens(index) = mot
              mot(7:11) = '    '
              Tens(index) = aimag( secdqia_m(i,j,k,initlr,ia) )
            endif
          end do
        end do
      end do
    endif
    if( E2E2 ) then
      mot = ' '
      mot(1:2) = 'Q_'
      do i = 1,3
        mot(3:3) = achar(i+119)
        do j = i,3
          mot(4:4) = achar(j+119)
          do k = 1,3
            mot(5:5) = achar(k+119)
            do l = k,3
              mot(6:6) = achar(l+119)
              index = index + 1
              if( tens_comp ) mot(7:8) = '_r'
              nomtens(index) = mot
              Tens(index) = real( secqqia(i,j,k,l,initlr,ia),db )
              if( tens_comp ) then
                index = index + 1
                mot(8:8) = 'i'
                nomtens(index) = mot
                Tens(index) = aimag( secqqia(i,j,k,l,initlr,ia) )
              endif
            end do
          end do
        end do
      end do
      do ijk = 1,6
        index = index + 1
        select case(ijk)
          case(1)
            i = 1; j = 2; k = 1; l = 2
          case(2)
            i = 1; j = 2; k = 1; l = 3
          case(3)
            i = 1; j = 3; k = 1; l = 3
          case(4)
            i = 1; j = 3; k = 2; l = 2
          case(5)
            i = 1; j = 3; k = 2; l = 3
          case(6)
            i = 2; j = 3; k = 2; l = 3
        end select
        mot(3:3) = achar(i+119)
        mot(4:4) = achar(j+119)
        mot(5:5) = achar(k+119)
        mot(6:6) = achar(l+119)
        if( tens_comp ) mot(7:8) = '_r'
        nomtens(index) = mot
        Tens(index) = real( secqqia(i,j,k,l,initlr,ia),db )
        if( tens_comp ) then
          index = index + 1
          mot(8:8) = 'i'
          nomtens(index) = mot
          Tens(index) = aimag( secqqia(i,j,k,l,initlr,ia) )
        endif
      end do
    endif
    if( M1M1 ) then
      mot = ' '
      mot(1:2) = 'M_'
      do i = 1,3
        mot(3:3) = achar(i+119)
        do j = i,3
          index = index + 1
          mot(4:4) = achar(j+119)
          if( tens_comp ) mot(5:6) = '_r'
          Tens(index) = real( secmdia(i,j,initlr,ia),db )
          nomtens(index) = mot
          if( tens_comp ) then
            index = index + 1
            mot(6:6) = 'i'
            nomtens(index) = mot
            Tens(index) = aimag( secmdia(i,j,initlr,ia) )
          endif
        end do
      end do
    endif
    n_tens = index

  end do

  zero_c = (0._db, 0._db)

  call Write_out(Abs_U_iso,rdum,rdum,zero_c,E_cut,Ephseuil,Epsii,Eseuil,First_E,.false.,rdum, &
           i_range,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,nomtens,1,0,0,0,nseuil,numat_abs, &
           cdum,cdum,Tens,V0muf,Core_resolved,0._db,Surface_ref,Volume_maille) 

  return
end

!***********************************************************************

subroutine Write_out(Abs_U_iso,angxyz,axyz,f_avantseuil,E_cut,Ephseuil,Epsii,Eseuil,First_E,Green_int,hkl_dafs, &
            i_range,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,Title,np,npp,nppa,npps,nseuil,numat, &
            ph1,ph2,Tens,V0muf,Core_resolved,natomsym_f,Surface_ref,Volume_maille)

  use declarations
  implicit none

  integer, parameter:: n_tens_max = 10000

  integer:: i, i_bulk_z, i_range, icor, ipl, ipr, jseuil, n, n_bulk_z, n_dim, n_tens, ninit1, ninitlr, np, npp, nppa, npps, &
            nseuil, numat 
  
  character(len=Length_name):: nomficht
  character(len=92):: mot1
  character(len=63):: mot2
  character(len=Length_word):: mot
  character(len=Length_word), dimension(n_dim):: title
  character(len=10+(n_tens-2*npp)*Length_word):: dummy

  complex(kind=db):: f_avantseuil
  complex(kind=db), dimension(np):: ph1, ph2

  logical:: Core_resolved, First_E, Gnuplot, Green_int, Truncature

  real(kind=db):: Abs_U_iso, E_cut, Ephseuil, Eseuil, f0_forward, fpp_avantseuil, Surface, Surface_ref, natomsym_f, V0muf, &
                  Volume, Volume_maille
  
  real(kind=db), dimension(nppa):: angxyz, axyz
  real(kind=db), dimension(ninitlr):: Epsii
  real(kind=db), dimension(n_dim):: Tens
  real(kind=db), dimension(3,npps):: hkl_dafs

  Gnuplot = .false.
  dummy = ' '

  if( Length_word < 11 .or. Length_word > 17 ) then
    call write_error
    do ipr = 6,9,3
      write(ipr,100) Length_word
    end do
    stop
  endif
    
! For 2D resonant diffraction, when an absorbing atom is in the bulk,
! hkl_dafs(2,1) is a fake number which when high means that hkl_dafs(i,ipl), contains Length_abs.  
  if( npps > 0 ) then
    Truncature = hkl_dafs(2,1) > 1000._db
    if( Truncature ) then
      hkl_dafs(2,1) = hkl_dafs(2,1) - 10000._db
      n_bulk_z = npps / np - 1
    endif
  else
    Truncature = .false.
  endif
  
  if( numat == -1 ) then
    ipr = 6
  elseif( numat < -1 ) then
    ipr = - numat
  else
    ipr = 4
  endif

  if( ipr == 6 ) then
    n = min(n_tens,4)
  else
    n = n_tens
  endif

  if( n_tens > n_tens_max ) then
    call write_error
    do ipr = 6,9,3
      write(ipr,105) n_tens
    end do
    stop
  endif

  fpp_avantseuil = aimag( f_avantseuil )
  f0_forward = real( f_avantseuil, db )

  if( First_E ) then
    if( numat >= 0 ) then
      if( i_range == 1 ) then !***JDB Jan 2018, added ELSE option
        open(ipr, file = nomficht)
      else
        open(ipr, file = nomficht, position='append')
      endif
    elseif( numat < -1 ) then
      open(ipr, status='scratch')
    endif
    if( numat > 0 ) then
      if( Core_resolved .and. nseuil /= 0 ) then
        icor = ninit1
      else
        icor = 1
      endif
      if( Green_int ) icor = - icor
      
      mot1 = ' = E_edge, Z, n_edge, j_edge, Abs_before_edge, VO_interstitial, E_cut, ninitl, ninit1, Epsii'
      mot2 = ', UnitCell_Volume, Surface_ref, f0_forward, natomsym_f, abs_u_iso'

      Volume = Volume_maille * bohr**3
      Surface = Surface_ref * bohr**2
                       
      if( i_range == 1 ) then !***JDB Jan 2018
        if( nseuil == 0 ) then
          write(ipr,120) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, 1, icor, Epsii(1)*rydb, Volume, Surface, f0_forward, natomsym_f, Abs_U_iso, mot1, mot2
        elseif( ninitlr == 1 ) then
          write(ipr,120) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Volume, Surface, f0_forward, natomsym_f, Abs_U_iso, mot1, mot2
        elseif( ninitlr == 2 ) then
          write(ipr,121) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Volume, Surface, f0_forward, natomsym_f, Abs_U_iso, mot1, ninitlr, mot2
        elseif( ninitlr == 4 ) then
          write(ipr,122) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Volume, Surface, f0_forward, natomsym_f, Abs_U_iso, mot1, ninitlr, mot2
        elseif( ninitlr == 6 ) then
          write(ipr,123) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Volume, Surface, f0_forward, natomsym_f, Abs_U_iso, mot1, ninitlr, mot2
        elseif( ninitlr == 8 ) then
          write(ipr,124) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Volume, Surface, f0_forward, natomsym_f, Abs_U_iso, mot1, ninitlr, mot2
        elseif( ninitlr == 10 ) then
          write(ipr,125) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Volume, Surface, f0_forward, natomsym_f, Abs_U_iso, mot1, ninitlr, mot2
        else     ! 14
          write(ipr,126) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Volume, Surface, f0_forward, natomsym_f, Abs_U_iso, mot1, ninitlr, mot2
        endif
      endif

    endif

    if( npp > 0 .and. numat >= 0 ) then

      select case(Length_word)
        case(11)
          write(ipr,131) dummy, ph2(1:npp)
          write(ipr,131) dummy, ph1(1:npp)
        case(12)
          write(ipr,132) dummy, ph2(1:npp)
          write(ipr,132) dummy, ph1(1:npp)
        case(13)
          write(ipr,133) dummy, ph2(1:npp)
          write(ipr,133) dummy, ph1(1:npp)
        case(14)
          write(ipr,134) dummy, ph2(1:npp)
          write(ipr,134) dummy, ph1(1:npp)
        case(15)
          write(ipr,135) dummy, ph2(1:npp)
          write(ipr,135) dummy, ph1(1:npp)
        case(16)
          write(ipr,136) dummy, ph2(1:npp)
          write(ipr,136) dummy, ph1(1:npp)
        case(17)
          write(ipr,137) dummy, ph2(1:npp)
          write(ipr,137) dummy, ph1(1:npp)
      end select

      if( nppa > 0 ) write(ipr,145) axyz(:)*bohr, angxyz(:) / radian, ( nint( hkl_dafs(1:3,i) ), i = 1,npp )

      if( Truncature) then
        do i = 1,3
          select case(Length_word)
            case(15)
              write(ipr,151) dummy, hkl_dafs(i,1:np)
            case(16)
              write(ipr,152) dummy, hkl_dafs(i,1:np)
            case(17)
              write(ipr,153) dummy, hkl_dafs(i,1:np)
            case default
              write(ipr,151) dummy, hkl_dafs(i,1:np)
          end select
        end do
        do i_bulk_z = 1,n_bulk_z
          select case(Length_word)
            case(15)
              write(ipr,135) dummy, ( hkl_dafs(2,ipl), hkl_dafs(3,ipl), ipl = np*i_bulk_z+1,np*(i_bulk_z+1) ) 
              write(ipr,151) dummy, ( hkl_dafs(1,ipl), ipl = np*i_bulk_z+1,np*(i_bulk_z+1) )
            case(16)
              write(ipr,136) dummy, ( hkl_dafs(2,ipl), hkl_dafs(3,ipl), ipl = np*i_bulk_z+1,np*(i_bulk_z+1) ) 
              write(ipr,152) dummy, ( hkl_dafs(1,ipl), ipl = np*i_bulk_z+1,np*(i_bulk_z+1) )
            case(17)
              write(ipr,137) dummy, ( hkl_dafs(2,ipl), hkl_dafs(3,ipl), ipl = np*i_bulk_z+1,np*(i_bulk_z+1) ) 
              write(ipr,153) dummy, ( hkl_dafs(1,ipl), ipl = np*i_bulk_z+1,np*(i_bulk_z+1) )
            case default
              write(ipr,135) dummy, ( hkl_dafs(2,ipl), hkl_dafs(3,ipl), ipl = np*i_bulk_z+1,np*(i_bulk_z+1) ) 
              write(ipr,151) dummy, ( hkl_dafs(1,ipl), ipl = np*i_bulk_z+1,np*(i_bulk_z+1) )
            end select
        end do
      endif
    endif
    do i = 1,n
      mot = title(i)
      call center_word( mot, Length_word )
      title(i) = mot
    end do
    if( i_range == 1 ) then !***JDB
      if( Gnuplot ) then
        write(ipr,158) title(1:n)
      else
        write(ipr,160) title(1:n)
      endif
    endif
  elseif( numat >= 0 ) then
    open(ipr, file = nomficht, position='append')
  endif

  if( abs( Ephseuil*rydb ) < 9.999995_db ) then
    select case(Length_word)
      case(11)
        write(ipr,161) Ephseuil*rydb, Tens(1:n)  
      case(12)
        write(ipr,162) Ephseuil*rydb, Tens(1:n)
      case(13)
        write(ipr,163) Ephseuil*rydb, Tens(1:n)
      case(14)
        write(ipr,164) Ephseuil*rydb, Tens(1:n)
      case(15)
        write(ipr,165) Ephseuil*rydb, Tens(1:n)
      case(16)
        write(ipr,166) Ephseuil*rydb, Tens(1:n)
      case(17)
        write(ipr,167) Ephseuil*rydb, Tens(1:n)
    end select
  elseif( abs( Ephseuil*rydb ) < 99.99995_db ) then
    select case(Length_word)
      case(11)
        write(ipr,171) Ephseuil*rydb, Tens(1:n)  
      case(12)
        write(ipr,172) Ephseuil*rydb, Tens(1:n)
      case(13)
        write(ipr,173) Ephseuil*rydb, Tens(1:n)
      case(14)
        write(ipr,174) Ephseuil*rydb, Tens(1:n)
      case(15)
        write(ipr,175) Ephseuil*rydb, Tens(1:n)
      case(16)
        write(ipr,176) Ephseuil*rydb, Tens(1:n)
      case(17)
        write(ipr,177) Ephseuil*rydb, Tens(1:n)
    end select
  elseif( abs( Ephseuil*rydb ) < 999.9995_db ) then
    select case(Length_word)
      case(11)
        write(ipr,181) Ephseuil*rydb, Tens(1:n)  
      case(12)
        write(ipr,182) Ephseuil*rydb, Tens(1:n)
      case(13)
        write(ipr,183) Ephseuil*rydb, Tens(1:n)
      case(14)
        write(ipr,184) Ephseuil*rydb, Tens(1:n)
      case(15)
        write(ipr,185) Ephseuil*rydb, Tens(1:n)
      case(16)
        write(ipr,186) Ephseuil*rydb, Tens(1:n)
      case(17)
        write(ipr,187) Ephseuil*rydb, Tens(1:n)
    end select
  else
    select case(Length_word)
      case(11)
        write(ipr,191) Ephseuil*rydb, Tens(1:n)  
      case(12)
        write(ipr,192) Ephseuil*rydb, Tens(1:n)
      case(13)
        write(ipr,193) Ephseuil*rydb, Tens(1:n)
      case(14)
        write(ipr,194) Ephseuil*rydb, Tens(1:n)
      case(15)
        write(ipr,195) Ephseuil*rydb, Tens(1:n)
      case(16)
        write(ipr,196) Ephseuil*rydb, Tens(1:n)
      case(17)
        write(ipr,197) Ephseuil*rydb, Tens(1:n)
    end select
  endif

  if( numat >= 0 ) close(ipr)

  return
  100 format(//' Length_word =',i3, ' This parameter must be set between 11 and 17 !'//)
  105 format(//' The number of column to be written is ',i5, // &
               ' This is greater than the maximum possible value given in the routine write_out in the file coabs.f !', / & 
               ' To change that you must modify the formats 121 up to 157 in this routine,', / &
               ' and increase to the same value the parameter n_tens_max.', / &
               ' Then you compile again.' //)

  120 format(f10.3,i5,2i3,1p,3e15.7,2i3,e15.7,5e15.7,a92,a63)
  121 format(f10.3,i5,2i3,1p,3e15.7,2i3,2e15.7,5e15.7,a92,'(1..',i1,')',a63)
  122 format(f10.3,i5,2i3,1p,3e15.7,2i3,4e15.7,5e15.7,a92,'(1..',i1,')',a63)
  123 format(f10.3,i5,2i3,1p,3e15.7,2i3,6e15.7,5e15.7,a92,'(1..',i1,')',a63)
  124 format(f10.3,i5,2i3,1p,3e15.7,2i3,8e15.7,5e15.7,a92,'(1..',i1,')',a63)
  125 format(f10.3,i5,2i3,1p,3e15.7,2i3,10e15.7,5e15.7,a92,'(1..',i2,')',a63)
  126 format(f10.3,i5,2i3,1p,3e15.7,2i3,14e15.7,5e15.7,a92,'(1..',i2,')',a63)

  131 format(A,1p,10000e11.3)
  132 format(A,1p,10000e12.4)
  133 format(A,1p,10000e13.5)
  134 format(A,1p,10000e14.6)
  135 format(A,1p,10000e15.7)
  136 format(A,1p,10000e16.8)
  137 format(A,1p,10000e17.9)
  145 format(3f10.5,3x,3f10.5,10000(14x,3i4))
  151 format(A,1p,10000(15x,e15.7))
  152 format(A,1p,10000(16x,e16.8))
  153 format(A,1p,10000(17x,e17.9))
  158 format('#   Energy',10000A)
  160 format('    Energy',10000A)
  161 format(f10.5,1p,10000e11.3)
  162 format(f10.5,1p,10000e12.4)
  163 format(f10.5,1p,10000e13.5)
  164 format(f10.5,1p,10000e14.6)
  165 format(f10.5,1p,10000e15.7)
  166 format(f10.5,1p,10000e16.8)
  167 format(f10.5,1p,10000e17.9)
  171 format(f10.4,1p,10000e11.3)
  172 format(f10.4,1p,10000e12.4)
  173 format(f10.4,1p,10000e13.5)
  174 format(f10.4,1p,10000e14.6)
  175 format(f10.4,1p,10000e15.7)
  176 format(f10.4,1p,10000e16.8)
  177 format(f10.4,1p,10000e17.9)
  181 format(f10.3,1p,10000e11.3)
  182 format(f10.3,1p,10000e12.4)
  183 format(f10.3,1p,10000e13.5)
  184 format(f10.3,1p,10000e14.6)
  185 format(f10.3,1p,10000e15.7)
  186 format(f10.3,1p,10000e16.8)
  187 format(f10.3,1p,10000e17.9)
  191 format(f10.3,1p,10000e11.3)
  192 format(f10.2,1p,10000e12.4)
  193 format(f10.2,1p,10000e13.5)
  194 format(f10.2,1p,10000e14.6)
  195 format(f10.2,1p,10000e15.7)
  196 format(f10.2,1p,10000e16.8)
  197 format(f10.2,1p,10000e17.9)
end

!***********************************************************************

subroutine center_word( mot, Length_word )

  character(len=*) mot
  character(len=Length_word) mot2

  mot2 = ' '
  mot = adjustl( mot )
  l = len_trim( mot )
  lshift = ( Length_word - l + 1 ) / 2
  lshift = max( 0, lshift )
  lm = min( l, Length_word )
  mot2(1+lshift:lm+lshift) = mot(1:lm)
  mot = mot2

  return
end

!***********************************************************************

subroutine ad_number(ib,nomfich,Length)

  character(len=Length) nomfich

  l = len_trim(nomfich)

  if( ib < 0 ) then
    l = l + 1
    if( l <= Length ) nomfich(l:l) = '-'
  endif

  i = abs(ib)

  do iu = 2,10
    if( i / 10**(iu-1) < 1 ) exit
  end do
  iumax = iu - 1

  do iu = iumax,1,-1

    ipuis = 10**(iu-1)

    in = i / ipuis

! S'il n'y a pas la place on ecrit que les derniers chiffres
    if( l + iu <= Length ) then
      l = l + 1
      nomfich(l:l) = achar(in+48)
    elseif( l + iu == Length-1 ) then
      if( nomfich(l:l) == '_' ) nomfich(l:l) = achar(in+48)
    endif

    i = i - ipuis * in

  end do

  return
end

!***********************************************************************

subroutine ad_number_r(r,nom)

  use declarations
  implicit none

  integer:: i, l, l2, l3

  character(len=1) a
  character(len=Length_word) mot, nom

  real(kind=db):: r
  
  if( Length_word >= 15 ) then 
    write(mot,'(f15.4)') r
  elseif( Length_word >= 13 ) then 
    write(mot,'(f13.3)') r
  elseif( Length_word >= 10 ) then 
    write(mot,'(f10.2)') r
  else
    write(mot,'(f5.2)') r
  endif
  
  mot = adjustl( mot )
  
  l2 = len_trim(mot) 
  l3 = l2

  do i = l3,2,-1
    a = mot(i:i)
    if( a /= '0' .and. a /= '.' ) exit
    mot(i:i) = ' '
    l2 = l2 - 1
    if( a == '.' ) exit
  end do
 
  l = len_trim( nom )

  do i = l+1,min(l+l2,Length_word) 
    nom(i:i) = mot(i-l:i-l)
  end do
      
  return
end

!***********************************************************************

subroutine ad_sufix(Core_resolved,initlr,jseuil,Length,mot,nseuil)

  use declarations
  implicit none

  integer:: initlr, jseuil, l, Length, nseuil
    
  character(len=Length):: mot
  
  Logical:: Core_resolved 

  l = len_trim(mot)

  if( Core_resolved ) then
    if( l < Length - 2 ) mot(l+1:l+1) = '_'
    call ad_number(initlr,mot,Length)
  else  
    if( l < Length - 3 ) then
      l = l + 1
      mot(l:l) = '_'
    endif
    if( l < Length - 2 ) then
      l = l + 1
      mot(l:l) = achar(nseuil+74)
    endif
    call ad_number(jseuil+initlr-1,mot,Length)
  endif

  return

end

!***********************************************************************

subroutine Spherical_tensor_cal(Abs_U_iso,Bragg_abs,ct_nelec,Core_resolved,Volume_maille,E_cut,E1E1,E1E2,E1M1,E2E2, &
            Energ,Ephseuil,Epsii,Eseuil,First_E,i_range,icheck,ie,Int_tens,jseuil,Magn_sens,Moyenne,n_rel,n_tens_max, &
            natomsym,nenerg,ninit1,ninitlr,nomfich_s,nphim,npldafs,nplr,nplrm,nplt,nseuil,numat_abs,pdp,phdf0t,phdt,pol, &
            Poldafse,Poldafss,secddia,secdqia,secmdia,secqqia,Spherical_signal,Surface_ref,Taux_eq,V0muf,Vec,Vecdafse,Vecdafss)

  use declarations
  implicit none

  integer, parameter:: n_tens_dd = 9 
  integer, parameter:: n_tens_dq = 15 
  integer, parameter:: n_tens_qq = 25 
  integer, parameter:: n_tens_dm = 9 

  integer:: i, i_range, ib, icheck, ie, initlr, ipl, ipldafs, jseuil, n_rel, n_tens_max, natomsym, nenerg, &
            ninit1, ninitlr, nphim, npldafs, nplr, nplrm, nplt, nseuil, numat_abs
  
  character(len=Length_name):: nomfich_s
  character(len=58), dimension(9):: Com_dd, com_dm, Com_dm_m
  character(len=58), dimension(15):: Com_dq, Com_dq_m
  character(len=58), dimension(25):: Com_qq

  complex(kind=db):: ph
  complex(kind=db), dimension(n_tens_dd,ninitlr):: Sph_tensor_dd
  complex(kind=db), dimension(n_tens_dq,ninitlr):: Sph_tensor_dq, Sph_tensor_dq_m
  complex(kind=db), dimension(n_tens_qq,ninitlr):: Sph_tensor_qq
  complex(kind=db), dimension(n_tens_dm,ninitlr):: Sph_tensor_dm, Sph_tensor_dm_m
  complex(kind=db), dimension(3,3,n_rel,ninitlr,0:natomsym):: secddia
  complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secmdia
  complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia
  complex(kind=db), dimension(3):: Pl_i, Pl_s
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secqqia
  complex(kind=db), dimension(3,npldafs,nphim):: Poldafse, Poldafss
  complex(kind=db), dimension(3,nplrm):: Pol
  complex(kind=db), dimension(npldafs):: phdf0t, phdt
  complex(kind=db), dimension(natomsym,npldafs):: Bragg_abs
  complex(kind=db), dimension(:,:), allocatable:: Tensor_pol_dd, Tensor_pol_dq, Tensor_pol_qq, Tensor_pol_dm

  logical:: Core_resolved, E1E1, E1E2, E1M1, E2E2, First_E, Magn_sens, Moyenne, Spherical_signal, Write_bav

  real(kind=db):: Abs_U_iso, Volume_maille, dph, E_cut, Ephseuil, Eseuil, Surface_ref, V0muf
  real(kind=db), dimension(3):: Vo_i, Vo_s
  real(kind=db), dimension(ninitlr):: ct_nelec, Epsii
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(natomsym) :: Taux_eq
  real(kind=db), dimension(3,nplrm):: Vec
  real(kind=db), dimension(nplrm,2) :: pdp
  real(kind=db), dimension(3,npldafs,nphim):: Vecdafse, Vecdafss
  real(kind=db), dimension(n_tens_max*ninitlr,0:natomsym):: Int_tens
  real(kind=db), dimension(n_tens_dd,0:natomsym,ninitlr):: Sph_tensor_ia_dd
  real(kind=db), dimension(n_tens_dq,0:natomsym,ninitlr):: Sph_tensor_ia_dq, Sph_tensor_ia_dq_m 
  real(kind=db), dimension(n_tens_qq,0:natomsym,ninitlr):: Sph_tensor_ia_qq
  real(kind=db), dimension(n_tens_dm,0:natomsym,ninitlr):: Sph_tensor_ia_dm, Sph_tensor_ia_dm_m

  data Com_dd/  ' l = 0, non-magnetic monopole:  D(00) = (Dxx+Dyy+Dzz)/3 =', &
                ' l = 1, magnetic dipole:          D(10)          = l_z  =', &
                '                                  D(11)          = l_x  =', &
                '                                 D(1-1)          = l_y  =', &
                ' l = 2, non-magnetic quadrupole:  D(20)         = d_z2  =', &
                '                                  D(21)         = d_xz  =', &
                '                                 D(2-1)         = d_yz  =', &
                '                                 D(2-2)         = d_xy  =', &
                '                                  D(22)      = d_x2-y2  ='/
  
  data Com_dq/  ' l = 1, dipole:           I(10)                 = d_z   =', &
                '                          I(11)                 = d_x   =', &
                '                         I(1-1)                 = d_y   =', &
                ' l = 2, quadrupole:       I(20)          = lz*Omega_z   =', &
                '                          I(21)          = (l,Omega)2   =', &
                '                         I(2-1)          = (l,Omega)2   =', &
                '                         I(2-2)          = (l,Omega)2   =', &
                '                          I(22)          = (l,Omega)2   =', &
                ' l = 3, octupole:         I(30)       d_z*(3lz2 - l2)   =', &
                '                          I(31)         = (n,(l,l)2)3   =', &
                '                         I(3-1)         = (n,(l,l)2)3   =', &
                '                         I(3-2)         = (n,(l,l)2)3   =', &
                '                          I(32)         = (n,(l,l)2)3   =', &
                '                         I(3-3)         = (n,(l,l)2)3   =', &
                '                          I(33)         = (n,(l,l)2)3   ='/
                                     
  data Com_dq_m/' l = 1, dipole:           I(10)               Omega_z   =', &
                '                          I(11)               Omega_x   =', &
                '                         I(1-1)               Omega_y   =', &
                ' l = 2, quadrupole:       I(20)               d_z*l_z   =', &
                '                          I(21)                (n,l)2   =', &
                '                         I(2-1)                (n,l)2   =', &
                '                         I(2-2)                (n,l)2   =', &
                '                          I(22)                (n,l)2   =', &
                ' l = 3, octupole:         I(30)   Omega_z*(3lz2 - l2)   =', &
                '                          I(31)       (Omega,(l,l)2)3   =', &
                '                         I(3-1)       (Omega,(l,l)2)3   =', &
                '                         I(3-2)       (Omega,(l,l)2)3   =', &
                '                          I(32)       (Omega,(l,l)2)3   =', &
                '                         I(3-3)       (Omega,(l,l)2)3   =', &
                '                          I(33)       (Omega,(l,l)2)3   ='/

  data Com_qq/  ' l = 0, non magn monopole:    Q(00)                     =', &
                ' l = 1, magnetic dipole:      Q(10)            l_qq_z   =', &
                '                              Q(11)            l_qq_x   =', &
                '                             Q(1-1)            l_qq_y   =', &
                ' l = 2, non magn quadrupole:  Q(20)                     =', &
                '                              Q(21)                     =', &
                '                             Q(2-1)                     =', &
                '                             Q(2-2)                     =', &
                '                              Q(22)                     =', &
                ' l = 3, magnetic octupole:    Q(30)                     =', &
                '                              Q(31)                     =', &
                '                             Q(3-1)                     =', &
                '                             Q(3-2)                     =', &
                '                              Q(32)                     =', &
                '                             Q(3-3)                     =', &
                '                              Q(33)                     =', &
                ' l = 4, non magn hexadecapol: Q(40)                     =', &
                '                              Q(41)                     =', &
                '                             Q(4-1)                     =', &
                '                             Q(4-2)                     =', &
                '                              Q(42)                     =', &
                '                             Q(4-3)                     =', &
                '                              Q(43)                     =', &
                '                             Q(4-4)                     =', &
                '                              Q(44)                     ='/

  data Com_dm/  ' l = 0, monopole:   T(00) = Im(Txx+Tyy+Tzz)/3           =', &
                ' l = 1, dipole:     T(10) = Im(Txy-Tyx )/2        = d_z =', &
                '                    T(11) = Im(Tyz-Tzy )/2        = d_x =', &
                '                    T(1-1) = Im(Tzx-Txz)/2        = d_y =', &
                ' l = 2, quadrupole: T(20) = Im(2Tzz-Txx-Tyy)/6    = dz2 =', &
                '                    T(21) = Im( Tzx + Txz )/2    = d_xz =', &
                '                    T(2-1) = Im( Tyz + Tzy )/2   = d_yz =', &
                '                    T(2-2) = Im( Txy + Tyx )/2   = d_xy =', &
                '                    T(22) = Im( Txx - Tyy )/2 = d_x2-y2 ='/

  data Com_dm_m/' l = 0, monopole:   T(00) = Re(Txx+Tyy+Tzz)/3           =', &
                ' l = 1, dipole:     T(10) = Re(Txy-Tyx )/2    = Omega_z =', &
                '                    T(11) = Re(Tyz-Tzy )/2    = Omega_x =', &
                '                    T(1-1) = Re(Tzx-Txz)/2    = Omega_y =', &
                ' l = 2, quadrupole: T(20) = Re(2Tzz-Txx-Tyy)/6          =', &
                '                    T(21) = Re( Tzx + Txz )/2           =', &
                '                    T(2-1) = Re( Tyz + Tzy )/2          =', &
                '                    T(2-2) = Re( Txy + Tyx )/2          =', &
                '                    T(22) = Re( Txx - Tyy )/2           ='/
  
! Calcul des tenseurs de l'absorption
  call Abs_Spherical_tensor(Com_dm,Com_dm_m,Com_dd,Com_dq,Com_dq_m,Com_qq,ct_nelec,E1E1,E1E2,E1M1,E2E2,First_E,icheck, &
       Magn_sens,n_rel,natomsym,ninitlr,secddia,secdqia,secmdia,secqqia, &
       Sph_tensor_ia_dd,Sph_tensor_ia_dq,Sph_tensor_ia_dq_m,Sph_tensor_ia_qq,Sph_tensor_ia_dm,Sph_tensor_ia_dm_m,Taux_eq)

  call Write_Spherical_tensor(Abs_U_iso,Core_resolved,Volume_maille,E_cut,E1E1,E1E2,E1M1,E2E2,Energ,Ephseuil,Epsii, &
       Eseuil,First_E,i_range,icheck,ie,Int_tens,jseuil,Magn_sens,n_tens_dd,n_tens_dm,n_tens_dq,n_tens_max,n_tens_qq, &
       natomsym,nenerg,ninit1,ninitlr,nomfich_s,nseuil,numat_abs,Sph_tensor_ia_dd, &
       Sph_tensor_ia_dq,Sph_tensor_ia_dq_m,Sph_tensor_ia_dm,Sph_tensor_ia_dm_m,Sph_tensor_ia_qq,Surface_ref,V0muf)

  if( .not. Spherical_signal ) return
  
 ! Calcul des tenseurs de la polarisation
 
  nplt = nplr + npldafs

  allocate( Tensor_pol_dd(0:nplt,n_tens_dd) )
  allocate( Tensor_pol_dq(0:nplt,2*n_tens_dq) ) 
  allocate( Tensor_pol_qq(0:nplt,n_tens_qq) )
  allocate( Tensor_pol_dm(0:nplt,2*n_tens_dm) ) 
  
  do ipl = 1,nplt
      
    if( ipl <= nplr ) then
      Pl_i(:) = Pol(:,ipl)
      Pl_s(:) = Pol(:,ipl)
      Vo_i(:) = Vec(:,ipl)
      Vo_s(:) = Vec(:,ipl)
    else
! Dans convolution on reprend le complexe conjugue de l'ensemble polarisation x diffusion
      Pl_i(:) = Conjg( Poldafse(:,ipl-nplr,1) )
      Pl_s(:) = Conjg( Poldafss(:,ipl-nplr,1) )
      Vo_i(:) = Vecdafse(:,ipl-nplr,1)
      Vo_s(:) = Vecdafss(:,ipl-nplr,1)
    endif

    if( E1E1 ) call Tensor_pol_dd_cal(ipl,n_tens_dd,nplt,Pl_i,Pl_s,Tensor_pol_dd)

    if( E1E2 ) call Tensor_pol_dq_cal(ipl,n_tens_dq,nplt,Pl_i,Pl_s,Vo_i,Vo_s,Tensor_pol_dq)

    if( E2E2 ) call Tensor_pol_qq_cal(ipl,n_tens_qq,nplt,Pl_i,Pl_s,Vo_i,Vo_s,Tensor_pol_qq)
  
    if( E1M1 ) call Tensor_pol_dm_cal(ipl,n_tens_dm,nplt,Pl_i,Pl_s,Vo_i,Vo_s,Tensor_pol_dm)
    
  end do

  if( Moyenne ) then
   if( E1E1 ) then
     do i = 1,n_tens_dd
       Tensor_pol_dd(0,i) = sum( pdp(1:nplr,1) * Tensor_pol_dd(1:nplr,i) )
     end do 
   endif
   if( E1E2 ) Tensor_pol_dq(0,:) = ( 0._db, 0._db )
   if( E2E2 ) then 
     do i = 1,n_tens_qq
       Tensor_pol_qq(0,i) = sum( pdp(1:nplr,2) * Tensor_pol_qq(1:nplr,i) )
     end do 
   endif
   if( E1M1 ) then
     Tensor_pol_dm(0,:) = ( 0._db, 0._db ) 
     Tensor_pol_dm(0,1) = ( 1._db, 0._db ) 
     if( Magn_sens ) Tensor_pol_dm(0,n_tens_dm+1) = ( 1._db, 0._db ) 
   endif
  endif
  
  if( ie == 1 .and. icheck > 0 ) then
    
    if( E1E1 ) then
      write(3,110) 'E1E1'
      write(3,115) ( 'Abs', i, i = 1,nplr ), ( 'Rxs', i, i = 1,npldafs )
      do i = 1,n_tens_dd
        write(3,120) i, Tensor_pol_dd(1:nplt,i)
      end do
    endif
    
    if( E1E2 ) then
      write(3,110) 'E1E2'
      write(3,115) ( 'Abs', i, i = 1,nplr ), ( 'Rxs', i, i = 1,npldafs )
      do i = 1,n_tens_dq
        write(3,120) i, Tensor_pol_dq(1:nplt,i)
      end do
      if( Magn_sens ) then
        do i = n_tens_dq+1,2*n_tens_dq
          write(3,120) i, Tensor_pol_dq(1:nplt,i)
        end do
      endif
    endif
    
    if( E2E2 ) then
      write(3,110) 'E2E2'
      write(3,115) ( 'Abs', i, i = 1,nplr ), ( 'Rxs', i, i = 1,npldafs )
      do i = 1,n_tens_qq
        write(3,120) i, Tensor_pol_qq(1:nplt,i)
      end do
    endif
    
    if( E1M1 ) then
      write(3,110) 'E1M1'
      write(3,115) ( 'Abs', i, i = 1,nplr ), ( 'Rxs', i, i = 1,npldafs )
      do i = 1,n_tens_dm
        write(3,120) i, Tensor_pol_dm(1:nplt,i)
      end do
      if( Magn_sens ) then
        do i = n_tens_dm+1,2*n_tens_dm
          write(3,120) i, Tensor_pol_dm(1:nplt,i)
        end do
      endif
    endif

  endif

  do ipl = 0,nplt
    if( .not. Moyenne .and. ipl == 0 ) cycle 
    ipldafs = ipl - nplr
      
    if( ipl <= nplr ) then
  
      do initlr = 1,ninitlr
        Sph_tensor_dd(:,initlr) = cmplx( Sph_tensor_ia_dd(:,0,initlr), 0._db ) 
        Sph_tensor_dq(:,initlr) = cmplx( Sph_tensor_ia_dq(:,0,initlr), 0._db ) 
        Sph_tensor_dq_m(:,initlr) = cmplx( 0._db, Sph_tensor_ia_dq_m(:,0,initlr) ) 
        Sph_tensor_qq(:,initlr) = cmplx( Sph_tensor_ia_qq(:,0,initlr), 0._db ) 
        Sph_tensor_dm(:,initlr) = cmplx( 0._db, Sph_tensor_ia_dm(:,0,initlr) ) 
        Sph_tensor_dm_m(:,initlr) = cmplx( Sph_tensor_ia_dm_m(:,0,initlr), 0._db )
      end do
      
    else
    
      Sph_tensor_dd(:,:) = cmplx( 0._db, 0._db ) 
      Sph_tensor_dq(:,:) = cmplx( 0._db, 0._db ) 
      Sph_tensor_dq_m(:,:) = cmplx( 0._db, 0._db ) 
      Sph_tensor_qq(:,:) = cmplx( 0._db, 0._db ) 
      Sph_tensor_dm(:,:) = cmplx( 0._db, 0._db ) 
      Sph_tensor_dm_m(:,:) = cmplx( 0._db, 0._db )

      do ib = 1,natomsym

        ph = conjg( Bragg_abs(ib,ipldafs) )

        do initlr = 1,ninitlr
          Sph_tensor_dd(:,initlr)   = Sph_tensor_dd(:,initlr)   + ph * Sph_tensor_ia_dd(:,ib,initlr) 
          Sph_tensor_dq(:,initlr)   = Sph_tensor_dq(:,initlr)   + ph * Sph_tensor_ia_dq(:,ib,initlr) 
          Sph_tensor_dq_m(:,initlr) = Sph_tensor_dq_m(:,initlr) + img * ph * Sph_tensor_ia_dq_m(:,ib,initlr) 
          Sph_tensor_qq(:,initlr)   = Sph_tensor_qq(:,initlr)   + ph * Sph_tensor_ia_qq(:,ib,initlr) 
          Sph_tensor_dm(:,initlr)   = Sph_tensor_dm(:,initlr)   + img * ph * Sph_tensor_ia_dm(:,ib,initlr) 
          Sph_tensor_dm_m(:,initlr) = Sph_tensor_dm_m(:,initlr) + ph * Sph_tensor_ia_dm_m(:,ib,initlr)
        end do
      end do
 
      if( natomsym > 1 .and. ie == 1 ) then
        Write_bav = .true.
        if( ipldafs > 1 ) then
          dph = sum( abs( Bragg_abs(:,ipldafs) - Bragg_abs(:,ipldafs-1) ) )
          if( dph < eps10 ) Write_bav = .false.
        endif
        if( Write_bav ) then
          write(3,125) ipldafs
          if( E1E1 ) then
            write(3,128) ( initlr, initlr = 1,ninitlr )
            do i = 1,n_tens_dd
              if( i == 1 .or. i == 2 .or. i == 5 ) write(3,*)
              write(3,130) Com_dd(i), Sph_tensor_dd(i,:)
            end do
          endif
          if( E1E2 ) then
            write(3,140) ( initlr, initlr = 1,ninitlr )
            do i = 1,n_tens_dq
              if( i == 1 .or. i == 4 .or. i == 9 ) write(3,*)
              write(3,130) Com_dq(i), Sph_tensor_dq(i,:)
            end do
            if( Magn_sens ) then
              write(3,150) ( initlr, initlr = 1,ninitlr )
              do i = 1,n_tens_dq
                if( i == 1 .or. i == 4 .or. i == 9 ) write(3,*)
                write(3,130) Com_dq_m(i), Sph_tensor_dq_m(i,:)
              end do
            endif
          endif
          if( E2E2 ) then
            write(3,160) ( initlr, initlr = 1,ninitlr )
            do i = 1,n_tens_qq
              if( i == 1 .or. i == 2 .or. i == 5  .or. i == 10 .or. i == 17 ) write(3,*)
              write(3,130) Com_qq(i), Sph_tensor_qq(i,:)
            end do
          endif
          if( E1M1 ) then
            write(3,170) ( initlr, initlr = 1,ninitlr )
            do i = 1,n_tens_dm
              if( i == 1 .or. i == 2 .or. i == 5 ) write(3,*)
              write(3,130) Com_dm(i), Sph_tensor_dm(i,:)
            end do
            if( Magn_sens ) then
              write(3,180) ( initlr, initlr = 1,ninitlr )
              do i = 1,n_tens_dm
                if( i == 1 .or. i == 2 .or. i == 5 ) write(3,*)
                write(3,130) Com_dm_m(i), Sph_tensor_dm_m(i,:)
              end do
            endif
          endif
        endif
      endif 
    endif 

    call Write_Signal(Abs_U_iso,ct_nelec,Core_resolved,Volume_maille,E_cut,E1E1,E1E2,E1M1,E2E2,Ephseuil,Epsii, &
      Eseuil,First_E,i_range,0,ipl,ipldafs,jseuil,Magn_sens,n_tens_dd,n_tens_dq,n_tens_dm, &
      n_tens_qq,ninit1,ninitlr,nomfich_s,npldafs,nplt,nseuil,numat_abs,phdf0t,phdt, &
      Sph_tensor_dd,Sph_tensor_dq,Sph_tensor_dq_m,Sph_tensor_dm,Sph_tensor_dm_m,Sph_tensor_qq, &
      Surface_ref,Tensor_pol_dd,Tensor_pol_dq,Tensor_pol_dm,Tensor_pol_qq,V0muf)

  end do
  
  deallocate( Tensor_pol_dd,  Tensor_pol_dq,  Tensor_pol_qq, Tensor_pol_dm ) 

  return
  110 format(/1x,A,' polarisation tensor')
  115 format(3x,25(8x,a3,i3,5x))
  120 format(i3,25(1x,2f9.5))
  125 format(/' Crystal spherical tensor (numb. of electron), reflexion',i3)
  128 format(/' E1E1 spherical tensor,                           initl =',10(14x,i2,11x))
  130 format(1p,A,10(2e13.5,1x))
  140 format(/' E1E2 spherical tensor, non-magnetic part,        initl =',10(14x,i2,11x))
  150 format(/' E1E2 spherical tensor, magnetic part,            initl =',10(14x,i2,11x))
  160 format(/' E2E2 spherical tensor,                           initl =',10(14x,i2,11x))
  170 format(/' E1M1 spherical tensor, non-magnetic part,        initl =',10(14x,i2,11x))
  180 format(/' E1M1 spherical tensor, magnetic part,            initl =',10(14x,i2,11x))
end

!***********************************************************************

subroutine Abs_Spherical_tensor(Com_dm,Com_dm_m,Com_dd,Com_dq,Com_dq_m,Com_qq,ct_nelec,E1E1,E1E2,E1M1,E2E2,First_E,icheck, &
       Magn_sens,n_rel,natomsym,ninitlr,secddia,secdqia,secmdia,secqqia, &
       Sph_tensor_ia_dd,Sph_tensor_ia_dq,Sph_tensor_ia_dq_m,Sph_tensor_ia_qq,Sph_tensor_ia_dm,Sph_tensor_ia_dm_m,Taux_eq)

  use declarations
  implicit none

  integer, parameter:: n_tens_dd = 9 
  integer, parameter:: n_tens_dq = 15 
  integer, parameter:: n_tens_qq = 25 
  integer, parameter:: n_tens_dm = 9 

  integer:: i, ia, ib, icheck, initlr, n_rel, natomsym, ninitlr
  
  character(len=58), dimension(9):: Com_dd, com_dm, Com_dm_m
  character(len=58), dimension(15):: Com_dq, Com_dq_m
  character(len=58), dimension(25):: Com_qq

  complex(kind=db), dimension(n_tens_dd):: Sph_tensor_dd
  complex(kind=db), dimension(n_tens_dq):: Sph_tensor_dq
  complex(kind=db), dimension(n_tens_qq):: Sph_tensor_qq
  complex(kind=db), dimension(n_tens_dm):: Sph_tensor_dm
  complex(kind=db), dimension(3,3):: sec_dd, sec_md
  complex(kind=db), dimension(3,3,3):: sec_dq
  complex(kind=db), dimension(3,3,3,3):: sec_qq
  complex(kind=db), dimension(3,3,n_rel,ninitlr,0:natomsym):: secddia
  complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secmdia
  complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secqqia

  logical:: E1E1, E1E2, E1M1, E2E2, First_E, Magn_sens

  real(kind=db), dimension(ninitlr):: ct_nelec
  real(kind=db), dimension(natomsym) :: Taux_eq
  real(kind=db), dimension(n_tens_dd,0:natomsym,ninitlr):: Sph_tensor_ia_dd
  real(kind=db), dimension(n_tens_dq,0:natomsym,ninitlr):: Sph_tensor_ia_dq, Sph_tensor_ia_dq_m 
  real(kind=db), dimension(n_tens_qq,0:natomsym,ninitlr):: Sph_tensor_ia_qq
  real(kind=db), dimension(n_tens_dm,0:natomsym,ninitlr):: Sph_tensor_ia_dm, Sph_tensor_ia_dm_m
  
  Sph_tensor_ia_dd(:,:,:) = 0._db
  Sph_tensor_ia_dq(:,:,:) = 0._db
  Sph_tensor_ia_dq_m(:,:,:) = 0._db
  Sph_tensor_ia_qq(:,:,:) = 0._db
  Sph_tensor_ia_dm(:,:,:) = 0._db
  Sph_tensor_ia_dm_m(:,:,:) = 0._db
  
  do initlr = 1,ninitlr

    do ia = 1,natomsym

      if( E1E1 ) then
        sec_dd(:,:) = secddia(:,:,1,initlr,ia)
        call Sph_tensor_dd_cal(n_tens_dd,sec_dd,Sph_tensor_dd)
 ! conversion en nombre d'electrons
        Sph_tensor_ia_dd(:,ia,initlr) = ct_nelec(initlr) * real( Sph_Tensor_dd(:), db )
        Sph_tensor_ia_dd(:,0,initlr) = Sph_tensor_ia_dd(:,0,initlr) + Taux_eq(ia) * Sph_tensor_ia_dd(:,ia,initlr) 
      endif

      if( E1E2 ) then
        sec_dq(:,:,:) = secdqia(:,:,:,initlr,ia)
        call Sph_tensor_dq_cal(n_tens_dq,sec_dq,Sph_tensor_dq)
        Sph_tensor_ia_dq(:,ia,initlr) = ct_nelec(initlr) * real( Sph_Tensor_dq(:), db )
        if( Magn_sens ) Sph_tensor_ia_dq_m(:,ia,initlr) = ct_nelec(initlr) * aimag( Sph_Tensor_dq(:) )
        Sph_tensor_ia_dq(:,0,initlr) = Sph_tensor_ia_dq(:,0,initlr) + Taux_eq(ia) * Sph_tensor_ia_dq(:,ia,initlr)
        if( Magn_sens ) Sph_tensor_ia_dq_m(:,0,initlr) = Sph_tensor_ia_dq_m(:,0,initlr) &
                                                       + Taux_eq(ia) * Sph_tensor_ia_dq_m(:,ia,initlr) 
      endif

      if( E2E2 ) then
        sec_qq(:,:,:,:) = secqqia(:,:,:,:,initlr,ia)
        call Sph_tensor_qq_cal(n_tens_qq,sec_qq,Sph_tensor_qq)
        Sph_tensor_ia_qq(:,ia,initlr) = ct_nelec(initlr) * real( Sph_Tensor_qq(:), db )
        Sph_tensor_ia_qq(:,0,initlr) = Sph_tensor_ia_qq(:,0,initlr) + Taux_eq(ia) * Sph_tensor_ia_qq(:,ia,initlr)
      endif
     
      if( E1M1 ) then  ! pour E1M1, la partie magnetique est le terme reel
        sec_md(:,:) = secmdia(:,:,initlr,ia)
        call Sph_tensor_dm_cal(n_tens_dm,sec_md,Sph_tensor_dm)
        Sph_tensor_ia_dm(:,ia,initlr) = ct_nelec(initlr) * aimag( Sph_tensor_dm(:) )
        if( Magn_sens ) Sph_tensor_ia_dm_m(:,ia,initlr) = ct_nelec(initlr) * real( Sph_tensor_dm(:), db )
        Sph_tensor_ia_dm(:,0,initlr) = Sph_tensor_ia_dm(:,0,initlr) + Taux_eq(ia) * Sph_tensor_ia_dm(:,ia,initlr)
        if( Magn_sens ) Sph_tensor_ia_dm_m(:,0,initlr) = Sph_tensor_ia_dm_m(:,0,initlr) &
                                                       + Taux_eq(ia) * Sph_tensor_ia_dm_m(:,ia,initlr) 
      endif

    end do  ! fin boucle ia

  end do ! fin boucle initlr

  do ib = 1,natomsym+1
    if( icheck == 0 .or. .not. First_E ) exit
    if( ib == natomsym+1 ) then
      if( natomsym == 1 ) exit
      ia = 0
      write(3,110)
    else
      ia = ib
      if( ia > 1 .and. icheck < 2 ) cycle
      write(3,120) ia
    endif

    if( E1E1 ) then
      write(3,128) ( initlr, initlr = 1,ninitlr )
      do i = 1,n_tens_dd
        if( i == 1 .or. i == 2 .or. i == 5 ) write(3,*)
        write(3,130) Com_dd(i), Sph_tensor_ia_dd(i,ia,:)
      end do
    endif
    if( E1E2 ) then
      write(3,140) ( initlr, initlr = 1,ninitlr )
      do i = 1,n_tens_dq
        if( i == 1 .or. i == 4 .or. i == 9 ) write(3,*)
        write(3,130) Com_dq(i), Sph_tensor_ia_dq(i,ia,:)
      end do
      if( Magn_sens ) then
        write(3,150) ( initlr, initlr = 1,ninitlr )
        do i = 1,n_tens_dq
          if( i == 1 .or. i == 4 .or. i == 9 ) write(3,*)
          write(3,130) Com_dq_m(i), Sph_tensor_ia_dq_m(i,ia,:)
        end do
      endif
    endif
    if( E2E2 ) then
      write(3,160) ( initlr, initlr = 1,ninitlr )
      do i = 1,n_tens_qq
        if( i == 1 .or. i == 2 .or. i == 5  .or. i == 10 .or. i == 17 ) write(3,*)
        write(3,130) Com_qq(i), Sph_tensor_ia_qq(i,ia,:)
      end do
    endif
    if( E1M1 ) then
      write(3,170) ( initlr, initlr = 1,ninitlr )
      do i = 1,n_tens_dm
        if( i == 1 .or. i == 2 .or. i == 5 ) write(3,*)
        write(3,130) Com_dm(i), Sph_tensor_ia_dm(i,ia,:)
      end do
      if( Magn_sens ) then
        write(3,180) ( initlr, initlr = 1,ninitlr )
        do i = 1,n_tens_dm
          if( i == 1 .or. i == 2 .or. i == 5 ) write(3,*)
          write(3,130) Com_dm_m(i), Sph_tensor_ia_dm_m(i,ia,:)
        end do
      endif
    endif

  end do

  return
  110 format(/' Crystal spherical tensor (numb. of electron)')
  120 format(/' Spherical tensor (numb. of electron), atome =',i3)
  128 format(/' E1E1 spherical tensor,                           initl =',10(7x,i2,5x))
  130 format(1p,A,10(e13.5,1x))
  140 format(/' E1E2 spherical tensor, non-magnetic part,        initl =',10(7x,i2,5x))
  150 format(/' E1E2 spherical tensor, magnetic part,            initl =',10(7x,i2,5x))
  160 format(/' E2E2 spherical tensor,                           initl =',10(7x,i2,5x))
  170 format(/' E1M1 spherical tensor, non-magnetic part,        initl =',10(7x,i2,5x))
  180 format(/' E1M1 spherical tensor, magnetic part,            initl =',10(7x,i2,5x))
end

!***********************************************************************

subroutine Sph_tensor_dd_cal(n_tens_dd,sec_dd,Sph_tensor_dd)

  use declarations
  implicit none
  
  integer:: n_tens_dd

  complex(kind=db), dimension(n_tens_dd):: Sph_tensor_dd
  complex(kind=db), dimension(3,3):: sec_dd

! l = 0

  Sph_tensor_dd(1) = ( 1 / 3._db ) * ( sec_dd(1,1) + sec_dd(2,2) + sec_dd(3,3) )

! l = 1

! Multiplie par - img, equivalent a prendre la partie imaginaire 
! lz = D(1,0)
  Sph_tensor_dd(2) = - img * ( sec_dd(1,2) - sec_dd(2,1) ) / 2

! lx = D(1,1)
  Sph_tensor_dd(3) = - img * ( sec_dd(2,3) - sec_dd(3,2) ) / 2

! ly = D(1,-1)
  Sph_tensor_dd(4) = - img * ( sec_dd(3,1) - sec_dd(1,3) ) / 2

! l = 2

! D02 = ( 2*Dzz - Dxx - Dyy ) / 6
  Sph_tensor_dd(5) =( 2*sec_dd(3,3) - sec_dd(1,1)  - sec_dd(2,2) ) / 6

! D(2,1)
  Sph_tensor_dd(6) = ( sec_dd(1,3) + sec_dd(3,1) ) / 2

! D(2,-1)
  Sph_tensor_dd(7) = ( sec_dd(2,3) + sec_dd(3,2) ) / 2

! D(2,-2)
  Sph_tensor_dd(8) = ( sec_dd(1,2) + sec_dd(2,1) ) / 2

! D(2,2)
  Sph_tensor_dd(9) = ( sec_dd(1,1) - sec_dd(2,2) ) / 2

  return
end

!***********************************************************************

subroutine Sph_tensor_dq_cal(n_tens_dq,sec_dq,Sph_tensor_dq)

  use declarations
  implicit none
  
  integer:: n_tens_dq

  complex(kind=db), dimension(3,3,3):: sec_dq
  complex(kind=db), dimension(n_tens_dq):: Sph_tensor_dq

  real(kind=db):: fac
  
! Tenseur 1

! I(10)
  Sph_tensor_dq(1) = ( 2 * sec_dq(3,3,3) - sec_dq(3,1,1) - sec_dq(3,2,2) + 3 * sec_dq(1,1,3) + 3 * sec_dq(2,2,3) ) / 15

! I(11)
  Sph_tensor_dq(2) = ( 2 * sec_dq(1,1,1) - sec_dq(1,2,2) - sec_dq(1,3,3) + 3 * sec_dq(2,1,2) + 3 * sec_dq(3,1,3) ) / 15

! I(1-1)
  Sph_tensor_dq(3) = ( 2 * sec_dq(2,2,2) - sec_dq(2,1,1) - sec_dq(2,3,3) + 3 * sec_dq(1,1,2) + 3 * sec_dq(3,2,3) ) / 15

! Tenseur 2

! -i*I(20)
  Sph_tensor_dq(4) = sec_dq(1,2,3) - sec_dq(2,1,3)

  fac = 1 / 3._db

! (-i/sqrt(2))*(I(21) - I(2-1))
  Sph_tensor_dq(5) = fac * ( sec_dq(2,1,1) - sec_dq(2,3,3) - sec_dq(1,1,2) + sec_dq(3,2,3) )

! (1/sqrt(2))*(I(21) + I(2-1))
  Sph_tensor_dq(6) = fac * ( sec_dq(1,2,2) - sec_dq(1,3,3) - sec_dq(2,1,2) + sec_dq(3,1,3) )

  Sph_tensor_dq(7) = fac * ( sec_dq(1,1,3) - sec_dq(2,2,3) - sec_dq(3,1,1) + sec_dq(3,2,2) )

! (-i/sqrt(2))*(I(22) + I(2-2))
  Sph_tensor_dq(8) = fac * ( sec_dq(1,2,3) + sec_dq(2,1,3) - sec_dq(3,1,2) - sec_dq(3,2,1) )

! Tenseur 3

  Sph_tensor_dq(9) = 1 / sqrt(10._db) * ( 2 * sec_dq(3,3,3) - 2 * sec_dq(1,1,3) - sec_dq(3,1,1) &
                   - 2 * sec_dq(2,2,3) - sec_dq(3,2,2) )

  fac = 1 / sqrt(60._db)

  Sph_tensor_dq(10) = fac * ( 3 * sec_dq(1,1,1) + sec_dq(1,2,2) + 2 * sec_dq(2,1,2) - 4 * sec_dq(1,3,3) - 8 * sec_dq(3,1,3) )

! (-i/sqrt(2))*(I(31) + I(3-1))
  Sph_tensor_dq(11) = fac * ( 3 * sec_dq(2,2,2) + sec_dq(2,1,1) + 2 * sec_dq(1,2,1) - 4 * sec_dq(2,3,3) - 8 * sec_dq(3,2,3) )

  fac = 1 / sqrt(6._db)

! (-i/sqrt(2))*(I(32) - I(3-2))
  Sph_tensor_dq(12) = 2 * fac * ( sec_dq(1,2,3) + sec_dq(2,1,3) + sec_dq(3,1,2) )

  Sph_tensor_dq(13) = fac * ( 2 * sec_dq(1,1,3) - 2 * sec_dq(2,2,3) + sec_dq(3,1,1) - sec_dq(3,2,2) )

  fac = 0.5_db

  Sph_tensor_dq(14) = fac * ( sec_dq(1,2,2) + 2 * sec_dq(2,1,2) - sec_dq(1,1,1) )

! (-i/sqrt(2))*(I(33) + I(3-3))
  Sph_tensor_dq(15) = - fac * ( sec_dq(2,1,1) + 2 * sec_dq(1,2,1) - sec_dq(2,2,2) )

  return
end

!***********************************************************************

subroutine Sph_tensor_qq_cal(n_tens_qq,sec_qq,Sph_tensor_qq)

  use declarations
  implicit none
  
  integer:: n_tens_qq

  complex(kind=db), dimension(3,3,3,3):: sec_qq
  complex(kind=db), dimension(n_tens_qq):: Sph_tensor_qq

  real(kind=db):: fac
  
! La multiplication par - img correspond a la partie imaginaire quand le
! tenseur n'est pas multiplie par le terme de Bragg.

! Tenseur 0, scalaire, signal isotropique quadrupolaire

  fac = 1 / 45._db
! Q(00)
  Sph_tensor_qq(1) = fac * ( 6 * ( sec_qq(1,3,1,3) + sec_qq(2,3,2,3) + sec_qq(1,2,1,2) ) &
                           + 2 * ( sec_qq(1,1,1,1) + sec_qq(2,2,2,2) + sec_qq(3,3,3,3) ) &
                 - ( sec_qq(1,1,2,2) + sec_qq(1,1,3,3) + sec_qq(3,3,2,2) + sec_qq(2,2,1,1) + sec_qq(3,3,1,1) + sec_qq(2,2,3,3) ) )

! Tenseur 1, vecteur, lz, lx, ly, magnetique

  fac = 2 / 10._db

! Q(10)
  Sph_tensor_qq(2) = - img * fac * ( sec_qq(2,1,1,1) - sec_qq(1,2,2,2) + sec_qq(2,3,1,3) &
                                   - sec_qq(1,1,2,1) + sec_qq(2,2,1,2) - sec_qq(1,3,2,3) )

! Q(11)
  Sph_tensor_qq(3) = - img * fac * ( sec_qq(3,2,2,2) - sec_qq(2,3,3,3) + sec_qq(3,1,2,1) &
                                   - sec_qq(2,2,3,2) + sec_qq(3,3,2,3) - sec_qq(2,1,3,1) )

! Q(1-1)
  Sph_tensor_qq(4) = - img * fac * ( sec_qq(1,3,3,3) - sec_qq(3,1,1,1) + sec_qq(1,2,3,2) &
                                   - sec_qq(3,3,1,3) + sec_qq(1,1,3,1) - sec_qq(3,2,1,2) )

! Tenseur 2 : quadrupole non magnetique

  fac = 2 / ( 3 * sqrt(14._db) )
! Q(20)
  Sph_tensor_qq(5) = fac * ( 6 * sec_qq(1,2,1,2) - 3 * sec_qq(1,3,1,3) - 3 * sec_qq(2,3,2,3) + sec_qq(1,1,1,1) + sec_qq(2,2,2,2) &
           - 2 * sec_qq(1,1,2,2) - 2 * sec_qq(2,2,1,1) - 2 * sec_qq(3,3,3,3) + sec_qq(1,1,3,3) + sec_qq(2,2,3,3) &
           + sec_qq(3,3,1,1) + sec_qq(3,3,2,2) )

  fac = 2 / sqrt(42._db)

! (1/sqrt(2)) * ( Q(21) - Q(2-1) )_comp = Q(21) reel
  Sph_tensor_qq(6) = fac * ( 3 * sec_qq(2,3,1,2) + sec_qq(3,3,1,3) + sec_qq(1,1,1,3) - 2 * sec_qq(2,2,1,3) &
                     + 3 * sec_qq(1,2,2,3) + sec_qq(1,3,3,3) + sec_qq(1,3,1,1) - 2 * sec_qq(1,3,2,2) )

! (1/sqrt(2)) * ( Q(21) + Q(2-1) )
  Sph_tensor_qq(7) = fac * ( 3 * sec_qq(1,3,1,2) + sec_qq(3,3,2,3) + sec_qq(2,2,2,3) - 2 * sec_qq(1,1,2,3) &
                     + 3 * sec_qq(1,2,1,3) + sec_qq(2,3,3,3) + sec_qq(2,3,2,2) - 2 * sec_qq(2,3,1,1) )

  fac = 2 / sqrt( 42._db )

! (1/sqrt(2)) * ( Q(22) - Q(2-2) )
  Sph_tensor_qq(8) = fac * ( 2 * sec_qq(3,3,1,2) - sec_qq(1,1,1,2) - sec_qq(2,2,2,1) - 3 * sec_qq(1,3,2,3) &
             + 2 * sec_qq(1,2,3,3) - sec_qq(1,2,1,1) - sec_qq(2,1,2,2) - 3 * sec_qq(2,3,1,3) )

! (1/sqrt(2)) * ( Q(22) + Q(2-2) )
  Sph_tensor_qq(9) = fac * ( sec_qq(3,3,1,1) - sec_qq(3,3,2,2) + sec_qq(1,1,3,3) - sec_qq(2,2,3,3) &
             + sec_qq(2,2,2,2) - sec_qq(1,1,1,1) + 3 * sec_qq(2,3,2,3) - 3 * sec_qq(1,3,1,3) )

! Tenseur 3 : octupole magnetique

  fac = 1 / sqrt(10._db)

! Q(30)
  Sph_tensor_qq(10) = - img * fac * ( sec_qq(1,2,1,1) - sec_qq(2,1,2,2) + 4 * sec_qq(1,3,2,3) &
           - sec_qq(1,1,1,2) + sec_qq(2,2,2,1) - 4 * sec_qq(2,3,1,3) )

  fac = 0.5_db / sqrt(15._db)

! (1/sqrt(2)) * ( Q(31) - Q(3-1) )
  Sph_tensor_qq(11) = - img *  fac * ( 6 * sec_qq(1,2,1,3) + 4 * sec_qq(3,3,2,3) - 5 * sec_qq(1,1,2,3) + sec_qq(2,2,2,3) &
                      - 6 * sec_qq(1,3,1,2) - 4 * sec_qq(2,3,3,3) + 5 * sec_qq(2,3,1,1) - sec_qq(2,3,2,2) )

! (1/sqrt(2)) * ( Q(31) + Q(3-1) )
  Sph_tensor_qq(12) = - img * fac * ( 6 * sec_qq(1,2,2,3) + 4 * sec_qq(3,3,1,3) - 5 * sec_qq(2,2,1,3) + sec_qq(1,1,1,3) &
                      - 6 * sec_qq(2,3,1,2) - 4 * sec_qq(1,3,3,3) + 5 * sec_qq(1,3,2,2) - sec_qq(1,3,1,1) )

  fac = 1 / sqrt(6._db)

! (1/sqrt(2)) * ( Q(32) - Q(3-2) )
  Sph_tensor_qq(13) = - img * fac * ( sec_qq(1,1,3,3) + sec_qq(2,2,1,1) + sec_qq(3,3,2,2) &
               - sec_qq(3,3,1,1) - sec_qq(1,1,2,2) - sec_qq(2,2,3,3) )

! (1/sqrt(2)) * ( Q(32) + Q(3-2) )
  Sph_tensor_qq(14) = - img * fac * ( 2 * sec_qq(1,2,3,3) - sec_qq(1,2,1,1) - sec_qq(1,2,2,2) &
           - 2 * sec_qq(3,3,1,2) + sec_qq(1,1,1,2) + sec_qq(2,2,1,2) )

  fac = 0.5_db

! (1/sqrt(2)) * ( Q(33) - Q(3-3) )
  Sph_tensor_qq(15) = - img * fac * ( sec_qq(2,3,1,1) - sec_qq(2,3,2,2) + 2 * sec_qq(1,3,1,2) &
           - sec_qq(1,1,2,3) + sec_qq(2,2,2,3) - 2 * sec_qq(1,2,1,3) )

! (1/sqrt(2)) * ( Q(33) + Q(3-3) )
  Sph_tensor_qq(16) = - img * fac * ( sec_qq(1,3,1,1) - sec_qq(1,3,2,2) + 2 * sec_qq(2,1,3,2) &
           - sec_qq(1,1,1,3) + sec_qq(2,2,1,3) - 2 * sec_qq(3,2,2,1) )

! Tenseur 4 : hexadecapole non magnetique

  fac = 1 / ( 2. * sqrt(70._db) )

! Q(40)
  Sph_tensor_qq(17) = fac * ( 3 * sec_qq(1,1,1,1) + 3 * sec_qq(2,2,2,2) &
           + 8 * sec_qq(3,3,3,3) + sec_qq(1,1,2,2) + sec_qq(2,2,1,1) - 4 * sec_qq(3,3,1,1) - 4 * sec_qq(3,3,2,2) &
           - 4 * sec_qq(1,1,3,3) - 4 * sec_qq(2,2,3,3) + 4 * sec_qq(1,2,1,2) - 16 * sec_qq(1,3,1,3) - 16 * sec_qq(2,3,2,3) )

  fac = 0.5_db / sqrt(7._db)

! (1/sqrt(2)) * ( Q(41) - Q(4-1) )
  Sph_tensor_qq(18) = fac * ( 3 * sec_qq(1,3,1,1) + sec_qq(1,3,2,2) - 4 * sec_qq(3,3,1,3) + 2 * sec_qq(2,3,1,2) &
          + 3 * sec_qq(1,1,1,3) + sec_qq(2,2,1,3) - 4 * sec_qq(1,3,3,3) + 2 * sec_qq(1,2,2,3) )

! (1/sqrt(2)) * ( Q(41) + Q(4-1) )
  Sph_tensor_qq(19) = fac * ( 3 * sec_qq(2,3,2,2) + sec_qq(2,3,1,1) - 4 * sec_qq(2,3,3,3) + 2 * sec_qq(1,3,1,2) &
          + 3 * sec_qq(2,2,2,3) + sec_qq(1,1,2,3) - 4 * sec_qq(3,3,2,3) + 2 * sec_qq(1,2,1,3) )

  fac = 0.5_db / sqrt(14._db)

! (1/sqrt(2)) * ( Q(42) - Q(4-2) )
  Sph_tensor_qq(20) = 2 * fac * ( 2 * sec_qq(3,3,1,2) - sec_qq(1,1,1,2) - sec_qq(2,2,2,1) + 4 * sec_qq(1,3,2,3) &
          + 2 * sec_qq(1,2,3,3) - sec_qq(1,2,1,1) - sec_qq(2,1,2,2) + 4 * sec_qq(2,3,1,3) )

! (1/sqrt(2)) * ( Q(42) + Q(4-2) )
  Sph_tensor_qq(21) = fac * ( 2 * sec_qq(3,3,1,1) - 2 * sec_qq(3,3,2,2) + 2 * sec_qq(1,1,3,3) - 2 * sec_qq(2,2,3,3) &
        - 2 * sec_qq(1,1,1,1) + 2 * sec_qq(2,2,2,2) + 8 * sec_qq(1,3,1,3) - 8 * sec_qq(2,3,2,3) )

  fac = 0.5_db

! (1/sqrt(2)) * ( Q(43) - Q(4-3) )
  Sph_tensor_qq(22) = fac * ( sec_qq(1,3,2,2) - sec_qq(1,3,1,1) + 2 * sec_qq(2,3,1,2) &
           + sec_qq(2,2,1,3) - sec_qq(1,1,1,3) + 2 * sec_qq(1,2,2,3) )

! (1/sqrt(2)) * ( Q(43) + Q(4-3) )
  Sph_tensor_qq(23) = fac * ( sec_qq(2,3,2,2) - sec_qq(2,3,1,1) - 2 * sec_qq(1,3,1,2) &
           + sec_qq(2,2,2,3) - sec_qq(1,1,2,3) - 2 * sec_qq(1,2,1,3) )

  fac = 1 / sqrt(2._db)

! (1/sqrt(2)) * ( Q(44) - Q(4-4) )
  Sph_tensor_qq(24) = fac * ( sec_qq(1,2,1,1) - sec_qq(1,2,2,2) + sec_qq(1,1,1,2) - sec_qq(2,2,1,2) )

! (1/sqrt(2)) * ( Q(44) + Q(4-4) )
  Sph_tensor_qq(25) = 0.5_db * fac * ( sec_qq(1,1,1,1) + sec_qq(2,2,2,2) - sec_qq(1,1,2,2) &
             - sec_qq(2,2,1,1) - 4 * sec_qq(1,2,1,2) )

  return
end

!***********************************************************************

subroutine Sph_tensor_dm_cal(n_tens_dm,sec_md,Sph_tensor_dm)

  use declarations
  implicit none
  
  integer:: n_tens_dm

  complex(kind=db), dimension(n_tens_dm):: Sph_tensor_dm
  complex(kind=db), dimension(3,3):: sec_md

  real(kind=db):: fac
  
! Tenseur 0

  Sph_tensor_dm(1) = ( 1 / 3._db ) * ( sec_md(1,1) + sec_md(2,2) + sec_md(3,3) )

! Tenseur 1
! Les composantes de ce tenseur sont en cas de seuil K : -lx, ly et lz.
  fac = 1 / 2._db

! Omega_z
  Sph_tensor_dm(2) = fac * ( sec_md(1,2) - sec_md(2,1) )

! Omega_x
  Sph_tensor_dm(3) = fac * ( sec_md(2,3) - sec_md(3,2) )

! Omega_y
  Sph_tensor_dm(4) = fac * ( sec_md(3,1) - sec_md(1,3) )

! Tenseur 2

! D02 = (1/sqrt(6))*(2*Dzz-Dxx-Dyy) = T_3z2-r2
  fac = 1 / 6._db
  Sph_tensor_dm(5) = fac * ( 2*sec_md(3,3) - sec_md(1,1)  - sec_md(2,2) )

  fac = 1 / 2._db

! (-1/sqrt(2))*(D(12) - D(-12)) = T_xz
  Sph_tensor_dm(6) = fac * ( sec_md(1,3) + sec_md(3,1) )

! (i/sqrt(2))*(D(12) + D(-12)) = T_yz
  Sph_tensor_dm(7) = fac * ( sec_md(2,3) + sec_md(3,2) )

! (-i/sqrt(2))*(D(22) - D(-22)) = T_xy
  Sph_tensor_dm(8) = fac * ( sec_md(1,2) + sec_md(2,1) )

! (1/sqrt(2))*(D(22) + D(-22)) = T_x2-y2
  Sph_tensor_dm(9) = fac * ( sec_md(1,1) - sec_md(2,2) )

! Les tenseurs cartesiens sont definis par conjg(D).M, il faut conjg(M).D 
  Sph_tensor_dm(2:4) = - Sph_tensor_dm(2:4)
  Sph_tensor_dm(:) = Conjg( Sph_tensor_dm(:) ) 

  return
end

!***********************************************************************

! Ecriture des tenseurs spheriques et de leur integrale

subroutine Write_Spherical_tensor(Abs_U_iso,Core_resolved,Volume_maille,E_cut,E1E1,E1E2,E1M1,E2E2,Energ,Ephseuil,Epsii, &
       Eseuil,First_E,i_range,icheck,ie,Int_tens,jseuil,Magn_sens,n_tens_dd,n_tens_dm,n_tens_dq,n_tens_max,n_tens_qq, &
       natomsym,nenerg,ninit1,ninitlr,nomfich_s,nseuil,numat_abs,Sph_tensor_ia_dd, &
       Sph_tensor_ia_dq,Sph_tensor_ia_dq_m,Sph_tensor_ia_dm,Sph_tensor_ia_dm_m,Sph_tensor_ia_qq,Surface_ref,V0muf)

  use declarations
  implicit none
  
  integer:: i, i_range, ia, ie, icheck, initlr, j, jseuil, long, n_tens, n_tens_dd, n_tens_dq, n_tens_max, n_tens_dm, &
      n_tens_qq, natomsym, nenerg, ninit1, ninitlr, nseuil, numat_abs

  character(len=Length_word):: mot
  character(len=Length_name):: nomfich_s, nomficht
  character(len=8), dimension(9):: Tens_name_D
  character(len=8), dimension(15):: Tens_name_I, Tens_name_I_m 
  character(len=8), dimension(25):: Tens_name_Q 
  character(len=8), dimension(9):: Tens_name_T, Tens_name_T_m 
  character(len=Length_word), dimension(n_tens_max*ninitlr):: nomten

  complex(kind=db):: zero_c
  complex(kind=db), dimension(1):: cdum

  logical:: Core_resolved, E1E1, E1E2, E1M1, E2E2, First_E, Magn_sens

  real(kind=db):: Abs_U_iso, de, Volume_maille, E_cut, Ephseuil, Eseuil, Surface_ref, V0muf
  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(n_tens_max*ninitlr,0:natomsym):: Int_tens
  real(kind=db), dimension(n_tens_max*ninitlr):: Int_tenst, Tens
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(ninitlr):: Epsii
  real(kind=db), dimension(n_tens_dd,0:natomsym,ninitlr):: Sph_tensor_ia_dd
  real(kind=db), dimension(n_tens_dq,0:natomsym,ninitlr):: Sph_tensor_ia_dq, Sph_tensor_ia_dq_m 
  real(kind=db), dimension(n_tens_qq,0:natomsym,ninitlr):: Sph_tensor_ia_qq
  real(kind=db), dimension(n_tens_dm,0:natomsym,ninitlr):: Sph_tensor_ia_dm, Sph_tensor_ia_dm_m

  data Tens_name_D/ '  D(00) ','  lz_dd ','  lx_dd ','  ly_dd ','  D_z2  ','   D_xz ','  D_yz  ','  D_xy  ',' D_x2-y2'/
  data Tens_name_I/ '   I_z  ','   I_x  ','    I_y ','  I(20) ',' I(21)  ',' I(2-1) ','  I(22) ',' I(2-2) ', &
                    '  I(30) ',' I(31)  ',' I(3-1) ',' I(32)  ',' I(3-2) ','  I(33) ',' I(3-3) '/
  data Tens_name_I_m/ '  I_z_m ',' I_x_m  ','  I_y_m ',' I(20)_m',' I(21)_m','I(2-1)_m',' I(22)_m','I(2-2)_m', &
                      ' I(30)_m',' I(31)_m','I(3-1)_m',' I(32)_m','I(3-2)_m',' I(33)_m','I(3-3)_m'/
  data Tens_name_Q/ '  Q(00) ','  lz_qq ','  lx_qq ','  ly_qq ','  Q(20) ','  Q(21) ',' Q(2-1) ','  Q(22) ',' Q(2-2) ', &
                    '  Q(30) ','  Q(31) ',' Q(3-1) ','  Q(32) ',' Q(3-2) ','  Q(33) ',' Q(3-3) ','  Q(40) ','  Q(41) ', &
                    ' Q(4-1) ','  Q(42) ',' Q(4-2) ','  Q(43) ',' Q(4-3) ','  Q(44) ',' Q(4-4) '/
  data Tens_name_T/   '  T(00) ','  T_z   ','   T_x  ','    T_y ','  T_z2  ','  T_xz  ','  T_yz  ','  T_yz  ','T_x2-y2 '/
  data Tens_name_T_m/ ' T(00)_m','  T_z_m ','  T_x_m ','  T_y_m ',' T_z2_m ',' T_xz_m ',' T_yz_m ',' T_yz_m ','Tx2-y2_m'/

  zero_c = (0._db, 0._db)

  do ia = 0,natomsym
 
    if( natomsym == 1 .and. ia == 0 ) cycle
    if( icheck <= 1 .and. ia > 1 ) exit
    
    j = 0
     
    do initlr = 1,ninitlr
  
      if( E1E1 ) then
        do i = 1,n_tens_dd
          j = j + 1
          select case(i)
            case(1)
! On divise la premiere composante du tenseur spherique par rac(3) (premiere composante du tenseur de polarisation) pour obtenir le
! terme de diffusion isotrope
              Tens(j) = Sph_tensor_ia_dd(i,ia,initlr)
            case(2,3,4)
! On divise les composantes 2, 3 et 4 du tenseur spherique par rac(2) (composantes du tenseur de polarisation) pour obtenir le
! moment magnetique
              Tens(j) = Sph_tensor_ia_dd(i,ia,initlr)
            case default
              Tens(j) = Sph_tensor_ia_dd(i,ia,initlr)
          end select
          mot = Tens_name_D(i) 
          if( ninitlr > 1 ) call ad_sufix(Core_resolved,initlr,jseuil,Length_word,mot,nseuil)
          nomten(j) = mot
        end do 
      endif

      if( E1E2 ) then
        do i = 1,n_tens_dq
          j = j + 1
          Tens(j) = Sph_tensor_ia_dq(i,ia,initlr)
          mot = Tens_name_I(i) 
          if( ninitlr > 1 ) call ad_sufix(Core_resolved,initlr,jseuil,Length_word,mot,nseuil)
          nomten(j) = mot
          if( Magn_sens ) then
            j = j+1
            Tens(j) = Sph_tensor_ia_dq_m(i,ia,initlr)
            mot = Tens_name_I_m(i) 
            if( ninitlr > 1 ) call ad_sufix(Core_resolved,initlr,jseuil,Length_word,mot,nseuil)
            nomten(j) = mot
          endif
        end do 
      endif
    
      if( E2E2 ) then
        do i = 1,n_tens_qq
          j = j + 1
          select case(i)
            case(1)
! On divise la premiere composante du tenseur spherique par rac(3) (premiere composante du tenseur de polarisation) pour obtenir le
! terme de diffusion isotrope
              Tens(j) = Sph_tensor_ia_qq(i,ia,initlr)
            case(2,3,4)
! On divise les composantes 2, 3 et 4 du tenseur spherique par rac(2) (composantes du tenseur de polarisation) pour obtenir le
! moment magnetique
              Tens(j) = Sph_tensor_ia_qq(i,ia,initlr)
            case default
              Tens(j) = Sph_tensor_ia_qq(i,ia,initlr)
          end select
          mot = Tens_name_Q(i) 
          if( ninitlr > 1 ) call ad_sufix(Core_resolved,initlr,jseuil,Length_word,mot,nseuil)
          nomten(j) = mot
       end do 
      endif

      if( E1M1 ) then
        do i = 1,n_tens_dm
          j = j + 1
          Tens(j) = Sph_tensor_ia_dm(i,ia,initlr)
          mot = Tens_name_T(i) 
          if( ninitlr > 1 ) call ad_sufix(Core_resolved,initlr,jseuil,Length_word,mot,nseuil)
          nomten(j) = mot
          if( Magn_sens ) then
            j = j+1
            Tens(j) = Sph_tensor_ia_dm_m(i,ia,initlr)
            mot = Tens_name_T_m(i) 
            if( ninitlr > 1 ) call ad_sufix(Core_resolved,initlr,jseuil,Length_word,mot,nseuil)
            nomten(j) = mot
          endif
        end do 
      endif
     
    end do

    n_tens = j
    
    nomficht = nomfich_s
    long = len_trim(nomficht)
    nomficht(long+1:long+4) = '_sph'
    if( ia > 0 ) then
      nomficht(long+5:long+9) = '_atom'
      call ad_number(ia,nomficht,Length_name)
    else
      nomficht(long+5:long+9) = '_xtal'
    endif
    long = len_trim(nomficht)
    nomficht(long+1:long+4) = '.txt'

    call write_out(Abs_U_iso,rdum,rdum,zero_c,E_cut,Ephseuil,Epsii,Eseuil,First_E,.false.,rdum, &
          i_range,jseuil,n_tens_max*ninitlr,n_tens,ninit1,ninitlr,nomficht,nomten,1,0,0,0,nseuil,numat_abs, &
          cdum,cdum,Tens,V0muf,Core_resolved,0._db,Surface_ref,Volume_maille)

    if( nenerg == 1 ) cycle

! Integrale
  
    nomficht = nomfich_s
    long = len_trim(nomficht)
    nomficht(long+1:long+4) = '_sph'
    if( ia > 0 ) then
      nomficht(long+5:long+9) = '_atom'
      call ad_number(ia,nomficht,Length_name)
    else
      nomficht(long+5:long+9) = '_xtal'
    endif
    long = len_trim(nomficht)
    nomficht(long+1:long+8) = '_int.txt'

    if( First_E ) then
      do i = 1,n_tens
        mot = nomten(i)
        mot = adjustl( mot )
        long = len_trim( mot )
        if( long < Length_word - 2 ) then
          mot(1:Length_word) = 'I_' // mot(1:Length_word-2)
        elseif( long == Length_word - 2 ) then
          mot(1:Length_word-1) = 'I' // mot(1:Length_word-1)
        endif
        nomten(i) = mot
      end do
      de = Energ(ie+1) - Energ(ie)
      Int_tens(1:n_tens,ia) = de * Tens(1:n_tens)
    else
      if( ie == nenerg ) then
        de = Energ(ie) - Energ(ie-1)
      else
        de = 0.5_db * ( Energ(ie+1) -  Energ(ie-1) )
      endif
      Int_tens(1:n_tens,ia) = Int_tens(1:n_tens,ia) + de * Tens(1:n_tens)
    endif

    Int_tenst(1:n_tens) = Int_tens(1:n_tens,ia)
    call write_out(Abs_U_iso,rdum,rdum,zero_c,E_cut,Ephseuil,Epsii,Eseuil,First_E,.false.,rdum, &
          i_range,jseuil,n_tens_max*ninitlr,n_tens,ninit1,ninitlr,nomficht,nomten,1,0,0,0,nseuil,numat_abs, &
          cdum,cdum,Int_tenst,v0muf,Core_resolved,0._db,Surface_ref,Volume_maille)

  end do

  return
end

!***********************************************************************

subroutine Tensor_pol_dd_cal(ipl,n_tens_dd,nplt,pe,ps,Tensor_pol_dd)

  use declarations
  implicit none
  
  integer:: ipl, n_tens_dd, nplt

  complex(kind=db) Px, Py, Pz, Qx, Qy, Qz
  complex(kind=db), dimension(3):: pe, ps
  complex(kind=db), dimension(n_tens_dd):: Tens
  complex(kind=db), dimension(0:nplt,n_tens_dd):: Tensor_pol_dd

  Px = Pe(1); Qx = conjg( Ps(1) )
  Py = Pe(2); Qy = conjg( Ps(2) )
  Pz = Pe(3); Qz = conjg( Ps(3) )

  Tens(1) = ( Qx*Px + Qy*Py + Qz*Pz )

  Tens(2) = - img * ( Qx*Py - Qy*Px )

  Tens(3) = - img * ( Qy*Pz - Qz*Py )

  Tens(4) = - img * ( Qz*Px - Qx*Pz )

  Tens(5) = ( 2*Qz*Pz - Qx*Px - Qy*Py )

  Tens(6) = ( Qx*Pz + Qz*Px )

  Tens(7) = ( Qy*Pz + Qz*Py )

  Tens(8) = ( Qx*Py + Qy*Px )

  Tens(9) = ( Qx*Px - Qy*Py )

  Tensor_pol_dd(ipl,:) = Tens(:)

  return
end

!***********************************************************************

subroutine Tensor_pol_dq_cal(ipl,n_tens_dq,nplt,pe,ps,ve,vs,Tensor_pol_dq)

  use declarations
  implicit none

  integer:: i, ipl, j, n_tens_dq, nplt
  
  complex(kind=db) Px, Py, Pz, Qx, Qy, Qz
  complex(kind=db), dimension(3):: pe, ps
  complex(kind=db), dimension(2*n_tens_dq):: Tens
  complex(kind=db), dimension(0:nplt,2*n_tens_dq):: Tensor_pol_dq

  real(kind=db):: fac, Vx, Vy, Vz, Wx, Wy, Wz
  
  real(kind=db), dimension(3):: ve, vs

  Px = Pe(1); Qx = conjg( Ps(1) )
  Py = Pe(2); Qy = conjg( Ps(2) )
  Pz = Pe(3); Qz = conjg( Ps(3) )

  Vx = Ve(1); Wx = Vs(1)
  Vy = Ve(2); Wy = Vs(2)
  Vz = Ve(3); Wz = Vs(3)

  Tens(:) = 0._db

  j = 0

  do i = 1,2

    if( i == 2 ) then
      Wx = - Wx; Wy = - Wy; Wz = - Wz
      j = n_tens_dq
    endif

! T(10)
    Tens(j+1) = ( 2*Qz*Pz + 1.5_db*Qx*Px + 1.5_db*Qy*Py ) * ( Vz - Wz ) &
              + ( 1.5_db*Qx*Pz - Qz*Px ) * Vx - ( 1.5_db*Qz*Px - Qx*Pz ) * Wx &
              + ( 1.5_db*Qy*Pz - Qz*Py ) * Vy - ( 1.5_db*Qz*Py - Qy*Pz ) * Wy

! (T(11)
    Tens(j+2) = ( 2*Qx*Px + 1.5_db*Qy*Py + 1.5_db*Qz*Pz ) * ( Vx - Wx ) &
              + ( 1.5_db*Qy*Px - Qx*Py ) * Vy - ( 1.5_db*Qx*Py - Qy*Px ) * Wy &
              + ( 1.5_db*Qz*Px - Qx*Pz ) * Vz - ( 1.5_db*Qx*Pz - Qz*Px ) * Wz

! T(1-1)
    Tens(j+3) = ( 2*Qy*Py + 1.5_db*Qx*Px + 1.5_db*Qz*Pz ) * ( Vy - Wy ) &
              + ( 1.5*Qx*Py - Qy*Px ) * Vx - ( 1.5_db*Qy*Px - Qx*Py ) * Wx &
              + ( 1.5*Qz*Py - Qy*Pz ) * Vz - ( 1.5_db*Qy*Pz - Qz*Py ) * Wz

! T(20)
    Tens(j+4) = 0.5_db * ( ( Qx*Py - Qy*Px ) * ( Vz + Wz ) - Qy*Pz*Vx + Qz*Py*Wx + Qx*Pz*Vy - Qz*Px*Wy )

    fac = 0.5_db

! (T(21)-T(2-1))/sqrt(2)
    Tens(j+5) = fac * ( ( Qz*Pz - Qx*Px ) * ( Vy - Wy ) + ( Qz*Py - 2*Qy*Pz ) * Vz - ( Qy*Pz - 2*Qz*Py ) * Wz &
            + ( 2*Qy*Px - Qx*Py ) * Vx - ( 2*Qx*Py - Qy*Px ) * Wx )

! (T(21)+T(2-1))/sqrt(2)
    Tens(j+6) = fac * ( ( Qz*Pz - Qy*Py ) * ( Vx - Wx ) + ( Qz*Px - 2*Qx*Pz ) * Vz - ( Qx*Pz - 2*Qz*Px ) * Wz &
            + ( 2*Qx*Py - Qy*Px ) * Vy - ( 2*Qy*Px - Qx*Py ) * Wy )

! (T(22)-T(2-2))/sqrt(2)
    Tens(j+7) = fac * ( ( Qx*Px - Qy*Py ) * ( Vz - Wz ) + ( Qx*Pz - 2*Qz*Px ) * Vx - ( Qz*Px - 2*Qx*Pz ) * Wx &
            + ( 2*Qz*Py - Qy*Pz ) * Vy - ( 2*Qy*Pz - Qz*Py ) * Wy )

! (T(22)+T(2-2))/sqrt(2)
    Tens(j+8) = fac * ( ( Qx*Py + Qy*Px ) * ( Vz - Wz ) + ( Qx*Pz - 2*Qz*Px ) * Vy - ( Qz*Px - 2*Qx*Pz ) * Wy &
            - ( 2*Qz*Py - Qy*Pz ) * Vx + ( 2*Qy*Pz - Qz*Py ) * Wx )

! T(30)
    Tens(j+9) = ( 1 / sqrt(10._db) ) * ( ( 2*Qz*Pz - Qx*Px - Qy*Py ) * ( Vz - Wz ) - ( Qx*Pz + Qz*Px ) * ( Vx - Wx ) &
                      - ( Qy*Pz + Qz*Py ) * ( Vy - Wy ) )

    fac = 1 / sqrt( 60._db )

! (T(31)-T(3-1))/sqrt(2)
    Tens(j+10) = fac * ( ( 3*Qx*Px + Qy*Py - 4*Qz*Pz ) * ( Vx - Wx ) + ( Qx*Py + Qy*Px ) * ( Vy - Wy ) &
                    - 4 * ( Qx*Pz + Qz*Px ) * ( Vz - Wz ) )

! (T(31)+T(3-1))/sqrt(2)
    Tens(j+11) = fac * ( ( 3*Qy*Py + Qx*Px - 4*Qz*Pz ) * ( Vy - Wy ) + ( Qx*Py + Qy*Px ) * ( Vx - Wx ) &
                    - 4 * ( Qy*Pz + Qz*Py ) * ( Vz - Wz ) )

    fac = 1 / sqrt( 6._db )

! (T(32)-T(3-2))/sqrt(2)
    Tens(j+12) = fac * ( ( Qx*Py + Qy*Px ) * ( Vz - Wz ) + ( Qx*Pz + Qz*Px ) * ( Vy - Wy ) &
                             + ( Qy*Pz + Qz*Py ) * ( Vx - Wx ) )

! (T(32)+T(3-2))/sqrt(2)
    Tens(j+13) = fac  * ( ( Qx*Px - Qy*Py ) * ( Vz - Wz ) + ( Qx*Pz + Qz*Px ) * ( Vx - Wx ) &
                        - ( Qy*Pz + Qz*Py ) * ( Vy - Wy ) )

! (T(33)-T(3-3))/sqrt(2)
    Tens(j+14) = 0.5_db * ( ( Qy*Py - Qx*Px ) * ( Vx - Wx ) + ( Qy*Px + Qx*Py ) * ( Vy - Wy ) )

! (T(33)+T(3-3))/sqrt(2)
    Tens(j+15) = 0.5_db * ( ( Qy*Py - Qx*Px ) * ( Vy - Wy ) - ( Qy*Px + Qx*Py ) * ( Vx - Wx ) )

  end do

  Tensor_pol_dq(ipl,:) = Tens(:)

  return
end

!***********************************************************************

subroutine Tensor_pol_qq_cal(ipl,n_tens_qq,nplt,pe,ps,ve,vs, Tensor_pol_qq)

  use declarations
  implicit none
  
  integer:: ipl, n_tens_qq, nplt

  complex(kind=db) Px, Py, Pz, Qx, Qy, Qz
  complex(kind=db), dimension(3):: pe, ps
  complex(kind=db), dimension(n_tens_qq):: Tens
  complex(kind=db), dimension(0:nplt,n_tens_qq):: Tensor_pol_qq

  real(kind=db):: fac, Vx, Vy, Vz, Wx, Wy, Wz
  
  real(kind=db), dimension(3):: ve, vs

  Px = Pe(1); Qx = conjg( Ps(1) )
  Py = Pe(2); Qy = conjg( Ps(2) )
  Pz = Pe(3); Qz = conjg( Ps(3) )

  Vx = Ve(1); Wx = Vs(1)
  Vy = Ve(2); Wy = Vs(2)
  Vz = Ve(3); Wz = Vs(3)

! Scalaire

  Tens(1) = 1.5_db * ( Qx*Wz*Px*Vz + Qx*Wz*Pz*Vx + Qz*Wx*Px*Vz + Qz*Wx*Pz*Vx + Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy &
                     + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy + Qx*Wy*Px*Vy + Qx*Wy*Py*Vx + Qy*Wx*Px*Vy + Qy*Wx*Py*Vx ) &
          + 2 * ( Qx*Wx*Px*Vx + Qy*Wy*Py*Vy + Qz*Wz*Pz*Vz ) &
          - ( Qx*Wx*Py*Vy + Qy*Wy*Px*Vx + Qx*Wx*Pz*Vz + Qz*Wz*Px*Vx + Qy*Wy*Pz*Vz + Qz*Wz*Py*Vy )

! Dipole magnetique

  Tens(2) = img * ( Qy*Wy*Px*Vy + Qy*Wy*Py*Vx + Qx*Wy*Px*Vx + Qy*Wx*Px*Vx &
                  - Qx*Wy*Py*Vy - Qy*Wx*Py*Vy - Qx*Wx*Px*Vy - Qx*Wx*Py*Vx ) &
          + 0.5_db * img * ( Qy*Wz*Px*Vz + Qy*Wz*Pz*Vx + Qz*Wy*Px*Vz + Qz*Wy*Pz*Vx &
                           - Qx*Wz*Py*Vz - Qx*Wz*Pz*Vy - Qz*Wx*Py*Vz - Qz*Wx*Pz*Vy )

  Tens(3) = img * ( Qz*Wz*Py*Vz + Qz*Wz*Pz*Vy + Qy*Wz*Py*Vy + Qz*Wy*Py*Vy &
                  - Qy*Wz*Pz*Vz - Qz*Wy*Pz*Vz - Qy*Wy*Py*Vz - Qy*Wy*Pz*Vy ) &
          + 0.5_db * img * ( Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx &
                           - Qx*Wy*Px*Vz - Qx*Wy*Pz*Vx - Qy*Wx*Px*Vz - Qy*Wx*Pz*Vx )

  Tens(4) = - img * ( Qz*Wz*Px*Vz + Qz*Wz*Pz*Vx + Qx*Wz*Px*Vx + Qz*Wx*Px*Vx &
                    - Qx*Wx*Px*Vz - Qx*Wx*Pz*Vx - Qx*Wz*Pz*Vz - Qz*Wx*Pz*Vz ) &
          + 0.5_db * img * ( Qx*Wy*Py*Vz + Qx*Wy*Pz*Vy + Qy*Wx*Py*Vz + Qy*Wx*Pz*Vy &
                           - Qy*Wz*Px*Vy - Qy*Wz*Py*Vx - Qz*Wy*Px*Vy - Qz*Wy*Py*Vx )

! Quadrupole non-magnetique

  fac = 1 / ( 3 * sqrt(14._db) )

  Tens(5) = 3 * fac * ( Qx*Wy*Px*Vy + Qx*Wy*Py*Vx + Qy*Wx*Px*Vy + Qy*Wx*Py*Vx ) - 1.5_db * fac &
     * ( Qx*Wz*Px*Vz + Qx*Wz*Pz*Vx + Qz*Wx*Px*Vz + Qz*Wx*Pz*Vx + Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy ) &
          + 2 * fac * ( Qx*Wx*Px*Vx + Qy*Wy*Py*Vy + Qx*Wx*Pz*Vz + Qz*Wz*Px*Vx + Qy*Wy*Pz*Vz + Qz*Wz*Py*Vy) - 4 * fac &
     * ( Qx*Wx*Py*Vy + Qy*Wy*Px*Vx + Qz*Wz*Pz*Vz )

  fac = 1 / sqrt( 42._db )

  Tens(6) = 1.5_db * fac * ( Qy*Wz*Px*Vy + Qz*Wy*Px*Vy + Qy*Wz*Py*Vx + Qz*Wy*Py*Vx &
       + Qx*Wy*Py*Vz + Qx*Wy*Pz*Vy + Qy*Wx*Py*Vz + Qy*Wx*Pz*Vy ) + fac &
     * ( Qz*Wz*Px*Vz + Qz*Wz*Pz*Vx + Qx*Wx*Px*Vz + Qx*Wx*Pz*Vx + Qx*Wz*Pz*Vz + Qz*Wx*Pz*Vz + Qx*Wz*Px*Vx + Qz*Wx*Px*Vx ) &
          - 2 * fac * ( Qy*Wy*Px*Vz + Qy*Wy*Pz*Vx + Qx*Wz*Py*Vy + Qz*Wx*Py*Vy )

  Tens(7) = 1.5_db * fac * ( Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx &
       + Qx*Wy*Px*Vz + Qx*Wy*Pz*Vx + Qy*Wx*Px*Vz + Qy*Wx*Pz*Vx ) + fac  &
     * ( Qz*Wz*Py*Vz + Qz*Wz*Pz*Vy + Qy*Wy*Py*Vz + Qy*Wy*Pz*Vy + Qy*Wz*Pz*Vz + Qz*Wy*Pz*Vz + Qy*Wz*Py*Vy + Qz*Wy*Py*Vy ) &
          - 2 * fac * ( Qx*Wx*Py*Vz + Qx*Wx*Pz*Vy + Qy*Wz*Px*Vx + Qz*Wy*Px*Vx )

  fac = 1 / sqrt( 42._db )

  Tens(8) = - fac * ( Qx*Wx*Px*Vy + Qx*Wx*Py*Vx + Qy*Wy*Px*Vy + Qy*Wy*Py*Vx &
       + Qx*Wy*Px*Vx + Qy*Wx*Px*Vx + Qx*Wy*Py*Vy + Qy*Wx*Py*Vy ) + 2 * fac &
     * ( Qz*Wz*Px*Vy + Qz*Wz*Py*Vx + Qx*Wy*Pz*Vz + Qy*Wx*Pz*Vz ) - 1.5_db * fac &
     * ( Qx*Wz*Py*Vz + Qx*Wz*Pz*Vy + Qz*Wx*Py*Vz + Qz*Wx*Pz*Vy + Qy*Wz*Px*Vz + Qy*Wz*Pz*Vx + Qz*Wy*Px*Vz + Qz*Wy*Pz*Vx )

  Tens(9) = 2 * fac * ( Qz*Wz*Px*Vx - Qz*Wz*Py*Vy + Qx*Wx*Pz*Vz - Qy*Wy*Pz*Vz + Qy*Wy*Py*Vy - Qx*Wx*Px*Vx ) + 1.5_db * fac &
     * ( Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy - Qx*Wz*Px*Vz - Qx*Wz*Pz*Vx - Qz*Wx*Px*Vz - Qz*Wx*Pz*Vx )

! Octupole magnetique

  fac = 1 / sqrt( 10._db )

  Tens(10) = 0.5_db * fac * img * ( Qx*Wy*Px*Vx - Qx*Wy*Py*Vy + Qy*Wx*Px*Vx - Qy*Wx*Py*Vy &
       - Qx*Wx*Px*Vy + Qy*Wy*Px*Vy - Qx*Wx*Py*Vx + Qy*Wy*Py*Vx ) + fac * img &
     * ( Qx*Wz*Py*Vz + Qx*Wz*Pz*Vy + Qz*Wx*Py*Vz + Qz*Wx*Pz*Vy - Qy*Wz*Px*Vz - Qz*Wy*Px*Vz - Qy*Wz*Pz*Vx - Qz*Wy*Pz*Vx )

  fac = 1 / ( 4 * sqrt( 15._db ) )

  Tens(11) = 3 * fac * img * ( Qx*Wy*Px*Vz + Qx*Wy*Pz*Vx + Qy*Wx*Px*Vz + Qy*Wx*Pz*Vx &
       - Qx*Wz*Px*Vy - Qx*Wz*Py*Vx - Qz*Wx*Px*Vy - Qz*Wx*Py*Vx ) + 4 * fac * img &
     * ( Qz*Wz*Py*Vz + Qz*Wz*Pz*Vy - Qy*Wz*Pz*Vz - Qz*Wy*Pz*Vz ) + 5 * fac * img &
     * ( Qy*Wz*Px*Vx + Qz*Wy*Px*Vx - Qx*Wx*Py*Vz - Qx*Wx*Pz*Vy ) + fac * img &
     * ( Qy*Wy*Py*Vz + Qy*Wy*Pz*Vy - Qy*Wz*Py*Vy - Qz*Wy*Py*Vy )

  Tens(12) = 3 * fac * img * ( Qx*Wy*Py*Vz + Qx*Wy*Pz*Vy + Qy*Wx*Py*Vz + Qy*Wx*Pz*Vy &
       - Qy*Wz*Px*Vy - Qy*Wz*Py*Vx - Qz*Wy*Px*Vy - Qz*Wy*Py*Vx ) + 4 * fac * img &
     * ( Qz*Wz*Px*Vz + Qz*Wz*Pz*Vx - Qx*Wz*Pz*Vz - Qz*Wx*Pz*Vz ) + 5 * fac * img &
     * ( Qx*Wz*Py*Vy + Qz*Wx*Py*Vy - Qy*Wy*Px*Vz - Qy*Wy*Pz*Vx ) + fac * img &
     * ( Qx*Wx*Px*Vz + Qx*Wx*Pz*Vx - Qx*Wz*Px*Vx - Qz*Wx*Px*Vx )

  fac = 1 / sqrt( 6._db )

  Tens(13) = fac * img * ( Qx*Wx*Pz*Vz - Qz*Wz*Px*Vx + Qy*Wy*Px*Vx - Qx*Wx*Py*Vy + Qz*Wz*Py*Vy - Qy*Wy*Pz*Vz )

  Tens(14) = fac * img * ( Qx*Wy*Pz*Vz + Qy*Wx*Pz*Vz - Qz*Wz*Px*Vy - Qz*Wz*Py*Vx ) + 0.5 * fac * img &
     * ( Qx*Wx*Px*Vy + Qx*Wx*Py*Vx + Qy*Wy*Px*Vy + Qy*Wy*Py*Vx - Qx*Wy*Px*Vx - Qy*Wx*Px*Vx - Qx*Wy*Py*Vy - Qy*Wx*Py*Vy )

  fac = 0.25_db

  Tens(15) = fac * img * ( Qy*Wz*Px*Vx - Qy*Wz*Py*Vy + Qz*Wy*Px*Vx - Qz*Wy*Py*Vy &
       + Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx - Qx*Wx*Py*Vz + Qy*Wy*Py*Vz - Qx*Wx*Pz*Vy + Qy*Wy*Pz*Vy &
       - Qx*Wy*Px*Vz - Qy*Wx*Px*Vz - Qx*Wy*Pz*Vx - Qy*Wx*Pz*Vx )

  Tens(16) = fac * img * ( Qx*Wz*Px*Vx - Qx*Wz*Py*Vy + Qz*Wx*Px*Vx - Qz*Wx*Py*Vy &
       - Qy*Wz*Px*Vy - Qy*Wz*Py*Vx - Qz*Wy*Px*Vy - Qz*Wy*Py*Vx - Qx*Wx*Px*Vz + Qy*Wy*Px*Vz - Qx*Wx*Pz*Vx + Qy*Wy*Pz*Vx &
       + Qx*Wy*Py*Vz + Qy*Wx*Py*Vz + Qx*Wy*Pz*Vy + Qy*Wx*Pz*Vy )

! Hexadecapole non-magnetique

  fac = 1 / ( 2 * sqrt( 70._db ) )

  Tens(17) = 3 * fac * ( Qx*Wx*Px*Vx + Qy*Wy*Py*Vy ) + fac * ( Qx*Wx*Py*Vy + Qy*Wy*Px*Vx + Qx*Wy*Px*Vy + Qx*Wy*Py*Vx &
       + Qy*Wx*Px*Vy + Qy*Wx*Py*Vx ) + 8 * fac * Qz*Wz*Pz*Vz - 4 * fac &
     * ( Qz*Wz*Px*Vx + Qz*Wz*Py*Vy + Qx*Wx*Pz*Vz + Qy*Wy*Pz*Vz + Qx*Wz*Px*Vz + Qx*Wz*Pz*Vx + Qz*Wx*Px*Vz + Qz*Wx*Pz*Vx &
       + Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy )

  fac = 1 / ( 4 * sqrt( 7._db ) )

  Tens(18) = 3 * fac * ( Qx*Wz*Px*Vx + Qz*Wx*Px*Vx + Qx*Wx*Px*Vz + Qx*Wx*Pz*Vx ) + fac &
     * ( Qx*Wz*Py*Vy + Qz*Wx*Py*Vy + Qy*Wy*Px*Vz + Qy*Wy*Pz*Vx + Qy*Wz*Px*Vy + Qy*Wz*Py*Vx + Qz*Wy*Px*Vy + Qz*Wy*Py*Vx &
       + Qy*Wx*Pz*Vy + Qx*Wy*Pz*Vy + Qx*Wy*Py*Vz + Qy*Wx*Py*Vz ) - 4 * fac &
     * ( Qz*Wz*Px*Vz + Qz*Wz*Pz*Vx + Qx*Wz*Pz*Vz + Qz*Wx*Pz*Vz )

  Tens(19) = 3 * fac * ( Qy*Wz*Py*Vy + Qz*Wy*Py*Vy + Qy*Wy*Py*Vz + Qy*Wy*Pz*Vy ) + fac &
     * ( Qy*Wz*Px*Vx + Qz*Wy*Px*Vx + Qx*Wx*Py*Vz + Qx*Wx*Pz*Vy + Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx &
       + Qx*Wy*Px*Vz + Qx*Wy*Pz*Vx + Qy*Wx*Px*Vz + Qy*Wx*Pz*Vx ) - 4 * fac &
     * ( Qy*Wz*Pz*Vz + Qz*Wy*Pz*Vz + Qz*Wz*Pz*Vy + Qz*Wz*Py*Vz )

  fac = 1 / sqrt( 14._db )

  Tens(20) = fac * ( Qz*Wz*Px*Vy + Qz*Wz*Py*Vx + Qx*Wy*Pz*Vz + Qy*Wx*Pz*Vz &
       + Qx*Wz*Py*Vz + Qx*Wz*Pz*Vy + Qz*Wx*Py*Vz + Qz*Wx*Pz*Vy + Qy*Wz*Px*Vz + Qy*Wz*Pz*Vx + Qz*Wy*Px*Vz + Qz*Wy*Pz*Vx ) &
           - 0.5_db * fac * ( Qx*Wx*Px*Vy + Qx*Wx*Py*Vx + Qy*Wy*Px*Vy + Qy*Wy*Py*Vx &
       + Qx*Wy*Px*Vx + Qy*Wx*Px*Vx + Qx*Wy*Py*Vy + Qy*Wx*Py*Vy )

  Tens(21) = fac * ( Qz*Wz*Px*Vx - Qz*Wz*Py*Vy - Qx*Wx*Px*Vx + Qy*Wy*Py*Vy &
       + Qx*Wx*Pz*Vz - Qy*Wy*Pz*Vz + Qx*Wz*Px*Vz + Qx*Wz*Pz*Vx + Qz*Wx*Px*Vz + Qz*Wx*Pz*Vx - Qy*Wz*Py*Vz - Qy*Wz*Pz*Vy &
       - Qz*Wy*Py*Vz - Qz*Wy*Pz*Vy )

  fac = 0.25_db

  Tens(22) = - fac * ( Qx*Wz*Px*Vx - Qx*Wz*Py*Vy + Qz*Wx*Px*Vx - Qz*Wx*Py*Vy &
       - Qy*Wz*Py*Vx - Qy*Wz*Px*Vy - Qz*Wy*Py*Vx - Qz*Wy*Px*Vy - Qx*Wy*Pz*Vy + Qx*Wx*Px*Vz - Qy*Wy*Px*Vz + Qx*Wx*Pz*Vx &
       - Qy*Wy*Pz*Vx - Qy*Wx*Py*Vz - Qx*Wy*Py*Vz - Qy*Wx*Pz*Vy )

  Tens(23) = - fac * ( Qy*Wz*Px*Vx - Qy*Wz*Py*Vy + Qz*Wy*Px*Vx - Qz*Wy*Py*Vy &
       + Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx + Qy*Wx*Pz*Vx + Qx*Wx*Py*Vz - Qy*Wy*Py*Vz + Qx*Wx*Pz*Vy &
       - Qy*Wy*Pz*Vy + Qx*Wy*Px*Vz + Qy*Wx*Px*Vz + Qx*Wy*Pz*Vx )

  fac = 1 / ( 2 * sqrt( 2._db ) )

  Tens(24) = fac * ( Qx*Wy*Px*Vx - Qx*Wy*Py*Vy + Qy*Wx*Px*Vx - Qy*Wx*Py*Vy &
       + Qx*Wx*Px*Vy + Qx*Wx*Py*Vx - Qy*Wy*Px*Vy - Qy*Wy*Py*Vx )

  Tens(25) = fac * ( Qx*Wx*Px*Vx + Qy*Wy*Py*Vy - Qx*Wx*Py*Vy - Qy*Wy*Px*Vx &
       - Qx*Wy*Px*Vy - Qx*Wy*Py*Vx - Qy*Wx*Px*Vy - Qy*Wx*Py*Vx )

  Tensor_pol_qq(ipl,:) = Tens(:)

  return
end

!***********************************************************************

subroutine Tensor_pol_dm_cal(ipl,n_tens_dm,nplt,Pl_i,Pl_s,Vo_i,Vo_s,Tensor_pol_dm)

  use declarations
  implicit none
  
  integer:: i, ipl, is, j, k, n_tens_dm, nplt 

  complex(kind=db), dimension(3):: Pl_i, Pl_s, U_i, U_s
  complex(kind=db), dimension(2*n_tens_dm):: Tens
  complex(kind=db), dimension(3,3):: T
  complex(kind=db), dimension(0:nplt,2*n_tens_dm):: Tensor_pol_dm

  real(kind=db), dimension(3):: Vo_i, Vo_s
  
  Tens(:) = ( 0._db, 0._db )

  U_i(1) = Vo_i(2) * Pl_i(3) - Vo_i(3) * Pl_i(2)
  U_i(2) = Vo_i(3) * Pl_i(1) - Vo_i(1) * Pl_i(3)
  U_i(3) = Vo_i(1) * Pl_i(2) - Vo_i(2) * Pl_i(1)
  U_s(1) = Vo_s(2) * Pl_s(3) - Vo_s(3) * Pl_s(2)
  U_s(2) = Vo_s(3) * Pl_s(1) - Vo_s(1) * Pl_s(3)
  U_s(3) = Vo_s(1) * Pl_s(2) - Vo_s(2) * Pl_s(1)

! is = -1: non magnetique; is = 1 : magnetique
  do is = -1,1,2

    do i = 1,3
      do j = 1,3
        T(i,j) = Conjg( Pl_s(i) ) * U_i(j) + is * Pl_i(i) * Conjg( U_s(j) ) 
      end do
    end do

    if( is == - 1 ) then
      k = 0
    else
      k = n_tens_dm
    endif

    Tens(k+1) = ( T(1,1) + T(2,2) + T(3,3) )

    Tens(k+2) = ( T(1,2) - T(2,1) )

    Tens(k+3) = ( T(2,3) - T(3,2) )

    Tens(k+4) = ( T(3,1) - T(1,3) )

    Tens(k+5) = ( 2*T(3,3) - T(1,1) - T(2,2) )

    Tens(k+6) = ( T(1,3) + T(3,1) )

    Tens(k+7) = ( T(2,3) + T(3,2) )

    Tens(k+8) = ( T(1,2) + T(2,1) )

    Tens(k+9) = ( T(1,1) - T(2,2) )
  
  end do

  Tensor_pol_dm(ipl,:) = Tens(:)

  return
end

!***********************************************************************

! Ecriture du signal par tenseur

subroutine Write_Signal(Abs_U_iso,ct_nelec,Core_resolved,Volume_maille,E_cut,E1E1,E1E2,E1M1,E2E2,Ephseuil,Epsii, &
      Eseuil,First_E,i_range,ia,ipl,ipldafs,jseuil,magn_sens,n_tens_dd,n_tens_dq,n_tens_dm, &
      n_tens_qq,ninit1,ninitlr,nomfich_s,npldafs,nplt,nseuil,numat_abs,phdf0t,phdt, &
      Sph_tensor_dd_ni,Sph_tensor_dq_ni,Sph_tensor_dq_m_ni,Sph_tensor_dm_ni,Sph_tensor_dm_m_ni,Sph_tensor_qq_ni, &
      Surface_ref,Tensor_pol_dd,Tensor_pol_dq,Tensor_pol_dm,Tensor_pol_qq,V0muf)

  use declarations
  implicit none
  
  integer:: i, i_range, ia, initlr, ipl, ipldafs, is, j, j0, jdd, jdd0, jdm, jdmm, jdm0, jdmm0, &
      jdq, jdq0, jdqm, jdqm0, jqq, jqq0, jseuil, jtot, jtot0, long, n_tens, n_tens_dd, n_tens_dq, n_tens_dm, &
      n_tens_qq, n_tens2, ninit1, ninitlr, npldafs, nplt, nseuil, numat_abs

  character(len=Length_name):: nomfich_s, nomficht
  character(len=Length_word), dimension(:), allocatable:: Tens_name

  complex(kind=db):: cf, cg, zero_c
  complex(kind=db), dimension(1):: cdum
  complex(kind=db), dimension(n_tens_dd):: Sph_tensor_dd
  complex(kind=db), dimension(n_tens_dq):: Sph_tensor_dq, Sph_tensor_dq_m
  complex(kind=db), dimension(n_tens_dm):: Sph_tensor_dm, Sph_tensor_dm_m
  complex(kind=db), dimension(n_tens_qq):: Sph_tensor_qq
  complex(kind=db), dimension(n_tens_dd,ninitlr):: Sph_tensor_dd_ni
  complex(kind=db), dimension(n_tens_dq,ninitlr):: Sph_tensor_dq_ni, Sph_tensor_dq_m_ni
  complex(kind=db), dimension(n_tens_dm,ninitlr):: Sph_tensor_dm_ni, Sph_tensor_dm_m_ni
  complex(kind=db), dimension(n_tens_qq,ninitlr):: Sph_tensor_qq_ni
  complex(kind=db), dimension(0:nplt,n_tens_dd):: Tensor_pol_dd
  complex(kind=db), dimension(0:nplt,2*n_tens_dq):: Tensor_pol_dq
  complex(kind=db), dimension(0:nplt,2*n_tens_dm):: Tensor_pol_dm
  complex(kind=db), dimension(0:nplt,n_tens_qq):: Tensor_pol_qq
  complex(kind=db), dimension(npldafs):: phdf0t, phdt
  complex(kind=db), dimension(:), allocatable:: ph0, phtem, Resul

  logical:: Core_resolved, Dafs, E1E1, E1E2, E1M1, E2E2, First_E, Magn_sens
  
  real(kind=db):: Abs_U_iso, E_cut, Ephseuil, Eseuil, fac, Surface_ref, V0muf, Volume_maille
  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(ninitlr):: ct_nelec, Epsii
  real(kind=db), dimension(:), allocatable:: Tens

  zero_c = (0._db, 0._db)
  Dafs = ipldafs > 0

! Nombre de colonnes sommes
  j0 = 0
  if( E1E1 ) j0 = j0 + 1
  if( E1E2 ) j0 = j0 + 1
  if( E1E2 .and. Magn_sens ) j0 = j0 + 1
  if( E2E2 ) j0 = j0 + 1
  if( E1M1 ) j0 = j0 + 1
  if( E1M1 .and. Magn_sens ) j0 = j0 + 1
  if( j0 > 1 ) j0 = j0 + 1 

  n_tens = j0
  
! Indices des colones sommes
  jtot0 = 1; jdd0 = 0; jdq0 = 0; jdqm0 = 0; jqq0 = 0; jdm0 = 0; jdmm0 = 0
  if( j0 == 1 ) then
    j = 0
  else
    j = 1
  endif
  if( E1E1 ) then
    j = j + 1
    jdd0 = j
    n_tens = n_tens + n_tens_dd
  endif
  if( E1E2 ) then
    j = j + 1
    jdq0 = j
    n_tens = n_tens + n_tens_dq
  endif
  if( E1E2 .and. Magn_sens) then
    j = j + 1
    jdqm0 = j
    n_tens = n_tens + n_tens_dq
  endif
  if( E2E2 ) then
    j = j + 1
    jqq0 = j
    n_tens = n_tens + n_tens_qq
  endif
  if( E1M1 ) then
    j = j + 1
    jdm0 = j
    n_tens = n_tens + n_tens_dm
  endif
  if( E1M1 .and. Magn_sens) then
    j = j + 1
    jdmm0 = j
    n_tens = n_tens + n_tens_dm
  endif

  n_tens = n_tens * ninitlr
  allocate( Resul( n_tens ) )
  allocate( Ph0( n_tens ) )
  allocate( phtem( n_tens ) )
  
  Resul(:) = (0._db,0._db)
  
  if( Dafs ) n_tens = 2 * n_tens
  
  allocate( Tens( n_tens ) )
  allocate( Tens_name( n_tens ) )
  
  j = 0
  
  do initlr = 1,ninitlr
      
    jtot = jtot0 + j
    jdd = jdd0 + j
    jdq = jdq0 + j
    jdqm = jdqm0 + j
    jqq = jqq0 + j
    jdm = jdm0 + j
    jdmm = jdmm0 + j
    j = j + j0
  
! On donne le Xanes en Megabarns
! Pour le DAFS, dans la convolution on multiplie par pi
    if( Dafs ) then
      fac = 1 / pi
    else
      fac = 1 / ct_nelec(initlr) 
    endif      

    Sph_tensor_dd(:) = fac * Sph_tensor_dd_ni(:,initlr) 
    Sph_tensor_dq(:) = fac * Sph_tensor_dq_ni(:,initlr) 
    Sph_tensor_dq_m(:) = fac * Sph_tensor_dq_m_ni(:,initlr) 
    Sph_tensor_qq(:) = fac * Sph_tensor_qq_ni(:,initlr)
    Sph_tensor_dm(:) = fac * Sph_tensor_dm_ni(:,initlr) 
    Sph_tensor_dm_m(:) = fac * Sph_tensor_dm_m_ni(:,initlr) 

    if( E1E1 ) then
      
      do i = 1,n_tens_dd

        j = j + 1
        Resul(j) = Tensor_pol_dd(ipl,i) * Sph_tensor_dd(i)

! Devant le produit, il faut mettre un signe -1 devant les tenseurs impairs (vient des 2 multiplications par -i)
        if( i >= 2 .and. i <= 4 ) Resul(j) = - Resul(j)

        call Fill_line( Core_resolved, 1, Dafs, i, initlr, j, jseuil, n_tens, nseuil, ninitlr, Resul(j), Tens, Tens_name )

        Resul(jdd) = Resul(jdd) + Resul(j)

      end do
 
      call Fill_line( Core_resolved, 1, Dafs, 0, initlr, jdd, jseuil, n_tens, nseuil, ninitlr, Resul(jdd), Tens, Tens_name )
      
    endif

    if( E1E2 ) then

      do is = 1,2 
        if( is == 2 .and. .not. Magn_sens ) exit
        do i = 1,n_tens_dq
        
          j = j + 1
! Les parties reelles et imaginaires sont a considerer avant la multiplication par le img eventuel.
          if( is == 1 ) then
            Resul(j) = Sph_tensor_dq(i) * Tensor_pol_dq(ipl,i)
          else
            Resul(j) = Sph_tensor_dq_m(i) * Tensor_pol_dq(ipl,i+n_tens_dq)
          endif

! On recupere le img exterieur au tenseur propre au dipole-quadrupole
          Resul(j) = img * Resul(j)

! Signe est oppose car complex conjugue fait dans conv
          if( Dafs ) Resul(j) = - Resul(j) 

          call Fill_line( Core_resolved, 1+is, Dafs, i, initlr, j, jseuil, n_tens, nseuil, ninitlr, Resul(j), Tens, Tens_name )

          if( is == 1 ) then
             Resul(jdq) = Resul(jdq) + Resul(j)
          else
             Resul(jdqm) = Resul(jdqm) + Resul(j)
          endif

        end do

        if( is == 1 ) then
          call Fill_line( Core_resolved, 2, Dafs, 0, initlr, jdq, jseuil, n_tens, nseuil, ninitlr, Resul(jdq), Tens, Tens_name )
        else
          call Fill_line( Core_resolved, 3, Dafs, 0, initlr, jdqm, jseuil, n_tens, nseuil, ninitlr, Resul(jdqm), Tens, Tens_name )
        endif
      
      end do
      
    endif
        
    if( E2E2 ) then

      do i = 1,n_tens_qq
      
        j = j + 1

        Resul(j) = Tensor_pol_qq(ipl,i) * Sph_tensor_qq(i)

        call Fill_line( Core_resolved, 4, Dafs, i, initlr, j, jseuil, n_tens, nseuil, ninitlr, Resul(j), Tens, Tens_name )

        Resul(jqq) = Resul(jqq) + Resul(j)

      end do

      call Fill_line( Core_resolved, 4, Dafs, 0, initlr, jqq, jseuil, n_tens, nseuil, ninitlr, Resul(jqq), Tens, Tens_name )
        
    endif
      
    if( E1M1 ) then

      do is = 1,2
        if( is == 2 .and. .not. Magn_sens ) exit
     
        do i = 1,n_tens_dm
        
          j = j + 1
! Modif pour que ca cole (je ne comprend pas vraiment), voir aussi fin de la routine Sph_tensor_dm_cal 
          if( is == 1 ) then
            Resul(j) = Sph_tensor_dm(i) * Tensor_pol_dm(ipl,i)
            if( Dafs .and. ( i >= 2 .and. i <= 4 )  ) Resul(j) = - Resul(j) 
          else
            Resul(j) = Sph_tensor_dm_m(i) * Tensor_pol_dm(ipl,i+n_tens_dm)
            if( i >= 2 .and. i <= 4 ) Resul(j) = - Resul(j) 
          endif

          call Fill_line( Core_resolved, 4+is, Dafs, i, initlr, j, jseuil, n_tens, nseuil, ninitlr, Resul(j), Tens, Tens_name )

          if( is == 1 ) then
             Resul(jdm) = Resul(jdm) + Resul(j)
          else
             Resul(jdmm) = Resul(jdmm) + Resul(j)
          endif

        end do

        if( is == 1 ) then
          call Fill_line( Core_resolved, 5, Dafs, 0, initlr, jdm, jseuil, n_tens, nseuil, ninitlr, Resul(jdm), Tens, Tens_name )
        else
          call Fill_line( Core_resolved, 6, Dafs, 0, initlr, jdmm, jseuil, n_tens, nseuil, ninitlr, Resul(jdmm), Tens, Tens_name )
        endif
      
      end do

    endif

    if( j0 > 1 ) then
      Resul(jtot) = sum( Resul(jtot+1:jtot+j0-1) )  
      call Fill_line( Core_resolved, 0, Dafs, 0, initlr, jtot, jseuil, n_tens, nseuil, ninitlr, Resul(jtot), Tens, Tens_name )
    endif    

  end do
  
! Ecriture des tenseurs

  nomficht = nomfich_s

  long = len_trim(nomficht)
  nomficht(long+1:long+11) = '_sph_signal'
  if( ia > 0 ) then
    nomficht(long+12:long+16) = '_atom'
    call ad_number(ia,nomficht,Length_name)
  endif
  long = len_trim(nomficht)
  if( ipldafs > 0 ) then
    nomficht(long+1:long+4) = '_rxs'
    call ad_number(ipldafs,nomficht,Length_name)
  elseif( ipl == 0 ) then
    nomficht(long+1:long+4) = '_xan'
  else
    nomficht(long+1:long+4) = '_pol'
    call ad_number(ipl,nomficht,Length_name)
  endif
  long = len_trim(nomficht)
  nomficht(long+1:long+4) = '.txt'

  if( Dafs ) then
    if( ia == 0 ) then
      cf = phdt(ipldafs)
      cg = phdf0t(ipldafs)
    else
      cf = (1._db,0._db)
      cg = (0._db,0._db)
    endif
    phtem(:) = cf
    ph0(:) = cg
    n_tens2 = n_tens / ( 2 * ninitlr )
    call write_out(Abs_U_iso,rdum,rdum,zero_c,E_cut,Ephseuil,Epsii,Eseuil,First_E,.false.,rdum, &
          i_range,jseuil,n_tens,n_tens,ninit1,ninitlr,nomficht,Tens_name,n_tens,n_tens2,0,0,nseuil,numat_abs, &
          phtem,ph0,Tens,v0muf,Core_resolved,0._db,Surface_ref,Volume_maille)
  else
    call write_out(Abs_U_iso,rdum,rdum,zero_c,E_cut,Ephseuil,Epsii,Eseuil,First_E,.false.,rdum, &
          i_range,jseuil,n_tens,n_tens,ninit1,ninitlr,nomficht,Tens_name,1,0,0,0,nseuil,numat_abs, &
          cdum,cdum,Tens,V0muf,Core_resolved,0._db,Surface_ref,Volume_maille)
  endif

  deallocate( Ph0, Phtem, Resul, Tens, Tens_name )

  return
end

!***************************************************************************************************************

! Remplissage valeur et nom du tenseur

subroutine Fill_line( Core_resolved, index, Dafs, i, initlr, j, jseuil, n_tens, nseuil, ninitlr, Resul, Tens, Tens_name )

  use declarations
  implicit none
  
  integer:: i, index, initlr, j, jj, jseuil, l, n_tens, ninitlr, nseuil
  
  character(len=Length_word):: mot
  character(len=8), dimension(0:0):: Tens_name_tot
  character(len=8), dimension(0:9):: Tens_name_D
  character(len=8), dimension(0:15):: Tens_name_I, Tens_name_I_m 
  character(len=8), dimension(0:25):: Tens_name_Q 
  character(len=8), dimension(0:9):: Tens_name_T, Tens_name_T_m 
  character(len=Length_word), dimension(n_tens):: Tens_name

  complex(kind=db):: Resul
  
  logical:: Core_resolved, Dafs 

  real(kind=db), dimension(n_tens):: Tens
  
  data Tens_name_tot/ ' Sum_tot'/
  data Tens_name_D/ ' Sum_dd ','  D(00) ','  lz_dd ','  lx_dd ','  ly_dd ','  D_z2  ','   D_xz ','  D_yz  ','  D_xy  ',' D_x2-y2'/
  data Tens_name_I/ ' Sum_dq ','   I_z  ','   I_x  ','    I_y ','  I(20) ',' I(21)  ',' I(2-1) ','  I(22) ',' I(2-2) ', &
                    '  I(30) ',' I(31)  ',' I(3-1) ',' I(32)  ',' I(3-2) ','  I(33) ',' I(3-3) '/
  data Tens_name_I_m/ 'Sum_dq_m','  I_z_m ',' I_x_m  ','  I_y_m ',' I(20)_m',' I(21)_m','I(2-1)_m',' I(22)_m','I(2-2)_m', &
                      ' I(30)_m',' I(31)_m','I(3-1)_m',' I(32)_m','I(3-2)_m',' I(33)_m','I(3-3)_m'/
  data Tens_name_Q/ ' Sum_qq ','  Q(00) ','  lz_qq ','  lx_qq ','  ly_qq ','  Q(20) ','  Q(21) ',' Q(2-1) ','  Q(22) ',' Q(2-2) ', &
                    '  Q(30) ','  Q(31) ',' Q(3-1) ','  Q(32) ',' Q(3-2) ','  Q(33) ',' Q(3-3) ','  Q(40) ','  Q(41) ', &
                    ' Q(4-1) ','  Q(42) ',' Q(4-2) ','  Q(43) ',' Q(4-3) ','  Q(44) ',' Q(4-4) '/
  data Tens_name_T/   ' Sum_dm ','  T(00) ','  T_z   ','   T_x  ','    T_y ','  T_z2  ','  T_xz  ','  T_yz  ','  T_yz  ','T_x2-y2 '/
  data Tens_name_T_m/ 'Sum_dm_m',' T(00)_m','  T_z_m ','  T_x_m ','  T_y_m ',' T_z2_m ',' T_xz_m ',' T_yz_m ',' T_yz_m ','Tx2-y2_m'/

  mot = ' '
  
  Select case(index)
    case(0)
      mot = Tens_name_tot(i)
    case(1)
      mot = Tens_name_D(i)
    case(2)
      mot = Tens_name_I(i)
    case(3)
      mot = Tens_name_I_m(i)
    case(4)
      mot = Tens_name_Q(i)
    case(5)
      mot = Tens_name_T(i)
    case(6)
      mot = Tens_name_T_m(i)
  end select
  
  if( Dafs ) then
    jj = 2 * j - 1
  else
    jj = j
  endif
  
  Tens(jj) = real( Resul, db )
  
  if( Dafs ) then

    l = len_trim( mot )
    mot(l+1:l+2) = '_r'
    if( ninitlr > 1 ) call ad_sufix(Core_resolved,initlr,jseuil,Length_word,mot,nseuil)
    Tens_name(jj) = mot

    jj = jj + 1
    Tens(jj) = aimag( Resul )
    mot(l+2:l+2) = 'i'
    Tens_name(jj) = mot

  else

    if( ninitlr > 1 ) call ad_sufix(Core_resolved,initlr,jseuil,Length_word,mot,nseuil)
    Tens_name(jj) = mot

  endif

  return
end