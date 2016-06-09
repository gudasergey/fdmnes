! FDMNES subroutine

! Calculate the absorption cross sections and the RXS amplitudes

subroutine write_coabs(Allsite,angxyz,axyz,Base_spin,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
            Densite_atom,E_cut,Energ,Energphot,Extract,Epsii,Eseuil,Final_tddft, &
            f_avantseuil,Full_self_abs,Green_int,hkl_dafs,iabsorig,icheck,ie,ie_computer, &
            Int_tens,isigpi,isymeq,jseuil,ltypcal,Moyenne,mpinodee,Multipole,n_multi_run,n_oo,n_tens_max,natomsym,nbseuil, &
            ncolm,ncolr,ncolt,nenerg,ninit1,ninitlr,nomabs,nomfich,nomfich_cal_convt,nomfich_s,nphi_dafs,npldafs, &
            nphim,nplr,nplrm,nseuil,nspin,numat_abs,nxanout,pdp,phdafs,phdf0t,phdt,pol,poldafse,poldafss,Rot_int, &
            sec_atom,secdd_a,secdd_m_a,secdq_a,secdq_m_a,secdo_a,secdo_m_a, &
            secmd_a,secmd_m_a,secmm_a,secmm_m_a,secoo_a,secoo_m_a,secqq_a,secqq_m_a, &
            Self_abs,Spherical_signal, &
            Spherical_tensor,Spinorbite_p,Taux_eq,V0muf,Vecdafse,Vecdafss,Vec,Volume_maille,Xan_atom)

  use declarations
  implicit none

  integer:: he, hs, i, ia, iabsorig, ib, ic1, ic2, icheck, icn1, icn2, id, ie, &
    ie_computer, ig, ii, ind_mu, initlr, ip, ipl, ipldafs, iseuil, &
    isym, ixandafs, j, j1, je, jhe, jhs, jpl, js, jseuil, ke, ks, l, &
    ll, long, mpinodee, n_dim, n_tens_max, n_multi_run, n_oo, n_tens, n1, n2, na, &
    natomsym, nb, nbseuil, nc, nccm, ncolm, ncolr, ncolt, nenerg, ninit1, ninitlr, nl, np, npldafs, &
    nphim, nplt, nplr, nplrm, nseuil, nspin, numat_abs, nxanout, nw

  character(len=length_word):: nomab
  character(len=length_word), dimension(ncolm):: nomabs
  character(len=length_word), dimension(ncolm*ninitlr):: title
  character(len=13), dimension(nplrm):: ltypcal
  character(len=132) nomfich, nomfich_cal_convt, nomfich_s, nomfichdafst, nomficht
  character(len=2310) mot

  complex(kind=db):: cf, f_avantseuil, ph, ph_m, sec
  complex(kind=db), dimension(3):: plae, plas, uae, uas
  complex(kind=db), dimension(8*ninitlr):: compnum
  complex(kind=db), dimension(3,3):: secdd, secmd, secmm
  complex(kind=db), dimension(3,3,3):: secdq
  complex(kind=db), dimension(3,3,3,3):: secdo, secqq
  complex(kind=db), dimension(3,3,3,3,3,3):: Mat6
  complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secddia, secddia_m, secmdia, secmdia_m, secmmia, secmmia_m
  complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia, secdqia_m
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secdoia, secdoia_m, secqqia, secqqia_m
  complex(kind=db), dimension(3,3,ninitlr,0:mpinodee-1):: secdd_a, secdd_m_a, secmd_a, secmd_m_a, secmm_a, secmm_m_a
  complex(kind=db), dimension(3,3,3,ninitlr,0:mpinodee-1):: secdq_a, secdq_m_a
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:mpinodee-1):: secdo_a, secdo_m_a, secqq_a, secqq_m_a
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:mpinodee-1):: secoo_a, secoo_m_a
  complex(kind=db), dimension(npldafs):: phdtem, phdf0t1, phdt1
  complex(kind=db), dimension(3,nplrm):: pol
  complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdt
  complex(kind=db), dimension(natomsym,npldafs):: phdafs
  complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss

  complex(kind=db), dimension(:,:,:,:), allocatable :: ampldafs, ampldafsdd, ampldafsdo, ampldafsdq, ampldafsmd, &
    ampldafsmm, ampldafsoo, ampldafsqq
  complex(kind=db), dimension(:,:,:,:,:), allocatable :: mu, mudd, mudo, mudq, mumd, mumm, muoo, muqq
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable :: secooia, secooia_m

  integer, dimension(natomsym):: isymeq
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(npldafs,2):: isigpi

  logical Allsite, Base_spin, Cartesian_tensor, Cor_abs, Core_resolved, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, E_vec, Dafs, &
    Dafs_bio, Final_tddft, Energphot, Extract, Full_self_abs, Green_int, Green_int_mag, idafs, M1M1, Magn_sens, &
    Moyenne, mu_cal, Self_abs, Spherical_signal, Spherical_tensor, Spinorbite_p, Tens_comp, Tensor_eval, Xan_atom

  logical, dimension(10):: Multipole

  real(kind=db):: ang, c_micro, c_milli, cst, conv_mbarn_nelec, dang, Densite_atom, E_cut, &
    eph2, Ephoton, Ephseuil, V0muf, Volume_maille

  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(ninitlr) :: ct_nelec, Epsii
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(3):: angxyz, axyz, voae, voas
  real(kind=db), dimension(3,3):: matopsym, Rot_int
  real(kind=db), dimension(nenerg) :: Energ
  real(kind=db), dimension(ninitlr) :: sec_atom
  real(kind=db), dimension(3,nplrm) :: vec
  real(kind=db), dimension(nplrm,2) :: pdp
  real(kind=db), dimension(3,npldafs):: hkl_dafs
  real(kind=db), dimension(ncolm*ninitlr) :: tens
  real(kind=db), dimension(natomsym) :: Taux_eq
  real(kind=db), dimension(n_tens_max*ninitlr,0:natomsym):: Int_tens
  real(kind=db), dimension(ncolr,ninitlr,0:natomsym):: secabs, secabsdd, secabsdq, secabsdo, secabsmd, secabsmm, secabsoo, &
        secabsqq
  real(kind=db), dimension(3,npldafs,nphim):: vecdafse, vecdafss

  if( icheck > 0 ) write(3,110)

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  if( E3E3 ) then
    allocate( secooia(3,n_oo,3,n_oo,ninitlr,0:natomsym) )
    allocate( secooia_m(3,n_oo,3,n_oo,ninitlr,0:natomsym) )
  endif

  if( E1E1 ) secabsdd(:,:,:) = ( 0._db, 0._db )
  if( E2E2 ) secabsqq(:,:,:) = ( 0._db, 0._db )
  if( E1E3 ) secabsdo(:,:,:) = ( 0._db, 0._db )
  if( E3E3 ) secabsoo(:,:,:) = ( 0._db, 0._db )
  if( M1M1 ) secabsmm(:,:,:) = ( 0._db, 0._db )
  if( E1M1 ) secabsmd(:,:,:) = ( 0._db, 0._db )
  if( E1E2 ) secabsdq(:,:,:) = ( 0._db, 0._db )

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

  if( Allsite ) then
    na = natomsym
    nb = natomsym
  else
    na = 0
    nb = 1
  endif

  if( dafs ) then

    allocate( ampldafs(npldafs,nphim,ninitlr,0:natomsym) )
    if( E1E1 ) allocate( ampldafsdd(npldafs,nphim,ninitlr,0:natomsym) )
    if( E1E2 ) allocate( ampldafsdq(npldafs,nphim,ninitlr,0:natomsym) )
    if( E2E2 ) allocate( ampldafsqq(npldafs,nphim,ninitlr,0:natomsym) )
    if( E1E3 ) allocate( ampldafsdo(npldafs,nphim,ninitlr,0:natomsym) )
    if( E3E3 ) allocate( ampldafsoo(npldafs,nphim,ninitlr,0:natomsym) )
    if( M1M1 ) allocate( ampldafsmm(npldafs,nphim,ninitlr,0:natomsym) )
    if( E1M1 ) allocate( ampldafsmd(npldafs,nphim,ninitlr,0:natomsym) )
    if( Cor_abs ) then
      allocate( mu(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E1E1 ) allocate( mudd(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E1E2 ) allocate( mudq(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E2E2 ) allocate( muqq(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E1E3 ) allocate( mudo(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E3E3 ) allocate( muoo(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( M1M1 ) allocate( mumm(npldafs,nphim,2,ninitlr,0:natomsym) )
      if( E1M1 ) allocate( mumd(npldafs,nphim,2,ninitlr,0:natomsym) )
    endif

  endif

  do initlr = 1,ninitlr       ! ----------> Boucle sur les seuils ou les etats initiaux

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
! Pour les seuils de tres basse Energie:
    Ephoton = max(0.001_db/Rydb, Ephoton)
    if( Energphot ) Ephseuil = Ephoton

    ct_nelec(initlr) = conv_mbarn_nelec(Ephoton)
    eph2 = 0.5 * Ephoton**2
! Pour avoir les tenseurs et sections efficace en Megabarn
    cst = eph2 / ct_nelec(initlr)

! Les tenseurs sont convertis en megabarn
    if( .not. Extract ) then
      if( xan_atom ) sec_atom(initlr) = sec_atom(initlr) * cst
      if( E1E1 ) secdd_a(:,:,initlr,ie_computer) = secdd_a(:,:,initlr,ie_computer) * cst
      if( E1E2 ) secdq_a(:,:,:,initlr,ie_computer) = secdq_a(:,:,:,initlr,ie_computer) * cst
      if( E2E2 ) secqq_a(:,:,:,:,initlr,ie_computer) = secqq_a(:,:,:,:,initlr,ie_computer) * cst
      if( E1E3 ) secdo_a(:,:,:,:,initlr,ie_computer) = secdo_a(:,:,:,:,initlr,ie_computer) * cst
      if( E3E3 ) secoo_a(:,:,:,:,initlr,ie_computer) = secoo_a(:,:,:,:,initlr,ie_computer) * cst
      if( E1M1 ) secmd_a(:,:,initlr,ie_computer) = - secmd_a(:,:,initlr,ie_computer) * cst
      if( M1M1 ) secmm_a(:,:,initlr,ie_computer) = secmm_a(:,:,initlr,ie_computer) * cst
      if( Green_int_mag ) then
        if( E1E1 ) secdd_m_a(:,:,initlr,ie_computer) = secdd_m_a(:,:,initlr,ie_computer) * cst
        if( E1E2 ) secdq_m_a(:,:,:,initlr,ie_computer) = secdq_m_a(:,:,:,initlr,ie_computer) * cst
        if( E2E2 ) secqq_m_a(:,:,:,:,initlr,ie_computer) = secqq_m_a(:,:,:,:,initlr,ie_computer) * cst
        if( E1E3 ) secdo_m_a(:,:,:,:,initlr,ie_computer) = secdo_m_a(:,:,:,:,initlr,ie_computer) * cst
        if( E3E3 ) secoo_m_a(:,:,:,:,initlr,ie_computer) = secoo_m_a(:,:,:,:,initlr,ie_computer) * cst
        if( E1M1 ) secmd_m_a(:,:,initlr,ie_computer) = - secmd_m_a(:,:,initlr,ie_computer) * cst
        if( M1M1 ) secmm_m_a(:,:,initlr,ie_computer) = secmm_m_a(:,:,initlr,ie_computer) * cst
      endif
    endif

    do ia = 1,natomsym

      isym = abs( isymeq(ia) )
      call opsym(isym,matopsym)
      if( base_spin ) then
        matopsym = matmul( matopsym, rot_int )
        matopsym = matmul( transpose(rot_int), matopsym )
      endif

      if( E1E1 ) then
        secdd(:,:) = secdd_a(:,:,initlr,ie_computer)
        if( isym /= 1 ) call rot_tensor_2( secdd, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) secdd(:,:) = conjg( secdd(:,:) )
        secddia(:,:,initlr,ia) = secdd(:,:)
        if( Green_int_mag ) then
          secdd(:,:) = secdd_m_a(:,:,initlr,ie_computer)
          if( isym /= 1 ) call rot_tensor_2( secdd, matopsym )
          if( isymeq(ia) < 0 ) secdd(:,:) = - secdd(:,:)
          secddia_m(:,:,initlr,ia) = secdd(:,:)
        endif
      endif

      if( E1E2 ) then
        secdq(:,:,:) = secdq_a(:,:,:,initlr,ie_computer)
        if( isym /= 1 ) call rot_tensor_3( secdq, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) secdq(:,:,:) = conjg( secdq(:,:,:) )
        secdqia(:,:,:,initlr,ia) = secdq(:,:,:)
        if( Green_int_mag ) then
          secdq(:,:,:) = secdq_m_a(:,:,:,initlr,ie_computer)
          if( isym /= 1 ) call rot_tensor_3( secdq, matopsym )
          if( isymeq(ia) < 0 ) secdq(:,:,:) = - secdq(:,:,:)
          secdqia_m(:,:,:,initlr,ia) = secdq(:,:,:)
        endif
      endif

      if( E2E2 ) then
        secqq(:,:,:,:) = secqq_a(:,:,:,:,initlr,ie_computer)
        if( isym /= 1 ) call rot_tensor_4( secqq, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) secqq(:,:,:,:) = conjg( secqq(:,:,:,:) )
        secqqia(:,:,:,:,initlr,ia) = secqq(:,:,:,:)
        if( Green_int_mag ) then
          secqq(:,:,:,:) = secqq_m_a(:,:,:,:,initlr,ie_computer)
          if( isym /= 1 ) call rot_tensor_4( secqq, matopsym )
          if( isymeq(ia) < 0 ) secqq(:,:,:,:) = - secqq(:,:,:,:)
          secqqia_m(:,:,:,:,initlr,ia) = secqq(:,:,:,:)
        endif
      endif

      if( E1E3 ) then
        secdo(:,:,:,:) = secdo_a(:,:,:,:,initlr,ie_computer)
        if( isym /= 1 ) call rot_tensor_4( secdo, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) secdo(:,:,:,:) = conjg( secdo(:,:,:,:) )
        secdoia(:,:,:,:,initlr,ia) = secdo(:,:,:,:)
        if( Green_int_mag ) then
          secdo(:,:,:,:) = secdo_m_a(:,:,:,:,initlr,ie_computer)
          if( isym /= 1 ) call rot_tensor_4( secdo, matopsym )
          if( isymeq(ia) < 0 ) secdo(:,:,:,:) = - secdo(:,:,:,:)
          secdoia_m(:,:,:,:,initlr,ia) = secdo(:,:,:,:)
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
                Mat6(:,je,he,:,js,hs) = secoo_a(:,jhe,:,jhs,initlr,ie_computer)
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
                  Mat6(:,je,he,:,js,hs) = secoo_m_a(:,jhe,:,jhs,initlr,ie_computer)
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
        secmd(:,:) = secmd_a(:,:,initlr,ie_computer)
        if( isym /= 1 ) call rot_tensor_2( secmd, matopsym )
   ! C'est la partie reelle qui porte le magnetisme
        if( isymeq(ia) < 0 .and. .not. Green_int ) secmd(:,:) = - conjg( secmd(:,:) )
        secmdia(:,:,initlr,ia) = secmd(:,:)
        if( Green_int_mag ) then
          secmd(:,:) = secmd_m_a(:,:,initlr,ie_computer)
          if( isym /= 1 ) call rot_tensor_2( secmd, matopsym )
          if( isymeq(ia) < 0 ) secmd(:,:) = - secmd(:,:)
          secmdia_m(:,:,initlr,ia) = secmd(:,:)
        endif
      endif

      if( M1M1 ) then
        secmm(:,:) = secmm_a(:,:,initlr,ie_computer)
        if( isym /= 1 ) call rot_tensor_2( secmm, matopsym )
        if( isymeq(ia) < 0 .and. .not. Green_int ) secmm(:,:) = conjg( secmm(:,:) )
        secmmia(:,:,initlr,ia) = secmm(:,:)
        if( Green_int_mag ) then
          secmm(:,:) = secmm_m_a(:,:,initlr,ie_computer)
          if( isym /= 1 ) call rot_tensor_2( secmm, matopsym )
          if( isymeq(ia) < 0 ) secmm(:,:) = - secmm(:,:)
          secmmia_m(:,:,initlr,ia) = secmm(:,:)
        endif
      endif

    end do ! fin boucle ia

  end do ! fin boucle initlr

  if( Dafs ) then
    phdt1(:) = phdt(:,1)
    phdf0t1(:) = phdf0t(:,1)
  endif

  if( Spherical_tensor ) call Spherical_tensor_cal(ct_nelec,Core_resolved,Densite_atom,E_cut,E1E1,E1E2,E1M1,E2E2, &
            Energ,Ephseuil,Epsii,Eseuil(nbseuil),icheck,ie,Int_tens,jseuil,Magn_sens,Moyenne,n_tens_max, &
            natomsym,nenerg,ninit1,ninitlr,nomfich_s,nphim,npldafs,nplr,nplrm,nplt,nseuil,numat_abs,pdp,phdafs,phdf0t1, &
            phdt1,pol,Poldafse,Poldafss,secddia,secdqia,secmdia,secqqia,Spherical_signal,Taux_eq,V0muf,Vec,Vecdafse,Vecdafss)

  E_vec = E1E2 .or. E2E2 .or. E1E3 .or. E3E3 .or. E1M1 .or. M1M1

  jpl = 0

  if( Cor_abs ) then
    nplt = nplr + 3*npldafs
  else
    nplt = nplr + npldafs
  endif
  
  do ixandafs = 1,2

    do ipl = 1,nplt

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

        secddia(:,:,:,0) = (0._db,0._db)
        secdqia(:,:,:,:,0) = (0._db,0._db)
        secdqia_m(:,:,:,:,0) = (0._db,0._db)
        secqqia(:,:,:,:,:,0) = (0._db,0._db)
        secdoia(:,:,:,:,:,0) = (0._db,0._db)
        if( E3E3 ) secooia(:,:,:,:,:,0) = (0._db,0._db)
        secmdia(:,:,:,0) = (0._db,0._db)
        secmdia_m(:,:,:,0) = (0._db,0._db)
        secmmia(:,:,:,0) = (0._db,0._db)
        if( Green_int_mag ) then
          secddia_m(:,:,:,0) = (0._db,0._db)
          secqqia_m(:,:,:,:,:,0) = (0._db,0._db)
          secdoia_m(:,:,:,:,:,0) = (0._db,0._db)
          if( E3E3 ) secooia_m(:,:,:,:,:,0) = (0._db,0._db)
          secmmia_m(:,:,:,0) = (0._db,0._db)
        endif

        do ia = 1,natomsym

          if( idafs ) then
! Le exp(iQr) est converti. On recupere le complexe conjugue dans convolution.
! Phdafs contient deja Taux_eq.
            ph = conjg( phdafs(ia,ipldafs) )
          else
            ph = (1._db, 0._db) * Taux_eq(ia)
          endif
          ph_m = img * ph

          if( E1E1 ) secddia(:,:,:,0) = secddia(:,:,:,0) + ph * secddia(:,:,:,ia)
          
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
          ! la partie magnetique est sur la composante reelle
            if( Green_int ) then
              secmdia(:,:,:,0) = secmdia(:,:,:,0) + ph * secmdia(:,:,:,ia)
            else
              secmdia(:,:,:,0) = secmdia(:,:,:,0) + ph_m * aimag( secmdia(:,:,:,ia) )
              if( Magn_sens ) secmdia_m(:,:,:,0) = secmdia_m(:,:,:,0) + ph * real( secmdia(:,:,:,ia), db )
            endif
          endif

          if( M1M1 ) secmmia(:,:,:,0) = secmmia(:,:,:,0) + ph * secmmia(:,:,:,ia)

          if( Green_int_mag ) then
            if( E1E1 ) secddia_m(:,:,:,0) = secddia_m(:,:,:,0) + ph * secddia_m(:,:,:,ia)
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
            call write_cartesian_tensor(Densite_atom,E_cut,E1E2,E2E2,Ephseuil,Epsii,Eseuil(nbseuil),ia,ie,ipldafs,jseuil, &
               M1M1,Magn_sens,natomsym,ninit1,ninitlr,nomfich_s,nseuil,numat_abs,secddia,secdqia, &
               secdqia_m,secqqia,secmdia,Tens_comp,V0muf,Core_resolved)
          end do
        endif

      endif

      if( ( icheck > 0 .and. ipl == 1 ) .or. ( icheck > 1 .and. idafs ) ) then

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

          if( Final_tddft ) then
            if( nbseuil == 2 ) then
              write(3,120) achar(nseuil+74) // achar(jseuil+iseuil+46) // achar(jseuil+iseuil+47)
            else
              write(3,120) achar(nseuil+74) // achar(jseuil+iseuil+46)
            endif
          elseif( nseuil == 0 ) then  ! optic
            if( Core_resolved ) then
              if( initlr == 1 ) then
                write(3,123) ninitlr
              else
                write(3,125) initlr-2, ninitlr
              endif
            else
              write(3,120) 'Opt'
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
              do ke = 1,3
                if( Green_int_mag ) then
                  write(3,150) secddia(ke,:,initlr,ia), secddia_m(ke,:,initlr,ia)
                elseif( Tens_comp ) then
                  write(3,150) secddia(ke,:,initlr,ia)
                else
                  write(3,150) real( secddia(ke,:,initlr,ia) )
                endif
              end do
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
                    write(3,150) real(secdqia(ke,ks,:,initlr,ia),db)
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
                      write(3,150) real( secqqia(ke,1:3,ks,js,initlr,ia) )
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
                elseif( tens_comp ) then
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

      endif

      if( idafs .or. mu_cal ) then
        np = nphi_dafs(ipldafs)
      else
        np = 1
      endif

      do ip = 1,np

        if( .not. ( idafs .or. mu_cal ) ) then

          plae(:) = pol(:,ipl)
          plas(:) = pol(:,ipl)
          if( E_vec ) voae(:) = vec(:,ipl)
          if( E_vec ) voas(:) = vec(:,ipl)

        elseif( idafs ) then

          plae(:) = poldafse(:,ipldafs,ip)
          plas(:) = poldafss(:,ipldafs,ip)
          if( E_vec ) voae(:) = vecdafse(:,ipldafs,ip)
          if( E_vec ) voas(:) = vecdafss(:,ipldafs,ip)

        else ! calcul mu

          if( Full_self_abs ) then

            select case(mod(ipldafs,4))
              case(1,0)
                if( ind_mu == 1 ) then   ! entrant
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
                if( ind_mu == 1 ) then   ! entrant
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
                if( ind_mu == 1 ) then   ! entrant
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

        do ia = 0,na
          do initlr = 1,ninitlr

            if( E1E1 ) then
              sec = (0._db,0._db)
              do ke = 1,3
                sec = sec + plae(ke) * sum( conjg(plas(:)) * secddia(:,ke,initlr,ia))
              end do
              if( Green_int_mag ) then
                do ke = 1,3
                  sec = sec + plae(ke) * sum( conjg(plas(:)) * secddia_m(:,ke,initlr,ia) )
                end do
              endif
! Il manque un facteur pi qui a deja ete pris en compte dans le calcul
! du tenseur dans la routine Tens_ab.
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                ampldafsdd(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                mudd(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                secabsdd(jpl,initlr,ia) = real( sec,db )
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

              sec = img * sec      ! C'est ici qu'on recupere le img
              if( idafs ) sec = - sec

              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                ampldafsdq(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                mudq(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                secabsdq(jpl,initlr,ia) = real( sec, db )
              endif
            endif

            if( E2E2 ) then
              sec = (0._db,0._db)
              do ke = 1,3
                do je = 1,3
                  do ks = 1,3
                    sec = sec + conjg( plas(ks) ) * plae(ke) * voae(je) * sum( voas(:) * secqqia(ks,:,ke,je,initlr,ia) )
                    if( .not. Green_int_mag ) cycle
                    sec = sec + conjg( plas(ks) ) * plae(ke) * voae(je) * sum( voas(:) * secqqia_m(ks,:,ke,je,initlr,ia) )
                  end do
                end do
              end do
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                ampldafsqq(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                muqq(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                secabsqq(jpl,initlr,ia) = real( sec,db )
              endif
            endif

            if( E1E3 ) then
              sec = (0._db,0._db)
              do ke = 1,3
                do ks = 1,3
                  do j1 = 1,3
                    sec = sec + conjg( plas(ks) ) * plae(ke) * ( voae(j1) * sum( voae(:) * secdoia(ks,ke,j1,:,initlr,ia) ) &
                        + voas(j1) * sum( voas(:) * conjg( secdoia(ke,ks,j1,:,initlr,ia) ) ) )
                    if( .not. Green_int_mag ) cycle
                    sec = sec + conjg( plas(ks) ) * plae(ke) * ( voae(j1) * sum( voae(:) &
                        * secdoia_m(ks,ke,j1,:,initlr,ia) ) + voas(j1) * sum( voas(:) &
                        * conjg( secdoia_m(ke,ks,j1,:,initlr,ia) )))
                  end do
                end do
              end do
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                ampldafsdo(ipldafs,ip,initlr,ia)=sec
              elseif( mu_cal ) then
                mudo(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                secabsdo(jpl,initlr,ia) = real( sec, db )
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
                          if( .not. Green_int_mag ) cycle
                          sec = sec +conjg( plas(ks) ) * voas(js) * voas(hs) * plae(ke) * voae(je) * voae(he) &
                            * secooia_m(ks,jhs,ke,jhe,initlr,ia)
                        end do
                      end do
                    end do
                  end do
                end do
              end do
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                ampldafsoo(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                muoo(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                secabsoo(jpl,initlr,ia) = real( sec,db )
              endif
            endif

            if( M1M1 ) then
              sec = (0._db,0._db)
              do ke = 1,3
                sec = sec + uae(ke) * sum( conjg(uas(:)) * secmmia(:,ke,initlr,ia) )
                if( .not. Green_int_mag ) cycle
                sec = sec + uae(ke) * sum( conjg(uas(:))*secmmia_m(:,ke,initlr,ia) )
              end do
! Il manque un facteur pi qui a deja ete pris en compte dans le calcul
! du tenseur dans la routine Tens_ab.
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                ampldafsmm(ipldafs,ip,initlr,ia)=sec
              elseif( mu_cal ) then
                mumm(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                secabsmm(jpl,initlr,ia) = real( sec, db )
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

! Il manque un facteur pi qui a deja ete pris en compte dans le calcul du tenseur dans la routine Tens_ab.
              if( idafs ) then
                if( Green_int ) sec = pi * img * conjg( sec )
                ampldafsmd(ipldafs,ip,initlr,ia) = sec
              elseif( mu_cal ) then
                mumd(ipldafs,ip,ind_mu,initlr,ia) = sec
              else
                secabsmd(jpl,initlr,ia) = real( sec, db )
              endif
            endif

          end do   ! Fin de la boucle sur initlr
        end do   ! Fin de la boucle sur les atomes

        if( ipl > nplr ) cycle
        if( ltypcal(ipl) == 'xanes circ d' ) jpl = jpl + 1

      end do

    end do   ! Fin de la boucle sur les polarisation
  end do

  if( Moyenne ) then
    if( Xan_atom ) then
      i = ncolr - 1
    else
      i = ncolr
    endif
    ipl = 0
    do j = 1,ncolr
      if( ipl > 1 ) then
        if( ltypcal(ipl) == 'xanes circ d' ) cycle
      endif
      ipl = ipl + 1
      if( E1E1 ) secabsdd(i,:,:) = secabsdd(i,:,:) + pdp(ipl,1) * secabsdd(j,:,:)
      if( E2E2 ) secabsqq(i,:,:) = secabsqq(i,:,:) + pdp(ipl,2) * secabsqq(j,:,:)
      if( E1E3 ) secabsdo(i,:,:) = secabsdo(i,:,:) + pdp(ipl,1) * secabsdo(j,:,:)
      if( E3E3 ) secabsoo(i,:,:) = secabsoo(i,:,:) + pdp(ipl,1) * secabsoo(j,:,:)
      if( M1M1 ) secabsmm(i,:,:) = secabsmm(i,:,:) + pdp(ipl,1) * secabsmm(j,:,:)
      if( E1M1 ) secabsmd(i,:,:) = secabsmd(i,:,:) + pdp(ipl,1) * secabsmd(j,:,:)
      if( ipl == nplr ) exit
    end do
    if( E1E2 ) secabsdq(i,:,:) = (0._db,0._db)
  endif

  jpl = 0
  do ipl = 1,nplr
    jpl = jpl + 1
    if( ltypcal(ipl) /= 'xanes circ d' ) cycle
    jpl = jpl + 1
    ig = jpl - 1
    id = jpl - 2
    if( E1E1 ) secabsdd(jpl,:,:) = secabsdd(id,:,:) - secabsdd(ig,:,:)
    if( E1E2 ) secabsdq(jpl,:,:) = secabsdq(id,:,:) - secabsdq(ig,:,:)
    if( E2E2 ) secabsqq(jpl,:,:) = secabsqq(id,:,:) - secabsqq(ig,:,:)
    if( E1E3 ) secabsdo(jpl,:,:) = secabsdo(id,:,:) - secabsdo(ig,:,:)
    if( E3E3 ) secabsoo(jpl,:,:) = secabsoo(id,:,:) - secabsoo(ig,:,:)
    if( M1M1 ) secabsmm(jpl,:,:) = secabsmm(id,:,:) - secabsmm(ig,:,:)
    if( E1M1 ) secabsmd(jpl,:,:) = secabsmd(id,:,:) - secabsmd(ig,:,:)
  end do

  if( xan_atom ) then
    secabsdd(ncolr,:,0) = sec_atom(:) * natomsym
    do ia = 1,na
      secabsdd(ncolr,:,ia) = sec_atom(:)
    end do
  endif

  secabs(:,:,:) = 0._db
  if( E1E1 ) secabs(:,:,:) = secabs(:,:,:) + secabsdd(:,:,:)
  if( E1E2 ) secabs(:,:,:) = secabs(:,:,:) + secabsdq(:,:,:)
  if( E2E2 ) secabs(:,:,:) = secabs(:,:,:) + secabsqq(:,:,:)
  if( E1E3 ) secabs(:,:,:) = secabs(:,:,:) + secabsdo(:,:,:)
  if( E3E3 ) secabs(:,:,:) = secabs(:,:,:) + secabsoo(:,:,:)
  if( M1M1 ) secabs(:,:,:) = secabs(:,:,:) + secabsmm(:,:,:)
  if( E1M1 ) secabs(:,:,:) = secabs(:,:,:) + secabsmd(:,:,:)

  if( Dafs ) then
    ampldafs(:,:,:,:) = (0._db,0._db)
    if( E1E1 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:) + ampldafsdd(:,:,:,:)
    if( E1E2 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:) + ampldafsdq(:,:,:,:)
    if( E2E2 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:) + ampldafsqq(:,:,:,:)
    if( E1E3 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:) + ampldafsdo(:,:,:,:)
    if( E3E3 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:) + ampldafsoo(:,:,:,:)
    if( M1M1 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:) + ampldafsmm(:,:,:,:)
    if( E1M1 ) ampldafs(:,:,:,:) = ampldafs(:,:,:,:) + ampldafsmd(:,:,:,:)
  endif

  if( Cor_abs ) then
    mu(:,:,:,:,:) = (0._db,0._db)
    if( E1E1 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mudd(:,:,:,:,:)
    if( E1E2 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mudq(:,:,:,:,:)
    if( E2E2 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + muqq(:,:,:,:,:)
    if( E1E3 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mudo(:,:,:,:,:)
    if( E3E3 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + muoo(:,:,:,:,:)
    if( M1M1 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mumm(:,:,:,:,:)
    if( E1M1 ) mu(:,:,:,:,:) = mu(:,:,:,:,:) + mumd(:,:,:,:,:)
  endif

! Conversion en nombre d'electrons
  if( dafs ) then
    do initlr = 1,ninitlr
      ampldafs(:,:,initlr,:) = ct_nelec(initlr) * ampldafs(:,:,initlr,:)
      if( E1E1 ) ampldafsdd(:,:,initlr,:) = ct_nelec(initlr) * ampldafsdd(:,:,initlr,:)
      if( E1E2 ) ampldafsdq(:,:,initlr,:) = ct_nelec(initlr) * ampldafsdq(:,:,initlr,:)
      if( E2E2 ) ampldafsqq(:,:,initlr,:) = ct_nelec(initlr) * ampldafsqq(:,:,initlr,:)
      if( M1M1 ) ampldafsmm(:,:,initlr,:) = ct_nelec(initlr) * ampldafsmm(:,:,initlr,:)
      if( E1M1 ) ampldafsmd(:,:,initlr,:) = ct_nelec(initlr) * ampldafsmd(:,:,initlr,:)
    end do
  endif

! Conversion en micrometres
  if( Cor_abs ) then
    c_micro = 100 / ( Volume_maille * bohr**3 )
    mu(:,:,:,:,:) = c_micro * mu(:,:,:,:,:)
    if( E1E1 ) mudd(:,:,:,:,:) = c_micro * mudd(:,:,:,:,:)
    if( E1E2 ) mudq(:,:,:,:,:) = c_micro * mudq(:,:,:,:,:)
    if( E2E2 ) muqq(:,:,:,:,:) = c_micro * muqq(:,:,:,:,:)
    if( E1E3 ) mudo(:,:,:,:,:) = c_micro * mudo(:,:,:,:,:)
    if( E3E3 ) muoo(:,:,:,:,:) = c_micro * muoo(:,:,:,:,:)
    if( M1M1 ) mumm(:,:,:,:,:) = c_micro * mumm(:,:,:,:,:)
    if( E1M1 ) mumd(:,:,:,:,:) = c_micro * mumd(:,:,:,:,:)
  endif
  if( nseuil == 0 ) then ! cas de l'optique: en millimetre^-1
    c_milli = 100000 / ( Volume_maille * bohr**3 )
    secabs(:,:,:) = c_milli * secabs(:,:,:)
    if( E1E1 ) secabsdd(:,:,:) = c_milli * secabsdd(:,:,:)
    if( E1E2 ) secabsdq(:,:,:) = c_milli * secabsdq(:,:,:)
    if( E2E2 ) secabsqq(:,:,:) = c_milli * secabsqq(:,:,:)
    if( E1E3 ) secabsdo(:,:,:) = c_milli * secabsdo(:,:,:)
    if( E3E3 ) secabsoo(:,:,:) = c_milli * secabsoo(:,:,:)
    if( M1M1 ) secabsmm(:,:,:) = c_milli * secabsmm(:,:,:)
    if( E1M1 ) secabsmd(:,:,:) = c_milli * secabsmd(:,:,:)
  endif

  if( icheck > 0 ) then
    do ia = 0,na
      if( ia == 0 ) then
        if( nseuil /= 0 ) then
          write(3,283) ct_nelec(:) * pi
        else
          write(3,284) c_milli
       endif
        write(3,285)
      else
        write(3,290) ia
      endif
      do ipl = 1,ncolt
        nomab = nomabs(ipl)
        call center_word(nomab,Length_word)
        nomabs(ipl) = nomab
      end do
      nccm = 36
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
          write(3,310) Ephseuil*rydb, secabs(ic1:ic2,initlr,ia)
          if( E1E1 .and. ( E1E2 .or. E2E2 .or. E1E3 .or. M1M1 .or. E3E3 .or. E1M1) ) write(3,320) secabsdd(ic1:ic2,initlr,ia)
          if( E1E2 .and. ( E1E1 .or. E2E2 .or. E1E3 .or. M1M1 .or. E3E3 .or. E1M1) ) write(3,330) secabsdq(ic1:ic2,initlr,ia)
          if( E2E2 .and. ( E1E1 .or. E2E2 .or. E1E3 .or. M1M1 .or. E3E3 .or. E1M1) ) write(3,340) secabsqq(ic1:ic2,initlr,ia)
          if( E1E3 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. M1M1 .or. E3E3 .or. E1M1) ) write(3,350) secabsdo(ic1:ic2,initlr,ia)
          if( E3E3 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. M1M1 .or. E1E3 .or. E1M1) ) write(3,351) secabsoo(ic1:ic2,initlr,ia)
          if( M1M1 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. E1E3 .or. E3E3 .or. E1M1) ) write(3,352) secabsmm(ic1:ic2,initlr,ia)
          if( E1M1 .and. ( E1E1 .or. E1E2 .or. E2E2 .or. E1E3 .or. E3E3 .or. M1M1) ) write(3,354) secabsmd(ic1:ic2,initlr,ia)
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
              write(3,370) (ampldafs(j,1,initlr,ia), Real(mu(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E1E1 .and. ( E1E2 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,320) (ampldafsdd(j,1,initlr,ia), Real(mudd(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E1E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,330) (ampldafsdq(j,1,initlr,ia), Real(mudq(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E2E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,340) (ampldafsqq(j,1,initlr,ia), Real(muqq(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E1E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) (ampldafsdo(j,1,initlr,ia), Real(mudo(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E3E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E1E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) (ampldafsoo(j,1,initlr,ia), Real(mudo(j,1,:,initlr,ia)), j = ic1,ic2)
              if( M1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. E1M1 ) ) &
                write(3,352) (ampldafsmm(j,1,initlr,ia), Real(mumm(j,1,:,initlr,ia)), j = ic1,ic2)
              if( E1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. M1M1 ) ) &
                write(3,354) (ampldafsmd(j,1,initlr,ia), Real(mumd(j,1,:,initlr,ia)), j = ic1,ic2)
            elseif(Full_self_abs ) then
              write(3,370) (ampldafs(j,1,initlr,ia), mu(j,1,:,initlr,ia), j = ic1,ic2)
              if( E1E1 .and. ( E1E2 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,320) (ampldafsdd(j,1,initlr,ia), mudd(j,1,:,initlr,ia), j = ic1,ic2)
              if( E1E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,330) (ampldafsdq(j,1,initlr,ia), mudq(j,1,:,initlr,ia), j = ic1,ic2)
              if( E2E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,340) (ampldafsqq(j,1,initlr,ia), muqq(j,1,:,initlr,ia), j = ic1,ic2)
              if( E1E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) (ampldafsdo(j,1,initlr,ia), mudo(j,1,:,initlr,ia), j = ic1,ic2)
              if( E3E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E1E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) (ampldafsoo(j,1,initlr,ia), mudo(j,1,:,initlr,ia), j = ic1,ic2)
              if( M1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. E1M1 ) ) &
                write(3,352) (ampldafsmm(j,1,initlr,ia), mumm(j,1,:,initlr,ia), j = ic1,ic2)
              if( E1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. M1M1 ) ) &
                write(3,354) (ampldafsmd(j,1,initlr,ia), mumd(j,1,:,initlr,ia), j = ic1,ic2)
            else
              write(3,370) ampldafs(ic1:ic2,1,initlr,ia)
              if( E1E1 .and. ( E1E2 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,320) ampldafsdd(ic1:ic2,1,initlr,ia)
              if( E1E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,330) ampldafsdq(ic1:ic2,1,initlr,ia)
              if( E2E2 .and. ( E1E1 .or.E2E2 .or. E1E3 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,340) ampldafsqq(ic1:ic2,1,initlr,ia)
              if( E1E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) ampldafsdo(ic1:ic2,1,initlr,ia)
              if( E3E3 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E1E3 .or. M1M1 .or. E1M1 ) ) &
                write(3,350) ampldafsoo(ic1:ic2,1,initlr,ia)
              if( M1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. E1M1 ) ) &
                write(3,352) ampldafsmm(ic1:ic2,1,initlr,ia)
              if( E1M1 .and. ( E1E1 .or.E1E2 .or. E2E2 .or. E3E3 .or. E1E3 .or. M1M1 ) ) &
                write(3,354) ampldafsmd(ic1:ic2,1,initlr,ia)
            endif
          end do
        endif
      end do
    end do
  endif

  if( ie == 1 ) write(6,392) nenerg
  n_dim = ncolm * ninitlr

  do ia = 0,na

    nomficht = nomfich
    nomfichdafst = nomfich

    if( ia > 0 ) then
      long = len_trim(nomficht)
      nomficht(long+1:long+5) = '_atom'
      call ad_number(ia,nomficht,132)
      nomfichdafst(long+1:long+5) = '_atom'
      call ad_number(ia,nomfichdafst,132)
    endif
    long = len_trim(nomficht)

    if( Final_tddft .and. .not. Extract ) then
      nomficht(long+1:long+6) = '_tddft'
      nomfichdafst(long+1:long+11) = '_tddft_scan'
    else
      nomfichdafst(long+1:long+5) = '_scan'
    end if

    if( n_multi_run > 1 ) then
      l = len_trim(nomficht)
      nomficht(l+1:l+1) = '_'
      call ad_number(iabsorig,nomficht,132)
      l = len_trim(nomfichdafst)
      nomfichdafst(l+1:l+1) = '_'
      call ad_number(iabsorig,nomfichdafst,132)
    endif
    l = len_trim(nomficht)
    nomficht(l+1:l+4) = '.txt'
    l = len_trim(nomfichdafst)
    nomfichdafst(l+1:l+4) = '.txt'

    if( ie == 1 .and. ia == 0 ) nomfich_cal_convt = nomficht

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
      Tens(n_tens+1:ipl) = secabs(nxanout:ncolr,initlr,ia)
      do ipldafs = 1,npldafs
        if( ia == 0 ) then
          cf = ampldafs(ipldafs,1,initlr,ia)
        else
! Le exp(iQr) est converti. On recupere le complexe conjugue dans convolution.
          cf = conjg( phdafs(ia,ipldafs) ) * ampldafs(ipldafs,1,initlr,ia)
        endif
        ipl = ipl + 1
        Tens(ipl) = real( cf,db )
        ipl = ipl + 1
        Tens(ipl) = aimag( cf )
        if( self_abs ) then
          do i = 1,2
            ipl = ipl + 1
            Tens(ipl) = real( mu(ipldafs,1,i,initlr,ia), db )
          end do
        elseif( Full_self_abs ) then
          do i = 1,2
            ipl = ipl + 1
            Tens(ipl) = real( mu(ipldafs,1,i,initlr,ia), db )
            ipl = ipl + 1
            Tens(ipl) = aimag( mu(ipldafs,1,i,initlr,ia) )
          end do
        endif
      end do

      n_tens = ipl

    end do

    if( ia == 0 ) then
      phdtem(:) = phdt(:,1)
    else
      phdtem(:) = phdafs(ia,:)
    endif

! Ecriture dans le fichier
    if( Full_self_abs .or. Self_abs ) then
      call write_out(angxyz,axyz,Densite_atom,f_avantseuil,E_cut, Ephseuil, &
            Epsii,Eseuil(nbseuil),Green_int,hkl_dafs,ie,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title, &
            npldafs,npldafs,3,npldafs,nseuil,numat_abs,phdtem,phdf0t1,tens,v0muf,Core_resolved,natomsym)
    else
      call write_out(rdum,rdum,Densite_atom,f_avantseuil,E_cut, Ephseuil, &
            Epsii,Eseuil(nbseuil),Green_int,rdum,ie,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title, &
            npldafs,npldafs,0,0,nseuil,numat_abs,phdtem,phdf0t1,tens,v0muf,Core_resolved,natomsym)
    endif

! Ecriture a l'ecran
    if( ia == 0 ) call write_out(rdum,rdum,Densite_atom,f_avantseuil,E_cut, Ephseuil, &
            Epsii,Eseuil(nbseuil),Green_int,rdum,ie,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title, &
            npldafs,npldafs,0,0,nseuil,-1,phdtem,phdf0t1,tens,v0muf,Core_resolved,natomsym)

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
          write(7,410) nint( hkl_dafs(:,ipl) ), isigpi(ipl,:)
        endif
        dang = 360._db / nphi_dafs(ipl)

        do ip = 1,nphi_dafs(ipl)
          ang = ( ip - 1 ) * dang
          if( ia == 0 ) then
            cf = (1._db,0._db)
          else
            cf = conjg( phdafs(ia,ipl) )
          endif
          nw = 0
          do initlr = 1,ninitlr
            nw = nw + 1
            compnum(nw) = cf * ampldafs(ipl,ip,initlr,ia)
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
          compnum(nw) = phdt(ipl,ip)
          if( Dafs_bio ) then
            write(7,420) compnum(1:nw)
          else
            write(7,430) ang, compnum(1:nw)
          endif
        end do
      end do
      close(7)
    endif

  end do  ! fin boucle sur atomes

  if( E3E3 ) deallocate( secooia, secooia_m ) 

  if( dafs ) then
    deallocate( ampldafs )
    if( E1E1 ) deallocate( ampldafsdd )
    if( E1E2 ) deallocate( ampldafsdq )
    if( E2E2 ) deallocate( ampldafsqq )
    if( E1E3 ) deallocate( ampldafsdo )
    if( E3E3 ) deallocate( ampldafsoo )
    if( M1M1 ) deallocate( ampldafsmm )
    if( E1M1 ) deallocate( ampldafsmd )
    if( Cor_abs ) then
      deallocate( mu )
      if( E1E1 ) deallocate( mudd )
      if( E1E2 ) deallocate( mudq )
      if( E2E2 ) deallocate( muqq )
      if( E1E3 ) deallocate( mudo )
      if( E3E3 ) deallocate( muoo )
      if( M1M1 ) deallocate( mumm )
      if( E1M1 ) deallocate( mumd )
    endif
  endif

  return
  110 format(/' ---- Coabs --------',100('-'))
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
  283 format(/' Conversion factor (numb. of electron/Mbarn) =',10f10.5)
  284 format(/' Output in mm^(-1), conversion factor (mm^(-1)/Mbarn) =',10f10.5)
  285 format(/'   Total signal')
  290 format(/'   Signal atom',i3)
  295 format(/'   Core state or edge',i3)
  297 format(/'   Contribution for valence state l =',i3)
  300 format(/4x,'Energy',320a13)
  310 format(f10.3,1p,320e13.5)
  320 format('     E1-E1',1p,320e13.5)
  330 format('     E1-E2',1p,320e13.5)
  340 format('     E2-E2',1p,320e13.5)
  350 format('     E1-E3',1p,320e13.5)
  351 format('     E3-E3',1p,320e13.5)
  352 format('     M1-M1',1p,320e13.5)
  354 format('     M1-E1',1p,320e13.5)
  360 format(/10x,1p,320a13)
  370 format('  Ampldafs',1p,320e13.5)
  392 format(' Number of Energies =',i5)
  400 format(i5,4x,' = Number of angles')
  405 format(f10.3,A)
  407 format(3i5,' = (h,k,l)')
  410 format(' (h,k,l) = ',3i3,', sigpi =',2i3)
  420 format(7x,1p,320e13.5)
  430 format(f7.1,1p,320e13.5)
end

!***********************************************************************

subroutine write_nrixs(All_nrixs,Allsite,Core_resolved, &
                  Densite_atom,E_cut,Energ,Energphot,Extract,Epsii,Eseuil,Final_tddft, &
                  f_avantseuil,Green_int,iabsorig,icheck,ie,ie_computer,l0_nrixs,lmax_nrixs,isymeq, &
                  jseuil,mpinodee,n_multi_run,natomsym,nbseuil,nenerg,ninit1,ninitlr,nomfich,nomfich_cal_convt, &
                  nq_nrixs,nseuil,nspinp,numat_abs,q_nrixs,S_nrixs,S_nrixs_l,S_nrixs_l_m,S_nrixs_m,Spinorbite,Taux_eq,v0muf)
                  
  use declarations
  implicit none

  integer:: i, ia, iabsorig, icheck, ie, ie_computer, initlr, iq, l0_nrixs, lmax_nrixs, jseuil, l, long, &
    mpinodee, n, n_dim, n_multi_run, n_tens, natomsym, nbseuil, nenerg, ninit1, ninitlr, npldafs, nq_nrixs, &
    nseuil, nspinp, numat_abs

  character(len=length_word):: nomab
  character(len=132) nomfich, nomfich_cal_convt, nomficht
  character(len=length_word), dimension(:), allocatable:: title

  complex(kind=db):: f_avantseuil
  complex(kind=db), dimension(0):: cdum

  integer, dimension(natomsym):: isymeq

  logical:: All_nrixs, Allsite, Core_resolved, Final_tddft, Energphot, Extract, Green_int, Green_int_mag, Magn_sens, Spinorbite

  real(kind=db):: Densite_atom, E_cut, Ephoton, Ephseuil, q, V0muf

  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(ninitlr) :: Epsii
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(natomsym):: Taux_eq
  real(kind=db), dimension(nq_nrixs,ninitlr,0:mpinodee-1):: S_nrixs, S_nrixs_m
  real(kind=db), dimension(nq_nrixs,l0_nrixs:lmax_nrixs,ninitlr,0:mpinodee-1):: S_nrixs_l, S_nrixs_l_m
  real(kind=db), dimension(nq_nrixs,ninitlr,0:natomsym):: S_nrixs_T
  real(kind=db), dimension(nq_nrixs,l0_nrixs:lmax_nrixs,ninitlr,0:natomsym):: S_nrixs_T_l
  real(kind=db), dimension(nq_nrixs):: q_nrixs
  real(kind=db), dimension(:), allocatable:: Tens

  if( icheck > 0 ) write(3,110)
  
  npldafs = 0

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

  do initlr = 1,ninitlr

    do ia = 1,natomsym
      S_nrixs_T(:,initlr,ia) = Taux_eq(ia) * S_nrixs(:,initlr,ie_computer)
      S_nrixs_T_l(:,:,initlr,ia) = Taux_eq(ia) * S_nrixs_l(:,:,initlr,ie_computer)
      
      if( Green_int_mag ) then
        if( isymeq(ia) < 0 ) then
          S_nrixs_T(:,initlr,ia) = S_nrixs_T(:,initlr,ia) - Taux_eq(ia) * S_nrixs_m(:,initlr,ie_computer)
          S_nrixs_T_l(:,:,initlr,ia) = S_nrixs_T_l(:,:,initlr,ia) - Taux_eq(ia) * S_nrixs_l_m(:,:,initlr,ie_computer)
        else
          S_nrixs_T(:,initlr,ia) = S_nrixs_T(:,initlr,ia) + Taux_eq(ia) * S_nrixs_m(:,initlr,ie_computer)
          S_nrixs_T_l(:,:,initlr,ia) = S_nrixs_T_l(:,:,initlr,ia) + Taux_eq(ia) * S_nrixs_l_m(:,:,initlr,ie_computer)
        endif
      endif
      
      S_nrixs_T(:,initlr,0) =  S_nrixs_T(:,initlr,0) + S_nrixs_T(:,initlr,ia)
      S_nrixs_T_l(:,:,initlr,0) =  S_nrixs_T_l(:,:,initlr,0) + S_nrixs_T_l(:,:,initlr,ia)
    end do

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
      call ad_number(ia,nomficht,132)
    endif
    long = len_trim(nomficht)

    if( Final_tddft .and. .not. Extract ) nomficht(long+1:long+6) = '_tddft'

    if( n_multi_run > 1 ) then
      long = len_trim(nomficht)
      nomficht(long+1:long+1) = '_'
      call ad_number(iabsorig,nomficht,132)
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
          q = q_nrixs(iq) / bohr + eps10
          n = int( q )   
          call ad_number(n,nomab,length_word)
          long = len_trim( nomab )
          long = long + 1
          nomab(long:long) = '.'
          if( n < 10 ) then
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

! Ecriture dans le fichier
    call write_out(rdum,rdum,Densite_atom,f_avantseuil,E_cut,Ephseuil, &
            Epsii,Eseuil(nbseuil),Green_int,rdum,ie,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,Title, &
            npldafs,npldafs,0,0,nseuil,numat_abs,cdum,cdum,Tens,V0muf,Core_resolved,natomsym)

! Ecriture a l'ecran
    if( ia == 0 ) call write_out(rdum,rdum,Densite_atom,f_avantseuil,E_cut,Ephseuil, &
            Epsii,Eseuil(nbseuil),Green_int,rdum,ie,jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,Title, &
            npldafs,npldafs,0,0,nseuil,-1,cdum,cdum,Tens,V0muf,Core_resolved,natomsym)

  end do
  
  deallocate( Tens, Title )
   
  return
  110 format(/' ---- Write_nrixs --------',100('-'))
end

!***********************************************************************

! Calcul du facteur de conversion Mbarn --> nbr. d'elec (divise par pi)

function conv_mbarn_nelec(Ephoton)

  use declarations
  implicit none
  
  real(kind=db):: conv_mbarn_nelec, cst, Eph2, Ephoton, ptrans  

! Calcul de la constante multiplicative.
! ptrans = S02 fixe a 1.
  ptrans = 1._db
! alfa_sf = e*e/(2*epsilon0*h*c) est la constante de structure fine.
  cst = quatre_pi * pi * ptrans * alfa_sf * Ephoton
! pour avoir le resultat en megabarn (10E-18 cm2)
  cst = 100 * bohr**2 * cst

  Eph2 = 0.5_db * Ephoton**2
! Constante multiplication pour avoir le resultat en nombre d'electron
  conv_mbarn_nelec = Eph2 / cst

  return
end

!***********************************************************************

subroutine write_cartesian_tensor(Densite_atom,E_cut,E1E2,E2E2,Ephseuil, Epsii,Eseuil,ia,ie,ipldafs,jseuil, &
                 M1M1,magn_sens,natomsym,ninit1,ninitlr,nomfich_s,nseuil,numat_abs,secddia,secdqia,secdqia_m, &
                 secqqia,secmdia,tens_comp,v0muf,Core_resolved)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( n_dim=10*(168+12) )

  character(len=132) nomficht, nomfich_s
  character(len=Length_word) mot
  character(len=Length_word), dimension(n_dim):: nomtens

  complex(kind=db):: zero_c
  complex(kind=db), dimension(1):: cdum
  complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secddia, secmdia
  complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia, secdqia_m
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secqqia

  logical:: Core_resolved, E1E2, E2E2, M1M1, magn_sens, tens_comp

  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(ninitlr):: Epsii
  real(kind=db), dimension(n_dim):: Tens

  nomficht = nomfich_s
  long = len_trim(nomficht)
  if( ia == 0 ) then
    nomficht(long+1:long+9) = '_car_xtal'
    if( ipldafs > 0 ) then
      nomficht(long+10:long+13) = '_rxs'
      call ad_number(ipldafs,nomficht,132)
    endif
  else
    nomficht(long+1:long+9) = '_car_atom'
    call ad_number(ia,nomficht,132)
  endif
  long = len_trim(nomficht)
  nomficht(long+1:long+4) = '.txt'

  index = 0

  do initlr = 1,ninitlr

    mot = ' '
    mot(1:2) = 'D_'
    do i = 1,3
      mot(3:3) = achar(i+119)
      do j = i,3
        index = index + 1
        mot(4:4) = achar(j+119)
        if( tens_comp ) mot(5:6) = '_r'
        Tens(index) = real( secddia(i,j,initlr,ia),db )
        nomtens(index) = mot
        if( tens_comp ) then
          index = index + 1
          mot(6:6) = 'i'
          nomtens(index) = mot
          Tens(index) = aimag( secddia(i,j,initlr,ia) )
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

  call write_out(rdum,rdum,Densite_atom,zero_c,E_cut,Ephseuil,Epsii,Eseuil,.false.,rdum,ie, &
           jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,nomtens,1,0,0,0,nseuil,numat_abs,cdum, &
           cdum,tens,v0muf,Core_resolved,0) 

  return
end

!***********************************************************************

subroutine write_out(angxyz,axyz,Densite_atom,f_avantseuil,E_cut,Ephseuil,Epsii,Eseuil,Green_int,hkl_dafs,ie, &
            jseuil,n_dim,n_tens,ninit1,ninitlr,nomficht,title,np,npp,nppa,npps,nseuil,numat,ph1, &
            ph2,tens,v0muf,Core_resolved,natomsym)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter(n_tens_max = 10000)

  character(len=132):: nomficht
  character(len=95):: mot1
  character(len=27):: mot2
  character(len=Length_word):: mot
  character(len=Length_word), dimension(n_dim):: title
  character(len=10+(n_tens-2*npp)*Length_word):: dummy

  complex(kind=db):: f_avantseuil
  complex(kind=db), dimension(np):: ph1, ph2

  logical:: Core_resolved, Gnuplot, Green_int

  real(kind=db):: f0_forward, fpp_avantseuil
  real(kind=db), dimension(nppa):: angxyz, axyz
  real(kind=db), dimension(ninitlr):: Epsii
  real(kind=db), dimension(n_dim):: Tens
  real(kind=db), dimension(3,npps):: hkl_dafs


  dummy = ' '

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

  if( ie == 1 ) then
    if( numat >= 0 ) then
      open(ipr, file = nomficht)
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
      
      mot1 = ' = E_edge, Z, n_edge, j_edge, Abs_before_edge, VO_interstitial, E_Fermi, ninitl, ninit1, Epsii'
      mot2 = ', Atom_density, f0_forward'

      Gnuplot = .false.

      if( Gnuplot ) then
        if( nseuil == 0 ) then
          write(ipr,110) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, 1, icor, Epsii(1)*rydb, Densite_atom, f0_forward, mot1, mot2
        elseif( ninitlr == 1 ) then
          write(ipr,110) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, mot2
        elseif( ninitlr == 2 ) then
          write(ipr,111) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        elseif( ninitlr == 4 ) then
          write(ipr,112) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        elseif( ninitlr == 6 ) then
          write(ipr,113) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        elseif( ninitlr == 8 ) then
          write(ipr,114) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        elseif( ninitlr == 10 ) then
          write(ipr,115) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        else     ! 14
          write(ipr,116) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        endif
      else
        if( nseuil == 0 ) then
          write(ipr,120) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, 1, icor, Epsii(1)*rydb, Densite_atom, f0_forward, mot1, mot2
        elseif( ninitlr == 1 ) then
          write(ipr,120) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, mot2
        elseif( ninitlr == 2 ) then
          write(ipr,121) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        elseif( ninitlr == 4 ) then
          write(ipr,122) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        elseif( ninitlr == 6 ) then
          write(ipr,123) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        elseif( ninitlr == 8 ) then
          write(ipr,124) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        elseif( ninitlr == 10 ) then
          write(ipr,125) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        else     ! 14
          write(ipr,126) Eseuil*rydb, numat, nseuil, jseuil, fpp_avantseuil, v0muf*rydb, &
            E_cut*rydb, ninitlr, icor, Epsii(:)*rydb, Densite_atom, f0_forward, mot1, ninitlr, mot2
        endif
      endif

    endif

    if( npp > 0 .and. numat >= 0 ) then
      select case(Length_word)
        case(11)
          write(ipr,141) dummy, ph2(1:npp)
          write(ipr,141) dummy, ph1(1:npp)
        case(12)
          write(ipr,142) dummy, ph2(1:npp)
          write(ipr,142) dummy, ph1(1:npp)
        case(13)
          write(ipr,143) dummy, ph2(1:npp)
          write(ipr,143) dummy, ph1(1:npp)
        case(14)
          write(ipr,144) dummy, ph2(1:npp)
          write(ipr,144) dummy, ph1(1:npp)
        case(15)
          write(ipr,145) dummy, ph2(1:npp)
          write(ipr,145) dummy, ph1(1:npp)
        case(16)
          write(ipr,146) dummy, ph2(1:npp)
          write(ipr,146) dummy, ph1(1:npp)
        case(17)
          write(ipr,147) dummy, ph2(1:npp)
          write(ipr,147) dummy, ph1(1:npp)
        case default
          call write_error
          do ipr = 6,9,3
            write(ipr,150) Length_word
          end do
          stop
      end select
      if( nppa > 0 ) write(ipr,155) natomsym, axyz(:)*bohr, angxyz(:), ( nint( hkl_dafs(:,i) ), i = 1,npp )
    endif
    do i = 1,n
      mot = title(i)
      call center_word( mot, Length_word )
      title(i) = mot
    end do
    write(ipr,160) title(1:n)
  elseif( numat >= 0 ) then
    open(ipr, file = nomficht, position='append')
  endif

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
    case default
      call write_error
      do ipr = 6,9,3
        write(ipr,150) Length_word
      end do
      stop
  end select

  if( numat >= 0 ) close(ipr)

  return
  105 format(//' The number of column to be written is ',i5, // &
               ' This is greater than the maximum possible value given in the routine write_out in the file coabs.f !', / & 
               ' To change that you must modify the formats 121 up to 157 in this routine,', / &
               ' and increase to the same value the parameter n_tens_max.', / &
               ' Then you compile again.' //)
  110 format('# ',f10.3,i5,2i3,3f12.5,2i3,2f13.3,f13.9,f12.5,a95,a27)
  111 format('# ',f10.3,i5,2i3,3f12.5,2i3,2f13.3,f13.9,f12.5,a95,'(1..',i1,')',a27)
  112 format('# ',f10.3,i5,2i3,3f12.5,2i3,4f13.3,f13.9,f12.5,a95,'(1..',i1,')',a27)
  113 format('# ',f10.3,i5,2i3,3f12.5,2i3,6f13.3,f13.9,f12.5,a95,'(1..',i1,')',a27)
  114 format('# ',f10.3,i5,2i3,3f12.5,2i3,8f13.3,f13.9,f12.5,a95,'(1..',i1,')',a27)
  115 format('# ',f10.3,i5,2i3,3f12.5,2i3,10f13.3,f13.9,f12.5,a95,'(1..',i2,')',a27)
  116 format('# ',f10.3,i5,2i3,3f12.5,2i3,14f13.3,f13.9,f12.5,a95,'(1..',i2,')',a27)

  120 format(f10.3,i5,2i3,3f12.5,2i3,f13.3,f13.9,f12.5,a95,a27)
  121 format(f10.3,i5,2i3,3f12.5,2i3,2f13.3,f13.9,f12.5,a95,'(1..',i1,')',a27)
  122 format(f10.3,i5,2i3,3f12.5,2i3,4f13.3,f13.9,f12.5,a95,'(1..',i1,')',a27)
  123 format(f10.3,i5,2i3,3f12.5,2i3,6f13.3,f13.9,f12.5,a95,'(1..',i1,')',a27)
  124 format(f10.3,i5,2i3,3f12.5,2i3,8f13.3,f13.9,f12.5,a95,'(1..',i1,')',a27)
  125 format(f10.3,i5,2i3,3f12.5,2i3,10f13.3,f13.9,f12.5,a95,'(1..',i2,')',a27)
  126 format(f10.3,i5,2i3,3f12.5,2i3,14f13.3,f13.9,f12.5,a95,'(1..',i2,')',a27)

  141 format(A,1p,10000e11.3)
  142 format(A,1p,10000e12.4)
  143 format(A,1p,10000e13.5)
  144 format(A,1p,10000e14.6)
  145 format(A,1p,10000e15.7)
  146 format(A,1p,10000e16.8)
  147 format(A,1p,10000e17.9)
  150 format(//' Length_word =',i3, ' This parameter must be set between 11 and 17 !'//)
  155 format(i4,3f10.5,3x,3f10.5,10000(14x,3i4))
  160 format('    Energy',10000A)
  161 format(f10.3,1p,10000e11.3)
  162 format(f10.3,1p,10000e12.4)
  163 format(f10.3,1p,10000e13.5)
  164 format(f10.3,1p,10000e14.6)
  165 format(f10.3,1p,10000e15.7)
  166 format(f10.3,1p,10000e16.8)
  167 format(f10.3,1p,10000e17.9)
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
    elseif( l + iu == Length-1 .and. nomfich(l:l) == '_' ) then
      nomfich(l:l) = achar(in+48)
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

subroutine Spherical_tensor_cal(ct_nelec,Core_resolved,Densite_atom,E_cut,E1E1,E1E2,E1M1,E2E2, &
            Energ,Ephseuil,Epsii,Eseuil,icheck,ie,Int_tens,jseuil,Magn_sens,Moyenne,n_tens_max, &
            natomsym,nenerg,ninit1,ninitlr,nomfich_s,nphim,npldafs,nplr,nplrm,nplt,nseuil,numat_abs,pdp,phdafs,phdf0t, &
            phdt,pol,Poldafse,Poldafss,secddia,secdqia,secmdia,secqqia,Spherical_signal,Taux_eq,V0muf,Vec,Vecdafse,Vecdafss)

  use declarations
  implicit none

  integer, parameter:: n_tens_dd = 9 
  integer, parameter:: n_tens_dq = 15 
  integer, parameter:: n_tens_qq = 25 
  integer, parameter:: n_tens_dm = 9 

  integer:: i, ib, icheck, ie, initlr, ipl, ipldafs, jseuil, n_tens_max, natomsym, nenerg, &
            ninit1, ninitlr, nphim, npldafs, nplr, nplrm, nplt, nseuil, numat_abs
  
  character(len=132):: nomfich_s
  character(len=58), dimension(9):: Com_dd, com_dm, Com_dm_m
  character(len=58), dimension(15):: Com_dq, Com_dq_m
  character(len=58), dimension(25):: Com_qq

  complex(kind=db):: ph
  complex(kind=db), dimension(n_tens_dd,ninitlr):: Sph_tensor_dd
  complex(kind=db), dimension(n_tens_dq,ninitlr):: Sph_tensor_dq, Sph_tensor_dq_m
  complex(kind=db), dimension(n_tens_qq,ninitlr):: Sph_tensor_qq
  complex(kind=db), dimension(n_tens_dm,ninitlr):: Sph_tensor_dm, Sph_tensor_dm_m
  complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secddia, secmdia
  complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia
  complex(kind=db), dimension(3):: Pl_i, Pl_s
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secqqia
  complex(kind=db), dimension(3,npldafs,nphim):: Poldafse, Poldafss
  complex(kind=db), dimension(3,nplrm):: Pol
  complex(kind=db), dimension(npldafs):: phdf0t, phdt
  complex(kind=db), dimension(natomsym,npldafs):: phdafs
  complex(kind=db), dimension(:,:), allocatable:: Tensor_pol_dd, Tensor_pol_dq, Tensor_pol_qq, Tensor_pol_dm

  logical:: Core_resolved, E1E1, E1E2, E1M1, E2E2, Magn_sens, Moyenne, Spherical_signal, Write_bav

  real(kind=db):: Densite_atom, dph, E_cut, Ephseuil, Eseuil, V0muf
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
                ' l = 4, non magn octupole:    Q(40)                     =', &
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
  call Abs_Spherical_tensor(Com_dm,Com_dm_m,Com_dd,Com_dq,Com_dq_m,Com_qq,ct_nelec,E1E1,E1E2,E1M1,E2E2,icheck,ie, &
       Magn_sens,natomsym,ninitlr,secddia,secdqia,secmdia,secqqia, &
       Sph_tensor_ia_dd,Sph_tensor_ia_dq,Sph_tensor_ia_dq_m,Sph_tensor_ia_qq,Sph_tensor_ia_dm,Sph_tensor_ia_dm_m,Taux_eq)

  call Write_Spherical_tensor(Core_resolved,Densite_atom,E_cut,E1E1,E1E2,E1M1,E2E2,Energ,Ephseuil,Epsii, &
       Eseuil,icheck,ie,Int_tens,jseuil,Magn_sens,n_tens_dd,n_tens_dm,n_tens_dq,n_tens_max,n_tens_qq, &
       natomsym,nenerg,ninit1,ninitlr,nomfich_s,nseuil,numat_abs, &
       Sph_tensor_ia_dd,Sph_tensor_ia_dq,Sph_tensor_ia_dq_m,Sph_tensor_ia_qq,Sph_tensor_ia_dm,Sph_tensor_ia_dm_m,V0muf)

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
  
  if( ie == 1 ) then
    
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

        ph = conjg( phdafs(ib,ipldafs) )

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
          dph = sum( abs( phdafs(:,ipldafs) - phdafs(:,ipldafs-1) ) )
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

    call Write_Signal(ct_nelec,Core_resolved,Densite_atom,E_cut,E1E1,E1E2,E1M1,E2E2,Ephseuil,Epsii, &
      Eseuil,0,ie,ipl,ipldafs,jseuil,Magn_sens,n_tens_dd,n_tens_dq,n_tens_dm, &
      n_tens_qq,ninit1,ninitlr,nomfich_s,npldafs,nplt,nseuil,numat_abs,phdf0t,phdt, &
      Sph_tensor_dd,Sph_tensor_dq,Sph_tensor_dq_m,Sph_tensor_dm,Sph_tensor_dm_m,Sph_tensor_qq, &
      Tensor_pol_dd,Tensor_pol_dq,Tensor_pol_dm,Tensor_pol_qq,V0muf)

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

subroutine Abs_Spherical_tensor(Com_dm,Com_dm_m,Com_dd,Com_dq,Com_dq_m,Com_qq,ct_nelec,E1E1,E1E2,E1M1,E2E2,icheck,ie, &
       Magn_sens,natomsym,ninitlr,secddia,secdqia,secmdia,secqqia, &
       Sph_tensor_ia_dd,Sph_tensor_ia_dq,Sph_tensor_ia_dq_m,Sph_tensor_ia_qq,Sph_tensor_ia_dm,Sph_tensor_ia_dm_m,Taux_eq)

  use declarations
  implicit none

  integer, parameter:: n_tens_dd = 9 
  integer, parameter:: n_tens_dq = 15 
  integer, parameter:: n_tens_qq = 25 
  integer, parameter:: n_tens_dm = 9 

  integer:: i, ia, ib, icheck, ie, initlr, natomsym, ninitlr
  
  character(len=58), dimension(9):: Com_dd, com_dm, Com_dm_m
  character(len=58), dimension(15):: Com_dq, Com_dq_m
  character(len=58), dimension(25):: Com_qq

  complex(kind=db), dimension(n_tens_dd):: Sph_tensor_dd
  complex(kind=db), dimension(n_tens_dq):: Sph_tensor_dq
  complex(kind=db), dimension(n_tens_qq):: Sph_tensor_qq
  complex(kind=db), dimension(n_tens_dm):: Sph_tensor_dm
  complex(kind=db), dimension(3,3):: secdd, secmd
  complex(kind=db), dimension(3,3,3):: secdq
  complex(kind=db), dimension(3,3,3,3):: secqq
  complex(kind=db), dimension(3,3,ninitlr,0:natomsym):: secddia, secmdia
  complex(kind=db), dimension(3,3,3,ninitlr,0:natomsym):: secdqia
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:natomsym):: secqqia

  logical:: E1E1, E1E2, E1M1, E2E2, Magn_sens

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
        secdd(:,:) = secddia(:,:,initlr,ia)
        call Sph_tensor_dd_cal(n_tens_dd,secdd,Sph_tensor_dd)
 ! conversion en nombre d'electrons
        Sph_tensor_ia_dd(:,ia,initlr) = ct_nelec(initlr) * pi * real( Sph_Tensor_dd(:), db )
        Sph_tensor_ia_dd(:,0,initlr) = Sph_tensor_ia_dd(:,0,initlr) + Taux_eq(ia) * Sph_tensor_ia_dd(:,ia,initlr) 
      endif

      if( E1E2 ) then
        secdq(:,:,:) = secdqia(:,:,:,initlr,ia)
        call Sph_tensor_dq_cal(n_tens_dq,secdq,Sph_tensor_dq)
        Sph_tensor_ia_dq(:,ia,initlr) = ct_nelec(initlr) * pi * real( Sph_Tensor_dq(:), db )
        if( Magn_sens ) Sph_tensor_ia_dq_m(:,ia,initlr) = ct_nelec(initlr) * pi * aimag( Sph_Tensor_dq(:) )
        Sph_tensor_ia_dq(:,0,initlr) = Sph_tensor_ia_dq(:,0,initlr) + Taux_eq(ia) * Sph_tensor_ia_dq(:,ia,initlr)
        if( Magn_sens ) Sph_tensor_ia_dq_m(:,0,initlr) = Sph_tensor_ia_dq_m(:,0,initlr) &
                                                       + Taux_eq(ia) * Sph_tensor_ia_dq_m(:,ia,initlr) 
      endif

      if( E2E2 ) then
        secqq(:,:,:,:) = secqqia(:,:,:,:,initlr,ia)
        call Sph_tensor_qq_cal(n_tens_qq,secqq,Sph_tensor_qq)
        Sph_tensor_ia_qq(:,ia,initlr) = ct_nelec(initlr) * pi * real( Sph_Tensor_qq(:), db )
        Sph_tensor_ia_qq(:,0,initlr) = Sph_tensor_ia_qq(:,0,initlr) + Taux_eq(ia) * Sph_tensor_ia_qq(:,ia,initlr)
      endif
     
      if( E1M1 ) then  ! pour E1M1, la partie magnetique est le terme reel
        secmd(:,:) = secmdia(:,:,initlr,ia)
        call Sph_tensor_dm_cal(n_tens_dm,secmd,Sph_tensor_dm)
        Sph_tensor_ia_dm(:,ia,initlr) = ct_nelec(initlr) * pi * aimag( Sph_tensor_dm(:) )
        if( Magn_sens ) Sph_tensor_ia_dm_m(:,ia,initlr) = ct_nelec(initlr) * pi * real( Sph_tensor_dm(:), db )
        Sph_tensor_ia_dm(:,0,initlr) = Sph_tensor_ia_dm(:,0,initlr) + Taux_eq(ia) * Sph_tensor_ia_dm(:,ia,initlr)
        if( Magn_sens ) Sph_tensor_ia_dm_m(:,0,initlr) = Sph_tensor_ia_dm_m(:,0,initlr) &
                                                       + Taux_eq(ia) * Sph_tensor_ia_dm_m(:,ia,initlr) 
      endif

    end do  ! fin boucle ia

  end do ! fin boucle initlr

  do ib = 1,natomsym+1
    if( icheck == 0 .or. ie > 1 ) exit
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

subroutine Sph_tensor_dd_cal(n_tens_dd,secdd,Sph_tensor_dd)

  use declarations
  implicit none
  
  integer:: n_tens_dd

  complex(kind=db), dimension(n_tens_dd):: Sph_tensor_dd
  complex(kind=db), dimension(3,3):: secdd

! l = 0

  Sph_tensor_dd(1) = ( 1 / 3._db ) * ( secdd(1,1) + secdd(2,2) + secdd(3,3) )

! l = 1

! Multiplie par - img, equivalent a prendre la partie imaginaire 
! lz = D(1,0)
  Sph_tensor_dd(2) = - img * ( secdd(1,2) - secdd(2,1) ) / 2

! lx = D(1,1)
  Sph_tensor_dd(3) = - img * ( secdd(2,3) - secdd(3,2) ) / 2

! ly = D(1,-1)
  Sph_tensor_dd(4) = - img * ( secdd(3,1) - secdd(1,3) ) / 2

! l = 2

! D02 = ( 2*Dzz - Dxx - Dyy ) / 6
  Sph_tensor_dd(5) =( 2*secdd(3,3) - secdd(1,1)  - secdd(2,2) ) / 6

! D(2,1)
  Sph_tensor_dd(6) = ( secdd(1,3) + secdd(3,1) ) / 2

! D(2,-1)
  Sph_tensor_dd(7) = ( secdd(2,3) + secdd(3,2) ) / 2

! D(2,-2)
  Sph_tensor_dd(8) = ( secdd(1,2) + secdd(2,1) ) / 2

! D(2,2)
  Sph_tensor_dd(9) = ( secdd(1,1) - secdd(2,2) ) / 2

  return
end

!***********************************************************************

subroutine Sph_tensor_dq_cal(n_tens_dq,secdq,Sph_tensor_dq)

  use declarations
  implicit none
  
  integer:: n_tens_dq

  complex(kind=db), dimension(3,3,3):: secdq
  complex(kind=db), dimension(n_tens_dq):: Sph_tensor_dq

  real(kind=db):: fac
  
! Tenseur 1

! I(10)
  Sph_tensor_dq(1) = ( 2 * secdq(3,3,3) - secdq(3,1,1) - secdq(3,2,2) + 3 * secdq(1,1,3) + 3 * secdq(2,2,3) ) / 15

! I(11)
  Sph_tensor_dq(2) = ( 2 * secdq(1,1,1) - secdq(1,2,2) - secdq(1,3,3) + 3 * secdq(2,1,2) + 3 * secdq(3,1,3) ) / 15

! I(1-1)
  Sph_tensor_dq(3) = ( 2 * secdq(2,2,2) - secdq(2,1,1) - secdq(2,3,3) + 3 * secdq(1,1,2) + 3 * secdq(3,2,3) ) / 15

! Tenseur 2

! -i*I(20)
  Sph_tensor_dq(4) = secdq(1,2,3) - secdq(2,1,3)

  fac = 1 / 3._db

! (-i/sqrt(2))*(I(21) - I(2-1))
  Sph_tensor_dq(5) = fac * ( secdq(2,1,1) - secdq(2,3,3) - secdq(1,1,2) + secdq(3,2,3) )

! (1/sqrt(2))*(I(21) + I(2-1))
  Sph_tensor_dq(6) = fac * ( secdq(1,2,2) - secdq(1,3,3) - secdq(2,1,2) + secdq(3,1,3) )

  Sph_tensor_dq(7) = fac * ( secdq(1,1,3) - secdq(2,2,3) - secdq(3,1,1) + secdq(3,2,2) )

! (-i/sqrt(2))*(I(22) + I(2-2))
  Sph_tensor_dq(8) = fac * ( secdq(1,2,3) + secdq(2,1,3) - secdq(3,1,2) - secdq(3,2,1) )

! Tenseur 3

  Sph_tensor_dq(9) = 1 / sqrt(10._db) * ( 2 * secdq(3,3,3) - 2 * secdq(1,1,3) - secdq(3,1,1) &
                   - 2 * secdq(2,2,3) - secdq(3,2,2) )

  fac = 1 / sqrt(60._db)

  Sph_tensor_dq(10) = fac * ( 3 * secdq(1,1,1) + secdq(1,2,2) + 2 * secdq(2,1,2) - 4 * secdq(1,3,3) - 8 * secdq(3,1,3) )

! (-i/sqrt(2))*(I(31) + I(3-1))
  Sph_tensor_dq(11) = fac * ( 3 * secdq(2,2,2) + secdq(2,1,1) + 2 * secdq(1,2,1) - 4 * secdq(2,3,3) - 8 * secdq(3,2,3) )

  fac = 1 / sqrt(6._db)

! (-i/sqrt(2))*(I(32) - I(3-2))
  Sph_tensor_dq(12) = 2 * fac * ( secdq(1,2,3) + secdq(2,1,3) + secdq(3,1,2) )

  Sph_tensor_dq(13) = fac * ( 2 * secdq(1,1,3) - 2 * secdq(2,2,3) + secdq(3,1,1) - secdq(3,2,2) )

  fac = 0.5_db

  Sph_tensor_dq(14) = fac * ( secdq(1,2,2) + 2 * secdq(2,1,2) - secdq(1,1,1) )

! (-i/sqrt(2))*(I(33) + I(3-3))
  Sph_tensor_dq(15) = - fac * ( secdq(2,1,1) + 2 * secdq(1,2,1) - secdq(2,2,2) )

  return
end

!***********************************************************************

subroutine Sph_tensor_qq_cal(n_tens_qq,secqq,Sph_tensor_qq)

  use declarations
  implicit none
  
  integer:: n_tens_qq

  complex(kind=db), dimension(3,3,3,3):: secqq
  complex(kind=db), dimension(n_tens_qq):: Sph_tensor_qq

  real(kind=db):: fac
  
! La multiplication par - img correspond a la partie imaginaire quand le
! tenseur n'est pas multiplie par le terme de Bragg.

! Tenseur 0, scalaire, signal isotropique quadrupolaire

  fac = 1 / 45._db
! Q(00)
  Sph_tensor_qq(1) = fac * ( 6 * ( secqq(1,3,1,3) + secqq(2,3,2,3) + secqq(1,2,1,2) ) &
                           + 2 * ( secqq(1,1,1,1) + secqq(2,2,2,2) + secqq(3,3,3,3) ) &
                   - ( secqq(1,1,2,2) + secqq(1,1,3,3) + secqq(3,3,2,2) + secqq(2,2,1,1) + secqq(3,3,1,1) + secqq(2,2,3,3) ) )

! Tenseur 1, vecteur, lz, lx, ly, magnetique

  fac = 2 / 10._db

! Q(10)
  Sph_tensor_qq(2) = - img * fac * ( secqq(2,1,1,1) - secqq(1,2,2,2) + secqq(2,3,1,3) &
                                   - secqq(1,1,2,1) + secqq(2,2,1,2) - secqq(1,3,2,3) )

! Q(11)
  Sph_tensor_qq(3) = - img * fac * ( secqq(3,2,2,2) - secqq(2,3,3,3) + secqq(3,1,2,1) &
                                   - secqq(2,2,3,2) + secqq(3,3,2,3) - secqq(2,1,3,1) )

! Q(1-1)
  Sph_tensor_qq(4) = - img * fac * ( secqq(1,3,3,3) - secqq(3,1,1,1) + secqq(1,2,3,2) &
                                   - secqq(3,3,1,3) + secqq(1,1,3,1) - secqq(3,2,1,2) )

! Tenseur 2 : quadrupole non magnetique

  fac = 2 / ( 3 * sqrt(14._db) )
! Q(20)
  Sph_tensor_qq(5) = fac * ( 6 * secqq(1,2,1,2) - 3 * secqq(1,3,1,3) - 3 * secqq(2,3,2,3) + secqq(1,1,1,1) + secqq(2,2,2,2) &
           - 2 * secqq(1,1,2,2) - 2 * secqq(2,2,1,1) - 2 * secqq(3,3,3,3) + secqq(1,1,3,3) + secqq(2,2,3,3) &
           + secqq(3,3,1,1) + secqq(3,3,2,2) )

  fac = 2 / sqrt(42._db)

! (1/sqrt(2)) * ( Q(21) - Q(2-1) )_comp = Q(21) reel
  Sph_tensor_qq(6) = fac * ( 3 * secqq(2,3,1,2) + secqq(3,3,1,3) + secqq(1,1,1,3) - 2 * secqq(2,2,1,3) &
                     + 3 * secqq(1,2,2,3) + secqq(1,3,3,3) + secqq(1,3,1,1) - 2 * secqq(1,3,2,2) )

! (1/sqrt(2)) * ( Q(21) + Q(2-1) )
  Sph_tensor_qq(7) = fac * ( 3 * secqq(1,3,1,2) + secqq(3,3,2,3) + secqq(2,2,2,3) - 2 * secqq(1,1,2,3) &
                     + 3 * secqq(1,2,1,3) + secqq(2,3,3,3) + secqq(2,3,2,2) - 2 * secqq(2,3,1,1) )

  fac = 2 / sqrt( 42._db )

! (1/sqrt(2)) * ( Q(22) - Q(2-2) )
  Sph_tensor_qq(8) = fac * ( 2 * secqq(3,3,1,2) - secqq(1,1,1,2) - secqq(2,2,2,1) - 3 * secqq(1,3,2,3) &
             + 2 * secqq(1,2,3,3) - secqq(1,2,1,1) - secqq(2,1,2,2) - 3 * secqq(2,3,1,3) )

! (1/sqrt(2)) * ( Q(22) + Q(2-2) )
  Sph_tensor_qq(9) = fac * ( secqq(3,3,1,1) - secqq(3,3,2,2) + secqq(1,1,3,3) - secqq(2,2,3,3) &
             + secqq(2,2,2,2) - secqq(1,1,1,1) + 3 * secqq(2,3,2,3) - 3 * secqq(1,3,1,3) )

! Tenseur 3 : octupole magnetique

  fac = 1 / sqrt(10._db)

! Q(30)
  Sph_tensor_qq(10) = - img * fac * ( secqq(1,2,1,1) - secqq(2,1,2,2) + 4 * secqq(1,3,2,3) &
           - secqq(1,1,1,2) + secqq(2,2,2,1) - 4 * secqq(2,3,1,3) )

  fac = 0.5_db / sqrt(15._db)

! (1/sqrt(2)) * ( Q(31) - Q(3-1) )
  Sph_tensor_qq(11) = - img *  fac * ( 6 * secqq(1,2,1,3) + 4 * secqq(3,3,2,3) - 5 * secqq(1,1,2,3) + secqq(2,2,2,3) &
                      - 6 * secqq(1,3,1,2) - 4 * secqq(2,3,3,3) + 5 * secqq(2,3,1,1) - secqq(2,3,2,2) )

! (1/sqrt(2)) * ( Q(31) + Q(3-1) )
  Sph_tensor_qq(12) = - img * fac * ( 6 * secqq(1,2,2,3) + 4 * secqq(3,3,1,3) - 5 * secqq(2,2,1,3) + secqq(1,1,1,3) &
                      - 6 * secqq(2,3,1,2) - 4 * secqq(1,3,3,3) + 5 * secqq(1,3,2,2) - secqq(1,3,1,1) )

  fac = 1 / sqrt(6._db)

! (1/sqrt(2)) * ( Q(32) - Q(3-2) )
  Sph_tensor_qq(13) = - img * fac * ( secqq(1,1,3,3) + secqq(2,2,1,1) + secqq(3,3,2,2) &
               - secqq(3,3,1,1) - secqq(1,1,2,2) - secqq(2,2,3,3) )

! (1/sqrt(2)) * ( Q(32) + Q(3-2) )
  Sph_tensor_qq(14) = - img * fac * ( 2 * secqq(1,2,3,3) - secqq(1,2,1,1) - secqq(1,2,2,2) &
           - 2 * secqq(3,3,1,2) + secqq(1,1,1,2) + secqq(2,2,1,2) )

  fac = 0.5_db

! (1/sqrt(2)) * ( Q(33) - Q(3-3) )
  Sph_tensor_qq(15) = - img * fac * ( secqq(2,3,1,1) - secqq(2,3,2,2) + 2 * secqq(1,3,1,2) &
           - secqq(1,1,2,3) + secqq(2,2,2,3) - 2 * secqq(1,2,1,3) )

! (1/sqrt(2)) * ( Q(33) + Q(3-3) )
  Sph_tensor_qq(16) = - img * fac * ( secqq(1,3,1,1) - secqq(1,3,2,2) + 2 * secqq(2,1,3,2) &
           - secqq(1,1,1,3) + secqq(2,2,1,3) - 2 * secqq(3,2,2,1) )

! Tenseur 4 : hexadecapole non magnetique

  fac = 1 / ( 2. * sqrt(70._db) )

! Q(40)
  Sph_tensor_qq(17) = fac * ( 3 * secqq(1,1,1,1) + 3 * secqq(2,2,2,2) &
           + 8 * secqq(3,3,3,3) + secqq(1,1,2,2) + secqq(2,2,1,1) - 4 * secqq(3,3,1,1) - 4 * secqq(3,3,2,2) &
           - 4 * secqq(1,1,3,3) - 4 * secqq(2,2,3,3) + 4 * secqq(1,2,1,2) - 16 * secqq(1,3,1,3) - 16 * secqq(2,3,2,3) )

  fac = 0.5_db / sqrt(7._db)

! (1/sqrt(2)) * ( Q(41) - Q(4-1) )
  Sph_tensor_qq(18) = fac * ( 3 * secqq(1,3,1,1) + secqq(1,3,2,2) - 4 * secqq(3,3,1,3) + 2 * secqq(2,3,1,2) &
          + 3 * secqq(1,1,1,3) + secqq(2,2,1,3) - 4 * secqq(1,3,3,3) + 2 * secqq(1,2,2,3) )

! (1/sqrt(2)) * ( Q(41) + Q(4-1) )
  Sph_tensor_qq(19) = fac * ( 3 * secqq(2,3,2,2) + secqq(2,3,1,1) - 4 * secqq(2,3,3,3) + 2 * secqq(1,3,1,2) &
          + 3 * secqq(2,2,2,3) + secqq(1,1,2,3) - 4 * secqq(3,3,2,3) + 2 * secqq(1,2,1,3) )

  fac = 0.5_db / sqrt(14._db)

! (1/sqrt(2)) * ( Q(42) - Q(4-2) )
  Sph_tensor_qq(20) = 2 * fac * ( 2 * secqq(3,3,1,2) - secqq(1,1,1,2) - secqq(2,2,2,1) + 4 * secqq(1,3,2,3) &
          + 2 * secqq(1,2,3,3) - secqq(1,2,1,1) - secqq(2,1,2,2) + 4 * secqq(2,3,1,3) )

! (1/sqrt(2)) * ( Q(42) + Q(4-2) )
  Sph_tensor_qq(21) = fac * ( 2 * secqq(3,3,1,1) - 2 * secqq(3,3,2,2) + 2 * secqq(1,1,3,3) - 2 * secqq(2,2,3,3) &
        - 2 * secqq(1,1,1,1) + 2 * secqq(2,2,2,2) + 8 * secqq(1,3,1,3) - 8 * secqq(2,3,2,3) )

  fac = 0.5_db

! (1/sqrt(2)) * ( Q(43) - Q(4-3) )
  Sph_tensor_qq(22) = fac * ( secqq(1,3,2,2) - secqq(1,3,1,1) + 2 * secqq(2,3,1,2) &
           + secqq(2,2,1,3) - secqq(1,1,1,3) + 2 * secqq(1,2,2,3) )

! (1/sqrt(2)) * ( Q(43) + Q(4-3) )
  Sph_tensor_qq(23) = fac * ( secqq(2,3,2,2) - secqq(2,3,1,1) - 2 * secqq(1,3,1,2) &
           + secqq(2,2,2,3) - secqq(1,1,2,3) - 2 * secqq(1,2,1,3) )

  fac = 1 / sqrt(2._db)

! (1/sqrt(2)) * ( Q(44) - Q(4-4) )
  Sph_tensor_qq(24) = fac * ( secqq(1,2,1,1) - secqq(1,2,2,2) + secqq(1,1,1,2) - secqq(2,2,1,2) )

! (1/sqrt(2)) * ( Q(44) + Q(4-4) )
  Sph_tensor_qq(25) = 0.5_db * fac * ( secqq(1,1,1,1) + secqq(2,2,2,2) - secqq(1,1,2,2) &
             - secqq(2,2,1,1) - 4 * secqq(1,2,1,2) )

  return
end

!***********************************************************************

subroutine Sph_tensor_dm_cal(n_tens_dm,secmd,Sph_tensor_dm)

  use declarations
  implicit none
  
  integer:: n_tens_dm

  complex(kind=db), dimension(n_tens_dm):: Sph_tensor_dm
  complex(kind=db), dimension(3,3):: secmd

  real(kind=db):: fac
  
! Tenseur 0

  Sph_tensor_dm(1) = ( 1 / 3._db ) * ( secmd(1,1) + secmd(2,2) + secmd(3,3) )

! Tenseur 1
! Les composantes de ce tenseur sont en cas de seuil K : -lx, ly et lz.
  fac = 1 / 2._db

! Omega_z
  Sph_tensor_dm(2) = fac * ( secmd(1,2) - secmd(2,1) )

! Omega_x
  Sph_tensor_dm(3) = fac * ( secmd(2,3) - secmd(3,2) )

! Omega_y
  Sph_tensor_dm(4) = fac * ( secmd(3,1) - secmd(1,3) )

! Tenseur 2

! D02 = (1/sqrt(6))*(2*Dzz-Dxx-Dyy) = T_3z2-r2
  fac = 1 / 6._db
  Sph_tensor_dm(5) = fac * ( 2*secmd(3,3) - secmd(1,1)  - secmd(2,2) )

  fac = 1 / 2._db

! (-1/sqrt(2))*(D(12) - D(-12)) = T_xz
  Sph_tensor_dm(6) = fac * ( secmd(1,3) + secmd(3,1) )

! (i/sqrt(2))*(D(12) + D(-12)) = T_yz
  Sph_tensor_dm(7) = fac * ( secmd(2,3) + secmd(3,2) )

! (-i/sqrt(2))*(D(22) - D(-22)) = T_xy
  Sph_tensor_dm(8) = fac * ( secmd(1,2) + secmd(2,1) )

! (1/sqrt(2))*(D(22) + D(-22)) = T_x2-y2
  Sph_tensor_dm(9) = fac * ( secmd(1,1) - secmd(2,2) )

! Les tenseurs cartesiens sont definis par conjg(D).M, il faut conjg(M).D 
  Sph_tensor_dm(2:4) = - Sph_tensor_dm(2:4)
  Sph_tensor_dm(:) = Conjg( Sph_tensor_dm(:) ) 

  return
end

!***********************************************************************

! Ecriture des tenseurs spheriques et de leur integrale

subroutine Write_Spherical_tensor(Core_resolved,Densite_atom,E_cut,E1E1,E1E2,E1M1,E2E2,Energ,Ephseuil,Epsii, &
       Eseuil,icheck,ie,Int_tens,jseuil,Magn_sens,n_tens_dd,n_tens_dm,n_tens_dq,n_tens_max,n_tens_qq, &
       natomsym,nenerg,ninit1,ninitlr,nomfich_s,nseuil,numat_abs, &
       Sph_tensor_ia_dd,Sph_tensor_ia_dq,Sph_tensor_ia_dq_m,Sph_tensor_ia_dm,Sph_tensor_ia_dm_m,Sph_tensor_ia_qq,V0muf)

  use declarations
  implicit none
  
  integer:: i, ia, ie, icheck, initlr, j, jseuil, long, n_tens, n_tens_dd, n_tens_dq, n_tens_max, n_tens_dm, &
      n_tens_qq, natomsym, nenerg, ninit1, ninitlr, nseuil, numat_abs

  character(len=Length_word):: mot
  character(len=132):: nomfich_s, nomficht
  character(len=8), dimension(9):: Tens_name_D
  character(len=8), dimension(15):: Tens_name_I, Tens_name_I_m 
  character(len=8), dimension(25):: Tens_name_Q 
  character(len=8), dimension(9):: Tens_name_T, Tens_name_T_m 
  character(len=Length_word), dimension(n_tens_max*ninitlr):: nomten

  complex(kind=db):: zero_c
  complex(kind=db), dimension(1):: cdum

  logical:: Core_resolved, E1E1, E1E2, E1M1, E2E2, Magn_sens

  real(kind=db):: de, Densite_atom, E_cut, Ephseuil, Eseuil, V0muf
  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(n_tens_max*ninitlr,0:natomsym):: Int_tens
  real(kind=db), dimension(n_tens_max*ninitlr):: Int_tenst, Tens
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(ninitlr):: Epsii
  real(kind=db), dimension(n_tens_dd,0:natomsym,ninitlr):: Sph_tensor_ia_dd
  real(kind=db), dimension(n_tens_dq,0:natomsym,ninitlr):: Sph_tensor_ia_dq, Sph_tensor_ia_dq_m 
  real(kind=db), dimension(n_tens_qq,0:natomsym,ninitlr):: Sph_tensor_ia_qq
  real(kind=db), dimension(n_tens_dm,0:natomsym,ninitlr):: Sph_tensor_ia_dm, Sph_tensor_ia_dm_m

  data Tens_name_D/ '  D(00) ','  lz_dd ','  lx_dd ','  ly_dd ','  D_z2  ','   D_xy ','  D_yx  ','  D_yz  ',' D_x2-y2'/
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
      call ad_number(ia,nomficht,132)
    else
      nomficht(long+5:long+9) = '_xtal'
    endif
    long = len_trim(nomficht)
    nomficht(long+1:long+4) = '.txt'

    call write_out(rdum,rdum,Densite_atom,zero_c,E_cut,Ephseuil,Epsii,Eseuil,.false.,rdum,ie, &
          jseuil,n_tens_max*ninitlr,n_tens,ninit1,ninitlr,nomficht,nomten,1,0,0,0,nseuil,numat_abs,cdum, &
          cdum,Tens,V0muf,Core_resolved,0)

    if( nenerg == 1 ) cycle

! Integrale
  
    nomficht = nomfich_s
    long = len_trim(nomficht)
    nomficht(long+1:long+4) = '_sph'
    if( ia > 0 ) then
      nomficht(long+5:long+9) = '_atom'
      call ad_number(ia,nomficht,132)
    else
      nomficht(long+5:long+9) = '_xtal'
    endif
    long = len_trim(nomficht)
    nomficht(long+1:long+8) = '_int.txt'

    if( ie == 1 ) then
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
      de = Energ(2) - Energ(1)
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
    call write_out(rdum,rdum,Densite_atom,zero_c,E_cut,Ephseuil,Epsii,Eseuil,.false.,rdum,ie, &
          jseuil,n_tens_max*ninitlr,n_tens,ninit1,ninitlr,nomficht,nomten,1,0,0,0,nseuil,numat_abs,cdum, &
          cdum,Int_tenst,v0muf,Core_resolved,0)

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

subroutine Tensor_pol_dq_cal(ipl,n_tens_dq,nplt,pe,ps,ve,vs, Tensor_pol_dq)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db) Px, Py, Pz, Qx, Qy, Qz
  complex(kind=db), dimension(3):: pe, ps
  complex(kind=db), dimension(2*n_tens_dq):: Tens
  complex(kind=db), dimension(0:nplt,2*n_tens_dq):: Tensor_pol_dq

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
    Tens(j+1) = ( 2*Qz*Pz + 1.5*Qx*Px + 1.5*Qy*Py ) * ( Vz - Wz ) &
              + ( 1.5*Qx*Pz - Qz*Px ) * Vx - ( 1.5*Qz*Px - Qx*Pz ) * Wx &
              + ( 1.5*Qy*Pz - Qz*Py ) * Vy - ( 1.5*Qz*Py - Qy*Pz ) * Wy

! (T(11)
    Tens(j+2) = ( 2*Qx*Px + 1.5*Qy*Py + 1.5*Qz*Pz ) * ( Vx - Wx ) &
              + ( 1.5*Qy*Px - Qx*Py ) * Vy - ( 1.5*Qx*Py - Qy*Px ) * Wy &
              + ( 1.5*Qz*Px - Qx*Pz ) * Vz - ( 1.5*Qx*Pz - Qz*Px ) * Wz

! T(1-1)
    Tens(j+3) = ( 2*Qy*Py + 1.5*Qx*Px + 1.5*Qz*Pz ) * ( Vy - Wy ) &
              + ( 1.5*Qx*Py - Qy*Px ) * Vx - ( 1.5*Qy*Px - Qx*Py ) * Wx &
              + ( 1.5*Qz*Py - Qy*Pz ) * Vz - ( 1.5*Qy*Pz - Qz*Py ) * Wz

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
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db) Px, Py, Pz, Qx, Qy, Qz
  complex(kind=db), dimension(3):: pe, ps
  complex(kind=db), dimension(n_tens_qq):: Tens
  complex(kind=db), dimension(0:nplt,n_tens_qq):: Tensor_pol_qq

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

  Tens(5) = 3 * fac * ( Qx*Wy*Px*Vy + Qx*Wy*Py*Vx + Qy*Wx*Px*Vy + Qy*Wx*Py*Vx ) - 1.5 * fac &
     * ( Qx*Wz*Px*Vz + Qx*Wz*Pz*Vx + Qz*Wx*Px*Vz + Qz*Wx*Pz*Vx + Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy ) &
          + 2 * fac * ( Qx*Wx*Px*Vx + Qy*Wy*Py*Vy + Qx*Wx*Pz*Vz + Qz*Wz*Px*Vx + Qy*Wy*Pz*Vz + Qz*Wz*Py*Vy) - 4 * fac &
     * ( Qx*Wx*Py*Vy + Qy*Wy*Px*Vx + Qz*Wz*Pz*Vz )

  fac = 1 / sqrt( 42._db )

  Tens(6) = 1.5 * fac * ( Qy*Wz*Px*Vy + Qz*Wy*Px*Vy + Qy*Wz*Py*Vx + Qz*Wy*Py*Vx &
       + Qx*Wy*Py*Vz + Qx*Wy*Pz*Vy + Qy*Wx*Py*Vz + Qy*Wx*Pz*Vy ) + fac &
     * ( Qz*Wz*Px*Vz + Qz*Wz*Pz*Vx + Qx*Wx*Px*Vz + Qx*Wx*Pz*Vx + Qx*Wz*Pz*Vz + Qz*Wx*Pz*Vz + Qx*Wz*Px*Vx + Qz*Wx*Px*Vx ) &
          - 2 * fac * ( Qy*Wy*Px*Vz + Qy*Wy*Pz*Vx + Qx*Wz*Py*Vy + Qz*Wx*Py*Vy )

  Tens(7) = 1.5 * fac * ( Qx*Wz*Px*Vy + Qx*Wz*Py*Vx + Qz*Wx*Px*Vy + Qz*Wx*Py*Vx &
       + Qx*Wy*Px*Vz + Qx*Wy*Pz*Vx + Qy*Wx*Px*Vz + Qy*Wx*Pz*Vx ) + fac  &
     * ( Qz*Wz*Py*Vz + Qz*Wz*Pz*Vy + Qy*Wy*Py*Vz + Qy*Wy*Pz*Vy + Qy*Wz*Pz*Vz + Qz*Wy*Pz*Vz + Qy*Wz*Py*Vy + Qz*Wy*Py*Vy ) &
          - 2 * fac * ( Qx*Wx*Py*Vz + Qx*Wx*Pz*Vy + Qy*Wz*Px*Vx + Qz*Wy*Px*Vx )

  fac = 1 / sqrt( 42._db )

  Tens(8) = - fac * ( Qx*Wx*Px*Vy + Qx*Wx*Py*Vx + Qy*Wy*Px*Vy + Qy*Wy*Py*Vx &
       + Qx*Wy*Px*Vx + Qy*Wx*Px*Vx + Qx*Wy*Py*Vy + Qy*Wx*Py*Vy ) + 2 * fac &
     * ( Qz*Wz*Px*Vy + Qz*Wz*Py*Vx + Qx*Wy*Pz*Vz + Qy*Wx*Pz*Vz ) - 1.5 * fac &
     * ( Qx*Wz*Py*Vz + Qx*Wz*Pz*Vy + Qz*Wx*Py*Vz + Qz*Wx*Pz*Vy + Qy*Wz*Px*Vz + Qy*Wz*Pz*Vx + Qz*Wy*Px*Vz + Qz*Wy*Pz*Vx )

  Tens(9) = 2 * fac * ( Qz*Wz*Px*Vx - Qz*Wz*Py*Vy + Qx*Wx*Pz*Vz - Qy*Wy*Pz*Vz + Qy*Wy*Py*Vy - Qx*Wx*Px*Vx ) + 1.5 * fac &
     * ( Qy*Wz*Py*Vz + Qy*Wz*Pz*Vy + Qz*Wy*Py*Vz + Qz*Wy*Pz*Vy - Qx*Wz*Px*Vz - Qx*Wz*Pz*Vx - Qz*Wx*Px*Vz - Qz*Wx*Pz*Vx )

! Octupole magnetique

  fac = 1 / sqrt( 10._db )

  Tens(10) = 0.5 * fac * img * ( Qx*Wy*Px*Vx - Qx*Wy*Py*Vy + Qy*Wx*Px*Vx - Qy*Wx*Py*Vy &
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
           - 0.5 * fac * ( Qx*Wx*Px*Vy + Qx*Wx*Py*Vx + Qy*Wy*Px*Vy + Qy*Wy*Py*Vx &
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

subroutine Write_Signal(ct_nelec,Core_resolved,Densite_atom,E_cut,E1E1,E1E2,E1M1,E2E2,Ephseuil,Epsii, &
      Eseuil,ia,ie,ipl,ipldafs,jseuil,magn_sens,n_tens_dd,n_tens_dq,n_tens_dm, &
      n_tens_qq,ninit1,ninitlr,nomfich_s,npldafs,nplt,nseuil,numat_abs,phdf0t,phdt, &
      Sph_tensor_dd_ni,Sph_tensor_dq_ni,Sph_tensor_dq_m_ni,Sph_tensor_dm_ni,Sph_tensor_dm_m_ni,Sph_tensor_qq_ni, &
      Tensor_pol_dd,Tensor_pol_dq,Tensor_pol_dm,Tensor_pol_qq,V0muf)

  use declarations
  implicit none
  
  integer:: i, ia, ie, initlr, ipl, ipldafs, is, j, j0, jdd, jdd0, jdm, jdmm, jdm0, jdmm0, &
      jdq, jdq0, jdqm, jdqm0, jqq, jqq0, jseuil, jtot, jtot0, long, n_tens, n_tens_dd, n_tens_dq, n_tens_dm, &
      n_tens_qq, n_tens2, ninit1, ninitlr, npldafs, nplt, nseuil, numat_abs

  character(len=132):: nomfich_s, nomficht
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

  logical:: Core_resolved, Dafs, E1E1, E1E2, E1M1, E2E2, Magn_sens
  
  real(kind=db):: Densite_atom, E_cut, Ephseuil, Eseuil, fac, V0muf
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
      fac = 1 / ( ct_nelec(initlr) * pi )
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
    call ad_number(ia,nomficht,132)
  endif
  long = len_trim(nomficht)
  if( ipldafs > 0 ) then
    nomficht(long+1:long+4) = '_rxs'
    call ad_number(ipldafs,nomficht,132)
  elseif( ipl == 0 ) then
    nomficht(long+1:long+4) = '_xan'
  else
    nomficht(long+1:long+4) = '_pol'
    call ad_number(ipl,nomficht,132)
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
    call write_out(rdum,rdum,Densite_atom,zero_c,E_cut,Ephseuil,Epsii,Eseuil,.false.,rdum,ie, &
          jseuil,n_tens,n_tens,ninit1,ninitlr,nomficht,Tens_name,n_tens,n_tens2,0,0,nseuil,numat_abs,phtem, &
          ph0,Tens,v0muf,Core_resolved,0)
  else
    call write_out(rdum,rdum,Densite_atom,zero_c,E_cut,Ephseuil,Epsii,Eseuil,.false.,rdum,ie, &
          jseuil,n_tens,n_tens,ninit1,ninitlr,nomficht,Tens_name,1,0,0,0,nseuil,numat_abs,cdum,cdum,Tens,v0muf, &
          Core_resolved,0)
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
  data Tens_name_D/ ' Sum_dd ','  D(00) ','  lz_dd ','  lx_dd ','  ly_dd ','  D_z2  ','   D_xy ','  D_yx  ','  D_yz  ',' D_x2-y2'/
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