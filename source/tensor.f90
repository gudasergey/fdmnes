! FDMNES subroutines

! Calculation of the cartesian tensors
! Par rapport a coabs, les indices "e" et "s" sont inverses
! C'est l'indicage de coabs qui est logique

subroutine tenseur_car(Base_spin,coef_g,Core_resolved,Ecinetic, &
                Eimag,Energ,Enervide,Eseuil,Final_optic,Final_tddft,Full_potential,Green_int,Green_plus,Hubb_a,Hubb_d, &
                icheck,ie,ip_max,ip0,is_g,lmax,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
                n_oo,nbseuil,ndim1,ndim2,nenerg_tddft,ninit1,ninitl,ninitlr,ninitlt,ninitlv,nlm_pot,nlma,nlma2,nlmamax, &
                nr,nrm,nspin,nspino,nspinp,numat,psii,r,Relativiste,Rmtg,Rmtsd,rof0,rot_atom_abs,Rot_int, &
                secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m, &
                secqq,secqq_m,Solsing,Solsing_only,Spinorbite,Taull,Tddft,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer:: he, hhe, hhs, hs, i, icheck, ie, initlr, initlt, ip_max, ip0, ipr, irang, irang1, is, isi, isol, &
    isp, j, je, jhe, jhs, jje, jjhe, jjhs, jjs, jrang, js, k, ke, kke, kks, ks, le, lm, lmax, lmax_pot, lme, &
    lmomax, lms, lomax, ls, lseuil, m_hubb, me, mpinodes, mpirank, mpirank0, &
    ms, n_oo, nbseuil, ndim1, ndim2, ne, nenerg_tddft, ninit1, ninitl, &
    ninitlr, ninitlt, ninitlv, nhe, nje, nhs, njs, nlm_pot, nlm1g, &
    nlm2g, nlma, nlma2, nlmamax, nr, nrang, nrm, nspin, nspino, nspinp, numat

  parameter( lomax = 3, lmomax = ( lomax + 1 )**2 )

  complex(kind=db), dimension(ninitlr):: Te, Te_m, Ten, Ten_m, Tens, Tens_m
  complex(kind=db), dimension(3,3) :: Mat2
  complex(kind=db), dimension(3,3,3) :: Mat3
  complex(kind=db), dimension(3,3,3,3) :: Mat4
  complex(kind=db), dimension(3,3,3,3,3,3) :: Mat6
  complex(kind=db), dimension(3,3,ninitlr,0:mpinodes-1):: secdd, secdd_m, secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitlr,0:mpinodes-1):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:mpinodes-1):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:mpinodes-1):: secoo, secoo_m
  complex(kind=db), dimension(lmomax,lmomax,ninitlr):: Tens_lm, Tens_lm_m
  complex(kind=db), dimension(nlma*nspino,nlma*nspino,2,2,ndim1,ndim2,ndim2):: Taull
  complex(kind=db), dimension(nlma,nspinp,ip0:ip_max,ip0:ip_max,ninitlv):: Singul
  complex(kind=db), dimension(nlma,nlma2,nspinp,nspino,ip0:ip_max, ninitlv):: rof
  complex(kind=db), dimension(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil):: rof0
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  integer, dimension(ninitl,2):: m_g
  integer, dimension(ninitl):: is_g
  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3):: msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi

  logical Base_spin, Core_resolved, E1E1, E1E2, E1E3, E1M1, &
    E2E2, E3E3, Final_optic, Final_tddft, Full_potential, Green_int, &
    Green_plus, Hubb_a, Hubb_d, lmoins1, lplus1, M_depend, M1M1, &
    No_diag, Relativiste, Solsing, Solsing_only, Spinorbite, Tddft, Ylm_comp

  logical, dimension(10):: Multipole

  real(kind=db):: c, c0, c1, c12, c120, c3, c5, c8, clme, clms, Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic_e, V0bd_e
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitlr):: Tensi, Tensr
  real(kind=db), dimension(0:lomax,3,3,3,lmomax):: clm
  real(kind=db), dimension(3,3):: rot_atom_abs, Rot_int, rot_tem
  real(kind=db), dimension(ninitl,2):: coef_g
  real(kind=db), dimension(nspin,ninitlt):: Ecinetic, V0bd
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato_e
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nlm_pot,nspin,ninitlt):: Vrato
  real(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: roff_ii, roff_ir, roff_ri, roff_rr

  if( icheck > 1 ) write(3,100)

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

! ninitlt = ninitl si TDDFT et Core_resolved_t = nbseuil sinon
!         = 1 en DFT

! Calcul des integrales radiales et de la solution singuliere
  if( Final_optic ) then

    No_diag = Full_potential .or. ( Hubb_a .and. .not. Hubb_d )
    M_depend = Spinorbite .or. Hubb_a .or. No_diag

    ne = ninitlt

    if( No_diag ) then
      nlm1g = ( lmax + 1 )**2
      nlm2g = nlm1g
    elseif( M_depend ) then
      nlm1g = ( lmax + 1 )**2
      nlm2g = 1
    else
      nlm1g = lmax + 1
      nlm2g = 1
    endif
    allocate( roff_ii(nlm1g,nlm1g,nlm2g,nlm2g,nspinp**2,nspino**2, ip0:ip_max) )
    allocate( roff_ir(nlm1g,nlm1g,nlm2g,nlm2g,nspinp**2,nspino**2, ip0:ip_max) )
    allocate( roff_ri(nlm1g,nlm1g,nlm2g,nlm2g,nspinp**2,nspino**2, ip0:ip_max) )
    allocate( roff_rr(nlm1g,nlm1g,nlm2g,nlm2g,nspinp**2,nspino**2, ip0:ip_max) )

    call radial_optic(Ecinetic,Eimag,Energ, Full_potential,Green_plus,Hubb_a,Hubb_d,icheck,ip_max,ip0, &
         lmax,lmax_pot,m_depend,m_hubb,ne,nlm_pot,nlm1g,nlm2g,No_diag,nr,nspin,nspino,nspinp,numat,r,Relativiste,Rmtg,Rmtsd, &
         roff_ii,roff_ir,roff_ri,roff_rr,Solsing,Spinorbite,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp)

  else

    rof(:,:,:,:,:,:) = (0._db, 0._db)
    Singul(:,:,:,:,:) = (0._db, 0._db)

    do initlt = 1,ninitlt
      if( Final_tddft ) then
        do isp = 1,nspin
          Ecinetic_e(isp) = Ecinetic(isp,initlt)
          V0bd_e(isp) = V0bd(isp,initlt)
          Vrato_e(1:nr,:,isp) = Vrato(1:nr,:,isp,initlt)
        end do
        Enervide = Ecinetic(1,initlt) + V0bd(1,initlt)
      else
        do isp = 1,nspin
          Ecinetic_e(isp) = Ecinetic(isp,1)
          V0bd_e(isp) = V0bd(isp,1)
          Vrato_e(1:nr,:,isp) = Vrato(1:nr,:,isp,1)
        end do
      endif

      call radial(Ecinetic_e,Eimag,Energ,Enervide,Eseuil,Final_tddft,Full_potential,Green_plus,Hubb_a,Hubb_d,icheck, &
         initlt,ip_max,ip0,lmax,lmax_pot,m_hubb,nbseuil,ninit1,ninitlv,nlm_pot,nlma,nlma2,nr,nrm,nspin,nspino,nspinp,numat,psii, &
         r,Relativiste,Rmtg,Rmtsd,rof,Singul,Solsing,Spinorbite,V_hubb,V_intmax,V0bd_e,Vrato_e,Ylm_comp)

    end do

  endif

  if( icheck > 1 ) write(3,110)

  if( Tddft .and. .not. Final_tddft ) then
    do isol = 1,nspino
      do isp = 1,nspinp
        do lm = 1,nlma
          do i = 1,nbseuil
            if( nlma2 == 1 ) then
              rof0(ie,lm,isp,isol,:) = rof(lm,1,isp,isol,0,i)
            else
! ce qui sert de base en tddft, ce sont les (ls,ms,s) et pas les (l,m) vrais.
              rof0(ie,lm,isp,isol,:) = sum( rof(1:nlma,lm,isp,isol,0,i) )
            endif
          end do
        end do
      end do
    end do
  endif

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

  if( nrang > lomax .and. mpirank0 == 0 ) then
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
        clm(irang,1,:,:,4) = - c
        clm(irang,2,:,:,2) = - c
        clm(irang,3,:,:,3) = c

      case(2)
! Quadrupole
        c0 = sqrt( 4 * pi ) / 3
        c = sqrt( 4 * pi / 15 )
        c3 = c / sqrt( 3._db )

        clm(irang,1,1,:,1) = c0
        clm(irang,1,1,:,7) = - c3;   clm(irang,1,1,:,9) = c
        clm(irang,1,2,:,5) = c
        clm(irang,1,3,:,db) = - c
        clm(irang,2,2,:,1) = c0
        clm(irang,2,2,:,7) = - c3;   clm(irang,2,2,:,9) = - c
        clm(irang,2,3,:,6) = - c
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

        clm(irang,1,1,1,4) = - 3 * c1
        clm(irang,1,1,1,14) = 3 * c120; clm(irang,1,1,1,16) = - c8
        clm(irang,2,2,2,2) = - 3 * c1
        clm(irang,2,2,2,12) = 3 * c120; clm(irang,2,2,2,10) = c8
        clm(irang,3,3,3,3) = 3 * c1
        clm(irang,3,3,3,13) = 2 * c5
        clm(irang,1,1,2,2) = - c1
        clm(irang,1,1,2,10) = - c8;      clm(irang,1,1,2,12) = c120
        clm(irang,1,2,2,4) = - c1
        clm(irang,1,2,2,16) = c8;        clm(irang,1,2,2,14) = c120
        clm(irang,1,1,3,3) = c1
        clm(irang,1,1,3,13) = - c5;      clm(irang,1,1,3,15) = c12
        clm(irang,2,2,3,3) = c1
        clm(irang,2,2,3,13) = - c5;      clm(irang,2,2,3,15) = - c12
        clm(irang,2,3,3,2) = - c1
        clm(irang,2,3,3,12) = - 4 * c120
        clm(irang,1,3,3,4) = - c1
        clm(irang,1,3,3,14) = - 4 * c120
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

  secdd(:,:,:,mpirank) = (0._db,0._db)
  secmm(:,:,:,mpirank) = (0._db,0._db)
  secmd(:,:,:,mpirank) = (0._db,0._db)
  secdq(:,:,:,:,mpirank) = (0._db,0._db)
  secqq(:,:,:,:,:,mpirank) = (0._db,0._db)
  secdo(:,:,:,:,:,mpirank) = (0._db,0._db)
  if( E3E3 ) secoo(:,:,:,:,:,mpirank) = (0._db,0._db)

  if( Green_int ) then
    secdd_m(:,:,:,mpirank) = (0._db,0._db)
    secmm_m(:,:,:,mpirank) = (0._db,0._db)
    secmd_m(:,:,:,mpirank) = (0._db,0._db)
    secdq_m(:,:,:,:,mpirank) = (0._db,0._db)
    secqq_m(:,:,:,:,:,mpirank) = (0._db,0._db)
    secdo_m(:,:,:,:,:,mpirank) = (0._db,0._db)
    if( E3E3 ) secoo_m(:,:,:,:,:,mpirank) = (0._db,0._db)
  endif
! Boucles sur le rang des tenseurs

  if( icheck > 2 .and. .not. Final_optic ) write(3,120) lseuil

  do irang = irang1,nrang
    do jrang = irang,nrang
      Tens_lm(:,:,:) = (0._db,0._db)
      Tens_lm_m(:,:,:) = (0._db,0._db)
      if( .not. E1E1 .and. irang == 1 .and. jrang == 1 ) cycle
      if( .not. E1E2 .and. irang == 1 .and. jrang == 2 ) cycle
      if( .not. E2E2 .and. irang == 2 .and. jrang == 2 ) cycle
      if( .not. E1E3 .and. irang == 1 .and. jrang == 3 ) cycle
      if( .not. M1M1 .and. irang == 0 .and. jrang == 0 ) cycle
      if( .not. E1M1 .and. irang == 0 .and. jrang == 1 ) cycle
      if( .not. E3E3 .and. irang == 3 .and. jrang == 3 ) cycle
      if( irang == 0 .and. jrang > 1 ) cycle
      if( jrang == 3 .and. irang == 2 ) cycle

      if( icheck > 1 ) write(3,130) irang, jrang

! Boucles sur les indices des tenseurs

      do ke = 1,3

        if( irang == 1 .and. ldip(ke) /= 1 ) cycle

        if( irang > 1 ) then
          nje = 3
        else
          nje = 1
        endif

        do je = 1,nje

          if( irang == 2 .and. lqua(ke,je) /= 1 ) cycle

          if( irang == 3 ) then
            nhe = 3
          else
            nhe = 1
          endif

          do he = 1,nhe

            if( irang == 3 .and. loct(ke,je,he) /= 1 ) cycle
            jhe = 3 * ( je - 1 ) + he

            do ks = 1,3

              if( jrang == 1 .and. ldip(ks) /= 1 ) cycle

              if( jrang > 1 ) then
                njs = 3
              else
                njs = 1
              endif

              do js = 1,njs

                if( jrang == 2 .and. lqua(ks,js) /= 1 ) cycle

                if( jrang == 3 ) then
                  nhs = 3
                else
                  nhs = 1
                endif

                do hs = 1,nhs

                  if( jrang == 3 .and. loct(ks,js,hs) /= 1 ) cycle
                  jhs = 3 * ( js - 1 ) + hs

                  if( irang == 0 .and. jrang == 0 ) then
                    if( ks < ke ) cycle
                  elseif( irang == 1 .and. jrang == 1 ) then
                    if( ks < ke ) cycle
                    if( msymdd(ke,ks) == 0 .and. msymddi(ke,ks) == 0 ) cycle
                    if( sum( abs( secdd(ke,ks,:,mpirank) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 1 .and. jrang == 2 ) then
                    if( msymdq(ke,ks,js) == 0 .and. msymdqi(ke,ks,js) == 0 ) cycle
                    if( sum( abs( secdq(ke,ks,js,:,mpirank) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 2 .and. jrang == 2 ) then
                    if( msymqq(ke,je,ks,js) == 0 .and. msymqqi(ke,je,ks,js) == 0 ) cycle
                    if( sum( abs( secqq(ke,je,ks,js,:,mpirank) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 1 .and. jrang == 3 ) then
                    if( msymdo(ke,ks,js,hs) == 0 .and. msymdoi(ke,ks,js,hs) == 0 ) cycle
                  elseif( irang == 3 .and. jrang == 3 ) then
                    if( msymoo(ke,jhe,ks,jhs) == 0 .and. msymooi(ke,jhe,ks,jhs) == 0 ) cycle
                    if( sum( abs( secoo(ke,jhe,ks,jhs,:,mpirank) ) ) > 1.e-15_db ) cycle
                  endif

                  Tens(:) = (0._db,0._db)
                  Tens_m(:) = (0._db,0._db)

                  if( icheck > 1 ) write(3,140) ke, je, he, ks, js, hs

! Boucles sur les composantes spheriques des tenseurs
                  lme = 0
                  do le = 0,max(irang,1)
                    do me = -le,le
                      lme = lme + 1

                      clme = clm(irang,ke,je,he,lme)
                      if( abs(clme) < eps10 ) cycle

                      lms = 0
                      do ls = 0,max(jrang,1)
                        do ms = -ls,ls
                          lms = lms + 1

                          clms = clm(jrang,ks,js,hs,lms)
                          if( abs(clms) < eps10 ) cycle

! Calcul de la composante du tenseur
                          if( sum( abs( Tens_lm(lme,lms,:) ) ) < 1.e-15_db ) then

                            if( icheck > 1 ) write(3,150) le, me, ls, ms, clme, clms

                            if( Final_optic ) then
                              call tens_op(Core_resolved,icheck,ip_max,ip0,irang,jrang,le,me,ls,ms,lmax,lmoins1,lplus1,M_depend, &
                                ndim1,ndim2,ninitlr,nlm1g,nlm2g,nlma,nlma2,nspinp,nspino,roff_ii,roff_ir,roff_ri,roff_rr, &
                                Spinorbite,Taull,Ten,Ten_m,Ylm_comp)
                            else
                              call tens_ab(coef_g,Core_resolved,Final_tddft,Green_int,icheck,ip_max,ip0,irang,is_g,jrang, &
                                le,me,ls,m_g,ms,lmax,lmoins1,lplus1,lseuil,ndim1,ndim2,ninit1, &
                                ninitl,ninitlv,ninitlr,nlma,nlma2,nspinp,nspino,rof,Singul,Solsing, &
                                Solsing_only,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)
                            endif

                            Tens_lm(lme,lms,:) = Ten(:)
                            Tens_lm_m(lme,lms,:) = Ten_m(:)
                          endif

                          Tens(:) = Tens(:) + clme * clms * Tens_lm(lme,lms,:)
                          Tens_m(:) = Tens_m(:) + clme * clms * Tens_lm_m(lme,lms,:)
                        end do
                      end do
                    end do
                  end do ! fin boucle le

! Remplissage de la valeur calculee dans tous les elements du tenseur
! equivalent par symetrie.

                  Tensr(:) = real( Tens(:),db )
                  Tensi(:) = aimag( Tens(:) )

! M1-M1 (dipole magnetique - dipole magnetique)
                  if( irang == 0 .and. jrang == 0 ) then
                    secmm(ke,ks,:,mpirank) = Tens(:)
                    if( Green_int ) then
                      secmm(ks,ke,:,mpirank) = Tens(:)
                      secmm_m(ke,ks,:,mpirank) = Tens_m(:)
                      secmm_m(ks,ke,:,mpirank) = - Tens_m(:)
                    else
                      secmm(ks,ke,:,mpirank) = conjg( Tens(:) )
                    endif
                  endif

! M1-E1 (dipole magnetique - dipole electrique)
                  if( irang == 0 .and. jrang == 1 ) then
                    secmd(ke,ks,:,mpirank) = Tens(:)
                    if( Green_int ) secmd_m(ke,ks,:,mpirank) = Tens_m(:)
                  endif

! E1-E1 (Dipole-dipole)
                  if( irang == 1 .and. jrang == 1 ) then

                    do kke = 1,3
                      do kks = kke,3
                        if( abs(msymdd(kke,kks)) /= abs(msymdd(ke,ks)) .or. abs(msymddi(kke,kks)) &
                           /= abs(msymddi(ke,ks)) ) cycle
                        if( msymdd(kke,kks) /= 0 ) then
                          is = msymdd(kke,kks) / msymdd(ke,ks)
                        else
                          is = 0
                        endif
                        if( msymddi(kke,kks) /= 0 ) then
                          isi = msymddi(kke,kks) / msymddi(ke,ks)
                        else
                          isi = 0
                        endif
                        if( Green_int ) then
                          secdd(kke,kks,:,mpirank) = is * Tens(:)
                          secdd(kks,kke,:,mpirank) = is * Tens(:)
                          secdd_m(kke,kks,:,mpirank) = isi*Tens_m(:)
                          secdd_m(kks,kke,:,mpirank) =-isi*Tens_m(:)
                        else
! tenseur hermitique
                          Te(:) = cmplx(is*Tensr(:),isi*Tensi(:),db)
                          secdd(kke,kks,:,mpirank) = Te(:)
                          secdd(kks,kke,:,mpirank) = conjg( Te(:) )
                        endif
                      end do
                    end do

! Dipole-quadrupole
                  elseif( irang == 1 .and. jrang == 2 ) then
                    do kke = 1,3
                      do kks = 1,3
                        do jjs = kks,3
                          if( abs(msymdq(kke,kks,jjs)) /= abs(msymdq(ke,ks,js)) .or. abs(msymdqi(kke,kks,jjs)) &
                             /= abs(msymdqi(ke,ks,js)) ) cycle
                          if( msymdq(kke,kks,jjs) /= 0 ) then
                            is = msymdq(kke,kks,jjs) / msymdq(ke,ks,js)
                          else
                            is = 0
                          endif
                          if( msymdqi(kke,kks,jjs) /= 0 ) then
                            isi = msymdqi(kke,kks,jjs) / msymdqi(ke,ks,js)
                          else
                            isi = 0
                          endif
                          if( Green_int ) then
                            secdq(kke,kks,jjs,:,mpirank) =is*Tens(:)
                            secdq(kke,jjs,kks,:,mpirank) =is*Tens(:)
                            secdq_m(kke,kks,jjs,:,mpirank) = isi*Tens_m(:)
                            secdq_m(kke,jjs,kks,:,mpirank) = isi*Tens_m(:)
                          else
                            Te(:)=cmplx(is*Tensr(:),isi*Tensi(:),db)
                            secdq(kke,kks,jjs,:,mpirank) = Te(:)
                            secdq(kke,jjs,kks,:,mpirank) = Te(:)
                          endif
                        end do
                      end do
                    end do

! Dipole-octupole
                  elseif( irang == 1 .and. jrang == 3 ) then
                    do kke = 1,3
                      do kks = 1,3
                        do jjs = 1,3
                          do hhs = 1,3
                            if( abs(msymdo(kke,kks,jjs,hhs)) /= abs(msymdo(ke,ks,js,hs)) .or. abs(msymdoi(kke,kks,jjs,hhs)) &
                               /= abs(msymdoi(ke,ks,js,hs))) cycle
                            if( msymdo(kke,kks,jjs,hhs)/=0 ) then
                              is = msymdo(kke,kks,jjs,hhs) / msymdo(ke,ks,js,hs)
                            else
                              is = 0
                            endif
                            if( msymdoi(kke,kks,jjs,hhs)/=0) then
                              isi = msymdoi(kke,kks,jjs,hhs) / msymdoi(ke,ks,js,hs)
                            else
                              isi = 0
                            endif
                            if( Green_int ) then
                              secdo(kke,kks,jjs,hhs,:,mpirank) = is * Tens(:)
                              secdo(kke,kks,hhs,jjs,:,mpirank) = is * Tens(:)
                              secdo(kke,jjs,kks,hhs,:,mpirank) = is * Tens(:)
                              secdo(kke,hhs,kks,jjs,:,mpirank) = is * Tens(:)
                              secdo(kke,jjs,hhs,kks,:,mpirank) = is * Tens(:)
                              secdo(kke,hhs,jjs,kks,:,mpirank) = is * Tens(:)
                              secdo_m(kke,kks,jjs,hhs,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kke,kks,hhs,jjs,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kke,jjs,kks,hhs,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kke,hhs,kks,jjs,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kke,jjs,hhs,kks,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kke,hhs,jjs,kks,:,mpirank) = isi * Tens_m(:)
                            else
                              Te(:) = cmplx(is*Tensr(:),isi*Tensi(:),db)
                              secdo(kke,kks,jjs,hhs,:,mpirank)=Te(:)
                              secdo(kke,kks,hhs,jjs,:,mpirank)=Te(:)
                              secdo(kke,jjs,kks,hhs,:,mpirank)=Te(:)
                              secdo(kke,hhs,kks,jjs,:,mpirank)=Te(:)
                              secdo(kke,jjs,hhs,kks,:,mpirank)=Te(:)
                              secdo(kke,hhs,jjs,kks,:,mpirank)=Te(:)
                            endif
                          end do
                        end do
                      end do
                    end do

! Quadrupole-quadrupole
                  elseif( irang == 2 .and. jrang == 2 ) then
                    do kke = 1,3
                      do jje = 1,3
                        do kks = 1,3
                          do jjs = 1,3
                            if( sum( abs( secqq(kke,jje,kks,jjs,:,mpirank) ) ) > 1.e-15_db ) cycle
                            if( abs(msymqq(kke,jje,kks,jjs)) /= abs(msymqq(ke,je,ks,js)) .or. abs(msymqqi(kke,jje,kks,jjs)) &
                               /= abs(msymqqi(ke,je,ks,js)) ) cycle
                            if( msymqq(kke,jje,kks,jjs) /= 0 ) then
                              is = msymqq(kke,jje,kks,jjs) / msymqq(ke,je,ks,js)
                            else
                              is = 0
                            endif
                            if( msymqqi(kke,jje,kks,jjs) /= 0 ) then
                              isi = msymqqi(kke,jje,kks,jjs) / msymqqi(ke,je,ks,js)
                            else
                              isi = 0
                            endif
                            if( Green_int ) then
                              Te(:) = is * Tens(:)
                              secqq(kke,jje,kks,jjs,:,mpirank)=Te(:)
                              secqq(jje,kke,kks,jjs,:,mpirank)=Te(:)
                              secqq(kke,jje,jjs,kks,:,mpirank)=Te(:)
                              secqq(jje,kke,jjs,kks,:,mpirank)=Te(:)
                              secqq(kks,jjs,kke,jje,:,mpirank)=Te(:)
                              secqq(jjs,kks,kke,jje,:,mpirank)=Te(:)
                              secqq(kks,jjs,jje,kke,:,mpirank)=Te(:)
                              secqq(jjs,kks,jje,kke,:,mpirank)=Te(:)
                              Te(:) = isi * Tens_m(:)
                              secqq_m(kke,jje,kks,jjs,:,mpirank) = Te(:)
                              secqq_m(jje,kke,kks,jjs,:,mpirank) = Te(:)
                              secqq_m(kke,jje,jjs,kks,:,mpirank) = Te(:)
                              secqq_m(jje,kke,jjs,kks,:,mpirank) = Te(:)
                              secqq_m(kks,jjs,kke,jje,:,mpirank) = - Te(:)
                              secqq_m(jjs,kks,kke,jje,:,mpirank) = - Te(:)
                              secqq_m(kks,jjs,jje,kke,:,mpirank) = - Te(:)
                              secqq_m(jjs,kks,jje,kke,:,mpirank) = - Te(:)
                            else
                              Te(:)= cmplx( is*Tensr(:), isi*Tensi(:),db)
                              secqq(kke,jje,kks,jjs,:,mpirank)=Te(:)
                              secqq(jje,kke,kks,jjs,:,mpirank)=Te(:)
                              secqq(kke,jje,jjs,kks,:,mpirank)=Te(:)
                              secqq(jje,kke,jjs,kks,:,mpirank)=Te(:)
                              Te(:) = Conjg( Te(:) )
                              secqq(kks,jjs,kke,jje,:,mpirank)=Te(:)
                              secqq(jjs,kks,kke,jje,:,mpirank)=Te(:)
                              secqq(kks,jjs,jje,kke,:,mpirank)=Te(:)
                              secqq(jjs,kks,jje,kke,:,mpirank)=Te(:)
                            endif
                          end do
                        end do
                      end do
                    end do
! Octupole-octupole
                  elseif( irang == 3 .and. jrang == 3 ) then

                    do kke = 1,3
                      do jje = 1,3
                        do hhe = 1,3
                          jjhe = 3 * ( jje - 1 ) + hhe
                          do kks = 1,3
                            do jjs = 1,3
                              do hhs = 1,3
                                jjhs = 3 * ( jjs - 1 ) + hhs
                                i = msymoo(kke,jjhe,kks,jjhs)
                                j = msymoo(ke,jhe,ks,jhs)
                                if( abs(i) /= abs(j) ) cycle
                                if( j /= 0 ) then
                                  is = i / j
                                else
                                  is = 0
                                endif
                                i = msymooi(kke,jjhe,kks,jjhs)
                                j = msymooi(ke,jhe,ks,jhs)
                                if( abs(i) /= abs(j) ) cycle
                                if( j /= 0 ) then
                                  isi = i / j
                                else
                                  isi = 0
                                endif
                                if( Green_int ) then
                                  Te(:) = is * Tens(:)
                                  Te_m(:) = isi * Tens_m(:)
                                else
                                  Te(:)= cmplx( is * Tensr(:), isi * Tensi(:), db )
                                  Te_m(:) = Conjg( Te(:) )
                                endif

                                secoo(kke,jjhe,kks,jjhs,:,mpirank) = Te(:)
                                if( Green_int ) secoo_m(kke,jjhe,kks,jjhs,: ,mpirank) = Te_m(:)

                              end do
                            end do
                          end do
                        end do
                      end do
                    end do

                  endif

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
      lme = 0
      do le = 0,max(irang,1)
        do me = -le,le
         lme = lme + 1
         lms = 0
         do ls = 0,max(jrang,1)
            do ms = -ls,ls
              lms = lms + 1
              if( sum( abs(Tens_lm(lme,lms,:)) ) > eps15 ) write(3,170) irang, jrang, le, me, ls, ms, &
                   ( Tens_lm(lme,lms,initlr), initlr = 1,ninitlr )
              if( .not. Green_int ) cycle
              if( sum( abs(Tens_lm_m(lme,lms,:)) ) > eps15 ) write(3,180) irang, jrang, le, me, ls, ms, &
                   ( Tens_lm_m(lme,lms,initlr) ,initlr = 1,ninitlr )
            end do
          end do
        end do
      end do

    end do
  end do

! Rotation pour avoir les tenseurs dans la base R1

  if( .not. Base_spin ) then

    rot_tem = matmul( rot_int, transpose(rot_atom_abs) )

    do initlr = 1,ninitlr

      if( E1E1 ) then
        mat2(:,:) = secdd(:,:,initlr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secdd(:,:,initlr,mpirank) = mat2(:,:)
      endif

      if( E1E2 ) then
        mat3(:,:,:) = secdq(:,:,:,initlr,mpirank)
        call rot_tensor_3( mat3, rot_tem )
        secdq(:,:,:,initlr,mpirank) = mat3(:,:,:)
      endif

      if( E2E2 ) then
        mat4(:,:,:,:) = secqq(:,:,:,:,initlr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secqq(:,:,:,:,initlr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1E3 ) then
        mat4(:,:,:,:) = secdo(:,:,:,:,initlr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secdo(:,:,:,:,initlr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1M1 ) then
        mat2(:,:) = secmd(:,:,initlr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmd(:,:,initlr,mpirank) = mat2(:,:)
      endif

      if( M1M1 ) then
        mat2(:,:) = secmm(:,:,initlr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmm(:,:,initlr,mpirank) = mat2(:,:)
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
                Mat6(:,je,he,:,js,hs) = secoo(:,jhe,:,jhs,initlr,mpirank)
              end do
            end do
          end do
        end do
        call rot_tensor_6( Mat6, Rot_tem )
        jhe = 0
        do je = 1,3
          do he = 1,3
            jhe = jhe + 1
            jhs = 0
            do js = 1,3
              do hs = 1,3
                jhs = jhs + 1
                secoo(:,jhe,:,jhs,initlr,mpirank) = Mat6(:,je,he,:,js,hs)
              end do
            end do
          end do
        end do
      endif

      if( .not. Green_int ) Cycle

      if( E1E1 ) then
        mat2(:,:) = secdd_m(:,:,initlr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secdd_m(:,:,initlr,mpirank) = mat2(:,:)
      endif

      if( E1E2 ) then
        mat3(:,:,:) = secdq_m(:,:,:,initlr,mpirank)
        call rot_tensor_3( mat3, rot_tem )
        secdq_m(:,:,:,initlr,mpirank) = mat3(:,:,:)
      endif

      if( E2E2 ) then
        mat4(:,:,:,:) = secqq_m(:,:,:,:,initlr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secqq_m(:,:,:,:,initlr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1E3 ) then
        mat4(:,:,:,:) = secdo_m(:,:,:,:,initlr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secdo_m(:,:,:,:,initlr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1M1 ) then
        mat2(:,:) = secmd_m(:,:,initlr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmd_m(:,:,initlr,mpirank) = mat2(:,:)
      endif

      if( M1M1 ) then
        mat2(:,:) = secmm_m(:,:,initlr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmm_m(:,:,initlr,mpirank) = mat2(:,:)
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
                Mat6(:,je,he,:,js,hs) = secoo_m(:,jhe,:,jhs,initlr,mpirank)
              end do
            end do
          end do
        end do
        call rot_tensor_6( Mat6, Rot_tem )
        jhe = 0
        do je = 1,3
          do he = 1,3
            jhe = jhe + 1
            jhs = 0
            do js = 1,3
              do hs = 1,3
                jhs = jhs + 1
                secoo_m(:,jhe,:,jhs,initlr,mpirank) = Mat6(:,je,he,:,js,hs)
              end do
            end do
          end do
        end do
      endif

    end do

  endif

  return
  100 format(/' ---- Radial ---------',100('-'))
  110 format(/' ---- Tens_ab --------',100('-'))
  120 format(/' lseuil =',i2)
  130 format(/' -- irang =',i2,', jrang =',i2,' --')
  140 format(/' ke, je, he =',3i3,',  ks, js, hs =',3i3)
  150 format(/' le, me =',2i3,',  ls, ms =',2i3,', Clme, Clms =',1p, 2e13.5)
  160 format(/' Tensor by harmonics (basis R4) :'/, ' ir jr  le me  ls ms      Tens(Ylm,i=1,ninitlr)')
  170 format(2i3,2(i4,i3),1p,20e15.7)
  180 format('i',i2,i3,2(i4,i3),1p,20e15.7)

end

!***********************************************************************

subroutine tens_ab(coef_g,Core_resolved, Final_tddft,Green_int,icheck,ip_max,ip0,irang,is_g,jrang, &
                              le,me,ls,m_g,ms,lmax,lmoins1,lplus1,lseuil,ndim1,ndim2,ninit1, &
                              ninitl,ninitlv,ninitlr,nlma,nlma2,nspinp,nspino,rof,Singul,Solsing, &
                              Solsing_only,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)

  use declarations
  implicit none

  integer, intent(in):: icheck, ip_max, ip0, irang, jrang, le, me, ls, ms, lseuil, ndim1, ndim2, ninit1, &
      ninitl, ninitlv, ninitlr, nlma, nlma2, nspinp, nspino
  integer, dimension(ninitl,2), intent(in):: m_g
  integer, dimension(ninitl), intent(in):: is_g

  integer:: i_g_1, i_g_2, initl1, initl2, initlr, is_r1, is_r2, iseuil1, iseuil2, iso1, iso2, ispf1, ispf2, ispg12, ispg21, &
    ispinf1, ispinf2, isping1, isping2, ispp1, ispp2, l1, l2, li, lm01, lm02, lm1, lm2, lmax, lmp01, lmp02, lmp1, lmp2, &
    lms1, lms2, lp1, lp2, m1, m2, mi1, mi2, mp1, mp2, mv1, mv2

  complex(kind=db) :: cfe, cfs, Cg, dfe, dfs, Tau_rad, Tau_rad_i

  complex(kind=db):: Gaunte, Gauntm, Gauntmag, Gaunts
  complex(kind=db), dimension(ninitlr):: Ten, Ten_m
  complex(kind=db), dimension(nlma,nlma2,nspinp,nspino,ip0:ip_max,ninitlv):: rof
  complex(kind=db), dimension(nlma*nspino,nlma*nspino,2,2,ndim1,ndim2,ndim2):: Taull
  complex(kind=db), dimension(nlma,nspinp,ip0:ip_max,ip0:ip_max,ninitlv):: Singul

  logical:: Core_resolved, Final_tddft, Green_int, lmoins1, lplus1, Solsing, Solsing_only, Spinorbite, Titre, Ylm_comp

  real(kind=db):: cfi, cfr, Ci_1, Ci_2, Ci2, J_initl1, J_initl2, Jz1, Jz2
  real(kind=db), dimension(ninitl,2):: coef_g

  li = lseuil

  Ten(:) = (0._db,0._db)
  Ten_m(:) = (0._db,0._db)

! Boucle sur les etats initiaux
  do i_g_1 = 1,ninitl

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
!        if(Final_tddft ) then
!          write(6,*) J_initl1, li, i_g_1, is_g(i_g_1)
!        endif

    if( i_g_1 <= ninit1 ) then
      iseuil1 = 1
    else
      iseuil1 = 2
    endif

    if( Final_tddft ) then
      initlr = 1
    elseif( Core_resolved ) then
      initlr = i_g_1
    else
      initlr = iseuil1
    endif

    if( ninitlv == ninitl .and. Final_tddft ) then
      is_r1 = i_g_1
    else
      is_r1 = iseuil1
    endif

! Boucle sur le spin
    do isping1 = 1,2  ! Spin de l'etat initial

      mi1 = m_g(i_g_1,isping1)
      Ci_1 = Coef_g(i_g_1,isping1)

      if( abs( Ci_1 ) < eps6 ) cycle

      do i_g_2 = 1,ninitl

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

        if( ninitlv == ninitl .and. Final_tddft ) then
          is_r2 = i_g_2
        else
          is_r2 = iseuil2
        endif

        do isping2 = 1,2  ! Spin de l'etat initial
!        if( isping1 /= isping2 ) cycle
          mi2 = m_g(i_g_2,isping2)
          Ci_2 = Coef_g(i_g_2,isping2)
          if( ndim1 > 1 ) then
            ispg12 = isping2 + 2 * ( isping1 - 1 )
            ispg21 = isping1 + 2 * ( isping2 - 1 )
          else
            ispg12 = 1
            ispg21 = 1
          endif

          if( abs( Ci_2 ) < eps6 ) cycle

          Ci2 = Ci_1 * Ci_2

! Boucle sur les harmoniques de l'etat final en entree
          do l1 = 0,lmax

            if( lplus1 .and. l1 < li + 1 ) cycle
            if( lmoins1 .and. l1 > li ) cycle
            if( l1 > li + irang .or. l1 < li - irang .or. mod(l1,2) /= mod(li+irang,2) ) cycle

            lm01 = l1**2 + l1 + 1
            do m1 = -l1,l1
              lm1 = lm01 + m1
              if( lm1 > nlma ) cycle

              do ispinf1 = 1,2  ! spin de l'etat final en entree
                ispf1 = min( ispinf1, nspinp )

                if( irang == 0 ) then
                  Gaunte = Gauntmag(isping1,ispinf1,me,l1,m1,li,mi1,Ylm_comp,.true.)
                else
                  Gaunte = Gauntm(l1,m1,le,me,li,mi1,Ylm_comp)
                endif
                if( abs(Gaunte) < eps10 ) cycle

! Il n'y a que pour le dipole magnetique qu'on peut avoir du spin-flip
                if( ispinf1 /= isping1 .and. irang /= 0 ) cycle

! Boucle sur les harmoniques de l'etat final en sortie
                do l2 = 0,lmax
                  if( lplus1 .and. l2 < li + 1 ) cycle
                  if( lmoins1 .and. l2 > li ) cycle
                  if( l2 > li + jrang .or. l2 < li - jrang .or. mod(l2,2) /= mod(li+jrang,2) ) cycle

                  lm02 = l2**2 + l2 + 1
                  do m2 = -l2,l2
                    lm2 = lm02 + m2
                    if( lm2 > nlma ) cycle

                    do ispinf2 = 1,2  ! spin etat final en sortie
                      ispf2 = min( ispinf2, nspinp )

                      if( ispinf2 /= isping2 .and. jrang /= 0 ) cycle
!                          if( .not. Spinorbite .and. ispinf2 /= ispinf1 ) cycle

                      if( (.not. Final_tddft) .and. ispinf2 /= ispinf1 ) cycle

                      if( jrang == 0 ) then
                        Gaunts = Gauntmag(isping2,ispinf2,ms,l2,m2,li,mi2,Ylm_comp,.true.)
                      else
                        Gaunts = Gauntm(l2,m2,ls,ms,li,mi2,Ylm_comp)
                      endif
                      if( abs(Gaunts) < eps10 ) cycle

                      Cg = - Ci2 * conjg( Gaunte ) * Gaunts

                      cfe = (0._db, 0._db)
                      cfs = (0._db, 0._db)
                      Titre = .true.

        if( .not. Solsing_only ) then

! Boucle sur les harmoniques de l'etat final en entree potentiel non spherique
          lmp1 = 0
          do lp1 = 0,lmax
            if( nlma2 == 1 .and. lp1 /= l1 ) cycle
            lmp01 = lp1**2 + lp1 + 1

            do mp1 = -lp1,lp1
              if( nlma2 == 1 .and. mp1 /= m1 ) cycle
              lmp1 = lmp1 + 1

! Boucle sur les solutions en entree
              do ispp1 = 1,nspinp
                if( .not. Spinorbite .and. ispp1 /= ispf1 ) cycle
                iso1 = min( ispp1, nspino )

                if( Spinorbite .and. nlma2 == 1 ) then
                  mv1 = mp1 - ispinf1 + iso1
                  if( mv1 > lp1 .or. mv1 < -lp1 ) cycle
                else
                  mv1 = mp1
                endif
                lms1 = nlma * ( iso1 - 1 ) + lmp01 + mv1

! Boucle sur les harmoniques de l'etat final en entree potentiel non spherique
                lmp2 = 0
                do lp2 = 0,lmax
                  if( nlma2 == 1 .and. lp2 /= l2 ) cycle
                  lmp02 = lp2**2 + lp2 + 1

                  do mp2 = -lp2,lp2
                    if( nlma2 == 1 .and. mp2 /= m2 ) cycle
                    lmp2 = lmp2 + 1

                    if( nlma2 == 1 ) then
                      lmp2 = 1
                    else
                      lmp2 = lmp02 + mp2
                    endif

                    do ispp2 = 1,nspinp
                      if( .not. Spinorbite .and. ispp2 /= ispf2 ) cycle
                      iso2 = min( ispp2, nspino )

                      if( Spinorbite .and. nlma2 == 1 ) then
                        mv2 = mp2 - ispinf2 + iso2
                        if( mv2 > lp2 .or. mv2 < -lp2 ) cycle
                      else
                        mv2 = mp2
                      endif
                      lms2 = nlma * (iso2 - 1) + lmp02 + mv2

                      dfe = rof(lm1,lmp1,ispf1,iso1,irang,is_r1) * rof(lm2,lmp2,ispf2,iso2,jrang,is_r2) &
                          * Taull(lms2,lms1,ispinf2,ispinf1,ispg12,initl2,initl1)
                      dfs = rof(lm2,lmp2,ispf2,iso2,jrang,is_r2) * rof(lm1,lmp1,ispf1,iso1,irang,is_r1) &
                          * Taull(lms1,lms2,ispinf1,ispinf2,ispg21,initl1,initl2)

                      cfe = cfe + dfe
                      cfs = cfs + dfs

                      if( icheck > 2 .and. ( abs( dfe ) > eps15 .or. abs( dfs ) > eps15) ) then
                        if( Titre ) then
                          Titre = .false.
                          write(3,120)
                          write(3,130) i_g_1, isping1, Jz1, mi1, Ci_1, i_g_2, isping2, Jz2, mi2, Ci_2, l1, &
                          m1, ispinf1, l2, m2, ispinf2
                          write(3,132)
                        endif
                        cfr = 0.5_db * real( dfe - dfs, db )
                        cfi = 0.5_db * aimag( dfe + dfs )
                        write(3,135) lp1, mp1, iso1, lp2, mp2, iso2, cfr, cfi, Gaunts, Gaunte, &
                          rof(lm1,lmp1,ispf1,iso1,irang,is_r1), rof(lm2,lmp2,ispf2,iso2,jrang,is_r2), &
                          Taull(lms2,lms1,ispinf2,ispinf1,ispg12,initl2,initl1), Taull(lms1,lms2,ispinf1,ispinf2, &
                                    ispg21,initl1,initl2)
                      endif

                    end do ! fin boucle sur les solutions en sortie
                  end do ! fin boucle mp2 sortie pot non spherique
                end do ! fin boucle lp2 sortie pot non spherique

              end do ! fin boucle isol1 solutions entree
            end do ! fin boucle mp1 entree pot non spherique
          end do ! fin boucle lp1 entree pot non spherique

        endif

! Comme on ne divise pas par pi, le resultat apparait comme multiplie
! par pi, si on calcule ensuite le facteur de structure. La
! normalisation par pi n'est donc pas a faire dans coabs sur l'amplitude
! dafs.
                      if( Solsing .and. i_g_1 == i_g_2 .and. lm1 == lm2 .and. ispinf1 == ispinf2 ) then
                        cfe = cfe + Singul(lm1,ispf1,irang,jrang,is_r1)
                        cfs = cfs + Singul(lm1,ispf1,irang,jrang,is_r1)
                      endif

                      if( Green_int ) then
                        Tau_rad = 0.5_db * ( cfe + cfs )
                        Tau_rad_i = 0.5_db * (cfe - cfs)
                        Ten(initlr) = Ten(initlr) + Real(Cg,db) * tau_rad + img * aimag(Cg)* tau_rad_i
                        Ten_m(initlr) = Ten_m(initlr) + Real(Cg,db) * tau_rad_i + img * aimag(Cg) * tau_rad
                      else
                        cfr = real( cfe - cfs, db )
                        cfi = aimag( cfe + cfs )
                        Tau_rad = 0.5_db * cmplx(cfr, cfi,db)
                        Ten(initlr) = Ten(initlr) + Cg * Tau_rad
                      endif

                      if( icheck > 2 .and.  abs( Tau_rad ) > eps15 ) write(3,140) Ten(initlr), Tau_rad, Cg, &
                                 Singul(lm1,ispf1,irang,jrang,is_r1)

                    end do ! fin boucle ispinf2, spin etat final sortie
                  end do ! fin boucle m2 etat final sortie
                end do ! fin boucle l2 etat final sortie

              end do  ! fin boucle sur spin d'etat final en entree
            end do  ! fin boucle m1 etat final entree
          end do ! fin boucle l1 etat final entree

        end do    ! fin boucle sur le spin d'etat initial sortie
      end do    ! fin boucle sur les etats initiaux sortie
    end do    ! fin boucle sur le spin d'etat initial entree
  end do   ! fin boucle sur les etats initiaux entree

! On multiplie par -i pour que ce soit la partie reelle du tenseur qui soit l'absorption.
  Ten(:) = - img * Ten(:)
  Ten_m(:) = - img * Ten_m(:)

  return
  120 format(/' ini1 isg1 Jz1 mi1  Ci1  ini2 isg2 Jz2 mi2  Ci2   l1  m1 isf1  l2  m2 isf2')
  130 format(2(2i4,f6.1,i3,f7.3),2(3i4,2x))
  132 format(7x,'lp1 mp1 is1 lp2 mp2 is2',12x,'Tau_rad',24x, 'Gaunts',25x,'Gaunte',26x,'rofe',26x,'rofs',27x,'Taue',27x,'Taus')
  135 format(5x,6i4,1p,2e15.7,2(1x,2e15.7),1x,4e15.7,8(1x,2e15.7))
  140 format(17x,'Ten',26x,'Tau_rad',26x,'Cg',27x,'Singul'/,1p,1x,4(1x,2e15.7))
end

!***********************************************************************

subroutine tens_op(Core_resolved,icheck,ip_max,ip0,irang,jrang,le,me,ls,ms,lmax,lmoins1,lplus1,M_depend, &
                                ndim1,ndim2,ninitlr,nlm1g,nlm2g,nlma,nlma2,nspinp,nspino,roff_ii,roff_ir,roff_ri,roff_rr, &
                                Spinorbite,Taull,Ten,Ten_m,Ylm_comp)

  use declarations
  implicit none

  integer, intent(in):: icheck, ip_max, ip0, irang, jrang, le, me, &
    ls, ms, ndim1, ndim2, nlm1g, nlm2g, nlma, nlma2, nspinp, nspino

  integer:: i, io1, io2, is1, is2, iso1, iso2, isoi1, isoi2, isp1, &
    isp2, ispi1, ispi2, ispin1, ispin2, ispini1, ispini2, ispip1, &
    ispip2, ispp1, ispp2, l1, l2, li1, li2, lip1, lip2, lm1, lm2, lmax, lmi1, lmi2, lmip1, lmip2, lmis1, lmis2, lmp1, lmp2, &
    lms1, lms2, lp1, lp2, m1, m2, mi1, mi2, mip1, mip2, miv1, miv2, mp1, mp2, mv1, mv2, ninitlr

  complex(kind=db):: cf, Gaunte, Gaunt_crc, Gauntmag, Gaunts, Sum_ff, Sum_ii, Tau_Tau_rad
  complex(kind=db), dimension(ninitlr):: Ten, Ten_m
  complex(kind=db), dimension(nlma*nspino,nlma*nspino,2,2,ndim1,ndim2,ndim2):: Taull

  logical:: Core_resolved, lmoins1, lplus1, M_depend, Spinorbite, Titre, Ylm_comp

  real(kind=db):: cfi, cfr, Ro
  real(kind=db), dimension(nlm1g,nlm1g,nlm2g,nlm2g,nspinp**2,nspino**2,ip0:ip_max):: roff_ii, roff_ir, roff_ri, roff_rr

  Titre = .true.

  Ten(:) = (0._db,0._db)
  Ten_m(:) = (0._db,0._db)

! Boucle harmoniques etat initiaux en entree
  do li1 = 0,lmax

    do mi1 = -li1,li1
      if( M_depend ) then
        lmi1 = li1**2 + li1 + 1 + mi1
      else
        lmi1 = li1 + 1
      endif

! Boucle spin etat initial, en entree
      do ispini1 = 1,2
        ispi1 = min( ispini1, nspinp )

! Boucle harmoniques etat initial, en sortie
        do li2 = 0,lmax

          do mi2 = -li2,li2
            if( M_depend ) then
              lmi2 = li2**2 + li2 + 1 + mi2
            else
              lmi2 = li2 + 1
            endif

! Boucle spin etat initial, en sortie
            do ispini2 = 1,2
              ispi2 = min( ispini2, nspinp )

              if( ispini2 /= ispini1 ) cycle

! Boucle harmoniques etat final en entree
              do l1 = 0,lmax

                if( lplus1 .and. l1 < li1 + 1 ) cycle
                if( lmoins1 .and. l1 > li1 ) cycle
                if( l1 > li1 + irang .or. l1 < li1 - irang .or. mod(l1,2) /= mod(li1+irang,2) ) cycle

                do m1 = -l1,l1
                  if( M_depend ) then
                    lm1 = l1**2 + l1 + 1 + m1
                  else
                    lm1 = l1 + 1
                  endif

! Boucle spin etat final, en entree
                  do ispin1 = 1,2

! Il n'y a que pour le dipole magnetique qu'on peut avoir du spin-flip
                    if( ispin1 /= ispini1 .and. irang /= 0 ) cycle

                    isp1 = min( ispin1, nspinp )

                    if( irang == 0 ) then
                      Gaunte = Gauntmag(ispin1,ispini1,me,li1,mi1,l1,m1,Ylm_comp,Ylm_comp)
                    else
                      Gaunte = Gaunt_crc(li1,mi1,le,me,l1,m1,Ylm_comp)
                    endif
                    if( abs(Gaunte) < eps10 ) cycle

! Boucle harmoniques etat final, en sortie
                    do l2 = 0,lmax
                      if( lplus1 .and. l2 < li2 + 1 ) cycle
                      if( lmoins1 .and. l2 > li2 ) cycle
                      if( l2 > li2 + jrang .or. l2 < li2 - jrang .or. mod(l2,2) /= mod(li2+jrang,2) ) cycle

                      do m2 = -l2,l2
                        if( M_depend ) then
                          lm2 = l2**2 + l2 + 1 + m2
                        else
                          lm2 = l2 + 1
                        endif

                        do ispin2 = 1,2  ! spin etat final en sortie
                          isp2 = min( ispin2, nspinp )

                          if( ispin2 /= ispini2 .and. jrang /= 0 ) cycle
!                          if( .not. Spinorbite .and. ispin2 /= ispin1 ) cycle

                          if( ispin2 /= ispin1 ) cycle

                          if( jrang == 0 ) then
                            Gaunts = Gauntmag(ispini2,ispin2,ms,l2,m2,li2,mi2,Ylm_comp,Ylm_comp)
                          else
                            Gaunts = Gaunt_crc(l2,m2,ls,ms,li2,mi2,Ylm_comp)
                          endif
                          if( abs(Gaunts) < eps10 ) cycle

                          Tau_Tau_rad = 0._db

                          if( Spinorbite .or. nlma2 > 1 ) Titre = .true.

! Boucle harmoniques pot non spherique, etat initial, en entree
    do lip1 = 0,lmax
      if( nlma2 == 1 .and. lip1 /= li1 ) cycle

      do mip1 = -lip1,lip1
        if( nlma2 == 1 .and. mip1 /= mi1 ) cycle
        lmip1 = min( lip1**2 + lip1 + 1 + mip1, nlma2 )

! Boucle solutions etat initial en entree
        do ispip1 = 1,nspinp
          if( .not. Spinorbite .and. ispip1 /= ispi1 ) cycle
          isoi1 = min( ispip1, nspino )

          if( Spinorbite .and. nlma2 == 1 ) then
            miv1 = mip1 - ispini1 + isoi1
            if( miv1 > lip1 .or. miv1 < -lip1 ) cycle
          else
            miv1 = mip1
          endif
          lmis1 = nlma * ( isoi1 - 1 ) + lip1**2 + lip1 + 1 + miv1

! Boucle harmoniques pot non spherique, etat initial, en sortie
          do lip2 = 0,lmax
            if( nlma2 == 1 .and. lip2 /= li2 ) cycle

            do mip2 = -lip2,lip2
              if( nlma2 == 1 .and. mip2 /= mi2 ) cycle
              lmip2 = min( lip2**2 + lip2 + 1 + mip2, nlma2 )

! Boucle solutions, etat initial, en sortie
              do ispip2 = 1,nspinp
                if( .not. Spinorbite .and. ispip2 /= ispi2 ) cycle
                isoi2 = min( ispip2, nspino )

                if( Spinorbite .and. nlma2 == 1 ) then
                  miv2 = mip2 - ispini2 + isoi2
                  if( miv2 > lip2 .or. miv2 < -lip2 ) cycle
                else
                  miv2 = mip2
                endif
                lmis2 = nlma * (isoi2 - 1) + lip2**2 + lip2 + 1+miv2

! Boucle harmoniques pot non spherique, etat final, en entree
                do lp1 = 0,lmax
                  if( nlma2 == 1 .and. lp1 /= l1 ) cycle

                  do mp1 = -lp1,lp1
                    if( nlma2 == 1 .and. mp1 /= m1 ) cycle

                    lmp1 = min( lp1**2 + lp1 + 1 + mp1, nlma2 )

! Boucle solutions, etat final, en entree
                    do ispp1 = 1,nspinp
                      if( .not. Spinorbite.and. ispp1 /= isp1) cycle
                      iso1 = min( ispp1, nspino )

                      if( Spinorbite .and. nlma2 == 1 ) then
                        mv1 = mp1 - ispin1 + iso1
                        if( mv1 > lp1 .or. mv1 < -lp1 ) cycle
                      else
                        mv1 = mp1
                      endif
                      lms1 = nlma * (iso1-1) + lp1**2 + lp1 + 1 +mv1

! Boucle harmoniques pot non spherique, etat final, en sortie
                      do lp2 = 0,lmax
                        if( nlma2 == 1 .and. lp2 /= l2 ) cycle

                        do mp2 = -lp2,lp2
                          if(nlma2 == 1 .and. mp2 /= m2 ) cycle

                          lmp2 = min( lp2**2 + lp2 + 1 + mp2, nlma2)

! Boucle solutions, etat final, en sortie
                          do ispp2 = 1,nspinp
                            if( .not. Spinorbite .and. ispp2 /= isp2 ) cycle
                            iso2 = min( ispp2, nspino )

                            if( Spinorbite .and. nlma2 == 1 ) then
                              mv2 = mp2 - ispin2 + iso2
                              if( mv2 > lp2 .or. mv2 < -lp2 ) cycle
                            else
                              mv2 = mp2
                            endif
                            lms2 = nlma * (iso2 - 1) + lp2**2 + lp2 + 1 + mv2

                            is1 = isp1 + (ispi1 - 1) * nspinp
                            is2 = isp2 + (ispi2 - 1) * nspinp
                            io1 = iso1 + (isoi1 - 1) * nspino
                            io2 = iso2 + (isoi2 - 1) * nspino

                            Ro = roff_rr(lmi1,lm1,lmip1,lmp1,is1,io1,irang) * roff_rr(lmi2,lm2,lmip2,lmp2,is2,io2,jrang)
                            
                            cfr = real( Taull(lms1,lms2,ispin1,ispin2,2,1,1) - Taull(lms2,lms1,ispin2,ispin1,2,1,1), db )  
                            cfi = aimag( Taull(lms1,lms2,ispin1,ispin2,2,1,1) + Taull(lms2,lms1,ispin2,ispin1,2,1,1) )
                            Sum_ff = 0.5_db * cmplx( cfr, cfi, db ) 

                            cfr = real( Taull(lmis2,lmis1,ispini2,ispini1,1,1,1) - Taull(lmis1,lmis2,ispini1,ispini2,1,1,1), db )  
                            cfi = aimag( Taull(lmis2,lmis1,ispini2,ispini1,1,1,1) + Taull(lmis1,lmis2,ispini1,ispini2,1,1,1) )
                            Sum_ii = 0.5_db * cmplx( cfr, cfi, db ) 
                            Tau_Tau_rad = Tau_Tau_rad - Ro * Sum_ff * Sum_ii

                    !        cfr = real( Taull(lms1,lms2,ispin1,ispin2,2,1,1) * Taull(lmis2,lmis1,ispini2,ispini1,1,1,1) &
                    !               + conjg( Taull(lms2,lms1,ispin2,ispin1,2,1,1) * Taull(lmis1,lmis2,ispini1,ispini2,1,1,1) ), db )  

                    !        cfi = aimag( Taull(lms1,lms2,ispin1,ispin2,2,1,1) * Taull(lmis2,lmis1,ispini2,ispini1,1,1,1) &
                    !                   - conjg( Taull(lms2,lms1,ispin2,ispin1,2,1,1) * Taull(lmis1,lmis2,ispini1,ispini2,1,1,1) ) )  

                    !        Tau_Tau_rad = Tau_Tau_rad - Ro * 0.5_db * cmplx( cfr, cfi, db )
                             
                            if( icheck < 2 .or. ( .not. Spinorbite .and. nlma2 == 1 ) .or. &
                                                abs(Ro * Sum_ff * Sum_ii) < 1.e-15_db ) cycle

                            if( nlma2 == 1 ) then
                              if( Titre ) then
                                Titre = .false.
                                write(3,110)
                                write(3,120) li1, mi1, ispini1, l1, m1, ispin1, l2, m2, ispin2, li2, mi2, &
                                ispini2, isoi1, iso1, iso2, isoi2, roff_rr(lmi1,lm1,lmip1,lmp1,is1,io1,irang), &
                                roff_rr(lmi2,lm2,lmip2,lmp2,is2,io2,jrang), Taull(lmis2,lmis1,ispini2,ispini1,1,1,1), &
                                Taull(lms1,lms2,ispin1,ispin2,2,1,1), Tau_Tau_rad
                              else
                                write(3,130) isoi1, iso1, iso2, isoi2, roff_rr(lmi1,lm1,lmip1,lmp1,is1,io1,irang), &
                                roff_rr(lmi2,lm2,lmip2,lmp2,is2,io2,jrang), Taull(lmis2,lmis1,ispini2,ispini1,1,1,1), &
                                Taull(lms1,lms2,ispin1,ispin2,2,1,1), Tau_Tau_rad
                              endif
                            else
                              if( Titre ) then
                                Titre = .false.
                                write(3,140)
                                write(3,150) li1, mi1, ispini1, l1, m1, ispin1, l2, m2, ispin2, li2, mi2, ispini2
                                write(3,160)
                              endif
                              write(3,170) lip1, mip1, isoi1, lp1, mp1, iso1, lp2, mp2, iso2,  lip2, mip2, isoi2, &
                                roff_rr(lmi1,lm1,lmip1,lmp1,is1,io1,irang), roff_rr(lmi2,lm2,lmip2,lmp2,is2,io2,jrang), &
                                Taull(lmis2,lmis1,ispini2,ispini1,1,1,1), Taull(lms1,lms2,ispin1,ispin2,2,1,1), Tau_Tau_rad
                            endif

                          end do ! fin boucle sur les solutions en sortie
                        end do ! fin boucle mp2 sortie pot non spherique
                      end do ! fin boucle lp2 sortie pot non spherique

                    end do ! fin boucle isol1 solutions entree
                  end do ! fin boucle mp1 entree pot non spherique
                end do ! fin boucle lp1 entree pot non spherique

              end do ! fin boucle sur les solutions etat initial, en sortie
            end do ! fin boucle mip2 sortie pot non spherique etat initial,
          end do ! fin boucle lip2 sortie pot non spherique etat initial

        end do ! fin boucle isoli1 solutions etat initial, entree
      end do ! fin boucle mip1 entree pot non spherique etat initial,
    end do ! fin boucle lip1 entree pot non spherique etat initial

! Bien que Tau apparaisse 2 fois, on divise une seule fois par pi.
! La normalisation par pi n'est donc pas a faire dans coabs sur l'amplitude dafs.

                          cf = Tau_Tau_rad * Gaunte * Gaunts / pi

                          Ten(1) = Ten(1) + cf
                          if( Core_resolved ) then
                            do i = 2,ninitlr
                              if( li1 == i-2 .and. li2 == i-2 ) Ten(i) = Ten(i) + cf
                            end do
                          endif

                          if( icheck < 2 ) cycle
                          if( Spinorbite .or. nlma2 > 1 ) then
                            if( Titre ) cycle
                            write(3,180) Ten(1), cf, Gaunte, Gaunts
                          elseif( abs(cf) > 1.e-15_db ) then
                            if( Titre ) write(3,190)
                            Titre = .false.
                            write(3,200) li1, mi1, ispini1, l1, m1, ispin1, l2, m2, ispin2, li2, mi2, &
                              ispini2, Ten(1), cf, Gaunte, Gaunts, &
                              roff_rr(lmi1,lm1,lmip1,lmp1,is1,io1,irang), roff_rr(lmi2,lm2,lmip2,lmp2,is2,io2,jrang), &
                              Taull(lmis1,lmis2,ispini1,ispini2,1,1,1), Taull(lms2,lms1,ispin2,ispin1,2,1,1)
                          endif

                        end do ! fin boucle ispin2, spin etat final sortie
                      end do ! fin boucle m2 etat final sortie
                    end do ! fin boucle l2 etat final sortie

                  end do  ! fin boucle sur spin d'etat final en entree
                end do  ! fin boucle m1 etat final entree
              end do ! fin boucle l1 etat final entree

            end do    ! fin boucle spin etat initial, sortie
          end do    ! fin boucle mi2, etat initial, sortie
        end do    ! fin boucle li2, etat initial, sortie

      end do    ! fin boucle spin, etat initial, entree
    end do    ! fin boucle mi1, etat initial, entree
  end do   ! fin boucle li1, etat initial, entree

  return
  110 format(/'  li1 mi1 si1  l1  m1 sp1  l2  m2 sp2 li2 mi2 si2 soi1 so1 so2 soi2    roff1',8x,'roff2',14x,'Tau_1',21x,'Tau_2', &
              18x,'Tau_Tau_rad')
  120 format(12i4,1x,4i4,1x,1p,28e13.5)
  130 format(49x,4i4,1x,1p,28e13.5)
  140 format(/'  li1 mi1 si1  l1  m1 sp1  l2  m2 sp2 li2 mi2 si2')
  150 format(12i4)
  160 format(3x,'lip1 mip1 isoi1 lp1 mp1 iso1  lp2  mp2 iso2 lip2 mip2 isoi2   roff1',8x,'roff2',14x,'Tau_1',21x,'Tau_2',18x, &
                'Tau_Tau_rad')
  170 format(12i5,1p,28e13.5)
  180 format(5x,'Ten =',1p,2e13.5,', dTen =',2e13.5,', Gaunte =',2e13.5,', Gaunts =',2e13.5)
  190 format(/'  li1 mi1 si1  l1  m1 sp1  l2  m2 sp2 li2 mi2 si2',11x,'Ten',22x,'dTen',21x,'Gaunte',20x,'Gaunts',14x, &
              'roff1',8x,'roff2',15x,'Tau_1',21x,'Tau_2')
  200 format(12i4,1p,28e13.5)
end

!***********************************************************************

! Calcule le coef de Gaunt avec harmoniques complexes = Int( Y(l1,m1)*Y(l2,m2)Y(l3,m3)dOmega )
! Formule de M. E. Rose p. 62

function Gauntcp(l1,m1,l2,m2,l3,m3)

  use declarations
  
  implicit none
  
  integer:: l1, l2, l3, m1, m2, m3
  
  real(kind=db):: cgc, fac, Gauntcp

  if( m1 /= m2 + m3 .or. mod(l1+l2+l3,2) == 1 .or. l1 < abs(l2-l3) .or. l1 > abs(l2+l3) ) then

    Gauntcp = 0._db

  else

    fac = sqrt( (2*l3 + 1._db) * (2*l2 + 1._db) / ( (2*l1 + 1._db) * quatre_pi ) )
    Gauntcp = fac * cgc(l3,l2,l1,m3,m2) * cgc(l3,l2,l1,0,0)

  endif

  return
end

!***********************************************************************

! Clebsch-gordan coefficient  eq. 3.18, Rose

function cgc(l1,l2,l3,m1,m2)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(5):: nd, num

  cgc = 0._db

  m3 = m1 + m2

! Arguments of factorials
  num( 1 ) = l3 + l1 - l2
  num( 2 ) = l3 - l1 + l2
  num( 3 ) = l1 + l2 - l3
  num( 4 ) = l3 + m3
  num( 5 ) = l3 - m3
  nd( 1 ) = l1 + l2 + l3 + 1
  nd( 2 ) = l1 - m1
  nd( 3 ) = l1 + m1
  nd( 4 ) = l2 - m2
  nd( 5 ) = l2 + m2

! Check triangle and projection conditions
  do i = 1, 5
    if( num(i) < 0 .or. nd(i) < 0 ) return
  end do

  nff = 0
  nmin = max( 0, l2 + m3 - l1 )
  ff = 1.0_db

! Two sets of factorial products
  n = 5
  do nfac = 1, 2
    n1 = n - 1
! Arrange arguments in descending order
    do i = 1, n1
      inum = i
      id = i
      i1 = i + 1
      do j = i1, n
        if( num( j ) > num( inum ) ) inum = j
        if( nd( j ) > nd( id ) ) id = j
      end do
      ntemp = num( i )
      num( i ) = num( inum )
      num( inum ) = ntemp
      ntemp = nd( i )
      nd( i ) = nd( id )
      nd( id ) = ntemp
    end do
! Compute factorial ratios
    do i = 1, n
      if( num( i ) < nd( i ) ) then
        jm = nd( i )
        if( jm == 1 ) cycle
        j0 = num( i ) + 1
        if( num( i ) == 0 ) j0 = 2
        do j = j0, jm
          if( abs ( ff ) < 1.0e-20_db ) then
            ff = ff * 1.e20_db
            nff = nff - 2
          endif
          ff = ff / j
        end do
      elseif( num( i ) > nd( i ) ) then
        jm = num( i )
        if( jm == 1 ) cycle
        j0 = nd( i ) + 1
        if( nd( i ) == 0 ) j0 = 2
        do j = j0, jm
          if( abs ( ff ) > 1.0e20_db ) then
            ff = ff / 1.0e20_db
            nff = nff + 2
          endif
          ff = j * ff
        end do
      endif
    end do

    if ( nfac == 2 ) exit

    nff = nff / 2
    ff = sqrt( ( 2 * l3 + 1 ) * ff )

! Second set of factorial arguments
    num( 1 ) = l2 + l3 + m1 - nmin
    num( 2 ) = l1 - m1 + nmin
    num( 3 ) = 0
    nd( 1 ) = nmin
    if( nmin == 0 ) nd( 1 ) = l1 - l2 - m3
    nd( 2 ) = l3 - l1 + l2 - nmin
    nd( 3 ) = l3 + m3 - nmin
    n = 3
  end do

  if( mod( nmin + l2 + m2, 2 ) /= 0 ) ff = - ff
  ff = ff * 1.0e10_db**nff
  cgcp = ff
  nmax = min( l3 - l1 + l2, l3 + m3 )

  do nu = nmin+1, nmax
    ff= - ff * ( (l1-m1+nu) * (l3-l1+l2-nu+1) * (l3+m3-nu+1) ) / real( nu * (nu+l1-l2-m3) * (l2+l3+m1-nu+1), db )
    cgcp = cgcp + ff
  end do

  cgc = cgcp

  return
end

!***********************************************************************

! This subroutine calculates the Gaunt coefficients for real spherical harmonics

function Gauntc( l1, m1, l2, m2, l3, m3 )

  use declarations
  implicit none

  integer:: itr, ityfct1, ityfct2, ityfct3, itygam, ityp, itypi, l1, l2, l3, ltsti, ltst, m, m1, m2, m3, mgam, mi, nexp, ns, nsf
 
  real(kind=db):: ab, cgc, pre, Gauntc, gf, sqr2i
  
  if( mod( l3+l1+l2, 2 ) /= 0 .or. abs( l1 - l3 ) > l2 .or. ( l1 + l3 ) < l2 ) then
    Gauntc = 0._db
    return
  endif

  mi = abs( m1 )
  mgam = abs( m2 )
  m = abs( m3 )

  if( mi+mgam /= m .and. mi+m /= mgam .and. m+mgam /= mi ) then
    gauntc = 0._db
    return
  endif

  sqr2i = 1 / sqrt( 2._db )

  ltsti = l1
  itr = l2
  ltst = l3

  pre = sqrt( ( 2*ltst+1 ) * ( 2*itr+1 ) / ( ( 2*ltsti +1 ) * quatre_pi ) ) * cgc(ltst,itr,ltsti,0,0)

  if( m1 == 0 ) then
    itypi = 0
  else
    itypi = mi / m1
  endif
  if( m2 == 0 ) then
    itygam = 0
  else
    itygam = mgam / m2
  endif
  if( m3 == 0 ) then
    ityp = 0
  else
    ityp = m / m3
  endif

  nexp = abs(itypi) + abs(itygam) + abs(ityp)
  ns = itypi + ityp + itygam
  if( ns == 0 ) ns = 1
  nsf = isign(1,ns)
  gf = nsf * (-1)**mi * sqr2i**nexp * pre

  if( m+mgam+mi == 0 ) then
    ab = cgc(ltst,itr,ltsti,0,0)
  else
    ityfct1 = itypi + itygam * ityp
    ityfct2 = itygam + itypi * ityp
    ityfct3 = ityp + itypi * itygam
    ab = 0._db
    if( m+mgam == mi .and. ityfct1 /= 0 ) ab = (-1)**mi * cgc(ltst,itr,ltsti,m,mgam) * ityfct1
    if( m+mi == mgam .and. ityfct2 /= 0 ) ab = ab + (-1)**mgam * cgc(ltst,itr,ltsti,-m,mgam) * ityfct2
    if( mgam+mi == m .and. ityfct3 /= 0 ) ab = ab + (-1)**m * cgc(ltst,itr,ltsti,m,-mgam) * ityfct3
  endif

  gauntc = gf * ab

  return
end

!***********************************************************************

! Calcule le coefficient de Gaunt avec Y(l,m) complexe ou reel, Y(lo,mo) reel et Y(li,mi) complexe. 

function Gauntm(l,m,lo,mo,li,mi,Ylmcomp)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  complex(kind=db) Gauntm

  logical Ylmcomp

  if( Ylmcomp ) then

    if( mo == 0 ) then
      gr = gauntcp(l,m,lo,mo,li,mi)
      gi = 0._db
    elseif( mo > 0 ) then
      gr = (   gauntcp(l,m,lo,mo,li,mi) + (-1)**mo * gauntcp(l,m,lo,-mo,li,mi) ) / sqrt( 2._db )
      gi = 0._db
    else
      gi = ( (-1)**mo * gauntcp(l,m,lo,mo,li,mi) - gauntcp(l,m,lo,-mo,li,mi) ) / sqrt( 2._db )
      gr = 0._db
    endif

  else

    if( mi == 0 ) then
      gr = gauntc(l,m,lo,mo,li,mi)
      gi = 0._db
    elseif( mi > 0 ) then
      gr = gauntc(l,m,lo,mo,li,mi) / sqrt( 2._db )
      gi = gauntc(l,m,lo,mo,li,-mi) / sqrt( 2._db )
    else
      gr = (-1)**mi * gauntc(l,m,lo,mo,li,-mi) / sqrt( 2._db )
      gi = - (-1)**mi * gauntc(l,m,lo,mo,li,mi) / sqrt( 2._db )
    endif

  endif

  Gauntm = cmplx( gr, gi,db )

  return
end

!***********************************************************************

! Calcule le coefficient de Gaunt avec Y(l1,m1) complexe, Y(l2,m2) reel et Y(l3,m3) complexe

function Gaunt_crc(l1,m1,l2,m2,l3,m3,Ylm_comp)

  use declarations
  implicit none
  
  integer:: l1, l2, l3, m1, m2, m3
  
  complex(kind=db):: Gaunt_crc

  logical Ylm_comp

  real(kind=db):: Gauntc, Gauntcp, Gtot

  if( Ylm_comp ) then

    if( m2 == 0 ) then
      Gaunt_crc = cmplx( gauntcp(l1,m1,l2,m2,l3,m3), 0._db, db )
    elseif( m2 > 0 ) then
      Gtot = ( Gauntcp(l1,m1,l2,m2,l3,m3) + (-1)**m2 * Gauntcp(l1,m1,l2,-m2,l3,m3) ) / sqrt( 2._db ) 
      Gaunt_crc = cmplx( Gtot, 0._db, db )
    else
      Gtot = ( (-1)**m2 * Gauntcp(l1,m1,l2,m2,l3,m3) - Gauntcp(l1,m1,l2,-m2,l3,m3) ) / sqrt( 2._db )
      Gaunt_crc = cmplx( 0._db, Gtot, db )
    endif

  else

    Gaunt_crc = cmplx( Gauntc(L1,m1,l2,m2,l3,m3), 0._db, db )

  endif

  return
end

!***********************************************************************

! Calcule le pseudo coefficient de Gaunt pour la transition dipole magnetique
! G = 1/hbar < (l,m) | L_o + 2 S_o | (li,mi) >

function Gauntmag(isping,ispinf,mo,l,m,li,mi,Ylm_comp,Ylm_compi)

  use declarations
  implicit none

  integer:: i, ispinf, isping, j, l, li, m, mi, mo, mp, mq

  complex(kind=db):: Gauntmag

  logical:: Ylm_comp, Ylm_compi

  real(kind=db):: f

  Gauntmag = ( 0._db, 0._db )

  if( l /= li ) return

! boucle si harmo reel a droite
  do j = 1,2

    if( j == 1 ) then
      mq = mi 
    else
      mq = - mi
    endif

! boucle si harmo reelle a gauche
    do i = 1,2

      if( i == 1 ) then
        mp = m 
      else
        mp = - m
      endif

      f = 0._db

      select case(mo)

        case(-1)     ! Ly + 2*Sy
          if( isping == ispinf ) then
            if( mp == mq + 1 ) then
              f = - 0.5_db * sqrt( li*(li+1._db) - mq*(mq+1) )
            elseif( mp == mq - 1 ) then
              f = 0.5_db * sqrt( li*(li+1._db) - mq*(mq-1) )
            endif
          elseif( mp == mq ) then
            if( isping == 1 ) then
              f = 1._db
            else
              f = - 1._db
            endif
          endif

        case(0)     ! Lz + 2*Sz
          if( isping == ispinf .and. mp == mq ) then
            if( isping == 1 ) then
              f = mq + 1._db
            else
              f = mq - 1._db
            endif
          endif

        case(1)     ! Lx + 2*Sx
          if( isping == ispinf ) then
            if( mp == mq + 1 ) then
              f = 0.5_db * sqrt( li*(li+1._db) - mq*(mq+1) )
            elseif( mp == mi - 1 ) then
              f = 0.5_db * sqrt( li*(li+1._db) - mq*(mq-1) )
            endif
          elseif( mp == mq ) then
            f = 1._db
          endif

      end select

      if( .not. ( Ylm_comp .or. m == 0 ) ) then
        f = f / sqrt( 2._db )
        if(  m > 0 ) then
          if( i == 2 .and. mod(abs(m),2) == 1 ) f = - f
        else
          if( ( i == 1 .and. mod(abs(m),2) == 1 ) .or. i == 2 ) f = -f
        endif
      endif

      if( .not. ( Ylm_compi .or. mi == 0 ) ) then
        f = f / sqrt( 2._db )
        if(  mi > 0 ) then
          if( j == 2 .and. mod(abs(mi),2) == 1 ) f = - f
        else
          if( ( j == 1 .and. mod(abs(mi),2) == 1 ) .or. j == 2 ) f = -f
        endif
      endif

      if( mo == -1 ) then
        Gauntmag = Gauntmag + cmplx( 0._db, f, db )
      else
        Gauntmag = Gauntmag + cmplx( f, 0._db, db )
      endif

      if( Ylm_comp .or. m == 0 ) exit

    end do

    if( Ylm_compi .or. mi == 0 ) exit
    
  end do
  
! Le signe moins est car c'est dans le <bra, donc complexe conjugue
  if( .not. Ylm_comp .and. m < 0 ) Gauntmag = - img * Gauntmag

  if( .not. Ylm_compi .and. mi < 0 ) Gauntmag = img * Gauntmag

  return
end

!***********************************************************************

subroutine rot_tensor_2( mat2, rot_int )

  use declarations
  implicit none

  integer:: i, j, k
  complex(kind=db):: cmat
  complex(kind=db), dimension(3,3):: mat, mat2

  real(kind=db), dimension(3,3):: rot_int

  mat(:,:) = (0._db, 0._db)
  do i = 1,3
    do j = 1,3
      do k = 1,3
        cmat = sum( rot_int(j,:) * mat2(k,:) )
        mat(i,j) = mat(i,j) + rot_int(i,k) * cmat
      end do
    end do
  end do
  mat2(:,:) = mat(:,:)

  return
end

!***********************************************************************

subroutine rot_tensor_3( mat3, rot_int )

  use declarations
  implicit none

  integer:: i, j, k, l, m
  complex(kind=db):: cmas, cmat
  complex(kind=db), dimension(3,3,3):: mat, mat3

  real(kind=db), dimension(3,3):: rot_int

  mat(:,:,:) = (0._db, 0._db)
  do i = 1,3
    do j = 1,3
      do k = 1,3
        do l = 1,3
          cmas = (0._db, 0._db)
          do m = 1,3
            cmat = sum( rot_int(k,:) * mat3(l,m,:) )
            cmas = cmas + rot_int(j,m) * cmat
          end do
          mat(i,j,k) = mat(i,j,k) + rot_int(i,l) * cmas
        end do
      end do
    end do
  end do
  mat3(:,:,:) = mat(:,:,:)

  return
end

!***********************************************************************

subroutine rot_tensor_4( Mat4, Rot_int )

  use declarations
  implicit none

  integer:: i, j, k, l, m, n, n1
  complex(kind=db):: cmar, cmas, cmat
  complex(kind=db), dimension(3,3,3,3):: mat, mat4

  real(kind=db), dimension(3,3):: rot_int

  mat(:,:,:,:) = (0._db, 0._db)

  do i = 1,3
    do j = 1,3
      do k = 1,3
        do l = 1,3
          do m = 1,3
            cmar = (0._db, 0._db)
            do n = 1,3
              cmas = (0._db, 0._db)
              do n1 = 1,3
                cmat = sum( rot_int(l,:) * mat4(m,n,n1,:) )
                cmas = cmas + rot_int(k,n1) * cmat
              end do
              cmar = cmar + rot_int(j,n) * cmas
            end do
            mat(i,j,k,l) = mat(i,j,k,l) + rot_int(i,m) * cmar
          end do
        end do
      end do
    end do
  end do

  mat4(:,:,:,:) = mat(:,:,:,:)

  return
end

!***********************************************************************

subroutine rot_tensor_6( Mat6, Rot_int )

  use declarations

  implicit none

  integer:: i1, i2, i3, i4, i5, i6, j1, j2, j3, j4, j5

  complex(kind=db):: Cm1, Cm2, Cm3, Cm4, Cm5, Cm6
  complex(kind=db), dimension(3,3,3,3,3,3):: Mat, Mat6

  real(kind=db), dimension(3,3):: Rot_int

  Mat(:,:,:,:,:,:) = (0._db, 0._db)

  do i1 = 1,3
    do i2 = 1,3
      do i3 = 1,3
        do i4 = 1,3
          do i5 = 1,3
            do i6 = 1,3

              Cm1 = (0._db, 0._db)
              do j1 = 1,3
                Cm2 = (0._db, 0._db)
                do j2 = 1,3
                  Cm3 = (0._db, 0._db)
                  do j3 = 1,3
                    Cm4 = (0._db, 0._db)
                    do j4 = 1,3
                      Cm5 = (0._db, 0._db)
                      do j5 = 1,3
                        Cm6 = sum( Rot_int(i6,:) * Mat6(j1,j2,j3,j4,j5,:) )
                        Cm5 = Cm5 + Rot_int(i5,j5) * Cm6
                      end do
                      Cm4 = Cm4 + Rot_int(i4,j4) * Cm5
                    end do
                    Cm3 = Cm3 + Rot_int(i3,j3) * Cm4
                  end do
                  Cm2 = Cm2 + Rot_int(i2,j2) * Cm3
                end do
                Cm1 = Cm1 + Rot_int(i1,j1) * Cm2
              end do

              Mat(i1,i2,i3,i4,i5,i6) = Mat(i1,i2,i3,i4,i5,i6) + Cm1

            end do
          end do
        end do
      end do
    end do
  end do

  Mat6(:,:,:,:,:,:) = Mat(:,:,:,:,:,:)

  return
end

!***********************************************************************

subroutine extract_coabs(Ang_rotsup,Core_resolved,Green_int,icheck,ie,isymext,multi_run,Multipole, &
              n_oo,nenerg,ninit1,ninitlr,nom_fich_extract,Rotsup,secdd,secdd_m,secdo,secdo_m,secdq, &
              secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Tensor_rot,Tddft)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer:: eof, he, hs, je, js, jhe, jhs, n_oo

  character(len=132) mot, nom_fich_extract

  complex(kind=db), dimension(3,3,ninitlr,0:0):: secdd, secdd_m, secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitlr,0:0):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:0):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:0):: secoo, secoo_m
  complex(kind=db), dimension(3,3):: secdd_t, secmd_t, secmm_t
  complex(kind=db), dimension(3,3,3):: secdq_t
  complex(kind=db), dimension(3,3,3,3):: secdo_t, secqq_t
  complex(kind=db), dimension(3,3,3,3,3,3):: Mat6

  logical:: Comp, Core_resolved, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Green_int, M1M1, Tddft, Tensor_rot

  logical, dimension(10):: Multipole

  real(kind=db), dimension(3):: Ang_rotsup
  real(kind=db), dimension(3,3):: rot_tem, Rotsup

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  if( ie == 1 ) then
    open(1, file = nom_fich_extract, status='old', iostat=istat)
    if( istat /= 0 ) call write_open_error(nom_fich_extract,istat,1)

! Si les tenseurs importes doivent subir une rotation
    Tensor_rot = .false.
    rotsup = 0._db
    do i = 1,3
      rotsup(i,i) = 1._db
    end do
    do i = 1,3
      if( abs( ang_rotsup(i) ) < eps6 ) cycle
      tensor_rot = .true.
      ang_rotsup(i) = ang_rotsup(i) * pi / 180
      cosa = cos( ang_rotsup(i) )
      sina = sin( ang_rotsup(i) )
      j = mod(i+2,3) + 1
      k = mod(i,3) + 1
      l = mod(i+1,3) + 1
      rot_tem(j,j) = cosa; rot_tem(j,k) = sina; rot_tem(j,l) = 0._db
      rot_tem(k,j) = -sina; rot_tem(k,k) = cosa;
      rot_tem(k,l) = 0._db; rot_tem(l,j) = 0._db
      rot_tem(l,k) = 0._db; rot_tem(l,l) = 1._db

      rotsup = matmul( rot_tem, rotsup )
    end do

    if( tensor_rot .and. icheck > 0 ) then
      write(3,110)
      write(3,120) ( rotsup(i,1:3), i = 1,3 )
    endif

    i = 0
    do
      read(1,'(A)' ) mot
      if( mot(2:15) /= 'Absorbing atom' ) cycle
      i = i + 1
      if( i == multi_run ) exit
    end do

    if( tddft ) then
      do
        read(1,'(A)',iostat=eof) mot
        if( eof /= 0 ) exit
        if( mot(2:12) == 'Cycle TDDFT' .or. mot(2:12) == 'Cycle Tddft') exit
      end do
    endif

  endif

  do
    read(1,'(A)') mot
    if( mot(7:11) == 'Coabs' ) exit
  end do

  call extract_tens(Comp,Core_resolved,E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,Green_int,isymext,M1M1,n_oo,ninit1,ninitlr,secdd, &
          secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m)

  if( tensor_rot ) then

    do i_g = 1,ninitlr

      if( E1E1 ) then
        secdd_t(:,:) = secdd(:,:,i_g,0)
        call rot_tensor_2( secdd_t, rotsup )
        secdd(:,:,i_g,0) = secdd_t(:,:)
      endif

      if( E1E2 ) then
        secdq_t(:,:,:) = secdq(:,:,:,i_g,0)
        call rot_tensor_3( secdq_t, rotsup )
        secdq(:,:,:,i_g,0) = secdq_t(:,:,:)
      endif

      if( E2E2 ) then
        secqq_t(:,:,:,:) = secqq(:,:,:,:,i_g,0)
        call rot_tensor_4( secqq_t, rotsup )
        secqq(:,:,:,:,i_g,0) = secqq_t(:,:,:,:)
      endif

      if( E1E3 ) then
        secdo_t(:,:,:,:) = secdo(:,:,:,:,i_g,0)
        call rot_tensor_4( secdo_t, rotsup )
        secdo(:,:,:,:,i_g,0) = secdo_t(:,:,:,:)
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
                Mat6(:,je,he,:,js,hs) = secoo(:,jhe,:,jhs,i_g,0)
              end do
            end do
          end do
        end do
        call rot_tensor_6( Mat6, Rotsup )
        jhe = 0
        do je = 1,3
          do he = 1,3
            jhe = jhe + 1
            jhs = 0
            do js = 1,3
              do hs = 1,3
                jhs = jhs + 1
                secoo(:,jhe,:,jhs,i_g,0) = Mat6(:,je,he,:,js,hs)
              end do
            end do
          end do
        end do
      endif

      if( E1M1 ) then
        secmd_t(:,:) = secmd(:,:,i_g,0)
        call rot_tensor_2( secmd_t, rotsup )
        secmd(:,:,i_g,0) = secmd_t(:,:)
      endif

      if( M1M1 ) then
        secmm_t(:,:) = secmm(:,:,i_g,0)
        call rot_tensor_2( secmm_t, rotsup )
        secmm(:,:,i_g,0) = secmm_t(:,:)
      endif

      if( .not. ( Comp .and. Green_int ) ) cycle

      if( E1E1 ) then
        secdd_t(:,:) = secdd_m(:,:,i_g,0)
        call rot_tensor_2( secdd_t, rotsup )
        secdd_m(:,:,i_g,0) = secdd_t(:,:)
      endif

      if( E1E2 ) then
        secdq_t(:,:,:) = secdq_m(:,:,:,i_g,0)
        call rot_tensor_3( secdq_t, rotsup )
        secdq_m(:,:,:,i_g,0) = secdq_t(:,:,:)
      endif

      if( E2E2 ) then
        secqq_t(:,:,:,:) = secqq_m(:,:,:,:,i_g,0)
        call rot_tensor_4( secqq_t, rotsup )
        secqq_m(:,:,:,:,i_g,0) = secqq_t(:,:,:,:)
      endif

      if( E1E3 ) then
        secdo_t(:,:,:,:) = secdo_m(:,:,:,:,i_g,0)
        call rot_tensor_4( secdo_t, rotsup )
        secdo_m(:,:,:,:,i_g,0) = secdo_t(:,:,:,:)
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
                Mat6(:,je,he,:,js,hs) = secoo_m(:,jhe,:,jhs,i_g,0)
              end do
            end do
          end do
        end do
        call rot_tensor_6( Mat6, Rotsup )
        jhe = 0
        do je = 1,3
          do he = 1,3
            jhe = jhe + 1
            jhs = 0
            do js = 1,3
              do hs = 1,3
                jhs = jhs + 1
                secoo_m(:,jhe,:,jhs,i_g,0) = Mat6(:,je,he,:,js,hs)
              end do
            end do
          end do
        end do
      endif

      if( E1M1 ) then
        secmd_t(:,:) = secmd_m(:,:,i_g,0)
        call rot_tensor_2( secmd_t, rotsup )
        secmd_m(:,:,i_g,0) = secmd_t(:,:)
      endif

      if( M1M1 ) then
        secmm_t(:,:) = secmm_m(:,:,i_g,0)
        call rot_tensor_2( secmm_t, rotsup )
        secmm_m(:,:,i_g,0) = secmm_t(:,:)
      endif

    end do

  endif

  if( ie == nenerg ) Close(1)

  return
  110 format(/' Matrix rotation for the extracted tensors, rot_sup :')
  120 format(3x,3f9.5)
end

!***********************************************************************

subroutine extract_tens(Comp,Core_resolved,E1E1,E1E2,E1E3,E1M1, E2E2,E3E3,Green_int,isymext,M1M1,n_oo,ninit1,ninitlr,secdd, &
      secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm, secmm_m,secoo,secoo_m,secqq,secqq_m)

  use declarations
  implicit none

  integer:: he, hs, i_g, ipr, istat, isym, isymext, j_g, j1, je, &
    jhe, jhs, js, k, ke, ks, l, Multipole, n, n_g, n_oo, ninit1, ninitlr, nnombre
  character(len=132) mot

  logical comp, Core_resolved_e, Core_resolved, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Green_int, M1M1

  real(kind=db), dimension(3):: vi, vr, wi, wr
  real(kind=db), dimension(3,3):: matopsym, voi, vor, woi, wor

  complex(kind=db), dimension(3,3):: mat2
  complex(kind=db), dimension(3,3,3):: mat3
  complex(kind=db), dimension(3,3,3,3):: mat4
  complex(kind=db), dimension(3,3,3,3,3,3):: mat6
  complex(kind=db), dimension(3,3,ninitlr):: secdd, secdd_m, secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitlr):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitlr):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr):: secoo, secoo_m

  if( E1E1 ) secdd(:,:,:) = (0._db,0._db)
  if( E1E2 ) secdq(:,:,:,:) = (0._db,0._db)
  if( E2E2 ) secqq(:,:,:,:,:) = (0._db,0._db)
  if( E1E3 ) secdo(:,:,:,:,:) = (0._db,0._db)
  if( E3E3 ) secoo(:,:,:,:,:) = (0._db,0._db)
  if( E1M1 ) secmd(:,:,:) = (0._db,0._db)
  if( M1M1 ) secmm(:,:,:) = (0._db,0._db)
  if( E1E1 ) secdd_m(:,:,:) = (0._db,0._db)
  if( E1E2 ) secdq_m(:,:,:,:) = (0._db,0._db)
  if( E2E2 ) secqq_m(:,:,:,:,:) = (0._db,0._db)
  if( E1E3 ) secdo_m(:,:,:,:,:) = (0._db,0._db)
  if( E3E3 ) secoo_m(:,:,:,:,:) = (0._db,0._db)
  if( E1M1 ) secmd_m(:,:,:) = (0._db,0._db)
  if( M1M1 ) secmm_m(:,:,:) = (0._db,0._db)

  do
    read(1,'(A)') mot
    if( mot(2:4) == 'sec' .or. mot(2:7) == 'Tensor' ) exit
  end do

  n_g = 1
  Core_resolved_e = .false.

  backspace(1)
  backspace(1)
  backspace(1)
  read(1,'(A)') mot
  if( mot(13:16) == 'tate' .or. mot(26:29) == 'Edge') then
    if( mot(13:16) == 'tate' ) Core_resolved_e = .true.
    backspace(1)
    read(1,'(24x,i3)') n_g
  else
    n_g = ninitlr
  endif

  do j_g = 1,n_g

    if( ( Core_resolved_e .and. Core_resolved ) .or. ( .not. Core_resolved_e .and. .not. Core_resolved ) ) then
      i_g = j_g
    elseif( j_g <= ninit1 ) then
      i_g = 1
    else
      i_g = 2
    endif

    Boucle_Multipole: do Multipole = 1,7

      if( Multipole == 1 .and. .not. E1E1 ) cycle
      if( Multipole == 2 .and. .not. E1E2 ) cycle
      if( Multipole == 3 .and. .not. E2E2 ) cycle
      if( Multipole == 4 .and. .not. E1E3 ) cycle
      if( Multipole == 5 .and. .not. E1M1 ) cycle
      if( Multipole == 6 .and. .not. M1M1 ) cycle
      if( Multipole == 7 .and. .not. E3E3 ) cycle

      do

        read(1,'(A)',iostat=istat) mot

        if( istat /= 0 ) then
          call write_error
          do ipr = 6,9,3
            select case(Multipole)
              Case(1)
                write(ipr,110) ' E1-E1 '
              Case(2)
                write(ipr,110) ' E1-E2 '
              Case(3)
                write(ipr,110) ' E2-E2 '
              Case(4)
                write(ipr,110) ' E1-E3 '
              Case(5)
                write(ipr,110) ' E1-M1 '
              Case(6)
                write(ipr,110) ' M1-M1 '
              Case(7)
                write(ipr,110) ' E3-E3 '
            end select
          end do
          stop
        endif

        select case(Multipole)
          Case(1)
            if( mot(2:10) == 'Tensor_dd' ) exit
          Case(2)
            if( mot(2:10) == 'Tensor_dq' ) exit
            if( ( mot(2:7) == 'Tensor' .and. mot(9:10) /= 'dd' ) .or. mot(2:4) == '---' ) then
               backspace(1)
               Cycle Boucle_Multipole
            endif
          Case(3)
            if( mot(2:10) == 'Tensor_qq' ) exit
          Case(4)
            if( mot(2:10) == 'Tensor_do' ) exit
          Case(5)
            if( mot(2:10) == 'Tensor_md' ) exit
          Case(6)
            if( mot(2:10) == 'Tensor_mm' ) exit
          Case(7)
            if( mot(2:10) == 'Tensor_oo' ) exit
        end select

      end do

      if( mot(43:56) == 'Green integral' ) Green_int = .true.

      n = nnombre(1,1320)
      if( n == 12 .or. ( n == 6 .and. .not. Green_int ) .or. n == 36 .or. ( n == 18 .and. .not. Green_int ) ) then
        Comp = .true.
      else
        Comp = .false.
        wi(:) = 0._db
      endif

      select case(Multipole)
        Case(1,5,6)

          do ke = 1,3

            if( Green_int .and. comp ) then
              read(1,*) ( wr(k), wi(k), k = 1,3 ), ( vr(k), vi(k), k = 1,3 )
            elseif( comp ) then
              read(1,*) ( wr(k), wi(k), k = 1,3 )
            else
              read(1,*) wr(1:3)
            endif

            if( Multipole == 1 ) then
              if( Green_int .and. comp ) secdd_m(ke,:,i_g) = cmplx( vr(:), vi(:),db )
              secdd(ke,:,i_g) = cmplx( wr(:), wi(:),db )
            elseif( Multipole == 5 ) then
              if( Green_int .and. comp ) secmd_m(ke,:,i_g) = cmplx( vr(:), vi(:),db )
              secmd(ke,:,i_g) = cmplx( wr(:), wi(:),db )
            elseif( Multipole == 6 ) then
              if( Green_int .and. comp ) secmm_m(ke,:,i_g) = cmplx( vr(:), vi(:),db )
              secmm(ke,:,i_g) = cmplx( wr(:), wi(:),db )
            endif

          end do

        Case(2)

          do ke = 1,3
            do ks = 1,3
              if( Green_int .and. comp ) then
                read(1,*) ( wr(k), wi(k), k = 1,3 ), ( vr(k), vi(k), k = 1,3 )
              elseif( comp ) then
                read(1,*) ( wr(k), wi(k), k = 1,3 )
              else
                read(1,*) wr(1:3)
              endif
              if( Green_int .and. comp ) secdq_m(ke,ks,:,i_g) = cmplx( vr(:), vi(:),db )
              secdq(ke,ks,:,i_g) = cmplx( wr(:), wi(:),db )
            end do
            read(1,*)
            read(1,*)
          end do
          Backspace(1)

        Case(3)

          do js = 1,3
            do ks = 1,3
              do ke = 1,3
                if( Green_int .and. comp ) then
                  read(1,*) ( wr(k), wi(k), k = 1,3 ), ( vr(k), vi(k), k = 1,3 )
                elseif( comp ) then
                  read(1,*) ( wr(k), wi(k), k = 1,3 )
                else
                  read(1,*) wr(1:3)
                endif
                if( Green_int .and. comp ) secqq_m(ke,:,ks,js,i_g) = cmplx( vr(:), vi(:),db )
                secqq(ke,:,ks,js,i_g) = cmplx( wr(:), wi(:),db )
              end do
              read(1,*)
              read(1,*)
            end do
          end do
          Backspace(1)

        Case(4)

          do ke = 1,3
            do ks = 1,3
              do j1 = 1,3
                if( Green_int .and. comp ) then
                  read(1,*) ( wr(k), wi(k), k = 1,3 ), ( vr(k), vi(k), k = 1,3 )
                elseif( comp ) then
                  read(1,*) ( wr(k), wi(k), k = 1,3 )
                else
                  read(1,*) wr(1:3)
                endif
                if( Green_int .and. comp ) secdo_m(ke,ks,j1,:,i_g) = cmplx( vr(:), vi(:),db )
                secdo(ke,ks,j1,:,i_g) = cmplx( wr(:), wi(:),db )
              end do
              read(1,*)
              read(1,*)
            end do
          end do
          Backspace(1)

        Case(7)

          do hs = 1,3
            do js = 1,3
              jhs = 3 * ( js - 1 ) + hs
              do he = 1,3
                do je = 1,3
                  jhe = 3 * ( je - 1 ) + he
                  if( Green_int .and. comp ) then
                    read(1,*) (( wor(k,ks), woi(k,ks), k = 1,3 ), ks = 1,3), (( vor(k,ks), voi(k,ks), k = 1,3 ), ks = 1,3)
                  elseif( comp ) then
                    read(1,*) (( wor(k,ks), woi(k,ks), k = 1,3 ), ks=1,3)
                  else
                    read(1,*) ( wor(1:3,ks), ks = 1,3)
                  endif
                  do ks = 1,3
                    if( Green_int .and. comp ) secoo_m(:,jhe,ks,jhs,i_g) = cmplx( vor(:,ks), voi(:,ks),db )
                    secoo(:,jhe,ks,jhs,i_g) = cmplx( wor(:,ks), woi(:,ks),db )
                  end do
                end do
              end do
              read(1,*)
              read(1,*)
            end do
          end do
          Backspace(1)

       end select

    end do Boucle_Multipole

  end do

  do i_g = 1,ninitlr

    if( isymext == 1 ) cycle
    isym = abs( isymext )
    call opsym(isym,matopsym)

    if( E1E1 ) then
      mat2(:,:) = secdd(:,:,i_g)
      call rot_tensor_2( mat2, matopsym )
      if( isymext < 0 ) mat2(:,:) = conjg( mat2(:,:) )
      secdd(:,:,i_g) = mat2(:,:)
    endif

    if( E1E2 ) then
      mat3(:,:,:) = secdq(:,:,:,i_g)
      call rot_tensor_3( mat3, matopsym )
      if( isymext < 0 ) mat3(:,:,:) = conjg( mat3(:,:,:) )
      secdq(:,:,:,i_g) = mat3(:,:,:)
    endif

    if( E2E2 ) then
      mat4(:,:,:,:) = secqq(:,:,:,:,i_g)
      call rot_tensor_4( mat4, matopsym )
      if( isymext < 0 ) mat4(:,:,:,:) = conjg( mat4(:,:,:,:) )
      secqq(:,:,:,:,i_g) = mat4(:,:,:,:)
    endif

    if( E1E3 ) then
      mat4(:,:,:,:) = secdo(:,:,:,:,i_g)
      call rot_tensor_4( mat4, matopsym )
      if( isymext < 0 ) mat4(:,:,:,:) = conjg( mat4(:,:,:,:) )
      secdo(:,:,:,:,i_g) = mat4(:,:,:,:)
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
              Mat6(:,je,he,:,js,hs) = secoo(:,jhe,:,jhs,i_g)
            end do
          end do
        end do
      end do
      call rot_tensor_6( Mat6, Matopsym )
      if( isymext < 0 ) mat6(:,:,:,:,:,:) = - mat6(:,:,:,:,:,:)
      jhe = 0
      do je = 1,3
        do he = 1,3
          jhe = jhe + 1
          jhs = 0
          do js = 1,3
            do hs = 1,3
              jhs = jhs + 1
              secoo(:,jhe,:,jhs,i_g) = Mat6(:,je,he,:,js,hs)
            end do
          end do
        end do
      end do
    endif

    if( E1M1 ) then
      mat2(:,:) = secmd(:,:,i_g)
      call rot_tensor_2( mat2, matopsym )
      if( isymext < 0 ) mat2(:,:) = conjg( mat2(:,:) )
      secmd(:,:,i_g) = mat2(:,:)
    endif

    if( M1M1 ) then
      mat2(:,:) = secmm(:,:,i_g)
      call rot_tensor_2( mat2, matopsym )
      if( isymext < 0 ) mat2(:,:) = conjg( mat2(:,:) )
      secmm(:,:,i_g) = mat2(:,:)
    endif

    if( .not. ( Green_int .and. Comp ) ) cycle

    if( E1E1 ) then
      mat2(:,:) = secdd_m(:,:,i_g)
      call rot_tensor_2( mat2, matopsym )
      if( isymext < 0 ) mat2(:,:) = - mat2(:,:)
      secdd_m(:,:,i_g) = mat2(:,:)
    endif

    if( E1E2 ) then
      mat3(:,:,:) = secdq_m(:,:,:,i_g)
      call rot_tensor_3( mat3, matopsym )
      if( isymext < 0 ) mat3(:,:,:) = - mat3(:,:,:)
      secdq_m(:,:,:,i_g) = mat3(:,:,:)
    endif

    if( E2E2 ) then
      mat4(:,:,:,:) = secqq_m(:,:,:,:,i_g)
      call rot_tensor_4( mat4, matopsym )
      if( isymext < 0 ) mat4(:,:,:,:) = - mat4(:,:,:,:)
      secqq_m(:,:,:,:,i_g) = mat4(:,:,:,:)
    endif

    if( E1E3 ) then
      mat4(:,:,:,:) = secdo_m(:,:,:,:,i_g)
      call rot_tensor_4( mat4, matopsym )
      if( isymext < 0 ) mat4(:,:,:,:) = - mat4(:,:,:,:)
      secdo_m(:,:,:,:,i_g) = mat4(:,:,:,:)
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
              Mat6(:,je,he,:,js,hs) = secoo_m(:,jhe,:,jhs,i_g)
            end do
          end do
        end do
      end do
      call rot_tensor_6( Mat6, Matopsym )
      if( isymext < 0 ) mat6(:,:,:,:,:,:) = - mat6(:,:,:,:,:,:)
      jhe = 0
      do je = 1,3
        do he = 1,3
          jhe = jhe + 1
          jhs = 0
          do js = 1,3
            do hs = 1,3
              jhs = jhs + 1
              secoo_m(:,jhe,:,jhs,i_g) = Mat6(:,je,he,:,js,hs)
            end do
          end do
        end do
      end do
    endif

    if( E1M1 ) then
      mat2(:,:) = secmd_m(:,:,i_g)
      call rot_tensor_2( mat2, matopsym )
      if( isymext < 0 ) mat2(:,:) = - mat2(:,:)
      secmd_m(:,:,i_g) = mat2(:,:)
    endif

    if( M1M1 ) then
      mat2(:,:) = secmm_m(:,:,i_g)
      call rot_tensor_2( mat2, matopsym )
      if( isymext < 0 ) mat2(:,:) = - mat2(:,:)
      secmm_m(:,:,i_g) = mat2(:,:)
    endif

  end do

  return
  110 format(//A,' not found in the extract file !'//)
end
