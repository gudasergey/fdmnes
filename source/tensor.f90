! FDMNES subroutines

! Calculation of the cartesian tensors

! < g_1 | o*( l_s m_s, k_s, j_s, l_s, irang ) | f_1 > < f_2 | o*( l_i, m_i, k_i, j_i, l_i, jrang  ) | g_2 >  

! n_Ec = ninit si TDDFT et XANES et Core_resolved
!      = nbseuil si TDDFT et XANES sans core_resolved
!      = 1 en DFT XANES
!      = 2 en Optic
! n_V = ninit si TDDFT et Core_resolved mais pas Optic
!     = nbseuil si TDDFT sans core_recolved mais pas Optic
!     = 1 en DFT ou optic

subroutine tenseur_car(coef_g,Core_resolved,Ecinetic, &
                Eimag,Energ,Enervide,Eseuil,FDM_comp,Final_optic,Final_tddft,Full_potential,Green,Green_int,Hubb_a,Hubb_d, &
                icheck,ie,ip_max,ip0,is_g,lmax,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
                n_Ec,n_oo,n_rel,n_RI,n_V,nbseuil,ns_dipmag,ndim2,nenerg_tddft,ninit1,ninit,ninitr,ninitv,nlm_pot,nlm_probe, &
                nlm_p_fp,nlmamax,nr,nrm,nspin,nspino,nspinp,numat,psii,r,Relativiste,Renorm, &
                Rmtg,Rmtsd,rof0,rot_atom_abs,Rot_int, &
                secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m, &
                secqq,secqq_m,Solsing,Solsing_only,Spinorbite,Tau_RI_abs,Taull,Tddft,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer:: h_s, hh_s, hh_i, h_i, i, icheck, ie, ief, initr, initlt, ip, ip_max, ip0, ipr, irang, irang1, is, isi, isol, &
    isp, isp1, isp2,isp3, isp4, ispfg, j, j_i, j_s, jh_i, jh_s, jj_i, jj_s, jjh_i, jjh_s, jrang, k, k_i, k_s, kk_i, kk_s, &
    l_i, l_s, lm, lm_i, lm_s, lmax, lmax_pot, lmomax, lomax, lseuil, m_hubb, m_i, m_s, mpinodes, mpirank, mpirank0, &
    n_Ec, n_oo, n_rel, n_RI, n_V, nbseuil, ndim2, nenerg_tddft, ninit1, ninit, &
    ninitr, ninitv, nh_i, nh_s, nj_i, nj_s, nlm_p_fp, nlm_pot, nlm_probe, nlm1g, &
    nlm2g, nlmamax, nr, nr_zet, nrang, nrm, ns_dipmag, nspin, nspino, nspinp, numat

  parameter( lomax = 3, lmomax = ( lomax + 1 )**2 )

  complex(kind=db), dimension(ninitr):: Te, Te_m, Ten, Ten_m, Tens, Tens_m
  complex(kind=db), dimension(3,3) :: Mat2
  complex(kind=db), dimension(3,3,3) :: Mat3
  complex(kind=db), dimension(3,3,3,3) :: Mat4
  complex(kind=db), dimension(3,3,3,3,3,3) :: Mat6
  complex(kind=db), dimension(3,3,ninitr,0:mpinodes-1):: secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitr,0:mpinodes-1):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitr,0:mpinodes-1):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,3,n_rel,ninitr,0:mpinodes-1):: secdd, secdd_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:mpinodes-1):: secoo, secoo_m
  complex(kind=db), dimension(lmomax,lmomax,n_rel,ninitr):: Tens_lm, Tens_lm_m
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull
  complex(kind=db), dimension(nlm_probe,nspinp,nlm_probe,nspinp,n_RI):: Tau_RI_abs
  complex(kind=db), dimension(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil):: rof0
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nlm_probe**2,nlm_p_fp**2,nspinp**2,nspino**2,ip_max-ip0+1,ip_max-ip0+1,nbseuil):: Rad_ioRIoi
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: Singul
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: rof

  integer, dimension(ninit,2):: m_g
  integer, dimension(ninit):: is_g
  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3):: msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi

  logical:: Core_resolved, Dip_rel, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, FDM_comp, Final_optic, &
    Final_tddft, Full_potential, Green, Green_int, Hubb_a, Hubb_d, lmoins1, lplus1, M_depend, M1M1, No_diag, NRIXS, Relativiste, &
    Renorm, Solsing, Solsing_only, Spinorbite, Tddft, Ylm_comp

  logical, dimension(10):: Multipole

  real(kind=db):: c, c0, c1, c12, c120, c3, c5, c8, clm_s, clm_i, Eimag, Energ, Enervide_t, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(1):: Energ_t
  real(kind=db), dimension(nspin):: Ecinetic_e, V0bd_e
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitr):: Tensi, Tensr
  real(kind=db), dimension(0:lomax,3,3,3,lmomax):: clm
  real(kind=db), dimension(3,3):: rot_atom_abs, Rot_int, rot_tem
  real(kind=db), dimension(ninit,2):: coef_g
  real(kind=db), dimension(nspin,n_Ec):: Ecinetic
  real(kind=db), dimension(n_Ec):: Eimag_t, Enervide
  real(kind=db), dimension(nspin,n_V):: V0bd
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato_e
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nlm_pot,nspin,n_V):: Vrato
  real(kind=db), dimension(nr,ip0:ip_max):: r_or_bess
  real(kind=db), dimension(:,:,:,:,:,:), allocatable:: zet
  real(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: roff_rr, roff0

  if( icheck > 1 ) write(3,100)

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  NRIXS = .false.
  Dip_rel = n_rel > 1

! Calculation of radial integrals for regular and irregular solutions
  if( Final_optic ) then

    No_diag = Full_potential .or. ( Hubb_a .and. .not. Hubb_d )
    M_depend = Spinorbite .or. Hubb_a .or. No_diag

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
    allocate( roff_rr(nlm1g,nlm1g,nlm2g,nlm2g,nspinp**2,nspino**2,ip0:ip_max) )
    nr_zet = 0
    allocate( zet(nr_zet,nlm1g,nlm2g,nspinp,nspino,n_Ec) )
    ief = 2
    allocate( roff0(ief-1,ief:n_Ec,nlm1g,nlm1g,nspinp,nspinp,nspino**2) )
    Eimag_t(:) = Eimag
    Energ_t(1) = Energ

    call radial_optic(Ecinetic,Eimag_t,Energ_t,Enervide,Full_potential,Hubb_a,Hubb_d,icheck,ief,ip_max,ip0, &
         lmax,lmax_pot,m_depend,m_hubb,1,n_Ec,n_V,nlm_pot,nlm1g,nlm2g,No_diag,nr,nr_zet,nspin,nspino,nspinp,numat,r,Relativiste, &
         Renorm,Rmtg,Rmtsd,roff_rr,roff0,Solsing,Spinorbite,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp,zet)
    
    deallocate( roff0, zet )

  else

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
  
    if( FDM_comp ) then

        Ecinetic_e(:) = Ecinetic(:,1)
        do isp = 1,nspin
          V0bd_e(isp) = V0bd(isp,1)
          Vrato_e(1:nr,:,isp) = Vrato(1:nr,:,isp,1)
        end do
        Enervide_t = Enervide(1)

      call Radial_integral_FDM_comp(Ecinetic_e,Eimag,Energ,Enervide_t,Eseuil,Full_potential,Hubb_a,Hubb_d,icheck, &
           ip_max,ip0,Lmax,Lmax_pot,m_hubb,n_RI,nbseuil,nlm_p_fp,nlm_pot,nlm_probe,nr,NRIXS,nrm,nspin,nspino, &
           nspinp,numat,psii,r,r_or_bess,Rad_ioRIoi,Relativiste,Renorm,Rmtg,Rmtsd,Spinorbite,Tau_RI_abs,V_hubb,V_intmax, &
           V0bd_e,Vrato_e,Ylm_comp)
 
! not used
      allocate( Singul(nlm_probe,nspinp,nlm_probe,nspinp,ip0:ip_max,ip0:ip_max,ninitv) )
      allocate( rof(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,ninitv) )
 
    else

      allocate( Singul(nlm_probe,nspinp,nlm_probe,nspinp,ip0:ip_max,ip0:ip_max,ninitv) )
      allocate( rof(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,ninitv) )

      rof(:,:,:,:,:,:) = (0._db, 0._db)
      Singul(:,:,:,:,:,:,:) = (0._db, 0._db)

      do initlt = 1,n_Ec
        if( Final_tddft ) then
          Ecinetic_e(:) = Ecinetic(:,initlt)
        else
          Ecinetic_e(:) = Ecinetic(:,1)
        endif
        if( Final_tddft ) then
          do isp = 1,nspin
            V0bd_e(isp) = V0bd(isp,initlt)
            Vrato_e(1:nr,:,isp) = Vrato(1:nr,:,isp,initlt)
          end do
        else
          do isp = 1,nspin
            V0bd_e(isp) = V0bd(isp,1)
            Vrato_e(1:nr,:,isp) = Vrato(1:nr,:,isp,1)
          end do
        endif
        Enervide_t = Enervide(initlt)

        call Radial_integral(Ecinetic_e,Eimag,Energ,Enervide_t,Eseuil,Final_tddft,Full_potential,Hubb_a,Hubb_d,icheck, &
           initlt,ip_max,ip0,lmax,lmax_pot,m_hubb,nbseuil,ninit1,ninitv,nlm_pot,nlm_probe,nlm_p_fp,nr,NRIXS,nrm, &
           nspin,nspino,nspinp,numat,psii,r,r_or_bess,Relativiste,Renorm,Rmtg,Rmtsd,rof,Singul,Solsing,Spinorbite,V_hubb,V_intmax, &
           V0bd_e,Vrato_e,Ylm_comp)
      
      end do

    endif

    if( Tddft .and. .not. Final_tddft ) then
      do isol = 1,nspino
        do isp = 1,nspinp
          do lm = 1,nlm_probe
            do i = 1,nbseuil
              if( nlm_p_fp == 1 ) then
                rof0(ie,lm,isp,isol,i) = rof(lm,1,isp,isol,0,i)
              else
! ce qui sert de base en tddft, ce sont les (l_i,m_i,s) et pas les (l,m) vrais.
                rof0(ie,lm,isp,isol,:) = sum( rof(1:nlm_probe,lm,isp,isol,0,i) )
              endif
            end do
          end do
        end do
      end do
    endif

  endif

  if( icheck > 1 ) write(3,110)

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

  secmm(:,:,:,mpirank) = (0._db,0._db)
  secmd(:,:,:,mpirank) = (0._db,0._db)
  secdq(:,:,:,:,mpirank) = (0._db,0._db)
  secqq(:,:,:,:,:,mpirank) = (0._db,0._db)
  secdo(:,:,:,:,:,mpirank) = (0._db,0._db)
  secdd(:,:,:,:,mpirank) = (0._db,0._db)

  if( E3E3 ) secoo(:,:,:,:,:,mpirank) = (0._db,0._db)

  if( Green_int ) then
    secmm_m(:,:,:,mpirank) = (0._db,0._db)
    secmd_m(:,:,:,mpirank) = (0._db,0._db)
    secdq_m(:,:,:,:,mpirank) = (0._db,0._db)
    secqq_m(:,:,:,:,:,mpirank) = (0._db,0._db)
    secdo_m(:,:,:,:,:,mpirank) = (0._db,0._db)
    secdd_m(:,:,:,:,mpirank) = (0._db,0._db)
    if( E3E3 ) secoo_m(:,:,:,:,:,mpirank) = (0._db,0._db)
  endif

  if( icheck > 1 .and. .not. Final_optic ) write(3,120) lseuil

! Loops over tensor rank
  do irang = irang1,nrang
    do jrang = irang,nrang
      Tens_lm(:,:,:,:) = (0._db,0._db)
      Tens_lm_m(:,:,:,:) = (0._db,0._db)
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

! Loops over tensor indices
      do k_s = 1,3

        if( irang == 1 .and. ldip(k_s) /= 1 ) cycle

        if( irang > 1 ) then
          nj_s = 3
        else
          nj_s = 1
        endif

        do j_s = 1,nj_s

          if( irang == 2 .and. lqua(k_s,j_s) /= 1 ) cycle

          if( irang == 3 ) then
            nh_s = 3
          else
            nh_s = 1
          endif

          do h_s = 1,nh_s

            if( irang == 3 .and. loct(k_s,j_s,h_s) /= 1 ) cycle
            jh_s = 3 * ( j_s - 1 ) + h_s

            do k_i = 1,3

              if( jrang == 1 .and. ldip(k_i) /= 1 ) cycle

              if( jrang > 1 ) then
                nj_i = 3
              else
                nj_i = 1
              endif

              do j_i = 1,nj_i

                if( jrang == 2 .and. lqua(k_i,j_i) /= 1 ) cycle

                if( jrang == 3 ) then
                  nh_i = 3
                else
                  nh_i = 1
                endif

                do h_i = 1,nh_i

                  if( jrang == 3 .and. loct(k_i,j_i,h_i) /= 1 ) cycle
                  jh_i = 3 * ( j_i - 1 ) + h_i

                  if( irang == 0 .and. jrang == 0 ) then
                    if( k_i < k_s ) cycle
                  elseif( irang == 1 .and. jrang == 1 ) then
                    if( k_i < k_s ) cycle
                    if( msymdd(k_s,k_i) == 0 .and. msymddi(k_s,k_i) == 0 ) cycle
                    if( sum( abs( secdd(k_s,k_i,1,:,mpirank) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 1 .and. jrang == 2 ) then
                    if( msymdq(k_s,k_i,j_i) == 0 .and. msymdqi(k_s,k_i,j_i) == 0 ) cycle
                    if( sum( abs( secdq(k_s,k_i,j_i,:,mpirank) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 2 .and. jrang == 2 ) then
                    if( msymqq(k_s,j_s,k_i,j_i) == 0 .and. msymqqi(k_s,j_s,k_i,j_i) == 0 ) cycle
                    if( sum( abs( secqq(k_s,j_s,k_i,j_i,:,mpirank) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 1 .and. jrang == 3 ) then
                    if( msymdo(k_s,k_i,j_i,h_i) == 0 .and. msymdoi(k_s,k_i,j_i,h_i) == 0 ) cycle
                  elseif( irang == 3 .and. jrang == 3 ) then
                    if( msymoo(k_s,jh_s,k_i,jh_i) == 0 .and. msymooi(k_s,jh_s,k_i,jh_i) == 0 ) cycle
                    if( sum( abs( secoo(k_s,jh_s,k_i,jh_i,:,mpirank) ) ) > 1.e-15_db ) cycle
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
                    Tens_m(:) = (0._db,0._db)

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
                            
                            if( Final_optic ) then
                              call tens_op(Core_resolved,Final_tddft,icheck,ip_max,ip0,irang,jrang,l_s,m_s,l_i,m_i,lmax,lmoins1, &
                                lplus1,M_depend,ns_dipmag,ndim2,ninitr,nlm1g,nlm2g,nlm_probe,nlm_p_fp,nspinp,nspino, &
                                roff_rr,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)
                            else
                              call tens_ab(coef_g,Core_resolved,Dip_rel,FDM_comp,Final_tddft,Full_potential,Green,Green_int, &
                                icheck,ip_max, &
                                ip0,irang,is_g,isp1,isp2,isp3,isp4,jrang,l_s,m_s,l_i,m_g,m_i,lmax,lmoins1,lplus1,lseuil,nbseuil, &
                                ns_dipmag,ndim2,ninit1,ninit,ninitv,ninitr,nlm_probe,nlm_p_fp,NRIXS,nspinp,nspino,Rad_ioRIoi,rof, &
                                Singul,Solsing,Solsing_only,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)
                              if( numat == 1) then
                           ! For hydrogen there is only 1 core state
                                Ten(:) = Ten(:) / 2
                                Ten_m(:) = Ten_m(:) / 2
                              endif
                            endif

                            Tens_lm(lm_s,lm_i,ispfg,:) = Ten(:)
                            Tens_lm_m(lm_s,lm_i,ispfg,:) = Ten_m(:)

                          endif

                          Tens(:) = Tens(:) + clm_s * clm_i * Tens_lm(lm_s,lm_i,ispfg,:)
                          Tens_m(:) = Tens_m(:) + clm_s * clm_i * Tens_lm_m(lm_s,lm_i,ispfg,:)
                        end do
                      end do
                    end do
                  end do ! fin boucle l_s

! Remplissage de la valeur calculee dans tous les elements du tenseur
! equivalent par symetrie.

                  Tensr(:) = real( Tens(:),db )
                  Tensi(:) = aimag( Tens(:) )

! M1-M1 (dipole magnetique - dipole magnetique)
                  if( irang == 0 .and. jrang == 0 ) then
                    secmm(k_s,k_i,:,mpirank) = Tens(:)
                    if( Green_int ) then
                      secmm(k_i,k_s,:,mpirank) = Tens(:)
                      secmm_m(k_s,k_i,:,mpirank) = Tens_m(:)
                      secmm_m(k_i,k_s,:,mpirank) = - Tens_m(:)
                    else
                      secmm(k_i,k_s,:,mpirank) = conjg( Tens(:) )
                    endif
                  endif

! M1-E1 (dipole magnetique - dipole electrique)
                  if( irang == 0 .and. jrang == 1 ) then
                    secmd(k_s,k_i,:,mpirank) = Tens(:)
                    if( Green_int ) secmd_m(k_s,k_i,:,mpirank) = Tens_m(:)
                  endif

! E1-E1 (Dipole-dipole)
                  if( irang == 1 .and. jrang == 1 ) then

                    do kk_s = 1,3
                      do kk_i = kk_s,3
                        if( abs(msymdd(kk_s,kk_i)) /= abs(msymdd(k_s,k_i)) .or. abs(msymddi(kk_s,kk_i)) &
                           /= abs(msymddi(k_s,k_i)) ) cycle
                        if( msymdd(kk_s,kk_i) /= 0 ) then
                          is = msymdd(kk_s,kk_i) / msymdd(k_s,k_i)
                        else
                          is = 0
                        endif
                        if( msymddi(kk_s,kk_i) /= 0 ) then
                          isi = msymddi(kk_s,kk_i) / msymddi(k_s,k_i)
                        else
                          isi = 0
                        endif
                        if( Green_int ) then
                          secdd(kk_s,kk_i,ispfg,:,mpirank) = is * Tens(:)
                          secdd(kk_i,kk_s,ispfg,:,mpirank) = is * Tens(:)
                          secdd_m(kk_s,kk_i,ispfg,:,mpirank) = isi * Tens_m(:)
                          secdd_m(kk_i,kk_s,ispfg,:,mpirank) = -isi * Tens_m(:)
                        else
! tenseur hermitique
                          Te(:) = cmplx( is*Tensr(:), isi*Tensi(:), db)
                          secdd(kk_s,kk_i,ispfg,:,mpirank) = Te(:)
                          secdd(kk_i,kk_s,ispfg,:,mpirank) = conjg( Te(:) )
                        endif
                      end do
                    end do

! Dipole-quadrupole
                  elseif( irang == 1 .and. jrang == 2 ) then
                    do kk_s = 1,3
                      do kk_i = 1,3
                        do jj_i = kk_i,3
                          if( abs(msymdq(kk_s,kk_i,jj_i)) /= abs(msymdq(k_s,k_i,j_i)) .or. abs(msymdqi(kk_s,kk_i,jj_i)) &
                             /= abs(msymdqi(k_s,k_i,j_i)) ) cycle
                          if( msymdq(kk_s,kk_i,jj_i) /= 0 ) then
                            is = msymdq(kk_s,kk_i,jj_i) / msymdq(k_s,k_i,j_i)
                          else
                            is = 0
                          endif
                          if( msymdqi(kk_s,kk_i,jj_i) /= 0 ) then
                            isi = msymdqi(kk_s,kk_i,jj_i) / msymdqi(k_s,k_i,j_i)
                          else
                            isi = 0
                          endif
                          if( Green_int ) then
                            secdq(kk_s,kk_i,jj_i,:,mpirank) =is*Tens(:)
                            secdq(kk_s,jj_i,kk_i,:,mpirank) =is*Tens(:)
                            secdq_m(kk_s,kk_i,jj_i,:,mpirank) = isi*Tens_m(:)
                            secdq_m(kk_s,jj_i,kk_i,:,mpirank) = isi*Tens_m(:)
                          else
                            Te(:)=cmplx(is*Tensr(:),isi*Tensi(:),db)
                            secdq(kk_s,kk_i,jj_i,:,mpirank) = Te(:)
                            secdq(kk_s,jj_i,kk_i,:,mpirank) = Te(:)
                          endif
                        end do
                      end do
                    end do

! Dipole-octupole
                  elseif( irang == 1 .and. jrang == 3 ) then
                    do kk_s = 1,3
                      do kk_i = 1,3
                        do jj_i = 1,3
                          do hh_i = 1,3
                            if( abs( msymdo(kk_s,kk_i,jj_i,hh_i) ) /= abs( msymdo(k_s,k_i,j_i,h_i) )  &
                              .or. abs( msymdoi(kk_s,kk_i,jj_i,hh_i) ) /= abs( msymdoi(k_s,k_i,j_i,h_i) ) ) cycle
                            if( msymdo(kk_s,kk_i,jj_i,hh_i) /= 0 ) then
                              is = msymdo(kk_s,kk_i,jj_i,hh_i) / msymdo(k_s,k_i,j_i,h_i)
                            else
                              is = 0
                            endif
                            if( msymdoi(kk_s,kk_i,jj_i,hh_i)/=0) then
                              isi = msymdoi(kk_s,kk_i,jj_i,hh_i) / msymdoi(k_s,k_i,j_i,h_i)
                            else
                              isi = 0
                            endif
                            if( Green_int ) then
                              secdo(kk_s,kk_i,jj_i,hh_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,kk_i,hh_i,jj_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,jj_i,kk_i,hh_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,hh_i,kk_i,jj_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,jj_i,hh_i,kk_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,hh_i,jj_i,kk_i,:,mpirank) = is * Tens(:)
                              secdo_m(kk_s,kk_i,jj_i,hh_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,kk_i,hh_i,jj_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,jj_i,kk_i,hh_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,hh_i,kk_i,jj_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,jj_i,hh_i,kk_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,hh_i,jj_i,kk_i,:,mpirank) = isi * Tens_m(:)
                            else
                              Te(:) = cmplx(is*Tensr(:),isi*Tensi(:),db)
                              secdo(kk_s,kk_i,jj_i,hh_i,:,mpirank) = Te(:)
                              secdo(kk_s,kk_i,hh_i,jj_i,:,mpirank) = Te(:)
                              secdo(kk_s,jj_i,kk_i,hh_i,:,mpirank) = Te(:)
                              secdo(kk_s,hh_i,kk_i,jj_i,:,mpirank) = Te(:)
                              secdo(kk_s,jj_i,hh_i,kk_i,:,mpirank) = Te(:)
                              secdo(kk_s,hh_i,jj_i,kk_i,:,mpirank) = Te(:)
                            endif
                          end do
                        end do
                      end do
                    end do

! Quadrupole-quadrupole
                  elseif( irang == 2 .and. jrang == 2 ) then
                    do kk_s = 1,3
                      do jj_s = 1,3
                        do kk_i = 1,3
                          do jj_i = 1,3
                            if( sum( abs( secqq(kk_s,jj_s,kk_i,jj_i,:,mpirank) ) ) > 1.e-15_db ) cycle
                            if( abs( msymqq(kk_s,jj_s,kk_i,jj_i) ) /= abs( msymqq(k_s,j_s,k_i,j_i) ) &
                                .or. abs( msymqqi(kk_s,jj_s,kk_i,jj_i) ) /= abs( msymqqi(k_s,j_s,k_i,j_i) ) ) cycle
                            if( msymqq(kk_s,jj_s,kk_i,jj_i) /= 0 ) then
                              is = msymqq(kk_s,jj_s,kk_i,jj_i) / msymqq(k_s,j_s,k_i,j_i)
                            else
                              is = 0
                            endif
                            if( msymqqi(kk_s,jj_s,kk_i,jj_i) /= 0 ) then
                              isi = msymqqi(kk_s,jj_s,kk_i,jj_i) / msymqqi(k_s,j_s,k_i,j_i)
                            else
                              isi = 0
                            endif
                            if( Green_int ) then
                              Te(:) = is * Tens(:)
                              secqq(kk_s,jj_s,kk_i,jj_i,:,mpirank)=Te(:)
                              secqq(jj_s,kk_s,kk_i,jj_i,:,mpirank)=Te(:)
                              secqq(kk_s,jj_s,jj_i,kk_i,:,mpirank)=Te(:)
                              secqq(jj_s,kk_s,jj_i,kk_i,:,mpirank)=Te(:)
                              secqq(kk_i,jj_i,kk_s,jj_s,:,mpirank)=Te(:)
                              secqq(jj_i,kk_i,kk_s,jj_s,:,mpirank)=Te(:)
                              secqq(kk_i,jj_i,jj_s,kk_s,:,mpirank)=Te(:)
                              secqq(jj_i,kk_i,jj_s,kk_s,:,mpirank)=Te(:)
                              Te(:) = isi * Tens_m(:)
                              secqq_m(kk_s,jj_s,kk_i,jj_i,:,mpirank) = Te(:)
                              secqq_m(jj_s,kk_s,kk_i,jj_i,:,mpirank) = Te(:)
                              secqq_m(kk_s,jj_s,jj_i,kk_i,:,mpirank) = Te(:)
                              secqq_m(jj_s,kk_s,jj_i,kk_i,:,mpirank) = Te(:)
                              secqq_m(kk_i,jj_i,kk_s,jj_s,:,mpirank) = - Te(:)
                              secqq_m(jj_i,kk_i,kk_s,jj_s,:,mpirank) = - Te(:)
                              secqq_m(kk_i,jj_i,jj_s,kk_s,:,mpirank) = - Te(:)
                              secqq_m(jj_i,kk_i,jj_s,kk_s,:,mpirank) = - Te(:)
                            else
                              Te(:)= cmplx( is*Tensr(:), isi*Tensi(:),db)
                              secqq(kk_s,jj_s,kk_i,jj_i,:,mpirank) = Te(:)
                              secqq(jj_s,kk_s,kk_i,jj_i,:,mpirank) = Te(:)
                              secqq(kk_s,jj_s,jj_i,kk_i,:,mpirank) = Te(:)
                              secqq(jj_s,kk_s,jj_i,kk_i,:,mpirank) = Te(:)
                              Te(:) = Conjg( Te(:) )
                              secqq(kk_i,jj_i,kk_s,jj_s,:,mpirank) = Te(:)
                              secqq(jj_i,kk_i,kk_s,jj_s,:,mpirank) = Te(:)
                              secqq(kk_i,jj_i,jj_s,kk_s,:,mpirank) = Te(:)
                              secqq(jj_i,kk_i,jj_s,kk_s,:,mpirank) = Te(:)
                            endif
                          end do
                        end do
                      end do
                    end do
! Octupole-octupole
                  elseif( irang == 3 .and. jrang == 3 ) then

                    do kk_s = 1,3
                      do jj_s = 1,3
                        do hh_s = 1,3
                          jjh_s = 3 * ( jj_s - 1 ) + hh_s
                          do kk_i = 1,3
                            do jj_i = 1,3
                              do hh_i = 1,3
                                jjh_i = 3 * ( jj_i - 1 ) + hh_i
                                i = msymoo(kk_s,jjh_s,kk_i,jjh_i)
                                j = msymoo(k_s,jh_s,k_i,jh_i)
                                if( abs(i) /= abs(j) ) cycle
                                if( j /= 0 ) then
                                  is = i / j
                                else
                                  is = 0
                                endif
                                i = msymooi(kk_s,jjh_s,kk_i,jjh_i)
                                j = msymooi(k_s,jh_s,k_i,jh_i)
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

                                secoo(kk_s,jjh_s,kk_i,jjh_i,:,mpirank) = Te(:)
                                if( Green_int ) secoo_m(kk_s,jjh_s,kk_i,jjh_i,:,mpirank) = Te_m(:)

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
              if( .not. Green_int ) cycle
              do ispfg = 1,n_rel
                if( ( irang /= 1 .or. jrang /= 1 ) .and. ispfg > 1 ) exit
                if( sum( abs(Tens_lm_m(lm_s,lm_i,ispfg,:)) ) > eps15 ) write(3,180) irang, jrang, l_s, m_s, l_i, m_i, &
                   ispfg, ( Tens_lm_m(lm_s,lm_i,ispfg,initr), initr = 1,ninitr )
              end do
            end do
          end do
        end do
      end do

    end do
  end do

! Rotation to get the tensors in R1 basis

    rot_tem = matmul( rot_int, transpose(rot_atom_abs) )

    do initr = 1,ninitr

      if( E1E1 ) then
        do ispfg = 1,n_rel
          mat2(:,:) = secdd(:,:,ispfg,initr,mpirank)
          call rot_tensor_2( mat2, rot_tem )
          secdd(:,:,ispfg,initr,mpirank) = mat2(:,:)
        end do
      endif

      if( E1E2 ) then
        mat3(:,:,:) = secdq(:,:,:,initr,mpirank)
        call rot_tensor_3( mat3, rot_tem )
        secdq(:,:,:,initr,mpirank) = mat3(:,:,:)
      endif

      if( E2E2 ) then
        mat4(:,:,:,:) = secqq(:,:,:,:,initr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secqq(:,:,:,:,initr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1E3 ) then
        mat4(:,:,:,:) = secdo(:,:,:,:,initr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secdo(:,:,:,:,initr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1M1 ) then
        mat2(:,:) = secmd(:,:,initr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmd(:,:,initr,mpirank) = mat2(:,:)
      endif

      if( M1M1 ) then
        mat2(:,:) = secmm(:,:,initr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmm(:,:,initr,mpirank) = mat2(:,:)
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
                Mat6(:,j_s,h_s,:,j_i,h_i) = secoo(:,jh_s,:,jh_i,initr,mpirank)
              end do
            end do
          end do
        end do
        call rot_tensor_6( Mat6, Rot_tem )
        jh_s = 0
        do j_s = 1,3
          do h_s = 1,3
            jh_s = jh_s + 1
            jh_i = 0
            do j_i = 1,3
              do h_i = 1,3
                jh_i = jh_i + 1
                secoo(:,jh_s,:,jh_i,initr,mpirank) = Mat6(:,j_s,h_s,:,j_i,h_i)
              end do
            end do
          end do
        end do
      endif

      if( .not. Green_int ) Cycle

      if( E1E1 ) then
        do ispfg = 1,n_rel
          mat2(:,:) = secdd_m(:,:,ispfg,initr,mpirank)
          call rot_tensor_2( mat2, rot_tem )
          secdd_m(:,:,ispfg,initr,mpirank) = mat2(:,:)
        end do
      endif

      if( E1E2 ) then
        mat3(:,:,:) = secdq_m(:,:,:,initr,mpirank)
        call rot_tensor_3( mat3, rot_tem )
        secdq_m(:,:,:,initr,mpirank) = mat3(:,:,:)
      endif

      if( E2E2 ) then
        mat4(:,:,:,:) = secqq_m(:,:,:,:,initr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secqq_m(:,:,:,:,initr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1E3 ) then
        mat4(:,:,:,:) = secdo_m(:,:,:,:,initr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secdo_m(:,:,:,:,initr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1M1 ) then
        mat2(:,:) = secmd_m(:,:,initr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmd_m(:,:,initr,mpirank) = mat2(:,:)
      endif

      if( M1M1 ) then
        mat2(:,:) = secmm_m(:,:,initr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmm_m(:,:,initr,mpirank) = mat2(:,:)
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
                Mat6(:,j_s,h_s,:,j_i,h_i) = secoo_m(:,jh_s,:,jh_i,initr,mpirank)
              end do
            end do
          end do
        end do
        call rot_tensor_6( Mat6, Rot_tem )
        jh_s = 0
        do j_s = 1,3
          do h_s = 1,3
            jh_s = jh_s + 1
            jh_i = 0
            do j_i = 1,3
              do h_i = 1,3
                jh_i = jh_i + 1
                secoo_m(:,jh_s,:,jh_i,initr,mpirank) = Mat6(:,j_s,h_s,:,j_i,h_i)
              end do
            end do
          end do
        end do
      endif

    end do

  if( Final_optic ) then
    deallocate( roff_rr )
  else
    deallocate( rof, Singul )
  endif
  
  return
  100 format(/' ---- Tenseur_car Radial ---------',100('-'))
  110 format(/' ---- Tenseur_car Tens_ab --------',100('-'))
  120 format(/' lseuil =',i2)
  130 format(/' -- irang =',i2,', jrang =',i2,' --')
  140 format(/' k_s, j_s, h_s =',3i3,',  k_i, j_i, h_i =',3i3,',  ispfg =',i2)
  150 format(/' l_s, m_s =',2i3,',  l_i, m_i =',2i3,', clm_s, clm_i =',1p, 2e13.5)
  160 format(/' Tensor by harmonics (basis R4) :'/, ' ir jr  l_s m_s  l_i m_i ispfg     Tens(Ylm,i=1,ninitr)')
  170 format(2i3,i4,i3,2x,i4,i3,i5,1p,20e15.7)
  180 format('i',i2,i3,i4,i3,2x,i4,i3,1p,20e15.7)

end

!***********************************************************************

subroutine tens_ab(coef_g,Core_resolved,Dip_rel,FDM_comp,Final_tddft,Full_potential,Green,Green_int,icheck,ip_max, &
                              ip0,irang,is_g,isp1,isp2,isp3,isp4,jrang,l_s,m_s,l_i,m_g,m_i,lmax,lmoins1,lplus1,lseuil,nbseuil, &
                              ns_dipmag,ndim2,ninit1,ninit,ninitv,ninitr,nlm_probe,nlm_p_fp,NRIXS,nspinp,nspino,Rad_ioRIoi,rof, &
                              Singul,Solsing,Solsing_only,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)

  use declarations
  implicit none

  integer, intent(in):: icheck, ip_max, ip0, irang, jrang, l_i, l_s, m_i, m_s, lseuil, nbseuil, ndim2, ninit1, &
      ninit, ninitv, ninitr, nlm_probe, nlm_p_fp, ns_dipmag, nspinp, nspino
  integer, dimension(ninit,2), intent(in):: m_g
  integer, dimension(ninit), intent(in):: is_g

  integer:: i_g_1, i_g_2, initl1, initl2, initr, is_dipmag, is_r1, is_r2, iseuil1, iseuil2, iso1, iso2, &
    isoo, isp1, isp2, isp3, isp4, ispf1, ispf2, ispinf1, ispinf2, isping1, isping2, &
    ispp, ispp_f1, ispp_f2, l_f1, l_f2, li, lm01, lm02, lm_f1, lm_f2, lmax, Lmm, lmp01, &
    lmp02, lmp_f1, lmp_f2, lms_f1, lms_f2, lp_f1, lp_f2, Lpp, m_f1, m_f2, mi1, mi2, mp_f1, mp_f2, mv1, mv2

  complex(kind=db):: cfe, cfs, Cg, rof_1, rof_2, Singul_e, Singul_s, Tau_rad, Tau_rad_i

  complex(kind=db):: Gaunt_i, Gaunt_nrixs, Gaunt_s, Gaunt_xrc, Gauntmag 
  complex(kind=db), dimension(ninitr):: Ten, Ten_m
  complex(kind=db), dimension(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,ninitv):: rof
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull
  complex(kind=db), dimension(nlm_probe,nspinp,nlm_probe,nspinp,ip0:ip_max,ip0:ip_max,ninitv):: Singul
  complex(kind=db), dimension(nlm_probe**2,nlm_p_fp**2,nspinp**2,nspino**2,ip_max-ip0+1,ip_max-ip0+1,nbseuil):: Rad_ioRIoi

  logical:: Core_resolved, Dip_rel, FDM_comp, Final_tddft, Full_potential, Green, Green_int, lmoins1, lplus1, NRIXS, &
            Solsing, Solsing_only, Spinorbite, Titre, Ylm_comp

  real(kind=db):: Ci_1, Ci_2, Ci2, J_initl1, J_initl2, Jz1, Jz2
  real(kind=db), dimension(ninit,2):: coef_g
  
  li = lseuil

  Ten(:) = (0._db,0._db)
  Ten_m(:) = (0._db,0._db)

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
        else
          is_r2 = iseuil2
        endif

! Loop over core state spin < ...  > < ... g_2 >
        do isping2 = 1,2  
        
!          if( .not. Final_tddft .and. isping1 /= isping2 ) cycle
          if( .not. ( Final_tddft .or. li /= 0 ) .and. isping1 /= isping2 ) cycle
          if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp4 /= isping2 ) cycle

          mi2 = m_g(i_g_2,isping2)
          Ci_2 = Coef_g(i_g_2,isping2)

          if( abs( Ci_2 ) < eps6 ) cycle

          Ci2 = Ci_1 * Ci_2

         if( icheck > 1 ) then
           write(3,110)
           write(3,120) i_g_1, isping1, Jz1, mi1, Ci_1, i_g_2, isping2, Jz2, mi2, Ci_2
         endif
         Titre = .true.

! Loop over spherical harmonics of final states < ... f_1 > < ... >
          do l_f1 = 0,lmax

            if( lplus1 .and. l_f1 < li + 1 ) cycle
            if( lmoins1 .and. l_f1 > li ) cycle
            if( l_f1 > li + irang .or. l_f1 < li - irang .or. mod(l_f1,2) /= mod(li+irang,2) ) cycle

            lm01 = l_f1**2 + l_f1 + 1
            do m_f1 = -l_f1,l_f1
              lm_f1 = lm01 + m_f1
              if( lm_f1 > nlm_probe ) cycle

! Loop over spin of final states < ... f_1 > < ... >
              do ispinf1 = 1,2

                if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp2 /= ispinf1 ) cycle

                ispf1 = min( ispinf1, nspinp )

! Il n'y a que pour l_s dipole magnetique qu'on peut avoir du spin-flip (ou pour l_s terme dipolaire relativiste)
                if( ispinf1 /= isping1 .and. ( NRIXS .or. irang > 1 .or. ( irang == 1 .and. &
                                                                           ( .not. Dip_rel .or. jrang /= 1 ) ) ) ) cycle

                if( NRIXS ) then
                  Gaunt_s = Gaunt_nrixs(l_f1,m_f1,l_s,m_s,li,mi1,Ylm_comp)
                elseif( irang == 0 ) then
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
                    if( lm_f2 > nlm_probe ) cycle

                    do ispinf2 = 1,2  ! spin etat final en sortie
                      ispf2 = min( ispinf2, nspinp )

                      if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp3 /= ispinf2 ) cycle

                      if( ispinf2 /= isping2 .and. ( NRIXS .or. jrang > 1 .or. ( jrang == 1 &
                                                                             .and. ( .not. Dip_rel .or. irang /= 1 ) ) ) ) cycle 

                      if( .not. ( Final_tddft .or. Spinorbite ) .and. ispinf2 /= ispinf1 ) cycle

                      if( NRIXS ) then
                        Gaunt_i = Gaunt_nrixs(l_f2,m_f2,l_i,m_i,li,mi2,Ylm_comp)
                      elseif( jrang == 0 ) then
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
              lmp_f1 = min( lmp01 + mp_f1, nlm_p_fp ) 

! Loop over spin-orbit solution index < ... f_1 > < ... >
              do ispp_f1 = 1,nspinp
                if( .not. Spinorbite .and. ispp_f1 /= ispf1 ) cycle
                iso1 = min( ispp_f1, nspino )

                if( Spinorbite .and. nlm_p_fp == 1 ) then
                  mv1 = mp_f1 - ispinf1 + iso1
                  if( mv1 > l_f1 .or. mv1 < -l_f1 ) cycle
                else
                  mv1 = mp_f1
                endif
                lms_f1 = ( lmp01 + mv1 - 1 ) * nspino + iso1

! Loop over harmonics of non spherical final states < ... > < f_2 ... >
                do lp_f2 = 0,lmax
                  if( .not. Full_potential .and. lp_f2 /= l_f2 ) cycle
                  lmp02 = lp_f2**2 + lp_f2 + 1

                  do mp_f2 = -lp_f2,lp_f2
                    if( nlm_p_fp == 1 .and. mp_f2 /= m_f2 ) cycle
                    lmp_f2 = min( lmp02 + mp_f2, nlm_p_fp ) 

! Loop over spin-orbit solution index < ... > < f_2 ... >
                    do ispp_f2 = 1,nspinp
                      if( .not. Spinorbite .and. ispp_f2 /= ispf2 ) cycle
                      iso2 = min( ispp_f2, nspino )

                      if( Spinorbite .and. nlm_p_fp == 1 ) then
                        mv2 = mp_f2 - ispinf2 + iso2
                        if( mv2 > l_f2 .or. mv2 < -l_f2 ) cycle
                      else
                        mv2 = mp_f2
                      endif
                      lms_f2 = ( lmp02 + mv2 - 1 ) * nspino + iso2

                      if( .not. Solsing_only ) then

                        cfe = ( 0._db, 0._db ) 
                          
                        if( FDM_comp ) then
                          ispp = (ispf1 - 1) * nspinp + ispf2   
                          isoo = (iso1 - 1) * nspino + iso2
                          Lmm = (L_f1**2 + L_f1 + m_f1) * nlm_probe + L_f2**2 + L_f2 + m_f2 + 1  
                          if( nlm_p_fp > 1 ) then
                            Lpp = (Lp_f1**2 + Lp_f1 + mp_f1) * nlm_p_fp + Lp_f2**2 + Lp_f2 + mp_f2 + 1
                          else
                            Lpp = 1
                          endif                         
                          cfe = Rad_ioRIoi(Lmm,Lpp,ispp,isoo,irang,jrang,is_r1)
  
                          ispp = (ispf2 - 1) * nspinp + ispf1   
                          isoo = (iso2 - 1) * nspino + iso1
                          Lmm = (L_f2**2 + L_f2 + m_f2) * nlm_probe + L_f1**2 + L_f1 + m_f1 + 1  
                          if( nlm_p_fp > 1 ) then
                            Lpp = (Lp_f2**2 + Lp_f2 + mp_f2) * nlm_p_fp + Lp_f1**2 + Lp_f1 + mp_f1 + 1
                          else
                            Lpp = 1
                          endif                         
                          cfs = Rad_ioRIoi(Lmm,Lpp,ispp,isoo,jrang,irang,is_r1)
                        
                        else

                          rof_1 = rof(lm_f1,lmp_f1,ispf1,iso1,irang,is_r1)
                          rof_2 = rof(lm_f2,lmp_f2,ispf2,iso2,jrang,is_r2)
 
                          if( Green ) then
                            cfe = rof_1 * rof_2 * Taull(lms_f1,lms_f2,initl1,initl2,ispinf1,ispinf2,is_dipmag)
                            cfs = rof_2 * rof_1 * Taull(lms_f2,lms_f1,initl2,initl1,ispinf2,ispinf1,is_dipmag)
                          else
                            cfe = conjg( rof_1 ) * rof_2 * Taull(lms_f1,lms_f2,initl1,initl2,ispinf1,ispinf2,is_dipmag)
                            cfs = conjg( rof_2 ) * rof_1 * Taull(lms_f2,lms_f1,initl2,initl1,ispinf2,ispinf1,is_dipmag)
                          endif
    
                        endif
                        
                      else

                        rof_1 = ( 0._db, 0._db )
                        rof_2 = ( 0._db, 0._db )
                        cfe = ( 0._db, 0._db )
                        cfs = ( 0._db, 0._db )

                      endif
                      
                      if( Solsing .and. i_g_1 == i_g_2 .and. lp_f1 == l_f1 .and. lp_f2 == l_f2  &
                                  .and. mp_f1 == m_f1 .and. mp_f2 == m_f2 .and. iso1 == iso2 .and. ispf1 == ispf2 &
                                  .and. ( .not. Spinorbite .or. iso1 == ispf1 ) ) then
                                                                                          ! index solution is diagonal
                        Singul_e = Singul(lm_f1+mv1-m_f1,ispf1,lm_f2+mv2-m_f2,ispf2,irang,jrang,is_r1)
                        Singul_s = Singul(lm_f2+mv2-m_f2,ispf2,lm_f1+mv1-m_f1,ispf1,jrang,irang,is_r1)
                        cfe = cfe + Singul_e
                        cfs = cfs + Singul_s
                      else
                        Singul_e = ( 0._db, 0._db ) 
                        Singul_s = ( 0._db, 0._db ) 
                      endif
                      
! Here one does not divide by pi, so convenient results seems multiplied by pi, when calculating atomic form factor.
! So, the normalization by pi must not be done in coabs to get the atomic form factor in DAFS.
                      if( Green_int ) then
                        Tau_rad = 0.5_db * ( cfe + cfs )
                        Tau_rad_i = 0.5_db * (cfe - cfs)
                        Ten(initr) = Ten(initr) + Real(Cg,db) * tau_rad + img * aimag(Cg)* tau_rad_i
                        Ten_m(initr) = Ten_m(initr) + Real(Cg,db) * tau_rad_i + img * aimag(Cg) * tau_rad
                      else
                        Tau_rad = 0.5_db * ( cfe - conjg( cfs ) )
                        Ten(initr) = Ten(initr) + Cg * Tau_rad
                      endif

                      if( icheck > 1 .and. abs( Tau_rad ) > eps15 ) then
                        if( Titre ) then
                          Titre = .false.
                          if( nlm_p_fp == 1 .and. FDM_comp ) then
                            write(3,125)
                          elseif( nlm_p_fp == 1 ) then
                            write(3,130)
                          elseif( FDM_comp ) then
                            write(3,132)
                          else
                            write(3,135)
                          endif
                        endif
                        if( nlm_p_fp == 1 .and. FDM_comp ) then
                          write(3,140) l_f1, m_f1, ispinf1, iso1, l_f2, m_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, cfe, cfs, Gaunt_i, Gaunt_s
                        elseif( nlm_p_fp == 1 ) then
                          write(3,140) l_f1, m_f1, ispinf1, iso1, l_f2, m_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, cfe, cfs, Gaunt_i, Gaunt_s, &
                            rof_1, rof_2, &
                            Taull(lms_f2,lms_f1,initl2,initl1,ispp_f2,ispp_f1,is_dipmag), &
                            Taull(lms_f1,lms_f2,initl1,initl2,ispp_f1,ispp_f2,is_dipmag), &
                            Singul_e, Singul_s
                        elseif( FDM_comp ) then
                          write(3,145) l_f1, m_f1, lp_f1, mp_f1, ispinf1, iso1, l_f2, m_f2, lp_f2, mp_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, cfe, cfs, Gaunt_i, Gaunt_s
                        else
                          write(3,145) l_f1, m_f1, lp_f1, mp_f1, ispinf1, iso1, l_f2, m_f2, lp_f2, mp_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, cfe, cfs, Gaunt_i, Gaunt_s, &
                            rof_1, rof_2, &
                            Taull(lms_f2,lms_f1,initl2,initl1,ispp_f2,ispp_f1,is_dipmag), &
                            Taull(lms_f1,lms_f2,initl1,initl2,ispp_f1,ispp_f2,is_dipmag), &
                            Singul_e, Singul_s
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
  Ten_m(:) = img * Ten_m(:)

  return
  110 format(/' ini1 isg1 Jz1 mi1  Ci1  ini2 isg2 Jz2 mi2  Ci2')
  120 format(2(2i4,f6.1,i3,f7.3))
  125 format('  l1  m1 is1 io1  l2  m2 is2 io2',15x,'Ten',27x,'dTen',24x,'Rad_ioRIoi',19x,'Rad_ioRIoi_t',23x, &
                'Gaunt_i',24x,'Gaunt_s')
  130 format('  l1  m1 is1 io1  l2  m2 is2 io2',15x,'Ten',27x,'dTen',28x,'cfe',27x,'cfs',27x, &
                'Gaunt_i',24x,'Gaunt_s',25x,'rof_s',26x,'rof_i',26x,'Tau_s',26x,'Tau_i',25x,'Singul_e',24x,'Singul_s')
  132 format('  l1  m1 lp1 mp1 is1 io1  l2  m2 lp2 mp2 is2 io2',15x,'Ten',27x,'dTen',24x,'Rad_ioRIoi',19x,'Rad_ioRIoi_t',23x, &
                'Gaunt_i',24x,'Gaunt_s')
  135 format('  l1  m1 lp1 mp1 is1 io1  l2  m2 lp2 mp2 is2 io2',15x,'Ten',27x,'dTen',28x,'cfe',27x,'cfs',27x, &
                'Gaunt_i',24x,'Gaunt_s',25x,'rof_s',26x,'rof_i',26x,'Tau_s',26x,'Tau_i',25x,'Singul_e',24x,'Singul_s')
  140 format(8i4,1p,16(1x,2e15.7))
  145 format(12i4,1p,16(1x,2e15.7))
end

!***********************************************************************

subroutine Tens_op(Core_resolved,Final_tddft,icheck,ip_max,ip0,irang,jrang,l_s,m_s,l_i,m_i,lmax,lmoins1, &
                                lplus1,M_depend,ns_dipmag,ndim2,ninitr,nlm1g,nlm2g,nlm_probe,nlm_p_fp,nspinp,nspino, &
                                roff_rr,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)

  use declarations
  implicit none

  integer, intent(in):: icheck, ip_max, ip0, irang, jrang, l_s, m_s, &
    l_i, m_i, ndim2, nlm1g, nlm2g, nlm_probe, nlm_p_fp, ns_dipmag, nspinp, nspino

  integer:: i, io1, io2, is1, is12, is2, is21, iso, iso_f1, iso_f2, iso_g1, iso_g2, isp, ispm, ispm_f1, ispm_f2, &
    ispm_g1, ispm_g2, isp_f1, isp_f2, isp_g1, isp_g2, l, l_f1, l_f2, l_g1, l_g2, lm, lm_f1, lm_f2, lm_g1, lm_g2, &
    lmax, lmp_f1, lmp_f2, lmp_g1, lmp_g2, lm_i, lms_f1, lms_f2, lms_g1, lmp, lmp0, lms_g2, lmss, lmss_f1, lmss_f2, lmss_g1, &
    lmss_g2, lp, m, m_f1, m_f2, m_g1, m_g2, mp, mv, ninitr, nlmss

  integer, dimension(:), allocatable:: iso_val, isp_val, ispm_val, l_val, lm_val, lmp_val, lms_val, m_val
 
  complex(kind=db):: cf, Gaunt_s, Gaunt_xrx, Gauntmag, Gaunt_i, Tau_rad
  complex(kind=db), dimension(ninitr):: Ten, Ten_m
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull

  logical:: Core_resolved, Final_tddft, lmoins1, lplus1, M_depend, Spinorbite, Ylm_comp

  real(kind=db):: Ro
  real(kind=db), dimension(nlm1g,nlm1g,nlm2g,nlm2g,nspinp**2,nspino**2,ip0:ip_max):: roff_rr

  Ten(:) = (0._db,0._db)
  Ten_m(:) = (0._db,0._db)

! Dimension calculation -----------------------------------

  nlmss = 0
  do isp = 1,2
    do l = 0,lmax
      do m = -l,l
! Loop over harmonics, non-spherical potential
        do lp = 0,lmax
          if( nlm_p_fp == 1 .and. lp /= l ) cycle
          do mp = -lp,lp
            if( nlm_p_fp == 1 .and. mp /= m ) cycle
! Loop over solutions
            do iso = 1,nspino
              if( Spinorbite ) then
                mv = mp - isp + iso
                if( mv > lp .or. mv < -lp ) cycle
              endif
              nlmss = nlmss + 1
            end do
          end do
        end do
      end do
    end do
  end do

  allocate( iso_val(nlmss) )
  allocate( isp_val(nlmss) )
  allocate( ispm_val(nlmss) )
  allocate( l_val(nlmss) )
  allocate( lm_val(nlmss) )
  allocate( lmp_val(nlmss) )
  allocate( lms_val(nlmss) )
  allocate( m_val(nlmss) )

! Building of quantum number table -------------------------------

  lmss = 0
! Loop over harmonics
  do l = 0,lmax

    do m = -l,l
      if( M_depend ) then
        lm = l**2 + l + 1 + m
      else
        lm = l + 1
      endif

! Loop over harmonics, non-spherical potential
      lmp = 0
      do lp = 0,lmax
        if( nlm_p_fp == 1 .and. lp /= l ) cycle
        lmp0 = nspino * (lp**2 + lp + 1)

        do mp = -lp,lp
          if( nlm_p_fp == 1 .and. mp /= m ) cycle
          lmp = lmp + 1

          do isp = 1,2
          ispm = min( isp, nspinp )

! Loop over solutions
            do iso = 1,nspino

              if( Spinorbite ) then
                mv = mp - isp + iso
                if( mv > lp .or. mv < -lp ) cycle
              else
                mv = mp
              endif
              lm_i = lmp0 + ( mv - 1 ) * nspino + iso

              lmss = lmss + 1

              iso_val(lmss) = iso
              isp_val(lmss) = isp 
              ispm_val(lmss) = ispm
              l_val(lmss) = l
              lm_val(lmss) = lm
              lmp_val(lmss) = lmp 
              lms_val(lmss) = lm_i 
              m_val(lmss) = m
              
            end do
          end do
        end do
      end do
    end do
  end do

  if( icheck > 1 ) then
    if( nlm_p_fp == 1 .and. Spinorbite .and. Final_tddft ) then
      write(3,110)
    elseif( nlm_p_fp == 1 .and. Final_tddft ) then
      write(3,112)
    elseif( nlm_p_fp == 1 .and. Spinorbite ) then
      write(3,115)
    elseif( nlm_p_fp == 1 ) then
      write(3,120)
    else
      write(3,130) 
    endif
  endif

! Calculation of the tensor components -------------------------------------------

! Loop over occupied states (outgoing)
  do lmss_g1 = 1,nlmss

    iso_g1 = iso_val(lmss_g1) 
    isp_g1 = isp_val(lmss_g1) 
    ispm_g1 = ispm_val(lmss_g1)
    l_g1 = l_val(lmss_g1)
    lm_g1 = lm_val(lmss_g1)
    lmp_g1 = lmp_val(lmss_g1) 
    lms_g1 = lms_val(lmss_g1) 
    m_g1 = m_val(lmss_g1) 

! Loop over non-occupied states ( outgoing )
    do lmss_f1 = 1,nlmss

      iso_f1 = iso_val(lmss_f1) 
      isp_f1 = isp_val(lmss_f1) 
      ispm_f1 = ispm_val(lmss_f1)
      l_f1 = l_val(lmss_f1)
      lm_f1 = lm_val(lmss_f1)
      lmp_f1 = lmp_val(lmss_f1) 
      lms_f1 = lms_val(lmss_f1) 
      m_f1 = m_val(lmss_f1) 

! spin-flip is possible only for magnetic dipole
      if( isp_f1 /= isp_g1 .and. irang /= 0 ) cycle

      if( ( lplus1 .and. l_f1 < l_g1 + 1 ) .or. ( lmoins1 .and. l_f1 > l_g1 )  ) cycle
      if( l_f1 > l_g1 + irang .or. l_f1 < l_g1 - irang .or. mod(l_f1,2) /= mod(l_g1+irang,2) ) cycle

! < g1 m_s f1 >< f2 m_i g2 >
      if( irang == 0 ) then
!        Gaunt_s = Gauntmag(l_f1,m_f1,isp_f1,m_s,l_g1,m_g1,isp_g1,Ylm_comp,Ylm_comp)
        Gaunt_s = conjg( Gauntmag(l_g1,m_g1,isp_g1,m_s,l_f1,m_f1,isp_f1,Ylm_comp,Ylm_comp) )
      else
        Gaunt_s = Gaunt_xrx(l_f1,m_f1,l_s,m_s,l_g1,m_g1,Ylm_comp)
      endif
      if( abs(Gaunt_s) < eps10 ) cycle

! Loop over non-occupied states ( incoming )
      do lmss_f2 = 1,nlmss

        iso_f2 = iso_val(lmss_f2) 
        isp_f2 = isp_val(lmss_f2) 
        ispm_f2 = ispm_val(lmss_f2)
        l_f2 = l_val(lmss_f2)
        lm_f2 = lm_val(lmss_f2)
        lmp_f2 = lmp_val(lmss_f2) 
        lms_f2 = lms_val(lmss_f2) 
        m_f2 = m_val(lmss_f2) 
    
! Loop over occupied states ( incoming )
        do lmss_g2 = 1,nlmss

          iso_g2 = iso_val(lmss_g2) 
          isp_g2 = isp_val(lmss_g2) 
          ispm_g2 = ispm_val(lmss_g2)
          l_g2 = l_val(lmss_g2)
          lm_g2 = lm_val(lmss_g2)
          lmp_g2 = lmp_val(lmss_g2) 
          lms_g2 = lms_val(lmss_g2) 
          m_g2 = m_val(lmss_g2) 

          if( isp_f2 /= isp_g2 .and. jrang /= 0 ) cycle
          
!          if( .not. ( Final_tddft .or. Spinorbite ) .and. ( isp_g2 /= isp_g1 .or. isp_f2 /= isp_f1 ) ) cycle
!          if( ( isp_g2 /= isp_g1 .or. isp_f2 /= isp_f1 ) ) cycle
          
          if( ( lplus1 .and. l_f2 < l_g2 + 1 ) .or. ( lmoins1 .and. l_f2 > l_g2 ) ) cycle
          if( l_f2 > l_g2 + jrang .or. l_f2 < l_g2 - jrang .or. mod(l_f2,2) /= mod(l_g2+jrang,2) ) cycle

! < g1 m_s f1 >< f2 m_i g2 >
          if( jrang == 0 ) then
            Gaunt_i = Gauntmag(l_f2,m_f2,isp_f2,m_i,l_g2,m_g2,isp_g2,Ylm_comp,Ylm_comp)
          else
            Gaunt_i = Gaunt_xrx(l_f2,m_f2,l_i,m_i,l_g2,m_g2,Ylm_comp)
          endif
          if( abs(Gaunt_i) < eps10 ) cycle

          is1 = ispm_f1 + (ispm_g1 - 1) * nspinp
          is2 = ispm_f2 + (ispm_g2 - 1) * nspinp
          io1 = iso_f1 + (iso_g1 - 1) * nspino
          io2 = iso_f2 + (iso_g2 - 1) * nspino

          Ro = roff_rr(lm_g1,lm_f1,lmp_g1,lmp_f1,is1,io1,irang) * roff_rr(lm_g2,lm_f2,lmp_g2,lmp_f2,is2,io2,jrang)

          if( Final_tddft ) then
            is12 = 2 * isp_g1 - 2 + isp_g2
            is21 = 2 * isp_g2 - 2 + isp_g1
     ! + conjg is because this Taull must be multiplied by - i to be classical tau like.
            Tau_rad = - 0.5 * Ro * ( Taull(lms_f1,lms_f2,lms_g1,lms_g2,isp_f1,isp_f2,is12) &
                               + conjg( Taull(lms_f2,lms_f1,lms_g2,lms_g1,isp_f2,isp_f1,is21) ) )
 !           Tau_rad = - Ro * cmplx( real( Taull(lms_f1,lms_f2,lms_g1,lms_g2,isp_f1,isp_f2,is12), db ), 0._db, db )
          else                       
            Tau_rad = - Ro * Taull(lms_f1,lms_f2,1,1,isp_f1,isp_f2,2) * Taull(lms_g2,lms_g1,1,1,isp_g2,isp_g1,1)
          endif

! Evenso Tau appear 2 times, one devide only 1 time by pi.
! Normalisation by pi must then do not to be done in coabs for dafs amplitude calculation.
          cf = Tau_rad * Conjg( Gaunt_s ) * Gaunt_i / pi
            
          Ten(1) = Ten(1) + cf
          if( Core_resolved ) then
            do i = 2,ninitr
              if( l_g1 == i-2 .and. l_g2 == i-2 ) Ten(i) = Ten(i) + cf
            end do
          endif

          if( icheck < 2  .or. abs( cf ) < 1.e-15_db ) cycle
          if( nlm_p_fp == 1 .and. Spinorbite .and. Final_tddft ) then
            write(3,150) l_g1, m_g1, isp_g1, iso_g1, l_f1, m_f1, isp_f1, iso_f1, &
                l_f2, m_f2, isp_f2, iso_f2, l_g2, m_g2, isp_g2, iso_g2, &
                Ten(1), cf, Taull(lms_f2,lms_f1,lms_g2,lms_g1,isp_f1,isp_f2,is12), Gaunt_s, Gaunt_i, & 
                roff_rr(lm_g1,lm_f1,lmp_g1,lmp_f1,is1,io1,irang), &
                roff_rr(lm_g2,lm_f2,lmp_g2,lmp_f2,is2,io2,jrang)
          elseif( nlm_p_fp == 1 .and. Final_tddft ) then
            write(3,155) l_g1, m_g1, isp_g1, l_f1, m_f1, isp_f1, &
                l_f2, m_f2, isp_f2, l_g2, m_g2, isp_g2, &
                Ten(1), cf, Taull(lms_f2,lms_f1,lms_g2,lms_g1,isp_f1,isp_f2,is12), Gaunt_s, Gaunt_i, & 
                roff_rr(lm_g1,lm_f1,lmp_g1,lmp_f1,is1,io1,irang), &
                roff_rr(lm_g2,lm_f2,lmp_g2,lmp_f2,is2,io2,jrang)
          elseif( nlm_p_fp == 1 .and. Spinorbite ) then
            write(3,150) l_g1, m_g1, isp_g1, iso_g1, l_f1, m_f1, isp_f1, iso_f1, &
                l_f2, m_f2, isp_f2, iso_f2, l_g2, m_g2, isp_g2, iso_g2, &
                Ten(1), cf, Taull(lms_g2,lms_g1,1,1,isp_g2,isp_g1,1), &
                Taull(lms_f1,lms_f2,1,1,isp_f1,isp_f2,2), Gaunt_s, Gaunt_i, & 
                roff_rr(lm_g1,lm_f1,lmp_g1,lmp_f1,is1,io1,irang), &
                roff_rr(lm_g2,lm_f2,lmp_g2,lmp_f2,is2,io2,jrang)
          elseif( nlm_p_fp == 1 ) then
            write(3,155) l_g1, m_g1, isp_g1, l_f1, m_f1, isp_f1, &
                l_f2, m_f2, isp_f2, l_g2, m_g2, isp_g2, &
                Ten(1), cf, Taull(lms_g2,lms_g1,1,1,isp_g2,isp_g1,1), &
                Taull(lms_f1,lms_f2,1,1,isp_f1,isp_f2,2), Gaunt_s, Gaunt_i, & 
                roff_rr(lm_g1,lm_f1,lmp_g1,lmp_f1,is1,io1,irang), &
                roff_rr(lm_g2,lm_f2,lmp_g2,lmp_f2,is2,io2,jrang)
          else
            write(3,160) l_g1, m_g1, lmp_g1, isp_g1, iso_g1, l_f1, m_f1, lmp_f1, isp_f1, iso_f1,&
                l_f2, m_f2, lmp_f2, isp_f2, iso_f2,  l_g2, m_g2, lmp_g2, isp_g2, iso_g2, &
                Ten(1), cf, Taull(lms_f1,lms_f2,lms_g2,lms_g1,isp_f1,isp_f2,1), Gaunt_s, Gaunt_i, & 
                roff_rr(lm_g1,lm_f1,lmp_g1,lmp_f1,is1,io1,irang), &
                roff_rr(lm_g2,lm_f2,lmp_g2,lmp_f2,is2,io2,jrang)
          endif

        end do    ! end of loop over non-occupied states, incoming
      end do    ! end of loop over non-occupied states, outgoing
    end do    ! end of loop over occupied states, incoming
  end do   ! end of loop over occupied states, outgoing

  deallocate( iso_val, isp_val, ispm_val, l_val, lm_val, lmp_val, lms_val, m_val )

  return
  110 format(/' lg1 mg1 sg1 og1 lf1 mf1 sf1 of1 lf2 mf2 sf2 of2 lg2 mg2 sg2 og2',11x,'Tens',21x,'d_Tens', &
              17x,'Tau_f1f2g1g2',14x,'Gaunte_g1f1',15x,'Gaunts_g2f2',10x,'rof_f1g1',5x,'rof_f2g2')
  112 format(/' lg1 mg1 sg1 lf1 mf1 sf1 lf2 mf2 sf2 lg2 mg2 sg2',11x,'Tens',21x,'d_Tens', &
              17x,'Tau_f1f2g1g2',14x,'Gaunte_g1f1',15x,'Gaunts_g2f2',10x,'rof_f1g1',5x,'rof_f2g2')
  115 format(/' lg1 mg1 sg1 og1 lf1 mf1 sf1 of1 lf2 mf2 sf2 of2 lg2 mg2 sg2 og22',11x,'Tens',21x,'d_Tens', &
              19x,'Tau_g1g2',18x,'Tau_f1f2',16x,'Gaunte_g1f1',15x,'Gaunts_g2f2',10x,'rof_f1g1',5x,'rof_f2g2')
  120 format(/' lg1 mg1 sg1 lf1 mf1 sf1 lf2 mf2 sf2 lg2 mg2 sg2',11x,'Tens',21x,'d_Tens', &
              19x,'Tau_g1g2',18x,'Tau_f1f2',16x,'Gaunte_g1f1',15x,'Gaunts_g2f2',10x,'rof_f1g1',5x,'rof_f2g2')
  130 format(/' lg1 mg1 pg1 sg1 og1 lf1 mf1 pf1 sf1 of1 lf2 mf2 pf2 sf2 of2 lg2 mg2 pg2 sg2 og2',11x,'Tens',21x,'d_Tens', &
              17x,'Tau_f1f2g1g2',14x,'Gaunte_g1f1',15x,'Gaunts_g2f2',10x,'rof_f1g1',3x,'rof_f2g2')
  150 format(i3,15i4,1p,14e13.5)
  155 format(i3,11i4,1p,14e13.5)
  160 format(i3,19i4,1p,14e13.5)
end

!***********************************************************************

! Gaunt coefficient calculation with complex harmonics = Int( Y(l1,m1)*Y(l2,m2)Y(l3,m3)dOmega )
! Formula from M. E. Rose p. 62

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

function Gaunt_r( l1, m1, l2, m2, l3, m3 )

  use declarations
  implicit none

  integer:: itr, ityfct1, ityfct2, ityfct3, itygam, ityp, itypi, l1, l2, l3, ltsti, ltst, m, m1, m2, m3, mgam, mi, nexp, ns, nsf
 
  real(kind=db):: ab, cgc, pre, Gaunt_r, gf, sqr2i
  
  if( mod( l3+l1+l2, 2 ) /= 0 .or. abs( l1 - l3 ) > l2 .or. ( l1 + l3 ) < l2 ) then
    Gaunt_r = 0._db
    return
  endif

  mi = abs( m1 )
  mgam = abs( m2 )
  m = abs( m3 )

  if( mi+mgam /= m .and. mi+m /= mgam .and. m+mgam /= mi ) then
    Gaunt_r = 0._db
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

  Gaunt_r = gf * ab

  return
end

!***********************************************************************

! Calcule le coefficient de Gaunt avec Y(l,m) complexe ou reel, Y(lo,mo) reel et Y(li,mi) complexe. 

function Gaunt_xrc(l,m,lo,mo,li,mi,Ylm_comp)

  use declarations
  implicit none
  
  integer:: l, li, lo, m, mi, mo
  
  complex(kind=db):: Gaunt_xrc

  logical:: Ylm_comp

  real(kind=db):: Gaunt_r, Gauntcp, gi, gr
  
  if( Ylm_comp ) then

    if( mo == 0 ) then
      gr = gauntcp(l,m,lo,mo,li,mi)
      gi = 0._db
    elseif( mo > 0 ) then
      gr = ( (-1)**mo * gauntcp(l,m,lo,mo,li,mi) + gauntcp(l,m,lo,-mo,li,mi) ) / sqrt( 2._db )
      gi = 0._db
    else
      gi = ( gauntcp(l,m,lo,mo,li,mi) - (-1)**mo * gauntcp(l,m,lo,-mo,li,mi) ) / sqrt( 2._db )
      gr = 0._db
    endif

  else

    if( mi == 0 ) then
      gr = Gaunt_r(l,m,lo,mo,li,mi)
      gi = 0._db
    elseif( mi > 0 ) then
      gr = (-1)**mi * Gaunt_r(l,m,lo,mo,li,mi) / sqrt( 2._db )
      gi = (-1)**mi * Gaunt_r(l,m,lo,mo,li,-mi) / sqrt( 2._db )
    else
      gr = Gaunt_r(l,m,lo,mo,li,-mi) / sqrt( 2._db )
      gi = - Gaunt_r(l,m,lo,mo,li,mi) / sqrt( 2._db )
    endif

  endif

  Gaunt_xrc = cmplx( gr, gi, db )

  return
end

!***********************************************************************

! Calcule le coefficient de Gaunt avec Y(l1,m1) complexe ou reel, Y(l2,m2) reel et Y(l3,m3) complexe ou reel

function Gaunt_xrx(l1,m1,l2,m2,l3,m3,Ylm_comp)

  use declarations
  implicit none
  
  integer:: l1, l2, l3, m1, m2, m3
  
  complex(kind=db):: Gaunt_xrx

  logical Ylm_comp

  real(kind=db):: Gaunt_r, Gauntcp, Gtot

  if( Ylm_comp ) then

    if( m2 == 0 ) then
      Gaunt_xrx = cmplx( gauntcp(l1,m1,l2,m2,l3,m3), 0._db, db )
    elseif( m2 > 0 ) then
      Gtot = ( Gauntcp(l1,m1,l2,m2,l3,m3) + Gauntcp(l1,m1,l2,-m2,l3,m3) ) / sqrt( 2._db ) 
      Gaunt_xrx = cmplx( Gtot, 0._db, db )
    else
      Gtot = ( (-1)**m2 * Gauntcp(l1,m1,l2,m2,l3,m3) - (-1)**m2 * Gauntcp(l1,m1,l2,-m2,l3,m3) ) / sqrt( 2._db )
      Gaunt_xrx = cmplx( 0._db, Gtot, db )
    endif

  else

    Gaunt_xrx = cmplx( Gaunt_r(l1,m1,l2,m2,l3,m3), 0._db, db )

  endif

  return
end

!***********************************************************************

! Calcule le coefficient de Gaunt modifie: Int ( Y_L1* Y_L2* Y_L3 dOmega )
! avec Y(l1,m1) complexe ou reel, Y(l2,m2) et Y(l3,m3) complexe
! Pour Y(l2,m2), on prend Y(l2,m2)* = (-1)^m2 Y(l2,-m2)  

function Gaunt_nrixs(l1,m1,l2,m2,l3,m3,Ylm_comp)

  use declarations
  implicit none
  
  integer:: l1, l2, l3, m1, m2, m3
  
  complex(kind=db):: Gaunt_nrixs

  logical Ylm_comp

  real(kind=db):: Gauntcp, Gtot

  if( Ylm_comp .or. m1 == 0 ) then
    Gaunt_nrixs = (-1)**m2 * cmplx( gauntcp(l1,m1,l2,-m2,l3,m3), 0._db, db )
  elseif( m1 > 0 ) then
    Gtot = (-1)**m2 * ( (-1)**m1 * Gauntcp(l1,m1,l2,-m2,l3,m3) + Gauntcp(l1,-m1,l2,-m2,l3,m3) ) / sqrt( 2._db ) 
    Gaunt_nrixs = cmplx( Gtot, 0._db, db )
  else
    Gtot = (-1)**m2 * ( Gauntcp(l1,m1,l2,-m2,l3,m3) - (-1)**m1 * Gauntcp(l1,-m1,l2,-m2,l3,m3) ) / sqrt( 2._db )
    Gaunt_nrixs = cmplx( 0._db, Gtot, db )
  endif

  return
end

!***********************************************************************

! Calcule le pseudo coefficient de Gaunt pour la transition dipole magnetique
! G = 1/hbar < (l,m,isp_f) | L_o + 2 S_o | (l_g,m_g,isp_g) >

function Gauntmag(l,m,isp_f,mo,l_g,m_g,isp_g,Ylm_comp,Ylm_comp_g)

  use declarations
  implicit none

  integer:: i, isp_f, isp_g, j, l, l_g, m, m_g, mo, m_cf, m_cg

  complex(kind=db):: Gauntmag

  logical:: Ylm_comp, Ylm_comp_g

  real(kind=db):: f

  Gauntmag = ( 0._db, 0._db )

  if( l /= l_g ) return

! Loop for real harmonics, right side
  do j = 1,2

    if( j == 1 ) then
      m_cg = m_g 
    else
      m_cg = - m_g
    endif

! Loop for real harmonics, left side
    do i = 1,2

      if( i == 1 ) then
        m_cf = m 
      else
        m_cf = - m
      endif

      f = 0._db

      select case(mo)

        case(-1)     ! Ly + 2*Sy
          if( isp_g == isp_f ) then
            if( m_cf == m_cg + 1 ) then
              f = - 0.5_db * sqrt( l_g * (l_g + 1._db) - m_cg * (m_cg + 1) )
            elseif( m_cf == m_cg - 1 ) then
              f = 0.5_db * sqrt( l_g * (l_g + 1._db) - m_cg * (m_cg - 1) )
            endif
          elseif( m_cf == m_cg ) then
            if( isp_g == 1 ) then
              f = 1._db
            else
              f = - 1._db
            endif
          endif

        case(0)     ! Lz + 2*Sz
          if( isp_g == isp_f .and. m_cf == m_cg ) then
            if( isp_g == 1 ) then
              f = m_cg + 1._db
            else
              f = m_cg - 1._db
            endif
          endif

        case(1)     ! Lx + 2*Sx
          if( isp_g == isp_f ) then
            if( m_cf == m_cg + 1 ) then
              f = 0.5_db * sqrt( l_g*(l_g + 1._db) - m_cg * (m_cg + 1) )
            elseif( m_cf == m_cg - 1 ) then
              f = 0.5_db * sqrt( l_g*(l_g + 1._db) - m_cg * (m_cg - 1) )
            endif
          elseif( m_cf == m_cg ) then
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

      if( .not. ( Ylm_comp_g .or. m_g == 0 ) ) then
        f = f / sqrt( 2._db )
        if(  m_g > 0 ) then
          if( j == 2 .and. mod(abs(m_g),2) == 1 ) f = - f
        else
          if( ( j == 1 .and. mod(abs(m_g),2) == 1 ) .or. j == 2 ) f = -f
        endif
      endif

      if( mo == -1 ) then
        Gauntmag = Gauntmag + cmplx( 0._db, f, db )
      else
        Gauntmag = Gauntmag + cmplx( f, 0._db, db )
      endif

      if( Ylm_comp .or. m == 0 ) exit

    end do

    if( Ylm_comp_g .or. m_g == 0 ) exit
    
  end do
  
! Le signe moins est car c'est dans le <bra, donc complexe conjugue
  if( .not. Ylm_comp .and. m < 0 ) Gauntmag = - img * Gauntmag

  if( .not. Ylm_comp_g .and. m_g < 0 ) Gauntmag = img * Gauntmag

  if( .not. Ylm_comp .and. mod(abs(m),2) == 1 ) Gauntmag = - Gauntmag
  if( .not. Ylm_comp_g .and. mod(abs(m_g),2) == 1 ) Gauntmag = - Gauntmag

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
              n_oo,n_rel,nenerg,ninit1,ninitr,nom_fich_extract,Rotsup,secdd,secdd_m,secdo,secdo_m,secdq, &
              secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Tensor_rot,Tddft)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer:: eof, he, hs, ispfg, je, js, jhe, jhs, n_oo, n_rel

  character(len=132) mot
  character(len=Length_name) nom_fich_extract

  complex(kind=db), dimension(3,3,n_rel,ninitr,0:0):: secdd, secdd_m
  complex(kind=db), dimension(3,3,ninitr,0:0):: secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitr,0:0):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitr,0:0):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:0):: secoo, secoo_m
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

  call extract_tens(Comp,Core_resolved,E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,Green_int,isymext,M1M1,n_oo,n_rel,ninit1, &
     ninitr,secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m)

  if( tensor_rot ) then

    do i_g = 1,ninitr

      if( E1E1 ) then
        do ispfg = 1,n_rel
          secdd_t(:,:) = secdd(:,:,ispfg,i_g,0)
          call rot_tensor_2( secdd_t, rotsup )
          secdd(:,:,ispfg,i_g,0) = secdd_t(:,:)
        end do
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
        do ispfg = 1,n_rel
          secdd_t(:,:) = secdd_m(:,:,ispfg,i_g,0)
          call rot_tensor_2( secdd_t, rotsup )
          secdd_m(:,:,ispfg,i_g,0) = secdd_t(:,:)
        end do
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

subroutine Extract_tens(Comp,Core_resolved,E1E1,E1E2,E1E3,E1M1,E2E2,E3E3,Green_int,isymext,M1M1,n_oo,n_rel,ninit1, &
      ninitr,secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m)

  use declarations
  implicit none

  integer:: he, hs, i_g, ipr, istat, isym, ispfg, isymext, j_g, j1, je, jhe, jhs, js, &
    k, ke, ks, Multipole, n, n_g, n_oo, n_rel, ninit1, ninitr, nnombre

  character(len=132):: mot

  logical:: comp, Core_resolved_e, Core_resolved, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Green_int, M1M1

  real(kind=db), dimension(3):: vi, vr, wi, wr
  real(kind=db), dimension(3,3):: matopsym, voi, vor, woi, wor

  complex(kind=db), dimension(3,3):: mat2
  complex(kind=db), dimension(3,3,3):: mat3
  complex(kind=db), dimension(3,3,3,3):: mat4
  complex(kind=db), dimension(3,3,3,3,3,3):: mat6
  complex(kind=db), dimension(3,3,n_rel,ninitr):: secdd, secdd_m
  complex(kind=db), dimension(3,3,ninitr):: secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitr):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitr):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr):: secoo, secoo_m

  if( E1E1 ) secdd(:,:,:,:) = (0._db,0._db)
  if( E1E2 ) secdq(:,:,:,:) = (0._db,0._db)
  if( E2E2 ) secqq(:,:,:,:,:) = (0._db,0._db)
  if( E1E3 ) secdo(:,:,:,:,:) = (0._db,0._db)
  if( E3E3 ) secoo(:,:,:,:,:) = (0._db,0._db)
  if( E1M1 ) secmd(:,:,:) = (0._db,0._db)
  if( M1M1 ) secmm(:,:,:) = (0._db,0._db)
  if( E1E1 ) secdd_m(:,:,:,:) = (0._db,0._db)
  if( E1E2 ) secdq_m(:,:,:,:) = (0._db,0._db)
  if( E2E2 ) secqq_m(:,:,:,:,:) = (0._db,0._db)
  if( E1E3 ) secdo_m(:,:,:,:,:) = (0._db,0._db)
  if( E3E3 ) secoo_m(:,:,:,:,:) = (0._db,0._db)
  if( E1M1 ) secmd_m(:,:,:) = (0._db,0._db)
  if( M1M1 ) secmm_m(:,:,:) = (0._db,0._db)

  Green_int = .false.

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
    n_g = ninitr
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

          do ispfg = 1,n_rel
          
            if( Multipole /= 1 .and. ispfg > 1 ) exit
            if( ispfg > 1 ) read(1,*)
              
            do ke = 1,3

              if( Green_int .and. comp ) then
                read(1,*) ( wr(k), wi(k), k = 1,3 ), ( vr(k), vi(k), k = 1,3 )
              elseif( comp ) then
                read(1,*) ( wr(k), wi(k), k = 1,3 )
              else
                read(1,*) wr(1:3)
              endif

              if( Multipole == 1 ) then
                if( Green_int .and. comp ) secdd_m(ke,:,ispfg,i_g) = cmplx( vr(:), vi(:),db )
                secdd(ke,:,ispfg,i_g) = cmplx( wr(:), wi(:),db )
              elseif( Multipole == 5 ) then
                if( Green_int .and. comp ) secmd_m(ke,:,i_g) = cmplx( vr(:), vi(:),db )
                secmd(ke,:,i_g) = cmplx( wr(:), wi(:),db )
              elseif( Multipole == 6 ) then
                if( Green_int .and. comp ) then
                  secmm_m(ke,:,i_g) = cmplx( vr(:), vi(:),db )
                elseif( Comp ) then
                  secmm(ke,:,i_g) = cmplx( wr(:), wi(:),db )
                else
                  secmm(ke,:,i_g) = cmplx( wi(:), wr(:),db )
                endif
              endif

            end do
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

  do i_g = 1,ninitr

    if( isymext == 1 ) cycle
    isym = abs( isymext )
    call opsym(isym,matopsym)

    if( E1E1 ) then
      do ispfg = 1,n_rel
        mat2(:,:) = secdd(:,:,ispfg,i_g)
        call rot_tensor_2( mat2, matopsym )
        if( isymext < 0 ) mat2(:,:) = conjg( mat2(:,:) )
        secdd(:,:,ispfg,i_g) = mat2(:,:)
      end do
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
      do ispfg = 1,n_rel
        mat2(:,:) = secdd_m(:,:,ispfg,i_g)
        call rot_tensor_2( mat2, matopsym )
        if( isymext < 0 ) mat2(:,:) = - mat2(:,:)
        secdd_m(:,:,ispfg,i_g) = mat2(:,:)
      end do
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

!*********************************************************************************************

! Calculation of S = Sum_if | <f| exp( iq.r ) |i> |^2
! exp(iqr) = 4 * pi * Sum_L ( i^l * j_l(qr) * Y_L^*(Omega_r) * Y_L^*(Omega_q) )

subroutine S_nrixs_cal(coef_g,Core_resolved,Ecinetic, &
                    Eimag,Energ,Enervide,Eseuil,FDM_comp,Final_tddft,Full_potential,Green,Green_int,Hubb_a, &
                    Hubb_d,icheck,l0_nrixs,lmax_nrixs,is_g,lmax,lmax_pot,lmoins1,lplus1,lseuil,m_g,m_hubb, &
                    mpinodes,mpirank,n_Ec,n_RI,n_V,nbseuil,ns_dipmag,ndim2, &
                    ninit1,ninit,ninitr,ninitv,nlm_pot,nlm_probe,nlm_p_fp,nq_nrixs,nr,nrm,nspin,nspino,nspinp, &
                    numat,psii,q_nrixs,r,Relativiste,Renorm,Rmtg,Rmtsd, &
                    S_nrixs,S_nrixs_m,Solsing,Solsing_only,Spinorbite,Tau_RI_abs,Taull, &
                    V_hubb,V_intmax,V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none
  
  integer, parameter::  ip_max = 0
  integer, parameter::  ip0 = 0

  integer:: icheck, initr, initlt, iq, isp, l_i, l_s, l0_nrixs, lmax, lmax_nrixs, lmax_pot, lm_i, lm_s, &
    lseuil, m_hubb, m_i, m_s, mpinodes, mpirank, &
    n_Ec ,n_RI, n_V, nbseuil, ndim2, ninit1, ninit, ninitr, ninitv, nlm_pot, nlm_probe, &
    nlm_p_fp, nq_nrixs, nr, nrm, ns_dipmag, nspin, nspino, nspinp, numat

  complex(kind=db):: dfac   
  complex(kind=db), dimension(ninitr):: Ten, Ten_m
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull
  complex(kind=db), dimension(nlm_probe,nspinp,nlm_probe,nspinp,n_RI):: Tau_RI_abs
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: Singul
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: rof
  complex(kind=db), dimension(nq_nrixs,(lmax_nrixs+1)**2,(lmax_nrixs+1)**2,ninitr,0:mpinodes-1):: S_nrixs, S_nrixs_m
  complex(kind=db), dimension(nlm_probe**2,nlm_p_fp**2,nspinp**2,nspino**2,Lmax_nrixs-L0_nrixs+1,lmax_nrixs-L0_nrixs+1,nbseuil):: &
                                                                                                  Rad_ioRIoi

  integer, dimension(ninit,2):: m_g
  integer, dimension(ninit):: is_g

  logical:: Core_resolved, FDM_comp, Final_tddft, Full_potential, Green, Green_int, &
    Hubb_a, Hubb_d, lmoins1, lplus1, Monocrystal, NRIXS, Relativiste, Renorm, RIXS, Solsing, Solsing_only, Spinorbite, &
    Ylm_comp

  real(kind=db):: Eimag, Energ, Enervide_t, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(3):: q_vec
  real(kind=db), dimension(nspin):: Ecinetic_e, V0bd_e
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninit,2):: coef_g
  real(kind=db), dimension(nspin,n_Ec):: Ecinetic
  real(kind=db), dimension(n_Ec):: Enervide
  real(kind=db), dimension(nspin,n_V):: V0bd
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato_e
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nlm_pot,nspin,n_V):: Vrato
  real(kind=db), dimension(nr,l0_nrixs:lmax_nrixs):: bessel
  real(kind=db), dimension(4,nq_nrixs):: q_nrixs

  NRIXS = .true.
  RIXS = .false.

  allocate( Singul(nlm_probe,nspinp,nlm_probe,nspinp,l0_nrixs:lmax_nrixs,l0_nrixs:lmax_nrixs,ninitv) )
  allocate( rof(nlm_probe,nlm_p_fp,nspinp,nspino,l0_nrixs:lmax_nrixs,ninitv) )
  S_nrixs(:,:,:,:,mpirank) = (0._db,0._db)
  S_nrixs_m(:,:,:,:,mpirank) = (0._db,0._db)

  do iq = 1,nq_nrixs

    if( icheck > 1 ) write(3,100) iq, q_nrixs(4,iq) / bohr
    
    rof(:,:,:,:,:,:) = (0._db, 0._db)
    Singul(:,:,:,:,:,:,:) = (0._db, 0._db)

    call cbessel(bessel,l0_nrixs,lmax_nrixs,nr,q_nrixs(4,iq),r)
  
    do initlt = 1,n_Ec
      if( Final_tddft ) then
        do isp = 1,nspin
          Ecinetic_e(isp) = Ecinetic(isp,initlt)
          V0bd_e(isp) = V0bd(isp,initlt)
          Vrato_e(1:nr,:,isp) = Vrato(1:nr,:,isp,initlt)
        end do
      else
        do isp = 1,nspin
          Ecinetic_e(isp) = Ecinetic(isp,1)
          V0bd_e(isp) = V0bd(isp,1)
          Vrato_e(1:nr,:,isp) = Vrato(1:nr,:,isp,1)
        end do
      endif
      Enervide_t = Enervide(initlt)

      call radial_integral(Ecinetic_e,Eimag,Energ,Enervide_t,Eseuil,Final_tddft,Full_potential,Hubb_a,Hubb_d,icheck, &
        initlt,lmax_nrixs,l0_nrixs,lmax,lmax_pot,m_hubb,nbseuil,ninit1,ninitv,nlm_pot,nlm_probe,nlm_p_fp,nr,NRIXS,nrm, &
        nspin,nspino,nspinp,numat,psii,r,bessel,Relativiste,Renorm,Rmtg,Rmtsd,rof,Singul,Solsing,Spinorbite,V_hubb,V_intmax, &
        V0bd_e,Vrato_e,Ylm_comp)

    end do

    if( icheck > 1 ) write(3,110) 
    
    q_vec(:) = q_nrixs(1:3,iq)
    Monocrystal = sum( abs(q_vec(:)) ) > eps10

    lm_s = 0   
    do l_s = l0_nrixs,lmax_nrixs
      do m_s = -l_s,l_s
        lm_s = lm_s + 1

        lm_i = 0
        do l_i = l0_nrixs,lmax_nrixs

          if( l_s == l_i ) then
            dfac = cmplx( quatre_pi**2, 0._db, db )
          else 
            dfac = img**(l_s - l_i) * quatre_pi**2
          endif

          do m_i = -l_i,l_i
            lm_i = lm_i + 1
    
            if( ( l_s /= l_i .or. m_s /= m_i ) .and. .not. Monocrystal ) cycle

            if( icheck > 1 ) write(3,150) l_s, m_s, l_i, m_i

            call tens_ab(coef_g,Core_resolved,.false.,FDM_comp,Final_tddft,Full_potential,Green,Green_int,icheck,lmax_nrixs, &
                                l0_nrixs,l_s,is_g,1,1,1,1,l_i,l_s,m_s,l_i,m_g,m_i,lmax,lmoins1,lplus1,lseuil,nbseuil,ns_dipmag, &
                                ndim2,ninit1,ninit,ninitv,ninitr,nlm_probe,nlm_p_fp,NRIXS,nspinp,nspino,Rad_ioRIoi,rof, &
                                Singul,Solsing,Solsing_only,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)

            S_nrixs(iq,lm_s,lm_i,:,mpirank) = Ten(:) * dfac
            if( Green_int ) S_nrixs_m(iq,lm_s,lm_i,:,mpirank) = Ten_m(:) * dfac

            if( icheck > 1 ) write(3,160) l_s, m_s, l_i, m_i, dfac, &
                                         ( S_nrixs(iq,lm_s,lm_i,initr,mpirank), initr = 1,ninitr )  

          end do
        end do
      end do
    end do ! end of loop over l_s

  end do

  deallocate( rof, Singul )
  
  return
  100 format(/' ---- S_nrixs_cal Radial ---------',100('-'),//'  iq =',i3,', q =',f6.3)
  110 format(/' ---- S_nrixs_cal Tens_ab --------',100('-'))
  120 format(/'  iq =',i3,', q =',f6.3,', q_vec =',3f10.5)
  130 format(/' Normalized q_vec in the absorbing atom local basis =',3f10.5)
  140 format(/'  iq =',i3,', q =',f6.3)
  150 format(/' l_s, m_s =',2i3,',  l_i, m_i =',2i3)
  160 format(/' l_s, m_s =',2i3,',  l_i, m_i =',2i3,', dfac =',1p,2e15.7,', S_nrixs(1..ninitr) =',28e15.7)

end


