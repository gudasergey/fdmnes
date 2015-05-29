! FDMNES subroutine

! Sous-ensemble de routines qui servent a la procedure optic.

!***********************************************************************

subroutine main_optic(alfpot,angxyz,Allsite,Atom_nonsph,axyz,Base_spin,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
          Densite_atom,dv0bdcF,E_cut,E_cut_imp,E_Fermi_man,Eclie,Eneg,Energ_t, &
          Extract,Eseuil,Full_atom,Full_potential,Full_self_abs,Green_int,Green_plus,hkl_dafs,Hubb_a,Hubb_d,icheck,iaabsi, &
          iabsorig,iapot,iaprotoi,ip_max,ip0,iprabs,iprabs_nonexc,isigpi,isymeq, &
          itabs,itypepr,jseuil,ldip,lmax_pot,lmax_probe,lmaxabs_t,lmaxat0,lmaxfree,lmoins1,loct,lplus1,lqua, &
          lseuil,ltypcal,m_hubb,Magnetic,Moyenne,MPI_host_num_for_mumps,mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq, &
          msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole,n_atom_0,n_atom_ind,n_atom_proto,n_multi_run, &
          n_oo,n_tens_max,natome,natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ninit1,ninitl,ninitlr,nlm_pot,nlm_probe, &
          nomabs,nomfich,nomfich_s,Nonexc_g,nphi_dafs,nphim,npldafs,nplr,nplrm, &
          npsom,nptmoy,nptmoy_out,nr,nrato,nrm,nseuil,nspin,nspino,nspinp,ntype, &
          numat,nxanout,pdp,phdafs,phdf0t,phdt,pol,Polarise,poldafse,poldafss,psii, &
          r,rato,Relativiste,Rmtg,Rmtsd,rot_atom_abs,Rot_int,rsato,rsbdc,rsbdc_out,Self_abs,Solsing_only, &
          Spherical_signal,Spherical_tensor,Spinorbite,Taull_opt,Taux_eq,tpt,V_intmax,V_hubb,V0muf, &
          Vcato,vec,vecdafse,vecdafss,Vhbdc,Vhbdc_out,Volume_maille,VxcbdcF,VxcbdcF_out,Vxcato,Workf,xyz,Ylm_comp)

  use declarations
  implicit none

  integer:: iaabsi, iabsorig, icheck_s, ie, ie_computer, ie_q, ie_s, ie_t, ip_max, ip0, iprabs_nonexc, &
    iprabs, iso1, iso2, isp, isp1, isp2, iss1, iss2, itabs, je, jseuil, lm1, lm1g, lm2, lm2g, lmax, lmax_pot, &
    lmax_probe, lmaxabs_t, lseuil, m_hubb, MPI_host_num_for_mumps, mpinodes, mpirank, mpirank0, &
    lmaxat0, n_atom_0, n_atom_ind, n_atom_proto, n_multi_run, n_oo, &
    n_tens_max, n_vr_0, n_vr_ind, natome, natomsym, ne_t, nbseuil, ncolm, ncolr, ncolt, &
    ndim1, ndim2, nenerg, nenerg_s, nenerg_tddft, nge, ninit1, ninitl, ninitlr, nlm_pot, nlm_probe, nlma, nlma2, &
    nlmamax, nphim, npldafs, nplr, nplrm, npoint, npoint_ns, npsom, nptmoy, nptmoy_out, nr, nrm, nseuil, nspin, &
    nspino, nspinp, ntype, numat, nxanout, Trace_k
  integer, dimension(30):: icheck
  integer, dimension(natome):: iaprotoi
  integer, dimension(natomsym):: isymeq
  integer, dimension(0:ntype):: nrato
  integer, dimension(0:n_atom_proto):: iapot, itypepr
  integer, dimension(ninitl):: is_g
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(ninitl,2):: m_g
  integer, dimension(npldafs,2):: isigpi
  integer, dimension(3,npldafs):: hkl_dafs
  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3)::  msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi

  parameter(ne_t = 2)
  parameter(ndim2 = 1)

  character(len=132):: nomfich, nomfich_s, nomfich_cal_convt
  character(len=13), dimension(nplrm):: ltypcal
  character(len=length_word), dimension(ncolm):: nomabs

  complex(kind=db):: f_avantseuil
  complex(kind=db), dimension(3,nplrm):: pol
  complex(kind=db), dimension(natomsym,npldafs):: phdafs
  complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdt
  complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss
  complex(kind=db), dimension(3,3,ninitlr,0:mpinodes-1):: secdd, &
    secdd_m, secmd, secmd_m, secmm, secmm_m, secdd_t, secdd_m_t, secmd_t, secmd_m_t, secmm_t, secmm_m_t
  complex(kind=db), dimension(3,3,3,ninitlr,0:mpinodes-1):: secdq, secdq_m, secdq_t, secdq_m_t
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:mpinodes-1):: secdo, &
    secdo_m, secqq, secqq_m, secdo_t, secdo_m_t, secqq_t, secqq_m_t
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:mpinodes-1):: secoo, secoo_m, secoo_t, secoo_m_t
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nenerg_s,nlm_probe,nlm_probe,nspinp,nspinp):: Taull_opt
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,2,2,ne_t,ndim2,ndim2):: Taull_abs
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: rof0

  logical:: Allsite, Atom_nonsph, Base_spin, Cartesian_tensor, Core_resolved, Dafs, Dafs_bio, E_Fermi_man, &
    E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Eneg, Energphot, Extract, Final_optic, Final_tddft, Full_atom, Full_potential, &
    Full_self_abs, Green_i, Green_int, Green_plus, Hubb_a, Hubb_d, lmaxfree, lmoins1, lplus1, M1M1, Magnetic, &
    Moyenne, Nonexc_g, Polarise, Relativiste, Self_abs, Solsing, &
    Solsing_only, Spherical_signal, Spherical_tensor, Spinorbite, Tddft, Xan_atom, Ylm_comp

  logical, dimension(10):: Multipole

  real(kind=sg):: time

  real(kind=db):: alfpot, delta_E, Densite_atom, E_cut, E_cut_imp, E_cut_optic, Ec_min, Eclie, Ecmax, &
    Enervide, p, Rmtsd, rsbdc, rsbdc_out, V_intmax, V0muf, Vhbdc, Vhbdc_out, Volume_maille, Workf

  real(kind=db), dimension(ne_t):: En
  real(kind=db), dimension(3):: angxyz, axyz
  real(kind=db), dimension(4):: tp
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(12):: tpt
  real(kind=db), dimension(nspin):: dv0bdcF, V0bdc, V0bdc_out, VxcbdcF, VxcbdcF_out
  real(kind=db), dimension(nenerg_s):: Energ_s, Energ_t
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitlr):: Epsii, sec_atom
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(natomsym) :: Taux_eq
  real(kind=db), dimension(0:n_atom_proto):: Rmtg
  real(kind=db), dimension(3,3):: rot_atom_abs, Rot_int
  real(kind=db), dimension(nplrm,2):: pdp
  real(kind=db), dimension(3,nplrm):: vec
  real(kind=db), dimension(ninitl,2):: coef_g
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(n_tens_max,0:natomsym):: Int_tens
  real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: rsato
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(3,npldafs,nphim):: vecdafse, vecdafss
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vxcato

  real(kind=db), dimension(:), allocatable:: Ecinetic_e, Eimag, Energ, rhons, rs, Vh, Vhns
  real(kind=db), dimension(:,:), allocatable:: Ecinetic, rho, V0bdc_e, Vr, Vxc
  real(kind=db), dimension(:,:,:,:), allocatable:: Vrato, Vrato_e

  if( icheck(1) > 0 ) write(3,100)

  if( E_Fermi_man ) then
    E_cut_optic = E_cut_imp
  else
    E_cut_optic = E_cut
  endif
  Energ_s(:) = Energ_t(:) + E_cut_optic

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  f_avantseuil = (0._db, 0._db)
  Energphot = .false.
  Final_optic = .true.
  Final_tddft = .false.
  Tddft = .false.
  sec_atom(:) = 0._db
  nlmamax = 0

  n_vr_0 = 1
  n_vr_ind = ne_t

  allocate( Ecinetic(nspin,ne_t) )
  allocate( Ecinetic_e(nspin) )
  allocate( V0bdc_e(nspin,ne_t) )

  Green_i = Green_int

! Elaboration de la grille optic
  call dim_grille_optic(E_cut_optic,Energ_s,mpirank0,nenerg,nenerg_s)

  allocate( Energ(nenerg) )
  allocate( Eimag(nenerg) )

  call grille_optic(E_cut_optic,Energ,Energ_s,icheck(29),nenerg,nenerg_s)

! On garde l'energie imaginaire nulle
  Eimag(:) = 0._db
  Solsing = .false.

  Xan_atom = .false.

  allocate( Vrato(0:nrm,nlm_pot,nspin,ne_t) )
  allocate( Vrato_e(nr,nlm_pot,nspin,ne_t) )

! Boucle sur l'energie
  nge = ( nenerg - 1 ) / mpinodes + 1

  boucle_energ: do je = 1,nge

    ie = ( je - 1 ) * mpinodes + mpirank + 1

    if( ie <= nenerg ) then

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

    do ie_s = 1,nenerg_s

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        tp(1) = real(time,db)
      endif

      En(1) = Energ_s(ie_s)
      En(2) = Energ_s(ie_s) + Energ(ie)
      if( En(1) > E_cut_optic - eps10 ) exit
      if( En(2) < E_cut_optic .or. En(2) > Energ_s(nenerg_s) ) cycle

      icheck_s = max( icheck(16), icheck(29) )
      if( icheck_s == 1 ) icheck_s = 0

      if( ie_s == nenerg_s ) then
        delta_E = Energ_s(ie_s) - Energ_s(ie_s-1)
      elseif( ie_s == 1 ) then
        delta_E = Energ_s(ie_s+1) - Energ_s(ie_s)
      elseif(  Energ_s(ie_s+1) > E_cut_optic ) then
        delta_E = E_cut_optic - 0.5_db*(Energ_s(ie_s) + Energ_s(ie_s-1))
      else
        delta_E = 0.5_db * ( Energ_s(ie_s+1) - Energ_s(ie_s-1) )
      endif

! Les variables ci-dessous ne sont pas utilisees
      npoint = 0
      allocate( rho(npoint,nspin) )
      allocate( rs(npoint) )
      allocate( Vh(npoint) )
      allocate( Vr(npoint,nspin) )
      allocate( Vxc(npoint,nspin) )
      npoint_ns = 0
      allocate( rhons(npoint_ns) )
      allocate( Vhns(npoint_ns) )
      Trace_k = 0
      Trace_p(:) = 0._db
 
      do ie_t = 1,ne_t

        Enervide = En(ie_t) - Workf

! Calcul du potentiel dans l'etat excite
        call potex(Atom_nonsph,axyz,alfpot,.true.,dv0bdcF,En(ie_t),Enervide,Final_optic,Full_atom, &
          iaabsi,iapot,iaprotoi,icheck_s,ie_t,iprabs_nonexc,iprabs,itabs,itypepr,Magnetic, &
          n_atom_0,n_atom_ind,n_atom_proto,n_vr_0,n_vr_ind,natome,nlm_pot,Nonexc_g,npoint,npoint_ns,npsom,nptmoy, &
          nptmoy_out,nrato,nrm,nspin,ntype,rato,rho,rhons,Rmtg,rs,rsato,rsbdc,rsbdc_out,Trace_k,Trace_p, &
          V_intmax,Vcato,Vh,Vhbdc,Vhbdc_out,Vhns,Vr,Vxc,VxcbdcF,VxcbdcF_out,Vrato,Vxcato,V0bdc,V0bdc_out,xyz)

        V0bdc_e(:,ie_t) = V0bdc(:)
        Vrato_e(1:nr,:,:,ie_t) = Vrato(1:nr,:,:,ie_t)
        Ecinetic(:,ie_t) = Enervide - V0bdc_e(:,ie_t)
        if( .not. Eneg ) Ecinetic(:,ie_t) = max( Ecinetic(:,ie_t), Eclie )

      end do

      deallocate( rho, rhons, rs, Vh, Vhns, Vr, Vxc )

      Ecmax = 0._db
      do ie_t = 1,ne_t
        do isp = 1,nspin
! L'energie cinetique minimale est prise egale a celle du niveau de Fermi.
!          Ec_min = E_cut_optic - Workf - V0bdc_e(isp,ie_t)
!          Ec_min = max( Ec_min, 1.00_db / Rydb )
         Ec_min = 1.00_db / Rydb
!!!
          Ecinetic(isp,ie_t) = max( Ecinetic(isp,ie_t), Ec_min )
        end do
        do isp = 1,nspin
          Ecmax = max( Ecmax, Ecinetic(isp,ie_t) )
        end do
      end do

! En Optic sur 2 energies, la deuxieme a plus haute energie peut
! avoir un lmax plus grand (et inutile).
      call clmax(Ecmax,Rmtg(iprabs),lmaxat0,lmax,numat,lmaxfree)
      lmax = min(lmax,lmaxabs_t)

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        tp(2) = real(time,db)
      endif

! Calcul des tenseurs cartesiens

      nlma = nlm_probe
      if( Hubb_a .or. Full_potential ) then
        nlma2 = nlma
      else
        nlma2 = 1
      endif

      icheck_s = max( icheck(29), icheck(20) )

      do ie_t = 1,ne_t
        do ie_q = 2,nenerg_s
          if( Energ_s(ie_q) > En(ie_t) - eps10 ) exit
        end do
        p = ( En(ie_t) - Energ_s(ie_q - 1) ) / ( Energ_s(ie_q) - Energ_s(ie_q - 1) )
        do lm1 = 1,nlma
          do iso1 = 1,nspino
            lm1g = nlma*(iso1-1) + lm1
            do lm2 = 1,nlma
              do iso2 = 1,nspino
                lm2g = nlma*(iso2-1) + lm2
                do isp1 = 1,2
                  do isp2 = 1,2
                    if( Spinorbite ) then
                      iss1 = iso1; iss2 = iso2
                    else
                      if( isp1 /= isp2 ) cycle
                      iss1 = min(isp1,nspinp)
                      iss2 = min(isp2,nspinp)
                    endif

                    Taull_abs(lm1g,lm2g,isp1,isp2,ie_t,1,1) = p * Taull_opt(ie_q,lm1,lm2,iss1,iss2) &
                                                            + (1 - p ) * Taull_opt(ie_q-1,lm1,lm2,iss1,iss2)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      ndim1 = ne_t
      nenerg_tddft = 0
      allocate( rof0(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil) )

      call tenseur_car(Base_spin,coef_g,Core_resolved,Ecinetic, &
                Eimag(ie),Energ(ie),Enervide,Eseuil,Final_optic,Final_tddft,Full_potential,Green_i,Green_plus,Hubb_a,Hubb_d, &
                icheck_s,ie,ip_max,ip0,is_g,lmax_probe,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
                n_oo,nbseuil,ndim1,ndim2,nenerg_tddft,ninitl,ninitl,ninitlr,ne_t,ninitlr,nlm_pot,nlma,nlma2,nlmamax, &
                nr,nrm,nspin,nspino,nspinp,numat,psii,r,Relativiste,Rmtg(iprabs),Rmtsd,rof0,rot_atom_abs,Rot_int, &
                secdd_t,secdd_m_t,secdo_t,secdo_m_t,secdq_t,secdq_m_t,secmd_t,secmd_m_t,secmm_t,secmm_m_t,secoo_t,secoo_m_t, &
                secqq_t,secqq_m_t,Solsing,Solsing_only,Spinorbite,Taull_abs,Tddft,V_hubb,V_intmax,V0bdc_e,Vrato_e,Ylm_comp)

      deallocate( rof0 )

      if( E1E1 ) secdd(:,:,:,mpirank) = secdd(:,:,:,mpirank) + delta_E * secdd_t(:,:,:,mpirank)
      if( M1M1 ) secmm(:,:,:,mpirank) = secmm(:,:,:,mpirank) + delta_E * secmm_t(:,:,:,mpirank)
      if( E1M1 ) secmd(:,:,:,mpirank) = secmd(:,:,:,mpirank) + delta_E * secmd_t(:,:,:,mpirank)
      if( E1E2 ) secdq(:,:,:,:,mpirank) = secdq(:,:,:,:,mpirank) + delta_E * secdq_t(:,:,:,:,mpirank)
      if( E2E2 ) secqq(:,:,:,:,:,mpirank) = secqq(:,:,:,:,:,mpirank) + delta_E * secqq_t(:,:,:,:,:,mpirank)
      if( E1E3 ) secdo(:,:,:,:,:,mpirank) = secdo(:,:,:,:,:,mpirank) + delta_E * secdo_t(:,:,:,:,:,mpirank)
      if( E3E3 ) secoo(:,:,:,:,:,mpirank) = secoo(:,:,:,:,:,mpirank) + delta_E * secoo_t(:,:,:,:,:,mpirank)

      if( Green_int ) then
        if( E1E1 ) secdd_m(:,:,:,mpirank) = secdd_m(:,:,:,mpirank) + delta_E * secdd_m_t(:,:,:,mpirank)
        if( M1M1 ) secmm_m(:,:,:,mpirank) = secmm_m(:,:,:,mpirank) + delta_E * secmm_m_t(:,:,:,mpirank)
        if( E1M1 ) secmd_m(:,:,:,mpirank) = secmd_m(:,:,:,mpirank) + delta_E * secmd_m_t(:,:,:,mpirank)
        if( E1E2 ) secdq_m(:,:,:,:,mpirank)=secdq_m(:,:,:,:,mpirank) + delta_E * secdq_m_t(:,:,:,:,mpirank)
        if( E2E2 ) secqq_m(:,:,:,:,:,mpirank) = secqq_m(:,:,:,:,:,mpirank) + delta_E * secqq_m_t(:,:,:,:,:,mpirank)
        if( E1E3 ) secdo_m(:,:,:,:,:,mpirank) = secdo_m(:,:,:,:,:,mpirank) + delta_E * secdo_m_t(:,:,:,:,:,mpirank)
        if( E3E3 ) secoo_m(:,:,:,:,:,mpirank) = secoo_m(:,:,:,:,:,mpirank) + delta_E * secoo_m_t(:,:,:,:,:,mpirank)
        endif

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        tp(3) = real(time,db)
        tpt(3) = tpt(3) + tp(2) - tp(1)
        tpt(8) = tpt(8) + tp(3) - tp(2)
      endif

    end do ! fin boucle sur energie des electrons

    endif  ! arrive ie > nenerg

    if( mpinodes > 1 ) then

      call MPI_RECV_all(MPI_host_num_for_mumps,mpinodes,mpirank,mpirank0,Multipole,n_oo,ninitlr,secdd,secdo, &
                       secdq,secmd,secmm,secoo,secqq)
      if( Green_i ) call MPI_RECV_all(MPI_host_num_for_mumps,mpinodes,mpirank,mpirank0,Multipole,n_oo,ninitlr,secdd_m,secdo_m, &
                                      secdq_m,secmd_m,secmm_m,secoo_m,secqq_m)

    endif

    if( mpirank0 /= 0 ) cycle

    icheck_s = max( icheck(29), icheck(21) )

    do ie_computer = 0,mpinodes-1

      ie = ( je - 1 ) * mpinodes + ie_computer + 1

      if( ie > nenerg ) exit

      call write_coabs(Allsite,angxyz,axyz,Base_spin, Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
            Densite_atom,E_cut_optic,Energ,Energphot,Extract,Epsii,Eseuil,Final_tddft, &
            f_avantseuil,Full_self_abs,Green_i,Green_plus,hkl_dafs,iabsorig,icheck_s,ie,ie_computer, &
            Int_tens,isigpi,isymeq,jseuil,ltypcal,Moyenne,mpinodes,Multipole,n_multi_run,n_oo,natomsym,nbseuil, &
            ncolm,ncolr,ncolt,nenerg,ninit1,ninitlr,nomabs,nomfich,nomfich_cal_convt,nomfich_s,nphi_dafs,npldafs, &
            nphim,nplr,nplrm,nseuil,nspinp,numat,nxanout,pdp,phdafs,phdf0t, &
            phdt,pol,Polarise,poldafse,poldafss,Rot_int,sec_atom,secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, &
            secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Self_abs,Spherical_signal, &
            Spherical_tensor,Spinorbite,Taux_eq,V0muf,vecdafse,vecdafss,vec,Volume_maille,Xan_atom)

    end do

    if( mpirank0 == 0 ) then
      call CPU_TIME(time)
      tp(4) = real(time,db)
      tpt(10) = tpt(10) + tp(4) - tp(3)
    endif

  end do boucle_energ

  deallocate( Ecinetic )
  deallocate( Ecinetic_e )
  deallocate( Energ )
  deallocate( Eimag )
  deallocate( V0bdc_e )
  deallocate( Vrato, Vrato_e )

  return
  100 format(/1x,120('-')//,' Cycle Optic')
end

!***********************************************************************

! Sous-programme qui definit les dimensions de la nouvelle gamme d'energie pour le calcul optic.

subroutine dim_grille_optic(E_cut,Energ_s,mpirank0,nenerg,nenerg_s)

  use declarations
  implicit none

  integer,intent(in):: mpirank0, nenerg_s
  integer,intent(out):: nenerg

  real(kind=db), intent(in):: E_cut
  real(kind=db), dimension(nenerg_s), intent(in):: Energ_s

  integer:: ie, ipr, ne

  if( E_cut < Energ_s(1) .and. mpirank0 == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,110) E_cut*rydb, Energ_s(1)*rydb
    end do
    stop
  endif

  if( E_cut > Energ_s(nenerg_s) .and. mpirank0 == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,120) E_cut*rydb, Energ_s(nenerg_s)*rydb
    end do
    stop
  endif

  do ie = 1,nenerg_s
    if( Energ_s(ie) > E_cut + eps10 ) exit
  end do
  ne = ie - 1

  nenerg = min( nenerg_s - ne, ne )

  return
  110 format(///' E_Fermi =',f7.3,' eV < First energy =',f7.3,' eV',// ' Increase the energy range of calculation !'//)
  120 format(///' E_Fermi =',f7.3,' eV > Last energy =',f7.3,' eV',// ' Increase the energy range of calculation !'//)
end

!***********************************************************************

! Sous-programme qui definit une nouvelle gamme d'energie pour le calcul optic.

subroutine grille_optic(E_cut,Energ,Energ_s,icheck,nenerg, nenerg_s)

  use declarations
  implicit none

  integer,intent(in):: icheck, nenerg, nenerg_s

  real(kind=db), intent(in):: E_cut
  real(kind=db), dimension(nenerg_s), intent(in):: Energ_s
  real(kind=db), dimension(nenerg), intent(out):: Energ

  integer:: ie, ne

  do ie = 1,nenerg_s
    if( Energ_s(ie) > E_cut + eps10 ) exit
  end do
  ne = ie - 1

  do ie = 1,nenerg
    Energ(ie) = Energ_s(ie+ne) - E_cut
  end do

  if( icheck > 2 ) then
    write(3,100)
    write(3,110)
    write(3,120) Energ(:)*rydb
  end if

  return
  100 format(/' ---- Grille_optic -------',100('-'))
  110 format(/'The energy grid for the optic calculation:')
  120 format(5f13.7)
end

