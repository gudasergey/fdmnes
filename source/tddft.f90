! FDMNES subroutine

! Subset of routines for TDDFT procedure.
! For each energy the kernel K is calculated.

!***********************************************************************

! nbseuil : nr of edges
! ninitl  : nbr of core states
! ninitl_out = ninitl in XANES core_resolved
!            = nbseuil in XANES not core_resolved
!            = number of l in optic core_resolved
!            = 1 in optic not core_resolved
! ninitlv = nbseuil in DFT
!         = ninitl_out in TDDFT
! n_Ec = nbr of electron energies a the photon energy
!       = 1 in DFT
!       = ninitlv in TDDFT
!       = 2 in TDDFT + Optic

subroutine main_tddft(alfpot,All_nrixs,angxyz,Allsite,Atomic_scr,axyz,Bulk_step,Classic_irreg, &
        coef_g,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio,Delta_edge,Delta_Eseuil,Dipmag, &
        dv0bdcF,Dyn_eg,Dyn_g,E_cut,E_cut_imp,E_Fermi,E_cut_man,Ecent,Eclie,Elarg,Eneg, &
        Energ_t,Energphot,Epsii_a,Epsii_moy,Eseuil,Estart,f_avantseuil,Full_potential,Full_self_abs, &
        Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,hkl_dafs,Hubb_a,Hubb_d,icheck, &
        iabsorig,iopsymc_25,is_g,isigpi,isymeq,jseuil,Kern_fac, &
        l0_nrixs,ldip,Length_abs,lmax_pot,lmax_nrixs,lmaxabs_t,lmaxat0,lmaxfree,lmoins1,loct,lplus1,lqua,lseuil, &
        ltypcal,m_g,m_hubb,Magnetic,Matper,Moyenne,mpinodes,mpirank,mpirank0,msymdd,msymddi, &
        msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
        n_multi_run,n_oo,n_rel,n_rout,n_tens_max, &
        natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,nenerg_tddft,ngamh,ninit1,ninitl,ninitl_out,ninitlv,nlm_pot,nlmamax, &
        nomabs,nomfich,nomfich_cal_tddft_conv,nomfich_s, &
        nphi_dafs,nphim,npldafs,nplr,nplrm,nq_nrixs,nr,NRIXS,nrm,nseuil,nspin,nspino,nspinp, &
        numat,nxanout,Octupole,Old_zero,pdp,phdafs,phdf0t,phdt,pol,poldafse,poldafss, &
        psii,q_nrixs,Quadrupole,r,Relativiste,Renorm,rhoato_abs,Rmtg,Rmtsd, &
        rof0,rot_atom_abs,Rot_int,RPALF,rsato,rsbdc,Self_abs,Solsing_only, &
        Spherical_signal,Spherical_tensor,Spinorbite,Taull_tdd,Taux_eq,Time_rout,V_intmax,V_hubb,V0muf, &
        Vcato,vec,vecdafse,vecdafss,Vhbdc,Volume_maille,VxcbdcF, &
        Vxcato,Workf,xsect_file,Ylm_comp_e)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: cal_nenerge, iabsorig, icheck_s, ie, ie_computer, ie_e, initl, iopsymc_25, ip_max, ip0, &
    ir, isp, je, jseuil, l, l0_nrixs, lmax, lmax_pot, lmax_probe, lmax_nrixs, lmaxabs_t, &
    lmaxat0, lseuil,m_hubb, mpinodes, mpirank, mpirank0, &
    n_Ec, n_multi_run, n_oo, n_rel, n_rout, n_tens_max, n_V, natomsym, nbseuil, &
    ncolm, ncolr, ncolt, nd3, nenerg, nenerg_s, nenerg_tddft, nenerge, ngamh, nge, ninit1, ninitl, ninitl_out, &
    ninitlv, nlm, nlm_fp, nlm_pot, nlm_probe, nlm_p_fp, nlmamax, nlms_f, nlms_g, nlmsm_f, &
    nphim, npldafs, nplr, nplrm, nq_nrixs, nr, nrm, ns_dipmag, &
    nseuil, nspin, nspino, nspinp, numat, nxanout

  integer, dimension(30):: icheck
  integer, dimension(natomsym):: isymeq
  integer, dimension(ninitl):: is_g
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(ninitl,2):: m_g
  integer, dimension(2,npldafs):: isigpi
  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3)::  msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi

  character(len=132):: nomfich, nomfich_s, nomfich_cal_convt, nomfich_cal_tddft_conv, xsect_file
  character(len=13), dimension(nplrm):: ltypcal
  character(len=length_word), dimension(ncolm):: nomabs

  complex(kind=db):: f_avantseuil
  complex(kind=db), dimension(3,nplrm):: pol
  complex(kind=db), dimension(natomsym,npldafs):: phdafs
  complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdt
  complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss
  complex(kind=db), dimension(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil):: rof0
  complex(kind=db), dimension(3,3,ninitl_out,0:mpinodes-1):: secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,n_rel,ninitl_out,0:mpinodes-1):: secdd, secdd_m
  complex(kind=db), dimension(3,3,3,ninitl_out,0:mpinodes-1):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitl_out,0:mpinodes-1):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitl_out,0:mpinodes-1):: secoo, secoo_m
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nlmamax,nspinp) :: Taull_tdd

  complex(kind=db), dimension(:,:), allocatable:: Trans, V
  complex(kind=db), dimension(:,:,:,:), allocatable:: rof_ph, V_hubb_t
  complex(kind=db), dimension(:,:,:,:,:), allocatable::  rof0_e
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: Chi_0
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: Chi

  logical:: All_nrixs, Allsite, Atomic_scr, Bulk_step, Cartesian_tensor, Classic_irreg, Core_resolved, Dafs, Dafs_bio, Dipmag, &
    Dyn_eg, Dyn_g, E_cut_man, Eneg, Energphot, FDM_comp, Final_optic, Final_tddft, First_E, &
    Full_potential, Full_self_abs, Gamma_hole_imp, Gamma_tddft, Green, Green_int, &
    Hubb_a, Hubb_d, lmaxfree, lmoins1, lplus1, Magnetic, Matper, &
    Moyenne, NRIXS, Octupole, Old_zero, Optic, Quadrupole, Radial_comp, &
    Relativiste, Renorm, RPALF, Self_abs, Solsing, Solsing_only, Spherical_signal, &
    Spherical_tensor, Spinorbite, Xan_atom, Ylm_comp, Ylm_comp_e

  logical, dimension(10):: Multipole

  real(kind=sg) time

  real(kind=db):: alfpot, Delta_edge, Delta_Eseuil, E_cut, E_cut_imp, E_cut_tddft, E_fermi, Eclie, &
     Ecent, Ecmax, EFermi_min, Elarg, Energ_tt, Enervide_t, Epsii_moy, Estart, &
     Gamma_max, Kern_fac, p, Rmtg, Rmtsd, rsbdc, V_intmax, V0muf, Vhbdc, Volume_maille, Workf

  real(kind=db), dimension(3):: angxyz, axyz
  real(kind=db), dimension(8):: Time_loc
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(n_rout):: Time_rout
  real(kind=db), dimension(nspin):: dv0bdcF, V0bdc_t, VxcbdcF
  real(kind=db), dimension(nenerg_s):: Energ_s, Energ_t
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitlv):: Epsii_a
  real(kind=db), dimension(ninitl_out):: Epsii, sec_atom
  real(kind=db), dimension(nr):: r, rsato
  real(kind=db), dimension(natomsym) :: Taux_eq
  real(kind=db), dimension(3,3):: rot_atom_abs, Rot_int
  real(kind=db), dimension(nplrm,2) :: pdp
  real(kind=db), dimension(3,nplrm) :: vec
  real(kind=db), dimension(ninitl,2):: coef_g
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nspin):: rhoato_abs
  real(kind=db), dimension(n_tens_max*ninitl_out,0:natomsym):: Int_tens
  real(kind=db), dimension(nr,nlm_pot):: Vcato
  real(kind=db), dimension(nr,2,2):: fxc
  real(kind=db), dimension(3,npldafs,nphim):: Vecdafse, Vecdafss
  real(kind=db), dimension(3,npldafs):: hkl_dafs
  real(kind=db), dimension(npldafs):: Length_abs
  real(kind=db), dimension(nr,nlm_pot,nspin) :: Vxcato
  real(kind=db), dimension(nq_nrixs):: q_nrixs
  real(kind=db), dimension(nq_nrixs,ninitl_out,0:mpinodes-1):: S_nrixs, S_nrixs_m
  real(kind=db), dimension(nq_nrixs,l0_nrixs:lmax_nrixs,ninitl_out,0:mpinodes-1):: S_nrixs_l, S_nrixs_l_m

  real(kind=db), dimension(:), allocatable:: Decal_initl, Ecinetic_e, Eimag, Energ, Enervide, rst, Vct, Vrt, Vxct
  real(kind=db), dimension(:,:), allocatable:: Ecinetic, V0bdc
  real(kind=db), dimension(:,:,:), allocatable::  Vrato_t
  real(kind=db), dimension(:,:,:,:), allocatable:: Vrato
  real(kind=db), dimension(:,:,:,:,:,:), allocatable:: fppn, zet
  real(kind=db), dimension(:,:,:,:,:,:,:), allocatable::  Kern

  if( icheck(1) > 0 ) write(3,100)

! Green complet avec partie reelle et imaginaire (en cas de Gamma non nul)
! Si Green_int, f' et f" sont directement calcule et l'integrale sur l'energie n'est pas a faire dans convolution.
  Green_int = Gamma_tddft

  Xan_atom = .false.
  Optic = .false.
  FDM_comp = .false.
  Final_optic = Optic
  Final_tddft = .true.
  Green = .true.

  if( Dipmag ) then
    ns_dipmag = 2
    ip0 = 0
  else
    ns_dipmag = 1
    ip0 = 1
  endif
  if( Octupole ) then
    ip_max = 3
  else
    ip_max = 2
  endif

  sec_atom(:) = 0._db

! Tous les calculs sont imposes en harmoniques complexes
! S'ils n'etaient pas en complexes, ils ne pouvaient pas etre en Spinorbite.
  Ylm_comp = .true.
  if( .not. Ylm_comp_e ) then
    Call Tr_Taull( icheck(22), nlmamax, nenerg_s, nspinp, Taull_tdd )

    if( Hubb_a ) then
      l = m_hubb
      allocate( V(-l:l,-l:l) )
      allocate( Trans(-l:l,-l:l) )
      Call Cal_trans_l(l,Trans)
      do isp = 1,nspinp
        V(-l:l,-l:l) = V_hubb(-l:l,-l:l,isp,isp)
        V = Matmul( Conjg( Transpose(Trans) ), Matmul( V, Trans ))
        V_hubb(-l:l,-l:l,isp,isp) = V(-l:l,-l:l)
      end do
      deallocate( Trans, V )
    endif
  endif

  ! Les Taull_tdd sont les amplitudes "efficaces". On multiplie par "i" pour que que la partie imaginaire soit l'absorption
  Taull_tdd(:,:,:,:,:) = img * Taull_tdd(:,:,:,:,:)

  allocate( V_hubb_t(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp) )
  V_hubb_t(:,:,:,:) = V_hubb(:,:,:,:)

  if( E_cut_man ) then
    E_cut_tddft = E_cut_imp
  else
    E_cut_tddft = E_cut
  endif

  Energ_s(:) = Energ_t(:)
  n_Ec = ninitlv
  n_V = ninitlv

  allocate( Ecinetic_e(nspin) )
  allocate( V0bdc(nspin,n_V) )

  allocate( Decal_initl(n_Ec) )
! Epsii_moy est la moyenne pour le 2eme seuil (si 2 seuils) par exemple L3
  Decal_initl(:) = Epsii_a(:) - Epsii_moy
  Epsii(1) = Epsii_moy

! Elaboration de la grille etendue pour la Tddft
  call dim_grille_tddft(Energ_s,Delta_Eseuil,Estart,nbseuil,nenerg,nenerg_s)
  allocate( Energ(nenerg) )
  allocate( Eimag(nenerg) )
  call grille_tddft(Energ,Energ_s,Delta_Eseuil,Estart,icheck(22),nbseuil,nenerg,nenerg_s)

! On garde l'energie imaginaire nulle
  Eimag(:) = 0._db
  Solsing = .false.

  if( Spinorbite ) then
    nlmsm_f = ( lmaxabs_t + 1 ) * ( 2 * lmaxabs_t + 1 )
  else
    nlmsm_f = ( lmaxabs_t + 1 )**2
  endif

  nlms_g = 0
  do initl = 1,ninitl
    if( abs( Coef_g(initl,1) ) > eps10 ) nlms_g = nlms_g + 1
  end do
  ! En xanes, la partie radiale change seulement entre les 2 seuils
  nd3 = ninitl    ! pour Chi envoye dans tenseur

  if( Octupole ) then
    lmax_probe = lseuil + 3
  elseif( Quadrupole ) then
    lmax_probe = lseuil + 2
  else
    lmax_probe = lseuil + 1
  endif

  nlm_probe = (lmax_probe + 1)**2
  if( Hubb_a .or. Full_potential ) then
    nlm_p_fp = nlm_probe
  else
    nlm_p_fp = 1
  endif

  allocate( Vrato(nr,nlm_pot,nspin,n_V) )

  if( alfpot > eps6 ) then

    do ie_e = 1,n_V
      do isp = 1,nspinp
        Vrato(1:nr,:,isp,ie_e) = Vcato(1:nr,:) + Vxcato(1:nr,:,isp)
        V0bdc(isp,ie_e) = Vhbdc + VxcbdcF(isp) + dV0bdcF(isp)
      end do
    end do

  endif

  nenerge = cal_nenerge(Energ_s,Eseuil(nbseuil),nenerg_s)
  allocate( fppn(nenerg_s-1:nenerge,nlmsm_f,nlmsm_f,nbseuil,nbseuil,nspinp) )
  call Extrap_fpp( Core_resolved, decal_initl, Delta_edge, Energ_s, Eseuil(nbseuil), icheck(23), E_cut_tddft, EFermi_min, fppn, &
                     lmaxabs_t, lseuil, nbseuil, n_Ec, nenerg_s, nenerge, ninit1, nlmamax, nlmsm_f, nspino, &
                     nspinp, numat, rof0, Spinorbite, Taull_tdd, xsect_file )
  E_cut = EFermi_min

! Valeur eventuellement decalee vers le bas pour ne rien couper a la convolution.

  allocate( Chi(nlm_probe*nspino,nlm_probe*nspino,nd3,nd3,2,2,ns_dipmag) )
  allocate( Ecinetic(nspin,n_Ec) )
  allocate( Enervide(n_Ec) )

 ! Calcul de la fonctionnelle d'echange-correlation
  if( .not. RPALF ) call fxcorr(alfpot,fxc,icheck(24),Magnetic,nr,nspin,r,rhoato_abs,rsato)

  nge = ( nenerg - 1 ) / mpinodes + 1

! Boucle sur l'energie des photons ---------------------------------------------------------------------------------------

  First_E = .true.

  boucle_energ: do je = 1,nge

    if( mpirank0 == 0 ) then
      call CPU_TIME(time)
      Time_loc(1) = real(time,db)
    endif

    ie = ( je - 1 ) * mpinodes + mpirank + 1

    if( ie <= nenerg ) then

      do ie_e = 1,n_Ec
        if( Old_zero ) then
          Enervide(ie_e) = Energ(ie) - Decal_initl(ie_e) - Workf
        else
          Enervide(ie_e) = Energ(ie) - Decal_initl(ie_e) + E_Fermi
        endif
      end do

! Calcul du potentiel excite
      if( alfpot < eps6 )  then

        allocate( rst(nr) )
        allocate( Vct(nr) )
        allocate( Vrt(nr) )
        allocate( Vxct(nr) )
        rst(1:nr) = rsato(1:nr)
        Vct(1:nr) = Vcato(1:nr,1)
        do isp = 1,nspin
          Vxct(1:nr) = Vxcato(1:nr,1,isp)
          do ie_e = 1,n_Ec
            call subpotex(nr,Vrt,Vct,Vxct,rst,Enervide(ie_e))
            Vrato(1:nr,1,isp,ie_e) = Vrt(1:nr)
            do ir = 1,nr
              p = ( Vrt(ir) - Vcato(ir,1) ) / Vxct(ir)
              Vrato(ir,2:nlm_pot,isp,ie_e) = Vcato(ir,2:nlm_pot) + p * Vxcato(ir,2:nlm_pot,isp)
            end do
          end do
        end do
        deallocate( rst, Vct, Vrt, Vxct )

        allocate( rst(1) )
        allocate( Vct(1) )
        allocate( Vrt(1) )
        allocate( Vxct(1) )
        rst(1) = rsbdc
        Vct(1) = Vhbdc
        do isp = 1,nspin
          Vxct(1) = VxcbdcF(isp)
          do ie_e = 1,n_Ec
            call subpotex(1,Vrt,Vct,Vxct,rst,Enervide(ie_e))
            V0bdc(isp,ie_e) = Vrt(1) + dV0bdcF(isp)
          end do
        end do
        deallocate( rst, Vct, Vrt, Vxct )

      endif

      do ie_e = 1,n_Ec
        Ecinetic(:,ie_e) = Enervide(ie_e) - V0bdc(:,min(n_V,ie_e))
        if( .not. Eneg ) Ecinetic(:,ie_e) = max( Ecinetic(:,ie_e), Eclie )
      end do

      Ecmax = 0._db
      do ie_e = 1,n_Ec
        do isp = 1,nspin
          Ecmax = max( Ecmax, Ecinetic(isp,ie_e) )
        end do
      end do
      call clmax(Ecmax,Rmtg,lmaxat0,lmax,numat,lmaxfree)
      lmax = min(lmax,lmaxabs_t)

      nlm = ( lmax + 1 )**2
      if( Full_potential .or. ( Hubb_a .and. .not. Hubb_d ) ) then
        nlm_fp = nlm
      else
        nlm_fp = 1
      endif

      if( Spinorbite ) then
        nlms_f = ( lmax + 1 ) * ( 2 * lmax + 1 )
      else
        nlms_f = ( lmax + 1 )**2
      endif

! On ne garde que la partie reelle de l'orbitale radiale zet. Verifier plus tard si correct
      Radial_comp = Eimag(ie) > eps10 .or. ( Hubb_a .and. Ylm_comp )

      allocate( zet(nr,nlm,nlm_fp,nspinp,nspino,n_Ec) )
      zet(:,:,:,:,:,:) = 0._db

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        Time_loc(2) = real(time,db)
      endif

      icheck_s = max( icheck(18), icheck(22) )

      allocate( rof_ph(nlm,nspinp,nspino,n_Ec) )
      rof_ph(:,:,:,:) = ( 0._db, 0._db )

      do ie_e = 1,n_Ec

        allocate( Vrato_t(nr,nlm_pot,nspin) )

        Ecinetic_e(:) = Ecinetic(:,ie_e)

        Energ_tt = Energ(ie) - Decal_initl(ie_e)
        Enervide_t = Enervide(ie_e)

        Vrato_t(1:nr,:,:) = Vrato(1:nr,:,:,ie_e)
        V0bdc_t(:) = V0bdc(:,ie_e)

        call Radial_wave(Ecinetic_e,Eimag(ie),Energ_tt,Enervide_t,Full_potential,Hubb_a,Hubb_d,icheck_s,ie_e,lmax,lmax_pot, &
          m_hubb,nbseuil,ninit1,n_Ec,nlm_pot,nlm,nlm_fp,nr,nrm,nspin,nspino,nspinp,numat,psii,r,Radial_comp, &
          Relativiste,Renorm,Rmtg,Rmtsd,rof_ph,Spinorbite,V_hubb_t,V_intmax,V0bdc_t,Vrato_t,Ylm_comp,zet)

        deallocate( Vrato_t )

      end do

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        Time_loc(3) = real(time,db)
      endif

      allocate( Kern(nlms_f,nlms_f,nlms_g,nlms_g,2,2,ns_dipmag) )
      allocate( Chi_0(nlms_f,nlms_f,nlms_g,nlms_g,2,ns_dipmag) )

! Calcul de Chi_0
      call Chi_0_int(Chi_0,Coef_g,Core_resolved,Decal_initl,Delta_edge,E_cut_tddft,Ecent,Elarg,Energ(ie),Energ_s, &
         Eseuil(nbseuil),fppn,Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,icheck(23),ie,jseuil,lmax, &
         nbseuil,nenerge,nenerg_s,ngamh,ninit1,ninitl,n_Ec,nlm,nlmamax,nlmsm_f,nlms_f,nlms_g,nomfich,ns_dipmag, &
         nseuil,nspino,nspinp,numat,rof_ph,rof0,Spinorbite,Taull_tdd)
      deallocate( rof_ph )

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        Time_loc(4) = real(time,db)
      endif

! Kernel calculation
      call kernel(Atomic_scr,coef_g,Core_resolved,Dipmag,Dyn_eg,Dyn_g,Energ(ie),fxc,icheck(24),ie,Kern,Kern_fac, &
                lmax,lseuil,m_g,nbseuil,ninit1,ninitl,n_Ec,nlm,nlm_fp,nlms_f,nlms_g,nomfich,nr,nrm, &
                ns_dipmag,nspino,nspinp,psii,r,Rmtsd,RPALF,Spinorbite,zet)

      deallocate( zet )

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        Time_loc(5) = real(time,db)
      endif

      call Cal_Chi(Chi, Chi_0, coef_g, Energ(ie), First_E, icheck(25), iopsymc_25, Kern, lmax_probe, lmax, &
              lseuil, nd3, ninitl, nlm_probe, nlm_fp, nlms_f, nlms_g, &
              nomfich, NRIXS, ns_dipmag, nspino, Optic, Quadrupole, Spinorbite)

      deallocate( Chi_0, Kern )

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        Time_loc(6) = real(time,db)
      endif

! Calcul des tenseurs cartesiens

      nenerg_tddft = 0

      allocate( rof0_e(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil) )

      icheck_s = icheck(20)

      if( NRIXS ) call S_nrixs_cal(Classic_irreg,coef_g,Core_resolved,Ecinetic, &
                      Eimag(ie),Energ(ie),Enervide,Eseuil,FDM_comp,Final_tddft,Full_potential,Green,Green_int,Hubb_a,Hubb_d, &
                      icheck_s,l0_nrixs,lmax_nrixs,is_g,lmax_probe,lmax_pot,lmoins1,lplus1,lseuil,m_g,m_hubb, &
                      mpinodes,mpirank, &
                      n_Ec,n_V,nbseuil,ns_dipmag,nd3,ninit1,ninitl,ninitl_out,ninitlv,nlm_pot,nlm_probe, &
                      nlm_p_fp,nq_nrixs,nr,nrm,nspin,nspino,nspinp,numat,psii,q_nrixs,r,Relativiste,Renorm,Rmtg, &
                      Rmtsd,S_nrixs,S_nrixs_l,S_nrixs_l_m,S_nrixs_m,Solsing,Solsing_only,Spinorbite,Chi, &
                      V_hubb_t,V_intmax,V0bdc,Vrato,Ylm_comp)

      call tenseur_car(Classic_irreg,coef_g,Core_resolved,Ecinetic, &
                Eimag(ie),Energ(ie),Enervide,Eseuil,FDM_comp,Final_optic,Final_tddft,Full_potential,Green,Green_int,Hubb_a, &
                Hubb_d,icheck_s,ie,ip_max,ip0,is_g,lmax_probe,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
                n_Ec,n_oo,n_rel,n_V,nbseuil,ns_dipmag,nd3,nenerg_tddft,ninit1,ninitl,ninitl_out,ninitlv,nlm_pot,nlm_probe, &
                nlm_p_fp,nlmamax,nr,nrm,nspin,nspino,nspinp,numat,psii,r,Relativiste,Renorm,Rmtg,Rmtsd,rof0_e,rot_atom_abs, &
                Rot_int,secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m, &
                secqq,secqq_m,Solsing,Solsing_only,Spinorbite,Chi,.true.,V_hubb_t,V_intmax,V0bdc,Vrato,Ylm_comp)

      deallocate( rof0_e )

    endif ! arriving when ie > nenerg

    if( mpirank0 == 0 ) then
      call CPU_TIME(time)
      Time_loc(7) = real(time,db)
    endif

    if( mpinodes > 1 ) then

      call MPI_RECV_all(l0_nrixs,lmax_nrixs,mpinodes,mpirank,mpirank0,Multipole,n_oo,n_rel,ninitl_out, &
                        nq_nrixs,S_nrixs,S_nrixs_l,secdd,secdo,secdq,secmd,secmm,secoo,secqq)
      if( Green_int ) call MPI_RECV_all(l0_nrixs,lmax_nrixs,mpinodes,mpirank,mpirank0,Multipole,n_oo, &
                        n_rel,ninitl_out,nq_nrixs,S_nrixs_m,S_nrixs_l_m,secdd_m,secdo_m,secdq_m,secmd_m,secmm_m,secoo_m,secqq_m)

    endif

    if( mpirank0 /= 0 ) cycle

    do ie_computer = 0,mpinodes-1

      ie = ( je - 1 ) * mpinodes + ie_computer + 1

      if( ie > nenerg ) exit

      call Write_coabs(Allsite,angxyz,axyz,Bulk_step,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
            E_cut,Energ,Energphot,.false.,Epsii,Eseuil,Final_tddft,First_E, &
            f_avantseuil,Full_self_abs,Green_int,hkl_dafs,iabsorig,icheck(21),ie,ie_computer,Int_tens, &
            isigpi,isymeq,jseuil,Length_abs,ltypcal,Matper,Moyenne,mpinodes,Multipole,n_multi_run,n_oo,n_rel,n_tens_max, &
            natomsym,nbseuil, &
            ncolm,ncolr,ncolt,nenerg,ninit1,ninitl_out,nomabs,nomfich,nomfich_cal_convt,nomfich_s,nphi_dafs,npldafs, &
            nphim,nplr,nplrm,nseuil,nspinp,numat,nxanout,pdp,phdafs,phdf0t, &
            phdt,pol,poldafse,poldafss,sec_atom,secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, &
            secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Self_abs,Spherical_signal, &
            Spherical_tensor,Spinorbite,Taux_eq,V0muf,Vecdafse,Vecdafss,Vec,Volume_maille,Xan_atom)

      if( NRIXS ) call Write_nrixs(All_nrixs,Allsite,Core_resolved,Volume_maille, &
                  E_cut,Energ,Energphot,.false.,Epsii,Eseuil,Final_tddft,First_E, &
                  f_avantseuil,Green_int,iabsorig,icheck(21),ie,ie_computer,l0_nrixs,lmax_nrixs,isymeq, &
                  jseuil,mpinodes,n_multi_run,natomsym,nbseuil,nenerg,ninit1,ninitl_out,nomfich,nomfich_cal_convt, &
                  nq_nrixs,nseuil,nspinp,numat,q_nrixs,S_nrixs,S_nrixs_l,S_nrixs_l_m,S_nrixs_m,Spinorbite,Taux_eq,V0muf)

      if( First_E ) nomfich_cal_tddft_conv = nomfich_cal_convt
      First_E = .false.

    end do

    if( mpirank0 == 0 ) then
      call CPU_TIME(time)
      Time_loc(8) = real(time,db)
      Time_rout(6) = Time_rout(6) + Time_loc(2) - Time_loc(1)  ! Potex

      Time_rout(14) = Time_rout(14) + Time_loc(3) - Time_loc(2)
      Time_rout(15) = Time_rout(15) + Time_loc(4) - Time_loc(3)
      Time_rout(16) = Time_rout(16) + Time_loc(5) - Time_loc(4)
      Time_rout(17) = Time_rout(17) + Time_loc(6) - Time_loc(5)

      Time_rout(11) = Time_rout(11) + Time_loc(7) - Time_loc(6) ! Tenseur
      Time_rout(13) = Time_rout(13) + Time_loc(8) - Time_loc(7) ! Coabs
    endif

  end do boucle_energ   ! Fin de la boucle sur l'energie.

  deallocate( Decal_initl )
  deallocate( fppn )
  deallocate( Chi )
  deallocate( Ecinetic, Enervide )
  deallocate( Ecinetic_e )
  deallocate( Energ )
  deallocate( Eimag )
  deallocate( V0bdc )
  deallocate( V_hubb_t )
  deallocate( Vrato )

  return
  100 format(/1x,120('-')//,' Cycle Tddft')
end

!***********************************************************************

! Sous-programme qui definit les dimensions de la nouvelle gamme
! d'energie pour le calcul tddft. Sert lorsqu'on a deux seuils qui se recouvrent

subroutine dim_grille_tddft(Energ_s,Delta_Eseuil,Estart,nbseuil,nenerg,nenerg_s)

  use declarations
  implicit none

  integer,intent(in):: nbseuil, nenerg_s
  integer,intent(out):: nenerg

  real(kind=db),intent(in):: Delta_Eseuil, Estart
  real(kind=db),dimension(nenerg_s),intent(in):: Energ_s

  integer ie, je0

  real(kind=db):: Pasdeb

  if( nbseuil == 1 ) then

    nenerg = nenerg_s

  else

    do ie = 1,nenerg_s
! Fonctionne meme si on est en energie de photon
      if( Energ_s(ie) - Energ_s(1) > Delta_Eseuil - eps10 ) exit
    end do
    je0 = ie - 1
    do ie = 1,nenerg_s
      if(  Delta_Eseuil + Energ_s(ie) > Energ_s(nenerg_s) ) exit
    end do
    nenerg = je0 + ie

  end if

  if( Energ_s(1) > Estart - 1.e-10_db ) then
    pasdeb = 0.5_db / rydb
    nenerg = nenerg + nint( ( Energ_s(1) - Estart ) / pasdeb )
  endif

  return
end

!***********************************************************************

! Sous-programme qui definit une nouvelle gamme d'energie pour le calcul
! tddft; elle sert lorsque deux seuils se recouvrent

subroutine grille_tddft(Energ,Energ_s,Delta_Eseuil,Estart,icheck, nbseuil,nenerg,nenerg_s)

  use declarations
  implicit none

  integer,intent(in):: icheck, nbseuil, nenerg, nenerg_s

  real(kind=db),intent(in):: Delta_Eseuil, Estart
  real(kind=db),dimension(nenerg_s),intent(in):: Energ_s

  real(kind=db),dimension(nenerg),intent(out):: Energ

  integer ie, je0, n

  real(kind=db):: Pasdeb

  if( Energ_s(1) > Estart - 1.e-10_db ) then
    pasdeb = 0.5_db / rydb
    n = nint( ( Energ_s(1) - Estart ) / pasdeb )
    do ie = 1,n
      Energ(ie) = Energ_s(1) + ( ie - n - 1 ) * Pasdeb
    end do
  else
    n = 0
  endif

  if( nbseuil == 1 ) then

    Energ(n+1:n+nenerg_s) = Energ_s(1:nenerg_s)

  else

    do ie = 1,nenerg_s
      Energ(n+ie) = Energ_s(ie)
! Fonctionne meme si on est en energie de photon
      if( Energ_s(ie) - Energ_s(1) > Delta_Eseuil - eps10 ) exit
    end do
    je0 = n + ie - 1
    do ie = 1,nenerg_s
      Energ( je0 + ie ) = Delta_Eseuil + Energ_s(ie)
      if(  Delta_Eseuil + Energ_s(ie) > Energ_s(nenerg_s) ) exit
    end do
  end if

  if( icheck > 2 ) then
    write(3,100)
    write(3,110)
    write(3,120) Energ(:)*rydb
  end if

  return
  100 format(/' ---- Grille_tddft -------',100('-'))
  110 format(/'The energy grid for the tddft calculation:')
  120 format(5f13.7)
end

!***************************************************************

function Cal_nenerge(Energ_s,Eseuil,nenerg_s)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: Cal_nenerge, ie, nenerg_s

  real(kind=db):: alfa, dde, def, Ephm, Ephoton, Eseuil
  real(kind=db), dimension(nenerg_s), intent(in):: Energ_s

  Ephm = Eseuil + 10000._db / rydb
  if( nenerg_s > 1 ) then
    def = Energ_s(nenerg_s) - Energ_s(nenerg_s-1)
  else
    def = 2 / rydb
  endif

  alfa = 1.02_db
  dde = def

  Ephoton = Eseuil + Energ_s(nenerg_s)
  do ie = 1,1000000
    dde = alfa * dde
    Ephoton = Ephoton + dde
    if( Ephoton > Ephm ) exit
  end do

  cal_nenerge = nenerg_s + ie - 1

  return
end

!***********************************************************************

subroutine Extrap_fpp( Core_resolved, decal_initl, Delta_edge, Energ_s, Eseuil, icheck, EFermi, EFermi_min, fppn, &
                       lmaxabs_t, lseuil, nbseuil, n_Ec, nenerg_s, nenerge, ninit1, nlmamax, nlmsm_f, nspino, &
                       nspinp, numat, rof0, Spinorbite, Taull_tdd, xsect_file )

  use declarations
  implicit none
  include 'mpif.h'

  integer:: icheck, ie, ie_saut, ie_ph, iseuil, iso1, iso2, isp, iss1, iss2, l1, l2, lm1, lm2, lmaxabs_t, &
            lms1, lms2, lmv1, lmv2, lseuil, m1, m2, mv1, mv2, nbseuil, nenerg_s, nenerge, ninit1, n_Ec, nlmamax, &
            nlmsm_f, numat, nspino, nspinp

  character(len=132):: xsect_file

  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nlmamax,nspinp) :: Taull_tdd
  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nspino,nbseuil):: rof0

  logical:: Core_resolved, Saut, Spinorbite

  real(kind=db):: alfa, dde, def, delta, Delta_edge, EFermi, EFermi_min, Ephoton, Eseuil, f_interp1, fp, fpp0, Pente, Rap, &
                  x, x1, x2, y,  y1, y2
  real(kind=db), dimension(n_Ec):: decal_initl
  real(kind=db), dimension(nenerge):: Energe
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(nenerg_s-1:nenerge):: fpp
  real(kind=db), dimension(nenerg_s-1:nenerge,nlmsm_f,nlmsm_f,nbseuil,nbseuil,nspinp):: fppn

  if( icheck > 1 ) write(3,100)

  fppn(:,:,:,:,:,:) = 0._db

  Ephoton = max( Eseuil - 2 / rydb, 0.001_db )
  call fprime(numat,Ephoton,fpp0,fp,xsect_file)  ! fpp0 = avantseuil de f"

  Energe(1:nenerg_s) = Energ_s(1:nenerg_s)

  if( nenerg_s > 1 ) then
    def = Energ_s(nenerg_s) - Energ_s(nenerg_s-1)
  else
    def = 2 / rydb
  endif

  alfa = 1.02_db
  dde = def

! Extrapolation su seuil avec suppression des nouveaux seuils eventuels
  Saut = .false.
  do ie = nenerg_s-1,nenerge
    if( ie > nenerg_s ) then
      dde = alfa * dde
      Energe(ie) = Energe(ie-1) + dde
    endif
    if( Saut ) cycle

    Ephoton = Energe(ie) + Eseuil
! fpp(ie) = f'' atomique pour une energie du spectre etendu
    if( numat > 2 ) then
      call fprime(numat,Ephoton,fpp(ie),fp,xsect_file)
    elseif( ie < nenerg_s+1 ) then
      fpp(ie) = Real( Taull_tdd(ie,3,1,3,1) * rof0(ie,3,1,1,1) * rof0(ie,3,1,1,1), db )
    else
      fpp(ie) = fpp(ie-1)*0.98_db
    endif
    fpp(ie) = fpp(ie) - fpp0
    if( ie <=  nenerg_s ) cycle
    Saut = fpp(ie) > fpp(ie-1)
    if( Saut ) ie_saut = ie
  end do

!Presence d'un seuil: on interpole
  if( Saut ) then
    Pente = ( fpp(ie_saut-2) - fpp(ie_saut-1) ) / ( Energe(ie_saut-2) - Energe(ie_saut-1) )
    do ie = ie_saut,nenerge
      fpp(ie)  = ( Energe(ie) - Energe(ie_saut-1) ) * pente + fpp(ie_saut-1)
      fpp(ie) = max( fpp(ie), 0._db )
    end do
  endif

  EFermi_min = EFermi
  do ie_ph = 1,n_Ec
    delta = decal_initl(ie_ph)
    if( abs(Delta_edge) > eps10 .and. ( ( Core_resolved .and. ie_ph <= ninit1 ) .or. &
        ( .not. Core_resolved .and. ie_ph < n_Ec ) ) ) delta = delta + Delta_edge
    EFermi_min = min( EFermi_min, EFermi + delta )
  end do

  do isp = 1,nspinp

    lm1 = 0
    lms1 = 0
    do l1 = 0,lmaxabs_t
      do m1 = -l1,l1
        lm1 = lm1 + 1
        do iso1 = 1,nspino
          if( Spinorbite ) then
            mv1 = m1 + isp - iso1
            if( abs(mv1) > l1 ) cycle
            iss1 = min( iso1, nspinp )
          else
            mv1 = m1
            iss1 = isp
          endif
          lms1 = lms1 + 1
          lmv1 = l1**2 + l1 + 1 + mv1

          lm2 = 0
          lms2 = 0
          do l2 = 0,lmaxabs_t
            do m2 = -l2,l2
              lm2 = lm2 + 1
              do iso2 = 1,nspino
                if( Spinorbite) then
                  mv2 = m2 + isp - iso2
                  if( abs(mv2) > l2 ) cycle
                  iss2 = min( iso2, nspinp )
                else
                  mv2 = m2
                  iss2 = isp
                endif
                lms2 = lms2 + 1
                lmv2 = l2**2 + l2 + 1 + mv2

                do iseuil = 1, nbseuil

                  if( l1 == lseuil+1 .and. lm1 == lm2 .and. iss1 == iss2 )  then

                    Rap = Real( Taull_tdd(nenerg_s,lm1,iss1,lm2,iss2) * rof0(nenerg_s,lmv1,isp,iso1,iseuil) &
                                                                      * rof0(nenerg_s,lmv2,isp,iso2,iseuil), db ) / fpp(nenerg_s)
                    do ie = nenerg_s+1,nenerge
                      fppn(ie,lms1,lms2,iseuil,iseuil,isp) = Rap * fpp(ie)
                    end do

                  else

                    y1 = Real( Taull_tdd(nenerg_s-1,lm1,iss1,lm2,iss2) * rof0(nenerg_s-1,lmv1,isp,iso1,iseuil) &
                                                                       * rof0(nenerg_s-1,lmv2,isp,iso2,iseuil), db )
                    y2 = Real( Taull_tdd(nenerg_s,lm1,iss1,lm2,iss2) * rof0(nenerg_s,lmv1,isp,iso1,iseuil) &
                                                                     * rof0(nenerg_s,lmv2,isp,iso2,iseuil), db )
                    if( y1 > y2 ) then
! Si fppn decroit on extrapole lineairement
                      x1 = Energe(nenerg_s-1)
                      x2 = Energe(nenerg_s)
                      do ie = nenerg_s,nenerge
                        x = Energe(ie)
                        y = f_interp1(x,x1,x2,y1,y2)
                        fppn(ie,lms1,lms2,iseuil,iseuil,isp) = max( y, 0._db )
                      end do
                    else
! Si fppn croit on prolonge par une constante
                      fppn(nenerg_s:nenerge,lms1,lms2,iseuil,iseuil,isp) = y2
                    endif
                  endif

                end do

              end do
            end do
          end do

        end do
      end do
    end do
  end do

  if( icheck > 2 ) then
    write(3,110) (((( lms1, lms2, iseuil, isp, lms1 = 1,nlmsm_f), lms2 = 1,nlmsm_f), iseuil = 1,nbseuil), isp = 1,nspinp)
    do ie = nenerg_s+1, nenerge
      write(3,120) Energe(ie)*rydb, (((( fppn(ie,lms1,lms2,iseuil,iseuil,isp), lms1 = 1,nlmsm_f), lms2 = 1,nlmsm_f), &
                                                                          iseuil = 1,nbseuil), isp = 1,nspinp)
    end do
  end if

  return
  100 format(/' ---- Extrap_fpp -------',100('-'))
  110 format(/' fpp(lms1,lms2,iseuil,isp)'/ '   Energy ',250(1x,i2,',',i2,',',i1,',',i1,1x) )
  120 format(f9.3,1p,500e11.3)
end

!***********************************************************************

! Calcul de Chi_0

subroutine Chi_0_int(Chi_0,Coef_g,Core_resolved,Decal_initl,Delta_edge,EFermi,Ecent,Elarg,Energ,Energ_s, &
       Eseuil,fppn,Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,icheck,ie,jseuil,lmax, &
       nbseuil,nenerge,nenerg_s,ngamh,ninit1,ninitl,n_Ec,nlm,nlmamax,nlmsm_f,nlms_f,nlms_g,nomfich,ns_dipmag, &
       nseuil,nspino,nspinp,numat,rof_ph,rof0,Spinorbite,Taull_tdd)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: icheck, ie, ie_ph, ief, initl, ipr, is_dipmag, iseuil, iso, iso1, iso2, isp, ispf, ispg, ispm, iss1, &
    iss2, je, jinitl, js, jseuil, l, l1, l2, lm, lm1, lm2, lmax, lms_g, lms1_g, lms2_g, lms1, lms2, lmv1, &
    lmv2, m, m1,  m2, mv1, mv2, nbseuil, nenerge, nenerg_s, ngamh, ninit1, ninitl, n_Ec, nlm, nlmamax, &
    nlmsm_f, nlms_f, nlms_g, ns_dipmag, nseuil, nspinp, nspino, numat

  integer, dimension(nlms_f,2):: i_val, l_val, m_val

  logical:: Core_resolved, Gamma_hole_imp, Gamma_tddft, Spinorbite

  character(len=2):: ch2
  character(len=132):: nomfich

  complex(kind=db):: dch, y1_c
  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nspino,nbseuil):: rof0
  complex(kind=db), dimension(nlm,nspinp,nspino,n_Ec):: rof_ph

  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nlmamax,nspinp) :: Taull_tdd
  complex(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g,2,ns_dipmag):: Chi_0

  real(kind=db):: alfa, Chi_0_r, Chi_0_i, d_i, d_r, dch_i, dch_r, dde, def, delta, EFermi_i, Energ, &
    f_interp1, param, pasmin, t1, t2, Tab_width, x, x1, x2, y1_i, y1_r

  real(kind=db):: Delta_edge, Ecent, EFermi, Elarg, Eseuil, Gamma_max
  real(kind=db), dimension(n_Ec):: Decal_initl
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(nenerge):: E1, E2, Energe, Gamma
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(ninitl,2):: coef_g
  real(kind=db), dimension(nenerg_s-1:nenerge,nlmsm_f,nlmsm_f,nbseuil,nbseuil,nspinp):: fppn

  Chi_0(:,:,:,:,:,:) = (0._db,0._db)

  if( .not. Gamma_tddft )  then
    pasmin = 1.0_db / rydb
! le pas minimum
    do je = 1,nenerg_s-1
      pasmin = min( pasmin, Energ_s(je+1)-Energ_s(je) )
    end do
    param = 0.5_db * pasmin
  end if

  if( ( icheck > 0 .and. ie == 1 ) .or. icheck > 1 ) write(3,100)

  if( nenerg_s > 1 ) then
    def = Energ_s(nenerg_s) - Energ_s(nenerg_s-1)
  else
    def = 2 / rydb
  endif

  alfa = 1.02_db
  dde = def

  Energe(1:nenerg_s) = Energ_s(1:nenerg_s)
  dde = def
  do je = nenerg_s-1,nenerge
    if( je > nenerg_s ) then
      dde = alfa * dde
      Energe(je) = Energe(je-1) + dde
    endif
  end do

  do ispg = 1,2 ! spin de l'etat de coeur

    lms1_g = 0
    lms2_g = 0
    boucle_seuil: do ie_ph = 1,n_Ec

      if( Core_resolved ) then
        initl = ie_ph
        if( abs( coef_g(initl,ispg) ) < eps10 ) cycle
        lms1_g = lms1_g + 1
        lms2_g = lms2_g + 1
        if( initl > ninit1 ) then
          iseuil = 2
        else
          iseuil = 1
        endif
      else
        iseuil = ie_ph
        if( iseuil == 1 ) then
          lms1_g = 1
          lms2_g = 0
          do jinitl = 1,ninit1
            if( abs( coef_g(jinitl,ispg) ) > eps10 ) lms2_g = lms2_g + 1
          end do
        else
          lms1_g = 1
          do jinitl = 1,ninit1
            if( abs( coef_g(jinitl,ispg) ) > eps10 ) lms1_g = lms1_g + 1
          end do
          lms2_g = 0
          do jinitl = 1,ninitl
            if( abs( coef_g(jinitl,ispg) ) > eps10 ) lms2_g = lms2_g + 1
          end do
       endif
      endif

! Energe est la gamme interne shiftee
      delta = decal_initl(ie_ph)
      if( abs(Delta_edge) > eps10 .and. ( ( Core_resolved .and. ie_ph <= ninit1 ) .or. &
          ( .not. Core_resolved .and. ie_ph < n_Ec ) ) ) delta = delta + Delta_edge

      EFermi_i = EFermi + delta

      Energe(:) = Energe(:) + delta

! Elargissement
      Gamma(:) = 0._db
      if( Gamma_tddft )  then

        if( Gamma_max > eps10 ) call gammarc(Ecent,Elarg,Gamma_max,EFermi_i,nenerge,Energe,Gamma)

        if( Gamma_hole_imp ) then
          if( ngamh == 1 ) then
            Gamma(:) = Gamma(:) + Gamma_hole(1)
          elseif( ngamh == n_Ec ) then
            Gamma(:) = Gamma(:) + Gamma_hole(ie_ph)
          elseif( ie_ph <= ninit1 ) then
            Gamma(:) = Gamma(:) + Gamma_hole(1)
          else
            Gamma(:) = Gamma(:) + Gamma_hole(2)
          endif
        else
          js = jseuil + iseuil - 1
          Gamma_hole(1) = Tab_width(Eseuil,js,nseuil,numat)
          Gamma(:) = Gamma(:) + Gamma_hole(1)
        end if
        if( ( icheck > 0 .and. ie == 1 ) .or. icheck > 1 ) write(3,105) ie_ph, iseuil, Gamma_hole(1) * rydb

! On prend la mi-largeur
        Gamma(:) = Gamma(:) / 2

      elseif( ie_ph == 1 ) then

        if( ( icheck > 0 .and. ie == 1 ) .or. icheck > 1 ) write(3,'(/A)' ) ' Gamma = 0'

      end if

! Indice du niveau de Fermi
      if( EFermi_i < Energe(1) .or. EFermi_i > Energe(nenerg_s) ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,500) EFermi*rydb, Energ_s(1)*rydb, Energ_s(nenerg_s)*rydb
        end do
        stop
      end if
      ief = 0
      do je = 1,nenerg_s-1
        if( Energe(je) > EFermi_i  ) exit
      end do
      ief = je

! Creation des tableaux bornes des intervales:
      do je = 2,nenerge-1
        e1(je) =  0.5_db * ( Energe(je-1) + Energe(je) )
        e2(je) =  0.5_db * ( Energe(je) + Energe(je+1) )
      end do

! Bornes de l'intervale:
      e1(1) = Energe(1) - 0.5_db * ( Energe(2) - Energe(1) )
      e2(1) = 0.5_db * ( Energe(1) + Energe(2) )
      e1(nenerge) =  0.5_db * ( Energe(nenerge-1) + Energe(nenerge) )
      e2(nenerge) = Energe(nenerge) + 0.5_db * ( Energe(nenerge-1) - Energe(nenerge) )

      e1(ief) = Efermi_i

      do is_dipmag = 1,ns_dipmag

! Spin etat final
        if( is_dipmag == 1 ) then
          ispf = ispg
        else
          ispf = 3 - ispg
        endif
        ispm = min( ispf, nspinp )

        lm1 = 0
        lms1 = 0
        do l1 = 0,lmax
          do m1 = -l1,l1
            lm1 = lm1 + 1
            do iso1 = 1,nspino
              if( Spinorbite ) then
                mv1 = m1 + ispf - iso1
                if( abs(mv1) > l1 ) cycle
                iss1 = min( iso1, nspinp)
              else
                mv1 = m1
                iss1 = min( ispf, nspinp )
              endif
              lmv1 = l1**2 + l1 + 1 + mv1
              lms1 = lms1 + 1
              l_val(lms1,ispf) = l1
              m_val(lms1,ispf) = m1
              i_val(lms1,ispf) = iso1

              lm2 = 0
              lms2 = 0
              do l2 = 0,lmax
                do m2 = -l2,l2
                  lm2 = lm2 + 1
                  do iso2 = 1,nspino
                    if( Spinorbite ) then
                      mv2 = m2 + ispf - iso2
                      if( abs(mv2) > l2 ) cycle
                      iss2 = min( iso2, nspinp )
                    else
                      mv2 = m2
                      iss2 = min( ispf, nspinp )
                    endif
                    lmv2 = l2**2 + l2 + 1 + mv2
                    lms2 = lms2 + 1

                    Chi_0_r = 0._db ; Chi_0_i = 0._db

                    do je = ief,nenerge

                      t1 = Energ - e1(je)
                      t2 = Energ - e2(je)

                      if( je > nenerg_s ) then
                        dch = cmplx( fppn(je,lms1,lms2,iseuil,iseuil,ispm), 0._db )
                      else
                        dch = Taull_tdd(je,lm1,iss1,lm2,iss2) * rof0(je,lmv1,ispm,iso1,iseuil) * rof0(je,lmv2,ispm,iso2,iseuil)
                      endif
                      dch = dch / ( rof_ph(lmv1,ispm,iso1,ie_ph) * rof_ph(lmv2,ispm,iso2,ie_ph) )
                      dch_r = real( dch, db )
                      dch_i = aimag( dch )

                      if( Gamma_tddft ) then

                        t1 = t1 / Gamma(je)
                        t2 = t2 / Gamma(je)
                        d_i = - ( atan(t2) - atan(t1) ) / pi
                        d_r = - 0.5_db * log( (1+t2**2)/(1+t1**2) ) / pi
                        Chi_0_r = Chi_0_r + dch_r * d_r - dch_i * d_i
                        Chi_0_i = Chi_0_i + dch_r * d_i + dch_i * d_r

                      else    ! gamma = 0

! Eviter les divergences qd les gammes se chevauchent
                        if( abs(t1) < param ) t1 =  param
                        if( abs(t2) < param ) t2 = -param
! On elimine la divergence en 1/e du f' (t1 -> -t1;
! si jamais |t1| = |t2| les contributions des deux cotes du point s'annulent.
! abs change le signe quand t1 et t2 sont de signes opposes.
                        if( abs(t1+t2) > eps10 ) then
                          Chi_0_r = Chi_0_r - dch_r * log( abs(t2/t1) ) / pi
                          Chi_0_i = Chi_0_i - dch_i * log( abs(t2/t1) ) / pi
                        endif

                        if( Energe(je) < Energ - eps10 .or. Energ < Efermi_i ) cycle
! Partie imaginaire, on interpole quand on depasse le point courant
                        if( abs( Energe(je) - Energ ) < eps10 .or. je == 1 ) then
                          Chi_0_i = Chi_0_i + dch_r
                        elseif( Energe(je-1) < Energ ) then
                          if( je-1 > nenerg_s ) then
                            y1_r = fppn(je-1,lms1,lms2,iseuil,iseuil,ispm)
                          else
                            y1_c = Taull_tdd(je-1,lm1,iss1,lm2,iss2) * rof0(je-1,lmv1,ispm,iso1,iseuil) &
                                                                     * rof0(je-1,lmv2,ispm,iso2,iseuil)
                            y1_r = real( y1_c, db )
                            y1_i = aimag( y1_c )
                          endif
                          x  = Energ
                          x1 = Energe(je-1)
                          x2 = Energe(je)
                          Chi_0_i = Chi_0_i + f_interp1(x,x1,x2,y1_r,dch_r)
                          if( je-1 <= nenerg_s ) Chi_0_r = Chi_0_r + f_interp1(x,x1,x2,y1_i,dch_i)
                        endif

                      end if

                    end do

                    do lms_g = lms1_g,lms2_g
                      Chi_0(lms1,lms2,lms_g,lms_g,ispf,is_dipmag) = cmplx( Chi_0_r, Chi_0_i, db )
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do

      end do  ! fin boucle is_dipmag

      Energe(:) = Energe(:) - delta

    end do boucle_seuil

  end do  ! fin boucle ispg

! Ecriture des termes diagonaux dans le fichier bavard:
  if( icheck > 1 ) then
    write(3,130) (((( lms1, lms1, lms_g, lms_g, ispf, is_dipmag, &
                                         lms1 = 1,nlms_f), lms_g = 1,nlms_g), ispf = 1,2 ), is_dipmag = 1,ns_dipmag)
    write(3,120) Energ*rydb, (((( Chi_0(lms1,lms1,lms_g,lms_g,ispf,is_dipmag), &
                                         lms1 = 1,nlms_f), lms_g = 1,nlms_g), ispf = 1,2 ), is_dipmag = 1,ns_dipmag)
  end if

  if( icheck > 2 ) then
    write(3,140) ((( l, ch2(m), isp, isp = 1,nspinp),  m = -l,l), l = 0,lmax )
    do je = 1,nenerg_s
      write(3,120) Energ_s(je)*rydb, (( Taull_tdd(je,lm,isp,lm,isp), isp = 1,nspinp), lm = 1,(1+lmax)**2 )
    end do

    write(3,150) (((( l, ch2(m), isp, iseuil, iseuil = 1,nbseuil), isp = 1,nspinp), m = -l,l), l = 0,lmax )
    do je = 1,nenerg_s
      write(3,120) Energ_s(je)*rydb, (( real( rof0(je,lm,isp,min(isp,nspino),1:nbseuil)), isp = 1,nspinp), lm = 1,(1+lmax)**2 )
    end do

    write(3,160) Energ*rydb, ((( isp ,iso, ie_ph, isp = 1,nspinp), iso = 1,nspino), ie_ph = 1,n_Ec)
    lm = 0
    do l = 0,lmax
      do m = -l, l
        lm = lm + 1
        write(3,170) l, m, (( real( rof_ph(lm,:,iso,ie_ph) ), iso = 1,nspino), ie_ph = 1,n_Ec)
      end do
    end do
  end if

  if( icheck > 1 ) call write_Chi_0(Chi_0,Energ,i_val,ie,l_val,m_val,nlms_f,nlms_g,nomfich,ns_dipmag)

  return
  100 format(/' ---- Chi_0_int -------',100('-'))
  102 format('   Chi_0 calculation, state =',i2,' on',i2)
  103 format('   Chi_0 calculation, state =',i2,' on',i2,', computer,', i3)
  105 format(/' ie_ph =',i2,', iseuil =',i2,', Gamma_hole =',f7.3,' eV')
  120 format(f9.3,1p,500e11.3)
  130 format(/' Chi_0(lms,lms,lms_g,lms_g,ispf,is_dipmag) for diagonal terms, ie_ph ',i2,/ &
              '   Energy ',250(2x,'(',i2,',',i1,',',i1,',',i1,')',3x,'Im',5x) )
  140 format(/' Taull_tdd(l,m,isp,l,m,isp) for diagonal terms',/ '   Energy ',250(8x,'(',i1,',',a2,',',i1,') ',6x) )
  150 format(/' Monopole rof0(l,m,isp,iso,iseuil) for iso = isp',/ '   Energy ',250('(',i1,',',a2,',',i1,',',i1,') ') )
  160 format(/' rof0 at the photon energy versus (isp,iso,ie_ph) at ',f10.3,' eV', / &
              '   l   m', 100(4x,3i2,3x))
  170 format(2i3,1p,100e13.3)
  500 format(//' Error: the Fermi energy is out of the energy calculation range !',/ 5x,'EFermi =',f10.3,' eV',/ &
            5x,'E_min  =',f10.3,' eV,  E_max  =',f10.3,' eV',// ' It is not possible in the TDDFT part !',/)
end

!**********************************************************************

subroutine write_Chi_0(Chi_0,Energ,i_val,je,l_val,m_val,nlms_f,nlms_g,nomfich,ns_dipmag)

  use declarations
  implicit none

  character(len=132):: chi_conv, nomfich
  character(len=1), dimension(nlms_f,2):: sgn

  integer:: is_dipmag, isp, je, l, lms, lms1, lms2, nlms_f, nlms_g, ns_dipmag
  integer, dimension(nlms_f,2):: i_val, l_val, m_val

  complex(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g,2,ns_dipmag):: Chi_0

  logical:: Diag

  real(kind=db) :: Energ

  Diag = .true.

! Nom du fichier de sortie
  chi_conv = nomfich
  l = len_trim( chi_conv )
  chi_conv(l+1:l+10) = '_Chi_0.txt'

  if( je == 1 ) then
    do isp = 1,2
      do lms1 = 1,nlms_f
        if( m_val(lms1,isp) < 0 ) then
          sgn(lms1,isp) = '-'
        else
          sgn(lms1,isp) = '+'
        endif
      end do
    end do
    open(31, file = chi_conv)
    if( Diag ) then
      isp = 1
      write(31,105) (( l_val(lms1,isp), sgn(lms1,isp), abs(m_val(lms1,isp)), i_val(lms1,isp), &
                       l_val(lms,isp), sgn(lms,isp), abs(m_val(lms,isp)), i_val(lms,isp), lms1 = 1,nlms_f), lms = 1,nlms_g)
    else
      write(31,110) ((((( l_val(lms1,isp), sgn(lms1,isp), abs(m_val(lms1,isp)), i_val(lms1,isp), l_val(lms2,isp), sgn(lms2,isp), &
                        abs(m_val(lms2,isp)), i_val(lms2,isp), lms, isp, is_dipmag, lms1 = 1,nlms_f), lms2 = 1,nlms_f), &
                        lms = 1,nlms_g), isp = 1,2), is_dipmag = 1,ns_dipmag)
    endif
  else
    open(31, file = chi_conv, position='append')
  endif

  if( Diag ) then
    isp = 1
    is_dipmag = 1
    write(31,120) Energ*rydb, (( Chi_0(lms1,lms1,lms,lms,isp,is_dipmag), lms1 = 1,nlms_f), lms = 1,nlms_g)
  else
    write(31,120) Energ*rydb, ((((( Chi_0(lms1,lms2,lms,lms,isp,is_dipmag), lms1 = 1,nlms_f), lms2 = 1,nlms_f), lms = 1,nlms_g), &
                                                                          isp = 1,2), is_dipmag = 1,ns_dipmag)
  endif
  Close(31)

  return
  105 format('   Energy ',5000(' (',i1,',',a1,i1,',',i1,',',i1,',',a1,i1,',',i1,')  Imag '))
  110 format('   Energy ',5000(' (',2(i1,',',a1,i1,',',i1,','),i1,',',i1,',',i1,')'))
  120 format(f10.3,1p,11664e11.3)
end

!***********************************************************************

character(len=2) function ch2(m)

  integer:: m

  if( m < 0 ) then
    ch2(1:1) = '-'
  elseif( m == 0 ) then
    ch2(1:1) = '.'
  else
    ch2(1:1) = '+'
  endif
  ch2(2:2) = achar( abs(m) + 48 )

  return
end

!***********************************************************************

! Sous-programme qui calcule le noyau pour le calcul TDDFT

! Le developement en harmoniques spheriques du potentiel coulombien est pris d'apres J.Phys.A:Math. Gen. 39(2006) 8613-8630
! On le calcule dans le cas harmoniques complexe.

subroutine kernel(Atomic_scr,coef_g,Core_resolved,Dipmag,Dyn_eg,Dyn_g,Energ,fxc,icheck,ie,Kern,Kern_fac, &
                lmax,lseuil,m_g,nbseuil,ninit1,ninitl,n_Ec,nlm,nlm_fp,nlms_f,nlms_g,nomfich,nr,nrm, &
                ns_dipmag,nspino,nspinp,psii,r,Rmtsd,RPALF,Spinorbite,zet)

  use declarations
  implicit none

  integer:: icheck, ie, lmax, lseuil, nbseuil, ninit1, ninitl, n_Ec, nlm, nlm_fp, nlms_f, nlms_g, &
     nr, nrm, ns_dipmag, nspinp, nspino
   integer:: is_dipmag, is1, is2, isf1, isf2, isg1, isg2, iso1, iso2, isp1, isp2, ist1, ist2, iz1, iz2, l0, l1, &
    l2, lcut, lg, lm1, lm2, lmg1, lmg2, lmp1, lmp2, lms1, lms2, lmv1, lmv2, lp1, lp2, m0, m1, m2, mg1, mg2, mp1, mp2, mv1, &
    mv2
  integer, dimension(ninitl,2):: m_g
  integer, dimension(nlms_f,2):: i_val, l_val, m_val

  character(len=132):: nomfich

  complex(kind=db):: Gaunt_xcc

  logical:: Atomic_scr, Core_resolved, Dipmag, Dyn_eg, Dyn_g, RPALF, Spinorbite, Ylm_comp

  real(kind=db), intent(in):: Energ, Kern_fac, Rmtsd
  real(kind=db), dimension(nr), intent(in):: r
  real(kind=db), dimension(nrm,nbseuil), intent(in):: psii
  real(kind=db), dimension(nr,2,2), intent(in):: fxc
  real(kind=db), dimension(ninitl,2), intent(in):: coef_g
  real(kind=db), dimension(nr,nlm,nlm_fp,nspinp,nspino,n_Ec), intent(in):: zet

  real(kind=db):: Angl, Angl1, Angl2, f_integr3, fac, Gaunt4Y, intrad_r
  real(kind=db), dimension(nr):: f, r1, r2, t1, t2, zet_1, zet_2
  real(kind=db),dimension(nlms_f,nlms_f,nlms_g,nlms_g,2,2,ns_dipmag):: Kern

  if( ( icheck > 1 .and. ie == 1 ) .or. icheck > 2 ) then
    write(3,100)
    write(3,110)
  endif

  Kern(:,:,:,:,:,:,:) = 0._db
  Ylm_comp = .true.

! l'etat initial est calcule en corestate:
  lg = lseuil

! L'indice 0 porte sur le developpement du potentiel coulombien, et
! l'indice g sur l'etat de coeur
! psii est reel

! le cut-off: voir les regles de selection pour les coef de gaunt

  lcut = lg + lmax

! Boucle sur les etats initiaux

  do isg1 = 1,2

    lmg1 = 0
    do ist1 = 1,ninitl

      if( abs( coef_g(ist1,isg1) ) < eps10 ) cycle
      lmg1 = lmg1 + 1
      mg1 = m_g(ist1,isg1)

      if( ist1 > ninit1 ) then
        is1 = 2
      else
        is1 = 1
      end if
      if( Core_resolved ) then
        iz1 = ist1
      else
        iz1 = is1
      endif

      do isg2 = 1,2

        lmg2 = 0
        do ist2 = 1,ninitl
          if( abs( coef_g(ist2,isg2) ) < eps10 ) cycle
          lmg2 = lmg2 + 1
          mg2 = m_g(ist2,isg2)

          if( ist2 > ninit1 ) then
            is2 = 2
          else
            is2 = 1
          endif
          if( Core_resolved ) then
            iz2 = ist2
          else
            iz2 = is2
          endif

          do isf1 = 1,2         ! Spin du premier etat final

            if( isf1 == isg1 ) then
              is_dipmag = 1
            else
              if( .not. dipmag ) cycle
              is_dipmag = 2
            endif
! Pour traiter le cas des calculs non magnetiques
            isp1 = min(isf1,nspinp)

            do isf2 = 1,2       ! Spin du deuxieme etat final

              if( isf2 /= isg2 .and. .not. dipmag ) cycle
              if( (isf1 == isf2 .and. isg1 /= isg2) .or. (isf1 /= isf2 .and. isg1 == isg2)) cycle

              isp2 = min(isf2,nspinp)

              lm1 = 0
              lms1 = 0
              do l1 = 0,lmax
                do m1 = -l1,l1
                  lm1 = lm1 + 1

! boucle cas non spherique
                  lmp1 = 0
                  do lp1 = 0,lmax

                    do mp1 = -lp1,lp1
                      if( nlm_fp == 1 .and. ( lp1 /= l1 .or. mp1 /= m1 ) ) cycle
                      lmp1 = lmp1 + 1

                      do iso1 = 1,nspino

                        if( Spinorbite .and. nlm_fp == 1 ) then
                          mv1 = m1 + isf1 - iso1
                          if( abs(mv1) > l1 ) cycle
                        else
                          mv1 = m1
                        endif
                        lmv1 = l1**2 + l1 + 1 + mv1
                        lms1 = lms1 + 1
                        l_val(lms1,isf1) = l1
                        m_val(lms1,isf1) = m1
                        i_val(lms1,isf1) = iso1

                        lm2 = 0
                        lms2 = 0
                        do l2 = 0,lmax

                          do m2 = -l2,l2
                            lm2 = lm2 + 1

                            lmp2 = 0
                            do lp2 = 0,lmax

                              do mp2 = -lp2,lp2
                                if( nlm_fp == 1 .and. ( lp2 /= l2 .or. mp2 /= m2 ) ) cycle
                                lmp2 = lmp2 + 1

                                do iso2 = 1,nspino

                                  if( Spinorbite .and. nlm_fp == 1 ) then
                                    mv2 = m2 + isf2 - iso2
                                    if( abs(mv2) > l2 ) cycle
                                  else
                                    mv2 = m2
                                  endif
                                  lmv2 = l2**2 + l2 + 1 + mv2
                                  lms2 = lms2 + 1
! Atomic screening
                                  if( Atomic_scr .and. l1 /= lg + 1 .and. l2 /= lg + 1 ) cycle

                                  zet_1(1:nr) = psii(1:nr,is1) * zet(1:nr,lmv1,lmp1,isp1,iso1,iz1)
                                  zet_2(1:nr) = psii(1:nr,is2) * zet(1:nr,lmv2,lmp2,isp2,iso2,iz2)
! Noyau coulombien
                                  do l0 = 0,lcut
                                    if( l0 < abs(l1 - lg) .or. l0 > l1 + lg ) cycle
                                    if( l0 < abs(l2 - lg) .or. l0 > l2 + lg ) cycle
                                    if( mod( l1 + l0 + lg, 2 ) /= 0 ) cycle
                                    if( mod( l2 + l0 + lg, 2 ) /= 0 ) cycle

                                    r1(1:nr) = r(1:nr)**(l0+1)
                                    r2(1:nr) = r(1:nr)**(-l0)

                                    f(1:nr) = r1(1:nr) * zet_2(1:nr)

                                    call ffintegr2_r(t1,f,r,nr,1,Rmtsd)

                                    f(1:nr) = r2(1:nr) * zet_2(1:nr)

                                    call ffintegr2_r(t2,f,r,nr,-1,Rmtsd)

                                    f(1:nr) = zet_1(1:nr) * ( r2(1:nr) * t1(1:nr) + r1(1:nr) * t2(1:nr) )

                                    intrad_r = f_integr3(r,f,1,nr,Rmtsd)
! le 2 de 2 / (r-r')
                                    fac = 8 * pi * intrad_r / ( 2 * l0 + 1 )

                                    do m0 = -l0,l0

                                      Angl1 = real( Gaunt_xcc(l1,mv1,l0,m0,lg,mg1,Ylm_comp), db)
                                      if( abs( Angl1 ) < eps10 ) cycle
                                      Angl2 = real( conjg( Gaunt_xcc(l2,mv2,l0,m0,lg,mg2,Ylm_comp) ), db )
                                      if( abs( Angl2 ) < eps10 ) cycle

                                      Kern(lms1,lms2,lmg1,lmg2,isf1,isf2,is_dipmag) &
                                                                = Kern(lms1,lms2,lmg1,lmg2,isf1,isf2,is_dipmag) &
                                                                + coef_g(ist1,isg1) * coef_g(ist2,isg2) * Angl1 * Angl2 * fac

                                    end do  ! fin boucle m0
                                  end do    ! fin boucle l0

! On rajoute la partie xc en LDA:
                                  if( .not. ( RPALF .or. (Dyn_g .and. ist1 /= ist2) .or. (Dyn_eg .and. is1 /= is2) ) ) then

! Les harmoniques sont supposees complexes
! Dans Gaunt4Y, ce sont les 2 premieres harmoniques qui sont complexe-conjuguees.
                                    Angl = Gaunt4Y(lg,mg2,l1,mv1,l2,mv2,lg,mg1)

                                    if( Abs( Angl ) > eps10 ) then

                                      if( isf1 == isg1 ) then
                                        f(1:nr) = fxc(1:nr,isf1,isf2) * zet_1(1:nr) * zet_2(1:nr)
                                      else
                                        f(1:nr) = fxc(1:nr,isf1,isg1) * zet_1(1:nr) * zet_2(1:nr)
                                      endif
                                      fac = f_integr3(r,f,1,nr,Rmtsd)

                                      Kern(lms1,lms2,lmg1,lmg2,isf1,isf2,is_dipmag) &
                                                                      = Kern(lms1,lms2,lmg1,lmg2,isf1,isf2,is_dipmag) &
                                                                      + coef_g(ist1,isg1) * coef_g(ist2,isg2) * Angl * fac

                                    endif
                                  endif

                                  if( ( ( icheck > 1 .and. ie == 1 ) .or. icheck > 2 ) &
                                        .and. abs(Kern(lms1,lms2,lmg1,lmg2,isf1,isf2,is_dipmag) ) > eps10 ) &
                                      write(3,120) isg1, ist1, isg2, ist2, isf1, isf2, l1, m1, iso1, l2, m2, iso2, is_dipmag, &
                                                   Kern(lms1,lms2,lmg1,lmg2,isf1,isf2,is_dipmag)

                                end do ! fin boucle iso2
                              end do   ! fin boucle mp2
                            end do     ! fin boucle lp2
                          end do       ! fin boucle m2
                        end do         ! fin boucle l2

                      end do ! fin boucle iso1
                    end do   ! fin boucle mp1
                  end do     ! fin boucle lp1
                end do       ! fin boucle m1
              end do         ! fin boucle l1

            end do   ! fin boucle isf2
          end do     ! fin boucle isf1

        end do   ! fin boucle ist2
      end do     ! fin boucle isg2

    end do   ! fin boucle ist1
  end do     ! fin boucle isg1

  if( abs( Kern_fac - 1._db ) > eps10 ) Kern(:,:,:,:,:,:,:) = Kern_fac * Kern(:,:,:,:,:,:,:)

  if( icheck > 1 ) call write_Kern(Kern,Energ,i_val,ie,l_val,m_val,nlms_f,nlms_g,nomfich,ns_dipmag)

  return
  100 format(/' ---- Kernel -------',100('-'))
  110 format(/' isg1 ist1 isg2 ist2 isf1 isf2   l1   m1 iso1   l2   m2 iso2 is_dipmag    Kern ')
  120 format(13i5,1x,1p,100e11.3)
end

!**********************************************************************

subroutine write_Kern(Kern,Energ,i_val,je,l_val,m_val,nlms_f,nlms_g,nomfich,ns_dipmag)

  use declarations
  implicit none

  character(len=132):: chi_conv, nomfich
  character(len=1), dimension(nlms_f,2):: sgn

  integer:: is_dipmag, isp, isp1, isp2, je, l, lms, lms1, lms2, nlms_f, nlms_g, ns_dipmag
  integer, dimension(nlms_f,2):: i_val, l_val, m_val

  real(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g,2,2,ns_dipmag):: Kern

  logical:: Diag

  real(kind=db) :: Energ

  Diag = .true.

! Outpur file name
  chi_conv = nomfich
  l = len_trim( chi_conv )
  chi_conv(l+1:l+11) = '_Kern.txt'

  if( je == 1 ) then
    do isp = 1,2
      do lms1 = 1,nlms_f
        if( m_val(lms1,isp) < 0 ) then
          sgn(lms1,isp) = '-'
        else
          sgn(lms1,isp) = '+'
        endif
      end do
    end do
    open(31, file = chi_conv)
    if( Diag ) then
      write(31,105) (((( l_val(lms1,isp), sgn(lms1,isp), abs(m_val(lms1,isp)), i_val(lms1,isp), &
                        lms, isp, is_dipmag, lms1 = 1,nlms_f), lms = 1,nlms_g), isp = 1,2), is_dipmag = 1,ns_dipmag)
    else
      write(31,110) ((((( l_val(lms1,isp), sgn(lms1,isp), abs(m_val(lms1,isp)), i_val(lms1,isp), l_val(lms2,isp), sgn(lms2,isp), &
                        abs(m_val(lms2,isp)), i_val(lms2,isp), lms, isp, is_dipmag, lms1 = 1,nlms_f), lms2 = 1,nlms_f), &
                        lms = 1,nlms_g), isp = 1,2), is_dipmag = 1,ns_dipmag)
    endif
  else
    open(31, file = chi_conv, position='append')
  endif

  if( Diag ) then
    write(31,120) Energ*rydb, (((( Kern(lms1,lms1,lms,lms,isp1,isp1,is_dipmag), lms1 = 1,nlms_f), lms = 1,nlms_g), &
                                                                          isp1 = 1,2), is_dipmag = 1,ns_dipmag)
  else
    write(31,120) Energ*rydb, (((((( Kern(lms1,lms2,lms,lms,isp1,isp1,is_dipmag), lms1 = 1,nlms_f), lms2 = 1,nlms_f), &
                                                       lms = 1,nlms_g), isp1 = 1,2), isp2 = 1,2), is_dipmag = 1,ns_dipmag)
  endif
  Close(31)

  return
  105 format('   Energy ',5000(' (',i1,',',a1,i1,',',i1,',',i1,',',i1,',',i1,')'))
  110 format('   Energy ',5000(' (',2(i1,',',a1,i1,',',i1,','),i1,',',i1,',',i1,')'))
  120 format(f10.3,1p,11664e11.3)
end

!***********************************************************************

! Transformation, si harmoniques reelles
! On ne peut pas etre en spinorbite

subroutine Tr_Taull( icheck, nlmamax, nenerg_s, nspinp, Taull_tdd )

  use declarations
  implicit none

  integer:: i, icheck, ie, isp, nenerg_s, nlmamax, nspinp

  complex(kind=db), dimension(nlmamax,nlmamax):: A, Trans
  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nlmamax,nspinp) :: Taull_tdd

  Call Cal_Trans(nlmamax,Trans)

  do ie = 1,nenerg_s

    do isp = 1,nspinp

      do i = 1,nlmamax
        A(:,i) = Taull_tdd(ie,:,isp,i,isp)
      end do

      A = matmul( A, Trans )
      A = matmul( conjg( transpose( Trans ) ), A )

      if( icheck > 2 ) then
        write(3,110) ie, isp
        do i = 1,nlmamax
          write(3,120) A(i,:)
        end do
      endif

      do i = 1,nlmamax
        Taull_tdd(ie,:,isp,i,isp) = A(:,i)
      end do

    end do

  end do

  return
  110 format(/' Taull in complex harmonic basis, ie =',i4,', isp =',i2)
  120 format(1p,100e13.5)
end

!***********************************************************************

! On calcule Chi suivant la formule : Chi = Chi0 + Chi0 * K * Chi
!                                     Chi = (1 - Chi0 * K )^(-1) * Chi0

! Chi_0 en entree est en convention cristallo
! Chi en sortie est en convention physique (compl. conj. de cristallo)

subroutine Cal_Chi(Chi, Chi_0, coef_g, Energ, First_E, icheck, iopsymc_25, Kern, lmax_probe, lmax, &
              lseuil, nd3, ninitl, nlm_probe, nlm_fp, nlms_f, nlms_g, &
              nomfich, NRIXS, ns_dipmag, nspino, Optic, Quadrupole, Spinorbite)

  use declarations
  implicit none

  character(len=132), intent(in):: nomfich

  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,nd3,nd3,2,2,ns_dipmag), intent(out):: Chi
  complex(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g,2,ns_dipmag):: Chi_0
  complex(kind=db), dimension(:), allocatable:: V
  complex(kind=db), dimension(:,:), allocatable:: A, B

  integer:: i, i1, i2, icheck, is_dipmag, iopsymc_25, imparite, isp, isp_f1, isp_f2, isp_g1, isp_g2, iso, &
    j, l, l_f1, l_f2, l_g1, l_g2, l1, l2, lm, lmax, lmax_probe, lms, lms_f1, lms_f2, lms_g1, lms_g2, &
    lms1, lms2, lmsg1,lmsg2, lseuil, m, m_f1, m_f2, nd3, ndim, ninitl, &
    nlm_fp, nlm_probe, nlms, nlms_f, nlms_g, nlms_gg, ns_dipmag, nspino

  integer, dimension(nlms_f,2):: i_val, l_val, m_val

  logical:: First_E, NRIXS, Optic, Quadrupole, Spinorbite, Stop_job

  real(kind=db), intent(in):: Energ
  real(kind=db), dimension(ninitl,2),intent(in):: coef_g
  real(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g,2,2,ns_dipmag), intent(in):: Kern
  real(kind=db), dimension(:,:), allocatable:: K

  if( icheck > 1 ) write(3,100)
  if( icheck > 1 ) write(3,110) Energ * rydb

  Chi(:,:,:,:,:,:,:) = ( 0._db, 0._db )
  Stop_job = .false.

  if( iopsymc_25 == 0 .or. Quadrupole .or. ns_dipmag == 2 .or. Optic .or. NRIXS ) then
    imparite = 2
  elseif( mod(lseuil,2) == 0 ) then
    imparite = 1
  else
    imparite = 0
  endif

  if( imparite /= 2 ) then
    ndim = 0
    do l = 0, lmax
      if( mod(l,2) /= imparite .and. imparite /= 2 ) cycle
      if( Spinorbite ) then
        ndim = ndim + 4 * l + 1
      else
        ndim = ndim + 2 * l + 1
      endif
    end do
  else
    ndim = nlms_f
  endif

  ndim = 2 * ndim * nlms_g

  do isp = 1,2
    lms = 0
    do l = 0,lmax
      do m = -l,l
        do iso = 1,nspino
          if( Spinorbite .and. nlm_fp == 1 ) then
            if( abs(m + isp - iso) > l ) cycle
          endif
          lms = lms + 1
          l_val(lms,isp) = l
          m_val(lms,isp) = m
          i_val(lms,isp) = iso
        end do
      end do
    end do
  end do

  nlms = lms

  do is_dipmag = 1,ns_dipmag

    do ! retour en cas de pb (pour faire ecriture bavarde)

    allocate( A(ndim,ndim), B(ndim,ndim), K(ndim,ndim) )

    A(:,:) = (0._db, 0._db)
    B(:,:) = (0._db, 0._db)
    K(:,:) = 0._db

    i1 = 0
    do isp_g1 = 1,2
      do lms_g1 = 1,nlms_g

        do isp_f1 = 1,2         ! Spin du premier etat final

          if( (isp_f1 /= isp_g1 .and. is_dipmag == 1) .or. (isp_f1 == isp_g1 .and. is_dipmag == 2) ) cycle

          do lms_f1 = 1,nlms
            l1 = l_val(lms_f1,isp_f1)
            if( mod(l1,2) /= imparite .and. imparite /= 2 ) cycle
            i1 = i1 + 1

            i2 = 0
            do isp_g2 = 1,2
              do lms_g2 = 1,nlms_g

                do isp_f2 = 1,2       ! Spin du deuxieme etat final

                  if( (isp_f2 /= isp_g2 .and. is_dipmag == 1) .or. (isp_f2 == isp_g2 .and. is_dipmag == 2) ) cycle
                  if( (isp_f1 == isp_f2 .and. isp_g1 /= isp_g2) .or. (isp_f1 /= isp_f2 .and. isp_g1 == isp_g2) ) cycle

                  do lms_f2 = 1,nlms
                    l2 = l_val(lms_f2,isp_f2)
                    if( mod(l2,2) /= imparite .and. imparite /= 2 ) cycle
                    i2 = i2 + 1

! Conjugue car convention physique
                    if( isp_g1 == isp_g2 .and. isp_f1 == isp_f2 ) then
                      A(i1,i2) = Conjg( Chi_0(lms_f1,lms_f2,lms_g1,lms_g2,isp_f1,is_dipmag) )
                    endif

                    K(i1,i2) = Kern(lms_f1,lms_f2,lms_g1,lms_g2,isp_f1,isp_f2,is_dipmag)

                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    if( icheck > 2 ) then
      write(3,'(/A)') '  Chi_0'
      do lm = 1,ndim
        write(3,140) A(lm,:)
      end do
      write(3,'(/A)') '  Kern'
      do lm = 1,ndim
        write(3,140) K(lm,:)
      end do
    endif

    do i = 1,ndim
      do j = 1,ndim
        B(i,j) = - sum( A(i,:) * K(:,j) )
      end do
    end do

    deallocate( K )

    do i = 1,ndim
      B(i,i) = 1 + B(i,i) ! B = (1 - chi0*K)
    end do

    if( icheck > 2 ) then
      write(3,'(/A)') '  (1 - Chi_0.K)(lm,iso,isg_ist)'
      do lm = 1,ndim
        write(3,140) B(lm,:)
      end do
    endif

    if( Stop_job ) stop

! B = (1 - chi0*K)**(-1)
    call invcomp(ndim,B,ndim,ndim,0,Stop_job)

    if( Stop_job ) then
      deallocate( A, B )
      icheck = 3
    else
      exit
    endif

    end do  ! boucle pour retour en case de pb

    if( icheck > 2 ) then
      write(3,'(/A)') '  (1 - Chi_0.K)**(-1)'
      do lm = 1,ndim
        write(3,140) B(lm,:)
      end do
    endif

! dChi = (1 - chi0*K)**(-1) * dChi0
! La recopie sur V fait gagner de l'espace memoire.
    allocate( V(ndim) )

    do i = 1,ndim
      V(:) = B(i,:)
      do j = 1,ndim
        B(i,j) = sum( V(:) * A(:,j) )
      end do
    end do

    deallocate( V )

    if( icheck > 2 ) then
      write(3,'(/A)') '  Chi'
      do lm = 1,ndim
        write(3,140) B(lm,:)
      end do
    endif

    if( icheck > 1 ) write(3,150)

    if( Optic ) then
      nlms_gg = nlms
    else
      nlms_gg = nd3
    endif

! Remplissage de Chi pour les lm utiles (jusqu'a nlm_probe).
    i1 = 0
    do isp_g1 = 1,2    ! Spin etat initial entree
      do lms_g1 = 1,nlms_gg
        if( Optic ) then
          l_g1 = l_val(lms_g1,isp_g1)
          m = m_val(lms_g1,isp_g1)
          iso = i_val(lms_g1,isp_g1)
          if( mod(l_g1,2) /= imparite .and. imparite /= 2 ) cycle
          lmsg1 = nspino *( l_g1**2 + l_g1 + m ) + iso
        else
          if( abs( coef_g(lms_g1,isp_g1) ) < eps10 ) cycle
          lmsg1 = lms_g1
          l_g1 = lseuil
        endif

        do isp_f1 = 1,2         ! Spin du premier etat final

          if( (isp_f1 /= isp_g1 .and. is_dipmag == 1) .or. (isp_f1 == isp_g1 .and. is_dipmag == 2) ) cycle

          do lms_f1 = 1,nlms
            l_f1 = l_val(lms_f1,isp_f1)
            m_f1 = m_val(lms_f1,isp_f1)
            iso = i_val(lms_f1,isp_f1)
            if( mod(l_f1,2) /= imparite .and. imparite /= 2 ) cycle
            lms1 = nspino *( l_f1**2 + l_f1 + m_f1 ) + iso
            i1 = i1 + 1

            if( l_g1 > lmax_probe .or. l_f1 > lmax_probe ) cycle

            i2 = 0
            do isp_g2 = 1,2
              do lms_g2 = 1,nlms_gg
                if( Optic ) then
                  l_g2 = l_val(lms_g2,isp_g2)
                  m = m_val(lms_g2,isp_g2)
                  iso = i_val(lms_g2,isp_g2)
                  if( mod(l_g2,2) /= imparite .and. imparite /= 2 ) cycle
                  lmsg2 = nspino *( l_g2**2 + l_g2 + m ) + iso
                else
                  if( abs( coef_g(lms_g2,isp_g2) ) < eps10 ) cycle
                  lmsg2 = lms_g2
                  l_g2 = lseuil
                endif

                do isp_f2 = 1,2       ! Spin du deuxieme etat final

                  if( (isp_f2 /= isp_g2 .and. is_dipmag == 1) .or. (isp_f2 == isp_g2 .and. is_dipmag == 2) ) cycle
                  if( (isp_f1 == isp_f2 .and. isp_g1 /= isp_g2) .or. (isp_f1 /= isp_f2 .and. isp_g1 == isp_g2) ) cycle

                  do lms_f2 = 1,nlms
                    l_f2 = l_val(lms_f2,isp_f2)
                    m_f2 = m_val(lms_f2,isp_f2)
                    iso = i_val(lms_f2,isp_f2)
                    if( mod(l_f2,2) /= imparite .and. imparite /= 2 ) cycle
                    i2 = i2 + 1
                    lms2 = nspino *( l_f2**2 + l_f2 + m_f2 ) + iso

                    if( l_g2 > lmax_probe .or. l_f2 > lmax_probe ) cycle

                    Chi(lms1,lms2,lmsg1,lmsg2,isp_f1,isp_f2,is_dipmag) = B(i1,i2)
                    if( icheck > 1 .and. abs( B(i1,i2) ) > eps10 ) write(3,160) l_f1, m_f1, lmsg1, isp_g1, l_f2, &
                                                                          m_f2, lmsg2, isp_g2, lms1, lms2, B(i1,i2)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    deallocate( A, B )

  end do  ! fin boucle is_dipmag

  if( icheck > 1 ) call write_Chi(Chi,Energ,First_E,ninitl,icheck,nlm_probe,nomfich,ns_dipmag,nspino)

  return
  100 format(/' ---- Cal_chi -------',100('-'))
  110 format(/' Energ =',f10.3,' eV')
  120 format(1p,9(1x,2e10.2))
  130 format(500(7x,'(',2(i2,','),i2,')',5x))
  135 format(500(1x,'(',2(i2,','),i2,')'))
  140 format(1p,1000e11.3)
  150 format(' l1 m1 g1 s1 l2 m2 g2 s2  lms1 lms2              chi')
  160 format(8i3,2i5,3x,1p,2e13.5)

end

!**********************************************************************

subroutine Write_Chi(Chi,Energ,First,ninitl,icheck,nlm_probe,nomfich,ns_dipmag,nspino)

  use declarations
  implicit none

  character(len=132),intent(in)::  nomfich

  integer,intent(in):: ninitl, icheck, nlm_probe, ns_dipmag, nspino

  complex(kind=db),dimension(nlm_probe*nspino,nlm_probe*nspino,ninitl,ninitl,2,2,ns_dipmag),intent(in):: Chi

  logical:: Diag, First

  real(kind=db),intent(in):: Energ

  character(len=132)::  chi_conv

  integer l, isp1, isp2, ist1, ist2, lm1, lm2

  if( icheck > 2 ) then
    Diag = .false.
  else
    Diag = .true.
  endif

! Output file name
  chi_conv = nomfich
  l = len_trim( chi_conv )
  chi_conv(l+1:l+8) = '_chi.txt'

  if( First ) then
    open(31, file = chi_conv)
    if( diag ) then
      write(31,110) ((( mod(lm1,10), isp1, mod(ist1,10) , lm1 = 1,nspino*nlm_probe), isp1 = 1,2), ist1 = 1,ninitl)
    else
      write(31,120) (((((( mod(lm1,10), mod(lm2,10), isp1, isp2, ist1, ist2, &
        lm2 = 1,nspino*nlm_probe), lm1 = 1,nspino*nlm_probe), isp2 = 1,2), isp1 = 1,2), ist2 = 1,ninitl),ist1 = 1,ninitl)
    endif
  else
    open(31, file = chi_conv, position='append')
  endif

  if( Diag ) then
    write(31,130) Energ*rydb, ((( chi(lm1,lm1,ist1,ist1,isp1,isp1,1), lm1 = 1,nspino*nlm_probe), &
         isp1 = 1,2), ist1 = 1,ninitl)
  else
    write(31,130) Energ*rydb, (((((( chi(lm1,lm2,ist1,ist2,isp1,isp2,1), &
        lm2 = 1,nspino*nlm_probe),lm1 = 1,nspino*nlm_probe), isp2 = 1,2), isp1 = 1,2), ist2 = 1,ninitl), ist1 = 1,ninitl)
  endif

  Close(31)

  return
  110 format('  Energy  ',5000('  (lm=',i1,',isp=',i1,',i=',i1,')  Im'))
  120 format('  Energy  ',5000('  (',i1,',',i1,',',i1,',',i1,',',i1, ',',i1,')    Im '))
  130 format(f10.3,1p,11664e11.3)
end


!***********************************************************************

! Calculation of Gaunt coefficient: Int( Y(l1,m1)*Y(l2,m2)Y(l3,m3) dOmega  )
! with Y(l2,m2) and Y(l3,m3) complex and  Y(l1,m1) real or complex

function Gaunt_xcc(l1,m1,l2,m2,l3,m3,Ylm_comp)

  use declarations
  implicit none

  complex(kind=db):: Gaunt_xcc

  integer, intent(in):: l1, m1, l2, m2, l3, m3
  logical, intent(in):: Ylm_comp

  real(kind=db) gr, gi, gauntcp

  if( Ylm_comp .or. m1 == 0 ) then
    Gaunt_xcc = cmplx( gauntcp(l1,m1,l2,m2,l3,m3), 0._db, db )
  else if( m1 > 0 ) then
! Gauntcp calculant Gaunt pour le complexe conjugue, on appele avec
! (-1)**m Y(l1,-m1) = Y(l1,m1)*
    gr = ( (-1)**m1 * gauntcp(l1,-m1,l2,m2,l3,m3) + gauntcp(l1,m1,l2,m2,l3,m3) ) / sqrt(2._db)
    gi = 0._db
    Gaunt_xcc = cmplx(gr, gi, db)
  else
! Gauntcp calculant Gaunt pour le complexe conjugue, on appele avec
! (-1)**m Y(l1,-m1) = Y(l1,m1)*
    gr = 0._db
    gi = ( gauntcp(l1,-m1,l2,m2,l3,m3) - (-1)**m1 * gauntcp(l1,m1,l2,m2,l3,m3) ) / sqrt(2._db)
    Gaunt_xcc = cmplx(gr, gi, db)
  end if

  return
end

!***********************************************************************

! Transformation harmo comp toward Harmo real
! The inverse transformation is the conjugate complex of the transpose

subroutine Cal_Trans(nlm,Trans)

  use declarations
  implicit none

  integer:: is, l1, l2, lm1, lm2, m1, m2, nlm

  complex(kind=db):: r2_r, r2_i
  complex(kind=db),dimension(nlm,nlm):: Trans

  real(kind=db):: r2

  Trans(:,:) = (0._db, 0._db)

  r2 = 1 / sqrt(2._db)
  r2_r = cmplx( r2,    0._db, db)
  r2_i = cmplx( 0._db, r2,    db)

  lm1 = 0
  boucle_l1: do l1 = 0,100
    do m1 = -l1,l1
      lm1 = lm1 + 1
      if( lm1 > nlm ) exit boucle_l1
      is = (-1)**m1

      lm2 = 0
      boucle_l2: do l2 = 0,100
        do m2 = -l2,l2
          lm2 = lm2 + 1
          if( l1 /= l2 ) cycle
          if( lm2 > nlm ) exit boucle_l2

          if( m1 == m2 ) then

            if( m1 == 0 ) then
              Trans(lm1,lm2) = (1._db,0._db)
            elseif( m1 > 0 ) then
              Trans(lm1,lm2) = r2_r
            else
              Trans(lm1,lm2) = - is * r2_i
            endif

          elseif( m1 == - m2 ) then

            if( m1 > 0 ) then
              Trans(lm1,lm2) = is * r2_r
            else
              Trans(lm1,lm2) = r2_i
            endif

          endif

        end do
      end do boucle_l2

    end do
  end do boucle_l1

  return
end

!***********************************************************************

! Only two spin indexes because on both side, occupied and non occupied states have same spin
! (isp_f1 = isp_g1 and isp_f2 = isp_g2)

subroutine fxcorr(alfpot,fxc,icheck,magnetic,nr,nspin,r,rhoato_abs,rsato)

  use declarations
  implicit none

  integer, intent(in):: icheck, nr, nspin

  logical, intent(in):: magnetic

  real(kind=db), intent(in):: alfpot
  real(kind=db),dimension(nr),intent(in):: r, rsato
  real(kind=db),dimension(nr,nspin),intent(in):: rhoato_abs
  real(kind=db),dimension(nr,2,2),intent(out):: fxc

  integer ir, isp, isp1, isp2

  real(kind=db):: f_vonbarth, fprime_vonbarth
  real(kind=db):: a, b, c, c_p, c_f, d, e, f, f1, f2, f3, f4, f5, fac, r_p, r_f, rsa, qtr, tr, x, xx

  if( icheck > 2 ) write(3,100)

  fxc(:,:,:) = 0._db

  if( alfpot > eps4 ) then

! Xalpha potential
    fac = - 2 * pi * alfpot / 3
    do isp = 1,2
      fxc(:,isp,isp) = fac * rsato(:)**2
    end do

  else if( alfpot < eps4 ) then

! Pour les valeurs de c_p et r_p on garde les valeurs non polarise
! choisie aussi par Moruzzi Janak et William (1978). Pour r_f et c_f
! on prend aussi leurs parametres plutot que les originaux de Von Barth
! qui sont : c_p = 0.0504, r_p = 30., c_f = 0.0254, r_f = 75.
    c_p =  0.045_db
    r_p = 21._db
    c_f = 0.0254_db
    r_f = 75._db

    tr  =  1._db / 3._db
    qtr = 4._db / 3._db

    f1 = ( 36._db / pi**2 )**tr
    f2 = 4 * pi / 9
    f3 = (2._db)**(-tr)
    f4 = qtr / ( 1 - f3 )
    f5 = 1 / ( 1 - f3 )

    do isp1 = 1,2
      do isp2 = 1,2
        do ir = 1,nr
          if( Magnetic ) then
            x = rhoato_abs(ir,isp2) / sum( rhoato_abs(ir,1:nspin) )
          else
            x = 0.5_db
          end if
          rsa = rsato(ir)

          if( isp1 == isp2 ) then
            xx = 1 - x
            a = f1 * ( x**(-2*tr) * xx + x**tr ) / rsa
          else
            xx = - x
            a = 0._db
          endif

          b = c_p * r_p / ( rsa + r_p)

          c = f4 * rsa * ( x**qtr + (1-x)**qtr - x**tr) * ( fprime_vonbarth(rsa/r_f) * (c_f/r_f) &
              - fprime_vonbarth(rsa/r_p) * (c_p/r_p) )

          d = f5 * ( x**qtr + (1-x)**qtr - f3 ) * ( c_f*r_f / ( r_f + rsa ) -  c_p*r_p / ( r_p + rsa ) )

          e = f4 * ( x**(-2*tr) + 4*x**tr - 4*(1-x)**tr ) * xx * ( c_p * f_vonbarth(rsa/r_p) - c_f * f_vonbarth(rsa/r_f) )

          f = 4 * f5 * xx * ( (1-x)**tr - x**tr ) * ( c_p * log(1 + r_p/rsa) - c_f * log(1 + r_f/rsa) )

          fxc(ir,isp1,isp2) = - ( a + b + c + d + e + f ) * f2 * rsa**3

        end do
      end do
    end do
  end if

  if( icheck > 2 ) then
    if( magnetic ) then
      write(3,110)
    else
      write(3,120)
    end if
    do ir = 1,nr
      write(3,130) r(ir)*bohr, rsato(ir)*bohr, ( rhoato_abs(ir,isp)*r(ir)**2, isp = 1,nspin ), fxc(ir,:,:)
    end do
  end if

  return
  100 format(/' ---- Fxcorr -------',100('-'))
  110 format(/'     Radius       Rsato   rhoato_up*r**2 rhoato_dn*r**2  fxc_uu      fxc_ud       fxc_du       fxc_dd')
  120 format(/'     Radius       Rsato     rhoato*r**2     fxc')
  130 format(1p,8e13.5)
end

!***********************************************************************

! Calculation of Integral( Y(l1,m1)* Y(l2,m2)* Y(l3,m3) Y(l4,m4) dOmega )

! One uses: Y(L2)* Y(L3) = Sum_L  Gaunt(L3,L2,L) Y(L)
! with L = (l,m)

! G4(L1,L2,L3,L4) = Sum_L  G3(L2,L3,L) x G3(L4,L1,L)

function Gaunt4Y(l1,m1,l2,m2,l3,m3,l4,m4)

  use declarations
  implicit none

  integer, intent(in):: l1, m1, l2, m2, l3, m3, l4, m4
  integer:: l, lmin, lmax, m

  real(kind=db):: Gaunt4Y, Gauntcp

  Gaunt4Y = 0._db

  lmin = max( abs(l2-l3), abs(l1-l4) )
  lmax = min( l2+l3, l1+l4 )

  do l = lmin,lmax,2
    do m = -l,l
      Gaunt4Y = Gaunt4Y + Gauntcp(l2,m2,l3,m3,l,m) * Gauntcp(l4,m4,l1,m1,l,m)
    end do
  end do

  return
end

!***********************************************************************

! Calculation of Integral( Y(l1,m1) Y(l2,m2) Y(l3,m3) Y(l4,m4) dOmega )
! with the 4 harmonics real

function Gaunt4Y_r(l1,m1,l2,m2,l3,m3,l4,m4)

  use declarations
  implicit none

  real(kind=db), parameter:: sqrt2_i = 1._db / sqrt( 2._db )

  complex(kind=db):: c1, c2, c3, c4

  integer:: l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4, s1, s2, s3, s4

  real(kind=db):: Gaunt4Y, Gaunt4Y_r

  Gaunt4Y_r = 0._db

  n1 = m1
  do s1 = 1,-1,-2

    if( s1 == 1 ) then
      if( m1 < 0 ) then
        if( mod(m1,2) == 0 ) then
          c1 = cmplx( 0._db, sqrt2_i, db )
        else
          c1 = cmplx( 0._db, - sqrt2_i, db )
        endif
      elseif( m1 == 0 ) then
        c1 = cmplx( 1._db, 0._db, db )
      else
        c1 = cmplx( sqrt2_i, 0._db, db )
      endif
    else
      if( m1 < 0 ) then
        c1 = cmplx( 0._db, - sqrt2_i, db )
      else
        if( mod(m1,2) == 0 ) then
          c1 = cmplx( sqrt2_i, 0._db, db )
        else
          c1 = cmplx( - sqrt2_i, 0._db, db )
        endif
      endif
    endif

    n2 = m2
    do s2 = 1,-1,-2

      if( s2 == 1 ) then
        if( m2 < 0 ) then
          if( mod(m2,2) == 0 ) then
            c2 = cmplx( 0._db, sqrt2_i, db )
          else
            c2 = cmplx( 0._db, - sqrt2_i, db )
          endif
        elseif( m2 == 0 ) then
          c2 = cmplx( 1._db, 0._db, db )
        else
          c2 = cmplx( sqrt2_i, 0._db, db )
        endif
      else
        if( m2 < 0 ) then
          c2 = cmplx( 0._db, - sqrt2_i, db )
        else
          if( mod(m2,2) == 0 ) then
            c2 = cmplx( sqrt2_i, 0._db, db )
          else
            c2 = cmplx( - sqrt2_i, 0._db, db )
          endif
        endif
      endif

      n3 = m3
      do s3 = 1,-1,-2

        if( s3 == 1 ) then
          if( m3 < 0 ) then
            if( mod(m3,2) == 0 ) then
              c3 = cmplx( 0._db, sqrt2_i, db )
            else
              c3 = cmplx( 0._db, - sqrt2_i, db )
            endif
          elseif( m3 == 0 ) then
            c3 = cmplx( 1._db, 0._db, db )
          else
            c3 = cmplx( sqrt2_i, 0._db, db )
          endif
        else
          if( m3 < 0 ) then
            c3 = cmplx( 0._db, - sqrt2_i, db )
          else
            if( mod(m3,2) == 0 ) then
              c3 = cmplx( sqrt2_i, 0._db, db )
            else
              c3 = cmplx( - sqrt2_i, 0._db, db )
            endif
          endif
        endif

        n4 = m4
        do s4 = 1,-1,-2

          if( s4 == 1 ) then
            if( m4 < 0 ) then
              if( mod(m4,2) == 0 ) then
                c4 = cmplx( 0._db, sqrt2_i, db )
              else
                c4 = cmplx( 0._db, - sqrt2_i, db )
              endif
            elseif( m4 == 0 ) then
              c4 = cmplx( 1._db, 0._db, db )
            else
              c4 = cmplx( sqrt2_i, 0._db, db )
            endif
          else
            if( m4 < 0 ) then
              c4 = cmplx( 0._db, - sqrt2_i, db )
            else
              if( mod(m4,2) == 0 ) then
                c4 = cmplx( sqrt2_i, 0._db, db )
              else
                c4 = cmplx( - sqrt2_i, 0._db, db )
              endif
            endif
          endif

          Gaunt4Y_r = Gaunt4Y_r + real( conjg( c1 * c2 ) * c3 * c4, db ) * Gaunt4Y(l1,n1,l2,n2,l3,n3,l4,n4)

          if( m4 == 0 ) exit
          n4 = - m4
        end do

        if( m3 == 0 ) exit
        n3 = - m3
      end do

      if( m2 == 0 ) exit
      n2 = - m2
    end do

    if( m1 == 0 ) exit
    n1 = - m1
  end do

  return
end

!***********************************************************************

function fprime_vonbarth(x)

  use declarations
  implicit none
  real(kind=db):: fprime_vonbarth, x

  fprime_vonbarth = 3*x**2 * log(1 + 1/x) - 1/x -3*x + 1.5_db

  return
end

!***********************************************************************

! Main for TDDFT in optics

! One calculates dChi(E_g,hbar.omega) = ( 1 - Chi_0(hbar.omega)*K(E_g,E_f) )^-1 dChi_0(E_g,hbar.omega)
! with hbar.omega = E_f - E_g
! Chi(hbar.omega) = Sum_ie_g ( dChi(E_g,hbar.omega) ) is not explicitely calculated

subroutine main_tddft_optic(alfpot,angxyz,Allsite,Atomic_scr,axyz,Classic_irreg,coef_g, &
        Cartesian_tensor,Core_resolved,Dafs,Dafs_bio,Dipmag, &
        dv0bdcF,E_cut,E_cut_imp,E_cut_man,Eclie,Eneg, &
        Energ_s,Energphot,Ephot_min,Eseuil,f_avantseuil,Full_potential,Full_self_abs, &
        Gamma_tddft,hkl_dafs,Hubb_a,Hubb_d,icheck,iabsorig,is_g,isigpi,isymeq, &
        jseuil,Kern_fac,Kern_fast,ldip,Length_abs,lmax_pot,lmaxabs_t,lmoins1,loct,lplus1,lqua,lseuil, &
        ltypcal,m_g,m_hubb,Magnetic,Matper,Moyenne,mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq, &
        msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
        n_multi_run,n_oo,n_rel,n_rout,n_tens_max, &
        natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,nenerg_tddft,ninit1,ninitl,ninitl_out,ninitlv,nlm_pot,nlmamax, &
        nomabs,nomfich,nomfich_cal_tddft_conv,nomfich_s, &
        nphi_dafs,nphim,npldafs,nplr,nplrm,nr,nrm,nseuil,nspin,nspino,nspinp, &
        numat,nxanout,Octupole,pdp,phdafs,phdf0t,phdt,pol,poldafse,poldafss, &
        psii,r,Relativiste,Renorm,rhoato_abs,Rmtg,Rmtsd, &
        rot_atom_abs,Rot_int,RPALF,rsato,Self_abs,Solsing_only, &
        Spherical_signal,Spherical_tensor,Spinorbite,Taull_tdd,Taux_eq,Time_rout,V_intmax,V_hubb,V0muf, &
        Vcato,Vec,Vecdafse,Vecdafss,Vhbdc,Volume_maille,VxcbdcF, &
        Vxcato,Workf,Ylm_comp_e)

  use declarations
  implicit none
  include 'mpif.h'

  integer, parameter:: nq_nrixs = 0
  integer, parameter:: l0_nrixs = 0
  integer, parameter:: lmax_nrixs = 0

  integer:: i, iabsorig, icheck_s, ie, ie_computer, ie_e, ie_g, ief, ip_max, ip0, &
    iso1, iso2, isp, isp1, isp2, j, je, jef, jseuil, lm1, lm2, lmax, lmax_pot, &
    lmax_probe, lmaxabs_t, lseuil, m_hubb, mpinodes, mpirank, mpirank0, &
    n_rel, n_Ec, n_multi_run, n_oo, n_rout, n_stk, n_tens_max, n_V, natomsym, nbseuil, &
    ncolm, ncolr, ncolt, nd3, nenerg, nenerg_s, nenerg_tddft, nge, ninit1, ninitl, ninitl_out, &
    ninitlv, nlm, nlm_fp, nlm_pot, nlm_probe, nlm_p_fp, nlmamax, nlms, nlms_g, nlms_f, &
    nphim, npldafs, nplr, nplrm, nr, nr_zet, nrm, ns_dipmag, &
    nseuil, nspin, nspino, nspinp, nspint, numat, nxanout

  integer, dimension(30):: icheck
  integer, dimension(natomsym):: isymeq
  integer, dimension(ninitl):: is_g
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(ninitl,2):: m_g
  integer, dimension(2,npldafs):: isigpi
  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3)::  msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi

  character(len=132):: nomfich, nomfich_s, nomfich_cal_convt, nomfich_cal_tddft_conv
  character(len=13), dimension(nplrm):: ltypcal
  character(len=length_word), dimension(ncolm):: nomabs

  complex(kind=db):: f_avantseuil
  complex(kind=db), dimension(3,nplrm):: pol
  complex(kind=db), dimension(natomsym,npldafs):: phdafs
  complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdt
  complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss
  complex(kind=db), dimension(3,3,n_rel,ninitl_out,0:mpinodes-1):: secdd, secdd_m, secdd_t, secdd_m_t
  complex(kind=db), dimension(3,3,ninitl_out,0:mpinodes-1):: secmd, secmd_m, secmm, secmm_m, &
                                                             secmd_t, secmd_m_t, secmm_t, secmm_m_t
  complex(kind=db), dimension(3,3,3,ninitl_out,0:mpinodes-1):: secdq, secdq_m, secdq_t, secdq_m_t
  complex(kind=db), dimension(3,3,3,3,ninitl_out,0:mpinodes-1):: secdo, secdo_m, secqq, secqq_m, secdo_t, secdo_m_t, secqq_t, &
                                                                 secqq_m_t
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitl_out,0:mpinodes-1):: secoo, secoo_m, secoo_t, secoo_m_t
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nlmamax,nspinp) :: Taull_tdd

  complex(kind=db), dimension(:,:), allocatable:: B_stk
  complex(kind=db), dimension(:,:,:,:), allocatable:: V_hubb_t
  complex(kind=db), dimension(:,:,:,:,:), allocatable::  rof0_e
  complex(kind=db), dimension(:,:,:,:), allocatable:: Chi_0
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: dChi_0, dChi_0_i
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: Chi

  logical:: Allsite, Atomic_scr, Cartesian_tensor, Classic_irreg, Core_resolved, Dafs, Dafs_bio, Dipmag, &
    E_cut_man, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Eneg, Energphot, First_E, First_K, FDM_comp, Final_optic, Final_tddft, &
    Full_potential, Full_self_abs, Gamma_tddft, Green, Green_int, &
    Hubb_a, Hubb_d, Kern_fast, lmoins1, lplus1, M_depend, M1M1, Magnetic, Matper, &
    Moyenne, No_diag, NRIXS, Octupole, Optic, Radial_comp, &
    Relativiste, Renorm, RPALF, Self_abs, Solsing, Solsing_only, Spherical_signal, &
    Spherical_tensor, Spinorbite, Xan_atom, Ylm_comp, Ylm_comp_e

  logical, dimension(10):: Multipole

  real(kind=sg) time

  real(kind=db):: alfpot, Delta_E, E_cut, E_cut_imp, E_cut_tddft, Eclie, Ephot_min, &
     Kern_fac, Rmtg, Rmtsd, V_intmax, V0muf, Vhbdc, Volume_maille, Workf

  real(kind=db), dimension(1):: Energ_u
  real(kind=db), dimension(3):: angxyz, axyz
  real(kind=db), dimension(8):: Time_loc
  real(kind=db), dimension(n_rout):: Time_rout
  real(kind=db), dimension(nspin):: dv0bdcF, VxcbdcF
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitl_out):: Epsii, sec_atom
  real(kind=db), dimension(nr):: r, rsato
  real(kind=db), dimension(natomsym) :: Taux_eq
  real(kind=db), dimension(3,3):: rot_atom_abs, Rot_int
  real(kind=db), dimension(nplrm,2) :: pdp
  real(kind=db), dimension(3,npldafs):: hkl_dafs
  real(kind=db), dimension(npldafs):: Length_abs
  real(kind=db), dimension(3,nplrm) :: vec
  real(kind=db), dimension(ninitl,2):: coef_g
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nspin):: rhoato_abs
  real(kind=db), dimension(n_tens_max,0:natomsym):: Int_tens
  real(kind=db), dimension(nr,nlm_pot):: Vcato
  real(kind=db), dimension(nr,2,2):: fxc
  real(kind=db), dimension(3,npldafs,nphim):: Vecdafse, Vecdafss
  real(kind=db), dimension(nr,nlm_pot,nspin) :: Vxcato
  real(kind=db), dimension(nq_nrixs,ninitl_out,0:mpinodes-1):: S_nrixs, S_nrixs_m
  real(kind=db), dimension(nq_nrixs,l0_nrixs:lmax_nrixs,ninitl_out,0:mpinodes-1):: S_nrixs_l, S_nrixs_l_m

  real(kind=db), dimension(:), allocatable:: Ecinetic_e, Eimag, Eimag_t, En, Energ, Enervide
  real(kind=db), dimension(:,:), allocatable:: Ecinetic, V0bdc
  real(kind=db), dimension(:,:,:,:), allocatable:: Kern, Vrato
  real(kind=db), dimension(:,:,:,:,:,:), allocatable:: roff0_ph, zet
  real(kind=db), dimension(:,:,:,:,:,:,:), allocatable::  roff_rr, roff0, roff0t

  if( icheck(1) > 0 ) write(3,100)

! Green complet avec partie reelle et imaginaire (en cas de Gamma non nul)
! Si Green_int, f' et f" sont directement calcule et l'integrale sur l'energie n'est pas a faire dans convolution.
  Green_int = Gamma_tddft

  Xan_atom = .false.
  Optic = .true.
  NRIXS = .false.
  Final_optic = Optic
  Final_tddft = .true.
  First_E = .true.
  First_K = .true.
  FDM_comp = .true.
  Green = .true.
  Epsii(:) = 0._db
  n_tens_max = 0

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  ns_dipmag = 4
  if( Dipmag ) then
    ip0 = 0
  else
    ip0 = 1
  endif
  if( Octupole ) then
    ip_max = 3
  else
    ip_max = 2
  endif

  sec_atom(:) = 0._db

! All calculations are now using complex harmonics
! For calculations with spinorbit, they where necessarily already with complex harmonics
!  Ylm_comp = .true.

  Ylm_comp = Ylm_comp_e

!  if( .not. Ylm_comp_e ) then
!    Call Tr_Taull( icheck(22), nlmamax, nenerg_s, nspinp, Taull_tdd )

!    if( Hubb_a ) then
!      l = m_hubb
!      allocate( V(-l:l,-l:l) )
!      allocate( Trans(-l:l,-l:l) )
!      Call Cal_trans_l(l,Trans)
!      do isp = 1,nspinp
!        V(-l:l,-l:l) = V_hubb(-l:l,-l:l,isp,isp)
!        V = Matmul( Conjg( Transpose(Trans) ), Matmul( V, Trans ))
!        V_hubb(-l:l,-l:l,isp,isp) = V(-l:l,-l:l)
!      end do
!      deallocate( Trans, V )
!    endif
!  endif

! The Taull_tdd are "efficient" scattering amplitudes. One multiplies by "i" in order the imaginary part be the absorption
  Taull_tdd(:,:,:,:,:) = img * Taull_tdd(:,:,:,:,:)

  allocate( V_hubb_t(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp) )
  V_hubb_t(:,:,:,:) = V_hubb(:,:,:,:)

  if( E_cut_man ) then
    E_cut_tddft = E_cut_imp
  else
    E_cut_tddft = E_cut
  endif

  n_Ec = 2  ! Energies at E_F +/- hbar*omega
  n_V = 1   ! Potential does not depend on energy in optics
  jef = n_Ec

  allocate( En(n_Ec) )

! Index of Fermi level (at this energy index, states are not occupied)
  do ie = 1,nenerg_s-1
    if( Energ_s(ie) > E_cut_tddft - eps10  ) exit
  end do
  ief = ie

  lmax = lmaxabs_t
  lmax_probe = lmaxabs_t
  nlm_probe = ( lmax_probe + 1 )**2
  nd3 = nlm_probe*nspino

  No_diag = Full_potential .or. ( Hubb_a .and. .not. Hubb_d )
  M_depend = Spinorbite .or. Hubb_a .or. No_diag

  if( No_diag ) then
    nlm = ( lmax + 1 )**2
    nlm_fp = nlm
    nlm_p_fp = nlm_probe
  elseif( M_depend ) then
    nlm = ( lmax + 1 )**2
    nlm_fp = 1
    nlm_p_fp = 1
  else
    nlm = lmax + 1
    nlm_fp = 1
    nlm_p_fp = 1
  endif

  if( Spinorbite ) then
    nlms = ( lmax + 1 ) * ( 2 * lmax + 1 )
  else
    nlms = ( lmax + 1 )**2
  endif

  nspint = nspin

  nlms_f = nspint * nlms
  nlms_g = nspint * nlms

  allocate( Ecinetic_e(nspin) )
  allocate( V0bdc(nspin,n_V) )

! Elaboration of the energy grid for Tddft
  call dim_grille_optic(E_cut_tddft,Energ_s,icheck(22),mpirank0,nenerg,nenerg_s)
  allocate( Energ(nenerg) )
  call grille_optic(E_cut_tddft,Energ,Energ_s,icheck(22),nenerg,nenerg_s)

  Solsing = .false.

  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    Time_loc(1) = real(time,db)
  endif

  allocate( Vrato(nr,nlm_pot,nspin,n_V) )

  do ie_e = 1,n_V
    do isp = 1,nspin
      Vrato(1:nr,:,isp,ie_e) = Vcato(1:nr,:) + Vxcato(1:nr,:,isp)
      V0bdc(isp,ie_e) = Vhbdc + VxcbdcF(isp) + dV0bdcF(isp)
    end do
  end do

  allocate( roff_rr(nlm,nlm,nlm_fp,nlm_fp,nspinp**2,nspino**2,0:0) )
  allocate( roff0(ief-1,ief:nenerg_s,nlm,nlm,nspinp,nspinp,nspino**2) )

! Calculation of radial integrals for Chi_0 ( roff0 ) ----------
  allocate( Ecinetic(nspin,nenerg_s) )
  allocate( Enervide(nenerg_s) )
! One keeps the imaginary energy as zero
  allocate( Eimag(nenerg_s) )
  Eimag(:) = 0._db
  nr_zet = 0
  allocate( zet(nr_zet,nlm,nlm_fp,nspinp,nspino,nenerg_s) )

  do ie = 1,nenerg_s
    Enervide(ie) = Energ_s(ie) - Workf
    Ecinetic(:,ie) = Enervide(ie) - V0bdc(:,1)
    if( .not. Eneg ) Ecinetic(:,ie) = max( Ecinetic(:,ie), Eclie )
  end do

  call radial_optic(Classic_irreg,Ecinetic,Eimag,Energ,Enervide,Full_potential,Hubb_a,Hubb_d,icheck(22),ief,0,0, &
       lmax,lmax_pot,m_depend,m_hubb,nenerg,nenerg_s,n_V,nlm_pot,nlm,nlm_fp,No_diag,nr,nr_zet,nspin,nspino,nspinp,numat,r, &
       Relativiste,Renorm,Rmtg,Rmtsd,roff_rr,roff0,Solsing,Spinorbite,V_hubb,V_intmax,V0bdc,Vrato,Ylm_comp,zet)

  deallocate( Ecinetic, Eimag, Enervide, zet )
! End of calculation of radial integrals for Chi0 --------------

  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    Time_loc(2) = real(time,db)
  endif
  Time_rout(14) = Time_rout(14) + Time_loc(2) - Time_loc(1) ! Radial

  allocate( Chi(nlm_probe*nspino,nlm_probe*nspino,nd3,nd3,2,2,ns_dipmag) )
  allocate( Chi_0(nlms_f,nlms_f,nlms_g,nlms_g) )
  allocate( dChi_0(ief-1,nlms_f,nlms_f,nlms_g,nlms_g) )
  allocate( dChi_0_i(ief-1,nlms_f,nlms_f,nlms_g,nlms_g) )
  allocate( Kern(nlms_f,nlms_f,nlms_g,nlms_g) )
  allocate( Ecinetic(nspin,n_Ec) )
  allocate( Eimag_t(n_Ec) )
  allocate( Enervide(n_Ec) )
  Eimag_t(:) = 0._db
  nr_zet = nr
  allocate( zet(nr_zet,nlm,nlm_fp,nspinp,nspino,nenerg_s) )

  if( Kern_fast ) then
    n_stk = nlms_f * nlms_g
  else
    n_stk = 0
  endif
  allocate( B_stk(n_stk,n_stk) )

! On ne garde que la partie reelle de l'orbitale radiale zet. Verifier plus tard si correct
  Radial_comp = Eimag_t(1) > eps10 .or. ( Hubb_a .and. Ylm_comp )

 ! Calculation of exchange-correlation functional
  if( .not. RPALF ) call fxcorr(alfpot,fxc,icheck(24),Magnetic,nr,nspin,r,rhoato_abs,rsato)

  nge = ( nenerg - 1 ) / mpinodes + 1

! Loop over photon energy ---------------------------------------------------------------------------------------

  boucle_energ: do je = 1,nge

    ie = ( je - 1 ) * mpinodes + mpirank + 1

    if( Energ(ie) < Ephot_min - eps10 ) cycle

    if( ie <= nenerg ) then

      dChi_0(:,:,:,:,:) = (0._db,0._db)
      dChi_0_i(:,:,:,:,:) = (0._db,0._db)
      Chi_0(:,:,:,:) = (0._db,0._db)

! Calculation of Chi_0
      do ie_g = 1,ief-1  ! Loop over occupied state electron

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(2) = real(time,db)
        endif

        En(1) = Energ_s(ie_g)
        En(2) = Energ_s(ie_g) + Energ(ie)
        if( En(1) > E_cut_tddft - eps10 ) exit
        if( En(2) < E_cut_tddft - eps10 .or. En(2) > Energ_s(nenerg_s) + eps10 ) cycle

        do ie_e = 1,n_Ec
          Enervide(ie_e) = En(ie_e) - Workf
          Ecinetic(:,ie_e) = Enervide(ie_e) - V0bdc(:,1)
          if( .not. Eneg ) Ecinetic(:,ie_e) = max( Ecinetic(:,ie_e), Eclie )
        end do

        zet(:,:,:,:,:,:) = 0._db

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(2) = real(time,db)
        endif

        allocate( roff0t(jef-1,jef:n_Ec,nlm,nlm,nspinp,nspinp,nspino**2) )
        allocate( roff0_ph(nlm,nlm,nspinp,nspinp,nspino,nspino) )
        roff0_ph(:,:,:,:,:,:) = 0._db
        Energ_u(1) = Energ(ie)

        call radial_optic(Classic_irreg,Ecinetic,Eimag_t,Energ_u,Enervide,Full_potential,Hubb_a,Hubb_d,icheck(22),jef,0,0, &
          lmax,lmax_pot,m_depend,m_hubb,1,n_Ec,n_V,nlm_pot,nlm,nlm_fp,No_diag,nr,nr_zet,nspin,nspino,nspinp,numat,r, &
          Relativiste,Renorm,Rmtg,Rmtsd,roff_rr,roff0t,Solsing,Spinorbite,V_hubb_t,V_intmax,V0bdc,Vrato,Ylm_comp,zet)

        deallocate( roff0t )

        i = 0
        do isp1 = 1,nspinp
          do isp2 = 1,nspinp
            i = i + 1
            j = 0
            do iso1 = 1,nspino
              do iso2 = 1,nspino
                j = j + 1
                if( nlm_fp == 1 ) then
                  roff0_ph(:,:,isp1,isp2,iso1,iso2) = roff_rr(:,:,1,1,i,j,0)
                else
                  do lm1 = 1,nlm
                    do lm2 = 1,nlm_fp
                      roff0_ph(lm1,lm2,isp1,isp2,iso1,iso2) = sum( roff_rr(:,:,lm1,lm2,i,j,0) )
                    end do
                  end do
                endif
              end do
            end do
          end do
        end do

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(3) = real(time,db)
        endif
        Time_rout(14) = Time_rout(14) + Time_loc(3) - Time_loc(2) ! Radial

        if( ie_g == nenerg_s ) then
          delta_E = Energ_s(ie_g) - Energ_s(ie_g-1)
        elseif( ie_g == 1 ) then
          delta_E = Energ_s(ie_g+1) - Energ_s(ie_g)
        elseif(  Energ_s(ie_g+1) > E_cut_tddft ) then
          delta_E = E_cut_tddft - 0.5_db * ( Energ_s(ie_g) + Energ_s(ie_g-1) )
        else
          delta_E = 0.5_db * ( Energ_s(ie_g+1) - Energ_s(ie_g-1) )
        endif

        call Chi_0_opt(Chi_0, dChi_0, dChi_0_i, delta_E, E_cut_tddft, Energ(ie), Energ_s, Gamma_tddft, icheck(23), ie_g, ief, &
           je, lmax, M_depend, nenerg_s, nlm, nlmamax, nlms_f, nlms_g, nomfich, &
           nspino, nspinp, nspint, roff0, roff0_ph, Spinorbite, Taull_tdd, Ylm_comp)

        deallocate( roff0_ph )

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(4) = real(time,db)
        endif
        Time_rout(15) = Time_rout(15) + Time_loc(4) - Time_loc(3) ! Chi_0

      end do

      secdd(:,:,:,:,mpirank) = (0._db,0._db)
      secmm(:,:,:,mpirank) = (0._db,0._db)
      secmd(:,:,:,mpirank) = (0._db,0._db)
      secdq(:,:,:,:,mpirank) = (0._db,0._db)
      secqq(:,:,:,:,:,mpirank) = (0._db,0._db)
      secdo(:,:,:,:,:,mpirank) = (0._db,0._db)
      secdd_m(:,:,:,:,mpirank) = (0._db,0._db)
      secmm_m(:,:,:,mpirank) = (0._db,0._db)
      secmd_m(:,:,:,mpirank) = (0._db,0._db)
      secdq_m(:,:,:,:,mpirank) = (0._db,0._db)
      secqq_m(:,:,:,:,:,mpirank) = (0._db,0._db)
      secdo_m(:,:,:,:,:,mpirank) = (0._db,0._db)
      if( E3E3 ) secoo(:,:,:,:,:,mpirank) = (0._db,0._db)
      if( E3E3 ) secoo_m(:,:,:,:,:,mpirank) = (0._db,0._db)

! Second loop on occupied states electron energy
      do ie_g = 1,ief-1

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(2) = real(time,db)
        endif

        En(1) = Energ_s(ie_g)
        En(2) = Energ_s(ie_g) + Energ(ie)
        if( En(1) > E_cut_tddft - eps10 ) exit
        if( En(2) < E_cut_tddft - eps10 .or. En(2) > Energ_s(nenerg_s) + eps10 ) cycle

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(2) = real(time,db)
        endif

        if( .not. Kern_fast .or. First_K ) then

          if( First_K .and. Kern_fast ) then
            En(1) = Energ_s(ief-1)
            En(2) = Energ_s(ief-1) + Energ(ie)
          endif
          do ie_e = 1,n_Ec
            Enervide(ie_e) = En(ie_e) - Workf
            Ecinetic(:,ie_e) = Enervide(ie_e) - V0bdc(:,1)
            if( .not. Eneg ) Ecinetic(:,ie_e) = max( Ecinetic(:,ie_e), Eclie )
          end do

          allocate( roff0t(jef-1,jef:n_Ec,nlm,nlm,nspinp,nspinp,nspino**2) )
          Energ_u(1) = Energ(ie)
          zet(:,:,:,:,:,:) = 0._db

          call radial_optic(Classic_irreg,Ecinetic,Eimag_t,Energ_u,Enervide,Full_potential,Hubb_a,Hubb_d,icheck(22),jef,0,0, &
            lmax,lmax_pot,m_depend,m_hubb,1,n_Ec,n_V,nlm_pot,nlm,nlm_fp,No_diag,nr,nr_zet,nspin,nspino,nspinp,numat,r, &
            Relativiste,Renorm,Rmtg,Rmtsd,roff_rr,roff0t,Solsing,Spinorbite,V_hubb_t,V_intmax,V0bdc,Vrato,Ylm_comp,zet)

          deallocate( roff0t )
        endif

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(3) = real(time,db)
        endif
        Time_rout(14) = Time_rout(14) + Time_loc(3) - Time_loc(2) ! Radial

! Kernel calculation
        if( .not. Kern_fast .or. First_K ) &
          call kernel_optic(Atomic_scr,Energ(ie),First_E,fxc,icheck(24),Kern,Kern_fac,lmax,M_depend,n_Ec,nlm,nlm_fp,nlms_f,nlms_g, &
                  nomfich,nr,nspino,nspinp,nspint,r,Rmtsd,RPALF,Spinorbite,Ylm_comp,zet)

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(5) = real(time,db)
        endif

! Chi is in fact dChi(ie_g) = (1 - Chi_0.K)^-1 dChi0(ie_g)
        call Cal_Chi_opt(B_stk, Chi, Chi_0, dChi_0, dChi_0_i, Energ(ie), First_E, First_K, icheck(25), ie_g, ief, &
              Kern, Kern_fast, lmax_probe, lmax, n_stk, nd3, ninitl, nlm_probe, nlm_fp, nlms_f, nlms_g, &
              nomfich, ns_dipmag, nspino, nspint, Spinorbite, Ylm_comp)

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(6) = real(time,db)
        endif

! Calculation of cartesian tensor
        nenerg_tddft = 0
        allocate( rof0_e(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil) )
        icheck_s = icheck(20)

        call tenseur_car(Classic_irreg,coef_g,Core_resolved,Ecinetic, &
                Eimag_t(1),Energ(ie),Enervide,Eseuil,FDM_comp,Final_optic,Final_tddft,Full_potential,Green,Green_int,Hubb_a, &
                Hubb_d,icheck_s,ie,ip_max,ip0,is_g,lmax_probe,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
                n_Ec,n_oo,n_rel,n_V,nbseuil,ns_dipmag,nd3,nenerg_tddft,ninit1,ninitl,ninitl_out,ninitlv,nlm_pot,nlm_probe, &
                nlm_p_fp,nlmamax,nr,nrm,nspin,nspino,nspinp,numat,psii,r,Relativiste,Renorm,Rmtg,Rmtsd,rof0_e,rot_atom_abs, &
                Rot_int,secdd_t,secdd_m_t,secdo_t,secdo_m_t,secdq_t,secdq_m_t,secmd_t,secmd_m_t,secmm_t,secmm_m_t,secoo_t, &
                secoo_m_t,secqq_t,secqq_m_t,Solsing,Solsing_only,Spinorbite,Chi,.true.,V_hubb_t,V_intmax,V0bdc,Vrato,Ylm_comp)

        deallocate( rof0_e )

        if( E1E1 ) secdd(:,:,:,:,mpirank) = secdd(:,:,:,:,mpirank) + secdd_t(:,:,:,:,mpirank)
        if( E1E1 ) secdd_m(:,:,:,:,mpirank) = secdd_m(:,:,:,:,mpirank) + secdd_m_t(:,:,:,:,mpirank)
        if( M1M1 ) secmm(:,:,:,mpirank) = secmm(:,:,:,mpirank) + secmm_t(:,:,:,mpirank)
        if( M1M1 ) secmm_m(:,:,:,mpirank) = secmm_m(:,:,:,mpirank) + secmm_m_t(:,:,:,mpirank)
        if( E1M1 ) secmd(:,:,:,mpirank) = secmd(:,:,:,mpirank) + secmd_t(:,:,:,mpirank)
        if( E1M1 ) secmd_m(:,:,:,mpirank) = secmd_m(:,:,:,mpirank) + secmd_m_t(:,:,:,mpirank)
        if( E1E2 ) secdq(:,:,:,:,mpirank) = secdq(:,:,:,:,mpirank) + secdq_t(:,:,:,:,mpirank)
        if( E1E2 ) secdq_m(:,:,:,:,mpirank) = secdq_m(:,:,:,:,mpirank) + secdq_m_t(:,:,:,:,mpirank)
        if( E2E2 ) secqq(:,:,:,:,:,mpirank) = secqq(:,:,:,:,:,mpirank) + secqq_t(:,:,:,:,:,mpirank)
        if( E2E2 ) secqq_m(:,:,:,:,:,mpirank) = secqq_m(:,:,:,:,:,mpirank) + secqq_m_t(:,:,:,:,:,mpirank)
        if( E1E3 ) secdo(:,:,:,:,:,mpirank) = secdo(:,:,:,:,:,mpirank) + secdo_t(:,:,:,:,:,mpirank)
        if( E1E3 ) secdo_m(:,:,:,:,:,mpirank) = secdo_m(:,:,:,:,:,mpirank) + secdo_m_t(:,:,:,:,:,mpirank)
        if( E3E3 ) secoo(:,:,:,:,:,mpirank) = secoo(:,:,:,:,:,mpirank) + secoo_t(:,:,:,:,:,mpirank)
        if( E3E3 ) secoo_m(:,:,:,:,:,mpirank) = secoo_m(:,:,:,:,:,mpirank) + secoo_m_t(:,:,:,:,:,mpirank)

        if( icheck(21) > 1 ) call Write_ten_opt(E1E1,E1E2,E2E2,E1M1,En,Energ(ie),M1M1,mpinodes,mpirank,n_rel,ninitl_out, &
                                                                   secdd,secdq,secmd,secmm,secqq)

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(7) = real(time,db)
          Time_rout(16) = Time_rout(16) + Time_loc(5) - Time_loc(3) ! Kern
          Time_rout(17) = Time_rout(17) + Time_loc(6) - Time_loc(5) ! Chi
          Time_rout(11) = Time_rout(11) + Time_loc(7) - Time_loc(6) ! Tenseur
        endif

      end do ! End of loop over occupied state electron energy

    endif ! arriving for parallel computing when ie > nenerg

    if( mpinodes > 1 ) then
      call MPI_RECV_all(l0_nrixs,lmax_nrixs,mpinodes,mpirank,mpirank0,Multipole,n_oo,n_rel,ninitl_out, &
                        nq_nrixs,S_nrixs,S_nrixs_l,secdd,secdo,secdq,secmd,secmm,secoo,secqq)
      if( Green_int ) call MPI_RECV_all(l0_nrixs,lmax_nrixs,mpinodes,mpirank,mpirank0,Multipole,n_oo, &
                        n_rel,ninitl_out,nq_nrixs,S_nrixs_m,S_nrixs_l_m,secdd_m,secdo_m,secdq_m,secmd_m,secmm_m,secoo_m,secqq_m)
    endif

    if( mpirank0 /= 0 ) cycle

    do ie_computer = 0,mpinodes-1

      ie = ( je - 1 ) * mpinodes + ie_computer + 1

      if( ie > nenerg ) exit

      call write_coabs(Allsite,angxyz,axyz,.false.,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
                  E_cut,Energ,Energphot,.false.,Epsii,Eseuil,Final_tddft,First_E, &
                  f_avantseuil,Full_self_abs,Green_int,hkl_dafs,iabsorig,icheck(21),ie,ie_computer,Int_tens, &
                  isigpi,isymeq,jseuil,Length_abs,ltypcal,Matper,Moyenne,mpinodes,Multipole,n_multi_run,n_oo,n_rel,n_tens_max, &
                  natomsym,nbseuil, &
                  ncolm,ncolr,ncolt,nenerg,ninit1,ninitl_out,nomabs,nomfich,nomfich_cal_convt,nomfich_s,nphi_dafs, &
                  npldafs,nphim,nplr,nplrm,nseuil,nspinp,numat,nxanout,pdp,phdafs,phdf0t, &
                  phdt,pol,poldafse,poldafss,sec_atom,secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, &
                  secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Self_abs,Spherical_signal, &
                  Spherical_tensor,Spinorbite,Taux_eq,V0muf,Vecdafse,Vecdafss,Vec,Volume_maille,Xan_atom)

      if( First_E ) nomfich_cal_tddft_conv = nomfich_cal_convt
      First_E = .false.

    end do

    if( mpirank0 == 0 ) then
      call CPU_TIME(time)
      Time_loc(8) = real(time,db)
      Time_rout(13) = Time_rout(13) + Time_loc(8) - Time_loc(7) ! Coabs
    endif

  end do boucle_energ   ! End of loop over photon energy

  deallocate( B_stk, Chi, Chi_0, dChi_0, dChi_0_i, Kern )
  deallocate( Ecinetic, Eimag_t, En, Enervide )
  deallocate( Ecinetic_e, Energ )
  deallocate( roff_rr, V0bdc, V_hubb_t, Vrato, zet )

  return
  100 format(/1x,120('-')//,' Cycle Tddft')
end

!***********************************************************************

subroutine Write_ten_opt(E1E1,E1E2,E2E2,E1M1,En,Energ,M1M1,mpinodes,mpirank,n_rel,ninitl_out, &
                                                                   secdd,secdq,secmd,secmm,secqq)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: i, j, k, mpinodes, mpirank, n_rel, ninitl_out

  complex(kind=db), dimension(3,3,n_rel,ninitl_out,0:mpinodes-1):: secdd
  complex(kind=db), dimension(3,3,ninitl_out,0:mpinodes-1):: secmd, secmm
  complex(kind=db), dimension(3,3,3,ninitl_out,0:mpinodes-1):: secdq
  complex(kind=db), dimension(3,3,3,3,ninitl_out,0:mpinodes-1):: secqq

  logical:: E1E1, E1E2, E2E2, E1M1, M1M1

  real(kind=db):: Energ
  real(kind=db), dimension(2):: En

  write(3,110) Energ*rydb, En(:)*rydb, mpirank

  if( E1E1 ) then
    write(3,'(/A)') ' E1E1'
    do i = 1,3
      write(3,120) secdd(i,:,1,1,mpirank)
    end do
  endif

  if( E1E2 ) then
    write(3,'(/A)') ' E1E2'
    do i = 1,3
      write(3,120) ( secdq(i,j,:,1,mpirank), j = 1,3 )
    end do
  endif

  if( E2E2 ) then
    write(3,'(/A)') ' E2E2'
    do i = 1,3
      write(3,120) ( ( secqq(i,j,k,:,1,mpirank), k = 1,3), j = 1,3 )
    end do
  endif

  if( E1M1 ) then
    write(3,'(/A)') ' E1M1'
    do i = 1,3
      write(3,120) secmd(i,:,1,mpirank)
    end do
  endif

  if( M1M1 ) then
    write(3,'(/A)') ' M1M1'
    do i = 1,3
      write(3,120) secmm(i,:,1,mpirank)
    end do
  endif

  return
  110 format(/ 50('-'), //' Photon energy =',f10.5,', Electron energies =',2f10.5,', mpirank =',i2)
  120 format(1p,3(2x,3(2x,3(2x,2e18.10))))
end

!***********************************************************************

! Calculation of Chi_0 for optics
! Chi_0 matrix is not symmetric because Taull_dd is using complex spherical harmonics and are not symmetric (and complex)

subroutine Chi_0_opt(Chi_0, dChi_0, dChi_0_i, delta_E, E_cut_tddft, Energ, Energ_s, Gamma_tddft, icheck, ie_g, ief, &
       je, lmax, M_depend, nenerg_s, nlm, nlmamax, nlms_f, nlms_g ,nomfich, &
       nspino, nspinp, nspint, roff0, roff0_ph, Spinorbite, Taull_tdd, Ylm_comp)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: icheck, ie, ie_g, ief, is1, is2, iso, iso_f1, iso_f2, iso_g1, iso_g2, isp, isp_f1, isp_f12, isp_f2, isp_g1, &
    isp_g12, isp_g2, ispm_f1, ispm_f2, ispm_g1, ispm_g2, iss, iss_f1, iss_f2, iss_g1, iss_g2, je, l, &
    lm, lm_f1, lm_f2, lm_g1, lm_g2, lmax, lms, lms_g1, lms_g2, lms_f1, lms_f2, lmv_f1, lmv_g1, lmv_f2, lmv_g2, m,&
    mv, nenerg_s, nlm, nlmamax, nlms, nlms_f, nlms_g, nspinp, nspino, nspint

  integer, dimension(nlms_f):: i_val, iss_val, l_val, lm_val, lmv_val, m_val, s_val

  logical:: Gamma_tddft, M_depend, Spinorbite, Ylm_comp

  character(len=132):: nomfich

  complex(kind=db):: dch, dch1
  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nlmamax,nspinp) :: Taull_tdd
  complex(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g):: Chi_0
  complex(kind=db), dimension(ief-1,nlms_f,nlms_f,nlms_g,nlms_g):: dChi_0, dChi_0_i

  real(kind=db):: Chi_0_r, Chi_0_r_im, Chi_0_i, Chi_0_i_im, d_i, d_r, dch_i, dch_r, delta_E, E_cut_tddft, Energ, &
    Gamma, p, param, pasmin, t1, t2

  real(kind=db), dimension(nenerg_s):: e1, e2, Energ_s
  real(kind=db), dimension(nlm,nlm,nspinp,nspinp,nspino,nspino):: roff0_ph
  real(kind=db), dimension(ief-1,ief:nenerg_s,nlm,nlm,nspinp,nspinp,nspino**2):: roff0

  Gamma = 0._db

  if( .not. Gamma_tddft )  then
    pasmin = 1.0_db / rydb
! le pas minimum
    do ie = 1,nenerg_s-1
      pasmin = min( pasmin, Energ_s(ie+1)-Energ_s(ie) )
    end do
    param = 0.5_db * pasmin
  end if

! Creation des tableaux bornes des intervales:
  do ie = 2,nenerg_s-1
    e1(ie) =  0.5_db * ( Energ_s(ie-1) + Energ_s(ie) )
    e2(ie) =  0.5_db * ( Energ_s(ie) + Energ_s(ie+1) )
  end do

! Bornes de l'intervale:
  e1(1) = Energ_s(1) - 0.5_db * ( Energ_s(2) - Energ_s(1) )
  e2(1) = 0.5_db * ( Energ_s(1) + Energ_s(2) )
  e1(nenerg_s) =  0.5_db * ( Energ_s(nenerg_s-1) + Energ_s(nenerg_s) )
  e2(nenerg_s) = Energ_s(nenerg_s) + 0.5_db * ( Energ_s(nenerg_s-1) - Energ_s(nenerg_s) )

  e1(ief) = E_cut_tddft

  lms = 0
  do isp = 1,nspint
    lm = 0
    do l = 0,lmax
      do m = -l,l
        lm = lm + 1
        do iso = 1,nspino
          if( Spinorbite ) then
            mv = m + isp - iso
            if( abs( mv ) > l ) cycle
            iss = min( iso, nspinp )
          else
            mv = m
            iss = min( isp, nspinp )
          endif
          lms = lms + 1
          s_val(lms) = isp
          l_val(lms) = l
          m_val(lms) = m
          i_val(lms) = iso
          iss_val(lms) = iss
          lm_val(lms) = lm
          if( M_depend ) then
            lmv_val(lms) = l**2 + l + 1 + mv
          else
            lmv_val(lms) = l + 1
          endif
        end do
      end do
    end do
  end do

  nlms = lms

  do lms_g1 = 1,nlms
    isp_g1 = s_val(lms_g1)
    ispm_g1 = min( isp_g1, nspinp )
    lm_g1 = lm_val(lms_g1)
    iso_g1 = i_val(lms_g1)
    iss_g1 = iss_val(lms_g1)
    lmv_g1 = lmv_val(lms_g1)

    do lms_g2 = 1,nlms
      isp_g2 = s_val(lms_g2)
      ispm_g2 = min( isp_g2, nspinp )
      lm_g2 = lm_val(lms_g2)
      iso_g2 = i_val(lms_g2)
      iss_g2 = iss_val(lms_g2)
      lmv_g2 = lmv_val(lms_g2)

      if( .not. Spinorbite .and. isp_g1 /= isp_g2 ) cycle

      isp_g12 = 2 * isp_g1 - 2 + isp_g2

      do lms_f1 = 1,nlms
        isp_f1 = s_val(lms_f1)
        ispm_f1 = min( isp_f1, nspinp )
        lm_f1 = lm_val(lms_f1)
        iso_f1 = i_val(lms_f1)
        iss_f1 = iss_val(lms_f1)
        lmv_f1 = lmv_val(lms_f1)

        is1 = ( iso_g1 - 1 )*nspino + iso_f1

        do lms_f2 = 1,nlms
          isp_f2 = s_val(lms_f2)
          ispm_f2 = min( isp_f2, nspinp )
          lm_f2 = lm_val(lms_f2)
          iso_f2 = i_val(lms_f2)
          iss_f2 = iss_val(lms_f2)
          lmv_f2 = lmv_val(lms_f2)

          if( .not. Spinorbite .and. isp_f1 /= isp_f2 ) cycle

          isp_f12 = 2 * isp_f1 - 2 + isp_f2

          is2 = ( iso_g2 - 1 )*nspino + iso_f2

          Chi_0_r = 0._db ; Chi_0_i = 0._db
          Chi_0_r_im = 0._db ; Chi_0_i_im = 0._db

          do ie = ief,nenerg_s

            t1 = Energ - e1(ie) + Energ_s(ie_g)
            t2 = Energ - e2(ie) + Energ_s(ie_g)
! < 1g 1 > < 2 2g >
            dch = Taull_tdd(ie,lm_f1,iss_f1,lm_f2,iss_f2) * Taull_tdd(ie_g,lm_g2,iss_g2,lm_g1,iss_g1) &
                * roff0(ie_g,ie,lmv_g1,lmv_f1,ispm_g1,ispm_f1,is1) * roff0(ie_g,ie,lmv_g2,lmv_f2,ispm_g2,ispm_f2,is2) * Delta_E

            dch = dch / ( roff0_ph(lmv_g1,lmv_f1,ispm_g1,ispm_f1,iso_g1,iso_f1) &
                        * roff0_ph(lmv_g2,lmv_f2,ispm_g2,ispm_f2,iso_g2,iso_f2) )

            dch_r = real( dch, db )
            dch_i = aimag( dch )

            if( Gamma_tddft ) then

              t1 = t1 / Gamma
              t2 = t2 / Gamma
              d_i = - ( atan(t2) - atan(t1) ) / pi
              d_r = - 0.5_db * log( (1+t2**2)/(1+t1**2) ) / pi
              if( Ylm_comp ) then
                Chi_0_r = Chi_0_r + dch_r * d_r
                Chi_0_i = Chi_0_i + dch_r * d_i
                Chi_0_r_im = Chi_0_r_im - dch_i * d_i
                Chi_0_i_im = Chi_0_i_im + dch_i * d_r
              else
                Chi_0_r = Chi_0_r + dch_r * d_r - dch_i * d_i
                Chi_0_i = Chi_0_i + dch_r * d_i + dch_i * d_r
              endif

            else    ! gamma = 0

! Avoid divergences when ranges overlap
              if( abs(t1) < param ) t1 =  param
              if( abs(t2) < param ) t2 = -param
! One eliminates divergence in 1/e of f' (t1 -> -t1;
! when |t1| = |t2| contributions from both sides of the point subsstract and become zero
! abs changes the signe when t1 and t2 are opposite in sign
              if( abs(t1+t2) > eps10 ) then
                Chi_0_r = Chi_0_r - dch_r * log( abs(t2/t1) ) / pi
                if( Ylm_comp ) then
                  Chi_0_i_im = Chi_0_i_im - dch_i * log( abs(t2/t1) ) / pi
                else
                  Chi_0_i = Chi_0_i - dch_i * log( abs(t2/t1) ) / pi
                endif
              endif

              if( abs( Energ - Energ_s(ie) + Energ_s(ie_g) ) < eps10 ) then
                Chi_0_i = Chi_0_i + dch_r
                if( Ylm_comp ) then
                  Chi_0_r_im = Chi_0_r_im + dch_i
                else
                  Chi_0_r = Chi_0_r + dch_i
                endif
              elseif( Energ > Energ_s(ie) - Energ_s(ie_g) .and. Energ < Energ_s(ie-1) - Energ_s(ie_g) ) then
                if( ie == ief ) then
                  dch1 = dch
                else
                  dch1 = Taull_tdd(ie-1,lm_f1,iss_f1,lm_f2,iss_f2) * Taull_tdd(ie_g,lm_g2,iss_g2,lm_g1,iss_g1) &
                       * roff0(ie_g,ie-1,lmv_g1,lmv_f1,ispm_g1,ispm_f1,is1) &
                       * roff0(ie_g,ie-1,lmv_g2,lmv_f2,ispm_g2,ispm_f2,is2)
                  dch1 = dch1 / ( roff0_ph(lmv_g1,lmv_f1,ispm_g1,ispm_f1,iso_g1,iso_f1) &
                                * roff0_ph(lmv_g2,lmv_f2,ispm_g2,ispm_f2,iso_g2,iso_f2) )
                endif
                p = ( Energ - Energ_s(ie-1) ) / ( Energ_s(ie) - Energ_s(ie-1) )
                dch = ( 1 - p ) * dch1 + p * dch
                dch_r = real( dch, db )
                dch_i = aimag( dch )
                Chi_0_i = Chi_0_i + dch_r
                if( Ylm_comp ) then
                  Chi_0_r_im = Chi_0_r_im + dch_i
                else
                  Chi_0_r = Chi_0_r + dch_i
                endif
              endif

            end if

          end do ! fin boucle energie etats inoccupes

          dChi_0(ie_g,lms_f1,lms_f2,lms_g1,lms_g2) = cmplx( Chi_0_r, Chi_0_i, db )
          if( Ylm_comp ) dChi_0_i(ie_g,lms_f1,lms_f2,lms_g1,lms_g2) = cmplx( Chi_0_r_im, Chi_0_i_im, db )

        end do ! end of loop over lms_f2
      end do ! end of loop over lms_f1
    end do
  end do

  Chi_0(:,:,:,:) = Chi_0(:,:,:,:) + dChi_0(ie_g,:,:,:,:)

  if( ie_g < ief-1 .or. icheck == 0 ) return

  if( icheck > 1 ) call Write_Chi_0_opt(Chi_0,Energ,Energ_s,i_val,icheck,je,l_val,lmax,m_val,nlms_f,nenerg_s,nlmamax,nlms_g, &
                           nomfich,nspinp,s_val,Taull_tdd)

  return
end

!**********************************************************************

subroutine Write_Chi_0_opt(Chi_0,Energ,Energ_s,i_val,icheck,je,l_val,lmax,m_val,nlms_f,nenerg_s,nlmamax,nlms_g, &
                           nomfich,nspinp,s_val,Taull_tdd)

  use declarations
  implicit none

  character(len=132):: chi_conv, nomfich
  character(len=1), dimension(nlms_f):: sgn

  integer:: icheck, ie, isp1, isp2, je, l, l1, l2, lm1, lm2, lmax, lms_f1, lms_f2, lms_g1, lms_g2, m1, m2, nenerg_s, &
            nlmamax, nlms_f, nlms_g, nspinp
  integer, dimension(nlms_f):: i_val, l_val, m_val, s_val

  complex(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g):: Chi_0
  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nlmamax,nspinp) :: Taull_tdd

  logical:: Diag

  real(kind=db):: Energ
  real(kind=db), dimension(nenerg_s):: Energ_s

  write(3,100)
  write(3,110) Energ*rydb, (( l_val(lms_f2), m_val(lms_f2), i_val(lms_f2), s_val(lms_f2), &
                  l_val(lms_g2), m_val(lms_g2), i_val(lms_g2), s_val(lms_g2), lms_f2 = 1,nlms_f), lms_g2 = 1,nlms_g)
  do lms_g1 = 1,nlms_g
    do lms_f1 = 1,nlms_f
      write(3,120) l_val(lms_g1), m_val(lms_g1), i_val(lms_g1), s_val(lms_g1), &
                   l_val(lms_f1), m_val(lms_f1), i_val(lms_f1), s_val(lms_f1), &
                   (( Chi_0(lms_f1,lms_f2,lms_g1,lms_g2), lms_f2 = 1,nlms_f), lms_g2 = 1,nlms_g)
    end do
  end do

  if( icheck > 2 .and. je == 1 ) then
    write(3,130) (((((( l1, m1, isp1, l2, m2, isp2, isp2 = 1,nspinp),  m2 = -l2,l2), l2 = 0,lmax ), isp1 = 1,nspinp), &
                                                                                         m1 = -l1,l1), l1 = 0,lmax )
    do ie = 1,nenerg_s
      write(3,140) Energ_s(ie)*rydb, (((( Taull_tdd(ie,lm1,isp1,lm2,isp2), isp2 = 1,nspinp), lm2 = 1,(1+lmax)**2 ), &
                                          isp1 = 1,nspinp), lm1 = 1,(1+lmax)**2 )
    end do
  end if

!------------------------------
  Diag = icheck < 3

! Nom du fichier de sortie
  chi_conv = nomfich
  l = len_trim( chi_conv )

  chi_conv(l+1:l+10) = '_Chi_0.txt'

  if( je == 1 ) then
    do lms_f1 = 1,nlms_f
      if( m_val(lms_f1) < 0 ) then
        sgn(lms_f1) = '-'
      else
        sgn(lms_f1) = '+'
      endif
    end do
    open(31, file = chi_conv)
    if( Diag ) then
      write(31,150) (( l_val(lms_g1), sgn(lms_g1), abs(m_val(lms_g1)), i_val(lms_g1), s_val(lms_g1), &
                       l_val(lms_f1), sgn(lms_f1), abs(m_val(lms_f1)), i_val(lms_f1), s_val(lms_f1), &
                       l_val(lms_f1), sgn(lms_f1), abs(m_val(lms_f1)), i_val(lms_f1), s_val(lms_f1), &
                       l_val(lms_g1), sgn(lms_g1), abs(m_val(lms_g1)), i_val(lms_g1), s_val(lms_g1), &
                       lms_f1 = 1,nlms_g ), lms_g1 = 1,nlms_f )
    else
      write(31,150) (((( l_val(lms_g1), sgn(lms_g1), abs(m_val(lms_g1)), i_val(lms_g1), s_val(lms_g1), &
                         l_val(lms_f1), sgn(lms_f1), abs(m_val(lms_f1)), i_val(lms_f1), s_val(lms_f1), &
                         l_val(lms_f2), sgn(lms_f2), abs(m_val(lms_f2)), i_val(lms_f2), s_val(lms_f2), &
                         l_val(lms_g2), sgn(lms_g2), abs(m_val(lms_g2)), i_val(lms_g2), s_val(lms_g2), &
                         lms_f1 = 1,nlms_f ), lms_f2 = 1,nlms_f ), lms_g1 = 1,nlms_g ), lms_g2 = 1,nlms_g )

    endif
  else
    open(31, file = chi_conv, position='append')
  endif

  if( Diag ) then
    write(31,160) Energ*rydb, (( Chi_0(lms_f1,lms_f1,lms_g1,lms_g1), lms_g1 = 1,nlms_g ), lms_f1 = 1,nlms_f )
  else
    write(31,160) Energ*rydb, (((( Chi_0(lms_f1,lms_f2,lms_g1,lms_g2), lms_g1 = 1,nlms_g), lms_g2 = 1,nlms_g), &
                                   lms_f1 = 1,nlms_f), lms_f2 = 1,nlms_f )
  endif
  Close(31)

  return
  100 format(/' ---- Chi_0_opt -------',100('-'))
  110 format(/' Chi_0, Energy =',f10.3,' eV',/ &
              32x,' l2f m2f i2f s2f l2g m2g i2g s2g',/ &
              ' l1g m1g i1g s1g l1f m1f i1f s1f',3x,10000(3x,'(',4i2,',',4i2,')',3x) )
  120 format(8i4,2x,1p,100000(1x,2e12.5))
  130 format(/' Taull_tdd(l1,m1,isp1,l2,m2,isp2)',/ &
              ' Energy', 250(5x,'(',3i2,',',3i2,')  ') )
  140 format(f7.3,1p,1000e11.3)
  150 format('   Energy ',100000('(',4(i1,a1,3i1),')'))
  160 format(f10.3,1p,200000e11.3)
end

!***********************************************************************

! Calculation of the kernel for optic in TDDFT

! The expansion in spherical harmonics of Coulombian potential is from J.Phys.A:Math. Gen. 39(2006) 8613-8630
! One calculates it using complex or real spherical harmonics

subroutine kernel_optic(Atomic_scr,Energ,First_E,fxc,icheck,Kern,Kern_fac,lmax,M_depend,n_Ec,nlm,nlm_fp,nlms_f,nlms_g, &
                nomfich,nr,nspino,nspint,nspinp,r,Rmtsd,RPALF,Spinorbite,Ylm_comp,zet)

  use declarations
  implicit none

  integer:: icheck, iso, iso_f1, iso_f2, iso_g1, iso_g2, &
     isp, isp_f1, isp_f2, isp_g1, isp_g2, ispv_f1, &
     ispv_f2, ispv_g1, ispv_g2, l, l_f1, l_f1_t, l_f2, l_f2_t, l_g1, l_g1_t, l_g2_t, l_g2, l0, lcut, lmax, lmp, lmp_f1, lmp_f2, &
     lmp_g1, lmp_g2, lms, lms_f1, lms_f1_t, lms_f2, lms_f2_t, lms_g1, lms_g1_t, lms_g2, lms_g2_t, &
     lmv_f1, lmv_f2, lmv_g1, lmv_g2, lp, m, m_f1, m_f2, m_g1, m_g2, m0, mp, mv, mv_f1, mv_f1_t, mv_f2, mv_f2_t, &
     mv_g1,  mv_g1_t, mv_g2, mv_g2_t, n_Ec, nlm, nlm_fp, nlms, nlms_f, nlms_g, nr, nspino, nspinp, nspint

  integer, dimension(nlms_f):: i_val, l_val, lmp_val, lmv_val, m_val, mv_val, s_val

  character(len=132):: nomfich

  complex(kind=db):: Gaunt_1, Gaunt_12, Gaunt_2, Gaunt_xcc, Gaunt_xrc

  logical:: Atomic_scr, First_E, M_depend, RPALF, Spinorbite, Ylm_comp

  real(kind=db), intent(in):: Energ, Kern_fac, Rmtsd
  real(kind=db), dimension(nr), intent(in):: r
  real(kind=db), dimension(nr,2,2), intent(in):: fxc
  real(kind=db), dimension(nr,nlm,nlm_fp,nspinp,nspino,n_Ec), intent(in):: zet

  real(kind=db):: Angl, f_integr3, fac, Gaunt4Y, Gaunt4Y_r, intrad_r
  real(kind=db), dimension(nr):: f, r1, r2, t1, t2, zet_1, zet_2
  real(kind=db),dimension(nlms_f,nlms_f,nlms_g,nlms_g):: Kern

  Kern(:,:,:,:) = 0._db

  lms = 0
  do isp = 1,nspint
    do l = 0,lmax
      do m = -l,l
! The (lp, mp) loops are for non spherical potential (not used in practice)
        lmp = 0
        do lp = 0,lmax
          do mp = -lp,lp
            if( nlm_fp == 1 .and. ( lp /= l .or. mp /= m ) ) cycle
            lmp = lmp + 1
            do iso = 1,nspino
              if( Spinorbite .and. nlm_fp == 1 ) then
                mv = m + isp - iso
                if( abs(mv) > l ) cycle
              else
                mv = m
              endif
              lms = lms + 1
              s_val(lms) = isp
              l_val(lms) = l
              m_val(lms) = m
              i_val(lms) = iso
              mv_val(lms) = mv
              if( M_depend ) then
               lmv_val(lms) = l**2 + l + 1 + mv
              else
               lmv_val(lms) = l + 1
              endif
              lmp_val(lms) = lmp
            end do
          end do
        end do
      end do
    end do
  end do

  nlms = lms

! 1) Coulomb kernel

! Cut-off: see the selction rules for Gaunt coefficients
  lcut = 2 * lmax

! Loop over occupied states (g)

  do lms_g1 = 1,nlms_g

    isp_g1 = s_val(lms_g1)
    ispv_g1 = min( isp_g1, nspinp )
    l_g1 = l_val(lms_g1)
    m_g1 = m_val(lms_g1)
    iso_g1 = i_val(lms_g1)
    mv_g1 = mv_val(lms_g1)
    lmv_g1 = lmv_val(lms_g1)
    lmp_g1 = lmp_val(lms_g1)

! When radial functions do not depend on m, there is a recopy further on
    if( .not. M_depend .and. m_g1 /= 0 ) cycle

    do lms_g2 = 1,nlms_g

      isp_g2 = s_val(lms_g2)
      ispv_g2 = min( isp_g2, nspinp )
      l_g2 = l_val(lms_g2)
      m_g2 = m_val(lms_g2)
      iso_g2 = i_val(lms_g2)
      mv_g2 = mv_val(lms_g2)
      lmv_g2 = lmv_val(lms_g2)
      lmp_g2 = lmp_val(lms_g2)

! When radial functions do not depend on m, there is a recopy further on
      if( .not. M_depend .and. m_g2 /= 0 ) cycle

      do lms_f1 = 1,nlms_f

        isp_f1 = s_val(lms_f1)
        ispv_f1 = min( isp_f1, nspinp )
        l_f1 = l_val(lms_f1)
        m_f1 = m_val(lms_f1)
        iso_f1 = i_val(lms_f1)
        mv_f1 = mv_val(lms_f1)
        lmv_f1 = lmv_val(lms_f1)
        lmp_f1 = lmp_val(lms_f1)

        if( isp_f1 /= isp_g1 ) cycle

! When radial functions do not depend on m, there is a recopy further on
        if( .not. M_depend .and. m_f1 /= 0 ) cycle

        do lms_f2 = 1,nlms_f

          isp_f2 = s_val(lms_f2)
          ispv_f2 = min( isp_f2, nspinp )
          l_f2 = l_val(lms_f2)
          m_f2 = m_val(lms_f2)
          iso_f2 = i_val(lms_f2)
          mv_f2 = mv_val(lms_f2)
          lmv_f2 = lmv_val(lms_f2)
          lmp_f2 = lmp_val(lms_f2)

          if( isp_f2 /= isp_g2 ) cycle

! When radial functions do not depend on m, there is a recopy further on
          if( .not. M_depend .and. m_f2 /= 0 ) cycle

! Atomic screening (most often set as false)
          if( Atomic_scr .and. l_f1 /= l_g1 + 1 .and. l_f2 /= l_g2 + 1 ) cycle

          zet_1(1:nr) = zet(1:nr,lmv_g1,lmp_g1,ispv_g1,iso_g1,1) * zet(1:nr,lmv_f1,lmp_f1,ispv_f1,iso_f1,2)
          zet_2(1:nr) = zet(1:nr,lmv_g2,lmp_g2,ispv_g2,iso_g2,1) * zet(1:nr,lmv_f2,lmp_f2,ispv_f2,iso_f2,2)

! Index 0 is on the expansion of the coulombian potential
          do l0 = 0,lcut
            if( l0 < abs(l_f1 - l_g1) .or. l0 > l_f1 + l_g1 ) cycle
            if( l0 < abs(l_f2 - l_g2) .or. l0 > l_f2 + l_g2 ) cycle
            if( mod( l_f1 + l0 + l_g1, 2 ) /= 0 ) cycle
            if( mod( l_f2 + l0 + l_g2, 2 ) /= 0 ) cycle

            select case(l0)
              case(0)
                r1(1:nr) = r(1:nr)**2
                r2(1:nr) = r(1:nr)
              case(1)
                r1(1:nr) = r(1:nr)**3
                r2(1:nr) = 1._db
              case(2)
                r1(1:nr) = r(1:nr)**4
                r2(1:nr) = 1 / r(1:nr)
              case default
                r1(1:nr) = r(1:nr)**(l0+2)
                r2(1:nr) = r(1:nr)**(-l0+1)
            end select

            f(1:nr) = r1(1:nr) * zet_2(1:nr)

            call ffintegr2_r(t1,f,r,nr,1,Rmtsd)

            f(1:nr) = r2(1:nr) * zet_2(1:nr)

            call ffintegr2_r(t2,f,r,nr,-1,Rmtsd)

            f(1:nr) = zet_1(1:nr) * ( r2(1:nr) * t1(1:nr) + r1(1:nr) * t2(1:nr) )

            intrad_r = f_integr3(r,f,1,nr,Rmtsd)
! le 2 de 2 / (r-r')
            fac = 8 * pi * intrad_r / ( 2 * l0 + 1 )

! These loops are to make the recopy when needed
            do lms_g1_t = 1,nlms_g
              if( M_depend .and. lms_g1_t /= lms_g1 ) cycle
              l_g1_t = l_val(lms_g1_t)
              if( l_g1 /= l_g1_t .or. isp_g1 /= s_val(lms_g1_t) ) cycle
              mv_g1_t = mv_val(lms_g1_t)
              do lms_g2_t = 1,nlms_g
                if( M_depend .and. lms_g2_t /= lms_g2 ) cycle
                l_g2_t = l_val(lms_g2_t)
                if( l_g2 /= l_g2_t .or. isp_g2 /= s_val(lms_g2_t) ) cycle
                mv_g2_t = mv_val(lms_g2_t)
                do lms_f1_t = 1,nlms_f
                  if( M_depend .and. lms_f1_t /= lms_f1 ) cycle
                  l_f1_t = l_val(lms_f1_t)
                  if( l_f1 /= l_f1_t .or. isp_f1 /= s_val(lms_f1_t) ) cycle
                  mv_f1_t = mv_val(lms_f1_t)
                  do lms_f2_t = 1,nlms_f
                    if( M_depend .and. lms_f2_t /= lms_f2 ) cycle
                    l_f2_t = l_val(lms_f2_t)
                    if( l_f2_t /= l_f2 .or. isp_f2 /= s_val(lms_f2_t) ) cycle
                    mv_f2_t = mv_val(lms_f2_t)

                    Gaunt_12 = ( 0._db, 0._db )

                    do m0 = -l0,l0

                      if( Ylm_comp ) then
                        Gaunt_1 = Gaunt_xcc(l_f1_t,mv_f1_t,l0,m0,l_g1_t,mv_g1_t,Ylm_comp)
                        Gaunt_2 = Gaunt_xcc(l_f2_t,mv_f2_t,l0,m0,l_g2_t,mv_g2_t,Ylm_comp)
                      else
   ! Only Y(l0,m0) is complex
                        Gaunt_1 = Gaunt_xrc(l_f1_t,mv_f1_t,l_g1_t,mv_g1_t,l0,m0,Ylm_comp)
                        Gaunt_2 = Gaunt_xrc(l_f2_t,mv_f2_t,l_g2_t,mv_g2_t,l0,m0,Ylm_comp)
                      endif

                      Gaunt_12 = Gaunt_12 + Gaunt_1 * Gaunt_2

                    end do

                    Kern(lms_f1_t,lms_f2_t,lms_g1_t,lms_g2_t) = Kern(lms_f1_t,lms_f2_t,lms_g1_t,lms_g2_t) &
                                                              + real( Gaunt_12, db ) * fac

                  end do
                end do
              end do
            end do

          end do    ! end of loop over l0

        end do   ! end of loop over lms_g2
      end do     ! end of loop over lms_g1

    end do  !  end of loop over lms_f2
  end do    !  end of loop over lms_f1

! Adding of the exchange-correlation part in LDA:
  if( .not. RPALF ) then

    do lms_g1 = 1,nlms_g

      isp_g1 = s_val(lms_g1)
      ispv_g1 = min( isp_g1, nspinp )
      l_g1 = l_val(lms_g1)
      m_g1 = m_val(lms_g1)
      iso_g1 = i_val(lms_g1)
      mv_g1 = mv_val(lms_g1)
      lmv_g1 = lmv_val(lms_g1)
      lmp_g1 = lmp_val(lms_g1)

      if( .not. M_depend .and. m_g1 /= 0 ) cycle

      do lms_g2 = 1,nlms_g

        isp_g2 = s_val(lms_g2)
        ispv_g2 = min( isp_g2, nspinp )
        l_g2 = l_val(lms_g2)
        m_g2 = m_val(lms_g2)
        iso_g2 = i_val(lms_g2)
        mv_g2 = mv_val(lms_g2)
        lmv_g2 = lmv_val(lms_g2)
        lmp_g2 = lmp_val(lms_g2)

        if( .not. M_depend .and. m_g2 /= 0 ) cycle

        do lms_f1 = 1,nlms_f

          isp_f1 = s_val(lms_f1)
          ispv_f1 = min( isp_f1, nspinp )
          l_f1 = l_val(lms_f1)
          m_f1 = m_val(lms_f1)
          iso_f1 = i_val(lms_f1)
          mv_f1 = mv_val(lms_f1)
          lmv_f1 = lmv_val(lms_f1)
          lmp_f1 = lmp_val(lms_f1)

          if( isp_f1 /= isp_g1 ) cycle

          if( .not. M_depend .and. m_f1 /= 0 ) cycle

          do lms_f2 = 1,nlms_f

            isp_f2 = s_val(lms_f2)
            ispv_f2 = min( isp_f2, nspinp )
            l_f2 = l_val(lms_f2)
            m_f2 = m_val(lms_f2)
            iso_f2 = i_val(lms_f2)
            mv_f2 = mv_val(lms_f2)
            lmv_f2 = lmv_val(lms_f2)
            lmp_f2 = lmp_val(lms_f2)

            if( isp_f2 /= isp_g2 ) cycle

            if( .not. M_depend .and. m_f2 /= 0 ) cycle

! Atomic screening
            if( Atomic_scr .and. l_f1 /= l_g1 + 1 .and. l_f2 /= l_g2 + 1 ) cycle

            zet_1(1:nr) = zet(1:nr,lmv_g1,lmp_g1,ispv_g1,iso_g1,1) * zet(1:nr,lmv_f1,lmp_f1,ispv_f1,iso_f1,2)
            zet_2(1:nr) = zet(1:nr,lmv_g2,lmp_g2,ispv_g2,iso_g2,1) * zet(1:nr,lmv_f2,lmp_f2,ispv_f2,iso_f2,2)

! isp_f1 = isp_f2 and isp_g1 = isp_g2
            f(1:nr) = fxc(1:nr,isp_f1,isp_f2) * zet_1(1:nr) * zet_2(1:nr)

            fac = f_integr3(r,f,1,nr,Rmtsd)

! These loops are to make the recopy when needed
            do lms_g1_t = 1,nlms_g
              if( M_depend .and. lms_g1_t /= lms_g1 ) cycle
              l_g1_t = l_val(lms_g1_t)
              if( l_g1 /= l_g1_t .or. isp_g1 /= s_val(lms_g1_t) ) cycle
              mv_g1_t = mv_val(lms_g1_t)
              do lms_g2_t = 1,nlms_g
                if( M_depend .and. lms_g2_t /= lms_g2 ) cycle
                l_g2_t = l_val(lms_g2_t)
                if( l_g2 /= l_g2_t .or. isp_g2 /= s_val(lms_g2_t) ) cycle
                mv_g2_t = mv_val(lms_g2_t)
                do lms_f1_t = 1,nlms_f
                  if( M_depend .and. lms_f1_t /= lms_f1 ) cycle
                  l_f1_t = l_val(lms_f1_t)
                  if( l_f1 /= l_f1_t .or. isp_f1 /= s_val(lms_f1_t) ) cycle
                  mv_f1_t = mv_val(lms_f1_t)
                  do lms_f2_t = 1,nlms_f
                    if( M_depend .and. lms_f2_t /= lms_f2 ) cycle
                    l_f2_t = l_val(lms_f2_t)
                    if( l_f2_t /= l_f2 .or. isp_f2 /= s_val(lms_f2_t) ) cycle
                    mv_f2_t = mv_val(lms_f2_t)

! In Gaunt4Y, the 2 first Ylm are complex-conjugate.
                    if( Ylm_comp ) then
                      Angl = Gaunt4Y(l_g2_t,mv_g2_t,l_f1_t,mv_f1_t,l_f2_t,mv_f2_t,l_g1_t,mv_g1_t)
                    else
                      Angl = Gaunt4Y_r(l_g2_t,mv_g2_t,l_f1_t,mv_f1_t,l_f2_t,mv_f2_t,l_g1_t,mv_g1_t)
                    endif

                    if( Abs( Angl ) < eps10 ) cycle

                    Kern(lms_f1_t,lms_f2_t,lms_g1_t,lms_g2_t) = Kern(lms_f1_t,lms_f2_t,lms_g1_t,lms_g2_t) + Angl * fac

                  end do
                end do
              end do
            end do

          end do
        end do
      end do
    end do

  endif

! To correct numerical differences
  do lms_f1 = 1,nlms_f
  do lms_g1 = 1,nlms_g
  do lms_f2 = 1,nlms_f
  do lms_g2 = 1,nlms_g
    Kern(lms_f1,lms_f2,lms_g1,lms_g2) = 0.5_db * ( Kern(lms_f1,lms_f2,lms_g1,lms_g2) + Kern(lms_f2,lms_f1,lms_g2,lms_g1) )
    Kern(lms_f2,lms_f1,lms_g2,lms_g1) = Kern(lms_f1,lms_f2,lms_g1,lms_g2)
  end do
  end do
  end do
  end do

  if( abs( Kern_fac - 1._db ) > eps10 ) Kern(:,:,:,:) = Kern_fac * Kern(:,:,:,:)

  if( icheck > 1 ) call Write_Kern_opt(icheck,Kern,Energ,First_E,i_val,l_val,m_val,nlms_f,nlms_g,nomfich,s_val)

  return
end

!**********************************************************************

subroutine Write_Kern_opt(icheck,Kern,Energ,First_E,i_val,l_val,m_val,nlms_f,nlms_g,nomfich,s_val)

  use declarations
  implicit none

  character(len=132):: chi_conv, nomfich
  character(len=1), dimension(nlms_f):: sgn

  integer:: icheck, l, lms_g1, lms_g2, lms_f1, lms_f2, nlms_f, nlms_g
  integer, dimension(nlms_f):: i_val, l_val, m_val, s_val

  real(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g):: Kern

  logical:: Diag, First_E

  real(kind=db) :: Energ

  if( First_E .or. icheck > 2 ) then
    write(3,100)
    write(3,110) Energ*rydb, (( l_val(lms_f2), m_val(lms_f2), i_val(lms_f2), s_val(lms_f2), &
                    l_val(lms_g2), m_val(lms_g2), i_val(lms_g2), s_val(lms_g2), lms_f2 = 1,nlms_f), lms_g2 = 1,nlms_g)
    do lms_g1 = 1,nlms_g
      do lms_f1 = 1,nlms_f

        write(3,120) l_val(lms_f1), m_val(lms_f1), i_val(lms_f1), s_val(lms_f1), &
                     l_val(lms_g1), m_val(lms_g1), i_val(lms_g1), s_val(lms_g1), &
                     (( Kern(lms_f1,lms_f2,lms_g1,lms_g2), lms_f2 = 1,nlms_f), lms_g2 = 1,nlms_g)
      end do
    end do
  end if

!-------------------------------

  Diag = icheck < 3

! Output file name
  chi_conv = nomfich
  l = len_trim( chi_conv )

  chi_conv(l+1:l+9) = '_Kern.txt'

  if( First_E ) then
    do lms_f1 = 1,nlms_f
      if( m_val(lms_f1) < 0 ) then
        sgn(lms_f1) = '-'
      else
        sgn(lms_f1) = '+'
      endif
    end do
    open(31, file = chi_conv)
    if( Diag ) then
      write(31,150) (( l_val(lms_f1), sgn(lms_f1), abs(m_val(lms_f1)), i_val(lms_f1), s_val(lms_f1), &
                         l_val(lms_f1), sgn(lms_f1), abs(m_val(lms_f1)), i_val(lms_f1), s_val(lms_f1), &
                         l_val(lms_g1), sgn(lms_g1), abs(m_val(lms_g1)), i_val(lms_g1), s_val(lms_g1), &
                         l_val(lms_g1), sgn(lms_g1), abs(m_val(lms_g1)), i_val(lms_g1), s_val(lms_g1), &
                         lms_g1 = 1,nlms_g ), lms_f1 = 1,nlms_f )
    else
      write(31,150) (((( l_val(lms_f1), sgn(lms_f1), abs(m_val(lms_f1)), i_val(lms_f1), s_val(lms_f1), &
                         l_val(lms_f2), sgn(lms_f2), abs(m_val(lms_f2)), i_val(lms_f2), s_val(lms_f2), &
                         l_val(lms_g1), sgn(lms_g1), abs(m_val(lms_g1)), i_val(lms_g1), s_val(lms_g1), &
                         l_val(lms_g2), sgn(lms_g2), abs(m_val(lms_g2)), i_val(lms_g2), s_val(lms_g2), &
                         lms_g1 = 1,nlms_g ), lms_g2 = 1,nlms_g ), lms_f1 = 1,nlms_f ), lms_f2 = 1,nlms_f )

    endif
  else
    open(31, file = chi_conv, position='append')
  endif

  if( Diag ) then
    write(31,160) Energ*rydb, (( Kern(lms_f1,lms_f1,lms_g1,lms_g1), lms_g1 = 1,nlms_g ), lms_f1 = 1,nlms_f )
  else
    write(31,160) Energ*rydb, (((( Kern(lms_f1,lms_f2,lms_g1,lms_g2), lms_g1 = 1,nlms_g), lms_g2 = 1,nlms_g), lms_f1 = 1,nlms_f), &
                                   lms_f2 = 1,nlms_f )
  endif
  Close(31)

  return
  100 format(/' ---- Kernel -------',100('-'))
  110 format(/' Kernel, Energy =',f10.3,' eV',// &
              ' l1f m1f i1f s1f l1g m1g i1g s1g',3x,10000('(',i1,i2,2i1,','i1,i2,2i1,')') )
  120 format(8i4,2x,1p,100000e13.5)
  150 format('   Energy ',100000('(',4(i1,a1,3i1),')'))
  160 format(f10.3,1p,100000(5x,e13.5,4x))
end

!***********************************************************************

! Calculation of Chi with :  dChi = (1 - Chi0 * K )^(-1) * dChi0

! Chi_0 in input is in crystallo convention (f" > 0)
! Chi in output is in physical convention (compl. conj. of crystallo)

! When not magnetic only the ( isp_g1 = isp_f1, isp_g2 = isp_f2 ) of chi are not zero
! Then, when isp_g1 == isp_g2  Chi = ( 1/2 + 1/2 /( 1 - 2.Chi_0.K ) ) Chi_0
!       when isp_g1 /= isp_g2  Chi = ( - 1/2 + 1/2 /( 1 - 2.Chi_0.K ) ) Chi_0

subroutine Cal_Chi_opt(B_stk, Chi, Chi_0, dChi_0, dchi_0_i, Energ, First_E, First_K, icheck, ie_g, ief, &
              Kern, Kern_fast, lmax_probe, lmax, n_stk, nd3, ninitl, nlm_probe, nlm_fp, nlms_f, nlms_g, &
              nomfich, ns_dipmag, nspino, nspint, Spinorbite, Ylm_comp)

  use declarations
  implicit none

  character(len=132), intent(in):: nomfich

  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,nd3,nd3,2,2,ns_dipmag), intent(out):: Chi
  complex(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g):: Chi_0
  complex(kind=db), dimension(ief-1,nlms_f,nlms_f,nlms_g,nlms_g):: dChi_0, dChi_0_i
  complex(kind=db), dimension(n_stk,n_stk):: B_stk
  complex(kind=db), dimension(:), allocatable:: V
  complex(kind=db), dimension(:,:), allocatable:: B, B2, B2_i, B_i, C, dC, dC_i

  integer:: i, i1, i2, icheck, ie_g, ief, is_dipmag, iso, iso_f1, iso_f2, iso_g1, iso_g2, isp, &
    isp_f1, isp_f1_t, isp_f2, isp_f2_t, isp_g1, isp_g1_t, isp_g2, isp_g2_t, j, l, l_f1, l_f2, l_g1, l_g2, lm, lmax, &
    lmax_probe, lms, lms_f1, lms_f2, lms_g1, lms_g2, lmsv_f1, lmsv_f2, lmsv_g1, lmsv_g2, m, &
    m_f1, m_f2, m_g1, m_g2, n_stk, nd3, ndim, ninitl, nlm_fp, nlm_probe, nlms, nlms_f, nlms_g, ns_dipmag, nspino, nspint

  integer, dimension(nlms_f):: i_val, l_val, lmsv_val, m_val, s_val

  logical:: First_E, First_K, Full_Kern, Kern_fast, Spinorbite, Stop_job, Ylm_comp

  real(kind=db):: Energ
  real(kind=db), dimension(nlms_f,nlms_f,nlms_g,nlms_g), intent(in):: Kern
  real(kind=db), dimension(:,:), allocatable:: K

  if( icheck > 1 ) write(3,100)
  if( icheck > 1 ) write(3,110) Energ * rydb

  Full_Kern = .true.

  Chi(:,:,:,:,:,:,:) = ( 0._db, 0._db )
  Stop_job = .false.

  ndim = nlms_f * nlms_g

  lms = 0
  do isp = 1,nspint
    do l = 0,lmax
      do m = -l,l
        do iso = 1,nspino
          if( Spinorbite .and. nlm_fp == 1 ) then
            if( abs(m + isp - iso) > l ) cycle
          endif
          lms = lms + 1
          s_val(lms) = isp
          l_val(lms) = l
          m_val(lms) = m
          i_val(lms) = iso
          lmsv_val(lms) = nspino * ( l**2 + l + m ) + iso
        end do
      end do
    end do
  end do

  nlms = lms

  allocate( B(ndim,ndim), dC(ndim,ndim) )
  if( Ylm_comp ) allocate( dC_i(ndim,ndim) )
  B(:,:) = (0._db, 0._db)
  dC(:,:) = (0._db, 0._db)
  if( Ylm_comp ) dC_i(:,:) = (0._db, 0._db)

  if( .not. Kern_fast .or. First_K ) then
    allocate( C(ndim,ndim), K(ndim,ndim) )
    C(:,:) = (0._db, 0._db)
    K(:,:) = 0._db
  endif

  i1 = 0
  do lms_g1 = 1,nlms_g
    isp_g1 = s_val(lms_g1)
    do lms_f1 = 1,nlms
      i1 = i1 + 1
      isp_f1 = s_val(lms_f1)

!      l1 = l_val(lms_f1)
!      if( mod(l1,2) /= imparite .and. imparite /= 2 ) cycle

      i2 = 0
      do lms_g2 = 1,nlms_g
        isp_g2 = s_val(lms_g2)
        do lms_f2 = 1,nlms
          i2 = i2 + 1
          isp_f2 = s_val(lms_f2)

 !         l2 = l_val(lms_f2)
 !         if( mod(l2,2) /= imparite .and. imparite /= 2 ) cycle

! Conjugate physical convention
          if( Spinorbite .or. ( isp_g1 == isp_g2 .and. isp_f1 == isp_f2 ) ) then
            dC(i1,i2) = Conjg( dChi_0(ie_g,lms_f1,lms_f2,lms_g1,lms_g2) )
            if( Ylm_comp ) dC_i(i1,i2) = Conjg( dChi_0_i(ie_g,lms_f1,lms_f2,lms_g1,lms_g2) )
            if( .not. Kern_fast .or. First_K ) C(i1,i2) = Conjg( Chi_0(lms_f1,lms_f2,lms_g1,lms_g2) )
          endif

          if( Full_Kern .or. ( isp_f1 == isp_g1 .and. isp_f2 == isp_g2 ) ) then
            if( .not. Kern_fast .or. First_K ) K(i1,i2) = Kern(lms_f1,lms_f2,lms_g1,lms_g2)
          endif

        end do
      end do
    end do
  end do

  if( .not. Kern_fast .or. First_K ) then

    do ! point to come bach when pb (for a talkative writing)

      if( icheck > 2 .or. Stop_job ) then
        write(3,'(/A)') '  dChi_0'
        do lm = 1,ndim
          write(3,140) dC(lm,:)
        end do
        write(3,'(/A)') '  Kern'
        do lm = 1,ndim
          write(3,145) K(lm,:)
        end do
      endif

      B = - Matmul( C, K )

! When not magnetic one calculates the inverse of 1 - 2 Chi_0.K
      if( nspint == 1 ) then
        if( Full_Kern ) then
          B(:,:) = 4 * B(:,:)
        else
          B(:,:) = 2 * B(:,:)
        endif
      endif

      do i = 1,ndim
        B(i,i) = 1 + B(i,i) ! B = (1 - chi0*K)
      end do

      if( icheck > 2 .or. Stop_job ) then
        write(3,'(/A)') '  (1 - Chi_0.K)'
        do lm = 1,ndim
          write(3,140) B(lm,:)
        end do
      endif

      if( Stop_job ) stop

! B = (1 - chi0*K)**(-1)
! The product of 2 symmetrical matrices is not necessarily symmetric, thus B is not symmetric
      call invcomp(ndim,B,ndim,ndim,0,Stop_job)

      if( .not. Stop_job ) exit

    end do  ! Loop for writing when problem

    if( icheck > 2 ) then
      write(3,'(/A)') '  (1 - Chi_0.K)**(-1)'
      do lm = 1,ndim
        write(3,140) B(lm,:)
      end do
    endif

  endif

  if( nspint == 1 ) then
    if( Full_Kern ) then
      B(:,:) = 0.25_db * B(:,:)
      do i = 1,ndim
        B(i,i) = 0.75_db + B(i,i) ! B = (1 - chi0*K)
      end do
    else
      B(:,:) = 0.5_db * B(:,:)
      do i = 1,ndim
        B(i,i) = 0.5_db + B(i,i) ! B = (1 - chi0*K)
      end do
    endif
  endif

  if( Kern_fast ) then
    if( First_K ) then
      B_stk(:,:) = B(:,:)
    else
      B(:,:) = B_stk(:,:)
    endif
  endif
  First_K = .false.

! dChi = (1 - chi0*K)**(-1) * dChi0
! Copy on V is to save memory space
  allocate( V(ndim) )
  if( Ylm_comp ) allocate( B_i(ndim,ndim) )

  do i = 1,ndim
    V(:) = B(i,:)
    do j = 1,ndim
      B(i,j) = sum( V(:) * dC(:,j) )
      if( Ylm_comp ) B_i(i,j) = sum( V(:) * dC_i(:,j) )
    end do
  end do

  deallocate( V )

  if( nspint == 1 ) then
    allocate( B2(ndim, ndim) )
    if( Ylm_comp ) allocate( B2_i(ndim, ndim) )
    B2(:,:) = B(:,:) - dC(:,:)
    if( Ylm_comp ) B2_i(:,:) = B_i(:,:) - dC_i(:,:)
  endif

  if( icheck > 1 ) then
    write(3,'(/A)') '  Chi'
    do lm = 1,ndim
      write(3,140) B(lm,:)
    end do
  endif

  if( icheck > 1 ) write(3,150)

! Filling of Chi for usefull lm (up to nlm_probe).
  i1 = 0
  do lms_g1 = 1,nlms
    isp_g1 = s_val(lms_g1)
    l_g1 = l_val(lms_g1)
    m_g1 = m_val(lms_g1)
    iso_g1 = i_val(lms_g1)
!    if( mod(l_g1,2) /= imparite .and. imparite /= 2 ) cycle
    lmsv_g1 = lmsv_val(lms_g1)

    do lms_f1 = 1,nlms
      i1 = i1 + 1
      isp_f1 = s_val(lms_f1)
      l_f1 = l_val(lms_f1)
      m_f1 = m_val(lms_f1)
      iso_f1 = i_val(lms_f1)
!      if( mod(l_f1,2) /= imparite .and. imparite /= 2 ) cycle
      lmsv_f1 = lmsv_val(lms_f1)

      if( l_g1 > lmax_probe .or. l_f1 > lmax_probe ) cycle

      i2 = 0
      do lms_g2 = 1,nlms
        isp_g2 = s_val(lms_g2)
        l_g2 = l_val(lms_g2)
        m_g2 = m_val(lms_g2)
        iso_g2 = i_val(lms_g2)
!        if( mod(l_g2,2) /= imparite .and. imparite /= 2 ) cycle
        lmsv_g2 = lmsv_val(lms_g2)

!            if( (isp_f1 == isp_f2 .and. isp_g1 /= isp_g2) .or. (isp_f1 /= isp_f2 .and. isp_g1 == isp_g2) ) cycle

        do lms_f2 = 1,nlms
          i2 = i2 + 1
          isp_f2 = s_val(lms_f2)
          l_f2 = l_val(lms_f2)
          m_f2 = m_val(lms_f2)
          iso_f2 = i_val(lms_f2)
!          if( mod(l_f2,2) /= imparite .and. imparite /= 2 ) cycle
          lmsv_f2 = lmsv_val(lms_f2)

          if( l_g2 > lmax_probe .or. l_f2 > lmax_probe ) cycle

          do isp_g1_t = 1,2
            if( nspint == 2 .and. isp_g1_t /= isp_g1 ) cycle
            do isp_f1_t = 1,2
              if( nspint == 2 .and. isp_f1_t /= isp_f1 ) cycle
              do isp_g2_t = 1,2
                if( nspint == 2 .and. isp_g2_t /= isp_g2 ) cycle
                do isp_f2_t = 1,2
                  if( nspint == 2 .and. isp_f2_t /= isp_f2 ) cycle

                  is_dipmag = 2 * isp_g1_t - 2 + isp_g2_t

                  if( nspint == 1 .and. ( ( isp_g1_t == isp_f1_t .and. isp_g2_t == isp_f2_t .and. isp_g1_t /= isp_g2_t ) &
                                   .or. ( Full_Kern .and. ( isp_g1_t /= isp_g2_t .or. isp_f1_t /= isp_f2_t ) ) ) ) then
                    if( Ylm_comp ) then
                      Chi(lmsv_f1,lmsv_f2,lmsv_g1,lmsv_g2,isp_f1_t,isp_f2_t,is_dipmag) &
                                                            = aimag( B2(i1,i2) ) - img * real( B2_i(i1,i2), db )
                    else
                      Chi(lmsv_f1,lmsv_f2,lmsv_g1,lmsv_g2,isp_f1_t,isp_f2_t,is_dipmag) = cmplx( aimag( B2(i1,i2) ), 0._db )
                    endif
                  elseif( ( nspint == 1 .and. .not. Full_Kern ) .and. isp_g1_t /= isp_f1_t .and. isp_g2_t /= isp_f2_t &
                                      .and. isp_g1_t == isp_g2_t .and. isp_f1_t == isp_f2_t ) then
                    if( Ylm_comp ) then
                      Chi(lmsv_f1,lmsv_f2,lmsv_g1,lmsv_g2,isp_f1_t,isp_f2_t,is_dipmag) &
                                                            = aimag( dC(i1,i2) ) - img * real( dC_i(i1,i2), db )
                    else
                      Chi(lmsv_f1,lmsv_f2,lmsv_g1,lmsv_g2,isp_f1_t,isp_f2_t,is_dipmag) = cmplx( aimag( dC(i1,i2) ), 0._db )
                    endif
                  elseif( .not. nspint == 1 .or. ( isp_g1_t == isp_f1_t .and. isp_g2_t == isp_f2_t  .and. isp_g1_t == isp_g2_t ) &
                                            .or. Full_Kern ) then
                    if( Ylm_comp ) then
                      Chi(lmsv_f1,lmsv_f2,lmsv_g1,lmsv_g2,isp_f1_t,isp_f2_t,is_dipmag) &
                                                               = aimag( B(i1,i2) ) - img * real( B_i(i1,i2), db )
                    else
                      Chi(lmsv_f1,lmsv_f2,lmsv_g1,lmsv_g2,isp_f1_t,isp_f2_t,is_dipmag) = cmplx( aimag( B(i1,i2) ), 0._db )
                    endif

                  endif

                  if( icheck > 1 .and. abs( B(i1,i2) ) > eps10 ) write(3,160) l_g1, m_g1, iso_g1, isp_g1_t, &
                     l_f1, m_f1, iso_f1, isp_f1_t, l_f2, m_f2, iso_f2, isp_f2_t, l_g2, m_g2, iso_g2, isp_g2_t, &
                     Chi(lmsv_f1,lmsv_f2,lmsv_g1,lmsv_g2,isp_f1_t,isp_f2_t,is_dipmag)

                end do
              end do
            end do
          end do

        end do
      end do
    end do
  end do

  deallocate( B, dC )
  if( Ylm_comp ) deallocate( B_i, dC_i )
  if( nspint == 1 ) deallocate( B2 )
  if( nspint == 1 .and. Ylm_comp ) deallocate( B2_i )
  if( .not. Kern_fast .or. First_K ) deallocate( C, K )

  if( icheck > 1 ) call Write_Chi(Chi,Energ,First_E,ninitl,icheck,nlm_probe,nomfich,ns_dipmag,nspino)

  return
  100 format(/' ---- Cal_chi -------',100('-'))
  110 format(/' Energ =',f10.3,' eV')
  120 format(1p,9(1x,2e10.2))
  130 format(500(7x,'(',2(i2,','),i2,')',5x))
  135 format(500(1x,'(',2(i2,','),i2,')'))
  140 format(1p,1000(1x,2e10.3))
  145 format(1p,1000e11.3)
  150 format(/'  l_g1 m_g1 i_g1 s_g1 l_f1 m_f1 i_f1 s_f1 l_f2 m_f2 i_f2 s_f2 l_g2 m_g2 i_g2 s_g2              chi')
  160 format(16i5,3x,1p,2e13.5)

end
