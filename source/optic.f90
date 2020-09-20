! FDMNES subroutine

! Set of subroutines for procedure optic

!***********************************************************************

subroutine main_optic(Abs_U_iso,angxyz,Allsite,axyz,Bragg_abs,Cartesian_tensor,Classic_irreg, &
          Core_resolved,Dafs,Dafs_bio,dv0bdcF,E_cut,E_cut_imp,E_cut_man,Eclie,Eneg,Energ_s,Ephot_min, &
          Extract_ten,Eseuil,Full_potential,Full_self_abs,Green,hkl_dafs,Hubb_a,Hubb_d,icheck, &
          iabsorig,ip_max,ip0,isigpi,isymeq, &
          jseuil,ldip,Length_abs,lmax_pot,lmax_probe,lmaxabs_t,lmaxat0,lmaxfree,lmoins1,loct,lplus1,lqua, &
          lseuil,ltypcal,m_hubb,Matper,Moyenne,mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq, &
          msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,multi_0,Multipole,n_bulk_sup,n_multi_run, &
          n_oo,n_rel,n_rout,n_tens_max,natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ninit1,ninit,ninitr,nlm_pot, &
          nlm_probe,nlmamax,nomabs,nomfich,nomfich_cal_conv,nomfich_s,nphi_dafs,nphim,npldafs,nplr,nplrm, &
          nr,nrm,nseuil,nspin,nspino,nspinp, &
          numat_abs,nxanout,pdp,phdf0t,phdt,pol,poldafse,poldafss,psii, &
          r,Relativiste,Renorm,Rmtg,Rmtsd,rot_atom_abs,Rot_int,Self_abs,Solsing_only,Spherical_signal, &
          Spherical_tensor,Spinorbite,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Taull_tdd,Taux_eq,Tddft,Time_rout,V_intmax, &
          V_hubb,V0muf,Vcato,vec,vecdafse,vecdafss,Vhbdc,Volume_maille,VxcbdcF,Vxcato,Workf,Ylm_comp)

  use declarations
  implicit none

  integer, parameter:: n_bulk_z = 0
  integer, parameter:: n_bulk_z_abs = 0
  integer, parameter:: n_bulk_z_max_abs = 0
  integer, parameter:: n_max = 1

  integer:: i_range, iabsorig, icheck_s, ie, ie_computer, ie_q, ie_s, ie_t, ip_max, ip0, &
    iso1, iso2, isp, isp1, isp2, iss1, iss2, je, jseuil, l0_nrixs, lm1, lm2, lmax, lmax_pot, &
    lmax_probe, lmaxabs_t, lms1, lms2, lseuil, m_hubb, mpinodes, mpirank, mpirank0, &
    lmax_nrixs, lmaxat0, multi_0, n_abs_rgh, n_bulk_sup, n_comp, n_Ec, n_multi_run, n_oo, n_rel, n_rout, &
    n_tens_max, n_V, natomsym, nbseuil, ncolm, ncolr, ncolt, &
    ndim2, nenerg, nenerg_s, nenerg_tddft, nge, ninit1, ninit, ninitr, nlm_pot, nlm_probe, nlm_p_fp, &
    nlmamax, nphim, npldafs, nplr, nplrm, nq_nrixs, nr, nrm, ns_dipmag, nseuil, nspin, &
    nspino, nspinp, numat_abs, nxanout

  integer, dimension(30):: icheck
  integer, dimension(natomsym):: isymeq
  integer, dimension(ninit):: is_g
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(ninit,2):: m_g
  integer, dimension(2,npldafs):: isigpi
  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3)::  msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi
  integer, dimension(n_bulk_z_abs):: n_bulk_zc_abs
  integer, dimension(n_bulk_z_max_abs,n_bulk_z_abs):: igr_bulk_z_abs

  character(len=Length_name):: nomfich, nomfich_s
  character(len=5), dimension(nplrm):: ltypcal
  character(len=length_word), dimension(ncolm):: nomabs
  character(len=Length_name), dimension(n_multi_run+n_bulk_sup):: nomfich_cal_conv

  complex(kind=db):: f_avantseuil
  complex(kind=db), dimension(3,nplrm):: pol
  complex(kind=db), dimension(natomsym,npldafs):: Bragg_abs
  complex(kind=db), dimension(npldafs,nphim):: phdf0t
  complex(kind=db), dimension(npldafs,nphim,n_max):: phdt
  complex(kind=db), dimension(npldafs,n_max):: Sum_Bragg_bulk_abs_f0
  complex(kind=db), dimension(npldafs,n_bulk_z):: Sum_Bragg_nonabs_f
  complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss
  complex(kind=db), dimension(3,3,n_rel,ninitr,0:mpinodes-1):: secdd, secdd_m, secdd_t, secdd_m_t
  complex(kind=db), dimension(3,3,ninitr,0:mpinodes-1):: secmd, secmd_m, secmm, secmm_m, secmd_t, secmd_m_t, secmm_t, secmm_m_t
  complex(kind=db), dimension(3,3,3,ninitr,0:mpinodes-1):: secdq, secdq_m, secdq_t, secdq_m_t
  complex(kind=db), dimension(3,3,3,3,ninitr,0:mpinodes-1):: secdo, &
    secdo_m, secqq, secqq_m, secdo_t, secdo_m_t, secqq_t, secqq_m_t
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:mpinodes-1):: secoo, secoo_m, secoo_t, secoo_m_t
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nenerg_s,nlmamax,nspinp,nlmamax,nspinp):: Taull_tdd
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: rof0, S_nrixs
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: Taull_abs

  logical:: Allsite, Cartesian_tensor, Classic_irreg, Core_resolved, Dafs, Dafs_bio, E_cut_man, &
    E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, Eneg, Energphot, Extract_ten, FDM_comp, Final_optic, Final_tddft, First_E, &
    Full_potential, Full_self_abs, Green, Green_int, Hubb_a, Hubb_d, lmaxfree, lmoins1, lplus1, M1M1, Matper, &
    Moyenne, Relativiste, Renorm, Self_abs, Solsing, &
    Solsing_only, Spherical_signal, Spherical_tensor, Spinorbite, Tddft, Xan_atom, Ylm_comp

  logical, dimension(10):: Multipole

  real(kind=sg):: time

  real(kind=db):: Abs_U_iso, delta_E, E_cut, E_cut_imp, E_cut_optic, Eclie, Ecmax, Ephot_min, &
    p, Rmtg, Rmtsd, Surface_ref, V_intmax, V0muf, Vhbdc, Volume_maille, Workf

  real(kind=db), dimension(3):: angxyz, axyz
  real(kind=db), dimension(3):: tp
  real(kind=db), dimension(n_rout):: Time_rout
  real(kind=db), dimension(nspin):: dv0bdcF, VxcbdcF
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitr):: Epsii, sec_atom
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(natomsym) :: Taux_eq
  real(kind=db), dimension(3,3):: rot_atom_abs, Rot_int
  real(kind=db), dimension(nplrm,2):: pdp
  real(kind=db), dimension(3,nplrm):: vec
  real(kind=db), dimension(3,npldafs):: hkl_dafs
  real(kind=db), dimension(npldafs):: Length_abs
  real(kind=db), dimension(n_bulk_z):: Length_rel
  real(kind=db), dimension(n_bulk_z_abs):: Length_rel_abs
  real(kind=db), dimension(ninit,2):: coef_g
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(n_tens_max*ninitr,0:natomsym):: Int_tens
  real(kind=db), dimension(nr,nlm_pot):: Vcato
  real(kind=db), dimension(3,npldafs,nphim):: vecdafse, vecdafss
  real(kind=db), dimension(nr,nlm_pot,nspin) :: Vxcato

  real(kind=db), dimension(:), allocatable:: Eimag, En, Energ, Enervide
  real(kind=db), dimension(:,:), allocatable:: Ecinetic, V0bdc
  real(kind=db), dimension(:,:,:,:), allocatable:: Vrato

  if( icheck(1) > 0 ) write(3,100)

  i_range = 1
  n_abs_rgh = 0
  Epsii(:) = 0._db
  Green_int = .false.
  FDM_comp = .false.
  n_V = 1
  n_Ec = 2
  ns_dipmag = 2  ! corresponds here at the 2, occupied and non-occupied state, energies, not the spin-flip in the E1M1 transition 
  ndim2 = 1
  Surface_ref = 0._db
  
  n_comp = 1
  
  if( Hubb_a .or. Full_potential ) then
    nlm_p_fp = nlm_probe
  else
    nlm_p_fp = 1
  endif

  allocate( Taull_abs(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag) )
  Taull_abs(:,:,:,:,:,:,:) = (0._db, 0._db)
  allocate( En(n_Ec) )
  allocate( V0bdc(nspin,n_V) )
  
  if( E_cut_man ) then
    E_cut_optic = E_cut_imp
  else
    E_cut_optic = E_cut
  endif

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  f_avantseuil = (0._db, 0._db)
  Energphot = .false.
  Final_optic = .true.
  Final_tddft = .false.
  sec_atom(:) = 0._db

  allocate( Ecinetic(nspin,n_Ec) )
  allocate( Enervide(n_Ec) )

! Elaboration of energy grid 
  call dim_grille_optic(E_cut_optic,Energ_s,icheck(29),mpirank0,nenerg,nenerg_s)

  allocate( Energ(nenerg) )
  allocate( Eimag(nenerg) )

  call grille_optic(E_cut_optic,Energ,Energ_s,icheck(29),nenerg,nenerg_s)

  Eimag(:) = 0._db
  Solsing = .false.

  Xan_atom = .false.

  allocate( Vrato(nr,nlm_pot,nspin,n_V) )

  do isp = 1,nspin
    Vrato(1:nr,:,isp,1) = Vcato(1:nr,:) + Vxcato(1:nr,:,isp)
    V0bdc(isp,1) = Vhbdc + VxcbdcF(isp) + dV0bdcF(isp)   
  end do

  nge = ( nenerg - 1 ) / mpinodes + 1

  First_E = .true.
  
! Loop over photon energy
  boucle_energ: do je = 1,nge

    ie = ( je - 1 ) * mpinodes + mpirank + 1
    
    if( ie <= nenerg ) then

      if( Energ(ie) < Ephot_min - eps10 ) cycle

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

! Loop over valence electron energy (below Fermi level)
      do ie_s = 1,nenerg_s

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          tp(1) = real(time,db)
        endif

        En(1) = Energ_s(ie_s)
        En(2) = Energ_s(ie_s) + Energ(ie)
        if( En(1) > E_cut_optic - eps10 ) exit
        if( En(2) < E_cut_optic - eps10 .or. En(2) > Energ_s(nenerg_s) + eps10 ) cycle

! n_Ec = 2
! ie_t = 1: index of valence level energy
! ie_t = 2: index of conduction level energy
        do ie_t = 1,n_Ec
          Enervide(ie_t) = En(ie_t) - Workf
          Ecinetic(:,ie_t) = Enervide(ie_t) - V0bdc(:,1)
          if( .not. Eneg ) Ecinetic(:,ie_t) = max( Ecinetic(:,ie_t), Eclie )
        end do

        Ecmax = maxval( Ecinetic )
        call clmax(Ecmax,Rmtg,lmaxat0,lmax,numat_abs,lmaxfree)
        lmax = min(lmax,lmaxabs_t)

! Calculation of cartesian tensor

        icheck_s = max( icheck(29), icheck(20) )

        do ie_t = 1,n_Ec
          do ie_q = 2,nenerg_s
            if( Energ_s(ie_q) > En(ie_t) - eps10 ) exit
          end do
          p = ( En(ie_t) - Energ_s(ie_q - 1) ) / ( Energ_s(ie_q) - Energ_s(ie_q - 1) )
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
                        if( isp1 /= isp2 ) cycle
                        iss1 = min(isp1,nspinp)
                        iss2 = min(isp2,nspinp)
                      endif
! In practice interpolation is only on the conduction electron (ie_t = 2 ) and when steps are different
                      Taull_abs(lms1,lms2,1,1,isp1,isp2,ie_t) = p * Taull_tdd(ie_q,lm1,iss1,lm2,iss2) &
                                                              + (1 - p ) * Taull_tdd(ie_q-1,lm1,iss1,lm2,iss2)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do

        if( ie_s == nenerg_s ) then
          delta_E = Energ_s(ie_s) - Energ_s(ie_s-1)
        elseif( ie_s == 1 ) then
          delta_E = Energ_s(ie_s+1) - Energ_s(ie_s)
        elseif(  Energ_s(ie_s+1) > E_cut_optic ) then
          delta_E = E_cut_optic - 0.5_db * ( Energ_s(ie_s) + Energ_s(ie_s-1) )
        else
          delta_E = 0.5_db * ( Energ_s(ie_s+1) - Energ_s(ie_s-1) )
        endif
! The square-root of the step energy is set to the square in subroutine tenseur_car
        Taull_abs(:,:,:,:,:,:,:) = sqrt( Delta_E ) * Taull_abs(:,:,:,:,:,:,:) 
 
        nenerg_tddft = 0
        allocate( rof0(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil) )

        call tenseur_car(Classic_irreg,coef_g,Core_resolved,Ecinetic, &
                Eimag(ie),Energ(ie),Enervide,Eseuil,FDM_comp,Final_optic,Final_tddft,Full_potential,Green,Green_int,Hubb_a, &
                Hubb_d,icheck_s,ie,ip_max,ip0,is_g,lmax_probe,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
                n_comp,n_Ec,n_oo,n_rel,n_V,nbseuil,ns_dipmag,ndim2,nenerg_tddft,ninit1,ninit,ninitr,ninitr,nlm_pot,nlm_probe, &
                nlm_p_fp,nlmamax,nr,nrm,nspin,nspino,nspinp,numat_abs,psii,r,Relativiste,Renorm, &  
                Rmtg,Rmtsd,rof0,rot_atom_abs, &
                Rot_int,secdd_t,secdd_m_t,secdo_t,secdo_m_t,secdq_t,secdq_m_t,secmd_t,secmd_m_t,secmm_t,secmm_m_t,secoo_t, &
                secoo_m_t,secqq_t,secqq_m_t,Solsing,Solsing_only,Spinorbite,Taull_abs,Tddft,V_hubb,V_intmax,V0bdc,Vrato,Ylm_comp)

        deallocate( rof0 )

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

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          tp(2) = real(time,db)
          Time_rout(11) = Time_rout(11) + tp(2) - tp(1)
        endif

      end do ! End of loop over electron energies

    endif  ! arrival when ie > nenerg

    if( mpinodes > 1 ) then
      l0_nrixs = 0; lmax_nrixs = 0; nq_nrixs = 0
      allocate( S_nrixs(nq_nrixs,(lmax_nrixs+1)**2,(lmax_nrixs+1)**2,ninitr,0:mpinodes-1) )
      call MPI_RECV_all(l0_nrixs,lmax_nrixs,mpinodes,mpirank,mpirank0,Multipole,n_oo,n_rel,ninitr, &
                       nq_nrixs,S_nrixs,secdd,secdo,secdq,secmd,secmm,secoo,secqq)
      deallocate( S_nrixs )
    endif

    if( mpirank0 /= 0 ) cycle
    call CPU_TIME(time)
    tp(2) = real(time,db)

    icheck_s = max( icheck(29), icheck(21) )

    do ie_computer = 0,mpinodes-1

      ie = ( je - 1 ) * mpinodes + ie_computer + 1

      if( ie > nenerg ) exit

      call Write_coabs(Abs_U_iso,Allsite,angxyz,axyz,Bragg_abs,.false.,Cartesian_tensor, &
             Core_resolved,Dafs,Dafs_bio,E_cut_optic,Energ,Energphot,Extract_ten,Epsii,Eseuil,Final_tddft,First_E, &
             f_avantseuil,Full_self_abs,Green_int,hkl_dafs,i_range,iabsorig,icheck_s,ie,ie_computer,igr_bulk_z_abs, &
             Int_tens,isigpi,isymeq,jseuil,Length_abs,Length_rel,Length_rel_abs,ltypcal,Matper,Moyenne,mpinodes, &
             multi_0,Multipole,n_abs_rgh,n_bulk_sup,n_multi_run, &
             n_bulk_z,n_bulk_z_abs,n_bulk_z_max_abs,n_bulk_zc_abs,n_max,n_oo,n_rel,n_tens_max,natomsym,nbseuil, &
             ncolm,ncolr,ncolt,nenerg,ninit1,ninitr,nomabs,nomfich,nomfich_cal_conv,nomfich_s,nphi_dafs,npldafs, &
             nphim,nplr,nplrm,nseuil,nspinp,numat_abs,nxanout,pdp,phdf0t,phdt,pol,poldafse,poldafss, &
             sec_atom,secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, &
             secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Self_abs,Spherical_signal, &
             Spherical_tensor,Spinorbite,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_ref,Taux_eq,V0muf,Vecdafse,Vecdafss, &
             Vec,Volume_maille,Xan_atom) 

      First_E = .false.
      
    end do

    call CPU_TIME(time)
    tp(3) = real(time,db)
    Time_rout(13) = Time_rout(13) + tp(3) - tp(2)

  end do boucle_energ

  deallocate( En, Taull_abs, V0bdc, Vrato )
  deallocate( Ecinetic, Eimag, Energ, Enervide )

  return
  100 format(/1x,120('-')//,' Cycle Optic')
end

!***********************************************************************

! Subroutine giving the dimension of the photon energy table in optics

subroutine dim_grille_optic(E_cut,Energ_s,icheck,mpirank0,nenerg,nenerg_s)

  use declarations
  implicit none

  integer:: icheck, mpirank0, nenerg, nenerg_s

  real(kind=db), intent(in):: E_cut
  real(kind=db), dimension(nenerg_s), intent(in):: Energ_s

  integer:: ie, ipr, ne

  if( E_cut < Energ_s(1) .and. mpirank0 == 0 ) then
    call write_error
    do ipr = 3,9,3
      if( ipr == 3 .and. icheck == 0 ) cycle
      write(ipr,110) E_cut*rydb, Energ_s(1)*rydb
    end do
    stop
  endif

  if( E_cut > Energ_s(nenerg_s) .and. mpirank0 == 0 ) then
    call write_error
    do ipr = 3,9,3
      if( ipr == 3 .and. icheck == 0 ) cycle
      write(ipr,120) E_cut*rydb, Energ_s(nenerg_s)*rydb
    end do
    stop
  endif

  do ie = 1,nenerg_s
    if( Energ_s(ie) > E_cut + eps10 ) exit
  end do
  ne = ie - 1

  nenerg = nenerg_s - ne

  return
  110 format(///' E_Fermi =',f7.3,' eV < First energy =',f7.3,' eV',// ' Increase the energy range of calculation !'//)
  120 format(///' E_Fermi =',f7.3,' eV > Last energy =',f7.3,' eV',// ' Increase the energy range of calculation !'//)
end

!***********************************************************************

! Subroutine giving the photon energy table for optics.

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
    Energ(ie) = Energ_s(ie+ne) - Energ_s(ne)
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

