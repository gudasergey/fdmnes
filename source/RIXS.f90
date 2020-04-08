! FDMNES subroutine

! Set of subroutines for procedure RIXS

!***********************************************************************

! Vhbdc = Potential interticiel de Hartree


subroutine main_RIXS(Ampl_abs,Classic_irreg,Core_resolved,dV0bdcF,E_cut,E_cut_imp,E_Fermi,E_cut_man,Eclie,Eneg,Energ_s, &
        Energphot,Epsii_a,Epsii_moy,Eseuil,FDM_comp_m,Full_potential,Hubb_a,Hubb_d, &
        icheck,ip0,ip_max,is_g,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb,msymdd,msymddi,msymdq,msymdqi, &
        msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole,n_oo,n_rel,nbseuil,nenerg_s,ninit1,ninitl,ninitlr,nlm_pot, &
        nlm_probe,nlmomax,nr,nrm,ns_dipmag,nspin,nspino,nspinp,numat,psii,r,Relativiste,Renorm,Rhoato,Rmtg,Rmtsd,Rot_atom_abs, &
        Rot_int,Spinorbite,V_intmax,V_hubb,V0muf,Vcato,Vhbdc,VxcbdcF,Vxcato,Ylm_comp)

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

  integer:: icheck, ie, ie_loss, ip0, ip_max, isp, je, lmax_pot, lmax_probe, lseuil, m_hubb, n_oo, n_rel, nbseuil, ne_loss, &
            nenerg, nenerg_s, ninit1, ninitl, ninitlr, ninitlv, nlm_p_fp, nlm_pot, nlm_probe, nlmomax, nr, nrm, ns_dipmag, nspin, &
            nspino, nspinp, numat

  integer, dimension(ninitl,2):: m_g
  integer, dimension(ninitl):: is_g
  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3):: msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi
  
  complex(kind=db), dimension(nenerg_s,nlm_probe,nspinp,nlmomax,nspinp):: Ampl_abs
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull
  complex(kind=db), dimension(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil):: rof0   ! empty
  complex(kind=db), dimension(3,3,ninitlr,0:mpinodes-1):: secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitlr,0:mpinodes-1):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:mpinodes-1):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,3,n_rel,ninitlr,0:mpinodes-1):: secdd, secdd_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:mpinodes-1):: secoo, secoo_m

  logical:: Classic_irreg, Core_resolved, E_cut_man, Eneg, Energphot, FDM_comp_m, Final_optic, Final_tddft, Full_potential, &
            Green, Green_int, Hubb_a, Hubb_d, lmoins1, lplus1, NRIXS, Relativiste, Renorm, RIXS, Solsing, Solsing_only, &
            Spinorbite, Tddft, Ylm_comp 

  logical, dimension(10):: Multipole

  real(kind=db):: E_cut, E_cut_imp, E_Fermi, E_cut_rixs, Eclie, Eimag, Epsii_moy, Rmtg, Rmtsd, V_intmax, Vhbdc, V0muf
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitlr):: Epsii_a
  real(kind=db), dimension(nr,nspin):: rhoato
  real(kind=db), dimension(ninitl,2):: coef_g
  real(kind=db), dimension(nspin):: dV0bdcF, VxcbdcF
  real(kind=db), dimension(nspin,n_V):: V0bdc
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(n_Ec):: Enervide
  real(kind=db), dimension(3,3):: Rot_atom_abs, Rot_int
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nspin,n_Ec):: Ecinetic
  real(kind=db), dimension(nr,nlm_pot):: Vcato
  real(kind=db), dimension(nr,nlm_pot,nspin,n_V):: Vrato
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vxcato

  real(kind=db), dimension(:), allocatable:: Decal_initl, E_loss, Energ

  Eimag = 0._db
  Final_optic = .false.
  Final_tddft = .false.
  Green = .false.
  Green_int = .false.
  ninitlv = n_Ec
  NRIXS = .false.
  RIXS = .true.
  Solsing = .false.
  Solsing_only = .false.
  Tddft = .false.

  allocate( Decal_initl(ninitlr) )
! Epsii_moy is the average for the second edge (when 2 edges) for example L3 in L23
  Decal_initl(:) = Epsii_a(:) - Epsii_moy

  if( E_cut_man ) then
    E_cut_rixs = E_cut_imp
  else
    E_cut_rixs = E_cut
  endif
  
  lmax_probe = nint( sqrt( real( nlm_probe, db ) ) ) - 1
  lmax_pot = nint( sqrt( real( nlm_pot, db ) ) ) - 1

  do ie = 1,nenerg_s
    if( Energ_s(ie) > E_cut_RIXS - eps10 ) exit
  end do
  nenerg = nenerg_s - ie + 1 
  ne_loss = nenerg_s - 1
  allocate( Energ(nenerg) )
  allocate( E_loss(ne_loss) )

  call Energy_Grid_RIXS(E_loss,Energ,Energ_s,icheck,ne_loss,nenerg,nenerg_s)

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
    
  do ie = 1,nenerg
    do ie_loss = 1,ne_loss

      call Call_Tau_loss(Ampl_abs,E_loss(ie_loss),Energ(ie),Energ_s,ie,ndim2,nenerg_s,nlm_probe,nlmomax,ns_dipmag,nspin,nspino, &
                         nspinp,Spinorbite,Taull)

      Enervide(1) = Energ_s(ie) + E_Fermi
      Enervide(2) = Enervide(2) - E_loss(ie_loss)
      do je = 1,n_Ec 
        Ecinetic(:,je) = Enervide(je) - V0bdc(:,1)
        if( .not. Eneg ) Ecinetic(:,je) = max( Ecinetic(:,je), Eclie )
      end do

      call Tenseur_car(Classic_irreg,coef_g,Core_resolved,Ecinetic, &
                Eimag,Energ(ie),Enervide,Eseuil,FDM_comp_m,Final_optic,Final_tddft,Full_potential,Green,Green_int,Hubb_a,Hubb_d, &
                icheck,ie,ip_max,ip0,is_g,lmax_probe,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
                n_Ec,n_oo,n_rel,n_V,nbseuil,ns_dipmag,ndim2,nenerg_tddft,ninit1,ninitl,ninitlr,ninitlv,nlm_pot,nlm_probe, &
                nlm_p_fp,nlmamax,nr,nrm,nspin,nspino,nspinp,numat,psii,r,Relativiste,Renorm,RIXS,Rmtg,Rmtsd,rof0,rot_atom_abs, &
                Rot_int,secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m, &
                secqq,secqq_m,Solsing,Solsing_only,Spinorbite,Taull,Tddft,V_hubb,V_intmax,V0bdc,Vrato,Ylm_comp)

    end do
  end do
  
  deallocate( Decal_initl, Energ, E_loss )
  
  return
end

!***********************************************************************

! Energy grid for RIXS

subroutine Energy_Grid_RIXS(E_loss,Energ,Energ_s,icheck,ne_loss,nenerg,nenerg_s)

  use declarations
  implicit none

  integer,intent(in):: icheck, ne_loss, nenerg, nenerg_s
  integer:: ie, je

  real(kind=db), dimension(nenerg_s), intent(in):: Energ_s

  real(kind=db), dimension(nenerg), intent(out):: Energ
  real(kind=db), dimension(ne_loss), intent(out):: E_loss

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

! Calculation of the inelastic scattering amplitude

subroutine Call_Tau_loss(Ampl_abs,E_loss,Energ,Energ_s,ie,ndim2,nenerg_s,nlm_probe,nlmomax,ns_dipmag,nspin,nspino,nspinp, &
                         Spinorbite,Taull)

  use declarations
  implicit none

  integer:: ie, ie_loss, iso1, iso2, isp1, isp2, ispf1, ispf2, iss1, iss2, je, je_loss, lm1, lmf1, lms1, lm2, lmf2, lms2, &
            ndim2, ne, nenerg_s, nlm_probe, nlmomax, ns_dipmag, nspin, nspino, nspinp

  complex(kind=db), dimension(nenerg_s,nlm_probe,nspinp,nlmomax,nspinp):: Ampl_abs
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull

  logical:: Spinorbite
  
  real(kind=db):: E_loss, E_val, Energ
  real(kind=db), dimension(2):: p
  real(kind=db), dimension(nenerg_s):: Energ_s
  
  E_val = Energ - E_loss
  
  do ie_loss = 1,nenerg_s
    if( Energ_s(ie_loss) > E_val - eps10 ) exit
  end do
  
  if( abs( Energ_s(ie_loss) - E_val ) > eps10 ) then
    ne = 2
    ie_loss = ie_loss - 1
    p(1) = ( Energ_s(ie_loss+1) - E_val ) / ( Energ_s(ie_loss+1) - Energ_s(ie_loss) )  
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
                if( isp1 /= isp2 ) cycle
                iss1 = min(isp1,nspin)
                iss2 = min(isp2,nspin)
              endif

              do je = 1,ne
                je_loss = ie_loss + je - 1
                do ispf1 = 1,nspinp
                  if( .not. Spinorbite .and. ispf1 /= isp1 ) cycle
                  do ispf2 = 1,nspinp
                    if( .not. Spinorbite .and. ispf2 /= isp2 ) cycle
                    do lmf1 = 1,nlmomax
                      do lmf2 = 1,nlmomax
   ! multiplication by "- img" in order it looks like an inelatic multiple scattering amplitude
                        Taull(lms1,lms2,1,1,isp1,isp2,1) = Taull(lms1,lms2,1,1,isp1,isp2,1) &
                                - img *  p(je) * Ampl_abs(je_loss,lm1,iss1,lmf1,ispf1) * Conjg( Ampl_abs(ie,lm2,iss2,lmf1,ispf2) )
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
  
  return
end

