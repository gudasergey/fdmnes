! Routines of the FDMNES package

!***********************************************************************

! Calculation of the starting energy for the self-consistency, to get at the end the Fermi energy

subroutine Starting_energ(Bulk_step,E_start,E_starta,Full_atom,iaprotoi,icheck,itypepr,lcoeur,n_atom_0, &
              n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,n_atom_proto_uc,natome,ncoeur,nenerg_coh,nlm_pot, &
              nrato,nrm,numat,nspin,ntype,Pas_SCF,Proto_calculated,psi_coeur,rato,Vcato,Vxcato, &
              Workf)

  use declarations
  implicit none

  integer:: i, iapr, icheck, ipr, ir, it, L, n, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, &
     n_atom_proto, n_atom_proto_uc, natome, nenerg_coh, nlm_pot, nr, nrm, nspin, ntype, Z

  integer, dimension(2,0:ntype):: lcoeur, ncoeur
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(natome):: iaprotoi
  integer, dimension(0:n_atom_proto):: itypepr

  character(len=3):: Atom_kind
  
  logical:: Bulk_step, Full_atom
  logical, dimension(0:n_atom_proto):: Proto_calculated

  real(kind=db):: E_KS_Dirac, E_marge, E_max_Fermi, E_start, Es, J, Pas_SCF, Workf

  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: E_starta
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato
  real(kind=db), dimension(0:nrm,2,0:ntype):: psi_coeur
  real(kind=db), dimension(2,n_atom_0_self:n_atom_ind_self):: E_KS
  real(kind=db), dimension(:), allocatable:: r, rho, psi
  real(kind=db), dimension(:,:), allocatable:: Vrs

! One selects the core level of highest energy and the valence level of lowest energy 
! The levels are ordered along the usual atomic order. Some could be inverted when the atom is in the cluster 
! E_KS is the storage of these energies

! Index 1 is for the top most core orbital
! Index 2 is for thelowest valence orbital

  if( icheck > 0 ) write (3,110)

  E_marge = 10._db / rydb
  E_start = 10000._db

  boucle_agr: do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
    else
      ipr = iapr
      if( Bulk_step .and. ipr <= n_atom_proto_uc ) cycle  
      if( .not. Proto_calculated(ipr) ) cycle
    endif
    it = itypepr(ipr)
    Z = numat(it)
    nr = nrato(it)
    allocate( r(nr) )
    allocate( Vrs(nr,nspin) )
    allocate( psi(nr) )
    allocate( rho(nr) )

    r(1:nr) = rato(1:nr,it)
  ! In the subroutine potential is in Hartree
    do ir = 1,nr
      Vrs(ir,:) = 0.5_db * r(ir) * ( Vcato(ir,1,iapr) + Vxcato(ir,1,:,iapr) )
      psi(ir) = psi_coeur(ir,1,it)
    end do

    do i = 1,2
      psi(1:nr) = psi_coeur(1:nr,i,it)
      n = ncoeur(i,it)
      L = Lcoeur(i,it)

      if( i == 1 .or. L == 0 ) then
        J = 0.5_db * ( 2 * L + 1 )
      else
        J = 0.5_db * ( 2 * L - 1 )
      endif

      E_KS(i,iapr) = E_KS_Dirac(.true.,icheck,J,L,n,nr,nspin,Z,Vrs,psi,r,rho)

      if( icheck > 2 ) then
        write(3,120) Z, n, L, nint(2*J)
        do ir = 1,nr
          write(3,130) r(ir)*bohr, psi(ir), psi_coeur(ir,i,it)
        end do
      end if

    end do

! In case of H or He there is one single core level
    if( E_KS(2,iapr) - E_KS(1,iapr) > 2*E_marge .or. numat(it) < 3 ) then
      E_starta(iapr) = E_KS(2,iapr) - E_marge
    else
      E_starta(iapr) = 0.5_db * ( E_KS(1,iapr) + E_KS(2,iapr) )
    endif

    E_start = min( E_start, E_starta(iapr) )

    E_starta(iapr) = E_starta(iapr) - eps10

    deallocate( psi, r, rho, Vrs )
    
  end do boucle_agr

  do iapr = n_atom_0_self,n_atom_ind_self
    if( .not. Full_atom ) then
      ipr = iapr
      if( Bulk_step .and. ipr <= n_atom_proto_uc ) cycle  
      if( .not. Proto_calculated(ipr) ) cycle
    endif
    if( E_start > E_KS(1,iapr) .and. E_starta(iapr) > E_start + eps10 ) then
      Es = max( E_start, E_KS(1,iapr) + E_marge )
      E_starta(iapr) = min( Es, E_starta(iapr) ) 
    endif    
  end do

! les potentiels sont references en fonction du niveau du vide alors que
! la gamme d'energie commence dessous (workf)
  E_start = E_start + Workf
  E_starta(:) = E_starta(:) + Workf
  
  E_max_Fermi = 30._db / rydb

! Evaluation du nombre de points pour la grille en energie
  nenerg_coh = nint( ( E_max_Fermi - E_start ) / Pas_SCF ) + 1

  if( icheck > 0 ) then
    if( Full_atom ) then
      Atom_kind = ' ia'
    else
      Atom_kind = 'ipr'
    endif
    write(3,140) Atom_kind
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        ipr = iapr
        if( Bulk_step .and. ipr <= n_atom_proto_uc ) cycle  
        if( .not. Proto_calculated(ipr) ) cycle
      endif
      it = itypepr(ipr)
      write(3,150) iapr, numat(it), ( ncoeur(i,it), lcoeur(i,it), E_KS(i,iapr)*rydb, i = 1,2), E_starta(iapr)*rydb
    end do
  end if

  if( icheck > 0 ) write(3,160) E_start * rydb

  return
  110 format(/' ---- En_dep --------',100('-'))
  120 format(/'  Z =',i4,', n =',i2,', L =',i2,', J =',i2,'/2',// &
              '     Radius      psi_new        psi')
  130 format(f13.7,1p,2e13.5)
  140 format(/1x,a3,'     Z    n  l     E_core  n  l     E_val    E_starta')
  150 format(i4,3x,i3,2x,2(i3,i3,f11.3),f11.3)
  160 format(/' Starting energy = ',f8.3,' eV')
end

!***********************************************************************

! Making of the energy grid for the SCF step

subroutine grille_coh(Eimag_coh,Energ_coh,E_start,Green,icheck,nenerg_coh,Pas_SCF)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer nenerg_coh

  logical Green

  real(kind=db), dimension(nenerg_coh) :: Energ_coh, Eimag_coh

  energ_coh(1) = E_start

  do ie = 2,nenerg_coh
    energ_coh(ie) = energ_coh(ie-1) + Pas_SCF
  end do

  if( Green ) then
    eimag_coh(1:nenerg_coh) = 2 * Pas_SCF
  else
    eimag_coh(1:nenerg_coh) = 0._db
  endif

  if( icheck > 1 ) then
    write(3,110)
    write(3,120)
    do ie = 1,nenerg_coh
      write(3,130) energ_coh(ie)*rydb, eimag_coh(ie)*rydb
    end do
  endif

  return
  110 format(/' ---- Grille_coh -',100('-'))
  120 format(/'    Energy   E_imag      (eV)')
  130 format(4f9.3)
end

!***********************************************************************

! Calculation of the charge inside the cluster used for SCF

subroutine Chg_agr(Bulk_step,Chargat,Chargat_init,Chg_coeur,Chg_reference,chg_open_val,Doping,Full_atom, &
                 iaprotoi,ipr_dop,iprabs,ispin_maj,itabs,icheck,itypepr,mpirank, &
                 n_atom_0_self,n_atom_ind_self,n_atom_proto,n_atom_proto_bulk,natome, &
                 nb_eq,nb_eq_2D,ngreq,nrato,nrm,nrm_self,nspin,ntype,numat,pop_open_val, &
                 Proto_calculated,psi_open_val,rato,rho_chg,rho_coeur,rhoato_init,rmtsd,SCF_mag_fix,SCF_mag_free,Sym_2D)

  use declarations
  implicit none

  integer:: i, iapr, icheck, ipr, iprabs, ipr_dop, iprint, ir, isp, ispin, it, itabs, mpirank, n, natome, n_atom_0_self, &
            n_atom_ind_self, n_atom_proto, n_atom_proto_bulk, nr, nrm, nrm_self, nspin, ntype, Sum_Z
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, nb_eq_2D, ngreq
  integer, dimension(natome):: iaprotoi, nb_eq
  integer, dimension(n_atom_0_self:n_atom_ind_self):: ispin_maj

  character(len=3):: Atom_kind
  
  logical:: Bulk_step, Doping, Full_atom, SCF_mag_fix, SCF_mag_free, Sym_2D
  logical, dimension(0:n_atom_proto):: Proto_calculated

  real(kind=db):: Charge_init, Chg, Chg_coeur_ia, Chg_sup, f, f_integr3, Rayint, Res
  
  real(kind=db), dimension(nspin):: Chg_reference
  real(kind=db), dimension(0:n_atom_proto):: chargat
  real(kind=db), dimension(0:nrm,0:ntype):: rato, rho_coeur
  real(kind=db), dimension(0:nrm_self,nspin, n_atom_0_self:n_atom_ind_self):: rho_chg, rhoato_init
  real(kind=db), dimension(0:n_atom_proto):: rmtsd
  real(kind=db), dimension(0:nrm)::  r, rh
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: ch_v, Chg_coeur, chargat_sup
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self,nspin):: chargat_init
  real(kind=db), dimension(nrm,2):: psi_open_val
  real(kind=db), dimension(2):: chg_open_val, pop_open_val

  Charge_init = 0._db
  Chg_reference(:) = 0._db
  chg_sup = 0._db
  Chg_coeur_ia = 0._db
  chargat_init(:,:) = 0._db
  chargat_sup(:) = 0._db

  if( icheck > 0 ) write(3,110)

  Sum_Z = 0
  SCF_mag_fix = .false.

  do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
      n = nb_eq(iapr)
    else
      ipr = iapr
      if( Bulk_step .and. ipr <= n_atom_proto - n_atom_proto_bulk ) cycle
      if( .not. Proto_calculated(ipr) ) cycle
      if( Sym_2D ) then
        n = nb_eq_2D(ipr)
      else
        n = ngreq(ipr)
      endif
      if( Doping .and. iapr == ipr_dop ) n = n - 1
    endif
    it = itypepr(ipr)
    nr = nrato(it)
    Rayint = rmtsd(ipr)
    r(0:nrm) = rato(0:nrm,it)

! Calcul de la charge totale de L'agregat: rhoato est la vraie densite
    do ispin = 1,nspin
      rh(0:nr) = rhoato_init(0:nr,ispin,iapr) * r(0:nr)**2
      res = quatre_pi * f_integr3(r,rh,0,nrm,rayint)
      chargat_init(iapr,ispin) = chargat_init(iapr,ispin) + res
    end do
    if( nspin == 2 .and. chargat_init(iapr,nspin) > chargat_init(iapr,1) + eps10 ) then
      ispin_maj(iapr) = 2
    else
      ispin_maj(iapr) = 1
    endif
    if( nspin == 2 .and. .not. SCF_mag_free .and. abs( chargat_init(iapr,nspin) - chargat_init(iapr,1) ) > eps10 ) &
      SCF_mag_fix = .true.

    isp = ispin_maj(iapr)
    Chg_reference(1) = Chg_reference(1) + chargat_init(iapr,isp) * n
    if( nspin == 2 ) then
      isp = 3 - ispin_maj(iapr)
      Chg_reference(2) = Chg_reference(2) + chargat_init(iapr,isp) * n
    endif

! Calcul de la charge totale de L'agregat: rho_chg est la densite
! venant des atomes exterieurs au petit agregat.
    do ispin = 1, nspin
      rh(0:nr) = rho_chg(0:nr,ispin,iapr) * r(0:nr)**2
      res = quatre_pi * f_integr3(r,rh,0,nrm,rayint)
      chargat_sup(iapr) = chargat_sup(iapr) + res
    end do
    chg_sup = chg_sup + chargat_sup(iapr) * n

! Calcul de la charge des orbitales de coeur: rho_coeur est la vraie densite
    rh(0:nr) = rho_coeur(0:nr,it) * r(0:nr)**2
    Chg_coeur(iapr) = quatre_pi * f_integr3(r,rh,0,nrm,rayint)

    Chg_coeur_ia = Chg_coeur_ia + Chg_coeur(iapr) * n
    Sum_Z = Sum_Z + n * numat(it)
    Charge_init = Charge_init + n * chargat(ipr)
  end do

! Calcul de la charge de valence:
  Chg = sum( Chg_reference(:) ) - Chg_coeur_ia

  if( icheck > 0 ) then
    write(3,120) Chg
    write(3,130) Chg_coeur_ia
    if( .not. SCF_Mag_fix ) then
      write(3,140) sum( Chg_reference(:) )
    else
      write(3,142) sum( Chg_reference(:) ), Chg_reference(:)
    endif
    write(3,145) Chg_sup
    write(3,150) Sum_Z
    write(3,160) Charge_init
    write(3,170) Sum_Z - Charge_init - sum( Chg_reference(:) )
  end if

  if( Full_atom ) then
    Atom_kind = ' ia'
  else
    Atom_kind = 'ipr'
  endif

  do iprint = 3,6,3
    if( iprint == 3 .and. icheck == 0 ) cycle
    if( iprint == 6 .and. mpirank /= 0 ) cycle
    if( nspin == 1 ) then
      write(iprint,300) Atom_kind
    else
      write(iprint,305) Atom_kind
    endif
    do iapr = n_atom_0_self,n_atom_ind_self
      ch_v(iapr) = sum(chargat_init(iapr,:)) - Chg_coeur(iapr)
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
        n = nb_eq(iapr)
      else
        ipr = iapr
        if( Sym_2D ) then
          n = nb_eq_2D(ipr)
        else
          n = ngreq(ipr)
        endif
        if( Bulk_step .and. ipr <= n_atom_proto - n_atom_proto_bulk ) cycle  
        if( .not. Proto_calculated(ipr) ) cycle
      endif
      
      it = itypepr(ipr)
      if( nspin == 1 ) then
        write(iprint,310) iapr, numat(it), n, ch_v(iapr), Chg_coeur(iapr), sum(chargat_init(iapr,:)), chargat_sup(iapr), &
                   numat(it)- sum(chargat_init(iapr,:))
      else
        write(iprint,310) iapr, numat(it), n, ch_v(iapr), Chg_coeur(iapr), sum(chargat_init(iapr,:)), &
                   chargat_init(iapr,1) - chargat_init(iapr,2), chargat_sup(iapr), numat(it)- sum(chargat_init(iapr,:))
      endif
    end do
  end do

  Sum_Z = Sum_Z - charge_init

  if( icheck > 1 ) then
    write(3,350)
    write(3,351) sum(ch_v(:))
    write(3,352) sum(Chg_coeur(:))
    write(3,353) sum(chargat_init(:,:))
  endif

  if( icheck > 2 ) then
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        ipr = iapr
        if( Bulk_step .and. ipr <= n_atom_proto - n_atom_proto_bulk ) cycle  
        if( .not. Proto_calculated(ipr) ) cycle
      endif
      it = itypepr(ipr)
      write(3,403) iapr
      if( nspin == 1 ) then
        write(3,401)
      else
        write(3,402)
      end if
      do ir = 1, nrato(it)
        write(3,405) rato(ir,it)*bohr, quatre_pi * rato(ir,it)**2 * rho_chg(ir,1:nspin,iapr)
      end do
    end do
  end if

  do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
    else
      ipr = iapr
      if( Bulk_step .and. ipr <= n_atom_proto - n_atom_proto_bulk ) cycle  
      if( .not. Proto_calculated(ipr) ) cycle
    endif
    it = itypepr(ipr)
    chargat_init(iapr,:) = real(numat(it),db) / nspin - chargat_init(iapr,:)
  end do

! coupure ad hoc au niveau de Fermi: indice 1 <=> nonexcited absorber
!                                    indice 2 <=>    excited absorber

  do i = 1,2
    if( i == 1 ) then
      it = itabs
      rayint = rmtsd(iprabs)
    else
      it = 0
      rayint = rmtsd(0)
    endif
    r(:) = rato(:,it)
    rh(0) = 0._db
    rh(1:nrm) = psi_open_val(1:nrm,i)**2
    nr = nrato(it)
    chg_open_val(i) = f_integr3(r,rh,0,nrm,rayint)
! psi_level_val etait normalise a L'unite:
    chg_open_val(i) = chg_open_val(i) * pop_open_val(i)
  end do

  if( icheck > 0 ) then
    write(3,415) chg_open_val(1), pop_open_val(1)
    write(3,416) chg_open_val(2), pop_open_val(2)
  end if

  if( icheck > 2 ) then
   write(3,500)
   do ir = 1,nrato(itabs)
     f = quatre_pi * rato(ir,itabs)**2
     write(3,405) rato(ir,itabs)*bohr, f * psi_open_val(ir,:)**2
   end do
  end if

  return
  110 format(/' ---- Chg_agr --------',100('-'))
  120 format(/' Number of valence electrons =',f9.3)
  130 format(' Number of core electrons    =',f9.3)
  140 format(' Total                       =',f9.3)
  142 format(' Total, spin maj, spin min   =',3f9.3)
  145 format(' Charge from outer sphere    =',f9.3)
  150 format(' Sum of atomic number        =',i5)
  160 format(' Initial charge              =',f9.3)
  170 format(' Cluster charge              =',f9.3)
  300 format(/1x,a3,'   Z  mult     ch_val    ch_core   ch_total     ch_out   Atom charge')
  305 format(/1x,a3,'   Z  mult     ch_val    ch_core   ch_total   ch_up-dn     ch_out     Charge')
  310 format(2i4,i5,1x,6f11.3)
  350 format(/' Reference charges in the symmetrised cluster ')
  351 format(10x,'valence electrons: ',f9.3)
  352 format(10x,'core electrons: ',f9.3)
  353 format(10x,'total number of electrons: ',f9.3)
  401 format(/' Density coming from atoms outside the cluster',/ 7x,'rato     4pi*r2*rho_chg')
  402 format(/' Density coming from atoms outside the cluster',/ 7x,'rato    4pi*r2*rho_chg(u)    4pi*r2*rho_chg(d)')
  403 format(/7x,'ia =  ',i4)
  405 format(1p,e13.5,2e17.5)
  415 format(/' Number of electron in the valence orbital of the absorbing atom:',/ 20x,'Integrated  Mulliken', &
         /'   non excited atom =',f7.3,f10.3)
  416 format('   excited atom     =',f7.3,f10.3)
  500 format(/6x,'rato   4pi*r2*psi_open_val_nonexc**2  4pi*r2*psi_open_val_exc**2')
end

!***********************************************************************

! Calculation and writing of the density of states and atomic charges
! Comparison with the reference charge and check if the Fermi level is reached.

! chargat_self = atomic charge at the current iteration
! chargat_self_s = atomic charge at the previous iteration
! ch_ia = total number of electrons including the core ones at the current iteration.

subroutine Search_Fermi(Bulk_atom_done,Bulk_step,Chg_reference,chg_open_val,chargat_self,Density,Doping,drho_self,E_cut, &
                E_Open_val,E_Open_val_exc,E_starta,Energ,E_Fermi,Energ_self,Fermi,Full_atom,hubb,Hubb_diag,iaabsi, &
                iaprotoi,i_self,icheck,ie,ie_computer,Int_statedens,ipr_dop,ispin_maj,itypei,itypepr,lamstdens, &
                Open_val,Open_val_exc,lla_state,lla2_state,lmaxat,m_hubb,mpinodes,n_atom_0,n_atom_0_self, &
                n_atom_ind,n_atom_ind_self,n_atom_proto,n_atom_proto_bulk,n_atom_proto_uc,natome,nb_eq,nb_eq_2D,nenerg,ngreq, &
                nlm_pot,nomfich_s,nrato, &
                nrm,nrm_self,nspin,nspinp,ntype,numat,occ_hubb,occ_hubb_i,pop_orb_val,Proto_calculated,rato,rho_self,rho_self_t, &
                Rmtsd,SCF_elecabs,SCF_mag_fix,Self_nonexc,Statedens,Statedens_i,Sym_2D,V_hubb,V_hubbard,Ylm_comp)

  use declarations
  implicit none

  integer:: i, i_self, ia, iaabsi, iapr, icheck, ie, ie_computer, imax, ipr, ipr_dop, iprabs, iprint, ir, isp, isp1, isp2, iss, &
    it, L, l_hubbard, l_level_val, la, lamstdens, ll, lla_state, lla2_state, lh, lm, lm0, lm1, lm2, lma, m, m_hubb,mpinodes, m1, &
    m2, n, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, n_atom_proto, n_atom_proto_bulk, n_atom_proto_uc, natome, nenerg, &
    nlm_pot, nr, nrm, nrm_self, nspin, nspinp, ntype, Numat_tot, Z

  character(len=3):: Atom_kind
  character(len=Length_name):: nomfich_s

  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, la_ipr, ll_ipr, lmaxat, nb_eq_2D, ngreq
  integer, dimension(natome):: iaprotoi, itypei, nb_eq
  integer, dimension(n_atom_0_self:n_atom_ind_self):: ispin_maj

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb

  logical, save:: Fermi_gen, Fermi_maj, Fermi_min
  logical:: Bulk_atom_done, Bulk_step, Density, Doping, Fermi, First_loop, Full_atom, Open_file, Open_val, Open_val_exc, &
    SCF_elecabs, SCF_mag_fix, Self_nonexc, Sym_2D, This_bulk_atom_done, Ylm_comp
  logical, dimension(0:n_atom_proto):: Proto_calculated, proto_done
  logical, dimension(0:ntype):: Hubb
  logical, dimension(n_atom_0_self:n_atom_ind_self):: Hubb_diag

  real(kind=db):: ch, ch_s, Charge, Charge_maj, Charge_min, Charge_val_abs, chg_lim, D_Energ, de, Dens, ds, E_cut, E_Fermi, &
    E_Open_val, E_Open_val_exc, En_f, f_integr3, poids

  real(kind=db), save:: Charge_maj_s, Charge_min_s, Charge_s, Charge_val_abs_e_s, Charge_val_abs_s, E_Fermi_maj, E_Fermi_min

  real(kind=db), dimension(2) :: chg_open_val
  real(kind=db), dimension(nspin):: Chg_reference, chg_ref
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(0:ntype):: V_hubbard
  real(kind=db), dimension(0:n_atom_proto):: Rmtsd
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: E_starta
  real(kind=db), dimension(nrm_self,nlm_pot,nspinp,n_atom_0_self:n_atom_ind_self,0:mpinodes-1):: drho_self
  real(kind=db), dimension(lla2_state,nspinp,n_atom_0:n_atom_ind):: Int_Statedens
  real(kind=db), dimension(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1):: Statedens, Statedens_i
  real(kind=db), dimension(0:lla_state,nspinp,n_atom_0:n_atom_ind):: Int_Statedens_l
  real(kind=db), dimension(0:lla_state,nspinp):: Statedens_l
  real(kind=db), dimension(nspinp,n_atom_0:n_atom_ind):: Int_Statedens_t
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: ch_ia_t, Energ_self, Energ_self_s
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self,nspin):: ch_ia, chargat_self, chargat_self_s, pop_orb_val_s, pop_orb_val
  real(kind=db), dimension(0:nrm):: r, r2, rh
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm_self,nlm_pot,n_atom_0_self:n_atom_ind_self):: rho_self_t
  real(kind=db), dimension(0:nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self):: rho_self, rho_self_s
  real(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: occ_hubb, occ_hubb_i, &
                                                                                                        occ_hubb_i_s, occ_hubb_s

  if( icheck > 1 )  then
    write(3,110)
    write(3,120) Energ(ie)*rydb
  endif

  if( Full_atom ) then
    Atom_kind = ' ia'
  else
    Atom_kind = 'ipr'
  endif
  
! Storage of the table from previous iteration. Used for the interpolation when Fermi level is reached
  if( ie == 1 ) then
    Fermi_gen = .false.
    Fermi_maj = .false.
    Fermi_min = .false.
    E_Fermi_maj = E_Fermi
    E_Fermi_min = E_Fermi
    Charge_s = 0._db; Charge_maj_s = 0._db; Charge_min_s = 0._db
    if( Bulk_atom_done ) then
      pop_orb_val(n_atom_0_self:n_atom_proto-n_atom_proto_bulk,:) = 0._db
    else
      pop_orb_val(:,:) = 0._db
    endif
    occ_hubb(:,:,:,:,:) = 0._db
    occ_hubb_i(:,:,:,:,:) = 0._db
    occ_hubb_i_s(:,:,:,:,:) = 0._db
    occ_hubb_s(:,:,:,:,:) = 0._db
  else
    chargat_self_s(:,:) = chargat_self(:,:)
    occ_hubb_i_s(:,:,:,:,:) = occ_hubb_i(:,:,:,:,:)
    occ_hubb_s(:,:,:,:,:) = occ_hubb(:,:,:,:,:)
  endif
  ch_ia(:,:) = 0._db
  ch_ia_t(:) = 0._db
  rho_self_s(:,:,:,:) = rho_self(:,:,:,:)
  Energ_self_s(:) = Energ_self(:)
  pop_orb_val_s(:,:) = pop_orb_val(:,:)

  proto_done(:) = .false.

  if( nenerg == 1 ) then
    de = 1._db
  elseif( ie == 1 ) then
    de = Energ(2) - Energ(1)
  elseif( ie == nenerg ) then
    de = Energ(ie) - Energ(ie-1)
  else
    de = 0.5_db * ( Energ(ie+1) -  Energ(ie-1) )
  endif

  do ipr = 0,n_atom_proto
    if( lamstdens > -1 ) then
      la = min(lamstdens,lmaxat(ipr))
    else
      Z = numat( itypepr(ipr) )
      if( Z < 19 ) then
        la = min(1,lmaxat(ipr))
        ll = min(2,lmaxat(ipr))
      elseif( Z > 18 .and. Z < 55 ) then
        la = min(2,lmaxat(ipr))
        ll = min(3,lmaxat(ipr))
      else
        la = min(3,lmaxat(ipr))
        ll = min(4,lmaxat(ipr))
      endif
    end if
    la_ipr(ipr) = la    ! maximum number of harmonics for the writing of the density of states
    ll_ipr(ipr) = ll    ! maximum number of harmonics for the calculation of the density of states
  end do

  boucle_ia: do ia = 1,natome

    ipr = iaprotoi(ia)
    la = la_ipr(ipr)
    ll = ll_ipr(ipr)

    This_bulk_atom_done = Bulk_atom_done .and. ipr > n_atom_proto - n_atom_proto_bulk

    if( ia == iaabsi ) iprabs = ipr

    if( Full_atom ) then
      iapr = ia
    else
      iapr = ipr
    endif

    it = itypei(ia)
    Z = numat(it)
    if( proto_done(ipr) .and. .not. Full_atom ) cycle boucle_ia
    if( Energ(ie) < E_starta(iapr) ) cycle
    Open_file = .false.
    if( ie == 1 ) then
      Open_file = .true.
    else
      if( Energ(ie-1) < E_starta(iapr) ) Open_file = .true.
    endif
     
    if( ie == 1 ) then
      First_loop = .true.
    else
      if( Energ(ie-1) < E_starta(iapr) ) then
        First_loop = .true.
      else
        First_loop = .false.
      endif
    endif

! Calculation of the density of state integrated up to the radius Rmstd of the atom

    do ir = 1,nrato(it)
      if( rato(ir,it) > Rmtsd(ipr) + eps10 ) exit
    end do
    nr = ir

    lma = (ll + 1)**2

    Statedens_l(:,:) = 0._db
    do L = 0,la
      lm1 = L**2 + 1
      lm2 = ( L + 1 )**2
      do isp = 1,nspinp
        do lm = lm1,lm2
          Statedens_l(L,isp) = Statedens_l(L,isp) + Statedens(lm,isp,lm,isp,iapr,ie_computer)
        end do
      end do
    end do
    if( nspinp == 1 ) Statedens_l(0:la,1) = 2 * Statedens_l(0:la,1)

! Integral of density of state
    do isp = 1,nspinp
      do lm = 1,lma
        Int_Statedens(lm,isp,iapr) = Int_Statedens(lm,isp,iapr) + de * Statedens(lm,isp,lm,isp,iapr,ie_computer)
      end do
    end do

    do L = 0,ll
      lm1 = L**2 + 1
      lm2 = ( L + 1 )**2
      do isp = 1,nspinp
        Int_Statedens_l(L,isp,iapr) = sum(Int_Statedens(lm1:lm2,isp,iapr))
      end do
      if( nspinp == 1 ) Int_Statedens_l(L,1,iapr) = 2 * Int_Statedens_l(L,1,iapr)
    end do

    Int_statedens_t(:,iapr) = 0._db
    do isp = 1,nspinp
      Int_Statedens_t(isp,iapr) = Int_Statedens_t(isp,iapr) + sum( Int_Statedens_l(0:la,isp,iapr) )
    end do

    if( icheck > 1 ) then
      write(3,130) Atom_kind, iapr
      lm = 0
      do L = 0,ll
        do m = -L,L
          lm = lm + 1
          do isp = 1,nspinp
            write(3,140) L, m, isp, Statedens(lm,isp,lm,isp,iapr,ie_computer), Int_Statedens(lm,isp,iapr)
          end do
        end do
      end do
      write(3,'(/A)') '    L   sum_m(Integral)'
      do L = 0,la
        write(3,150) L, Int_Statedens_l(L,1:nspinp,iapr)
      end do
      write(3,160) Int_Statedens_t(1:nspinp,iapr)
    endif

    if( ia == iaabsi .and. Density .and. icheck > 1 ) then
      call Write_stdens(.false.,.false.,Ylm_comp,Energ(ie),.false.,i_self,iapr,ie_computer,la, &
           Int_Statedens_l,lla_state,lla2_state,mpinodes,n_atom_0,n_atom_ind,nomfich_s,nspin,nspinp, &
           Open_file,Statedens,Statedens_l,Ylm_comp)
    endif

! Occupation matrix for Hubbard
    if( hubb(it) ) then
      if( nspinp == 1 ) then
        ds = 2 * de
      else
        ds = de
      endif
      lh = l_hubbard( Z )
      lm0 = lh**2 + lh + 1
      do isp1 = 1, nspinp
        do isp2 = 1, nspinp
          do m1 = -lh,lh
            do m2 = -lh,lh
              occ_hubb(m1,m2,isp1,isp2,iapr) = occ_hubb(m1,m2,isp1,isp2,iapr) &
                                             + ds * Statedens(lm0+m1,isp1,lm0+m2,isp2,iapr,ie_computer)
              occ_hubb_i(m1,m2,isp1,isp2,iapr) = occ_hubb_i(m1,m2,isp1,isp2,iapr) &
                                               + ds * Statedens_i(lm0+m1,isp1,lm0+m2,isp2,iapr,ie_computer)
            end do
          end do
        end do
      end do
    endif

! Energy of the cluster (or of the unit cell)
    if( nspinp == 1 ) then
      ds = 2 * de
    else
      ds = de
    endif
    do isp = 1,nspinp
      if( ( isp == ispin_maj(iapr) .and. Fermi_maj ) .or. ( isp /= ispin_maj(iapr) .and. Fermi_min ) ) cycle

      Dens = 0._db
      do lm = 1,lma
        Dens = Dens + Statedens(lm,isp,lm,isp,iapr,ie_computer)
      end do
      D_Energ = ds * Energ(ie) * Dens
      Energ_self(iapr) =  Energ_self(iapr) + D_Energ
    end do

    if( SCF_mag_fix .and. First_loop ) rho_self_t(:,:,iapr) = rho_self(:,:,1,iapr) + rho_self(:,:,nspin,iapr)

     r(:) = rato(:,it)
     r2(:) = rato(:,it)**2

    if( .not. This_bulk_atom_done ) then
    
      L = min( l_level_val(Z), lmaxat(ipr) )
      do isp = 1,nspin
        if( SCF_mag_fix .and. ( ( isp == ispin_maj(iapr) .and. Fermi_maj ) .or. &
            ( isp /= ispin_maj(iapr) .and. Fermi_min ) ) ) cycle
        if( nspinp /= nspin ) then
          pop_orb_val(iapr,isp) = sum( Int_statedens_l(L,:,iapr) )
        else
          pop_orb_val(iapr,isp) = Int_statedens_l(L,isp,iapr)
        endif
      end do

      if( First_loop ) then
        do isp = 1, nspin
          rh(0:nr) = rho_self(0:nr,1,isp,iapr) * r2(0:nr)
          chargat_self(iapr,isp) = Real( Z, db ) / nspin - quatre_pi * f_integr3(r,rh,0,nrm,Rmtsd(ipr))
          chargat_self_s(iapr,isp) = chargat_self(iapr,isp)
        end do

        if( icheck > 1 ) write(3,170) Atom_kind, iapr, chargat_self(iapr,:)
        if( icheck > 2 ) then
          write(3,180) Atom_kind, iapr, Z
          do ir = 1,nr
            write(3,190) rato(ir,it)*bohr, ( quatre_pi * r2(ir) * rho_self(ir,lm,:,iapr), lm = 1,nlm_pot )
          end do
        end if
      endif

      if( nspin == 1 .and. .not. nspinp == 2 ) then
        ds = 2 * de
      else
        ds = de
      endif
   ! drho is in "solution" basis, here we make the approximation that it is like "spin" basis even if it is a bad quantum number
      do isp = 1,nspinp
        iss = min(isp,nspin)
        if( SCF_mag_fix .and. .not. Fermi_gen ) rho_self_t(1:nr,:,iapr) = rho_self_t(1:nr,:,iapr) &
                      + ds * drho_self(1:nr,:,isp,iapr,ie_computer)
        if( ( isp == ispin_maj(iapr) .and. Fermi_maj ) .or. ( isp /= ispin_maj(iapr) .and. Fermi_min ) ) cycle
        rho_self(1:nr,:,iss,iapr) = rho_self(1:nr,:,iss,iapr) + ds * drho_self(1:nr,:,isp,iapr,ie_computer)
      end do

    endif

! Calculation of atomic charge
    do isp = 1,nspin
      rh(0:nr) = rho_self(0:nr,1,isp,iapr) * r2(0:nr)
      ch_ia(iapr,isp) = quatre_pi * f_integr3(r,rh,0,nrm,Rmtsd(ipr))
      if( .not. This_bulk_atom_done ) chargat_self(iapr,isp) = Real( Z, db ) / nspin - ch_ia(iapr,isp)
    end do
    if( SCF_mag_fix .and. .not. Fermi_gen ) then
      rh(0:nr) = rho_self_t(0:nr,1,iapr) * r2(0:nr)
      ch_ia_t(iapr) = quatre_pi * f_integr3(r,rh,0,nrm,Rmtsd(ipr))
    endif

    if( icheck > 1 ) then
      if( nspin == 1 ) then
        write(3,200) Atom_kind, iapr, chargat_self(iapr,1), ch_ia(iapr,1), Int_statedens_t(1,iapr)
      else
        write(3,210) Atom_kind, iapr, ( isp, chargat_self(iapr,isp), ch_ia(iapr,isp), &
                     Int_statedens_t(min(isp,nspin),iapr), isp = 1,nspin )
      endif
    endif
    if( icheck > 2 ) then
      write(3,180) Atom_kind, iapr, Z
      do ir = 1,nr
        write(3,190) rato(ir,it)*bohr, ( quatre_pi * r2(ir) * rho_self(ir,lm,:,iapr), lm = 1,nlm_pot )
      end do
    end if

    proto_done(ipr) = .true.

  end do boucle_ia

! Interpolation:
! E_f = E_i*(ch_ref - ch_i-1)/(ch_i-ch_i-1) + E_i-1*(ch_i - ch_ref)/(ch_i-ch_i-1)
! ch_s = total number of electrons at previous iteration
! ch = total number of electrons
! Chg_reference = charge to compare with the reference

! Fermi level evaluation

  Z = numat( itypei(iaabsi) )
  L = l_level_val(Z)
  Charge = 0._db
  Charge_maj = 0._db
  Charge_min = 0._db
  chg_ref(:) = Chg_reference(:)
  Numat_tot = 0
  if( Full_atom ) then
    Charge_val_abs = sum( pop_orb_val(iaabsi,:) )
  else
    Charge_val_abs = sum( pop_orb_val(iprabs,:) )
  endif

  do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      n = nb_eq(iapr)
      Z = numat( itypei(iapr) )
    else
      if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
      if( .not. Proto_calculated(iapr) ) cycle
      if( Sym_2D ) then
        n = nb_eq_2D(iapr)
      else
        n = ngreq(iapr)
      endif
      if( Doping .and. iapr == ipr_dop ) n = n - 1
      Z = numat( itypepr(iapr) )
    endif
    if( SCF_mag_fix ) then
      isp = ispin_maj(iapr)
      Charge_maj = Charge_maj + n * ch_ia(iapr,isp)
      isp = 3 - ispin_maj(iapr)
      Charge_min = Charge_min + n * ch_ia(iapr,isp)
      Charge = Charge + n * ch_ia_t(iapr)
    else
      Charge = Charge + n * sum( ch_ia(iapr,:) )
    endif
    Numat_tot = Numat_tot + n * Z
  end do

  if( icheck > 1 ) then
    if( nspin == 1 ) then
      write(3,220) Charge, chg_ref(1), Numat_tot
    else
      write(3,230) Charge, sum( chg_ref(:) ), chg_ref(:), Numat_tot
    endif
  endif

  if( icheck > 1 ) then
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Energ(ie) < E_starta(iapr) ) cycle
      if( Full_atom) then
        it = itypei(iapr)
      else
        if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
        if( .not. Proto_calculated(iapr) ) cycle
        it = itypepr(iapr)
      endif
      if( .not. hubb(it) ) cycle
      Z = numat(it)
      L = l_hubbard( Z )
      call Write_occ_hubb(Atom_kind,iapr,L,m_hubb,n_atom_0_self,n_atom_ind_self,nspinp,occ_hubb,occ_hubb_i,Z)
    end do
  endif

  if( SCF_mag_fix ) then
    imax = 4
  else
    imax = 2
  endif

  do i = 0,imax

    select case(i)
      case(0)
        if( Open_val ) cycle
        ch = Charge_val_abs
        chg_lim = chg_open_val(1)
      case(1)
        if( Open_val_exc .or. self_nonexc ) cycle
        ch = Charge_val_abs
        chg_lim = chg_open_val(2)
      case(2)
        if( Fermi_gen ) Cycle
        ch = Charge
        chg_lim = sum( chg_ref(:) )
      case(3)
        if( Fermi_maj ) Cycle
        ch = Charge_maj
        chg_lim = chg_ref(1)
      case(4)
        if( Fermi_min ) Cycle
        ch = Charge_min
        chg_lim = chg_ref(nspin)
    end select

    if( ch < chg_lim ) then

      select case(i)
        case(0)
          Charge_val_abs_s = ch
        case(1)
          Charge_val_abs_e_s = ch
        case(2)
          Charge_s = ch
        case(3)
          Charge_maj_s = ch
        case(4)
          Charge_min_s = ch
      end select

! We have reached the Fermi level
    else

      select case(i)
        case(0)
          Open_val = .true.
          if( self_nonexc ) Open_val_exc = .true.
          ch_s = Charge_val_abs_s
        case(1)
          Open_val_exc = .true.
          ch_s = Charge_val_abs_e_s
        case(2)
          Fermi_gen = .true.
          if( .not. SCF_mag_fix ) then
            Fermi_maj = .true.
            Fermi_min = .true.
          endif
          ch_s = Charge_s
        case(3)
          Fermi_maj = .true.
          ch_s = Charge_maj_s
        case(4)
          Fermi_min = .true.
          ch_s = Charge_min_s
      end select
      Fermi = Fermi_gen .and. Fermi_maj .and. Fermi_min

! Interpolation, une fois qu'on a atteint le niveau de Fermi:
      poids = ( chg_lim - ch_s ) / ( ch - ch_s )

      if( ie == 1 ) then
        En_f = Energ(ie)
      else
        En_f = Energ(ie) * poids + Energ(ie-1) * ( 1 - poids )
      end if

      select case(i)
        case(0)
          E_Open_val = En_f
          if( Self_nonexc ) E_Open_val_exc = En_f
          Charge_val_abs_s = ch
          cycle
        case(1)
          E_Open_val_exc = En_f
          Charge_val_abs_e_s = ch
          cycle
        case(2)
          E_Fermi = En_f
 ! On demand or when "calculation and SCF" excited and L = 2 or 3.
          if( Open_val_exc .and. scf_elecabs ) then
            E_cut = E_Open_val_exc
            E_cut = max( E_cut, En_f - 2._db/Rydb ) ! when E_Open_val_exc too far from E_Fermi
          else
            E_cut = E_Fermi
          endif
        case(3)
          E_Fermi_maj = En_f
        case(4)
          E_Fermi_min = En_f
      end select

      do iapr = n_atom_0_self,n_atom_ind_self
        if( .not. Full_atom ) then
          if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
          if( .not. Proto_calculated(iapr) ) cycle
        endif
        do isp = 1,nspin
          if( i == 2 .and. SCF_mag_fix ) cycle
          if( i == 3 .and. isp /= ispin_maj(iapr) ) cycle
          if( i == 4 .and. isp == ispin_maj(iapr) ) cycle
          chargat_self(iapr,isp) = chargat_self(iapr,isp) * poids + chargat_self_s(iapr,isp) * ( 1 - poids )
          rho_self(:,:,isp,iapr) = rho_self(:,:,isp,iapr) * poids + rho_self_s(:,:,isp,iapr) * ( 1 - poids )
          pop_orb_val(iapr,isp) = pop_orb_val(iapr,isp) * poids + pop_orb_val_s(iapr,isp) * ( 1 - poids )
        end do
      end do
      do iapr = n_atom_0_self,n_atom_ind_self
        if( .not. Full_atom ) then
          if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
          if( .not. Proto_calculated(iapr) ) cycle
        endif
        if( i == 2 .and. SCF_mag_fix ) cycle
        occ_hubb(:,:,:,:,iapr) = occ_hubb(:,:,:,:,iapr) * poids + occ_hubb_s(:,:,:,:,iapr) * ( 1 - poids )
        occ_hubb_i(:,:,:,:,iapr) = occ_hubb_i(:,:,:,:,iapr) * poids + occ_hubb_i_s(:,:,:,:,iapr) * ( 1 - poids )
      end do

      if( .not. Fermi ) cycle

      Energ_self(:) = Energ_self(:) * poids + Energ_self_s(:) * ( 1 - poids )

      if( Full_atom ) then
        Charge_val_abs = sum( pop_orb_val(iaabsi,:) )
      else
        Charge_val_abs = sum( pop_orb_val(iprabs,:) )
      endif

      if( icheck == 1 .and. i_self == 1 ) write(3,110)

      do iprint = 3,6,3
        if( icheck == 0 .and. iprint == 3 ) cycle
        if( SCF_mag_fix ) then
          write(iprint,289) i_self, E_Fermi*rydb, E_Fermi_maj*rydb, E_Fermi_min*rydb, Charge_val_abs 
        else
          write(iprint,290) i_self, E_Fermi*rydb, Charge_val_abs
        endif
        if( icheck > 1 .and. iprint == 3 ) then
          if( Open_val_exc .and. .not. self_nonexc ) then
            write(iprint,292) E_Open_val_exc*rydb
            write(iprint,293) E_Open_val*rydb
          elseif( Open_val ) then
            write(iprint,293) E_Open_val*rydb
          endif
        endif
        
        if( icheck < 2 .and. iprint == 3 ) cycle
        if( nspin == 1 ) then
          write(iprint,300) Atom_kind
        else
          write(iprint,305) Atom_kind
        endif
        do iapr = n_atom_0_self,n_atom_ind_self
          if( Full_atom ) then
            ipr = iaprotoi(iapr)
          else
            if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
            ipr = iapr
            if( .not. Proto_calculated(ipr) ) cycle
          endif
          Z = numat( itypepr(ipr) )
          L = min( l_level_val(Z), lmaxat(ipr) )
          if( nspin == 2 ) then
            write(iprint,310) iapr, Z, Energ_self(iapr)*rydb, sum(chargat_self(iapr,:)), &
                 chargat_self(iapr,2) - chargat_self(iapr,1), sum(pop_orb_val(iapr,:)), L, Rmtsd(ipr)*bohr
          else
            write(iprint,320) iapr, Z, Energ_self(iapr)*rydb, sum(chargat_self(iapr,:)), sum(pop_orb_val(iapr,:)), &
                L, Rmtsd(ipr)*bohr
          endif
        end do
      end do
      if( icheck > 1 ) then
        do iapr = n_atom_0_self,n_atom_ind_self
          if( Full_atom ) then
            ipr = iaprotoi(iapr)
          else
           if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
           ipr = iapr
           if( .not. Proto_calculated(ipr) ) cycle
          endif
          it = itypepr(ipr)
          Z = numat( it )

          write(3,180) Atom_kind, iapr, Z
          do ir = 1,nrato(it)
            if( rato(ir,it) > Rmtsd(ipr) + eps10 ) exit
            write(3,190) rato(ir,it)*bohr, ( quatre_pi * ( rato(ir,it)**2 ) * rho_self(ir,lm,:,iapr), lm = 1,nlm_pot )
          end do
        end do
      endif

! Matrice de Hubbard

      Hubb_diag(:) = .true.
      do iapr = n_atom_0_self,n_atom_ind_self

        if( Full_atom ) then
          ipr = iaprotoi(iapr)
        else
          if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
          ipr = iapr
          if( .not. Proto_calculated(ipr) ) cycle
        endif
        it = itypepr(ipr)
        Z = numat( it )
        if( .not. hubb(it) ) cycle
        L = l_hubbard( Z )

! Occupancy counts for exchange, it is thus per spin
        if( nspinp == 1 ) then
          occ_hubb(:,:,:,:,iapr) = 0.5_db * occ_hubb(:,:,:,:,iapr)
          occ_hubb_i(:,:,:,:,iapr) = 0.5_db * occ_hubb_i(:,:,:,:,iapr)
        endif

        V_hubb(:,:,:,:,iapr) = - V_hubbard(it) * cmplx( occ_hubb(:,:,:,:,iapr), occ_hubb_i(:,:,:,:,iapr) )
        do isp = 1,nspinp
          do m = -L,L
            V_hubb(m,m,isp,isp,iapr) = V_hubb(m,m,isp,isp,iapr) + 0.5_db * V_hubbard(it)
          end do
        end do

        if( icheck > 1 ) then
          call Write_occ_hubb(Atom_kind,iapr,L,m_hubb,n_atom_0_self,n_atom_ind_self,nspinp,occ_hubb,occ_hubb_i,Z)
          call Write_V_hubb(Atom_kind,iapr,it,L,m_hubb,n_atom_0_self,n_atom_ind_self,nspinp,ntype,V_hubb,V_hubbard,Z)
        endif

        boucle_m: do m1 = -L,L
          do isp1 = 1,nspinp
            do m2 = -L,L
              do isp2 = 1,nspinp
                if( m1 == m2 .and. isp1 == isp2 ) cycle
                if( abs( V_hubb(m1,m2,isp1,isp2,iapr) ) < eps10 ) cycle
                Hubb_diag(iapr) = .false.
                exit boucle_m
              end do
            end do
          end do
        end do boucle_m

      end do

    endif

  end do

  return
  110 format(/' ---- Search_Fermi ----',100('-'))
  120 format(/' Energy =',f10.3,' eV')
  130 format(/'  L  m is  Density of states   Integral   ',a3,' =',i3)
  140 format(3i3,3f15.7)
  150 format(i5,2f15.7)
  160 format(' Total =',f12.7,f15.7)
  170 format(/'  Before integration: ',a3,' = ',i3,'  charge_self = ', 2f10.5)
  180 format(/1x,a3,' =',i3,', Z =',i3/ '   Radius_(A) 4*pi*r2*Rho_self')
  190 format(1p,30e13.5)
  200 format(/1x,a3,' =',i3,', Charge_self =',f10.5,', ch_ia =',f10.5,', Int_state_t =',f10.5)
  210 format(/1x,a3,' =',i3,2(', isp =',i2,', Charge_self =',f10.5,', ch_ia =',f10.5,', Int_state_t =',f10.5))
  220 format(/' Total charge at the current iteration =',f10.3,', Reference charge =',f10.3,', Sum of atomic numbers =',i5)
  230 format(/' Total charge at the current iteration =',f10.3,', Reference charge =',f10.3, &
              ', spin maj, spin min =',2f10.3,', Sum of atomic numbers =',i5)
  289 format(/' Cycle',i4,',  Fermi Energy =',f8.3,' eV,  maj, min =',2f8.3,' eV, occupancy val abs =',f7.3)
  290 format(/' Cycle',i4,',  Fermi energy =',f8.3,' eV, Occupancy val abs =',f7.3)
  292 format(9x,'Level val excite =',f8.3,' eV')
  293 format(9x,'Level val absorb =',f8.3,' eV')
  300 format(/6x,a3,'   Z   Energy_KS      Charge  pop_orb_val(L)   L    Radius')
  305 format(/6x,a3,'   Z   Energy_KS      Charge       up-dn  pop_orb_val(L)   L    Radius')
  310 format(5x,2i4,4f12.3,i8,f10.5)
  320 format(5x,2i4,3f12.3,i8,f10.5)
end

!********************************************************************************

subroutine Write_occ_hubb(Atom_kind,iapr,L,m_hubb,n_atom_0_self,n_atom_ind_self,nspinp,occ_hubb,occ_hubb_i,Z)

  use declarations
  implicit none

  character(len=3):: Atom_kind

  integer:: iapr, isp, isp1, isp2, L, m, m_hubb, m1, m2, n_atom_0_self, n_atom_ind_self, nspinp, Z

  real(kind=db):: Pop 
  real(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: occ_hubb, occ_hubb_i

  Pop = 0._db
  do isp = 1,nspinp
    do m = -L,L
      Pop = Pop + occ_hubb(m,m,isp,isp,iapr)
    end do
  end do
! Occupancy count for exchange, it is thus per spin. Pop is the total occupancy
  if( nspinp == 1 ) Pop = 2._db * Pop

  if( nspinp == 1 ) then
    write(3,266) Atom_kind, iapr, Z, ( m2, m2 = -L,L ) 
  else
    write(3,267) Atom_kind, iapr, Z, ( ( m2, isp2, m2 = -L,L ), isp2 = 1,nspinp )
  endif
  do isp1 = 1,nspinp
    do m1 = -L,L
      if( nspinp == 1 ) then
        write(3,268) m1, (( occ_hubb(m1,m2,isp1,isp2,iapr), occ_hubb_i(m1,m2,isp1,isp2,iapr), m2 = -L,L ), &
                                isp2 = 1,nspinp )
      else
        write(3,269) m1, isp1, (( occ_hubb(m1,m2,isp1,isp2,iapr), occ_hubb_i(m1,m2,isp1,isp2,iapr), m2 = -L,L ), &
                                  isp2 = 1,nspinp )
      endif 
    end do
  end do
  write(3,350) Pop

  return
  266 format(/' Hubbard occupation matrix for ',a3,' =',i3,', Z =',i3,/'  m',2(1x,7(10x,i3,10x)))
  267 format(/' Hubbard occupation matrix for ',a3,' =',i3,', Z =',i3,/' m1 s1',2(1x,7(9x,2i3,8x)))
  268 format(i3,7(1x,2f11.7))
  269 format(2i3,2(1x,7(1x,2f11.7)))
  350 format(/' Total number of electron, summed over spin:',f9.5)
end

!********************************************************************************

subroutine Write_V_hubb(Atom_kind,iapr,it,L,m_hubb,n_atom_0_self,n_atom_ind_self,nspinp,ntype,V_hubb,V_hubbard,Z)

  use declarations
  implicit none

  integer:: iapr, isp1, isp2, it, L, m_hubb, m1, m2, n_atom_0_self, n_atom_ind_self, nspinp, ntype, Z

  character(len=3):: Atom_kind

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb

  real(kind=db), dimension(0:ntype):: V_hubbard

  if( nspinp == 1 ) then
    write(3,360) Atom_kind, iapr, Z, V_hubbard(it) * Rydb, ( m2, m2 = -L,L )
  else
    write(3,361) Atom_kind, iapr, Z, V_hubbard(it) * Rydb, ( ( m2, isp2, m2 = -L,L ), isp2 = 1,nspinp )
  endif
  do isp1 = 1,nspinp
    do m1 = -L,L
      if( nspinp == 1 ) then
        write(3,268) m1, (( V_hubb(m1,m2,isp1,isp2,iapr) * Rydb, m2 = -L,L ), isp2 = 1,nspinp )
      else
        write(3,269) m1, isp1, (( V_hubb(m1,m2,isp1,isp2,iapr) * Rydb, m2 = -L,L ), isp2 = 1,nspinp ) 
      endif
    end do
  end do

  return
  268 format(i3,7(1x,2f11.7))
  269 format(2i3,2(1x,7(1x,2f11.7)))
  360 format(/' Hubbard potential for ',a3,' =',i3,', Z =',i3,', Hubbard parameter U - J =',f5.2,' eV',/ &
              '  m',2(1x,7(10x,i3,10x)))
  361 format(/' Hubbard potential for ',a3,' =',i3,', Z =',i3,', Hubbard parameter U - J =',f5.2,' eV',/ &
              ' m1 s1',2(1x,7(9x,2i3,8x)))
end

!********************************************************************************

! Routine calcculating the Kohn Sham energy of the top most core orbital;
! One considers that the wave function of this orbital is not modified along the SCF cycles. 
! Nevertheless the potential and thus its energy is changing.
! The Approximation is that thie shift is the same for all the core levels.

! Result will be used to calculate the cluster (or unit cell) energy
! Chg_coeur is the total charge of the core, that is the number of electrons which are not valence

subroutine Energ_KS_core(Bulk_step,E_KS_core,Full_atom,iaprotoi,icheck,itypepr,Jseuil,lcoeur,Lseuil,n_atom_0, &
              n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,n_atom_proto_uc,natome,ncoeur,nlm_pot,nrato, &
              nrm,nspin,nseuil,ntype,numat,Proto_calculated,psi_coeur,rato,Vcato,Vxcato)

  use declarations
  implicit none

  integer i_j, iapr, icheck, ipr, ir, it, Jseuil, L, Lseuil, n, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, &
          n_atom_proto, n_atom_proto_uc, natome, nlm_pot, nr, nrm, nseuil, nspin, nspino, nspinp, ntype, Z
  integer,dimension(2,0:ntype):: lcoeur, ncoeur
  integer,dimension(natome):: iaprotoi
  integer,dimension(0:n_atom_proto):: itypepr
  integer,dimension(0:ntype):: nrato, numat

  logical:: Bulk_step, Full_atom, Ylm_comp
  logical, dimension(0:n_atom_proto):: Proto_calculated
            
  real(kind=db):: E_KS_Dirac, E_KS, J, Level_occupancy 
  real(kind=db), dimension(nspin):: V0bd
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: Core_occupancy, E_KS_core, E_KS_core_s
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato
  real(kind=db), dimension(0:nrm,2,0:ntype):: psi_coeur
  real(kind=db), dimension(:), allocatable:: r, psi, psi_s, rho
  real(kind=db), dimension(:,:), allocatable:: Vrs

! Only the spherical part of the potential
  Ylm_comp = .false.
  V0bd(:) = 0._db
  nspino = 2
  nspinp = 2
  
  E_KS_core_s(:) = E_KS_core(:)
  E_KS_core(:) = 0._db
  Core_occupancy(:) = 0._db

  do iapr = n_atom_0_self,n_atom_ind_self

    if( Full_atom ) then
      ipr = iaprotoi(iapr)
    else
      if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
      ipr = iapr
      if( .not. Proto_calculated(ipr) ) cycle
    endif
    it = itypepr(ipr)
    Z = numat(it)
    nr = nrato(it)
    allocate( r(nr) )
    allocate( Vrs(nr,nspin) )
    allocate( psi(nr) )
    allocate( psi_s(nr) )
    allocate( rho(nr) )
    
    r(1:nr) = rato(1:nr,it)

  ! In the subroutine potential is in Hartree
    do ir = 1,nr
      Vrs(ir,:) = 0.5_db * r(ir) * ( Vcato(ir,1,iapr) + Vxcato(ir,1,:,iapr) )
      psi(ir) = psi_coeur(ir,1,it)
    end do

    do n = 1,ncoeur(1,it)
      do L = 0,n-1
      
        do i_j = 1,2  

          if( L == 0 .and. i_j == 1 ) cycle

          if( i_j == 1 ) then
            J = 0.5_db * ( 2 * L - 1 )
          else
            J = 0.5_db * ( 2 * L + 1 )
          endif 

          if( it == 0 .and. L == Lseuil .and. n == nseuil .and. &
             ( ( Jseuil == 1 .or. Jseuil == 2 .or. Jseuil == 4 .or. Jseuil == 6 ) .and. i_j == 1 ) ) then
            Level_occupancy = 2*J    ! missing electron 
          else 
            Level_occupancy = 2*J + 1
          endif

          Core_occupancy(iapr) = Core_occupancy(iapr) + Level_Occupancy
          
          psi_s(:) = psi(:)
        
          E_KS =  E_KS_Dirac(.true.,icheck,J,L,n,nr,nspin,Z,Vrs,psi,r,rho)

          E_KS_core(iapr) = E_KS_core(iapr) + Level_occupancy * E_KS
          if( icheck > 1 ) then
            write(3,*) 
            write(3,110) Z, n, L, nint(2*J), E_KS * rydb 
          endif
          if( icheck > 2 ) then
            write(3,'(/A)') '    Radius      old_psi      new_psi' 
            do ir = 1,nr
              write(3,'(f11.7,1p,2e13.5)') r(ir)*bohr, psi_s(ir), psi(ir) 
            end do
          endif        
        
        end do
     
        if( n == ncoeur(1,it) .and. L == Lcoeur(1,it) ) exit
      end do
    end do
    
    deallocate( psi, psi_s, r, rho, Vrs )
    
  end do ! en of loop over atoms

  if( icheck > 1 ) then
    write(3,'(/A)') ' iapr  Z  Occ.     E_KS      Delta E_KS'
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
        ipr = iapr
        if( .not. Proto_calculated(ipr) ) cycle
      endif
      it = itypepr(ipr)
      Z = numat(it)
      write(3,120) iapr, Z, nint(Core_occupancy(iapr)), E_KS_core(iapr) * rydb, ( E_KS_core(iapr) - E_KS_core_s(iapr) ) * rydb
    end do
  endif
  
  return
  110 format(' Z =',i3,' n =',i2,', L =',i2,', J =',i2,'/2, E =',f15.5,' eV')
  120 format(2i4,i5,2f13.3)
end

!********************************************************************************

! Routine calcculating the Kohn Sham energy of the top most core orbital;
! One considers that the wave function of this orbital is not modified along the SCF cycles. 
! Nevertheless the potential and thus its energy is changing.
! The Approximation is that thie shift is the same for all the core levels.

! Result will be used to calculate the cluster (or unit cell) energy
! Chg_coeur is the total charge of the core, that is the number of electrons which are not valence

subroutine Eps_coeur(Bulk_step,Chg_coeur,E_coeur,E_coeur_s,Full_atom,iaprotoi,icheck,itypepr,lcoeur,n_atom_0,n_atom_0_self, &
              n_atom_ind,n_atom_ind_self,n_atom_proto,n_atom_proto_uc,natome,ncoeur,nlm_pot,nrato,nrm,nspin,ntype,numat, &
              Proto_calculated,psi_coeur,rato,Relativiste,Rmtg,Rmtsd,V_intmax,Vcato,Vxcato)

  use declarations
  implicit none

  integer iapr, icheck, ipr, ir, it, lmax_pot_loc, natome, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, &
          n_atom_proto, n_atom_proto_uc, nlm_pot, nlm_pot_loc, nr, nrm, nspin, ntype

  logical:: Bulk_step, Full_atom, Relativiste, Ylm_comp
  logical, dimension(0:n_atom_proto):: Proto_calculated

  integer,dimension(2,0:ntype):: lcoeur, ncoeur
  integer,dimension(natome):: iaprotoi
  integer,dimension(0:n_atom_proto):: itypepr
  integer,dimension(0:ntype):: nrato, numat

  real(kind=db):: psiHpsi, res, V_intmax, Vxc
  real(kind=db), dimension(nspin):: V0bd
  real(kind=db), dimension(0:n_atom_proto):: Rmtg, Rmtsd
  real(kind=db),dimension(n_atom_0_self:n_atom_ind_self):: Chg_coeur
  real(kind=db),dimension(n_atom_0_self:n_atom_ind_self):: E_coeur, E_coeur_s
  real(kind=db),dimension(0:nrm,0:ntype):: rato
  real(kind=db),dimension(nrm):: r, psi
  real(kind=db),dimension(nrm,1,1):: pot
  real(kind=db),dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db),dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato
  real(kind=db),dimension(0:nrm,2,0:ntype):: psi_coeur

! Only the spherical part of the potential
  Ylm_comp = .false.
  lmax_pot_loc = 0
  nlm_pot_loc = 1
  V0bd(:) = 0._db

  E_coeur(:) = 0._db

  do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
    else
      if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
      ipr = iapr
      if( .not. Proto_calculated(ipr) ) cycle
    endif
    it = itypepr(ipr)
    nr = nrato(it)
    r(1:nrm) = rato(1:nrm,it)
! Average on both spins if needed:
    do ir = 1,nrm
      Vxc = sum( Vxcato(ir,1,1:nspin,iapr) ) / nspin
      pot(ir,1,1) = Vcato(ir,1,iapr) + Vxc
      psi(ir) = psi_coeur(ir,1,it)
    end do
    res =  psiHpsi(.false.,icheck,1,lcoeur(1,it),lmax_pot_loc,0,ncoeur(1,it),nlm_pot_loc,nr,nrm,1,1, &
              1,1,numat(it),V0bd,pot,psi,r,Relativiste,Rmtg(ipr),Rmtsd(ipr),.false.,V_intmax,Ylm_comp)
    res = res - E_coeur_s(iapr)
    E_coeur(iapr) = res * nint( Chg_coeur(iapr) )

    if( icheck > 1 ) then
      if( iapr == 1 ) write(3,'(/A)') ' iapr      E_KS      Delta E_KS'
      write(3,110) iapr, ( E_coeur(iapr) + E_coeur_s(iapr)*nint( Chg_coeur(iapr) ) ) * rydb, E_coeur(iapr) * rydb
    endif

  end do

  return
  110 format(i4,2f13.3)
end

!***********************************************************************

! Calculation of the total energy

subroutine Energ_DFT(Bulk_step,Doping,En_cluster,Energ_self,E_KS_core,Excato,Full_atom,Hubb,iaprotoi,icheck,ipr_dop, &
           itypepr,m_hubb,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,n_atom_proto_uc,natome,nb_eq,nb_eq_2D, &
           ngreq,nlm_pot,nrm,nrm_self,nrato,nspin,nspinp,ntype,numat,occ_hubb,Proto_calculated,rato,rho_self,rmtsd,Sym_2D, &
           V_hubb,Vcato,Vxcato)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: iapr, icheck, ipr, ipr_dop, ir, isp, ispin, it, m, m_hubb, n, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, &
            n_atom_proto, n_atom_proto_uc, natome, nlm_pot, nr, nrm, nrm_self, nspin, nspinp, ntype

  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(natome):: iaprotoi, nb_eq
  integer, dimension(0:n_atom_proto):: itypepr, nb_eq_2D, ngreq

  character(len=3):: Atom_kind
  
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb

  logical:: Bulk_step, Doping, Full_atom, Sym_2D
  logical, dimension(0:ntype):: Hubb
  logical, dimension(0:n_atom_proto):: Proto_calculated

  real(kind=db):: E_elec_agr, E_exc_agr, E_hubbard_agr, E_kin_agr, E_KS_core_agr, E_Vxc_agr, En_cluster, &
                  Energ_self_KS_agr, f_integr3, fac, U_agr

  real(kind=db), dimension(0:nrm_self,nlm_pot,nspin, n_atom_0_self:n_atom_ind_self):: rho_self
  real(kind=db), dimension(nrm,n_atom_0:n_atom_ind):: Excato
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: E_KS_core, Energ_self, Energ_self_KS, &
                     E_exc, E_elec, E_hubbard, E_kin, E_Vxc, U
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vxcato
  real(kind=db), dimension(nrm):: Vhartree, V_Z
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind)::Vcato
  real(kind=db), dimension(0:nrm):: fct1, fct2, r, r2
  real(kind=db), dimension(0:n_atom_proto):: rmtsd
  real(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: occ_hubb

  if( nrm /= nrm_self ) return  ! because rho_self is not calculated

  Energ_self_KS(:) = Energ_self(:)

  boucle_ia: do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
    else
      ipr = iapr
      if( Bulk_step .and. ipr <= n_atom_proto_uc ) cycle  
      if( .not. Proto_calculated(ipr) ) cycle
    endif
    it = itypepr(ipr)
    nr = nrato(it)
    r(:) = rato(:,it)
    r2(1:nr) = r(1:nr)**2
    r2(0) = 0._db

    fct1(:) = 0._db; fct2(:) = 0._db
    do ir = 1,nr
      fct1(ir) = sum( rho_self(ir,1,1:nspin,iapr) ) * r2(ir)
    end do

!  Hartree potential is the coulonmb potential of just the electrons, not the nucleus.
    V_Z(1:nr) = - 2 * numat(it) / r(1:nr)
    Vhartree(1:nr) = Vcato(1:nr,1,iapr) - V_Z(1:nr)
    
    fct2(1:nr) = fct1(1:nr) * Excato(1:nr,iapr)
    E_exc(iapr) = quatre_pi * f_integr3(r,fct2,0,nrm,Rmtsd(ipr))

    fct2(1:nr) = fct1(1:nr) * Vhartree(1:nr)
    E_elec(iapr) = 0.5_db * quatre_pi * f_integr3(r,fct2,0,nrm,Rmtsd(ipr))
   
    fct2(1:nr) = fct1(1:nr) * ( Vcato(1:nr,1,iapr) + V_Z(1:nr) )
    U(iapr) = 0.5_db * quatre_pi * f_integr3(r,fct2,0,nrm,Rmtsd(ipr)) 

    E_Vxc(iapr) = 0._db
    do ispin = 1, nspin
      fct2(1:nr) = Vxcato(1:nr,1,ispin,iapr) * rho_self(1:nr,1,ispin,iapr) * r2(1:nr)
      E_Vxc(iapr) = E_Vxc(iapr) + quatre_pi * f_integr3(r,fct2,0,nrm,Rmtsd(ipr))
    end do

    E_hubbard(iapr) = 0._db
    fac = 1._db / ( 3 - nspin )
    if( Hubb(it) ) then
      do isp = 1,nspinp
        do m = -m_hubb,m_hubb
          E_hubbard(iapr) = E_hubbard(iapr) + occ_hubb(m,m,isp,isp,iapr) * real( V_hubb(m,m,isp,isp,iapr), db )
        end do
      end do
    endif

! Becareful to the signs
    Energ_self(iapr) = Energ_self_KS(iapr) + E_KS_core(iapr) + E_exc(iapr) - E_elec(iapr) - E_Vxc(iapr) + E_hubbard(iapr)
    E_kin(iapr) = Energ_self(iapr) - U(iapr) - E_exc(iapr) 

  end do boucle_ia

  Energ_self_KS_agr = 0._db
  E_exc_agr = 0._db
  E_elec_agr = 0._db
  E_KS_core_agr = 0._db
  E_Vxc_agr = 0._db
  En_cluster = 0._db
  E_hubbard_agr = 0._db
  E_kin_agr = 0._db
  U_agr = 0._db

  do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      n = nb_eq(iapr)
    else
      if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
      if( .not. Proto_calculated(iapr) ) cycle
      if( Sym_2D ) then
        n = nb_eq_2D(iapr)
      else
        n = ngreq(iapr)
      endif
      if( Doping .and. iapr == ipr_dop ) n = n - 1
    endif
    Energ_self_KS_agr = Energ_self_KS_agr + n * Energ_self_KS(iapr)
    E_exc_agr = E_exc_agr + n * E_exc(iapr)
    E_elec_agr = E_elec_agr + n * E_elec(iapr)
    E_Vxc_agr = E_Vxc_agr + n * E_Vxc(iapr)
    E_hubbard_agr = E_hubbard_agr + n * E_hubbard(iapr)
    En_cluster = En_cluster + n * Energ_self(iapr)
    E_kin_agr = E_kin_agr + n * E_kin(iapr)
    E_KS_core_agr = E_KS_core_agr + n * E_KS_core(iapr)
    U_agr = U_agr + n * U(iapr)
  end do

  if( icheck > 1 ) then
    write(3,100)
    if( Full_atom ) then
      Atom_kind = ' ia'
    else
      Atom_kind = 'ipr'
    endif
    write(3,500) Atom_kind
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        ipr = iapr
        if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
        if( .not. Proto_calculated(ipr) ) cycle
     endif
      it = itypepr( ipr )
      write(3,510) iapr, numat(it), Energ_self(iapr) * rydb, Energ_self_KS(iapr) * rydb, E_KS_core(iapr) * rydb, &
         E_elec(iapr) * rydb, E_exc(iapr) * rydb, E_Vxc(iapr) * rydb, E_hubbard(iapr) * rydb, E_kin(iapr) * rydb, &
         U(iapr) * rydb
    end do
    write(3,520) En_cluster * rydb, Energ_self_KS_agr * rydb, E_KS_core_agr * rydb, &
         E_elec_agr * rydb, E_exc_agr * rydb, E_Vxc_agr * rydb, E_hubbard_agr * rydb, E_kin_agr * rydb, U_agr * rydb
  end if

  return

  100 format(/' ---- Energ_DFT ---',100('-'))
  500 format(/1x,a3,'   Z    Energ_atom      E_KS_val     E_KS_core        E_elec        E_exc        E_Vxc', &
    '        E_hubbard       E_kin    E_elec + E_Z')
  510 format(2i4,9f14.3)
  520 format(/' Total =',9f14.3)
end

!***********************************************************************

! Preparation of next iteration

subroutine prep_next_iter(Bulk_step,chargat_self,chargat_self_s,Convergence,Delta_En_conv,Delta_energ,Delta_energ_s, &
               Doping,En_cluster,En_cluster_s,Energ_self,Energ_self_s,Fermi,Fermi_first,Full_atom,Hubb, &
               Hubbard,i_self,iaprotoi,icheck,ipr_dop,itypepr,Lmaxat,m_hubb,mpirank,n_atom_0_self,n_atom_ind_self, &
               n_atom_proto,n_atom_proto_uc,n_devide,natome,nb_eq,nb_eq_2D,ngreq,nlm_pot,nrm_self,nself,nspin,nspinp, &
               ntype,numat,Occ_hubb,Occ_hubb_i,p_self,p_self_max,p_self0,pop_orb_val,Proto_calculated,rho_self,rho_self_s, &
               Rmtsd,Sym_2D,V_Hubb,V_Hubb_s,V_hubbard) 

  use declarations
  implicit none

  integer:: i_self, iapr, icheck, ipr, ipr_dop, it, L, L_hubbard, L_level_val, m_hubb, mpirank, n, n_atom_0_self, &
    n_atom_ind_self, n_atom_proto, n_atom_proto_uc, natome, n_devide, nlm_pot, nrm_self, nself, nspin, nspinp, ntype, Z
  integer, save:: Mod_p
  
  integer, dimension(0:n_atom_proto):: itypepr, Lmaxat, nb_eq_2D, ngreq
  integer, dimension(natome):: iaprotoi, nb_eq
  integer, dimension(0:ntype):: numat

  character(len=3):: Atom_kind

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb, V_hubb_s

  logical:: Bulk_step, Convergence, Doping, Fermi, Fermi_first, Full_atom, Hubbard, Sym_2D
  logical, dimension(0:n_atom_proto):: Proto_calculated
  logical, dimension(0:ntype):: Hubb

  real(kind=db):: Delta_En_conv, Delta_energ, Delta_energ_s, Delta_lim, En_cluster, En_cluster_s, &
    p_self, p_self_max, p_self_s, p_self0, pds
  real(kind=db), save:: Delta_energ_min
  real(kind=db), dimension(0:n_atom_proto):: Rmtsd
  real(kind=db), dimension(0:ntype):: V_hubbard
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: Energ_self, Energ_self_s
  real(kind=db), dimension(0:nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self):: rho_self, rho_self_s
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self,nspin):: chargat_self, chargat_self_s, pop_orb_val
  real(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: occ_hubb, occ_hubb_i

  if( .not. Fermi .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      if( icheck > 0 .or. ipr > 3 ) write(ipr,110)
    end do
    stop
  endif

! Interpolation between 2 last iterations
  if( i_self == 1 ) then

    Delta_energ = 1000000._db
    Delta_energ_min = Delta_energ
    p_self = p_self0 / 2
    En_cluster_s = En_cluster
    Mod_p = 0
    
  else

    if( icheck > 1 ) write(3,140)

! Test convergence: on the energy and on the charge of the absorbing atom
! To to before interpolation

    Delta_lim = 0._db
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        n = nb_eq(iapr)
        Z = numat( itypepr( iaprotoi(iapr) ) )
      else
        if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
        if( .not. Proto_calculated(iapr) ) cycle
        if( Sym_2D ) then
          n = nb_eq_2D(iapr)
        else
          n = ngreq(iapr)
        endif
        if( Doping .and. iapr == ipr_dop ) n = n - 1
        Z = numat( itypepr( iapr ) )
      endif
      if( Z <= 30 ) then
        pds = 1._db
      else
      ! had hoc formula to have beyond Z = 30 a limit increasing as Z^2
        pds = ( Z / 30._db )**2
      endif
      Delta_lim = Delta_lim + n * pds
    end do
    Delta_lim = Delta_lim * Delta_En_conv

    Delta_energ = 0._db
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        n = nb_eq(iapr)
      else
        if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
        if( .not. Proto_calculated(iapr) ) cycle
        if( Sym_2D ) then
          n = nb_eq_2D(iapr)
        else
          n = ngreq(iapr)
        endif
        if( Doping .and. iapr == ipr_dop ) n = n - 1
      endif
      Delta_energ = Delta_energ + n * abs( Energ_self(iapr) - Energ_self_s(iapr) )
    end do

    if( Delta_energ < Delta_lim ) Convergence = .true.

    if( Convergence ) then

      if( icheck > 0 ) then
        if( icheck > 1 ) write(3,*)
        write(3,150) En_cluster*rydb, Delta_energ*rydb, '<', Delta_lim*rydb, p_self
      endif

      if( mpirank == 0 ) write(6,155) i_self, En_cluster*rydb, Delta_energ*rydb, '<', Delta_lim*rydb, p_self
      
      if( Fermi_first ) then
        Convergence = .false.
        Fermi_first = .false.
        if( icheck > 0 ) write(3,160) n_devide
        if( mpirank == 0 ) write(6,160) n_devide
        i_self = min( i_self, nself - 1 )
      endif
        
    else

      if( icheck > 0 ) then
        if( icheck > 1 ) write(3,*)
        write(3,162) En_cluster*rydb, Delta_energ*rydb, '>', Delta_lim*rydb, p_self, Mod_p
      endif
      if( mpirank == 0 ) write(6,165) i_self, En_cluster*rydb, Delta_energ*rydb, '>', Delta_lim*rydb, p_self, Mod_p

      if( i_self == nself .and. mpirank == 0 ) then
        if( icheck > 0 ) write(3,'(/A)') ' Calculation has not converged !'
        write(6,'(/A)') ' Calculation has not converged !'
      else
        if( i_self == 1 .or. ( En_cluster - En_cluster_s ) > Delta_lim / 10 ) Mod_p = Mod_p + 1  
! When calculation is diverging, or with too much beating, one decreases the weight.
        p_self_s = p_self
        if( ( En_cluster - En_cluster_s  < - eps10 ) .or. Delta_energ > 3._db * abs( En_cluster - En_cluster_s ) ) then
          if( p_self0 > 0.39_db .or. i_self > 50 ) then
            p_self = max( p_self / 2, p_self0 / 8 )
          else
            p_self = max( p_self / 2, p_self0 / 4 )
          endif
          Mod_p = 0
        elseif( Mod_p >= 4 .or. ( Mod_p >= 3 .and. i_self <= 50 ) ) then
          p_self = min( 2 * p_self, 4 * p_self0 )
          p_self = min( p_self, p_self_max )
          Mod_p = 0
        endif
        if( abs( p_self_s - p_self ) > eps10 ) then
          if( icheck > 0 ) write(3,170) p_self, Mod_p
          if( mpirank == 0 ) write(6,170) p_self, Mod_p
        endif
      end if
    endif
  endif

  if( p_self > eps10 ) then
! In the other case, no self-consistency, chargat_self is just calculated for the fermi evaluation and not kept for potential calculation in xanes step 
    rho_self(:,:,:,:) = p_self * rho_self(:,:,:,:) + ( 1 - p_self ) * rho_self_s(:,:,:,:)
    chargat_self(:,:) = p_self * chargat_self(:,:) + ( 1 - p_self ) * chargat_self_s(:,:)
  endif

  if( Hubbard .and. i_self > 1 ) V_hubb(:,:,:,:,:) =  p_self * V_hubb(:,:,:,:,:) + ( 1 - p_self ) * V_hubb_s(:,:,:,:,:)

! One keeps the values of the current iteration, to use them in the next one
  if( .not. convergence .and. i_self /= nself ) then
    Delta_energ_s = Delta_energ
    En_cluster_s = En_cluster
    Energ_self_s(:) = Energ_self(:)
    chargat_self_s(:,:) = chargat_self(:,:)
    rho_self_s(:,:,:,:) = rho_self(:,:,:,:)
    if( Hubbard ) V_hubb_s(:,:,:,:,:) = V_hubb(:,:,:,:,:)
    Delta_energ_min = min( Delta_energ, Delta_energ_min )
  endif

  if( icheck == 1 .and. ( Convergence .or. i_self == nself ) ) then
    if( Full_atom ) then
      Atom_kind = ' ia'
    else
      Atom_kind = 'ipr'
    endif
    if( nspin == 1 ) then
      write(3,300) Atom_kind
    else
      write(3,305) Atom_kind
    endif
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
        ipr = iapr
        if( .not. Proto_calculated(ipr) ) cycle
      endif
      Z = numat( itypepr(ipr) )
      L = min( l_level_val(Z), lmaxat(ipr) )
      if( nspin == 2 ) then
        write(3,310) iapr, Z, Energ_self(iapr)*rydb, sum(chargat_self(iapr,:)), &
             chargat_self(iapr,2) - chargat_self(iapr,1), sum(pop_orb_val(iapr,:)), L, Rmtsd(ipr)*bohr
      else
        write(3,320) iapr, Z, Energ_self(iapr)*rydb, sum(chargat_self(iapr,:)), sum(pop_orb_val(iapr,:)), &
            L, Rmtsd(ipr)*bohr
      endif
    end do

    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        if( Bulk_step .and. iapr <= n_atom_proto_uc ) cycle  
        ipr = iapr
        if( .not. Proto_calculated(ipr) ) cycle
      endif
      it = itypepr(ipr)
      if( .not. Hubb(it) ) cycle
      Z = numat( itypepr(ipr) )
      L = L_hubbard( Z )
      call Write_occ_hubb(Atom_kind,iapr,L,m_hubb,n_atom_0_self,n_atom_ind_self,nspinp,occ_hubb,occ_hubb_i,Z)
      call Write_V_hubb(Atom_kind,iapr,it,L,m_hubb,n_atom_0_self,n_atom_ind_self,nspinp,ntype,V_hubb,V_hubbard,Z)
    end do

  endif

  return
  110 format(/' The Fermi level was not reached ! ')
  120 format(/' Cycle',i4,', Total Cluster energy =',f14.3,' eV')
  130 format(/' Cycle',i4,', Total Unit cell energy =',f14.3,' eV')
  140 format(/'----- prep_next_iter ', 100('-'))
  150 format(12x,' Total energy =',f13.3,' eV, Delta energy =',f11.3,' ',a1,f8.3,' eV,  Weight =',f8.5,'. Converged !')
  155 format(/' Cycle',i4,': Total energy =',f13.3,' eV, Delta energy =',f11.3,' ',a1,f8.3,' eV,  Weight =',f8.5,'. Converged !')
  160 format(/' New Fermi optimization with imaginary energy divided by',i2)
  162 format(12x,' Total energy =',f13.3,' eV, Delta energy =',f11.3,' ',a1,f8.3,' eV,  Weight =',f8.5,', Mod_p =',i2)
  165 format(/' Cycle',i4,': Total energy =',f13.3,' eV, Delta energy =',f11.3,' ',a1,f8.3,' eV,  Weight =',f8.5,', Mod_p =',i2)
  170 format(79x,' Next weight =',f8.5,', Mod_p =',i2)
  300 format(/6x,a3,'   Z       Energy       Charge  pop_orb_val(L)   L    Radius')
  305 format(/6x,a3,'   Z       Energy       Charge       up-dn  pop_orb_val(L)   L    Radius')
  310 format(5x,2i4,f14.3,3f12.3,i8,f10.5)
  320 format(5x,2i4,f14.3,2f12.3,i8,f10.5)
end

!***********************************************************************

! Transformation base reelle vers base complexe

subroutine Trans_rc(L,V)

  use declarations
  implicit none

  integer:: L, m1, m2

  complex(kind=db), dimension(-L:L,-L:L):: V, W
  complex(kind=db):: pm1_1, pm1_2, pm2_1, pm2_2

  real(kind=db):: r2

  r2 = 1 / sqrt(2._db)

  do m1 = -L,L
    if( m1 < 0 ) then
      pm1_1 = cmplx( 0._db, -r2 * (-1)**m1, db )
      pm1_2 = cmplx( r2 * (-1)**m1, 0._db, db )
    elseif( m1 == 0 ) then
      pm1_1 = cmplx( 1._db, 0._db, db )
      pm1_2 = cmplx( 0._db, 0._db, db )
    else
      pm1_1 = cmplx( r2, 0._db, db )
      pm1_2 = cmplx( 0._db, r2, db )
    endif
    do m2 = -L,L
      if( m2 < 0 ) then
        pm2_1 = cmplx( 0._db, -r2 * (-1)**m2, db )
        pm2_2 = cmplx( r2 * (-1)**m2, 0._db, db )
      elseif( m2 == 0 ) then
        pm2_1 = cmplx( 1._db, 0._db, db )
        pm2_2 = cmplx( 0._db, 0._db, db )
      else
        pm2_1 = cmplx( r2, 0._db, db )
        pm2_2 = cmplx( 0._db, r2, db )
      endif

      W(m1,m2) = Conjg( pm1_1 ) * pm2_1 * V(m1,m2) + Conjg( pm1_1 ) * pm2_2 * V(m1,-m2) &
               + Conjg( pm1_2 ) * pm2_1 * V(-m1,m2) + Conjg( pm1_2 ) * pm2_2 * V(-m1,-m2)

    end do
  end do

  V(:,:) = W(:,:)

end

!***********************************************************************

! Writing of the DOS

subroutine Write_stdens(Abs_exc,Cal_xanes,Density_comp,Energ,Harm_cubic,i_self,iapr,ie_computer,la, &
           Int_statedens_l,lla_state,lla2_state,mpinodes,n_atom_0,n_atom_ind,nomfich_s,nspin,nspinp, &
           Open_file,Statedens,statedens_l,Ylm_comp)

  use declarations
  implicit none

  integer:: iapr, i_self, ie_computer, index, isp, L, la, lla_state, lla2_state, lm, long, length, m, mpinodes, n_atom_0, &
         n_atom_ind, nspin, nspinp

  character(len=1), dimension(0:4):: Orb_L
  character(len=9), dimension(16):: Orb_Lm, Orb_Lm_D, Orb_Lm_K, Orb_Lm_Z
  character(len=13):: mot13
  character(len=13), dimension((lla2_state+2*lla_state)*nspinp+2):: Col_name
  character(len=Length_name) nomfich_s, nomficht

  logical:: Abs_exc, Cal_xanes, Density_comp, Harm_cubic, Open_file, Spin_out, Ylm_comp

  real(kind=db):: Energ, Int_rho, Rho
  real(kind=db), dimension(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1):: Statedens
  real(kind=db), dimension(0:lla_state,nspinp,n_atom_0:n_atom_ind):: Int_statedens_l
  real(kind=db), dimension(0:lla_state,nspinp):: statedens_l

  data Orb_L/'s','p','d','f','g'/
  data Orb_Lm_K/'s','px','py','pz','dx2-y2','dz2','dyz','dxz','dxy','fxyz','fx3','fy3','fz3','fx(y2-z2)','fy(z2-x2)','fz(x2-y2)'/
  data Orb_Lm_Z/'s','py','pz','px','dxy','dyz','dz2','dxz','dx2-y2','(3,-3)','fxyz','(3,-1)','fz3','(3,1)','fz(x2-y2)','(3,3)'/
  data Orb_Lm_D/'(0,0)','(1,-1)','(1,0)','(1,1)','(2,-2)','(2,-1)','(2,0)','(2,1)','(2,2)','(3,-3)','(3,-2)','(3,-1)','(3,0)', &
                '(3,1)','(3,2)','(3,3)'/

  if( Cal_xanes ) then
    if( Density_comp ) then
      Orb_Lm(:) = Orb_Lm_D(:)
    elseif( Harm_cubic ) then
      Orb_Lm(:) = Orb_Lm_K(:)
    else
      Orb_Lm(:) = Orb_Lm_Z(:)
    endif
  else
    if( Ylm_comp ) then
      Orb_Lm(:) = Orb_Lm_D(:)
    else
      Orb_Lm(:) = Orb_Lm_Z(:)
    endif
  endif

  Spin_out = nspin == 2 .or. ( nspinp == 2 .and. ( Density_comp .and. Cal_xanes ) .or. .not. Cal_xanes )
    
  nomficht = nomfich_s
  long = len_trim(nomfich_s)
  nomficht(long+1:long+3) = '_sd'

  if( Abs_exc ) then
    call ad_number(0,nomficht,Length_name)
  elseif( Cal_xanes ) then
    call ad_number(iapr,nomficht,Length_name)
  else
    call ad_number(i_self,nomficht,Length_name)
  endif
  long = len_trim(nomficht)
  nomficht(long+1:long+4) = '.txt'

  if( Open_file ) then

    open(4, file = nomficht )

    index = 0
    lm = 0    
    do L = 0,la
      if( .not. ( Spin_out .and. L == 0 ) ) then
        do m = -L,L
          mot13 = ' '
          lm = lm + 1
          index = index + 1
          mot13 = Orb_lm(lm)
          if( Spin_out ) then
            length = len_trim(mot13)
            mot13(length+1:length+3) = '_up'
            call Center_word(mot13,13) 
            Col_name(index) = mot13
    
            index = index + 1
            mot13 = adjustl(mot13)
            mot13(length+2:length+3) = 'dn'
            call Center_word(mot13,13) 
            Col_name(index) = mot13
          else
            call Center_word(mot13,13) 
            Col_name(index) = mot13
          endif
        end do
      else
        lm = lm + 1
      endif
      if( L /= 0 .or. Spin_out ) then
        mot13 = Orb_L(L)
        if( Spin_out ) mot13(2:4) = '_up'
        index = index + 1
        call Center_word(mot13,13) 
        Col_name(index) = mot13
      endif
      mot13 = 'Int( )'
      mot13(5:5) = Orb_L(L)
      if( Spin_out ) mot13(7:9) = '_up'
      index = index + 1
      call Center_word(mot13,13) 
      Col_name(index) = mot13
      if( Spin_out ) then
        mot13 = Orb_L(L)
        mot13(2:4) = '_dn'
        index = index + 1
        call Center_word(mot13,13) 
        Col_name(index) = mot13
        mot13 = 'Int( )'
        mot13(5:5) = Orb_L(L)
        mot13(7:9) = '_dn'
        index = index + 1
        call Center_word(mot13,13) 
        Col_name(index) = mot13
      endif
    end do
    mot13 = 'total'
    index = index + 1
    call Center_word(mot13,13) 
    Col_name(index) = mot13

    mot13 = 'Int(total)'
    index = index + 1
    call Center_word(mot13,13) 
    Col_name(index) = mot13
   
    write(4,110) Col_name(1:index)

  else

    open(4, file = nomficht, position='append')

  endif

  Rho = sum( Statedens_l(0:la,:) )
  Int_Rho = sum( Int_statedens_l(0:la,:,iapr) )

  if( Spin_out ) then
     write(4,120) Energ*rydb, ( Statedens_l(0,isp), Int_statedens_l(0,isp,iapr), isp = 1,nspinp ), &
      ((( Statedens(L**2+L+1+m,isp,L**2+L+1+m,isp,iapr,ie_computer), isp = 1,nspinp), m = -L,L ), &
        ( Statedens_l(L,isp), Int_statedens_l(L,isp,iapr), isp = 1,nspinp ), L = 1,la ), Rho, Int_rho
  else   
     write(4,120) Energ*rydb, sum( Statedens_l(0,:) ), sum( Int_statedens_l(0,:,iapr) ), &
      (( 2*Statedens(L**2+L+1+m,1,L**2+L+1+m,1,iapr,ie_computer), m = -L,L ), &
         sum( Statedens_l(L,:) ), sum( Int_statedens_l(L,:,iapr) ), L = 1,la ), Rho, Int_rho
  endif

  close(4)

  return
  110 format(4x,'Energy',320a13)
  120 format(f10.4,1p,320e13.5)
end
