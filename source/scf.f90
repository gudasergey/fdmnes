! Routines of the FDMNES package

!***********************************************************************

! Sousprogramme elaborant la grille en energie pour le calcul du niveau
! de Fermi

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
  120 format(/'   energie    eimag      en eV')
  130 format(4f9.3)
end

!***********************************************************************

! Calcul de la charge du petit agregat

subroutine chg_agr(chargat,chargat_init,ch_c,chg_cluster,chg_open_val,Doping,Full_atom,iaprotoi,ipr_dop,iprabs, &
                 ispin_maj,itabs,icheck,itypepr,mpirank, natome,n_atom_0_self,n_atom_ind_self, &
                 n_atom_proto,nb_eq,ngreq,nrato,nrm, nrm_self,nspin,ntype,numat,pop_open_val, &
                 psi_open_val,rato,rho_chg,rho_coeur, rhoato_init,rmtsd,SCF_mag_fix,SCF_mag_free)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer:: ipr_dop, itabs, Sum_Z
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, ngreq
  integer, dimension(natome):: iaprotoi, nb_eq
  integer, dimension(n_atom_0_self:n_atom_ind_self):: ispin_maj

  logical:: Doping, Full_atom, SCF_mag_fix, SCF_mag_free

  real(kind=db), dimension(nspin):: chg_cluster
  real(kind=db), dimension(0:n_atom_proto):: chargat
  real(kind=db), dimension(0:nrm,0:ntype):: rato, rho_coeur
  real(kind=db), dimension(0:nrm_self,nspin, n_atom_0_self:n_atom_ind_self):: rho_chg, rhoato_init
  real(kind=db), dimension(0:n_atom_proto):: rmtsd
  real(kind=db), dimension(0:nrm)::  r, rh
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: ch_v, ch_c, chargat_sup
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self,nspin):: chargat_init
  real(kind=db), dimension(nrm,2):: psi_open_val
  real(kind=db), dimension(2):: chg_open_val, pop_open_val

  charge_init = 0._db
  chg_cluster(:) = 0._db
  chg_sup = 0._db
  chg_coeur = 0._db
  chargat_init(:,:) = 0._db
  chargat_sup(:) = 0._db

  if( icheck > 0 ) write(3,110)

  sum_Z = 0
  SCF_mag_fix = .false.

  do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
      n = nb_eq(iapr)
    else
      ipr = iapr
      n = ngreq(ipr)
      if( Doping .and. iapr == ipr_dop ) n = n - 1
    endif
    it = itypepr(ipr)
    nr = nrato(it)
    rayint = rmtsd(ipr)
    r(0:nrm) = rato(0:nrm,it)

! Calcul de la charge totale de l'agregat: rhoato est la vraie densite
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
    chg_cluster(1) = chg_cluster(1) + chargat_init(iapr,isp) * n
    if( nspin == 2 ) then
      isp = 3 - ispin_maj(iapr)
      chg_cluster(2) = chg_cluster(2) + chargat_init(iapr,isp) * n
    endif

! Calcul de la charge totale de l'agregat: rho_chg est la densite
! venant des atomes exterieurs au petit agregat.
    do ispin = 1, nspin
      rh(0:nr) = rho_chg(0:nr,ispin,iapr) * r(0:nr)**2
      res = quatre_pi * f_integr3(r,rh,0,nrm,rayint)
      chargat_sup(iapr) = chargat_sup(iapr) + res
    end do
    chg_sup = chg_sup + chargat_sup(iapr) * n

! Calcul de la charge des orbitales de coeur: rho_coeur est la vraie
! densite
    rh(0:nr) = rho_coeur(0:nr,it) * r(0:nr)**2
    ch_c(iapr) = quatre_pi * f_integr3(r,rh,0,nrm,rayint)

    chg_coeur = chg_coeur + ch_c(iapr) * n
    Sum_Z = Sum_Z + n * numat(it)
    charge_init = charge_init + n * chargat(ipr)
  end do

! Calcul de la charge de valence:
  chg = sum( chg_cluster(:) ) - chg_coeur

  if( icheck > 0 ) then
     write(3,120) chg
     write(3,130) chg_coeur
     if( .not. SCF_Mag_fix ) then
       write(3,140) sum( chg_cluster(:) )
     else
       write(3,142) sum( chg_cluster(:) ), chg_cluster(:)
     endif
     write(3,145) chg_sup
     write(3,150) Sum_Z
     write(3,160) charge_init
     write(3,170) Sum_Z - charge_init - sum( chg_cluster(:) )
  end if

  do iprint = 3,6,3
    if( iprint == 3 .and. icheck == 0 ) cycle
    if( iprint == 6 .and. mpirank /= 0 ) cycle
    if( nspin == 1 ) then
      write(iprint,300)
    else
      write(iprint,305)
    endif
    do iapr = n_atom_0_self,n_atom_ind_self
      ch_v(iapr) = sum(chargat_init(iapr,:)) - ch_c(iapr)
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        ipr = iapr
      endif
      it = itypepr(ipr)
      if( nspin == 1 ) then
        write(iprint,310) iapr, numat(it), ch_v(iapr), ch_c(iapr), sum(chargat_init(iapr,:)), chargat_sup(iapr), &
                   numat(it)- sum(chargat_init(iapr,:))
      else
        write(iprint,310) iapr, numat(it), ch_v(iapr), ch_c(iapr), sum(chargat_init(iapr,:)), &
                   chargat_init(iapr,1) - chargat_init(iapr,2), chargat_sup(iapr), numat(it)- sum(chargat_init(iapr,:))
      endif
    end do
  end do

  Sum_Z = Sum_Z - charge_init

  if( icheck > 1 ) then
    write(3,350)
    write(3,351) sum(ch_v(:))
    write(3,352) sum(ch_c(:))
    write(3,353) sum(chargat_init(:,:))
  endif

  if( icheck > 2 ) then
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        ipr = iapr
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
! psi_level_val etait normalise a l'unite:
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
  300 format(/'  ia   Z     ch_val    ch_core   ch_total     ch_out   Atom charge')
  305 format(/'  ia   Z     ch_val    ch_core   ch_total   ch_up-dn     ch_out     Charge')
  310 format(i4,i4,6f11.3)
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

! Calcul du point du depart en energie, pour la determination du niveau
! de Fermi

subroutine En_dep(E_coeur_s,E_start,E_starta,Full_atom,iaprotoi,icheck,itypepr,lcoeur,n_atom_0,n_atom_0_self, &
                n_atom_ind,n_atom_ind_self,n_atom_proto,natome,ncoeur,nenerg_coh,nlm_pot,nrato,nrm,numat,nspin, &
                ntype,Pas_SCF,psi_coeur,relativiste,rato,Rmtg,Rmtsd,V_intmax,Vcato,Vxcato,workf)

  use declarations
  implicit none

  integer:: i, iapr, icheck, ipr, ir, it, lmax_pot_loc, n_atom_0, n_atom_0_self, &
     n_atom_ind, n_atom_ind_self, n_atom_proto, natome, nenerg_coh, nlm_pot, nlm_pot_loc, nr, nrm, nspin, ntype

  integer, dimension(2,0:ntype):: lcoeur, ncoeur
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(natome):: iaprotoi
  integer, dimension(0:n_atom_proto):: itypepr

  logical:: Full_atom, Relativiste, Ylm_comp

  real(kind=db):: E_marge, E_max_Fermi, E_start, Pas_SCF, psiHpsi, V_intmax, Vxc, Workf

  real(kind=db), dimension(1):: V0bdm
  real(kind=db), dimension(0:n_atom_proto):: Rmtg, Rmtsd
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(nrm):: r, psi
  real(kind=db), dimension(nrm,1,1):: pot
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: E_coeur_s, E_starta
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vxcato
  real(kind=db), dimension(0:nrm,2,0:ntype):: psi_coeur
  real(kind=db), dimension(2,n_atom_0_self:n_atom_ind_self):: E_psi

! On selectione l'orbitale de coeur de plus grande energie et celle de
! valence de plus basse energie (valable pour les niveaux atomiques
! ordones selon la regle d'occupation. Il est possible qu'ils soient
! inverses dans le cas d'un agregat)
! ces energies sont stoquees en E_psi

! l'indice 1 correspond a la derniere orbitale de coeur
! l'indice 2 correspond a la premiere orbitale de valence

  if( icheck > 0 ) write (3,110)

! Potentiel ici garde spherique
  Ylm_comp = .false.
  lmax_pot_loc = 0
  nlm_pot_loc = 1

  E_marge = 10._db / rydb
  E_start = 10000._db

  boucle_agr: do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
    else
      ipr = iapr
    endif
    it = itypepr(ipr)
    r(1:nrm) = rato(1:nrm,it)
    nr = nrato(it)

    do i = 1,2
      psi(1:nrm) = psi_coeur(1:nrm,i,it)

! On fait une moyenne sur les deux spins, si besoin :
     do ir = 1,nrm
       Vxc = sum( Vxcato(ir,1,:,iapr) ) / nspin
       pot(ir,1,1) = Vcato(ir,1,iapr) + Vxc
     end do
     V0bdm(1) = 0._db
     E_psi(i,iapr) =  psiHpsi(.false.,icheck,1,lcoeur(i,it),lmax_pot_loc,0,ncoeur(i,it),nlm_pot_loc,nr,nrm,1,1, &
            1,1,numat(it),V0bdm,pot,psi,r,Relativiste,Rmtg(ipr),Rmtsd(ipr),.false.,V_intmax,Ylm_comp)
    end do

    if( E_psi(2,iapr) - E_psi(1,iapr) > 2*E_marge ) then
      E_starta(iapr) = E_psi(2,iapr) - E_marge
    else
      E_starta(iapr) = 0.5_db * ( E_psi(1,iapr) + E_psi(2,iapr) )
    endif

    E_start = min( E_start, E_starta(iapr) )

    E_starta(iapr) = E_starta(iapr) - eps10

  end do boucle_agr

! les potentiels sont references en fonction du niveau du vide alors que
! la gamme d'energie commence dessous (workf)
  E_start = E_start + Workf
  E_starta(:) = E_starta(:) + Workf
  
  E_coeur_s(:) = E_psi(1,:)

  E_max_Fermi = 30._db / rydb

! Evaluation du nombre de points pour la grille en energie
  nenerg_coh = nint( ( E_max_Fermi - E_start ) / Pas_SCF ) + 1

  if( icheck > 0 ) then
    write(3,120)
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        ipr = iapr
      endif
      it = itypepr(ipr)
      write(3,130) iapr, numat(it), ( ncoeur(i,it), lcoeur(i,it), E_psi(i,iapr)*rydb, i = 1,2), E_starta(iapr)*rydb
    end do
  end if

  if( icheck > 2 ) then
    do it = 0,ntype
      write(3,135) numat(it), (ncoeur(i,it), lcoeur(i,it), i = 1,2)
      do ir = 1,nrato(it)
        write(3,137) ( psi_coeur(ir,i,it), i = 1,2 )
      end do
    end do
  end if

  if( icheck > 0 ) write (3,150) E_start * rydb, - Workf *rydb

  return
  110 format(/' ---- En_dep --------',100('-'))
  120 format(/'  ia     Z    n  l     E_core  n  l     E_val    E_starta')
  130 format(i4,3x,i3,2x,2(i3,i3,f11.3),f11.3)
  135 format(/'  Z = ',i3,3x,'n_coeur = ',i1,3x,'l_coeur = ',i1,3x, 'n_val = ',i1,3x,'l_val = ',i1// &
            4x,'psi_coeur    ','psi_val    '/)
  137 format(1p,2e14.6)
  150 format(/' Starting energy = ',f8.3,' eV,    zero at',f7.2,' eV, E_core, E_val zero at infinity.')
end

!***********************************************************************

! Sousroutine qui corrige l'energie sortie du calcul DFT
! les energies sont en rydb

subroutine Energ_DFT(Doping,En_cluster,Energ_self,E_coeur,excato,Full_atom,Hubb,iaprotoi,icheck,ipr_dop,itypepr, &
           m_hubb,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nb_eq,ngreq,nlm_pot,nrm,nrm_self, &
           nrato,nspin,nspinp,ntype,numat,rato,rho_self,rmtsd,V_hubb,Vcato,Vxcato)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: nlm_pot, nspin, nspinp

  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(natome):: iaprotoi, nb_eq
  integer, dimension(0:n_atom_proto):: itypepr, ngreq

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb

  logical:: Doping, Full_atom
  logical, dimension(0:ntype):: Hubb

  real(kind=db), dimension(0:nrm_self,nlm_pot,nspin, n_atom_0_self:n_atom_ind_self):: rho_self
  real(kind=db), dimension(nrm,n_atom_0:n_atom_ind):: excato
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: E_coeur, Energ_self, Energ_self_KS, &
                     E_exc, E_coul, E_hubbard, E_Vxc
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vxcato
  real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: Vhartree
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind)::Vcato
  real(kind=db), dimension(0:nrm):: fct1, fct2, fct3, r, r2
  real(kind=db), dimension(0:n_atom_proto):: rmtsd

! Enragr: la somme des valeurs propres pour les orbitales Kohn Sham occupees

  if( nrm /= nrm_self ) return  ! car rho_self n'est pas calcule

  Energ_self_KS(:) = Energ_self(:)

  boucle_ia: do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
    else
      ipr = iapr
    endif
    it = itypepr(ipr)
    nr = nrato(it)
    r(:) = rato(:,it)
    r2(1:nr) = r(1:nr)**2
    r2(0) = 0._db

    fct1(:) = 0._db; fct2(:) = 0._db; fct3(:) = 0._db

! le potentiel dont on a besoin pour la correction est le vrai Hartree,
!  c'est a dire il faut enlever la contribution des noyaux

    Vhartree(1:nr,iapr) = Vcato(1:nr,1,iapr) + 2 * numat(it)/r(1:nr)

    do ir = 1,nr
      fac = sum( rho_self(ir,1,1:nspin,iapr) ) * r2(ir)
      fct1(ir) = excato(ir,iapr) * fac
      fct2(ir) = Vhartree(ir,iapr) * fac
    end do

    E_exc(iapr) = quatre_pi * f_integr3(r,fct1,0,nrm,Rmtsd(ipr))
    E_coul(iapr) = quatre_pi * f_integr3(r,fct2,0,nrm,Rmtsd(ipr))
    E_coul(iapr) = 0.5_db * E_coul(iapr)

    E_Vxc(iapr) = 0._db
    do ispin = 1, nspin
      do ir = 1,nr
        fct3(ir) = Vxcato(ir,1,ispin,iapr) * rho_self(ir,1,ispin,iapr) * r2(ir)
      end do
      E_Vxc(iapr) = E_Vxc(iapr) + quatre_pi * f_integr3(r,fct3,0,nrm,Rmtsd(ipr))
    end do

    E_hubbard(iapr) = 0._db
    fac = 1._db / ( 3 - nspin )
    if( Hubb(it) ) then
      do isp = 1,nspinp
        do m = -m_hubb,m_hubb
          E_hubbard(iapr) = E_hubbard(iapr) + real( V_hubb(m,m,isp,isp,iapr), db )
        end do
      end do
    endif

! E_coeur: energie des orbitales KS correspondant aux etats de coeur, pour chaque atome;

! ATTENTION aux signes
    Energ_self(iapr) = Energ_self_KS(iapr) + E_coeur(iapr) + E_exc(iapr) - E_coul(iapr) - E_Vxc(iapr) + E_hubbard(iapr)
  end do boucle_ia

  Energ_self_KS_agr = 0._db
  Delta_E_coeur_agr = 0._db
  E_exc_agr = 0._db
  E_coul_agr = 0._db
  E_Vxc_agr = 0._db
  En_cluster = 0._db
  E_hubbard_agr = 0._db

  do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      n = nb_eq(iapr)
    else
      n = ngreq(iapr)
      if( Doping .and. iapr == ipr_dop ) n = n - 1
    endif
    Energ_self_KS_agr = Energ_self_KS_agr + n * Energ_self_KS(iapr)
    Delta_E_coeur_agr = Delta_E_coeur_agr + n * E_coeur(iapr)
    E_exc_agr = E_exc_agr + n * E_exc(iapr)
    E_coul_agr = E_coul_agr + n * E_coul(iapr)
    E_Vxc_agr = E_Vxc_agr + n * E_Vxc(iapr)
    E_hubbard_agr = E_hubbard_agr + n * E_hubbard(iapr)
    En_cluster = En_cluster + n * Energ_self(iapr)
  end do

  if( icheck > 0 ) then
    write(3,100)
    if( Full_atom ) then
      write(3,110) En_cluster * rydb
    else
      write(3,120) En_cluster * rydb
    endif
    write(3,500)
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        ipr = iapr
      endif
      it = itypepr( ipr )
      write(3,510) iapr, numat(it), Energ_self(iapr) * rydb, Energ_self_KS(iapr) * rydb, E_coeur(iapr) * rydb, &
         E_exc(iapr) * rydb, E_coul(iapr) * rydb, E_Vxc(iapr) * rydb, E_hubbard(iapr) * rydb
    end do
    write(3,520) En_cluster * rydb, Energ_self_KS_agr * rydb, Delta_E_coeur_agr * rydb, &
         E_exc_agr * rydb, E_coul_agr * rydb, E_Vxc_agr * rydb, E_hubbard_agr * rydb
  end if

  return

  100 format(/' ---- En_DFT ------',100('-'))
  110 format(/' Cluster energy: ', f15.3,' eV')
  120 format(/' Unit cell energy: ', f15.3,' eV')
  500 format(/'   ia   Z    Energ_atom      Energ_KS    Delta_E_coeur     E_exc        E_coul        E_Vxc', &
    '        E_hubbard')
  510 format(i5,i4,7f14.3)
  520 format(/'  Total:',7f14.3)
end

!********************************************************************************

! Sousprogramme qui calcule l'energie de la derniere orbitales de coeur Kohn Sham;
! On considere que pendant les iterations la fonction d'onde de cette orbitale ne
! se modifie pas, par contre le potentiel et donc son energie bougent. On fait
! l'approximation que ce deplacement est le meme pour tous les niveaux de coeur

! le resultat sera utilise pour le calcul de l'energie de l'agregat

subroutine eps_coeur(ch_coeur,E_coeur,E_coeur_s,Full_atom,iaprotoi,icheck,itypepr,lcoeur,n_atom_0,n_atom_0_self, &
              n_atom_ind,n_atom_ind_self,n_atom_proto,natome,ncoeur,nlm_pot,nrato,nrm,nspin,ntype,numat,psi_coeur, &
              rato,Relativiste,Rmtg,Rmtsd,V_intmax,Vcato,Vxcato)

  use declarations
  implicit none

  integer iapr, icheck, ipr, ir, it, lmax_pot_loc, natome, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, &
          n_atom_proto, nlm_pot, nlm_pot_loc, nr, nrm, nspin, ntype

  logical:: Full_atom, Relativiste, Ylm_comp

  integer,dimension(2,0:ntype):: lcoeur, ncoeur
  integer,dimension(natome):: iaprotoi
  integer,dimension(0:n_atom_proto):: itypepr
  integer,dimension(0:ntype):: nrato, numat

  real(kind=db):: psiHpsi, res, V_intmax, Vxc
  real(kind=db), dimension(nspin):: V0bd
  real(kind=db), dimension(0:n_atom_proto):: Rmtg, Rmtsd
  real(kind=db),dimension(n_atom_0_self:n_atom_ind_self):: ch_coeur
  real(kind=db),dimension(n_atom_0_self:n_atom_ind_self):: E_coeur, E_coeur_s
  real(kind=db),dimension(0:nrm,0:ntype):: rato
  real(kind=db),dimension(nrm):: r, psi
  real(kind=db),dimension(nrm,1,1):: pot
  real(kind=db),dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db),dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato
  real(kind=db),dimension(0:nrm,2,0:ntype):: psi_coeur

! Potentiel ici garde spherique
  Ylm_comp = .false.
  lmax_pot_loc = 0
  nlm_pot_loc = 1
  V0bd(:) = 0._db

  E_coeur(:) = 0._db

  do iapr = n_atom_0_self,n_atom_ind_self
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
    else
      ipr = iapr
    endif
    it = itypepr(ipr)
    nr = nrato(it)
    r(1:nrm) = rato(1:nrm,it)
! On fait une moyenne sur les deux spins, si besoin :
    do ir = 1,nrm
      Vxc = sum( Vxcato(ir,1,1:nspin,iapr) ) / nspin
      pot(ir,1,1) = Vcato(ir,1,iapr) + Vxc
      psi(ir) = psi_coeur(ir,1,it)
    end do
    res =  psiHpsi(.false.,icheck,1,lcoeur(1,it),lmax_pot_loc,0,ncoeur(1,it),nlm_pot_loc,nr,nrm,1,1, &
              1,1,numat(it),V0bd,pot,psi,r,Relativiste,Rmtg(ipr),Rmtsd(ipr),.false.,V_intmax,Ylm_comp)
    res = res - E_coeur_s(iapr)
    E_coeur(iapr) = res * nint( ch_coeur(iapr) )
  end do

  return
end

!***********************************************************************

! Preparation de l'iteration suivante

subroutine prep_next_iter(chargat_self,chargat_self_s,Convergence,Delta_En_conv,Delta_energ,Delta_energ_s, &
               Delta_energ_t,Doping,En_cluster,En_cluster_s,Energ_self,Energ_self_s,Fermi,Fermi_first,Full_atom, &
               Hubbard,i_self,icheck,ipr_dop,m_hubb,mpirank,n_atom_0_self,n_atom_ind_self,n_atom_proto,n_devide,natome,natomeq, &
               nb_eq,ngreq,nlm_pot,nrm_self,nself,nspin,nspinp,p_self,p_self_max,p_self0,rho_self,rho_self_s, &
               V_Hubb,V_Hubb_s) 

  use declarations
  implicit none

  integer:: i_self, iapr, icheck, ipr, ipr_dop, m_hubb, mpirank, n, n_atom_0_self, n_atom_ind_self, n_atom_proto, &
    natome, natomeq, n_devide, nlm_pot, nrm_self, nself, nspin, nspinp
  integer, save:: Mod_p
  
  integer, dimension(0:n_atom_proto):: ngreq
  integer, dimension(natome):: nb_eq

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb, V_hubb_s

  logical:: Convergence, Doping, Fermi, Fermi_first, Full_atom, Hubbard

  real(kind=db):: Delta_En_conv, Delta_energ, Delta_energ_s, Delta_energ_t, Delta_lim, En_cluster, En_cluster_s, &
    p_self, p_self_max, p_self0
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: Energ_self, Energ_self_s
  real(kind=db), dimension(0:nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self):: rho_self, rho_self_s
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self,nspin):: chargat_self, chargat_self_s

  if( .not. Fermi .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      if( icheck > 0 .or. ipr > 3 ) write(ipr,110)
    end do
    stop
  endif

  if( mpirank == 0 ) write(6,120) i_self, En_cluster * rydb

! Interpolation d'une iteration de la boucle coherente a l'autre
  if( i_self == 1 ) then

    Delta_energ = 1000000._db
    p_self = p_self0 / 2
    En_cluster_s = En_cluster
    Mod_p = 0
    
  else

    if( icheck > 0 ) write(3,125)

! Test convergence: sur l'energie et sur la charge de l'atome central
! a faire avant l'interpolation
    Delta_energ = 0._db
    if( Full_atom ) then
      Delta_lim = Delta_En_conv * natomeq
    else
      Delta_lim = Delta_En_conv * sum( ngreq(1:n_atom_proto) )
    endif

    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        n = nb_eq(iapr)
      else
        n = ngreq(iapr)
        if( Doping .and. iapr == ipr_dop ) n = n - 1
      endif
      Delta_energ = Delta_energ + n * abs( Energ_self(iapr) - Energ_self_s(iapr) )
    end do

    if( Delta_energ < Delta_lim ) Convergence = .true.

    if( Convergence ) then

      if( icheck > 0 ) then
        write(3,*)
        write(3,130) Delta_energ*rydb, '<', Delta_lim*rydb, p_self
      endif
      if( mpirank == 0 ) write(6,130) Delta_energ*rydb, '<', Delta_lim*rydb, p_self
      
      if( Fermi_first ) then
        Convergence = .false.
        Fermi_first = .false.
        if( icheck > 0 ) write(3,135) n_devide
        if( mpirank == 0 ) write(6,135) n_devide
        i_self = min( i_self, nself - 1 )
      endif
        
    else

      if( icheck > 0 ) then
        write(3,*)
        write(3,130) Delta_energ*rydb, '>', Delta_lim*rydb, p_self
      endif
      if( mpirank == 0 ) write(6,130) Delta_energ*rydb, '>', Delta_lim*rydb, p_self

    endif
    
    if( .not. Convergence ) then
      if( i_self == nself .and. mpirank == 0 ) then
        if( icheck > 0 ) write(3,'(/A)') ' Calculation has not converged !'
        write(6,'(/A)') ' Calculation has not converged !'
      else
        if( i_self == 1 .or. ( En_cluster - En_cluster_s ) > Delta_lim / 10 ) Mod_p = Mod_p + 1  
! Si calcul divergent, ou avec trop de batement, on diminue le poids.
        if( ( En_cluster - En_cluster_s  < - eps10 ) .or. Delta_energ > 3._db * abs( En_cluster - En_cluster_s ) ) then
          if( p_self0 > 0.39_db .or. i_self > 50 ) then
            p_self = max( p_self / 2, p_self0 / 8 )
          else
            p_self = max( p_self / 2, p_self0 / 4 )
          endif
          Mod_p = 0
        elseif( Mod_p >= 4 .or. ( Mod_p >= 3 .and. i_self <= 50 ) ) then
          if( p_self0 > 0.39_db .or. i_self > 50 ) then
            p_self = min( 2 * p_self, 2 * p_self0 )
          else
            p_self = min( 2 * p_self, 4 * p_self0 )
          endif
          p_self = min( p_self, p_self_max )
          Mod_p = 0
        endif
      end if
    endif
  endif

  rho_self(:,:,:,:) = p_self * rho_self(:,:,:,:) + ( 1 - p_self ) * rho_self_s(:,:,:,:)
  chargat_self(:,:) = p_self * chargat_self(:,:) + ( 1 - p_self ) * chargat_self_s(:,:)

  if( Hubbard .and. i_self > 1 ) V_hubb(:,:,:,:,:) =  p_self * V_hubb(:,:,:,:,:) + ( 1 - p_self ) * V_hubb_s(:,:,:,:,:)

! On stoque les valeurs de l'iteration courrante, pour les injecter dans la suivante
  if( .not. convergence .and. i_self /= nself ) then
    Delta_energ_t = Delta_energ_s
    Delta_energ_s = Delta_energ
    En_cluster_s = En_cluster
    Energ_self_s(:) = Energ_self(:)
    chargat_self_s(:,:) = chargat_self(:,:)
    rho_self_s(:,:,:,:) = rho_self(:,:,:,:)
    if( Hubbard ) V_hubb_s(:,:,:,:,:) = V_hubb(:,:,:,:,:)
  endif

  return
  110 format(/' The Fermi level was not reached ! ')
  120 format(/' Cycle',i4,', Total Cluster energy =',f14.3,' eV')
  125 format(/'----- prep_next_iter ', 70('-'))
  130 format(12x,'Delta_energ =',f11.3, ' eV ',a1,' Delta =',f8.3,' eV,  Weight =',f8.5)
  135 format(/' New Fermi optimization with imaginary energy divided by',i2)
end

!***********************************************************************

! Calcul et ecriture de la densite d'etat et calcul du niveau de Fermi

! chargat_self = charge atomique a l'iteration courrante
! chargat_self_s = charge atomique a l'iteration precedente
! ch_ia = nombre total d'electrons (y compris de coeur) a
!         l'iteration courrante, pour l'atome ia

subroutine Cal_State(chg_cluster,chg_open_val,Cal_xanes,chargat_self,Density,Doping,drho_self,E_cut,E_Open_val, &
                E_Open_val_exc,E_starta,Energ,E_Fermi,Enragr,Energ_self,Fermi,Full_atom,hubb,Hubb_diag,iaabsi,iaprotoi, &
                i_self,icheck,ie,ie_computer,Int_statedens,ipr_dop,ispin_maj,itypei,itypepr,lamstdens, &
                Open_val,Open_val_exc,lla_state,lla2_state,lmaxat,m_hubb,mpinodes,n_atom_0,n_atom_0_self, &
                n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nb_eq,nenerg,ngreq,nlm_pot,nomfich_s,nonexc_g,nrato, &
                nrm,nrm_self,nspin,nspinp,ntype,numat,occ_hubb,occ_hubb_i,pop_orb_val,rato,rho_self,rho_self_t,rmtsd,SCF_elecabs, &
                SCF_mag_fix,Self_nonexc,State_all,Statedens,Statedens_i,V_hubb,V_hubbard)

  use declarations
  implicit none

  integer:: i, i_self, ia, iaabsi, iapr, icheck, ie, ie_computer, imax, ipr, ipr_dop, iprabs, iprint, ir, isp, isp1, isp2, it, l, &
    l_hubbard, l_level_val, la, lamstdens, ll, lla_state, lla2_state, lh, lm, lm0, lm1, lm2, lma, m, m_hubb,mpinodes, m1, m2, n, &
    n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, n_atom_proto, natome, nenerg, nlm_pot, nr, nrm, nrm_self, nspin, &
    nspinp, ntype, Numat_tot, Z

  character(len=132) nomfich_s

  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, la_ipr, ll_ipr, lmaxat, ngreq
  integer, dimension(natome):: iaprotoi, itypei, nb_eq
  integer, dimension(n_atom_0_self:n_atom_ind_self):: ispin_maj

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb

  logical, save:: Fermi_gen, Fermi_maj, Fermi_min
  logical::Abs_exc, Cal_xanes, Density, Doping, Fermi, First_loop, Full_atom, nonexc_g, Open_file, Open_val, Open_val_exc, &
    SCF_elecabs, SCF_mag_fix, Self_nonexc, State_all
  logical, dimension(0:n_atom_proto):: proto_done
  logical, dimension(0:ntype):: Hubb
  logical, dimension(n_atom_0_self:n_atom_ind_self):: Hubb_diag

  real(kind=db):: ch, ch_s, Charge, Charge_maj, Charge_min, Charge_val_abs, chg_lim, D_Energ, de, Dens, ds, E_cut, E_Fermi, &
    E_Open_val, E_Open_val_exc, En_f, Enragr, Enragr_s, f_integr3, poids, Pop

  real(kind=db), save:: Charge_maj_s, Charge_min_s, Charge_s, Charge_val_abs_e_s, Charge_val_abs_s, E_Fermi_maj, E_Fermi_min

  real(kind=db), dimension(2) :: chg_open_val
  real(kind=db), dimension(nspin):: chg_cluster, chg_ref
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(0:ntype):: V_hubbard
  real(kind=db), dimension(0:n_atom_proto):: rmtsd
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self):: E_starta
  real(kind=db), dimension(nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self,0:mpinodes-1):: drho_self
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

  if( ( icheck > 1 .and. .not. Cal_xanes ) .or. icheck > 2 )  then
    write(3,110)
    write(3,270) Energ(ie)*rydb
  endif

! Stockage des tableaux de l'iteration precedante, utilises pour
! l'interpolation:
  if( .not. cal_xanes ) then
    if( ie == 1 ) then
      Fermi_gen = .false.
      Fermi_maj = .false.
      Fermi_min = .false.
      E_Fermi_maj = E_Fermi
      E_Fermi_min = E_Fermi
      Charge_s = 0._db; Charge_maj_s = 0._db; Charge_min_s = 0._db
      pop_orb_val(:,:) = 0._db
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
    Enragr_s = Enragr
    Energ_self_s(:) = Energ_self(:)
    pop_orb_val_s(:,:) = pop_orb_val(:,:)
  endif

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
!            la = min(2,lmaxat(ipr))  ! Version Hiroshi
        ll = min(2,lmaxat(ipr))
      elseif( Z > 18 .and. Z < 55 ) then
        la = min(2,lmaxat(ipr))
        ll = min(3,lmaxat(ipr))
      else
        la = min(3,lmaxat(ipr))
        ll = min(4,lmaxat(ipr))
      endif
    end if
    la_ipr(ipr) = la    ! nombre maximal d'harmoniques pour l'ecriture de la densite d'etat
    ll_ipr(ipr) = ll    ! nombre maximal d'harmoniques pour le calcul  de la densite d'etat
  end do

  boucle_ia: do ia = 1,natome

    ipr = iaprotoi(ia)
    la = la_ipr(ipr)
    ll = ll_ipr(ipr)

    if( ia == iaabsi ) iprabs = ipr

    if( Full_atom ) then
      iapr = ia
    else
      iapr = ipr
    endif

    it = itypei(ia)
    Z = numat(it)
    if( proto_done(ipr) .and. .not. Full_atom ) cycle boucle_ia
    if( Cal_xanes ) then
      if( .not. ( ( ia == iaabsi .and. Density ) .or. State_all ) ) cycle
    else
      if( Energ(ie) < E_starta(iapr) ) cycle
    endif
    Open_file = .false.
    if( ie == 1 ) then
      Open_file = .true.
    elseif( .not. Cal_xanes ) then
      if( Energ(ie-1) < E_starta(iapr) ) Open_file = .true.
    endif
     
    if( ie == 1 ) then
      First_loop = .true.
    elseif( .not. Cal_xanes ) then
      if( Energ(ie-1) < E_starta(iapr) ) then
        First_loop = .true.
      else
        First_loop = .false.
      endif
    else
      First_loop = .false.
    endif

! Calcul de la densite d'etat integree jusqu'au rayon Rmstd de l'atome

    do ir = 1,nrato(it)
      if( rato(ir,it) > Rmtsd(ipr) + eps10 ) exit
    end do
    nr = ir

    lma = (ll + 1)**2

    Statedens_l(:,:) = 0._db
    do l = 0,la
      lm1 = l**2 + 1
      lm2 = ( l + 1 )**2
      do isp = 1,nspinp
        do lm = lm1,lm2
          Statedens_l(l,isp) = Statedens_l(l,isp) + Statedens(lm,isp,lm,isp,iapr,ie_computer)
        end do
      end do
    end do
    if( nspinp == 1 ) Statedens_l(0:la,1) = 2 * Statedens_l(0:la,1)

! Integrale de la densite d'etats
    do isp = 1,nspinp
      do lm = 1,lma
        Int_Statedens(lm,isp,iapr) = Int_Statedens(lm,isp,iapr) + de * Statedens(lm,isp,lm,isp,iapr,ie_computer)
      end do
    end do

    do l = 0,ll
      lm1 = l**2 + 1
      lm2 = ( l + 1 )**2
      do isp = 1,nspinp
        Int_Statedens_l(l,isp,iapr) = sum(Int_Statedens(lm1:lm2,isp,iapr))
      end do
      if( nspinp == 1 ) Int_Statedens_l(l,1,iapr) = 2 * Int_Statedens_l(l,1,iapr)
    end do

    Int_statedens_t(:,iapr) = 0._db
    do isp = 1,nspinp
      Int_Statedens_t(isp,iapr) = Int_Statedens_t(isp,iapr) + sum( Int_Statedens_l(0:la,isp,iapr) )
    end do

    if( icheck > 1 ) then
      write(3,180) iapr
      lm = 0
      do l = 0,ll
        do m = -l,l
          lm = lm + 1
          do isp = 1,nspinp
            write(3,190) l, m, isp, Statedens(lm,isp,lm,isp,iapr,ie_computer), Int_Statedens(lm,isp,iapr)
          end do
        end do
      end do
      write(3,200)
      do l = 0,la
        write(3,210) l, Int_Statedens_l(l,1:nspinp,iapr)
      end do
      write(3,220) Int_Statedens_t(1:nspinp,iapr)
    endif

    if(    ( Cal_xanes .and. ( ( ia == iaabsi .and. Density ) .or. state_all ) ) &
      .or. ( .not. Cal_xanes .and. ia == iaabsi .and. Density .and. icheck > 1 ) ) then
      Abs_exc = Full_atom .and. ( ia == iaabsi ) .and. Cal_xanes .and. .not. nonexc_g
      call write_stdens(Abs_exc,Cal_xanes,Energ(ie),i_self,iapr,ie_computer,la,Int_Statedens,Int_Statedens_t, &
           Int_Statedens_l,lla_state,lla2_state,mpinodes,n_atom_0,n_atom_ind,nomfich_s,nspinp,Open_file,Statedens, &
           Statedens_l)
    endif

! Matrice d'occupation pour Hubbard
    if( .not. Cal_xanes .and. hubb(it) ) then
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

    if( .not. cal_xanes ) then

! Energie de l'agregat (ou de la cellule elementaire)
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
        if( Full_atom ) then
          Enragr = Enragr + nb_eq(iapr) * D_Energ
        else
          Enragr = Enragr + ngreq(iapr) * D_Energ
        endif
      end do

      l = min( l_level_val(Z), lmaxat(ipr) )
      do isp = 1,nspin
        if( SCF_mag_fix .and. ( ( isp == ispin_maj(iapr) .and. Fermi_maj ) .or. &
            ( isp /= ispin_maj(iapr) .and. Fermi_min ) ) ) cycle
        if( nspinp /= nspin ) then
          pop_orb_val(iapr,isp) = sum( Int_statedens_l(l,:,iapr) )
        else
          pop_orb_val(iapr,isp) = Int_statedens_l(l,isp,iapr)
        endif
      end do

      r(:) = rato(:,it)
      r2(:) = rato(:,it)**2

      if( First_loop ) then
        do isp = 1, nspin
          rh(0:nr) = rho_self(0:nr,1,isp,iapr) * r2(0:nr)
          chargat_self(iapr,isp) = Real( Z,db ) / nspin - quatre_pi * f_integr3(r,rh,0,nrm,Rmtsd(ipr))
          chargat_self_s(iapr,isp) = chargat_self(iapr,isp)
        end do

        if( SCF_mag_fix ) rho_self_t(:,:,:) = rho_self(:,:,1,:) + rho_self(:,:,nspin,:)

        if( icheck > 1 ) write(3,150) iapr, chargat_self(iapr,:)
        if( icheck > 2 ) then
          write(3,160) iapr, Z
          do ir = 1,nr
            write(3,170) rato(ir,it)*bohr, ( quatre_pi * r2(ir) * rho_self(ir,lm,:,iapr), lm = 1,nlm_pot )
          end do
        end if
      endif

      if( nspinp == 1 ) then
        ds = 2 * de
      else
        ds = de
      endif
      do isp = 1,nspin
        if( SCF_mag_fix .and. .not. Fermi_gen ) rho_self_t(1:nr,:,iapr) = rho_self_t(1:nr,:,iapr) &
                      + ds * drho_self(1:nr,:,isp,iapr,ie_computer)
        if( ( isp == ispin_maj(iapr) .and. Fermi_maj ) .or. ( isp /= ispin_maj(iapr) .and. Fermi_min ) ) cycle
        rho_self(1:nr,:,isp,iapr) = rho_self(1:nr,:,isp,iapr) + ds * drho_self(1:nr,:,isp,iapr,ie_computer)
      end do

! Calcul de la charge
      do isp = 1,nspin
        rh(0:nr) = rho_self(0:nr,1,isp,iapr) * r2(0:nr)
        ch_ia(iapr,isp) = quatre_pi * f_integr3(r,rh,0,nrm,Rmtsd(ipr))
        chargat_self(iapr,isp) = Real( Z, db ) / nspin - ch_ia(iapr,isp)
      end do
      if( SCF_mag_fix .and. .not. Fermi_gen ) then
        rh(0:nr) = rho_self_t(0:nr,1,iapr) * r2(0:nr)
        ch_ia_t(iapr) = quatre_pi * f_integr3(r,rh,0,nrm,Rmtsd(ipr))
      endif

      if( icheck > 1 ) write(3,230) iapr, ( isp, chargat_self(iapr,isp), Int_statedens_t(min(isp,nspin),iapr), isp = 1,nspin )
      if( icheck > 2 ) then
        write(3,160) iapr, Z
        do ir = 1,nr
          write(3,170) rato(ir,it)*bohr, ( quatre_pi * r2(ir) * rho_self(ir,lm,:,iapr), lm = 1,nlm_pot )
        end do
      end if
    endif

    proto_done(ipr) = .true.

  end do boucle_ia

  if( Cal_xanes ) return

! Interpolation:
! E_f = E_i*(ch_ref - ch_i-1)/(ch_i-ch_i-1) + E_i-1*(ch_i - ch_ref)/(ch_i-ch_i-1)
! ch_s, enragr = nombre total d'electrons et l'Energie de l'agregat a l'iteration precedente
! ch, enr = nombre total d'electrons et Energie courrantes
! chg_cluster = la charge qu'on va comparer avec la reference

! Evaluation du niveau de Fermi

  Z = numat( itypei(iaabsi) )
  l = l_level_val(Z)
  Charge = 0._db
  Charge_maj = 0._db
  Charge_min = 0._db
  chg_ref(:) = chg_cluster(:)
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
      n = ngreq(iapr)
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

  if( icheck > 1 .and. ie == 1 ) then
    if( nspin == 1 ) then
      write(3,260) Numat_tot, chg_ref(1)
    else
      write(3,265) Numat_tot, sum( chg_ref(:) ), chg_ref(:)
    endif
  endif

  if( icheck > 1 ) then
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom) then
        it = itypei(iapr)
      else
        it = itypepr(iapr)
      endif
      if( .not. hubb(it) ) cycle
      Z = numat(it)
      l = l_hubbard( Z )
      if( nspinp == 1 ) then
        write(3,266) iapr, Z, ( m2, m2 = -l,l ) 
      else
        write(3,267) iapr, Z, ( ( m2, isp2, m2 = -l,l ), isp2 = 1,nspinp )
      endif
      do isp1 = 1,nspinp
        do m1 = -l,l
          if( nspinp == 1 ) then
            write(3,268) m1, (( occ_hubb(m1,m2,isp1,isp2,iapr), occ_hubb_i(m1,m2,isp1,isp2,iapr), m2 = -l,l ), &
                                    isp2 = 1,nspinp )
          else
            write(3,269) m1, isp1, (( occ_hubb(m1,m2,isp1,isp2,iapr), occ_hubb_i(m1,m2,isp1,isp2,iapr), m2 = -l,l ), &
                                      isp2 = 1,nspinp )
          endif 
        end do
      end do
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

! On a atteind le niveau de Fermi
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
 ! On demand or when "calculation and SCF" excited and l = 2 or 3.
          if( Open_val_exc .and. scf_elecabs ) then
            E_cut = E_Open_val_exc
          else
            E_cut = E_Fermi
          endif
        case(3)
          E_Fermi_maj = En_f
        case(4)
          E_Fermi_min = En_f
      end select

      do iapr = n_atom_0_self,n_atom_ind_self
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
        if( i == 2 .and. SCF_mag_fix ) cycle
        occ_hubb(:,:,:,:,iapr) = occ_hubb(:,:,:,:,iapr) * poids + occ_hubb_s(:,:,:,:,iapr) * ( 1 - poids )
        occ_hubb_i(:,:,:,:,iapr) = occ_hubb_i(:,:,:,:,iapr) * poids + occ_hubb_i_s(:,:,:,:,iapr) * ( 1 - poids )
      end do

      if( .not. Fermi ) cycle

      Enragr = Enragr * poids + Enragr_s * ( 1 - poids )
      Energ_self(:) = Energ_self(:) * poids + Energ_self_s(:) * ( 1 - poids )

      if( Full_atom ) then
        Charge_val_abs = sum( pop_orb_val(iaabsi,:) )
      else
        Charge_val_abs = sum( pop_orb_val(iprabs,:) )
      endif

      do iprint = 3,6,3
        if( icheck == 0 .and. iprint == 3 ) cycle
        if( SCF_mag_fix ) then
          write(iprint,289) i_self, E_Fermi*rydb, E_Fermi_maj*rydb, E_Fermi_min*rydb, Enragr*rydb
        else
          write(iprint,290) i_self, E_Fermi*rydb, Enragr*rydb
        endif
        if( Open_val_exc .and. .not. self_nonexc ) then
          write(iprint,292) E_Open_val_exc*rydb
          write(iprint,293) E_Open_val*rydb
        elseif( Open_val ) then
          write(iprint,293) E_Open_val*rydb
        else
          write(iprint,294) Charge_val_abs
        endif
        if( nspin == 1 ) then
          write(iprint,300)
        else
          write(iprint,305)
        endif
        do iapr = n_atom_0_self,n_atom_ind_self
          if( Full_atom ) then
            ipr = iaprotoi(iapr)
          else
            ipr = iapr
          endif
          Z = numat( itypepr(ipr) )
          l = min( l_level_val(Z), lmaxat(ipr) )
          if( nspin == 2 ) then
            write(iprint,310) iapr, Z, Energ_self(iapr)*rydb, sum(chargat_self(iapr,:)), &
                 chargat_self(iapr,2) - chargat_self(iapr,1), sum(pop_orb_val(iapr,:)), l, Rmtsd(ipr)*bohr
          else
            write(iprint,320) iapr, Z, Energ_self(iapr)*rydb, sum(chargat_self(iapr,:)), sum(pop_orb_val(iapr,:)), &
                l, Rmtsd(ipr)*bohr
          endif
        end do
      end do
      if( icheck > 1 ) then
        do iapr = n_atom_0_self,n_atom_ind_self
          if( Full_atom ) then
            ipr = iaprotoi(iapr)
          else
            ipr = iapr
          endif
          it = itypepr(ipr)
          Z = numat( it )

          write(3,160) iapr, Z
          do ir = 1,nrato(it)
            if( rato(ir,it) > Rmtsd(ipr) + eps10 ) exit
            write(3,170) rato(ir,it)*bohr, ( quatre_pi * ( rato(ir,it)**2 ) * rho_self(ir,lm,:,iapr), lm = 1,nlm_pot )
          end do
        end do
      endif

! Matrice de Hubbard

      Hubb_diag(:) = .true.
      do iapr = n_atom_0_self,n_atom_ind_self

        if( Full_atom ) then
          ipr = iaprotoi(iapr)
        else
          ipr = iapr
        endif
        it = itypepr(ipr)
        Z = numat( it )
        if( .not. hubb(it) ) cycle
        l = l_hubbard( Z )

        Pop = 0._db
        do isp = 1,nspinp
          do m = -l,l
            Pop = Pop + occ_hubb(m,m,isp,isp,iapr)
          end do
        end do

        if( nspinp == 1 ) then
          occ_hubb(:,:,:,:,iapr) = 0.5_db * occ_hubb(:,:,:,:,iapr)
          occ_hubb_i(:,:,:,:,iapr) = 0.5_db * occ_hubb_i(:,:,:,:,iapr)
        endif
        if( icheck > 0 ) then
          if( nspinp == 1 ) then
            write(3,266) iapr, Z, ( m2, m2 = -l,l ) 
          else
            write(3,267) iapr, Z, ( ( m2, isp2, m2 = -l,l ), isp2 = 1,nspinp )
          endif
          do isp1 = 1,nspinp
            do m1 = -l,l
              if( nspinp == 1 ) then
                write(3,268) m1, (( occ_hubb(m1,m2,isp1,isp2,iapr), occ_hubb_i(m1,m2,isp1,isp2,iapr), m2 = -l,l ), &
                                        isp2 = 1,nspinp )
              else
                write(3,269) m1, isp1, (( occ_hubb(m1,m2,isp1,isp2,iapr), occ_hubb_i(m1,m2,isp1,isp2,iapr), m2 = -l,l ), &
                                          isp2 = 1,nspinp )
              endif 
            end do
          end do
          write(3,350) Pop
        endif

        occ_hubb(:,:,:,:,iapr) = - occ_hubb(:,:,:,:,iapr)
        do isp = 1,nspinp
          do m = -l,l
            occ_hubb(m,m,isp,isp,iapr) = occ_hubb(m,m,isp,isp,iapr) + 0.5_db
          end do
        end do

        V_hubb(:,:,:,:,iapr) = V_hubbard(it) * cmplx( occ_hubb(:,:,:,:,iapr), occ_hubb_i(:,:,:,:,iapr) )

        if( icheck > 0 ) then
          if( nspinp == 1 ) then
            write(3,360) iapr, Z, V_hubbard(it) * Rydb, ( m2, m2 = -l,l )
          else
            write(3,361) iapr, Z, V_hubbard(it) * Rydb, ( ( m2, isp2, m2 = -l,l ), isp2 = 1,nspin )
          endif
          do isp1 = 1,nspinp
            do m1 = -l,l
              if( nspinp == 1 ) then
                write(3,268) m1, (( V_hubb(m1,m2,isp1,isp2,iapr) * Rydb, m2 = -l,l ), isp2 = 1,nspinp )
              else
                write(3,269) m1, isp1, (( V_hubb(m1,m2,isp1,isp2,iapr) * Rydb, m2 = -l,l ), isp2 = 1,nspinp ) 
              endif
            end do
          end do
        endif

        boucle_m: do m1 = -l,l
          do isp1 = 1,nspinp
            do m2 = -l,l
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
  110 format(/' ---- Cal_state --------',100('-'))
  150 format(/'  Before integration:  ia = ',i3,'  charge_self = ', 2f10.5)
  160 format(/' ia =',i3,', Z =',i3/ '   Radius_(A) 4*pi*r2*Rho_self')
  170 format(1p,30e13.5)
  180 format(/'  l  m is  Density of states   Integral    ia =',i3)
  190 format(3i3,3f15.7)
  200 format(/'    l   sum_m(Integral)')
  210 format(i5,2f15.7)
  220 format(' Total =',f12.7,f15.7)
  230 format(/' Atome =',i3,2(', isp =',i2,', Charge_self =',f10.5, ', Int_state_t =',f10.5))
  260 format(/' Sum of atomic numbers =',i5,', Reference charge =', f10.3)
  265 format(/' Sum of atomic numbers =',i5,', Reference charge =', f10.3,', spin maj, spin min =',2f10.3)
  266 format(/' Hubbard occupation matrix for ia =',i3,', Z =',i3,/'  m',2(1x,7(10x,i3,10x)))
  267 format(/' Hubbard occupation matrix for ia =',i3,', Z =',i3,/' m1 s1',2(1x,7(9x,2i3,8x)))
  268 format(i3,7(1x,2f11.7))
  269 format(2i3,2(1x,7(1x,2f11.7)))
  270 format(15x,' Energy =',f10.3,' eV')
  289 format(/' Cycle',i4,',  Fermi Energy =',f8.3,' eV,  maj, min =',2f8.3,' eV', &
                      ',  Cluster Energy_KS =',f11.3,' eV')
  290 format(/' Cycle',i4,',  Fermi Energy =',f8.3,' eV,  Cluster Energy_KS =',f11.3,' eV')
  292 format(9x,'Level val excite =',f8.3,' eV')
  293 format(9x,'Level val absorb =',f8.3,' eV')
  294 format(9x,'Popul val absorb =',f8.3)
  300 format(/5x,'  ia   Z   Energy_KS      Charge  pop_orb_val(l)   l    Radius')
  305 format(/5x,'  ia   Z   Energy_KS      Charge       up-dn  pop_orb_val(l)   l    Radius')
  310 format(5x,2i4,4f12.3,i8,f10.5)
  320 format(5x,2i4,3f12.3,i8,f10.5)
  350 format(/' Total number of electron, summed over spin:',f9.5)
  360 format(/' Hubbard potential for ia =',i3,', Z =',i3,', Hubbard parameter U - J =',f5.2,' eV',/ &
              '  m',2(1x,7(10x,i3,10x)))
  361 format(/' Hubbard potential for ia =',i3,', Z =',i3,', Hubbard parameter U - J =',f5.2,' eV',/ &
              ' m1 s1',2(1x,7(9x,2i3,8x)))
end

!***********************************************************************

! Transformation base reelle vers base complexe

subroutine Trans_rc(l,V)

  use declarations
  implicit none

  integer:: l, m1, m2

  complex(kind=db), dimension(-l:l,-l:l):: V, W
  complex(kind=db):: pm1_1, pm1_2, pm2_1, pm2_2

  real(kind=db):: r2

  r2 = 1 / sqrt(2._db)

  do m1 = -l,l
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
    do m2 = -l,l
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

! Ecriture de la densite d'etat.

subroutine write_stdens(Abs_exc,Cal_xanes,Energ,i_self,iapr,ie_computer,la,Int_statedens,Int_statedens_t, &
           Int_statedens_l,lla_state,lla2_state,mpinodes,n_atom_0,n_atom_ind,nomfich_s,nspinp,Open_file,Statedens, &
           statedens_l)

  use declarations
  implicit none

  integer:: i, iapr, i_self, ie_computer, isp, l, la, lla_state, lla2_state, long, m, mpinodes, n_atom_0, &
         n_atom_ind, nspinp

  character(len=13):: mot13
  character(len=13), dimension(nspinp):: nomtIn
  character(len=13), dimension(0:lla_state,nspinp):: nomln, nomlIn
  character(len=13), dimension(0:lla_state,-lla_state:lla_state,nspinp):: nomn, nomIn
  character(len=132) nomfich_s, nomficht

  logical:: Abs_exc, Cal_xanes, Open_file

  real(kind=db):: Energ
  real(kind=db), dimension(lla2_state,nspinp,n_atom_0:n_atom_ind):: Int_statedens
  real(kind=db), dimension(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1):: Statedens
  real(kind=db), dimension(0:lla_state,nspinp,n_atom_0:n_atom_ind):: Int_statedens_l
  real(kind=db), dimension(0:lla_state,nspinp):: statedens_l
  real(kind=db), dimension(nspinp,n_atom_0:n_atom_ind):: Int_statedens_t

  nomficht = nomfich_s
  long = len_trim(nomfich_s)
  nomficht(long+1:long+3) = '_sd'

  if( Abs_exc ) then
    call ad_number(0,nomficht,132)
  elseif( Cal_xanes ) then
    call ad_number(iapr,nomficht,132)
  else
    call ad_number(i_self,nomficht,132)
  endif
  long = len_trim(nomficht)
  nomficht(long+1:long+4) = '.txt'

  if( Open_file ) then

    open(4, file = nomficht )

    if( nspinp == 1 ) then
      nomtIn(1) = '     Int_t   '
    else
      nomtIn(1) = '   Int_t(u)  '
      nomtIn(nspinp) = '   Int_t(d)  '
    endif
    do l = 0,la
      mot13 = 'n_l('
      mot13(5:5) = achar(l+48)
      if( nspinp == 1 ) then
        mot13(6:6) = ')'
        nomln(l,1) = '   ' // mot13(1:6)
        nomlIn(l,1) = '  Int' // mot13(1:6)
      else
        mot13(6:8) = ',u)'
        nomln(l,1) = '   ' // mot13(1:8)
        nomlIn(l,1) = '  Int' // mot13(1:8)
        mot13(6:8) = ',d)'
        nomln(l,nspinp) = '   ' // mot13(1:8)
        nomlIn(l,nspinp) = '  Int' // mot13(1:8)
      endif
      do m = -l,l
        mot13 = 'n('
        call ad_number(l,mot13,13)
        mot13(4:4) = ','
        call ad_number(m,mot13,13)
        i = len_trim(mot13)
        if( nspinp == 1 ) then
          mot13(i+1:i+1) = ')'
          nomn(l,m,1) = '    ' // mot13(1:7)
          nomIn(l,m,1) = '   Int' // mot13(1:7)
        else
          mot13(i+1:i+3) = ',u)'
          nomn(l,m,1) = '   ' // mot13(1:9)
          nomIn(l,m,1) = ' Int' // mot13(1:9)
          mot13(i+2:i+2) = 'd'
          nomn(l,m,nspinp) = '   ' // mot13(1:9)
          nomIn(l,m,nspinp) = ' Int' // mot13(1:9)
        endif
      end do
    end do

    write(4,110) (nomtIn(isp), isp = 1,nspinp), ( ( ( nomn(l,m,isp), nomIn(l,m,isp), isp = 1,nspinp), m = -l,l ), &
      ( nomln(l,isp), nomlIn(l,isp), isp = 1,nspinp ), l = 0,la )

  else

    open(4, file = nomficht, position='append')

  endif

  write(4,120) Energ*rydb,(Int_statedens_t(isp,iapr),isp = 1,nspinp), &
    (((Statedens(l**2+l+1+m,isp,l**2+l+1+m,isp,iapr,ie_computer), Int_statedens(l**2+l+1+m,isp,iapr), isp = 1,nspinp), &
    m = -l,l ),(Statedens_l(l,isp), Int_statedens_l(l,isp,iapr),isp = 1,nspinp ), l = 0,la )
  close(4)

  return
  110 format(4x,'Energy',320a13)
  120 format(f10.4,1p,320e13.5)
end
