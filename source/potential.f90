! FDMNES subroutines
! Calculation of the potential

! Potsup makes the superposition of the atomic potential and density

subroutine Potsup(alfpot,Axe_atom_gr,Bulk_atom_done,Cal_xanes,Chargat,chargat_init, &
            chargat_self,Delta_helm,dExc_ex_nex,drho_ex_nex,dv_ex_nex,dvc_ex_nex,Exc_abs_i,Excato,Full_atom,Helm_cos,Hybrid, &
            i_self,ia_eq_inv,ia_eq_inv_self,iaabs,iaproto,iaprotoi,iapot,icheck,igreq,igroup,iprabs_nonexc,ipr1,itab, &
            itypei,itypep,itypepr,lmax_pot,lvval,Magnetic,mpirank,n_atom_0,n_atom_0_self,n_atom_ind, &
            n_atom_ind_self,n_atom_proto,n_atom_proto_bulk,natome,natome_self,natomeq,natomeq_self,natomp,neqm, &
            ngroup_m,ngroup_nonsph, &
            nhybm,nlat,nlatm,nlm_pot,Nonexc,Nonsph,norbv,normrmt,npoint,npoint_ns,npsom,nrato,nrm,nrm_self,nspin, &
            ntype,numat,Overlap,pop_nonsph,popatm,popatv,pos,posi,posi_self,Proto_calculated,psival,Rato,rho,rho_chg, &
            rho_self,rhoato,rhoato_init,rhoit,rhons,Rmtg,Rmtimp,Rmtg0,Rmtsd,Rot_Atom_gr,Rot_int,rs,rsato,Rsort, &
            SCF,Self_nonexc,Sym_2D,V_abs_i,V_helm,V_intmax,Vc_abs_i,Vcato,Vcato_init,Vh,Vhns,Vsphere,Vxc,Vxcato, &
            V0bdcFimp,Width_helm,xyz)

  use declarations
  implicit none

  integer:: i_self, ia, iaabs, iapr, iapr0, iaprabs, iprabs_nonexc, iaprex, ipr, ipr1, ir, ispin, it, itab, &
    japr, lmax_pot, mpirank, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, n_atom_proto, n_atom_proto_bulk, &
    n_iapr, natome,natome_self, &
    natomeq, natomeq_self, natomp, neqm, ngroup_m, ngroup_nonsph, nhybm, nlatm, nlm_pot, normrmt, npoint, npoint_ns, &
    npsom, nr, nrm, nrm_self, nspin, ntype

  integer, dimension(30):: icheck
  integer, dimension(natomp):: iaproto, igroup, itypep
  integer, dimension(natome):: iaprotoi, itypei
  integer, dimension(natomeq):: ia_eq_inv
  integer, dimension(natomeq_self):: ia_eq_inv_self
  integer, dimension(0:ntype):: nlat, nrato, numat
  integer, dimension(0:ngroup_nonsph):: norbv
  integer, dimension(0:n_atom_proto):: iapot, itypepr
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(0:ntype,nlatm):: lvval

  complex(kind=db), dimension(nhybm,16,ngroup_nonsph) :: hybrid

  logical:: Bulk_atom_done, Cal_xanes, Do_init, Full_atom, Helm_cos, Magnetic, Nonexc, Nonsph, SCF, Self_nonexc, &
            Sym_2D, This_bulk_atom_done
  logical, dimension(0:n_atom_proto):: Proto_calculated

  real(kind=db):: alfpot, Delta_helm, Dist, f_integr3, Overlap, Orig_helm, Rayint, Rsort, V_helm, V_intmax, V0bdcFimp, &
                  Vsphere, Width_helm, z_bulk, z_surf
  real(kind=db), dimension(3):: V_surf
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self,nspin):: chargat_init, chargat_self
  real(kind=db), dimension(npoint):: rs, Vh
  real(kind=db), dimension(npoint_ns):: rhons, Vhns
  real(kind=db), dimension(0:ntype):: Rmtimp
  real(kind=db), dimension(0:n_atom_proto):: Chargat, rhonspr, Rmtg, Rmtg0, Rmtsd, Vcmft, Vhnspr
  real(kind=db), dimension(3,natome_self):: posi_self
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_gr
  real(kind=db), dimension(0:nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self):: rho_self
  real(kind=db), dimension(0:nrm_self,nspin,n_atom_0_self:n_atom_ind_self):: rho_chg, rhoato_init
  real(kind=db), dimension(nrm):: dExc_ex_nex, dvc_ex_nex, Exc_abs_i, Vc_abs_i
  real(kind=db), dimension(nrm,nspin):: drho_ex_nex, dv_ex_nex, V_abs_i
  real(kind=db), dimension(0:nrm_self,n_atom_0_self:n_atom_ind_self):: Vcato_init
  real(kind=db), dimension(npoint,nspin):: Vxc, rho
  real(kind=db), dimension(nrm):: exc
  real(kind=db), dimension(3,3):: Rot_int
  real(kind=db), dimension(nrm,n_atom_0:n_atom_ind):: Excato
  real(kind=db), dimension(nhybm,ngroup_nonsph) :: pop_nonsph
  real(kind=db), dimension(0:ntype,nlatm) :: popatv
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(3,3,ngroup_m):: Rot_Atom_gr
  real(kind=db), dimension(0:nrm):: rsato_e
  real(kind=db), dimension(0:nrm,nlm_pot):: Vcato_e
  real(kind=db), dimension(0:nrm,nspin):: rho_chg_e, rho_no_sup_e, rhoato_e, rhoato_init_e
  real(kind=db), dimension(0:nrm,nlm_pot,nspin):: vxcato_e
  real(kind=db), dimension(0:nrm,0:n_atom_proto):: vato
  real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: dvcato, drhoato, rsato
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nspin,0:n_atom_proto):: rhoigr
  real(kind=db), dimension(0:nrm,nspin,n_atom_0:n_atom_ind):: rhoato, rho_no_sup
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato
  real(kind=db), dimension(0:n_atom_proto,nspin):: rhomft, Vxcmft
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm
  real(kind=db), dimension(0:nrm,0:ntype):: Rato, rhoit
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(0:nrm):: drhoato_e, dvcato_e, r, rhr2, Vcato_init_e
  real(kind=db), dimension(n_atom_0:n_atom_ind):: ch

  drhoato_e(:)= 0._db

  do it = 0,ntype
    call potato(icheck(10),it,itypepr,n_atom_proto,nlat,nlatm,nrato,nrm,nspin,ntype,numat, &
        popatm,popatv,psival,Rato,rhoigr,rhoit,vato)
  end do

  if( i_self == 1 ) then
    dExc_ex_nex(:) = 0._db
    drho_ex_nex(:,:) = 0._db
    dvc_ex_nex(:) = 0._db
    dv_ex_nex(:,:) = 0._db
  endif
  rhonspr(:) = 0._db
  Vhnspr(:) = 0._db

  if( Nonsph ) call orbval(Hybrid,iaproto,iapot,icheck(11), &
        igreq,igroup,itypepr,lvval,mpirank,n_atom_proto,natomeq,natomp,neqm,ngroup_m,ngroup_nonsph,nhybm,nlat,nlatm,norbv, &
        npoint,npoint_ns,npsom,nrato,nrm,ntype,pop_nonsph,pos,psival,Rato,rhons,rhonspr,Rot_Atom_gr,Rot_int,Vhns,Vhnspr,xyz)

! Some preliminary works on Helmholtz potential
  Orig_helm = 0._db
  z_bulk = - 1000._db
  z_surf = - 1000._db
  if( abs( V_helm ) > eps10 ) then
    Orig_helm = -1000000._db
! V_surf is the direction perpendicular to the surface in the rotated cluster basis
    V_surf(1:3) = Rot_int(3,1:3)
    if( Sym_2D ) then
      do ia = 1,natomp
        Dist = sum( V_surf(1:3)*pos(1:3,ia) )
        Orig_helm = max( Orig_helm, Dist + Delta_helm )
      end do
    else
      do ia = 1,natomp
        Orig_helm = max( Orig_helm, sqrt( sum( pos(:,ia)**2 ) ) + Delta_helm )
      end do
    endif
    if( icheck(13) > 0 .and. i_self == 1 ) write(3,'(/a12,f10.5,a2)') ' Orig_helm =',Orig_helm*bohr,' A'
! Looking for the bulk atom, the closest from the suface
    do ia = 1,natomp
      ipr = iaproto( ia )
      Dist = sum( V_surf(1:3) * pos(1:3,ia) )
      if( ipr > n_atom_proto - n_atom_proto_bulk ) then
        z_bulk = max( z_bulk, Dist )
      else
        z_surf = max( z_surf, Dist )
      endif
    end do
  endif

! Calculation of the superposed potential of the atoms in the cluster
! It is calculated for all non equivalent atoms 1) of the ponctual group when "Full_atom"
!                                               2) of the space group when not "Full_atom"
! The absorbing atom is calculated in both excited and not excited states

  if( Cal_xanes ) then
    iapr0 = n_atom_0
    n_iapr = n_atom_ind
  else
    iapr0 = n_atom_0_self
    n_iapr = n_atom_ind_self
  endif

  if( Full_atom ) then
    if( ( Cal_xanes .and. .not. nonexc ) .or. ( .not. Cal_xanes .and. .not. self_nonexc ) ) then
      iaprabs = 0
      do ia = 1,n_atom_ind
        if( iaprotoi(ia) == 0 ) exit
      end do
      iaprex = ia  ! index of the absorbing atom in the cluster (before excitation)
    else
      do iaprabs = 1,n_atom_ind
        if( iaprotoi(iaprabs) == iprabs_nonexc ) exit
      end do
      iaprex = 0   ! when the absorbing atom is not excited, iaprex is the index of the excited one
    endif
  else
    iaprabs = iprabs_nonexc
    iaprex = 0
  endif

! At the first iteration, one calculates both the excited and non excited absorbing atom potential
  if( i_self == 1 ) iapr0 = 0   ! index for excited atom

  if( icheck(13) > 2 ) write(3,110) i_self, Full_atom, iapr0, n_iapr

  Do_init = i_self == 1 .and. .not. Cal_xanes

! Loop over all the atoms
  do iapr = iapr0,n_iapr

    if( Full_atom ) then
      if( iapr == 0 ) then
        if( ( Cal_xanes .and. nonexc ) .or. ( .not. Cal_xanes .and. self_nonexc ) ) then
          ipr = 0
        else
          ipr = iprabs_nonexc
        endif
      else
        ipr = iaprotoi( iapr )
      endif
    else
      ipr = iapr
    endif

    it = itypepr(ipr)
    nr = nrato(it)
    if( iapr > 0 ) then
      This_bulk_atom_done = Bulk_atom_done .and. ipr > n_atom_proto - n_atom_proto_bulk
    else
      This_bulk_atom_done = .false.
    endif

    if( ( i_self > 1 .and. iapr <= n_atom_ind_self ) .or. This_bulk_atom_done ) then
      if( Cal_xanes .and. Full_atom ) then
! In the XANES step calculation, it is possible to work with a ponctual group higher in symmetry than when making the SCF calculation
        do japr = 1,n_atom_ind_self
          if( sum( abs( posi(:,iapr) - posi_self(:,japr) ) ) < eps10 ) exit
        end do
      elseif( cal_xanes .and. iapr == 0 .and. Self_nonexc ) then
        japr = iprabs_nonexc
      else
        japr = iapr
      endif
      rho_chg_e(:,:) = rho_chg(:,:,japr)
      if( This_bulk_atom_done ) then
        Vcato_e(:,:) = Vcato(:,:,japr)
      else
        Vcato_init_e(:) = Vcato_init(:,japr)
      endif
      rhoato_init_e(:,:) = rhoato_init(:,:,japr)
    else
      rho_chg_e(:,:) = 0._db
      rhoato_init_e(:,:) = 0._db
    endif

! Calculation of the potential inside the atoms
    call Pot0muffin(alfpot,Cal_xanes,Chargat,chargat_init,chargat_self, &
        drho_ex_nex,drhoato_e,dvc_ex_nex,dvcato_e,Exc,Full_atom,Helm_cos,i_self,ia_eq_inv_self,iaproto,iapot,iapr, &
        icheck(13),ipr,iprabs_nonexc,itypep,itypepr,lmax_pot,Magnetic,n_atom_0,n_atom_0_self, &
        n_atom_ind_self,n_atom_proto,natome,natome_self,natomeq,natomeq_self,natomp,nlm_pot,nonexc,nrato,nrm,nrm_self, &
        nspin,ntype,numat,Orig_helm,pos,posi,Rato,rho_chg_e,rho_no_sup_e,rho_self,rhoato_e,rhoato_init_e,rhoigr, &
        rhonspr(ipr),Rmtsd(ipr),rsato_e,SCF,Self_nonexc,Sym_2D,This_bulk_atom_done,V_helm,V_surf,Vato,Vcato_e, &
        Vcato_init_e,Vhnspr(ipr),Vsphere,Vxcato_e,Width_helm,z_bulk,z_surf)

    if( iapr >= n_atom_0 ) then
      rsato(:,iapr) = rsato_e(:)
      Excato(1:nr,iapr) = Exc(1:nr)
      Vcato(:,:,iapr) = Vcato_e(:,:)
      Vxcato(:,:,:,iapr) = Vxcato_e(:,:,:)
      rhoato(:,:,iapr) = rhoato_e(:,:)
      rho_no_sup(:,:,iapr) = rho_no_sup_e(:,:)
      dvcato(:,iapr) = dvcato_e(:)
      drhoato(:,iapr) = drhoato_e(:)
    endif

    if( Do_init .and. iapr >= n_atom_0_self ) then
      Vcato_init(:,iapr) = Vcato_e(:,1)
      rhoato_init(:,:,iapr) = rhoato_e(:,:)
      rho_chg(:,:,iapr) = rho_chg_e(:,:)
    endif

    if( i_self == 1 ) then
! Calcul de la difference excite - non excite
      if( iapr == iaprex ) then
        drho_ex_nex(1:nrm,:) = drho_ex_nex(1:nrm,:) + rhoato_e(1:nrm,:)
        dvc_ex_nex(1:nrm) = dvc_ex_nex(1:nrm) + Vcato_e(1:nrm,1)
        do ispin = 1,nspin
          dv_ex_nex(1:nrm,ispin) = dv_ex_nex(1:nrm,ispin) + Vcato_e(1:nrm,1) + Vxcato_e(1:nrm,1,ispin)
        end do
        dExc_ex_nex(1:nrm) = dExc_ex_nex(1:nrm) + Exc(1:nrm)
      elseif( iapr == iaprabs ) then
        drho_ex_nex(1:nrm,:) = drho_ex_nex(1:nrm,:) - rhoato_e(1:nrm,:)
        dvc_ex_nex(1:nrm) = dvc_ex_nex(1:nrm) - Vcato_e(1:nrm,1)
         do ispin = 1,nspin
          dv_ex_nex(1:nrm,ispin) = dv_ex_nex(1:nrm,ispin) - Vcato_e(1:nrm,1) - Vxcato_e(1:nrm,1,ispin)
        end do
        dExc_ex_nex(1:nrm) = dExc_ex_nex(1:nrm) - Exc(1:nrm)
      endif
    endif

    if( ( ( .not. Cal_xanes .and. Self_nonexc ) .or. i_self == 1 ) .and. iapr == iaprabs ) then
      do ispin = 1,nspin
        V_abs_i(1:nrm,ispin) = Vcato_e(1:nrm,1) + Vxcato_e(1:nrm,1,ispin)
      end do
      Vc_abs_i(1:nrm) = Vcato_e(1:nrm,1)
      Exc_abs_i(1:nrm) = Exc(1:nrm)
    elseif( .not. Cal_xanes .and. .not. Self_nonexc .and. iapr == iaprex ) then
      do ispin = 1,nspin
        V_abs_i(1:nrm,ispin) = Vcato_e(1:nrm,1) + Vxcato_e(1:nrm,1,ispin) - dv_ex_nex(1:nrm,ispin)
      end do
      Vc_abs_i(1:nrm) = Vcato_e(1:nrm,1) - dvc_ex_nex(1:nrm)
      Exc_abs_i(1:nrm) = Exc(1:nrm) - dExc_ex_nex(1:nrm)
    endif

  end do ! end of loop on the atoms

  if( Full_atom ) call Force_pot_eq(Do_init,drhoato,dVcato,Excato,iaprotoi,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self, &
                        natome,nlm_pot,nrm,nrm_self,nspin,posi,rho_chg,rho_no_sup,rhoato,rhoato_init,rsato,Vcato,Vcato_init,Vxcato)

  if( i_self == 1 ) then
    it = itab
    if( icheck(13) > 2 ) then
      if( nspin == 1 ) then
        write(3,120) iaprabs
      else
        write(3,125) iaprabs
      endif
      do ir = 1,nrato(it)
        write(3,130) Rato(ir,it)*bohr, quatre_pi * Rato(ir,it)**2 * drho_ex_nex(ir,:), dv_ex_nex(ir,:)*Rydb
      end do
    endif
  endif

! Calculation of the muffin-tin atomic radius
  call Raymuf(Cal_xanes,Chargat,Full_atom,i_self,iapot,iaproto,iaprotoi,icheck(13),iprabs_nonexc, &
        ipr1,itab,itypei,itypep,itypepr,n_atom_0,n_atom_ind,n_atom_proto,natome,natomeq,natomp,nlm_pot, &
        normrmt,nrato,nrm,nspin,ntype,numat,Overlap,pos,Rato,rhoato,Rmtg,Rmtg0,Rmtimp, &
        Rmtsd,Rsort,V0bdcFimp,Vcato,Vxcato)

  call Potrmt(Full_atom,iaprotoi,ipr1,iprabs_nonexc,itypepr,n_atom_0,n_atom_ind,n_atom_proto,natome,nlm_pot,nrato,nrm,nspin,ntype, &
                  numat,Rato,rhoato,rhomft,Rmtg0,Vcato,Vcmft,Vxcato,Vxcmft,iapot)

! Calculation of the intersitial potential, on the FDM grid of points
  call Pot0(alfpot,Axe_Atom_gr,Bulk_atom_done,Chargat,drhoato,dvcato,Full_atom,Helm_cos,i_self, &
        ia_eq_inv,iaabs,iaproto,icheck(12),igreq,igroup,itypep,Magnetic,n_atom_0,n_atom_ind,n_atom_proto,n_atom_proto_bulk, &
        natomeq,natomp, &
        neqm,ngroup_m,Nonsph,npoint,npoint_ns,npsom,nrato,nrm,nspin,ntype,Orig_helm,pos,Rato,rho,rhons,rs,rhoigr,rhomft, &
        Rmtg0,Rmtsd,SCF,Sym_2D,V_helm,V_intmax,V_surf,Vato,Vcmft,Vh,Vhns,Vsphere,Vxc,Width_helm,xyz,z_bulk,z_surf)

! Writing: the atom charge is integrated up to Rmtsd, here before any superposition
  if( i_self == 1 .and. icheck(13) > 2 ) then
    ch(:) = 0._db
    write(3,140); write(3,150)
    do iapr = n_atom_0,n_atom_ind
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        ipr = iapr
        if( ( ipr == 0 .and. ipr1 == 1 ) .or. .not. Proto_calculated(ipr) ) cycle
      endif
      it = itypepr(ipr)
      r(:) = Rato(:,it)
      rayint = Rmtsd(ipr)
      do ispin = 1, nspin
        rhr2(:) = rho_no_sup(:,ispin,iapr) * r(:)**2
        ch(iapr) = ch(iapr) + quatre_pi * f_integr3(r,rhr2,0,nrm,rayint)
      end do
      write(3,160) iapr, ch(iapr)
    end do
  end if

  return
  110 format(/' i_self =',i3,', Full_atom =',l2,/ ' iapr0 =',i2,', n_iapr =',i3)
  120 format(/' Difference of charge and potential between excited and non excited atom (iapr =',i3,') :',/ &
          '       r      4*pi*r2*drho   dv_ex_nex')
  125 format(/' Difference of charge and potential between excited and non excited atom (iapr =',i3,') :',/ &
          '       r    4*pi*r2*drho(up) 4*pi*r2*drho(dn)  dv_ex_nex(up)  dv_ex_nex(dn)')
  130 format(1p,10e13.5)
  140 format(/' Before superposing the electronic densities')
  150 format(/'    ia      charge ')
  160 format(2x,i3,f9.3)
end

!***********************************************************************

! Calculation of the atomic potential

subroutine potato(icheck,it,itypepr,n_atom_proto,nlat,nlatm,nrato,nrm,nspin,ntype,numat, &
        popatm,popatv,psival,Rato,rhoigr,rhoit,Vato)

  use declarations
  implicit none

  integer:: DeuxZ, icalcul, icheck, ipr, ir, ispin, it, l, n_atom_proto, nl, nlatm, nr, nrm, nspin, ntype

  integer, dimension(0:ntype):: nlat, nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr

  real(kind=db):: dp, facspin, p1, p2, qpi

  real(kind=db), dimension(0:nrm):: r, rho, Vh, vhato
  real(kind=db), dimension(0:nrm,nlatm):: rhoval, vhval
  real(kind=db), dimension(0:nrm,0:n_atom_proto):: vato
  real(kind=db), dimension(0:ntype,nlatm):: popatv
  real(kind=db), dimension(0:nrm,nspin,0:n_atom_proto):: rhoigr
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm
  real(kind=db), dimension(0:nrm,0:ntype):: Rato, rhoit
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival

  if( icheck > 1 .and. it == 0 ) write(3,110)

! Cas des spheres vides
  if( numat(it) == 0 ) then
    nrato(it) = nrato(0)
    nr = nrato(0)
    Rato(0:nr,it) = Rato(0:nr,0)
    do ipr = 0,n_atom_proto
      if( itypepr(ipr) /= it ) cycle
      rhoigr(0:nr,:,ipr) = 0._db
      vato(:,ipr) = 0._db
    end do
    return
  endif

  nr = nrato(it)
  nl = nlat(it)
  r(0:nr) = Rato(0:nr,it)
  r(0) = 0._db

! Extrapolation at the center of the atom
  p1 = r(2) / ( r(2) - r(1) )
  p2 = 1 - p1
  rhoit(0,it) = p1 * rhoit(1,it) + p2 * rhoit(2,it)
  if( nl > 0 ) psival(0,1:nl,it) = 0._db

! Wave functions psival are in fact: sqrt(4*pi)*r*psi.
! rhoval is the density of the valence orbitals
  qpi = 0.25_db / pi
  do l = 1,nl
    rhoval(1:nr,l) = qpi * ( psival(1:nr,l,it) / r(1:nr ) )**2
    rhoval(0,l) = p1 * rhoval(1,l) + p2 * rhoval(2,l)
  end do

! Calculation of the atomic Hartree potential
  do icalcul = 1,1+nl
    l = icalcul - 1

    if( icalcul == 1 ) then
      rho(0:nr) = rhoit(0:nr,it)
    else
      rho(0:nr) = rhoval(0:nr,l)
    endif

    call Poisson(nr,nrm,r,rho,Vh)

    if( icalcul == 1 ) then
      DeuxZ = 2 * numat(it)
      Vhato(0) = - numat(it) * 100000._db / rydb
      Vhato(1:nr) = Vh(1:nr) - DeuxZ / r(1:nr)
    else
      Vhval(0:nr,l) = Vh(0:nr)
    endif

    if( icheck > 2 ) then
      if( icalcul == 1 ) then
        write(3,120) it, numat(it), icalcul
        do ir = 1,nr
          write(3,130) r(ir)*bohr, vhato(ir)*rydb
        end do
      else
        write(3,140) it, numat(it), icalcul
        do ir = 1,nr
          write(3,130) r(ir)*bohr, vhval(ir,l)*rydb
        end do
      endif
    endif

  end do

  facspin = 1._db / nspin

  do ipr = 0,n_atom_proto

    if( itypepr(ipr) /= it ) cycle

    vato(0:nr,ipr) = vhato(0:nr)
    do ispin = 1,nspin
      rhoigr(0:nr,ispin,ipr) = rhoit(0:nr,it) * facspin
    end do
! rhoit: the real density

! Here, dp take eventualy into account the magnetism for the calculation of the potential
    do l = 1,nl
      do ispin = 1,nspin
        dp = popatm(ipr,l,ispin) - popatv(it,l) * facspin
        Vato(0:nr,ipr) = Vato(0:nr,ipr) + dp * vhval(0:nr,l)
        rhoigr(0:nr,ispin,ipr) = rhoigr(0:nr,ispin,ipr) + dp * rhoval(0:nr,l)
      end do
    end do

    if( icheck > 1 ) then
      write(3,150) it, ipr
      do ir = 0,nr
        write(3,130) r(ir)*bohr, vato(ir,ipr)*rydb, rhoigr(ir,1:nspin,ipr)
      end do
    endif

  end do

  return
  110 format(/' ---- Potato -------',100('-'))
  120 format(/5x,'Rato        vhato(it)        it =',i3,', Z =',i3,', icalcul =',i3)
  130 format(1p,9e13.5)
  140 format(/5x,'Rato          vhval          it =',i3,', Z =',i3,', icalcul =',i3)
  150 format(/5x,'Rato        vato         rhoato(ispin=1,nspin)  it =',i3,', ipr =',i3)
end

!********************************************************************

! Le potentiel de Hartree est calcule a l'aide de la formulation integrale de l'equation de Poisson.
! Il y a 2 integrales a calculer pour chaque rayon r.
! Toutes deux sont obtenues par integration du troisieme degree.
! Resolution de l'equation de Poisson spherique
!       Laplacien(V) = - 8 * pi * rho (en u.a.).

subroutine Poisson(nr,nrm,r,rho,Vh)

  use declarations
  implicit none

  integer:: ir, j, nr, nrm

  real(kind=db):: prim
  real(kind=db), dimension(nrm):: a, b, c, d
  real(kind=db), dimension(0:nrm):: ch, cx, r, r2, rds, rho, Vh

  r2(0:nr) = r(0:nr)**2
  rds(0:nr) = huit_pi * r2(0:nr) * rho(0:nr)

! Calcul de la charge interieure en fonction du rayon.
  ch(0) = 0._db
  call coefpol3(r,rds,a,b,c,d,nr,nrm)
  do ir = 1,nr
    if( ir == 1 ) then
      j = ir
    elseif( ir < nr ) then
      j = ir - 1
    else
      j = ir - 2
    endif
    prim = a(j) * ( r(ir) + r(ir-1) )*( r2(ir-1) + r2(ir) ) / 4 + b(j) * ( r2(ir) + r(ir)*r(ir-1) + r2(ir-1) ) / 3
    prim = ( prim +  c(j) * ( r(ir) + r(ir-1) ) / 2 + d(j) ) * ( r(ir) - r(ir-1) )
    ch(ir) = ch(ir-1) + prim
  end do

! Calcul de l'integrale exterieure en fonction du rayon.

  rds(1:nr) = huit_pi * r(1:nr) * rho(1:nr)

  cx(nr) = 0._db
  call coefpol3(r,rds,a,b,c,d,nr,nrm)
  do ir = nr-1,0,-1
    if( ir == nr-1 ) then
      j = nr - 2
    elseif( ir == 0 ) then
      j = 1
    else
      j = ir
    endif
    prim = a(j) * ( r(ir) + r(ir+1) ) * ( r2(ir) + r2(ir+1) ) / 4 + b(j) * ( r2(ir) + r(ir) * r(ir+1) + r2(ir+1) ) / 3
    prim = ( prim + c(j) * ( r(ir) + r(ir+1) ) / 2 + d(j) ) * ( r(ir+1) - r(ir) )
    cx(ir) = cx(ir+1) + prim
  end do

  Vh(0) = cx(0)
  Vh(1:nr) = ch(1:nr) / r(1:nr) + cx(1:nr)

  return
end

!********************************************************************

! Resolution de l'equation de Poisson pour l /= 0
! On resoud d2/dr2 U_L - l(l+1)/r2 U_L = - 8 pi r rho_L

subroutine Poisson_lm(lmax_pot,nlm_pot,nr,r,rho_lm,Rmtsd,Vc_lm)

  use declarations
  implicit none

  integer:: icheck, ir, l, l0, lm, lmax_pot, lmp, m, nlm_pot, nm, nr, nrmtsd

  real(kind=db):: dr, r0, rap, rm, Rmtsd, rp
  real(kind=db), dimension(nr):: clapl0, claplm, claplp, diag, f2, r
  real(kind=db), dimension(nr,nlm_pot):: rho_lm, Vc_lm
  real(kind=db), dimension(:,:), allocatable:: Sm, Sol

  icheck = 1

! Coefficients du Laplacien
  do ir = 1,nr
    if( ir == 1 ) then
      rm = 0._db
    else
      rm = r(ir-1)
    endif
    r0 = r(ir)
    if( ir == nr ) then
      rp = 2 * r0 - rm
    else
      rp = r(ir+1)
    endif
    dr = 0.5_db * ( rp - rm )
    claplm(ir) = 1 / ( ( r0 - rm ) * dr )
    claplp(ir) = 1 / ( ( rp - r0 ) * dr )
    clapl0(ir) = - claplm(ir) - claplp(ir)
  end do

  f2(1:nr) = 1 / r(1:nr)**2

  if( icheck > 2 ) then
    write(3,110)
    do ir = 1,nr
      write(3,120) r(ir)*bohr, claplm(ir), clapl0(ir), claplp(ir), f2(ir)
    end do
  endif

  if( icheck > 2 ) then
    l0 = 0
  else
    l0 = 1
  endif

  do ir = 1,nr
    if( r(ir) > Rmtsd ) exit
  end do
  nrmtsd = min( ir, nr )

  do l = l0,lmax_pot

    nm = 2*l + 1

    Diag(1:nr) = - l * ( l + 1 ) * f2(1:nr) + clapl0(1:nr)

    allocate( Sm(nr,nm) )
    allocate( Sol(nr,nm) )

    do m = -l,l
      lm = l + 1 + m
      lmp = l**2 + l + 1 + m
      Sm(1:nrmtsd,lm) = - huit_pi * r(1:nrmtsd)*rho_lm(1:nrmtsd,lmp)
    end do
    Sm(nrmtsd+1:nr,:) = 0._db

    if( icheck > 2 ) then
      write(3,130) l
      do ir = 1,nr
        write(3,120) r(ir)*bohr, Diag(ir), Sm(ir,:)
      end do
    endif

! Triangularisation
    do ir = nr-1,1,-1
      rap = claplp(ir) / diag(ir+1)
      Diag(ir) = Diag(ir) - rap * claplm(ir+1)
      Sm(ir,:) = Sm(ir,:) - rap * Sm(ir+1,:)
    end do

! Resolution
    Sol(1,:) = Sm(1,:) / Diag(1)
    do ir = 2,nr
      Sol(ir,:) = (Sm(ir,:) - claplm(ir) * Sol(ir-1,:)) / Diag(ir)
    end do

    if( l > 0 ) then
      do m = -l,l
        lm = l + 1 + m
        lmp = l**2 + l + 1 + m
        Vc_lm(1:nr,lmp) = Sol(1:nr,lm) / r(1:nr)
      end do
    endif

    if( icheck > 2 ) then
      write(3,140)
      do ir = 1,nr
        write(3,120) r(ir)*bohr, Rydb * Sol(ir,:) / r(ir)
      end do
    endif

    Deallocate( Sol, Sm )

  end do

  110 format(/' Poisson_lm',//5x, 'Radius       claplm       clapl0       claplp        f2')
  120 format(1p,100e13.5)
  130 format(/' l =',i2,/ 5x,'Radius',7x,'Diag',8x,'Sm(1)',8x,'Sm(2)',8x,'Sm(3)...')
  140 format(/5x,'Radius',7x,'Vc(1)',8x,'Vc(2)',8x,'Vc(3)...')

  return
end

!********************************************************************

subroutine coefpol3(x,y,a,b,c,d,nr,nrm)

  use declarations
  implicit none

  integer:: i, nr, nrm

  real(kind=db):: deter, xa, xb, xba, xc, xca, xcb, xd, xda, xdb, xdc, y1, y2, y3, ya, yb, yc, yd

  real(kind=db), dimension(0:nrm) :: x, y
  real(kind=db), dimension(nrm) :: a, b, c, d

  do i = 1,nr-2

    xa = x(i-1)
    xb = x(i)
    xc = x(i+1)
    xd = x(i+2)
    ya = y(i-1)
    yb = y(i)
    yc = y(i+1)
    yd = y(i+2)
    xdc = xd - xc
    xdb = xd - xb
    xda = xd - xa
    xcb = xc - xb
    xca = xc - xa
    xba = xb - xa

    deter = xdc * xdb * xda * xcb * xca * xba
    a(i) = - ya * xdc * xdb * xcb + yb * xdc * xda * xca - yc * xdb * xda * xba + yd * xcb * xca * xba
    a(i) = a(i) / deter
    y1 = ya - a(i) * xa**3
    y2 = yb - a(i) * xb**3
    y3 = yc - a(i) * xc**3
    b(i) = y1 / ( xba * xca ) - y2 / ( xba * xcb ) + y3 / ( xcb * xca )
    y1 = y1 - b(i) * xa**2
    y2 = y2 - b(i) * xb**2
    c(i) = ( y2 - y1 ) / xba
    d(i) = y1 - c(i) * xa

  end do

  return
end

!*********************************************************************

! Calculation of the electronic density on the FDM grid of point coming from the non spherical valence orbital

subroutine orbval(Hybrid,iaproto,iapot,icheck,igreq,igroup,itypepr,lvval,mpirank,n_atom_proto,natomeq, &
        natomp,neqm,ngroup_m,ngroup_nonsph,nhybm,nlat,nlatm,norbv,npoint,npoint_ns,npsom,nrato,nrm,ntype,pop_nonsph,pos,psival, &
        Rato,rhons,rhonspr,Rot_Atom_gr,Rot_int,Vhns,Vhnspr,xyz)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( lmaxval=3, lmx2=(2*lmaxval+1)**2 )

  complex(kind=db), dimension(nhybm,16,ngroup_nonsph) :: hybrid
  complex(kind=db), dimension(:), allocatable :: ylmc

  integer, dimension(0:ngroup_nonsph) :: norbv
  integer, dimension(natomp):: iaproto, igroup
  integer, dimension(0:n_atom_proto):: iapot, itypepr
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(0:ntype):: nlat, nrato
  integer, dimension(0:ntype,nlatm):: lvval

  logical, dimension(lmx2):: rho_nul

  real(kind=db), dimension(3):: v
  real(kind=db), dimension(npoint):: dist
  real(kind=db), dimension(npoint_ns):: rhons, Vhns
  real(kind=db), dimension(0:n_atom_proto):: rhonspr, vhnspr
  real(kind=db), dimension(0:nrm):: cdiag, claplm, claplp, clapl0, com, cop, f2, r, vdiag
  real(kind=db), dimension(3,3):: rot, Rot_int
  real(kind=db), dimension(3,3,ngroup_m):: Rot_Atom_gr
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(nhybm,ngroup_nonsph) :: pop_nonsph
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(0:nrm,lmx2):: rholm, vlm
  real(kind=db), dimension(0:nrm,0:ntype):: Rato
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(:), allocatable :: ylmr
  real(kind=db), dimension(:,:), allocatable :: ylm

  if(icheck > 0) write(3,110)

  rhons(:) = 0._db
  vhns(:) = 0._db

  do ipr = 0,n_atom_proto

    rholm(:,:) = 0._db
    rho_nul(:) = .true.

    it = itypepr(ipr)
    igr = igreq(ipr,1)
    if( norbv(igr) == 0 ) cycle

    nr = nrato(it)
    r(:) = Rato(:,it)
    do ir = 1,nr
      rm = r(ir-1)
      r0 = r(ir)
      if( ir == nr ) then
        rp = 2 * r0 - rm
      else
        rp = r(ir+1)
      endif
      dr = 0.5_db * ( rp - rm )
      claplm(ir) = 1 / ( ( r0 - rm ) * dr )
      claplp(ir) = 1 / ( ( rp - r0 ) * dr )
      clapl0(ir) = - claplm(ir) - claplp(ir)
      f2(ir) = 1 / r(ir)**2
    end do

    lm = 1
    do l = 1,2*lmaxval
      do m = -l,l
        lm = lm + 1
! Calcul de la densite electronique radiale atomique

        do io = 1,norbv(igr)

          lmp = 0
          do lp = 0,lmaxval
            do mp = -lp,lp
              lmp = lmp + 1
              if( abs(hybrid(io,lmp,igr)) < eps10 ) cycle
              do np = 1,nlat(it)
                if( lvval(it,np) == lp ) exit
              end do
              if( np == nlat(it) + 1 .and. mpirank == 0 ) then
                call write_error
                do iprt = 3,9,3
                  write(iprt,120) lp
                end do
                stop
              endif
              lmpp = 0
              do lpp = 0,lmaxval
                do mpp = -lpp,lpp
                  lmpp = lmpp + 1
                  if( abs(hybrid(io,lmpp,igr)) < eps10 ) cycle

                  g = Gaunt_r(l,m,lp,mp,lpp,mpp)
                  if( abs( g ) < eps10 ) cycle

                  do npp = 1,nlat(it)
                    if( lvval(it,npp) == lpp ) exit
                  end do
                  if( npp == nlat(it) + 1 .and. mpirank == 0 ) then
                    call write_error
                    do iprt = 3,9,3
                      write(iprt,120) lpp
                    end do
                    stop
                  endif
                  rho_nul(lm) = .false.
                  aap = g * pop_nonsph(io,igr) * real( conjg(hybrid(io,lmp,igr)) * hybrid(io,lmpp,igr), db )

                  rholm(1:nr,lm) = rholm(1:nr,lm) + aap * psival(1:nr,np,it) * psival(1:nr,npp,it) / r(1:nr)**2
                end do
              end do

            end do
          end do

        end do

! Resolution de l'equation de Poisson spherique
!       Laplacien(V) = - 8*pi*rho (en u.a.).

        if( rho_nul(lm) ) cycle

        cdiag(:) = ( clapl0(:) - l * ( l + 1 ) * f2(:) )
        com(:) = claplm(:) / cdiag(:)
        cop(:) = claplp(:) / cdiag(:)
        vdiag(:) = - 8 * pi * rholm(:,lm) / cdiag(:)

        vlm(0:nr,lm) = 0._db
        do iter = 1,10*nr
          do ir = nr-1,1,-1
            vlm(ir,lm) = vdiag(ir) - com(ir) * vlm(ir-1,lm) - cop(ir) * vlm(ir+1,lm)
          end do
        end do

        if( icheck > 1 ) then
          write(3,130) ipr, l, m
          do ir = 0,nr
            write(3,140) r(ir)*bohr, rholm(ir,lm), vlm(ir,lm)*rydb
          end do
        endif

      end do
    end do

! Calcul du potentiel et de la densite totale non spherique

! Calcul des densites electroniques

    lm = 1
    lmaxv = 0
    do l = 1,2*lmaxval
      do m = -l,l
        lm = lm + 1
        if( .not. rho_nul(lm) ) lmaxv = l
      end do
    end do

    nlmr = ( lmaxv + 1 )**2
    nlmc = ( ( lmaxv + 1 ) * ( lmaxv + 2 ) ) / 2

    allocate( ylmc(nlmc) )
    allocate( ylmr(nlmr) )
    allocate( ylm(npoint,nlmr) )

    boucle_ia: do ia = 1,natomeq

      if( iaproto(ia) /= ipr ) cycle

      igr = igroup(ia)
      rot(:,:) = Rot_Atom_gr(:,:,igr)
      rot = matmul( transpose( rot_int ), rot )

      if( icheck > 0 ) then
        write(3,150) ipr
        write(3,160) ( rot(i,1:3), i = 1,3 )
      endif
! On prend la rotation inverse pour amener les points vers l'orbitale
      rot = transpose( rot )

      do i = 1,npoint
        v(1:3) = xyz(1:3,i) - pos(1:3,ia)
        v = matmul( rot, v )
        dist(i) = sqrt( sum( v(:)**2 ) )
        call cylm(lmaxv,v,dist(i),ylmc,nlmc)
        call ylmcr(lmaxv,nlmc,nlmr,ylmc,ylmr)
        ylm(i,:) = ylmr(:)
      end do

      do i = 1,npoint

        do ir = 1,nrato(it)
          if( Rato(ir,it) > dist(i) ) exit
        end do
        if( ir > nrato(it) ) cycle

        p1 = ( Rato(ir,it) - dist(i) ) / ( Rato(ir,it) - Rato(ir-1,it) )
        p2 = 1 - p1

        lm = 1
        do l = 1,lmaxv
          do m = -l,l
            lm = lm + 1
            if( rho_nul(lm) ) cycle

            rhons(i) = rhons(i) + ylm(i,lm) * ( p1 * rholm(ir-1,lm) + p2 * rholm(ir,lm) )
            vhns(i) = vhns(i) + ylm(i,lm) * ( p1 * vlm(ir-1,lm) + p2 * vlm(ir,lm) )

          end do
        end do

      end do

! On fait la meme chose pour les atomes pour calculer leur eventuel shift en energie

      i = 1
      do iprb = 0,n_atom_proto

        if( iapot(iprb) == 0 ) cycle
        ib = iapot(iprb)
        v(1:3) = pos(1:3,ib) - pos(1:3,ia)
        v = matmul( rot, v )
        dist(i) = sqrt( sum( v(:)**2 ) )
        call cylm(lmaxv,v,dist(i),ylmc,nlmc)
        call ylmcr(lmaxv,nlmc,nlmr,ylmc,ylmr)
        ylm(i,:) = ylmr(:)

        do ir = 1,nrato(it)
          if( Rato(ir,it) > dist(i) ) exit
        end do
        if( ir > nrato(it) ) cycle

        p1 = ( Rato(ir,it) - dist(i) ) / ( Rato(ir,it) - Rato(ir-1,it) )
        p2 = 1 - p1

        lm = 1
        do l = 1,lmaxv
          do m = -l,l
            lm = lm + 1
            if( rho_nul(lm) ) cycle

            rhonspr(iprb) = rhonspr(iprb) + ylm(i,lm) * ( p1 * rholm(ir-1,lm) + p2 * rholm(ir,lm) )
            vhnspr(iprb) = vhnspr(iprb) + ylm(i,lm) * ( p1 * vlm(ir-1,lm) + p2 * vlm(ir,lm) )

          end do
        end do
      end do

    end do boucle_ia

    deallocate( ylmc )
    deallocate( ylmr )
    deallocate( ylm )

  end do

  if( icheck > 0 ) then
    write(3,170)
    do ipr = 0,n_atom_proto
      if( iapot(ipr) == 0 ) cycle
      write(3,180) ipr, pos(1:3,iapot(ipr))*bohr, rhonspr(ipr), vhnspr(ipr)*rydb
    end do
  endif

  if( icheck > 1 ) then
    write(3,190)
    do i = 1,npoint
      v(1:3) = abs( xyz(1:3,i) )
      if( icheck > 2 .or. ( (v(1) < eps10 .and. v(2) < eps10 ) .or. (v(1) < eps10 .and. v(3) < eps10 ) .or. &
           (v(2) < eps10 .and. v(3) < eps10 ) ) ) write(3,200) xyz(1:3,i)*bohr, rhons(i), vhns(i)*rydb
    end do
  endif

  return
  110 format(/' ---- Orbval -------',100('-'))
  120 format(///' The orbital l =',i3,' is not defined under atom !'///)
  130 format(/'   igr =',i3,',  l =',i2,', m =',i2,/ '   Radius       rholm         Vlm')
  140 format(f10.7,1p,2e13.5)
  150 format(/' Matrix rotation of the atom',i2,' :')
  160 format(3f8.4)
  170 format(/'  ipr    posx    posy    posz     rhonspr        Vnspr')
  180 format(i4,2x,3f8.4,1p,2e13.5)
  190 format(/'    x       y       z        rhons         Vns')
  200 format(3f8.3,1p,2e13.5)
end

!***********************************************************************

! Routine making the superposition to calculate the Hartree (Coulomb) potential, the exchange-correlation potential and
! the electronic density in the ground state

subroutine Pot0muffin(alfpot,Cal_xanes,Chargat,chargat_init,chargat_self, &
        drho_ex_nex,drhoato,dvc_ex_nex,dvcato,Exc,Full_atom,Helm_cos,i_self,ia_eq_inv_self,iaproto,iapot,iapr, &
        icheck,ipr,iprabs,itypep,itypepr,lmax_pot,Magnetic,n_atom_0,n_atom_0_self, &
        n_atom_ind_self,n_atom_proto,natome,natome_self,natomeq,natomeq_self,natomp,nlm_pot,nonexc,nrato,nrm,nrm_self, &
        nspin,ntype,numat,Orig_helm,pos,posi,Rato,rho_chg,rho_no_sup,rho_self,rhoato,rhoato_init,rhoigr, &
        rhonspr,Rmtsd,rsato,SCF,Self_nonexc,Sym_2D,This_bulk_atom_done,V_helm,V_surf,Vato,Vcato,Vcato_init,Vhnspr, &
        Vsphere,Vxcato,Width_helm,z_bulk,z_surf)

  use declarations
  implicit none

! number of theta values for the integration in the sphere
  integer, parameter:: nthetam = 60
  integer, parameter:: ntpm = 2 * nthetam**2

  integer:: i, i_self, ia, iai, iapr, ib, ibb, icheck, ipr, iprabs, iprb, ir, ir0, ispin, it, itia, l, lm, lmax_pot, m, &
    n_atom_0, n_atom_0_self, n_atom_ind_self, n_atom_proto, natome, &
    natome_self, natomeq, natomeq_self, natomp, nlm_pot, nr, nrm, nrm_self, nrmin, nspin, ntheta, ntheta1, ntype

  integer, dimension(natomeq_self):: ia_eq_inv_self
  integer, dimension(natomp):: iaproto, itypep
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: iapot, itypepr

  logical:: Atom_self, This_bulk_atom_done, Cal_xanes, Full_atom, Helm_cos, Magnetic, Nonexc, Self_nonexc, SCF, Sym_2D

  real(kind=db):: alfa, alfpot, Cos_theta, Dab, dab2, dabr, dch, Delta_V_helm, dist, drhm, dt2, dtheta, dvcm, f, &
    fac, Orig_helm, p1, p2, ra1, ra2, ra3, Rayop, rhonspr, Rmtsd, rpotmin, Sin_theta, Smooth, Theta, Tiers, V_helm, Vrop, &
    Vsphere, Vhnspr, Width_helm, Width_p, x, z_bulk, z_surf

  real(kind=db), dimension(3):: p, V_surf
  real(kind=db), dimension(n_atom_0_self:n_atom_ind_self,nspin):: chargat_init, chargat_self
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(nthetam):: dvc
  real(kind=db), dimension(nthetam,0:lmax_pot):: YlmSintdt
  real(kind=db), dimension(nspin):: rhoigrop
  real(kind=db), dimension(nlm_pot):: Dlm0
  real(kind=db), dimension(0:n_atom_proto):: Chargat
  real(kind=db), dimension(nthetam,nspin):: drh
  real(kind=db), dimension(0:nrm,0:n_atom_proto):: vato
  real(kind=db), dimension(0:nrm):: drhoato, dvcato, r, r2, rsato, Vcato_init
  real(kind=db), dimension(0:nrm,nlm_pot):: Vcato
  real(kind=db), dimension(nrm):: Exc, dvc_ex_nex
  real(kind=db), dimension(nrm,nspin):: drho_ex_nex
  real(kind=db), dimension(0:nrm,nspin):: rho_chg, rho_no_sup, rhoato, rhoato_init
  real(kind=db), dimension(0:nrm,nlm_pot,nspin):: Vxcato
  real(kind=db), dimension(0:nrm,nspin,0:n_atom_proto):: rhoigr
  real(kind=db), dimension(0:nrm_self,nlm_pot,nspin, n_atom_0_self:n_atom_ind_self):: rho_self
  real(kind=db), dimension(0:nrm,0:ntype):: Rato

  real(kind=db), dimension(:), allocatable:: ra, rst, Vr, Vrp
  real(kind=db), dimension(:,:), allocatable :: rho_lm, rhot, Vc_lm, Vxct

  if( .not. This_bulk_atom_done ) Vcato(:,:) = 0._db
  Vxcato(:,:,:) = 0._db

  if( icheck > 1 .and. ( iapr < n_atom_0 .or. ( i_self > 1 .and. iapr == n_atom_0 ) ) ) write(3,110)
  if( iapr < n_atom_0 ) then
    ia = iapot(0)
  elseif( Full_atom ) then
    if( iapr == 0 ) then
      do ia = 1,natomp
        if( sum( abs( posi(:,iapot(0)) - pos(:,ia) ) ) < eps10 ) exit
      end do
    else
      do ia = 1,natomp
        if( sum( abs( posi(:,iapr) - pos(:,ia) ) ) < eps10 ) exit
      end do
    endif
  else
    ia = iapot(iapr)
    if( ia == 0 .and. ipr == iprabs ) ia = iapot(0)
  endif
  if( ia == 0 ) return
  itia = itypepr(ipr)

  nr = nrato(itia)

  r(0:nr) = Rato(0:nr,itia)
  r2(0:nr) = r(0:nr)**2

  Atom_self = .false.

  if( SCF .and. ( ( i_self > 1 .and. ia <= natomeq_self ) .or. This_bulk_atom_done ) ) then
    if( Full_atom .and. cal_xanes .and. ia <= natomeq .and. ( ( ia /= iapot(0) .and. ipr /= 0 ) .or. ipr == 0 ) ) then
      do iai = 1,natome_self
        if( ia_eq_inv_self(ia) /= iai ) cycle
        Atom_self = .true.
        exit
      end do
    else
      if( Cal_xanes .and. iapr == 0 .and. self_nonexc .and. .not. Full_atom ) then
        iai = iprabs
      else
        iai = iapr
      endif
      Atom_self = .true.
    endif
  endif

  if( Atom_self ) then

! Calculation of the Hartree potential by resolution of the Poisson equation

    do ispin = 1,nspin
      rhoato(0:nr,ispin) = rho_self(0:nr,1,ispin,iai)
    end do

    if( .not. This_bulk_atom_done ) then

      do ir = 0,nr
        drhoato(ir) = sum( rhoato(ir,:) - rhoato_init(ir,:) )
      end do

      call Poisson(nr,nrm,r,drhoato,dvcato)

      Vcato(0:nr,1) = dvcato(0:nr) + Vcato_init(0:nr)

      if( Cal_xanes .and. ipr == 0 .and. self_nonexc .and. .not. nonexc ) then
        rhoato(1:nr,:) = rhoato(1:nr,:) + drho_ex_nex(1:nr,:)
        Vcato(1:nr,1) = Vcato(1:nr,1) + dvc_ex_nex(1:nr)
      endif

      if( nlm_pot > 1 ) then

        allocate( rho_lm(nr,nlm_pot) )
        allocate( Vc_lm(nr,nlm_pot) )
        allocate( ra(nr) )

        ra(1:nr) = r(1:nr)

        do lm = 1,nlm_pot
          do ir = 1,nr
            rho_lm(ir,lm) = sum( rho_self(ir,lm,:,iai) )
          end do
        end do

        call Poisson_lm(lmax_pot,nlm_pot,nr,ra,rho_lm,Rmtsd,Vc_lm)

        do lm = 2,nlm_pot
          Vcato(0,lm) = 0._db
          Vcato(1:nr,lm) = Vc_lm(1:nr,lm)
        end do

        deallocate( ra, rho_lm, Vc_lm )

      endif

    endif

  else

    dvcato(:) = 0._db

    Vcato(0:nr,1) = Vato(0:nr,ipr) + Vsphere + Vhnspr

    if( abs( V_helm ) > eps10 ) then
      if( .not. Helm_cos ) then
    ! alfa is such that sqrt(pi/2) * erf( alfa/sqrt(2) ) / alfa = 0.5
        alfa = 2.47304793497921_db
        Width_p = Width_helm / ( sqrt(2._db) * alfa )
        fac = 0.5_db * sqrt(pi) * Width_p * V_helm
      endif
      if( Sym_2D ) then
        Dist = sum( V_surf(1:3) * pos(1:3,ia) )
        if( z_bulk > -999._db .and. z_surf > -999._db ) then
          Smooth = ( Dist - z_bulk ) / ( z_surf - z_bulk )
          Smooth = max( Smooth, 0._db )
        else
          Smooth = 1._db
        endif
      else
        Dist = sqrt( sum( pos(:,ia)**2 ) )
        Smooth = 1._db
      endif
      if( Helm_cos ) then
        if( dist > Orig_helm - Width_helm .and. dist < Orig_helm + Width_helm ) then
          Delta_V_helm = Smooth * 0.5_db * V_helm * ( 1 + cos( pi * ( dist - Orig_helm ) / Width_helm ) )
          Vcato(0:nr,1) = Vcato(0:nr,1) + Delta_V_helm
          if( icheck > 1 ) write(3,*)
          if( icheck > 0 ) write(3,115) ia, ipr, numat(itia), Smooth, Delta_V_helm * rydb
        endif
      else
        x = ( dist - Orig_helm ) / Width_p
        if( abs( x ) < eps10 ) then
          Delta_V_helm = V_helm
        else
          Delta_V_helm = Smooth * fac * ( 1 - erfc( x ) ) / ( dist - Orig_helm )
        endif
        Vcato(0:nr,1) = Vcato(0:nr,1) + Delta_V_helm
        if( icheck > 1 ) write(3,*)
        if( icheck > 0 ) write(3,115) ia, ipr, numat(itia), Smooth, Delta_V_helm *rydb
      endif
    endif

    do ispin = 1,nspin
      rhoato(0:nr,ispin) = rhoigr(0:nr,ispin,ipr) + rhonspr
    end do

! Storage of the electronic density before superposition
    if( i_self == 1 ) rho_no_sup(:,:) = rhoato(:,:)

  endif

  if( icheck > 2 ) then
    write(3,120) ia, iapr, ipr, numat(itia)
    if( nspin == 1 ) then
      write(3,140)
    else
      write(3,150)
    endif
    do ir = 1,nr
      write(3,160) r(ir)*bohr, quatre_pi * rhoato(ir,1:nspin) * r2(ir), Vcato(ir,1:nlm_pot)*rydb
    end do
  endif

  rpotmin = 0.1_db / bohr

  do ir = 1,nr-1
    if( Rato(ir,itia) > rpotmin ) exit
  end do
  nrmin = ir

  dtheta = pi / nthetam
  dt2 = dtheta / 2
  ra1 = sqrt(3 / quatre_pi)
  ra2 = 0.5_db * sqrt(45 / quatre_pi)
  ra3 = 0.5_db * sqrt(175 / quatre_pi)
  Tiers = 1 / 3._db

  do i = 1,nthetam
    Theta = ( i - 0.5_db ) * dtheta
    sin_theta = sin( theta )
    if( lmax_pot > 0 ) cos_theta = cos( theta )
    do l = 0,lmax_pot
      select case(l)
        case(0)
! A peu pres 0.5*sin(theta)*dtheta  ( pour l = 0, on fait la simple moyenne )
          YlmSintdt(i,l) = sin_theta * sin(dt2)
        case(1)
! Y(l,0)*Sin(theta)*dtheta
          YlmSintdt(i,l) = ra1 * cos_theta * sin_theta * dtheta
        case(2)
          YlmSintdt(i,l) = ra2 * ( cos_theta**2 - Tiers ) * sin_theta * dtheta
        case(3)
          YlmSintdt(i,l) = ra3 * ( cos_theta**2 - 0.6_db ) * cos_theta * sin_theta * dtheta
        case default
          YlmSintdt(i,l) = ( sqrt(4*l**2 - 1._db) / l ) * ( cos_theta * YlmSintdt(i,l-1) &
               - ( ( l - 1._db ) / sqrt((2*l-1._db)*(2*l-3._db)) ) * YlmSintdt(i,l-2) )
      end select
    end do
  end do

! Superposition

  do ib = 1,natomp

    if( This_bulk_atom_done ) cycle

    it = itypep(ib)
    iprb = iaproto(ib)
    if( iprb == 0 .and. self_nonexc ) iprb = iprabs

! on exclut la superposition de l'atome avec lui même...

    p(1:3) = pos(1:3,ia) - pos(1:3,ib)
    Dab = sqrt( sum( p(:)**2 ) )
    if( Dab < epspos ) then
      if( i_self == 1 .and. ib == natomeq .and. .not. This_bulk_atom_done ) rho_chg(:,:) = rhoato(:,:)
      cycle
    endif

    dab2 = Dab**2
! Si l'atome est a plus de 10 ua 3.17 A, on fait un calcul moins precis
    if( Dab > 10._db ) then
      ntheta = 1
    else
      ntheta = nthetam
    endif

    if( Atom_self ) then

! Density is in this case already superposed
      if( ib <= natomeq_self ) then
        if( ia_eq_inv_self(ib) <= natome_self ) then
          if( Full_atom ) then
            ibb = ia_eq_inv_self(ib)
          else
            ibb = iprb
          endif
          dch = sum( chargat_self(ibb,:) - chargat_init(ibb,:) )
          dvcm = - 2 * dch / Dab
          do ir = 0,nrato(it)
            Vcato(ir,1) = Vcato(ir,1) + dvcm
          end do
        endif
      else
        exit
      endif

    elseif( Dab > Rato(nrato(it),it) ) then

      dvcm = - 2 * Chargat(iaproto(ib)) / Dab
      Vcato(0:nr,1) = Vcato(0:nr,1) + dvcm
      if( i_self == 1 .and. ib == natomeq .and. .not. This_bulk_atom_done ) rho_chg(:,:) = rhoato(:,:)

    else

! Determination du rayon en dessous duquel on considere le potentiel constant.
      do ir = 1,nrato(it)
        if( vato(ir,iprb) > -1._db ) exit
      end do
      ir = min(ir,nrato(it))
      Rayop = Rato(ir,it)
      vrop = vato(ir,iprb)
      rhoigrop(1:nspin) = rhoigr(ir,1:nspin,iprb)

! Calculation of Wigner matrices
     if( nlm_pot > 1 .and. ntheta > 1 ) call Cal_Dlm0(Dab,Dlm0,lmax_pot,nlm_pot,p)

      do ir0 = 0,nr

        dabr = 2 * Dab * r(ir0)
        if( ir0 < nrmin ) then
          ntheta1 = 1
        else
          ntheta1 = ntheta
        endif
        do i = 1,ntheta1
          if( ntheta1 == 1 ) then
            dist = sqrt( dab2 + r2(ir0) )
          else
            theta = ( i - 0.5_db ) * dtheta
            dist = sqrt( dab2 + r2(ir0) - dabr * cos( theta ) )
          endif

          if( dist > Rato(nrato(it),it) ) then
            dvc(i) = - 2 * Chargat(iprb) / dist
            drh(i,1:nspin) = 0._db
          elseif( dist < Rayop ) then
            dvc(i) = vrop
            drh(i,1:nspin) = rhoigrop(1:nspin)
          else
            do ir = 1,nrato(it)
              if( Rato(ir,it) > dist ) exit
            end do
            p1 = (dist-Rato(ir-1,it)) / (Rato(ir,it)-Rato(ir-1,it))
            p2 = 1 - p1
            dvc(i) = p1 * vato(ir,iprb) + p2 * vato(ir-1,iprb)
            do ispin = 1,nspin
              drh(i,ispin) = p1 * rhoigr(ir,ispin,iprb) + p2 * rhoigr(ir-1,ispin,iprb)
            end do
          endif
          if( Magnetic ) then
            drh(i,1) = 0.5_db * sum( drh(i,1:nspin) )
            drh(i,nspin) = drh(i,1)
          endif

        end do
        do l = 0,lmax_pot
          if( ntheta1 == 1 ) then
            if( l > 0 ) cycle
            dvcm = dvc(1)
          else
            dvcm = sum( dvc(1:ntheta1) * YlmSintdt(1:ntheta1,l) )
          endif
          do m = -l,l
            if( l == 0 ) then
              Vcato(ir0,1) = Vcato(ir0,1) + dvcm
            else
              lm = l**2 + l + 1 + m
              Vcato(ir0,lm) = Vcato(ir0,lm) + Dlm0(lm) * dvcm
            endif
          end do
        end do

        do ispin = 1,nspin
          if( ntheta1 == 1 ) then
            drhm = drh(1,ispin)
          else
            drhm = sum( drh(1:ntheta1,ispin) * YlmSintdt(1:ntheta1,0))
          endif
          rhoato(ir0,ispin) = rhoato(ir0,ispin) + drhm

          end do
        end do

! On fait sortir la densité electronique pour le calcul de la charge de
! l'agregat et la distance pour l'atome le plus eloigne du petit agregat:
! natome = nombre d'atomes dans l'agregat symetrise
! natomeq = nombre d'atomes dans le petit agregat
      if( i_self == 1 .and. ib == natomeq ) then
        rho_chg(:,:) = rhoato(:,:)
      endif

    endif

  end do  ! end of superposition

! Densite venant de l'exterieur de l'agregat par superposition
  if( i_self == 1 .and. .not. This_bulk_atom_done ) rho_chg(:,:) = rhoato(:,:) - rho_chg(:,:)

  do ispin = 1,nspin
    do ir = 0,nr
      rhoato(ir,ispin) = max( rhoato(ir,ispin), 1.e-15_db )
    end do
  end do

! Calcul du rayon de Fermi rs.
  f = 0.75_db / pi
  do ir = 0,nr
    rsato(ir) = ( f / sum( rhoato(ir,1:nspin) ) )**tiers
  end do

  allocate( rhot(nr,nspin) )
  allocate( Vxct(nr,nspin) )
  allocate( rst(nr) )

  do ispin = 1,nspin
    rhot(1:nr,ispin) = rhoato(1:nr,ispin)
  end do
  rst(1:nr) = rsato(1:nr)
  call potxc(Magnetic,nr,nrm,nspin,alfpot,rhot,rst,Vxct,exc)
  do ispin = 1,nspin
    Vxcato(1:nr,1,ispin) = Vxct(1:nr,ispin)
    if( .not. Atom_self ) cycle
    do lm = 2,nlm_pot
      Vxcato(1:nr,lm,ispin) = Vxct(1:nr,ispin) * rho_self(1:nr,lm,ispin,iai) / ( 3 * rhoato(1:nr,ispin) )
    end do
  end do

  deallocate( rhot )
  deallocate( Vxct )
  deallocate( rst )

! Le potentiel doit etre toujours croissant.
  if( numat(itia ) /= 0 ) then
    Allocate( Vrp(nspin), Vr(nspin) )
    Vrp(:) = Vcato(1,1) + Vxcato(1,1,:)
    do ir = 2,nr
      Vr(:) = Vcato(ir,1) + Vxcato(ir,1,:)
      if( Vr(1) < Vrp(1) .or. Vr(nspin) < Vrp(nspin) ) then
        Vcato(ir,1) = Vcato(ir-1,1)
        Vxcato(ir,1,:) = Vxcato(ir-1,1,:)
        rhoato(ir,:) = rhoato(ir-1,:)
        rsato(ir) = rsato(ir-1)
        drhoato(ir) = drhoato(ir-1)
      endif
      Vrp(:) = Vcato(ir,1) + Vxcato(ir,1,:)
    end do
    Deallocate( Vrp, Vr )
  endif

  if( icheck > 1 ) then
    write(3,132) ia, iapr
    if( Atom_self ) then
      if( nspin == 1 .and. nlm_pot == 1 ) then
        write(3,170)
      elseif( nspin == 1 .and. nlm_pot > 1 ) then
        write(3,175)
      elseif( nspin == 2 .and. nlm_pot == 1 ) then
        write(3,180)
      else
        write(3,185)
      endif
      do ir = 1,nr
        dvc(1:nspin) = Vcato(ir,1) + Vxcato(ir,1,1:nspin)

        write(3,160) r(ir)*bohr, dvc(1:nspin)*rydb, Vcato(ir,1)*rydb, Vxcato(ir,1,1:nspin)*rydb, &
            quatre_pi * rhoato(ir,1:nspin) * r2(ir), rsato(ir), quatre_pi * rho_chg(ir,1:nspin) * r2(ir), &
            quatre_pi * drhoato(ir) * r2(ir), dvcato(ir)*rydb, Vcato(ir,2:nlm_pot)*rydb, Vxcato(ir,2:nlm_pot,1:nspin)*rydb
      end do
    else
      if( nspin == 1 .and. nlm_pot == 1 ) then
        write(3,190)
      elseif( nspin == 1 .and. nlm_pot > 1 ) then
        write(3,195)
      elseif( nspin == 2 .and. nlm_pot == 1 ) then
        write(3,200)
      else
        write(3,205)
      endif
      do ir = 1,nr
        dvc(1:nspin) = Vcato(ir,1) + Vxcato(ir,1,1:nspin)

        write(3,160) r(ir)*bohr, dvc(1:nspin)*rydb, Vcato(ir,1)*rydb, Vxcato(ir,1,1:nspin)*rydb, &
            quatre_pi * rhoato(ir,1:nspin) * r2(ir), rsato(ir), quatre_pi * rho_chg(ir,1:nspin) * r2(ir), &
            Vcato(ir,2:nlm_pot)*rydb
      end do
    endif
  endif

  return
  110 format(/' ---- Pot0muffin ',98('-'))
  115 format(' ia =',i3,', ipr =',i3,', Z =',i3,', Smooth =',f6.3,', Delta_V_helmholtz =',f10.5,' eV')
  120 format(/' ia =',i3,', iapr =',i3,', ipr =',i3,', Z =',i3, '. Avant superposition')
  132 format(/' ia =',i3,', iapr =',i3)
  140 format(4x,'Radius     4pi*r2*rho    Vcato(lm=1,..nlm_pot)')
  150 format(4x,'Radius    4pi*r2*rho_u 4pi*r2*rho_d   Vcato(lm=1,..nlm_pot)')
  160 format(1p,50e13.5)
  170 format(6x,'Rato        Vato        Vcato        Vxcato     4pi*r2*rho      rsato   4pi*r2*rho_chg 4pi*r2*d_rho   d_vcato')
  175 format(6x,'Rato        Vato        Vcato        Vxcato     4pi*r2*rho      rsato   4pi*r2*rho_chg 4pi*r2*d_rho   d_vcato', &
                '     Vcato(lm=2,nlm_pot)',104x,'Vxcato(lm=2,nlm_pot)')
  180 format(6x,'Rato       Vato(u)      Vato(d)        Vcato      Vxcato(u)    Vxcato(d)  4pi*r2*rho(u) 4pi*r2*rho(d)', &
                '   rsato 4pi*r2*rho_chg(u) 4pi*r2*rho_chg(d) 4pi*r2*d_rho   d_vcato')
  185 format(6x,'Rato       Vato(u)      Vato(d)        Vcato      Vxcato(u)    Vxcato(d)  4pi*r2*rho(u) 4pi*r2*rho(d)', &
                '   rsato 4pi*r2*rho_chg(u) 4pi*r2*rho_chg(d) 4pi*r2*d_rho   d_vcato     Vcato(lm=2,..nlm_pot)')
  190 format(6x,'Rato        Vato        Vcato        Vxcato     4pi*r2*rho      rsato   4pi*r2*rho_chg')
  195 format(6x,'Rato        Vato        Vcato        Vxcato     4pi*r2*rho      rsato   4pi*r2*rho_chg     Vcato(lm=2,nlm_pot)')
  200 format(6x,'Rato       Vato(u)      Vato(d)       Vcato      Vxcato(u)    Vxcato(d)  4pi*r2*rho(u) 4pi*r2*rho(d)', &
                '   rsato 4pi*r2*rho_chg(u) 4pi*r2*rho_chg(d)')
  205 format(6x,'Rato       Vato(u)      Vato(d)       Vcato      Vxcato(u)    Vxcato(d)  4pi*r2*rho(u) 4pi*r2*rho(d)', &
                '   rsato 4pi*r2*rho_chg(u) 4pi*r2*rho_chg(d)     Vcato(lm=2,nlm_pot)')
end

!***********************************************************************

! Calcul des matrices de rotation de Wigner quand le premier indice est nul
! Dlm0 = racine(4pi/(l+1))*Ylm

subroutine Cal_Dlm0(Dab,Dlm0,lmax_pot,nlm_pot,p)

  use declarations
  implicit none

  integer:: l, lm, lmax_pot, m, nlm_pot, nlmc

  complex(kind=db), dimension(:), allocatable:: Ylmc

  real(kind=db):: Dab, r
  real(kind=db), dimension(3):: p
  real(kind=db), dimension(nlm_pot):: Dlm0, Ylmr

  nlmc = ( ( lmax_pot + 1 ) * ( lmax_pot + 2 ) ) / 2

  allocate( Ylmc(nlmc) )

  call cylm(lmax_pot,p,Dab,Ylmc,nlmc)
  call ylmcr(lmax_pot,nlmc,nlm_pot,Ylmc,Ylmr)

  deallocate( Ylmc )

  lm = 0
  do l = 0,lmax_pot
    r = sqrt( quatre_pi / (2*l + 1 ) )
    do m = -l,l
      lm = lm + 1
      Dlm0(lm) = r * Ylmr(lm)
    end do
  end do

  return
end

!***********************************************************************

! Force the atom equivalent in the space group and at the same distance from the center of the sphere of calculation
! to have the same potential

subroutine Force_pot_eq(Do_init,drhoato,dVcato,Excato,iaprotoi,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self, &
                        natome,nlm_pot,nrm,nrm_self,nspin,posi,rho_chg,rho_no_sup,rhoato,rhoato_init,rsato,Vcato,Vcato_init,Vxcato)

  use declarations
  implicit none

  integer:: i, ia, ib, n_at_eq, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, natome, nlm_pot, nrm, nrm_self, nspin

  integer, dimension(natome):: iaprotoi, index_at_eq

  logical:: Do_init, First_couple

  logical, dimension(natome):: Done

  real(kind=db):: Dist
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(nrm):: Exc
  real(kind=db), dimension(0:nrm):: dVc, drh, rsa
  real(kind=db), dimension(0:nrm,nlm_pot):: Vca
  real(kind=db), dimension(0:nrm,nspin):: rho, rho_no
  real(kind=db), dimension(0:nrm,nlm_pot,nspin):: Vxc

  real(kind=db), dimension(nrm,n_atom_0:n_atom_ind):: Excato
  real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: dVcato, drhoato, rsato
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nspin,n_atom_0:n_atom_ind):: rhoato, rho_no_sup
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato

  real(kind=db), dimension(0:nrm_self):: Vca_init
  real(kind=db), dimension(0:nrm_self,nspin):: rho_ch, rho_init
  real(kind=db), dimension(0:nrm_self,n_atom_0_self:n_atom_ind_self):: Vcato_init
  real(kind=db), dimension(0:nrm_self,nspin,n_atom_0_self:n_atom_ind_self):: rho_chg, rhoato_init

! We are necessarily with  n_atom_0 = 1
!                          n_atom_ind = natome

  Done(:) = .false.
  index_at_eq(:) = 0

  do ia = 1,natome

    if( Done(ia) ) cycle

    Done(ia) = .true.
    Dist = sum( posi(:,ia)**2 )
    First_couple = .true.
    n_at_eq = 0

    do ib = ia + 1,natome
      if( Done(ib) ) cycle
      if( iaprotoi(ia) /= iaprotoi(ib) ) cycle
      if( abs( sum( posi(:,ib)**2 ) - Dist ) > eps10 ) cycle
      Done(ib) = .true.
      if( First_couple ) then
        First_couple = .false.
        n_at_eq = 1
        index_at_eq(1) = ia
      endif
      n_at_eq = n_at_eq + 1
      index_at_eq(n_at_eq) = ib
    end do

    if( n_at_eq == 0 ) cycle

    Exc(:) = 0._db
    dVc(:) = 0._db
    drh(:) = 0._db
    rsa(:) = 0._db
    Vca(:,:) = 0._db
    rho(:,:) = 0._db
    rho_no(:,:) = 0._db
    Vxc(:,:,:) = 0._db
    Vca_init(:) = 0._db
    rho_init(:,:) = 0._db
    rho_ch(:,:) = 0._db
    do i = 1,n_at_eq
      ib = index_at_eq(i)
      Exc(:) = Exc(:) + Excato(:,ib) / n_at_eq
      dVc(:) = dVc(:) + dVcato(:,ib) / n_at_eq
      drh(:) = dVc(:) + drhoato(:,ib) / n_at_eq
      rsa(:) = rsa(:) + rsato(:,ib) / n_at_eq
      Vca(:,:) = Vca(:,:) + Vcato(:,:,ib) / n_at_eq
      rho(:,:) = rho(:,:) + rhoato(:,:,ib) / n_at_eq
      rho_no(:,:) = rho_no(:,:) + rho_no_sup(:,:,ib) / n_at_eq
      Vxc(:,:,:) = Vxc(:,:,:) + Vxcato(:,:,:,ib) / n_at_eq
      if( Do_init ) then
        Vca_init(:) = Vca_init(:) + Vcato_init(:,ib) / n_at_eq
        rho_init(:,:) = rho_init(:,:) + rhoato_init(:,:,ib) / n_at_eq
        rho_ch(:,:) = rho_ch(:,:) + rho_chg(:,:,ib) / n_at_eq
      endif
    end do
    do i = 1,n_at_eq
      ib = index_at_eq(i)
      Excato(:,ib) = Exc(:)
      dVcato(:,ib) = dVc(:)
      drhoato(:,ib) = drh(:)
      rsato(:,ib) = rsa(:)
      Vcato(:,:,ib) = Vca(:,:)
      rhoato(:,:,ib) = rho(:,:)
      rho_no_sup(:,:,ib) = rho_no(:,:)
      Vxcato(:,:,:,ib) = Vxc(:,:,:)
      if( Do_init ) then
        Vcato_init(:,ib) = Vca_init(:)
        rhoato_init(:,:,ib) = rho_init(:,:)
        rho_chg(:,:,ib) = rho_ch(:,:)
      endif
    end do

  end do

  return
end

!***********************************************************************

subroutine Raymuf(Cal_xanes,Chargat,Full_atom,i_self,iapot,iaproto,iaprotoi,icheck,iprabs, &
        ipr1,itab,itypei,itypep,itypepr,n_atom_0,n_atom_ind,n_atom_proto,natome,natomeq,natomp,nlm_pot, &
        normrmt,nrato,nrm,nspin,ntype,numat,Overlap,pos,Rato,rhoato,Rmtg,Rmtg0,Rmtimp, &
        Rmtsd,Rsort,V0bdcFimp,Vcato,Vxcato)

  use declarations
  implicit none

  integer:: i_self, ia, iapr, iaprb, ib, icheck, ipr, ipr1, ipra, iprabs, iprb, ir, ira, irb, it, ita, itab, itb, &
    n_atom_0, n_atom_ind, n_atom_proto, natome, natomeq, natomp, nlm_pot, normrmt, nr, nra, nrb, nrm, nspin, ntype, Z

  integer, dimension(natomp):: iaproto, itypep
  integer, dimension(natome):: iaprotoi, itypei
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: iapot, iaproxp, itypepr

  logical:: Cal_xanes, Full_atom

  real(kind=db):: a1, a2, b1, b2, chtot, dab_ov, dist, Overlap, p1, rb, rm, Rsort, V0bdcFimp, Vr, Vr1, Vrop, &
                  Vropmax

  real(kind=db), dimension(3):: ps
  real(kind=db), dimension(0:nrm):: r, rhr2, vra, vrb
  real(kind=db), dimension(0:ntype):: Rmtimp
  real(kind=db), dimension(0:n_atom_proto):: Chargat, Dab, Rayop, Rdem, Rmtg, Rmtg0, Rmtsd, Rn, Rnorm, RV0
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nspin,n_atom_0:n_atom_ind):: rhoato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vxcato
  real(kind=db), dimension(0:nrm,0:ntype):: Rato
  real(kind=db), dimension(3,natomp):: pos

  if( Cal_xanes .or. i_self == 1 ) then

  if( icheck > 0 ) write(3,110)

  iaproxp(:) = 0
  Rmtg(:) = 0._db
  Rmtg0(:) = 0._db

  if( natomp == 1 ) then
    if( normrmt == 4 ) then
      rm = min( Rsort, Rmtimp( 1 ) )
    else
      rm = min( Rsort, 2.5_db / bohr )
    endif
    Dab(:) = 0._db
    Rdem(:) = rm
    Rnorm(:) = rm
    Rn(:) = rm
    Rayop(:) = rm
    Rmtg(:) = rm
    Rmtg0(:) = rm
    RV0(:) = rm
    Vrop = 0._db

  else

! Norman radius
    do ipr = ipr1,n_atom_proto
      Rnorm(ipr) = 0._db
      if( iapot(ipr) == 0 ) cycle
      it = itypepr(ipr)
      if( numat(it) == 0 ) cycle
      if( Full_atom ) then
        do ia = n_atom_0,n_atom_ind
          if( iaprotoi(ia) == ipr ) exit
        end do
        if( ia > n_atom_ind ) then
          do ib = 1,n_atom_ind
            if( numat( itypei(ib) ) == numat( it ) ) exit
          end do
          if( ib > n_atom_ind ) then
            Rnorm(ipr) = 0._db
            cycle
          else
            iapr = ib
          endif
        else
          iapr = ia
        endif
      else
        iapr = ipr
      endif
      nr = nrato(it)
      r(0:nr) = Rato(0:nr,it)
      do ir = 0,nr
        rhr2(ir) = sum( rhoato(ir,1:nspin,iapr) ) * r(ir)**2
      end do
      chtot = numat(it) - Chargat(ipr)
      call rnorman(r,rhr2,nr,Rnorm(ipr),chtot,nrm)
    end do

    Rn(:) = 0._db

    do ipr = ipr1, n_atom_proto
      if( natomp < 2 ) then
        Rn(ipr) = Rsort
      else
        do ia = 1,natomp
          if( iaproto(ia) == ipr ) exit
        end do
        if( ia == natomp + 1 ) then
          iaproxp(ipr) = 0
          Dab(ipr) = Dab(0)
          Rdem(ipr) = Rdem(0)
          Rn(ipr) = Rn(0)
        else
          Dab(ipr) = 100000._db
          do ib = 1,natomp
            if( numat(itypep(ib)) == 0 ) cycle
            ps(1:3) = pos(1:3,ia) - pos(1:3,ib)
            dist = sqrt( sum( ps(:)**2 ) )
            if( dist < Dab(ipr) - epspos .and. dist > epspos ) then
              iaproxp(ipr) = ib
              Dab(ipr) = dist
              Rdem(ipr) = 0.5_db * Dab(ipr)
            endif
          end do
          ib = iaproxp(ipr)
          iprb = iaproto(ib)
          if( Rnorm(ipr) > eps10 .and. Rnorm( iprb ) > eps10 ) then
            Rn(ipr) = min( 0.75_db * Dab(ipr), Dab(ipr) / ( 1 + Rnorm(iprb) / Rnorm(ipr) ) )
          else
            Rn(ipr) = min( 0.75_db * Dab(ipr), Rnorm(ipr) )
          endif
        endif
      endif
    end do

    do ipr = ipr1,n_atom_proto
      if( iaproxp(ipr) /= 0 ) cycle
      iaproxp(ipr) = natomp + 1
      Dab(ipr) = Dab(iprabs)
      Rdem(ipr) = Rdem(iprabs)
      Rn(ipr) = Rn(iprabs)
    end do

! Optimized radius

    Vropmax = 0._db

    if( n_atom_proto == 1 ) then

      Rayop(ipr1:1) = 0.5_db * ( 1 + Overlap ) * Dab(1)
      Vrop = 0._db

    else

      vra(:) = 0._db; vrb(:) = 0._db

! Calculation of the cutting potential
      if( ipr1 == 0 ) then
        ipra = ipr1
      else
        ipra = iprabs
      endif
      if( Full_atom ) then
        do ia = n_atom_0,n_atom_ind
          if( itypei(ia) == itab ) exit
        end do
        iapr = ia
      else
        iapr = ipra
      endif
      ita = itab
      nra = nrato( ita )
      do ir = 1,nra
        Vra(ir) = Vcato(ir,1,iapr) + sum(Vxcato(ir,1,1:nspin,iapr)) / nspin
        if( ir == 1 ) cycle
        if( Vra(ir) < Vra(ir-1) ) Vropmax = vra(ir-1) - eps6
      end do
      iprb = iaproto( iaproxp(ipra) )
      if( Full_atom ) then
        do ia = n_atom_0,n_atom_ind
          if( iaprotoi(ia) == iprb ) exit
        end do
        if( ia > n_atom_ind ) then ! case with only 1 atom
          iaprb = iapr
        else
          iaprb = ia
        endif
      else
        iaprb = iprb
      endif
      itb = itypepr(iprb)
      nrb = nrato( itb )
      do ir = 1,nrb
        Vrb(ir) = Vcato(ir,1,iaprb) + sum( Vxcato(ir,1,1:nspin,iaprb) ) / nspin
        if( vrb(ir) < vrb(ir-1) ) Vropmax = min( Vrb(ir-1)-eps6, Vropmax )
      end do

      dab_ov = ( 1 + Overlap ) * Dab(ipra)

      do ira = 2,nra-1
        rb = dab_ov - Rato(ira,ita)
        do irb = nrb,2,-1
          if( Rato(irb,itb) < rb ) exit
        end do
        if( Vrb( irb ) < Vra( ira ) ) exit
      end do
      a1 = ( Vra(ira+1) - Vra(ira) ) / ( Rato(ira+1,ita) - Rato(ira,ita) )
      a2 = ( Vrb(irb) - Vrb(irb-1) ) / ( Rato(irb-1,itb) - Rato(irb,itb) )
      b1 = Vra(ira) - a1 * Rato(ira,ita)
      b2 = Vrb(irb) - a2 * ( dab_ov - Rato(irb,itb) )
      if( abs(a1) < eps10 ) then
        Vrop = vra(ira)
      else
        Vrop = ( b1*a2 - b2*a1 ) / ( a2 - a1 )
      endif
      Vrop = min( Vrop, Vropmax )

      if( icheck > 2 ) then
        write(3,120) dab_ov*bohr, numat(ita), ita, iapr, numat(itb), itb, iaproxp(ipra), iaprb, Vropmax*rydb, ira, irb
        write(3,125)
        do ir = 1,min(nra,nrb)
          write(3,130) Rato(ir,ita)*bohr, Vra(ir)*rydb, (dab_ov - Rato(ir,itb))*bohr, Vrb(ir)*rydb
        end do
      endif

      boucle_ia: do ipr = ipr1, n_atom_proto
        Rayop(ipr) = 0._db
        if( iapot(ipr) == 0 ) cycle
        it = itypepr(ipr)
        if( Full_atom ) then
          do ia = n_atom_0,n_atom_ind
            if( iaprotoi(ia) == ipr ) exit
          end do
          if( ia > n_atom_ind ) then
            do ib = 1,n_atom_ind
              if( numat( itypei(ib) ) == numat( it ) ) exit
            end do
            if( ib > n_atom_ind ) then
              iapr = 1
            else
              iapr = ib
            endif
          else
            iapr = ia
          endif
        else
          iapr = ipr
        endif
        if( numat(it) == 0 ) cycle
        vra(1) = Vcato(1,1,iapr) + sum( Vxcato(1,1,1:nspin,iapr) ) / nspin
        do ir = 2,nrato( it )
          vra(ir) = Vcato(ir,1,iapr) + sum(Vxcato(ir,1,1:nspin,iapr)) / nspin
          if( vra(ir) > vrop - eps10 ) exit
          if( vra(ir) < vra(ir-1) + eps10 ) then
            Rayop(ipr) = Rato(ir-1,it)
            cycle boucle_ia
          endif
        end do
        p1 = ( vra(ir) - vrop ) / ( vra(ir) - vra(ir-1) )
        Rayop(ipr) = p1 * Rato(ir-1,it) + ( 1 - p1 ) * Rato(ir,it)
        Rayop(ipr) = min( Rayop(ipr), Dab(ipr) )
      end do boucle_ia

    endif

! Hydrogen case
    do ipr = ipr1,n_atom_proto
      if( iapot(ipr) == 0 .or. Rn(ipr) < eps10 ) cycle
      Z = numat( itypepr(ipr)  )
      if( Z /= 1 ) cycle
      if( Rayop(ipr) > 0.45_db*Rdem(ipr) .and. Rayop(ipr) < 1.75_db*Rdem(ipr) ) cycle
      Rayop(ipr) = min( (1 + Overlap) * Rdem(ipr), 0.8_db / bohr )
    end do

    do ipr = ipr1,n_atom_proto
      RV0(ipr) = 0._db
      if( iapot(ipr) == 0 ) cycle
      it = itypepr(ipr)
      if( Full_atom ) then
        do ia = n_atom_0,n_atom_ind
          if( iaprotoi(ia) == ipr ) exit
        end do
        if( ia > n_atom_ind ) then
          do ib = 1,n_atom_ind
            if( numat( itypei(ib) ) == numat( it ) ) exit
          end do
          if( ib > n_atom_ind ) then
            iapr = 1
          else
            iapr = ib
          endif
        else
          iapr = ia
        endif
      else
        iapr = ipr
      endif
      Vr1 = Vcato(1,1,iapr) + sum( Vxcato(1,1,1:nspin,iapr) ) / nspin
      do ir = 2,nrato( it )
        Vr = Vcato(ir,1,iapr) + sum( Vxcato(ir,1,1:nspin,iapr) )/nspin
        if( ir == nrato(it) ) then
          RV0(ipr) = Rato(ir,it)
          exit
        elseif( Rato(ir,it) > Dab(ipr) .and. natomp > 1 ) then
          RV0(ipr) = Dab(ipr)
          exit
        elseif( Vr < Vr1 .and. numat(it) /= 0 ) then
          RV0(ipr) = Rato(ir-1,it)
          exit
        elseif( ( Vr - V0bdcFimp )*( vr1 - V0bdcFimp ) <= 0._db ) then
          p1 = ( Vr - V0bdcFimp ) / ( Vr - vr1 )
          RV0(ipr) = p1 * Rato(ir-1,it) + ( 1 - p1 ) * Rato(ir,it)
          exit
        endif
        Vr1 = Vr
      end do
    end do

    do ipr = ipr1,n_atom_proto
      if( iapot(ipr) == 0 .or. ( numat(itypepr(ipr)) == 0 .and. normrmt /= 4 ) ) then
        Rmtg0(ipr) = Rdem(ipr)
        Rmtg(ipr) = (1 + Overlap) * Rmtg0(ipr)
        cycle
      elseif( Rn(ipr) < eps10 ) then
        Rmtg0(ipr) = 0._db
        Rmtg(ipr) = 0._db
        cycle
      endif
      select case(normrmt)
        case(1)
          Rmtg(ipr) = Rayop(ipr)
          Rmtg0(ipr) = Rmtg(ipr) / ( 1 + Overlap )
        case(2)
          Rmtg0(ipr) = Rn(ipr)
          Rmtg(ipr) = (1 + Overlap) * Rmtg0(ipr)
        case(3)
          ! Suppressed
        case(4)
          if( Rmtimp( itypepr(ipr) ) < eps10 ) then
            Rmtg(ipr) = Rayop(ipr)
            Rmtg0(ipr) = Rmtg(ipr) / ( 1 + Overlap )
          else
            Rmtg0(ipr) = Rmtimp( itypepr(ipr) )
            Rmtg(ipr) = Rmtg0(ipr)
          endif
        case(5)
          Rmtg0(ipr) = RV0(ipr)
          Rmtg(ipr) = Rmtg0(ipr)
      end select
    end do

    do ipr = ipr1,n_atom_proto
      it = itypepr(ipr)
      Rmtg(ipr) = min( Rmtg(ipr), Rato(nrato(it),it) )
      Rmtg0(ipr) = min( Rmtg0(ipr), Rato(nrato(it),it) )
    end do

! When optimized radius are not ok, one uses Normann criterion
    if( normrmt == 1 ) then
      do ipr = ipr1,n_atom_proto
        if( iapot(ipr) == 0 .or. Rn(ipr) < eps10 ) cycle
        Z = numat( itypepr(ipr)  )
        if( Z < 2 ) cycle
        ia = iaproxp(ipr)
        if( ia > natomp ) cycle
        Z = numat( itypep(ia) )
        if( ( Z == 1 .and. ( Rayop(ipr) > 0.99_db*Rdem(ipr) .and. Rayop(ipr) < 1.75_db*Rdem(ipr) ) ) &
            .or. ( Z /= 1 .and. ( Rayop(ipr) > 0.25_db*Rdem(ipr) .and. Rayop(ipr) < 1.5_db*Rdem(ipr) ) ) ) then
          cycle
        elseif( ( Z == 1 .and. ( Rn(ipr) > 0.99_db*Rdem(ipr) .and. Rn(ipr) < 1.75_db*Rdem(ipr) ) ) &
            .or. ( Z /= 1 .and. ( Rn(ipr) > 0.25_db*Rdem(ipr) .and. Rn(ipr) < 1.5_db*Rdem(ipr) ) ) ) then
          if( iapot(ipr) <= natomeq ) then
            Rmtg0(:) = Rn(:)
            Rmtg(:) = (1 + Overlap) * Rmtg0(:)
          else
            Rmtg0(ipr) = Rn(ipr)
            Rmtg(ipr) = (1 + Overlap) * Rmtg0(ipr)
          endif
        else
          Rmtg0(ipr) = Rdem(ipr)
          Rmtg(ipr) = (1 + Overlap) * Rmtg0(ipr)
        endif
      end do
    endif

  endif ! arriving when natomp = 1

! Case of atoms outside the sphere
! Give an effect for the potential calculation inside the MT radius of the corresponding atoms overlaping with the outersphere.
! One artificially chooses the potential of the nonexcited absorbing atom.
  do ipr = ipr1,n_atom_proto
    if( Rmtg(ipr) > eps10 ) cycle
    Rmtg(ipr) = Rmtg(iprabs)
    Rmtg0(ipr) = Rmtg0(iprabs)
  end do

  Rmtsd(:) = Rmtg(:)
  if( ipr1 == 1 ) Rmtsd(0) = Rmtsd(iprabs)

  if( icheck > 0 ) then
    if( .not. n_atom_proto == 1 ) write(3,140) Vrop*rydb
    write(3,150)
    write(3,160)
    do ipr = ipr1,n_atom_proto
      if( iapot(ipr) == 0 .or. Rmtg0(ipr) < eps10 ) cycle
      Z = numat( itypepr(ipr) )
      write(3,170) ipr, iapot(ipr), Z, Rn(ipr)*bohr, Rnorm(ipr)*bohr, Dab(ipr)*bohr, Rdem(ipr)*bohr, iaproxp(ipr), &
              Rayop(ipr)*bohr, RV0(ipr)*bohr, Rmtsd(ipr)*bohr, Rmtg(ipr)*bohr
    end do
  endif

  endif  ! arriving when i_self /= 1 and not cal_xanes

  return
  110 format(/' ---- Raymuf -------',100('-'))
  120 format(/' Vrop calculation: dab_ov =',f6.3,/ &
              '  Central atom: Z =',i3,', ita =',i2,', iapr  =',i2,/ &
              ' Neighbor atom: Z =',i3,', itb =',i2,', iaprb =',i3,', ib =',i3, / &
              ' Vrop_max =',f8.3,' eV, ira =',i4,', irb =',i4)
  125 format(/'     Rato          vra      dab_ov-Rato     vrb')
  130 format(1p,4e13.5)
  140 format(/' Vrop = ',f10.3,' eV',/)
  150 format(' Rmtg : muffin-tin radius',/ &
             ' Rmtsd : Radius for the density of state calculation'/)
  160 format('  ipr ia   Z     Rn     Rnorm     Dab     Rdem  iaprox  Rayop     RV0     Rmtsd     Rmtg')
  170 format(3i4,4f9.5,i5,4f9.5)
end

!***********************************************************************

! Calculation of the potential and density at the muffin-tin radius without Overlap

subroutine Potrmt(Full_atom,iaprotoi,ipr1,iprabs,itypepr,n_atom_0,n_atom_ind,n_atom_proto,natome,nlm_pot,nrato,nrm,nspin,ntype, &
                  numat,Rato,rhoato,rhomft,Rmtg0,Vcato,Vcmft,Vxcato,Vxcmft,iapot)

  use declarations
  implicit none

  integer:: iapr, ipr, ipr1, iprabs, iprb, ir, it, n_atom_0, n_atom_ind, n_atom_proto, natome, nlm_pot, nrm, nrmtg0, nspin, ntype

  logical:: Full_atom

  integer, dimension(natome):: iaprotoi
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: iapot, itypepr

  real(kind=db):: p1
  real(kind=db), dimension(0:n_atom_proto):: Vcmft
  real(kind=db), dimension(0:n_atom_proto,nspin):: rhomft, Vxcmft
  real(kind=db), dimension(0:n_atom_proto):: Rmtg0
  real(kind=db), dimension(0:nrm,0:ntype):: Rato
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nspin,n_atom_0:n_atom_ind):: rhoato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato

  Vcmft(:) = 0._db
  rhomft(:,:) = 0._db
  vxcmft(:,:) = 0._db

  do ipr = ipr1,n_atom_proto

    if( Full_atom ) then
      do iapr = n_atom_0,n_atom_ind
        if( ipr == iaprotoi(iapr) ) exit
      end do
      if( iapr > n_atom_ind ) cycle
    else
      iapr = ipr
    endif

    it = itypepr(ipr)
    do ir = 2,nrato(it)-1
      if( Rato(ir,it) > Rmtg0(ipr) - eps10 ) exit
    end do
    nrmtg0 = ir

    p1 = (Rmtg0(ipr) - Rato(nrmtg0-1,it)) / ( Rato(nrmtg0,it) - Rato(nrmtg0-1,it) )
    Vcmft(ipr) = p1 * Vcato(nrmtg0,1,iapr) + ( 1-p1 ) * Vcato(nrmtg0-1,1,iapr)
    Vxcmft(ipr,1:nspin) = p1 * Vxcato(nrmtg0,1,1:nspin,iapr) + ( 1 - p1 ) * Vxcato(nrmtg0-1,1,1:nspin,iapr)
    rhomft(ipr,1:nspin) = p1 * rhoato(nrmtg0,1:nspin,iapr) + ( 1 - p1 ) * rhoato(nrmtg0-1,1:nspin,iapr)

  end do

  do ipr = ipr1,n_atom_proto

    if( abs( Vcmft(ipr) ) > eps10 ) cycle

! Search of an atom of same Z
    if( Full_atom .and. .not. iapot(ipr) == 0 ) then
      do iapr = n_atom_0,n_atom_ind
        iprb = iaprotoi(iapr)
        if( ipr == iprb ) exit
      end do
      if( iapr > n_atom_ind ) then
        do iapr = n_atom_0,n_atom_ind
          iprb = iaprotoi(iapr)
          if( numat( itypepr( ipr ) ) == numat( itypepr( iprb ) ) ) exit
        end do
      endif
      if( iprb <= n_atom_ind ) then
        Vcmft(ipr) = Vcmft(iprb)
        Vxcmft(ipr,1:nspin) = Vxcmft(iprb,1:nspin)
        rhomft(ipr,1:nspin) = rhomft(iprb,1:nspin)
        cycle
      endif
   endif

    Vcmft(ipr) = Vcmft(iprabs)
    Vxcmft(ipr,1:nspin) = Vxcmft(iprabs,1:nspin)
    rhomft(ipr,1:nspin) = rhomft(iprabs,1:nspin)

  end do

  return
end

!***********************************************************************

! Calculation of Norman radius such that
!       integral(0..Rnorm)(4*pi*r2*rho*dr) = Z

subroutine Rnorman(r,rh,nr,Rnorm,z,nrm)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db), dimension(nrm):: r, rh

  charge = 0._db
  zs4pi = z / quatre_pi

  do ir = 2,nr-1
    x1 = ( rh(ir+1) - rh(ir) ) / ( r(ir+1) - r(ir) )
    x2 = ( rh(ir-1) - rh(ir) ) / ( r(ir-1) - r(ir) )
    a = (  x1 - x2  ) / ( r(ir+1) - r(ir-1) )
    b = x1 - a * ( r(ir+1) + r(ir) )
    c = rh(ir) - a*r(ir)**2 - b*r(ir)
    if(ir == 2) then
      r2 = 0.5_db * ( r(ir+1) + r(ir) )
      r1 = 0._db
    elseif(ir == nr-1) then
      r2 = r(nr)
      r1 = 0.5_db * ( r(ir) + r(ir-1) )
    else
      r2 = 0.5_db * ( r(ir+1) + r(ir) )
      r1 = 0.5_db * ( r(ir) + r(ir-1) )
    endif
    dcharge = (a/3) * (r2**3 - r1**3) + (b/2) * (r2**2 - r1**2) + c * (r2 - r1)
    charge = charge + dcharge
    if( charge > zs4pi ) exit
  end do
  if( ir < nr-1 ) then
    ch2 = quatre_pi * charge
    ch1 = quatre_pi * ( charge - dcharge )
    p = ( z - ch1 ) / ( ch2 - ch1 )
    Rnorm = p * r(ir) + ( 1 - p ) * r(ir-1)
  else
    Rnorm = r(nr)
  endif

  return
end

!***********************************************************************

! Routine de superposition calculant le potentiel de Hartree total, le potentiel d'echange-correlation total, et la densite electronique
! totale dans l'etat fondamental.
! Mise a part la contribution venant de orbval, on neglige l'interference entre orbitales que ce soit entre atomes ou a
! l'interieur d'un meme atome.
! Dans chaque atome les orbitales sont toutes a symetrie spherique.

subroutine Pot0(alfpot,Axe_Atom_gr,Bulk_atom_done,Chargat,drhoato,dvcato,Full_atom,Helm_cos,i_self, &
        ia_eq_inv,iaabs,iaproto,icheck,igreq,igroup,itypep,Magnetic,n_atom_0,n_atom_ind,n_atom_proto,n_atom_proto_bulk, &
        natomeq,natomp, &
        neqm,ngroup_m,Nonsph,npoint,npoint_ns,npsom,nrato,nrm,nspin,ntype,Orig_helm,pos,Rato,rho,rhons,rs,rhoigr,rhomft, &
        Rmtg0,Rmtsd,SCF,Sym_2D,V_helm,V_intmax,V_surf,Vato,Vcmft,Vh,Vhns,Vsphere,Vxc,Width_helm,xyz,z_bulk,z_surf)

  use declarations
  implicit none

  integer:: i, i_self, ia, iaabs, iae, icheck, iga, igr, ipr, ir, isp, ispin, it, itm, kgr, n_atom_0, n_atom_ind, n_atom_proto, &
            n_atom_proto_bulk, natomeq, natomp, neqm, ngroup_m, npoint, npoint_ns, npsom, nrm, nspin, ntype

  integer, dimension(npoint):: ia_close
  integer, dimension(nspin):: ispp
  integer, dimension(natomeq):: ia_eq_inv
  integer, dimension(natomp):: iaproto, igroup, itypep
  integer, dimension(0:ntype):: nrato
  integer, dimension(0:n_atom_proto,neqm):: igreq

  logical:: Bulk_atom_done, Full_atom, Helm_cos, Magnetic, Nonsph, SCF, Sym_2D, This_bulk_atom_done
  logical, dimension(npoint):: iok

  real(kind=db):: alfa, alfpot, cosang, Delta_V_helm, drho, dist, dist_min, dv, f, fac, Orig_helm, p1, p2, Smooth, Tiers, &
    V_helm, V_intmax, Vsphere, Width_helm, Width_p, x, z_bulk, z_surf

  real(kind=db), dimension(3):: V, V_surf
  real(kind=db), dimension(npoint):: rs, Vh
  real(kind=db), dimension(npoint_ns):: rhons, Vhns
  real(kind=db), dimension(npoint,nspin):: Vxc, rho
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_gr
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(0:n_atom_proto):: Chargat, Rmtg0, Rmtsd, Vcmft
  real(kind=db), dimension(0:n_atom_proto,nspin):: rhomft
  real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: drhoato, dvcato
  real(kind=db), dimension(0:nrm,0:n_atom_proto):: Vato
  real(kind=db), dimension(0:nrm,nspin,0:n_atom_proto):: rhoigr
  real(kind=db), dimension(0:nrm,0:ntype):: Rato
  real(kind=db), dimension(:), allocatable:: rst
  real(kind=db), dimension(:,:), allocatable:: rhot, Vxct
  real(kind=db), dimension(nrm):: temp ! on ne se sert pas de ce tableau

  if( icheck > 2 ) write(3,110)

  Vh(:) = 0._db
  Vxc(:,:) = 0._db
  rs(:) = 1.0e+05_db
  rho(:,:) = 0._db
  iok(:) = .false.

! When SCF, looking for the closest atom
  ia_close(:) = 0
  if( i_self > 1 .and. SCF ) then
    do i = 1,npoint
      dist_min = 10000._db
      do ia = 1,natomeq
        ipr = iaproto(ia)
        v(1:3) = xyz(1:3,i) - pos(1:3,ia)
        dist = sqrt( sum( v(:)**2 ) )
        if( dist > dist_min .or. dist > Rmtsd(ipr) ) cycle
        ia_close(i) = ia
      end do
    end do
  endif

! Calculation of Hartree potential and electronic density

  do ia = 1,natomp

    it = itypep(ia)
    ipr = iaproto(ia)

    This_bulk_atom_done = Bulk_atom_done .and. ipr > n_atom_proto - n_atom_proto_bulk

    if( ia <= natomeq .and. ( i_self > 1 .or. This_bulk_atom_done ) .and. SCF ) then
      if( Full_atom ) then
        iae = ia_eq_inv(ia)
      else
        iae = ipr
      endif
    endif

    itm = 1
    if( nspin == 2 ) then
      igr = igroup(ia)
      iga = igroup(iaabs)
      kgr = igreq(ipr,1)

      cosang = sum( Axe_Atom_gr(:,kgr) * Axe_Atom_gr(:,igr) )
      cosang = cosang * abs( sum( Axe_Atom_gr(:,iga) * Axe_Atom_gr(:,igr) ) )
      if( abs(cosang - 1) < eps4 ) then
        ispp(1) = 1
        ispp(nspin) = nspin
      elseif( abs(cosang + 1) < eps4 ) then
        ispp(1) = nspin
        ispp(nspin) = 1
      else
        itm = 0
      endif
    else
      ispp(1) = 1
    endif

    do i = 1,npoint
      if( iok(i) ) cycle
      v(1:3) = xyz(1:3,i) - pos(1:3,ia)
      dist = sqrt( sum( v(:)**2 ) )
! Si on tombe dans un atome chevauchant la frontiere exterieure,
! on prend le potentiel au niveau de son rayon muffin-tin.
      if( ia > natomeq ) then
        if( dist < Rmtg0(ipr) ) then
          vh(i) = Vcmft(ipr)
          if( itm == 1 ) then
            rho(i,1:nspin) = rhomft(ipr,ispp(1:nspin))
          else
            rho(i,1:nspin) = 0.5_db * ( rhomft(ipr,1 )+ rhomft(ipr,nspin) )
          endif
          iok(i) = .true.
          cycle
        endif
      endif
      do ir = 1,nrato(it)
        if( Rato(ir,it) > dist ) exit
      end do
      if( ir > nrato(it) ) then
        Vh(i) = Vh(i) - 2 * Chargat(ipr) / dist
        cycle
      endif
      p1 = ( dist - Rato(ir-1,it) ) / ( Rato(ir,it) - Rato(ir-1,it))
      p2 = 1 - p1
      vh(i) = vh(i) + p1*vato(ir,ipr) + p2*vato(ir-1,ipr)
      if( itm == 1 ) then
        do ispin = 1,nspin
          isp = ispp(ispin)
          rho(i,ispin) = rho(i,ispin) + p1*rhoigr(ir,isp,ipr) + p2*rhoigr(ir-1,isp,ipr)
        end do
      else
        rho(i,:) = rho(i,:) + 0.5_db * ( p1*rhoigr(ir,1,ipr) + p2*rhoigr(ir-1,1,ipr) &
                                       + p1*rhoigr(ir,nspin,ipr) + p2*rhoigr(ir-1,nspin,ipr) )
      endif

! Self-consistence :
      if( ia <= natomeq .and. ( i_self > 1 .or. This_bulk_atom_done ) .and. SCF ) then
        Vh(i) = Vh(i) + p1*dvcato(ir,iae) + p2*dvcato(ir-1,iae)
! On ne superpose pas les densites a l'interieur de l'atome le plus proche
        if( ia_close(i) == ia ) then
          drho = ( p1*drhoato(ir,iae) + p2*drhoato(ir-1,iae) ) / nspin
          rho(i,:) = rho(i,:) + drho
        endif
      endif

    end do
  end do

! On ajoute eventuellement la contribution non spherique des orbitales de valence.
  if( Nonsph ) then
    do i = 1,npoint
      if( iok(i) ) cycle
      rho(i,1:nspin) = rho(i,1:nspin) + rhons(i)
      vh(i) = vh(i) + vhns(i)
    end do
  endif

! One adds the contribution from the outer-sphere
  do i = 1,npoint
    if( iok(i) ) cycle
    vh(i) = vh(i) + Vsphere
  end do

  if( abs( V_helm ) > eps10 ) then
    if( .not. Helm_cos ) then
    ! alfa is such that sqrt(pi/2) * erf( alfa/sqrt(2) ) / alfa = 0.5
      alfa = 2.47304793497921_db
      Width_p = Width_helm / ( sqrt(2._db) * alfa )
      fac = 0.5_db * sqrt(pi) * Width_p * V_helm
    endif
    do i = 1,npoint
      if( Sym_2D ) then
        dist = sum( V_surf(1:3) * xyz(1:3,i) )
        if( z_bulk > -999._db .and. z_surf > -999._db .and. Dist < z_surf ) then
          Smooth = ( Dist - z_bulk ) / ( z_surf - z_bulk )
          Smooth = max( Smooth, 0._db )
        else
          Smooth = 1._db
        endif
      else
        Smooth = 1._db
        dist = xyz(4,i)
      endif
      if( Helm_cos ) then
        if( dist > Orig_helm - Width_helm .and. dist < Orig_helm + Width_helm ) then
          Delta_V_helm = Smooth * 0.5_db * V_helm * ( 1 + cos( pi * ( dist - Orig_helm ) / Width_helm ) )
          Vh(i) = Vh(i) + Delta_V_helm
        endif
      else
        x = ( dist - Orig_helm ) / Width_p
        if( abs( x ) < eps10 ) then
          Delta_V_helm = V_helm
        else
          Delta_V_helm = Smooth * fac * ( 1 - erfc( x ) ) / ( dist - Orig_helm )
        endif
        Vh(i) = Vh(i) + Delta_V_helm
      endif
    end do
  endif

! Calculation of the Fermi radius, rs
  f = 0.75_db / pi
  do isp = 1,nspin
    do i = 1,npoint
      rho(i,isp) = max( rho(i,isp), eps10 )
    end do
  end do
  tiers = 1.0_db / 3.0_db
  do i = 1,npoint
    rs(i) = ( f / sum( rho(i,1:nspin) ) )**tiers
  end do

  allocate( rhot(npoint,nspin) )
  allocate( Vxct(npoint,nspin) )
  allocate( rst(npoint) )

  do isp = 1,nspin
    rhot(1:npoint,isp) = rho(1:npoint,isp)
  end do
  rst(1:npoint) = rs(1:npoint)

! Calculation of the exchange-correlation potential in the ground state.
  call potxc(Magnetic,npoint,nrm,nspin,alfpot,rhot,rst,Vxct,temp)

  do isp = 1,nspin
    Vxc(1:npoint,isp) = Vxct(1:npoint,isp)
  end do

  deallocate( rhot )
  deallocate( Vxct )
  deallocate( rst )

  if( V_intmax < 1000._db ) then
    do isp = 1,nspin
      do i = 1,npoint
        dv = V_intmax - Vh(i) - Vxc(i,isp)
        if( dv < 0._db ) Vh(i) = Vh(i) + dv
      end do
    end do
  endif

  if( icheck > 2 ) then
    if( Nonsph ) then
      write(3,120)
      do i = 1,npoint
        write(3,130) i, xyz(1:3,i)*bohr, Vh(i)*rydb, Vhns(i)*rydb, rhons(i), rs(i), Vxc(i,1:nspin)*rydb
      end do
    else
      if( nspin == 1 ) then
        write(3,140)
      else
        write(3,145)
      endif
      do i = 1,npoint
        write(3,130) i, xyz(1:3,i)*bohr, Vh(i)*rydb, rs(i), Vxc(i,1:nspin)*rydb
      end do
    endif
  endif

  return
  110 format(/' ---- Pot0 ----------',100('-'))
  120 format(/4x,'i    x      y      z      Vh_(eV)     Vhns_(eV)     rhons       rs_(ua)',8x,'Vxc(1,nspin)_(eV)')
  130 format(i5,3f7.3,1p,6e13.5)
  140 format(/4x,'i    x      y      z      Vh_(eV)      rs_(ua)      Vxc_(eV)')
  145 format(/4x,'i    x      y      z      vh_(eV)      rs_(ua)    Vxc_up_(eV)  Vxc_dn_(eV)')
end

!***********************************************************************

subroutine potxc(Magnetic,np,nrm,nspin,alfpot,rhot,rst,Vxct,exct)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( nA = 6 )

  logical:: Magnetic

  real(kind=db):: mu_cf, mu_cp, mu_xp, nu_c
  real(kind=db), dimension(np):: rst
  real(kind=db), dimension(np,nspin):: rhot, Vxct
  real(kind=db), dimension(nrm):: exct

! Vxct = potentiel d'echange correlation
! exct = energie d echange correlation par particule

  tiers = 1._db/3._db
  qtiers = 4 * tiers
  fac = ( (18*pi)**tiers ) / pi
  alpha0 = -0.75_db * (18/pi**2)**tiers

  if( alfpot >= eps4 ) then

! X_alpha potential
    f = - 1.5_db * alfpot * fac
    Vxct(1:np,1) = f / rst(1:np)
    if( Magnetic ) Vxct(1:np,nspin) = Vxct(1:np,1)
    exct(1:min(np,nrm)) = 0.75_db * Vxct(1:min(np,nrm),1)

  elseif( abs(alfpot) < eps4 ) then
! Hedin and Lundqvist potential,
! J. Phys. C: Solid State Phys., 4, 2064 (1971)
! Von Barth and Hedin, J. Phys. C: Solid State Phys., 5, 1629 (1972)

    a = (1 / 2._db)**tiers
    gamma = qtiers * a / ( 1 - a )
! Pour les valeurs de c_p et r_p on garde les valeurs non polarise
! choisie aussi par Moruzzi Janak et William (1978). Pour r_f et c_f
! on prend aussi leurs parametres plutot que les originaux de Von Barth
! qui sont : c_p = 0.0504, r_p = 30., c_f = 0.0254, r_f = 75.
    c_p = 0.045_db
    r_p = 21._db
    c_f = c_p / 2
    r_f = 2._db**qtiers * r_p   ! = 52.9166841

    do ip = 1,np

      mu_xp = - fac / rst(ip)
      mu_cp = - c_p * log( 1 + r_p / rst(ip) )
      e_cp = - c_p * f_vonbarth( rst(ip) / r_p )

      if( nspin == 1 ) then

        Vxct(ip,1) = mu_xp + mu_cp

        if( ip <= nrm ) then  ! juste pour ne pas depasser les tableaux; les appels de pot0 sont pas signifiants
          ex = alpha0 / rst(ip)
          ect = e_cp
          exct(ip) = ex + ect
        end if

      else

        e_cf = - c_f * f_vonbarth( rst(ip) / r_f )
        nu_c = gamma * ( e_cf - e_cp )
        mu_cf = - c_f * log( 1 + r_f / rst(ip) )
        tau_c = mu_cf - mu_cp - qtiers * ( e_cf - e_cp )

        x = rhot(ip,1) / sum( rhot(ip,1:nspin) )
        fx = ( x**qtiers + (1-x)**qtiers - a ) / ( 1 - a )

        f1 = mu_xp + nu_c
        f2 = mu_cp  - nu_c + tau_c * fx
        Vxct(ip,1) = f1 * ( 2 * x )**tiers + f2
        Vxct(ip,nspin) = f1 * ( 2 - 2 * x )**tiers + f2
        if( ip <= nrm ) then
          ex = alpha0 / rst(ip) + mu_xp * fx / gamma
          ect = e_cp + fx * nu_c / gamma
          exct(ip) = ex + ect
        end if
      endif

    end do

  else

! Potentiel de correlation de Perdew et Wang
    do ip = 1,np

      if( Magnetic ) then
        dzeta = ( rhot(ip,1) - rhot(ip,nspin) ) / sum( rhot(ip,:) )
      else
        dzeta = 0._db
      endif
      Rs = rst(ip)
      call cor_perdew(Rs,dzeta,Vcup,Vcdn)

! On ajoute l'echange
      mu_xp = - fac / rst(ip)
      if( nspin == 1 ) then
        Vxct(ip,1) = Vcup + mu_xp
      else
        x = rhot(ip,1) / sum( rhot(ip,1:nspin) )
        Vxct(ip,1) = Vcup + mu_xp * ( 2 * x )**tiers
        Vxct(ip,nspin) = Vcdn + mu_xp * ( 2 - 2 * x )**tiers
      endif

! Mis en attendant (relation approchee prise du Xalpha)
      if( ip <= nrm ) exct(ip) = 0.75_db * Vxct(ip,1)

    end do

  endif

  return
end

!***********************************************************************

function f_vonbarth(x)

  use declarations
  implicit none
  real(kind=db):: f_vonbarth, x

  f_vonbarth = (1 + x**3) * log(1 + 1/x) + x / 2 - x**2 - 1/3._db

  return
end

!***********************************************************************

! Perdew and Wang correlation potential, PRB, 45, 13244 (1992-I)

subroutine cor_perdew(Rs,dzeta,Vcup,Vcdn)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  data gamma, fs0 / 0.5198421_db, 1.709921_db /
  data Tiers, QTiers / 0.333333333333333_db, 1.333333333333333_db /

  f = ( (1 + dzeta)**QTiers + (1 - dzeta)**QTiers - 2 ) / gamma

  call G_cor(0.0310907_db,0.21370_db,7.5957_db,3.5876_db,1.6382_db, 0.49294_db,1.00_db,Rs,Eu,dEu_rs)
  call G_cor(0.01554535_db,0.20548_db,14.1189_db,6.1977_db, 3.3662_db,0.62517_db,1.00_db,Rs,Ep,dEp_rs)
  call G_cor(0.0168869_db,0.11125_db,10.357_db,3.6231_db,0.88026_db, 0.49671_db,1.00_db,Rs,alf,dalf)

  alfac = - alf
  dalfac = - dalf

  dzeta4 = dzeta**4
  Ec = Eu * ( 1 - f*dzeta4 ) + Ep * f * dzeta4 + alfac * f * ( 1 - dzeta4 ) / fs0

  dEc_rs = dEu_rs * ( 1 - f * dzeta4 ) + dEp_rs * f * dzeta4 + dalfac * f * ( 1 - dzeta4 ) / fs0
  df = QTiers * ( (1+dzeta)**Tiers - (1-dzeta)**Tiers ) / gamma
  dEc_dz = 4 * (dzeta**3) * f * ( Ep - Eu - alfac / fs0 ) + df * ( dzeta4*Ep - dzeta4*Eu + ( 1 - dzeta4 )*alfac/fs0)
  Vcom = Ec - Rs * dEc_rs / 3 - dzeta * dEc_dz
  Vcup = Vcom + dEc_dz
  Vcdn = Vcom - dEc_dz

! Pour convertir de Hartree en Rydberg
  Vcup = 2 * Vcup
  Vcdn = 2 * Vcdn

  return
end

!***********************************************************************

subroutine G_cor(A,A1,B1,B2,B3,B4,P,Rs,GG,G_Rs)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  P1 = P + 1
  Q0 = - 2 * A * ( 1 + A1 * Rs )
  Rs12 = sqrt( Rs )
  Rs32 = Rs12**3
  RsP = Rs**P
  Q1 = 2 * A * ( B1 * Rs12 + B2 * Rs + B3 * Rs32 + B4 * Rs * RsP )
  Q2 = log( 1 + 1 / Q1 )
  GG = Q0 * Q2
  Q3 = A * ( B1 / Rs12 + 2 * B2 + 3 * B3 * Rs12 + 2 * B4 * P1 * RsP)

  G_Rs = - 2 * A * A1 * Q2 - Q0 * Q3 / ( Q1**2 + Q1 )

  return
end

!***********************************************************************

! Routine of projection of the potentials coming from FLAPW Wien2k, on the FDM grid and in the atomic sphere.

subroutine potlapw(axyz,Chargat,Coupelapw,deccent,Flapw_new,Full_atom,iapot, &
            iaproto,iaprotoi,icheck,igroup,iprabs,ipr1,itabs,its_lapw,itypei,itypep,itypepr,Magnetic,mpinodes0, &
            mpirank,n_atom_0,n_atom_ind,n_atom_proto,natome,natomeq,natomp,ngroup,ngroup_lapw,nklapw,nlm_pot, &
            nlmlapwm,nmatsym,normrmt,npoint,npsom,nrato,nrato_lapw,nrm,nslapwm,nspin,ntype, &
            numat,Orthmat,Overlap,pos,Rato,rho,rhoato,rlapw,Rmtg,Rmtg0,Rmtimp,Rmtsd,Rot_int,rotloc_lapw,rs,rsato,Rsort, &
            Trace_format_wien,Trace_k,Trace_p,V0bdcFimp,Vcato,Vh,Vxc,Vxcato,Wien_file,Wien_matsym, &
            Wien_save,Wien_taulap,xyz)

  use declarations
  implicit none
  include 'mpif.h'

  integer, parameter:: ndir = 98

  integer:: i, ia, iap, icheck, idir, igr, ij1, ij2, ik, iprabs, ipr1, ir, is, isp, ispin, istat, it, itabs, j, j1, &
    j2, mpierr, mpinodes0, mpirank, n, n_atom_0, n_atom_ind, n_atom_proto, natome, natomeq, natomp, n1, n2, ndim, ngroup, &
    ngroup_lapw, nklapw, nlm_pot, nlmlapwm, nmatsym, normrmt, np, npoint, npsom, nr, nrm, ns, nslapwm, nspin, ntype, &
    Trace_k, Wien_save

  character(len=Length_name), dimension(9):: Wien_file

  complex(kind=db), dimension(nslapwm) :: taupp
  complex(kind=db), dimension(nklapw,nmatsym) :: tauk
  complex(kind=db), dimension(:), allocatable :: vcklapw
  complex(kind=db), dimension(:,:), allocatable :: rhoklapw, vxklapw

  integer, dimension(3):: kzz
  integer, dimension(0:ntype):: nlmlapw, nrato, nrato_lapw, numat
  integer, dimension(0:n_atom_proto):: iapot, itypepr
  integer, dimension(3,nslapwm):: kkk
  integer, dimension(natomp):: iaproto, igroup, itypep
  integer, dimension(natome):: iaprotoi, itypei
  integer, dimension(ngroup_lapw):: its_lapw
  integer, dimension(3,3,nslapwm):: Wien_matsym
  integer, dimension(:), allocatable:: nksym
  integer, dimension(:,:), allocatable:: llapw, mlapw
  integer, dimension(:,:,:), allocatable:: kxyz

  logical:: Coupelapw, Flapw_new, Full_atom, Magnetic, Trace_format_wien

  real(kind=db):: ctrace, f, hj1, hj2, Overlap, pd, r4pi, Rsort, Tiers, V0bdcFimp, Vht, Vol, xlength, ylength
  real(kind=db), dimension(3):: axyz, deccent, p, ptrace, v, vectrace, vx, vy, vz, w, wx, wy, wz
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(3,3):: Orthmat, Rot_int, rottem
  real(kind=db), dimension(3,ndir):: vdir
  real(kind=db), dimension(nspin):: rhot, Vxct
  real(kind=db), dimension(npoint):: rs, Vh
  real(kind=db), dimension(0:ntype):: Rmtimp
  real(kind=db), dimension(0:n_atom_proto):: Chargat, Rmtg, Rmtg0, Rmtsd, Vcmft
  real(kind=db), dimension(0:n_atom_proto,nspin):: rhomft, Vxcmft
  real(kind=db), dimension(3,nslapwm):: Wien_taulap
  real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: rsato
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nspin,n_atom_0:n_atom_ind):: rhoato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vxcato
  real(kind=db), dimension(npoint,nspin):: Vxc, rho
  real(kind=db), dimension(3,3,ngroup_lapw) :: rotloc_lapw
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(0:ntype):: rlapw
  real(kind=db), dimension(0:nrm,0:ntype):: Rato
  real(kind=db), dimension(:), allocatable :: vh_plot
  real(kind=db), dimension(:,:), allocatable :: rho_plot, vxc_plot, xxx
  real(kind=db), dimension(:,:,:), allocatable :: qxyz, vclapw
  real(kind=db), dimension(:,:,:,:), allocatable :: rholapw, vxlapw

  Vcmft(:) = 0._db
  rhomft(:,:) = 0._db
  vxcmft(:,:) = 0._db

  vectrace(1:3) = Trace_p(1:3)
  ptrace(1:3) = Trace_p(4:6)

  if( Wien_save == - 1 .and. mpirank == 0 ) then

    open(8, file = Wien_file(9), status='old', iostat=istat)
    if( istat /= 0 ) call write_open_error(Wien_file(9),istat,1)
    do i = 1,npoint
      read(8,*) Vh(i), Vxc(i,1:nspin), rho(i,1:nspin), rs(i)
    end do
    do ia = 1,n_atom_proto
      do ir = 0,nrato( itypepr(ia) )
        read(8,*) Vcato(ir,:,ia), Vxcato(ir,:,1:nspin,ia), rhoato(ir,1:nspin,ia), rsato(ir,ia)
      end do
    end do
    close(8)

  elseif( mpirank == 0 ) then

    allocate( kxyz(3,nklapw,nslapwm) )
    allocate( llapw(nlmlapwm,0:ntype) )
    allocate( mlapw(nlmlapwm,0:ntype) )
    allocate( qxyz(3,nklapw,nslapwm) )
    allocate( rhoklapw(nklapw,2*nspin-1) )
    allocate( rholapw(nrm,nlmlapwm,0:ntype,2*nspin-1) )
    allocate( vcklapw(nklapw) )
    allocate( vclapw(nrm,nlmlapwm,0:ntype) )
    allocate( vxklapw(nklapw,nspin) )
    allocate( vxlapw(nrm,nlmlapwm,0:ntype,2*nspin-1) )

    call lect_pot_lapw(flapw_new,kxyz,llapw,magnetic,mlapw,nklapw,nlmlapw,nlmlapwm, &
        nrato_lapw,nrm,nslapwm,nspin,ntype,Rato,rhoklapw,rholapw,vcklapw,vclapw,vxklapw,vxlapw,Wien_file)

    vh(:) = 0._db
    Vxc(:,:) = 0._db
    rs(:) = 100000._db
    rho(:,:) = 0._db
    ns = 1 + 2 * (nspin - 1 )

    do igr = 1,ngroup
      rottem(:,:) = rotloc_lapw(:,:,igr)
      do i = 1,3
        do j = 1,3
          rotloc_lapw(i,j,igr) = sum( rottem(i,:) * rot_int(:,j) )
        end do
      end do
    end do

    allocate( nksym(nklapw) )

    do ik = 1,nklapw
      kzz(1:3) = kxyz(1:3,ik,1)
      call stern(nslapwm,nksym(ik),nmatsym,Wien_matsym,kzz,Wien_taulap,kkk,taupp)
      do is = 1,nksym(ik)
        tauk(ik,is) = taupp(is)
        kxyz(1:3,ik,is) = kkk(1:3,is)
      end do
    end do

! Wave vector in the orthonomal internal basis
    vx(:) = orthmat(:,1)
    vy(:) = orthmat(:,2)
    vz(:) = orthmat(:,3)

! wx, wy, wz : basis of the reciprocal unit cell
    call prodvec(wx,vy,vz)

    vol = sum( wx(:) * vx(:) )
    wx(:) = wx(:) / vol
    call prodvec(wy,vz,vx)
    wy(:) = wy(:) / vol
    call prodvec(wz,vx,vy)
    wz(:) = wz(:) / vol

    do ik = 1,nklapw
      do is = 1,nksym(ik)
        qxyz(:,ik,is) = deux_pi * ( kxyz(1,ik,is) * wx(:) + kxyz(2,ik,is) * wy(:) + kxyz(3,ik,is) * wz(:) )
      end do
    end do

    call cal_vdir(ndir,vdir)
    pd = 1._db / ndir

    tiers = 1._db / 3._db
    f = 0.75_db / pi
    r4pi = 1 / sqrt( quatre_pi )

! Atomic part
    n = min(nlm_pot,nlmlapwm)

    boucle_ia: do ia = 1,n_atom_proto

      it = itypepr(ia)
      iap = iapot(ia)

      do ir = 1,nrato(it)
        if( Rato(ir,it) > rlapw(it) ) exit
        Vcato(ir,1:n,ia) = Vclapw(ir,1:n,it)
        Vxcato(ir,1:n,1:nspin,ia) = Vxlapw(ir,1:n,it,1:nspin)

        Vcato(ir,1,ia) = Vcato(ir,1,ia) * r4pi
        Vxcato(ir,1,1:nspin,ia) = Vxcato(ir,1,1:nspin,ia) * r4pi

        rhoato(ir,1:nspin,ia) = rholapw(ir,1,it,1:nspin)
        if( nspin == 1 ) then
          rhoato(ir,1,ia) = rholapw(ir,1,it,1)
        else
          do ispin = 1,nspin
            if( ispin == 1 ) then
              is = 1
            else
              is = - 1
            endif
            rhoato(ir,ispin,ia) = 0.5_db * ( rholapw(ir,1,it,1) + is * ( rholapw(ir,1,it,nspin) - rholapw(ir,1,it,ns) ) )

          end do
        endif
      end do
      rhoato(0,1:nspin,ia) = rhoato(1,1:nspin,ia)
      Vcato(0,:,ia) = Vcato(1,:,ia)
      Vxcato(0,:,1:nspin,ia) = Vxcato(1,:,1:nspin,ia)
      nr = ir
      do ir = nr,nrato(it)

        Vcato(ir,:,ia) = 0._db
        Vxcato(ir,:,1:nspin,ia) = 0._db
        rhoato(ir,1:nspin,ia) = 0._db
        do idir = 1,ndir
          p(1:3) = pos(1:3,iap) + vdir(1:3,idir) * Rato(ir,it)
          call calpot(axyz,deccent,iaproto,igroup,its_lapw,itypep,qxyz,llapw,mlapw,n_atom_proto, &
            natomeq,natomp,ngroup_lapw,nklapw,nksym,nlmlapw,nlmlapwm,nmatsym,nrato,nrm,nslapwm,nspin,ntype, &
            Orthmat,p,pos,Rato,rhoklapw,rholapw,rhomft,rhot,rlapw,Rmtg0,rotloc_lapw,tauk,vcklapw,vclapw,Vcmft,vht,vxcmft,Vxct, &
            vxklapw,vxlapw,.false.)
          Vcato(ir,1,ia) = Vcato(ir,1,ia) + pd * Vht
          do ispin = 1,nspin
            Vxcato(ir,1,ispin,ia) = Vxcato(ir,1,ispin,ia) + pd * Vxct(ispin)
          end do
          rhoato(ir,1:nspin,ia) = rhoato(ir,1:nspin,ia) + pd * rhot(1:nspin)
        end do

      end do

      do ir = 1,nrato(it)
        do ispin = 1,nspin
          rhoato(ir,ispin,ia) = max( eps10, rhoato(ir,ispin,ia) )
        end do
        rsato(ir,ia) = ( f / sum( rhoato(ir,1:nspin,ia) ) )**tiers
      end do
      rsato(0,ia) = rsato(1,ia)

    end do boucle_ia

  endif

  if( mpinodes0 > 1 ) then
    ndim = ( nrm + 1 ) * ( n_atom_ind - n_atom_0 + 1 )

    call MPI_Bcast(rsato,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    ndim = ndim*nlm_pot
    call MPI_Bcast(Vcato,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    ndim = ( nrm + 1 ) * ( n_atom_ind - n_atom_0 + 1 ) * nspin

    call MPI_Bcast(rhoato,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    ndim = ndim*nlm_pot
    call MPI_Bcast(Vxcato,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  endif

  call Raymuf(.true.,Chargat,Full_atom,1,iapot,iaproto,iaprotoi,icheck,iprabs, &
        ipr1,itabs,itypei,itypep,itypepr,n_atom_0,n_atom_ind,n_atom_proto,natome,natomeq,natomp,nlm_pot, &
        normrmt,nrato,nrm,nspin,ntype,numat,Overlap,pos,Rato,rhoato,Rmtg,Rmtg0,Rmtimp, &
        Rmtsd,Rsort,V0bdcFimp,Vcato,Vxcato)

  call Potrmt(Full_atom,iaprotoi,ipr1,iprabs,itypepr,n_atom_0,n_atom_ind,n_atom_proto,natome,nlm_pot,nrato,nrm,nspin,ntype, &
                  numat,Rato,rhoato,rhomft,Rmtg0,Vcato,Vcmft,Vxcato,Vxcmft,iapot)

  if( Wien_save /= - 1 .and. mpirank == 0 ) then

    if( icheck > 0 ) write(3,110)

    do   ! Loop for the cuting

      if( coupelapw ) then
        select case(Trace_k)
          case(1)
            n1 = 501
            n2 = 1
            np = n1 * n2
            allocate( xxx(np,3) )
            w(1:3) = axyz(1:3) * vectrace(1:3) / ( np - 1 )
            v(1:3) = axyz(1:3) * ptrace(1:3)
            do i = 1,np
              xxx(i,1:3) = ( i - 1 ) * w(1:3) + v(1:3)
            end do
          case(2)
            n1 = 501
            n2 = 21
            np = n1 * n2
            ctrace = ptrace(1) / bohr
            allocate( xxx(np,3) )
            do j = 3,1,-1
              if( abs( vectrace(j) ) < eps4 ) cycle
              if( j == 3 ) then
                j1 = 1
                j2 = 2
              elseif( j == 2 ) then
                j1 = 1
                j2 = 3
              else
                j1 = 2
                j2 = 3
              endif
              hj1 = axyz(j1) / ( n1 - 1 )
              hj2 = axyz(j2) / ( n2 - 1 )
              i = 0
              do ij2 = 1,n2
                do ij1 = 1,n1
                  i = i + 1
                  xxx(i,j1) = ( ij1 - 1 ) * hj1
                  xxx(i,j2) = ( ij2 - 1 ) * hj2
                  xxx(i,j) = ( ctrace - vectrace(j1) * xxx(i,j1) - vectrace(j2) * xxx(i,j2) ) / vectrace(j)
                end do
              end do
            end do
            np = i
        end select
      else
        np = npoint
        allocate( xxx(np,3) )
        do i = 1,npoint
          xxx(i,1:3) = xyz(1:3,i)
        end do
      endif

      if( Coupelapw ) then
        if( trace_format_wien ) then
          p(:) = xxx(n1,:) - xxx(1,:)
          ylength = sqrt( sum( p(:)**2 ) )
          p(:) = xxx(np,:) - xxx(n1,:)
          xlength = sqrt( sum( p(:)**2 ) )
          allocate( vh_plot(np) )
          allocate( vxc_plot(np,nspin) )
          allocate( rho_plot(np,nspin) )
        endif
        if( nspin == 1 ) then
          write(3,112)
        else
          write(3,113)
        endif
      endif

      boucle_point: do i = 1,np

        p(:) = xxx(i,:)
        call calpot(axyz,deccent,iaproto,igroup,its_lapw,itypep,qxyz,llapw,mlapw,n_atom_proto, &
            natomeq,natomp,ngroup_lapw,nklapw,nksym,nlmlapw,nlmlapwm,nmatsym,nrato,nrm,nslapwm,nspin,ntype, &
            Orthmat,p,pos,Rato,rhoklapw,rholapw,rhomft,rhot,rlapw,Rmtg0,rotloc_lapw,tauk,vcklapw,vclapw,Vcmft,vht,Vxcmft,Vxct, &
            vxklapw,vxlapw,.true.)

        if( coupelapw ) then
          write(3,115) xxx(i,1:3)*bohr, vht*rydb, Vxct(1:nspin)*rydb, rhot(1:nspin)
          if( trace_format_wien ) then
            vh_plot(i) = vht
            vxc_plot(i,1:nspin) = Vxct(1:nspin)
            rho_plot(i,1:nspin) = rhot(1:nspin)
          endif
        else
          vh(i) = vht
          Vxc(i,1:nspin) = Vxct(1:nspin)
          rho(i,1:nspin) = rhot(1:nspin)
        endif

      end do boucle_point     ! End of loop on the points

      if( coupelapw .and. trace_format_wien ) then
        write(3,116) n2, n1, xlength, ylength
        write(3,117) ( vh_plot(i), i = 1,np )
        do isp = 1,nspin
          write(3,118) isp, n2, n1, xlength, ylength
          write(3,117) ( vxc_plot(i,isp), i = 1,np )
          write(3,119) isp, n2, n1, xlength, ylength
          write(3,117) ( rho_plot(i,isp), i = 1,np )
        end do
        deallocate( vh_plot )
        deallocate( vxc_plot )
        deallocate( rho_plot )
      endif

      deallocate( xxx )

      if( Coupelapw ) then
        Coupelapw = .false.
      else
        exit
      endif

    end do  ! End of loop Coupe

    deallocate( nksym )

! Calculation of the Fermi radius, rs

    do ispin = 1,nspin
      do i = 1,npoint
        rho(i,ispin) = max( rho(i,ispin), eps10 )
      end do
    end do

    tiers = 1._db / 3._db
    f = 0.75_db / pi

    do i = 1,npoint
      rs(i) = ( f / sum( rho(i,1:nspin) ) )**tiers
    end do

    deallocate( kxyz, llapw, mlapw, qxyz, rhoklapw, rholapw, vcklapw, vclapw, vxklapw, vxlapw )

  endif   ! arriving Wien_save == - 1 .or. mpirank /= 0

  if( mpinodes0 > 1 ) then

    call MPI_Bcast(vh,npoint,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(rs,npoint,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    ndim = npoint * nspin

    call MPI_Bcast(rho,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Vxc,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  endif

  if( icheck > 1 ) then
    if( nspin == 1 ) then
      write(3,120)
    else
      write(3,125)
    endif
    do i = 1,npoint
      write(3,130) i, vh(i)*rydb, Vxc(i,1:nspin)*rydb, rho(i,1:nspin), rs(i)
    end do
    do  ia = 1,n_atom_proto
      if( nspin == 1 ) then
        write(3,140) ia
      else
        write(3,150) ia
      endif
      it = itypepr(ia)
      do ir = 1,nrato(it)
        write(3,160) Rato(ir,it)*bohr, quatre_pi * rhoato(ir,1:nspin,ia) * Rato(ir,it)**2, &
          rsato(ir,ia), ( Vxcato(ir,i,1:nspin,ia)*rydb, Vcato(ir,i,ia)*rydb, i = 1,nlm_pot )
      end do
    end do
  endif

  if( Wien_save == 1 .and. mpirank == 0 ) then
    open(8,file = Wien_file(9))
    do i = 1,npoint
      write(8,*) vh(i), Vxc(i,1:nspin), rho(i,1:nspin), rs(i)
    end do
    do ia = 1,n_atom_proto
      do ir = 0,nrato( itypepr(ia) )
        write(8,*) Vcato(ir,:,ia), Vxcato(ir,:,1:nspin,ia), rhoato(ir,1:nspin,ia), rsato(ir,ia)
      end do
    end do
    close(8)
  endif

  return
  110 format(/' ---- Potlapw ------',100('-'))
  112 format(/6x,'x',12x,'y',12x,'z',10x,'vh(eV)',6x,'Vxc(eV)',8x,'rho')
  113 format(/6x,'x',12x,'y',12x,'z',9x,'vh(eV)      Vxc(up)',9x, 'Vxc(dn)',9x,'rho(up)',8x,'rho(dn)')
  115 format(1p,9e13.5)
  116 format(/' Format WIEN :',//, ' Potential coulombien (rydb) :'/,2i5,2f10.5)
  117 format(5e16.8)
  118 format(/' Potential d echange (rydb), ispin =',i2,' :'/, 2i5,2f10.5)
  119 format(/' Densite electronique, ispin =',i2,' :'/,2i5,2f10.5)
  120 format(/4x,'i     vh(eV)      Vxc(eV)',8x,'rho           rs(ua)')
  125 format(/4x,'i     vh(eV)    Vxc(up)(eV) Vxc(down)(eV)   rho(up)     rho(down)       rs(ua)')
  130 format(i5,1p,6e13.5)
  140 format(/'    Rato  4*pi*r2*rhoato   rsato     vxcato      vato    ia =',i3)
  150 format(/'    Rato    rhoato(up) rhoato(down)  rsato   vxcato(up) vxcato(down)   vato    ia =',i3)
  160 format(1p,70e11.3)
end

!***********************************************************************

! Calcul des vecteurs directions pour le calcul du potentiel moyen radial

subroutine cal_vdir(ndir,vdir)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db), dimension(3,ndir):: vdir

  r3 = 1 / sqrt( 3._db )
  rac2 = 1 / sqrt( 2._db )
  vdir(:,1:3) = 0._db
  vdir(:,4:7) = r3
  vdir(:,8:13) = rac2
  do idir = 1,3
    vdir(idir,idir) = 1._db
  end do
  do idir = 5,7
    k = idir - 4
    vdir(k,idir) = - vdir(k,idir)
  end do
  do idir = 8,10
    k = idir - 7
    vdir(k,idir) = 0._db
    jdir = idir + 3
    vdir(k,jdir) = 0._db
    k = mod(k,3) + 1
    vdir(k,jdir) = - vdir(k,jdir)
  end do
  dr5 = 2 / sqrt(5._db)
  ur5 = 1 / sqrt(5._db)
  vdir(1,14) = dr5;  vdir(2,14) = ur5;  vdir(3,14) = 0._db
  vdir(1,15) = ur5;  vdir(2,15) = dr5;  vdir(3,15) = 0._db
  vdir(1,16) = -dr5;  vdir(2,16) = ur5;  vdir(3,16) = 0._db
  vdir(1,17) = -ur5;  vdir(2,17) = dr5;  vdir(3,17) = 0._db
  do idir = 18,21
    vdir(1,idir) = vdir(3,idir-4)
    vdir(2,idir) = vdir(1,idir-4)
    vdir(3,idir) = vdir(2,idir-4)
  end do
  do idir = 22,25
    vdir(1,idir) = vdir(3,idir-4)
    vdir(2,idir) = vdir(1,idir-4)
    vdir(3,idir) = vdir(2,idir-4)
  end do
  dr6 = 2 / sqrt(6._db)
  ur6 = 1 / sqrt(6._db)
  vdir(1,26) = dr6;  vdir(2,26) = ur6;  vdir(3,26) = ur6
  vdir(1,27) = ur6;  vdir(2,27) = dr6;  vdir(3,27) = ur6
  vdir(1,28) = ur6;  vdir(2,28) = ur6;  vdir(3,28) = dr6
  dr9 = 2 / 3._db
  ur9 = 1 / 3._db
  vdir(1,29) = ur9;  vdir(2,29) = dr9;  vdir(3,29) = dr9
  vdir(1,30) = dr9;  vdir(2,30) = ur9;  vdir(3,30) = dr9
  vdir(1,31) = dr9;  vdir(2,31) = dr9;  vdir(3,31) = ur9
  do idir = 26,31
    vdir(1,idir+6) = -vdir(1,idir)
    vdir(2,idir+6) = vdir(2,idir)
    vdir(3,idir+6) = vdir(3,idir)
    vdir(1,idir+12) = vdir(1,idir)
    vdir(2,idir+12) = -vdir(2,idir)
    vdir(3,idir+12) = vdir(3,idir)
    vdir(1,idir+18) = vdir(1,idir)
    vdir(2,idir+18) = vdir(2,idir)
    vdir(3,idir+18) = -vdir(3,idir)
  end do
  n = ndir / 2
  do idir = n+1,ndir
    vdir(:,idir) = - vdir(:,idir-n)
  end do

  return
end

!***********************************************************************

! Reading routine for the potential and electronic densities coming from WIEN2k

subroutine lect_pot_lapw(flapw_new,kxyz,llapw,magnetic,mlapw,nklapw,nlmlapw,nlmlapwm, &
        nrato_lapw,nrm,nslapwm,nspin,ntype,Rato,rhoklapw,rholapw,vcklapw,vclapw,vxklapw,vxlapw,Wien_file)

!     nlmlapwm : nombre max de termes (l,m)
!     nslapwm : nombre max d'operations de symetrie ponctuelle
!     nklapw : nombre d'ondes planes

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer:: eof, nslapwm

  character(len=Length_name):: nomvcoul, nomr2v, nomr2vdn, nomclm(2*nspin-1)
  character(len=Length_name), dimension(9):: Wien_file

! modif delphine 8/06/01 pour les structures sans centre d'inversion
! one then has to read the imaginary part of Fourrier coefficients
  complex(kind=db), dimension(nklapw) :: vcklapw
  complex(kind=db), dimension(nklapw,nspin) :: vxklapw
  complex(kind=db), dimension(nklapw,2*nspin-1) :: rhoklapw

  integer, dimension(0:ntype):: nlmlapw, nrato_lapw
  integer, dimension(nlmlapwm,0:ntype):: llapw, mlapw
  integer, dimension(3,nklapw,nslapwm):: kxyz

  logical flapw_new, magnetic

  real(kind=db), dimension(nrm):: r2
  real(kind=db), dimension(nrm,nlmlapwm,0:ntype):: vclapw
  real(kind=db), dimension(nrm,nlmlapwm,0:ntype,2*nspin-1):: rholapw, vxlapw
  real(kind=db), dimension(0:nrm,0:ntype):: Rato

! kxyz           = vecteurs du reseau reciproque
! ll             = nombre de (l,m) par atome
! ntype          = nombre d'atomes inequivalents
! nmatsym        = nombre d'op. de symetrie

! Lecture du potentiel coulombien (case.vcoul), du potentiel d'echange
! et de correlation (case.r2v) et de la densite de charge (case.clmsum)

  nomvcoul = Wien_file(2)
  nomr2v = Wien_file(3)
  nomr2vdn = Wien_file(4)
  nomclm(1) = Wien_file(5)
  if( nspin == 2 ) nomclm(2:3) = Wien_file(6:7)

! 1- atomic spheres

  if( magnetic .and. flapw_new ) then
    nfich = 6
  elseif( magnetic ) then
    nfich = 5
  else
    nfich = 3
  endif

  do ifich = 1,nfich
    jfich = ifich
    if( magnetic .and. flapw_new ) then
      if( ifich == 2 ) then
        isp = 1
      elseif( ifich == 3 ) then
        jfich = 6
        isp = 2
      else
        isp = ifich - 3
      endif
    else
      isp = ifich - 2
    endif
    select case(jfich)
      case(1)
        open(8, file = nomvcoul, status='old', iostat=istat)
        if( istat /= 0 ) call write_open_error(nomvcoul,istat,1)
      case(2)
        open(8, file = nomr2v, status='old', iostat=istat)
        if( istat /= 0 ) call write_open_error(nomr2v,istat,1)
      case(3,4,5)
        open(8, file = nomclm(isp), status='old', iostat=istat)
        if( istat /= 0 ) call write_open_error(nomclm(isp),istat,1)
      case(6)
        open(8, file = nomr2vdn, status='old', iostat=istat)
        if( istat /= 0 ) call write_open_error(nomr2vdn,istat,1)
    end select

    if( .not. (magnetic .and. Flapw_new) .and. ifich == 2 ) then
      ns = nspin
    else
      ns = 1
    endif

    do it = 1,ntype

      do i = 1,3
        read(8,*)
      end do
      read(8,'(16x,i2)') natom
      read(8,'(16x,i2)',iostat=eof) ll
      if( eof /= 0 ) then
        backspace(8)
        read(8,'(17x,i2)') ll
      endif

      nlmlapw(it) = ll

      do ispin = 1,ns
        do l = 1,ll
          read(8,*)
          read(8,*)
          read(8,'(16x,i2,5x,i2)') llapw(l,it), mlapw(l,it)
          read(8,*)

          do j = 1,nrato_lapw(it),4
            jcard = min(nrato_lapw(it),j+3)
            select case(jfich)
              case(1)
                read(8,'(3X,4E19.12)') vclapw(j:jcard,l,it)
              case(2)
                read(8,'(3X,4E19.12)') vxlapw(j:jcard,l,it,ispin)
              case(3,4,5)
                read(8,'(3X,4E19.12)') rholapw(j:jcard,l,it,ispin)
              case(6)
                read(8,'(3X,4E19.12)') vxlapw(j:jcard,l,it,nspin)
            end select
          end do
        end do
      end do

      do i = 1,3
        read(8,*)
      end do

    end do ! end of loop over atoms

! 2- interstitial reegion

    do i = 1,3
      read(8,*)
    end do

    do ispin = 1,ns
      read(8,*)
      read(8,*)
      read(8,*)
      do ik = 1,nklapw
        select case(jfich)
          case(1)
            read(8,'(3X,3I5,2E19.12)') kxyz(1:3,ik,1), cfr, cfi
            vcklapw(ik) = cmplx( cfr, cfi,db )
          case(2)
            read(8,'(3X,3I5,2E19.12)') kxyz(1:3,ik,1), cfr, cfi
            vxklapw(ik,ispin) = cmplx( cfr, cfi,db )
          case(3,4,5)
            read(8,'(3X,3I5,2E19.12)') kxyz(1:3,ik,1), cfr, cfi
            rhoklapw(ik,isp) = cmplx( cfr, cfi,db )
          case(6)
            read(8,'(3X,3I5,2E19.12)') kxyz(1:3,ik,1), cfr, cfi
            vxklapw(ik,nspin) = cmplx( cfr, cfi,db )
        end select
      end do

    end do

    close(8)
  end do

! Division par r2 et division par 4pi pour le terme lm=(0,0) de rholapw

  ns = 1 + 2 * ( nspin - 1 )
  do it = 1,ntype
    r2(1:nrato_lapw(it)) = Rato(1:nrato_lapw(it),it)**2
    do j = 1,nrato_lapw(it)
      do l = 1,nlmlapw(it)
        vclapw(j,l,it) = vclapw(j,l,it) / r2(j)
        vxlapw(j,l,it,1:nspin) = vxlapw(j,l,it,1:nspin) / r2(j)
        rholapw(j,l,it,1:ns) = rholapw(j,l,it,1:ns) / r2(j)
      enddo
      rholapw(j,1,it,1:ns) = rholapw(j,1,it,1:ns) / quatre_pi
    enddo
  enddo

  return
end

!***********************************************************************

subroutine calpot(axyz,deccent,iaproto,igroup,its_lapw,itypep,qxyz,llapw,mlapw,n_atom_proto, &
            natomeq,natomp,ngroup_lapw,nklapw,nksym,nlmlapw,nlmlapwm,nmatsym,nrato,nrm,nslapwm,nspin,ntype, &
            Orthmat,p,pos,Rato,rhoklapw,rholapw,rhomft,rhot,rlapw,Rmtg0,rotloc_lapw,tauk,vcklapw,vclapw,Vcmft,vht,vxcmft,Vxct, &
            vxklapw,vxlapw,centat)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer:: nslapwm

  complex(kind=db):: coef, imag1, rhoc, tauk(nklapw,nmatsym)
! modif delphine 8/06/01 pour les structures sans centre d'inversion
! il faut alors absolument lire la partie imaginaire des coeff de Fourier.
  complex(kind=db), dimension(nklapw):: vcklapw
  complex(kind=db), dimension(nklapw,nspin):: vxklapw
  complex(kind=db), dimension(nklapw,2*nspin-1):: rhoklapw
! fin modif delphine
  complex(kind=db), dimension(:), allocatable:: ylm

  integer, dimension(natomp):: iaproto, igroup, itypep
  integer, dimension(0:ntype):: nlmlapw, nrato
  integer, dimension(ngroup_lapw):: its_lapw
  integer, dimension(nklapw):: nksym
  integer, dimension(nlmlapwm,0:ntype):: llapw, mlapw

  logical:: centat

  real(kind=db), dimension(2):: vc
  real(kind=db), dimension(3):: axyz, deccent, p, v, w
  real(kind=db), dimension(2,3):: rh(2,3)
  real(kind=db), dimension(3,3):: Orthmat
  real(kind=db), dimension(nspin):: rhot, Vxct
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(3,3,ngroup_lapw) :: rotloc_lapw
  real(kind=db), dimension(3,nklapw,nslapwm):: qxyz
  real(kind=db), dimension(nrm,nlmlapwm,0:ntype) :: vclapw
  real(kind=db), dimension(nrm,nlmlapwm,0:ntype,2*nspin-1) :: rholapw, vxlapw
  real(kind=db), dimension(0:n_atom_proto):: Rmtg0, Vcmft
  real(kind=db), dimension(0:n_atom_proto,nspin):: rhomft, vxcmft
  real(kind=db), dimension(0:ntype):: rlapw
  real(kind=db), dimension(0:nrm,0:ntype):: Rato

  vht = 0._db
  Vxct(1:nspin) = 0._db
  rhot(1:nspin) = 0._db
  rac_2 = sqrt( 2._db )
  sqin2 = 1._db / rac_2
  ns = 1 + 2 * (nspin - 1 )

! Dans le premier cas, c'est pour calculer le potentiel atomique
! pour un rayon grand, il faut eviter la divergence si on s'approche
! d'un autre atome.
  if( centat ) then

! On regarde si on est dans un atome
  do ia = 1,natomp

    ipr = iaproto(ia)
    it = itypep(ia)
    igr = igroup(ia)
    v(1:3) = p(1:3) - pos(1:3,ia)
    dist = sqrt( sum( v(:)**2 ) )
    dist = max(dist,Rato(1,it))

! Si on tombe dans un atome chevauchant la frontiere exterieur,
! on prend le potentiel au niveau de son rayon muffin-tin.
! On ne peut pas le faire au debut car Rmtg0 n'est pas encore calcule.
    if( ia > natomeq .and. dist < Rmtg0(ipr) .and. centat ) then
      vht = Vcmft(ipr)
      Vxct(1:nspin) = vxcmft(ipr,1:nspin)
      rhot(1:nspin) = rhomft(ipr,1:nspin)
      return
    endif

    if( dist >= rlapw(it) ) cycle
! On est dans un atome FLAPW.

    do ir = 2,nrato(it)
      if(Rato(ir,it) > dist) exit
    end do
    p1 = ( dist - Rato(ir-1,it) ) / ( Rato(ir,it)-Rato(ir-1,it) )
    p2 = 1 - p1

    do j = 1,3
      w(j) = sum( rotloc_lapw(j,:,igr) * v(:) )
    end do

    l = abs( llapw(nlmlapw(it),it) )
    nlm = ( l + 1 )**2
    allocate( ylm(nlm) )
    call ylmlapw(l,nlm,w,ylm)

    do lm = 1,nlmlapw(it)
      l = abs( llapw(lm,it) )
      m = mlapw(lm,it)
      lm1 = l*(l+1) + 1 + abs(m)
      minu = 1
      imag1 = (1._db,0._db)
      if( llapw(lm,it) < 0 ) then
        imag1 = - img
        minu = - 1
      endif
      if( mod(m,2) == 1 ) then
        imag1 = - imag1
        minu = - minu
      endif
      idx = l*(l+1) + 1
      if( m == 0 ) then
        ylmr = real( ylm(idx),db )
      else
        idp = idx + m
        idm = idx - m
        ylmr = sqin2 * real( (ylm(idp) + minu * ylm(idm))*imag1,db )
      endif
      if( l == 0 ) then
        ylm1 = 1._db
      else
        ylm1 = ylmr
      endif

      if( its_lapw(igr) > 0  .and. lm > 1 .and. lm < 6 ) then
        if( lm == 2 .or. lm == 4 ) then
          lm1 = lm
          lm2 = lm + 1
        else
          lm1 = lm - 1
          lm2 = lm
        endif
      endif

      do k = ir-1,ir
        j = k - ir + 2
        if( its_lapw(igr) > 0 .and. lm > 1 .and. lm < 6 ) then
          call cascfc(vc(j), vclapw(k,lm1,it),vclapw(k,lm2,it),lm)
        else
          vc(j) = vclapw(k,lm,it)
        endif
      end do
      vht = vht + ( p1 * vc(2) + p2 * vc(1) ) * ylmr

      do isp = 1,nspin
        do k = ir-1,ir
          j = k - ir + 2
          if( its_lapw(igr) > 0  .and. lm > 1 .and. lm < 6 ) then
            call cascfc(vc(j),vxlapw(k,lm1,it,isp), vxlapw(k,lm2,it,isp),lm)
          else
            vc(j) = vxlapw(k,lm,it,isp)
          endif
        end do
        Vxct(isp) = Vxct(isp) + ( p1 * vc(2) + p2 * vc(1) ) * ylmr
      end do

      do isp = 1,ns
        do k = ir-1,ir
          j = k - ir + 2
          if( its_lapw(igr) > 0 .and. lm > 1 .and. lm < 6 ) then
            call cascfc(rh(j,isp),rholapw(k,lm1,it,isp), rholapw(k,lm2,it,isp),lm)
          else
            rh(j,isp) = rholapw(k,lm,it,isp)
          endif
        end do
      end do

      if( nspin == 1 ) then
        rhot(1) = rhot(1) + ( p1 * rh(2,1) + p2 * rh(1,1) ) * ylm1
      else
        do ispin = 1,nspin
          if( ispin == 1 ) then
            rho1 = 0.5_db * ( rh(2,1) + rh(2,nspin) - rh(2,ns) )
            rho2 = 0.5_db * ( rh(1,1) + rh(1,nspin) - rh(1,ns) )
          else
            rho1 = 0.5_db * ( rh(2,1) - rh(2,nspin) + rh(2,ns) )
            rho2 = 0.5_db * ( rh(1,1) - rh(1,nspin) + rh(1,ns) )
          endif
          rhot(ispin) = rhot(ispin) + ( p1*rho1 + p2*rho2 ) * ylm1
        end do
      endif

    end do  ! fin de la boucle sur les harmoniques

    deallocate( ylm )

    return
  end do

! On est dans la zone interstitielle FLAPW

  endif  ! fin test centat

  v(1:3) = deccent(1:3) / axyz(1:3)
  v = matmul( orthmat, v )

  v(1:3) = p(1:3) + v(1:3)

  do ik = 1,nklapw
    coef = (0._db,0._db)
    do is = 1,nksym(ik)
      arg = sum( qxyz(:,ik,is) * v(:) )
      coef = coef + cmplx( cos(arg), sin(arg),db ) * tauk(ik,is)
    end do
    coef  = coef / nksym(ik)

    vht = vht + real( vcklapw(ik) * coef, db )

    Vxct(1:nspin) = Vxct(1:nspin) + real( vxklapw(ik,1:nspin) * coef, db )

    if( nspin == 1 ) then
      rhot(1) = rhot(1) + real( rhoklapw(ik,1) * coef, db )
    else
      do ispin = 1,nspin
        if( ispin == 1 ) then
          is = 1
        else
          is = - 1
        endif
        rhoc = 0.5_db * ( rhoklapw(ik,1) + is * ( rhoklapw(ik,nspin) - rhoklapw(ik,ns) ) )
        rhot(ispin) = rhot(ispin) + real( rhoc * coef,db )
      end do
    endif
  end do

  return
end

!***********************************************************************

! Routine adaptee de la routine sul dans lapw5

subroutine cascfc(vo,v1,v2,lm)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db), dimension(2) :: cf, vi

!  NOTE THIS SUBROUTINE IS ONLY FOR CUBIC FCC, NACL, CSCL OR
! ANTIFLUORIT - STRUCTUR
! K0(R)=1/SQRT(4PI)
! K4(R)=SQRT(7/12)Y40(R) + SQRT(5/24)*(Y44(R) + Y4-4(R))
! K6(R)=SQRT(2)/4*Y60(R) - SQRT(7)/4* (Y64(R) + Y6-4(R))
! K7(R)=(-I)*(Y32(R)-Y3-2(R))  FOR ME IN ANFL STRUCTUR
! CSO=SQRT(2)/4
! CFO=SQRT(7/12)
! CFF=SQRT(5/12)
! CSO=-SQRT(14)/4

  vi(1) = v1
  vi(2) = v2
  if( lm == 2 .or. lm == 3 ) then
    cf(1) = .763762616
    cf(2) = .645497224
  else
    cf(1) = .35355339
    cf(2) =-.93541435
  endif
  i = mod(lm,2) + 1
  vo = cf(i) * sum( cf(:) * vi(:) )

  return
end

!***********************************************************************

subroutine ylmlapw(lomax,nlm,v,y)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db) y(nlm)
  real(kind=db) v(3), p(8,8)

! Calcul des sinus et cosinus du vecteur v

  fpi = 4._db * acos(-1.d0)
  xy = v(1)**2 + v(2)**2
  xyz = xy + v(3)**2
  if( xy > eps10 ) then
    xy = sqrt(xy)
    xyz = sqrt(xyz)
    cth = v(3) / xyz
    sth = xy / xyz
    cfi = v(1) / xy
    sfi = v(2) / xy
  else
    cth = 1._db
    if( v(3) < 0._db )  cth = - cth
    sth = 0._db
    cfi = 1._db
    sfi = 0._db
  endif

  yr = 1 / sqrt(fpi)
  y(1) = cmplx(yr,0._db,db)

  i = 1
  p(1,1) = 1._db
  p(2,1) = cth
  c2l = cth
  tcth = cth + cth

  do l = 2,lomax+1
    i = i + l - 1
    idwn = i + 2
    l1 = l + 1
    l2 = l - 1
    lm = l2
    lm2 = l
    cmfi = 1._db
    smfi = 0._db
    cd = 1._db
    c2l = c2l + tcth
    sgnm = 1._db

    do m = 1,l

      if( m < l ) then
        m1 = m + 1
        p(l1,m) = ( c2l*p(l,m) - lm*p(l2,m) ) / lm2
        c1l = ( lm + 1 ) * cth
        if( abs( sth ) < eps10 ) then
          p(l,m1) = 0._db
        else
          p(l,m1) = ( c1l*p(l,m) - lm2*p(l1,m) ) / sth
        endif
      endif

      i = i + 1
      idwn = idwn - 1
      csr = sqrt( (2*l-1._db) / (fpi*cd) )
      cyp = sgnm * csr * p(l,m)
      yr = cyp * cmfi
      yi = cyp * smfi
      y(i) = cmplx(yr,yi,db)
      if( idwn /= i ) y(idwn) = sgnm * cmplx(yr,-yi,db)

      cn = cmfi
      cmfi = cn*cfi - sfi*smfi
      smfi = sfi*cn + smfi*cfi
      lm2 = lm2-1
      lm = lm+1
      cd = cd*lm*lm2
      sgnm = - sgnm

    end do

  end do

  return
end

!***********************************************************************

! Stern generates the star of rec lattice vector kzz(i).
!       The star vectors are stored in kkk, the star-size in nst,
!       imat contains the symmetry-matrices.

subroutine stern(nslapwm,nst,iord,imat,kzz,Wien_taulap,kkk,taupp)

  use declarations
  implicit real(kind=db) (a-h,o-z)

!        arguments

  integer:: iord, nslapwm, nst
  integer:: kzz(3), kkk(3,nslapwm), imat(3,3,nslapwm)

  real(kind=db):: Wien_taulap(3,nslapwm)

  complex(kind=db):: taupp(nslapwm)

!        locals

  integer:: i, j, k, l, m
  integer:: g(3), indexx(nslapwm)

  real(kind=db):: tpi, tk

!        intrinsic functions

  intrinsic atan, exp

  tpi = 8.0d+0*atan(1.0d+0)
  g(1:3) = kzz(1:3)
  nst = 0

!         start loop over all symmetry operations

  boucle_ext: do i = 1, iord
     tk = 0._db
     do j = 1, 3
       tk = tk + Wien_taulap(j,i)*g(j)*tpi
       k = 0
       do l = 1, 3
         k = imat(j,l,i)*g(l) + k
       end do
       kkk(j,i) = k
     end do
     if( nst /= 0 ) then

!        proof, if the vector kkk(j,i) is a new starmember or not

       boucle_m: do m = 1, nst
         do j = 1,3
           if( kkk(j,m) /= kkk(j,i) ) cycle boucle_m
         end do

!        kkk(j,i) is not a new starmember, it is equiv to kkk(j,m).
!        but the tauphase of kkk(j,i) can be new.  therefore the
!        already defined phase taupp(m) is averaged with the phase
!        of kkk(j,m).

         taupp(m) = taupp(m) + exp( img * tk )
         indexx(m) = indexx(m) + 1
         cycle boucle_ext
       end do boucle_m
     endif

!        new vector found

     nst = nst + 1
     do j = 1,3
       kkk(j,nst) = kkk(j,i)
     end do
     taupp(nst) = exp( img * tk )
     indexx(nst) = 1
  end do boucle_ext

  do i = 1,nst
    taupp(i) = taupp(i) / indexx(i)
  end do

  return
end

!***********************************************************************

! Complementary calculations on the potential
! Calculation of the average interstitial potential and of the the maximum kinetic energy

subroutine Potential_comp(Bulk_atom_done,Cal_xanes,distai,dV0bdcF,Full_atom, &
            Green,i_range,i_self,iaabs,iapot,iaproto,iaprotoi,icheck,imoy,imoy_out,iopsymr,ipr1,isrt,itypepr,korigimp,Magnetic, &
            Moy_loc,mpirank,n_atom_0,n_atom_ind,n_atom_proto,n_atom_proto_bulk,natome,natomp,nim,nlm_pot,npoint,npsom,nptmoy, &
            nptmoy_out,nrato,nrm,nsortf,nspin,nstm,ntype,numat,poidsov,poidsov_out,pos, &
            Rato,Rchimp,rhoato,Rmtg,Rmtg0,Rmtsd,rs,rsbdc,rsbdc_out,Rsort,rvol,SCF,Sym_2D,V_intmax,V0bdcF,V0bdcF_out,V0bdcFimp, &
            V0muf,Vcato,Vh,Vhbdc,Vhbdc_init,Vhbdc_out,Vr,Vxc,Vxcato,VxcbdcF,VxcbdcF_out,xyz,Workf)

  use declarations
  implicit none

  integer:: i_range, i_self, iaabs, icheck, icheck_loc, ipr1, ispin, mpirank, n_atom_0, n_atom_ind, n_atom_proto, &
            n_atom_proto_bulk, natome, natomp, nim, nlm_pot, npoint, npsom, nptmoy, nptmoy_out, nrm, nsortf, nspin, nstm, ntype

  integer, dimension(natomp):: iaproto
  integer, dimension(natome):: iaprotoi
  integer, dimension(npoint):: imoy, imoy_out
  integer, dimension(nstm):: isrt
  integer, dimension(nopsm):: iopsymr
  integer, dimension(0:n_atom_proto):: iapot, itypepr
  integer, dimension(0:ntype):: nrato, numat

  logical:: Bulk_atom_done, Cal_xanes, Full_atom, Green, korigimp, Magnetic, Moy_loc, SCF, Sym_2D

  real(kind=db):: distai, rsbdc, rsbdc_out, Rsort, V_intmax, V0muf, Vhbdc, &
    Vhbdc_init, Vhbdc_out, Workf

  real(kind=db), dimension(nspin):: dV0bdcF, V0bdcF, V0bdcFimp, V0bdcF_out, VmoyF, VmoyF_out, VxcbdcF, VxcbdcF_out
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(0:n_atom_proto):: Rchrg, Rmtg, Rmtg0, Rmtsd
  real(kind=db), dimension(0:ntype):: Rchimp
  real(kind=db), dimension(nim):: rvol
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(npoint):: poidsov, poidsov_out, rs, Vh
  real(kind=db), dimension(npoint,nspin):: Vxc, Vr
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,0:ntype):: Rato
  real(kind=db), dimension(0:nrm,nspin,n_atom_0:n_atom_ind):: rhoato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato

  if( i_self == 1 .or. Cal_xanes ) then
    icheck_loc = icheck
  else
    icheck_loc = icheck - 1
  endif

  if( icheck_loc > 0 ) write(3,110)

  call Ptmoy(distai,Green,iaabs,iaproto,icheck,imoy,imoy_out,iopsymr,isrt,Moy_loc, &
             n_atom_proto,natomp,nim,npoint,nsortf,npsom,nptmoy,nptmoy_out,nstm,poidsov,poidsov_out, &
             pos,Rmtg0,Rsort,rvol,Sym_2D,xyz)

  if( npoint > 0 ) then
    do ispin = 1,nspin
      Vr(1:npoint,ispin) = Vh(1:npoint) + Vxc(1:npoint,ispin)
    end do
  endif

  call Cal_Vmoy(Cal_xanes,i_range,i_self,icheck,imoy,imoy_out,korigimp,Magnetic,mpirank,npoint,nptmoy,nptmoy_out,nspin, &
             poidsov,poidsov_out,rs,rsbdc,rsbdc_out,SCF,V0bdcFimp,Vh,Vhbdc,Vhbdc_out,VmoyF,VmoyF_out,Vr,VxcbdcF,VxcbdcF_out)

! One imposes to the interstial coulomb potential to remain the same to avoid Madelung constant potential problem.
  if( i_self == 1 ) then
    Vhbdc_init = Vhbdc
  elseif( SCF ) then
    call Shift_V(Bulk_atom_done,Cal_xanes,i_range,icheck_loc,korigimp,mpirank,n_atom_0,n_atom_ind,n_atom_proto,n_atom_proto_bulk, &
                 nlm_pot,npoint,nptmoy_out,nrm,nspin,rsbdc, &
                 rsbdc_out,V0bdcFimp,Vcato,Vh,Vhbdc,Vhbdc_init,Vhbdc_out,VmoyF,VmoyF_out,Vr)
  endif

  call Writing_atom_carac(Cal_xanes,Full_atom,iapot,icheck_loc,ipr1,iaprotoi,itypepr,mpirank,n_atom_0,n_atom_ind, &
        n_atom_proto,natome,nlm_pot,nrato,nrm,nspin,ntype,numat,Rato,Rchimp,Rchrg,rhoato, &
        Rmtg,Rmtg0,Rmtsd,V_intmax,Vcato,Vxcato,i_range)

  if( korigimp ) then
    V0bdcF(1:nspin) = V0bdcFimp(1:nspin)
    dV0bdcF(1:nspin) = V0bdcFimp(1:nspin) - VmoyF(1:nspin)
  else
    V0bdcF(1:nspin) = VmoyF(1:nspin)
    dV0bdcF(1:nspin) = 0._db
  endif
  if( Green ) then
    V0bdcF_out(:) = VmoyF(:)
    rsbdc_out = rsbdc
    VxcbdcF_out(:) = VxcbdcF(:)
    Vhbdc_out = Vhbdc
  else
    V0bdcF_out(1:nspin) = VmoyF_out(1:nspin)
  endif
! Seulement utilise dans la convolution pour calculer le libre parcours moyen
  V0muf = Workf + sum( V0bdcF(1:nspin) ) / nspin

  return
  110 format(/' ---- Potent_comp ',100('-'))
end

!***********************************************************************

! Routine making the selection of the FDM points for the calculation of the average potential in the material (muffin-tin zero)
! and at the outer sphere

subroutine Ptmoy(distai,Green,iaabs,iaproto,icheck,imoy,imoy_out,iopsymr,isrt,Moy_loc, &
                 n_atom_proto,natomp,nim,npoint,nsortf,npsom,nptmoy,nptmoy_out,nstm,poidsov,poidsov_out, &
                 pos,Rmtg0,Rsort,rvol,Sym_2D,xyz)

  use declarations
  implicit none

  integer:: i, i1, i2, ia, iaabs, ib, icheck, io, is, j, k, k1, k2, n_atom_proto, natomp, nim, npoint, ns, nsortf, npsom, nptmoy, &
            nptmoy_out, nstm

  integer, dimension(natomp):: iaproto
  integer, dimension(npoint):: imoy, imoy_out
  integer, dimension(nstm):: isrt
  integer, dimension(nopsm):: iopsymr

  logical:: Green, Moy_loc, Sym_2D

  real(kind=db):: D_max, dcour, dist, dist_min, distai, pdmoy, pdmoy_out, R_centre, Rsort, Rvmmax, z_max

  real(kind=db), dimension(3):: p, ps
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(natomp):: ray
  real(kind=db), dimension(npoint):: poidsov, poidsov_out
  real(kind=db), dimension(0:n_atom_proto):: Rmtg0
  real(kind=db), dimension(nim):: rvol
  real(kind=db), dimension(4,npsom):: xyz

  D_max = 0._db
  do ia = 1,natomp
    ray(ia) = Rmtg0( iaproto(ia) )
    D_max = max( Ray(ia), D_max )
  end do
  D_max = 2 * D_max

  do  ! boucle reduction

  nptmoy_out = 0

  if( .not. Green ) then

    pdmoy_out = 0._db
    j = 0

    boucle_ib: do ib = 1,nsortf
      i = isrt(ib)
      do ia = 1,natomp
        ps(1:3) = pos(1:3,ia) - xyz(1:3,i)
        dist = sqrt( sum( ps(:)**2 ) )
        if( dist < ray(ia) - eps6 ) cycle boucle_ib
      end do
      j = j + 1
      imoy_out(j) = i
      pdmoy_out = pdmoy_out + rvol(i)
    end do boucle_ib

    nptmoy_out = j
    do i = 1,nptmoy_out
      poidsov_out(i) = rvol(imoy_out(i))
    end do

  endif

  z_max = -1000000000._db
  if( Sym_2D ) then
    do ia = 1,natomp
      z_max = max( z_max, pos(3,ia) )
    end do
    z_max = z_max + eps6
  endif
 ! Moy_loc is usual false (see fdm.f)
  if( Moy_loc ) then
    dcour = Rsort
    do ia = 1,natomp
      if( ia == iaabs ) cycle
      ps(:) = pos(:,ia) - pos(:,iaabs)
      dcour = min(dcour,sqrt( sum( ps(:)**2 ) ))
      exit
    end do
    Rvmmax = min(dcour,Rsort) + eps6
  else
    Rvmmax = Rsort + eps6
! Distai is the distance from the center for the farest atom
    if( Distai > eps6 ) then
      R_centre = Distai + eps6
    else
      R_centre = Rsort + eps6
    endif
  endif

  ns = 1
  ib = 0
  pdmoy = 0._db
  do i = 1,npoint
    boucle_is: do is = 1,ns
      p(1:3) = xyz(1:3,i)
      if( is > 1  ) then
        select case(is)
          case(2,3,4)
            k = is - 1
            io = 39 + k
            if( iopsymr(io) == 1 .and. abs(p(k)) > epspos ) then
              p(k) = - p(k)
            else
              cycle
            endif
          case(5,6,7)
            k = is - 4
            k1 = 1 + mod(k,3)
            k2 = 1 + mod(k+1,3)
            i1 = 39 + k1
            i2 = 39 + k2
            if( iopsymr(i1) == 1 .and. abs(p(k1)) > epspos .and. iopsymr(i2) == 1 .and. abs(p(k2)) > epspos ) then
              p(k1) = - p(k1)
              p(k2) = - p(k2)
            else
              cycle
            endif
          case(8)
            if( iopsymr(40) == 1 .and. abs(p(1)) > epspos .and. iopsymr(41) == 1 .and. abs(p(2)) > epspos &
          .and. iopsymr(42) == 1 .and. abs(p(3)) > epspos ) then
              p(1:3) = - p(1:3)
            else
              cycle
            endif
        end select
      endif

      if( Sym_2D .and. p(3) > z_max ) cycle ! One eliminates the points beyond the surface layer
      if( Moy_loc ) then
        ps(:) = p(:) - pos(:,iaabs)
        dist = sqrt( sum( ps(:)**2 ) )
      else
        ps(:) = p(:)
        dist = sqrt( sum( ps(:)**2 ) )
        if( dist > R_centre ) cycle
      endif
      if( dist > rvmmax ) cycle
      dist_min = 10000000._db
      do ia = 1,natomp
        ps(:) = pos(:,ia) - p(:)
        dist = sqrt( sum( ps(:)**2 ) )
        if( dist < ray(ia) - eps6  ) cycle boucle_is
        dist_min = min( dist_min, dist )
      end do
      if( dist_min > D_max ) cycle
      ib = ib + 1
      imoy(ib) = i
      pdmoy = pdmoy + rvol(i)
    end do boucle_is
  end do
  nptmoy = ib

  do i = 1,nptmoy
    poidsov(i) = rvol(imoy(i))
  end do

  if( .not. ( nptmoy == 0 .or. ( nptmoy_out == 0 .and. .not. Green ) ) ) exit

! One starts again at the begining of the routine
    ray(:) = 0.9_db * ray(:)

  end do

  if( icheck > 2 ) then
    write(3,110) nptmoy, pdmoy
    write(3,120) (imoy(ib), xyz(1:4,imoy(ib))*bohr, poidsov(ib), ib = 1,nptmoy)
    if( .not. Green ) then
      write(3,130)
      write(3,110) nptmoy_out, pdmoy_out
      write(3,120) (imoy_out(ib), xyz(1:4,imoy_out(ib))*bohr, poidsov_out(ib), ib = 1,nptmoy_out)
    endif
  endif
  poidsov(1:nptmoy) = poidsov(1:nptmoy) / pdmoy
  if( .not. Green ) then
   poidsov_out(1:nptmoy_out) = poidsov_out(1:nptmoy_out) / pdmoy_out
  end if

  return
  110 format(/' nptmoy =',i4,'  pdmoy =',f10.5,/ &
              '  imoy       x         y        z         r      poidsov')
  120 format(i6,5f10.5)
  130 format(/' For the outer sphere :')
end

!***********************************************************************

! Calcul du potentiel moyen dans la zone interstitielle.

subroutine Cal_Vmoy(Cal_xanes,i_range,i_self,icheck,imoy,imoy_out,korigimp,Magnetic,mpirank,npoint,nptmoy,nptmoy_out,nspin, &
            poidsov,poidsov_out,rs,rsbdc,rsbdc_out,SCF,V0bdcFimp,Vh,Vhbdc,Vhbdc_out,VmoyF,VmoyF_out,Vr,VxcbdcF,VxcbdcF_out)

  use declarations
  implicit none

  integer:: i_self,icheck, ipr, mpirank, npoint, nptmoy, nptmoy_out, ispin, nspin, i_range

  integer, dimension(npoint):: imoy, imoy_out

  logical:: Cal_xanes, korigimp, Magnetic, SCF

  real(kind=db):: rsbdc, rsbdc_out, Vhbdc, Vhbdc_out

  real(kind=db), dimension(nspin):: VmoyF, VmoyF_out, V0bdcFimp, VxcbdcF, VxcbdcF_out
  real(kind=db), dimension(npoint):: poidsov, poidsov_out, rs, Vh
  real(kind=db), dimension(npoint,nspin):: Vr

! Calcul du potentiel moyen :
  if( nptmoy > 0 ) then
    do ispin = 1,nspin
      VmoyF(ispin) = sum( Vr(imoy(1:nptmoy),ispin) * poidsov(1:nptmoy) )
    end do

    if( Magnetic ) then
      VmoyF(1) = 0.5_db * ( vmoyF(1) + vmoyF(nspin) )
      VmoyF(nspin) = vmoyF(1)
    endif

    Vhbdc = sum( Vh(imoy(1:nptmoy)) * poidsov(1:nptmoy) )

    VxcbdcF(:) = VmoyF(:) - Vhbdc

! La moyenne est faite non pas sur le rayon rs, mais sur EF = ckf / rs**2, ckf constante.
    rsbdc = sum( (1/rs(imoy(1:nptmoy))**2) * poidsov(1:nptmoy) )
    rsbdc = 1 / sqrt( rsbdc )

    if( mpirank == 0 .and. ( i_self == 1 .or. .not. SCF ) ) then
      do ipr = 3,6,3
        if( ipr == 6 .and. i_range /= 1 ) cycle
        if( icheck == 0 .and. ipr == 3 ) cycle
        if( .not. Cal_xanes .and. ipr == 6 ) cycle
        if( korigimp ) then
          write(ipr,130) V0bdcFimp(1) * rydb, Vhbdc * rydb, rsbdc
        else
          write(ipr,130) VmoyF(1) * rydb, Vhbdc * rydb, rsbdc
        endif
      end do
    endif

  endif

  if( nptmoy_out > 0 ) then
    do ispin = 1,nspin
      VmoyF_out(ispin) = sum( Vr(imoy_out(1:nptmoy_out),ispin) * poidsov_out(1:nptmoy_out) )
    end do

    if( Magnetic ) then
      VmoyF_out(1) = 0.5_db * ( vmoyF_out(1) + vmoyF_out(nspin) )
      VmoyF_out(nspin) = vmoyF_out(1)
    endif

    Vhbdc_out = sum( Vh(imoy_out(1:nptmoy_out)) * poidsov_out(1:nptmoy_out) )

    VxcbdcF_out(:) = VmoyF_out(:) - Vhbdc_out

! La moyenne est faite non pas sur le rayon rs, mais sur EF = ckf / rs**2, ckf constante.
    rsbdc_out = sum( (1/rs(imoy_out(1:nptmoy_out))**2) * poidsov_out(1:nptmoy_out) )
    rsbdc_out = 1 / sqrt( rsbdc_out )

    if( mpirank == 0 .and. ( i_self == 1 .or. .not. SCF ) ) then
      do ipr = 3,6,3
        if( ipr == 6 .and. i_range /= 1 ) cycle
        if( icheck == 0 .and. ipr == 3 ) cycle
        if( .not. Cal_xanes .and. ipr == 6 ) cycle
        write(ipr,150) VmoyF_out(1) * rydb,  Vhbdc_out * rydb, rsbdc_out
      end do
    endif

  endif

  return
  130 format(/'     V0bdcF =',f12.5,' eV,     Vhbdc =',f10.3,' eV,     rsbdc =',f8.3,' a.u.')
  150 format(' V0bdcF_out =',f12.5,' eV, Vhbdc_out =',f10.3,' eV, rsbdc_out =',f8.3,' a.u.')
end

!***********************************************************************

  subroutine Shift_V(Bulk_atom_done,Cal_xanes,i_range,icheck,korigimp,mpirank,n_atom_0,n_atom_ind,n_atom_proto,n_atom_proto_bulk, &
                 nlm_pot,npoint,nptmoy_out,nrm,nspin,rsbdc, &
                 rsbdc_out,V0bdcFimp,Vcato,Vh,Vhbdc,Vhbdc_init,Vhbdc_out,VmoyF,VmoyF_out,Vr)

  use declarations
  implicit none

  integer:: i_range, icheck, ipr, mpirank, n_atom_0, n_atom_ind, n_atom_proto, n_atom_proto_bulk, nlm_pot, npoint, nptmoy_out, &
            nrm, nspin

  logical:: Bulk_atom_done, Cal_xanes, korigimp

  real(kind=db):: Delta_V, rsbdc, rsbdc_out, Vhbdc, Vhbdc_init, Vhbdc_out

  real(kind=db), dimension(nspin):: V0bdcFimp, VmoyF, VmoyF_out
  real(kind=db), dimension(npoint):: Vh
  real(kind=db), dimension(npoint,nspin):: Vr
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato

  Delta_V = Vhbdc_init - Vhbdc

  Vhbdc = Vhbdc_init
  Vhbdc_out = Vhbdc_out + Delta_V
  VmoyF(:) = VmoyF(:) + Delta_V
  VmoyF_out(:) = VmoyF_out(:) + Delta_V

  Vr(:,:) = Vr(:,:) + Delta_V
  Vh(:) = Vh(:) + Delta_V
  do ipr = n_atom_0,n_atom_ind
! bulk atom are not shifted, and we are necessarily in non Full_atom_mode in this case
    if( Bulk_atom_done .and. ipr > n_atom_proto - n_atom_proto_bulk ) exit
    Vcato(:,:,ipr) = Vcato(:,:,ipr) + Delta_V
  end do

  if( icheck > 0 ) write(3,110) Delta_V * rydb

  if( mpirank == 0 ) then
    do ipr = 3,6,3
      if( ipr == 6 .and. i_range /= 1 ) cycle
      if( icheck == 0 .and. ipr == 3 ) cycle
      if( .not. Cal_xanes .and. ipr == 6 ) cycle
      if( korigimp ) then
        write(ipr,130) V0bdcFimp(1) * rydb, Vhbdc * rydb, rsbdc
      else
        write(ipr,130) VmoyF(1) * rydb, Vhbdc * rydb, rsbdc
      endif
      if( nptmoy_out > 0 ) write(ipr,150) VmoyF_out(1) * rydb,  Vhbdc_out * rydb, rsbdc_out
    end do
  endif

  return
  110 format(/' To keep the interstitial potential constant, potentials are shifted by:',f10.5,' eV')
  130 format(/'     V0bdcF =',f10.3,' eV,     Vhbdc =',f10.3,' eV,     rsbdc =',f8.3,' a.u.')
  150 format(' V0bdcF_out =',f10.3,' eV, Vhbdc_out =',f10.3,' eV, rsbdc_out =',f8.3,' a.u.')
end

!***********************************************************************

subroutine Writing_atom_carac(Cal_xanes,Full_atom,iapot,icheck,ipr1,iaprotoi,itypepr,mpirank,n_atom_0,n_atom_ind, &
        n_atom_proto,natome,nlm_pot,nrato,nrm,nspin,ntype,numat,Rato,Rchimp,Rchrg,rhoato, &
        Rmtg,Rmtg0,Rmtsd,V_intmax,Vcato,Vxcato,i_range)

  use declarations
  implicit none

  integer:: i_range, iapr, icheck, ipr, ipr1, iprt, ir, ispin, it, mpirank, n, n_atom_0, n_atom_ind, n_atom_proto, &
            natome, nlm_pot, nr, nrm, nspin, ntype, Z

  integer, dimension(natome):: iaprotoi
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: iapot, itypepr, nrmtg

  character(len=1):: Star

  logical:: Cal_xanes, Full_atom

  real(kind=db):: f_integr3, p1, Rayion, V_intmax

  real(kind=db), dimension(:), allocatable:: r, rhr2
  real(kind=db), dimension(0:ntype):: Rchimp
  real(kind=db), dimension(0:n_atom_proto):: charge_ion, Radius_ion, Rchrg, Rmtg, Rmtg0, Rmtsd
  real(kind=db), dimension(0:n_atom_proto,nspin):: charge, VmftF
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind)::Vcato
  real(kind=db), dimension(0:nrm,nspin,n_atom_0:n_atom_ind):: rhoato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vxcato
  real(kind=db), dimension(0:nrm,0:ntype):: Rato

  do ipr = ipr1,n_atom_proto

    if( iapot(ipr) == 0 .or. Rmtg0(ipr) < eps10 ) cycle

    it = itypepr(ipr)
    Z = numat(it)
    nr = nrato(it)

    if( abs( Rchimp(it) ) > eps10 ) then
      Rchrg(ipr) = Rchimp(it)
    else
      Rchrg(ipr) = Rmtsd(ipr)
    endif

    if( Z > 0 ) then
      Radius_ion(ipr) = Rayion( Z ) / bohr
    else
      Radius_ion(ipr) = Rchrg(ipr)
    endif

    do ir = 2,nr-1
      if( Rato(ir,it) > Rmtg(ipr) - eps10 ) exit
    end do
    nrmtg(ipr) = ir

    charge_ion(ipr) = real( numat(it),db)
    charge(ipr,:) = real( numat(it),db) / nspin

    if( Full_atom ) then
      do iapr = n_atom_0,n_atom_ind
        if( ipr == iaprotoi(iapr) ) exit
      end do
      if( iapr > n_atom_ind ) cycle
    else
      iapr = ipr
    endif

    allocate( r(0:nr) )
    allocate( rhr2(0:nr) )
    r(0:nr) = Rato(0:nr,it)

    Charge_ion(ipr) = 0._db

    do ispin = 1,nspin
      rhr2(0:nr) = rhoato(0:nr,ispin,iapr) * r(0:nr)**2
      charge(ipr,ispin) = quatre_pi * f_integr3(r,rhr2,0,nr,Rchrg(ipr))
      charge(ipr,ispin) = real(Z, db) / nspin - charge(ipr,ispin)

      Charge_ion(ipr) = Charge_ion(ipr) + quatre_pi * f_integr3(r,rhr2,0,nr,Radius_ion(ipr))
    end do
    Charge_ion(ipr) = Z - Charge_ion(ipr)

! Calcul du potentiel au rayon muffintin :
    n = nrmtg(ipr)
    p1 = ( Rmtg(ipr) - r(n-1) ) / ( r(n) - r(n-1) )
    VmftF(ipr,1:nspin) = p1 * ( Vcato(n,1,iapr) + Vxcato(n,1,1:nspin,iapr) ) &
                       + ( 1 - p1 ) * ( Vcato(n-1,1,iapr) + Vxcato(n-1,1,1:nspin,iapr) )
    VmftF(ipr,:) = min( VmftF(ipr,:), V_intmax )

    deallocate( r )
    deallocate( rhr2 )

  end do

! The part below is suppressed. It was in order to have full charge in the atoms and to get neutrality in the unit cell
! It is confusing

 ! do ipr = 1,n_atom_proto
 !   if( iapot( ipr ) /= 0 .and. Rmtg0(ipr) > eps10 ) cycle
 !   All_found = .false.
 !   exit
 ! end do
 ! if( All_found ) then
 !   n_elec_tot = 0
 !   ch_ion_tot = 0._db
 !   charge_tot = 0._db
 !   do ipr = 1,n_atom_proto
 !     it = itypepr( ipr )
 !     n_elec_tot = n_elec_tot + ngreq(ipr) * numat( it )
 !     charge_tot = charge_tot + ngreq(ipr) * sum( charge(ipr,:) )
 !     ch_ion_tot = ch_ion_tot + ngreq(ipr) * charge_ion(ipr)
 !   end do
 !   rap_ch = n_elec_tot / charge_tot
 !   rap_io = n_elec_tot / ch_ion_tot
 !   do ipr = ipr1,n_atom_proto
 !     it = itypepr( ipr )
 !     if( abs( Rchimp(it) ) < eps10 ) charge(ipr,:) = rap_ch * charge(ipr,:)
 !     charge_ion(ipr) = rap_io * charge_ion(ipr)
 !   end do
 ! endif
 ! do ipr = ipr1,n_atom_proto
 !   if( iapot(ipr) == 0 .or. Rmtg0(ipr) < eps10 ) cycle
 !   charge_ion(ipr) = numat( itypepr(ipr) ) - charge_ion(ipr)
 ! end do

  if( mpirank == 0 ) then
    do iprt = 3,6,3
      if( iprt == 6 .and. i_range /= 1 ) cycle
      if( icheck == 0 .and. iprt == 3 ) cycle
      if( .not. Cal_xanes .and. iprt == 6 ) cycle
      if( nspin == 1 ) then
        write(iprt,120)
      else
        write(iprt,130)
      endif
      do ipr = ipr1,n_atom_proto
        if( iapot(ipr) == 0 .or. Rmtg0(ipr) < eps10 ) cycle
        if( Full_atom ) then
          do iapr = n_atom_0,n_atom_ind
            if( ipr == iaprotoi(iapr) ) exit
          end do
          if( iapr > n_atom_ind ) cycle
        endif
        it = itypepr(ipr)
        Z = numat( it )
        if( it > 0 ) then
          Star = ' '
        else
          Star = '*'
        endif
        write(iprt,140) ipr, Z, Star, Charge(ipr,1:nspin), Rchrg(ipr) * bohr, charge_ion(ipr),  Radius_ion(ipr) * bohr, &
                          VmftF(ipr,1:nspin) * rydb
      end do
    end do
  endif

  return
  120 format(/' ipr   Z     Atom-Charge  Atom-Radius Ionic-Charge Ionic-Radius     Vmft')
  130 format(/' ipr   Z     Atom-Ch(up)  Atom-Ch(dn)  Atom-Radius   Ionic-Ch   Ionic-Radius   Vmft(up)     Vmft(dn)')
  140 format(i3,i5,a1,7f13.3)

end

!***********************************************************************

subroutine Check_Ecmax(Ecineticmax,Ecineticmax_out,Eclie,Eclie_out,Eneg,Energ_max,mpirank,nspin,V0bdcF,V0bdcF_out,Workf)

  use declarations
  implicit none

  integer:: ipr, mpirank, nspin

  logical:: Eneg

  real(kind=db):: Ecineticmax, Ecineticmax_out, Eclie, Eclie_out, Energ_max, Workf
  real(kind=db), dimension(nspin):: V0bdcF, V0bdcF_out

  Ecineticmax = Energ_max - Workf - min( V0bdcF(1), V0bdcF(nspin) )
  Ecineticmax_out = Energ_max - Workf - min( V0bdcF_out(1), V0bdcF_out(nspin) )

! Pour tenir compte de l'eventuel baisse de Vbd apres le niveau de Fermi, on ajoute 1 eV :
  Ecineticmax = Ecineticmax + 1._db / rydb

  if( .not. Eneg ) then
    Ecineticmax = max( Ecineticmax, Eclie )
    Ecineticmax_out = max( Ecineticmax_out, Eclie_out )
  endif
  if( ( Ecineticmax < eps10 .or. Ecineticmax_out < eps10 ) .and. .not. Eneg .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,110) Ecineticmax * rydb, Ecineticmax_out * rydb
    end do
    close(9)
    stop
  endif

  return
  110 format(/' E_kinetic_max =',f7.3,' eV  or E_kinetic_max_out =',f7.3,' eV < 0.',/ &
              ' Start the calculation at higher energy !'///)
end

!***********************************************************************

! Calcul du potentiel dans l'etat excite.
! Amene une modification de Vr.

subroutine potex(Nonsph,axyz,alfpot,dv0bdcF,Energ,Enervide,Full_atom, &
          iaabsi,iapot,iaprotoi,icheck,iprabs,iprabs_reel,itab,itypepr,Magnetic, &
          n_atom_0,n_atom_ind,n_atom_proto,n_vr_0,n_vr_ind,natome,nlm_pot,Nonexc_g,npoint,npoint_ns,npsom,nptmoy, &
          nptmoy_out,nrato,nrm,nspin,ntype,Optic,Rato,rho,rhons,Rmtg,rs,rsato,rsbdc,rsbdc_out,Trace_k,Trace_p, &
          V_intmax,Vcato,Vh,Vhbdc,Vhbdc_out,Vhns,Vr,Vxc,VxcbdcF,VxcbdcF_out,Vrato,Vxcato,V0bdc,V0bdc_out,xyz)

  use declarations
  implicit none

  integer:: i, ia, iaabsi, iapr, icheck, ipr, iprabs, iprabs_reel, ir, ispin, it, itab, n_atom_0, n_atom_ind, &
    n_atom_proto, n_vr_0, n_vr_ind, natome, nlm_pot, npoint, npoint_ns, npsom, nptmoy, nptmoy_out, nr, nrm, nspin, ntype, Trace_k

  integer, dimension(0:ntype):: nrato
  integer, dimension(natome):: iaprotoi
  integer, dimension(0:n_atom_proto):: iapot, itypepr

  logical:: Optic, Full_atom, magnetic, nonexc_g, Nonsph

  real(kind=db):: alfpot, Energ, Enervide, p, p1, rsbdc, rsbdc_out, V_intmax, Vhbdc, Vhbdc_out
  real(kind=db), dimension(3):: axyz
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(nspin):: dV0bdcF, V0bdc, V0bdc_out, Vmftabs, VxcbdcF, VxcbdcF_out
  real(kind=db), dimension(npoint):: rs, Vh
  real(kind=db), dimension(npoint_ns):: rhons, Vhns
  real(kind=db), dimension(npoint,nspin):: rho, Vxc, Vr, Vr_td
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(0:n_atom_proto):: Rmtg
  real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: rsato
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_vr_0:n_vr_ind):: Vrato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato
  real(kind=db), dimension(0:nrm,0:ntype):: Rato
  real(kind=db), dimension(:), allocatable:: rst, Vct, Vrt, Vxct

  if( icheck > 1 ) then
    write(3,100)
    write(3,110) Energ * rydb
  endif

  do ispin = 1,nspin

! Calcul des potentiels atomiques dans l'etat excitee.

    do ia = n_atom_0,n_atom_ind
      if( Full_atom ) then
        ipr = iaprotoi(ia)
      else
        ipr = ia
        if( nonexc_g .and. ipr == 0 ) cycle
      endif
      if( iapot(ipr) == 0 ) cycle

      nr = nrato( itypepr(ipr) )

      if( alfpot < eps6 .and. .not. Optic ) then

        allocate( rst(nr) )
        allocate( Vct(nr) )
        allocate( Vrt(nr) )
        allocate( Vxct(nr) )

        Vct(1:nr) = Vcato(1:nr,1,ia)
        Vxct(1:nr) = Vxcato(1:nr,1,ispin,ia)
        rst(1:nr) = rsato(1:nr,ia)

        call subpotex(nr,Vrt,Vct,Vxct,rst,Enervide)

        Vrato(1:nr,1,ispin,ia) = Vrt(1:nr)

        do ir = 1,nr
          p = ( Vrt(ir) - Vcato(ir,1,ia) ) / Vxct(ir)
          Vrato(ir,2:nlm_pot,ispin,ia) = Vcato(ir,2:nlm_pot,ia) + p * Vxcato(ir,2:nlm_pot,ispin,ia)
        end do

        deallocate( rst )
        deallocate( Vct )
        deallocate( Vrt )
        deallocate( Vxct )

      else

        Vrato(1:nr,1:nlm_pot,ispin,ia) = Vcato(1:nr,1:nlm_pot,ia) + Vxcato(1:nr,1:nlm_pot,ispin,ia)

      endif

    end do

    if( npoint > 0 ) then
      allocate( rst(npoint) )
      allocate( Vct(npoint) )
      allocate( Vrt(npoint) )
      allocate( Vxct(npoint) )

      Vct(1:npoint) = Vh(1:npoint)
      Vxct(1:npoint) = Vxc(1:npoint,ispin)
      rst(1:npoint) = rs(1:npoint)

      if( alfpot < eps6 .and. .not. Optic ) then
        call subpotex(npoint,Vrt,Vct,Vxct,rst,Enervide)
        Vr_td(1:npoint,ispin) = Vrt(1:npoint)
      else
        Vr_td(1:npoint,ispin) = Vct(1:npoint) + Vxct(1:npoint)
      endif

      deallocate( rst )
      deallocate( Vct )
      deallocate( Vrt )
      deallocate( Vxct )

      if( V_intmax < 1000._db ) then
        do i = 1,npoint
          Vr_td(1:npoint,ispin) = min( Vr_td(1:npoint,ispin), V_intmax )
        end do
      endif

    endif

! Pour le potentiel interstitiel, on fait l'approximation que la moyenne du potentiel excite est egale
! au potentiel fonction de la moyenne de rs (rsbdc) et VxcF (VxcbdcF).
! En fait on a pris rsbc = 1/sqrt( <1/rs**2> ), c'est donc une moyenne sur EF.
    if( nptmoy > 0 ) then
      allocate( rst(1) )
      allocate( Vct(1) )
      allocate( Vrt(1) )
      allocate( Vxct(1) )

      Vct(1) = Vhbdc
      Vxct(1) = VxcbdcF(ispin)
      rst(1) = rsbdc

      if( alfpot < eps6 .and. .not. Optic ) then
        call subpotex(1,Vrt,Vct,Vxct,rst,Enervide)
        V0bdc(ispin) = Vrt(1)
      else
        V0bdc(ispin) = Vct(1) + Vxct(1)
      endif

      deallocate( rst )
      deallocate( Vct )
      deallocate( Vrt )
      deallocate( Vxct )

      if( V_intmax < 1000._db ) V0bdc(ispin) = min( V0bdc(ispin), V_intmax )
    endif

    if( nptmoy_out > 0 ) then
      allocate( rst(1) )
      allocate( Vct(1) )
      allocate( Vrt(1) )
      allocate( Vxct(1) )

      Vct(1) = Vhbdc_out
      Vxct(1) = VxcbdcF_out(ispin)
      rst(1) = rsbdc_out

      if( alfpot < eps6 .and. .not. Optic ) then
        call subpotex(1,Vrt,Vct,Vxct,rst,Enervide)
        V0bdc_out(ispin) = Vrt(1)
      else
        V0bdc_out(ispin) = Vct(1) + Vxct(1)
      endif

      deallocate( rst )
      deallocate( Vct )
      deallocate( Vrt )
      deallocate( Vxct )

      if( V_intmax < 1000._db ) V0bdc_out(ispin) = min( V0bdc_out(ispin), V_intmax )

    else
      V0bdc_out(ispin) = V0bdc(ispin)
    endif

    it = itab
    if( iprabs_reel == 0 .and. nonexc_g ) then
      ipr = iprabs
    else
      ipr = iprabs_reel
    endif
    do ir = 1,nrato(it)
      if( Rato(ir,it) > Rmtg(ipr) ) exit
    end do
    p1 = ( Rmtg(ipr) - Rato(ir-1,it) ) / ( Rato(ir,it) - Rato(ir-1,it) )
    if( Full_atom ) then
      iapr = iaabsi
    else
      iapr = ipr
    endif
    Vmftabs(ispin) = p1 * Vrato(ir,1,ispin,iapr) + ( 1 - p1 ) * Vrato(ir-1,1,ispin,iapr)
    if( V_intmax < 1000._db ) Vmftabs(ispin) = min( Vmftabs(ispin), V_intmax )

  end do

  V0bdc(1:nspin) = V0bdc(1:nspin) + dV0bdcF(1:nspin)

  if( Magnetic ) then
    V0bdc(1) = 0.5_db * ( V0bdc(1) + V0bdc(nspin) )
    V0bdc(nspin) = V0bdc(1)
    V0bdc_out(1) = 0.5_db * ( V0bdc_out(1) + V0bdc_out(nspin) )
    V0bdc_out(nspin) = V0bdc_out(1)
  endif

  Vr(:,:) = Vr_td(:,:)

  if( icheck > 1 ) then
    if( nspin == 1 ) then
      write(3,120) V0bdc(1) * rydb, Vmftabs(1:nspin) * rydb
    else
      write(3,125) V0bdc(1) * rydb, Vmftabs(1:nspin) * rydb
    endif
    if( nptmoy_out > 0 ) write(3,127) V0bdc_out(1) * rydb

    if( npoint > 0 ) then
      do ispin = 1,nspin
        write(3,130) ispin
        write(3,140) (i, Vr(i,ispin)*rydb, i = 1,npoint)
      end do
    endif
    do iapr = n_vr_0,n_vr_ind
      if( Full_atom ) then
        ipr = iaprotoi( iapr )
      else
        ipr = iapr
        if( nonexc_g .and. ipr == 0 ) cycle
        if( iapot(ipr) == 0 ) cycle
      endif
      it = itypepr(ipr)
      nr = nrato(it)
      if(nspin == 1 ) then
        write(3,150) iapr
      else
        write(3,155) iapr
      endif
      do ir = 1,nr
        write(3,160) Rato(ir,it)*bohr, Vrato(ir,1,:,iapr)*rydb, Vcato(ir,1,iapr)*rydb, Vxcato(ir,1,:,iapr)*rydb
      end do
    end do
  endif

  if( Trace_k /= 0 .and. icheck > 0 ) then
    do ispin = 1,nspin
      call Trace(Nonsph,axyz,npoint,npoint_ns,ispin,npsom,xyz,rhons,Trace_k,Trace_p,Vhns,nspin,Vr,rho)
    end do
  endif

  return
  100 format(/' ---- Potex --------',100('-'))
  110 format(/' Energ =',f10.5,' eV')
  120 format(/' V0bdc =',f10.5,' eV, Vmftabs =',f10.5,' eV',f10.5)
  125 format(/' V0bdc =',f10.5,' eV, Vmftabs(up) =',f10.5,' eV, Vmftabs(dn) =',f10.5,' eV')
  127 format(/' V0bdc_out =',f10.5,' eV')
  130 format(/4x,'i     Vr_(eV)     ispin = ',i2)
  140 format(5(i5,e15.5))
  150 format('  ipr =',i3,/7x,'rato_(A)    Vrato_(eV)     Vcato(eV)     Vxcato(eV)')
  155 format('  ipr =',i3,/7x,'rato_(A)   Vrato(up)_(eV) Vrato(dn)_(eV)     Vcato(eV) Vxcato(up)_(eV) Vxcato(dn)_(eV)')
  160 format(f15.6,1p,5e15.5)
end

!***********************************************************************

! Calculation of atomic excited potential (when extracted)

subroutine Potex_abs(alfpot,dV0bdcF,Enervide,icheck,nlm_pot,nrato_abs,nspin,Optic,rato_abs,rsato_abs,rsbdc,V_intmax, &
                     V0bdc,Vcato,Vhbdc,Vrato,Vxcato,VxcbdcF)

  use declarations
  implicit none

  integer:: icheck,ir, isp, lm, nrato_abs, nlm_pot, nspin

  logical:: Optic

  real(kind=db):: alfpot, rsbdc, Enervide, V_intmax, Vhbdc
  real(kind=db), dimension(1):: rs, Vc, Vr, Vxc
  real(kind=db), dimension(nspin):: dV0bdcF,  V0bdc, VxcbdcF
  real(kind=db), dimension(nrato_abs):: rato_abs, rsato_abs, Vct, Vrt, Vxct
  real(kind=db), dimension(nrato_abs,nlm_pot):: Vcato
  real(kind=db), dimension(nrato_abs,nlm_pot,nspin):: Vrato, Vxcato

  if( alfpot < eps6 .and. .not. Optic )  then

    Vct(:) = Vcato(:,1)
    do isp = 1,nspin
      Vxct(:) = Vxcato(:,1,isp)
      call subpotex(nrato_abs,Vrt,Vct,Vxct,rsato_abs,Enervide)
      Vrato(:,1,isp) = Vrt(:)
      do lm = 2,nlm_pot
        Vrato(:,lm,isp) = Vcato(:,lm) + ( ( Vrt(:) - Vcato(:,1) ) / Vxct(:) ) * Vxcato(:,lm,isp)
      end do
    end do

    rs(1) = rsbdc
    Vc(1) = Vhbdc
    do isp = 1,nspin
      Vxc(1) = VxcbdcF(isp)
      call subpotex(1,Vr,Vc,Vxc,rs,Enervide)
      V0bdc(isp) = Vr(1) + dV0bdcF(isp)
    end do

  else

    do isp = 1,nspin
      do lm = 1,nlm_pot
        Vrato(:,lm,isp) = Vcato(:,lm) + Vxcato(:,lm,isp)
      end do
      V0bdc(isp) = Vhbdc + VxcbdcF(isp) + dV0bdcF(isp)
    end do

  endif

  if( V_intmax < 1000._db ) V0bdc(:) = min( V0bdc(:), V_intmax )

  if( icheck > 2 ) then
    if( nspin == 1 ) then
      write(3,'(A)') '        Radius      Vrato_abs      Vcato_abs     Vxcato_abs'
    else
      write(3,'(A)') '        Radius      Vrato_abs_up   Vrato_abs_dn   Vcato_abs     Vxcato_abs_up   Vxcato_abs_dn'
    endif
    do ir = 1,nrato_abs
      write(3,120) rato_abs(ir)*bohr, ( Vrato(ir,1,isp)*rydb, isp = 1,nspin ), Vcato(ir,1)*rydb, &
                   ( Vxcato(ir,1,isp)*rydb, isp = 1,nspin )
    end do
  endif

  return
  120 format(f15.6,1p,5e15.5)

end

!***********************************************************************

! Calcul du potentiel dans l'etat excite selon Hedin et Lundqvist.
! Amene une modification de Vr.

subroutine subpotex(np,Vrt,Vct,Vxct,rst,Enervide)

  use declarations
  implicit none

  integer, parameter:: npkfm = 16
  integer, parameter:: nraym = 17

  real(kind=db), parameter:: fac = ( 4._db/(9._db*pi) )**(1._db/3._db)

  integer:: i, j, jj, np

  real(kind=db):: Enervide, p1, p2, pskf, Vhli, Vhlm, Vhlp

  real(kind=db), dimension(nraym,npkfm):: vhl
  real(kind=db), dimension(nraym,npkfm/2):: vhl1
  real(kind=db), dimension(nraym,npkfm/2+1:npkfm):: vhl2
  real(kind=db), dimension(npkfm):: pkf
  real(kind=db), dimension(nraym):: ray
  real(kind=db), dimension(np):: Vct, Vxct, Vrt, rst

! Tableau Von Barth (Hedin et Lundqvist) normalise a la valeur du niveau de Fermi.
  data pkf/ 0.00_db, 0.20_db, 0.40_db, 0.60_db, 0.80_db, 1.00_db, 1.20_db, 1.40_db, 1.60_db, 1.80_db, &
            2.00_db, 2.20_db, 2.40_db, 2.60_db, 2.80_db, 3.00_db/

  data ray/ 0.00_db, 0.01_db, 0.02_db, 0.05_db, 0.10_db, 0.20_db, 0.30_db, 0.40_db, 0.50_db, 0.70_db, &
            1.00_db, 1.50_db, 2.00_db, 3.00_db, 4.00_db, 5.00_db, 6.00_db/

  data vhl1/ 2.0000_db,1.8564_db,1.7983_db,1.7124_db,1.5881_db,1.4398_db,1.3486_db,1.2939_db,1.2501_db,1.1624_db, &
             1.0528_db,1.0134_db,0.9716_db,0.9541_db,0.9227_db,0.9525_db,0.9545_db, &
             1.9732_db,1.8353_db,1.7687_db,1.6943_db,1.5704_db,1.4241_db,1.3333_db,1.2814_db,1.2379_db,1.1508_db, &
             1.0555_db,1.0066_db,0.9767_db,0.9595_db,0.9552_db,0.9569_db,0.9581_db, &
             1.8897_db,1.7545_db,1.7012_db,1.6225_db,1.5058_db,1.3668_db,1.2836_db,1.2371_db,1.1993_db,1.1237_db, &
             1.0528_db,1.0100_db,0.9839_db,0.9689_db,0.9666_db,0.9654_db,0.9687_db, &
             1.7395_db,1.6139_db,1.5625_db,1.4909_db,1.3851_db,1.2679_db,1.2020_db,1.1697_db,1.1427_db,1.0888_db, &
             1.0422_db,1.0134_db,0.9939_db,0.9814_db,0.9777_db,0.9785_db,0.9790_db, &
             1.4944_db,1.3774_db,1.3341_db,1.2754_db,1.2057_db,1.1376_db,1.1026_db,1.0884_db,1.0757_db,1.0504_db, &
             1.0262_db,1.0118_db,1.0000_db,0.9922_db,0.9889_db,0.9915_db,0.9896_db, &
             1.0000_db,1.0000_db,1.0000_db,1.0000_db,1.0000_db,1.0000_db,1.0000_db,1.0000_db,1.0000_db,1.0000_db, &
             1.0000_db,1.0000_db,1.0000_db,1.0000_db,1.0000_db,1.0000_db,1.0000_db, &
             0.5604_db,0.6617_db,0.7030_db,0.7904_db,0.8324_db,0.8697_db,0.8940_db,0.9207_db,0.9304_db,0.9499_db, &
             0.9762_db,0.9913_db,0.9987_db,1.0069_db,1.0055_db,1.0132_db,1.0106_db, &
             0.3856_db,0.4730_db,0.5088_db,0.5868_db,0.6645_db,0.7708_db,0.8497_db,0.8982_db,0.9142_db,0.9460_db, &
             0.9735_db,0.9948_db,1.0084_db,1.0208_db,1.0168_db,1.0261_db,1.0279_db/

  data vhl2/ 0.2852_db,0.3624_db,0.3926_db,0.4731_db,0.5294_db,0.6219_db,0.6911_db,0.7629_db,0.8086_db,0.8998_db, &
             1.0262_db,1.0236_db,1.0294_db,1.0272_db,1.0446_db,1.0478_db,1.0419_db, &
             0.2204_db,0.2876_db,0.3165_db,0.3772_db,0.4352_db,0.5172_db,0.5813_db,0.6407_db,0.6807_db,0.7607_db, &
             0.8945_db,1.0168_db,1.0516_db,1.0699_db,1.0669_db,1.0651_db,1.0595_db, &
             0.1761_db,0.2365_db,0.2628_db,0.3174_db,0.3677_db,0.4412_db,0.4986_db,0.5514_db,0.5877_db,0.6602_db, &
             0.7730_db,0.9239_db,1.0110_db,1.1157_db,1.0949_db,1.0953_db,1.0804_db, &
             0.1441_db,0.2007_db,0.2224_db,0.2694_db,0.3175_db,0.3839_db,0.4356_db,0.4835_db,0.5167_db,0.5830_db, &
             0.6675_db,0.7988_db,0.8967_db,1.0356_db,1.1340_db,1.1689_db,1.1190_db, &
             0.1200_db,0.1707_db,0.1924_db,0.2334_db,0.2794_db,0.3394_db,0.3859_db,0.4287_db,0.4582_db,0.5174_db, &
             0.6015_db,0.7230_db,0.8133_db,0.9524_db,1.0557_db,1.1343_db,1.1749_db, &
             0.1016_db,0.1466_db,0.1684_db,0.2033_db,0.2496_db,0.3031_db,0.3448_db,0.3840_db,0.4103_db,0.4629_db, &
             0.5460_db,0.6608_db,0.7438_db,0.8899_db,0.9824_db,1.0713_db,1.1073_db, &
             0.0874_db,0.1286_db,0.1492_db,0.1795_db,0.2231_db,0.2731_db,0.3121_db,0.3474_db,0.3718_db,0.4206_db, &
             0.5011_db,0.6068_db,0.6860_db,0.8269_db,0.9099_db,0.9935_db,1.0376_db, &
             0.0757_db,0.1139_db,0.1343_db,0.1614_db,0.2024_db,0.2487_db,0.2851_db,0.3172_db,0.3388_db,0.3818_db, &
             0.4643_db,0.5627_db,0.6356_db,0.7732_db,0.8373_db,0.9199_db,0.9784_db/

  do i = 1,npkfm/2
    vhl(:,i) = vhl1(:,i)
  end do
  do i = npkfm/2+1,npkfm
    vhl(:,i) = vhl2(:,i)
  end do

  Vrt(1:np) = Vct(1:np) + Vxct(1:np)

  do i = 1,np
! Si on est au dela du tableau pour le rayon rs on prend la valeur du bord du tableau.
    if( ray(nraym) <= rst(i) ) then
      j = nraym
      rst(i) = ray(nraym)
    else
      do j = 2,nraym
        if( ray(j) > rst(i) ) exit
      end do
    endif

    if( Enervide  <=  Vrt(i) ) then
      pskf = 0._db
    else
  !    Ef = ( ( 9*pi/4 )**(1/3.)  / rs )**2
  !    pskf = sqrt( (Enervide - Vr) / Ef ) = rs * (4/(9*pi))**tiers * sqrt( (Enervide - Vr) / Ef )
      pskf = fac * rst(i) * sqrt( Enervide - Vrt(i) )
    endif

    p1 = ( ray(j) - rst(i) ) / ( ray(j) - ray(j-1) )
    p2 = 1 - p1
! Au dessus de la valeur max de pskf du tableau on prend le potentiel
! Vxc en 1/k.
    if( pskf > pkf(npkfm) ) then
      Vhli = p1 * Vhl(j-1,npkfm) + p2 * Vhl(j,npkfm)
      Vhli = Vhli * pkf(npkfm) / pskf
! En dessous de cette valeur on interpole a 2 dimensions dans le tableau
    else
      do jj = 1,npkfm
        if( pkf(jj) > pskf ) exit
      end do
      Vhlm = p1 * Vhl(j-1,jj-1) + p2 * Vhl(j,jj-1)
      Vhlp = p1 * Vhl(j-1,jj) + p2 * Vhl(j,jj)
      Vhli = ( (pkf(jj) - pskf) * Vhlm + (pskf - pkf(jj-1)) * Vhlp ) / ( pkf(jj) - pkf(jj-1) )
    endif
    Vrt(i) = Vct(i) + Vhli * Vxct(i)
  end do

  return
end

!***********************************************************************

! Trace en sortie des coupes du potentiel.

subroutine Trace(Nonsph,axyz,npoint,npoint_ns,ispin,npsom,xyz,rhons,Trace_k,Trace_p,vhns,nspin,Vr,rho)

  use declarations
  implicit none

  integer:: i, ispin, j, j_min, j1, j2, jdir, k, npoint, npoint_ns, npsom, nspin, ntrace, Trace_k

  integer, dimension(npoint):: itrace
  integer, dimension(:), allocatable:: jtrace

  logical:: Nonsph
  logical, dimension(:), allocatable:: ok

  real(kind=db):: ctrace, fac, p, x1, x2, z_min
  real(kind=db), dimension(3):: axyz, Ptrace, v, Vectrace
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(npoint_ns):: rhons, Vhns
  real(kind=db), dimension(npoint,nspin):: rho, Vr
  real(kind=db), dimension(4,npsom):: xyz

  Vectrace(1:3) = Trace_p(1:3)
  Ptrace(1:3) = Trace_p(4:6)

! Calcul des points ou on ecrit la fonction d'onde

  k = 0
  if( Trace_k == 1 ) then
    do j = 1,3
      if( abs(Vectrace(j)) > eps10 ) exit
    end do
    jdir = j
    if( jdir == 1 ) then
      j1 = 2
      j2 = 3
    elseif( jdir == 2 ) then
      j1 = 1
      j2 = 3
    else
      j1 = 1
      j2 = 2
    endif
    v(1:3) = Ptrace(1:3) * axyz(1:3)
    do  i = 1,npoint
      p = ( xyz(jdir,i) - v(jdir) ) / Vectrace(jdir)
      x1 = v(j1) + p * Vectrace(j1)
      x2 = v(j2) + p * Vectrace(j2)
      if( abs(xyz(j1,i)-x1) < eps6 .and. abs(xyz(j2,i)-x2) < eps6 ) then
        k = k+1
        itrace(k) = i
      endif
    end do
  elseif( Trace_k == 2 ) then
    ctrace = Ptrace(1) / bohr
    do i = 1,npoint
      fac = sum( Vectrace(1:3) * xyz(1:3,i) )
      if( abs( fac - ctrace ) < eps6 ) then
        k = k+1
        itrace(k) = i
      endif
    end do
  elseif( Trace_k == 3 ) then
    do i = 1,npoint
      itrace(i) = i
    end do
    k = npoint
  endif
  ntrace = k

  if( Trace_k == 1 ) then
    write(3,110) vectrace(1:3), ptrace(1:3)
  elseif( Trace_k == 2 ) then
    write(3,120) vectrace(1:3), ptrace(1)
  endif
  if( Nonsph ) then
    write(3,130) ntrace
  else
    write(3,135) ntrace
  endif
  if( Trace_k == 1 ) then
    do jdir = 1,3
      if( abs( abs( Vectrace(jdir) ) - 1 ) < eps10 ) exit
    end do
    if( jdir < 4 ) then
      allocate( ok(ntrace) )
      allocate( jtrace(ntrace) )
      ok(:) = .false.
      do i = 1,ntrace
        z_min = 1000000._db
        do j = 1,ntrace
          if( ok(j) ) cycle
          if( xyz(jdir,itrace(j)) > z_min ) cycle
          j_min = j
          z_min = xyz(jdir,itrace(j))
        end do
        ok(j_min) = .true.
        jtrace(i) = itrace(j_min)
      end do
      itrace(1:ntrace) = jtrace(1:ntrace)
      deallocate( jtrace, ok )
    endif
  endif
  do j = 1,ntrace
    k = itrace(j)
    if( Nonsph ) then
      write(3,140) xyz(1:3,k)*bohr, Vr(k,ispin)*rydb, vhns(k)*rydb, rho(k,ispin), rhons(k)
    else
      write(3,140) xyz(1:3,k)*bohr, Vr(k,ispin)*rydb, rho(k,ispin)
    endif
  end do

  return
  110 format(/' Cut along the line of unitary vector =',3f7.3,/ '  and passing through  =',3f7.3)
  120 format(/'Cut along the plane : ax + by + cz = d,',/ ' with : a =',f5.2,'  b = ',f5.2,'  c = ',f5.2,'  d =',f5.2)
  130 format(/i5,'  / ntrace',/3x, 'xval    yval    zval      Vr        Vhns       rho       rhons')
  135 format(/i5,'  / ntrace',/3x, 'xval    yval    zval      Vr       rho')
  140 format(3f8.4,1p,6e11.3)
end

!***********************************************************************

! On rend le potentiel muffintin.

subroutine mdfmuf(Nonsph,Axe_Atom_gr,axyz,Full_atom,iaabs,iaproto,iaprotoi,icheck, &
              igreq,igroup,ispin,itypep,n_atom_0,n_atom_ind,n_atom_proto,natome,natomp,neqm,ngroup_m,nlm_pot, &
              npoint,npoint_ns,npsom,nrato,nrm,nspin,ntype,pos,Rato,rho,rhons,Rmtg,Trace_k,Trace_p,Vhns,Vm,Vr,Vrato,xyz)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer:: Trace_k
  integer, dimension(natome):: iaprotoi
  integer, dimension(natomp):: iaproto, igroup, itypep
  integer, dimension(0:ntype):: nrato
  integer, dimension(0:n_atom_proto,neqm):: igreq

  logical:: Full_atom, Nonsph

  real(kind=db), dimension(3):: axyz, v
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_gr
  real(kind=db), dimension(npoint_ns):: rhons, vhns
  real(kind=db), dimension(npoint,nspin):: rho, Vr
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(0:n_atom_proto) :: Rmtg
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vrato
  real(kind=db), dimension(0:nrm,0:ntype):: Rato

  boucle_i: do i = 1,npoint
    do ia = 1,natomp

      if( nspin == 2 ) then
        igr = igroup(ia)
        iga = igroup(iaabs)
        kgr = igreq(ipr,1)

        cosang = sum( Axe_Atom_gr(:,kgr) * Axe_Atom_gr(:,igr) )
        cosang = cosang * abs( sum( Axe_Atom_gr(:,iga) * Axe_Atom_gr(:,igr) ) )
        if( abs(cosang - 1) < eps4 ) then
          itm = 1
          isp = ispin
        elseif( abs(cosang + 1) < eps4 ) then
          itm = - 1
          isp = 3 - ispin
        else
          itm = 0
        endif
      else
        itm = 1
        isp = ispin
      endif

      it = itypep(ia)
      ipr = iaproto(ia)
      v(1:3) = xyz(1:3,i) - pos(1:3,ia)
      r = sqrt( sum( v(:)**2 ) )
      if( r < Rmtg(ipr) ) then
        if( Full_atom ) then
          do iapr = 1,natome
            if( iaprotoi(iapr) == ipr ) exit
          end do
          if( iapr > natome ) cycle boucle_i
        else
          iapr = ipr
        endif
        do ir = 1,nrato(it)
          if( Rato(ir,it) > r ) exit
        end do
        p1 = ( r - Rato(ir-1,it) ) / ( Rato(ir,it) - Rato(ir-1,it) )
        if( itm == 0 ) then
          Vr(i,ispin) = 0.5_db * ( ( 1 - p1 ) * sum( Vrato(ir-1,1,:,iapr) ) + p1 * sum( Vrato(ir,1,:,iapr) ) )
        else
          Vr(i,ispin) = ( 1 - p1 ) * Vrato(ir-1,1,isp,iapr) + p1 * Vrato(ir,1,isp,iapr)
        endif
        cycle boucle_i
      endif
    end do

    Vr(i,ispin) = vm

  end do boucle_i

  if( icheck > 1 ) then
    write(3,130)
    write(3,140) (i, Vr(i,ispin)*rydb, i = 1,npoint)
  endif

  if( Trace_k /= 0 .and. icheck > 0 ) then
    do ispin = 1,nspin
      call Trace(Nonsph,axyz,npoint,npoint_ns,ispin,npsom,xyz,rhons,Trace_k,Trace_p,vhns,nspin,Vr,rho)
    end do
  endif

  return
  130 format(/4x,'i     Vr (eV)  apres muffin-tin')
  140 format(5(i5,f10.3))
end

!***********************************************************************

! Potential modification to make it muffin-tin, even in FDM (for testing).

subroutine modmuf(Full_atom,iaprotoi,itypepr,icheck, ispin,n_atom_0,n_atom_ind,n_atom_proto,natome,nlm_pot, &
            nrato,nrm,nspin,ntype,Rato,Rmtg,vm,Vrato)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(natome):: iaprotoi
  integer, dimension(0:ntype):: nrato
  integer, dimension(0:n_atom_proto):: itypepr

  logical:: Full_atom

  real(kind=db), dimension(0:n_atom_proto):: Rmtg
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vrato
  real(kind=db), dimension(0:nrm,0:ntype):: Rato

  do iapr = n_atom_0,n_atom_ind
    if( Full_atom ) then
      ipr = iaprotoi(iapr)
    else
      ipr = iapr
    endif
    it = itypepr(ipr)
    rm = 0.9_db * Rmtg(ipr)
    do ir = 1,nrato(it)
      if( Rato(ir,it) < rm ) cycle
      if( Rato(ir,it) > Rmtg(ipr) ) then
        Vrato(ir,1,ispin,iapr) = vm
        Vrato(ir,2:nlm_pot,ispin,iapr) = 0._db
      else
        p2 = ( Rato(ir,it) - rm ) / ( Rmtg(ipr) - rm )
        p1 = 1 - p2
        Vrato(ir,1,ispin,iapr) = p1 * Vrato(ir,1,ispin,iapr) + p2*vm
      endif
    end do
  end do

  if( icheck > 1 ) then
    write(3,110)
    do iapr = n_atom_0,n_atom_ind
      if( Full_atom ) then
        ipr = iaprotoi(iapr)
      else
        ipr = iapr
      endif
      it = itypepr(ipr)
      write(3,120) iapr, ispin
      do ir = 1,nrato(it)
        write(3,130) Rato(ir,it)*bohr, Vrato(ir,1,ispin,iapr)*rydb
      end do
    end do
  endif

  return
  110 format(/' ---- Modmuf -------',100('-'))
  120 format('     Rato (A)    Vrato (eV)    iapr =',i3,', ispin =',i2)
  130 format(f15.6,1p,e11.3)
end

!***********************************************************************

subroutine gradpot(Base_hexa,cgrad,gradvr,icheck,iord,ispin,ivois,nicm,npoint,npsom,nspin,nvois,Vr)

  use declarations
  implicit none

  integer:: i, icheck, io, iord, is, ispin, iv, iw, j, k, nicm, npoint, npsom, nspin, nvois

  integer, dimension(npsom,nvois):: ivois

  logical Base_hexa

  real(kind=db):: fac
  real(kind=db), dimension(nvois):: cgrad
  real(kind=db), dimension(npoint,nspin):: Vr
  real(kind=db), dimension(nicm,3,nspin):: gradvr

  gradvr(:,:,ispin) = 0._db

  iv = 0
  do io = 1,iord / 2
    do k = 1,3
      do is = 1,2
        iv = iv + 1
        do i = 1,npoint
          j = ivois(i,iv)
          if( j == 0 .or. j > npoint ) j = i

! Pour le calcul du gradient, on neglige la rotation eventuelle de l'axe
! de spin (y compris quand elle est due a une symetrie).
          gradvr(i,k,ispin) = gradvr(i,k,ispin) + cgrad(iv) * Vr(j,ispin)

          if( Base_hexa .and. k == 2 ) then
            iw = iv + 4
            j = ivois(i,iw)
            if( j == 0 .or. j > npoint ) j = i
            gradvr(i,k,ispin) = gradvr(i,k,ispin) + cgrad(iw) * Vr(j,ispin)
          endif

        end do
      end do
    end do
  end do

  if( icheck > 1 ) then
    write(3,130)
    fac = rydb / bohr
    write(3,140) (i, gradvr(i,1:3,ispin)*fac, i = 1,npoint)
  endif

  return
  130 format(/4x,'i   dvr/dx    dvr/dy    dvr/dz   (eV/Angstroem)')
  140 format(i5,3f10.5)
end

!***********************************************************************

subroutine Potential_writing(Delta_Eseuil,dV0bdcf,E_cut,E_Fermi,Ecineticmax,Epsii, &
                             Epsii_moy,Hubb_abs,Hubb_diag_abs,lmax_abs,m_hubb,ninitr,nlm_pot,nrato_abs,nspin,nspinp, &
                             rato_abs,rhoato_abs,Rmtg_abs,Rmtsd_abs,rsato_abs,rsbdc,V_hubb_abs,V0muf,Vcato_abs,Vhbdc,Vxcato_abs, &
                             VxcbdcF)

  use declarations
  implicit none

  integer:: i, ir, isp, jsp, lmax_abs, m, m_hubb, ninitr, nlm_pot, nrato_abs, nspin, nspinp

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb_abs

  logical:: Hubb_abs, Hubb_diag_abs

  real(kind=db):: Delta_Eseuil, E_cut, E_Fermi, Ecineticmax, Epsii_moy, Rmtg_abs, Rmtsd_abs, rsbdc, V0muf, Vhbdc
  real(kind=db), dimension(nspin):: dV0bdcf, VxcbdcF
  real(kind=db), dimension(ninitr):: Epsii
  real(kind=db), dimension(nrato_abs):: rato_abs, rsato_abs
  real(kind=db), dimension(nrato_abs,nspin):: rhoato_abs
  real(kind=db), dimension(nrato_abs,nlm_pot):: Vcato_abs
  real(kind=db), dimension(nrato_abs,nlm_pot,nspin):: Vxcato_abs

  write(3,110)
  write(3,'(/A)') '   nspin nspinp nrato_abs nlm_pot m_hubb lmax_abs ninitr'
  write(3,'(2i7,i9,4i7)') nspin, nspinp, nrato_abs, nlm_pot, m_hubb, lmax_abs, ninitr

  if( nspin == 2 ) then
    write(3,120)
  else
    write(3,130)
  endif
  write(3,'(1p,10e19.11)') E_cut*Rydb, E_Fermi*Rydb, V0muf*Rydb, Vhbdc*Rydb, VxcbdcF(:)*Rydb, dV0bdcf(:)*Rydb, rsbdc, &
                         Ecineticmax*Rydb

  write(3,'(/A)') '  Hubb_abs  Hubb_abs_diag'
  write(3,*) Hubb_abs, Hubb_diag_abs
  if( Hubb_abs ) then
    write(3,'(/A)') '  V_Hubb_abs'
    do isp = 1,nspinp
      do jsp = 1,nspinp
        do m = -m_hubb,m_hubb
        ! it is the transpose but also in the reading !
          write(3,'(1p,28e19.11)') V_hubb_abs(-m_hubb:m_hubb,m,isp,jsp)*Rydb
        end do
      end do
    end do
  endif

  write(3,'(/A)') '       Rmtg_abs          Rmtsd_abs'
  write(3,'(1p,5e19.11)') Rmtg_abs*bohr, Rmtsd_abs*bohr

  write(3,140) ( 'Epsii_', i, i = 1,min(ninitr,9) ), ( 'Epsii_', i, i = 10,ninitr )
  write(3,'(1p,8e19.11)') Epsii_moy*Rydb, Delta_Eseuil*Rydb, Epsii(:)*Rydb

  if( nspin == 2 ) then
    write(3,'(/7x,A)') &
    'Radius             Vcato            Vxcato_up          Vxcato_dn            rhoato_up          rhoato_dn            Rsato'
  else
    write(3,'(/A)') '       Radius             Vcato              Vxcato              rhoato             Rsato'
  endif

  do ir = 1,nrato_abs
    write(3,'(1p,250e19.11)') rato_abs(ir)*bohr, Vcato_abs(ir,:)*Rydb, Vxcato_abs(ir,:,:)*Rydb, rhoato_abs(ir,:), &
                              rsato_abs(ir)
  end do

  return
  110 format(/' ---- Potential data --',100('-'))
  120 format(/'        E_cut             E_Fermi            V0muf              Vhbdc            VxbdcF_up          VxcbdcF_dn', &
              '         dV0bdcF_up         dV0bdcF_dn            rsbdc          E_kin_max')
  130 format(/'        E_cut             E_Fermi            V0muf              Vhbdc             VxcbdcF             dV0bdcF', &
              '             rsbdc           E_kin_max')
  140 format(/6x,'Epsii_moy',9x,'Delta_Eseuil',3x,9(6x,a6,i1,6x),20(5x,a6,i2,6x))
end

