! FDMNES subroutines

!***********************************************************************

! Sousprogramme qui donne la dimension des etats initiaux.
! ninitl : nombre total d'etats initiaux
! ninit1 : nombre d'etat initiaux du premier seuil

subroutine dim_init(jseuil,lseuil,nbseuil,ninit1,ninitl)

  use declarations
  implicit none

  integer,intent(in):: jseuil, lseuil, nbseuil
  integer,intent(out):: ninitl, ninit1

  if( nbseuil == 1 ) then

    select case( jseuil )
      case( 1 )           ! seuil K, L1, M1 ...
        ninit1 = 2
      case( 2 )           ! seuil L2, M2
        ninit1 = 2
      case( 3 )           ! seuil L3, M3
        ninit1 = 4
      case( 4 )           !       M4, N4
        ninit1 = 4
      case( 5 )           !       M5
        ninit1 = 6
      case( 6 )           !       N6
        ninit1 = 6
      case( 7 )           !       N7
        ninit1 = 8
    end select
    ninitl = ninit1

  else

    select case( lseuil )
      case( 0 )           ! K ou L1
        ninit1 = 2
        ninitl = 2
      case( 1 )           ! L23, M23
        ninit1 = 2
        ninitl = 6
      case( 2 )           ! M45
        ninit1 = 4
        ninitl = 10
      case( 3 )           ! N67
        ninit1 = 6
        ninitl = 14
    end select

  end if

  return
end

!***********************************************************************

! Sous-routine qui calcule les m et les coefficients qui multiplient les
! Ylm pour un certain etat |j,mj>

subroutine coef_init(coef_g,is_g,ninitl,jseuil,lseuil,m_g,nbseuil)

  use declarations
  implicit none

  integer,intent(in):: ninitl, jseuil, lseuil, nbseuil
  integer,dimension(ninitl),intent(out):: is_g
  integer,dimension(ninitl,2),intent(out):: m_g

  integer i, initl, isinitl, iseuil, is, isp, li, mi, ni

  real(kind=db) Ci2, J_initl, Jz, spin
  real(kind=db),dimension(ninitl,2),intent(out):: coef_g

! il existe une redondance entre la variable locale ninitl et celle
! externe ninit
  m_g(:,:) = 1000
  coef_g(:,:) = 0._db

  initl = 0

  select case(jseuil)
    case(1,3,5,7)
      isinitl = 1
    case default
      isinitl = -1
  end select

  do iseuil = 1, nbseuil

    if( iseuil == 2 ) isinitl = - isinitl
    li = lseuil
    J_initl = li + 0.5_db * isinitl
    ni = nint( 2*J_initl  + 1 )

! Boucle sur les etats initiaux |j,mj>, dont la multiplicite est 2*j+1

    do i = 1,ni

      initl = initl + 1
      Jz = - J_initl + i - 1
      is_g(initl) = isinitl

      do isp = 1,2

        if( isp == 1 ) then
          spin = 0.5_db
          is = isinitl
        else
          spin = -0.5_db
          is = 1
        endif

        mi = nint( Jz - spin )
        if( abs(mi) > li ) cycle

        Ci2 = ( li + 0.5_db + 2*spin*isinitl*Jz ) / ( 2*li + 1 )
        if( Ci2 < eps6 ) cycle

        coef_g(initl,isp) = is * sqrt( Ci2 )
        m_g(initl,isp) = mi

      end do

    end do

  end do

  return
end

!*********************************************************************

! Recherche des atomes equivalents a l'absorbeur et calcul de leur symetrie relative.

subroutine symsite(absauto,angxyz,Atom_with_axe,Atom_nonsph,Atom_nsph_e,Axe_atom_gr,axyz,Base_ortho,Cubmat,dcosxyz, &
        Doping,Extract,Flapw,iabsm,icheck,igr_i,igr_is,igr_proto,itype, &
        Magnetic,Matper,Memory_save,n_atom_proto,n_multi_run_e,neqm, &
        ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,numat_abs,popats,posn,Struct)

  use declarations
  implicit none

  integer:: i, ia, ib, icheck, igr, igra, igrb, io, irev, is, isp, it, ita, itb, jgr, js, long, n_atom_proto, n_multi_run_e, na, &
     nabs, nb, nla, nlb, neq, neqm, ngroup, ngroup_d, ngroup_m, nlatm, nspin, ntype, numat_abs, numat_it, numat_jgr

  character(len=5):: struct
  character(len=9):: nomsym
  character(len=11):: nomt

  integer, dimension(0:ntype):: nlat, numat
  integer, dimension(ngroup):: igr_i, igr_is, igr_proto, itype
  integer, dimension(n_multi_run_e):: iabsm
  integer, dimension(ngroup):: igreq, itypegen
  integer, dimension(ntype):: itequ
  integer, dimension(:), allocatable:: iabsmm, isymq

  logical:: absauto, Atom_nonsph, base_ortho, Doping, Extract, Flapw, Magnetic, matper, Memory_save, mspinor
  logical, dimension(ngroup):: Atom_nsph_e, ok
  logical, dimension(0:ngroup_m):: Atom_with_axe
  logical, dimension(:), allocatable:: Far_atom

  real(kind=db):: dist, distm_neq, dpop, pp1, pp2, rad, vnorme

  real(kind=db), dimension(3):: angxyz, Axe_atom_c, axyz, dcosxyz, ps, spini, vspin, wspin
  real(kind=db), dimension(3,3):: Cubmat, Cubmati, matopsym
  real(kind=db), dimension(3,ngroup):: posn
  real(kind=sg), dimension(3,ngroup):: pos, posg, poss
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_gr, Axe_atom_s
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats

  call invermat(cubmat,cubmati)

  allocate( Far_atom(ngroup) )

  if( Doping ) then
    ngroup_d = ngroup - 1
  else
    ngroup_d = ngroup
  endif

  neqm = 0

  rad = pi / 180
  dcosxyz(:) = 2 * cos( angxyz(:) * rad )
  if( abs( dcosxyz(1) ) < eps10 .and. abs( dcosxyz(2) ) < eps10 .and. abs( dcosxyz(3) ) < eps10 ) then
    Base_ortho = .true.
  else
    Base_ortho = .false.
  endif

  ok(:) = .false.

  do it = 1,ntype
    itequ(it) = it
  end do

  do ita = 1,ntype
    na = numat(ita)
    nla = nlat(ita)
    do igra = 1,ngroup_d
      if( abs(itype(igra)) == ita ) exit
    end do
    if( igra > ngroup_d ) cycle
    do itb = ita+1,ntype
!!
      if( itequ(itb) /= itb ) cycle
      if( FLAPW .and. .not. Extract ) cycle
!!
      nb = numat(itb)
      if( nb /= na ) cycle
      nlb = nlat(itb)
      if( nlb /= nla ) cycle
      do igrb = 1,ngroup_d
        if( abs(itype(igrb)) == itb ) exit
      end do
      if( igrb > ngroup_d ) cycle

      dpop = 0._db
      do io = 1,nlb
        do isp = 1,nspin
          dpop = dpop + abs( popats(igra,io,isp) - popats(igrb,io,isp) )
        end do
      end do
      if( dpop < eps6 ) then
        itequ(itb) = ita
      elseif( nspin > 1 ) then
        dpop = 0._db
        do io = 1,nlb
          do isp = 1,nspin
            dpop = dpop + abs( popats(igra,io,isp) - popats(igrb,io,nspin-isp+1) )
          end do
        end do
        if( dpop < eps6 ) itequ(itb) = -ita
      endif

    end do
  end do

  do igr = 1,ngroup_d
    it = abs( itype(igr) )
    itypegen(igr) = itequ(it)
  end do

  Atom_with_axe(1:ngroup_m) = Atom_nsph_e(1:ngroup_m)
  if( Magnetic ) then
    if( Flapw ) then
      Atom_with_axe(:) = .true.
    else
      do igr = 1,ngroup_d
        it = abs( itype(igr) )
        if( nlat(it) > 0 ) then
          dpop = sum( abs( popats(igr,1:nlat(it),nspin) - popats(igr,1:nlat(it),1) ) )
          if( dpop > eps6 ) Atom_with_axe(igr) = .true.
        endif
      end do
    endif
  endif

  n_atom_proto = 0

! Liste des atomes absorbeurs
  if( Absauto ) then
    nabs = 0
    do igr = 1,ngroup_d
      if( numat( abs( itype(igr) ) ) == numat_abs ) nabs = nabs+1
    end do
    allocate( iabsmm(nabs) )
    ia = 0
    do igr = 1,ngroup_d
      if( numat( abs( itype(igr) ) ) /= numat_abs ) cycle
      ia = ia + 1
      iabsmm(ia) = igr
    end do
  else
    nabs = n_multi_run_e
    allocate( iabsmm(n_multi_run_e) )
    iabsmm(:) = iabsm(:)
  endif

! Parametres servant a limiter le calcul du nombre d'atomes non
! equivalents
  if( Memory_save ) then

!        distm_neq = 2.6_db / bohr
    distm_neq = 1.0_db / bohr

    do igr = 1,ngroup_d

      Far_atom(igr) = .true.
      do ia = 1,nabs
        ps(:) = posn(:,iabsmm(ia)) - posn(:,igr)
        if( Matper ) then
          where( ps > 0.5_db ) ps = ps - 1._db
          where( ps < -0.5_db ) ps = ps + 1._db
        endif
        ps(:) = ps(:) * axyz(:)
        dist = vnorme(Base_ortho,dcosxyz,ps)
        if( dist <= distm_neq ) then
          Far_atom(igr) =  .false.
          exit
        endif
      end do

    end do

  else

    Far_atom(:) = .false.

  endif

  allocate( isymq(ngroup_d) )

  boucle_jgr: do jgr = 1,ngroup_d

    if( ok(jgr) .or. Far_atom(jgr) ) cycle

    numat_jgr = numat( abs( itype(jgr) ) )

    n_atom_proto = n_atom_proto + 1

    ps(:) = posn(:,jgr)
    do ia = 1,ngroup_d
      posg(:,ia) = posn(:,ia) - ps(:)
    end do

    if( matper ) call inmesh_sg(posg,ngroup_d)

    neq = 1
    isymq(1) = 1
    igreq(1) = jgr
    ok(jgr) = .true.

    boucle_igr: do igr = 1,ngroup_d

      if( ok(igr) .or. Far_atom(igr) ) cycle boucle_igr

      if( abs(itypegen(igr)) /= abs(itypegen(jgr))) cycle boucle_igr

      do ia = 1,ngroup_d
        pos(:,ia) = posg(:,ia) - posg(:,igr)
      end do
      if( matper ) call inmesh_sg(pos,ngroup_d)

! Recherche de la  symetrie qui rend les atomes equivalents.
      boucle_is: do is = 1,nopsm

        call opsym(is,matopsym)
        if( struct /= 'cubic' ) then
          matopsym = matmul(matopsym,cubmat)
          matopsym = matmul(cubmati,matopsym)
        endif

        do ia = 1,ngroup_d
          ps(:) = posg(:,ia) * axyz(:)
          ps = matmul( matopsym, ps )
          poss(:,ia) = ps(:) / axyz(:)
        end do
        if( matper ) call inmesh_sg(poss,ngroup_d)

        if( Magnetic .or. Atom_nonsph ) then
          call opsym(is,matopsym)
          do ia = 1,ngroup_d
            Axe_atom_c(:) = axyz(:) * Axe_atom_gr(:,ia)
            Axe_atom_c =  matmul( cubmat, Axe_atom_c )
            if( abs(Axe_atom_c(1)) < eps6 ) then
              wspin(1) = 1._db;  wspin(2:3) = 0._db
            else
              wspin(3) = 1._db;  wspin(1:2) = 0._db
            endif
            call prodvec( vspin, Axe_atom_c, wspin )
            vspin(:) = vspin(:) / sum( vspin(:)**2 )
            call prodvec( wspin, Axe_atom_c, vspin )
            wspin(:) = wspin(:) / sum( wspin(:)**2 )
            vspin = matmul( matopsym, vspin )
            wspin = matmul( matopsym, wspin )
            call prodvec(spini,vspin,wspin)
            spini =  matmul( cubmati, spini )
            spini(:) = spini(:) / axyz(:)

            Axe_atom_s(:,ia) = spini(:)

         end do
        endif

        boucle_irev: do irev = 1,nspin    ! 1 identite,
! 2 renversement du temps
          boucle_ia: do ia = 1,ngroup_d

            boucle_ib: do ib = 1,ngroup_d

              if( abs( pos(1,ia) - poss(1,ib) ) > epspos .or. abs( pos(2,ia) - poss(2,ib) ) > epspos .or. &
                  abs( pos(3,ia) - poss(3,ib) ) > epspos ) cycle boucle_ib

              if( abs(itypegen(ia)) /= abs(itypegen(ib)) ) cycle boucle_ib

              if( .not. ( Magnetic .or. Atom_nonsph ) ) cycle boucle_ia

              if( .not. ( Atom_with_axe(ia) .or. Atom_with_axe(ib))) cycle boucle_ia

              if( ( Atom_with_axe(ia) .and. .not. Atom_with_axe(ib)) .or. &
                ( .not. Atom_with_axe(ia) .and. Atom_with_axe(ib) )) cycle boucle_is

              mspinor = itypegen(ia) == itypegen(ib)
              if( Atom_with_axe(ia) ) then

                if( ( irev == 1 .and. mspinor ) .or. ( irev == 2 .and. .not. mspinor ) ) then
                  if(     abs( Axe_atom_gr(1,ia) - Axe_atom_s(1,ib) ) > epspos  &
                     .or. abs( Axe_atom_gr(2,ia) - Axe_atom_s(2,ib) ) > epspos  &
                     .or. abs( Axe_atom_gr(3,ia) - Axe_atom_s(3,ib) ) > epspos ) cycle boucle_irev
                else
                  if(     abs( Axe_atom_gr(1,ia) + Axe_atom_s(1,ib) ) > epspos  &
                     .or. abs( Axe_atom_gr(2,ia) + Axe_atom_s(2,ib) ) > epspos  &
                     .or. abs( Axe_atom_gr(3,ia) + Axe_atom_s(3,ib) ) > epspos ) cycle boucle_irev
                endif

              else

                 pp1 = sum( abs( Axe_atom_gr(:,ia) - Axe_atom_s(:,ib) ) )
                 pp2 = sum( abs( Axe_atom_gr(:,ia) + Axe_atom_s(:,ib) ) )
                 if( pp1 > epspos .and. pp2 > epspos ) cycle boucle_irev
              endif

              js = ia

              cycle boucle_ia

            end do boucle_ib

            cycle boucle_is

          end do boucle_ia

! Cette symetrie est la bonne
          if( irev == 2 ) then
            js = - is
          else
            js = is
          endif

          neq = neq + 1
          ok(igr) = .true.
          igreq(neq) = igr
          isymq(neq) = js

          exit boucle_is

        end do boucle_irev

      end do boucle_is

    end do boucle_igr

! On prend les axes de spins en reference au premier atome du groupe
    if( Magnetic .or. Atom_nonsph ) then

      do i = 2,neq

        is = abs( isymq(i) )
        call opsym(is,matopsym)

        Axe_atom_c(:) = axyz(:) * Axe_atom_gr(:,igreq(1))
        Axe_atom_c =  matmul( cubmat, Axe_atom_c )
        if( abs(Axe_atom_c(1)) < eps6 ) then
          wspin(1) = 1._db;  wspin(2:3) = 0._db
        else
          wspin(3) = 1._db;  wspin(1:2) = 0._db
        endif
        call prodvec( vspin, Axe_atom_c, wspin )
        vspin(:) = vspin(:) / sum( vspin(:)**2 )
        call prodvec( wspin, Axe_atom_c, vspin )
        wspin(:) = wspin(:) / sum( wspin(:)**2 )
        vspin = matmul( matopsym, vspin )
        wspin = matmul( matopsym, wspin )
        call prodvec(spini,vspin,wspin)
        spini =  matmul( cubmati, spini )
        spini(:) = spini(:) / axyz(:)

        ia = igreq(i)
        Axe_atom_gr(:,ia) = (isymq(i) / is ) * spini(:)

      end do

    endif

    neqm = max(neqm,neq)
    do i = 1,neq
      igr_proto(igreq(i)) = n_atom_proto
      igr_i(igreq(i)) = i
      igr_is(igreq(i)) = isymq(i)
    end do

    if( icheck > 0 ) then
      if( n_atom_proto == 1 ) write(3,110)
      write(3,120) n_atom_proto, numat_jgr, neq
      do ia = 1,neq
        nomt = '         '
        nomt(3:11) = nomsym( abs(isymq(ia)) )
        if(isymq(ia) < 0 ) then
          long = len_trim( adjustl( nomsym(abs(isymq(ia))) ) )
          nomt(10-long:11-long) = 'T.'
        endif
        igra = igreq(ia)
        write(3,130) igreq(ia), posn(:,igra), nomt, isymq(ia)
      end do
    endif

  end do boucle_jgr

  deallocate( Far_atom )
  deallocate( isymq )

  if( Memory_save ) then

    do jgr = 1,ngroup_d

      if( ok(jgr) ) cycle

      ok(jgr) = .true.
      n_atom_proto = n_atom_proto + 1
      neq = 1

      igr_proto(jgr) = n_atom_proto
      igreq(1) = jgr
      igr_i(jgr) = neq
      igr_is(jgr) = 1

      do igr = jgr+1,ngroup_d
        if( ok(igr ) .or. itypegen(igr) /= itypegen(jgr) ) cycle

        ok(igr) = .true.
        neq = neq + 1

        igreq(neq) = igr
        igr_proto(igr) = n_atom_proto
        igr_i(igr) = neq
        igr_is(igr) = 1
      end do
      neqm = max(neqm,neq)

      if( icheck > 0 ) then
        numat_it = numat( abs( itype(jgr) ) )
        write(3,115) n_atom_proto, numat_it, neq
        write(3,125) igreq(1:neq)
      endif

    end do

  endif

  if( Doping ) then
    igr_proto(ngroup) = n_atom_proto + 1
  endif

  deallocate( iabsmm )

  return
  110 format(/' ---- Symsite ------',100('-'))
  115 format(/3x,'ipr =',i3,', Z =',i3,', natomsym =',i6//,'  igr...')
  120 format(/3x,'ipr =',i3,', Z =',i3,', natomsym =',i6//, '  igr      posx     posy     posz          sym    code')
  125 format(20i6)
  130 format(i5,2x,3f9.5,2x,a11,3x,129i3)
end

!*********************************************************************

subroutine inmesh(pos,ngroup)

  use declarations
  implicit none

  integer:: igr, k, m, ngroup

  real(kind=db) pos(3,ngroup)

  do igr = 1,ngroup
    do k = 1,3
      if( pos(k,igr) > 1._db - eps10 ) then
        m = int( pos(k,igr) + eps10 )
        pos(k,igr) = pos(k,igr) - m
      elseif( pos(k,igr) < - eps10 ) then
        m = 1 - int( pos(k,igr) + eps10 )
        pos(k,igr) = pos(k,igr) + m
      endif
    end do
  end do

  return
end

!*********************************************************************

subroutine inmesh_sg(pos,ngroup)

  use declarations
  implicit none

  integer:: igr, k, m, ngroup

  real(kind=sg) pos(3,ngroup)

  do igr = 1,ngroup
    do k = 1,3
      if( pos(k,igr) > 1._sg - eps10 ) then
        m = int( pos(k,igr) + eps10 )
        pos(k,igr) = pos(k,igr) - m
      elseif( pos(k,igr) < - eps10 ) then
        m = 1 - int( pos(k,igr) + eps10 )
        pos(k,igr) = pos(k,igr) + m
      endif
    end do
  end do

  return
end

!*********************************************************************

! Sousprogramme elaborant la grille en energie pour le calcul de XANES

subroutine grille_xanes(eeient,eimag,eimagent,egamme,energ,icheck,lin_gam,ngamme,neimagent,nenerg)

  use declarations
  implicit none

  integer:: icheck, i, ie, lin_gam, ngamme, neimagent, nenerg, ngc

  real(kind=db):: de, def, p1, r
  real(kind=db), dimension(neimagent) :: eeient, eimagent
  real(kind=db), dimension(nenerg) :: energ, eimag
  real(kind=db), dimension(ngamme) :: egamme

! Calcul de la grille en energie.
  energ(1) = egamme(1)

  if( lin_gam == 1 ) then            ! 'rangel'
    def = 10._db / rydb
    do ie = 2,nenerg
      r = 1 + energ(ie-1) / def
      r = max( r, 0.25_db )
      de = sqrt( r ) * egamme(2)
      energ(ie) = energ(ie-1) + de
    end do
  else
    ngc = 2
    do ie = 2,nenerg
      energ(ie) = energ(ie-1) + egamme(ngc)
      if( energ(ie) >= egamme(ngc+1) - eps10 ) ngc = ngc + 2
    end do
  endif

  eimag(1:nenerg) = 0._db
  if( neimagent > 0 ) then
    do ie = 1,nenerg
      if( energ(ie) >= eeient(neimagent) ) then
        eimag(ie) = eimagent(neimagent)
      elseif( energ(ie) <= eeient(1) ) then
        eimag(ie) = eimagent(1)
      else
        do i = 2,neimagent
          if( eeient(i) >= energ(ie) ) exit
        end do
        p1 = ( energ(ie) - eeient(i-1) ) / ( eeient(i) - eeient(i-1) )
        eimag(ie) = p1 * eimagent(i) + (1 - p1) * eimagent(i-1)
      endif
    end do
  endif

  if( icheck > 1 .and. neimagent > 0 ) then
    write(3,110)
    write(3,120)
    do ie = 1,nenerg
      write(3,130) energ(ie)*rydb, eimag(ie)*rydb
    end do
  endif

  return
  110 format(/' ---- Grille_xanes -',100('-'))
  120 format(/'   energie    eimag      en eV')
  130 format(4f9.3)
end

!***********************************************************************

function extract_green(nom_fich_extract)

  use declarations
  implicit none

  integer:: istat, l

  character(len=132) mot, nom_fich_extract

  logical green_s, extract_green

  open(1, file = nom_fich_extract, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(nom_fich_extract,istat,1)

  do l = 1,100000
    read(1,'(A)' ) mot
    if( mot(2:20) == 'Multiple scattering' ) then
      green_s = .true.
      exit
    elseif( mot(2:25) == 'Finite difference method' ) then
      green_s = .false.
      exit
    endif
  end do

  Close(1)

  extract_green = green_s

  return
end

!***********************************************************************

subroutine init_run(Chargat,Charge_free,Chargm,Clementi,Com,Doping,Ecrantage,Flapw,Force_ecr,Hubb,iabsorbeur, &
      iabsorig,icheck,icom,igreq,iprabs,iprabs_nonexc,itabs,itype,itype_dop,itypepr,jseuil,lcoeur,lecrantage, &
      lseuil,lvval,mpinodes,mpirank,n_atom_proto,n_multi_run,n_orbexc,nbseuil,ncoeur,necrantage,neqm,ngreq,ngroup,nlat,nlatm, &
      nnlm,nomfich_s,Nonexc,nrato,nrato_dirac,nrato_lapw,nrm,nseuil,nspin,ntype,numat,numat_abs,nvval,pop_open_val, &
      popatm,popats,popatv,popval,psi_coeur,psii,psi_open_val,psival,r0_lapw,rato,rchimp,Relativiste,rho_coeur, &
      rhoit,rlapw,Rmt,Rmtimp,V_hubbard,nompsii)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: i, iabsorig, iabsorbeur, icheck, igr, ipr, iprabs, iprabs_nonexc, itabs, itype_dop, jseuil, l, lecrantage, lseuil, &
      mpinodes, mpirank, n_atom_proto, n_multi_run, n_orbexc, nbseuil, necrantage, neqm, ngroup, nlatm, nnlm, nrato_dirac, nrm, &
      nseuil, nspin, ntype, numat_abs

  character(len=132):: nomfich_s, nompsii
  character(len=35), dimension(0:ntype):: com

  integer, dimension(ngroup):: itype
  integer, dimension(nnlm):: lqnexc, nqnexc
  integer, dimension(0:n_atom_proto):: itypepr, ngreq
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(0:ntype):: icom, nlat, nrato, nrato_lapw, numat
  integer, dimension(2,0:ntype):: lcoeur, ncoeur
  integer, dimension(0:ntype,nlatm):: lvval, nvval

  logical:: Charge_free, Clementi, Doping, Flapw, Force_ecr, nonexc, relativiste
  logical, dimension(0:ntype):: hubb

  real(kind=db):: Chargm
  real(kind=db), dimension(0:n_atom_proto):: chargat
  real(kind=db), dimension(ngroup):: chargatg
  real(kind=db), dimension(0:ntype):: popatc, r0_lapw, rchimp, rlapw, rmt, rmtimp, V_hubbard
  real(kind=db), dimension(nspin):: ecrantage
  real(kind=db), dimension(nnlm,nspin):: popexc
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm
  real(kind=db), dimension(0:ntype,nlatm):: popatv
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats
  real(kind=db), dimension(0:ntype,nlatm,nspin):: popval
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(0:nrm,0:ntype):: rato, rhoit, rho_coeur
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(nrm,2):: psi_open_val
  real(kind=db), dimension(2):: pop_open_val
  real(kind=db), dimension(0:nrm,2,0:ntype):: psi_coeur

  if( n_multi_run > 1 .and. mpirank == 0 ) then
    l = len_trim(nomfich_s)
    nomfich_s(l+1:l+1) = '_'
    call ad_number(iabsorig,nomfich_s,132)
    write(6,'(/3x,A)') nomfich_s
  endif

  if( icheck > 0 ) then
    write(3,100) iabsorig
    write(3,110)
  endif

  boucle_1: do iprabs_nonexc = 1,n_atom_proto
    do i = 1,ngreq(iprabs_nonexc)
      if( abs(igreq(iprabs_nonexc,i)) == iabsorbeur ) exit boucle_1
    end do
  end do boucle_1
  if( Nonexc ) then
    iprabs = iprabs_nonexc
  else
    iprabs = 0
  endif

! Atomic electron density
  call atom(Clementi,com,icheck,icom,itype,jseuil,lcoeur,lseuil,lvval,mpirank,nbseuil,ncoeur,ngroup,nlat, &
        nlatm,nnlm,Nonexc,nrato,nrato_dirac,nrato_lapw,nrm,nseuil,nspin,ntype,numat,nvval,popatc,popats, &
        popatv,popexc,popval,psi_coeur,psii,psival,r0_lapw,rato,Relativiste,rho_coeur,rhoit,rlapw)

  call pop_group(Chargatg,Charge_free,Chargm,Flapw,icheck,itype,mpirank,ngroup,nlat,nlatm,nspin,ntype,numat,popatc,popats)

  do ipr = 1,n_atom_proto
    igr = igreq(ipr,1)
    Chargat(ipr) = Chargatg(igr)
    do l = 1,nlat(abs(itype(igr)))
      popatm(ipr,l,1:nspin) = popats(igr,l,1:nspin)
    end do
  end do

! Absorbing atom electron density
  call type_work(com,Doping,Ecrantage,force_ecr,hubb,iabsorbeur,icheck,icom,itabs,itype,itype_dop,jseuil,lcoeur,lecrantage, &
        lqnexc,lseuil,lvval,mpinodes,mpirank,n_orbexc,nbseuil,ncoeur,necrantage,ngroup,nlat,nlatm,nnlm,nompsii,Nonexc,nqnexc, &
        nrato,nrato_dirac,nrm,nseuil,nspin,ntype,numat,nvval,pop_open_val,popatc,popats,popatv,popexc,popval, &
        psi_coeur,psii,psi_open_val,psival,rato,rchimp,Relativiste,rho_coeur,rhoit,rmt,rmtimp,V_hubbard)

  call Pop_mod(Chargat,Chargm,Doping,Flapw,iabsorbeur,icheck,igreq,itabs,iprabs_nonexc,itype,itype_dop,itypepr,lqnexc,lvval, &
        n_atom_proto,n_orbexc,neqm,ngreq,ngroup,nlat,nlatm,nnlm,nqnexc,nspin,ntype,numat_abs,nvval,popatm,popexc)

  return
  100 format(/' ----------',110('-'),//' Absorbing atom',i4)
  110 format(/' ---- Init_run ------',100('-'))
end

!***********************************************************************

subroutine atom(Clementi,Com,icheck,icom,itype,jseuil,lcoeur,lseuil,lvval,mpirank,nbseuil,ncoeur,ngroup,nlat, &
        nlatm,nnlm,Nonexc,nrato,nrato_dirac,nrato_lapw,nrm,nseuil,nspin,ntype,numat,nvval,popatc,popats, &
        popatv,popexc,popval,psi_coeur,psii,psival,r0_lapw,rato,Relativiste,rho_coeur,rhoit,rlapw)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: icheck, igr, ipr, it, j, jseuil, k, l, lseuil, mpirank, n, n_orbexc, &
    nbseuil, ngroup, nlatm, nnlm, nr, nrato_dirac, nrm, nseuil, nspin, ntype

  real(kind=db):: dr, drmax, dx, rmax

  character(len=35), dimension(0:ntype) :: com

  integer, dimension(nnlm):: lqnexc, nqnexc
  integer, dimension(ngroup):: itype
  integer, dimension(0:ntype):: icom, nlat, nrato, nrato_lapw, numat
  integer, dimension(2,0:ntype):: lcoeur, ncoeur
  integer, dimension(0:ntype,nlatm):: lvval, nvval

  logical:: Clementi, Nonexc, Relativiste

  real(kind=db), dimension(0:ntype):: popatc, r0_lapw, rlapw
  real(kind=db), dimension(nrm):: ra
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats
  real(kind=db), dimension(nnlm,nspin):: popexc
  real(kind=db), dimension(0:ntype,nlatm) :: popatv
  real(kind=db), dimension(0:ntype,nlatm,nspin):: popval
  real(kind=db), dimension(0:nrm,0:ntype):: rato, rhoit, rho_coeur
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(nrm,2):: psi_open_val
  real(kind=db), dimension(2):: pop_open_val
  real(kind=db), dimension(0:nrm,2,0:ntype):: psi_coeur

  if( clementi ) then
    do it = 1,ntype
      if( icom(it) == 1 .and. numat(it) < 55 ) then
        com(it) = ' Clementi and Roetti'
        icom(it) = 2
      endif
    end do
  endif

  do it = 1,ntype
    if( numat(it) == 0 ) cycle
    select case( icom(it) )

      case(1)

! Valeurs pour l'atome excite non encore utilises
        lqnexc(:) = 0
        nqnexc(:) = 0

        call dirgen(icheck,it,0,jseuil,lcoeur,lqnexc,lseuil,lvval,mpirank,n_orbexc,nbseuil,ncoeur,nlat,nlatm, &
          nnlm,nonexc,nqnexc,nr,nrato_dirac,nrm,nseuil,nspin,ntype,nvval,pop_open_val,popatc,popatv, &
          popexc,popval,psi_coeur,psii,psi_open_val,psival,ra,rho_coeur,rhoit,numat(it),Relativiste)

        nrato(it) = nr
        rato(1:nr,it) = ra(1:nr)
        rato(0,it) = 0._db

        do igr = 1,ngroup
          if( abs( itype(igr) ) /= it ) cycle
          do l = 1,nlat(it)
            popats(igr,l,1:nspin) = popval(it,l,1:nspin)
          end do
        end do

      case(2)

        call clem(icheck,it,0,lseuil,lvval,mpirank,nbseuil,nlat,nlatm,nonexc,nr,nrm,nseuil,ntype,nvval,popatc, &
               popatv,psii,psival,ra,rhoit,numat(it))
        nrato(it) = nr
        rato(1:nr,it) = ra(1:nr)
        rato(0,it) = 0._db

      case(3)

        nrato(it) = nrato_lapw(it)
        if( nrato(it) > nrm .and. mpirank == 0 ) then
          call write_error
          do ipr = 3,9,3
            write(ipr,110) nrato(it), nrm
          end do
          stop
        endif
        dx = log( rlapw(it) / r0_lapw(it) ) / ( nrato_lapw(it) - 1 )
        do j = 1,nrato_lapw(it)
          rato(j,it) = r0_lapw(it) * exp( ( j - 1 ) * dx )
        enddo
        rato(0,it) = 0._db
! on prolonge
        rmax = 2.5_db / bohr
        drmax = 0.1_db / bohr

        boucle_r: do j = nrato_lapw(it)+1,nrm
          rato(j,it) = r0_lapw(it) * exp( ( j - 1 ) * dx )
          if( rato(j,it) > rmax ) then
            nrato(it) = j
            exit boucle_r
          endif
          dr = rato(j,it) - rato(j-1,it)
          if( dr > drmax ) then
            do k = j,nrm
              rato(k,it) = rato(k-1,it) + drmax
              if( rato(k,it) > rmax ) then
                nrato(it) = k
                exit boucle_r
              endif
            end do
          endif
        end do boucle_r

        if( rato(nrato(it),it) < 2.5_db/bohr ) then
          n = 1 + Int( ( rmax - rato(nrato(it),it) ) / drmax )
          call write_error
          do ipr = 3,9,3
            write(ipr,120) nrm + n, nrm
          end do
          stop
        endif

    end select

  end do     ! fin de la boucle sur les especes chimiques

  return
  110 format(/' nrato =',i6,' > nrm =',i6,' !', /' Increase the value of nrm using keyword nrato')
  120 format(/' The number of radius needs to be increased',/ ' nrato =',i6,' > nrm =',i6,' after extrapolation !', &
         /' Increase the value of nrm using keyword nrato')
end

!***********************************************************************

! Calcul des atomes selon Clementi et Roetti

!  References :
!   HS  F. Herman and S. Skillman "Atomic Structure Calculations"
!       Prentice Hall 1963
!   CR  E. Clementi and C. Roetti "Atomic Data and Nuclear Data Tables"
!       Vol.14,(3-4) [1974] p.177-478

!   NPRIN Principal quantum number
!   NAZIM Azimuthal quantum number
!   NCONF Number of electrons occupying this orbital
!   NTERM Number of radial basis functions
!   NEXP  Ni for the basis functions (See Clementi Roetti book)
!   ESP   Exponential constant for the basis function
!   CON   Coefficient which multiplies the basis function
!
subroutine clem(icheck,it,itabs,lseuil,lvval,mpirank,nbseuil,nlat,nlatm,Nonexc,nrato,nrm,nseuil,ntype,nvval,popatc, &
              popatv,psii,psival,rato,rhoit,Z)

  use declarations
  implicit none

  integer, parameter:: norbm = 11
  integer, parameter:: ntrm = 11

  integer:: i, icalc, icheck, io, io0, ipr, ir, it, itabs, iv, j, lseuil, mpirank, nbseuil, ncalc, nlatm, norb, norb0, nr, &
    nrato, nrm, nseuil, ntype, Z, Zp

  integer, dimension(norbm):: nazim, nazim0, nconf, nconf0, nprin, nprin0, nterm, nterm0
  integer, dimension(norbm,ntrm) :: nexp, nexp0
  integer, dimension(0:ntype):: nlat
  integer, dimension(0:ntype,nlatm):: lvval, nvval

  logical:: Nonexc

  real(kind=db):: alfa, arg, ch, delta, f, phi, r0, rmax, ra

  real(kind=db), dimension(6):: nfact
  real(kind=db), dimension(nrm):: rato
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(norbm,ntrm):: esp, esp0, con, con0, onor
  real(kind=db), dimension(0:ntype):: popatc
  real(kind=db), dimension(0:ntype,nlatm):: popatv
  real(kind=db), dimension(0:nrm,0:ntype):: rhoit
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival

  if( icheck > 1 ) write(3,100) it, Z

! Rayons
  r0 = 0.00005_db / bohr
  rmax = 10._db / bohr
  delta = 1.02_db
  rato(1) = r0
  ra = r0
  f = log( delta )
  do ir = 2,1000000
   ra = ra + r0 * exp( ( ir - 1) * f )
   if( ir > nrm ) cycle
   rato(ir) = ra
   if( rato(ir) > rmax ) exit
  end do
  nrato = ir
  nr = nrato
  if( nr > nrm .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,105) nr, nrm
    end do
    stop
  endif

! Factoriels
  Nfact(1) = sqrt( 2._db )
  do i = 2,6
    Nfact(i) = Nfact(i-1) * sqrt( 2._db * i * ( 2 * i - 1 ) )
  end do

! Pour l'absorbeur, on calcule l'atome excite (Z+1) et l'atome non
! excite
  if( it == itabs .and. .not. nonexc ) then
    ncalc = 2
  else
    ncalc = 1
  endif

  do icalc = 1,ncalc

    Zp = Z + icalc - 1 ! Approximation Z + 1

    if( Zp < 30 ) then
      call clempar1(Zp,norb,nprin,nazim,nconf,nterm,nexp,esp,con)
    elseif( Z < 50 ) then
      call clempar2(Zp,norb,nprin,nazim,nconf,nterm,nexp,esp,con)
    else
      call clempar3(Zp,norb,nprin,nazim,nconf,nterm,nexp,esp,con)
    endif

    if( it == itabs .and. icalc == 1 ) then
      if( nseuil == 0 ) then
        psii(:,:) = 0._db
      else
        ch = 0._db
        do io = 1,norb
          if( nprin(io) == nseuil .and. nazim(io) == lseuil ) exit
        end do
        do i = 1,Nterm(io)
          alfa =  ( 2 * esp(io,i) )**( nexp(io,i) + 0.5_db )
          onor(io,i) =  alfa / Nfact(nexp(io,i))
        end do
        do ir = 1,nr
          phi = 0._db
          do j = 1,nterm(io)
            arg = esp(io,j) * rato(ir)
            if( arg > 100._db ) cycle
            phi = phi + onor(io,j) * con(io,j) * rato(ir)**nexp(io,j) * exp( - arg )
          end do
          psii(ir,1) = phi
          psii(ir,nbseuil) = phi
          ch = ch + phi**2 * ( rato(ir+1) - rato(ir) )
        end do
        if( icheck > 1 ) write(3,120) io, nprin(io), nazim(io), nconf(io), ch
      endif
    endif

    if( it == itabs .and. .not. nonexc ) then
      if( icalc == 1 ) then
        norb0 = norb
        nprin0(1:norb) = nprin(1:norb)
        nazim0(1:norb) = nazim(1:norb)
        nconf0(1:norb) = nconf(1:norb)
        nterm0(1:norb) = nterm(1:norb)
        nexp0(1:norb,:) = nexp(1:norb,:)
        esp0(1:norb,:) = esp(1:norb,:)
        con0(1:norb,:) = con(1:norb,:)
      else
        do io = 1,norb
          if( nprin(io) == nseuil .and. nazim(io) == lseuil ) nconf(io) = nconf(io) - 1
          if( nprin(io) > nseuil .or. ( nprin(io) == nseuil .and. nazim(io) >= lseuil ) ) cycle
          do io0 = 1,norb0
            if( nprin0(io0) /= nprin(io) .or. nazim0(io0) /= nazim(io) ) cycle
            nconf(io) = nconf0(io0)
            nterm(io) = nterm0(io0)
            nexp(io,:) = nexp0(io0,:)
            esp(io,:) = esp0(io0,:)
            con(io,:) = con0(io0,:)
            exit
          end do
        end do
      endif
    endif

  end do

  do io = 1,nlat(it)
    do i = 1,norb
      if( nprin(i) == nvval(it,io) .and. nazim(i) == lvval(it,io)) exit
    end do
    if( i == norb+1 .and. mpirank == 0 ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,110) it, nvval(it,io), lvval(it,io)
      end do
      stop
    endif
    popatv(it,io) = 1._db * nconf(i)
  end do

  popatc(it) = 1._db * sum( nconf(1:norb) )
  do io = 1,nlat(it)
    popatc(it) = popatc(it) - popatv(it,io)
  end do

!  Calculate normalization coefficients
!  (see CR p.180, equation (6))

  do io = 1,norb
    do i = 1,Nterm(io)
      alfa =  ( 2 * esp(io,i) )**( nexp(io,i) + 0.5_db )
      onor(io,i) =  alfa / Nfact(nexp(io,i))
    end do
  end do

  rhoit(:,it) = 0._db
  do io = 1,norb
    do iv = 1,nlat(it)
      if( nprin(io) == nvval(it,iv) .and. nazim(io) == lvval(it,iv) ) exit
    end do
    if( iv == nlat(it)+1 ) iv = 0
    ch = 0._db
    do ir = 1,nr
      phi = 0._db
      do j = 1,nterm(io)
        arg = esp(io,j) * rato(ir)
        if( arg > 100._db ) cycle
        phi = phi + onor(io,j) * con(io,j) * rato(ir)**nexp(io,j) * exp( - arg )
      end do
      if( iv > 0 ) psival(ir,iv,it) = phi
      rhoit(ir,it) = rhoit(ir,it) + nconf(io) * ( phi / rato(ir) )**2
      ch = ch + phi**2 * ( rato(ir+1) - rato(ir) )
    end do
    if( icheck > 1 ) write(3,120) io, nprin(io), nazim(io), nconf(io), ch
  end do
  rhoit(1:nr,it) = rhoit(1:nr,it) / quatre_pi

  return
  100 format(/' it =',i2,'  Z =',i4)
  105 format(/' Number of radius =',i6,' > nrm =',i6,// ' Increase nrm using keyword nrato !'/)
  110 format(///' Orbital not found in the Clementi and Roetti bases :',/'   it =',i2,',  n =',i2,'  l =',i2)
  120 format(' io, nprin, nazim, nconf =',4i3,' ch =',f10.7)
end

!***********************************************************************

! Charge et population de chaque atome de la maille

subroutine pop_group(chargatg,Charge_free,Chargm,Flapw,icheck,itype,mpirank,ngroup,nlat,nlatm,nspin,ntype,numat,popatc,popats)

  use declarations
  implicit none

  integer:: icheck, igr, it, l, mpirank, ngroup, nlatm, nspin, ntype
  integer, dimension(0:ntype) :: nlat, numat
  integer, dimension(ngroup):: itype

  logical:: Charge_free, Flapw

  real(kind=db):: Chargm
  real(kind=db), dimension(ngroup):: chargatg
  real(kind=db), dimension(0:ntype):: popatc
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats

  if( icheck > 1 ) write(3,110)

  if( Flapw ) then
    chargatg(:) = 0
    return
  endif

  do igr = 1,ngroup
    it = abs( itype(igr) )
    if( nlat(it) > 0 ) then
        chargatg(igr) = numat(it) - popatc(it) - sum( popats(igr,1:nlat(it),1:nspin) )
    else
      chargatg(igr) = numat(it) - popatc(it)
    endif
  end do

  if( icheck > 1 ) then
    if( nspin == 1 ) then
      write(3,130)
    else
      write(3,140)
    endif
    do igr = 1,ngroup
      it = abs( itype(igr) )
      write(3,150) igr, chargatg(igr), ( popats(igr,l,1:nspin), l = 1,nlat(it) )
    end do
  endif

! Test sur la neutralite de la maille ou de la molecule
  Chargm = sum( chargatg( 1:ngroup ) )
  if( abs(Chargm) > 0.0001 .and. mpirank == 0 ) then
    write(3,120) Chargm
    write(6,120) Chargm
    if( .not. Charge_free ) then
      call write_error
      write(9,120) Chargm
      stop
    endif
  endif

  return
  110 format(/' ---- Pop_group ----',100('-'))
  120 format(/' Unit mesh or molecule charge =',f8.4)
  130 format(/'  igr    charge   popats(1)  popats(2)')
  140 format(/'  igr    charge   popats(1,up)  popats(1,dn) ...')
  150 format(i4,11f11.5)
end

!*********************************************************************

subroutine type_work(com,Doping,Ecrantage,force_ecr,hubb,iabsorbeur,icheck,icom,itabs,itype,itype_dop,jseuil,lcoeur,lecrantage, &
        lqnexc,lseuil,lvval,mpinodes,mpirank,n_orbexc,nbseuil,ncoeur,necrantage,ngroup,nlat,nlatm,nnlm,nompsii,nonexc,nqnexc, &
        nrato,nrato_dirac,nrm,nseuil,nspin,ntype,numat,nvval,pop_open_val,popatc,popats,popatv,popexc,popval, &
        psi_coeur,psii,psi_open_val,psival,rato,rchimp,relativiste,rho_coeur,rhoit,rmt,rmtimp,V_hubbard)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  character(len=35), dimension(0:ntype) :: com
  character(len=132) :: nompsii

  integer, dimension(ngroup):: itype
  integer, dimension(nnlm):: lqnexc, nqnexc
  integer, dimension(0:ntype):: icom, nlat, nrato, numat
  integer, dimension(0:ntype,nlatm):: lvval, nvval
  integer, dimension(2,0:ntype):: lcoeur, ncoeur
  logical, dimension(0:ntype):: hubb

  logical:: Doping, Force_ecr, nonexc, relativiste

  real(kind=db), dimension(0:ntype):: popatc, rchimp, rmt, rmtimp, V_hubbard
  real(kind=db), dimension(nspin):: ecrantage
  real(kind=db), dimension(nrm):: rr
  real(kind=db), dimension(nrm,nbseuil):: psii, psiit
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats
  real(kind=db), dimension(nnlm,nspin):: popexc
  real(kind=db), dimension(0:ntype,nlatm) :: popatv
  real(kind=db), dimension(0:ntype,nlatm,nspin):: popval
  real(kind=db), dimension(0:nrm,0:ntype):: rato, rhoit, rho_coeur
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(nrm,2):: psi_open_val
  real(kind=db), dimension(2):: pop_open_val
  real(kind=db), dimension(0:nrm,2,0:ntype):: psi_coeur

  if( Doping ) then
    itabs = itype_dop
  else
    itabs = abs( itype(iabsorbeur) )
  endif

  if( numat(itabs) == 0 .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,110) iabsorbeur
    end do
    stop
  endif

  itabs_t = itabs

  rmt(0) = rmt(itabs)
  numat(0) = numat( itabs )
  nlat(0) = nlat( itabs )
  nvval(0,:) = nvval(itabs,:)
  ncoeur(:,0) = ncoeur(:,itabs)
  lcoeur(:,0) = lcoeur(:,itabs)
  lvval(0,:) = lvval(itabs,:)
  popval(0,:,:) = popval(itabs,:,:)
  com(0) = com(itabs)
  icom(0) = icom(itabs)
  rchimp(0) = rchimp(itabs)
  rmtimp(0) = rmtimp(itabs)
  Hubb(0) = Hubb(itabs)
  V_hubbard(0) = V_hubbard(itabs)

  call Screendef(ecrantage,force_ecr,icheck,lecrantage,lqnexc,lseuil,lvval,mpirank,n_orbexc,necrantage, &
        nlat,nlatm,nnlm,nqnexc,nseuil,nspin,ntype,nvval,popexc,popval,numat(itabs))

! itabs est remis a la valeur initiale a la fin de la routine si on est en nonexc
  itabs = 0

! Il faut recalculer cette partie meme en nonexc pour avoir psii.
  psii(:,:) = 0._db

  if( icom(itabs) == 2 .or. nompsii == 'clementi' ) then

    call clem(icheck,itabs,itabs,lseuil,lvval,mpirank,nbseuil,nlat,nlatm,nonexc,nr,nrm,nseuil,ntype,nvval,popatc, &
          popatv,psii,psival,rr,rhoit,numat(itabs))

  elseif( icom(itabs) == 1 .or. nompsii == 'dirac' ) then

    call dirgen(icheck,itabs,itabs,jseuil,lcoeur,lqnexc, lseuil,lvval,mpirank,n_orbexc,nbseuil,ncoeur,nlat,nlatm, &
         nnlm,nonexc,nqnexc,nr,nrato_dirac,nrm,nseuil,nspin,ntype,nvval,pop_open_val,popatc,popatv, &
         popexc,popval,psi_coeur,psii,psi_open_val,psival, rr,rho_coeur,rhoit,numat(itabs),relativiste)

! popats pour l'absorbeur
    do igr = 1,ngroup
      if( abs( itype(igr) ) /= itabs ) cycle
      do l = 1,nlat(itabs)
        popats(igr,l,1:nspin) = popval(itabs,l,1:spin)
      end do
    end do

  elseif( icom(itabs) == 3 .and. nompsii /= 'clementi' .and. nompsii /= 'dirac' ) then
    if( mpirank == 0 ) then
      open(8, file = nompsii, status='old', iostat=istat)
      if( istat /= 0 ) call write_open_error(nompsii,istat,1)
      read(8,*) nr
      do ir = 1,nr
        read(8,*) rr(ir), psii(ir,1)
      end do
      if( nbseuil == 2 ) psii(:,nbseuil) = psii(:,1)
      close(8)
    endif
    if( mpinodes > 1 ) then
      call MPI_Bcast(nr,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(rr,nrm,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(psii,nrm*nbseuil,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    endif
  endif

! ici on fait les atributions pour l'absorbeur rato(:,0)

  if( icom(itabs) /= 3 ) then
    nrato(itabs) = nr
    rato(0,itabs) = 0._db
    rato(1:nr,itabs) = rr(1:nr)
  endif

  if( icom(itabs) == 3 .and. ( nompsii == 'clementi' .or. nompsii == 'dirac' ) ) then
    psiit(1:nr,:) = psii(1:nr,:)
! Interpolation pour avoir la fct d'onde initiale de coeur dans les rayons FLAPW
    psii(:,:) = 0._db
    do ir = 1,nrato(itabs_t)
      do irt = 2,nr
        if( rr(irt) > rato(ir,itabs_t) ) exit
      end do
      if( irt > nr ) exit
      p1 = ( rato(ir,itabs_t) -  rr(irt-1) ) / (rr(irt) - rr(irt-1))
      p2 = 1 - p1
      psii(ir,:) = p1*psiit(irt,:) + p2*psiit(irt-1,:)
    end do
  endif

  if( icheck > 0 ) then
    write(3,120)
    do it = 0,ntype
      if( nonexc .and. it == 0 ) cycle
      write(3,130) it, com(it), numat(it), ( nvval(it,l), lvval(it,l), popatv(it,l), l = 1,nlat(it))
    end do
    if( .not. nonexc ) write(3,140) necrantage, lecrantage, ecrantage(1:nspin)
  endif
  if( icheck > 1 ) then
    write(3,150)
    do ir = 1,nrato(itabs)
      write(3,160) rato(ir,itabs)*bohr, psii(ir,:)
    end do
  endif

  if( nonexc ) itabs = itabs_t

  return
  110 format(//' The absorbing atom cannot be the number',i3,' because it is an empty sphere !'//)
  120 format(/' Atom type',21x,'Z  n  l  popatv')
  130 format(i3,1x,a25,i3,8(2i3,f6.2))
  140 format(/' Default or imposed orbital screening :',/ &
              '      When default is used, if the screening orbital is full,',/ &
              '      it is the next one which is filled',/ &
              '   n_screening_orbital =',i2,/ &
              '   l_screening_orbital =',i2,/ &
              '   Screening =',2f6.3)
  150 format(/' Core wave function after interpolation' )
  160 format(f12.8,2f11.7)
end

!***********************************************************************

! Calcul de la configuration electronique de l'atome excite

subroutine Screendef(Ecrantage,Force_ecr,icheck,lecrantage,lqnexc,lseuil,lvval,mpirank,n_orbexc,necrantage, &
        nlat,nlatm,nnlm,nqnexc,nseuil,nspin,ntype,nvval,popexc,popval,Z)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer Z
  integer, dimension(nnlm):: lqnexc, nqnexc
  integer, dimension(0:ntype):: nlat
  integer, dimension(0:ntype,nlatm):: lvval, nvval

  logical force_ecr

  real(kind=db), dimension(nnlm):: nel, rqn
  real(kind=db), dimension(0:ntype,nlatm,nspin):: popval
  real(kind=db), dimension(nspin) :: dp, dpp, Ecrantage
  real(kind=db), dimension(nnlm,nspin):: popexc

  if( necrantage == 0 .or. .not. Force_ecr ) then
    select case( Z )
      case(1,2,3)
        necrantage = 2
        lecrantage = 0
      case(4,5,6,7,8,9)
        necrantage = 2
        lecrantage = 1
      case(10,11)
        necrantage = 3
        lecrantage = 0
      case(12,13,14,15,16,17)
        necrantage = 3
        lecrantage = 1
      case(18,19,20)
        necrantage = 4
        lecrantage = 0
      case(21,22,23,24,25,26,27,28,29,30)
        necrantage = 3
        lecrantage = 2
      case(31,32,33,34,35)
        necrantage = 4
        lecrantage = 1
      case(36,37,38)
        necrantage = 5
        lecrantage = 0
      case(39,40,41,42,43,44,45,46,47,48)
        necrantage = 4
        lecrantage = 2
      case(49,50,51,52,53)
        necrantage = 5
        lecrantage = 1
      case(54,55,56)
        necrantage = 6
        lecrantage = 0
      case(57,58,59,60,61,62,63,64,65,66,67,68,69,70,71)
        necrantage = 4
        lecrantage = 3
      case(72,73,74,75,76,77,78,79,80)
        necrantage = 5
        lecrantage = 2
      case(81,82,83,84,85)
        necrantage = 6
        lecrantage = 1
      case(86,87,88)
        necrantage = 7
        lecrantage = 0
      case(89,90,91,92,93,94,95,96,97,98,99,100,101)
        necrantage = 5
        lecrantage = 3
    end select
  endif

  irel = 0
  call config(Z,irel,n_coeur,n_orbexc,nnlm,nqnexc,lqnexc,rqn,nel)

  do ispin = 1,nspin
    popexc(1:n_orbexc,ispin) = nel(1:n_orbexc) / nspin
  end do

  if( Force_ecr ) then
    do l = 1,nlat(0)
      if( nvval(0,l) == necrantage .and. lvval(0,l) == lecrantage ) exit
    end do
    if( l > nlat(0) ) then
      nlat(0) = l
      nvval(0,l) = necrantage
      lvval(0,l) = lecrantage

      popval(0,l,1:nspin) = 0._db
      do ip = 1,n_orbexc
        if( nqnexc(ip) /= necrantage .or. lqnexc(ip) /= lecrantage ) cycle
        popval(0,l,1:nspin) = popexc(ip,1:nspin)
        exit
      end do

    endif
  endif

  it = 0

! On modifie la configuration electronique, compte tenu des entrees
  boucle_io: do io = 1,nlat(it)
    do ip = 1,n_orbexc
      if( nqnexc(ip) /= nvval(it,io) .or. lqnexc(ip) /= lvval(it,io)) cycle
      popexc(ip,1:nspin) = popval(it,io,1:nspin)
      cycle boucle_io
    end do
! On doit constuire une orbitale supplementaire
    n_orbexc = n_orbexc + 1
    if( n_orbexc > nnlm .and. mpirank == 0 ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,103)
      end do
      close(9)
      stop
    endif
    nqnexc(n_orbexc) = nvval(it,io)
    lqnexc(n_orbexc) = lvval(it,io)
    popexc(n_orbexc,1:nspin) = popval(it,io,1:nspin)
  end do boucle_io

! On enleve l'electron de coeur
  do io = 1,n_orbexc
    if( nqnexc(io) /= nseuil .or. lqnexc(io) /= lseuil ) cycle
    do ispin = 1,nspin
      popexc(io,ispin) = popexc(io,ispin) - 1._db / nspin
    end do
    exit
  end do

! On ajoute l'ecrantage

  if( Force_ecr ) then

    do io = 1,n_orbexc
      if( necrantage /= nqnexc(io) .or. lecrantage /= lqnexc(io) ) cycle
      popexc(io,1:nspin) = popexc(io,1:nspin) + ecrantage(1:nspin)
      exit
    end do

    if( io > n_orbexc ) then
      n_orbexc = n_orbexc + 1
      nqnexc(n_orbexc) = necrantage
      lqnexc(n_orbexc) = lecrantage
      popexc(n_orbexc,1:nspin) = ecrantage(1:nspin)
    endif

  else

    dp(1:nspin) = ecrantage(1:nspin)

    boucle_i: do i = 1,2
      do io = 1,n_orbexc
        if( nqnexc(io) < nseuil .or. ( nqnexc(io) == nseuil .and. lqnexc(io) <= lseuil ) ) cycle

        if( i == 1 .and. ( necrantage /= nqnexc(io) .or. lecrantage /= lqnexc(io) ) ) cycle
        if( i == 2 .and. ( necrantage == nqnexc(io) .and. lecrantage == lqnexc(io) ) ) cycle

        elmax = ( 2 + 4. * lqnexc(io) ) / nspin
        do ispin = 1,nspin
          dpp(ispin) = min(elmax - popexc(io,ispin), dp(ispin) )
          popexc(io,ispin) = popexc(io,ispin) + dpp(ispin)
          dp(ispin) = dp(ispin) - dpp(ispin)
        end do
        if( sum( dp(:) ) < eps10 ) exit boucle_i
        if( sum( popexc(io,:) ) < nspin * elmax - eps10 ) then
          dpp(1) = min(elmax - popexc(io,1), dp(nspin) )
          dpp(nspin) = min(elmax - popexc(io,nspin), dp(1) )
          popexc(io,:) = popexc(io,:) + dpp(:)
          dp(1) = dp(1) - dpp(nspin)
          dp(nspin) = dp(nspin) - dpp(1)
        endif
        if( sum( dp(:) ) < eps10 ) exit boucle_i
        if( i == 1 ) exit
      end do
    end do boucle_i

    if( dp(1) > eps10 .or. dp(nspin) > eps10 ) then
      nmax = 0
      do io = 1,n_orbexc
        nmax = max(nqnexc(io),nmax)
      end do
      if( nmax == 1 ) then
        iaug = 1
      else
        iaug = 0
        do io = 1,n_orbexc
          if( nqnexc(io) == nmax .and. lqnexc(io) == 1 ) then
            iaug = 1
            exit
          endif
        end do
      endif
      n_orbexc = n_orbexc + 1
      if( iaug == 1 ) then
        nqnexc(n_orbexc) = nmax + 1
        lqnexc(n_orbexc) = 0
      elseif( Z == 20 .or. Z == 38 .or.Z == 56 .or. Z == 88 ) then
        nqnexc(n_orbexc) = nmax - 1
        lqnexc(n_orbexc) = 2
      else
        nqnexc(n_orbexc) = nmax
        lqnexc(n_orbexc) = 1
      endif
      do ispin = 1,nspin
        if( dp(ispin) > eps10 ) then
          popexc(n_orbexc,ispin) = dp(ispin)
        else
          popexc(n_orbexc,ispin) = 0._db
        endif
      end do
    endif
  endif

  if( icheck > 1 ) then
    write(3,120)
    do io = 1,n_orbexc
      write(3,130) nqnexc(io), lqnexc(io), popexc(io,1:nspin)
    end do
  endif

  return
  103 format(///'   n_orb < nnlm ',// ' Increase nnlm in routine lecture.f !')
  120 format(/'  n  l popexc')
  130 format(2i3,2f7.3)
end

!***********************************************************************

! Allocation des tableaux itypepr et popatm pour l'excite

subroutine Pop_mod(chargat,chargm,Doping,Flapw,iabsorbeur,icheck,igreq,itabs,iprabs,itype,itype_dop,itypepr,lqnexc,lvval, &
        n_atom_proto,n_orbexc,neqm,ngreq,ngroup,nlat,nlatm,nnlm,nqnexc,nspin,ntype,numat_abs,nvval,popatm,popexc)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(nnlm):: lqnexc, nqnexc
  integer, dimension(ngroup):: itype
  integer, dimension(0:n_atom_proto):: itypepr, ngreq
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(0:ntype):: nlat
  integer, dimension(0:ntype,nlatm):: lvval, nvval

  logical:: Doping, flapw

  real(kind=db), dimension(0:n_atom_proto):: chargat
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm
  real(kind=db), dimension(nnlm,nspin):: popexc

  boucle_iprabs: do iprabs = 1,n_atom_proto
    do igrabs = 1,ngreq(iprabs)
      if( igreq(iprabs,igrabs) == iabsorbeur ) exit boucle_iprabs
    end do
  end do boucle_iprabs

  itypepr(:) = 0

  if( flapw ) then

    do ipr = 1,n_atom_proto
      itypepr(ipr) = abs( itype( igreq(ipr,1) ) )
    end do

  else

    itypepr(0) = 0
    do ipr = 1,n_atom_proto
      if( Doping .and. ipr == n_atom_proto ) then
        itypepr(ipr) = abs( itype_dop )
      else
        itypepr(ipr) = abs( itype( igreq(ipr,1) ) )
      endif
    end do

    igreq(0,1) = igreq(iprabs,igrabs)
    ngreq(0) = 1

    chargat(0) = numat_abs - sum( popexc(1:n_orbexc,1:nspin) )

    do io = 1,nlat(itabs)
      do jo = 1,n_orbexc
        if( nvval(itabs,io) /= nqnexc(jo) .or. lvval(itabs,io) /= lqnexc(jo) ) cycle
        popatm(0,io,1:nspin) = popexc(jo,1:nspin)
        exit
      end do
    end do

    chargm = chargm + chargat(0) - chargat(iprabs)

    if( icheck > 2 ) then
      write(3,120) iprabs
      if( nspin == 1 ) then
        write(3,130)
      else
        write(3,140)
      endif
      do ipr = 0,n_atom_proto
        nl = nlat( itypepr(ipr) )
        write(3,150) ipr, itypepr(ipr), chargat(ipr), ( popatm(ipr,l,1:nspin), l = 1,nl )
      end do
    endif

  endif

  return
  120 format(/' ---- Pop_mod ------',100('-')//,' iprabs =', i3)
  130 format(/'  ipr  it    charge   popatm(1)  popatm(2)')
  140 format(/'  ipr  it    charge   popatm(1,up)  popatm(1,dn) ...')
  150 format(2i4,11f11.5)
end

!***********************************************************************

function extract_E_cut(multi_run,nom_fich_extract)

  use declarations
  implicit none

  integer:: eof, i, istat, multi_run

  character(len=132) mot, nom_fich_extract

  real(kind=db):: E_cut, extract_E_cut

  open(1, file = nom_fich_extract, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(nom_fich_extract,istat,1)

  E_cut = -5._db / Rydb

  i = 0
  do
    read(1,'(A)' ) mot
    if( mot(2:15) /= 'Absorbing atom' ) cycle
    i = i + 1
    if( i == multi_run ) exit
  end do

  do
    read(1,'(A)',iostat=eof) mot
    if( eof /= 0 ) exit
    if( mot(2:6) /= 'E_cut' ) cycle
    backspace(1)
    read(1,'(8x,f11.5)') E_cut
    exit
  end do

  extract_E_cut = E_cut / Rydb

  return
end

!***********************************************************************

function extract_v0bdcF()

  use declarations
  implicit none

  character(len=132):: mot

  real(kind=db):: extract_v0bdcF, v0bdcF

  do
    read(1,'(A)' ) mot
    if( mot(3:7) == 'VmoyF' ) then
      backspace(1)
      read(1,'(9x,f10.3)') V0bdcF
      exit
    elseif( mot(6:10) == 'VmoyF' ) then
      backspace(1)
      read(1,'(12x,f10.3)') V0bdcF
      exit
    elseif( mot(6:11) == 'V0bdcF' ) then
      backspace(1)
      read(1,'(13x,f10.3)') V0bdcF
      exit
    endif
  end do

  extract_v0bdcF = v0bdcF / rydb

  return
end

!***********************************************************************

subroutine extract_Epsii(Core_resolved,Delta_Epsii,Delta_Eseuil, Epsii,Epsii_moy,icheck,nbseuil,ninit1,ninitl,ninitlr)

  use declarations
  implicit none

  integer:: i, icheck, nbseuil, ninit1, ninitl, ninitlr

  character(len=132) mot

  logical:: Core_resolved, Found

  real(kind=db):: Delta, Delta_Epsii, Delta_Eseuil, E, E_up, E_dn, Ep_moy, Epsii_moy
  real(kind=db), dimension(ninitlr):: Epsii

  Found = .false.
  Epsii(:) = 0._db

  if( icheck > 0 ) write(3,*)

  do
    read(1,'(A)' ) mot
    if( mot(5:9) == 'Epsii' ) then
      backspace(1)
      do i = 1,ninitl
        read(1,'(11x,f10.3)') E
        if( icheck > 0 ) write(3,110) E, i
        if( Core_resolved ) then
          Epsii(i) = E
        else
          if( i <= ninit1 ) then
            Epsii(1) = Epsii(1) + E / ninit1
          else
            Epsii(2) = Epsii(2) + E / ( ninitl - ninit1 )
          endif
        endif
      end do
      exit
    else

! Anciennes versions
      if( mot(3:7) == 'Epsii' ) then
        backspace(1)
        read(1,'(9x,f10.3)') E
        Found = .true.
      elseif( mot(6:12) == 'Epsii_up' ) then
        backspace(1)
        read(1,'(15x,f10.3,17x,f10.3)') E_up, E_dn
        E = ( E_up + E_dn ) / 2
        Found = .true.
      elseif( mot(6:10) == 'Epsii' ) then
        backspace(1)
        read(1,'(12x,f10.3)') E
        Found = .true.
      endif
      if( Found ) then
        Epsii(1:ninitlr) = E
        do i = 1,ninitl
          if( icheck > 0 ) write(3,110) E, i
        end do
        exit
      endif

    endif
  end do

  Close(1)

  Epsii(:) = Epsii(:) / rydb

  if( Core_resolved ) then
    if( ninit1 /= ninitl ) then
      Epsii_moy = sum( Epsii(ninit1+1:ninitl) ) / ( ninitl - ninit1 )
      Ep_moy = sum( Epsii(1:ninit1) ) / ninit1
    else
      Epsii_moy = sum( Epsii(1:ninitl) ) / ninitl
      Ep_moy = Epsii_moy
    endif
  else
    Epsii_moy = Epsii(ninitlr)
    Ep_moy = Epsii(1)
  endif
  Delta_Eseuil = Ep_moy - Epsii_moy

  if( Delta_Epsii < 100000._db .and. nbseuil > 1 ) then
    if( Core_resolved ) then
      Delta = Delta_Epsii - Delta_Eseuil
      Epsii(1:ninit1) = Epsii(1:ninit1) + Delta
    else
      Epsii(1) = Epsii(ninitlr) + Delta_Epsii
    endif
    Delta_Eseuil = Delta_Epsii
  endif

  if( .not. Core_resolved .and. icheck > 0 ) then
    if( ninitlr == 1 ) then
      Write(3,120) Epsii(:) * Rydb
    else
      Write(3,130) Epsii(:) * Rydb
    endif
  endif

  Close(1)

  return
  110 format(4x,'Epsii =',f10.3,' eV,  initial state',i3)
  120 format(/7x,'Epsii used =',f11.3,' eV')
  130 format(/7x,'Epsii used =',f11.3,'  and',f11.3,' eV')
end

!***********************************************************************

function extract_nenerg(multi_run,nom_fich_extract,Tddft)

  use declarations
  implicit none

  integer:: extract_nenerg, eof, i, ipr, istat, multi_run, nenerg

  character(len=132):: mot, nom_fich_extract

  logical:: Tddft

  open(1, file = nom_fich_extract, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(nom_fich_extract,istat,1)

  i = 0
  do
    read(1,'(A)' ) mot
    if( mot(2:15) /= 'Absorbing atom' ) cycle
    i = i + 1
    if( i == multi_run ) exit
  end do

  if( tddft ) then
    do while(eof==0)
      read(1,'(A)',iostat=eof) mot
      if( mot(2:12) == 'Cycle TDDFT' .or. mot(2:12) == 'Cycle Tddft') exit
    end do

    if( eof /= 0 ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,110)
      end do
      stop
    endif

  endif

  nenerg = 0
  boucle_je: do
    do
      read(1,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit boucle_je
      if( mot(2:12) == 'Cycle TDDFT' .or. mot(2:12) == 'Cycle Tddft' .or. mot(2:15) == 'Absorbing atom' ) exit boucle_je
      if( mot(7:11) == 'Coabs' ) nenerg = nenerg + 1
    end do
  end do boucle_je

  Close(1)

  extract_nenerg = nenerg

  return
  110 format(//' Error in the the indata file:',// ' An exctraction of tensors calculated using TDDFT is asked.', &
           /' The corresponding file does not contain TDDFT outputs !!')
end

!***********************************************************************

subroutine extract_energ(Energ_s,Eseuil,multi_run,nbbseuil,nbseuil,nenerg_s,nom_fich_extract,tddft)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  character(len=132) mot, nom_fich_extract

  logical tddft

  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(nenerg_s):: energ_s

  open(1, file = nom_fich_extract, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(nom_fich_extract,istat,1)

  i = 0
  do
    read(1,'(A)' ) mot
    if( mot(2:15) /= 'Absorbing atom' ) cycle
    i = i + 1
    if( i == multi_run ) exit
  end do

  if( tddft ) then
    do
      read(1,'(A)') mot
      if( mot(2:12) == 'Cycle TDDFT' .or. mot(2:12) == 'Cycle Tddft' ) exit
    end do
  endif

  boucle_ie: do ie = 1,nenerg_s

    do
      read(1,'(A)') mot
      if( mot(7:11) == 'Coabs' ) exit
    end do

    do
      read(1,'(A)') mot
      if( mot(6:10) == 'energ' .or. mot(5:10) == 'Energy' ) then
        read(1,*) energ_s(ie)
        energ_s(ie) = energ_s(ie) / rydb
        cycle boucle_ie
      endif
    end do

  end do boucle_ie

  if( energ_s(nenerg_s) > Eseuil(min(nbbseuil,nbseuil)) .and. Eseuil(nbbseuil) > energ_s(nenerg_s) - energ_s(1) ) then
    do ie = 1,nenerg_s
      energ_s(ie) = energ_s(ie) - Eseuil(min(nbbseuil,nbseuil))
    end do
  endif

  Close(1)

  return
end

!*********************************************************************

! Elaboration de l'agregat (celui qui sert au calcul du potentiel)
! Evaluation de la dimension des tableaux qu'on va utiliser dans le sousprogramme agregat

subroutine natomp_cal(angxyz,ATA,axyz,Base_ortho,chargat,d_ecrant,dcosxyz,deccent,Doping,dpos,Flapw,iabsorbeur,iabsfirst, &
      icheck,igr_dop,igreq,itabs,itype,Kgroup,Matper,mpirank,multrmax,n_atom_proto,n_radius,natomeq_s,natomeq_coh, &
      natomp,neqm,ngreq,ngroup,ngroup_pdb,ngroup_taux,noncentre,posn,Proto_all,r_self,rsorte_s,rmax,rpotmax, &
      self_cons,Taux,Taux_oc,Test_dist_min)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer:: i_radius, iabsfirst, igr_dop, itabs, n_atom_proto, mpirank, multrmax, n_radius, natomeq_coh, natomp, ngroup, &
    ngroup_taux

  integer, dimension(ngroup):: itype
  integer, dimension(n_radius):: natomeq_s
  integer, dimension(ngroup_pdb):: Kgroup
  integer, dimension(0:n_atom_proto):: ngreq
  integer, dimension(0:n_atom_proto,neqm):: igreq

  logical:: ATA, Base_ortho, Doping, flapw, matper, noncentre, Proto_all, self_cons, Taux
  logical, dimension(n_atom_proto):: ipr_ok

  real(kind=db):: r_self, rmax, rpotmax
  real(kind=db), dimension(n_radius):: rsorte_s
  real(kind=db), dimension(3):: angxyz, axyz, dcosxyz, deccent, dpos, ps, v
  real(kind=db), dimension(ngroup_taux):: Taux_oc
  real(kind=db), dimension(3,ngroup):: posg, posn
  real(kind=db), dimension(0:n_atom_proto):: chargat

  rad = pi / 180
  dcosxyz(:) = 2 * cos( angxyz(:) * rad )
  if( abs( dcosxyz(1) ) < eps10 .and. abs( dcosxyz(2) ) < eps10 .and. abs( dcosxyz(3) ) < eps10 ) then
    Base_ortho = .true.
  else
    Base_ortho = .false.
  endif

! L'atome absorbeur est mis au centre
! (en cas de One_run, c'est le premier absorbeur)
  if( noncentre ) then
    deccent(:) = 0._db
  else
    deccent(:) = posn(:,iabsfirst) * axyz(:)
  endif

  do igr = 1,ngroup
    posg(1:3,igr) = posn(1:3,igr) * axyz(1:3) - deccent(1:3)
  end do

  if( rpotmax > eps4 ) then
    rmax = rpotmax
  else
    rmax = max(rsorte_s(1),r_self) + multrmax * 2.5 / bohr
  endif

  if( matper ) then
    nxmaille = nint( rmax / axyz(1) ) + 3
    nymaille = nint( rmax / axyz(2) ) + 3
    nzmaille = nint( rmax / axyz(3) ) + 3
  else
    nxmaille = 0
    nymaille = 0
    nzmaille = 0
  endif

  natomeq_s(:) = 0
  natomeq_coh = 0
  natomp = 0
  ipr_ok(:) = .false.

  do ix = -nxmaille,nxmaille
    v(1) = ix * axyz(1)
    do iy = -nymaille,nymaille
      v(2) = iy * axyz(2)
      do iz = -nzmaille,nzmaille
        v(3) = iz * axyz(3)
        boucle_igr: do igr = 1,ngroup
          if( Doping .and. igr == ngroup .and. .not. ( ix == 0 .and. iy == 0 .and. iz == 0 ) ) cycle
          if( Doping .and. igr == igr_dop .and. ( ix == 0 .and. iy == 0 .and. iz == 0 ) ) cycle
          if( ngroup_pdb > 0 ) then
            if( Kgroup(iabsorbeur) /= Kgroup(igr) .and. Kgroup(iabsorbeur) /= 0 .and. Kgroup(igr) /= 0 ) cycle
          endif
          if( ( ix == 0 .and. iy == 0 .and. iz == 0 ) .and. itype(igr) /= itype(iabsorbeur) &
              .and. sum( abs(posg(:,igr) - posg(:,iabsorbeur)) ) < eps6 ) cycle
          ps(1:3) = posg(1:3,igr) + v(1:3)
          if( .not. ATA .and. .not. ( ix == 0 .and. iy == 0 .and. iz == 0 .and. igr == iabsorbeur ) ) then
            do jgr = 1,igr-1
              if( sum( abs( posg(:,igr) - posg(:,jgr) ) ) < eps6 )  cycle boucle_igr
            end do
          endif
          dist = vnorme(Base_ortho,dcosxyz,ps)
          if( dist <= r_self + eps10 ) then
            natomeq_coh = natomeq_coh + 1
            boucle_ipr: do ipr = 1,n_atom_proto
              do jgr = 1,ngreq(ipr)
                if( igreq(ipr,jgr) /= igr ) cycle
                ipr_ok(ipr) = .true.
                exit boucle_ipr
              end do
            end do boucle_ipr
          endif
          do i_radius = 1,n_radius
            if( dist <= rsorte_s(i_radius) + eps10 ) natomeq_s(i_radius) = natomeq_s(i_radius) + 1
          end do
          if( dist <= rmax + eps10 ) natomp = natomp + 1
        end do boucle_igr
      end do
    end do
  end do

  Proto_all = .true.
  do ipr = 1,n_atom_proto
    if( ipr_ok(ipr) ) cycle
    Proto_all = .false.
    exit
  end do

  if( .not. flapw .and. matper ) then
    rsortm = max( r_self, rsorte_s(1) )
    call reduc_natomp(ATA,axyz,Base_ortho,chargat,d_ecrant,dcosxyz,deccent,Doping,dpos,iabsorbeur,icheck,igr_dop,igreq,itabs, &
            itype,Kgroup,matper,mpirank,n_atom_proto,natomp,natomr,ngreq,ngroup,ngroup_pdb,ngroup_taux,noncentre,posn,rmax, &
            rsortm,Taux_oc)
    natomp = natomr
  endif

  if( mpirank == 0 ) then

    if( icheck > 0 ) write(3,110)

! Test sur les atomes trop proches
    istop = 0
    do igr1 = 1,ngroup
      if( Doping .and. igr1 == ngroup ) cycle
      do igr2 = igr1+1,ngroup
        if( Doping .and. igr2 == ngroup ) cycle
        if( ngroup_pdb > 0 ) then
          if( Kgroup(igr1) /= Kgroup(igr2) .and. Kgroup(igr1) /= 0 .and. Kgroup(igr2) /= 0 ) cycle
        endif
        ps(1:3) = posg(1:3,igr1) - posg(1:3,igr2)
        r = vnorme(Base_ortho,dcosxyz,ps)
        if( r > Test_dist_min ) cycle
        if( Taux .and. r < eps10 .and. .not. Doping ) then
          if( Taux_oc(igr1) + Taux_oc(igr2) < 1._db + eps10 ) cycle
        endif
        if( istop == 0 ) call write_error
        do iprint = 3,9,3
          if( istop == 0 ) write(iprint,120)
          write(iprint,130) igr1, posn(1:3,igr1), igr2, posn(1:3,igr2), r*bohr
        end do
        istop = 1
      end do
    end do
    if( istop == 1 ) then
      do iprint = 3,9,3
        write(iprint,140)
      end do
      stop
    endif

    do ipr = 3,6,3
      if( icheck == 0 .and. ipr == 3 ) cycle
      if( .not. self_cons .or. abs(rsorte_s(1) - r_self) < eps10 ) then
        write(ipr,150) ( rsorte_s(i_radius)*bohr, natomeq_s(i_radius), i_radius = 1,n_radius )
      else
        do i_radius = 1,n_radius
          if( i_radius == 1 ) then
            write(ipr,160) rsorte_s(i_radius)*bohr, natomeq_s(i_radius)
          else
            write(ipr,165) rsorte_s(i_radius)*bohr, natomeq_s(i_radius)
          endif
        end do
        if( self_cons ) write(ipr,170) r_self*bohr, natomeq_coh
      endif
      if( ipr == 3 ) write(ipr,180) rmax*bohr, natomp
    end do
  endif

  return
  110 format(/' ---- Natomp_cal ----',100('-'),/)
  120 format(//' Error in the indata file, the following atoms are too close:',/)
  130 format(' Atoms',i5,' at p =',3f11.7,/'   and',i5,' at p =',3f11.7, ', distance =',f11.7,' A')
  140 format(/' This can come from:'/, &
              '  - Two atoms set at the same position,'/, &
              '  - Two atoms too close,'/, &
              '  - When using Spgroup keyword, the position of an atom standing on',/ &
              '    a symmetry plane or on a rotation axis must follow precisely the convenient',/ &
              '    ratio between x, y and z of the site (for example (h,2h,0) ),'/, &
              '  - An unsufficient number of digits for the atom position when it is 1/3, 2/3...'/ &
              '    One must write them with at least 10 digits,'/, &
              '    for example 0.3333333333 and not 0.3333 !'//)
  150 format(' Cluster radius =',f5.2,' A, nb. of atom =',i4)
  160 format(' Absorption calculation   : cluster radius =',f5.2,' A, nb. of atom =',i4)
  165 format('                            cluster radius =',f5.2,' A, nb. of atom =',i4)
  170 format(' Fermi energy calculation : cluster radius =',f5.2,' A, nb. of atom =',i4)
  180 format(' Potential sup calculation: cluster radius =',f5.2,' A, nb. of atom =',i4)
   end

!*********************************************************************

subroutine reduc_natomp(ATA,axyz,Base_ortho,chargat,d_ecrant,dcosxyz,deccent,Doping,dpos,iabsorbeur,icheck,igr_dop,igreq,itabs, &
            itype,Kgroup,matper,mpirank,n_atom_proto,natomp,natomr,ngreq,ngroup,ngroup_pdb,ngroup_taux,noncentre,posn,rmax, &
            rsortm,Taux_oc)

  use declarations
  implicit none

  integer:: ia, iaabs, iaabsfirst, iabsorbeur, icheck, igr, igr_dop, ipr, itabs, mpirank, n_atom_proto, na, natomp, &
            natomq, natomr, ngroup, ngroup_pdb, ngroup_taux
  integer, dimension(natomp):: iaproto, igroup, itypep
  integer, dimension(ngroup_pdb):: Kgroup
  integer, dimension(ngroup):: itype
  integer, dimension(0:n_atom_proto):: ngreq
  integer, dimension(0:n_atom_proto,ngroup):: igreq

  logical:: ATA, Base_ortho, Doping, Matper, noncentre, OK

  real(kind=db):: ch, ch_min, ch_test, chagreg, chg, d_ecrant, rmax, rsortm
  real(kind=db), dimension(3):: axyz, dcosxyz, deccent, dpos
  real(kind=db), dimension(natomp):: dista
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(3,ngroup):: posn
  real(kind=db), dimension(ngroup_taux):: Taux_oc
  real(kind=db), dimension(0:n_atom_proto):: chargat

  matper = .true.
  call clust(ATA,axyz,Base_ortho,dcosxyz,deccent,dista,Doping,dpos,iaabs,iaabsfirst,iabsorbeur,igr_dop,igroup,itabs,itype, &
              itypep,Kgroup,matper,mpirank,1,natomp,ngroup,ngroup_pdb,ngroup_taux,noncentre,.false.,pos,posn,rmax,Taux_oc)

  do ia = 1,natomp
    boucle_i: do ipr = 1,n_atom_proto
      do igr = 1,ngreq(ipr)
        if( igroup(ia) == igreq(ipr,igr) ) exit boucle_i
      end do
    end do boucle_i
    iaproto(ia) = ipr
  end do

! Calcul de la charge de l'agregat :
  chagreg = 0._db
  do ia = 1,natomp
    ipr = iaproto(ia)
    chagreg = chagreg + chargat( ipr )
  end do

  if( abs( chagreg ) < eps10 ) then
    natomr = natomp
    return
  endif

  do ia = 1,natomp-1
    if( dista(ia+1) > rsortm ) exit
  end do
  natomq = ia

  na = natomq

  OK = .false.

  ch = chagreg
  do ia = natomp-1,natomq,-1
    ipr = iaproto(ia+1)
    ch = ch - chargat(ipr)
    if( abs( dista(ia+1) - dista(ia) ) > eps10 ) then
      if( abs( ch - d_ecrant ) > eps10  ) cycle
      na = ia
      chg = ch
      OK = .true.
      exit
    endif
  end do

  if( .not. OK ) then
    na = natomp
    ch_min = chagreg - d_ecrant
    chg = chagreg
    do ia = natomp-1,natomq,-1
      ipr = iaproto(ia+1)
      chg = chg - chargat(ipr)
      if( abs( dista(ia+1) - dista(ia) ) < eps10 ) cycle
      ch_test = chg - d_ecrant
      if( abs(ch_test) > abs(ch_min) - eps10 ) cycle
      ch_min = ch_test
      na = ia
    end do
  endif

  if( na > natomq .and. na < natomp ) then
    natomr = na
    chagreg = chg
    if( icheck > 0 ) write(3,110) na, dista(na)*bohr
    rmax = dista(natomr) + eps10
  else
    natomr = natomp
  endif

  return
  110 format(/' natomp reduced at =',i4,', dista(natomp) =',f11.5)
end

!*********************************************************************

subroutine clust(ATA,axyz,Base_ortho,dcosxyz,deccent,dista,Doping,dpos,iaabs,iaabsfirst,iabsorbeur,igr_dop,igroup,itab,itype, &
              itypep,Kgroup,matper,mpirank,multi_run,natomp,ngroup,ngroup_pdb,ngroup_taux,noncentre,One_run,pos,posn,rmax,Taux_oc)

  use declarations
  implicit none

  integer:: ia, ia1, ia2, iaabs, iaabsfirst, iabsorbeur, ib, igr, igr_dop, igr12, ipr, itab, ity12, ix, iy, iz, mpirank, &
    multi_run, natomp, ngroup, ngroup_pdb, ngroup_taux, nxmaille, nymaille, nzmaille
  integer, dimension(natomp):: igroup, itypep
  integer, dimension(ngroup):: itype
  integer, dimension(ngroup_pdb):: Kgroup

  logical:: ATA, Base_ortho, Doping, Matper, Noncentre, One_run

  real(kind=db):: dist, dist12, rmax, vnorme
  real(kind=db), dimension(3):: axyz, dcosxyz, deccent, dpos, ps, v
  real(kind=db), dimension(ngroup_taux):: Taux_oc
  real(kind=db), dimension(natomp):: dista
  real(kind=db), dimension(3,ngroup):: posg, posn
  real(kind=db), dimension(3,natomp):: pos

  do igr = 1,ngroup
    posg(1:3,igr) = posn(1:3,igr) * axyz(1:3) - deccent(1:3)
  end do

  if( matper ) then
    nxmaille = nint( rmax / axyz(1) ) + 3
    nymaille = nint( rmax / axyz(2) ) + 3
    nzmaille = nint( rmax / axyz(3) ) + 3
  else
    nxmaille = 0
    nymaille = 0
    nzmaille = 0
  endif

  iaabs = 0
  ia = 0

  do ix = -nxmaille,nxmaille
    v(1) = ix * axyz(1)
    do iy = -nymaille,nymaille
      v(2) = iy * axyz(2)
      do iz = -nzmaille,nzmaille
        v(3) = iz * axyz(3)
        boucle_igr: do igr = 1,ngroup
          if( Doping .and. igr == ngroup .and. .not. ( ix == 0 .and. iy == 0 .and. iz == 0 ) ) cycle
          if( Doping .and. igr == igr_dop .and. ( ix == 0 .and. iy == 0 .and. iz == 0 ) ) cycle
          if( ngroup_pdb > 0 ) then
            if( Kgroup(iabsorbeur) /= Kgroup(igr) .and. Kgroup(iabsorbeur) /= 0 .and. Kgroup(igr) /= 0 ) cycle
          endif
          if( ( ix == 0 .and. iy == 0 .and. iz == 0 ) .and. itype(igr) /= itype(iabsorbeur) &
              .and. sum( abs(posg(:,igr) - posg(:,iabsorbeur)) ) < eps6 ) cycle
          ps(1:3) = posg(1:3,igr) + v(1:3)
          if( .not. ATA .and. .not. ( ix == 0 .and. iy == 0 .and. iz == 0 .and. igr == iabsorbeur ) ) then
            do ib = 1,ia
              if( sum( abs( pos(:,ib) - ps(:) ) ) < eps6 ) then
                if( Taux_oc( igroup(ib) ) <= Taux_oc( igr ) - eps10 ) itypep(ib) = abs( itype(igr) )
                cycle boucle_igr
              endif
            end do
          endif
          dist = vnorme(Base_ortho,dcosxyz,ps)
          if( dist > rmax + eps10 ) cycle
          ia = ia + 1
          pos(1:3,ia) = ps(1:3)
          igroup(ia) = igr
          if( igr == iabsorbeur .and. ix == 0 .and. iy == 0 .and. iz == 0 ) then
            itypep(ia) = itab
            iaabs = ia
          else
            itypep(ia) = abs( itype(igr) )
          endif
        end do boucle_igr
      end do
    end do
  end do

  if( iaabs == 0 .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,110)
    end do
    stop
  endif

  if( .not. One_run .or. ( One_run .and. multi_run == 1 ) ) iaabsfirst = iaabs

  if( noncentre ) then
    pos(:,iaabsfirst) = pos(:,iaabsfirst) + dpos(:)
  else
    do ia = 1,natomp
      if( ia /= iaabsfirst ) pos(:,ia) = pos(:,ia) - dpos(:)
    end do
  endif

! Mise en ordre par distance croissante au centre
  do ia = 1,natomp
    v(:) = pos(:,ia)
    dista(ia) = vnorme(Base_ortho,dcosxyz,v)
  end do

  do ia1 = 1,natomp
    do ia2 = ia1+1,natomp
      if( dista(ia1) < dista(ia2)+eps10 ) cycle
      dist12 = dista(ia1)
      dista(ia1) = dista(ia2)
      dista(ia2) = dist12
      igr12 = igroup(ia1)
      igroup(ia1) = igroup(ia2)
      igroup(ia2) = igr12
      ity12 = itypep(ia1)
      itypep(ia1) = itypep(ia2)
      itypep(ia2) = ity12
      v(:) = pos(:,ia1)
      pos(:,ia1) = pos(:,ia2)
      pos(:,ia2) = v(:)
      if( ia1 == iaabs ) then
        iaabs = ia2
      elseif( ia2 == iaabs ) then
        iaabs = ia1
      endif
    end do
  end do

  return
  110 format(/' There is no absorbing atom in the calculation sphere !')
end

!*********************************************************************

! Elaboration de l'agregat.

subroutine agregat(angxyz,ATA,Atom_with_axe,Atom_nonsph,Axe_atom_clu,Axe_atom_gr,Axe_atom_grn,axyz,Base_hexa,Base_ortho, &
          chargat,chargm,Cubmat,dcosxyz,deccent,dista,Doping,dpos,Flapw,iaabs,iaabsfirst,iabsorbeur,iaproto,iapot, &
          icheck,igr_dop,igreq,igroup,igrpt_nomag,igrpt0,iopsymc,iopsymr,itab,itype,itypep,karact,Kgroup,Magnetic,Matper, &
          mpirank,multi_run,n_atom_proto,natomp,nb_rep,nb_sym_op,neqm,ngreq,ngroup,ngroup_m,ngroup_pdb,ngroup_taux,nlat, &
          nlatm,noncentre,nonexc_g,nspin,ntype,numat,One_run,Orthmat,Orthmati,PointGroup,PointGroup_Auto, &
          popats,pos,posn,rmax,Rot_int,Self_nonexc,Spinorbite,Rot_Atom_gr,sym_4,Struct,Sym_cubic,Symmol,Taux,Taux_oc,Vsphere)

  use declarations
  implicit none

  integer:: i, ia, iaabs, iaabsfirst, iabsorbeur, ib, icheck, igr, igr_dop, igrpt, igrpt_nomag, igrpt_sg, igrpt_sg_cmp, &
    igrpt_sg_cal, igrpt_sg_so, igrpt0, igrptn, ipr, ired, it, itab, k, mpirank, multi_run, n_atom_proto, natomp, nb_rep, &
    nb_sym_op, neqm, ngroup, ngroup_m, ngroup_pdb, ngroup_taux, nlatm, npr1, nspin, ntype

  integer, dimension(natomp):: iaproto, igroup, itypep
  integer, dimension(nopsm):: iopsymc, iopsymr
  integer, dimension(0:ntype):: nlat, numat
  integer, dimension(0:n_atom_proto):: iapot, ngreq
  integer, dimension(ngroup):: itype
  integer, dimension(ngroup_pdb):: Kgroup
  integer, dimension(0:n_atom_proto,neqm):: igreq

  character(len=8):: PointGroup, PointGroup_name, PointGroup_nomag_name, PointGroup_Sch, &
                     PointSubGroup_name, PointSubGroup_Sch
  character(len=5):: Struct

  logical:: ATA, Atom_mag_cal, Atom_nonsph, Base_hexa, Base_ortho, Doping, Flapw, Magnetic, Matper, Noncentre, Nonexc_g, &
    One_run, PointGroup_Auto, Self_nonexc, Spinorbite, Sym_4, Sym_cubic, Symmol, Taux
  logical, dimension(0:ngroup_m):: Atom_with_axe

  real(kind=db):: Angle, Chagreg, Chargm, Dist_max, Distm, Dsa, Dsb, Prod, rad, Rmax, Vsphere

  real(kind=db), dimension(3):: angxyz, Axe_spin, Axe_spin_0, axyz, dcosxyz, deccent, dpos, pa, pb, v
  real(kind=db), dimension(3,3):: Cubmat, Orthmat, Orthmati, Orthmatt, Rot_un, Rotmat, Rot_int, Rot_so, Rot_tem, Rotspin
  real(kind=db), dimension(natomp):: dista
  real(kind=db), dimension(3,natomp):: Axe_atom_clu, pos
  real(kind=db), dimension(3,ngroup):: posn
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_gr, Axe_atom_grn
  real(kind=db), dimension(ngroup,nlatm,nspin) :: popats
  real(kind=db), dimension(0:n_atom_proto):: chargat
  real(kind=db), dimension(3,3,ngroup_m):: Rot_Atom_gr
  real(kind=db), dimension(ngroup_taux):: taux_oc

  complex(kind=db), dimension(nopsm,nrepm):: karact

  common/PointGroup_name/ PointGroup_name, PointGroup_nomag_name
  common/PointGroup_Sch/ PointGroup_Sch
  common/PointSubGroup_name/ PointSubGroup_name, PointSubGroup_Sch

  if( icheck > 0 ) write(3,110)

  rad = pi / 180
  dcosxyz(:) = 2 * cos( angxyz(:) * rad )
  if( abs( dcosxyz(1) ) < eps10 .and. abs( dcosxyz(2) ) < eps10 .and. abs( dcosxyz(3) ) < eps10 ) then
    base_ortho = .true.
  else
    base_ortho = .false.
  endif

  call clust(ATA,axyz,Base_ortho,dcosxyz,deccent,dista,Doping,dpos,iaabs,iaabsfirst,iabsorbeur,igr_dop,igroup,itab,itype, &
              itypep,Kgroup,matper,mpirank,multi_run,natomp,ngroup,ngroup_pdb,ngroup_taux,noncentre,One_run,pos,posn,rmax,Taux_oc)

  if( icheck > 0 ) write(3,130) iaabs

  call cal_iaeqrmt(iaabs,iaproto,iapot,igreq,igroup,n_atom_proto,natomp,neqm,ngreq,nonexc_g,self_nonexc)

  if( magnetic .or. Atom_nonsph ) then
    do igr = 1,ngroup
      v(:) = Axe_atom_gr(:,igr) * axyz(:)
      Axe_atom_grn(:,igr) = v(:)
      do ia = 1,natomp
        if( igroup(ia) == igr ) Axe_atom_clu(:,ia) = v(:)
      end do
    end do
  else
    Axe_atom_clu(:,:) = 0._db
  endif

! On passe en repere orthogonal
  if( Struct /= 'cubic' ) then
    do ia = 1,natomp
      v(:) = pos(:,ia)
      do k = 1,3
        pos(k,ia) = sum( cubmat(k,1:3) * v(1:3) )
      end do
    end do
    if( Magnetic .or. Atom_nonsph ) then
      do ia = 1,ngroup
        v(:) = Axe_atom_gr(:,ia) * axyz(:)
        do k = 1,3
          Axe_atom_grn(k,ia) = sum( cubmat(k,1:3) * v(1:3) )
        end do
      end do
      do ia = 1,natomp
        v(:) = Axe_atom_clu(:,ia)
        do k = 1,3
          Axe_atom_clu(k,ia) = sum( cubmat(k,1:3) * v(1:3) )
        end do
      end do
    endif
    dcosxyz(:) = 0._db
    Base_ortho = .true.
  endif

  Orthmatt(:,:) = 0._db
  do k = 1,3
    orthmatt(k,k) = axyz(k)
  end do
  Orthmatt = matmul( cubmat, orthmatt )

  Orthmat = Orthmatt

  rot_un(:,:) = 0._db
  do i = 1,3
    rot_un(i,i) = 1._db
  end do
  rotspin(:,:) = rot_un(:,:)
  rotmat(:,:) = rot_un(:,:)

! On prend comme axe z, l'axe du spin de l'absorbeur
! Faux en cas de one_run avec des spins non alignes.
  Axe_spin(:) = 0._db
  Atom_with_axe(0) = .false.
  if( Magnetic .or. Atom_nonsph ) then
    do ia = 1,natomp
      if( ia /= iaabs ) cycle
      igr = igroup(ia)
      if( .not. Atom_with_axe(igr) ) cycle
      Atom_with_axe(0) = Atom_with_axe(igr)
      Axe_spin(:) = Axe_atom_clu(:,ia)
      rotspin(:,:) = Rot_Atom_gr(:,:,igr)
      rotspin = transpose( rotspin )
      exit
    end do
  endif

! Rotation de la base interne en cas d'axe de spin n'etant pas selon Oz
  if( Atom_with_axe(0) ) then

    orthmat = matmul( rotspin, orthmat )

    do ia = 1,natomp
      v(:) = pos(:,ia)
      v = matmul( rotspin, v )
      pos(:,ia) = v(:)
    end do

    do ia = 1,natomp
      v(:) = Axe_atom_clu(:,ia)
      v = matmul( rotspin, v )
      Axe_atom_clu(:,ia) = v(:)
    end do

    Axe_spin = matmul( rotspin, Axe_spin )

  endif

! Determination des symetries du groupe ponctuel (iopsymc)

  if( PointGroup_Auto ) then

    rot_so = rot_un
    Axe_spin_0 = Axe_spin

    do i = 1,2

      call sym_cluster(Atom_with_axe,Axe_atom_clu,Base_ortho,dcosxyz,iaabs,igroup,iopsymc,itype,itypep, &
             natomp,ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,pos,Symmol)

      call cluster_rot(iopsymc,rot_tem)

      if( sum( abs( rot_tem(:,:) - rot_un(:,:) ) ) < eps10  ) exit

      orthmat = matmul( rot_tem, orthmat )
      rotmat = matmul( rot_tem, rotmat )

      do ia = 1,natomp
        v(:) = pos(:,ia)
        v = matmul( rot_tem, v )
        pos(:,ia) = v(:)
      end do
      if( magnetic .or. Atom_nonsph ) then
        do ia = 1,natomp
          v(:) = Axe_atom_clu(:,ia)
          v = matmul( rot_tem, v )
          Axe_atom_clu(:,ia) = v(:)
        end do
        Axe_spin = matmul( rot_tem, Axe_spin )
      endif
      if( sum( abs( axe_spin(:) - axe_spin_0(:) ) ) > eps10 ) rot_so = matmul( rot_so, transpose(rot_tem) )

    end do

    call sym_cluster(Atom_with_axe,Axe_atom_clu,Base_ortho,dcosxyz,iaabs,igroup,iopsymc,itype,itypep, &
             natomp,ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,pos,Symmol)

    call numgrpt(iopsymc,igrpt,igrpt_nomag,mpirank,.false.)

!        if( spinorbite) then
      orthmat = matmul( rot_so, orthmat )
      rotmat = matmul( rot_so, rotmat )

      do ia = 1,natomp
        v(:) = pos(:,ia)
        v = matmul( rot_so, v )
        pos(:,ia) = v(:)
      end do
      if( magnetic .or. Atom_nonsph ) then
        do ia = 1,natomp
          v(:) = Axe_atom_clu(:,ia)
          v = matmul( rot_so, v )
          Axe_atom_clu(:,ia) = v(:)
        end do
        Axe_spin = matmul( rot_so, Axe_spin )
      endif
      call sym_cluster(Atom_with_axe,Axe_atom_clu,Base_ortho,dcosxyz,iaabs,igroup,iopsymc,itype,itypep, &
             natomp,ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,pos,Symmol)
!        endif

  else

    call grp_opsym_imp(iopsymc,igrpt,igrpt_nomag,mpirank,PointGroup)

  endif

! iopsymc est la symetrie vraie.
! iopsymr est la symetrie correspondant utilisee pour le calcul de la stucture
! electronique. Il correspond au sous-groupe pour lequel les
! representations a calculer sont unidimensionnelles.

  if( icheck > 0 ) then
    npr1 = 3
  else
    npr1 = 6
  endif

  if( mpirank == 0 ) then
    do ipr = npr1,6,3
      if( igrpt > 32 ) then
        write(ipr,140) PointGroup_name, igrpt
        write(ipr,150) PointGroup_nomag_name, PointGroup_Sch
        if( PointGroup_Auto ) write(ipr,160) PointSubgroup_name, PointSubgroup_Sch
      else
        write(ipr,170) PointGroup_name, PointGroup_Sch
      endif
    end do
    if( icheck > 0 ) then
      write(3,175)
      call write_iopsym(iopsymc,3)
    endif
  endif

! On prend un sous-groupe.
  igrptn = igrpt
  if( PointGroup_Auto ) then
    call numgrpt(iopsymc,igrpt,igrpt_nomag,mpirank,.true.)
    igrpt_sg = igrpt_sg_cal(.false.,igrpt,igrpt_nomag, Spinorbite)
  else
    igrpt_sg = igrpt
  endif

  sym_cubic = .false.
  sym_4 = .false.
  if( igrpt_sg /= 1 .and. .not. spinorbite ) then
    if( igrpt_nomag > 27 ) then
      sym_cubic = .true.
    elseif( igrpt_nomag >= 12 .and. igrpt_nomag <= 15 ) then
      sym_4 = .true.
    endif
  endif

  if( igrpt_sg /= igrpt_nomag .or. igrptn /= igrpt ) then
    PointGroup = 'SousGr'
    igrpt_nomag = igrpt_sg
    iopsymr(:) = iopsymc(:)
    call grp_opsym_imp(iopsymr,igrpt,igrpt_nomag,mpirank,PointGroup)
    ired = 1
  else
    iopsymr(:) = iopsymc(:)
    ired = 0
  endif
  if( Spinorbite .or. igrpt_sg /= 8 ) then
    igrpt_sg = igrpt_sg_so(igrpt,igrpt_nomag)
    if( igrpt_sg /= igrpt_nomag ) then
      PointGroup = 'SousGr'
      igrpt_nomag = igrpt_sg
      iopsymr(:) = iopsymc(:)
      call grp_opsym_imp(iopsymr,igrpt,igrpt_nomag,mpirank,PointGroup)
      ired = 1
    endif
    igrpt_sg = igrpt_sg_cmp(igrpt,igrpt_nomag)
    if( igrpt_sg /= igrpt_nomag ) then
      PointGroup = 'SousGr'
      igrpt_nomag = igrpt_sg
      iopsymr(:) = iopsymc(:)
      call grp_opsym_imp(iopsymr,igrpt,igrpt_nomag,mpirank,PointGroup)
      ired = 1
    end if
  endif

  Atom_with_axe(0) = Atom_mag_cal(igrpt)

! igrpt0 est transporte dans la routine point_group_atom
  igrpt0 = igrpt

  if( ired == 1 .and. mpirank == 0 ) then
    do ipr = npr1,6,3
      if( igrpt > 32 ) then
        write(ipr,180) PointGroup_name, igrpt
        write(ipr,150) PointGroup_nomag_name, PointGroup_Sch
        if( PointGroup_Auto ) write(ipr,160) PointSubgroup_name, PointSubgroup_Sch
      else
        write(ipr,182) Pointgroup_name, Pointgroup_Sch
      endif
    end do
    if( icheck > 0 ) then
      write(3,183)
      call write_iopsym(iopsymr,3)
    endif
  endif

  call character_table(icheck,igrpt_nomag,karact,nb_rep)

  call invermat(orthmat,orthmati)

  rotmat = matmul( rotmat, rotspin )
  call invermat(rotmat,Rot_int)

! Nombre d'operations de symetrie du groupe utilise.
  nb_sym_op = sum( abs( iopsymr(:) ) )

  if( icheck > 0 ) then
    write(3,185)
    write(3,187)
    write(3,200) ( orthmatt(i,1:3), i = 1,3 )
    write(3,190)
    write(3,200) ( cubmat(i,1:3), i = 1,3 )
    write(3,202)
    write(3,200) ( rotmat(i,1:3), i = 1,3 )
    write(3,203)
    write(3,200) ( Rot_int(i,1:3), i = 1,3 )
    write(3,206)
    write(3,200) ( orthmat(i,1:3), i = 1,3 )
    write(3,207)
    write(3,200) ( orthmati(i,1:3), i = 1,3 )
    if( Taux ) then
      write(3,210)
    else
      write(3,215)
    endif
    do ia = 1,natomp
      ipr = iaproto(ia)
      it = itypep(ia)
      if( Taux ) then
        write(3,220) numat(it), pos(:,ia)*bohr, dista(ia)*bohr, Taux_oc( igroup(ia) ), ia, igroup(ia), it, ipr, chargat( ipr )
      else
        write(3,225) numat(it), pos(:,ia)*bohr, dista(ia)*bohr, ia, igroup(ia), it, ipr, chargat( ipr )
      endif
    end do
    write(3,230)
    do ipr = 0,n_atom_proto
      write(3,240) ipr, iapot(ipr)
    end do

! Angle des liasons autour de l'absorbeur
    write(3,'(/A)') ' ia ib Za Zb   Angle(a,O,b)'
    rad = 180 / pi
    Distm = 3.2_db / bohr
    Dist_max = Dista(iaabs) + Distm
    do ia = 1,natomp
      if( Dista(ia) > Dist_max ) exit
      if( ia == iaabs ) cycle
      pa(:) = pos(:,ia) - pos(:,iaabs)
      Dsa = sqrt( sum( pa(:)**2 ) )
      if( Dsa > Distm ) cycle
      do ib = ia + 1, natomp
        if( dista(ib) > Dist_max ) exit
        if( ib == iaabs ) cycle
        pb(:) = pos(:,ib) - pos(:,iaabs)
        Dsb = sqrt( sum( pb(:)**2 ) )
        if( Dsb > Distm ) cycle
        Prod = sum( pa(:) * pb(:) ) / ( Dsa * Dsb )
        if( Prod < -1._db + eps10 ) then
          Angle = 180._db
        elseif( Prod > 1._db - eps10 ) then
          Angle = 0._db
        else
          Angle = rad * acos( Prod )
        endif
        write(3,250) ia, ib, numat(itypep(ia)), numat(itypep(ib)), Angle
      end do
    end do
  endif

! En cas d'axe 3, maillage hexagonal.
  if( abs(iopsymc(49)) == 1 .or. abs(iopsymc(53)) == 1 ) then
    Base_hexa = .true.
  else
    Base_hexa = .false.
  endif

! Calcul de la charge de l'agregat :
  if( .not. flapw .and. matper ) then
    chagreg = 0._db
    do ia = 1,natomp
      ipr = iaproto(ia)
      chagreg = chagreg + chargat( ipr )
    end do
! On met une charge opposee sur la sphere a dista(natomp) du centre
    if( natomp > 1 ) then
      Vsphere = 2 * ( Chagreg - Chargm ) / dista(natomp)
    else
      Vsphere = 0._db
    endif
    if( icheck > 0 ) write(3,260) Chagreg, Vsphere * rydb
  else
    Vsphere = 0._db
  endif

  return
  110 format(/' ---- Agregat ------',100('-'))
  130 format(/' Index of the absorbing atom, iaabs =',i3)
  140 format(/' Point group : ',a8,'     number =',i3)
  150 format( '   Non magnetic corresponding point group   : ',a8,' (',a5,')')
  160 format( '   Subgroup not multiplied by time reversal : ',a8,' (',a5,')')
  170 format(/' Point group : ',a8,' (',a5,')')
  175 format(/' iopsymc =')
  180 format(/' Point group used : ',a8,'     number =',i3)
  182 format(/' Point group used : ',a8,' (',a5,')')
  183 format(/' iopsymr =')
  185 format(/' Transformation matrices :',/ '    R1 : internal orthonormal bases with z along c, x along b x c',/ &
  '         but for trigonal symmetry where z is along the hexagonal axis,', &
  /'         used for the tensorial expansion.',/ '    R2 : internal orthonormal bases used for the electronic', &
  ' structure calculation')
  187 format(/' Transformation Crystal Bases - Bases R1 :')
  190 format(/' Transformation Crystal Normalized Bases - Bases R1 :')
  200 format(3x,3f11.5)
  202 format(/' Rotation Bases R1 - Bases R2 :')
  203 format(/' Inverse matrix :')
  206 format(/' Transformation Crystal Bases - Bases R2 :')
  207 format(/' Inverse matrix :')
  210 format(/' Cluster: atom positions in order, in the internal R2 bases',/ &
   '   Z      posx       posy       posz          dista     Taux_oc  ia   igr ity ipr  chargat')
  215 format(/' Cluster: atom positions in order, in the internal R2 bases',/ &
   '   Z      posx       posy       posz          dista    ia   igr ity ipr  chargat')
  220 format(i4,3f11.6,'   ! ',f11.6,f10.5,i4,i6,2i4,f9.5)
  225 format(i4,3f11.6,'   ! ',f11.6,i4,i6,2i4,f9.5)
  230 format(/' ipr   iapot')
  240 format(i4,i6)
  250 format(4i3,f13.3)
  260 format(/' Cluster charge =',f11.5,'  Vsphere =',f11.5,' eV')
end

! **********************************************************************

! Determination du sous groupe ayant toutes les representations utiles
! de dimension 1.

! En cas de green on prend dans certains cas une symetrie plus basse
! a cause d'une erreur non trouvee dans la symetrisation (peut-etre
! calcul de Cmat dans la routine mat).
! (En cas de hubbard et fdm on prend un sous groupe a representations
! reelles) --> supprime

  integer function igrpt_sg_cal(Hubbard,igrpt,igrpt_nomag, Spinorbite)

  implicit none

  integer:: igrpt, igrpt_nomag, igrpt_sg

  logical:: Hubbard, Spinorbite

  igrpt_sg = igrpt_nomag

!      if( hubbard .and. .not. green ) then
  if( hubbard ) then

    select case(igrpt_nomag)
      case(9,10,22,25,26)
        igrpt_sg = 4
      case(11,23,27)
        igrpt_sg = 5
      case(12)
        igrpt_sg = 6
      case(13,14,28,30,31)
        igrpt_sg = 7   ! D2
      case(15,29)
        if( igrpt /= igrpt_nomag ) then
          igrpt_sg = 5
        else
          igrpt_sg = 8
        endif
      case(16,18,19)
        igrpt_sg = 1
      case(17,20)
        igrpt_sg = 2
      case(21,24)
        igrpt_sg = 3
      case(32)
        igrpt_sg = 8   ! D2h

    end select

  else

    select case(igrpt_nomag)
      case(12)
        igrpt_sg = 9    ! C4 a la place de C2v (mm2)
      case(13)
        igrpt_sg = 10   ! S4 a la place de D2 (222)
      case(14)
        igrpt_sg = 9    ! C4 a la place de D2 (222)
      case(28)
        igrpt_sg = 7
      case(15,29)
        if( igrpt /= igrpt_nomag ) then
          igrpt_sg = 11
        else
          igrpt_sg = 8
        endif
      case(18,19)
        igrpt_sg = 16
      case(20)
        igrpt_sg = 17
      case(24)
        igrpt_sg = 21
      case(25,26)
        igrpt_sg = 22
      case(27)
        igrpt_sg = 23
      case(30)
!            if( l_val_max > 2 .or. green ) then
          igrpt_sg = 10  ! S4 a la place de D2 (222)
!            else
!              igrpt_sg = 13  ! D2d
!            endif
      case(31)
!            if( l_val_max > 2 .or. green ) then
          igrpt_sg = 9   ! C4 a la place de D2 (222)
!            else
!              igrpt_sg = 14  ! D4
!            endif
      case(32)
!            if( l_val_max > 2 .or. green .or. Spinorbite ) then
        if( Spinorbite ) then
          igrpt_sg = 11   ! C4h a la place de D2h
        else
          igrpt_sg = 8  ! D2h = mmm
        endif
    end select

  endif

  igrpt_sg_cal = igrpt_sg

  return
end

!*********************************************************************

subroutine write_iopsym(iopsym,ipr)

  use declarations

  character(len=9):: mot9, nomsym
  character(len=42):: mot

  integer, dimension(nopsm):: iopsym

  i1 = 1
  do ligne = 1,20
    select case(ligne)
      case(1,9)
        ni = 1
      case(4,5,6,7,8,10,11,14,15,16)
        ni = 3
      case(2,3,12,13,17,18,19,20)
        ni = 4
    end select
    i2 = i1 + ni - 1
    mot = ' '
    do i = i1,i2
      l = len_trim(mot)
      l = l + 1
      mot9 = adjustl( nomsym(i) )
      ln = len_trim( mot9 )
      mot(l+1:l+ln) = mot9(1:ln)
      if( i < i2 ) mot(l+ln+1:l+ln+1) = ','
    end do
    write(ipr,110) adjustr(mot),iopsym(i1:i2)
    i1 = i2 + 1
  end do

  return
  110 format(a,' : ',8i2)
end

!*********************************************************************

  real(kind=db) function vnorme(Base_ortho,dcosxyz,v)

  use declarations
  implicit none

  logical Base_ortho

  real(kind=db), dimension(3):: dcosxyz, v

  if( Base_ortho ) then
    vnorme = sqrt( sum( v(:)**2 ) )
  else
    vnorme = sqrt( sum( v(:)**2 ) + v(1)*v(2)*dcosxyz(3) + v(1)*v(3)*dcosxyz(2) + v(2)*v(3)*dcosxyz(1) )
  endif

  return
end

!*********************************************************************

! Sous-programme calculant la rotation a effectuer pour se ramener a un groupe ponctuel standard.

subroutine cluster_rot(iopsym,rotmat)

  use declarations
  implicit none

  integer:: i, i2, is, n_opsym, n_opsym_2, n_opsym_3, n_opsym_s
  integer, dimension(nopsm):: iops, iopsym

  real(kind=db):: cs, r2, sn, vn, wn
  real(kind=db), dimension(3):: k, v, w, x
  real(kind=db), dimension(3,3):: rotmat

  iops(:) = abs( iopsym(:) )

  rotmat(:,:) = 0._db
  rotmat(1,1) = 1._db; rotmat(2,2) = 1._db; rotmat(3,3) = 1._db;

  n_opsym = sum( iops(:) )

  n_opsym_3 = sum( iops(2:9) )
  if( n_opsym_3 == 8 ) return

  r2 = 1 / sqrt( 2._db )
  v(:) = 0._db
  k(1) = 0._db; k(2) = 0._db; k(3) = 1._db

  n_opsym_2 = sum( iops(10:15) ) + sum( iops(22:24) ) + sum( iops(58:64:2) )
  if( n_opsym_2 == 1 ) then
    do i = 10,15
      if( iops(i) == 1 ) i2 = i
    end do
    do i = 22,24
      if( iops(i) == 1 ) i2 = i
    end do
    do i = 58,64,2
      if( iops(i) == 1 ) i2 = i
    end do
  else
    i2 = 0
  endif

  n_opsym_s = sum( iops(40:48) ) + sum( iops(57:63:2) )
  if( n_opsym_s == 1 ) then
    do i = 40,48
      if( iops(i) == 1 ) is = i
    end do
    do i = 57,63,2
      if( iops(i) == 1 ) is = i
    end do
  else
    is = 0
  endif

  if( n_opsym_2 == 1 .or. n_opsym_s == 1 ) then

    if( i2 == 10 .or. is == 48 ) then
      v(1) = 1._db; v(2) = 1._db; v(3) = 0._db
    elseif( i2 == 11 .or. is == 45 ) then
      v(1) = -1._db; v(2) = 1._db; v(3) = 0._db
    elseif( i2 == 12 .or. is == 47 ) then
      v(1) = 1._db; v(2) = 0._db; v(3) = 1._db
    elseif( i2 == 13 .or. is == 44 ) then
      v(1) = -1._db; v(2) = 0._db; v(3) = 1._db
    elseif( i2 == 14 .or. is == 46 ) then
      v(1) = 0._db; v(2) = 1._db; v(3) = 1._db
    elseif( i2 == 15 .or. is == 43 ) then
      v(1) = 0._db; v(2) = -1._db; v(3) = 1._db
    elseif( i2 == 22 .or. is == 40 )then
      v(1) = 1._db; v(2) =  0._db; v(3) = 0._db
    elseif( i2 == 23 .or. is == 41 )then
      v(1) = 0._db; v(2) =  1._db; v(3) = 0._db
    elseif( i2 == 58 .or. is == 61 ) then
      v(1) = 0.866025403784439_db; v(2) = 0.5_db; v(3) = 0._db
    elseif( i2 == 60 .or. is == 63 ) then
      v(1) = 0.5_db; v(2) = 0.866025403784439_db; v(3) = 0._db
    elseif( i2 == 62 .or. is == 57 ) then
      v(1) = -0.5_db; v(2) = 0.866025403784439_db; v(3) = 0._db
    elseif( i2 == 64 .or. is == 59 ) then
      v(1) = -0.866025403784439_db; v(2) = 0.5_db; v(3) = 0._db
    endif

! Axes 4 horizontaux
  elseif( ( iops(18) == 0 .and. iops(28) == 0 ) .and. ( iops(16) == 1 .or. iops(26) == 1 ) ) then
    v(1) = 1._db; v(2) =  0._db; v(3) = 0._db

  elseif( ( iops(18) == 0 .and. iops(28) == 0 ) .and. ( iops(17) == 1 .or. iops(27) == 1 ) ) then
    v(1) = 0._db; v(2) =  1._db; v(3) = 0._db

! Abaissement depuis cubique
  elseif( iops(24) == 0 .and. ( iops(2) == 1 .or. iops(32) == 1 ) ) then
    v(1) = 1._db; v(2) = 1._db; v(3) = 1._db

  elseif( iops(24) == 0 .and. ( iops(4) == 1 .or. iops(34) == 1 ) ) then
    v(1) = 1._db; v(2) = -1._db; v(3) = 1._db

  elseif( iops(24) == 0 .and. ( iops(6) == 1 .or. iops(36) == 1 ) ) then
    v(1) = -1._db; v(2) = 1._db; v(3) = 1._db

  elseif( iops(24) == 0 .and. ( iops(8) == 1 .or. iops(38) == 1 ) ) then
    v(1) = 1._db; v(2) = 1._db; v(3) = -1._db

! Autour de Ox
  elseif( ( ( iops(41) == 0 .and. iops(42) == 0 ) .and. ( iops(43) == 1 .or. iops(46) == 1 ) ) .or. &
      ( ( iops(23) == 0 .and. iops(24) == 0 ) .and. ( iops(14) == 1 .or. iops(15) == 1 ) ) ) then
      v(1) = 0._db; v(2) = 1._db; v(3) = 1._db

! Autour de Oy
  elseif( ( ( iops(42) == 0 .and. iops(40) == 0 ) .and. ( iops(44) == 1 .or. iops(47) == 1 ) ) .or. &
      ( ( iops(22) == 0 .and. iops(24) == 0 ) .and. ( iops(12) == 1 .or. iops(13) == 1 ) ) ) then
      v(1) = 1._db; v(2) = 0._db; v(3) = 1._db

  endif

  vn = sqrt( sum(v(:)**2) )

  if( vn > eps10 ) then

    v(:) = v(:) / vn
    call prodvec(w,v,k)
    wn = sqrt( sum(w(:)**2) )
    w(:) = w(:) / wn
    call prodvec(x,v,w)
    rotmat(1,:) = w(:)
    rotmat(2,:) = x(:)
    rotmat(3,:) = v(:)

! Rotations autour de Oz

! Rotation de + 45 degres en cas D2d ou D2h
  elseif( iops(49) == 0 .and. iops(2) == 0 .and. ( ( iops(40) == 1 .and. iops(41) == 1 .and. &
                iops(45) == 0 .and. iops(48) == 0 .and. iops(22) == 0 .and. iops(23) == 0 .and. &
                iops(10) == 1 .and. iops(11) == 1  ) .or. ( iops(40) == 0 .and. iops(41) == 0 .and. &
                iops(22) == 0 .and. iops(23) == 0 .and. ( iops(45) == 1 .or. iops(48) == 1 .or. &
                iops(10) == 1 .or. iops(11) == 1 ) ) ) ) then
    rotmat(1,1) =   r2;  rotmat(1,2) =   r2; rotmat(1,3) = 0._db
    rotmat(2,1) =  -r2;  rotmat(2,2) =   r2; rotmat(2,3) = 0._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 1._db

! Rotation de 30 degres en cas d'axe 3 selon Oz
  elseif( ( iops(58) == 1 .and. iops(22) == 0 .and. iops(57) == 0 ) .or. &
      ( iops(49) == 1 .and. iops(41) == 1 .and. iops(22) == 1 .and. iops(57) == 0 .and. iops(58) == 0 ) .or. &
      ( iops(57) == 0 .and. iops(61) == 0 .and. iops(59) == 1 .and. iops(63) == 1 ) .or. &
      ( iops(58) == 0 .and. iops(62) == 0 .and. iops(60) == 1 .and. iops(64) == 1 ) ) then

    cs = cos( pi / 6 )
    sn = sin( pi / 6 )
    rotmat(1,1) =   cs;  rotmat(1,2) =  -sn; rotmat(1,3) = 0._db
    rotmat(2,1) =   sn;  rotmat(2,2) =   cs; rotmat(2,3) = 0._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 1._db

  elseif( ( iops(58) == 1 .and. iops(62) == 1 .and. iops(60) == 0 .and. iops(64) == 1 ) .or. &
          ( iops(57) == 1 .and. iops(61) == 1 .and. iops(59) == 0 .and. iops(63) == 0 ) ) then

    cs = cos( pi / 6 )
    sn = - sin( pi / 6 )
    rotmat(1,1) =   cs;  rotmat(1,2) =  -sn; rotmat(1,3) = 0._db
    rotmat(2,1) =   sn;  rotmat(2,2) =   cs; rotmat(2,3) = 0._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 1._db

  elseif( n_opsym_s == 3 .and. iops(49) == 1 .and. iops(59) == 1) then

    cs = cos( pi / 6 )
    sn = sin( pi / 6 )
    rotmat(1,1) =   cs;  rotmat(1,2) =  -sn; rotmat(1,3) = 0._db
    rotmat(2,1) =   sn;  rotmat(2,2) =   cs; rotmat(2,3) = 0._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 1._db

! Groupes magnetiques

! Rotation de + 45 degres en cas 4'/mm'm' --> 4'/mmm
  elseif(  iopsym(18) == -1 .and. iopsym(25) == 1 &
     .and. iopsym(40) == -1 .and. iopsym(41) == -1 .and. iopsym(45) == 1 ) then
    rotmat(1,1) =   r2;  rotmat(1,2) =   r2; rotmat(1,3) = 0._db
    rotmat(2,1) =  -r2;  rotmat(2,2) =   r2; rotmat(2,3) = 0._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 1._db

! Rotation de + 45 degres en cas 4'm'm --> 4'mm'   (le deuxieme "m" est le plan diagonal)
  elseif( n_opsym == 8 .and. iopsym(18) == -1  &
     .and. iopsym(40) == -1 .and. iopsym(41) == -1 .and. iopsym(45) == 1 ) then
    rotmat(1,1) =   r2;  rotmat(1,2) =   r2; rotmat(1,3) = 0._db
    rotmat(2,1) =  -r2;  rotmat(2,2) =   r2; rotmat(2,3) = 0._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 1._db

! mm'm' --> m'm'm
  elseif(  iopsym(40) == 1 .and. iopsym(41) == -1 .and. iopsym(42) == -1 .and. iopsym(52) == 0 ) then
    rotmat(1,1) = 0._db;  rotmat(1,2) = 0._db; rotmat(1,3) = -1._db
    rotmat(2,1) = 0._db;  rotmat(2,2) = 1._db; rotmat(2,3) = 0._db
    rotmat(3,1) = 1._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 0._db

! m'mm' --> m'm'm
  elseif(  iopsym(40) == 1 .and. iopsym(41) == -1 .and. iopsym(42) == -1 .and. iopsym(52) == 0 ) then
    rotmat(1,1) = 1._db;  rotmat(1,2) = 0._db; rotmat(1,3) = 0._db
    rotmat(2,1) = 0._db;  rotmat(2,2) = 0._db; rotmat(2,3) = -1._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 1._db; rotmat(3,3) = 0._db

! m'mm --> mmm'
  elseif(  iopsym(40) == -1 .and. iopsym(41) == 1 .and. iopsym(42) == 1 .and. iopsym(52) == 0 ) then
    rotmat(1,1) = 0._db;  rotmat(1,2) = 0._db; rotmat(1,3) = -1._db
    rotmat(2,1) = 0._db;  rotmat(2,2) = 1._db; rotmat(2,3) = 0._db
    rotmat(3,1) = 1._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 0._db

! mm'm --> mmm'
  elseif(  iopsym(40) == 1 .and. iopsym(41) == -1 .and. iopsym(42) == 1 .and. iopsym(52) == 0 ) then
    rotmat(1,1) = 1._db;  rotmat(1,2) = 0._db; rotmat(1,3) = 0._db
    rotmat(2,1) = 0._db;  rotmat(2,2) = 0._db; rotmat(2,3) = -1._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 1._db; rotmat(3,3) = 0._db

! 2'mm' --> 2'm'm
  elseif( n_opsym == 4 .and. iopsym(24) == -1 .and. iopsym(40) == 1 .and. iopsym(41) == -1 ) then
    rotmat(1,1) = 0._db;  rotmat(1,2) = 1._db; rotmat(1,3) = 0._db
    rotmat(2,1) = -1._db;  rotmat(2,2) = 0._db; rotmat(2,3) = 0._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 1._db

! 2'2'2 --> 22'2'   (ordre z,y,x)
  elseif( n_opsym == 4 .and. iopsym(22) == 1 .and. iopsym(23) == -1 .and. iopsym(24) == -1 ) then
    rotmat(1,1) = 0._db;  rotmat(1,2) = 0._db; rotmat(1,3) = -1._db
    rotmat(2,1) = 0._db;  rotmat(2,2) = 1._db; rotmat(2,3) = 0._db
    rotmat(3,1) = 1._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 0._db

! 2'22' --> 22'2'   (ordre z,y,x)
  elseif( n_opsym == 4 .and. iopsym(22) == -1 .and. iopsym(23) == 1 .and. iopsym(24) == -1 ) then
    rotmat(1,1) = 1._db;  rotmat(1,2) = 0._db; rotmat(1,3) = 0._db
    rotmat(2,1) = 0._db;  rotmat(2,2) = 0._db; rotmat(2,3) = -1._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 1._db; rotmat(3,3) = 0._db

! Rotation de 30 degres (6'mm --> 6'mm' et famille)
  elseif( iopsym(40) == -1 .and. iopsym(41) == 1 .and. iopsym(51) == -1 ) then
    cs = cos( pi / 6 )
    sn = sin( pi / 6 )
    rotmat(1,1) =   cs;  rotmat(1,2) =  -sn; rotmat(1,3) = 0._db
    rotmat(2,1) =   sn;  rotmat(2,2) =   cs; rotmat(2,3) = 0._db
    rotmat(3,1) = 0._db;  rotmat(3,2) = 0._db; rotmat(3,3) = 1._db

  endif

  return
end

!***********************************************************************

! Sous-programme calculant la symetrie ponctuelle de l'agregat

subroutine sym_cluster(Atom_with_axe,Axe_atom_clu,Base_ortho,dcosxyz,iaabs,igroup,iopsym,itype,itypep, &
             natomp,ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,pos,Symmol)

  use declarations
  implicit none

  integer:: ia, iaabs, ib, igr, igra, igrb, io, iok, is, isg, isp, ispin, it, ita, itb, na, natomp, nb, ngroup, ngroup_m, nla, &
    nlatm, nlb, nspin, ntype

  integer, dimension(natomp):: igroup, itypep
  integer, dimension(ngroup):: itype
  integer, dimension(0:ntype):: itequ, nlat, numat
  integer, dimension(ngroup):: itypegen
  integer, dimension(nopsm):: iopsym

  logical:: Base_ortho, Symmol
  logical, dimension(0:ngroup_m):: Atom_with_axe

  real(kind=db):: dpop, vnorme

  real(kind=db), dimension(3):: Axe_atom_c, Axe_atom_s, dcosxyz, ps, pt, vspin, vspini, wspin, wspini
  real(kind=db), dimension(3,3):: matopsym
  real(kind=db), dimension(ngroup,nlatm,nspin) :: popats
  real(kind=db), dimension(3,natomp):: Axe_atom_clu, pos

  do it = 0,ntype
    itequ(it) = it
  end do

  do ita = 0,ntype
    na = numat(ita)
    nla = nlat(ita)
    do igra = 1,ngroup
      if( abs(itype(igra)) == ita ) exit
    end do
    if( igra > ngroup ) cycle
    do itb = ita+1,ntype
      if( itequ(itb) /= itb ) cycle
      nb = numat(itb)
      if( nb /= na ) cycle
      nlb = nlat(itb)
      if( nlb /= nla ) cycle
      do igrb = 1,ngroup
        if( abs(itype(igrb)) == itb ) exit
      end do
      if( igrb > ngroup ) cycle

      dpop = 0._db
      do io = 1,nlb
        do isp = 1,nspin
          dpop = dpop + abs( popats(igra,io,isp) - popats(igrb,io,isp) )
        end do
      end do
      if( dpop < eps6 ) then
        itequ(itb) = ita
      elseif( nspin > 1 ) then
        dpop = 0._db
        do io = 1,nlb
          do isp = 1,nspin
            dpop = dpop + abs( popats(igra,io,isp) - popats(igrb,io,nspin-isp+1) )
          end do
        end do
        if( dpop < eps6 ) itequ(itb) = -ita
      endif

    end do
  end do

  do igr = 1,ngroup
    it = abs( itype(igr) )
    itypegen(igr) = itequ(it)
  end do

  iopsym(1) = 1; iopsym(2:nopsm) = 0

  boucle_sym: do is = 2,nopsm

    call opsym(is,matopsym)

    boucle_exter: do ispin = 1,nspin

      if( ispin == 1 ) then
        isg = 1
      else
        isg = - 1
      endif

      do ia = 1,natomp
        ita = itypep(ia)
        igra = igroup(ia)
        ps(:) = pos(:,ia)
        pt = matmul( matopsym, ps )
        iok = 0
        if( ngroup_m > 0 ) then
          if( Atom_with_axe(igra) ) then
            Axe_atom_c(:) = Axe_atom_clu(:,ia)
            if( abs(Axe_atom_c(1)) < eps6 ) then
              wspin(1) = 1._db;  wspin(2:3) = 0._db
            else
              wspin(3) = 1._db;  wspin(1:2) = 0._db
            endif
            call prodvec( vspin, Axe_atom_c, wspin )
            vspin = vspin / vnorme(Base_ortho,dcosxyz,vspin)

            call prodvec( wspin, Axe_atom_c, vspin )

            vspini = matmul( matopsym, vspin )
            wspini = matmul( matopsym, wspin )
            call prodvec(Axe_atom_s,vspini,wspini)
            if( ispin == 2 ) Axe_atom_s(:) = - Axe_atom_s(:)

          endif
        endif

        do ib = 1,natomp
          if( .not. Symmol ) then
            if( ( ia == iaabs .and. ib /= iaabs ) .or. ( ia /= iaabs .and. ib == iaabs ) ) cycle
          endif
          if( abs( pos(1,ib) - pt(1) ) > epspos .or. abs( pos(2,ib) - pt(2) ) > epspos .or. &
              abs( pos(3,ib) - pt(3) ) > epspos ) cycle
          igrb = igroup(ib)
          if( abs(itypegen(igrb)) /= abs(itypegen(igra)) ) cycle
!     &                                               cycle boucle_sym

          if( ngroup_m == 0 ) then
            iok = 1
            exit
          endif
          if( .not. ( Atom_with_axe(igra) .or. Atom_with_axe(igrb) ) ) then
            iok = 1
            exit
          endif

          if( ( Atom_with_axe(igra) .and. .not. Atom_with_axe(igrb)) &
          .or.(.not. Atom_with_axe(igra) .and. Atom_with_axe(igrb))) cycle boucle_exter

          if( Atom_with_axe(igra) ) then
! L'axe de l'atome a ete renverse dans symsite si l'atome est antiferro
!            if( itypegen(igra) == itypegen(igrb) ) then
              if( abs( Axe_atom_clu(1,ib) - Axe_atom_s(1) ) > epspos .or. abs( Axe_atom_clu(2,ib) - Axe_atom_s(2) ) &
                                                    > epspos .or. abs( Axe_atom_clu(3,ib) - Axe_atom_s(3) ) &
                                    > epspos ) cycle boucle_exter
!            else

!              if( abs( Axe_atom_clu(1,ib) + Axe_atom_s(1) ) > epspos .or. abs( Axe_atom_clu(2,ib) + Axe_atom_s(2) ) &
!                                                    > epspos .or. abs( Axe_atom_clu(3,ib) + Axe_atom_s(3) ) &
!                                    > epspos ) cycle boucle_exter
!            endif

          else

            if(   ( abs( Axe_atom_clu(1,ib) - Axe_atom_s(1) ) > epspos &
             .and.  abs( Axe_atom_clu(1,ib) + Axe_atom_s(1) ) > epspos ) &
             .or. ( abs( Axe_atom_clu(2,ib) - Axe_atom_s(2) ) > epspos &
             .and.  abs( Axe_atom_clu(2,ib) + Axe_atom_s(2) ) > epspos ) &
             .or. ( abs( Axe_atom_clu(3,ib) - Axe_atom_s(3) ) > epspos &
             .and.  abs( Axe_atom_clu(3,ib) + Axe_atom_s(3) ) > epspos ) ) cycle boucle_sym

          endif

          iok = 1
          exit

        end do

        if( iok == 0 ) cycle boucle_exter

      end do

      iopsym(is) = isg
      exit

    end do boucle_exter

  end do boucle_sym

  return
end

!*********************************************************************

! Sousprogramme contenant les operations de symetries.

subroutine opsym(is,matopsym)

  use declarations
  implicit none

  integer:: is

  real(kind=db):: matsym1(3,3,6), matsym2(3,3,7:12), matsym3(3,3,13:18), matsym4(3,3,19:24), &
                 matsym5(3,3,25:30), matsym6(3,3,31:36), matsym7(3,3,37:42), matsym8(3,3,43:48), &
                 matsym9(3,3,49:54), matsym10(3,3,55:60), matsym11(3,3,61:nopsm), matopsym(3,3)

!      1  : identite
!      2 : rot 2*pi/3 autour de (1,1,1)
!      3 : rot 4*pi/3 autour de (1,1,1)
!      4 : rot 2*pi/3 autour de (1,-1,1)
!      5 : rot 4*pi/3 autour de (1,-1,1)
!      6 : rot 2*pi/3 autour de (-1,1,1)
  data matsym1/ 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 1.0_db, &
                1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, &
                1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 1.0_db, 0.0_db, 0.0_db, &
                0.0_db, 0.0_db, 1.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, &
                0.0_db, 0.0_db, 1.0_db, -1.0_db, 0.0_db, 0.0_db/
!      7 : rot 4*pi/3 autour de (-1,1,1)
!      8 : rot 2*pi/3 autour de (1,1,-1)
!      9 : rot 4*pi/3 autour de (1,1,-1)
!      10 : rot 2*pi/2 autour de (1,1,0)
!      11 : rot 2*pi/2 autour de (-1,1,0)
!      12 : rot 2*pi/2 autour de (1,0,1)
  data matsym2/ 0.0_db, 0.0_db,-1.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, &
                0.0_db, 0.0_db,-1.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 1.0_db, 0.0_db, 0.0_db, &
                0.0_db,-1.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, &
                0.0_db,-1.0_db, 0.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db, 1.0_db, &
                0.0_db,-1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db/
!      13 : rot 2*pi/2 autour de (-1,0,1)
!      14 : rot 2*pi/2 autour de (0,1,1)
!      15 : rot 2*pi/2 autour de (0,-1,1)
!      16 : C4x, rot 2*pi/4 selon 0x
!      17 : C4y, rot 2*pi/4 selon 0y
!      18 : C4z, rot 2*pi/4 selon 0z
  data matsym3/ 0.0_db, 0.0_db,-1.0_db, 0.0_db,-1.0_db, 0.0_db, -1.0_db, 0.0_db, 0.0_db, -1.0_db, 0.0_db, 0.0_db, &
                0.0_db, 0.0_db, 1.0_db, 0.0_db, 1.0_db, 0.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, &
                0.0_db,-1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 1.0_db, 0.0_db, &
                0.0_db, 0.0_db, 1.0_db, 0.0_db, 1.0_db, 0.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, &
                1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db/
!      19 : -C4x, rot -2*pi/4 selon 0x
!      20 : -C4y, rot -2*pi/4 selon 0y
!      21 : -C4z, rot -2*pi/4 selon 0z
!      22 : rot 2*pi/2 selon 0x
!      23 : rot 2*pi/2 selon 0y
!      24 : rot 2*pi/2 selon 0z
  data matsym4/ 1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, &
                0.0_db, 1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, -1.0_db, 0.0_db, 0.0_db, &
                0.0_db, 0.0_db, 1.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, &
               -1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, -1.0_db, 0.0_db, 0.0_db, &
                0.0_db,-1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db/
!      25 : inversion
!      26 : S4x, rot 2*pi/4 selon 0x et miroir
!      27 : S4y, rot 2*pi/4 selon 0y et miroir
!      28 : S4z, rot 2*pi/4 selon 0z et miroir
!      29 : -S4x, rot -2*pi/4 selon 0x et miroir
!      30 : -S4y, rot -2*pi/4 selon 0y et miroir
  data matsym5/ -1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, -1.0_db, 0.0_db, 0.0_db, &
                 0.0_db, 0.0_db,-1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db,-1.0_db, 0.0_db, &
                -1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, &
                -1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, &
                 0.0_db,-1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db/
!      31 : -S4z, rot -2*pi/4 selon 0z et miroir
!      32 : rot 2*pi/3 autour de (1,1,1) et inversion
!      33 : rot 4*pi/3 autour de (1,1,1) et inversion
!      34 : rot 2*pi/3 autour de (1,-1,1) et inversion
!      35 : rot 4*pi/3 autour de (1,-1,1) et inversion
!      36 : rot 2*pi/3 autour de (-1,1,1) et inversion
  data matsym6/ 0.0_db, 1.0_db, 0.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db,-1.0_db, &
               -1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, &
               -1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, -1.0_db, 0.0_db, 0.0_db, &
                0.0_db, 0.0_db,-1.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, &
                0.0_db, 0.0_db,-1.0_db, 1.0_db, 0.0_db, 0.0_db/
!      37 : rot 4*pi/3 autour de (-1,1,1) et inversion
!      38 : rot 2*pi/3 autour de (1,1,-1) et inversion
!      39 : rot 4*pi/3 autour de (1,1,-1) et inversion
!      40 : plan perpendiculaire a 0x
!      41 : plan perpendiculaire a 0y
!      42 : plan perpendiculaire a 0z
  data matsym7/ 0.0_db, 0.0_db, 1.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, &
                0.0_db, 0.0_db, 1.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, -1.0_db, 0.0_db, 0.0_db, &
                0.0_db, 1.0_db, 0.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, &
                1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 1.0_db, 0.0_db, 0.0_db, &
                0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db/
!      43 : plan diagonal y = z contenant 0x
!      44 : plan diagonal x = z contenant 0y
!      45 : plan diagonal x = y contenant 0z
!      46 : plan diagonal y = -z contenant 0x
!      47 : plan diagonal x = -z contenant 0y
!      48 : plan diagonal x = -y contenant Oz,
  data matsym8/ 1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, &
                0.0_db, 1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db, 0.0_db, 1.0_db, 0.0_db, 0.0_db, &
                0.0_db, 0.0_db, 1.0_db, 1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db,-1.0_db, 0.0_db, &
                0.0_db, 0.0_db,-1.0_db, 0.0_db, 1.0_db, 0.0_db, -1.0_db, 0.0_db, 0.0_db, 0.0_db,-1.0_db, 0.0_db, &
               -1.0_db, 0.0_db, 0.0_db, 0.0_db, 0.0_db, 1.0_db/
!      49 : rot 2*pi/3 autour de 0z
!      50 : rot 4*pi/3 autour de 0z
!      51 : rot 2*pi/6 autour de 0z
!      52 : rot 10*pi/6 autour de 0z
!      53 : rot 2*pi/3 autour de 0z, axe negatif
!      54 : rot 4*pi/3 autour de 0z, axe negatif
  data matsym9/ -0.5_db,               -0.866025403784439_db, 0.0_db, 0.866025403784439_db, -0.5_db,               0.0_db, &
      0.0_db,                0.0_db,               1.0_db,-0.5_db,                0.866025403784439_db, 0.0_db, &
     -0.866025403784439_db, -0.5_db,               0.0_db, 0.0_db,                0.0_db,               1.0_db, &
      0.5_db,               -0.866025403784439_db, 0.0_db, 0.866025403784439_db,  0.5_db,               0.0_db, &
      0.0_db,                0.0_db,               1.0_db, 0.5_db,                0.866025403784439_db, 0.0_db, &
     -0.866025403784439_db,  0.5_db,               0.0_db, 0.0_db,                0.0_db,               1.0_db, &
     -0.5_db,               -0.866025403784439_db, 0.0_db, 0.866025403784439_db, -0.5_db,               0.0_db, &
      0.0_db,                0.0_db,              -1.0_db,-0.5_db,                0.866025403784439_db, 0.0_db, &
     -0.866025403784439_db, -0.5_db,               0.0_db, 0.0_db,                0.0_db,              -1.0_db/
!      55 : rot 2*pi/6 autour de 0z, axe negatif
!      56 : rot 10*pi/6 autour de 0z, axe negatif
!      57 : plan de sym contenant Oz, a phi = 30 degres
!      58 : Axe 2 perpendiculaire a Oz, a phi = 30 degres
!      59 : plan de sym contenant Oz, a phi = 60 degres
!      60 : Axe 2 perpendiculaire a Oz, a phi = 60 degres
  data matsym10/ 0.5_db,               -0.866025403784439_db, 0.0_db, 0.866025403784439_db,  0.5_db,               0.0_db, &
      0.0_db,                0.0_db,              -1.0_db, 0.5_db,                0.866025403784439_db, 0.0_db, &
     -0.866025403784439_db,  0.5_db,               0.0_db, 0.0_db,                0.0_db,              -1.0_db, &
      0.5_db,                0.866025403784439_db, 0.0_db, 0.866025403784439_db, -0.5_db,               0.0_db, &
      0.0_db,                0.0_db,               1.0_db, 0.5_db,                0.866025403784439_db, 0.0_db, &
      0.866025403784439_db, -0.5_db,               0.0_db, 0.0_db,                0.0_db,              -1.0_db, &
     -0.5_db,                0.866025403784439_db, 0.0_db, 0.866025403784439_db,  0.5_db,               0.0_db, &
      0.0_db,                0.0_db,               1.0_db,-0.5_db,                0.866025403784439_db, 0.0_db, &
      0.866025403784439_db,  0.5_db,               0.0_db, 0.0_db,                0.0_db,              -1.0_db/
!      61 : plan de sym contenant Oz, a phi = 120 degres
!      62 : Axe 2 perpendiculaire a Oz, a phi = 120 degres
!      63 : plan de sym contenant Oz, a phi = 150 degres
!      64 : Axe 2 perpendiculaire a Oz, a phi = 150 degres
  data matsym11/ -0.5_db,               -0.866025403784439_db, 0.0_db, &
     -0.866025403784439_db,  0.5_db,               0.0_db, 0.0_db,                0.0_db,               1.0_db, &
     -0.5_db,               -0.866025403784439_db, 0.0_db,-0.866025403784439_db,  0.5_db,               0.0_db, &
      0.0_db,                0.0_db,              -1.0_db, 0.5_db,               -0.866025403784439_db, 0.0_db, &
     -0.866025403784439_db, -0.5_db,               0.0_db, 0.0_db,                0.0_db,               1.0_db, &
      0.5_db,               -0.866025403784439_db, 0.0_db,-0.866025403784439_db, -0.5_db,               0.0_db, &
      0.0_db,                0.0_db,              -1.0_db/

  if( is <= 6 ) then
    matopsym(:,:) = matsym1(:,:,is)
  elseif( is <= 12 ) then
    matopsym(:,:) = matsym2(:,:,is)
  elseif( is <= 18 ) then
    matopsym(:,:) = matsym3(:,:,is)
  elseif( is <= 24 ) then
    matopsym(:,:) = matsym4(:,:,is)
  elseif( is <= 30 ) then
    matopsym(:,:) = matsym5(:,:,is)
  elseif( is <= 36 ) then
    matopsym(:,:) = matsym6(:,:,is)
  elseif( is <= 42 ) then
    matopsym(:,:) = matsym7(:,:,is)
  elseif( is <= 48 ) then
    matopsym(:,:) = matsym8(:,:,is)
  elseif( is <= 54 ) then
    matopsym(:,:) = matsym9(:,:,is)
  elseif( is <= 60 ) then
    matopsym(:,:) = matsym10(:,:,is)
  else
    matopsym(:,:) = matsym11(:,:,is)
  endif

! Car les data du fortran boucle sur les colonnes a l'interieur des
! lignes
  matopsym = transpose( matopsym)

  return
end


!***********************************************************************

subroutine numgrpt(iopsym,igrpt,igrpt_nomag,mpirank,Sgrp)

  use declarations
  implicit none

  integer, parameter:: ngrptm = 32
  integer, parameter:: ngrptmagm = 90

  integer:: i, igr_sg, igrpt, igrpt_nomag, ipr, is, k, mpirank, nb_ord, ni

  integer, dimension(nopsm):: iopsymc, iopsym

  character(len=8) PointGroup_name, PointGroup_nomag_name, PointGroup_Sch, PointSubGroup_name, &
                   PointSubGroup_Sch, ptgrname_int, ptgrname_int_nomag, ptgrname_sch

  logical:: Sgrp

  common/PointGroup_name/ PointGroup_name, PointGroup_nomag_name
  common/PointGroup_Sch/ PointGroup_Sch
  common/PointSubGroup_name/ PointSubGroup_name, PointSubGroup_Sch

! La sphere est transformee en m3m
  if( sum( abs(iopsym(:)) ) == 64 ) iopsym(49:64) = 0
! Le cylindre est transformee en 6/mmm
  if( abs(iopsym(51)) == 1 .and. abs(iopsym(45)) == 1 ) then
    iopsym(10:11) = 0
    iopsym(18:21:3) = 0
    iopsym(28:31:3) = 0
    iopsym(45:48:3) = 0
  endif

  if( sgrp ) then
    ni = 2
  else
    ni = 1
  endif

  do i = 1,ni

    boucle_grpt: do igrpt = 1,ngrptmagm

      call grp_opsym(igr_sg,igrpt,igrpt_nomag,iopsymc,nb_ord)

      do k = 1,nopsm
        if( iopsym(k) /= iopsymc(k) ) cycle boucle_grpt
      end do

      PointGroup_name = ptgrname_int(igrpt)
      if( igrpt > ngrptm ) then
        PointSubGroup_name = ptgrname_int_nomag(igr_sg)
        PointSubGroup_Sch = ptgrname_Sch(igr_sg)
      endif
      PointGroup_nomag_name = ptgrname_int_nomag(igrpt_nomag)
      PointGroup_Sch = ptgrname_sch(igrpt_nomag)

      return

    end do boucle_grpt

! Cas des groupes magnetiques ou on prend un sous-groupe pour avoir
! l'axe du spin conforme.
    if( Sgrp .and. i == 1 ) then
      is = nopsm
      if( abs( iopsym(25) ) == 1 ) then
        is = 25
      elseif( iopsym(42) == 1 ) then
        is = 42
      elseif( iopsym(24) == 1 ) then
        is = 24
      endif
      iopsym(2:is-1) = 0
      iopsym(is+1:nopsm) = 0
    endif

  end do

  if( mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,'(/A)') ' Point group not found !'
      write(ipr,170)
      call write_iopsym(iopsym,ipr)
    end do
    stop
  endif

  return
  170 format(/' iopsymc =')
end

!***********************************************************************

subroutine grp_opsym_imp(iopsymc,igrpt,igrpt_nomag,mpirank,PointGroup)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( ngrptm = 32, ngrptmagm = 90)

  character(len=8) PointGroup, PointGroup_name, PointGroup_nomag_name, PointGroup_Sch, ptgrname_int, ptgrname_sch

  integer, dimension(nopsm):: iopsymc, iopsymt

  common/PointGroup_name/ PointGroup_name, PointGroup_nomag_name
  common/PointGroup_Sch/ PointGroup_Sch

  if( PointGroup == 'SousGr' ) then

    boucle_grpt1: do igrptt = 1,ngrptmagm

      call grp_opsym(igr_sg,igrptt,igrpt_nomagt,iopsymt,nb_ord)

      if( igrpt_nomagt /= igrpt_nomag ) cycle

      do k = 1,nopsm
        if( iopsymt(k) == 0 ) cycle
        if( iopsymt(k) /= iopsymc(k) ) cycle boucle_grpt1
      end do

      exit

    end do boucle_grpt1

    igrpt = igrptt
    igrpt_nomag = igrpt_nomagt
    iopsymc(:) = iopsymt(:)

    PointGroup_name = ptgrname_int(igrpt)
    PointGroup_nomag_name = ptgrname_int(igrpt_nomag)
    PointGroup_Sch = ptgrname_sch(igrpt_nomag)

    return

  else

    boucle_grpt: do igrpt = 1,ngrptmagm

      if( PointGroup /= ptgrname_int(igrpt) ) then
        if( igrpt > ngrptm ) cycle
        if( PointGroup /= ptgrname_sch(igrpt) ) cycle
      endif

      call grp_opsym(igr_sg,igrpt,igrpt_nomag,iopsymc,nb_ord)

      PointGroup_name = ptgrname_int(igrpt)
      PointGroup_nomag_name = ptgrname_int(igrpt_nomag)
      PointGroup_Sch = ptgrname_sch(igrpt_nomag)

      return

    end do boucle_grpt

  endif

  if( mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,'(/A)') ' Point group not found !'
    end do
    stop
  endif

  return
end

!***********************************************************************

subroutine grp_opsym(igr_sg,igrpt,igrpt_nomag,iopsymc,nb_ord)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( ngrptm = 32, ngrpt_compm = 11, ngrptmagm = 90, ngrptmag_compm = 10 )

  integer, dimension(ngrptm+1:ngrptmagm+ngrptmag_compm):: ngref, ngnmg
  integer, dimension(nopsm):: iopsymc

! Numero de groupe non magnetique reference pour les groupes magnetiques
  data ngref/ 2, 4, 3, 5, 5, 5, 7, 6, 6, 8, 8, 8,19,18,21,24,24,24, 9,10,14,14,11,11,11,12,12,13,13,13, &
             15,15,15,15,15,22,17,20,20,20,26,26,23,23,23,25,25,27,27,27,27,27,29,30,31,32,32,32, &
              7, 7, 6, 8, 8, 8, 8,12,26,25/

! Numero du sous-groupe non magnetique non multiplie par le renversement
! du temps pour les groupes magnetiques
  data ngnmg/ 1, 1, 1, 2, 3, 4, 4, 4,34, 7, 6, 5,16,16,16,21,18,41, 4, 4, 9, 7, 9,10, 5, 9, 6,10, 7,42, &
             14,12, 8,13,11,16,16,19,18,17,22,19,22,17,21,22,18,24,20,23,26,25,28,28,28,31,30,29, &
             35,36,33,39,40,37,38,42,41,43/

  iopsymc(:) = 0

  if( igrpt > ngrptm ) then
    igrpt_nomag = ngref(igrpt)
  else
    igrpt_nomag = igrpt
  endif

  nb_ord = numbops( igrpt_nomag )

  ideb = 0
  do is = 1,igrpt_nomag - 1
    ideb = ideb + numbops(is)
  end do

  do is = 1,nb_ord
    iopsymc( ios(ideb+is) ) = 1
  end do

! Groupes magnetiques : on met un signe negatif pour les operations
! multipliees par le renversement du temps.

  if( igrpt > ngrptm ) then

! Sous-groupe contenant les operations non multiplies par le
! renversement du temps
    igr_sg = ngnmg(igrpt)
    nb_ord_sg = numbops(igr_sg)
    jdeb = 0
    do is = 1,igr_sg - 1
      jdeb = jdeb + numbops(is)
    end do
    boucle_is: do is = 1,nb_ord
      do js = 1,nb_ord_sg
        if( ios(ideb+is) == ios(jdeb+js) ) cycle boucle_is
      end do
      iopsymc( ios(ideb+is) ) = - iopsymc( ios(ideb+is) )
    end do boucle_is

  endif

  return
end

!***********************************************************************

subroutine character_table(icheck,igrpt_nomag,karact,nb_rep)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( nclasm=12, ngrptm = 32, ngrpt_compm = 11 )

  character(len=8):: PointGroup_name, PointGroup_nomag_name
  character(len=9):: nomsym
  character(len=36), dimension(12):: ch
  character(len=144):: mot

  complex(kind=db), dimension(nclasm,nrepm):: kar
  complex(kind=db), dimension(nopsm,nrepm):: karact

  integer, dimension(ngrptm):: nb_cla_table, nb_rep_table
  integer, dimension(nopsm):: nb_sympcl(nclasm)

  logical inversion, jexp, kexp

  common/PointGroup_name/ PointGroup_name, PointGroup_nomag_name

  data nb_cla_table/1,2,2,2,4,4,4,8,4,4, 8,5,5,5,10,3,6,3,3,6, 6,6,12,6,6,6,12,4,6,5, 5,10/

  data nb_rep_table/1,2,2,2,4,4,4,8,3,3, 6,5,5,5,10,2,4,3,3,6, 4,4,8,6,6,6,12,3,6,5, 5,10/

  nb_cla = nb_cla_table(igrpt_nomag)
  nb_rep = nb_rep_table(igrpt_nomag)

  inversion = .false.

  select case(PointGroup_nomag_name)

    case('1')                      ! C1 1
      ch(1) = '  1'

    case('-1','m','2')             ! Ci, Cs, C2, 2 3 4
      inversion = .true.
      ch(1) = '  1'

    case('2/m','mm','222')         ! C2h, C2v, D2 5 6 7
      inversion = .true.
      ch(1) = '  1  1'
      ch(2) = '  1 -1'

    case('mmm')                    ! D2h 8
      inversion = .true.
      ch(1) = '  1  1  1  1'
      ch(2) = '  1  1 -1 -1'
      ch(3) = '  1 -1  1 -1'
      ch(4) = '  1 -1 -1  1'

    case('4','-4')                 ! C4, S4 9 10
      ch(1) = '  1  1  1  1'
      ch(2) = '  1  1 -1 -1'
      ch(3) = '  1 -1  i -i'

    case('4/m')                    ! C4h 11
      inversion = .true.
      ch(1) = '  1  1  1  1'
      ch(2) = '  1 -1  1 -1'
      ch(3) = '  1  i -1 -i'

    case('4mm','-42m','422')       ! C4v, D2d, D4   12 13 14
      ch(1) = '  1  1  1  1  1'
      ch(2) = '  1  1  1 -1 -1'
      ch(3) = '  1  1 -1  1 -1'
      ch(4) = '  1  1 -1 -1  1'
      ch(5) = '  2 -2  0  0  0'

    case('4/mmm')                  ! D4h  15
      inversion = .true.
      ch(1) = '  1  1  1  1  1'
      ch(2) = '  1  1  1 -1 -1'
      ch(3) = '  1 -1  1  1 -1'
      ch(4) = '  1 -1  1 -1  1'
      ch(5) = '  2  0 -2  0  0'

    case('3')                      ! C3  16
      ch(1) = '  1  1  1'
      ch(2) = '  1  j j2'

    case('-3','-6','6')            ! S6, C3h, C6   17 21 22
      inversion = .true.
      ch(1) = '  1  1  1'
      ch(2) = '  1  j j2'

    case('3m','32')                ! C3v D3  18 19
      ch(1) = '  1  1  1'
      ch(2) = '  1  1 -1'
      ch(3) = '  2 -1  0'

    case('-3m','-6m2','6mm','622') ! D3d D3h C6v D6  20 24 25 26
      inversion = .true.
      ch(1) = '  1  1  1'
      ch(2) = '  1  1 -1'
      ch(3) = '  2 -1  0'

    case('6/m')                    ! C6h  23
      inversion = .true.
      ch(1) = '  1  1  1  1  1  1'
      ch(2) = '  1 -1  1 -1  1 -1'
      ch(3) = '  1  k  j -1 j2 k5'
      ch(4) = '  1  j j2  1  j j2'

    case('6/mmm')                  ! D6h 27
      inversion = .true.
      ch(1) = '  1  1  1  1  1  1'
      ch(2) = '  1  1  1  1 -1 -1'
      ch(3) = '  1 -1  1 -1  1 -1'
      ch(4) = '  1 -1  1 -1 -1  1'
      ch(5) = '  2  1 -1 -2  0  0'
      ch(6) = '  2 -1 -1  2  0  0'

    case('23')                     ! T 28
      ch(1) = '  1  1  1  1'
      ch(2) = '  1  1  j j2'
      ch(3) = '  3 -1  0  0'

    case('m3')                     ! Th 29
      inversion = .true.
      ch(1) = '  1  1  1  1'
      ch(2) = '  1  j j2  1'
      ch(3) = '  3  0  0 -1'

    case('-43m','432')             ! Td, O  30 31
      ch(1) = '  1  1  1  1  1'
      ch(2) = '  1  1  1 -1 -1'
      ch(3) = '  2 -1  2  0  0'
      ch(4) = '  3  0 -1 -1  1'
      ch(5) = '  3  0 -1  1 -1'

    case('m3m')                    ! Oh 32
      inversion = .true.
      ch(1) = '  1  1  1  1  1'
      ch(2) = '  1  1 -1 -1  1'
      ch(3) = '  2 -1  0  0  2'
      ch(4) = '  3  0 -1  1 -1'
      ch(5) = '  3  0  1 -1 -1'

  end select

  if( inversion ) then
    nc = nb_cla / 2
    nr = nb_rep / 2
  else
    nc = nb_cla
    nr = nb_rep
  endif

  if( inversion ) then
    do irep = 1,nr
      mot = ch(irep)
      l = len_trim(mot)
      mot(l+1:2*l) = mot(1:l)
      ch(irep) = mot
      do ic = nc+1, nb_cla
        j = 3*ic-1
        select case( mot(j:j) )
          case(' ')
            mot(j:j) = '-'
          case('-')
            mot(j:j) = ' '
          case('j')
            mot(j:j+1) = ' k'
          case('k')
            mot(j:j+1) = ' j'
        end select
      end do
      ch(irep+nr) = mot
    end do
  endif

  r = sqrt(3._db)/2

  do j = 1,nr
    jc = -3
    mot = ch(j)
    do i = 1,nc
      jc = jc + 3

      select case(mot(jc+2:jc+3))
        case(' 0')
          kar(i,j) = (0._db,0._db)
        case(' 1')
          kar(i,j) = (1._db,0._db)
        case('-1')
          kar(i,j) = (-1._db,0._db)
        case(' 2')
          kar(i,j) = (2._db,0._db)
        case('-2')
          kar(i,j) = (-2._db,0._db)
        case(' 3')
          kar(i,j) = (3._db,0._db)
        case('-3')
          kar(i,j) = (-3._db,0._db)
        case(' i')
          kar(i,j) = (0._db,1._db)
        case('-i')
          kar(i,j) = (0._db,-1._db)
        case(' j')
          kar(i,j) = cmplx(-0.5_db,r,db)
        case('j2')
          kar(i,j) = cmplx(-0.5_db,-r,db)
        case(' k')
          kar(i,j) = cmplx(0.5_db,r,db)
        case('k5')
          kar(i,j) = cmplx(0.5_db,-r,db)
      end select
    end do
  end do

  if( inversion ) then
    do j = 1,nc
      kar(j+nc,1:nr) = kar(j,1:nr)
    end do
    do i = 1,nr
      kar(1:nc,i+nr) = kar(1:nc,i)
      kar(nc+1:nb_cla,i+nr) = - kar(nc+1:nb_cla,i)
    end do
  endif

  nb_sympcl(1:nc) = 1

  select case(PointGroup_nomag_name)

    case('4mm  ','-42m ','422  ')
      nb_sympcl(3:5) = 2

    case('4/mmm')
      nb_sympcl(2) = 2; nb_sympcl(4:5) = 2

    case('3m','32','-3m','-6m2','6mm','622')
      nb_sympcl(2) = 2; nb_sympcl(3) = 3

    case('23')
      nb_sympcl(2) = 3; nb_sympcl(3) = 4; nb_sympcl(4) = 3

    case('m3')
      nb_sympcl(2) = 4; nb_sympcl(3) = 4; nb_sympcl(4) = 3

    case('-43m','432')
      nb_sympcl(2) = 8; nb_sympcl(3) = 3; nb_sympcl(4:5) = 6

    case('m3m')
      nb_sympcl(2) = 8; nb_sympcl(3:4) = 6; nb_sympcl(5) = 3

  end select

  if( inversion ) then
    do j = 1,nc
      nb_sympcl(j+nc) = nb_sympcl(j)
    end do
  endif

  ideb = 0
  do is = 1,igrpt_nomag - 1
    ideb = ideb + numbops(is)
  end do

  karact(:,:) = (0._db,0._db)

  boucle_is: do is = 1,numbops( igrpt_nomag )

    k= 0
    do i = 1,nb_cla
      do j = 1,nb_sympcl(i)
        k = k + 1
        if( k == is ) then
          karact(ios(ideb+is),1:nb_rep) = kar(i,1:nb_rep)
          cycle boucle_is
        endif
      end do
    end do

  end do boucle_is

  if( icheck > 0 ) then
    ifin = 0
    do is = 1,igrpt_nomag
      ifin = ifin + numbops(is)
    end do
    write(3,110)
    write(3,'(5x,48i3)') ( ios(is), is = ideb+1,ifin )
    do irep = 1,nb_rep
      mot = ch(irep)
      ic = 0
      do i = 1,nb_cla
        ic = ic + 3
        if( nb_sympcl(i) > 1 ) then
          id = 3 * nb_sympcl(i) - 3
          l = len_trim(mot)
          mot(ic+1+id:l+id) = mot(ic+1:l)
          mot(ic+1:ic+id) = ' '
          ic = ic + id
        endif
      end do
      write(3,'(5x,a)') mot
    end do
    jexp = .false.; kexp = .false.
    do i = 1,len_trim(mot)
      if( mot(i:i) == 'j' ) jexp = .true.
      if( mot(i:i) == 'k' ) kexp = .true.
    end do
    if( jexp .or. kexp ) then
      write(3,'(/A)') ' Symmetry code      Character code'
    else
      write(3,'(/A)') ' Symmetry code:'
    endif
    do is = ideb+1,ifin
      if( is == ideb+1 .and. jexp ) then
        write(3,120) ios(is), nomsym(ios(is))
      elseif( is == ideb+1 .and. kexp ) then
        write(3,130) ios(is), nomsym(ios(is))
      elseif( is == ideb+2 .and. jexp .and. kexp ) then
        write(3,130) ios(is), nomsym(ios(is))
      else
        write(3,140) ios(is), nomsym(ios(is))
      endif
    end do
  endif

  return
  110 format(/' Character table :')
  120 format(2x,i3,a9,'     j = exp(2.i.pi/3)')
  130 format(2x,i3,a9,'     k = exp(2.i.pi/6)')
  140 format(2x,i3,a9)
end

!*********************************************************************

! Sousprogramme calculant la forme des tenseurs

subroutine Tensor_shape(Atom_with_axe,Atom_nonsph, Axe_atom_clu,Base_ortho,dcosxyz,Dipmag, &
               E1E2e,Green,iaabs,icheck,igroup,igrpt0,iopsymc,iopsymr,itype,itypep,ldip,loct,lqua,lseuil,magnetic, &
               msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo, msymooi,msymqq,msymqqi,Multipole,n_oo,natomp,ngroup, &
               ngroup_m,nlat,nlatm,nspin,ntype,numat, octupole,popats,pos,quadrupole,rot_atom_abs,Rot_int, &
               Spinorbite,State_all,Symmol,Tensor_imp)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer:: he, hs, i, icheck, j, k
  integer, dimension(3):: ldip, ldipimp
  integer, dimension(12):: Tensor_imp
  integer, dimension(3,3):: lqua, lquaimp, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3):: msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi
  integer, dimension(nopsm):: irotiops, iopsymc, iopsymr, iopsymt
  integer, dimension(natomp):: igroup, itypep
  integer, dimension(ngroup):: itype
  integer, dimension(0:ntype):: nlat, numat

  logical:: Atom_comp, Atom_comp_cal, Atom_mag, Atom_mag_cal, Atom_nonsph, Base_ortho, Dipmag, E1E1, E1E2, E1E2e, E1E3, &
     E2E2, E3E3, Green, magnet, Magnetic, Octupole, Quadrupole, Spinorbite, State_all, Symmol
  logical, dimension(10):: Multipole
  logical, dimension(0:ngroup_m):: Atom_with_axe

  real(kind=db), dimension(3):: dcosxyz, ps
  real(kind=db), dimension(3,3):: rot_atom_abs, Rot_int, rot_tem
  real(kind=db), dimension(3,natomp):: Axe_atom_clu, pos
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats

  common/iopsym_abs/ iopsym_abs(nopsm)

  if( icheck > 0 ) write(3,100)

  E1E1 = Multipole(1); E1E3 = Multipole(3);
  E2E2 = Multipole(6); E3E3 = Multipole(7);

  msymdd(:,:) = 0
  msymddi(:,:) = 0

  ldipimp(1:3) = Tensor_imp(1:3)
  k = 3
  do i = 1,3
    do j = 1,3
      k = k + 1
      lquaimp(i,j) = Tensor_imp(k)
    end do
  end do

! Evaluation de la symetrie de l'atome absorbeur
  if( Magnetic ) then
    Atom_mag = Atom_mag_cal(igrpt0)
  else
    Atom_mag = .false.
  endif
  Atom_comp = Atom_comp_cal(igrpt0)

  ps(:) = pos(:,iaabs)
  if( nspin == 1 .and. sum( abs( pos(:,iaabs) ) ) < eps10 ) then
    iopsym_abs(:) = iopsymc(:)
  else
    if( spinorbite ) then
      call point_group_atom(Atom_comp,Atom_mag,Atom_with_axe,Atom_nonsph,Axe_atom_clu,Base_ortho,dcosxyz, &
        iaabs,iaabs,icheck,igroup,igroup(iaabs),igrpt,igrpt0,iopsym_abs,iopsymr,itype,itypep,Magnetic,mpirank, &
        natomp,ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,pos,ps,rot_atom_abs,Spinorbite,Symmol)
    else
! On prend la symetrie complete (iopsymc)
      call point_group_atom(Atom_comp,Atom_mag,Atom_with_axe,Atom_nonsph,Axe_atom_clu,Base_ortho,dcosxyz, &
        iaabs,iaabs,icheck,igroup,igroup(iaabs),igrpt,igrpt0,iopsym_abs,iopsymc,itype,itypep,Magnetic,mpirank, &
        natomp,ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,pos,ps,rot_atom_abs,Spinorbite,Symmol)
    endif
  endif
! La rotation est calcule avec la symetrie utilisee (iopsymr)
  call point_group_atom(Atom_comp,Atom_mag,Atom_with_axe,Atom_nonsph,Axe_atom_clu,Base_ortho,dcosxyz, &
      iaabs,iaabs,icheck,igroup,igroup(iaabs),igrpt,igrpt0,iopsymt,iopsymr,itype,itypep,Magnetic,mpirank, &
      natomp,ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,pos,ps,rot_atom_abs,Spinorbite,Symmol)

  if( spinorbite .or. ( nspin == 2 .and. lseuil > 0 ) ) then
    magnet = .true.
  else
    magnet = .false.
  endif

  do ical = 1,2
    if( State_all ) then
      irotiops(:) = 0
      irotiops(1) = 1
    else
      if( ical == 1 ) then
        rot_tem = transpose(matmul(Rot_int,transpose(rot_atom_abs)))
        if( abs( rot_tem(1,1) - 1 ) < eps10 .and. abs( rot_tem(2,2) - 1 ) < eps10 .and. &
            abs( rot_tem(3,3) - 1 ) < eps10 ) cycle
        call iop_rot(icheck,irotiops,rot_tem)
        where( irotiops /= 0 ) irotiops = iopsym_abs(irotiops)
      else
        irotiops = iopsym_abs
      endif
    endif

! Calcul des tenseurs
    if( E1E1 .or. dipmag ) call tensdd(irotiops,magnet,msymdd,msymddi)
    if( E1E2e ) call tensdq(irotiops,magnet,msymdq,msymdqi)
    if( E2E2 ) call tensqq(irotiops,magnet,msymqq,msymqqi)
    if( E1E3 ) call tensdo(irotiops,magnet,msymdo,msymdoi)
    if( E3E3 ) call tensoo(irotiops,magnet,msymoo,msymooi,n_oo)

! Terme d'interference dipole - quadrupole
    E1E2 = E1E2e
    if( E1E2e ) then
      n = maxval( abs(msymdq) )
      n = max( n, maxval( abs(msymdqi) ) )
      if( n == 0 ) E1E2 = .false.
    endif
    Multipole(2) = E1E2

    if( icheck > 0 ) then
      if( ical == 1 ) then
        write(3,102)
      else
        write(3,104)
      endif
      if( magnet ) then
        if( E1E1 .or. dipmag ) then
          write(3,110)
          do i = 1,3
            write(3,120) ( msymdd(i,k), msymddi(i,k), k = 1,3)
          end do
        endif
        if( E1E2e ) then
          write(3,130)
          do j = 1,3
            write(3,120) ( ( msymdq(i,j,k), msymdqi(i,j,k), k =1,3), i = 1,3 )
          end do
        endif
        if( E2E2 ) then
          write(3,140)
          do i = 1,3
            write(3,150) i, i, i
            do k = 1,3
              write(3,120) ( ( msymqq(i,j,k,l), msymqqi(i,j,k,l), l = 1,3), j = 1,3 )
            end do
          end do
        endif
        if( E1E3 ) then
          write(3,160)
          do i = 1,3
            write(3,150) i, i, i
            do k = 1,3
              write(3,120) ( (msymdo(i,j,k,l), msymdoi(i,j,k,l), l = 1,3 ), j = 1,3 )
            end do
          end do
        endif
        if( E3E3 ) then
          write(3,165)
          do hs = 1,3
            do js = 1,3
              jhs = 3 * ( js - 1 ) + hs
              write(3,155) (k, js, hs, k = 1,3)
              do he = 1,3
                do je = 1,3
                  jhe = 3 * ( je - 1 ) + he
                  write(3,167) je, he, ( ( msymoo(i,jhe,k,jhs), msymooi(i,jhe,k,jhs), i = 1,3), k =1,3)
                end do
              end do
            end do
          end do
        endif
      else
        if( E1E1 .or. dipmag ) then
          write(3,110)
          do i = 1,3
            write(3,170) msymdd(i,:)
          end do
        endif
        if( E1E2e ) then
          write(3,180)
          do j = 1,3
            write(3,170) (msymdq(i,j,:), i = 1,3)
          end do
        endif
        if( E2E2 ) then
          write(3,140)
          do i = 1,3
            write(3,190) i, i, i
            do k = 1,3
              write(3,170) (msymqq(i,j,k,:), j = 1,3)
            end do
          end do
        endif
        if( E1E3 ) then
          write(3,160)
          do i = 1,3
            write(3,190) i, i, i
            do k = 1,3
              write(3,170) (msymdo(i,j,k,:), j = 1,3)
            end do
          end do
        endif
        if( E3E3 ) then
          write(3,165)
          do hs = 1,3
            do js = 1,3
              jhs = 3 * ( js - 1 ) + hs
              write(3,195) (k, js, hs, k = 1,3)
              do he = 1,3
                do je = 1,3
                  jhe = 3 * ( je - 1 ) + he
                  write(3,168) je, he, (msymoo(:,jhe,k,jhs), k =1,3)
                end do
              end do
            end do
          end do
        endif
      endif
    endif

  end do

! Determination des operateurs a calculer
  ldip = 0
  lqua = 0
  loct = 0

  if( E1E1 .or. Dipmag ) then
    n = maxval( abs(msymdd) )
    n = max( n, maxval( abs(msymddi) ) )
    do i = 1,n
      boucle_dip: do ke = 3,1,-1
        do ks = 3,1,-1
          if( abs(msymdd(ke,ks)) == i .or. abs(msymddi(ke,ks)) == i ) exit boucle_dip
        end do
      end do boucle_dip
      ldip(ke) = 1
      ldip(ks) = 1
    end do
  endif

  if( E2E2 ) then
    n = maxval(abs(msymqq))
    n = max( n, maxval( abs(msymqqi) ) )
    do i = 1,n
      boucle_qua: do ke = 1,3
        do ks = 1,3
          do je = 1,3
            do js = 1,3
              if( abs(msymqq(ke,je,ks,js)) == i .or. abs(msymqqi(ke,je,ks,js)) == i) exit boucle_qua
            end do
          end do
        end do
      end do boucle_qua
      kk = min(ke,je)
      jj = max(ke,je)
      lqua(kk,jj) = 1
      kk = min(ks,js)
      jj = max(ks,js)
      lqua(kk,jj) = 1
    end do
  endif

! Terme d'interference dipole - quadrupole
  if( E1E2 ) then
    n = maxval(abs(msymdq))
    n = max( n, maxval( abs(msymdqi) ) )
    do i = 1,n
      boucle_int: do ke = 3,1,-1
        do ks = 3,1,-1
          do js = 3,1,-1
            if( abs(msymdq(ke,ks,js)) == i .or. abs(msymdqi(ke,ks,js)) == i ) exit boucle_int
          end do
        end do
      end do boucle_int
      ldip(ke) = 1
      kk = min(ks,js)
      jj = max(ks,js)
      lqua(kk,jj) = 1
    end do
  endif

  if( E1E3 ) then
    n = maxval(abs(msymdo))
    n = max( n, maxval( abs(msymdoi) ) )
    do i = 1,n
      boucle_oct: do ke = 3,1,-1
        do ks = 3,1,-1
          do j1 = 3,1,-1
            do j2 = 3,1,-1
              if( abs(msymdo(ke,ks,j1,j2)) == i .or. abs(msymdoi(ke,ks,j1,j2)) == i ) exit boucle_oct
            end do
          end do
        end do
      end do boucle_oct
      ldip(ke) = 1
      kk = min(ks,j1,j2)
      jj2 = max(ks,j1,j2)
      if( (ks == kk .and. j1 == jj2) .or. (j1 == kk .and. ks == jj2) ) then
        jj1 = j2
      elseif( (ks == kk .and. j2 == jj2) .or. (j2 == kk .and. ks == jj2) ) then
        jj1 = j1
      else
        jj1 = ks
      endif
      loct(kk,jj1,jj2) = 1
    end do
  endif

  if( E3E3 ) then
    n = maxval(abs(msymoo))
    n = max( n, maxval( abs(msymooi) ) )
    do i = 1,n
      boucle_oo: do ke = 3,1,-1
        do je = 3,1,-1
          do he = 3,1,-1
            jhe = 3*(je - 1) + he
            do ks = 3,1,-1
              do js = 3,1,-1
                do hs = 3,1,-1
                  jhs = 3*(js - 1) + hs
                  if( abs(msymoo(ke,jhe,ks,jhs)) == i .or. abs(msymooi(ke,jhe,ks,jhs)) == i ) exit boucle_oo
                end do
              end do
            end do
          end do
        end do
      end do boucle_oo

   end do
     loct(:,:,:) = 1
  endif

  if( icheck > 0 ) write(3,210) ldip(1:3)

  if( ( E1E1 .or. Dipmag ) .and. ldipimp(1) /= -1 .and. .not. green ) then
    ldip(:) = ldipimp(:)
    if( icheck > 0 ) write(3,215) ldip(1:3)
  endif

  if( Quadrupole .or. Octupole ) then
    if( icheck > 0 ) write(3,220) (lqua(i,1:3), i = 1,3)

    if( lquaimp(1,1) /= -1  .and. .not. green ) then
      lqua(:,:) = 0
      do i = 1,3
        do j = 1,3
          if( lquaimp(i,j) == 0 ) cycle
          ii = min(i,j)
          jj = max(i,j)
          lqua(ii,jj) = 1
        end do
      end do
      if( icheck > 0 ) write(3,225) (lqua(i,1:3), i = 1,3)
    endif
  endif

  if( icheck > 0 .and. Octupole ) write(3,230) ( (loct(i,j,1:3), j = 1,3 ), i = 1,3 )

  return
  100 format(//' ---- Tensor_shape --------',100('-'))
  102 format(/'  Internal basis, z along c :')
  104 format(/'  Absorbing atom basis :')
  110 format(/'  Dipole-dipole matrix shape :'/)
  120 format(3(3x,3(i4,i3)))
  130 format(/'  Dipole-quadrupole matrix shape :'/,/ 11x,'(1,j,k)',17x,'(2,j,k)',17x,'(3,j,k)')
  140 format(/'  Quadrupole-quadrupole matrix shape :'/)
  150 format( 10x,'(',i1,',1,k,l)',15x,'(',i1,',2,k,l)',15x,'(',i1,',3,k,l)')
  155 format(/' je he',3(13x,'(:,je,he',3(',',i1),')'))
  160 format(/'  Dipole-Octupole matrix shape :'/)
  165 format(/'  Octupole-Octupole matrix shape :'/)
  167 format(2i3,3(6x,3(i4,i3)))
  168 format(2i3,5x,3(3i4,8x))
  170 format(3(3x,3i3))
  180 format(/'  Dipole-quadrupole matrix shape :'/,/ 5x,'(1,j,k)     (2,j,k)     (3,j,k)')
  190 format(/'    (',i1,',1,k,l)   (',i1,',2,k,l)   (',i1,',3,k,l)')
  195 format(/' je he',3(5x,'(:,je,he',3(',',i1),')'))
  210 format(/' ldip = (',3i2,')')
  215 format(/' ldip = (',3i2,'), apres imposition')
  220 format(/' lqua = (',3i2,')',2(/8x,'(',3i2,')'))
  225 format(/' lqua = (',3i2,'), apres imposition',2(/8x,'(',3i2,')'))
  230 format(/' loct = ',3('(',3i2,')  '),2(/8x,3('(',3i2,')  ')))
end

!*********************************************************************

! Calcul du tenseur dipole-dipole

subroutine tensdd(iopsymt,magnet,msymdd,msymddi)

  use declarations

  integer im(3), iopsymt(nopsm)
  integer, dimension(3,3):: msymdd, msymddi

  logical magnet
  logical, dimension(9):: impos, imposi

  real(kind=db):: matopsym(3,3)

  msymdd(:,:) = 0;  msymddi(:,:) = 0;
  impos(:) = .false.
  if( magnet ) then
   imposi(:) = .false.
  else
   imposi(:) = .true.
  endif

  kkk = 0
  do i = 1,3
    do j = i,3
      if( msymdd(i,j) /= 0 .or. msymddi(i,j) /= 0 ) cycle
      kkk = kkk + 1
      msymdd(i,j) = kkk
      msymdd(j,i) = kkk
      if( i /= j ) then
        msymddi(i,j) = kkk
        msymddi(j,i) = - kkk
      endif

      do iss = 2,min(48,nopsm)

        is = iordresym(iss)
        if( iopsymt(is) == 0 ) cycle
        isgmag = iopsymt(is)
        call opsym(is,matopsym)
        do ii = 1,3
          do jj = 1,3
            if( abs(matopsym(jj,ii)) < 0.0001_db ) cycle
            im(ii) = jj * nint( matopsym(jj,ii) )
            exit
          end do
        end do
        ii = abs( im(i) )
        iis = ii / im(i)
        jj = abs( im(j) )
        jjs = jj / im(j)
        isg = iis * jjs
        isgi = isg * isgmag

        if( ( isg > 0 .and. msymdd(ii,jj) == - kkk ) .or. ( isg < 0 .and. msymdd(ii,jj) == kkk ) ) impos(kkk) = .true.
        if( .not. impos(kkk) ) then
          msymdd(ii,jj) = isg * kkk
          msymdd(jj,ii) = isg * kkk
        endif
        if( ii == jj ) cycle
        if( is /= 51 .and. is /= 52 .and. is < 55 ) then
          if( ( isgi > 0 .and. msymddi(ii,jj) == - kkk ) .or. ( isgi < 0 .and. msymddi(ii,jj) == kkk ) ) &
                                                 imposi(kkk) =.true.
          if( .not. imposi(kkk) .and. isgmag /= 0 ) then
            msymddi(ii,jj) = isgi * kkk
            msymddi(jj,ii) = - isgi * kkk
          endif
        endif

      end do

    end do
  end do

  kkkmax = kkk
  kk = 0
  do kkk = 1,kkkmax
    kk = kk + 1
    if( impos(kkk) ) where( abs(msymdd) == kk ) msymdd = 0
    if( imposi(kkk) ) where( abs(msymddi) == kk ) msymddi = 0
    if( impos(kkk) .and. imposi(kkk) ) then
      where( abs(msymdd) > kk ) msymdd = msymdd - abs(msymdd)/msymdd
      where( abs(msymddi) > kk ) msymddi = msymddi - abs(msymddi) / msymddi
      kk = kk - 1
    endif
  end do

  n = maxval( abs(msymdd) )
  boucle_i: do i = 1,n
    do ke = 1,3
      do ks = 1,3
        if( abs(msymdd(ke,ks)) == i .or. abs(msymddi(ke,ks)) == i ) cycle boucle_i
      end do
    end do
    where( abs(msymdd) > i ) msymdd = msymdd - abs(msymdd) / msymdd
    where( abs(msymddi) > i ) msymddi = msymddi - abs(msymddi) / msymddi
  end do boucle_i

  return
end

!*********************************************************************

function ismag(is)

  select case(is)

    case(1,18,21,24,25,28,31,42,49,50,51,52,53,54,55,56)
      ismag = 1

    case(2,3,4,5,6,7,8,9,12,13,14,15,16,17,19,20,26,27,29,30,32,33, 34,35,36,37,38,39,43,44,46,47)
      ismag = 0

    case(10,11,22,23,40,41,45,48,57,58,59,60,61,62,63,64)
      ismag = -1

  end select

  return
end

!*********************************************************************

subroutine tensdq(iopsymt,magnet,msymdq,msymdqi)

  use declarations
  integer im(3), iopsymt(nopsm)
  integer, dimension(3,3,3):: msymdq, msymdqi

  logical magnet
  logical, dimension(27):: impos, imposi

  real(kind=db):: matopsym(3,3)

  msymdq(:,:,:) = 0; msymdqi(:,:,:) = 0
  impos(:) = .false.
  if( magnet ) then
   imposi(:) = .false.
  else
   imposi(:) = .true.
  endif

  kkk = 0
  do i = 1,3
    do j = 1,3
      do k = j,3
        if( msymdq(i,j,k) /= 0 .or. msymdqi(i,j,k) /= 0 ) cycle
        kkk = kkk + 1
        msymdq(i,j,k) = kkk
        msymdq(i,k,j) = kkk
        msymdqi(i,j,k) = kkk
        msymdqi(i,k,j) = kkk

        do iss = 2,min(48,nopsm)

          is = iordresym(iss)
          if( is > 48 ) cycle
          if( iopsymt(is) == 0 ) cycle
          isgmag = iopsymt(is)
          call opsym(is,matopsym)
          do ii = 1,3
            do jj = 1,3
              if( abs(matopsym(jj,ii)) < 0.0001_db ) cycle
              im(ii) = jj * nint( matopsym(jj,ii) )
              exit
            end do
          end do
          ii = abs( im(i) )
          iis = ii / im(i)
          jj = abs( im(j) )
          jjs = jj / im(j)
          kk = abs( im(k) )
          kks = kk / im(k)
          isg = iis * jjs * kks
          isgi = isg * isgmag

          if( ( isg > 0 .and. msymdq(ii,jj,kk) == - kkk ) .or. ( isg < 0 .and. msymdq(ii,jj,kk) == kkk ) ) &
                                                 impos(kkk) = .true.
          if( ( isgi > 0 .and. msymdqi(ii,jj,kk) == - kkk ).or. ( isgi < 0 .and. msymdqi(ii,jj,kk) == kkk ) ) &
                                                 imposi(kkk) =.true.

          if( .not. impos(kkk) ) then
            msymdq(ii,jj,kk) = isg * kkk
            msymdq(ii,kk,jj) = isg * kkk
          endif
          if( .not. imposi(kkk) .and. isgmag /= 0 ) then
            msymdqi(ii,jj,kk) = isgi * kkk
            msymdqi(ii,kk,jj) = isgi * kkk
          endif

        end do

      end do
    end do
  end do

  kkkmax = kkk

  kk = 0
  do kkk = 1,kkkmax
    kk = kk + 1
    if( impos(kkk) ) where( abs(msymdq) == kk ) msymdq = 0
    if( imposi(kkk) ) where( abs(msymdqi) == kk ) msymdqi = 0
    if( impos(kkk) .and. imposi(kkk) ) then
      where( abs(msymdq) > kk ) msymdq = msymdq - abs(msymdq)/msymdq
      where( abs(msymdqi) > kk ) msymdqi = msymdqi - abs(msymdqi) / msymdqi
      kk = kk - 1
    endif
  end do

  n = maxval( abs(msymdq) )
  boucle_i: do i = 1,n
    do k = 1,3
      do l = 1,3
        do m = 1,3
          if( abs(msymdq(k,l,m)) == i .or. abs(msymdqi(k,l,m)) == i ) cycle boucle_i
        end do
      end do
    end do
    where( abs(msymdq) > i ) msymdq = msymdq - abs(msymdq) / msymdq
    where( abs(msymdqi) > i ) msymdqi = msymdqi - abs(msymdqi) / msymdqi
  end do boucle_i

  return
end

!*********************************************************************

! Calcul du tenseur quadrupole-quadrupole

subroutine tensqq(iopsymt,magnet,msymqq,msymqqi)

  use declarations
  integer im(3), iopsymt(nopsm)
  integer, dimension(3,3,3,3):: msymqq, msymqqi

  logical magnet
  logical, dimension(81):: impos, imposi

  real(kind=db):: matopsym(3,3)

  msymqq(:,:,:,:) = 0; msymqqi(:,:,:,:) = 0
  impos(:) = .false.
  if( magnet ) then
   imposi(:) = .false.
  else
   imposi(:) = .true.
  endif

  kkk = 0
  do i = 1,3
    do j = i,3
      do k = 1,3
        do l = k,3
          if( msymqq(i,j,k,l) /= 0 .or. msymqqi(i,j,k,l) /= 0) cycle
          kkk = kkk + 1
          msymqq(i,j,k,l) = kkk
          msymqq(j,i,k,l) = kkk
          msymqq(i,j,l,k) = kkk
          msymqq(j,i,l,k) = kkk
          msymqq(k,l,i,j) = kkk
          msymqq(l,k,i,j) = kkk
          msymqq(k,l,j,i) = kkk
          msymqq(l,k,j,i) = kkk
          if( i /= k .or. j /= l ) then
            msymqqi(i,j,k,l) = kkk
            msymqqi(j,i,k,l) = kkk
            msymqqi(i,j,l,k) = kkk
            msymqqi(j,i,l,k) = kkk
            msymqqi(k,l,i,j) = - kkk
            msymqqi(l,k,i,j) = - kkk
            msymqqi(k,l,j,i) = - kkk
            msymqqi(l,k,j,i) = - kkk
          endif

          do iss = 2,min(48,nopsm)

            is = iordresym(iss)
            if( is > 48 ) cycle
            if( iopsymt(is) == 0 ) cycle
            isgmag = iopsymt(is)
            call opsym(is,matopsym)
            do ii = 1,3
              do jj = 1,3
                if( abs(matopsym(jj,ii)) < 0.0001_db ) cycle
                im(ii) = jj * nint( matopsym(jj,ii) )
                exit
              end do
            end do
            ii = abs( im(i) )
            iis = ii / im(i)
            jj = abs( im(j) )
            jjs = jj / im(j)
            kk = abs( im(k) )
            kks = kk / im(k)
            ll = abs( im(l) )
            lls = ll / im(l)
            isg = iis * jjs * kks * lls
            isgi = isg * isgmag

            if( ( isg > 0 .and. msymqq(ii,jj,kk,ll) == - kkk ) .or. ( isg < 0 .and. msymqq(ii,jj,kk,ll) == kkk ) ) &
                                                 impos(kkk) = .true.
            if( ( isgi > 0 .and. msymqqi(ii,jj,kk,ll) == - kkk ).or. ( isgi < 0 .and. msymqqi(ii,jj,kk,ll) == kkk ) ) &
                                                 imposi(kkk) =.true.

            if( .not. impos(kkk) ) then
              msymqq(ii,jj,kk,ll) = isg * kkk
              msymqq(jj,ii,kk,ll) = isg * kkk
              msymqq(ii,jj,ll,kk) = isg * kkk
              msymqq(jj,ii,ll,kk) = isg * kkk
              msymqq(kk,ll,ii,jj) = isg * kkk
              msymqq(ll,kk,ii,jj) = isg * kkk
              msymqq(kk,ll,jj,ii) = isg * kkk
              msymqq(ll,kk,jj,ii) = isg * kkk
            endif

            if( .not. imposi(kkk) .and. isgmag /= 0 ) then
              msymqqi(ii,jj,kk,ll) = isgi * kkk
              msymqqi(jj,ii,kk,ll) = isgi * kkk
              msymqqi(ii,jj,ll,kk) = isgi * kkk
              msymqqi(jj,ii,ll,kk) = isgi * kkk
              msymqqi(kk,ll,ii,jj) = - isgi * kkk
              msymqqi(ll,kk,ii,jj) = - isgi * kkk
              msymqqi(kk,ll,jj,ii) = - isgi * kkk
              msymqqi(ll,kk,jj,ii) = - isgi * kkk
            endif

          end do

        end do
      end do
    end do
  end do

  kkkmax = kkk
  kk = 0
  do kkk = 1,kkkmax
    kk = kk + 1
    if( impos(kkk) ) where( abs(msymqq) == kk ) msymqq = 0
    if( imposi(kkk) ) where( abs(msymqqi) == kk ) msymqqi = 0
    if( impos(kkk) .and. imposi(kkk) ) then
      where( abs(msymqq) > kk ) msymqq = msymqq - abs(msymqq)/msymqq
      where( abs(msymqqi) > kk ) msymqqi = msymqqi - abs(msymqqi) / msymqqi
      kk = kk - 1
    endif
  end do

  n = maxval( abs(msymqq) )
  boucle_i: do i = 1,n
    do j = 1,3
      do k = 1,3
        do l = 1,3
          do m = 1,3
            if( abs(msymqq(j,k,l,m)) == i .or. abs(msymqqi(j,k,l,m)) == i ) cycle boucle_i
          end do
        end do
      end do
    end do
    where( abs(msymqq) > i ) msymqq = msymqq - abs(msymqq) / msymqq
    where( abs(msymqqi) > i ) msymqqi = msymqqi - abs(msymqqi) / msymqqi
  end do boucle_i

  return
end

!*********************************************************************

! Calcul du tenseur dipole-octupole

subroutine tensdo(iopsymt,magnet,msymdo,msymdoi)

  use declarations

  integer im(3), iopsymt(nopsm)
  integer, dimension(3,3,3,3):: msymdo, msymdoi

  logical magnet
  logical, dimension(81):: impos, imposi

  real(kind=db):: matopsym(3,3)

  msymdo(:,:,:,:) = 0; msymdoi(:,:,:,:) = 0
  impos(:) = .false.
  if( magnet ) then
   imposi(:) = .false.
  else
   imposi(:) = .true.
  endif

  kkk = 0
  do i = 1,3
    do j = 1,3
      do k = j,3
        do l = k,3
          if( msymdo(i,j,k,l) /= 0 .or. msymdoi(i,j,k,l) /= 0 )cycle
          kkk = kkk + 1
          msymdo(i,j,k,l) = kkk
          msymdo(i,j,l,k) = kkk
          msymdo(i,k,j,l) = kkk
          msymdo(i,l,j,k) = kkk
          msymdo(i,k,l,j) = kkk
          msymdo(i,l,k,j) = kkk
          msymdoi(i,j,k,l) = kkk
          msymdoi(i,j,l,k) = kkk
          msymdoi(i,k,j,l) = kkk
          msymdoi(i,l,j,k) = kkk
          msymdoi(i,k,l,j) = kkk
          msymdoi(i,l,k,j) = kkk

          do iss = 2,min(48,nopsm)

            is = iordresym(iss)

            if( iopsymt(is) == 0 ) cycle
            isgmag = iopsymt(is)
            call opsym(is,matopsym)
            do ii = 1,3
              do jj = 1,3
                if( abs(matopsym(jj,ii)) < 0.0001_db ) cycle
                im(ii) = jj * nint( matopsym(jj,ii) )
                exit
              end do
            end do
            ii = abs( im(i) )
            iis = ii / im(i)
            jj = abs( im(j) )
            jjs = jj / im(j)
            kk = abs( im(k) )
            kks = kk / im(k)
            ll = abs( im(l) )
            lls = ll / im(l)
            isg = iis * jjs * kks * lls
            isgi = isg * isgmag

            if( ( isg > 0 .and. msymdo(ii,jj,kk,ll) == - kkk ) .or. ( isg < 0 .and. msymdo(ii,jj,kk,ll) == kkk ) ) &
                                                 impos(kkk) = .true.
            if( ( isgi > 0 .and. msymdoi(ii,jj,kk,ll) == - kkk ).or. ( isgi < 0 .and. msymdoi(ii,jj,kk,ll) == kkk ) ) &
                                                 imposi(kkk) =.true.


            if( .not. impos(kkk) ) then
              msymdo(ii,jj,kk,ll) = isg * kkk
              msymdo(ii,jj,ll,kk) = isg * kkk
              msymdo(ii,kk,jj,ll) = isg * kkk
              msymdo(ii,ll,jj,kk) = isg * kkk
              msymdo(ii,kk,ll,jj) = isg * kkk
              msymdo(ii,ll,kk,jj) = isg * kkk
            endif

            if( .not. impos(kkk) .and. isgmag /= 0 ) then
              msymdoi(ii,jj,kk,ll) = isgi * kkk
              msymdoi(ii,jj,ll,kk) = isgi * kkk
              msymdoi(ii,kk,jj,ll) = isgi * kkk
              msymdoi(ii,ll,jj,kk) = isgi * kkk
              msymdoi(ii,kk,ll,jj) = isgi * kkk
              msymdoi(ii,ll,kk,jj) = isgi * kkk
            endif

          end do

        end do
      end do
    end do
  end do

  kkkmax = kkk
  kk = 0
  do kkk = 1,kkkmax
    kk = kk + 1
    if( impos(kkk) ) where( abs(msymdo) == kk ) msymdo = 0
    if( imposi(kkk) ) where( abs(msymdoi) == kk ) msymdoi = 0
    if( impos(kkk) .and. imposi(kkk) ) then
      where( abs(msymdo) > kk ) msymdo = msymdo - abs(msymdo)/msymdo
      where( abs(msymdoi) > kk ) msymdoi = msymdoi - abs(msymdoi) / msymdoi
      kk = kk - 1
    endif
  end do

  n = maxval( abs(msymdo) )
  boucle_i: do i = 1,n
    do j = 1,3
      do k = 1,3
        do l = 1,3
          do m = 1,3
            if( abs(msymdo(j,k,l,m)) == i .or. abs(msymdoi(j,k,l,m)) == i ) cycle boucle_i
          end do
        end do
      end do
    end do
    where( abs(msymdo) > i ) msymdo = msymdo - abs(msymdo) / msymdo
    where( abs(msymdoi) > i ) msymdoi = msymdoi - abs(msymdoi) / msymdoi
  end do boucle_i

  return
end

!*********************************************************************

! Calcul du tenseur E3E3

subroutine tensoo(iopsymt,magnet,msymoo,msymooi,n_oo)

  use declarations
  implicit none

  integer:: he, hhe, hhe_s, hhs, hhs_s, hs, i, i1, i2, i3, ii, iie, &
    iis, is, isg, isgi, iordresym, isgmag, iss, j, j1, j2, j3, je, jhe, jhs, jj, jje, jje_s, jjhe, jjhs, &
    jjs, jjs_s, js, ke, kk, kke, kke_s, kkk, kkkmax, kks, kks_s, ks, n, n_oo
  integer, dimension(3):: im(3)
  integer, dimension(nopsm):: iopsymt
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi

  logical:: Magnet
  logical, dimension(9*n_oo**2):: impos, imposi

  real(kind=db):: matopsym(3,3)

  msymoo(:,:,:,:) = 0; msymooi(:,:,:,:) = 0
  impos(:) = .false.
  if( Magnet ) then
   imposi(:) = .false.
  else
   imposi(:) = .true.
  endif

  kkk = 0
  do ke = 1,3
    do je = 1,3
      do he = 1,3
        jhe = 3 * ( je - 1 ) + he
        do ks = 1,3
          do js = 1,3
            do hs = 1,3
              jhs = 3 * ( js - 1 ) + hs

              if( msymoo(ke,jhe,ks,jhs) /= 0 .or. msymooi(ke,jhe,ks,jhs) /= 0 ) cycle

              kkk = kkk + 1

              do iie = 1,6
                select case(iie)
                  case(1)
                    i1 = ke; i2 = je; i3 = he
                  case(2)
                    i1 = ke; i2 = he; i3 = je
                    if( i2 == i3 ) cycle
                  case(3)
                    i1 = je; i2 = ke; i3 = he
                    if( i1 == i2 ) cycle
                  case(4)
                    i1 = he; i2 = ke; i3 = je
                    if( i1 == i3 ) cycle
                  case(5)
                    i1 = je; i2 = he; i3 = ke
                    if( i2 == i3 ) cycle
                  case(6)
                    i1 = he; i2 = je; i3 = ke
                    if( i1 == i3 ) cycle
                end select
                i = 3 * ( i2 - 1 ) + i3

                do iis = 1,6
                  select case(iie)
                    case(1)
                      j1 = ks; j2 = js; j3 = hs
                    case(2)
                      j1 = ks; j2 = hs; j3 = js
                      if( j2 == j3 ) cycle
                    case(3)
                      j1 = js; j2 = ks; j3 = hs
                      if( j1 == j2 ) cycle
                    case(4)
                      j1 = hs; j2 = ks; j3 = js
                      if( j1 == j3 ) cycle
                    case(5)
                      j1 = js; j2 = hs; j3 = ks
                      if( j2 == j3 ) cycle
                    case(6)
                      j1 = hs; j2 = js; j3 = ks
                      if( j1 == j3 ) cycle
                  end select
                  j = 3 * ( j2 - 1 ) + j3

                  msymoo(i1,i,j1,j) = kkk
                  msymoo(j1,j,i1,i) = kkk

                  if( ke /= ks .or. je /= js .or. he /= hs ) then
                    msymooi(i1,i,j1,j) = kkk
                    msymooi(j1,j,i1,i) = - kkk
                  endif

                end do
              end do  ! fin boucle iie

              do iss = 2,min(48,nopsm)

                is = iordresym(iss)
                if( is > 48 ) cycle
                if( iopsymt(is) == 0 ) cycle
                isgmag = iopsymt(is)
                call opsym(is,matopsym)
                do ii = 1,3
                  do jj = 1,3
                    if( abs(matopsym(jj,ii)) < 0.0001_db ) cycle
                    im(ii) = jj * nint( matopsym(jj,ii) )
                    exit
                  end do
                end do
                kke = abs( im(ke) )
                kke_s =  kke / im(ke)
                jje = abs( im(je) )
                jje_s =  jje / im(je)
                hhe = abs( im(he) )
                hhe_s =  hhe / im(he)
                kks = abs( im(ks) )
                kks_s =  kks / im(ks)
                jjs = abs( im(js) )
                jjs_s =  jjs / im(js)
                hhs = abs( im(hs) )
                hhs_s =  hhs / im(hs)
                isg = kke_s * jje_s * hhe_s * kks_s * jjs_s * hhs_s
                isgi = isg * isgmag

                jjhe = 3 * ( jje - 1 ) + hhe
                jjhs = 3 * ( jjs - 1 ) + hhs

                if( ( isg > 0 .and. msymoo(kke,jjhe,kks,jjhs) == - kkk ) .or. ( isg < 0 .and. &
                          msymoo(kke,jjhe,kks,jjhs) == kkk ) ) impos(kkk) = .true.
                if( ( isgi > 0 .and. msymooi(kke,jjhe,kks,jjhs) == - kkk ).or. ( isgi < 0 .and. &
                          msymooi(kke,jjhe,kks,jjhs) == kkk ) ) imposi(kkk) =.true.

                if( impos(kkk) .and. imposi(kkk) ) cycle

                do iie = 1,6
                  select case(iie)
                    case(1)
                      i1 = kke; i2 = jje; i3 = hhe
                    case(2)
                      i1 = kke; i2 = hhe; i3 = jje
                      if( i2 == i3 ) cycle
                    case(3)
                      i1 = jje; i2 = kke; i3 = hhe
                      if( i1 == i2 ) cycle
                    case(4)
                      i1 = hhe; i2 = kke; i3 = jje
                      if( i1 == i3 ) cycle
                    case(5)
                      i1 = jje; i2 = hhe; i3 = kke
                      if( i2 == i3 ) cycle
                    case(6)
                      i1 = hhe; i2 = jje; i3 = kke
                      if( i1 == i3 ) cycle
                  end select
                  i = 3 * ( i2 - 1 ) + i3

                  do iis = 1,6
                    select case(iie)
                      case(1)
                        j1 = kks; j2 = jjs; j3 = hhs
                      case(2)
                        j1 = kks; j2 = hhs; j3 = jjs
                        if( j2 == j3 ) cycle
                      case(3)
                        j1 = jjs; j2 = kks; j3 = hhs
                        if( j1 == j2 ) cycle
                      case(4)
                        j1 = hhs; j2 = kks; j3 = jjs
                        if( j1 == j3 ) cycle
                      case(5)
                        j1 = jjs; j2 = hhs; j3 = kks
                        if( j2 == j3 ) cycle
                      case(6)
                        j1 = hhs; j2 = jjs; j3 = kks
                        if( j1 == j3 ) cycle
                    end select
                    j = 3 * ( j2 - 1 ) + j3

                    if( .not. impos(kkk) ) then
                      msymoo(i1,i,j1,j) = isg * kkk
                      msymoo(j1,j,i1,i) = isg * kkk
                    endif

                    if( .not. imposi(kkk) .and. isgmag /= 0 .and. (kke /= kks .or. jje /= jjs .or. hhe /= hhs) ) then
                      msymooi(i1,i,j1,j) = isgi * kkk
                      msymooi(j1,j,i1,i) = - isgi * kkk
                    endif

                  end do
                end do  ! fin boucle iie

              end do  ! boucle_iss

            end do
          end do
        end do
      end do
    end do
  end do

  kkkmax = kkk
  kk = 0
  do kkk = 1,kkkmax
    kk = kk + 1
    if( impos(kkk) ) where( abs(msymoo) == kk ) msymoo = 0
    if( imposi(kkk) ) where( abs(msymooi) == kk ) msymooi = 0
    if( impos(kkk) .and. imposi(kkk) ) then
      where( abs(msymoo) > kk ) msymoo = msymoo - abs(msymoo)/msymoo
      where( abs(msymooi) > kk ) msymooi = msymooi - abs(msymooi) / msymooi
      kk = kk - 1
    endif
  end do

  n = maxval( abs(msymoo) )
  boucle_i: do i = 1,n
    do ke = 1,3
      do jhe = 1,n_oo
        do ks = 1,3
          do jhs = 1,n_oo
            if( abs(msymoo(ke,jhe,ks,jhs)) == i .or. abs(msymooi(ke,jhe,ks,jhs)) == i ) cycle boucle_i
          end do
        end do
      end do
    end do
    where( abs(msymoo) > i ) msymoo = msymoo - abs(msymoo) / msymoo
    where( abs(msymooi) > i ) msymooi = msymooi - abs(msymooi) / msymooi
  end do boucle_i

  return
end

!***********************************************************************

subroutine prodvec(u,v,w)

  use declarations
  real(kind=db), dimension(3):: u, v, w

  u(1) = v(2) * w(3) - v(3) * w(2)
  u(2) = v(3) * w(1) - v(1) * w(3)
  u(3) = v(1) * w(2) - v(2) * w(1)

  return
end

!***********************************************************************

subroutine prodvec_cp(u,v,w)

  use declarations
  complex(kind=db), dimension(3):: u, v, w

  u(1) = v(2) * w(3) - v(3) * w(2)
  u(2) = v(3) * w(1) - v(1) * w(3)
  u(3) = v(1) * w(2) - v(2) * w(1)

  return
end

!***********************************************************************

! Calcul direct de l'inverse d'une matrice 3 x 3

subroutine invermat(a,b)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db), dimension(3,3):: a, b, cofac
  real(kind=db), dimension(2,2):: c

  n = 3

  do i = 1,n
    do j = 1,n

      do i1 = 1,n
        do j1 = 1,n
          if( i1 == i .or. j1 == j) cycle
          if( i1 < i ) then
            i2 = i1
          else
            i2 = i1 - 1
          endif
          if( j1 < j ) then
            j2 = j1
          else
            j2 = j1 - 1
          endif
          c(i2,j2) = a(i1,j1)
        end do
      end do
      cofac(i,j) = ( c(1,1)*c(2,2) - c(1,2)*c(2,1) ) * (-1)**(i+j)

    end do
  end do

  det = sum( a(1,:) * cofac(1,:) )

  b = transpose( cofac ) / det

  return
end

!*********************************************************************

subroutine cal_iaeqrmt(iaabs,iaproto,iapot,igreq,igroup,n_atom_proto,natomp,neqm,ngreq,nonexc_g,self_nonexc)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(natomp):: igroup, iaproto
  integer, dimension(0:n_atom_proto):: iapot, ngreq
  integer, dimension(0:n_atom_proto,neqm):: igreq

  logical nonexc_g, self_nonexc

  iapot(:) = 0
! Choix des atomes prototypiques dont le potentiel sera calcule
  iapot(0) = iaabs
  do ipr = 1,n_atom_proto
    iapot(ipr) = 0
    boucle_ia: do ia = 1,natomp
      if( ia == iaabs .and. .not. self_nonexc ) cycle
      do igr = 1,ngreq(ipr)
        if( igroup(ia) == igreq(ipr,igr) ) then
          iapot(ipr) = ia
          exit boucle_ia
        endif
      end do
    end do boucle_ia
  end do

  do ia = 1,natomp
    boucle_i: do ipr = 1,n_atom_proto
      do igr = 1,ngreq(ipr)
        if( igroup(ia) == igreq(ipr,igr) ) exit boucle_i
      end do
    end do boucle_i
    iaproto(ia) = ipr
  end do

  if( .not. nonexc_g ) iaproto(iaabs) = 0

  return
end

!*********************************************************************

! Calcul de l'absorption avant seuil

  complex(kind=db) function f_cal(Doping,Eseuil,Full_self_abs,icheck, itypepr,n_atom_proto,nbseuil,ngreq,ntype,numat,Z_abs, &
          self_abs,Taux_ipr,Volume_maille,xsect_file)

  use declarations
  implicit none

  character(len=4):: elemv
  character(len=2):: Chemical_symbol
  character(len=132):: xsect_file

  integer icheck, ipr, n_atom_proto, nbseuil, ntype, Z, Z_abs
  integer, dimension(0:ntype):: numat
  integer, dimension(0:n_atom_proto):: itypepr, ngreq

  logical:: Doping, Full_self_abs, self_abs

  real(kind=sg):: getf0

  real(kind=db):: Conv_mbarn_nelec, Ea, fpp_avantseuil, fp_avantseuil, Volume_maille
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(n_atom_proto):: f0, fpa, fppa, Taux_ipr

  fpp_avantseuil = 0._db
  fp_avantseuil = 0._db
  fppa(:) = 0._db
  fpa(:) = 0._db

  if( Eseuil(nbseuil) > eps10 ) then ! Sinon on est en optic

    do ipr = 1,n_atom_proto
      if( Doping .and. ipr == n_atom_proto ) cycle
      Z = numat( itypepr( ipr ) )
      if( Z == 0 ) cycle
      Ea = Eseuil(nbseuil) - 1._db
      Ea = max( Ea,  0.5_db )
      call fprime(Z,Ea,fppa(ipr),fpa(ipr),xsect_file)
      fpp_avantseuil = fpp_avantseuil + Taux_ipr(ipr) * ngreq(ipr) * fppa(ipr)

      if( Z /= Z_abs ) fp_avantseuil = fp_avantseuil + Taux_ipr(ipr) * ngreq(ipr) * fpa(ipr)

      elemv = ' '
      elemv(1:2) = Chemical_Symbol(Z)
      f0(ipr) = getf0(elemv,0._sg)
      fp_avantseuil = fp_avantseuil + Taux_ipr(ipr) * ngreq(ipr) * f0(ipr)
    end do

! Conversion en Megabarn (= 10^-18 cm2 = 10^-22 m2 = 10^-2 A2)
    fpp_avantseuil = fpp_avantseuil / ( conv_mbarn_nelec(Ea) * pi )
    fp_avantseuil = fp_avantseuil / ( conv_mbarn_nelec(Ea) * pi )
! Conversion en coefficient d'absorption lineaire en micrometre^-1
    fpp_avantseuil = 100 * fpp_avantseuil /(Volume_maille * bohr**3)
    fp_avantseuil = 100 * fp_avantseuil /(Volume_maille * bohr**3)

    if( icheck > 0 ) then
      write(3,310)
      do ipr = 1,n_atom_proto
        if( Doping .and. ipr == n_atom_proto ) cycle
        write(3,320) ipr, numat( itypepr(ipr) ), f0(ipr), fpa(ipr), fppa(ipr)
      end do
      write(3,330) fpp_avantseuil, fp_avantseuil
    endif

  endif

  if( self_abs ) then
    f_cal = cmplx( fp_avantseuil, fpp_avantseuil )
  elseif( Full_self_abs ) then
    f_cal = cmplx( fp_avantseuil, - fpp_avantseuil )
  else
    f_cal = (0._db, 0._db)
  endif

  return

  310 format(/' Absorption before the edge :',/ &
          ' Site  Z      f0           fp          fpp    (per atom in nbr of electron)')
  320 format(2i4,1p,3e13.5)
  330 format(/' Complex linear absorption coefficient :',1p,2e13.5,' micrometer^(-1)')
end

!*********************************************************************

! Sousprogramme effectuant certaines preparations pour le DAFS

subroutine prepdafs(Angle_or,angpoldafs,angxyz,Axe_atom_gr,axyz,Base_spin,Bormann,Dafs_bio,Eseuil,f_no_res,hkl_dafs, &
            icheck,igreq,iprabs,isigpi,itabs,itypepr,lvval,Magnetic,Mat_or,mpirank,n_atom_proto,natomsym,nbseuil, &
            neqm,ngreq,ngrm,ngroup,ngroup_m,ngroup_taux,ngroup_temp,nlat,nlatm,nphi_dafs,nphim,npldafs, &
            nrato,nrm,nspin,ntype,numat,Orthmat,Orthmati,phdafs,phdf0t,phdt,poldafse,poldafsem,poldafss,poldafssm, &
            popatm,posn,psival,rato,Rot_int,Taux,Taux_oc,temp,Temp_coef,Temperature,Vec_orig,vecdafse,vecdafsem, &
            vecdafss,vecdafssm,xsect_file)

  use declarations
  implicit none

  integer:: i, icheck, igr, ip, ipl, ipr, iprabs, istop, it, itabs, &
    iwrite, j, jgr, kgr, mpirank, n, n_atom_proto, natomsym, &
    nbseuil, neqm, ngrm, ngroup, ngroup_m, ngroup_taux, ngroup_temp, nlatm, nphim, npldafs, nrm, nspin, ntype

  character(len=4) elemv, mot4
  character(len=2) Chemical_symbol
  character(len=132):: xsect_file

  complex(kind=db):: cfac, ph, ph_cjg
  complex(kind=db), dimension(3):: vec_a, vec_b, pe, ps
  complex(kind=db), dimension(npldafs):: phdabs, v_vec_a, v_vec_b
  complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdt
  complex(kind=db), dimension(natomsym,npldafs):: phdafs
  complex(kind=db), dimension(n_atom_proto,npldafs):: phd, phd_f0, phd_f0_cjg, phd_fan_cjg, phd_fan, phd_fmag
  complex(kind=db), dimension(3,npldafs):: poldafsem, poldafssm
  complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss
  complex(kind=db), dimension(:,:), allocatable:: Bragg
  complex(kind=db), dimension(:,:,:), allocatable:: phd_fmo,phd_fms

  integer, dimension(0:ntype):: nlat, nrato, numat
  integer, dimension(0:ntype,nlatm):: lvval
  integer, dimension(3,npldafs):: hkl_dafs
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(npldafs,2):: isigpi
  integer, dimension(0:n_atom_proto):: itypepr, ngreq
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(neqm):: igr_tem

  logical:: Base_spin, Bormann, Dafs_bio, Debye, Magnetic, Taux, Temperature

  real(kind=sg):: getf0, s

  real(kind=db):: Angle, arg, cos_pe, cos_ps, cosb, cosp, Deb, Delta_2, deltak_A, detmat, dhkl, dp, dpdeg, DW, &
    Emc2, fac, fmo, fms, konde, pp, psi, qkn, rad, rap_lsur2s, sin_pe, sin_ps, sinb, sinp, Temp, Tempt, Thetabragg, vol, wn, wsn

  real(kind=db), dimension(2):: f_no_res
  real(kind=db), dimension(3):: ang, angxyz, axyz, cosdir, hklred, &
     qi, qj, qk, v, Vec_orig, vpie, vpe, vps, vpis, vsig, vx, vy, vz, w, we, ws, wx, wy, wz
  real(kind=db), dimension(npldafs):: Angle_or, deltak
  real(kind=db), dimension(ngroup_taux):: Taux_oc
  real(kind=db), dimension(ngroup_temp):: Temp_coef
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(0:3):: det
  real(kind=db), dimension(3,3):: Mat, Mat_or, Mat_ori, Orthmat, Orthmati, Rot_int
  real(kind=db), dimension(3,ngroup):: posn
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_g, Axe_atom_gr
  real(kind=db), dimension(3,npldafs):: angpoldafs, vecdafsem, vecdafssm
  real(kind=db), dimension(3,npldafs,nphim):: vecdafse, vecdafss
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(n_atom_proto):: fp, fpp
  real(kind=db), dimension(n_atom_proto,npldafs):: f_ms, f_mo, f0
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm

  if( icheck > 0 ) write(3,110)

! lambda = 2 * pi / k = 2 * d * sintheta
! En S.I. vecond = k = E*alfa*4*pi*epsilon0 / (e*e)
! En ua et rydb : k = 0.5 * alfa * E
  konde = 0.5_db * alfa_sf * eseuil(nbseuil)   ! exact que pour nbseuil = 1

  istop = 0
  cosdir(:) = cos( angxyz(:) * pi / 180 )

  if( Bormann ) then
    poldafse(:,:,:) = ( 0._db, 0._db )
    poldafss(:,:,:) = ( 0._db, 0._db )
    vecdafse(:,:,:) = 0._db
    vecdafss(:,:,:) = 0._db
  endif

  do ipl = 1,npldafs

    if( hkl_dafs(1,ipl) == 0 .and. hkl_dafs(2,ipl) == 0 .and. hkl_dafs(3,ipl) == 0 ) then

      dhkl = 0._db
      thetabragg = 0._db

    else

      hklred(:) = hkl_dafs(:,ipl) / axyz(:)
      do i = 0,3
        do j = 1,3
          mat(j,j) = 1._db
        end do
        mat(1,2) = cosdir(3)
        mat(1,3) = cosdir(2)
        mat(2,3) = cosdir(1)
        mat(2,1) = mat(1,2)
        mat(3,1) = mat(1,3)
        mat(3,2) = mat(2,3)
        if( i > 0 ) mat(:,i) = hklred(:)
        det(i) = detmat(mat)
      end do
! Distance interplan
      dhkl = sqrt( det(0) / sum( hklred(1:3) * det(1:3) ) )
      fac = abs( pi / ( konde * dhkl ) )

      if( fac > 1._db .and. mpirank == 0 ) then
        if( istop == 0 ) call write_error
        do ipr = 3,9,3
          if( ipr == 3 .and. icheck == 0 ) cycle
          write(ipr,120) ipl, hkl_dafs(:,ipl)
        end do
        istop = 1
        cycle
      endif

      thetabragg = asin( pi / ( konde * dhkl ) )

    endif

! Quand angpoldafs(3,ipl) = 10000, polarisation et veconde sont imposes
    if( angpoldafs(3,ipl) > 9999._db ) then

      poldafse(:,ipl,1) = poldafsem(:,ipl)
      poldafss(:,ipl,1) = poldafssm(:,ipl)
      vecdafse(:,ipl,1) = vecdafsem(:,ipl)
      vecdafss(:,ipl,1) = vecdafssm(:,ipl)

! On passe en base interne (orthonormee)
      v(:) = real( poldafse(:,ipl,1), db )
      if( sum( v(:)**2 )  > eps6 ) then
        call trvec(mpirank,Orthmat,v,w)
        poldafse(:,ipl,1) = cmplx( w(:), 0._db )
      endif

      v(:) = real( poldafss(:,ipl,1), db )
      if( sum( v(:)**2 )  > eps6 ) then
        call trvec(mpirank,Orthmat,v,w)
        poldafss(:,ipl,1) = cmplx( w(:), 0._db )
      endif

      v(:) = vecdafse(:,ipl,1)
      if( sum( v(:)**2 )  > eps6 ) then
        call trvec(mpirank,Orthmat,v,w)
        vecdafse(:,ipl,1) = w(:)
      endif

      v(:) = vecdafss(:,ipl,1)
      if( sum( v(:)**2 )  > eps6 ) then
        call trvec(mpirank,Orthmat,v,w)
        vecdafss(:,ipl,1) = w(:)
      endif

    elseif( Bormann ) then

      select case(ipl)
        case(1,10,19,28)
          poldafse(1,ipl,1) = ( 1._db, 0._db )
          poldafss(1,ipl,1) = ( 1._db, 0._db )
        case(2,11,20,29)
          poldafse(1,ipl,1) = ( 1._db, 0._db )
          poldafss(2,ipl,1) = ( 1._db, 0._db )
        case(3,12,21,30)
          poldafse(1,ipl,1) = ( 1._db, 0._db )
          poldafss(3,ipl,1) = ( 1._db, 0._db )
        case(4,13,22,31)
          poldafse(2,ipl,1) = ( 1._db, 0._db )
          poldafss(1,ipl,1) = ( 1._db, 0._db )
        case(5,14,23,32)
          poldafse(2,ipl,1) = ( 1._db, 0._db )
          poldafss(2,ipl,1) = ( 1._db, 0._db )
        case(6,15,24,33)
          poldafse(2,ipl,1) = ( 1._db, 0._db )
          poldafss(3,ipl,1) = ( 1._db, 0._db )
        case(7,16,25,34)
          poldafse(3,ipl,1) = ( 1._db, 0._db )
          poldafss(1,ipl,1) = ( 1._db, 0._db )
        case(8,17,26,35)
          poldafse(3,ipl,1) = ( 1._db, 0._db )
          poldafss(2,ipl,1) = ( 1._db, 0._db )
        case(9,18,27,36)
          poldafse(3,ipl,1) = ( 1._db, 0._db )
          poldafss(3,ipl,1) = ( 1._db, 0._db )
      end select

    else

      sinb = sin( thetabragg )
      cosb = cos( thetabragg )

! Vecteur diffraction dans la base orthogonale interne
      vx(:) = orthmat(:,1)
      vy(:) = orthmat(:,2)
      vz(:) = orthmat(:,3)

! wx, wy, wz : base du reseau reciproque
      call prodvec(wx,vy,vz)

      vol = sum( wx(:) * vx(:) )
      wx(:) = wx(:) / vol
      call prodvec(wy,vz,vx)
      wy(:) = wy(:) / vol
      call prodvec(wz,vx,vy)
      wz(:) = wz(:) / vol

      qk(:) = hkl_dafs(1,ipl) * wx(:) + hkl_dafs(2,ipl) * wy(:) + hkl_dafs(3,ipl) * wz(:)
      qkn = sqrt( sum( qk(:)**2 ) )
      if( abs(qkn) > eps10 ) then
        qk(:) = qk(:) / qkn
      else
! Cas speculaire
        qk(1:2) = 0._db; qk(3) = 1._db
      endif
! On prend l'origine de l'azimut selon le plan (Q,Vec_orig)
! si pas possible, on prend selon (Q,Ox) ou (Q,Oz).
      do i = 1,3
        select case(i)
          case(1)
            v(:) = Vec_orig(:)
          case(2)
            v(1) = 1._db; v(2:3) = 0._db
          case(3)
            v(1:2) = 0._db; v(3) = 1._db
        end select
        v = matmul( orthmat, v )
        call prodvec(w,qk,v)
        wn = sqrt( dot_product(w,w) )
        if( wn > eps10 ) exit
      end do
      w(:) = w(:) / wn
      call prodvec(qi,w,qk)
      call prodvec(qj,qk,qi)

      if( icheck > 1 ) then
        write(3,130) hkl_dafs(:,ipl)
        write(3,140) ( wx(i)/bohr, wy(i)/bohr, wz(i)/bohr, i = 1,3 )
        write(3,150) ( qi(i), qj(i), qk(i), i = 1,3 )
      endif

      if( Dafs_bio ) then

        if( ipl == 1 ) then
          call invermat(Mat_or,Mat_ori)
          Mat_ori = bohr * Transpose( Mat_ori )
        endif

! Matrice de rotation inverse
        Angle = Angle_or(ipl) * pi / 180
        Mat(1,1) = 1._db; Mat(1,2:3) = 0._db
        Mat(2,1) = 0._db;  Mat(2,2) = cos( Angle )
                                           Mat(2,3) = sin( Angle )
        Mat(3,1) = 0._db; Mat(3,1) = - Mat(2,3); Mat(3,3) = Mat(2,2)

        Mat = Matmul( Mat_ori, Mat )
        Mat = Matmul( Orthmat, Mat )
! La polarisation est selon Ox
        v(:) = Mat(:,1)
! Le vecteur d'onde est selon Oz
        we(:) =  Mat(:,3)

        ws(:) = qk(:) - we(:)
        wsn = sqrt( sum( ws(:)**2 ) )
        ws(:) = ws(:) / wsn

        call prodvec(vsig,qk,we)
        wsn = sqrt( sum( vsig(:)**2 ) )
        vsig(:) = vsig(:) / wsn

        call prodvec(vpis,vsig,ws)

        poldafse(:,ipl,1) = cmplx( v(:), 0._db )
        vecdafse(:,ipl,1) = we(:)
        poldafss(:,ipl,1) = cmplx( vsig(:), 0._db )
        vecdafss(:,ipl,1) = ws(:)

        poldafse(:,ipl,2) = cmplx( v(:), 0._db )
        vecdafse(:,ipl,2) = we(:)
        poldafss(:,ipl,2) = cmplx( vpis(:), 0._db )
        vecdafss(:,ipl,2) = ws(:)

! Reflection (-Q)
        poldafse(:,ipl,3) = cmplx( v(:), 0._db )
        vecdafse(:,ipl,3) = - we(:)
        poldafss(:,ipl,3) = cmplx( vsig(:), 0._db )
        vecdafss(:,ipl,3) = - ws(:)

        poldafse(:,ipl,4) = cmplx( v(:), 0._db )
        vecdafse(:,ipl,4) = - we(:)
        poldafss(:,ipl,4) = cmplx( - vpis(:), 0._db )
        vecdafss(:,ipl,4) = - ws(:)

      else

        dp = 2 * pi / nphi_dafs(ipl)

        if( angpoldafs(1,ipl) > -9999._db ) then
          cos_pe = cos( angpoldafs(1,ipl) )
          sin_pe = sin( angpoldafs(1,ipl) )
        endif
        if( angpoldafs(2,ipl) > -9999._db ) then
          cos_ps = cos( angpoldafs(2,ipl) )
          sin_ps = sin( angpoldafs(2,ipl) )
        endif

        do ip = 1,nphi_dafs(ipl)

          if( angpoldafs(1,ipl) < -9999._db ) then
            Angle = ( ip - 1 ) * dp
            cos_pe = cos( angle )
            sin_pe = sin( angle )
          elseif( angpoldafs(2,ipl) < -9999._db ) then
            Angle = ( ip - 1 ) * dp
            cos_ps = cos( angle )
            sin_ps = sin( angle )
          endif
          if( angpoldafs(3,ipl) < -9999._db ) then
            psi = ( ip - 1 ) * dp
          else
            psi = angpoldafs(3,ipl)
          endif

          if( ip == 1 .or. angpoldafs(3,ipl) < -9999._db ) then
            sinp = sin( psi )
            cosp = cos( psi )
            v(1) = cosb * cosp
            v(2) = cosb * sinp
            v(3) = - sinb
            we(:) = v(1) * qi(:) - v(2) * qj(:) + v(3) * qk(:)
            ws(:) = v(1) * qi(:) - v(2) * qj(:) - v(3) * qk(:)
            vsig(:) = sinp * qi(:) + cosp * qj(:)
            call prodvec(vpie,we,vsig)
            call prodvec(vpis,ws,vsig)
          endif

          vpe(:) = cos_pe * vsig(:) + sin_pe * vpie(:)
          vps(:) = cos_ps * vsig(:) + sin_ps * vpis(:)

! En diffraction tout est defini avec le complexe conjuge
! donc les polarisations circulaires aussi
          select case( isigpi(ipl,1) )
            case(3)
              poldafse(:,ipl,ip) = cmplx( vsig(:),-vpie(:), db) / sqrt(2._db)
            case(4)
              poldafse(:,ipl,ip) = cmplx( vsig(:), vpie(:), db) / sqrt(2._db)
            case default
              poldafse(:,ipl,ip) = cmplx( vpe(:), 0._db, db)
          end select
          vecdafse(:,ipl,ip) = we(:)

          select case( isigpi(ipl,2) )
            case(3)
              poldafss(:,ipl,ip) = cmplx( vsig(:),-vpis(:), db) / sqrt(2._db)
            case(4)
              poldafss(:,ipl,ip) = cmplx( vsig(:), vpis(:), db) / sqrt(2._db)
            case default
              poldafss(:,ipl,ip) = cmplx( vps(:), 0._db, db)
          end select
          vecdafss(:,ipl,ip) = ws(:)

        end do

      endif
    endif

    if( hkl_dafs(1,ipl) == 0 .and. hkl_dafs(2,ipl) == 0 .and. hkl_dafs(3,ipl) == 0 ) then

      if( Bormann ) then

        deltak(ipl) = 0._db

      else

        v(:) = vecdafse(:,ipl,1)
        w(:) = vecdafss(:,ipl,1)
        if( mpirank == 0 .and. ( sum(abs(v(:))) < eps10 .or. sum(abs(w(:))) < eps10 ) ) then
          call write_error
          do ipr = 3,9,3
            write(ipr,160)
          end do
          stop
        endif

! ThetaBragg is here the scattering angle
        pp = dot_product(v,w)
        if( abs( pp - 1._db ) < eps10 ) then
          Thetabragg = 0._db
        else
          Thetabragg = 0.5_db * acos( pp )
        endif
        deltak(ipl) = 2 * sin( thetabragg ) * konde
      endif

    else

      deltak(ipl) = 2 * pi / dhkl

    endif

! Facteur de structure atomique
    s = deltak(ipl) / ( 4 * pi )      ! = sin(theta) / lambda
    s = s / bohr                 ! on le veut en angstroem - 1

    do ipr = 1,n_atom_proto

      it = itypepr( ipr )
      elemv = ' '
      elemv(1:2) = Chemical_Symbol(numat(it))
      f0(ipr,ipl) = getf0(elemv,s)

! Diffusion magnetique non resonante
! Le moment orbital est pris parralele a l'axe de spin...
      if( Magnetic ) then
        call get_fmag(deltak(ipl),fmo,fms,ipr,it,lvval,n_atom_proto,nlat(it),nlatm,nrato(it),nrm,nspin,ntype, &
              popatm,psival,rato)
! Emc2 = ( 9.1093897 * 2.99792458**2 / 1.602) * 10000 / Rydb
        Emc2 = 37561.88242_db
        f_ms(ipr,ipl) = - ( eseuil(nbseuil) / Emc2 ) * f_no_res(1) * fms
        if( f_no_res(2) < -99._db ) then
          f_mo(ipr,ipl) = - ( eseuil(nbseuil) / Emc2 ) * f_no_res(1) * fmo
        else
! Facteur 2 car f_ms correspond a S et non a 2S.
          f_mo(ipr,ipl) = 2 * f_no_res(2) * f_ms(ipr,ipl)
        endif
      endif
    end do

    if( icheck > 0 ) then
      if( ipl == 1 ) write(3,170)
      rad = 180 / pi
      do i = 1,3
        if( abs( angpoldafs(i,ipl) ) > 9999._db ) then
          Ang(i) = 0._db
        else
          Ang(i) = angpoldafs(i,ipl)
        endif
      end do
      pe(:) = poldafse(:,ipl,1)
      ps(:) = poldafss(:,ipl,1)
      we(:) = vecdafse(:,ipl,1)
      ws(:) = vecdafss(:,ipl,1)
      pe = matmul( orthmati, pe )
      ps = matmul( orthmati, ps )
      we = matmul( orthmati, we )
      ws = matmul( orthmati, ws )
      write(3,180) hkl_dafs(:,ipl), Thetabragg * rad, dhkl * bohr, Ang(1:3)*rad, ( Real(pe(i)), aimag(pe(i)), we(i), &
                 Real(ps(i)), aimag(ps(i)), ws(i), i = 1,3)
    endif

  end do  ! fin boucle ipl

  if( istop == 1 ) stop

  if( Magnetic ) then
    allocate( phd_fmo(n_atom_proto,ngrm,npldafs) )
    allocate( phd_fms(n_atom_proto,ngrm,npldafs) )
    phd_fms(:,:,:) = (0._db,0._db)
    phd_fmo(:,:,:) = (0._db,0._db)
  endif

  if( Temp > 0.00001_db ) then
    Debye = .true.
    Tempt = max (1._db, temp)      ! le cas defaut: T = 1 K
  else
    Debye = .false.
  endif

  if( Temperature ) write(3,'(/A)') ' ( h, k, l)   Z   Debye-Waller Attenuation'

  do ipr = 1,n_atom_proto

! Calcul des termes de Bragg
    it = itypepr( ipr )
    if( numat( it ) == 0 ) cycle
    igr_tem(:) = igreq(ipr,:)

    allocate( Bragg(ngreq(ipr),npldafs) )

    do ipl = 1,npldafs

      if( Debye ) then
        deltak_A = deltak(ipl) / bohr   ! on le veut en angstroem - 1
        Deb = DW(deltak_A,numat(it),tempt)
        write(3,165) hkl_dafs(1:3,ipl), numat(it), Deb
      elseif( Temperature ) then
! Delta_2 = ( Sin(Theta_Bragg)/Lambda )^2
! Temp_coef = 8*pi^2 * <u>^2 est en Angtroem^2
        Delta_2 = ( deltak(ipl) / ( quatre_pi * bohr ) )**2
      endif

      do igr = 1,ngreq(ipr)
        jgr = igr_tem(igr)
        arg = deux_pi * sum( posn(:,jgr) * hkl_dafs(:,ipl) )
        Bragg(igr,ipl) = cmplx( cos(arg), sin(arg), db )
        if( Taux ) Bragg(igr,ipl) = Bragg(igr,ipl) * Taux_oc(jgr)
        if( Temperature ) Deb = exp( - Temp_coef(igr) * Delta_2 )
        if( Debye .or. Temperature) Bragg(igr,ipl) = Bragg(igr,ipl) * Deb
        if( ipr == iprabs ) phdafs(igr,ipl) = Bragg(igr,ipl)
      end do
      phd(ipr,ipl) = sum( Bragg(1:ngreq(ipr),ipl) )
      if( ipr == iprabs ) phdabs(ipl) = phd(ipr,ipl)

    end do

    n = numat( itypepr( ipr ) )
    if( n == 0 ) cycle
! On ne tient pas compte de l'anomale de l'atome absorbeur puisqu'il est
! calcule dans le programme !
    if( n /= numat(itabs) ) then
      call fprime(n,eseuil(nbseuil),fpp(ipr),fp(ipr),xsect_file)
    else
      fp(:) = 0._db; fpp(ipr) = 0._db
    endif

    do ipl = 1,npldafs

      phd_f0(ipr,ipl) = phd(ipr,ipl) * f0(ipr,ipl)
      phd_fan(ipr,ipl) = phd(ipr,ipl) * cmplx(fp(ipr), fpp(ipr), db)

      if( Dafs_bio ) then
        phd_f0_cjg(ipr,ipl) = conjg( phd(ipr,ipl) ) * f0(ipr,ipl)
        phd_fan_cjg(ipr,ipl) = conjg( phd(ipr,ipl) ) * cmplx(fp(ipr), fpp(ipr), db)
      endif

      if( Magnetic ) then

        if( abs(f_ms(ipr,ipl)) > eps10 ) then
          kgr = igreq(ipr,1)
          do igr = 1,ngreq(ipr)
            jgr = igreq(ipr,igr)
! Les signes moins sont pour transformer vers la convention cristallo (Green moins), les formules de Gibbs.
! Ces donnees sont utilisees dans la routine convolution a un endroit ou on est en Green -.
! Bragg est deja en Green -.
            phd_fms(ipr,igr,ipl) = - img * Bragg(igr,ipl) * f_ms(ipr,ipl)
            phd_fmo(ipr,igr,ipl) = - img * Bragg(igr,ipl) * f_mo(ipr,ipl)
          end do
        endif

      endif

    end do

    deallocate( Bragg )

  end do  ! fin boucle ipr

  if( icheck > 1 .and. magnetic ) then
    write(3,'(/A)') ' ipr ipl      f_ms         f_mo      phd_fms(ipr,ipl,1..ngr)'
    do ipr = 1,n_atom_proto
      do ipl = 1,npldafs
        if( abs(f_ms(ipr,ipl)) < eps10 ) cycle
        write(3,'(2i4,1p,50e13.5)') ipr, ipl, f_ms(ipr,ipl), f_mo(ipr,ipl), phd_fms(ipr,1:ngreq(ipr),ipl)
      end do
    end do
  endif

! Multplication par epsilon_e * espilon_s
  do ipl = npldafs,1,-1
    ph = sum( phd_f0(1:n_atom_proto,ipl) ) + sum( phd_fan(1:n_atom_proto,ipl) )
    if( Dafs_bio ) ph_cjg = sum( phd_f0_cjg(1:n_atom_proto,ipl) ) + sum( phd_fan_cjg(1:n_atom_proto,ipl) )
    do ip = nphi_dafs(ipl),1,-1
      pe(:) = poldafse(:,ipl,ip)
      ps(:) = poldafss(:,ipl,ip)
      cfac = sum( conjg( ps(:) ) * pe(:) )
      if( Dafs_bio .and. ip > 2 ) then
        phdf0t(ipl,ip) = cfac * ph_cjg
        phdt(ipl,ip) = cfac * conjg( phdabs(ipl) )
      else
        phdf0t(ipl,ip) = cfac * ph
        phdt(ipl,ip) = cfac * phdabs(ipl)
      endif
    end do
    phd_f0(1:n_atom_proto,ipl) = cfac * phd_f0(1:n_atom_proto,ipl)
    phd_fan(1:n_atom_proto,ipl) = cfac * phd_fan(1:n_atom_proto,ipl)
  end do

! Voir M. Blume and Doon Gibbs, PRB, 37, 1779 (1988)
! Formule modifiee pour vect A'' car erreur trouve dans l'article
  phd_fmag(:,:) = (0._db,0._db)
  if( Magnetic ) then
    Axe_atom_g = matmul( orthmat, Axe_atom_gr )
    do ipl = 1,npldafs
      do ip = nphi_dafs(ipl),1,-1
        pe(:) = poldafse(:,ipl,ip)
        ps(:) = poldafss(:,ipl,ip)
        we(:) = vecdafse(:,ipl,ip)
        ws(:) = vecdafss(:,ipl,ip)
! Approximation pour le moment orbital pris oriente parallelement a lui.
! Voir aussi appendix de G. T. Trammell, PR 92, 1387 (1953)
        call get_vec_b(pe,ps,vec_b,we,ws)
        call get_vec_a(pe,ps,vec_a,we,ws)
        do ipr = 1,n_atom_proto

          do igr = 1,ngreq(ipr)
            jgr = igreq(ipr,igr)
            v(:) = Axe_atom_g(:,jgr)
            v_vec_b(ipl) = dot_product(v,vec_b)
            v_vec_a(ipl) = dot_product(v,vec_a)
            if( icheck > 1 ) then
              write(3,182) ipl
              do i = 1,3
                write(3,185) v(i), vec_b(i), vec_a(i)
              end do
              write(3,187) v_vec_b(ipl), v_vec_a(ipl)
            endif

            if( ip == 1 ) &
              phd_fmag(ipr,ipl) = phd_fmag(ipr,ipl) + v_vec_a(ipl) * phd_fmo(ipr,igr,ipl) + v_vec_b(ipl) * phd_fms(ipr,igr,ipl)

            phdf0t(ipl,ip) = phdf0t(ipl,ip) + v_vec_a(ipl) * phd_fmo(ipr,igr,ipl) + v_vec_b(ipl) * phd_fms(ipr,igr,ipl)

          end do
        end do
      end do

    end do

  endif

  if( magnetic ) then
    deallocate( phd_fmo )
    deallocate( phd_fms )
  endif

  if( icheck > 0 ) then

    if( icheck > 1 ) then
      do ipl = 1,npldafs
        if( nphi_dafs(ipl) == 1 ) cycle
        dpdeg = 360._db / nphi_dafs(ipl)
        if( Dafs_bio ) then
          write(3,200)
        else
          write(3,210) hkl_dafs(:,ipl), isigpi(ipl,:)
        endif
        do ip = 1,nphi_dafs(ipl)
          pe(:) = poldafse(:,ipl,ip)
          ps(:) = poldafss(:,ipl,ip)
          we(:) = vecdafse(:,ipl,ip)
          ws(:) = vecdafss(:,ipl,ip)

          pe = matmul( orthmati, pe )
          ps = matmul( orthmati, ps )
          we = matmul( orthmati, we )
          ws = matmul( orthmati, ws )
          if( Dafs_bio ) then
            if( ip == 1 .or. ip == 3 ) then
              mot4 = ' sig'
            else
              mot4 = '  pi'
            endif
            if( ip == 1 .or. ip == 2 ) then
              write(3,215) hkl_dafs(:,ipl), mot4, real(pe(:)), real(ps(:)), we(:), ws(:)
            else
              write(3,215) -hkl_dafs(:,ipl), mot4, real(pe(:)), real(ps(:)), we(:), ws(:)
            endif
          else
            write(3,220)  ( ip - 1 ) * dpdeg, real(pe(:)), aimag(pe(:)), real(ps(:)), aimag(ps(:)), we(:), ws(:)
          endif
        end do
      end do
    endif

    write(3,230)
    do ipl = 1,npldafs
      write(3,240) hkl_dafs(:,ipl), phdafs(1:natomsym,ipl)
    end do

    if( Magnetic ) then
      write(3,245) f_no_res(1)
      ipl = 1  ! le rapport ne depend pas de la reflexion.
      do ipr = 1,n_atom_proto
        if( abs( f_ms(ipr,ipl) ) > eps10 ) then
          rap_lsur2s = 0.5_db * f_mo(ipr,ipl) / f_ms(ipr,ipl)
        else
          rap_lsur2s = 0._db
        endif
        write(3,246) ipr, rap_lsur2s
      end do
    endif

    write(3,250)
    do ipr = 1,n_atom_proto
      if( nspin == 1 ) then
        write(3,260) ipr, numat( itypepr(ipr) )
        do ipl = 1,npldafs
          write(3,270) hkl_dafs(:,ipl), f0(ipr,ipl), fp(ipr), fpp(ipr), phd(ipr,ipl), phd_f0(ipr,ipl), phd_fan(ipr,ipl)
        end do
      else
        write(3,280) ipr, numat( itypepr(ipr) )
        do ipl = 1,npldafs
          write(3,290) hkl_dafs(:,ipl), f0(ipr,ipl), fp(ipr), fpp(ipr), f_ms(ipr,ipl), f_mo(ipr,ipl), phd(ipr,ipl), &
            phd_f0(ipr,ipl), phd_fan(ipr,ipl), phd_fmag(ipr,ipl), v_vec_b(ipl), v_vec_a(ipl)
        end do
      endif
    end do

    write(3,300)
    do ipl = 1,npldafs
      write(3,270) hkl_dafs(:,ipl), phdf0t(ipl,1)
    end do

  endif

! On passe en base R1
  pp  = abs( Rot_int(1,1) - 1._db ) + abs( Rot_int(2,2) - 1._db ) + abs( Rot_int(3,3) - 1._db )
  if( .not. Base_spin .and.  pp > eps10 ) then
    do ipl = 1,npldafs
      do ip = 1,nphi_dafs(ipl)
        v(:) = real( poldafse(:,ipl,ip), db )
        v = matmul( Rot_int, v )
        w(:) = aimag( poldafse(:,ipl,ip) )
        w = matmul( rot_int, w )
        poldafse(:,ipl,ip) = cmplx( v(:), w(:), db )

        v(:) = real( poldafss(:,ipl,ip), db )
        v = matmul( rot_int, v )
        w(:) = aimag( poldafss(:,ipl,ip) )
        w = matmul( rot_int, w )
        poldafss(:,ipl,ip) = cmplx( v(:), w(:), db )

        v(:) = vecdafse(:,ipl,ip)
        v = matmul( rot_int, v )
        vecdafse(:,ipl,ip) = v(:)

        v(:) = vecdafss(:,ipl,ip)
        v = matmul( rot_int, v )
        vecdafss(:,ipl,ip) = v(:)
      end do
    end do
  endif

  if( icheck > 1 ) then
    iwrite = 1
    do ipl = 1,npldafs
      if( nphi_dafs(ipl) == 1 ) cycle
      if( iwrite == 1 ) then
        if( Base_spin ) then
          write(3,370)
        else
          write(3,375)
        endif
        iwrite = 0
      endif
      write(3,376) ipl
      do ip = 1,nphi_dafs(ipl)
        write(3,*)
        do i = 1,3
         write(3,380) poldafse(i,ipl,ip), vecdafse(i,ipl,ip), poldafss(i,ipl,ip), vecdafss(i,ipl,ip)
        end do
      end do
    end do
  endif

  return
  110 format(/' ---- Prepdafs -----',100('-'))
  120 format(//' The reflection number',i3,' : (h,k,l) = (',3i3,') does not exist at this energy !'//)
  130 format(/' (h,k,l) = (',3i3,')')
  140 format('  Reciprocal mesh base (A-1) :',/5x,'X',8x,'Y',8x,'Z', 3(/3f9.5))
  150 format('  Local base (A-1) :',/5x,'I',8x,'J',8x,'Q', 3(/3f9.5))
  160 format(' Calculations with (h,k,l) = (0,0,0)  need non zero values for the incoming and outgoing'/, &
        ' wave vector in order to calculate the scattering angle !')
  165 format(1x,3i3,i5,10x,f9.5)
  170 format(/'  (h,k,l) ThetaBragg  d_hkl (A)    Polarization angles        poldafse',5x, &
      'vecdafse',7x,'poldafss',5x,'vecdafss'/)
  180 format(3i3,2f10.3,2x,3f8.3,2(2x,2f8.5,2x,f8.5)/, 2(55x,2(2x,2f8.5,2x,f8.5)/))
  182 format(/' ipl =',i3,//'  Spin_axis          vec_b',17x,'vec_a')
  185 format(f10.5,2(2x,2f10.5))
  187 format(/' V.Vec_b =',2f10.5,',  V.Vec_a =',2f10.5)
  200 format(/'  (h, k, l)   pol',11x,'Poldafse',20x, 'Poldafss',20x,'Vecdafse',20x,'Vecdafss')
  210 format(/' (h,k,l) = (',3i3,'),  isigpi =',2i3,/' angle', &
  8x,'real(Poldafse)',13x,'imag(Poldafse)',11x,'real(Poldafss)',12x, 'imag(Poldafse)',15x,'Vecdafse',15x,'Vecdafss')
  215 format(3i4,1x,a4, 6(1x,3f9.5))
  220 format(f6.1, 6(1x,3f9.5))
  230 format(/'  (h,k,l)  exp(i*Q.R_ia) (ia = 1,natomsym)')
  240 format(3i3,1p,48(1x,2e13.5))
  245 format(/' Attenuation factor for non resonant magnetic structure', ' factor :'/,'      General : fma =',f7.3,/ &
    '      Orbital : Site    L/2S     (f_orb = (L/S)*f_spin)')
  246 format(10x,i9,f9.3)
  250 format(/' Value of the structure factors per atom site :')
  260 format(/' Site =',i3,', Z =',i3,//, '  (h,k,l)',6x,'f0',14x,'fp',10x,'fpp',10x,'Ph = Sum(exp(iQR))', &
    10x,'Ph*f0*(eps_e.eps_s)     Ph*(fp+i*fpp)*(eps_e.eps_s)')
  270 format(3i3,1p,e13.5,2x,2e13.5,5(2x,2e13.5))
  280 format(/' Site =',i3,', Z =',i3,//, '  (h,k,l)',6x,'f0',14x,'fp',10x,'fpp',11x,'f_spin',9x, &
    'f_orb',10x,'Ph = Sum(exp(iQR))',10x,'Ph*f0*(eps_e.eps_s)    Ph*(fp+i*fpp)*(eps_e.eps_s)', &
    ' i*Ph_m*fma*(f_spin*vs+f_orb*vo)   vs = spin.vec_b',13x, 'vo = spin.vec_a')
  290 format(3i3,1p,e13.5,2x,2e13.5,2(2x,e13.5),6(2x,2e13.5))
  300 format(/'  (h,k,l)   Total Structure factor')
  370 format(/'  Polarization and wave vectors in the internal basis R4 ( orthogonal basis, z along spin direction )',// &
              '         pol            vec           pls',12x,'vos')
  375 format(/'  Polarization and wave vectors in the internal basis R1 ( orthogonal basis, z along c crystal )',// &
              '         pol            vec           pls',12x,'vos')
  376 format(/' ipl = ',i5)
  380 format(2(2f9.5,1x,f9.5,1x))

end

!***********************************************************************

subroutine col_dafs_name(angpoldafs,Bormann,Full_self_abs, hkl_dafs,isigpi,mpirank,ncolm,ncolr,ncolt, &
         nomabs,npldafs,Self_abs)

  use declarations
  implicit none

  integer:: i, icol, index, ipldafs, ipr, j, k, mpirank, ncolm, ncolr, ncolt, npldafs
  integer, dimension(3):: hkl
  integer, dimension(3,npldafs):: hkl_dafs
  integer, dimension(npldafs,2):: isigpi

  character(len=Length_word) nomab
  character(len=Length_word), dimension(ncolm):: nomabs

  Logical:: Bormann, Full_self_abs, Self_abs

  real(kind=db):: a
  real(kind=db), dimension(3,npldafs):: angpoldafs

  icol = ncolr

  do ipldafs = 1,npldafs
    icol = icol + 1
    hkl(:) = hkl_dafs(1:3,ipldafs)
    nomab = ' r('
    j = 3
    call trnom(j,hkl,Length_word,nomab)
    j = j + 1
    if( j < Length_word-1 ) nomab(j:j) = ')'
    if( Bormann ) then
      j = j + 1
      if( j > Length_word ) exit
      select case( ipldafs )
        case(1,2,3,10,11,12,19,20,21,28,29,30)
          nomab(j:j) = 'x'
        case(4,5,6,13,14,15,22,23,24,31,32,33)
          nomab(j:j) = 'y'
        case(7,8,9,16,17,18,25,26,27,34,35,36)
          nomab(j:j) = 'z'
      end select
      j = j + 1
      if( j > Length_word ) exit
      select case( ipldafs )
        case(1,4,7,10,13,16,19,22,25,28,31,34)
          nomab(j:j) = 'x'
        case(2,5,8,11,14,17,20,23,26,29,32,35)
          nomab(j:j) = 'y'
        case(3,6,9,12,15,18,21,24,27,30,33,36)
          nomab(j:j) = 'z'
      end select
    else
      do i = 1,2
        j = j + 1
        if( i == 1 ) k = j
        if( j > Length_word ) exit
        select case( isigpi(ipldafs,i) )
          case(1)
            nomab(j:j) = 's'
          case(2)
            nomab(j:j) = 'p'
          case(3)
            nomab(j:j) = 'r'
          case(4)
            nomab(j:j) = 'l'
          case default
            nomab(j:j) = 'a'
        end select
      end do
    endif
    j = j + 1
    if( j < Length_word-1 ) then
      nomab(j:j) = '_'
      if( abs( angpoldafs(3,ipldafs) ) < 9999._db ) then
        a = angpoldafs(3,ipldafs) * 180 / pi
      else
        a = 0._db
      endif
      index = nint( a )
      call ad_number(index,nomab,Length_word)
    endif
    nomabs(icol) = nomab

    icol = icol + 1
    nomab(2:2) = 'i'
    nomabs(icol) = nomab

    if( Full_self_abs ) then
      nomab = ' '
      icol = icol + 1
      select case( mod(ipldafs,4) )
        case(1)
          nomab(2:11) = 'mu_ss_in_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ss_in_i'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ss_ou_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ss_ou_i'
          nomabs(icol) = nomab
        case(2)
          nomab(2:11) = 'mu_sp_in_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_sp_in_i'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_sp_ou_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_sp_ou_i'
          nomabs(icol) = nomab
        case(3)
          nomab(2:11) = 'mu_ps_in_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ps_in_i'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ps_ou_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ps_ou_i'
          nomabs(icol) = nomab
        case(0)
          nomab(2:11) = 'mu_pp_in_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_pp_in_i'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_pp_ou_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_pp_ou_i'
          nomabs(icol) = nomab
      end select
    elseif( self_abs ) then
      icol = icol + 1
      nomab(2:2) = 'A'
      if( k < Length_word ) nomab(k:k+1) = 'in'
      if( k == Length_word ) nomab(k:k) = 'i'
      nomabs(icol) = nomab
      icol = icol + 1
      if( k < Length_word ) nomab(k:k+1) = 'ou'
      if( k == Length_word ) nomab(k:k) = 'o'
      nomabs(icol) = nomab
    endif
  end do

! ncolt correspond au nombre de colonnes total (d'absorption reelle ou
! virtuelle)
  ncolt = icol
  if( ncolt > ncolm .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,105) ncolt, ncolm
    end do
    stop
  endif

  105 format(/' ncolt =',i3,' > ncolm = ',i3, /' Bug in the program !')

  return
end

!***********************************************************************

  real(kind=db) function detmat(mat)

  use declarations
  implicit none

  real(kind=db), dimension(3,3):: mat

  detmat = mat(1,1) * ( mat(2,2)*mat(3,3) - mat(3,2)*mat(2,3) ) - mat(2,1) * ( mat(1,2)*mat(3,3) - mat(3,2)*mat(1,3) ) &
         + mat(3,1) * ( mat(1,2)*mat(2,3) - mat(2,2)*mat(1,3) )

  return
end

!***********************************************************************

function nomat(numat)

  character(len=4) nomat
  character(len=2) Chemical_Symbol

  nomat = ' '
  nomat(1:2) = Chemical_Symbol(numat)

  return
end

!***********************************************************************

subroutine get_fmag(deltak,fmo,fms,ipr,it,lvval,n_atom_proto,nl,nlatm,nr,nrm,nspin,ntype,popatm,psival,rato)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(0:ntype,nlatm):: lvval

  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival

  fms = 0._db
  fmo = 0._db

  do io = 1,nl

    popt = popatm(ipr,io,1) + popatm(ipr,io,nspin)
    spin = 0.5_db * ( popatm(ipr,io,1) - popatm(ipr,io,nspin) )
    pop = abs( popatm(ipr,io,1) - popatm(ipr,io,nspin) )

    if( abs(spin) < eps10 ) cycle

    l = lvval(it,io)
    hund_mo = ( pop * l - pop * ( pop - 1 ) / 2 ) / (2 * abs(spin))
    if( popt < 2 * l + 1 ) hund_mo = - hund_mo

    f = 0._db
! psival est en fait sqrt(4*pi)*r*psi
    do ir = 2,nr-1
      dr = rato(ir+1,it) - rato(ir-1,it)
      f = f + ( psival(ir,io,it)**2 / rato(ir,it) ) * sin( deltak * rato(ir,it) ) * dr
    end do
    f = 0.5_db * f / deltak

    fms = fms + spin * f
    fmo = fmo + hund_mo * 2 * spin * f

  end do

! Facteur empirique
  fmo = 0.2_db * fmo

  return
end

!*********************************************************************

subroutine get_vec_b(pe,ps,vec_b,we,ws)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db), dimension(3):: vec_b, kve, kvs, p, pe, ps, q, ve, vs

  real(kind=db), dimension(3):: we, ws

  ve(:) = cmplx( we(:), 0._db, db )
  vs(:) = cmplx( ws(:), 0._db, db )

  call prodvec_cp( kve, ve, pe )
  call prodvec_cp( kvs, vs, ps )
  call prodvec_cp( p, ps, pe )
  call prodvec_cp( q, kvs, kve )

  vec_b(:) = p(:) + dot_product(vs,pe) * kvs(:) - dot_product(ve,ps) * kve(:) - q(:)

  return
end

!*********************************************************************

! Calcul de A" dans la formule de Blume et Gibbs, PRB 37,1779 (1988).
! en fait modifie a cause erreur

subroutine get_vec_a(pe,ps,vec_a,we,ws)

  use declarations
  implicit none

  complex(kind=db):: vv
  complex(kind=db), dimension(3):: vec_a, p, pe, ps, ve, ves, vs
!      complex(kind=db), dimension(3):: kve, kvs

  real(kind=db), dimension(3):: we, ws

  ve(:) = cmplx( we(:), 0._db, db )
  vs(:) = cmplx( ws(:), 0._db, db )

  vv = dot_product(ve,vs)

!      call prodvec_cp( kve, ve, pe )
!      call prodvec_cp( kvs, vs, ps )
  call prodvec_cp( p, ps, pe )

! Formule Blume
!      vec_a(:) = 2 * ( 1 - vv ) * p(:) - dot_product(ve,ps) * kve(:)
!     &         + dot_product(vs,pe) * kvs(:)

! Formule modifiee
! La premiere partie de la formule en haut de la page 1781 est bonne
! A" = A' - (A'.K)K, mais le suite est fausse.
!      vec_a(:) = - 2 * ( 1 - vv )
!     &         * ( p(:) + ( dot_product(kve,ps) + dot_product(kvs,pe) )
!     &                  * ( ve(:) - vs(:) ) )

! Remodifie grace a Stephane
  ves(:) = ve(:) - vs(:)

  vec_a(:) = - 2 * ( 1 - vv ) * p(:) - dot_product( ves, p ) * ves(:)

  return
end

!*********************************************************************

! Oana Bunau 15 avril 2007

function DW(q,Z,temp)

  use declarations
  implicit none

  integer Z, nassm
  parameter( nassm=103 )

  real(kind=db), dimension(nassm):: TD

! source: Ashcroft and Mermin; the Debye temperatures correspond to the
! free element
! Unknown Debye temperatures are taken as 1

  data TD /  110,  26, 400,1000,1250,1860,  79,  46,   1,  63, 150, 318, 394, 625,   1,   1,   1,  85, 100, 230, &
             359, 380, 390, 460, 400, 420, 385, 375, 315, 234, 240, 360, 285, 150,   1,  73,  56, 147, 256, 250, &
             275, 380,   1, 382, 350, 275, 215, 120, 129, 170, 200, 139,   1,  55,  40, 110, 132, 139, 152, 157, &
               1, 160, 107, 176, 188, 186, 191, 196, 200, 118, 207,   1, 225, 310, 416, 400, 430, 230, 170, 100, &
              96,  88, 120,   1,   1,   1,   1,   1,   1, 100, 1, 210, 188, 150,   1,   1,   1,   1,   1,   1, 1,   1,   1/

! fonction de Debye, facteur de Debye Waller
  real(kind=db) Debye_function, DW, FD, Mass_atom
! vecteur transfert d'impuls, en 1/A
  real(kind=db) q
  real(kind=db) temp, x

  FD = Debye_function(temp,TD(Z))

! 1 uam = 1.660538921*10**(-27)kg
  x = - 300 * (q*hbar)**2 * FD / ( Mass_atom(Z) * 1.660538921 * k_Boltzmann * TD(Z) )

  DW = exp(x)

! Ecriture dans le fichier d'erreur
  if( abs(TD(Z)) < eps10 ) then
    call write_error
    write(9,*) ' Debye temperature is not available for the element with Z =', Z
    stop
  end if

end

!*********************************************************************

function Debye_function(m,n)

  use declarations
  real(kind=db) m,n,dx,a,x,Debye_function
  integer i,imax
  Debye_function = 0._db

  a = n/m
  dx = 0.0001_db
  imax = Nint(a/dx)
  dx = a/imax
  x = dx/2

  do i = 1,imax
   x = x + dx
   Debye_function = Debye_function + ( dx*x /( EXP(x)-1 ) )
  end do

  Debye_function = (m/n)**2 * Debye_function + 0.25_db

end

!*********************************************************************

! Sousprogramme evaluant les orbitales sondees par les regles de selection des transitions

subroutine etafin(Ylm_comp,icheck,iopsymr,irep_util, jseuil,karact,ldip,lmoins1,loct,lplus1,lqua,lseuil, &
            mpirank,nb_rep,nbseuil,ngrph,nspino,Spinorbite, State_all,Sym_cubic)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db):: gnt

  integer, dimension(2):: ninitl
  integer, dimension(6):: loperat, moperat
  integer, dimension(nopsm):: iopsymr
  integer, dimension(nrepm,2):: irep_util(nrepm,2)
  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua
  integer, dimension(3,3,3):: loct
  integer, dimension(:), allocatable:: iselec, lselec, mselec

  logical:: lmoins1, lplus1, spinorbite, state_all, sym_cubic, Ylm_comp

  real(kind=db):: jinitl(2), jzinitl
  complex(kind=db), dimension(nopsm,nrepm):: karact

  if( icheck > 0 )  write(3,110)

  linitl = lseuil

! Calcul du (l,m) = (linitl,minitl) de l'etat initial
  select case(jseuil)
    case(1,3,5,7)
      isinitl = 1
    case default
      isinitl = -1
  end select

  do iseuil = 1,nbseuil
    if( iseuil == 2 ) isinitl = - isinitl
    jinitl(iseuil) = linitl + 0.5_db * isinitl
    ninitl(iseuil) = nint( 2 * jinitl(iseuil) + 1 )
  end do

  if( Ylm_comp ) then
    nsm = 2
  else
    nsm = 1
  endif

  nlfm0 = 52

  nlfm = nlfm0

  do     ! Boucle dimension

  allocate( iselec(nlfm) )
  allocate( lselec(nlfm) )
  allocate( mselec(nlfm) )

! Calcul des (l,m) = (lselec,mselec) des etats d'arrivee
  lf = 0
! Boucle sur toutes les transitions possibles, on est dans la base de
! l'atome absorbeur.
  do kv2 = 0,3
  do kvo = 0,3
  do kpl = 1,3
    if( kvo == 0 ) then
      if( ldip(kpl) == 0 ) cycle
    else
      if( kvo < kpl ) cycle
    endif
    if( kv2 == 0 ) then
      if( kvo > 0 ) then
        if( lqua(kpl,kvo) == 0 ) cycle
      endif
    else
      if( kv2 < kvo .or. kvo == 0 ) cycle
      if( loct(kpl,kvo,kv2) == 0 ) cycle
    endif

! Calcul du (l,m) = (loperat,moperat) de la transition
    call lmtrans(kpl,kvo,kv2,noperat,loperat,moperat,spinorbite)
    if( icheck > 1 ) then
      if( kvo == 0 ) then
        write(3,130) kpl
      elseif( kv2 == 0 ) then
        write(3,140) kpl, kvo
      else
        write(3,150) kpl, kvo, kv2
      endif
      do iop = 1,noperat
        write(3,160) loperat(iop), moperat(iop)
      end do
    endif

    do iseuil = 1,nbseuil

      do initl = 1,ninitl(iseuil)
        jzinitl = - jinitl(iseuil) + initl - 1

        do ispin = 1,2
          m_initl = nint( jzinitl + ispin - 1.5_db )
          jspin = min(ispin,nspino)
          if( m_initl > linitl .or. m_initl < -linitl ) cycle

          do iop = 1,noperat
            lo = loperat(iop)
            mo = moperat(iop)

            lf1 = abs( linitl - lo )
            lf2 = linitl + lo
            if( lplus1 ) lf1 = lf2
            if( lmoins1 ) lf2 = lf1

            m = mo + m_initl

            do ispo = 1,nspino
              if( ispo == 1 ) then
                kspin = jspin
              else
                kspin = 3 - jspin
                if( jspin == 1 ) then
                  m = m + 1
                else
                  m = m - 1
                endif
              endif

            boucle_ll: do l = lf1,lf2,2

              if( m > l .or. m < -l ) cycle

              if( Ylm_comp ) then

                boucle_ism: do ism = 1,nsm
                  if( ism == 1 ) then
                    mm = m
                  else
                    if( m == 0 ) cycle
                    mm = - m
                  endif
                  do lg = 1,min(lf,nlfm)
                    if( lselec(lg) == l .and. mselec(lg) == mm .and. iselec(lg) == kspin ) cycle boucle_ism
                  end do
                  lf = lf + 1
                  if( lf > nlfm ) cycle
                  lselec(lf) = l
                  mselec(lf) = mm
                  iselec(lf) = kspin
                end do boucle_ism

              else  ! on passe en base reelle

                do ms = -1,1,2
                  if( m_initl == 0 .and. ms == 1 ) exit
                  mi = ms * abs( m_initl )
                  gnt = gauntc(l,m,lo,mo,linitl,mi)
                  if( abs(gnt) < eps6 ) cycle
                  do lg = 1,min(lf,nlfm)
                    if( lselec(lg) == l .and. mselec(lg) == m ) cycle boucle_ll
                  end do
                  lf = lf + 1
                  if( lf > nlfm ) cycle
                  lselec(lf) = l
                  mselec(lf) = m
                  iselec(lf) = kspin
                end do

              endif

            end do boucle_ll
          end do
          end do

        end do  ! fin de la boucle sur les operateurs

       end do
    end do
  end do
  end do
  end do
  nselec = lf

  if( sym_cubic .and. .not. Ylm_comp ) then
    do lf = 1,nselec
      if( lselec(lf) == 2 ) exit
    end do
    if( lf < nselec + 1 ) then
      do lg = lf,min(nselec,nlfm)
        if( lselec(lg) == 2 .and. mselec(lg) == - 2 ) exit
      end do
      if( lg == nselec + 1 ) then
        nselec = nselec + 1
        if( nselec <= nlfm ) then
         lselec(nselec) = 2
         mselec(nselec) = -2
         iselec(nselec) = 1
        endif
      endif
    endif
  endif

    if( nselec <= nlfm ) exit

    nlfm = 2 * nlfm
    deallocate( iselec, lselec, mselec )

  end do   ! fin boucle dimension

  if( mpirank == 0 .and. nselec > nlfm0 ) write(3,170) nselec

  if( state_all ) then
    call irep_util_all(icheck,iopsymr,irep_util,karact,nb_rep, ngrph,nspino,Spinorbite)
  else
    call cal_irep(icheck,iopsymr,irep_util,iselec,karact,lselec, mselec,nb_rep,ngrph,nlfm,nselec,nspino,Spinorbite)
  endif

  deallocate( iselec, lselec, mselec )

  return
  110 format(/' ---- Etafin -------',100('-'))
  130 format(/' Dipole, kpl =',i2)
  140 format(/' Quadrupole, kpl =',i2,', kvo =',i2)
  150 format(/' Octupole, kpl =',i2,', kvo =',i2,', kv2 =',i2)
  160 format('   (loperat,moperat) = (',2i2,')')
  170 format(/' nselec =',i2)
end

!*********************************************************************

subroutine cal_irep(icheck,iopsymr,irep_util,iselec,karact,lselec, mselec,nb_rep,ngrph,nlfm,nselec,nspino,spinorbite)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  character(len=9):: nomsym

  complex(kind=db), dimension(nopsm,nrepm):: karact
  complex(kind=db), dimension(nopsm):: kopsymo

  integer, dimension(nopsm):: iopsymr
  integer, dimension(nrepm):: nlmrep
  integer, dimension(nselec):: irepo, lass
  integer, dimension(nrepm,2):: irep_util(nrepm,2)
  integer, dimension(nrepm,nlfm):: lmrep
  integer, dimension(nselec):: iselec, lselec, mselec

  logical spinorbite
  logical, dimension(nselec):: fait

  if( icheck > 1 ) then
    write(3,110)
    write(3,120)
    call write_iopsym(iopsymr,3)
  endif

  irepo(:) = 0
  lass(:) = 0
  fait(:) = .false.
  do lf = 1,nselec

    if( fait(lf) ) cycle

    l = lselec(lf)
    m = mselec(lf)
    ispin = iselec(lf)
    fait(lf) = .true.

    do isp = 1,nspino

      if( isp == 2 ) then
        if( ispin == 1 ) then
          m = m + 1
        else
          m = m - 1
        endif
        ispin = 3 - ispin
        if( m > l .or. m < -l ) cycle
        do lg = 1,nselec
          if( fait(lg) ) cycle
          if( lselec(lg) == l .and. mselec(lg) == m .and. iselec(lg) == ispin ) then
            fait(lg) = .true.
            lass(lf) = lg
            lass(lg) = lf
            exit
          endif
        end do
      endif

! Calcul des characteres des orbitales d'arrivee :
! Dans la base de l'atome absorbeur
      call symorb(l,m,kopsymo)

! Recherche de la representation a laquelle appartient l'orbitale
      boucle_irep: do irep = 1,nb_rep

        do is = 1,nopsm
          if( iopsymr(is) == 0 ) cycle
          if( abs( karact(is,irep) - kopsymo(is) ) > eps10 ) cycle boucle_irep
        end do
        if( isp == 1 ) then
          irepo(lf) = irep
        else
          irepo(lg) = irep
        endif
        exit

      end do boucle_irep

    end do
  end do

! Cas des representations conjuguees
  if( spinorbite ) then
    do lf = 1,nselec
      if( irepo(lf) /= 0 ) cycle
      l = lselec(lf)
      m = - mselec(lf)
      do lg = 1,nselec
        if( l == lselec(lg) .and. m == mselec(lg) ) exit
      end do
      if( lg == nselec + 1 ) then
        call symorb(l,m,kopsymo)
        boucle_irep2: do irep = 1,nb_rep
          do is = 1,nopsm
            if( iopsymr(is) == 0 ) cycle
            if( abs( karact(is,irep) - kopsymo(is) ) > eps10 ) cycle boucle_irep2
          end do
          irepo(lf) = - irep
        end do boucle_irep2
      else
        irepo(lf) = - irepo(lg)
      endif
    end do
  endif

  if(icheck > 1 .and. lf <= nlfm ) then
    write(3,'(/A)') '  lselec mselec iselec  irep   lass'
    do lf = 1,nselec
      write(3,190) lselec(lf), mselec(lf), iselec(lf), irepo(lf), lass(lf)
    end do
  endif

  ngrph = 0
  irep_util(:,:) = 0

  if( nb_rep == 1 ) then
    ngrph = 1
    irep_util(ngrph,:) = 1
  else
    boucle_l: do lf = 1,nselec
      if( irepo(lf) == 0 ) cycle
      lg = lass(lf)
      if( lg == 0 ) cycle
      ispf = iselec(lf)
      ispg = nspino + 1 - ispf
      do igrph = 1,ngrph
        if( irep_util(igrph,ispf) == irepo(lf) .and. irep_util(igrph,ispg) == irepo(lg) ) cycle boucle_l
      end do
      ngrph = ngrph + 1
      irep_util(ngrph,ispf) = irepo(lf)
      irep_util(ngrph,ispg) = irepo(lg)
    end do boucle_l

    boucle_l2: do lf = 1,nselec
      if( irepo(lf) == 0 ) cycle
      lg = lass(lf)
      if( lg /= 0 ) cycle
      ispf = iselec(lf)
      do igrph = 1,ngrph
        if( irep_util(igrph,ispf) == irepo(lf) ) cycle boucle_l2
      end do
      ngrph = ngrph + 1
      irep_util(ngrph,ispf) = irepo(lf)
    end do boucle_l2
  endif

  nlmrep(:) = 0
  boucle_igrph: do igrph = 1,ngrph
    do lf = 1,nselec
      if( irepo(lf) /= irep_util(igrph,iselec(lf)) ) cycle
      nlmrep(igrph) = nlmrep(igrph) + 1
      lmrep(igrph,nlmrep(igrph)) = lf
    end do
  end do boucle_igrph

! Cas ou aucune fonction n'a ete trouve pour une des 2 representations
  if( spinorbite ) then
    do igrph = 1,ngrph
      if( irep_util(igrph,1) /= 0 .and. irep_util(igrph,nspino) /= 0 ) cycle
      lf = lmrep(igrph,1)
      isp = iselec(lf)
      l = lselec(lf)
      m = mselec(lf)

      if( isp == 1 ) then
        m = m + 1
      else
        m = m - 1
      endif
      l = l + 4

      boucle_isg: do isg = 1,-1,-2
        m = isg * m
        call symorb(l,m,kopsymo)

        boucle_irep3: do irep = 1,nb_rep

          do is = 1,nopsm
            if( iopsymr(is) == 0 ) cycle
            if( abs( karact(is,irep) - kopsymo(is) ) > eps10 ) cycle boucle_irep3
          end do
          irep_util(igrph,nspino+1-isp) = isg * irep
          exit boucle_isg

        end do boucle_irep3
      end do boucle_isg
    end do
  endif

  if( icheck > 0 ) then
    do igrph = 1,ngrph
      write(3,220) igrph
      if( spinorbite ) then
        write(3,230) irep_util(igrph,1:nspino)
      else
        write(3,240) irep_util(igrph,1:nspino)
      endif
      write(3,'(/A)') '  Orbitals belonging to the representation :'
      if( spinorbite ) then
        write(3,250) ( lselec(lmrep(igrph,lm)), mselec(lmrep(igrph,lm)), iselec(lmrep(igrph,lm)), lm = 1,nlmrep(igrph) )
      else
        write(3,260) ( lselec(lmrep(igrph,lm)), mselec(lmrep(igrph,lm)), lm = 1,nlmrep(igrph) )
      endif
      write(3,'(/A)') '     is  Symmetry    Character'
      do isp = 1,nspino
        irep = abs(irep_util(igrph,isp))
        if( irep == 0 ) cycle
        do is = 1,nopsm
          if( abs( karact(is,irep) ) < eps10 ) cycle
          if( irep_util(igrph,isp) > 0 ) then
            write(3,270) is, nomsym(is), karact(is,irep)
          else
            write(3,270) is, nomsym(is), conjg( karact(is,irep) )
          endif
        end do
      end do
    end do
  endif

  return
  110 format(/' ---- Cal_irep -------',100('-'))
  120 format(/' iopsymr =')
  190 format(5i7)
  220 format(/' Useful representation number',i3,' :')
  230 format(/'   Number of the corresponding representation :',i2, ' x',i2)
  240 format(/'   Number of the corresponding representation :',i2)
  250 format(30(2x,3i3))
  260 format(30(2x,2i3))
  270 format(i7,a9,13(1x,2f7.3))
end

!*********************************************************************

! Stockage dans irep_util de toutes les representations
! Utilisee dans le cas d'un calcul autocoherent car pour calculer
! la densite totale, il faut avoir toutes les representations.

subroutine irep_util_all(icheck,iopsymr,irep_util,karact,nb_rep, ngrph,nspino,spinorbite)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(nopsm):: iopsymr
  integer, dimension(:), allocatable:: iselec, lselec, mselec
  integer, dimension(nrepm,2):: irep_util(nrepm,2)

  complex(kind=db), dimension(nopsm,nrepm):: karact

  logical:: Spinorbite

  nlfm = 52

  allocate( iselec(nlfm) )
  allocate( lselec(nlfm) )
  allocate( mselec(nlfm) )

  nselec = 0

  do l = 0,3
    do m = -l,l
      do isp = 1,nspino
        nselec = nselec + 1
        lselec(nselec) = l
        mselec(nselec) = m
        iselec(nselec) = isp
      end do
    end do
  end do

  call cal_irep(icheck,iopsymr,irep_util,iselec,karact,lselec, mselec,nb_rep,ngrph,nlfm,nselec,nspino,Spinorbite)

  deallocate( iselec, lselec, mselec )

  return
end

!*********************************************************************

! Sousprogramme calculant les valeurs des vecteurs polarisation et onde

subroutine polond(Base_spin,Dipmag,icheck,ltypcal,moyenne,mpirank,msymdd, &
        msymqq,ncolm,ncolr,nomabs,nple,nplr,nplrm,nxanout,Octupole,Orthmat,Orthmati,pdp, &
        pdpolar,pol,polar,Polarise,Quadrupole,Rot_int,Veconde,Vec,Xan_atom)

  use declarations
  implicit none

  integer:: i, icheck, ii, ipl, ipr, j, jj, jpl, k, mpirank, ncolm, ncolr, nj, nple, nplr, nplrm, nxanout

  character(len=Length_word):: nomab
  character(len=Length_word), dimension(ncolm):: nomabs
  character(len=13), dimension(nplrm):: ltypcal

  complex(kind=db), dimension(3,nplrm):: pol

  integer, dimension(3):: ipol, kpol
  integer, dimension(3,3):: msymdd
  integer, dimension(3,3,3,3):: msymqq

  logical:: Base_spin, Dipmag, Moyenne, Octupole, Polarise, Quadrupole, Xan_atom

  real(kind=db):: pdt, plmin, pp, pp1, pp2, pv, r, rac_2, vomin, vv
  real(kind=db), dimension(3):: p, pl, v, v1, v2, vo, w
  real(kind=db), dimension(3,3):: Orthmat, Orthmati, Rot_int
  real(kind=db), dimension(3,nple):: polar, veconde
  real(kind=db), dimension(nple,2):: pdpolar
  real(kind=db), dimension(3,nplrm):: vec
  real(kind=db), dimension(nplrm,2):: pdp

  if( icheck > 0 ) write(3,100)

  pol(:,:) = 0._db
  vec(:,:) = 0._db
  pdp(:,:) = 0._db

  rac_2 = sqrt( 2._db )
  jpl = 0

! Si aucune polarisation n'est definie en entree, on en construit par
! defaut
! On les definit dans le repere orthonorme interne.
  if( nple == 0 ) then

    if( Quadrupole .or. Octupole ) then

      boucle_ij: do i = 1,3
        j = mod(i,3) + 1
        do ii = 1,i-1
          jj = mod(ii,3) + 1
          if( abs(msymdd(i,i)) == abs(msymdd(ii,ii)) .and. abs(msymqq(i,j,i,j)) == abs(msymqq(ii,jj,ii,jj) ) ) &
               cycle boucle_ij
        end do
        jpl = jpl + 1
        pol(i,jpl) = 1._db
        vec(j,jpl) = 1._db
        jpl = jpl + 1
        pol(i,jpl) = 1._db / rac_2
        vec(i,jpl) = 1._db / rac_2
        pol(j,jpl) = 1._db / rac_2
        vec(j,jpl) = - 1._db / rac_2
        do ii = 1,3
          jj = mod(ii,3) + 1
          if( abs(msymdd(i,i)) == abs(msymdd(ii,ii)) .and. abs(msymqq(i,j,i,j)) == abs(msymqq(ii,jj,ii,jj)) ) then
            pdp(jpl-1,1) = pdp(jpl-1,1) + 1._db
            pdp(jpl-1,2) = pdp(jpl-1,2) + 1.5_db
            pdp(jpl,2) = pdp(jpl,2) + 1._db
          endif
        end do
      end do boucle_ij

    else

      boucle_exter: do i = 1,3
        do ii = 1,i-1
          if( abs( msymdd(i,i) ) == abs( msymdd(ii,ii)) ) cycle boucle_exter
        end do
        jpl = jpl + 1
        pol(i,jpl) = 1._db
        do ii = 1,3
          if( abs(msymdd(i,i)) == abs( msymdd(ii,ii)) ) pdp(jpl,1) = pdp(jpl,1) + 1._db
        end do
      end do boucle_exter

      if( dipmag ) then
        select case(jpl)
          case(1)
            if( abs( pol(1,1) ) < 0.99_db ) then
              v(2:3) = 0._db; v(1) = 1._db
            else
              v(1:2) = 0._db; v(3) = 1._db
            endif
            p(1:3) = real( pol(1:3,1), db )
            call prodvec(w,p,v)
            w = w / sqrt( sum( w(:)**2 ) )
            call prodvec(v,w,p)
            vec(1:3,1) = v(1:3)
          case(2)
            p(1:3) = real( pol(1:3,1), db )
            w(1:3) = real( pol(1:3,2), db )
            call prodvec(v,p,w)
            vec(1:3,1) = v
            vec(1:3,2) = real( pol(1:3,1) )
          case(3)
            vec(1:3,1) = real( pol(1:3,2), db )
            vec(1:3,2) = real( pol(1:3,3), db )
            vec(1:3,3) = real( pol(1:3,1), db )
        end select
      endif

    endif

    nplr = jpl

    pdt = sum( pdp(1:nplr,1) )
    pdp(1:nplr,1) = pdp(1:nplr,1) / pdt
    if( Quadrupole .or. Octupole ) then
      pdt = sum( pdp(1:nplr,2) )
      pdp(1:nplr,2) = pdp(1:nplr,2) / pdt
    endif

    ltypcal(:) = 'xanes rectil'

  else

    do ipl = 1,nple

      pl(:) = 0._db
      vo(:) = 0._db

! On passe en base interne (orthonormee)
      p(:) = polar(:,ipl)
      pp = sum( p(:)**2 )
      if( abs(pp) > eps6 ) call trvec(mpirank,Orthmat,p,pl)

      v(:) = veconde(:,ipl)
      vv = sum( v(:)**2 )
      if( vv > eps6 ) call trvec(mpirank,Orthmat,v,vo)

! Test sur l'orthogonalite
      if( Quadrupole .or. Octupole ) then
        pv = abs( sum( vo(:) * pl(:) ) )
        if( pv > eps4 ) then
          if( mpirank == 0 ) then
            call write_error
            do ipr = 3,9,3
              write(ipr,110) ipl, p(:), v(:)
            end do
          endif
          stop
        endif
      endif

      jpl = jpl + 1

      if( pp > eps6 ) then

        pol(:,jpl) = cmplx( pl(1:3), 0._db, db)
        vec(:,jpl) = vo(:)
        pdp(jpl,:) = pdpolar(ipl,:)
        ltypcal(jpl) = 'xanes rectil'

      else

! Si le vecteur polarisation est nul, la polarisation est circulaire
        do i = 1,3
          if( abs( vo(i) ) > eps4 ) exit
        end do
        j = i - 1
        if( j == 0 ) j = 3
        v(:) = 0._db
        v(j) = 1._db
        call prodvec(v1,v,vo)
        r = sqrt( sum( v1(:)**2 ) )
        v1(:) = v1(:) / r
        call prodvec(v2,vo,v1)
        pol(:,jpl) = cmplx( v1(:), v2(:), db ) / rac_2
        vec(:,jpl) = vo(:)
        pdp(jpl,:) = 0.5_db*pdpolar(ipl,:)
        ltypcal(jpl) = 'xanes circ g'

        jpl = jpl + 1
        pol(:,jpl) = cmplx( v1(:), - v2(:), db ) / rac_2
        vec(:,jpl) = vo(:)
        pdp(jpl,:) = 0.5_db*pdpolar(ipl,:)
        ltypcal(jpl) = 'xanes circ d'

      endif

    end do

    nplr = jpl

  endif

! Determination des noms des colonnes de resultats dans les fichiers de
! sortie

  jpl = 0   ! indice de colonne

  do ipl = 1,nplr

    jpl = jpl + 1

    if( ltypcal(ipl) == 'xanes circ g' ) then

      nomabs(jpl) = '   left_pol  '

    elseif( ltypcal(ipl) == 'xanes circ d' ) then

      nomabs(jpl) = '  right_pol  '

! On ajoute la colonne difference
      jpl = jpl + 1
      vo(:) = vec(:,ipl)
      vo = matmul( orthmati, vo)
      vomin = 1._db
      do k = 1,3
        pp = abs( vo(k) )
        if( pp > eps4 ) vomin = min( vomin, pp)
      end do
      kpol(:) = nint( vo(:) / vomin )

      nomab = ' dic('
      j = 5
      call trnom(j,kpol,Length_word,nomab)
      j = j + 1
      if( j <= Length_word ) nomab(j:j) = ')'
      nomabs(jpl) = nomab

    else

! On se place dans la base cristallographique
      pl(:) = pol(:,ipl)
      pl = matmul( orthmati, pl)
      vo(:) = vec(:,ipl)
      vo = matmul( orthmati, vo)

      plmin = 1._db
      do k = 1,3
        pp = abs( pl(k) )
        if( pp > eps4 ) plmin = min( plmin, pp )
      end do
      ipol(:) = nint( pl(:) / plmin )
      vomin = 1._db
      do k = 1,3
        pp = abs( vo(k) )
        if( pp > eps4 ) vomin = min( vomin, pp)
      end do
      kpol(:) = nint( vo(:) / vomin )

      j = 0
      nomab = ' '
      call trnom(j,ipol,Length_word,nomab)
      if( ( kpol(1) /= 0 .or. kpol(2) /= 0 .or. kpol(3) /= 0 ) .and. ( Quadrupole .or. Octupole .or. Dipmag ) ) then
        j = j + 1
        nomab(j:j) = ','
        call trnom(j,kpol,Length_word,nomab)
      endif
      nj = j
      if( nj < Length_word-2 ) nomab(nj+1:nj+1) = ')'
      if( nj < Length_word-1 ) then
        nomab(1:Length_word) = '(' // nomab(1:Length_word-1)
      else
        nomab(1:Length_word) = ' ' // nomab(1:Length_word-1)
      endif
      nomabs(jpl) = nomab
    endif
  end do

  pp1 = sum( pdp(1:nplr,1) )
  pp2 = sum( pdp(1:nplr,2) )
  pp = sum( abs(pdp(1:nplr,:)) )
  if(.not. Polarise ) jpl = nplr
  if( nplr > 1 .and. pp > eps4 ) then
    if( abs(pp1) > eps4 ) pdp(1:nplr,1) = pdp(1:nplr,1) / pp1
    if( ( Quadrupole .or. Octupole ) .and. abs(pp2) > eps4 ) pdp(1:nplr,2) = pdp(1:nplr,2) / pp2
    Moyenne = .true.
    jpl = jpl + 1
    if( abs(pp1+pp2) > eps4 ) then
      nomabs(jpl) = '    <xanes>  '
    else
      nomabs(jpl) = '    dichroi  '
    endif
  else
    Moyenne = .false.
  endif

! nxanout est l'indice de colonne a partir de laquelle on imprime
  if( Polarise ) then
    nxanout = 1
  else
    if( nplr == 1 ) nomabs(1) = '    <xanes>  '
    nxanout = jpl
  endif

  if( xan_atom ) then
    jpl = jpl + 1
    nomabs(jpl) = '  XANES_atom '
  endif

! ncolr correspond au nombre de colonnes ecrites pour le xanes
  ncolr = jpl

! On passe en base R1
  pp  = abs( rot_int(1,1) - 1._db ) + abs( rot_int(2,2) - 1._db ) + abs( rot_int(3,3) - 1._db )
  if( .not. Base_spin .and.  pp > eps10 ) then
    do ipl = 1,nplr

      v(:) = real( pol(:,ipl), db )
      v = matmul( rot_int, v )

      w(:) = aimag( pol(:,ipl) )
      w = matmul( rot_int, w )
      pol(:,ipl) = cmplx( v(:), w(:), db )

      v(:) = vec(:,ipl)
      v = matmul( rot_int, v )
      vec(:,ipl) = v(:)

    end do

  endif

  if( icheck > 0 ) then
    if( Base_spin ) then
      write(3,150)
    else
      write(3,155)
    endif
    do ipl = 1,nplr
      write(3,160)
      write(3,160) pol(1,ipl), vec(1,ipl), ltypcal(ipl), pdp(ipl,:)
      do i= 2,3
        write(3,160) pol(i,ipl), vec(i,ipl)
      end do
    end do
  endif

  return
  100 format(/' ---- Polond -----',100('-'))
  110 format(/' The incoming wave vector and the polarization must be perpendicular !'/, &
              ' ipl =',i3,'   pol =',3f10.5,/ 12x,'vec =',3f10.5)
  150 format(/'  Polarization and wave vectors in the internal basis R4 ( orthogonal basis, z along spin direction )',// &
              '         pol           vec        type      weight_d weight_q')
  155 format(/'  Polarization and wave vectors in the internal basis R1 ( orthogonal basis, z along c crystal )',// &
              '         pol           vec        type      weight_d weight_q')
  160 format(2f9.5,1x,f9.5,3x,a12,2f9.5)
end

!***********************************************************************

subroutine trnom(j,ind,Length_word,nomab)

  implicit none

  integer:: j, k, Length_word

  character(len=Length_word):: nomab, nomac

  integer, dimension(3):: ind

  nomac = ' '
  nomac(1:Length_word-1) = nomab(1:Length_word-1)
  do k = 1,3
    call ad_number(ind(k),nomac,Length_word)
  end do
  nomab(1:Length_word-1) = nomac(1:Length_word-1)
  j = len_trim(nomab)

  return
end

!***********************************************************************

subroutine trvec(mpirank,Orthmat,ve,vs)

  use declarations
  implicit none

  integer:: ipr, mpirank

  real(kind=db):: pp
  real(kind=db), dimension(3):: ve, vs, w
  real(kind=db), dimension(3,3):: Orthmat

! On se place en base orthonormee
  w = matmul( Orthmat, ve )

! Normalisation
  pp = sqrt( sum( w(:)**2 ) )
  if( pp > eps4 ) then
    w(:) = w(:) / pp
  elseif( mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,110)
    end do
    stop
  endif

  vs(:) = w(:)

  return
  110 format(//' The polarisation vector is zero in routine trvec !'/)
end

!***********************************************************************

! Calcul du (l,m) = (loperat,moperat) de la transition

subroutine lmtrans(kpl,kvo,kv2,noperat,loperat,moperat,spinorbite)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(6) :: loperat, moperat
  logical spinorbite

! dipole
  if( kvo == 0 ) then

    noperat = 1
    loperat(1) = 1
    select case(kpl)
      case(1)
        moperat(1) = 1
      case(2)
        moperat(1) = -1
      case(3)
        moperat(1) = 0
    end select

! quadrupole
  elseif( kv2 == 0 ) then

    kprod = kpl * kvo

    select case(kprod)
      case(2)
        noperat = 1
        loperat(1) = 2
        moperat(1) = -2
      case(6)
        noperat = 1
        loperat(1) = 2
        moperat(1) = -1
      case(3)
        noperat = 1
        loperat(1) = 2
        moperat(1) = 1
      case(1,4)
        noperat = 3
        loperat(1) = 0
        moperat(1) = 0
        loperat(2:3) = 2
        moperat(2) = 0
        moperat(3) = 2
      case(9)
        noperat = 2
        loperat(1) = 0
        moperat(1) = 0
        loperat(2) = 2
        moperat(2) = 0
    end select

! octupole
  else

    kprod = kpl * kvo * kv2

    select case(kprod)
      case(6)
        noperat = 1
        loperat(1) = 3
        moperat(1) = -2
      case(2)
        noperat = 3
        loperat(1) = 1
        moperat(1) = -1
        loperat(2:3) = 3
        moperat(2) = -3
        moperat(3) = -1
      case(4)
        noperat = 3
        loperat(1) = 1
        moperat(1) = 1
        loperat(2:3) = 3
        moperat(2) = 3
        moperat(3) = 1
      case(3,12)
        noperat = 3
        loperat(1) = 1
        moperat(1) = 0
        loperat(2:3) = 3
        moperat(2) = 0
        moperat(3) = 2
      case(18)
        noperat = 2
        loperat(1) = 1
        moperat(1) = -1
        loperat(2) = 3
        moperat(2) = -1
      case(9)
        noperat = 2
        loperat(1) = 1
        moperat(1) = 1
        loperat(2) = 3
        moperat(2) = 1
      case(1)
        noperat = 3
        loperat(1) = 1
        moperat(1) = 1
        loperat(2:3) = 3
        moperat(2) = 3
        moperat(3) = 1
      case(8)
        noperat = 3
        loperat(1) = 1
        moperat(1) = -1
        loperat(2:3) = 3
        moperat(2) = -3
        moperat(3) = -1
      case(27)
        noperat = 2
        loperat(1) = 1
        moperat(1) = 0
        loperat(2) = 3
        moperat(2) = 0
    end select

  endif

  if( spinorbite ) then
    jop = noperat
    do iop = 1,noperat
      if( moperat(iop) == 0 ) cycle
      jop = jop + 1
      loperat(jop) = loperat(iop)
      moperat(jop) = - moperat(iop)
    end do
    noperat = jop
  endif

  return
end

!***********************************************************************

subroutine symorb(l,m,kopsymo)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db), dimension(nopsm):: kopsymo

  integer, dimension(nopsm):: iopsymo

  logical lpair, pair

  if( l == 0 ) then
    kopsymo(:) = (1._db, 0._db)
    return
  endif
  pair = .false.
  iopsymo(:) = 0
  iopsymo(1) = 1
  lpair = mod(l,2) == 0

! Axes 2
  im = abs( m )
  do is = 22,24
    select case(is)
      case(22)   ! Axe 2 Ox
        mm = mod(l+im,2)
        pair = (mm == 0 .and. m >= 0) .or. (mm == 1 .and. m < 0)
      case(23)   ! Axe 2 Oy
        mm = mod(l+2*im,2)
        pair = (mm == 0 .and. m >= 0) .or. (mm == 1 .and. m < 0)
      case(24)   ! Axe 2 Oz
        mm = mod(im,2)
        pair = mm == 0
    end select
    if( pair ) then
      iopsymo(is) = 1
    else
      iopsymo(is) = -1
    endif
  end do

! Axes 2 selon (110) et (1-10)
  m4 = mod(im,4)
  mlm = mod(l+im,2)
  if( ( m >= 0 .and. ( ( m4 == 0 .and. mlm == 0 ) .or. ( m4 == 2 .and. mlm == 1 ) ) ) &
   .or. ( m < 0 .and. ( ( m4 == 2 .and. mlm == 0 ) .or. ( m4 == 0 .and. mlm == 1 ) ) ) ) then
    iopsymo(10:11) = 1
  elseif( ( m >= 0 .and. ( ( m4 == 0 .and. mlm == 1 ) .or. ( m4 == 2 .and. mlm == 0 ) ) ) &
   .or. ( m < 0 .and. ( ( m4 == 2 .and. mlm == 1 ) .or. ( m4 == 0 .and. mlm == 0 ) ) ) ) then
    iopsymo(10:11) = -1
  endif

! Axe 3 Oz
! Quand la valeur est differente de 1 ou - 1, le chiffre represente n
! pour le caractere egal a exp(n*i*pi/12).
  mm = m
  do
    if( mm >= 0 ) exit
    mm = mm + 3
  end do
  mm = mod(mm,3)
  select case(mm)
    case(0)
      iopsymo(49:50) = 1
    case(1)
      iopsymo(49) = 4; iopsymo(50) = 8
    case(2)
      iopsymo(49) = 8; iopsymo(50) = 4
  end select
! Axes 3 negatif
  if( mod(l+m,2) == 0 ) then
    iopsymo(53:54) = iopsymo(49:50)
  else
    select case(mm)
      case(0)
        iopsymo(53:54) = -1
      case(1)
        iopsymo(53) = 10; iopsymo(54) = 2
      case(2)
        iopsymo(53) = 2;  iopsymo(54) = 10
    end select
  endif

! Axe 4 Oz
  mm = m
  do
    if( mm >= 0 ) exit
    mm = mm + 4
  end do
  mm = mod(mm,4)
  select case(mm)
    case(0)
      iopsymo(18:21:3) = 1
    case(1)
      iopsymo(18) = 3; iopsymo(21) = 9
    case(2)
      iopsymo(18:21:3) = -1
    case(3)
      iopsymo(18) = 9; iopsymo(21) = 3
  end select
! Axes S4
  if( mod(l+m,2) == 0 ) then
    iopsymo(28:31:3) = iopsymo(18:21:3)
  else
    select case(mm)
      case(0)
        iopsymo(28:31:3) = -1
      case(1)
        iopsymo(28) = 9; iopsymo(31) = 3
      case(2)
        iopsymo(28:31:3) = 1
      case(3)
        iopsymo(28) = 3; iopsymo(31) = 9
      end select
  endif

! Axe 6 Oz
  mm = m
  do
    if( mm >= 0 ) exit
    mm = mm + 6
  end do
  mm = mod(mm,6)
  select case(mm)
    case(0)
      iopsymo(51:52) = 1
    case(1)
      iopsymo(51) = 2; iopsymo(52) = 10
    case(2)
      iopsymo(51) = 4; iopsymo(52) = 8
    case(3)
      iopsymo(51:52) = -1
    case(4)
      iopsymo(51) = 8; iopsymo(52) = 4
    case(5)
      iopsymo(51) = 10; iopsymo(52) = 2
  end select
! Axes 6 negatif
  if( mod(l+m,2) == 0 ) then
    iopsymo(55:56) = iopsymo(51:52)
  else
    select case(mm)
      case(0,3)
        iopsymo(55:56) = - iopsymo(51:52)
      case default
        do is = 51,52
          if( iopsymo(is) < 6 ) then
            iopsymo(is+4) = iopsymo(is) + 6
          else
            iopsymo(is+4) = iopsymo(is) - 6
          endif
        end do
    end select
  endif

! Centrosymmetie
  if( lpair ) then
    iopsymo(25) = 1
  else
    iopsymo(25) = - 1
  endif

! Plans de symetrie
  m2 = mod(im,2)

  if( (m > 0 .and. m2 == 1) .or. (m < 0 .and. m2 == 0) ) then
    iopsymo(40) = -1
  else
    iopsymo(40) = 1
  endif
  if( m < 0 ) then
    iopsymo(41) = -1
  else
    iopsymo(41) = 1
  endif
  if( mod(l+m,2) == 1 ) then
    iopsymo(42) = -1
  else
    iopsymo(42) = 1
  endif

! Plans diagonaux contenant Oz
  m4 = mod(abs(m),4)
  if( (m4 == 0 .and. m >= 0) .or. (m4 == 2 .and. m < 0) ) then
    iopsymo(45:48:3) = 1
  elseif((m4 == 0 .and. m < 0) .or. (m4 == 2 .and. m >= 0)) then
    iopsymo(45:48:3) = -1
  else
    iopsymo(45:48:3) = 0
  endif

! Quand la valeur est differente de 1 ou - 1, le chiffre represente n
! pour le caractere egal a exp(n*i*pi/12).
  do is = 1,nopsm
    select case( iopsymo(is) )
      case(0)
        kopsymo(is) = (0._db, 0._db)
      case(1)
        kopsymo(is) = (1._db, 0._db)
      case(-1)
        kopsymo(is) = (-1._db, 0._db)
      case default
        kopsymo(is) = exp( iopsymo(is) * img * pi / 6 )
    end select
  end do

  return
end

!***********************************************************************

function natome_cal(igrpt_nomag,iopsymr,mpirank0,natomeq,natomp,pos)

  use declarations
  implicit none

  integer:: ia, igrpt_nomag, isym, mpirank0, natome, natome_cal, natomeq, natomp

  integer, dimension(nopsm):: iopsymr

  real(kind=db), dimension(3):: dp, ps
  real(kind=db), dimension(3,natomp):: pos

  natome = 0
  do ia = 1,natomeq
    ps(:) = pos(:,ia)
    call posequiv(mpirank0,ps,iopsymr,isym,igrpt_nomag)
    dp(:) = abs( ps(:) - pos(:,ia) )
    if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos ) cycle
    natome = natome + 1
  end do
  natome_cal = natome

  return
end

!***********************************************************************

! Selection des atomes du petit agregat
! Evaluation de leur groupe ponctuel dans cet agregat.

subroutine Atom_selec(adimp,Atom_axe,Atom_with_axe,Atom_nonsph,Atom_occ_mat,Axe_atom_clu,Base_ortho,dcosxyz, &
           dista,distai,Full_atom,Green,hubbard,i_self,ia_eq,ia_eq_inv,ia_rep,iaabs,iaabsi,iaproto,iaprotoi,icheck,igreq, &
           igroup,igroupi,igrpt_nomag,igrpt0,iopsym_atom,iopsymr,iord,is_eq,itype,itypei,itypep,itypepr,magnetic,m_hubb, &
           m_hubb_e,mpirank,natome,n_atom_0_self,n_atom_ind_self,n_atom_proto,natomeq,natomp,nb_eq,nb_rpr, &
           nb_rep_t,nb_sym_op,neqm,ngroup,ngroup_hubb,ngroup_m,nlat,nlatm,nspin,nspinp,ntype,numat,nx,occ_hubb_e,overad,popats, &
           pos,posi,rmt,rot_atom,roverad,rsort,rsorte,Spinorbite,Symmol,V_hubb,V_hubbard,Ylm_comp)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  character(len=8):: ptgrname_int

  integer, dimension(nopsm):: iopsyma, iopsymr
  integer, dimension(nopsm,natome):: iopsym_atom
  integer, dimension(natome):: iaprotoi, iatomp, igroupi, igrpt, itypei, nb_eq
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(natomeq):: ia_eq_inv
  integer, dimension(nb_sym_op,natome):: ia_eq, ia_eq_sym, is_eq
  integer, dimension(nb_sym_op,natome,natome):: ia_rep, nb_rep_t
  integer, dimension(natome,natome):: nb_rpr
  integer, dimension(natomp):: iaproto, igroup, itypep
  integer, dimension(ngroup):: itype
  integer, dimension(0:ntype):: nlat, numat
  integer, dimension(0:n_atom_proto):: itypepr

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb

  logical:: Atom_comp_cal, Atom_mag_cal, Atom_nonsph, Atom_occ_mat, Base_ortho, Full_atom, Green, Hubbard, magnetic, overad, &
          Spinorbite, Symmol, Ylm_comp
  logical, dimension(nb_sym_op):: Fait
  logical, dimension(natome):: Atom_comp, Atom_axe
  logical, dimension(0:natome):: Atom_mag
  logical, dimension(0:ngroup_m):: Atom_with_axe

  real(kind=db), dimension(3):: dcosxyz, dp, ps, v
  real(kind=db), dimension(3,3):: matopsym, rot_a
  real(kind=db), dimension(0:ntype) :: rmt, V_hubbard
  real(kind=db), dimension(3,natome):: Axe_atom_clui, posi
  real(kind=db), dimension(natome):: distai
  real(kind=db), dimension(natomp):: dista
  real(kind=db), dimension(3,natomp):: Axe_atom_clu, pos
  real(kind=db), dimension(3,3,natome):: rot_atom
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats
  real(kind=db), dimension(-m_hubb_e:m_hubb_e,-m_hubb_e:m_hubb_e,nspin,ngroup_hubb):: occ_hubb_e

  common/iopsym_abs/ iopsym_abs(nopsm)

  if( icheck > 0) write(3,110)

! Selection des atomes ou on effectue un developpement en harmoniques
! spheriques.

  iaabsi = 0

  ib = 0
  do ia = 1,natomeq
    ps(:) = pos(:,ia)
    call posequiv(mpirank,ps,iopsymr,isym,igrpt_nomag)
    dp(:) = abs( ps(:) - pos(:,ia) )
    if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos ) cycle
    ib = ib + 1
    iatomp(ib) = ia
    iaprotoi(ib) = iaproto(ia)
    itypei(ib) = itypep(ia)
    igroupi(ib) = igroup(ia)
    posi(:,ib) = pos(:,ia)
    Axe_atom_clui(:,ib) = Axe_atom_clu(:,ia)
    distai(ib) = dista(ia)
    if( ia == iaabs ) iaabsi = ib
  end do
  if( iaabsi == 0 ) then
    ps(:) = pos(:,iaabs)
    call posequiv(mpirank,ps,iopsymr,isym,igrpt_nomag)
    do ia = 1,natome
      dp(:) = abs( ps(:) - posi(:,ia) )
      if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos ) cycle
      iaabsi = ia
      exit
    end do
  endif

  if( Atom_occ_mat .and. hubbard .and. i_self == 1 ) then
    V_hubb(:,:,:,:,:) = ( 0._db, 0._db )
    do iapr = n_atom_0_self,n_atom_ind_self
      if( Full_atom ) then
        igr = igroupi(iapr)
        it = itypei(iapr)
      else
        igr = igreq(iapr,1)
        ipr = iaprotoi(igr)
        it = itypepr(ipr)
      end if
      l = l_hubbard( numat(it) )
      do m = -l,l
        do mp = -l,l
          do isp = 1,nspin
            if( m == mp ) V_hubb(m,mp,isp,isp,iapr) = cmplx( 0.5_db * V_hubbard(it), 0._db )
            V_hubb(m,mp,isp,isp,iapr) = V_hubb(m,mp,isp,isp,iapr) - cmplx( V_hubbard(it) * occ_hubb_e(m,mp,isp,igr), 0._db )
            if( nspinp > nspin ) V_hubb(-mp,-m,nspinp,nspinp,iapr) = V_hubb(m,mp,isp,isp,iapr)
          end do
        end do
      end do
    end do
  end if

! Elaboration de la liste d'atomes equivalents
! On met en premier l'atome representant.
  nb_eq(:) = 0
  do ia = 1,natomeq
    ps(:) = pos(:,ia)
    call posequiv(mpirank,ps,iopsymr,isym,igrpt_nomag)
    if( abs(isym) /= 1 ) cycle
    do ib = 1,natome
      dp(:) = abs( ps(:) - posi(:,ib) )
      if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos .or. itypep(ia) /= itypei(ib) ) cycle
      nb_eq(ib) = nb_eq(ib) + 1
      ia_eq( nb_eq(ib), ib ) = ia
      is_eq( nb_eq(ib), ib ) = isym
      exit
    end do
  end do
  do ia = 1,natomeq
    ps(:) = pos(:,ia)
    call posequiv(mpirank,ps,iopsymr,isym,igrpt_nomag)
    if( abs(isym) == 1 ) cycle
    do ib = 1,natome
      dp(:) = abs( ps(:) - posi(:,ib) )
      if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos .or. itypep(ia) /= itypei(ib) ) cycle
      nb_eq(ib) = nb_eq(ib) + 1
      ia_eq( nb_eq(ib), ib ) = ia
      is_eq( nb_eq(ib), ib ) = isym
      exit
    end do
  end do

  do ib = 1,natome
    do i = 1,nb_eq(ib)
      ia_eq_inv( ia_eq( i, ib ) ) = ib
    end do
  end do

  do ia = 1,natome
    js = 0
    do is = 1,nopsm
      if( iopsymr(is) == 0 ) cycle
      js = js + 1
      call opsym(is,matopsym)
      ps(:) = posi(:,ia)
      v = matmul( matopsym, ps )
      do i = 1,nb_eq(ia)
        ie = ia_eq(i,ia)
        dp(:) = abs( v(:) - pos(:,ie) )
        if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos .or. itypei(ia) /= itypep(ie) ) cycle
        ia_eq_sym(js,ia) = i
        exit
      end do
    end do
  end do

! Calcul du nombre d'atomes n' tels que S(n0) = n0 et S(n1) = n'
! et evaluation des n1
  nb_rpr(:,:) = 0
  nb_rep_t(:,:,:) = 0
  do ia = 1,natome
    do ib = 1,natome

      Fait(:) = .false.
      ind_rep = 0

      do i = 1,nb_eq(ib)

!          Fait(:) = .false.
        ie = ia_eq(i,ib)
        if( Fait(i) ) cycle  ! i est deja represente ou representant

        ind_rep = ind_rep + 1
        ia_rep(ind_rep,ia,ib) = ie

        do is = 1,nopsm
          if( iopsymr(is) == 0 ) cycle
          call opsym(is,matopsym)
          ps(:) = posi(:,ia)
          v = matmul( matopsym, ps )
          dp(:) = abs( v(:) - ps(:) )
          if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos ) cycle         ! S(ia) /= ia

          ps(:) = pos(:,ie)
          ps = matmul( matopsym, ps )

          do j = 1,nb_eq(ib)   ! recherche du n' = S(n1)
            if( Fait(j) ) cycle
            ig = ia_eq(j,ib)
            dp(:) = abs( pos(:,ig) - ps(:) )
            if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos  .or. itypep(ie) /= itypep(ig) ) cycle
            nb_rep_t(ind_rep,ia,ib) = nb_rep_t(ind_rep,ia,ib) + 1
            Fait(j) = .true.
          end do

        end do

      end do

      nb_rpr(ia,ib) = ind_rep
    end do
  end do

  if( magnetic ) then
    Atom_mag(0) = Atom_mag_cal(igrpt0)
  else
    Atom_mag(0) = .false.
  endif

! Evaluation de la symetrie locale
  do ia = 1,natome
    if( icheck > 1 ) write(3,130) ia
    ps(:) = posi(:,ia)
    call point_group_atom(Atom_comp(ia),Atom_mag(ia),Atom_with_axe, Atom_nonsph,Axe_atom_clu,Base_ortho,dcosxyz, &
        iaabs,iatomp(ia),icheck,igroup,igroupi(ia),igrpt(ia),igrpt0,iopsyma,iopsymr,itype,itypep,magnetic,mpirank,natomp, &
        ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,pos,ps,rot_a,Spinorbite,Symmol)

    Atom_comp(ia) = Ylm_comp

    if( magnetic ) then
      Axe_atom_clui(:,ia) = rot_a(:,3)
    else
      Atom_mag(ia) = .false.
    endif
    iopsym_atom(:,ia) = iopsyma(:)
    rot_atom(:,:,ia) = rot_a(:,:)
  end do

! Atomes se trouvant sur l'axe de l'agregat
  do ia = 1,natome
    if( sum( abs( posi(1:2,ia) ) ) < eps10 ) then
      Atom_axe(ia) = .true.
    else
      Atom_axe(ia) = .false.
    endif
  end do

! Parametres du maillage
  if( Green ) then
    rsort = min( rsorte, distai(natome) ) + epspos
  else
    it = itypei(natome)
    if( overad ) then
      rm = distai(natome) + rmt(it) + roverad
    else
      rm = distai(natome) + rmt(it) + adimp
    endif
    rsort = max( rsorte, rm + adimp )
  endif
  rsort = max( rsort, rmt(itypei(1)) + adimp )
  rmax = ( rsort / adimp + sqrt(2._db) * ( iord / 2 ) + epspos ) / cos( pi / 6 )
  nx = nint( rmax )

  if( icheck > 0) then
    write(3,140) rsort*bohr
    write(3,150) nx
    write(3,160) natome, igrpt0, Ylm_comp, Atom_mag(0)
    if( Full_atom ) then
      write(3,'(A)') ' Full_atom mode'
    else
      write(3,'(A)') ' No Full_atom mode'
    endif
    if( Magnetic .or. Atom_nonsph ) then
      write(3,170)
    else
      write(3,180)
    endif
    do ia = 1,natome
      it = itypei(ia)
      if( atom_mag(ia) ) then
        write(3,190) ia, numat(it), it, igroupi(ia), iaprotoi(ia), iatomp(ia), posi(1:3,ia)*bohr, igrpt(ia), &
              ptgrname_int(igrpt(ia)), Ylm_comp, Atom_axe(ia), Atom_mag(ia), Axe_atom_clui(:,ia)
      else
        write(3,190) ia, numat(it), it, igroupi(ia), iaprotoi(ia), iatomp(ia), posi(1:3,ia)*bohr, igrpt(ia), &
              ptgrname_int(igrpt(ia)), Ylm_comp, Atom_axe(ia), Atom_mag(ia)
      endif
    end do
    write(3,'(/A)') ' Atom rotation matrices :'
    na_ligne = 4
    do iga = 1,(natome-1)/na_ligne + 1
      na1 = na_ligne * iga - na_ligne + 1
      na2 = min(na1+na_ligne-1,natome)
      write(3,200) (ia, ia = na1,na2)
      do i = 1,3
        write(3,210) ( rot_atom(i,:,ia), ia = na1,na2)
      end do
    end do

    write(3,220)
    do ia = 1,natome
      write(3,230) ia, ( ia_eq(i,ia), is_eq(i,ia), i = 1,nb_eq(ia) )
    end do

  endif

  if( icheck > 1) then
    write(3,232)
    do ia = 1,natomeq
      write(3,230) ia, ia_eq_inv(ia)
    end do
    write(3,240) ( ia, ia = 1,natome )
    js = 0
    do is = 1,nopsm
      if( iopsymr(is) == 0 ) cycle
      js = js + 1
      write(3,250) is, ( ia_eq_sym(js,ia), ia = 1,natome )
    end do
    write(3,260)
    do ia = 1,natome
      do ib = 1,natome
      write(3,270) ia, ib, ( ia_rep(i,ia,ib), nb_rep_t(i,ia,ib), i = 1,nb_rpr(ia,ib) )
      end do
    end do
  endif

! L'elaboration du reseau s'effectuant autour de chaque atome, en
! definissant la zone atomique par la sphere muffin-tin, il est
! necessaire d'avoir des rayons muffin-tin laissant l'espace a au moins
! un point entre les atomes.

  if( .not. Green ) then
    do ia1 = 1,natome
      it1 = itypei(ia1)
      do ia2 = ia1+1,natome
        it2 = itypei(ia2)
        v(:) = posi(:,ia1) - posi(:,ia2)
        dist = vnorme(Base_ortho,dcosxyz,v)
        if( dist < rmt(it1) + rmt(it2) - adimp .and. mpirank == 0) then
          call write_error
          do ipr = 3,9,3
            write(ipr,280) ia1, numat(itypei(ia1)), igroupi(ia1), ia2, numat(itypei(ia2)), igroupi(ia2), dist*bohr, &
                   rmt(it1)*bohr, rmt(it2)*bohr, adimp*bohr
          end do
          stop
        endif
      end do
    end do
  endif

  if( iaabsi == 0 .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,290)
    end do
    stop
  endif

  return
  110 format(/' ---- Atom_selec --',100('-'))
  130 format(/' Atome',i4,', with working cluster sub-point group')
  140 format(/'  Rsort = ',f8.3,' A')
  150 format('  nx =',i4)
  160 format('  natome =',i4,', igrpt =',i4,', Cluster_comp =',l2,', Cluster_mag =',l2)
  170 format(/'  ia   Z  it  igr ipr iap     posx      posy      posz   igrpt PtGrName Comp  Axe  Mag',9x,'Axe_atom')
  180 format(/'  ia   Z  it  igr ipr iap     posx      posy      posz   igrpt PtGrName Comp  Axe  Mag')
  190 format(i4,2i4,i5,2i4,3f10.5,i5,6x,a8,l1,2l5,3f10.5)
  200 format(/'  Atom :',4x,i4,3(21x,i4))
  210 format(4(3x,3f7.3))
  220 format(/' Atom  ia_eq  is_eq')
  230 format(i4,3x,48(1x,2i4))
  232 format(/' Atom  ia_eq_inv')
  240 format(/' ia_sym_eq (index of the representative atom)', /' Sym \ ia ',1x,48i4)
  250 format(i4,7x,48i4)
  260 format(/' Atom Atom ia_rep, nb_rep_t ')
  270 format(i4,i5,3x,96i4)
  280 format(//' When using the FDM mode, the minimum distance between the atoms is the sum',/ &
               ' of their Rmt radius plus the interpoint distance.',// &
               '   The atoms',i4,' (Z =',i3,', igr =',i4,') and',i4, ' (Z =',i3,', igr =',i4,') are too close.',// &
               '   Their interatomic distance is',f8.5,' Angstroem.',/ &
               '   Their Rmt radius are ',f8.5,' and ',f8.5,' Angstroem.',/ &
               '   The interpoint distance is ',f8.5,' Angstroem.')
  290 format(//' The absorbing atom is outside the calculation sphere !' //)
end

!***********************************************************************

! Cherche le point equivalent a l'interieur de la zone de calcul
! compte tenu des symmetries.

subroutine posequiv(mpirank,pos,iopsym,isym,igrpt_nomag)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(nopsm):: iopsym

  logical:: d6, d6p, dpn12, dp13, d3p12, d3pn12, Found, gdiag, symok
  logical, dimension(3):: pnul, ppn, ppos
  logical, dimension(2):: p3nul, p3pn, p3pos

  real(kind=db), dimension(3):: p, pos, v
  real(kind=db), dimension(3,3):: matopsym

  isym = 1

  do k = 1,3
    pnul(k) = abs(pos(k)) < epspos
  end do
  if( pnul(1) .and. pnul(2) .and. pnul(3) ) return

  boucle: do is = 1,nopsm
    if( iopsym(is) == 0 ) cycle

    call opsym(is,matopsym)
    v = matmul( matopsym, pos )
    do k = 1,3
      pnul(k) = abs(v(k)) < epspos
      ppos(k) = v(k) > epspos
      ppn(k) = ppos(k) .or. pnul(k)
    end do
    dpn12 = v(1) > v(2) - epspos
    dp13 = v(1) > v(3) + epspos
    gdiag = ppn(1) .and. ( abs(v(1)-v(2)) < epspos ) .and. ( abs(v(1)-v(3)) < epspos )

! Cas des groupes trigo et hexagonal

    d3p12 = .false.; d3pn12 = .false.; d6 = .false.; d6p = .false.
    if( igrpt_nomag > 15 .and. igrpt_nomag < 28 ) then
      r3 = sqrt(3._db)
      f = 1 / r3
      p(1) = v(1) + f * v(2)
      p(2) = 2 * f * v(2)
      p(3) = v(3)
      do k = 1,2
        p3nul(k) = abs(p(k)) < epspos
        p3pos(k) = p(k) > epspos
        p3pn(k) = p3pos(k) .or. p3nul(k)
      end do
      d3p12 = (p(1) > p(2) + epspos) .or. ( pnul(1) .and. pnul(2) )
      d3pn12 = p(1) > p(2) - epspos
      d6 = v(1) > abs( r3 * v(2) ) - epspos
      d6p = ( ( v(1) > - sqrt(3._db) * v(2) + epspos ) .and. (  v(1) > r3 * v(2) - epspos ) ) .or. &
             ( pnul(1) .and. pnul(2) )
    endif

    Found = .false.

    select case(igrpt_nomag)

      case(1)
        Found = .true.

      case(2)
        if( ppos(3) .or. ( pnul(3) .and.  ppos(1) ) .or. ( pnul(3) .and.  pnul(1) .and. ppn(2) ) ) Found = .true.

      case(3)
        if( ppn(3) ) Found = .true.

      case(4)
        if( ppos(1) .or. ( pnul(1) .and.  ppn(2) ) ) Found = .true.

      case(5)
        if( ppn(3) .and. ( ppos(1) .or. ( pnul(1) .and.  ppn(2) ) ) ) Found = .true.

      case(6)
        if( ppn(1) .and. ppn(2) ) Found = .true.

      case(7)
        if( ( ppos(1) .and. ppos(2) ) .or. ( pnul(1) .and. ppn(2) .and. ppn(3) ) .or. &
            ( pnul(2) .and. ppn(1) .and. ppn(3) ) ) Found = .true.

      case(8)
        if( ppn(1) .and. ppn(2) .and. ppn(3) ) Found = .true.

      case(9)
        if( ( ppos(1) .and. ppn(2) ) .or. ( pnul(1) .and. pnul(2) ) ) Found = .true.

      case(10)
        if( ( ppos(1) .and. ppn(2) ) .or. ( pnul(1) .and. pnul(2) .and. ppn(3) ) ) Found = .true.

      case(11)
        if( ( ppn(3) .and. ppos(1) .and. ppn(2) ) .or. ( ppn(3) .and. pnul(1) .and. pnul(2) ) ) Found = .true.

      case(12)
        if( ppn(2) .and. dpn12 ) Found = .true.

      case(13)
        if( ( ppos(2) .and. dpn12 ) .or. ( pnul(2) .and. ppn(1) .and. ppn(3) ) ) Found = .true.

      case(14)
        if( ( ppos(2) .and. dpn12 ) .or. ( pnul(2) .and. dpn12  .and. ppn(3) ) ) Found = .true.

      case(15)
        if( ppn(2) .and. dpn12  .and. ppn(3)) Found = .true.

      case(16)  ! C3
        if( ( p3pos(1) .and. p3pn(2) ) .or. ( p3nul(1) .and. p3nul(2) ) ) Found = .true.

      case(17)  ! S6
        if( ( p3pos(1) .and. p3pn(2) .and. ppos(3) ) .or. ( p3nul(1) .and. p3nul(2).and. ppn(3) ) .or. &
            ( d3p12 .and. p3pn(1) .and. p3pn(2) .and. pnul(3) ) ) Found = .true.

      case(18)  ! C3v
        if( d6 ) Found = .true.

      case(19)
        if( ( (  ( p3pos(1) .and. p3pn(2) ) .or. ( pnul(1) .and. pnul(2) )  ) .and. ppos(3) )  .or. &
           ( pnul(3) .and. p3pos(1) .and. d3pn12   )) Found = .true.

      case(20)
        if( d6p ) Found = .true.

      case(21)  ! C3h
        if( ( p3pos(1) .and. p3pn(2) .and. ppn(3) ) .or. ( p3nul(1) .and. p3nul(2).and. ppn(3) ) ) Found = .true.

      case(22)   ! C6
        if( p3pn(1) .and. p3pn(2) .and. d3p12 ) Found = .true.

      case(23)   ! C6h
        if( p3pn(1) .and. p3pn(2) .and. d3p12 .and. ppn(3) ) Found = .true.

      case(24)
        if( d6 .and. ppn(3) ) Found = .true.

      case(25)
        if( d6p .and. ppn(2) ) Found = .true.

      case(26)
        if( ( ppn(2) .and. d3p12 .and. ppos(3) ).or. ( d6 .and. ppn(2) .and. pnul(3) ) ) Found = .true.

      case(27)
        if( d6 .and. ppn(2) .and. ppn(3) ) Found = .true.

      case(28)
        a3 = abs( v(3) )
        if( ppn(2) .and. (( dpn12 .and. v(1) > a3+epspos ) .or. ( abs(v(1)-v(2)) < epspos .and. abs(v(1)-a3) < epspos ) )) &
                                                      Found = .true.

      case(29) ! m3
        if( ppn(2) .and. ppn(3) .and. ( ( dpn12 .and. dp13 ) .or. gdiag ) ) Found = .true.

      case(30) ! -43m
        a3 = abs( v(3) )
        if( dpn12 .and. v(2) > a3-epspos ) Found = .true.

      case(31) ! 432
        if( ( ppn(3) .and. ppos(2) .and. dpn12 .and. dp13 ) .or. ( ppn(1) .and. pnul(2) .and. pnul(3) ) &
             .or. gdiag ) Found = .true.

      case(32) ! m3m
        a3 = abs( v(3) )
        if( ppn(3) .and. ( dpn12 .and. v(2) > a3-epspos ) ) Found = .true.

      case(35)
        if( ppn(1) .and. ( ppos(3) .or. ( pnul(3) .and. ppn(2) ) ) ) Found = .true.

      case default

        symok = .false.

        select case(is)

          case(10)        ! Axe 2 diagonal (110)
            if( ppos(3) .or. (pnul(3) .and. dpn12) ) symok = .true.

          case(11)        ! Axe 2 diagonal (1-10)
            if( ppos(3) .or. ( pnul(3) .and. pos(2) < - pos(1) + epspos ) ) symok = .true.

          case(22,23,24)  ! Axes 2
            k = is - 21
            k2 = 1 + mod(k,3)
            k3 = 1 + mod(k+1,3)
            if( ppos(k2) .or. ( pnul(k2) .and. ppn(k3) ) ) symok = .true.

          case(25)        ! Centrosymetrie
            if( ppos(3) .or. ( pnul(3) .and.  ppos(1) ) .or. ( pnul(3) .and.  pnul(1) .and. ppn(2) ) ) symok = .true.

          case(40,41,42)  ! Plans de symetrie principaux
            k = is - 39
            if( ppn(k) ) symok = .true.

          case(45)        ! Plan diagonal
            if( dpn12 ) symok = .true.

          case(48)        ! Plan diagonal
            if( pos(2) < -pos(1)+epspos ) symok = .true.

        end select

      if( symok ) pos = v

    end select

    if( Found ) exit boucle

  end do boucle

  if( is > nopsm .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,120) igrpt_nomag, pos(1:3)
      call write_iopsym(iopsym,ipr)
    end do
    stop
  endif

  pos = v

  isym = iopsym(is) * isym_inv(is)

  return
  120 format(/' Probleme dans posequiv : symetrie mal programmee !'/ ' igrpt_nomag =',i4,', pos =',3f7.3,//' iopsym =')
end

!***********************************************************************

! Determination de l'operation de symmetrie inverse.

function isym_inv(is)

  select case(is)
    case(2,4,6,8,32,34,36,38,49,51,53,55)
      isym_inv = is + 1
    case(3,5,7,9,33,35,37,39,50,52,54,56)
      isym_inv = is - 1
    case(16,17,18,26,27,28)
      isym_inv = is + 3
    case(19,20,21,29,30,31)
      isym_inv = is - 3
    case default
      isym_inv = is
  end select

  return
end

!***********************************************************************

subroutine point_group_atom(Ylm_comp,Atom_mag,Atom_with_axe,Atom_nonsph,Axe_atom_clu,Base_ortho,dcosxyz, &
          iaabs,iap,icheck,igroup,igra,igrpt,igrpt0,iopsym,iopsymr,itype,itypep,magnetic,mpirank,natomp, &
          ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,pos,posi,rot_atom,Spinorbite,Symmol)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(natomp):: igroup, itypep
  integer, dimension(nopsm):: iopsym, iopsymr, irotiops
  integer, dimension(ngroup):: itype
  integer, dimension(0:ntype):: nlat, numat

  logical:: Ylm_comp, Atom_comp_cal, Atom_mag, Atom_mag_cal, Atom_nonsph, Base_ortho, Magnetic, Spinorbite, Symmol
  logical, dimension(0:ngroup_m):: Atom_with_axe

  real(kind=db), dimension(3):: dcosxyz, posi, px, py, pz, v
  real(kind=db), dimension(3,3):: rot_atom, rot_tem, rot_un
  real(kind=db), dimension(3,natomp):: Axe_atom_clu, Axe_atom_clut, p, pos
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats

  rot_un = 0._db
  do i = 1,3
    rot_un(i,i) = 1._db
  end do

  rot_atom(:,:) = rot_un(:,:)
  Axe_atom_clut(:,1:natomp) = Axe_atom_clu(:,1:natomp)

  if( abs( posi(1) ) < eps10 .and. abs( posi(2) ) < eps10 .and. abs( posi(3) ) < eps10 ) then

    iopsym = iopsymr(:)
    igrpt = igrpt0

  else

    do ib = 1,natomp
      p(:,ib) = pos(:,ib) - posi(:)
    end do

! Determination des symmetries du groupe ponctuel (iopsymc)
    do i = 1,2

      call sym_cluster(Atom_with_axe,Axe_atom_clut,Base_ortho,dcosxyz,iaabs,igroup,iopsym,itype,itypep, &
             natomp,ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,popats,p,Symmol)

! Le groupe ponctuel de l'atome doit etre un sous-groupe du groupe de l'agregat.
      rot_tem = transpose(rot_atom)
      call iop_rot(icheck,irotiops,rot_tem)
      where( irotiops /= 0 ) irotiops = iopsymr(irotiops)
      where( irotiops == 0 ) iopsym = 0

      call cluster_rot(iopsym,rot_tem)

      if( sum( abs( rot_tem(:,:) - rot_un(:,:) ) ) < eps10  ) exit

      rot_atom(:,:) = matmul( rot_tem, rot_atom )

      do ib = 1,natomp
        v(:) = p(:,ib)
        v = matmul( rot_tem, v )
        p(:,ib) = v(:)
      end do
      if( magnetic .or. Atom_nonsph ) then
        do ib = 1,natomp
          v(:) = Axe_atom_clut(:,ib)
          v = matmul( rot_tem, v )
          Axe_atom_clut(:,ib) = v(:)
        end do
      endif

    end do

    call numgrpt(iopsym,igrpt,igrpt_nomag,mpirank,.false.)

  endif

  if( magnetic ) then
    Atom_mag = Atom_mag_cal(igrpt)
  else
    Atom_mag = .false.
  endif

  if( Spinorbite .and. Magnetic ) then

    if( Atom_mag .and. .not. Atom_with_axe(igra) ) then
      call Axe_mag_cal(igrpt,pz)
    else
      pz = Axe_atom_clut(:,iap)
    endif
    if( abs( pz(2) - 1._db ) < epspos ) then
      px(1:2) = 0._db
      px(3) = - 1._db
      call prodvec(py,pz,px)
    elseif( abs( pz(2) + 1._db ) < epspos ) then
      px(1:2) = 0._db
      px(3) = 1._db
      call prodvec(py,pz,px)
    else
      py(1:3:2) = 0._db
      py(2) = 1._db
      call prodvec(px,py,pz)
      pp = sqrt( sum( px(:)**2 ) )
      px(:) = px(:) / pp
      call prodvec(py,pz,px)
    endif

    rot_tem(1,1:3) = px(1:3)
    rot_tem(2,1:3) = py(1:3)
    rot_tem(3,1:3) = pz(1:3)
    rot_atom(:,:) = matmul( rot_tem, rot_atom )

    rot_tem = transpose(rot_tem)

    call iop_rot(icheck,irotiops,rot_tem)
    where( irotiops /= 0 ) iopsym = iopsym(irotiops)

  endif

  Ylm_comp = Atom_comp_cal(igrpt)

  if( icheck > 1 ) then
    if( Atom_mag ) then
      write(3,110) igrpt
    else
      write(3,120) igrpt
    endif
    write(3,130)
    call write_iopsym(iopsym,3)
    write(3,150) rot_atom(:,:)
  endif

  return
  110 format(/' Local point group number',i3,', magnetic atom.')
  120 format(/' Local point group number',i3,', non-magnetic atom.')
  130 format(/' iopsym =')
  150 format(/' Atomic rotation matrix :',3(/,3f10.5))
end

!***********************************************************************

! Pour certains groupes magnetique, l'atome central ne peut pas etre magnetique

function Atom_mag_cal(igrpt)

  implicit none

  integer:: igrpt
  logical:: Atom_mag_cal

  select case(igrpt)
    case(1,2,3,4,5,9,10,11,16,17,21,22,23,35,39,40,41,44,45,46,48, 53,58,60,67,70,73,78,84)
      Atom_mag_cal = .true.
    case default
      Atom_mag_cal = .false.
  end select

  return
end

!***********************************************************************

! Pour certains groupes magnetique, l'atome central ne peut pas etre magnetique

subroutine Axe_mag_cal(igrpt,Axe)

  use declarations
  real(kind=db), dimension(3):: Axe

  Axe(1:3) = 0._db

  select case(igrpt)
    case(34,35,41)
      Axe(2) = 1._db
    case default
      Axe(3) = 1._db
  end select

  return
end

!***********************************************************************

! Pour certains groupes la base commode est la base des Y(l,m) complexes

function Atom_comp_cal(igrpt)

  logical Atom_comp_cal

  select case(igrpt)
    case(9,10,11,16,17,21,22,23,47,51,52,55,56,57,68,69,75,76,77)
      Atom_comp_cal = .true.
    case default
      Atom_comp_cal = .false.
  end select

  return
end

!***********************************************************************

! Calcul du nombre de point du maillage.

subroutine nbpoint(Adimp,Base_hexa,Base_ortho,D_max_pot,dcosxyz,Green,iaabs,igrpt_nomag,iopsymr,iord,Moy_loc, &
               mpirank,natomp,npoint,npso,nvois,nx,pos,rsort)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(nopsm):: iopsymr

  logical Base_hexa, Base_ortho, Green, Moy_loc

  real(kind=db), dimension(3):: dcosxyz, p, v, w
  real(kind=db), dimension(3,natomp):: pos

  f = sqrt(3._db)/2

  if( Moy_loc ) then
    dcour = rsort
    do ia = 1,natomp
      if( ia == iaabs ) cycle
      p(:) = pos(:,ia) - pos(:,iaabs)
      dcour = min(dcour,vnorme(Base_ortho,dcosxyz,p))
    end do
    rmax = min(dcour,rsort) / adimp + epspos
  else
    rmax = rsort / adimp + sqrt(2._db) * iord / 2 + epspos
  endif

  p(:) = pos(:,iaabs) / adimp
  npoint = 0
  npso = 0
  do ix = -nx,nx
    do iy = -nx,nx
      if( Base_hexa ) then
        v(1) = ix - 0.5_db * iy
        v(2) = f * iy
      else
        v(1) = 1._db * ix
        v(2) = 1._db * iy
      endif
      do iz = -nx,nx
        v(3) = 1._db * iz
        if( Moy_loc ) then
          w(:) = v(:) - p(:)
          dist = vnorme(Base_ortho,dcosxyz,w)
        else
          dist = vnorme(Base_ortho,dcosxyz,v)
        endif
        if( dist > rmax ) cycle
        w(:) = v(:)
        call posequiv(mpirank,w,iopsymr,isym,igrpt_nomag)
        w(:) = abs( (v(:) - w(:) ) )
        if( w(1) > epspos .or. w(2) > epspos .or. w(3) > epspos ) cycle
        if( Green .and. .not. Moy_loc ) then
! Pour la moyenne du potentiel, on ne prend pas les points loin de tout
! atome
          do ia = 1,natomp
            w(:) = pos(:,ia) / adimp - v(:)
            dista = vnorme(Base_ortho,dcosxyz,w)
            if( dista < D_max_Pot / adimp ) exit
          end do
          if( ia > natomp ) cycle
        endif
        npso = npso + 1
        if( dist*adimp <= rsort + epspos ) npoint = npoint + 1
     end do
    end do
  end do

  if( Base_hexa ) then
    nvois = 4 * iord
  else
    nvois = 3 * iord
  endif

  return
end

!***********************************************************************

! Elaboration du maillage.

subroutine reseau(Adimp,Base_hexa,Base_ortho,D_max_pot,dcosxyz,Green,iaabs,icheck, &
               igrpt_nomag,indice,iopsymr,iord,itypei,Moy_loc,mpirank,mpres,natome,natomp,nim, &
               npoint,npr,npso,npsom,ntype,numia,nx,pos,posi,rmt,rsort,rvol,xyz)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(natome):: itypei
  integer, dimension(nopsm):: iopsymr
  integer, dimension(npsom,3):: indice
  integer, dimension(npsom):: numia
  integer, dimension(-nx:nx,-nx:nx,-nx:nx):: mpres

  logical:: Base_hexa, Base_ortho, Green, Moy_loc

  real(kind=db), dimension(3):: dcosxyz, p, x, v, w, w1
  real(kind=db), dimension(nim):: rvol
  real(kind=db), dimension(0:ntype) :: rmt
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(4,npsom):: xyz

  if( icheck > 0 ) write(3,100)

  mpres(:,:,:) = 0
  f = sqrt(3._db) / 2

  if( Moy_loc ) then
    dcour = rsort
    do ia = 1,natomp
      if( ia == iaabs ) cycle
      p(:) = pos(:,ia) - pos(:,iaabs)
      dcour = min(dcour,vnorme(Base_ortho,dcosxyz,p))
    end do
    rmax = min(dcour,rsort) / adimp + epspos
  else
    rmax = rsort / adimp + sqrt(2._db) * iord / 2 + epspos
  endif

  p(:) = pos(:,iaabs) / adimp
  i = 0
  do ix = -nx,nx
    do iy = -nx,nx
      if( base_hexa ) then
        v(1) = ix - 0.5_db * iy
        v(2) = f * iy
      else
        v(1) = 1._db * ix
        v(2) = 1._db * iy
      endif
      do iz = -nx,nx
        v(3) = 1._db * iz
        if( Moy_loc ) then
          w(:) = v(:) - p(:)
          dist = vnorme(Base_ortho,dcosxyz,w)
        else
          dist = vnorme(Base_ortho,dcosxyz,v)
        endif
        if( dist > rmax ) cycle
        w(:) = v(:)
        call posequiv(mpirank,w,iopsymr,isym,igrpt_nomag)
        if( abs(v(1)-w(1)) > epspos .or. abs(v(2)-w(2)) > epspos .or. abs(v(3)-w(3)) > epspos ) cycle
        if( Green .and. .not. Moy_loc ) then
! Pour la moyenne du potentiel, on ne prend pas les points loin de tout
! atome
          do ia = 1,natomp
            w(:) = pos(:,ia) / adimp - v(:)
            dista = vnorme(Base_ortho,dcosxyz,w)
            if( dista < D_max_Pot / adimp ) exit
          end do
          if( ia > natomp ) cycle
        endif
        i = i + 1
        indice(i,1) = ix
        indice(i,2) = iy
        indice(i,3) = iz
        xyz(4,i) = adimp * dist
      end do
    end do
  end do
  npso = i

! Mise en ordre par rapport a la distance au centre de l'agregat.
  do i = 1,npso
    do j = i+1,npso
      if( xyz(4,i) < xyz(4,j) + epspos ) cycle
      xyzij = xyz(4,i)
      xyz(4,i) = xyz(4,j)
      xyz(4,j) = xyzij
      indi1 = indice(i,1)
      indi2 = indice(i,2)
      indi3 = indice(i,3)
      indice(i,1:3) = indice(j,1:3)
      indice(j,1) = indi1
      indice(j,2) = indi2
      indice(j,3) = indi3
    end do
  end do

  do i = 1,npso
    if( base_hexa ) then
      xyz(1,i) = ( indice(i,1) - 0.5_db * indice(i,2) ) * adimp
      xyz(2,i) = f * indice(i,2) * adimp
      xyz(3,i) = indice(i,3) * adimp
    else
      xyz(1:3,i) = indice(i,1:3) * adimp
    endif
  end do

! Nombre de points a l'interieur de la zone de calcul :
  do i = 1,npso
    if( xyz(4,i) > rsort + epspos ) exit
  end do
  npoint = i - 1
  numia(1:npoint) = 0
  if( npoint < npso ) numia(npoint+1:npso) = -2

  do i = 1,npso
    mpres(indice(i,1),indice(i,2),indice(i,3)) = i
  end do

! Liste des points qui sont dans un atome.
  do ia = 1,natome
    it = itypei(ia)
    p(1:3) = posi(1:3,ia)
    do i = 1,npoint
      x(1:3) = xyz(1:3,i)
      call posrel(Base_ortho,dcosxyz,iopsymr,x,p,v,dr,isym)
      if( dr < rmt(it) + epspos ) numia(i) = ia
    end do
  end do
  npr = npoint
  do i = 1,npoint
    if( numia(i) /= 0 ) npr = npr - 1
  end do

! Calcul du volume des points.

  f3 = 1 / sqrt(3._db)
  df3 = 2 * f3
  istop = 0
  p(:) = pos(:,iaabs) / adimp
  rvol(:) = 0._db
  boucle_ix: do ix = -nx,nx
    do iy = -nx,nx
      if( base_hexa ) then
        v(1) = ix - 0.5_db * iy
        v(2) = f * iy
      else
        v(1) = 1._db * ix
        v(2) = 1._db * iy
      endif
      do iz = -nx,nx
        v(3) = 1._db * iz
        if( Moy_loc ) then
          w(:) = v(:) - p(:)
          dist = vnorme(Base_ortho,dcosxyz,w)
        else
          dist = vnorme(Base_ortho,dcosxyz,v)
        endif
        if( dist > rmax ) cycle
        w(:) = v(:)
        call posequiv(mpirank,w,iopsymr,isym,igrpt_nomag)

        if( Green .and. .not. Moy_loc ) then
! Pour la moyenne du potentiel, on ne prend pas les points loin de tout
! atome
          do ia = 1,natomp
            w1(:) = pos(:,ia) / adimp - w(:)
            dista = vnorme(Base_ortho,dcosxyz,w1)
            if( dista < D_max_Pot /adimp ) exit
          end do
          if( ia > natomp ) cycle
        endif

        if( base_hexa ) then
          ix1 = nint( w(1) + f3 * w(2) )
          iy1 = nint( df3 * w(2) )
        else
          ix1 = nint( w(1) )
          iy1 = nint( w(2) )
        endif
        iz1 = nint( w(3) )
        i = mpres(ix1,iy1,iz1)
        if( i > npoint ) cycle
        if( i == 0 .and. mpirank == 0 ) then
          call write_error
          do ipr = 3,9,3
            write(ipr,130) ix, iy, iz, ix1, iy1, iz1
          end do
          istop = 1
          exit boucle_ix
        endif
        rvol(i) = rvol(i) + 1
      end do
    end do
  end do boucle_ix

  if( istop == 0 ) then
    rvmax = 1._db
    do i = 1,npoint
      rvmax = max(rvmax,rvol(i))
    end do
    rvol(1:npoint) = rvol(1:npoint) / rvmax
  endif

  if( icheck > 0 ) then
    write(3,150) npr, npoint, npso
    write(3,160) natome
  endif
  if( icheck > 1 .or. istop == 1 ) then
    if( base_hexa ) then
      write(3,170)
    else
      write(3,175)
    endif
    do i = 1,npoint
      write(3,180) i, indice(i,1:3), numia(i), xyz(1:4,i)*bohr, rvol(i)
    end do
    do i = npoint+1,npso
      write(3,180) i, indice(i,1:3), numia(i), xyz(1:4,i)*bohr
    end do
    if( istop == 1 ) stop
  endif

  return
  100 format(/' ---- Reseau -------',100('-'))
  130 format(//' Pour ix, iy, iz =',3i3,', point equivalent en ix1, iy1, iz1 =',3i3,/' non trouve dans reseau !'/ )
  150 format(/' npr =',i6,', npoint =',i6,', npso =',i6)
  160 format(' natome =',i4)
  170 format(/ ' Point index along hexagonal mesh',// '    i  ix  iy  iz  ia   xval    yval    zval    rval    rvol')
  175 format(/'    i  ix  iy  iz  ia    xval    yval    zval    rval    rvol')
  180 format(i5,4i4,5f8.4)
end

!***********************************************************************

subroutine posrel(Base_ortho,dcosxyz,iopsymr,x,a,u,dr,isym)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(nopsm):: iopsymr

  logical:: Base_ortho

  real(kind=db), dimension(3):: a, b, dcosxyz, x, u, v, w
  real(kind=db), dimension(3,3):: matopsym

  u(:) = x(:) - a(:)
  dr = vnorme(Base_ortho,dcosxyz,u)
  isym = 1

  do is = 2,nopsm
    select case(is)
      case(2,4,6,8,32,34,36,38,49,51,53,55)
        isi = is + 1
      case(3,5,7,9,33,35,37,39,50,52,54,56)
        isi = is - 1
      case(16,17,18,26,27,28)
        isi = is + 3
      case(19,20,21,29,30,31)
        isi = is - 3
      case default
        isi = is
    end select
    if( iopsymr(is) == 0 ) cycle
    call opsym(is,matopsym)
    b = matmul( matopsym, a )
    w(1:3) = x(1:3) - b(1:3)
    call opsym(isi,matopsym)
    v = matmul( matopsym, w )
    dr1 = vnorme(Base_ortho,dcosxyz,v)
    if( dr1 < dr-epspos ) then
      u(:) = v(:)
      dr = dr1
      isym = is * iopsymr(is)   ! pour avoir le signe
    endif
  end do

  return
end

!***********************************************************************

! Routine calculant les coefficients du laplacien.

subroutine laplac(Adimp,Base_hexa,cgrad,clapl,icheck, igrpt_nomag,indice,iopsymr,iord,ivois,isvois, &
                  mpirank,mpres,npso,npsom,nvois,nx)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(3):: la
  integer, dimension(0:nvois,3):: indpro
  integer, dimension(nopsm):: iopsymr
  integer, dimension(npsom,3):: indice
  integer, dimension(npsom,nvois):: ivois, isvois
  integer, dimension(-nx:nx,-nx:nx,-nx:nx):: mpres

  logical Base_hexa

  real(kind=db), dimension(2):: coeford
  real(kind=db), dimension(3):: v
  real(kind=db), dimension(nvois):: cgrad
  real(kind=db), dimension(0:nvois):: clapl

  if( icheck > 0 ) write(3,110)

! Coefficients du laplacien :

  indpro = 0

  ad2 = adimp**2
  coeford(1) = 1._db
  if( iord == 2 ) then
    cl = 1._db / ad2
  else
    cl = ( 4._db / 3._db ) / ad2
    coeford(2) = - 1._db / 16
  endif

  d3 = 2._db / 3
  iv = 0
  do io = 1,iord/2

    do i = 1,3
      do is = -1,1,2
        iv = iv + 1
        clapl(iv) = coeford(io) * cl
        if( base_hexa .and. ( i == 1 .or. i == 2 ) ) clapl(iv) = clapl(iv) * d3
        indpro(iv,i) = is * io
      end do
    end do

    if( base_hexa ) then
      do is = -1,1,2
        iv = iv + 1
        clapl(iv) = coeford(io) * cl * d3
        indpro(iv,1:2) = is * io
      end do
    endif

  end do

  clapl(0) = - sum( clapl(1:nvois) )

! Coefficients du gradient

  coeford(1) = 1._db
  if( iord == 2 ) then
    cd = 1 / ( 2 * adimp )
  else
    cd = (4._db/3._db) / ( 2 * adimp )
    coeford(2) = - 1._db/8
  endif
  f = 1 /sqrt(3._db)

  iv = 0
  do io = 1,iord/2
    if( base_hexa ) then
      do i = 1,4
        do is = -1,1,2
          iv = iv + 1
          cgrad(iv) = is * coeford(io) * cd
! i = 4 correspond a la direction (110)
          if( i == 2 .or. i == 4 ) cgrad(iv) = cgrad(iv) * f
        end do
      end do
    else
      do i = 1,3
        do is = -1,1,2
          iv = iv + 1
          cgrad(iv) = is * coeford(io) * cd
        end do
      end do
    endif
  end do

  if( icheck > 0 ) then
    write(3,130) nvois
    write(3,140)
    do iv = 0,nvois
      write(3,150) iv, indpro(iv,:), clapl(iv) * ad2
    end do
    write(3,160)
    do iv = 1,nvois
      write(3,150) iv, indpro(iv,:), cgrad(iv) * adimp
    end do
  endif

  f = sqrt(3._db) / 2
  f3 = 1 / sqrt(3._db)
  df3 = 2 * f3
  do i = 1,npso
    do  iv = 1,nvois

      la(1:3) = indice(i,1:3) + indpro(iv,1:3)
      if( base_hexa ) then
        v(1) = la(1) - 0.5_db * la(2)
        v(2) = f * la(2)
        v(3) = 1._db * la(3)
      else
        v(1:3) = 1._db * la(1:3)
      endif

! Si le voisin tombe en dehors de la maille, on le ramene dedans par
! la symetrie eventuelle.
      call posequiv(mpirank,v,iopsymr,isym,igrpt_nomag)
      if( base_hexa ) then
        la(1) = nint( v(1) + f3 * v(2) )
        la(2) = nint( df3 * v(2) )
        la(3) = nint( v(3) )
      else
        la(1:3) = nint( v(1:3) )
      endif

      if( abs(la(1)) <= nx .and. abs(la(2)) <= nx .and. abs(la(3)) <= nx ) then
        ivois(i,iv) = mpres(la(1),la(2),la(3))
        isvois(i,iv) = isym
      else
        ivois(i,iv) = 0
        isvois(i,iv) = 0
      endif

    end do
  end do

  if( icheck > 1 ) then
    write(3,170)
    do i = 1,npso
      write(3,180) i, ( ivois(i,iv), isvois(i,iv), iv = 1,nvois )
    end do
  endif

  return
  110 format(/' ---- Laplac -------',100('-'))
  130 format(/' nvois =',i3)
  140 format(/' iv  indpro  clapl * ad2')
  150 format(4i3,f12.6)
  160 format(/' iv  indpro  cgrad * adimp')
  170 format(/'    i  ivois, isvois')
  180 format(i6,16(i6,i4))
end

!***********************************************************************

! Routine effectuant la selection des points en bordure des spheres muffin-tin, et calculant les distances relatives de ces points au
! centre de l'atome correspondant.

subroutine bordure(Base_ortho,dcosxyz,Green,icheck,iopsymr,iord,iscratch,ivois,mpirank,natome,nbm,nbtm,nim, &
             npoint,npso,npsom,nsm,nstm,numia,nvois,posi,rvol,xyz)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(natome):: ibp, nbord, nbordf
  integer, dimension(nopsm):: iopsymr
  integer, dimension(npsom):: numia
  integer, dimension(npsom,nvois):: ivois
  integer, dimension(:), allocatable:: isrt
  integer, dimension(:,:), allocatable:: ibord, isbord

  logical:: Base_ortho, Green

  real(kind=db), dimension(3):: dcosxyz, p, v, x
  real(kind=db), dimension(natome):: pdharm
  real(kind=db), dimension(nim):: rvol
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(:), allocatable :: poidso
  real(kind=db), dimension(:,:), allocatable :: poidsa

! De ib = 1 a nbordf, on stocke les points frontiere exterieurs,
! de ib = nbordf+1 a nbord, les points interieurs.
  nbordf(1:natome) = 0
  nsortf = 0
  ibp(1:natome) = 0
  iss = 0
! Calcul des dimensions

  nv = 2 * nvois / iord

  do i = 1,npoint
    if( numia(i) /= 0 ) cycle
    do iv = 1,nv
      j = ivois(i,iv)
      if( j == 0 ) cycle
      ia = numia(j)
      if( ia == -2 ) then
        if( iss /= i ) nsortf = nsortf + 1
        iss = i
      elseif( ia > 0 .and. .not. green ) then
        if( ibp(ia) /= i ) nbordf(ia) = nbordf(ia) + 1
        ibp(ia) = i
      endif
    end do
  end do

  nbord(1:natome) = nbordf(1:natome)
  if( .not. green ) then
    do i = 1,npoint
      ia = numia(i)
      if( ia < 1 ) cycle
      do iv = 1,nvois
        j = ivois(i,iv)
        if( j == 0 .and. mpirank == 0 ) then
          call write_error
          do ipr = 3,9,3
            write(ipr,110) i, ia, iv
          end do
          stop
        endif
        if( numia(j) /= 0 ) cycle
        nbord(ia) = nbord(ia) + 1
        exit
      end do
    end do
  endif

  nsort = nsortf
  do i = npoint+1,npso
    do iv = 1,nvois
      j = ivois(i,iv)
      if( j == 0 ) cycle
      if( numia(j) /= 0 ) cycle
      nsort = nsort + 1
      exit
    end do
  end do

  do ia = 1,natome
    nbtm = max( nbtm, nbord(ia) )
    nbm = max( nbm, nbordf(ia) )
  end do
  nstm = max( nstm, nsort )
  nsm = max( nsm, nsortf )

  allocate( ibord(nbtm,natome) )
  allocate( isbord(nbtm,natome) )
  allocate( isrt(nsort) )
  allocate( poidsa(nbm,natome) )
  allocate( poidso(nsortf) )

  nbordf(1:natome) = 0
  nsortf = 0
  pdsort = 0._db
  pdharm(1:natome) = 0._db
  ibp(1:natome) = 0

! Les points en bordures sont seulement les points voisins dans les
! directions cristallographiques.
  do i = 1,npoint
    if( numia(i) /= 0 ) cycle
    do iv = 1,nv
      j = ivois(i,iv)
      if( j == 0 ) cycle
      ia = numia(j)
      if( ia == -2 ) then
        if( iss /= i ) then
          pdsort = pdsort + rvol(i)
          nsortf = nsortf + 1
          isrt(nsortf) = i
        endif
        iss = i
      elseif( ia > 0 .and. .not. green ) then
        if( ibp(ia) /= i ) then
          pdharm(ia) = pdharm(ia) + rvol(i)
          nbordf(ia) = nbordf(ia) + 1
          p(1:3) = posi(1:3,ia)
          x(1:3) = xyz(1:3,i)
          call posrel(Base_ortho,dcosxyz,iopsymr,x,p,v,dr,isym)
          ibord(nbordf(ia),ia) = i
          isbord(nbordf(ia),ia) = isym
        endif
        ibp(ia) = i
      endif
    end do
  end do

  nbord(1:natome) = nbordf(1:natome)
  if( .not. green ) then
    do i = 1,npoint
      ia = numia(i)
      if( ia < 1 ) cycle
      do iv = 1,nvois
        j = ivois(i,iv)
        if( j == 0 .and. mpirank == 0 ) then
          call write_error
          do ipr = 3,9,3
            write(ipr,110) i, ia, iv
          end do
          stop
        endif
        if( numia(j) /= 0 ) cycle

        nbord(ia) = nbord(ia) + 1
        p(1:3) = posi(1:3,ia)
        x(1:3) = xyz(1:3,i)
        call posrel(Base_ortho,dcosxyz,iopsymr,x,p,v,dr,isym)
        ibord(nbord(ia),ia) = i
        isbord(nbord(ia),ia) = isym
        exit
      end do
    end do
  endif

  nsort = nsortf
  do i = npoint+1,npso
    do iv = 1,nvois
      j = ivois(i,iv)
      if( j == 0 ) cycle
      if( numia(j) /= 0 ) cycle
      nsort = nsort + 1
      isrt(nsort) = i
      exit
    end do
  end do

  do ia = 1,natome
    if( nbordf(ia) > 0 ) fac = quatre_pi / pdharm(ia)
    do ib = 1,nbordf(ia)
      poidsa(ib,ia) = fac * rvol( ibord(ib,ia) )
    end do
  end do
  if( nsortf > 0 ) then
    fac = quatre_pi / pdsort
    do ib = 1,nsortf
      poidso(ib) = fac * rvol( isrt(ib) )
    end do
  endif

  if( mpirank == 0 ) then
    open(iscratch, status='SCRATCH')
    do ia = 1,natome
      write(iscratch,*) nbordf(ia), nbord(ia)
      if( nbord(ia) > 0 ) then
        write(iscratch,*) ibord(1:nbord(ia),ia)
        write(iscratch,*) isbord(1:nbord(ia),ia)
      endif
      if( nbordf(ia) > 0 ) write(iscratch,*) poidsa(1:nbordf(ia),ia)
    end do
    write(iscratch,*) nsortf, nsort
    if( nsort > 0 ) write(iscratch,*) isrt(1:nsort)
    if( nsortf > 0 ) write(iscratch,*) poidso(1:nsortf)
  endif

  if( icheck > 0 ) then
    write(3,120)
    if( .not. green ) then
      do ia = 1,natome
        write(3,130) ia, nbordf(ia), nbord(ia)
        if( nbordf(ia) > 0 .and. icheck > 1 ) then
          write(3,140)
          do ib = 1,nbordf(ia)
            write(3,150) ibord(ib,ia), isbord(ib,ia), poidsa(ib,ia)
          end do
        endif
        if( nbord(ia) > nbordf(ia) .and. icheck > 1) then
          do ib = nbordf(ia)+1,nbord(ia)
            write(3,150) ibord(ib,ia), isbord(ib,ia)
          end do
        endif
      end do
    endif
    write(3,160) nsortf, nsort
    if( nsortf > 0 .and. icheck > 1 ) then
      write(3,165)
      write(3,170) (isrt(ib), poidso(ib), ib = 1,nsortf)
      if( nsort > nsortf ) write(3,180) (isrt(ib), ib = nsortf+1,nsort)
    endif
  endif

  deallocate( ibord )
  deallocate( isbord )
  deallocate( isrt )
  deallocate( poidsa )
  deallocate( poidso )

  return
  110 format(//' Erreur dans bordure pour i, ia, iv =',3i6/)
  120 format(/' ---- Bordure ------',100('-')/)
  130 format('     Atome',i3,', nbordf =',i5,', nbord =',i5)
  140 format('  ibord  isbord   poidsa')
  150 format(2i6,f12.5)
  160 format(' Outer sphere, nsortf =',i5,', nsort =',i5)
  165 format('   isrt  poidso')
  170 format(10(i6,f10.5))
  180 format(10(i6,10x))
end

!***********************************************************************

subroutine recup_bordure(ibord,isbord,iscratch,isrt,mpinodes,mpirank,natome,nbord,nbordf,nbm,nbtm,nsm, &
        nsort,nsortf,nstm,poidsa,poidso)

  use declarations

  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer, dimension(nstm):: isrt
  integer, dimension(natome):: nbord, nbordf
  integer, dimension(nbtm,natome):: ibord, isbord

  real(kind=db), dimension(nsm):: poidso
  real(kind=db), dimension(nbm,natome):: poidsa

  if( mpirank == 0 ) then
    rewind(iscratch)
    do ia = 1,natome
      read(iscratch,*) nbordf(ia), nbord(ia)
      if( nbord(ia) > 0 ) then
        read(iscratch,*) ibord(1:nbord(ia),ia)
        read(iscratch,*) isbord(1:nbord(ia),ia)
      endif
      if( nbordf(ia) > 0 ) read(iscratch,*) poidsa(1:nbordf(ia),ia)
    end do
    read(iscratch,*) nsortf, nsort
    if( nsort > 0 ) read(iscratch,*) isrt(1:nsort)
    if( nsortf > 0 ) read(iscratch,*) poidso(1:nsortf)
    Close(iscratch)
  endif
  if( mpinodes > 1 ) then
    ndim = nbtm * natome
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nbordf,natome,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(nbord,natome,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ibord,ndim,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(isbord,ndim,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(poidsa,nbm*natome,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(nsortf,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nsort,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(isrt,nstm,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(poidso,nsm,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  endif

  return
end

!***********************************************************************

! Calcul de l'energie du niveau de coeur initial.

subroutine Energseuil(Core_resolved,Delta_Epsii,Delta_Eseuil,Epsii,Epsii_moy,Eseuil,icheck,is_g, &
           itabs_nonexc,lseuil,m_g,mpirank,nbseuil,ninit1,ninitl,ninitlr,nr,nrm,nseuil,nspine,ntype,numat,psii,rato,Rmtg, &
           Rmtsd,V_abs_i,V_intmax,V0bde)

  use declarations
  implicit none

  integer icheck, initl, ipr, isol, isp, itabs_nonexc, jsp, lmax_pot_loc, lseuil, m, mpirank, nbseuil, ninit1, ninitl, &
    ninitlr, nlm_pot_loc, nm1, nm2, nr, nrm, nseuil, nsol, nspin, nspine, nspino, nspinp, ntype, numat

  parameter(nspin = 2, nspino = 2, nspinp = 2)

  logical:: Core_resolved, Relativiste, Spinorbite, Ylm_comp

  integer, dimension(ninitl):: is_g
  integer, dimension(ninitl,2):: m_g

  real(kind=db):: Delta, Delta_Epsii, Delta_Eseuil, Ep_moy, Epsii_moy, j, mj, psiHpsi, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nbseuil):: Epsii_m, Eseuil
  real(kind=db), dimension(ninitlr):: Epsii
  real(kind=db), dimension(ninitl):: dEpsii, Epsi
  real(kind=db), dimension(nspine):: V0bde
  real(kind=db), dimension(nspin):: V0bd
  real(kind=db), dimension(nr):: psi, r
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nrm,nspine):: V_abs_i
  real(kind=db), dimension(nr,1,nspin):: V_abs
  real(kind=db), dimension(0:nrm,0:ntype):: rato

  if( icheck > 0 ) write(3,110)

  Spinorbite = .true.
  Relativiste = .true.

! La contribution non diagonale du potentiel est negligee
  Ylm_comp = .false.
  lmax_pot_loc = 0
  nlm_pot_loc = 1
  do isp = 1,nspin
    jsp = min(isp,nspine)
    V0bd(isp) = V0bde(jsp)
    V_abs(1:nr,1,isp) = V_abs_i(1:nr,jsp)
  end do

  r(1:nr) = rato(1:nr,itabs_nonexc)

! Le calcul de l'energie des niveaux initiaux est toujours avec
! spinorbite
  nm1 = - lseuil - 1
  nm2 = lseuil

  do initl = 1,ninitl

    if( initl <= ninit1 ) then
      psi(1:nr) = psii(1:nr,1)
    else
      psi(1:nr) = psii(1:nr,2)
    endif
    if( m_g(initl,1) == 1000 ) then
      m = nm1
      nsol = 1
      isol = 2
    elseif( m_g(initl,2) == 1000 ) then
      m = nm2
      nsol = 1
      isol = 1
    else
      m = m_g(initl,1)
      nsol = 2
      if( is_g(initl) == 1 ) then
        isol = 1
!        isol = 2
      else
        isol = 2
!        isol = 1
      endif
    endif

    if( initl == 1 .and. icheck > 0 ) write(3,*)

    Epsi(initl) = - psiHpsi(.true.,icheck,isol,lseuil,lmax_pot_loc,m,nseuil,nlm_pot_loc,nr,nr,nsol,nspin, &
               nspino,nspinp,numat,V0bd,V_abs,psi,r,Relativiste,Rmtg,Rmtsd,Spinorbite,V_intmax,Ylm_comp)
  end do

  Epsii_m(1) = sum( Epsi(1:ninit1) ) / ninit1

  if( ninit1 < ninitl ) then
    Epsii_m(2) = sum( Epsi(ninit1+1:ninitl) ) / ( ninitl - ninit1 )
    Delta = ( Eseuil(1) - Eseuil(2) ) - ( Epsii_m(1) - Epsii_m(2) )
    dEpsii(1:ninit1) = Epsi(1:ninit1)
    dEpsii(ninit1+1:ninitl) = Epsi(ninit1+1:ninitl)
  else
    dEpsii(1:ninitl) = Epsi(1:ninitl)
  endif

  if( Core_resolved ) then
    Epsii(1:ninitl) = dEpsii(1:ninitl)
  else
    Epsii(1:ninitlr) = Epsii_m(1:ninitlr)
  endif

  if( Core_resolved ) then
    if( ninit1 /= ninitl ) then
      Epsii_moy = sum( Epsii(ninit1+1:ninitl) ) / ( ninitl - ninit1 )
      Ep_moy = sum( Epsii(1:ninit1) ) / ninit1
    else
      Epsii_moy = sum( Epsii(1:ninitl) ) / ninitl
      Ep_moy = Epsii_moy
    endif
  else
    Epsii_moy = Epsii(ninitlr)
    Ep_moy = Epsii(1)
  endif
  Delta_Eseuil = Ep_moy - Epsii_moy

  if( Delta_Epsii < 100000._db .and. nbseuil > 1 ) then
    if( Core_resolved ) then
      Delta = Delta_Epsii - Delta_Eseuil
      Epsii(1:ninit1) = Epsii(1:ninit1) + Delta
    else
      Epsii(1) = Epsii(ninitlr) + Delta_Epsii
    endif
    Delta_Eseuil = Delta_Epsii
  endif

  do ipr = 3,6,3

    if( lseuil == 0 .or. ninitl == 2*lseuil + 2 ) then
      j = lseuil + 0.5_db
    else
      j = lseuil - 0.5_db
    endif

    if( ( ipr == 3 .and. icheck > 0 ) .or. ( ipr == 6 .and. Core_resolved .and. mpirank == 0 ) ) then
      if( ipr == 3 ) write(ipr,115)
      if( ipr == 3 .and. ninit1 /= ninitl ) write(ipr,120) Delta * rydb
      mj = - j - 1
      do initl = 1,ninitl
        mj = mj + 1
        if( initl == ninit1 + 1 ) then
          j = j + 1
          mj = - lseuil - 0.5_db
        endif
        if( Delta_Epsii < 100000._db .and. nbseuil > 1 ) then
          write(ipr,130) dEpsii(initl)*rydb, Epsi(initl)*rydb, initl, j, mj
        else
          write(ipr,135) dEpsii(initl)*rydb, initl, j, mj
        endif
      end do
    endif

    if( .not. Core_resolved .and. ( ( ipr == 3 .and. icheck > 0) .or. ( ipr == 6 .and. mpirank == 0 ) ) ) then
      if( ipr == 3 ) then
        if( ninitlr == 1 ) then
          Write(ipr,140) Epsii(:) * Rydb
        else
          Write(ipr,150) Epsii(:) * Rydb
        endif
      else
        if( ninitlr == 1 ) then
          Write(ipr,160) Epsii(:) * Rydb
        else
          Write(ipr,170) Epsii(:) * Rydb
        endif
      endif
    endif

  end do

  return
  110 format(/' ---- Energseuil ---',100('-'))
  115 format(/' Kohn-Scham energy of the initial states:',/)
  120 format(4x,'Edge(1) - E_edge(2) - ( <E_KS(1)> - <E_KS(2)> ) =', f6.3,' eV',/)
  130 format(4x,'Epsii =',f10.3,' eV,  E_KS =',f10.3, ' eV,  State',i3,', j =',f5.1,', mj =',f5.1)
  135 format(4x,'Epsii =',f10.3, ' eV,  State',i3,', j =',f5.1,', mj =',f5.1)
  140 format(/7x,'Epsii used =',f10.3,' eV')
  150 format(/7x,'Epsii used =',f10.3,'  and',f10.3,' eV')
  160 format(' Epsii used =',f10.3,' eV')
  170 format(' Epsii used =',f10.3,'  and',f10.3,' eV')
end

!***********************************************************************

! Calcul de l'energie d'un niveau lie par E = <psi|H|psi> / <psi|psi>
! Dans le cas relativiste, il faut plusieurs cycles car l'Hamiltonien
! depend de l'energie calculee.

function psiHpsi(Cal_psi,icheck,isol,l,lmax_pot,m,n,nlm_pot,nr,nrm,nsol,nspin, &
              nspino,nspinp,numat,V0bd,Vrato,psi,r,Relativiste,Rmtg,Rmtsd,Spinorbite,V_intmax,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, ie_cycle, im, ip, ip_cycle, ir, isol, isp, l, l2, lmax, lmax_pot, m, m_hubb, m1, m2, n, nie_cycle, &
    nip_cycle, nlm, nlm_pot, nodes, nr, nr_max, nr_pos, nrm, nrmtg, nrmtsd, nsol, nspin, nspino, nspinp, numat

  parameter(m_hubb = 0)

  logical:: Cal_psi, Ecomp, Full_potential, Hubb_a, Hubb_d, Radial_comp, Relativiste, Renorm, Spinorbite, Ylm_comp

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(2*l+1,nspin,2*l+1,nspin):: Tau
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  real(kind=db):: E, Eimag, Ep, Epp, fac, f_integr, f_so, Hpsi, &
    k_onde, php, pp, psiHpsi, Rmtg, Rmtsd, spp, td, V_intmax, vme
  real(kind=db), dimension(nspin)::  Ecinetic, V0bd
  real(kind=db), dimension(nrm):: dr, f2, psi, psi2, r
  real(kind=db), dimension(nrm,nspin):: g0, gm, gp, gpi, ur
  real(kind=db), dimension(nrm,nlm_pot,nspin):: Vrato
  real(kind=db), dimension(nrm,nspino):: gso
  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: uis, urs

  if( icheck > 1 ) then
    if( Spinorbite ) then
      write(3,110) n, l, m
    else
      write(3,115) n, l
    endif
  endif

  Eimag = 0._db
  Ecomp = .false.
! Pour le calcul de l'energie des orbitales de coeur on neglige
! l'aspect non spherique
  Full_potential = .false.
  Hubb_a = .false.
  Hubb_d = .true.
  m2 = 1
  m1 = 2 * l + 1
  lmax = l
  Renorm = .false.
  Radial_comp = .false.

  if( Full_potential ) then
    nlm = ( lmax + 1 )**2
  else
    nlm = 1
  endif
  allocate( V(nr,nlm,nlm,nspin) )

  call mod_V(icheck,lmax,lmax_pot,nlm,nlm_pot,nrm,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  dr(1) = 0.5_db * ( r(2) + r(1) )
  do ir = 2,nrmtg-1
    dr(ir) = 0.5_db * ( r(ir+1) - r(ir-1) )
  end do
  dr(nrmtg) = Rmtg - 0.5 * ( r(nrmtg-2) + r(nrmtg-1) )

  nip_cycle = 200

  if( Relativiste .or. Spinorbite ) then
    nie_cycle = 10
  else
    nie_cycle = 1
  endif

  Epp = 0._db

  l2 = l**2 + l

  ur(:,:) = 0._db

  if( Spinorbite ) then
    if( m == -l-1 ) then
      isp = 2
      ur(:,2) = psi(:)
    elseif( m == l ) then
      isp = 1
      ur(:,1) = psi(:)
    else
      ur(:,1) = psi(:) / sqrt(2._db)
      ur(:,2) = psi(:) / sqrt(2._db)
    endif
    f_so = sqrt( (l - m) * ( l + m + 1._db) )
  else
    isp = 1
    ur(:,1) = psi(:)
  endif
  psi2(1:nrm) = psi(1:nrm)**2
  pp = f_integr(r,psi2,nrmtg,1,nrm,Rmtg)
  ur(:,:) = ur(:,:) / sqrt( pp )

  E = 0._db

! cycle pour calcul de la fonction d'onde
  do ip_cycle = 1,nip_cycle

    do ie_cycle = 1,nie_cycle  ! cycle pour le calcul de l'energie

      if( ip_cycle == 1 .and. ie_cycle == 1 ) then
        call coef_sch_rad(E,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,.false.,.false.,V)
        gpi(:,:) = 1 / gp(:,:)
      elseif( Relativiste .or. Spinorbite ) then
        call coef_sch_rad(E,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)
        gpi(:,:) = 1 / gp(:,:)
      endif

      php = 0._db

      do isp = 1,nspin
        if( Spinorbite ) then
          if( m == -l-1 .and. isp == 1 ) cycle
          if( m == l .and. isp == 2 ) cycle
        endif
        do ir = 2,nrmtg
! Il faut ajouter E car g0 contient V - E
          td = g0(ir,isp) + l2 * f2(ir) + E
          if( Spinorbite .and. (ip_cycle > 1 .or. ie_cycle >1)) then
            if( isp == 1 ) then
              td = td + m * gso(ir,isp)
            else
              td = td - ( m + 1 ) * gso(ir,isp)
            endif
          endif

          Hpsi = gm(ir,isp) * ur(ir-1,isp) +  td * ur(ir,isp) + gp(ir,isp) * ur(ir+1,isp)
          if( nsol > 1 .and. (ip_cycle > 1 .or. ie_cycle >1) ) Hpsi = Hpsi + f_so * gso(ir,isp) * ur(ir,3-isp)
          php = php + ur(ir,isp) * Hpsi * dr(ir)
        end do
      end do

      Ep = E
      E =  php

      if( icheck > 1 ) write(3,120) ip_cycle, ie_cycle, E*Rydb, pp

      if( ie_cycle == 1 ) cycle

      if( abs( Ep - E ) < 0.001_db ) exit

    end do

    if( .not. cal_psi ) exit

    if( abs( Epp - E ) < 0.001_db ) exit

    Epp = E

    Ecinetic(:) = E - V0bd(:)

    do

! Konde ne sert pas, car on ne calcule pas Tau
      konde(:) = sqrt( cmplx(E - V0bd(:), Eimag,db) )

      allocate( urs(nr,m1,m2,nspin,nspino) )
      allocate( uis(nr,m1,m2,nspin,nspino) )
      call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde, &
         l,lmax,m_hubb,nlm,m1,m2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau, &
         uis,urs,V,V_hubb)

      if( Spinorbite ) then
        if( m == -l-1 ) then
          ur(:,1) = 0._db
          ur(:,2) = urs(:,1,m2,2,nspino)
        elseif( m == l ) then
          ur(:,1) = urs(:,m1,m2,1,1)
          ur(:,2) = 0._db
        else
          ur(:,1) = urs(:,l+1+m,m2,1,isol)
          ur(:,2) = urs(:,l+1+m+1,m2,2,isol)
        endif
      else
        ur(:,:) = urs(:,min(l+1+m,m1),m2,:,1)
      endif
      deallocate( urs, uis )

      do ir = nrmtg,1,-1
        if( V(ir,1,1,1) < E .or. V(ir,1,1,nspin) < E ) exit
      end do
      nr_pos = ir
      nr_pos = min(nr_pos,nr-2)

      nodes = 0
      if( m == -l-1 ) then
        isp = 2
      else
        isp = 1
      endif
      do ir = 2,min(nr_pos,nrmtg)
        if( ur(ir,isp)*ur(ir-1,isp) < 0._db ) nodes = nodes + 1
      end do
      if( nodes == n - l - 1 .or. E > V(nrmtg,1,1,1) ) exit

      if( nodes > n - l - 1 ) then
        E = E - 0.1_db * abs(E)
      else
        E = E + 0.1_db * abs(E)
      endif
      if( icheck > 1 ) write(3,125) nodes, n-l-1, E*Rydb
    end do

    if( E < V(nrmtg,1,1,1) ) then

! Inward integration
      do ir = nr,1,-1
        if( ( V(ir,1,1,1) - E ) * r(ir)**2 < 625._db ) exit
      end do
      nr_max = ir

      do ir = nr_max+1,nr
        ur(ir,:) = 0._db
      end do

      do isp = 1,nspin
        if( Spinorbite ) then
          if( m == -l-1 .and. isp == 1 ) cycle
          if( m == l .and. isp == 2 ) cycle
        endif

        do ir = nr_max,nr_max-1,-1
          VmE = V(ir,1,1,isp) - E
          VmE = Max( VmE, 0.000001_db )
          k_onde = sqrt( VmE )
          ur(ir,isp) = exp( - k_onde * r(ir) ) / r(ir)**(l+1)
        end do

        fac = ur(nr_pos+1,isp)
        do ir = nr_max-1,nr_pos+2,-1
          im = ir - 1
          ip = ir + 1
          td = g0(ir,isp) + l2 * f2(ir)
          if( Spinorbite ) then
            if( isp == 1 ) then
              td = td + m * gso(ir,isp)
            else
              td = td - ( m + 1 ) * gso(ir,isp)
            endif
          endif
          ur(im,isp) = - ( td * ur(ir,isp) + gp(ir,isp) * ur(ip,isp) ) / gm(ir,isp)
        end do
        fac = fac / ur(nr_pos+1,isp)
        do ir = nr_pos+1,nr_max
          ur(ir,isp) = fac * ur(ir,isp)
        end do

      end do

    else

      nr_max = nrmtg

    endif


    psi2(:) = 0._db
    do isp = 1,nspin
      psi2(1:nrm) = psi2(1:nrm) + ur(1:nrm,isp)**2
    end do

    pp = f_integr(r,psi2,nr,1,nrm,r(nr_max))
    spp = sqrt( pp )

    do isp = 1,nspin
      if( Spinorbite ) then
        if( m == -l-1 .and. isp == 1 ) cycle
        if( m == l .and. isp == 2 ) cycle
      endif
      ur(:,isp) = ur(:,isp) / spp
    end do

  end do

  if( icheck > 2 ) then
    if( nspin == 2 ) then
      write(3,130)
    else
      write(3,135)
    endif
    do ir = 1,nr
      write(3,140) r(ir)*bohr, psi(ir), ur(ir,:), V(ir,1,1,:)*Rydb
    end do
  endif

  psiHpsi =  E  ! les valeurs sont negatives

  deallocate( V )

  return
  110 format(/' Core wave function : n =',i2,', l =',i2,', m =',i2/)
  115 format(/' Core wave function : n =',i2,', l =',i2/)
  120 format(' Cycle',2i4,',  E = ',f10.3,' eV,  <psi.psi> =',f10.7)
  125 format('   Nodes =',i2,' >',i2,' -->  E = ',f10.3,' eV')
  130 format('    rato        psi      psi_new(up)   psi_new(dn)       Vr_up       Vr_dn')
  135 format('    rato        psi        psi_new        V')
  140 format(f10.5,1p,5e13.5)
end

!***********************************************************************

! Calcul des (l,m) de chaque representation

subroutine lmrep(Green,iaprotoi,iato,icheck,iopsym_atom,iopsymr,irep_util,iso,itypei,karact,lato, &
             lmaxat,lmaxso,lso,mato,mpirank,mso,n_atom_proto,natome,nbordf,ngrph,nlmsa0,nlmsam,nlmso0, &
             nso1,nsortf,nspino,ntype,numat,Orthmati,posi,rot_atom)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db), dimension(nopsm):: kopsymo
  complex(kind=db), dimension(nopsm,nrepm):: karact

  integer, dimension(0:ntype) :: numat
  integer, dimension(natome):: iaprotoi, itypei
  integer, dimension(ngrph):: nlmso0
  integer, dimension(0:n_atom_proto):: lmaxat
  integer, dimension(nopsm):: iopsyma, iopsymr, irotiops
  integer, dimension(nopsm,natome):: iopsym_atom
  integer, dimension(nrepm,2):: irep_util
  integer, dimension(natome,ngrph):: nlmsa0
  integer, dimension(natome):: nbordf
  integer, dimension(nso1,ngrph):: iso, lso, mso
  integer, dimension(nlmsam,natome,ngrph):: iato, lato, mato

  logical green

  real(kind=db), dimension(3):: p
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(3,3):: Orthmati, rot_a
  real(kind=db), dimension(3,3,natome):: rot_atom

  if( icheck > 1 ) write(3,110)

  nlmso0(:) = 0
  nlmsa0(:,:) = 0

  do ia = 0,natome

    if( ia == 0 ) then

      if( green ) cycle
      lmax = lmaxso
      iopsyma(:) = iopsymr(:)
      do is = 1,nopsm
        irotiops(is) = is
      end do

    else

      ipr = iaprotoi(ia)
      lmax = lmaxat(ipr)
      if( icheck > 1 .or. (ia == 1 .and. icheck == 2) ) then
        do k = 1,3
          p(k) = sum( orthmati(k,1:3) * posi(1:3,ia) )
        end do
        write(3,130) ia, p(1:3)
      endif
      p(1:3) = posi(1:3,ia)
      iopsyma(:) = iopsym_atom(:,ia)
      rot_a(:,:) = rot_atom(:,:,ia)
      call iop_rot(icheck,irotiops,rot_a)
    endif

    do l = 0,lmax
      do m = -l,l

        do ispin = 1,nspino

          call symorb(l,m,kopsymo)

! Recherche de la representation a laquelle appartient l'orbitale
          boucle_irep: do igrph = 1,ngrph
!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if( ia == 2 .and. (igrph == 3 .or. igrph == 5 ) .and. l == 4 .and. abs(m) == 4 ) cycle
            irep = abs( irep_util(igrph,ispin) )
            if( irep == 0 ) cycle
            do is = 1,nopsm
              js = irotiops(is)
              if( js == 0 ) cycle
              if( iopsymr(is) == 0 .or. iopsyma(js) == 0 ) cycle
              if( irep_util(igrph,ispin) > 0 ) then
                if( abs( karact(is,irep) - kopsymo(js) ) > eps10 ) cycle boucle_irep
              else
                if( abs( conjg(karact(is,irep)) - kopsymo(js) ) > eps10 ) cycle boucle_irep
              endif
            end do

            if( ia == 0 ) then
              nlmso0(igrph) = nlmso0(igrph) + 1
              iso( nlmso0(igrph), igrph ) = ispin
              lso( nlmso0(igrph), igrph ) = l
              mso( nlmso0(igrph), igrph ) = m
            else
              nlmsa0(ia,igrph) = nlmsa0(ia,igrph) + 1
              iato( nlmsa0(ia,igrph), ia, igrph ) = ispin
              lato( nlmsa0(ia,igrph), ia, igrph ) = l
              mato( nlmsa0(ia,igrph), ia, igrph ) = m
            endif

          end do boucle_irep

        end do
      end do
    end do

    if( icheck > 1 ) then

      if( ia == 0 ) then

        do igrph = 1,ngrph
          write(3,170) lmax, nlmso0(igrph), igrph
          write(3,180)
          write(3,190) (lm, lso(lm,igrph), mso(lm,igrph), iso(lm,igrph), lm = 1,nlmso0(igrph))
        end do

      else

        do igrph = 1,ngrph
          write(3,175) ia, lmax, nlmsa0(ia,igrph), igrph
          write(3,180)
          write(3,190) (lm, lato(lm,ia,igrph), mato(lm,ia,igrph), iato(lm,ia,igrph), lm = 1,nlmsa0(ia,igrph))
        end do

      endif

    endif

  end do

  if( .not. green ) then

    istop = 0
    do igrph = 1,ngrph

      nsp = nspino * nsortf
      if( nlmso0(igrph) > nsp .and. mpirank == 0 ) then
        if( istop == 0 ) call write_error
        do ipr = 3,9,3
          write(ipr,200) nlmso0(igrph), nsp
        end do
        istop = 1
      endif

      do ia = 1,natome
        if( numat( itypei(ia) ) == 0 ) cycle
        nsp = nspino * nbordf(ia)
        if( nlmsa0(ia,igrph) > nsp .and. mpirank == 0 ) then
          if( istop == 0 ) call write_error
          do ipr = 3,9,3
            write(ipr,210) ia, nlmsa0(ia,igrph), nsp
          end do
          istop = 1
        endif
      end do

    end do

    if( istop == 1 ) stop

  endif

  return
  110 format(/' ---- lmrep --------',100('-'))
  130 format(/' Atome ia =',i3,', position =',3f7.3)
  170 format(/' lmax =',i3,', nlm =',i3,', igrph =',i2)
  175 format(/'  ia = ',i3,', lmax =',i3,', nlm =',i3,', igrph =',i2)
  180 format(/'   lm    l    m    is')
  190 format(4i5)
  200 format(///'   nlmso0 =',i4,' > nspino*nsortf =',i4, // &
    ' Solution :',/'  Reduce the maximum energy or ',/ &
    '  Reduce the outersphere lmax using the lmaxso keyword or',/ &
    '  Reduce the interpoint distance using the keyword adimp'//)
  210 format(///'   nlmsa0(ia=',i2,') =',i4,' > nspino*nbordf =',i4, // &
    ' Solution :',/'  Reduce the maximum energy or ',/ &
    '  Reduce the atomic lmax using the lmax keyword or',/ '  Reduce the interpoint distance using the keyword adimp'//)
end

!***********************************************************************

! Reduction de la symetrie compte tenu de la position de l'atome
! Rotation de la base locale tel que Oz inchange et Ox direction radiale.

subroutine iop_rot(icheck,irotiops,rot_atom)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(nopsm):: irotiops

  real(kind=db), dimension(3,3):: matopsym, rot_atom, rot_tem

  irotiops(:) = 0

  do is = 1,nopsm
    call opsym(is,matopsym)
    rot_tem = matmul( transpose(rot_atom), matmul( matopsym, rot_atom ) )
    do js = 1,nopsm
      call opsym(js,matopsym)
      if( sum( abs( matopsym(:,:) - rot_tem(:,:) ) ) > eps6 ) cycle
      irotiops(js) = is
      exit
    end do
  end do

  if( icheck > 1 ) write(3,110) irotiops(:)

  return
  110 format(/' irotiops =',/ 4(1x,5i3) )
end

!***********************************************************************

! Calcul du nombre d'harmoniques pour l'energie courante.

subroutine clmax(energ,r,lmax0,lmax,Z,lmaxfree)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer Z

  logical lmaxfree

  if( energ > 0._db .and. r > eps10 ) then
    rk = sqrt( energ ) * r
  else
    rk = 0._db
  endif
  rk2 = rk**2
  if( lmax0 < 0 ) then
    lmax = int( 0.5_db + 0.5_db * sqrt(1 + 4*rk2) ) - lmax0
  else
    lmax = int( 0.5_db + 0.5_db * sqrt(1 + 4*rk2) ) + 2
  endif
  if( rk < 1._db ) then
    lmax = min(6,lmax)
  elseif( rk < 1.7_db ) then
    lmax = min( nint(4.30*rk+0.30),lmax )
  else
    lmax = min( nint(3.04*rk+3.83),lmax )
  endif
  if( Z == 1 ) then
    lmax = 2
  elseif( Z > 54 ) then
    lmax = max( 3,lmax )
  else
    lmax = max( 2,lmax )
  endif

  if( .not. lmaxfree ) then
    if( Z < 3 ) then
      lmax = min(lmax,3)
    elseif( Z < 19 ) then
      lmax = min(lmax,4)
    elseif( Z < 37 ) then
      lmax = min(lmax,5)
    elseif( Z < 55 ) then
      lmax = min(lmax,6)
    elseif( Z < 87 ) then
      lmax = min(lmax,7)
    else
      lmax = min(lmax,8)
    endif
  endif

  if( lmax0 > 0 ) lmax = min( lmax, lmax0 )

  return
end

!***********************************************************************

! Calcul des Ylm et des distances au centre de l'atome de chaque point.

subroutine ylmpt(Base_ortho,dcosxyz,iaprotoi,ibord,icheck, iopsymr,isrt,lmaxat,lmaxso,n_atom_proto,natome, &
             nbord,nbtm,nlmmax,nlmomax,nsort,nstm,npsom,posi, rot_atom,rot_atom_abs,xyz,ylmato,ylmso)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db), dimension(:), allocatable :: ylmc

  integer, dimension(0:n_atom_proto):: lmaxat
  integer, dimension(nopsm):: iopsymr
  integer, dimension(nstm):: isrt
  integer, dimension(natome):: iaprotoi, nbord
  integer, dimension(nbtm,natome):: ibord

  logical:: Base_ortho

  real(kind=db), dimension(3):: dcosxyz, p, v, w
  real(kind=db), dimension(3,3):: rot_a, rot_atom_abs
  real(kind=db), dimension(3,3,natome):: rot_atom
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(nbtm,nlmmax,natome):: ylmato
  real(kind=db), dimension(nsort,nlmomax):: ylmso
  real(kind=db), dimension(:), allocatable :: ylmr

  if( icheck > 1 ) write(3,110)

  do ia = 0,natome

    if( ia > 0 ) then
      ipr = iaprotoi(ia)
      lmax = lmaxat(ipr)
      np = nbord(ia)
      p(1:3) = posi(1:3,ia)
! Changement de base :
      rot_a(:,:) = rot_atom(:,:,ia)
    else
      p(:) = 0._db
      lmax = lmaxso
      np = nsort
      rot_a(:,:) = rot_atom_abs(:,:)
    endif
    nlmr = ( lmax + 1 )**2
    nlmc = ( ( lmax + 1 ) * ( lmax + 2 ) ) / 2
    allocate( ylmc(nlmc) )
    allocate( ylmr(nlmr) )

    do ib = 1,np
      if( ia == 0 ) then
        i = isrt(ib)
      else
        i = ibord(ib,ia)
      endif
      w(1:3) = xyz(1:3,i)
      call posrel(Base_ortho,dcosxyz,iopsymr,w,p,v,r,isym)
      w = matmul( rot_a, v )
      call cylm(lmax,w,r,ylmc,nlmc)
! Conversion en Ylm reels
      call ylmcr(lmax,nlmc,nlmr,ylmc,ylmr)
      if( ia == 0 ) then
        ylmso(ib,1:nlmr) = ylmr(1:nlmr)
      else
        ylmato(ib,1:nlmr,ia) = ylmr(1:nlmr)
      endif
    end do

    if( icheck > 2 ) then
      write(3,120) ia
      lm = 0
      do l = 0,lmax
        do m = -l,l
          write(3,130) l, m
          lm = lm + 1
          if( ia == 0 ) then
            write(3,140) ylmso(1:np,lm)
          else
            write(3,140) ylmato(1:np,lm,ia)
          endif
        end do
      end do
    endif

    deallocate( ylmc )
    deallocate( ylmr )

  end do

  return
  110 format(/' ---- Ylmpt --------',100('-'))
  120 format(/'    Ylm,  Atom',i3)
  130 format(/' l, m =',2i3)
  140 format(1p,8f10.3)
end

!***********************************************************************

! Calcul des Ylm complexes

subroutine cylm(lmax,v,r,ylmc,nlm)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db):: ylmc(nlm), exphi, sint_exphi

  real(kind=db), dimension(3):: v

  if( r > eps6 ) then
    cost = v(3) / r
  else
    cost = 1._db
  endif

  den = sqrt( v(1)**2 + v(2)**2 )
  if( den > eps6 ) then
    exphi = cmplx( v(1), v(2), db ) / den
    sint = sqrt( 1 - cost**2 )
    sint_exphi = cmplx( v(1), v(2), db ) / r
    cott = cost / sint
  else
    exphi = (1._db,0._db)
    sint = 0._db
    sint_exphi = (0._db,0._db)
  endif

! Calcul de Y(0,0) :
  ylmc(1) = 1 / sqrt( quatre_pi )

! Calcul des Y(l,l) :
  do l = 1,lmax
    lm = ((l+1)*(l+2)) / 2
    lm1 = (l*(l+1)) / 2
    f = - sqrt(1 + 0.5_db/l)
    ylmc(lm) = f * sint_exphi * ylmc(lm1)
  end do

! Calcul de Y(1,0) :
  if( lmax > 0 ) ylmc(2) = sqrt( 3 / quatre_pi ) * cost

! Calcul des Y(l,l-1) :
  do l = 2,lmax
    lm = ( ( l + 1 ) * ( l + 2 ) ) / 2 - 1
    lm1 = ( l**2 + l ) / 2 - 1
    f = - sqrt( (2*l + 1._db) / (2*l - 2) )
    ylmc(lm) = f * sint_exphi * ylmc(lm1)
  end do

! Calcul des Y(l,m) :
  exphi = conjg( exphi )

  do l = 2,lmax
    lm0 = ( l**2 + l ) / 2
    do m = l-2,1,-1
      lm = lm0 + m + 1
      f = - 2 * (m + 1._db) / sqrt( l*(l+1._db) - m*(m+1._db) )
      g = - sqrt( ( l*(l+1._db) - (m+1._db)*(m+2._db) ) / ( l*(l+1._db) - m*(m+1._db) ) )
      if( abs(sint) > eps10 ) then
        ylmc(lm) = exphi*( cott*f*ylmc(lm+1) + exphi*g*ylmc(lm+2) )
      else
        ylmc(lm) = (0._db,0._db)
      endif
    end do
  end do

! Calcul des Y(l,0) :
  do l = 2,lmax
    lm = ( l * (l+1) ) / 2 + 1
    lm1 = ( (l-1)*l ) / 2 + 1
    lm2 = ( (l-2)*(l-1) ) / 2 + 1
    f = sqrt( (2*l + 1._db) / (2*l - 1._db) ) * (2*l - 1._db) / l
    g = - sqrt( (2*l + 1._db) / (2*l - 3._db) ) * (l - 1._db) / l
    ylmc(lm) = f * cost * ylmc(lm1) + g * ylmc(lm2)
  end do

  return
end

!***********************************************************************

! Conversion des Ylm complexes en Ylm reels.

subroutine ylmcr(lmax,nlmc,nlmr,ylmc,ylmr)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db), dimension(nlmc):: ylmc
  real(kind=db), dimension(nlmr):: ylmr

  rac_2 = sqrt(2._db)

  lm = 0
  do l = 0,lmax

    l2 = ( l**2 + l ) / 2 + 1

    do m = -l,l

      lm = lm + 1
      lm0 = l2 + abs(m)

      if( m < 0 ) then
        ylmr(lm) = rac_2 * aimag( ylmc(lm0) )
      elseif( m == 0 ) then
        ylmr(lm) = real( ylmc(lm0), db )
      else
        ylmr(lm) = rac_2 * real( ylmc(lm0), db )
      endif

    end do
  end do

  return
end

!***********************************************************************

! Selection des harmoniques spheriques compte tenu des symetries

subroutine ylmsym(iopsym,lmax,lv,mpirank,mv,nlms,nlmtot)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer iops(nopsm), iopsym(nopsm), lv(nlms), mv(nlms)
  logical imparite, parite, parrot, parsym

  iops(:) = abs( iopsym(:) )
  parsym = .false.

  lm = 0

  do l = 0,lmax
    do m = -l,l

! Centrosymetrie
      if( iops(25) == 1 ) then
        parsym = iopsym(25) == 1
        parite = mod(l+2*abs(m),2) == 0
        if( (parsym .and. .not. parite) .or. (.not. parsym .and. parite)  ) cycle
      endif

! Plan Ox
      if( iops(40) == 1 ) then
        parsym = iopsym(40) == 1
        m2 = mod(abs(m),2)
        parite = (m >= 0 .and. m2 == 0) .or. (m < 0 .and. m2 == 1)
        if( (parsym .and. .not. parite) .or. (.not. parsym .and. parite)  ) cycle
      endif

! Plan Oy
      if( iops(41) == 1 ) then
        parsym = iopsym(41) == 1
        parite = m >= 0
        if( (parsym .and. .not. parite) .or. (.not. parsym .and. parite)  ) cycle
      endif

! Plan Oz
      if( iops(42) == 1 ) then
        parsym = iopsym(42) == 1
        parite = mod(l+abs(m),2) == 0
        if( (parsym .and. .not. parite) .or. (.not. parsym .and. parite)  ) cycle
      endif

! Plan diagonal 0z
      if( iops(44) == 1 .or. iops(45) == 1 ) then
        parsym = ( iopsym(44) == 1 .or. iopsym(45) == 1 )
        m4 = mod(abs(m),4)
        parite = (m4 == 0 .and. m >= 0) .or. (m4 == 2 .and. m < 0)
        imparite = (m4 == 2 .and. m >= 0) .or. (m4 == 0 .and. m < 0)
        if( ( (parsym .and. .not. parite) .or. (.not. parsym .and. .not. imparite) ) &
             .and. .not. ( (m4 == 1 .or. m4 == 3) .and. (m > 0) .and. iops(18) /= 1 ) ) cycle
      endif

! Axes de rotation selon Oz
      if( iops(51) == 1 ) then
        irot = 6
        parsym = iopsym(51) == 1
      elseif( iops(55) == 1 ) then
        irot = 6
        parsym = iopsym(55) == 1
      elseif( iops(18) == 1 ) then
        irot = 4
        parsym = iopsym(18) == 1
      elseif( iops(28) == 1 ) then
        irot = 4
        parsym = iopsym(28) == 1
      elseif( iops(49) == 1 ) then
        irot = 3
        parsym = iopsym(49) == 1
      elseif( iops(53) == 1 ) then
        irot = 3
        parsym = iopsym(53) == 1
      elseif( iops(24) == 1 ) then
        irot = 2
        parsym = iopsym(24) == 1
      else
        irot = 1
      endif

      if( irot > 1) then
        mr = mod(abs(m),irot)
        ml = mod(l+abs(m),2)
        parrot = mod(irot,2) == 0
        if( abs(iopsym(28)) == 1 .and. iopsym(18) == 0 ) then
          parite = ( mr == 0 .and. ml == 0 ) .or. ( mr == 2 .and. ml == 1 )
          imparite = ( mr == 2 .and. ml == 0 ).or. ( mr == 0 .and. ml == 1 )
        else
          parite = mr == 0
          imparite = parrot .and. mr == irot/2
        endif
        if( (parsym .and. .not. parite) .or. (.not. parsym .and. .not. imparite)  ) cycle
      endif

! Axe 2 selon 0x
      if( iops(22) == 1 ) then
        parsym = iopsym(22) > 0
        parite = (m >= 0 .and. mod(l+m,2) == 0) .or. (m < 0 .and. mod(l-m,2) == 1)
        if( (parsym .and. .not. parite) .or. (.not. parsym .and. parite)  ) cycle
      endif

! Axe 2 selon 0y
      if( iops(23) == 1 ) then
        parsym = iopsym(23) > 0
        parite = (m >= 0 .and. mod(l+2*m,2) == 0) .or. (m < 0 .and. mod(l-2*m,2) == 1)
        if( (parsym .and. .not. parite) .or. (.not. parsym .and. parite)  ) cycle
      endif

! Axes 4 selon 0x ou 0y
      if( mpirank == 0 .and. ( iops(16) == 1 .or. iops(17) == 1 .or. iops(26) == 1 .or. iops(27) == 1 ) ) then
        call write_error
        do ipr = 3,9,3
          write(ipr,110)
        end do
        stop
      endif

      lm = lm + 1
      if( lm > nlms .and. mpirank == 0 ) then
        call write_error
        do ipr = 3,9,3
          write(ipr,120) lm, nlms
        end do
        stop
      endif
      lv(lm) = l
      mv(lm) = m

    end do
  end do

  nlmtot = lm

  return
  110 format(/' Les axes 4 selon Ox et Oy ne sont pas programmes !')
  120 format(/' lm =',i4,' > nlms =',i4)
end

!*********************************************************************

! Calcul de la racine de la densite d'etats avec energie complexe.
! Prise comme la convolution par une lorenzienne de la densite d'etats du vide.

function cal_norm(Energ,Eimag)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db):: k

  Eimag2 = Eimag**2
  de_int = abs( Eimag ) / 20
  Delta_E = 500 * abs( Eimag )
  E1 = Energ - Delta_E
  E2 = Energ + Delta_E
  je1 = nint( E1 / de_int )
  je1 = max( 1, je1 )
  je2 = nint( E2 / de_int )
  je2 = max( je2, 10 )

  rint = 0._db
  rho = 0._db

  do je = je1,je2
    e = ( je - 0.5 ) * de_int
    k = sqrt( e )
    fac = de_int / ( ( Energ - e )**2 + Eimag2 )
    rho = rho + k * fac
    rint = rint + fac
  end do
  rho = rho * abs(Eimag) / pi**2

  rint = rint * abs(Eimag) / pi

  cal_norm = sqrt( rho )

end

!***********************************************************************

! Calcul des fonctions de bessel et neuman.

subroutine cbess(fnorm,z,lmax,lmaxm,bessel)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db):: bessel(0:lmaxm), fnorm, z

  if( abs( z ) < eps10 ) then
    bessel(0) = fnorm
  else
    bessel(0) = fnorm * sin(z) / z
  endif

  if( lmax == 0 ) return

  if( abs( z ) < eps10 ) then

    bessel(1:lmax) = (0._db, 0._db)

  else

    bessel(1) = bessel(0) / z - fnorm * cos(z) / z

    do l = 2,lmax
      l1 = 2*l - 1
      bessel(l) = l1 * bessel(l-1) / z - bessel(l-2)
    end do

  endif

  return
end

!***********************************************************************

! Calcul des fonctions de bessel et neuman.

subroutine cbessneu(fnorm,z,lmax,lmaxm,bessel,neuman)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db):: bessel(0:lmaxm), fnorm, neuman(0:lmaxm), z

  neuman(0) = - fnorm * cos(z) / z
  bessel(0) = fnorm * sin(z) / z

  if( lmax == 0 ) return

  neuman(1) = neuman(0) / z - bessel(0)
  bessel(1) = bessel(0) / z + neuman(0)

  do l = 2,lmax
    l1 = 2*l - 1
    neuman(l) = l1 * neuman(l-1) / z - neuman(l-2)
    bessel(l) = l1 * bessel(l-1) / z - bessel(l-2)
  end do

  return
end

!***********************************************************************

! Calcul des fonctions de bessel et neuman.

subroutine cbessneur(fnorm,z,lmax,lmaxm,bessel,neuman)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db):: bessel(0:lmaxm), neuman(0:lmaxm), z

  neuman(0) = - fnorm * cos(z) / z
  bessel(0) = fnorm * sin(z) / z

  neuman(1) = neuman(0) / z - bessel(0)
  bessel(1) = bessel(0) / z + neuman(0)

  do l = 2,lmax
    l1 = 2*l - 1
    neuman(l) = l1 * neuman(l-1) / z - neuman(l-2)
    bessel(l) = l1 * bessel(l-1) / z - bessel(l-2)
  end do

  return
end

!***********************************************************************

! Calcul de la fonction de Bessel

subroutine cbessel(bess,ip0,lmax,nr,q,r)

  use declarations
  implicit none

  integer:: ip0, l, l1, lmax, nr

  real(kind=db):: q
  real(kind=db), dimension(nr):: r, z
  real(kind=db), dimension(nr,0:lmax):: bessel
  real(kind=db), dimension(nr,ip0:lmax):: bess

  z(:) = q * r(:)
  bessel(:,0) = sin( z(:) ) / z(:)
  if( lmax > 0 ) bessel(:,1) = ( bessel(:,0) - cos( z(:) ) ) / z(:)

  do l = 2,lmax
    l1 = 2*l - 1
    bessel(:,l) = l1 * bessel(:,l-1) / z(:) - bessel(:,l-2)
  end do

  do l = ip0,lmax
     bess(:,l) = bessel(:,l)
  end do

  return
end

!***********************************************************************

subroutine cnlmmax(lmax,lv,nlmam,nlm)

  integer lv(nlmam)

  do lm = nlmam,1,-1
    if( lv(lm) <= lmax ) exit
  end do
  nlm = lm

  return
end

!***********************************************************************

subroutine irep_comp(Atom_comp0,green,iopsymr,irep_util, karact,ngrph,nspino,Repres_comp)

  use declarations
  implicit none

  integer, intent(in):: ngrph, nspino
  logical, intent(in):: Atom_comp0, green

  integer, dimension(nopsm), intent(in):: iopsymr
  integer, dimension(nrepm,2), intent(in):: irep_util

  complex(kind=db), dimension(nopsm,nrepm), intent(in):: karact
  logical, dimension(ngrph), intent(out):: Repres_comp

  integer igrph, irep, is, isp

  Repres_comp(:) = .false.

  do igrph = 1, ngrph
    do isp = 1,nspino
      irep = abs( irep_util(igrph,isp) )
      if( green ) then
        if( irep == 0 ) cycle
      else
! Il n'y a pas de couplage up-down dans le cas suivant
        if( irep == 0 ) irep = abs( irep_util(igrph,3-isp) )
        if( .not. Atom_comp0 ) cycle
      end if
      do is = 1,nopsm
        if( green .and. iopsymr(is) == 0 ) cycle
        if( abs( aimag(karact(is,irep)) ) > eps10 ) Repres_comp(igrph) = .true.
      end do
    end do
  end do

  return
end

!***********************************************************************

! Calcul de l'integrale de 0 a Radius de f.
! L'integrale est calculee avec un polynome d'interpolation d'ordre deux
! r: contient les valeurs des rayons pour les points qu'on integre
! nr: nombre de points qu'on integre

  real(kind=db) function f_integr(r,fct,nr,ir0,nrm,Radius)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db), dimension(ir0:nrm):: fct, r

  tiers = 1._db / 3

  f_integr = 0._db

  do i = 2,nr-1
    rm = r(i-1)
    r0 = r(i)
    rp = r(i+1)
    if( i == 2 ) then
      x1 = 0._db
    else
      x1 = 0.5_db * ( rm + r0 )
    endif
    if( rp > Radius ) then
      x2 = Radius
    else
      if( i == nr - 1 ) then
        x2 = rp
      else
        x2 = 0.5_db * ( rp + r0 )
      endif
      fm = fct(i-1)
      f0 = fct(i)
      fp = fct(i+1)
      if( abs(fp) < 1e-20_db .and. abs(f0) < 1e-20_db .and. abs(fm) < 1e-20_db ) cycle

      dp0 = rp - r0
      dpm = rp - rm
      d0m = r0 - rm
      a = ( fm * dp0 - f0 * dpm + fp * d0m ) / ( d0m * dp0 * dpm )
      b = ( f0 - fm ) / d0m - a * ( r0 + rm )
      c = f0 - a * r0**2 - b * r0
    endif

    f_integr = f_integr + ( tiers * a * ( x1**2 + x1 * x2  + x2**2 ) + 0.5_db * b * ( x1 + x2 ) + c ) * ( x2 - x1 )

    if( rp > Radius ) exit

  end do

  return
end

!***********************************************************************

! Calcul de l'integrale de 0 a Radius de f.
! L'integrale est calculee avec un polynome d'interpolation d'ordre 3
! r: contient les valeurs des rayons pour les points qu'on integre
! On considere que seule les valeurs jusqu'à l'indice n tel que
! r(i) > Radius > Radius(i-1) sont sures.

  real(kind=db) function f_integr3(r,fct,ir0,nrm,Radius)

  use declarations
  implicit none

  integer:: i, ir0, nrm

  logical:: This_is_the_end

  real(kind=db):: a, b, c, d, d0m, dp0, dpm, f0, f1, f2, f3, f4, fac1, fac2, fm, fp, r0, Radius, r1, r2, r3, r4, rap14, rap24, &
                  rap34, rm, rp, Tiers
  real(kind=db), dimension(ir0:nrm):: fct, r

  Tiers = 1._db / 3
  This_is_the_end = .false.
  f_integr3 = 0._db

  do i = 1,nrm-1
    if( i == 1 ) then
      if( ir0 == 1 ) then
 ! Construction de l'ordonnée par extrapolation du second degre
        rm = r(1)
        r0 = r(2)
        rp = r(3)
        fm = fct(1)
        f0 = fct(2)
        fp = fct(3)
        if( abs(fp) < 1e-20_db .and. abs(f0) < 1e-20_db .and. abs(fm) < 1e-20_db ) cycle
        dp0 = rp - r0
        dpm = rp - rm
        d0m = r0 - rm
        a = ( fm * dp0 - f0 * dpm + fp * d0m ) / ( d0m * dp0 * dpm )
        b = ( f0 - fm ) / d0m - a * ( r0 + rm )
        c = f0 - a * r0**2 - b * r0
        f1 = c
      else
        f1 = fct(i-1)
      endif
      r1 = 0._db
      r2 = r(i)
      f2 = fct(i)
      r3 = r(i+1)
      f3 = fct(i+1)
      r4 = r(i+2)
      f4 = fct(i+2)
      rm = 0._db
      rp = r3
    elseif( i == nrm-1 .or. r(i+1) > Radius ) then
      rm = r(i)
      rp = r(nrm)
    else
      r1 = r(i-1)
      f1 = fct(i-1)
      r2 = r(i)
      f2 = fct(i)
      r3 = r(i+1)
      f3 = fct(i+1)
      r4 = r(i+2)
      f4 = fct(i+2)
      rm = r2
      rp = r3
    endif
    if( rp > Radius ) then
      rp = Radius
      This_is_the_end = .true.
    endif
    if( abs(f2) < 1e-20_db .and. abs(f3) < 1e-20_db ) cycle

    rap14 = (f1 - f4) / (r1 - r4)
    rap24 = (f2 - f4) / (r2 - r4)
    rap34 = (f3 - f4) / (r3 - r4)
    fac1 = ( rap14 - rap34 ) / (r1 - r3)
    fac2 = ( rap24 - rap34 ) / (r2 - r3)
    a = ( fac1 - fac2 ) / ( r1 - r2 )
    b = fac1 - a * ( r1 + r3 + r4 )
    c = rap14 - a * ( r1**2 + r1 * r4 + r4**2 ) - b * ( r1 + r4 )
    d = f1 - ( a * r1**2 + b * r1 + c ) * r1

    f_integr3 = f_integr3 + ( 0.25_db * ( rm + rp ) * ( rm**2 + rp**2 ) * a + tiers * ( rm**2 + rm*rp + rp**2 ) * b &
              + 0.5_db * ( rm + rp ) * c + d ) * ( rp - rm )

    if( This_is_the_end ) exit

  end do

  return
end

!***********************************************************************

! Calcul de l'integrale de 0 a Radius de f. f(0) = 0.
! L'integrale est calculee avec un polynome d'interpolation d'ordre 3
! r: contient les valeurs des rayons pour les points qu'on integre
! nr: nombre de points qu'on integre

  real(kind=db) function f_integr32(r,fct,ir0,nr,R0,Radius)

  use declarations
  implicit none

  integer:: i, ir0, ir1, nr

  logical This_is_the_end

  real(kind=db):: a, b, c, d, f1, f2, f3, f4, fac1, fac2, R0, r1, r2, r3, r4, rap14, rap24, rap34, rm, rp, Radius, tiers
  real(kind=db), dimension(ir0:nr):: fct, r

  tiers = 1._db / 3
  This_is_the_end = .false.
  f_integr32 = 0._db
  ir1 = ir0 + 1

  do i = ir1,nr-2
    if( r(i) < R0 ) cycle
    r1 = r(i-1)
    f1 = fct(i-1)
    r2 = r(i)
    f2 = fct(i)
    r3 = r(i+1)
    f3 = fct(i+1)
    r4 = r(i+2)
    f4 = fct(i+2)
    if( i == ir1 .or. r1 < R0 ) then
      rm = R0
    else
      rm = r2
    endif
    if( i == nr-2 .or. r4 > Radius ) then
      rp = Radius
      This_is_the_end = .true.
    else
      rp = r3
    endif

    if( abs(f2) < 1e-20_db .and. abs(f3) < 1e-20_db ) cycle

    rap14 = (f1 - f4) / (r1 - r4)
    rap24 = (f2 - f4) / (r2 - r4)
    rap34 = (f3 - f4) / (r3 - r4)
    fac1 = ( rap14 - rap34 ) / (r1 - r3)
    fac2 = ( rap24 - rap34 ) / (r2 - r3)
    a = ( fac1 - fac2 ) / ( r1 - r2 )
    b = fac1 - a * ( r1 + r3 + r4 )
    c = rap14 - a * ( r1**2 + r1 * r4 + r4**2 ) - b * ( r1 + r4 )
    d = f1 - ( a * r1**2 + b * r1 + c ) * r1

    f_integr32 = f_integr32 + ( 0.25_db * ( rm + rp ) * ( rm**2 + rp**2 ) * a + tiers * ( rm**2 + rm*rp + rp**2 ) * b &
               + 0.5_db * ( rm + rp ) * c + d ) * ( rp - rm )

    if( This_is_the_end ) exit

  end do

  return
end

!***********************************************************************

  real(kind=db) function f_interp1(x,x1,x2,y1,y2)

  use declarations
  implicit none

  real(kind=db),intent(in):: x, x1, x2, y1, y2
  real(kind=db):: p

  p = (x - x1) / (x2 - x1)
  f_interp1 = p*y2 + (1-p)*y1

  return
end

!***********************************************************************

  real(kind=db) function f_interp2(r,rm,r0,rp,fm,f0,fp)

  use declarations
  implicit none

  real(kind=db):: a, b, c, d0m, dp0, dpm, f0, fm, fp, r0, rm, rp, r

  dp0 = rp - r0
  dpm = rp - rm
  d0m = r0 - rm

  a = ( fm * dp0 - f0 * dpm + fp * d0m ) / ( d0m * dp0 * dpm )
  b = ( f0 - fm ) / d0m - a * ( r0 + rm )
  c = f0 - a * r0**2 - b * r0

  f_interp2 = a * r**2 + b * r + c

  return
end

!**********************************************************************

! Interpolation d'ordre 3, pour obtenir fonction et derivee

subroutine interp3(f,fp,r,r1,r2,r3,r4,f1,f2,f3,f4)

  use declarations
  implicit none

  real(kind=db):: a, b, c, d, f, f1, f2, f3, f4, fac1, fac2, fp, r, r1, r2, r3, r4, rap14, rap24, rap34

  rap14 = (f1 - f4) / (r1 - r4)
  rap24 = (f2 - f4) / (r2 - r4)
  rap34 = (f3 - f4) / (r3 - r4)
  fac1 = ( rap14 - rap34 ) / (r1 - r3)
  fac2 = ( rap24 - rap34 ) / (r2 - r3)
  a = ( fac1 - fac2 ) / ( r1 - r2 )
  b = fac1 - a * ( r1 + r3 + r4 )
  c = rap14 - a * ( r1**2 + r1 * r4 + r4**2 ) - b * ( r1 + r4 )
  d = f1 - ( a * r1**2 + b * r1 + c ) * r1

  f = a * r**3 + b * r**2 + c * r + d
  fp = 3 * a * r**2 + 2 * b * r + c

  return
end

!**********************************************************************

! Interpolation d'ordre 3

  real(kind=db) function f_interp3(r,r1,r2,r3,r4,f1,f2,f3,f4)

  use declarations
  implicit none

  real(kind=db):: a, b, c, d, f1, f2, f3, f4, fac1, fac2, r, r1, r2, r3, r4, rap14, rap24, rap34

  rap14 = (f1 - f4) / (r1 - r4)
  rap24 = (f2 - f4) / (r2 - r4)
  rap34 = (f3 - f4) / (r3 - r4)
  fac1 = ( rap14 - rap34 ) / (r1 - r3)
  fac2 = ( rap24 - rap34 ) / (r2 - r3)
  a = ( fac1 - fac2 ) / ( r1 - r2 )
  b = fac1 - a * ( r1 + r3 + r4 )
  c = rap14 - a * ( r1**2 + r1 * r4 + r4**2 ) - b * ( r1 + r4 )
  d = f1 - ( a * r1**2 + b * r1 + c ) * r1

  f_interp3 = a * r**3 + b * r**2 + c * r + d

  return
end

!**********************************************************************

! Calcul du volume de la maille elementaire

function Cal_Volume_maille(axyz,angxyz)

  use declarations
  implicit none

  real(kind=db):: Cal_Volume_maille, cosa, cosb, cosc
  real(kind=db), dimension(3):: angxyz, axyz

  cosa = cos( pi * angxyz(1) / 180 )
  cosb = cos( pi * angxyz(2) / 180 )
  cosc = cos( pi * angxyz(3) / 180 )
  Cal_Volume_maille = axyz(1)*axyz(2)*axyz(3) * sqrt( 1 - cosa**2 - cosb**2 - cosc**2 + 2*cosa*cosb*cosc )

  return
end



