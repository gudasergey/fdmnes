! Fdmnes subroutine
! Calculation of the atomic electronic densities using dirac-slater or
! Hartree-Fock Slater.

subroutine dirgen(icheck,it,itabs,jseuil,lcoeur,lqnexc,lseuil,lvval,mpirank,n_orbexc,nbseuil,ncoeur,nlat,nlatm, &
          nnlm,nonexc,nqnexc,n_ray,nrato_dirac, nrm,nseuil,nspin,ntype,nvval,pop_open_val,popatc,popatv, &
          popexc,popval,psi_coeur,psii,psi_open_val,psival,rato,rho_coeur,rhoit,Z,Relativiste)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: i, iaug, ibav, icheck, io, ip, ipr, ir, irel, it, itabs, j, jo, jseuil, ko, l, l_level_val, lseuil, mpirank, &
            mseuil, n_coeur, n_orb, n_orbexc, nbseuil, nlatm, nmax, nnlm, n_ray, nrato_dirac, nrm, nseuil, nspin, ntype, Z
  
  integer, dimension(nnlm):: lqn, lqnexc, nqn, nqnexc
  integer, dimension(2,0:ntype):: lcoeur, ncoeur
  integer, dimension(0:ntype):: nlat
  integer, dimension(0:ntype,nlatm):: lvval, nvval

  logical:: nonexc, relativiste

  real(kind=db):: Charge, dp, E_total, E_total_exc, elmax, h_ray, p, p1, p2, Popt, ppp, Ptot, Ray_max
  
  real(kind=db), dimension(nnlm):: nel, pop, rqn
  real(kind=db), dimension(nrm):: rato
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(0:ntype):: popatc
  real(kind=db), dimension(0:ntype,nlatm):: popatv
  real(kind=db), dimension(0:ntype,nlatm,nspin):: popval
  real(kind=db), dimension(nnlm,nspin):: popexc
  real(kind=db), dimension(0:nrm,0:ntype):: rhoit, rho_coeur
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(nrm,2):: psi_open_val
  real(kind=db), dimension(2):: pop_open_val
  real(kind=db), dimension(0:nrm,2,0:ntype):: psi_coeur

  real(kind=db), dimension(:), allocatable:: ray, rho
  real(kind=db), dimension(:,:), allocatable:: psi, psi_small

  if( icheck > 1 ) write(3,110) it, Z
  
  if( Z <= 0 .or. Z > 103 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,115)
    end do
    stop
  endif

  n_ray = nrato_dirac ! Number of radius in the radial mesh
  ray_max = 20._db   ! Maximum radius (au)
  h_ray = n_ray * 54._db / 600      ! log step

  if( n_ray > nrm .and. mpirank == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,120) n_ray, nrm
    end do
    stop
  endif
  ibav = 3

  if( it == itabs ) then
    select case(jseuil)
      case(1,2,4,6)
        mseuil = 1
      case default
        mseuil = 2
    end select
  endif

  if( Z > 2 .or. relativiste ) then
    irel = 1
  else
    irel = 0
    rqn(:) = 0._db
  endif
  if( irel == 0 ) mseuil = 1

  call config(Z,irel,n_coeur,n_orb,nnlm,nqn,lqn,rqn,nel)

  pop(1:n_orb) = nel(1:n_orb)

! On remplace les population par defaut par celles qui sont donnees en
! entree
  boucle_io: do io = 1,nlat(it)
    Popt = sum( popval(it,io,1:nspin) )
    do ip = 1,n_orb
      if(nqn(ip) /= nvval(it,io) .or. lqn(ip) /= lvval(it,io)) cycle
      if( irel == 0 .or. lqn(ip) == 0 ) then
        pop(ip) = popt
      else
        pop(ip) = ( lqn(ip) / ( 2*lqn(ip) + 1._db ) ) * popt
        pop(ip+1) = popt - pop(ip)
      endif
      cycle boucle_io
    end do
! Si en entree on trouve une orbitale qui normalement n'existe pas:
    n_orb = n_orb + 1
    if( n_orb > nnlm ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,130)
      end do
      stop
    endif
    nqn(n_orb) = nvval(it,io)
    lqn(n_orb) = lvval(it,io)
    if( lvval(it,io) == 0 ) then
      rqn(n_orb) = 0.5_db
    else
      rqn(n_orb) = lvval(it,io) - 0.5_db
    endif
    if( irel == 0 ) then
      pop(n_orb) = popt
    else
      pop(n_orb) = ( lqn(n_orb) / ( 2*lqn(n_orb)+1._db ) ) * popt
      if( lqn(n_orb) /= 0 ) then
! Si relativiste la nouvelle orbitale trouvee a un splitting
        n_orb = n_orb + 1
        if( n_orb > nnlm ) then
          call write_error
          do ipr = 3,9,3
            write(ipr,130)
          end do
          stop
        endif
        nqn(n_orb) = nvval(it,io)
        lqn(n_orb) = lvval(it,io)
        rqn(n_orb) = lvval(it,io) + 0.5_db
        pop(n_orb) = popt - pop(n_orb-1)
      endif
    endif
  end do boucle_io

! construction de l'atome neutre (le calcul des niveaux atomiques se
! fait avec des atomes neutres)

  Charge = Z - sum( pop(1:n_orb) )
  if( charge > eps6 ) then
    dp = charge
    do io = 1,n_orb
      if( ( nqn(io) <= nseuil .and. it == itabs ) .or. lqn(io) > 1 ) cycle
      if( irel == 0 .or. lqn(io) == 0 ) then
        elmax = 2._db + 4._db * lqn(io)
        if( pop(io) > elmax - eps6 ) cycle

        do j = 1,nlat(it)
          if(nqn(io) /= nvval(it,j) .or. lqn(io)/=lvval(it,j)) cycle
          exit
        end do
        if( j == nlat(it)+1 ) then
          nlat(it) = nlat(it) + 1
          if( nlat(it) > nlatm ) then
            call write_error
            do ipr = 3,9,3
              write(ipr,140) nlat(it), nlatm
            end do
            stop
          endif
          nvval(it,nlat(it)) = nqn(io)
          lvval(it,nlat(it)) = lqn(io)
          popval(it,nlat(it),1:nspin) = pop(io) / nspin
        endif

        pop(io) = pop(io) + dp
        dp = pop(io) - elmax
        if( dp < eps6 ) exit
        pop(io) = elmax
      else
        if( lqn(io-1) == lqn(io) ) cycle
        elmax = 2._db + 4._db * lqn(io)
        if( sum(pop(io:io+1)) > elmax - eps6 ) cycle

        do j = 1,nlat(it)
          if(nqn(io) /= nvval(it,j) .or. lqn(io)/=lvval(it,j)) cycle
          exit
        end do
        if( j == nlat(it)+1 ) then
          nlat(it) = nlat(it) + 1
          if( nlat(it) > nlatm ) then
            call write_error
            do ipr = 3,9,3
              write(ipr,140) nlat(it), nlatm
            end do
            stop
          endif
          nvval(it,nlat(it)) = nqn(io)
          lvval(it,nlat(it)) = lqn(io)
          popval(it,nlat(it),1:nspin) = sum(pop(io:io+1)) / nspin
        endif

        p = lqn(io) / ( 2 * lqn(io) + 1._db )
        pop(io) = pop(io) + p * dp
        pop(io+1) = pop(io+1) + ( 1 - p ) * dp
        dp = sum( pop(io:io+1) ) - elmax
        if( dp < eps6 ) exit
        pop(io) = 2._db * lqn(io)
        pop(io+1) = 2 + 2._db * lqn(io)
      endif
    end do
    if( dp > eps6 ) then
      nmax = 0
      do io = 1,n_orb
        nmax = max(nqn(io),nmax)
      end do

      if( nmax == 1 ) then
        iaug = 1
      else
        iaug = 0
        do io = 1,n_orb
          if( nqn(io) == nmax .and. lqn(io) == 1 ) then
            iaug = 1
            exit
          endif
        end do
      endif
      n_orb = n_orb + 1
      if( n_orb > nnlm ) then
        call write_error
        do ipr = 3,9,3
          write(ipr,130)
        end do
        stop
      endif
      if( iaug == 1 ) then
        nqn(n_orb) = nmax + 1
        lqn(n_orb) = 0
        pop(n_orb) = dp
        if( irel == 1 ) rqn(n_orb) = 0.5_db
      else
        nqn(n_orb) = nmax
        lqn(n_orb) = 1
        if( irel == 0 ) then
          pop(n_orb) = dp
        else
          rqn(n_orb) = 0.5_db
          pop(n_orb) = dp / 3
          n_orb = n_orb + 1
          if( n_orb > nnlm ) then
            call write_error
            do ipr = 3,9,3
              write(ipr,130)
            end do
            stop
          endif
          nqn(n_orb) = nmax
          lqn(n_orb) = 1
          rqn(n_orb) = 1.5_db
          pop(n_orb) = 2 * dp / 3
        endif
      endif
      nlat(it) = nlat(it) + 1
      if( nlat(it) > nlatm ) then
        call write_error
        do ipr = 3,9,3
          write(ipr,140) nlat(it), nlatm
        end do
        stop
      endif
      nvval(it,nlat(it)) = nqn(n_orb)
      lvval(it,nlat(it)) = lqn(n_orb)
      popval(it,nlat(it),1:nspin) = 0._db
    endif

  elseif( charge < - eps6 ) then

    do io = n_orb,1,-1
      dp = charge
      if( ( nqn(io) <= nseuil .and. it == itabs ) .or. lqn(io) > 1 ) cycle

      do j = 1,nlat(it)
        if(nqn(io) /= nvval(it,j) .or. lqn(io)/=lvval(it,j)) cycle
        exit
      end do
      if( j == nlat(it)+1 ) then
        nlat(it) = nlat(it) + 1
        if( nlat(it) > nlatm ) then
          call write_error
          do ipr = 3,9,3
            write(ipr,140) nlat(it), nlatm
          end do
          stop
        endif
        nvval(it,nlat(it)) = nqn(io)
        lvval(it,nlat(it)) = lqn(io)
        if( irel == 0 .or. lqn(io) == 0 ) then
          popval(it,nlat(it),1:nspin) = pop(io) / nspin
        else
          popval(it,nlat(it),1:nspin) = sum(pop(io-1:io)) / nspin
        endif
      endif

      if( irel == 0 .or. lqn(io) == 0 ) then
        pop(io) = pop(io) + dp
      else
        pop(io) = pop(io) + 2 * dp / 3
        pop(io-1) = pop(io-1) + dp / 3
      endif
      exit
    end do

  endif

  allocate( ray(n_ray) )
  allocate( rho(n_ray) )
  allocate( psi(n_ray,nnlm) )
  allocate( psi_small(n_ray,nnlm) )

! Remplissage de popatv: ces populations correspondent a l'atome neutre

  if( it /= itabs .or. nonexc ) then
    do io = 1,nlat(it)
      do j = 1,n_orb
        if(nqn(j) /= nvval(it,io) .or. lqn(j) /= lvval(it,io)) cycle
        if( irel == 0 .or. lqn(j) == 0 ) then
          popatv(it,io) = pop(j)
        else
          popatv(it,io) = pop(j) + pop(j+1)
        endif
        exit
      end do
    end do
  endif

  if( icheck > 0 .and. it /= 0 ) then
    if( nlat(it) /= 0 ) then
      write(3,'(/A/A/A)') ' The occupancies below are the one used for the self-consistent calculation', &
                      ' of the atomic radial wave functions of the occupied levels', &
       ' The real occupancy corresponding to the demand in the indata file can be different and is applied afterwards'
    endif
      
    if( irel == 0 ) then
      write(3,150) it, Z
      do io = 1,n_orb
        if( Z > 18 .and. nqn(io) < 3 ) cycle
        if( Z > 36 .and. nqn(io) < 4 ) cycle
        if( Z > 86 .and. nqn(io) < 5 ) cycle
        write(3,160) nqn(io), lqn(io), pop(io)
      end do
    else
      write(3,170) it, Z
      do io = 1,n_orb
        if( Z > 18 .and. nqn(io) < 3 ) cycle
        if( Z > 36 .and. nqn(io) < 4 ) cycle
        if( Z > 86 .and. nqn(io) < 5 ) cycle
        write(3,180) nqn(io), lqn(io), rqn(io), pop(io)
      end do
    endif
  endif

! Influence par ce qu'on lui indique en entree, mais
! neanmoins pop etaient modifiees de ce que l'atome soit neutre

  call dirac(E_total,h_ray,icheck,ibav,irel,lqn,n_orb,n_ray,nnlm,nqn,pop,psi,psi_small,ray,ray_max,rho,rqn,Z)

! Construction de l'atome excite (appel en init_run):

  if( it == itabs ) then

    l = l_level_val(Z)

! Recuperation de la fonction d'onde de l'orbitale de valence potentiellement occupee
    do io = n_orb,1,-1
      if( lqn(io) /= l ) cycle
      if( irel == 1 .and. lqn(io) /= 0 ) then
        psi_open_val(1:n_ray,1) = 0.5_db * ( psi(1:n_ray,io) + psi(1:n_ray,io-1) )
        pop_open_val(1) = pop(io) + pop(io-1)
      else
        psi_open_val(1:n_ray,1) = psi(1:n_ray,io)
        pop_open_val(1) = pop(io)
      end if
      exit
    end do

    if( nseuil == 0 ) then ! Cas de l'optic
      psii(:,:) = 0._db
    else
! Recuperation de la fonction d'onde de coeur dite initiale
      do j = 1,n_orb
        if( nqn(j) /= nseuil .or. lqn(j) /= lseuil ) cycle
        if( irel == 0 .or. lqn(j) == 0 ) then
          psii(1:n_ray,1) = psi(1:n_ray,j)
          if( nbseuil == 2 ) psii(1:n_ray,nbseuil) = psi(1:n_ray,j)
        elseif( mseuil == 1 ) then
          psii(1:n_ray,1:nbseuil) = psi(1:n_ray,j:j+nbseuil-1)
        else
          psii(1:n_ray,1) = psi(1:n_ray,j+1)
        endif
        exit
      end do
    endif

    jo = 0
    ko = 0

    do io = 1,n_orbexc

      if( jo < nnlm ) jo = jo + 1
      ko = ko + 1

      nqn(jo) = nqnexc(io)
      lqn(jo) = lqnexc(io)
      Ptot = sum( popexc(io,1:nspin) )

      if( irel == 0 .or. lqn(jo) == 0 ) then
        pop(jo) = ptot
        if( irel == 1 ) rqn(jo) = 0.5_db
      else
        if( nqn(jo) == nseuil .and. lqn(jo) == lseuil ) ptot = ptot + 1._db
        p = lqn(jo) / ( 2 * lqn(jo) + 1._db )
        pop(jo) = p * ptot
        if( nqn(jo) == nseuil .and. lqn(jo) == lseuil .and. mseuil == 1 ) pop(jo) = pop(jo) - 1._db
        nqn(jo) = nqnexc(io)
        lqn(jo) = lqnexc(io)
        rqn(jo) = lqn(jo) - 0.5_db
        jo = jo + 1
        nqn(jo) = nqnexc(io)
        lqn(jo) = lqnexc(io)
        rqn(jo) = lqn(jo) + 0.5_db
        pop(jo) = ( 1 - p ) * ptot
        if( nqn(jo) == nseuil .and. lqn(jo) == lseuil .and. mseuil == 2 ) pop(jo) = pop(jo) - 1._db
      endif

    end do

    n_orb = jo        ! norbexc relativist

    if( ko > nnlm ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,130)
      end do
      stop
    endif

    charge = Z - sum( pop(1:n_orb) )
    if( charge > eps6 ) then
      dp = charge
      do io = 2,n_orb
        if( nqn(io) <= nseuil .or. lqn(io) > 1 .or. lqn(io-1) == lqn(io) ) cycle
        if( irel == 0 .or. lqn(io) == 0 ) then
          elmax = 2 + 4. * lqn(io)
          if( pop(io) > elmax - eps6 ) cycle
          ppp = pop(io)
          pop(io) = min( pop(io) + dp, elmax )
          dp = dp - pop(io) + ppp
        else
          elmax = 2 + 4. * lqn(io)
          if( sum(pop(io:io+1)) > elmax - eps6 ) cycle
          p = lqn(io) / ( 2 * lqn(io) + 1. )
          ppp = sum( pop(io:io+1) )
          pop(io) = min( pop(io) + p * dp, elmax*p )
          pop(io+1) = min( pop(io+1) + ( 1 - p ) * dp, elmax*(1-p) )
          dp = dp - sum( pop(io:io+1) ) + ppp
        endif
        if( dp < eps6 ) exit
      end do
      if( dp > eps6 ) then
        nmax = 0
        do io = 1,n_orb
          nmax = max(nqn(io),nmax)
        end do
        iaug = 0
        do io = 1,n_orb
          if( nqn(io) == nmax .and. lqn(io) == 1 ) then
            iaug = 1
            exit
          endif
        end do
        if( iaug == 1 ) then
          n_orb = n_orb + 1
          if( n_orb > nnlm ) then
            call write_error
            do ipr = 3,9,3
              write(ipr,130)
            end do
            stop
          endif
          nqn(n_orb) = nmax + 1
          lqn(n_orb) = 0
          pop(n_orb) = dp
          if( irel == 1 ) rqn(n_orb) = 0.5_db
        else
          n_orb = n_orb + 1
          nqn(n_orb) = nmax
          lqn(n_orb) = 1
          if( irel == 0 ) then
            pop(n_orb) = dp
          else
            rqn(n_orb) = 0.5_db
            pop(n_orb) = dp / 3
            n_orb = n_orb + 1
            if( n_orb > nnlm ) then
              call write_error
              do ipr = 3,9,3
                write(ipr,130)
              end do
              stop
            endif
            nqn(n_orb) = nmax
            lqn(n_orb) = 1
            rqn(n_orb) = 1.5_db
            pop(n_orb) = 2 * dp / 3
          endif
        endif
      endif

    elseif( charge < - eps6 ) then

      do io = n_orb,2,-1
        dp = charge
        if( nqn(io) <= nseuil .or. lqn(io) > 1 ) cycle
        if( irel == 0 .or. lqn(io) == 0 ) then
          pop(io) = pop(io) + dp
        else
          pop(io) = pop(io) + 2 * dp / 3
          pop(io-1) = pop(io-1) + dp / 3
        endif
        exit
      end do

    endif

    do io = 1,nlat(it)
      do j = 1,n_orb
        if(nqn(j) /= nvval(it,io) .or. lqn(j) /= lvval(it,io)) cycle
        if( irel == 0 .or. lqn(j) == 0 ) then
          popatv(it,io) = pop(j)
        else
          popatv(it,io) = pop(j) + pop(j+1)
        endif
        exit
      end do
    end do

    if( icheck > 0 ) then
      if( irel == 0 ) then
        write(3,165) Z
        do io = 1,n_orb
          if( Z > 18 .and. nqn(io) < 3 ) cycle
          if( Z > 36 .and. nqn(io) < 4 ) cycle
          if( Z > 86 .and. nqn(io) < 5 ) cycle
          write(3,160) nqn(io), lqn(io), pop(io)
        end do
      else
        write(3,175) Z
        do io = 1,n_orb
          if( Z > 18 .and. nqn(io) < 3 ) cycle
          if( Z > 36 .and. nqn(io) < 4 ) cycle
          if( Z > 86 .and. nqn(io) < 5 ) cycle
          write(3,180) nqn(io), lqn(io), rqn(io), pop(io)
        end do
      endif
    endif

    call dirac(E_total_exc,h_ray,icheck,ibav,irel,lqn,n_orb,n_ray,nnlm,nqn,pop,psi,psi_small,ray,ray_max,rho,rqn,Z)

! Recuperation de la fonction d'onde de l'orbitale de valence excite
    do io = 1, n_orb
      if( io <= n_coeur .or. lqn(io) /= l ) cycle
      if( irel == 1 .and. lqn(io) /= 0 ) then
        psi_open_val(1:n_ray,2) = 0.5_db * ( psi(1:n_ray,io) + psi(1:n_ray,io-1) )
        pop_open_val(2) = pop(io) + pop(io-1)
      else
        psi_open_val(1:n_ray,2) = psi(1:n_ray,io)
        pop_open_val(2) = pop(io)
      end if
    end do

    if( icheck > 0 ) write(3,190) E_total_exc * 2 * Rydb, E_total * 2 * Rydb, ( E_total_exc - E_total ) * 2 * Rydb

  endif   ! fin de la partie appel type_work pour atome excite

  rato(1:n_ray) = ray(1:n_ray)
  rhoit(1:n_ray,it) = rho(1:n_ray)    ! densite

! Les fonctions d'onde psi sont multipliees par r
!  rho_coeur(1:n_ray,it) = 0._db
!  do io = 1,n_coeur
!    rho_coeur(1:n_ray,it) = rho_coeur(1:n_ray,it) + pop(io) * psi(1:n_ray,io)**2
!  end do

!  rho_coeur(1:n_ray,it) = rho_coeur(1:n_ray,it) / ( quatre_pi * ray(1:n_ray)**2 )

  rho_coeur(1:n_ray,it) = 0._db
  do io = n_coeur+1,n_orb
    if( irel == 1 ) then
      rho_coeur(1:n_ray,it) = rho_coeur(1:n_ray,it) - pop(io) * ( psi(1:n_ray,io)**2 + psi_small(1:n_ray,io)**2 )
    else
      rho_coeur(1:n_ray,it) = rho_coeur(1:n_ray,it) - pop(io) * psi(1:n_ray,io)**2
    endif
  end do
  rho_coeur(1:n_ray,it) = rho_coeur(1:n_ray,it) / ( quatre_pi * ray(1:n_ray)**2 )
  rho_coeur(1:n_ray,it) = rho_coeur(1:n_ray,it) + rho(1:n_ray) 

! Extrapolation au centre de l'atome.
  p1 = rato(2) / ( rato(2) - rato(1) )
  p2 = 1 - p1
  rho_coeur(0,it) = p1 * rho_coeur(1,it) + p2 * rho_coeur(2,it)

! Oana: tous les tableaux que je fais sortir vers le main devraient etre remplis
! apres la construction de l'atome excite; si l'appel se fait a partir de atom cela
! ne gene pas, car on ne construit pas l'excite; si l'appel se fait a partir de type_work
! on aurait deja eu toutes les quantites pour it = itabs

! ici on stoque les fonctions d'onde (nonrelativistes, on tient pas compte du spin) pour les orbitales de valence

  do io = 1,nlat(it)
    do j = 1,n_orb
      if( nqn(j) /= nvval(it,io) .or. lqn(j) /= lvval(it,io) ) cycle
      if( irel == 0 .or. lqn(j) == 0 ) then
        psival(1:n_ray,io,it) = psi(1:n_ray,j)
      else
        do ir = 1,n_ray
          psival(ir,io,it) = 0.5_db * sum( psi(ir,j:j+1) )
        end do
      endif
      exit
    end do
  end do

! Recuperation des fonctions d'onde de la derniere orbitale de coeur
! et de la premiere orbitale de valence

  select case (Z)
    case(1,2)
       do i = 1,2
         psi_coeur(1:n_ray,i,it) = psi(1:n_ray,1)
         lcoeur(i,it) = lqn(1)
         ncoeur(i,it) = nqn(1)
       end do
    case default
       do i = 1,2
         if( i == 1 ) then
           io = n_coeur
         else
           io = n_coeur + 1
           if( irel == 1 .and. lqn(io) /= 0 ) io = io + 1
         endif
         lcoeur(i,it) = lqn(io)
         ncoeur(i,it) = nqn(io)
         if( irel == 1 .and. lqn(io) /=0 ) then
           psi_coeur(1:n_ray,i,it) = 0.5_db * ( psi(1:n_ray,io) + psi(1:n_ray,io-1) )
         else
           psi_coeur(1:n_ray,i,it) = psi(1:n_ray,io)
         end if
         psi_coeur(0,i,it) = 0._db
       end do
  end select

  if( icheck > 2 ) then
    write(3,135) it, ( ncoeur(i,it), lcoeur(i,it), i = 1,2 )
     do ir = 1, n_ray
        write(3,137) ray(ir)*bohr, ( psi_coeur(ir,i,it),i = 1,2 )
     end do
  end if
! on prend en compte la charge eventuelle de l'atome

  popatc(it) = 1._db * Z
  do io = 1,nlat(it)
    popatc(it) = popatc(it) - popatv(it,io)
  end do

  if( icheck > 1 ) then
     write(3,210) popatc(it)
     do l = 1,nlat(it)
       if( l == 1 ) write(3,215)
       write(3,220) nvval(it,l), lvval(it,l), popatv(it,l)
     end do

     if( it == itabs ) then
       write(3,230) ( nvval(it,l), lvval(it,l), l = 1,nlat(it) ), nseuil, lseuil
       do ir = 1,n_ray
         p = quatre_pi * rato(ir)**2
         if( nlat(it) > 0 ) then
           write(3,240) rato(ir) * bohr, p * rhoit(ir,it), p * rho_coeur(ir,it), ( psival(ir,l,it), l = 1,nlat(it) ), psii(ir,:)
         else
           write(3,240) rato(ir) * bohr, p * rhoit(ir,it), p * rho_coeur(ir,it), psii(ir,:)
         endif
       end do
     else
       write(3,230) ( nvval(it,l), lvval(it,l), l = 1,nlat(it) )
       do ir = 1,n_ray
         p = quatre_pi * rato(ir)**2
         if( nlat(it) > 0 ) then
           write(3,240) rato(ir) * bohr, p * rhoit(ir,it), p * rho_coeur(ir,it), ( psival(ir,l,it), l = 1,nlat(it) )
         else
           write(3,240) rato(ir) * bohr, p * rhoit(ir,it), p * rho_coeur(ir,it)
         endif
       end do
       
    endif
  endif

  deallocate( Ray )
  deallocate( rho )
  deallocate( psi, psi_small )

  return
  110 format(/' ---- Dirac --------',100('-')//' it =',i2,'  Z =',i4)
  115 format(//' The atomic number must be between Z = 1 and Z = 103', / &
               ' In n_orb_coeur, it is Z =',i12, ' !'//) 
  120 format(///'   n_ray =',i5,' > nrm =',i4,/ ' Use keyword n_ray in the indata file to decrease its value !')
  130 format(///'   n_orb > nnlm ',// ' Increase nnlm in lecture.f !')
  135 format(/'  it =',i3,// '      rato    psi_coeur(n=',i1,',l=',i1,')  ', 'psi_val(n=',i1,',l=',i1,')')
  137 format(1p,3e15.5)
  140 format(//' nlat =',i3,' > nlatm =',i3)
  150 format(/' Atom type',i3,',  Z =',i3, '   Non relativistic atomic calculation',/'   n  l     pop')
  160 format(i4,i3,f9.3)
  165 format(/' Excited atom, type 0',',  Z =',i3, '   Non relativistic atomic calculation',/'   n  l     pop')
  170 format(/' Atom type',i3,',  Z =',i3, '   Relativistic atomic calculation',/'   n  l   j      pop')
  175 format(/' Excited atom, type 0',',  Z =',i3, '   Relativistic atomic calculation',/'   n  l   j      pop')
  180 format(i4,i3,f5.1,2f9.3)
  190 format(/' Excited atom total energy =',1p,e13.5,' eV',/ &
              '         atom total energy =',e13.5,' eV',/&
              '                Difference =',e13.5,' eV') 
  210 format(/'  Popatc =',f9.4)
  215 format(/'  n  l  Popatv')
  220 format(2i3,f8.4)
  230 format(/'     rato      4*pi*r2*rho 4*pi*r2*rho_c',10(7x,2i3))
  240 format(1p,12e13.5)
end

!***********************************************************************

subroutine config(Z,irel,n_coeur,n_orb,nnlm,nqn,lqn,rqn,nel)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter(nzm=103)

  integer Z, n_orb, n_coeur
  integer, dimension(nnlm):: lqn, nqn

  real(kind=db), dimension(nnlm):: nel, rqn

  if( Z > nzm .or. Z <= 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,110) nzm, Z
    end do
    stop
  endif

  n_coeur = n_orb_coeur(Z)
  n_orb = n_orb_base(Z)

  do io = 1,n_orb
    ip = 0
    boucle_n: do n = 1,7
      do l = 0,n-1
        ip = ip + 1
        if( ip == io ) then
          nqn(io) = n
          lqn(io) = l
          exit boucle_n
        endif
      end do
    end do boucle_n
  end do

  Select Case(Z)
    Case(19,20)
      nqn(n_orb) = 4;   lqn(n_orb) = 0
    Case(21,22,23,24,25,26,27,28,29,30)
      nqn(n_orb-1) = 4; lqn(n_orb-1) = 0
      nqn(n_orb) = 3;   lqn(n_orb) = 2
    Case(37,38)
      nqn(n_orb) = 5; lqn(n_orb) = 0
    Case(39,40,41,42,43,44,45,46,47,48)
      nqn(n_orb-1) = 5; lqn(n_orb-1) = 0
      nqn(n_orb) = 4;   lqn(n_orb) = 2
    Case(49,50,51,52,53,54)
      nqn(n_orb-1) = 5; lqn(n_orb-1) = 0
      nqn(n_orb) = 5;   lqn(n_orb) = 1
    Case(55,56)
      nqn(n_orb-2) = 5; lqn(n_orb-2) = 0
      nqn(n_orb-1) = 5; lqn(n_orb-1) = 1
      nqn(n_orb) = 6;   lqn(n_orb) = 0
    Case(57) ! La ne converge pas avec 4f1
      nqn(n_orb-3) = 5; lqn(n_orb-3) = 0
      nqn(n_orb-2) = 5; lqn(n_orb-2) = 1
      nqn(n_orb-1) = 6; lqn(n_orb-1) = 0
      nqn(n_orb) = 5;   lqn(n_orb) = 2
    Case(58,59,60,61,62,63,64,65,66,67,68,69,70,71) ! converge pas
      nqn(n_orb-4) = 5; lqn(n_orb-4) = 0            ! avec 5d0
      nqn(n_orb-3) = 5; lqn(n_orb-3) = 1
      nqn(n_orb-2) = 6; lqn(n_orb-2) = 0
      nqn(n_orb-1) = 5; lqn(n_orb-1) = 2
      nqn(n_orb) = 4; lqn(n_orb) = 3
    Case(72,73,74,75,76,77,78,79,80)
      nqn(n_orb-1) = 6; lqn(n_orb-1) = 0
      nqn(n_orb) = 5;   lqn(n_orb) = 2
    Case(81,82,83,84,85,86)
      nqn(n_orb-1) = 6; lqn(n_orb-1) = 0
      nqn(n_orb) = 6;   lqn(n_orb) = 1
    Case(87,88)
      nqn(n_orb-2) = 6; lqn(n_orb-2) = 0
      nqn(n_orb-1) = 6; lqn(n_orb-1) = 1
      nqn(n_orb) = 7;   lqn(n_orb) = 0
    Case(89)
      nqn(n_orb-3) = 6; lqn(n_orb-3) = 0
      nqn(n_orb-2) = 6; lqn(n_orb-2) = 1
      nqn(n_orb-1) = 7; lqn(n_orb-1) = 0
      nqn(n_orb) = 6; lqn(n_orb) = 2
    Case(90,91,92,93,94,95,96,97,98,99,100,101,102,103)
      nqn(n_orb-4) = 6; lqn(n_orb-4) = 0
      nqn(n_orb-3) = 6; lqn(n_orb-3) = 1
      nqn(n_orb-2) = 7; lqn(n_orb-2) = 0
      nqn(n_orb-1) = 6; lqn(n_orb-1) = 2
      nqn(n_orb) = 5;   lqn(n_orb) = 3
  end select

  nel(1:n_orb) = 4 * lqn(1:n_orb) + 2._db  ! orbitale pleine
  nel(n_orb) = nel(n_orb) + Z - sum( nel(1:n_orb) )

  if( Z > 57 .and. Z < 67 ) then
    nel(n_orb) = Z - 57._db
    nel(n_orb-1) = 1._db
  elseif( Z > 66 .and. Z < 72 ) then
! 58:  5s5p6s ocupees + 2 el sur 5d
    nel(n_orb) = Z - 58._db
    nel(n_orb-1) = 2._db
  elseif( Z > 89 .and. Z <= nzm ) then
    nel(n_orb) = Z - 89._db
    nel(n_orb-1) = 1._db
  endif

  if( irel == 1 ) then
    n_orb_rel = n_orb

    do io = 1,n_orb
      if( lqn(io) /= 0 ) n_orb_rel = n_orb_rel + 1
    end do

! Cas calcul atomique relativiste
    n_coeur_rel = n_coeur
    do io = 1,n_coeur
      if( lqn(io) /= 0 ) n_coeur_rel = n_coeur_rel + 1
    end do
    n_coeur = n_coeur_rel
! la boucle doit etre placee avant le changement des configurations

! passage de la base l,s,ml,ms a l,s,j,mj

    jo = n_orb_rel + 1
    do io = n_orb,1,-1
      jo = jo - 1
      nqn(jo) = nqn(io)
      lqn(jo) = lqn(io)
      rqn(jo) = lqn(jo) + 0.5_db     ! mj
      if( lqn(jo) == 0 ) then
        nel(jo) = nel(io)
      elseif( abs( nel(io) - 4 * lqn(io) - 2 ) < eps10 ) then   ! cas d'une couche pleine
        nel(jo) = 2 * ( lqn(jo) + 1._db )
      else
        nel(jo) = ( ( lqn(jo) + 1._db ) / ( 2*lqn(jo) + 1._db ) ) * nel(io)
      endif
      if( lqn(jo) == 0 ) cycle
      jo = jo - 1
      nqn(jo) = nqn(io)
      lqn(jo) = lqn(io)
      rqn(jo) = lqn(jo) - 0.5_db     ! mj
      if( abs( nel(io) - 4 * lqn(io) - 2 ) < eps10 ) then
        nel(jo) = 2._db * lqn(jo)
      else
        nel(jo) = nel(io) - nel(jo+1)
      endif
    end do
    n_orb = n_orb_rel
  endif

  return
  110 format(//' The atomic number must be between Z = 1 and Z = ',i3, / &
               ' In subroutine config, it is Z =',i12, ' !'//) 
end


!***********************************************************************

! PROGRAM DIRAC

! DIRAC, A COMPUTER PROGRAM TO CARRY OUT BOTH NONRELATIVISTIC AND
! RELATIVISTIC SELF CONSISTENT FIELD CALCULATIONS FOR ATOMS AND IONS

! THE PROGRAM CONSISTS OUT OF THE FOLLOWING ROUTINES

!       DIRAC - CONTROL PROGRAM
!       DINPT - INPUT ROUTINE
!       DIPOT - GENERATES POTENTIAL
!       DIDIF - SOLVES DIRAC-SLATER DIFFERENTIAL EQUATIONS
!       DIDF1 - INTEGRATES INNER PART OF DIRAC SLATER EQUATIONS
!       HSDIF - SOLVES HARTREE-FOCK-SLATER DIFFERENTIAL EQUATIONS
!       DIOUT - PRINTED OUTPUT ROUTINE
!       DIADL - ROUTINE FOR NUMERICAL INTEGRATION

! WRITTEN AND/OR MODIFIED FROM LOS ALAMOS PROGRAM BY
! P. ROS AND D.E. ELLIS, CHEMISTRY DEPARTMENT
! FREE UNIVERSITY OF AMSTERDAM, HOLLAND.
! MOD.BY.ELLIS AND JANSEN(1977)FOR EFG AND DIPOLE INTEGRALS.
! MOD. BY G.A. BENESH FOR CONTINUUM ORBITALS..OCT78
! MODIFIED TO GENERATE BASIS FUNCTIONS AND CHARGE DENSITIES FOR THE
! NON-RELATIVISTIC AND RELATIVISTIC MOLECULAR PROGRAMS
! Modified by Tapio T Rantala to include different exchange and
! correlation potentials by implementing the routine XC(     ).
! 1985-11-06. The changes are done in DIPOT.
! A. Rosen found that it was some error for spinpolarized calc.
! This was corrected by Bengt Lindgren and A. Rosen Dec 1987.
! Modified by Y Joly 2000-2016

subroutine dirac(E_total,h_ray,icheck,ibav,irel,lqn,n_orb,n_ray,nnlm,nqn,pop,psi,psi_small,ray,ray_max,rho,rqn,Z)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  parameter( ahz = 1.e-10_db )

  integer:: Z
  integer, dimension(nnlm):: lqn, nqn
  integer, dimension(n_orb):: xl, xn

  logical:: last

  real(kind=db), dimension(n_ray,nnlm):: psi, psi_small
  real(kind=db), dimension(nnlm):: pop, rqn
  real(kind=db), dimension(n_orb):: xe, xj, xz
  real(kind=db), dimension(n_ray):: a, b, ray, rh, rho, rj, vr, Vrf, vs, y
  real(kind=db), dimension(5):: da, dbb, voc
  real(kind=db), dimension(9):: ha
  real(kind=db), dimension(15):: a0, b0

  psi(:,:) = 0._db
  psi_small(:,:) = 0._db

  Zn = 1._db * Z ! ATOMIC NUMBER
! IDIRC= ZERO FOR HARTREE-FOCK-SLATER, NONZERO FOR DIRAC-SLATER
  idirc = irel

  last = .false.

  call dinpt(Delta_rv,h_ray,ibav,icheck,idirc,lqn,n_orb,n_ray,nnlm, nqn,ph,pop,ray,ray_max,rh,rj,rqn,vr,vs, &
         xe,xj,xl,xn,xz,Zn,nc1,jspn,ha, rn,h,xion,phi,eps,del,xalph,rnuc,voc,anuc,sumel, rbar,rba2,vbar,h3)

! ITERATION TIE POINT HERE...
  do ncycl = nc1,0,-1

! RH=ELECTR.DENS. RJ=SPIN.DENS. VR=R*POT.
    call dipot(a,b,Convr,Delta_rv,E_v,E_Vxc,E_xc,E_Z,ibav, icheck,idirc,1,n_orb,n_ray,Q,ray,rh,rj,vr,vs,xte,y,Zn, &
          jspn,h,del,xalph,rnuc,da,dbb,voc,sumel,rbar,rba2,vbar,h3)

    if( Convr < Delta_rv .or. ncycl <= 0 ) then
      last = .true.
      Vrf(:) = Vr(:)
    endif

    Summa = 0._db
    ps = 1 - ph
    if( icheck > 2 .or. (icheck > 1 .and. last ) ) write(ibav,20) ncycl, ph
    rh(1:n_ray) = ps * rh(1:n_ray)

    do i = 1, n_orb

      nq = xn(i)
      lq = xl(i)
      fj = xj(i)
      e = xe(i)

      if( i == jspn + 1 ) then
        call dipot(a,b,Convr,Delta_rv,E_v,E_Vxc,E_xc,E_Z,ibav, icheck,idirc,2,n_orb,n_ray,Q,ray,rh,rj,vr,vs, xte,y,Zn, &
          jspn,h,del,xalph,rnuc,da,dbb,voc,sumel,rbar,rba2,vbar,h3)
        rj(1:n_ray) = ps * rj(1:n_ray)
      endif

! SKIP UNOCCUPIED STATES UNTIL LAST CYCLE.
      if ( (xz(i) < ahz) .and. .not. last ) cycle

      if( idirc == 0 ) then
        call hsdif(a,b,e,ibav,icheck,last,lq,n_ray,nq,ray,vr,y,Zn, npts,ha,rn,h,del,v0,da,dbb,a0,b0,voc,rbar,vbar)
        y(1:n_ray) = ph * xz(i) * a(1:n_ray)**2
      else
        call didif(a,b,e,ibav,icheck,last,lq,n_ray,nq,ray,vr,y,Zn, npts,ha,rn,h,eps,del,fj,v0,da,dbb,a0,b0,voc, &
                   vbar,fk,cs,g,q11,q22,tcs)
        y(1:n_ray) = ph * xz(i) * ( a(1:n_ray)**2 + b(1:n_ray)**2 )
      endif

      rh(1:n_ray) = rh(1:n_ray) + y(1:n_ray)
      if( i > jspn ) rj(1:n_ray) = rj(1:n_ray) + y(1:n_ray)
      Summa = Summa + xz(i) * e
      xe(i) = e

      if( .not. last) cycle

      psi(1:npts,i) = a(1:npts)
      if( idirc /= 0 ) psi_small(1:npts,i) = b(1:npts)

      if( icheck > 2 ) then
        write(ibav,24) v0, ( a0(k), k = 1,5 )
        if( idirc > 0 ) write(ibav,24) v0, (b0(k), k = 1,5)
      endif

    end do

    E_total = Summa - 0.5_db * (E_v - E_Z) - E_Vxc + E_xc

    if( icheck > 2 ) write(ibav,23) E_Vxc, E_v, E_Z, Summa, rh(1), E_total
    if( ncycl < 21 ) ph = 0.25_db * phi + 0.75_db * ph

    if( last ) exit

  end do

  call diout(E_total,E_v,E_Vxc,E_xc,E_Z,ibav,icheck, idirc,n_orb,n_ray,nnlm,psi,psi_small,Q,ray,rh,rho, &
             Summa,Vrf,xe,xj,xl,xn,xte,xz,Z)

  return
   20 format(' CYCLE ',i5,'  DENSITY AVERAGING ',f10.7)
   23 format(' E_Vxc, E_V, E_Z ,E_sum, Rho(1), E_total =',1p,7e14.6)
   24 format(1p,6e20.8)
end

! *********************************************************************

! INPUT ROUTINE FOR DIRAC PROGRAM
! FOR INPUT SPECIFICATIONS, SEE MAIN PROGRAM DIRAC
! n_ray : NUMBER OF POINTS IN RADIAL MESH
! n_orb : NUMBER OF SHELLS

subroutine dinpt(Delta_rv,h_ray,ibav,icheck,idirc,lqn,n_orb,n_ray, nnlm,nqn,ph,pop,ray,ray_max,rh,rj,rqn,vr,vs, &
               xe,xj,xl,xn,xz,Zn,nc1,jspn,ha, rn,h,xion,phi,eps,del,xalph,rnuc,voc,anuc,sumel, rbar,rba2,vbar,h3)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( adz = 0.3_db, aez = .25e4_db, afz = .1e-6_db, agz = .1e-4_db, arz = .208e-4_db)

  integer, dimension(nnlm):: lqn, nqn
  integer, dimension(n_orb):: xl, xn

  real(kind=db), dimension(nnlm):: pop, rqn
  real(kind=db), dimension(n_orb):: a, b, xe, xj, xz
  real(kind=db), dimension(n_ray):: ray, rh, rj, vr, vs
  real(kind=db), dimension(5):: voc
  real(kind=db), dimension(9):: ha

  nc1 = 199     ! MAXIMUM NUMBER OF CYCLES
  rn = ray_max  ! MAXIMUM RADIUS
  h = h_ray     ! INTERVAL OF LOGARITHMIC MESH
! ATOMIC WEIGHT, ANUC MUST BE SPECIFIED FOR PROPER VOLUME
! AVERAGING(DIDEN)
  anuc = 0._db
  xion = 0._db   ! IONIC CHARGE

  voc(2) = 0._db

! OPTION TO INSERT POTL WELL - VBAR NEGATIVE
  ph = 0.3_db  ! INITIAL DENSITY ITERATION AVERAGING FACTORS
  phi = 0.3_db ! FINAL DENSITY ITERATION AVERAGING FACTORS
  eps = 0._db  ! PARAMETER USED TO DETERMINE THE PRACTICAL INFINETY
! OF EACH ORBITAL ( EXP(- SQRT(EPS)) IS SET TO ZERO)
  del = 0._db   ! ABSOLUTE ACCURACY CRITERIUM FOR EIGENVALUES
  Delta_rv = 0.00001_db ! ABSOLUTE ACCURACY CRITERIUM FOR POTENTIAL

! PARAMETERS TO GENERATE EXTENDED BASIS SETS
  rbar = 0._db  ! RADIUS FOR THE POTENTIAL WELL  (ATOMIC UNITS )
  vbar = 0._db  ! VBAR = POTENTIAL WELL DEPTH ( ATOMIC UNITS )
  rba2 = 0._db
! RBA2 IS OUTER LIMIT OF FUNNEL WELL,DEFAULT IS RN
! SET DEFAULT VALUES OF  MISSING PARAMETERS

  if( rn < eps10 ) rn = 60._db
  if( phi < eps10 ) then
    phi = adz
    ph = phi
  endif
  if( eps < eps10 ) eps = aez
  if( del < eps10 ) del = afz
  if( rba2 < eps10 ) rba2 = rn
! THE ABOVE VALUES OF n_ray, RN,H, ETC. SUFFICE FOR MOST ATOMS AND IONS.
  if( h < agz ) h = 33._db - 0.1_db * Zn
  if( icheck > 1 ) then
    write(ibav,110)
    write(ibav,120) n_ray, n_orb, rn, h, Zn, xion, anuc
  endif
  if( icheck > 2 ) then
    write(ibav,130)
    write(ibav,140) ph, phi, eps, del, Delta_rv
    write(ibav,150) rbar, vbar, rba2
  endif

  rh(1:n_ray) = 0._db
  rj(1:n_ray) = 0._db
  jspn = 0
  if( jspn <= 0 ) jspn = n_orb
  if( nc1 == 0 ) nc1 = 30
  if( icheck > 2 ) then
    write(ibav,160)
    write(ibav,170) idirc, nc1, jspn
  endif

! XALPH = SLATER EXCHANGE COEFFICIENT
! XALPH = 2 GIVES HEDIN-LUNDQUIST EXCHANGE
  xalph = 2.01

! DEFINE H-CONSTANTS FOR DIFFER  START1,2  DIFF1,2
  h = 1 / h
  h3 = h / 3
  ha(1) = 8 * h3
  ha(2) = 1.125_db * h3
  ha(3) = 2.375_db * h3
  ha(4) = 0.625_db * h3
  ha(5) = 0.125_db * h3
  ha(6) = h3 / 240
  ha(7) = h3 / 120
  ha(8) = h3 / 80
  ha(9) = h3 / 60

! RNUC = NUCLEAR RADIUS IN BOHR RADII, IF SET TO 1. OR LARGER
!        RNUC IS COMPUTED AS  .0000208*ANUC**(1/3) Fermis
  rnuc = 0._db
  if( rnuc > 0._db ) then
    if( rnuc >= 1._db ) rnuc = arz * ( anuc**(1/3._db) )
    r = rn * exp(- (n_ray - 3) * h)
    if( rnuc < r ) then
      if( icheck > 2 ) write(ibav,180)
    endif
  endif
  sum = 0._db
  if( icheck > 2 ) then
    write(ibav,190) xalph, rnuc
    write(ibav,200)
  endif
  zbar = Zn

  do i = 1, n_orb
! XN = RADIAL QUANTUM NUMBER
! XL = ORBITAL ANGULAR MOMENTUM (OF MAJOR COMPONENT IN DIRAC-S
! XJ = TOTAL ANGULAR MOMENTUM ( NOT USED IN H-F-S, READ ZERO)
! XE = ORBITAL ENERGY ( OR ESTIMATE) IN ATOMIC UNITS
! XZ = NUMBER OF ELECTRONS IN SHELL
    xn(i) = nqn(i)
    xl(i) = lqn(i)
    xj(i) = rqn(i)
    xe(i) = 0._db
    xz(i) = pop(i)
    zbar = max( adz, zbar - xz(i) )
    if( abs(xe(i)) < afz ) xe(i) = - 0.5_db * (zbar / xn(i))**2 + vbar
    if( idirc == 0 ) xj(i) = 0._db
    if( icheck > 2 ) write(ibav,210) xn(i), xl(i), xj(i), xe(i), xz(i)
    if( xn(i) < 1-eps6 ) exit
    sum = sum + xz(i)
    if( xn(i) == xl(i) ) exit
    if( idirc == 0 ) cycle
    if( abs( abs(xl(i) - xj(i)) - 0.5_db ) > eps6 ) exit
    if( 2*xj(i) + 1 < xz(i)-eps6 ) exit
  end do

  if( i <= n_orb ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,220) i
      write(ipr,225) xn(i), xl(i)
      if( idirc /= 0 ) write(ipr,227) xj(i), xz(i)
    end do
    stop
  endif

  sumel = sum
  sum = Zn - sum
  if( icheck > 2 ) write(ibav,230) xion, sum
  xion = sum

  d = exp(h)
  r = rn / (d**n_ray)
  do k = 1, n_ray
    r = r * d
    ray(k) = r
    vr(k) = 0._db
    vs(k) = 0._db
  end do

  if( icheck > 2 ) write(ibav,240) ray(1), d, rnuc

! ANALYTIC APPROXIMATION TO CHARGE DENSITY
  if( rh(1) < eps10 ) then
    a(1:n_orb) = 0._db
    njm = 0
    do k = 1, n_orb
      b(k) = sqrt( 8 * abs(xe(k)) )
      njm = 2 * ( xn(k) - xl(k) )
      dfac = 1._db
      do nk = 1, njm
        fnk = nk
        dfac = dfac * fnk
      end do
      a(k) = ( b(k)**(njm + 1) ) / dfac
    end do
    do i = 1, n_ray
      rv = ray(i)
      sum = 0._db
      dfac = 0._db
      do k = 1, n_orb
        val = rv * b(k)
        if( val >= 32._db ) cycle
        njm = 2 * (xn(k) - xl(k))
        fnk = xz(k) * a(k) * exp(- val) * (rv**njm)
        if( k > jspn ) dfac = dfac + fnk
        sum = sum + fnk
      end do
      rj(i) = dfac
      rh(i) = sum
    end do
  endif

! SET STARTING SPIN DENSITY = .5*CHARGE DENS.
  if( abs( rj(1) ) < eps10 ) rj(1:n_ray) = 0.5_db * rh(1:n_ray)

  if( icheck > 3 ) then
    write(ibav,250)
    write(ibav,260) (rh(i),i = 1, n_ray)
    if( n_orb <= jspn ) return
    write(ibav,250)
    write(ibav,260) (rj(i),i = 1, n_ray)
  endif

  return
  110 format(/'   N    J    RN      H     Zn    XION   ANUC')
  120 format(2i5,5f7.3)
  130 format('  PH   PHI     EPS        DEL       Delta_rv')
  140 format(2f6.2,f12.2,2f12.8)
  150 format(/' PARAMTERS TO GENRT ADDTNL BASIS FNS ', ' RADIUS OF POTL WELL',f12.8,' ITS DEPTH',f12.8/,' CUTOFF',f12.8)
  160 format('  DIRC CYCL' )
  170 format(5i5)
  180 format(/' LESS THAN THREE MESH POINTS INSIDE NUCLEUS')
  190 format('  XALPHA = ',e20.8,', RNUC = ',e20.8)
  200 format('   XN   XL   XJ      XE          XZ')
  210 format(2i5,f5.1,2f10.3)
  220 format(//' Problem in Dirac routine ! :'//, ' ERROR ON EIGENVALUE CARD',i5)
  225 format(/' Check the values :'/,3x,' n = ',i3,/3x,' l = ',i3)
  227 format(/3x,' j = ',f7.3,/3x,' number of electron = ',f7.3)
  230 format(' XION was = ',f12.5,',and is = ',f12.5)
  240 format(' R(1) =',e15.7,', D =',e15.7,', RNUC =',e15.7)
  250 format(' CHARGE DENSITY')
  260 format(1x,8e15.7)
end

!***********************************************************************

! DIPOT COMPUTES THE SELF-CONSISTENT FIELD POTENTIAL FUNCTION FROM
! THE CHARGE DENSITY
! MODIFIED FOR RELATIVISTIC EXCHANGE...SEE 'SQK','QRELX'
! POINT ION ONLY OPTION FOR RHO(1).LT.0

subroutine dipot(a,b,Convr,Delta_rv,E_v,E_Vxc,E_xc,E_Z,ibav, icheck,idirc,isw,n_orb,n_ray,Q,ray,rh,rj,vr,vs,xte,y,Zn, &
       jspn,h,del,xalph,rnuc,da,dbb,voc,sumel,rbar,rba2,vbar,h3)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( pi4 = 4 * pi)

  real(kind=db), dimension(5):: da, dbb, voc
  real(kind=db), dimension(15):: dc
  real(kind=db), dimension(30):: uerr
  real(kind=db), dimension(n_ray):: a, b, c, ray, rh, rhd, rj, vr, vs, y, y2

  if( rh(1) <= 0._db ) then

    if( rh(2) <= 0._db ) then

      do k = 1, n_ray
!   FUNNEL = 0 IF NO FUNNEL POTENTIAL
        vr(k) = - Zn + funnel(ray(k),rbar,rba2,vbar)
      end do
    else
      vr(1:n_ray) = rh(1:n_ray)
    endif
    Convr = 0._db
    E_v = 0._db
    E_Vxc = 0._db
    E_Z = 0._db
!     NORMAL ENTRY

  else

    nerro = 0
    do k = 1, n_ray
      s2 = vr(k)
      vr(k) = vs(k)
      vs(k) = s2
    end do

    if( isw > 1 ) then

      y(1:n_ray) = y2(1:n_ray)
      E_Vxc = E_Vxc + diadl(rj,y,n_ray,h,2,xp)

    else

      s2 = log( rh(2) / rh(1) ) / h
      s3 = s2 + 1
      do k = 1,2
        r = ray(k)
        a(k) = (rh(k) * r) / s3
        b(k) = rh(k) / s2
        dbb(k) = h3 * rh(k)
        da(k) = dbb(k) * r
      end do
      nc = 0
      do k = 3,n_ray
        r = ray(k)
        dbb(3) = h3 * rh(k)
        da(3) = dbb(3) * r
        a(k) = ((a(k - 2) + da(3)) + (4 * da(2))) + da(1)
        b(k) = ((b(k - 2) + dbb(3)) + (4 * dbb(2))) + dbb(1)
        if( r <= rnuc ) nc = k
        da(1) = da(2)
        dbb(1) = dbb(2)
        da(2) = da(3)
        dbb(2) = dbb(3)
      end do

      nc = nc + 1
!...RENORMALIZE TOTAL CHARGE
      Q = a(n_ray)
      ren = sumel / Q
!...FORM AND SAVE COULOMB POTENTIAL...
      bn = b(n_ray)
      do k = 1, n_ray
        r = ray(k)
        if(r < rnuc ) then
          x = r / rnuc
          rvn = - Zn * x * (1.5_db - 0.5_db * x**2)
        else
          rvn = - Zn
          rh(k) = ren * rh(k)
        endif
        rhd(k) = ren * ( a(k) + r * ( bn - b(k) ) ) + rvn
      end do

      E_v = diadl(rh,rhd,n_ray,h,2,xp)
      E_Z = - Zn * bn

! OPTION FOR FUNNEL POTENTIAL
! FUNNEL=0 IF NO FUNNEL POTENTIAL
      do k = 1, n_ray
        rshft = funnel(ray(k),rbar,rba2,vbar)
        rhd(k) = rhd(k) + rshft
        a(k) = rshft
      end do
      xte = diadl(rh,a,n_ray,h,1,xp)

! INTERACTION POTENTIAL INSIDE NUCLEUS
      if( rnuc > 0._db ) then
        kk = nc + 2
        s5 = s2 + 3
        do k = 1, 2
          rr3 = rh(k) * (ray(k) ** 3)
          c(k) = rr3 / s5
          dc(k) = h3 * rr3
        end do
        do k = 3, kk
          dc(3) = (h3 * rh(k)) * (ray(k) ** 3)
          c(k) = ((c(k - 2) + dc(3)) + (4 * dc(2))) + dc(1)
          dc(1) = dc(2)
          dc(2) = dc(3)
        end do
!...CUBIC INTERPOLATION
        anc = 0._db
        bnc = 0._db
        cnc = 0._db
        do k = nc, kk
          rr3 = 1._db
          do i = nc, kk
            if( k /= i ) rr3 = rr3*(rnuc - ray(i)) / (ray(k)-ray(i))
          end do
          anc = anc + a(k) * rr3
          bnc = bnc + b(k) * rr3
          cnc = cnc + c(k) * rr3
        end do

        E_Z = E_Z - Zn * ( ( ( 1.5_db * anc / rnuc) - bnc ) - 0.5_db * cnc / rnuc**3 )
      endif

      do k = 1, n_ray
        pi4r2 = pi4 * ray(k)**2
        rhok = rh(k) / pi4r2
        if( n_orb <= jspn ) then
          spnk = 0.d0
        else
          spnk = (rh(k) - 2 * rj(k)) / pi4r2
        end if
        call xc(xalph, rhok, spnk, vxc1, vxc2, epsilon_xc)
        y(k) = vxc1 * ray(k)
        y2(k) = vxc2 * ray(k)
        a(k) = rh(k) - rj(k)
        b(k) = epsilon_xc * ray(k)
      end do

      if( n_orb <= jspn ) then
        E_Vxc = diadl(rh,y,n_ray,h,2,xp)
      else
        E_Vxc = diadl(a,y,n_ray,h,2,xp)
      end if
      E_xc = diadl(rh,b,n_ray,h,2,xp)

    end if

!     ASSEMBLE POTENTIAL , LOCATE MAX ERROR
    epslo = 0._db
    do k = 1, n_ray
      rv = rhd(k) + y(k)
      error = abs(rv - vr(k))
      iflag = (k / 20) - ((k - 1) / 20)
      if ( icheck > 3 .and. iflag == 1 ) then
        iflag1 = k / 20
        uerr(iflag1) = rv - vr(k)
      endif
      if( error > epslo ) then
        epslo = error
        nerro = k
      endif
      vr(k) = rv
    end do

    Convr = abs( epslo )

    if( icheck > 2 ) write(ibav,20) epslo, nerro
    if ( icheck > 3 ) then
      write(ibav,22) ( uerr(k), k = 2,15 )
      write(ibav,22) ( vr(20*k), k = 2,15 )
      if ( Convr < Delta_rv ) then
        write(ibav,23) rhd
        write(ibav,23) y
      endif
    endif

! COMPUTE POTENTIAL EXPANSION COEFFICIENTS FOR SMALL R

    del = min( epslo * 0.01_db, 0.1_db )

  endif

  ra = ray(1)
  rb = ray(2) / ra
  rc = ray(3) / ra

!     POINT NUCLEUS
  if( rnuc <= 0._db ) then
    voc(1) = - Zn
    ta = vr(1) + Zn
    tb = (vr(2) + Zn) / rb
    tc = (vr(3) + Zn) / rc
!     FINITE NUCLEUS
  else
    voc(1) = 0._db
    ta = vr(1)
    tb = vr(2) / rb
    tc = vr(3) / rc
  endif

  ra = 1._db
  deta = rb * rc * (rc - rb)
  detb = ra * rc * (ra - rc)
  detc = ra * rb * (rb - ra)
  det = deta + detb + detc
!     TEMPORARILY STORED IN VOC(5) UNTIL E CAN BE ADDED
  voc(5) = ( ta * deta + tb * detb + tc * detc ) / det
  deta = ra * ra
  detb = rb * rb
  detc = rc * rc
  voc(3) = ( ta * (detb - detc) + tb * (detc - deta) + tc * (deta - detb) ) / det
  voc(4) = (ta * (rc - rb) + tb * (ra - rc) + tc * (rb - ra)) / det
  if( idirc /= 0 ) voc(1:5) = - voc(1:5) * alfa_sf

  if( icheck > 2 ) write(ibav,21) voc(1), voc(5), voc(3:4)

  return
   20 format(' Max error in r.V(r) is ',e10.2,' at the ',i4,'th point')
   21 format(//' Potential expansion coefficients',/' V0C(1)= ',e18.8, &
  ', V0C(2)= ',e18.8,', V0C(3)= ',e18.8,', V0C(4)= ',e18.8/)
   22 format(1x,15e8.1)
   23 format(6x,8e14.7)
end

!***********************************************************************

! CALCULATES FUNNEL POTENTIAL R*SHIFT

function funnel(r, rbar, rba2, vbar)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter(abz = 1.5_db, agz = 0.5_db)

  funnel = 0._db
  if( (rbar <= 0._db) .or. (r > rba2) ) return
  funnel = vbar
  if( r > rbar ) funnel = ( vbar * (r - rba2) * (r - rba2) * ( r - abz * rbar + agz * rba2 ) ) &
           / ( ( rbar - rba2 ) * ( rbar - rba2 ) * ( agz * rba2 - agz * rbar ) )
  funnel = r * funnel

  return
end

!***********************************************************************

! INTEGRATES THE PRODUCT X*Y, OPTION TO CALC.END CORRECTION
! SIMPSON'S RULE INTEGRATION FROM R(1) TO INFINITY
! WEIGHTS H*(1 4 2 4 2 4 2 4 2 1)/3,INTEGRAND ASSUMED 0 PAST R(n_ray)
! INTEGRAL 0 TO R(1) MAY BE OMITTED OR CALCULATED BY 1 TERM.
! EXPONENT MAY BE SPECIFIED,OR COMPUTED BY FIT HERE.

function diadl(x, y, n_ray, h, nmx, xp )

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db), dimension(n_ray):: x, y

  select case(nmx)
    case(1)
      beg = 0._db
    case(2)
      dx = log( x(2) * y(2) / (x(1) * y(1)) ) / h
      beg = x(1) * y(1) / dx
    case(3)
      beg = x(1) * y(1) / xp
  end select

  sa = 0._db
  sb = 0._db
  do k = 2, n_ray, 2
    sa = sa + x(k-1) * y(k-1)
    sb = sb + x(k) * y(k)
  end do
  diadl = h * ( 2 * sa + 4 * sb - x(1) * y(1) ) / 3 + beg

  return
end

!***********************************************************************

subroutine xc(alfa, rho, spn, vxc1, vxc2, epsilon_xc)

!       RHO =TOTAL ELECTRON DENSITY, SPN =SPIN DENSITY (RHO(1)-RHO(2))
!       VXC1,VXC2 = EXCHANGE POTENTIALS
!       epsilon_xc = EXCHANGE ENERGY DENSITY, TOTAL ENERGY CONTRIBUTION IS
!             OBTAINED BY THE INTEGRAL (RHO*epsilon_xc)
!       CODE FOR SOME EXCHANGE-CORRELATION POTENTIALS AND ENERGIES
!         0 .LE. ALFA .LE.  1.   XALFA FOR ALFA
!                ALFA .EQ.  2.   VON BARTH - HEDIN
!                ALFA .EQ.  3.   VOSKO ET AL. - MANNINEN
!                ALFA .EQ. -1.   EXCHANGE ONLY
!                ALFA .LE. -2.   CORRELATION ONLY FOR ABS(ALFA)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  third = 1 / 3._db
  pi75 = 0.75_db / pi
  a = ( 4 / ( 9 * pi ) )**third
  pia = pi * a
  aa = 0.5_db**third
  aa1 = 1 / aa
  gamma = 4 * aa / ( 3 * (1 - aa ) )

!......VON BARTH-HEDIN - HEDIN-LUNDQVIST PARAMETERS:
  cp = 0.0225_db
  rp = 21._db
  cf = cp / 2
  rf = 2 * rp * aa1

  absxa = dabs(alfa)
  epsilon_xc = 0.d0
  vxc1 = 0.d0
  vxc2 = 0.d0
  if( rho <= 0 ) return
  s = spn / rho

! * * * EXCHANGE CONTRIBUTION * * *
! * * * A. PARAMAGNETIC

  rs = (pi75 / rho) ** third
  vxc1 = - 1 / ( pia * rs )
  vxc2 = vxc1

! * * * B. FERROMAGNETIC
  epsilon_xc = 0.75_db * vxc1

  if( abs( s )  > eps10 ) then
    s = max( s, -1._db )
    s = min( s,  1._db )
    s1 = 1 + s
    s2 = 1 - s
    if( s1 <= 0._db .or. s2 <= 0._db ) return
    s13 = s1**third
    s23 = s2**third
    gs = ((s1 * s13) + (s2 * s23)) * 0.5_db
    vxc1 = vxc1 * s13
    vxc2 = vxc2 * s23
    epsilon_xc = epsilon_xc * gs
!   *  EXCHANGE ONLY OPTION
  end if

! * * * CORRELATION CONTRIBUTION * * *
!   *  CORRELATION ONLY OPTION
  if( (alfa < 0._db) .and. (absxa < 2._db) ) return
  if( alfa < -1._db ) then
    vxc1 = 0._db
    vxc2 = 0._db
    epsilon_xc = 0._db
  end if

!   *  XALFA "CORRELATION"
  if (absxa < 2._db ) then
    xalfa = 1.5 * alfa
    vxc1 = xalfa * vxc1
    vxc2 = xalfa * vxc2
    epsilon_xc = xalfa * epsilon_xc
    return
  endif

! * * *  VON BARTH - HEDIN CORRELATION
! * * * A. PARAMAGNETIC

  if( absxa < 3._db ) then

    ecp = - cp * f_vonbarth(rs / rp)
    ucp = - cp * log(1 + (rp / rs))
    epsilon_xc = epsilon_xc + ecp
    vxc1 = vxc1 + ucp

! * * * B. FERROMAGNETIC

    vxc2 = vxc2 + ucp
    if( abs( s ) > eps10 ) then
      ecf = - cf * f_vonbarth(rs / rf)
      ucf = - cf * log(1 + (rf / rs))
      decf = ecf - ecp
      vc = gamma * decf
      tc = ucf - ucp - 4 * decf / 3._db
      fx = (gs - 1) / (aa1 - 1)
      epsilon_xc = epsilon_xc + decf * fx
      vxc1 = vxc1 + vc * (s13 - 1) + tc * fx
      vxc2 = vxc2 + vc * (s23 - 1) + tc * fx
    end if

    return
  endif

! * * *  VOSKO ET AL. - MANNINEN CORRELATION

  call ucor(rs, s, uc1, uc2, ec)
  epsilon_xc = epsilon_xc + ec * 0.5_db
  vxc1 = vxc1 + uc1 * 0.5
  vxc2 = vxc2 + uc2 * 0.5

  return
end

!***********************************************************************

subroutine ucor(rs, s, uc1, uc2, ec)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  data ap / .0621814d0 /
  data xp0 / -.10498d0 /
  data bp / 3.72744d0 /
  data cp / 12.9352d0 /
  data qp / 6.1519908d0 /
  data cp1 / 1.2117833d0 /
  data cp2 / 1.1435257d0 /
  data cp3 / -.031167608d0 /
  data af / .0310907d0 /
  data xf0 / -.325d0 /
  data bf / 7.06042d0 /
  data cf / 18.0578d0 /
  data qf / 4.7309269d0 /
  data cf1 / 2.9847935d0 /
  data cf2 / 2.7100059d0 /
  data cf3 / -.1446006d0 /

  third = 1 / 3._db
  fothi = 1 + third

  x = sqrt(rs)
  xpx = x * x + bp * x + cp
  xfx = x * x + bf * x + cf
  s4 = s**4 - 1
  fs = ((((1 + s)**fothi) + ((1 - s)**fothi)) - 2) / (2._db**fothi - 2)
  beta = 1._db / ( 2.74208_db + 3.182_db*x + 0.09873_db*x**2 + 0.18268_db*x**3)
  dfs = (fothi * (((1 + s)**third) - ((1 - s)**third))) / (2._db**fothi - 2)
  dbeta = - ((.27402_db * x + .09873_db + 1.591_db / x) * beta**2 )
  atnp = datan(qp / (2 * x + bp))
  atnf = datan(qf / (2 * x + bf))
  ecp = ap * ((log((x * x) / xpx) + (cp1 * atnp)) - (cp3 * (log(((x - xp0)**2) / xpx) + cp2 * atnp)))
  ecf = af * ((log((x * x) / xfx) + (cf1 * atnf)) - (cf3 * (log(((x - xf0)**2) / xfx) + cf2 * atnf)))
  ec = ecp + ((fs * (ecf - ecp)) * (1._db + s4 * beta))
  tp1 = (x * x + bp * x) / xpx
  tf1 = (x * x + bf * x) / xfx
  ucp = ecp - ((ap / 3) * (1._db - tp1 - (cp3 * ((x / (x - xp0) - tp1) - ((xp0 * x) / xpx)))))
  ucf = ecf - ((af / 3) * (1._db - tf1 - (cf3 * ((x / (x - xf0) - tf1) - xf0 * x / xfx))))
  uc0 = ucp + (ucf - ucp) * fs
  uc10 = uc0 - (ecf - ecp) * (s - 1) * dfs
  uc20 = uc0 - (ecf - ecp) * (s + 1) * dfs
  duc = (ucf - ucp) * beta * s4 * fs + (((((ecf - ecp) * (- rs / 3)) * dbeta) * s4) * fs)
  duc1 = duc - (ecf - ecp) * beta * (s - 1) * (4*(s**3)*fs + s4 * dfs)
  duc2 = duc - (ecf - ecp) * beta * (s + 1) * (4*(s**3)*fs + s4 * dfs)
  uc1 = uc10 + duc1
  uc2 = uc20 + duc2

  return
end

!**********************************************************************

!     SUBROUTINE FOR THE INTEGRATION OF THE HARTREE-FOCK-SLATER EQUATION
!     FINDS EIGENVALUES. NORMALIZES ORBITAL FUNCTIONS.
!     USES HAMMING'S PREDICTOR-CORRECTOR INTEGRATION SCHEME (MODIFIED)
!     MODIFIED FOR CONTINUUM STATES 1978....GREG BENESH...

subroutine hsdif(a,b,e,ibav,icheck,last,lq,n_ray,nq,ray,vr,y,Zn, npts,ha,rn,h,del,v0,da,dbb,a0,b0,voc,rbar,vbar)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter(delt = 0.01_db, delk = 0.01_db, acut = 1.e-10_db, &
     acz = 0.1e-4_db, adz = 1.2_db, aiz = 1.125_db, ajz = 0.125_db, &
     abb = 1.02_db, fmodc = .92561983_db, fcorc = .07438017_db, epx = 625._db, smset = 1.e20_db)

  Logical Continuum, last

  real(kind=db), dimension(5):: chix, ax
  real(kind=db), dimension(n_ray):: a, b, chi, ray, vr, y
  real(kind=db), dimension(5):: da, dbb, voc
  real(kind=db), dimension(9):: ha
  real(kind=db), dimension(15):: a0, b0

  de = 0._db
  eold = e
  smin = smset
!.....EXTRA SHIFT ACCOMODATES HYDROGENIC CASE.....
  emin = - 0.5_db * ( (Zn / nq)**2 ) + vbar - 0.5_db
  emax = - acz
  lqp = lq + 1
  lqpo = lqp + lq
  lqp2 = lqp + lqp
  nomax = nq - lqp

  if( e > emax ) then
    Continuum = .true.
  else
    Continuum = .false.
    if( abs(e - emax) < eps10 .or. e < emin ) e = 0.5_db * ( emax + emin )
  endif

  a0(1) = 1._db
  dg = exp(h * lqp)
!===>BEGIN OUTWARD INTEGRATION....
  rst = ( rn**lqp ) / ( dg**n_ray )

  do niter = 1,10

    ra = ray(1)
    iflag = 0
    voc(2) = voc(5) - e * ra

    b0(1:4) = 2 * voc(1:4)

    fis = 1._db
    lqpp = lqp2
    do is = 1,10
      suma = 0._db
      ismx = min0(is,4)
      do it = 1, ismx
        itm = is - it + 1
        suma = suma + b0(it) * a0(itm)
      end do
      a0(is+1) = ra * suma / (fis * lqpp)
      fis = fis + 1
      lqpp = lqpp + 1
    end do

! vr = r * V(r)
    do i = 1, n_ray
      r = ray(i)
      b(i) = 2 * r * ( vr(i) - r*e )
      if( Continuum ) chi(i) = - ( lq*lqp + b(i) )
    end do

    do i = 1,4
      r = ray(i)
      ta = r / ra
      suma = 0._db
      sumb = 0._db
      fis = 0._db
      tpow = 1._db
      do is = 1, 11
        suma = suma + a0(is) * tpow
        sumb = sumb + fis * a0(is) * tpow
        fis = fis + 1._db
        tpow = tpow * ta
      end do
      a(i) = suma
      da(i) = sumb
      dbb(i) = - lqpo * da(i) + b(i) * a(i)
    end do

    nodes = 0

    m = n_ray - 10
    do k = 11, m
      km = n_ray - k
      if( b(km) < 0._db ) exit
    end do

    if( Continuum ) then

      nm3 = n_ray - 3
      do is = 1, nm3
        if (chi(is) > 0._db) exit
      end do
      do i = is, n_ray
        if( chi(i) >= 0._db ) then
          chi(i) = sqrt( chi(i) )
          cycle
        endif
        if( icheck > 3 ) write(ibav,3000) i, chi(i)
        chi(i) = 0._db
      end do

      if( icheck > 3 ) write(ibav,15) is, chi(is), ray(is)
      km = n_ray - 1

    endif

    ki = km + 1
    adif = 0._db
    bdif = 0._db

    do k = 5, ki
      r = ray(k)
      bpred = da(1) + ha(1) * ( dbb(4) - 0.5_db * dbb(3) + dbb(2) )
      apred = a(k-4) + ha(1) * ( da(4) - 0.5_db * da(3) + da(2) )
      bmode = bpred + fmodc * bdif
      amode = apred + fmodc * adif
      fmode = - lqpo * bmode + b(k) * amode
      bcorr = aiz * da(4) - ajz * da(2) + ha(2) * ( fmode + 2 *  dbb(4) - dbb(3) )
      acorr = aiz * a(k-1) - ajz * a(k-3) + ha(2) * ( bmode + 2 * da(4) - da(3) )
      adif = acorr - apred
      bdif = bcorr - bpred
      a(k) = acorr - fcorc * adif
      da(5) = bcorr - fcorc * bdif
      dbb(5) = - lqpo * da(5) + b(k) * a(k)
      da(1) = da(2)
      da(2) = da(3)
      da(3) = da(4)
      da(4) = da(5)
      dbb(2) = dbb(3)
      dbb(3) = dbb(4)
      dbb(4) = dbb(5)
      asq = a(k) * a(k-1)
      if( asq < 0._db ) nodes = nodes + 1
    end do

!===>CONTINUUM STATE BY WKB ,FITTED TO OUTWARD INTEGRATION...
    if( Continuum ) then

      b(is) = 0._db
      b(is + 1) = (h / 24) * (((9 * chi(is) + 19 * chi(is + 1)) - (5 * chi(is + 2))) + chi(is + 3))

      in = is + 2
      do i = in, n_ray
        b(i) = ( h / 3 ) * ( chi(i-2) + 4*chi(i-1) + chi(i) ) + b(i-2)
      end do

      do i = is, n_ray
        if( chi(i) <= 0._db ) then
          if( icheck > 2 ) write(ibav,3001) i, chi(i)
          chi(i) = 0._db
        else
          chi(i) = sqrt(chi(i) / ray(i))
        endif
      end do

      nm4 = n_ray - 4
      do i = is, nm4
        r1 = ray(i)**lqp * a(i) * chi(i)
        r2 = ray(i+2)**lqp * a(i+2) * chi(i+2)
        dp = b(i+2) - b(i)
        arg = atan( ( 1 / tan(dp) ) - (r2 / r1) / sin(dp) )
        beth = b(i) - arg
        facts = r1 / cos(b(i) - beth)
        i4 = i + 4
        sumsq = 0._db
        rg = rst * dg**i
        do k = i, i4
          chix(k+1-i) = (facts / chi(k)) * cos( b(k) - beth ) / ray(k)
          ax(k+1-i) = rg * a(k) / ray(k)
          rg = rg * dg
          sumsq = sumsq +  ( chix(k+1-i) - ax(k+1-i) )**2
        end do
        if( i /= is .and. sumsq >= smin ) cycle
        smin = sumsq
        fact = facts
        isi = i
        beta = beth
        argg = arg
        dpp = dp
        r11 = r1
        r21 = r2
      end do

      if( icheck > 3 ) then
        write(ibav,11) r11, r21
        write(ibav,12) fact, argg, dpp
        write(ibav,14) beta, isi
      endif

      rg = rst * dg**isi
      do k = isi, n_ray
        a(k) = (fact / chi(k)) * cos( b(k) - beta )
      end do
      kj = isi - 1

    else

      if( icheck > 3 ) write(ibav,18) nq, lq, emin, e, emax

!     TOO MANY NODES
      if( nomax < nodes ) then
        if( e < emax ) emax = e
        e = e + 0.5_db * max(e + e,emin - e)

        dl1 = abs(emax - emin) / (abs(e) + 1._db)
!===>ERROR,TIGHT BOUNDS,BUT WRONG NR OF NODES....
        if( dl1 >= delk ) cycle

        if( icheck > 2 ) write(ibav,20) nodes, nomax
        if( icheck > 2 ) write(ibav,18) nq, lq, emin, e, emax
        emin = adz * emin
        cycle

!     TOO FEW NODES
      elseif( nomax > nodes ) then

        if( e > emin ) emin = e
        e = 0.5_db * (e + emax)
        dl1 = abs( emax - emin ) / (abs(e) + 1._db)
        if( dl1 >= delk ) cycle

        if( icheck > 2 ) write(ibav,20) nodes, nomax
        if( icheck > 2 ) write(ibav,18) nq, lq, emin, e, emax
        if( e >= - acz ) then
          if( icheck > 2 ) write(ibav,19) nq, lq
          exit
        endif
        e = emax
        cycle

      endif

! CORRECT NO.OF NODES. BEGIN INWARD INTEGRATION.

      ra = a(ki)
      rb = da(5)
      do k = n_ray, ki+5,-1
        if( epx > b(k) ) exit
        a(k) = 0._db
      end do
      kj = k

! AYMPTOTIC SOLUTION FOR OUTER 4 POINTS...
!  CORR OCT79,NOTE R**(L+1) SCALING ON A(K) AND DP/DT, WHERE T=LN(R)...
      do i = 1, 4
        if( b(k) < 0._db ) then
          if( icheck > 2 ) write(ibav,3003) k, b(k)
          b(k) = 0._db
        endif
        r = ray(k)
        alfr = sqrt( b(k) )
        a(k) = exp(- alfr) / (r**lqp)
        da(i) = - ( lqp + alfr * r ) * a(k)
        dbb(i) = - lqpo * da(i) + b(k) * a(k)
        k = k - 1
      end do
      adif = 0._db
      bdif = 0._db

      kk = k
      do k = kk,ki,-1
        bpred = da(1) - ha(1) * ( dbb(4) - 0.5_db * dbb(3) + dbb(2))
        apred = a(k+4) - ha(1) * ( da(4) - 0.5_db * da(3) + da(2) )
        bmode = bpred + fmodc * bdif
        amode = apred + fmodc * adif
        fmode = - lqpo * bmode + b(k) * amode
        bcorr = aiz * da(4) - ajz * da(2) - ha(2) * ( fmode + 2 * dbb(4) - dbb(3) )
        acorr = aiz * a(k+1) - ajz * a(k+3) - ha(2) * ( bmode + 2 * da(4) - da(3) )
        adif = acorr - apred
        bdif = bcorr - bpred
        a(k) = acorr - fcorc * adif
        da(5) = bcorr - fcorc * bdif
        dbb(5) = - lqpo * da(5) + b(k) * a(k)
        if( k == ki ) exit
        da(1) = da(2)
        da(2) = da(3)
        da(3) = da(4)
        da(4) = da(5)
        dbb(2) = dbb(3)
        dbb(3) = dbb(4)
        dbb(4) = dbb(5)
      end do

      ra = ra / a(ki)
      a(ki:kj) = a(ki:kj) * ra
      rc = da(5) * ra

    endif

    rg = rst
    do k = 1, kj
      rg = rg * dg
      a(k) = a(k) * rg
    end do
    y(1:n_ray) = a(1:n_ray)**2

    xp = 2*lqp + 1._db

!===>CONTINUUM STATE NORMALIZED ON SPHERE RADIUS RBAR....
    if( Continuum) then
      do k = 1, n_ray
        if( ray(k) < rbar ) cycle
        nrb = k
        exit
      end do
      w = diadl(ray,y,nrb,h,3,xp)
      exit  ! On sort de la boucle sur les iterations
    endif

    w = diadl(ray,y,n_ray,h,3,xp)

!===>SET NEW ENERGY AND ADJUST BOUNDS....
    de = - 0.5_db * (ray(ki)**lq) * a(ki) * (rc - rb) / w

    if( de > 0._db ) then
      emin = e
    elseif( de < 0._db ) then
      emax = e
    endif

    chkbd = abs(de) + ajz * (e - acz)
    if( chkbd > 0._db ) then
      iflag = 1
      if( icheck > 3 )  write(ibav,21) nq, lq
    endif

!===>TEXT FOR ENERGY CONVERGENCE....
    e = e + de
    dell = del
    if( .not. last ) dell = max( dell, delt * abs(eold - e) )
    if( icheck > 3 ) write(ibav,18) nq, lq, emin, e, emax, de, dell

    if( e >= emax ) then
      e = 0.5_db * ( e - de + emax)
      if( iflag == 1 ) e = ajz * e + (1 - ajz) * emax
    else
      if( e <= emin ) e = 0.5_db * ( e - de + emin )
    endif

    if( ( dell > abs(de) ) .and. ( iflag == 0) ) exit

  end do ! fin de la boucle sur les iterations

  if( icheck > 2 .and. niter > 10 ) write(ibav,13) nq, lq, emin, e, emax

  if (w <= 0._db) then
    if( icheck > 2 ) write(ibav,3005) w
    w = 1._db
  endif

  x = 1 / sqrt( w )
  npts = 0

  do k = 1, n_ray
    a(k) = a(k) * x
    if( abs(a(k)) > acut ) npts = npts + 1
  end do

  v0 = x

  return
   11 format(' R1 AND R2 ARE',2f10.5)
   12 format(//' FACTS,ARG,DP ARE',3f10.5)
   13 format(' 10 E TRYALS EXCEEDED, PROCEEDING  ',2i5,4x,4e15.5)
   14 format(' THE FIRST BETH IS',f10.5,' AT I EQUAL TO',i5)
   15 format(' IS,CHI,ray ARE',i5,2f10.7)
   18 format(2i4,4x,e15.5,e19.9,e15.5,e13.3,e13.3)
   19 format(2i4,'...LEVEL NOT BOUND,PROCEEDING..BEWARE')
   20 format(' EMAX-EMIN TOO SMALL, NUMBER OF NODES FOUND=',i5, ' NUMBER  NODES REQUIRED=',i5/)
   21 format(' *** DE LARGE ***   ',2i4)
 3000 format(' ]]]]] TP1 ',i5,1p,e15.5)
 3001 format(' ]]]]] TP2 ',i5,1p,e15.5)
 3003 format(' ]]]]] TP3 ',i5,1p,e15.5)
 3005 format(' ]]]]] TP4 ',1p,e15.5)
end

!***********************************************************************

!  CONTROL SUBROUTINE FOR THE INTEGRATION OF THE DIRAC-SLATER EQUATIONS.
!  FINDS EIGENVALUES. NORMALIZES ORBITAL FUNCTIONS.

subroutine didif(a,b,e,ibav,icheck,last,lq,n_ray,nq,ray,vr,y,Zn, npts,ha,rn,h,eps,del,fj,v0,da,dbb,a0,b0,voc, &
                   vbar,fk,cs,g,q11,q22,tcs)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter( ga1 = 251._db, ga2 = 646._db, ga3 = 264._db, ga4 = 106._db, ga5 = 19._db, ga6 = 116._db, ga7 = 496._db, &
   ga8 = 96._db, ga9 = 16._db, ga10 = 4._db, ga11 = 81._db, ga12 = 306._db, ga13 = 216._db, ga14 = 126._db, ga15 = 9._db, &
   ga16 = 56._db, ga17 = 256._db, ga18 = 96._db, ga19 = 256._db, ga20 = 56._db )

  parameter( acut = 1e-10_db )
  parameter( delt = 0.1_db, delk = 0.01_db, efac = 0.6_db )
  parameter( aaz = 1.2_db, acz = .1e-4_db, aiz = 1.e-33_db )

  logical last

  real(kind=db), dimension(5):: da, dbb, voc
  real(kind=db), dimension(9):: ha
  real(kind=db), dimension(15):: a0, b0
  real(kind=db), dimension(n_ray):: a, b, ray, vr, y

  c = 1 / alfa_sf

  eold = e
  ntest = 0
  etest = del * (1 - e)
  dels = etest
  de = 0._db
  emin = - 0.5_db * ( 1 + Zn * alfa_sf ) * (Zn / nq)**2 + vbar
  emax = - acz
  if( e > emax .or. e < emin ) e = 0.5_db * (emax + emin)
!   PREPARE CONSTANTS USED FOR THIS ORBITAL
  fk = 2 * (lq - fj) * (fj + 0.5_db)
  cs = c
  tcs = 2 * c
  q21 = voc(1)
  g = sqrt( fk**2 - q21**2 )
  q11 = - g - fk
  q22 = - g + fk
  if( fk < 0._db ) then
    a0(1) = g - fk
    b0(1) = q21
  else
    a0(1) = - q21
    b0(1) = g + fk
  endif
  dg = exp(h * g)

  rst = (rn**g) / (dg**n_ray)

  do

    call didf1(a,b,e,ki,lq,n_ray,nb_nodes,ray,vr, ha,da,dbb,a0,b0,voc,fk,cs,g,q11,q22,tcs)

!  nb_nodes = Number of nodes + l + 1. nb_nodes SHOULD EQUAL nq.
    if( icheck > 3 ) write(ibav,10) nq, lq, fj, emin, e, emax

!  TOO MANY NODES.
    if( nq < nb_nodes) then
      if( e < emax ) emax = ( 2 * emax + e) / 3
      dele = max(e,emin - e)
      e = e + (0.5_db * dele)
      dl1 = abs(emax - emin) / (abs(e) + 1)

      if( dl1 < delk ) then
        if( icheck > 2 ) write(ibav,15) nb_nodes, nq, lq, fj, emin, e, emax
        emin = aaz * emin
      endif

      cycle

    elseif( nq > nb_nodes ) then
!  TOO FEW NODES.
      if( e > emin ) emin = e
      e = 0.5_db * (e + emax)
      dl1 = abs(emax - emin) / (abs(e) + 1)

      if( dl1 >= delk ) cycle

      if( icheck > 2 ) write(ibav,50) nb_nodes, nq, lq, fj, emin, e, emax

      if( ntest /= 0 ) then
        call write_error
        do ipr = 3,9,3
          write(ipr,60)
        end do
      stop
      endif

      ntest = 1
      emax = - acz
      e = 0.5_db * (e + emax)

      cycle
    endif

!  Correct number of nodes
    ra = a(ki)
    rb = b(ki)
    do i = 1, ki
      a(i) = a(i) / ra
      b(i) = b(i) / ra
    end do

!  STARTS INWARD INTEGRATION. OUTER BOUNDARY CONDITION IS SET HERE.
    kj = ki + 4
    kk1 = ki + 4

    if( kk1 > n_ray ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,70) ki, n_ray
      end do
      stop
    endif

    do k = kk1, n_ray
      r = ray(k)
      if( eps < r * ( vr(k) - e * r ) ) exit
    end do

!  ISOLATED ATOM BOUNDARY CONDITIONS.
    kj = k - 1
    rz = - vr(kj) / r

    do k = kj, n_ray
      a(k) = 0._db
      b(k) = 0._db
    end do
    rl = (lq + 0.5_db) / r
    rk = fk / r
    p = - 2 * ( e + rz ) + rl * rl
    if( p < 0._db ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,19)
      end do
      stop
    endif
    apa = - sqrt(p) + 0.5_db * (rz - rl * rl) / (r * p)
    a(kj) = 1._db
    b(kj) = cs * (apa + rk) / (- cs * tcs - e - rz)
    do l = 1,4
      k = kj - l
      a(k) = a(kj)
      b(k) = b(kj)
    end do

    do i = 1, 4
      k = kj + 1
      do l = 1, 5
        k = k - 1
        rp21 = ( e * ray(k) - vr(k) ) / cs
        rp12 = - rp21 - tcs * ray(k)
        da(l) = q11 * a(k) + rp12 * b(k)
        dbb(l) = q22 * b(k) + rp21 * a(k)
      end do
      a(kj-1) = a(kj) - ((ga1 * da(1) + ga2 * da(2) - ga3 * da(3) + ga4 * da(4) - ga5 * da(5)) * ha(6))
      b(kj-1) = b(kj) - ((ga1 * dbb(1) + ga2 * dbb(2) - ga3 * dbb(3) + ga4 * dbb(4) - ga5 * dbb(5)) * ha(6))
      a(kj-2) = a(kj) - ((ga6 * da(1) + ga7 * da(2) + ga8 * da(3) + ga9 * da(4) - ga10 * da(5)) * ha(7))
      b(kj-2) = b(kj) - ((ga6 * dbb(1) + ga7 * dbb(2) + ga8 * dbb(3) + ga9 * dbb(4) - ga10 * dbb(5)) * ha(7))
      a(kj-3) = a(kj) - ((ga11 * da(1) + ga12 * da(2) + ga13 * da(3) + ga14 * da(4) - ga15 * da(5)) * ha(8))
      b(kj-3) = b(kj)-((ga11 * dbb(1) + ga12 * dbb(2) + ga13 *dbb(3) + ga14 * dbb(4) - ga15 * dbb(5)) * ha(8))
      a(kj-4) = a(kj) - ((ga16 * da(1) + ga17 * da(2) + ga18 * da(3) + ga19 * da(4) + ga20 * da(5)) * ha(9))
      b(kj-4) = b(kj)-((ga16 * dbb(1) + ga17 * dbb(2) + ga18 *dbb(3) + ga19 * dbb(4) + ga20 * dbb(5)) * ha(9))
    end do

    k = kj - 3
    da(1) = da(2)
    dbb(1) = dbb(2)
    da(2) = da(3)
    dbb(2) = dbb(3)
    da(3) = da(4)
    dbb(3) = dbb(4)

    do k = kj-4,ki,-1
      rp21 = ( e * ray(k) - vr(k) ) / cs
      rp12 = - rp21 - tcs * ray(k)
      akk = a(k+4) - ha(1) * (da(3) - 0.5_db * da(2) + da(1))
      bkk = b(k+4) - ha(1) * (dbb(3) - 0.5_db * dbb(2) + dbb(1))
      da(4) = q11 * akk + rp12 * bkk
      dbb(4) = q22 * bkk + rp21 * akk
      a(k) = a(k+1) - ha(2) * da(4) - ha(3) * da(3) + ha(4) * da(2) - ha(5) * da(1)
      b(k) = b(k+1) - ha(2) * dbb(4) - ha(3) * dbb(3) + ha(4) * dbb(2) - ha(5) * dbb(1)
      da(1) = da(2)
      dbb(1) = dbb(2)
      da(2) = da(3)
      dbb(2) = dbb(3)
      da(3) = (q11 * a(k)) + (rp12 * b(k))
      dbb(3) = (q22 * b(k)) + (rp21 * a(k))
    end do

    rc = a(ki)
    rb = rb / b(ki)
    do k = ki, kj
      a(k) = a(k) / rc
      b(k) = b(k) / rc
    end do
    rg = rst
    do k = 1, kj
      rg = rg * dg
      a(k) = a(k) * rg
      b(k) = b(k) * rg
    end do
    do k = 1, n_ray
      y(k) = a(k)**2 + b(k)**2
    end do

! ESTIMATE CHANGE IN ENERGY FOR CONTINUOUS DERIVATIVE......
    w = diadl(ray,y,n_ray,h,3,2*g+1._db)
    de = cs * a(ki) * b(ki) * ( 1 - (rb * rc / ra) ) / w

    do

      dl1 = abs( emax - emin ) / ( abs(e) + 1 )
      dl = abs( de )
      dll = min( dl, - (0.5_db * e) )

      if( dl <= delk .or. dl1 >= delk ) exit

! Reduce 'de' and try again
      if( icheck > 2 ) write(ibav,40) nq, lq, fj, emin, e, emax, de
      de = 0.5_db * de

    end do

    if( de > 0._db ) then
      emin = e
      de = dll
    elseif( de < 0._db ) then
      emax = ( 2 * e + emax) / 3
      de = - dll
    endif

    e = e + de
    dell = del
    if( .not. last ) dell = max(dell,delt * abs(eold - e),dels)
    if( icheck > 3 ) write(ibav,10) nq, lq, fj, emin, e,emax, de, dell

    if( dl >= dell ) then
      if( e > emax ) e = 0.5_db * ( e - de + emax )
      if( e < emin ) e = 0.5_db * ( e - de + emin )
      cycle
    endif

    exit

  end do

  x = 1 / sqrt(w)
  npts = 0
  do k = 1, n_ray
    a(k) = a(k) * x
    if( abs( a(k) ) > acut ) npts = npts + 1
    b(k) = b(k) * x
  end do
  v0 = x / ra

  return

   10 format(2i5,f5.1,e15.5,e19.9,e15.5,2e13.3)
   15 format(' Too many nodes',i4,2i4,f4.1,3e15.7)
   19 format(//' Dirac routine:'/, ' THE TURNING POINT IS BEYOND RN. A LARGER VALUE OF RN', ' MAY BE NEEDED.')
   40 format(//' Dirac routine:'/,' DE too large',2i4,f4.1,4e15.7)
   50 format(//' Dirac routine:'/,' Too few nodes',i4,2i4,f4.1,3e15.7)
   60 format(//' Dirac routine:'/,' ntest /= 0')
   70 format(//' Dirac routine:'/,' NO ROOM FOR TAIL..',2i5)
end

!***********************************************************************

!  INTEGRATES OUTWARD TO CLASSICAL TURNING POINT, COUNTS NODES IN THE
!  MAJOR COMPONENT OF THE RADIAL WAVE FUNCTION.

subroutine didf1(a,b,e,ki,lq,n_ray,nb_nodes,ray,vr, ha,da,dbb,a0,b0,voc,fk,cs,g,q11,q22,tcs)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db), dimension(5):: da, dbb, voc
  real(kind=db), dimension(9):: ha
  real(kind=db), dimension(15):: a0, b0
  real(kind=db), dimension(n_ray):: a, b, ray, vr

!  STARTS OUTWARD INTEGRATION WITH A POWER SERIES.
  ra = ray(1)
  tcr = tcs * ra
  voc(2) = voc(5) + e * ra / cs
  fis = 0._db
  do is = 2, 11
    fis = fis + 1
    suma = 0._db
    sumb = 0._db
    ismx = min0(is,4)
    do it = 2, ismx
      itm = (is - it) + 1
      suma = suma - (voc(it) * ((b0(itm) * (fis + g - fk)) - (voc(1) * a0(itm))))
      sumb = sumb + (voc(it) * ((a0(itm) * (fis + g + fk)) - (voc(1) * b0(itm))))
    end do
    fac = fis * (fis + g + g)
    a0(is) = ((- (tcr * (fis + g - fk) * b0(is - 1))) + suma) / fac
    b0(is) = ((- (tcr * voc(1) * b0(is - 1))) + sumb) / fac
  end do
  do k = 1, 5
    r = ray(k)
    ta = r / ra
    suma = 0._db
    sumb = 0._db
    tpow = 1._db
    do is = 1, 11
      suma = suma + a0(is) * tpow
      sumb = sumb + b0(is) * tpow
      tpow = tpow * ta
    end do
    a(k) = suma
    b(k) = sumb
    rp21 = (e * r - vr(k)) / cs
    rp12 = - rp21 - tcs * r
    da(k) = q11 * a(k) + rp12 * b(k)
    dbb(k) = q22 * b(k) + rp21 * a(k)
  end do

  da(1) = da(2)
  dbb(1) = dbb(2)
  da(2) = da(3)
  dbb(2) = dbb(3)
  da(3) = da(4)
  dbb(3) = dbb(4)
  nb_nodes = 1 + lq
  m = n_ray - 10
  do k = 11,m
    km = n_ray - k
    if( e*ray(km) > vr(km) ) exit
  end do
  ki = km + 1
  do k = 5, ki
    rp21 = (e * ray(k) - vr(k)) / cs
    rp12 = - rp21 - tcs * ray(k)
    akk = a(k-4) + ha(1) * (da(3) - 0.5_db * da(2) + da(1))
    bkk = b(k-4) + ha(1) * (dbb(3) - 0.5_db * dbb(2) + dbb(1))
    da(4) = q11 * akk + rp12 * bkk
    dbb(4) = q22 * bkk + rp21 * akk
    a(k) = a(k-1) + ha(2)*da(4) + ha(3)*da(3) - ha(4)*da(2) + ha(5) * da(1)
    b(k) = b(k-1) + ha(2)*dbb(4) + ha(3)*dbb(3) - ha(4)*dbb(2) + ha(5)*dbb(1)
    da(1) = da(2)
    dbb(1) = dbb(2)
    da(2) = da(3)
    dbb(2) = dbb(3)
    da(3) = q11 * a(k) + rp12 * b(k)
    dbb(3) = q22 * b(k) + rp21 * a(k)
    if( a(k)*a(k-1) < 0._db ) nb_nodes = nb_nodes + 1
  end do

  return
end

!***********************************************************************

!   DIOUT PRINTS DATA FROM THE CALCULATION.

subroutine diout(E_total,E_v,E_Vxc,E_xc,E_Z,ibav,icheck, idirc,n_orb,n_ray,nnlm,psi,psi_small,Q,ray,rh, &
                   rho,Summa,Vrf,xe,xj,xl,xn,xte,xz,Z)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter(abz = 1.501_db, ahz = 51.001_db)

  character(len=6), dimension(80):: nlj

  integer:: Z
  integer, dimension(n_orb):: nist
  integer, dimension(n_orb):: xl, xn

  real(kind=db), dimension(n_orb):: xe, xj, xz
  real(kind=db), dimension(n_ray):: ray, rh, rho, Vrf
  real(kind=db), dimension(n_ray,nnlm):: psi, psi_small

  data nlj / ' 1S1/2', ' 2S1/2', ' 2P1/2', ' 2P3/2', ' 3S1/2', ' 3P1/2', ' 3P3/2', ' 3D3/2', ' 3D5/2', ' 4S1/2', ' 4P1/2', &
   ' 4P3/2', ' 4D3/2', ' 4D5/2', ' 4F5/2', ' 4F7/2', ' 5S1/2', ' 5P1/2', ' 5P3/2', ' 5D3/2', ' 5D5/2', ' 5F5/2', ' 5F7/2', &
   ' 5G7/2', ' 5G9/2', ' 6S1/2', ' 6P1/2', ' 6P3/2' , ' 6D3/2', ' 6D5/2', ' 6F5/2', ' 6F7/2', ' 6G7/2', &
   ' 6G9/2', ' 6H9/2', '6H11/2', ' 7S1/2', ' 7P1/2', ' 7P3/2', ' 7D3/2', ' 7D5/2', ' 7F5/2', ' 7F7/2', ' 7G7/2', ' 7G9/2', &
   ' 7H9/2', '7H11/2', '7I11/2', '7I13/2', ' 8S1/2', '1S', '2S', '2P', '3S', '3P', '3D', '4S', '4P', '4D', &
   '4F', '5S', '5P', '5D', '5F', '5G', '6S', '6P', '6D', '6F', '6G', &
   '6H', '7S', '7P', '7D', '7F', '7G', '7H', '7I', '8S', '8P' /

  p4 = 1 / ( 4 * pi )

  do i = 1, n_orb
    if( idirc == 0 )  then
      nist(i) = nint( 0.5_db * xn(i) * (xn(i) - 1) + xl(i) + ahz )
    else
      nist(i) = nint( xn(i) * (xn(i) - 2) + xl(i) + xj(i) + abz )
    endif
  end do

! E_v = Integrale( (Z/r + V_hartree)*rho*d3r )
! E_Z = Integrale( (Z/r)*rho*d3r )
  Ekin = Summa - E_v - E_Vxc - xte
  U = 0.5_db * ( E_v - E_Z )

  Hartree = 2 * Rydb

  if( icheck > 1 ) then
    write(ibav,12)
    do i = 1, n_orb
      ni = nist(i)
      write(ibav,14) nlj(ni), xn(i), xl(i), xj(i), xz(i), Hartree * xe(i)
    end do
    write(ibav,15) Z, Q, Hartree * E_total, Hartree * (E_Z + U), Hartree * E_xc, Hartree * Ekin, &
                   Hartree * Summa, Hartree * U, Hartree * E_Vxc
  endif

  rho(1:n_ray) = p4 * rh(1:n_ray) / ray(1:n_ray)**2

  if( icheck > 1 ) then
    write(ibav,110) nlj( nist(1:n_orb) )
    do i = 1,n_ray
      f = quatre_pi * ray(i)**2
      write(ibav,120) bohr*ray(i), Hartree*Vrf(i)/ray(i), f*rho(i), psi(i,1:n_orb)
    end do
    if( idirc/= 0 ) then
      write(ibav,130) nlj( nist(1:n_orb) )
      do i = 1,n_ray
        f = quatre_pi * ray(i)**2
        write(ibav,120) bohr*ray(i), Hartree*vrf(i)/ray(i), f*rho(i), psi_small(i,1:n_orb)
      end do
    endif
  endif

   12 format(/' Orbital  n  l  j    Electron  Eigenvalue(eV)')
   14 format(a6,i5,i3,f4.1,f10.3,3f15.3)
   15 format(/' Nuclear charge                  =',i12/, ' Integral of charge density      =',f12.3/, &
          ' Total energy = Ue + U + Exc + T =',f12.3/, ' Potential energy = U + Ue       =',f12.3/, &
          ' Exchange energy = Exc           =',f12.3/, ' Kinetic energy = T              =',f12.3/, &
          ' Sum of the energy eigenvalues   =',f12.3/, ' 0.5*Integral(Vh.rho.dr)         =',f12.3/, &
          ' Integral(Vxc.rho.dr)            =',f12.3)
  110 format(/'      Radius    Potential  4*pi*r2*Rho ',36(4x,a6,3x))
  120 format(f13.6,1p,38e13.5)
  130 format(/'  Small component:',// '      Radius    Potential  4*pi*r2*Rho ',36(4x,a6,3x))
end

