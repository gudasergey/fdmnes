! FDMNES subroutines

! Calculation of the minimum of metrics in a grid of points.

subroutine minim(fdmfit_out,index_Met_Fit,itape_minim,Minim_fdm_ok,minimok,ncali,ndm,ngroup_par,ngroup_par_conv,nmetric, &
               nmetricm,Nom_Met,nparam,npm,par,par_op,typepar)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  character(len=132) :: fdmfit_out
  character(len=2), dimension(nmetricm) :: Nom_Met
  character(len=9), dimension(ngroup_par) :: typepar

  integer, dimension(ngroup_par) :: ic, nparam
  integer, dimension(ngroup_par,nmetric) :: ia

  logical ellipse, Minim_fdm_ok, minimok, polynome, polynomefdm

  real(kind=db), dimension(nmetric) :: distmin_ell, distmin_pol
  real(kind=db), dimension(ngroup_par) :: a, b, parm_ext, par_op, parp_ext, sol, sm
  real(kind=db), dimension(ngroup_par,nmetric) :: distm, par_pol, par_ell
  real(kind=db), dimension(ngroup_par,ngroup_par) :: ab, mat
  real(kind=db), dimension(ngroup_par,npm) :: par
  real(kind=db), dimension(ngroup_par,npm,nmetric) :: dist_cut

  itape = itape_minim
  ellipse = .true.
  polynome = .false.
  polynomefdm = .false.

  ncal = ncali * ndm

  ngr_par = ngroup_par
  do igr = 1,ngroup_par
    if( nparam(igr) < 3 ) ngr_par = ngr_par - 1
  end do

! Valeurs extremes possibles des parametres
  do igr = 1,ngroup_par
    if( nparam(igr) < 2 ) then
      parm_ext(igr) = par(igr,1)
      parp_ext(igr) = par(igr,1)
    else
      pp1 = 2 * par(igr,1) - par(igr,2)
      pp2 = 2 * par(igr,nparam(igr)) - par(igr,nparam(igr)-1)
      if( par(igr,1) < par(igr,2) ) then
        parm_ext(igr) = pp1
        parp_ext(igr) = pp2
      else
        parm_ext(igr) = pp2
        parp_ext(igr) = pp1
      endif
    endif
  end do

! Boucle sur les distances metriques.
  boucle_im: do im = 1,nmetric

! Recherche du minima des distances metriques dans la grille d'entrees.
    dmin = 10000000.
    Rewind(itape)
    do ical = 1,ncal
      read(itape,*) ( dist, i = 1,im )
      if( dist > dmin ) cycle
      dmin = dist
      ical_min = ical
    end do
    index_read = ncal

    ndem = 1
    do igr = 1,ngroup_par
      j = ( ical_min - 1 ) / ndem
      ndem = ndem * nparam(igr)
      ia(igr,im) = mod( j, nparam(igr) ) + 1
    end do

! Recherche des points pour les coupes autour des minimums
!    if( ncali > 1 ) then
      do igr = 1,ngroup_par
        if( nparam(igr) < 2 ) cycle
        do ip = 1,nparam(igr)
          ic(:) = ia(:,im)
          ic(igr) = ip
          iligne = icalical(ic,nparam,ngroup_par)
          dist_cut(igr,ip,im) =read_dist(iligne,index_read,itape,im)
          index_read = iligne
        end do
      end do
!    endif

! Polynome du second degre.

    do igr = 1,ngroup_par

      if( nparam(igr) < 3 ) then
        par_pol(igr,im) = par(igr,ia(igr,im))
        cycle
      endif

      ic(:) = ia(:,im)
      if( ia(igr,im) == 1 ) then
        ic(igr) = ia(igr,im) + 1
      elseif( ia(igr,im) == nparam(igr) ) then
        ic(igr) = ia(igr,im) - 2
      else
        ic(igr) = ia(igr,im) - 1
      endif
      ical = icalical(ic,nparam,ngroup_par)
      dist = read_dist(ical,index_read,itape,im)
      index_read = ical

      dm = dist - dmin
      xm = par(igr,ic(igr)) - par(igr,ia(igr,im))
      if( ia(igr,im) == 1 ) then
        ic(igr) = ia(igr,im) + 2
      elseif( ia(igr,im) == nparam(igr) ) then
        ic(igr) = ia(igr,im) - 1
      else
        ic(igr) = ia(igr,im) + 1
      endif
      ical = icalical(ic,nparam,ngroup_par)
      dist = read_dist(ical,index_read,itape,im)
      index_read = ical
      dp = dist - dmin
      xp = par(igr,ic(igr)) - par(igr,ia(igr,im))
      a(igr) = ( dm / xm - dp / xp ) / ( xm - xp )
      b(igr) = dm / xm - a(igr) * xm

      if( a(igr) <= 0._db ) then
! Cas de la courbe convexe.
        distm(igr,im) = dmin
        par_pol(igr,im) = par(igr,ia(igr,im))
      else
        par_pol(igr,im) = -0.5 * b(igr) / a(igr)

        p = par_pol(igr,im) + par(igr,ia(igr,im))
        if( p < parm_ext(igr) ) then
          par_pol(igr,im) = parm_ext(igr) - par(igr,ia(igr,im))
        elseif( p > parp_ext(igr) ) then
          par_pol(igr,im) = parp_ext(igr) - par(igr,ia(igr,im))
        endif
        distm(igr,im) = a(igr) * par_pol(igr,im)**2 + b(igr) * par_pol(igr,im) + dmin
        if( distm(igr,im) < 0.5 * dmin ) then
          distm(igr,im) = dmin
          par_pol(igr,im) = par(igr,ia(igr,im))
        else
          if( igr > 1 .and. im == 1 ) polynome = .true.
          if( igr > ngroup_par_conv .and. im == 1 ) polynomefdm = .true.
          par_pol(igr,im) = par_pol(igr,im) + par(igr,ia(igr,im))
        endif

      endif

    end do

    if( ngroup_par < 2 ) cycle

! Calcul de l'hyperellipsoide.

    ab(:,:) = 0._db
    do igr = 1,ngroup_par
      if( nparam(igr) < 3 ) cycle
      do jgr = 1,igr - 1
        if( nparam(jgr) < 3 ) cycle
        ifac = 0
        do i = -1,1,2
          if( ( ia(igr,im) == 1 .and. i == - 1 ) .or. ( ia(igr,im) == nparam(igr) .and. i == 1 ) ) cycle
          do j = -1,1,2
            if( ( ia(jgr,im) == 1 .and. j == - 1 ) .or. ( ia(jgr,im) == nparam(jgr) .and. j == 1 ) ) cycle
            ifac = ifac + 1
            ic(:) = ia(:,im)
            ic(igr) = ia(igr,im) + i
            ic(jgr) = ia(jgr,im) + j
            ical = icalical(ic,nparam,ngroup_par)
            dist = read_dist(ical,index_read,itape,im)
            index_read = ical
            dd = dist - dmin
            xi = par(igr,ic(igr)) - par(igr,ia(igr,im))
            xj = par(jgr,ic(jgr)) - par(jgr,ia(jgr,im))
            ab(igr,jgr) = ab(igr,jgr) + ( dd - a(igr)*xi**2 - a(jgr)*xj**2 - b(igr)*xi - b(jgr)*xj ) / ( xi * xj )
          end do
        end do
        ab(igr,jgr) = 0.5 * ab(igr,jgr) / ifac
        ab(jgr,igr) = ab(igr,jgr)
      end do
    end do

! Calcul des parametres optimises.

! Remplissage de la matrice

    i = 0
    do igr = 1,ngroup_par
      if( nparam(igr) < 3 ) cycle
      i = i + 1
      j = 0
      do jgr = 1,ngroup_par
        if( nparam(jgr) < 3 ) cycle
        j = j + 1
        if( igr /= jgr ) then
          mat(i,j) = ab(igr,jgr)
        else
          mat(i,j) = 2 * a(igr)
        endif
      end do
      sm(i) = - b(igr)
    end do

! Triangularisation
    l = ngr_par
    if( abs( mat(l,l) ) < 0.0000000001_db ) then
      if( im == 1 ) ellipse = .false.
      distmin_ell(im) = dmin
      do igr = 1,ngroup_par
        par_ell(igr,im) = par(igr,ia(igr,im))
      end do
      cycle boucle_im
    endif

    fac = 1 / mat(l,l)
    do k = 1,l
      mat(l,k) = fac * mat(l,k)
    end do
    sm(l) = fac * sm(l)

    do l = ngr_par-1,1,-1
      do j = ngr_par,l+1,-1
        do k = j,1,-1
          mat(l,k) = mat(l,k) - mat(l,j) * mat(j,k)
        end do
        sm(l) = sm(l) - mat(l,j) * sm(j)
      end do

      if( abs( mat(l,l) ) < 0.0000000001_db ) then
        if( im == 1 ) ellipse = .false.
        distmin_ell(im) = dmin
        do igr = 1,ngroup_par
          par_ell(igr,im) = par(igr,ia(igr,im))
        end do
        cycle boucle_im
      endif

      fac = 1 / mat(l,l)
      do k = 1,l
        mat(l,k) = fac * mat(l,k)
      end do
      sm(l) = fac * sm(l)

    end do

! Resolution
    do i = 1,ngr_par
      sol(i) = sm(i)
      do j = 1,i-1
        sol(i) = sol(i) - mat(i,j) * sol(j)
      end do
    end do

! Test pour voir si on est trop loin de la grille
    i = 0
    do igr = 1,ngroup_par
      if( nparam(igr) < 3 ) cycle
      i = i + 1
      p = sol(i) + par(igr,ia(igr,im))
      if( p < parm_ext(igr) ) then
        sol(i) = parm_ext(igr) - par(igr,ia(igr,im))
      elseif( p > parp_ext(igr) ) then
        sol(i) = parp_ext(igr) - par(igr,ia(igr,im))
      endif
    end do

! Metrique minimum
    distmin_ell(im) = dmin
    i = 0
    do igr = 1,ngroup_par
      if( nparam(igr) < 3 ) cycle
      i = i + 1
      j = 0
      distmin_ell(im) = distmin_ell(im) + a(igr) * sol(i)**2 + b(igr)*sol(i)
      do jgr = 1,ngroup_par
        if( nparam(jgr) < 3 .or. igr == jgr ) cycle
        j = j + 1
        distmin_ell(im) = distmin_ell(im) + ab(igr,jgr) * sol(i) * sol(j)
      end do
    end do

    i = 0
    do igr = 1,ngroup_par
      if( nparam(igr) < 3 ) then
        par_ell(igr,im) = par(igr,ia(igr,im))
      else
        i = i + 1
        par_ell(igr,im) = sol(i) + par(igr,ia(igr,im))
      endif
    end do

    if( distmin_ell(im) > dmin + 0.00000001_db .or. distmin_ell(im) < 0.5 * dmin ) then
      if( im == 1 ) ellipse = .false.
      distmin_ell(im) = dmin
      do igr = 1,ngroup_par
        par_ell(igr,im) = par(igr,ia(igr,im))
      end do
    endif

  end do boucle_im

  do im = 1,nmetric
    distmin_pol(im) = 0._db
    do igr = 1,ngroup_par
      if( nparam(igr) < 3 ) cycle
      distmin_pol(im) = max( distmin_pol(im), distm(igr,im) )
    end do
  end do

  if( ellipse ) then
    par_op(2:ngroup_par) = par_ell(2:ngroup_par,index_Met_Fit)
  else
    par_op(2:ngroup_par) = par_pol(2:ngroup_par,index_Met_Fit)
  endif
  minimok = ellipse .or. polynome
  Minim_fdm_ok = ellipse .or. polynomefdm

  open( 4, file = fdmfit_out, position='append')

  do ipr = 4,6,2

    write(ipr,110) Nom_Met(1:nmetric)
    write(ipr,115) distmin_pol(1:nmetric)
    do igr = 1,ngroup_par
      write(ipr,120) adjustr(typepar(igr)), par_pol(igr,:)
    end do

    if( ngroup_par > 1 ) then
      write(ipr,130) Nom_Met(1:nmetric)
      write(ipr,135) distmin_ell(1:nmetric)
      do igr = 1,ngroup_par
        write(ipr,120) adjustr(typepar(igr)), par_ell(igr,:)
      end do
    endif

  end do

! Ecriture des coupes autour des minima de la grille.

!  if( ncali > 1 ) then
    do igr = 1,ngroup_par
      if( nparam(igr) < 2 ) cycle
      write(4,140) adjustl(typepar(igr)), adjustr(typepar(igr)), Nom_Met(1:nmetric)
      do ip = 1,nparam(igr)
        write(4,150) par(igr,ip), ( dist_cut(igr,ip,im), im = 1,nmetric )
      end do
      if( nmetric == 3 ) then
        write(4,160) distm(igr,1:nmetric)
        write(4,165) par_pol(igr,1:nmetric)
      else
        write(4,166) distm(igr,1:nmetric)
        write(4,167) par_pol(igr,1:nmetric)
      endif
      if( ngroup_par > 1 ) write(4,'(A)') '      with the other parameters :'
      do jgr = 1,ngroup_par
        if( jgr == igr ) cycle
        write(4,170) adjustr(typepar(jgr)), par(jgr,ia(jgr,:))
      end do
    end do
!  endif

  Close(4)

  return
  110 format(//'  Minimum for the second order polynomia :', /18x,4(3x,a2,'min',3x))
  115 format(16x,2f11.5,2f11.6)
  120 format(5x,a9,' =',4f11.5)
  130 format(/'  Minimum for the hyperellipsoid :', /18x,4(3x,a2,'min',3x))
  135 format(16x,2f11.5,2f11.6)
  140 format(/'  Cut along the parameter ',a9/,7x,a9,4(7x,a2,2x))
  150 format(1x,f15.5,2f11.5,2f11.6)
  160 format(/16x,2f11.5,f11.6,' : minimum')
  165 format(16x,3f11.5,' : corresponding parameter value')
  166 format(/16x,2f11.5,2f11.6,' : minimum')
  167 format(16x,4f11.5,' : corresponding parameter value')
  170 format(5x,a9,' =',4f11.5)
end

!***********************************************************************

function icalical(ib,nparam,ngroup_par)

  integer, dimension(ngroup_par) :: ib, nparam

  ical = ib(1)
  do jgr = 2,ngroup_par
    nn = nparam(1)
    do kgr = 2,jgr-1
      nn = nn * nparam(kgr)
    end do
    ical = ical + nn * ( ib(jgr) - 1 )
  end do

  icalical = ical

  return
end

!***********************************************************************

function read_dist(iligne,index_read,itape,icol)

  use declarations
  implicit none

  integer:: i, icol, iligne, index_read, itape, ligne, ligne_1

  real(kind=db):: dist, read_dist

  if( iligne > index_read ) then
    ligne_1 = index_read + 1
  else
    ligne_1 = 1
    Rewind(itape)
  endif

  do ligne = ligne_1,iligne-1
    read(itape,*)
  end do

  read(itape,*) ( dist, i = 1,icol )

  read_dist = dist

  return
end



