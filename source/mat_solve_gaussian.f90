
! Routine solving the system of linear equations by Gaussian elimination

subroutine mat_solve(Base_hexa, Basereel, Bessel, Besselr, Cal_comp, cgrad, clapl, E_comp, Eimag, Enervide, gradvr, &
        ianew, iato, ibord, icheck, ie, igrph, ii, isbord, iso, ispinin, isrt, isvois, ivois, Kar, Kari, lato, &
        lb1, lb2, lmaxso, lso, mato, MPI_host_num_for_mumps, mpirank0, mso, natome, nbm, nbord, nbordf, nbtm, Neuman, Neumanr, &
        new, newinv, ngrph, nicm, nim, nligne, nligne_i, nligneso, nlmsam,  nlmagm, nlmmax, nlmomax, nlmsa, nlmso, nlmso_i, &
        nphiato1, nphiato7, npoint, npsom, nsm, nso1, nsort, nsort_c, nsort_r, nsortf, nspin, nspino, nspinp, nspinr, nstm, &
        numia, nvois, phiato, poidsa, poidso, Relativiste, Repres_comp, rvol, smi, smr, Spinorbite, Time_fill, Time_tria, Vr, &
        Ylm_comp, Ylmato, Ylmso)
  
  use declarations
  implicit none

  integer:: i, i_newind, ia, ib, icheck, ie, igrph, ii, ipr, isp, ispin, ispinin, iv, j, jj, k, lb1i, lb1r, lb2i, lb2r, lm, &
    lmaxso, lms, MPI_host_num_for_mumps, mpirank0, natome, nbm, nbtm, ngrph, nicm, nim, nligne, nligne_i, nligneso, nlmagm, &
    nlmmax, nlmomax, nlmsam, nlmso, nlmso_i, nphiato1, nphiato7, npoint, & 
    npsom, nsm, nso1, nsort, nsort_c, nsort_r, nsortf, nspin, nspino, nspinp, nspinr, nstm, nvois

  integer(kind=db):: lk, nb_not_zero_tot
  integer(kind=db), dimension(nligne):: lbz
!  integer:: lk, nb_not_zero_tot
!  integer, dimension(nligne):: lbz
   
  integer, dimension(0:npoint):: new
  integer, dimension(natome):: ianew, nbord, nbordf, nlmsa
  integer, dimension(nstm):: isrt
  integer, dimension(npsom):: numia
  integer, dimension(nso1,ngrph):: lso, mso, iso
  integer, dimension(nlmsam,natome,ngrph):: lato, mato, iato
  integer, dimension(nbtm,natome):: ibord, isbord
  integer, dimension(npsom,nvois):: ivois, isvois
  integer, dimension(nligne):: lb1, lb2, nb_not_zero, newinv

  integer, dimension(:), allocatable:: IndColNotZero, IndColNotZero_inv

  character(len=2), dimension(nligne):: mletl

  complex(kind=db), dimension(nsort_c,0:lmaxso,nspinr):: Bessel, Neuman

  logical:: Base_hexa, Basereel, Cal_comp, E_comp, Relativiste, Repres_comp, Second, Spinorbite, Ylm_comp
  logical, dimension(nligne):: ColNotZero 

  real(kind=sg):: time

  real(kind=db):: den, Enervide, Eimag, fac, faci, facr, fi, fr, tp1, tp2, tp3, Time_fill, Time_tria, vnormali, vnormalr
  
  real(kind=db), dimension(nopsm,nspino):: Kar, Kari
  real(kind=db), dimension(nvois):: cgrad
  real(kind=db), dimension(0:nvois):: clapl
  real(kind=db), dimension(npoint,nspin):: Vr
  real(kind=db), dimension(nbm,natome):: poidsa
  real(kind=db), dimension(nim):: rvol
  real(kind=db), dimension(nsm):: poidso
  real(kind=db), dimension(nsort,nlmomax):: Ylmso
  real(kind=db), dimension(nicm,3,nspino):: gradvr
  real(kind=db), dimension(nbtm,nlmmax,natome):: Ylmato
  real(kind=db), dimension(nphiato1,nlmagm,nspinp,nspino,natome,nphiato7):: phiato
  real(kind=db), dimension(nlmso,nligne):: smr
  real(kind=db), dimension(nlmso_i,nligne_i):: smi
  real(kind=db), dimension(nsort_r,0:lmaxso,nspinr):: Besselr, Neumanr

  real(kind=db), dimension(:), allocatable:: abi, abr, abvi, abvr

! Evaluation indicag des colonnes non nulles apres triangularisation
  do i_newind = 1,2
! 1 = Dimensionnement  

    if( i_newind == 1 ) then
      Second = .false.
    else
      Second = .true.
      allocate( IndColNotZero_inv(nb_not_zero_tot) ) 
    endif
      
    nb_not_zero_tot = 0
  
    ColNotZero(:) = .false.
    nb_not_zero(:) = 0
 
    do ii = nligne,1,-1
 
      i = newinv(ii)

! Points du reseau de base
      if( i > 0 ) then

        do iv = 1,nvois

          j = ivois(i,iv)
          ia = numia(j)

! Le voisin est dans le reseau de base :
          if( ia == 0 ) then
            ColNotZero(new(j)-nspino+1:new(j)) = .true.
          elseif( ia > 0 ) then
            ColNotZero(ianew(ia)+1:ianew(ia)+nlmsa(ia)) = .true.
          elseif( ia == -2 ) then
            ColNotZero(nligneso+1:nligneso+nlmso) = .true.
          endif
        
        end do

! Developpement en sortie
      elseif( i == 0 ) then

        ColNotZero(nligneso+1:nligneso+nlmso) = .true.

        lm = ii - nligneso
        isp = iso(lm,igrph)
 
        do ib = 1,nsortf
          j = isrt(ib)
          jj = new(j) - nspino + isp
          ColNotZero(jj) = .true.
        end do
    
! Developpement dans un atome
      else

        ia = -i

        ColNotZero(ianew(ia)+1:ianew(ia)+nlmsa(ia)) = .true.

        lm = ii - ianew(ia)
        isp = iato(lm,ia,igrph)

        do ib = 1,nbordf(ia)
          j = ibord(ib,ia)
          jj = new(j) - nspino + isp
          ColNotZero(jj) = .true.
        end do

      endif

      do j = ii-1,1,-1
        if( .not. ColNotZero(j) ) cycle
        nb_not_zero_tot = nb_not_zero_tot + 1
        if( Second ) then 
          nb_not_zero(ii) = nb_not_zero(ii) + 1
          IndColNotZero_inv(nb_not_zero_tot) = j
        endif  
      end do
    
    end do

  end do ! fin boucle i_newind

  allocate( IndColNotZero(nb_not_zero_tot) ) 
  
  do lk = 1,nb_not_zero_tot
    IndColNotZero( lk ) = IndColNotZero_inv( nb_not_zero_tot + 1 - lk )
  end do

  deallocate( IndColNotZero_inv )

  lbz(1) = 0
  do ii = 2,nligne
    lbz(ii) = lbz(ii-1) + nb_not_zero(ii-1)
  end do 
      
  if( icheck > 2 ) then
    write(3,110) nb_not_zero_tot
    do ii = 1,nligne
      write(3,120) ii, nb_not_zero(ii), lbz(ii)
      write(3,130) IndColNotZero(lbz(ii)+1:lbz(ii)+nb_not_zero(ii))
    end do
  endif
    
! Debut remplissage

  allocate( abr(nb_not_zero_tot) )
  abr(:) = 0._db

  if( Cal_comp ) then
    allocate( abi(nb_not_zero_tot) )
    abi(:) = 0._db
  endif  

  if( Spinorbite ) then
    ispin = 1
  else
    ispin = ispinin
  endif

  boucle_ii: do ii = nligne,1,-1

    call CPU_TIME(time)
    tp1 = real(time,db)

    lb1r = lb1(ii)
    lb2r = lb2(ii)
    if( Cal_comp ) then
      lb1i = lb1(ii)
      lb2i = lb2(ii)
    else
      lb1i = 0
      lb2i = 0
    endif
    allocate( abvr(lb1r:lb2r) )
    abvr(:) = 0._db
    allocate( abvi(lb1i:lb2i) )
    abvi(:) = 0._db

    call calcMatRow( abvr, abvi, Base_hexa, Basereel, Bessel, Besselr, Cal_comp, cgrad, clapl, E_comp, Eimag, &
    Enervide, gradvr, ianew, iato, ibord, icheck, igrph, ii, isbord, iso, ispin, isrt, isvois, ivois, Kar, Kari, &
    lato, lb1i, lb1r, lb2i, lb2r, lmaxso, lso, mato, mletl, mso, natome, nbm, nbord, &
    nbordf, nbtm, Neuman, Neumanr, new, newinv, ngrph, nicm, nim, nligne, nligne_i, &
    nligneso, nlmagm, nlmmax, nlmomax, nlmsa, nlmsam, nlmso, nlmso_i, nphiato1, nphiato7, npoint, npsom, nsm, nso1, &
    nsort, nsort_c, nsort_r, nsortf, nspin, nspino, nspinp, nspinr, nstm, &
    numia, nvois, phiato, poidsa, poidso, Relativiste, Repres_comp, rvol, Spinorbite,  &
    smi, smr, Vr, Ylm_comp, Ylmato, Ylmso )

    call CPU_TIME(time)
    tp2 = real(time,db)

    Time_fill = Time_fill + tp2 - tp1

! Triangularisation :

    do j = lb2(ii),ii+1,-1
      if( Cal_comp ) then

        facr = abvr(j)
        faci = abvi(j)
        if( abs( facr ) > eps10 .or. abs( faci ) > eps10 ) then
!CDIR NODEP
          do lk =  lbz(j)+1, lbz(j)+nb_not_zero(j)
            k = IndColNotZero(lk)
            abvr(k) = abvr(k) - facr * abr(lk) + faci * abi(lk)
            abvi(k) = abvi(k) - facr * abi(lk) - faci * abr(lk)
          end do

          smr(1:nlmso,ii) = smr(1:nlmso,ii) - facr * smr(1:nlmso,j) + faci * smi(1:nlmso,j)
          smi(1:nlmso,ii) = smi(1:nlmso,ii) - facr * smi(1:nlmso,j) - faci * smr(1:nlmso,j)
        endif

      else

        fac = abvr(j)
        if( abs( fac ) > eps10 ) then
!CDIR NODEP
          do lk =  lbz(j)+1, lbz(j)+nb_not_zero(j)
            k = IndColNotZero(lk)
            abvr(k) = abvr(k) - fac * abr(lk)
          end do
          smr(1:nlmso,ii) = smr(1:nlmso,ii) - fac * smr(1:nlmso,j)
        endif

      endif
    end do

! Normalisation :
    if( Cal_comp ) then
      den = abvr(ii)**2 + abvi(ii)**2
      if( den < 1.e-20_db .and. mpirank0 == 0 ) then
        call write_error
        do ipr = 3,9,3
          write(ipr,140) ii, i, ispin, ia, lm
        end do
        stop
      endif
      vnormalr = abvr(ii) / den
      vnormali = - abvi(ii) / den
      do lk =  lbz(ii)+1, lbz(ii)+nb_not_zero(ii)
        k = IndColNotZero(lk)
        abr(lk) = abvr(k) * vnormalr - abvi(k) * vnormali
        abi(lk) = abvi(k) * vnormalr + abvr(k) * vnormali
      end do
!CDIR NODEP
      do lms = 1,nlmso
        fr = smr(lms,ii) * vnormalr - smi(lms,ii) * vnormali
        fi = smr(lms,ii) * vnormali + smi(lms,ii) * vnormalr
        smr(lms,ii) = fr
        smi(lms,ii) = fi
      end do
    else
      if( abs( abvr(ii) ) < 1.e-10_db .and. mpirank0 == 0 ) then
        call write_error
        do ipr = 3,9,3
          write(ipr,140) ii, i, ispin, ia, lm
        end do
        stop
      endif
      vnormalr = 1 / abvr(ii)
      do lk =  lbz(ii)+1, lbz(ii)+nb_not_zero(ii)
        k = IndColNotZero(lk)
        abr(lk) =  abvr(k) * vnormalr
      end do
      smr(1:nlmso,ii) = smr(1:nlmso,ii) * vnormalr
    endif

    deallocate( abvr )
    deallocate( abvi )

    call CPU_TIME(time)
    tp3 = real(time,db)

    Time_tria = Time_tria + tp3 - tp2

  end do boucle_ii      ! fin de la boucle sur les lignes

  if( icheck > 2 ) then
    write(3,'(/A/)') '   IndcolNotZero, abr, after triangularisation'
    do ii = 2,nligne
      Write(3,'(a5,i6)') ' ii =', ii
      write(3,'(10(i6,1p,e13.5))') ( IndcolNotZero(lk), abr(lk), lk = lbz(ii)+1, lbz(ii)+nb_not_zero(ii) )
    end do

    write(3,'(/A)') '   ii   sm(lmf = 1,nlmso) after triangularisation'
    do ii = 1,nligne
      if( Cal_comp ) then
        write(3,150) ii, ( smr(lms,ii), smi(lms,ii), lms = 1,nlmso )
      else
        write(3,160) ii, smr(1:nlmso,ii)
      endif
    end do
  endif

! Resolution :

  do ii = 2,nligne
    do lk =  lbz(ii)+1, lbz(ii)+nb_not_zero(ii)
      j = IndColNotZero(lk)
      if( Cal_comp ) then
        smr(1:nlmso,ii) = smr(1:nlmso,ii) - abr(lk)*smr(1:nlmso,j) + abi(lk) * smi(1:nlmso,j)
        smi(1:nlmso,ii) = smi(1:nlmso,ii) - abr(lk)*smi(1:nlmso,j) - abi(lk) * smr(1:nlmso,j)
      else
        smr(1:nlmso,ii) = smr(1:nlmso,ii) - abr(lk)*smr(1:nlmso,j)
      endif
    end do
  end do

  deallocate( abr )
  if( Cal_comp ) deallocate( abi )
  deallocate( IndColNotZero )

  return
  110 format(/' nb_not_zero_tot =',i9)
  120 format(/' Line =',i7,', nb_not_zero =',i6,', lbz =',i9,' / IndColNotZero')
  130 format(20i6)
  140  format(/' Division par zero dans mat',// &
               '  ii =',i6,', i =',i6,', ispin =',i3,', ia =',i3,', lm =',i3)
  150 format(i6,1p,250(1x,2e11.3))
  160 format(i6,1p,500e11.3)
end
      
!***********************************************************************
        
subroutine getSolverParams(MPI_host_num_for_mumps,mpinodes0,Solver)

  use declarations
  implicit none

  integer:: MPI_host_num_for_mumps, mpinodes0
  
  Character(len=5):: Solver
  
  Solver = 'Gauss'
  MPI_host_num_for_mumps = 1

  return      
end
      