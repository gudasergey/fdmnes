
! Routine solving the system of linear equations using sparse matrix technique and MUMPS
! From Alexander Guda et al, Rostov, Russia. 

subroutine mat_solve(Base_hexa, Basereel, Bessel, Besselr, Cal_comp, cgrad, clapl, E_comp, Eimag, Enervide, gradvr, &
        ianew, iato, ibord, icheck, igrph, ii, isbord, iso, ispinin, isrt, isvois, ivois, Kar, Kari, lato, &
        lb1, lb2, lmaxso, lso, mato, MPI_host_num_for_mumps, mpirank0, mso, natome, nbm, nbord, nbordf, nbtm, Neuman, Neumanr, &
        new, newinv, ngrph, nicm, nim, nligne, nligne_i, nligneso, nlmsam,  nlmagm, nlmmax, nlmomax, nlmsa, nlmso, nlmso_i, &
        nphiato1, nphiato7, npoint, npsom, nsm, nso1, nsort, nsort_c, nsort_r, nsortf, nspin, nspino, nspinp, nspinr, nstm, &
        numia, nvois, phiato, poidsa, poidso, Relativiste, Repres_comp, rvol, smi, smr, Spinorbite, tpt1, tpt2, Vr, Ylmato, Ylmso)

  use declarations
  implicit none

  integer:: icheck, igrph, ii, isp, ispin, ispinin, j, lb1i, lb1r, lb2i, lb2r, &
    lmaxso, MPI_host_num_for_mumps, mpirank_in_mumps_group, mpirank0, natome, nbm, nbtm, ngrph, nicm, nim, nligne, &
    nligne_i, nligneso, nlmagm, nlmmax, nlmomax, nlmsam, nlmso, nlmso_i, nphiato1, nphiato7, npoint, & 
    npsom, nsm, nso1, nsort, nsort_c, nsort_r, nsortf, nspin, nspino, nspinp, nspinr, nstm, nvois
  
  integer, dimension(0:npoint):: new
  integer, dimension(natome):: ianew, nbord, nbordf, nlmsa
  integer, dimension(nstm):: isrt
  integer, dimension(npsom):: numia
  integer, dimension(nso1,ngrph):: lso, mso, iso
  integer, dimension(nlmsam,natome,ngrph):: lato, mato, iato
  integer, dimension(nbtm,natome):: ibord, isbord
  integer, dimension(npsom,nvois):: ivois, isvois
  integer, dimension(nligne):: lb1, lb2, newinv

  integer(kind=8):: inz, nz, nz_i, nligne8, newSize, oldSize
  integer(kind=8), dimension(:), allocatable:: rowIndexes, columnIndexes, tempRI, tempCI

  character(len=2), dimension(nligne):: mletl

  complex(kind=db), dimension(nsort_c,0:lmaxso,nspinr):: Bessel, Neuman

  logical:: Base_hexa, Basereel, Cal_comp, E_comp, Relativiste, Repres_comp, Spinorbite

  real(kind=sg):: time

  real(kind=db):: Enervide, Eimag, p2, tp1, tp2, tp3, tpt1, tpt2
  
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

  real(kind=db), dimension(:), allocatable:: A, A_im, tempA, tempAi, abvi, abvr

  mpirank_in_mumps_group = mod( mpirank0, MPI_host_num_for_mumps )

!  Define problem on the host (processor 0)      
  if( mpirank_in_mumps_group == 0 ) then
  
    call CPU_TIME(time)
    tp1 = real(time,db)

    if( Spinorbite ) then
      ispin = 1
    else
      ispin = ispinin
    endif

    nz = nligne/100*nligne
    allocate( rowIndexes ( nz ) )
    allocate( columnIndexes ( nz ) )
    allocate( A( nz ) )
    if ( Cal_comp ) then
      nz_i = nz
    else
      nz_i = 0
    endif
    allocate( A_im( nz_i ) )
    
    smr(:,:) = 0._db

    if( Cal_comp ) smi(:,:) = 0._db  
    
    inz = 0
    do ii = nligne,1,-1
      
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
      smi, smr, Vr, Ylmato, Ylmso )

! expand arrays            
      oldSize = size(A)
      if ( oldSize < inz+lb2r-lb1r+1 ) then
        allocate(tempRI(oldSize))
        allocate(tempCI(oldSize))
        allocate(tempA(oldSize))
        if( Cal_comp ) allocate(tempAi(oldSize))
        tempRI = rowIndexes
        tempCI = columnIndexes
        tempA = A
        if( Cal_comp ) tempAi = A_im
        deallocate(rowIndexes)
        deallocate(columnIndexes)
        deallocate(A)
        if( Cal_comp ) deallocate(A_im)
        newSize = inz+lb2r-lb1r+1
        if (newSize < oldSize*2) newSize = oldSize*2
        allocate(rowIndexes(newSize))
        allocate(columnIndexes(newSize))
        allocate(A(newSize))
        if( Cal_comp ) allocate(A_im(newSize))
        rowIndexes(1:oldSize) = tempRI
        columnIndexes(1:oldSize) = tempCI
        A(1:oldSize) = tempA
        if ( Cal_comp ) A_im(1:oldSize) = tempAi
        deallocate(tempRI)
        deallocate(tempCI)
        deallocate(tempA)
        if ( Cal_comp ) deallocate(tempAi)
      endif
      
! Fill matrises
      do j = lb1(ii), lb2(ii)
        if( Cal_comp ) then
          if( abvr(j) == 0 .and. abvi(j) == 0 ) cycle
        else
          if( abvr(j) == 0 ) cycle
        endif

        inz = inz+1
        rowIndexes(inz) = ii
        columnIndexes(inz) = j
        A(inz) = abvr(j)
        if( Cal_comp ) A_im(inz) = abvi(j)
      end do
 
      deallocate( abvi, abvr )
     
    end do      ! end of cycle by lines
    nz = inz
    
    if ( icheck > 0 ) write(3,100) nligne, nz
    if ( icheck > 1 ) then
      nligne8 = nligne
      write(6,'(" Sizes of linear equation system:")')
      write(6,'(" nligne   =",I24)') nligne
      write(6,'(" nligne^2 =",I24)') nligne8**2
      p2 = 100 * real(nz, db) / ( nligne8**2 )
      write(6,'(" not zero =",I24,5X,F6.3," %")') nz, p2
    endif

  else

    nz = 0
    nz_i = 0
    allocate( A( nz ), A_im(nz_i), rowIndexes(nz), columnIndexes(nz) )

  endif

  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    tp2 = real(time,db)
    tpt1 = tpt1 + tp2 - tp1
  endif
      
! run solver
  call mat_solver(A, A_im, rowIndexes, columnIndexes, smr, smi, nligne, nligne_i, nz, nz_i, nlmso, nlmso_i, Cal_comp, &
                  mpirank0, icheck)
  
  if ( mpirank_in_mumps_group == 0) then
    deallocate( rowIndexes )
    deallocate( columnIndexes )
    deallocate( A )
    deallocate( A_im )
  endif

  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    tp3 = real(time,db)
    tpt2 = tpt2 + tp3 - tp2
  endif
  
  return
  100 format(/' FDM matrix: number of line =',i7, / &
              '             number of not zero terms =',i8)
end

!**************************************************************************************************************

subroutine mat_solver(A, A_im, rowIndexes, columnIndexes, b, b_im, nligne, nligne_i, nz, nz_i, nlmso, nlmso_i, Cal_comp, &
                      mpirank0, icheck)
 
  use declarations
  implicit none
 
  INCLUDE 'dmumps_struc.h'
  INCLUDE 'zmumps_struc.h'
  
  integer:: nligne, nligne_i, nlmso, nlmso_i, MYID, mpirank0, i, j, icheck
  integer*8:: nz, nz_i
  integer*8, dimension(nz):: rowIndexes, columnIndexes
  
  logical:: Cal_comp

  real(kind=db), dimension(nz):: A
  real(kind=db), dimension(nz_i):: A_im
  real(kind=db), dimension(nlmso,nligne):: b
  real(kind=db), dimension(nlmso_i,nligne_i):: b_im
  
  TYPE (DMUMPS_STRUC) dmumps_par
  TYPE (ZMUMPS_STRUC) zmumps_par
  
  if ( Cal_comp ) then
! Complex case ========================================================

    zmumps_par%COMM = MPI_COMM_MUMPS
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
    zmumps_par%JOB = -1
    zmumps_par%SYM = 0
    zmumps_par%PAR = 1
    CALL ZMUMPS(zmumps_par)
    zmumps_par%ICNTL(2:3) = -1
    zmumps_par%ICNTL(4) = 2 !printLevel
    zmumps_par%ICNTL(7) = 5 !ordering
    zmumps_par%ICNTL(14) = 100 !memIncrease
    MYID = zmumps_par%MYID
!  Define problem on the host (processor 0)      
    IF ( MYID == 0 ) THEN
      zmumps_par%N = nligne
      zmumps_par%LRHS = nligne
      zmumps_par%NRHS = nlmso
      zmumps_par%NZ = nz
      ALLOCATE( zmumps_par%IRN ( nz ) )
      zmumps_par%IRN = rowIndexes
      ALLOCATE( zmumps_par%JCN ( nz ) )
      zmumps_par%JCN = columnIndexes
      ALLOCATE( zmumps_par%A( nz ) )
      zmumps_par%A = CMPLX(A, A_im, db)
      ALLOCATE( zmumps_par%RHS ( nligne*nlmso ) )
      do i = 1, nligne
        do j = 1, nlmso
          zmumps_par%RHS(i+nligne*(j-1)) = CMPLX( b(j,i), b_im(j,i), db )
        end do
      end do
    endif
!  Call package for solution
    zmumps_par%JOB = 6
    CALL ZMUMPS(zmumps_par)
    if ( zmumps_par%INFOG(1) /= 0 ) stop
!  Solution has been assembled on the host
    IF ( MYID == 0 ) THEN
      do j = 1, nlmso
        b(j,1:nligne) = real(zmumps_par%RHS(1+nligne*(j-1) : nligne*j), db)
        b_im(j,1:nligne) = aimag(zmumps_par%RHS(1+nligne*(j-1) : nligne*j))
      end do
!  Deallocate user data
      DEALLOCATE( zmumps_par%IRN )
      DEALLOCATE( zmumps_par%JCN )
      DEALLOCATE( zmumps_par%A   )
      DEALLOCATE( zmumps_par%RHS )
    END IF
!  Destroy the instance (deallocate internal data structures)
    zmumps_par%JOB = -2
    CALL ZMUMPS(zmumps_par)
    
  else
! Real case ===========================================================
    dmumps_par%COMM = MPI_COMM_MUMPS
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
    dmumps_par%JOB = -1
    dmumps_par%SYM = 0
    dmumps_par%PAR = 1
    CALL DMUMPS(dmumps_par)
    dmumps_par%ICNTL(2:3) = -1
    dmumps_par%ICNTL(4) = 2 !printLevel
    dmumps_par%ICNTL(7) = 5 !ordering
    dmumps_par%ICNTL(14) = 100 !memIncrease
    MYID = dmumps_par%MYID
!  Define problem on the host (processor 0)      
    IF ( MYID == 0 ) THEN
      dmumps_par%N = nligne
      dmumps_par%LRHS = nligne
      dmumps_par%NRHS = nlmso
      dmumps_par%NZ = nz
      ALLOCATE( dmumps_par%IRN ( nz ) )
      dmumps_par%IRN = rowIndexes
      ALLOCATE( dmumps_par%JCN ( nz ) )
      dmumps_par%JCN = columnIndexes
      ALLOCATE( dmumps_par%A( nz ) )
      dmumps_par%A = A
      ALLOCATE( dmumps_par%RHS ( nligne*nlmso ) )
      do i = 1, nligne
        do j = 1, nlmso
          dmumps_par%RHS(i+nligne*(j-1)) = b(j,i)
        end do
      end do
    endif
!  Call package for solution
    dmumps_par%JOB = 6
    CALL DMUMPS(dmumps_par)
    if ( dmumps_par%INFOG(1) /= 0 ) stop
!  Solution has been assembled on the host
    IF ( MYID == 0 ) THEN
      do j = 1, nlmso
        b(j,1:nligne) = dmumps_par%RHS(1+nligne*(j-1) : nligne*j)
      end do
!  Deallocate user data
      DEALLOCATE( dmumps_par%IRN )
      DEALLOCATE( dmumps_par%JCN )
      DEALLOCATE( dmumps_par%A   )
      DEALLOCATE( dmumps_par%RHS )
    END IF
!  Destroy the instance (deallocate internal data structures)
    dmumps_par%JOB = -2
    CALL DMUMPS(dmumps_par)

  endif
  
end
    
!**************************************************************************************        
        
subroutine getSolverParams(MPI_host_num_for_mumps,mpinodes0,Solver)

  use declarations
  implicit none

  integer:: ipr, MPI_host_num_for_mumps, mpi_split, mpinodes0
  
  character(len=5):: Solver
  character(len=255):: str
  
  Solver = 'MUMPS'
  
  CALL getenv("HOST_NUM_FOR_MUMPS", str)
  if ( LEN_TRIM(str) > 0 ) then
    read( str, '(i10)' ) MPI_host_num_for_mumps
  else
    MPI_host_num_for_mumps = 1
  endif

  mpi_split = mod( mpinodes0, MPI_host_num_for_mumps ) 
  if( mpi_split /= 0 ) then
    call write_error
    do ipr = 3,6,3
      write(6,100) mpinodes0, MPI_host_num_for_mumps
    end do
    stop
  endif

  return
  100 format(//' The total number of nodes =',i3,' must be a multiple of the number of nodes for MUMPS =',i3, &
               ' It is not the case !'//) 
end
