
! Routine solving the system of linear equations using sparse matrix technique and MUMPS
! From Alexander Guda et al, Rostov, Russia. 

subroutine mat_solve(Base_hexa, Basereel, Bessel, Besselr, Cal_comp, cgrad, clapl, E_comp, Eimag, Enervide, gradvr, &
        ianew, iato, ibord, icheck, ie, igrph, ii, isbord, iso, ispinin, isrt, isvois, ivois, Kar, Kari, lato, &
        lb1, lb2, lmaxso, lso, mato, MPI_host_num_for_mumps, mpirank0, mso, natome, nbm, nbord, nbordf, nbtm, Neuman, Neumanr, &
        new, newinv, ngrph, nicm, nim, nligne, nligne_i, nligneso, nlmsam,  nlmagm, nlmmax, nlmomax, nlmsa, nlmso, nlmso_i, &
        nphiato1, nphiato7, npoint, npsom, nsm, nso1, nsort, nsort_c, nsort_r, nsortf, nspin, nspino, nspinp, nspinr, nstm, &
        numia, nvois, phiato, poidsa, poidso, Relativiste, Repres_comp, rvol, smi, smr, Spinorbite, Time_fill, Time_tria, &
        Vr, Ylm_comp, Ylmato, Ylmso)

  use declarations
  implicit none
  include 'mpif.h'
  
  INTERFACE 
    subroutine expandArrayINT(x, n, minsz)
      integer(kind=8):: n, minsz
      integer, dimension(:), pointer:: x
    end subroutine
    
    subroutine expandArrayREAL(x, n, minsz)
      use declarations
      integer(kind=8):: n,minsz
      real(kind=db), dimension(:), pointer:: x
    end subroutine
    
    subroutine expandArrayCOMPLEX(x, n, minsz)
      use declarations
      integer(kind=8):: n,minsz
      complex(kind=db), dimension(:), pointer:: x
    end subroutine
    
    subroutine gatherINT(xGath, nzSum, x, nz, nzGather, mpirank, mpinodes)
      use declarations
      integer:: mpirank, mpinodes
      integer*8:: nzSum, nz
      integer*8, dimension(0:mpinodes-1):: nzGather
      integer, dimension(:), pointer:: xGath, x
    end subroutine
  
    subroutine gatherREAL(xGath, nzSum, x, nz, nzGather, mpirank, mpinodes)
      use declarations
      integer:: mpirank, mpinodes
      integer*8:: nzSum, nz
      integer*8, dimension(0:mpinodes-1):: nzGather
      real(kind=db), dimension(:), pointer:: xGath, x
    end subroutine
  
    subroutine gatherCOMPLEX(xGath, nzSum, x, nz, nzGather, mpirank, mpinodes)
      use declarations
      integer:: mpirank, mpinodes
      integer*8:: nzSum, nz
      integer*8, dimension(0:mpinodes-1):: nzGather
      complex(kind=db), dimension(:), pointer:: xGath, x
    end subroutine
  
    subroutine mat_solver(A, AZ, rowIndexes, columnIndexes, b, b_im, nligne, nligne_i, nz, nlmso, nlmso_i, Cal_comp, &
                          mpirank0, mpirank_in_mumps_group, MPI_host_num_for_mumps, icheck)
      use declarations
      integer:: nligne, nligne_i, nlmso, nlmso_i, mpirank0, mpirank_in_mumps_group, icheck
      integer*8:: nz
      integer, dimension(:), pointer:: rowIndexes, columnIndexes
      logical:: Cal_comp
      real(kind=db), dimension(:), pointer:: A
      complex(kind=db), dimension(:), pointer:: AZ
      real(kind=db), dimension(nlmso,nligne):: b
      real(kind=db), dimension(nlmso_i,nligne_i):: b_im
    end subroutine
  END INTERFACE

  integer:: icheck, ie, igrph, ii, isp, ispin, ispinin, i, j, lb1i, lb1r, lb2i, lb2r, &
    lmaxso, MPI_host_num_for_mumps, mpirank_in_mumps_group, mpirank0, natome, nbm, nbtm, ngrph, nicm, nim, nligne, &
    nligne_i, nligneso, nlmagm, nlmmax, nlmomax, nlmsam, nlmso, nlmso_i, nphiato1, nphiato7, npoint, & 
    npsom, nsm, nso1, nsort, nsort_c, nsort_r, nsortf, nspin, nspino, nspinp, nspinr, nstm, nvois, mpierr
  
  integer, dimension(0:npoint):: new
  integer, dimension(natome):: ianew, nbord, nbordf, nlmsa
  integer, dimension(nstm):: isrt
  integer, dimension(npsom):: numia
  integer, dimension(nso1,ngrph):: lso, mso, iso
  integer, dimension(nlmsam,natome,ngrph):: lato, mato, iato
  integer, dimension(nbtm,natome):: ibord, isbord
  integer, dimension(npsom,nvois):: ivois, isvois
  integer, dimension(nligne):: lb1, lb2, newinv

  integer(kind=8):: inz, nz, nzSum, nligne8, oldSize
  integer(kind=8), dimension(0:MPI_host_num_for_mumps-1):: nzGather
  integer, dimension(0:MPI_host_num_for_mumps-1):: displs, nsmrGather
  integer, dimension(:), pointer:: rowIndexes, rowIndexesGather, columnIndexes, columnIndexesGather

  character(len=2), dimension(nligne):: mletl

  complex(kind=db), dimension(nsort_c,0:lmaxso,nspinr):: Bessel, Neuman

  logical:: Base_hexa, Basereel, Cal_comp, E_comp, Relativiste, Repres_comp, Spinorbite, Ylm_comp

  real(kind=sg):: time

  real(kind=db):: Enervide, Eimag, p2, Precision, tp1, tp2, tp3, Time_fill, Time_tria
  
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

  real(kind=db), dimension(:), allocatable:: abvi, abvr
  real(kind=db), dimension(:), pointer:: A, AGather 
  complex(kind=db), dimension(:), pointer:: AZ, AZGather

  mpirank_in_mumps_group = mod( mpirank0, MPI_host_num_for_mumps )

!  Define problem in parallel by mumps group
  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    tp1 = real(time,db)
  endif
  
  if( Spinorbite ) then
    ispin = 1
  else
    ispin = ispinin
  endif

  nz = nligne/1000*nligne/MPI_host_num_for_mumps + 1000000 !start value for nz
  allocate( rowIndexes ( nz ) )
  allocate( columnIndexes ( nz ) )
  if ( Cal_comp ) then
    allocate( AZ( nz ) )
  else
    allocate( A( nz ) )
  endif
  smr(:,:) = 0._db
  if( Cal_comp ) smi(:,:) = 0._db  
  
  inz = 0
  do ii = mpirank_in_mumps_group+1, nligne, MPI_host_num_for_mumps
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
    smi, smr, Vr, Ylm_comp, Ylmato, Ylmso, .true. )
    
! expand arrays            
    oldSize = size(rowIndexes)
    if ( oldSize < inz+lb2r-lb1r+1 ) then
      call expandArrayINT(rowIndexes, oldSize, inz+lb2r-lb1r+1)
      call expandArrayINT(columnIndexes, oldSize, inz+lb2r-lb1r+1)
      if( Cal_comp ) then
        call expandArrayCOMPLEX(AZ, oldSize, inz+lb2r-lb1r+1)
      else
        call expandArrayREAL(A, oldSize, inz+lb2r-lb1r+1)
      endif
    endif

    Precision = 1.e-20_db
        
! Fill matrises
    do j = lb1(ii), lb2(ii)
      if( Cal_comp ) then
        if( abs( abvr(j) ) < Precision .and. abs( abvi(j) ) < Precision ) cycle
      else
        if( abs( abvr(j) ) < Precision ) cycle
      endif
!      if( Cal_comp ) then
!        if( abvr(j) == 0 .and. abvi(j) == 0 ) cycle
!      else
!        if( abvr(j) == 0 ) cycle
!      endif

      inz = inz+1
      rowIndexes(inz) = ii
      columnIndexes(inz) = j
      if( Cal_comp ) then
        AZ(inz) = CMPLX(abvr(j), abvi(j), db)
      else
        A(inz) = abvr(j)
      endif
    end do
    deallocate( abvi, abvr )
  end do      ! end of cycle by lines
  nz = inz
  
  if ( MPI_host_num_for_mumps > 1 ) then
    call MPI_GATHER( nz, 1, MPI_INTEGER8, nzGather, 1, MPI_INTEGER8, 0, MPI_COMM_MUMPS, mpierr)
    nzSum = sum(nzGather)
    call gatherINT(rowIndexesGather, nzSum, rowIndexes, nz, nzGather, mpirank_in_mumps_group, MPI_host_num_for_mumps)
    call gatherINT(columnIndexesGather, nzSum, columnIndexes, nz, nzGather, mpirank_in_mumps_group, MPI_host_num_for_mumps)
    if( Cal_comp ) then
      call gatherCOMPLEX(AZGather, nzSum, AZ, nz, nzGather, mpirank_in_mumps_group, MPI_host_num_for_mumps)
    else
      call gatherREAL(AGather, nzSum, A, nz, nzGather, mpirank_in_mumps_group, MPI_host_num_for_mumps)
    endif
    call gather_sm(smr,smi,nlmso,nligne,nlmso_i,nligne_i,MPI_host_num_for_mumps,mpirank_in_mumps_group,Cal_comp)
    nz = nzSum
  endif
  if ( mpirank0 == 0 ) then
    if ( icheck > 1 .or. ( icheck == 1 .and. ie == 1 ) ) then
      if( ngrph == 1 ) then
        write(3,100) nligne, nz
      else
        write(3,110) igrph, nligne, nz
      endif
    endif
    if ( icheck > 2 ) then
      nligne8 = nligne
      write(6,'(" Sizes of linear equation system:")')
      write(6,'(" nligne   =",I24)') nligne
      write(6,'(" nligne^2 =",I24)') nligne8**2
      p2 = 100 * real(nz, db) / ( nligne8**2 )
      write(6,'(" not zero =",I24,5X,F6.3," %")') nz, p2
    endif
  endif

  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    tp2 = real(time,db)
    Time_fill = Time_fill + tp2 - tp1
  endif
      
! run solver
  call mat_solver(A, AZ, rowIndexes, columnIndexes, smr, smi, nligne, nligne_i, &
                nz, nlmso, nlmso_i, Cal_comp, mpirank0, mpirank_in_mumps_group, MPI_host_num_for_mumps, icheck)
 
  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    tp3 = real(time,db)
    Time_tria = Time_tria + tp3 - tp2
  endif
  
  return
  100 format(/' FDM matrix: number of line =',i7,',  number of not zero terms =',i8)
  110 format(/' FDM matrix, igrph =',i2,' : number of line =',i7,',  number of not zero terms =',i8)
end

!**************************************************************************************************************

subroutine expandArrayINT(x,n,minsz)
  use declarations
  implicit none
  integer(kind=8):: n,minsz,newsz
  integer, dimension(:), pointer:: x, x_new
  newsz = n*2
  if (newsz < minsz) newsz=minsz
  allocate(x_new(newsz))
  x_new(1:n) = x(1:n)
  deallocate(x)
  x => x_new
end

!**************************************************************************************************************

subroutine expandArrayREAL(x,n,minsz)
  use declarations
  implicit none
  integer(kind=8):: n,minsz,newsz
  real(kind=db), dimension(:), pointer:: x, x_new
  newsz = n*2
  if (newsz < minsz) newsz=minsz
  allocate(x_new(newsz))
  x_new(1:n) = x(1:n)
  deallocate(x)
  x => x_new
end

!**************************************************************************************************************

subroutine expandArrayCOMPLEX(x,n,minsz)
  use declarations
  implicit none
  integer(kind=8):: n,minsz,newsz
  complex(kind=db), dimension(:), pointer:: x, x_new
  newsz = n*2
  if (newsz < minsz) newsz=minsz
  allocate(x_new(newsz))
  x_new(1:n) = x(1:n)
  deallocate(x)
  x => x_new
end

!**************************************************************************************************************

subroutine gatherINT(xGath, nzSum, x, nz, nzGather, mpirank, mpinodes)
  use declarations
  implicit none
  include 'mpif.h'
  integer:: mpirank, mpinodes, sz, n, sendCount, isend, mpierr, j, rank
  integer*8:: nzSum, nz
  integer*8, dimension(0:mpinodes-1):: nzGather,displs, nzSent
  integer, dimension(:), pointer:: xGath,x
  integer status(MPI_STATUS_SIZE)
  
  if (mpirank == 0) allocate(xGath(nzSum))
  sz = 25000000 !max size for mpi send/recv
  nzSent(:) = 0
  
  displs(0)=0
  do j = 1,mpinodes-1
    displs(j) = displs(j-1)+nzGather(j-1); 
  end do
! gather array on root node of MUMPS group
! MPI_GATHERV can't gather more than 2*10^9 elements, so we do it by parts
  
  if ( mpirank == 0 ) then
    xGath(1:nz) = x(1:nz)
    do rank = 1,mpinodes-1
      sendCount = (nzGather(rank)+sz-1)/sz
      do isend = 1,sendCount
        n = MIN(sz, nzGather(rank)-nzSent(rank))
        call MPI_recv(xGath(1+displs(rank)+nzSent(rank)), n, MPI_INTEGER, rank, 123, MPI_COMM_MUMPS, status, mpierr)
        nzSent(rank) = nzSent(rank) + n
      end do
    end do
  else
    rank = mpirank
    sendCount = (nz+sz-1)/sz
    do isend = 1,sendCount
      n = MIN(sz, nz-nzSent(rank))
      call MPI_send(x(1+nzSent(rank)), n, MPI_INTEGER, 0, 123, MPI_COMM_MUMPS, mpierr)
      nzSent(rank) = nzSent(rank) + n
    end do
  endif
  deallocate(x)
  if (mpirank == 0) x => xGath
  return
end

!**************************************************************************************************************

subroutine gatherREAL(xGath, nzSum, x, nz, nzGather, mpirank, mpinodes)
  use declarations
  implicit none
  include 'mpif.h'
  
  integer:: mpirank, mpinodes, sz, n, sendCount, isend, mpierr, j, rank
  integer*8:: nzSum, nz
  integer*8, dimension(0:mpinodes-1):: nzGather,displs, nzSent
  real(kind=db), dimension(:), pointer:: xGath,x
  integer status(MPI_STATUS_SIZE)
  
  if (mpirank == 0) allocate(xGath(nzSum))
  sz = 25000000 !max size for mpi send/recv
  nzSent(:) = 0
  
  displs(0)=0
  do j = 1,mpinodes-1
    displs(j) = displs(j-1)+nzGather(j-1); 
  end do
  
! gather array on root node of MUMPS group
! MPI_GATHERV can't gather more than 2*10^9 elements, so we do it by parts
  
  if ( mpirank == 0 ) then
    xGath(1:nz) = x(1:nz)
    do rank = 1,mpinodes-1
      sendCount = (nzGather(rank)+sz-1)/sz
      do isend = 1,sendCount
        n = MIN(sz, nzGather(rank)-nzSent(rank))
        call MPI_recv(xGath(1+displs(rank)+nzSent(rank)), n, MPI_DOUBLE_PRECISION, rank, 123, MPI_COMM_MUMPS, status, mpierr)
        nzSent(rank) = nzSent(rank) + n
      end do
    end do
  else
    rank = mpirank
    sendCount = (nz+sz-1)/sz
    do isend = 1,sendCount
      n = MIN(sz, nz-nzSent(rank))
      call MPI_send(x(1+nzSent(rank)), n, MPI_DOUBLE_PRECISION, 0, 123, MPI_COMM_MUMPS, mpierr)
      nzSent(rank) = nzSent(rank) + n
    end do
  endif
  deallocate(x)
  if (mpirank == 0) x => xGath
  return
end

!**************************************************************************************************************

subroutine gatherCOMPLEX(xGath, nzSum, x, nz, nzGather, mpirank, mpinodes)
  use declarations
  implicit none
  include 'mpif.h'
  
  integer:: mpirank, mpinodes, sz, n, sendCount, isend, mpierr, j, rank
  integer*8:: nzSum, nz
  integer*8, dimension(0:mpinodes-1):: nzGather,displs, nzSent
  complex(kind=db), dimension(:), pointer:: xGath,x
  integer status(MPI_STATUS_SIZE)
  
  if (mpirank == 0) allocate(xGath(nzSum))
  sz = 25000000 !max size for mpi send/recv
  nzSent(:) = 0
  
  displs(0)=0
  do j = 1,mpinodes-1
    displs(j) = displs(j-1)+nzGather(j-1); 
  end do
  
! gather array on root node of MUMPS group
! MPI_GATHERV can't gather more than 2*10^9 elements, so we do it by parts
  
  if ( mpirank == 0 ) then
    xGath(1:nz) = x(1:nz)
    do rank = 1,mpinodes-1
      sendCount = (nzGather(rank)+sz-1)/sz
      do isend = 1,sendCount
        n = MIN(sz, nzGather(rank)-nzSent(rank))
        call MPI_recv(xGath(1+displs(rank)+nzSent(rank)), n, MPI_DOUBLE_COMPLEX, rank, 123, MPI_COMM_MUMPS, status, mpierr)
        nzSent(rank) = nzSent(rank) + n
      end do
    end do
  else
    rank = mpirank
    sendCount = (nz+sz-1)/sz
    do isend = 1,sendCount
      n = MIN(sz, nz-nzSent(rank))
      call MPI_send(x(1+nzSent(rank)), n, MPI_DOUBLE_COMPLEX, 0, 123, MPI_COMM_MUMPS, mpierr)
      nzSent(rank) = nzSent(rank) + n
    end do
  endif
  deallocate(x)
  if (mpirank == 0) x => xGath
  return
end

!**************************************************************************************************************

subroutine gather_sm(smr,smi,nlmso,nligne,nlmso_i,nligne_i,mpinodes,mpirank,Cal_comp)
  use declarations
  implicit none
  include 'mpif.h'
  
  integer:: mpirank, mpinodes, nlmso, nlmso_i, nligne, nligne_i, i, j, rank, mpierr
  real(kind=db), dimension(nlmso,nligne):: smr
  real(kind=db), dimension(:,:), allocatable:: recvBuf
  real(kind=db), dimension(nlmso_i,nligne_i):: smi
  integer status(MPI_STATUS_SIZE)
  logical:: Cal_comp
  
  if ( mpinodes == 1 ) return
  allocate(recvBuf(nlmso,nligne))
  if ( mpirank == 0 ) then
    do rank = 1,mpinodes-1
      call MPI_recv(recvBuf, nlmso*nligne, MPI_REAL8, rank, 100, MPI_COMM_MUMPS, status, mpierr)
      smr(:,rank+1:nligne:mpinodes) = recvBuf(:,rank+1:nligne:mpinodes)
      if ( Cal_comp ) then
        call MPI_recv(recvBuf, nlmso*nligne, MPI_REAL8, rank, 100, MPI_COMM_MUMPS, status, mpierr)
        smi(:,rank+1:nligne:mpinodes) = recvBuf(:,rank+1:nligne:mpinodes)
      endif
    end do
  else
    call MPI_send(smr, nlmso*nligne, MPI_REAL8, 0, 100, MPI_COMM_MUMPS, mpierr)
    if ( Cal_comp ) call MPI_send(smi, nlmso*nligne, MPI_REAL8, 0, 100, MPI_COMM_MUMPS, mpierr)
  endif
  deallocate(recvBuf)
  return
end

!**************************************************************************************************************

subroutine mat_solver(A, AZ, rowIndexes, columnIndexes, b, b_im, nligne, nligne_i, nz, nlmso, nlmso_i, Cal_comp, &
                      mpirank0, mpirank_in_mumps_group, MPI_host_num_for_mumps, icheck)
 
  use declarations
  implicit none
  
  !DEC$ OPTIONS /NOWARN
  INCLUDE 'dmumps_struc.h'
  INCLUDE 'zmumps_struc.h'
  !DEC$ END OPTIONS
  
  integer:: nligne, nligne_i, nlmso, nlmso_i, MYID, mpirank0, mpirank_in_mumps_group, i, j, icheck, par, &
            MPI_host_num_for_mumps, ipr
  integer*8:: nz
  integer, dimension(:), pointer:: rowIndexes, columnIndexes
  
  logical:: Cal_comp

  real(kind=db), dimension(:), pointer:: A
  complex(kind=db), dimension(:), pointer:: AZ
  real(kind=db), dimension(nlmso,nligne):: b
  real(kind=db), dimension(nlmso_i,nligne_i):: b_im
  
  TYPE (DMUMPS_STRUC) dmumps_par
  TYPE (ZMUMPS_STRUC) zmumps_par
  
  character(len=255):: str
  
  par = 1
! for many nodes release root node from calculation to save memory
  if ( MPI_host_num_for_mumps > 4 ) par = mpirank_in_mumps_group
  CALL getenv("MUMPS_CENTRAL_WORK", str)
  if ( LEN_TRIM(str) > 0 ) then
    read( str, '(i10)' ) i
    if ( i == 0 ) par = mpirank_in_mumps_group
  endif
  
  
  if ( Cal_comp ) then
! Complex case ========================================================

    zmumps_par%COMM = MPI_COMM_MUMPS
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
    zmumps_par%JOB = -1
    zmumps_par%SYM = 0
    zmumps_par%PAR = par
    CALL ZMUMPS(zmumps_par)
    zmumps_par%ICNTL(2:3) = -1
    zmumps_par%ICNTL(4) = 2 !printLevel
    zmumps_par%ICNTL(7) = 5 !ordering (meaningles when ICNTL(28)=2)
    zmumps_par%ICNTL(14) = 100 !memIncrease
    zmumps_par%ICNTL(28) = 1 ! 2 - use parallel ordering, 1 - sequential
    zmumps_par%ICNTL(29) = 2 ! parallel ordering 
    MYID = zmumps_par%MYID
!  Define problem on the host (processor 0)      
    IF ( MYID == 0 ) THEN
      zmumps_par%N = nligne
      zmumps_par%LRHS = nligne
      zmumps_par%NRHS = nlmso
      zmumps_par%NZ = nz
      zmumps_par%IRN => rowIndexes
      zmumps_par%JCN => columnIndexes
      zmumps_par%A => AZ
    endif
!  Call package for solution
    zmumps_par%JOB = 4
    CALL ZMUMPS(zmumps_par)
    if ( zmumps_par%INFOG(1) /= 0 ) then
      if ( zmumps_par%INFOG(1) == -13 ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,*) 'Not enough memory'
        end do
      endif
      stop 1
    endif
    IF ( MYID == 0 ) THEN
      DEALLOCATE( zmumps_par%IRN )
      DEALLOCATE( zmumps_par%JCN )
      DEALLOCATE( zmumps_par%A   )
      ALLOCATE( zmumps_par%RHS ( nligne*nlmso ) )
      do i = 1, nligne
        do j = 1, nlmso
          zmumps_par%RHS(i+nligne*(j-1)) = CMPLX( b(j,i), b_im(j,i), db )
        end do
      end do
    endif
    zmumps_par%JOB = 3
    CALL ZMUMPS(zmumps_par)
!  Solution has been assembled on the host
    IF ( MYID == 0 ) THEN
      do j = 1, nlmso
        b(j,1:nligne) = real(zmumps_par%RHS(1+nligne*(j-1) : nligne*j), db)
        b_im(j,1:nligne) = aimag(zmumps_par%RHS(1+nligne*(j-1) : nligne*j))
      end do
!  Deallocate user data
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
    dmumps_par%PAR = par
    CALL DMUMPS(dmumps_par)
    dmumps_par%ICNTL(2:3) = -1
    dmumps_par%ICNTL(4) = 2 !printLevel
    dmumps_par%ICNTL(7) = 5 !ordering (meaningles when ICNTL(28)=2)
    dmumps_par%ICNTL(14) = 100 !memIncrease
    dmumps_par%ICNTL(28) = 1 ! 2 - use parallel ordering, 1 - sequential
    dmumps_par%ICNTL(29) = 2 ! parallel ordering 
    MYID = dmumps_par%MYID
!  Define problem on the host (processor 0)      
    IF ( MYID == 0 ) THEN
      dmumps_par%N = nligne
      dmumps_par%LRHS = nligne
      dmumps_par%NRHS = nlmso
      dmumps_par%NZ = nz
      dmumps_par%IRN => rowIndexes
      dmumps_par%JCN => columnIndexes
      dmumps_par%A => A
    endif
!  Call package for solution
    dmumps_par%JOB = 4
    CALL DMUMPS(dmumps_par)
    if ( dmumps_par%INFOG(1) /= 0 ) then
      if ( dmumps_par%INFOG(1) == -13 ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,*) 'Not enough memory'
        end do
      endif
      stop 1
    endif
    if ( MYID == 0 ) then
      DEALLOCATE( dmumps_par%IRN )
      DEALLOCATE( dmumps_par%JCN )
      DEALLOCATE( dmumps_par%A   )
      ALLOCATE( dmumps_par%RHS ( nligne*nlmso ) )
      do i = 1, nligne
        do j = 1, nlmso
          dmumps_par%RHS(i+nligne*(j-1)) = b(j,i)
        end do
      end do
    endif
    dmumps_par%JOB = 3
    CALL DMUMPS(dmumps_par)
!  Solution has been assembled on the host
    IF ( MYID == 0 ) THEN
      do j = 1, nlmso
        b(j,1:nligne) = dmumps_par%RHS(1+nligne*(j-1) : nligne*j)
      end do
!  Deallocate user data
      DEALLOCATE( dmumps_par%RHS )
    END IF
!  Destroy the instance (deallocate internal data structures)
    dmumps_par%JOB = -2
    CALL DMUMPS(dmumps_par)
  endif
  return
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
      write(ipr,100) mpinodes0, MPI_host_num_for_mumps
    end do
    stop 1
  endif

  return
  100 format(//' The total number of nodes =',i3,' must be a multiple of the number of nodes for MUMPS =',i3, &
               ' It is not the case !'//) 
end


