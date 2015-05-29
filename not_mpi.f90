! FDMNES subroutines
! These are fake MPI routines that return trivial values when
! running the non-parallel version of fdmnes.

subroutine MPI_Init(mpierr)
  mpierr = 0
  return
end
!----------------------------------------------------------------------
subroutine MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierr)
  mpirank = 0
  mpierr = 0
  return
end
!----------------------------------------------------------------------
subroutine MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  mpierr = 0
  return
end
!----------------------------------------------------------------------
subroutine MPI_BCAST(ivar,idim,MPI_INTEGER,ii,MPI_COMM_WORLD,mpierr)
  dimension ivar(idim)
  mpierr = 0
  return
end
!----------------------------------------------------------------------
subroutine MPI_Finalize(mpierr)
  mpierr = 0
  return
end
!----------------------------------------------------------------------
subroutine MPI_Comm_Size(MPI_COMM_WORLD,mpinodes,mpierr)
  mpinodes = 1
  mpierr = 0
  return
end
!----------------------------------------------------------------------
subroutine MPI_SEND(ivar,idim,MPI_INTEGER,mpirank,i,MPI_COMM_WORLD,mpierr)
  dimension ivar(idim)
  mpierr = 0
  i = 0
  return
end
!----------------------------------------------------------------------
subroutine MPI_RECV(ivar,idim,MPI_INTEGER,mpirank,i,MPI_COMM_WORLD,mstatus,mpierr)
  dimension ivar(idim)
  mpierr = 0
  mstatus = 0
  i = 0
  return
end
!----------------------------------------------------------------------
subroutine MPI_GATHER(MPI_IN_PLACE,idim,MPI_INTEGER,var,idim1,MPI_REAL8,i,MPI_COMM_WORLD,mpierr)
  dimension var(idim)
  mpierr = 0
  i = 0
  return
end
!----------------------------------------------------------------------
subroutine MPI_COMM_SPLIT(MPI_COMM_WORLD,i,j,MPI_COMM_MUMPS,mpierr)
  mpierr = 0
  return
end


