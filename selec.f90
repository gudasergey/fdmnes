! FDMNES subroutines
! Subprogram selecting spectra or azimuthal scan from the output scan
! file after convolution.

subroutine selec(itape5)

  use declarations
  implicit none

  integer:: eof, i, ia, iaa, iang, ibeam, ie, igrdat, indm, iref, istat, itape5, n, n_ang, n1, nbeam, ne, nf, nnombre

  character(len=9) grdat, grdat1
  character(len=Length_word), dimension(:), allocatable:: nombeam
  character(len=132) file_in, file_out, identmot, mots

  integer,dimension(1000):: n_angl
  integer, dimension(:), allocatable:: ind

  logical scan_true

  real(kind=db):: aa, aa1, e, ff, ff1, p, x
  real(kind=db), dimension(:), allocatable:: ang, energ
  real(kind=db), dimension(:,:,:), allocatable:: f

  ne = 0
  nf = 0
  n_ang = 0
  scan_true = .true.

  Rewind(itape5)

  do igrdat = 1,1000

    read(itape5,'(A)',iostat=eof) mots
    if( eof /= 0 ) exit

    grdat = identmot(mots,9)
    if( grdat(1:1) /= ' ' ) write(6,'(3x,A)') grdat

    select case(grdat)

      case('fin','end')
        exit

      case('selec_inp')
        n = nnombre(itape5,132)
        read(itape5,'(A)') file_in
! Start ESRF Changes
        if( file_in(1:1) == ' ' ) file_in = adjustl(file_in)
! End ESRF Changes

      case('selec_out')
        n = nnombre(itape5,132)
        read(itape5,'(A)') file_out
! Start ESRF Changes
        if( file_out(1:1) == ' ' ) file_out = adjustl(file_out)
! End ESRF Changes

      case('energy')
        do i = 1,1000
          n = nnombre(itape5,132)
          if( n == 0 ) exit
          read(itape5,*)
          ne = ne + n
        end do
        if( ne > 0 ) then
          allocate( energ(ne) )
          rewind(itape5)
          do i = 1,1000
            read(itape5,'(A)') mots
            grdat1 = identmot(mots,9)
            if( grdat1 == 'energy' ) exit
          end do
          n1 = 0
          do i = 1,ne
            n = nnombre(itape5,132)
            read(itape5,*) energ(n1+1:n1+n)
            n1 = n1 + n
            if( n1 == ne ) exit
          end do
        endif

      case('reflectio')
        do i = 1,1000
          n = nnombre(itape5,132)
          if( n == 0 ) exit
          read(itape5,*)
          nf = nf + n
        end do
        if( nf > 0 ) then
          allocate( ind(nf) )
          allocate( nombeam(nf) )
          rewind(itape5)
          do i = 1,1000
            read(itape5,'(A)') mots
            grdat1 = identmot(mots,9)
            if( grdat1 == 'reflectio' ) exit
          end do
          n1 = 0
          do iref = 1,nf
            n = nnombre(itape5,132)
            read(itape5,*) ind(n1+1:n1+n)
            n1 = n1 + n
            if( n1 == nf ) exit
          end do
        endif

      case('azimuth')
        do i = 1,1000
          n = nnombre(itape5,132)
          if( n == 0 ) exit
          read(itape5,*)
          n_ang = n_ang + n
        end do
        if( n_ang == 0 ) then
          scan_true = .true.
        else
          scan_true = .false.
          allocate( ang(n_ang) )
          rewind(itape5)
          do i = 1,1000
            read(itape5,'(A)') mots
            grdat1 = identmot(mots,9)
            if( grdat1 == 'azimuth' ) exit
          end do
          n1 = 0
          do iang = 1,n_ang
            n = nnombre(itape5,132)
            read(itape5,*) ang(n1+1:n1+n)
            n1 = n1 + n
            if( n1 == n_ang ) exit
          end do
        endif

    end select

  end do

  Close(itape5)

  if( ne == 0 .and. n_ang == 0 ) then
    call write_error
    write(6,100)
    write(9,100)
  stop
  endif
  if( ne /= 0 .and. n_ang > 0 ) then
    call write_error
    write(6,105)
    write(9,105)
  stop
  endif

  open(1, file = file_in, status = 'old', iostat=istat )
  if( istat /= 0 ) call write_open_error(file_in,istat,1)
  open(2, file = file_out, iostat=istat)
  if( istat /= 0 ) call write_open_error(file_out,istat,1)

! Determination of the number of reflections and number of angles
  n = nnombre(1,132)
  read(1,*) x
  read(1,'(A)') mots
  boucle_i: do ibeam = 1,10000
    boucle_ia: do ia = 1,100000

     n = nnombre(1,132)
     select case(n)
        case(2)
          read(1,*) x
        case(1)
          nbeam = ibeam
          n_angl(ibeam) = ia - 1
          exit boucle_i
        case(0)
          read(1,*) mots
          n_angl(ibeam) = ia - 1
          if( mots == ' ' ) exit boucle_i
          exit boucle_ia

      end select
    end do boucle_ia
  end do boucle_i

  if( nf == 0 ) then ! par defaut, on prend toutes les reflections
    do i = 1,nbeam
      if( n_angl(i) == 1 ) cycle
      nf = nf + 1
    end do
    allocate( ind(nf) )
    allocate( nombeam(nf) )
    iref = 0
    do i = 1,nbeam
      if( n_angl(i) == 1 ) cycle
      iref = iref + 1
      ind(iref) = i
    end do
  else
    indm = 1
    do iref = 1,nf
      indm = max( indm, ind(iref) )
    end do
    if( indm > nbeam ) then
      call write_error
      write(6,106) indm, nbeam
      write(9,106) indm, nbeam
      stop
    endif
  endif

  write(6,107) nbeam
  do ibeam = 1,nbeam
    write(6,108) ibeam, n_angl(ibeam)
  end do

  if( scan_true ) then
    n_ang = 1
    do iref = 1,nbeam
      if( n_angl(iref) > 1 ) then
        n_ang = n_angl(iref)
        exit
      endif
    end do
    do iref = 1,nf
      if( n_angl(ind(iref)) == n_ang ) cycle
      call write_error
      write(6,109) ind(iref), n_angl(ind(iref)), n_ang
      write(9,109) ind(iref), n_angl(ind(iref)), n_ang
      stop
    end do
  endif

  rewind(1)

  if( scan_true ) then
    allocate( ang(n_ang) )
  else
    boucle_ie: do ie = 1,10000
      read(1,*,iostat=eof)
      if( eof /= 0 ) exit boucle_ie
      read(1,*,iostat=eof) e
      if( eof /= 0 ) exit boucle_ie
      do ibeam = 1,nbeam
        read(1,*,iostat=eof)
        if( eof /= 0 ) exit boucle_ie
        do i = 1,n_angl(ibeam)
          read(1,*,iostat=eof)
          if( eof /= 0 ) exit boucle_ie
        end do
      end do
    end do boucle_ie
    ne = ie - 1
    allocate( energ(ne) )
    rewind(1)
  endif
  allocate( f(n_ang,nf,ne) )

! Lecture
  do ie = 1,ne

    do i = 1,100000
      n = nnombre(1,132)
      read(1,*,iostat=eof) e
      if( eof /= 0 ) exit
      if( scan_true ) then
        if( e > energ(ie) - 0.0001_db ) exit
      else
        exit
      endif
      do iref = 1,nbeam
        read(1,*)
        do ia = 1,n_angl(iref)
          read(1,*)
        end do
      end do
    end do

    energ(ie) = e

    do ibeam = 1,nbeam

      do iref = 1,nf
        if( ind(iref) == ibeam ) exit
      end do

      if( iref > nf ) then
        read(1,*)
        do ia = 1,n_angl(ibeam)
          read(1,*)
        end do
        cycle
      endif

      read(1,'(A13)') nombeam(iref)
      if( scan_true ) then
      do ia = 1,n_ang
          read(1,*) ang(ia), f(ia,iref,ie)
        end do
      else
        aa1 = 10000.
        do iaa = 1,n_angl(ibeam)
          read(1,*) aa, ff
          do ia = 1,n_ang
            if( abs( aa - ang(ia) ) < 0.0001 ) then
              f(ia,iref,ie) = ff
            elseif( ang(ia) < aa .and. ang(ia) > aa1 + 0.0001 ) then
              p = ( ang(ia) - aa1 ) / ( aa - aa1 )
              f(ia,iref,ie) = p * ff + (1 - p ) * ff1
            endif
          end do
          aa1 = aa
          ff1 = ff
        end do
      endif

    end do

  end do

  do iref = 1,nf
    mots = nombeam(iref)
    call center_word( mots, Length_word )
    nombeam(iref) = mots
  end do

  if( scan_true ) then
    if( nf == 1 ) then
      write(2,130) energ(1:ne)
      do i = -1,0
        do ia = 1,n_ang
          write(2,140) ang(ia)+i*360, f(ia,1,1:ne)
        end do
      end do
    else
      write(2,150) (nombeam(iref), iref = 1,nf)
      do i = -1,0
        do ia = 1,n_ang
          write(2,140) ang(ia)+i*360, f(ia,1:nf,1)
        end do
      end do
    endif
  else        ! Spectre a plusieurs angles
    if( n_ang == 1 ) then
      write(2,160) (nombeam(iref), iref = 1,nf)
    else
      write(2,170) ang(1:n_ang)
    endif
    do ie = 1,ne
      write(2,140) energ(ie), (f(1:n_ang,iref,ie), iref = 1,nf)
    end do
  endif

  deallocate( ang )
  deallocate( energ )
  deallocate( f )
  deallocate( ind )
  deallocate( nombeam )

  Close(1)
  Close(2)

  100 format(///' One has to choose azimuthal scan or spectra !'/, ' Number of energy or number of angle must not be zero', &
         ' together !'///)
  105 format(///' One has to choose azimuthal scan or spectra !'/, &
         ' Number of energy and number of angle cannot be non zero', ' together !'///)
  106 format(///' Feflection index =',i3,' > Reflection number =',i3///)
  107 format(/' nbeam =',i4)
  108 format(' ibeam =',i4,',   n_angle =',i4)
  109 format(///' nb_angle(',i3,') =',i4,' /= n_ang =',i4///)
  130 format('     Angle ',100f13.3)
  140 format(f11.3,1p,100e13.5)
  150 format('     Angle ',100a13)
  160 format('   Energy ',100a13)
  170 format('   Energy ',100f13.3)
end

