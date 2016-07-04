! FDMNES subroutines

! Program giving all the atoms of the mesh from the non equivalent atoms and from the space group.
! A big part comes from Ch. Brouder.

!   Space_Group      Name of the symmetry group
!   NMAXOP      Maximum number of symmetry operations

subroutine spgroup(Cif,Do_exp,neq,ngroup,ngroup_neq,itype,posn,posout,Space_file,Space_group)

  use declarations
  implicit none

  integer, parameter:: nmaxop = 192

  integer:: i, ia, is, j, ja, js, k, ks, ngroup, ngroup_neq, nsym
  
  character(len=1):: SGTrans
  character(len=13):: Space_Group
  character(len=132):: Space_file

  integer, dimension(ngroup):: itype
  integer, dimension(ngroup_neq):: neq
  integer, dimension(:), allocatable:: itypeq

  logical:: Check, Cif, Do_exp

  real(kind=db):: eps
  real(kind=db), dimension(3):: Along_XY, Along_YZ, Along_XZ, Along_XYZ, q
  real(kind=db), dimension(3,3,nmaxop):: Mat
  real(kind=db), dimension(3,nmaxop):: Trans
  real(kind=db), dimension(3,ngroup_neq):: posn
  real(kind=db), dimension(3,ngroup):: posout
  real(kind=db), dimension(:,:), allocatable:: qq

  Along_XY(1) = 1._db; Along_XY(2) = 1._db; Along_XY(3) = 0._db;
  Along_YZ(1) = 0._db; Along_YZ(2) = 1._db; Along_YZ(3) = 1._db;
  Along_XZ(1) = 1._db; Along_XZ(2) = 0._db; Along_XZ(3) = 1._db;
  Along_XYZ(1) = 1._db; Along_XYZ(2) = 1._db; Along_XYZ(3) = 1._db;

  Mat(:,:,:) = 0._db
  Trans(:,:) = 0._db
  Check = .false.

  call symgrp(Cif,Space_Group,Mat,Trans,nsym,nmaxop,SGTrans,Space_file)

  select case(SGTrans)

    case('A')

      do is = 1,nsym
        js = is + nsym
        Mat(:,:,js) = Mat(:,:,is)
        Trans(:,js) = Trans(:,is) + 0.5_db * Along_YZ(:)
      end do
      nsym = 2 * nsym

    case('B')

      do is = 1,nsym
        js = is + nsym
        Mat(:,:,js) = Mat(:,:,is)
        Trans(:,js) = Trans(:,is) + 0.5_db * Along_XZ(:)
      end do
      nsym = 2 * nsym

    case('C')

      do is = 1,nsym
        js = is + nsym
        Mat(:,:,js) = Mat(:,:,is)
        Trans(:,js) = Trans(:,is) + 0.5_db * Along_XY(:)
      end do
      nsym = 2 * nsym

    case('F')

      do is = 1,nsym
        do k = 2,4
          js = is + ( k - 1 ) * nsym
          Mat(:,:,js) = Mat(:,:,is)
          select case(k)
            Case(2)
              Trans(:,js) = Trans(:,is) + 0.5_db * Along_YZ(:)
            Case(3)
              Trans(:,js) = Trans(:,is) + 0.5_db * Along_XZ(:)
            Case(4)
              Trans(:,js) = Trans(:,is) + 0.5_db * Along_XY(:)
          end select
        end do
      end do
      nsym = 4 * nsym

    case('I')

      do is = 1,nsym
        js = is + nsym
        Mat(:,:,js) = Mat(:,:,is)
        Trans(:,js) = Trans(:,is) + 0.5_db * Along_XYZ(:)
      end do
      nsym = 2 * nsym

    case('H')

      do is = 1,nsym
        do k = 2,3
          js = is + ( k - 1 ) * nsym
          Mat(:,:,js) = Mat(:,:,is)
          Trans(:,js) = Trans(:,is)
          select case(k)
            Case(2)
              Trans(1,js) = Trans(1,js) + 2._db / 3
              Trans(2:3,js) = Trans(2:3,js) + 1._db / 3
            Case(3)
              Trans(1,js) = Trans(1,js) + 1._db / 3
              Trans(2:3,js) = Trans(2:3,js) + 2._db / 3
          end select
        end do
      end do
      nsym = 3 * nsym

  end select

  if( check ) then
    write(3,'(/A)') '        Matrix           Trans'
    do is = 1,nsym
      write(3,'(A,i3)') '  is =', is
      do i = 1,3
        write(3,'(3f7.3,3x,f7.3)') Mat(i,:,is), Trans(i,is)
      end do
    end do
  endif

  eps = 0.0000001_db
  do j = 1,2
    where( posn > 1._db - eps ) posn = posn - 1
    where( posn < - eps ) posn = posn + 1
  end do

  allocate( qq(3,nmaxop*ngroup_neq) )
  allocate( itypeq(nmaxop*ngroup_neq) )

  js = 0
  do ia = 1,ngroup_neq
    ja = 0
    boucle_is: do is = 1,nsym
      do i = 1,3
        q(i) = sum( Mat(i,:,is) * posn(:,ia) ) + Trans(i,is)
      end do
      do j = 1,2
        where( q > 1._db - eps ) q = q - 1
        where( q < - eps ) q = q + 1
      end do
      do ks = 1,js
        if( sum( abs(qq(:,ks)-q(:)) ) < 0.000001_db .and. itypeq(ks) == itype(ia) ) cycle boucle_is
      end do
      js = js + 1
      ja = ja + 1
      qq(:,js) = q(:)
      itypeq(js) = itype(ia)
      if( Do_exp ) posout(:,js ) = qq(:,js)
    end do boucle_is
    if( Do_exp ) neq(ia) = ja
  end do

  deallocate( itypeq, qq )
  if( .not. Do_exp ) ngroup = js

  return
end

!*********************************************************************

subroutine symgrp(Cif,Space_Group,Mat,Trans,nbsyop,nmaxop,SGTrans,Space_file)

! This subroutine looks for the space group whose name is in Space_Group.
! If it is found, it outputs the number of symmetry operations, and builds the matrices for these operations.

! Variable description :
!   Space_Group Name of the symmetry group
!   NBSYOP      Number of symmetry operations
!   NMAXOP      Maximum number of symmetry operations
!   sgnb        Space group number

  use declarations
  implicit none

  integer:: eof, nmaxop

  character(len=1):: SGTrans
  character(len=10):: sgnbcar, sgnbcar0
  character(len=13):: sgschoenfliess, Space_Group
  character(len=27):: sgHMshort, sgHMlong
  character(len=80):: line, mot
  character(len=132):: Space_file
  character(len=80), dimension(nmaxop):: lines

  integer:: i, i1, i2, ipr, istat, itape, j, nbsyop, sgnb

  logical:: Cif, pareil

  real(kind=db), dimension(3,3,nmaxop):: Mat
  real(kind=db), dimension(3,nmaxop):: Trans
  real(kind=db), dimension(3,4):: Matrix(3,4)

  pareil = .false.
  itape = 7

  Open(itape, file = Space_file, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(Space_file,istat,1)

  if( Cif ) then

    do
     read(itape,'(A)') mot
     if( mot(1:26) == '_symmetry_equiv_pos_as_xyz' .or. mot(1:32) == '_space_group_symop_operation_xyz' ) exit 
    end do
    
    boucle_i: do i = 1,1000

      read(itape,'(A)' ) mot
      line = ' '
      
      do j = 1,80
        if( mot(j:j) == ',' ) exit
      end do
      if( j > 80 ) exit boucle_i

      do j = 1,80
        if( mot(j:j) == "'" ) mot(j:j) = ' '
      end do

      lines(i) = adjustl(mot)
        
    end do boucle_i
    
    nbsyop = i - 1

    SGTrans = ' ' ! Space_group(1:1)
        
  else
  
! Ask for exact definition of space group.
! sgnbcar0 is the detailed number of the spacegroup, specifying axis and origin conventions

    call locateSG(itape,Space_file,Space_Group,sgnbcar0)

    Rewind(itape)

! In spacegroup.txt the name of the symmetry group follows a *<space>
! look for it. If a * is found, check that the following string
! is the name of the desired symmetry group.

    do i = 1,10000

      Read(itape,'(a80)',iostat=istat) line
      if( istat /= 0 ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,100) Space_group, Trim( Space_file )
        end do
        stop
      endif
      if (line(1:1) /= '*') cycle

! Analyse the line giving space group name(s)
      call analysename(line,sgnb,sgnbcar,sgschoenfliess,sgHMshort,sgHMlong)
      Pareil = sgnbcar == sgnbcar0
      SGTrans = sgHMlong(1:1)
      if( index(sgHMlong,'H') /= 0 ) SGTrans = 'H'
      if( pareil ) exit

    end do

! Look for nbsyop
    do i = 1,1000
      read(itape,'(a80)',iostat=eof) line
      if( eof /= 0 ) exit
      if( line(1:1) == '*' .or. line(1:1) == ' ' ) exit
      lines(i) = line
    end do
    nbsyop = i - 1

    if( index(sgnbcar,'R') /= 0 ) then
      call write_error
      do ipr = 6,9,3
        write(ipr,110) Space_file
      end do
      stop
    end if

  endif
   
! Read the NBSYOP symmetry operations, and build the corresponding
! transformation matrix.

  do i = 1,nbsyop
    call findop(lines(i),Matrix)
    do i1=1,3
      do i2=1,3
        Mat(i1,i2,i) = Matrix(i1,i2)
      end do
      Trans(i1,i) = Matrix(i1,4)
    end do
  end do

  Close(itape)

  return
  100 format(//' Space group name, ',a13,', not found in the file ',A//)
  110 format(//' Rhombohedral axes are not implemented in the file ',A,/ ' Please, convert to hexagonal axes'//)
end

!***********************************************************************

subroutine findop(line,matrix)

! This subroutine takes the line LINE coming from file spacegroup.txt
! and builds the matrix corresponding to the symmetry operation
! written in the line.
! The symmetry operation is written in the line as e.g. -y,x,-z

!   line    Line red from the input file
!   Matrix  Matrix of the symmetry operation (output) in the
!           crystal axes

  use declarations
  implicit none

  character(len=80) line

  integer i, j, ibegin, iend, ifin, ipr, ncar

  real(kind=db), dimension(3,4):: matrix

! Initialize Matrix
  do i=1,3
    do j=1,4
      Matrix(i,j) = 0._db
    end do
  end do
  ibegin = 1
  iend = len_trim(line)

! The symmetry operation is written as a line of 3 words
! separated by commas e.g. -y,x,-z
!   NCAR is the number of characters in each word
!   IBEGIN the place of the beginning of the word
!   IFIN   the place of the end of the word
!   IEND   the place of the end of the line

  do i = 1,3

    ncar = Index(line(ibegin:iend),',') - 1
    ifin = ibegin + ncar - 1
    if( ncar == -1 ) ifin = iend

    select case( line(ibegin:ifin) )

      case('x')
        Matrix(i,1) = 1._db

      case('-x')
        Matrix(i,1) = -1._db

      case('y')
        Matrix(i,2) = 1._db

      case('-y')
        Matrix(i,2) = -1._db

      case('z')
        Matrix(i,3) = 1._db

      case('-z')
        Matrix(i,3) = -1._db

      case('x-y')
        Matrix(i,1) = 1._db
        Matrix(i,2) = -1._db

      case('-y+x')
        Matrix(i,1) = 1._db
        Matrix(i,2) = -1._db

      case('y-x')
        Matrix(i,1) = -1._db
        Matrix(i,2) = 1._db

      case('-x+y')
        Matrix(i,1) = -1._db
        Matrix(i,2) = 1._db

      case('1/2+x')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 0.5_db

      case('x+1/2')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 0.5_db

      case('1/2-x')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 0.5_db

      case('-x+1/2')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 0.5_db

      case('1/2+y')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 0.5_db

      case('y+1/2')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 0.5_db

      case('1/2-y')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 0.5_db

      case('-y+1/2')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 0.5_db

      case('1/2+z')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 0.5_db

      case('z+1/2')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 0.5_db

      case('1/2-z')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 0.5_db

      case('-z+1/2')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 0.5_db

      case('1/4+x')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 0.25_db

      case('x+1/4')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 0.25_db

      case('1/4-x')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 0.25_db

      case('-x+1/4')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 0.25_db

      case('1/4+y')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 0.25_db

      case('y+1/4')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 0.25_db

      case('1/4-y')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 0.25_db

      case('-y+1/4')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 0.25_db

      case('1/4+z')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 0.25_db

      case('z+1/4')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 0.25_db

      case('1/4-z')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 0.25_db

      case('-z+1/4')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 0.25_db

      case('3/4+x')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 0.75_db

      case('x+3/4')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 0.75_db

      case('3/4-x')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 0.75_db

      case('-x+3/4')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 0.75_db

      case('3/4+y')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 0.75_db

      case('y+3/4')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 0.75_db

      case('3/4-y')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 0.75_db

      case('-y+3/4')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 0.75_db

      case('3/4+z')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 0.75_db

      case('z+3/4')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 0.75_db

      case('3/4-z')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 0.75_db

      case('-z+3/4')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 0.75_db

      case('z+1/6')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 1/6._db

      case('z+1/3')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 1/3._db

      case('-z+1/3')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 1/3._db

      case('z+2/3')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 2/3._db

      case('-z+2/3')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 2/3._db

      case('z+5/6')
        Matrix(i,3) = 1.
        Matrix(i,4) = 5/6._db

      case default
        call write_error
        do ipr = 6,9,3
          write(ipr,100) line
        end do
        stop

    end select

    ibegin = ifin + 2

  end do

  return
  100 format(//' Sorry, an operation is not known in the line',/, &
      1x,A,/,' Please add it to subroutine Findop in spgroup.f'/)
end

!***********************************************************************

subroutine analysename(line,sgnb,sgnbcar,sgschoenfliess, sgHMshort,sgHMlong)

! This program analyses the line containing various names of a space
! group. This line was generated with the space group program SGinfo

!  line       Line of data
!  sgnb       Space group number
!  sgnbcar    Space group number (and eventually axis choice
!             or origin choice) in characters
!  sgschoenfliess Space group name in Schoenfliess notation
!  sgHMshort  Space group name in Hermann-Mauguin short notation
!  sgHMlong   Space group name in Hermann-Mauguin long notation

  use declarations
  implicit none

  character(len=10):: sgnbcar
  character(len=13):: sgschoenfliess
  character(len=80):: line
  character(len=27):: sgHMshort,sgHMlong

  integer sgnb,i0,i

!   Read space group number
  read(line,'(1x,i3)') sgnb
  read(line,'(1x,a10)') sgnbcar
  read(line,'(12x,a13)') sgschoenfliess
  read(line,'(26x,a26)') sgHMshort

!    Find Hermann-Mauguin long name
  i0 = index(sgHMshort,'=')
  if(i0.eq.0) then
    sgHMlong = sgHMshort
  else
    do i = i0+2,26
      sgHMlong(i-i0-1:i-i0-1) = sgHMshort(i:i)
    end do
    do i = 26-i0,26
      sgHMlong(i:i) = ' '
    end do
    do i = i0,26
      sgHMshort(i:i) = ' '
    end do
  end if

!    For Hermann-Mauguin short name, strip additional characters
  i0 = index(sgHMshort,':')
  if(i0.ne.0) then
    do i = i0,len(sgHMshort)
      sgHMshort(i:i) = ' '
    end do
  end if

  return
end

!***********************************************************************

subroutine locateSG(itape,Space_file,Space_Group,sgnbcar0)

!    This program locates all space groups whose names
!    look like Space_Group and asks to choose the right one
!  line       Line of data
!  sgnb       Space group number
!  sgnbcar    Space group number (and eventually axis choice
!             or origin choice) in characters
!  sgschoenfliess Space group name in Schoenfliess notation
!  sgHMshort  Space group name in Hermann-Mauguin short notation
!  sgHMlong   Space group name in Hermann-Mauguin long notation

  use declarations
  implicit none

  character(len=10) sgnbcar, sgnbcar0, sgnbcar1
  character(len=13) sgHMlong13, sgHMshort13, sgnbcar13, sgschoenfliess, sgschoenfliess1, Space_Group
  character(len=27) sgHMshort, sgHMshort1, sgHMlong, sgHMlong1
  character(len=80) line
  character(len=132) Space_file

  integer i, ipr, istat, itape, sgnb, sgnb1, nbsol

  logical pareil

  nbsol = 0

  do i = 1,10000

    read(itape,'(A)',iostat=istat) line
    if( istat /= 0 ) exit
    if( line(1:1) /= '*' ) cycle

    call analysename(line,sgnb,sgnbcar,sgschoenfliess,sgHMshort,sgHMlong)

    sgHMlong = Adjustl( sgHMlong )
    sgHMlong13(1:13) = sgHMlong(1:13)
    sgHMshort = Adjustl( sgHMshort )
    sgHMshort13(1:13) = sgHMshort(1:13)
    sgnbcar13 = ' '
    sgnbcar13(1:10) = sgnbcar(1:10)
    sgnbcar13 = Adjustl( sgnbcar13 )

    Pareil = ( sgnbcar13 == Space_group ) .or. ( sgschoenfliess == Space_group ) .or. ( sgHMshort13 == Space_group ) &
           .or. ( sgHMlong13 == Space_group )

    if( pareil ) then
      nbsol = nbsol + 1
      sgnbcar0 = sgnbcar
      if( nbsol == 2 ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,110) Space_Group
          write(ipr,120)
          write(ipr,130) sgnb1, sgnbcar1, sgschoenfliess1, sgHMshort1, sgHMlong1
        end do
      end if
      if( nbsol == 1 ) then
        sgnb1 = sgnb; sgnbcar1 = sgnbcar
        sgschoenfliess1 = sgschoenfliess
        sgHMshort1 = sgHMshort
        sgHMlong1 = sgHMlong
      else
        do ipr = 6,9,3
          write(ipr,130) sgnb, sgnbcar, sgschoenfliess, sgHMshort, sgHMlong
        end do
      endif
    end if

  end do

  if( nbsol > 1 ) then
    do ipr = 6,9,3
      write(ipr,140) Space_file
    end do
    stop
  endif

  return
  110 format(/' Space group is ',a13)
  120 format(/' Nb  Full Nb', ' Schoenfliess Hermann-Mauguin Long Hermann-Mauguin')
  130 format(i3,2x,a11,1x,a10,3x,a10,9x,a10)
  140 format(/'   There are more than one definition of the', ' group operations !',/, &
    ' Please enter the full number (Full Nb in the above list)',/ ' of the set of operations that you desire.',/ &
    ' See the International Tables or the file ',a10,' for more', ' detail.'/)
end
