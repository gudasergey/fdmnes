! FDMNES subroutines

! Program giving all the atoms of the mesh from the non equivalent atoms and from the space group.
! A big part comes from Ch. Brouder.

! Space_Group Name of the symmetry group
! NMAXOP      Maximum number of symmetry operations

subroutine spgroup(Cif,Cif_file,Do_exp,neq,ngroup,ngroup_neq,itype,posn,posout,Space_group)

  use declarations
  implicit none

  integer, parameter:: nmaxop = 192

  integer:: i, ia, is, j, ja, js, k, ks, ngroup, ngroup_neq, nsym
  
  character(len=1):: SGTrans
  character(len=13):: Space_Group
  character(len=Length_name):: Cif_file

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

  call symgrp(Cif,Cif_file,Space_Group,Mat,Trans,nsym,nmaxop,SGTrans)

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
!        write(6,'(6i3,1p,7e13.5)') ia,is,ks,js,itypeq(ks), itype(ia), sum( abs(qq(:,ks)-q(:)) ), q(:), qq(:,ks) 
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

subroutine symgrp(Cif,Cif_file,Space_Group,Mat,Trans,nbsyop,nmaxop,SGTrans)

! This subroutine looks for the space group whose name is in Space_Group.
! If it is found, it outputs the number of symmetry operations, and builds the matrices for these operations.

! Variable description :
!   Space_Group Name of the symmetry group
!   NBSYOP      Number of symmetry operations
!   NMAXOP      Maximum number of symmetry operations
!   sgnb        Space group number

  use declarations
  implicit none

  integer, parameter:: nline_spgr = 5046
  integer:: Length, nmaxop

  character(len=1):: SGTrans
  character(len=10):: sgnbcar, sgnbcar0
  character(len=13):: sgschoenfliess, Space_Group
  character(len=27):: sgHMshort, sgHMlong
  character(len=80):: line, mot, motb
  character(len=Length_name):: Cif_file
  character(len=80), dimension(nmaxop):: lines

  character(len=80), dimension(nline_spgr):: spacegroupdata
  common /spacegroupdata/ spacegroupdata

  integer:: i, i1, i2, ipr, istat, itape, j, k, l, nbsyop, sgnb, currentLineNum

  logical:: Cif, index_sym, pareil

  real(kind=db), dimension(3,3,nmaxop):: Mat
  real(kind=db), dimension(3,nmaxop):: Trans
  real(kind=db), dimension(3,4):: Matrix(3,4)

  pareil = .false.
  index_sym = .false.

  if( Cif ) then
    itape = 7

    Open(itape, file = Cif_file, status='old', iostat=istat)
    if( istat /= 0 ) call write_open_error(Cif_file,istat,1)

    do
      read(itape,'(A)') mot
      Length = len_trim(mot)
      do i = 1,Length
       if( mot(i:i) == char(9) ) mot(i:i) = ' '
      end do
      mot = adjustl(mot)
      if( mot(1:27) == '_symmetry_equiv_pos_site_id' ) index_sym = .true.
      if( mot(1:26) == '_symmetry_equiv_pos_as_xyz' .or. mot(1:32) == '_space_group_symop_operation_xyz' ) exit 
    end do
    
    boucle_i: do i = 1,1000

      read(itape,'(A)' ) mot
      do j = 1,80
        if( mot(j:j) == char(9) ) mot(j:j) = ' '
      end do
      line = ' '

      mot = adjustl(mot)
      
      if( index_sym ) then
        do j = 1,80
          if( mot(j:j) == ' ' ) exit
          mot(j:j) = ' '
        end do
      endif
            
      do j = 1,80
        if( mot(j:j) == ',' ) exit
      end do
      if( j > 80 ) exit boucle_i

      boucle_j: do j = 1,80
        if( mot(j:j) == "'" ) then
          do k = 1,j
            mot(k:k) = ' '
          end do
          do k = j+1,80
            if( mot(k:k) == "'" ) then
              mot(k:k) = ' '
              exit boucle_j
            endif
          end do
        endif 
      end do boucle_j

      do k = 1,2
        motb = adjustl( mot )
        l = len_trim( motb )
        do j = 1,l-1
          if( motb(j:j) /= ' ' ) cycle
          mot = ' '
          mot(1:j-1) = motb(1:j-1)
          mot(j:l-1) = motb(j+1:l)  
        end do
      end do
      
      lines(i) = adjustl(mot)
        
    end do boucle_i
    
    nbsyop = i - 1

    SGTrans = ' ' ! Space_group(1:1)

    Close(itape)
  else
  
! Ask for exact definition of space group.
! sgnbcar0 is the detailed number of the spacegroup, specifying axis and origin conventions

    call locateSG(Space_Group,sgnbcar0)

! In spacegroupdata the name of the symmetry group follows a *<space> look for it.
! If a * is found, check that the following string is the name of the desired symmetry group.

    do i = 1,nline_spgr

      line = spacegroupdata(i)
      if( line(1:1) /= '*' ) cycle

! Analyse the line giving space group name(s)
      call analysename(line,sgnb,sgnbcar,sgschoenfliess,sgHMshort,sgHMlong)
      Pareil = sgnbcar == sgnbcar0
      SGTrans = sgHMlong(1:1)
      if( index(sgHMlong,'H') /= 0 ) SGTrans = 'H'
      if( pareil ) exit

    end do
    currentLineNum = i
! Look for nbsyop
    do i = 1,1000 
      j = currentLineNum + i
      if( j > nline_spgr ) exit 
      line = spacegroupdata(j)
      if( line(1:1) ==   '*' ) exit
      lines(i) = line
    end do
    nbsyop = i - 1

    if( index(sgnbcar,'R') /= 0 ) then
      call write_error
      do ipr = 6,9,3
        write(ipr,110)
      end do
      stop
    end if

  endif
   
! Read the NBSYOP symmetry operations, and build the corresponding transformation matrix.

  do i = 1,nbsyop
    call findop(lines(i),Matrix)
    do i1=1,3
      do i2=1,3
        Mat(i1,i2,i) = Matrix(i1,i2)
      end do
      Trans(i1,i) = Matrix(i1,4)
    end do
  end do

  return
  110 format(//' Rhombohedral axes are not implemented ! Please, convert to hexagonal axes'//)
end

!***********************************************************************

subroutine findop(line,matrix)

! This subroutine takes the line LINE coming from spacegroupdata
! and builds the matrix corresponding to the symmetry operation written in the line.
! The symmetry operation is written in the line as e.g. -y,x,-z

!   line    Line red from the input file
!   Matrix  Matrix of the symmetry operation (output) in the crystal axes

  use declarations
  implicit none

  integer, parameter:: nline_spgr = 5046

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

      case('x','+x')
        Matrix(i,1) = 1._db

      case('-x')
        Matrix(i,1) = -1._db

      case('y','+y')
        Matrix(i,2) = 1._db

      case('-y')
        Matrix(i,2) = -1._db

      case('z','+z')
        Matrix(i,3) = 1._db

      case('-z')
        Matrix(i,3) = -1._db

      case('x-y','+x-y','-y+x')
        Matrix(i,1) = 1._db
        Matrix(i,2) = -1._db

      case('y-x','+y-x','-x+y')
        Matrix(i,1) = -1._db
        Matrix(i,2) = 1._db

      case('1/2+x','x+1/2','+x+1/2')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 0.5_db

      case('x-1/2','+x-1/2','-1/2+x')
        Matrix(i,1) = 1._db
        Matrix(i,4) = -0.5_db

      case('1/2-x','-x+1/2')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 0.5_db

      case('-1/2-x','-x-1/2')
        Matrix(i,1) = -1._db
        Matrix(i,4) = -0.5_db

      case('1/2+y','y+1/2','+y+1/2')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 0.5_db

      case('1/2-y','-y+1/2')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 0.5_db

      case('-1/2+y','y-1/2','+y-1/2')
        Matrix(i,2) =  1._db
        Matrix(i,4) = -0.5_db

      case('-y-1/2','-1/2-y')
        Matrix(i,2) = -1._db
        Matrix(i,4) = -0.5_db

      case('z+1/2','+z+1/2','1/2+z')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 0.5_db

      case('1/2-z','-z+1/2')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 0.5_db

      case('-1/2+z','z-1/2','+z-1/2')
        Matrix(i,3) = 1._db
        Matrix(i,4) = -0.5_db

      case('-z-1/2','-1/2-z')
        Matrix(i,3) = -1._db
        Matrix(i,4) = -0.5_db

      case('1/4+x','x+1/4','+x+1/4')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 0.25_db

      case('1/4-x','-x+1/4')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 0.25_db

      case('1/4+y','y+1/4','+y+1/4')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 0.25_db

      case('1/4-y','-y+1/4')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 0.25_db

      case('1/4+z','z+1/4','+z+1/4')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 0.25_db

      case('1/4-z','+1/4-z','-z+1/4')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 0.25_db

      case('3/4+x','x+3/4','+x+3/4')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 0.75_db

      case('3/4-x','+3/4-x','-x+3/4')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 0.75_db

      case('3/4+y','+3/4+y','y+3/4','+y+3/4')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 0.75_db

      case('3/4-y','+3/4-y','-y+3/4')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 0.75_db

      case('3/4+z','+3/4+z','z+3/4','+z+3/4')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 0.75_db

      case('3/4-z','+3/4-z','-z+3/4')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 0.75_db

      case('z+1/6','+z+1/6','1/6+z')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 1/6._db

      case('-z+1/6','+1/6-z','1/6-z')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 1/6._db

      case('z+1/3','+z+1/3','1/3+z')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 1/3._db

      case('-z+1/3','1/3-z','+1/3-z')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 1/3._db

      case('-z+5/6','5/6-z','+5/6-z')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 5/6._db

      case('z+2/3','+z+2/3','2/3+z')
        Matrix(i,3) = 1._db
        Matrix(i,4) = 2/3._db

      case('-z+2/3','2/3-z')
        Matrix(i,3) = -1._db
        Matrix(i,4) = 2/3._db

      case('z+5/6','+z+5/6','5/6+z','+5/6+z')
        Matrix(i,3) = 1.
        Matrix(i,4) = 5/6._db

      case('x+1/3','+x+1/3','1/3+x','+1/3+x')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 1/3._db

      case('x+2/3','+x+2/3','2/3+x')
        Matrix(i,1) = 1._db
        Matrix(i,4) = 2/3._db

      case('-x+1/3','1/3-x')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 1/3._db

      case('-x+2/3','2/3-x')
        Matrix(i,1) = -1._db
        Matrix(i,4) = 2/3._db

      case('y+1/3','+y+1/3','1/3+y')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 1/3._db

      case('y+2/3','+y+2/3','2/3+y')
        Matrix(i,2) = 1._db
        Matrix(i,4) = 2/3._db

      case('-y+1/3','1/3-y')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 1/3._db

      case('-y+2/3','2/3-y')
        Matrix(i,2) = -1._db
        Matrix(i,4) = 2/3._db

      case('x-y+1/3','+x-y+1/3','-y+x+1/3','1/3+x-y')
        Matrix(i,1) = 1._db
        Matrix(i,2) = -1._db
        Matrix(i,4) = 1/3._db

      case('x-y+2/3','+x-y+2/3','-y+x+2/3','2/3+x-y')
        Matrix(i,1) = 1._db
        Matrix(i,2) = -1._db
        Matrix(i,4) = 2/3._db

      case('-x+y+1/3','y-x+1/3','+y-x+1/3','1/3-x+y')
        Matrix(i,1) = -1._db
        Matrix(i,2) = 1._db
        Matrix(i,4) = 1/3._db

      case('-x+y+2/3','y-x+2/3','+y-x+2/3','2/3-x+y')
        Matrix(i,1) = -1._db
        Matrix(i,2) = 1._db
        Matrix(i,4) = 2/3._db

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

subroutine analysename(line,sgnb,sgnbcar,sgschoenfliess,sgHMshort,sgHMlong)

! This program analyses the line containing various names of a space group. 
! This line was generated with the space group program SGinfo

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

subroutine locateSG(Space_Group,sgnbcar0)

!    This program locates all space groups whose names
!    look like Space_Group and asks to choose the right one
!  line       Line of data
!  sgnb       Space group number
!  sgnbcar    Space group number (and eventually axis choice or origin choice) in characters
!  sgschoenfliess Space group name in Schoenfliess notation
!  sgHMshort  Space group name in Hermann-Mauguin short notation
!  sgHMlong   Space group name in Hermann-Mauguin long notation

  use declarations
  implicit none

  integer, parameter:: nline_spgr = 5046

  character(len=80), dimension(nline_spgr):: spacegroupdata
  common /spacegroupdata/ spacegroupdata

  character(len=10) sgnbcar, sgnbcar0, sgnbcar1
  character(len=13) sgHMlong13, sgHMshort13, sgnbcar13, sgschoenfliess, sgschoenfliess1, Space_Group
  character(len=27) sgHMshort, sgHMshort1, sgHMlong, sgHMlong1
  character(len=80) line

  integer i, ipr, sgnb, sgnb1, nbsol

  logical pareil

  nbsol = 0

  do i = 1,nline_spgr

    line = spacegroupdata(i)
    if( line(1:1) /=   '*' ) cycle

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

    if( Pareil ) then
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
    endif
    
  end do

  if( nbsol == 0 ) then 
    call write_error
    do ipr = 6,9,3
      write(ipr,140) Space_group
    end do
    stop
  elseif( nbsol > 1 ) then 
    do ipr = 6,9,3
      write(ipr,150)
    end do
    stop
  endif
  
  return
  110 format(/' Space group is ',a13)
  120 format(/' Nb  Full Nb', ' Schoenfliess Hermann-Mauguin Long Hermann-Mauguin')
  130 format(i3,2x,a11,1x,a10,3x,a10,9x,a10)
  140 format(//' Space-group name: ',a13, / &
               ' not found in the fdmnes basis',// &
               ' Check:', / &
               ' 1) The spelling is correct and corresponds to the Internationnal Tables convention.',/ &
               ' 2) Some space group can be defined in hexagonal or rhombohedral unit cell and must have',/ &
               '    respectively :H or :R extention after the general name, (example R-3C:H or 167:H)',/ &
               '    (groups 146, 148, 155, 160, 161, 166, 167).',/ &
               ' 3) A number of other space-groups can be defined with different center or axis definition',/ &
               '    which must be defined with :1, :2 (example Pn-3:1) or some others extensions', // &
               ' See the International Tables for more details.', // &
               ' Without answer after checking all this, contact the author of the code',//)
  150 format(/'   There are more than one possible definition of the group operations !',/, &
              ' To define completely the space group in the indata file, chose in the list above,'/ &
              ' one of the "Full Nb" or of the "Long Hermann-Mauguin" symbol.',// &
              ' See the International Tables for more detail.'/)
end

!----------------------------------------------------------------------------------------------------------------

! Spacegroup.data
 
! From Ch. Brouder 4-dec-95

!  This is list of all space groups, with all choice of axes and all origins 
!  as defined in the Handbook of Crystallography, Vol.A.
!  For each space group, the first line gives the space group number,
!  the Schoenfliess symbol, the Hermann-Mauguin symbol (short and long),
!  and the Hall symbol (S.R. Hall; Space-Group Notation with an Explicit
!  Origin, Acta Cryst. (1981). A37, 517-525, or International Tables
!  Volume B 1994, Section 1.4. Symmetry in reciprocal space (Sydney R. Hall,
!  Crystallography Centre, University of Western Australia. (syd@crystal.uwa.edu.au).
!  This file was generated automatically using the space group analysis
!  program "sginfo", written by Ralf W. Grosse-Kunstleve, Laboratory of
!  Crystallography, ETH Zurich, Switzerland (ralf@kristall.erdw.ethz.ch).
!  In this list all space group operations are given except for the
!  translations given by the lattice type (R,I,F)
!  Group C-1 added by F. Farges.

block data spacegroupdatasubroutine

  integer, parameter:: nline_spgr = 5046

  character(len=80), dimension(nline_spgr):: spacegroupdata
  common /spacegroupdata/ spacegroupdata

 data spacegroupdata(   1:   2)/ &                              ! group index =   1, number of symmetry op =  1
  '*  1        C1^1          P1                           P 1', &                   
  'x,y,z' /                                                                                                                         

 data spacegroupdata(   3:   5)/ &                              ! group index =   2, number of symmetry op =  2
  '*  2:a      Ci^1          P-1                         -P 1', &                   
  'x,y,z','-x,-y,-z' /                                                                                                              

 data spacegroupdata(   6:  10)/ &                              ! group index =   3, number of symmetry op =  4
  '*  2:b      Ci^1          C-1                         -C 1', &                   
  'x,y,z','-x,-y,-z','x+1/2,y+1/2,z','-x+1/2,-y+1/2,-z' /                                                                           

 data spacegroupdata(  11:  13)/ &                              ! group index =   4, number of symmetry op =  2
  '*  3:b      C2^1          P2:b = P121                  P 2y', &                  
  'x,y,z','-x,y,-z' /                                                                                                               

 data spacegroupdata(  14:  16)/ &                              ! group index =   5, number of symmetry op =  2
  '*  3:c      C2^1          P2:c = P112                  P 2', &                   
  'x,y,z','-x,-y,z' /                                                                                                               

 data spacegroupdata(  17:  19)/ &                              ! group index =   6, number of symmetry op =  2
  '*  3:a      C2^1          P2:a = P211                  P 2x', &                  
  'x,y,z','x,-y,-z' /                                                                                                               

 data spacegroupdata(  20:  22)/ &                              ! group index =   7, number of symmetry op =  2
  '*  4:b      C2^2          P21:b = P1211                P 2yb', &                 
  'x,y,z','-x,y+1/2,-z' /                                                                                                           

 data spacegroupdata(  23:  25)/ &                              ! group index =   8, number of symmetry op =  2
  '*  4:c      C2^2          P21:c = P1121                P 2c', &                  
  'x,y,z','-x,-y,z+1/2' /                                                                                                           

 data spacegroupdata(  26:  28)/ &                              ! group index =   9, number of symmetry op =  2
  '*  4:a      C2^2          P21:a = P2111                P 2xa', &                 
  'x,y,z','x+1/2,-y,-z' /                                                                                                           

 data spacegroupdata(  29:  31)/ &                              ! group index =  10, number of symmetry op =  2
  '*  5:b1     C2^3          C2:b1 = C121                 C 2y', &                  
  'x,y,z','-x,y,-z' /                                                                                                               

 data spacegroupdata(  32:  34)/ &                              ! group index =  11, number of symmetry op =  2
  '*  5:b2     C2^3          C2:b2 = A121                 A 2y', &                  
  'x,y,z','-x,y,-z' /                                                                                                               

 data spacegroupdata(  35:  37)/ &                              ! group index =  12, number of symmetry op =  2
  '*  5:b3     C2^3          C2:b3 = I121                 I 2y', &                  
  'x,y,z','-x,y,-z' /                                                                                                               

 data spacegroupdata(  38:  40)/ &                              ! group index =  13, number of symmetry op =  2
  '*  5:c1     C2^3          C2:c1 = A112                 A 2', &                   
  'x,y,z','-x,-y,z' /                                                                                                               

 data spacegroupdata(  41:  43)/ &                              ! group index =  14, number of symmetry op =  2
  '*  5:c2     C2^3          C2:c2 = B112                 B 2', &                   
  'x,y,z','-x,-y,z' /                                                                                                               

 data spacegroupdata(  44:  46)/ &                              ! group index =  15, number of symmetry op =  2
  '*  5:c3     C2^3          C2:c3 = I112                 I 2', &                   
  'x,y,z','-x,-y,z' /                                                                                                               

 data spacegroupdata(  47:  49)/ &                              ! group index =  16, number of symmetry op =  2
  '*  5:a1     C2^3          C2:a1 = B211                 B 2x', &                  
  'x,y,z','x,-y,-z' /                                                                                                               

 data spacegroupdata(  50:  52)/ &                              ! group index =  17, number of symmetry op =  2
  '*  5:a2     C2^3          C2:a2 = C211                 C 2x', &                  
  'x,y,z','x,-y,-z' /                                                                                                               

 data spacegroupdata(  53:  55)/ &                              ! group index =  18, number of symmetry op =  2
  '*  5:a3     C2^3          C2:a3 = I211                 I 2x', &                  
  'x,y,z','x,-y,-z' /                                                                                                               

 data spacegroupdata(  56:  58)/ &                              ! group index =  19, number of symmetry op =  2
  '*  6:b      Cs^1          Pm:b = P1m1                  P -2y', &                 
  'x,y,z','x,-y,z' /                                                                                                                

 data spacegroupdata(  59:  61)/ &                              ! group index =  20, number of symmetry op =  2
  '*  6:c      Cs^1          Pm:c = P11m                  P -2', &                  
  'x,y,z','x,y,-z' /                                                                                                                

 data spacegroupdata(  62:  64)/ &                              ! group index =  21, number of symmetry op =  2
  '*  6:a      Cs^1          Pm:a = Pm11                  P -2x', &                 
  'x,y,z','-x,y,z' /                                                                                                                

 data spacegroupdata(  65:  67)/ &                              ! group index =  22, number of symmetry op =  2
  '*  7:b1     Cs^2          Pc:b1 = P1c1                 P -2yc', &                
  'x,y,z','x,-y,z+1/2' /                                                                                                            

 data spacegroupdata(  68:  70)/ &                              ! group index =  23, number of symmetry op =  2
  '*  7:b2     Cs^2          Pc:b2 = P1n1                 P -2yac', &               
  'x,y,z','x+1/2,-y,z+1/2' /                                                                                                        

 data spacegroupdata(  71:  73)/ &                              ! group index =  24, number of symmetry op =  2
  '*  7:b3     Cs^2          Pc:b3 = P1a1                 P -2ya', &                
  'x,y,z','x+1/2,-y,z' /                                                                                                            

 data spacegroupdata(  74:  76)/ &                              ! group index =  25, number of symmetry op =  2
  '*  7:c1     Cs^2          Pc:c1 = P11a                 P -2a', &                 
  'x,y,z','x+1/2,y,-z' /                                                                                                            

 data spacegroupdata(  77:  79)/ &                              ! group index =  26, number of symmetry op =  2
  '*  7:c2     Cs^2          Pc:c2 = P11n                 P -2ab', &                
  'x,y,z','x+1/2,y+1/2,-z' /                                                                                                        

 data spacegroupdata(  80:  82)/ &                              ! group index =  27, number of symmetry op =  2
  '*  7:c3     Cs^2          Pc:c3 = P11b                 P -2b', &                 
  'x,y,z','x,y+1/2,-z' /                                                                                                            

 data spacegroupdata(  83:  85)/ &                              ! group index =  28, number of symmetry op =  2
  '*  7:a1     Cs^2          Pc:a1 = Pb11                 P -2xb', &                
  'x,y,z','-x,y+1/2,z' /                                                                                                            

 data spacegroupdata(  86:  88)/ &                              ! group index =  29, number of symmetry op =  2
  '*  7:a2     Cs^2          Pc:a2 = Pn11                 P -2xbc', &               
  'x,y,z','-x,y+1/2,z+1/2' /                                                                                                        

 data spacegroupdata(  89:  91)/ &                              ! group index =  30, number of symmetry op =  2
  '*  7:a3     Cs^2          Pc:a3 = Pc11                 P -2xc', &                
  'x,y,z','-x,y,z+1/2' /                                                                                                            

 data spacegroupdata(  92:  94)/ &                              ! group index =  31, number of symmetry op =  2
  '*  8:b1     Cs^3          Cm:b1 = C1m1                 C -2y', &                 
  'x,y,z','x,-y,z' /                                                                                                                

 data spacegroupdata(  95:  97)/ &                              ! group index =  32, number of symmetry op =  2
  '*  8:b2     Cs^3          Cm:b2 = A1m1                 A -2y', &                 
  'x,y,z','x,-y,z' /                                                                                                                

 data spacegroupdata(  98: 100)/ &                              ! group index =  33, number of symmetry op =  2
  '*  8:b3     Cs^3          Cm:b3 = I1m1                 I -2y', &                 
  'x,y,z','x,-y,z' /                                                                                                                

 data spacegroupdata( 101: 103)/ &                              ! group index =  34, number of symmetry op =  2
  '*  8:c1     Cs^3          Cm:c1 = A11m                 A -2', &                  
  'x,y,z','x,y,-z' /                                                                                                                

 data spacegroupdata( 104: 106)/ &                              ! group index =  35, number of symmetry op =  2
  '*  8:c2     Cs^3          Cm:c2 = B11m                 B -2', &                  
  'x,y,z','x,y,-z' /                                                                                                                

 data spacegroupdata( 107: 109)/ &                              ! group index =  36, number of symmetry op =  2
  '*  8:c3     Cs^3          Cm:c3 = I11m                 I -2', &                  
  'x,y,z','x,y,-z' /                                                                                                                

 data spacegroupdata( 110: 112)/ &                              ! group index =  37, number of symmetry op =  2
  '*  8:a1     Cs^3          Cm:a1 = Bm11                 B -2x', &                 
  'x,y,z','-x,y,z' /                                                                                                                

 data spacegroupdata( 113: 115)/ &                              ! group index =  38, number of symmetry op =  2
  '*  8:a2     Cs^3          Cm:a2 = Cm11                 C -2x', &                 
  'x,y,z','-x,y,z' /                                                                                                                

 data spacegroupdata( 116: 118)/ &                              ! group index =  39, number of symmetry op =  2
  '*  8:a3     Cs^3          Cm:a3 = Im11                 I -2x', &                 
  'x,y,z','-x,y,z' /                                                                                                                

 data spacegroupdata( 119: 121)/ &                              ! group index =  40, number of symmetry op =  2
  '*  9:b1     Cs^4          Cc:b1 = C1c1                 C -2yc', &                
  'x,y,z','x,-y,z+1/2' /                                                                                                            

 data spacegroupdata( 122: 124)/ &                              ! group index =  41, number of symmetry op =  2
  '*  9:b2     Cs^4          Cc:b2 = A1n1                 A -2yac', &               
  'x,y,z','x+1/2,-y,z+1/2' /                                                                                                        

 data spacegroupdata( 125: 127)/ &                              ! group index =  42, number of symmetry op =  2
  '*  9:b3     Cs^4          Cc:b3 = I1a1                 I -2ya', &                
  'x,y,z','x+1/2,-y,z' /                                                                                                            

 data spacegroupdata( 128: 130)/ &                              ! group index =  43, number of symmetry op =  2
  '*  9:-b1    Cs^4          Cc:-b1 = A1a1                A -2ya', &                
  'x,y,z','x+1/2,-y,z' /                                                                                                            

 data spacegroupdata( 131: 133)/ &                              ! group index =  44, number of symmetry op =  2
  '*  9:-b2    Cs^4          Cc:-b2 = C1n1                C -2ybc', &               
  'x,y,z','x,-y+1/2,z+1/2' /                                                                                                        

 data spacegroupdata( 134: 136)/ &                              ! group index =  45, number of symmetry op =  2
  '*  9:-b3    Cs^4          Cc:-b3 = I1c1                I -2yc', &                
  'x,y,z','x,-y,z+1/2' /                                                                                                            

 data spacegroupdata( 137: 139)/ &                              ! group index =  46, number of symmetry op =  2
  '*  9:c1     Cs^4          Cc:c1 = A11a                 A -2a', &                 
  'x,y,z','x+1/2,y,-z' /                                                                                                            

 data spacegroupdata( 140: 142)/ &                              ! group index =  47, number of symmetry op =  2
  '*  9:c2     Cs^4          Cc:c2 = B11n                 B -2bc', &                
  'x,y,z','x,y+1/2,-z+1/2' /                                                                                                        

 data spacegroupdata( 143: 145)/ &                              ! group index =  48, number of symmetry op =  2
  '*  9:c3     Cs^4          Cc:c3 = I11b                 I -2b', &                 
  'x,y,z','x,y+1/2,-z' /                                                                                                            

 data spacegroupdata( 146: 148)/ &                              ! group index =  49, number of symmetry op =  2
  '*  9:-c1    Cs^4          Cc:-c1 = B11b                B -2b', &                 
  'x,y,z','x,y+1/2,-z' /                                                                                                            

 data spacegroupdata( 149: 151)/ &                              ! group index =  50, number of symmetry op =  2
  '*  9:-c2    Cs^4          Cc:-c2 = A11n                A -2ac', &                
  'x,y,z','x+1/2,y,-z+1/2' /                                                                                                        

 data spacegroupdata( 152: 154)/ &                              ! group index =  51, number of symmetry op =  2
  '*  9:-c3    Cs^4          Cc:-c3 = I11a                I -2a', &                 
  'x,y,z','x+1/2,y,-z' /                                                                                                            

 data spacegroupdata( 155: 157)/ &                              ! group index =  52, number of symmetry op =  2
  '*  9:a1     Cs^4          Cc:a1 = Bb11                 B -2xb', &                
  'x,y,z','-x,y+1/2,z' /                                                                                                            

 data spacegroupdata( 158: 160)/ &                              ! group index =  53, number of symmetry op =  2
  '*  9:a2     Cs^4          Cc:a2 = Cn11                 C -2xbc', &               
  'x,y,z','-x,y+1/2,z+1/2' /                                                                                                        

 data spacegroupdata( 161: 163)/ &                              ! group index =  54, number of symmetry op =  2
  '*  9:a3     Cs^4          Cc:a3 = Ic11                 I -2xc', &                
  'x,y,z','-x,y,z+1/2' /                                                                                                            

 data spacegroupdata( 164: 166)/ &                              ! group index =  55, number of symmetry op =  2
  '*  9:-a1    Cs^4          Cc:-a1 = Cc11                C -2xc', &                
  'x,y,z','-x,y,z+1/2' /                                                                                                            

 data spacegroupdata( 167: 169)/ &                              ! group index =  56, number of symmetry op =  2
  '*  9:-a2    Cs^4          Cc:-a2 = Bn11                B -2xbc', &               
  'x,y,z','-x,y+1/2,z+1/2' /                                                                                                        

 data spacegroupdata( 170: 172)/ &                              ! group index =  57, number of symmetry op =  2
  '*  9:-a3    Cs^4          Cc:-a3 = Ib11                I -2xb', &                
  'x,y,z','-x,y+1/2,z' /                                                                                                            

 data spacegroupdata( 173: 177)/ &                              ! group index =  58, number of symmetry op =  4
  '* 10:b      C2h^1         P2/m:b = P12/m1             -P 2y', &                  
  'x,y,z','-x,y,-z','-x,-y,-z','x,-y,z' /                                                                                           

 data spacegroupdata( 178: 182)/ &                              ! group index =  59, number of symmetry op =  4
  '* 10:c      C2h^1         P2/m:c = P112/m             -P 2', &                   
  'x,y,z','-x,-y,z','-x,-y,-z','x,y,-z' /                                                                                           

 data spacegroupdata( 183: 187)/ &                              ! group index =  60, number of symmetry op =  4
  '* 10:a      C2h^1         P2/m:a = P2/m11             -P 2x', &                  
  'x,y,z','x,-y,-z','-x,-y,-z','-x,y,z' /                                                                                           

 data spacegroupdata( 188: 192)/ &                              ! group index =  61, number of symmetry op =  4
  '* 11:b      C2h^2         P21/m:b = P121/m1           -P 2yb', &                 
  'x,y,z','-x,y+1/2,-z','-x,-y,-z','x,-y+1/2,z' /                                                                                   

 data spacegroupdata( 193: 197)/ &                              ! group index =  62, number of symmetry op =  4
  '* 11:c      C2h^2         P21/m:c = P1121/m           -P 2c', &                  
  'x,y,z','-x,-y,z+1/2','-x,-y,-z','x,y,-z+1/2' /                                                                                   

 data spacegroupdata( 198: 202)/ &                              ! group index =  63, number of symmetry op =  4
  '* 11:a      C2h^2         P21/m:a = P21/m11           -P 2xa', &                 
  'x,y,z','x+1/2,-y,-z','-x,-y,-z','-x+1/2,y,z' /                                                                                   

 data spacegroupdata( 203: 207)/ &                              ! group index =  64, number of symmetry op =  4
  '* 12:b1     C2h^3         C2/m:b1 = C12/m1            -C 2y', &                  
  'x,y,z','-x,y,-z','-x,-y,-z','x,-y,z' /                                                                                           

 data spacegroupdata( 208: 212)/ &                              ! group index =  65, number of symmetry op =  4
  '* 12:b2     C2h^3         C2/m:b2 = A12/m1            -A 2y', &                  
  'x,y,z','-x,y,-z','-x,-y,-z','x,-y,z' /                                                                                           

 data spacegroupdata( 213: 217)/ &                              ! group index =  66, number of symmetry op =  4
  '* 12:b3     C2h^3         C2/m:b3 = I12/m1            -I 2y', &                  
  'x,y,z','-x,y,-z','-x,-y,-z','x,-y,z' /                                                                                           

 data spacegroupdata( 218: 222)/ &                              ! group index =  67, number of symmetry op =  4
  '* 12:c1     C2h^3         C2/m:c1 = A112/m            -A 2', &                   
  'x,y,z','-x,-y,z','-x,-y,-z','x,y,-z' /                                                                                           

 data spacegroupdata( 223: 227)/ &                              ! group index =  68, number of symmetry op =  4
  '* 12:c2     C2h^3         C2/m:c2 = B112/m            -B 2', &                   
  'x,y,z','-x,-y,z','-x,-y,-z','x,y,-z' /                                                                                           

 data spacegroupdata( 228: 232)/ &                              ! group index =  69, number of symmetry op =  4
  '* 12:c3     C2h^3         C2/m:c3 = I112/m            -I 2', &                   
  'x,y,z','-x,-y,z','-x,-y,-z','x,y,-z' /                                                                                           

 data spacegroupdata( 233: 237)/ &                              ! group index =  70, number of symmetry op =  4
  '* 12:a1     C2h^3         C2/m:a1 = B2/m11            -B 2x', &                  
  'x,y,z','x,-y,-z','-x,-y,-z','-x,y,z' /                                                                                           

 data spacegroupdata( 238: 242)/ &                              ! group index =  71, number of symmetry op =  4
  '* 12:a2     C2h^3         C2/m:a2 = C2/m11            -C 2x', &                  
  'x,y,z','x,-y,-z','-x,-y,-z','-x,y,z' /                                                                                           

 data spacegroupdata( 243: 247)/ &                              ! group index =  72, number of symmetry op =  4
  '* 12:a3     C2h^3         C2/m:a3 = I2/m11            -I 2x', &                  
  'x,y,z','x,-y,-z','-x,-y,-z','-x,y,z' /                                                                                           

 data spacegroupdata( 248: 252)/ &                              ! group index =  73, number of symmetry op =  4
  '* 13:b1     C2h^4         P2/c:b1 = P12/c1            -P 2yc', &                 
  'x,y,z','-x,y,-z+1/2','-x,-y,-z','x,-y,z+1/2' /                                                                                   

 data spacegroupdata( 253: 257)/ &                              ! group index =  74, number of symmetry op =  4
  '* 13:b2     C2h^4         P2/c:b2 = P12/n1            -P 2yac', &                
  'x,y,z','-x+1/2,y,-z+1/2','-x,-y,-z','x+1/2,-y,z+1/2' /                                                                           

 data spacegroupdata( 258: 262)/ &                              ! group index =  75, number of symmetry op =  4
  '* 13:b3     C2h^4         P2/c:b3 = P12/a1            -P 2ya', &                 
  'x,y,z','-x+1/2,y,-z','-x,-y,-z','x+1/2,-y,z' /                                                                                   

 data spacegroupdata( 263: 267)/ &                              ! group index =  76, number of symmetry op =  4
  '* 13:c1     C2h^4         P2/c:c1 = P112/a            -P 2a', &                  
  'x,y,z','-x+1/2,-y,z','-x,-y,-z','x+1/2,y,-z' /                                                                                   

 data spacegroupdata( 268: 272)/ &                              ! group index =  77, number of symmetry op =  4
  '* 13:c2     C2h^4         P2/c:c2 = P112/n            -P 2ab', &                 
  'x,y,z','-x+1/2,-y+1/2,z','-x,-y,-z','x+1/2,y+1/2,-z' /                                                                           

 data spacegroupdata( 273: 277)/ &                              ! group index =  78, number of symmetry op =  4
  '* 13:c3     C2h^4         P2/c:c3 = P112/b            -P 2b', &                  
  'x,y,z','-x,-y+1/2,z','-x,-y,-z','x,y+1/2,-z' /                                                                                   

 data spacegroupdata( 278: 282)/ &                              ! group index =  79, number of symmetry op =  4
  '* 13:a1     C2h^4         P2/c:a1 = P2/b11            -P 2xb', &                 
  'x,y,z','x,-y+1/2,-z','-x,-y,-z','-x,y+1/2,z' /                                                                                   

 data spacegroupdata( 283: 287)/ &                              ! group index =  80, number of symmetry op =  4
  '* 13:a2     C2h^4         P2/c:a2 = P2/n11            -P 2xbc', &                
  'x,y,z','x,-y+1/2,-z+1/2','-x,-y,-z','-x,y+1/2,z+1/2' /                                                                           

 data spacegroupdata( 288: 292)/ &                              ! group index =  81, number of symmetry op =  4
  '* 13:a3     C2h^4         P2/c:a3 = P2/c11            -P 2xc', &                 
  'x,y,z','x,-y,-z+1/2','-x,-y,-z','-x,y,z+1/2' /                                                                                   

 data spacegroupdata( 293: 297)/ &                              ! group index =  82, number of symmetry op =  4
  '* 14:b1     C2h^5         P21/c:b1 = P121/c1          -P 2ybc', &                
  'x,y,z','-x,y+1/2,-z+1/2','-x,-y,-z','x,-y+1/2,z+1/2' /                                                                           

 data spacegroupdata( 298: 302)/ &                              ! group index =  83, number of symmetry op =  4
  '* 14:b2     C2h^5         P21/c:b2 = P121/n1          -P 2yn', &                 
  'x,y,z','-x+1/2,y+1/2,-z+1/2','-x,-y,-z','x+1/2,-y+1/2,z+1/2' /                                                                   

 data spacegroupdata( 303: 307)/ &                              ! group index =  84, number of symmetry op =  4
  '* 14:b3     C2h^5         P21/c:b3 = P121/a1          -P 2yab', &                
  'x,y,z','-x+1/2,y+1/2,-z','-x,-y,-z','x+1/2,-y+1/2,z' /                                                                           

 data spacegroupdata( 308: 312)/ &                              ! group index =  85, number of symmetry op =  4
  '* 14:c1     C2h^5         P21/c:c1 = P1121/a          -P 2ac', &                 
  'x,y,z','-x+1/2,-y,z+1/2','-x,-y,-z','x+1/2,y,-z+1/2' /                                                                           

 data spacegroupdata( 313: 317)/ &                              ! group index =  86, number of symmetry op =  4
  '* 14:c2     C2h^5         P21/c:c2 = P1121/n          -P 2n', &                  
  'x,y,z','-x+1/2,-y+1/2,z+1/2','-x,-y,-z','x+1/2,y+1/2,-z+1/2' /                                                                   

 data spacegroupdata( 318: 322)/ &                              ! group index =  87, number of symmetry op =  4
  '* 14:c3     C2h^5         P21/c:c3 = P1121/b          -P 2bc', &                 
  'x,y,z','-x,-y+1/2,z+1/2','-x,-y,-z','x,y+1/2,-z+1/2' /                                                                           

 data spacegroupdata( 323: 327)/ &                              ! group index =  88, number of symmetry op =  4
  '* 14:a1     C2h^5         P21/c:a1 = P21/b11          -P 2xab', &                
  'x,y,z','x+1/2,-y+1/2,-z','-x,-y,-z','-x+1/2,y+1/2,z' /                                                                           

 data spacegroupdata( 328: 332)/ &                              ! group index =  89, number of symmetry op =  4
  '* 14:a2     C2h^5         P21/c:a2 = P21/n11          -P 2xn', &                 
  'x,y,z','x+1/2,-y+1/2,-z+1/2','-x,-y,-z','-x+1/2,y+1/2,z+1/2' /                                                                   

 data spacegroupdata( 333: 337)/ &                              ! group index =  90, number of symmetry op =  4
  '* 14:a3     C2h^5         P21/c:a3 = P21/c11          -P 2xac', &                
  'x,y,z','x+1/2,-y,-z+1/2','-x,-y,-z','-x+1/2,y,z+1/2' /                                                                           

 data spacegroupdata( 338: 342)/ &                              ! group index =  91, number of symmetry op =  4
  '* 15:b1     C2h^6         C2/c:b1 = C12/c1            -C 2yc', &                 
  'x,y,z','-x,y,-z+1/2','-x,-y,-z','x,-y,z+1/2' /                                                                                   

 data spacegroupdata( 343: 347)/ &                              ! group index =  92, number of symmetry op =  4
  '* 15:b2     C2h^6         C2/c:b2 = A12/n1            -A 2yac', &                
  'x,y,z','-x+1/2,y,-z+1/2','-x,-y,-z','x+1/2,-y,z+1/2' /                                                                           

 data spacegroupdata( 348: 352)/ &                              ! group index =  93, number of symmetry op =  4
  '* 15:b3     C2h^6         C2/c:b3 = I12/a1            -I 2ya', &                 
  'x,y,z','-x+1/2,y,-z','-x,-y,-z','x+1/2,-y,z' /                                                                                   

 data spacegroupdata( 353: 357)/ &                              ! group index =  94, number of symmetry op =  4
  '* 15:-b1    C2h^6         C2/c:-b1 = A12/a1           -A 2ya', &                 
  'x,y,z','-x+1/2,y,-z','-x,-y,-z','x+1/2,-y,z' /                                                                                   

 data spacegroupdata( 358: 362)/ &                              ! group index =  95, number of symmetry op =  4
  '* 15:-b2    C2h^6         C2/c:-b2 = C12/n1           -C 2ybc', &                
  'x,y,z','-x,y+1/2,-z+1/2','-x,-y,-z','x,-y+1/2,z+1/2' /                                                                           

 data spacegroupdata( 363: 367)/ &                              ! group index =  96, number of symmetry op =  4
  '* 15:-b3    C2h^6         C2/c:-b3 = I12/c1           -I 2yc', &                 
  'x,y,z','-x,y,-z+1/2','-x,-y,-z','x,-y,z+1/2' /                                                                                   

 data spacegroupdata( 368: 372)/ &                              ! group index =  97, number of symmetry op =  4
  '* 15:c1     C2h^6         C2/c:c1 = A112/a            -A 2a', &                  
  'x,y,z','-x+1/2,-y,z','-x,-y,-z','x+1/2,y,-z' /                                                                                   

 data spacegroupdata( 373: 377)/ &                              ! group index =  98, number of symmetry op =  4
  '* 15:c2     C2h^6         C2/c:c2 = B112/n            -B 2bc', &                 
  'x,y,z','-x,-y+1/2,z+1/2','-x,-y,-z','x,y+1/2,-z+1/2' /                                                                           

 data spacegroupdata( 378: 382)/ &                              ! group index =  99, number of symmetry op =  4
  '* 15:c3     C2h^6         C2/c:c3 = I112/b            -I 2b', &                  
  'x,y,z','-x,-y+1/2,z','-x,-y,-z','x,y+1/2,-z' /                                                                                   

 data spacegroupdata( 383: 387)/ &                              ! group index = 100, number of symmetry op =  4
  '* 15:-c1    C2h^6         C2/c:-c1 = B112/b           -B 2b', &                  
  'x,y,z','-x,-y+1/2,z','-x,-y,-z','x,y+1/2,-z' /                                                                                   

 data spacegroupdata( 388: 392)/ &                              ! group index = 101, number of symmetry op =  4
  '* 15:-c2    C2h^6         C2/c:-c2 = A112/n           -A 2ac', &                 
  'x,y,z','-x+1/2,-y,z+1/2','-x,-y,-z','x+1/2,y,-z+1/2' /                                                                           

 data spacegroupdata( 393: 397)/ &                              ! group index = 102, number of symmetry op =  4
  '* 15:-c3    C2h^6         C2/c:-c3 = I112/a           -I 2a', &                  
  'x,y,z','-x+1/2,-y,z','-x,-y,-z','x+1/2,y,-z' /                                                                                   

 data spacegroupdata( 398: 402)/ &                              ! group index = 103, number of symmetry op =  4
  '* 15:a1     C2h^6         C2/c:a1 = B2/b11            -B 2xb', &                 
  'x,y,z','x,-y+1/2,-z','-x,-y,-z','-x,y+1/2,z' /                                                                                   

 data spacegroupdata( 403: 407)/ &                              ! group index = 104, number of symmetry op =  4
  '* 15:a2     C2h^6         C2/c:a2 = C2/n11            -C 2xbc', &                
  'x,y,z','x,-y+1/2,-z+1/2','-x,-y,-z','-x,y+1/2,z+1/2' /                                                                           

 data spacegroupdata( 408: 412)/ &                              ! group index = 105, number of symmetry op =  4
  '* 15:a3     C2h^6         C2/c:a3 = I2/c11            -I 2xc', &                 
  'x,y,z','x,-y,-z+1/2','-x,-y,-z','-x,y,z+1/2' /                                                                                   

 data spacegroupdata( 413: 417)/ &                              ! group index = 106, number of symmetry op =  4
  '* 15:-a1    C2h^6         C2/c:-a1 = C2/c11           -C 2xc', &                 
  'x,y,z','x,-y,-z+1/2','-x,-y,-z','-x,y,z+1/2' /                                                                                   

 data spacegroupdata( 418: 422)/ &                              ! group index = 107, number of symmetry op =  4
  '* 15:-a2    C2h^6         C2/c:-a2 = B2/n11           -B 2xbc', &                
  'x,y,z','x,-y+1/2,-z+1/2','-x,-y,-z','-x,y+1/2,z+1/2' /                                                                           

 data spacegroupdata( 423: 427)/ &                              ! group index = 108, number of symmetry op =  4
  '* 15:-a3    C2h^6         C2/c:-a3 = I2/b11           -I 2xb', &                 
  'x,y,z','x,-y+1/2,-z','-x,-y,-z','-x,y+1/2,z' /                                                                                   

 data spacegroupdata( 428: 432)/ &                              ! group index = 109, number of symmetry op =  4
  '* 16        D2^1          P222                         P 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z' /                                                                                           

 data spacegroupdata( 433: 437)/ &                              ! group index = 110, number of symmetry op =  4
  '* 17        D2^2          P2221                        P 2c 2', &                
  'x,y,z','-x,-y,z+1/2','x,-y,-z','-x,y,-z+1/2' /                                                                                   

 data spacegroupdata( 438: 442)/ &                              ! group index = 111, number of symmetry op =  4
  '* 17:cab    D2^2          P2122                        P 2a 2a', &               
  'x,y,z','-x+1/2,-y,z','x+1/2,-y,-z','-x,y,-z' /                                                                                   

 data spacegroupdata( 443: 447)/ &                              ! group index = 112, number of symmetry op =  4
  '* 17:bca    D2^2          P2212                        P 2 2b', &                
  'x,y,z','-x,-y,z','x,-y+1/2,-z','-x,y+1/2,-z' /                                                                                   

 data spacegroupdata( 448: 452)/ &                              ! group index = 113, number of symmetry op =  4
  '* 18        D2^3          P21212                       P 2 2ab', &               
  'x,y,z','-x,-y,z','x+1/2,-y+1/2,-z','-x+1/2,y+1/2,-z' /                                                                           

 data spacegroupdata( 453: 457)/ &                              ! group index = 114, number of symmetry op =  4
  '* 18:cab    D2^3          P22121                       P 2bc 2', &               
  'x,y,z','-x,-y+1/2,z+1/2','x,-y,-z','-x,y+1/2,-z+1/2' /                                                                           

 data spacegroupdata( 458: 462)/ &                              ! group index = 115, number of symmetry op =  4
  '* 18:bca    D2^3          P21221                       P 2ac 2ac', &             
  'x,y,z','-x+1/2,-y,z+1/2','x+1/2,-y,-z+1/2','-x,y,-z' /                                                                           

 data spacegroupdata( 463: 467)/ &                              ! group index = 116, number of symmetry op =  4
  '* 19        D2^4          P212121                      P 2ac 2ab', &             
  'x,y,z','-x+1/2,-y,z+1/2','x+1/2,-y+1/2,-z','-x,y+1/2,-z+1/2' /                                                                   

 data spacegroupdata( 468: 472)/ &                              ! group index = 117, number of symmetry op =  4
  '* 20        D2^5          C2221                        C 2c 2', &                
  'x,y,z','-x,-y,z+1/2','x,-y,-z','-x,y,-z+1/2' /                                                                                   

 data spacegroupdata( 473: 477)/ &                              ! group index = 118, number of symmetry op =  4
  '* 20:cab    D2^5          A2122                        A 2a 2a', &               
  'x,y,z','-x+1/2,-y,z','x+1/2,-y,-z','-x,y,-z' /                                                                                   

 data spacegroupdata( 478: 482)/ &                              ! group index = 119, number of symmetry op =  4
  '* 20:bca    D2^5          B2212                        B 2 2b', &                
  'x,y,z','-x,-y,z','x,-y+1/2,-z','-x,y+1/2,-z' /                                                                                   

 data spacegroupdata( 483: 487)/ &                              ! group index = 120, number of symmetry op =  4
  '* 21        D2^6          C222                         C 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z' /                                                                                           

 data spacegroupdata( 488: 492)/ &                              ! group index = 121, number of symmetry op =  4
  '* 21:cab    D2^6          A222                         A 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z' /                                                                                           

 data spacegroupdata( 493: 497)/ &                              ! group index = 122, number of symmetry op =  4
  '* 21:bca    D2^6          B222                         B 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z' /                                                                                           

 data spacegroupdata( 498: 502)/ &                              ! group index = 123, number of symmetry op =  4
  '* 22        D2^7          F222                         F 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z' /                                                                                           

 data spacegroupdata( 503: 507)/ &                              ! group index = 124, number of symmetry op =  4
  '* 23        D2^8          I222                         I 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z' /                                                                                           

 data spacegroupdata( 508: 512)/ &                              ! group index = 125, number of symmetry op =  4
  '* 24        D2^9          I212121                      I 2b 2c', &               
  'x,y,z','-x,-y+1/2,z','x,-y,-z+1/2','-x+1/2,y,-z' /                                                                               

 data spacegroupdata( 513: 517)/ &                              ! group index = 126, number of symmetry op =  4
  '* 25        C2v^1         Pmm2                         P 2 -2', &                
  'x,y,z','-x,-y,z','-x,y,z','x,-y,z' /                                                                                             

 data spacegroupdata( 518: 522)/ &                              ! group index = 127, number of symmetry op =  4
  '* 25:cab    C2v^1         P2mm                         P -2 2', &                
  'x,y,z','x,-y,-z','x,y,-z','x,-y,z' /                                                                                             

 data spacegroupdata( 523: 527)/ &                              ! group index = 128, number of symmetry op =  4
  '* 25:bca    C2v^1         Pm2m                         P -2 -2', &               
  'x,y,z','-x,y,-z','x,y,-z','-x,y,z' /                                                                                             

 data spacegroupdata( 528: 532)/ &                              ! group index = 129, number of symmetry op =  4
  '* 26        C2v^2         Pmc21                        P 2c -2', &               
  'x,y,z','-x,-y,z+1/2','-x,y,z','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 533: 537)/ &                              ! group index = 130, number of symmetry op =  4
  '* 26:ba-c   C2v^2         Pcm21                        P 2c -2c', &              
  'x,y,z','-x,-y,z+1/2','-x,y,z+1/2','x,-y,z' /                                                                                     

 data spacegroupdata( 538: 542)/ &                              ! group index = 131, number of symmetry op =  4
  '* 26:cab    C2v^2         P21ma                        P -2a 2a', &              
  'x,y,z','x+1/2,-y,-z','x+1/2,y,-z','x,-y,z' /                                                                                     

 data spacegroupdata( 543: 547)/ &                              ! group index = 132, number of symmetry op =  4
  '* 26:-cba   C2v^2         P21am                        P -2 2a', &               
  'x,y,z','x+1/2,-y,-z','x,y,-z','x+1/2,-y,z' /                                                                                     

 data spacegroupdata( 548: 552)/ &                              ! group index = 133, number of symmetry op =  4
  '* 26:bca    C2v^2         Pb21m                        P -2 -2b', &              
  'x,y,z','-x,y+1/2,-z','x,y,-z','-x,y+1/2,z' /                                                                                     

 data spacegroupdata( 553: 557)/ &                              ! group index = 134, number of symmetry op =  4
  '* 26:a-cb   C2v^2         Pm21b                        P -2b -2', &              
  'x,y,z','-x,y+1/2,-z','x,y+1/2,-z','-x,y,z' /                                                                                     

 data spacegroupdata( 558: 562)/ &                              ! group index = 135, number of symmetry op =  4
  '* 27        C2v^3         Pcc2                         P 2 -2c', &               
  'x,y,z','-x,-y,z','-x,y,z+1/2','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 563: 567)/ &                              ! group index = 136, number of symmetry op =  4
  '* 27:cab    C2v^3         P2aa                         P -2a 2', &               
  'x,y,z','x,-y,-z','x+1/2,y,-z','x+1/2,-y,z' /                                                                                     

 data spacegroupdata( 568: 572)/ &                              ! group index = 137, number of symmetry op =  4
  '* 27:bca    C2v^3         Pb2b                         P -2b -2b', &             
  'x,y,z','-x,y,-z','x,y+1/2,-z','-x,y+1/2,z' /                                                                                     

 data spacegroupdata( 573: 577)/ &                              ! group index = 138, number of symmetry op =  4
  '* 28        C2v^4         Pma2                         P 2 -2a', &               
  'x,y,z','-x,-y,z','-x+1/2,y,z','x+1/2,-y,z' /                                                                                     

 data spacegroupdata( 578: 582)/ &                              ! group index = 139, number of symmetry op =  4
  '* 28:ba-c   C2v^4         Pbm2                         P 2 -2b', &               
  'x,y,z','-x,-y,z','-x,y+1/2,z','x,-y+1/2,z' /                                                                                     

 data spacegroupdata( 583: 587)/ &                              ! group index = 140, number of symmetry op =  4
  '* 28:cab    C2v^4         P2mb                         P -2b 2', &               
  'x,y,z','x,-y,-z','x,y+1/2,-z','x,-y+1/2,z' /                                                                                     

 data spacegroupdata( 588: 592)/ &                              ! group index = 141, number of symmetry op =  4
  '* 28:-cba   C2v^4         P2cm                         P -2c 2', &               
  'x,y,z','x,-y,-z','x,y,-z+1/2','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 593: 597)/ &                              ! group index = 142, number of symmetry op =  4
  '* 28:bca    C2v^4         Pc2m                         P -2c -2c', &             
  'x,y,z','-x,y,-z','x,y,-z+1/2','-x,y,z+1/2' /                                                                                     

 data spacegroupdata( 598: 602)/ &                              ! group index = 143, number of symmetry op =  4
  '* 28:a-cb   C2v^4         Pm2a                         P -2a -2a', &             
  'x,y,z','-x,y,-z','x+1/2,y,-z','-x+1/2,y,z' /                                                                                     

 data spacegroupdata( 603: 607)/ &                              ! group index = 144, number of symmetry op =  4
  '* 29        C2v^5         Pca21                        P 2c -2ac', &             
  'x,y,z','-x,-y,z+1/2','-x+1/2,y,z+1/2','x+1/2,-y,z' /                                                                             

 data spacegroupdata( 608: 612)/ &                              ! group index = 145, number of symmetry op =  4
  '* 29:ba-c   C2v^5         Pbc21                        P 2c -2b', &              
  'x,y,z','-x,-y,z+1/2','-x,y+1/2,z','x,-y+1/2,z+1/2' /                                                                             

 data spacegroupdata( 613: 617)/ &                              ! group index = 146, number of symmetry op =  4
  '* 29:cab    C2v^5         P21ab                        P -2b 2a', &              
  'x,y,z','x+1/2,-y,-z','x,y+1/2,-z','x+1/2,-y+1/2,z' /                                                                             

 data spacegroupdata( 618: 622)/ &                              ! group index = 147, number of symmetry op =  4
  '* 29:-cba   C2v^5         P21ca                        P -2ac 2a', &             
  'x,y,z','x+1/2,-y,-z','x+1/2,y,-z+1/2','x,-y,z+1/2' /                                                                             

 data spacegroupdata( 623: 627)/ &                              ! group index = 148, number of symmetry op =  4
  '* 29:bca    C2v^5         Pc21b                        P -2bc -2c', &            
  'x,y,z','-x,y+1/2,-z','x,y+1/2,-z+1/2','-x,y,z+1/2' /                                                                             

 data spacegroupdata( 628: 632)/ &                              ! group index = 149, number of symmetry op =  4
  '* 29:a-cb   C2v^5         Pb21a                        P -2a -2ab', &            
  'x,y,z','-x,y+1/2,-z','x+1/2,y,-z','-x+1/2,y+1/2,z' /                                                                             

 data spacegroupdata( 633: 637)/ &                              ! group index = 150, number of symmetry op =  4
  '* 30        C2v^6         Pnc2                         P 2 -2bc', &              
  'x,y,z','-x,-y,z','-x,y+1/2,z+1/2','x,-y+1/2,z+1/2' /                                                                             

 data spacegroupdata( 638: 642)/ &                              ! group index = 151, number of symmetry op =  4
  '* 30:ba-c   C2v^6         Pcn2                         P 2 -2ac', &              
  'x,y,z','-x,-y,z','-x+1/2,y,z+1/2','x+1/2,-y,z+1/2' /                                                                             

 data spacegroupdata( 643: 647)/ &                              ! group index = 152, number of symmetry op =  4
  '* 30:cab    C2v^6         P2na                         P -2ac 2', &              
  'x,y,z','x,-y,-z','x+1/2,y,-z+1/2','x+1/2,-y,z+1/2' /                                                                             

 data spacegroupdata( 648: 652)/ &                              ! group index = 153, number of symmetry op =  4
  '* 30:-cba   C2v^6         P2an                         P -2ab 2', &              
  'x,y,z','x,-y,-z','x+1/2,y+1/2,-z','x+1/2,-y+1/2,z' /                                                                             

 data spacegroupdata( 653: 657)/ &                              ! group index = 154, number of symmetry op =  4
  '* 30:bca    C2v^6         Pb2n                         P -2ab -2ab', &           
  'x,y,z','-x,y,-z','x+1/2,y+1/2,-z','-x+1/2,y+1/2,z' /                                                                             

 data spacegroupdata( 658: 662)/ &                              ! group index = 155, number of symmetry op =  4
  '* 30:a-cb   C2v^6         Pn2b                         P -2bc -2bc', &           
  'x,y,z','-x,y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z+1/2' /                                                                             

 data spacegroupdata( 663: 667)/ &                              ! group index = 156, number of symmetry op =  4
  '* 31        C2v^7         Pmn21                        P 2ac -2', &              
  'x,y,z','-x+1/2,-y,z+1/2','-x,y,z','x+1/2,-y,z+1/2' /                                                                             

 data spacegroupdata( 668: 672)/ &                              ! group index = 157, number of symmetry op =  4
  '* 31:ba-c   C2v^7         Pnm21                        P 2bc -2bc', &            
  'x,y,z','-x,-y+1/2,z+1/2','-x,y+1/2,z+1/2','x,-y,z' /                                                                             

 data spacegroupdata( 673: 677)/ &                              ! group index = 158, number of symmetry op =  4
  '* 31:cab    C2v^7         P21mn                        P -2ab 2ab', &            
  'x,y,z','x+1/2,-y+1/2,-z','x+1/2,y+1/2,-z','x,-y,z' /                                                                             

 data spacegroupdata( 678: 682)/ &                              ! group index = 159, number of symmetry op =  4
  '* 31:-cba   C2v^7         P21nm                        P -2 2ac', &              
  'x,y,z','x+1/2,-y,-z+1/2','x,y,-z','x+1/2,-y,z+1/2' /                                                                             

 data spacegroupdata( 683: 687)/ &                              ! group index = 160, number of symmetry op =  4
  '* 31:bca    C2v^7         Pn21m                        P -2 -2bc', &             
  'x,y,z','-x,y+1/2,-z+1/2','x,y,-z','-x,y+1/2,z+1/2' /                                                                             

 data spacegroupdata( 688: 692)/ &                              ! group index = 161, number of symmetry op =  4
  '* 31:a-cb   C2v^7         Pm21n                        P -2ab -2', &             
  'x,y,z','-x+1/2,y+1/2,-z','x+1/2,y+1/2,-z','-x,y,z' /                                                                             

 data spacegroupdata( 693: 697)/ &                              ! group index = 162, number of symmetry op =  4
  '* 32        C2v^8         Pba2                         P 2 -2ab', &              
  'x,y,z','-x,-y,z','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z' /                                                                             

 data spacegroupdata( 698: 702)/ &                              ! group index = 163, number of symmetry op =  4
  '* 32:cab    C2v^8         P2cb                         P -2bc 2', &              
  'x,y,z','x,-y,-z','x,y+1/2,-z+1/2','x,-y+1/2,z+1/2' /                                                                             

 data spacegroupdata( 703: 707)/ &                              ! group index = 164, number of symmetry op =  4
  '* 32:bca    C2v^8         Pc2a                         P -2ac -2ac', &           
  'x,y,z','-x,y,-z','x+1/2,y,-z+1/2','-x+1/2,y,z+1/2' /                                                                             

 data spacegroupdata( 708: 712)/ &                              ! group index = 165, number of symmetry op =  4
  '* 33        C2v^9         Pna21                        P 2c -2n', &              
  'x,y,z','-x,-y,z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z' /                                                                     

 data spacegroupdata( 713: 717)/ &                              ! group index = 166, number of symmetry op =  4
  '* 33:ba-c   C2v^9         Pbn21                        P 2c -2ab', &             
  'x,y,z','-x,-y,z+1/2','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z+1/2' /                                                                     

 data spacegroupdata( 718: 722)/ &                              ! group index = 167, number of symmetry op =  4
  '* 33:cab    C2v^9         P21nb                        P -2bc 2a', &             
  'x,y,z','x+1/2,-y,-z','x,y+1/2,-z+1/2','x+1/2,-y+1/2,z+1/2' /                                                                     

 data spacegroupdata( 723: 727)/ &                              ! group index = 168, number of symmetry op =  4
  '* 33:-cba   C2v^9         P21cn                        P -2n 2a', &              
  'x,y,z','x+1/2,-y,-z','x+1/2,y+1/2,-z+1/2','x,-y+1/2,z+1/2' /                                                                     

 data spacegroupdata( 728: 732)/ &                              ! group index = 169, number of symmetry op =  4
  '* 33:bca    C2v^9         Pc21n                        P -2n -2ac', &            
  'x,y,z','-x,y+1/2,-z','x+1/2,y+1/2,-z+1/2','-x+1/2,y,z+1/2' /                                                                     

 data spacegroupdata( 733: 737)/ &                              ! group index = 170, number of symmetry op =  4
  '* 33:a-cb   C2v^9         Pn21a                        P -2ac -2n', &            
  'x,y,z','-x,y+1/2,-z','x+1/2,y,-z+1/2','-x+1/2,y+1/2,z+1/2' /                                                                     

 data spacegroupdata( 738: 742)/ &                              ! group index = 171, number of symmetry op =  4
  '* 34        C2v^10        Pnn2                         P 2 -2n', &               
  'x,y,z','-x,-y,z','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2' /                                                                     

 data spacegroupdata( 743: 747)/ &                              ! group index = 172, number of symmetry op =  4
  '* 34:cab    C2v^10        P2nn                         P -2n 2', &               
  'x,y,z','x,-y,-z','x+1/2,y+1/2,-z+1/2','x+1/2,-y+1/2,z+1/2' /                                                                     

 data spacegroupdata( 748: 752)/ &                              ! group index = 173, number of symmetry op =  4
  '* 34:bca    C2v^10        Pn2n                         P -2n -2n', &             
  'x,y,z','-x,y,-z','x+1/2,y+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2' /                                                                     

 data spacegroupdata( 753: 757)/ &                              ! group index = 174, number of symmetry op =  4
  '* 35        C2v^11        Cmm2                         C 2 -2', &                
  'x,y,z','-x,-y,z','-x,y,z','x,-y,z' /                                                                                             

 data spacegroupdata( 758: 762)/ &                              ! group index = 175, number of symmetry op =  4
  '* 35:cab    C2v^11        A2mm                         A -2 2', &                
  'x,y,z','x,-y,-z','x,y,-z','x,-y,z' /                                                                                             

 data spacegroupdata( 763: 767)/ &                              ! group index = 176, number of symmetry op =  4
  '* 35:bca    C2v^11        Bm2m                         B -2 -2', &               
  'x,y,z','-x,y,-z','x,y,-z','-x,y,z' /                                                                                             

 data spacegroupdata( 768: 772)/ &                              ! group index = 177, number of symmetry op =  4
  '* 36        C2v^12        Cmc21                        C 2c -2', &               
  'x,y,z','-x,-y,z+1/2','-x,y,z','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 773: 777)/ &                              ! group index = 178, number of symmetry op =  4
  '* 36:ba-c   C2v^12        Ccm21                        C 2c -2c', &              
  'x,y,z','-x,-y,z+1/2','-x,y,z+1/2','x,-y,z' /                                                                                     

 data spacegroupdata( 778: 782)/ &                              ! group index = 179, number of symmetry op =  4
  '* 36:cab    C2v^12        A21ma                        A -2a 2a', &              
  'x,y,z','x+1/2,-y,-z','x+1/2,y,-z','x,-y,z' /                                                                                     

 data spacegroupdata( 783: 787)/ &                              ! group index = 180, number of symmetry op =  4
  '* 36:-cba   C2v^12        A21am                        A -2 2a', &               
  'x,y,z','x+1/2,-y,-z','x,y,-z','x+1/2,-y,z' /                                                                                     

 data spacegroupdata( 788: 792)/ &                              ! group index = 181, number of symmetry op =  4
  '* 36:bca    C2v^12        Bb21m                        B -2 -2b', &              
  'x,y,z','-x,y+1/2,-z','x,y,-z','-x,y+1/2,z' /                                                                                     

 data spacegroupdata( 793: 797)/ &                              ! group index = 182, number of symmetry op =  4
  '* 36:a-cb   C2v^12        Bm21b                        B -2b -2', &              
  'x,y,z','-x,y+1/2,-z','x,y+1/2,-z','-x,y,z' /                                                                                     

 data spacegroupdata( 798: 802)/ &                              ! group index = 183, number of symmetry op =  4
  '* 37        C2v^13        Ccc2                         C 2 -2c', &               
  'x,y,z','-x,-y,z','-x,y,z+1/2','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 803: 807)/ &                              ! group index = 184, number of symmetry op =  4
  '* 37:cab    C2v^13        A2aa                         A -2a 2', &               
  'x,y,z','x,-y,-z','x+1/2,y,-z','x+1/2,-y,z' /                                                                                     

 data spacegroupdata( 808: 812)/ &                              ! group index = 185, number of symmetry op =  4
  '* 37:bca    C2v^13        Bb2b                         B -2b -2b', &             
  'x,y,z','-x,y,-z','x,y+1/2,-z','-x,y+1/2,z' /                                                                                     

 data spacegroupdata( 813: 817)/ &                              ! group index = 186, number of symmetry op =  4
  '* 38        C2v^14        Amm2                         A 2 -2', &                
  'x,y,z','-x,-y,z','-x,y,z','x,-y,z' /                                                                                             

 data spacegroupdata( 818: 822)/ &                              ! group index = 187, number of symmetry op =  4
  '* 38:ba-c   C2v^14        Bmm2                         B 2 -2', &                
  'x,y,z','-x,-y,z','-x,y,z','x,-y,z' /                                                                                             

 data spacegroupdata( 823: 827)/ &                              ! group index = 188, number of symmetry op =  4
  '* 38:cab    C2v^14        B2mm                         B -2 2', &                
  'x,y,z','x,-y,-z','x,y,-z','x,-y,z' /                                                                                             

 data spacegroupdata( 828: 832)/ &                              ! group index = 189, number of symmetry op =  4
  '* 38:-cba   C2v^14        C2mm                         C -2 2', &                
  'x,y,z','x,-y,-z','x,y,-z','x,-y,z' /                                                                                             

 data spacegroupdata( 833: 837)/ &                              ! group index = 190, number of symmetry op =  4
  '* 38:bca    C2v^14        Cm2m                         C -2 -2', &               
  'x,y,z','-x,y,-z','x,y,-z','-x,y,z' /                                                                                             

 data spacegroupdata( 838: 842)/ &                              ! group index = 191, number of symmetry op =  4
  '* 38:a-cb   C2v^14        Am2m                         A -2 -2', &               
  'x,y,z','-x,y,-z','x,y,-z','-x,y,z' /                                                                                             

 data spacegroupdata( 843: 847)/ &                              ! group index = 192, number of symmetry op =  4
  '* 39        C2v^15        Abm2                         A 2 -2c', &               
  'x,y,z','-x,-y,z','-x,y,z+1/2','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 848: 852)/ &                              ! group index = 193, number of symmetry op =  4
  '* 39:ba-c   C2v^15        Bma2                         B 2 -2c', &               
  'x,y,z','-x,-y,z','-x,y,z+1/2','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 853: 857)/ &                              ! group index = 194, number of symmetry op =  4
  '* 39:cab    C2v^15        B2cm                         B -2c 2', &               
  'x,y,z','x,-y,-z','x,y,-z+1/2','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 858: 862)/ &                              ! group index = 195, number of symmetry op =  4
  '* 39:-cba   C2v^15        C2mb                         C -2b 2', &               
  'x,y,z','x,-y,-z','x,y+1/2,-z','x,-y+1/2,z' /                                                                                     

 data spacegroupdata( 863: 867)/ &                              ! group index = 196, number of symmetry op =  4
  '* 39:bca    C2v^15        Cm2a                         C -2b -2b', &             
  'x,y,z','-x,y,-z','x,y+1/2,-z','-x,y+1/2,z' /                                                                                     

 data spacegroupdata( 868: 872)/ &                              ! group index = 197, number of symmetry op =  4
  '* 39:a-cb   C2v^15        Ac2m                         A -2c -2c', &             
  'x,y,z','-x,y,-z','x,y,-z+1/2','-x,y,z+1/2' /                                                                                     

 data spacegroupdata( 873: 877)/ &                              ! group index = 198, number of symmetry op =  4
  '* 40        C2v^16        Ama2                         A 2 -2a', &               
  'x,y,z','-x,-y,z','-x+1/2,y,z','x+1/2,-y,z' /                                                                                     

 data spacegroupdata( 878: 882)/ &                              ! group index = 199, number of symmetry op =  4
  '* 40:ba-c   C2v^16        Bbm2                         B 2 -2b', &               
  'x,y,z','-x,-y,z','-x,y+1/2,z','x,-y+1/2,z' /                                                                                     

 data spacegroupdata( 883: 887)/ &                              ! group index = 200, number of symmetry op =  4
  '* 40:cab    C2v^16        B2mb                         B -2b 2', &               
  'x,y,z','x,-y,-z','x,y+1/2,-z','x,-y+1/2,z' /                                                                                     

 data spacegroupdata( 888: 892)/ &                              ! group index = 201, number of symmetry op =  4
  '* 40:-cba   C2v^16        C2cm                         C -2c 2', &               
  'x,y,z','x,-y,-z','x,y,-z+1/2','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 893: 897)/ &                              ! group index = 202, number of symmetry op =  4
  '* 40:bca    C2v^16        Cc2m                         C -2c -2c', &             
  'x,y,z','-x,y,-z','x,y,-z+1/2','-x,y,z+1/2' /                                                                                     

 data spacegroupdata( 898: 902)/ &                              ! group index = 203, number of symmetry op =  4
  '* 40:a-cb   C2v^16        Am2a                         A -2a -2a', &             
  'x,y,z','-x,y,-z','x+1/2,y,-z','-x+1/2,y,z' /                                                                                     

 data spacegroupdata( 903: 907)/ &                              ! group index = 204, number of symmetry op =  4
  '* 41        C2v^17        Aba2                         A 2 -2ac', &              
  'x,y,z','-x,-y,z','-x+1/2,y,z+1/2','x+1/2,-y,z+1/2' /                                                                             

 data spacegroupdata( 908: 912)/ &                              ! group index = 205, number of symmetry op =  4
  '* 41:ba-c   C2v^17        Bba2                         B 2 -2bc', &              
  'x,y,z','-x,-y,z','-x,y+1/2,z+1/2','x,-y+1/2,z+1/2' /                                                                             

 data spacegroupdata( 913: 917)/ &                              ! group index = 206, number of symmetry op =  4
  '* 41:cab    C2v^17        B2cb                         B -2bc 2', &              
  'x,y,z','x,-y,-z','x,y+1/2,-z+1/2','x,-y+1/2,z+1/2' /                                                                             

 data spacegroupdata( 918: 922)/ &                              ! group index = 207, number of symmetry op =  4
  '* 41:-cba   C2v^17        C2cb                         C -2bc 2', &              
  'x,y,z','x,-y,-z','x,y+1/2,-z+1/2','x,-y+1/2,z+1/2' /                                                                             

 data spacegroupdata( 923: 927)/ &                              ! group index = 208, number of symmetry op =  4
  '* 41:bca    C2v^17        Cc2a                         C -2bc -2bc', &           
  'x,y,z','-x,y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z+1/2' /                                                                             

 data spacegroupdata( 928: 932)/ &                              ! group index = 209, number of symmetry op =  4
  '* 41:a-cb   C2v^17        Ac2a                         A -2ac -2ac', &           
  'x,y,z','-x,y,-z','x+1/2,y,-z+1/2','-x+1/2,y,z+1/2' /                                                                             

 data spacegroupdata( 933: 937)/ &                              ! group index = 210, number of symmetry op =  4
  '* 42        C2v^18        Fmm2                         F 2 -2', &                
  'x,y,z','-x,-y,z','-x,y,z','x,-y,z' /                                                                                             

 data spacegroupdata( 938: 942)/ &                              ! group index = 211, number of symmetry op =  4
  '* 42:cab    C2v^18        F2mm                         F -2 2', &                
  'x,y,z','x,-y,-z','x,y,-z','x,-y,z' /                                                                                             

 data spacegroupdata( 943: 947)/ &                              ! group index = 212, number of symmetry op =  4
  '* 42:bca    C2v^18        Fm2m                         F -2 -2', &               
  'x,y,z','-x,y,-z','x,y,-z','-x,y,z' /                                                                                             

 data spacegroupdata( 948: 952)/ &                              ! group index = 213, number of symmetry op =  4
  '* 43        C2v^19        Fdd2                         F 2 -2d', &               
  'x,y,z','-x,-y,z','-x+1/4,y+1/4,z+1/4','x+1/4,-y+1/4,z+1/4' /                                                                     

 data spacegroupdata( 953: 957)/ &                              ! group index = 214, number of symmetry op =  4
  '* 43:cab    C2v^19        F2dd                         F -2d 2', &               
  'x,y,z','x,-y,-z','x+1/4,y+1/4,-z+1/4','x+1/4,-y+1/4,z+1/4' /                                                                     

 data spacegroupdata( 958: 962)/ &                              ! group index = 215, number of symmetry op =  4
  '* 43:bca    C2v^19        Fd2d                         F -2d -2d', &             
  'x,y,z','-x,y,-z','x+1/4,y+1/4,-z+1/4','-x+1/4,y+1/4,z+1/4' /                                                                     

 data spacegroupdata( 963: 967)/ &                              ! group index = 216, number of symmetry op =  4
  '* 44        C2v^20        Imm2                         I 2 -2', &                
  'x,y,z','-x,-y,z','-x,y,z','x,-y,z' /                                                                                             

 data spacegroupdata( 968: 972)/ &                              ! group index = 217, number of symmetry op =  4
  '* 44:cab    C2v^20        I2mm                         I -2 2', &                
  'x,y,z','x,-y,-z','x,y,-z','x,-y,z' /                                                                                             

 data spacegroupdata( 973: 977)/ &                              ! group index = 218, number of symmetry op =  4
  '* 44:bca    C2v^20        Im2m                         I -2 -2', &               
  'x,y,z','-x,y,-z','x,y,-z','-x,y,z' /                                                                                             

 data spacegroupdata( 978: 982)/ &                              ! group index = 219, number of symmetry op =  4
  '* 45        C2v^21        Iba2                         I 2 -2c', &               
  'x,y,z','-x,-y,z','-x,y,z+1/2','x,-y,z+1/2' /                                                                                     

 data spacegroupdata( 983: 987)/ &                              ! group index = 220, number of symmetry op =  4
  '* 45:cab    C2v^21        I2cb                         I -2a 2', &               
  'x,y,z','x,-y,-z','x+1/2,y,-z','x+1/2,-y,z' /                                                                                     

 data spacegroupdata( 988: 992)/ &                              ! group index = 221, number of symmetry op =  4
  '* 45:bca    C2v^21        Ic2a                         I -2b -2b', &             
  'x,y,z','-x,y,-z','x,y+1/2,-z','-x,y+1/2,z' /                                                                                     

 data spacegroupdata( 993: 997)/ &                              ! group index = 222, number of symmetry op =  4
  '* 46        C2v^22        Ima2                         I 2 -2a', &               
  'x,y,z','-x,-y,z','-x+1/2,y,z','x+1/2,-y,z' /                                                                                     

 data spacegroupdata( 998:1002)/ &                              ! group index = 223, number of symmetry op =  4
  '* 46:ba-c   C2v^22        Ibm2                         I 2 -2b', &               
  'x,y,z','-x,-y,z','-x,y+1/2,z','x,-y+1/2,z' /                                                                                     

 data spacegroupdata(1003:1007)/ &                              ! group index = 224, number of symmetry op =  4
  '* 46:cab    C2v^22        I2mb                         I -2b 2', &               
  'x,y,z','x,-y,-z','x,y+1/2,-z','x,-y+1/2,z' /                                                                                     

 data spacegroupdata(1008:1012)/ &                              ! group index = 225, number of symmetry op =  4
  '* 46:-cba   C2v^22        I2cm                         I -2c 2', &               
  'x,y,z','x,-y,-z','x,y,-z+1/2','x,-y,z+1/2' /                                                                                     

 data spacegroupdata(1013:1017)/ &                              ! group index = 226, number of symmetry op =  4
  '* 46:bca    C2v^22        Ic2m                         I -2c -2c', &             
  'x,y,z','-x,y,-z','x,y,-z+1/2','-x,y,z+1/2' /                                                                                     

 data spacegroupdata(1018:1022)/ &                              ! group index = 227, number of symmetry op =  4
  '* 46:a-cb   C2v^22        Im2a                         I -2a -2a', &             
  'x,y,z','-x,y,-z','x+1/2,y,-z','-x+1/2,y,z' /                                                                                     

 data spacegroupdata(1023:1031)/ &                              ! group index = 228, number of symmetry op =  8
  '* 47        D2h^1         Pmmm                        -P 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z','-x,-y,-z','x,y,-z','-x,y,z','x,-y,z' /                                                     

 data spacegroupdata(1032:1040)/ &                              ! group index = 229, number of symmetry op =  8
  '* 48:1      D2h^2         Pnnn:1                       P 2 2 -1n', &             
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x+1/2,y+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(1041:1049)/ &                              ! group index = 230, number of symmetry op =  8
  '* 48:2      D2h^2         Pnnn:2                      -P 2ab 2bc', &             
  'x,y,z','-x+1/2,-y+1/2,z','x,-y+1/2,-z+1/2','-x+1/2,y,-z+1/2','-x,-y,-z','x+1/2,y+1/2,-z','-x,y+1/2,z+1/2','x+1/2,-y,z+1/2' /     

 data spacegroupdata(1050:1058)/ &                              ! group index = 231, number of symmetry op =  8
  '* 49        D2h^3         Pccm                        -P 2 2c', &                
  'x,y,z','-x,-y,z','x,-y,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x,y,-z','-x,y,z+1/2','x,-y,z+1/2' /                                     

 data spacegroupdata(1059:1067)/ &                              ! group index = 232, number of symmetry op =  8
  '* 49:cab    D2h^3         Pmaa                        -P 2a 2', &                
  'x,y,z','-x+1/2,-y,z','x,-y,-z','-x+1/2,y,-z','-x,-y,-z','x+1/2,y,-z','-x,y,z','x+1/2,-y,z' /                                     

 data spacegroupdata(1068:1076)/ &                              ! group index = 233, number of symmetry op =  8
  '* 49:bca    D2h^3         Pbmb                        -P 2b 2b', &               
  'x,y,z','-x,-y+1/2,z','x,-y+1/2,-z','-x,y,-z','-x,-y,-z','x,y+1/2,-z','-x,y+1/2,z','x,-y,z' /                                     

 data spacegroupdata(1077:1085)/ &                              ! group index = 234, number of symmetry op =  8
  '* 50:1      D2h^4         Pban:1                       P 2 2 -1ab', &            
  'x,y,z','-x+1/2,-y+1/2,-z','-x,-y,z','x,-y,-z','-x,y,-z','x+1/2,y+1/2,-z','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z' /                     

 data spacegroupdata(1086:1094)/ &                              ! group index = 235, number of symmetry op =  8
  '* 50:2      D2h^4         Pban:2                      -P 2ab 2b', &              
  'x,y,z','-x+1/2,-y+1/2,z','x,-y+1/2,-z','-x+1/2,y,-z','-x,-y,-z','x+1/2,y+1/2,-z','-x,y+1/2,z','x+1/2,-y,z' /                     

 data spacegroupdata(1095:1103)/ &                              ! group index = 236, number of symmetry op =  8
  '* 50:1cab   D2h^4         Pncb:1                       P 2 2 -1bc', &            
  'x,y,z','-x,-y+1/2,-z+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z+1/2','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1104:1112)/ &                              ! group index = 237, number of symmetry op =  8
  '* 50:2cab   D2h^4         Pncb:2                      -P 2b 2bc', &              
  'x,y,z','-x,-y+1/2,z','x,-y+1/2,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x,y+1/2,-z','-x,y+1/2,z+1/2','x,-y,z+1/2' /                     

 data spacegroupdata(1113:1121)/ &                              ! group index = 238, number of symmetry op =  8
  '* 50:1bca   D2h^4         Pcna:1                       P 2 2 -1ac', &            
  'x,y,z','-x+1/2,-y,-z+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x+1/2,y,-z+1/2','-x+1/2,y,z+1/2','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1122:1130)/ &                              ! group index = 239, number of symmetry op =  8
  '* 50:2bca   D2h^4         Pcna:2                      -P 2a 2c', &               
  'x,y,z','-x+1/2,-y,z','x,-y,-z+1/2','-x+1/2,y,-z+1/2','-x,-y,-z','x+1/2,y,-z','-x,y,z+1/2','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1131:1139)/ &                              ! group index = 240, number of symmetry op =  8
  '* 51        D2h^5         Pmma                        -P 2a 2a', &               
  'x,y,z','-x+1/2,-y,z','x+1/2,-y,-z','-x,y,-z','-x,-y,-z','x+1/2,y,-z','-x+1/2,y,z','x,-y,z' /                                     

 data spacegroupdata(1140:1148)/ &                              ! group index = 241, number of symmetry op =  8
  '* 51:ba-c   D2h^5         Pmmb                        -P 2b 2', &                
  'x,y,z','-x,-y+1/2,z','x,-y,-z','-x,y+1/2,-z','-x,-y,-z','x,y+1/2,-z','-x,y,z','x,-y+1/2,z' /                                     

 data spacegroupdata(1149:1157)/ &                              ! group index = 242, number of symmetry op =  8
  '* 51:cab    D2h^5         Pbmm                        -P 2 2b', &                
  'x,y,z','-x,-y,z','x,-y+1/2,-z','-x,y+1/2,-z','-x,-y,-z','x,y,-z','-x,y+1/2,z','x,-y+1/2,z' /                                     

 data spacegroupdata(1158:1166)/ &                              ! group index = 243, number of symmetry op =  8
  '* 51:-cba   D2h^5         Pcmm                        -P 2c 2c', &               
  'x,y,z','-x,-y,z+1/2','x,-y,-z+1/2','-x,y,-z','-x,-y,-z','x,y,-z+1/2','-x,y,z+1/2','x,-y,z' /                                     

 data spacegroupdata(1167:1175)/ &                              ! group index = 244, number of symmetry op =  8
  '* 51:bca    D2h^5         Pmcm                        -P 2c 2', &                
  'x,y,z','-x,-y,z+1/2','x,-y,-z','-x,y,-z+1/2','-x,-y,-z','x,y,-z+1/2','-x,y,z','x,-y,z+1/2' /                                     

 data spacegroupdata(1176:1184)/ &                              ! group index = 245, number of symmetry op =  8
  '* 51:a-cb   D2h^5         Pmam                        -P 2 2a', &                
  'x,y,z','-x,-y,z','x+1/2,-y,-z','-x+1/2,y,-z','-x,-y,-z','x,y,-z','-x+1/2,y,z','x+1/2,-y,z' /                                     

 data spacegroupdata(1185:1193)/ &                              ! group index = 246, number of symmetry op =  8
  '* 52        D2h^6         Pnna                        -P 2a 2bc', &              
  'x,y,z','-x+1/2,-y,z','x,-y+1/2,-z+1/2','-x+1/2,y+1/2,-z+1/2','-x,-y,-z','x+1/2,y,-z','-x,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(1194:1202)/ &                              ! group index = 247, number of symmetry op =  8
  '* 52:ba-c   D2h^6         Pnnb                        -P 2b 2n', &               
  'x,y,z','-x,-y+1/2,z','x+1/2,-y+1/2,-z+1/2','-x+1/2,y,-z+1/2','-x,-y,-z','x,y+1/2,-z','-x+1/2,y+1/2,z+1/2','x+1/2,-y,z+1/2' /     

 data spacegroupdata(1203:1211)/ &                              ! group index = 248, number of symmetry op =  8
  '* 52:cab    D2h^6         Pbnn                        -P 2n 2b', &               
  'x,y,z','-x+1/2,-y+1/2,z+1/2','x,-y+1/2,-z','-x+1/2,y,-z+1/2','-x,-y,-z','x+1/2,y+1/2,-z+1/2','-x,y+1/2,z','x+1/2,-y,z+1/2' /     

 data spacegroupdata(1212:1220)/ &                              ! group index = 249, number of symmetry op =  8
  '* 52:-cba   D2h^6         Pcnn                        -P 2ab 2c', &              
  'x,y,z','-x+1/2,-y+1/2,z','x,-y,-z+1/2','-x+1/2,y+1/2,-z+1/2','-x,-y,-z','x+1/2,y+1/2,-z','-x,y,z+1/2','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(1221:1229)/ &                              ! group index = 250, number of symmetry op =  8
  '* 52:bca    D2h^6         Pncn                        -P 2ab 2n', &              
  'x,y,z','-x+1/2,-y+1/2,z','x+1/2,-y+1/2,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x+1/2,y+1/2,-z','-x+1/2,y+1/2,z+1/2','x,-y,z+1/2' /     

 data spacegroupdata(1230:1238)/ &                              ! group index = 251, number of symmetry op =  8
  '* 52:a-cb   D2h^6         Pnan                        -P 2n 2bc', &              
  'x,y,z','-x+1/2,-y+1/2,z+1/2','x,-y+1/2,-z+1/2','-x+1/2,y,-z','-x,-y,-z','x+1/2,y+1/2,-z+1/2','-x,y+1/2,z+1/2','x+1/2,-y,z' /     

 data spacegroupdata(1239:1247)/ &                              ! group index = 252, number of symmetry op =  8
  '* 53        D2h^7         Pmna                        -P 2ac 2', &               
  'x,y,z','-x+1/2,-y,z+1/2','x,-y,-z','-x+1/2,y,-z+1/2','-x,-y,-z','x+1/2,y,-z+1/2','-x,y,z','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1248:1256)/ &                              ! group index = 253, number of symmetry op =  8
  '* 53:ba-c   D2h^7         Pnmb                        -P 2bc 2bc', &             
  'x,y,z','-x,-y+1/2,z+1/2','x,-y+1/2,-z+1/2','-x,y,-z','-x,-y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z+1/2','x,-y,z' /                     

 data spacegroupdata(1257:1265)/ &                              ! group index = 254, number of symmetry op =  8
  '* 53:cab    D2h^7         Pbmn                        -P 2ab 2ab', &             
  'x,y,z','-x+1/2,-y+1/2,z','x+1/2,-y+1/2,-z','-x,y,-z','-x,-y,-z','x+1/2,y+1/2,-z','-x+1/2,y+1/2,z','x,-y,z' /                     

 data spacegroupdata(1266:1274)/ &                              ! group index = 255, number of symmetry op =  8
  '* 53:-cba   D2h^7         Pcnm                        -P 2 2ac', &               
  'x,y,z','-x,-y,z','x+1/2,-y,-z+1/2','-x+1/2,y,-z+1/2','-x,-y,-z','x,y,-z','-x+1/2,y,z+1/2','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1275:1283)/ &                              ! group index = 256, number of symmetry op =  8
  '* 53:bca    D2h^7         Pncm                        -P 2 2bc', &               
  'x,y,z','-x,-y,z','x,-y+1/2,-z+1/2','-x,y+1/2,-z+1/2','-x,-y,-z','x,y,-z','-x,y+1/2,z+1/2','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1284:1292)/ &                              ! group index = 257, number of symmetry op =  8
  '* 53:a-cb   D2h^7         Pman                        -P 2ab 2', &               
  'x,y,z','-x+1/2,-y+1/2,z','x,-y,-z','-x+1/2,y+1/2,-z','-x,-y,-z','x+1/2,y+1/2,-z','-x,y,z','x+1/2,-y+1/2,z' /                     

 data spacegroupdata(1293:1301)/ &                              ! group index = 258, number of symmetry op =  8
  '* 54        D2h^8         Pcca                        -P 2a 2ac', &              
  'x,y,z','-x+1/2,-y,z','x+1/2,-y,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x+1/2,y,-z','-x+1/2,y,z+1/2','x,-y,z+1/2' /                     

 data spacegroupdata(1302:1310)/ &                              ! group index = 259, number of symmetry op =  8
  '* 54:ba-c   D2h^8         Pccb                        -P 2b 2c', &               
  'x,y,z','-x,-y+1/2,z','x,-y,-z+1/2','-x,y+1/2,-z+1/2','-x,-y,-z','x,y+1/2,-z','-x,y,z+1/2','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1311:1319)/ &                              ! group index = 260, number of symmetry op =  8
  '* 54:cab    D2h^8         Pbaa                        -P 2a 2b', &               
  'x,y,z','-x+1/2,-y,z','x,-y+1/2,-z','-x+1/2,y+1/2,-z','-x,-y,-z','x+1/2,y,-z','-x,y+1/2,z','x+1/2,-y+1/2,z' /                     

 data spacegroupdata(1320:1328)/ &                              ! group index = 261, number of symmetry op =  8
  '* 54:-cba   D2h^8         Pcaa                        -P 2ac 2c', &              
  'x,y,z','-x+1/2,-y,z+1/2','x,-y,-z+1/2','-x+1/2,y,-z','-x,-y,-z','x+1/2,y,-z+1/2','-x,y,z+1/2','x+1/2,-y,z' /                     

 data spacegroupdata(1329:1337)/ &                              ! group index = 262, number of symmetry op =  8
  '* 54:bca    D2h^8         Pbcb                        -P 2bc 2b', &              
  'x,y,z','-x,-y+1/2,z+1/2','x,-y+1/2,-z','-x,y,-z+1/2','-x,-y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z','x,-y,z+1/2' /                     

 data spacegroupdata(1338:1346)/ &                              ! group index = 263, number of symmetry op =  8
  '* 54:a-cb   D2h^8         Pbab                        -P 2b 2ab', &              
  'x,y,z','-x,-y+1/2,z','x+1/2,-y+1/2,-z','-x+1/2,y,-z','-x,-y,-z','x,y+1/2,-z','-x+1/2,y+1/2,z','x+1/2,-y,z' /                     

 data spacegroupdata(1347:1355)/ &                              ! group index = 264, number of symmetry op =  8
  '* 55        D2h^9         Pbam                        -P 2 2ab', &               
  'x,y,z','-x,-y,z','x+1/2,-y+1/2,-z','-x+1/2,y+1/2,-z','-x,-y,-z','x,y,-z','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z' /                     

 data spacegroupdata(1356:1364)/ &                              ! group index = 265, number of symmetry op =  8
  '* 55:cab    D2h^9         Pmcb                        -P 2bc 2', &               
  'x,y,z','-x,-y+1/2,z+1/2','x,-y,-z','-x,y+1/2,-z+1/2','-x,-y,-z','x,y+1/2,-z+1/2','-x,y,z','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1365:1373)/ &                              ! group index = 266, number of symmetry op =  8
  '* 55:bca    D2h^9         Pcma                        -P 2ac 2ac', &             
  'x,y,z','-x+1/2,-y,z+1/2','x+1/2,-y,-z+1/2','-x,y,-z','-x,-y,-z','x+1/2,y,-z+1/2','-x+1/2,y,z+1/2','x,-y,z' /                     

 data spacegroupdata(1374:1382)/ &                              ! group index = 267, number of symmetry op =  8
  '* 56        D2h^10        Pccn                        -P 2ab 2ac', &             
  'x,y,z','-x+1/2,-y+1/2,z','x+1/2,-y,-z+1/2','-x,y+1/2,-z+1/2','-x,-y,-z','x+1/2,y+1/2,-z','-x+1/2,y,z+1/2','x,-y+1/2,z+1/2' /     

 data spacegroupdata(1383:1391)/ &                              ! group index = 268, number of symmetry op =  8
  '* 56:cab    D2h^10        Pnaa                        -P 2ac 2bc', &             
  'x,y,z','-x+1/2,-y,z+1/2','x,-y+1/2,-z+1/2','-x+1/2,y+1/2,-z','-x,-y,-z','x+1/2,y,-z+1/2','-x,y+1/2,z+1/2','x+1/2,-y+1/2,z' /     

 data spacegroupdata(1392:1400)/ &                              ! group index = 269, number of symmetry op =  8
  '* 56:bca    D2h^10        Pbnb                        -P 2bc 2ab', &             
  'x,y,z','-x,-y+1/2,z+1/2','x+1/2,-y+1/2,-z','-x+1/2,y,-z+1/2','-x,-y,-z','x,y+1/2,-z+1/2','-x+1/2,y+1/2,z','x+1/2,-y,z+1/2' /     

 data spacegroupdata(1401:1409)/ &                              ! group index = 270, number of symmetry op =  8
  '* 57        D2h^11        Pbcm                        -P 2c 2b', &               
  'x,y,z','-x,-y,z+1/2','x,-y+1/2,-z','-x,y+1/2,-z+1/2','-x,-y,-z','x,y,-z+1/2','-x,y+1/2,z','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1410:1418)/ &                              ! group index = 271, number of symmetry op =  8
  '* 57:ba-c   D2h^11        Pcam                        -P 2c 2ac', &              
  'x,y,z','-x,-y,z+1/2','x+1/2,-y,-z+1/2','-x+1/2,y,-z','-x,-y,-z','x,y,-z+1/2','-x+1/2,y,z+1/2','x+1/2,-y,z' /                     

 data spacegroupdata(1419:1427)/ &                              ! group index = 272, number of symmetry op =  8
  '* 57:cab    D2h^11        Pmca                        -P 2ac 2a', &              
  'x,y,z','-x+1/2,-y,z+1/2','x+1/2,-y,-z','-x,y,-z+1/2','-x,-y,-z','x+1/2,y,-z+1/2','-x+1/2,y,z','x,-y,z+1/2' /                     

 data spacegroupdata(1428:1436)/ &                              ! group index = 273, number of symmetry op =  8
  '* 57:-cba   D2h^11        Pmab                        -P 2b 2a', &               
  'x,y,z','-x,-y+1/2,z','x+1/2,-y,-z','-x+1/2,y+1/2,-z','-x,-y,-z','x,y+1/2,-z','-x+1/2,y,z','x+1/2,-y+1/2,z' /                     

 data spacegroupdata(1437:1445)/ &                              ! group index = 274, number of symmetry op =  8
  '* 57:bca    D2h^11        Pbma                        -P 2a 2ab', &              
  'x,y,z','-x+1/2,-y,z','x+1/2,-y+1/2,-z','-x,y+1/2,-z','-x,-y,-z','x+1/2,y,-z','-x+1/2,y+1/2,z','x,-y+1/2,z' /                     

 data spacegroupdata(1446:1454)/ &                              ! group index = 275, number of symmetry op =  8
  '* 57:a-cb   D2h^11        Pcmb                        -P 2bc 2c', &              
  'x,y,z','-x,-y+1/2,z+1/2','x,-y,-z+1/2','-x,y+1/2,-z','-x,-y,-z','x,y+1/2,-z+1/2','-x,y,z+1/2','x,-y+1/2,z' /                     

 data spacegroupdata(1455:1463)/ &                              ! group index = 276, number of symmetry op =  8
  '* 58        D2h^12        Pnnm                        -P 2 2n', &                
  'x,y,z','-x,-y,z','x+1/2,-y+1/2,-z+1/2','-x+1/2,y+1/2,-z+1/2','-x,-y,-z','x,y,-z','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(1464:1472)/ &                              ! group index = 277, number of symmetry op =  8
  '* 58:cab    D2h^12        Pmnn                        -P 2n 2', &                
  'x,y,z','-x+1/2,-y+1/2,z+1/2','x,-y,-z','-x+1/2,y+1/2,-z+1/2','-x,-y,-z','x+1/2,y+1/2,-z+1/2','-x,y,z','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(1473:1481)/ &                              ! group index = 278, number of symmetry op =  8
  '* 58:bca    D2h^12        Pnmn                        -P 2n 2n', &               
  'x,y,z','-x+1/2,-y+1/2,z+1/2','x+1/2,-y+1/2,-z+1/2','-x,y,-z','-x,-y,-z','x+1/2,y+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2','x,-y,z' /     

 data spacegroupdata(1482:1490)/ &                              ! group index = 279, number of symmetry op =  8
  '* 59:1      D2h^13        Pmmn:1                       P 2 2ab -1ab', &          
  'x,y,z','-x+1/2,-y+1/2,-z','-x,-y,z','x+1/2,-y+1/2,-z','-x+1/2,y+1/2,-z','x+1/2,y+1/2,-z','-x,y,z','x,-y,z' /                     

 data spacegroupdata(1491:1499)/ &                              ! group index = 280, number of symmetry op =  8
  '* 59:2      D2h^13        Pmmn:2                      -P 2ab 2a', &              
  'x,y,z','-x+1/2,-y+1/2,z','x+1/2,-y,-z','-x,y+1/2,-z','-x,-y,-z','x+1/2,y+1/2,-z','-x+1/2,y,z','x,-y+1/2,z' /                     

 data spacegroupdata(1500:1508)/ &                              ! group index = 281, number of symmetry op =  8
  '* 59:1cab   D2h^13        Pnmm:1                       P 2bc 2 -1bc', &          
  'x,y,z','-x,-y+1/2,-z+1/2','-x,-y+1/2,z+1/2','x,-y,-z','-x,y+1/2,-z+1/2','x,y,-z','-x,y+1/2,z+1/2','x,-y,z' /                     

 data spacegroupdata(1509:1517)/ &                              ! group index = 282, number of symmetry op =  8
  '* 59:2cab   D2h^13        Pnmm:2                      -P 2c 2bc', &              
  'x,y,z','-x,-y,z+1/2','x,-y+1/2,-z+1/2','-x,y+1/2,-z','-x,-y,-z','x,y,-z+1/2','-x,y+1/2,z+1/2','x,-y+1/2,z' /                     

 data spacegroupdata(1518:1526)/ &                              ! group index = 283, number of symmetry op =  8
  '* 59:1bca   D2h^13        Pmnm:1                       P 2ac 2ac -1ac', &        
  'x,y,z','-x+1/2,-y,-z+1/2','-x+1/2,-y,z+1/2','x+1/2,-y,-z+1/2','-x,y,-z','x,y,-z','-x,y,z','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1527:1535)/ &                              ! group index = 284, number of symmetry op =  8
  '* 59:2bca   D2h^13        Pmnm:2                      -P 2c 2a', &               
  'x,y,z','-x,-y,z+1/2','x+1/2,-y,-z','-x+1/2,y,-z+1/2','-x,-y,-z','x,y,-z+1/2','-x+1/2,y,z','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1536:1544)/ &                              ! group index = 285, number of symmetry op =  8
  '* 60        D2h^14        Pbcn                        -P 2n 2ab', &              
  'x,y,z','-x+1/2,-y+1/2,z+1/2','x+1/2,-y+1/2,-z','-x,y,-z+1/2','-x,-y,-z','x+1/2,y+1/2,-z+1/2','-x+1/2,y+1/2,z','x,-y,z+1/2' /     

 data spacegroupdata(1545:1553)/ &                              ! group index = 286, number of symmetry op =  8
  '* 60:ba-c   D2h^14        Pcan                        -P 2n 2c', &               
  'x,y,z','-x+1/2,-y+1/2,z+1/2','x,-y,-z+1/2','-x+1/2,y+1/2,-z','-x,-y,-z','x+1/2,y+1/2,-z+1/2','-x,y,z+1/2','x+1/2,-y+1/2,z' /     

 data spacegroupdata(1554:1562)/ &                              ! group index = 287, number of symmetry op =  8
  '* 60:cab    D2h^14        Pnca                        -P 2a 2n', &               
  'x,y,z','-x+1/2,-y,z','x+1/2,-y+1/2,-z+1/2','-x,y+1/2,-z+1/2','-x,-y,-z','x+1/2,y,-z','-x+1/2,y+1/2,z+1/2','x,-y+1/2,z+1/2' /     

 data spacegroupdata(1563:1571)/ &                              ! group index = 288, number of symmetry op =  8
  '* 60:-cba   D2h^14        Pnab                        -P 2bc 2n', &              
  'x,y,z','-x,-y+1/2,z+1/2','x+1/2,-y+1/2,-z+1/2','-x+1/2,y,-z','-x,-y,-z','x,y+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y,z' /     

 data spacegroupdata(1572:1580)/ &                              ! group index = 289, number of symmetry op =  8
  '* 60:bca    D2h^14        Pbna                        -P 2ac 2b', &              
  'x,y,z','-x+1/2,-y,z+1/2','x,-y+1/2,-z','-x+1/2,y+1/2,-z+1/2','-x,-y,-z','x+1/2,y,-z+1/2','-x,y+1/2,z','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(1581:1589)/ &                              ! group index = 290, number of symmetry op =  8
  '* 60:a-cb   D2h^14        Pcnb                        -P 2b 2ac', &              
  'x,y,z','-x,-y+1/2,z','x+1/2,-y,-z+1/2','-x+1/2,y+1/2,-z+1/2','-x,-y,-z','x,y+1/2,-z','-x+1/2,y,z+1/2','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(1590:1598)/ &                              ! group index = 291, number of symmetry op =  8
  '* 61        D2h^15        Pbca                        -P 2ac 2ab', &             
  'x,y,z','-x+1/2,-y,z+1/2','x+1/2,-y+1/2,-z','-x,y+1/2,-z+1/2','-x,-y,-z','x+1/2,y,-z+1/2','-x+1/2,y+1/2,z','x,-y+1/2,z+1/2' /     

 data spacegroupdata(1599:1607)/ &                              ! group index = 292, number of symmetry op =  8
  '* 61:ba-c   D2h^15        Pcab                        -P 2bc 2ac', &             
  'x,y,z','-x,-y+1/2,z+1/2','x+1/2,-y,-z+1/2','-x+1/2,y+1/2,-z','-x,-y,-z','x,y+1/2,-z+1/2','-x+1/2,y,z+1/2','x+1/2,-y+1/2,z' /     

 data spacegroupdata(1608:1616)/ &                              ! group index = 293, number of symmetry op =  8
  '* 62        D2h^16        Pnma                        -P 2ac 2n', &              
  'x,y,z','-x+1/2,-y,z+1/2','x+1/2,-y+1/2,-z+1/2','-x,y+1/2,-z','-x,-y,-z','x+1/2,y,-z+1/2','-x+1/2,y+1/2,z+1/2','x,-y+1/2,z' /     

 data spacegroupdata(1617:1625)/ &                              ! group index = 294, number of symmetry op =  8
  '* 62:ba-c   D2h^16        Pmnb                        -P 2bc 2a', &              
  'x,y,z','-x,-y+1/2,z+1/2','x+1/2,-y,-z','-x+1/2,y+1/2,-z+1/2','-x,-y,-z','x,y+1/2,-z+1/2','-x+1/2,y,z','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(1626:1634)/ &                              ! group index = 295, number of symmetry op =  8
  '* 62:cab    D2h^16        Pbnm                        -P 2c 2ab', &              
  'x,y,z','-x,-y,z+1/2','x+1/2,-y+1/2,-z','-x+1/2,y+1/2,-z+1/2','-x,-y,-z','x,y,-z+1/2','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(1635:1643)/ &                              ! group index = 296, number of symmetry op =  8
  '* 62:-cba   D2h^16        Pcmn                        -P 2n 2ac', &              
  'x,y,z','-x+1/2,-y+1/2,z+1/2','x+1/2,-y,-z+1/2','-x,y+1/2,-z','-x,-y,-z','x+1/2,y+1/2,-z+1/2','-x+1/2,y,z+1/2','x,-y+1/2,z' /     

 data spacegroupdata(1644:1652)/ &                              ! group index = 297, number of symmetry op =  8
  '* 62:bca    D2h^16        Pmcn                        -P 2n 2a', &               
  'x,y,z','-x+1/2,-y+1/2,z+1/2','x+1/2,-y,-z','-x,y+1/2,-z+1/2','-x,-y,-z','x+1/2,y+1/2,-z+1/2','-x+1/2,y,z','x,-y+1/2,z+1/2' /     

 data spacegroupdata(1653:1661)/ &                              ! group index = 298, number of symmetry op =  8
  '* 62:a-cb   D2h^16        Pnam                        -P 2c 2n', &               
  'x,y,z','-x,-y,z+1/2','x+1/2,-y+1/2,-z+1/2','-x+1/2,y+1/2,-z','-x,-y,-z','x,y,-z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z' /     

 data spacegroupdata(1662:1670)/ &                              ! group index = 299, number of symmetry op =  8
  '* 63        D2h^17        Cmcm                        -C 2c 2', &                
  'x,y,z','-x,-y,z+1/2','x,-y,-z','-x,y,-z+1/2','-x,-y,-z','x,y,-z+1/2','-x,y,z','x,-y,z+1/2' /                                     

 data spacegroupdata(1671:1679)/ &                              ! group index = 300, number of symmetry op =  8
  '* 63:ba-c   D2h^17        Ccmm                        -C 2c 2c', &               
  'x,y,z','-x,-y,z+1/2','x,-y,-z+1/2','-x,y,-z','-x,-y,-z','x,y,-z+1/2','-x,y,z+1/2','x,-y,z' /                                     

 data spacegroupdata(1680:1688)/ &                              ! group index = 301, number of symmetry op =  8
  '* 63:cab    D2h^17        Amma                        -A 2a 2a', &               
  'x,y,z','-x+1/2,-y,z','x+1/2,-y,-z','-x,y,-z','-x,-y,-z','x+1/2,y,-z','-x+1/2,y,z','x,-y,z' /                                     

 data spacegroupdata(1689:1697)/ &                              ! group index = 302, number of symmetry op =  8
  '* 63:-cba   D2h^17        Amam                        -A 2 2a', &                
  'x,y,z','-x,-y,z','x+1/2,-y,-z','-x+1/2,y,-z','-x,-y,-z','x,y,-z','-x+1/2,y,z','x+1/2,-y,z' /                                     

 data spacegroupdata(1698:1706)/ &                              ! group index = 303, number of symmetry op =  8
  '* 63:bca    D2h^17        Bbmm                        -B 2 2b', &                
  'x,y,z','-x,-y,z','x,-y+1/2,-z','-x,y+1/2,-z','-x,-y,-z','x,y,-z','-x,y+1/2,z','x,-y+1/2,z' /                                     

 data spacegroupdata(1707:1715)/ &                              ! group index = 304, number of symmetry op =  8
  '* 63:a-cb   D2h^17        Bmmb                        -B 2b 2', &                
  'x,y,z','-x,-y+1/2,z','x,-y,-z','-x,y+1/2,-z','-x,-y,-z','x,y+1/2,-z','-x,y,z','x,-y+1/2,z' /                                     

 data spacegroupdata(1716:1724)/ &                              ! group index = 305, number of symmetry op =  8
  '* 64        D2h^18        Cmca                        -C 2bc 2', &               
  'x,y,z','-x,-y+1/2,z+1/2','x,-y,-z','-x,y+1/2,-z+1/2','-x,-y,-z','x,y+1/2,-z+1/2','-x,y,z','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1725:1733)/ &                              ! group index = 306, number of symmetry op =  8
  '* 64:ba-c   D2h^18        Ccmb                        -C 2bc 2bc', &             
  'x,y,z','-x,-y+1/2,z+1/2','x,-y+1/2,-z+1/2','-x,y,-z','-x,-y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z+1/2','x,-y,z' /                     

 data spacegroupdata(1734:1742)/ &                              ! group index = 307, number of symmetry op =  8
  '* 64:cab    D2h^18        Abma                        -A 2ac 2ac', &             
  'x,y,z','-x+1/2,-y,z+1/2','x+1/2,-y,-z+1/2','-x,y,-z','-x,-y,-z','x+1/2,y,-z+1/2','-x+1/2,y,z+1/2','x,-y,z' /                     

 data spacegroupdata(1743:1751)/ &                              ! group index = 308, number of symmetry op =  8
  '* 64:-cba   D2h^18        Acam                        -A 2 2ac', &               
  'x,y,z','-x,-y,z','x+1/2,-y,-z+1/2','-x+1/2,y,-z+1/2','-x,-y,-z','x,y,-z','-x+1/2,y,z+1/2','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1752:1760)/ &                              ! group index = 309, number of symmetry op =  8
  '* 64:bca    D2h^18        Bbcm                        -B 2 2bc', &               
  'x,y,z','-x,-y,z','x,-y+1/2,-z+1/2','-x,y+1/2,-z+1/2','-x,-y,-z','x,y,-z','-x,y+1/2,z+1/2','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1761:1769)/ &                              ! group index = 310, number of symmetry op =  8
  '* 64:a-cb   D2h^18        Bmab                        -B 2bc 2', &               
  'x,y,z','-x,-y+1/2,z+1/2','x,-y,-z','-x,y+1/2,-z+1/2','-x,-y,-z','x,y+1/2,-z+1/2','-x,y,z','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1770:1778)/ &                              ! group index = 311, number of symmetry op =  8
  '* 65        D2h^19        Cmmm                        -C 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z','-x,-y,-z','x,y,-z','-x,y,z','x,-y,z' /                                                     

 data spacegroupdata(1779:1787)/ &                              ! group index = 312, number of symmetry op =  8
  '* 65:cab    D2h^19        Ammm                        -A 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z','-x,-y,-z','x,y,-z','-x,y,z','x,-y,z' /                                                     

 data spacegroupdata(1788:1796)/ &                              ! group index = 313, number of symmetry op =  8
  '* 65:bca    D2h^19        Bmmm                        -B 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z','-x,-y,-z','x,y,-z','-x,y,z','x,-y,z' /                                                     

 data spacegroupdata(1797:1805)/ &                              ! group index = 314, number of symmetry op =  8
  '* 66        D2h^20        Cccm                        -C 2 2c', &                
  'x,y,z','-x,-y,z','x,-y,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x,y,-z','-x,y,z+1/2','x,-y,z+1/2' /                                     

 data spacegroupdata(1806:1814)/ &                              ! group index = 315, number of symmetry op =  8
  '* 66:cab    D2h^20        Amaa                        -A 2a 2', &                
  'x,y,z','-x+1/2,-y,z','x,-y,-z','-x+1/2,y,-z','-x,-y,-z','x+1/2,y,-z','-x,y,z','x+1/2,-y,z' /                                     

 data spacegroupdata(1815:1823)/ &                              ! group index = 316, number of symmetry op =  8
  '* 66:bca    D2h^20        Bbmb                        -B 2b 2b', &               
  'x,y,z','-x,-y+1/2,z','x,-y+1/2,-z','-x,y,-z','-x,-y,-z','x,y+1/2,-z','-x,y+1/2,z','x,-y,z' /                                     

 data spacegroupdata(1824:1832)/ &                              ! group index = 317, number of symmetry op =  8
  '* 67        D2h^21        Cmma                        -C 2b 2', &                
  'x,y,z','-x,-y+1/2,z','x,-y,-z','-x,y+1/2,-z','-x,-y,-z','x,y+1/2,-z','-x,y,z','x,-y+1/2,z' /                                     

 data spacegroupdata(1833:1841)/ &                              ! group index = 318, number of symmetry op =  8
  '* 67:ba-c   D2h^21        Cmmb                        -C 2b 2b', &               
  'x,y,z','-x,-y+1/2,z','x,-y+1/2,-z','-x,y,-z','-x,-y,-z','x,y+1/2,-z','-x,y+1/2,z','x,-y,z' /                                     

 data spacegroupdata(1842:1850)/ &                              ! group index = 319, number of symmetry op =  8
  '* 67:cab    D2h^21        Abmm                        -A 2c 2c', &               
  'x,y,z','-x,-y,z+1/2','x,-y,-z+1/2','-x,y,-z','-x,-y,-z','x,y,-z+1/2','-x,y,z+1/2','x,-y,z' /                                     

 data spacegroupdata(1851:1859)/ &                              ! group index = 320, number of symmetry op =  8
  '* 67:-cba   D2h^21        Acmm                        -A 2 2c', &                
  'x,y,z','-x,-y,z','x,-y,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x,y,-z','-x,y,z+1/2','x,-y,z+1/2' /                                     

 data spacegroupdata(1860:1868)/ &                              ! group index = 321, number of symmetry op =  8
  '* 67:bca    D2h^21        Bmcm                        -B 2 2c', &                
  'x,y,z','-x,-y,z','x,-y,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x,y,-z','-x,y,z+1/2','x,-y,z+1/2' /                                     

 data spacegroupdata(1869:1877)/ &                              ! group index = 322, number of symmetry op =  8
  '* 67:a-cb   D2h^21        Bmam                        -B 2c 2', &                
  'x,y,z','-x,-y,z+1/2','x,-y,-z','-x,y,-z+1/2','-x,-y,-z','x,y,-z+1/2','-x,y,z','x,-y,z+1/2' /                                     

 data spacegroupdata(1878:1886)/ &                              ! group index = 323, number of symmetry op =  8
  '* 68:1      D2h^22        Ccca:1                       C 2 2 -1bc', &            
  'x,y,z','-x,-y+1/2,-z+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z+1/2','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1887:1895)/ &                              ! group index = 324, number of symmetry op =  8
  '* 68:2      D2h^22        Ccca:2                      -C 2b 2bc', &              
  'x,y,z','-x,-y+1/2,z','x,-y+1/2,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x,y+1/2,-z','-x,y+1/2,z+1/2','x,-y,z+1/2' /                     

 data spacegroupdata(1896:1904)/ &                              ! group index = 325, number of symmetry op =  8
  '* 68:1ba-c  D2h^22        Cccb:1                       C 2 2 -1bc', &            
  'x,y,z','-x,-y+1/2,-z+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z+1/2','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1905:1913)/ &                              ! group index = 326, number of symmetry op =  8
  '* 68:2ba-c  D2h^22        Cccb:2                      -C 2b 2c', &               
  'x,y,z','-x,-y+1/2,z','x,-y,-z+1/2','-x,y+1/2,-z+1/2','-x,-y,-z','x,y+1/2,-z','-x,y,z+1/2','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1914:1922)/ &                              ! group index = 327, number of symmetry op =  8
  '* 68:1cab   D2h^22        Abaa:1                       A 2 2 -1ac', &            
  'x,y,z','-x+1/2,-y,-z+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x+1/2,y,-z+1/2','-x+1/2,y,z+1/2','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1923:1931)/ &                              ! group index = 328, number of symmetry op =  8
  '* 68:2cab   D2h^22        Abaa:2                      -A 2a 2c', &               
  'x,y,z','-x+1/2,-y,z','x,-y,-z+1/2','-x+1/2,y,-z+1/2','-x,-y,-z','x+1/2,y,-z','-x,y,z+1/2','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1932:1940)/ &                              ! group index = 329, number of symmetry op =  8
  '* 68:1-cba  D2h^22        Acaa:1                       A 2 2 -1ac', &            
  'x,y,z','-x+1/2,-y,-z+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x+1/2,y,-z+1/2','-x+1/2,y,z+1/2','x+1/2,-y,z+1/2' /                     

 data spacegroupdata(1941:1949)/ &                              ! group index = 330, number of symmetry op =  8
  '* 68:2-cba  D2h^22        Acaa:2                      -A 2ac 2c', &              
  'x,y,z','-x+1/2,-y,z+1/2','x,-y,-z+1/2','-x+1/2,y,-z','-x,-y,-z','x+1/2,y,-z+1/2','-x,y,z+1/2','x+1/2,-y,z' /                     

 data spacegroupdata(1950:1958)/ &                              ! group index = 331, number of symmetry op =  8
  '* 68:1bca   D2h^22        Bbcb:1                       B 2 2 -1bc', &            
  'x,y,z','-x,-y+1/2,-z+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z+1/2','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1959:1967)/ &                              ! group index = 332, number of symmetry op =  8
  '* 68:2bca   D2h^22        Bbcb:2                      -B 2bc 2b', &              
  'x,y,z','-x,-y+1/2,z+1/2','x,-y+1/2,-z','-x,y,-z+1/2','-x,-y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z','x,-y,z+1/2' /                     

 data spacegroupdata(1968:1976)/ &                              ! group index = 333, number of symmetry op =  8
  '* 68:1a-cb  D2h^22        Bbab:1                       B 2 2 -1bc', &            
  'x,y,z','-x,-y+1/2,-z+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x,y+1/2,-z+1/2','-x,y+1/2,z+1/2','x,-y+1/2,z+1/2' /                     

 data spacegroupdata(1977:1985)/ &                              ! group index = 334, number of symmetry op =  8
  '* 68:2a-cb  D2h^22        Bbab:2                      -B 2b 2bc', &              
  'x,y,z','-x,-y+1/2,z','x,-y+1/2,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x,y+1/2,-z','-x,y+1/2,z+1/2','x,-y,z+1/2' /                     

 data spacegroupdata(1986:1994)/ &                              ! group index = 335, number of symmetry op =  8
  '* 69        D2h^23        Fmmm                        -F 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z','-x,-y,-z','x,y,-z','-x,y,z','x,-y,z' /                                                     

 data spacegroupdata(1995:2003)/ &                              ! group index = 336, number of symmetry op =  8
  '* 70:1      D2h^24        Fddd:1                       F 2 2 -1d', &             
  'x,y,z','-x+1/4,-y+1/4,-z+1/4','-x,-y,z','x,-y,-z','-x,y,-z','x+1/4,y+1/4,-z+1/4','-x+1/4,y+1/4,z+1/4','x+1/4,-y+1/4,z+1/4' /     

 data spacegroupdata(2004:2012)/ &                              ! group index = 337, number of symmetry op =  8
  '* 70:2      D2h^24        Fddd:2                      -F 2uv 2vw', &             
  'x,y,z','-x+1/4,-y+1/4,z','x,-y+1/4,-z+1/4','-x+1/4,y,-z+1/4','-x,-y,-z','x+3/4,y+3/4,-z','-x,y+3/4,z+3/4','x+3/4,-y,z+3/4' /     

 data spacegroupdata(2013:2021)/ &                              ! group index = 338, number of symmetry op =  8
  '* 71        D2h^25        Immm                        -I 2 2', &                 
  'x,y,z','-x,-y,z','x,-y,-z','-x,y,-z','-x,-y,-z','x,y,-z','-x,y,z','x,-y,z' /                                                     

 data spacegroupdata(2022:2030)/ &                              ! group index = 339, number of symmetry op =  8
  '* 72        D2h^26        Ibam                        -I 2 2c', &                
  'x,y,z','-x,-y,z','x,-y,-z+1/2','-x,y,-z+1/2','-x,-y,-z','x,y,-z','-x,y,z+1/2','x,-y,z+1/2' /                                     

 data spacegroupdata(2031:2039)/ &                              ! group index = 340, number of symmetry op =  8
  '* 72:cab    D2h^26        Imcb                        -I 2a 2', &                
  'x,y,z','-x+1/2,-y,z','x,-y,-z','-x+1/2,y,-z','-x,-y,-z','x+1/2,y,-z','-x,y,z','x+1/2,-y,z' /                                     

 data spacegroupdata(2040:2048)/ &                              ! group index = 341, number of symmetry op =  8
  '* 72:bca    D2h^26        Icma                        -I 2b 2b', &               
  'x,y,z','-x,-y+1/2,z','x,-y+1/2,-z','-x,y,-z','-x,-y,-z','x,y+1/2,-z','-x,y+1/2,z','x,-y,z' /                                     

 data spacegroupdata(2049:2057)/ &                              ! group index = 342, number of symmetry op =  8
  '* 73        D2h^27        Ibca                        -I 2b 2c', &               
  'x,y,z','-x,-y+1/2,z','x,-y,-z+1/2','-x+1/2,y,-z','-x,-y,-z','x,y+1/2,-z','-x,y,z+1/2','x+1/2,-y,z' /                             

 data spacegroupdata(2058:2066)/ &                              ! group index = 343, number of symmetry op =  8
  '* 73:ba-c   D2h^27        Icab                        -I 2a 2b', &               
  'x,y,z','-x+1/2,-y,z','x,-y+1/2,-z','-x,y,-z+1/2','-x,-y,-z','x+1/2,y,-z','-x,y+1/2,z','x,-y,z+1/2' /                             

 data spacegroupdata(2067:2075)/ &                              ! group index = 344, number of symmetry op =  8
  '* 74        D2h^28        Imma                        -I 2b 2', &                
  'x,y,z','-x,-y+1/2,z','x,-y,-z','-x,y+1/2,-z','-x,-y,-z','x,y+1/2,-z','-x,y,z','x,-y+1/2,z' /                                     

 data spacegroupdata(2076:2084)/ &                              ! group index = 345, number of symmetry op =  8
  '* 74:ba-c   D2h^28        Immb                        -I 2a 2a', &               
  'x,y,z','-x+1/2,-y,z','x+1/2,-y,-z','-x,y,-z','-x,-y,-z','x+1/2,y,-z','-x+1/2,y,z','x,-y,z' /                                     

 data spacegroupdata(2085:2093)/ &                              ! group index = 346, number of symmetry op =  8
  '* 74:cab    D2h^28        Ibmm                        -I 2c 2c', &               
  'x,y,z','-x,-y,z+1/2','x,-y,-z+1/2','-x,y,-z','-x,-y,-z','x,y,-z+1/2','-x,y,z+1/2','x,-y,z' /                                     

 data spacegroupdata(2094:2102)/ &                              ! group index = 347, number of symmetry op =  8
  '* 74:-cba   D2h^28        Icmm                        -I 2 2b', &                
  'x,y,z','-x,-y,z','x,-y+1/2,-z','-x,y+1/2,-z','-x,-y,-z','x,y,-z','-x,y+1/2,z','x,-y+1/2,z' /                                     

 data spacegroupdata(2103:2111)/ &                              ! group index = 348, number of symmetry op =  8
  '* 74:bca    D2h^28        Imcm                        -I 2 2a', &                
  'x,y,z','-x,-y,z','x+1/2,-y,-z','-x+1/2,y,-z','-x,-y,-z','x,y,-z','-x+1/2,y,z','x+1/2,-y,z' /                                     

 data spacegroupdata(2112:2120)/ &                              ! group index = 349, number of symmetry op =  8
  '* 74:a-cb   D2h^28        Imam                        -I 2c 2', &                
  'x,y,z','-x,-y,z+1/2','x,-y,-z','-x,y,-z+1/2','-x,-y,-z','x,y,-z+1/2','-x,y,z','x,-y,z+1/2' /                                     

 data spacegroupdata(2121:2125)/ &                              ! group index = 350, number of symmetry op =  4
  '* 75        C4^1          P4                           P 4', &                   
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z' /                                                                                             

 data spacegroupdata(2126:2130)/ &                              ! group index = 351, number of symmetry op =  4
  '* 76        C4^2          P41                          P 4w', &                  
  'x,y,z','-y,x,z+1/4','-x,-y,z+1/2','y,-x,z+3/4' /                                                                                 

 data spacegroupdata(2131:2135)/ &                              ! group index = 352, number of symmetry op =  4
  '* 77        C4^3          P42                          P 4c', &                  
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2' /                                                                                     

 data spacegroupdata(2136:2140)/ &                              ! group index = 353, number of symmetry op =  4
  '* 78        C4^4          P43                          P 4cw', &                 
  'x,y,z','-y,x,z+3/4','-x,-y,z+1/2','y,-x,z+1/4' /                                                                                 

 data spacegroupdata(2141:2145)/ &                              ! group index = 354, number of symmetry op =  4
  '* 79        C4^5          I4                           I 4', &                   
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z' /                                                                                             

 data spacegroupdata(2146:2150)/ &                              ! group index = 355, number of symmetry op =  4
  '* 80        C4^6          I41                          I 4bw', &                 
  'x,y,z','-y,x+1/2,z+1/4','-x,-y,z','y,-x+1/2,z+1/4' /                                                                             

 data spacegroupdata(2151:2155)/ &                              ! group index = 356, number of symmetry op =  4
  '* 81        S4^1          P-4                          P -4', &                  
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z' /                                                                                           

 data spacegroupdata(2156:2160)/ &                              ! group index = 357, number of symmetry op =  4
  '* 82        S4^2          I-4                          I -4', &                  
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z' /                                                                                           

 data spacegroupdata(2161:2169)/ &                              ! group index = 358, number of symmetry op =  8
  '* 83        C4h^1         P4/m                        -P 4', &                   
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','-x,-y,-z','y,-x,-z','x,y,-z','-y,x,-z' /                                                     

 data spacegroupdata(2170:2178)/ &                              ! group index = 359, number of symmetry op =  8
  '* 84        C4h^2         P42/m                       -P 4c', &                  
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2','-x,-y,-z','y,-x,-z+1/2','x,y,-z','-y,x,-z+1/2' /                                     

 data spacegroupdata(2179:2187)/ &                              ! group index = 360, number of symmetry op =  8
  '* 85:1      C4h^3         P4/n:1                       P 4ab -1ab', &            
  'x,y,z','-x+1/2,-y+1/2,-z','-y+1/2,x+1/2,z','-x,-y,z','y+1/2,-x+1/2,z','y,-x,-z','-y,x,-z','x+1/2,y+1/2,-z' /                     

 data spacegroupdata(2188:2196)/ &                              ! group index = 361, number of symmetry op =  8
  '* 85:2      C4h^3         P4/n:2                      -P 4a', &                  
  'x,y,z','-y+1/2,x,z','-x+1/2,-y+1/2,z','y,-x+1/2,z','-x,-y,-z','y+1/2,-x,-z','x+1/2,y+1/2,-z','-y,x+1/2,-z' /                     

 data spacegroupdata(2197:2205)/ &                              ! group index = 362, number of symmetry op =  8
  '* 86:1      C4h^4         P42/n:1                      P 4n -1n', &              
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','y,-x,-z','-y,x,-z','x+1/2,y+1/2,-z+1/2' /     

 data spacegroupdata(2206:2214)/ &                              ! group index = 363, number of symmetry op =  8
  '* 86:2      C4h^4         P42/n:2                     -P 4bc', &                 
  'x,y,z','-y,x+1/2,z+1/2','-x+1/2,-y+1/2,z','y+1/2,-x,z+1/2','-x,-y,-z','y,-x+1/2,-z+1/2','x+1/2,y+1/2,-z','-y+1/2,x,-z+1/2' /     

 data spacegroupdata(2215:2223)/ &                              ! group index = 364, number of symmetry op =  8
  '* 87        C4h^5         I4/m                        -I 4', &                   
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','-x,-y,-z','y,-x,-z','x,y,-z','-y,x,-z' /                                                     

 data spacegroupdata(2224:2232)/ &                              ! group index = 365, number of symmetry op =  8
  '* 88:1      C4h^6         I41/a:1                      I 4bw -1bw', &            
  'x,y,z','-x,-y+1/2,-z+1/4','-y,x+1/2,z+1/4','-x,-y,z','y,-x+1/2,z+1/4','y,-x,-z','-y,x,-z','x,y+1/2,-z+1/4' /                     

 data spacegroupdata(2233:2241)/ &                              ! group index = 366, number of symmetry op =  8
  '* 88:2      C4h^6         I41/a:2                     -I 4ad', &                 
  'x,y,z','-y+3/4,x+1/4,z+1/4','-x,-y+1/2,z','y+1/4,-x+1/4,z+1/4','-x,-y,-z','y+1/4,-x+3/4,-z+3/4','x,y+1/2,-z', &                  
  '-y+3/4,x+3/4,-z+3/4' /                                                                                                           

 data spacegroupdata(2242:2250)/ &                              ! group index = 367, number of symmetry op =  8
  '* 89        D4^1          P422                         P 4 2', &                 
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-y,-z','-x,y,-z','y,x,-z','-y,-x,-z' /                                                     

 data spacegroupdata(2251:2259)/ &                              ! group index = 368, number of symmetry op =  8
  '* 90        D4^2          P4212                        P 4ab 2ab', &             
  'x,y,z','-y+1/2,x+1/2,z','-x,-y,z','y+1/2,-x+1/2,z','x+1/2,-y+1/2,-z','-x+1/2,y+1/2,-z','y,x,-z','-y,-x,-z' /                     

 data spacegroupdata(2260:2268)/ &                              ! group index = 369, number of symmetry op =  8
  '* 91        D4^3          P4122                        P 4w 2c', &               
  'x,y,z','-y,x,z+1/4','-x,-y,z+1/2','y,-x,z+3/4','x,-y,-z+1/2','-x,y,-z','y,x,-z+3/4','-y,-x,-z+1/4' /                             

 data spacegroupdata(2269:2277)/ &                              ! group index = 370, number of symmetry op =  8
  '* 92        D4^4          P41212                       P 4abw 2nw', &            
  'x,y,z','-y+1/2,x+1/2,z+1/4','-x,-y,z+1/2','y+1/2,-x+1/2,z+3/4','x+1/2,-y+1/2,-z+3/4','-x+1/2,y+1/2,-z+1/4','y,x,-z', &           
  '-y,-x,-z+1/2' /                                                                                                                  

 data spacegroupdata(2278:2286)/ &                              ! group index = 371, number of symmetry op =  8
  '* 93        D4^5          P4222                        P 4c 2', &                
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2','x,-y,-z','-x,y,-z','y,x,-z+1/2','-y,-x,-z+1/2' /                                     

 data spacegroupdata(2287:2295)/ &                              ! group index = 372, number of symmetry op =  8
  '* 94        D4^6          P42212                       P 4n 2n', &               
  'x,y,z','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','x+1/2,-y+1/2,-z+1/2','-x+1/2,y+1/2,-z+1/2','y,x,-z','-y,-x,-z' /     

 data spacegroupdata(2296:2304)/ &                              ! group index = 373, number of symmetry op =  8
  '* 95        D4^7          P4322                        P 4cw 2c', &              
  'x,y,z','-y,x,z+3/4','-x,-y,z+1/2','y,-x,z+1/4','x,-y,-z+1/2','-x,y,-z','y,x,-z+1/4','-y,-x,-z+3/4' /                             

 data spacegroupdata(2305:2313)/ &                              ! group index = 374, number of symmetry op =  8
  '* 96        D4^8          P43212                       P 4nw 2abw', &            
  'x,y,z','-y+1/2,x+1/2,z+3/4','-x,-y,z+1/2','y+1/2,-x+1/2,z+1/4','x+1/2,-y+1/2,-z+1/4','-x+1/2,y+1/2,-z+3/4','y,x,-z', &           
  '-y,-x,-z+1/2' /                                                                                                                  

 data spacegroupdata(2314:2322)/ &                              ! group index = 375, number of symmetry op =  8
  '* 97        D4^9          I422                         I 4 2', &                 
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-y,-z','-x,y,-z','y,x,-z','-y,-x,-z' /                                                     

 data spacegroupdata(2323:2331)/ &                              ! group index = 376, number of symmetry op =  8
  '* 98        D4^10         I4122                        I 4bw 2bw', &             
  'x,y,z','-y,x+1/2,z+1/4','-x,-y,z','y,-x+1/2,z+1/4','x,-y+1/2,-z+1/4','-x,y+1/2,-z+1/4','y,x,-z','-y,-x,-z' /                     

 data spacegroupdata(2332:2340)/ &                              ! group index = 377, number of symmetry op =  8
  '* 99        C4v^1         P4mm                         P 4 -2', &                
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','-x,y,z','x,-y,z','-y,-x,z','y,x,z' /                                                         

 data spacegroupdata(2341:2349)/ &                              ! group index = 378, number of symmetry op =  8
  '*100        C4v^2         P4bm                         P 4 -2ab', &              
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z','-y+1/2,-x+1/2,z','y+1/2,x+1/2,z' /                         

 data spacegroupdata(2350:2358)/ &                              ! group index = 379, number of symmetry op =  8
  '*101        C4v^3         P42cm                        P 4c -2c', &              
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2','-x,y,z+1/2','x,-y,z+1/2','-y,-x,z','y,x,z' /                                         

 data spacegroupdata(2359:2367)/ &                              ! group index = 380, number of symmetry op =  8
  '*102        C4v^4         P42nm                        P 4n -2n', &              
  'x,y,z','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2','-y,-x,z','y,x,z' /         

 data spacegroupdata(2368:2376)/ &                              ! group index = 381, number of symmetry op =  8
  '*103        C4v^5         P4cc                         P 4 -2c', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','-x,y,z+1/2','x,-y,z+1/2','-y,-x,z+1/2','y,x,z+1/2' /                                         

 data spacegroupdata(2377:2385)/ &                              ! group index = 382, number of symmetry op =  8
  '*104        C4v^6         P4nc                         P 4 -2n', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2' /         

 data spacegroupdata(2386:2394)/ &                              ! group index = 383, number of symmetry op =  8
  '*105        C4v^7         P42mc                        P 4c -2', &               
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2','-x,y,z','x,-y,z','-y,-x,z+1/2','y,x,z+1/2' /                                         

 data spacegroupdata(2395:2403)/ &                              ! group index = 384, number of symmetry op =  8
  '*106        C4v^8         P42bc                        P 4c -2ab', &             
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2' /         

 data spacegroupdata(2404:2412)/ &                              ! group index = 385, number of symmetry op =  8
  '*107        C4v^9         I4mm                         I 4 -2', &                
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','-x,y,z','x,-y,z','-y,-x,z','y,x,z' /                                                         

 data spacegroupdata(2413:2421)/ &                              ! group index = 386, number of symmetry op =  8
  '*108        C4v^10        I4cm                         I 4 -2c', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','-x,y,z+1/2','x,-y,z+1/2','-y,-x,z+1/2','y,x,z+1/2' /                                         

 data spacegroupdata(2422:2430)/ &                              ! group index = 387, number of symmetry op =  8
  '*109        C4v^11        I41md                        I 4bw -2', &              
  'x,y,z','-y,x+1/2,z+1/4','-x,-y,z','y,-x+1/2,z+1/4','-x,y,z','x,-y,z','-y,-x+1/2,z+1/4','y,x+1/2,z+1/4' /                         

 data spacegroupdata(2431:2439)/ &                              ! group index = 388, number of symmetry op =  8
  '*110        C4v^12        I41cd                        I 4bw -2c', &             
  'x,y,z','-y,x+1/2,z+1/4','-x,-y,z','y,-x+1/2,z+1/4','-x,y,z+1/2','x,-y,z+1/2','-y+1/2,-x,z+1/4','y+1/2,x,z+1/4' /                 

 data spacegroupdata(2440:2448)/ &                              ! group index = 389, number of symmetry op =  8
  '*111        D2d^1         P-42m                        P -4 2', &                
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','x,-y,-z','-x,y,-z','-y,-x,z','y,x,z' /                                                     

 data spacegroupdata(2449:2457)/ &                              ! group index = 390, number of symmetry op =  8
  '*112        D2d^2         P-42c                        P -4 2c', &               
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','x,-y,-z+1/2','-x,y,-z+1/2','-y,-x,z+1/2','y,x,z+1/2' /                                     

 data spacegroupdata(2458:2466)/ &                              ! group index = 391, number of symmetry op =  8
  '*113        D2d^3         P-421m                       P -4 2ab', &              
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','x+1/2,-y+1/2,-z','-x+1/2,y+1/2,-z','-y+1/2,-x+1/2,z','y+1/2,x+1/2,z' /                     

 data spacegroupdata(2467:2475)/ &                              ! group index = 392, number of symmetry op =  8
  '*114        D2d^4         P-421c                       P -4 2n', &               
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','x+1/2,-y+1/2,-z+1/2','-x+1/2,y+1/2,-z+1/2','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2' /     

 data spacegroupdata(2476:2484)/ &                              ! group index = 393, number of symmetry op =  8
  '*115        D2d^5         P-4m2                        P -4 -2', &               
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','y,x,-z','-y,-x,-z','-x,y,z','x,-y,z' /                                                     

 data spacegroupdata(2485:2493)/ &                              ! group index = 394, number of symmetry op =  8
  '*116        D2d^6         P-4c2                        P -4 -2c', &              
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','y,x,-z+1/2','-y,-x,-z+1/2','-x,y,z+1/2','x,-y,z+1/2' /                                     

 data spacegroupdata(2494:2502)/ &                              ! group index = 395, number of symmetry op =  8
  '*117        D2d^7         P-4b2                        P -4 -2ab', &             
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','y+1/2,x+1/2,-z','-y+1/2,-x+1/2,-z','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z' /                     

 data spacegroupdata(2503:2511)/ &                              ! group index = 396, number of symmetry op =  8
  '*118        D2d^8         P-4n2                        P -4 -2n', &              
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','y+1/2,x+1/2,-z+1/2','-y+1/2,-x+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2' /     

 data spacegroupdata(2512:2520)/ &                              ! group index = 397, number of symmetry op =  8
  '*119        D2d^9         I-4m2                        I -4 -2', &               
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','y,x,-z','-y,-x,-z','-x,y,z','x,-y,z' /                                                     

 data spacegroupdata(2521:2529)/ &                              ! group index = 398, number of symmetry op =  8
  '*120        D2d^10        I-4c2                        I -4 -2c', &              
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','y,x,-z+1/2','-y,-x,-z+1/2','-x,y,z+1/2','x,-y,z+1/2' /                                     

 data spacegroupdata(2530:2538)/ &                              ! group index = 399, number of symmetry op =  8
  '*121        D2d^11        I-42m                        I -4 2', &                
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','x,-y,-z','-x,y,-z','-y,-x,z','y,x,z' /                                                     

 data spacegroupdata(2539:2547)/ &                              ! group index = 400, number of symmetry op =  8
  '*122        D2d^12        I-42d                        I -4 2bw', &              
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','x,-y+1/2,-z+1/4','-x,y+1/2,-z+1/4','-y,-x+1/2,z+1/4','y,x+1/2,z+1/4' /                     

 data spacegroupdata(2548:2564)/ &                              ! group index = 401, number of symmetry op = 16
  '*123        D4h^1         P4/mmm                      -P 4 2', &                 
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-y,-z','-x,y,-z','y,x,-z','-y,-x,-z','-x,-y,-z','y,-x,-z','x,y,-z','-y,x,-z','-x,y,z', &   
  'x,-y,z','-y,-x,z','y,x,z' /                                                                                                      

 data spacegroupdata(2565:2581)/ &                              ! group index = 402, number of symmetry op = 16
  '*124        D4h^2         P4/mcc                      -P 4 2c', &                
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-y,-z+1/2','-x,y,-z+1/2','y,x,-z+1/2','-y,-x,-z+1/2','-x,-y,-z','y,-x,-z','x,y,-z', &      
  '-y,x,-z','-x,y,z+1/2','x,-y,z+1/2','-y,-x,z+1/2','y,x,z+1/2' /                                                                   

 data spacegroupdata(2582:2598)/ &                              ! group index = 403, number of symmetry op = 16
  '*125:1      D4h^3         P4/nbm:1                     P 4 2 -1ab', &            
  'x,y,z','-x+1/2,-y+1/2,-z','-y,x,z','-x,-y,z','y,-x,z','y+1/2,-x+1/2,-z','-y+1/2,x+1/2,-z','x,-y,-z','-x,y,-z','y,x,-z', &        
  '-y,-x,-z','x+1/2,y+1/2,-z','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z','-y+1/2,-x+1/2,z','y+1/2,x+1/2,z' /                                 

 data spacegroupdata(2599:2615)/ &                              ! group index = 404, number of symmetry op = 16
  '*125:2      D4h^3         P4/nbm:2                    -P 4a 2b', &               
  'x,y,z','-y+1/2,x,z','-x+1/2,-y+1/2,z','y,-x+1/2,z','x,-y+1/2,-z','-x+1/2,y,-z','y,x,-z','-y+1/2,-x+1/2,-z','-x,-y,-z', &         
  'y+1/2,-x,-z','x+1/2,y+1/2,-z','-y,x+1/2,-z','-x,y+1/2,z','x+1/2,-y,z','-y,-x,z','y+1/2,x+1/2,z' /                                

 data spacegroupdata(2616:2632)/ &                              ! group index = 405, number of symmetry op = 16
  '*126:1      D4h^4         P4/nnc:1                     P 4 2 -1n', &             
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','-y,x,z','-x,-y,z','y,-x,z','y+1/2,-x+1/2,-z+1/2','-y+1/2,x+1/2,-z+1/2','x,-y,-z','-x,y,-z', &     
  'y,x,-z','-y,-x,-z','x+1/2,y+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2' /    

 data spacegroupdata(2633:2649)/ &                              ! group index = 406, number of symmetry op = 16
  '*126:2      D4h^4         P4/nnc:2                    -P 4a 2bc', &              
  'x,y,z','-y+1/2,x,z','-x+1/2,-y+1/2,z','y,-x+1/2,z','x,-y+1/2,-z+1/2','-x+1/2,y,-z+1/2','y,x,-z+1/2','-y+1/2,-x+1/2,-z+1/2', &    
  '-x,-y,-z','y+1/2,-x,-z','x+1/2,y+1/2,-z','-y,x+1/2,-z','-x,y+1/2,z+1/2','x+1/2,-y,z+1/2','-y,-x,z+1/2','y+1/2,x+1/2,z+1/2' /     

 data spacegroupdata(2650:2666)/ &                              ! group index = 407, number of symmetry op = 16
  '*127        D4h^5         P4/mbm                      -P 4 2ab', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x+1/2,-y+1/2,-z','-x+1/2,y+1/2,-z','y+1/2,x+1/2,-z','-y+1/2,-x+1/2,-z','-x,-y,-z', &         
  'y,-x,-z','x,y,-z','-y,x,-z','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z','-y+1/2,-x+1/2,z','y+1/2,x+1/2,z' /                                

 data spacegroupdata(2667:2683)/ &                              ! group index = 408, number of symmetry op = 16
  '*128        D4h^6         P4/mnc                      -P 4 2n', &                
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x+1/2,-y+1/2,-z+1/2','-x+1/2,y+1/2,-z+1/2','y+1/2,x+1/2,-z+1/2','-y+1/2,-x+1/2,-z+1/2', &    
  '-x,-y,-z','y,-x,-z','x,y,-z','-y,x,-z','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2' /     

 data spacegroupdata(2684:2700)/ &                              ! group index = 409, number of symmetry op = 16
  '*129:1      D4h^7         P4/nmm:1                     P 4ab 2ab -1ab', &        
  'x,y,z','-x+1/2,-y+1/2,-z','-y+1/2,x+1/2,z','-x,-y,z','y+1/2,-x+1/2,z','y,-x,-z','-y,x,-z','x+1/2,-y+1/2,-z','-x+1/2,y+1/2,-z', & 
  'y,x,-z','-y,-x,-z','x+1/2,y+1/2,-z','-x,y,z','x,-y,z','-y+1/2,-x+1/2,z','y+1/2,x+1/2,z' /                                        

 data spacegroupdata(2701:2717)/ &                              ! group index = 410, number of symmetry op = 16
  '*129:2      D4h^7         P4/nmm:2                    -P 4a 2a', &               
  'x,y,z','-y+1/2,x,z','-x+1/2,-y+1/2,z','y,-x+1/2,z','x+1/2,-y,-z','-x,y+1/2,-z','y+1/2,x+1/2,-z','-y,-x,-z','-x,-y,-z', &         
  'y+1/2,-x,-z','x+1/2,y+1/2,-z','-y,x+1/2,-z','-x+1/2,y,z','x,-y+1/2,z','-y+1/2,-x+1/2,z','y,x,z' /                                

 data spacegroupdata(2718:2734)/ &                              ! group index = 411, number of symmetry op = 16
  '*130:1      D4h^8         P4/ncc:1                     P 4ab 2n -1ab', &         
  'x,y,z','-x+1/2,-y+1/2,-z','-y+1/2,x+1/2,z','-x,-y,z','y+1/2,-x+1/2,z','y,-x,-z','-y,x,-z','x+1/2,-y+1/2,-z+1/2', &               
  '-x+1/2,y+1/2,-z+1/2','y,x,-z+1/2','-y,-x,-z+1/2','x+1/2,y+1/2,-z','-x,y,z+1/2','x,-y,z+1/2','-y+1/2,-x+1/2,z+1/2', &             
  'y+1/2,x+1/2,z+1/2' /                                                                                                             

 data spacegroupdata(2735:2751)/ &                              ! group index = 412, number of symmetry op = 16
  '*130:2      D4h^8         P4/ncc:2                    -P 4a 2ac', &              
  'x,y,z','-y+1/2,x,z','-x+1/2,-y+1/2,z','y,-x+1/2,z','x+1/2,-y,-z+1/2','-x,y+1/2,-z+1/2','y+1/2,x+1/2,-z+1/2','-y,-x,-z+1/2', &    
  '-x,-y,-z','y+1/2,-x,-z','x+1/2,y+1/2,-z','-y,x+1/2,-z','-x+1/2,y,z+1/2','x,-y+1/2,z+1/2','-y+1/2,-x+1/2,z+1/2','y,x,z+1/2' /     

 data spacegroupdata(2752:2768)/ &                              ! group index = 413, number of symmetry op = 16
  '*131        D4h^9         P42/mmc                     -P 4c 2', &                
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2','x,-y,-z','-x,y,-z','y,x,-z+1/2','-y,-x,-z+1/2','-x,-y,-z','y,-x,-z+1/2','x,y,-z', &  
  '-y,x,-z+1/2','-x,y,z','x,-y,z','-y,-x,z+1/2','y,x,z+1/2' /                                                                       

 data spacegroupdata(2769:2785)/ &                              ! group index = 414, number of symmetry op = 16
  '*132        D4h^10        P42/mcm                     -P 4c 2c', &               
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2','x,-y,-z+1/2','-x,y,-z+1/2','y,x,-z','-y,-x,-z','-x,-y,-z','y,-x,-z+1/2','x,y,-z', &  
  '-y,x,-z+1/2','-x,y,z+1/2','x,-y,z+1/2','-y,-x,z','y,x,z' /                                                                       

 data spacegroupdata(2786:2802)/ &                              ! group index = 415, number of symmetry op = 16
  '*133:1      D4h^11        P42/nbc:1                    P 4n 2c -1n', &           
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','y,-x,-z','-y,x,-z','x,-y,-z+1/2', &           
  '-x,y,-z+1/2','y+1/2,x+1/2,-z','-y+1/2,-x+1/2,-z','x+1/2,y+1/2,-z+1/2','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z','-y,-x,z+1/2', &         
  'y,x,z+1/2' /                                                                                                                     

 data spacegroupdata(2803:2819)/ &                              ! group index = 416, number of symmetry op = 16
  '*133:2      D4h^11        P42/nbc:2                   -P 4ac 2b', &              
  'x,y,z','-y+1/2,x,z+1/2','-x+1/2,-y+1/2,z','y,-x+1/2,z+1/2','x,-y+1/2,-z','-x+1/2,y,-z','y,x,-z+1/2','-y+1/2,-x+1/2,-z+1/2', &    
  '-x,-y,-z','y+1/2,-x,-z+1/2','x+1/2,y+1/2,-z','-y,x+1/2,-z+1/2','-x,y+1/2,z','x+1/2,-y,z','-y,-x,z+1/2','y+1/2,x+1/2,z+1/2' /     

 data spacegroupdata(2820:2836)/ &                              ! group index = 417, number of symmetry op = 16
  '*134:1      D4h^12        P42/nnm:1                    P 4n 2 -1n', &            
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','y,-x,-z','-y,x,-z','x,-y,-z','-x,y,-z', &     
  'y+1/2,x+1/2,-z+1/2','-y+1/2,-x+1/2,-z+1/2','x+1/2,y+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2','-y,-x,z','y,x,z' /    

 data spacegroupdata(2837:2853)/ &                              ! group index = 418, number of symmetry op = 16
  '*134:2      D4h^12        P42/nnm:2                   -P 4ac 2bc', &             
  'x,y,z','-y+1/2,x,z+1/2','-x+1/2,-y+1/2,z','y,-x+1/2,z+1/2','x,-y+1/2,-z+1/2','-x+1/2,y,-z+1/2','y,x,-z','-y+1/2,-x+1/2,-z', &    
  '-x,-y,-z','y+1/2,-x,-z+1/2','x+1/2,y+1/2,-z','-y,x+1/2,-z+1/2','-x,y+1/2,z+1/2','x+1/2,-y,z+1/2','-y,-x,z','y+1/2,x+1/2,z' /     

 data spacegroupdata(2854:2870)/ &                              ! group index = 419, number of symmetry op = 16
  '*135        D4h^13        P42/mbc                     -P 4c 2ab', &              
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2','x+1/2,-y+1/2,-z','-x+1/2,y+1/2,-z','y+1/2,x+1/2,-z+1/2','-y+1/2,-x+1/2,-z+1/2', &    
  '-x,-y,-z','y,-x,-z+1/2','x,y,-z','-y,x,-z+1/2','-x+1/2,y+1/2,z','x+1/2,-y+1/2,z','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2' /     

 data spacegroupdata(2871:2887)/ &                              ! group index = 420, number of symmetry op = 16
  '*136        D4h^14        P42/mnm                     -P 4n 2n', &               
  'x,y,z','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','x+1/2,-y+1/2,-z+1/2','-x+1/2,y+1/2,-z+1/2','y,x,-z','-y,-x,-z', &    
  '-x,-y,-z','y+1/2,-x+1/2,-z+1/2','x,y,-z','-y+1/2,x+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2','-y,-x,z','y,x,z' /     

 data spacegroupdata(2888:2904)/ &                              ! group index = 421, number of symmetry op = 16
  '*137:1      D4h^15        P42/nmc:1                    P 4n 2n -1n', &           
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','y,-x,-z','-y,x,-z','x+1/2,-y+1/2,-z+1/2', &   
  '-x+1/2,y+1/2,-z+1/2','y,x,-z','-y,-x,-z','x+1/2,y+1/2,-z+1/2','-x,y,z','x,-y,z','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2' /      

 data spacegroupdata(2905:2921)/ &                              ! group index = 422, number of symmetry op = 16
  '*137:2      D4h^15        P42/nmc:2                   -P 4ac 2a', &              
  'x,y,z','-y+1/2,x,z+1/2','-x+1/2,-y+1/2,z','y,-x+1/2,z+1/2','x+1/2,-y,-z','-x,y+1/2,-z','y+1/2,x+1/2,-z+1/2','-y,-x,-z+1/2', &    
  '-x,-y,-z','y+1/2,-x,-z+1/2','x+1/2,y+1/2,-z','-y,x+1/2,-z+1/2','-x+1/2,y,z','x,-y+1/2,z','-y+1/2,-x+1/2,z+1/2','y,x,z+1/2' /     

 data spacegroupdata(2922:2938)/ &                              ! group index = 423, number of symmetry op = 16
  '*138:1      D4h^16        P42/ncm:1                    P 4n 2ab -1n', &          
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','y,-x,-z','-y,x,-z','x+1/2,-y+1/2,-z', &       
  '-x+1/2,y+1/2,-z','y,x,-z+1/2','-y,-x,-z+1/2','x+1/2,y+1/2,-z+1/2','-x,y,z+1/2','x,-y,z+1/2','-y+1/2,-x+1/2,z','y+1/2,x+1/2,z' /  

 data spacegroupdata(2939:2955)/ &                              ! group index = 424, number of symmetry op = 16
  '*138:2      D4h^16        P42/ncm:2                   -P 4ac 2ac', &             
  'x,y,z','-y+1/2,x,z+1/2','-x+1/2,-y+1/2,z','y,-x+1/2,z+1/2','x+1/2,-y,-z+1/2','-x,y+1/2,-z+1/2','y+1/2,x+1/2,-z','-y,-x,-z', &    
  '-x,-y,-z','y+1/2,-x,-z+1/2','x+1/2,y+1/2,-z','-y,x+1/2,-z+1/2','-x+1/2,y,z+1/2','x,-y+1/2,z+1/2','-y+1/2,-x+1/2,z','y,x,z' /     

 data spacegroupdata(2956:2972)/ &                              ! group index = 425, number of symmetry op = 16
  '*139        D4h^17        I4/mmm                      -I 4 2', &                 
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-y,-z','-x,y,-z','y,x,-z','-y,-x,-z','-x,-y,-z','y,-x,-z','x,y,-z','-y,x,-z','-x,y,z', &   
  'x,-y,z','-y,-x,z','y,x,z' /                                                                                                      

 data spacegroupdata(2973:2989)/ &                              ! group index = 426, number of symmetry op = 16
  '*140        D4h^18        I4/mcm                      -I 4 2c', &                
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-y,-z+1/2','-x,y,-z+1/2','y,x,-z+1/2','-y,-x,-z+1/2','-x,-y,-z','y,-x,-z','x,y,-z', &      
  '-y,x,-z','-x,y,z+1/2','x,-y,z+1/2','-y,-x,z+1/2','y,x,z+1/2' /                                                                   

 data spacegroupdata(2990:3006)/ &                              ! group index = 427, number of symmetry op = 16
  '*141:1      D4h^19        I41/amd:1                    I 4bw 2bw -1bw', &        
  'x,y,z','-x,-y+1/2,-z+1/4','-y,x+1/2,z+1/4','-x,-y,z','y,-x+1/2,z+1/4','y,-x,-z','-y,x,-z','x,-y+1/2,-z+1/4','-x,y+1/2,-z+1/4', & 
  'y,x,-z','-y,-x,-z','x,y+1/2,-z+1/4','-x,y,z','x,-y,z','-y,-x+1/2,z+1/4','y,x+1/2,z+1/4' /                                        

 data spacegroupdata(3007:3023)/ &                              ! group index = 428, number of symmetry op = 16
  '*141:2      D4h^19        I41/amd:2                   -I 4bd 2', &               
  'x,y,z','-y+1/4,x+3/4,z+1/4','-x,-y+1/2,z','y+1/4,-x+1/4,z+3/4','x,-y,-z','-x,y+1/2,-z','y+1/4,x+3/4,-z+1/4', &                   
  '-y+1/4,-x+1/4,-z+3/4','-x,-y,-z','y+3/4,-x+1/4,-z+3/4','x,y+1/2,-z','-y+3/4,x+3/4,-z+1/4','-x,y,z','x,-y+1/2,z', &               
  '-y+3/4,-x+1/4,z+3/4','y+3/4,x+3/4,z+1/4' /                                                                                       

 data spacegroupdata(3024:3040)/ &                              ! group index = 429, number of symmetry op = 16
  '*142:1      D4h^20        I41/acd:1                    I 4bw 2aw -1bw', &        
  'x,y,z','-x,-y+1/2,-z+1/4','-y,x+1/2,z+1/4','-x,-y,z','y,-x+1/2,z+1/4','y,-x,-z','-y,x,-z','x+1/2,-y,-z+1/4','-x+1/2,y,-z+1/4', & 
  'y,x,-z+1/2','-y,-x,-z+1/2','x,y+1/2,-z+1/4','-x,y,z+1/2','x,-y,z+1/2','-y+1/2,-x,z+1/4','y+1/2,x,z+1/4' /                        

 data spacegroupdata(3041:3057)/ &                              ! group index = 430, number of symmetry op = 16
  '*142:2      D4h^20        I41/acd:2                   -I 4bd 2c', &              
  'x,y,z','-y+1/4,x+3/4,z+1/4','-x,-y+1/2,z','y+1/4,-x+1/4,z+3/4','x,-y,-z+1/2','-x+1/2,y,-z','y+3/4,x+1/4,-z+1/4', &               
  '-y+1/4,-x+1/4,-z+1/4','-x,-y,-z','y+3/4,-x+1/4,-z+3/4','x,y+1/2,-z','-y+3/4,x+3/4,-z+1/4','-x,y,z+1/2','x+1/2,-y,z', &           
  '-y+1/4,-x+3/4,z+3/4','y+3/4,x+3/4,z+3/4' /                                                                                       

 data spacegroupdata(3058:3061)/ &                              ! group index = 431, number of symmetry op =  3
  '*143        C3^1          P3                           P 3', &                   
  'x,y,z','-y,x-y,z','-x+y,-x,z' /                                                                                                  

 data spacegroupdata(3062:3065)/ &                              ! group index = 432, number of symmetry op =  3
  '*144        C3^2          P31                          P 31', &                  
  'x,y,z','-y,x-y,z+1/3','-x+y,-x,z+2/3' /                                                                                          

 data spacegroupdata(3066:3069)/ &                              ! group index = 433, number of symmetry op =  3
  '*145        C3^3          P32                          P 32', &                  
  'x,y,z','-y,x-y,z+2/3','-x+y,-x,z+1/3' /                                                                                          

 data spacegroupdata(3070:3073)/ &                              ! group index = 434, number of symmetry op =  3
  '*146:H      C3^4          R3:H                         R 3', &                   
  'x,y,z','-y,x-y,z','-x+y,-x,z' /                                                                                                  

 data spacegroupdata(3074:3077)/ &                              ! group index = 435, number of symmetry op =  3
  '*146:R      C3^4          R3:R                         P 3*', &                  
  'x,y,z','z,x,y','y,z,x' /                                                                                                         

 data spacegroupdata(3078:3084)/ &                              ! group index = 436, number of symmetry op =  6
  '*147        C3i^1         P-3                         -P 3', &                   
  'x,y,z','-y,x-y,z','-x+y,-x,z','-x,-y,-z','y,-x+y,-z','x-y,x,-z' /                                                                

 data spacegroupdata(3085:3091)/ &                              ! group index = 437, number of symmetry op =  6
  '*148:H      C3i^2         R-3:H                       -R 3', &                   
  'x,y,z','-y,x-y,z','-x+y,-x,z','-x,-y,-z','y,-x+y,-z','x-y,x,-z' /                                                                

 data spacegroupdata(3092:3098)/ &                              ! group index = 438, number of symmetry op =  6
  '*148:R      C3i^2         R-3:R                       -P 3*', &                  
  'x,y,z','z,x,y','y,z,x','-x,-y,-z','-z,-x,-y','-y,-z,-x' /                                                                        

 data spacegroupdata(3099:3105)/ &                              ! group index = 439, number of symmetry op =  6
  '*149        D3^1          P312                         P 3 2', &                 
  'x,y,z','-y,x-y,z','-x+y,-x,z','-y,-x,-z','-x+y,y,-z','x,x-y,-z' /                                                                

 data spacegroupdata(3106:3112)/ &                              ! group index = 440, number of symmetry op =  6
  '*150        D3^2          P321                         P 3 2"', &                
  'x,y,z','-y,x-y,z','-x+y,-x,z','x-y,-y,-z','-x,-x+y,-z','y,x,-z' /                                                                

 data spacegroupdata(3113:3119)/ &                              ! group index = 441, number of symmetry op =  6
  '*151        D3^3          P3112                        P 31 2c (0 0 1)', &       
  'x,y,z','-y,x-y,z+1/3','-x+y,-x,z+2/3','-y,-x,-z+2/3','-x+y,y,-z+1/3','x,x-y,-z' /                                                

 data spacegroupdata(3120:3126)/ &                              ! group index = 442, number of symmetry op =  6
  '*152        D3^4          P3121                        P 31 2"', &               
  'x,y,z','-y,x-y,z+1/3','-x+y,-x,z+2/3','x-y,-y,-z+2/3','-x,-x+y,-z+1/3','y,x,-z' /                                                

 data spacegroupdata(3127:3133)/ &                              ! group index = 443, number of symmetry op =  6
  '*153        D3^5          P3212                        P 32 2c (0 0 -1)', &      
  'x,y,z','-y,x-y,z+2/3','-x+y,-x,z+1/3','-y,-x,-z+1/3','-x+y,y,-z+2/3','x,x-y,-z' /                                                

 data spacegroupdata(3134:3140)/ &                              ! group index = 444, number of symmetry op =  6
  '*154        D3^6          P3221                        P 32 2"', &               
  'x,y,z','-y,x-y,z+2/3','-x+y,-x,z+1/3','x-y,-y,-z+1/3','-x,-x+y,-z+2/3','y,x,-z' /                                                

 data spacegroupdata(3141:3147)/ &                              ! group index = 445, number of symmetry op =  6
  '*155:H      D3^7          R32:H                        R 3 2"', &                
  'x,y,z','-y,x-y,z','-x+y,-x,z','x-y,-y,-z','-x,-x+y,-z','y,x,-z' /                                                                

 data spacegroupdata(3148:3154)/ &                              ! group index = 446, number of symmetry op =  6
  '*155:R      D3^7          R32:R                        P 3* 2', &                
  'x,y,z','z,x,y','y,z,x','-y,-x,-z','-x,-z,-y','-z,-y,-x' /                                                                        

 data spacegroupdata(3155:3161)/ &                              ! group index = 447, number of symmetry op =  6
  '*156        C3v^1         P3m1                         P 3 -2"', &               
  'x,y,z','-y,x-y,z','-x+y,-x,z','-x+y,y,z','x,x-y,z','-y,-x,z' /                                                                   

 data spacegroupdata(3162:3168)/ &                              ! group index = 448, number of symmetry op =  6
  '*157        C3v^2         P31m                         P 3 -2', &                
  'x,y,z','-y,x-y,z','-x+y,-x,z','y,x,z','x-y,-y,z','-x,-x+y,z' /                                                                   

 data spacegroupdata(3169:3175)/ &                              ! group index = 449, number of symmetry op =  6
  '*158        C3v^3         P3c1                         P 3 -2"c', &              
  'x,y,z','-y,x-y,z','-x+y,-x,z','-x+y,y,z+1/2','x,x-y,z+1/2','-y,-x,z+1/2' /                                                       

 data spacegroupdata(3176:3182)/ &                              ! group index = 450, number of symmetry op =  6
  '*159        C3v^4         P31c                         P 3 -2c', &               
  'x,y,z','-y,x-y,z','-x+y,-x,z','y,x,z+1/2','x-y,-y,z+1/2','-x,-x+y,z+1/2' /                                                       

 data spacegroupdata(3183:3189)/ &                              ! group index = 451, number of symmetry op =  6
  '*160:H      C3v^5         R3m:H                        R 3 -2"', &               
  'x,y,z','-y,x-y,z','-x+y,-x,z','-x+y,y,z','x,x-y,z','-y,-x,z' /                                                                   

 data spacegroupdata(3190:3196)/ &                              ! group index = 452, number of symmetry op =  6
  '*160:R      C3v^5         R3m:R                        P 3* -2', &               
  'x,y,z','z,x,y','y,z,x','y,x,z','x,z,y','z,y,x' /                                                                                 

 data spacegroupdata(3197:3203)/ &                              ! group index = 453, number of symmetry op =  6
  '*161:H      C3v^6         R3c:H                        R 3 -2"c', &              
  'x,y,z','-y,x-y,z','-x+y,-x,z','-x+y,y,z+1/2','x,x-y,z+1/2','-y,-x,z+1/2' /                                                       

 data spacegroupdata(3204:3210)/ &                              ! group index = 454, number of symmetry op =  6
  '*161:R      C3v^6         R3c:R                        P 3* -2n', &              
  'x,y,z','z,x,y','y,z,x','y+1/2,x+1/2,z+1/2','x+1/2,z+1/2,y+1/2','z+1/2,y+1/2,x+1/2' /                                             

 data spacegroupdata(3211:3223)/ &                              ! group index = 455, number of symmetry op = 12
  '*162        D3d^1         P-31m                       -P 3 2', &                 
  'x,y,z','-y,x-y,z','-x+y,-x,z','-y,-x,-z','-x+y,y,-z','x,x-y,-z','-x,-y,-z','y,-x+y,-z','x-y,x,-z','y,x,z','x-y,-y,z', &          
  '-x,-x+y,z' /                                                                                                                     

 data spacegroupdata(3224:3236)/ &                              ! group index = 456, number of symmetry op = 12
  '*163        D3d^2         P-31c                       -P 3 2c', &                
  'x,y,z','-y,x-y,z','-x+y,-x,z','-y,-x,-z+1/2','-x+y,y,-z+1/2','x,x-y,-z+1/2','-x,-y,-z','y,-x+y,-z','x-y,x,-z','y,x,z+1/2', &     
  'x-y,-y,z+1/2','-x,-x+y,z+1/2' /                                                                                                  

 data spacegroupdata(3237:3249)/ &                              ! group index = 457, number of symmetry op = 12
  '*164        D3d^3         P-3m1                       -P 3 2"', &                
  'x,y,z','-y,x-y,z','-x+y,-x,z','x-y,-y,-z','-x,-x+y,-z','y,x,-z','-x,-y,-z','y,-x+y,-z','x-y,x,-z','-x+y,y,z','x,x-y,z', &        
  '-y,-x,z' /                                                                                                                       

 data spacegroupdata(3250:3262)/ &                              ! group index = 458, number of symmetry op = 12
  '*165        D3d^4         P-3c1                       -P 3 2"c', &               
  'x,y,z','-y,x-y,z','-x+y,-x,z','x-y,-y,-z+1/2','-x,-x+y,-z+1/2','y,x,-z+1/2','-x,-y,-z','y,-x+y,-z','x-y,x,-z','-x+y,y,z+1/2', &  
  'x,x-y,z+1/2','-y,-x,z+1/2' /                                                                                                     

 data spacegroupdata(3263:3275)/ &                              ! group index = 459, number of symmetry op = 12
  '*166:H      D3d^5         R-3m:H                      -R 3 2"', &                
  'x,y,z','-y,x-y,z','-x+y,-x,z','x-y,-y,-z','-x,-x+y,-z','y,x,-z','-x,-y,-z','y,-x+y,-z','x-y,x,-z','-x+y,y,z','x,x-y,z', &        
  '-y,-x,z' /                                                                                                                       

 data spacegroupdata(3276:3288)/ &                              ! group index = 460, number of symmetry op = 12
  '*166:R      D3d^5         R-3m:R                      -P 3* 2', &                
  'x,y,z','z,x,y','y,z,x','-y,-x,-z','-x,-z,-y','-z,-y,-x','-x,-y,-z','-z,-x,-y','-y,-z,-x','y,x,z','x,z,y','z,y,x' /               

 data spacegroupdata(3289:3301)/ &                              ! group index = 461, number of symmetry op = 12
  '*167:H      D3d^6         R-3c:H                      -R 3 2"c', &               
  'x,y,z','-y,x-y,z','-x+y,-x,z','x-y,-y,-z+1/2','-x,-x+y,-z+1/2','y,x,-z+1/2','-x,-y,-z','y,-x+y,-z','x-y,x,-z','-x+y,y,z+1/2', &  
  'x,x-y,z+1/2','-y,-x,z+1/2' /                                                                                                     

 data spacegroupdata(3302:3314)/ &                              ! group index = 462, number of symmetry op = 12
  '*167:R      D3d^6         R-3c:R                      -P 3* 2n', &               
  'x,y,z','z,x,y','y,z,x','-y+1/2,-x+1/2,-z+1/2','-x+1/2,-z+1/2,-y+1/2','-z+1/2,-y+1/2,-x+1/2','-x,-y,-z','-z,-x,-y','-y,-z,-x', &  
  'y+1/2,x+1/2,z+1/2','x+1/2,z+1/2,y+1/2','z+1/2,y+1/2,x+1/2' /                                                                     

 data spacegroupdata(3315:3321)/ &                              ! group index = 463, number of symmetry op =  6
  '*168        C6^1          P6                           P 6', &                   
  'x,y,z','x-y,x,z','-y,x-y,z','-x,-y,z','-x+y,-x,z','y,-x+y,z' /                                                                   

 data spacegroupdata(3322:3328)/ &                              ! group index = 464, number of symmetry op =  6
  '*169        C6^2          P61                          P 61', &                  
  'x,y,z','x-y,x,z+1/6','-y,x-y,z+1/3','-x,-y,z+1/2','-x+y,-x,z+2/3','y,-x+y,z+5/6' /                                               

 data spacegroupdata(3329:3335)/ &                              ! group index = 465, number of symmetry op =  6
  '*170        C6^3          P65                          P 65', &                  
  'x,y,z','x-y,x,z+5/6','-y,x-y,z+2/3','-x,-y,z+1/2','-x+y,-x,z+1/3','y,-x+y,z+1/6' /                                               

 data spacegroupdata(3336:3342)/ &                              ! group index = 466, number of symmetry op =  6
  '*171        C6^4          P62                          P 62', &                  
  'x,y,z','x-y,x,z+1/3','-y,x-y,z+2/3','-x,-y,z','-x+y,-x,z+1/3','y,-x+y,z+2/3' /                                                   

 data spacegroupdata(3343:3349)/ &                              ! group index = 467, number of symmetry op =  6
  '*172        C6^5          P64                          P 64', &                  
  'x,y,z','x-y,x,z+2/3','-y,x-y,z+1/3','-x,-y,z','-x+y,-x,z+2/3','y,-x+y,z+1/3' /                                                   

 data spacegroupdata(3350:3356)/ &                              ! group index = 468, number of symmetry op =  6
  '*173        C6^6          P63                          P 6c', &                  
  'x,y,z','x-y,x,z+1/2','-y,x-y,z','-x,-y,z+1/2','-x+y,-x,z','y,-x+y,z+1/2' /                                                       

 data spacegroupdata(3357:3363)/ &                              ! group index = 469, number of symmetry op =  6
  '*174        C3h^1         P-6                          P -6', &                  
  'x,y,z','-x+y,-x,-z','-y,x-y,z','x,y,-z','-x+y,-x,z','-y,x-y,-z' /                                                                

 data spacegroupdata(3364:3376)/ &                              ! group index = 470, number of symmetry op = 12
  '*175        C6h^1         P6/m                        -P 6', &                   
  'x,y,z','x-y,x,z','-y,x-y,z','-x,-y,z','-x+y,-x,z','y,-x+y,z','-x,-y,-z','-x+y,-x,-z','y,-x+y,-z','x,y,-z','x-y,x,-z', &          
  '-y,x-y,-z' /                                                                                                                     

 data spacegroupdata(3377:3389)/ &                              ! group index = 471, number of symmetry op = 12
  '*176        C6h^2         P63/m                       -P 6c', &                  
  'x,y,z','x-y,x,z+1/2','-y,x-y,z','-x,-y,z+1/2','-x+y,-x,z','y,-x+y,z+1/2','-x,-y,-z','-x+y,-x,-z+1/2','y,-x+y,-z','x,y,-z+1/2', & 
  'x-y,x,-z','-y,x-y,-z+1/2' /                                                                                                      

 data spacegroupdata(3390:3402)/ &                              ! group index = 472, number of symmetry op = 12
  '*177        D6^1          P622                         P 6 2', &                 
  'x,y,z','x-y,x,z','-y,x-y,z','-x,-y,z','-x+y,-x,z','y,-x+y,z','x-y,-y,-z','-x,-x+y,-z','y,x,-z','-y,-x,-z','-x+y,y,-z', &         
  'x,x-y,-z' /                                                                                                                      

 data spacegroupdata(3403:3415)/ &                              ! group index = 473, number of symmetry op = 12
  '*178        D6^2          P6122                        P 61 2 (0 0 -1)', &       
  'x,y,z','x-y,x,z+1/6','-y,x-y,z+1/3','-x,-y,z+1/2','-x+y,-x,z+2/3','y,-x+y,z+5/6','x-y,-y,-z','-x,-x+y,-z+2/3','y,x,-z+1/3', &    
  '-y,-x,-z+5/6','-x+y,y,-z+1/2','x,x-y,-z+1/6' /                                                                                   

 data spacegroupdata(3416:3428)/ &                              ! group index = 474, number of symmetry op = 12
  '*179        D6^3          P6522                        P 65 2 (0 0 1)', &        
  'x,y,z','x-y,x,z+5/6','-y,x-y,z+2/3','-x,-y,z+1/2','-x+y,-x,z+1/3','y,-x+y,z+1/6','x-y,-y,-z','-x,-x+y,-z+1/3','y,x,-z+2/3', &    
  '-y,-x,-z+1/6','-x+y,y,-z+1/2','x,x-y,-z+5/6' /                                                                                   

 data spacegroupdata(3429:3441)/ &                              ! group index = 475, number of symmetry op = 12
  '*180        D6^4          P6222                        P 62 2c (0 0 1)', &       
  'x,y,z','x-y,x,z+1/3','-y,x-y,z+2/3','-x,-y,z','-x+y,-x,z+1/3','y,-x+y,z+2/3','x-y,-y,-z','-x,-x+y,-z+1/3','y,x,-z+2/3', &        
  '-y,-x,-z+2/3','-x+y,y,-z','x,x-y,-z+1/3' /                                                                                       

 data spacegroupdata(3442:3454)/ &                              ! group index = 476, number of symmetry op = 12
  '*181        D6^5          P6422                        P 64 2c (0 0 -1)', &      
  'x,y,z','x-y,x,z+2/3','-y,x-y,z+1/3','-x,-y,z','-x+y,-x,z+2/3','y,-x+y,z+1/3','x-y,-y,-z','-x,-x+y,-z+2/3','y,x,-z+1/3', &        
  '-y,-x,-z+1/3','-x+y,y,-z','x,x-y,-z+2/3' /                                                                                       

 data spacegroupdata(3455:3467)/ &                              ! group index = 477, number of symmetry op = 12
  '*182        D6^6          P6322                        P 6c 2c', &               
  'x,y,z','x-y,x,z+1/2','-y,x-y,z','-x,-y,z+1/2','-x+y,-x,z','y,-x+y,z+1/2','x-y,-y,-z','-x,-x+y,-z','y,x,-z','-y,-x,-z+1/2', &     
  '-x+y,y,-z+1/2','x,x-y,-z+1/2' /                                                                                                  

 data spacegroupdata(3468:3480)/ &                              ! group index = 478, number of symmetry op = 12
  '*183        C6v^1         P6mm                         P 6 -2', &                
  'x,y,z','x-y,x,z','-y,x-y,z','-x,-y,z','-x+y,-x,z','y,-x+y,z','-x+y,y,z','x,x-y,z','-y,-x,z','y,x,z','x-y,-y,z','-x,-x+y,z' /     

 data spacegroupdata(3481:3493)/ &                              ! group index = 479, number of symmetry op = 12
  '*184        C6v^2         P6cc                         P 6 -2c', &               
  'x,y,z','x-y,x,z','-y,x-y,z','-x,-y,z','-x+y,-x,z','y,-x+y,z','-x+y,y,z+1/2','x,x-y,z+1/2','-y,-x,z+1/2','y,x,z+1/2', &           
  'x-y,-y,z+1/2','-x,-x+y,z+1/2' /                                                                                                  

 data spacegroupdata(3494:3506)/ &                              ! group index = 480, number of symmetry op = 12
  '*185        C6v^3         P63cm                        P 6c -2', &               
  'x,y,z','x-y,x,z+1/2','-y,x-y,z','-x,-y,z+1/2','-x+y,-x,z','y,-x+y,z+1/2','-x+y,y,z+1/2','x,x-y,z+1/2','-y,-x,z+1/2','y,x,z', &   
  'x-y,-y,z','-x,-x+y,z' /                                                                                                          

 data spacegroupdata(3507:3519)/ &                              ! group index = 481, number of symmetry op = 12
  '*186        C6v^4         P63mc                        P 6c -2c', &              
  'x,y,z','x-y,x,z+1/2','-y,x-y,z','-x,-y,z+1/2','-x+y,-x,z','y,-x+y,z+1/2','-x+y,y,z','x,x-y,z','-y,-x,z','y,x,z+1/2', &           
  'x-y,-y,z+1/2','-x,-x+y,z+1/2' /                                                                                                  

 data spacegroupdata(3520:3532)/ &                              ! group index = 482, number of symmetry op = 12
  '*187        D3h^1         P-6m2                        P -6 2', &                
  'x,y,z','-x+y,-x,-z','-y,x-y,z','x,y,-z','-x+y,-x,z','-y,x-y,-z','-y,-x,-z','-x+y,y,-z','x,x-y,-z','-x+y,y,z','x,x-y,z', &        
  '-y,-x,z' /                                                                                                                       

 data spacegroupdata(3533:3545)/ &                              ! group index = 483, number of symmetry op = 12
  '*188        D3h^2         P-6c2                        P -6c 2', &               
  'x,y,z','-x+y,-x,-z+1/2','-y,x-y,z','x,y,-z+1/2','-x+y,-x,z','-y,x-y,-z+1/2','-y,-x,-z','-x+y,y,-z','x,x-y,-z','-x+y,y,z+1/2', &  
  'x,x-y,z+1/2','-y,-x,z+1/2' /                                                                                                     

 data spacegroupdata(3546:3558)/ &                              ! group index = 484, number of symmetry op = 12
  '*189        D3h^3         P-62m                        P -6 -2', &               
  'x,y,z','-x+y,-x,-z','-y,x-y,z','x,y,-z','-x+y,-x,z','-y,x-y,-z','x-y,-y,-z','-x,-x+y,-z','y,x,-z','y,x,z','x-y,-y,z', &          
  '-x,-x+y,z' /                                                                                                                     

 data spacegroupdata(3559:3571)/ &                              ! group index = 485, number of symmetry op = 12
  '*190        D3h^4         P-62c                        P -6c -2c', &             
  'x,y,z','-x+y,-x,-z+1/2','-y,x-y,z','x,y,-z+1/2','-x+y,-x,z','-y,x-y,-z+1/2','x-y,-y,-z','-x,-x+y,-z','y,x,-z','y,x,z+1/2', &     
  'x-y,-y,z+1/2','-x,-x+y,z+1/2' /                                                                                                  

 data spacegroupdata(3572:3596)/ &                              ! group index = 486, number of symmetry op = 24
  '*191        D6h^1         P6/mmm                      -P 6 2', &                 
  'x,y,z','x-y,x,z','-y,x-y,z','-x,-y,z','-x+y,-x,z','y,-x+y,z','x-y,-y,-z','-x,-x+y,-z','y,x,-z','-y,-x,-z','-x+y,y,-z', &         
  'x,x-y,-z','-x,-y,-z','-x+y,-x,-z','y,-x+y,-z','x,y,-z','x-y,x,-z','-y,x-y,-z','-x+y,y,z','x,x-y,z','-y,-x,z','y,x,z', &          
  'x-y,-y,z','-x,-x+y,z' /                                                                                                          

 data spacegroupdata(3597:3621)/ &                              ! group index = 487, number of symmetry op = 24
  '*192        D6h^2         P6/mcc                      -P 6 2c', &                
  'x,y,z','x-y,x,z','-y,x-y,z','-x,-y,z','-x+y,-x,z','y,-x+y,z','x-y,-y,-z+1/2','-x,-x+y,-z+1/2','y,x,-z+1/2','-y,-x,-z+1/2', &     
  '-x+y,y,-z+1/2','x,x-y,-z+1/2','-x,-y,-z','-x+y,-x,-z','y,-x+y,-z','x,y,-z','x-y,x,-z','-y,x-y,-z','-x+y,y,z+1/2','x,x-y,z+1/2', &
  '-y,-x,z+1/2','y,x,z+1/2','x-y,-y,z+1/2','-x,-x+y,z+1/2' /                                                                        

 data spacegroupdata(3622:3646)/ &                              ! group index = 488, number of symmetry op = 24
  '*193        D6h^3         P63/mcm                     -P 6c 2', &                
  'x,y,z','x-y,x,z+1/2','-y,x-y,z','-x,-y,z+1/2','-x+y,-x,z','y,-x+y,z+1/2','x-y,-y,-z+1/2','-x,-x+y,-z+1/2','y,x,-z+1/2', &        
  '-y,-x,-z','-x+y,y,-z','x,x-y,-z','-x,-y,-z','-x+y,-x,-z+1/2','y,-x+y,-z','x,y,-z+1/2','x-y,x,-z','-y,x-y,-z+1/2', &              
  '-x+y,y,z+1/2','x,x-y,z+1/2','-y,-x,z+1/2','y,x,z','x-y,-y,z','-x,-x+y,z' /                                                       

 data spacegroupdata(3647:3671)/ &                              ! group index = 489, number of symmetry op = 24
  '*194        D6h^4         P63/mmc                     -P 6c 2c', &               
  'x,y,z','x-y,x,z+1/2','-y,x-y,z','-x,-y,z+1/2','-x+y,-x,z','y,-x+y,z+1/2','x-y,-y,-z','-x,-x+y,-z','y,x,-z','-y,-x,-z+1/2', &     
  '-x+y,y,-z+1/2','x,x-y,-z+1/2','-x,-y,-z','-x+y,-x,-z+1/2','y,-x+y,-z','x,y,-z+1/2','x-y,x,-z','-y,x-y,-z+1/2','-x+y,y,z', &      
  'x,x-y,z','-y,-x,z','y,x,z+1/2','x-y,-y,z+1/2','-x,-x+y,z+1/2' /                                                                  

 data spacegroupdata(3672:3684)/ &                              ! group index = 490, number of symmetry op = 12
  '*195        T^1           P23                          P 2 2 3', &               
  'x,y,z','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-x,-y,z','x,-y,-z','-x,y,-z' /               

 data spacegroupdata(3685:3697)/ &                              ! group index = 491, number of symmetry op = 12
  '*196        T^2           F23                          F 2 2 3', &               
  'x,y,z','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-x,-y,z','x,-y,-z','-x,y,-z' /               

 data spacegroupdata(3698:3710)/ &                              ! group index = 492, number of symmetry op = 12
  '*197        T^3           I23                          I 2 2 3', &               
  'x,y,z','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-x,-y,z','x,-y,-z','-x,y,-z' /               

 data spacegroupdata(3711:3723)/ &                              ! group index = 493, number of symmetry op = 12
  '*198        T^4           P213                         P 2ac 2ab 3', &           
  'x,y,z','z,x,y','y,z,x','-y+1/2,-z,x+1/2','z+1/2,-x+1/2,-y','-y,z+1/2,-x+1/2','-z+1/2,-x,y+1/2','-z,x+1/2,-y+1/2', &              
  'y+1/2,-z+1/2,-x','-x+1/2,-y,z+1/2','x+1/2,-y+1/2,-z','-x,y+1/2,-z+1/2' /                                                         

 data spacegroupdata(3724:3736)/ &                              ! group index = 494, number of symmetry op = 12
  '*199        T^5           I213                         I 2b 2c 3', &             
  'x,y,z','z,x,y','y,z,x','-y,-z+1/2,x','z,-x,-y+1/2','-y+1/2,z,-x','-z,-x+1/2,y','-z+1/2,x,-y','y,-z,-x+1/2','-x,-y+1/2,z', &      
  'x,-y,-z+1/2','-x+1/2,y,-z' /                                                                                                     

 data spacegroupdata(3737:3761)/ &                              ! group index = 495, number of symmetry op = 24
  '*200        Th^1          Pm-3                        -P 2 2 3', &               
  'x,y,z','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-x,-y,z','x,-y,-z','-x,y,-z','-x,-y,-z', &   
  '-z,-x,-y','-y,-z,-x','y,z,-x','-z,x,y','y,-z,x','z,x,-y','z,-x,y','-y,z,x','x,y,-z','-x,y,z','x,-y,z' /                          

 data spacegroupdata(3762:3786)/ &                              ! group index = 496, number of symmetry op = 24
  '*201:1      Th^2          Pn-3:1                       P 2 2 3 -1n', &           
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x', &                     
  '-z+1/2,-x+1/2,-y+1/2','-y+1/2,-z+1/2,-x+1/2','y+1/2,z+1/2,-x+1/2','-z+1/2,x+1/2,y+1/2','y+1/2,-z+1/2,x+1/2', &                   
  'z+1/2,x+1/2,-y+1/2','z+1/2,-x+1/2,y+1/2','-y+1/2,z+1/2,x+1/2','-x,-y,z','x,-y,-z','-x,y,-z','x+1/2,y+1/2,-z+1/2', &              
  '-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2' /                                                                                       

 data spacegroupdata(3787:3811)/ &                              ! group index = 497, number of symmetry op = 24
  '*201:2      Th^2          Pn-3:2                      -P 2ab 2bc 3', &           
  'x,y,z','z,x,y','y,z,x','-y+1/2,-z+1/2,x','z,-x+1/2,-y+1/2','-y+1/2,z,-x+1/2','-z+1/2,-x+1/2,y','-z+1/2,x,-y+1/2', &              
  'y,-z+1/2,-x+1/2','-x+1/2,-y+1/2,z','x,-y+1/2,-z+1/2','-x+1/2,y,-z+1/2','-x,-y,-z','-z,-x,-y','-y,-z,-x','y+1/2,z+1/2,-x', &      
  '-z,x+1/2,y+1/2','y+1/2,-z,x+1/2','z+1/2,x+1/2,-y','z+1/2,-x,y+1/2','-y,z+1/2,x+1/2','x+1/2,y+1/2,-z','-x,y+1/2,z+1/2', &         
  'x+1/2,-y,z+1/2' /                                                                                                                

 data spacegroupdata(3812:3836)/ &                              ! group index = 498, number of symmetry op = 24
  '*202        Th^3          Fm-3                        -F 2 2 3', &               
  'x,y,z','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-x,-y,z','x,-y,-z','-x,y,-z','-x,-y,-z', &   
  '-z,-x,-y','-y,-z,-x','y,z,-x','-z,x,y','y,-z,x','z,x,-y','z,-x,y','-y,z,x','x,y,-z','-x,y,z','x,-y,z' /                          

 data spacegroupdata(3837:3861)/ &                              ! group index = 499, number of symmetry op = 24
  '*203:1      Th^4          Fd-3:1                       F 2 2 3 -1d', &           
  'x,y,z','-x+1/4,-y+1/4,-z+1/4','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x', &                     
  '-z+1/4,-x+1/4,-y+1/4','-y+1/4,-z+1/4,-x+1/4','y+1/4,z+1/4,-x+1/4','-z+1/4,x+1/4,y+1/4','y+1/4,-z+1/4,x+1/4', &                   
  'z+1/4,x+1/4,-y+1/4','z+1/4,-x+1/4,y+1/4','-y+1/4,z+1/4,x+1/4','-x,-y,z','x,-y,-z','-x,y,-z','x+1/4,y+1/4,-z+1/4', &              
  '-x+1/4,y+1/4,z+1/4','x+1/4,-y+1/4,z+1/4' /                                                                                       

 data spacegroupdata(3862:3886)/ &                              ! group index = 500, number of symmetry op = 24
  '*203:2      Th^4          Fd-3:2                      -F 2uv 2vw 3', &           
  'x,y,z','z,x,y','y,z,x','-y+1/4,-z+1/4,x','z,-x+1/4,-y+1/4','-y+1/4,z,-x+1/4','-z+1/4,-x+1/4,y','-z+1/4,x,-y+1/4', &              
  'y,-z+1/4,-x+1/4','-x+1/4,-y+1/4,z','x,-y+1/4,-z+1/4','-x+1/4,y,-z+1/4','-x,-y,-z','-z,-x,-y','-y,-z,-x','y+3/4,z+3/4,-x', &      
  '-z,x+3/4,y+3/4','y+3/4,-z,x+3/4','z+3/4,x+3/4,-y','z+3/4,-x,y+3/4','-y,z+3/4,x+3/4','x+3/4,y+3/4,-z','-x,y+3/4,z+3/4', &         
  'x+3/4,-y,z+3/4' /                                                                                                                

 data spacegroupdata(3887:3911)/ &                              ! group index = 501, number of symmetry op = 24
  '*204        Th^5          Im-3                        -I 2 2 3', &               
  'x,y,z','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-x,-y,z','x,-y,-z','-x,y,-z','-x,-y,-z', &   
  '-z,-x,-y','-y,-z,-x','y,z,-x','-z,x,y','y,-z,x','z,x,-y','z,-x,y','-y,z,x','x,y,-z','-x,y,z','x,-y,z' /                          

 data spacegroupdata(3912:3936)/ &                              ! group index = 502, number of symmetry op = 24
  '*205        Th^6          Pa-3                        -P 2ac 2ab 3', &           
  'x,y,z','z,x,y','y,z,x','-y+1/2,-z,x+1/2','z+1/2,-x+1/2,-y','-y,z+1/2,-x+1/2','-z+1/2,-x,y+1/2','-z,x+1/2,-y+1/2', &              
  'y+1/2,-z+1/2,-x','-x+1/2,-y,z+1/2','x+1/2,-y+1/2,-z','-x,y+1/2,-z+1/2','-x,-y,-z','-z,-x,-y','-y,-z,-x','y+1/2,z,-x+1/2', &      
  '-z+1/2,x+1/2,y','y,-z+1/2,x+1/2','z+1/2,x,-y+1/2','z,-x+1/2,y+1/2','-y+1/2,z+1/2,x','x+1/2,y,-z+1/2','-x+1/2,y+1/2,z', &         
  'x,-y+1/2,z+1/2' /                                                                                                                

 data spacegroupdata(3937:3961)/ &                              ! group index = 503, number of symmetry op = 24
  '*206        Th^7          Ia-3                        -I 2b 2c 3', &             
  'x,y,z','z,x,y','y,z,x','-y,-z+1/2,x','z,-x,-y+1/2','-y+1/2,z,-x','-z,-x+1/2,y','-z+1/2,x,-y','y,-z,-x+1/2','-x,-y+1/2,z', &      
  'x,-y,-z+1/2','-x+1/2,y,-z','-x,-y,-z','-z,-x,-y','-y,-z,-x','y,z+1/2,-x','-z,x,y+1/2','y+1/2,-z,x','z,x+1/2,-y','z+1/2,-x,y', &  
  '-y,z,x+1/2','x,y+1/2,-z','-x,y,z+1/2','x+1/2,-y,z' /                                                                             

 data spacegroupdata(3962:3986)/ &                              ! group index = 504, number of symmetry op = 24
  '*207        O^1           P432                         P 4 2 3', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-z,y','x,-y,-z','x,z,-y','z,y,-x','-x,y,-z','-z,y,x','z,x,y','y,z,x','-y,-z,x','z,-x,-y', &
  '-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','y,x,-z','-y,-x,-z','-x,z,y','-x,-z,-y','z,-y,x','-z,-y,-x' /                             

 data spacegroupdata(3987:4011)/ &                              ! group index = 505, number of symmetry op = 24
  '*208        O^2           P4232                        P 4n 2 3', &              
  'x,y,z','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','x+1/2,-z+1/2,y+1/2','x,-y,-z','x+1/2,z+1/2,-y+1/2', &                
  'z+1/2,y+1/2,-x+1/2','-x,y,-z','-z+1/2,y+1/2,x+1/2','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x', &
  'y+1/2,x+1/2,-z+1/2','-y+1/2,-x+1/2,-z+1/2','-x+1/2,z+1/2,y+1/2','-x+1/2,-z+1/2,-y+1/2','z+1/2,-y+1/2,x+1/2', &                   
  '-z+1/2,-y+1/2,-x+1/2' /                                                                                                          

 data spacegroupdata(4012:4036)/ &                              ! group index = 506, number of symmetry op = 24
  '*209        O^3           F432                         F 4 2 3', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-z,y','x,-y,-z','x,z,-y','z,y,-x','-x,y,-z','-z,y,x','z,x,y','y,z,x','-y,-z,x','z,-x,-y', &
  '-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','y,x,-z','-y,-x,-z','-x,z,y','-x,-z,-y','z,-y,x','-z,-y,-x' /                             

 data spacegroupdata(4037:4061)/ &                              ! group index = 507, number of symmetry op = 24
  '*210        O^4           F4132                        F 4d 2 3', &              
  'x,y,z','-y+1/4,x+1/4,z+1/4','-x,-y,z','y+1/4,-x+1/4,z+1/4','x+1/4,-z+1/4,y+1/4','x,-y,-z','x+1/4,z+1/4,-y+1/4', &                
  'z+1/4,y+1/4,-x+1/4','-x,y,-z','-z+1/4,y+1/4,x+1/4','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x', &
  'y+1/4,x+1/4,-z+1/4','-y+1/4,-x+1/4,-z+1/4','-x+1/4,z+1/4,y+1/4','-x+1/4,-z+1/4,-y+1/4','z+1/4,-y+1/4,x+1/4', &                   
  '-z+1/4,-y+1/4,-x+1/4' /                                                                                                          

 data spacegroupdata(4062:4086)/ &                              ! group index = 508, number of symmetry op = 24
  '*211        O^5           I432                         I 4 2 3', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-z,y','x,-y,-z','x,z,-y','z,y,-x','-x,y,-z','-z,y,x','z,x,y','y,z,x','-y,-z,x','z,-x,-y', &
  '-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','y,x,-z','-y,-x,-z','-x,z,y','-x,-z,-y','z,-y,x','-z,-y,-x' /                             

 data spacegroupdata(4087:4111)/ &                              ! group index = 509, number of symmetry op = 24
  '*212        O^6           P4332                        P 4acd 2ab 3', &          
  'x,y,z','-y+3/4,x+1/4,z+3/4','-x+1/2,-y,z+1/2','y+3/4,-x+3/4,z+1/4','x+3/4,-z+3/4,y+1/4','x+1/2,-y+1/2,-z','x+1/4,z+3/4,-y+3/4', &
  'z+1/4,y+3/4,-x+3/4','-x,y+1/2,-z+1/2','-z+3/4,y+1/4,x+3/4','z,x,y','y,z,x','-y+1/2,-z,x+1/2','z+1/2,-x+1/2,-y', &                
  '-y,z+1/2,-x+1/2','-z+1/2,-x,y+1/2','-z,x+1/2,-y+1/2','y+1/2,-z+1/2,-x','y+1/4,x+3/4,-z+3/4','-y+1/4,-x+1/4,-z+1/4', &            
  '-x+3/4,z+1/4,y+3/4','-x+1/4,-z+1/4,-y+1/4','z+3/4,-y+3/4,x+1/4','-z+1/4,-y+1/4,-x+1/4' /                                         

 data spacegroupdata(4112:4136)/ &                              ! group index = 510, number of symmetry op = 24
  '*213        O^7           P4132                        P 4bd 2ab 3', &           
  'x,y,z','-y+1/4,x+3/4,z+1/4','-x+1/2,-y,z+1/2','y+1/4,-x+1/4,z+3/4','x+1/4,-z+1/4,y+3/4','x+1/2,-y+1/2,-z','x+3/4,z+1/4,-y+1/4', &
  'z+3/4,y+1/4,-x+1/4','-x,y+1/2,-z+1/2','-z+1/4,y+3/4,x+1/4','z,x,y','y,z,x','-y+1/2,-z,x+1/2','z+1/2,-x+1/2,-y', &                
  '-y,z+1/2,-x+1/2','-z+1/2,-x,y+1/2','-z,x+1/2,-y+1/2','y+1/2,-z+1/2,-x','y+3/4,x+1/4,-z+1/4','-y+3/4,-x+3/4,-z+3/4', &            
  '-x+1/4,z+3/4,y+1/4','-x+3/4,-z+3/4,-y+3/4','z+1/4,-y+1/4,x+3/4','-z+3/4,-y+3/4,-x+3/4' /                                         

 data spacegroupdata(4137:4161)/ &                              ! group index = 511, number of symmetry op = 24
  '*214        O^8           I4132                        I 4bd 2c 3', &            
  'x,y,z','-y+1/4,x+3/4,z+1/4','-x,-y+1/2,z','y+1/4,-x+1/4,z+3/4','x+1/4,-z+1/4,y+3/4','x,-y,-z+1/2','x+3/4,z+1/4,-y+1/4', &        
  'z+3/4,y+1/4,-x+1/4','-x+1/2,y,-z','-z+1/4,y+3/4,x+1/4','z,x,y','y,z,x','-y,-z+1/2,x','z,-x,-y+1/2','-y+1/2,z,-x','-z,-x+1/2,y', &
  '-z+1/2,x,-y','y,-z,-x+1/2','y+3/4,x+1/4,-z+1/4','-y+1/4,-x+1/4,-z+1/4','-x+1/4,z+3/4,y+1/4','-x+1/4,-z+1/4,-y+1/4', &            
  'z+1/4,-y+1/4,x+3/4','-z+1/4,-y+1/4,-x+1/4' /                                                                                     

 data spacegroupdata(4162:4186)/ &                              ! group index = 512, number of symmetry op = 24
  '*215        Td^1          P-43m                        P -4 2 3', &              
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','-x,z,-y','x,-y,-z','-x,-z,y','-z,-y,x','-x,y,-z','z,-y,-x','z,x,y','y,z,x','-y,-z,x', &    
  'z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-y,-x,z','y,x,z','x,-z,-y','x,z,y','-z,y,-x','z,y,x' /                         

 data spacegroupdata(4187:4211)/ &                              ! group index = 513, number of symmetry op = 24
  '*216        Td^2          F-43m                        F -4 2 3', &              
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','-x,z,-y','x,-y,-z','-x,-z,y','-z,-y,x','-x,y,-z','z,-y,-x','z,x,y','y,z,x','-y,-z,x', &    
  'z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-y,-x,z','y,x,z','x,-z,-y','x,z,y','-z,y,-x','z,y,x' /                         

 data spacegroupdata(4212:4236)/ &                              ! group index = 514, number of symmetry op = 24
  '*217        Td^3          I-43m                        I -4 2 3', &              
  'x,y,z','y,-x,-z','-x,-y,z','-y,x,-z','-x,z,-y','x,-y,-z','-x,-z,y','-z,-y,x','-x,y,-z','z,-y,-x','z,x,y','y,z,x','-y,-z,x', &    
  'z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-y,-x,z','y,x,z','x,-z,-y','x,z,y','-z,y,-x','z,y,x' /                         

 data spacegroupdata(4237:4261)/ &                              ! group index = 515, number of symmetry op = 24
  '*218        Td^4          P-43n                        P -4n 2 3', &             
  'x,y,z','y+1/2,-x+1/2,-z+1/2','-x,-y,z','-y+1/2,x+1/2,-z+1/2','-x+1/2,z+1/2,-y+1/2','x,-y,-z','-x+1/2,-z+1/2,y+1/2', &            
  '-z+1/2,-y+1/2,x+1/2','-x,y,-z','z+1/2,-y+1/2,-x+1/2','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y', &        
  'y,-z,-x','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2','x+1/2,-z+1/2,-y+1/2','x+1/2,z+1/2,y+1/2','-z+1/2,y+1/2,-x+1/2', &            
  'z+1/2,y+1/2,x+1/2' /                                                                                                             

 data spacegroupdata(4262:4286)/ &                              ! group index = 516, number of symmetry op = 24
  '*219        Td^5          F-43c                        F -4c 2 3', &             
  'x,y,z','y,-x,-z+1/2','-x,-y,z','-y,x,-z+1/2','-x,z,-y+1/2','x,-y,-z','-x,-z,y+1/2','-z,-y,x+1/2','-x,y,-z','z,-y,-x+1/2', &      
  'z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-y,-x,z+1/2','y,x,z+1/2','x,-z,-y+1/2','x,z,y+1/2', &
  '-z,y,-x+1/2','z,y,x+1/2' /                                                                                                       

 data spacegroupdata(4287:4311)/ &                              ! group index = 517, number of symmetry op = 24
  '*220        Td^6          I-43d                        I -4bd 2c 3', &           
  'x,y,z','y+1/4,-x+3/4,-z+1/4','-x,-y+1/2,z','-y+1/4,x+1/4,-z+3/4','-x+1/4,z+1/4,-y+3/4','x,-y,-z+1/2','-x+3/4,-z+1/4,y+1/4', &    
  '-z+3/4,-y+1/4,x+1/4','-x+1/2,y,-z','z+1/4,-y+3/4,-x+1/4','z,x,y','y,z,x','-y,-z+1/2,x','z,-x,-y+1/2','-y+1/2,z,-x', &            
  '-z,-x+1/2,y','-z+1/2,x,-y','y,-z,-x+1/2','-y+3/4,-x+1/4,z+1/4','y+1/4,x+1/4,z+1/4','x+1/4,-z+3/4,-y+1/4','x+1/4,z+1/4,y+1/4', &  
  '-z+1/4,y+1/4,-x+3/4','z+1/4,y+1/4,x+1/4' /                                                                                       

 data spacegroupdata(4312:4360)/ &                              ! group index = 518, number of symmetry op = 48
  '*221        Oh^1          Pm-3m                       -P 4 2 3', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-z,y','x,-y,-z','x,z,-y','z,y,-x','-x,y,-z','-z,y,x','z,x,y','y,z,x','-y,-z,x','z,-x,-y', &
  '-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','y,x,-z','-y,-x,-z','-x,z,y','-x,-z,-y','z,-y,x','-z,-y,-x','-x,-y,-z','y,-x,-z', &       
  'x,y,-z','-y,x,-z','-x,z,-y','-x,y,z','-x,-z,y','-z,-y,x','x,-y,z','z,-y,-x','-z,-x,-y','-y,-z,-x','y,z,-x','-z,x,y','y,-z,x', &  
  'z,x,-y','z,-x,y','-y,z,x','-y,-x,z','y,x,z','x,-z,-y','x,z,y','-z,y,-x','z,y,x' /                                                

 data spacegroupdata(4361:4409)/ &                              ! group index = 519, number of symmetry op = 48
  '*222:1      Oh^2          Pn-3n:1                      P 4 2 3 -1n', &           
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','-y,x,z','-x,-y,z','y,-x,z','x,-z,y','x,-y,-z','x,z,-y','z,y,-x','-x,y,-z','-z,y,x', &             
  'y+1/2,-x+1/2,-z+1/2','-y+1/2,x+1/2,-z+1/2','-x+1/2,z+1/2,-y+1/2','-x+1/2,-z+1/2,y+1/2','-z+1/2,-y+1/2,x+1/2', &                  
  'z+1/2,-y+1/2,-x+1/2','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-z+1/2,-x+1/2,-y+1/2', &       
  '-y+1/2,-z+1/2,-x+1/2','y+1/2,z+1/2,-x+1/2','-z+1/2,x+1/2,y+1/2','y+1/2,-z+1/2,x+1/2','z+1/2,x+1/2,-y+1/2','z+1/2,-x+1/2,y+1/2', &
  '-y+1/2,z+1/2,x+1/2','y,x,-z','-y,-x,-z','-x,z,y','-x,-z,-y','z,-y,x','-z,-y,-x','x+1/2,y+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2', &     
  'x+1/2,-y+1/2,z+1/2','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2','x+1/2,-z+1/2,-y+1/2','x+1/2,z+1/2,y+1/2','-z+1/2,y+1/2,-x+1/2', & 
  'z+1/2,y+1/2,x+1/2' /                                                                                                             

 data spacegroupdata(4410:4458)/ &                              ! group index = 520, number of symmetry op = 48
  '*222:2      Oh^2          Pn-3n:2                     -P 4a 2bc 3', &            
  'x,y,z','-y+1/2,x,z','-x+1/2,-y+1/2,z','y,-x+1/2,z','x,-z+1/2,y','x,-y+1/2,-z+1/2','x,z,-y+1/2','z,y,-x+1/2','-x+1/2,y,-z+1/2', & 
  '-z+1/2,y,x','z,x,y','y,z,x','-y+1/2,-z+1/2,x','z,-x+1/2,-y+1/2','-y+1/2,z,-x+1/2','-z+1/2,-x+1/2,y','-z+1/2,x,-y+1/2', &         
  'y,-z+1/2,-x+1/2','y,x,-z+1/2','-y+1/2,-x+1/2,-z+1/2','-x+1/2,z,y','-x+1/2,-z+1/2,-y+1/2','z,-y+1/2,x','-z+1/2,-y+1/2,-x+1/2', &  
  '-x,-y,-z','y+1/2,-x,-z','x+1/2,y+1/2,-z','-y,x+1/2,-z','-x,z+1/2,-y','-x,y+1/2,z+1/2','-x,-z,y+1/2','-z,-y,x+1/2', &             
  'x+1/2,-y,z+1/2','z+1/2,-y,-x','-z,-x,-y','-y,-z,-x','y+1/2,z+1/2,-x','-z,x+1/2,y+1/2','y+1/2,-z,x+1/2','z+1/2,x+1/2,-y', &       
  'z+1/2,-x,y+1/2','-y,z+1/2,x+1/2','-y,-x,z+1/2','y+1/2,x+1/2,z+1/2','x+1/2,-z,-y','x+1/2,z+1/2,y+1/2','-z,y+1/2,-x', &            
  'z+1/2,y+1/2,x+1/2' /                                                                                                             

 data spacegroupdata(4459:4507)/ &                              ! group index = 521, number of symmetry op = 48
  '*223        Oh^3          Pm-3n                       -P 4n 2 3', &              
  'x,y,z','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','x+1/2,-z+1/2,y+1/2','x,-y,-z','x+1/2,z+1/2,-y+1/2', &                
  'z+1/2,y+1/2,-x+1/2','-x,y,-z','-z+1/2,y+1/2,x+1/2','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x', &
  'y+1/2,x+1/2,-z+1/2','-y+1/2,-x+1/2,-z+1/2','-x+1/2,z+1/2,y+1/2','-x+1/2,-z+1/2,-y+1/2','z+1/2,-y+1/2,x+1/2', &                   
  '-z+1/2,-y+1/2,-x+1/2','-x,-y,-z','y+1/2,-x+1/2,-z+1/2','x,y,-z','-y+1/2,x+1/2,-z+1/2','-x+1/2,z+1/2,-y+1/2','-x,y,z', &          
  '-x+1/2,-z+1/2,y+1/2','-z+1/2,-y+1/2,x+1/2','x,-y,z','z+1/2,-y+1/2,-x+1/2','-z,-x,-y','-y,-z,-x','y,z,-x','-z,x,y','y,-z,x', &    
  'z,x,-y','z,-x,y','-y,z,x','-y+1/2,-x+1/2,z+1/2','y+1/2,x+1/2,z+1/2','x+1/2,-z+1/2,-y+1/2','x+1/2,z+1/2,y+1/2', &                 
  '-z+1/2,y+1/2,-x+1/2','z+1/2,y+1/2,x+1/2' /                                                                                       

 data spacegroupdata(4508:4556)/ &                              ! group index = 522, number of symmetry op = 48
  '*224:1      Oh^4          Pn-3m:1                      P 4n 2 3 -1n', &          
  'x,y,z','-x+1/2,-y+1/2,-z+1/2','-y+1/2,x+1/2,z+1/2','-x,-y,z','y+1/2,-x+1/2,z+1/2','x+1/2,-z+1/2,y+1/2','x,-y,-z', &              
  'x+1/2,z+1/2,-y+1/2','z+1/2,y+1/2,-x+1/2','-x,y,-z','-z+1/2,y+1/2,x+1/2','y,-x,-z','-y,x,-z','-x,z,-y','-x,-z,y','-z,-y,x', &     
  'z,-y,-x','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-z+1/2,-x+1/2,-y+1/2', &                   
  '-y+1/2,-z+1/2,-x+1/2','y+1/2,z+1/2,-x+1/2','-z+1/2,x+1/2,y+1/2','y+1/2,-z+1/2,x+1/2','z+1/2,x+1/2,-y+1/2','z+1/2,-x+1/2,y+1/2', &
  '-y+1/2,z+1/2,x+1/2','y+1/2,x+1/2,-z+1/2','-y+1/2,-x+1/2,-z+1/2','-x+1/2,z+1/2,y+1/2','-x+1/2,-z+1/2,-y+1/2', &                   
  'z+1/2,-y+1/2,x+1/2','-z+1/2,-y+1/2,-x+1/2','x+1/2,y+1/2,-z+1/2','-x+1/2,y+1/2,z+1/2','x+1/2,-y+1/2,z+1/2','-y,-x,z','y,x,z', &   
  'x,-z,-y','x,z,y','-z,y,-x','z,y,x' /                                                                                             

 data spacegroupdata(4557:4605)/ &                              ! group index = 523, number of symmetry op = 48
  '*224:2      Oh^4          Pn-3m:2                     -P 4bc 2bc 3', &           
  'x,y,z','-y,x+1/2,z+1/2','-x+1/2,-y+1/2,z','y+1/2,-x,z+1/2','x+1/2,-z,y+1/2','x,-y+1/2,-z+1/2','x+1/2,z+1/2,-y', &                
  'z+1/2,y+1/2,-x','-x+1/2,y,-z+1/2','-z,y+1/2,x+1/2','z,x,y','y,z,x','-y+1/2,-z+1/2,x','z,-x+1/2,-y+1/2','-y+1/2,z,-x+1/2', &      
  '-z+1/2,-x+1/2,y','-z+1/2,x,-y+1/2','y,-z+1/2,-x+1/2','y+1/2,x+1/2,-z','-y,-x,-z','-x,z+1/2,y+1/2','-x,-z,-y','z+1/2,-y,x+1/2', & 
  '-z,-y,-x','-x,-y,-z','y,-x+1/2,-z+1/2','x+1/2,y+1/2,-z','-y+1/2,x,-z+1/2','-x+1/2,z,-y+1/2','-x,y+1/2,z+1/2','-x+1/2,-z+1/2,y', &
  '-z+1/2,-y+1/2,x','x+1/2,-y,z+1/2','z,-y+1/2,-x+1/2','-z,-x,-y','-y,-z,-x','y+1/2,z+1/2,-x','-z,x+1/2,y+1/2','y+1/2,-z,x+1/2', &  
  'z+1/2,x+1/2,-y','z+1/2,-x,y+1/2','-y,z+1/2,x+1/2','-y+1/2,-x+1/2,z','y,x,z','x,-z+1/2,-y+1/2','x,z,y','-z+1/2,y,-x+1/2', &       
  'z,y,x' /                                                                                                                         

 data spacegroupdata(4606:4654)/ &                              ! group index = 524, number of symmetry op = 48
  '*225        Oh^5          Fm-3m                       -F 4 2 3', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-z,y','x,-y,-z','x,z,-y','z,y,-x','-x,y,-z','-z,y,x','z,x,y','y,z,x','-y,-z,x','z,-x,-y', &
  '-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','y,x,-z','-y,-x,-z','-x,z,y','-x,-z,-y','z,-y,x','-z,-y,-x','-x,-y,-z','y,-x,-z', &       
  'x,y,-z','-y,x,-z','-x,z,-y','-x,y,z','-x,-z,y','-z,-y,x','x,-y,z','z,-y,-x','-z,-x,-y','-y,-z,-x','y,z,-x','-z,x,y','y,-z,x', &  
  'z,x,-y','z,-x,y','-y,z,x','-y,-x,z','y,x,z','x,-z,-y','x,z,y','-z,y,-x','z,y,x' /                                                

 data spacegroupdata(4655:4703)/ &                              ! group index = 525, number of symmetry op = 48
  '*226        Oh^6          Fm-3c                       -F 4c 2 3', &              
  'x,y,z','-y,x,z+1/2','-x,-y,z','y,-x,z+1/2','x,-z,y+1/2','x,-y,-z','x,z,-y+1/2','z,y,-x+1/2','-x,y,-z','-z,y,x+1/2','z,x,y', &    
  'y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','y,x,-z+1/2','-y,-x,-z+1/2','-x,z,y+1/2','-x,-z,-y+1/2', &    
  'z,-y,x+1/2','-z,-y,-x+1/2','-x,-y,-z','y,-x,-z+1/2','x,y,-z','-y,x,-z+1/2','-x,z,-y+1/2','-x,y,z','-x,-z,y+1/2','-z,-y,x+1/2', & 
  'x,-y,z','z,-y,-x+1/2','-z,-x,-y','-y,-z,-x','y,z,-x','-z,x,y','y,-z,x','z,x,-y','z,-x,y','-y,z,x','-y,-x,z+1/2','y,x,z+1/2', &   
  'x,-z,-y+1/2','x,z,y+1/2','-z,y,-x+1/2','z,y,x+1/2' /                                                                             

 data spacegroupdata(4704:4752)/ &                              ! group index = 526, number of symmetry op = 48
  '*227:1      Oh^7          Fd-3m:1                      F 4d 2 3 -1d', &          
  'x,y,z','-x+1/4,-y+1/4,-z+1/4','-y+1/4,x+1/4,z+1/4','-x,-y,z','y+1/4,-x+1/4,z+1/4','x+1/4,-z+1/4,y+1/4','x,-y,-z', &              
  'x+1/4,z+1/4,-y+1/4','z+1/4,y+1/4,-x+1/4','-x,y,-z','-z+1/4,y+1/4,x+1/4','y,-x,-z','-y,x,-z','-x,z,-y','-x,-z,y','-z,-y,x', &     
  'z,-y,-x','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','-z+1/4,-x+1/4,-y+1/4', &                   
  '-y+1/4,-z+1/4,-x+1/4','y+1/4,z+1/4,-x+1/4','-z+1/4,x+1/4,y+1/4','y+1/4,-z+1/4,x+1/4','z+1/4,x+1/4,-y+1/4','z+1/4,-x+1/4,y+1/4', &
  '-y+1/4,z+1/4,x+1/4','y+1/4,x+1/4,-z+1/4','-y+1/4,-x+1/4,-z+1/4','-x+1/4,z+1/4,y+1/4','-x+1/4,-z+1/4,-y+1/4', &                   
  'z+1/4,-y+1/4,x+1/4','-z+1/4,-y+1/4,-x+1/4','x+1/4,y+1/4,-z+1/4','-x+1/4,y+1/4,z+1/4','x+1/4,-y+1/4,z+1/4','-y,-x,z','y,x,z', &   
  'x,-z,-y','x,z,y','-z,y,-x','z,y,x' /                                                                                             

 data spacegroupdata(4753:4801)/ &                              ! group index = 527, number of symmetry op = 48
  '*227:2      Oh^7          Fd-3m:2                     -F 4vw 2vw 3', &           
  'x,y,z','-y,x+1/4,z+1/4','-x+1/4,-y+1/4,z','y+1/4,-x,z+1/4','x+1/4,-z,y+1/4','x,-y+1/4,-z+1/4','x+1/4,z+1/4,-y', &                
  'z+1/4,y+1/4,-x','-x+1/4,y,-z+1/4','-z,y+1/4,x+1/4','z,x,y','y,z,x','-y+1/4,-z+1/4,x','z,-x+1/4,-y+1/4','-y+1/4,z,-x+1/4', &      
  '-z+1/4,-x+1/4,y','-z+1/4,x,-y+1/4','y,-z+1/4,-x+1/4','y+1/4,x+1/4,-z','-y,-x,-z','-x,z+1/4,y+1/4','-x,-z,-y','z+1/4,-y,x+1/4', & 
  '-z,-y,-x','-x,-y,-z','y,-x+3/4,-z+3/4','x+3/4,y+3/4,-z','-y+3/4,x,-z+3/4','-x+3/4,z,-y+3/4','-x,y+3/4,z+3/4','-x+3/4,-z+3/4,y', &
  '-z+3/4,-y+3/4,x','x+3/4,-y,z+3/4','z,-y+3/4,-x+3/4','-z,-x,-y','-y,-z,-x','y+3/4,z+3/4,-x','-z,x+3/4,y+3/4','y+3/4,-z,x+3/4', &  
  'z+3/4,x+3/4,-y','z+3/4,-x,y+3/4','-y,z+3/4,x+3/4','-y+3/4,-x+3/4,z','y,x,z','x,-z+3/4,-y+3/4','x,z,y','-z+3/4,y,-x+3/4', &       
  'z,y,x' /                                                                                                                         

 data spacegroupdata(4802:4850)/ &                              ! group index = 528, number of symmetry op = 48
  '*227:S      Oh^7          Fd-3mS                      -F d - 3 m S', &           
  'z+1/4,y+1/4,-x+1/4','y+1/4,x+1/4,-z+1/4','x+1/4,z+1/4,-y+1/4','z+1/4,x+1/4,-y+1/4','y+1/4,z+1/4,-x+1/4','x+1/4,y+1/4,-z+1/4', &  
  'z+1/4,-y+1/4,x+1/4','y+1/4,-x+1/4,z+1/4','x+1/4,-z+1/4,y+1/4','z+1/4,-x+1/4,y+1/4','y+1/4,-z+1/4,x+1/4','x+1/4,-y+1/4,z+1/4', &  
  '-z+1/4,y+1/4,x+1/4','-y+1/4,x+1/4,z+1/4','-x+1/4,z+1/4,y+1/4','-z+1/4,x+1/4,y+1/4','-y+1/4,z+1/4,x+1/4','-x+1/4,y+1/4,z+1/4', &  
  '-z+1/4,-y+1/4,-x+1/4','-y+1/4,-x+1/4,-z+1/4','-x+1/4,-z+1/4,-y+1/4','-z+1/4,-x+1/4,-y+1/4','-y+1/4,-z+1/4,-x+1/4', &             
  '-x+1/4,-y+1/4,-z+1/4','-z,-y,x','-y,-x,z','-x,-z,y','-z,-x,y','-y,-z,x','-x,-y,z','-z,y,-x','-y,x,-z','-x,z,-y','-z,x,-y', &     
  '-y,z,-x','-x,y,-z','z,-y,-x','y,-x,-z','x,-z,-y','z,-x,-y','y,-z,-x','x,-y,-z','z,y,x','y,x,z','x,z,y','z,x,y','y,z,x','x,y,z' / 

 data spacegroupdata(4851:4899)/ &                              ! group index = 529, number of symmetry op = 48
  '*228:1      Oh^8          Fd-3c:1                      F 4d 2 3 -1cd', &         
  'x,y,z','-x+1/4,-y+1/4,-z+3/4','-y+1/4,x+1/4,z+1/4','-x,-y,z','y+1/4,-x+1/4,z+1/4','x+1/4,-z+1/4,y+1/4','x,-y,-z', &              
  'x+1/4,z+1/4,-y+1/4','z+1/4,y+1/4,-x+1/4','-x,y,-z','-z+1/4,y+1/4,x+1/4','y,-x,-z+1/2','-y,x,-z+1/2','-x,z,-y+1/2', &             
  '-x,-z,y+1/2','-z,-y,x+1/2','z,-y,-x+1/2','z,x,y','y,z,x','-y,-z,x','z,-x,-y','-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x', &          
  '-z+1/4,-x+1/4,-y+3/4','-y+1/4,-z+1/4,-x+3/4','y+1/4,z+1/4,-x+3/4','-z+1/4,x+1/4,y+3/4','y+1/4,-z+1/4,x+3/4', &                   
  'z+1/4,x+1/4,-y+3/4','z+1/4,-x+1/4,y+3/4','-y+1/4,z+1/4,x+3/4','y+1/4,x+1/4,-z+1/4','-y+1/4,-x+1/4,-z+1/4','-x+1/4,z+1/4,y+1/4', &
  '-x+1/4,-z+1/4,-y+1/4','z+1/4,-y+1/4,x+1/4','-z+1/4,-y+1/4,-x+1/4','x+1/4,y+1/4,-z+3/4','-x+1/4,y+1/4,z+3/4', &                   
  'x+1/4,-y+1/4,z+3/4','-y,-x,z+1/2','y,x,z+1/2','x,-z,-y+1/2','x,z,y+1/2','-z,y,-x+1/2','z,y,x+1/2' /                              

 data spacegroupdata(4900:4948)/ &                              ! group index = 530, number of symmetry op = 48
  '*228:2      Oh^8          Fd-3c:2                     -F 4cvw 2vw 3', &          
  'x,y,z','-y,x+1/4,z+3/4','-x+1/4,-y+1/4,z','y+1/4,-x,z+3/4','x+1/4,-z,y+3/4','x,-y+1/4,-z+1/4','x+1/4,z+3/4,-y', &                
  'z+1/4,y+3/4,-x','-x+1/4,y,-z+1/4','-z,y+1/4,x+3/4','z,x,y','y,z,x','-y+1/4,-z+1/4,x','z,-x+1/4,-y+1/4','-y+1/4,z,-x+1/4', &      
  '-z+1/4,-x+1/4,y','-z+1/4,x,-y+1/4','y,-z+1/4,-x+1/4','y+1/4,x+3/4,-z','-y,-x,-z+1/2','-x,z+1/4,y+3/4','-x,-z,-y+1/2', &          
  'z+1/4,-y,x+3/4','-z,-y,-x+1/2','-x,-y,-z','y,-x+3/4,-z+1/4','x+3/4,y+3/4,-z','-y+3/4,x,-z+1/4','-x+3/4,z,-y+1/4', &              
  '-x,y+3/4,z+3/4','-x+3/4,-z+1/4,y','-z+3/4,-y+1/4,x','x+3/4,-y,z+3/4','z,-y+3/4,-x+1/4','-z,-x,-y','-y,-z,-x','y+3/4,z+3/4,-x', & 
  '-z,x+3/4,y+3/4','y+3/4,-z,x+3/4','z+3/4,x+3/4,-y','z+3/4,-x,y+3/4','-y,z+3/4,x+3/4','-y+3/4,-x+1/4,z','y,x,z+1/2', &             
  'x,-z+3/4,-y+1/4','x,z,y+1/2','-z+3/4,y,-x+1/4','z,y,x+1/2' /                                                                     

 data spacegroupdata(4949:4997)/ &                              ! group index = 531, number of symmetry op = 48
  '*229        Oh^9          Im-3m                       -I 4 2 3', &               
  'x,y,z','-y,x,z','-x,-y,z','y,-x,z','x,-z,y','x,-y,-z','x,z,-y','z,y,-x','-x,y,-z','-z,y,x','z,x,y','y,z,x','-y,-z,x','z,-x,-y', &
  '-y,z,-x','-z,-x,y','-z,x,-y','y,-z,-x','y,x,-z','-y,-x,-z','-x,z,y','-x,-z,-y','z,-y,x','-z,-y,-x','-x,-y,-z','y,-x,-z', &       
  'x,y,-z','-y,x,-z','-x,z,-y','-x,y,z','-x,-z,y','-z,-y,x','x,-y,z','z,-y,-x','-z,-x,-y','-y,-z,-x','y,z,-x','-z,x,y','y,-z,x', &  
  'z,x,-y','z,-x,y','-y,z,x','-y,-x,z','y,x,z','x,-z,-y','x,z,y','-z,y,-x','z,y,x' /                                                

 data spacegroupdata(4998:5046)/ &                              ! group index = 532, number of symmetry op = 48
  '*230        Oh^10         Ia-3d                       -I 4bd 2c 3', &            
  'x,y,z','-y+1/4,x+3/4,z+1/4','-x,-y+1/2,z','y+1/4,-x+1/4,z+3/4','x+1/4,-z+1/4,y+3/4','x,-y,-z+1/2','x+3/4,z+1/4,-y+1/4', &        
  'z+3/4,y+1/4,-x+1/4','-x+1/2,y,-z','-z+1/4,y+3/4,x+1/4','z,x,y','y,z,x','-y,-z+1/2,x','z,-x,-y+1/2','-y+1/2,z,-x','-z,-x+1/2,y', &
  '-z+1/2,x,-y','y,-z,-x+1/2','y+3/4,x+1/4,-z+1/4','-y+1/4,-x+1/4,-z+1/4','-x+1/4,z+3/4,y+1/4','-x+1/4,-z+1/4,-y+1/4', &            
  'z+1/4,-y+1/4,x+3/4','-z+1/4,-y+1/4,-x+1/4','-x,-y,-z','y+3/4,-x+1/4,-z+3/4','x,y+1/2,-z','-y+3/4,x+3/4,-z+1/4', &                
  '-x+3/4,z+3/4,-y+1/4','-x,y,z+1/2','-x+1/4,-z+3/4,y+3/4','-z+1/4,-y+3/4,x+3/4','x+1/2,-y,z','z+3/4,-y+1/4,-x+3/4','-z,-x,-y', &   
  '-y,-z,-x','y,z+1/2,-x','-z,x,y+1/2','y+1/2,-z,x','z,x+1/2,-y','z+1/2,-x,y','-y,z,x+1/2','-y+1/4,-x+3/4,z+3/4', &                 
  'y+3/4,x+3/4,z+3/4','x+3/4,-z+1/4,-y+3/4','x+3/4,z+3/4,y+3/4','-z+3/4,y+3/4,-x+1/4','z+3/4,y+3/4,x+3/4'/                          

end