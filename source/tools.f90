! Subroutines for FDMNES
! Contains different tools

!***********************************************************************

subroutine ad_number(ib,nomfich,Length)

  character(len=Length) nomfich

  l = len_trim(nomfich)

  if( ib < 0 ) then
    l = l + 1
    if( l <= Length ) nomfich(l:l) = '-'
  endif

  i = abs(ib)

  do iu = 2,10
    if( i / 10**(iu-1) < 1 ) exit
  end do
  iumax = iu - 1

  do iu = iumax,1,-1

    ipuis = 10**(iu-1)

    in = i / ipuis

! S'il n'y a pas la place on ecrit que les derniers chiffres
    if( l + iu <= Length ) then
      l = l + 1
      nomfich(l:l) = achar(in+48)
    elseif( l + iu == Length-1 ) then
      if( nomfich(l:l) == '_' ) nomfich(l:l) = achar(in+48)
    endif

    i = i - ipuis * in

  end do

  return
end

!***********************************************************************

subroutine ad_number_r(r,nom)

  use declarations
  implicit none

  integer:: i, l, l2, l3

  character(len=1) a
  character(len=Length_word) mot, nom

  real(kind=db):: r
  
  if( Length_word >= 15 ) then 
    write(mot,'(f15.4)') r
  elseif( Length_word >= 13 ) then 
    write(mot,'(f13.3)') r
  elseif( Length_word >= 10 ) then 
    write(mot,'(f10.2)') r
  else
    write(mot,'(f5.2)') r
  endif
  
  mot = adjustl( mot )
  
  l2 = len_trim(mot) 
  l3 = l2

  do i = l3,2,-1
    a = mot(i:i)
    if( a /= '0' .and. a /= '.' ) exit
    mot(i:i) = ' '
    l2 = l2 - 1
    if( a == '.' ) exit
  end do
 
  l = len_trim( nom )

  do i = l+1,min(l+l2,Length_word) 
    nom(i:i) = mot(i-l:i-l)
  end do
      
  return
end

!***********************************************************************

subroutine ad_sufix(Core_resolved,initr,jseuil,Length,mot,nseuil)

  use declarations
  implicit none

  integer:: initr, jseuil, l, Length, nseuil
    
  character(len=Length):: mot
  
  Logical:: Core_resolved 

  l = len_trim(mot)

  if( Core_resolved ) then
    if( l < Length - 2 ) mot(l+1:l+1) = '_'
    call ad_number(initr,mot,Length)
  else  
    if( l < Length - 3 ) then
      l = l + 1
      mot(l:l) = '_'
    endif
    if( l < Length - 2 ) then
      l = l + 1
      mot(l:l) = achar(nseuil+74)
    endif
    call ad_number(jseuil+initr-1,mot,Length)
  endif

  return

end

!***********************************************************************

subroutine center_word( mot, Length_word )

  character(len=*) mot
  character(len=Length_word) mot2

  mot2 = ' '
  mot = adjustl( mot )
  l = len_trim( mot )
  lshift = ( Length_word - l + 1 ) / 2
  lshift = max( 0, lshift )
  lm = min( l, Length_word )
  mot2(1+lshift:lm+lshift) = mot(1:lm)
  mot = mot2

  return
end

!***********************************************************************

! Fonction giving the number of number in the next non empty line et setting at the biginning of this line.
! With line starting with a character, nnombre = 0

function nnombre(irec,Length)

  use declarations
  implicit none

  integer, parameter:: nlet = 53

  integer:: eof, i, icol, irec, j, k, l, Length, ligne, nmots, n, nnombre

  character(len=Length):: mot, test
  character(len=1), dimension(nlet):: let

  real(kind=db):: x

  data let /'d','e','D','E','a','A','b','B','c','C','f','F','g','G','h','H','i','I','j','J', &
            'k','K','l','L','m','M','n','N', 'o','O','p','P','q','Q','r','R','s','T','t','T','u','U','v','V', &
            'w','W','x','X','y','Y','z','Z','/'/

  do i = 1,Length
    mot(i:i) = 'P'
  end do

  nmots = 0

  boucle_ligne: do ligne = 1,1000000

    read(irec,'(A)',iostat=eof) mot

! An "END=" condition leaves the file position after the endfile record,
! thus backspace to be just before the endfile record. This way a
! subsequent READ with another "END=" clause can execute correctly.
    if( eof /= 0 ) then
      nnombre = 0
      backspace(irec)
      return
    endif

    do i = 1,Length
      if( mot(i:i) /= ' ' .and. mot(i:i) /= char(9) ) exit boucle_ligne
    end do

  end do boucle_ligne

  backspace(irec)

  Open( 16, status='SCRATCH' )

  do i = 1,Length
    if( mot(i:i) == char(9) ) mot(i:i) = ' '
  end do

  n = 0
  i = 0
  do icol = 1,Length
    i = i + 1
    if( i > Length ) exit
    if( mot(i:i) == ' ' ) cycle
    do j = i+1,Length
      if( mot(j:j) == ' ' ) exit
    end do
    j = j - 1
    n = n + 1
    write(16,'(A)') mot(i:j)
    i = j
  end do

  Rewind(16)

  boucle_i: do i = 1,n
    read(16,*,iostat=eof) x
    if( eof /= 0 ) exit boucle_i
    backspace(16)
    read(16,'(A)',iostat=eof) test
    l = len_trim(test)
    do j = 1,l
      do k = 5,nlet
        if( test(j:j) == let(k) ) exit boucle_i
      end do
    end do
    do k = 1,4
      if( test(1:1) == let(k) .or. test(l:l)==let(k) ) exit boucle_i
    end do
    nmots = nmots + 1
  end do boucle_i
  Close(16)


  nnombre = nmots
  return

end

!***********************************************************************

! Calculation of the integral of f(r), from 0 to Radius
! Integral is calculated using a second order polynomia
! nr: number of points

  function f_integr(r,fct,nr,ir0,nrm,Radius)

  use declarations
  implicit none

  integer:: i, ir0, nr, nrm

  real(kind=db):: a, b, c, d0m, dp0, dpm, f_integr, f0, fm, fp, r0, Radius, rm, rp, tiers, x1, x2
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

! Calcul de L'integrale de 0 a Radius de f.
! L'integrale est calculee avec un polynome d'interpolation d'ordre 3
! r: contient les valeurs des rayons pour les points qu'on integre
! On considere que seule les valeurs jusqu'à L'indice n tel que
! r(i) > Radius > Radius(i-1) sont sures.

  real(kind=db) function f_integr3(r,fct,ir0,nrm,Radius)

  use declarations
  implicit none

  integer:: i, ipr, ir0, nrm

  logical:: This_is_the_end

  real(kind=db):: a, b, c, d, d0m, dp0, dpm, f0, f1, f2, f3, f4, fac1, fac2, fm, fp, r0, Radius, r1, r2, r3, r4, rap14, rap24, &
                  rap34, rm, rp, Tiers
  real(kind=db), dimension(ir0:nrm):: fct, r

  f_integr3 = 0._db

  if( Radius < eps10 ) return

  if( nrm < 4 ) then
    call write_error
    do ipr = 6,9,3
      write(ipr,110) nrm
    end do
    stop
  endif

  Tiers = 1._db / 3
  This_is_the_end = .false.

  do i = 1,nrm-1
    if( i == 1 ) then
      if( ir0 == 1 ) then
 ! Construction de L'ordonnée par extrapolation du second degre
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
      if( i < 3 ) then
        r1 = r(i-1)
        f1 = fct(i-1)
        r2 = r(i)
        f2 = fct(i)
        r3 = r(i+1)
        f3 = fct(i+1)
        r4 = r(i+2)
        f4 = fct(i+2)
      endif
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
  110 format(/' Call to function f_integr3 needs a number of point higher than 3, here it is =',i2,' !')
end

!***********************************************************************

! Calcul de L'integrale de 0 a Radius de f. f(0) = 0.
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

!***********************************************************************

subroutine mult_cell(itape,File_out)

  use declarations
  implicit none

  integer, parameter:: nam = 100000
  integer, parameter:: ntypem = 100

  integer:: eof, eoff, i, igrdat, ipr, itape, ix, iy, iz, j, l, n, na, nnombre, ntype

  real(kind=db):: Thickness
  real(kind=db), dimension(3):: A, abc, abc_new, Angle, Angle_new, B, C, CosAngle, q, V
  real(kind=db), dimension(3,3):: Mat, Mat_i
  real(kind=db), dimension(3,nam):: p

  integer, dimension(3):: nm
  integer, dimension(nam):: itype, Z
  integer, dimension(ntypem):: Numat

  Logical:: Ang, Mat_mul, Surface, Typ

  character(len=2):: Chemical_Symbol
  character(len=9):: keyword
  character(len=132):: identmot, mot
  character(len=Length_name):: File_out

  Ang = .false.
  Mat_mul = .false.
  Surface = .false.
  Typ = .false.
  Thickness = 0._db
  Numat(:) = 0
  nm(:) = 1

  l = len_trim(File_out)
  if( File_out(l-3:l) /= '.txt' ) File_out(l+1:l+4) = '.txt'
  open(2,file = File_out)

! Lecture

  Rewind(itape)

  do igrdat = 1,100000

    read(itape,'(A)',iostat=eoff) mot
    if( eoff /= 0 ) exit

    keyword = identmot(mot,9)
    if( keyword(1:1) == '!' ) cycle

    select case(keyword)

      case('end')
        exit

      case('mult_cell')
        n = nnombre(itape,132)
        if( n == 0 ) call write_err_form(itape,keyword)
        n = min(n,3)
        read(itape,*) nm(1:n)

      case('mat_mul')
        Mat_mul = .true.
        do i  = 1,3
          n = nnombre(itape,132)
          if( n < 3 ) call write_err_form(itape,keyword)
          read(itape,*) Mat(i,:)
        end do
        
      case('surf_cell')
        Surface = .true.
        n = nnombre(itape,132)
        if( n == 0 ) call write_err_form(itape,keyword)
        n = min(n,3)
        if( n == 1 ) then
          read(itape,*) Thickness
        else
          read(itape,*) nm(1:2), Thickness
        endif

      case('atomic_nu')
        Typ = .true.
        ntype = nnombre(itape,132)
        if( ntype > ntypem ) then
          call write_error
            do ipr = 6,9,3
              write(ipr,110) ntype, ntypem
            end do
          stop
        endif
        read(itape,*,iostat=eoff)  Numat(1:ntype)
        if( eoff /= 0 ) call write_err_form(itape,keyword)

      case('unit_cell')
        n = nnombre(itape,132)
        if( n == 3 .or. n == 4 .or. n == 5 ) then
          read(itape,*) abc(:)
        elseif( n >= 6 ) then
          Ang = .true.
          read(itape,*) abc(:), Angle(:)
        else
          call write_err_form(itape,keyword)
        endif

        i = 0
        do j = 1,100000
          if( j <= nam ) i = i + 1
          read(itape,*,iostat=eof) itype(i), p(:,i)
          if( eof /= 0 ) exit
          if( Typ ) then
            Z(i) = Numat(itype(i))
          else
            Z(i) = itype(i)
          endif
        end do

! ESRF changes start
! An "END=" condition leaves the file position after the endfile record,
! thus backspace to be just before the endfile record. This way a
! subsequent READ with another "END=" clause can execute correctly.
        backspace(itape)
! ESRF changes end

        if( j > nam ) then
          call write_error
          do ipr = 6,9,3
            write(ipr,130) j, nam
          end do
          stop
        endif

        na = i-1

      case default

        if( keyword(1:1) /= ' ' ) then
          call write_error
          do ipr = 6,9,3
            write(ipr,140)
            write(ipr,150) mot
          end do
          stop
        endif

    end select

  end do

  if( Mat_mul ) then
    do i = 1,3
      A(i) = Mat(1,i) * abc(i) 
      B(i) = Mat(2,i) * abc(i) 
      C(i) = Mat(3,i) * abc(i)
    end do
    if( Ang ) then
      CosAngle(:) = cos( Angle(:) * pi / 180._db )
      abc_new(1) = sqrt( sum( A(:)**2 ) + 2 * A(1)*A(2) * cosAngle(3) + 2 * A(1)*A(3) * cosAngle(2) &
                                        + 2 * A(2)*A(3) * cosAngle(1) ) 
      abc_new(2) = sqrt( sum( B(:)**2 ) + 2 * B(1)*B(2) * cosAngle(3) + 2 * B(1)*B(3) * cosAngle(2) &
                                        + 2 * B(2)*B(3) * cosAngle(1) ) 
      abc_new(3) = sqrt( sum( C(:)**2 ) + 2 * C(1)*C(2) * cosAngle(3) + 2 * C(1)*C(3) * cosAngle(2) &
                                        + 2 * C(2)*C(3) * cosAngle(1) )
      Angle_new(1) = sum( B(:)*C(:) ) + ( B(1)*C(2) + B(2)*C(1) ) * cosAngle(3) + ( B(1)*C(3) + B(3)*C(1) ) * cosAngle(2) &
                                      + ( B(2)*C(3) + B(3)*C(2) ) * cosAngle(1)
      Angle_new(1) = Angle_new(1) / ( abc_new(2) * abc_new(3) )                                           

      Angle_new(2) = sum( A(:)*C(:) ) + ( A(1)*C(2) + A(2)*C(1) ) * cosAngle(3) + ( A(1)*C(3) + A(3)*C(1) ) * cosAngle(2) &
                                      + ( A(2)*C(3) + A(3)*C(2) ) * cosAngle(1)
      Angle_new(2) = Angle_new(2) / ( abc_new(1) * abc_new(3) )                                           

      Angle_new(3) = sum( A(:)*B(:) ) + ( A(1)*B(2) + A(2)*B(1) ) * cosAngle(3) + ( A(1)*B(3) + A(3)*B(1) ) * cosAngle(2) &
                                      + ( A(2)*B(3) + A(3)*B(2) ) * cosAngle(1)
      Angle_new(3) = Angle_new(3) / ( abc_new(1) * abc_new(2) )
      
      Angle_new(:) = acos( Angle_new(:) ) * 180._db / pi                                           
    else
      abc_new(1) = sqrt( sum( A(:)**2 ) ) 
      abc_new(2) = sqrt( sum( B(:)**2 ) )
      abc_new(3) = sqrt( sum( C(:)**2 ) )
    endif
    
  else
  
    if( Surface ) nm(3) = 1
    abc_new(:) = abc(:) * nm(:)
    Angle_new(:) = Angle(:)
 
  endif
  
  if( Ang ) then
    write(2,120) abc_new(:), Angle_new(:)
    write(6,120) abc_new(:), Angle_new(:)
  else
    write(2,120) abc_new(:)
    write(6,120) abc_new(:)
  endif
    
  if( Mat_mul ) then

    call invermat(Mat,Mat_i)
    
    j = 0
    do i = 1,na

      do ix = -10,10
        v(1) = p(1,i) + ix
        do iy = -10,10
          v(2) = p(2,i) + iy
          do iz = -10,10
            v(3) = p(3,i) + iz
          
            q = Matmul( Mat_i, v )
            
            if( q(1) > 1._db - eps10 .or. q(2) > 1._db - eps10 .or. q(3) > 1._db - eps10 &
             .or. q(1) < - eps10 .or. q(2) < - eps10 .or. q(3) < - eps10 ) cycle
         
            j = j + 1
               
            write(2,160) itype(i), q(:), j, Chemical_Symbol(Z(i))
            write(6,160) itype(i), q(:), j, Chemical_Symbol(Z(i))
            
          end do
        end do
      end do
    end do
    
  else
  
    if( Surface ) nm(3) = int( Thickness / abc(3) ) + 1
    j = 0
    do ix = 0,nm(1)-1
      do iy = 0,nm(2)-1
        do iz = 0,nm(3)-1
          write(2,*)
          do i = 1,na
            j = j + 1
            q(1) = ( p(1,i) + ix ) / nm(1)
            q(2) = ( p(2,i) + iy ) / nm(2)
            if( Surface ) then
              q(3) = p(3,i) + iz
            else
              q(3) = ( p(3,i) + iz ) / nm(3)
            endif
            if( Surface .and. q(3) > Thickness ) cycle
            write(2,160) itype(i), q(:), j, Chemical_Symbol(Z(i))
            write(6,160) itype(i), q(:), j, Chemical_Symbol(Z(i))
          end do

        end do
      end do
    end do

  endif
  
  Close(itape)
  Close(2)

  return

  110 format(//' number of type =',i5,' > ntypem =',i4,// ' Change the parameter ntypem in the code !'//)
  120 format(/' Crystal',/5x,3f16.12,3f12.5)
  130 format(//' number of atoms =',i6,' > nam =',i6,// ' Change the parameter nam in the code !'//)
  140 format(//'  Error in the indata file :')
  150 format(//' The following line is not understood :',/A,// &
          ' If it is a keyword, check the spelling.'/, &
          ' If the line is not supposed to be a keyword but contains numbers, check:'/ &
          5x,' - How many numbers must be in the line ?'/, &
          5x,' - Are there spaces between the numbers ?'/, &
          5x,' - Tabulations are forbidden !'//)
  160 format(i5,3f16.12,'  ! ',i5,3x,a2)
end
