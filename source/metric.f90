! FDMNES subroutines

! Programme de calcul des distances metriques D1 et D2 et du R-factor Rx
! entre deux jeux de spectres.

! Calcule :
!  * le rapport entre les integrales des spectres : "rap",
!  * l'ecart type entre les faisceaux de ces rapports : "Ecart",
!  * les distances metriques D1, D2 et Rx de chaque spectre,
!  * les distances metriques totales D1t et D2t et Rxt.

subroutine metric(comt,convolution_out,Dafs_bio,Dist_min,Dist_min_g,fdmfit_out,fit_cal,Gen_Shift_min,ical, &
             ical_Met_min,index_Met_Fit,iscratchconv,itape_minim,itape2,length_line,Length_block,nb_datafile, &
             ncal,ndm,ng,ngroup_par,nmetric,nmetricm,Nom_Met,npar,nparam,nparm,param,parmax,parmin, &
             RapIntegrT_min_g,typepar)

  use declarations
  implicit none

  integer:: eof, eog, eoh, i, i_data, ical, id, ifich, ig, ig0, igr, ii, im, immax, index_D2, index_Met_Fit, index_Rx, &
    index_Rxg, ipar, ipr, iscratchconv, istat, itape_minim, itape2, j, j0, jpas, k, kk, l, length_line, n, nb_datafile, &
    ncal, ncolm, ndem, ndm, ng, ngroup_par, nmetric, nmetricm, nnombre, nparm, npm, npp

  character(len=8):: dat
  character(len=10):: tim
  character(len=9):: keyword
  character(len=50):: com_date, com_time
  character(len=132):: comt, convolution_out, fdmfit_out, identmot, mot
  character(len=length_line):: motl
  character(len=2), dimension(nmetricm) :: Nom_Met
  character(len=9), dimension(ngroup_par,nparm) :: typepar
  character(len=132), dimension(nb_datafile,2):: File_dat

  integer, dimension(nmetricm) :: ical_Met_min, i_Shift_Met, indparMet, jMet
  integer, dimension(ngroup_par) :: Length_block, npar, nparam
  integer, dimension(:), allocatable:: np, numcol2
  integer, dimension(:,:), allocatable:: npf, npfile, numcol, numlincom

  logical:: Cal_D2, Cal_Rx, Cal_Rxg, Dafs_bio, detail, fit_cal, Fit_Rx, kev, Met_min_g, Simple_comp, Print_all
  logical, dimension(nmetricm) :: Met_min

  real(kind=db):: a, b, c, ci, ci00, de, E0, Emaxdec, Emindec, &
    Energ, fac, p1, p2, pas, pas_shift, Rap_init, Rapg_c, Rapg_min, Rxii, Xanes
  real(kind=db), dimension(2):: Integrt
  real(kind=db), dimension(-50:50) :: Rxi
  real(kind=db), dimension(ng) :: E1, E2, EMetrMax, EMetrMin, RapIntegrT_min_g, Siexp2
  real(kind=db), dimension(nmetricm) :: Dist_min, Dist_min_c, Gen_Shift_min, Gen_Shift_min_c, paramMet
  real(kind=db), dimension(ng,nmetricm) :: Dist_min_g
  real(kind=db), dimension(ngroup_par,nparm) :: param, parmax, parmin
  real(kind=db), dimension(:), allocatable :: decalE, E, EcartType, Poidt, Rapg, RapMoy, Yl
  real(kind=db), dimension(:,:), allocatable :: Dist_Met, Integr, Poids, RapIntegrT, weig
  real(kind=db), dimension(:,:,:), allocatable :: Ef, Met, Y, Yf

!--- Entrees -----------------------------------------------------

  ndm = 1
  Emindec = 0._db; Emaxdec = 0._db; Pas_shift = 0._db
  index_Rxg = 0
  Cal_D2 = .false.
  Cal_Rx = .true.
  Fit_Rx = .false.
  Cal_Rxg = .false.

  kev = .false.
  detail = .false.
  Simple_comp = .false.

  Rewind(itape2)

  allocate( np(ng) )
  allocate( npf(ng,2) )
  allocate( npfile(nb_datafile,2) )
  allocate( numcol(ng,2) )
  allocate( numcol2(ng) )
  allocate( numlincom(ng,2) )
  allocate( weig(ng,2) )
  numcol2(:) = 0

  numlincom(:,:) = 0
  EMetrMin(:) = -1000000.
  EMetrMax(:) = 1000000.

  do ii = 1,1000

    n = nnombre(itape2,132)
    read(itape2,'(A)',iostat=eof) mot
    if( eof /= 0 ) exit

    keyword = identmot(mot,9)

    select case(keyword)

      case('experimen')

        do i_data = 1,nb_datafile
          n = nnombre(itape2,132)
          mot = ' '
          read(itape2,'(A)') mot
          if( mot(1:1) == ' ' ) mot = adjustl( mot )
          File_dat(i_data,1) = mot
          File_dat(i_data,2) = convolution_out
          n = nnombre(itape2,132)
          if( Dafs_bio ) then
            if( n /= 0 ) read(itape2,*)
          else
            ig = i_data
            numcol(ig,:) = 2
            select case(n)
              case(1)
                read(itape2,*) numcol(ig,2)
              case(2)
                read(itape2,*) numcol(ig,:)
              case(3)
                read(itape2,*,iostat=eof) numcol(ig,:), numlincom(ig,1)
                if( eof/= 0 ) call write_err_form(itape2,keyword)
              case(4)
                read(itape2,*) numcol(ig,2), weig(ig,1),numcol2(ig), weig(ig,2)
            end select
            if( numcol(ig,1) < 2 .or. numcol(ig,2) < 2 ) then
              call write_error
              do ipr = 6,9,3
                write(ipr,110)
              end do
              stop
            endif
          endif
        end do

      case('file_met')
        Simple_comp = .true.
        do i = 1,2
          mot = ' '
          read(itape2,'(A)') mot
          if( mot(1:1) == ' ' ) mot = adjustl( mot )
          l = len_trim(mot)
          if( l > 4 ) then
           if( mot(l-3:l-3) /= '.' ) mot(l+1:l+4) = '.txt'
          endif
          File_dat(1,i) = mot
          if( i == 2 ) then
            open(99,file = File_dat(1,i),status = 'old', iostat = istat)
            if( istat /= 0 ) call write_open_error(File_dat(1,i),istat,1)
            read(99,*)
            n = nnombre(99,132000) - 1
            if( n /= ng ) then
              call write_error
              do ipr = 6,9,3
                write(ipr,115) n, ng
              end do
              stop
            endif
            Close(99)
          endif
          do ig = 1,ng
            numcol(ig,:) = 1 + ig
          end do
        end do

      case('gen_shift')

        n = nnombre(itape2,132)
        read(itape2,*,iostat=eof) Emindec, Emaxdec, ndm
        if( eof/= 0 ) call write_err_form(itape2,keyword)

      case('emin')

        n = nnombre(itape2,132)
        if( n == 1 ) then
          Read(itape2,*,iostat=eof) EMetrMin(1)
          do ig = 2,ng
            EMetrMin(ig) = EMetrMin(1)
          end do
        else
          Read(itape2,*,iostat=eof) EMetrMin(1:ng)
        endif
        if( eof/= 0 ) call write_err_form(itape2,keyword)

      case('emax')

        n = nnombre(itape2,132)
        if( n == 1 ) then
          Read(itape2,*,iostat=eof) EMetrMax(1)
          do ig = 2,ng
            EMetrMax(ig) = EMetrMax(1)
          end do
        else
          Read(itape2,*,iostat=eof) EMetrMax(1:ng)
        endif
        if( eof/= 0 ) call write_err_form(itape2,keyword)

      case('kev')

        kev = .true.

      case('rx')

        Cal_Rx = .true.
        Fit_Rx = .true.

      case('detail')

        detail = .true.

      case('fit_out','metric_ou')

        n = nnombre(itape2,132)
        read(itape2,'(A)') fdmfit_out
        fdmfit_out = adjustl( fdmfit_out )
        l = len_trim(fdmfit_out)
        if( l > 4 ) then
          if( fdmfit_out(l-3:l-3) /= '.' ) fdmfit_out(l+1:l+4) = '.txt'
        endif 

      case('rxg')
        Cal_Rxg = .true.

      case('d2')
        Cal_D2 = .true.

      case default

        call write_error
        do ipr = 6,9,3
          write(ipr,120) mot
        end do
        stop

     end select

  end do

  nmetric = 1
  Nom_Met(1) = 'D1'

  if( Cal_D2 ) then
    nmetric = nmetric + 1
    index_D2 = nmetric
    Nom_Met(nmetric) = 'D2'
  endif

  nmetric = nmetric + 1
  index_Rx = nmetric
  Nom_Met(nmetric) = 'Rx'

  if( Cal_Rxg ) then
    nmetric = nmetric + 1
    index_Rxg = nmetric
    Nom_Met(nmetric) = 'Rg'
  endif

  if( Fit_Rx ) then
    index_Met_Fit = index_Rx
  else
    index_Met_Fit = 1
  endif

  if( Dafs_bio ) then
    ncolm = 0
  else
    ncolm = 2
    do ifich = 1,2
      do ig = 1,ng
        ncolm = max( ncolm, numcol(ig,ifich) )
      end do
    end do
    do ig = 1,ng
      ncolm = max( ncolm, numcol2(ig) )
    end do
  endif

  allocate( Dist_Met(ndm,nmetric) )
  allocate( decalE(ndm) )
  allocate( EcartType(ndm) )
  allocate( Met(ng,ndm,nmetric) )
  allocate( Poids(ng,ndm) )
  allocate( Poidt(ndm) )
  allocate( RapIntegrT(ng,ndm) )
  allocate( Rapg(ndm) )
  allocate( RapMoy(ndm) )
  allocate( Yl(ncolm) )

  RapIntegrT(:,:) = 0
  RapMoy(:) = 0

  if( ndm > 1 ) pas_shift = ( Emaxdec - Emindec ) / (ndm - 1 )
  do id = 1,ndm
    decalE(id) = Emindec + ( id - 1 ) * pas_shift
  end do

!--- Lecture -----------------------------------------------------

  npm = 0
  do i_data = 1,nb_datafile

    do ifich = 1,2    ! 1: experience, 2: calcul
      if( fit_cal .and. ifich == 2 ) then
        ipr = iscratchconv
        Rewind(ipr)
      else
        ipr = 2
        open(ipr,file = File_dat(i_data,ifich),status = 'old', iostat = istat)
        if( istat /= 0 ) call write_open_error(File_dat(i_data,ifich),istat,1)
      endif
      if( Dafs_bio .and. ifich == 1 ) then
        read(ipr,*)
        read(ipr,*)
        do i = 1,10000
          Read(ipr,*,iostat=eog) Energ
          if( eog /= 0 ) exit
        end do
        npfile(i_data,ifich) = i - 1
        npm = max( npm, npfile(i_data,ifich) )
      else
        do ii = 1,1000
          n = nnombre(ipr,132)
          if( n < 1 ) then
            read(ipr,*)
          else
            exit
          endif
        end do
        do i = 1,10000
          Read(ipr,*,iostat=eoh) Energ
          if( eoh /= 0 ) exit
        end do
        if( Dafs_bio ) then
          npfile(i_data,ifich) = i - 1
          npm = max( npm, npfile(i_data,ifich) )
        else
          ig = i_data
          if( Simple_comp ) then
            npf(:,ifich) = i - 1
            npm = max( npm, i - 1 )
          else
            npf(ig,ifich) = i - 1
            npm = max( npm, npf(ig,ifich) )
          endif
        endif
      endif
      if( .not. (fit_cal .and. ifich == 2) ) close(ipr)
    end do
  end do

  allocate( E(npm) )
  allocate( Ef(npm,ng,2) )
  allocate( Integr(npm,2) )
  allocate( Y(npm,ng,2) )
  allocate( Yf(npm,ng,2) )

  if( Dafs_bio ) then

! Fichiers experiences
    ig0 = 0
    do i_data = 1,nb_datafile
      ipr = 2
      open(ipr, file=File_dat(i_data,1), status='old', iostat=istat)
      n = nnombre(ipr,100000)
      n = n / 3
      Read(ipr,*)
      Read(ipr,*)
      do i = 1,npfile(i_data,1)
        Read(ipr,*) Energ, Yf(i,ig0+1:ig0+n,1)
        Ef(i,ig0+1:ig0+n,1) = Energ
      end do
      npf(ig0+1:ig0+n,1) = npfile(i_data,1)
      ig0 = ig0 + n
      close(ipr)
    end do

! Fichier calcul
    if( Fit_cal ) then
      ipr = iscratchconv
      Rewind(ipr)
    else
      ipr = 2
      open(ipr, file = File_dat(1,2), status='old', iostat=istat)
    endif
    Read(ipr,*)
    do i = 1,npfile(1,2)
      Read(ipr,*) Energ, xanes, Yf(i,1:ng,2)
      Ef(i,:,2) = Energ
    end do
    npf(:,2) = npfile(1,2)
    close(ipr)

  else

    do ig = 1,ng
      do ifich = 1,2
        if( fit_cal .and. ifich == 2 ) then
          ipr = iscratchconv
          Rewind(ipr)
        else
          ipr = 2
          if( Simple_comp ) then
            open(ipr,file = File_dat(1,ifich),status='old', iostat=istat)
          else
            open(ipr,file = File_dat(ig,ifich),status='old', iostat=istat)
          endif
          if( istat /= 0 ) call write_open_error(File_dat(ig,ifich),istat,1)
        endif
        do i = 1, numlincom(ig,ifich)
          read(ipr,*)
        end do
        do ii = 1,1000
          n = nnombre(ipr,132)
          if( n < 1 ) then
            read(ipr,*)
          else
            exit
          endif
        end do
        do i = 1,npf(ig,ifich)
          Read(ipr,*) Ef(i,ig,ifich), Yl(2:numcol(ig,ifich))
          if( ifich == 2 .and. numcol2(ig) > 0 ) then
            Yf(i,ig,ifich) = weig(ig,1) * Yl(numcol(ig,ifich))
            backspace(ipr)
            Read(ipr,*) Ef(i,ig,ifich), Yl(2:numcol2(ig))
            Yf(i,ig,ifich) = Yf(i,ig,ifich) + weig(ig,2) * Yl(numcol2(ig))
          else
            Yf(i,ig,ifich) = Yl(numcol(ig,ifich))
          endif
        end do
        if( .not. (fit_cal .and. ifich == 2) ) close(ipr)
      end do
    end do

  endif

  if( kev ) Ef(:,:,1) = 1000 * Ef(:,:,1)

  de = emindec - pas_shift
  Ef(:,:,2) = Ef(:,:,2) + de

!-----------------------------------------------------------------------

  do id = 1,ndm

    do ig = 1,ng

      Ef(:,ig,2) = Ef(:,ig,2) + pas_shift
! Bornes de la metrique.
      E1(ig) = max(EMetrMin(ig),Ef(1,ig,1),Ef(1,ig,2))
      E2(ig) = min(EMetrMax(ig),Ef(npf(ig,1),ig,1), Ef(npf(ig,2),ig,2))
      if (E2(ig) <= E1(ig)) then
        Poids(ig,id) = 0._db
      else
        Poids(ig,id) = E2(ig) - E1(ig)
      endif
    end do

    Poidt(id) = sum( Poids(1:ng,id) )
    Poids(1:ng,id) = Poids(1:ng,id) / Poidt(id)

    if( Poidt(id) < 1e-08_db ) then
      call write_error
      do ipr = 6,9,3
        write(ipr,125)
        do ig = 1,ng
          write(ipr,126) ig, Ef(1,ig,1), Ef(1,ig,2), Ef(npf(ig,1),ig,1), Ef(npf(ig,2),ig,2), EMetrMin(ig), EMetrMax(ig)
        end do
      end do
      stop
    endif

    do ig = 1,ng

      if ( E2(ig) < E1(ig) ) then
        Met(ig,id,:) = 0._db
        cycle
      endif

! Elaboration de la grille en energie pour le calcul de la metrique.
! On prend un nombre de point egal a la moyenne entre les nombres
! de points des 2 courbes inclus dans la gamme de calcul.

      np(ig) = nint( 0.5 * ( E2(ig) - E1(ig) ) * ( npf(ig,1) / ( Ef(npf(ig,1),ig,1) - Ef(1,ig,1) ) &
               + npf(ig,2) / ( Ef(npf(ig,2),ig,2) - Ef(1,ig,2) ) ) )
      np(ig) = max(np(ig),4)
      dE = ( E2(ig) - E1(ig) ) / ( np(ig) - 1 )
      E(1) = E1(ig)
      do i = 2,np(ig)
        E(i) = E(i-1) + dE
      end do

! Calcul des intensites sur la grille par interpolation.
      do ifich = 1,2
        do i = 1,np(ig)
          do j = 2,npf(ig,ifich)
            if( ( Ef(j,ig,ifich) >= E(i) .and. Ef(j,ig,ifich) > Ef(j-1,ig,ifich) + 1.e-10_db ) &
             .or. (j == npf(ig,ifich) ) ) exit
          end do
          if( abs( Ef(j,ig,ifich) - Ef(j-1,ig,ifich) ) < eps10 ) then
            p2 = 0.5_db
          else 
            p2 = ( E(i) - Ef(j-1,ig,ifich) ) / ( Ef(j,ig,ifich) - Ef(j-1,ig,ifich) )
          endif
          p1 = 1 - p2
          Y(i,ig,ifich) = p1*Yf(j-1,ig,ifich) + p2*Yf(j,ig,ifich)
       end do
      end do

! Normalisation
      E0 = E(1)
      E(:) = ( E(:) - E0 ) / ( E2(ig) - E1(ig) )

      de = 1.0 / (np(ig) - 1)
      do ifich = 1,2
        Integr(1,ifich) = 0.5 * de * abs( Y(1,ig,ifich) )
        do i = 2,np(ig)-1
          Integr(i,ifich) = Integr(i-1,ifich) + de * abs( Y(i,ig,ifich) )
        end do
        i = np(ig)
        Integr(i,ifich) = Integr(i-1,ifich) + 0.5 * de * abs( Y(i,ig,ifich) )
        IntegrT(ifich) = Integr(i,ifich)
      end do

      RapIntegrT(ig,id) = IntegrT(2) / IntegrT(1)
      RapMoy(id) = RapMoy(id) + Poids(ig,id) * RapIntegrT(ig,id)

! Calcul de Rx., Voir Horsky et al. Phys. Rev. B 46, 7011 (1992) }
      if( Cal_Rx ) then
        siexp2(ig) = sum( Y(1:np(ig),ig,2)**2 )

        Rxi(-3:0) = 0
        j = - 2
        jpas = 1
        pas = 0.1_db * RapIntegrT(ig,id)
        do k = 1,1000
          j0 = j
          j = j + jpas
          ci = RapIntegrT(ig,id) + j * pas
          Rxi(j) = sum((ci*Y(1:np(ig),ig,1) - Y(1:np(ig),ig,2) )**2)
          Rxi(j) = Rxi(j) / siexp2(ig)
          if( j == 0 .and. Rxi(0) > Rxi(-1) ) then
            j = -1
            jpas = -1
          endif
          if( k < 3 ) cycle
          if( (Rxi(j0+1) > Rxi(j0)).and.(Rxi(j0-1) > Rxi(j0)) ) exit
        end do

        ci00 = RapIntegrT(ig,id) + j0 * pas
        a = 0.5*(Rxi(j0+1) + Rxi(j0-1) - 2*Rxi(j0)) / pas**2
        b = ( Rxi(j0+1) - Rxi(j0-1) ) / (2*pas) - 2 * a * ci00
        c = Rxi(j0) - a*ci00**2 - b*ci00
        ci = -0.5*b/a
        Met(ig,id,index_Rx) = a*ci**2 + b*ci + c
        Met(ig,id,index_Rx) = (1.5 + 2.0/3) * Met(ig,id,index_Rx)
      endif

! 3) Normalisation.
      do i = 1,np(ig)
        Y(i,ig,:) = Y(i,ig,:) / IntegrT(:)
        Integr(i,:) = Integr(i,:) / IntegrT(:)
      end do

! Calcul de la distance metrique D1
      Met(ig,id,1) = 0.5 * abs( Y(1,ig,1) - Y(1,ig,2) ) + sum( abs( Y(2:np(ig)-1,ig,1) - Y(2:np(ig)-1,ig,2) ) ) &
               + 0.5 * abs( Y(np(ig),ig,1) - Y(np(ig),ig,2) )
      Met(ig,id,1) = 50 * ( E(2) - E(1) ) * Met(ig,id,1)

      if( Cal_D2 ) then
! Calcul de la distance metrique D2
        Met(ig,id,index_D2) = 0.5 * abs( Integr(np(ig),1) - Integr(np(ig),2) ) &
         + sum( abs( Integr(2:np(ig)-1,1) - Integr(2:np(ig)-1,2) ) )
        Met(ig,id,index_D2) = 100 * ( E(2) - E(1) ) * Met(ig,id,index_D2)
      endif

! 4) Denormalisation.
      do i = 1,np(ig)
        Y(i,ig,:) = Y(i,ig,:) * IntegrT(:)
      end do

    end do

! Calcul de Rxg., Voir Horsky et al. Phys. Rev. B 46, 7011 (1992)
    if( Cal_Rxg ) then
      Rap_init = 10 * Rapmoy(id)
      boucle_kk: do kk = 1,10
        Rap_init = 0.1 * Rap_init
        Rxi(-3:0) = 0
        j = - 2
        jpas = 1
        pas = 0.1_db * Rap_init
        do k = 1,10000
          j0 = j
          j = j + jpas
          ci = Rap_init + j * pas
          Rxi(j) = 0
          do ig = 1,ng
            Rxii = sum((ci*Y(1:np(ig),ig,1) - Y(1:np(ig),ig,2))**2)
            Rxii = Rxii / siexp2(ig)
            Rxi(j) = Rxi(j) + Rxii * Poids(ig,id)
          end do
          if( (j == 0) .and. ( Rxi(0) > Rxi(-1) ) ) then
            j = -1
            jpas = -1
          endif
          if( k < 3 ) cycle
          if((Rxi(j0+1) > Rxi(j0)) .and. (Rxi(j0-1) > Rxi(j0))) exit boucle_kk
          if( j == -9 ) cycle boucle_kk
        end do
      end do boucle_kk

      ci00 = Rap_init + j0 * pas
      a = 0.5*(Rxi(j0+1) + Rxi(j0-1) - 2*Rxi(j0)) / pas**2
      b = ( Rxi(j0+1) - Rxi(j0-1) ) / (2*pas) - 2*a*ci00
      c = Rxi(j0) - a*ci00**2 - b*ci00
      ci = -0.5 * b / a
      Rapg(id) = ci
      Dist_Met(id,index_Rxg) = a*ci**2 + b*ci + c
      Dist_Met(id,index_Rxg) = (1.5/ng + 2.0/3) * Dist_Met(id,index_Rxg)
    endif

  end do

! Calcul des metriques totales et de l'ecart type.
  do id = 1,ndm
    EcartType(id) = sum( ( RapIntegrT(1:ng,id) - RapMoy(id) )**2 )
    EcartType(id) = sqrt( EcartType(id) / ng )
    do im = 1,nmetric
      if( im == index_Rxg ) cycle
      Dist_Met(id,im) = sum( Poids(1:ng,id) * Met(1:ng,id,im) )
    end do
    if( Cal_Rx ) Dist_Met(id,index_Rx) = (1.5/ng + 2./3) * Dist_Met(id,index_Rx) / (1.5 + 2./3)
  end do

! ---- Calcul du minimum et ecriture -----------------------------------

  Met_min(:) = .false.

  do im = 1,nmetric

    if( ical == 1 ) Dist_min(im) = 1000000000._db
    Dist_Min_c(im) = 1000000000._db

    do id = 1,ndm
      if( Dist_Met(id,im) < Dist_min_c(im) ) then
        i_Shift_Met(im) = id
        Dist_min_c(im) = Dist_Met(id,im)
        Gen_Shift_min_c(im) = decalE(id)
        if( im == index_Rxg ) Rapg_c = Rapg(id)
      endif
      if( Dist_Met(id,im) < Dist_min(im) ) then
        Met_min(im) = .true.
        ical_Met_min(im) = ical
        Dist_min(im) = Dist_Met(id,im)
        Gen_Shift_min(im) = decalE(id)
        Dist_min_g(:,im) = Met(:,id,im)
        RapIntegrT_min_g(:) = RapIntegrT(:,id)
        if( im == index_Rxg ) Rapg_min = Rapg(id)
      endif
    end do

  end do

  if( fit_cal .and. Met_min(index_Met_Fit) ) then
    Rewind(iscratchconv)
    open(2, file = convolution_out, iostat=istat )
    if( istat /= 0 ) call write_open_error(convolution_out,istat,1)
    do i = 1,100000
      read(iscratchconv,'(A)',iostat=eof) motl
      if( eof /= 0 ) exit
      l = len_trim(motl)
      write(2,'(A)') motl(1:l)
    end do
    Close(2)
  endif
  if( fit_cal ) Close(iscratchconv)

  if( ical == 1 ) then

    open(4, file = fdmfit_out )
    call date_and_time( date = dat, time = tim )
    com_date = '   Date = ' // dat(7:8) // ' ' // dat(5:6) // ' ' // dat(1:4)
    com_time = '   Time = ' // tim(1:2)  // ' h ' // tim(3:4) // ' mn ' // tim(5:6) // ' s'
    write(4,'(1x,A/A/A)') Revision, com_date, com_time
    if( comt /= ' ') write(4,'(A)') comt

    write(4,130)
    if( fit_cal ) then
      write(4,144)
      do ig = 1,ng
        write(4,140) File_dat(ig,1)
      end do
    else
      do ifich = 1,2
        do i_data = 1,nb_datafile
          write(4,140) File_dat(i_data,ifich)
        end do
        if( ifich == 1 ) write(4,150)
      end do
    endif

    if( ncal > 1 ) then
      write(4,145) ngroup_par - 1
      do igr = 2,ngroup_par
        write(4,146) igr-1, ( adjustr( typepar(igr,ipar) ), ipar = 1,npar(igr) )
        if( Length_block(igr) == 0 ) then
          npp = npar(igr)
        else
          npp = 1
        endif
        write(4,147) ( parmin(igr,ipar), ipar = 1,npp )
        write(4,147) ( parmax(igr,ipar), ipar = 1,npp )
        write(4,148) ( nparam(igr), ipar = 1,npp )
      end do
    endif

  else

    open( 4, file = fdmfit_out, position='append')

  endif

  Met_min_g = .false.
  do im = 1,nmetric
    if( Met_min(im) ) Met_min_g = .true.
  end do
  Print_all = ncal <= 3600 .or. Met_min_g .or. detail

  do ipr = 4,6,3

    if( ncal > 1 ) then
      if( ical == ncal+1 ) then
        write(ipr,155)
      else
        write(ipr,160) ical, ncal
      endif
      do igr = 2,ngroup_par
        if( nparam(igr) < 2 ) cycle
        if( ipr == 6 ) then
          npp = 1
        else
          npp = npar(igr)
        endif
        write(ipr,180) ( adjustr(typepar(igr,ipar)), param(igr,ipar), ipar = 1,npp )
      end do
    endif

    write(ipr,*)
    do im = 1,nmetric
      if( Met_min(im) .and. .not. ncal == 1 ) then
        write(ipr,190) Nom_Met(im), Dist_min_c(im), Gen_Shift_min_c(im)
      else
        write(ipr,200) Nom_Met(im), Dist_min_c(im), Gen_Shift_min_c(im)
      endif
      if( im == index_Rxg ) write(ipr,205) Rapg_c
    end do

  end do

  if( Cal_Rxg ) then
    immax = nmetric - 1
  else
    immax = nmetric
  endif
  if( ng > 1 .and. Print_all ) then

    if( immax == 3 ) then
      write(4,210) ( Nom_Met(im), decalE(i_Shift_Met(im)), im = 1,immax )
      do ig = 1,ng
        write(4,220) ig, ( Met(ig,i_Shift_Met(im),im), im = 1,immax), RapIntegrT(ig,i_Shift_Met(1))
      end do
    else
      write(4,215) ( Nom_Met(im), decalE(i_Shift_Met(im)), im = 1,immax )
      do ig = 1,ng
        write(4,225) ig, ( Met(ig,i_Shift_Met(im),im), im = 1,immax), RapIntegrT(ig,i_Shift_Met(1))
      end do
    endif
  endif

  if( detail ) then
    if( Cal_Rxg ) then
      if( Cal_D2 ) then
        write(4,310) Nom_Met(1:nmetric)
      else
        write(4,315) Nom_Met(1:nmetric)
      endif
      do id = 1,ndm
        if( Cal_D2 ) then
          write(4,320) decalE(id), Dist_Met(id,1:nmetric), Poidt(id), RapMoy(id), Rapg(id), EcartType(id)
        else
          write(4,325) decalE(id), Dist_Met(id,1:nmetric), Poidt(id), RapMoy(id), Rapg(id), EcartType(id)
        endif
      end do
    else
      if( Cal_D2 ) then
        write(4,330) Nom_Met(1:nmetric)
      else
        write(4,335) Nom_Met(1:nmetric)
      endif
      do id = 1,ndm
        if( Cal_D2 ) then
          write(4,340) decalE(id), Dist_Met(id,1:nmetric), Poidt(id), RapMoy(id), EcartType(id)
        else
          write(4,345) decalE(id), Dist_Met(id,1:nmetric), Poidt(id), RapMoy(id), EcartType(id)
        endif
      end do
    endif
  endif

  if( ng > 1 .and. detail ) then
    do ig = 1,ng
      write(4,350) ig
      if( Dafs_bio ) then
        if( ig == 1 ) then
          i_data = 1
          do ifich = 1,2
            write(4,140) File_dat(i_data,ifich)
          end do
        endif
      else
        i_data = ig
        do ifich = 1,2
          write(4,140) File_dat(i_data,ifich)
        end do
      endif
      if( Cal_D2 ) then
        write(4,360) Nom_Met(1:immax)
      else
        write(4,365) Nom_Met(1:immax)
      endif
      do id = 1,ndm
        if( Cal_D2 ) then
          write(4,370) decalE(id), Met(ig,id,1:immax), Poids(ig,id), RapIntegrT(ig,id)
        else
          write(4,375) decalE(id), Met(ig,id,1:immax), Poids(ig,id), RapIntegrT(ig,id)
        endif
      end do
    end do
  end if

  if( ncal > 1 .and. ical == ncal ) then

    do ipr = 4,6,3
      write(ipr,450) Nom_Met(1:nmetric)
      if( Cal_Rxg ) then
        if( Cal_D2 ) then
          write(ipr,460) Dist_min(1:nmetric), Rapg_min
        else
          write(ipr,462) Dist_min(1:nmetric), Rapg_min
        endif
      else
        if( Cal_D2 ) then
          write(ipr,465) Dist_min(1:nmetric)
        else
          write(ipr,467) Dist_min(1:nmetric)
        endif
      endif
      write(ipr,470) adjustr(typepar(1,1)), Gen_Shift_min(1:nmetric)
    end do
    ndem = 1
    do igr = 2,ngroup_par
      jMet(1:nmetric) = ( ical_Met_min(1:nmetric) - 1 ) / ndem
      ndem = ndem * nparam(igr)
      indparMet(1:nmetric) = mod( jMet(1:nmetric), nparam(igr) ) + 1
      if( Length_block(igr) == 0 ) then
        npp = npar(igr)
      else
        npp = 1
      endif
      do ipar = 1,npp
        if( nparam(igr) < 2 ) then
          paramMet(1:nmetric) = parmin(igr,ipar)
        else
          fac = ( parmax(igr,ipar) - parmin(igr,ipar) ) / ( nparam(igr) - 1 )
          paramMet(1:nmetric) = parmin(igr,ipar) + fac * (indparMet(1:nmetric) - 1)
        endif
        do ipr = 4,6,3
          write(ipr,470) adjustr(typepar(igr,ipar)), paramMet(1:nmetric)
        end do

      end do
    end do

    if( ng > 1 ) then
      do ipr = 4,6,3
        if( immax == 3 ) then
          write(ipr,210) ( Nom_Met(im), Gen_Shift_min(im), im = 1,immax )
          do ig = 1,ng
            write(ipr,220) ig, ( Dist_min_g(ig,im), im = 1,immax), RapIntegrT_min_g(ig)
          end do
        else
          write(ipr,215) ( Nom_Met(im), Gen_Shift_min(im), im = 1,immax )
          do ig = 1,ng
            write(ipr,225) ig, ( Dist_min_g(ig,im), im = 1,immax), RapIntegrT_min_g(ig)
          end do
        endif
      end do
    endif

  endif

  do id = 1,ndm
    write(itape_minim,'(4f19.13)') Dist_Met(id,1:nmetric)
  end do

  Close(4)

  deallocate( decalE )
  deallocate( Dist_Met )
  deallocate( E )
  deallocate( Ef )
  deallocate( EcartType )
  deallocate( Integr )
  deallocate( Met )
  deallocate( np )
  deallocate( npf, npfile )
  deallocate( numcol )
  deallocate( numcol2 )
  deallocate( numlincom )
  deallocate( Poids )
  deallocate( Poidt )
  deallocate( RapIntegrT )
  deallocate( Rapg )
  deallocate( RapMoy )
  deallocate( weig )
  deallocate( Y )
  deallocate( Yf )
  deallocate( Yl )

  return
  110 format(///' Number of column must be greater than 1, first is', ' the energy'///)
  115 format(///' There is not the same number of columns in the 2 files:',i4,' and',i4,' !'///)
  120 format(///' Unknown keyword in the indata file :'//1x,A)
  125 format(///' In then metric calculation, there is no overlap',/ &
  ' between the energy ranges of the experiment and the', ' calculation !'/,' Check the keywords Emin, Emax, the energy', &
  ' range (photon or photoelectron energy), the unit (eV, keV', ' ...'//,' ig   Emin(1)   Emin(2)', &
  '   Emax(1)   Emax(2)   E_met_min    E_met_max ')
  126 format(i3,4f10.3,2f13.3)
  130 format(/' Metric distance calculation between :')
  140 format(a132)
  144 format('  the actual calculation and the experimental files :')
  145 format(/' Number of parameter =',i3)
  146 format(/' Parameter',i2,' / first value / last value /', ' number of value',/,3x,100(2x,a9))
  147 format(3x,100f11.5)
  148 format(3x,100i11)
  150 format('  and :')
  155 format(//80('-')//' Calculation with optimized parameters :')
  160 format(/3x,'Calculation',i6,' /',i6,10x,' Parameters :')
  180 format(3x,100(2x,a9,' =',f11.5,','))
  190 format(3x,a2,' =',f9.5,',  general shift =',f9.2,' eV,', '   ... up to now, best value !')
  200 format(3x,a2,' =',f9.5,',  general shift =',f9.2,' eV')
  205 format(17x,'  General ratio =',1p,e13.5)
  210 format(/'   Spectra  ',3(2x,a2,'(',f9.2,')'),7x,'Rap')
  215 format(/'   Spectra  ',2(2x,a2,'(',f9.2,')'),7x,'Rap')
  220 format(i8,2f15.5,f15.7,1p,e17.3)
  225 format(i8,f15.5,f15.7,1p,e17.3)
  310 format(/4x,'Gen_Shift',4(5x,a2,3x), '    Range     <Rap>       Rapg      Ecart')
  315 format(/4x,'Gen_Shift',3(5x,a2,3x), '    Range     <Rap>       Rapg      Ecart')
  320 format(f13.2,2f10.5,2f10.6,f10.2,1p,3e11.2)
  325 format(f13.2,f10.5,2f10.6,f10.2,1p,3e11.2)
  330 format(/4x,'Gen_Shift',3(5x,a2,3x), '    Range     <Rap>      Ecart')
  335 format(/4x,'Gen_Shift',2(5x,a2,3x), '    Range     <Rap>      Ecart')
  340 format(f13.2,2f10.5,f10.6,f10.2,1p,2e11.2)
  345 format(f13.2,f10.5,f10.6,f10.2,1p,2e11.2)
  350 format(/'  ig =',i3,' :')
  360 format(/'    Gen_Shift',5x,a2,6x,a3,5x,a2,'   Weights    Rap')
  365 format(/'    Gen_Shift',5x,a2,6x,a2,'   Weights    Rap')
  370 format(f13.2,f8.2,f8.3,f8.4,f8.3,1p,e11.2)
  375 format(f13.2,f8.2,f8.4,f8.3,1p,e11.2)
  450 format(//80('-')//,'  Minimum in the grid of parameters :',// 18x,4(4x,a2,'min',3x))
  460 format(16x,2f12.5,2f12.6,'   Rapg =',1p,e13.5/)
  462 format(16x,f12.5,2f12.6,'   Rapg =',1p,e13.5/)
  465 format(16x,2f12.5,f12.6/)
  467 format(16x,f12.5,f12.6/)
  470 format(5x,a9,' =',4f12.5)
end
