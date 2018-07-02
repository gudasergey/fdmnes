! FDMNES subroutines

! Post treatment of the spectra calculated by FDMNES
! 1) Sommation with energy shift of varius spectra
! 2) Convolution by a lorentzienne for XANES
!              L(x) = (1/(pi*b)) * 1 / ( 1 + ( (x-a)/b )**2 )
!    Integration over energy for Dafs.

subroutine Convolution(bav_open,Bormann,Conv_done,convolution_out,Delta_edge,E_cut_imp,E_cut_man,Ecent,Elarg,Estart, &
        Fit_cal,Gamma_hole,Gamma_hole_imp,Gamma_max,ical, &
        icheck,indice_par,iscratchconv,itape1,kw_conv,Length_line,n_col_max, &
        ngamh,ngroup_par,nkw_conv,nomfich,nomfichbav,npar,nparm,param,Scan_a,typepar,ncal)

  use declarations
  implicit none

  integer:: Dafs_exp_type, eof, i, i_col, i_conv, i_hk, i_Trunc, i1, ical, icheck, ie, ie1, ie2, ifich, igr, ii, &
    index_hk, initl, ip, ipar, ipas, ipl, ipr, ipr1, ipr2, is, iscr, iscratchconv, istop, istat, itape1, j, &
    jfich, jpl, js, jseuil, k, kpl, l, Length_line, ll, mfich, n, n_col, n_col_max, n_energ_tr, &
    n_mat_pol, n_selec_core, n_signal, n_Stokes, n_Trunc, &
    ncal, ne2, nef, nelor, nen2, nenerg, nenerge, nes, nes_in, nfich, &
    ngamh, ngroup_par, ninit, ninit1, ninitlm, nkw_conv, nnombre, np_stokes, nparm, nphim, npldafs, npldafs_b, &
    npldafs_t, npldafs_th, nseuil, nxan, nw

! njp : Points beyond the energy range to make the border effect less strong
  integer, parameter:: njp = 500

  integer, dimension(3):: hkl_S
  integer, dimension(10):: num_core
  integer, dimension(ngroup_par):: npar
  integer, dimension(ngroup_par,nparm):: indice_par
  integer, dimension(:), allocatable:: i_done, indf, ne, ninitl, nphi, numat
  integer, dimension(:,:), allocatable:: hkl_dafs, nsup, ne_initl
  integer, dimension(:), allocatable:: n_index_hk, n_div_fpp

  character(len=6):: mot6, mot6_b
  character(len=9):: keyword, mot9, Traduction
  character(len=15):: mot15
  character(len=132):: identmot, mot, mots
  character(len=Length_word):: nomab   
  character(len=Length_name):: chemin, convolution_out, fichscanout, nomfich, nomfichbav
  character(len=9), dimension(nkw_conv) :: kw_conv
  character(len=9), dimension(ngroup_par,nparm) :: typepar
  character(len=Length_word), dimension(:), allocatable:: nom_col
  character(len=Length_name), dimension(:), allocatable:: Convolution_out_all, fichin, fichscanin, fichscanout_all
  character(len=13), dimension(:), allocatable:: Stokes_name

  complex(kind=db):: cf, zero_c
  complex(kind=db), dimension(1):: cdum
  complex(kind=db), dimension(:), allocatable:: dampl, dph, dph_t, dpht, f0, f0_bulk, f0_th
  complex(kind=db), dimension(:,:), allocatable:: f0scan, f0scan_bulk, phdtscan, Trs
  complex(kind=db), dimension(:,:,:), allocatable:: Ad, Adafs, As, Mu_m, Mu_mat_comp, Mu_t
  complex(kind=db), dimension(:,:,:,:), allocatable:: As_bulk, mu, mus

  logical:: Abs_before, Abs_in_bulk, Analyzer, Another_one, Arc, bav_open, Bormann, Dafs, Dafs_bio, Check_conv, chem, Circular, &
    Conv_done, Cor_abs, decferm, Deuxieme, Double_cor, E_cut_man, Energphot, Epsii_ref_man, Extrap, Fermip, First_E, Fit_cal, &
    Forbidden, fprim, fprime_atom, Full_self_abs, Gamma, Gamma_hole_imp, Gamma_var, Gaussian_default, Green_int, Just_total, &
    Magn, no_extrap, nxan_lib, Photoemission, Scan_a, scan_true, Seah, Self_abs, &
    Signal_Sph, Stokes, Stokes_Dafs, Stokes_xan, Sup_sufix, Tenseur, Tenseur_car, Thomson, Transpose_file, U_iso_inp
  logical, dimension(:), allocatable:: run_done, Skip_run, Trunc

  real(kind=db):: a, a1, a2, a3, a4, Abs_U_iso_inp, alambda, Asea, b, b1, b2, b3, b4, bba, bbb, c, c_micro, &
    conv_mbarn_nelec, ct_epsilon, ct_nelec, d, d_dead, de_obj, de1, de2, Delta_edge, Deltar, &
    E, E_cut_imp, E_obj, e1m, Ecent, E_cut, E_cut_orig, Eintmax, Elarg, Eph, Epsii_ref, Esmin, &
    Estart, f0_forward, fac, fpp0, Gamm, Gamma_h, Gamma_max, &
    Im_pi, Im_sig, Ip_pi, Ip_sig, mu_0_bulk, natomsym, p1, p2, pasdeb, Pdt, Pdt_bulk, S0_2, Sample_thickness, Shift_U_iso,  &
    Surface_ref, Tab_width, Vibration, Volume_maille, Volume_maille_bulk

  real(kind=db), dimension(0):: rdum
  real(kind=db), dimension(3):: angxyz, axyz
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(ngroup_par,nparm) :: param
  real(kind=db), dimension(:), allocatable:: Abs_U_iso, angle, bb, betalor, decal, e1, e2, Efermip, Elor, En_fermi, Energ, &
        Energe, Energ_tr, Ep, Eph1, Ephm, Ephoton, Es, Es_temp, Eseuil, fpp_avantseuil, fi, fr, l_dafs, Length_abs, lori, & 
        lorix, lorr, lorrx, natomsym_f, Pds, p1f, p2f, Tens, V0muf, Ts, Yr, Yi
  real(kind=db), dimension(:,:), allocatable:: decal_initl, Epsii, Length_rel, Mu_tt, mua_r, mua_i, Signal, Stokes_param, Xa, &
                                               Xanes, Xs
  real(kind=db), dimension(:,:,:), allocatable:: Icirc, Icirccor, Icor, Icircdcor, Idcor, Mu_mat, Mus_mat

! Put to .true. when TDDFT
  Another_one = .false.

  do  ! end of loop at the end of the routine

  Abs_before = .false.
  Abs_in_bulk = .false.
  Abs_U_iso_inp = 0._db
  Analyzer = .true.
  Arc = .true.
  Asea = 0.2_db  ! Slope of Gamma at the origin in Seah Dench model
  chem = .false.
  Check_conv = .false.
  Circular = .false.
  Convolution_out = ' '
  d_dead = 0._db
  Dafs = .false.
  Dafs_exp_type = 4
  decferm = .false.
  Deltar = 0._db
  Volume_maille = 0._db
  Volume_maille_bulk = 0._db
  Deuxieme = .false.
  Double_cor = .false.
  E_cut = E_cut_imp * rydb
  Epsii_ref_man = .false.
  eintmax = 1000000._db
  f0_forward = 0._db
  Fermip = .false.
  fichscanout = ' '
  Forbidden = .false.
  fprim = .false.
  fprime_atom = .false.
  Full_self_abs = .false.
  Gamma_var = .false.
  Gaussian_default = .false.
  Green_int = .false.
  hkl_S(:) = 0
  jseuil = 1
  Just_total = .true.
  Magn = .false.
  n_Stokes = 0
  nelor = 0
  nfich = 0
  npldafs = 0
  npldafs_t = 0
  no_extrap = .false. 
  nseuil = -1
  do i = 1,10
    num_core(i) = i
  end do
  nxan_lib = .false.
  Photoemission = .false.
  S0_2 = 1._db
  Sample_thickness = 0._db
  scan_true = .false.
  seah = .false.
  self_abs = .false.
  Signal_Sph = .false.
  Shift_U_iso = 0._db
  Stokes = .false.
  Stokes_dafs = .false.
  Stokes_xan = .false.
  Surface_ref = 0._db
  Tenseur = .false.
  Tenseur_car = .false.
  Thomson = .false.
  Transpose_file = .false.
  U_iso_inp = .false.
  Vibration = 0._db

  if( bav_open .or. icheck > 1 ) Check_conv = .true.
  
! -- Lecture --------------------------------------------

  Rewind(itape1)

  boucle_ii: do ii = 1,1000

    n = nnombre(itape1,132)
    read(itape1,'(A)',iostat=eof) mot
    if( eof /= 0 ) exit boucle_ii

    keyword = identmot(mot,9)

    if( keyword == 'calculati' ) then
      nfich = 0
      boucle_i: do i = 1,1000
        n = nnombre(itape1,132)
        read(itape1,'(A)',iostat=eof) mots
        if( eof /= 0 ) exit boucle_ii
        mot9 = identmot(mots,9)
        do j = 1,nkw_conv
          if( mot9 /= kw_conv(j) ) cycle
          backspace(itape1)
          exit boucle_i
        end do
        if( mot9 == 'run_done' ) exit boucle_ii
        if( n == 0 ) nfich = nfich + 1
      end do boucle_i
    endif

    if( keyword == 'stokes' ) then
      do i = 1,1000
        n = nnombre(itape1,132)
        if( n == 0 ) exit
        n_stokes = n_stokes + 1
        read(itape1,*)
      end do
    endif

    if( keyword == 'circular' ) n_stokes = n_stokes + 4

  end do boucle_ii

  rewind(itape1)

  allocate( run_done(nfich) )
  allocate( skip_run(nfich) )
  allocate( i_done(nfich) )
  allocate( Stokes_param(5,n_Stokes) )
  allocate( Stokes_name(n_Stokes) )
  allocate( V0muf(nfich) )

  Stokes_param(1:3,:) = 0._db
! Analyzer: when the 2 angles are zero, no analyzer
! Rotation angle de of the analyzer
  Stokes_param(4,:) = 0._db
! Bragg angle of the analyzer
  Stokes_param(5,:) = 0._db
  Stokes_name(:) = 'no_name'

  run_done(:) = .false.
  Skip_run(:) = .false.
  do ifich = 1,nfich
    i_done(ifich) = ifich
  end do
  mfich = nfich

  boucle_jj: do ii = 1,1000
    n = nnombre(itape1,132)
    read(itape1,'(A)',iostat=eof) mot
    if( eof /= 0 ) exit boucle_jj
    keyword = identmot(mot,9)
    k = 0
    if( keyword == 'run_done' ) then
      do i = 1,mfich
        read(itape1,*,iostat=eof) j, is
        if( eof /= 0 ) exit boucle_jj
        if( j == 0 ) then
          run_done(i) = .true.
        else
          k = k + 1
          if( is == 0 ) Skip_run(k) = .true.
        endif
        if( run_done(i) ) then
          nfich = nfich - 1
          do j = i+1,mfich
            i_done(j) = i_done(j) - 1
          end do
        endif
      end do
      exit
    endif
  end do boucle_jj

  rewind(itape1)

  if( nfich == 0 ) then
    nfich = 1
    mfich = 1
    allocate( fichin(nfich) )
    mot = ' '
    mot = nomfich
    l = len_trim(mot)
    mot(l+1:l+4) = '.txt'
    fichin(1) = mot
  else
    allocate( fichin(nfich) )
  endif

  allocate( decal(nfich) )
  allocate( Efermip(nfich) )
  allocate( En_fermi(nfich) )
  allocate( Eseuil(nfich) )
  allocate( fichscanin(nfich) )
  allocate( fpp_avantseuil(nfich) )
  allocate( n_div_fpp(nfich) )
  allocate( ne(nfich) )
  allocate( ninitl(nfich) )
  allocate( numat(nfich) )
  allocate( Pds(nfich) )

  Pds(1) = 1._db
  decal(:) = 0._db
  Efermip(1) = 0._db
  ninitl(:) = 1
  ninitlm = 1
  n_selec_core = 0

  if( ncal > 1 .and. ical == 1 ) then
    do jfich = 1,mfich
      iscr = 100 + jfich
      open( iscr, status = 'scratch' )
    end do
  endif

  boucle_lect: do ii = 1,1000

    n = nnombre(itape1,132)
    read(itape1,'(A)',iostat=eof) mot
    if( eof /= 0 ) exit boucle_lect

    keyword = identmot(mot,9)
    keyword = Traduction(keyword)

    select case(keyword)

      case('run_done')
        exit

      case('abs_befor')
        Abs_before = .true.

      case('all_conv')
        Just_total = .false.

      case('calculati')

        ifich = 0
        do i = 1,mfich
          n = nnombre(itape1,132)
          if( run_done(i) ) then
            read(itape1,*)
            n = nnombre(itape1,132)
            if( n /= 0 ) read(itape1,*)
            cycle
          endif
          ifich = ifich + 1
          read(itape1,'(A)') fichin(ifich)
          fichin(ifich) = adjustl( fichin(ifich) )
          n = nnombre(itape1,132)
          if( n == 0 ) then
            Pds(ifich) = 1._db
            decal(ifich) = 0._db
          elseif( n == 1 ) then
            read(itape1,*) Pds(ifich)
            decal(ifich) = 0._db
          elseif( n == 2 ) then
            read(itape1,*) Pds(ifich), decal(ifich)
          else
            Fermip = .true.
            read(itape1,*) Pds(ifich), decal(ifich), Efermip(ifich)
          endif
        end do

      case('cal_tddft')

        ifich = 0
        do i = 1,mfich
          n = nnombre(itape1,132)
          if( run_done(i) ) then
            read(itape1,*)
            n = nnombre(itape1,132)
            if( n /= 0 ) read(itape1,*)
            cycle
          endif
          ifich = ifich + 1
          if( Another_one ) then
            read(itape1,'(A)') fichin(ifich)
            fichin(ifich) = adjustl( fichin(ifich) )
            Deuxieme = .true.
          else
            read(itape1,'(A)') mot
          endif
          Pds(ifich) = 1._db
          decal(ifich) = 0._db
        end do
        Another_one = .not. Another_one

      case('conv_out')

        n = nnombre(itape1,132)
        read(itape1,'(A)') convolution_out
        convolution_out = adjustl( convolution_out )

      case('dafs_exp_')

        n = nnombre(itape1,132)
        read(itape1,*) Dafs_exp_type

      case('gaussian')

        n = nnombre(itape1,132)
        if( n == 1 ) then
          read(itape1,*) deltar
          vibration = 0._db
        elseif( n == 2 ) then
          read(itape1,*) deltar, vibration
        else
          Gaussian_default = .true.
        endif

      case('fprime')

        fprim = .true.

      case('forbidden')

        Forbidden = .true.

      case('scan_file')

        scan_true = .true.
        do ifich = 1,nfich
          n = nnombre(itape1,132)
          read(itape1,'(A)') fichscanin(ifich)
          fichscanin(ifich) = adjustl( fichscanin(ifich) )
        end do

      case('scan')

        Scan_a = .true.

      case('scan_conv')

        n = nnombre(itape1,132)
        read(itape1,'(A)') fichscanout
        fichscanout = adjustl( fichscanout )

      case('directory')

        chem = .true.
        n = nnombre(itape1,132)
        read(itape1,'(A)') chemin

      case('abs_b_iso')
        U_iso_inp = .true.
        n = nnombre(itape1,132)
        if( n == 1 ) then
          read(itape1,*) Abs_U_iso_inp
        else
          read(itape1,*) Abs_U_iso_inp, Shift_U_iso          
        endif
        Abs_U_iso_inp = Abs_U_iso_inp / ( 8 * pi**2 ) 

      case('abs_u_iso')
        U_iso_inp = .true.
        n = nnombre(itape1,132)
        if( n == 1 ) then
          read(itape1,*) Abs_U_iso_inp
        else
          read(itape1,*) Abs_U_iso_inp, Shift_U_iso          
        endif

      case('seah')

        Arc = .false.
        seah = .true.
        n = nnombre(itape1,132)
        select case(n)
          case(0)
            continue
          case(1)
            read(itape1,*) asea
          case(2)
            read(itape1,*) asea, Gamma_max
            Gamma_max = Gamma_max / rydb
          case(3)
            read(itape1,*) asea, Gamma_max, Gamma_hole(1)
            Gamma_hole(1) = Gamma_hole(1) / rydb
            Gamma_max = Gamma_max / rydb
            Gamma_hole_imp = .true.
          case default
            read(itape1,*) asea, Gamma_max, Gamma_hole(1), E_cut
            Gamma_hole(1) = Gamma_hole(1) / rydb
            Gamma_max = Gamma_max / rydb
            Gamma_hole_imp = .true.
            E_cut_man = .true.
        end select

      case('convoluti')

        Arc = .true.
        n = nnombre(itape1,132)
        select case(n)
          case(0)
            continue
          case(1)
            read(itape1,*) Ecent
            Ecent = Ecent / rydb
            Elarg = Elarg / rydb
          case(2)
            read(itape1,*) Ecent, Elarg
            Ecent = Ecent / rydb
            Elarg = Elarg / rydb
          case(3)
            read(itape1,*) Ecent, Elarg, Gamma_max
            Ecent = Ecent / rydb
            Elarg = Elarg / rydb
            Gamma_max = Gamma_max / rydb
          case(4)
            read(itape1,*) Ecent, Elarg, Gamma_max, Gamma_hole(1)
            Ecent = Ecent / rydb
            Elarg = Elarg / rydb
            Gamma_hole(1) = Gamma_hole(1) / rydb
            Gamma_max = Gamma_max / rydb
            Gamma_hole_imp = .true.
          case default
            read(itape1,*) Ecent, Elarg, Gamma_max, Gamma_hole(1), E_cut
            Ecent = Ecent / rydb
            Elarg = Elarg / rydb
            Gamma_hole(1) = Gamma_hole(1) / rydb
            Gamma_max = Gamma_max / rydb
            Gamma_hole_imp = .true.
            E_cut_man = .true.
        end select

      case('table')

        Arc = .false.
        n = nnombre(itape1,132)
        if( n == 1 ) read(itape1,*)
        nelor = 0
        do ie = 1,10000
          n = nnombre(itape1,132)
          if( n == 0 ) exit
          read(itape1,*)
          nelor = nelor + 1
        end do
        rewind(itape1)
        do i = 1,10000
          read(itape1,'(A)') mots
          mot9 = identmot(mots,9)
          if( mot9 == 'table' ) exit
        end do
        n = nnombre(itape1,132)
        if( n == 1 ) then
          read(itape1,*) E_cut
          E_cut_man = .true.
        endif
        allocate( Elor(nelor) )
        allocate( betalor(nelor) )
        do ie = 1,nelor
          n = nnombre(itape1,132)
          read(itape1,*,iostat=eof) Elor(ie), betalor(ie)
          if( eof > 0 ) call write_err_form(itape1,keyword)
        end do
        Elor(:) = Elor(:) / rydb
        betalor(:) = betalor(:) / rydb

      case('eintmax')

        n = nnombre(itape1,132)
        read(itape1,*,iostat=eof) eintmax
        if( eof > 0 ) call write_err_form(itape1,keyword)

      case('gamma_fix')
        gamma_var = .false.

      case('gamma_var')
        gamma_var = .true.

! To make the shift before the convolution.
      case('dec')

        decferm = .true.

      case('no_extrap')

        no_extrap = .true.

      case('nxan_lib')

        nxan_lib = .true.

      case('thomson')

        Thomson = .true.
        npldafs_th = nnombre(itape1,Length_line) / 2
        allocate( fr(npldafs_th) )
        allocate( fi(npldafs_th) )
        allocate( f0_th(npldafs_th) )
        read(itape1,*,iostat=eof) (fr(ipl), fi(ipl),ipl = 1,npldafs_th)
        if( eof > 0 ) call write_err_form(itape1,keyword)
        f0_th(:) = cmplx( fr(:), fi(:),db)
        deallocate( fr )
        deallocate( fi )

      case('transpose')
        Transpose_file = .true.
        n_energ_tr = nnombre(itape1,132)
        if( n_energ_tr > 0 ) then 
          allocate( Energ_tr(n_energ_tr) )
          read(itape1,*) Energ_tr(:)
        endif
         
      case('s0_2')

        read(itape1,*,iostat=eof) S0_2

      case('photo_emi')

        Photoemission = .true.

      case('epsii')
        n = nnombre(itape1,132)
        read(itape1,*,iostat=eof) Epsii_ref
        if( eof > 0 ) call write_err_form(itape1,keyword)
        Epsii_ref_man = .true.

      case('selec_cor')
        n_selec_core = nnombre(itape1,132)
        read(itape1,*,iostat=eof) num_core(1:n_selec_core)
        if( eof > 0 ) call write_err_form(itape1,keyword)

      case('surface_p')
        n = nnombre(itape1,132)
        read(itape1,*,iostat=eof) hkl_S(:)
        if( eof > 0 ) call write_err_form(itape1,keyword)

      case('circular')
        Circular = .true.

      case('sample_th')
        n = nnombre(itape1,132)
        read(itape1,*,iostat=eof) Sample_thickness

      case('stokes')
        Stokes = .true.
        do i = 1,n_stokes
          n = nnombre(itape1,132)
          if( n == 0 ) exit
          read(itape1,*,iostat=eof) Stokes_param(1:min(5,n),i)
          if( eof > 0 ) call write_err_form(itape1,keyword)
!Si l'angle de rotation de l'analyseur est defini mais pas son angle de Bragg,
!ce dernier est supposé = 45 degres correspondant a un analyseur "parfait"
          if( n == 4 ) Stokes_param(5,i) = 45._db
          Stokes_param(4:5,i) = Stokes_param(4:5,i) * pi / 180
        end do

      case('stokes_na')
        boucle_n: do i = 1,n_stokes
          read(itape1,'(A)',iostat=eof) mots
          if( eof /= 0 ) exit boucle_lect
          mot9 = identmot(mots,9)
          do j = 1,nkw_conv
            if( mot9 /= kw_conv(j) ) cycle
            backspace(itape1)
            exit boucle_n
          end do
          mots = adjustl( mots )
          Stokes_name(i) = mots(1:13)
          if( Stokes_name(i) == 'noname       ') Stokes_name(i) = 'no_name'
        end do boucle_n

      case('dead_laye')
        n = nnombre(itape1,132)
        read(itape1,*,iostat=eof) d_dead
        if( eof > 0 ) call write_err_form(itape1,keyword)
        d_dead = d_dead / 10000

      case('double_co')
        Double_cor = .true.

      case('no_analyz')
        Analyzer = .false.

      case default

        call write_error
        do ipr = 6,9,3
          write(ipr,110) mot
        end do
        stop

    end select

  end do boucle_lect

! Test on input file names
  do ifich = 1,nfich
    open(2, file = fichin(ifich), status='old', iostat=istat)
    if( istat /= 0 ) then
      mot = fichin(ifich)
      l = len_trim(mot)
      if( mot(l-3:l) /= '.txt' ) then
        mot(l+1:l+4) = '.txt'
        Close(2) 
        open(2, file = mot, status='old', iostat=istat)
        if( istat /= 0 ) then
          l = len_trim(mot)
          mot(l-3:l+2) = '_1.txt'
          open(2, file = mot, status='old', iostat=istat)
          if( istat == 0 ) then
           call write_error
           do ipr = 6,9,3
             write(ipr,112) fichin(ifich), mot
           end do
           stop
         else
           call write_open_error(fichin(ifich),istat,1)
         endif
        endif
        fichin(ifich) = mot
      else
       l = len_trim(mot)
       mot(l-3:l+2) = '_1.txt'
       open(2, file = mot, status='old', iostat=istat)
       if( istat == 0 ) then
         call write_error
         do ipr = 6,9,3
           write(ipr,112) fichin(ifich), mot
         end do
         stop
       else
         call write_open_error(fichin(ifich),istat,1)
       endif
      endif
      Close(2)
    endif
  end do

  allocate( fichscanout_all(nfich) ) 
  allocate( Convolution_out_all(nfich) )
  Convolution_out_all(:) = ' ' 
  fichscanout_all(:) = ' '

  call Conv_out_name(bav_open,check_conv,Chem,Chemin,convolution_out,Convolution_out_all,Dafs_bio,Deuxieme,fichin, &
                        fichscanin,fichscanout,fichscanout_all,Length_line,nfich,nomfichbav,Photoemission,Scan_a,Scan_true)

  do ifich = 1,nfich
    if( convolution_out /= fichin(ifich) ) cycle
    call write_error
    do ipr = 6,9,3
      write(ipr,120) fichin(ifich), mot
    end do
    stop
  end do

  if( .not. ( Seah .or. Arc ) .and. nelor == 0 ) then
    nelor = 1
    allocate( Elor(nelor) )
    allocate( betalor(nelor) )
    Elor(1) = 0._db
    betalor(1) = 0._db
  endif

! -- Table dimensions -------------------------------------

  ninitlm = 1
  do ifich = 1,nfich
    open(2, file = fichin(ifich), status='old', iostat=istat)
    n = nnombre(2,Length_line)
    if( n > 8 ) read(2,*) Eseuil(ifich), numat(ifich), nseuil, jseuil, fpp_avantseuil(ifich), V0muf(ifich), En_fermi(ifich), ninit
    if( n /= ninit + 11 .and. n /= ninit + 12 .and. n /= ninit + 13 .and. n /= ninit + 14 ) then
      call write_error
      do ipr = 6,9,3
        write(ipr,'(///,A/)') ' Old format in the first line of not convoluted file !'
      end do
      stop
    endif
    ninitlm = max( ninitlm, ninit )
    Close(2)
  end do

  call n_div_fpp_avanseuil(Eseuil,n_div_fpp,nfich)

  allocate( Abs_U_iso(nfich) )
  allocate( Eph1(nfich) )
  allocate( Ephm(nfich) )
  allocate( Epsii(ninitlm,nfich) )
  allocate( decal_initl(ninitlm,nfich) )
  allocate( natomsym_f(nfich) )
  allocate( nsup(ninitlm,nfich) )
  allocate( ne_initl(ninitlm,nfich) )
  allocate( Trunc(nfich) )

  Abs_U_iso(:) = 0._db
  Epsii(:,:) = 0._db
  decal_initl(:,:) = 0._db
  Trunc(:) = .false.

  call Dimension_file(Abs_in_bulk,Abs_U_iso,Cor_abs,Eintmax,En_Fermi,Eph1,Ephm,Epsii,f0_forward,Fichin,Fichscanin, &
          fprim,Full_self_abs,Green_int,jseuil,Length_line,Magn,mu_0_bulk,n_col_max,n_mat_pol,n_Trunc,natomsym_f,ne,nfich, &
          ninit1,ninitl,ninitlm,nphim,npldafs,nseuil,nxan,nxan_lib,Sample_thickness,Scan_true,Self_abs,Signal_Sph, &
          Surface_ref,Tenseur,Tenseur_car,Trunc,V0muf,Volume_maille,Volume_maille_bulk)

  if( .not. Cor_abs ) Double_cor = .false.

  if( jseuil >= 4 ) no_extrap = .true.

! Modification when fitting
  if( Fit_cal ) then

    do igr = 2,ngroup_par
      istop = 0
      do ipar = 1,npar(igr)
        if( typepar(igr,ipar) /= 'shift' .and. typepar(igr,ipar) /= 'weight' ) cycle
        if( indice_par(igr,ipar) > nfich ) then
          call write_error
          do ipr = 6,9,3
            write(ipr,145) typepar(igr,ipar), indice_par(igr,ipar), nfich
          end do
          istop = 1
         endif
      end do
    end do
    if( istop == 1 ) stop

    do igr = 2,ngroup_par
      do ipar = 1,npar(igr)
        select case( typepar(igr,ipar) )
          case('aseah')
            asea = param(igr,ipar)
          case('ecent')
            Ecent = param(igr,ipar)
            Ecent = Ecent / rydb
          case('gaussian')
            deltar = param(igr,ipar)
          case('vibr')
            vibration = param(igr,ipar)
          case('elarg')
            Elarg = param(igr,ipar)
            Elarg = Elarg / rydb
          case('gamma_max')
            Gamma_max = param(igr,ipar)
            Gamma_max = Gamma_max / rydb
          case('gamma_hol')
            Gamma_hole(1) = param(igr,ipar)
            Gamma_hole(1) = Gamma_hole(1) / rydb
            Gamma_hole_imp = .true.
          case('E_cut','EFermi')
            E_cut = param(igr,ipar)
            E_cut_man = .true.
          case('shift')
            if( .not. run_done( indice_par(igr,ipar) ) ) then
              ifich = i_done( indice_par(igr,ipar) )
              decal( ifich ) = param(igr,ipar)
            endif
          case('weight')
            if( .not. run_done( indice_par(igr,ipar) ) ) then
              ifich = i_done( indice_par(igr,ipar) )
              Pds( ifich ) = param(igr,ipar)
            endif
        end select
      end do
    end do
  endif

  deallocate( i_done )

  if( Cor_abs .or. Dafs_bio ) then
    npldafs_b = npldafs
  else
    npldafs_b = 0
  endif
  if( Analyzer ) then
    npldafs_t = npldafs
  else
    npldafs_t = npldafs / 2
  endif
  allocate( hkl_dafs(3,npldafs_b) )

  Dafs = npldafs > 0
   
  allocate( dph(npldafs) )
  allocate( dpht(npldafs) )
  dpht(:) = (0._db,0._db)
  allocate( fi(npldafs) ); allocate( fr(npldafs) )
  allocate( f0(npldafs) ) 
  f0(:) = (0._db, 0._db )
  allocate( nphi(npldafs) )
  if( Abs_in_bulk ) then
    allocate( f0_bulk(npldafs) )
    f0_bulk(:) = (0._db, 0._db )
    allocate( l_dafs(npldafs) )
    allocate( Length_rel(npldafs,n_Trunc) )
    allocate( Length_abs(npldafs) )
  endif

! Reading of Thomson part
  if( npldafs > 0 ) then

    if( Thomson ) then
      n = min( npldafs, npldafs_th )
      f0(1:n) = f0_th(1:n)
      deallocate( f0_th )
    else
      Pdt = 0._db; Pdt_bulk = 0._db
      do ifich = 1,nfich
        if( Skip_run(ifich) ) cycle
        open(2, file = fichin(ifich), status='old', iostat=istat)
        if( istat /= 0) call write_open_error(fichin(ifich),istat,1)
        read(2,*)
        n = nnombre(2,Length_line)
        if( n > 0 ) then
          read(2,*) ( fr(ipl), fi(ipl), ipl = 1,npldafs )
          if( Trunc(ifich) ) then
            f0_bulk(:) = f0_bulk(:) + Pds(ifich) * cmplx( fr(:), fi(:), db )
            Pdt_bulk = Pdt_bulk + Pds(ifich)
          else
            f0(:) = f0(:) + Pds(ifich) * cmplx( fr(:), fi(:), db )
            Pdt = Pdt + Pds(ifich)
          endif
        endif
        Close(2)
      end do
      if( abs(Pdt) > 1e-10_db ) f0(:) = f0(:) / Pdt
      if( Abs_in_bulk .and. abs(Pdt_bulk) > 1e-10_db ) f0_bulk(:) = f0_bulk(:) / Pdt_bulk
    endif

  endif

  if( E_cut_man ) En_Fermi(:) = E_cut

  if( Gaussian_default ) deltar = max( Eseuil(1) / 10000, 0.0001_db )

  deltar = deltar / sqrt( 8 * log(2._db) )
  vibration = vibration / sqrt( 8 * log(2._db) )

  decal(:) = decal(:) / rydb
  En_fermi(:) = En_fermi(:) / rydb
  Eph1(:) = Eph1(:) / rydb
  Ephm(:) = Ephm(:) / rydb
  Epsii(:,:) = Epsii(:,:) / rydb
  Eseuil(:) = Eseuil(:) / rydb
  V0muf(:) = V0muf(:) / rydb
  deltar = deltar / rydb
  Shift_U_iso = Shift_U_iso / rydb
  if( Fermip ) Efermip(1:nfich) = Efermip(1:nfich) / rydb
  if( .not. Arc .and. .not. seah ) Elor(:) = Elor(:) - En_fermi(1)

! Elaboration of energy grid

  call Shift_eval(decal,decal_initl,delta_edge,Energphot,Eph1,Ephm,Epsii,Epsii_ref,Epsii_ref_man,Eseuil,Esmin,nfich, &
                        ninit1,ninitl,ninitlm)
 
  pasdeb = 0.5_db / rydb
  nes_in = 100000
  allocate( Es_temp(nes_in) )
   
  call Output_Energy_Grid(decal_initl,Energphot,Eph1,Es_temp,Eseuil,Esmin,Estart,fichin,Length_line,ne,ne_initl,nes,nes_in, &
                              nfich,ninitl,ninitlm,nsup,pasdeb)  
  allocate( Es(nes) )
  Es(1:nes) = Es_temp(1:nes)
  deallocate( Es_temp )

  np_stokes = n_stokes * npldafs / 4
  allocate( angle(nphim) )
  allocate( As(nes,nphim,npldafs) )
  allocate( Icor(nes,nphim,npldafs) )
  allocate( Icirccor(nes,nphim,np_stokes) )
  allocate( Icircdcor(nes,nphim,np_stokes) )
  allocate( Icirc(nes,nphim,np_stokes) )
  allocate( Idcor(nes,nphim,npldafs) )
  allocate( indf(nes) )
  allocate( f0scan(nphim,npldafs) )
  allocate( p1f(nes) )
  allocate( p2f(nes) )
  allocate( phdtscan(nphim,npldafs) )
  allocate( Xs(nes,nxan) )
  allocate( Mus_mat(nes,n_stokes,n_mat_pol) )
  allocate( Mu_mat_comp(nes,4,n_mat_pol) )
  if( Abs_in_bulk ) then
    allocate( As_bulk(nes,nphim,npldafs,n_Trunc) )
    allocate( f0scan_bulk(nphim,npldafs) )
    allocate( Ts(nes) )
  endif
  if( Cor_abs ) allocate( mus(nes,nphim,npldafs,2) )

  Stokes = n_stokes > 0
  Stokes_Dafs = Circular .or. Full_self_abs
  Stokes_xan = Stokes .and. n_mat_pol > 0

  if( Circular ) then
! The analysor is supposed perfect
    i = n_Stokes - 3
! Sigma - Sigma
    Stokes_param(3,i) = 1._db; Stokes_param(4,i) = 0._db
    Stokes_param(5,i) = pi / 4
    i = i + 1
! Sigma - Pi
    Stokes_param(3,i) = 1._db; Stokes_param(4,i) = pi / 2
    Stokes_param(5,i) = pi / 4
    i = i + 1
    Stokes_param(5,i) = pi / 4
    Stokes_param(3,i) = -1._db; Stokes_param(4,i) = 0._db
    i = i + 1
    Stokes_param(3,i) = -1._db; Stokes_param(4,i) = pi / 2
    Stokes_param(5,i) = pi / 4
  endif

  if( Tenseur ) then
    n_col = 2*npldafs
  elseif( Bormann ) then
    n_col = nxan + n_mat_pol * n_stokes + 2*npldafs
  elseif( .not. Analyzer ) then
    n_col = nxan + n_mat_pol * ( n_stokes + 2 ) + npldafs / 2
  else
    n_col = nxan + n_mat_pol * ( n_stokes + 2 ) + npldafs
    if( Dafs_bio ) then
      if( icheck > 1 ) then
         n_col = n_col + 4*npldafs
         nw = 5
      else
         nw = 1
      endif
    else
      if( fprim ) n_col = n_col + 2*npldafs
      if( Self_abs ) n_col = n_col + 3*npldafs
      if( Full_self_abs ) n_col = n_col + 5*npldafs
      if( Double_cor ) n_col = n_col + npldafs
      if( Stokes_Dafs ) n_col = n_col + np_stokes
      if( Stokes .and. Full_self_abs ) n_col = n_col + np_stokes
      if( Stokes .and. Double_cor ) n_col = n_col + np_stokes
    endif
  endif

  allocate( nom_col(n_col) )

! Preparation is done --------------------------------------

  Sup_sufix = ninitl(1) > 1
   
  call Col_name(Analyzer,Bormann,Cor_abs,Dafs_bio,Double_cor,fichin(1),Fichscanin(1),fprim,Full_self_abs,hkl_dafs, &
      Length_line,n_col,n_mat_pol,n_stokes,nom_col,npldafs,npldafs_b,nxan,Photoemission,Self_abs,Signal_sph,Stokes, &
      Stokes_name,Stokes_param,Sup_sufix,Tenseur)

! This loop is for the case of output for all indata files  
  boucle_conv_file: do i_conv = 0,nfich
    if( i_conv > 0 ) then
      if ( Just_total .or. nfich == 1 ) exit
      Convolution_out = Convolution_out_all(i_conv)
      fichscanout = fichscanout_all(i_conv)
    endif

  Xs(:,:) = 0._db
  Mus_mat(:,:,:) = 0._db
  Mu_mat_comp(:,:,:) = (0._db,0._db)
  if( Dafs ) As(:,:,:) = (0._db,0._db)
  if( Cor_abs ) mus(:,:,:,:) = (0._db,0._db)
  if( Abs_in_bulk ) then
    As_bulk(:,:,:,:) = (0._db,0._db)
    Ts(:) = 0._db
    f0scan_bulk(:,:) = ( 0._db, 0._db )
  endif
  f0scan(:,:) = ( 0._db, 0._db )
  angle(:) = 0._db
  phdtscan(:,:) = ( 0._db, 0._db ) 
 
  ifich = 0
  i_Trunc = 0
  do jfich = 1,mfich
    if( Trunc(jfich) ) i_Trunc = i_Trunc + 1
    if( run_done(jfich) ) cycle
    ifich = ifich + 1

    if( .not. Just_total .and. i_conv > 0 .and. i_conv /= ifich ) cycle
     
    if( Fermip ) then
      E_cut_orig = Efermip(ifich)
    else
      E_cut_orig = En_fermi(ifich)
    endif

    do initl = 1,ninitl(ifich)

      if( n_selec_core /= 0 ) then
        do i = 1,n_selec_core
          if( initl == num_core(i) ) exit
        end do
        if( i > n_selec_core ) cycle
      endif

      nenerg = ne_initl(initl,ifich)

      allocate( Adafs(nenerg,nphim,npldafs) )
      Adafs(:,:,:) = (0._db, 0._db)
      if( Cor_abs ) then
        allocate( mu(nenerg,nphim,npldafs,2) )
        mu(:,:,:,:) = (0._db, 0._db)
      endif
      allocate( Xanes(nenerg,nxan) )
      allocate( Mu_mat(nenerg,6,n_mat_pol) )
      Xanes(:,:) = 0._db
      Mu_mat(:,:,:) = 0._db

! -- Lecture -----------------------------------------------------------

      open(2, file = fichin(ifich), status='old', iostat=istat)

      read(2,*)

      natomsym = natomsym_f(ifich)
      
      n = nnombre(2,Length_line)
      if( n > 0 ) then
        read(2,*)      ! This line has already been red
        read(2,*) ( fr(ipl), fi(ipl), ipl = 1,npldafs )
        dph(:) = cmplx( fr(:), fi(:),db )
        dpht(:) = dpht(:) + Pds(ifich) * dph(:)
        if( Cor_abs ) then
          read(2,*) axyz(:), angxyz(:), ( hkl_dafs(:,ipl), ipl = 1,npldafs)
          angxyz(:) = angxyz(:) * radian
          axyz(:) = axyz(:) / bohr
        elseif( Trunc(ifich) ) then
          read(2,*) Length_abs(:)
          read(2,*) Length_rel(:,i_Trunc)
          read(2,*) l_dafs(:)
        endif
      elseif( Tenseur ) then
        f0(:) = ( 0._db, 0._db )
        dph(:) = ( 0._db, 0._db )
        dph(1) = ( 1._db, 0._db )
        if( Tenseur_car ) then
          dph(4) = ( 1._db, 0._db )
          dph(6) = ( 1._db, 0._db )
        endif
        dpht(:) = dpht(:) + Pds(ifich) * dph(:)
      endif

      read(2,*)

      allocate( Energ(nenerg) )

      if( Dafs ) then

        nphi(:) = 1

        if( Cor_abs ) allocate( mua_r(npldafs,2) )
        if( Full_self_abs ) allocate( mua_i(npldafs,2) )

        do ie = nsup(initl,ifich)+1,nenerg
          if( Tenseur .and. .not. Magn ) then
            Read(2,*) Energ(ie), ((fr(ipl), ipl = 1,npldafs), j = 1,initl)
            if( .not. scan_true ) Adafs(ie,1,:) = cmplx(fr(:), 0._db,db)
          elseif( nxan > 0 .or. n_mat_pol > 0 ) then
            if( Full_self_abs ) then
              Read(2,*) Energ(ie),( ( Xanes(ie,ipl), ipl = 1,nxan ), ( Mu_mat(ie,:,ipl), ipl = 1,n_mat_pol ), &
                                   ( fr(ipl), fi(ipl), (mua_r(ipl,i), mua_i(ipl,i), i = 1,2), ipl = 1,npldafs ), j = 1,initl)
            elseif( self_abs ) then
              Read(2,*) Energ(ie),(( Xanes(ie,ipl), ipl = 1,nxan ), ( Mu_mat(ie,:,ipl), ipl = 1,n_mat_pol ), &
                                  ( fr(ipl), fi(ipl), (mua_r(ipl,i), i = 1,2), ipl = 1,npldafs ), j = 1,initl)
            else
              Read(2,*) Energ(ie), (( Xanes(ie,ipl), ipl = 1,nxan ), ( Mu_mat(ie,:,ipl), ipl = 1,n_mat_pol ), &
                                  ( fr(ipl), fi(ipl), ipl = 1,npldafs ), j = 1,initl )
            endif
          else
            if( Full_self_abs ) then
              Read(2,*) Energ(ie), ( ( fr(ipl), fi(ipl), ( mua_r(ipl,i), mua_i(ipl,i), &
                i = 1,2), ipl = 1,npldafs ), j = 1,initl)
            elseif( self_abs ) then
              Read(2,*) Energ(ie), ( ( fr(ipl), fi(ipl), ( mua_r(ipl,i), i = 1,2 ), ipl = 1,npldafs ), j = 1,initl)
            else
              Read(2,*) Energ(ie), ( ( fr(ipl), fi(ipl), ipl = 1,npldafs), j = 1,initl )
            endif
          endif
          if( .not. scan_true ) then
            Adafs(ie,1,:) = cmplx( fr(:), fi(:), db)
            if( Full_self_abs ) mu(ie,1,:,:) = cmplx( mua_r(:,:), mua_i(:,:), db)
            if( Self_abs ) mu(ie,1,:,:) = cmplx( mua_r(:,:), 0._db, db)
          endif
        end do

        if( Cor_abs ) deallocate( mua_r)
        if( Full_self_abs ) deallocate( mua_i )

      else

        do ie = nsup(initl,ifich)+1,nenerg
          Read(2,*) Energ(ie), (( Xanes(ie,ipl), ipl = 1,nxan ), ( Mu_mat(ie,:,ipl), ipl = 1,n_mat_pol ), j = 1,initl)
        end do

      endif

      close(2)

      if( Scan_true ) then
        open(2, file = fichscanin(ifich), status='old',iostat=istat)

        do ipl = 1,npldafs
          read(2,*) nphi(ipl)
        end do
        
        do ie = nsup(initl,ifich)+1,nenerg
          read(2,*) mots
          do ipl = 1,npldafs
            if( Dafs_bio ) then
              read(2,*) hkl_dafs(:,ipl)
            else
              read(2,*) mots
            endif
            do i = 1,nphi(ipl)
              if( Self_abs ) then
                read(2,*) angle(i), ( a, b, a3, a4, j = 1,initl ), ( c, d, c, d, j = initl+1,ninitl(ifich) ), &
                          a1, b1, a2, b2
                mu(ie,i,ipl,1) = cmplx( a3, 0._db, db )
                mu(ie,i,ipl,2) = cmplx( a4, 0._db, db )
              elseif( Full_self_abs ) then
                read(2,*) angle(i), ( a, b, a3, b3, a4, b4, j = 1,initl ), ( c, d, c, d, c, d, j = initl+1,ninitl(ifich) ), &
                          a1, b1, a2, b2
                mu(ie,i,ipl,1) = cmplx( a3, b3, db )
                mu(ie,i,ipl,2) = cmplx( a4, b4, db )
              else
                if( Dafs_bio ) then
                  read(2,*) ( a, b, j = 1,initl ), ( c, d, j = initl+1,ninitl(ifich) ), a1, b1, a2, b2
                else
                  read(2,*) angle(i), ( a, b, j = 1,initl ), ( c, d, j = initl+1,ninitl(ifich) ), a1, b1, a2, b2
                endif
              endif
              Adafs(ie,i,ipl) = cmplx( a, b, db)
              if( Trunc(ifich) ) then
                f0scan_bulk(i,ipl) = cmplx( a1, b1, db)
              else
                f0scan(i,ipl) = cmplx( a1, b1, db)
              endif
              phdtscan(i,ipl) = cmplx( a2, b2, db)
            end do
          end do
        end do
        close(2)

      endif

      Energ(:) = Energ(:) / rydb

      do ie = nsup(initl,ifich),1,-1
        Energ(ie) = Energ(ie+1) - pasdeb
      end do

      if( Decferm ) then
      ! Fermi does not follow the shift of the energy
        if( Energphot ) then
          E_cut = E_cut_orig - decal_initl(initl,ifich) + Eseuil(ifich) 
        else
          E_cut = E_cut_orig - decal_initl(initl,ifich) + Eseuil(ifich) - Esmin 
        endif
      elseif( Energphot ) then
        E_cut = E_cut_orig + Eseuil(ifich) 
      else
        E_cut = E_cut_orig 
      endif

      allocate( Mu_m(nenerg,4,n_mat_pol) )

      if( .not. ( Conv_done .or. Skip_run(ifich) ) ) then

        nenerge = nenerg + njp

        if( Cor_abs ) then
          allocate( Ad(nenerg,nphim,3*npldafs) )
        else
          allocate( Ad(nenerg,nphim,npldafs) )
        endif
        allocate( Xa(nenerge,nxan) )

        if( seah .or. Arc ) then
          nelor = nenerg
          allocate( Elor(nelor) )
          allocate( betalor(nelor) )
          Elor(:) = Energ(:)
        endif

        if( Photoemission .or. Green_int ) then
          betalor(:) = 0._db
        elseif( seah ) then
          call seahdench(asea,E_cut,Gamma_max,nelor,Elor,betalor)
        elseif( Arc ) then
          call gammarc(Ecent,Elarg,Gamma_max,E_cut,nelor,Elor,betalor)
        endif

        if( .not. Gamma_hole_imp ) then
          if( ninitl(ifich) == 1 .or. ninit1 == ninitl(ifich) .or. initl <= ninit1 ) then
            js = jseuil
          else
            js = jseuil + 1
          endif
          Gamma_h = Tab_width(Eseuil(ifich),js,nseuil,numat(ifich))
        elseif( ngamh == 1 .or. ninitl(ifich) == 1 ) then
          Gamma_h = Gamma_hole(1)
        elseif( ngamh == ninitl(ifich) ) then
          Gamma_h = Gamma_hole(initl)
        elseif( initl <= ninit1 ) then
          Gamma_h = Gamma_hole(1)
        else
          Gamma_h = Gamma_hole(2)
        endif
        betalor(:) = betalor(:) + Gamma_h

        if( Check_conv ) then
          ipr1 = 3
        else
          ipr1 = 6
        endif
        ipr2 = 6
        do ipr = ipr1,ipr2,3
          if( ifich == 1 .and. ( initl == 1 .or. initl == num_core(1) ) )  then
            if( seah ) then
              Write(ipr,'(/A)') ' Seah-Dench model'
              Write(ipr,150) Gamma_max*rydb, asea
            elseif( Arc ) then
              Write(ipr,'(/A)') ' Arctangent model'
              Write(ipr,160) Gamma_max*rydb, Ecent*rydb, Elarg*rydb
            endif
          endif
          if( nfich == 1 .and. ninitl(ifich) == 1 ) then
            Write(ipr,170) Gamma_h*rydb, E_cut_orig*rydb, decal_initl(initl,ifich)*rydb
          elseif( nfich == 1 ) then
            Write(ipr,172) Gamma_h*rydb, E_cut_orig*rydb, decal_initl(initl,ifich)*rydb, initl
          elseif( ninitl(ifich) == 1 ) then
            Write(ipr,174) Gamma_h*rydb, E_cut_orig*rydb, decal_initl(initl,ifich)*rydb, ifich
          else
            Write(ipr,176) Gamma_h*rydb, E_cut_orig*rydb, decal_initl(initl,ifich)*rydb, ifich, initl
          endif
        end do

        if( i_conv == 0 .and. ( ifich == 1 .or. Esmin /= Eseuil(ifich) ) &
                .and. ( initl == 1 .or. initl == num_core(1) ) .and. Gamma_max > eps10 ) then
          de_obj = ( Elor(nelor) - Elor(1) ) / 15
          E_obj = Elor(1) - 0.00001_db
          do ie = 1,nelor
            Gamm = betalor(ie) - Gamma_h
            if( Gamm > 0._db ) exit
          end do
          n = ie - 1
          do ie = 1,nelor
            if( ie > 1 .and. ie < n ) cycle
            if( Energphot ) then
              E = Elor(ie) - Eseuil(ifich) - V0muf(ifich)
            else
              E = Elor(ie) - V0muf(ifich)
            endif
            Gamm = betalor(ie) - Gamma_h
            if( Gamm > 0._db ) then
              alambda = sqrt( 2 / ( sqrt( E**2 + Gamm**2 ) - E ) )
            else
              alambda = 0._db
            endif
            E = E + V0muf(ifich)
            do ipr = ipr1,ipr2,3
              if( ie == 1 ) then
                if( Energphot ) then
                  Write(ipr,'(A)') '    E - E_edge   Width_(eV) lambda_(A)'
                else
                  Write(ipr,'(A)') '      E_(eV)     Width_(eV) lambda_(A)'
                endif
              endif
              if( ipr == 6 ) then
                if( Elor(ie) < E_obj .and. ie /= nelor ) cycle
                E_obj = E_obj + de_obj
              endif
              if( alambda * bohr < 10000000 ) then
                Write(ipr,180) E * Rydb, betalor(ie) * Rydb, alambda * bohr
              else
                Write(ipr,185) E * Rydb, betalor(ie) * Rydb, alambda * bohr
              endif
            end do
          end do
        endif
  
        allocate( dampl(nenerg) )
        allocate( bb(nenerge) )
        allocate( Energe(nenerge) )
        allocate( Ephoton(nenerg) )
        allocate( e1(nenerge) )
        allocate( e2(nenerge) )
        Ephoton(:) = 0._db
        allocate( lori(nenerge) )
        allocate( lorix(nenerge) )
        allocate( lorrx(nenerge) )
        if( Dafs .or. Stokes_xan ) allocate( lorr(nenerge) )

! Calculation of the coefficients of the lorentzian
        call cflor(bb,betalor,E_cut,Energe,Elor,Energ,ie1,ie2,nef,nelor,nenerg,nenerge,Photoemission)
        if( Energphot ) then
          Ephoton(:) = Energ(:) + decal_initl(initl,ifich)
        else
          Ephoton(:) = Energ(:) + decal_initl(initl,ifich) + Esmin
        endif
  
  ! Temperature effect in the Debye model
        if( ( U_iso_inp .and. Abs_U_iso_inp > eps10 ) .or. Abs_U_iso(ifich) > eps10   ) then

! natomsym is real
          fac = natomsym / ninitl(ifich)

          call Debye_effect(Abs_U_iso,Abs_U_iso_inp,Energ,Energphot,Ephoton,Eseuil,fac,ifich,n_col,nenerg,nfich, &
                        nom_col,numat,nxan,Shift_U_iso,U_iso_inp,V0muf,Xanes)

          fac = 1 / ninitl(ifich)

          if( Dafs ) call Debye_effect_a(Abs_U_iso,Abs_U_iso_inp,Adafs,dph,Energ,Energphot,Ephoton,Eseuil,fac,ifich, &
                        nenerg,nfich,.true.,nphim,npldafs,numat,Shift_U_iso,U_iso_inp,V0muf) 

          fac = natomsym * 100 / ( ninitl(ifich) * Volume_maille )
  
          if( Cor_abs ) then
            allocate( Mu_t(nenerg,nphim,npldafs) )
            allocate( dph_t(npldafs) )
            dph_t(:) = ( 1._db, 0._db )
            do i = 1,2
              Mu_t(:,:,:) = Mu(:,:,:,i) 
              call Debye_effect_a(Abs_U_iso,Abs_U_iso_inp,Mu_t,dph_t,Energ,Energphot,Ephoton,Eseuil,fac,ifich, &
                        nenerg,nfich,.false.,nphim,npldafs,numat,Shift_U_iso,U_iso_inp,V0muf)
              Mu(:,:,:,i) = Mu_t(:,:,:) 
            end do
            deallocate( dph_t, Mu_t )
          endif 

          if( Stokes_xan ) then
            allocate( Mu_tt(nenerg,npldafs) )
            do i = 1,6,5
              Mu_tt(:,:) = Mu_mat(:,i,:) 
              call Debye_effect(Abs_U_iso,Abs_U_iso_inp,Energ,Energphot,Ephoton,Eseuil,fac,ifich,n_col,nenerg,nfich, &
                        nom_col,numat,npldafs,Shift_U_iso,U_iso_inp,V0muf,Mu_tt)
              Mu_mat(:,i,:) = Mu_tt(:,:)  
            end do  
            deallocate( Mu_tt )
          endif

        endif

        if( no_extrap ) then
          Extrap = .false.
        else
          Extrap = Stokes_xan .or. Cor_abs
          if( Dafs .and. .not. Extrap ) then
            do ipl = 1,npldafs
              if( abs( dph(ipl) ) < 1.e-10_db ) cycle
              Extrap = .true.
              exit
            end do
          endif
        endif

        fprime_atom = icheck > 2 .and. ifich == 1 .and. i_conv == 0 .and. ( initl == 1 .or. initl == num_core(1) )
        if( Extrap .or. fprime_atom .or. Cor_abs .or. Stokes_xan ) then
          call Extrapat(bb(nenerg),dampl,Ephoton,Eseuil(ifich),Extrap,fpp0,fprime_atom,icheck,nenerg,numat(ifich))
           dampl(:) = dampl(:) / ninitl(ifich)
           fpp0 = fpp0 / ninitl(ifich)
        endif

! Convolution par la lorentzienne
        if( Dafs ) Ad(1:nenerg,:,1:npldafs) = Adafs(1:nenerg,:,1:npldafs)
        if( Cor_abs ) then
          do i = 1,2
            do ipl = 1,npldafs
              Ad(1:nenerg,:,i*npldafs+ipl) = mu(1:nenerg,:,ipl,i)
            end do
          end do
        endif

        do ipl = 1,nxan
          Xa(1:nenerg,ipl) = Xanes(1:nenerg,ipl)
! Extrapolation
          do ie = nenerg+1,nenerge
            Xa(ie,ipl) = Xa(nenerg,ipl)
          end do
        end do

        gamma = .false.
        do ie = nef,nenerge
          if( abs( bb(ie) ) > 1.e-10_db ) then
            gamma = .true.
            exit
          endif
        end do

        ne2 = ie2

        do ie = ie1,ie2
          if( ie == ie1 ) then
            if( Photoemission ) then
              e1(ie) = 1.5 * Energe(1) - 0.5 * Energe(2)
            else
              e1(ie) = E_cut
              if( ie == 1 ) then
                e1m = 1.5 * Energe(1) - 0.5 * Energe(2)
                e1(ie) = max( e1(ie), e1m )
              endif
            endif
          else
            e1(ie) = 0.5 * ( Energe(ie) + Energe(ie-1) )
          endif
          if( ie == ne2 ) then
            e2(ie) = 1.5 * Energe(ne2) - 0.5 * Energe(ne2-1)
          else
            e2(ie) = 0.5 * ( Energe(ie) + Energe(ie+1) )
          endif
        end do

        if( Photoemission ) then
          nen2 = ie2
        else
          nen2 = nenerg
        endif
   
        if( Gamma .and. .not. Green_int ) then
   
          do ie = 1,nenerg
   
            do j = ie1,ne2
   
              if( .not. Gamma_var ) then
                bba = bb(ie)
                if( bba >= 0._db ) then
                  bbb = max( bba, 1.E-08_db )
                else
                  bbb = bba
                endif
                de2 = ( e2(j) - Energe(ie) ) / bbb
                de1 = ( e1(j) - Energe(ie) ) / bbb
                lorix(j) = atan( de1 ) - atan( de2 )
                if( ( Dafs .or. Stokes_xan ) .and. j <= nenerg ) lorrx(j) = 0.5 * log( (1 + de1**2) / (1 + de2**2) )
              endif
              if( Dafs .or. Gamma_var ) then
                bba = bb(j)
                if( bba >= 0._db ) then
                  bbb = max( bba, 1.E-08_db )
                else
                  bbb = bba
                endif
                de2 = ( e2(j) - Energe(ie) ) / bbb
                de1 = ( e1(j) - Energe(ie) ) / bbb
                lori(j) = atan( de1 ) - atan( de2 )
                if( Gamma_var ) then
                  lorix(j) = lori(j)
                  lorrx(j) = lorr(j)
                endif
                if( ( Dafs .or. Stokes_xan ) .and. j <= nenerg ) lorr(j) = 0.5 * log( (1 + de1**2) / (1 + de2**2) )
              endif
            end do

            do ipl = 1,nxan
              Xanes(ie,ipl) = - sum( lorix(ie1:ne2) * Xa(ie1:ne2,ipl) ) / pi
            end do
     
            if( Stokes_xan ) then
              do ipl = 1,n_mat_pol
                Mu_m(ie,1,ipl) = sum( cmplx( lorrx(ie1:nen2), lorix(ie1:nen2), db) * Mu_mat(ie1:nen2,1,ipl) ) / pi
                Mu_m(ie,2,ipl) = sum( cmplx( lorrx(ie1:nen2), lorix(ie1:nen2), db) &
                                    * cmplx( Mu_mat(ie1:nen2,2,ipl), Mu_mat(ie1:nen2,3,ipl), db ) ) / pi
                Mu_m(ie,3,ipl) = sum( cmplx( lorrx(ie1:nen2), lorix(ie1:nen2), db) &
                                    * cmplx( Mu_mat(ie1:nen2,4,ipl), Mu_mat(ie1:nen2,5,ipl), db ) ) / pi
                Mu_m(ie,4,ipl) = sum( cmplx( lorrx(ie1:nen2), lorix(ie1:nen2), db) * Mu_mat(ie1:nen2,6,ipl) ) / pi
              end do
            endif
     
            do ipl = 1,npldafs
              do ip = 1,nphi(ipl)

                if( Tenseur ) then

 ! Here Ad is in number of electron, the convolution needs division by pi
                  Adafs(ie,ip,ipl) = sum( cmplx( lorr(ie1:nen2), lori(ie1:nen2),db ) * real( Ad(ie1:nen2,ip,ipl),db ) ) / pi

                else

                  Adafs(ie,ip,ipl) = sum( cmplx(lorr(ie1:nen2),lori(ie1:nen2),db) * Ad(ie1:nen2,ip,ipl) )

! Here the imaginary part becomes the Kramers Koenig of the real part
                  if( Cor_abs ) then
                    do i = 1,2
                      jpl = ipl + i*npldafs
                      mu(ie,ip,ipl,i) = sum( cmplx( lorrx(ie1:nen2), lorix(ie1:nen2), db) * Ad(ie1:nen2,ip,jpl) ) / pi
                    end do
                  endif

                endif

              end do
            end do

          end do

        else

          if( Dafs .or. Stokes_xan ) then

            do ie = 1,nenerg

              do j = ie1,nen2
                de2 = e2(j) - Energe(ie)
                de1 = e1(j) - Energe(ie)
                if( j == ie ) then
                  if( abs( de2 + de1 ) < 1.e-10_db ) then
                    lorr(j) = 0._db
                  elseif( abs( de1 ) < 1.e-10_db ) then
                    lorr(j) = - log( 4._db )
                  else
                    lorr(j) = log( - de1 / de2 )
                  endif
                else
                  lorr(j) = log( de1 / de2 )
                endif
              end do

              do ipl = 1,npldafs
                do ip = 1,nphi(ipl)

                  if( Tenseur .and. .not. Green_int ) then
 ! Here Ad is in number of electron, the convolution needs division by pi
                    Adafs(ie,ip,ipl) = sum( lorr(ie1:nen2) * real( Ad(ie1:nen2,ip,ipl),db ) ) / pi
                    if( ie >= ie1 ) Adafs(ie,ip,ipl) = Adafs(ie,ip,ipl) - img * real( Ad(ie,ip,ipl),db )
                  else
                    if( .not. Green_int ) then
                      Adafs(ie,ip,ipl) = sum( lorr(ie1:nen2) * Ad(ie1:nen2,ip,ipl) )
                      if( ie >= ie1 ) Adafs(ie,ip,ipl) = Adafs(ie,ip,ipl) + img * pi * Ad(ie,ip,ipl)
                    endif
! For the correction by absorption, the second xanes becomes the Kramers Koenig of the frst one
                    if( Cor_abs ) then
                      do i = 1,2
                        jpl = ipl + i*npldafs
                        mu(ie,ip,ipl,i) = sum( lorr(ie1:nen2) * Ad(ie1:nen2,ip,jpl)) / pi
                        if( ie >= ie1 ) mu(ie,ip,ipl,i) = mu(ie,ip,ipl,i) - img * Ad(ie,ip,jpl)
                      end do
                    endif
                  endif
                end do
              end do

              if( Stokes_xan ) then
                do ipl = 1,n_mat_pol
                  Mu_m(ie,1,ipl) = sum( lorr(ie1:nen2) * Mu_mat(ie1:nen2,1,ipl) ) / pi
                  if( ie >= ie1 ) Mu_m(ie,1,ipl) = Mu_m(ie,1,ipl) - img * Mu_mat(ie,1,ipl)
                  Mu_m(ie,2,ipl) = sum( lorr(ie1:nen2) &
                             * cmplx( Mu_mat(ie1:nen2,2,ipl), Mu_mat(ie1:nen2,3,ipl), db ) ) / pi
                  if( ie >= ie1 ) Mu_m(ie,2,ipl) = Mu_m(ie,2,ipl) - img * cmplx( Mu_mat(ie,2,ipl), Mu_mat(ie,3,ipl), db )
                  Mu_m(ie,3,ipl) = sum( lorr(ie1:nen2) &
                             * cmplx( Mu_mat(ie1:nen2,4,ipl), Mu_mat(ie1:nen2,5,ipl), db ) ) / pi
                  if( ie >= ie1 ) Mu_m(ie,3,ipl) = Mu_m(ie,3,ipl) - img * cmplx( Mu_mat(ie,4,ipl), Mu_mat(ie,5,ipl), db )
                  Mu_m(ie,4,ipl) = sum( lorr(ie1:nen2) * Mu_mat(ie1:nen2,6,ipl) ) / pi
                  if( ie >= ie1 ) Mu_m(ie,4,ipl) = Mu_m(ie,4,ipl) - img * Mu_mat(ie,6,ipl)
                end do
              endif
          
            end do

          endif

          if( Photoemission ) then
            do ipl = 1,nxan
              Xanes(nef+1:nenerg,ipl) = 0._db
            end do
          elseif( .not. Green_int ) then
            do ipl = 1,nxan
              Xanes(1:nef-1,ipl) = 0._db
            end do
          endif
          if( Stokes_xan ) then
            do ipl = 1,n_mat_pol
              Mu_m(1:nef-1,1:4,ipl) = 0._db
            end do
          endif
     
        endif

! Transformation to crystallo convention with f" > 0
        if( Dafs ) then
          Adafs(:,:,:) = conjg( Adafs(:,:,:) )
          if( Cor_abs ) mu(:,:,:,:) = Conjg( mu(:,:,:,:) )
        endif
        if( Stokes_xan ) Mu_m(:,:,:) = Conjg( Mu_m(:,:,:) )

! One wants f' and f" in number of electrons
        if( Tenseur_car ) then
          do ie = 1,nenerg
            Eph = Ephoton(ie)
            ct_nelec = conv_mbarn_nelec(Eph) / pi
            adafs(ie,:,:) = ct_nelec * adafs(ie,:,:)
          end do
        endif

        if( Extrap ) then
          do ipl = 1,npldafs
            if( Signal_Sph .and. ipl > 2 ) exit
            do ip = 1,nphi(ipl)
              if( ip == 1 ) then
                cf = dph(ipl)
              else
                cf = phdtscan(ip,ipl)
              endif
              if( Green_int ) then
                do ie = 1,nenerg
                  Adafs(ie,ip,ipl) = Adafs(ie,ip,ipl) + img * fpp0 * cf
                end do
              else
                do ie = 1,nenerg
                  Adafs(ie,ip,ipl) = Adafs(ie,ip,ipl) + ( dampl(ie) + img * fpp0 ) * cf
                end do
              endif
            end do
          end do
        endif

! One transforms to crystallo convention with f" > 0
        if( Bormann ) Adafs(:,:,:) = conjg( Adafs(:,:,:) )

        if( Cor_abs ) then
          do ipl = 1,npldafs
            if( Full_self_abs .and. ( mod(ipl,4) == 2 .or. mod(ipl,4) == 3 ) ) cycle
            do ie = 1,nenerg ! Unit cell volume is in A^3
              fac = natomsym * 100 / ( Volume_maille * Conv_mbarn_nelec(Ephoton(ie)) )
              mu(ie,:,ipl,:) = mu(ie,:,ipl,:) + fac * dampl(ie)
            end do
          end do
        endif
      
        if( Stokes_xan ) then
          do ie = 1,nenerg
            fac = natomsym * 100 / ( Volume_maille * Conv_mbarn_nelec(Ephoton(ie)) )
            Mu_m(ie,1,:) = Mu_m(ie,1,:) + fac * dampl(ie)
            Mu_m(ie,4,:) = Mu_m(ie,4,:) + fac * dampl(ie)
          end do
! Real part is absorption
        endif

! Real part is absorption
        if( Cor_abs ) mu(:,:,:,:) = img * Conjg( mu(:,:,:,:) )
        if( Stokes_xan ) Mu_m(:,:,:) = img * Conjg( Mu_m(:,:,:) )
        
! Addition of the constant term corresponding to the edges of lower energy (but on dichroic terms)
        if( i_conv == 0 .and. Abs_before .and. initl == 1 ) then
          do i = 1,nxan
            nomab = nom_col(i)
            nomab = adjustl(nomab)
            if( nomab(1:3) == 'dic' .or. nomab(1:3) == 'Dic' ) cycle
            Xanes(:,i) = Xanes(:,i) + fpp_avantseuil(ifich) / n_div_fpp(ifich)
          end do
        endif
    
        if( ncal > 1 ) then
          iscr = 100 + jfich
          rewind( iscr )
          do ie = 1,nenerg
            write(iscr,*) Xanes(ie,1:nxan)
          end do
          if( Dafs ) then
            do ie = 1,nenerg
              write(iscr,*) ( ( Adafs(ie,ip,ipl), ip = 1,nphi(ipl) ), ipl = 1,npldafs )
            end do
          endif
        endif

        deallocate( Ad )
        deallocate( bb )
        deallocate( dampl )
        deallocate( Energe )
        deallocate( Ephoton )
        deallocate( e1 )
        deallocate( e2 )
        deallocate( lori )
        deallocate( lorix, lorrx )
        if( Dafs .or. Stokes_xan ) deallocate( lorr )
        deallocate( Xa )
        if( Seah .or. Arc ) then
          deallocate( Elor )
          deallocate( betalor )
        endif
 
      else

        iscr = 100 + jfich
        rewind( iscr )
        do ie = 1,nenerg
          read(iscr,*) Xanes(ie,1:nxan)
        end do
        if( Dafs ) then
          do ie = 1,nenerg
            read(iscr,*) ( ( Adafs(ie,ip,ipl), ip = 1,nphi(ipl) ), ipl = 1,npldafs )
          end do
        endif

      endif
      
      Energ(1:nenerg) = Energ(1:nenerg) + decal_initl(initl,ifich)
      
      fac = 1._db ! (Rydb)
      
      do ie = 1,nes
        E = Es(ie)
        E = max(E ,Energ(1) )
        E = min(E ,Energ(nenerg) )
        do i = 2,nenerg
          if( Energ(i) > E - eps10 ) exit
        end do
        p1 = ( E - Energ(i-1) ) / ( Energ(i) - Energ(i-1) )
        p2 = 1 - p1
        if( Es(ie) < Energ(1) - fac .or. Es(ie) > Energ(nenerg) + fac &
                                    .or. ( Es(ie) > Energ(nenerg) + eps10 .and. Es(nes) - Es(ie) > fac ) ) then
          p1f(ie) = 0._db
          p2f(ie) = 0._db
        else  
          p1f(ie) = Pds(ifich) * p1
          p2f(ie) = Pds(ifich) * p2
        endif
        indf(ie) = i
      end do

      do ie = 1,nes
        i = indf(ie)
        if( Dafs .and. .not. Trunc(ifich) ) As(ie,:,:) = As(ie,:,:) + p2f(ie) * Adafs(i-1,:,:) + p1f(ie) * Adafs(i,:,:)
        if( Cor_abs ) mus(ie,:,:,:) = mus(ie,:,:,:) + p2f(ie) * mu(i-1,:,:,:) + p1f(ie) * mu(i,:,:,:)

        if( Stokes_xan ) Mu_mat_comp(ie,:,:) = Mu_mat_comp(ie,:,:) + p2f(ie) * Mu_m(i-1,:,:) + p1f(ie) * Mu_m(i,:,:)
        if( Trunc(ifich) ) then
          Ts(ie) = Ts(ie) + p2f(ie) * Xanes(i-1,nxan) + p1f(ie) * Xanes(i,nxan)
          As_bulk(ie,:,:,i_Trunc) = As_bulk(ie,:,:,i_Trunc) + p2f(ie) * Adafs(i-1,:,:) + p1f(ie) * Adafs(i,:,:)
        else
          Xs(ie,:) = Xs(ie,:) + p2f(ie) * Xanes(i-1,:) + p1f(ie) * Xanes(i,:)
        endif
      end do

      deallocate( Adafs )
      deallocate( Energ )
      deallocate( Mu_mat )
      deallocate( Mu_m )
      deallocate( Xanes )
      if( Cor_abs ) deallocate( mu )
      
    end do ! end of loop over initial states (initl)

  end do ! end of loop over files

  if( .not. ( seah .or. Arc ) ) then
    deallocate( Elor )
    deallocate( betalor )
  endif

! One applies the truncation
  if( Abs_in_bulk ) then
    allocate( Trs(nes,npldafs) )
! Addition of pre-edge absorption and conversion in micrometer^-1
    c_micro = 100 / Volume_maille_bulk
    Ts(:) = c_micro * ( Ts(:) + mu_0_bulk )
    
    do ipl = 1,npldafs 
      Trs(:,ipl) = 1 / ( 1 - exp( - Ts(:) * Length_abs(ipl) - img * 2 * pi * l_dafs(ipl) ) )
      do i_Trunc = 1,n_Trunc 
        do ip = 1,nphi(ipl)
          As(:,ip,ipl) = As(:,ip,ipl) + exp( - Ts(:) * Length_rel(ipl,i_Trunc) ) * Trs(:,ipl) * As_bulk(:,ip,ipl,i_Trunc) 
        end do
      end do
    end do
    
    if( icheck > 1 ) then
      write(3,'(/a17,1p,e13.5,a14)') ' mu_0_bulk      =', c_micro*mu_0_bulk,' micrometer^-1'
      write(3,'(A/,23x,1p,120e13.5)') ' Length_abs(1,11,...) (micrometer) =', ( Length_abs(i), i = 1,min(npldafs,201), 10 )
      do i_Trunc = 1,n_Trunc
        write(3,'(A/,23x,1p,120e13.5)') ' Length_rel(1,11,...) (micrometer) =', &
                                                                      ( Length_rel(i,i_Trunc), i = 1,min(npldafs,201), 10 )
      end do
      write(3,'(/A)') '    Energy      Ts       1 - exp(-Ts*Length_abs(i)), i = 1,11,...'
      do ie = 1,nes
        write(3,'(f10.3,1p,121e13.5)') Es(ie)*rydb, Ts(ie), &
                                   ( 1 - exp( - Ts(ie) * Length_abs(i)), i = 1,min(npldafs,201), 10 ) 
      end do
    endif
    
  endif
  
  if( abs( S0_2 - 1._db ) > eps10 ) then
    Xs(:,:) = S0_2 * Xs(:,:)
    if( Dafs ) As(:,:,:) = S0_2 * As(:,:,:)
  endif

  if( .not. Forbidden .and. Dafs .and. .not. Tenseur ) then
    do ipl = 1,npldafs
      do ip = 1,nphi(ipl)
        if( nphi(ipl) == 1 ) then
          if( Bormann ) then
            As(:,ip,ipl) = As(:,ip,ipl) + conjg( f0(ipl) )
          else
            As(:,ip,ipl) = As(:,ip,ipl) + f0(ipl)
            if( Abs_in_bulk ) As(:,ip,ipl) = As(:,ip,ipl) + Trs(:,ipl) * f0_bulk(ipl)
          endif
        else
          if( Bormann ) then
            As(:,ip,ipl) = As(:,ip,ipl) + conjg( f0scan(ip,ipl) )
          else
            As(:,ip,ipl) = As(:,ip,ipl) + f0scan(ip,ipl)
            if( Abs_in_bulk ) As(:,ip,ipl) = As(:,ip,ipl) + Trs(:,ipl) * f0scan_bulk(ip,ipl)
          endif
        endif
      end do
    end do
  endif

  if( Abs_in_bulk ) deallocate ( Trs )
  
! In fact one calculates epsilon - 1 = -4*pi*Atom_density*r0/k2
! We are in atomic units r0 --> r0 / a0 = alfa^2
! k = 1/2 alfa E
! Volume_maille = unit cell volume in A^3
! ct = -4*pi * alfa^2  * 4  * bohr**3 / ( alfa^2 * E^2 * Volume_maille )
!    = -16*pi  * bohr**3 / ( E^2 * Volume_maille )
  if( Bormann ) then
    do ie = 1,nes
      if( Energphot ) then
        Eph = Es(ie) + Esmin
      else
        Eph = Es(ie)
      endif
      ct_epsilon = - 16 * pi * bohr**3 / ( Volume_maille * Eph**2 )
      As(ie,:,:) = ct_epsilon * As(ie,:,:)
    end do
  endif

  if( Cor_abs ) then
! Conversion from Megabarn to micrometer^-1
    c_micro = 100 / Volume_maille
! mus is already in micrometre -1
    do ipl = 1,npldafs
      if( Full_self_abs .and. ( mod(ipl,4) == 2 .or. mod(ipl,4) == 3 ) ) cycle
      mus(:,:,ipl,:) = mus(:,:,ipl,:) + c_micro * cmplx( fpp_avantseuil(ifich), f0_forward, db )
    end do
  endif

  if( Stokes_xan ) then
! Conversion from Megabarn to micrometer^-1
    c_micro = 100 / Volume_maille
! Mu_mat_comp is already in micrometre -1
    do i = 1,4,3
      Mu_mat_comp(:,i,:) = Mu_mat_comp(:,i,:) + c_micro * cmplx( fpp_avantseuil(ifich), f0_forward, db)
    end do
  endif

! Convolution by a Gaussian
! la vibration ne fonctionne pas !
  if( abs(deltar) > 1.e-10_db .or. abs(vibration) > 1.e-10_db ) then

    allocate( Yr(nes) )
    allocate( Yi(nes) )

    do ipl = 1,nxan
      Yr(:) = Xs(:,ipl)
      call gaussi(deltar,E_cut,Es,nes,vibration,Yr)
      Xs(:,ipl) = Yr(:)
    end do

    do ipl = 1,npldafs
      do ip = 1,nphi(ipl)
        Yr(:) = real( As(:,ip,ipl), db )
        call gaussi(deltar,E_cut,Es,nes,vibration,Yr)
        Yi(:) = aimag( As(:,ip,ipl) )
        call gaussi(deltar,E_cut,Es,nes,vibration,Yi)
        As(:,ip,ipl) = cmplx( Yr(:), Yi(:),db )

        if( Cor_abs ) then
          do i = 1,2
            Yr(:) = real( mus(:,ip,ipl,i), db )
            call gaussi(deltar,E_cut,Es,nes,vibration,Yr)
            Yi(:) = aimag( mus(:,ip,ipl,i) )
            call gaussi(deltar,E_cut,Es,nes,vibration,Yi)
            mus(:,ip,ipl,i) = cmplx( Yr(:), Yi(:),db )
          end do
        endif

      end do
    end do

    deallocate( Yr )
    deallocate( Yi )

  endif

  if( Stokes .and. Full_Self_abs ) call Cal_stokes(As,Icirc,nes,n_stokes,np_stokes,nphi,nphim,npldafs,Stokes_param)

  if( Stokes_xan ) then
    if( icheck > 1 ) then
      write(3,190) (( ipl, i = 1,8 ), ipl = 1,n_mat_pol)  
      do ie = 1,nes
        write(3,'(f10.3,1p,121e13.5)') Es(ie)*rydb, ( Mu_mat_comp(ie,:,ipl), ipl = 1,n_mat_pol ) 
      end do
    endif
    call Cal_MuTrans(icheck,Mu_mat_comp,Mus_mat,n_stokes,n_mat_pol,nes,Sample_thickness,Stokes_param)
  endif

  if( Cor_abs ) &
    call Corr_abs(angxyz,As,axyz,c_micro,d_dead,Double_cor,Eseuil(1),fpp_avantseuil(1),hkl_dafs,hkl_S,Icor,Icirccor, &
          Icircdcor,Idcor,mus,n_stokes,nes,np_stokes,nphi,nphim,npldafs,Self_abs,Stokes_param)

!---- Writing -------------------------------------------------

  allocate( Tens(n_col) )

  if( fprim ) then
    do ipl = 1,npldafs
      if( abs( dpht(ipl) ) > 1.e-10_db ) cycle
      dpht(ipl) = (1._db,0._db)
    end do
  endif

  zero_c = (0._db, 0._db)

  if( Double_cor ) then
    n_signal = 3 * nes
    i1 = n_col - npldafs_t - 4*npldafs + 1
    ipas = 5
  elseif( Cor_abs ) then
    n_signal = 2*nes
    i1 = n_col - npldafs_t - 3*npldafs + 1
    ipas = 4
  else
    n_signal = nes
    i1 = n_col - npldafs_t + 1
    ipas = 1
  endif
  if( Transpose_file ) allocate( Signal(npldafs_t,n_signal) )

  if( .not. Analyzer ) then

    do j = 1,2
    
      i_hk = 0
      index_hk = 0
      mot6_b = ' '
      do i_col = i1,n_col,ipas
        mot15 = ' '
        mot15 = nom_col(i_col)
        do i = 1,15
          if( mot15(i:i) /= '(' ) cycle
          mot6 = ' '
          if( mot15(i+1:i+1) == '-' .and. mot15(i+3:i+3) == '-' ) then
            mot6(1:4) = mot15(i+1:i+4)
          elseif( mot15(i+1:i+1) == '-' .or. mot15(i+2:i+2) == '-' ) then
            mot6(1:3) = mot15(i+1:i+3)
          else
            mot6(1:2) = mot15(i+1:i+2)  
          endif
          exit
        end do

        do i = 1,14
          if( mot15(i:i) /= ')' ) cycle
          ll = len_trim( mot6 )
          if( mot15(i+1:i+1) /= ' ' ) mot6(ll+1:ll+1) = mot15(i+1:i+1)
          if( i >= 14 ) exit
          if( mot15(i+2:i+2) /= ' ' ) mot6(ll+2:ll+2) = mot15(i+2:i+2)
          exit
        end do

        if( mot6 == mot6_b ) then    
          index_hk = index_hk + 1
          if( j == 2 ) n_index_hk(i_hk) = index_hk 
        else  
          i_hk = i_hk + 1
          index_hk = 1
          mot6_b = mot6 
        endif
      end do
      
      if( j == 1 ) then
        allocate( n_index_hk(i_hk) )
        n_index_hk(:) = 1
      endif
      
    end do
  
  endif
  
  First_E = .true.
  
  do ie = 1,nes

    do ipl = 1,nxan
      Tens(ipl) = Xs(ie,ipl)
    end do

    jpl = nxan

    do ipl = 1,n_mat_pol
      jpl = jpl + 1
      Tens(jpl) = real( Mu_mat_comp(ie,1,ipl), db )   ! sigma-sigma
      jpl = jpl + 1
      Tens(jpl) = real( Mu_mat_comp(ie,4,ipl), db )   ! pi-pi
      do i = 1,n_stokes
        jpl = jpl + 1
        Tens(jpl) = Mus_mat(ie,i,ipl)
      end do
    end do

    if( Dafs_bio ) then
! There are 4 polarization states corresponding to Q and -Q, both with sigma out and pi out.
! The intensity is the difference between I(Q) = I(Q)_sigma + I(Q)_pi, and I(-Q)
      do ipl = 1,npldafs
        Ip_sig = real(As(ie,1,ipl),db)**2 + aimag(As(ie,1,ipl))**2
        Ip_pi = real(As(ie,2,ipl),db)**2 + aimag(As(ie,2,ipl))**2
        Im_sig = real(As(ie,3,ipl),db)**2 + aimag(As(ie,3,ipl))**2
        Im_pi = real(As(ie,4,ipl),db)**2 + aimag(As(ie,4,ipl))**2
        do i = 1,nw
          jpl = jpl + 1
          select case(i)
            case(1)
              select case( Dafs_exp_type )
                case(1)
                  Tens(jpl) = Ip_sig + Ip_pi
                case(2)
                  Tens(jpl) = sqrt( Ip_sig + Ip_pi )
                case(3)
                  Tens(jpl) = Ip_sig + Ip_pi - ( Im_sig + Im_pi )
                case default
                  Tens(jpl) = sqrt( Ip_sig + Ip_pi ) - sqrt( Im_sig + Im_pi )
              end select
            case(2)
              Tens(jpl) = Ip_sig
            case(3)
              Tens(jpl) = Ip_pi
            case(4)
              Tens(jpl) = Im_sig
            case(5)
              Tens(jpl) = Im_pi
          end select
        end do
      end do

    else

      i_hk = 1
      kpl = 0      
      do ipl = 1,npldafs
        if( .not. ( Tenseur .or. Bormann ) ) then
          if( Analyzer ) then
           jpl = jpl + 1
           if( Surface_ref > eps10 ) then
              Tens(jpl) = abs( As(ie,1,ipl) / Surface_ref )**2  ! 2D diffraction case
            else
              Tens(jpl) = abs( As(ie,1,ipl) )**2
            endif
            if( Transpose_file ) Signal(ipl,ie) = Tens(jpl)
          else
            if( ipl > 2 * sum( n_index_hk(1:i_hk) ) ) i_hk = i_hk + 1
            n = 2 * sum( n_index_hk(1:i_hk-1) )
            if( ipl > n .and. ipl <= n + n_index_hk(i_hk) ) then
              jpl = jpl + 1
              kpl = kpl + 1
              j = ipl + n_index_hk(i_hk) 
              if( Surface_ref > eps10 ) then
                Tens(jpl) = abs( As(ie,1,ipl) / Surface_ref )**2 + abs( As(ie,1,j) / Surface_ref )**2  ! 2D diffraction case
              else
                Tens(jpl) = abs( As(ie,1,ipl) )**2 + abs( As(ie,1,j) )**2
              endif
              if( Transpose_file ) Signal(kpl,ie) = Tens(jpl)
            endif
          endif
        endif
        if( fprim .or. Tenseur .or. Bormann ) then
          jpl = jpl + 1
          Tens(jpl) = real( As(ie,1,ipl), db )
          jpl = jpl + 1
          Tens(jpl) = aimag( As(ie,1,ipl) )
          if( Tenseur .or. Bormann ) cycle
        endif
        if( Cor_abs ) then
          jpl = jpl + 1
          Tens(jpl) = Icor(ie,1,ipl)
          if( Transpose_file ) Signal(ipl,nes+ie) = Tens(jpl)
        endif
        if( Double_cor ) then
          jpl = jpl + 1
          Tens(jpl) = Idcor(ie,1,ipl)
          if( Transpose_file ) Signal(ipl,2*nes+ie) = Tens(jpl)
        endif
        if( Cor_abs ) then
          do i = 1,2
            jpl = jpl + 1
            Tens(jpl) = real( mus(ie,1,ipl,i), db )
            if( Full_self_abs ) then
              jpl = jpl + 1
              Tens(jpl) = aimag( mus(ie,1,ipl,i) )
            endif
          end do
        endif
        if( Stokes .and. mod(ipl,4) == 0 ) then
          j = ( ipl/4 - 1 )*n_stokes
          do i = 1,n_stokes
            jpl = jpl + 1
            Tens(jpl) = Icirc(ie,1,j+i)
            if( Cor_abs ) then
              jpl = jpl + 1
              Tens(jpl) = Icirccor(ie,1,j+i)
            endif
            if( Double_cor ) then
              jpl = jpl + 1
              Tens(jpl) = Icircdcor(ie,1,j+i)
            endif
          end do
        endif
      end do

    endif

    if( Fit_cal ) then
! Command the writing in a temporary file
      n = -iscratchconv
    else
      n = 0
    endif
    allocate( Ep(1) )
    Ep(:) = 0._db
    call Write_out(Abs_U_iso_inp,rdum,rdum,zero_c,E_cut,Es(ie),Ep,0._db,First_E,.false.,rdum, &
                   1,0,n_col,jpl,0,1,Convolution_out,nom_col,1,0,0,0,0,n,cdum,cdum,Tens,V0muf(1),.false.,0._db, &
                   Surface_ref,Volume_maille)
    deallocate( Ep )
    First_E = .false.
  end do

  deallocate( Tens )
  
  if( Transpose_file ) then
    if( n_energ_tr == 0 ) then
      n_energ_tr = nes
      allocate( Energ_tr(n_energ_tr) )
      Energ_tr(:) = Es(:)
    else
      Energ_tr(:) = Energ_tr(:) / rydb
    endif
    call Write_transpose(Convolution_out,Energ_tr,Es,n_col,n_energ_tr,n_signal,nes,nom_col,npldafs_t,Signal)
    deallocate( Energ_tr, Signal )
  endif

  if( .not. Analyzer ) deallocate( n_index_hk )
  
  if( Scan_true .and. .not. Dafs_bio ) then

    Open(7, file = fichscanout)

! One skips the '_0'
    do i = 1,n_col
      nomab = adjustl( nom_col(i) )
      if( nomab(1:1) /= 'I' ) cycle
      l = len_trim(nomab)
      if( nomab(l-1:l) == '_0' ) nomab(l-1:l) = '  '
      nom_col(i) = nomab
    end do

    do ie = 1,nes
      Write(7,*)
      Write(7,200) Es(ie) * rydb
      jpl = nxan

      do ipl = 1,npldafs

        jpl = jpl + 1
        Write(7,*) nom_col(jpl)
        do ip = 1,nphi(ipl)
          Write(7,210) angle(ip), real(As(ie,ip,ipl),db)**2 + aimag(As(ie,ip,ipl))**2
        end do

        if( fprim ) jpl = jpl + 2

        if( Cor_abs ) then
          jpl = jpl + 1
          Write(7,*) nom_col(jpl)
          do ip = 1,nphi(ipl)
            Write(7,210) angle(ip), Icor(ie,ip,ipl)
          end do
          if( Double_cor  ) then
            jpl = jpl + 1
            Write(7,*) nom_col(jpl)
            do ip = 1,nphi(ipl)
              Write(7,210) angle(ip), Idcor(ie,ip,ipl)
            end do
          endif
          jpl = jpl + 2  ! pour sauter les mu
          if( Full_self_abs ) jpl = jpl + 2  ! pour sauter les mu
        endif

        if( Stokes .and. mod(ipl,4) == 0 ) then
          j = ( ipl/4 - 1 )*n_stokes
          do i = 1,n_stokes
            jpl = jpl + 1
            Write(7,*) nom_col(jpl)
            do ip = 1,nphi(ipl)
              Write(7,210) angle(ip), Icirc(ie,ip,j+i)
            end do
            jpl = jpl + 1
            Write(7,*) nom_col(jpl)
            do ip = 1,nphi(ipl)
              Write(7,210) angle(ip), Icirccor(ie,ip,j+i)
            end do
            if( Double_cor  ) then
              jpl = jpl + 1
              Write(7,*) nom_col(jpl)
              do ip = 1,nphi(ipl)
                Write(7,210) angle(ip), Icircdcor(ie,ip,j+i)
              end do
            endif
          end do
        endif

      end do
    end do

  endif

  if( Scan_true .and. .not. Dafs_bio ) Close(7)

  end do boucle_conv_file
 
  deallocate( Abs_U_iso, angle, As, Convolution_out_all, decal, decal_initl, dph, dpht )
  deallocate( Efermip, En_fermi, Eph1, Ephm, Epsii, Es, Eseuil )
  deallocate( f0, f0scan, fi, Fichin, Fichscanin, Fichscanout_all, fpp_avantseuil, fr )
  deallocate( hkl_dafs )
  deallocate( Icirc, Icirccor, Icircdcor, Icor, Idcor, indf )
  deallocate( Mu_mat_comp, Mus_mat, n_div_fpp, natomsym_f, ne, ne_initl, ninitl, nom_col, nphi, nsup, numat )
  deallocate( p1f, p2f, Pds, phdtscan )
  deallocate( Run_done, Stokes_name, Stokes_param, Skip_run )
  deallocate( Trunc, V0muf, Xs )
  if( Abs_in_bulk ) deallocate( As_bulk, f0_bulk, f0scan_bulk, l_dafs, Length_abs, Length_rel, Ts  )
  if( Cor_abs ) deallocate( mus )

  if( ncal > 1 .and. ical == ncal+1 ) then
    do jfich = 1,mfich
      iscr = 100 + jfich
      Close( iscr )
    end do
  endif

  if( .not. Another_one ) exit

  end do

  return

  110 format(///' Line not understood in the indata file :'//1x,A//, &
          ' If it is supposed to be a keyword check the spelling.',/ &
          ' If it is a line containing numbers, check:',/ &
          '   - How many numbers are supposed to be there.',/ &
          '   - Numbers must be separated by spaces.',/ &
          '   - Are there extra characters (points, tabulations...)'//)
  112 format(//' Error opening the file:',//3x,A,//  &
           5x,'It does not exist but the following file is found existing: ',//3x,A,// &
           5x,'Modify the input file and eventualy add other indata files with the convenient indexes !',//)
  120 format(///' The output file has the same name than one of the input files:',/ &
                3x,A,// &
                ' Choose another one with keyword "Conv_out" !',//)
  145 format(///' Error under the keyword Par_',a6,'in the indata file:',// &
                ' The wanted file is the number',i3,' !',/ &
                ' There are only',i3,' files in the job !'//)
  150 format('    Gamma_max  =',f7.2,',  Aseah =',f7.2)
  160 format('    Gamma_max  =',f7.2,',  Ecent =',f7.2,', Elarg =',f7.2)
  170 format('    Gamma_hole =',f7.2,', E_cut =',f8.3,', Shift =',f8.3,' eV')
  172 format('    Gamma_hole =',f7.2,', E_cut =',f8.3,', Shift =',f8.3,' eV,  initl =',i3)
  174 format('    Gamma_hole =',f7.2,', E_cut =',f8.3,', Shift =',f8.3,' eV,  site', i3)
  176 format('    Gamma_hole =',f7.2,', E_cut =',f8.3,', Shift =',f8.3,' eV,  site', i3,', initl =',i3)
  178 format(/'      E_(eV)      Width_(eV) lambda_(A)')
  179 format(/' E_(eV) - E_edge  Width_(eV) lambda_(A)')
  180 format(3f12.3)
  185 format(2f12.3,1p,e12.4)
  190 format(/'   Energy   ',9(' Mu_ss_r_',i1,'    Mu_ss_i_',i1,'    Mu_sp_r_',i1,'    Mu_sp_i_',i1,'    Mu_ps_r_',i1, &
                              '    Mu_ps_i_',i1,'    Mu_pp_r_',i1,'    Mu_pp_i_',i1,'   '))
  200 format(f10.3,1p,240e13.5)
  210 format(f7.1,1p,3e13.5)
end

!***********************************************************************************************************************

subroutine n_div_fpp_avanseuil(Eseuil,n_div_fpp,nfich)

  use declarations
  implicit none
  
  integer:: i, ifich, jfich, nfich

  integer, dimension(nfich):: n_div_fpp
  
  logical, dimension(nfich):: ok
  
  real(kind=db), dimension(nfich):: Eseuil
   
  ok(:) = .false.
  do ifich = 1,nfich
    if( ok(ifich) ) cycle
    i = 1
    do jfich = ifich+1,nfich
      if( ok(jfich) ) cycle
      if( abs( Eseuil(ifich) - Eseuil(jfich) ) > 1._db ) cycle
      i = i + 1
    end do 
    do jfich = ifich,nfich
      if( ok(jfich) ) cycle
      if( abs( Eseuil(ifich) - Eseuil(jfich) ) > 1._db ) cycle
      n_div_fpp(jfich) =  i
      ok(jfich) = .true.
    end do 
    
  end do
  return
end

!********************************************************************************************************************

! Evaluation of relative shifts between input data

  subroutine Shift_eval(decal,decal_initl,delta_edge,Energphot,Eph1,Ephm,Epsii,Epsii_ref,Epsii_ref_man,Eseuil,Esmin,nfich, &
                        ninit1,ninitl,ninitlm)
  
  use declarations
  implicit none

  integer:: ifich, jfich, n, nfich, ninit1, ninitlm
  
  integer, dimension(nfich):: ninitl
  
  logical:: Energphot, Epsii_ref_man
  logical, dimension(:), allocatable:: ok
  
  real(kind=db):: decal_min, delta_edge, Epsii_moy, Epsii_ref, Esmin
  real(kind=db), dimension(nfich):: decal, decal_in, Eph1, Ephm, Eseuil 
  real(kind=db), dimension(ninitlm,nfich):: decal_initl, Epsii 

! At this stage contains only the imposed value.  
  decal_in(:) = decal(:)
  decal(:) = 0._db
  
  Esmin = Eseuil(1)
  do ifich = 2,nfich
    Esmin = min( Esmin, Eseuil(ifich) ) 
  end do

  Energphot = ( Ephm(1) > Eseuil(1) ) .and. ( Eseuil(1) > Ephm(1) - Eph1(1) )

  if( .not. Energphot ) then
    do ifich = 1,nfich
      decal(ifich) = decal(ifich) + Eseuil(ifich) - Esmin
    end do
  endif

  allocate( ok(nfich) )
  ok(:) = .false.
 
  do ifich = 1,nfich
    if( ok(ifich) ) cycle
    ok(ifich) = .true.

    Epsii_moy = sum( Epsii(1:ninitl(ifich),ifich) )
    n = ninitl(ifich) 

    do jfich = ifich+1,nfich
      if( abs( Eseuil(jfich) - Eseuil(ifich) ) > 1._db ) cycle
      n = n + ninitl(jfich) 
      ok(jfich) = .true.
      Epsii_moy = Epsii_moy + sum( Epsii(1:ninitl(jfich),jfich) )
    end do
    Epsii_moy = Epsii_moy / n

    if( Epsii_ref_man ) Epsii_moy = Epsii_ref

    decal_min = 0._db
    do jfich = ifich,nfich
      if( abs( Eseuil(jfich) - Eseuil(ifich) ) > 1._db ) cycle
      decal(jfich) = decal(jfich) - Epsii_moy + sum( Epsii(1:ninitl(jfich),jfich) ) / ninitl(jfich)
      decal_min = min( decal(jfich), decal_min ) 
    end do
    do jfich = ifich,nfich
      if( abs( Eseuil(jfich) - Eseuil(ifich) ) > 1._db ) cycle
      decal(jfich) = decal(jfich) - decal_min
    end do
    
  end do
 
  deallocate( ok )
  
  decal(:) = decal(:) + decal_in(:)
  
  decal_initl(:,:) = 0._db
  do ifich = 1,nfich
    select case( ninitl(ifich) )
      case(1)
        Epsii_moy = Epsii(1,ifich)
      case(2)
         if( ninit1 == 2 ) then
           Epsii_moy = sum( Epsii(1:ninitl(ifich),ifich) ) / ninitl(ifich)
         else
           Epsii_moy = Epsii(2,ifich)
         endif
      case(4,6,10)
         if( ninit1 /= ninitl(ifich) ) then
           Epsii_moy = sum( Epsii(ninit1+1:ninitl(ifich),ifich) ) / ( ninitl(ifich) - ninit1 )
         else
           Epsii_moy = sum( Epsii(1:ninitl(ifich),ifich) ) / ninitl(ifich)
         endif
    end select
    decal_initl(1:ninitl(ifich),ifich) = Epsii(1:ninitl(ifich),ifich) - Epsii_moy + decal(ifich)
    if( abs(Delta_edge) > eps10 .and. ( ninitl(ifich) /= ninit1 ) ) decal_initl(1:ninit1,ifich) &
               = decal_initl(1:ninit1,ifich) + Delta_edge * Rydb
  end do

  return
end

!********************************************************************************************************************

subroutine Output_Energy_Grid(decal_initl,Energphot,Eph1,Es_temp,Eseuil,Esmin,Estart,fichin,Length_line,ne,ne_initl,nes,nes_in, &
                              nfich,ninitl,ninitlm,nsup,pasdeb)

  use declarations
  implicit none
  
  integer:: i, ie, ifich, ifichref, initl, initlref, ipr, je, istat, Length_line, n, &
            nemax, nes, nes_in, nfich, ninitlm, nnombre 

  integer, dimension(nfich):: ne, ninitl
  integer, dimension(ninitlm,nfich):: ne_initl, nsup
  
  character(len=Length_name), dimension(nfich):: Fichin
  
  logical:: Energphot, File_change
  logical, dimension(:,:), allocatable:: Fichdone 
  
  real(kind=db):: d, dmin, E, Emax, Emin, Esmin, Estart, Estart_l, pasdeb
  real(kind=db), dimension(nes_in):: Es_temp 
  real(kind=db), dimension(nfich):: Eph1, Eseuil 
  real(kind=db), dimension(ninitlm,nfich):: decal_initl 
  real(kind=db), dimension(:), allocatable:: Energ 
  real(kind=db), dimension(:,:,:), allocatable:: Ef 
  
  Estart_l = Estart
  do ifich = 1,nfich
    do initl = 1,ninitl(ifich)
      if( Energphot ) then
        E = Eph1(ifich) - Eseuil(ifich) + decal_initl(initl,ifich)
      else
        E = Eph1(ifich) + decal_initl(initl,ifich) - Eseuil(ifich) + Esmin
      endif
      Estart_l = Min( Estart_l, E )
    end do
  end do

  do ifich = 1,nfich
    do initl = 1,ninitl(ifich)
      if( Energphot ) then
        E = Eph1(ifich) - Eseuil(ifich) + decal_initl(initl,ifich)
      else
        E = Eph1(ifich) + decal_initl(initl,ifich) - Eseuil(ifich) + Esmin
      endif
      if( E > Estart_l - eps10 ) then
        nsup(initl,ifich) = nint( ( E - Estart_l ) / pasdeb )
        ne_initl(initl,ifich) = ne(ifich) + nsup(initl,ifich)
      else
        ne_initl(initl,ifich) = ne(ifich)
        nsup(initl,ifich) = 0
      endif
    end do
  end do

  nemax = maxval( ne_initl )

  n = sum( ninitl(:) )
  allocate( Ef(nemax,ninitlm,nfich) )

  do ifich = 1,nfich

    allocate( Energ( ne(ifich) ) )
    
    open(2, file = fichin(ifich), status='old', iostat=istat)

    read(2,*)
    do i = 1,6
      n = nnombre(2,Length_line)
      if( n == 0 ) then
        read(2,*)
        exit
      endif
      read(2,*)
    end do

    do ie = 1,ne(ifich)
      Read(2,*) Energ(ie)
    end do

    Close(2)
    
    Energ(:) = Energ(:) / Rydb
    
    do initl = 1,ninitl(ifich)
      do ie = 1,ne(ifich)
        je = ie + nsup(initl,ifich)
        Ef(je,initl,ifich) = Energ(ie)
      end do
    end do

    do initl = 1,ninitl(ifich)
      do ie = nsup(initl,ifich),1,-1
        Ef(ie,initl,ifich) = Ef(ie+1,initl,ifich) - pasdeb
      end do
    end do

    do initl = 1,ninitl(ifich)
      n = ne_initl(initl,ifich)
      Ef(1:n,initl,ifich) = Ef(1:n,initl,ifich) + decal_initl(initl,ifich)
    end do

    deallocate( Energ )
    
  end do

 ! When there are absorbing atom with different Z, 
 ! the maximum energy of the convoluted spectra corresponds to the one of maximum edge energy  
  E = Eseuil(1)
  i = 0
  do ifich = 2,nfich
    if( abs( Eseuil(ifich) - E ) < 1._db ) cycle
    if( ifich == 2 ) i = 1
    if( Eseuil(ifich) < E ) cycle
    E = Eseuil(ifich)
    i = ifich
  end do

  Emin = Ef(1,1,1)
  if( i == 0 ) then
    Emax = Ef(ne_initl(1,1),1,1)
  else
    Emax = Ef(ne_initl(1,i),1,i)
  endif

  Emin = Estart_l
  
  do ifich = 1,nfich
    do initl = 1,ninitl(ifich)
!      Emin = min( Ef(1,initl,ifich), Emin )
      if( i == 0 .or. abs( Eseuil(ifich) - Eseuil(max(i,1)) ) < 1._db ) then
        if( Ef(ne_initl(initl,ifich),initl,ifich) - Emax > Emax - Emin ) then
          Emax = Ef(ne_initl(initl,ifich),initl,ifich)
        elseif( .not. Emax - Ef(ne_initl(initl,ifich),initl,ifich) > Ef(ne_initl(initl,ifich),initl,ifich) - Emin ) then  
          Emax = min( Ef(ne_initl(initl,ifich),initl,ifich), Emax )
        endif
      endif
    end do
  end do

! Grid elaboration

  allocate( Fichdone(ninitlm,nfich) )
  Fichdone(:,:) = .false.

  dmin = 1000000000._db
  do ifich = 1,nfich
    do initl = 1,ninitl(ifich)
      if( decal_initl(initl,ifich) > dmin - eps10 ) cycle
      ifichref = ifich
      initlref = initl
      dmin = decal_initl(initlref,ifichref) 
    end do
  end do
  Fichdone(initlref,ifichref) = .true.

  je = 0
  ie = 0
  do
    je = je + 1
    if( je <= ne_initl(initlref,ifichref) ) then 
      E = Ef(je,initlref,ifichref)
      ie = ie + 1

      if( ie > nes_in ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,110)
        end do
        stop
      endif

      Es_temp(ie) = E
      if( E > Emax - eps10 ) exit
    endif

    if( ie == 1 ) cycle 

    d = decal_initl(initlref,ifichref)
    dmin = 1000000000._db
    File_change = .false.
    do ifich = 1,nfich
      do initl = 1,ninitl(ifich)
        if( Fichdone(initl,ifich) ) cycle
        if( decal_initl(initl,ifich) < d + eps10 ) then
          Fichdone(initl,ifich) = .true.
          cycle
        endif  
        if( decal_initl(initl,ifich) > d + eps10 .and. decal_initl(initl,ifich) < dmin &
           .and. ( E > Ef(nsup(initl,ifich)+1,initl,ifich) - eps10 .or. je > ne_initl(initlref,ifichref) ) ) then
          File_change = .true.
          initlref = initl
          ifichref = ifich
          dmin = decal_initl(initl,ifich)
        endif
      end do
    end do
    if( File_change ) then
      Fichdone(initlref,ifichref) = .true.
      do je = 1,ne_initl(initlref,ifichref)
        if( Ef(je,initlref,ifichref) > Es_temp(ie-1) + eps10 ) exit
      end do
      if( je > ne_initl(initlref,ifichref) ) exit
      Es_temp(ie) = Ef(je,initlref,ifichref)  
    endif

  end do
  
  if( E > Emax + eps10 ) then
    nes = ie - 1
  else
    nes = ie
  endif

  deallocate( Fichdone )
  deallocate( Ef )

  if( nes == 0 ) then
    call write_error
    do ipr = 6,9,3
      write(ipr,120)
    end do
    stop
  endif

  return
  110 format(///' Problem in the energy grid elaboration in convolution',/ &
                ' Too much points !' //)
  120 format(///' Taking into account the energy shifts, there is no overlap between',/ &
                ' the energy ranges of the different absorbing atoms !',/ &
                ' No summation and convolution are possible !' //)
end

!********************************************************************************************************************

! Name of the output file

  subroutine Conv_out_name(bav_open,check_conv,Chem,Chemin,convolution_out,Convolution_out_all,Dafs_bio,Deuxieme,fichin, &
                           fichscanin,fichscanout,fichscanout_all,Length_line,nfich,nomfichbav,Photoemission,Scan_a,Scan_true)

  use declarations
  implicit none

  integer:: ifich, istat, l, Length_line, long, longf, n, nfich, nnombre, ns
  
  character(len=8):: dat
  character(len=10):: tim
  character(len=50):: com_date, com_time
  character(len=Length_name):: chemin, convolution_out, fichscanout, mot, mots, nomfichbav
  character(len=Length_name), dimension(nfich):: convolution_out_all, fichin, fichscanin, fichscanout_all

  logical:: bav_open, check_conv, Chem, Dafs_bio, Deuxieme, Photoemission, Scan_a, Scan_true
  
  if( convolution_out == ' ' ) then
 
    do ifich = 0,nfich
      if( nfich == 1 .and. ifich > 0 ) exit
 
      mot = fichin( max(1,ifich) )
      l = len_trim( mot )
      if( l > 4 ) then
        if( mot(l-3:l) == '.txt' )  mot(l-3:l) = '    '     
      endif
      
      l = len_trim( mot)
      if( nfich > 1 .and. ifich == 0 .and. l > 2 ) then
        if( mot(l-1:l-1) == '_' ) then
          mot(l-1:l) = '  '
        elseif( mot(l-2:l-2) == '_' ) then
          mot(l-2:l) = '   '
        endif
      endif

      l = len_trim( mot) + 1
      if( Photoemission ) then
        mot(l:l+12) = '_xes_conv.txt'
      else
        mot(l:l+8) = '_conv.txt'
      endif
      
      if( ifich == 0 ) then
        convolution_out = mot
      else
        convolution_out_all(ifich) = mot
      endif

    end do      
  
  else
  
    l = len_trim( convolution_out )

    if( l > 5 ) then
      if( convolution_out(l-3:l) /= '.txt' ) convolution_out(l+1:l+4) = '.txt'  
    endif

    do ifich = 1,nfich
      if( nfich == 1 ) exit
      mot = convolution_out
      l = len_trim( convolution_out )
      if( convolution_out(l-8:l) == '_conv.txt' ) then
        mot(l-7:l) = '        '
        call ad_number(ifich,mot,Length_name)
        l = len_trim( mot )
        if( Photoemission ) then
          mot(l+1:l+16) = '_xes_conv.txt'
        else
          mot(l+1:l+9) = '_conv.txt'
        endif 
      else
        if( Photoemission ) then
          mot(l-3:l+1) = '_xes_'
        else
          mot(l-3:l) = '_   '
        endif
        call ad_number(ifich,mot,Length_name)
        l = len_trim( mot )
        mot(l+1:l+4) = '.txt'
      endif 
      convolution_out_all(ifich) = mot
    end do
    
  endif

  if( .not. bav_open .and. Check_conv ) then
    l = len_trim(convolution_out)
    nomfichbav = ' '
    nomfichbav(1:l-4) = convolution_out(1:l-4)
    nomfichbav(l-3:l+4) = '_bav.txt'
    open(3, file = nomfichbav, status='unknown',iostat=istat)
    if( istat /= 0 ) call write_open_error(nomfichbav,istat,1)
    bav_open = .true.
    call date_and_time( date = dat, time = tim )
    com_date = '   Date = ' // dat(7:8) // ' ' // dat(5:6) // ' ' // dat(1:4)
    com_time = '   Time = ' // tim(1:2)  // ' h ' // tim(3:4) // ' mn ' // tim(5:6) // ' s'
    write(3,'(1x,A/A/A)') Revision, com_date, com_time
  endif

  if( Check_conv .and. .not. Deuxieme ) write(3,110)

  if( Scan_a .and. .not. Scan_true ) then
    do ifich = 1,nfich
      mots = ' '
      mot = fichin(ifich)
      l = len_trim(mot) - 4
      if( nfich == 1 ) then
        ns = 0
      else
        if( mot(l-2:l-2) == '_' ) then
          ns = 3
        elseif( mot(l-1:l-1) == '_' ) then
          ns = 2
        else
          ns = 0
        endif
      endif
      mots(1:l-ns) = mot(1:l-ns)
      mots(l-ns+1:l-ns+5) = '_scan'
      if( ns > 0 ) mots(l-ns+6:l+5) = mot(l-ns+1:l)
      mots(l+6:l+9) = '.txt'
      fichscanin(ifich) = mots
    end do
    Scan_true = .true.
  endif

  Dafs_bio = .false.
  if( Scan_true ) then
    open(2, file = fichscanin(1), status='old',iostat=istat)
    if( istat /= 0 ) call write_open_error(fichscanin(1),istat,1)
    n = nnombre(2,Length_line)
    read(2,*) n
    if( n == 4 ) Dafs_bio = .true.
    close(2)
  endif

  if( Scan_true .and. fichscanout == ' ' .and. .not. Dafs_bio ) then
    mot = convolution_out
    l = len_trim(mot)
    mot(l-7:l+5) = 'scan_conv.txt'
    fichscanout = mot
    do ifich = 1,nfich
      if( nfich == 1 ) exit
      mot = convolution_out_all(ifich)
      l = len_trim(mot)
      mot(l-7:l+5) = 'scan_conv.txt'
      fichscanout_all(ifich) = mot
    end do
  endif

  if( Chem ) then
    long = len_trim(chemin)
    do ifich = 1,nfich
      mot = fichin(ifich)
      longf = len_trim(mot)
      fichin(ifich) = chemin(1:long) // mot(1:longf)
    end do
    longf = len_trim(convolution_out)
    convolution_out = chemin(1:long) // convolution_out(1:longf)
    do ifich = 1,nfich
      if( nfich == 1 ) exit
      longf = len_trim(convolution_out_all(ifich))
      mot = convolution_out_all(ifich)
      convolution_out_all(ifich) = chemin(1:long) // mot(1:longf)
    end do
    if( Scan_true ) then
      do ifich = 1,nfich
        mot = fichscanin(ifich)
        longf = len_trim(mot)
        fichscanin(ifich) = chemin(1:long) // mot(1:longf)
      end do
      longf = len_trim(fichscanout)
      fichscanout = chemin(1:long) // fichscanout(1:longf)
       do ifich = 1,nfich
        if( nfich == 1 ) exit
        longf = len_trim(fichscanout_all(ifich))
        mot = fichscanout_all(ifich)
        fichscanout_all(ifich) = chemin(1:long) // mot(1:longf)
      end do
   endif
  endif

  return
  110 format(/' ---- Convolution ------',100('-'))
end

!********************************************************************************************************************

! Reading of dimensions

subroutine Dimension_file(Abs_in_bulk,Abs_U_iso,Cor_abs,Eintmax,En_Fermi,Eph1,Ephm, Epsii,f0_forward,Fichin,Fichscanin, &
          fprim,Full_self_abs,Green_int,jseuil,Length_line,Magn,mu_0_bulk,n_col_max,n_mat_pol,n_Trunc,natomsym_f,ne,nfich, &
          ninit1,ninitl,ninitlm,nphim,npldafs,nseuil,nxan,nxan_lib,Sample_thickness,Scan_true,Self_abs,Signal_Sph, &
          Surface_ref,Tenseur,Tenseur_car,Trunc,V0muf,Volume_maille,Volume_maille_bulk)

  use declarations
  implicit none

  integer:: eof, i, ie, ifich, ipl, ipr, istat, jseuil, l, Length_line, n, n_col_max, n_mat_pol, n_mat_pol1, n_Trunc, nfich, &
            ninit1, ninitlm, nnombre, nphi, nphim, nphim1, npldafs, npldafs1, nseuil, numat, nxan, nxan1
  integer, dimension(nfich):: ne, ninitl

  character(len=15):: nomab
  character(len=132):: mot
  character(len=Length_name), dimension(nfich):: fichin, Fichscanin
  character(len=15), dimension(:), allocatable:: nomcol

  logical:: Cor_abs, fprim, Full_self_abs, Green_int, Magn, nxan_lib, Scan_true, Self_abs, Signal_Sph, Tenseur, Tenseur_car, &
            Abs_in_bulk
  logical, dimension(nfich):: Trunc

  real(kind=db):: Eintmax, Eph, Eseuil, f0_forward, fpp_avantseuil, mu_0_bulk, Sample_thickness, Surface_ref, Volume, &
                  Volume_maille, Volume_maille_bulk  
  real(kind=db), dimension(ninitlm,nfich):: Epsii
  real(kind=db), dimension(nfich):: Abs_U_iso, En_Fermi, Eph1, Ephm, natomsym_f, V0muf

  n_Trunc = 0
  natomsym_f(:) = 0
 
  do ifich = 1,nfich

    n_mat_pol = 0

    open(2, file = fichin(ifich), status='old', iostat=istat)

    n = nnombre(2,Length_line)
    if( n == 11 + ninitl(ifich) ) then
      read(2,*) Eseuil, numat, nseuil, jseuil, fpp_avantseuil, V0muf(ifich), En_fermi(ifich), ninitl(ifich), &
        ninit1, Epsii(1:ninitl(ifich),ifich), Volume, f0_forward
    elseif( n == 12 + ninitl(ifich) ) then
      read(2,*) Eseuil, numat, nseuil, jseuil, fpp_avantseuil, V0muf(ifich), En_fermi(ifich), ninitl(ifich), &
        ninit1, Epsii(1:ninitl(ifich),ifich), Volume, f0_forward, natomsym_f(ifich)
    elseif( n == 13 + ninitl(ifich) ) then
      read(2,*) Eseuil, numat, nseuil, jseuil, fpp_avantseuil, V0muf(ifich), En_fermi(ifich), ninitl(ifich), &
        ninit1, Epsii(1:ninitl(ifich),ifich), Volume, Surface_ref, f0_forward, natomsym_f(ifich)
    else
      read(2,*) Eseuil, numat, nseuil, jseuil, fpp_avantseuil, V0muf(ifich), En_fermi(ifich), ninitl(ifich), &
        ninit1, Epsii(1:ninitl(ifich),ifich), Volume, Surface_ref, f0_forward, natomsym_f(ifich), Abs_U_iso(ifich)
    endif

    if( ninit1 < 0 ) then
      Green_int = .true.
      ninit1 = abs( ninit1 )
    else
      Green_int = .false.
    endif
    
    n = nnombre(2,Length_line)
    
    if( n >= n_col_max ) then
      call write_error
      do ipr = 6,9,3
        write(ipr,110) n_mat_pol
      end do
      stop
    endif
    
    npldafs = n / 2
    if( n /= 0 ) then    
      read(2,*)
      read(2,*)
    endif

    n = nnombre(2,Length_line)
    if( n /= 0 ) then
      read(2,*)
      n = nnombre(2,Length_line)
      if( n == 0 ) then
        read(2,'(A)') mot
        do i = 1,131
          if( mot(i:i+1) == 'mu' ) then
            Full_self_abs = .true.
            exit
          endif
        end do
        Self_abs = i > 131
        backspace(2)
      else
        read(2,*); read(2,*) 
        Trunc(ifich) = .true.
        n_Trunc = n_Trunc + 1
        Abs_in_bulk = .true.
        mu_0_bulk = fpp_avantseuil
      endif
    endif            
    Cor_abs = Self_abs .or. Full_self_abs

    if( Trunc(ifich) ) then
      Volume_maille_bulk = Volume
    else
      Volume_maille = Volume
    endif

    read(2,'(10x,a15)') nomab
    nomab = adjustl( nomab )

    n = ( nnombre(2,Length_line) - 1 ) / ninitl(ifich)
    if( npldafs > 0 ) then
      if( Full_self_abs ) then
        nxan = n - 6 * npldafs
      elseif( Self_abs ) then
        nxan = n - 4 * npldafs
      else
        nxan = n - 2 * npldafs
      endif
      if( nomab(1:8) == 'Sum_dd_r' .or. nomab(1:8) == 'Sum_tot_' ) Signal_Sph = .true.
    elseif( fprim .and. nomab(1:1) == 'D' ) then
      Tenseur = .true.
      nxan = 0
      if( nomab(1:4) == 'D_xx' ) Tenseur_car = .true.
      Magn = nomab(1:6) == 'D_xx_r'
      if( Magn ) then
        npldafs = n / 2
      else
        npldafs = n
      endif
    else
      nxan = n
    endif

! number of polarization matrix
    backspace( 2)
    allocate( nomcol(nxan) )
    read(2,'(10x,10000a15)') nomcol(:)
    do ipl = 1,nxan
      nomab = nomcol(ipl)
      l = len_trim( nomab )
      if( l <= 3 ) cycle
      do i = 1,l-2
        if( nomab(i:i+2) /= '_ss' ) cycle
        n_mat_pol = n_mat_pol + 1
        exit
      end do
    end do
    deallocate( nomcol )
    
    if( n_mat_pol /= 0 .and. Sample_thickness < eps10 ) then
      call write_error
      do ipr = 6,9,3
        if( n_mat_pol == 1 ) then
          write(ipr,120) n_mat_pol
        else
          write(ipr,122) n_mat_pol
        endif 
      end do
      stop
    endif

    nxan = nxan - 6 * n_mat_pol
     
    if( ifich == 1 ) then
      n_mat_pol1 = n_mat_pol
      nxan1 = nxan
      npldafs1 = npldafs
    else
      if( nxan_lib ) nxan = min(nxan1,nxan)
      if( nxan1 /= nxan .or. npldafs1 /= npldafs ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,130) nxan1, npldafs1, n_mat_pol1, ifich, nxan, npldafs, n_mat_pol
        end do
        stop
      endif
    endif

    do ie = 1,10000
      Read(2,*,iostat=eof) Eph
      if( eof /= 0 ) exit
      if( ie == 1 ) eph1(ifich) = eph
      Ephm(ifich) = eph
      if( eph > Eintmax ) exit
    end do
    ne(ifich) = ie - 1

    if( ne(ifich) < 2 ) then
      call write_error
      do ipr = 6,9,3
        write(ipr,140) fichin(ifich)
      end do
      stop
    endif

    Close(2)
  end do

  if( npldafs > 0 ) then 
    nphim = 1
  else
    nphim = 0
  endif
  if( Scan_true ) then
    do ifich = 1,nfich
      open(2, file = fichscanin(ifich), status='old',iostat=istat)
      if( istat /= 0 ) call write_open_error(fichscanin(ifich),istat,1)
      n = nnombre(2,Length_line)
      do ipl = 1,npldafs
        read(2,*) nphi
        nphim = max( nphim, nphi )
      end do
      Close(2)
      
      if( ifich == 1 ) then
        nphim1 = nphim
      elseif( nphim /= nphim1 ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,150) nphim1, nphim
        end do
        stop
      endif
 
    end do
  endif

  return
  110 format(///' The maximum of number of columns =',i6,', is reached in the files to convolute.', & 
               ' Use the keyword Length_line to increase this number.', &
               ' Some compilers such as ftn95 does not allow higher number.',// ) 
  120 format(///'     In convolution, there is', i2,' polarization matrix to calculate',/&
                '     A sample tickness must be given in micrometer with keyword "Sample_tickness" !',// )
  122 format(///'     In convolution, there are', i3,' polarization matrix to calculate',/&
                '     A sample tickness must be given in micrometer with keyword "Sample_tickness" !',// )
  130 format(///'     In convolution, indata files are different,',/ &
                '  file   1 : number of xanes =',i4,', number of reflexions =',i4,', number of polarization matrix =',i4,/ &
                '  file',i4,' : number of xanes =',i4,', number of reflexions =',i4,', number of polarization matrix =',i4)
  140 format(//' For the convolution, the number of energy must be greater than one !'/, &
               ' It is not the case in the file:'//A//)
  150 format(///'     Input files are different,',/ &
                '  file 1 : nphim =',i4,'  file ',i2,' nphim =',i4)
end

!***********************************************************************

subroutine col_name(Analyzer,Bormann,Cor_abs,Dafs_bio,Double_cor,fichin,Fichscanin,fprim,Full_self_abs,hkl_dafs, &
      Length_line,n_col,n_mat_pol,n_stokes,nom_col,npldafs,npldafs_b,nxan,Photoemission,Self_abs,Signal_sph,Stokes, &
      Stokes_name,Stokes_param,Sup_sufix,Tenseur)

  use declarations
  implicit none

  integer:: i, i_hk, ii, ipl, istat, j, k, kk, l, Length_line, m, n, n_col, n_col_o, n_col_in, n_index_hk, n_mat_pol, &
    n_stokes, nc, nnombre, npldafs, npldafs_b, nxan
  integer, dimension(npldafs):: nphi

  character(len=Length_name):: fichin, Fichscanin
  character(len=Length_word):: nomab, nomac
  character(len=length_line):: motl
  character(len=13), dimension(n_stokes):: Stokes_name
  character(len=Length_word), dimension(n_col):: nom_col
  character(len=Length_word), dimension(2*n_col):: nom_col_o
  character(len=Length_word), dimension(:), allocatable:: nom_col_in

  integer, dimension(3,npldafs_b):: hkl_dafs

  logical:: Analyzer, Bormann, Cor_abs, Dafs_bio, Double_cor, fprim, Full_self_abs, Photoemission, Self_abs, Signal_sph, &
            Stokes, Sup_sufix, Tenseur

  real(kind=db), dimension(5,n_stokes):: Stokes_param

  open(2, file = fichin, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(fichin,istat,1)

  do i = 1,5
    n = nnombre(2,Length_line)
    if( n == 0 ) exit
    read(2,*)
  end do

  if( Dafs_bio ) then
    n_col_in = 1 + nxan
  else
    n_col_in = 1 + nxan + 6 * n_mat_pol + 2 * npldafs
    if( Self_abs ) n_col_in = n_col_in + 2*npldafs
    if( Full_self_abs ) n_col_in = n_col_in + 4*npldafs
  endif

  allocate( nom_col_in(n_col_in) )

  read(2,'(A)') motl

  call Extract_word(Length_line,motl,nom_col_in,n_col_in,n_col_in)

  if( Sup_sufix .or. nxan /= 0 .or. n_mat_pol /= 0 ) then
    if( nxan == 0 .and. n_mat_pol == 0 ) then
      nc = n_col_in - 1
    else
      nc = nxan + 6 * n_mat_pol 
    endif 
    i = 0
    j = 1   ! Energy
    do ipl = 1,nc
      i = i + 1
      j = j + 1
      nomab = nom_col_in(j)
      l = len_trim( nomab )
      if( nomab(l-1:l-1) == '_' ) nomab(l-1:l) = '  '
      if( nomab(l-2:l-1) == '_L' ) nomab(l-2:l) = '   '
      if( nomab(l-2:l-1) == '_M' ) nomab(l-2:l) = '   '
      if( nomab(l-2:l-1) == '_N' ) nomab(l-2:l) = '   '
      if( nomab(l-2:l-1) == '_O' ) nomab(l-2:l) = '   '

      if( j > nxan + 1 .and. j <= nxan + 6 * n_mat_pol + 1 ) then

        if( mod(j-nxan-1,6) /= 1 ) then
          i = i - 1
          cycle
        endif
        l = len_trim( nomab )
        do k = 1,l
          if( nomab(k:k) == '_' ) exit
        end do
        do m = k+1,l
          nomab(m:m) = ' '
        end do
        nomac = nomab
        l = len_trim( nomab )
        nomab(l+1:l+2) = 'ss'
        nom_col_o(i) = nomab
        nomab(l+1:l+2) = 'pp'
        i = i + 1
        nom_col_o(i) = nomab
        do n = 1,n_stokes
          if( Stokes_name(n) /= 'no_name' ) then
            nomab = Stokes_name(n)
          else
           nomab = nomac
           call ad_number(n,nomab,Length_word)
          endif
          i = i + 1
          nom_col_o(i) = nomab
        end do
      else 
        if( l > 4 ) then
          if( nomab(l-3:l-2) == '_r' .or. nomab(l-3:l-2) == '_i' ) nomab(l-3:l-2) = '  '
        endif
        if( nxan == 0 ) then
           nom_col_in(j) = nomab
        else
          nom_col_o(i) = nomab
        endif
      endif
      nomab = adjustl( nom_col_o(i) )
      if( Photoemission .and. nomab(1:7) == '<xanes>' ) then
        nomab(1:7) = '  <xes>'
        call center_word(nomab,Length_word)
        nom_col_o(i) = nomab
      endif
       
    end do
  endif

  i = nxan + n_stokes * n_mat_pol
  j = nxan + 6 * n_mat_pol + 1

  if( Dafs_bio ) then

! When using Dafs_bio, there is necessarily a scan file with 4 polarizations value
! We use this file to get the hkl_dafs indexes
    open(2, file = fichscanin, status='old',iostat=istat)

    do ipl = 1,npldafs
      read(2,*) nphi(ipl)
    end do
        
    read(2,*)
    do ipl = 1,npldafs
      read(2,*) hkl_dafs(:,ipl)
      do k = 1,nphi(ipl)
        read(2,*) 
      end do
    end do
 
    Close(2)
 
    if( n_col == nxan + 6 * n_mat_pol + npldafs ) then
      n = 1
    else
      n = 5
    endif

    do ipl = 1,npldafs
      do j = 1,n
        i = i + 1
        nomab = '('
        do k = 1,3
          if( j == 4 .or. j == 5 ) then
            call ad_number(-hkl_dafs(k,ipl),nomab,Length_word)
          else
            call ad_number(hkl_dafs(k,ipl),nomab,Length_word)
          endif
          l = len_trim(nomab) + 1
          if( l > Length_word ) exit
          if( k == 3 ) then
            nomab(l:l) = ')'
          else
            nomab(l:l) = ','
          endif
        end do
        if( l+1 <= Length_word ) then 
          if( j == 2 .or. j == 4 ) then
            nomab(l+1:l+1) = 's'
          elseif( j == 3 .or. j == 5 ) then
            nomab(l+1:l+1) = 'p'
          endif
        endif
        call center_word(nomab,Length_word)
        nom_col_o(i) = nomab
      end do
    end do

  else

    do ipl = 1,npldafs

      if( Tenseur ) then

        i = i + 1
        j = j + 1
        nomab = nom_col_in(j)
        l = len_trim( nomab )
        nomab(l+1:l+1) = 'p'
        nom_col_o(i) = nomab

        i = i + 1
        nomab(l+1:l+1) = 's'
        nom_col_o(i) = nomab

        cycle

      elseif( Bormann ) then

        i = i + 1
        j = j + 2
        nomab = nom_col_in(j)
        nomab(1:1) = 'e'
        nom_col_o(i) = nomab

        i = i + 1
        l = min( len_trim( nomab ) + 1, Length_word )
        nomab(l:l) = 'i'
        nom_col_o(i) = nomab

        cycle

      endif

      j = j + 2
      nomab = nom_col_in(j)
      if( Signal_sph ) then
        l = len_trim( nomab )
        nomab(2:l+1) = nomab(1:l)
      endif
      nomab(1:1) = 'I'
      i = i + 1
      nom_col_o(i) = nomab

      if( fprim ) then
        i = i + 1
        nomab = nom_col_in(j)
        l = len_trim( nomab )
        do k = 1,l-1
          nomab(k:k) = nomab(k+1:k+1)
        end do
        nomab(l:l) = 'p'
        nom_col_o(i) = nomab

        i = i + 1
        nomab(l:l) = 's'
        nom_col_o(i) = nomab
      endif

      if( Cor_abs ) then
        i = i + 1
        nomab = nom_col_in(j)
        l = len_trim( nomab )
        do k = min(l+1,Length_word),2,-1
          nomab(k:k) = nomab(k-1:k-1)
        end do
        nomab(1:2) = 'Ic'
        nom_col_o(i) = nomab
      endif

      if( Double_cor ) then
        i = i + 1
        nomab = nom_col_o(i-1)
        nomab(2:2) = 'd'
        nom_col_o(i) = nomab
      endif

      if( Cor_abs ) then
        do ii = 1,4
          if( .not. Full_self_abs .and. ii > 2 ) exit
          i = i + 1
          j = j + 1
          nom_col_o(i) = nom_col_in(j)
        end do
      endif

      if( Stokes .and. mod(ipl,4) == 0 ) then
        do ii = 1,n_stokes

          i = i + 1

          nomab = ' '

          if( Stokes_name(ii) /= 'no_name' ) then

            nomab = Stokes_name(ii)
            l = len_trim( nomab )
            nomab(2:l+1) = nomab(1:l)
            nomab(1:1) = 'I'
            l = l + 1

          else

            nomab = nom_col_in(j-4)
            nomab(1:1) = 'I'
            l = len_trim( nomab )
            do k = 2,l-1
              if( nomab(k:k+1) == 'pp' ) exit
            end do
            if( abs( Stokes_param(1,ii) - 1 ) < 0.0001 ) then
              nomab(k:k) = 's'
            elseif( abs( Stokes_param(1,ii) + 1 ) < 0.0001 ) then
              nomab(k:k) = 'p'
            elseif( abs( Stokes_param(3,ii) - 1 ) < 0.0001 ) then
              nomab(k:k) = 'r'
            elseif( abs( Stokes_param(3,ii) + 1 ) < 0.0001 ) then
              nomab(k:k) = 'l'
            else
              nomab(k:k) = 'a'
            endif

            k = k + 1

            if( abs( Stokes_param(5,ii) ) < 0.0001 ) then
              nomab(k:l-1) = nomab(k+1:l)
              nomab(l:l) = ' '
            elseif( abs( Stokes_param(4,ii) ) < 0.0001 ) then
              nomab(k:k) = 's'
            elseif( abs( Stokes_param(4,ii) - pi/2 ) < 0.0001 ) then
              nomab(k:k) = 'p'
            else
              nomab(k:k) = 'a'
            endif

          endif

          nom_col_o(i) = nomab

          if( Cor_abs ) then
            i = i + 1
            do k = min(l+1,Length_word),2,-1
              nomab(k:k) = nomab(k-1:k-1)
            end do
            nomab(1:2) = 'Ic'
            nom_col_o(i) = nomab
          endif

          if( Double_cor ) then
            i = i + 1
            nomab(2:2) = 'd'
            nom_col_o(i) = nomab
          endif

        end do
      endif

    end do

  endif

  Close(2)
  
  deallocate( nom_col_in )

  if( Analyzer ) then
   
    nom_col(1:n_col) = nom_col_o(1:n_col)

  else
    n_col_o = i
  
    i_hk = 0
    k = 0
    kk = 0
    do i = 1,n_col_o
    
      if( i <= k ) cycle
      
      nomab = adjustl( nom_col_o(i) )
      if( nomab(1:2) /= 'I(' ) then
        k = k + 1
        kk = kk + 1 
        nom_col(kk) = nom_col_o(i)  
        cycle
      endif
      l = len_trim( nomab ) - 1
    
      do j = i+1,n_col_o
        nomac = adjustl( nom_col_o(j) )
        if( nomab(1:l) /= nomac(1:l) ) cycle
        i_hk = i_hk + 1
        n_index_hk = j - i
        k = k + 2 * n_index_hk 
        exit
      end do
      
      if( j > n_col_o ) then
        kk = kk + 1 
        nom_col(kk) = nom_col_o(i)  
      else
      
        do j = i, i + n_index_hk - 1
          kk = kk + 1  
          nomab = nom_col_o(j)
          l = len_trim( nomab )
          nomab(l-1:l) = '  '
          nom_col(kk) = nomab
        end do
      
      endif  

    end do

  endif
   
  return
end

!***********************************************************************

! Tables donnant les largeurs des niveaux de coeur 
! Seuils K, L1, L2 et L3: vient de M. O. Krause et J. H. Oliver, J. Phys. Chem. Ref. Data 8, 329 (1979)

function Tab_width(Eseuil,jseuil,nseuil,numat)

  use declarations
  implicit none
  
  integer, parameter:: nam = 110
  
  integer:: jseuil, nseuil, numat
  
  real(kind=db):: Es, Eseuil, Gamma_hole, Shift_edge, Tab_width
  real(kind=db), dimension(10:nam):: GH_K1, GH_L1, GH_L2, GH_L3

  data GH_K1/                                                                                  0.24_db, &
     0.30_db,  0.36_db,  0.42_db,  0.48_db,  0.53_db,  0.59_db,  0.64_db,  0.68_db,  0.74_db,  0.81_db, &
     0.86_db,  0.94_db,  1.01_db,  1.08_db,  1.16_db,  1.25_db,  1.33_db,  1.44_db,  1.55_db,  1.67_db, &
     1.82_db,  1.96_db,  2.14_db,  2.33_db,  2.52_db,  2.75_db,  2.99_db,  3.25_db,  3.52_db,  3.84_db, &
     4.14_db,  4.52_db,  4.91_db,  5.33_db,  5.77_db,  6.24_db,  6.75_db,  7.28_db,  7.91_db,  8.49_db, &
     9.16_db,  9.89_db, 10.6_db,  11.4_db,  12.3_db,  13.2_db,  14.1_db,  15.1_db,  16.2_db,  17.3_db,  &
    18.5_db,  19.7_db,  21.0_db,  22.3_db,  23.8_db,  25.2_db,  26.8_db,  28.4_db,  30.1_db,  31.9_db,  &
    33.7_db,  35.7_db,  37.7_db,  39.9_db,  42.1_db,  44.4_db,  46.8_db,  49.3_db,  52.0_db,  54.6_db,  &
    57.4_db,  60.4_db,  63.4_db,  66.6_db,  69.8_db,  73.3_db,  76.8_db,  80.4_db,  84.1_db,  88.0_db,  &
    91.9_db,  96.1_db, 100._db,  105._db,  109._db,  114._db,  119._db,  124._db,  129._db,  135._db,   &
   140._db,  145._db,  150._db,  156._db,  162._db,  168._db,  174._db,  181._db,  187._db,  193._db    /

  data GH_L1/                                                                                  0.1_db,  &
     0.2_db,   0.41_db,  0.73_db,  1.03_db,  1.26_db,  1.49_db,  1.58_db,  1.63_db,  1.92_db,  2.07_db, &
     2.21_db,  2.34_db,  2.41_db,  2.54_db,  2.62_db,  2.76_db,  2.79_db,  2.89_db,  3.06_db,  3.28_db, &
     3.38_db,  3.53_db,  3.79_db,  3.94_db,  4.11_db,  4.28_db,  4.44_db,  4.67_db,  4.71_db,  4.78_db, &
     3.94_db,  4.25_db,  4.36_db,  4.58_db,  4.73_db,  4.93_db,  4.88_db,  4.87_db,  5.00_db,  2.97_db, &
     3.13_db,  3.32_db,  3.46_db,  3.64_db,  3.78_db,  3.92_db,  4.06_db,  4.21_db,  4.34_db,  4.52_db, &
     4.67_db,  4.80_db,  4.91_db,  5.05_db,  5.19_db,  5.25_db,  5.33_db,  5.43_db,  5.47_db,  5.53_db, &
     5.54_db,  5.63_db,  5.58_db,  5.61_db,  6.18_db,  7.25_db,  8.30_db,  9.39_db, 10.5_db,  11.3_db,  &
    12.0_db,  12.2_db,  12.4_db,  12.6_db,  12.8_db,  13.1_db,  13.3_db,  13.4_db,  13.6_db,  13.7_db,  &
    14.3_db,  14.0_db,  14.0_db,  13.5_db,  13.3_db,  13.6_db,  13.8_db,  14.0_db,  14.3_db,  14.4_db,  &
    14.8_db,  15.1_db,  15.9_db,  16.2_db,  16.5_db,  16.8_db,  17.0_db,  17.2_db,  17.6_db,  18.1_db   /

  data GH_L2/                                                                                  0.0_db,  &
     0.0_db,  0.001_db, 0.004_db, 0.015_db, 0.032_db, 0.054_db, 0.083_db, 0.126_db, 0.152_db,  0.17_db, &
     0.19_db,  0.24_db,  0.26_db,  0.29_db,  0.34_db,  0.37_db,  0.43_db,  0.52_db,  0.62_db,  0.72_db, &
     0.83_db,  0.95_db,  1.03_db,  1.13_db,  1.21_db,  1.31_db,  1.43_db,  1.54_db,  1.65_db,  1.78_db, &
     1.87_db,  1.97_db,  2.08_db,  2.23_db,  2.35_db,  2.43_db,  2.57_db,  2.62_db,  2.72_db,  2.84_db, &
     3.00_db,  3.12_db,  3.25_db,  3.40_db,  3.51_db,  3.57_db,  3.68_db,  3.80_db,  3.89_db,  3.97_db, &
     4.06_db,  4.15_db,  4.23_db,  4.32_db,  4.43_db,  4.55_db,  4.66_db,  4.73_db,  4.79_db,  4.82_db, &
     4.92_db,  5.02_db,  5.15_db,  5.33_db,  5.48_db,  5.59_db,  5.69_db,  5.86_db,  6.00_db,  6.17_db, &
     6.32_db,  6.48_db,  6.67_db,  6.83_db,  7.01_db,  7.20_db,  7.47_db,  7.68_db,  7.95_db,  8.18_db, &
     8.75_db,  9.32_db,  9.91_db, 10.5_db,  10.9_db,  11.4_db,  11.8_db,  12.2_db,  12.7_db,  13.1_db,  &
    13.6_db,  14.0_db,  14.4_db,  14.9_db,  15.5_db,  15.9_db,  16.4_db,  17.0_db,  17.5_db,  18.1_db   /

  data GH_L3/                                                                                  0.0_db,  &
     0.0_db,  0.001_db, 0.004_db, 0.014_db, 0.033_db, 0.054_db, 0.087_db, 0.128_db, 0.156_db,  0.17_db, &
     0.19_db,  0.22_db,  0.24_db,  0.27_db,  0.32_db,  0.36_db,  0.43_db,  0.48_db,  0.56_db,  0.65_db, &
     0.76_db,  0.82_db,  0.94_db,  1.00_db,  1.08_db,  1.17_db,  1.27_db,  1.39_db,  1.50_db,  1.57_db, &
     1.66_db,  1.78_db,  1.91_db,  2.00_db,  2.13_db,  2.25_db,  2.40_db,  2.50_db,  2.65_db,  2.75_db, &
     2.87_db,  2.95_db,  3.08_db,  3.13_db,  3.25_db,  3.32_db,  3.41_db,  3.48_db,  3.60_db,  3.65_db, &
     3.75_db,  3.86_db,  3.91_db,  4.01_db,  4.12_db,  4.17_db,  4.26_db,  4.35_db,  4.48_db,  4.60_db, &
     4.68_db,  4.80_db,  4.88_db,  4.98_db,  5.04_db,  5.16_db,  5.25_db,  5.31_db,  5.41_db,  5.50_db, &
     5.65_db,  5.81_db,  5.98_db,  6.13_db,  6.29_db,  6.41_db,  6.65_db,  6.82_db,  6.98_db,  7.13_db, &
     7.33_db,  7.43_db,  7.59_db,  7.82_db,  8.04_db,  8.26_db,  8.55_db,  8.75_db,  9.04_db,  9.33_db, &
     9.61_db,  9.90_db, 10.1_db,  10.5_db,  10.8_db,  11.2_db,  11.6_db,  12.0_db,  12.3_db,  12.7_db    /

  select case(nseuil)

    case(1)   ! Seuil K
    
      if( numat < 10 .or. numat > nam ) then
   ! Formule ajustee sur les data de Krause et Oliver (correspond donc a une extrapolation)
   ! We substract the shift
        Es = ( Eseuil - Shift_edge(0,numat) ) * Rydb
        Gamma_hole = 0.13007_db + 0.00011544_db * Es + 9.2262E-09_db * Es**2 - 3.2344E-14_db * Es**3
      else
        Gamma_hole = GH_K1( numat )
      endif
      
    case(2)

      select case(jseuil)

        case(1)   ! Seuil L1
          if( numat < 10 ) then
            Gamma_hole = 0._db
          elseif( numat >= 10 .and. numat <= nam ) then
            Gamma_hole = GH_L1( numat )
          else
            Gamma_hole = 20._db
          endif

        case(2)   ! Seuil L2
          if( numat < 10 ) then
            Gamma_hole = 0._db
          elseif( numat >= 10 .and. numat <= nam ) then
            Gamma_hole = GH_L2( numat )
          else
            Gamma_hole = 19._db
          endif
! Anciennes valeurs avant 13/2/2015
!    Z < 26     Gamma_hole = 1.699*p  + 0.5001, p = ( numat - 1 ) / 25.    
!     Fe        Gamma_hole = 1.5     
!     Co        Gamma_hole = 1.3    
!     Ni        Gamma_hole = 1.1
!  Z = 29,30,31 Gamma_hole = -0.9*p + 2.2,  p = ( numat - 26 ) / 5.     
!  31 < Z < 61  Gamma_hole = 2.7*p + 1.3,   p = ( numat - 31 ) / 29.     
!  60 < Z < 81  Gamma_hole = 1.5*p + 4.0,   p = ( numat - 60 ) / 20.
!  Z > 80       Gamma_hole = 5*p + 5.5,     p = ( numat - 80 ) / 15.

        case(3)   ! Seuil L3
          if( numat < 10 ) then
            Gamma_hole = 0._db
          elseif( numat >= 10 .and. numat <= nam ) then
            Gamma_hole = GH_L3( numat )
          else
            Gamma_hole = 13._db
          endif
! Anciennes valeurs avant 13/2/2015, comme L2 - 0.5
      end select

! Pour les seuils M, N, O, les valeurs ne sont pas optimisees
    case default

      select case(jseuil)

        case(0) ! Optic

          Gamma_hole = 0._db

        case(1,2,3)

          Gamma_hole = 0.2_db

        case(4,5,6,7)

          Gamma_hole = 0.1_db

    end select

  end select

  Tab_width = Gamma_hole / Rydb

  return
end

!***********************************************************************

subroutine Seahdench(A,E_cut,Gamma_max,nelor,Elor,betalor)

  use declarations
  implicit none
  
  integer:: ie, nelor

  real(kind=db):: A, E, E_cut, Ep, lambda, Gamma_max
  real(kind=db), dimension(nelor):: Elor, betalor

  do ie = 1,nelor
    E = Elor(ie) - E_cut
    Ep = E
    if( E > 0._db .and. Ep > 0._db .and. ( A > 0._db .or. Gamma_max > 0._db ) ) then
      lambda = 0._db
      if( A > 0._db ) lambda = lambda + 1 / ( A * sqrt(Ep) )
      if( Gamma_max > 0._db ) lambda = lambda + sqrt(Ep) / Gamma_max
      betalor(ie) = sqrt(E) / lambda
    else
      betalor(ie) = 0._db
    endif
  end do

  return
end

!***********************************************************************

subroutine gammarc(Ecent,Elarg,Gamma_max,E_cut,nelor, Elor,betalor)

  use declarations
  implicit none

  integer ie, nelor

  real(kind=db):: E, Ec, Ecent, E_cut, El, Elarg, Gamma_max, p
  real(kind=db), dimension(nelor):: Elor, betalor

  Ec = max( Ecent, 1.E-10_db )
  El = max( Elarg, 1.E-10_db )
  p = ( pi / 3 ) * Gamma_max / El

  do ie = 1,nelor
    E = Elor(ie) - E_cut
    if ( E <= 0._db ) then
      betalor(ie) = 0._db
    else
      betalor(ie) = Gamma_max * ( 0.5 + atan( p*(E/Ec - (Ec/E)**2)) / pi)
    endif
  end do

  return
end

!***********************************************************************

! Convolution by a gaussian

subroutine Gaussi(deltar,E_cut,Energ,nenerg,vibration,Y)

  use declarations
  implicit none

! nj : Points au dela de la gamme pour diminuer les effets de bord de la convolution
  integer, parameter:: nj = 10
  integer, parameter:: n1m = 1 - nj
  
  integer:: i, ie, je, n, ne1, ne2, nenerg
  
  real(kind=db):: b, def, deltar, E, E_cut, Efac, fac, p, pas, pdt, vib, vibration, Yint
  real(kind=db), dimension(nenerg):: Energ, Y
  real(kind=db), dimension(n1m:nenerg+nj):: de, Ef, Gaus(n1m:nenerg+nj), Xa(n1m:nenerg+nj)

  Ef(1:nenerg) = Energ(1:nenerg)
  Xa(1:nenerg) = Y(1:nenerg)

! Creation des points au dela de la gamme
  ne1 = 1 - nj
!      ne2 = nenerg + nj
  def = Ef(2) - Ef(1)
  do ie = 0,ne1,-1
    Ef(ie) = Ef(ie+1) - def
    Xa(ie) = Xa(1)
  end do
  ne2 = nenerg

  de(ne1) = Ef(ne1+1) - Ef(ne1)
  do ie = ne1+1,ne2-1
    de(ie) = 0.5 * ( Ef(ie+1) - Ef(ie-1) )
  end do
  de(ne2) = ( Ef(ne2) - Ef(ne2-1) )

  if( abs(deltar) < 1.e-10_db .and. abs(vibration) < 1.e-10_db ) return

  do ie = 1,nenerg

    vib = 2 * vibration * ( Ef(ie) - E_cut + 0.5_db )
    vib = max( 0._db, vib )
    b = deltar + vib
    if( abs(b) < 1.e-10_db ) then
      Y(ie) = Xa(ie)
      cycle
    endif

    gaus(:) = 0._db
    Pdt = 0._db
    do je = ne1,ne2
      n = max( int( 10 * de(je) / b ), 1 )
      pas = de(je) / ( n + 1 )
      if ( je == ne1 ) then
        E = Ef(je) - 0.5 * ( Ef(je+1) -  Ef(je) )
      else
        E = Ef(je) - 0.5 * ( Ef(je) -  Ef(je-1) )
      endif
      do i = 1,n
        E = E + pas
        if( ( E < Ef(je) .and. je /= ne1 ) .or. je == ne2 ) then
          p = ( E - Ef(je-1) ) / ( Ef(je) - Ef(je-1) )
          Yint = (1 - p ) * Xa(je-1) + p * Xa(je)
        else
          p = ( E - Ef(je) ) / ( Ef(je+1) - Ef(je) )
          Yint = (1 - p ) * Xa(je) + p * Xa(je+1)
        endif
        fac = -0.5 * ( ( E - Ef(ie) ) / b )**2
        if( fac > -600._db ) then
          efac = exp( fac )
          gaus(je) = gaus(je) + efac * Yint
          Pdt = Pdt + efac * de(je) / n
        endif
      end do
      gaus(je) = ( gaus(je) / n ) * de(je)
    end do
    Y(ie) = sum( gaus(ne1:ne2) ) / Pdt

  end do

  return
end

!***********************************************************************

! Calculation of the coefficients of the lorentzian

subroutine cflor(bb,betalor,E_cut,Eph,Elor,Energ,ie1,ie2,nef,nelor,nenerg,nenerge,Photoemission)

  use declarations
  implicit none

  integer:: ie, ie1, ie2, j, nef, nelor, nenerg, nenerge
  
  logical:: Photoemission

  real(kind=db):: def, E_cut, p
  real(kind=db), dimension(nenerge):: bb, Eph
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(nelor):: betalor, Elor

  Eph(1:nenerg) = Energ(1:nenerg)

! Creation des points au dela de la gamme
  def = Eph(2) - Eph(1)
  def = Eph(nenerg) - Eph(nenerg-1)
  do ie = nenerg+1,nenerge
    Eph(ie) = Eph(ie-1) + def
  end do

! Les etats en dessous de Fermi sont occupes
  do ie = 1,nenerg - 1
    if( Eph(ie) > E_cut - 1.e-10_db ) exit
  end do
  nef = ie

  if( Photoemission ) nef = max(1, nef-1)

  if( Photoemission ) then
    ie1 = 1
    ie2 = nef
  else
    ie1 = nef
    ie2 = nenerge
  endif

  do ie = 1,nenerge
    if( Eph(ie) <= Elor(1) ) then
      bb(ie) = betalor(1)
    elseif( Eph(ie) >= Elor(nelor) ) then
      bb(ie) = betalor(nelor)
    else
      do j = 2,nelor
        if( Elor(j) >= Eph(ie) ) exit
      end do
      p = ( Eph(ie) - Elor(j-1) ) / ( Elor(j) - Elor(j-1) )
      bb(ie) = p * betalor(j) + ( 1 - p ) * betalor(j-1)
    endif
    bb(ie) = 0.5_db * bb(ie) ! bb est en fait Gamma / 2
  end do

  return
end

!***********************************************************************

! Calculation of the extrapolation after the grid for the atomic scattering

subroutine Extrapat(bb_nenerg,dampl,Ephoton,Eseuil,Extrap,fpp0,fprime_atom,icheck,nenerg,numat)

  use declarations
  implicit none

  integer, parameter:: nim = 10000

  integer:: i, icheck, ie, nenerg, ni, numat
  
  character(len=2):: Chemical_symbol

  complex(kind=db):: dampl(nenerg)

  logical:: Extrap, fprime_atom

  real(kind=db):: alfa, bb_nenerg, bb2, dde, E1, Eph, Ephi, Ephm, Eseuil, f, fp, fpp, fpp0
  real(kind=db):: ei(nim), ei1(nim), ei2(nim), Ephoton(nenerg), fppn(nim)

  if( fprime_atom ) then
    Eph = max( Eseuil - 200 / rydb, 10 / rydb )
    dde = 0.5 /rydb
    Ephm = Eseuil + 5000. / rydb
    Ephi = Eseuil + 500. / rydb
    write(3,110) Chemical_symbol(numat)
    write(3,'(/A)') '   Energy_eV       f_prim      f_second'
    do i = 1,nim
      if( Eph < Ephi ) then
        Eph = Eph + dde
      else
        Eph = Eph + 10 * dde
      endif
      if( Eph > Ephm ) exit
      call fprime(numat,Eph,fpp,fp)
      write(3,'(3f13.6)') Eph*rydb, fp, fpp
    end do
    if( .not. Extrap ) return
  endif

  Eph = Eseuil - 2 / rydb
  call fprime(numat,Eph,fpp0,fp)

  Eph = Ephoton(nenerg)
  Ephm = Eseuil + 10000._db / rydb
  if( nenerg > 1 ) then
    dde = Ephoton(nenerg) - Ephoton(nenerg-1)
  else
    dde = 2 / rydb
  endif

  alfa = 1.02_db
  f = 0.5_db * alfa
  do i = 1,nim
    dde = alfa * dde
    Eph = Eph + dde
    if( Eph > ephm ) exit
    call fprime(numat,Eph,fpp,fp)
    Ei(i) = Eph
    Ei2(i) = Ei(i) + f * dde
    Ei1(i) = Ei(i) - 0.5_db * dde
    fppn(i) = - ( fpp - fpp0 ) / pi
  end do
  ni = min(i-1,nim)

  bb2 = bb_nenerg**2

  do ie = 1,nenerg
    E1 = Ephoton(ie)
    if( abs(bb_nenerg) < 1.e-10_db ) then
      dampl(ie) = sum( fppn(1:ni) * log( ( Ei2(1:ni) - E1 ) / (Ei1(1:ni) - E1 ) ) )
    else
      dampl(ie) = sum( fppn(1:ni) * ( 0.5*log( ( ( Ei2(1:ni) - E1 )**2 + bb2 ) / ( ( Ei1(1:ni) - E1 )**2 + bb2 ) ) &
              - img * ( atan( ( Ei2(1:ni) - E1 ) / bb_nenerg ) - atan( ( Ei1(1:ni) - E1 ) / bb_nenerg ) ) ) )
    endif
  end do

  if( icheck > 1 ) then
    write(3,'(/A)') '  Energy      dampl_r     dampl_i'
    do ie = 1,nenerg
      write(3,'(f10.3,1p,2e13.5)') Ephoton(ie)*rydb, dampl(ie)
    end do
  endif

  return
  110 format(/1x,a2,' atomic scattering scattering amplitude')
end

!***********************************************************************

subroutine extract_word(Length_line,motl,word,n_word,n_dim)

  use declarations
  
  implicit none
  
  integer:: i, i_word, ii, j, Length_line, n_dim, n_word

  character(len=Length_line):: motl
  character(len=length_word):: mot
  character(len=length_word), dimension(n_dim):: word

  i_word = 0
  i = 0
  do ii = 1,Length_line
    i = i + 1
    if( i > Length_line ) exit
    if( motl(i:i) == ' ' ) cycle
    i_word = i_word + 1
    mot = ' '
    do j = i,Length_line
      if( motl(j:j) == ' ' ) exit
      mot(j-i+1:j-i+1) = motl(j:j)
    end do
    i = j
    word(i_word) = mot
    if( i_word == n_word ) exit
  end do

  return
end

!***********************************************************************

! Calculation of transmission

subroutine Cal_MuTrans(icheck,Mu_mat_comp,Mus_mat,n_stokes,n_mat_pol,nes,Sample_thickness,Stokes_param)

  use declarations
  implicit none

  integer:: i, icheck, indp, ie, ipl, nes, n_mat_pol, n_stokes

  complex(kind=db):: Ch, diag, exp_TS, exp_TS_i, mu_pp, mu_ps, mu_sp, mu_ss, Sh, Tau, Trace

  complex(kind=db), dimension(2,2):: AM, Mu, Mat, Mat_pol
  complex(kind=db), dimension(nes,4,n_mat_pol):: Mu_mat_comp

  real(kind=db):: Cos_eta, Cos_2Bra, P1, P2, P3, Sample_thickness, Sin_eta

  real(kind=db), dimension(2,2):: Analyseur
  real(kind=db), dimension(5,n_stokes):: Stokes_param
  real(kind=db),  dimension(nes,n_stokes,n_mat_pol):: Mus_mat

  if( icheck > 1 ) write(3,100)
  
  do ipl = 1,n_mat_pol

    do indp = 1,n_stokes   ! Boucle sur la polarisation

! Parametres de Stokes
      P1 = Stokes_param(1,indp)
      P2 = Stokes_param(2,indp)
      P3 = Stokes_param(3,indp)

! Matrice polarisation
      Mat_pol(1,1) = 1._db + P1;  Mat_pol(1,2) = P2 - img*P3
      Mat_pol(2,1) = P2 + img*P3; Mat_pol(2,2) = 1._db - P1
      Mat_pol = Mat_pol / 2

! Tout est avec le complexe conjugue (Bragg, F, ..)
      Mat_pol = Conjg( Mat_pol)

! Angle de rotation de l'analyseur
      Cos_eta = Cos( Stokes_param(4,indp) )
      Sin_eta = Sin( Stokes_param(4,indp) )
! Angle de Bragg de l'analyseur
      Cos_2Bra = Cos( 2 * Stokes_param(5,indp) )

! Matrice Analyseur
      Analyseur(1,1) = Cos_eta
      Analyseur(1,2) = - Sin_eta
      Analyseur(2,1) = Cos_2Bra * Sin_eta
      Analyseur(2,2) = Cos_2Bra * Cos_eta

      if( icheck > 1 ) then
        write(3,110) ipl, indp, P1, P2, P3, Stokes_param(4:5,indp) * 180 / pi
        do i = 1,2
          write(3,120) Mat_pol(i,:), Analyseur(i,:)
        end do 
      endif
      
      do ie = 1,nes

        mu_ss = Mu_mat_comp(ie,1,ipl) 
        mu_sp = Mu_mat_comp(ie,2,ipl) 
        mu_ps = Mu_mat_comp(ie,3,ipl) 
        mu_pp = Mu_mat_comp(ie,4,ipl) 
        if( abs( mu_pp - mu_ss ) < 1.e-20_db .and. abs( mu_sp * mu_ps ) < 1.e-20_db ) then
          Mu(1,1) = 1._db
          Mu(2,1) = 0._db
          Mu(1,2) = 0._db
          Mu(2,2) = 1._db
        else 
          Tau = 0.25_db * sqrt( (mu_pp - mu_ss)**2 + 4 * mu_sp * mu_ps )
          exp_TS = exp( Tau * Sample_thickness )
          exp_TS_i = 1 / exp_TS
          Ch  = 0.5_db * ( Exp_TS + Exp_Ts_i )  
          Sh = 0.5_db * ( Exp_TS - Exp_Ts_i ) 
          diag = 0.25_db * Sh * ( mu_pp - mu_ss ) / Tau
          Mu(1,1) = ch + diag
          Mu(2,1) = - 0.5_db * Sh * mu_ps / Tau
          Mu(1,2) = - 0.5_db * Sh * mu_sp / Tau
          Mu(2,2) = ch - diag
        endif
        
        Mu(:,:) =  exp( -0.25_db * ( mu_pp + mu_ss ) * Sample_thickness ) * Mu(:,:)

        if( icheck > 1 .and. ie == 1 ) then
          write(3,'(/A)') '   Mu Matrix'
          do i = 1,2
            write(3,130) Mu(i,:)
          end do 
        endif
      
        AM = Matmul( Analyseur, Mu )

        if( icheck > 1 .and. ie == 1 ) then
          write(3,'(/A)') '   A.Mu Matrix'
          do i = 1,2
            write(3,130) AM(i,:)
          end do 
        endif

        Mat = Transpose( Conjg( AM ) )
        Mat = Matmul( AM, Matmul( Mat_pol, Mat ) )
        Trace = real( Mat(1,1) + Mat(2,2), db )

        if( icheck > 1 .and. ie == 1 ) then
          write(3,'(/A)') '   A.Mu.P.Tr(Conj(A.Mu)) Matrix'
          do i = 1,2
            write(3,130) Mat(i,:)
          end do 
        endif

        Mus_mat(ie,indp,ipl) = Trace
!        Mus_mat(ie,indp,ipl) = - log( Trace ) / Sample_thickness

      end do
    end do
  end do

  return
 100 format(/' ---- Transmission calculation and polarization rotation --------',60('-')) 
 110 format(/' Wave vector direction',i2,', Stoke parameters',i3, / &
             ' P1, P2, P3 =',3f9.5,', Analyzer angles =',2f11.5, // & 
             '                  Polarization matrix                    Analyzer Matrix')
 120 format(5x,2(1x,2f10.7),5x,2f11.7)
 130 format(1p,5x,2(1x,2e13.5))
end

!***********************************************************************
! Calcul de l'intensite en fonction des parametres de Stokes

subroutine Cal_stokes(As,Icirc,nes,n_stokes,np_stokes, nphi,nphim,npldafs,Stokes_param)

  use declarations
  implicit none

  integer:: ie, indp, ip, ipl, is, nes, n_stokes, np_stokes, nphim, npldafs
  integer, dimension(npldafs):: nphi

  complex(kind = db):: fps, fpp, fsp, fss, Trace

  complex(kind=db), dimension(2,2):: AF, F, Mat, Mat_pol
  complex(kind=db), dimension(nes,nphim,npldafs):: As

  real(kind=db):: Cos_eta, Cos_2Bra, P1, P2, P3, Sin_eta

  real(kind=db), dimension(2,2):: Analyseur
  real(kind=db), dimension(5,n_stokes):: Stokes_param
  real(kind=db), dimension(nes,nphim,np_stokes):: Icirc

  is = 0

  do ipl = 1,npldafs,4

    do indp = 1,n_stokes   ! Boucle sur la polarisation

      is = is + 1

! Parametres de Stokes
      P1 = Stokes_param(1,indp)
      P2 = Stokes_param(2,indp)
      P3 = Stokes_param(3,indp)

! Matrice polarisation
      Mat_pol(1,1) = 1._db + P1;  Mat_pol(1,2) = P2 - img*P3
      Mat_pol(2,1) = P2 + img*P3; Mat_pol(2,2) = 1._db - P1
      Mat_pol = Mat_pol / 2

! Tout est avec le complexe conjugue (Bragg, F, ..)
      Mat_pol = Conjg( Mat_pol)

! Angle de rotation de l'analyseur
      Cos_eta = Cos( Stokes_param(4,indp) )
      Sin_eta = Sin( Stokes_param(4,indp) )
! Angle de Bragg de l'analyseur
      Cos_2Bra = Cos( 2 * Stokes_param(5,indp) )

! Matrice Analyseur
      Analyseur(1,1) = Cos_eta
      Analyseur(1,2) = - Sin_eta
      Analyseur(2,1) = Cos_2Bra * Sin_eta
      Analyseur(2,2) = Cos_2Bra * Cos_eta

      do ip = 1,nphi(ipl)

        do ie = 1,nes

          fss = As(ie,ip,ipl)
          fsp = As(ie,ip,ipl+1)
          fps = As(ie,ip,ipl+2)
          fpp = As(ie,ip,ipl+3)

          F(1,1) = fss; F(1,2) = fps
          F(2,1) = fsp; F(2,2) = fpp

          AF = Matmul( Analyseur, F )

          Mat = Transpose( Conjg( AF ) )
          Mat = Matmul( AF, Matmul( Mat_pol, Mat ) )
          Trace = Mat(1,1) + Mat(2,2)

          Icirc(ie,ip,is) = Real ( Trace , db )

        end do

      end do
    end do
  end do

  return
end

!***********************************************************************

! Calcul de l'intensite corrige d'absorption en tenant compte de la birefringence

subroutine Corr_abs(angxyz,As,axyz,c_micro,d_dead,Double_cor,Eseuil,fpp_avantseuil,hkl_dafs,hkl_S,Icor,Icirccor, &
          Icircdcor,Idcor,mus,n_stokes,nes,np_stokes,nphi,nphim, npldafs,Self_abs,Stokes_param)

  use declarations
  implicit none

  integer:: i, ie, indp, ip, ipl, ipl_pas, is, j, k, n_stokes, nes, np_stokes, nphim, npldafs
  integer, dimension(3):: hkl, hkl_S
  integer, dimension(npldafs):: nphi
  integer, dimension(3,npldafs):: hkl_dafs

  logical:: Check, Double_cor, Self_abs

  complex(kind = db):: dmu_i, dmu_s, Expo, fps, fpp, fsp, fss, mu_i, mu_pp_i, mu_pp_s, mu_ps_i, mu_ps_s, mu_s, &
       mu_sp_i, mu_sp_s, mu_ss_i, mu_ss_s, mups_i, mups_s, musp_i, musp_s, t_i, t_s, tau_i, tau_s, Trace

  complex(kind=db), dimension(2):: Ex_i, Ex_s
  complex(kind=db), dimension(4):: Ex
  complex(kind=db), dimension(2,2):: F, M1, M2, M3, Mat_pol
  complex(kind=db), dimension(2,2,2):: Tr_i, Tr_s
  complex(kind=db), dimension(2,2,4):: ATFT, TFT
  complex(kind=db), dimension(nes,nphim,npldafs):: As
  complex(kind=db), dimension(nes,nphim,npldafs,2):: mus

  real(kind=db):: alfa_S, Attenuation, Brag_A, c_micro, corr_ref, Correction, Cos_eta, Cos_2Bra, csi, css, d_dead, Eseuil, eta, &
    fpp_avantseuil, mu_i_m, mu_s_m, mum, mum_0, P1, P1A, P2, P3, Resultat, Sin_eta, theta_B

  real(kind=db), dimension(3):: angxyz, axyz
  real(kind=db), dimension(2,2):: Analyseur
  real(kind=db), dimension(5,n_stokes):: Stokes_param
  real(kind=db), dimension(nes,nphim,npldafs):: Icor, Idcor
  real(kind=db), dimension(nes,nphim,np_stokes):: Icirccor, Icircdcor

  Check = .false.

  if( Self_abs ) then
    ipl_pas = 1
  else
    ipl_pas = 4
  endif

  do ipl = 1,npldafs,ipl_pas

    hkl(:) = hkl_dafs(:,ipl)
    call Angle_Bragg(alfa_S,angxyz,axyz,Eseuil,hkl,hkl_S,theta_B)

    csi = 1 / sin( theta_B + alfa_S )
    css = 1 / sin( theta_B - alfa_S )
    mum_0 = c_micro * fpp_avantseuil * ( csi + css )
    corr_ref = mum_0 / exp( - mum_0 * d_dead )

    do ip = 1,nphi(ipl)

      do ie = 1,nes

        if( Self_abs .or. Double_cor ) mum = real( mus(ie,ip,ipl,1), db )*csi + real( mus(ie,ip,ipl,2), db )*css
        if( Double_cor ) Correction = mum_0 / mum

        if( self_abs ) then
          Attenuation = exp( - mum * d_dead ) * corr_ref
          Icor(ie,ip,ipl) = ( Attenuation / mum ) * abs( As(ie,ip,ipl) )**2
          if( Double_cor ) Idcor(ie,ip,ipl) = Icor(ie,ip,ipl) * Correction
          cycle
        endif

        mu_ss_i = mus(ie,ip,ipl,1)
        mu_ss_s = mus(ie,ip,ipl,2)
        mu_sp_i = mus(ie,ip,ipl+1,1)
        mu_sp_s = mus(ie,ip,ipl+1,2)
        mu_ps_i = mus(ie,ip,ipl+2,1)
        mu_ps_s = mus(ie,ip,ipl+2,2)
        mu_pp_i = mus(ie,ip,ipl+3,1)
        mu_pp_s = mus(ie,ip,ipl+3,2)

        fss = As(ie,ip,ipl)
        fsp = As(ie,ip,ipl+1)
        fps = As(ie,ip,ipl+2)
        fpp = As(ie,ip,ipl+3)

        mu_i = ( mu_ss_i + mu_pp_i ) * csi / 4
        mu_s = ( mu_ss_s + mu_pp_s ) * css / 4

        mum = 0.5_db*real( mu_ss_i + mu_pp_i, db )*csi + 0.5_db*real( mu_ss_s + mu_pp_s, db )*css

        Attenuation = exp( - mum * d_dead ) * corr_ref

        tau_i = sqrt( ( mu_pp_i - mu_ss_i )**2 + 4*mu_sp_i*mu_ps_i )
        tau_s = sqrt( ( mu_pp_s - mu_ss_s )**2 + 4*mu_sp_s*mu_ps_s )

        t_i = csi * tau_i / 4
        t_s = css * tau_s / 4

        if( abs(tau_i) > eps10 ) then
          mups_i = 2*mu_ps_i / tau_i
          musp_i = 2*mu_sp_i / tau_i
          dmu_i = ( mu_pp_i - mu_ss_i ) / tau_i
        else
          mups_i = (0._db, 0._db)
          musp_i = (0._db, 0._db)
          dmu_i = (0._db, 0._db)
        endif

        if( abs(tau_s) > eps10 ) then
          mups_s = 2*mu_ps_s / tau_s
          musp_s = 2*mu_sp_s / tau_s
          dmu_s = ( mu_pp_s - mu_ss_s ) / tau_s
        else
          mups_s = (0._db, 0._db)
          musp_s = (0._db, 0._db)
          dmu_s = (0._db, 0._db)
        endif

! Matrices transmittance
! Elles sont decomposee en 2 parties. La matrice transmittance totale
! est donnee par T = 0.5*(exp(Tau*l)*Tr(1) + exp(-Tau*l)*Tr(2) )
        Tr_i(1,1,1) = 1 + dmu_i; Tr_i(1,2,1) = - musp_i
        Tr_i(2,1,1) = - mups_i;  Tr_i(2,2,1) = 1 - dmu_i

        Tr_i(1,1,2) = 1 - dmu_i; Tr_i(1,2,2) = musp_i
        Tr_i(2,1,2) = mups_i;    Tr_i(2,2,2) = 1 + dmu_i

        Tr_i(:,:,:) = Tr_i(:,:,:) / 2

        Tr_s(1,1,1) = 1 + dmu_s; Tr_s(1,2,1) = - musp_s
        Tr_s(2,1,1) = - mups_s;  Tr_s(2,2,1) = 1 - dmu_s

        Tr_s(1,1,2) = 1 - dmu_s; Tr_s(1,2,2) = musp_s
        Tr_s(2,1,2) = mups_s;    Tr_s(2,2,2) = 1 + dmu_s

        Tr_s(:,:,:) = Tr_s(:,:,:) / 2

! Exposant des exponentielles

        Ex_i(1) = - mu_i + t_i;  Ex_i(2) = - mu_i - t_i

        Ex_s(1) = - mu_s + t_s;  Ex_s(2) = - mu_s - t_s

! Amplitude de diffusion

        F(1,1) = fss; F(1,2) = fps
        F(2,1) = fsp; F(2,2) = fpp

! Produit
        k = 0
        do i = 1,2
          do j = 1,2
            k = k + 1
            M1(:,:) = Tr_s(:,:,i)
            M2(:,:) = Tr_i(:,:,j)
            M3 = Matmul( M1, Matmul( F, M2 ) )
            TFT(:,:,k) = M3
            Ex(k) = Ex_s(i) + Ex_i(j)
          end do
        end do

        if( check ) then
          Write(6,'(A)') ' tau_i tau_s'
          Write(6,'(1p,8e13.5)') tau_i, tau_s, t_i, t_s
          Write(6,'(A)') ' Tr_i'
          do i = 1,2
            Write(6,'(1p,8e13.5)') Tr_i(i,:,1), Tr_i(i,:,2)
          end do
          Write(6,'(A)') ' Tr_s'
          do i = 1,2
            Write(6,'(1p,8e13.5)') Tr_s(i,:,1), Tr_s(i,:,2)
          end do
          Write(6,'(A)') ' Ex_i'
            Write(6,'(1p,8e13.5)') Ex_i(:)
          Write(6,'(A)') ' Ex_s'
            Write(6,'(1p,8e13.5)') Ex_s(:)
          Write(6,'(A)') ' F'
          do i = 1,2
            Write(6,'(1p,8e13.5)') F(i,:)
          end do
          do k = 1,4
            Write(6,'(a9,i2)') ' TFT, k =', k
            Write(6,'(1p,8e13.5)') TFT(1,:,k), Ex(k)
            Write(6,'(1p,8e13.5)') TFT(2,:,k)
          end do
        endif

        do indp = 1,4+n_stokes  ! Boucle sur polarisation et analyse

          select case(indp)

            case(1)
              P1 = 1._db; P2 = 0._db; P3 = 0._db
              eta = 0._db; Brag_A = pi / 4

            case(2)
              P1 = 1._db; P2 = 0._db; P3 = 0._db
              eta = pi / 2; Brag_A = pi / 4

            case(3)
              P1 = -1._db; P2 = 0._db; P3 = 0._db
              eta = 0._db; Brag_A = pi / 4

            case(4)
              P1 = -1._db; P2 = 0._db; P3 = 0._db
              eta = pi / 2; Brag_A = pi / 4

            case default

              is = ipl + indp - 5

              eta = Stokes_param(4,indp-4)
              Brag_A = Stokes_param(5,indp-4)

! Parametres de Stokes
              P1 = Stokes_param(1,indp-4)
              P2 = Stokes_param(2,indp-4)
              P3 = Stokes_param(3,indp-4)

! Parametre de Stokes P1 equivalent pour l'analyseur
              P1A = Cos( 2 * eta ) * Sin( 2 * Brag_A )

          end select

! Matrice polarisation
          Mat_pol(1,1) = 1._db + P1;  Mat_pol(1,2) = P2 - img*P3
          Mat_pol(2,1) = P2 + img*P3; Mat_pol(2,2) = 1._db - P1
          Mat_pol = Mat_pol / 2

! Tout est avec le complexe conjugue (Bragg, F, ..)
          Mat_pol = Conjg( Mat_pol)

! Angle de rotation de l'analyseur
          Cos_eta = Cos( eta )
          Sin_eta = Sin( eta )
! Angle de Bragg de l'analyseur
          Cos_2Bra = Cos( 2 * Brag_A )

! Matrice Analyseur
          Analyseur(1,1) = Cos_eta
          Analyseur(1,2) = - Sin_eta
          Analyseur(2,1) = Cos_2Bra * Sin_eta
          Analyseur(2,2) = Cos_2Bra * Cos_eta

          if( Check ) then
            Write(6,'(/A)') ' Analyseur'
            Write(6,'(2f8.5)') Analyseur(1,:)
            Write(6,'(2f8.5)') Analyseur(2,:)
          endif

          do i = 1,4
            M1(:,:) = TFT(:,:,i)
            M2 = Matmul( Analyseur, M1 )
            ATFT(:,:,i) = M2
          end do

          Resultat = 0._db
          do i = 1,4
            do j = 1,4
              M1(:,:) = ATFT(:,:,i)
              M2(:,:) = ATFT(:,:,j)
              M2 = Transpose( Conjg( M2 ) )
              M3 = Matmul( M1, Matmul( Mat_pol, M2 ) )
              Trace = M3(1,1) + M3(2,2)
              Expo = Ex(i) + Conjg( Ex(j) )
! Signe moins car integrale de 0 a l'infini = - Primitive(0)
              Resultat = Resultat - Real ( Trace / Expo, db )
            end do
          end do
          Resultat = Resultat * Attenuation

          if( check ) then
            Write(6,'(/A)') ' Polar'
            Write(6,'(4f8.3)') Mat_pol(1,:)
            Write(6,'(4f8.3)') Mat_pol(2,:)
            Write(6,'(/A)') ' Attenuation, Resultat'
            Write(6,'(1p,5e13.5)') Attenuation, Resultat
          endif

          select case(indp)

            case(1)  ! sigma - sigma
              Icor(ie,ip,ipl) = Resultat
              if( Double_cor ) Idcor(ie,ip,ipl) = Resultat * Corr_ref / Real( mu_ss_i*csi + mu_ss_s*css, db)
            case(2)  ! pi - sigma
              Icor(ie,ip,ipl+2) = Resultat
              if( Double_cor ) Idcor(ie,ip,ipl+2) = Resultat * Corr_ref / Real( mu_pp_i*csi + mu_ss_s*css, db)
            case(3)  ! sigma - pi
              Icor(ie,ip,ipl+1) = Resultat
              if( Double_cor ) Idcor(ie,ip,ipl+1) = Resultat * Corr_ref / Real( mu_ss_i*csi + mu_pp_s*css, db)
            case(4)  ! pi - pi
              Icor(ie,ip,ipl+3) = Resultat
                if( Double_cor ) Idcor(ie,ip,ipl+3) = Resultat * Corr_ref / Real( mu_pp_i*csi + mu_pp_s*css, db)

            case default

              Icirccor(ie,ip,is) = Resultat

              if( Double_cor ) then
                mu_i_m = ( 1 + P1 ) * Real( mu_ss_i, db ) + ( 1 - P1 ) * Real( mu_pp_i, db )

                mu_s_m = ( 1 + P1A ) * Real( mu_ss_s, db ) + ( 1 - P1A ) * Real( mu_pp_s, db )

                Icircdcor(ie,ip,is) = Resultat * Corr_ref * 2 / ( mu_i_m * csi + mu_s_m * css )
              endif

          end select

        end do ! Boucle indp

      end do
    end do
  end do

  return
end

!***********************************************************************

! Effect of temperature using the Debye model

subroutine Debye_effect(Abs_U_iso,Abs_U_iso_inp,Energ,Energphot,Ephoton,Eseuil,fac,ifich,n_col,nenerg,nfich, &
                        nom_col,numat,nxan,Shift_U_iso,U_iso_inp,V0muf,Xanes)

  use declarations
  implicit none

  character(len=Length_word):: nomab
  character(len=Length_word), dimension(n_col):: nom_col

  integer:: i, ie, ifich, n_col, nenerg, nfich, nxan, Z
  integer, dimension(nfich):: numat
  
  logical:: Energphot, U_iso_inp

  real(kind=db):: Abs_U_iso_inp, c, Conv_Mbarn_nelec, E, Eph, f, fac, fp, fpp, fpp0, p, Shift_U_iso
  real(kind=db), dimension(nenerg):: Energ, Ephoton, Xanes_atom 
  real(kind=db), dimension(nenerg,nxan):: Xanes 
  real(kind=db), dimension(nfich):: Abs_U_iso, Eseuil, V0muf 

  Z = numat(ifich) 
  
!  Abs_U_iso = <u^2> is in A^2
  c = 0.04 * Rydb * ( m_electron * e_electron / hbar**2 )

  do i = 1,nxan
    nomab = nom_col(i)
    nomab = adjustl(nomab)
    if( nomab /= 'XANES_atom' ) cycle
    Xanes_atom(:) = Xanes(:,i)
    exit   
  end do

 ! Calculation of atomic spectra when not existing
  if( i > nxan ) then

    Eph = Eseuil(ifich) - 1
    Eph = max( Eph,  0.5_db )
    call fprime(Z,Eph,fpp0,fp)

    do ie = 1,nenerg
      Eph = Ephoton(ie) - Shift_U_iso
      call fprime(Z,Eph,fpp,fp)
      f = fac / conv_mbarn_nelec(Eph)
      Xanes_atom(ie) = f * ( fpp - fpp0 )
    end do
  
  endif  

  do ie = 1,nenerg
    if( Energphot ) then
      E = Energ(ie) - Eseuil(ifich) - V0muf(ifich)
    else
      E = Energ(ie) - V0muf(ifich)
    endif
    if( E < eps10 ) cycle

    if( U_iso_inp ) then 
      p = exp( - c * Abs_U_iso_inp * E )
    else
      p = exp( - c * Abs_U_iso(ifich) * E )
    endif
    do i = 1,nxan
      nomab = nom_col(i)
      nomab = adjustl(nomab)
      if( nomab(1:3) == 'dic' .or. nomab(1:3) == 'Dic' ) then
        Xanes(ie,i) = p * Xanes(ie,i)
      else
        Xanes(ie,i) = p * Xanes(ie,i) + (1 - p) * Xanes_atom(ie)
      endif
    end do
  end do
           
  return
end

!***********************************************************************

subroutine Debye_effect_a(Abs_U_iso,Abs_U_iso_inp,Adafs,dph,Energ,Energphot,Ephoton,Eseuil,fac,ifich, &
                        nenerg,nfich,Dafs_cal,nphim,npldafs,numat,Shift_U_iso,U_iso_inp,V0muf) 
  use declarations
  implicit none

  integer:: i, icheck, ie, ifich, nenerg, nfich, nphim, npldafs, Z
  integer, dimension(nfich):: numat

  complex(kind=db), dimension(nenerg,nphim,npldafs):: Adafs 
  complex(kind=db), dimension(npldafs):: dph 
  complex(kind=db), dimension(nenerg):: Dafs_atom 
  
  logical:: Dafs_cal, Energphot, U_iso_inp

  real(kind=db):: Abs_U_iso_inp, c, Conv_Mbarn_nelec, E, Eph, f, fac, fp, fpp, fpp0, p, Shift_U_iso
  real(kind=db), dimension(nenerg):: Energ, Ephoton 
  real(kind=db), dimension(nfich):: Abs_U_iso, Eseuil, V0muf 

  Z = numat(ifich) 
  
!  Abs_U_iso = <u^2> is in A^2
  c = 0.04 * Rydb * ( m_electron * e_electron / hbar**2 )

 ! Calculation of atomic spectra when not existing

  Eph = Eseuil(ifich) - 1
  Eph = max( Eph,  0.5_db )
  call fprime(Z,Eph,fpp0,fp)

  f = fac / pi
  
  do ie = 1,nenerg
    Eph = Ephoton(ie) - Shift_U_iso
    call fprime(Z,Eph,fpp,fp)
    if( .not. Dafs_cal ) f = fac / ( pi * Conv_Mbarn_nelec(Eph) ) 
    Dafs_atom(ie) = f * cmplx( fpp - fpp0, 0._db, db )
  end do

  icheck = 1
  if( icheck > 1 ) then
    open(99, file = 'Debye_out.txt')
    write(99,'(A)') '   Energy     Dafs_r       Dafs_i      Dafs_a_r    Dafs_a_i'
    do ie = 1,nenerg
      write(99,'(f10.3,1p,4e13.5)') Energ(ie)*Rydb, Adafs(ie,1,5), Dafs_atom(ie)*conjg(dph(5))
    end do
    stop
  endif

  do ie = 1,nenerg
    if( Energphot ) then
      E = Energ(ie) - Eseuil(ifich) - V0muf(ifich)
    else
      E = Energ(ie) - V0muf(ifich)
    endif
    if( E < eps10 ) cycle
    
    if( U_iso_inp ) then 
      p = exp( - c * Abs_U_iso_inp * E )
    else
      p = exp( - c * Abs_U_iso(ifich) * E )
    endif
    do i = 1,npldafs
      Adafs(ie,:,i) = p * Adafs(ie,:,i) + (1 - p) * conjg( dph(i) ) * Dafs_atom(ie)
    end do
  end do
   
  return
end

!***********************************************************************

subroutine Angle_Bragg(alfa_S,angxyz,axyz,Eseuil,hkl,hkl_S, theta_B)

  use declarations
  implicit none

  integer:: i, j
  integer, dimension(3):: hkl, hkl_S

  real(kind=db):: alfa_S, detmat, dhkl, Eseuil, fac, konde, qkn, qkn_S, theta_B, vol
  real(kind=db), dimension(3):: angxyz, axyz, cosdir, hklred, qk, qk_S, vx, vy, vz, wx, wy, wz
  real(kind=db), dimension(0:3):: det
  real(kind=db), dimension(3,3):: Cubmat, mat

  character(len=5):: Struct

  call cal_cubmat(angxyz,Cubmat,struct)

  konde = 0.5_db * alfa_sf * Eseuil

  cosdir(:) = cos( angxyz(:) )

! Base du reseau direct, exprimee dans une base orthonormee
  vx(:) = Cubmat(:,1) * axyz(:)
  vy(:) = Cubmat(:,2) * axyz(:)
  vz(:) = Cubmat(:,3) * axyz(:)

! wx, wy, wz : base du reseau reciproque exprimee dans la meme base
  call prodvec(wx,vy,vz)

  vol = sum( wx(:) * vx(:) )
  wx(:) = wx(:) / vol
  call prodvec(wy,vz,vx)
  wy(:) = wy(:) / vol
  call prodvec(wz,vx,vy)
  wz(:) = wz(:) / vol

  qk_S(:) = hkl_S(1) * wx(:) + hkl_S(2) * wy(:) + hkl_S(3) * wz(:)
  qkn_S = sqrt( sum( qk_S(:)**2 ) )

  if( abs(qkn_S) > eps10 ) qk_S(:) = qk_S(:) / qkn

  if( hkl(1) == 0 .and. hkl(2) == 0 .and. hkl(3) == 0 ) then

    dhkl = 0._db
    theta_B = 0._db

  else

    hklred(:) = hkl(:) / axyz(:)
    do i = 0,3
      do j = 1,3
        mat(j,j) = 1._db
      end do
      mat(1,2) = cosdir(3)
      mat(1,3) = cosdir(2)
      mat(2,3) = cosdir(1)
      mat(2,1) = mat(1,2)
      mat(3,1) = mat(1,3)
      mat(3,2) = mat(2,3)
      if( i > 0 ) mat(:,i) = hklred(:)
      det(i) = detmat(mat)
    end do
! Distance interplan
    dhkl = sqrt( det(0) / sum( hklred(1:3) * det(1:3) ) )
    fac = abs( pi / ( konde * dhkl ) )

    if( fac > 1._db ) then
      call write_error
      write(6,120) hkl(:)
      stop
    endif

    theta_B = asin( pi / ( konde * dhkl ) )

  endif

! Calcul angle avec plan de surface

  if( qkn_S > eps10 ) then

    qk(:) = hkl(1) * wx(:) + hkl(2) * wy(:) + hkl(3) * wz(:)
    qkn = sqrt( sum( qk(:)**2 ) )
    if( abs(qkn) > eps10 ) then
      qk(:) = qk(:) / qkn
      alfa_S = acos( sum( qk(:)*qk_S(:) ) )
    else
      alfa_S = 0._db
    endif
  else
    alfa_S = 0._db
  endif

  return
  120 format(//' The reflection number',i3,' : (h,k,l) = (',3i3,') does not exist at this energy !'//)
end

!*******************************************************************************************************************

subroutine Write_transpose(Convolution_out,Energ_tr,Es,n_col,n_energ_tr,n_signal,nes,nom_col,npldafs,Signal)

  use declarations
  implicit none
  
  integer:: i, i_col, i_cor, i_hk, i1, i2, ie, index, index_hk, index_max, ipas, ipr, j, l, l2, l3, le, length_hk, ll, n_col, &
            n_cor, n_hk_tr, n_energ_tr, n_signal, nes, npldafs
  integer, dimension(:), allocatable:: hk_length, n_index_hk
  integer, dimension(n_energ_tr):: index_ie

  logical:: Double_cor
  
  character(len=6):: mot6, mot6_b
  character(len=6), dimension(:), allocatable:: hk
  character(len=15):: mot15, mot15_b, mot15_c, mot15_d
  character(len=15), dimension(n_col):: nom_col
  character(len=15), dimension(:,:,:), allocatable:: E_string, l_value
  character(len=15), dimension(:,:,:,:), allocatable:: Signal_out
  character(len=Length_name):: Convolution_out, Convolution_tr

  real(kind=db):: E
  real(kind=db), dimension(nes):: Es
  real(kind=db), dimension(n_energ_tr):: Energ_tr, p1
  real(kind=db), dimension(npldafs,n_signal):: Signal
 
  Convolution_tr = Convolution_out
  l = len_trim(Convolution_tr)
  Convolution_tr(l-3:l+3) = '_tr.txt'
  
  ipr = 4
  
  open( ipr, file = Convolution_tr )

  if( n_signal == 3*nes ) then
    Double_cor = .true.
    i1 = n_col - 5*npldafs + 1
    ipas = 5
    n_cor = 3
  elseif( n_signal == 2*nes ) then
    Double_cor = .false.
    i1 = n_col - 4*npldafs + 1
    ipas = 4
    n_cor = 2
  else
    Double_cor = .false.
    i1 = n_col - npldafs + 1
    ipas = 1
    n_cor = 1
  endif

  do ie = 1,n_energ_tr
    do i = 2,nes
      if( Es(i) > Energ_tr(ie) - eps10 ) exit
    end do
    i = min( i, nes )
    index_ie(ie) = i 
    p1(ie) = ( Energ_tr(ie) - Es(i-1) ) / ( Es(i) - Es(i-1) )
  end do

  do j = 1,2
  
    i_hk = 0
    index_hk = 0
    mot6_b = ' '
    do i_col = i1,n_col,ipas
      mot15 = ' '
      mot15 = nom_col(i_col)
      do i = 1,15
        if( mot15(i:i) /= '(' ) cycle
        mot6 = ' '
        if( mot15(i+1:i+1) == '-' .and. mot15(i+3:i+3) == '-' ) then
          mot6(1:4) = mot15(i+1:i+4)
        elseif( mot15(i+1:i+1) == '-' .or. mot15(i+2:i+2) == '-' ) then
          mot6(1:3) = mot15(i+1:i+3)
        else
          mot6(1:2) = mot15(i+1:i+2)  
        endif
        length_hk = len_trim( mot6 )
        exit
      end do

      do i = 1,14
        if( mot15(i:i) /= ')' ) cycle
        ll = len_trim( mot6 )
        if( mot15(i+1:i+1) /= ' ' ) mot6(ll+1:ll+1) = mot15(i+1:i+1)
        if( i >= 14 ) exit
        if( mot15(i+2:i+2) /= ' ' ) mot6(ll+2:ll+2) = mot15(i+2:i+2)
        exit
      end do

      if( mot6 == mot6_b ) then    
        index_hk = index_hk + 1
        if( j == 2 ) n_index_hk(i_hk) = index_hk 
      else  
        i_hk = i_hk + 1
        index_hk = 1
        if( j == 2 ) then
          hk(i_hk) = mot6
          hk_length(i_hk) = length_hk
        endif
        mot6_b = mot6 
      endif
    end do
    
    if( j == 1 ) then
      n_hk_tr = i_hk
      allocate( n_index_hk(n_hk_tr) )
      n_index_hk(:) = 1
      allocate( hk(n_hk_tr) )
      allocate( hk_length(n_hk_tr) )
    endif
    
  end do
  
  index_max = 0  
  do i_hk = 1,n_hk_tr
    index_max = max( index_max, n_index_hk(i_hk) )
  end do
  
  allocate( l_value(index_max,n_cor,n_hk_tr) )
  allocate( Signal_out(index_max,n_energ_tr,n_cor,n_hk_tr) )
  allocate( E_string(n_energ_tr,n_cor,n_hk_tr) )
  Signal_out(:,:,:,:) = ' '
  l_value(:,:,:) = ' '

  index = 0
  index_hk = 0
  i_hk = 1

  do i_col = i1,n_col,ipas
    index = index + 1
    index_hk = index_hk + 1
    mot15 = ' '
    mot15 = nom_col(i_col)
    i2 = 16
    do i = 1,15
      if( mot15(i:i) == '(' ) then
        i1 = i + hk_length(i_hk) + 1
      elseif( mot15(i:i) == ')' ) then
        i2 = i-1
        exit
      endif
    end do
    mot15_b = ' '
    mot15_b(1:i2-i1+1) = mot15(i1:i2)
    l_value(index_hk,:,i_hk) = adjustr( mot15_b )

    do i_cor = 1,n_cor
      do ie = 1,n_energ_tr
        j = ( i_cor - 1 ) * nes + index_ie(ie)
        mot15_b = ' '
        write(mot15_b,'(1p,e15.7)') ( 1 - p1(ie) ) * Signal(index,j-1) + p1(ie) * Signal(index,j)
        Signal_out(index_hk,ie,i_cor,i_hk) = mot15_b
      end do
    end do 
    
    if( index_hk == n_index_hk(i_hk) ) then
      index_hk = 0
      i_hk = i_hk + 1
      if( i_hk > n_hk_tr ) exit
    endif   
  end do 

  l = 1
  mot15_b = ''
  mot15_b(1:1) = 'I'
  do i_cor = 1,n_cor
    if( i_cor == 2 ) then
      mot15_b(2:2) = 'c'
      l = 2
    elseif( i_cor == 3 ) then
      mot15_b(2:2) = 'd'
      l = 2
    endif
    do i_hk = 1,n_hk_tr
      mot15_d = mot15_b
      mot6 = hk(i_hk)
      le = len_trim( mot6 ) 
      mot15_d(l+1:l+1) = '_'
      mot15_d(l+2:l+1+le) = mot6(1:le) 
      l2 = len_trim(mot15_d) + 1
      mot15_d(l2:l2) = '('
      do i = 1,n_energ_tr
        E = Energ_tr(i)*rydb
        mot15_c = mot15_d
        mot15 = ' '
        if( abs( nint( E ) - E ) < eps10 ) then
          write(mot15,'(i6)') nint( E )
        elseif( abs( nint( 10*E ) - 10*E ) < eps10 ) then
          write(mot15,'(f15.1)') Energ_tr(i)*rydb
        elseif( abs( nint( 100*E ) - 100*E ) < eps10 ) then
          write(mot15,'(f15.2)') Energ_tr(i)*rydb
        elseif( abs( nint( 1000*E ) - 1000*E ) < eps10 ) then
          write(mot15,'(f15.3)') Energ_tr(i)*rydb
        elseif( abs( nint( 1000*E ) - 1000*E ) < eps10 ) then
          write(mot15,'(f15.4)') Energ_tr(i)*rydb
        endif
        mot15 = adjustl(mot15)
        l3 = len_trim(mot15)
        mot15_c(l2+1:l2+l3) = mot15(1:l3)
        if( l2+l3 < 14 ) mot15_c(l2+l3+1:l2+l3+1) = ')'
        mot15 = adjustr( mot15_c )
        E_string(i,i_cor,i_hk) = mot15
      end do
    end do

  end do
  
  write(ipr,110) (( '         l_'//hk(i_hk), E_string(:,i_cor,i_hk), i_cor = 1,n_cor ), i_hk = 1,n_hk_tr )

  do index = 1,index_max
    write(ipr,110) ( ( l_value(index,i_cor,i_hk), Signal_out(index,:,i_cor,i_hk), i_cor = 1,n_cor ), i_hk = 1,n_hk_tr )
  end do 
  
  Close(ipr)
  
  deallocate( E_string, hk, hk_length, l_value, n_index_hk, Signal_out )

  return
  110 format(10000a15)
end

!***********************************************************************

! Convolution by a gaussian

subroutine Main_gaussian(File_in,itape,File_out)

  use declarations
  implicit none

  integer:: eof, i, ier, ipr, istat, itape, l, ligne, n, n_text, ncol, nnombre, nx

  character(len=9):: keyword
  character(len=132):: identmot, mot, Title, Traduction
  character(len=Length_name):: File_in, File_out

  real(kind=db):: a, sigma
  real(kind=db), dimension(:), allocatable:: x
  real(kind=db), dimension(:,:), allocatable:: y, z
 
  Rewind(itape)
  
  i = 0
  do ligne = 1,1000

    n = nnombre(itape,132)
    read(itape,'(A)',iostat=eof) mot
    if( eof /= 0 ) exit

    keyword = identmot(mot,9)
    keyword = Traduction(keyword)
    
    if( keyword == 'end' ) exit

    if( keyword == 'conv_gaus' ) then
      n = nnombre(itape,132)
      read(itape,*,iostat=ier) sigma
      if( ier > 0 ) call write_err_form(itape,keyword)
      i = i + 1
    endif
    
  end do
  
  if( i /= 1 .or. File_in == ' ' ) then
    call write_error
    do ipr = 6,9,3
      write(ipr,'(A)') ' Both "File_in" and "Conv_gaus" keywords must be in the indata file !'
    end do
    stop
  endif

  if( sigma < -eps10 ) then
    call write_error
    do ipr = 6,9,3
      write(ipr,110) sigma
    end do
    stop
  endif
  
! Dimensions

  Open(20, file = File_in, status='old', iostat=istat )
      
  do i = 1,10000
    n = nnombre(20,10000)
    if( n == 0 ) then
      read(20,'(A)') mot
    else
      exit
    endif
  end do
  n_text = i - 1
  
  ncol = n - 1

  do i = 1,1000000
    read(20,*,iostat=istat) a
    if( istat /= 0 ) exit
  end do
  nx = i - 1

  allocate( x(nx) ); allocate( y(nx,ncol) ); allocate( z(nx,ncol) )

! Reading

  Rewind(20)

  do i = 1, n_text - 1
    read(20,*)
  end do
  if( n_text > 0 ) read(20,'(A)') Title

  do i = 1,nx
    read(20,*,iostat=istat) x(i), y(i,:)
    if( istat /= 0 ) then
      call write_error
      do ipr = 6,9,3
        write(ipr,120) File_in
      end do
      stop
    endif
  end do

  Close(20) 
  
  call Conv_gaussian(ncol,nx,sigma,x,y,z)

! Writing
  
  l = len_trim( File_out )
  if( l > 4 ) then
    if( File_out(l-3:l-3) /= '.' ) File_out(l+1:l+4) = '.txt'
  endif

  Open(20, file = File_out )

  if( n_text > 0 ) write(20,'(1x,A)') Title
  
  do i = 1,nx
    write(20,130) x(i), z(i,:)
  end do
  
  Close(20)
    
  deallocate( x, y, z )
    
  return
  110 format(//' In the gaussian broadening, the width =',1p,e11.3,' is negative, what is not possible !'//)
  120 format(//' In the file:',/A,/'  the columns must have the same number of elements',//)
  130 format(f15.7,1p,1000e13.5) 
end

!***********************************************************************

subroutine Conv_gaussian(ncol,nx,sigma,x,y,z)

  use declarations
  implicit none

  integer:: i, j, ncol, nx
  
  real(kind=db):: alfa, f, g, g_sum, mu, sigma
  real(kind=db), dimension(nx):: dx, x
  real(kind=db), dimension(nx,ncol):: y, z

  if( abs( sigma ) < eps10 ) then
    z(:,:) = y(:,:)
    return
  endif
   
  alfa = 0.5_db / sigma**2
  f = 1 / ( sigma * sqrt( 2._db * pi ) )
  
  z(:,:) = 0._db

  dx(1) = ( x(2) - x(1) ) / 2
  do i = 2,nx-1
    dx(i) = ( x(i+1) - x(i-1) ) / 2
  end do
  dx(nx) = ( x(nx) - x(nx-1) ) / 2
    
  do i = 1,nx
  
    mu = x(i)
    
    g_sum = 0
    
    do j = 1,nx
    
      g = f * exp( - alfa * ( x(j) - mu )**2 ) * dx(j)
      g_sum = g_sum + g
      
      z(i,:) = z(i,:) + g * y(j,:)
      
    end do
    
    z(i,:) = z(i,:) / g_sum
    
  end do

  return
end

