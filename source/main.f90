! FDMNES II program, Yves Joly, Oana Bunau, Yvonne Soldo-Olivier, 19th of January 2018, 30 Nivose, An 226
!                 Institut Neel, CNRS - Universite Grenoble Alpes, Grenoble, France.
! MUMPS solver inclusion by S. Guda, A. Guda, M. Soldatov et al., University of Rostov-on-Don, Russia
! FDMX extension by J. Bourke and Ch. Chantler, University of Melbourne, Australia

! Program performing calculations of x-ray spectroscopies, XANES, RXD, dichroism, Raman X-ray scattering, SRXRD
! Work using the finite difference method or the multiple scattering theory
! Monoelectronic approach or TD-DFT (LSDA)

! Main routines of the FDMNES package
! Need also :
!   clemf0.f90, coabs.f90, convolution.f90, diffraction.f90, dirac.f90, fdm.f90, fprime.f90, frime_data.f90, general.f90,
!   lecture.f90, mat.f90, metric.f90, minim.f90, optic.f90, potential.f90, scf.f90 selec.f90,
!   spgroup.f90, sphere.f90, tab_data.f90, tensor.f90, tddft.f90

! When using the MUMPS library, one also needs:
!   mat_solve_MUMPS.f and MUMPS (and associated SCOTCH and METIS), BLAS and LAPACK libraries.
!   Inclusion of MUMPS library is due to S. Guda, A. Guda, M. Soldatov et al, University of Rostov-on-Don, Russia
!   MUMPS is from P. Amestoy et al (http://mumps-solver.org/)

! When using gaussian solver, (what is less efficient) one also needs:
!   mat_solve_gaussian.f and the BLAS and LAPACK libraries.
!   The BLAS and LAPACK library can be replaced by sub_util.f

! When working on sequential mode one also needs not_mpi.f (but with MUMPS which produces its own fake MPI files)

!***********************************************************************

! Declarations imported in most routines by instruction "use".

module declarations

  implicit none

  integer, parameter:: db = selected_real_kind(15)
  integer, parameter:: sg = selected_real_kind(6)

  integer:: MPI_COMM_MUMPS, MPI_COMM_GATHER

  integer, parameter:: Length_word = 15 ! number of character in words

  integer, parameter:: nassm = 103   ! Number of chemical elements in the Mendeleiev table
  integer, parameter:: nrepm = 12    ! Max number of representation
  integer, parameter:: nopsm = 64    ! Number of symmetry operation

  character(len=50), parameter:: Revision = 'FDMNES II program, Revision 19th of January 2018'
  character(len=16), parameter:: fdmnes_error = 'fdmnes_error.txt'

  complex(kind=db), parameter:: img = ( 0._db, 1._db )

! Physical constant
  real(kind=db), parameter:: alfa_sf = 0.0072973525664_db  ! Fine struct. cst. = 0.007 297 352 5664(17) = 1/137.035999139 CODATA 2014
  real(kind=db), parameter:: bohr = 0.52917721067_db       ! Bohr radius       = 0.529 177 210 67(12) x 10^(-10) m        CODATA 2014
  real(kind=db), parameter:: hbar = 1.054571800_db         ! Planck cst / 2 pi = 1.054 571 800(13) x 10^(-34) J s         CODATA 2014
  real(kind=db), parameter:: k_Boltzmann = 1.38064852_db   ! Boltzman constant = 1.380 648 52(79) x 10^(-23) J/K          CODATA 2014
  real(kind=db), parameter:: e_electron = 1.6021766208_db  ! Electron charge   = 1.602 176 6208(98) x  10^(-19) C         CODATA 2014
  real(kind=db), parameter:: m_electron = 9.10938356_db    ! Electron mass     = 9.109 383 56(11) x 10^(31) kg            CODATA 2014
  real(kind=db), parameter:: atom_mu = 1.660539040_db      ! Atomic mass unit  = 1.660 539 040(20) x 10^(-27) kg          CODATA 2014
  real(kind=db), parameter:: c_light = 2.99792458_db       ! Speed of light    = 2.997 924 58 x 10**8 m/s                 Exact

! Rydb = R_Rydb_lim * c * h / e = m (c.alfa)**2 / 2 e = 13.6056930129052 = convertion factor from Rydberg to eV
  real(kind=db), parameter:: Rydb = 5000 * m_electron * c_light**2 * alfa_sf**2 / e_electron

  real(kind=db), parameter:: pi = 3.1415926535897932384626433832_db
  real(kind=db), parameter:: deux_pi = 2 * pi
  real(kind=db), parameter:: quatre_pi = 4 * pi
  real(kind=db), parameter:: huit_pi = 8 * pi
  real(kind=db), parameter:: radian = pi / 180._db
  real(kind=db), parameter:: sqrt3s2 = sqrt( 3._db ) / 2
  real(kind=db), parameter:: msqrt3s2 = - sqrt( 3._db ) / 2

  real(kind=db), parameter:: eps4 = 1.e-4_db
  real(kind=db), parameter:: eps6 = 1.e-6_db
  real(kind=db), parameter:: eps10 = 1.e-10_db
  real(kind=db), parameter:: eps15 = 1.e-15_db
  real(kind=db), parameter:: epspos = 1.e-4_db     ! precision on the atom and point positions

end module declarations

!***********************************************************************

! Main of the FDMNES program

program fdmnes

  use declarations
  implicit none
  include 'mpif.h'

  integer:: i, ipr, istat, j, k, l, mpierr, mpinodes0, mpirank0, n, &
            ncalcul, nnombre

  character(len=1):: mot1
  character(len=8):: dat
  character(len=10):: tim
  character(len=11):: fdmfile
  character(len=50):: com_date, com_time
  character(len=132):: identmot, mot
  character(len=132), dimension(:), allocatable:: fdmnes_inp

! mpinodes0 is the total number of nodes for parallel computing
! mpinodes is the number of nodes for the energy loop for parallel computing
  call MPI_Init( mpierr )
  call MPI_Comm_Size(MPI_COMM_WORLD,mpinodes0,mpierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD,mpirank0,mpierr)

  fdmfile = 'fdmfile.txt'

  if( mpirank0 == 0 ) then

    call date_and_time( date = dat, time = tim )
    com_date = '   Date = ' // dat(7:8) // ' ' // dat(5:6) // ' ' // dat(1:4)
    com_time = '   Time = ' // tim(1:2)  // ' h ' // tim(3:4) // ' mn ' // tim(5:6) // ' s'

    write(6,'(1x,A/A/A)') Revision, com_date, com_time

    open(1, file = fdmfile, status='old', iostat=istat)
    if( istat /= 0 ) call write_open_error(fdmfile,istat,1)

    do i = 1,1000
      n = nnombre(1,132)
      if( n > 0 ) exit
      read(1,*)
    end do

    read(1,*) ncalcul
    allocate( fdmnes_inp(ncalcul) )

    j = 0
    do i = 1,ncalcul
      do k = 1,100
        n = nnombre(1,132)
        read(1,'(A)',iostat=istat) mot
        if( istat /= 0 ) then
          call write_error
          do ipr = 6,9,3
            write(ipr,110) i, ncalcul
          end do
          stop
        endif
        mot1 = identmot(mot,1)
        if( mot1 == ' ' ) mot = adjustl( mot )
        mot1 = identmot(mot,1)
        if( mot1 == '!' ) cycle
        exit
      end do
      Open(2,file=mot,status='old',iostat=istat)
      if( istat /= 0 ) then
        l = len_trim(mot)
        if( l > 4 ) then
          if( mot(l-3:l-3) /= '.') then
            mot(l+1:l+4) = '.txt'
            Open(2,file=mot,status='old',iostat=istat)
            if( istat /= 0 ) mot(l+1:l+4) = '    '
          endif
        endif
      endif
      if( istat /= 0 ) then
        call write_open_error(mot,istat,0)
        close(9)
      else
        j = j + 1
        fdmnes_inp(j) = mot
      endif
      Close(2)
    end do
    Close(1)
    ncalcul = j
  endif

  if( mpinodes0 > 1 ) then
! One does not need to broadcast fdmnes_inp because the reading of the indata file is done only by mpirank0 == 0
    call MPI_Bcast(ncalcul,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  endif

  if( mpirank0 /= 0 ) allocate( fdmnes_inp(ncalcul) )

  do i = 1,ncalcul
    call fit(fdmnes_inp(i),mpirank0,mpinodes0)
  end do

  deallocate( fdmnes_inp )

  call MPI_FINALIZE(mpierr)

  110 format(//' Problem when reading the file names in fdmfile.txt !'/, &
               '   The problem is for the file number',i4,',',/ &
               '   the total number of files is supposed to be',i4,'.'//)
end

!***********************************************************************

subroutine Write_error

  use declarations
  implicit none

  character(len=8):: dat
  character(len=10):: tim
  character(len=50):: com_date, com_time

  call date_and_time( date = dat, time = tim )
  com_date = '   Date = ' // dat(7:8) // ' ' // dat(5:6) // ' ' // dat(1:4)
  com_time = '   Time = ' // tim(1:2)  // ' h ' // tim(3:4) // ' mn ' // tim(5:6) // ' s'

  open(9, file = fdmnes_error)
  write(9,'(1x,A/A/A)') Revision, com_date, com_time

  return
end

!***********************************************************************

subroutine Write_open_error(File_name,istat,istop)

  implicit none

  integer:: ipr, istat, istop

  character(len=*):: File_name

  call write_error
  do ipr = 6,9,3
    write(ipr,110) File_name, istat
  end do
  if( istop == 1 ) stop

  return
  110 format(//' Error opening the file:',//3x,A,// 10x,'Status =',i5,// &
           10x,'It does not exist or ',/ &
           10x,'it is not in the good directory or',/ &
           10x,'there is some bad spelling or',/ &
           10x,'the line contains some extra hidden character !'//)
end

!***********************************************************************

! Contain the different steps of option of the calculation:
! 1) Calculation of not convoluted signal
! 2) convolution
! 3) comparison with data using metric distances
! 4) fitting
! 5) Multiplication of unit cells
! 6) Selection of scan or spectra in scan DAFS files

subroutine Fit(fdmnes_inp,mpirank0,mpinodes0)

  use declarations
  implicit none
  include 'mpif.h'

  integer, parameter:: nkw_all = 38
  integer, parameter:: nkw_fdm = 202
  integer, parameter:: nkw_conv = 36
  integer, parameter:: nkw_fit = 1
  integer, parameter:: nkw_gaus = 1
  integer, parameter:: nkw_metric = 12
  integer, parameter:: nkw_mult = 4
  integer, parameter:: nkw_selec = 5
  integer, parameter:: nmetricm = 4
  integer, parameter:: nparam_conv = 11
  integer, parameter:: nparam_fdm = 17

  integer, parameter:: nparam_tot = nparam_conv + nparam_fdm

  integer:: eof, i, i1, ibl, ical, ifdm, igr, index_Met_Fit, indpar, &
    inotskip, ip, ipar, ipbl, ipr, istat, istop, itape, itape_minim, itape1, itape2, itape3, itape4, itape5, itape6, itape7, &
    itph, itpj, itpm, itps, iscratch, iscratchconv, j, jgr, jpar, k, l, Length_line, &
    ligne, ligne2, m, mpierr, mpirank0, mpinodes0, n, &
    n_atom_proto_p, n_shift, n1, n2, nb_datafile, nblock, ncal, ncal_nonfdm, ndem, ndm, ng, ngamh, ngroup_par, ngroup_par_conv, &
    nmetric, nn, nnombre, nnotskip, nnotskipm, nparm, npm, mermrank

  character(len=9):: keyword, mot9, Traduction
  character(len=13):: mot13
  character(len=132):: comt, convolution_out, elf_infile, File_in, fdmfit_out, fdmnes_inp, Fichier, identmot, imfp_infile, &
    mot, nomfich, nomfichbav
  character(len=1320):: mot1320
  character(len=2), dimension(nmetricm):: Nom_Met
  character(len=9), dimension(nkw_all):: kw_all
  character(len=9), dimension(nkw_fdm):: kw_fdm
  character(len=9), dimension(nkw_conv):: kw_conv
  character(len=9), dimension(nkw_fit):: kw_fit
  character(len=9), dimension(nkw_gaus):: kw_gaus
  character(len=9), dimension(nkw_metric):: kw_metric
  character(len=9), dimension(nkw_selec):: kw_selec
  character(len=9), dimension(nkw_mult):: kw_mult
  character(len=9), dimension(nparam_tot):: param_conv
  character(len=9), dimension(:,:), allocatable:: typepar, typeparc
  character(len=9), dimension(:), allocatable:: typeparg

  integer, dimension(3):: hkl_borm
  integer, dimension(30):: icheck
  integer, dimension(nmetricm):: ical_Met_min
  integer, dimension(:), allocatable:: ifile_notskip, indparp, Length_block, npar, npbl, nparam
  integer, dimension(:,:), allocatable:: indice_par

  logical:: bav_open, Bormann, Case_fdm, Check_file, Conv_done, &
    Convolution_cal, Dafs_bio, E_cut_man, Fdmnes_cal, Fit_cal, Gamma_hole_imp, Gamma_tddft, Gaus_cal, Metric_cal, &
    Minim_fdm_ok, minimok, Mult_cal, Scan_a, Selec_cal, &
    Use_FDMX, FDMX_only, cm2g, nobg, nohole, nodw, noimfp, imfp_inp, elf_inp, dwfactor_inp, tdebye_inp, tmeas_inp, &
    expntl, victoreen

  logical, dimension(:), allocatable:: block_sum

  real(kind=db):: Ang_borm, Delta_edge, E_cut_imp, e1, e2, Ecent, &
    Elarg, Estart, Gamma_max, fac, Gmin, param_dep, prop, rtph, tp_deb, tp_fin, tpt, x, &
    dwfactor, tdebye, tmeas, expntlA, expntlB, victA, victB
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(nmetricm):: Dist_Min, Gen_Shift_min
  real(kind=db), dimension(:), allocatable:: par_op, parsum, RapIntegrT_min_g
  real(kind=db), dimension(:,:), allocatable:: Dist_min_g, par, param, parmax, parmin

  real(kind=sg) time

! Keywords of the indata file
  data kw_all /  'bormann  ','check    ','check_all','check_coa','check_con', &
     'check_pot','check_mat','check_sph','check_tdd','check_ten','comment  ','delta_edg','ecent    ','e_cut    ','elarg    ', &
     'estart   ','file_in  ','filout   ','fprime_at','gamma_hol','gamma_max','length_li','no_check ','imfpin   ','elfin    ', &
     'dwfactor ','tdebye   ','tmeas    ','expntl   ','victoreen','mermin   ', &
     'fdmx     ','fdmx_proc','cm2g     ','nobg     ','nohole   ','nodw     ','noimfp   '/

  data kw_conv / 'abs_befor','all_conv ','no_analyz','cal_tddft','calculati','circular ','conv_out ','convoluti','dafs_exp_', &
     'dead_laye','dec      ','directory','double_co','eintmax  ','epsii    ','forbidden','fprime   ', &
     'gamma_fix','gamma_var','gaussian ','no_extrap','nxan_lib ','photo_emi','s0_2     ','selec_cor','sample_th','scan     ', &
     'scan_conv','scan_file','seah     ','stokes   ','stokes_na','surface_p','table    ','thomson  ','transpose'/

  data kw_fdm/  &
     'absorbeur','adimp    ','all_nrixs','allsite  ','ata      ','atom     ','atom_conf','atom_nsph','ang_spin ','atomic_sc', &
     'axe_spin ','base_comp','base_reel','bond     ','bulk     ','bulk_roug','cap_disor','cap_rough','cap_layer','cap_shift', &
     'cap_thick','cartesian','center   ','center_ab','chlib    ','cif_file ','classic_i','clementi ','core_reso','crystal  ', &
     'crystal_c','crystal_t','d_max_pot','dafs     ','dafs_2d  ','dafs_exp ','debye    ','delta_en_','dip_rel  ','e1e1     ', &
     'delta_eps','density  ','density_a','density_c','dilatorb ','dipmag   ','doping   ','dpos     ','dyn_g    ','dyn_eg   ', &
     'edge     ','e1e2     ','e1e3     ','e1m1     ','e1m2     ','e2e2     ','e3e3     ','eimag    ','eneg     ','energphot', &
     'ephot_min','etatlie  ','excited  ','extract  ','extract_t','extractpo','extractsy','fdm_comp ','film     ','film_cif_', &
     'film_pdb_','film_t   ','film_roug','film_shif','film_zero','flapw    ','flapw_n  ','flapw_n_p','flapw_psi','flapw_r  ', &
     'flapw_s  ','flapw_s_p','full_atom','full_pote','full_self','gamma_tdd','green    ','green_int','hedin    ','helmholtz', &
     'hkl_film ','hubbard  ','iord     ','kern_fac ','kern_fast', &
     'lmax     ','lmax_nrix','lmax_tddf','lmaxfree ','lmaxso   ','lmaxstden','ldipimp  ','lmoins1  ','lplus1   ','mat_ub   ', &
     'memory_sa','lquaimp  ','m1m1     ','m1m2     ','m2m2     ','magnetism','mat_polar', &
     'molecule ','molecule_','muffintin','multrmax ','n_self   ','nchemin  ','new_zero ','no_core_r','no_dft   ','no_e1e1  ', &
     'no_e1e2  ','no_e1e3  ','no_e2e2  ','no_e3e3  ','no_fermi ','no_renorm','no_res_ma','no_res_mo','no_solsin','normaltau', &
     'norman   ','noncentre','non_relat','nonexc   ','not_eneg ','nrato    ','nrixs    ','occupancy','octupole ','old_zero ', &
     'one_run  ','optic    ','over_rad ','overlap  ','p_self   ','p_self_ma','pdb_file ','perdew   ','pointgrou','polarized', &
     'quadmag  ','quadrupol','radius   ','range    ','rangel   ','raydem   ','rchimp   ','readfast ','relativis', &
     'rmt      ','rmtg     ','rmtv0    ','rot_sup  ','rpalf    ','rpotmax  ','r_self   ','rydberg  ', &
     'self_abs ','scf      ','scf_abs  ','scf_exc  ','scf_mag_f','scf_non_e','scf_step ', &
     'screening','setaz    ','solsing  ','spgroup  ','sphere_al','spherical','spinorbit','step_azim','supermuf ','surface  ', &
     'surface_s','surface_t','symmol   ','symsite  ','tddft    ','temperatu','test_dist','trace    ','vmax     ','v0imp    ', &
     'xalpha   ','xan_atom ','ylm_comp ','z_absorbe','z_nospino','zero_azim'/

  data kw_fit / 'parameter'/

  data kw_gaus / 'conv_gaus'/

  data kw_metric/ 'd2       ','detail   ','emin     ','emax     ','experimen','file_met ','fit_out  ','gen_shift','kev      ', &
                  'metric_ou','rx       ','rxg      '/

  data kw_selec/ 'selec_inp','selec_out','energy   ','azimuth  ', 'reflectio'/

  data kw_mult / 'mult_cell','unit_cell','atomic_nu','surf_cell'/

! First the convolution parameters, then the others
  data param_conv / &
     'ecent    ','e_cut    ','elarg    ','gamma_hol','gamma_max','gaussian ','shift    ','aseah    ','bseah    ','vibr     ', &
     'weight   ','a        ','abc      ','anga     ','angb     ','angc     ','b        ','c        ','dposx    ','dposy    ', &
     'dposz    ','occup    ','phi      ','poporb   ','posx     ','posy     ','posz     ','theta    '/

  call CPU_TIME(time)
  tp_deb = real(time,db)

  Ang_borm = 0._db
  Bormann = .false.
  Check_file = .false.
  hkl_borm(:) = 0
  icheck(:) = 1
! icheck : 1: Lecture        2: Atom et Dirac    3: Symsite
!         10: Pot0          13: Potentiel
!         18: Sphere        19: Mat et MSM      20: Tensor
!         21: Coabs         22: tddft - Sphere  23: tddft - Chi_0   24: tddft - Kernel
!         25: tddft - Chi   26: Hubbard         27: SCF             28: State           29: Optic           30: convolution
  Ecent = 30._db / Rydb
  Elarg = 30._db / Rydb
  Estart = 100000._db / Rydb
  Gamma_hole(:) = 0._db
  Gamma_hole_imp = .false.
  Gamma_max = 15._db / Rydb
  Gamma_tddft = .false.
  Gaus_cal = .false.
  Length_line = 10 + 10001 * Length_word
  Minim_fdm_ok = .false.
  Minimok = .false.
  ngamh = 1
  nomfich = 'fdmnes_out'
  Scan_a = .false.
!*** JDB
  imfp_inp = .false.
  elf_inp = .false.
  dwfactor_inp = .false.
  tdebye_inp = .false.
  tmeas_inp = .false.
  expntl = .false.
  victoreen = .false.
  mermrank = 0
  Use_FDMX = .false.
  FDMX_only = .false.
  cm2g = .false.
  nobg = .false.
  nohole = .false.
  nodw = .false.
  noimfp = .false.
!*** JDB

  itape = 6
  itape1 = 11  ! Convolution
  itape2 = 12  ! Metric calculation
  itape3 = 13  ! Fitting
  itape4 = 14  ! FDM calculation
  itape5 = 17  ! Selection in DAFS scan
  itape6 = 18  ! Multiplication of unit cell
  itape7 = 19  ! Convolation by a gaussian
  itape_minim = 10
  iscratch = 15
  iscratchconv = 35

  comt = ' '
  Convolution_cal = .false.
  Delta_edge = 0._db
  Fdmnes_cal = .false.
  E_cut_imp = -5._db / rydb
  E_cut_man = .false.
  e1 = 0._db; e2 = 0._db
  Fit_cal = .false.
  Metric_cal = .false.
  Mult_cal = .false.
  Selec_cal = .false.
  ngroup_par = 0
  nparm = 1
  n_shift = 1
  nb_datafile = 0
  ng = 0

  if( mpirank0 == 0 ) then

  open(1, file=fdmnes_inp, status = 'old')

  boucle_ligne: do ligne = 1,100000

    read(1,'(A)',iostat=eof) mot

    if( eof /= 0 ) exit boucle_ligne

    keyword = identmot(mot,9)
    keyword = traduction(keyword)

    if( keyword(1:1) == '!' .or. keyword(1:1) == ' ' .or. keyword == 'endjump' ) cycle

    if( keyword == 'end' ) exit

    if( keyword == 'jump' ) then
      do i = 1,100000
        read(1,'(A)',iostat=eof) mot
        if( eof /= 0 ) exit boucle_ligne
        keyword = identmot(mot,9)
        keyword = traduction(keyword)
        if( keyword == 'endjump' ) exit
      end do
      read(1,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit boucle_ligne
      keyword = identmot(mot,9)
      keyword = traduction(keyword)
    endif

    boucle_k: do k = 1,8
      select case(k)

        case(1)
          do i = 1,nkw_fdm
            if( keyword /= kw_fdm(i) ) cycle
            itape = itape4
            if( .not. Fdmnes_cal ) then
              if( Check_file ) then
                Open( itape )
              else
                Open( itape, status='SCRATCH' )
              endif
              Fdmnes_cal = .true.
            endif
            exit boucle_k
          end do

        case(2)
          do i = 1,nkw_conv
            if( keyword /= kw_conv(i) ) cycle
            itape = itape1
            if( .not. Convolution_cal ) then
              Convolution_cal = .true.
              if( Check_file ) then
                Open( itape )
              else
                Open( itape, status='SCRATCH' )
              endif
            endif
            exit boucle_k
          end do

        case(3)
          do i = 1,nkw_metric
            if( keyword /= kw_metric(i) ) cycle
            itape = itape2
            if( .not. metric_cal ) then
              metric_cal = .true.
              if( Check_file ) then
                Open( itape )
              else
                Open( itape, status='SCRATCH' )
              endif
            endif
            exit boucle_k
          end do

        case(4)
          do i = 1,nkw_fit
            if( keyword /= kw_fit(i) ) cycle
            itape = itape3
            if( .not. Fit_cal ) then
              Fit_cal = .true.
              if( Check_file ) then
                Open( itape )
              else
                Open( itape, status='SCRATCH' )
              endif
            endif
            exit boucle_k
          end do

        case(5)
          do i = 1,nkw_selec
            if( keyword /= kw_selec(i) ) cycle
            itape = itape5
            if( .not. selec_cal ) then
              selec_cal = .true.
              if( Check_file ) then
                Open( itape )
              else
                Open( itape, status='SCRATCH' )
              endif
            endif
            exit boucle_k
          end do

        case(6)
          do i = 1,nkw_mult
            if( keyword /= kw_mult(i) ) cycle
            itape = itape6
            if( .not. Mult_cal ) then
              Mult_cal = .true.
              if( Check_file ) then
                Open( itape )
              else
                Open( itape, status='SCRATCH' )
              endif
            endif
            exit boucle_k
          end do

        case(7)
          do i = 1,nkw_gaus
            if( keyword /= kw_gaus(i) ) cycle
            itape = itape7
            if( .not. Gaus_cal ) then
              Gaus_cal = .true.
              if( Check_file ) then
                Open( itape )
              else
                Open( itape, status='SCRATCH' )
              endif
            endif
            exit boucle_k
          end do

        case(8)
          do i = 1,nkw_all
            if( keyword /= kw_all(i) ) cycle
            exit boucle_k
          end do

      end select
    end do boucle_k

    if( k /= 8 ) then
      write(itape,'(A)') keyword
    else

      select case(keyword)

        case('check')
          n1 = 1
          do i = 1,30
            n = nnombre(1,132)
            n2 = min( n1+n-1, 30)
            if( n2 >= n1 ) then
              read(1,*) icheck(n1:n2)
              n1 = n2 + 1
            else
              exit
            endif
          end do

        case('check_pot')
          n = nnombre(1,132)
          if( n == 0 ) then
            n = 3
          else
           read(1,*) n
         endif
         icheck(10:13) = n; icheck(16) = n

        case('check_sph')
          n = nnombre(1,132)
          if( n == 0 ) then
            n = 3
          else
           read(1,*) n
          endif
          icheck(18) = n

        case('check_mat')
          n = nnombre(1,132)
          if( n == 0 ) then
            n = 3
          else
           read(1,*) n
          endif
          icheck(19) = n

        case('check_tdd')
          n = nnombre(1,132)
          if( n == 0 ) then
            n = 3
          else
           read(1,*) n
          endif
          icheck(22:25) = n

        case('check_ten')
          n = nnombre(1,132)
          if( n == 0 ) then
            n = 3
          else
           read(1,*) n
          endif
          icheck(20) = n

        case('check_coa')
          n = nnombre(1,132)
          if( n == 0 ) then
            n = 3
          else
           read(1,*) n
          endif
          icheck(21) = n

        case('check_con')
          n = nnombre(1,132)
          if( n == 0 ) then
            n = 3
          else
           read(1,*) n
          endif
          icheck(30) = n

        case('fprime_at')
          n = nnombre(1,132)
          if( n == 0 ) then
            n = 3
          else
           read(1,*) n
          endif
          icheck(30) = n

        case('check_all')
          n = nnombre(1,132)
          if( n == 0 ) then
            n = 3
          else
           read(1,*) n
          endif
          icheck(:) = n

         case('no_check')
          icheck(:) = 0

       case('e_cut','efermi')
          read(1,*) E_cut_imp
          E_cut_imp = E_cut_imp /rydb
          E_cut_man = .true.

        case('length_li')
          n = nnombre(1,132)
          read(1,*) Length_line
          n = 10 + 201 * Length_word
          Length_line = max( Length_line, n )

        case('comment')
          n = nnombre(1,132)
          read(1,'(A)') comt

        case('file_in')
          n = nnombre(1,132)
          read(1,'(A)') mot
          File_in = adjustl( mot )
          mot = File_in
          open(20, file = mot, status='old', iostat=istat)
          if( istat /= 0 ) then
            l = len_trim(mot)
            if( mot(l-3:l) /= '.txt' ) then
              mot(l+1:l+4) = '.txt'
              open(20, file = mot, status='old', iostat=istat)
              if( istat /= 0 ) then
                mot(l+2:l+4) = '.dat'
                open(20, file = mot, status='old', iostat=istat)
              endif
            endif
            if( istat /= 0 ) call write_open_error(File_in,istat,1)
          endif
          Close(20)
          File_in = mot

        case('filout')
          n = nnombre(1,132)
          read(1,'(A)') mot
          nomfich = adjustl( mot )
          l = len_trim(nomfich)
          if( l > 3 ) then
            if( nomfich(l-3:l) == '.txt') nomfich(l-3:l) = '    '
          endif

        case('ecent')
          n = nnombre(1,132)
          read(1,*) Ecent
          Ecent = Ecent / rydb

        case('elarg')
          n = nnombre(1,132)
          read(1,*) Elarg
          Elarg = Elarg / rydb

        case('estart')
          n = nnombre(1,132)
          read(1,*) Estart
          Estart = Estart / rydb

        case('gamma_hol')
          n = nnombre(1,132)
          gamma_hole_imp = .true.
          ngamh = min(n,10)
          read(1,*) Gamma_hole(1:ngamh)
          Gamma_hole(1:ngamh) = Gamma_hole(1:ngamh) / rydb
          Gmin = Gamma_hole(1)
          do i = 2,ngamh
            Gmin = min( Gmin, Gamma_hole(i) )
          end do
          do i = 1,ngamh
            if( abs(Gmin - Gamma_hole(i)) < eps10 ) cycle
            Gamma_tddft = .true.
            exit
          end do

        case('gamma_max')
          n = nnombre(1,132)
          read(1,*) Gamma_max
          Gamma_max = Gamma_max / rydb

        case('delta_edg')
          n = nnombre(1,132)
          read(1,*) Delta_edge
          Delta_edge = Delta_edge / rydb

        case('bormann')
          Bormann = .true.
          n = nnombre(1,132)
          read(1,*) hkl_borm(:), Ang_borm

!*** JDB
        case('imfpin')
          imfp_inp = .true.
          read(1,'(A)') imfp_infile

        case('elfin')
          elf_inp = .true.
          read(1,'(A)') elf_infile

        case('dwfactor')
          dwfactor_inp = .true.
          read(1,*) dwfactor

        case('tdebye')
          tdebye_inp = .true.
          read(1,*) tdebye

        case('tmeas')
          tmeas_inp = .true.
          read(1,*) tmeas

        case('expntl')
          expntl = .true.
          read(1,*) expntlA, expntlB

        case('victoreen')
          victoreen = .true.
          read(1,*) victA, victB

        case('mermin')
          read(1,*) mermrank

        case('fdmx')
          Use_FDMX = .true.
        case('fdmx_proc')
          FDMX_only = .true.
        case('cm2g')
          cm2g = .true.
        case('nobg')
          nobg = .true.
        case('nohole')
          nohole = .true.
        case('nodw')
          nodw = .true.
        case('noimfp')
          noimfp = .true.

      end select
      cycle
    endif

    boucle_j: do j = 1,1000000

      read(1,'(A)',iostat=eof) mot1320
      if( eof /= 0 ) exit boucle_ligne
      backspace(1)
      read(1,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit boucle_ligne
      mot9 = identmot(mot,9)

      if( mot9(1:1) == '!' .or. mot9(1:1) == ' ') cycle
      mot9 = traduction(mot9)
      if( mot9 == 'end' ) exit boucle_ligne
      if( mot9 == 'jump' ) then
        do i = 1,100000
          read(1,'(A)',iostat=eof) mot
          if( eof /= 0 ) exit boucle_ligne
          mot9 = identmot(mot,9)
          mot9 = traduction(mot9)
          if( mot9 == 'endjump' ) exit
        end do
        read(1,'(A)',iostat=eof) mot1320
        if( eof /= 0 ) exit boucle_ligne
        backspace(1)
        read(1,'(A)',iostat=eof) mot
        if( eof /= 0 ) exit boucle_ligne
        mot9 = identmot(mot,9)
        mot9 = traduction(mot9)
        if( mot9(1:1) == '!' .or. mot9(1:1) == ' ' ) exit
        if( mot9 == 'end' ) exit
      endif

      do k = 1,nkw_fdm
        if( mot9 == kw_fdm(k) ) exit boucle_j
      end do
      do k = 1,nkw_conv
        if( mot9 == kw_conv(k) ) exit boucle_j
      end do
      do k = 1,nkw_metric
        if( mot9 == kw_metric(k) ) exit boucle_j
      end do
      do k = 1,nkw_fit
        if( mot9 /= kw_fit(k) ) cycle
        if( mot9 == 'experimen' ) then
          mot13 = identmot(mot,13)
          if(mot13 /= 'experiment' .and. mot13 /= 'experimen') cycle
        endif
        exit boucle_j
      end do
      do k = 1,nkw_selec
        if( mot9 == kw_selec(k) ) exit boucle_j
      end do
      do k = 1,nkw_mult
        if( mot9 == kw_mult(k) ) exit boucle_j
      end do
      do k = 1,nkw_gaus
        if( mot9 == kw_gaus(k) ) exit boucle_j
      end do
      do k = 1,nkw_all
        if( mot9 == kw_all(k) ) exit boucle_j
      end do

      if( mot9(1:4) == 'par_' ) then
        mot13 = identmot(mot,13)
        mot9 = mot13(5:13)
        mot9 = traduction(mot9)
        mot1320 = ' '
        mot1320(2:10) = mot9
      endif

      write(itape,'(A)') mot1320

    end do boucle_j

    backspace(1)

  end do boucle_ligne

  Close(1)

  Dafs_bio = .false.
  if( Metric_cal ) then
    Rewind(itape2)
    ngroup_par = 1
    boucle_m: do ligne = 1,100000
      read(itape2,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit boucle_m
      mot9 = identmot(mot,9)
      if( mot9 == 'gen_shift' ) then
        n = nnombre(itape2,132)
        read(itape2,*) e1, e2, n_shift
      elseif( mot9 == 'experimen' ) then
        boucle_i: do i = 1,100000
          n = nnombre(itape2,132)
          read(itape2,'(A)',iostat=eof) mot
          if( eof /= 0 ) exit boucle_m
          mot9 = identmot(mot,9)
          do j = 1,nkw_fdm
            if( mot9 == kw_fdm(j) ) exit boucle_i
          end do
          do j = 1,nkw_conv
            if( mot9 == kw_conv(j) ) exit boucle_i
          end do
          do j = 1,nkw_metric
            if( mot9 == kw_metric(j) ) exit boucle_i
          end do
          do j = 1,nkw_fit
            if( mot9 /= kw_fit(j) ) cycle
            if( mot9 == 'experimen' ) then
              mot13 = identmot(mot,13)
              if( mot13 /= 'experiment' .and. mot13 /= 'experimen' ) cycle
            endif
            exit boucle_i
          end do
          if( n == 0 ) then
            nb_datafile = nb_datafile + 1  ! number of files
            Fichier = Adjustl( mot )
            l = len_trim(Fichier)
            if( l > 4 ) then
             if( Fichier(l-3:l-3) /= '.' ) Fichier(l+1:l+4) = '.txt'
            endif
            open(99, file=Fichier, status='old',iostat=istat)
            if( istat /= 0 ) call write_open_error(Fichier,istat,1)
            n = nnombre(99,100000)
            if( n == 0 ) then
              ng = ng + 1
            else
              ng = ng + n / 3
              Dafs_bio = .true.
            endif
          endif
          Close(99)
        end do boucle_i
        backspace(itape2)
      elseif( mot9 == 'file_met' ) then
        nb_datafile = 1
        read(itape2,'(A)') mot
        Fichier = Adjustl( mot )
        l = len_trim(Fichier)
        if( l > 4 ) then
         if( Fichier(l-3:l-3) /= '.' ) Fichier(l+1:l+4) = '.txt'
        endif
        open(99, file=Fichier, status='old',iostat=istat)
        if( istat /= 0 ) call write_open_error(Fichier,istat,1)
        read(99,*)
        ng = nnombre(99,132000) - 1
      endif
    end do boucle_m
  endif

  nblock = 0
  if( Fit_cal ) then
    Rewind(itape3)
    boucle_f: do ligne = 1,100000
      read(itape3,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit boucle_f
      mot9 = identmot(mot,9)
      if( mot9 == 'parameter' ) then
        ngroup_par = ngroup_par + 1
        nblock = nblock + 1
        nn = 0

        do ligne2 = 1,100000
          n = nnombre(itape3,132)
          read(itape3,'(A)',iostat=eof) mot
          if( eof /= 0 ) exit boucle_f
          if( mot == ' ' ) cycle
          mot9 = identmot(mot,9)
          if( mot9 == 'parameter' ) then
            backspace(itape3)
            exit
          else
            if( n == 0 ) then
              cycle
            elseif( n > 2 ) then
              nn = nn + 1
            else
              ngroup_par = ngroup_par + nn - 1
              exit
            endif
          endif
        end do
      endif

    end do boucle_f
    Rewind(itape3)
  endif

  allocate( npar(ngroup_par) )
  allocate( npbl(nblock) )
  allocate( block_sum(nblock) )
  if( ngroup_par > 0 ) then
    block_sum(:) = .false.
    npar(:) = 0
  endif

  if( Fit_cal ) then

    Rewind(itape3)
    boucle_g: do ibl = 1,nblock
      n = nnombre(itape3,132)
      read(itape3,'(A)') mot
      mot9 = identmot(mot,9)
      if( mot9 == 'parameter' ) then
        do ligne = 1,100000
          n = nnombre(itape3,132)
          read(itape3,'(A)',iostat=eof) mot
          if( eof /= 0 ) exit boucle_g
          if( mot == ' ' ) cycle
          mot9 = identmot(mot,9)
          if( mot9 == 'parameter' ) then
            backspace(itape3)
            exit
          else
            if( n == 1 .or. n == 2 ) then
              block_sum(ibl) = .true.
              exit
            else
              cycle
            endif
          endif
        end do
      endif

    end do boucle_g
  endif

  if( Metric_cal ) then
    nparm = n_shift
    npar(1) = 1
  endif

  if( Fit_cal ) then

    Rewind(itape3)
    igr = 1
    npbl(:) = 0
    do ibl = 1,nblock
      i1 = 1
      n = nnombre(itape3,132)
      read(itape3,'(A)')
      if( .not. Block_sum(ibl) ) igr = igr + 1
      boucle_h: do ligne = 1,100000
        n = nnombre(itape3,132)
        read(itape3,'(A)',iostat=eof) mot
        if( eof /= 0 ) exit boucle_h
        if( n > 0 ) cycle
        mot9 = identmot(mot,9)
        if( mot9 == 'parameter' ) then
          backspace(itape3)
          exit
        else
          istop = 0
          do i = 1,nparam_tot
            if( mot9 /= param_conv(i) ) cycle
            if( .not. Fdmnes_cal .and. i > nparam_conv ) then
              if( istop == 0 ) call write_error
              do ipr = 6,9,3
                write(ipr,110) mot9
              end do
              istop = 1
            endif
            if( block_sum(ibl) ) then
              n = nnombre(itape3,132)
              if( n > 2 ) then
                igr = igr + 1
                if( i1 == 1 ) then
                  npar(igr) = 2
                  i1 = 2     ! to avoid a double addition in subroutine "Lecture"
                else
                  npar(igr) = 1
                endif
              endif
            else
              npar(igr) = npar(igr) + 1
            endif
            npbl(ibl) = npbl(ibl) + 1
            exit
          end do
          if( istop == 1 ) then
            close(9)
            stop
          endif
          if( i > nparam_tot .and. mpirank0 == 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,120) mot9
            end do
            close(9)
            stop
          endif
        endif
      end do boucle_h
      nparm = max( nparm, npar(igr) )
    end do

  endif

  endif   ! Arriving for parallel computing

  if( mpinodes0 > 1 ) then
    call MPI_Bcast(Mult_cal,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)

  endif

  if( Mult_cal ) then
! Multiplication of unit cell
    if( mpirank0 == 0 ) call mult_cell(itape6,nomfich)
    return
  elseif( Gaus_cal ) then
! Concolution by a gaussian
    if( mpirank0 == 0 ) call Main_gaussian(File_in,itape7,nomfich)
    return
  endif

  if( mpinodes0 > 1 ) then

    call MPI_Bcast(Ang_borm,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Bormann,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E_cut_imp,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E_cut_man,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(nblock,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ngroup_par,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Length_line,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(hkl_borm,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(e1,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(e2,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ecent,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Elarg,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Estart,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_shift,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nparm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ng,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nb_datafile,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Fit_cal,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Gamma_hole,10,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Gamma_hole_imp,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Gamma_max,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Gamma_tddft,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Fdmnes_cal,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ngamh,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  endif

  if( mpirank0 > 0 ) then
    allocate( npar(ngroup_par) )
    allocate( npbl(nblock) )
    allocate( block_sum(nblock) )
  endif
  allocate( Length_block(ngroup_par) )
  allocate( param(ngroup_par,nparm) )
  allocate( parmax(ngroup_par,nparm) )
  allocate( parmin(ngroup_par,nparm) )
  allocate( parsum(ngroup_par) )
  allocate( nparam(ngroup_par) )
  allocate( indice_par(ngroup_par,nparm) )
  allocate( typepar(ngroup_par,nparm) )
  allocate( typeparc(ngroup_par,nparm) )
  allocate( typeparg(ngroup_par) )
  if( ngroup_par > 0 ) then
    indice_par(:,:) = 0
    Length_block(:) = 0
    parmin(1,1) = e1
    parmax(1,1) = e2
    nparam(1) = n_shift
    typepar(1,1) = 'Gen_Shift'
    typeparg(1) = 'Gen_Shift'
  endif

  if( mpirank0 == 0 ) then
    allocate( Dist_min_g(ng,nmetricm) )
    allocate( RapIntegrT_min_g(ng) )
  endif

  if( Fit_cal .and. mpirank0 == 0 ) then

    rewind(itape3)
    igr = 1
    do ibl = 1,nblock
      n = nnombre(itape3,132)
      read(itape3,'(A)') mot
      mot9 = identmot(mot,9)
      if( .not. Block_sum(ibl) ) igr = igr + 1
      do ipbl = 1,npbl(ibl)
        if( .not. Block_sum(ibl) ) then
          ipar = ipbl
        else
          if( ipbl == npbl(ibl) ) then
            ipar = 2
          else
            ipar = 1
            igr = igr + 1
          endif
        endif
        if( Block_sum(ibl) ) then
          if( ipbl < npbl(ibl) - 1 ) then
            Length_block(igr) = - 1
          elseif( ipbl == npbl(ibl) - 1 ) then
            Length_block(igr) = npbl(ibl) - 1
          endif
        endif
        n = nnombre(itape3,132)
        read(itape3,'(A)') mot
        mot9 = identmot(mot,9)
        if( Block_sum(ibl) .and. ipar == 2 ) then
          typepar(igr-npbl(ibl)+2:igr,ipar) = mot9
        else
          typepar(igr,ipar) = mot9
        endif
        n = nnombre(itape3,132)
        select case( n )
          case(1)
            read(itape3,*) x
            parsum(igr-npbl(ibl)+2:igr) = x
          case(2)
            read(itape3,*) x, i
            parsum(igr-npbl(ibl)+2:igr) = x
            indice_par(igr-npbl(ibl)+2:igr,ipar) = i
          case(3)
            read(itape3,*) parmin(igr,ipar), parmax(igr,ipar), nparam(igr)
          case(4)
            read(itape3,*) parmin(igr,ipar), parmax(igr,ipar), nparam(igr), indice_par(igr,ipar)
        end select
      end do
    end do
    Close(itape3)

  endif

  if( mpinodes0 > 1 ) then
    if( ngroup_par > 0 ) then
      call MPI_Bcast(Length_block,ngroup_par,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(npar,ngroup_par,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(npbl,nblock,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    endif

    if( Fit_cal ) then
      call MPI_Bcast(nparam,ngroup_par,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(indice_par,nparm*ngroup_par,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(parmin,nparm*ngroup_par,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(parmax,nparm*ngroup_par,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(parsum,ngroup_par,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      do igr = 1,ngroup_par
        do ipar = 1,npar(igr)
          if( mpirank0 /= 0 ) mot = ' '
          if( mpirank0 == 0 ) mot = typepar(igr,ipar)
          do i = 1,9
            if( mpirank0 == 0 ) j = iachar( mot(i:i) )
            call MPI_Bcast(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

            if( mpirank0 /= 0 ) mot(i:i) = achar( j )
          end do
          if( mpirank0 /= 0 ) typepar(igr,ipar) = mot
        end do
      end do
    endif
  endif

  ncal = 1
  if( Fit_cal ) then

    do igr = 2,ngroup_par
      if( abs(parmin(igr,1) - parmax(igr,1)) < 0.00000001_db ) nparam(igr) = 1
    end do

    do igr = 2,ngroup_par
      nparam(igr) = max( nparam(igr), 1 )
      ncal = ncal * nparam(igr)
    end do

  endif

! Parameter ordering : first shift, then conv, then fdm.
  boucle_igr: do igr = 2,ngroup_par
    do ipar = 1,npar(igr)
      do i = 1,nparam_conv
        if( typepar(igr,ipar) == param_conv(i) ) exit
      end do
      if( i > nparam_conv ) exit
    end do
    Case_fdm = .true.
    if( i <= nparam_conv ) then
! Test if only shift or gaus
      do ipar = 1,npar(igr)
        if( typepar(igr,ipar) /= 'shift' .and. typepar(igr,ipar) /= 'weight' .and. typepar(igr,ipar) /= 'gaussian' ) exit
      end do
      if( ipar > npar(igr) ) cycle
      Case_fdm = .false.
    endif

    do jgr = igr+1,ngroup_par
      if( Case_fdm ) then
        do jpar = 1,npar(jgr)
          do j = 1,nparam_conv
            if( typepar(jgr,jpar) == param_conv(j) ) exit
          end do
          if( j > nparam_conv ) exit
        end do
        if( j > nparam_conv ) cycle
      else
        do jpar = 1,npar(jgr)
          if( typepar(jgr,jpar) /= 'shift' .and. typepar(jgr,jpar) /= 'weight' .and. typepar(jgr,jpar) /= 'gaussian' ) exit
        end do
        if( jpar <= npar(jgr) ) cycle
      endif

      n = Length_block(igr)
      Length_block(igr) = Length_block(jgr)
      Length_block(jgr) = n
      n = nparam(igr)
      nparam(igr) = nparam(jgr)
      nparam(jgr) = n
      n = npar(igr)
      npar(igr) = npar(jgr)
      npar(jgr) = n
      n = min(npar(igr),npar(jgr))
      do ipar = 1,n
        x = parmin(igr,ipar)
        parmin(igr,ipar) = parmin(jgr,ipar)
        parmin(jgr,ipar) = x
        x = parmax(igr,ipar)
        parmax(igr,ipar) = parmax(jgr,ipar)
        parmax(jgr,ipar) = x
        mot9 = typepar(igr,ipar)
        typepar(igr,ipar) = typepar(jgr,ipar)
        typepar(jgr,ipar) = mot9
        m = indice_par(igr,ipar)
        indice_par(igr,ipar) = indice_par(jgr,ipar)
        indice_par(jgr,ipar) = m
      end do
      do ipar = n+1,npar(igr)
        parmin(igr,ipar) = parmin(jgr,ipar)
        parmin(jgr,ipar) = 0._db
        parmax(igr,ipar) = parmax(jgr,ipar)
        parmax(jgr,ipar) = 0._db
        typepar(igr,ipar) = typepar(jgr,ipar)
        typepar(jgr,ipar) = '         '
        indice_par(igr,ipar) = indice_par(jgr,ipar)
        indice_par(jgr,ipar) = 0
      end do
      do ipar = n+1,npar(jgr)
        parmin(jgr,ipar) = parmin(igr,ipar)
        parmin(igr,ipar) = 0._db
        parmax(jgr,ipar) = parmax(igr,ipar)
        parmax(igr,ipar) = 0._db
        typepar(jgr,ipar) = typepar(igr,ipar)
        typepar(igr,ipar) = '         '
        indice_par(jgr,ipar) = indice_par(igr,ipar)
        indice_par(igr,ipar) = 0
      end do
      x = parsum(igr)
      parsum(igr) = parsum(jgr)
      parsum(jgr) = x

    end do

  end do boucle_igr

  ncal_nonfdm = 1
  boucle_ext: do igr = 2,ngroup_par
    do ipar = 1,npar(igr)
      do i = 1,nparam_conv
        if( typepar(igr,ipar) == param_conv(i) ) exit
      end do
      if( i > nparam_conv ) exit boucle_ext
    end do
    ncal_nonfdm = ncal_nonfdm * nparam(igr)
  end do boucle_ext
  ngroup_par_conv = igr - 1

  do igr = 1,ngroup_par
    do ipar = 1,npar(igr)
      typeparc(igr,ipar) = typepar(igr,ipar)
      if( indice_par(igr,ipar) /= 0 ) then
        mot9 = typepar(igr,ipar)
        l = len_trim(mot9) + 1
        mot9(l:l) = '_'
        call ad_number(indice_par(igr,ipar),mot9,9)
        typeparc(igr,ipar) = mot9
      endif
    end do
  end do

  if( ngroup_par > 1 ) then
    allocate( indparp(ngroup_par) )
    indparp(:) = 0
  endif

  nnotskipm = 0
  do igr = 2,ngroup_par
    nnotskipm = nnotskipm + npar(igr)
  end do
  allocate( ifile_notskip( nnotskipm ) )

  n_atom_proto_p = 0

  bav_open = .false.

  if( Metric_cal .and. mpirank0 == 0 ) open(itape_minim, status = 'SCRATCH')

  do ical = 1,ncal

    if( ical == 1 ) then
      Conv_done = .false.
    else
      Conv_done = .true.
    endif

    ndem = 1
    inotskip = 0
    do igr = 2,ngroup_par

      j = ( ical - 1 ) / ndem
      indpar = mod( j, nparam(igr) ) + 1
      ndem = ndem * nparam(igr)

      if( nparam(igr) < 2 ) then
        param(igr,1:npar(igr)) = parmin(igr,1:npar(igr))
      else
        do ipar = 1,npar(igr)
          if( Length_block(igr) /= 0 .and. ipar > 1 ) exit
          param(igr,ipar) = parmin(igr,ipar) + ( indpar - 1 ) * ( parmax(igr,ipar) - parmin(igr,ipar) ) &
                          / ( nparam(igr) - 1 )
        end do
      endif

      if( Length_block(igr) > 0 ) then ! on est a la fin du block
        param_dep = parsum(igr) - sum( param(igr-Length_block(igr)+1:igr,1) )
        param(igr-Length_block(igr)+1:igr,2) = param_dep
      endif

      if( indpar /= indparp(igr) ) then

        do ipar = 1,npar(igr)
          if( typepar(igr,ipar) /= 'shift' .and. typepar(igr,ipar) /= 'weight' .and. &
              typepar(igr,ipar) /= 'gaussian' ) Conv_done = .false.
        end do

        if( ical > 1 ) then
          do ipar = 1,npar(igr)
            if( typepar(igr,ipar) /= 'dposx' .and. typepar(igr,ipar) /= 'dposy' .and. typepar(igr,ipar) /= 'dposz' .and. &
                typepar(igr,ipar) /= 'posx' .and. typepar(igr,ipar) /= 'posy' .and. typepar(igr,ipar) /= 'posz' .and. &
                typepar(igr,ipar) /= 'poporb' ) cycle
            inotskip = inotskip + 1
            ifile_notskip(inotskip) = indice_par(igr,ipar)
          end do
        endif

      endif
      indparp(igr) = indpar
    end do
    nnotskip = inotskip

    if( ncal_nonfdm == 1 ) then
      ifdm = 1
    else
      ifdm = mod(ical,ncal_nonfdm)
    endif

    if( Fdmnes_cal .and. ifdm == 1 ) then
      if( ical > 1 ) Close(3)

      call fdm(Ang_borm,Bormann,comt,Convolution_cal,Delta_edge,E_cut_imp,E_cut_man,Ecent,Elarg,Estart,Fit_cal, &
          Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,hkl_borm,icheck,ifile_notskip,indice_par,iscratch, &
          itape1,itape4,mpinodes0,mpirank0,n_atom_proto_p,ngamh,ngroup_par,nnotskip,nnotskipm, &
          nomfich,nomfichbav,npar,nparm,param,Scan_a,typepar,Use_FDMX,FDMX_only, &
          fdmnes_inp,cm2g,nobg,nohole,nodw,noimfp,imfp_inp,imfp_infile,elf_inp,elf_infile,dwfactor_inp,dwfactor,tdebye_inp, &
          tdebye,tmeas_inp,tmeas,expntl,expntlA,expntlB,victoreen,victA,victB,mermrank)

      if( sum(icheck(1:27)) > 0 ) then
        bav_open = .true.
      else
        bav_open = .false.
      endif
    endif

    if( mpirank0 /= 0 ) cycle

    if( Convolution_cal ) call convolution(bav_open,Bormann,Conv_done, &
        convolution_out,Delta_edge,E_cut_imp,E_cut_man,Ecent,Elarg,Estart,Fit_cal,Gamma_hole,Gamma_hole_imp,Gamma_max, &
        ical,icheck(30),indice_par,iscratchconv, itape1,kw_conv,length_line, &
        ngamh,ngroup_par,nkw_conv,nomfich,nomfichbav, npar,nparm,param,Scan_a,typepar,ncal)

    if( Metric_cal ) then

      if( ical == 1 ) then
        if( nomfich == 'fdmnes_out' ) then
          fdmfit_out = 'fdmfit_out.txt'
        else
          mot = nomfich
          l = len_trim(mot)
          mot(l+1:l+8) = '_fit.txt'
          fdmfit_out = mot
        endif
      endif

      call metric(comt,convolution_out,Dafs_bio,Dist_min,Dist_min_g,fdmfit_out,Fit_cal,Gen_Shift_min,ical, &
             ical_Met_min,index_Met_Fit,iscratchconv,itape_minim,itape2,length_line,Length_block,nb_datafile,ncal,ndm, &
             ng,ngroup_par,nmetric,nmetricm,Nom_Met,npar,nparam,nparm,param,parmax,parmin, RapIntegrT_min_g,typeparc)

    endif

  end do

  if( ngroup_par > 1 ) deallocate( indparp )

  allocate( par_op(ngroup_par) )

  if( Metric_cal .and. mpirank0 == 0 ) then

    npm = nparam(1)
    do igr = 2,ngroup_par
      npm = max( npm, nparam(igr) )
    end do

    allocate( par(ngroup_par,npm) )
    do igr = 1,ngroup_par
      typeparg(igr) = typeparc(igr,1)
      do ip = 1,nparam(igr)
        if( nparam(igr) < 2 ) then
          par(igr,ip) = parmin(igr,1)
        else
          fac = ( parmax(igr,1) - parmin(igr,1) ) / ( nparam(igr) - 1 )
          par(igr,ip) = parmin(igr,1) +  fac * ( ip - 1 )
        endif
      end do
    end do

    call minim(fdmfit_out,index_Met_Fit,itape_minim,Minim_fdm_ok,minimok,ncal,ndm,ngroup_par,ngroup_par_conv,nmetric, &
               nmetricm,Nom_Met,nparam,npm,par,par_op,typeparg)

    deallocate( par )

  endif

  if( Fit_cal .and. mpinodes0 > 1 ) then
    call MPI_Bcast(Minim_fdm_ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(minimok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)

  endif

  if( Fit_cal .and. minimok .and. ncal > 1 ) then

    if( mpinodes0 > 1 ) call MPI_Bcast(par_op,ngroup_par, MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    ndem = 1
    inotskip = 0
    do igr = 2,ngroup_par

      j = ( ical - 1 ) / ndem
      indpar = mod( j, nparam(igr) ) + 1
      ndem = ndem * nparam(igr)

      if( nparam(igr) > 1 ) then
        prop = ( par_op(igr) - parmin(igr,1) ) / ( parmax(igr,1) - parmin(igr,1) )
        do ipar = 1,npar(igr)
          if( Length_block(igr) /= 0 .and. ipar > 1 ) exit
          param(igr,ipar) = ( 1 - prop ) * parmin(igr,ipar) + prop * parmax(igr,ipar)
        end do
      endif

      if( Length_block(igr) > 0 ) then ! on est a la fin du block
        param_dep = parsum(igr) - sum( param(igr-Length_block(igr)+1:igr,1) )
        param(igr-Length_block(igr)+1:igr,2) = param_dep
      endif

      do ipar = 1,npar(igr)
        if( typepar(igr,ipar) /= 'shift' .and. typepar(igr,ipar) /= 'weight' .and. &
            typepar(igr,ipar) /= 'gaussian' ) Conv_done = .false.
      end do

      do ipar = 1,npar(igr)
        if( typepar(igr,ipar) /= 'dposx' .and. typepar(igr,ipar) /= 'dposy' .and. typepar(igr,ipar) /= 'dposz' .and. &
            typepar(igr,ipar) /= 'posx' .and. typepar(igr,ipar) /= 'posy' .and. typepar(igr,ipar) /= 'posz' .and. &
            typepar(igr,ipar) /= 'poporb') cycle
        inotskip = inotskip + 1
        ifile_notskip(inotskip) = indice_par(igr,ipar)
      end do

    end do
    nnotskip = inotskip

    if( Fdmnes_cal .and. Minim_fdm_ok .and. ncal /= ncal_nonfdm ) then
      Close(3)
      call fdm(Ang_borm,Bormann,comt,Convolution_cal,Delta_edge,E_cut_imp,E_cut_man,Ecent,Elarg,Estart,Fit_cal, &
        Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,hkl_borm,icheck,ifile_notskip,indice_par,iscratch, &
        itape1,itape4,mpinodes0,mpirank0,n_atom_proto_p,ngamh,ngroup_par,nnotskip,nnotskipm, &
        nomfich,nomfichbav,npar,nparm,param,Scan_a,typepar,Use_FDMX,FDMX_only, &
        fdmnes_inp,cm2g,nobg,nohole,nodw,noimfp,imfp_inp,imfp_infile,elf_inp,elf_infile,dwfactor_inp,dwfactor,tdebye_inp, &
        tdebye,tmeas_inp,tmeas,expntl,expntlA,expntlB,victoreen,victA,victB,mermrank)

      if( sum(icheck(1:27)) > 0 ) then
        bav_open = .true.
      else
        bav_open = .false.
      endif
    endif

    if( mpirank0 == 0 ) then
      if( Convolution_cal ) call convolution(bav_open,Bormann, .false., &
        convolution_out,Delta_edge,E_cut_imp,E_cut_man,Ecent, Elarg,Estart,Fit_cal,Gamma_hole,Gamma_hole_imp,Gamma_max, &
        ical,icheck(30),indice_par,iscratchconv, itape1,kw_conv,length_line, &
        ngamh,ngroup_par,nkw_conv,nomfich,nomfichbav,npar,nparm,param,Scan_a,typepar,ncal)

      call metric(comt,convolution_out,Dafs_bio,Dist_min, Dist_min_g,fdmfit_out,Fit_cal,Gen_Shift_min,ical, &
             ical_Met_min,index_Met_Fit,iscratchconv,itape_minim,itape2,length_line,Length_block,nb_datafile,ncal,ndm, &
             ng,ngroup_par,nmetric,nmetricm,Nom_Met,npar,nparam,nparm,param,parmax,parmin, RapIntegrT_min_g,typeparc)
    endif

  endif

  if( Metric_cal .and. mpirank0 == 0 ) close(itape_minim)

  deallocate( par_op )
  deallocate( ifile_notskip )

  deallocate( block_sum )
  deallocate( npbl )
  deallocate( Length_block )
  deallocate( npar )
  deallocate( param )
  deallocate( parmax )
  deallocate( parmin )
  deallocate( parsum )
  deallocate( nparam )
  deallocate( indice_par )
  deallocate( typepar )
  deallocate( typeparc )
  deallocate( typeparg )

  if( mpirank0 /= 0 ) return

  deallocate( Dist_min_g )
  deallocate( RapIntegrT_min_g )

! Selection of DAFS scan
  if( selec_cal ) call selec(itape5)

  if( Convolution_cal ) Close(itape1)
  if( Metric_cal ) Close(itape2)
  if( Fdmnes_cal ) Close(itape4)
  if( selec_cal ) Close(itape5)

  if( bav_open ) then
    call CPU_TIME(time)
    tp_fin = real(time,db)
    tpt = tp_fin - tp_deb
    itpj = int( tpt / 86400 )
    itph = int( tpt / 3600 )
    rtph = tpt - itph * 3600
    itpm = int( rtph / 60 )
    itps = nint( rtph - itpm * 60 )
    write(3,130) tpt
    if( itph > 0 .or. itpm > 0 ) write(3,140) itph, itpm, itps
    if( itpj > 1 ) write(3,150) itpj, itph-24*itpj, itpm, itps
    write(3,160)
    Close(3)
  endif

  return
  110 format(///' The parameter called Par_',a9,/ &
         ' needs a complete calculation. Not only the convolution part !',/ &
         ' Check your indata file under the keyword parameter, or',/ &
         ' rewrite your indata file with a complete calculation !'//)
  120 format(///' The parameter called Par_',a9,' does not exist !'/, &
                ' Check your indata file under the keyword parameter !'//)
  130 format(/' ',121('-'),//' Total time =',f10.1,' sCPU')
  140 format('            =',i4,' h,',i3,' min,',i3,' sCPU')
  150 format('            =',i4,' d,',i3,' h,',i3,' min,',i3,' sCPU')
  160 format(/'    Have a beautiful day !')

end

!*********************************************************************

function identmot(mot,longueur)

  use declarations
  implicit none

  integer:: i, idebut, ifin, l, longueur

  character(len=1) let(52)
  character(len=132) identmot, mot

  data let/'a','b','c','d','e','f','g','h','i','j','k','l','m', &
   'n','o','p','q','r','s','t','u','v','w','x','y','z','A','B','C', &
   'D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

  do i = 1,131-longueur
    if( mot(i:i) == char(9) ) then
      mot(i:i) = ' '
      cycle
    endif
  end do

  do i = 1,131-longueur
    if(mot(i:i) /= ' ') exit
  end do
  idebut = i

  do i = idebut+1,idebut+longueur
    if(mot(i:i) == ' ') exit
  end do
  ifin = i-1
  identmot = mot(idebut:ifin)

! Convertion from capital to lower case letter
  do i = 1,longueur
    do l = 27,52
      if( identmot(i:i) == let(l) ) then
        identmot(i:i) = let(l-26)
        exit
      endif
    end do
  end do

  return
end

!*********************************************************************

function Traduction(keyword)

  use declarations
  implicit none

  character(len=9) keyword, traduction

  traduction = keyword

  select case(keyword)
    case('iabsorbeu','absorbor','assorbito','absorber')
      traduction = 'absorbeur'
    case('interpoin','inter_poi')
      traduction = 'adimp'
    case('allnrixs','allnrix','allxraman','all_xrama')
      traduction = 'all_nrixs'
    case('all_site')
      traduction = 'allsite'
    case('angspin','spin_ang')
      traduction = 'ang_spin'
    case('atome','atomes','atoms','atomo','atomi')
      traduction = 'atom'
    case('atome_ch','atomch','atomech','atomoch','atomo_ch', 'atomich','atomi_ch')
      traduction = 'atom_ch'
    case('axespin','spin_axe','spin_axes','spinaxe','spinaxes', 'spin_axis','spinaxis')
      traduction = 'axe_spin'
    case('basecomp')
      traduction = 'base_comp'
    case('basereel','real_basi','real_base')
      traduction = 'base_reel'
    case('borman')
      traduction = 'bormann'
    case('cartesien')
      traduction = 'cartesian'
    case('chlibre','freech')
      traduction = 'chlib'
    case('cristal_c','crystal_c','ciffile')
      traduction = 'cif_file'
    case('icheck')
      traduction = 'check'
    case('centre')
      traduction = 'center'
    case('centre_ab','centreabs','centerabs')
      traduction = 'center_ab'
    case('checkall','allcheck','all_check')
      traduction = 'check_all'
    case('checkcoab')
      traduction = 'check_coa'
    case('checkmat')
      traduction = 'check_mat'
    case('checkpot')
      traduction = 'check_pot'
    case('checksph','checksphe')
      traduction = 'check_sph'
    case('spinresol','spin_reso','coreresol')
      traduction = 'core_reso'
    case('crist','cryst','cristallo','cristal')
      traduction = 'crystal'
    case('dmaxpot','dmax_pot','d_maxpot','distmaxpo','dist_maxp')
      traduction = 'd_max_pot'
    case('dafsscan','dafscan','rxs','rxd')
      traduction = 'dafs'
     case('rxs_bio','rxd_bio','rxsbio','rxdbio','dafsbio','rxs_exp', 'rxd_exp')
      traduction = 'dafs_exp'
    case('rxs_2d','rxd_2d','dafs_surf','rxs_surfa','rxd_surfa')
      traduction = 'dafs_2d'
    case('state_den','statedens','densite')
      traduction = 'density'
    case('stateall','state_all','densityal')
      traduction = 'density_a'
    case('converge','deltae','delta_e_c','delta_en','delta_e_s')
      traduction = 'delta_en_'
    case('dilat','dilatati','dilat_or')
      traduction = 'dilatorb'
   case('dip_relat')
      traduction = 'dip_rel'
    case('magdip','dip_mag','mag_dip')
      traduction = 'dipmag'
    case('magquad','dip_quad')
      traduction = 'dipquad'
    case('dynamical')
      traduction = 'dyn_g'
    case('dipole')
      traduction = 'e1e1'
    case('e3e1')
      traduction = 'e1e3'
    case('seuil','threshold','soglia')
      traduction = 'edge'
    case('potimag  ')
      traduction = 'eimag'
    case('fin','fine')
      traduction = 'end'
    case('end_jump')
      traduction = 'endjump'
    case('eclie')
      traduction = 'etatlie'
    case('extract_s')
      traduction = 'extracsy'
    case('enrgpsii','epsiia')
      traduction = 'epsii'
    case('enrgphot','energpho')
      traduction = 'energphot'
    case('extractio')
      traduction = 'extrac'
    case('lapw','wien')
      traduction = 'flapw'
    case('lapw_s','wien_s','flapw_sau','lapw_sauv','wien_sauv', 'flapw_sav','lapw_save','wien_save')
      traduction = 'flapw_s'
    case('lapw_r','wien_r','flapw_rec','lapw_recu','wien_recu')
      traduction = 'flapw_r'
    case('fprim_ato','fprimatom','fprimeato')
      traduction = 'fprime_at'
    case('hedin_lun','hedinlund')
      traduction = 'hedin'
    case('helm','helmoltz','helmhotz','helmotz')
      traduction = 'helmholtz'
    case('hubard')
      traduction = 'hubbard'
    case('lmaxat','lmaxat0','lmax_at','lmax_atom')
      traduction = 'lmax'
    case('lmaxnrixs','lmax_xram','lmaxxrama')
      traduction = 'lmax_nrix'
    case('lamstdens')
      traduction = 'lmaxstden'
    case('lmaxso','lmaxso0','lmax_so','lmax_sor','lmaxsort')
      traduction = 'lmaxso'
    case('lminus1')
      traduction = 'lmoins1'
    case('lquadimp')
      traduction = 'lquaimp'
    case('m1_m2','m2m1','m2_m1')
      traduction = 'm1m2'
    case('m1_m1')
      traduction = 'm1m1'
    case('magnetiqu','magnetic')
      traduction = 'magnetism'
    case('biologie','biology','save_memo','savememor','memorysav')
      traduction = 'memory_sa'
    case('molec','molecola')
      traduction = 'molecule'
    case('muffin_ti')
      traduction = 'muffintin'
    case('n_scf','nself','nscf')
      traduction = 'n_self'
    case('chemin')
      traduction = 'nchemin'
    case('notdipole','nondipole','nodipole','no_dipole','noe1e1')
      traduction = 'no_e1e1'
    case('notinterf','noninterf','nointerf','no_dipqua','nodipquad', 'nondipqua','no_interf','noe1e2')
      traduction = 'no_e1e2'
    case('noe1e3','noe3e1','no_e3e1')
      traduction = 'no_e1e3'
    case('noe2e2')
      traduction = 'no_e2e2'
    case('not_core_','not_corer','notcorere','notcorer_','no_core_r','no_corere','nocoreres','nocore_re')
      traduction = 'no_core_r'
    case('nofermi')
      traduction = 'no_fermi'
    case('nonrelat','nonrelati')
      traduction = 'non_relat'
    case('non_exc','nonexcite','non_excit','no_hole','nohole')
      traduction = 'nonexc'
    case('normal_ta')
      traduction = 'normaltau'
    case('noteneg')
      traduction = 'not_eneg'
    case('nrato_dir','nratodira','nr_ato')
      traduction = 'nrato'
    case('nixs','xraman','xramman')
      traduction = 'nrixs'
    case('dipole_oc','dip_oct')
      traduction = 'octupole'
    case('oldzero')
      traduction = 'old_zero'
    case('overad','overrad')
      traduction = 'over_rad'
    case('optics','optique')
      traduction = 'optic'
    case('p_self0','pself','pself0','p0_self','p_scf','p0_scf', 'pscf','pscf0','p0scf')
      traduction = 'p_self'
    case('p_scf_max','pselfmax','pscfmax','ppscf_max','p_scfmax','p_selfmax', 'pself_max')
      traduction = 'p_self_max'
    case('cristal_p','crystal_p','pdbfile')
      traduction = 'pdb_file'
    case('polarise','polarised','polarisat','polarizat','polar', 'polarize')
      traduction = 'polarized'
    case('quadripol')
      traduction = 'quadrupol'
    case('r_selfcon','rself','rselfcons','r_scf','rscf')
      traduction = 'r_self'
    case('rayon','rsort','rsorte','raggio')
      traduction = 'radius'
    case('gamme')
      traduction = 'range'
    case('gammel')
      traduction = 'rangel'
    case('raychimp')
      traduction = 'rchimp'
    case('relat')
      traduction = 'relativis'
    case('rmtimp','rmt_imp','rmtgimp','rmtg_imp')
      traduction = 'rmtg'
    case('rmtvo')
      traduction = 'rmtv0'
    case('ang_rotsu','angrotsup','rotsup')
      traduction = 'rot_sup'
    case('rpa','rpa_lf')
      traduction = 'rpalf'
    case('selfcons','self_cons')
      traduction = 'scf'
    case('selfexc','self_exc','scfexc')
      traduction = 'scf_exc'
    case('selfnonex','self_non_','self_none','scfnonexc', 'scf_nonex','scf_nohol','scf_no_ho')
      traduction = 'scf_non_e'
    case('step_scf','pas_scf')
      traduction = 'scf_step'
    case('ecrantage')
      traduction = 'screening'
    case('singsol','solsing_o')
      traduction = 'solsing'
    case('spin_orbi')
      traduction = 'spinorbit'
    case('stepazim','stepazimu','azim_step')
      traduction = 'step_azim'
    case('selfabsor','self_abso','selfabs')
      traduction = 'self_abs'
    case('sym','point_gro')
      traduction = 'pointgrou'
    case('sym_mol')
      traduction = 'symmol'
    case('fileout','file_out')
      traduction = 'filout'
    case('sphere_si')
      traduction = 'sphere_al'
    case('xcorr','tdlda','td_lda','td_dft')
      traduction = 'tddft'
    case('testdistm','dist_min','distmin','distamin')
      traduction = 'test_dist'
    case('voimp','v0bdcfim','vmoyf','korigimp')
      traduction = 'v0imp'
    case('v_intmax','v_max')
      traduction = 'vmax'
    case('xalfa','alfpot')
      traduction = 'xalpha'
    case('ylm_compl','ylmcomp','ylmcomple','harmo_com','harmocomp')
      traduction = 'ylm_comp'
    case('zabsorber','zabsorbeu','numatabs','numat_abs')
      traduction = 'z_absorbe'

! Convolution
    case('conv_all')
      traduction = 'all_conv'
    case('arc')
      traduction = 'convoluti'
    case('resolutio','gaus','gauss')
      traduction = 'gaussian'
    case('fprim')
      traduction = 'fprime'
    case('seah_denc')
      traduction = 'seah'
    case('noextrap')
      traduction = 'no_extrap'
    case('photoemis','photo')
      traduction = 'photo_emi'
    case('s02','so2','so_2')
      traduction = 's0_2'
    case('elarge','e_large','ewidth','e_width')
      traduction = 'elarg'
    case('e_cent','ecenter','ecentre','e_center','e_centre')
      traduction = 'ecent'
    case('scan_out','scanout','scanconv')
      traduction = 'scan_conv'
    case('select_co','selectcor','seleccore','sel_core','selcore')
      traduction = 'selec_cor'
    case('stoke')
      traduction = 'stokes   '
    case('stokesnam','stokename','stoke_nam')
      traduction = 'stokes_na'
    case('circulair')
      traduction = 'circular'
    case('no_analys','noanalyzo','noanalyze','noanalyso','noanalyse')
      traduction = 'no_analyz'

! Fit
    case('pop_orb')
      traduction = 'poporb'
    case('fit_rx')
      traduction = 'rx'

! Conv_gaus
    case('filein','filin')
      traduction = 'file_in'

! Selec
    case('reflexion')
      traduction = 'reflectio'
    case('select_in','selec_in','selecin','selecinp','selec_ind')
      traduction = 'selec_inp'
    case('select_ou','selecout','selec_ou')
      traduction = 'selec_out'

! General
    case('ecut','e_cut_imp','ecut_imp','ecutimp','e_cutimp','e_fermi','efermi')
      traduction = 'e_cut'
    case('e_start')
      traduction = 'estart'

  end select

  return
end


!***********************************************************************

! Fonction giving the number of number in the next non empty line et setting at the biginning of this line.
! With line starting with a character, nnombre = 0

function nnombre(irec,length_line)

  use declarations
  implicit none

  integer, parameter:: nlet = 53

  integer:: eof, i, icol, irec, j, k, l, length_line, ligne, nmots, n, nnombre

  character(len=length_line):: mot, test
  character(len=1), dimension(nlet):: let

  real(kind=db):: x

  data let /'d','e','D','E','a','A','b','B','c','C','f','F','g','G','h','H','i','I','j','J', &
            'k','K','l','L','m','M','n','N', 'o','O','p','P','q','Q','r','R','s','T','t','T','u','U','v','V', &
            'w','W','x','X','y','Y','z','Z','/'/

  do i = 1,length_line
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

    do i = 1,length_line
      if( mot(i:i) /= ' ' .and. mot(i:i) /= char(9) ) exit boucle_ligne
    end do

  end do boucle_ligne

  backspace(irec)

  Open( 16, status='SCRATCH' )

  do i = 1,length_line
    if( mot(i:i) == char(9) ) mot(i:i) = ' '
  end do

  n = 0
  i = 0
  do icol = 1,length_line
    i = i + 1
    if( i > length_line ) exit
    if( mot(i:i) == ' ' ) cycle
    do j = i+1,length_line
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

subroutine mult_cell(itape,File_out)

  use declarations
  implicit none

  integer, parameter:: nam = 100000
  integer, parameter:: ntypem = 100

  integer:: eof, eoff, i, igrdat, ipr, itape, ix, iy, iz, j, l, n, nnombre, ntype

  real(kind=db):: alfa, ax, ay, az, beta, gamma, Thickness
  real(kind=db), dimension(3):: q
  real(kind=db), dimension(3,nam):: p

  integer, dimension(3):: nm
  integer, dimension(nam):: itype, Z
  integer, dimension(ntypem):: Numat

  Logical:: Ang, Surface, Typ

  character(len=2):: Chemical_Symbol
  character(len=9):: keyword
  character(len=132):: File_out, identmot, mot

  Ang = .false.
  Surface = .false.
  Typ = .false.
  Thickness = 0._db
  Numat(:) = 0
  nm(:) = 1

  mot = File_out
  l = len_trim(mot)
  if( mot(l-3:l) /= '.txt' ) then
    mot(l+1:l+4) = '.txt'
    File_out = mot
  endif
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
          read(itape,*) ax, ay, az
        elseif( n >= 6 ) then
          Ang = .true.
          read(itape,*) ax, ay, az, alfa, beta, gamma
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

        n = i-1

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

  if( Surface ) then
    if( Ang ) then
      write(2,120) ax*nm(1), ay*nm(2), az, alfa, beta, gamma
      write(6,120) ax*nm(1), ay*nm(2), az, alfa, beta, gamma
    else
      write(2,120) ax*nm(1), ay*nm(2), az
      write(6,120) ax*nm(1), ay*nm(2), az
    endif
  else
    if( Ang ) then
      write(2,120) ax*nm(1), ay*nm(2), az*nm(3), alfa, beta, gamma
      write(6,120) ax*nm(1), ay*nm(2), az*nm(3), alfa, beta, gamma
    else
      write(2,120) ax*nm(1), ay*nm(2), az*nm(3)
      write(6,120) ax*nm(1), ay*nm(2), az*nm(3)
    endif
  endif

  if( Surface ) nm(3) = int( Thickness / az ) + 1
  j = 0
  do ix = 0,nm(1)-1
    do iy = 0,nm(2)-1
      do iz = 0,nm(3)-1
        write(2,*)
        do i = 1,n
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

  Close(itape)
  Close(2)

  return

  110 format(//' number of type =',i5,' > ntypem =',i4,// ' Change the parameter ntypem in the code !'//)
  120 format(/' Crystal',/5x,3f14.10,3f12.5)
  130 format(//' number of atoms =',i6,' > nam =',i6,// ' Change the parameter nam in the code !'//)
  140 format(//'  Error in the indata file :')
  150 format(//' The following line is not understood :',/A,// &
          ' If it is a keyword, check the spelling.'/, &
          ' If the line is not supposed to be a keyword but contains numbers, check:'/ &
          5x,' - How many numbers must be in the line ?'/, &
          5x,' - Are there spaces between the numbers ?'/, &
          5x,' - Tabulations are forbidden !'//)
  160 format(i5,3f14.10,'  ! ',i5,3x,a2)
end

