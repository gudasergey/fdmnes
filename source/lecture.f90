! FDMNES subroutines related to the indata files reading

!-------------------------------------------------------------------------------------------------------------------------------

! Reading subroutine giving the dimensions of the tables, necessary for the rest of the code, including the main reading routine 'lecture'

subroutine lectdim(Absauto,Atom_occ_hubb,Atom_nonsph,Axe_loc,Bormann,Bulk,Cap_layer,Dafs_bio,Doping,Extract,Extract_ten,Film, &
    Flapw,Full_self_abs,Hubbard,itape4,Length_line,Magnetic,Memory_save,mpinodes0,mpirank0,n_adimp,n_atom,n_atom_coop,  &
    n_file_dafs_exp,n_mat_polar,n_multi_run_e,n_q_rixs,n_radius,n_range,n_theta_in,n_two_theta,n_Z_abs,nb_atom_conf_m,ncolm, &
    neimagent,nenerg_in_rixs,nenerg_s,ngamme, &
    ngroup,nhybm,nklapw,nlatm,nlmlapwm,nmatsym,norbdil,npldafs,npldafs_2d,npldafs_e,nple,nplrm,nq_nrixs, &
    NRIXS,nspin,nspino,nspinp,ntype,ntype_bulk,ntype_conf,Occupancy_first,Pdb,Readfast,Self_abs,Taux,Temp_B_iso,Use_FDMX,Xan_atom)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: eof, eoff, i, iabsm_1, ie, ier, igamme, igc, igr, igrdat, io, ipr, ipl, istat, it, itape4, itype_dop, j, jgr, &
    jpl, k, l, Length, Length_line, lin_gam, lmax_nrixs, mpierr, mpinodes0, mpirank0, n, n_adimp, n_atom_bulk, n_atom_cap, &
    n_atom_coop, n_atom_int, n_atom_per, n_atom_per_neq, n_atom_sur, n_atom_uc, n_dic, n_file_dafs_exp, &
    n_fract_x, n_fract_y, n_fract_z, n_label, n_mat_polar, n_mu, n_multi_run_e, n_q_rixs, n_radius, n_range, n_symbol, &
    n_theta_in, n_two_theta, n_Z_abs, &
    nb, nb_atom_conf_m, ncolm, neimagent, nenerg_in_rixs, nenerg_s, ngamme, ngc, ngroup, n_atom_neq, nhybm, nklapw, nl, &
    nlatm, nlmlapwm, nmatsym, nn, nnombre, norbdil, norbv, npldafs, npldafs_2d, npldafs_e, npldafs_t, &
    nple, nplm, nplrm, nq_nrixs, nspin, nspino, nspinp, ntype, ntype_bulk, numat_abs, ntype_conf, Wien_save, Z

  integer, dimension(8):: n_atom
  integer, dimension(:), allocatable :: Z_bulk

  character(len=4):: mot4
  character(len=6):: mot6
  character(len=9):: keyword
  character(len=13):: Space_Group, Spgr
  character(len=132):: identmot, mot132, mot
  character(len=Length_name):: Fcif_file, Fichier, Fichier_pdb
  character(len=Length_name), dimension(9):: Wien_file

  logical:: Absauto, adimpin, Apostrophe, Atom_conf, Atom_conf_dop, Atom_nonsph, Atom_occ_hubb, Axe_loc, Bormann, Bulk, Fcif, &
     Cap_layer, Dafs_bio, Doping, Extract, Extract_ten, Film, Flapw, Full_self_abs, Hubbard, &
     Magnetic, Matper, Memory_save, NRIXS, Occupancy_first, Pdb, Pol_dafs_in, Quadrupole, Readfast, RIXS_case, Screening, &
     Self_abs, Taux, Temp_B_iso, Use_FDMX, Very_fast, Xan_atom

  real(kind=db):: Adimp, Angz, de, def, number_from_text, E, r, Rsorte_s, Step, Theta_max, Theta_min
  real(kind=db), dimension(3):: p
  real(kind=db), dimension(:), allocatable:: E_adimp, E_radius, Egamme
  real(kind=db), dimension(:,:), allocatable:: pop

  Absauto = .true.
  Angz = 0._db
  Atom_conf = .false.
  Atom_conf_dop = .false.
  Atom_nonsph = .false.
  Atom_occ_hubb = .false.
  Axe_loc = .false.
  Bulk = .false.
  Cap_layer = .false.
  Fcif = .false.
  Dafs_bio = .false.
  Doping = .false.
  Extract = .false.
  Extract_ten = .false.
  Film = .false.
  Flapw = .false.
  Full_self_abs = .false.
  Hubbard = .false.
  iabsm_1 = 1
  lmax_nrixs = 2
  Matper = .false.
  Memory_save = .false.
  n_adimp = 1
  n_atom_bulk = 0
  n_atom_cap = 0
  n_atom_coop = 0
  n_atom_int = 0
  n_atom_per = 0
  n_atom_per_neq = 0
  n_atom_sur = 0
  n_atom_uc = 0
  n_dic = 0
  n_file_dafs_exp = 0
  n_mat_polar = 0
  n_mu = 0
  n_multi_run_e = 1
  n_q_rixs = 0
  n_range = 1
  n_radius = 1
  n_theta_in = 0
  n_two_theta = 0
  n_Z_abs = 1
  nb_atom_conf_m = 0
  neimagent = 0
  nenerg_in_rixs = 0
  nenerg_s = 125
  ngamme = 3
  nhybm = 0
  nklapw = 1
  nlmlapwm = 1
  nlatm = 0
  nmatsym = 1
  norbdil = 0
  nple = 0
  npldafs = 0
  npldafs_2d = 0
  npldafs_e = 0
  nq_nrixs = 0
  NRIXS = .false.
  nspin = 1
  nspino = 1
  nspinp = 1 ! pour la diffusion ( = 2 si spinorbite meme potentiel non magnetique )
  ntype = 0
  ntype_bulk = 0
  ntype_conf = 0
  numat_abs = 0
  Occupancy_first = .true.
  Pdb = .false.
  Pol_dafs_in = .false.
  Quadrupole = .false.
  Readfast = .false.
  Screening = .false.
  Self_abs = .false.
  Space_Group = ' '
  Taux = .false.
  Temp_B_iso = .false.
  Xan_atom = .false.
  Wien_save = 0

  adimpin = .false.

  if( Bormann ) npldafs = 36

  if( mpirank0 == 0 ) then

    Rewind(itape4)

    boucle_gen: do

      read(itape4,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit boucle_gen

      keyword = identmot(mot,9)

      if( keyword(1:1) == '!' ) cycle

      select case(keyword)

        case('bulk')
          Bulk = .true.
          read(itape4,*,iostat=ier) p(:)
          do igr = 1,100000
            read(itape4,*,iostat=ier) i, p(:)
            if( ier /= 0 ) exit
          end do
          n_atom_bulk = igr - 1

          do igr = 1,n_atom_bulk+1
            Backspace(itape4)
          end do
          allocate( Z_bulk(n_atom_bulk) )
          do igr = 1,n_atom_bulk
            read(itape4,*,iostat=ier) Z_bulk(igr)
          end do
          ntype_bulk = 1
          do igr = 2,n_atom_bulk
            do jgr = 1,igr-1
              if( Z_bulk(igr) == Z_bulk(jgr) ) exit
            end do
            if( jgr >= igr ) ntype_bulk = ntype_bulk + 1
          end do

        case('cap_layer')
          Cap_layer = .true.
          read(itape4,*,iostat=ier) p(:)
          do igr = 1,100000
            read(itape4,*,iostat=ier) i, p(:)
            if( ier /= 0 ) exit
          end do
          Backspace(itape4)
          n_atom_cap = igr - 1

        case('coop')
          n_atom_coop = nnombre(itape4,132)

        case('doping')
          Doping = .true.
          read(itape4,*) itype_dop
          Absauto = .false.
          n_multi_run_e = 1

        case('absorbeur')
          Absauto = .false.
          n = nnombre(itape4,132)
          if( n > 0 ) n_multi_run_e = 0
          do i = 1,100000
            n = nnombre(itape4,132)
            if( n < 1 ) exit
            n_multi_run_e = n_multi_run_e + n
            if( i == 1 ) then
              read(itape4,*) iabsm_1
            else
              read(itape4,*)
            endif
          end do

        case('z_absorbe')
          n = nnombre(itape4,132)
          n_Z_abs = 0
          do i = 1,100000
            n = nnombre(itape4,132)
            if( n < 1 ) exit
            n_Z_abs = n_Z_abs + n
            if( i == 1 ) then
              read(itape4,*) numat_abs
            else
              read(itape4,*)
            endif
          end do
          n_Z_abs = max( 1, n_Z_abs )

        case('end')
          exit

        case('extract')
          Extract = .true.

        case('extract_t')
          Extract_ten = .true.
          Extract = .true.

        case('hubbard','hubbard_z')
          Hubbard = .true.

        case('lmax_nrix')
          read(itape4,*,iostat=ier) lmax_nrixs

        case('rixs')
          do i = 1,100000
            read(itape4,*,iostat=ier) r
            if( ier > 0 ) then
              backspace(itape4)
              exit
            endif
          end do
          n_q_rixs = i - 1

        case('nrixs','nrixs_mon')
          NRIXS = .true.
          nq_nrixs = 0
          do
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            if( keyword == 'nrixs_mon' ) then
              nq_nrixs = nq_nrixs + 1
            else
              nq_nrixs = nq_nrixs + n
            endif
            read(itape4,*)
          end do
          if( nq_nrixs == 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,'(/A/)') ' After keyword "NRIXS" there is no q values'
            end do
            stop
          endif

        case('spgroup')
          n = nnombre(itape4,132)
          read(itape4,'(A)') mot
          if( mot(1:1) == ' ' ) mot = adjustl(mot)
          Space_group = mot(1:13)

        case('atom_b_is')
          Temp_B_iso = .true.

        case('atom_u_is')
          Temp_B_iso = .true.

        case('occupancy')
          Taux = .true.
          Occupancy_first = .not. Temp_B_iso

        case('xan_atom')
          Xan_atom = .true.

        case('radius')
          n = nnombre(itape4,132)
          if( mod(n,2) == 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,115) n
            end do
            stop
          endif
          if( n > 1 ) then
            n_radius = ( n + 1 ) / 2
            allocate( E_radius( n_radius - 1  ) )
            read(itape4,*,iostat=ier) ( rsorte_s, E_radius(i), i = 1,n_radius-1 )
          endif

        case('adimp')
          n = nnombre(itape4,132)
          adimpin = .true.
          if( mod(n,2) == 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,117) n
            end do
            stop
          endif
          if( n > 1 ) then
            n_adimp = ( n + 1 ) / 2
            allocate( E_adimp( n_adimp - 1 ) )
            read(itape4,*,iostat=ier) ( adimp, E_adimp(i), i = 1,n_adimp-1 )
          endif

        case('range','rangel','rixs_ener')
          RIXS_case = keyword(1:4) == 'rixs'

          n = nnombre(itape4,Length_line)
          if( .not. RIXS_case ) ngamme = n

          if( mod(n,2) == 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,120) keyword, n
            end do
            stop
          endif
          if( keyword == 'rangel' ) then
            lin_gam = 1
            if( n /= 1 ) then
              ngamme = 3
              n = 3
            endif
          else
            lin_gam = 0
          endif
          allocate( Egamme(n) )
          read(itape4,*,iostat=ier) Egamme(1:n)
          if( ier > 0 ) call write_err_form(itape4,keyword)

          do igamme = 2,n,2
            if( Egamme(igamme) > eps6 ) cycle
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,130) Keyword
            end do
            stop
          end do

          E = Egamme(1)
          if( n == 1 ) then
            if( RIXS_case ) then
              nenerg_in_rixs = 1
            else
              nenerg_s = 1
            endif
          elseif( Egamme(3) <= Egamme(1) ) then
            if( RIXS_case ) then
              nenerg_in_rixs = 1
            else
              nenerg_s = 1
            endif
          elseif( lin_gam == 1 ) then
            def = 10 / rydb
            do ie = 2,10000000
              r = 1 + ( E/rydb ) / def
              r = max( r, 0.25_db )
              de = sqrt( r ) * Egamme(2)
              e = e + de
              if( e > Egamme(n) + eps10 ) then
                nenerg_s = ie - 1
                exit
              endif
            end do
          else
            ngc = 2
            do ie = 2,1000000
              E = E + Egamme(ngc)
              if( E > Egamme(n) - eps10 ) then
                if( RIXS_case ) then
                  nenerg_in_rixs = ie
                else
                  nenerg_s = ie
                endif
                exit
              else
                do igc = ngc, n - 3, 2
                  if( E >= Egamme(igc+1) - eps10 ) then
                    ngc = ngc + 2
                    E = min( Egamme(ngc-1), E )
                  else
                    exit
                  endif
                end do
              endif
            end do
          endif
          deallocate( Egamme )

        case('eimag')
          do ie = 1,100000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            read(itape4,*)
          end do
          neimagent = ie - 1

        case('quadrupol')
          Quadrupole = .true.

        case('e1e2')
          Quadrupole = .true.

        case('e2e2')
          Quadrupole = .true.

        case('spinorbit')
          nspinp = 2
          nspino = 2

        case('two_theta')
          n = nnombre(itape4,132)
          n = min(n,3)
          if( n == 1 ) then
            n_two_theta = 1
          elseif( n == 3 ) then
            read(itape4,*) Theta_min, step, Theta_max
            if( abs(Step) < eps10 ) then
              n_two_theta = 1
            else
              n_two_theta = nint( ( Theta_max - Theta_min ) / Step ) + 1
              n_two_theta = max( n_two_theta, 1)
            endif
          elseif( n /= 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,'(//A//)') ' Under keyword "Two_theta" only 1 or 3 numbers (2theta_min, step, 2theta_max) are possible !'
            end do
            stop
          endif

        case('theta_in')
          n = nnombre(itape4,132)
          n = min(n,3)
          if( n == 1 ) then
            n_theta_in = 1
          elseif( n == 3 ) then
            read(itape4,*) Theta_min, step, Theta_max
            if( abs(Step) < eps10 ) then
              n_theta_in = 1
            else
              n_theta_in = nint( ( Theta_max - Theta_min ) / Step ) + 1
              n_theta_in = max( n_theta_in, 1)
            endif
          elseif( n /= 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,'(//A//)') ' Under keyword "Theta_in" only 1 or 3 numbers (theta_min, step, theta_max) are possible !'
            end do
            stop
          endif

        case('magnetism')
          nspinp = 2

        case('dilatorb')
          do i = 1,100000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            norbdil = norbdil + 1
            read(itape4,*)
          end do

        case('polarized')
          do jpl = 1,100000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            read(itape4,*,iostat=ier) p(:)
            if( ier > 0 ) call write_err_form(itape4,keyword)
            if( sum( p(:) )**2 < eps10 ) n_dic = n_dic + 1
          end do
          nple = jpl - 1

        case('mat_polar')
          do jpl = 1,100000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            read(itape4,*,iostat=ier) p(:), p(:)
            if( ier > 0 ) call write_err_form(itape4,keyword)
            n_mat_polar = n_mat_polar + 1
          end do

        case('dafs','dafs_2d')
          do ipl = 1,100000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            select case(n)
              case(3)
                read(itape4,*)
                nn = nnombre(itape4,132)
                select case(nn)
                  case(2,3,4,5)
                    read(itape4,*)
                  case(6,9)
                    Pol_dafs_in = .true.
                    read(itape4,*)
                    read(itape4,*)
                  case default
                    call write_error
                    do ipr = 6,9,3
                      write(ipr,140) ipl
                    end do
                    stop
                end select
              case(5,6,7,8,9,10,11,12,13)
                read(itape4,*)
              case default
                call write_error
                do ipr = 6,9,3
                   write(ipr,140) ipl
                end do
                stop
            end select
          end do
          npldafs = ipl - 1
          if( keyword == 'dafs_2d' ) npldafs_2d = npldafs

        case('dafs_exp')
          Dafs_bio = .true.
          do i = 1,3
            read(itape4,*)
          end do
          npldafs = 0
          do i = 1,1000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            read(itape4,'(/A)') Fichier
            Fichier = Adjustl(Fichier)
            Length = len_trim(Fichier)
            if( Length > 4 ) then
              if( Fichier(Length-3:Length-3) /= '.' ) Fichier(Length+1:Length+4) = '.txt'
            endif
            open(99, file=Fichier, status='old',iostat=istat)
            if( istat /= 0 ) call write_open_error(Fichier,istat,1)
            n = nnombre(99,Length_line)
            Close(99)
            npldafs = npldafs + n / 3
          end do
          n_file_dafs_exp = i - 1

        case('readfast')
          Readfast = .true.

        case('self_abs')
          Self_abs = .true.

        case('full_self')
          Full_self_abs = .true.

        case('atom_conf')

          Atom_conf = .true.
          do it = 1,100000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            ntype_conf = ntype_conf + 1
            read(itape4,*) nb, nn
            if( nn == 0 ) Atom_conf_dop = .true.
            nb_atom_conf_m = max( nb, nb_atom_conf_m )
            backspace(itape4)
            read(itape4,*,iostat=ier) ( nl, i = 1,nb+2 )
            if( ier > 0 ) call write_err_form(itape4,keyword)
            nlatm = max( nlatm, nl )
            if( n == nb + 2 + 4*nl .and. nl > 0 .and. nspin == 1 ) then
              backspace(itape4)
              allocate( pop(nl,2) )
              read(itape4,*,iostat=ier) ( i, l = 1,nb+1 ), nl, ( i, i, pop(l,:), l = 1,nl )
              if( ier > 0 ) call write_err_form(itape4,keyword)
              do l = 1,nl
                if( abs( pop(l,1) - pop(l,2) ) < eps10 ) cycle
                nspin = 2
                exit
              end do
              deallocate( pop )
            endif
          end do

        case('atom_nsph')
          do
            n = nnombre(itape4,132)
            if( n == 0 ) then
              exit
            elseif( n < 3 ) then
              call write_error
              do ipr = 6,9,3
                write(ipr,135)
              end do
              stop
            endif
            read(itape4,*) ( norbv, i = 1,n )
            nhybm = max (nhybm, norbv )
            do io = 1,norbv
              read(itape4,*)
            end do
          end do
          if( nhybm > 0 ) Atom_nonsph = .true.

        case('atom')

          do it = 1,100000
            n = nnombre(itape4,132)
            select case(n)
              case(0)
                exit
              case(1)
                ntype = ntype + 1
                read(itape4,*)
              case(3)
                cycle
              case default
                ntype = ntype + 1
                read(itape4,*) Z, nl
                nlatm = max( nlatm, nl )
                if( n == 2 + 4*nl .and. nl > 0 .and. nspin == 1 ) then
                  backspace(itape4)
                  allocate( pop(nl,2) )
                  read(itape4,*) Z, nl, ( i, i, pop(l,:), l = 1,nl )
                  do l = 1,nl
                    if( abs( pop(l,1) - pop(l,2) ) < eps10 ) cycle
                    nspin = 2
                    exit
                  end do
                  deallocate( pop )
                endif
            end select
          end do

! Cluster, film or unit cell description
        case('crystal','molecule','crystal_t','molecule_','film','film_t','surface','surface_t','interface','interfac_')
          if( keyword(6:6) == 't' .or.keyword(9:9) == 't' .or. keyword(9:9) == '_' ) Taux = .true.

          if( keyword(1:8) /= 'molecule' ) Matper = .true.
          if( keyword(1:1) == 'f' .or. keyword(1:1) == 's' .or. keyword(1:1) == 'i' ) Film = .true.

          n = nnombre(itape4,132)
          if( n == 0 ) then
            call write_error
            read(itape4,'(A)') mot132
            do ipr = 6,9,3
              write(ipr,150) keyword
              write(ipr,'(A)') mot132
              write(ipr,160)
            end do
            stop
          endif
          if( n > 5 ) then
            read(itape4,*,iostat=ier) ( Angz, i = 1,6 )
          else
            read(itape4,*)
          endif

          if( Readfast .or. Taux .or. Temp_B_iso ) then
            do igr = 1,100000
              read(itape4,*,iostat=ier) i, p(:)
              if( ier /= 0 ) exit
            end do
            Backspace(itape4)
          else
            do igr = 1,100000
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              if( n == 2 .or. n == 3 ) then
                Axe_loc = .true.
                read(itape4,*)
                n = nnombre(itape4,132)
              endif
              norbv = 0
              select case(n)
                case(4)
                  read(itape4,*)
                case(5)
                  read(itape4,*,iostat=ier) i, p(:), norbv
                  if( ier > 0 ) call write_err_form(itape4,keyword)
                  if( norbv < 0 ) then
                    nhybm = max( nhybm, - norbv - 1 )
                  else
                    nhybm = max( nhybm, norbv )
                  endif
                case default
                  call write_error
                  read(itape4,'(A)') mot132
                  do ipr = 6,9,3
                    write(ipr,150) keyword
                    write(ipr,'(A)') mot132
                    write(ipr,160)
                  end do
                  stop
              end select
              if( norbv == 0 ) cycle
              do io = 1,abs(norbv)
                read(itape4,*)
              end do
              if( norbv /= -1 ) Atom_nonsph = .true.
              if( norbv < 0 ) then
                Atom_occ_hubb = .true.
                Hubbard = .true.
              endif
            end do
          endif
          if( keyword(1:1) == 'c' .or. keyword(1:1) == 'm' ) then
            n_atom_uc = igr - 1
            n_atom_per = n_atom_uc
          elseif( keyword(1:1) == 'f' ) then
            n_atom_per = igr - 1
          elseif( keyword(1:1) == 's' ) then
            n_atom_sur = igr - 1
          elseif( keyword(1:1) == 'i' ) then
            n_atom_int = igr - 1
          endif

! Description de l'agregat :
        case('cif_file','film_cif_')

          Fcif = .true.
          Matper = .true.
          Film = keyword == 'film_cif_'

          Fcif_file = ' '
          read(itape4,'(A)') Fcif_file
          Fcif_file = Adjustl(Fcif_file)
          Length = len_trim(Fcif_file)
          if( Length > 4 ) then
            if( Fcif_file(Length-3:Length-3) /= '.' ) Fcif_file(Length+1:Length+4) = '.cif'
          endif
          open(8, file = Fcif_file, status='old', iostat=istat)
          if( istat /= 0 ) call write_open_error(Fcif_file,istat,1)

          do

            read(8,'(A)',iostat=eoff) mot

            if( eoff /= 0 ) exit

            Length = len_trim(mot)
            do i = 1,Length
             if( mot(i:i) == char(9) ) mot(i:i) = ' '
            end do
            mot = adjustl(mot)

            if( mot(1:30) == '_symmetry_space_group_name_H-M' .or. mot(1:25) == '_space_group_name_H-M_alt') then

              Length = len_trim(mot)

              Apostrophe = .true.
              do i = 31,Length
               if( mot(i:i) == '''' ) exit
               if( mot(i:i) /= ' ' ) then
                 Apostrophe = .false.
                 exit
               endif
              end do
              if( .not. Apostrophe ) i = i - 1
              k = 0
              j = i + 1
              do i = j,Length
                if( mot(i:i) == ' ' ) then
                  if( Apostrophe ) then
                    cycle
                  else
                    exit
                  endif
                endif
                if( mot(i:i) == '''' ) exit
                k = k + 1
                Space_group(k:k) = mot(i:i)
              end do

            elseif( mot(1:17) == '_cell_angle_gamma' ) then
              Angz = number_from_text(1,mot)

            elseif( mot(2:16) == 'atom_site_label' .or. mot(2:22) == 'atom_site_symbol' .or. &
                    mot(2:16) == 'atom_site_fract' ) then

              backspace(8)
              n_fract_x = -1; n_fract_y = -1; n_fract_z = -1
              n_label = -1; n_symbol = -1

              do n = 1,100000
                read(8,'(A)') mot
                mot = adjustl(mot)
                if( mot(1:1) /= '_' ) then
                  backspace(8)
                  exit
                elseif( mot(2:16) == 'atom_site_label') then
                  n_label = n - 1
                elseif( mot(2:22) == 'atom_site_type_symbol') then
                  n_symbol = n - 1
                elseif( mot(2:18) == 'atom_site_fract_x') then
                  n_fract_x = n - 1
                elseif( mot(2:18) == 'atom_site_fract_y') then
                  n_fract_y = n - 1
                elseif( mot(2:18) == 'atom_site_fract_z') then
                  n_fract_z = n - 1
                elseif( mot(2:20) == 'atom_site_occupancy') then
                  Taux = .true.
                elseif( mot(2:25) == 'atom_site_U_iso_or_equiv' .or. mot(2:25) == 'atom_site_B_iso_or_equiv' ) then
                  Temp_B_iso = .true.
                endif
              end do

              if( n_fract_x == -1 .or. n_fract_y == -1 .or. n_fract_z == -1 .or. ( n_label == -1 .and. n_symbol == -1 ) ) then
                call write_error
                do ipr = 6,9,3
                  write(ipr,110)
                  write(ipr,165) Fcif_file
                end do
                stop
              endif

              igr = 0
              do
                read(8,'(a4)',iostat=eoff) mot
                mot = adjustl(mot)
                mot4(1:4) = mot(1:4)
                if( eoff /= 0 .or. mot4(1:1) == '_' .or. mot4 == 'loop' .or. mot4 == ' ' .or. mot4(1:1) == '#' ) exit
                igr = igr + 1
              end do
              exit

            endif

          end do

          Close(8)

          if( keyword(1:1) == 'c' ) then
            n_atom_uc = igr
          elseif( keyword(1:1) == 'f' ) then
            n_atom_per = igr
          endif

          if( Space_group == ' ' ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,166) Fcif_file
            end do
            stop
          endif

          if( n_atom_uc == 0 .and. n_atom_per == 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,167) Fcif_file
            end do
            stop
          endif

         case('pdb_file','film_pdb_')

          Pdb = .true.
          Taux = .true.
          Matper = .true.
          Film = keyword == 'film_pdb_'

          Fichier_pdb = ' '
          read(itape4,'(A)') Fichier_pdb
          Fichier_pdb = Adjustl(Fichier_pdb)
          Length = len_trim(Fichier_pdb)
          if( Length > 4 ) then
            if( Fichier_pdb(Length-3:Length-3) /= '.' ) Fichier_pdb(Length+1:Length+4) = '.pdb'
          endif
          open(8, file = Fichier_pdb, status='old', iostat=istat)
          if( istat /= 0 ) call write_open_error(Fichier_pdb,istat,1)

          igr = 0

          do

            read(8,'(a6)',iostat=eoff) mot6

            if( eoff /= 0 ) exit

            select case(mot6)

              case('END   ')

                exit

              case('CRYST1' )

                backspace(8)
                read(8,'(47x,f7.2,a13)') Angz, Spgr

                Length = len_trim(Spgr)
                j = 0
                do i = 1,Length
                  if( Spgr(i:i) == ' ' ) cycle
                  j = j + 1
                  Space_group(j:j) = Spgr(i:i)
                end do

              case('ATOM  ','HETATM')

                igr = igr + 1

              case default

                cycle

              end select

          end do

          Close(8)

          if( keyword(1:1) == 'p' ) then
            n_atom_uc = igr
          else
            n_atom_per = igr
          endif

        case('flapw','flapw_s','flapw_r','flapw_s_p','flapw_psi','flapw_n','flapw_n_p','flapw_s_n')

          Flapw = .true.
          Matper = .true.

          if( keyword(6:7) == '_s' ) then
            n = nnombre(itape4,132)
            read(itape4,*)
          elseif( keyword(6:7) == '_r' ) then
            Wien_save = - 1
            n = nnombre(itape4,132)
            read(itape4,*)
          endif

          do i = 1,2
            if( i == 2 .and. Wien_save == - 1 ) exit
            n = nnombre(itape4,132)
            read(itape4,'(A)') Wien_file(i)
            if( Wien_file(i)(1:1) == ' ' ) Wien_file(i) = Adjustl( Wien_file(i) )
          end do

          call lectpot_dim(n_atom_uc,nklapw,nlmlapwm,nmatsym,Wien_file(1),Wien_file(2),ntype,Wien_save)

        case('memory_sa')
          Memory_save = .true.

        case('screening')
          Screening = .true.

      end select

    end do boucle_gen

!*** JDB
    if( .NOT. adimpin .AND. Use_FDMX ) then
      allocate( E_adimp(4) )
      E_adimp(1:4) = (/ 100.0, 250.0, 400.0, 500.0 /)
    endif
!*** JDB

    nspin = min( nspin, nspinp )

    if( Extract ) Xan_atom = .false.
    if( Flapw .or. Extract_ten ) Hubbard = .false.
    if( Flapw ) nspin = nspinp
    if( nspin == 2  ) then
      Magnetic = .true.
    else
      Magnetic = .false.
    endif
    if( npldafs == 0 ) Self_abs = .false.
    if( npldafs == 0 ) Full_self_abs = .false.
    if( Self_abs ) n_mu = 2 * npldafs
    if( Full_self_abs ) n_mu = 4 * npldafs
    if( .not. Matper ) Space_Group = ' '
    if( Pdb ) Temp_B_iso = .true.

    if( ( npldafs /= 0 .or. npldafs_2d /= 0 ) .and. .not. Dafs_bio ) then
      Rewind(itape4)
      do igrdat = 1,100000
        read(itape4,'(A)') mot
        keyword = identmot(mot,9)
        if( keyword(1:4) /= 'dafs' ) cycle
        npldafs_t = 0
        do ipl = 1,npldafs
          n = nnombre(itape4,132)
          select case(n)
            case(3)
              read(itape4,*)
              npldafs_t = npldafs_t + 1
              nn = nnombre(itape4,132)
              select case(nn)
                case(2,3,4,5)
                  read(itape4,*)
                case(6)
                  read(itape4,*)
                  read(itape4,*)
              end select
            case(5,6,7)
              read(itape4,*)
              npldafs_t = npldafs_t + 1
            case(8,9,10,11)
              if( keyword(1:7) == 'dafs_2d' ) then
                npldafs_t = npldafs_t + 1
                read(itape4,*) r, r, r, i, j
                if( i < 0 .or. j < 0 ) npldafs_t = npldafs_t + 1
              else
                read(itape4,*) r, r, p(:)
                if( abs( p(2) ) < eps10 ) then
                  npldafs_t = npldafs_t + 1
                else
                  k = nint( ( p(3) - p(1) ) / p(2) ) + 1
                  k = max( 1, k )
                  npldafs_t = npldafs_t + k
                endif
              endif
            case(12,13)
              read(itape4,*) r, r, p(:), i, j
              if( abs( p(2) ) < eps10 ) then
                k = 1
              else
                k = nint( ( p(3) - p(1) ) / p(2) ) + 1
                k = max( 1, k )
              endif
              npldafs_t = npldafs_t + k
              if( keyword(1:7) == 'dafs_2d' .and. ( i < 0 .or. j < 0 ) ) npldafs_t = npldafs_t + k
          end select
        end do
        npldafs = npldafs_t
        exit
      end do
      if( keyword == 'dafs_2d' ) npldafs =  2 * npldafs
    endif

    if( Pol_dafs_in ) npldafs_e = npldafs

    if( n_adimp > 1 .and. n_radius > 1 .and. .not. Extract ) then

      n_range = n_adimp + n_radius - 1

      do i = 1,n_adimp - 1
        do j = 1,n_radius - 1
          if( abs( E_adimp(i) - E_radius(j) ) < 0.00001 ) then
            n_range = n_range - 1
            exit
          endif
        end do
      end do
      deallocate( E_adimp, E_radius )

    elseif( n_adimp > 1 .and. .not. Extract  ) then

      n_range = n_adimp
      deallocate( E_adimp )

    elseif( n_radius > 1 .and. .not. Extract ) then

      n_range = n_radius
      deallocate( E_radius )

    else

      n_range = 1

    endif

    if( Bulk ) deallocate( Z_bulk ) ! allocation done during the main reading

    if( n_atom_per == 0 ) n_atom_per = n_atom_uc
    n_atom_neq = n_atom_per + n_atom_sur + n_atom_int
    n_atom_per_neq = n_atom_per

    if( Space_Group /= ' ' .or. ntype == 0 .or. Atom_conf ) then
      Very_fast = Readfast .or. Taux .or. Temp_B_iso .or. .not. ( Atom_nonsph .or. Atom_occ_hubb .or. Axe_loc)
      call Dim_reading(Angz,Atom_conf,Atom_conf_dop,Fcif,Doping,Fcif_file,Fichier_pdb,itape4,itype_dop,n_atom_bulk,n_atom_int, &
                       n_atom_per,n_atom_per_neq,n_atom_sur,n_atom_neq,n_fract_x,n_fract_y,n_fract_z,n_label,n_symbol,ntype, &
                       ntype_conf,Pdb,Space_Group,Very_fast)
    endif
    n_atom_uc = n_atom_per + n_atom_sur + n_atom_int

! Pour le cas des atomes charges ou il faut ajouter une orbitale:
    if( nlatm > 0 .or. Screening ) nlatm = nlatm + 1

! ntype is the total number of atom types including the bulk atoms when existing
! ngroup is the total number of atoms including the bulk atoms when existing
    ngroup = n_atom_uc + n_atom_bulk
    if( Doping ) ngroup = ngroup + 1
    if( Doping .and. Atom_conf .and. .not. Atom_conf_dop ) ntype = ntype + 1

    n_atom(1) = n_atom_bulk; n_atom(2) = n_atom_cap; n_atom(3) = n_atom_int; n_atom(4) = n_atom_neq; n_atom(5) = n_atom_per
    n_atom(6) = n_atom_per_neq; n_atom(7) = n_atom_sur; n_atom(8) = n_atom_uc

  endif   ! arrive mpirank0 /= 0

  if( mpinodes0 > 1 ) then
    call MPI_Bcast(Absauto,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Atom_nonsph,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Atom_occ_hubb,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Axe_loc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Bulk,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Cap_layer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Dafs_bio,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Doping,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Memory_save,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Extract,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Extract_ten,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Film,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Flapw,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Full_self_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Hubbard,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(itype_dop,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Magnetic,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_adimp,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_atom,8,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_atom_coop,8,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_dic,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_file_dafs_exp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_mat_polar,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_mu,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_multi_run_e,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_q_rixs,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_radius,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_range,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_theta_in,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_two_theta,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_Z_abs,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nb_atom_conf_m,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(neimagent,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nenerg_in_rixs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nenerg_s,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ngamme,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ngroup,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nhybm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nklapw,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nlatm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nlmlapwm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nmatsym,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(norbdil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(npldafs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(npldafs_2d,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(npldafs_e,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nple,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nq_nrixs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nrixs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nspin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nspino,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nspinp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ntype,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ntype_bulk,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(ntype_conf,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Pdb,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Quadrupole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Self_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Taux,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Temp_B_iso,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Xan_atom,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Wien_save,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)


  endif

  if( nq_nrixs >= 1 .and. lmax_nrixs > 1 ) Quadrupole = .true.

  if( nple == 0 ) then
    if( Quadrupole ) then
      nplm = 6
    else
      nplm = 3
    endif
  else
    nplm = nple
  endif
  nple = nple + n_mat_polar

  nplrm = nplm + n_dic + n_mu + 4 * n_mat_polar
  ncolm = nplm + n_dic + 2 * npldafs + 2 * n_mu + 6 * n_mat_polar + 1
  if( xan_atom ) ncolm = ncolm + 1
  if( self_abs ) ncolm = ncolm  + 2 * npldafs
  if( Full_self_abs ) ncolm = ncolm  + 4 * npldafs

  return
  110 format(//'  Error in the indata file :')
  115 format(//' After keyword "Radius", the number of value given for the radius and energies',/ &
               ' is even (',i2,'), it must be odd !',/ &
               ' If it is not the case, check the presence of extra characters or tabulations.',/ &
               ' They are forbidden !'//)
  117 format(//' After keyword "Adimp", the number of value given for interpoint distances and energies',/ &
               ' is even (',i2,'), it must be odd !',/ &
               ' If it is not the case, check the presence of extra characters or tabulations.',/ &
               ' They are forbidden !'//)
  120 format(//' After keyword "',a9,'", the number of value given for the energies and steps',/ &
               ' is even (',i2,'), it must be odd !',/ &
               ' If it is not the case, check the presence of extra characters or tabulations.',/ &
               ' They are forbidden !'//)
  130 format(//' Energy step is zero or negative after keyword "',a9,'", forbidden !'//)
  135 format(///' Just after the keyword "Atom_nsph", the number of digit must be >= 3 !',// &
                ' ( They are: nb of atoms, indexes of the atoms, nb of non-spherical orbitals' )
  140 format(//' After the keyword "Dafs", for the reflection number', i3,',',/ &
               ' the polarization or the indexes are not well set.',/ &
               ' Check the format !'//)
  150 format(///' Just after the keyword "',A,'"', ', the following line is read :',/)
  160 format(/' It must not be there or it contains unwanted characters !'//)
  165 format(///' Error in the lecture of the Cif_file :',A,/ &
                '   The position of the atoms are not found',/ &
                '   (with tag _atom_site_fract_x, _atom_site_fract_y and _atom_site_fract_z) !' //)
  166 format(///' Error in the lecture of the Cif_file :',A,/ &
                '   No space group found ( with tag _symmetry_space_group_name_H-M ) !' //)
  167 format(///' Error in the lecture of the Cif_file :',A,/ &
                '   No list of atom position found !' //)
end

!*********************************************************************

subroutine write_err_form(irec,keyword)

  implicit none

  integer:: ipr, irec, pb_line

  character(len=9):: keyword
  character(len=132):: mot

  call write_error
  backspace(irec)
  read(irec,'(A)',iostat=pb_line) mot

  do ipr = 6,9,3
    write(ipr,110) keyword
    if( pb_line > 0 ) then
      write(ipr,'(//5x,A//)') ' Check if the line is terminated by a cariage return !'
    else
      write(ipr,'(//A/)') ' The line is:'
      write(ipr,'(A)') mot
      write(ipr,120)
    endif
  end do
  stop

  return
  110 format(//' Format error when reading in the indata file under the keyword:',//,5x,A)
  120 format(//' Check :',/ '  - How many numbers must be in the line ?',/ &
                            '  - Are there spaces between the numbers ?',/ &
                            '  - Are there unwanted characters, extra  points... ?',//)
end

!***********************************************************************

! Lecture subroutine of potentials and electronic densities comming from Wien2k

subroutine lectpot_dim(n_atom_uc,nklapw,nlmlapwm,nmatsym,nomstruct,nomvcoul,ntype,Wien_save)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  character(len=1) Trans
  character(len=Length_name) nomstruct, nomvcoul

  integer:: Wien_save
  integer, dimension(:), allocatable :: nrato_lapw

  open(8, file = nomstruct, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(nomstruct,istat,1)

  read(8,*)
  read(8,'(a1,27x,i2)') Trans, ntype
  select case(Trans)
    case('F')
      ntrans = 3
    case('B','C')
      ntrans = 1
    case default
      ntrans = 0
  end select
  read(8,*)
  read(8,*)

  allocate( nrato_lapw(ntype) )

  index = 0
  do ia = 1,ntype
    index = index + ntrans + 1
    read(8,*)
    read(8,'(15x,i2)') mult
    do mu = 1,mult-1
      index = index + ntrans + 1
      read(8,*)
    end do
    read(8,'(15x,i5)') nrato_lapw(ia)
    read(8,*)
    read(8,*)
    read(8,*)
  end do

  read(8,'(i4)') nmatsym

  close(8)

  n_atom_uc = index

  if( Wien_save /= - 1 ) then

    open(8, file = nomvcoul, status='old', iostat=istat)
    if( istat /= 0 ) call write_open_error(nomvcoul,istat,1)

    do ia = 1,ntype

      do i = 1,4
        read(8,*)
      end do
      read(8,'(16x,i2)') ll
      nlmlapwm = max( abs(ll), nlmlapwm )

      do l = 1,ll
        do i = 1,4
          read(8,*)
        end do
        do j = 1,nrato_lapw(ia),4
          read(8,*)
        end do
      end do

      do i = 1,3
        read(8,*)
      end do

    end do ! end of loop over atoms

    do i = 1,5
      read(8,*)
    end do
    read(8,'(14x,i5)') nklapw

    close(8)

  endif

  deallocate ( nrato_lapw )

  return
end

!*********************************************************************

subroutine Dim_reading(Angz,Atom_conf,Atom_conf_dop,Fcif,Doping,Fcif_file,Fichier_pdb,itape4,itype_dop,n_atom_bulk,n_atom_int, &
                       n_atom_per,n_atom_per_neq,n_atom_sur,n_atom_neq,n_fract_x,n_fract_y,n_fract_z,n_label,n_symbol,ntype, &
                       ntype_conf,Pdb,Space_Group,Very_fast)

  use declarations
  implicit none

  integer:: eof, i, ier, igr, igrdat, io, ipr, istat, it, itape4, itype_dop, jgr, kgr, ligne, n, n_atom_bulk, n_atom_neq, &
    n_atom_int, n_atom_per, n_atom_per_neq, n_atom_sur, n_label, n_fract_x, n_fract_y, n_fract_z, na, n_symbol, nnombre, norbv, &
    ntype, ntype_conf

  integer, dimension(n_atom_per):: neq, numat
  integer, dimension(n_atom_neq+n_atom_bulk):: igra, itype

  character(len=2):: Chemical_Symbol, Chemical_Symbol_c, Symbol
  character(len=6):: mot6
  character(len=9):: keyword
  character(len=13):: Space_Group
  character(len=132):: identmot, mot, mot132, Word_from_text
  character(len=Length_name):: Fcif_file, Fichier_pdb

  logical:: Atom_conf, Atom_conf_dop, Fcif, Doping, Pdb, Very_fast

  real(kind=db):: Angz, number_from_text
  real(kind=db), dimension(3):: p
  real(kind=db), dimension(3,3):: Mat
  real(kind=db), dimension(3,n_atom_per_neq):: posn, posout

  if( Pdb ) then

    open(8, file = Fichier_pdb, status='old', iostat=istat)

    igr = 0

    do ligne = 1,100000

      read(8,'(a6)') mot6

      select case(mot6)

        case('END   ')

          exit

        case('SCALE1' )

          backspace(8)
          do i = 1,3
            read(8,'(10x,3f10.6)') Mat(i,:)
          end do

        case('ATOM  ','HETATM')

          backspace(8)
          read(8,'(30x,3f8.3,22x,a2)') p(:), Symbol

          igr = igr + 1
          p = Matmul( Mat, p )
          posn(:,igr) = p(:)

          Symbol = adjustl(Symbol)
          do i = 1,103
            if( Chemical_Symbol_c(i) /= Symbol .and. Chemical_Symbol(i) /= Symbol ) cycle
            numat(igr) = i
            exit
          end do
          if( i == 104 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,170) igr, Symbol
            end do
            stop
          endif
          if( igr == n_atom_per ) exit

        case default

          cycle

      end select

    end do

    Close(8)

  elseif( Fcif ) then

    open(8, file = Fcif_file, status='old', iostat=istat)

    igr = 0

    do
      read(8,'(A)') mot
      mot = adjustl(mot)
      if( mot(2:16) == 'atom_site_label' .or. mot(2:22) == 'atom_site_type_symbol' .or. mot(2:16) == 'atom_site_fract' ) exit
    end do
    backspace(8)

    n = 0
    do
      n = n + 1
      read(8,'(A)') mot
      mot = adjustl(mot)
      if( mot(1:1) /= '_' ) then
        backspace(8)
        exit
      endif
    end do

    do igr = 1,n_atom_per

      read(8,'(A)') mot
      if( n_symbol == 0 ) then
        mot132 = word_from_text(n_label,mot)
      else
        mot132 = word_from_text(n_symbol,mot)
      endif

! Lecture element chimique
      Symbol = mot132(1:2)
      if( Symbol(2:2) == '1' .or. Symbol(2:2) == '2' .or. Symbol(2:2) == '3' .or. Symbol(2:2) == '4' .or.  &
          Symbol(2:2) == '5' .or. Symbol(2:2) == '6' .or. Symbol(2:2) == '7' .or. Symbol(2:2) == '8' .or.  &
          Symbol(2:2) == '9' .or. Symbol(2:2) == '+' .or. Symbol(2:2) == '-' ) Symbol(2:2) = ' '

      do i = 1,103
        if( Chemical_Symbol_c(i) /= Symbol .and. Chemical_Symbol(i) /= Symbol ) cycle
        numat(igr) = i
        exit
      end do
      if( i == 104 ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,170) igr, Symbol
        end do
        stop
      endif

      posn(1,igr) = number_from_text(n_fract_x,mot)
      posn(2,igr) = number_from_text(n_fract_y,mot)
      posn(3,igr) = number_from_text(n_fract_z,mot)

    end do

    Close(8)

  elseif( Space_Group /= ' ' .or. ntype == 0 ) then

    Rewind(itape4)

    do igrdat = 1,100000
      read(itape4,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit
      keyword = identmot(mot,9)
      if( keyword(1:7) == 'crystal' .or. keyword(1:4) == 'film'  .or. keyword(1:7) == 'molecul' ) exit
    end do

    if( eof == 0 ) then ! can corresponds to the case, there is only "surface" or "interface"

      n = nnombre(itape4,132)
      read(itape4,*)

      do igr = 1,n_atom_per_neq
        if( Very_fast ) then
          read(itape4,*) numat(igr), posn(:,igr)
        else
          n = nnombre(itape4,132)
          if( n == 0 ) exit
          if( n == 2 .or. n == 3 ) then
            read(itape4,*)
            n = nnombre(itape4,132)
          endif
          norbv = 0
          select case(n)
            case(4)
              read(itape4,*) numat(igr), posn(:,igr)
            case default
              read(itape4,*,iostat=ier) numat(igr), posn(:,igr), norbv
              if( ier > 0 ) call write_err_form(itape4,keyword)
          end select
          if( norbv == 0 ) cycle
          do io = 1,abs(norbv)
            read(itape4,*)
          end do
        endif

      end do

    endif

  endif

  if( Space_Group /= ' ' ) then

! Atom defined with few digits
    if( abs( Angz - 120._db ) < 0.0001_db ) then
      do i = 1,11
        if( i == 3 .or. i == 6 .or. i == 9 ) cycle
        Where( abs( Posn - i / 12._db ) < 0.00005_db ) Posn = i / 12._db
      end do
    endif

    call spgroup(Fcif,Fcif_file,.false.,neq,n_atom_per,n_atom_per_neq,numat,posn,posout,Space_group)

  endif

  itype(1:n_atom_per_neq) = numat(1:n_atom_per_neq)

  if( ntype == 0 ) then

    ntype = ntype_conf

    na = 0
    igra(:) = 0
    if( Atom_conf ) then
      igra(:) = 0
      Rewind(itape4)

      do igrdat = 1,100000
        read(itape4,'(A)') mot
        keyword = identmot(mot,9)
        if( keyword /= 'atom_conf' ) cycle

        do it = 1,ntype_conf
          read(itape4,*) n, igra(na+1:na+n)
          na = na + n
        end do
        exit
      end do
    endif

    if( n_atom_int > 0 ) then

      Rewind(itape4)

      do igrdat = 1,100000
        read(itape4,'(A)') mot
        keyword = identmot(mot,9)
        if( keyword(1:7) == 'interfa' ) exit
      end do

      read(itape4,*)

      do igr = n_atom_per_neq + 1, n_atom_per_neq + n_atom_int
        read(itape4,*) itype(igr)
      end do

    endif

    if( n_atom_sur > 0 ) then

      Rewind(itape4)

      do igrdat = 1,100000
        read(itape4,'(A)') mot
        keyword = identmot(mot,9)
        if( keyword(1:7) == 'surface' ) exit
      end do

      read(itape4,*)

      do igr = n_atom_per_neq + n_atom_int+ 1, n_atom_neq
        read(itape4,*) itype(igr)
      end do

    endif

    if( n_atom_bulk > 0 ) then

      Rewind(itape4)

      do igrdat = 1,100000
        read(itape4,'(A)') mot
        keyword = identmot(mot,9)
        if( keyword(1:4) == 'bulk' ) exit
      end do

      read(itape4,*)

      do igr = n_atom_neq + 1, n_atom_neq + n_atom_bulk
        read(itape4,*) itype(igr)
      end do

    endif

    boucle_1: do igr = 1, n_atom_neq + n_atom_bulk
      do kgr = 1,na
        if( igra(kgr) == igr ) cycle boucle_1
      end do
      boucle_jgr: do jgr = 1,igr-1
        do kgr = 1,na
          if( igra(kgr) == jgr ) cycle boucle_jgr
        end do
        if( itype(igr) == itype(jgr) ) cycle boucle_1
      end do boucle_jgr
      ntype = ntype + 1
    end do boucle_1
    if( Doping .and. .not. Atom_conf_dop ) then
      do igr = 1,n_atom_neq
        if( itype(igr) == itype_dop ) exit
      end do
      if( igr > n_atom_neq ) ntype = ntype + 1
    endif

  endif

  return
  170 format(/' Error in the Pdb file :',// &
              ' The chemical symbol of atom number',i6,' is: ', a2,', what is not known!'//)
end

!*********************************************************************

! Reading subroutine
! All the indata are in Angstroem and eV.
! Angle are in degrees
! They are converted in this routine in atomic units (Bohr radius and Rydberg).

subroutine lecture(Absauto,adimp,alfpot,All_nrixs,All_site_rixs,Allsite,Ampl_rixs,Ang_borm,Ang_rotsup,Angle_mode,Angle_or, &
    Angpoldafs,Angxyz,Angxyz_bulk, &
    Angxyz_cap,Angxyz_int,Angxyz_sur,ATA,Atom_occ_hubb,Atom_nonsph,Atom_nsph_e,Atomic_scr,Axe_atom_gr,Axe_loc,axyz,axyz_bulk, &
    axyz_cap,axyz_int,axyz_sur,Basereel,Bormann,Bulk,Bulk_roughness,Cap_B_iso,Cap_layer,Cap_roughness, &
    Cap_shift,Cap_thickness,Cartesian_tensor,cdil,Center_s,Centre,Charge_free,Check_extract,Classic_irreg,Clementi,Com,Comt, &
    COOP,COOP_z_along_bond,Core_energ_tot,Core_resolved,Coupelapw,D_max_pot,Dafs,Dafs_bio,Delta_En_conv,Delta_Epsii, &
    Delta_helm,Density, &
    Density_comp,Dip_rel,Dipmag,Dist_coop,Doping,dpos,dyn_eg,dyn_g,E_adimp,E_radius,E_max_range,Eclie,Eclie_out,Ecrantage,Eeient, &
    Egamme,Eimagent,Eneg_i,Eneg_n_i,Energ_rixs,Energphot,Ephot_min,Extract,Extract_ten,f_no_res,FDM_comp,FDMX_only,Film, &
    Film_roughness,Film_shift,Film_thickness,Fit_cal,Flapw,Flapw_new,Force_ecr,Full_atom_e,Full_potential,Full_self_abs, &
    Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,Green_bulk,Green_int,Green_s,Green_self,Harm_cubic,Helm_cos, &
    hkl_borm,hkl_dafs,hkl_film,hkl_ref,Hubb, &
    Hubbard,Hybrid,iabsm,iabsorig,icheck,icom,igr_coop,igr_dop,indice_par,Interface_shift,iscratch,isigpi,itdil,its_lapw,iord, &
    itape4,itype,itype_dop,jseuil,Kern_fac,Kern_fast,Kgroup,korigimp,lmax_nrixs,lamstdens,ldil,lecrantage,Length_line,lin_gam, &
    lmax_pot,lmax_tddft_inp,lmaxfree,lmaxso_max,lmaxso0,lmaxat0,lmoins1,lplus1,Lseuil,lvval,m_hubb_e,Magnetic,Mat_or,Mat_UB, &
    Matper,Moment_conservation,mpinodes,mpinodes0,mpirank,mpirank0,Muffintin,Multipole, &
    multrmax,n_adimp,n_atom,n_atom_bulk,n_atom_cap,n_atom_coop,n_atom_uc,n_atom_proto,n_devide,n_file_dafs_exp,n_mat_polar, &
    n_multi_run_e,n_q_rixs,n_radius,n_range,n_theta_in,n_two_theta,n_Z_abs,nb_atom_conf_m,nbseuil,nchemin,necrantage, &
    neimagent,nenerg_in_rixs,nenerg_s,ngamh,ngamme,ngroup,ngroup_hubb,ngroup_lapw,ngroup_m,ngroup_nonsph,ngroup_par,ngroup_pdb, &
    ngroup_taux,ngroup_temp,nhybm,nlat,nlatm,No_DFT,No_solsing,nom_fich_Extract, &
    nomfich,nomfichbav,Noncentre, &
    Nonexc,norbdil,norbv,Normaltau,normrmt,npar,nparm,nphi_dafs, &
    nphim,npl_2d,npldafs,npldafs_2d,npldafs_e,npldafs_f,nple,nposextract,nq_nrixs,nrato,nrato_dirac,nrato_lapw,nrm, &
    nself,nself_bulk,nseuil,nslapwm,nspin,nsymextract,ntype,ntype_bulk,ntype_conf,ntype_mod,numat,numat_abs, &
    nvval,occ_hubb_e,Occupancy_first,Octupole,Old_zero,One_run,One_SCF,Operation_mode,Operation_mode_used,Optic,Overad,Overlap, &
    p_self_max,p_self0,p_self0_bulk,param,Pas_SCF,pdpolar,phi_0,PointGroup,PointGroup_Auto,polar,Polarise,poldafsem,poldafssm, &
    pop_nonsph,popats,popval,posn,posn_bulk,posn_cap,Powder,q_nrixs,q_rixs,Quadmag,Quadrupole,R_rydb,R_self,R_self_bulk, &
    r0_lapw,Ray_max_dirac,rchimp,Readfast,Relativiste,Renorm,RIXS,RIXS_core,RIXS_fast,Rlapw,Rmt,Rmtimp,Rot_Atom_gr,rotloc_lapw, &
    roverad,RPALF,rpotmax,rydberg,rsorte_s,SCF_log,Self_abs, &
    Solsing_s,Solsing_only,Spherical_signal,Spherical_tensor,Spin_channel,Spinorbite,State_all,State_all_out, &
    Step_loss,Sum_rixs,Supermuf,Surface_shift,Symauto,Symmol,Taux,Taux_cap,Taux_oc,Tddft,Temp,Temp_coef,Temp_B_iso, &
    Tensor_imp,Test_dist_min,Theta_in,Trace_format_wien,Trace_k,Trace_p,Two_theta,Typepar,Use_FDMX,V_helm,V_hubbard,V_intmax, &
    Vec_orig,Vecdafsem,Vecdafssm,Veconde,V0bdcFimp,Width_helm,Wien_file,Wien_matsym,Wien_save,Wien_taulap,Write_modQ, &
    Ylm_comp_inp,Z_bulk,Z_cap,Z_nospinorbite)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: eof, eoff, i, i_range, ier, ia, ie, igr, igr_dop, igrdat, io, iord, ip, ipar, ipl, ipl0, ipr, ipr0, &
    iscratch, isp, ispin, istat, istop, isymeq, it, itape4, itype_dop, j, jgr, jpl, jseuil, jt, k, kgr, kt, l, l_hubbard, &
    l_level_val, l1, l2, lamstdens, lecrantage, Length, Length_line, lin_gam, lmax_nrixs, lmax_pot, lmax_pot_default, &
    lmax_tddft_inp, lmaxat0, lmaxso0, lmaxso_max, long, Lseuil, m, m_hubb_e, MPI_host_num_for_mumps, mpierr, mpinodes, &
    mpinodes0, mpirank, mpirank0, &
    mpirank_in_mumps_group, multi_run, multrmax, n, n_adp_type, n_adimp, n_atom_bulk, n_atom_cap, n_atom_coop, n_atom_int, &
    n_atom_neq,  n_atom_neq_min, n_atom_per, n_atom_per_neq, n_atom_sur, n_atom_proto, n_atom_uc, n_B_iso, n_rchimp_z, n_devide, &
    n_file_dafs_exp, n_fract_x, n_fract_y, n_fract_z, &
    n_hubbard_Z, n_label, n_mat_polar, n_multi_run_e, n_occupancy, n_q_rixs, n_radius, n_range, n_rmtg_Z, n_symbol, n_theta_in, &
    n_two_theta, n_Z_abs, n1, n2, natomsym, nb_atom_conf_m, nb_atom_nsph, nbseuil, nchemin, &
    necrantage, neimagent, nenerg_in_rixs, nenerg_s, ngamme, ngamh, ngroup, ngroup_hubb, ngroup_lapw, ngroup_m, ngroup_nonsph, &
    ngroup_par, ngroup_pdb, ngroup_taux, ngroup_temp, nhybm, nlatm, nn, nnombre, non_relat, norb, norbdil, normrmt, &
    nparm, nphim, nphimt, npldafs, npldafs_2d, npldafs_e, npldafs_f, nple, nq_nrixs, nrato_dirac, nrm, nscan, nself, nself_bulk, &
    nseuil, nslapwm, nspin, ntype, ntype_bulk, ntype_conf, ntype_mod, Trace_k, Wien_save, Z_nospinorbite, Z

  character(len=1):: Let
  character(len=2):: Chemical_Symbol, Chemical_Symbol_c, Symbol
  character(len=3):: Seuil
  character(len=5):: Solver, Struct, Struct_bulk
  character(len=6):: mot6
  character(len=8):: dat, PointGroup
  character(len=9):: keyword
  character(len=10):: tim
  character(len=11):: motpol
  character(len=13):: Chemical_Name, mot13, Space_Group, Spgr
  character(len=50):: com_date, com_time
  character(len=132):: comt, Error_message, identmot, mot, mot132, Word_from_text
  character(len=Length_name):: Fichier, Fcif_file, nomfich, nom_fich_Extract, nomfichbav
  character(len=35), dimension(0:ntype):: com
  character(len=Length_name), dimension(9):: Wien_file
  character(len=Length_name), dimension(n_file_dafs_exp):: File_dafs_exp
  character(len=9), dimension(ngroup_par,nparm):: typepar

  integer, dimension(2):: Mult_bulk, Mult_film
  integer, dimension(3):: hkl_borm, hkl_ref, ldipimp
  integer, dimension(8):: n_atom
  integer, dimension(12):: Tensor_imp
  integer, dimension(30):: icheck
  integer, dimension(3,3):: lquaimp
  integer, dimension(n_atom_coop):: igr_coop
  integer, dimension(n_multi_run_e):: iabsm, iabsorig, nposextract, nsymextract
  integer, dimension(ngroup_par):: npar
  integer, dimension(n_Z_abs):: numat_abs
  integer, dimension(ngroup_lapw):: its_lapw
  integer, dimension(norbdil):: itdil, ldil
  integer, dimension(0:ntype):: icom, nlat, nrato, nrato_lapw, numat
  integer, dimension(ngroup):: itype
  integer, dimension(ngroup_pdb):: Kgroup
  integer, dimension(n_atom_bulk):: Z_bulk
  integer, dimension(n_atom_cap):: Z_cap
  integer, dimension(Z_Mendeleiev_max):: Z_rchimp, Z_hubbard, Z_rmtg

  integer, dimension(0:ngroup_nonsph):: norbv
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(npldafs_2d):: npl_2d
  integer, dimension(2,npldafs):: isigpi
  integer, dimension(5,npldafs_2d):: Operation_mode
  integer, dimension(ntype_conf):: nb_atom_conf
  integer, dimension(nb_atom_conf_m,ntype_conf):: igra
  integer, dimension(ngroup_par,nparm):: indice_par
  integer, dimension(0:ntype,nlatm):: lvval, nvval
  integer, dimension(3,3,nslapwm):: Wien_matsym

  integer, dimension(:), allocatable :: igr_nsph, itZ, neq

  complex(kind=db), dimension(3,npldafs_e):: Poldafsem, Poldafssm
  complex(kind=db), dimension(nhybm,16,ngroup_nonsph):: Hybrid

  logical:: Absauto, Absorber, adimpin, All_nrixs, All_site_rixs, Allsite, Ampl_rixs, Apostrophe, ATA, Atom, Atom_B_iso, &
    Atom_conf, Atom_nonsph, &
    Atom_occ_hubb, Atomic_scr, Atom_U_iso, Avoid, Axe_loc, Basereel, Bormann, Bulk, Cartesian_tensor, Charge_free, Centre_auto, &
    Centre_auto_abs, Center_s, Check_extract, Classic_irreg, Clementi, COOP, Core_energ_tot, Core_resolved, COOP_z_along_bond, &
    Core_resolved_e, Coupelapw, Cap_layer, Cylindre, Dafs, Dafs_bio, Density, Density_comp, Diagonal, &
    Dip_rel, Dipmag, Doping, dyn_eg, dyn_g, E1E1, E1E2e, E1E3, E1M1, E1M2, E2E2, E3E3, eneg_i, eneg_n_i, Energphot, Fcif, Film, &
    exc_imp, Extract, Extract_ten, FDM_comp, FDMX_only, Fermi_auto, Fit_cal, Flapw, Flapw_new, Force_ecr, Full_atom_e, &
    Full_potential, Full_self_abs, Gamma_hole_imp, Gamma_tddft, Green_bulk, Green_s, Green_self, Green_int, Harm_cubic, Hedin, &
    Helm_cos, hkl_film, Hubbard, Interface_shift_given, Kern_fac_default, Kern_fast, korigimp, lmaxfree, lmoins1, lplus1, M1M1, &
    M1M2, M2M2, Magnetic, Matper, Moment_conservation, Monocrystal, muffintin, No_DFT, no_core_resolved, no_dipquad, &
    no_e1e3, no_e2e2, no_e3e3, No_renorm, No_SCF_bulk, No_solsing, noncentre, Nonexc, nonexc_imp, &
    normaltau, Occupancy_first, Octupole, Old_zero, One_run, One_SCF, Operation_mode_used, Optic, Overad, &
    Pas_SCF_imp, Pdb, Perdew, PointGroup_Auto, Polarise, Powder, quadmag, Quadrupole, r_self_imp, Readfast, Relativiste, &
    Renorm, RIXS, RIXS_core, RIXS_fast, rpalf, rydberg, SCF, SCF_bulk, SCF_elecabs, SCF_mag_free, Self_abs, Self_cons, &
    Self_cons_bulk, SCF_exc_imp, Self_nonexc, Self_nonexc_imp, &
    solsing_only, solsing_s, spherical_signal, spherical_tensor, spherique, Spin_channel, Spinorbite, &
    State_all, State_all_out, Sum_rixs, Supermuf, Surface_shift_given, Symauto, Symmol, Taux, Tddft, Temp_B_iso, &
    Trace_format_wien, Use_FDMX, Write_modQ, Ylm_comp_inp

  logical, dimension(7):: SCF_log
  logical, dimension(10):: Multipole
  logical, dimension(ngroup):: Atom_nsph_e
  logical, dimension(0:ntype):: Hubb
  logical, dimension(ntype):: Suppress

  real(kind=db):: Alfpot, Ang_borm, Bulk_roughness, Cap_B_iso, Cap_roughness, Cap_shift, Cap_thickness, Cap_U_iso, &
    D_max_pot, Delta_En_conv, Delta_Epsii, Delta_helm, Eclie, Eclie_out, Ephot_min, Film_roughness, Film_thickness, &
    Film_zero, g1, g2, Gamma_max, Kern_fac, number_from_text, Overlap, p_self_max, p_self0, p_self0_bulk, Pas_SCF, phi, phi_0, &
    pop_nsph, pp, q, r, r_self, R_self_bulk, R_rydb, Ray_max_dirac, rn, Roverad, Rpotmax, Step, Step_azim, Step_loss, t, tc, &
    Temp, Test_dist_min, Theta, Theta_min, V_helm, V_intmax, vv, Width_helm, z_min

  real(kind=db), dimension(4,nq_nrixs):: q_nrixs
  real(kind=db), dimension(2):: Dist_coop, f_no_res
  real(kind=db), dimension(3):: Ang, Ang_rotsup, Ang_spin, Angxyz, Angxyz_bulk, Angxyz_cap, Angxyz_int, Angxyz_sur, Axe, Axe_spin, &
                                axyz, axyz_bulk, axyz_cap, axyz_int, axyz_sur, Centre, dpos, p, Vec_orig
  real(kind=db), dimension(4):: Film_shift, Interface_shift, Surface_shift
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(3,3):: Cubmat, Cubmat_bulk, Cubmati, Cubmati_bulk, Mat, Mat_or, Mat_UB, Rot, Rot_gen
  real(kind=db), dimension(norbdil):: cdil
  real(kind=db), dimension(0:ntype):: r0_lapw, rchimp, rlapw, rmt, rmtimp, V_hubbard
  real(kind=db), dimension(neimagent):: eeient, eimagent
  real(kind=db), dimension(ngamme):: Egamme
  real(kind=db), dimension(nspin):: ecrantage, V0bdcFimp
  real(kind=db), dimension(n_adimp):: Adimp
  real(kind=db), dimension(n_radius):: Rsorte_s
  real(kind=db), dimension(n_adimp):: E_adimp
  real(kind=db), dimension(n_radius):: E_radius
  real(kind=db), dimension(n_range):: E_max_range
  real(kind=db), dimension(nenerg_in_rixs):: Energ_rixs
  real(kind=db), dimension(n_theta_in):: Theta_in
  real(kind=db), dimension(n_two_theta):: Two_theta
  real(kind=db), dimension(Z_Mendeleiev_max):: Rchimp_Z, Hubbard_Z, Rmtg_Z
  real(kind=db), dimension(3,nple):: polar, veconde
  real(kind=db), dimension(3,npldafs):: hkl_dafs
  real(kind=db), dimension(4,n_q_rixs):: q_rixs
  real(kind=db), dimension(nple,2):: pdpolar
  real(kind=db), dimension(3,0:ntype):: Ang_base_loc
  real(kind=db), dimension(n_file_dafs_exp):: Angle_dafs_exp
  real(kind=db), dimension(n_atom_cap):: Taux_cap
  real(kind=db), dimension(ngroup_taux):: Taux_oc
  real(kind=db), dimension(ngroup_temp):: Temp_coef
  real(kind=db), dimension(3,n_atom_uc):: posn
  real(kind=db), dimension(3,ngroup):: Ang_base_loc_gr
  real(kind=db), dimension(3,n_atom_bulk):: posn_bulk
  real(kind=db), dimension(3,n_atom_cap):: posn_cap
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_gr
  real(kind=db), dimension(-m_hubb_e:m_hubb_e,-m_hubb_e:m_hubb_e,nspin,ngroup_hubb):: occ_hubb_e
  real(kind=db), dimension(3,3,ngroup_m):: Rot_atom_gr
  real(kind=db), dimension(3,3,ngroup_lapw):: rotloc_lapw
  real(kind=db), dimension(npldafs_f):: Angle_or
  real(kind=db), dimension(3,npldafs):: angpoldafs
  real(kind=db), dimension(3,npldafs_e):: Poldafse, Poldafsei, Poldafss, Poldafssi, Vecdafsem, Vecdafssm
  real(kind=db), dimension(3,npldafs_2d):: Angle_mode
  real(kind=db), dimension(ngroup_par,nparm):: param
  real(kind=db), dimension(nhybm,ngroup_nonsph):: pop_nonsph
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats
  real(kind=db), dimension(0:ntype,nlatm,nspin):: popval
  real(kind=db), dimension(3,nslapwm):: Wien_taulap

  real(kind=db), dimension(:), allocatable:: Eg, Hyb, x
  real(kind=db), dimension(:,:), allocatable:: pos
  real(kind=db), dimension(:,:,:), allocatable:: Hybrid_i, Hybrid_r

  n_atom_int = n_atom(3); n_atom_neq = n_atom(4); n_atom_per = n_atom(5)
  n_atom_per_neq = n_atom(6); n_atom_sur = n_atom(7)

  n_atom_neq_min = min( n_atom_neq, n_atom_uc ) ! because of the fake cif_file from Vista...

! Default values of parametres
  Absorber = .false.
  adimp(:) = 0.25_db
  adimpin = .false. !*** JDB
  alfpot = 0._db
  All_nrixs = .false.
  All_site_rixs = .false.
  Allsite = .false.
  Ampl_rixs = .false.
  Ang_base_loc(1,:) = -10000._db; Ang_base_loc(2:3,:) = 0._db
  Ang_base_loc_gr(1,:) = -10000._db; Ang_base_loc_gr(2:3,:) = 0._db
  Ang_rotsup(:) = 0._db
  Angle_mode(:,:) = 0._db
  Angpoldafs(:,:) = 10000._db
  Angxyz_bulk(:) = 0._db
  angxyz_cap(:) = 0._db
  angxyz_int(:) = 0._db
  angxyz_sur(:) = 0._db
  ATA = .false.
  Atom = .false.
  Atom_B_iso = .false.
  Atom_conf = .false.
  Atom_U_iso = .false.
  Atomic_scr = .false.
  Ang_spin(:) = 0._db
  Axe_spin(1) = 0._db; Axe_spin(2) = 0._db; Axe_spin(3) = 1._db
  axyz(:) = 0._db
  axyz_int(:) = 0._db
  axyz_sur(:) = 0._db
  axyz_bulk(1:3) = 0._db
  axyz_cap(1:3) = 0._db
  Basereel = .true.
  Bulk_roughness = 0._db
  Cap_B_iso = 0._db
  Cap_roughness = 0._db
  Cap_shift = -1000._db
  Cap_thickness = -1000._db
  Cap_U_iso = 0._db
  Cartesian_tensor = .false.
  cdil(:) = 0._db
  Centre(:) = 0._db
  Centre_auto = .false.
  Centre_auto_abs = .false.
  Center_s = .false.
 ! The irregular solution is taken by continuity with bessel if false.
  Classic_irreg = .true.
  Clementi = .false.
  Fcif = .false.
  com(:) = ' Dirac'
  COOP = .false.
  COOP_z_along_bond = .true.
  Core_energ_tot = .false.
  Core_resolved_e = .false.
  Coupelapw = .false.
  Cylindre = .false.
  Delta_Helm = 0._db
  D_max_pot = 2.5_db
  Dafs = .false.
  Delta_En_conv = 0.1_db
  Delta_Epsii = 1000000._db * Rydb
  Density = .false.
  Density_comp = .false.
  Dip_rel = .false.
  dipmag = .false.
  Dist_coop(:) = 0._db
  dyn_eg = .false.
  dyn_g = .false.
  dpos(:) = 0._db
  E_adimp(:) = 100000000._db
  E_max_range(:) = 100000000._db
  E_radius(:) = 100000000._db
  E1E1 = .true.;  E1E2e = .false.; E1E3 = .false.
  E1M1 = .false.
  E1M2 = .true.
  E2E2 = .false.
  E3E3 = .false.
  Ecrantage(:) = 1._db / nspin
  Eclie = 0.2_db
  Eclie_out = 1._db
  Eneg_i = .false.
  Eneg_n_i = .false.
  Energphot = .false.
  Ephot_min = - 1000000._db
  Exc_imp = .false.
  f_no_res(1) = 1._db     ! mag
  f_no_res(2) = -100._db  ! mom
  Fdm_comp = .false.
  Fermi_auto = .true.
  Film_roughness = 0._db
  Film_shift(:) = 0._db
  Film_thickness = 0._db
  Film_zero = -1000._db
  Flapw_new = .false.
  Force_ecr = .false.
  Full_atom_e = .false.
  Full_potential = .false.
  Gamma_tddft = .false.
  Green_bulk = .false.
  Green_int = .false.
  Green_s = .false.
  Green_self = .true.
  Harm_cubic = .true.
  Hedin = .false.
  Helm_cos = .false.
  hkl_film = .false.
  hkl_ref(1:2) = 0; hkl_ref(3) = 1
  Hybrid(:,:,:) = ( 0._db, 0._db )
  iabsm(1) = 1
  ier = 0
  Interface_shift(:) = 0._db
  Interface_shift_given = .false.
  Charge_free = .false.
  icom(:) = 1
  igr_dop = 0
  igr_coop(:) = 0
  iord = 4
  isigpi(:,:) = 0
  istop = 0
  Hubb(:) = .false.
  Kern_fac_default = .true.
  Kern_fast = .false.
  korigimp = .false.
  lamstdens = -1
  ldipimp(:) = -1
  lecrantage = 0
  lquaimp(:,:) = -1
  lin_gam = -1
  lmax_nrixs = 2
  lmax_tddft_inp = -1
  lmaxso0 = -5
  lmaxso_max = 1000000
  lmaxat0 = -1
  lmaxfree = .false.
  lmax_pot = 0
  lmax_pot_default = 1
  lmoins1 = .false.
  lplus1 = .false.
  M1M1 = .false.
  M1M2 = .false.
  M2M2 = .false.
  Mat_UB(:,:) = 0._db
  Matper = .false.
  Moment_conservation = .false.
  Monocrystal = .false.
  muffintin = .false.
  multrmax = 1
  n_rchimp_Z = 0
  n_devide = 2
  n_hubbard_Z = 0
  n_rmtg_Z = 0
  nchemin = - 1
  necrantage = 0
  nlat(:) = 0
  no_core_resolved = .false.
  No_DFT = .false.
  no_dipquad = .false.
  no_e1e3 = .false.
  no_e3e3 = .false.
  no_e2e2 = .false.
  No_renorm = .false.
  No_SCF_bulk = .false.
  No_solsing = .false.
  Nonexc = .false.
  nonexc_imp = .false.
  noncentre = .false.
  non_relat = 0
  norbv(:) = 0
  normaltau = .false.
  normrmt = 1
  nphim = 180
  nrato(:) = 0
  nrato_dirac = 600
  nrato_lapw(:) = 0
  nrm = 0
  nself = 0
  nself_bulk = 0
  nsymextract(:) = 1
  do i = 1,n_multi_run_e
    nposextract(i) = i
  end do
  numat_abs(:) = 0
  if( Atom_occ_hubb ) occ_hubb_e(:,:,:,:) = 0._db
  Octupole = .false.
  Old_zero = .false.
  One_run = .false.
  One_SCF = .false.
  Operation_mode(:,:) = 0
  Operation_mode_used = .false.
  Optic = .false.
  Overad = .false.
  Overlap = 0.1_db
  p_self_max = 1._db
  p_self0 = 0.1_db
  p_self0_bulk = 0.1_db
  Pas_SCF_imp = .false.
  Pdb = .false.
  pdpolar(:,:) = 0._db
  Width_helm = 0._db
  perdew = .false.
  phi_0 = 0._db
  PointGroup = ' '
  PointGroup_Auto = .true.
  polar(:,:) = 0._db
  polarise = .false.
  popats(:,:,:) = 0._db
  Powder = .false.
  q_nrixs(:,:) = 0._db
  Quadmag = .false.
  Quadrupole = .false.
  R_self_imp = .false.
  r0_lapw(:) = 0._db
  Rchimp(:) = 0._db
  Ray_max_dirac = 0._db
  Relativiste = .false.
  Renorm = .false.
  RIXS = .false.
  RIXS_core = .false.
  RIXS_fast = .false.
  rlapw(:) = 0._db
  Rmt(:) = 0._db
  Rmtimp(:) = 0._db
  rotloc_lapw(:,:,:) = 0._db
  Roverad = 0._db
  RPALF = .false.
  Rpotmax = 0._db
  Rsorte_s(:) = 3._db
  Rydberg = .false.
  R_rydb = 1._db
  SCF = .false.
  SCF_bulk = .false.
  SCF_elecabs = .false.
  SCF_mag_free = .false.
  Self_cons = .false.
  Self_cons_bulk = .false.
  SCF_exc_imp = .false.
  Self_nonexc = .true.
  Self_nonexc_imp = .false.
  Seuil = 'K1'
  solsing_s = .false.
  solsing_only = .false.
  Space_Group = ' '
  Spherical_tensor = .false.
  Spherical_signal = .false.
  Spherique = .false.
  Spinorbite = .false.
  State_all = .false.
  State_all_out = .false.
  Step_loss = -1._db
  Spin_channel = .false.
  Sum_rixs = .false.
  Supermuf = .false.
  Surface_shift(:) = 0._db
  Surface_shift_given = .false.
  Symauto = .true.
  Symmol = .false.
  Taux_cap(:) = 1._db
  Taux_oc(:) = 1._db
  Temp_coef(:) = 0._db
  Tddft = .false.
  Temp = 0._db
  Test_dist_min = 0.7_db * bohr ! minimum distance between 2 atoms
  Trace_format_wien = .false.
  Trace_k = 0
  Trace_p(:) = 0._db
  V_intmax = 1000000 * rydb
  V_helm = 0._db
  V_hubbard(:) = 0._db
  V0bdcFimp(:) = 0._db
  Vec_orig(1:2) = 0._db; Vec_orig(3) = 1._db
  Veconde(:,:) = 0._db
  Ylm_comp_inp = .false.
  Write_modQ = .false.
  Z_nospinorbite = 0

!Wien_file(1): struct, (2): vcoul, (3): r2v, (4): r2vdn, (5,6,7): clm(isp), (8) psii, (9) Save_potlapw
  Wien_file(:) = ' '
  Wien_file(8) = 'dirac'
  Wien_save = 0

  Multipole(1) = E1E1; Multipole(2) = E1E2e; Multipole(3) = E1E3
  Multipole(4) = E1M1; Multipole(5) = E1M2; Multipole(6) = E2E2
  Multipole(7) = E3E3; Multipole(8) = M1M1; Multipole(9) = M1M2
  Multipole(10) = M2M2

  if( mpirank0 == 0 ) then

    Rewind(itape4)

    write(6,*)

    boucle_lect: do igrdat = 1,100000

      read(itape4,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit boucle_lect

      keyword = identmot(mot,9)
      if( keyword(1:1) /= ' ' ) write(6,'(3x,A)') keyword

      select case(keyword)

        case('all_site_')
          All_site_rixs = .true.

        case('rixs_ampl')
          Ampl_rixs = .true.

        case('ata')
          ATA = .true.

        case('atom_b_is')
          Atom_B_iso = .true.

        case('atom_u_is')
          Atom_U_iso = .true.

        case('bulk')
          n = nnombre(itape4,132)
          if( n == 3 ) then
            read(itape4,*) axyz_bulk(1:3)
            angxyz_bulk(1:3) = 90._db
          else
            read(itape4,*,iostat=ier) axyz_bulk(1:3), angxyz_bulk(1:3)
          endif
          if( ier > 0 ) call write_err_form(itape4,keyword)
          do igr = 1,n_atom_bulk
            n = nnombre(itape4,132)
            if( n == 4 .or. ( .not. Temp_B_iso .and. .not. Taux ) ) then
              read(itape4,*) Z_bulk(igr), posn_bulk(:,igr)
            elseif( Temp_B_iso .and. Taux ) then
              if( Occupancy_first ) then
                if( n > 5 ) then
                  read(itape4,*) Z_bulk(igr), posn_bulk(:,igr), Taux_oc(ngroup-n_atom_bulk+igr), Temp_coef(ngroup-n_atom_bulk+igr)
                else
                  read(itape4,*) Z_bulk(igr), posn_bulk(:,igr), Taux_oc(ngroup-n_atom_bulk+igr)
                endif
              else
                if( n > 5 ) then
                  read(itape4,*) Z_bulk(igr), posn_bulk(:,igr), Temp_coef(ngroup-n_atom_bulk+igr), Taux_oc(ngroup-n_atom_bulk+igr)
                else
                  read(itape4,*) Z_bulk(igr), posn_bulk(:,igr), Temp_coef(ngroup-n_atom_bulk+igr)
                endif
              endif
            elseif( Temp_B_iso ) then
              read(itape4,*) Z_bulk(igr), posn_bulk(:,igr), Temp_coef(ngroup-n_atom_bulk+igr)
            else
              read(itape4,*) Z_bulk(igr), posn_bulk(:,igr), Taux_oc(ngroup-n_atom_bulk+igr)
            endif
          end do

       case('bulk_roug')
         read(itape4,*,iostat=ier) bulk_roughness
         if( ier > 0 ) call write_err_form(itape4,keyword)

        case('cap_layer')
          n = nnombre(itape4,132)
          if( n == 3 ) then
            read(itape4,*) axyz_cap(1:3)
            angxyz_cap(1:3) = 90._db
          else
            read(itape4,*,iostat=ier) axyz_cap(1:3), angxyz_cap(1:3)
          endif
          if( ier > 0 ) call write_err_form(itape4,keyword)
          do igr = 1,n_atom_cap
            n = nnombre(itape4,132)
            if( n == 4 ) then
              read(itape4,*) Z_cap(igr), posn_cap(:,igr)
            else
              read(itape4,*) Z_cap(igr), posn_cap(:,igr), Taux_cap(igr)
            endif
          end do

       case('cap_b_iso')
         read(itape4,*,iostat=ier) Cap_B_iso
         if( ier > 0 ) call write_err_form(itape4,keyword)

       case('cap_u_iso')
         read(itape4,*,iostat=ier) Cap_U_iso
         if( ier > 0 ) call write_err_form(itape4,keyword)
         Cap_B_iso = 8 * pi**2 * Cap_U_iso

       case('cap_rough')
         read(itape4,*,iostat=ier) Cap_roughness
         if( ier > 0 ) call write_err_form(itape4,keyword)

       case('cap_thick')
         read(itape4,*,iostat=ier) Cap_thickness
         if( ier > 0 ) call write_err_form(itape4,keyword)

        case('cap_shift')
          read(itape4,*,iostat=ier) Cap_shift
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('classic_i')
          Classic_irreg = .false.

        case('coop')
          n = nnombre(itape4,132)
          COOP = .true.
          if( n > 0 ) read(itape4,*,iostat=ier) igr_coop(1:n)

        case('coop_dist')
          n = nnombre(itape4,132)
          if( n == 1 ) then
            read(itape4,*,iostat=ier) Dist_coop(2)
          elseif( n > 1 ) then
            read(itape4,*,iostat=ier) Dist_coop(:)
            if( Dist_coop(2) < Dist_coop(1) ) then
              r = Dist_coop(1)
              Dist_coop(1) = Dist_coop(2)
              Dist_coop(2) = r
            endif
          endif

        case('coop_z_ax')
          COOP_z_along_bond = .false.

        case('core_ener')
          Core_energ_tot = .true.

        case('doping')
          read(itape4,*,iostat=ier) itype_dop, igr_dop
          if( ier > 0 ) call write_err_form(itape4,keyword)

       case('film_roug')
         read(itape4,*,iostat=ier) Film_roughness
         if( ier > 0 ) call write_err_form(itape4,keyword)

        case('film_shif')
          n = nnombre(itape4,132)
          if( n == 1 ) then
            read(itape4,*,iostat=ier) Film_shift(3)
          else
            n = min(4,n)
            read(itape4,*,iostat=ier) Film_shift(1:n)
          endif
          if(ier > 0 ) call write_err_form(itape4,keyword)

        case('film_zero')
          read(itape4,*,iostat=ier) Film_zero
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('film_thic')
          read(itape4,*,iostat=ier) Film_thickness
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('harm_tess')
          Harm_cubic = .false.

        case('hkl_film')
          hkl_film = .true.

        case('inter_shi')
          n = nnombre(itape4,132)
          Interface_shift_given = .true.
          if( n == 1 ) then
            read(itape4,*,iostat=ier) Interface_shift(3)
          elseif( n == 3 ) then
            read(itape4,*,iostat=ier) Interface_shift(1:3)
          else
            read(itape4,*,iostat=ier) Interface_shift(1:4)
          endif
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('surface_s')
          n = nnombre(itape4,132)
          Surface_shift_given = .true.
          if( n == 1 ) then
            read(itape4,*,iostat=ier) Surface_shift(3)
          elseif( n == 3 ) then
            read(itape4,*,iostat=ier) Surface_shift(1:3)
          else
            read(itape4,*,iostat=ier) Surface_shift(1:4)
          endif
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('dip_rel')
          Dip_rel = .true.

        case('fdm_comp')
          FDM_comp = .true.

        case('full_pote')
          Full_potential = .true.
          n = nnombre(itape4,132)
          if( n > 0 ) then
            read(itape4,*,iostat=ier) lmax_pot
            if( ier > 0 ) call write_err_form(itape4,keyword)
          else
            lmax_pot = lmax_pot_default
          endif

        case('helmholtz','helm_cos')
          if( keyword(6:8) == 'cos' ) Helm_cos = .true.
          n = nnombre(itape4,132)
          if( n < 2 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,'(/A/)') ' After keyword "Helmholtz", 2 or 3 numbers must be given:'
              write(ipr,'(/A/)') '   Helmoltz layer potential, Helmoltz layer distance, and optionnaly Helmoltz potential extension'
            end do
            stop
          endif
          if( n == 2 ) then
            read(itape4,*,iostat=ier) V_helm, Delta_helm
            Width_helm = 2 * Delta_helm
          else
            read(itape4,*,iostat=ier) V_helm, Delta_helm, Width_helm
          endif

        case('optic')
          Optic = .true.

        case('clementi')
          clementi = .true.

        case('scf_mag_f')
          SCF_mag_free = .true.

        case('xan_atom')

        case('spherical')
          spherical_tensor = .true.

        case('sphere_al')
          spherical_tensor = .true.
          spherical_signal = .true.

        case('cartesian')
          cartesian_tensor = .true.

        case('extract','extract_t')
          n = nnombre(itape4,132)
          read(itape4,'(A)') nom_fich_Extract
          if( nom_fich_Extract(1:1) == ' ' ) nom_fich_Extract = adjustl( nom_fich_Extract )
          Length = len_trim(nom_fich_Extract)
          if( nom_fich_Extract(Length-3:Length) == '_bav' ) nom_fich_Extract(Length+1:Length+4) = '.txt'

        case('extractsy')
          n = nnombre(itape4,132)
          n = min(n,n_multi_run_e)
          read(itape4,*) nsymextract(1:n)

        case('extractpo')
          n = nnombre(itape4,132)
          n = min(n,n_multi_run_e)
          read(itape4,*) nposextract(1:n)

        case('trace')
          n = nnombre(itape4,132)
          Coupelapw = .true.
          read(itape4,*,iostat=ier) Trace_k, Trace_p(:)
          if( ier > 0 ) call write_err_form(itape4,keyword)
          if( keyword == 'trace_for' .or. keyword == 'trace_wie' ) trace_format_wien = .true.

        case('range','rangel')
          if( keyword == 'rangel' ) then
            lin_gam = 1
          else
            lin_gam = 0
          endif
          read(itape4,*,iostat=ier) Egamme(1:ngamme)

        case('rixs_ener')

          n = nnombre(itape4,Length_line)
          allocate( Eg(n) )
          read(itape4,*,iostat=ier) Eg(1:n)
          do i = 3,n,2
            if( Eg(i) < Eg(i-2) - eps10 ) exit
          end do
          n = min(i,n)

          Energ_rixs(1) = Eg(1)
          if( nenerg_in_rixs > 1 ) then
            j = 2
            do ie = 2,nenerg_in_rixs
              Energ_rixs(ie) = Energ_rixs(ie-1) + Eg(j)
              if( ie < nenerg_in_rixs ) then
                do i = j, n - 3, 2
                  if( Energ_rixs(ie) > Eg(i+1) - eps10 ) then
                    j = j + 2
                    Energ_rixs(ie) = min( Eg(j-1), Energ_rixs(ie) )
                  else
                    exit
                  endif
                end do
              endif
            end do
          endif
          deallocate( Eg )

        case('energphot')
          energphot = .true.

        case('density')
          Density = .true.

        case('density_c')
          Density = .true.
          Density_comp = .true.

        case('density_a')
          Density = .true.
          State_all = .true.
          State_all_out = .true.

        case('spin_chan')
          Spin_channel = .true.

        case('supermuf')
          supermuf = .true.

        case('two_theta')
          n = nnombre(itape4,132)
          n = min(n,3)
          if( n_two_theta == 1 ) then
            read(itape4,*) Two_theta(1)
          else
            read(itape4,*) Theta_min, Step
            Two_theta(1) = Theta_min
            do i = 2,n_two_theta
              Two_theta(i) = Two_theta(i-1) + Step
            end do
          endif

        case('theta_in')
          n = nnombre(itape4,132)
          n = min(n,3)
          if( n_theta_in == 1 ) then
            read(itape4,*) Theta_in(1)
          else
            read(itape4,*) Theta_min, Step
            Theta_in(1) = Theta_min
            do i = 2,n_theta_in
              Theta_in(i) = Theta_in(i-1) + Step
            end do
          endif

        case('powder')
          Powder = .true.

        case('old_zero')
          Old_zero = .true.

        case('new_zero')
          Old_zero = .false.

        case('e_out_min')
          n = nnombre(itape4,132)
          if( n == 1 ) then
            read(itape4,*,iostat=ier) Eclie_out
            Eclie = Eclie_out
          else
            read(itape4,*,iostat=ier) Eclie, Eclie_out
          endif
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('eneg')
          Eneg_i = .true.
          Eclie = 0._db
          Eclie_out = 0._db

        case('not_eneg')
          Eneg_n_i = .true.

        case('rydberg')
          Rydberg = .true.
          n = nnombre(itape4,132)
          if( n > 0 ) read(itape4,*) R_rydb

        case('nchemin')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) nchemin
          if( ier > 0 ) call write_err_form(itape4,keyword)
          normaltau = .true.

        case('lmoins1')
          lmoins1 = .true.

        case('lplus1')
          lplus1 = .true.

        case('base_reel')
          basereel = .true.

        case('base_comp')
          basereel = .false.

        case('spinorbit')
          Spinorbite = .true.

        case('core_reso')
          core_resolved_e = .true.

        case('no_core_r')
          no_core_resolved = .true.

        case('ang_spin')
          n = nnombre(itape4,132)
          n = min(3,n)
          read(itape4,*,iostat=ier) Ang_spin(1:n)
          if( ier > 0 ) call write_err_form(itape4,keyword)
          Ang_spin(1:n) = Ang_spin(1:n) * pi / 180

        case('axe_spin')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Axe_spin(1:3)
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('rot_sup')
          n = nnombre(itape4,132)
          n = min(3,n)
          read(itape4,*,iostat=ier) Ang_rotsup(1:n)
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('relativis')
          relativiste = .true.

        case('non_relat')
          non_relat = 1

        case('no_dft')
          No_DFT = .true.

        case('no_renorm')
          No_renorm = .true.

        case('renorm')
          Renorm = .true.

        case('no_res_ma')
          read(itape4,*,iostat=ier) f_no_res(1)
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('no_res_mo')
          read(itape4,*,iostat=ier) f_no_res(2)
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('z_nospino')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Z_nospinorbite
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('debye')
          read(itape4,*,iostat=ier) temp
          if( ier > 0 ) call write_err_form(itape4,keyword)

          if(temp > 10000-eps10) then
            call write_error
            do ipr = 6,9,3
              write(ipr,'(/A/)') ' The temperature must be less than 10000K!'
            end do
            stop
          elseif(temp <= - eps10) then
            call write_error
            do ipr = 6,9,3
              write(ipr,'(//A//)') ' Negative temperatures are not allowed!'
            end do
            stop
          end if

        case('ephot_min')
          read(itape4,*,iostat=ier) Ephot_min
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('radius')
          n = nnombre(itape4,132)
          if( n == 1 ) then
            read(itape4,*,iostat=ier) rsorte_s(1)
          else
            read(itape4,*,iostat=ier) rsorte_s(1), ( E_radius(i-1), rsorte_s(i), i = 2,n_radius )
          endif
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('nrato')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) nrato_dirac
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('multrmax')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) multrmax
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('rpotmax')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) rpotmax
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('over_rad')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) roverad
          if( ier > 0 ) call write_err_form(itape4,keyword)
          overad = .true.

        case('screening')
          n = nnombre(itape4,132)
          Select case(n)
            case(1,2)
              read(itape4,*,iostat=ier) ecrantage(1:min(nspin,n))
              necrantage = 0
              lecrantage = 0
            case(3,4)
              n = n - 2
              read(itape4,*,iostat=ier) necrantage, lecrantage, ecrantage(1:min(nspin,n))
          end select
          Force_ecr = .true.
          if( n == 1 .and. Magnetic ) then
            ecrantage(nspin) = ecrantage(1) / nspin
            ecrantage(1) = ecrantage(nspin)
          endif

        case('rixs_sum')
          Sum_rixs = .true.

        case('tddft')
          tddft = .true.

        case('rpalf')
          tddft = .true.
          rpalf = .true.

        case('dyn_eg')   ! noyau TDDFT fxc diag sur seuils
          tddft = .true.
          dyn_eg = .true.

        case('dyn_g')  ! noyau TDDFT fxc diag sur etats initiaux
          tddft = .true.
          dyn_g = .true.

        case('kern_fac')
          Kern_fac_default = .false.
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Kern_fac
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('kern_fast')
          Kern_fast = .true.

        case('gamma_tdd')
          Gamma_tddft = .true.

        case('atomic_sc')
          Atomic_scr = .true.

        case('edge')
          n = nnombre(itape4,132)
          read(itape4,'(A)') mot132
          Seuil = identmot(mot132,3)
          select case( Seuil(1:1) )
            case('k')
              Seuil = 'K1'
            case('l')
              Seuil(1:1) = 'L'
            case('m')
              Seuil(1:1) = 'M'
            case('n')
              Seuil(1:1) = 'N'
            case('o')
              Seuil(1:1) = 'O'
            case('p')
              Seuil(1:1) = 'P'
          end select

        case('quadrupol')
          Quadrupole = .true.

        case('octupole','dipole_oc','dip_oct')
          Octupole = .true.

        case('e1e1')
          E1E1 = .true.

        case('e1e2')
          E1E2e = .true.

        case('e1e3')
          E1E3 = .true.

        case('e2e2')
          E2E2 = .true.

        case('e1m1')
          E1M1 = .true.

        case('e1m2')
          E1M2 = .true.

        case('e3e3')
          E3E3 = .true.

        case('m1m1')
          M1M1 = .true.

        case('m2m2')
          M2M2 = .true.

        case('dipmag')
          dipmag = .true.
          E1M1 = .true.
          M1M1 = .true.

        case('quadmag')
          quadmag = .true.
          E1M2 = .true.
          M2M2 = .true.
          M1M2 = .true.

        case('no_e1e1')
          E1E1 = .false.

        case('no_e2e2')
          no_e2e2 = .true.

        case('no_e1e3')
          no_e1e3 = .true.

        case('no_e1e2')
          no_dipquad = .true.

        case('no_e3e3')
          no_e3e3 = .true.

        case('ldipimp')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) ldipimp(1:3)
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('lquaimp')
          do i = 1,3
            n = nnombre(itape4,132)
            read(itape4,*,iostat=ier) lquaimp(i,1:3)
            if( ier > 0 ) call write_err_form(itape4,keyword)
          end do

        case('normaltau')
          normaltau = .true.

        case('noncentre')
          noncentre = .true.
          State_all = .true.

        case('center')
          noncentre = .true.
          State_all = .true.
          n = nnombre(itape4,132)
          select case(n)
            case(0)
              Centre_auto = .true.
            case(1,2)
              Centre_auto = .true.
              read(itape4,*)
            case default
              read(itape4,*,iostat=ier) Centre(1:3)
              if( ier > 0 ) call write_err_form(itape4,keyword)
          end select

        case('center_s')
          Center_s = .true.
          read(itape4,*,iostat=ier) Centre(1:2)

        case('center_ab')
          noncentre = .true.
          State_all = .true.
          n = nnombre(itape4,132)
          select case(n)
            case(0)
              Centre_auto = .true.
              Centre_auto_abs = .true.
            case(1,2)
              Centre_auto = .true.
              Centre_auto_abs = .true.
              read(itape4,*)
            case default
              read(itape4,*,iostat=ier) Centre(1:3)
          end select
          Centre_auto = .true.

        case('polarized')
          polarise = .true.
          do ipl = 1,nple - n_mat_polar
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            select case(n)
              case(3)
                read(itape4,*,iostat=ier) polar(:,ipl)
              case(4)
                read(itape4,*,iostat=ier) polar(:,ipl), pdpolar(ipl,1)
              case(6)
                read(itape4,*,iostat=ier) polar(:,ipl), veconde(:,ipl)
              case(7)
                read(itape4,*,iostat=ier) polar(:,ipl),veconde(:,ipl), pdpolar(ipl,1)
                pdpolar(ipl,2) = pdpolar(ipl,1)
              case(8)
                read(itape4,*,iostat=ier) polar(:,ipl),veconde(:,ipl), pdpolar(ipl,:)
              case default
                Error_message = ' Wrong number of elements under the card Polarized !'
                call write_error_message(Error_message,6,0)
                stop
            end select
            if( ier > 0 ) call write_err_form(itape4,keyword)
          end do

        case('mat_polar')
          do ipl = nple-n_mat_polar+1,nple
            read(itape4,*,iostat=ier) polar(:,ipl), veconde(:,ipl)
          end do

        case('allsite')
          allsite = .true.

        case('all_nrixs')
          All_nrixs = .true.

        case('lmax_nrix')
          read(itape4,*,iostat=ier) lmax_nrixs

        case('lmax_tddf')
          read(itape4,*,iostat=ier) lmax_tddft_inp

        case('nrixs')
          i = 0
          do
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            read(itape4,*,iostat=ier) q_nrixs(4,i+1:i+n)
            i = i + n
          end do

        case('nrixs_mon')
          do i = 1,nq_nrixs
            n = nnombre(itape4,132)
            if( n == 1 ) then
              read(itape4,*,iostat=ier) q_nrixs(4,i)
            elseif( n == 2 .or. n == 3 ) then
              call write_error
              do ipr = 6,9,3
                write(ipr,100)
                write(ipr,'(/A)') ' After keyword "nrixs_mono" 4 numbers per line are needed to define modulus and direction of q !'
              end do
              stop
            else
              Monocrystal = .true.
              read(itape4,*,iostat=ier) q_nrixs(4,i), q_nrixs(1:3,i)
            endif
          end do

        case('step_azim')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Step_azim
          if( ier > 0 ) call write_err_form(itape4,keyword)
          if( abs( Step_azim ) < eps10 ) then
            nphim = 1
          else
            nphim = nint( 360._db / abs( Step_azim ) )
          endif

        case('symsite')
          symauto = .false.
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) n_atom_proto
          if( ier > 0 ) call write_err_form(itape4,keyword)
          igr = 0
          do ipr = 1,n_atom_proto
            n = nnombre(itape4,132)
            read(itape4,*,iostat=ier) natomsym
            if( ier > 0 ) call write_err_form(itape4,keyword)
            write(iscratch,*) natomsym
            do i = 1,natomsym
              n = nnombre(itape4,132)
              igr = igr + 1
              if( n == 1 ) then
                read(itape4,*,iostat=ier) isymeq
                p(:) = 0._db
              else
                read(itape4,*,iostat=ier) isymeq, p(:)
              endif
              if( ier > 0 ) call write_err_form(itape4,keyword)
              write(iscratch,*) igr, p(:), isymeq
            end do
          end do

        case('dafs')
          Dafs = .true.
          ipl = 0
          do jpl = 1,1000000000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            ipl = ipl + 1
            select case(n)
              case(3)
                read(itape4,*) hkl_dafs(:,ipl)
                nn = nnombre(itape4,132)
                select case(nn)
                  case(6)
                    read(itape4,*,iostat=ier) poldafse(1:3,ipl), vecdafsem(1:3,ipl)
                    read(itape4,*,iostat=ier) poldafss(1:3,ipl), vecdafssm(1:3,ipl)
                  case(9)
                    read(itape4,*,iostat=ier) (poldafse(i,ipl), poldafsei(i,ipl), i = 1,3), vecdafsem(1:3,ipl)
                    read(itape4,*,iostat=ier) (poldafss(1:3,ipl), poldafssi(i,ipl), i = 1,3), vecdafssm(1:3,ipl)
                  case default
                    call write_error
                    do ipr = 6,9,3
                      write(ipr,110) ipl
                    end do
                    stop
                end select
                if( ier > 0 ) call write_err_form(itape4,keyword)
              case(5)
                read(itape4,*,iostat=ier) hkl_dafs(:,ipl), isigpi(:,ipl)
                angpoldafs(3,ipl) = - 10000._db
                do i = 1,2
                  if( isigpi(i,ipl) == 1 ) then
                    angpoldafs(i,ipl) = 0._db
                  elseif( isigpi(i,ipl) == 2 ) then
                    angpoldafs(i,ipl) = 90._db
                  endif
                end do
              case(6)
                read(itape4,*,iostat=ier) hkl_dafs(:,ipl), isigpi(:,ipl), angpoldafs(3,ipl)
                do i = 1,2
                  if( isigpi(i,ipl) == 1 ) then
                    angpoldafs(i,ipl) = 0._db
                  elseif( isigpi(i,ipl) == 2 ) then
                    angpoldafs(i,ipl) = 90._db
                  endif
                end do
              case(7)
                read(itape4,*,iostat=ier) hkl_dafs(:,ipl), (isigpi(i,ipl), angpoldafs(i,ipl), i = 1,2)
                angpoldafs(3,ipl) = - 10000._db
              case(8)
                if( Film ) then
                  read(itape4,*,iostat=ier) hkl_dafs(1:2,ipl), p(1:3), isigpi(:,ipl), angpoldafs(3,ipl)
                  do i = 1,2
                    if( isigpi(i,ipl) == 1 ) then
                      angpoldafs(i,ipl) = 0._db
                    elseif( isigpi(i,ipl) == 2 ) then
                      angpoldafs(i,ipl) = 90._db
                    endif
                  end do
                  hkl_dafs(3,ipl) = p(1)
                  if( abs( p(2) ) < eps10 ) cycle
                  n = nint( ( p(3) - p(1) ) / p(2) )
                  n = max( 0, n )
                  if( n == 0 ) cycle
                  t = ( p(3) - p(1) ) / n
                  do i = 1,n
                    hkl_dafs(3,ipl+i) = hkl_dafs(3,ipl) + i * t
                    hkl_dafs(1:2,ipl+i) = hkl_dafs(1:2,ipl)
                    isigpi(:,ipl+i) = isigpi(:,ipl)
                    angpoldafs(1:3,ipl+i) = angpoldafs(1:3,ipl)
                  end do
                  ipl = ipl + n
                else
                  read(itape4,*,iostat=ier) hkl_dafs(:,ipl), (isigpi(i,ipl), angpoldafs(i,ipl), i = 1,2), angpoldafs(3,ipl)
                endif
              case default
                call write_error
                do ipr = 6,9,3
                   write(ipr,100)
                   write(ipr,120) ipl
                end do
                stop
            end select
            if( ier > 0 ) call write_err_form(itape4,keyword)
          end do

          do ipl = 1,npldafs
            do i = 1,2
              if( isigpi(i,ipl) == 3 .or. isigpi(i,ipl) == 4 .or. isigpi(i,ipl) == 10 ) cycle
              if( abs( angpoldafs(i,ipl) ) < eps10 ) then
                isigpi(i,ipl) = 1
              elseif( abs( angpoldafs(i,ipl) - 90 ) < eps10 ) then
                isigpi(i,ipl) = 2
              else
                isigpi(i,ipl) = 5
              endif
            end do
            do i = 1,2
              if( isigpi(i,ipl) == 10 ) angpoldafs(i,ipl) = -10000._db
            end do
            nscan = 0
            do i = 1,3
              if( angpoldafs(i,ipl) < - 9999._db ) nscan = nscan + 1
            end do
            if( nscan > 1 ) then
              call write_error
              do ipr = 6,9,3
                 write(ipr,100)
                 write(ipr,110) ipl
              end do
              stop
            endif
          end do

        case('phi_0')
          read(itape4,*,iostat=ier) phi_0
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('dafs_2d')
          Operation_mode_used = .true.
          Dafs = .true.
          ipl = 0
          do jpl = 1,npldafs_2d
            n = nnombre(itape4,132)
            n = min(n,13)
            ipl = ipl + 1
            npl_2d(jpl) = 1
            select case(n)
              case(10,11)
                read(itape4,*,iostat=ier) hkl_dafs(:,ipl), Operation_mode(:,jpl), Angle_mode(1:n-8,jpl)
                ipl = ipl + 1
                hkl_dafs(:,ipl) = hkl_dafs(:,ipl-1)
                if( Operation_mode(1,jpl) < 0 .or. Operation_mode(2,jpl) < 0 ) then
                  ipl = ipl + 1
                  hkl_dafs(:,ipl) = hkl_dafs(:,ipl-1)
                  ipl = ipl + 1
                  hkl_dafs(:,ipl) = hkl_dafs(:,ipl-1)
                endif
              case(12,13)
                read(itape4,*,iostat=ier) hkl_dafs(1:2,ipl), p(1:3), Operation_mode(:,jpl), Angle_mode(1:n-10,jpl)
                hkl_dafs(3,ipl) = p(1)
                if( abs( p(2) ) > eps10 ) then
                  n = nint( ( p(3) - p(1) ) / p(2) )
                  n = max( 0, n )
                  if( n > 0 ) then
                    t = ( p(3) - p(1) ) / n
                    ipl0 = ipl
                    do i = 1,n
                      ipl = ipl + 1
                      hkl_dafs(3,ipl) = hkl_dafs(3,ipl0) + i * t
                      hkl_dafs(1:2,ipl) = hkl_dafs(1:2,ipl0)
                    end do
                    npl_2d(jpl) = 1 + n
                  endif
                endif
! Then the same index for the other outgoing polarization
                do i = 1,npl_2d(jpl)
                  ipl = ipl + 1
                  hkl_dafs(:,ipl) = hkl_dafs(:,ipl-npl_2d(jpl))
                end do
                if( Operation_mode(1,jpl) < 0 .or. Operation_mode(2,jpl) < 0 ) then
                  do j = 1,2
                    do i = 1,npl_2d(jpl)
                      ipl = ipl + 1
                      hkl_dafs(:,ipl) = hkl_dafs(:,ipl-npl_2d(jpl))
                    end do
                  end do
                endif
              case default
                call write_error
                do ipr = 6,9,3
                   write(ipr,100)
                   write(ipr,120) ipl
                end do
                stop
            end select
            if( ier > 0 ) call write_err_form(itape4,keyword)
          end do

        case('dafs_exp')
          Dafs = .true.
          do i = 1,3
            read(itape4,*,iostat=ier) Mat_or(i,:)
            if( ier > 0 ) call write_err_form(itape4,keyword)
          end do
          ipl0 = 0
          do i = 1,n_file_dafs_exp
            read(itape4,*,iostat=ier) Angle_dafs_exp(i)
            if( ier > 0 ) call write_err_form(itape4,keyword)
            n = nnombre(itape4,132)
            read(itape4,'(A)') Fichier
            Fichier = Adjustl(Fichier)
            Length = len_trim(Fichier)
            if( Length > 4 ) then
              if( Fichier(Length-3:Length-3) /= '.' ) Fichier(Length+1:Length+4) = '.txt'
            endif
            open(99, file=Fichier, status='old',iostat=istat)
            File_dafs_exp(i) = Fichier
            n = nnombre(99,Length_line)
            n = n / 3
            read(99,*,iostat=ier) (hkl_dafs(:,ipl),ipl=ipl0+1,ipl0+n)
            if( ier > 0 ) call write_err_form(99,keyword)
            Close(99)
            Angle_or(ipl0+1:ipl0+n) = Angle_dafs_exp(i)
            ipl0 = ipl0 + n
          end do

        case('green')
          Green_s = .true.

        case('green_bul')
          Green_bulk = .true.

        case('green_int')
          Green_s = .true.
          Green_int = .true.

        case('zero_azim')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Vec_orig(:)
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('solsing')
          solsing_only = .true.
          solsing_s = .true.

        case('no_solsin')
          No_solsing = .true.

        case('norman')
          normrmt = 2

        case('ray_max_d')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Ray_max_dirac
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('rcharge')
          do it = 1,100000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            if( it <= ntype ) then
              read(itape4,*,iostat=ier) rchimp(it)
              if( ier > 0 ) call write_err_form(itape4,keyword)
            else
              read(itape4,*)
            endif
          end do

        case('rcharge_z')
          do i = 1,Z_Mendeleiev_max
            read(itape4,*,iostat=ier) Z_rchimp(i), Rchimp_Z(i)
            if( ier > 0 ) then
              backspace(itape4)
              exit
            endif
          end do
          n_rchimp_Z = i - 1

        case('rmtg')
          normrmt = 4
          n1 = 1
          do i = 1,ntype
            n = nnombre(itape4,132)
            n2 = min( n1+n-1, ntype )
            if( n2 >= n1 ) then
              read(itape4,*,iostat=ier) rmtimp(n1:n2)
              if( ier > 0 ) call write_err_form(itape4,keyword)
              n1 = n2 + 1
            else
              exit
            endif
          end do
          if( n2 < ntype ) normrmt = - n2
          overlap = 0._db

        case('rmtg_z')
          normrmt = 4
          do i = 1,Z_Mendeleiev_max
            read(itape4,*,iostat=ier) Z_rmtg(i), Rmtg_Z(i)
            if( ier > 0 ) then
              backspace(itape4)
              exit
            endif
          end do
          n_rmtg_Z = i - 1

        case('rmtv0')
          normrmt = 5
          n = nnombre(itape4,132)
          read(itape4,*) v0bdcFimp(1:min(n,nspin))
          if( nspin > n ) v0bdcFimp(nspin) = v0bdcFimp(1)
          korigimp = .true.
          overlap = 0._db

        case('overlap')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) overlap
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('muffintin')
          muffintin = .true.

        case('iord')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) iord
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('adimp')
          n = nnombre(itape4,132)
          adimpin = .true.
          if( n == 1 ) then
            read(itape4,*,iostat=ier) adimp(1)
          else
            read(itape4,*,iostat=ier) adimp(1), ( E_adimp(i-1), adimp(i), i = 2,n_adimp )
          endif
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('rmt')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) rmt(1:min(ntype,n))
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('lmax')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) lmaxat0
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('lmaxfree')
          lmaxfree = .true.

        case('lmaxstden')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) lamstdens
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('lmaxso')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) lmaxso0
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('lmaxso_ma')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) lmaxso_max
          if( ier > 0 ) lmaxso_max = 28

        case('chfree')
          Charge_free = .true.

        case('hedin')
          hedin = .true.

        case('perdew')
          perdew = .true.

        case('xalpha')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) alfpot
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('flapw','flapw_s','flapw_r','flapw_s_p','flapw_psi', 'flapw_n','flapw_n_p','flapw_s_n')

          Matper = .true.
          if( keyword(6:7) == '_n' ) Flapw_new = .true.
          if( keyword(6:7) == '_s' ) then
            Wien_save = 1
            n = nnombre(itape4,132)
            read(itape4,'(A)') Wien_file(9)
          elseif( keyword(6:7) == '_r' ) then
            Wien_save = - 1
            n = nnombre(itape4,132)
            read(itape4,'(A)') Wien_file(9)
          endif
          com(:) = ' Come from FLAPW'
          icom(:) = 3

          do i = 1,8
            if( i == 2 .and. Wien_save == - 1 ) exit
            if( i == 4 .and. .not. ( Flapw_new .and. nspin == 2 ) ) cycle
            if( ( i == 6 .or. i == 7 ) .and. nspin == 1 ) cycle
            if( i == 8 .and. ( keyword == 'flapw_s_p' .or. keyword == 'flapw_psi' .or. keyword == 'flapw_n_p' ) ) cycle
            n = nnombre(itape4,132)
            read(itape4,'(A)') Wien_file(i)
            if( Wien_file(i)(1:1) == ' ' ) Wien_file(i) = Adjustl( Wien_file(i) )
          end do

        case('delta_en_')
          read(itape4,*,iostat=ier) Delta_En_conv
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('scf')
          ier = 0
          SCF = .true.
          n = nnombre(itape4,132)
          if( n == 1 ) then
            read(itape4,*,iostat=ier) nself
          elseif( n == 2 ) then
            read(itape4,*,iostat=ier) nself, p_self0
          elseif( n > 2 ) then
            read(itape4,*,iostat=ier) nself, p_self0, Delta_En_conv
          endif
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('no_scf_bu')
          No_SCF_bulk = .true.

        case('scf_abs')
          scf_elecabs = .true.

        case('p_self')
          read(itape4,*,iostat=ier) p_self0
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('p_self_ma')
          read(itape4,*,iostat=ier) p_self_max
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('n_self')
          read(itape4,*,iostat=ier) nself
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('r_self')
          read(itape4,*,iostat=ier) r_self
          if( ier > 0 ) call write_err_form(itape4,keyword)
          r_self_imp = .true.

        case('scf_exc')
          SCF_exc_imp = .true.
          SCF = .true.

        case('scf_non_e')
          Self_nonexc_imp = .true.
          SCF = .true.

        case('nonexc')
          Nonexc = .true.
          nonexc_imp = .true.

        case('excited')
          Exc_imp = .true.

        case('no_fermi')
          Fermi_auto = .false.

        case('hubbard')
          n1 = 1
          do i = 1,100
            n = nnombre(itape4,132)
            n2 = min( n1+n-1, ntype)
            if( n2 >= n1 ) then
              read(itape4,*,iostat=ier) V_hubbard(n1:n2)
              if( ier > 0 ) call write_err_form(itape4,keyword)
              do j = n1,n2
               if( abs(V_hubbard(j)) > eps10 ) Hubb(j) = .true.
              end do
              n1 = n2 + 1
            else
              exit
            endif
          end do

        case('hubbard_z')
          do i = 1,Z_Mendeleiev_max
            read(itape4,*,iostat=ier) Z_hubbard(i), Hubbard_Z(i)
            if( ier /= 0 ) then
              backspace(itape4)
              exit
            endif
          end do
          n_hubbard_Z = i - 1

        case('full_atom')
          Full_atom_e = .true.

        case('absorbeur')
          Absorber = .true.
          k = 0
          do i = 1,100000
            n = nnombre(itape4,132)
            if( n < 1 ) exit
            read(itape4,*,iostat=ier) iabsm(k+1:k+n)
            if( ier > 0 ) call write_err_form(itape4,keyword)
            k = k + n
          end do

        case('atom_conf')

          Atom_conf = .true.

          do it = 1,ntype_conf
            read(itape4,*) nb_atom_conf(it), igra(1:nb_atom_conf(it),it), nlat(it)
            backspace(itape4)
            n = nnombre(itape4,132)
            if( n == nb_atom_conf(it) + 2 + 3*nlat(it) ) then
              read(itape4,*,iostat=ier) nb_atom_conf(it), igra(1:nb_atom_conf(it),it), nlat(it), ( nvval(it,l), lvval(it,l), &
                             popval(it,l,1), l = 1,nlat(it) )
              if( ier > 0 ) call write_err_form(itape4,keyword)
              if( nspin == 2 ) then
                do l = 1,nlat(it)
                  popval(it,l,1) = 0.5_db * popval(it,l,1)
                  popval(it,l,2) = popval(it,l,1)
                end do
              endif
            elseif( n == nb_atom_conf(it) + 2 + 4*nlat(it) ) then
              if( nspin == 1 ) then
                allocate( x(nlat(it)) )
                read(itape4,*,iostat=ier) nb_atom_conf(it), igra(1:nb_atom_conf(it),it), nlat(it),  &
                                         ( nvval(it,l), lvval(it,l), popval(it,l,1), x(l), l = 1,nlat(it) )
                if( ier > 0 ) call write_err_form(itape4,keyword)
                do l = 1,nlat(it)
                  popval(it,l,1) = popval(it,l,1) + x(l)
                end do
                deallocate( x )
              else
                read(itape4,*,iostat=ier) nb_atom_conf(it), igra(1:nb_atom_conf(it),it), nlat(it), &
                                         ( nvval(it,l), lvval(it,l), popval(it,l,:), l = 1,nlat(it) )
                if( ier > 0 ) call write_err_form(itape4,keyword)
              endif
            else
              call write_error
              do ipr = 6,9,3
                write(ipr,100)
                write(ipr,125)  it
              end do
              stop
            endif
          end do

        case('atom')
          atom = .true.
          it = 0
          do jt = 1,100000
            n = nnombre(itape4,132)

            if( n == 0 ) exit

            if( n == 3 ) then
              read(itape4,*,iostat=ier) Ang_base_loc(:,it)
              if( ier > 0 ) call write_err_form(itape4,keyword)
              n = nnombre(itape4,132)
            endif

            nrm = max( nrm, nrato_dirac )
            if( jt <= ntype ) it = it + 1
            if( n == 1 ) then
              read(itape4,*,iostat=ier) numat(it)
              if( ier > 0 ) call write_err_form(itape4,keyword)
              nlat(it) = 0
            else
              read(itape4,*,iostat=ier) numat(it), nlat(it)
              if( ier > 0 ) call write_err_form(itape4,keyword)
              if( nlat(it) > 0 ) then
                backspace(itape4)
                n = nnombre(itape4,132)
                if( n == 2 + 3*nlat(it) ) then
                  read(itape4,*,iostat=ier) numat(it), nlat(it), ( nvval(it,l), lvval(it,l), &
                                 popval(it,l,1), l = 1,nlat(it) )
                  if( nspin == 2 ) then
                    do l = 1,nlat(it)
                      popval(it,l,1) = 0.5_db * popval(it,l,1)
                      popval(it,l,2) = popval(it,l,1)
                    end do
                  endif
                elseif( n == 2 + 4*nlat(it) ) then
                  if( nspin == 1 ) then
                    allocate( x(nlat(it)) )
                    read(itape4,*,iostat=ier) numat(it),nlat(it), ( nvval(it,l), lvval(it,l), popval(it,l,1), x(l), &
                                                                    l = 1,nlat(it) )
                    do l = 1,nlat(it)
                      popval(it,l,1) = popval(it,l,1) + x(l)
                    end do
                    deallocate( x )
                  else
                    read(itape4,*,iostat=ier) numat(it), nlat(it), ( nvval(it,l), &
                       lvval(it,l), popval(it,l,:), l = 1,nlat(it) )
                  endif
                else
                  call write_error
                  do ipr = 6,9,3
                    write(ipr,100)
                    write(ipr,130) it
                  end do
                  stop
                endif
                if( ier > 0 ) call write_err_form(itape4,keyword)
              endif
            endif

          end do

        case('atom_nsph')

          do
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            n = n - 2
            allocate( igr_nsph(n) )
            read(itape4,*) nb_atom_nsph, igr_nsph(:), norb
            norbv( igr_nsph(:) ) = norb
            do io = 1,norb
              n = nnombre(itape4,132)
              n = n - 1
              allocate( Hyb(n) )
              read(itape4,*) Hyb(1:n), pp
              do jgr = 1,nb_atom_nsph
                igr = igr_nsph(jgr)
                pop_nonsph(io,igr) = pp
                select case(n)
                  Case(1,2,3,4,5,6,7,8,9,16)
                    Hybrid(io,1:n,igr) = cmplx( Hyb(1:n), 0._db, db )
                  Case default
                    do i = 1,n-1,2
                      Hybrid(io,(i+1)/2,igr) = cmplx( Hyb(i), Hyb(i+1), db )
                    end do
                end select
              end do
              deallocate ( Hyb )
            end do
            deallocate( igr_nsph )
          end do

        case('dilatorb')
          do i = 1,norbdil
            n = nnombre(itape4,132)
            read(itape4,*,iostat=ier) itdil(i), ldil(i), cdil(i)
            if( ier > 0 ) call write_err_form(itape4,keyword)
          end do

        case('v0imp')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) v0bdcFimp(1:min(n,nspin))
          if( ier > 0 ) call write_err_form(itape4,keyword)
          if( nspin > n ) v0bdcFimp(nspin) = v0bdcFimp(1)
          korigimp = .true.

        case('vmax')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) v_intmax
          if( ier > 0 ) call write_err_form(itape4,keyword)

        case('eimag')
          do ie = 1,neimagent
            n = nnombre(itape4,132)
            if( n == 1 ) then
              read(itape4,*,iostat=ier) eimagent(ie)
              eeient(ie) = 0._db
            else
              read(itape4,*,iostat=ier) eeient(ie), eimagent(ie)
            endif
            if( ier > 0 ) call write_err_form(itape4,keyword)
          end do

        case('pointgrou')
          n = nnombre(itape4,132)
          read(itape4,'(A)') mot132
          PointGroup = identmot(mot132,8)
          select case( PointGroup(1:1) )
            Case('_')
              PointGroup(1:1) = '-'
            Case('c')
              PointGroup(1:1) = 'C'
            Case('d')
              PointGroup(1:1) = 'D'
            Case('s')
              PointGroup(1:1) = 'S'
            Case('o')
              PointGroup(1:1) = 'O'
            Case('t')
              PointGroup(1:1) = 'T'
          end select
          PointGroup_Auto = .false.

        case('symmol')
          Symmol = .true.

        case('spgroup')
          n = nnombre(itape4,132)
          read(itape4,'(A)') mot132
          if( mot132(1:1) == ' ' ) mot132 = adjustl(mot132)
          Space_group = mot132(1:13)

! Cluster, film or unit cell description
        case('crystal','molecule','crystal_t','molecule_','film','film_t','surface','surface_t','interface','interfac_')

          if( keyword(1:8) /= 'molecule' ) Matper = .true.

          n = nnombre(itape4,132)
          if( n == 3 ) then
            read(itape4,*) axe(1:3)
            ang(1:3) = 90._db
          elseif( n == 2 ) then
            cylindre = .true.
            read(itape4,*) axe(1:3:2)
            axe(2) = axe(1)
            ang(1:3) = 90._db
          elseif( n == 1 ) then
            spherique = .true.
            read(itape4,*) axe(1)
            axe(2) = axe(1)
            axe(3) = axe(1)
            ang(1:3) = 90._db
          else
            read(itape4,*,iostat=ier) axe(1:3), ang(1:3)
          endif
          if( ier > 0 ) call write_err_form(itape4,keyword)
          if( keyword(1:1) == 's' ) then
            angxyz_sur(:) = ang(:)
            axyz_sur(:) = axe(:)
          elseif( keyword(1:1) == 'i' ) then
            angxyz_int(:) = ang(:)
            axyz_int(:) = axe(:)
          else
            angxyz(:) = ang(:)
            axyz(:) = axe(:)
          endif
          do igr = 1,n_atom_neq
            if( keyword(1:1) == 'i' ) then
              if( igr <= n_atom_per_neq .or. igr > n_atom_per_neq + n_atom_int) cycle
            elseif( keyword(1:1) == 's' ) then
              if( igr <= n_atom_per_neq + n_atom_int ) cycle
            elseif( igr > n_atom_per_neq ) then
              cycle
            endif
            if( Temp_B_iso .and. Taux ) then
              n = nnombre(itape4,132)
              if( n > 5 ) then
                if( Occupancy_first ) then
                  read(itape4,*) itype(igr), posn(:,igr), Taux_oc(igr), Temp_coef(igr)
                else
                  read(itape4,*) itype(igr), posn(:,igr), Temp_coef(igr), Taux_oc(igr)
                endif
              elseif( n == 5 ) then
                if( Occupancy_first ) then
                  read(itape4,*) itype(igr), posn(:,igr), Taux_oc(igr)
                else
                  read(itape4,*) itype(igr), posn(:,igr), Temp_coef(igr)
                endif
              else
                read(itape4,*) itype(igr), posn(:,igr)
              endif
            elseif( Temp_B_iso ) then
              n = nnombre(itape4,132)
              if( n > 4 ) then
                read(itape4,*) itype(igr), posn(:,igr), Temp_coef(igr)
              else
                read(itape4,*) itype(igr), posn(:,igr)
              endif
            elseif( Taux ) then
              n = nnombre(itape4,132)
              if( n > 4 ) then
                read(itape4,*) itype(igr), posn(:,igr), Taux_oc(igr)
              else
                read(itape4,*) itype(igr), posn(:,igr)
              endif
            elseif( Readfast .or. .not.( Atom_nonsph .or. Atom_occ_hubb .or. Axe_loc) ) then
              read(itape4,*) itype(igr), posn(:,igr)
            else
              n = nnombre(itape4,132)
              if( n == 0 ) exit
              if( n == 2 .or. n == 3 ) then
                read(itape4,*) Ang_base_loc_gr(1:n,igr)
                Ang_base_loc_gr(1:n,igr) = Ang_base_loc_gr(1:n,igr) * pi / 180
                n = nnombre(itape4,132)
              endif
              select case(n)
                case(4)
                  read(itape4,*) itype(igr), posn(:,igr)
                  norb = 0
                case(5)
                  read(itape4,*) itype(igr), posn(:,igr), norb
              end select
              if( norb == 0 ) cycle
              if( norb /= -1 ) then
                if( norb < 0 ) then
                  norbv(igr) = - norb - 1
                else
                  norbv(igr) = norb
                endif
                do io = 1,norbv(igr)
                  n = nnombre(itape4,132)
                  n = n - 1
                  allocate( Hyb(n) )
                  read(itape4,*) Hyb(1:n), pop_nonsph(io,igr)
                  Hybrid(io,1:n,igr) = cmplx( Hyb(1:n), 0._db, db )
                  deallocate( Hyb )
                end do
              endif
              if( norb < 0 ) then
                n = nnombre(itape4,132)
                l = (n/nspin - 1)/2
                read(itape4,*) ((occ_hubb_e(m,m,isp,igr), m = -l,l ), isp = 1,nspin)
              endif
            endif
          end do

! Unit cell reading
        case('pdb_file','film_pdb_')

          Pdb = .true.
          Fichier = ' '
          read(itape4,'(A)') Fichier
          Fichier = Adjustl(Fichier)
          Length = len_trim(Fichier)
          if( Length > 4 ) then
            if( Fichier(Length-3:Length-3) /= '.' ) Fichier(Length+1:Length+4) = '.pdb'
          endif
          open(8, file = Fichier, status='old', iostat=istat)
          if( istat /= 0 ) call write_open_error(Fichier,istat,1)

          Matper = .true.

          igr = 0

          do

            read(8,'(a6)',iostat=eoff) mot6
            if( eoff /= 0 ) exit

            select case(mot6)

              case('END   ')

                exit

              case('SCALE1' )

                backspace(8)
                do i = 1,3
                  read(8,'(10x,3f10.6)') Mat(i,:)
                end do

              case('CRYST1' )

                backspace(8)
                read(8,'(6x,3f9.3,3f7.2,a13)') axyz(1:3), angxyz(1:3), Spgr

                Length = len_trim(Spgr)
                j = 0
                do i = 1,Length
                  if( Spgr(i:i) == ' ' ) cycle
                  j = j + 1
                  Space_group(j:j) = Spgr(i:i)
                end do

              case('ATOM  ','HETATM')

                backspace(8)
                read(8,'(16x,a1,13x,3f8.3,2f6.2,10x,a2)') Let ,p(:), t, tc,  Symbol

                igr = igr + 1
                p = Matmul( Mat, p )
                posn(:,igr) = p(:)
                if( Taux ) Taux_oc(igr) = t
                if( Temp_B_iso ) Temp_coef(igr) = tc

                select case(Let)
                  case(' ')
                    Kgroup(igr) = 0
                  case('A')
                    Kgroup(igr) = 1
                  case('B')
                    Kgroup(igr) = 2
                  case('C')
                    Kgroup(igr) = 3
                  case('D')
                    Kgroup(igr) = 4
                  case default
                    Kgroup(igr) = 5
                end select

                Symbol = adjustl( Symbol )
                do i = 1,103
                  if( Chemical_Symbol_c(i) /= Symbol .and. Chemical_Symbol(i) /= Symbol ) cycle
                  itype(igr) = i
                  exit
                end do
                if( i == 104 ) then
                  call write_error
                  do ipr = 6,9,3
                    write(ipr,140) igr, Symbol
                  end do
                  stop
                endif

              case default

                cycle

              end select

          end do

          Close(8)

! Unit cell reading
        case('cif_file','film_cif_')

          Matper = .true.
          Fcif = .true.

          Fcif_file = ' '
          read(itape4,'(A)') Fcif_file
          Fcif_file = Adjustl(Fcif_file)
          Length = len_trim(Fcif_file)
          if( Length > 4 ) then
            if( Fcif_file(Length-3:Length-3) /= '.' ) Fcif_file(Length+1:Length+4) = '.cif'
          endif
          open(8, file = Fcif_file, status='old', iostat=istat)
          if( istat /= 0 ) call write_open_error(Fcif_file,istat,1)

          do

            read(8,'(A)',iostat=eof) mot

            if( eof /= 0 ) exit

            Length = len_trim(mot)
            do i = 1,Length
              if( mot(i:i) == char(9) ) mot(i:i) = ' '
            end do
            mot = adjustl(mot)

            if( mot(2:30) == 'symmetry_space_group_name_H-M' .or. mot(2:25) == 'space_group_name_H-M_alt' ) then

              Length = len_trim(mot)

              Apostrophe = .true.
              do i = 31,Length
               if( mot(i:i) == '''' ) exit
               if( mot(i:i) /= ' ' ) then
                 Apostrophe = .false.
                 exit
               endif
              end do
              if( .not. Apostrophe ) i = i - 1
              k = 0
              j = i + 1
              do i = j,Length
                if( mot(i:i) == ' ' ) then
                  if( Apostrophe ) then
                    cycle
                  else
                    exit
                  endif
                endif
                if( mot(i:i) == '''' ) exit
                k = k + 1
                Space_group(k:k) = mot(i:i)
              end do

            elseif( mot(2:17) == 'cell_angle_alpha' ) then
              angxyz(1) = number_from_text(1,mot)

            elseif( mot(2:16) == 'cell_angle_beta' ) then
              angxyz(2) = number_from_text(1,mot)

            elseif( mot(2:17) == 'cell_angle_gamma' ) then
              angxyz(3) = number_from_text(1,mot)

            elseif( mot(2:14) == 'cell_length_a') then
              axyz(1) = number_from_text(1,mot)

            elseif( mot(2:14) == 'cell_length_b') then
              axyz(2) = number_from_text(1,mot)

            elseif( mot(2:14) == 'cell_length_c') then
              axyz(3) = number_from_text(1,mot)

            elseif( mot(2:16) == 'atom_site_label' .or. mot(2:22) == 'atom_site_type_symbol' &
               .or. mot(2:16) == 'atom_site_fract' ) then

              backspace(8)
              n_adp_type = 0
              n_B_iso = 0
              n_fract_x = 0; n_fract_y = 0; n_fract_z = 0
              n_label = 0
              n_symbol = 0

              do n = 1,100000
                read(8,'(A)') mot
                mot = adjustl(mot)
                if( mot(1:1) /= '_' ) then
                  backspace(8)
                  exit
                elseif( mot(2:16) == 'atom_site_label') then
                  n_label = n - 1
                elseif( mot(2:22) == 'atom_site_type_symbol') then
                  n_symbol = n - 1
                elseif( mot(2:22) == 'atom_site_adp_type') then
                  n_adp_type = n - 1
                elseif( mot(2:18) == 'atom_site_fract_x') then
                  n_fract_x = n - 1
                elseif( mot(2:18) == 'atom_site_fract_y') then
                  n_fract_y = n - 1
                elseif( mot(2:18) == 'atom_site_fract_z') then
                  n_fract_z = n - 1
                elseif( mot(2:20) == 'atom_site_occupancy') then
                  n_occupancy = n - 1
                elseif( mot(2:25) == 'atom_site_U_iso_or_equiv') then
                  Atom_U_iso = .true.
                  n_B_iso = n - 1
                elseif( mot(2:25) == 'atom_site_B_iso_or_equiv') then
                  Atom_B_iso = .true.
                  n_B_iso = n - 1
                endif
              end do

              igr = 0
              Loop_kgr: do kgr = 1,n_atom_per_neq

                read(8,'(A)') mot
                mot = adjustl(mot)

! Reading chemical element
                if( n_symbol == 0 ) then
                  mot132 = word_from_text(n_label,mot)
                else
                  mot132 = word_from_text(n_symbol,mot)
                endif
                Symbol = mot132(1:2)
                if( Symbol(2:2) == '1' .or. Symbol(2:2) == '2' .or. Symbol(2:2) == '3' .or. Symbol(2:2) == '4' .or. &
                    Symbol(2:2) == '5' .or. Symbol(2:2) == '6' .or. Symbol(2:2) == '7' .or. Symbol(2:2) == '8' .or. &
                    Symbol(2:2) == '9' .or. Symbol(2:2) == '-' .or. Symbol(2:2) == '+' ) Symbol(2:2) = ' '

                do i = 1,103
                  if( Chemical_Symbol_c(i) /= Symbol .and. Chemical_Symbol(i) /= Symbol ) cycle
                  it = i
                  exit
                end do
                if( i == 104 ) then
                  call write_error
                  do ipr = 6,9,3
                    write(ipr,145) kgr, Symbol
                  end do
                  stop
                endif

                p(1) = number_from_text(n_fract_x,mot)
                p(2) = number_from_text(n_fract_y,mot)
                p(3) = number_from_text(n_fract_z,mot)

                do i = 1,3
                  if( p(i) < - eps10 ) p(i) = p(i) + 1._db
                  if( p(i) > 1 - eps10 ) p(i) = p(i) - 1._db
                end do
! Checking for the fake vista cif_file
                do jgr = 1,igr
                  if( sum( abs( posn(:,jgr) - p(:) ) ) < eps10 .and. it == itype(jgr) ) cycle loop_kgr
                end do

                igr = igr + 1
                itype(igr) = it
                posn(:,igr) = p(:)
                if( Taux ) Taux_oc(igr) = number_from_text(n_occupancy,mot)
                if( Temp_B_iso ) Temp_coef(igr) = number_from_text(n_B_iso,mot)

              end do loop_kgr

              exit

            endif
          end do
          Close(8)

        case('dpos')
          n = nnombre(itape4,132)
          read(itape4,*) dpos(1:3)

        case('test_dist')
          n = nnombre(itape4,132)
          read(itape4,*) Test_dist_min

        case('ylm_comp')
          Ylm_comp_inp = .true.

        case('delta_eps')
          n = nnombre(itape4,132)
          read(itape4,*) Delta_Epsii

        case('one_run')
          One_run = .true.
          State_all = .true.

        case('one_scf')
          One_scf = .true.

        case('mat_ub')
          m = 0
          allocate( x(9) )
          do
            n = nnombre(itape4,132)
            if( n == 0 ) call write_err_form(itape4,keyword)
            n = min( m + n, 9 )
            read(itape4,*,iostat=ier) x(m+1:n)
            if( ier > 0 ) call write_err_form(itape4,keyword)
            if( n == 9 ) exit
            m = m + 3
          end do
          Mat_UB(1,1:3) = x(1:3)
          Mat_UB(2,1:3) = x(4:6)
          Mat_UB(3,1:3) = x(7:9)
          deallocate( x )

        case('moment_co')
          Moment_conservation = .true.

        case('rixs')
          RIXS = .true.
          do i = 1,n_q_rixs
            n = nnombre(itape4,132)
            if( n < 3 ) then
              call write_err_form(itape4,keyword)
            elseif( n == 3 ) then
              read(itape4,*) q_rixs(1:3,i)
              q_rixs(4,i) = 0._db
            else
              read(itape4,*) q_rixs(:,i)
            endif
          end do

        case('rixs_core')
          RIXS_core = .true.

        case('rixs_fast')
          RIXS_fast = .true.

        case('step_loss')
          n = nnombre(itape4,132)
          if( n > 0 ) read(itape4,*) Step_loss

        case('setaz')
          read(itape4,*) hkl_ref(:)

        case('z_absorbe')
          n = nnombre(itape4,132)
          if( Absauto ) then
            j = 0
            do i = 1,100000
              n = nnombre(itape4,132)
              if( n < 1 ) exit
              read(itape4,*) numat_abs(j+1:j+n)
              j = j + n
            end do
          else
            read(itape4,*)
          endif

        case('d_max_pot')
          n = nnombre(itape4,132)
          read(itape4,*) D_max_pot

        case('scf_step')
          n = nnombre(itape4,132)
          Pas_SCF_imp = .true.
          read(itape4,*) Pas_SCF

        case('write_mod')
          Write_modQ = .true.

! Parameters already red in lectdim
        case('full_self','magnetism','memory_sa','occupancy','self_abs','readfast')

        case default

          if( igrdat == 1 ) then
            comt = mot
          elseif( keyword(1:1) /= ' ' ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,100)
              write(ipr,150) mot
            end do
            stop
          endif

      end select

    end do boucle_lect

! End of lecture.

! Talkative ("bav" for bavard in French !) file opening

    Length = len_trim(nomfich)
    write(6,'(/a9,A)') ' Filout: ',nomfich(1:Length)
    if( Length > 4 ) then
      if( nomfich(Length-3:Length) == '.txt' .or. nomfich(Length-3:Length) == '.dat') nomfich(Length-3:Length) = '    '
    endif
    nomfichbav = nomfich
    long = len_trim(nomfich)
    nomfichbav(long+1:long+8) = '_bav.txt'

    if( Extract .and. .not. Check_extract ) then
      do i = 1,30
        if( icheck(i) > 1 ) exit
      end do
      if( i > 30 ) icheck(1:30) = 0
    endif

    i = sum( icheck(1:28) )
    if( i > 0 ) then
      if (FDMX_only) then
        open(3, file = nomfichbav, status='old',iostat=istat) !*** JDB March 2018 fix
      else
        open(3, file = nomfichbav, status='unknown',iostat=istat)
      endif
      if( istat /= 0 ) call write_open_error(nomfichbav,istat,1)
    endif

    if( icheck(1) > 0 ) then
      call date_and_time( date = dat, time = tim )
      com_date = '   Date = ' // dat(7:8) // ' ' // dat(5:6) // ' ' // dat(1:4)
      com_time = '   Time = ' // tim(1:2)  // ' h ' // tim(3:4) // ' mn ' // tim(5:6) // ' s'
      write(3,'(1x,A/A/A)') Revision, com_date, com_time
      if( comt /= ' ') write(3,'(/A)') comt
      ipr0 = 3
    else
      ipr0 = 6
    endif

! For potential coming from Wien2k, the structure is read just here
    if( Flapw ) then
      nrm = nrato_dirac
      call lect_struct_lapw(angxyz,axyz,iabsm,iabsorig,icheck(1),its_lapw,itype,n_multi_run_e,n_atom_per,ngroup_lapw,Wien_file(1), &
                            nrato_lapw,nrm,nslapwm,ntype,ntype_bulk,numat,posn,r0_lapw,rlapw,rotloc_lapw,Wien_matsym,Wien_taulap)
      if( normrmt == 1 ) normrmt = 2
    elseif( nrm == 0 ) then
      nrm = nrato_dirac
    endif

!*** JDB
    if( .not. adimpin .and. Use_FDMX ) then
      adimp(1:5) = (/ 0.24_db, 0.20_db, 0.16_db, 0.12_db, 0.08_db /)
      E_adimp(1:4) = (/ 100.0_db, 250.0_db, 400.0_db, 500.0_db /)
    endif
!*** JDB

    if( .not. Extract ) then

      if( n_range > 1 ) then
        i = 1
        j = 1
        do i_range = 1,n_range-1
          if( i < n_adimp .and. j < n_radius ) then
            E_max_range(i_range) = min( E_adimp(i), E_radius(j) )
            if( abs( E_adimp(i) - E_radius(j) ) < eps10 ) then
              i = i + 1
              j = j + 1
            elseif( E_adimp(i) < E_radius(j) ) then
              i = i + 1
            else
              j = j + 1
            endif
          elseif( i < n_adimp ) then
            E_max_range(i_range) = E_adimp(i)
            i = i + 1
          elseif( j < n_radius ) then
            E_max_range(i_range) = E_radius(j)
            j = j + 1
          endif
        end do
      endif
    endif

    if( Doping ) then
      n = ngroup - n_atom_bulk ! the doping element is set at the end of crystal atoms but before the bulk atoms when existing
      Ang_base_loc_gr(:,n) = Ang_base_loc_gr(:,igr_dop)
      Atom_nsph_e(n) = Atom_nsph_e(igr_dop)
      if( Atom .or. Atom_conf ) then
        itype(n) = itype_dop
        if( Atom .and. itype_dop > ntype ) then
          call write_error
          do ipr = ipr0,9,3
            write(ipr,100)
            write(ipr,155) itype_dop, ntype
          end do
          stop
        endif
        if( Atom ) then
          popats(n,:,:) =  popval(itype_dop,:,:)
        else
          popats(n,:,:) =  popval(1,:,:)
        endif
      else
        itype(n) = ntype
        popats(n,:,:) =  popval(itype(n),:,:)
        if( .not. Atom_conf ) itype(n) = itype_dop
      endif
      if( ngroup_nonsph > 1 ) then
        norbv(n) = norbv(igr_dop)
        Hybrid(:,:,n) = Hybrid(:,:,igr_dop)
        pop_nonsph(:,n) = pop_nonsph(:,igr_dop)
      endif
      if( ngroup_m > 0 ) then
        Axe_atom_gr(:,n) = Axe_atom_gr(:,igr_dop)
        Rot_atom_gr(:,:,n) = Rot_atom_gr(:,:,igr_dop)
      endif
      if( ngroup_taux > 0 ) Taux_oc(n) = Taux_oc(igr_dop)
      if( ngroup_hubb > 0 ) occ_hubb_e(:,:,:,n) = occ_hubb_e(:,:,:,igr_dop)

      if( ngroup_pdb > 0 ) Kgroup(n) = Kgroup(igr_dop)
      if( ngroup_temp > 0 ) Temp_coef(n) = Temp_coef(igr_dop)
    endif

! Atom defined with few digits
    if( abs( Angxyz(3) - 120._db ) < 0.0001_db ) then
      do i = 1,11
        if( i == 3 .or. i == 6 .or. i == 9 ) cycle
        Where( abs( Posn - i / 12._db ) < 0.00005_db ) Posn = i / 12._db
      end do
    endif
    if( Bulk ) then
      if( abs( Angxyz_bulk(3) - 120._db ) < 0.0001_db ) then
        do i = 1,11
          if( i == 3 .or. i == 6 .or. i == 9 ) cycle
          Where( abs( Posn_bulk - i / 12._db ) < 0.00005_db ) Posn_bulk = i / 12._db
        end do
      endif
    endif

! Modification when fitting
    if( fit_cal ) then
      do igr = 2,ngroup_par
        istop = 0
        do ipar = 1,npar(igr)
          if( typepar(igr,ipar) /= 'dposx' .and. typepar(igr,ipar) /= 'dposy' .and. typepar(igr,ipar) /= 'dposz' .and. &
              typepar(igr,ipar) /= 'posx' .and. typepar(igr,ipar) /= 'posy' .and. typepar(igr,ipar) /= 'posz' .and. &
              typepar(igr,ipar) /= 'theta' .and. typepar(igr,ipar) /= 'phi' .and. typepar(igr,ipar) /= 'occup'  ) cycle
          if( indice_par(igr,ipar) > ngroup - n_atom_bulk ) then
            call write_error
            do ipr = ipr0,9,3
              write(ipr,100)
              write(ipr,160) typepar(igr,ipar), indice_par(igr,ipar), ngroup - n_atom_bulk
            end do
            istop = 1
           endif
        end do
      end do
      if( istop == 1 ) stop

      do i = 2,ngroup_par
        do ip = 1,npar(i)
          select case( typepar(i,ip) )
            case('occup')
              Taux_oc(indice_par(i,ip)) = param(i,ip)
            case('dposx')
              posn(1,indice_par(i,ip)) = posn(1,indice_par(i,ip)) + param(i,ip)
            case('dposy')
              posn(2,indice_par(i,ip)) = posn(2,indice_par(i,ip)) + param(i,ip)
            case('dposz')
              posn(3,indice_par(i,ip)) = posn(3,indice_par(i,ip)) + param(i,ip)
            case('posx')
              posn(1,indice_par(i,ip)) = param(i,ip)
            case('posy','theta')
              posn(2,indice_par(i,ip)) = param(i,ip)
            case('posz','phi')
              posn(3,indice_par(i,ip)) = param(i,ip)
            case('abc')
              axyz(1:3) = axyz(1:3) * (1 + 0.01*param(i,ip))
            case('a')
              axyz(1) = axyz(1) * (1 + 0.01*param(i,ip))
            case('b')
              axyz(2) = axyz(2) * (1 + 0.01*param(i,ip))
            case('c')
              axyz(3) = axyz(3) * (1 + 0.01*param(i,ip))
            case('angx','anga')
              angxyz(1) = param(i,ip)
            case('angy','angb')
              angxyz(2) = param(i,ip)
            case('angz','angc')
              angxyz(3) = param(i,ip)
            case('poporb')
              io = 0
              boucle_it: do it = 1,ntype
                do l = 1,nlat(it)
                  do ispin = 1,nspin
                    io = io + 1
                    if( io /= indice_par(i,ip) ) cycle
                    popval(it,l,ispin) = param(i,ip)
                    exit boucle_it
                  end do
                end do
              end do boucle_it
              if( it > ntype ) then
                call write_error
                do ipr = ipr0,9,3
                  write(ipr,100)
                  write(ipr,170) typepar(i,ip), indice_par(i,ip)
                end do
                stop
              endif
          end select
        end do
      end do
    endif

    if( cylindre ) then
      do igr = 1,n_atom_uc
        r = posn(1,igr)
        theta = pi * posn(2,igr) / 180
        posn(1,igr) = r * cos( theta )
        posn(2,igr) = r * sin( theta )
      end do
    elseif( spherique ) then
      do igr = 1,n_atom_uc
        r = posn(1,igr)
        theta = pi * posn(2,igr) / 180
        phi = pi * posn(3,igr) / 180
        posn(1,igr) = r * sin( theta ) * cos( phi)
        posn(2,igr) = r * sin( theta ) * sin( phi)
        posn(3,igr) = r * cos( theta )
      end do
    endif

    if( Atom_conf ) then

      do it = 1,ntype_conf
        if( igra(1,it) == 0 ) then
          numat(it) = itype(ngroup)    ! Case of atom_conf on doping atom
          itype_dop = 1
        else
          numat(it) = itype(igra(1,it))
        endif
      end do

      allocate( itZ(103) )

      jt = ntype_conf

      boucle_2: do igr = 1,n_atom_neq_min+1
        if( ( .not. Doping .or. ( Doping .and. itype_dop == 1 ) ) .and. igr == n_atom_neq_min+1 ) cycle
        do it = 1,ntype_conf
          do kgr = 1,nb_atom_conf(it)
            if( igr == igra(kgr,it) ) cycle boucle_2
          end do
        end do
        boucle_3: do jgr = 1,igr-1
          do it = 1,ntype_conf
            do kgr = 1,nb_atom_conf(it)
              if( jgr == igra(kgr,it) ) cycle boucle_3
            end do
          end do
          if( itype(igr) == itype(jgr) ) cycle boucle_2
        end do boucle_3
        jt = jt + 1
        if( Doping .and. igr == n_atom_neq_min+1 ) then
          numat(jt) = itype_dop
        else
          numat(jt) = itype(igr)
        endif
        itZ( numat(jt) ) = jt
      end do boucle_2

      boucle_igr: do igr = 1,n_atom_neq_min+1
        if( .not. Doping .and. igr == n_atom_neq_min+1 ) cycle
        do it = 1,ntype_conf
          do kgr = 1,nb_atom_conf(it)
            if( igr /= igra(kgr,it) ) cycle
            itype(igr) = it
            cycle boucle_igr
          end do
        end do
        if( Doping .and. igr == n_atom_neq_min+1 ) then
          if( itype_dop == 1 ) then
            itype(igr) = itype_dop
          else
            itype(igr) = itZ( itype_dop )
          endif
        else
          itype(igr) = itZ( itype(igr) )
        endif
      end do boucle_igr

      deallocate( itZ )

      boucle_ig: do igr = 1,n_atom_bulk
        do it = ntype_conf+1,jt
          if( Z_bulk(igr) == numat(it) ) cycle boucle_ig
        end do

        jt = jt + 1
        numat(jt) = Z_bulk(igr)
      end do boucle_ig

      if( itype_dop /= 1 ) itype_dop = ntype

    endif

    if( .not. Matper .or. n_atom_per == 0 ) Space_group = ' '

    iabsorig(:) = iabsm(:)

    if( Matper .and. Space_group /= ' ' ) then

! Recopy of surface and inter-layer atoms in film
      n = n_atom_int + n_atom_sur
      do i = n_atom_int + n_atom_sur,1,-1
        ia = n_atom_uc - n + i
        igr = n_atom_per_neq + i
        if( igr_dop == igr ) igr_dop = ia
        do j = 1,n_atom_coop
          if( igr_coop(j) == igr ) igr_coop(j) = ia
        end do
        posn(:,ia) = posn(:,igr)
        itype(ia) = itype(igr)
        if( Taux ) Taux_oc(ia) = Taux_oc(igr)
        if( Temp_B_iso ) Temp_coef(ia) = Temp_coef(igr)
        if( ngroup_pdb > 0 ) Kgroup(ia) = Kgroup(igr)

        Ang_base_loc_gr(:,ia) = Ang_base_loc_gr(:,igr )

        if( ia <= ngroup_nonsph ) then
          norbv(ia) = norbv(igr)
          if( norbv(ia) /= 0 ) then
            hybrid(:,:,ia) = hybrid(:,:,igr)
            pop_nonsph(:,ia) = pop_nonsph(:,igr)
          endif
        endif
        if( Atom_occ_hubb ) occ_hubb_e(:,:,:,ia) = occ_hubb_e(:,:,:,igr)

        do multi_run = 1,n_multi_run_e
          if( iabsorig(multi_run) == igr ) iabsm(multi_run) = ia
        end do
      end do

      iabsorig(:) = iabsm(:)

      nn = min( n_atom_per_neq, n_atom_uc )
      allocate( neq(nn) )
      allocate( pos(3,nn) )
      do igr = 1,nn
        pos(:,igr) = posn(:,igr)
      end do
      if( n_atom_uc > n_atom_per_neq ) then
        call spgroup(Fcif,Fcif_file,.true.,neq,n_atom_uc,nn,itype,pos,posn,Space_group)
      else
        neq(:) = 1
      endif

      if( Doping .and. igr_dop <= nn ) then
        do igr = 1,n_atom_uc - n
          if( sum( abs( posn(:,igr) - pos(:,igr_dop) ) ) > eps10 ) cycle
          igr_dop = igr
          exit
        end do
      endif

      do j = 1,n_atom_coop
        if( igr_coop(j) > nn ) cycle
        do igr = 1,n_atom_uc - n
          if( sum( abs( posn(:,igr) - pos(:,igr_coop(j)) ) ) > eps10 ) cycle
          igr_coop(j) = igr
          exit
        end do
      end do

      ia = n_atom_uc - n + 1
      do igr = nn,1,-1
        do i = 1,neq(igr)
          ia = ia - 1
          itype(ia) = itype(igr)
          if( Taux ) Taux_oc(ia) = Taux_oc(igr)
          if( Temp_B_iso ) Temp_coef(ia) = Temp_coef(igr)
          if( ngroup_pdb > 0 ) Kgroup(ia) = Kgroup(igr)

          Ang_base_loc_gr(:,ia) = Ang_base_loc_gr(:,igr )

          if( ia > ngroup_nonsph ) cycle

          norbv(ia) = norbv(igr)
          if( norbv(ia) /= 0 ) then
            hybrid(:,:,ia) = hybrid(:,:,igr)
            pop_nonsph(:,ia) = pop_nonsph(:,igr)
          endif
          if( Atom_occ_hubb ) occ_hubb_e(:,:,:,ia) = occ_hubb_e(:,:,:,igr)

        end do
        do multi_run = 1,n_multi_run_e
          if( iabsorig(multi_run) == igr ) iabsm(multi_run) = ia
        end do
      end do

      deallocate( neq, pos )

    endif

    if( Doping ) iabsm(1) = n_atom_uc + 1
    if( .not. Taux ) ATA = .FALSE.
    if( COOP ) Normaltau = .false. ! Normaltau is used in MSM. I do not remember well when it can be usefull...

    if( RIXS ) then
      Green_s = .false.
      Green_int = .false.
    endif

! With spinorbit, the renormalization is mandatory. I have not checked why.
    if( Green_s ) Green_bulk = .true.
    if( Green_s .or. Renorm .or. Spinorbite .or. .not. ( Optic  .or. No_renorm ) ) then
      Renorm = .true.
    else
      Renorm = .false.
    endif

    if( Pdb ) One_run = .false.

    if( numat_abs(1) /= 0 ) then
      do i = 1,n_Z_abs
        do igr = 1,ngroup
          if( Atom .or. Atom_conf .or. Flapw ) then
            if( abs( itype( igr ) ) > ntype ) then
              call write_error
              do ipr = ipr0,9,3
                write(ipr,174) numat_abs(i)
              end do
              stop
            endif
            if( numat( abs( itype( igr ) ) ) == numat_abs(i) ) exit
          else
            if( itype( igr ) == numat_abs(i) ) exit
          endif
        end do
        if( igr > ngroup ) then
          call write_error
          do ipr = ipr0,9,3
            write(ipr,175) numat_abs(i)
          end do
          stop
        endif
      end do
    else
      if( Doping ) then
        if( Atom .or. Atom_conf ) then
          numat_abs(1) = numat( abs( itype_dop ) )
        else
          numat_abs(1) = itype_dop
        endif
      else
        if( Atom .or. Atom_conf .or. Flapw ) then
          it = abs( itype( iabsm(1) ) )
          if( it > ntype ) then
            call write_error
            do ipr = ipr0,9,3
              write(ipr,176) it, ntype
            end do
            stop
          endif
          numat_abs(1) = numat( it )
        else
          numat_abs(1) = itype( iabsm(1) )
        endif
      endif
    endif

    if( Clementi ) then
      Wien_file(8) = 'clementi'
      nrm = nrato_dirac
      Fermi_auto = .false.
    endif

    if( Spinorbite ) COOP_z_along_bond = .false.  ! because the rotation does not work with spinorbit

    if( Atom_occ_hubb ) Fermi_auto = .false.
    if( Flapw .or. Extract .or. .not. Fermi_auto ) then
      SCF = .false.
      Fermi_auto = .false.
      nself = 0
      nself_bulk = 0
    end if
    if( .not. r_self_imp ) r_self = rsorte_s(1)
    if( SCF ) then
      Self_cons = .true.
      if( nself == 0 ) nself = 100
    elseif( ( Fermi_auto .or. ( Hubbard .and. .not. ( Atom_occ_hubb .or. Extract ) ) ) .and. nself == 0 ) then
      Self_cons = .true.
      p_self0 = 0._db
      nself = 1
      if( .not. ( r_self_imp .or. One_run ) ) r_self = min( rsorte_s(1), 3.5_db )
    endif

    if( No_SCF_bulk .and. SCF ) then
      SCF_bulk = .false.
      if( ( Fermi_auto .or. ( Hubbard .and. .not. ( Atom_occ_hubb .or. Extract ) ) ) .and. nself_bulk == 0 ) then
        Self_cons_bulk = .true.
        nself_bulk = 1
      else
        Self_cons_bulk = .false.
        nself_bulk = 0
      endif
      p_self0_bulk = 0._db
      R_self_bulk = min( rsorte_s(1), 3.5_db )
    else
      SCF_bulk = SCF
      Self_cons_bulk = Self_cons
      p_self0_bulk = p_self0
      nself_bulk = nself
      R_self_bulk = R_self
    endif

    if( Quadrupole ) then
      if( .not. no_dipquad ) E1E2e = .true.
      if( .not. no_e2e2 ) E2E2 = .true.
    endif
    if( Octupole .and. .not. no_e1e3 ) E1E3 = .true.
    if( Octupole .and. .not. no_e3e3 ) E3E3 = .true.
    if( Dipmag ) then
      E1M1 = .true.
      M1M1 = .true.
    endif
    if( E1E2e .or. E2E2 ) Quadrupole = .true.
    if( E1M1 .or. M1M1 ) Dipmag = .true.
    if( E1M2 .or. M2M2 ) Quadmag = .true.
    if( E1E3 .or. E3E3 ) Octupole = .true.

    if( Density_comp ) Harm_cubic = .false.

    if( .not. Pas_SCF_imp ) then
      if( Green_self ) then
        Pas_SCF = 0.1_db
      else
        Pas_SCF = 0.01_db
      endif
    endif

! Film preparation
    if( Film ) then
      if( .not. Bulk ) hkl_film = .true.
    else
      Bulk = .false.
      Cap_layer = .false.
    endif
    if( .not. Film .and. Matper ) V_helm = 0._db

    if( n_atom_per == 0 ) Film_shift(:) = 0._db
    if( n_atom_int == 0 ) Interface_shift(:) = 0._db
    if( n_atom_sur == 0 ) Surface_shift(:) = 0._db
    if( Film .and. n_atom_per == 0 .and.  n_atom_sur == 0 ) then
      n_atom_sur = n_atom_int
      n_atom_int = 0
      Surface_shift(:) = Interface_shift(:)
      Interface_shift(:) = 0._db
      Surface_shift_given = Interface_shift_given
      Interface_shift_given = .false.
      axyz_sur(:) = axyz_int(:)
      angxyz_sur(:) = angxyz_int(:)
    endif
    if( Film .and. n_atom_per == 0 ) then
      axyz(:) = axyz_sur(:)
      angxyz(:) = angxyz_sur(:)
    endif

    if( Film .and. n_atom_per /= 0 ) then
      if( Film_zero > eps10 ) then
        do jgr = 1,n_atom_per
          posn(3,jgr) = posn(3,jgr) - Film_zero
          if( posn(3,jgr) < - eps10 ) posn(3,jgr) = posn(3,jgr) + 1._db
        end do
      endif
      z_min = 10000._db
      do jgr = 1,n_atom_per
        z_min = min( z_min, posn(3,jgr) )
      end do
      do jgr = 1,n_atom_per
        posn(3,jgr) = posn(3,jgr) - z_min
      end do
    endif
    if( Film .and. n_atom_int /= 0 ) then
      z_min = 10000._db
      do jgr = n_atom_per+1, n_atom_per+n_atom_int
        z_min = min( z_min, posn(3,jgr) )
      end do
      do jgr = n_atom_per+1, n_atom_per+n_atom_int
        posn(3,jgr) = posn(3,jgr) - z_min
      end do
    endif
    if( Film .and. n_atom_sur /= 0 ) then
      z_min = 10000._db
      do jgr = n_atom_per+n_atom_int+1, n_atom_uc
        z_min = min( z_min, posn(3,jgr) )
      end do
      if( .not. Surface_shift_given ) Surface_shift(3) = z_min * axyz_sur(3)
      do jgr = n_atom_per+n_atom_int+1, n_atom_uc
        posn(3,jgr) = posn(3,jgr) - z_min
      end do
      if( Noncentre ) Centre(3) = Centre(3) - z_min
    endif
    if( Cap_layer ) then
      z_min = 10000._db
      do jgr = 1,n_atom_cap
        z_min = min( z_min, posn_cap(3,jgr) )
      end do
      do jgr = 1,n_atom_cap
        posn_cap(3,jgr) = posn_cap(3,jgr) - z_min
      end do
    endif
    if( Bulk ) then
      z_min = 10000._db
      do jgr = 1,n_atom_bulk
        z_min = min( z_min, posn_bulk(3,jgr) )
      end do
      do jgr = 1,n_atom_bulk
        posn_bulk(3,jgr) = posn_bulk(3,jgr) - z_min
      end do
    endif

    if( Operation_mode_used ) then
      Bormann = .false.
      Self_abs = .false.
      Full_self_abs = .false.
    endif

    if( Extract_ten ) then
      Density = .false.
      if( Quadrupole ) then
        open(1, file = nom_fich_Extract, status='old', iostat=istat)
        do i = 1,100000
          read(1,'(A)') mot
          if( mot(2:6) == 'Dipol'  ) then
            read(1,'(A)') mot
            if( mot(2:14) /= 'E1E1 includes' ) then
              Dip_rel = .false.
              backspace(1)
            endif
            read(1,'(A)') mot
            if( mot(2:6) /= 'Quadr' ) Quadrupole = .false.
            read(1,'(A)') mot
            if( mot(2:6) /= 'Octup' ) octupole = .false.
            exit
          endif
        end do
        Close(1)
     endif
    endif

    if( Extract ) then
      State_all = .false.
      State_all_out = .false.
      COOP = .false.
      call Extract_log(Core_resolved,Green_bulk,Green_s,nom_fich_extract,Optic,Renorm,Seuil)
    endif

    if( Spinorbite .and. non_relat == 0 ) Relativiste = .true.

    if( Optic .and. .not. Hedin .and. abs(alfpot) < eps6 ) Perdew = .true.

    if( .not. Flapw .and. abs(alfpot) < eps6 .and. .not. Perdew ) Hedin = .true.
    if( Hedin ) then
      alfpot = 0._db
    elseif( Flapw ) then
      alfpot = 0.33_db
    elseif( Perdew ) then
      alfpot = -1._db
    endif

    if( Kern_fac_default ) then
      if( Optic ) then
        Kern_fac = 1._db
      else
        Kern_fac = 2._db
      endif
    endif

    if( Optic ) then
      Seuil = 'Opt'
      Dip_rel = .false.
      nseuil = 0; Lseuil = 0; jseuil = 0; nbseuil = 1
    else
      select case( Seuil(1:1) )
        case('K')
          nseuil = 1
        case('L')
          nseuil = 2
        case('M')
          nseuil = 3
        case('N')
          nseuil = 4
        case('O')
          nseuil = 5
        case('P')
          nseuil = 6
      end select
      select case( Seuil(2:2) )
        case('1')
          Lseuil = 0; jseuil = 1
        case('2')
          Lseuil = 1; jseuil = 2
        case('3')
          Lseuil = 1; jseuil = 3
        case('4')
          Lseuil = 2; jseuil = 4
        case('5')
          Lseuil = 2; jseuil = 5
        case('6')
          Lseuil = 3; jseuil = 6
        case('7')
          Lseuil = 3; jseuil = 7
      end select
      nbseuil = len_trim( Seuil ) - 1
    endif

    if( Flapw .or. Lseuil > 0 .or. Optic ) Nonexc = .true.
    if( Exc_imp .and. .not. Flapw ) Nonexc = .false.
    if( Nonexc ) Symmol = .true.
    if( Self_nonexc_imp ) Self_nonexc = .true.
    if( SCF_exc_imp .and. .not. Optic ) Self_nonexc = .false.
    if( SCF_exc_imp .and. .not. Optic .and. .not. nonexc_imp ) Nonexc = .false.
    if( Nonexc ) Self_nonexc = .true.
    if( .not. Nonexc .and. .not. SCF ) Self_nonexc = .false.
    L = L_level_val(numat_abs(1))
! Below Self_nonexc would be more logical than ( Self_nonexc .or. .not. SCF ), but we do not mind
    if( .not. Nonexc .and. ( Self_nonexc .or. .not. SCF ) .and. ( L == 2 .or. L == 3 ) ) Scf_elecabs = .true.

    if( .not. ( Flapw .or. Atom .or. Atom_conf ) ) then
! In these case, up to here, itype was the atomic number
      jt = 0
      boucle_1: do igr = 1,ngroup - n_atom_bulk
        n = itype(igr)
        do it = 1,jt
          if( numat(it) /= n ) cycle
          itype(igr) = it
          if( Doping .and. igr == ngroup - n_atom_bulk ) itype_dop = jt
          cycle boucle_1
        end do
        jt = jt + 1
        if( Doping .and. igr == ngroup - n_atom_bulk ) itype_dop = jt
        numat(jt) = n
        itype(igr) = jt
      end do boucle_1
      nlat(1:jt) = 0
    endif

    Loop_bulk: do igr = 1,n_atom_bulk
      do it = 1,jt
        if( Z_bulk(igr) == numat(it) .and. nlat(it) == 0 ) then
          itype(ngroup - n_atom_bulk + igr) = it
          cycle Loop_bulk
        endif
      end do
      jt = jt + 1
      numat(jt) = Z_bulk(igr)
      itype( ngroup - n_atom_bulk + igr ) = jt
    end do Loop_bulk

 !   if( ntype_bulk > 0 ) then
 !     it = ntype - ntype_bulk + 1
 !     itype( ngroup - n_atom_bulk + 1 ) = it
 !     numat(it) = Z_bulk( 1 )
 !     do igr = 2,n_atom_bulk
 !       do jgr = 1,igr - 1
 !         if( Z_bulk(igr) == Z_bulk(jgr) ) then
 !           itype( ngroup - n_atom_bulk + igr ) = itype( ngroup - n_atom_bulk + jgr )
 !           exit
 !         endif
 !       end do
 !       if( jgr > igr - 1 ) then
 !         it = it + 1
 !         itype( ngroup - n_atom_bulk + igr ) = it
 !         numat(it) = Z_bulk(igr)
 !       endif
 !     end do
 !   endif

    do igr = 1,ngroup
      it = abs( itype(igr) )
      do l = 1,nlat(it)
        do ispin = 1,nspin
          popats(igr,l,ispin) = popval(it,l,ispin)
        end do
      end do
    end do

    do it = 1,ntype
      if( numat(it) == 0 ) then
        com(it) = ' Empty sphere'
        icom(it) = 5
      endif
    end do

    if( Extract_ten ) then
      if( no_core_resolved ) Core_resolved = .false.
    elseif( ( nspin == 1 .or. no_core_resolved ) .and. .not. Core_resolved_e ) then
      Core_resolved = .false.
    else
      Core_resolved = .true.
    endif
! When there are more than 1 edge, Gamma must be taken in the Chi_0 calculation
    if( nbseuil > 1 .or. ngamh > 1 ) then
      if( Gamma_max > eps10 ) Gamma_tddft = .true.
      if( Gamma_hole_imp ) then
        g1 = Gamma_hole(1)
        g2 = Gamma_hole(1)
        do i = 2,ngamh
          g1 = min( g1, Gamma_hole(i) )
          g2 = max( g2, Gamma_hole(i) )
        end do
        if( abs( g1 - g2 ) > eps10 ) Gamma_tddft = .true.
      else
        if( Lseuil == 1 .and. nseuil == 2 ) Gamma_tddft = .true.
      endif
    endif

    if( Bormann ) then
      Dafs = .true.
      Angpoldafs(3,:) = Ang_borm
      do ipl = 1,npldafs/4
       hkl_dafs(:,ipl) = real( hkl_borm(:), db )
      end do
      do ipl = npldafs/4 + 1,3*npldafs/4
        hkl_dafs(:,ipl) = 0._db
      end do
      do ipl = 3*npldafs/4 + 1,npldafs
        hkl_dafs(:,ipl) = - hkl_borm(:)
      end do
    endif

! To avoid singularities
    if( Film .and. Bulk ) then
      Avoid = .false.
      do igr = 1,n_atom_bulk
        do i = 1,n_Z_abs
          if( Z_bulk(igr) /= numat_abs(i) ) cycle
          Avoid = .true.
          exit
        end do
      end do
      Avoid = .false.
      if( Avoid ) then
        do ipl = 1,npldafs
          if( abs( nint( hkl_dafs(3,ipl) ) - hkl_dafs(3,ipl) ) < 0.004_db - eps10 ) then
            if( hkl_dafs(3,ipl) > nint( hkl_dafs(3,ipl) ) - eps10 ) then
              hkl_dafs(3,ipl) =  nint( hkl_dafs(3,ipl) ) + 0.004_db
            else
              hkl_dafs(3,ipl) =  nint( hkl_dafs(3,ipl) ) - 0.004_db
            endif
          endif
        end do
      endif
    endif

! Check of the inputs
    if( Seuil /= 'K1' .and. Seuil /= 'L1' .and. Seuil /= 'L2' .and. &
        Seuil /= 'L3' .and. Seuil /= 'M1' .and. Seuil /= 'M2' .and. &
        Seuil /= 'M3' .and. Seuil /= 'M4' .and. Seuil /= 'M5' .and. &
        Seuil /= 'N1' .and. Seuil /= 'N2' .and. Seuil /= 'N3' .and. &
        Seuil /= 'N4' .and. Seuil /= 'N5' .and. Seuil /= 'N6' .and. &
        Seuil /= 'N7' .and. Seuil /= 'O1' .and. Seuil /= 'O2' .and. &
        Seuil /= 'O3' .and. Seuil /= 'O4' .and. Seuil /= 'O5' .and. &
        Seuil /= 'P1' .and. Seuil /= 'P2' .and. Seuil /= 'P3' .and. &
        Seuil /= 'L23' .and. Seuil /= 'M23' .and. Seuil /= 'M45' .and. &
        Seuil /= 'N23' .and. Seuil /= 'N45' .and. Seuil /= 'N67' .and. &
        Seuil /= 'O23' .and. Seuil /= 'O45' .and. Seuil /= 'P23' .and. Seuil /= 'Opt' ) then
      call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,180) Seuil
      end do
      stop
    endif

    if( Doping .and. Absorber ) then
      Error_message = ' Absorber keyword not authorized when Doping keyword is used !'
      call write_error_message(Error_message,6,istop)
      stop
    endif

    if( One_run .and. Absorber ) then
      Error_message = ' Absorber keyword not authorized when One_run keyword is used !'
      call write_error_message(Error_message,6,istop)
      stop
    endif

    if( Atom ) then
      istop = 0
      do igr = 1,ngroup
        it = itype(igr)
        if( it <= ntype ) cycle
        if( istop == 0 ) then
          call write_error
          do ipr = ipr0,9,3
            write(ipr,100)
          end do
        endif
        do ipr = ipr0,9,3
          if( igr > ngroup ) then
            write(ipr,190) it, ntype
          else
            write(ipr,195) igr, it, ntype
          endif
        end do
        istop = 1
      end do
      if( istop == 1 ) then
        do ipr = ipr0,9,3
          write(ipr,200)
        end do
        stop
      endif
    endif

    if( Atom .or. Atom_conf ) then
      do it = 0,ntype
        do l = 1,nlat(it)
          q = sum( popval(it,l,:) )
          n = 2 * ( 2 * lvval(it,l) + 1 )
          if( q <=  n + eps10 ) cycle
          if( istop == 0 ) call write_error
          do ipr = ipr0,9,3
            write(ipr,100)
            write(ipr,210) it, l, lvval(it,l), q, n
          end do
          istop = 1
        end do
      end do
      if( istop == 1 ) stop
    endif

    if( Atom .and. Atom_conf ) then
      Error_message = &
      ' "Atom" .and. "Atom_conf" keywords cannot be used together !'
      call write_error_message(Error_message,ipr0,istop)
    endif

    if( .not. ( E1E1 .or. E1E2e .or. E1E3 .or. E2E2 .or. E1M1 .or. E1M2 .or. M1M1 .or. M1M2 .or. M2M2 ) ) then
      Error_message = &
      ' When keyword keyword "No_dipole" is used, at least one other transition channel must be selected (E1E2, E2E2, E1M1...)'
      call write_error_message(Error_message,ipr0,istop)
    endif

    if( ngroup == 0 ) then
      Error_message = '  There is no atom in the calculation ! A mandatary keywords as Molecule, Film or Crystal is missing'
      call write_error_message(Error_message,ipr0,istop)
    endif
    if( ntype == 0 .and. ngroup /= 0 ) then
      Error_message = '  There is no chemical specie specified in your calculation !'
      call write_error_message(Error_message,ipr0,istop)
    endif
    if( iord /= 2 .and. iord /= 4 ) then
      if( istop == 0 ) call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,220) iord
      end do
      istop = 1
    endif
    do igr = 1,ngroup
      if( abs(itype(igr)) <= ntype ) cycle
      if( istop == 0 ) call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,230) igr, itype(igr), ntype
      end do
      istop = 1
    end do
    do ipl = 1,nple
      pp = sum( polar(1:3,ipl)**2 )
      q = sum( veconde(1:3,ipl)**2 )
      if( pp < eps10 .and. q < eps10 ) then
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,240)
        end do
        istop = 1
      endif
    end do
    if( normrmt < 0 ) then
      if( istop == 0 ) call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,250) -normrmt, ntype
      end do
      istop = 1
    endif
    if( Full_self_abs ) then
      do ipl = 1,npldafs
        jpl = mod(ipl-1,4) + 1
        select case(jpl)
          case(1)
            if( isigpi(1,ipl) == 1 .and. isigpi(2,ipl) == 1 ) cycle
          case(2)
            if( isigpi(1,ipl) == 1 .and. isigpi(2,ipl) == 2 ) cycle
          case(3)
            if( isigpi(1,ipl) == 2 .and. isigpi(2,ipl) == 1 ) cycle
          case(4)
            if( isigpi(1,ipl) == 2 .and. isigpi(2,ipl) == 2 ) cycle
        end select
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,260) ipl, nint( hkl_dafs(:,ipl) ), isigpi(:,ipl)
        end do
        istop = 1
        exit
      end do
    endif

    if( npldafs > 0 ) then

      if( Film .and. Bulk ) then

        angxyz_bulk(:) = angxyz_bulk(:) * radian
        if( n_atom_per > 0 ) then
          Film_shift(4) = Film_shift(4) * radian
          call mult_film_bulk(angxyz_bulk,axyz,axyz_bulk,Film_shift,Mult_bulk,Mult_film)
          Film_shift(4) = Film_shift(4) / radian
          Diagonal = abs( angxyz_bulk(3) - 90 ) < 1._db .and. &
                  ( abs( Film_shift(4) - 45 ) <  1._db .or. abs( Film_shift(4) + 45 ) < 1._db )
        elseif( n_atom_sur > 0 ) then
          Surface_shift(4) = Surface_shift(4) * radian
          call mult_film_bulk(angxyz_bulk,axyz_sur,axyz_bulk,Surface_shift,Mult_bulk,Mult_film)
          Surface_shift(4) = Surface_shift(4) / radian
          Diagonal = abs( angxyz_bulk(3) - 90 ) < 1._db .and. &
                  ( abs( Surface_shift(4) - 45 ) <  1._db .or. abs( Surface_shift(4) + 45 ) < 1._db )
        endif
        angxyz_bulk(:) = angxyz_bulk(:) / radian

        i = 0
        do ipl = 1,npldafs
          if( Diagonal ) then
            if( hkl_film ) then
              if( sum( abs( nint( Mult_film(1:2)*hkl_dafs(1:2,ipl) ) - Mult_film(1:2)*hkl_dafs(1:2,ipl) ) ) < eps10 ) cycle
            else
    ! sum of index must be integer
              if( abs( nint( sum( Mult_film(1:2)*hkl_dafs(1:2,ipl) ) ) - sum( Mult_film(1:2)*hkl_dafs(1:2,ipl) ) ) < eps10 ) cycle
            endif
          else
            if( hkl_film ) then
              if( sum( abs( nint( Mult_film(1:2)*hkl_dafs(1:2,ipl) ) - Mult_film(1:2)*hkl_dafs(1:2,ipl) ) ) < eps10 ) cycle
            else
              if( sum( abs( nint( Mult_bulk(1:2)*hkl_dafs(1:2,ipl) ) - Mult_bulk(1:2)*hkl_dafs(1:2,ipl) ) ) < eps10 ) cycle
            endif
          endif
          if( istop == 0 .and. i == 0 ) call write_error
          do ipr = ipr0,9,3
            if( istop == 0 .and. i == 0 ) write(ipr,100)
            if( hkl_film .and. i == 0 ) then
              write(ipr,265) 'film', Mult_film(1:2)
            elseif( i == 0 ) then
              if( Diagonal ) then
                write(ipr,263) 'bulk', Mult_bulk(1:2)
              else
                write(ipr,265) 'bulk', Mult_bulk(1:2)
              endif
            endif
            write(ipr,270) ipl, hkl_dafs(:,ipl)
          end do
          i = 1
          istop = 1
        end do

      else

        if( Film ) then
          n = 2
        else
          n = 3
        endif

        i = 0
        do ipl = 1,npldafs
          if( sum( abs( nint(hkl_dafs(1:n,ipl)) - hkl_dafs(1:n,ipl) ) ) < eps10 ) cycle
          if( istop == 0 .and. i == 0 ) call write_error
          do ipr = ipr0,9,3
            if( istop == 0 .and. i == 0 ) write(ipr,100)
            if( .not. Film .and. i == 0 ) then
              write(ipr,'(/A)') ' Reflexion indexes must be integer !'
            elseif( i == 0 ) then
              write(ipr,'(/A)') ' For a film without bulk, the 2 first reflexion indexes must be integer !'
            endif
            write(ipr,280) ipl, hkl_dafs(:,ipl)
          end do
          i = 1
          istop = 1
        end do

      endif

    endif

    if( n_radius > 1 ) then
      do i = 2,n_radius
        if( Rsorte_s(i-1) > Rsorte_s(i) - 0.0001_db ) cycle
        Error_message = ' Under keyword "Radius", the radius values must be decreasing !'
        call write_error_message(Error_message,ipr0,istop)
      end do
      do i = 2,n_radius
        if( E_radius(i-1) < E_radius(i) - 0.0001_db ) cycle
        Error_message = ' Under keyword "Radius", the energy values must be increasing !'
        call write_error_message(Error_message,ipr0,istop)
      end do
    endif

    if( n_adimp > 1 ) then
      do i = 2,n_adimp
        if( E_adimp(i-1) < E_adimp(i) - 0.0001_db ) cycle
        Error_message = ' Under keyword "Adimp", the energy values must be increasing !'
        call write_error_message(Error_message,ipr0,istop)
      end do
    endif

    if( One_run .and. n_range > 1 ) then
      Error_message = ' It is not possible to have "One_run" keyword and several values of adimp or radius !'
      call write_error_message(Error_message,ipr0,istop)
    endif

    if( istop == 1 ) stop
! end of inputs checking

! Normalisation of hybridation parameters
    if( Atom_nonsph ) then
      do igr = 1,ngroup
        if( norbv(igr) == 0 ) cycle
        do io = 1,norbv(igr)
          rn = sqrt( sum( abs(hybrid(io,:,igr))**2 ) )
          hybrid(io,:,igr) = hybrid(io,:,igr) / rn
        end do
      end do
    endif

    if( Tddft .or. Optic ) Eimagent(:) = 0._db
    if( Spinorbite ) Ylm_comp_inp = .true.

! Default energy grid
    if( lin_gam == - 1 ) then
      lin_gam = 0
      Egamme(1) = -2.0_db; Egamme(2) =  0.5_db
      Egamme(3) = Egamme(1) + ( nenerg_s - 1 ) * Egamme(2)
    endif

    if( Green_s ) then
      FDM_comp = .false.
    elseif( .not. FDM_comp .and. neimagent > 0 ) then
      Eimagent(:) = 0._db
    endif

! Reduction of the number of atom type
    ntype_mod = ntype
    if( .not. Flapw ) then
      Suppress(:) = .false.
      do it = 1,ntype
        boucle_jt: do jt = it + 1,ntype
          if( Suppress(jt) ) cycle
          if( numat(it) /= numat(jt) ) cycle
          if( nlat(it) /= nlat(jt) ) cycle
          do l = 1,nlat(it)
            if( sum ( abs( popval(it,l,:) - popval(jt,l,:) ) ) > eps10 ) cycle boucle_jt
          end do
          ntype_mod = ntype_mod - 1
          Suppress(jt) = .true.
          where( itype == jt ) itype = it
        end do boucle_jt
      end do
      do it = 1,ntype
        if( .not. Suppress(it) ) cycle
        n = 0
        do jt = it,ntype-1
          if( Suppress( jt+1 ) ) then
            n = n + 1
            cycle
          endif
          kt = jt - n
          where( itype == jt+1 ) itype = itype - n - 1
          Ang_base_loc(:,kt) = Ang_base_loc(:,jt+1)
          com(kt) = com(jt+1)
          Hubb(kt) = Hubb(jt+1)
          icom(kt) = icom(jt+1)
          lvval(kt,:) = lvval(jt+1,:)
          nlat(kt) = nlat(jt+1)
          nrato(kt) = nrato(jt+1)
          nrato_lapw(kt) = nrato_lapw(jt+1)
          numat(kt) = numat(jt+1)
          nvval(kt,:) = nvval(jt+1,:)
          popval(kt,:,:) = popval(jt+1,:,:)
          r0_lapw(kt) = r0_lapw(jt+1)
          rchimp(kt) = rchimp(jt+1)
          rlapw(kt) = rlapw(jt+1)
          rmt(kt) = rmt(jt+1)
          rmtimp(kt) = rmtimp(jt+1)
          V_hubbard(kt) = V_hubbard(jt+1)
        end do
      end do
    endif

    if( n_rchimp_Z > 0 ) then
     do it = 1,ntype
       do i = 1,n_rchimp_Z
         if( numat(it) /= Z_rchimp(i) ) cycle
         rchimp(it) = Rchimp_Z(i)
         exit
       end do
     end do
   endif

   if( n_rmtg_Z > 0 ) then
     do it = 1,ntype
       do i = 1,n_rmtg_Z
         if( numat(it) /= Z_rmtg(i) ) cycle
         rmtimp(it) = Rmtg_Z(i)
         exit
       end do
     end do
   endif

   if( n_hubbard_Z > 0 ) then
     do it = 1,ntype
       do i = 1,n_hubbard_Z
         if( numat(it) /= Z_hubbard(i) ) cycle
         V_hubbard(it) = Hubbard_Z(i)
         Hubb(it) = abs( V_hubbard(it) ) > eps10
         exit
       end do
     end do
   endif

! Writing of indata:
    if( icheck(1) > 0 ) then

      if( Extract ) then
        Length = len_trim(nom_fich_Extract)
        if( Extract_ten ) then
          write(3,300) nom_fich_Extract(1:Length)
          write(6,300) nom_fich_Extract(1:Length)
        else
          write(3,305) nom_fich_Extract(1:Length)
          write(6,305) nom_fich_Extract(1:Length)
        endif
        open(1, file = nom_fich_Extract, status='old',iostat=istat)
        if( istat /= 0 ) call write_open_error(nom_fich_Extract,istat,1)
        do i = 1,100000
          read(1,'(A)') mot
          write(3,'(A)') mot
          if( mot(2:6) == 'Dipol'  ) exit
        end do
      else
        if( RIXS ) then
          write(3,'(/A)') ' RIXS calculation'
          do i = 1,n_q_rixs
            write(3,'(a6,3f7.3,a9,f8.3)') '   q =', q_rixs(1:3,i), ', Angle =', q_rixs(4,i)
          end do
          if( RIXS_fast ) write(3,'(A)') '   Fast calculation mode'
          if( Moment_conservation )then
            write(3,'(A)') '   Conservation of the moment'
          else
            write(3,'(A)') '   No moment conservation'
          endif
        endif
        do i = 1,n_Z_abs
          l2 = len_trim(Seuil)
          mot = ' '
          mot = ' Threshold:'
          j = 12
          mot13 = Chemical_Name(numat_abs(i))
          l1 = len_trim(mot13)
          mot(j+1:j+l1) = mot13(1:l1)
          j = j + l1 + 1
          mot(j+1:j+l2) = Seuil(1:l2)
          j = j + l2 + 1
          mot(j+1:j+4) = 'edge'
          j = j + 4
          write(3,'(/A)') mot
          write(6,'(A)') mot(1:j)
        end do
        if( n_radius == 1 ) then
          write(3,310) Rsorte_s(1)
        else
          write(3,315) Rsorte_s(1), ( E_radius(i-1), Rsorte_s(i), i = 2,n_radius )
        endif
        if( .not. Green_s .and. overad ) write(3,320) roverad
        if( lin_gam == 1 ) write(3,340)
        if( abs( Egamme(1) ) < 999.999_db .and. abs( Egamme(ngamme) ) < 999.999_db ) then
          write(3,350) Egamme(1:ngamme)
        else
          write(3,360) Egamme(1:ngamme)
        endif
        write(3,'(A)') ' Dipole component'
      endif
      if( Dip_rel ) write(3,'(A)') ' E1E1 includes the spin-orbit relativistic transition channel'
      if( Quadrupole ) write(3,'(A)') ' Quadrupole component'
      if( Octupole ) write(3,'(A)') ' Octupole component'
      if( Dipmag ) write(3,'(A)') ' Magnetic dipole component'
      if( .not. E1E1 ) write(3,'(A)') ' ... but E1E1 component neglected in the output'
      if( no_E2E2 ) write(3,'(A)') ' ... but E2E2 component neglected in the output'
      if( no_E1E3 ) write(3,'(A)') ' ... but E1E3 component neglected in the output'
      if( no_E3E3 ) write(3,'(A)') ' ... but E3E3 component neglected in the output'
      if( no_dipquad ) write(3,'(A)') ' ... but E1E2 component neglected in the output'
      if( tddft ) then
        if( rpalf ) then
          write(3,'(A)') ' TDDFT calculation using the RPA LF approximation'
        else
          write(3,'(A)') ' TDDFT calculation'
        endif
        if( Gamma_tddft ) then
          write(3,'(A)') '   Broadening in the Chi_0 calculation'
        else
          write(3,'(A)') '   No broadening in the Chi_0 calculation'
        endif
        write(3,'(a13,f5.2)') '   Kern_fac =', Kern_fac
        if( Optic .and. Kern_fast ) write(3,'(A)') '   Single calculation of Kernel with wave functions at Fermi energy'
      endif
      if( No_DFT ) write(3,'(A)') '   No DFT output'
      if( ( Optic .or. Tddft ) .and. lmax_tddft_inp >= 0 ) write(3,325) lmax_tddft_inp
      if( Core_resolved ) write(3,'(A)') ' Core resolved in outputs'

      if( Extract_ten ) then
        do i = 1,100000
          read(1,'(A)') mot
          if( mot(2:9) /= 'Molecule' .and. mot(2:8) /= 'Crystal' ) cycle
          write(3,'(A)') mot
          exit
        end do
      else
        if( Relativiste ) then
          write(3,'(A)') ' Relativistic calculation'
        else
          write(3,'(A)') ' Non-relativistic calculation'
        endif
        if( Magnetic .and. Spinorbite ) then
          write(3,'(A)') ' Magnetic calculation with spin-orbit interaction'
        elseif( Spinorbite ) then
          write(3,'(A)') ' Non-Magnetic calculation with spin-orbit interaction'
        elseif( Magnetic ) then
          write(3,'(A)') ' Magnetic calculation without spin-orbit interaction'
        else
          write(3,'(A)') ' Non-magnetic calculation'
        endif
        if( Z_nospinorbite /= 0 ) write(3,400) Z_nospinorbite
        if( lmoins1 ) write(3,'(A)') ' Approximation l-1'
        if( lplus1 ) write(3,'(A)') ' Approximation l+1'
        if( Basereel ) then
          write(3,'(A)') ' Real bases'
        else
          write(3,'(A)') ' Complex bases'
        endif
        if( Green_s ) then
          write(3,'(A)') ' Multiple scattering calculation (Green)'
          if( Green_int ) write(3,'(A)') '    Full Green mode'
          if( supermuf ) write(3,'(A)') '    Continuous potential (Supermuf)'
          if( nchemin /= - 1 ) write(3,410) nchemin
          if( normaltau ) write(3,'(A)') '    Tau normalization'
          select case(normrmt)
            case(1)
              write(3,'(A)') '   Optimized type muffin-tin radius'
            case(2)
              write(3,'(A)') '   Norman type muffin-tin radius'
            case(3)
              ! suppressed
            case(4)
              write(3,420) rmtimp(1:ntype_mod)
            case(5)
              write(3,'(A)') '   Potential imposed type muffin-tin radius'
          end select
          if( abs(overlap) > eps6 ) write(3,430) overlap
          if( lmaxfree ) then
            write(3,'(A)') '   No limitation on the maximum value of l'
          else
            write(3,'(A)') '   Limitation on the maximum value of l'
          endif
        else
          write(3,'(A)') ' Finite difference method calculation'
          if( Bulk .and. Green_bulk ) write(3,'(3x,A)') ' but in the bulk where multiple scattering theory (Green) is used'
          if( n_adimp == 1 ) then
            write(3,440) iord, adimp(1)
          else
            write(3,445) iord, adimp(1), ( E_adimp(i-1), adimp(i), i = 2,n_adimp )
          endif
          write(3,450) lmaxso0
          if( lmaxso_max < 1000000 ) write(3,455) lmaxso_max
          if( Muffintin ) write(3,'(A)') '   Muffin-tin potential'
          if( Rydberg ) write(3,460) R_rydb
          if( .not. Eneg_i ) write(3,480) Eclie, Eclie_out
          if( V_intmax < 100000._db ) write(3,490) V_intmax
        endif
      endif
      if( Renorm ) then
        write(3,'(A)') ' Normalization of the atomic radial wave functions to one state per Rydberg'
      else
        write(3,'(A)') ' No normalization of the atomic radial wave functions'
      endif
      if( Atom_nonsph ) write(3,'(/A)') ' Calculation with non spherical orbitals'
      if( Temp > eps10 ) write(3,500) Temp
      if( Ray_max_dirac > eps6 ) write(3,502) Ray_max_dirac

      if( nple > 0 ) then
        if( sum( abs( pdpolar(:,:) ) ) > eps10 ) then
          write(3,'(/A)') ' XANES :    Polarization             Wave vector     Weight_dip Weight_quad'
          do ipl = 1,nple
            write(3,505) polar(1:3,ipl), veconde(1:3,ipl), pdpolar(ipl,:)
          end do
        else
          write(3,'(/A)') ' XANES :    Polarization             Wave vector'
          do ipl = 1,nple - n_mat_polar
            write(3,505) polar(1:3,ipl), veconde(1:3,ipl)
          end do
        endif
        do ipl = nple - n_mat_polar + 1, nple
          write(3,510) polar(1:3,ipl), veconde(1:3,ipl)
        end do
      endif
      if( nq_nrixs >= 1 ) then
        write(3,'(/A)') ' NRIXS (X-ray Raman) calculation'
        write(3,'(a15,i2)') '   lmax_nrixs =' ,lmax_nrixs
        if( Monocrystal ) then
          do i = 1,nq_nrixs
            write(3,'(a6,f10.5,2x,3f10.5)') '   q =' , q_nrixs(4,i), q_nrixs(1:3,i)
          end do
        else
          write(3,'(a6,1000f10.5)') '   q =' , q_nrixs(4,:)
        endif
      endif

      if( Dafs_bio ) then

        write(3,'(1x,A)') 'DAFS reflections selected using the experimental files:'
        do i = 1,n_file_dafs_exp
          write(3,'(3x,A)') File_dafs_exp(i)
          write(3,'(5x,a12,f10.3)') 'with angle =', Angle_dafs_exp(i)
        end do

        write(3,'(A)') ' Experimental orientation matrix:'
        do i = 1,3
          write(3,'(3x,3f12.6)') Mat_or(i,:)
        end do

        Write(3,'(a23,i6)') ' Number of reflexions =', npldafs

      elseif( Bormann ) then

        write(3,511) nint( hkl_dafs(:,1) ), Angpoldafs(3,1)

      elseif( Operation_mode_used ) then

        write(3,'(/A)') ' DAFS :     h      k     l1   step     lf   Operation mode   Angle_1  Angle_2  Angle_3'
        jpl = 0
        do ipl = 1,npldafs_2d
          jpl = jpl + 2*npl_2d(ipl)
          if( npl_2d(ipl) > 1 ) then
            r = ( hkl_dafs(3,jpl) - hkl_dafs(3,jpl-npl_2d(ipl)+1) ) / ( npl_2d(ipl) - 1 )
          else
            r = 0._db
          endif
          write(3,512) hkl_dafs(1:3,jpl-2*npl_2d(ipl)+1), r, hkl_dafs(3,jpl), Operation_mode(:,ipl), Angle_mode(:,ipl)
        end do
        if( sum( abs( Mat_UB(:,:) ) ) > eps10 ) then
          write(3,'(/a14,3f12.7,2(/14x,3f12.7))') '   UB matrix =', ( Mat_UB(i,:), i = 1,3 )
        else
          write(3,'(/a36,f8.3)') '   Azimuth in the orientation matrix =', phi_0
        endif
        write(3,'(a24,3i4)')  '   Reference axis hkl =', hkl_ref(:)

      else

        do ipl = 1,npldafs
          motpol = ' '
          l = 0
          do i = 1,2
            if( isigpi(i,ipl) == 1 ) then
              motpol(l+1:l+5) = 'sigma'
            elseif( isigpi(i,ipl) == 2 ) then
              motpol(l+1:l+2) = 'pi'
            elseif( isigpi(i,ipl) == 3 ) then
              motpol(l+1:l+5) = 'right'
            elseif( isigpi(i,ipl) == 4 ) then
              motpol(l+1:l+4) = 'left'
            elseif( isigpi(i,ipl) == 5 .or. isigpi(i,ipl) == 10 ) then
              motpol(l+1:l+5) = 'recti'
            endif
            l = l + 6
          end do

          if( angpoldafs(3,ipl) > 9999._db ) then
            if( ipl == 1 ) write(3,'(/A)') ' DAFS :       Polarization                Wave vector'
            write(3,515) hkl_dafs(:,ipl)
            write(3,513) poldafse(1:3,ipl), vecdafsem(1:3,ipl), ' incoming'
            write(3,513) poldafss(1:3,ipl), vecdafssm(1:3,ipl), ' outgoing'
          else
            if( Film ) then
              if( ipl == 1 ) write(3,'(/A)') ' DAFS : (    h,     k,     l)  Polarization   Angle_i   Angle_o  Azimuth'
              if( npldafs > 20 .and. mod(ipl,10) /= 1 ) cycle
              if( angpoldafs(3,ipl) < -9999._db ) then
                write(3,515) hkl_dafs(:,ipl), motpol, angpoldafs(1:2,ipl)
              elseif( isigpi(1,ipl) == 10 ) then
                write(3,516) hkl_dafs(:,ipl), motpol, angpoldafs(2:3,ipl)
              elseif( isigpi(2,ipl) == 10 ) then
                write(3,517) hkl_dafs(:,ipl), motpol, angpoldafs(1:3:2,ipl)
              else
                write(3,518) hkl_dafs(:,ipl), motpol, angpoldafs(1:3,ipl)
              endif
            else
              if( ipl == 1 ) write(3,'(/A)') ' DAFS : (h, k, l)  Polarization   Angle_i   Angle_o  Azimuth'
              if( angpoldafs(3,ipl) < -9999._db ) then
                write(3,519) nint( hkl_dafs(:,ipl) ), motpol, angpoldafs(1:2,ipl)
              elseif( isigpi(1,ipl) == 10 ) then
                write(3,520) nint( hkl_dafs(:,ipl) ), motpol, angpoldafs(2:3,ipl)
              elseif( isigpi(2,ipl) == 10 ) then
                write(3,521) nint( hkl_dafs(:,ipl) ), motpol, angpoldafs(1:3:2,ipl)
              else
                write(3,522) nint( hkl_dafs(:,ipl) ), motpol, angpoldafs(1:3,ipl)
              endif
            endif
          endif

        end do

      endif

      if( Self_abs ) write(3,'(/A)') ' XANES for the RXS polarizations'
      if( Full_self_abs ) write(3,'(/A)') ' Self absorption taken into account with birefringence effect'

      if( Extract ) then

        Close(1)

      else

        write(3,*)
        write(3,542) ngroup, ntype_mod

        if( nlatm > 0 ) then
          write(3,'(A)') '   Typ  Z   n  l  Popul.'
          do it = 1,ntype_mod
            if( nspin == 1 ) then
              write(3,543) it, numat(it), ( nvval(it,l), lvval(it,l), popval(it,l,:), l = 1,nlat(it) )
            else
              write(3,544) it, numat(it), ( nvval(it,l), lvval(it,l), popval(it,l,:), l = 1,nlat(it) )
            endif
          end do
        endif
        if( norbdil > 0 ) write(3,548)
        do io = 1,norbdil
          write(3,549) itdil(io), ldil(io), cdil(io)
        end do
        if( Nonexc ) write(3,'(A)') '   Non excited absorbing atom'
        if( .not. PointGroup_Auto ) write(3,552) PointGroup
        if( Symmol .and. .not. Nonexc ) write(3,'(A)') '   Absorbing atom taken as not excited for the symmetry calculation'

      endif

      write(3,*)
      if( Film ) then
        write(3,'(a55,4i5)') ' Film,  n_atom_uc, n_atom_int, n_atom_sur, n_atom_per =', n_atom_uc, n_atom_int, &
                             n_atom_sur, n_atom_per
      elseif( Matper ) then
        write(3,'(a22,i5)') ' Crystal,  n_atom_uc =', n_atom_uc
      else
        write(3,'(a22,i5)') ' Molecule, n_atom_uc =', n_atom_uc
      endif
      if( Space_group /= ' ' .and. Matper ) write(3,554) Space_group
      if( Pdb .and. Taux .and. Atom_B_iso ) then
        mot = 'Kgroup  Occupancy    B_iso'
      elseif( Pdb .and. Taux .and. Atom_U_iso ) then
        mot = 'Kgroup  Occupancy    U_iso'
      elseif( Pdb .and. Taux ) then
        mot = 'Kgroup  Occupancy'
      elseif( Pdb .and. Temp_B_iso ) then
        mot = 'Kgroup    B_iso'
      elseif( Pdb ) then
        mot = 'Kgroup'
      elseif( Atom_nonsph ) then
        if( Taux .and. Atom_B_iso ) then
          mot = ' Occupancy    B_iso    norbv   popats'
        elseif( Taux .and. Atom_U_iso ) then
          mot = ' Occupancy    U_iso    norbv   popats'
        elseif( Taux ) then
          mot = ' Occupancy  norbv   popats'
        elseif( Atom_B_iso ) then
          mot = '   B_iso   norbv   popats'
        elseif( Atom_U_iso ) then
          mot = '   U_iso   norbv   popats'
        else
          mot = '  norbv  popats'
        endif
      else
        if( Taux .and. Atom_B_iso ) then
          mot = ' Occupancy   B_iso'
        elseif( Taux .and. Atom_U_iso ) then
          mot = ' Occupancy   U_iso'
        elseif( Taux ) then
          mot = ' Occupancy'
        elseif( Atom_B_iso ) then
          mot = '   B_iso'
        elseif( Atom_U_iso ) then
          mot = '   U_iso'
        else
          mot = ' '
        endif
      endif
      do jgr = 1,ngroup - n_atom_bulk
        if( Doping .and. jgr == ngroup - n_atom_bulk ) then
          write(3,'(/A)') '   Doping element :'
          igr = igr_dop
        else
          igr = jgr
        endif
        if( jgr == 1 .and. n_atom_per > 0 ) then
          write(3,555) axyz(1:3)
          r = sum( abs( real( nint( 1000_db * angxyz(:) ), db ) / 1000_db  - angxyz(:) ) )
          if( r < eps10 ) then
            write(3,556) angxyz(1:3)
          else
            write(3,557) angxyz(1:3)
          endif
        endif
        if( jgr == n_atom_per + 1 .and. n_atom_int > 0 ) then
          write(3,'(/A)') '   Interface elements :'
          write(3,555) axyz_int(1:3)
          r = sum( abs( real( nint( 1000_db * angxyz_int(:) ), db ) / 1000_db  - angxyz_int(:) ) )
          if( r < eps10 ) then
            write(3,556) angxyz_int(1:3)
          else
            write(3,557) angxyz_int(1:3)
          endif
        elseif( jgr == n_atom_per + n_atom_int + 1 .and. n_atom_sur > 0 ) then
          write(3,'(/A)') '   Surface elements :'
          write(3,555) axyz_sur(1:3)
          r = sum( abs( real( nint( 1000_db * angxyz_sur(:) ), db ) / 1000_db  - angxyz_sur(:) ) )
          if( r < eps10 ) then
            write(3,556) angxyz_sur(1:3)
          else
            write(3,557) angxyz_sur(1:3)
          endif
        endif
        if( jgr == 1 ) write(3,560) mot

        it = itype(jgr)
        Z = numat( abs(it) )
        if( Pdb .and. Taux .and. Temp_B_iso ) then
          write(3,565) Z, posn(:,igr), it, Kgroup(jgr), Taux_oc(jgr), Temp_coef(jgr)
        elseif( Pdb .and. Taux ) then
          write(3,565) Z, posn(:,igr), it, Kgroup(jgr), Taux_oc(jgr)
        elseif( Pdb .and. Temp_B_iso ) then
          write(3,570) Z, posn(:,igr), it, Kgroup(jgr), Temp_coef(jgr)
        elseif( Pdb ) then
          write(3,570) Z, posn(:,igr), it, Kgroup(jgr)
        elseif( norbv( min(jgr,ngroup_nonsph) ) == 0 ) then
          if( Taux .and. Temp_B_iso ) then
            write(3,571) Z, posn(:,igr), it, Taux_oc(jgr), Temp_coef(jgr)
          elseif( Taux ) then
            write(3,571) Z, posn(:,igr), it, Taux_oc(jgr)
          elseif( Temp_B_iso ) then
            write(3,572) Z, posn(:,igr), it, Temp_coef(jgr)
          else
            write(3,570) Z, posn(:,igr), it
          endif
        else
          if( Taux .and. Temp_B_iso ) then
            write(3,575) Z, posn(:,igr), it, Taux_oc(jgr), Temp_coef(jgr), norbv(jgr), &
                                                                     ( popats(jgr,l,1:nspin), l = 1,nlat( abs(it) ) )
          elseif( Taux ) then
            write(3,576) Z, posn(:,igr), it, Taux_oc(jgr), norbv(jgr), ( popats(jgr,l,1:nspin), l = 1,nlat( abs(it) ) )
          elseif( Temp_B_iso ) then
            write(3,576) Z, posn(:,igr), it, Temp_coef(jgr), norbv(jgr), ( popats(jgr,l,1:nspin), l = 1,nlat( abs(it) ) )
          else
            write(3,580) Z, posn(:,igr), it, norbv(jgr), (popats(jgr,l,1:nspin), l = 1,nlat( abs(it)) )
          endif
          do io = 1,norbv(jgr)
            if( sum( abs( hybrid(io,10:16,jgr) ) ) < eps10 ) then
              write(3,600) io, pop_nonsph(io,jgr), real( hybrid(io,1:9,jgr) )
            else
              write(3,600) io, pop_nonsph(io,jgr), real( hybrid(io,:,jgr) )
            endif
          end do
        endif
        if( Atom_occ_hubb .and. Hubb( abs(it) ) ) then
          l = l_hubbard( Z )
          write(3,610) ( ( occ_hubb_e(m,m,isp,jgr), m = -l,l ), isp = 1,nspin )
        endif
      end do

      if( abs( dpos(1) ) > epspos .or. abs( dpos(2) ) > epspos .or. abs( dpos(3) ) > epspos ) write(3,545) dpos(1:3)

      if( Film ) then
        write(3,*)
        if( n_atom_per > 0 ) then
          write(3,'(a17,f10.5)') ' Film thickness =', Film_thickness
          write(3,'(a17,3f10.5)') ' Film shift     =', Film_shift(1:3)
          if( abs( Film_shift(4) ) > eps10 ) write(3,'(a17,f10.5)') ' Film rotat     =', Film_shift(4)
          if( Film_zero > -100._db )  write(3,'(a17,f10.5)') ' Film zero      =', Film_zero
        endif
        if( Film_roughness > eps10 ) write(3,'(a17,f10.5)') ' Film roughness =', Film_roughness
        if( n_atom_int > 0 ) then
          write(3,'(a17,3f10.5)') ' Interface shift=', Interface_shift(1:3)
          if( abs( Interface_shift(4) ) > eps10 ) write(3,'(a17,f10.5)') ' Interface rotat=', Interface_shift(4)
        endif
        if( n_atom_sur > 0 ) then
          write(3,'(a17,3f10.5)') ' Surface shift  =', Surface_shift(1:3)
          if( abs( Surface_shift(4) ) > eps10 ) write(3,'(a17,f10.5)') ' Surface rotat  =', Interface_shift(4)
        endif
        if( hkl_film ) then
          write(3,'(A)') ' (h,k,l) corresponding to the reciprocal space of the film'
        else
          write(3,'(A)') ' (h,k,l) corresponding to the reciprocal space of the bulk'
        endif
      endif
      if( Bulk ) then
        if( Bulk_roughness > eps10 ) write(3,'(a17,f10.5)') ' Bulk roughness =', Bulk_roughness
        write(3,620) ' Bulk', axyz_bulk(:), Angxyz_bulk(:)
        if( Taux .and. Atom_B_iso ) then
          mot = ' Occupancy   B_iso'
        elseif( Taux .and. Atom_U_iso ) then
          mot = ' Occupancy   U_iso'
        elseif( Taux ) then
          mot = ' Occupancy'
        elseif( Atom_B_iso ) then
          mot = '   B_iso'
        elseif( Atom_U_iso ) then
          mot = '   U_iso'
        else
          mot = ' '
        endif

        write(3,560) mot
        do igr = 1,n_atom_bulk
          it = itype( ngroup - n_atom_bulk + igr )
          if( Taux .and. Atom_B_iso ) then
            write(3,625) Z_bulk(igr), posn_bulk(:,igr), it, Taux_oc(ngroup-n_atom_bulk+igr), Temp_coef(ngroup-n_atom_bulk+igr)
          elseif( Taux ) then
            write(3,625) Z_bulk(igr), posn_bulk(:,igr), it, Taux_oc(ngroup-n_atom_bulk+igr)
          elseif( Temp_B_iso ) then
            write(3,625) Z_bulk(igr), posn_bulk(:,igr), it, Temp_coef(ngroup-n_atom_bulk+igr)
          else
            write(3,625) Z_bulk(igr), posn_bulk(:,igr), it
          endif
        end do
      endif
      if( Cap_layer ) then
        write(3,620) '  Cap', axyz_cap(:), Angxyz_cap(:)
        write(3,'(A)') '    Z         x              y              z           Occupancy'
        do igr = 1,n_atom_cap
          write(3,630) Z_cap(igr), posn_cap(:,igr), Taux_cap(igr)
        end do
        if( Cap_shift < -100._db ) then
          write(3,'(A)') ' Default cap shift, equal to the sum of top film and bottom cap atom radii'
        else
          write(3,640) 'shift    ', Cap_shift
        endif
        if( Cap_thickness < -100._db ) then
          write(3,'(A)') ' Default cap thickness, equal to one layer'
        else
          write(3,640) 'thickness', Cap_thickness
        endif
        if( Cap_roughness > eps10 ) write(3,640) 'roughness', Cap_roughness
        if( Cap_U_iso > eps10 ) then
          write(3,640) 'U_iso    ', Cap_U_iso
        elseif( Cap_B_iso > eps10 ) then
          write(3,640) 'B_iso    ', Cap_B_iso
        endif
      endif

      if( Center_s ) then
        write(3,'(/A)') '   No surface centering'
        write(3,'(a11,2f10.3)') '   Center =', Centre(1:2)
      elseif( Noncentre ) then
        write(3,'(/A)') '   Cluster not automatically centered on the absorbing atom'
        write(3,'(a11,2f10.3)') '   Center =', Centre(:)
      endif

! About potential
      if( abs( V_helm ) > eps10 ) then
        write(3,650) V_helm, Delta_helm, Width_helm
        if( Helm_cos ) write(3,'(A)') '    with cosine model'
      endif
      if( Flapw ) then
        if( Hedin .or. Perdew ) then
          write(3,'(/A)') ' FLAPW potential energy dependant'
        else
          write(3,'(/A)') ' FLAPW potential not energy dependant'
        endif
      else
        if( Hedin ) then
          write(3,'(/A)') ' Hedin and Lundqvist exchange-correlation potential'
        elseif( Perdew ) then
          write(3,'(/A)') ' Perdew and Zunger exchange-correlation potential'
        else
          write(3,660) alfpot
        endif
      endif
      if( Full_potential ) write(3,670) lmax_pot
      if( neimagent == 1 ) then
        write(3,680) eimagent(1)
      elseif( neimagent > 1 ) then
        write(3,'(A)') '   Energy   E_imag    (eV)'
        write(3,690) (eeient(ie), eimagent(ie), ie = 2,neimagent)
      endif
      if( multrmax /= 1 ) write(3,702) multrmax
      if( abs( rpotmax ) > eps4 ) write(3,703) Rpotmax
      write(3,704) D_max_pot

      if( Hubbard ) then
        write(3,'(/A)') ' Hubbard calculation '
        write(3,'(A)') ' Type  Z  Hubbard parameter (eV)'
        do it = 1,ntype_mod
          if( Hubb(it) ) write(3,'(2i4,f12.3)') it, numat(it), V_hubbard(it)
        end do
      endif
      if( ATA ) write(3,'(/A)') ' Average T-matrix Approximation'

      if( SCF .or. Fermi_auto ) then
        if( Fermi_auto .and. abs(p_self0) < eps10 .and. nself == 1 ) then
          write(3,'(/A)') ' One cycle for the Fermi level calculation'
        else
          write(3,'(/A)') ' Self consistent calculation'
          if( Bulk .and. .not. SCF_bulk ) write(3,'(3x,A)') ' but in the bulk'
          write(3,708) nself, p_self0, p_self_max, Delta_En_conv
        endif
        write(3,709) r_self
        if( Scf_elecabs ) write(3,'(A)') '   Cuting energy criterium using the number of electron in the absorbing atom'
        if( Self_nonexc ) then
          write(3,'(A)') '   Non excited absorbing atom in this part'
        else
          write(3,'(A)') '   Excited absorbing atom in this part'
        endif
      elseif( .not. Extract ) then
        write(3,'(/A)') ' No Fermi level calculation. It is set at its default value = -5 eV'
      end if

    endif

  endif ! arrival point when mpirank0 /= 0

  if( mpinodes0 > 1 ) then
    call MPI_Bcast(Green_bulk,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Green_s,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Extract,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  endif
  call getSolverParams(MPI_host_num_for_mumps,mpinodes0,Solver)
  if( Green_s .or. Extract ) MPI_host_num_for_mumps = 1
  mpirank = mpirank0 / MPI_host_num_for_mumps
  mpinodes = mpinodes0 / MPI_host_num_for_mumps
  mpirank_in_mumps_group = mod( mpirank0, MPI_host_num_for_mumps )
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,mpirank,mpirank_in_mumps_group,MPI_COMM_MUMPS,mpierr)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,mpirank_in_mumps_group,mpirank,MPI_COMM_GATHER,mpierr)

  if( mpirank0 == 0 ) then
    do ipr = 3,6,3
      if( icheck(1) < 1 .and. ipr == 3 ) cycle
      if( Green_s .or. Extract ) then
        if( mpinodes0 > 1 ) then
          write(ipr,'(/A)') ' Parallel calculation'
          write(ipr,720) mpinodes0
        else
          write(ipr,'(/A)') ' Sequential calculation'
        endif
      else
        if( Solver == 'MUMPS' ) then
          if( mpinodes0 > 1 ) then
            write(ipr,'(/A)') ' Parallel calculation with MUMPS Solver'
            write(ipr,710) mpinodes0, mpinodes, MPI_host_num_for_mumps
          else
            write(ipr,'(/A)') ' Sequential calculation with MUMPS Solver'
          endif
        else
          if( mpinodes0 > 1 ) then
            write(ipr,'(/A)') ' Parallel calculation with Gaussian Solver'
            write(ipr,720) mpinodes0
          else
            write(ipr,'(/A)') ' Sequential calculation with Gaussian Solver'
          endif
        endif
      endif
    end do

  endif

  if( mpinodes0 > 1 ) then
    call MPI_Bcast(adimp,n_adimp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(alfpot,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(All_nrixs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(All_site_rixs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(allsite,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ampl_rixs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Angle_mode,3*npldafs_2d,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Angle_or,npldafs_f,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    if( npldafs > 0 ) call MPI_Bcast(angpoldafs,3*npldafs, MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(angxyz,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(angxyz_int,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(angxyz_sur,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ang_base_loc,3*(ntype+1),MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ang_base_loc_gr,3*ngroup,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ang_rotsup,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ang_spin,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Atom_B_iso,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Atom_U_iso,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Atomic_scr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Axe_spin,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(axyz,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(axyz_int,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(axyz_sur,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Basereel,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Bormann,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Bulk_roughness,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Cap_B_iso,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Cap_roughness,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Cap_shift,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Cap_thickness,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Cartesian_tensor,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    if( norbdil > 0 ) call MPI_Bcast(cdil,norbdil,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Centre,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Centre_auto,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Centre_auto_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Center_s,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Charge_free,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Classic_irreg,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Clementi,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(COOP,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(COOP_z_along_bond,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Core_energ_tot,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Core_resolved,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Coupelapw,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(D_max_pot,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Dip_rel,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(dyn_eg,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(dyn_g,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Dafs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Delta_Epsii,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Delta_En_conv,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Delta_helm,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Density,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Density_comp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Dipmag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Dist_coop,2,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Doping,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(dpos,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E_adimp,n_adimp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E_max_range,n_range,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E_radius,n_radius,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E1E1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E1E2e,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E1E3,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E1M1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E1M2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E2E2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(E3E3,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Eclie,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Eclie_out,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ecrantage,nspin,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    if( neimagent > 0 ) call MPI_Bcast(Eeient,neimagent,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Egamme,ngamme,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    if( neimagent > 0 ) call MPI_Bcast(eimagent,neimagent,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Eneg_i,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Eneg_n_i,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Energ_rixs,nenerg_in_rixs,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Energphot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ephot_min,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(f_no_res,2,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(FDM_comp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Film_roughness,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Film_thickness,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Film_shift,4,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Force_ecr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Full_atom_e,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Full_potential,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Gamma_tddft,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Green_int,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Green_self,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Harm_cubic,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Hedin,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    if( npldafs > 0 ) call MPI_Bcast(hkl_dafs,3*npldafs, MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Helm_cos,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(hkl_film,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(hkl_ref,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if( ngroup_nonsph > 0 ) then
      allocate( Hybrid_r(nhybm,16,ngroup_nonsph) )
      allocate( Hybrid_i(nhybm,16,ngroup_nonsph) )
      if( mpirank0 == 0 ) then
        hybrid_r(:,:,:) = real( hybrid(:,:,:),db )
        hybrid_i(:,:,:) = aimag( hybrid(:,:,:) )
      endif
      call MPI_Bcast(Hybrid_r,nhybm*16*ngroup_nonsph,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Hybrid_i,nhybm*16*ngroup_nonsph,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

      if( mpirank0 /= 0 ) Hybrid(:,:,:) = cmplx( Hybrid_r(:,:,:), Hybrid_i(:,:,:),db )
      deallocate( Hybrid_i, Hybrid_r )
    endif
    call MPI_Bcast(iabsm,n_multi_run_e,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(iabsorig,n_multi_run_e,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(igr_coop,n_atom_coop,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(igr_dop,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Interface_shift,4,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(iord,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(isigpi,npldafs*2, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(itdil,norbdil, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(itype,ngroup,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(itype_dop,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(jseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Kern_fac,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Kern_fast,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Kgroup,ngroup_pdb, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(korigimp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lamstdens,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ldil,norbdil,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ldipimp,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lecrantage,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lin_gam,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmax_nrixs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmax_tddft_inp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmax_pot,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmaxat0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmaxfree,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmaxso0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmaxso_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmoins1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lplus1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lquaimp,9,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Lseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lvval,(ntype+1)*nlatm, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(M1M1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(M1M2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(M2M2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Mat_or,9,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Mat_UB,9,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(matper,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Moment_conservation,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(muffintin,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(multrmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_atom_proto,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(n_devide,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(nbseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nchemin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(necrantage,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(No_dft,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(No_solsing,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(noncentre,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Nonexc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    if( ngroup_nonsph > 0 ) call MPI_Bcast(norbv,ngroup_nonsph+1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(normaltau,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(normrmt,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nphim,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(npl_2d,npldafs_2d,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nself,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nself_bulk,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nsymextract,n_multi_run_e,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nposextract,n_multi_run_e,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ntype_mod,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if( nlatm > 0 ) call MPI_Bcast(nvval,(ntype+1)*nlatm, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(numat_abs,n_Z_abs,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if( Atom_occ_hubb ) call MPI_Bcast(occ_hubb_e,nspin*ngroup_hubb *(2*m_hubb_e+1)**2,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Octupole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Old_zero,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(One_run,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(One_SCF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Operation_mode,5*npldafs_2d,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Operation_mode_used,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Optic,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Overad,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Overlap,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(p_self_max,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(p_self0,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(p_self0_bulk,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(pdpolar,nple*2,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(polar,3*nple,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(polarise,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(poldafse,3*npldafs_e,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(poldafsei,3*npldafs_e,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(poldafss,3*npldafs_e,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(poldafssi,3*npldafs_e,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    if( ngroup_nonsph > 0 ) call MPI_Bcast(pop_nonsph, nhybm*ngroup_nonsph,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(popats,ngroup*nlatm*nspin,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(popval,(ntype+1)*nlatm*nspin,MPI_REAL8, 0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(phi_0,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    do i = 1,8
      if( mpirank0 == 0 ) j = iachar( PointGroup(i:i) )
      call MPI_Bcast(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

      if( mpirank0 /= 0 ) PointGroup(i:i) = achar( j )
    end do
    call MPI_Bcast(posn,3*n_atom_uc,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Powder,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(q_rixs,4*n_q_rixs,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Quadmag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Quadrupole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ray_max_dirac,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(R_self,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(R_self_bulk,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Relativiste,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Renorm,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(roverad,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(rpotmax,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(R_rydb,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Rsorte_s,n_radius,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Rydberg,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(SCF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(SCF_bulk,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(SCF_elecabs,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Self_cons,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Self_cons_bulk,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Self_nonexc,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    do i = 1,3
      if( mpirank0 == 0 ) j = iachar( Seuil(i:i) )
      call MPI_Bcast(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

      if( mpirank0 /= 0 ) Seuil(i:i) = achar( j )
    end do
    call MPI_Bcast(Solsing_s,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Solsing_only,1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Spherical_tensor,1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Spherical_signal,1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Spin_channel,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Spinorbite,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(State_all,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(State_all_out,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Step_loss,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Supermuf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Surface_shift,4,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Temp_coef,ngroup_temp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Trace_k,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Trace_p,6,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Theta_in,n_theta_in,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Two_theta,n_two_theta,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Vecdafsem,3*npldafs_e,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Vecdafssm,3*npldafs_e,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(PointGroup_Auto,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Symauto,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Symmol,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Tddft,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Test_dist_min,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(q_nrixs,4*nq_nrixs,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(RIXS,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(RIXS_core,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(RIXS_fast,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Rpalf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Sum_rixs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Temp,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Taux_cap,n_atom_cap,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Taux_oc,ngroup_taux,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Vec_orig,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(V_intmax,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(V_helm,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(V_hubbard,ntype+1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(V0bdcFimp,nspin,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    if( nple > 0 ) call MPI_Bcast(Veconde,3*nple,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Width_helm,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ylm_comp_inp,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Z_nospinorbite,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    if( Bulk ) then
    call MPI_Bcast(axyz_bulk,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(angxyz_bulk,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(posn_bulk,3*n_atom_bulk,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Z_bulk,n_atom_bulk,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    endif
    if( Cap_layer ) then
    call MPI_Bcast(axyz_cap,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(angxyz_cap,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(posn_cap,3*n_atom_cap,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Z_cap,n_atom_cap,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    endif

    n = ntype + 1

    call MPI_Bcast(icom,n,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nlat,n,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Hubb,n,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nrato,n,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nrato_dirac,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(nrm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(numat,n,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Pas_SCF,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(rchimp,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(rmt,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(rmtimp,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Write_modQ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)


    if( Flapw ) then
      call MPI_Bcast(its_lapw,ngroup_lapw,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(nrato_lapw,n,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(r0_lapw,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(rlapw,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(rotloc_lapw,9*ngroup_lapw,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Wien_matsym,9*nslapwm,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Wien_taulap,3*nslapwm,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      do i = 1,Length_name
        if( mpirank0 == 0 ) j = iachar( Wien_file(8)(i:i) )
        call MPI_Bcast(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

        if( mpirank0 /= 0 ) Wien_file(8)(i:i) = achar( j )
      end do
    endif


    if( mpirank0 /= 0 ) icheck(:) = 0

  endif

  if( mpirank0 == 0 ) then
    istop = 0
    do k = 1,3
      if( abs(axyz(k)) < epspos ) then
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,740) k
        end do
        istop = 1
      endif
    end do

!    if( One_run .and. Taux ) then
!      do igr = 1,ngroup_taux
!        if( abs( Taux_oc(igr) - 1 ) > 0.000001_db ) exit
!      end do
!      if( igr <= ngroup_taux ) then
!        Error_message = ' It is not possible to have at the same time One_run mode and occupancy different from 1 !'
!        call write_error_message(Error_message,ipr0,istop)
!      endif
!    endif

    if( Film_roughness > eps10 ) then
      do igr = n_atom_per + 1, n_atom_uc
        if( posn(3,igr) < 1 + eps10 ) cycle
        Error_message = ' It is not possible to have at the same time film roughness and overlayer !'
        call write_error_message(Error_message,ipr0,istop)
        exit
      end do
    endif

    if( istop == 1 ) stop
  endif

! Conversion in atomic units (bohr and Rydberg) and in radian
  adimp(:) = adimp(:) / bohr
  Angle_mode(:,:) = Angle_mode(:,:) * radian
  angxyz(:) = angxyz(:) * radian
  angxyz_bulk(:) = angxyz_bulk(:) * radian
  angxyz_cap(:) = angxyz_cap(:) * radian
  angxyz_int(:) = angxyz_int(:) * radian
  angxyz_sur(:) = angxyz_sur(:) * radian
  where( abs(angpoldafs) < 9999._db ) angpoldafs = angpoldafs * radian
  axyz(1:3) = axyz(1:3) / bohr
  axyz_int(1:3) = axyz_int(1:3) / bohr
  axyz_sur(1:3) = axyz_sur(1:3) / bohr
  Bulk_roughness = Bulk_roughness / bohr
  if( Bulk ) axyz_bulk(1:3) = axyz_bulk(1:3) / bohr
  if( Cap_layer ) axyz_cap(1:3) = axyz_cap(1:3) / bohr
  Cap_roughness = Cap_roughness / bohr
  Cap_thickness = Cap_thickness / bohr
  Cap_shift = Cap_shift / bohr
  D_max_pot = D_max_pot / bohr
  Delta_En_conv = Delta_En_conv / Rydb
  Delta_Epsii = Delta_Epsii / Rydb
  Delta_helm = Delta_helm / bohr
  Dist_coop(:) = Dist_coop(:) / bohr
  dpos(:) = dpos(:) / bohr
  E_adimp(:) = E_adimp(:) / Rydb
  E_max_range(:) = E_max_range(:) / Rydb
  E_radius(:) = E_radius(:) / Rydb
  Eclie = Eclie / Rydb
  Eclie_out = Eclie_out / Rydb
  Eeient(:) = Eeient(:) / Rydb
  Egamme(:) = Egamme(:) / Rydb
  Eimagent(:) = Eimagent(:) / Rydb
  Energ_rixs(:) = Energ_rixs(:) / Rydb
  Ephot_min = Ephot_min / Rydb
  Film_roughness = Film_roughness / bohr
  Film_thickness = Film_thickness / bohr
  Film_shift(1:3) = Film_shift(1:3) / bohr
  Film_shift(4) = Film_shift(4) * radian
  Interface_shift(1:3) = Interface_shift(1:3) / bohr
  Interface_shift(4) = Interface_shift(4) * radian
  Mat_UB(:,:) = Mat_UB(:,:) * bohr  ! Mat_UB in A^-1 in inpout
  Pas_SCF = Pas_SCF / rydb
  phi_0 = phi_0 * radian
  q_nrixs(4,:) = q_nrixs(4,:) * bohr  ! q_nrxis in A^-1 in inpout
  q_rixs(4,:) = q_rixs(4,:) * radian
  R_rydb = R_rydb / bohr
  R_self = R_self / bohr
  R_self_bulk = R_self_bulk / bohr
  Ray_max_dirac = Ray_max_dirac / bohr
  rchimp(:) = rchimp(:) / bohr
  rmtimp(:) = rmtimp(:) / bohr
  Roverad = Roverad / bohr
  Rpotmax = Rpotmax / bohr
  Rsorte_s(:) = Rsorte_s(:) / bohr
  Step_loss = Step_loss / Rydb
  Surface_shift(1:3) = Surface_shift(1:3) / bohr
  Surface_shift(4) = Surface_shift(4) * radian
  Test_dist_min = Test_dist_min / bohr
  Two_theta(:) = Two_theta(:) * pi / 180_db
  Theta_in(:) = Theta_in(:) * pi / 180_db
  V_helm = V_helm / Rydb
  V_intmax = V_intmax / Rydb
  if( Hubbard ) V_hubbard(:) = V_hubbard(:) / Rydb
  V0bdcFimp(:) = V0bdcFimp(:) / Rydb
  Width_helm = Width_helm / bohr

  if( Atom_U_iso ) Temp_coef(:) = 8 * pi**2 * Temp_coef(:)

  call cal_cubmat(angxyz,Cubmat,Struct)

  if( Magnetic .or. Atom_nonsph ) then

    call invermat(Cubmat,Cubmati)

    if( abs( Axe_spin(1) ) < eps10 .and. abs( Axe_spin(2) ) < eps10 .and. abs( Axe_spin(3) - 1 ) < eps10 ) then
      Axe_spin(1:2) = 0._db; Axe_spin(3) = 1._db
      Axe_spin = matmul( Cubmat, Axe_spin )
      call mat_euler( Ang_spin, Rot_gen )
      Axe_spin = matmul( Rot_gen, Axe_spin )
    else
      Axe_spin(:) = Axe_spin(:) * axyz(:) * bohr
      Axe_spin = matmul( Cubmat, Axe_spin )
      vv = sqrt( sum( Axe_spin(:)**2 ) )
      Axe_spin = Axe_spin / vv
      if( abs( Axe_spin(3) ) < 1 - eps10 ) then
        Ang_spin(1) = datan2( Axe_spin(2), Axe_spin(1) )
        Ang_spin(2) = acos( Axe_spin(3) )
        Ang_spin(3) = 0._db
      else
        Ang_spin(:) = 0._db
      endif
      call mat_euler( Ang_spin, Rot_gen )
    endif
    Axe_spin = matmul( Cubmati, Axe_spin )
    Axe_spin(:) = Axe_spin(:) / axyz(:)

    do igr = 1,ngroup
      if( Ang_base_loc_gr(1,igr) > - 1000._db ) cycle
      if( Ang_base_loc(1, abs(itype(igr)) ) < -1000._db ) cycle
      Ang_base_loc_gr(:,igr) = Ang_base_loc(:, abs(itype(igr)) )
    end do
    do igr = 1,ngroup - n_atom_bulk

      if( igr == n_atom_per + n_atom_int + 1 .and. n_atom_per + n_atom_int /= 0 .and. n_atom_sur /= 0 ) then
        call cal_cubmat(angxyz_sur,Cubmat,Struct)
        call invermat(Cubmat,Cubmati)
      elseif( igr == n_atom_per + 1 .and. n_atom_int /= 0 ) then
        call cal_cubmat(angxyz_int,Cubmat,Struct)
        call invermat(Cubmat,Cubmati)
      endif

      if( Ang_base_loc_gr(1,igr) < - 1000._db ) Ang_base_loc_gr(:,igr) = Ang_spin(:)
      Ang(:) = Ang_base_loc_gr(:,igr)
      call mat_euler( Ang, Rot )
      Rot_Atom_gr(:,:,igr) = Rot(:,:)
      if( Struct == 'trigo' ) then
        Axe(1:3) = 1._db / sqrt(3._db)
      else
        Axe(1:2) = 0._db; Axe(3) = 1._db
      endif
      Axe = matmul( Cubmat, Axe )
      Axe = matmul( Rot, Axe )
      Axe = matmul( Cubmati, Axe )
      Axe_atom_gr(:,igr) = Axe(:) / axyz(:)
      if( itype(igr) < 0 ) Axe_atom_gr(:,igr) = - Axe_atom_gr(:,igr)

      p(:) = Axe_atom_gr(:,igr) * axyz(:)
      p = matmul( cubmat, p )
      vv = sum( p(:)**2 )
      if( vv > eps10 ) p(:) = p(:) / sqrt( vv )
      p = matmul( Cubmati, p )
! Axe_atom_gr is now normalized to 1
      Axe_atom_gr(:,igr) = p(:) / axyz(:)
    end do

    if( Bulk ) then
      call cal_cubmat(angxyz_bulk,Cubmat_bulk,Struct_bulk)
      call invermat(Cubmat_bulk,Cubmati_bulk)

      do igr = ngroup - n_atom_uc + 1, ngroup
        if( Ang_base_loc_gr(1,igr) < - 1000._db ) Ang_base_loc_gr(:,igr) = Ang_spin(:)
        Ang(:) = Ang_base_loc_gr(:,igr)
        call mat_euler( Ang, Rot )
        Rot_Atom_gr(:,:,igr) = Rot(:,:)
        if( Struct_bulk == 'trigo' ) then
          Axe(1:3) = 1._db / sqrt(3._db)
        else
          Axe(1:2) = 0._db; Axe(3) = 1._db
        endif
        Axe = matmul( Cubmat_bulk, Axe )
        Axe = matmul( Rot, Axe )
        Axe = matmul( Cubmati_bulk, Axe )
        Axe_atom_gr(:,igr) = Axe(:) / axyz_bulk(:)
        if( itype(igr) < 0 ) Axe_atom_gr(:,igr) = - Axe_atom_gr(:,igr)

        p(:) = Axe_atom_gr(:,igr) * axyz_bulk(:)
        p = matmul( Cubmat_bulk, p )
        vv = sum( p(:)**2 )
        if( vv > eps10 ) p(:) = p(:) / sqrt( vv )
        p = matmul( Cubmati_bulk, p )
! Axe_atom_gr is now normalized to 1
        Axe_atom_gr(:,igr) = p(:) / axyz_bulk(:)
      end do

    endif

    if( icheck(1) > 0 ) then
      write(3,750) Axe_spin(1:3) * axyz(1:3)
      write(3,760) Ang_spin(1:3) * 180 / pi
      do igr = 1,ngroup - n_atom_bulk
        write(3,770) igr
        do i = 1,3
          write(3,780) Rot_atom_gr(i,:,igr), Axe_atom_gr(i,igr) * axyz(i)
        end do
      end do
      do igr = ngroup - n_atom_bulk + 1, ngroup
        write(3,770) igr
        do i = 1,3
          write(3,780) Rot_atom_gr(i,:,igr), Axe_atom_gr(i,igr) * axyz_bulk(i)
        end do
      end do
    endif
  endif

  Atom_nsph_e(:) = .false.
  if( Atom_nonsph ) then
    do igr = 1,ngroup
      if( norbv(igr) == 0 ) cycle
      pop_nsph = sum( pop_nonsph(1:norbv(igr),igr) )
      if( pop_nsph > eps10 ) Atom_nsph_e(igr) = .true.
    end do
  endif

  if( nq_nrixs >= 1 .and. lmax_nrixs > 1 ) Quadrupole = .true.
  if( nq_nrixs >= 1 .and. lmax_nrixs > 2 ) Octupole = .true.

  if( Dafs_bio ) then
    nphim = 4
    nphi_dafs(:) = nphim
  else
    nphimt = 1
    nphi_dafs(:) = 1
    do ipl = 1,npldafs
      if( angpoldafs(1,ipl) > -9999._db .and. angpoldafs(2,ipl) > -9999._db .and. angpoldafs(3,ipl) > -9999._db ) cycle
      nphimt = nphim
      nphi_dafs(ipl) = nphim
    end do
    nphim = nphimt
  endif

  if( Dafs ) then
    poldafsem(:,:) = cmplx( poldafse(:,:), poldafsei(:,:), db )
    poldafssm(:,:) = cmplx( poldafss(:,:), poldafssi(:,:), db )
  endif

  if( Centre_auto ) call Auto_center(axyz,Centre,Centre_auto_abs,Cubmat,icheck(1),itype,n_atom_uc,n_Z_abs,ngroup,ntype,numat, &
                                     numat_abs,posn,Rsorte_s(1),Struct)

  Tensor_imp(1:3) = ldipimp(1:3)
  k = 3
  do i = 1,3
    do j = 1,3
      k = k + 1
      Tensor_imp(k) = lquaimp(i,j)
    end do
  end do

  call Rmt_fix(icheck(1),ntype,numat,Rmt)

! Stocking for shorter transmission
  Multipole(1) = E1E1; Multipole(2) = E1E2e; Multipole(3) = E1E3;
  Multipole(4) = E1M1; Multipole(5) = E1M2; Multipole(6) = E2E2;
  Multipole(7) = E3E3; Multipole(8) = M1M1; Multipole(9) = M1M2;
  Multipole(10) = M2M2

  SCF_log(1) = SCF_elecabs
  SCF_log(2) = SCF_mag_free
  SCF_log(3) = Self_cons
  SCF_log(4) = Self_nonexc
  SCF_log(5) = SCF
  SCF_log(6) = Self_cons_bulk
  SCF_log(7) = SCF_bulk

  return
  100 format(//'  Error in the indata file :')
  110 format(//' After the keyword Dafs (or RXS),',/ &
               '    for the reflexion number',i3,/ &
               '    there are at least 2 directions for an angular scan.',/ &
               ' Just one is allowed !'//)
  120 format(//' After the keyword Dafs (or RXS),',/ &
               '    the polarization is not defined for the reflexion number',i3,'!')
  125 format(/' After the keyword Atom_conf, for the atom number',i2,',',/ &
              ' check how is written the corresponding electronic configuration !'//)
  130 format(/' After the keyword Atom, for the atom number',i2,',',/  &
              ' check how is written the corresponding electronic configuration !'//)
  140 format(/' Error in the Pdb file :',// &
              ' The chemical symbol of atom number',i6,' is: ',a2,', what is not known!'//)
  145 format(/' Error in the Cif file :',// &
              ' The chemical symbol of atom number',i6,' is: ',a2,', what is not known!'//)
  150 format(//' The following line is not understood :',/A,// &
               ' If it is a keyword, check the spelling.'/, &
               ' If the line is not supposed to be a keyword but contains numbers, check:'/ &
               '      - How many numbers must be in the line ?'/, &
               '      - Are there spaces between the numbers ?'//)
  155 format(/'  When "Doping" and "Atom" keywords are used together,',/ &
              '  the first number after "Doping" must be the atom type index and not the atomic number !',/ &
              '  This index is the line number in the list under "Atom".',// &
              '  In the indata file, this index is ',i4,/ &
              '  which is higher that the total number of atomic configurations =',i2,', listed after "Atom"')
  160 format(//' Error under the keyword Par_',a6,/ &
                ' The wanted atom is the number',i3,' !',/ &
                ' There are only',i3,' atoms in the job !'//)
  170 format(//'  A parameter index for the fit is not possible !'/, &
               '  Check your indata file under the keyword ',a9,/ &
               '  The index is',i4//)
  174 format(//' The atomic number given under keyword Z_absorber =',i4,' is not in the list of surface atoms !'//)
  175 format(//' The atomic number given under keyword Z_absorber =',i4,' is not in the list of atoms !'//)
  176 format(//' The first type number defined under "Crystal", "Film" or "Molecule" =',i4,',',/ &
               ' is higher than the number of atom types defined under "Atom" or "Atom_conf" =',i4,' !'//)
  180 format(/' Edge = ',a3,' not programmed !'//)
  190 format(/' Doping atom is indexed with',i3,' what is more than',/ &
              ' the number of atom type =',i2,' declared under keyword Atom !')
  195 format(/' Atom number',i4,' is indexed with',i3,' what is more than',/ &
              ' the number of atom type =',i2,' declared under keyword Atom !')
  200 format(//' Under keyword Crystal or Molecule, when the index is',/ &
               ' supposed to be the atomic number, the simultaneous',/ &
               ' use of the keyword Atom is forbidden !'//)
  210 format(//' Under keyword "Atom" or "Atom_conf", the orbital occupancy is defined higher than what is possible:',/&
               '   This concerns the atom type',i2,', the orbital number',i2, /&
               '   with l =',i2,' and occupancy =',f8.5,' >',i3,' !'//)
  220 format(//' iord =',i2,' must be equal to 2 or 4 !'//)
  230 format(/' itype(igr=',i2,') =',i2,' > ntype =',i2,', forbidden !')
  240 format(//' Bad value in the indata file !'/, &
               ' Under the keyword Polarise, polarisation and wave vector (if specified) cannot be both zero !')
  250 format(//' Bad number of values in the indata file !'/, &
               ' Under the keyword Rmtg, the number of radius values must be equal to the number of atom type !',/ &
               ' In the indata file number of radius =',i3,/ &
               '                 number of atom type =',i3,/)
  260 format(//' Bad polarization definition for DAFS when using Full_self_abs option !'//, &
               ' For reflection number',i3,', (h,k,l) = (',3i3,')',/ &
               ' polarization indexes are',i2,' and',i2,' !'//, &
               ' For each reflection (h,k,l) one must have in order',/ &
              2x,' sigma-sigma, sigma-pi, pi-sigma and pi-pi,',/, &
              2x,' that is respectively indexes 1 1, 1 2, 2 1 and 2 2',/)
  263 format(//' This reflexion does not exist !',/ &
               ' The supercell from the ',a4,' unit cell is ',i3,'.sqrt(2) x ',i3,'.sqrt(2), R(45)')
  265 format(//' This reflexion does not exist !',/ &
               ' The supercell from the ',a4,' unit cell is ',i3,' x ',i3)
  270 format('   This is not compatible with reflection',i4,' : (h, k, l) =',3f10.5)
  280 format('   It is not the case for reflection',i4,' : (h, k, l) =',3f10.5)
  300 format(/' Tensors extracted from the file:'/,A,/)
  305 format(/' Electronic structure extracted from the file:'/,A,/)
  310 format(' Radius =',f6.2)
  315 format(' Radius, E_radius =',100(f6.2,f6.1))
  320 format('   Roverad =',f6.2)
  325 format('   In optics and tddft parts, expansion in spherical harmonics limited to lmax =',i2)
  340 format(' Linear range :')
  350 format(' Range =',40f8.3,5(/8x,40f8.3))
  360 format(' Range =',40f10.3,5(/8x,40f10.3))
  400 format(' Spin-orbit not taken into account for the atomic number', i3)
  410 format('   Path expansion, n =',i3)
  420 format('   Imposed type muffin-tin radius, Rmtimp =',10f6.3,/9x,10f6.3)
  430 format('   Overlap of the muffin-tin radius  =',f6.2)
  440 format('   iord =',i2,', adimp =',f6.2)
  445 format('   iord =',i2,', adimp, E_adimp =',100( f6.2, f6.1) )
  450 format('   lmaxso0 =',i3)
  455 format('   lmaxso_max =',i5)
  460 format('   R_rydb =',f7.3,' A')
  480 format('   E_ato_min, E_out_min =',2f7.3,' eV')
  490 format('   V_intmax =',f7.3,' eV')
  500 format(/' Temperature =',f6.1,' K')
  502 format(' Ray max Dirac  =',f6.2,' A')
  505 format(6x,3f7.3,3x,3f7.3,3x,f7.3,4x,f7.3)
  510 format(6x,3f7.3,3x,3f7.3,3x,'Matrix')
  511 format(/' Bormann',/ '  (h, k, l) = (',i3,',',i3,',',i3,')   Azimuth =',f7.2)
  512 format(7x,5f7.3,3x,5i2,4x,3f9.3)
  513 format(6x,3f9.5,3x,3f9.5,3x,A)
  515 format(7x,3f7.3,3x,a11,2f10.3,'      scan')
  516 format(7x,3f7.3,3x,a11,'      scan',2f10.3)
  517 format(7x,3f7.3,3x,a11,f10.3,'          scan',f10.3)
  518 format(7x,3f7.3,3x,a11,3f10.3)
  519 format(7x,3i3,3x,a11,2f10.3,'      scan')
  520 format(7x,3i3,3x,a11,'      scan',2f10.3)
  521 format(7x,3i3,3x,a11,f10.3,'          scan',f10.3)
  522 format(7x,3i3,3x,a11,3f10.3)
  542 format('   ngroup =',i5,', ntype =',i2)
  543 format(1x,2i4,10(i4,i3,f7.3))
  544 format(1x,2i4,10(i4,i3,2f7.3))
  545 format('   dpos =',3f7.3)
  548 format(' Orbital dilatation :',/'   it    l   cdil')
  549 format(2i5,f7.3)
  552 format('   Point Group = ',a8)
  554 format('   Space Group = ',a13)
  555 format('   a, b, c =',3f12.7)
  556 format('   alfa, beta, gamma =',3f9.3)
  557 format('   alfa, beta, gamma =',3f11.5)
  560 format('    Z         x              y              z      Typ ',A)
  565 format(i5,3f15.10,i4,i6,2f10.5)
  570 format(i5,3f15.10,i4,i6,f10.5)
  571 format(i5,3f15.10,i4,2f10.5)
  572 format(i5,3f15.10,i4,f10.5)
  575 format(i5,3f15.10,i4,2f10.5,2x,12f8.4)
  576 format(i5,3f15.10,i4,f10.5,2x,12f8.4)
  580 format(i5,3f15.10,i4,i5,2x,12f8.4)
  600 format(9x,'Pop_nsph(',i1,') =',f7.3,', Hybrid =',f8.3,1x,3f7.3,1x,5f7.3,1x,7f7.3,1x,100f7.3)
  610 format('    Occ. matrix :', 14f5.2)
  620 format(/1x,a5,' : '/,'   a, b, c =',3f12.7,/'   alfa, beta, gamma =',3f9.3)
  625 format(i5,3f15.10,i4,f10.5)
  630 format(i5,3f15.10,f14.5)
  640 format(' Cap ',a9,' = ',f10.5,' A')
  650 format(/' V_helm =',f10.5,' eV,  Delta_helm =',f10.5,' A,  Width_helm =',f10.5,' A')
  660 format(/' Xalfa potential , Xalfa =',f8.5)
  670 format(' Full potential inside the atomic spheres with lmax =',i2)
  680 format('   E_imag =',f9.3,' eV')
  690 format(2f9.3)
  702 format('   multrmax =',i2)
  703 format('   Rpotmax =',f7.3,' A')
  704 format('   D_max_pot =',f7.3,' A')
  708 format('   Maximum number of iteration = ',i3,/ '   Weight =',f6.3,',   Weight max =',f6.3/ &
             '   Delta energy for convergence =',f7.3,' eV / atom')
  709 format('   Cluster radius used for this part = ',f7.3,' A')
  710 format('  with',i3,' processors including',i3,' for the energy loop and',i3,' for MUMPS')
  720 format('  with',i3,' processors')
  740 format(//' The unit cell parameter',i2,' is zero !'//)
  750 format(/' General Z axis =',3f9.5)
  760 format(' Euler angles   =',3f9.3)
  770 format(/6x,' local matrix rotation',7x,'local Z axis   Atom =',i3)
  780 format(3x,3f9.5,5x,f9.5)

end

!***********************************************************************

subroutine write_error_message(Error_message,ipr0,istop)

  implicit none

  character(len=132):: Error_message
  integer:: ipr, ipr0, istop

  if( istop == 0 ) call write_error
  do ipr = ipr0,9,3
    write(ipr,100)
    write(ipr,'(/A/)') Error_message
  end do
  istop = 1

  return
  100 format(//'  Error in the indata file :')
end

!***********************************************************************

subroutine Extract_log(Core_resolved,Green_bulk,Green_s,nom_fich_extract,Optic,Renorm,Seuil)

  use declarations
  implicit none

  integer:: i, istat, Length, Line

  character(len=3):: Seuil
  character(len=132):: mot
  character(len=Length_name):: nom_fich_extract

  logical:: Core_resolved, Green_bulk, Green_s, Optic, Renorm

  open(1, file = nom_fich_extract, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(nom_fich_extract,istat,1)

  do i = 1,100000
    read(1,'(A)') mot
    if( mot(2:10) == 'Threshold' ) then
      Length = len_trim(mot)
      if( mot(Length:Length) == 'e' ) then
        Seuil = mot(Length-7:Length-5)
        Seuil = adjustl( Seuil )
      elseif( mot(Length:Length) == 's' ) then
        Seuil = mot(Length-8:Length-6)
        Seuil = adjustl( Seuil )
      else
        Seuil = mot(14:16)
      endif
      exit
    endif
    if( Seuil == 'Opt' ) Optic = .true.
  end do

  Rewind(1)
  Core_resolved = .false.
  do i = 1,100
    read(1,'(A)') mot
    if( mot(2:5) /= 'Core'  ) cycle
    Core_resolved = .true.
    exit
  end do

  Rewind(1)

  do Line = 1,100000
    read(1,'(A)' ) mot
    if( mot(2:20) == 'Multiple scattering' ) then
      Green_s = .true.
      exit
    elseif( mot(2:25) == 'Finite difference method' ) then
      Green_s = .false.
      read(1,'(A)' ) mot
      if( mot(5:19) == 'but in the bulk' ) then
        Green_bulk = .true.
      else
        backspace(1)
      endif
      exit
    endif
  end do

  do Line = 1,100000
    read(1,'(A)' ) mot
    if( mot(2:50) == 'Normalization of the atomic radial wave functions' ) then
      Renorm = .true.
      exit
    elseif( mot(2:53) == 'No normalization of the atomic radial wave functions' ) then
      Renorm = .false.
      exit
    endif
  end do

  Close(1)

  return
end

!***********************************************************************

! Muffin-tin radius fixation for FDM method

subroutine Rmt_fix(icheck,ntype,numat,Rmt)

  use declarations
  implicit none

  integer:: icheck, it, ntype
  integer, dimension(0:ntype):: numat

  real(kind=db):: p1, Rmt_H, Rmt_Ti
  real(kind=db), dimension(0:ntype):: Rmt

  Rmt_H = 0.3_db
  Rmt_Ti = 1.0_db

  do it = 1,ntype
    if( abs(rmt(it)) > eps10 ) cycle
    if( numat(it) == 0 ) then
      Rmt(it) = 0._db
    elseif( numat(it) == 1 ) then
      Rmt(it) = rmt_H
    elseif( numat(it) > 21 ) then
      Rmt(it) = rmt_Ti
    else
      p1 = ( numat(it) - 1._db ) / 21._db
      Rmt(it) = p1 * Rmt_Ti + (1 - p1) * Rmt_H
    endif
  end do

  if( icheck > 0 ) then
    write(3,110)
    do it = 1,ntype
      write(3,120) it, numat(it), Rmt(it)
    end do
  endif
  Rmt(:) = Rmt(:) / bohr

  return
  110 format(/' FDM atom radius: Type   Z    R')
  120 format(16x,2i5,f7.3)
end

!***********************************************************************

function number_from_text(n_skip,mot_in)

  use declarations
  implicit none

  character(len=132):: mot, mot_in

  integer:: i, ier, j, length, n, n_skip

  real(kind=db):: number_from_text

  mot = ' '
  mot = mot_in

  do n = 1,n_skip
    mot = adjustl(mot)
    length = len_trim(mot)

    boucle_i: do i = 1,length
      if( mot(i:i) == ' ' ) then
        do j = 1,i-1
          mot(j:j) = ' '
        end do
        exit boucle_i
      endif
    end do boucle_i
    if( n == n_skip ) exit
  end do

  mot = adjustl(mot)
  length = len_trim(mot)

  do i = 1,length
    if( mot(i:i) == ' ' .or. mot(i:i) == '(' ) exit
  end do

  open(9, status='SCRATCH')
  write(9,*) mot(1:i-1)
  backspace(9)
  read(9,*,iostat=ier) number_from_text
  if( ier /= 0 ) number_from_text = 0._db
  Close(9)

  return
end

!***********************************************************************

function word_from_text(n_skip,mot_in)

  use declarations
  implicit none

  character(len=132):: mot, mot_in, word_from_text

  integer:: i, j, length, n, n_skip

  word_from_text = ' '
  mot = ' '
  mot = mot_in

  do n = 1,n_skip
    mot = adjustl(mot)
    length = len_trim(mot)

    boucle_i: do i = 1,length
      if( mot(i:i) == ' ' ) then
        do j = 1,i-1
          mot(j:j) = ' '
        end do
        exit boucle_i
      endif
    end do boucle_i
    if( n == n_skip ) exit
  end do

  mot = adjustl(mot)
  length = len_trim(mot)

  do i = 1,length
    if( mot(i:i) == ' ' .or. mot(i:i) == '(' ) exit
  end do

  word_from_text(1:i-1) = mot(1:i-1)

  return
end

!***********************************************************************

! Lecture of the structure coming from WIEN-2k

subroutine lect_struct_lapw(angxyz,axyz,iabsm,iabsorig,icheck,its_lapw,itype,n_multi_run_e,n_atom_per,ngroup_lapw,nomstruct, &
        nrato_lapw,nrm,nslapwm,ntype,ntype_bulk,numat,posn,r0_lapw,rlapw,rotloc_lapw,Wien_matsym,Wien_taulap)

  use declarations
  implicit none

  character(len=1) Trans
  character(len=2) Plan
  character(len=Length_name) nomstruct

  integer:: i, ia, icheck, igr, index, ipr, iprot, is, istat, it, itr, its, ittt, j, jatom, mu, mult, multi_run, n_multi_run_e, &
    n_atom_per, ngroup_lapw, nrm, nmatsym, nslapwm, nt, ntrans, ntype, ntype_bulk

  integer, dimension(n_multi_run_e):: iabsm, iabsorig
  integer, dimension(n_atom_per):: numprot
  integer, dimension(0:ntype):: nrato_lapw, numat
  integer, dimension(ngroup_lapw):: its_lapw
  integer, dimension(n_atom_per):: itype
  integer, dimension(3,3,nslapwm):: Wien_matsym

  real(kind=db):: demi, rZ
  real(kind=db), dimension(3):: angxyz, axyz, dp, p
  real(kind=db), dimension(3,3):: Mas, Mat, ptrans, Rot
  real(kind=db), dimension(3,n_atom_per):: posn
  real(kind=db), dimension(3,3,ntype):: rotloc
  real(kind=db), dimension(3,3,ngroup_lapw):: rotloc_lapw
  real(kind=db), dimension(0:ntype):: r0_lapw, rlapw
  real(kind=db), dimension(3,nslapwm):: Wien_taulap

! ll                = number of (l,m) per atom
! ntype-ntype_bulk  = number of inequivalent atoms
! nmatsym           = number of symmetriy operations
!
! Unite 8 : 'case.struct'  (donnees structurales)

  open(8, file = nomstruct, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(nomstruct,istat,1)

! Lecture of xxxx.struct
! Read and formats partly taken in main1.f routine from Wien97/SRC_lapw5

  read(8,*)
  read(8,'(a1,a2)') Trans, Plan
  demi = 0.5_db
  ptrans(:,:) = 0._db
  select case(Trans)
    case('F')
      ntrans = 3
      ptrans(1,1) = demi;    ptrans(2,1) = demi;    ptrans(3,1) = 0._db
      ptrans(1,2) = demi;    ptrans(2,2) = 0._db;   ptrans(3,2) = demi
      ptrans(1,3) = 0._db;   ptrans(2,3) = demi;    ptrans(3,3) = demi
    case('B')
      ntrans = 1
      ptrans(1:3,1) = demi
    case('C')
      ntrans = 1
      if( Plan == 'XY' ) then
        ptrans(1:2,1) = demi
      elseif( Plan == 'YZ' ) then
        ptrans(1:3:2,1) = demi
      else
        ptrans(2:3,1) = demi
      endif
  case default
      ntrans = 0
      ptrans(1:3,1) = 0._db
  end select

  if( Trans == 'H' ) then
    Mat(1,1) = 1._db; Mat(1,2) = -1._db / 2;      Mat(1,3) = 0._db
    Mat(2,1) = 0._db; Mat(2,2) = sqrt(3._db) / 2; Mat(2,3) = 0._db
    Mat(3,1) = 0._db; Mat(3,2) = 0._db;           Mat(3,3) = 1._db
    Mas(1,1) = 1._db; Mas(1,2) = 1 / sqrt(3._db); Mas(1,3) = 0._db
    Mas(2,1) = 0._db; Mas(2,2) = 2 / sqrt(3._db); Mas(2,3) = 0._db
    Mas(3,1) = 0._db; Mas(3,2) = 0._db;           Mas(3,3) = 1._db
  endif

  read(8,*)
  read(8,'(6f10.7)') axyz(1:3), angxyz(1:3)
  axyz(1:3) = axyz(1:3) * bohr

  index = 0
  do jatom = 1,ntype - ntype_bulk

    index = index + 1
    read(8,'(5x,i3,1x,3(3x,f10.7))') its, posn(:,index)
    it = abs( its )
    itype(index) = it
    its_lapw(index) = its
    numprot(index) = index
    iprot = index

    do multi_run = 1,n_multi_run_e
      if( iabsorig(multi_run) == jatom ) iabsm(multi_run) = index
    end do

    do itr = 1,ntrans
      index = index + 1
      posn(1:3,index) = posn(1:3,iprot) + ptrans(1:3,itr)
      itype(index) = it
      its_lapw(index) = its
      numprot(index) = iprot
    end do

    read(8,'(15x,i2)') mult
    do mu = 1,mult-1
      index = index + 1
      read(8,'(5x,i3,1x,3(3x,f10.7))') ittt, posn(1:3,index)
      itype(index) = it
      its_lapw(index) = its
      numprot(index) = iprot
      do itr = 1,ntrans
        index = index + 1
        posn(1:3,index) = posn(1:3,index-itr) + ptrans(1:3,itr)
        itype(index) = it
        its_lapw(index) = its
        numprot(index) = iprot
      end do
    end do

! Maillage radial dans les spheres atomiques
    read(8,'(15x,i5,2(5x,f10.5),5x,f5.2)') nrato_lapw(it), r0_lapw(it), rlapw(it), rZ
    nrm = max( nrm, nrato_lapw(it) + 200 )

    numat(it) = nint( rZ )

! Matrice de rotation
    read(8,'(20x,3f10.7)') (rotloc(i,1:3,it), i = 1,3)

  end do

  where( posn > 1._db - eps10 ) posn = posn - 1._db

! Operations de symetrie
  read(8,'(i4)') nmatsym
  do is = 1,nmatsym
    read(8,'(3(3i2,f10.5,/))') (Wien_matsym(i,1:3,is), Wien_taulap(i,is), i = 1,3)
  end do
  close(8)

! Recherche de l'operation de symetrie qui renvoit a l'atome prototypique
  nt = ntrans + 1

  boucle_exter: do igr = 1,n_atom_per
    it = itype(igr)
    iprot = numprot(igr)

    if( igr == iprot ) then
      do i = 1,3
        rotloc_lapw(i,1:3,igr) = rotloc(1:3,i,it)
      end do
      cycle
    elseif( ntrans > 0 .and. mod(igr,nt) /= 1 ) then
      do i = 1,3
        rotloc_lapw(i,1:3,igr) = rotloc_lapw(i,1:3,igr-1)
      end do
      cycle
    endif

    do is = 1,nmatsym
      do j = 1,3
        p(j) = sum( Wien_matsym(j,:,is) * (posn(:,igr) - Wien_taulap(:,is)) )
        if( p(j) < -epspos ) then
          p(j) = p(j) + 1._db
        elseif( p(j) >= 1._db-epspos ) then
          p(j) = p(j) - 1._db
        endif
      end do

      dp(1:3) = abs( p(1:3) - posn(1:3,iprot) )
      if( dp(1) < epspos .and. dp(2) < epspos .and. dp(3) < epspos ) then
        Rot(:,:) = Wien_matsym(:,:,is)
        if( Trans == 'H' ) Rot = Matmul( Mat, Matmul( Rot, Mas ) )
        do i = 1,3
          do j = 1,3
            Rotloc_lapw(j,i,igr) = sum( Rot(i,1:3) * Rotloc(1:3,j,it) )
          end do
        end do
        cycle boucle_exter
      endif
    end do

    call write_error
    do ipr = 3,9,3
      if( icheck == 0 .and. ipr == 3 ) cycle
      write(ipr,110)
    end do
    stop

  end do boucle_exter

  if( icheck > 1 ) then
    do ia = 1,n_atom_per
      write(3,130) ia
      write(3,'(3f7.3)') (rotloc_lapw(i,1:3,ia), i = 1,3)
    end do
  endif

  return
  110 format(/' Atome symetrique non trouve !')
  130 format(/' rotloc(ia=',i2,')')
end

!***********************************************************************

! Calcul de la matrice de changement de repere maille - orthogonale

subroutine cal_cubmat(angxyz,Cubmat,Struct)

  use declarations
  implicit none

  character(len=5) struct

  logical ang(3), ange(3)

  real(kind=db):: a, alfa, b, beta, cosa, cosb, cosg, gamma, sina, sinb
  real(kind=db), dimension(3):: angxyz(3)
  real(kind=db), dimension(3,3):: cubmat

! Matrice de changement de repere cristallo, cubique
  ang(:) = abs( angxyz(:) - pi / 2 ) < eps4
  ange(1) = abs( angxyz(2) - angxyz(3) ) < eps4
  ange(2) = abs( angxyz(3) - angxyz(1) ) < eps4
  ange(3) = abs( angxyz(1) - angxyz(2) ) < eps4
  if( ange(1) .and. ange(2) .and. ange(3) ) then
    if( ang(1) ) then
      struct = 'cubic'
    else
      struct = 'trigo'
    endif
  elseif( ( abs( angxyz(3) - 2 * pi / 3 ) < eps4 ) .and. ang(1) .and. ang(2) ) then
    struct = 'hexag'
  else
    struct = 'autre'
  endif

  if( struct /= 'cubic' ) then

    if( struct == 'trigo' ) then
      alfa = angxyz(1)
      cosa = sqrt( ( 1._db + 2 * cos( alfa ) ) / 3._db )
      sina = sqrt( 1._db - cosa**2 )
      cubmat(1,1) = sina;  cubmat(1,2:3) = -0.5_db * sina
      cubmat(2,1) = 0._db;  cubmat(2,2) = sqrt(3._db) * sina / 2
      cubmat(2,3) = - cubmat(2,2)
      cubmat(3,1:3) = cosa
    else
      alfa = angxyz(1)
      beta = angxyz(2)
      gamma = angxyz(3)
      sina = sin( alfa )
      cosa = cos( alfa )
      sinb = sin( beta )
      cosb = cos( beta )
      cosg = cos( gamma )
      a = ( cosg - cosa*cosb ) / sinb
      b = sqrt( sina**2 - a**2 )
      cubmat(1,1) = sinb;  cubmat(1,2) = a;     cubmat(1,3) = 0._db
      cubmat(2,1) = 0._db; cubmat(2,2) = b;     cubmat(2,3) = 0._db
      cubmat(3,1) = cosb;  cubmat(3,2) = cosa;  cubmat(3,3) = 1._db
    endif

  else

     cubmat(:,:) = 0._db
     cubmat(1,1) = 1._db; cubmat(2,2) = 1._db; cubmat(3,3) = 1._db

  endif

  return
end

!***********************************************************************

! Calculation of the rotation matrix versus the Euler angles
! First rotation about z
! second rotation about y' (new y)  ( and not x' )
! third rotation about z' (new z)

subroutine Mat_Euler(Ang,Rot)

  use declarations
  implicit none

  real(kind=db):: c1, c2, c3, s1, s2, s3
  real(kind=db), dimension(3):: Ang
  real(kind=db), dimension(3,3):: Rot

  c1 = cos( Ang(1) )
  s1 = sin( Ang(1) )
  c2 = cos( Ang(2) )
  s2 = sin( Ang(2) )
  c3 = cos( Ang(3) )
  s3 = sin( Ang(3) )

  Rot(1,1) = c1*c2*c3 - s1*s3; Rot(1,2) = - c1*c2*s3 - s1*c3; Rot(1,3) = c1*s2
  Rot(2,1) = s1*c2*c3 + c1*s3; Rot(2,2) = - s1*c2*s3 + c1*c3; Rot(2,3) = s1*s2
  Rot(3,1) = - s2*c3;          Rot(3,2) = s2*s3;              Rot(3,3) = c2

  return
end

!***********************************************************************

! Calculation of the center of the cluster

subroutine Auto_center(axyz,Centre,Centre_auto_abs,Cubmat,icheck,itype,n_atom_uc,n_Z_abs,ngroup,ntype,numat, &
                       numat_abs,posn,Rsorte_s,Struct)

  use declarations
  implicit none

  integer:: i, ia1, ia2, ia3, ia4, icheck, igr, ipr, jgr, ngr, n_atom_uc, n_Z_abs, ngroup, ntype, Z
  integer, dimension(0:ntype):: numat
  integer, dimension(ngroup):: itype
  integer, dimension(n_Z_abs):: numat_abs

  character(len=2):: Chemical_Symbol
  character(len=5) struct

  logical:: Centre_auto_abs

  real(kind=db):: dist, dist_max, Radius, Rsorte_s
  real(kind=db), dimension(3):: axyz, b, Centre, p, q, v
  real(kind=db), dimension(3,3):: Cubmat, Mat, Mati
  real(kind=db), dimension(3,n_atom_uc):: posn, pos

  dist_max = 0._db

  jgr = 0
  do i = 1,n_Z_abs
    do igr = 1,n_atom_uc
      Z = numat( abs( itype(igr) ) )
      if( Centre_auto_abs .and. Z /= numat_abs(i) ) cycle
      jgr = jgr + 1
      if( struct /= 'cubic' ) then
        p(:) = posn(:,igr)
        p = matmul( Cubmat, p )
        pos(:,jgr) = p(:) * axyz(:)
      else
        pos(:,jgr) = posn(:,igr) * axyz(:)
      endif
    end do
  end do

  ngr = jgr

  do igr = 1,ngr
    Z = numat( abs( itype(igr) ) )
    if( igr == 1 ) Centre(:) = pos(:,igr)
    do jgr = igr+1,ngr
      p(:) = pos(:,igr) - pos(:,jgr)
      dist = sqrt( sum( ( p(:) )**2 ) )
      if( dist < dist_max ) cycle
      dist_max = dist
      Centre(:) = 0.5_db * ( pos(:,igr) + pos(:,jgr) )
      ia1 = igr
      ia2 = jgr
    end do
  end do

  Radius = dist_max / 2

  dist_max = Radius
  ia3 = 0
  do igr = 1,ngr
    p(:) = pos(:,igr) - Centre(:)
    dist = sqrt( sum( ( p(:) )**2 ) )
    if( dist < dist_max + eps10 ) cycle
    dist_max = dist
    ia3 = igr
  end do

  if( ia3 /= 0 ) then
! Recherche du cercle circonscrit

! Plan hauteur 1
    v(:) = pos(:,ia2) - pos(:,ia1)
    Mat(1,:) = v(:)
    b(1) = sum( v(:)*Centre(:) )

! Plan hauteur 2
    v(:) = pos(:,ia3) - pos(:,ia1)
    Mat(2,:) = v(:)
    b(2) = 0.5_db * sum( v(:) * ( pos(:,ia3) + pos(:,ia1) ) )

! Plan du triangle
    p(:) = pos(:,ia2) - pos(:,ia1)
    q(:) = pos(:,ia3) - pos(:,ia1)
    call prodvec(v,p,q)
    Mat(3,:) = v(:)
    b(3) = sum( v(:)*Centre(:) )

    call invermat(Mat,Mati)

    Centre = Matmul( Mati, b )

    v(:) = pos(:,ia3) - Centre(:)
    Radius = sqrt( sum( v(:)**2 ) )

    dist_max = Radius
    ia4 = 0
    do igr = 1,ngr
      p(:) = pos(:,igr) - Centre(:)
      dist = sqrt( sum( ( p(:) )**2 ) )
      if( dist < dist_max + eps10 ) cycle
      dist_max = dist
      ia4 = igr
    end do

    if( ia4 /= 0 ) then
! Recherche de la sphere circonscrite a 4 points

! Plan hauteur 3
      v(:) = pos(:,ia4) - pos(:,ia1)
      Mat(3,:) = v(:)
      b(3) = 0.5_db * sum( v(:) * ( pos(:,ia4) + pos(:,ia1) ) )

      call invermat(Mat,Mati)

      Centre = Matmul( Mati, b )

      v(:) = pos(:,ia4) - Centre(:)
      Radius = sqrt( sum( v(:)**2 ) )

    endif

  endif

  call invermat(Cubmat,Mati)

  Centre = Matmul( Mati, Centre )
  Centre(:) = Centre(:) / axyz(:)

  if( icheck > 0 ) then
    write(3,110) Centre(:)
    if( Centre_auto_abs ) then
      if( n_Z_abs == 1 ) then
         write(3,120) Chemical_Symbol(numat_abs(1)), Radius * bohr
      elseif( n_Z_abs == 2 ) then
         write(3,130) Chemical_Symbol(numat_abs(1)), Chemical_Symbol(numat_abs(2)), Radius * bohr
      else
        write(3,140) ( Chemical_Symbol(numat_abs(i)), i = 1,n_Z_abs )
        write(3,150) Radius * bohr
      endif
    else
      write(3,160) Radius * bohr
    endif
  endif

  if( Centre_auto_abs .and. Radius > Rsorte_s ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,170) Radius*bohr, Rsorte_s*bohr
    end do
    stop
  endif

  return
  110 format(/' Center set at  =',3f9.5,'  (cell unit)')
  120 format(' Farthest ',a2,' atom at =',f7.3,' A, from the center of the cluster')
  130 format(' Farthest ',a2,' or ',a2,' atom at =',f7.3,' A, from the center of the cluster')
  140 format(' Farthest atom',10(1x,a2))
  150 format(10x,'at =',f9.5,' A, from the center of the cluster')
  160 format(' Farthest atom at =',f7.3,' A, from the center of the cluster')
  170 format(//' The farther absorbing atom is at R =',f7.3,' A  > Cluster radius =',f7.3,' A',// &
               ' Increase the cluster radius !',/)
end
