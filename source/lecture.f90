! FDMNES subroutines

! Sousprogramme de lecture des fichiers d'entree permettant d'etablir
! les dimensions des differents tableaux

subroutine lectdim(Absauto,Atom_occ_hubb,Atom_nonsph,Axe_loc,Bormann,Doping,Extract,Flapw,Full_self_abs,Hubbard,itape4, &
    Magnetic,Memory_save,mpinodes0,mpirank0,n_file_dafs_exp,n_multi_run_e,nb_atom_conf_m,ncolm,neimagent,nenerg,ngamme,ngroup, &
    ngroup_neq,nhybm,nklapw,nlatm,nlmlapwm,nmatsym,norbdil,npldafs,nple,nplrm,n_adimp,n_radius,n_range,nq_nrixs,NRIXS,nspin, &
    nspino,nspinp,ntype,ntype_conf,Pdb,Readfast,Self_abs,Space_file,Taux,Temperature,Use_FDMX,Xan_atom)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: eof, eoff, i, ie, ier, igamme, igr, igrdat, io, ipr, ipl, istat, it, itape4, itype_dop, j,&
    jgr, jpl, k, kgr, l, ligne, lin_gam, mpinodes0, mpirank0, n, n_adimp, n_radius, n_range, n_dic, n_file_dafs_exp, &
    n_fract_x, n_fract_y, n_fract_z, n_multi_run_e, n_skip, na, nb, nb_atom_conf_m, ncolm, neimagent, &
    nenerg, ngamme, ngc, ngroup, ngroup_neq, nhybm, nklapw, nl, &
    nlatm, nlmlapwm, nmatsym, nn, nnombre, norbdil, norbv, mpierr, &
    npldafs, nple, nplm, nplrm, nq_nrixs, nspin, nspino, nspinp, ntype, ntype_conf, Wien_save, Z

  character(len=2):: Chemical_Symbol, Chemical_Symbol_c, Symbol
  character(len=4):: mot4
  character(len=6):: mot6
  character(len=9):: grdat
  character(len=13):: Space_Group, Spgr
  character(len=132):: Fichier, Fichier_cif, Fichier_pdb, identmot, mot, motsb, Space_file
  character(len=132), dimension(9):: Wien_file

  integer, dimension(:), allocatable :: igra, neq, numat

  logical:: Absauto, adimpin, Atom_conf, Atom_nonsph, Atom_occ_hubb, Axe_loc, Bormann, Cif, Cylindre, Doping, &
     Extract, Flapw, Full_self_abs, Hubbard, Magnetic, Matper, Memory_save, NRIXS, Pdb, Quadrupole, Readfast, &
     Screening, Self_abs, Spherique, Taux, Temperature, Use_FDMX, Xan_atom

  real(kind=db):: Adimp, Angz, de, def, number_from_text, E, phi, q_nrixs_first, q_nrixs, q_nrixs_step, q_nrixs_last, r, &
                  Rsorte_s, Theta
  real(kind=db), dimension(3):: p
  real(kind=db), dimension(3,3):: Mat
  real(kind=db), dimension(:), allocatable:: E_adimp, E_radius, Egamme
  real(kind=db), dimension(:,:), allocatable:: pop, posn, posout

  Absauto = .true.
  Angz = 0._db
  Atom_conf = .false.
  Atom_nonsph = .false.
  Atom_occ_hubb = .false.
  Axe_loc = .false.
  Cif = .false.
  Cylindre = .false.
  Doping = .false.
  Extract = .false.
  Flapw = .false.
  Full_self_abs = .false.
  Hubbard = .false.
  Matper = .false.
  Memory_save = .false.
  n_adimp = 1
  n_dic = 0
  n_file_dafs_exp = 0
  n_multi_run_e = 1
  n_range = 1
  n_radius = 1
  nb_atom_conf_m = 0
  neimagent = 0
  nenerg = 131
  ngamme = 3
  nhybm = 1
  nklapw = 1
  nlmlapwm = 1
  nlatm = 0
  nmatsym = 1
  norbdil = 0
  nple = 0
  nq_nrixs = 0
  NRIXS = .false.
  nspin = 1
  nspino = 1
  nspinp = 1 ! pour la diffusion ( = 2 si spinorbite meme potentiel non magnetique )
  ntype = 0
  ntype_conf = 0
  Pdb = .false.
  Quadrupole = .false.
  Readfast = .false.
  Screening = .false.
  Self_abs = .false.
  Space_Group = ' '
  Spherique = .false.
  Taux = .false.
  Temperature = .false.
  Xan_atom = .false.
  Wien_save = 0

  adimpin = .false.
  if( Use_FDMX ) n_adimp = 5

  if( Bormann ) then
    npldafs = 36
  else
    npldafs = 0
  endif

  if( mpirank0 == 0 ) then

    Rewind(itape4)

    boucle_gen: do

      read(itape4,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit boucle_gen

      grdat = identmot(mot,9)

      if( grdat(1:1) == '!' ) cycle

      select case(grdat)

        case('doping')
          Doping = .true.
          read(itape4,*) itype_dop
          absauto = .false.
          n_multi_run_e = 1

        case('absorbeur')
          absauto = .false.
          n = nnombre(itape4,132)
          if( n > 0 ) n_multi_run_e = 0
          do i = 1,100000
            n = nnombre(itape4,132)
            if( n < 1 ) exit
            n_multi_run_e = n_multi_run_e + n
            read(itape4,*)
          end do

        case('end')
          exit

        case('extract')
          Extract = .true.

        case('hubbard')
          Hubbard = .true.

        case('nrixs')
          NRIXS = .true.
          n = nnombre(itape4,132)
          if( n /= 1 .and. n /= 3 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,'(/A/)') ' After keyword "NRIXS" one must have 1 or 3 numbers : q_first, q_step, q_last !'
            end do
            stop
          endif
          if( n == 1 ) then
            read(itape4,*,iostat=ier) q_nrixs_first
            nq_nrixs = 1
          else
            read(itape4,*,iostat=ier) q_nrixs_first, q_nrixs_step, q_nrixs_last
          endif
          if( ier > 0 ) call write_err_form(itape4,grdat)
          if( n == 3 ) then
            q_nrixs = q_nrixs_first
            do nq_nrixs = 2,100000
              q_nrixs = q_nrixs + q_nrixs_step
              if( q_nrixs > q_nrixs_last + eps10 ) exit
            end do
            nq_nrixs = nq_nrixs - 1
          endif

        case('spgroup')
          n = nnombre(itape4,132)
          read(itape4,'(A)') mot
          if( mot(1:1) == ' ' ) mot = adjustl(mot)
          Space_group = mot(1:13)

        case('temperatu')
          Temperature = .true.

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

        case('range','rangel')
          ngamme = nnombre(itape4,132)
          if( mod(ngamme,2) == 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,120) ngamme
            end do
            stop
          endif
          if( grdat == 'rangel' ) then
            lin_gam = 1
            if( ngamme /= 1 ) ngamme = 3
          else
            lin_gam = 0
          endif
          allocate( egamme(ngamme) )
          read(itape4,*,iostat=ier) egamme(1:ngamme)
          if( ier > 0 ) call write_err_form(itape4,grdat)

          do igamme = 2,ngamme,2
            if( egamme(igamme) > eps6 ) cycle
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,130)
            end do
            stop
          end do

          E = egamme(1)
          if( ngamme == 1 ) then
            nenerg = 1
          elseif( egamme(3) <= egamme(1) ) then
            nenerg = 1
          elseif( lin_gam == 1 ) then
            def = 10 / rydb
            do ie = 2,10000000
              r = 1 + ( E/rydb ) / def
              r = max( r, 0.25_db )
              de = sqrt( r ) * egamme(2)
              e = e + de
              if( e > egamme(ngamme) + eps10 ) then
                nenerg = ie - 1
                exit
              endif
            end do
          else
            ngc = 2
            do ie = 2,1000000
              e = e + egamme(ngc)
              if( e > egamme(ngamme) + eps10 ) then
                nenerg = ie - 1
                exit
              elseif( e >= egamme(ngc+1) - eps10 ) then
                if( ngc+1 == ngamme ) then
                  nenerg = ie
                  exit
                endif
                if( egamme(ngc+3) <= egamme(ngc+1) ) then
                  nenerg = ie
                  exit
                endif
                ngc = ngc + 2
              endif
            end do
          endif
          deallocate( egamme )

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
            if( ier > 0 ) call write_err_form(itape4,grdat)
            if( sum( p(:) )**2 < eps10 ) n_dic = n_dic + 1
          end do
          nple = jpl - 1

        case('dafs')
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
                  case(6)
                    read(itape4,*)
                    read(itape4,*)
                  case default
                    call write_error
                    do ipr = 6,9,3
                      write(ipr,140) ipl
                    end do
                    stop
                end select
              case(5,6,7,8)
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

        case('dafs_exp')
          do i = 1,3
            read(itape4,*)
          end do
          npldafs = 0
          do i = 1,1000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            read(itape4,'(/A)') Fichier
            Fichier = Adjustl(Fichier)
            l = len_trim(Fichier)
            if( l > 4 ) then
              if( Fichier(l-3:l-3) /= '.' ) Fichier(l+1:l+4) = '.txt'
            endif
            open(99, file=Fichier, status='old',iostat=istat)
            if( istat /= 0 ) call write_open_error(Fichier,istat,1)
            n = nnombre(99,100000)
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
            read(itape4,*) nb
            nb_atom_conf_m = max( nb, nb_atom_conf_m )
            backspace(itape4)
            read(itape4,*) ( nl, i = 1,nb+2 )
            nlatm = max( nlatm, nl )
            if( n == nb + 2 + 4*nl .and. nl > 0 .and. nspin == 1 ) then
              backspace(itape4)
              allocate( pop(nl,2) )
              read(itape4,*) ( i, l = 1,nb+1 ), nl, ( i, i, pop(l,:), l = 1,nl )
              do l = 1,nl
                if( abs( pop(l,1) - pop(l,2) ) < eps10 ) cycle
                nspin = 2
                exit
              end do
              deallocate( pop )
            endif
          end do
          ntype = ntype_conf

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

! Description de l'agregat :
        case('crystal','molecule','crystal_t','molecule_')
          if( grdat == 'crystal_t' .or. grdat == 'molecule_' ) Taux = .true.

          Matper = grdat(1:5) == 'cryst'

          n = nnombre(itape4,132)
          if( n == 0 ) then
            call write_error
            read(itape4,'(A)') motsb
            do ipr = 6,9,3
              write(ipr,150) grdat
              write(ipr,'(A)') motsb
              write(ipr,160)
            end do
            stop
          endif
          if( n > 5 ) then
            read(itape4,*,iostat=ier) ( Angz, i = 1,6 )
          else
            read(itape4,*)
          endif

          if( Readfast .or. Taux ) then
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
                  if( ier > 0 ) call write_err_form(itape4,grdat)
                  if( norbv < 0 ) then
                    nhybm = max( nhybm, - norbv - 1 )
                  else
                    nhybm = max( nhybm, norbv )
                  endif
                case default
                  call write_error
                  read(itape4,'(A)') motsb
                  do ipr = 6,9,3
                    write(ipr,150) grdat
                    write(ipr,'(A)') motsb
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
          ngroup = igr - 1

! Description de l'agregat :
        case('cif_file')

          Cif = .true.
          Matper = .true.

          Fichier_cif = ' '
          read(itape4,'(A)') Fichier_cif
          Fichier_cif = Adjustl(Fichier_cif)
          l = len_trim(Fichier_cif)
          if( l > 4 ) then
            if( Fichier_cif(l-3:l-3) /= '.' ) Fichier_cif(l+1:l+4) = '.cif'
          endif
          open(8, file = Fichier_cif, status='old', iostat=istat)
          if( istat /= 0 ) call write_open_error(Fichier_cif,istat,1)

          do

            read(8,'(A)',iostat=eoff) mot

            if( eoff /= 0 ) exit

            if( mot(1:30) == '_symmetry_space_group_name_H-M' ) then

              l = len_trim(mot)

              do i = 31,l
                if( mot(i:i) == '''' ) exit
              end do
              k = 0
              j = i + 1
              do i = j,l
               if( mot(i:i) == ' ' ) cycle
               if( mot(i:i) == '''' ) exit
               k = k + 1
               Space_group(k:k) = mot(i:i)
              end do

            elseif( mot(1:17) == '_cell_angle_gamma' ) then
              Angz = number_from_text(1,mot)

            elseif( mot(1:16) == '_atom_site_label') then

              n_fract_x = 0; n_fract_y = 0; n_fract_z = 0

              do n = 1,100000
                read(8,'(A)') mot
                if( mot(1:1) /= '_' ) then
                  backspace(8)
                  exit
                elseif( mot(2:18) == 'atom_site_fract_x') then
                  n_fract_x = n
                elseif( mot(2:18) == 'atom_site_fract_y') then
                  n_fract_y = n
                elseif( mot(2:18) == 'atom_site_fract_z') then
                  n_fract_z = n
                elseif( mot(2:20) == 'atom_site_occupancy') then
                  Taux = .true.
                endif
              end do

              if( n_fract_x == 0 .or. n_fract_y == 0 .or. n_fract_z == 0 ) then
                call write_error
                do ipr = 6,9,3
                  write(ipr,110)
                  write(ipr,165) Fichier_cif
                end do
                stop
              endif

              igr = 0
              do
                read(8,'(a4)',iostat=eoff) mot4
                if( eoff /= 0 .or. mot4(1:1) == '_' .or. mot4 == 'loop' .or. mot4 == ' ' ) exit
                igr = igr + 1
              end do
              exit

            endif

          end do

          Close(8)

          ngroup = igr

          if( Space_group == ' ' ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,166) Fichier_cif
            end do
            stop
          endif

          if( ngroup == 0 ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,110)
              write(ipr,167) Fichier_cif
            end do
            stop
          endif

        case('pdb_file')

          Pdb = .true.
          Taux = .true.
          Matper = .true.

          Fichier_pdb = ' '
          read(itape4,'(A)') Fichier_pdb
          Fichier_pdb = Adjustl(Fichier_pdb)
          l = len_trim(Fichier_pdb)
          if( l > 4 ) then
            if( Fichier_pdb(l-3:l-3) /= '.' ) Fichier_pdb(l+1:l+4) = '.pdb'
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

                l = len_trim(Spgr)
                j = 0
                do i = 1,l
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
          ngroup = igr

        case('flapw','flapw_s','flapw_r','flapw_s_p','flapw_psi','flapw_n','flapw_n_p','flapw_s_n')

          Flapw = .true.

          if( grdat(6:7) == '_s' ) then
            n = nnombre(itape4,132)
            read(itape4,*)
          elseif( grdat(6:7) == '_r' ) then
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

          call lectpot_dim(ngroup,nklapw,nlmlapwm,nmatsym,Wien_file(1),Wien_file(2),ntype,Wien_save)

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

    if( Extract ) xan_atom = .false.
    if( Flapw .or. Extract ) Hubbard = .false.
    if( Flapw ) nspin = nspinp
    if( nspin == 2  ) then
      Magnetic = .true.
    else
      Magnetic = .false.
    endif
    if( npldafs == 0 ) Self_abs = .false.
    if( npldafs == 0 ) Full_self_abs = .false.
    if( Self_abs ) n_dic = n_dic + 2 * npldafs
    if( Full_self_abs ) n_dic = n_dic + 4 * npldafs
    if( .not. Matper ) Space_Group = ' '
    if( Pdb ) Temperature = .true.

    ngroup_neq = ngroup

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

    elseif( n_adimp > 1 .and. .not. Extract ) then

      n_range = n_adimp
      deallocate( E_adimp )

    elseif( n_radius > 1 .and. .not. Extract ) then

      n_range = n_radius
      deallocate( E_radius )

    else

      n_range = 1

    endif

    if( Space_Group /= ' ' .or. ntype == 0 .or. Atom_conf ) then

      allocate( neq(ngroup_neq) )
      allocate( posn(3,ngroup_neq) )
      allocate( posout(3,ngroup) )
      allocate( numat(ngroup) )

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
              if( igr == ngroup_neq ) exit

            case default

              cycle

          end select

        end do

        Close(8)

      elseif( Cif ) then

        open(8, file = Fichier_cif, status='old', iostat=istat)

        igr = 0

        do
          read(8,'(A)') mot
          if( mot(1:16) == '_atom_site_label' ) exit
        end do
        backspace(8)

        n = 0
        do
          n = n + 1
          read(8,'(A)') mot
          if( mot(1:1) /= '_' ) then
            backspace(8)
            exit
          elseif( mot(1:18) == '_atom_site_fract_x' ) then
            n_skip = n - 1
          elseif( mot(1:20) == '_atom_site_occupancy ') then
            Taux = .true.
          endif
        end do

        do igr = 1,ngroup_neq

          read(8,'(A)') mot

! Lecture element chimique
          Symbol = mot(1:2)
          if( Symbol(2:2) == '1' .or. Symbol(2:2) == '2' .or. Symbol(2:2) == '3' .or. Symbol(2:2) == '4' .or.  &
              Symbol(2:2) == '5' .or. Symbol(2:2) == '6' .or. Symbol(2:2) == '7' .or. Symbol(2:2) == '8' .or.  &
              Symbol(2:2) == '9' ) Symbol(2:2) = ' '

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

          posn(1,igr) = number_from_text(n_skip,mot)
          posn(2,igr) = number_from_text(n_skip+1,mot)
          posn(3,igr) = number_from_text(n_skip+2,mot)

        end do

        Close(8)

      else

        Rewind(itape4)

        do igrdat = 1,100000
          read(itape4,'(A)') mot
          grdat = identmot(mot,9)
          if( grdat(1:7) /= 'crystal' .and. grdat(1:7) /= 'molecul' ) cycle

          n = nnombre(itape4,132)
          if( n == 2 ) then
            cylindre = .true.
          elseif( n == 1 ) then
            spherique = .true.
          endif
          read(itape4,*)
          do igr = 1,ngroup_neq
            if( Readfast .or. Taux .or. .not. ( Atom_nonsph .or. Atom_occ_hubb .or. Axe_loc) ) then
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
                  if( ier > 0 ) call write_err_form(itape4,grdat)
              end select
              if( norbv == 0 ) cycle
              do io = 1,abs(norbv)
                read(itape4,*)
              end do
            endif
          end do

          exit
        end do

      endif

      if( ntype == 0 ) then

        ntype = 1
        boucle_1: do igr = 2,ngroup_neq
          do jgr = 1,igr-1
            if( numat(igr) == numat(jgr) ) cycle boucle_1
          end do
          ntype = ntype + 1
        end do boucle_1
        if( Doping ) then
          do igr = 1,ngroup_neq
            if( numat(igr) == itype_dop ) exit
          end do
          if( igr > ngroup_neq ) ntype = ntype + 1
        endif

      elseif( Atom_conf ) then

        allocate( igra(ngroup_neq) )

        Rewind(itape4)

        do igrdat = 1,100000
          read(itape4,'(A)') mot
          grdat = identmot(mot,9)
          if( grdat /= 'atom_conf' ) cycle

          na = 0
          do it = 1,ntype
            read(itape4,*) n, igra(na+1:na+n)
            na = na + n
          end do
          exit
        end do

        boucle_2: do igr = 1,ngroup_neq
          do kgr = 1,na
            if( igr == igra(kgr) ) cycle boucle_2
          end do
          boucle_3: do jgr = 1,igr-1
            do kgr = 1,na
              if( jgr == igra(kgr) ) cycle boucle_3
            end do
            if( numat(igr) == numat(jgr) ) cycle boucle_2
          end do boucle_3
          ntype = ntype + 1
        end do boucle_2

        deallocate( igra )

      endif

      if( Space_Group /= ' ' ) then

! Atom defined with few digits
        if( abs( Angz - 120._db ) < 0.0001_db ) then
          do i = 1,11
            if( i == 3 .or. i == 6 .or. i == 9 ) cycle
            Where( abs( Posn - i / 12._db ) < 0.00005_db ) Posn = i / 12._db
          end do
        endif

        do igr = 1,ngroup_neq
          if( cylindre ) then
            r = posn(1,igr)
            Theta = pi * posn(2,igr) / 180
            posn(1,igr) = r * cos( theta )
            posn(2,igr) = r * sin( theta )
          elseif( spherique ) then
            r = posn(1,igr)
            theta = pi * posn(2,igr) / 180
            phi = pi * posn(3,igr) / 180
            posn(1,igr) = r * sin( theta ) * cos( phi)
            posn(2,igr) = r * sin( theta ) * sin( phi)
            posn(3,igr) = r * cos( theta )
          endif
        end do

        call spgroup(.false.,neq,ngroup,ngroup_neq,numat,posn,posout,Space_file,space_group)

      endif

      deallocate( posn )
      deallocate( posout )
      deallocate( neq )
      deallocate( numat )

    endif

! Pour le cas des atomes charges ou il faut ajouter une orbitale:
    if( nlatm > 0 .or. Screening ) nlatm = nlatm + 1

  endif   ! arrive mpirank0 /= 0

  if( mpinodes0 > 1 ) then
    call MPI_Bcast(Absauto,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Atom_nonsph,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Atom_occ_hubb,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Axe_loc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Doping,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Memory_save,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Extract,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Flapw,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Full_self_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Hubbard,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Magnetic,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_adimp,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_dic,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_multi_run_e,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_radius,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_range,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nb_atom_conf_m,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(neimagent,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nenerg,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ngamme,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ngroup,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ngroup_neq,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nhybm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nklapw,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nlatm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nlmlapwm,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nmatsym,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(norbdil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(npldafs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nple,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nq_nrixs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nrixs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nspin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nspino,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nspinp,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ntype,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ntype_conf,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Pdb,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Quadrupole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Self_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Taux,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Temperature,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Xan_atom,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Wien_save,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
  endif

  if( nple == 0 ) then
    if( Quadrupole ) then
      nplm = 6
    else
      nplm = 3
    endif
  else
    nplm = nple
  endif

  nplrm = nplm + n_dic
  ncolm = nplm + 2*npldafs + 2*n_dic + 1
  if( xan_atom ) ncolm = ncolm + 1
  if( self_abs ) ncolm = ncolm  + 2*npldafs
  if( Full_self_abs ) ncolm = ncolm  + 4*npldafs

  if( Doping ) ngroup = ngroup + 1

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
  120 format(//' After keyword "Range", the number of value given for the energies and steps',/ &
               ' is even (',i2,'), it must be odd !',/ &
               ' If it is not the case, check the presence of extra characters or tabulations.',/ &
               ' They are forbidden !'//)
  130 format(//' Energy step is zero or negative after keyword "Range", forbidden !'//)
  140 format(//' After the keyword "Dafs", for the reflection number', i3,',',/ &
               ' the polarization or the indexes are not well set.',/ &
               ' Check the format !'//)
  150 format(///' Just after the keyword "',A,'"', ', the following line is read :',/)
  160 format(/' It must not be there or it contains unwanted characters !'//)
  165 format(///' Error in the lecture of the cif_file :',A,/ &
                '   The position of the atoms are not found',/ &
                '   (with tag _atom_site_fract_x, _atom_site_fract_y and _atom_site_fract_z) !' //)
  166 format(///' Error in the lecture of the cif_file :',A,/ &
                '   No space group found ( with tag _symmetry_space_group_name_H-M ) !' //)
  167 format(///' Error in the lecture of the cif_file :',A,/ &
                '   No list of atom position found !' //)
  170 format(/' Error in the Pdb file :',// &
              ' The chemical symbol of atom number',i6,' is: ', a2,', what is not known!'//)
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
  110 format(//' Format error when reading in the indata file under', ' the keyword:',//,5x,A)
  120 format(//' Check :',/ '  - How many numbers must be in the line ?',/ '  - Are there spaces between the numbers ?',/ &
           '  - Are there unwanted characters, extra  points... ?',/ '  - Tabulations are forbidden !'//)
end

!***********************************************************************

! Routine de lecture des potentiels et densites electroniques venant de WIEN

subroutine lectpot_dim(ngroup,nklapw,nlmlapwm,nmatsym,nomstruct,nomvcoul,ntype,Wien_save)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  character(len=1) Trans
  character(len=132) nomstruct, nomvcoul

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

  ngroup = index

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

    end do ! fin de la boucle sur les atomes

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

!   Sousprogramme de lecture.
! Toutes les entrees sauf les densites electroniques et les fonctions d'onde sont en Angstroem et en eV.
! Elles sont converties, pour tout le programme, en unites atomiques et Rydberg dans ce sous-programme.

subroutine lecture(Absauto,adimp,alfpot,All_nrixs,Allsite,Ang_borm,Ang_rotsup,Angle_or,Angpoldafs,Angxyz,ATA,Atom_occ_hubb, &
    Atom_nonsph,Atom_nsph_e,Atomic_scr,Axe_atom_gr,Axe_loc,axyz,Base_spin,basereel,Bormann,Cartesian_tensor,Charge_free, &
    Clementi,com,comt,Core_resolved,Coupelapw,Cubmat,D_max_pot,Dafs,Dafs_bio,Delta_En_conv,Delta_Epsii,Density,Density_comp, &
    Dipmag,Doping,dpos,dyn_eg,dyn_g,E_adimp,E_radius,E_max_range,Eclie,Eclie_out,Ecrantage,Eeient,Egamme,Eimagent, &
    Eneg_i,Eneg_n_i,Energphot,Extract,f_no_res,FDM_comp,Fit_cal,Flapw,Flapw_new, &
    Force_ecr,Full_atom_e,Full_potential,Full_self_abs,Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,Green_int,Green_s, &
    Green_self,hkl_borm,hkl_dafs,Hubb,Hubbard,Hybrid,iabsm,iabsorig,icheck,icom,igr_dop,indice_par,iscratch, &
    isigpi,itdil,its_lapw,iord,itape4,itype,itype_dop,jseuil,Kern_fac,Kgroup,korigimp,lmax_nrixs,l_selec_max, &
    lamstdens,ldil,lecrantage,lin_gam,lmax_pot,lmaxfree,lmaxso0,lmaxat0,lmoins1,lplus1,lseuil,lvval,m_hubb_e, &
    Magnetic,Mat_or,Matper,MPI_host_num_for_mumps,mpinodes,mpinodes0,mpirank0,Muffintin, &
    Multipole,multrmax,n_adimp,n_atom_proto,n_devide,n_file_dafs_exp,n_multi_run_e,n_radius,n_range, &
    nb_atom_conf_m,nbseuil,nchemin,necrantage,neimagent,nenerg,ngamh,ngamme,ngroup,ngroup_hubb,ngroup_lapw,ngroup_m, &
    ngroup_neq,ngroup_nonsph,ngroup_par,ngroup_pdb,ngroup_taux, &
    ngroup_temp,nhybm,nlat,nlatm,nnlm,No_solsing,nom_fich_Extract, &
    nomfich,nomfich_optic_data,nomfich_tddft_data,nomfichbav,Noncentre, &
    Nonexc,norbdil,norbv,Normaltau,normrmt,npar,nparm,nphi_dafs, &
    nphim,npldafs,nple,nposextract,nq_nrixs,nrato,nrato_dirac,nrato_lapw,nrm, &
    nself,nseuil,nslapwm,nspin,nsymextract,ntype,ntype_conf,numat,numat_abs, &
    nvval,occ_hubb_e,Octupole,Old_reference,One_run,Optic,Overad,Overlap,p_self_max,p_self0, &
    param,Pas_SCF,pdpolar,PointGroup,PointGroup_Auto,polar,Polarise,poldafsem,poldafssm, &
    pop_nonsph,popats,popval,posn,q_nrixs,Quadmag,Quadrupole,R_rydb, &
    r0_lapw,rchimp,Readfast,Recup_optic,Recup_tddft_data,Relativiste,r_self,rlapw,rmt,rmtimp,Rot_Atom_gr,rotloc_lapw, &
    roverad,RPALF,rpotmax,rydberg,rsorte_s,Save_optic,Save_tddft_data,SCF_log,Self_abs, &
    Solsing_s,Solsing_only,Solver,Space_file,Spherical_signal,Spherical_tensor, &
    Spinorbite,state_all,state_all_out,Struct,Supermuf,Symauto,Symmol,Taux,Taux_oc,Tddft,Temp,Temp_coef, &
    Temperature,Tensor_imp,Test_dist_min,Trace_format_wien,Trace_k,Trace_p,Typepar,Use_FDMX,V_hubbard,V_intmax,Vec_orig, &
    Vecdafsem,Vecdafssm,Veconde,V0bdcFimp,Wien_file,Wien_matsym,Wien_save,Wien_taulap,Ylm_comp_inp,Z_nospinorbite)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: eof, eoff, i, i_range, ier, ia, ie, igr, igr_dop, igrdat, io, iord, ip, ipar, ipl, ipl0, ipr, ipr0, iscratch, isp, &
    ispin, istat, istop, isymeq, it, itape4, itype_dop, j, jgr, jpl, jseuil, jt, k, kgr, l, l_hubbard, l_level_val, &
    l_selec_max, l1, l2, lamstdens, lecrantage, lin_gam, lmax_nrixs, lmax_pot, lmax_pot_default, lmaxat0, lmaxso0, long, &
    lseuil, m,  m_hubb_e, MPI_host_num_for_mumps, mpierr, mpinodes, mpinodes0, mpirank0, &
    multi_run, multrmax, n, n_adimp, n_atom_proto, n_devide, n_file_dafs_exp, n_fract_x, n_fract_y, n_fract_z, n_multi_run_e, &
    n_occupancy, n_orb_rel, n_radius, n_range, n1, n2, natomsym, nb_atom_conf_m, nbseuil, nchemin, &
    necrantage, neimagent, nenerg, ngamme, ngamh, ngroup, ngroup_hubb, ngroup_lapw, ngroup_m, ngroup_neq, ngroup_nonsph, &
    ngroup_par, ngroup_pdb, ngroup_taux, ngroup_temp, nhybm, nlatm, nn, nnlm, nnombre, non_relat, norb, norbdil, normrmt, &
    nparm, nphim, nphimt, npldafs, nple, nq_nrixs, nrato_dirac, nrm, nscan, nself, nseuil, nslapwm, nspin, ntype, ntype_conf, &
    numat_abs, Trace_k, Wien_save, Z_nospinorbite, Z

  character(len=1):: Let
  character(len=2):: Chemical_Symbol, Chemical_Symbol_c, Symbol
  character(len=3):: Seuil
  character(len=5):: Solver, Struct
  character(len=6):: mot6
  character(len=8):: PointGroup
  character(len=9):: grdat
  character(len=11):: motpol
  character(len=13):: Chemical_Name, mot13, Space_Group, Spgr
  character(len=132):: comt, Fichier, identmot, mot, motsb, nomfich, nomfich_optic_data, nomfich_tddft_data, &
    nom_fich_Extract, nomfichbav, Space_file
  character(len=35), dimension(0:ntype):: com
  character(len=132), dimension(9):: Wien_file
  character(len=132), dimension(n_file_dafs_exp):: File_dafs_exp
  character(len=9), dimension(ngroup_par,nparm):: typepar

  integer, dimension(3):: hkl_borm, ldipimp
  integer, dimension(12):: Tensor_imp
  integer, dimension(30):: icheck
  integer, dimension(3,3):: lquaimp
  integer, dimension(n_multi_run_e):: iabsm, iabsorig, nposextract, nsymextract
  integer, dimension(ngroup_par):: npar
  integer, dimension(ngroup_lapw):: its_lapw
  integer, dimension(norbdil):: itdil, ldil
  integer, dimension(0:ntype):: icom, nlat, nrato, nrato_lapw, numat
  integer, dimension(ngroup):: itype
  integer, dimension(ngroup_pdb):: Kgroup

  integer, dimension(0:ngroup_nonsph):: norbv
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(3,npldafs):: hkl_dafs
  integer, dimension(npldafs,2):: isigpi
  integer, dimension(ntype_conf):: nb_atom_conf
  integer, dimension(nb_atom_conf_m,ntype_conf):: igra
  integer, dimension(ngroup_par,nparm):: indice_par
  integer, dimension(0:ntype,nlatm):: lvval, nvval
  integer, dimension(3,3,nslapwm):: Wien_matsym

  integer, dimension(:), allocatable :: itZ, neq

  complex(kind=db), dimension(3,npldafs):: poldafsem, poldafssm
  complex(kind=db), dimension(nhybm,16,ngroup_nonsph):: Hybrid

  logical:: Absauto, adimpin, All_nrixs, Allsite, ATA, Atom, Atom_conf, Atom_nonsph, Atom_occ_hubb, Atomic_scr, Axe_loc, &
    Base_spin, Basereel, Bormann, Cartesian_tensor, Charge_free, Centre_auto, Centre_auto_abs, Clementi, Core_resolved, &
    Core_resolved_e, Coupelapw, Cylindre, Dafs, Dafs_bio, Density, Density_comp, Dipmag, Doping, dyn_eg, dyn_g, &
    E1E1, E1E2e, E1E3, E1M1, E1M2, E2E2, E3E3, eneg_i, eneg_n_i, Energphot, &
    exc_imp, Extract, FDM_comp, Fermi_auto, Fit_cal, Flapw, Flapw_new, Force_ecr, Full_atom_e, Full_potential, Full_self_abs, &
    Gamma_hole_imp, Gamma_tddft, Green_s, Green_self, Green_int, Hedin, Hubbard, korigimp, lmaxfree, lmoins1, lplus1, M1M1, &
    M1M2, M2M2, Magnetic, matper, muffintin, noncentre, nonexc, nonexc_imp, normaltau, no_core_resolved, no_dipquad, no_e1e3, &
    no_e2e2, no_e3e3, No_solsing, Octupole, Old_reference, One_run, Optic, Overad, Pas_SCF_imp, Pdb, Perdew, &
    PointGroup_Auto, Polarise, quadmag, Quadrupole, r_self_imp, Readfast, Recup_optic, Recup_tddft_data, Relativiste, &
    rpalf, rydberg, Save_optic, Save_tddft_data, SCF, SCF_elecabs, SCF_mag_free, Self_abs, self_cons, self_exc_imp, self_nonexc, &
    self_nonexc_imp, solsing_only, solsing_s, spherical_signal, spherical_tensor, spherique, &
    Spinorbite, State_all, State_all_out, Supermuf, Symauto, Symmol, Taux, Tddft, Temperature, &
    Trace_format_wien, Use_FDMX, Ylm_comp_inp

  logical, dimension(5):: SCF_log
  logical, dimension(10):: Multipole
  logical, dimension(ngroup):: Atom_nsph_e
  logical, dimension(0:ntype):: Hubb

  real(kind=db):: Alfpot, Ang_borm, D_max_pot, Delta_En_conv, Delta_Epsii, Eclie, Eclie_out, g1, g2, Gamma_max, &
    number_from_text, Overlap, p_self_max, p_self0, Pas_SCF, phi, pop_nsph, pp, q, r, R_rydb, rad, Rmtt, rn, Roverad, Rpotmax, &
    Step_azim, t, tc, Temp, Test_dist_min, Theta, V_intmax, vv

  real(kind=db):: Kern_fac, q_nrixs_first, q_nrixs_step, q_nrixs_last, r_self
  real(kind=db), dimension(nq_nrixs):: q_nrixs
  real(kind=db), dimension(2):: f_no_res
  real(kind=db), dimension(3):: Ang, Ang_rotsup, Ang_spin, angxyz, Axe, Axe_spin, axyz, Centre, dpos, p, Vec_orig
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(3,3):: Cubmat, Cubmati, Mat, Mat_or, Rot, Rot_gen
  real(kind=db), dimension(norbdil):: cdil
  real(kind=db), dimension(0:ntype):: r0_lapw, rchimp, rlapw, rmt, rmtimp, V_hubbard
  real(kind=db), dimension(neimagent):: eeient, eimagent
  real(kind=db), dimension(ngamme):: egamme
  real(kind=db), dimension(nspin):: ecrantage, V0bdcFimp
  real(kind=db), dimension(n_adimp):: Adimp
  real(kind=db), dimension(n_radius):: Rsorte_s
  real(kind=db), dimension(n_adimp):: E_adimp
  real(kind=db), dimension(n_radius):: E_radius
  real(kind=db), dimension(n_range):: E_max_range
  real(kind=db), dimension(3,nple):: polar, veconde
  real(kind=db), dimension(nple,2):: pdpolar
  real(kind=db), dimension(3,0:ntype) :: Ang_base_loc
  real(kind=db), dimension(n_file_dafs_exp):: Angle_dafs_exp
  real(kind=db), dimension(ngroup_taux):: Taux_oc
  real(kind=db), dimension(ngroup_temp):: Temp_coef
  real(kind=db), dimension(3,ngroup):: Ang_base_loc_gr, posn
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_gr
  real(kind=db), dimension(-m_hubb_e:m_hubb_e,-m_hubb_e:m_hubb_e,nspin,ngroup_taux):: occ_hubb_e
  real(kind=db), dimension(3,3,ngroup_m):: Rot_atom_gr
  real(kind=db), dimension(3,3,ngroup_lapw):: rotloc_lapw
  real(kind=db), dimension(npldafs):: Angle_or
  real(kind=db), dimension(3,npldafs):: angpoldafs, poldafse, poldafsei, poldafss, poldafssi, vecdafsem, vecdafssm
  real(kind=db), dimension(ngroup_par,nparm):: param
  real(kind=db), dimension(nhybm,ngroup_nonsph):: pop_nonsph
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats
  real(kind=db), dimension(0:ntype,nlatm,nspin):: popval
  real(kind=db), dimension(3,nslapwm):: Wien_taulap

  real(kind=db), dimension(:), allocatable:: Hyb, x
  real(kind=db), dimension(:,:), allocatable:: pos
  real(kind=db), dimension(:,:,:), allocatable:: Hybrid_i, Hybrid_r

! Parametres par defaut
  adimp(:) = 0.25_db
  adimpin = .false. !*** JDB
  alfpot = 0._db
  All_nrixs = .false.
  Allsite = .false.
  Ang_base_loc(1,:) = -10000._db; Ang_base_loc(2:3,:) = 0._db
  Ang_base_loc_gr(1,:) = -10000._db; Ang_base_loc_gr(2:3,:) = 0._db
  Ang_rotsup(:) = 0._db
  Angpoldafs(:,:) = 0._db
  ATA = .false.
  Atom = .false.
  Atom_conf = .false.
  Atomic_scr = .false.
  Ang_spin(:) = 0._db
  Axe_spin(1) = 0._db; Axe_spin(2) = 0._db; Axe_spin(3) = 1._db
  Base_spin = .false.
  Basereel = .true.
  cdil(:) = 0._db
  Cartesian_tensor = .false.
  Centre(:) = 0._db
  Centre_auto = .false.
  Centre_auto_abs = .false.
  Clementi = .false.
  com(:) = ' Dirac'
  E1E1 = .true.
  E1E2e = .false.
  E1E3 = .false.
  E1M1 = .false.
  E1M2 = .true.
  E2E2 = .false.
  E3E3 = .false.
  M1M1 = .false.
  M1M2 = .false.
  M2M2 = .false.
  Multipole(1) = E1E1; Multipole(2) = E1E2e; Multipole(3) = E1E3;
  Multipole(4) = E1M1; Multipole(5) = E1M2; Multipole(6) = E2E2;
  Multipole(7) = E3E3; Multipole(8) = M1M1; Multipole(9) = M1M2;
  Multipole(10) = M2M2
  Coupelapw = .false.
  Cylindre = .false.
  D_max_pot = 2.5_db
  Dafs = .false.
  Dafs_bio = .false.
  dipmag = .false.
  dyn_eg = .false.
  dyn_g = .false.
  dpos(:) = 0._db
  E_adimp(:) = 100000000._db
  E_max_range(:) = 100000000._db
  E_radius(:) = 100000000._db
  Ecrantage(:) = 1._db / nspin
  Eclie = 0.2_db
  Eclie_out = 1._db
  Delta_En_conv = 0.1_db
  Delta_Epsii = 1000000._db * Rydb
  Density = .false.
  Density_comp = .false.
  Eneg_i = .false.
  Eneg_n_i = .false.
  Energphot = .false.
  Exc_imp = .false.
  f_no_res(1) = 1._db     ! mag
  f_no_res(2) = -100._db  ! mom
  Fdm_comp = .false.
  Fermi_auto = .true.
  Flapw_new = .false.
  Force_ecr = .false.
  Full_atom_e = .false.
  Full_potential = .false.
  Gamma_tddft = .false.
  Green_int = .false.
  Green_s = .false.
  Green_self = .true.
  Hedin = .false.
  iabsm(1) = 1
  ier = 0
  Charge_free = .false.
  icom(:) = 1
  igr_dop = 0
  iord = 4
  isigpi(:,:) = 0
  istop = 0
  Hubb(:) = .false.
  Kern_fac = 2._db
  korigimp = .false.
  lamstdens = -1
  ldipimp(:) = -1
  lecrantage = 0
  lquaimp(:,:) = -1
  lin_gam = -1
  lmax_nrixs = 3
  lmaxso0 = -5
  lmaxat0 = -1
  lmaxfree = .false.
  lmax_pot = 0
  lmax_pot_default = 1
  lmoins1 = .false.
  lplus1 = .false.
  matper = .true.
  muffintin = .false.
  multrmax = 1
  n_devide = 2
  nchemin = - 1
  necrantage = 0
  nlat(:) = 0
  no_dipquad = .false.
  no_e1e3 = .false.
  no_e3e3 = .false.
  no_e2e2 = .false.
  No_solsing = .false.
  nonexc = .false.
  nonexc_imp = .false.
  noncentre = .false.
  non_relat = 0
  norbv(:) = 0
  normaltau = .false.
  normrmt = 1
  nphim = 180
  nrato(:) = 0
  nrato_dirac = 600
  nrm = 0
  nself = 0
  nsymextract(:) = 1
  do i = 1,n_multi_run_e
    nposextract(i) = i
  end do
  numat_abs = 0
  if( Atom_occ_hubb ) occ_hubb_e(:,:,:,:) = 0._db
  Octupole = .false.
  old_reference = .true.
  One_run = .false.
  Optic = .false.
  overad = .false.
  overlap = 0.1_db
  p_self_max = 1._db
  p_self0 = 0.1_db
  Pas_SCF_imp = .false.
  Pdb = .false.
  pdpolar(:,:) = 0._db
  perdew = .false.
  polar(:,:) = 0._db
  polarise = .false.
  popats(:,:,:) = 0._db
  quadmag = .false.
  Quadrupole = .false.
  Recup_optic = .false.
  Recup_tddft_data = .false.
  relativiste = .false.
  rchimp(:) = 0._db
  r_self_imp = .false.
  rmt(:) = 0._db
  rmtimp(:) = 0._db
  roverad = 0._db
  rpotmax = 0._db
  rsorte_s(:) = 3._db
  Rydberg = .false.
  R_rydb = 1._db
  Save_optic = .false.
  Save_tddft_data = .false.
  SCF = .false.
  SCF_elecabs = .false.
  SCF_mag_free = .false.
  self_cons = .false.
  self_exc_imp = .false.
  self_nonexc = .true.
  self_nonexc_imp = .false.
  seuil = 'K1'
  solsing_s = .false.
  solsing_only = .false.
  Space_Group = ' '
  spherical_tensor = .false.
  spherical_signal = .false.
  spherique = .false.
  Spinorbite = .false.
  core_resolved_e = .false.
  no_core_resolved = .false.
  State_all = .false.
  State_all_out = .false.
  supermuf = .false.
  PointGroup = ' '
  PointGroup_Auto = .true.
  RPALF = .false.
  Symauto = .true.
  Symmol = .false.
  if( Taux ) Taux_oc(:) = 1._db
  if( Temperature ) Temp_coef(:) = 0._db
  Tddft = .false.
  Temp = 0._db
  Test_dist_min = 0.7_db * bohr ! distance minimum entre 2 atomes
  trace_format_wien = .false.
  Trace_k = 0
  Trace_p(:) = 0._db
  veconde(:,:) = 0._db
  v0bdcFimp(:) = 0._db
  V_hubbard(:) = 0._db
  v_intmax = 1000000 * rydb
  Vec_orig(1:2) = 0._db; Vec_orig(3) = 1._db
  Ylm_comp_inp = .false.
  Z_nospinorbite = 0

!Wien_file(1): struct, (2): vcoul, (3): r2v, (4): r2vdn, (5,6,7): clm(isp), (8) psii, (9) Save_potlapw
  Wien_file(:) = ' '
  Wien_file(8) = 'dirac'
  Wien_save = 0

  if( mpirank0 == 0 ) then

    Rewind(itape4)

    write(6,*)

    boucle_lect: do igrdat = 1,100000

      read(itape4,'(A)',iostat=eof) mot
      if( eof /= 0 ) exit boucle_lect

      grdat = identmot(mot,9)
      if( grdat(1:1) /= ' ' ) write(6,'(3x,A)') grdat

      select case(grdat)

        case('ata')
          ATA = .true.

        case('doping')
          read(itape4,*,iostat=ier) itype_dop, igr_dop
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('fdm_comp')
          FDM_comp = .true.

        case('full_pote')

          Full_potential = .true.
          n = nnombre(itape4,132)
          if( n > 0 ) then
            read(itape4,*,iostat=ier) lmax_pot
            if( ier > 0 ) call write_err_form(itape4,grdat)
          else
            lmax_pot = lmax_pot_default
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

        case('extract')
          n = nnombre(itape4,132)
          read(itape4,'(A)') nom_fich_Extract
          if( nom_fich_Extract(1:1) == ' ' ) nom_fich_Extract = adjustl( nom_fich_Extract )
          l = len_trim(nom_fich_Extract)
          if( nom_fich_Extract(l-3:l) == '_bav' ) nom_fich_Extract(l+1:l+4) = '.txt'

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
          coupelapw = .true.
          read(itape4,*,iostat=ier) Trace_k, Trace_p(:)
          if( ier > 0 ) call write_err_form(itape4,grdat)
          if( grdat == 'trace_for' .or. grdat == 'trace_wie' ) trace_format_wien = .true.

        case('range','rangel')
          n = nnombre(itape4,132)
          if( grdat == 'rangel' ) then
            lin_gam = 1
          else
            lin_gam = 0
          endif
          read(itape4,*,iostat=ier) egamme(1:n)

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
          state_all_out = .true.

        case('supermuf')
          supermuf = .true.

        case('old_refer')
          old_reference = .true.

        case('new_refer')
          old_reference = .false.

        case('etatlie')
          n = nnombre(itape4,132)
          if( n == 1 ) then
            read(itape4,*,iostat=ier) Eclie_out
            Eclie = Eclie_out
          else
            read(itape4,*,iostat=ier) Eclie, Eclie_out
          endif
          if( ier > 0 ) call write_err_form(itape4,grdat)

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
          if( ier > 0 ) call write_err_form(itape4,grdat)
          normaltau = .true.

        case('lmoins1')
          lmoins1 = .true.

        case('lplus1')
          lplus1 = .true.

        case('base_reel')
          basereel = .true.

        case('base_comp')
          basereel = .false.

        case('base_spin')
          base_spin = .true.

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
          if( ier > 0 ) call write_err_form(itape4,grdat)
          Ang_spin(1:n) = Ang_spin(1:n) * pi / 180

        case('axe_spin')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Axe_spin(1:3)
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('rot_sup')
          n = nnombre(itape4,132)
          n = min(3,n)
          read(itape4,*,iostat=ier) Ang_rotsup(1:n)
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('relativis')
          relativiste = .true.

        case('non_relat')
          non_relat = 1

        case('no_res_ma')
          read(itape4,*,iostat=ier) f_no_res(1)
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('no_res_mo')
          read(itape4,*,iostat=ier) f_no_res(2)
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('z_nospino')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Z_nospinorbite
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('debye')
          read(itape4,*,iostat=ier) temp
          if( ier > 0 ) call write_err_form(itape4,grdat)

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

        case('radius')
          n = nnombre(itape4,132)
          if( n == 1 ) then
            read(itape4,*,iostat=ier) rsorte_s(1)
          else
            read(itape4,*,iostat=ier) rsorte_s(1), ( E_radius(i-1), rsorte_s(i), i = 2,n_radius )
          endif
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('nrato')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) nrato_dirac
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('multrmax')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) multrmax
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('rpotmax')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) rpotmax
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('over_rad')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) roverad
          if( ier > 0 ) call write_err_form(itape4,grdat)
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
          if( n == 1 .and. Magnetic ) then
            ecrantage(nspin) = ecrantage(1) / nspin
            ecrantage(1) = ecrantage(nspin)
          endif
          Force_ecr = .true.
          Charge_free = .true.

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
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Kern_fac
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('gamma_tdd')
          Gamma_tddft = .true.

        case('tddft_dat')
          Recup_tddft_data = .true.
          read(itape4,'(A)',iostat=ier) nomfich_tddft_data
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('optic_dat')
          Recup_optic = .true.
          read(itape4,'(A)',iostat=ier) nomfich_optic_data
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('save_tddf')
          Save_tddft_data = .true.

        case('save_opti')
          Save_optic = .true.

        case('atomic_sc')
          Atomic_scr = .true.

        case('edge')
          n = nnombre(itape4,132)
          read(itape4,'(A)') motsb
          seuil = identmot(motsb,3)
          select case( seuil(1:1) )
            case('k')
              seuil = 'K1'
            case('l')
              seuil(1:1) = 'L'
            case('m')
              seuil(1:1) = 'M'
            case('n')
              seuil(1:1) = 'N'
            case('o')
              seuil(1:1) = 'O'
            case('p')
              seuil(1:1) = 'P'
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
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('lquaimp')
          do i = 1,3
            n = nnombre(itape4,132)
            read(itape4,*,iostat=ier) lquaimp(i,1:3)
            if( ier > 0 ) call write_err_form(itape4,grdat)
          end do

        case('normaltau')
          normaltau = .true.

        case('noncentre')
          noncentre = .true.
          state_all = .true.

        case('center')
          noncentre = .true.
          state_all = .true.
          n = nnombre(itape4,132)
          select case(n)
            case(0)
              Centre_auto = .true.
            case(1,2)
              Centre_auto = .true.
              read(itape4,*)
            case default
              read(itape4,*,iostat=ier) Centre(1:3)
              if( ier > 0 ) call write_err_form(itape4,grdat)
          end select

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
          do ipl = 1,nple
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

                call write_error
                do ipr = 6,9,3
                  write(ipr,100)
                  write(ipr,'(//A/)') ' Wrong number of elements under the card Polarized !'
                end do
                stop
            end select
            if( ier > 0 ) call write_err_form(itape4,grdat)
          end do

        case('allsite')
          allsite = .true.

        case('all_nrixs')
          All_nrixs = .true.

        case('lmax_nrix')
          read(itape4,*,iostat=ier) lmax_nrixs

        case('nrixs')
          n = nnombre(itape4,132)
          if( n == 1 ) then
            read(itape4,*,iostat=ier) q_nrixs_first
            q_nrixs_step = 1._db
            q_nrixs_last = q_nrixs_first
          else
            read(itape4,*,iostat=ier) q_nrixs_first, q_nrixs_step, q_nrixs_last
          endif
          if( ier > 0 ) call write_err_form(itape4,grdat)
          q_nrixs(1) = q_nrixs_first
          do i = 2,nq_nrixs
            q_nrixs(i) = q_nrixs(i-1) + q_nrixs_step
          end do

        case('step_azim')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Step_azim
          if( ier > 0 ) call write_err_form(itape4,grdat)
          if( abs( Step_azim ) < eps10 ) then
            nphim = 1
          else
            nphim = nint( 360._db / abs( Step_azim ) )
          endif

        case('symsite')
          symauto = .false.
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) n_atom_proto
          if( ier > 0 ) call write_err_form(itape4,grdat)
          igr = 0
          do ipr = 1,n_atom_proto
            n = nnombre(itape4,132)
            read(itape4,*,iostat=ier) natomsym
            if( ier > 0 ) call write_err_form(itape4,grdat)
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
              if( ier > 0 ) call write_err_form(itape4,grdat)
              write(iscratch,*) igr, p(:), isymeq
            end do
          end do

        case('dafs')
          Dafs = .true.
          do ipl = 1,npldafs
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            select case(n)
              case(3)
                read(itape4,*) hkl_dafs(:,ipl)
                nn = nnombre(itape4,132)
                select case(nn)
                  case(6)
                    read(itape4,*,iostat=ier) poldafse(1:3,ipl), vecdafsem(1:3,ipl)
                    read(itape4,*,iostat=ier) poldafss(1:3,ipl), vecdafssm(1:3,ipl)
                    angpoldafs(3,ipl) = 10000._db
                  case(9)
                    read(itape4,*,iostat=ier) (poldafse(i,ipl), poldafsei(i,ipl), i = 1,3), vecdafsem(1:3,ipl)
                    read(itape4,*,iostat=ier) (poldafss(1:3,ipl), poldafssi(i,ipl), i = 1,3), vecdafssm(1:3,ipl)
                    angpoldafs(3,ipl) = 10000._db
                  case default
                    call write_error
                    do ipr = 6,9,3
                      write(ipr,110) ipl
                    end do
                    stop
                end select
                if( ier > 0 ) call write_err_form(itape4,grdat)
              case(5)
                read(itape4,*,iostat=ier) hkl_dafs(:,ipl), isigpi(ipl,1:2)
                angpoldafs(3,ipl) = - 10000._db
                do i = 1,2
                  if( isigpi(ipl,i) == 1 ) then
                    angpoldafs(i,ipl) = 0._db
                  elseif( isigpi(ipl,i) == 2 ) then
                    angpoldafs(i,ipl) = 90._db
                  endif
                end do
              case(6)
                read(itape4,*,iostat=ier) hkl_dafs(:,ipl), isigpi(ipl,1:2), angpoldafs(3,ipl)
                do i = 1,2
                  if( isigpi(ipl,i) == 1 ) then
                    angpoldafs(i,ipl) = 0._db
                  elseif( isigpi(ipl,i) == 2 ) then
                    angpoldafs(i,ipl) = 90._db
                  endif
                end do
              case(7)
                read(itape4,*,iostat=ier) hkl_dafs(:,ipl), (isigpi(ipl,i), angpoldafs(i,ipl), i = 1,2)
                angpoldafs(3,ipl) = - 10000._db
              case(8)
                read(itape4,*,iostat=ier) hkl_dafs(:,ipl), (isigpi(ipl,i), angpoldafs(i,ipl), i = 1,2), angpoldafs(3,ipl)
              case default
                call write_error
                do ipr = 6,9,3
                   write(ipr,100)
                   write(ipr,120) ipl
                end do
                stop
            end select
            if( ier > 0 ) call write_err_form(itape4,grdat)
          end do
          do ipl = 1,npldafs
            do i = 1,2
              if( isigpi(ipl,i) == 3 .or. isigpi(ipl,i) == 4 .or. isigpi(ipl,i) == 10 ) cycle
              if( abs( angpoldafs(i,ipl) ) < eps10 ) then
                isigpi(ipl,i) = 1
              elseif( abs( angpoldafs(i,ipl) - 90 ) < eps10 ) then
                isigpi(ipl,i) = 2
              else
                isigpi(ipl,i) = 5
              endif
            end do
            do i = 1,2
              if( isigpi(ipl,i) == 10 ) angpoldafs(i,ipl) = -10000._db
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

        case('dafs_exp')
          Dafs_bio = .true.
          Dafs = .true.
          do i = 1,3
            read(itape4,*,iostat=ier) Mat_or(i,:)
            if( ier > 0 ) call write_err_form(itape4,grdat)
          end do
          ipl0 = 0
          do i = 1,n_file_dafs_exp
            read(itape4,*,iostat=ier) Angle_dafs_exp(i)
            if( ier > 0 ) call write_err_form(itape4,grdat)
            n = nnombre(itape4,132)
            read(itape4,'(A)') Fichier
            Fichier = Adjustl(Fichier)
            l = len_trim(Fichier)
            if( l > 4 ) then
              if( Fichier(l-3:l-3) /= '.' ) Fichier(l+1:l+4) = '.txt'
            endif
            open(99, file=Fichier, status='old',iostat=istat)
            File_dafs_exp(i) = Fichier
            n = nnombre(99,100000)
            n = n / 3
            read(99,*,iostat=ier) (hkl_dafs(:,ipl),ipl=ipl0+1,ipl0+n)
            if( ier > 0 ) call write_err_form(99,grdat)
            Close(99)
            Angle_or(ipl0+1:ipl0+n) = Angle_dafs_exp(i)
            ipl0 = ipl0 + n
          end do

        case('green')
          Green_s = .true.

        case('green_int')
          Green_s = .true.
          Green_int = .true.

        case('zero_azim')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) Vec_orig(:)
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('solsing')
          solsing_only = .true.
          solsing_s = .true.

        case('no_solsin')
          No_solsing = .true.

        case('norman')
          normrmt = 2

        case('rchimp')
          do it = 1,100000
            n = nnombre(itape4,132)
            if( n == 0 ) exit
            if( it <= ntype ) then
              read(itape4,*,iostat=ier) rchimp(it)
              if( ier > 0 ) call write_err_form(itape4,grdat)
            else
              read(itape4,*)
            endif
          end do

        case('raydem')
          normrmt = 3

        case('rmtg')
          normrmt = 4
          n1 = 1
          do i = 1,ntype
            n = nnombre(itape4,132)
            n2 = min( n1+n-1, ntype )
            if( n2 >= n1 ) then
              read(itape4,*,iostat=ier) rmtimp(n1:n2)
              if( ier > 0 ) call write_err_form(itape4,grdat)
              n1 = n2 + 1
            else
              exit
            endif
          end do
          if( n2 < ntype ) normrmt = - n2
          overlap = 0._db

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
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('muffintin')
          muffintin = .true.

        case('iord')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) iord
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('adimp')
          n = nnombre(itape4,132)
          adimpin = .true.
          if( n == 1 ) then
            read(itape4,*,iostat=ier) adimp(1)
          else
            read(itape4,*,iostat=ier) adimp(1), ( E_adimp(i-1), adimp(i), i = 2,n_adimp )
          endif
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('rmt')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) rmtt
          if( ier > 0 ) call write_err_form(itape4,grdat)
          rmt(:) = rmtt

        case('lmax')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) lmaxat0
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('lmaxfree')
          lmaxfree = .true.

        case('lmaxstden')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) lamstdens
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('lmaxso')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) lmaxso0
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('chlib')
          Charge_free = .true.

        case('hedin')
          hedin = .true.

        case('perdew')
          perdew = .true.

        case('xalpha')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) alfpot
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('flapw','flapw_s','flapw_r','flapw_s_p','flapw_psi', 'flapw_n','flapw_n_p','flapw_s_n')

          if( grdat(6:7) == '_n' ) Flapw_new = .true.
          if( grdat(6:7) == '_s' ) then
            Wien_save = 1
            n = nnombre(itape4,132)
            read(itape4,'(A)') Wien_file(9)
          elseif( grdat(6:7) == '_r' ) then
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
            if( i == 8 .and. ( grdat == 'flapw_s_p' .or. grdat == 'flapw_psi' .or. grdat == 'flapw_n_p' ) ) cycle
            n = nnombre(itape4,132)
            read(itape4,'(A)') Wien_file(i)
            if( Wien_file(i)(1:1) == ' ' ) Wien_file(i) = Adjustl( Wien_file(i) )
          end do

        case('delta_en_')
          read(itape4,*,iostat=ier) Delta_En_conv
          if( ier > 0 ) call write_err_form(itape4,grdat)

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
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('scf_abs')
          scf_elecabs = .true.

        case('p_self')
          read(itape4,*,iostat=ier) p_self0
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('p_self_ma')
          read(itape4,*,iostat=ier) p_self_max
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('n_self')
          read(itape4,*,iostat=ier) nself
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('r_self')
          read(itape4,*,iostat=ier) r_self
          if( ier > 0 ) call write_err_form(itape4,grdat)
          r_self_imp = .true.

        case('scf_exc')
          self_exc_imp = .true.

        case('scf_non_e')
          self_nonexc_imp = .true.

        case('nonexc')
          nonexc = .true.
          nonexc_imp = .true.

        case('excited')
          exc_imp = .true.

        case('no_fermi')
          fermi_auto = .false.

        case('hubbard')
          n1 = 1
          do i = 1,100
            n = nnombre(itape4,132)
            n2 = min( n1+n-1, ntype)
            if( n2 >= n1 ) then
              read(itape4,*,iostat=ier) V_hubbard(n1:n2)
              if( ier > 0 ) call write_err_form(itape4,grdat)
              do j = n1,n2
               if( abs(V_hubbard(j)) > eps10 ) Hubb(j) = .true.
              end do
              n1 = n2 + 1
            else
              exit
            endif
          end do

        case('full_atom')
          Full_atom_e = .true.

        case('absorbeur')
          if( Doping ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,100)
              write(ipr,'(//A//)') ' Absorber keyword not authorized when Doping keyword present !'
            end do
            stop
          endif
          k = 0
          do i = 1,100000
            n = nnombre(itape4,132)
            if( n < 1 ) exit
            read(itape4,*,iostat=ier) iabsm(k+1:k+n)
            if( ier > 0 ) call write_err_form(itape4,grdat)
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
              if( ier > 0 ) call write_err_form(itape4,grdat)
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
                if( ier > 0 ) call write_err_form(itape4,grdat)
                do l = 1,nlat(it)
                  popval(it,l,1) = popval(it,l,1) + x(l)
                end do
                deallocate( x )
              else
                read(itape4,*,iostat=ier) nb_atom_conf(it), igra(1:nb_atom_conf(it),it), nlat(it), &
                                         ( nvval(it,l), lvval(it,l), popval(it,l,:), l = 1,nlat(it) )
                if( ier > 0 ) call write_err_form(itape4,grdat)
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
              if( ier > 0 ) call write_err_form(itape4,grdat)
              n = nnombre(itape4,132)
            endif

            nrm = max( nrm, nrato_dirac )
            if( jt <= ntype ) it = it + 1
            if( n == 1 ) then
              read(itape4,*,iostat=ier) numat(it)
              if( ier > 0 ) call write_err_form(itape4,grdat)
              nlat(it) = 0
            else
              read(itape4,*,iostat=ier) numat(it), nlat(it)
              if( ier > 0 ) call write_err_form(itape4,grdat)
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
                    read(itape4,*,iostat=ier) numat(it),nlat(it), ( nvval(it,l), lvval(it,l), popval(it,l,1), x(l), l = 1,nlat(it) )
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
                if( ier > 0 ) call write_err_form(itape4,grdat)
              endif
            endif

          end do

        case('dilatorb')
          do i = 1,norbdil
            n = nnombre(itape4,132)
            read(itape4,*,iostat=ier) itdil(i), ldil(i), cdil(i)
            if( ier > 0 ) call write_err_form(itape4,grdat)
          end do

        case('v0imp')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) v0bdcFimp(1:min(n,nspin))
          if( ier > 0 ) call write_err_form(itape4,grdat)
          if( nspin > n ) v0bdcFimp(nspin) = v0bdcFimp(1)
          korigimp = .true.

        case('vmax')
          n = nnombre(itape4,132)
          read(itape4,*,iostat=ier) v_intmax
          if( ier > 0 ) call write_err_form(itape4,grdat)

        case('eimag')
          do ie = 1,neimagent
            n = nnombre(itape4,132)
            if( n == 1 ) then
              read(itape4,*,iostat=ier) eimagent(ie)
              eeient(ie) = 0._db
            else
              read(itape4,*,iostat=ier) eeient(ie), eimagent(ie)
            endif
            if( ier > 0 ) call write_err_form(itape4,grdat)
          end do

        case('pointgrou')
          n = nnombre(itape4,132)
          read(itape4,'(A)') motsb
          PointGroup = identmot(motsb,8)
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
          read(itape4,'(A)') motsb
          if( motsb(1:1) == ' ' ) motsb = adjustl(motsb)
          Space_group = motsb(1:13)

        case('crystal','molecule','crystal_t','molecule_')
          if( grdat(1:8) == 'molecule') then
            matper = .false.
          else
            matper = .true.
          endif
          n = nnombre(itape4,132)
          if( n == 3 ) then
            read(itape4,*) axyz(1:3)
            angxyz(1:3) = 90._db
          elseif( n == 2 ) then
            cylindre = .true.
            read(itape4,*) axyz(1:3:2)
            axyz(2) = axyz(1)
            angxyz(1:3) = 90._db
          elseif( n == 1 ) then
            spherique = .true.
            read(itape4,*) axyz(1)
            axyz(2) = axyz(1)
            axyz(3) = axyz(1)
            angxyz(1:3) = 90._db
          else
            read(itape4,*,iostat=ier) axyz(1:3), angxyz(1:3)
          endif
          if( ier > 0 ) call write_err_form(itape4,grdat)
          do igr = 1,ngroup_neq
            if( Temperature .and. Taux ) then
              n = nnombre(itape4,132)
              if( n > 5 ) then
                read(itape4,*) itype(igr), posn(:,igr), Taux_oc(igr), Temp_coef(igr)
              else
                read(itape4,*) itype(igr), posn(:,igr)
              endif
            elseif( Temperature ) then
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
                  hybrid(io,:,igr) = 0._db
                  n = nnombre(itape4,132)
                  allocate( Hyb(n-1) )
                  read(itape4,*) Hyb(1:n-1), pop_nonsph(io,igr)
                  select case(n)
                    Case(5,6,8,17)
                      Hybrid(io,1:n-1,igr) = cmplx( Hyb(1:n-1), 0._db, db )
                    Case(9,11,15,33)
                      do i = 1,n-2,2
                        Hybrid(io,(i+1)/2,igr) = cmplx( Hyb(i), Hyb(i+1), db )
                      end do
                  end select
                  deallocate ( Hyb )
                end do
              endif
              if( norb < 0 ) then
                n = nnombre(itape4,132)
                l = (n/nspin - 1)/2
                read(itape4,*) ((occ_hubb_e(m,m,isp,igr), m = -l,l ), isp = 1,nspin)
              endif
            endif
          end do

        case('pdb_file')

          Pdb = .true.
          Fichier = ' '
          read(itape4,'(A)') Fichier
          Fichier = Adjustl(Fichier)
          l = len_trim(Fichier)
          if( l > 4 ) then
            if( Fichier(l-3:l-3) /= '.' ) Fichier(l+1:l+4) = '.pdb'
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

                l = len_trim(Spgr)
                j = 0
                do i = 1,l
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
                if( Temperature ) Temp_coef(igr) = tc

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

! Description de l'agregat :
        case('cif_file')

          Matper = .true.

          Fichier = ' '
          read(itape4,'(A)') Fichier
          Fichier = Adjustl(Fichier)
          l = len_trim(Fichier)
          if( l > 4 ) then
            if( Fichier(l-3:l-3) /= '.' ) Fichier(l+1:l+4) = '.cif'
          endif
          open(8, file = Fichier, status='old', iostat=istat)
          if( istat /= 0 ) call write_open_error(Fichier,istat,1)

          do

            read(8,'(A)',iostat=eof) mot

            if( eof /= 0 ) exit

            if( mot(1:30) == '_symmetry_space_group_name_H-M' ) then

            l = len_trim(mot)

            do i = 31,l
              if( mot(i:i) == '''' ) exit
            end do
            k = 0
            j = i + 1
            do i = j,l
              if( mot(i:i) == ' ' ) cycle
              if( mot(i:i) == '''' ) exit
              k = k + 1
              Space_group(k:k) = mot(i:i)
            end do

            elseif( mot(1:17) == '_cell_angle_alpha' ) then
              angxyz(1) = number_from_text(1,mot)

            elseif( mot(1:16) == '_cell_angle_beta' ) then
              angxyz(2) = number_from_text(1,mot)

            elseif( mot(1:17) == '_cell_angle_gamma' ) then
              angxyz(3) = number_from_text(1,mot)

            elseif( mot(1:14) == '_cell_length_a') then
              axyz(1) = number_from_text(1,mot)

            elseif( mot(1:14) == '_cell_length_b') then
              axyz(2) = number_from_text(1,mot)

            elseif( mot(1:14) == '_cell_length_c') then
              axyz(3) = number_from_text(1,mot)

            elseif( mot(1:16) == '_atom_site_label') then

              n_fract_x = 0; n_fract_y = 0; n_fract_z = 0

              do n = 1,100000
                read(8,'(A)') mot
                if( mot(1:1) /= '_' ) then
                  backspace(8)
                  exit
                elseif( mot(2:18) == 'atom_site_fract_x') then
                  n_fract_x = n
                elseif( mot(2:18) == 'atom_site_fract_y') then
                  n_fract_y = n
                elseif( mot(2:18) == 'atom_site_fract_z') then
                  n_fract_z = n
                elseif( mot(2:20) == 'atom_site_occupancy') then
                  Taux = .true.
                  n_occupancy = n
                endif
              end do

              do igr = 1,ngroup_neq

                read(8,'(A)') mot

! Lecture element chimique
                Symbol = mot(1:2)
                if( Symbol(2:2) == '1' .or. Symbol(2:2) == '2' .or. Symbol(2:2) == '3' .or. Symbol(2:2) == '4' .or. &
                    Symbol(2:2) == '5' .or. Symbol(2:2) == '6' .or. Symbol(2:2) == '7' .or. Symbol(2:2) == '8' .or. &
                    Symbol(2:2) == '9' ) Symbol(2:2) = ' '

                do i = 1,103
                  if( Chemical_Symbol_c(i) /= Symbol .and. Chemical_Symbol(i) /= Symbol ) cycle
                  itype(igr) = i
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
                if( Taux ) Taux_oc(igr) = number_from_text(n_occupancy,mot)

              end do

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

        case('z_absorbe')
          n = nnombre(itape4,132)
          if( Absauto ) then
            read(itape4,*) numat_abs
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

! Parametres deja lus dans lectdim
        case('full_self','magnetism','memory_sa','self_abs','readfast','temperatu')

        case default

          if( igrdat == 1 ) then
            comt = mot
          elseif( grdat(1:1) /= ' ' ) then
            call write_error
            do ipr = 6,9,3
              write(ipr,100)
              write(ipr,150) mot
            end do
            stop
          endif

      end select

    end do boucle_lect

! Fin de la lecture.

!*** JDB
    if (.NOT. adimpin .AND. Use_FDMX) then
      adimp(1:5) = (/ 0.24, 0.20, 0.16, 0.12, 0.08 /)
      E_adimp(1:4) = (/ 100.0, 250.0, 400.0, 500.0 /)
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

    l = len_trim(nomfich)
    write(6,'(/a9,A)') ' Filout: ',nomfich(1:l)
    if( l > 4 ) then
      if( nomfich(l-3:l) == '.txt' .or. nomfich(l-3:l) == '.dat') nomfich(l-3:l) = '    '
    endif
    nomfichbav = nomfich
    long = len_trim(nomfich)
    nomfichbav(long+1:long+8) = '_bav.txt'

    do i = 1,28
      if( icheck(i) > 1 ) exit
    end do
    if( Extract .and. i > 28 ) icheck(1:28) = 0

    i = sum( icheck(1:28) )
    if( i > 0 ) then
      open(3, file = nomfichbav, status='unknown',iostat=istat)
      if( istat /= 0 ) call write_open_error(nomfichbav,istat,1)
    endif

    if( icheck(1) > 0 ) then
      write(3,'(A/A/A)') Revision, com_date, com_time
      if( comt /= ' ') write(3,'(/A)') comt
      ipr0 = 3
    else
      ipr0 = 6
    endif

    if( Doping ) then
      Ang_base_loc_gr(:,ngroup) = Ang_base_loc_gr(:,igr_dop)
      Atom_nsph_e(ngroup) = Atom_nsph_e(igr_dop)
      itype(ngroup) = itype_dop
      popats(ngroup,:,:) =  popval(itype_dop,:,:)
      posn(:,ngroup) = posn(:,igr_dop)
      if( ngroup_nonsph > 1 ) then
        norbv(ngroup) = norbv(igr_dop)
        Hybrid(:,:,ngroup) = Hybrid(:,:,igr_dop)
        pop_nonsph(:,ngroup) = pop_nonsph(:,igr_dop)
      endif
      if( ngroup_m > 0 ) then
        Axe_atom_gr(:,ngroup) = Axe_atom_gr(:,igr_dop)
        Rot_atom_gr(:,:,ngroup) = Rot_atom_gr(:,:,igr_dop)
      endif
      if( ngroup_taux > 0 ) then
        Taux_oc(ngroup) = Taux_oc(igr_dop)
        occ_hubb_e(:,:,:,ngroup) = occ_hubb_e(:,:,:,igr_dop)
      endif
      if( ngroup_pdb > 0 ) Kgroup(ngroup) = Kgroup(igr_dop)
      if( ngroup_temp > 0 ) Temp_coef(ngroup) = Temp_coef(igr_dop)
    endif

! Atom defined with few digits
    if( abs( Angxyz(3) - 120._db ) < 0.0001_db ) then
      do i = 1,11
        if( i == 3 .or. i == 6 .or. i == 9 ) cycle
        Where( abs( Posn - i / 12._db ) < 0.00005_db ) Posn = i / 12._db
      end do
    endif

! Modification en cas de fit.
    if( fit_cal ) then
      do igr = 2,ngroup_par
        istop = 0
        do ipar = 1,npar(igr)
          if( typepar(igr,ipar) /= 'dposx' .and. typepar(igr,ipar) /= 'dposy' .and. typepar(igr,ipar) /= 'dposz' .and. &
              typepar(igr,ipar) /= 'posx' .and. typepar(igr,ipar) /= 'posy' .and. typepar(igr,ipar) /= 'posz' .and. &
              typepar(igr,ipar) /= 'theta' .and. typepar(igr,ipar) /= 'phi' .and. typepar(igr,ipar) /= 'occup'  ) cycle
          if( indice_par(igr,ipar) > ngroup ) then
            call write_error
            do ipr = ipr0,9,3
              write(ipr,100)
              write(ipr,160) typepar(igr,ipar), indice_par(igr,ipar), ngroup
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

    do igr = 1,ngroup
      if( cylindre ) then
        r = posn(1,igr)
        theta = pi * posn(2,igr) / 180
        posn(1,igr) = r * cos( theta )
        posn(2,igr) = r * sin( theta )
      elseif( spherique ) then
        r = posn(1,igr)
        theta = pi * posn(2,igr) / 180
        phi = pi * posn(3,igr) / 180
        posn(1,igr) = r * sin( theta ) * cos( phi)
        posn(2,igr) = r * sin( theta ) * sin( phi)
        posn(3,igr) = r * cos( theta )
      endif
    end do

    if( Atom_conf ) then

      do it = 1,ntype_conf
        numat(it) = itype(igra(1,it))
      end do

      allocate( itZ(103) )

      jt = ntype_conf

      boucle_2: do igr = 1,ngroup_neq+1
        if( .not. Doping .and. igr == ngroup_neq+1 ) cycle
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
        numat(jt) = itype(igr)
        itZ( numat(jt) ) = jt
      end do boucle_2

      boucle_igr: do igr = 1,ngroup_neq+1
        if( .not. Doping .and. igr == ngroup_neq+1 ) cycle
        do it = 1,ntype_conf
          do kgr = 1,nb_atom_conf(it)
            if( igr /= igra(kgr,it) ) cycle
            itype(igr) = it
            cycle boucle_igr
          end do
        end do
        itype(igr) = itZ( itype(igr) )
      end do boucle_igr

      deallocate( itZ )

      Atom = .true.

    endif

    if( .not. Matper ) Space_group = ' '

    iabsorig(:) = iabsm(:)
    if( Matper .and. Space_group /= ' ' ) then
      allocate( neq(ngroup_neq) )
      allocate( pos(3,ngroup_neq) )
      do igr = 1,ngroup_neq
        pos(:,igr) = posn(:,igr)
      end do
      call spgroup(.true.,neq,ngroup,ngroup_neq,itype,pos,posn,Space_file,Space_group)
      ia = ngroup + 1
      if( Doping ) then
        ia = ia - 1
        do igr = 1,ngroup-1
          if( sum( abs( posn(:,igr) - pos(:,igr_dop) ) ) > eps10 ) cycle
          igr_dop = igr
          exit
        end do
      endif
      do igr = ngroup_neq,1,-1
        do i = 1,neq(igr)
          ia = ia - 1
          itype(ia) = itype(igr)
          if( Taux ) Taux_oc(ia) = Taux_oc(igr)
          if( Temperature ) Temp_coef(ia) = Temp_coef(igr)
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
      deallocate( pos )
      deallocate( neq )
    endif

    if( Doping ) iabsm(1) = ngroup
    if( .not. Taux ) ATA = .FALSE.

    if( Flapw ) then
      nrm = nrato_dirac
      call lect_struct_lapw(angxyz,axyz,icheck(1),its_lapw,itype,ngroup,ngroup_lapw,Wien_file(1),nrato_lapw,nrm,nslapwm,ntype, &
          numat,posn,r0_lapw,rlapw,rotloc_lapw,Wien_matsym,Wien_taulap)
      if( normrmt == 1 ) normrmt = 2
    elseif( nrm == 0 ) then
      nrm = nrato_dirac
    endif

    if( Pdb ) One_run = .false.

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
            write(ipr,171) it, ntype
          else
            write(ipr,172) igr, it, ntype
          endif
        end do
        istop = 1
      end do
      if( istop == 1 ) then
        do ipr = ipr0,9,3
          write(ipr,174)
        end do
        stop
      endif
    endif

    if( numat_abs /= 0 ) then
      do igr = 1,ngroup
        if( Atom .or. Flapw ) then
          if( numat( abs( itype( igr ) ) ) == numat_abs ) exit
        else
          if( itype( igr ) == numat_abs ) exit
        endif
      end do
      if( igr > ngroup ) then
        call write_error
        do ipr = ipr0,9,3
          write(ipr,175) numat_abs
        end do
        stop
      endif
    else
      if( Doping ) then
        if( Atom ) then
          numat_abs = numat( abs( itype_dop ) )
        else
          numat_abs = itype_dop
        endif
      else
        if( Atom .or. Flapw ) then
          numat_abs = numat( abs( itype( iabsm(1) ) ) )
        else
          numat_abs = itype( iabsm(1) )
        endif
      endif
    endif

    if( Clementi ) then
      Wien_file(8) = 'clementi'
      nrm = nrato_dirac
      Fermi_auto = .false.
    endif

    if( Atom_occ_hubb ) Fermi_auto = .false.
    if( Flapw .or. Extract .or. .not. Fermi_auto ) then
      SCF = .false.
      Fermi_auto = .false.
      nself = 0
    end if
    if( .not. r_self_imp ) r_self = rsorte_s(1)
    if( SCF ) then
      Self_cons = .true.
      if( nself == 0 ) nself = 100
    elseif( ( Fermi_auto .or. ( Hubbard .and. .not. Atom_occ_hubb) ) .and. nself == 0 ) then
      Self_cons = .true.
      p_self0 = 0._db
      nself = 1
      if( .not. ( r_self_imp .or. One_run ) ) r_self = min( rsorte_s(1), 3.5_db )
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

    if( Recup_tddft_data ) Save_tddft_data = .false.
    if( Recup_optic ) Save_optic = .false.

    if( .not. Pas_SCF_imp ) then
      if( Green_self ) then
        Pas_SCF = 0.1_db
      else
        Pas_SCF = 0.01_db
      endif
    endif

    if( Extract ) then
      Density = .false.
      state_all = .false.
      state_all_out = .false.

      open(1, file = nom_fich_Extract, status='old', iostat=istat)
      if( istat /= 0 ) call write_open_error(nom_fich_Extract,istat,1)
      do i = 1,100000
        read(1,'(A)') mot
        if( mot(2:10) == 'Threshold' ) then
          l = len_trim(mot)
          if( mot(l:l) == 'e' ) then
            seuil = mot(l-7:l-5)
            seuil = adjustl( seuil )
          else
            seuil = mot(14:16)
          endif
          exit
        endif
        if( Seuil == 'Opt' ) Optic = .true.
      end do
      if( Quadrupole ) then
        do i = 1,100000
          read(1,'(A)') mot
          if( mot(2:6) == 'Dipol'  ) then
            read(1,'(A)') mot
            if( mot(2:6) /= 'Quadr' ) Quadrupole = .false.
            read(1,'(A)') mot
            if( mot(2:6) /= 'Octup' ) octupole = .false.
            exit
          endif
        end do
      endif
      Rewind(1)
      Core_resolved = .false.
      do i = 1,100
        read(1,'(A)') mot
        if( mot(2:5) /= 'Core'  ) cycle
        Core_resolved = .true.
        exit
      end do
      Close(1)
    endif

    if( Spinorbite .and. non_relat == 0 ) relativiste = .true.

    if( .not. Flapw .and. abs(alfpot) < eps6 .and. .not. Perdew ) Hedin = .true.
    if( Hedin ) then
      alfpot = 0._db
    elseif( Flapw ) then
      alfpot = 0.33_db
    elseif( Perdew ) then
      alfpot = -1._db
    endif

    if( Optic ) Seuil = 'Opt'

    if( seuil /= 'K1' .and. seuil /= 'L1' .and. seuil /= 'L2' .and. &
        seuil /= 'L3' .and. seuil /= 'M1' .and. seuil /= 'M2' .and. &
        seuil /= 'M3' .and. seuil /= 'M4' .and. seuil /= 'M5' .and. &
        seuil /= 'N1' .and. seuil /= 'N2' .and. seuil /= 'N3' .and. &
        seuil /= 'N4' .and. seuil /= 'N5' .and. seuil /= 'N6' .and. &
        seuil /= 'N7' .and. seuil /= 'O1' .and. seuil /= 'O2' .and. &
        seuil /= 'O3' .and. seuil /= 'O4' .and. seuil /= 'O5' .and. &
        seuil /= 'P1' .and. seuil /= 'P2' .and. seuil /= 'P3' .and. &
        seuil /= 'L23' .and. seuil /= 'M23' .and. seuil /= 'M45' .and. &
        seuil /= 'N23' .and. seuil /= 'N45' .and. seuil /= 'N67' .and. &
        seuil /= 'O23' .and. seuil /= 'O45' .and. seuil /= 'P23' .and. seuil /= 'Opt' ) then
      call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,180) seuil
      end do
      stop
    endif

    if( Optic ) then
      nseuil = 0; lseuil = 0; jseuil = 0; nbseuil = 1
    else
      select case( seuil(1:1) )
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
      select case( seuil(2:2) )
        case('1')
          lseuil = 0; jseuil = 1
        case('2')
          lseuil = 1; jseuil = 2
        case('3')
          lseuil = 1; jseuil = 3
        case('4')
          lseuil = 2; jseuil = 4
        case('5')
          lseuil = 2; jseuil = 5
        case('6')
          lseuil = 3; jseuil = 6
        case('7')
          lseuil = 3; jseuil = 7
      end select
      nbseuil = len_trim( seuil ) - 1
    endif

    if( Flapw .or. lseuil > 0 .or. Optic ) nonexc = .true.
    if( exc_imp .and. .not. flapw ) nonexc = .false.
    if( nonexc ) Symmol = .true.
    if( self_nonexc_imp ) self_nonexc = .true.
    if( self_exc_imp .and. .not. Optic ) self_nonexc = .false.
    if( self_exc_imp .and. .not. Optic .and. .not. nonexc_imp ) nonexc = .false.
    if( nonexc ) self_nonexc = .true.
    l = l_level_val(numat_abs)
    if( .not. ( nonexc .and. self_nonexc ) .and. ( l == 2 .or. l == 3 ) ) scf_elecabs = .true.

    if( .not. Flapw ) then

      if( .not. atom ) then
! Dans ce cas itype est pour l'instant le numero atomique
        jt = 0
        boucle_1: do igr = 1,ngroup
          n = itype(igr)
          do it = 1,jt
            if( numat(it) /= n ) cycle
            itype(igr) = it
            if( Doping .and. igr == ngroup ) itype_dop = jt
            cycle boucle_1
          end do
          jt = jt + 1
          if( Doping .and. igr == ngroup ) itype_dop = jt
          itype(igr) = jt
          numat(jt) = n
        end do boucle_1
        nlat(1:jt) = 0
      endif

      do igr = 1,ngroup
        it = abs( itype(igr) )
        do l = 1,nlat(it)
          do ispin = 1,nspin
            popats(igr,l,ispin) = popval(it,l,ispin)
          end do
        end do
      end do

    endif

    do it = 1,ntype
      if( numat(it) == 0 ) then
        com(it) = ' Empty sphere'
        icom(it) = 5
      endif
    end do

    if( Extract ) then
      if( no_core_resolved ) Core_resolved = .false.
    elseif( ( nspin == 1 .or. no_core_resolved ) .and. .not. Core_resolved_e ) then
      Core_resolved = .false.
    else
      Core_resolved = .true.
    endif
! S'il y a plusieurs seuils, le Gamma doit etre pris dans le calcul de Chi_0.
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
        if( lseuil == 1 .and. nseuil == 2 ) Gamma_tddft = .true.
      endif
    endif

! Verification des entrees :

    if( .not. ( E1E1 .or. E1E2e .or. E1E3 .or. E2E2 .or. E1M1 .or. E1M2 .or. M1M1 .or. M1M2 .or. M2M2 ) ) then
      if( istop == 0 ) call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,190)
      end do
      istop = 1
    endif

    if( Bormann ) then
      Dafs = .true.
      Angpoldafs(3,:) = Ang_borm
      do ipl = 1,npldafs/4
       hkl_dafs(:,ipl) = hkl_borm(:)
      end do
      do ipl = npldafs/4 + 1,3*npldafs/4
        hkl_dafs(:,ipl) = 0
      end do
      do ipl = 3*npldafs/4 + 1,npldafs
        hkl_dafs(:,ipl) = - hkl_borm(:)
      end do
    endif

    if( ngroup == 0 ) then
      if( istop == 0 ) call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,200)
      end do
      istop = 1
    endif
    if( ntype == 0 .and. ngroup /= 0 ) then
      if( istop == 0 ) call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,210)
      end do
      istop = 1
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
            if( isigpi(ipl,1) == 1 .and. isigpi(ipl,2) == 1 ) cycle
          case(2)
            if( isigpi(ipl,1) == 1 .and. isigpi(ipl,2) == 2 ) cycle
          case(3)
            if( isigpi(ipl,1) == 2 .and. isigpi(ipl,2) == 1 ) cycle
          case(4)
            if( isigpi(ipl,1) == 2 .and. isigpi(ipl,2) == 2 ) cycle
        end select
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,260) ipl, hkl_dafs(:,ipl), isigpi(ipl,:)
        end do
        istop = 1
        exit
      end do
    endif

    if( n_radius > 1 ) then
      do i = 2,n_radius
        if( Rsorte_s(i-1) > Rsorte_s(i) - 0.0001_db ) cycle
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,'(/A/)') ' Under keyword "Radius", the radius values must be decreasing !'
        end do
        istop = 1
      end do
      do i = 2,n_radius
        if( E_radius(i-1) < E_radius(i) - 0.0001_db ) cycle
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,'(/A/)') ' Under keyword "Radius", the energy values must be increasing !'
        end do
        istop = 1
      end do
    endif

    if( n_adimp > 1 ) then
      do i = 2,n_adimp
        if( E_adimp(i-1) < E_adimp(i) - 0.0001_db ) cycle
        if( istop == 0 ) call write_error
        do ipr = ipr0,9,3
          write(ipr,100)
          write(ipr,'(/A/)') ' Under keyword "Adimp", the energy values must be increasing !'
        end do
        istop = 1
      end do
    endif

    if( One_run .and. n_range > 1 ) then
      if( istop == 0 ) call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,'(/A/)') ' It is not possible to have "One_run" keyword and several values of adimp or radius !'
      end do
      istop = 1
    endif

    if( istop == 1 ) stop

! Normalisation des vecteurs hybridations :
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

    if( lin_gam == - 1 ) then
      lin_gam = 1
      egamme(1) = -5.0_db; egamme(2) =  0.5_db;
      egamme(3) = egamme(1) + ( nenerg - 1 ) * egamme(2)
    endif
    if( .not. green_s .and. ( .not. FDM_comp ) .and. neimagent > 0 ) eimagent(:) = 0._db

! Ecriture des entrees :
    if( icheck(1) > 0 ) then

      if( Extract ) then
        write(3,300) nom_fich_Extract
        write(6,300) nom_fich_Extract
        open(1, file = nom_fich_Extract, status='old',iostat=istat)
        if( istat /= 0 ) call write_open_error(nom_fich_Extract,istat,1)
        do i = 1,100000
          read(1,'(A)') mot
          write(3,'(A)') mot
          if( mot(2:6) == 'Dipol'  ) exit
        end do
      else
        mot13 = Chemical_Name(numat_abs)
        l1 = len_trim(mot13)
        l2 = len_trim(seuil)
        mot = ' '
        mot =' Threshold:'
        mot(13:12+l1) = mot13(1:l1)
        mot(14+l1:13+l1+l2) = seuil(1:l2)
        mot(15+l1+l2:18+l1+l2) = 'edge'
        write(3,'(/A)') mot
        write(6,'(A)') mot(12:18+l1+l2)
        if( n_radius == 1 ) then
          write(3,310) Rsorte_s(1)
        else
          write(3,315) Rsorte_s(1), ( E_radius(i-1), Rsorte_s(i), i = 2,n_radius )
        endif
        if( .not. green_s .and. overad ) write(3,320) roverad
        write(3,330) icheck(:)
        if( lin_gam == 1 ) write(3,340)
        write(3,350) egamme(1:ngamme)
        write(3,'(A)') ' Dipole component'
      endif
      if( Quadrupole ) write(3,'(A)') ' Quadrupole component'
      if( Octupole ) write(3,'(A)') ' Octupole component'
      if( Dipmag ) write(3,'(A)') ' Magnetic dipole component'
      if( .not. E1E1 ) write(3,'(A)') ' ... but E1E2 component neglected in the output'
      if( no_e2e2 ) write(3,'(A)') ' ... but E2E2 component neglected in the output'
      if( no_e1e3 ) write(3,'(A)') ' ... but E1E3 component neglected in the output'
      if( no_e3e3 ) write(3,'(A)') ' ... but E3E3 component neglected in the output'
      if( no_dipquad ) write(3,'(A)') ' ... but E1E2 component neglected in the output'
      if( tddft ) then
        if( rpalf ) then
          write(3,'(A)') ' TDDFT calculation using the RPA LF approximation'
        else
          write(3,'(A)') ' TDDFT calculation'
        endif
        if( Gamma_tddft ) then
          write(3,'(A)') '    Broadening in the Chi_0 calculation'
        else
          write(3,'(A)') '    No broadening in the Chi_0 calculation'
        endif
      endif
      if( Core_resolved ) write(3,'(A)') ' Core resolved in outputs'

      if( Extract ) then
        do i = 1,100000
          read(1,'(A)') mot
          if( mot(2:6) /= 'Relat' .and. mot(6:10) /= 'relat' ) cycle
          write(3,'(A)') mot
          exit
        end do
        do i = 1,100000
          read(1,'(A)') mot
          if(  mot(2:7) == 'ngroup' .or. mot(2:6) == 'XANES' .or. mot(2:5) == 'DAFS' ) then
            backspace(1)
            exit
          endif
          write(3,'(A)') mot
        end do
      else
        if( relativiste ) then
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
        if( basereel ) then
          write(3,'(A)') ' Real bases'
        else
          write(3,'(A)') ' Complex bases'
        endif
        if( green_s ) then
          write(3,'(A)') ' Multiple scattering calculation (Green)'
          if( Green_int ) write(3,'(A)') '    Full Green mode'
          if( supermuf ) write(3,'(A)') '    Continuous potential (Supermuf)'
          if( nchemin /= - 1 ) write(3,410) nchemin
          if( normaltau ) write(3,'(A)') '    Tau normalization'
          select case(normrmt)
            case(1)
              write(3,'(A)') '    Optimized type muffin-tin radius'
            case(2)
              write(3,'(A)') '    Norman type muffin-tin radius'
            case(3)
              write(3,'(A)') '    Half interatomic distance type muffin-tin radius'
            case(4)
              write(3,420) rmtimp(1:ntype)
            case(5)
              write(3,'(A)') '    Potential imposed type muffin-tin radius'
          end select
          if( abs(overlap) > eps6 ) write(3,430) overlap
          if( lmaxfree ) then
            write(3,'(A)') '    No limitation on the maximum value of l'
          else
            write(3,'(A)') '    Limitation on the maximum value of l'
          endif
        else
          write(3,'(A)') ' Finite difference method calculation'
          if( n_adimp == 1 ) then
            write(3,440) iord, adimp(1)
          else
            write(3,445) iord, adimp(1), ( E_adimp(i-1), adimp(i), i = 2,n_adimp )
          endif
          write(3,450) lmaxso0
          if( muffintin ) write(3,'(A)') '    Muffin-tin potential'
          if( rydberg ) write(3,460) R_rydb
          if( noncentre ) write(3,'(A)') '    Non centered absorbing atom'
          if( sum( abs(Centre(:)) ) > epspos ) write(3,470) Centre(:)
          if( .not. Eneg_i ) write(3,480) Eclie, Eclie_out
          if( v_intmax < 100000._db ) write(3,490) V_intmax
        endif
      endif
      if( Atom_nonsph ) write(3,'(/A)') ' Calculation with non spherical orbitals'
      if( temp > eps10 ) write(3,500) Temp
      if( Temperature ) write(3,'(/A)') ' Temperature coefficients taken into account for diffraction'

      if( nple > 0 ) then
        if( sum( abs( pdpolar(:,:) ) ) > eps10 ) then
          write(3,'(/A)') ' XANES :    Polarization             Wave vector     Weight_dip Weight_quad'
          do ipl = 1,nple
            write(3,510) polar(1:3,ipl), veconde(1:3,ipl), pdpolar(ipl,:)
          end do
        else
          write(3,'(/A)') ' XANES :    Polarization             Wave vector'
          do ipl = 1,nple
            write(3,510) polar(1:3,ipl), veconde(1:3,ipl)
          end do
        endif
      endif
      if( nq_nrixs >= 1 ) then
        write(3,'(/A)') ' NRIXS (X-ray Raman) calculation'
        write(3,'(a15,i2)') '   lmax_nrixs =' ,lmax_nrixs
        write(3,'(a28,3f8.5,a6,i4)') '   q_first, q_step, q_last =' ,q_nrixs_first, q_nrixs_step, q_nrixs_last, ', nq =',nq_nrixs
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

        write(3,511) hkl_dafs(:,1), Angpoldafs(3,1)

      else

        do ipl = 1,npldafs
          motpol = ' '
          l = 0
          do i = 1,2
            if( isigpi(ipl,i) == 1 ) then
              motpol(l+1:l+5) = 'sigma'
            elseif( isigpi(ipl,i) == 2 ) then
              motpol(l+1:l+2) = 'pi'
            elseif( isigpi(ipl,i) == 3 ) then
              motpol(l+1:l+5) = 'right'
            elseif( isigpi(ipl,i) == 4 ) then
              motpol(l+1:l+4) = 'left'
            elseif( isigpi(ipl,i) == 5 .or. isigpi(ipl,i) == 10 ) then
              motpol(l+1:l+5) = 'recti'
            endif
            l = l + 6
          end do

          if( angpoldafs(3,ipl) > 9999._db ) then
            if( ipl == 1 ) write(3,'(/A)') ' DAFS :     Polarization             Wave vector '
            write(3,512) poldafse(1:3,ipl), vecdafsem(1:3,ipl), ' incoming'
            write(3,512) poldafss(1:3,ipl), vecdafssm(1:3,ipl), ' outgoing'
          else
            if( ipl == 1 ) write(3,'(/A)') ' DAFS : (h, k, l)  Polarization   Angle_i   Angle_o  Azimuth'
            if( angpoldafs(3,ipl) < -9999._db ) then
              write(3,514) hkl_dafs(:,ipl), motpol, angpoldafs(1:2,ipl)
            elseif( isigpi(ipl,1) == 10 ) then
              write(3,517) hkl_dafs(:,ipl), motpol, angpoldafs(2:3,ipl)
            elseif( isigpi(ipl,2) == 10 ) then
              write(3,518) hkl_dafs(:,ipl), motpol, angpoldafs(1:3:2,ipl)
            else
              write(3,519) hkl_dafs(:,ipl), motpol, angpoldafs(1:3,ipl)
            endif
          endif

        end do

      endif

      if( Self_abs ) write(3,'(/A)') ' XANES for the RXS polarizations'
      if( Full_self_abs ) write(3,'(/A)') ' Self absorption taken into account with birefringence effect'

      if( Extract ) then

        Close(1)

      else

        if( Matper ) then
          write(3,'(/A)') ' Crystal'
        else
          write(3,'(/A)') ' Molecule'
        endif
        write(3,542) ngroup, ntype

        if( nlatm > 0 ) then
          write(3,'(A)') '   Typ  Z   n  l  Popul.'
          do it = 1,ntype
            if( nspin == 1 ) then
              write(3,543) it, numat(it), ( nvval(it,l), lvval(it,l), popval(it,l,:), l = 1,nlat(it) )
            else
              write(3,544) it, numat(it), ( nvval(it,l), lvval(it,l), popval(it,l,:), l = 1,nlat(it) )
            endif
          end do
        endif

        if( abs( dpos(1) ) > epspos .or. abs( dpos(2) ) > epspos .or. abs( dpos(3) ) > epspos ) write(3,545) dpos(1:3)
        if( norbdil > 0 ) write(3,548)
        do io = 1,norbdil
          write(3,549) itdil(io), ldil(io), cdil(io)
        end do
        if( nonexc ) write(3,'(A)') '   Non excited absorbing atom'
        if( .not. PointGroup_Auto ) write(3,552) PointGroup
        if( Symmol ) write(3,'(A)') '   Absorbing atom taken as not excited for the symmetry calculation'
        if( Space_group /= ' ' .and. matper ) write(3,554) Space_group

        write(3,555) axyz(1:3)
        write(3,556) angxyz(1:3)
        if( Pdb .and. Taux .and. Temperature ) then
          write(3,560) 'Kgroup  Occupancy  Temp_cf'
        elseif( Pdb .and. Taux ) then
          write(3,560) 'Kgroup  Occupancy'
        elseif( Pdb .and. Temperature ) then
          write(3,560) 'Kgroup   Temp_cf'
        elseif( Pdb ) then
          write(3,560) 'Kgroup'
        elseif( Taux .and. Temperature ) then
          write(3,560) ' Occupancy  Temp_cf'
        elseif( Taux ) then
          write(3,560) ' Occupancy'
        elseif( Temperature ) then
          write(3,560) ' Temp_cf'
        elseif( Atom_nonsph ) then
          write(3,560) 'norbv   popats'
        else
          write(3,560)
        endif
        do igr = 1,ngroup
          if( Doping .and. igr == ngroup ) write(3,'(/A)') '  Doping element :'
          it = itype(igr)
          Z = numat( abs(it) )
          if( Pdb .and. Taux .and. Temperature ) then
            write(3,565) Z, it, posn(:,igr), Kgroup(igr), Taux_oc(igr), Temp_coef(igr)
          elseif( Pdb .and. Taux ) then
            write(3,565) Z, it, posn(:,igr), Kgroup(igr), Taux_oc(igr)
          elseif( Pdb .and. Temperature ) then
            write(3,570) Z, it, posn(:,igr), Kgroup(igr), Temp_coef(igr)
          elseif( Pdb ) then
            write(3,570) Z, it, posn(:,igr), Kgroup(igr)
          elseif( Taux .and. Temperature ) then
            write(3,571) Z, it, posn(:,igr), Taux_oc(igr), Temp_coef(igr)
          elseif( Taux ) then
            write(3,571) Z, it, posn(:,igr), Taux_oc(igr)
          elseif( Temperature ) then
            write(3,572) Z, it, posn(:,igr), Temp_coef(igr)
          elseif( Atom_nonsph ) then
            write(3,580) Z, it, posn(:,igr), norbv(igr), (popats(igr,l,1:nspin), l = 1,nlat( abs(it)) )
            if( norbv(igr) /= 0 ) write(3,600) (hybrid(io,:,igr), pop_nonsph(io,igr), io, io = 1,norbv(igr))
          else
            write(3,570) Z, it, posn(:,igr)
          endif
          if( Atom_occ_hubb .and. Hubb( abs(it) ) ) then
            l = l_hubbard( Z )
            write(3,610) ( ( occ_hubb_e(m,m,isp,igr), m = -l,l ), isp = 1,nspin )
          endif
        end do
        write(3,*)

! About potental
        if( Flapw ) then
          if( Hedin .or. Perdew ) then
            write(3,'(A)') ' FLAPW potential energy dependant'
          else
            write(3,'(A)') ' FLAPW potential not energy dependant'
          endif
        else
          if( hedin ) then
            write(3,'(A)') ' Hedin and Lundqvist exchange-correlation potential'
          elseif( perdew ) then
            write(3,'(A)') ' Perdew and Zunger exchange-correlation potential'
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
        if( abs( rpotmax ) > eps4 ) write(3,703) rpotmax
        write(3,704) D_max_pot
      endif

      if( Hubbard ) then
        write(3,'(/A)') ' Hubbard calculation '
        write(3,'(/A)') ' Type  Z  Hubbard parameter (eV)'
        do it = 1,ntype
          if( Hubb(it) ) write(3,'(2i4,f12.3)') it, numat(it), V_hubbard(it)
        end do
      endif
      if( ATA ) write(3,'(/A)') ' Average T-matrix Approximation'

      if( SCF .or. Fermi_auto ) then
        if( Fermi_auto .and. abs(p_self0) < eps10 .and. nself == 1 ) then
          write(3,'(/A)')' One cycle for the Fermi level calculation'
        else
          write(3,'(/A)') ' Self consistent calculation'
          write(3,708) nself, p_self0, p_self_max, Delta_En_conv
        endif
        write(3,709) r_self
        if( scf_elecabs ) write(3,'(A)') '   Cuting energy criterium using the number of electron in the absorbing atom'
        if( self_nonexc ) then
          write(3,'(A)') '   Non excited absorbing atom in this part'
        else
          write(3,'(A)') '   Excited absorbing atom in this part'
        endif
      end if

    endif

    do ipr = 3,6,3
      if( icheck(1) < 1 .and. ipr == 3 ) cycle
      if( Green_s .or. Extract ) then
        if( mpinodes0 > 1 ) then
          write(ipr,'(A)') ' Parallel calculation'
          write(ipr,720) mpinodes0
        else
          write(ipr,'(A)') ' Sequential calculation'
        endif
      else
        if( Solver == 'MUMPS' ) then
          if( mpinodes0 > 1 ) then
            write(ipr,'(A)') ' Parallel calculation with MUMPS Solver'
            write(ipr,710) mpinodes0, mpinodes, MPI_host_num_for_mumps
          else
            write(ipr,'(A)') ' Sequential calculation with MUMPS Solver'
          endif
        else
          if( mpinodes0 > 1 ) then
            write(ipr,'(A)') ' Parallel calculation with Gaussian Solver'
            write(ipr,720) mpinodes0
          else
            write(ipr,'(A)') ' Sequential calculation with Gaussian Solver'
          endif
        endif
      endif
    end do

  endif  ! Point d'arrivee mpirank0 /= 0

  if( mpinodes0 > 1 ) then
    call MPI_Bcast(adimp,n_adimp,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(alfpot,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(allsite,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    if( npldafs > 0 ) call MPI_Bcast(angpoldafs,3*npldafs, MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(angxyz,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ang_base_loc,3*(ntype+1),MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ang_base_loc_gr,3*ngroup,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ang_rotsup,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ang_spin,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Axe_spin,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(axyz,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Basereel,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Base_spin,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(cartesian_tensor,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    if( norbdil > 0 ) call MPI_Bcast(cdil,norbdil,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Centre,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Centre_auto,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Centre_auto_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Clementi,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Core_resolved,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(D_max_pot,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(dyn_eg,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(dyn_g,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Bormann,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Charge_free,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Dafs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Dafs_bio,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Delta_Epsii,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Delta_En_conv,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Density,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Density_comp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Dipmag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
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
    call MPI_Bcast(Ecrantage,nspin,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    if( neimagent > 0 ) call MPI_Bcast(Eeient,neimagent,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Egamme,ngamme,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    if( neimagent > 0 ) call MPI_Bcast(eimagent,neimagent,MPI_REAL8, 0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Eneg_i,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Eneg_n_i,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Energphot,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(f_no_res,2,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(FDM_comp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Force_ecr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Full_atom_e,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Full_potential,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Gamma_tddft,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Green_int,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Green_s,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Green_self,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Hedin,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    if( npldafs > 0 ) call MPI_Bcast(hkl_dafs,3*npldafs, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if( ngroup_nonsph > 0 ) then
      allocate( Hybrid_r(nhybm,16,ngroup_nonsph) )
      allocate( Hybrid_i(nhybm,16,ngroup_nonsph) )
      if( mpirank0 == 0 ) then
        hybrid_r(:,:,:) = real( hybrid(:,:,:),db )
        hybrid_i(:,:,:) = aimag( hybrid(:,:,:) )
      endif
      call MPI_Bcast(Hybrid_r, nhybm*16*ngroup_nonsph,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Hybrid_i, nhybm*16*ngroup_nonsph,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
      if( mpirank0 /= 0 ) Hybrid(:,:,:) = cmplx( Hybrid_r(:,:,:), Hybrid_i(:,:,:),db )
      deallocate( Hybrid_i, Hybrid_r )
    endif
    call MPI_Bcast(iabsm,n_multi_run_e,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(iabsorig,n_multi_run_e,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(iord,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if( npldafs > 0 ) call MPI_Bcast(isigpi,npldafs*2, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if( norbdil > 0 ) call MPI_Bcast(itdil,norbdil, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(itype,ngroup,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(jseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Kern_fac,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    if( ngroup_pdb > 0 ) call MPI_Bcast(Kgroup,ngroup_pdb, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(korigimp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lamstdens,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if( norbdil > 0 ) call MPI_Bcast(ldil,norbdil,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(ldipimp,3,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lecrantage,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lin_gam,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmax_nrixs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmax_pot,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmaxat0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmaxfree,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmaxso0,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lmoins1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lplus1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lquaimp,9,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(lseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if( nlatm > 0 ) call MPI_Bcast(lvval,(ntype+1)*nlatm, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(M1M1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(M1M2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(M2M2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(matper,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(muffintin,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(multrmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(n_atom_proto,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(n_devide,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(nbseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nchemin,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(necrantage,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(No_solsing,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(noncentre,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nonexc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    if( ngroup_nonsph > 0 ) call MPI_Bcast(norbv,ngroup_nonsph+1, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(normaltau,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(normrmt,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nphim,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nself,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nseuil,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nsymextract,n_multi_run_e,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(nposextract,n_multi_run_e,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
    if( nlatm > 0 ) call MPI_Bcast(nvval,(ntype+1)*nlatm, MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(numat_abs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if( Atom_occ_hubb ) call MPI_Bcast(occ_hubb_e,nspin*ngroup_hubb *(2*m_hubb_e+1)**2,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Octupole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Old_reference,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(One_run,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Optic,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Overad,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Overlap,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(p_self_max,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(p_self0,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    if( nple > 0 ) then
      call MPI_Bcast(pdpolar,nple*2,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(polar,3*nple,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    endif
    call MPI_Bcast(polarise,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    if( npldafs > 0 ) then
      call MPI_Bcast(poldafse,3*npldafs,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(poldafsei,3*npldafs,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(poldafss,3*npldafs,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(poldafssi,3*npldafs,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(vecdafsem,3*npldafs,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(vecdafssm,3*npldafs,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    endif
    if( ngroup_nonsph > 0 ) call MPI_Bcast(pop_nonsph, nhybm*ngroup_nonsph,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    if( nlatm > 0 ) then
      call MPI_Bcast(popats,ngroup*nlatm*nspin,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(popval,(ntype+1)*nlatm*nspin,MPI_REAL8, 0,MPI_COMM_WORLD,mpierr)
    endif
    call MPI_Bcast(posn,3*ngroup,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Quadmag,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Quadrupole,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(r_self,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Recup_tddft_data,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Recup_optic,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Relativiste,1,MPI_LOGICAL,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(roverad,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(rpotmax,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(R_rydb,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Rsorte_s,n_radius,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Rydberg,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Save_optic,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Save_tddft_data,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(SCF,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(SCF_elecabs,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Self_cons,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Self_nonexc,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    do i = 1,3
      if( mpirank0 == 0 ) j = iachar( seuil(i:i) )
      call MPI_Bcast(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
      if( mpirank0 /= 0 ) seuil(i:i) = achar( j )
    end do
    call MPI_Bcast(Solsing_s,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Solsing_only,1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Spherical_tensor,1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Spherical_signal,1, MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Spinorbite,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(State_all,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(State_all_out,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Supermuf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    do i = 1,8
      if( mpirank0 == 0 ) j = iachar( PointGroup(i:i) )
      call MPI_Bcast(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
      if( mpirank0 /= 0 ) PointGroup(i:i) = achar( j )
    end do
    call MPI_Bcast(PointGroup_Auto,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Symauto,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Symmol,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Tddft,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Atomic_scr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Test_dist_min,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(q_nrixs,nq_nrixs,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Rpalf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Temp,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Taux_oc,ngroup_taux,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Vec_orig,3,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(V_intmax,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(V_hubbard,ntype+1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(V0bdcFimp,nspin,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    if( nple > 0 ) call MPI_Bcast(Veconde,3*nple,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Ylm_comp_inp,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Z_nospinorbite,1,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)

    n = ntype + 1
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
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
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if( flapw ) then
      call MPI_Bcast(its_lapw,ngroup_lapw,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(nrato_lapw,n,MPI_INTEGER,0, MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(r0_lapw,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(rlapw,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(rotloc_lapw,9*ngroup_lapw,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
      do i = 1,132
        if( mpirank0 == 0 ) j = iachar( Wien_file(8)(i:i) )
        call MPI_Bcast(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
        if( mpirank0 /= 0 ) Wien_file(8)(i:i) = achar( j )
      end do
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if( mpirank0 /= 0 ) then
      if( Extract ) return
      icheck(:) = 0
    endif

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

    if( One_run .and. Taux ) then
      if( istop == 0 ) call write_error
      do ipr = ipr0,9,3
        write(ipr,100)
        write(ipr,'(A/)') ' It is not possible to have at the same time One_run and Crystal_t keywords !'
      end do
      istop = 1
    endif

    if( istop == 1 ) stop
  endif

! Conversion en unites atomiques.
  rad = pi / 180

  adimp(:) = adimp(:) / bohr

  where( abs(angpoldafs) < 9999._db ) angpoldafs = angpoldafs * rad
  axyz(1:3) = axyz(1:3) / bohr

  D_max_pot = D_max_pot / bohr
  Delta_En_conv = Delta_En_conv / Rydb
  Delta_Epsii = Delta_Epsii / Rydb
  dpos(:) = dpos(:) / bohr
  E_adimp(:) = E_adimp(:) / Rydb
  E_max_range(:) = E_max_range(:) / Rydb
  E_radius(:) = E_radius(:) / Rydb
  Eclie = Eclie / Rydb
  Eclie_out = Eclie_out / Rydb
  eeient(:) = eeient(:) / Rydb
  egamme(:) = egamme(:) / Rydb
  eimagent(:) = eimagent(:) / Rydb
  if( Hubbard ) V_hubbard(:) = V_hubbard(:) / Rydb
  Pas_SCF = Pas_SCF / rydb
  q_nrixs(:) = q_nrixs(:) * bohr  ! q_nrxis en A^-1 en entree
  r_self = r_self / bohr
  rchimp(:) = rchimp(:) / bohr
  rmtimp(:) = rmtimp(:) / bohr
  roverad = roverad / bohr
  rpotmax = rpotmax / bohr
  R_rydb = R_rydb / bohr
  rsorte_s(:) = rsorte_s(:) / bohr
  Test_dist_min = Test_dist_min / bohr
  v_intmax = v_intmax / Rydb
  v0bdcFimp(:) = v0bdcFimp(:) / Rydb

! Nombre d'orbitales atomiques maximum
  nnlm = 0
  do it = 1,ntype
    n = n_orb_rel( numat(it) )
    if( numat(it) == numat_abs )then
      n = n + 2 * max( 1, nlat(it) )
    elseif( nlat(it) > 0 ) then
      n = n + 2 * nlat(it)
    endif
    nnlm = max( nnlm, n )
  end do

  call cal_cubmat(angxyz,cubmat,struct)

  if( Magnetic .or. Atom_nonsph ) then
    call invermat(cubmat,cubmati)

    if( abs( Axe_spin(1) ) < eps10 .and. abs( Axe_spin(2) ) < eps10 .and. abs( Axe_spin(3) - 1 ) < eps10 ) then
      Axe_spin(1:2) = 0._db; Axe_spin(3) = 1._db
      Axe_spin = matmul(cubmat, Axe_spin)
      call mat_euler( Ang_spin, Rot_gen )
      Axe_spin = matmul( Rot_gen, Axe_spin )
    else
      Axe_spin(:) = Axe_spin(:) * axyz(:) * bohr
      Axe_spin = matmul( cubmat, Axe_spin )
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
    Axe_spin = matmul( cubmati, Axe_spin )
    Axe_spin(:) = Axe_spin(:) / axyz(:)

    do igr = 1,ngroup
      if( Ang_base_loc_gr(1,igr) > - 1000._db ) cycle
      if( Ang_base_loc(1, abs(itype(igr)) ) < -1000._db ) cycle
      Ang_base_loc_gr(:,igr) = Ang_base_loc(:, abs(itype(igr)) )
    end do
    do igr = 1,ngroup
      if( Ang_base_loc_gr(1,igr) < - 1000._db ) Ang_base_loc_gr(:,igr) = Ang_spin(:)
      Ang(:) = Ang_base_loc_gr(:,igr)
      call mat_euler( Ang, Rot )
      Rot_Atom_gr(:,:,igr) = Rot(:,:)
      if( struct == 'trigo' ) then
        Axe(1:3) = 1._db / sqrt(3._db)
      else
        Axe(1:2) = 0._db; Axe(3) = 1._db
      endif
      Axe = matmul(cubmat, Axe)
      Axe = matmul(Rot, Axe)
      Axe = matmul( cubmati, Axe )
      Axe_atom_gr(:,igr) = Axe(:) / axyz(:)
      if( itype(igr) < 0 ) Axe_atom_gr(:,igr) = - Axe_atom_gr(:,igr)

      p(:) = Axe_atom_gr(:,igr) * axyz(:)
      p = matmul( cubmat, p )
      vv = sum( p(:)**2 )
      if( vv > eps10 ) p(:) = p(:) / sqrt( vv )
      p = matmul( cubmati, p )
! Axe_atom_gr est maintenant normalise a 1.
      Axe_atom_gr(:,igr) = p(:) / axyz(:)
    end do

    if( icheck(1) > 0 ) then
      write(3,750) Axe_spin(1:3) * axyz(1:3)
      write(3,760) Ang_spin(1:3) * 180 / pi
      do igr = 1,ngroup
        write(3,770) igr
        do i = 1,3
          write(3,780) Rot_atom_gr(i,:,igr), Axe_atom_gr(i,igr) * axyz(i)
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

  if( octupole) then
    l_selec_max = lseuil + 3
  elseif( Quadrupole) then
    l_selec_max = lseuil + 2
  else
    l_selec_max = lseuil + 1
  endif

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

  if( Centre_auto ) call Auto_center(axyz,Centre, Centre_auto_abs,cubmat,icheck(1),itype, &
                     ngroup,ntype,numat,numat_abs,posn,struct)

  do igr = 1,ngroup
    posn(:,igr) = posn(:,igr) - Centre(:)
  end do

  Tensor_imp(1:3) = ldipimp(1:3)
  k = 3
  do i = 1,3
    do j = 1,3
      k = k + 1
      Tensor_imp(k) = lquaimp(i,j)
    end do
  end do

  call Rmt_fix(icheck(1),ntype,numat,Rmt)

! Stockage pour transmission plus courte

  Multipole(1) = E1E1; Multipole(2) = E1E2e; Multipole(3) = E1E3;
  Multipole(4) = E1M1; Multipole(5) = E1M2; Multipole(6) = E2E2;
  Multipole(7) = E3E3; Multipole(8) = M1M1; Multipole(9) = M1M2;
  Multipole(10) = M2M2

  SCF_log(1) = SCF_elecabs
  SCF_log(2) = SCF_mag_free
  SCF_log(3) = Self_cons
  SCF_log(4) = Self_nonexc
  SCF_log(5) = SCF

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
  150 format(//' The following line is not understood :',/A,// &
               ' If it is a keyword, check the spelling.'/, &
               ' If the line is not supposed to be a keyword but contains numbers, check:'/ &
               '      - How many numbers must be in the line ?'/, &
               '      - Are there spaces between the numbers ?'//)
  160 format(///' Error under the keyword Par_',a6,/ &
                ' The wanted atom is the number',i3,' !',/ &
                ' There are only',i3,' atoms in the job !'//)
  170 format(//'  A parameter index for the fit is not possible !'/, &
               '  Check your indata file under the keyword ',a9,/ &
               '  The index is',i4//)
  171 format(/' Doping atom is indexed with',i3,' what is more than',/ &
              ' the number of atom type =',i2,' declared under keyword Atom !')
  172 format(/' Atom number',i4,' is indexed with',i3,' what is more than',/ &
              ' the number of atom type =',i2,' declared under keyword Atom !')
  174 format(//' Under keyword Crystal or Molecule, when the index is',/ &
               ' supposed to be the atomic number, the simultaneous',/ &
               ' use of the keyword Atom is forbidden !'//)
  175 format(//' The atomic number given under keyword Z_absorber =',i4,' is not in the list of atoms !'//)
  180 format(/' Edge = ',a3,' not programmed !'//)
  190 format(//'  No transitions are allowed in your calculation !'/, &
               ' Check the keywords Quadrupole, Octupole, No_dipole, No_dipquad, No_E2E2, Dipmag, Quadmag, ...')
  200 format(//'  There is no atom in your calculation !'/,  &
               ' Some necessary keywords as molecule or crystal could be missing'//)
  210 format(//'  There is no chemical species specified in your calculation !'//)
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
  300 format(/' Tensors Extracted from the file :'/,A,/)
  310 format(' Radius =',f6.2)
  315 format(' Radius, E_radius =',100(f6.2,f6.1))
  320 format('    Roverad =',f6.2)
  330 format(' icheck =',30i2)
  340 format(' Linear range :')
  350 format(' Range =',9f8.3,5(/9x,9f8.3))
  400 format(' Spin-orbit not taken into account for the atomic number', i3)
  410 format('    Path expansion, n =',i3)
  420 format('    Imposed type muffin-tin radius, Rmtimp =',10f6.3,/9x,10f6.3)
  430 format('    Overlap of the muffin-tin radius  =',f6.2)
  440 format('    iord =',i2,', adimp =',f6.2)
  445 format('    iord =',i2,', adimp, E_adimp =',100( f6.2, f6.1) )
  450 format('    lmaxso0 =',i3)
  460 format('    R_rydb =',f7.3,' A')
  470 format('      Center =',3f7.3)
  480 format('    Eclie, Eclie_out =',2f7.3,' eV')
  490 format('    V_intmax =',f7.3,' eV')
  500 format(/' Temperature =',f6.1,' K')
  510 format(6x,3f7.3,3x,3f7.3,3x,f7.3,4x,f7.3)
  511 format(/' Bormann',/ '  (h, k, l) = (',i3,',',i3,',',i3,')   Azimuth =',f7.2)
  512 format(6x,3f7.3,3x,3f7.3,3x,A)
  514 format(7x,3i3,3x,a11,2f10.3,'      scan')
  517 format(7x,3i3,3x,a11,'      scan',2f10.3)
  518 format(7x,3i3,3x,a11,f10.3,'          scan',f10.3)
  519 format(7x,3i3,3x,a11,3f10.3)
  542 format('   ngroup =',i5,', ntype =',i2)
  543 format(1x,2i4,10(i4,i3,f7.3))
  544 format(1x,2i4,10(i4,i3,2f7.3))
  545 format('    dpos =',3f7.3)
  548 format(' Orbital dilatation :',/'   it    l   cdil')
  549 format(2i5,f7.3)
  552 format('   Point Group = ',a8)
  554 format('   Space Group = ',a13)
  555 format('   a, b, c =',3f12.7)
  556 format('   alfa, beta, gamma =',3f9.3)
  560 format('    Z  Typ       posx           posy           posz    ',A)
  565 format(i5,i4,3f15.10,i6,f10.5,f11.2)
  570 format(i5,i4,3f15.10,i6,f10.2)
  571 format(i5,i4,3f15.10,f10.5,f11.2)
  572 format(i5,i4,3f15.10,f10.2)
  580 format(i5,i4,3f15.10,i5,2x,12f8.4)
  600 format(4x,17(1x,2f7.3),1x,'= hybrid, pop_nonsph(',i1,')')
  610 format('    Occ. matrix :', 14f5.2)
  660 format(' Xalfa potential , Xalfa =',f8.5)
  670 format(' Full potential inside the atomic spheres with lmax =',i2)
  680 format(' E_imag =',f9.3,' eV')
  690 format(2f9.3)
  702 format('   multrmax =',i2)
  703 format('   Rpotmax =',f7.3,' A')
  704 format('   D_max_pot =',f7.3,' A')
  708 format('   Maximum number of iteration = ',i3,/ '   Weight =',f6.3,',   Weight max =',f6.3/ &
             '   Delta energy for convergence =',f7.3,' eV / atom')
  709 format('   Cluster radius used for this part = ',f7.3,' A')
  710 format('  with',i3,' processors including',i3,' for the energy loop and',i3,' for MUMPS')
  720 format('  with',i3,' processors')
  740 format(//' The mesh parameter',i2,' is zero !'//)
  750 format(/' General Z axis =',3f9.5)
  760 format(' Euler angles   =',3f9.3)
  770 format(/6x,' local matrix rotation',7x,'local Z axis   Atom =',i3)
  780 format(3x,3f9.5,5x,f9.5)

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

  if( icheck > 0 ) write(3,110) Rmt(1:ntype)
  Rmt(:) = Rmt(:) / bohr

  return
  110 format(/' FDM atom radius =',10f6.3)
end

!***********************************************************************

function number_from_text(n_skip,mot_in)

  use declarations
  implicit none

  character(len=126):: mot, mot_in

  integer:: i, j, length, n, n_skip

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
  read(9,*) number_from_text
  Close(9)

  return
end

!***********************************************************************

! Routine de lecture de la structure venant de WIEN

subroutine lect_struct_lapw(angxyz,axyz,icheck,its_lapw,itype,ngroup,ngroup_lapw,nomstruct,nrato_lapw,nrm,nslapwm,ntype, &
        numat,posn,r0_lapw,rlapw,rotloc_lapw,Wien_matsym,Wien_taulap)

  use declarations
  implicit none

  character(len=1) Trans
  character(len=2) Plan
  character(len=132) nomstruct

  integer:: i, ia, icheck, igr, index, ipr, iprot, is, istat, it, itr, its, ittt, j, jatom, mu, mult, ngroup, ngroup_lapw, nrm, &
    nmatsym, nslapwm, nt, ntrans, ntype

  integer, dimension(ngroup) :: numprot
  integer, dimension(0:ntype) :: nrato_lapw, numat
  integer, dimension(ngroup_lapw) :: its_lapw
  integer, dimension(ngroup) :: itype
  integer, dimension(3,3,nslapwm):: Wien_matsym

  real(kind=db):: demi, rZ
  real(kind=db), dimension(3):: angxyz, axyz, dp, p
  real(kind=db), dimension(3,3):: Mas, Mat, ptrans, Rot
  real(kind=db), dimension(3,ngroup):: posn
  real(kind=db), dimension(3,3,ntype):: rotloc
  real(kind=db), dimension(3,3,ngroup_lapw):: rotloc_lapw
  real(kind=db), dimension(0:ntype):: r0_lapw, rlapw
  real(kind=db), dimension(3,nslapwm):: Wien_taulap

! ll             = nombre de (l,m) par atome
! ntype          = nombre d'atomes inequivalents
! nmatsym        = nombre d'op. de symetrie
!
! Unite 8 : 'case.struct'  (donnees structurales)

  open(8, file = nomstruct, status='old', iostat=istat)
  if( istat /= 0 ) call write_open_error(nomstruct,istat,1)

! Lecture de xxxx.struct
! Read et formats pris en partie dans la routine main1.f de
! Wien97/SRC_lapw5

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
  do jatom = 1,ntype

    index = index + 1
    read(8,'(5x,i3,1x,3(3x,f10.7))') its, posn(:,index)
    it = abs( its )
    itype(index) = it
    its_lapw(index) = its
    numprot(index) = index
    iprot = index
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

! Recherche de l'operation de symetrie qui renvoit a l'atome
! prototypique
  nt = ntrans + 1

  boucle_exter: do igr = 1,ngroup
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
    do ia = 1,ngroup
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

subroutine cal_cubmat(angxyz,cubmat,struct)

  use declarations
  implicit none
  include 'mpif.h'

  character(len=5) struct

  logical ang(3), ange(3)

  real(kind=db):: a, alfa, b, beta, cosa, cosb, cosg, gamma, rad, sina, sinb
  real(kind=db), dimension(3):: angxyz(3)
  real(kind=db), dimension(3,3):: cubmat

! Matrice de changement de repere cristallo, cubique
  ang(:) = abs( angxyz(:) - 90._db ) < eps4
  ange(1) = abs( angxyz(2) - angxyz(3) ) < eps4
  ange(2) = abs( angxyz(3) - angxyz(1) ) < eps4
  ange(3) = abs( angxyz(1) - angxyz(2) ) < eps4
  if( ange(1) .and. ange(2) .and. ange(3) ) then
    if( ang(1) ) then
      struct = 'cubic'
    else
      struct = 'trigo'
    endif
  elseif( ( abs( angxyz(3) - 120. ) < eps4 ) .and. ang(1) .and. ang(2) ) then
    struct = 'hexag'
  else
    struct = 'autre'
  endif

  if( struct /= 'cubic' ) then

    rad = pi / 180
    if( struct == 'trigo' ) then
      alfa = angxyz(1) * rad
      cosa = sqrt( ( 1._db + 2 * cos( alfa ) ) / 3. )
      sina = sqrt( 1._db - cosa**2 )
      cubmat(1,1) = sina;  cubmat(1,2:3) = -0.5_db * sina
      cubmat(2,1) = 0._db;  cubmat(2,2) = sqrt(3._db) * sina / 2
      cubmat(2,3) = - cubmat(2,2)
      cubmat(3,1:3) = cosa
    else
      alfa = angxyz(1) * rad
      beta = angxyz(2) * rad
      gamma = angxyz(3) * rad
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

! Calcule la matrice de rotation en fonction des angles d'Euler.

subroutine mat_euler(Ang,Rot)

  use declarations
  implicit none

  integer:: i, j, k, l

  real(kind=db):: Angr, cs, ss
  real(kind=db), dimension(3):: Ang
  real(kind=db), dimension(3,3):: Mat, Rot

  do l = 1,3

    Angr = Ang(4-l)

    cs = cos( Angr )
    ss = sin( Angr )
    i = mod(l,3) + 1
    j = mod(l+1,3) + 1
    k = mod(l-1,3) + 1
    mat(i,i) = cs;    mat(i,j) = -ss;   mat(i,k) = 0._db
    mat(j,i) = ss;    mat(j,j) = cs;    mat(j,k) = 0._db
    mat(k,i) = 0._db;  mat(k,j) = 0._db;  mat(k,k) = 1._db
    if( l == 1 ) then
      rot = mat
    else
      rot = matmul( mat, rot )
    endif
  end do

  return
end

!***********************************************************************

! Calcul du centre de l'agregat.

subroutine Auto_center(axyz,Centre,Centre_auto_abs, cubmat,icheck,itype,ngroup,ntype, numat,numat_abs,posn,struct)

  use declarations
  implicit none

  integer:: ia1, ia2, ia3, ia4, icheck, igr, jgr, ngr, ngroup, ntype, numat_abs, Z
  integer, dimension(0:ntype):: numat
  integer, dimension(ngroup):: itype

  character(len=2):: Chemical_Symbol
  character(len=5) struct

  logical:: Centre_auto_abs

  real(kind=db):: dist, dist_max, Radius
  real(kind=db), dimension(3):: axyz, b, Centre, p, q, v
  real(kind=db), dimension(3,3):: Cubmat, Mat, Mati
  real(kind=db), dimension(3,ngroup):: posn, pos

  dist_max = 0._db

  jgr = 0
  do igr = 1,ngroup
    Z = numat( abs( itype(igr) ) )
    if( Centre_auto_abs .and. Z /= numat_abs ) cycle
    jgr = jgr + 1
    if( struct /= 'cubic' ) then
      p(:) = posn(:,igr)
      p = matmul( Cubmat, p )
      pos(:,jgr) = p(:) * axyz(:)
    else
      pos(:,jgr) = posn(:,igr) * axyz(:)
    endif
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
      write(3,120) Chemical_Symbol(numat_abs), Radius * bohr
    else
      write(3,130) Radius * bohr
    endif
  endif

  return
  110 format(/' Center set at  =',3f9.5,'  (cell unit)')
  120 format(' Farest ',a2,' atom at =',f9.5,' A')
  130 format(' Farest atom at =',f7.3,' A')
end



