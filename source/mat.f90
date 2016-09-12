! FDMNES subroutines

! Remplissage de la matrice FDM et resolution du systeme d'equations lineaires.

subroutine mat(Adimp,Atom_axe,Axe_atom_grn,Base_hexa,Base_ortho,Basereelt,Cal_xanes,cgrad, &
                    clapl,dcosxyz,distai,E_comp,Ecinetic_out,Eclie_out,Eimag,Eneg,Enervide,FDM_comp,Full_atom, &
                    gradvr,iaabsi,iaprotoi,iato,ibord,icheck,ie,igreq,igroupi,igrph,irep_util,isbord,iso,ispinin,isrt,ivois, &
                    isvois,karact,lato,lmaxa,lmaxso,lso,mato,MPI_host_num_for_mumps,mpirank0,mso, &
                    natome,n_atom_0,n_atom_ind,n_atom_proto,nbm,nbord,nbordf,nbtm,neqm,ngroup_m,ngrph,nim,nicm, &
                    nlmagm,nlmmax,nlmomax,nlmsa,nlmsam,nlmso,nphiato1,nphiato7,npoint,npr,npsom,nsm, &
                    nsort,nsortf,nso1,nspin,nspino,nspinp,nstm,numia,nvois,phiato,poidsa,poidso,R_rydb,Recop,Relativiste, &
                    Repres_comp,Rsort,rvol,Rydberg,Solsing,Spinorbite,State_all_r,Sym_cubic,Tau_ato,Taull, &
                    Time_fill,Time_tria,V0bd_out,Vr,xyz,Ylm_comp,Ylmato,Ylmso)

  use declarations
  implicit none
  
  integer:: i, i1, i2, ia, iaabsi, icheck, idir1, idir2, idir3, ie, igrph, ii, ipas, irep, is1, is2, isens, isg, isp, ispinin, &
    jj, jsp, l, l1, l2, lm, lm01, lm01c, lm02, lm02c, lm1, lm2, lmaxso, lmf, lms, lmw, m, m1, m2, &
    MPI_host_num_for_mumps, mpirank0, natome, n_atom_0, n_atom_ind, n_atom_proto, nbm, nbtm, neqm, ngroup_m, ngrph, nim, nicm,&
    nligne, nligne_i, nligneso, nlmagm, nlmmax ,nlmomax, nlmsam, nlmso, nlmso_i, nlmw, nphiato1, nphiato7, npoint, npr, &
    npsom, nsm, nsort, nsort_c, nsort_r, nsortf ,nso1, nspin, nspino, nspinp, nspinr, nstm, nvois

  integer, dimension(0:npoint):: new
  integer, dimension(natome):: ianew, iaprotoi, igroupi, lmaxa, nbord, nbordf, nlmsa
  integer, dimension(nstm):: isrt
  integer, dimension(npsom):: numia
  integer, dimension(nso1,ngrph):: lso, mso, iso
  integer, dimension(nlmsam,natome,ngrph):: lato, mato, iato
  integer, dimension(nbtm,natome):: ibord, isbord
  integer, dimension(nrepm,2):: irep_util
  integer, dimension(npsom,nvois):: ivois, isvois
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(:), allocatable:: lb1, lb2, newinv

  character(len=1) mot1

  complex(kind=db):: ampl1, ampl2, cfac
  complex(kind=db), dimension(nopsm,nrepm):: karact
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,natome):: taull
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,n_atom_0:n_atom_ind):: Tau_ato
  complex(kind=db), dimension(:), allocatable :: smc
  complex(kind=db), dimension(:,:), allocatable :: norm,norm_t
  complex(kind=db), dimension(:,:,:), allocatable :: Bessel, Neuman
  complex(kind=db), dimension(:,:,:,:), allocatable :: taull_tem

  logical:: Base_hexa, Base_ortho, Basereel, Basereelt, Cal_xanes, Cal_comp, E_comp, Eneg, FDM_comp, Full_atom, recop, &
    Relativiste, Repres_comp, Rydberg, Spinorbite, Solsing, State_all_r, Stop_job, Sym_cubic, Ylm_comp
  logical, dimension(natome):: Atom_axe

  real(kind=db):: adimp, Eclie_out, Eimag, Enervide, R_rydb, Rsort, Time_fill, Time_tria 

  real(kind=db), dimension(3):: dcosxyz
  real(kind=db), dimension(nspin):: Ecinetic_out, V0bd_out
  real(kind=db), dimension(nopsm,nspino):: Kar, Kari
  real(kind=db), dimension(nvois):: cgrad
  real(kind=db), dimension(0:nvois):: clapl
  real(kind=db), dimension(npoint,nspin):: Vr
  real(kind=db), dimension(nbm,natome):: poidsa
  real(kind=db), dimension(nim):: rvol
  real(kind=db), dimension(nsm):: poidso
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(nsort,nlmomax):: Ylmso
  real(kind=db), dimension(natome):: distai
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_grn
  real(kind=db), dimension(nicm,3,nspin):: gradvr
  real(kind=db), dimension(nbtm,nlmmax,natome):: Ylmato
  real(kind=db), dimension(nphiato1,nlmagm,nspinp,nspino,natome,nphiato7):: phiato

  real(kind=db), dimension(:,:), allocatable :: smi, smr
  real(kind=db), dimension(:,:,:), allocatable :: Besselr, Neumanr

  Stop_job = .true.

  do isp = 1,nspino
    irep = abs( irep_util(igrph,isp) )
! Il n'y a pas de couplage up-down dans le cas suivant
    if( irep == 0 ) irep = abs( irep_util(igrph,3-isp) )
    Kar(:,isp) = real( karact(:,irep), db )
    if( Ylm_comp ) then
      if( irep_util(igrph,isp) > 0 ) then
        Kari(:,isp) = aimag( karact(:,irep) )
      else
        Kari(:,isp) = - aimag( karact(:,irep) )
      endif
    else
      Kari(:,isp) = 0._db
    endif
  end do
  Basereel = Basereelt
  if( E_comp .or. Ylm_comp ) then
    Cal_comp = .true.
    Basereel = .false.
  else
    Cal_comp = .false.
  endif

  if( .not. Basereel ) Cal_comp = .true.

  if( E_comp ) then
    nsort_c = nsort
    nsort_r = 0
  else
    nsort_c = 0
    nsort_r = nsort
  endif
  
  if( Spinorbite ) then
    nspinr = nspin
  else
    nspinr = 1
  endif
  
  allocate( Bessel(nsort_c,0:lmaxso,nspinr) )
  allocate( Neuman(nsort_c,0:lmaxso,nspinr) )
  allocate( Besselr(nsort_r,0:lmaxso,nspinr) )
  allocate( Neumanr(nsort_r,0:lmaxso,nspinr) )

  jsp = 0
  do isp = 1,nspin
    if( .not. Spinorbite .and. isp /= ispinin ) cycle
    jsp = jsp + 1
    call phiso(Adimp,Base_ortho,Bessel,Besselr,dcosxyz,E_comp,Ecinetic_out(isp),Eclie_out,Eimag, &
       Eneg,Enervide,icheck,jsp,isrt,lmaxso,mpirank0,Neuman,Neumanr,npsom,nsort,nsort_c,nsort_r,nspinr,nstm, &
       R_Rydb,Rsort,Rydberg,V0bd_out(isp),xyz)
  end do

  nligne = nspino * npr + nlmso + sum( nlmsa(1:natome) )

  allocate( lb1(nligne) ); allocate( lb2(nligne) )
  allocate( newinv(nligne) )

  call newind(distai,ianew,ibord,icheck,isrt,ivois,lb1,lb2,mpirank0,natome,nbordf,nbtm,new,newinv, &
              nligne,nlmsa,nlmso,npoint,nsortf,nspino,npsom,nstm,numia,nvois,xyz)

  if( Cal_comp ) then
    nligne_i = nligne
    nlmso_i = nlmso
  else
    nligne_i = 0
    nlmso_i = 0
  endif
  allocate( smi(nlmso_i,nligne_i) )
  allocate( smr(nlmso,nligne) )
  smi(:,:) = 0._db
  smr(:,:) = 0._db

  if( icheck > 1 ) write(3,110)

  nligneso = nligne - nlmso

  call mat_solve(Base_hexa, Basereel, Bessel, Besselr, Cal_comp, cgrad, clapl, E_comp, Eimag, Enervide, gradvr, &
        ianew, iato, ibord, icheck, ie, igrph, ii, isbord, iso, ispinin, isrt, isvois, ivois, Kar, Kari, lato, &
        lb1, lb2, lmaxso, lso, mato, MPI_host_num_for_mumps, mpirank0, mso, natome, nbm, nbord, nbordf, nbtm, Neuman, Neumanr, &
        new, newinv, ngrph, nicm, nim, nligne, nligne_i, nligneso, nlmsam,  nlmagm, nlmmax, nlmomax, nlmsa, nlmso, nlmso_i, &
        nphiato1, nphiato7, npoint, npsom, nsm, nso1, nsort, nsort_c, nsort_r, nsortf, nspin, nspino, nspinp, nspinr, nstm, &
        numia, nvois, phiato, poidsa, poidso, Relativiste, Repres_comp, rvol, smi, smr, Spinorbite, Time_fill, Time_tria, Vr, &
        Ylmato, Ylmso)

  deallocate( Bessel, Besselr, Neuman, Neumanr )

  if( icheck > 2 ) then
    write(3,120)
    do ii = 1,nligne
      if( Cal_comp ) then
        write(3,130) ii, ( smr(lms,ii), smi(lms,ii), lms = 1,nlmso )
      else
        write(3,140) ii, smr(1:nlmso,ii)
      endif
    end do
  endif

! Valeur des amplitudes des harmoniques spheriques

  ia = iaabsi

  if( Spinorbite ) then
    is1 = 1
    is2 = 2
  else
    is1 = ispinin
    is2 = ispinin
  endif

  if( Basereel ) then
! Si le calcul a ete fait sur la base reelle, il faut renormaliser.

    allocate( norm(nlmso,nlmso) )
    allocate( norm_t(nlmso,nlmso) )

    do lm = 1,nlmso
      ii = nligneso + lm
      if( Cal_comp ) then
        do lm2 = 1,nlmso
          norm(lm2,lm) = cmplx( -smi(lm2,ii), smr(lm2,ii), db )
        end do
      else
        do lm2 = 1,nlmso
          norm(lm2,lm) = cmplx( 0._db, smr(lm2,ii), db )
        end do
      endif
    end do
    do lm = 1,nlmso
      norm(lm,lm) = 1._db + norm(lm,lm)
    end do

    if( icheck > 2 ) then
      write(3,'(/A)') '   lm  norm(lm,1..nlmso)'
      do lm = 1,nlmso
        write(3,130) lm, norm(lm,:)
      end do
    endif

! Inversion
    call invcomp(nlmso,norm,nlmso,nlmso,0,Stop_job)
! for speedup access to norm
    norm_t = TRANSPOSE(norm)

  endif

! Amplitude en sortie.
  do ia = 1,natome

    if( ia /= iaabsi .and. .not. State_all_r .and. Cal_xanes ) cycle

    do lm1 = 1,nlmsa(ia)
      l1 = lato(lm1,ia,igrph)
      m = mato(lm1,ia,igrph)
      lm01 = l1**2 + l1 + 1 + m
      if( Spinorbite ) then
        is1 = iato(lm1,ia,igrph)
      else
        is1 = ispinin
      endif
      ii = ianew(ia) + lm1

      do lm2 = 1,nlmsa(ia)

        jj = ianew(ia) + lm2

        l2 = lato(lm2,ia,igrph)
        m = mato(lm2,ia,igrph)
        lm02 = l2**2 + l2 + 1 + m
        if( nspino == 2 ) then
          is2 = iato(lm2,ia,igrph)
        else
          is2 = ispinin
        endif

        cfac = ( 0._db, 0._db )

        do lmf = 1,nlmso

          if( Basereel ) then
            if( Cal_comp ) then
              ampl1 = sum( norm_t(1:nlmso,lmf) * cmplx( smr(1:nlmso,ii), smi(1:nlmso,ii), db ))
              ampl2 = sum( norm_t(1:nlmso,lmf) * cmplx( smr(1:nlmso,jj), smi(1:nlmso,jj), db ))
            else
              ampl1 = sum( norm_t(1:nlmso,lmf) * smr(1:nlmso,ii) )
              ampl2 = sum( norm_t(1:nlmso,lmf) * smr(1:nlmso,jj) )
            endif
          else
            ampl1 = cmplx(smr(lmf,ii),smi(lmf,ii),db)
            ampl2 = cmplx(smr(lmf,jj),smi(lmf,jj),db)
          endif

          if( FDM_comp ) then
            cfac = cfac + ampl1 * ampl2
          else
            cfac = cfac + conjg( ampl1 ) * ampl2
          endif
          
        end do

        if( Repres_comp .and. .not. Spinorbite .and. .not. Atom_axe(ia) ) cfac = cfac + Conjg(cfac)

! On multiplie par -img pour que ce soit -aimag(taull) qui soit le XANES
        cfac = - img * cfac

        taull(lm01,is1,lm02,is2,ia) = taull(lm01,is1,lm02,is2,ia) + cfac

      end do
    end do

  end do

! Representation conjugue
  if( Repres_comp .and. .not. Spinorbite ) then
    do ia = 1,natome

      if( .not. Ylm_comp  .or. .not. Atom_axe(ia) ) cycle

      do lm1 = 1,nlmsa(ia)
        l1 = lato(lm1,ia,igrph)
        m1 = mato(lm1,ia,igrph)
        lm01 = l1**2 + l1 + 1 + m1
        lm01c = lm01 - 2 * m1
        is1 = ispinin

        do lm2 = 1,nlmsa(ia)

          l2 = lato(lm2,ia,igrph)
          m2 = mato(lm2,ia,igrph)
          lm02 = l2**2 + l2 + 1 + m2
          lm02c = lm02 - 2 * m2
          is2 = ispinin
          isg = (-1)**(m1+m2)
          taull(lm02c,is1,lm01c,is2,ia) = taull(lm02c,is1,lm01c,is2,ia) + isg * taull(lm01,is1,lm02,is2,ia)
          
        end do
      end do
    end do

  endif

! Si State_all, recop est false.
  if( recop .and. igrph == ngrph .and. ( Spinorbite .or. ispinin == nspin ) ) then
    allocate( taull_tem(nlmagm,nspinp,nlmagm,nspinp) )
    do ia = 1,natome
      if( ia /= iaabsi .and. .not. Atom_axe(ia) ) Cycle
      taull_tem(:,:,:,:) = taull(:,:,:,:,ia)
      call recop_taull(nlmagm,nspinp,Sym_cubic,Taull_tem,Ylm_comp)
      taull(:,:,:,:,ia) = taull_tem(:,:,:,:)
    end do
    deallocate( taull_tem )
  endif

  if( icheck > 1 ) then

    if( Spinorbite ) then
      write(3,'(/A)') ' Wave function, versus (l_out,m_out,s_out):'
    else
      write(3,'(/A)') ' Wave function, versus (l_out,m_out):'
    endif
    select case(mso(1,igrph))
      case(-1,-2,-3)
        idir1 = 1; idir2 = 3; idir3 = 2
        mot1 = 'y'
      case(0)
        idir1 = 1; idir2 = 2; idir3 = 3
        mot1 = 'z'
      case(1,2,3)
        idir1 = 2; idir2 = 3; idir3 = 1
        mot1 = 'x'
    end select
    nlmw = min(nlmso,9)
    if( Cal_comp .or. Basereel ) then
      if( Spinorbite ) then
        write(3,150) mot1, (lso(lm,igrph), mso(lm,igrph), iso(lm,igrph), lm = 1,nlmw )
      else
        write(3,160) mot1, (lso(lm,igrph), mso(lm,igrph), lm = 1,nlmw )
      endif
    else
      if( Spinorbite ) then
        write(3,170) mot1, (lso(lm,igrph), mso(lm,igrph), iso(lm,igrph), lm = 1,nlmw )
      else
        write(3,180) mot1, (lso(lm,igrph), mso(lm,igrph), lm = 1,nlmw )
      endif
    endif
    if( Basereel ) allocate( smc(nlmso) )
   
    do jsp = 1,nspino
      if( Spinorbite ) write(3,'(a6,i2)') ' isp =', jsp
      do isens = 1,2
        if( isens == 1 ) then
          i1 = nligne  
          i2 = 1
          ipas = - 1
        else
          i1 = 1
          i2 = nligne
          ipas = 1
        endif
        do ii = i1,i2,ipas
          i = newinv(ii)
          if( i <= 0 ) cycle
          if( Spinorbite ) then
            isp = mod(ii+1,2) + 1
            if( isp /= jsp ) cycle
            if( irep_util(igrph,isp) == 0 ) cycle
          endif
          if( abs(xyz(idir1,i)) > eps10 .or. abs(xyz(idir2,i)) > eps10 ) cycle
          if( ( isens == 1 .and. xyz(idir3,i) > - eps10 ) .or. ( isens == 2 .and. xyz(idir3,i) < - eps10 ) ) cycle
          if( Cal_comp ) then
            if( Basereel ) then
              do lmf = 1,nlmw
                smc(lmf) = sum( norm(lmf,1:nlmso) * cmplx( smr(1:nlmso,ii), smi(1:nlmso,ii), db ))
              end do
              write(3,190) xyz(idir3,i)*bohr, smc(1:nlmw) / rvol(i)
            else
              write(3,200) xyz(idir3,i)*bohr, ( smr(lms,ii) / rvol(i), smi(lms,ii) / rvol(i), lms = 1,nlmw )
            endif
          else
            if( Basereel ) then
              do lmf = 1,nlmw
                smc(lmf) = sum( norm(lmf,1:nlmso) * smr(1:nlmso,ii) )
              end do
              write(3,190) xyz(idir3,i)*bohr, smc(1:nlmw) / rvol(i)
            else
              write(3,200) xyz(idir3,i)*bohr, smr(1:nlmw,ii) / rvol(i)
            endif
          endif
        end do
      end do
      
    end do
    
    if( Basereel ) deallocate(smc)
  endif

  if( Basereel ) then
    deallocate( norm ) 
    deallocate( norm_t )
  endif

  deallocate( smi, smr )

  deallocate( lb1 ); deallocate( lb2 )
  deallocate( newinv )

  if( Solsing .and. igrph == ngrph ) call soustract_tl(Axe_Atom_grn,Cal_xanes,Full_atom, &
                iaabsi,iaprotoi,igreq,igroupi,ispinin,lmaxa,n_atom_0,n_atom_ind,n_atom_proto,natome,neqm,ngroup_m,nlmagm, &
                nspin,nspino,nspinp,Spinorbite,State_all_r,Tau_ato,Taull)

  if( icheck > 1 .and. igrph == ngrph .and. ( Spinorbite .or. ispinin == nspin ) ) then
    do ia = 1,natome
      if( ia /= iaabsi .and. .not. State_all_r .and. Cal_xanes ) cycle
      if( ia == iaabsi ) then
        write(3,210)
      else
        write(3,220) ia
      endif
      lmw = min(3,lmaxa(ia))
      nlmw = min(nlmagm,(lmw+1)**2)
      if( nspinp == 2 ) then
        write(3,230) ((( l, m, i, m = -l,l), l = 0,lmw ), i = 1,2)
      else
        write(3,240) (( l, m, m = -l,l), l = 0,lmw )
      endif
      do is1 = 1,nspinp
        lm1 = 0
        do l1 = 0,lmw
          do m1 = -l1,l1
            lm1 = lm1 + 1
            write(3,250) l1, m1, is1, (( Taull(lm1,is1,lm2,is2,ia), lm2 = 1,nlmw), is2 = 1,nspinp )
          end do
        end do
      end do
    end do
  endif

! Ecriture pour extract
  if( icheck > 0 .and. igrph == ngrph .and. ( Spinorbite .or. ispinin == nspin ) .and. Cal_xanes ) then
    if( icheck == 1 ) write(3,110)
    ia = iaabsi
    ! si lmaxa > 10 il faut changer le format d'ecriture
    lmw = min( 10, lmaxa(ia) )
    nlmw = min( nlmagm, (lmw+1)**2 )
    write(3,260)
    if( nspinp == 1 ) then
      write(3,270) lmw, nspinp
      write(3,280) (( l, m, m = -l,l ), l = 0,lmw )
    else
      write(3,290) lmw, nspinp
      write(3,300) ((( l, m, i, m = -l,l ), l = 0,lmw ), i = 1,2)
    endif
    do is1 = 1,nspinp
      do lm1 = 1,nlmw
        write(3,310) ( ( Taull(lm1,is1,lm2,is2,ia), lm2 = 1,nlmw ), is2 = 1,nspinp )
      end do
    end do
  endif

  return
  110 format(/' ---- Mat ---------',100('-'))
  120 format(/'   ii   Solution(i,lmf = 1,nlmso)')
  130 format(i6,1p,250(1x,2e11.3))
  140 format(i6,1p,500e11.3)
  150 format(4x,a1,2x,18(7x,3i3,7x))
  160 format(4x,a1,2x,18(9x,2i3,8x))
  170 format(4x,a1,2x,18(3x,3i3,3x))
  180 format(4x,a1,2x,18(2x,2i3,3x))
  190 format(f7.3,1p,18(1x,2e11.3))
  200 format(f7.3,1p,36e11.3)
  210 format(/' Multiple scattering amplitude, Taull:')
  220 format(/' Scattering amplitude of the atom ia =',i3,' (Taull):')
  230 format('( l, m, s)',18(7x,3i3,6x))
  240 format('( l, m, s)',18(10x,2i3,7x))
  250 format(3i3,2x,1p,18(1x,2e11.3))
  260 format(/' Taull, Multiple scattering amplitude')
  270 format(2i4,' = lmax, nspinp / ( l, m)')
  280 format(242(14x,2i3,11x))
  290 format(2i4,' = lmax, nspinp / ( l, m, s)')
  300 format(242(12x,3i3,10x))
  310 format(1p,242(1x,2e15.7))

end

!***********************************************************************

! Calcul des fonctions d'onde emergeantes en sortie.

subroutine phiso(Adimp,Base_ortho,Bessel,Besselr,dcosxyz,E_comp,Ecinetic_out,Eclie_out,Eimag, &
       Eneg,Enervide,icheck,isp,isrt,lmaxso,mpirank0,Neuman,Neumanr,npsom,nsort,nsort_c,nsort_r,nspinr,nstm, &
       R_Rydb,Rsort,Rydberg,V0bd_out,xyz)

  use declarations
  implicit none

  integer, parameter:: nptrm = 2000
  
  integer:: i, ib, icheck, ipr, ir, isp, j, l, l2, lmaxso, mpirank0, npsom, nr, nsort, nsort_c, nsort_r, nspinr, nstm

  integer, dimension(nstm):: isrt

  complex(kind=db):: bs(0:lmaxso), fnormc, konde, nm(0:lmaxso), z
  complex(kind=db), dimension(nsort_c,0:lmaxso,nspinr):: Bessel, Neuman
  complex(kind=db), dimension(:), allocatable:: f12
  complex(kind=db), dimension(:,:), allocatable:: rbs, rnm, u

  logical:: Base_ortho, E_comp, Eneg, Rydberg

  real(kind=db):: Adimp, cal_norm, clapl, clapl0, deltar, dr, Ecinetic_out, Ec, Eclie_out, ee, Eimag, Enervide, fnorm, konder, &
                  p1, p2, pp, R_rydb, Rmax, rr, Rsort, V0bd_out, vnorme, zr

  real(kind=db), dimension(3):: dcosxyz, p
  real(kind=db), dimension(0:lmaxso):: bsr, nmr
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(nsort_r,0:lmaxso,nspinr):: Besselr, Neumanr
  real(kind=db), dimension(:), allocatable:: f1, f2, f12r, r, v
  real(kind=db), dimension(:,:), allocatable:: rbsr, rnmr, ur

  if( icheck > 1 ) write(3,110) Ecinetic_out * rydb, lmaxso

  if( E_comp ) then
    konde = sqrt( cmplx( Ecinetic_out, Eimag, db ) )
    if( abs( konde ) < eps10 ) konde = cmplx( eps10, 0._db, db )
    konder = real( konde, db )
!    fnorm = cal_norm(Ecinetic_out,eimag)
!    fnormc = cmplx( fnorm, 0._db, db)
    fnormc = sqrt( konde / pi )
    
    fnorm = real( fnormc, db )
  else
    konder = sqrt( Ecinetic_out )
    konde = cmplx( konder, 0._db, db )
    if( konder < eps10 .and. mpirank0 == 0 ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,120)
      end do
      stop
    endif
    fnorm = sqrt( konder / pi )
    fnormc = sqrt( konde / pi )
  endif

  if( Rydberg ) then

    if( abs( Enervide ) < eps10 ) then
      ee = eps10
    else
      ee = Enervide
    endif
    if( ee < eps10 ) then
      Rmax = - 2 / ee + R_rydb
    else
      Rmax = 100. / bohr
    endif

  endif

  if( Rydberg .and. Rmax >= Rsort ) then
    deltar = rmax - rsort + adimp + eps10
    dr = 0.05 / bohr
    nr = nint( deltar / dr ) + 1

    allocate( r(0:nr) )
    allocate( v(0:nr) )
    allocate( f1(0:nr) )
    allocate( f2(0:nr) )

    dr = deltar / ( nr - 1 )
    do ir = nr,0,-1
      r(ir) = rmax - ( ir - 1 ) * dr
    end do
    v(0:nr) = - 2 / r(0:nr)
    do ir = nr,0,-1
      v(ir) = max( v(ir), V0bd_out )
    end do

! Pour eviter la discontinuite
    if( ee >= eps10 ) then
      pp = - v(0) / ( r(0) - r(nr) )
      do ir = 0,nr
        v(ir) = v(ir) + ( r(ir) - r(nr) ) * pp
      end do
    endif

    if( E_comp ) then
      allocate( f12(0:nr) )
      allocate( u(-1:nr,2) )
      allocate( rnm(2,0:lmaxso) )
      allocate( rbs(2,0:lmaxso) )
    else
      allocate( f12r(0:nr) )
      allocate( ur(-1:nr,2) )
      allocate( rnmr(2,0:lmaxso) )
      allocate( rbsr(2,0:lmaxso) )
    endif

    Ec = Enervide - v(0)
    if( .not. Eneg ) Ec = max( Ec, Eclie_out )
    if( E_comp ) then
      konde = sqrt( cmplx( Ec, eimag, db ) )
      konder = real( konde, db )
    else
      konder = sqrt( ec )
      konde = cmplx( konder, 0._db, db )
    endif
    if( konder < eps10 .and. mpirank0 == 0 ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,120)
      end do
      stop
    endif
    if( E_comp ) then
      fnorm = cal_norm(Ecinetic_out,eimag)
      fnormc = cmplx( fnorm, 0._db, db)
    else
      fnorm = sqrt( konder / pi )
      fnormc = sqrt( konde / pi )
    endif

    do i = 1,2
      rr = rmax + i * dr
      if( E_comp ) then
        z = konde * rr
        call cbessneu(fnormc,z,lmaxso,lmaxso,bs,nm)
        rnm(i,0:lmaxso) = nm(0:lmaxso) * rr
        rbs(i,0:lmaxso) = bs(0:lmaxso) * rr
      else
        zr = konder * rr
        call cbessneur(fnorm,zr,lmaxso,lmaxso,bsr,nmr)
        rnmr(i,0:lmaxso) = nmr(0:lmaxso) * rr
        rbsr(i,0:lmaxso) = bsr(0:lmaxso) * rr
      endif
    end do

    clapl = 1 / dr**2
    clapl0 = - 2 * clapl
    f2(0:nr) = 1 / r(0:nr)**2
    f1(0:nr) = - clapl0 + v(0:nr) - Enervide

    do l = 0,lmaxso
      l2 = l**2 + l

      if( E_comp ) then
        f12(0:nr) = cmplx( f1(0:nr) + l2 * f2(0:nr), -eimag, db ) / clapl
        u(0,1) = rnm(1,l)
        u(-1,1) = rnm(2,l)
        u(0,2) = rbs(1,l)
        u(-1,2) = rbs(2,l)
        do j = 1,2
          do ir = 0,nr-1
            u(ir+1,j) = f12(ir) * u(ir,j) - u(ir-1,j)
          end do
          u(0:nr,j) = u(0:nr,j) / r(0:nr)
        end do
      else
        f12r(0:nr) = ( f1(0:nr) + l2 * f2(0:nr) ) / clapl
        ur(0,1) = rnmr(1,l)
        ur(-1,1) = rnmr(2,l)
        ur(0,2) = rbsr(1,l)
        ur(-1,2) = rbsr(2,l)
        do j = 1,2
          do ir = 0,nr-1
            ur(ir+1,j) = f12r(ir) * ur(ir,j) - ur(ir-1,j)
          end do
          ur(0:nr,j) = ur(0:nr,j) / r(0:nr)
        end do
      endif

      do ib = 1,nsort
        i = isrt(ib)
        p(1:3) = xyz(1:3,i)
        rr = vnorme(Base_ortho,dcosxyz,p)
        do ir = nr-1,2,-1
          if( r(ir) > rr ) exit
        end do
        p1 = ( rr - r(ir) ) / ( r(ir+1) - r(ir) )
        p2 = 1 - p1

        if( E_comp ) then
          Neuman(ib,l,isp) = p1 * u(ir+1,1) + p2 * u(ir,1)
          Bessel(ib,l,isp) = p1 * u(ir+1,2) + p2 * u(ir,2)
        else
          Neumanr(ib,l,isp) = p1 * ur(ir+1,1) + p2 * ur(ir,1)
          Besselr(ib,l,isp) = p1 * ur(ir+1,2) + p2 * ur(ir,2)
        endif

      end do

      if( icheck > 2 ) then
        write(3,126) l
        do ir = nr,0,-1
          if( E_comp ) then
            write(3,127) r(ir)*bohr, v(ir)*rydb, u(ir,1:2)
          else
            write(3,127) r(ir)*bohr, v(ir)*rydb, ur(ir,1:2)
          endif
        end do
      endif

    end do

    deallocate( r )
    deallocate( v )
    deallocate( f1 )
    deallocate( f2 )
    if( E_comp ) then
      deallocate( f12 )
      deallocate( u )
      deallocate( rnm )
      deallocate( rbs )
    else
      deallocate( f12r )
      deallocate( ur )
      deallocate( rnmr )
      deallocate( rbsr )
    endif

  else

! Calcul des fonctions de Bessel et Neuman.
    do ib = 1,nsort
      i = isrt(ib)
      p(1:3) = xyz(1:3,i)
      if( E_comp ) then
        z = konde * vnorme(Base_ortho,dcosxyz,p)
        call cbessneu(fnormc,z,lmaxso,lmaxso,bs,nm)
        Neuman(ib,0:lmaxso,isp) = nm(0:lmaxso)
        Bessel(ib,0:lmaxso,isp) = bs(0:lmaxso)
      else
        zr = konder * vnorme(Base_ortho,dcosxyz,p)
        call cbessneur(fnorm,zr,lmaxso,lmaxso,bsr,nmr)
        Neumanr(ib,0:lmaxso,isp) = nmr(0:lmaxso)
        Besselr(ib,0:lmaxso,isp) = bsr(0:lmaxso)
      endif
    end do

  endif

  if( icheck > 2 ) then
    do l = 0,lmaxso
      if( E_comp ) then
        write(3,130) l, isp, konde
        do ib = 1,nsort
          i = isrt(ib)
          p(1:3) = xyz(1:3,i)
          rr = vnorme(Base_ortho,dcosxyz,p)
          z = konde * rr
          write(3,140) isrt(ib), rr*bohr, Bessel(ib,l,isp), Neuman(ib,l,isp)
        end do
      else
        write(3,150) l, isp, konder
        do ib = 1,nsort
          i = isrt(ib)
          p(1:3) = xyz(1:3,i)
          rr = vnorme(Base_ortho,dcosxyz,p)
          zr = konder * rr
          write(3,140) isrt(ib), rr*bohr, Besselr(ib,l,isp), Neumanr(ib,l,isp)
        end do
      endif
    end do
  endif

  return
  110 format(/' ---- Phiso --------',100('-')//, 20x,' Ecinetic_out =',f10.5,' eV,  lmaxso =',i3)
  120 format(//' The wave vector is zero, what is forbidden !'/)
  126 format(/'    radius      v            bess-coul             neuman-coul     l =',i3)
  127 format(2f10.3,1p,4e12.4)
  130 format(/' isrt     r',13x,'bessel',18x,'neuman',9x, 'l =',i3,', isp =',i2,', konde =',1p,2e12.4)
  140 format(i5,f8.3,1x,1p,4e12.4)
  150 format(/' isrt     r',7x,'bessel      neuman',5x, 'l =',i3,', isp =',i2,', konder =',1p,e12.4)
end

!***********************************************************************

! Calculate one row abvr, abvi and smr, smi

subroutine calcMatRow( abvr, abvi, Base_hexa, Basereel, Bessel, Besselr, Cal_comp, cgrad, clapl, E_comp, Eimag, &
    Enervide, gradvr, ianew, iato, ibord, icheck, igrph, ii, isbord, iso, ispin0, isrt, isvois, ivois, Kar, Kari, &
    lato, lb1i, lb1r, lb2i, lb2r, lmaxso, lso, mato, mletl, mso, natome, nbm, nbord, &
    nbordf, nbtm, Neuman, Neumanr, new, newinv, ngrph, nicm, nim, nligne, nligne_i, &
    nligneso, nlmagm, nlmmax, nlmomax, nlmsa, nlmsam, nlmso, nlmso_i, nphiato1, nphiato7, npoint, npsom, nsm, nso1, &
    nsort, nsort_c, nsort_r, nsortf, nspin, nspino, nspinp, nspinr, nstm, &
    numia, nvois, phiato, poidsa, poidso, Relativiste, Repres_comp, rvol, Spinorbite,  &
    smi, smr, Vr, Ylmato, Ylmso )

  use declarations
  implicit none
  
  integer, parameter:: nletm = 52

  integer:: i, ii1, ia, ib, icheck, ii, ispin, ispin0, igrph, indg, isol, isp, ispint, ispo, ispp, ispq, ispt, isym, &
    iv, j, jj, jj1, jjj, k, l, lb1i, lb1r, lb2i, lb2r, ljj, lm, lm0, lmaxso, lmf, lmf0, lmp, lmp0, lms, lp, &
    m, mf, mp, n, natome, nbm, nbtm, ngrph, nicm, nim, nligne, nligne_i, nligneso, nlmagm, nlmso, nlmso_i, & 
    nphiato1, nphiato7, npoint, npsom, nsm, nsort, nsort_c, nsortf, nsort_r, nspino, nspinp, nstm, nso1, nlmsam, nspin, nspinr, &
    nlmmax, nlmomax, nvois
 
  integer, dimension(0:npoint):: new
  integer, dimension(natome):: ianew, nbord, nbordf, nlmsa
  integer, dimension(nstm):: isrt
  integer, dimension(npsom):: numia
  integer, dimension(nso1,ngrph):: lso, mso, iso
  integer, dimension(nlmsam,natome,ngrph):: lato, mato, iato
  integer, dimension(nbtm,natome):: ibord, isbord
  integer, dimension(npsom,nvois):: ivois, isvois
  integer, dimension(nligne):: newinv
  integer, dimension(:), allocatable:: jnzero

  integer, save:: letd, letu, nlet

  character(len=1):: let(0:nletm)
  character(len=2):: mlet(nletm**2)
  character(len=2), dimension(nligne):: mletl

  complex(kind=db) cfac, cfad, charmc, yc, ycomp
  complex(kind=db), dimension(nspino):: cfv_comp
  complex(kind=db), dimension(nsort_c,0:lmaxso,nspino):: Bessel, Neuman

  logical:: Base_hexa, Basereel, Cal_comp, E_comp, Relativiste, Repres_comp, Spinorbite

  real(kind=db):: a2s4, bder, charmr, charmi, crelat, Enervide, Eimag, f, fac, fad, fr, vme
  
  real(kind=db), dimension(2):: Charac_i, Charac_r, csor, csoi, cfvs_r, cfvs_i, v
  real(kind=db), dimension(nopsm,nspino):: Kar, Kari
  real(kind=db), dimension(nvois):: cgrad
  real(kind=db), dimension(0:nvois):: clapl
  real(kind=db), dimension(npoint,nspin):: Vr
  real(kind=db), dimension(nbm,natome):: poidsa
  real(kind=db), dimension(nim):: rvol
  real(kind=db), dimension(nsm):: poidso
  real(kind=db), dimension(nsort,nlmomax):: Ylmso
  real(kind=db), dimension(nicm,3,nspin):: gradvr
  real(kind=db), dimension(nbtm,nlmmax,natome):: Ylmato
  real(kind=db), dimension(nphiato1,nlmagm,nspinp,nspino,natome,nphiato7):: phiato
  real(kind=db), dimension(lb1i:lb2i):: abvi
  real(kind=db), dimension(lb1r:lb2r):: abvr
  real(kind=db), dimension(nlmso,nligne):: smr
  real(kind=db), dimension(nlmso_i,nligne_i):: smi
  real(kind=db), dimension(nsort_r,0:lmaxso,nspino):: Besselr, Neumanr
  
  real(kind=db), save, dimension(nletm**2):: vletr, vleti

  data let/ ' ','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z', &
                'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

  ispin = ispin0
  if( icheck > 2 .and. ii == nligne ) then
    write(3,120)
    nlet = 0
    letu = 0
    letd = 0
  endif

! alfa_sf constante de structure fine.
  a2s4 = 0.25_db * alfa_sf**2

! Remplissage de la ligne et des seconds membres.
  i = newinv(ii)

  if( Spinorbite .and. i > 0 ) then
    do ia = natome,1,-1
      if( ii > ianew(ia) ) exit
    end do
    ii1 = ii - ( ianew(ia) + nlmsa(ia) )
    ispin = 2 - mod(ii1,2)
  endif 
  
! Points du reseau de base
  if( i > 0 ) then

    ispp = min( ispin, nspin )
    indg = 3 - ispin
    ispq = min( indg, nspin )

    isp = min( ispin, nspino )

! Terme diagonal :
    VmE = Vr( i, ispp ) - Enervide
    abvr(ii) = VmE - clapl(0)
    if( E_comp ) abvi(ii) = - Eimag

    if( Relativiste ) abvr(ii) = abvr(ii) - a2s4 * VmE**2
    if( Relativiste .or. Spinorbite ) bder = - a2s4 / ( 1 - a2s4 * VmE )

    do iv = 1,nvois

! Coefficient reliant le point central au voisin
      j = ivois(i,iv)
      ia = numia(j)
      cfvs_r(isp) = - clapl(iv)

      Charac_r(1:nspino) = Kar( abs( isvois(i,iv) ), 1:nspino )
      Charac_i(1:nspino) = Kari( abs( isvois(i,iv) ), 1:nspino )

      if( Spinorbite .or. Relativiste ) then
        if( Base_hexa ) then
          k = 1 + mod((iv-1)/2,4)
        else
          k = 1 + mod((iv-1)/2,3)
        endif
      endif

      if (Relativiste ) then
        if( k == 4 ) then
          crelat = bder * cgrad(iv) * gradvr(i,2,ispp)
        else
          crelat = bder * cgrad(iv) * gradvr(i,k,ispp)
        endif
        cfvs_r(isp) = cfvs_r(isp) + crelat
      endif

      if( Spinorbite ) then
! rotationnel
        csor(:) = 0._db
        csoi(:) = 0._db
        select case(k)
          case(1)
            csoi(isp) = - gradvr(i,2,ispp)
            csor(indg) = gradvr(i,3,ispq)
          case(2,4)
            csoi(isp) = gradvr(i,1,ispp)
            csoi(indg) = - gradvr(i,3,ispq)
          case(3)
            csor(indg) = - gradvr(i,1,ispq) + gradvr(i,2,ispq)
        end select
        csor(:) = bder * cgrad(iv) * csor(:)
        csoi(:) = bder * cgrad(iv) * csoi(:)

        cfvs_r(isp) = cfvs_r(isp) + csor(isp)
        cfvs_r(indg) = csor(indg)
        cfvs_i(:) = csoi(:)
        if( Repres_comp ) then
          v(:) = Charac_r(:) * cfvs_r(:) - Charac_i(:) * cfvs_i(:)
          cfvs_i(:) = Charac_i(:) * cfvs_r(:) + Charac_r(:) * cfvs_i(:)
          cfvs_r(:) = v(:)
        else
          cfvs_i(:) = Charac_r(:) * cfvs_i(:)
          cfvs_r(:) = Charac_r(:) * cfvs_r(:)
        endif
      else
        if( Repres_comp ) cfvs_i(isp) = Charac_i(isp)*cfvs_r(isp)
        cfvs_r(isp) = Charac_r(isp) * cfvs_r(isp)
      endif

      cfvs_r(1:nspino) = cfvs_r(1:nspino) * rvol(i)
      if( Spinorbite .or. Repres_comp ) cfvs_i(1:nspino) = cfvs_i(1:nspino) * rvol(i)

! Le voisin est dans le reseau de base :
      if( ia == 0 ) then

        jj = new(j) - nspino + isp
        do ispo = 1,nspino
          jj1 = jj + ispo - isp
          abvr( jj1 ) = abvr( jj1 ) + cfvs_r(ispo) / rvol(j)
          if( Repres_comp .or. Spinorbite ) abvi( jj1 ) = abvi( jj1 ) + cfvs_i(ispo) / rvol(j)
        end do

      elseif( ia > 0 ) then

! Le voisin est dans l'atome ia :
        do ib = nbordf(ia)+1,nbord(ia)
          if( j == ibord(ib,ia) ) exit
        end do
        isym = abs( isbord(ib,ia) )
        if( isym /= 1 ) then
          do ispo = 1,nspino
            if( Repres_comp ) then
              cfac = cmplx( Kar(isym,ispo), Kari(isym,ispo), db ) * cmplx( cfvs_r(ispo), cfvs_i(ispo), db )
              cfvs_r(ispo) = real( cfac, db)
              cfvs_i(ispo) = aimag( cfac )
            elseif( Spinorbite ) then
              cfvs_r(ispo) = Kar(isym,ispo) * cfvs_r(ispo)
              cfvs_i(ispo) = Kar(isym,ispo) * cfvs_i(ispo)
            else
              cfvs_r(ispo) = Kar(isym,ispo) * cfvs_r(ispo)
            endif
          end do
        endif

        do lm = 1,nlmsa(ia)
          jjj = ianew(ia) + lm
          l = lato(lm,ia,igrph)
          isol = iato(lm,ia,igrph)
          do ispo = 1,nspino       ! c'est en fait le spin
            m = mato(lm,ia,igrph) + ispo - isol  ! m reel
            if( m > l .or. m < -l ) cycle
            lm0 = l**2 + l + 1 + m
            if( Spinorbite ) then
              cfac = cmplx( cfvs_r(ispo), cfvs_i(ispo), db ) * cmplx( phiato(ib,lm0,ispo,isol,ia,1), &
                            phiato(ib,lm0,ispo,isol,ia,2), db )
              abvr( jjj ) = abvr( jjj ) + Real( cfac, db )
              abvi( jjj ) = abvi( jjj ) + aimag( cfac )
            elseif( nphiato7 == 2 ) then
              cfac = cmplx( cfvs_r(isp), cfvs_i(isp), db ) * cmplx( phiato(ib,lm0,ispp,isol,ia,1), &
                            phiato(ib,lm0,ispp,isol,ia,2), db )
              abvr( jjj ) = abvr( jjj ) + Real( cfac, db )
              abvi( jjj ) = abvi( jjj ) + aimag( cfac )
            else
              abvr( jjj ) = abvr( jjj ) + cfvs_r(isp) * phiato(ib,lm0,ispp,isol,ia,1)
            endif
          end do

        end do

! Le voisin est a l'exterieur de la sphere :
      elseif( ia == -2 ) then

        do ib = nsortf+1,nsort
          if( j == isrt(ib) ) exit
        end do

        do lm = 1,nlmso
          l = lso(lm,igrph)
          m = mso(lm,igrph)
          lm0 = l**2 + l + 1 + m
          ispq = min( iso(lm,igrph), nspinr )
          jjj = nligneso + lm

          if( Cal_comp ) then
            if( Repres_comp ) then
              ycomp = yc(m,Ylmso(ib,lm0),Ylmso(ib,lm0-2*m))
            else
              ycomp = cmplx( Ylmso(ib,lm0), 0._db, db )
            endif
            if( Basereel ) then
              if( E_comp ) then
                cfac = Neuman(ib,l,ispq)
              else
                cfac = Neumanr(ib,l,ispq)
              endif
            else
              if( E_comp ) then
                cfac = Bessel(ib,l,ispq) + img*Neuman(ib,l,ispq)
              else
                cfac = cmplx( Besselr(ib,l,ispq), Neumanr(ib,l,ispq), db )
              endif
            endif
            cfv_comp(1:nspino) = cfac * ycomp * cmplx( cfvs_r(1:nspino), cfvs_i(1:nspino),db)
            if( isp == iso(lm,igrph) ) then
              abvr( jjj ) = abvr( jjj ) + real( cfv_comp(isp), db)
              abvi( jjj ) = abvi( jjj ) + aimag( cfv_comp(isp) )
            else
      ! on est forcement en spinorbite
              abvr( jjj ) = abvr( jjj ) + real( cfv_comp(indg),db)
              abvi( jjj ) = abvi( jjj ) + aimag( cfv_comp(indg) )
            endif
          else
            fr = Ylmso(ib,lm0) * Neumanr(ib,l,ispq)
            abvr( jjj ) = abvr( jjj ) + cfvs_r(isp) * fr
          endif

        end do

        do lm = 1,nlmso
          l = lso(lm,igrph)
          m = mso(lm,igrph)
          lm0 = l**2 + l + 1 + m
          ispq = min( iso(lm,igrph), nspinr )

          if( Cal_comp ) then
            if( Repres_comp ) then
              Ycomp = yc(m,Ylmso(ib,lm0),Ylmso(ib,lm0-2*m))
            else
              Ycomp = cmplx( Ylmso(ib,lm0), 0._db, db )
            endif
            if( E_comp ) then
              cfv_comp(1:nspino) = Ycomp * Bessel(ib,l,ispq) * cmplx( cfvs_r(1:nspino), cfvs_i(1:nspino), db )
            else
              cfv_comp(1:nspino) = Ycomp * Besselr(ib,l,ispq) * cmplx( cfvs_r(1:nspino), cfvs_i(1:nspino), db )
            endif
            if( isp == iso(lm,igrph) ) then
              smr(lm,ii) = smr(lm,ii) - real( cfv_comp(isp), db )
              smi(lm,ii) = smi(lm,ii) - aimag( cfv_comp(isp) )
            else
              smr(lm,ii) = smr(lm,ii) - real( cfv_comp(indg), db)
              smi(lm,ii) = smi(lm,ii) - aimag( cfv_comp(indg) )
            endif
          else
            fr = Ylmso(ib,lm0) * Besselr(ib,l,ispq)
            smr(lm,ii) = smr(lm,ii) - cfvs_r(isp) * fr
          endif
        end do

      endif

    end do ! fin de la boucle sur les voisins

! Developpement en sortie
  elseif( i == 0 ) then

    lm = ii - nligneso
    m = mso(lm,igrph)
    lm0 = lso(lm,igrph)**2 + lso(lm,igrph) + 1 + m
    isp = iso(lm,igrph)

    do lmp = 1,nlmso

      if( iso(lmp,igrph) /= isp ) cycle
      l = lso(lmp,igrph)
      mp = mso(lmp,igrph)
      lmp0 = l**2 + l + 1 + mp
      ispq = min( iso(lmp,igrph), nspinr )
      ljj = nligneso + lmp
      if( Cal_comp ) then
        charmc = ( 0._db, 0._db )
      else
        charmr = 0._db
      endif

      do ib = 1,nsortf
        if( Cal_comp ) then
          if( Repres_comp ) then
            cfac = conjg( yc(m,Ylmso(ib,lm0),Ylmso(ib,lm0-2*m)) ) * yc(mp,Ylmso(ib,lmp0), Ylmso(ib,lmp0-2*mp) )
          else
            cfac = cmplx(Ylmso(ib,lm0) * Ylmso(ib,lmp0), 0._db, db)
          endif
          cfac = poidso(ib) * cfac
          if( Basereel ) then
            if( E_comp ) then
              charmc = charmc + cfac * Neuman(ib,l,ispq)
            else
              charmc = charmc + cfac * Neumanr(ib,l,ispq)
            endif
          else
            if( E_comp ) then
              charmc = charmc + cfac * ( Bessel(ib,l,ispq) + img * Neuman(ib,l,ispq) )
            else
              charmc = charmc + cfac * cmplx( Besselr(ib,l,ispq), Neumanr(ib,l,ispq), db )
            endif
          endif
        else
          fac = poidso(ib) * Ylmso(ib,lm0) * Ylmso(ib,lmp0)
          charmr = charmr + fac * Neumanr(ib,l,ispq)
        endif
      end do
      if( Cal_comp ) then
        abvr(ljj) = abvr(ljj) + real( charmc, db )
        abvi(ljj) = abvi(ljj) + aimag( charmc )
      else
        abvr(ljj) = abvr(ljj) + charmr
      endif
      
    end do

    do ib = 1,nsortf

      j = isrt(ib)
      jj = new(j) - nspino + isp

      if( Repres_comp ) then
        cfac = poidso(ib) * conjg( yc(m,Ylmso(ib,lm0),Ylmso(ib,lm0-2*m)) )
        abvr( jj ) = abvr( jj ) - real( cfac, db ) / rvol(j)
        abvi( jj ) = abvi( jj ) - aimag( cfac ) / rvol(j)
      else
        fac = poidso(ib) * Ylmso(ib,lm0)
        abvr( jj ) = abvr( jj ) - fac / rvol(j)
      endif

      do lmf = 1,nlmso
        ispq = min( iso(lmf,igrph), nspinr )
        if( isp /= iso(lmf,igrph) ) cycle
        l = lso(lmf,igrph)
        mf = mso(lmf,igrph)
        lmf0 = l**2 + l + 1 + mf
        if( Cal_comp ) then
          if( Repres_comp ) then
            cfad = cfac * Yc(mf,Ylmso(ib,lmf0),Ylmso(ib,lmf0-2*mf))
          else
            cfad = fac * cmplx(Ylmso(ib,lmf0),0._db,db)
          endif
          if( E_comp ) then
            cfad = cfad * Bessel(ib,l,ispq)
          else
            cfad = cfad * Besselr(ib,l,ispq)
          endif
          smr(lmf,ii) = smr(lmf,ii) - real( cfad, db )
          smi(lmf,ii) = smi(lmf,ii) - aimag( cfad )
        else
          fad = fac * Ylmso(ib,lmf0)
          smr(lmf,ii) = smr(lmf,ii) - fad * Besselr(ib,l,ispq)
        endif

      end do
    end do

! Developpement dans un atome
! A ce niveau isp devient l'indice solution
  else
    ia = -i
    lm = ii - ianew(ia)
    l = lato(lm,ia,igrph)
    ispt = iato(lm,ia,igrph)
    m = mato(lm,ia,igrph)       ! on developpe sur le spin = isol
    lm0 = l**2 + l + 1 + m      ! c'est le m de l'harmonique

    do ib = 1,nbordf(ia)
      j = ibord(ib,ia)
      isym = abs( isbord(ib,ia) )
      jj = new(j) - nspino + ispt
      if( Cal_comp ) then
        if( m /= 0 ) then
          ycomp = yc(m,Ylmato(ib,lm0,ia),Ylmato(ib,lm0-2*m,ia))
          if( abs( Kari(isym,ispt) ) < eps10 ) then
            cfac = poidsa(ib,ia) * conjg( ycomp ) * Kar(isym,ispt) / rvol(j)
          else
            cfac = poidsa(ib,ia) * conjg( ycomp ) * cmplx( Kar(isym,ispt), -Kari(isym,ispt),db) / rvol(j)
          endif
        else
          cfac = poidsa(ib,ia) * Ylmato(ib,lm0,ia) * cmplx( Kar(isym,ispt), -Kari(isym,ispt),db) / rvol(j)
        endif
        abvr( jj ) = abvr( jj ) - real( cfac, db )
        abvi( jj ) = abvi( jj ) - aimag( cfac )
      else
        f = poidsa(ib,ia) * Kar(isym,ispt) / rvol(j)
        abvr( jj ) = abvr( jj ) - f * Ylmato(ib,lm0,ia)
      endif
    end do

    do lmp = 1,nlmsa(ia)
      lp = lato(lmp,ia,igrph)

      isol = iato(lmp,ia,igrph)
      mp = mato(lmp,ia,igrph) + ispt - isol
      if( mp > lp .or. mp < -lp ) cycle
      lmp0 = lp**2 + lp + 1 + mp

      charmr = 0._db
      if( Cal_comp ) charmi = 0._db
      if( Spinorbite ) then
        ispint = min( ispt, nspinp )  ! en fait nspinp = 2 dans ce cas
      else
        ispint = min( ispin, nspinp )
      endif
      do ib = 1,nbordf(ia)
        if( nphiato7 == 2 ) then
          cfac = conjg( Yc(m,Ylmato(ib,lm0,ia),Ylmato(ib,lm0-2*m,ia)) ) * cmplx( phiato(ib,lmp0,ispint,isol,ia,1), &
                        phiato(ib,lmp0,ispint,isol,ia,2), db )
          charmi = charmi + poidsa(ib,ia) * aimag( cfac )
          charmr = charmr + poidsa(ib,ia) * real( cfac, db )
        else
          fr = poidsa(ib,ia) * Ylmato(ib,lm0,ia)
          charmr = charmr + poidsa(ib,ia) * Ylmato(ib,lm0,ia) * phiato(ib,lmp0,ispint,isol,ia,1)
        endif
      end do
      ljj = ianew(ia) + lmp
      abvr(ljj) = abvr(ljj) + charmr
      if( Cal_comp ) abvi(ljj) = abvi(ljj) + charmi
    end do
  endif

  if( icheck > 2 ) then
    boucle_jj: do jj = lb1r,lb2r
      j = jj - lb1r + 1
      if( Cal_comp ) then
        if( ( abs(abvr(jj)) < eps6 ) .and. ( abs(abvi(jj)) < eps6 ) ) then
          mletl(j) = ' .'
          cycle
        endif
      else
        if( abs(abvr(jj))  <  eps6 ) then
          mletl(j) = ' .'
          cycle
        endif
      endif
      if( abs( abvr(jj) - 1 ) < eps6 ) then
        mletl(j) = ' 1'
      elseif( abs(abvr(jj) + 1) < eps6 ) then
        mletl(j) = '-1'
      else
        do k = 1,nlet
          if( Cal_comp ) then
            if( abs(abvr(jj) - vletr(k)) > eps6 .or. abs(abvi(jj) - vleti(k)) > eps6 ) cycle
          else
            if( abs(abvr(jj) - vletr(k)) > eps6 ) cycle
          endif
          mletl(j) = mlet(k)
          cycle boucle_jj
        end do
        if( nlet == nletm**2 ) then
          mletl(j) = ' ?'
        else
          nlet = nlet + 1
          letu = letu + 1
          if( letu > nletm ) then
            letu = 1
            letd = letd + 1
          endif
          mlet(nlet) = let(letd) // let(letu)
          mletl(j) = mlet(nlet)
          vletr(nlet) = abvr(jj)
          if( Cal_comp ) vleti(nlet) = abvi(jj)
        endif
      endif
    end do boucle_jj
    write(3,130) ii, i, lb1r, lb2r
    write(3,140) (mletl(j), j = 1,lb2r - lb1r + 1)
    if( icheck > 3 ) then
      allocate( jnzero(lb2r - lb1r + 1) )
      n = 0
      do j = lb1r,lb2r
        if( Cal_comp ) then
          if( abs( abvr(j) ) +  abs( abvi(j) ) < eps10 ) cycle
         else
          if( abs( abvr(j) ) < eps10 ) cycle
        endif
         n = n + 1
         jnzero(n) = j 
      end do        
      if( Cal_comp ) then
        write(3,'(300(i5,i4,2e12.4))') ( jnzero(j), newinv(jnzero(j)), abvr(jnzero(j)), abvi(jnzero(j)), j = 1,n )
      else
        write(3,'(300(i5,i4,e12.4))'), ( jnzero(j), newinv(jnzero(j)), abvr(jnzero(j)), j = 1,n )
      endif
      deallocate( jnzero )
    endif
    if( sum( abs(smr(:,ii)) ) > eps10 ) then
      write(3,'(A)') ' sm ='
      if( Cal_comp ) then
        write(3,150) ( smr(lms,ii), smi(lms,ii), lms = 1,nlmso )
      else
        write(3,155) ( smr(lms,ii), lms = 1,nlmso )
      endif
    endif
  endif

  if( icheck > 2 .and. ii == 1 ) then
    write(3,*) 
    if( .not. Cal_comp ) then
      write(3,160) (mlet(k), vletr(k), k = 1,nlet)
    else
      write(3,165) (mlet(k), vletr(k), vleti(k), k = 1,nlet)
    endif
  endif
    
  return
  120 format(/'  FDM Matrix :',/'  ii     i   lb1   lb2   / abr(i = lb1,lb2)')
  130 format(4i6)
  140 format(500a2)
  150 format(1p,250(1x,2e11.3))
  155 format(1p,500e11.3)
  160 format(8(a2,' =',1pe10.3,1x))
  165 format(1p,4(a2,' =',2e10.3,1x))
  end

!***********************************************************************

! Fonction donnant l'harmonique complexe en fonction des harmonique reelles

function Yc(m,Ylm1,Ylm2)

  use declarations
  implicit none

  integer:: m
  
  complex(kind=db):: Yc

  real(kind=db):: Ylm1, Ylm2

  if( m == 0 ) then
    Yc = cmplx( Ylm1, 0._db,db)
  else
    if( m < 0 ) then
      Yc = cmplx( Ylm2,-Ylm1,db) * ( (-1)**m ) / sqrt(2._db)
    else
      Yc = cmplx( Ylm1, Ylm2,db) / sqrt(2._db)
    endif
  endif

  return
end

!***********************************************************************

! Calcul du nouvel indicage des points, c'est-a-dire de leur numero de ligne dans la matrice generale.

subroutine newind(distai,ianew,ibord,icheck,isrt,ivois,lb1,lb2,mpirank0,natome,nbordf,nbtm,new,newinv, &
               nligne,nlmsa,nlmso,npoint,nsortf,nspino,npsom,nstm,numia,nvois,xyz)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(natome):: ianew, nlmsa
  integer, dimension(nligne):: lb1, lb2, newinv
  integer, dimension(nstm):: isrt
  integer, dimension(natome):: nbordf
  integer, dimension(npsom):: numia
  integer, dimension(0:npoint):: new
  integer, dimension(nbtm,natome):: ibord
  integer, dimension(npsom,nvois):: ivois

  real(kind=db), dimension(natome):: distai
  real(kind=db), dimension(4,npsom):: xyz

  ii = 0
  do i = 1,npoint
    if( numia(i) == 0 ) then
      ii = ii + nspino
      new(i) = ii
    else
      new(i) = 0
    endif
  end do

  npf = 1
  npf1 = 0
  do ia = 1,natome
    do i = npf,npoint
      if( numia(i) == 0 .and. xyz(4,i) >= distai(ia) ) exit
    end do
    npf = i
    if( npf1 == npf ) then
      ianew(ia) = ianew(ia-1) + nlmsa(ia-1)
    else
      ianew(ia) = new(npf) - nspino
    endif
    npf1 = npf
    do i = npf,npoint
      if( numia(i) == 0 ) new(i) = new(i) + nlmsa(ia)
    end do
  end do

! Calcul des positions des coefficients dans la matrice.
  do ii = 1,nligne
    lb1(ii) = ii
    lb2(ii) = ii
  end do

  n = nspino - 1
  do i = 1,npoint
    if( numia(i) /= 0 ) cycle
    ii = new(i)
    do iv = 1,nvois
      j = ivois(i,iv)
      if( j == 0 .and. mpirank0 == 0 ) then
        call write_error
        do ipr = 3,9,3
          write(ipr,105) i, iv
        end do
        stop
      endif
      ia = numia(j)
      if( ia == 0 ) then
        lb1(ii) = min(lb1(ii),new(j)-n)
        lb2(ii) = max(lb2(ii),new(j))
      elseif( ia > 0 ) then
        lb1(ii) = min(lb1(ii),ianew(ia) + 1)
        lb2(ii) = max(lb2(ii),ianew(ia) + nlmsa(ia))
      elseif( ia == -2 ) then
        lb2(ii) = nligne
      endif
    end do
    if( nspino == 2 ) then
      lb2(ii-1) = lb2(ii)
      lb1(ii-1) = lb1(ii)
    endif
  end do

  do ia = 1,natome
    lb11 = ianew(ia) + 1
    lb22 = ianew(ia) + nlmsa(ia)
    do ib = 1,nbordf(ia)
      lb11 = min( lb11, new( ibord(ib,ia) ) - n )
      lb22 = max( lb22, new( ibord(ib,ia) ) )
    end do
    do l = 1,nlmsa(ia)
      ii = ianew(ia) + l
      lb1(ii) = lb11
      lb2(ii) = lb22
    end do
  end do

  lb11 = nligne - nlmso + 1
  do ib = 1,nsortf
    lb11 = min(lb11,new(isrt(ib))-n)
  end do
  do ii = nligne-nlmso+1,nligne
    lb1(ii) = lb11
    lb2(ii) = nligne
  end do

  do ii = nligne-1,1,-1
    lb1(ii) = min(lb1(ii),lb1(ii+1))
  end do

  newinv = 0
  do i = 1,npoint
    ii = new(i)
    if( ii /= 0 ) then
      newinv(ii) = i
      if( nspino == 2 ) newinv(ii-1) = i
    endif
  end do
  do ia = 1,natome
    do lm = 1,nlmsa(ia)
      newinv(ianew(ia)+lm) = - ia
    end do
  end do

  if( icheck > 1 ) then
    write(3,110)
    write(3,120) nligne
    write(3,130)
    write(3,140) ianew(1:natome)
    write(3,145)
    write(3,140) nlmsa(1:natome)
    write(3,146)
    write(3,140) nlmso
  endif
  if( icheck > 2 ) then
    write(3,150)
    write(3,155) (i, new(i), i = 1,npoint)
    write(3,160)
    write(3,170) (ii, lb1(ii), lb2(ii), newinv(ii), ii = 1,nligne)
  endif

  return
  105 format(//' Error in newind for i, iv =',2i6)
  110 format(/' ---- Newind -------',100('-'))
  120 format(/' Number of line =',i6)
  130 format(' ianew(ia) =')
  140 format(12i6)
  145 format(' nlmsa(ia) =')
  146 format(' nlmso =')
  150 format(10('    i  new(i)'))
  155 format(10(i7,i6))
  160 format(/2('    ii   lb1   lb2   newinv'))
  170 format(2(3i6,i9))
end

!***********************************************************************

! Sous programme d'inversion de matrices generales complexes

subroutine invcomp(n,mat,nm,lwork,is,Stop_job)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer:: ipiv(nm)

  complex(kind=db):: mat(nm,nm), work(lwork)

  logical:: Stop_job

  if( is == 0 ) then
!cray      call cgetrf(n,n,mat,nm,ipiv,info)
    call zgetrf(n,n,mat,nm,ipiv,info)
  elseif( is == 1 ) then
!cray      call csytrf('u',n,mat,nm,ipiv,work,lwork,info)
    call zsytrf('u',n,mat,nm,ipiv,work,lwork,info)
  else
    call zhetrf('u',n,mat,nm,ipiv,work,lwork,info)
  endif

  if( info /= 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,110) info
    end do
    if( Stop_job ) stop
    Stop_job = .true.
  endif

  if( is == 0 ) then
!cray        call cgetri(n,mat,nm,ipiv,work,lwork,info)
    call zgetri(n,mat,nm,ipiv,work,lwork,info)
  elseif( is == 1 ) then
!cray      call csytri('u',n,mat,nm,ipiv,work,info)
    call zsytri('u',n,mat,nm,ipiv,work,info)
    do i = 1,nm
      do j = i + 1,nm
        mat(j,i) = mat(i,j)
     end do
    end do
  else
    call zhetri('u',n,mat,nm,ipiv,work,info)
    do i = 1,nm
      do j = i + 1,nm
        mat(j,i) = conjg( mat(i,j) )
     end do
    end do
  endif

  if( info /= 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,120) info
    end do
    if( Stop_job ) stop
    Stop_job = .true.
  endif

  return
  110 format(/' CGETRF : info =',i5)
  120 format(/' CGETRI : info =',i5)
end

!***********************************************************************

! Recopie des amplitudes pour les (l,m) non calculees

subroutine recop_taull(nlmagm,nspinp,Sym_cubic,taull,Ylm_comp)

  use declarations
  implicit none

  integer:: nlmagm, nspinp

  logical:: Sym_cubic, Ylm_comp

  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp):: taull

  if( Sym_cubic ) then
! l == 1
    if( nlmagm > 1 ) then
      taull(2,:,2,:) = taull(3,:,3,:)
      taull(4,:,4,:) = taull(3,:,3,:)
    endif

! l == 2
    if( nlmagm > 4 ) then
      if( Ylm_comp ) then
        if( abs( taull(8,1,8,1) ) > eps10 ) then
          taull(5,:,5,:) = 0.5_db*( taull(7,:,7,:) + taull(8,:,8,:))
          taull(9,:,9,:) = taull(5,:,5,:)
          taull(5,:,9,:) = 0.5_db*( taull(7,:,7,:) - taull(8,:,8,:))
          taull(9,:,5,:) = taull(5,:,9,:)
        else
          taull(8,:,8,:) = 2 * taull(5,:,5,:) - taull(7,:,7,:)
          taull(6,:,6,:) = taull(8,:,8,:)
        endif
      else
        taull(9,:,9,:) = taull(7,:,7,:)
        taull(6,:,6,:) = taull(5,:,5,:)
        taull(8,:,8,:) = taull(5,:,5,:)
      endif
    endif

  elseif( .not. Ylm_comp ) then
! sym_4
    if( nlmagm > 1 ) taull(4,:,4,:) = taull(2,:,2,:)
    if( nlmagm > 4 ) taull(6,:,6,:) = taull(8,:,8,:)

  endif

  return
end

!***********************************************************************

subroutine soustract_tl(Axe_Atom_grn,Cal_xanes,Full_atom,iaabsi,iaprotoi,igreq,igroupi,ispinin,lmaxa,n_atom_0, &
                n_atom_ind,n_atom_proto,natome,neqm,ngroup_m,nlmagm,nspin,nspino,nspinp,Spinorbite,State_all_r,Tau_ato,Taull)

  use declarations
  implicit none

  integer:: i1, i2, ia, iaabsi, iapr, iga, igr, ipr, is1, is2, isp1, isp2, ispinin, its, kgr, l1, l2, &
     lm10, lm20, lm1, lm2, m1, m2, n_atom_0, n_atom_ind, n_atom_proto, natome, neqm, ngroup_m, nlmagm, nspin, nspino, nspinp

  integer, dimension(natome):: iaprotoi, igroupi, lmaxa
  integer, dimension(0:n_atom_proto,neqm):: igreq

  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,n_atom_0:n_atom_ind):: tau_ato
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,natome):: taull

  logical:: Cal_xanes, Full_atom, Spinorbite, State_all_r

  real(kind=db):: cosang
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_grn

  do ia = 1,natome

    if( ia /= iaabsi .and. .not. State_all_r .and. Cal_xanes ) cycle

    if( nspin == 2 ) then
      igr = igroupi(ia)
      ipr = iaprotoi(ia)
      kgr = igreq(ipr,1)
      iga = igroupi(iaabsi)

      cosang = sum( Axe_Atom_grn(:,iga) * Axe_Atom_grn(:,igr) )

      if( abs( cosang - 1 ) < eps4 ) then
        its = 1
      elseif( abs( cosang + 1 ) < eps4 ) then
        its = -1
      else
        its = 0
      endif
    else
      its = 1
    endif

    if( Full_atom ) then
      iapr = ia
    else
      iapr = iaprotoi(ia)
    endif

    lm1 = 0
    do l1 = 0,lmaxa(ia)
      do m1 = -l1,l1
        lm1 = lm1 + 1
        do is1 = 1,nspino
          if( Spinorbite ) then
            isp1 = is1
          else
            isp1 = ispinin
          endif
          if( its == - 1 ) then
            i1 = 3 - isp1
            lm10 = l1**2 + l1 + 1 - m1
          else
            i1 = isp1
            lm10 = l1**2 + l1 + 1 + m1
          endif
          lm2 = 0
          do l2 = 0,lmaxa(ia)
! Meme en Hubbard, les l differents ne se melangent pas
            do m2 = -l2,l2
              lm2 = lm2 + 1
              if( l2 /= l1 ) cycle
              do is2 = 1,nspino
!                    if( is1 /= is2 ) cycle
                if( Spinorbite ) then
                  isp2 = is2
                else
                  isp2 = ispinin
                endif
                if( its == - 1 ) then
                  i2 = 3 - isp2
                  lm20 = l2**2 + l2 + 1 - m2
                else
                  i2 = isp2
                  lm20 = l2**2 + l2 + 1 + m2
                endif
                if( its == 0 ) then
                  if( i1 == i2 ) taull(lm1,isp1,lm2,isp2,ia) = taull(lm1,isp1,lm2,isp2,ia) - 0.5_db &
                     * ( tau_ato(lm10,1,lm20,1,iapr) + tau_ato(lm10,nspinp,lm20,nspinp,iapr) )
                else
                  taull(lm1,isp1,lm2,isp2,ia) = taull(lm1,isp1,lm2,isp2,ia) - tau_ato(lm10,i1,lm20,i2,iapr)
                end if
              end do
            end do
          end do
        end do
      end do
    end do

  end do

  return
end

!***********************************************************************

! Calcul des amplitudes de diffusions selon la theorie de la diffusion multiple.

subroutine msm(Axe_atom_grn,Base_ortho,Cal_xanes,dcosxyz,ecinetic,Eimag,Full_atom,ia_eq,ia_rep,iaabsi,iaprotoi, &
                    iato,icheck,igreq,igroupi,igrph,iopsymr,irep_util,is_eq,ispin,karact,lato,lmaxa,lmaxg,mato,n_atom_0, &
                    n_atom_ind,n_atom_proto,natome,natomp,nb_eq,nb_rpr,nb_rep_t,nb_sym_op,nchemin,neqm,ngroup_m, &
                    ngrph,nlmagm,nlmsa,nlmsam,nlmsamax,Normaltau,nspin,nspino,nspinp,pos,posi,recop,Repres_comp, &
                    rot_atom,Solsing,Spinorbite,State_all_r,Sym_cubic,Tau_ato,Tau_nondiag,Taull,Time_fill,Time_tria,Ylm_comp)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  parameter(nletm=52)

  character(len=1) let(0:nletm)
  character(len=3), dimension(:), allocatable:: motval
  character(len=3), dimension(:,:,:,:), allocatable:: mot
  
  complex(kind=db):: fnormc, gmat, gmatc, gmats, kr
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nopsm,nrepm):: karact
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,n_atom_0:n_atom_ind)::tau_ato
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,natome):: taull
  complex(kind=db), dimension(:), allocatable:: bessel, neuman, val, ylmc
  complex(kind=db), dimension(:,:), allocatable:: hankelm, mat, matp, taullp, taullq
  complex(kind=db), dimension(:,:,:), allocatable:: Cmatr, Dlmm_ia, taup, taup_inv
  complex(kind=db), dimension(:,:,:,:), allocatable :: taull_tem
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: Cmat

  integer, dimension(natome):: iaprotoi, igroupi, lmaxa, nb_eq, nlmsa
  integer, dimension(nopsm):: iopsymr
  integer, dimension(natome,natome):: nb_rpr
  integer, dimension(nb_sym_op,natome,natome):: ia_rep, nb_rep_t
  integer, dimension(nb_sym_op,natome):: ia_eq, is_eq
  integer, dimension(nrepm,2):: irep_util
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(nlmsam,natome,ngrph):: lato, mato, iato

  logical:: Base_ortho, Brouder, Cal_xanes, Ereel, Full_atom, Normaltau, Solsing, &
    Recop, Repres_comp, Spinorbite, State_all_r, Stop_job, Sym_cubic, Tau_nondiag, Ylm_comp

  real(kind=db):: tp1, tp2, tp3, Time_fill, Time_tria

  real(kind=db), dimension(nspin):: ecinetic
  real(kind=db), dimension(nspin):: konder
  real(kind=db), dimension(3):: dcosxyz, w
  real(kind=db), dimension(3,3):: Mat_rot
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_grn
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(3,natomp):: pos
  real(kind=db), dimension(3,3,natome):: rot_atom
  real(kind=db), dimension(:), allocatable :: besselr, neumanr, Ylm

  real(kind=sg) time

  data let/' ','a','b','c','d','e','f','g','h','i','j','k','l','m', &
   'n','o','p','q','r','s','t','u','v','w','x','y','z','A','B','C', &
   'D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

  Stop_job = .true.

! Pour l'instant toutes les representations sont de dimension 1.
  I_rep_dim = 1

! Correspond a Ch. Brouder et al. PRB 54, 7334 (1996) pour la matrice
! de diffusion multiple plus stable. Non programme avec spin-orbite.
  Brouder = .false.

! Remplissage de la matrice
  call CPU_TIME(time)
  tp1 = real(time,db)

  ndim = sum( nlmsa(1:natome) )
  nlmabs = nlmsa(iaabsi)
  nlmsmax = 0
  do ia = 1,natome
    nlmsmax = max( nlmsmax, nlmsa(ia) )
  end do

  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag, db) )

  if( icheck > 1 ) then
    write(3,110)
    write(3,115) igrph, ndim
  endif

  if( ndim /= 0 ) then

  allocate( taup(nlmsmax,nlmsmax,natome) )
  allocate( taullp(nlmsmax,nlmsmax) )
  allocate( taullq(nlmsmax,nlmsmax) )

  if( nchemin == 2 ) then
    nlmch = nlmabs
  else
    nlmch = ndim
  endif
  Taup(:,:,:) = (0._db, 0._db)

  if( Tau_nondiag ) allocate( Taup_inv(nlmsmax,nlmsmax,natome) )

  boucle_at: do ia = 1,natome

    if( nspin == 2 ) then
      igr = igroupi(ia)
      iga = igroupi(iaabsi)
      ipr = iaprotoi(ia)
      kgr = igreq(ipr,1)

      cosang = sum( Axe_Atom_grn(:,iga) * Axe_Atom_grn(:,igr) )
 !     cosang = sum( Axe_Atom_grn(:,kgr) * Axe_Atom_grn(:,igr) )
 !     cosang = cosang * abs( sum( Axe_Atom_grn(:,iga) * Axe_Atom_grn(:,igr) ) )
      if( abs( cosang - 1 ) < eps4 ) then
        its = 1
      elseif( abs( cosang + 1 ) < eps4 ) then
        its = -1
      else
        its = 0
      endif
    else
      its = 1
    endif

    if( Full_atom ) then
      iapr = ia
    else
      iapr = iaprotoi(ia)
    endif

    do lm1 = 1,nlmsa(ia)

      l1 = lato(lm1,ia,igrph)
      m1 = mato(lm1,ia,igrph)
      if( nspino == 2 ) then
        i1 = iato(lm1,ia,igrph)
      else
        i1 = ispin
      endif
      if( its == - 1 ) then
        i1 = 3 - i1
        m1 = - m1
      endif
      lm10 = l1**2 + l1 + 1 + m1
      do lm2 = 1,nlmsa(ia)
        l2 = lato(lm2,ia,igrph)
        m2 = mato(lm2,ia,igrph)
        if( nspino == 2 ) then
          i2 = iato(lm2,ia,igrph)
        else
          i2 = ispin
        endif
        if( its == - 1 ) then
          i2 = 3 - i2
          m2 = - m2
        endif
        lm20 = l2**2 + l2 + 1 + m2
        if( its == 0 ) then
          if( i1 == i2 ) taup(lm1,lm2,ia) = 0.5_db * ( tau_ato(lm10,1,lm20,1,iapr) + tau_ato(lm10,nspinp,lm20,nspinp,iapr) )
        else
          Taup(lm1,lm2,ia) = Tau_ato(lm10,i1,lm20,i2,iapr)
        endif
      end do
    end do
    if( Tau_nondiag .and. nlmsa(ia) > 0 ) then
      nlm = nlmsa(ia)
      allocate( Mat(nlm,nlm) )
      Mat(1:nlm,1:nlm) = Taup(1:nlm,1:nlm,ia)
      Mat = Transpose( Mat )
      if( Ylm_comp ) then
        call invcomp(nlm,mat,nlm,nlm,0,Stop_job)
      else
        call invcomp(nlm,mat,nlm,nlm,1,Stop_job)
      endif
      Taup_inv(1:nlm,1:nlm,ia) = Mat(1:nlm,1:nlm)
      deallocate( Mat )
    endif
  end do boucle_at

  allocate( mat(ndim,nlmch) )
  mat(:,:) = (0._db, 0._db)

  konder(1:nspin) = real( konde(1:nspin), db )
  if( abs(Eimag) < eps10 .and. Ecinetic(1) > eps10 .and. Ecinetic(nspin) > eps10 ) then
    Ereel = .true.
  else
    Ereel = .false.
  endif

  if( natome > 1 .or. nb_sym_op > 1 ) then
    allocate( Cmat(natome,nlmsamax,nb_sym_op,-lmaxg:lmaxg,nspino) )
    call Cmat_cal(Cmat,iato,icheck,igrph,iopsymr,irep_util,is_eq,karact,lato,lmaxg,mato,natome,nb_eq,nb_sym_op,ngrph, &
                nlmsamax,nlmsa,nlmsam,nspino,rot_atom,Ylm_comp)
  endif

  fnorm = 1._db
  fnormc = (1._db,0._db)

  iligne = 0
  icolone = 0
  if( nchemin == 2 ) then
    na = 1
  else
    na = natome
  endif

  do ia = 1,na   ! boucle sur la ligne
    iligne0 = sum( nlmsa(1:ia-1) )

    if( nb_sym_op > 1 .or. Normaltau .or. Ylm_comp ) then
      ib1 = 1
    else
      ib1 = ia + 1
    endif

    do ib = ib1,natome   ! boucle sur la colonne
      icolone0 = sum( nlmsa(1:ib-1) )
      lmax = lmaxa(ia) + lmaxa(ib)
      if( Ereel ) then
        allocate( neumanr(0:lmax) )
        allocate( besselr(0:lmax) )
      else
        allocate( neuman(0:lmax) )
        allocate( bessel(0:lmax) )
      endif
      nspinr = min( nspin, nspino )
      allocate( hankelm(0:lmax,nspinr) )
      nlmr = (lmax + 1)**2
      nlmc = ( (lmax + 1) * (lmax + 2) ) / 2
      allocate( ylmc(nlmc) )
      if( .not. Ylm_comp ) allocate( ylm(nlmr) )

! Boucle sur les atomes n1 tels qu'il existe une symetrie S avec :
! S(n0) = n0 et S(n1) = n', ou n' appartient a p'.
      do n1 = 1,nb_rpr(ia,ib)

! ibb : indice cluster non symetrise
        ibb = ia_rep(n1,ia,ib)

! L'atome ne diffuse pas vers lui-meme.
        w(1:3) = posi(1:3,ia) - pos(1:3,ibb)
        r = vnorme(Base_ortho,dcosxyz,w)
        if( r < eps10 ) cycle

! Recherche de l'indice "ind" correspondant au representant n1.
        do ind = 1,nb_eq(ib)
          if( ibb == ia_eq(ind,ib) ) exit
        end do

        fac = real( nb_eq(ia) * nb_rep_t(n1,ia,ib), db ) / I_rep_dim

! Calcul des fonctions de hankel
        do isp = 1,nspinr
          if( Spinorbite ) then
            jsp = isp
          else
            jsp = ispin
          endif
          if( Ereel ) then
            rkr = konder(jsp) * r
            call cbessneur(fnorm,rkr,lmax,lmax,besselr,neumanr)
            hankelm(:,isp) = cmplx(besselr(:),neumanr(:),db)
          else
            kr = konde(jsp) * r
            call cbessneu(fnormc,kr,lmax,lmax,bessel,neuman)
            hankelm(0:lmax,isp) = bessel(0:lmax) + img*neuman(0:lmax)
          endif
          do l = 1,lmax
            hankelm(l,isp) = img**l * hankelm(l,isp)
          end do
        end do

! Calcul des Ylm
        call cylm(lmax,w,r,Ylmc,nlmc)
        if( .not. Ylm_comp ) call ylmcr(lmax,nlmc,nlmr,Ylmc,Ylm)

        do lm1 = 1,nlmsa(ia)

          iligne = iligne0 + lm1
          l1 = lato(lm1,ia,igrph)
          m1 = mato(lm1,ia,igrph)
          is1 = iato(lm1,ia,igrph)
          if( Spinorbite ) then
            isp1 = min( is1, nspin )
          else
            isp1 = 1
          endif

          do lm2 = 1,nlmsa(ib)

            icolone = icolone0 + lm2
            l2 = lato(lm2,ib,igrph)
            m2 = mato(lm2,ib,igrph)
            is2 = iato(lm2,ib,igrph)
            if( is2 /= is1 ) cycle

            lma = abs( l1 - l2 )
            lmb = l1 + l2

! Boucles sur m et m' associes a la symetrie
            gmats = (0._db, 0._db)
            do ms1 = -l1,l1
              if( nb_sym_op == 1 ) then
                if( ms1 /= m1 ) cycle
              else
                if( abs( Cmat(ia,lm1,1,ms1,is1) ) < eps10 ) cycle
              endif
              do ms2 = -l2,l2

                if( nb_sym_op == 1 ) then
                  if( ms2 /= m2 ) cycle
                else
                  if( abs( Cmat(ib,lm2,ind,ms2,is2) ) < eps10) cycle
                endif

                gmat = (0._db,0._db)
                do l = lma,lmb,2
                  ll = l**2 + l + 1
                  lm0 = ( l**2 + l ) / 2 + 1
                  if( Ylm_comp ) then
                    gmatc = (0._db,0._db)
                  else
                    gmatr = 0._db
                  endif
                  do m = -l,l
                    if( Ylm_comp ) then
! c'est Y(l2,m2) qui est complexe conjugue
                      g = gauntcp(l2,ms2,l1,ms1,l,m)
                      if( abs( g ) < eps10 ) cycle
                      lmc = lm0 + abs(m)
                      if( m >= 0 ) then
                        gmatc = gmatc + ylmc(lmc) * g
                      else
                        gmatc = gmatc + (-1)**m * conjg(ylmc(lmc)) * g
                      endif
                    else
                      g = gauntc(l1,ms1,l2,ms2,l,m)
                      if( abs( g ) < eps10 ) cycle
                      lm = ll + m
                      gmatr = gmatr + ylm(lm) * g
                    endif
                  end do
                  if( Ylm_comp ) then
                    gmat = gmat + gmatc * hankelm(l,isp1)
                  else
                    gmat = gmat + gmatr * hankelm(l,isp1)
                  endif

                end do

                gmat = - quatre_pi * gmat * img**(l1-l2+1)

                if( nb_sym_op > 1 ) then
                  gmats = gmats + Cmat(ia,lm1,1,ms1,is1) * gmat * Conjg( Cmat(ib,lm2,ind,ms2,is2) )
                else
                  gmats = gmats + gmat
                endif

              end do
            end do  ! Fin des boucles sur m, m' liees a la symetrie

            mat(iligne,icolone) = mat(iligne,icolone) - fac * gmats
            if( .not. ( nb_sym_op > 1 .or. normaltau .or. Ylm_comp ) ) mat(icolone,iligne) = mat(iligne,icolone)

          end do
        end do  ! fin boucle sur l,m

      end do ! fin de la boucle sur les atomes equivalents

      if( Ereel ) then
        deallocate( besselr )
        deallocate( neumanr )
      else
        deallocate( bessel )
        deallocate( neuman )
      endif

      deallocate( hankelm )
      if( .not. Ylm_comp ) deallocate( Ylm )
      deallocate( ylmc )
    end do
  end do

! diagonale
  if( normaltau .or. Brouder ) then
    if( nchemin == - 1 ) then
      do iligne = 1,ndim
        mat(iligne,iligne) = mat(iligne,iligne) + (1._db,0._db)
      end do
    endif
  else
    iligne = 0
    do ia = 1,natome
      jligne0 = iligne
      do lm = 1,nlmsa(ia)
        iligne = iligne + 1
        if( tau_nondiag ) then
          do lm2 = 1,nlmsa(ia)
            jligne = jligne0 + lm2
            mat(iligne,jligne) = mat(iligne,jligne) + Taup_inv(lm,lm2,ia)
          end do
        else
          mat(iligne,iligne) = mat(iligne,iligne) + 1/Taup(lm,lm,ia)
        endif
      end do
    end do
  endif

  if( icheck > 2 ) then
    write(3,120)
    iligne = 0
    do ia = 1,natome
      do lm = 1,nlmsa(ia)
        iligne = iligne + 1
        if( Spinorbite ) then
          write(3,130) iligne, ia, lato(lm,ia,igrph), mato(lm,ia,igrph), iato(lm,ia,igrph)
        else
          write(3,130) iligne, ia, lato(lm,ia,igrph), mato(lm,ia,igrph), ispin
        endif
        write(3,140) mat(iligne,1:ndim)
      end do
    end do
  endif

  if( ( normaltau .or. Brouder ) .and. nchemin /= 2 ) then
    i0 = 0
    do ia = 1,natome
      j0 = 0
      do ib = 1,natome
        if( ia == ib ) then
          j0 = j0 + nlmsa(ib)
          cycle
        endif

        allocate( matp(nlmsa(ia),nlmsa(ib)) )
        do lma = 1,nlmsa(ia)
          i = i0 + lma
          do lmb = 1,nlmsa(ib)
            matp(lma,lmb) = mat(i,j0+lmb)
          end do
        end do

        if( Brouder ) then

          if( Spinorbite ) then
            write(6,'(//3x,A)') 'Brouder non programme avec Spinorbite !'
            stop
          else
            do lma = 1,nlmsa(ia)
              i = i0 + lma
              do lmb = 1,nlmsa(ib)
                mat(i,j0+lmb) = sqrt( Taup(lma,lma,ia) ) * matp(lma,lmb) * sqrt( Taup(lmb,lmb,ib) )
              end do
            end do
          endif

        else

          do lma = 1,nlmsa(ia)
            i = i0 + lma
            do lmb = 1,nlmsa(ib)
              mat(i,j0+lmb) = sum( taup(lma,1:nlmsa(ia),ia) * matp(1:nlmsa(ia),lmb) )
            end do
          end do

        endif

        deallocate( matp )
        j0 = j0 + nlmsa(ib)
      end do
      i0 = i0 + nlmsa(ia)
    end do
  endif

  call CPU_TIME(time)
  tp2 = real(time,db)

  Time_fill = tp2 - tp1

  if( natome == 1 ) then

    ia = 1
    do lm1 = 1,nlmsa(ia)
      l = lato(lm1,ia,igrph)
      m1 = mato(lm1,ia,igrph)
      lm01 = l**2 + l + 1 + m1
      if( Spinorbite ) then
        is1 = iato(lm1,ia,igrph)
      else
        is1 = ispin
      endif
      do lm2 = 1,nlmsa(ia)
        l = lato(lm2,ia,igrph)
        m2 = mato(lm2,ia,igrph)
        lm02 = l**2 + l + 1 + m2
        if( Spinorbite ) then
          is2 = iato(lm2,ia,igrph)
        else
          is2 = ispin
        endif
        Taull(lm01,is1,lm02,is2,ia) = Taull(lm01,is1,lm02,is2,ia) + Taup(lm1,lm2,ia)
        if( .not. Spinorbite .and. Repres_comp ) then
          lm01c = lm01 - 2 * m1
          lm02c = lm02 - 2 * m2
          isg = (-1)**(m1+m2)
          Taull(lm02c,is2,lm01c,is1,ia) = Taull(lm02c,is2,lm01c,is1,ia) + isg * Taup(lm1,lm2,ia)
        endif
      end do
    end do

! Le remplissage de la matrice est termine, suit la resolution.

  elseif( nchemin > -1 ) then

! Developpement en chemin
    call dev_chemin(iaabsi,iato,igrph,lato,mat,mato,natome,nchemin,ndim,ngrph,nlmabs,nlmch,nlmagm, &
             nlmsa,nlmsam,nlmsmax,nspinp,taull,taup)

  else

! Full multiple scattering, inversion de la matrice

    if( normaltau .or. nb_sym_op > 1 .or. Ylm_comp ) then
      call invcomp(ndim,mat,ndim,ndim,0,Stop_job)
    else
      call invcomp(ndim,mat,ndim,ndim,1,Stop_job)
    endif

    if( Brouder ) then

      do lm = 1,ndim
        mat(lm,lm) = mat(lm,lm) - (1._db,0._db)
      end do

      i0 = 0
      do ia = 1,natome
        j0 = 0
        do ib = 1,natome

          if( Spinorbite ) then
          else
            do lma = 1,nlmsa(ia)
              i = i0 + lma
              do lmb = 1,nlmsa(ib)
                mat(i,j0+lmb) = sqrt( taup(lma,lma,ia) ) * mat(i,j0+lmb) * sqrt( taup(lmb,lmb,ib) )
              end do
              if( ia == ib ) mat(i,i) = mat(i,i) + Taup(lma,lma,ia)
            end do
          endif

          j0 = j0 + nlmsa(ib)
        end do
        i0 = i0 + nlmsa(ia)
      end do

    endif

    lm0 = 0

    do ia = 1,natome

      if( ia > 1 ) lm0 = lm0 + nlmsa(ia-1)

      if( ia /= iaabsi .and. .not. State_all_r .and. Cal_xanes ) cycle

      if( normaltau ) then

        do lm1 = 1,nlmsa(ia)
          lmd1 = lm0 + lm1
          do lm2 = 1,nlmsa(ia)
            lmd2 = lm0 + lm2
            Taullp(lm2,lm1) =  sum( mat(lmd2,lm0+1:lm0+nlmsa(ia)) * Taup(1:nlmsa(ia),lm1,ia) )
!               taullp(lm1,lm2) = taullp(lm2,lm1)
          end do
        end do

      else

        do lm1 = 1,nlmsa(ia)
          lmd1 = lm0 + lm1
          taullp(lm1,lm1) = mat(lmd1,lmd1)
! A cause du cas harmonique complexe, la matrice n'est plus symetrique
          do lm2 = 1,nlmsa(ia)
            lmd2 = lm0 + lm2
            taullp(lm1,lm2) = mat(lmd1,lmd2)
          end do
        end do

      endif

      if( nb_sym_op > 1 ) then

        allocate( Cmatr(nlmsamax,-lmaxg:lmaxg,nspino) )

        Mat_rot(:,:) = rot_atom(:,:,ia)
        Mat_rot = Transpose( Mat_rot )
        allocate( Dlmm_ia(-lmaxa(ia):lmaxa(ia),-lmaxa(ia):lmaxa(ia), 0:lmaxa(ia)) )
        call Rotation_mat(Dlmm_ia,icheck,Mat_rot,lmaxa(ia), Ylm_comp)
        do lm = 1,nlmsa(ia)
          l = lato(lm,ia,igrph)
          do isp = 1,nspino
            do m = -l,l
              Cmatr(lm,m,isp) = sum( Cmat(ia,lm,1,-l:l,isp) * Dlmm_ia(-l:l,m,l) )
            end do
          end do
        end do
        deallocate( Dlmm_ia )

        taullq(:,:) = taullp(:,:)
        taullp(:,:) = (0._db,0._db)

        do lm1 = 1,nlmsa(ia)
          l1 = lato(lm1,ia,igrph)
          m1 = mato(lm1,ia,igrph)
          is1 = iato(lm1,ia,igrph)
          do lm2 = 1,nlmsa(ia)
            l2 = lato(lm2,ia,igrph)
            m2 = mato(lm2,ia,igrph)
            is2 = iato(lm2,ia,igrph)

            do lm1p = 1,nlmsa(ia)
              if( lato(lm1p,ia,igrph) /= l1 ) cycle
              if( abs( Cmatr(lm1p,m1,is1) ) < eps10 ) cycle
              do lm2p = 1,nlmsa(ia)
                if( lato(lm2p,ia,igrph) /= l2 ) cycle
                if( abs( Cmatr(lm2p,m2,is2) ) < eps10 ) cycle

                taullp(lm1,lm2) = taullp(lm1,lm2) + Cmatr(lm1p,m1,is1) * taullq(lm1p,lm2p) * Conjg( Cmatr(lm2p,m2,is2) )

              end do
            end do

          end do
        end do

        deallocate( Cmatr )

      endif

      do lm1 = 1,nlmsa(ia)
        l = lato(lm1,ia,igrph)
        m1 = mato(lm1,ia,igrph)
        lm01 = l**2 + l + 1 + m1
        if( Spinorbite ) then
          is1 = iato(lm1,ia,igrph)
        else
          is1 = ispin
        endif
        do lm2 = 1,nlmsa(ia)
          l = lato(lm2,ia,igrph)
          m2 = mato(lm2,ia,igrph)
          lm02 = l**2 + l + 1 + m2
          if( Spinorbite ) then
            is2 = iato(lm2,ia,igrph)
          else
            is2 = ispin
          endif
          taull(lm01,is1,lm02,is2,ia) = taull(lm01,is1,lm02,is2,ia) + taullp(lm1,lm2)
          if( .not. Spinorbite .and. Repres_comp ) then
            lm01c = lm01 - 2 * m1
            lm02c = lm02 - 2 * m2
            isg = (-1)**(m1+m2)
            taull(lm02c,is2,lm01c,is1,ia) = taull(lm02c,is2,lm01c,is1,ia) + isg * taullp(lm1,lm2)
          endif
        end do
      end do

    end do ! fin boucle sur ia

  endif

  if( natome > 1 .or. nb_sym_op > 1 ) deallocate( Cmat )

  if( icheck > 2 ) then
    write(3,150)
    iligne = 0
    do ia = 1,natome
      do lm = 1,nlmsa(ia)
        iligne = iligne + 1
        if( Spinorbite ) then
          write(3,130) iligne, ia, lato(lm,ia,igrph), mato(lm,ia,igrph), iato(lm,ia,igrph)
        else
          write(3,130) iligne, ia, lato(lm,ia,igrph), mato(lm,ia,igrph), ispin
        endif
        write(3,140) mat(iligne,1:ndim)
      end do
    end do
  endif

  deallocate( mat )
  deallocate( taullp )
  deallocate( taullq )
  deallocate( taup )
  if( tau_nondiag ) deallocate( taup_inv )

  endif            ! arrivee en cas de ndim = 0

  if( recop .and. igrph == ngrph .and. ( Spinorbite .or. ispin == nspin ) ) then
    allocate( Taull_tem(nlmagm,nspinp,nlmagm,nspinp) )
    Taull_tem(:,:,:,:) = Taull(:,:,:,:,iaabsi)
    call recop_taull(nlmagm,nspinp,Sym_cubic,Taull_tem,Ylm_comp)
    Taull(:,:,:,:,iaabsi) = Taull_tem(:,:,:,:)
    deallocate( taull_tem )
  endif

! Soustraction de l'amplitude de diffusion atomique
  if( Solsing .and. igrph == ngrph ) then
    if( icheck > 2 ) then
      nlmw = min(nlmagm,9)
      write(3,161)
      lm1 = 0
      do l1 = 0,5
        do m1 = -l1,l1
          lm1 = lm1 + 1
          if( lm1 > nlmw ) exit
          do is1 = 1,nspinp
            write(3,218) l1, m1, is1, (( Taull(lm1,is1,lm2,is2,iaabsi), is2 = 1,nspinp ), lm2 = 1,nlmw )
          end do
        end do
      end do
    end if
    call soustract_tl(Axe_Atom_grn,Cal_xanes,Full_atom,iaabsi,iaprotoi,igreq,igroupi,ispin,lmaxa,n_atom_0, &
                n_atom_ind,n_atom_proto,natome,neqm,ngroup_m,nlmagm,nspin,nspino,nspinp,Spinorbite,State_all_r,Tau_ato,Taull)
  end if

  if( icheck > 1 .and. igrph == ngrph .and. ( Spinorbite .or. ispin == nspin ) ) then

    do ia = 1,natome
      if(  ia /= iaabsi .and. .not. State_all_r .and. Cal_xanes ) cycle
      if( ia == iaabsi ) then
        write(3,160)
      else
        write(3,162) ia
      endif
      lmw = min( 3, lmaxa(ia) )
      nlmw = min( nlmagm, (lmw+1)**2 )

      allocate( val(nletm**2) )
      allocate( motval(nletm**2) )
      allocate( mot(nlmagm,nspinp,nlmagm,nspinp) )
      val(:) = (0._db, 0._db)
      motval(:) = ' '
      mot(:,:,:,:) = ' '
      nlet = 0
      letu = 0 
      letd = 0
      do is1 = 1,nspinp
        lm1 = 0
        do l1 = 0,lmw
          do m1 = -l1,l1
            lm1 = lm1 + 1
            do is2 = 1,nspinp
              lm2 = 0
              do l2 = 0,lmw
                boucle_m2: do m2 = -l2,l2
                  lm2 = lm2 + 1
                  if( abs( Taull(lm1,is1,lm2,is2,ia) ) < eps15 ) then
                    mot(lm1,is1,lm2,is2) = '  .'
                  else
                    do ival = 1,nlet
                      if( abs( Taull(lm1,is1,lm2,is2,ia) - val(ival) ) < eps10 ) then
                        mot(lm1,is1,lm2,is2) = motval(ival)
                        cycle boucle_m2
                      elseif( abs( Taull(lm1,is1,lm2,is2,ia) + val(ival) ) < eps10 ) then
                        mot(lm1,is1,lm2,is2) = motval(ival)
                        l = len_trim( adjustl( motval(ival) ) )
                        mot(lm1,is1,lm2,is2)(3-l:3-l) = '-'
                        cycle boucle_m2
                      endif
                    end do

                    if( nlet == nletm**2 ) then
                      mot(lm1,is1,lm2,is2) = '  ?'
                    else
                      nlet = nlet + 1
                      letu = letu + 1
                      if( letu > nletm ) then
                        letu = 1
                        letd = letd + 1
                      endif
                      val(nlet) = Taull(lm1,is1,lm2,is2,ia)
                      motval(nlet)(2:3) = let(letd) // let(letu) 
                      mot(lm1,is1,lm2,is2) =  motval(nlet) 
                    endif
                  endif  
                end do boucle_m2
              end do
            end do
          end do
        end do
      end do

      write(3,220) 'l', ((( l, m = -l,l), l = 0,lmw ), i = 1,nspinp )
      write(3,220) 'm', ((( m, m = -l,l), l = 0,lmw ), i = 1,nspinp )
      if( nspinp == 2 ) then
        write(3,220) 's', ((( i, m = -l,l), l = 0,lmw ), i = 1,nspinp )
        write(3,'(A)') '  l  m  s'
      else
        write(3,'(A)') '     l  m'
      endif
      do is1 = 1,nspinp
        lm1 = 0
        do l1 = 0,lmw
          do m1 = -l1,l1
            lm1 = lm1 + 1
            if( nspinp == 1 ) then
              write(3,230) l1, m1, (( mot(lm1,is1,lm2,is2), lm2 = 1,nlmw), is2 = 1,nspinp )
            else
              write(3,240) l1, m1, is1, (( mot(lm1,is1,lm2,is2), lm2 = 1,nlmw), is2 = 1,nspinp )
           endif
          end do
        end do
      end do
      
      write(3,*)
      write(3,250) ( motval(i), val(i), i = 1,nlet )
      deallocate( mot, motval, val)
    end do
  endif

! Ecriture pour extract
  if( icheck > 0 .and. igrph == ngrph .and. ( Spinorbite .or. ispin == nspin ) .and. Cal_xanes ) then
    if( icheck == 1 ) write(3,110)
    ia = iaabsi
    ! si lmaxa > 10 il faut changer le format d'ecriture
    lmw = min( 10, lmaxa(ia) )
    nlmw = min( nlmagm, (lmw+1)**2 )
    write(3,260)
    if( nspinp == 1 ) then
      write(3,270) lmw, nspinp
      write(3,280) (( l, m, m = -l,l ), l = 0,lmw )
    else
      write(3,290) lmw, nspinp
      write(3,300) ((( l, m, i, m = -l,l ), l = 0,lmw ), i = 1,2)
    endif
    do is1 = 1,nspinp
      do lm1 = 1,nlmw
        write(3,310) ( ( Taull(lm1,is1,lm2,is2,ia), lm2 = 1,nlmw ), is2 = 1,nspinp )
      end do
    end do
  endif

  if( ndim == 0 ) then
    Time_tria = 0._db
  else
    call CPU_TIME(time)
    tp3 = real(time,db)
    Time_tria = tp3 - tp2
  endif

  return
  110 format(/' ---- MSM ---------',100('-'))
  115 format(/ ' igrph =',i2,', Matrix dimension =',i7)
  120 format(/' Multiple scattering matrix'/)
  130 format(' iligne =',i3,', ia =',i2,', (l,m,isp) = (',3i2,')' )
  140 format(1p,12(1x,2e11.4))
  150 format(/' Inverse of the multiple scattering matrix'/)
  160 format(/' Scattering amplitude of the central atom (taull):')
  161 format(/' Scattering amplitude of the central atom (taull) before eliminating the contribution of the singular solution:')
  162 format(/' Scattering amplitude of the atom ia =',i3,' (taull):')
  210 format('( l, m, s)',34(8x,3i3,6x))
  212 format('( l, m, s)',34(10x,2i3,7x))
  218 format(3i3,2x,1p,34(1x,2e11.3))
  220 format(8x,a1,2x,102i3)
  230 format(3x,2i3,2x,102a3)
  240 format(3i3,2x,102a3)
  250 format(1p,5(1x,a3,' = ',2e15.7))    
  260 format(/' Taull, Multiple scattering amplitude')
  270 format(2i4,' = lmax, nspinp / ( l, m)')
  280 format(242(13x,2i3,12x))
  290 format(2i4,' = lmax, nspinp / ( l, m, s)')
  300 format(242(12x,3i3,10x))
  310 format(1p,242(1x,2e15.7))

end

!***********************************************************************

! Calcul des coefficients de matrice dus a la projection vers la base symetrisee.
! Il y a peut-etre une erreur dans cette routine --> resultat different
! pour le cuivre cfc en symetrie 4/mmm et mmm (ou 1). A cause de ca, le
! groupe ponctuel est pris de symetrie plus basse en cas de green dans
! igrpt_sg_cal appele par la routine agregat.

subroutine Cmat_cal(Cmat,iato,icheck,igrph,iopsymr,irep_util,is_eq,karact,lato,lmaxg,mato,natome,nb_eq,nb_sym_op,ngrph, &
                nlmsamax,nlmsa,nlmsam,nspino,rot_atom,Ylm_comp)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db), dimension(nb_sym_op,nspino):: C_Gamma
  complex(kind=db), dimension(-lmaxg:lmaxg):: Ctem
  complex(kind=db), dimension(natome,nlmsamax,nb_sym_op,-lmaxg:lmaxg,nspino):: Cmat
  complex(kind=db), dimension(nopsm,nrepm):: karact
  complex(kind=db), dimension(:,:,:), allocatable:: Dlmm_ia, Dlmm_ib,Dlmm_is
  complex(kind=db), dimension(:,:), allocatable:: Mat_ia, Mat_is

  integer, dimension(natome):: nb_eq, nlmsa
  integer, dimension(nopsm):: iopsymr
  integer, dimension(nlmsam,natome,ngrph):: lato, mato, iato
  integer, dimension(nb_sym_op,natome):: is_eq
  integer, dimension(nrepm,2):: irep_util

  logical:: Rotation, Sym_rot, Ylm_comp

  real(kind=db), dimension(3,3):: Mat_rot
  real(kind=db), dimension(3,3,natome):: rot_atom

  do isp = 1,nspino
    js = 0
    irep = abs( irep_util(igrph,isp) )
    if( irep == 0 ) cycle
    do is = 1,nopsm
      if( iopsymr(is) == 0 ) cycle
      js = js + 1
      if( irep_util(igrph,isp) > 0 ) then
        C_Gamma(js,isp) = Conjg( karact(is,irep) )
      else
        C_Gamma(js,isp) = karact(is,irep)
      endif
    end do
  end do
  Cmat(:,:,:,:,:) = (0._db, 0._db)

  do ia = 1,natome
    if( nlmsa(ia) == 0 ) cycle
    lmax = lato(nlmsa(ia),ia,igrph)
    allocate( Dlmm_ia(-lmax:lmax,-lmax:lmax,0:lmax) )
    allocate( Dlmm_ib(-lmax:lmax,-lmax:lmax,0:lmax) )
    allocate( Dlmm_is(-lmax:lmax,-lmax:lmax,0:lmax) )

    if( abs( rot_atom(1,1,ia) - 1._db ) > eps10 .or. abs( rot_atom(2,2,ia) - 1._db ) > eps10 .or. &
        abs( rot_atom(3,3,ia) - 1._db ) > eps10 ) then
      rotation = .true.
    else
      rotation = .false.
    endif
    Mat_rot(:,:) = rot_atom(:,:,ia)
    call Rotation_mat(Dlmm_ia,icheck,Mat_rot,lmax,Ylm_comp)

    if( icheck > 3 ) then
      write(3,100) ia
      do l = 0,lmax
        write(3,102) l
        do m = -l,l
          write(3,104) Dlmm_ia(m,-l:l,l)
        end do
      end do
    endif


! Boucle sur les atomes equivalents
    do ind = 1,nb_eq(ia)
      is = abs( is_eq(ind,ia) )

      js = 0
      do ks = 1,nopsm
        if( iopsymr(ks) == 0 ) cycle
        js = js + 1
        if( ks == is ) exit
      end do

      if( is == 1 .or. is == 25 ) then
        sym_rot = .false.
      else
        sym_rot = .true.
      endif

      call Dlmm_cal_is(icheck,is,lmax,Dlmm_is,Ylm_comp)

! La matrice de la symetrie est mise dans la base de l'atome ia.
! On prend tous les atomes equivalents selon la meme base.
      if( rotation .and. sym_rot ) then
        do l = 0,lmax
          allocate( Mat_is(-l:l,-l:l) )
          allocate( Mat_ia(-l:l,-l:l) )
          Mat_is(-l:l,-l:l) = Dlmm_is(-l:l,-l:l,l)
          Mat_ia(-l:l,-l:l) = Dlmm_ia(-l:l,-l:l,l)
          Mat_is = Matmul( Mat_is, Conjg(Transpose(Mat_ia)) )
          Mat_is = Matmul( Mat_ia, Mat_is )
          Dlmm_is(-l:l,-l:l,l) = Mat_is(-l:l,-l:l)
          deallocate( Mat_is )
          deallocate( Mat_ia )
        end do
      endif

      if( icheck > 3 ) then
        write(3,106) is, ind
        do l = 0,lmax
          write(3,102) l
          do m = -l,l
            write(3,104) Dlmm_is(m,-l:l,l)
          end do
        end do
      endif

      do lm = 1,nlmsa(ia)
        l = lato(lm,ia,igrph)
        m = mato(lm,ia,igrph)
        isp = iato(lm,ia,igrph)
        do mp = -l,l
          Cmat(ia,lm,ind,mp,isp) = Cmat(ia,lm,ind,mp,isp) + C_Gamma(js,isp) * conjg(Dlmm_is(m,mp,l))
        end do
      end do
    end do

! Normalisation
    do lm = 1,nlmsa(ia)
      do isp = 1,nspino
        C_norm = sum( abs( Cmat(ia,lm,:,:,isp) )**2 )
        if( C_norm < eps10 ) cycle
        Cmat(ia,lm,:,:,isp) = Cmat(ia,lm,:,:,isp) / sqrt( C_norm )
      end do
    end do

    if( rotation ) then
      do ind = 1,nb_eq(ia)

! Supprime car tous les atomes 'ind' sont dans la base de 'ia'
!            isym = abs( is_eq(ind,ia) )
!            if( isym /= 1 ) then
!              call Dlmm_cal_is(icheck,isym,lmax,Dlmm_is,Ylm_comp)
!              do l = 0,lmax
!                allocate( Mat_is(-l:l,-l:l) )
!                allocate( Mat_ia(-l:l,-l:l) )
!                Mat_is(-l:l,-l:l) = Dlmm_is(-l:l,-l:l,l)
!                Mat_ia(-l:l,-l:l) = Dlmm_ia(-l:l,-l:l,l)
!                Mat_ia = Matmul( Mat_is, Mat_ia )
!                Dlmm_ib(-l:l,-l:l,l) = Mat_ia(-l:l,-l:l)
!                deallocate( Mat_is )
!                deallocate( Mat_ia )
!              end do
!            else
!              Dlmm_ib(:,:,:) = Dlmm_ia(:,:,:)
!            endif

! Rotation pour se mettre dans la base du cluster (et de la diffusion).
        do lm = 1,nlmsa(ia)
          l = lato(lm,ia,igrph)
          do isp = 1,nspino
            Ctem(-l:l) = Cmat(ia,lm,ind,-l:l,isp)
            do m = -l,l
              Cmat(ia,lm,ind,m,isp) = sum( Ctem(-l:l) * Dlmm_ia(-l:l,m,l) )
            end do
          end do
        end do

      end do

    endif

    deallocate( Dlmm_ia )
    deallocate( Dlmm_ib )
    deallocate( Dlmm_is )

  end do

  if( icheck > 2 ) then
    write(3,110) igrph
    do ia = 1,natome
      write(3,*)
      do lm = 1,nlmsa(ia)
          l = lato(lm,ia,igrph)
          do ind = 1,nb_eq(ia)
            do m = -l,l
              if( sum( abs(Cmat(ia,lm,ind,m,:))) < eps10 ) cycle
              write(3,120) ia, lm, ind, l, m, Cmat(ia,lm,ind,m,:)
            end do
          end do
      end do
    end do
  endif

  return
  100 format(/' Dlmm_ia (Atom =',i2,') :')
!  100 format(/' Atom =',i2,', ind =',i2,', isym =',i3,//' Dlmm_ia :')
  102 format(/' l =',i2)
  104 format(7(1x,2f8.5))
  106 format(/' Dlmm_is ( is =',i3,', ind =',i2,' ) :')
  110 format(/' ia lm  ind  l  m   Cmat = <rep,ia,lm=(l,s)|ind,l,m>,    rep =',i2)
  120 format(2i3,2x,3i3,4f10.5)
end

!***********************************************************************

subroutine Dlmm_cal_is(icheck,is,lmax,Dlmm_is,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, is, js, l, lmax, m, mp

  complex(kind=db), dimension(-lmax:lmax,-lmax:lmax,0:lmax):: Dlmm_is

  logical:: Ylm_comp

  real(kind=db):: Dlmm_cal
  real(kind=db), dimension(3,3):: matopsym

  select case(is)

    Case(1,22,23,24,25,40,41,42)

       do l = 0,lmax
         do m = -l,l
           do mp = -l,l
             Dlmm_is(m,mp,l) = cmplx( Dlmm_cal(is,l,m,mp), 0._db,db)
           end do
         end do
       end do

    Case default

      select case(is)

        case(26,27,28)
          js = is - 7
        case(29,30,31)
          js = is - 13
        case(32,33,34,35,36,37,38,39)
          js = is - 30
        case(43)
          js = 15
        case(44)
          js = 13
        case(45)
          js = 11
        case(46)
          js = 14
        case(47)
          js = 12
        case(48)
          js = 10
        case(53)
          js = 52
        case(54)
          js = 51
        case(55)
          js = 50
        case(56)
          js = 49
        case(57,59,61,63)
          js = is + 1
        case default
          js = is
      end select

      call opsym(js,matopsym)

      call Rotation_mat(Dlmm_is,icheck,matopsym,lmax,Ylm_comp)

! On ajoute l'inversion (changement de signe pour l impair).
      if( is /= js ) then
       do l = 1,lmax,2
         do m = -l,l
           do mp = -l,l
             Dlmm_is(m,mp,l) = - Dlmm_is(m,mp,l)
           end do
         end do
       end do
      endif

  end select

  return
end

!***********************************************************************

function Dlmm_cal(isym,l,m,mp)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  logical:: pair

  Dlmm_cal = 0._db

  if( m /= mp ) return

  pair = .false.
  im = abs( m )

  select case(isym)

    case(1)
      pair = .true.

    case(22)
      mm = mod(l+im,2)
      pair = (mm == 0 .and. m >= 0) .or. (mm == 1 .and. m < 0)

    case(23)
      mm = mod(l+2*im,2)
      pair = (mm == 0 .and. m >= 0) .or. (mm == 1 .and. m < 0)

    case(24)
      mm = mod(im,2)
      pair = mm == 0

    case(25)
      pair = mod(l,2) == 0

    case(40)
      m2 = mod(im,2)
      pair = ( m <= 0 .or. m2 == 0 ) .and. ( m >= 0 .or. m2 == 1 )

    case(41)
      pair = m >= 0

    case(42)
      pair = mod(l+m,2) == 0

  end select

  if( pair ) then
    Dlmm_cal = 1._db
  else
    Dlmm_cal = -1._db
  endif

  return
end

!***********************************************************************

! Developpement en chemin

subroutine dev_chemin(iaabsi,iato,igrph,lato,mat,mato,natome,nchemin,ndim,ngrph,nlmabs,nlmch,nlmagm, &
             nlmsa,nlmsam,nlmsmax,nspinp,taull,taup)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db), dimension(ndim,nlmch):: mat
  complex(kind=db), dimension(nlmch,nlmch):: matp
  complex(kind=db), dimension(nlmsmax,nlmsmax,natome):: taup
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,natome):: taull

  integer, dimension(natome):: nlmsa
  integer, dimension(nlmsam,natome,ngrph):: lato, mato, iato

  iab = iaabsi

  lm0 = 0
  do ia = 1,natome
    if( ia == iaabsi ) exit
    lm0 = lm0 + nlmsa(ia)
  end do

  do lm1 = 1,nlmabs
    l = lato(lm1,iab,igrph)
    m = mato(lm1,iab,igrph)
    lm01 = l**2 + l + 1 + m
    is1 = iato(lm1,iab,igrph)
    do lm2 = 1,nlmabs
      l = lato(lm2,iab,igrph)
      m = mato(lm2,iab,igrph)
      lm02 = l**2 + l + 1 + m
      is2 = iato(lm2,iab,igrph)
      taull(lm01,is1,lm02,is2,iab) = taup(lm1,lm2,iab)
    end do
  end do

  if( nchemin /= 2 ) matp = - mat

  do ichemin = 2,nchemin
    if( nchemin == 2 ) then
      do lm2 = 1,nlmabs
        do lm1 = 1,lm2
          i = nlmabs
          do ia = 2,natome
            do lm3 = 1,nlmsa(ia)    ! boucle sur la ligne
              i = i + 1
              matp(lm1,lm2) = matp(lm1,lm2) + taup(lm3,lm3,iab) * mat(i,lm1) * mat(i,lm2)
            end do
          end do
          matp(lm1,lm2) = taup(lm1,lm1,iab) * matp(lm1,lm2)
        end do
      end do
    else
      matp = - matmul(mat,matp)
    endif

    do lm2 = 1,nlmabs
      l = lato(lm2,iab,igrph)
      m = mato(lm2,iab,igrph)
      lm02 = l**2 + l + 1 + m
      is2 = iato(lm2,iab,igrph)
      lmd2 = lm0 + lm2
      taull(lm02,is2,lm02,is2,iab) = taull(lm02,is2,lm02,is2,iab) + matp(lmd2,lmd2) * taup(lm2,lm2,iab)
      do lm1 = 1,lm2-1
        l = lato(lm1,iab,igrph)
        m = mato(lm1,iab,igrph)
        lm01 = l**2 + l + 1 + m
        is1 = iato(lm1,iab,igrph)
        lmd1 = lm0 + lm1
        taull(lm01,is1,lm02,is2,iab) = taull(lm01,is1,lm02,is2,iab) + matp(lmd1,lmd2) * taup(lm2,lm2,iab)
        taull(lm02,is2,lm01,is1,iab) = taull(lm01,is1,lm02,is2,iab)
      end do
    end do

  end do

  return
end

!***********************************************************************

! Calcul des Dlmm en fonction des matrices de rotations cartesiennes.
! Vient de Didier Sebillaud
! Modifie par Yves Joly en Fevrier 2008

subroutine Rotation_mat(Dlmm,icheck,Mat_rot,lmax,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, l, lmax, m, mp

  real(kind=db):: alfa, beta, gamma, r2

!  Number of N! values computed : Nfac

  complex(kind=db), dimension(-lmax:lmax,-lmax:lmax):: Mat, Trans, Transi
  complex(kind=db), dimension(-lmax:lmax,-lmax:lmax,0:lmax):: Dlmm

  Logical:: Ylm_comp

  real(kind=db), dimension(3,3):: Mat_rot
  real(kind=db), dimension(-lmax:lmax,-lmax:lmax,0:lmax):: R

  call Euler_mat(icheck,Mat_rot,alfa,beta,gamma)

  call DJMN(beta,R,lmax)

  do l = 0,lmax
    do m = -l,l
      do mp = -l,l
        Dlmm(m,mp,l) = exp(-img * ( m*alfa + mp*gamma) ) * R(m,mp,l)
      end do
    end do
  end do

  if( Ylm_comp ) return

! Cas des harmoniques reelles :

  Trans(:,:) = (0._db, 0._db)
  Transi(:,:) = (0._db, 0._db)
  r2 = 1 / sqrt(2._db)
  do m = -lmax,lmax
    if( m > 0 ) then
      Trans(m,m) = cmplx(r2, 0._db, db)
      Trans(m,-m) = cmplx(0._db, r2, db)
    elseif( m == 0 ) then
      Trans(0,0) = (1._db, 0._db)
    else
      if( mod(m,2) == 0 ) then
        Trans(m,-m) = cmplx(r2, 0._db, db)
        Trans(m,m) = cmplx(0._db, -r2, db)
      else
        Trans(m,-m) = cmplx(-r2, 0._db, db)
        Trans(m,m) = cmplx(0._db, r2, db)
      endif
    endif
  end do
  do m = -lmax,lmax
    if( m > 0 ) then
      Transi(m,m) = cmplx(r2, 0._db, db)
      if( mod(m,2) == 0 ) then
        Transi(m,-m) = cmplx(r2, 0._db, db)
      else
        Transi(m,-m) = cmplx(-r2, 0._db, db)
      endif
    elseif( m == 0 ) then
      Transi(0,0) = (1._db, 0._db)
    else
      if( mod(m,2) == 0 ) then
        Transi(m,m) = cmplx(0._db, r2, db)
      else
        Transi(m,m) = cmplx(0._db, -r2, db)
      endif
      Transi(m,-m) = - cmplx(0._db, r2, db)
    endif
  end do

  do l = 0,lmax

    do m = -l,l
      do mp = -l,l
        Mat(m,mp) = sum( Dlmm(m,-l:l,l) * Trans(-l:l,mp) )
      end do
    end do

    do m = -l,l
      do mp = -l,l
        Dlmm(m,mp,l) = sum( Transi(m,-l:l) * Mat(-l:l,mp) )
      end do
    end do

  end do

  return
end

!***********************************************************************

!  This routine calculates Wigner rotation matrices R^{L}_{M1 M2} up to
!             order LMAX, following Messiah's convention.
!                   They are stored as R(M2,M1,L).

Subroutine DJMN(beta,R,LMAX)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  Integer Nfac, EPS0

  Parameter( Nfac = 125 )

  real(kind = db), dimension(-lmax:lmax,-lmax:lmax,0:lmax):: R
  real(kind = db), dimension(0:2*lmax,0:2*lmax,0:2*lmax):: CF
  real(kind = db), dimension(0:2*lmax,0:2*lmax):: expr
  real(kind = db), dimension(0:Nfac):: GLG

  data SMALL,SQR2 /0.001_db,1.414213562_db/

!  Logarithm of the Gamma function : Log(N!)
  GLG(1) = 0._db
  do I = 2,Nfac
    J = I - 1
    GLG(I) = GLG(J) + Log( real(J, db) )
  end do

  NL_M = lmax + 1

!  Coefficients appearing in Wigner's rotation matrices

  expr(0,0) = 1._db
  do L = 1,2*NL_M-2
    do M1 = 0,L
      expr(M1,L) =  exp(0.5_db*(GLG(L+L+1)-GLG(L+M1+1)-GLG(L-M1+1)))
      if( L < NL_M ) then
        do M2 = 0,L
          CF(L,M1,M2) = sqrt( real( (L*L-M1*M1)*(L*L-M2*M2), db ) ) / L
        end do
      endif
    end do
  end do

  C = COS( beta ) * 0.5_db
  S = SIN( beta ) * 0.5_db
  CC = C + C
  CMUL = -1._db
  if( ABS(S) < SMALL ) then
    if( C > 0._db ) EPS0  =  1
    if( C < 0._db ) EPS0  =  -1
    do L = 0,LMAX
      do M1 = -L,L
        do M2 = -L,L
          if( M1 /=  M2 * EPS0 ) then
            R(M2,M1,L) = 0._db
          else
            if( EPS0 == 1 ) then
              R(M2,M1,L) = 1._db
            else
              if( MOD(L+M1,2) == 0 ) then
                R(M2,M1,L) = 1._db
              else
                R(M2,M1,L) = -1._db
              endif
            endif
          endif
        end do
      end do
    end do
  else
    S1 = S * SQR2
    C1 = 0.5_db + C
    R(0,0,0) = 1._db
    R(-1,-1,1) = C1
    R(0,-1,1) = S1
    R(1,-1,1) = 1._db - C1
    R(-1,0,1) = -S1
    R(0,0,1) = CC
    R(1,0,1) = S1
    R(-1,1,1) = 1._db - C1
    R(0,1,1) = -S1
    R(1,1,1) = C1

    PRODL = -S
    COEF = -S / C1
    CL = -1._db
    do L = 2,LMAX
      CL = -CL
      L1 = L - 1
      FLL1 = CC * ( L + L1 )
      FLL2 = 1._db / ( L * L1 * CC )
      PRODL = - PRODL * S

!  Case M = 0

      R_1 = expr(0,L) * PRODL
      R(-L,0,L) = R_1

      R(L,0,L) = R_1 * CL
      R(0,-L,L) = R_1 * CL

      R(0,L,L) = R_1

      CM2 = CL
      do M2 = -L1,-1
        CM2 = CM2 * CMUL
        CF1 = CF(L1,0,-M2) / FLL1
        CF2 = FLL1 / CF(L,0,-M2)
        if( -M2 < L1 ) then
          R_A = CF2 * ( R(M2,0,L1) - R(M2,0,L-2) * CF1 )
        else
          R_A = CF2 * R(M2,0,L1)
        endif

        R(M2,0,L) = R_A

        R(-M2,0,L) = R_A * CM2
        R(0,M2,L) = R_A * CM2

        R(0,-M2,L) = R_A

      end do

      R(0,0,L) = FLL1 * R(0,0,L1) / CF(L,0,0) - R(0,0,L-2) * CF(L1,0,0) / CF(L,0,0)

!  Case M > 0

      PRODM = 1._db
      CM = CL
      FLLM = 0._db
      do M = 1,L1
        CM = -CM
        PRODM = PRODM * COEF
        FLLM = FLLM + FLL2

        R_1 = expr(M,L) * PRODL * PRODM
        R_2 = R_1 / ( PRODM * PRODM )

        R(-L,M,L) = R_1
        R(-L,-M,L) = R_2

        R(L,-M,L) = R_1 * CM
        R(M,-L,L) = R_1 * CM
        R(L,M,L) = R_2 * CM
        R(-M,-L,L) = R_2 * CM

        R(-M,L,L) = R_1
        R(M,L,L) = R_2

        CM2 = CM
        do M2 = -L1,-M
          CM2 = -CM2
          D0 = M2 * FLLM
          CF1 = CF(L1,M,-M2) / FLL1
          CF2 = FLL1 / CF(L,M,-M2)
          if( (M < L1) .AND. (-M2 < L1) ) then
            R_A = CF2 * ( (1._db-D0)*R(M2,M,L1) - R(M2,M,L-2)*CF1 )
            R_B = CF2 * ( (1._db+D0)*R(M2,-M,L1) - R(M2,-M,L-2)*CF1)
          else
            R_A = CF2 * ( 1._db - D0 ) * R(M2,M,L1)
            R_B = CF2 * ( 1._db + D0 ) * R(M2,-M,L1)
          endif

          R(M2,M,L) = R_A
          R(M2,-M,L) = R_B

          R(-M2,-M,L) = R_A * CM2
          R(M,M2,L) = R_A * CM2
          R(-M,M2,L) = R_B * CM2
          R(-M2,M,L) = R_B * CM2

          R(-M,-M2,L) = R_A
          R(M,-M2,L) = R_B

        end do
      end do

      PRODM = PRODM * COEF
      R_1 = PRODL * PRODM
      R_2 = PRODL / PRODM
      R(-L,L,L) = R_1
      R(L,-L,L) = R_1
      R(L,L,L) = R_2
      R(-L,-L,L) = R_2

    end do
  endif

  return
end

! **********************************************************************

! Calcul des angles d'Euler a partir de la matrice de rotation

subroutine Euler_mat(icheck,Mat,alfa,beta,gamma)

  use declarations
  implicit none

  integer:: icheck

  real(kind=db):: alfa, beta, cos_a, cos_b, cos_g, gamma, sin_a, sin_b, sin_g
  real(kind=db), dimension(3,3):: Mat

  cos_b = Mat(3,3)
  sin_b = sqrt( abs(1 - cos_b**2) )    ! beta entre 0 et pi

  if( sin_b > 0.00001_db ) then

    cos_a = Mat(1,3) / sin_b
    sin_a = Mat(2,3) / sin_b
    cos_g = - Mat(3,1) / sin_b
    sin_g = Mat(3,2) / sin_b

  else

    cos_a = Mat(2,2)
    sin_a = - Mat(1,2)
    cos_g = 1._db
    sin_g = 0._db

  endif

  cos_a = min(cos_a, 1._db);  cos_a = max(cos_a, -1._db)
  cos_b = min(cos_b, 1._db);  cos_b = max(cos_b, -1._db)
  cos_g = min(cos_g, 1._db);  cos_g = max(cos_g, -1._db)

  alfa = acos( cos_a )
  if( sin_a < 0._db ) alfa = - alfa
  beta = acos( cos_b )
  if( sin_b < 0._db ) beta = - beta
  gamma = acos( cos_g )
  if( sin_g < 0._db ) gamma = - gamma

  if( icheck > 3 ) write(3,100) alfa / radian, beta / radian, gamma / radian

  return
  100 format(/' alfa, beta, gamma = ',3f10.5)
end

!***********************************************************************

! Ecriture ou lecture des amplitudes de diffusion en cas de calcul en mode One run.

subroutine Data_one_run(iabsm,iaprotoi,icheck,igreq,index_e,igroupi,ipr0,lmax_probe,lmaxa, &
                mpinodes,multi_run,n_atom_proto,n_multi_run,natome,neqm,nge,ngreq,nlm_probe,nlmagm,nspinp, &
                Rot_atom,Rot_int,Spinorbite,taull,taull_stk,Ylm_comp)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: i, ia, icheck, igr, index_e, ipr, ipr0, is1, is2, l, l1, l2, lm1, lm2, lmax_probe, lmw, lmx, m, m1, m2, mpinodes, &
            mu, multi_run, n_atom_proto, n_multi_run, natome, neqm, nge, nlm_probe, nlmagm, nlmw, nspinp
  integer, dimension(natome):: iaprotoi, igroupi, lmaxa
  integer, dimension(n_multi_run):: iabsm
  integer, dimension(0:n_atom_proto):: ngreq
  integer, dimension(0:n_atom_proto,neqm):: igreq

  complex(kind=db), dimension(nlm_probe,nspinp,nlm_probe,nspinp,2:n_multi_run,nge):: taull_stk

  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,natome):: taull
  complex(kind=db), dimension(:,:), allocatable:: Mat, Mat_1, Mat_2
  complex(kind=db), dimension(:,:,:), allocatable:: Dlmm

  logical:: Found, Spinorbite, Ylm_comp

  real(kind=db), dimension(3,3):: Mat_rot, Rot_int
  real(kind=db), dimension(3,3,natome):: Rot_atom

  if( index_e == 1 .and. mpinodes == 1 ) then
    if( multi_run == 1 ) then
      Open(50, status = 'SCRATCH')
    else
      Rewind(50)
    endif
  endif

  do mu = 2,n_multi_run
    Found = .false.

    do ipr = ipr0,n_atom_proto
      if( igreq(ipr,1) == iabsm(mu) ) exit
    end do

    do ia = 1,natome

      if( iaprotoi(ia) /= ipr ) cycle

      do igr = 1,ngreq(ipr)
        if( igreq(ipr,igr) == igroupi(ia) ) then
          Found = .true.
          exit
        endif
      end do
      if( .not. Found ) cycle

      lmx = min(lmaxa(ia), lmax_probe)

      if( abs( rot_atom(1,1,ia) - 1._db ) > eps10 .or. abs( rot_atom(2,2,ia) - 1._db ) > eps10 .or. &
        abs( rot_atom(3,3,ia) - 1._db ) > eps10 ) then
        Mat_rot(:,:) = rot_atom(:,:,ia)
        Mat_rot = Matmul( Rot_int, Mat_rot )
      else
        Mat_rot = Rot_int
      endif

! En ecriture, on applique la rotation inverse --> Transpose
      if( multi_run == 1 ) Mat_rot = Transpose( Mat_rot )

      Allocate( Dlmm(-lmx:lmx,-lmx:lmx,0:lmx) )
      call Rotation_mat(Dlmm,icheck,Mat_rot,lmx,Ylm_comp)

      do is1 = 1,nspinp
        do is2 = 1,nspinp
          if( .not. spinorbite .and. is1 /= is2 ) cycle

          do l1 = 0,lmx
            lm1 = l1**2 + l1 + 1
            Allocate( Mat_1(-l1:l1,-l1:l1) )
            Mat_1(-l1:l1,-l1:l1) = Dlmm(-l1:l1,-l1:l1,l1)
            Mat_1 = Transpose( Mat_1 )

            do l2 = 0,lmx
              lm2 = l2**2 + l2 + 1
              Allocate( Mat_2(-l2:l2,-l2:l2) )
              Mat_2(-l2:l2,-l2:l2) = Dlmm(-l2:l2,-l2:l2,l2)
              Allocate( Mat(-l1:l1,-l2:l2) )

              if( multi_run == 1 ) then

                do m1 = -l1,l1
                  do m2 = -l2,l2
                    Mat(m1,m2) = taull(lm1+m1,is1,lm2+m2,is2,ia)
                  end do
                end do
                Mat = Matmul( Mat_1, Matmul( Mat, Mat_2) )
                if( mpinodes == 1 ) then
                  write(50,*) Mat(:,:)
                else
                  do m1 = -l1,l1
                    do m2 = -l2,l2
                      taull_stk(lm1+m1,is1,lm2+m2,is2,mu,index_e) = Mat(m1,m2)
                    end do
                  end do
                endif

              else

                if( mpinodes == 1 ) then
                  Read(50,*) Mat(:,:)
                else
                  do m1 = -l1,l1
                    do m2 = -l2,l2
                      Mat(m1,m2) = taull_stk(lm1+m1,is1,lm2+m2,is2,mu,index_e)
                    end do
                  end do
                endif
                Mat = Matmul( Mat_1, Matmul( Mat, Mat_2) )
                do m1 = -l1,l1
                  do m2 = -l2,l2
                    taull(lm1+m1,is1,lm2+m2,is2,ia) = Mat(m1,m2)
                  end do
                end do

              endif

              Deallocate( Mat, Mat_2 )
            end do

            Deallocate( Mat_1 )
          end do

        end do
      end do

      deallocate( Dlmm )

      exit
    end do

    if( .not. Found ) then
      call write_error
      do ipr = 3,9,3
        write(ipr,100) mu
      end do
      stop
    endif

  end do

  if( icheck > 1 .and. multi_run > 1 ) then

    do ia = 1,natome
      if( igroupi(ia) /= iabsm(multi_run) ) cycle

      nlmw = min(nlmagm,9)
      write(3,110)
      write(3,120)
      if( nspinp == 2 ) then
        write(3,210) ((( l, m, i, i=1,2), m = -l,l), l = 0,2 )
      else
        write(3,212) (( l, m, m = -l,l), l = 0,2 )
      endif

      lm1 = 0
      do l1 = 0,lmx
        do m1 = -l1,l1
          lm1 = lm1 + 1
          do is1 = 1,nspinp
             write(3,218) l1, m1, is1, (( taull(lm1,is1,lm2,is2,ia), is2 = 1,nspinp ), lm2 = 1,nlmw)
          end do
        end do
      end do

    end do

  endif

! Ecriture pour extract
  if( icheck > 0 .and. multi_run > 1 ) then
    do ia = 1,natome
      if( igroupi(ia) /= iabsm(multi_run) ) cycle
    
    ! si lmaxa > 10 il faut changer le format d'ecriture
      lmw = min( 10, lmaxa(ia) )
      nlmw = min( nlmagm, (lmw+1)**2 )
      write(3,300)
      do is1 = 1,nspinp
        do lm1 = 1,nlmw
          write(3,310) ( ( Taull(lm1,is1,lm2,is2,ia), lm2 = 1,nlmw ), is2 = 1,nspinp )
        end do
      end do
    end do
  endif

  return
  100 format(//' Prototypical atom number ',i3,' not in the sphere of calculation !'//&
               ' Increase the cluster radius !'/)
  110 format(/' ---- Data_One_Run ',100('-'))
  120 format(/' Multiple scattering amplitude, Taull:')
  210 format('( l, m, s)',50(7x,3i3,6x))
  212 format('( l, m, s)',50(10x,2i3,7x))
  218 format(3i3,2x,1p,50(1x,2e11.3))
  300 format(/' Taull, Multiple scattering amplitude')
  310 format(1p,242(1x,2e15.7))
end
