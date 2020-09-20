! FDMNES subroutines
! File containing all concerning the solutions in the atomic sphere

!***********************************************************************

! Resolution of the spherical Schrodinger equation in the atoms

subroutine Sphere(Axe_Atom_grn,Ecinetic,Eimag,Energ,Enervide,Full_atom, &
            Full_potential,Green,Hubb_a,Hubb_d,iaabsi,iapr,iaprotoi, &
            ibord,icheck,igreq,igroupi,iopsymr,lmax,lmax_pot,m_hubb, n_atom_0, &
            n_atom_ind,n_atom_proto,n_comp,natome,nbord,nbtm,nbtm_fdm,neqm,ngroup_m,nlm_pot,nlmagm, &
            nlmmax,nphiato1,nphiato7,npsom,nr,nspin,nspino,nspinp,numat,phiato,posi,r,Relativiste,Renorm,Rmtg, &
            Rmtsd,Spinorbite,Tau_ato,V_hubb,V_intmax,V0bd,Vrato,xyz,Ylm_comp,Ylmato)

  use declarations
  implicit none

  integer:: ia, iaabsi, iang, iapr, ib, icheck, iga, igr, kgr, L, &
    l_hubbard, lfin, lm, lm0, lmax, lmax_pot, lmp, m, m_hubb, mp, n, &
    n_atom_0, n_atom_ind, n_atom_proto, n_comp, natome, nbtm, nbtm_fdm, neqm, ngroup_m, nlm, nlm_pot, &
    nlm1, nlm2, nlmagm, nlmmax, np, nphiato1, nphiato7, npsom, nr, nrmtg, nrmtsd, nspin, nspino, nspinp, numat
  integer, dimension(natome):: iaprotoi, igroupi, nbord
  integer, dimension(nopsm):: iopsymr
  integer, dimension(nbtm,natome):: ibord
  integer, dimension(0:n_atom_proto,neqm):: igreq

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,n_atom_0:n_atom_ind):: Tau_ato
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  logical:: Classic_irreg, Ecomp, Full_atom, Full_potential, Green, Hubb_a, Hubb_d, Hubb_m, Radial_comp, Relativiste, Renorm, &
    Spinorbite, Ylm_comp

  real(kind=db):: cosang, Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_grn
  real(kind=db), dimension(nbtm_fdm,nlmmax,natome):: Ylmato

  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nr,nspin):: g0, gm, gmi, gp, gpi
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(nphiato1,n_comp*nlmagm,nspinp,nspino,natome,nphiato7):: phiato

  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur, usi, usr

  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag, db) )

  if( icheck > 1 ) write(3,110) iapr, numat, lmax
  if( icheck > 1 ) write(3,120) Energ*rydb, Enervide*rydb
  if( icheck > 2 ) then
    write(3,130) Ecinetic(:)*rydb
    write(3,140) V0bd(:)*rydb
    write(3,150) konde(:)
  endif

  if( Full_potential ) then
    nlm = ( lmax + 1 )**2
  else
    nlm = 1
  endif
  allocate( V(nr,nlm,nlm,nspin) )

  call mod_V(icheck,lmax,lmax_pot,nlm,nlm_pot,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)

  gpi(1:nr-1,:) = 1 / gp(1:nr-1,:)

  if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10 .or. Ecinetic(nspin) < eps10 ) then
    Ecomp = .true.
  else
    Ecomp = .false.
  endif
  Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

  if( Full_potential ) then
    lfin = 0
  else
    lfin = lmax
  endif

  do L = 0,lfin

    if( Hubb_a .and. L == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = nlm   ! = (lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*L + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*L + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,L,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

! Copy
    if( Full_potential ) then
      do lm = 1,nlm1
        do lmp = 1,nlm2
          Tau_ato(lm,:,lmp,:,iapr) = Tau(lm,:,lmp,:)
        end do
      end do
    else
      lm0 = L**2 + L + 1
       do m = -L,L
        if( nlm1 == 1 ) then
          n = 1
        else
          n = L + 1 + m
        endif
        do mp = -L,L
          if( nlm1 == 1 ) then
            if( m /= mp ) cycle
            np = 1
          else
            np = L + 1 + mp
          endif
          Tau_ato(lm0+m,:,lm0+mp,:,iapr) = Tau(n,:,np,:)
        end do
      end do
    endif

    deallocate( Tau )

 ! Calculation of the radial functions phiato on the FDM grid points
    if( .not. Green ) then

      allocate( usi(nr,nlm1,nlm2,nspinp,nspino) )
      allocate( usr(nr,nlm1,nlm2,nspinp,nspino) )

! gm and gp reverse in subroutine
      if( n_comp == 2 ) then
        Classic_irreg = .false.
        gmi(1:nr-1,:) = 1 / gm(1:nr-1,:)
        call Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gmi,gp,gso,Hubb_a,Hubb_d,icheck,konde, &
          L,Lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtg,nr,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtg,Spinorbite,Tau,usi,usr,V,V_hubb)
      endif

      do ib = 1,natome
        if( Full_atom ) then
          if( ib > 1 ) exit
          ia = iapr
        else
          ia = ib
          if( iaprotoi(ia) /= iapr ) cycle
        endif

        if( nspin == 2 ) then
          igr = igroupi(ia)
          iga = igroupi(iaabsi)
          kgr = igreq(iaprotoi(ia),1)

          cosang = sum( Axe_Atom_grn(:,kgr) * Axe_Atom_grn(:,igr) )
          cosang = cosang * abs( sum( Axe_Atom_grn(:,iga) * Axe_Atom_grn(:,igr) ) )
          if( abs(cosang - 1) < eps4 ) then
            iang = 1
          elseif( abs(cosang + 1) < eps4 ) then
            iang = - 1
          else
            iang = 0
          endif
        else
          iang = 1
        endif

        call cal_phiato(Full_potential,ia,iang,ibord,icheck,iopsymr,L,lmax,n_comp,nlm1,nlm2,natome, &
             nbord(ia),nbtm,nbtm_fdm,nlmagm,nlmmax,nphiato1,nphiato7,npsom,nr,nspinp,nspino,phiato,posi,r, &
             Radial_comp,Spinorbite,ui,ur,usi,usr,xyz,Ylm_comp,Ylmato)
      end do
      
      deallocate( usi, usr )
      
    endif

    deallocate( ui, ur )

  end do   ! end of loop over L

  deallocate( V )

  return
  110 format(/' iapr =',i3,', Z =',i3,', lmax =',i2)
  120 format(' Energ =',f10.3,' eV,  Enervide =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format(' V0bd  =',2f10.3)
  150 format(' konde =',2f12.5)
end

!***********************************************************************

Subroutine mod_V(icheck,lmax,lmax_pot,nlm,nlm_pot,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, ir, ispin, L1, L2, l3, lm, lm1, lm2, lm3, lmax, &
    lmax_pot, m1, m2, m3, nlm, nlm_pot, nlma, nr, nrmtg, nrmtsd, nspin

  logical:: Ylm_comp

  real(kind=db):: g, Gaunt_r, gauntcp, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nspin):: V0bd
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato
  real(kind=db), dimension(nr,nlm,nlm,nspin):: V

  V(:,:,:,:) = 0._db

  do ispin = 1,nspin
    do lm = 1,nlm
      V(1:nr,lm,lm,ispin) = Vrato(1:nr,1,ispin)
    end do
  end do

  if( V_intmax < 1000._db ) then
    do ispin = 1,nspin
      do lm = 1,nlm
        do ir = 1,nr
          V(ir,lm,lm,ispin) = min( V(ir,lm,lm,ispin), V_intmax )
        end do
      end do
    end do
  endif

  do ir = 1,nr-1
    if( r(ir) > Rmtg + eps10 ) exit
  end do
  nrmtg = ir
  do ir = 1,nr-1
    if( r(ir) > Rmtsd + eps10 ) exit
  end do
  nrmtsd = ir

! No effect for FDM
  if( abs( V0bd(1) ) > eps10 ) then
    do ir = 1,nr-1
      if( r(ir) > Rmtg - eps10 ) exit
    end do
    do ispin = 1,nspin
      do lm = 1,nlm
        V(ir:nr,lm,lm,ispin) = V0bd(ispin)
      end do
    end do
  endif

  lm1 = 0
  do L1 = 0,lmax
    do m1 = -L1,L1
      lm1 = lm1 + 1
      lm2 = 0
      do L2 = 0,lmax
        do m2 = -L2,L2
          lm2 = lm2 + 1
          if( lm1 == lm2 .or. nlm == 1 ) cycle

          lm3 = 1
          do l3 = 1,lmax_pot
            do m3 = -l3,l3
              lm3 = lm3 + 1

              if( Ylm_comp ) then
                g = gauntcp(L1,m1,L2,m2,l3,m3)
              else
                g = Gaunt_r(L1,m1,L2,m2,l3,m3)
              endif

              V(1:nr,lm1,lm2,:) = V(1:nr,lm1,lm2,:) + g * Vrato(1:nr,lm3,:)
            end do
          end do

        end do
      end do
    end do
  end do

  if( icheck > 2 ) then
    nlma = min( nlm, 9 )
    if( nspin == 2 .and. nlm_pot == 1 ) then
      write(3,'(/A)') '    Radius     V(up)      V(dn)'
    elseif( nspin == 2 .and. nlm_pot > 1 ) then
      write(3,'(/A,A)') '    Radius  V(up,lm)   V(dn,lm)    lm = 1,nlm_pot'
    elseif( nspin == 1 .and. nlm_pot == 1 ) then
      write(3,'(/A)') '    Radius      V'
    else
      write(3,100) ( ( lm1, lm2, lm2 = 1,nlma ), lm1 = 1,nlma )
    endif
    do ir = 1,nr
      write(3,110) r(ir)*bohr, ( ( V(ir,lm1,lm2,:)*rydb, lm2 = 1,nlma ), lm1 = 1,nlma )
    end do
  endif

  return
  100 format('    Radius    ',1000('V(',i2,',',i2,')',3x))
  110 format(1p,1000e11.3)
end

!***********************************************************************

subroutine coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)

  use declarations
  implicit none

  integer:: ir, ispin, nlm, nr, nspin, nspino, numat

  logical:: Relativiste, Spinorbite

  real(kind=db):: a2s4, bder, dr, dvr, Enervide, fac, r0, rm, rp, Vme

  real(kind=db), dimension(nr):: cgrad0, cgradm, cgradp, clapl0, claplm, claplp, f2, r
  real(kind=db), dimension(nr,nspin):: g0, gm, gp
  real(kind=db), dimension(nr,nlm,nlm,nspin):: V
  real(kind=db), dimension(nr,nspino):: gso

  gso(:,:) = 0._db

  do ir = 1,nr
    if( ir == nr ) then
      rp = r(nr)**2 / r(nr-1)
    else
      rp = r(ir+1)
    endif
    if( ir == 1 ) then
      rm = r(1)**2 / r(2)
    else
      rm = r(ir-1)
    endif
    r0 = r(ir)
    dr = 0.5_db * ( rp - rm )
    claplm(ir) = 1 / ( ( r0 - rm ) * dr )
    claplp(ir) = 1 / ( ( rp - r0 ) * dr )
    clapl0(ir) = - claplm(ir) - claplp(ir)
    if( Spinorbite .or. Relativiste ) then
      cgradm(ir) = ( rp - r0 ) / ( ( rm - r0 ) * ( rp - rm ) )
      cgradp(ir) = ( rm - r0 ) / ( ( rp - r0 ) * ( rm - rp ) )
      cgrad0(ir) = - cgradm(ir) - cgradp(ir)
    endif
    f2(ir) = 1 / r(ir)**2
  end do

! alfa_sf = fine structure constant
  a2s4 = 0.25_db * alfa_sf**2

  do ir = 1,nr
    do ispin = 1,nspin

      Vme = V(ir,1,1,ispin) - Enervide
      g0(ir,ispin) = - clapl0(ir) + Vme
      gm(ir,ispin) = - claplm(ir)
      gp(ir,ispin) = - claplp(ir)

      if( Relativiste .or. Spinorbite ) then

        bder = 1 / ( 1 - a2s4 * Vme )
        if( ir == 1 ) then
          dvr = 2 * numat / r(1)**2
        elseif( ir == nr ) then
          dvr = cgradm(ir) * V(ir-1,1,1,ispin) + cgrad0(ir) * V(ir,1,1,ispin) + cgradp(ir) * V(ir,1,1,ispin)
        else
          dvr = cgradm(ir) * V(ir-1,1,1,ispin) + cgrad0(ir) * V(ir,1,1,ispin) + cgradp(ir) * V(ir+1,1,1,ispin)
        endif
        fac = a2s4 * bder * dvr

        if( Relativiste ) then
          g0(ir,ispin) = g0(ir,ispin) - a2s4 * Vme**2 - fac * ( cgrad0(ir) - 1 / r(ir) )
          gm(ir,ispin) = gm(ir,ispin) - fac * cgradm(ir)
          gp(ir,ispin) = gp(ir,ispin) - fac * cgradp(ir)
        endif
        if( Spinorbite ) gso(ir,ispin) = fac / r(ir)

      endif

    end do
  end do

  if( nspin < nspino ) gso(:,nspino) = gso(:,nspin)

  return
end

!***********************************************************************

! Resolution of the radial Schrodinger equation. Include the case of the non spherical potential
! Called by Sphere, radial, radial_sd
! In ur(ir,m,mp,isp,isol) : (mp, isol) define the basis state
!                           (m, isp) are the real moment and spin

subroutine Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,LL,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

  use declarations
  implicit none

  integer:: icheck, im, ip, ir, isol, isp, isq, isr, iss, L, l_hubbard, L2, li, LL, lf, lmax, m, m_hubb, mp, ms, &
    n, nlm, nlm1, nlm2, np, nr, nrmtg, nsol, nspin, nspino, nspinp, numat

  logical:: Ecomp, Failed, Full_potential, Full_potential_loc, Hubb_a, Hubb_d, Hubb_m, Hubb_nd, &
    Radial_comp, Renorm, Relativiste, Spino_simple, Spinorbite

  complex(kind=db):: V_h
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlm1,nspinp,nlm1,nspinp):: Tau
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  real(kind=db):: br, ci, Eimag, fac, p, Rmtg, td
  real(kind=db), dimension(nspin):: cr, Ecinetic
  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nr,nspin):: g0, gm, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr,nlm,nlm,nspin):: V
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur
  real(kind=db), dimension(-LL:LL,2):: fac_m

  Full_potential_loc = Full_potential
  Failed = .false.

! Point de retour si la normalisation en cas de full-potential a les elements non diagonaux plus grands que les diagonaux
  do

    ur(:,:,:,:,:) = 0._db
    ui(:,:,:,:,:) = 0._db

    Spino_simple = Spinorbite .and. nlm2 == 1

    if( Full_potential ) then
      li = 0  
      lf = lmax
    else
      li = LL
      lf = LL
    endif

  do L = li,lf

    L2 = L * ( L + 1 )

    if( Hubb_a .and. L == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif
    Hubb_nd = Hubb_m .and. .not. Hubb_d

    ci = - Eimag / ( 4 * L + 6 )

    br = - numat / ( L + 1._db )
    cr(:) = - ( 2 * numat * br + Ecinetic(:) ) / ( 4 * L + 6 )

! Je fais L'hypothese que le terme de Hubbard ne croise pas les solutions 1 et 2 en cas de spinorbite.

! Development at origin
    do isol = 1,nspino
      do isp = 1,nspinp
        isr = min( isp, nspin )
        do m = -L,L

          if( Full_potential ) then
            n = L**2 + L + 1 + m
          elseif( Hubb_m .or. Spinorbite ) then
            n = L + 1 + m
          else
            n = 1
          endif

          if( (m == L .and. isp == 1) .or. (m == -L .and. isp == 2 ) ) then
            nsol = 1
          else
            nsol = 2
          endif
          if( Spinorbite .and. nsol == 1 .and. isp /= isol ) cycle

          if( Spinorbite .and. Relativiste ) then
            if( nsol == 1  .or. isol == 1 ) then
              p = sqrt( ( L + 1 )**2 - ( alfa_sf * numat )**2 )
            else
              p = sqrt( L**2 - ( alfa_sf * numat )**2 )
            endif
          elseif( Spinorbite ) then
            if( nsol == 1 ) then
              p = L + 1._db
            elseif( isol == 1 ) then
              p = 0.5_db + 0.5_db * sqrt( 1._db + 4*(L**2) + 8*L )
            else
              if( L == 1 ) then
! In fact, there is no solution for L=1 with spin-orbit and not relativistic !
!                      p = 1._db * L
                p = L + 1._db
              else
                p = 0.5_db + 0.5_db * sqrt( -5._db + 4*(L**2) )
              endif
            endif
          elseif( Relativiste ) then
            p = sqrt( L**2 + L + 1 - ( alfa_sf * numat )**2 )
          else
            p = L + 1._db
          endif

 ! Becareful, it is from the point of vew of the m of the orbital
          if( Spinorbite .and. nsol /= 1 ) then
            if( isol == 1 .and. isp == 1 ) then
              fac = sqrt( ( L + m + 1 ) / ( 2*L + 1._db ) )
            elseif( isol == 1 .and. isp == 2 ) then
              fac = sqrt( ( L - m + 1 ) / ( 2*L + 1._db ) )
            elseif( isol == 2 .and. isp == 1 ) then
              fac = sqrt( ( L - m ) / ( 2*L + 1._db ) )
            else
              fac = - sqrt( ( L + m ) / ( 2*L + 1._db ) )
            endif
          else
            fac = 1._db
          endif
          if( nlm2 == 1 ) then
            np = 1
          elseif( .not. Spinorbite .or. isp == isol ) then
            np = n
          elseif( isp < isol ) then
            np = n + 1
          else
            np = n - 1
          endif
          ur(1:2,n,np,isp,isol) = fac * ( 1 + br * r(1:2) + cr(isr) * r(1:2)**2 ) * r(1:2)**p
          if( Radial_comp ) ui(1:2,n,np,isp,isol) = fac * ci * r(1:2)**(p+2)
        end do
      end do
    end do
  end do

  if( Spinorbite .and. .not. Full_potential ) then
    L = LL
    do m = -L,L
      fac_m(m,1) = sqrt( ( L - m ) * ( L + m + 1._db ) )
      fac_m(m,2) = sqrt( ( L + m ) * ( L - m + 1._db ) )
    end do
  endif

  do ir = 2,nr-1
    im = ir - 1
    ip = ir + 1

    do isol = 1,nspino
      do isp = 1,nspinp
        isr = min( isp, nspin )
        isq = 3 - isp
        if( Spinorbite ) then
          iss = isol
        else
          iss = isp
        endif
        
        do L = li,lf

          L2 = L * ( L + 1 )

          if( Hubb_a .and. L == l_hubbard( numat ) )  then
            Hubb_m = .true.
          else
            Hubb_m = .false.
          endif

          do m = -L,L

            if( Spino_simple .and. ( ( m == L .and. isp == 1 ) .or. ( m == -L .and. isp == 2 ) ) .and. isp /= isol ) cycle
  
            if( Full_potential ) then
              n = L**2 + L + 1 + m
            elseif( Hubb_m .or. Spinorbite ) then
              n = L + 1 + m
            else
              n = 1
              if( m /= 0 ) cycle
            endif

            td = g0(ir,isr) + L2 * f2(ir)

            if( Spinorbite ) then
              if( isp == 1 ) then
                ms = m
                mp = m + 1   ! m de L'autre spin
              else
                ms = - m
                mp = m - 1
              endif
              if( Full_potential ) then
                fac = sqrt( ( L - ms ) * ( L + ms + 1._db ) )
              else
                fac = fac_m(m,isp)
              endif
              td = td + ms * gso(ir,isp)
            else
              mp = L + 1 ! pour que le if en dessous soit faux
            endif
            if( Full_potential ) then
              np = L**2 + L + 1 + mp
            elseif( Hubb_nd .or. Spinorbite ) then
              np = L + 1 + mp
            else
              np = 1
            endif

! The second index "m" is the non zero component at origin
            ur(ip,n,:,isp,isol) = td * ur(ir,n,:,isp,isol) + gm(ir,isr) * ur(im,n,:,isp,isol)
            if( Radial_comp ) then
              ui(ip,n,:,isp,isol) = td * ui(ir,n,:,isp,isol) + gm(ir,isr) * ui(im,n,:,isp,isol)
              if( Ecomp ) then 
                ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol) + Eimag * ui(ir,n,:,isp,isol)
                ui(ip,n,:,isp,isol) = ui(ip,n,:,isp,isol) - Eimag * ur(ir,n,:,isp,isol)
              endif
            endif

            if( abs(mp) <= L ) then
              ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol) + fac * gso(ir,isp) * ur(ir,np,:,isq,isol)
              if( Radial_comp ) ui(ip,n,:,isp,isol) = ui(ip,n,:,isp,isol) + fac * gso(ir,isp) * ui(ir,np,:,isq,isol)
            endif

            if( Hubb_m ) then
              do mp = -L,L
                if( Hubb_d .and. m /= mp ) cycle
                if( Full_potential ) then
                  np = L**2 + L + 1 + mp
                else
                  np = L + 1 + mp
                endif
                V_h = V_hubb(m,mp,iss,iss) 
                if( Radial_comp ) then
                  ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol) + real( V_h, db) * ur(ir,np,:,isp,isol) &
                                                            - aimag( V_h ) * ui(ir,np,:,isp,isol)
                  ui(ip,n,:,isp,isol) = ui(ip,n,:,isp,isol) + real( V_h, db) * ui(ir,np,:,isp,isol) &
                                                            + aimag( V_h ) * ur(ir,np,:,isp,isol)
                else
                  ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol) + real( V_h, db) * ur(ir,np,:,isp,isol)
                endif
              end do
            endif

            if( Full_potential_loc ) then
              do np = 1,nlm
                if( n == np ) cycle
                ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol) + V(ir,n,np,isr) * ur(ir,np,:,isp,isol)
                if( Radial_comp ) ui(ip,n,:,isp,isol) = ui(ip,n,:,isp,isol) + V(ir,n,np,isr) * ui(ir,np,:,isp,isol)
              end do
            endif

            ur(ip,n,:,isp,isol) = - ur(ip,n,:,isp,isol) * gp(ir,isr)
            if( Radial_comp )  ui(ip,n,:,isp,isol) = - ui(ip,n,:,isp,isol) * gp(ir,isr)

          end do
        end do
      end do
    end do
  end do

  if( icheck > 3 ) then
    call write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspinp,nspino,numat,r,Radial_comp,10000._db,ui,ur,1)
  elseif( icheck > 2 ) then 
    call write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspinp,nspino,numat,r,Radial_comp,Rmtg,ui,ur,1)
  endif

  do n = 1,nlm1
    do np = 1,nlm2
      do isp = 1,nspinp
        do isol = 1 ,nspino
          ur(:,n,np,isp,isol) = ur(:,n,np,isp,isol) / r(:)
          if( Radial_comp )  ui(:,n,np,isp,isol) = ui(:,n,np,isp,isol) / r(:)
        end do
      end do
    end do
  end do

  if( .not. Renorm ) return

  call Renormal(Radial_comp,Full_potential,icheck,konde,LL,lmax,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,r,Rmtg, &
                Tau,ui,ur,Failed)

    if( Failed ) then
      Full_potential_loc = .false.
      Failed = .false.
    else
      exit
    endif

  end do ! retour au debut

  if( icheck > 2 ) then
    do ir = 1,nr
      ur(ir,:,:,:,:) = ur(ir,:,:,:,:) * r(ir)
      if( Radial_comp ) ui(ir,:,:,:,:) = ui(ir,:,:,:,:) * r(ir)
    end do
    if( icheck > 3 ) then
      call write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspinp,nspino,numat,r,Radial_comp,10000._db,ui,ur,2)
    else 
      call write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspinp,nspino,numat,r,Radial_comp,Rmtg,ui,ur,2)
    endif
    do ir = 1,nr
      ur(ir,:,:,:,:) = ur(ir,:,:,:,:) / r(ir)
      if( Radial_comp ) ui(ir,:,:,:,:) = ui(ir,:,:,:,:) / r(ir)
    end do
  endif

  return
end

!***********************************************************************

! Writing of the radial function

subroutine write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspinp,nspinpo,numat,r,Radial_comp,Rmtg,ui,ur,icom)

  use declarations
  implicit none

  integer:: i, icom, ir, isp, isol, k, kmax, L, li, lf, lp, m, nlm1, nlm2, mp, n, np, nr, nspinp, nspinpo, numat

  Character(len=14):: mot, mot14
  Character(len=14), dimension(196):: Label, Labeli

  logical:: Full_potential, Radial_comp

  real(kind=db):: Rmtg
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspinpo):: ui, ur

  if( nlm1 == 1 .and. nspinp == 1 ) then
    Label(1) = '     ur       '
    Labeli(1) = '     ui       '
    kmax = 1
  else
    k = 0
    do isol = 1,nspinpo
      np = 0
      do lp = li,lf
        do mp = -lp,lp
          if( nlm2 == 1 .and. mp /= 0 ) cycle
          do isp = 1,nspinp
            k = k + 1
            mot14 = ' ur(          '
            i = 5
            if( nspinp == 2 ) then
              if( isp == 1 ) then
                mot14(i:i) = 'u'
              else
                mot14(i:i) = 'd'
              endif
              i = i + 1
              if( nlm2 > 1 ) then
                mot14(i:i) = ','
                i = i + 1
              endif
            endif
            if( nlm2 > 1 ) then
              if( Full_potential ) then
                call ad_number(abs(lp),mot14,14)
                i = i + 1
                mot14(i:i) = ','
                i = i + 1
              endif
              if( mp < 0 ) then
                mot14(i:i) = '-'
                i = i + 1
              endif
              call ad_number(abs(mp),mot14,14)
              i = i + 1
            endif
            if( nspinpo == 2 .and. nlm2 == 1 ) then
              mot14(i:i) = ','
              i = i + 1
              if( isol == 1 ) then
                mot14(i:i) = 'u'
              else
                mot14(i:i) = 'd'
              endif
              i = i + 1
            endif
            mot14(i:i) = ')'
            label(k) = mot14
            if( Radial_comp ) then
              mot14(3:3) = 'i'
              labeli(k) = mot14
            endif
          end do
        end do
      end do
      if( nlm2 > 1 ) exit
    end do
    kmax = k
    do k = 1,kmax
      mot = label(k)
      Call Center_word( mot, 14 )
      label(k) = mot
      if( Radial_comp ) then
        mot = labeli(k)
        Call Center_word( mot, 14 )
        labeli(k) = mot
      endif
   end do
  endif

  n = 0
  np = 0
  do L = li,lf

    if( nlm2 == 1 ) then
      do m = -L,L
        if( nlm1 == 1 .and. m /= 0 ) cycle
        if( nlm1 == 1 ) then
          if( icom == 1 ) then
            write(3,110) numat, L
          elseif( icom == 2 ) then
            write(3,120) numat, L
          else
            write(3,130) numat, L
          endif
        else
          if( icom == 1 ) then
            write(3,140) numat, L, m
          elseif( icom == 2 ) then
            write(3,150) numat, L, m
          else
            write(3,160) numat, L, m
          endif
        endif
        n = n + 1
        if( nspinpo == 2 ) then
          if( Radial_comp ) then
            write(3,170)
          else
            write(3,180)
          endif
        endif
        if( Radial_comp ) then
          write(3,190) ( Label(k), Labeli(k), k = 1,kmax )
          do ir = 1,nr
            write(3,200) r(ir)*bohr, ( ( ( ur(ir,n,np,isp,isol), ui(ir,n,np,isp,isol), isp = 1,nspinp), np = 1,nlm2), &
                                             isol = 1,nspinpo )
            if( r(ir) > Rmtg + eps10 ) exit
          end do
        else
          write(3,190) ( Label(k), k = 1,kmax )
          do ir = 1,nr
            write(3,200) r(ir)*bohr, ( ( ( ur(ir,n,np,isp,isol), isp = 1,nspinp), np = 1,nlm2), isol = 1,nspinpo )
            if( r(ir) > Rmtg + eps10 ) exit
          end do
        endif
      end do
    else
      do mp = -L,L
        np = np + 1
        do isol = 1,nspinpo
          if( nspinpo == 1 ) then
            if( icom == 1 ) then
              write(3,140) numat, L, mp
            elseif( icom == 2 ) then
              write(3,150) numat, L, mp
            else
              write(3,160) numat, L, mp
            endif
          else
            if( icom == 1 ) then
              write(3,210) numat, L, mp, isol
            elseif( icom == 2 ) then
              write(3,220) numat, L, mp, isol
            else
              write(3,230) numat, L, mp, isol
            endif
          endif
          if( Radial_comp ) then
            write(3,190) ( Label(k), Labeli(k), k = 1,kmax )
            do ir = 1,nr
              write(3,200) r(ir)*bohr, ( ( ur(ir,n,np,isp,isol), ui(ir,n,np,isp,isol), isp = 1,nspinp), n = 1,nlm1 )
              if( r(ir) > Rmtg + eps10 ) exit
            end do
          else
            write(3,190) ( Label(k), k = 1,kmax )
            do ir = 1,nr
              write(3,200) r(ir)*bohr, ( ( ur(ir,n,np,isp,isol), isp = 1,nspinp), n = 1,nlm1 )
              if( r(ir) > Rmtg + eps10 ) exit
            end do
          endif
        end do
      end do
    endif
  end do

  return
  110 format(/' Radial wave function time r:  Z =',i3,', L =',i2,/)
  120 format(/' Radial wave function, time r, after normalization:', '  Z =',i3,', L =',i2,/)
  130 format(/' Radial singular function time r:  Z =',i3,', L =',i2,/)
  140 format(/' Radial wave function time r:  Z =',i3,', L =',i2, ', m =',i2,/)
  150 format(/' Radial wave function, time r, after normalization:', '  Z =',i3,', L =',i2,', m =',i2,/)
  160 format(/' Radial singular function time r:  Z =',i3,', L =',i2, ', m =',i2,/)
  170 format(37x,' Solution 1',46x,' Solution 2')
  180 format(23x,' Solution 1',18x,' Solution 2')
  190 format('     Radius    ',392a14)
  200 format(1p,491e14.6)
  210 format(/' Radial wave function time r:  Z =',i3,', L =',i2, ', m =',i2,', Solution',i2,/)
  220 format(/' Radial wave function, time r, after normalization:', '  Z =',i3,', L =',i2,', m =',i2,', Solution',i2,/)
  230 format(/' Radial singular function time r:  Z =',i3,', L =',i2, ', m =',i2,', Solution',i2,/)
end

!***********************************************************************

! Normalisation of radial functions using continuity at muffin-tin radius with solutions in vaccum
! Called by Sch_radial

subroutine Renormal(Radial_comp,Full_potential,icheck,konde,LL,lmax,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,r, &
                  Rmtg,Tau,ui,ur,Failed)

  use declarations
  implicit none

  integer icheck, i, inr, ir, isol, isp, isr, L, LL, lmax, lp, m, nlm1, nlm2, mp, n, np, nr, nrmtg, nspin, nspino, nspinp

  complex(kind=db):: z
  complex(kind=db), dimension(0:lmax):: bess, d_hank, d_bess, hank
  complex(kind=db), dimension(nspin):: konde, s1, s2
  complex(kind=db), dimension(nlm1,nspinp,nlm1,nspinp)::  Tau
  complex(kind=db), dimension(nlm1,nlm2,nspinp,nspino):: Ampl, u1, u2, Wronske, Wronsks
  complex(kind=db), dimension(2,0:lmax,nspin):: bs, hs

  logical:: Failed, Full_potential, Radial_comp

  real(kind=db):: a, ai, b, bi, c, ci, d0, d01, d02, d1, d2, d12, Rmtg, Wronskout
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur

  inr = min( nrmtg, nr )
  d12 = r(inr-1) - r(inr-2)
  d01 = r(inr) - r(inr-1)
  d02 = r(inr) - r(inr-2)
  d2 = 1 / ( d12 * d02 )
  d1 = 1 / ( d12 * d01 )
  d0 = 1 / ( d01 * d02 )

  Wronske(:,:,:,:) = (0._db, 0._db)
  Wronsks(:,:,:,:) = (0._db, 0._db)
  u1(:,:,:,:) = (0._db, 0._db)
  u2(:,:,:,:) = (0._db, 0._db)

! Hankel function and its derivative (versus r, not z) at Rmtg
  do isp = 1,nspin
    call Cal_Hankel(d_hank,hank,konde(isp),lmax,Rmtg,.false.)
    hs(1,0:lmax,isp) = hank(0:lmax)
    hs(2,0:lmax,isp) = d_hank(0:lmax)
    call Cal_Bessel(d_bess,bess,konde(isp),lmax,Rmtg)
    bs(1,0:lmax,isp) = bess(0:lmax)
    bs(2,0:lmax,isp) = d_bess(0:lmax)
  end do

  do isp = 1,nspinp
    do isol = 1,nspino
      do m = 1,nlm1
        do mp = 1,nlm2

          a = ur(inr-2,m,mp,isp,isol) * d2 - ur(inr-1,m,mp,isp,isol) * d1 + ur(inr,m,mp,isp,isol) * d0
          b = ( ur(inr-1,m,mp,isp,isol) - ur(inr-2,m,mp,isp,isol)) / d12 - a * ( r(inr-1) + r(inr-2) )
          c = ur(inr,m,mp,isp,isol) - a * r(inr)**2 - b * r(inr)

          if( Radial_comp ) then
            ai = ui(inr-2,m,mp,isp,isol) * d2 - ui(inr-1,m,mp,isp,isol) * d1 + ui(inr,m,mp,isp,isol) * d0
            bi = (ui(inr-1,m,mp,isp,isol) - ui(inr-2,m,mp,isp,isol)) / d12  - ai * ( r(inr-1) + r(inr-2) )
            ci = ui(inr,m,mp,isp,isol) - ai * r(inr)**2 - bi *r(inr)
          endif
          
          if( Radial_comp ) then
            u1(m,mp,isp,isol) = cmplx( a * Rmtg**2 + b * Rmtg + c,  ai * Rmtg**2 + bi * Rmtg + ci, db )
            u2(m,mp,isp,isol) = cmplx( 2 * a * Rmtg + b, 2 * ai * Rmtg + bi, db )
          else
            u1(m,mp,isp,isol) = cmplx( a * Rmtg**2 + b * Rmtg + c, 0._db, db )
            u2(m,mp,isp,isol) = cmplx( 2 * a * Rmtg + b, 0._db, db )
          endif

        end do
      end do
    end do
  end do

! The following formula for Wronskout comes from normalisation by sqrt(k/pi) of Bessel and Hankel functions
  Wronskout = 1 / ( pi * Rmtg**2 )

  np = 0
  do L = 0,lmax
    if( .not. Full_potential .and. L /= LL ) cycle
    s1(:) = - img * hs(1,L,:)
    s2(:) = - img * hs(2,L,:)
    do mp = -L,L
      if( nlm2 == 1 .and. mp > -L ) exit
      np = np + 1
      do isp = 1,nspinp
        isr = min( isp, nspin )
        if( nlm2 == 1 ) then
          Wronsks(:,1,isp,:) = u1(:,1,isp,:) * s2(isr) - u2(:,1,isp,:) * s1(isr)
          Wronske(:,1,isp,:) = u1(:,1,isp,:) * bs(2,L,isr) - u2(:,1,isp,:) * bs(1,L,isr)
        else
! Dans ce cas le deuxieme indice de "u" est celui de L'attaque
! c'est aussi le premier indice de Wronsk
! et L est celui sortant
          Wronsks(np,:,isp,:) = u1(np,:,isp,:) * s2(isr) - u2(np,:,isp,:) * s1(isr)
          Wronske(np,:,isp,:) = u1(np,:,isp,:) * bs(2,L,isr) - u2(np,:,isp,:) * bs(1,L,isr)
        endif
      end do
    end do
  end do

  if( icheck > 2 ) then
    n = 0
    do L = 0,lmax
      if( .not. Full_potential .and. L /= LL ) cycle
      write(3,110) L
      write(3,120) - img * hs(1,L,:)
      write(3,130) - img * hs(2,L,:)
      write(3,140) bs(1,L,:)
      write(3,150) bs(2,L,:)
      write(3,160) Wronskout
      if( nspino == 2 ) then
        write(3,162)
      elseif( nspinp == 2 ) then
        write(3,164)
      else
        write(3,166)
      endif
      do m = -L,L
        if( nlm1 == 1 .and. m /= 0 ) cycle
        n = n + 1
        np = 0
        do lp = 0,lmax
          if( .not. Full_potential .and. lp /= LL ) cycle
          do mp = -lp,lp
            if( nlm2 == 1 .and. mp /= 0 ) cycle
            np = np + 1
            write(3,170) L, m, lp, mp, ' Wronsks   =', ( Wronsks(n,np,:,isol), isol = 1,nspino )
            write(3,170) L, m, lp, mp, ' u1        =', ( u1(n,np,:,isol), isol = 1,nspino )
            write(3,170) L, m, lp, mp, ' u2        =', ( u2(n,np,:,isol), isol = 1,nspino )
          end do
        end do
      end do
    end do
  endif

  call cal_ampl(Ampl,Full_potential,icheck,LL,lmax,nlm1,nlm2,nspino,nspinp,Tau,Wronske,Wronsks,Wronskout,Failed)

  if( Failed ) return

  do ir = 1,nr

    if( Radial_comp ) then
      u1(:,:,:,:) =  cmplx( ur(ir,:,:,:,:), ui(ir,:,:,:,:), db )
      ur(ir,:,:,:,:) = 0._db
      ui(ir,:,:,:,:) = 0._db
    else
      u1(:,:,:,:) =  cmplx( ur(ir,:,:,:,:), 0._db, db )
      ur(ir,:,:,:,:) = 0._db
    endif

    if( nlm2 == 1 .and. nspino == 1 ) then

      do isp = 1,nspinp
        do n = 1,nlm1
          z = Ampl(n,1,isp,1) * u1(n,1,isp,1)
          ur(ir,n,1,isp,1) = real( z, db )
          if( Radial_comp ) ui(ir,n,1,isp,1) = aimag( z )
        end do
      end do

    elseif( nlm2 == 1 .and. nspino == 2 ) then

      do isp = 1,nspinp
        n = 0
        do m = -LL,LL
          n = n + 1
          do isol = 1,nspino ! le spin d'attaque devient L'indice de solution

            mp = m + isol - isp
            if( mp > LL .or. mp < -LL ) cycle
            np = n + isol - isp

            z = sum( Ampl(np,1,isol,:) * u1(n,1,isp,:) )

            ur(ir,n,1,isp,isol) = real( z, db )
            if( Radial_comp ) ui(ir,n,1,isp,isol) = aimag( z )

          end do
        end do
      end do

    else

      do isp = 1,nspinp
        do n = 1,nlm1
          do isol = 1,nspino ! le spin d'attaque devient L'indice de solution
            do np = 1,nlm2

! np and n indeed in that way !
              if( nspino == 1 ) then
                z = sum( Ampl(np,:,isp,isol) * u1(n,:,isp,isol) )
              else
                z = ( 0._db, 0._db )
                do i = 1,nspino
                  z = z + sum( Ampl(np,:,isol,i) * u1(n,:,isp,i) )
                end do
              endif

              ur(ir,n,np,isp,isol) = real( z, db )
              if( Radial_comp )  ui(ir,n,np,isp,isol) = aimag( z )

            end do
          end do
        end do
      end do

    endif

  end do

  return
  110 format(/' Renormalisation parameters, L =',i2)
  120 format(/' -i.hankel =',1p,8e11.3)
  130 format('     deriv =',1p,8e11.3)
  140 format('    bessel =',1p,8e11.3)
  150 format('     deriv =',1p,8e11.3)
  160 format(' Wronskout =',1p,8e11.3)
  162 format('  L  m lp mp',21x,'up up',17x,'dn up',17x,'up dn',17x, 'dn dn')
  164 format('  L  m lp mp',  23x,'up',19x,'dn')
  166 format('  L  m lp mp')
  170 format(4i3,a12,1p,16e11.3)
end

!***********************************************************************

! Calculation of the amplitudes of the radial functions

subroutine cal_ampl(Ampl,Full_potential,icheck,LL,lmax,nlm1,nlm2,nspino,nspinp,Tau,Wronske,Wronsks,Wronskout,Failed)

  use declarations
  implicit none

  integer:: i, icheck, isol, isp, ispp, j, L, LL, lmax, lp, m, nlm1, nlm2, mp, mpv, mq, n, nd, ndim, np, nspino, nspinp, nu

  complex(kind=db):: Det
  complex(kind=db), dimension(nlm1,nspinp,nlm1,nspinp):: Tau
  complex(kind=db), dimension(nlm1,nlm2,nspinp,nspino):: Amp, Ampl, Wronske, Wronsks
  complex(kind=db), dimension(:,:), allocatable:: Mat_A, Mat_e, Mat_s, Mat_T

  logical:: Failed, Full_potential, Stop_job

  real(kind=db):: Wronskout

  Failed = .false.
  Stop_job = .true.

  Amp(:,:,:,:) = (0._db, 0._db)
  Ampl(:,:,:,:) = (0._db, 0._db)
  Tau(:,:,:,:) = (0._db, 0._db)

  if( nlm2 == 1 ) then
    if( nspino == 1 ) then
      do isp = 1,nspinp
        do n = 1,nlm1
          Ampl(n,1,isp,1) = Wronskout / Wronske(n,1,isp,1)
          Tau(n,isp,n,isp) = - Wronske(n,1,isp,1) / Wronsks(n,1,isp,1)
        end do
      end do
    else
! Pour Ampl, isp est le spin d'attaque (donc L'indice de solution apres la renormalisation)
! Le m de Amp est m + 1/2 - is = m + isol - 1
      do m = -LL-1,LL
        if( m == -LL-1 ) then
          isp = 2
          isol = 2
          mq = m + isp - 1
          n = LL + 1 + mq
          Ampl(n,1,isp,isol) = Wronskout / Wronske(n,1,isp,isol)
          Tau(n,isp,n,isp) = - Wronske(n,1,isp,isol) / Wronsks(n,1,isp,isol)
        elseif( m == LL ) then
          isp = 1
          isol = 1
          mq = m + isp - 1
          n = LL + 1 + mq
          Ampl(n,1,isp,isol) = Wronskout / Wronske(n,1,isp,isol)
          Tau(n,isp,n,isp) = - Wronske(n,1,isp,isol) / Wronsks(n,1,isp,isol)
        else
          nu = LL + 1 + m
          nd = LL + 1 + m + 1
          Det = Wronsks(nu,1,1,1) * Wronsks(nd,1,2,2) - Wronsks(nd,1,2,1) * Wronsks(nu,1,1,2)
          Amp(nu,1,1,1) = Wronsks(nd,1,2,2) / Det
          Amp(nu,1,1,2) = - Wronsks(nd,1,2,1) / Det
          Amp(nd,1,2,1) = - Wronsks(nu,1,1,2) / Det
          Amp(nd,1,2,2) = Wronsks(nu,1,1,1) / Det

          do isp = 1,nspinp
            n = LL + 1 + m + isp - 1
            do ispp = 1,nspinp
              np = LL + 1 + m + ispp - 1
              Tau(n,isp,np,ispp) = (0._db,0._db)
              do isol = 1,nspino
                Tau(n,isp,np,ispp) = Tau(n,isp,np,ispp) - Amp(n,1,isp,isol) * Wronske(np,1,ispp,isol)
              end do
            end do
          end do

          Det = Wronske(nu,1,1,1) * Wronske(nd,1,2,2) - Wronske(nd,1,2,1) * Wronske(nu,1,1,2)
          Ampl(nu,1,1,1) = Wronskout * Wronske(nd,1,2,2) / Det
          Ampl(nu,1,1,2) = - Wronskout * Wronske(nd,1,2,1) / Det
          Ampl(nd,1,2,1) = - Wronskout * Wronske(nu,1,1,2) / Det
          Ampl(nd,1,2,2) = Wronskout * Wronske(nu,1,1,1) / Det
        endif
      end do
    endif

  else

    ndim = nspino * nlm2
    allocate( Mat_e(ndim,ndim) )
    allocate( Mat_s(ndim,ndim) )
    allocate( Mat_A(ndim,ndim) )
    allocate( Mat_T(ndim,ndim) )

    do isp = 1,nspinp ! spin d'attaque
      if( nspino == 1 .or. isp == 1 ) then
        Mat_e(:,:) = ( 0._db, 0._db )
        Mat_s(:,:) = ( 0._db, 0._db )
        Mat_A(:,:) = ( 0._db, 0._db )
        Mat_T(:,:) = ( 0._db, 0._db )
        i = 0
      endif

      n = 0
      do L = 0,lmax
        if( .not. Full_potential .and. L /= LL ) cycle
        do m = -L,L
          n = n + 1
          i = i + 1

          j = 0
          do ispp = 1,nspinp  ! solution
            if( nspino == 1 .and. ispp /= isp ) cycle
            isol = min(nspino,ispp)
            do lp = 0,lmax
              if( .not. Full_potential .and. lp /= LL ) cycle
              do mp = -lp,lp
!                if( nspino == 2 ) then
!                  mpv = mp - isp + isol
!                  if( abs(mpv) > lp ) cycle
!                else 
                  mpv = mp
!                endif
                if( Full_potential ) then
                  np = 1 + lp + lp**2 + mpv
                else
                  np = lp + 1 + mpv
                endif  
                j = j + 1 
! colonne harmonique reelle, ligne harmonique attaque
!                Mat_e(i,j) = Wronske(n,np,isp,isol)
!                Mat_s(i,j) = Wronsks(n,np,isp,isol)

                Mat_e(j,i) = Wronske(n,np,isp,isol)
                Mat_s(j,i) = Wronsks(n,np,isp,isol)
              end do
            end do
          end do
        end do
      end do

      if( nspino == 2 .and. isp == 1 ) cycle

      if( icheck > 3 ) then
        if(nspino == 1 ) then
          write(3,110) 'incoming', isp
        else
          write(3,120) 'incoming'
        endif
        do i = 1,ndim
          write(3,130) i, Mat_e(i,:)
        end do
        if(nspino == 1 ) then
          write(3,110) 'outgoing', isp
        else
          write(3,120) 'outgoing'
        endif
        do i = 1,ndim
          write(3,130) i, Mat_s(i,:)
        end do
      endif

      call invcomp(ndim,Mat_s,ndim,ndim,0,Stop_job)

      Mat_T = - matmul( Mat_s, Mat_e )

      call invcomp(ndim,Mat_e,ndim,ndim,0,Stop_job)
! le signe est plus car Wronske est en fait "- Wronskien"
      Mat_A = Wronskout * Mat_e

      if( nspino == 1 ) then
        Ampl(:,:,isp,1) = Mat_A(:,:)
        Tau(:,isp,:,isp) = Mat_T(:,:)
      endif

    end do

    if( nspino == 2 ) then

      i = 0
      do isp = 1,nspinp
        n = 0
        do L = 0,lmax
          if( .not. Full_potential .and. L /= LL ) cycle
          do m = -L,L
            n = n + 1
            i = i + 1

            j = 0
            do isol = 1,nspino
              np = 0
              do lp = 0,lmax
                if( .not. Full_potential .and. lp /= LL ) cycle
                do mp = -lp,lp
                  np = np + 1
                  j = j + 1
                  Ampl(n,np,isp,isol) = Mat_A(i,j)
                  Tau(n,isp,np,isol) = Mat_T(i,j)
                end do
              end do
            end do
          end do
        end do
      end do

    endif

    deallocate( Mat_A, Mat_e, Mat_s, Mat_T )

  endif

  if( icheck > 2 ) then
    if( LL == 0 .or. icheck > 1 ) then
      if( nlm1 == 1 .and. nspinp == 1 ) then
        write(3,140)
      elseif( nlm1 == 1 ) then
        write(3,150)
      elseif( nlm2 == 1 .and. nspino == 2 ) then
        write(3,155)
      elseif( nspinp == 1 ) then
        write(3,160)
      else
        write(3,170)
      endif
    endif
    do isp = 1,nspinp
      n = 0
      do L = 0,lmax
        if( .not. Full_potential .and. L /= LL ) cycle
        do m = -L,L
          if( nlm1 == 1 .and. m /= 0 ) cycle
          n = n + 1
          if( nlm1 == 1 .and. nspinp == 1 ) then
            write(3,220) L, Ampl(n,1,isp,min(isp,nspino))
          elseif( nlm1 == 1 ) then
            write(3,230) L, isp, Ampl(n,1,isp,min(isp,nspino))
          elseif( nlm2 == 1 .and. nspinp == 1 ) then
            write(3,230) L, m, Ampl(n,1,isp,:)
          elseif( nlm2 == 1 ) then
            write(3,240) L, m, isp, Ampl(n,1,isp,:)
          elseif( nspinp == 1 ) then
            write(3,230) L, m, ( Ampl(n,:,isp,ispp), ispp = 1,nspino )
          else
            write(3,240) L, m, isp, ( Ampl(n,:,isp,ispp), ispp = 1,nspino )
          endif
        end do
        if( nlm2 == 1 .and. nspino == 1 ) exit
      end do
    end do
  endif

  if( icheck > 1 ) then
    if( nlm1 == 1 .and. nspinp == 1 ) then
      write(3,180)
    elseif( nlm1 == 1 ) then
      write(3,190)
    elseif( nspinp == 1 ) then
      write(3,200)
    else
      write(3,210)
    endif
    do isp = 1,nspinp
      n = 0
      do L = 0,lmax
        if( .not. Full_potential .and. L /= LL ) cycle
        do m = -L,L
          if( nlm1 == 1 .and. m /= 0 ) cycle
          n = n + 1
          if( nlm1 == 1 .and. nspinp == 1 ) then
            write(3,220) L, Tau(n,isp,n,isp)
          elseif( nlm1 == 1 ) then
            write(3,230) L, isp, Tau(n,isp,n,isp)
          elseif( nlm2 == 1 .and. nspinp == 1 ) then
            write(3,230) L, m, Tau(n,isp,n,isp)
          elseif( nlm2 == 1 .and. nspino == 1 ) then
            write(3,240) L, m, isp, Tau(n,isp,n,isp)
          else
            write(3,240) L, m, isp, ( Tau(n,isp,:,ispp), ispp = 1,nspinp )
          endif
        end do
        if( nlm2 == 1 .and. nspino == 1 ) exit
      end do
    end do
  endif

! Test on the calculation
  if( Full_potential ) then
    do n = 1,nlm1
      do np = n+1,nlm1
        do isp = 1,nspinp
          if( abs( Tau(n,isp,np,isp) ) < abs( Tau(n,isp,n,isp) ) .and. abs( Tau(n,isp,np,isp) ) < abs( Tau(n,isp,n,isp) ) ) &
            cycle
          if( icheck > 0 ) write(3,250) n, np
          Failed = .true.
        end do
      end do
    end do
  endif

  return
  110 format(/' Wronskian ',a8,' matrix, isp =',i2)
  120 format(/' Wronskian ',a8,' matrix')
  130 format(i3,1p,100(2e11.3,1x))
  140 format(/' Amplitude',/'  L',11x,'Ampl')
  150 format(/' Amplitude',/'  L isp',9x,'Ampl')
  155 format(/' Amplitude',/'  L  m isp',6x,'Ampl(Sol 1)',12x, 'Ampl(Sol 2)')
  160 format(/' Amplitude',/'  L  m',10x,'Ampl')
  170 format(/' Amplitude',/'  L  m isp',10x,'Ampl')
  180 format(/'  L',11x,'Tau')
  190 format(/'  L isp',10x,'Tau')
  200 format(/'  L  m',11x,'Tau')
  210 format(/'  L  m isp',10x,'Tau')
  220 format(i3,1p,16(1x,2e13.5))
  230 format(2i3,1p,16(1x,2e13.5))
  240 format(3i3,1p,16(1x,2e13.5))
  250 format(/'  Off-digonal Tau higher than diagonal for', ' lm =',i3,', lmp =',i3,'. Off-diagonal potential neglected !')
end

!**********************************************************************

! Calculation of radial integrals and irregular solution
! Called by Tenseur_car and S_nrixs_cal

subroutine Radial_integral(Classic_irreg,Ecinetic,Eimag,Energ,Enervide,Eseuil,Final_tddft,Full_potential,Hubb_a,Hubb_d,icheck, &
         initlv,ip_max,ip0,lmax,lmax_pot,m_hubb,n_comp,nbseuil,ninit1,ninitv,nlm_pot,nlma,nlma2,nr,NRIXS,nrm,nspin,nspino, &
         nspinp,numat,psii,r,r_or_bess,Relativiste,Renorm,Rmtg,Rmtsd,rof,Singul,Solsing,Spinorbite,V_hubb,V_intmax, &
         V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer:: i, icheck, initl, initlv, ip, ip_max, ip0, ip1, ip2, iseuil, isp, L, l_hubbard, L1, L2, lfin, lm, lm1, &
    lm2, lmax, lmax_pot, lmp, lp, m, m1, m2, m_hubb, &
    n_comp, nlm1, nlm2, mp, nbseuil, ninit1, ninitv, nlm, nlm_pot, nlma, nlma2, nr, nrm, &
    nrmtsd, nrmtg, nspin, nspino, nspinp, numat

  character(len=2):: mot1, mot2
  character(len=104):: mot

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlma,nspinp,nlma,nspinp,ip0:ip_max,ip0:ip_max,ninitv):: Singul
  complex(kind=db), dimension(n_comp*nlma,nlma2,nspinp,nspino,ip0:ip_max,ninitv):: rof
  complex(kind=db), dimension(nlma,nlma2,nspinp,nspino,ip0:ip_max,ninitv):: rofb
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau

  logical:: Classic_irreg, Ecomp, Final_tddft, Full_potential, Hubb_a, &
    Hubb_d, Hubb_m, NRIXS, Radial_comp, Relativiste, Renorm, Solsing, Spinorbite, Ylm_comp

  real(kind=db):: Eimag, Energ, Enervide, Ephoton, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato

  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nbseuil):: Eseuil, Vecond
  real(kind=db), dimension(nr,nspin):: g0, gm, gmi, gpi, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,ip0:ip_max):: r_or_bess

  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur, usi, usr

  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

  if( nbseuil == 2 ) then
    if( Final_tddft .and. nbseuil /= ninitv ) then
      ninit1 = ( ninitv - 2 ) / 2
      if( initlv <= ninit1 ) then
        iseuil = 1
      else
        iseuil = 2
      endif
    else ! not used in DFT (only for TDDFT)
      iseuil = initlv
    endif
  else
    iseuil = 1
  endif

  if( icheck > 1 ) then
    write(3,110) initlv, iseuil
    write(3,120) Energ*rydb
    write(3,130) Ecinetic(:)*rydb
    write(3,140) V0bd(:)*rydb
    write(3,150) konde(:)
  endif

  if( Full_potential ) then
    nlm = ( lmax + 1 )**2
  else
    nlm = 1
  endif
  allocate( V(nr,nlm,nlm,nspin) )

  do i = 1,nbseuil
    if( Final_tddft ) then
      if( i /= nbseuil ) cycle
      Ephoton = Energ + Eseuil(nbseuil)
    else
      Ephoton = Energ + Eseuil(i)
    endif
    Ephoton = max( 0.001_db, Ephoton ) ! Optic case
! Multiplicative term for quadrupolar transitions
! In S.I. vecond = k = E*alfa_sf*4*pi*epsilon0 / (e*e)
! In a.u. and Rydberg : k = 0.5_db * alfa_sf * E
    Vecond(i) = 0.5_db * alfa_sf * Ephoton
  end do

  call mod_V(icheck,lmax,lmax_pot,nlm,nlm_pot,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)

  gpi(:,:) = 1 / gp(:,:)
  if( Solsing .or. n_comp == 2 ) gmi(1:nr-1,:) = 1 / gm(1:nr-1,:)

  if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10 .or. Ecinetic(nspin) < eps10 ) then
    Ecomp = .true.
  else
    Ecomp = .false.
  endif
!  Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

  if( Full_potential ) then
    lfin = 0
  else
    lfin = lmax
  endif

  do L = 0,lfin

    if( Hubb_a .and. L == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif
    Radial_comp = Ecomp .or. ( Hubb_m .and. Ylm_comp )

    if( Full_potential ) then
      nlm1 = nlm    ! = ( lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*L + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*L + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2, Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde, &
         L,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

! Radial integral for the golden rule
    call radial_matrix(Final_tddft,initlv,ip_max,ip0,iseuil,L,n_comp,nlm1,nlm2,nbseuil, &
           ninitv,nlma,nlma2,nr,NRIXS,nrm,nspino,nspinp,psii,r,r_or_bess,Radial_comp,Rmtsd,rof,ui,ur,Vecond)

! Calculation of irregular solution
    if( Solsing ) call Cal_Solsing(Classic_irreg,Ecomp,Eimag,f2,Final_tddft,Full_potential,g0,gmi,gp,gso,Hubb_a,Hubb_d, &
         icheck,initlv,ip_max,ip0,iseuil,konde,L,lmax,m_hubb,nlm, &
         nlm1,nlm2,nbseuil,ninitv,nlma,nr,NRIXS,nrm,nrmtsd,nspin,nspino,nspinp, &
         numat,psii,r,r_or_bess,Radial_comp,Rmtsd,Singul,Spinorbite,Tau,ui,ur,V,V_hubb,Vecond)

    allocate( usi(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( usr(nr,nlm1,nlm2,nspinp,nspino) )

! gm and gp reverse in subroutine
    if( n_comp == 2 ) then
      Classic_irreg = .false.
      call Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gmi,gp,gso,Hubb_a,Hubb_d,icheck,konde, &
          L,Lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtg,nr,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtg,Spinorbite,Tau,usi,usr,V,V_hubb)

      call radial_matrix(Final_tddft,initlv,ip_max,ip0,iseuil,L,1,nlm1,nlm2,nbseuil, &
           ninitv,nlma,nlma2,nr,NRIXS,nrm,nspino,nspinp,psii,r,r_or_bess,Radial_comp,Rmtg,rofb,usi,usr,Vecond)
      do lm = 1,nlma
        rof(lm+nlma,:,:,:,:,:) = rofb(nlma,:,:,:,:,:) 
      end do

    endif

    deallocate( Tau, ui, ur, usi, usr )

  end do   ! end of loop over L

  if( icheck > 1 ) then
    mot = ' '
    if( Ecomp .and. nspino == 2 ) then
      mot(11:17) = 'up sol1'
      mot(37:43) = 'up sol2'
      mot(63:69) = 'dn sol1'
      mot(89:95) = 'dn sol2'
    elseif( nspino == 2 ) then
      mot(5:50) = 'up sol1      up sol2      dn sol1      dn sol2'
    elseif( nspinp == 2 ) then
      mot = '      up           dn'
    else
      mot = ' '
    endif

    do i = 1,nbseuil
      if( Final_tddft ) then
        if( i /= nbseuil ) cycle
        initl = initlv
      else
        initl = i
        if( nbseuil > 1 ) write(3,155) i
      endif

      do ip = ip0,ip_max
        select case(ip)
          case(0)
            write(3,160) '  monopole'
          case(1)
            write(3,160) '    dipole'
          case(2)
            write(3,160) 'quadrupole'
          case(3)
            write(3,160) '  octupole'
        end select
        if( Full_potential .or. Hubb_a ) then
          write(3,170) '  L  m lp mp', mot
        else
          write(3,175) '  L  m', mot
        endif
        lm = 0
        do L = 0,lmax
          do m = -L,L
            lm = lm + 1
            if( .not. ( Hubb_a .or. Full_potential .or. Spinorbite ) .and. m /= 0 ) cycle
            lmp = 0
            do lp = 0,lmax
              if( .not. Full_potential .and. L /= lp ) cycle
              do mp = -lp,lp
                if( .not. ( Hubb_a .or. Full_potential .or. m == mp ) ) cycle
                lmp = lp**2 + lp + 1 + mp
                lmp = min(nlma2,lmp)
                if( Full_potential .or. Hubb_a ) then
                  if( Ecomp) then
                    write(3,180) L, m, lp, mp, ( rof(lm,lmp,isp,:,ip,initl), isp = 1,nspinp)
                  else
                    write(3,180) L, m, lp, mp, ( real(rof(lm,lmp,isp,:,ip,initl)), isp = 1,nspinp)
                  endif
                else
                  if( Ecomp) then
                    write(3,185) L, m, ( rof(lm,lmp,isp,:,ip,initl), isp = 1,nspinp )
                  else
                    write(3,185) L, m, ( real(rof(lm,lmp,isp,:,ip,initl)), isp = 1,nspinp )
                  endif
                endif
              end do
            end do
          end do
        end do
      end do

      if( .not. Solsing ) cycle

      do ip1 = ip0,ip_max
        select case(ip1)
          case(0)
            mot1 = 'M1'
          case(1)
            mot1 = 'E1'
          case(2)
            mot1 = 'E2'
          case(3)
            mot1 = 'E3'
        end select

        do ip2 = ip0,ip_max
          if( .not. Full_potential .and. mod(ip1+ip2,2) == 1 ) cycle

          select case(ip2)
            case(0)
              mot2 = 'M1'
            case(1)
              mot2 = 'E1'
            case(2)
              mot2 = 'E2'
            case(3)
              mot2 = 'E3'
          end select

          write(3,190) mot1, mot2
          if( Hubb_m .or. Full_potential .and. nspinp == 1 ) then 
            write(3,200) ' L1 m1 L2 m2','  Sing '
          elseif( Hubb_m .or. Full_potential .or. Spinorbite ) then 
            write(3,200) ' L1 m1 L2 m2','  up up','  up dn','  dn up','  dn dn'
          elseif( nspinp == 2 ) then
            write(3,210) '  L  m','   up  ','   dn  '
          else
            write(3,210) '  L  m','  Sing '
          endif

          lm1 = 0
          do L1 = 0,lmax
            do m1 = -L1,L1
              lm1 = lm1 + 1
              lm2 = 0
              do L2 = 0,lmax
                do m2 = -L2,L2
                  lm2 = lm2 + 1
                  if( Hubb_m .or. Full_potential .or. Spinorbite ) then
                    if( sum( abs( Singul(lm1,:,lm2,:,ip1,ip2,initl) ) ) < eps10 ) cycle
                    write(3,220) L1, m1, L2, m2, ( Singul(lm1,isp,lm2,:,ip1,ip2,initl), isp = 1,nspinp )
                  elseif( lm1 /= lm2 .or. m1 /= 0 ) then
                    cycle
                  else
                    write(3,230) L1, m1, ( Singul(lm1,isp,lm2,isp,ip1,ip2,initl), isp = 1,nspinp )
                  endif
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  endif

  deallocate( V )

  return
  110 format(/' initl =',i2,', iseuil =',i2)
  120 format(' Energ    =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format(' V0bd     =',2f10.3)
  150 format(' konde    =',2f12.5)
  155 format(/' iseuil =',i2)
  160 format(/' Radial matrix, ',a10,1x,'term:')
  170 format(a12,a104)
  175 format(a6,a104)
  180 format(4i3,1p,8e13.5)
  185 format(2i3,1p,8e13.5)
  190 format(/' Radial integral of singular solution, ',a2,'-',a2, ' term:')
  200 format(a12,10x,a7,4(20x,a7))
  210 format(a6,10x,a7,2(20x,a7))
  220 format(4i3,1p,4(1x,2e13.5)) 
  230 format(2i3,1p,4(1x,2e13.5))
end

!**********************************************************************

! Calculation of the radial integral in the (Fermi) Golden rule.
! psii, u and ur are wave functions time r.

subroutine radial_matrix(Final_tddft,initlv,ip_max,ip0,iseuil,L,n_comp,nlm1,nlm2,nbseuil, &
           ninitv,nlma,nlma2,nr,NRIXS,nrm,nspino,nspinp,psii,r,r_or_bess,Radial_comp,Rmtsd,rof,ui,ur,Vecond)

  use declarations
  implicit none

  integer initlv, ip, ip_max, ip0, iseuil, is, isol, isp, iss, L, lm, lm1, lm2, n_comp, nlm1, nlm2, n1, n2, nbseuil, &
    ninitv, nlma, nlma2, nr, nrm, ns1, ns2, nspino, nspinp

  complex(kind=db), dimension(n_comp*nlma,nlma2,nspinp,nspino,ip0:ip_max,ninitv):: rof

  logical:: Final_tddft, NRIXS, Radial_comp

  real(kind=db):: f_integr3, fac, radlr, radli, Rmtsd
  real(kind=db), dimension(nbseuil):: Vecond
  real(kind=db), dimension(nr):: r, fct
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,ip0:ip_max):: r_or_bess
  real(kind=db), dimension(nr,nbseuil):: psiir
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur

  if( Final_tddft ) then
    ns1 = initlv
    ns2 = initlv
  else
    ns1 = 1
    ns2 = nbseuil
  endif

  do is = 1,nbseuil
    psiir(1:nr,is) =  psii(1:nr,is) * r(1:nr) 
  end do

  lm1 = L**2
  lm2 = (L + 1)**2

  do ip = ip0,ip_max
  
    do is = ns1,ns2

      if( Final_tddft ) then
        iss = nbseuil   ! It is the second edge which is the energy origin
      else
        iss = is
        iseuil = is
      endif

      if( NRIXS ) then
        fac = 1._db
      else
        select case(ip)
          case(0)
            fac = 0.5_db * alfa_sf
          case(1)
            fac = 1._db
          case(2)
            fac = 0.5_db * Vecond(iss)
          case(3)
            fac = - ( 1._db / 6 ) * Vecond(iss)**2
        end select
      endif

      do isp = 1,nspinp
        do n1 = 1,nlm1
          do isol = 1,nspino
            do n2 = 1,nlm2

              fct(1:nr) = psiir(1:nr,iseuil) * ur(1:nr,n1,n2,isp,isol) * r_or_bess(1:nr,ip)
              radlr = fac * f_integr3(r,fct,1,nr,Rmtsd)
              if( Radial_comp ) then
                fct(1:nr) = psiir(1:nr,iseuil) * ui(1:nr,n1,n2,isp,isol) * r_or_bess(1:nr,ip)
                radli = fac * f_integr3(r,fct,1,nr,Rmtsd)
              else
                radli = 0._db
              endif

              if( nlm1 /= 1 .and. nlm2 /= 1 ) then
                rof(lm1+n1,lm1+n2,isp,isol,ip,is) = cmplx( radlr, radli,db)
              elseif( nlm1 /= 1 .and. nlma2 /= 1 ) then
                rof(lm1+n1,lm1+n1,isp,isol,ip,is) = cmplx( radlr, radli,db)
              elseif( nlm1 /= 1 ) then
                rof(lm1+n1,n2,isp,isol,ip,is) = cmplx( radlr, radli,db)
              elseif( nlm1 == 1 .and. nlma2 /= 1 ) then
                do lm = lm1+1,lm2
                  rof(lm,lm,isp,isol,ip,is) = cmplx(radlr, radli,db)
                end do
              else
                do lm = lm1+1,lm2
                  rof(lm,n2,isp,isol,ip,is) = cmplx(radlr, radli,db)
                end do
              endif

            end do
          end do
        end do
      end do
    end do
  end do

  return
end

!**********************************************************************

! Calculation of the radial integral of the (Fermi) golden rule
! ui et ur are the wave functions times r.

subroutine radial_matrix_optic(ip_max,ip0,ne,nlm1g,nlm_fp,nr,nrmtsd,nspinp,nspino,r,Radial_comp,Rmtsd,roff_rr,ui,ur,Vecond)

  use declarations
  implicit none

  integer:: ip, ip_max, ip0, ipp, ir, isol, isoli, isolf, isp, ispi, ispf, nlm1g, nlm_fp, n1i, n1f, n2i, n2f, ne, &
    nr, nrmtsd, nspinp, nspino

  logical:: Radial_comp

  real(kind=db):: f_integr3, fac, Rmtsd, Vecond
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nrmtsd):: rp, fct
  real(kind=db), dimension(nr,nlm1g,nlm_fp,nspinp,nspino,ne):: ui, ur
  real(kind=db), dimension(nlm1g,nlm1g,nlm_fp,nlm_fp,nspinp**2,nspino**2,ip0:ip_max):: roff_ii, roff_ir, roff_ri, roff_rr

  do ip = ip0,ip_max

    ipp = ip + 1

    rp(1:nrmtsd) = r(1:nrmtsd)**ipp

    select case(ip)
      case(0)
        fac = 0.5_db * alfa_sf
      case(1)
        fac = 1._db
      case(2)
        fac = 0.5_db * Vecond
      case(3)
      fac = - ( 1._db / 6 ) * Vecond**2
    end select

    do ispi = 1,nspinp
      do n1i = 1,nlm1g
        do isoli = 1,nspino
          do n2i = 1,nlm_fp

            do ispf = 1,nspinp
              isp = ispf + ( ispi - 1 ) * nspinp
              do n1f = 1,nlm1g
                do isolf = 1,nspino
                  isol = isolf + ( isoli - 1 ) * nspino
                  do n2f = 1,nlm_fp

                    do ir = 1,nrmtsd
                      fct(ir) = ur(ir,n1i,n2i,ispi,isoli,1) * ur(ir,n1f,n2f,ispf,isolf,2) * rp(ir)
                    end do

                    roff_rr(n1i,n1f,n2i,n2f,isp,isol,ip) = fac * f_integr3(rp,fct,1,nrmtsd,Rmtsd)

                    if( Radial_comp ) then
                      do ir = 1,nrmtsd
                        fct(ir) = ur(ir,n1i,n2i,ispi,isoli,1) * ui(ir,n1f,n2f,ispf,isolf,2) * rp(ir)
                      end do

                      roff_ri(n1i,n1f,n2i,n2f,isp,isol,ip) = fac * f_integr3(rp,fct,1,nrmtsd,Rmtsd)

                      do ir = 1,nrmtsd
                        fct(ir) = ui(ir,n1i,n2i,ispi,isoli,1) * ur(ir,n1f,n2f,ispf,isolf,2) * rp(ir)
                      end do

                      roff_ir(n1i,n1f,n2i,n2f,isp,isol,ip) = fac * f_integr3(rp,fct,1,nrmtsd,Rmtsd)

                      do ir = 1,nrmtsd
                        fct(ir) = ui(ir,n1i,n2i,ispi,isoli,1) * ui(ir,n1f,n2f,ispf,isolf,2) * rp(ir)
                      end do

                      roff_ii(n1i,n1f,n2i,n2f,isp,isol,ip) = fac * f_integr3(rp,fct,1,nrmtsd,Rmtsd)

                    endif

                  end do
                end do
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

! Calculation of the singular term coming from the irregular solution for the absorption
! Called by Radial

Subroutine Cal_Solsing(Classic_irreg,Ecomp,Eimag,f2,Final_tddft,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d, &
         icheck,initlv,ip_max,ip0,iseuil,konde,LL,lmax,m_hubb,nlm, &
         nlm1,nlm2,nbseuil,ninitv,nlma,nr,NRIXS,nrm,nrmtsd,nspin,nspino,nspinp, &
         numat,psii,r,r_or_bess,Radial_comp,Rmtsd,Singul,Spinorbite,Tau,ui,ur,V,V_hubb,Vecond)

  use declarations
  implicit none

  integer:: icheck, initlv, ip_max, ip0, ip1, ip2, is, ise, iseuil, isol, isp1, isp2, iss, L, l_hubbard, LL, lmax, m_hubb, &
    m1, m2, ms1, ms2, mv, n0, n1, n2, nbseuil, ninitv, nlm, nlm1, nlm2, np, nlma, nr, nrm, nrmax, nrmtsd, ns1, ns2, nspin, &
    nspino, nspinp, numat, nv, nv1, nv2

  complex(kind=db):: integr_sing, Sing
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlma,nspinp,nlma,nspinp,ip0:ip_max,ip0:ip_max,ninitv):: Singul
  complex(kind=db), dimension(nrmtsd):: f_reg, f_irg
  complex(kind=db), dimension(nlm1,nspinp,nlm1,nspinp):: Tau
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  logical:: Classic_irreg, Ecomp, Final_tddft, Full_potential, Hubb_a, Hubb_m, Hubb_d, NRIXS, Radial_comp, Spinorbite

  real(kind=db):: Eimag, fac1, fac2, Rmtsd
  real(kind=db), dimension(nbseuil):: Vecond
  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nrmtsd):: rr, phi1, phi2
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nspin):: g0, gm, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr,nlm,nlm,nspin):: V
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur
  real(kind=db), dimension(nr,ip0:ip_max):: r_or_bess

  real(kind=db), dimension(:,:,:,:,:), allocatable:: usi, usr

  rr(1:nrmtsd) = r(1:nrmtsd)

  if( Final_tddft ) then
    ns1 = initlv
    ns2 = initlv
  else
    ns1 = 1
    ns2 = nbseuil
  endif

  nrmax = nrmtsd+1 
  allocate( usi(nrmax,nlm1,nlm2,nspinp,nspino) )
  allocate( usr(nrmax,nlm1,nlm2,nspinp,nspino) )

! gm et gp inverse in subroutine
  call Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde, &
       LL,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtsd,nrmax,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Spinorbite,Tau,usi,usr,V,V_hubb)

  do is = ns1,ns2

    if( Final_tddft ) then
      iss = nbseuil   ! it is the second edge which is the origin (because at lower energy)
      ise = iseuil
    else
      iss = is
      ise = is
    endif

    do ip1 = ip0,ip_max  ! boucle sur dipole, quadrupole, Octupole

      if( NRIXS ) then
        fac1 = 1._db
      else
        select case(ip1)
          case(0)
            fac1 = - 0.5_db * alfa_sf
          case(1)
            fac1 = 1._db
          case(2)
            fac1 = 0.5_db * Vecond(iss)
          case(3)
            fac1 = - ( 1._db / 6 ) * Vecond(iss)**2
        end select
      endif

      phi1(1:nrmtsd) = psii(1:nrmtsd,ise) * r_or_bess(1:nrmtsd,ip1)

      do ip2 = ip0,ip_max  ! loop over multipole

        if( NRIXS ) then
          fac2 = 1._db
        else
          select case(ip2)
            case(0)
              fac2 = - 0.5_db * alfa_sf
            case(1)
              fac2 = 1._db
            case(2)
              fac2 = 0.5_db * Vecond(iss)
            case(3)
              fac2 = - ( 1._db / 6 ) * Vecond(iss)**2
          end select
        endif

        phi2(1:nrmtsd) = psii(1:nrmtsd,ise) * r_or_bess(1:nrmtsd,ip2)

        n1 = 0
        do L = 0,lmax
          if( .not. Full_potential .and. L /= LL ) cycle
          Hubb_m = Hubb_a .and. L == l_hubbard( numat )

          n0 = L**2 + L + 1
          do m1 = -L,L
            if( nlm1 == 1 .and. m1 /= 0 ) cycle
            n1 = n1 + 1
            do m2 = -L,L
              if( nlm1 == 1 .and. m2 /= 0 ) cycle
              if( .not. ( Full_potential .or. ( Hubb_m .and. .not. Hubb_d ) ) .and. abs(m2-m1) > 1 ) cycle
              n2 = n1 + m2 - m1
        
              do isp1 = 1,nspinp
                do isp2 = 1,nspinp
                  if( .not. spinorbite .and. isp1 /= isp2 ) cycle
                  if( nlm1 == 1 .and. isp1 /= isp2 ) cycle
                  if( .not. ( Full_potential .or. ( Hubb_m .and. .not. Hubb_d ) ) .and. isp2 - isp1 /= m2 - m1 ) cycle
            
                  do np = 1,nlm2
            
! Index solution is diagonal
                    do isol = 1,nspino
                      if( Spinorbite .and. nlm2 == 1 ) then
                        ms1 = m1 + isol - isp1
                        if( ms1 > L .or. ms1 < -L ) cycle
                        ms2 = m2 + isol - isp2
                        if( ms2 > L .or. ms2 < -L ) cycle
                      else
                        ms1 = m1
                        ms2 = m2
                      endif
            
                      if( Radial_comp ) then
                        f_reg(1:nrmtsd) = r(1:nrmtsd) * cmplx( ur(1:nrmtsd,n1,np,isp1,isol), ui(1:nrmtsd,n1,np,isp1,isol), db )
                      else
                        f_reg(1:nrmtsd) = r(1:nrmtsd) * cmplx( ur(1:nrmtsd,n1,np,isp1,isol), 0._db, db )
                      endif
            
                      f_irg(1:nrmtsd) = cmplx( usr(1:nrmtsd,n2,np,isp2,isol), usi(1:nrmtsd,n2,np,isp2,isol), db)
            
                      Sing = fac1 * fac2 * integr_sing(nrmtsd,phi1,phi2,f_reg,f_irg,Rmtsd,rr,icheck)
            
                      if( nlm1 == 1 ) then
                        do mv = -L,L
                          nv = n0 + mv
                          Singul(nv,isp1,nv,isp2,ip1,ip2,is) = Singul(nv,isp1,nv,isp2,ip1,ip2,is) + Sing
                        end do
                      else
                        nv1 = n0 + m1
                        nv2 = n0 + m2
                        Singul(nv1,isp1,nv2,isp2,ip1,ip2,is) = Singul(nv1,isp1,nv2,isp2,ip1,ip2,is) + Sing
                      endif
            
                    end do
                  end do
            
                end do    ! end of loop isp2
              end do    ! end of loop isp1
            end do   ! end of loop m2
          end do   ! end of loop m1
        end do  ! end of loop L

      end do ! end of loop dipole, quadrupole
    end do ! end of loop dipole, quadrupole

  end do ! end of loop seuil

  deallocate( usi, usr )
  
  return
end

!***********************************************************************

! Resolution of the radial Schrodinger equation for the irregular solution 

! Copy of Sch_radial: gm et gp are reversed in the call

subroutine Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gp,gm,gso,Hubb_a,Hubb_d,icheck,konde, &
       LL,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtsd,nrmax,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Spinorbite,Tau,usi,usr,V,V_hubb)

  use declarations
  implicit none

  integer:: icheck, im, ip, ir, isol, isp, isq, isr, L, l_hubbard, L2, li, LL, lf, lmax, lp, m, m_hubb, mp, ms, n, &
    nlm, nlm1, nlm2, np, nr, nrmax, nrmtsd, ns, nspin, nspino, nspinp, numat

  logical:: Classic_irreg, Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m, Hubb_nd, Radial_comp, Spinorbite

  complex(kind=db):: fnormc, V_h, z
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(0:lmax):: bess, neum
  complex(kind=db), dimension(nlm1,nspinp,nlm1,nspinp):: Tau
  complex(kind=db), dimension(nrmtsd:nrmax,0:lmax,nspin):: Bessel, Hankel
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: us

  real(kind=db):: Eimag, fac, fnorm, konder, Rmtsd, td, zr
  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(0:lmax):: bessr, neumr
  real(kind=db), dimension(nr,nspin):: g0, gm, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr,nlm,nlm,nspin):: V
  real(kind=db), dimension(:), allocatable:: rr
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur

  real(kind=db), dimension(nrmax,nlm1,nlm2,nspinp,nspino):: usi, usr
  
  allocate( us(nrmax,nlm1,nlm2,nspinp,nspino)  )

  us(:,:,:,:,:) = ( 0._db, 0._db )

! Hankel is in fact -i.hankel+
    do ir = nrmtsd,nrmax
      do isp = 1,nspin
        if( Ecomp ) then
          fnormc = sqrt( konde(isp) / pi )
          z = konde(isp) * r(ir)
          call cbessneu(fnormc,z,lmax,lmax,bess,neum)
          if( Classic_irreg ) then
            Hankel(ir,0:lmax,isp) = neum(0:lmax) - img * bess(0:lmax)
          else
 !  The minus is to change the sign of the irregular solution corresponding to the atomic solution
            Bessel(ir,0:lmax,isp) = - bess(0:lmax)
          endif
        else
          konder = Real( konde(isp), db )
          fnorm = sqrt( konder / pi )
          zr = konder * r(ir)
          call cbessneur(fnorm,zr,lmax,lmax,bessr,neumr)
          if( Classic_irreg ) then
            Hankel(ir,0:lmax,isp) = neumr(0:lmax) - img * bessr(0:lmax)
          else
            Bessel(ir,0:lmax,isp) = - cmplx( bessr(0:lmax), 0._db, db )
          endif
        endif
      end do
    end do

  if( Full_potential ) then
    li = 0
    lf = lmax
  else
    li = LL
    lf = LL
  endif

  n = 0
  do L = li,lf

    do m = -L,L
      if( nlm1 == 1 .and. m /= 0 ) cycle
      n = n + 1

      do isp = 1,nspinp  ! Spin reel
        isr = min( isp, nspin )
        
        np = 0
        do lp = li,lf
          if( .not. Full_potential .and. L /= lp ) cycle
          do mp = -lp,lp
            if( nlm2 == 1 .and. mp /= m ) cycle
            np = np  + 1

            do isq = 1,nspinp  ! spin d'attaque ( indice solution )
              isol = min( isq, nspino )
              if( nspino == 1 .and. isp /= isq ) cycle
              if( nlm2 == 1 ) then
                ms = m + isq - isp
                if( ms > L .or. ms < -L ) cycle
                ns = n + isq - isp
              else
                ns = np
              endif

              do ir = nrmtsd,nrmax
                if( Classic_irreg ) then
                  us(ir,n,np,isp,isol) = Tau(ns,isq,n,isp) * r(ir) * Hankel(ir,L,isr)
                else
                  us(ir,n,np,isp,isol) = r(ir) * Bessel(ir,L,isr)
                endif
              end do
            end do
          end do
        end do
      end do

    end do
  end do

  do ir = nrmtsd,2,-1
    ip = ir - 1      ! inverse par rapport a Sch_radial
    im = ir + 1

    do isp = 1,nspinp
      isr = min( isp, nspin )
      isq = 3 - isp

      do L = li,lf

        L2 = L * ( L + 1 )

        if( Hubb_a .and. L == l_hubbard(numat) ) then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif
        Hubb_nd = Hubb_m .and. .not. Hubb_d

        do m = -L,L

          if( Full_potential ) then
            n = L**2 + L + 1 + m
          elseif( Hubb_m .or. Spinorbite ) then
            n = L + 1 + m
          else
            n = 1
            if( m /= 0 ) cycle
          endif

          td = g0(ir,isr) + L2 * f2(ir)

          if( Spinorbite ) then
            if( isp == 1 ) then
              ms = m
              mp = m + 1   ! m de L'autre spin
            else
              ms = - m
              mp = m - 1
            endif
            fac = sqrt( ( L - ms ) * ( L + ms + 1._db ) )
            td = td + ms * gso(ir,isp)
          else
            mp = L + 1 ! pour que le if en dessous soit faux
          endif
          if( Full_potential ) then
            np = L**2 + L + 1 + mp
          elseif( Hubb_nd .or. Spinorbite ) then
            np = L + 1 + mp
          else
            np = 1
          endif

! the second index "m" is the non zero component at origin
          us(ip,n,:,isp,:) = td * us(ir,n,:,isp,:) + gm(ir,isr) * us(im,n,:,isp,:)
          if( Ecomp ) us(ip,n,:,isp,:) = us(ip,n,:,isp,:) - img * Eimag * us(ir,n,:,isp,:)

          if( abs(mp) <= L ) us(ip,n,:,isp,:) = us(ip,n,:,isp,:) + fac * gso(ir,isr) * us(ir,np,:,isq,:)

          if( Hubb_m ) then
            do mp = -L,L
              if( Hubb_d .and. m /= mp ) cycle
              if( Full_potential ) then
                np = L**2 + L + 1 + mp
              else
                np = L + 1 + mp
              endif
              do isol = 1,nspino
                if( Spinorbite ) then
                  V_h = V_hubb(m,mp,isol,isol) 
                else
                  V_h = V_hubb(m,mp,isp,isp) 
                endif
                us(ip,n,:,isp,isol) = us(ip,n,:,isp,isol) + V_h * us(ir,np,:,isp,isol)
              end do
            end do
          endif

          if( Full_potential ) then
            do np = 1,nlm
              if( n == np ) cycle
              us(ip,n,:,isp,:) = us(ip,n,:,isp,:) + V(ir,n,np,isr) * us(ir,np,:,isp,:)
            end do
          endif

          us(ip,n,:,isp,:) = - us(ip,n,:,isp,:) * gp(ir,isr)

        end do
      end do
    end do
  end do

  if( icheck > 2 ) then
    n = nrmax
    allocate( ui(n,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(n,nlm1,nlm2,nspinp,nspino) )
    allocate( rr(n) )
    ur(1:n,:,:,:,:) = Real( us(1:n,:,:,:,:), db )
    ui(1:n,:,:,:,:) = aimag( us(1:n,:,:,:,:) )
    rr(1:n) = r(1:n)
    call write_ur(Full_potential,li,lf,nlm1,nlm2,n,nspinp,nspino,numat,rr,Radial_comp,Rmtsd,ui,ur,3)
    deallocate( rr, ui, ur )
  endif

  usr(:,:,:,:,:) = real( us(:,:,:,:,:), db )
  usi(:,:,:,:,:) = aimag( us(:,:,:,:,:) )
  
  deallocate( us )
  
  return
end

!**********************************************************************

! Calcul des integrales radiales dans le cas optique
! Appele par tenseur_car

subroutine radial_optic(Classic_irreg,Ecinetic_t,Eimag_t,Energ,Enervide,Full_potential,Hubb_a,Hubb_d,icheck,ief,ip_max,ip0, &
         lmax,lmax_pot,m_depend,m_hubb,n_ph,n_Ec,n_V,nlm_pot,nlm1g,nlm_fp,No_diag,nr,nr_zet,nspin,nspino,nspinp,numat,r, &
         Relativiste,Renorm,Rmtg,Rmtsd,roff_rr,roff0,Solsing,Spinorbite,V_hubb,V_intmax,V0bd_t,Vrato_t,Ylm_comp,zet)

  use declarations
  implicit none

  integer:: i, icheck, ie, ie1, ie2, ief, ip, ip_max, ip0, iso1, iso2, isp, isp1, isp2, j, L, lf, l_hubbard, lfin, lm, lmf, &
    lm1, lm2, lmax, lmax_pot, lmp, lmpf, lp, lpf, m, m_hubb, mf, mp, mpf, n_ph, nlm, nlm_pot, nlm1, nlm1g, nlm2, nlm_fp, &
    n_Ec, n_V, nr, nr_zet, nrmax, nrmtsd, nrmtg, nspin, nspino, nspinp, numat

  character(len=104):: mot

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: us
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: us_t

  logical:: Classic_irreg, Ecomp, Diag, Full_potential, Hubb_a, Hubb_d, Hubb_m, m_depend, No_diag, Radial_comp, Relativiste, &
    Renorm, Solsing, Spinorbite, Ylm_comp

  real(kind=db):: Eimag, Rmtg, Rmtsd, V_intmax, Vecond
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato
  real(kind=db), dimension(n_Ec):: Eimag_t, Enervide
  real(kind=db), dimension(nspin,n_Ec):: Ecinetic_t
  real(kind=db), dimension(n_ph):: Energ
  real(kind=db), dimension(nspin,n_V):: V0bd_t
  real(kind=db), dimension(nr,nlm_pot,nspin,n_V):: Vrato_t

  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nr,nspin):: g0, gm, gmi, gpi, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr_zet,nlm1g,nlm_fp,nspinp,nspino,n_Ec):: zet
  real(kind=db), dimension(nlm1g,nlm1g,nlm_fp,nlm_fp,nspinp**2,nspino**2,ip0:ip_max):: roff_rr
  real(kind=db), dimension(ief-1,ief:n_Ec,nlm1g,nlm1g,nspinp,nspinp,nspino**2):: roff0

  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur, usi, usr
  real(kind=db), dimension(:,:,:,:,:,:), allocatable:: ui_t, ui1, ur_t, ur1

  roff_rr(:,:,:,:,:,:,:) = 0._db

! Terme multiplicatif pour les transitions quadrupolaires
! En S.I. vecond = k = E*alfa_sf*4*pi*epsilon0 / (e*e)
! En ua et rydb : k = 0.5_db * alfa_sf * E
  Vecond = 0.5_db * alfa_sf * Energ(1)

  if( Full_potential ) then
    nlm = ( lmax + 1 )**2
  else
    nlm = 1
  endif
  allocate( V(nr,nlm,nlm,nspin) )

  Vrato(:,:,:) = Vrato_t(:,:,:,1)
  V0bd(:) = V0bd_t(:,1)

  call mod_V(icheck,lmax,lmax_pot,nlm,nlm_pot,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  Ecomp = .false.
  if( abs(Eimag_t(1)) > eps10 ) then
    Ecomp = .true.
  else
    do isp = 1,nspin
      do ie = 1,n_Ec
        if( Ecinetic_t(isp,ie) > eps10 ) cycle
        Ecomp = .true.
        exit
      end do
    end do
  endif
  Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

  if( Full_potential ) then
    lfin = 0
  else
    lfin = lmax
  endif

  allocate( ur_t(nr,nlm1g,nlm_fp,nspinp,nspino,n_Ec) )
  allocate( ui_t(nr,nlm1g,nlm_fp,nspinp,nspino,n_Ec) )

  do L = 0,lfin

    if( Hubb_a .and. L == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = ( lmax + 1 )**2
      nlm2 = nlm1
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*L + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*L + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    nrmax = nrmtsd + 1 
    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )
    allocate( us(nrmax,nlm1,nlm2,nspinp,nspino) )
    allocate( usi(nrmax,nlm1,nlm2,nspinp,nspino) )
    allocate( usr(nrmax,nlm1,nlm2,nspinp,nspino) )

    if( L == 0 ) allocate( us_t(nrmax,nlm1g,nlm_fp,nspinp,nspino,n_Ec) )

    do ie = 1,n_Ec

      if( icheck > 1 ) write(3,110) ie, Enervide(ie)*rydb, Ecinetic_t(:,ie)*rydb 
      if( icheck > 1 .and. n_V <= 2 ) write(3,120) V0bd_t(1,:)*rydb 
     
      Ecinetic(:) = Ecinetic_t(:,ie)
      Eimag = Eimag_t(ie)
      konde(:) = sqrt( cmplx(Ecinetic(:), Eimag, db) )

      call coef_sch_rad(Enervide(ie),f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)

      gpi(:,:) = 1 / gp(:,:)

      call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,L,lmax,m_hubb, &
        nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
        ur,V,V_hubb)

! Calcul de la solution singuliere
! gm et gp inverse dans le sousprogramme
      if( Solsing ) then
        gmi(:,:) =  1 / gm(:,:)
        call Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gmi,gp,gso,Hubb_a,Hubb_d,icheck,konde, &
           L,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtsd,nrmax,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Spinorbite,Tau,usi,usr,V, &
           V_hubb)
      endif

      us(:,:,:,:,:) = cmplx( usr(:,:,:,:,:), usi(:,:,:,:,:), db )
      
      if( Radial_comp ) then
        if( Full_potential ) then
          ur_t(:,:,:,:,:,ie) = ur(:,:,:,:,:)
          ui_t(:,:,:,:,:,ie) = ui(:,:,:,:,:)
          us_t(:,:,:,:,:,ie) = us(:,:,:,:,:)
        elseif( Hubb_m .and. .not. Hubb_d ) then
          ur_t(:,L**2+1:L**2+nlm1,L**2+1:L**2+nlm1,:,:,ie) = ur(:,1:nlm1,1:nlm1,:,:)
          ui_t(:,L**2+1:L**2+nlm1,L**2+1:L**2+nlm1,:,:,ie) = ui(:,1:nlm1,1:nlm1,:,:)
          us_t(:,L**2+1:L**2+nlm1,L**2+1:L**2+nlm1,:,:,ie) = us(:,1:nlm1,1:nlm1,:,:)
        elseif( No_diag ) then
          do lm1 = L**2+1,(L+1)**2
            lm = min( lm1 - L**2, nlm1 )
            ur_t(:,lm1,lm1,:,:,ie) = ur(:,lm,1,:,:)
            ui_t(:,lm1,lm1,:,:,ie) = ui(:,lm,1,:,:)
            us_t(:,lm1,lm1,:,:,ie) = us(:,lm,1,:,:)
          end do
        elseif( Spinorbite .or. Hubb_m ) then
           do lm1 = L**2+1,(L+1)**2
             lm = lm1 - L**2
             ur_t(:,lm1,1,:,:,ie) = ur(:,lm,1,:,:)
             ui_t(:,lm1,1,:,:,ie) = ui(:,lm,1,:,:)
             us_t(:,lm1,1,:,:,ie) = us(:,lm,1,:,:)
           end do
        elseif( m_depend ) then
          do lm1 = L**2+1,(L+1)**2
            ur_t(:,lm1,1,:,:,ie) = ur(:,1,1,:,:)
            ui_t(:,lm1,1,:,:,ie) = ui(:,1,1,:,:)
            us_t(:,lm1,1,:,:,ie) = us(:,1,1,:,:)
          end do
        else
          ur_t(:,L,1,:,:,ie) = ur(:,1,1,:,:)
          ui_t(:,L,1,:,:,ie) = ui(:,1,1,:,:)
          us_t(:,L,1,:,:,ie) = us(:,1,1,:,:)
        endif
      else
        if( Full_potential ) then
          ur_t(:,:,:,:,:,ie) = ur(:,:,:,:,:)
        elseif( Hubb_m .and. .not. Hubb_d ) then
          ur_t(:,L**2+1:L**2+nlm1,L**2+1:L**2+nlm1,:,:,ie) = ur(:,1:nlm1,1:nlm1,:,:)
        elseif( No_diag ) then
          do lm1 = L**2+1,(L+1)**2
            lm = min( lm1 - L**2, nlm1 )
            ur_t(:,lm1,lm1,:,:,ie) = ur(:,lm,1,:,:)
          end do
        elseif( Spinorbite .or. Hubb_m ) then
           do lm1 = L**2+1,(L+1)**2
             lm = lm1 - L**2
             ur_t(:,lm1,1,:,:,ie) = ur(:,lm,1,:,:)
           end do
        elseif( m_depend ) then
          do lm1 = L**2+1,(L+1)**2
            ur_t(:,lm1,1,:,:,ie) = ur(:,1,1,:,:)
          end do
        else
          ur_t(:,L+1,1,:,:,ie) = ur(:,1,1,:,:)
        endif
      endif

    end do

    deallocate( Tau, ui, ur, us, usi, usr )

  end do   ! fin de la boucle sur L

  if( nr_zet > 0 ) zet(:,:,:,:,:,:) = ur_t(:,:,:,:,:,:) 

  allocate( ur1(nr,nlm1g,nlm_fp,nspinp,nspino,2) )
  allocate( ui1(nr,nlm1g,nlm_fp,nspinp,nspino,2) )

  do ie1 = 1,ief-1

    ur1(:,:,:,:,:,1) = ur_t(:,:,:,:,:,ie1)
    ui1(:,:,:,:,:,1) = ui_t(:,:,:,:,:,ie1)

    do ie2 = ief, n_Ec
    
      ur1(:,:,:,:,:,2) = ur_t(:,:,:,:,:,ie2)
      ui1(:,:,:,:,:,2) = ui_t(:,:,:,:,:,ie2)
     
      call radial_matrix_optic(ip_max,ip0,n_Ec,nlm1g,nlm_fp,nr,nrmtsd,nspinp,nspino,r,Radial_comp,Rmtsd,roff_rr,ui1,ur1,Vecond)

      if( ief > 1 .and. ip0 == 0 ) then
        i = 0
        do isp1 = 1,nspinp
          do isp2 = 1,nspinp
            i = i + 1
            j = 0
            do iso1 = 1,nspino
              do iso2 = 1,nspino
                j = j + 1
                if( nlm_fp == 1 ) then
                  roff0(ie1,ie2,:,:,isp1,isp2,j) = roff_rr(:,:,1,1,i,j,0)
                else
                  roff0(ie1,ie2,:,:,isp1,isp2,j) = ( 0._db, 0._db )
                  do lm1 = 1,nlm_fp
                    do lm2 = 1,nlm_fp
                      roff0(ie1,ie2,lm1,lm2,isp1,isp2,j) = sum( roff_rr(:,:,lm1,lm2,i,j,0) )
                    end do
                  end do 
                endif
              end do
            end do
          end do
        end do

      endif

    end do
  end do

  deallocate( ui_t, ui1, ur_t, ur1, us_t, V )

  if( icheck > 1 ) then

    Diag = .true.
    
    if( Diag ) then
      write(3,130) (((( lm, isp, iso1, iso2, iso2 = 1,nspino), iso1 = 1,nspino), isp = 1,nspinp), lm = 1,nlm1g)
      do ie1 = 1,ief-1
        do ie2 = ief, n_Ec
          write(3,140) Enervide(ie1)*rydb, Enervide(ie2)*rydb, (( roff0(ie1,ie2,lm,lm,isp,isp,:), isp = 1,nspinp), &
                                                                            lm = 1,nlm1g )
        end do
      end do
    else
      write(3,145) (((((( lm1, lm2, isp1, isp2, iso1, iso2, iso2 = 1,nspino), iso1 = 1,nspino), isp2 = 1,nspinp), &
                                                            isp1 = 1,nspinp), lm2 = 1,nlm_fp), lm1 = 1,nlm1g )
      do ie1 = 1,ief-1
        do ie2 = ief, n_Ec
          write(3,140) Enervide(ie1)*rydb, Enervide(ie2)*rydb, (((( roff0(ie1,ie2,lm1,lm2,isp1,isp2,:), &
                                          isp2 = 1,nspinp), isp1 = 1,nspinp), lm2 = 1,nlm_fp), lm1 = 1,nlm1g )
        end do
      end do
    endif
    
    mot = ' '
    if( Ecomp .and. nspino == 2 ) then
      mot(11:17) = 'up sol1'
      mot(37:43) = 'up sol2'
      mot(63:69) = 'dn sol1'
      mot(89:95) = 'dn sol2'
    elseif( nspino == 2 ) then
      mot(5:50) = 'up sol1      up sol2      dn sol1      dn sol2'
    elseif( nspinp == 2 ) then
      mot = '      up           dn'
    else
      mot = ' '
    endif

    do ip = ip0,ip_max
      select case(ip)
        case(0)
          write(3,160) '  monopole'
        case(1)
          write(3,160) '    dipole'
        case(2)
          write(3,160) 'quadrupole'
        case(3)
          write(3,160) '  octupole'
      end select
      if( Full_potential .or. Hubb_a ) then
        write(3,170) ' L1 m1 lp mp L2 m2 lp mp', mot
      else
        write(3,175) ' L1 m1 L2 m2', mot
      endif

      do L = 0,lmax
        do m = -L,L
          if( m_depend ) then
            lm = L**2 + L + 1 + m
          else
            if( m /= 0 ) cycle
            lm = L + 1
          endif
          do lp = 0,lmax
            if( .not. Full_potential .and. L /= lp ) cycle
            do mp = -lp,lp
              if( .not. ( Hubb_a .or. Full_potential .or. m == mp ) ) cycle
              lmp = min( lp**2 + lp + 1 + mp, nlm_fp )

              do lf = 0,lmax
                do mf = -lf,lf
                  if( m_depend ) then
                    lmf = lf**2 + lf + 1 + mf
                  else
                    if( mf /= 0 ) cycle
                    lmf = lf + 1
                  endif
                  lmpf = 0
                  do lpf = 0,lmax
                    if( .not. Full_potential .and. lf /= lpf ) cycle
                    do mpf = -lpf,lpf
                      if( .not. ( Hubb_a .or. Full_potential .or. mf == mpf ) ) cycle
                      lmpf = lpf**2 + lpf + 1 + mpf
                      lmpf = min( nlm_fp, lmpf )

                      if( Full_potential .or. Hubb_a ) then
                        write(3,180) L, m, lp, mp, lf, mf, lpf, mpf, (roff_rr(lm,lmf,lmp,lmpf,isp,:,ip), isp = 1,nspinp**2)
                      else
                        write(3,185) L, m, lf, mf, (roff_rr(lm,lmf,lmp,lmpf,isp,:,ip), isp = 1,nspinp**2)
                      endif
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end if

  return
  110 format(/' Energy index =',i4,',  Enervide =',f10.3,',  Ecinetic =',2f10.3,' eV')
  120 format(22x,'V0bd     =',4f10.3)
  130 format(/' Roff0',/'   Enervide_1   Enervide_2', 1000(4x,i3,3i2) ) 
  140 format(2f13.5,1p,1000e13.5)
  145 format(/' Roff0',/'   Enervide_1   Enervide_2', 1000(3x,2i3,4i1) ) 
  160 format(/' Radial matrix (Roff_rr), ',a10,1x,'term:')
  170 format(a24,a104)
  175 format(a12,a104)
  180 format(8i3,1p,16e13.5)
  185 format(4i3,1p,16e13.5)
end

!***********************************************************************

! Calculation of basis functions on the FDM grid points 

subroutine cal_phiato(Full_potential,ia,iang,ibord,icheck,iopsymr,LL,lmax,n_comp,nlm1,nlm2,natome, &
             nbord,nbtm,nbtm_fdm,nlmagm,nlmmax,nphiato1,nphiato7,npsom,nr,nspinp,nspino,phiato,posi,r, &
             Radial_comp,Spinorbite,ui,ur,usi,usr,xyz,Ylm_comp,Ylmato)

  use declarations
  implicit none

  integer:: i, i_comp, ia, iang, ib, icheck, isol, isp, isym, jsol, jsp, L, LL, lm, lm0, lmax, lmv, Lmw, &
    lp, m, mp, mpp, mv, n, n_comp, natome, nbord, nbtm, nbtm_fdm, nlm1, nlm2, nlmagm, nlmmax, np, &
    nphiato1, nphiato7, npsom, nr, nspinp, nspino
  integer, dimension(nopsm):: iopsymr
  integer, dimension(nbtm,natome):: ibord

  complex(kind=db):: phic, phid, yc

  logical:: Full_potential, Radial_comp, Spinorbite, Ylm_comp

  real(kind=db):: f_interp2, phi, phii, r0, rm, rp, rrel, u0, uc0, ucm, ucp, um, up
  real(kind=db), dimension(3):: ps, x, w
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: usi, usr
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(nbtm_fdm,nlmmax,natome):: Ylmato
  real(kind=db), dimension(nphiato1,n_comp*nlmagm,nspinp,nspino,natome,nphiato7):: phiato

  ps(1:3) = posi(1:3,ia)

  do ib = 1,nbord

    i = ibord(ib,ia)
    x(1:3) = xyz(1:3,i)

    call posrel(iopsymr,x,ps,w,rrel,isym)

    do i = 2,nr
      if( r(i) > rrel ) exit
    end do
    i = min(i,nr)   
    i = max(i,3)

    rm = r(i-2)
    r0 = r(i-1)
    rp = r(i)

    n = 0
    do L = 0,lmax
      if( .not. Full_potential .and. L /= LL ) cycle
      lm0 = L**2 + L + 1
      do m = -L,L
        if( nlm1 == 1 .and. m /= 0 ) cycle
        n = n + 1
        do isp = 1,nspinp
          if( iang == 1 ) then
            jsp = isp
          elseif( iang == - 1 ) then
            jsp = nspinp + 1 - isp
          endif

          np = 0
          do lp = 0,lmax
            if( .not. Full_potential .and. lp /= L ) cycle
            do mp = -lp,lp
              if( nlm2 == 1 .and. mp /= m ) cycle
              np = np + 1
              Loop_isol: do isol = 1,nspino
                if( iang == 1 ) then
                  jsol = isol
                elseif( iang == - 1 ) then
                  jsol = nspino + 1 - isol
                endif

                if( nlm2 > 1 ) then
                  lmv = lp**2 + lp + 1 + mp
                elseif( Spinorbite ) then
                  mv = m + isol - isp  ! m  faux
                  if( mv > L .or. mv < -L ) cycle
                  lmv = lm0 + m
                else
                  lmv = lm0 + m
                endif
                mv = m

                do i_comp = 1,n_comp

                  if( i_comp == 1 ) then
                                   
                    if( iang == 0 ) then
                      if( nspino == 2 .and. isp /= isol ) cycle
                      um = 0.5_db * ( ur(i-2,n,np,1,1) + ur(i-2,n,np,min(2,nspinp),min(2,nspino)) )
                      u0 = 0.5_db * ( ur(i-1,n,np,1,1) + ur(i-1,n,np,min(2,nspinp),min(2,nspino)) )
                      up = 0.5_db * ( ur(i,n,np,1,1) + ur(i,n,np,min(2,nspinp),min(2,nspino)) )
                    else
                      um = ur(i-2,n,np,jsp,jsol)
                      u0 = ur(i-1,n,np,jsp,jsol)
                      up = ur(i,n,np,jsp,jsol)
                    endif
               
                    if( Radial_comp ) then
                      if( iang == 0 ) then
                        if( nspino == 2 .and. isp /= isol ) cycle
                        ucm = 0.5_db * ( ui(i-2,n,np,1,1) + ui(i-2,n,np,min(2,nspinp),min(2,nspino)) )
                        uc0 = 0.5_db * ( ui(i-1,n,np,1,1) + ui(i-1,n,np,min(2,nspinp),min(2,nspino)) )
                        ucp = 0.5_db * ( ui(i,n,np,1,1) + ui(i,n,np,min(2,nspinp),min(2,nspino)) )
                      else
                        ucm = ui(i-2,n,np,jsp,jsol)
                        uc0 = ui(i-1,n,np,jsp,jsol)
                        ucp = ui(i,n,np,jsp,jsol)
                      endif
                    endif
                    
                  else
                                   
                    if( iang == 0 ) then
                      if( nspino == 2 .and. isp /= isol ) cycle
                      um = 0.5_db * ( usr(i-2,n,np,1,1) + usr(i-2,n,np,min(2,nspinp),min(2,nspino)) )
                      u0 = 0.5_db * ( usr(i-1,n,np,1,1) + usr(i-1,n,np,min(2,nspinp),min(2,nspino)) )
                      up = 0.5_db * ( usr(i,n,np,1,1) + usr(i,n,np,min(2,nspinp),min(2,nspino)) )
                    else
                      um = usr(i-2,n,np,jsp,jsol)
                      u0 = usr(i-1,n,np,jsp,jsol)
                      up = usr(i,n,np,jsp,jsol)
                    endif
               
                    if( Radial_comp ) then
                      if( iang == 0 ) then
                        if( nspino == 2 .and. isp /= isol ) cycle Loop_isol
                        ucm = 0.5_db * ( usi(i-2,n,np,1,1) + usi(i-2,n,np,min(2,nspinp),min(2,nspino)) )
                        uc0 = 0.5_db * ( usi(i-1,n,np,1,1) + usi(i-1,n,np,min(2,nspinp),min(2,nspino)) )
                        ucp = 0.5_db * ( usi(i,n,np,1,1) + usi(i,n,np,min(2,nspinp),min(2,nspino)) )
                      else
                        ucm = usi(i-2,n,np,jsp,jsol)
                        uc0 = usi(i-1,n,np,jsp,jsol)
                        ucp = usi(i,n,np,jsp,jsol)
                      endif
                    endif
                    
                  endif
             
                  phi = f_interp2(rrel,rm,r0,rp,um,u0,up)
                  if( Radial_comp ) then
                    phii = f_interp2(rrel,rm,r0,rp,ucm,uc0,ucp)
                  else
                    phii = 0._db
                  endif
                  phic = cmplx( phi, phii, db )
             
                  do mpp = -L,L
                    if( nlm1 /= 1 .and. mpp /= m ) cycle
                    lm = lm0 + mpp
                    if( nlm1 == 1 ) then
                      lmv = lm
                      mv = mpp
                    endif
                    if( i_comp == 1 ) then
                      lmw = lmv
                    else
                      lmw = lmv + nlmagm
                    endif
                    if( mpp == 0 .or. .not. Ylm_comp ) then
                      phiato(ib,lmw,isp,isol,ia,1) = phiato(ib,lmw,isp,isol,ia,1) + phi * Ylmato(ib,lm,ia)
                      if( Radial_comp ) phiato(ib,lmw,isp,isol,ia,2) = phiato(ib,lmw,isp,isol,ia,2) + phii * Ylmato(ib,lm,ia)
                    else
                      phid = phic * yc(mv,Ylmato(ib,lm,ia),Ylmato(ib,lm-2*mv,ia))
                      phiato(ib,lmw,isp,isol,ia,1) = phiato(ib,lmw,isp,isol,ia,1) + Real( phid, db)
                      phiato(ib,lmw,isp,isol,ia,2) = phiato(ib,lmw,isp,isol,ia,2) +  aimag( phid )
                    endif

                  end do             
 
                end do
              end do Loop_isol
            end do
          end do

        end do
      end do
    end do
  end do

  if( icheck > 2 .and. ( LL == lmax .or. Full_potential ) ) then
    do i_comp = 1,n_comp
      if( i_comp == 1 ) then
        write(3,110) ia
      else
        write(3,115) ia
      endif
      do L = 0,lmax
        lm0 = L**2 + L + 1
        do mp = -L,L
          if( nphiato7 == 1 ) then
            if( Spinorbite ) then
              write(3,120) L, mp
            else
              write(3,130) L, mp
            endif
          else
            if( Spinorbite ) then
              write(3,140) L, mp
            else
              write(3,150) L, mp
            endif
          endif
          lm = lm0 + mp
          if( i_comp == 2 ) lm = lm + nlmagm
          do ib = 1,nbord
            write(3,160) ibord(ib,ia), ( ( phiato(ib,lm,isp,isol,ia,:), isp = 1,nspinp ), isol = 1,nspino )
          end do
        end do
      end do
    end do
  endif

  return
  110 format(/' ia =',i3)
  115 format(/' ia =',i3,', irregular solution')
  120 format(/' ibord phiato(u,1)  phiato(d,1)  phiato(u,2)  ', 'phiato(d,2)   L =',i3,', m =',i3)
  130 format(/' ibord  phiato(u)    phiato(d)   L =',i3,', m =',i3)
  140 format(/' ibord',6x,' phiato(u,1)',13x,'  phiato(d,1)',13x, '  phiato(u,2)  ',13x,'phiato(d,2)   L =',i3,', m =',i3)
  150 format(/' ibord',6x,'  phiato(u)',13x,'    phiato(d)   L =',i3, ', m =',i3)
  160 format(i5,1p,8e13.4)
end

!***********************************************************************

! Calculation of density of state

subroutine Cal_dens(Cal_xanes,Classic_irreg,Density_comp,drho_self,Ecinetic,Eimag,Energ,Enervide,Full_atom, &
            Full_potential,Harm_cubic,Hubb,Hubb_diag,Hubb_diag_abs,iaabsi,iaprotoi,icheck,itypei,itypepr,lla2_state, &
            lmax_pot,lmaxat, &
            m_hubb,mpinodes,mpirank,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nlmagm,nlm_pot, &
            nrato,nrm,nrm_self,nspin,nspino,nspinp,ntype,numat,rato,Relativiste,Renorm,Rmtg,Rmtsd,Solsing,Solsing_only, &
            Spinorbite,State_all,Statedens,Statedens_i,Taull,V_hubb,V_hubb_abs,V_intmax,V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer ia, iaabsi, iapr, icheck, ipr, ir, isp, isp1, isp2, it, lla2_state, lm, lmax, lmax_pot, &
    m_hubb, mpinodes, mpirank, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, &
    n_atom_proto, natome, nlm_pot, nlma, nlmagm, nr, nrs, nrm, nrm_self, nspin, nspino, nspinp, ntype, Z

  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,natome):: Taull
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb_abs
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb
  complex(kind=db), dimension(:,:,:,:), allocatable:: V_hubb_t
  complex(kind=db), dimension(:,:,:,:), allocatable:: Taulla

  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, lmaxat
  integer, dimension(natome):: iaprotoi, itypei

  logical:: Absorbeur, Cal_xanes, Classic_irreg, Density_comp, Full_atom, Full_potential, Harm_cubic, Hubb_a, Hubb_d, &
    Hubb_diag_abs, Relativiste, Renorm, Self, Solsing, Solsing_only, Spinorbite, State_all, Ylm_comp, Comp_to_real
  logical, dimension(0:ntype):: Hubb
  logical, dimension(n_atom_0:n_atom_ind):: iapr_done
  logical, dimension(n_atom_0_self:n_atom_ind_self):: Hubb_diag

  real(kind=db):: Eimag, Energ, Enervide, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(0:n_atom_proto):: Rmtg, Rmtsd
  real(kind=db), dimension(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1):: Statedens, Statedens_i
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vrato
  real(kind=db), dimension(nrm_self,nlm_pot,nspinp,n_atom_0_self:n_atom_ind_self,0:mpinodes-1):: drho_self
  real(kind=db), dimension(:), allocatable:: r
  real(kind=db), dimension(:,:,:), allocatable:: drho, Vrato_t
  real(kind=db), dimension(:,:,:,:), allocatable:: State, State_i

  Self = .not. Cal_xanes

  if( Self ) drho_self(:,:,:,:,mpirank) = 0._db

  Statedens(:,:,:,:,:,mpirank) = 0._db
  Statedens_i(:,:,:,:,:,mpirank) = 0._db

  iapr_done(:) = .false.

  boucle_ia: do ia = 1,natome

    Absorbeur = ( ia == iaabsi .and. Cal_xanes )
    if( Cal_xanes .and. .not. ( State_all .or. Absorbeur ) ) cycle

    ipr = iaprotoi(ia)
    it = itypepr(ipr)
    Hubb_a = Hubb(it)

    Z = numat( it )
    if( Z < 19 ) then
      lmax = min(2,lmaxat(ipr))
    elseif( Z > 18 .and. Z < 55 ) then
      lmax = min(3,lmaxat(ipr))
    else
      lmax = min(4,lmaxat(ipr))
    endif
    nlma = (lmax + 1)**2

    if( Full_atom ) then
      iapr = ia
    else
      iapr = ipr
      if( iapr_done(iapr) ) cycle
    endif
    it = itypei(ia)

    do ir = 1,nrato(it)
      if( rato(ir,it) > Rmtsd(ipr) + eps10 ) exit
    end do
    nr = ir + 1

    if( Self ) then
      nrs = nr
    else
      nrs = 0
    endif
    allocate( r(nr) )
    allocate( drho(nrs,nlm_pot,nspinp) )
    allocate( State(nlma,nspinp,nlma,nspinp))
    allocate( State_i(nlma,nspinp,nlma,nspinp))
    allocate( Vrato_t(nr,nlm_pot,nspin) )
    allocate( V_hubb_t(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp) )
    allocate( Taulla(nlma,nspinp,nlma,nspinp) )
    State(:,:,:,:) = 0._db
    State_i(:,:,:,:) = 0._db
    drho(:,:,:) = 0._db
    r(1:nr) = rato(1:nr,it)
    Vrato_t(1:nr,:,:) = Vrato(1:nr,:,:,iapr)
    if( Hubb(it) .and. ( iapr <= n_atom_ind_self .or. Absorbeur ) ) then 
   ! ".or. Absorber" is because in case of Extract, n_atom_ind_self == 0, and the test becomes false
      Hubb_a = .true.
      if( Absorbeur ) then
        V_Hubb_t(:,:,:,:) = V_Hubb_abs(:,:,:,:)
        Hubb_d = Hubb_diag_abs
      else
        V_Hubb_t(:,:,:,:) = V_Hubb(:,:,:,:,iapr)
        Hubb_d = Hubb_diag(iapr)
      endif
    else
      Hubb_a = .false.
      Hubb_d = .true.
    end if

! Comp_to_real = .true. corresponds to calculation with complex harmonics that we project to real harmonics
    if( ( Cal_xanes .and. Density_comp ) .or. .not. Cal_xanes ) then
      Comp_to_real = .false.
    else
      Comp_to_real = Ylm_comp
    endif

    do isp1 = 1,nspinp
      do lm = 1,nlma
        do isp2 = 1,nspinp
          Taulla(1:nlma,isp1,lm,isp2) = Taull(1:nlma,isp1,lm,isp2,ia)
        end do
      end do
    end do

    call Radial_sd(Classic_irreg,drho,Ecinetic,Eimag,Energ,Enervide,Full_potential,Harm_cubic,Hubb_a,Hubb_d,ia,icheck,lmax, &
         lmax_pot,m_hubb,nlm_pot,nlma,nr,nrs,nspin,nspino,nspinp,Z,r,Relativiste,Renorm,Rmtg(ipr), &
         Rmtsd(ipr),Self,Solsing,Solsing_only,Spinorbite,State,State_i,Taulla,V_hubb_t,V_intmax,V0bd,Vrato_t,Ylm_comp,Comp_to_real)

    deallocate( r, Taulla, V_hubb_t, Vrato_t )

    do isp1 = 1,nspinp
      do lm = 1,nlma
        do isp2 = 1,nspinp
          Statedens(1:nlma,isp1,lm,isp2,iapr,mpirank) = State(1:nlma,isp1,lm,isp2)
          Statedens_i(1:nlma,isp1,lm,isp2,iapr,mpirank) = State_i(1:nlma,isp1,lm,isp2)
        end do
      end do
    end do

    if( Self ) then
      do isp = 1,nspinp
        do lm = 1,nlm_pot
          drho_self(1:nr,lm,isp,iapr,mpirank) = drho(1:nr,lm,isp)
        end do
      end do
    endif

    deallocate( State, State_i, drho )

    iapr_done(iapr) = .true.

  end do boucle_ia

  return
end

!**********************************************************************

! Calculation of radial integrals and irregular solution for the density of states
! Comp_to_real = .true. corresponds to calculation with complex harmonics that we project to real harmonics
! Called by Cal_dens

subroutine Radial_sd(Classic_irreg,drho,Ecinetic,Eimag,Energ,Enervide,Full_potential,Harm_cubic,Hubb_a,Hubb_d,ia,icheck,lmax, &
         lmax_pot,m_hubb,nlm_pot,nlma,nr,nrs,nspin,nspino,nspinp,numat,r,Relativiste,Renorm,Rmtg, &
         Rmtsd,Self,Solsing,Solsing_only,Spinorbite,State,State_i,Taull,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp,Comp_to_real)

  use declarations
  implicit none

  integer:: ia, icheck, isp1, isp2, L, l_hubbard, lfin, lmax, lmax_pot, m, m_hubb, m1, m2, nlm, nlm_pot, nlm1, nlm2, nlma, nr, &
    nrs, nrmtg, nrmtsd, nspin, nspino, nspinp, numat

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlma,nspinp,nlma,nspinp):: Taull, Taull_r
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb, V_hubb_r 
  complex(kind=db), dimension(:,:), allocatable:: Mat_c, Mat_r 
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau

  logical:: Classic_irreg, Diagonale, Ecomp, Full_potential, Harm_cubic, Hubb_a, Hubb_d, Hubb_m, Radial_comp, Relativiste, Renorm, &
    Self, Solsing, Solsing_only, Spinorbite, Ylm_comp, Comp_to_real

  real(kind=db):: Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nr,nspin):: g0, gm, gmi, gp, gpi
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato
  real(kind=db), dimension(nrs,nlm_pot,nspinp):: drho
  real(kind=db), dimension(nlma,nspinp,nlma,nspinp):: State, State_i
  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur

  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

  if( icheck > 1 ) then
    write(3,110)
    write(3,120) ia, numat, Rmtg*bohr, Rmtsd*bohr, Energ*rydb
    write(3,130) Ecinetic(:)*rydb
    write(3,140) V0bd(:)*rydb
    write(3,150) konde(:)
  endif

  if( Full_potential ) then
    nlm = ( lmax + 1 )**2
  else
    nlm = 1
  endif
  allocate( V(nr,nlm,nlm,nspin) )

  call mod_V(icheck,lmax,lmax_pot,nlm,nlm_pot,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)

  gpi(:,:) = 1 / gp(:,:)
  if( Solsing ) gmi(:,:) = 1 / gm(:,:)

  if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10 .or. Ecinetic(nspin) < eps10 ) then
    Ecomp = .true.
  else
    Ecomp = .false.
  endif
  
!  Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )
  Radial_comp = Ecomp

  if( Full_potential ) then
    lfin = 0
  else
    lfin = lmax
  endif

  do L = 0,lfin

    if( Hubb_a .and. L == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = nlm ! = ( lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*L + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*L + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif
    Diagonale = nlm2 == 1

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    V_hubb_r(:,:,:,:) = V_hubb(:,:,:,:)
    Taull_r(:,:,:,:) = Taull(:,:,:,:)

    if( .not. Full_potential .and. Hubb_m .and. .not. Hubb_d .and. Ylm_comp ) then
    
      allocate( Mat_c(-L:L,-L:L) )
      allocate( Mat_r(-L:L,-L:L) )
 
      do isp1 = 1,nspinp
        do isp2 = 1,nspinp
          Mat_c(:,:) = V_hubb(:,:,isp1,isp2)
          Call Trans_HarmoC_to_R(L,Mat_c,Mat_r)
          V_hubb_r(:,:,isp1,isp2) = Mat_r(:,:)
        end do
      end do

      do isp1 = 1,nspinp
        do isp2 = 1,nspinp
          do m = -L,L
            Mat_c(:,m) = Taull(L**2+1:(L+1)**2,isp1,L**2+L+1+m,isp2)
          end do
          Call Trans_HarmoC_to_R(L,Mat_c,Mat_r)
          do m1 = -L,L
            do m2 = -L,L
              if( ( m1 >= 0 .and. m2 < 0 ) .or. ( m1 < 0  .and. m2 >= 0 ) ) then 
                Taull_r(L**2+L+1+m1,isp1,L**2+L+1+m2,isp2) = - Mat_r(m1,m2)
              else
                Taull_r(L**2+L+1+m1,isp1,L**2+L+1+m2,isp2) = Mat_r(m1,m2)
              endif
            end do
          end do
        end do
      end do

      deallocate( Mat_c, Mat_r )

      call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,L,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb_r)

! Calculation of the radial integral
      if( .not. Solsing_only ) &
        call Radial_matrix_sd(Diagonale,drho,Full_potential,Harm_cubic,icheck,L,lmax,nlm_pot,nlm1,nlm2,nlma,nr,nrs,nrmtsd, &
           nspin,nspino,nspinp,r,Radial_comp,Rmtsd,Spinorbite,State,State_i,Self,Taull_r,ui,ur,.false.,.false.)
          
    else

      call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,L,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

! Calculation of the radial integral
      if( .not. Solsing_only ) &
        call Radial_matrix_sd(Diagonale,drho,Full_potential,Harm_cubic,icheck,L,lmax,nlm_pot,nlm1,nlm2,nlma,nr,nrs,nrmtsd, &
           nspin,nspino,nspinp,r,Radial_comp,Rmtsd,Spinorbite,State,State_i,Self,Taull_r,ui,ur,Ylm_comp,Comp_to_real)

    endif
    
! Calculation of irregular solution 
    if( Solsing ) call Cal_Solsing_sd(Classic_irreg,drho,Ecomp,Eimag,f2,Full_potential,g0,gmi,gp,gso,Harm_cubic,Hubb_a,Hubb_d, &
          icheck,konde,L,lmax,m_hubb,nlm,nlm_pot,nlm1,nlm2,nlma,nr,nrmtsd,nrs,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Self, &
          Spinorbite,State,Tau,ui,ur,V,V_hubb,Comp_to_real)

    deallocate( Tau )
    deallocate( ui, ur )

  end do   ! end of loop over L

  deallocate( V )

  return
  110 format(/' ---- Radial_sd ----',100('-'))
  120 format(/'     Atom =',i3,/'        Z =',i3,/ '     Rmtg =',f10.5,' A,  Rmtsd =',f10.5,' A',/ &
              '    Energ =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format('    V0bd  =',2f10.3)
  150 format('    konde =',2f12.5)
end

!********************************************************************************************************************************

subroutine Trans_HarmoC_to_R(L,Mat_c,Mat_r)

  use declarations
  implicit none

  integer:: is, L, m
  
  complex(kind=db), dimension(-L:L,-L:L):: Mat_c, Mat_r, Trans, Transi 
      
  Trans(:,:) = (0._db, 0._db)
  Transi(:,:) = (0._db, 0._db)

  do m = -L,L
    if( m > 0 ) then
      is = (-1)**m
      Trans(m,m) = cmplx( is * sqrt_1o2, 0._db, db)
      Trans(m,-m) = cmplx( 0._db, is * sqrt_1o2, db)
    elseif( m == 0 ) then
      Trans(0,0) = (1._db, 0._db)
    else
      Trans(m,-m) = cmplx(sqrt_1o2, 0._db, db)
      Trans(m,m) = cmplx(0._db, -sqrt_1o2, db)
    endif
  end do
  
  Transi = Conjg( Transpose( Trans ) )

  Mat_r = Matmul( Transi, Matmul( Mat_c, Trans ) )

  return
end

!********************************************************************************************************************************

! Calculation of the radial integral for the density of state
! Comp_to_real = .true. corresponds to the calculation with complex harmo that we project to real harmo 

subroutine Radial_matrix_sd(Diagonale,drho,Full_potential,Harm_cubic,icheck,L,lmax,nlm_pot,nlm1,nlm2,nlma,nr,nrs,nrmtsd, &
           nspin,nspino,nspinp,r,Radial_comp,Rmtsd,Spinorbite,State,State_i,Self,Taull,ui,ur,Ylm_comp,Comp_to_real)

  use declarations
  implicit none

  integer:: i_k1, i_k2, icheck, ir, is1, is2, isg1, isg2, iso, iso1, iso2, isp, isp1, isp2, iss1, iss2, &
    L, L1, L2, l3, lm, lm0, &
    lm1, lm2, lma1, lma2, lmax, Lp1, Lp2, m1, m2, m3, m_k1, m_k2, mp1, mp2, m_r1, m_r2, mv1, mv2, n_rk1, n_rk2, n1, n2, &
    nlm_pot, nlm1, nlm2, nlma, np1, np2, nr, nrmtsd, nrs, nspin, nspino, nspinp

! The volatile statetement prevent for some compilation optimization.
! On linux, at times without it, one gets segmentation fault erros message during the running with magnetic system (YJ 2018/01/25)
  integer, volatile:: n_k1, n_k2
!  integer:: n_k1, n_k2
  
  complex(kind=db):: c_harm, c_harm1, c_harm2, rof_sd, Sta_e, Sta_s
  complex(kind=db), dimension(nspinp):: rof_sd_l
  complex(kind=db), dimension(nrmtsd):: fc
  complex(kind=db), dimension(nrmtsd,nspinp):: fc_l
  complex(kind=db), dimension(nlma,nspinp,nlma,nspinp):: Taull

  logical:: Diagonale, Full_potential, Harm_cubic, Radial_comp, Self, Spinorbite, Ylm_comp, Comp_to_real

  real(kind=db), parameter:: sqrt_4pi = sqrt( 4*pi )

  real(kind=db):: c_cubic1, c_cubic2, f_integr3, g, Gaunt_r, gauntcp, radlr, radli, Rmtsd, Sta_i, Sta_r
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nrmtsd):: r2, rr, fcr, fci, ri2
  real(kind=db), dimension(nrs,nlm_pot,nspinp):: drho
  real(kind=db), dimension(nlma,nspinp,nlma,nspinp):: State, State_i
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur

  if( Full_potential ) then
    lm0 = 0
  else
    lm0 = L**2
  endif

  rr(1:nrmtsd) = r(1:nrmtsd)
  r2(1:nrmtsd) = r(1:nrmtsd)**2
  ri2(1:nrmtsd) = 1 / ( quatre_pi * r2(1:nrmtsd) )

  if( nlm1 == 1 ) then
    do isp = 1,nspinp
      iso = min( isp, nspino )
      if( Radial_comp ) then
        fcr(1:nrmtsd) = r2(1:nrmtsd) * ( ur(1:nrmtsd,1,1,isp,iso) * ur(1:nrmtsd,1,1,isp,iso) &
                                       - ui(1:nrmtsd,1,1,isp,iso) * ui(1:nrmtsd,1,1,isp,iso) )
        radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
        fci(1:nrmtsd) = r2(1:nrmtsd) * ( ur(1:nrmtsd,1,1,isp,iso) * ui(1:nrmtsd,1,1,isp,iso) &
                                       + ui(1:nrmtsd,1,1,isp,iso) * ur(1:nrmtsd,1,1,isp,iso) )
        radli = f_integr3(rr,fci,1,nrmtsd,Rmtsd)
      else
        fcr(1:nrmtsd) = r2(1:nrmtsd) * ur(1:nrmtsd,1,1,isp,iso) * ur(1:nrmtsd,1,1,isp,iso)
        radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
        fci(:) = 0._db
        radli = 0._db
      endif
      rof_sd_l(isp) = cmplx( radlr, radli, db )
      if( Self ) fc_l(1:nrmtsd,isp) = ri2(1:nrmtsd) * cmplx( fcr(1:nrmtsd), fci(1:nrmtsd), db )

    end do
  endif

  if( icheck > 2 ) then
    if( nlm2 > 1 ) then
      write(3,110)
    elseif( Spinorbite ) then
      write(3,120)
    else
      write(3,130)
    endif
  endif

  n_k1 = 0
  do L1 = 0,lmax
    if( .not. Full_potential .and. L /= L1 ) cycle
    do m_k1 = -L1,L1
      n_k1 = n_k1 + 1
      lm1 = lm0 + n_k1

      if( .not. Harm_cubic .or. L1 /= 3 .or. ( m_k1 == -3 .or. m_k1 == 0 .or. m_k1 == 3 ) ) then
        n_rk1 = 1
      else
        n_rk1 = 2
      endif

      do isp1 = 1,nspinp
      
        n_k2 = 0
        do L2 = 0,lmax
          if( .not. Full_potential .and. L /= L2 ) cycle
          do m_k2 = -L2,L2
            n_k2 = n_k2 + 1
            lm2 = lm0 + n_k2

            if( .not. Harm_cubic .or. L2 /= 3 .or. ( m_k2 == -3 .or. m_k2 == 0 .or. m_k2 == 3 ) ) then
              n_rk2 = 1
            else
              n_rk2 = 2
            endif
        
            do isp2 = 1,nspinp
        
              do i_k1 = 1,n_rk1

                call Trans_TtoK(c_cubic1,Harm_cubic,i_k1,L1,m_k1,m_r1)

                do isg1 = -1,1,2
 
                  if( ( .not. Comp_to_real .or. m_r1 == 0 ) .and. isg1 == -1 ) cycle
       
                  call Trans_RtoC(c_harm1,isg1,m1,m_r1,Comp_to_real)

                  np1 = 0
                  do Lp1 = 0,lmax
                    if( .not. Full_potential .and. Lp1 /= L1 ) cycle

                    do mp1 = -Lp1, Lp1
                      if( nlm2 == 1 .and. mp1 /= m1 ) cycle
                      if( nlm2 == 1 ) then
                        np1 = 1
                      else
                        np1 = np1 + 1
                      endif

                      do iso1 = 1,nspino
  
                        if( Spinorbite ) then
                          is1 = iso1
                        else
                          is1 = isp1
                        endif
                        if( Diagonale ) then
                          iss1 = isp1
                        else
                          iss1 = isp1
                        endif

                        if( Spinorbite .and. Diagonale ) then
!                          mv1 = m1 + iso1 - isp1
                          mv1 = mp1 + iso1 - isp1
                          if( mv1 > L1 .or. mv1 < - L1 ) cycle
                        else
!                          mv1 = m1
                          mv1 = mp1
                        endif                  

                        if( Full_potential ) then
                          n1 = lma1
                        elseif( nlm1 == 1 ) then
                          n1 = 1
                        else
                          n1 = L1 + 1 + m1 
                        endif

                        lma1 = L1**2 + L1 + 1 + mv1

                        do i_k2 = 1,n_rk2
        
                          call Trans_TtoK(c_cubic2,Harm_cubic,i_k2,L2,m_k2,m_r2)
        
                          do isg2 = -1,1,2
        
                            if( ( .not. Comp_to_real .or. m_r2 == 0 ) .and. isg2 == -1 ) cycle
                   
                            call Trans_CtoR(c_harm2,isg2,m2,m_r2,Comp_to_real)
        
!                            c_harm = c_harm1 * c_cubic1 * conjg( c_harm2 * c_cubic2 )
                            c_harm = c_harm1 * c_cubic1 * c_harm2 * c_cubic2
        
                            np2 = 0
                            do Lp2 = 0,lmax
                              if( .not. Full_potential .and. Lp2 /= L2 ) cycle
                              if( .not. Full_potential .and. Lp2 /= Lp1 ) cycle
        
                              do mp2 = -Lp2,Lp2
                                if( nlm2 == 1 .and. mp2 /= m2 ) cycle
                                if( nlm2 == 1 ) then
                                  np2 = 1
                                else
                                  np2 = np2 + 1
                                endif
        
                                do iso2 = 1,nspino
                                  if( isp1 /= isp2 ) cycle ! in a scalar product there is no spin-crossing (but there is solution crossing)
                                  if( .not. Spinorbite .and. isp1 /= isp2 ) cycle
 
                                  if( Spinorbite ) then
                                    is2 = iso2
                                  else
                                    is2 = isp2
                                  endif
                                  if( Diagonale ) then
                                    iss2 = isp2
                                  else
                                    iss2 = isp2
                                  endif

                                  if( Spinorbite .and. Diagonale ) then
                                    mv2 = mp2 + iso2 - isp2
!                                    mv2 = m2 + iso2 - isp2
                                    if( mv2 > L2 .or. mv2 < - L2 ) cycle
                                  else
!                                    mv2 = m2
                                    mv2 = mp2
                                  endif                  

                                  if( Full_potential ) then
                                    n2 = lma2
                                  elseif( nlm1 == 1 ) then
                                    n2 = 1
                                  else
                                    n2 = L2 + 1 + m2
                                  endif
                                  
                                  lma2 = L2**2 + L2 + 1 + mv2

                                  if( nlm1 == 1 ) then
        
                                    rof_sd = rof_sd_l(isp1)
                                    if( Self .and. isp1 == isp2 .and. lm1 == lm2 ) fc(:) = fc_l(:,isp1)
        
                                  else
                                    
                                    if( nlm2 == 1 ) then
        
                                      if( Radial_comp ) then
                                        fcr(1:nrmtsd) = r2(1:nrmtsd) &
                                                      * ( ur(1:nrmtsd,n1,np1,isp1,iso1) * ur(1:nrmtsd,n2,np2,isp2,iso2) &
                                                        - ui(1:nrmtsd,n1,np1,isp1,iso1) * ui(1:nrmtsd,n2,np2,isp2,iso2) )
                                        radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
                                        fci(1:nrmtsd) = r2(1:nrmtsd) &
                                                      * ( ur(1:nrmtsd,n1,np1,isp1,iso1) * ui(1:nrmtsd,n2,np2,isp2,iso2) &
                                                        + ui(1:nrmtsd,n1,np1,isp1,iso1) * ur(1:nrmtsd,n2,np2,isp2,iso2) )
                                        radli = f_integr3(rr,fci,1,nrmtsd,Rmtsd)
                                      else
                                        fcr(1:nrmtsd) = r2(1:nrmtsd) * ur(1:nrmtsd,n1,np1,isp1,iso1) &
                                                                     * ur(1:nrmtsd,n2,np2,isp2,iso2)
                                        radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
                                        fci(:) = 0._db
                                        radli = 0._db
                                      endif
                                      
                                    else ! seems not useful because u(m,mp) = conjg( u(mp,m) ) at least when Eimag = 0

                                      if( Radial_comp ) then
                                        fcr(1:nrmtsd) = r2(1:nrmtsd) &
                                                      * ( ur(1:nrmtsd,np1,n1,isp1,iso1) * ur(1:nrmtsd,np2,n2,isp2,iso2) &
                                                        - ui(1:nrmtsd,np1,n1,isp1,iso1) * ui(1:nrmtsd,np2,n2,isp2,iso2) )
                                        radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
                                        fci(1:nrmtsd) = r2(1:nrmtsd) &
                                                      * ( ur(1:nrmtsd,np1,n1,isp1,iso1) * ui(1:nrmtsd,np2,n2,isp2,iso2) &
                                                        + ui(1:nrmtsd,np1,n1,isp1,iso1) * ur(1:nrmtsd,np2,n2,isp2,iso2) )
                                        radli = f_integr3(rr,fci,1,nrmtsd,Rmtsd)
                                      else
                                        fcr(1:nrmtsd) = r2(1:nrmtsd) * ur(1:nrmtsd,np1,n1,isp1,iso1) &
                                                                     * ur(1:nrmtsd,np2,n2,isp2,iso2)
                                        radlr = f_integr3(rr,fcr,1,nrmtsd, Rmtsd)
                                        fci(:) = 0._db
                                        radli = 0._db
                                      endif
                                    
                                    endif
        
                                    rof_sd = cmplx( radlr, radli, db)
        
                                    if( Self .and. isp1 == isp2 ) fc(1:nrmtsd) &
                                                                       = cmplx( fcr(1:nrmtsd),fci(1:nrmtsd),db) * ri2(1:nrmtsd )
        
                                  endif

                                  Sta_e = c_harm * rof_sd * Taull(lma1,is1,lma2,is2)
                                  Sta_s = conjg( c_harm ) * rof_sd * Taull(lma2,is2,lma1,is1)
        
                                  Sta_r = 0.5_db * real( Sta_e - Sta_s, db)
                                  Sta_i = 0.5_db * aimag( Sta_e + Sta_s )
        
                                  State(lm1,iss1,lm2,iss2) = State(lm1,iss1,lm2,iss2) - Sta_i
                                  State_i(lm1,iss1,lm2,iss2) = State_i(lm1,iss1,lm2,iss2) + Sta_r

                                  if( Self .and. iss1 == iss2 ) then

                                    lm = 0
                                    boucle_L3: do L3 = 0,lmax
                                      if( .not. Full_potential .and. l3 > 0 ) cycle                                  
!                                    if( ( nlm1 == 1 .or. Diagonale .or. Comp_to_real ) .and. L3 /= 0 ) cycle
!                                    if( .not. Full_potential .and. L3 /= L1 .and. L3 /= 0 ) cycle

                                      do m3 = -L3,L3
                                        lm = lm + 1
!                                        if( Comp_to_real ) then
 !                                         g = 1._db
                                        if( Ylm_comp ) then
                                          g = sqrt_4pi * gauntcp(L1,m1,L2,m2,L3,m3)
                                        else
                                          g = sqrt_4pi * Gaunt_r(L1,m1,L2,m2,L3,m3)
                                        endif

                                        if( abs( g ) < eps10 ) cycle 
                                      
                                        do ir = 1,nrmtsd
                                          Sta_e = c_harm * fc(ir) * Taull(lma1,is1,lma2,is2)
                                          Sta_s = conjg(c_harm) * fc(ir) * Taull(lma2,is2,lma1,is1)
                                          Sta_i = 0.5_db * aimag( Sta_e + Sta_s)
                                          drho(ir,lm,iss1) = drho(ir,lm,iss1) - g * Sta_i        
                                        end do
                                      
                                      end do
        
                                    end do boucle_L3

                                  endif
        
                                  if( icheck > 2 .and. abs( Taull(lma1,is1,lma2,is2) ) > 1.e-15_db ) then
                                    if( nlm2 > 1 ) then
                                      write(3,140) L1, m_k1, m_r1, m1, isp1, L2, m_k2, m_r2, m2, isp2, Lp1, mp1, &
                                        iso1, Lp2, mp2, iso2, rof_sd, Taull(lma1,is1,lma2,is2), State(lm1,iss1,lm2,iss2), &
                                        State_i(lm1,iss1,lm2,iss2)
                                    elseif( Spinorbite ) then
                                      write(3,150) lm1, iss1, lm2, iss2, lma1, is1, lma2, is2, L1, m_k1, m_r1, m1, isp1, &
                                        L2, m_k2, m_r2, m2, isp2, iso1, iso2, rof_sd, Taull(lma1,is1,lma2,is2), &
                                        Taull(lma2,is2,lma1,is1), State(lm1,iss1,lm2,iss2), State_i(lm1,iss1,lm2,iss2)
                                    else
                                      write(3,160) L1, m_k1, m_r1, m1, isp1, L2, m_k2, m_r2, m2, isp2, rof_sd, &
                                        Taull(lma1,is1,lma2,is2), State(lm1,iss1,lm2,iss2), State_i(lm1,iss1,lm2,iss2)
                                    endif
                                  endif
        
                                end do  ! boucle iso2
                              end do  ! boucle mp2
                            end do  ! boucle Lp2
                          end do  ! boucle isg2
                        end do  ! boucle ik_2

                      end do  ! boucle i_k1
                    end do  ! boucle iso1
                  end do  ! boucle mp1
                end do  ! boucle isg1
              end do  ! boucle i_k1

            end do  ! boucle isp2
          end do  ! boucle m_k2
        end do  ! boucle L2

      end do  ! boucle isp1
    end do  ! boucle m_k1
  end do  ! boucle L1

  if( icheck > 1 ) then
    if( .not. Full_potential ) then
      write(3,165) (( L, m1, isp, m1 = -L,L), isp = 1,nspinp )
    else
      write(3,166)
    endif
    do isp1 = 1,nspinp
      n1 = 0
      do L1 = 0,lmax
        if( .not. Full_potential .and. L1 /= L ) cycle
        do m1 = -L1,L1
          n1 = n1 + 1
          lm1 = lm0 + n1
          if( .not. Full_potential ) then
            write(3,180) L1, m1, isp1, ((Taull(lm1,isp1,lm2,isp2), lm2 = lm0+1,lm0+2*L1+1), isp2 = 1,nspinp )
          else
            write(3,180) L1, m1, isp1, (( Taull(lm1,isp1,lm2,isp2), lm2 = 1,nlma), isp2 = 1,nspinp )
          endif
        end do
      end do
    end do

    if( .not. Full_potential ) then
      write(3,170) (( L, m1, isp, m1 = -L,L), isp = 1,nspinp )
    else
      write(3,175)
    endif
    do isp1 = 1,nspinp
      n1 = 0
      do L1 = 0,lmax
        if( .not. Full_potential .and. L1 /= L ) cycle
        do m1 = -L1,L1
          n1 = n1 + 1
          lm1 = lm0 + n1
          if( .not. Full_potential ) then
            write(3,180) L1, m1, isp1, (( State(lm1,isp1,lm2,isp2), State_i(lm1,isp1,lm2,isp2), lm2 = lm0+1,lm0+2*L1+1), &
                                                                                                isp2 = 1,nspinp )
          else
            write(3,180) L1, m1, isp1, (( State(lm1,isp1,lm2,isp2), State_i(lm1,isp1,lm2,isp2), lm2 = 1,nlma), isp2 = 1,nspinp )
          endif
        end do
      end do
    end do
  endif

  if( icheck > 2 .and. Self ) then
    if( nspin == 1 ) then
      write(3,190) ( lm, lm = 1,nlm_pot )
    else
      write(3,200) ( lm, lm, lm = 1,nlm_pot )
    endif
    do ir = 1,nrmtsd
      write(3,210) r(ir)*bohr, ( drho(ir,lm,:), lm = 1,nlm_pot )
    end do
  endif

  return
  110 format(/' Radial integral for the density of state:',/ &
              ' L1 mk1 mr1 m1 isp1 L2 mk2 mr2 m2 isp2 Lp1 mp1 iso1 Lp2 mp2 iso2',12x,'rof_sd',23x,'Taull',17x,'State',7x,'State_i')
  120 format(/' Radial integral for the density of state:',/ &
              ' lm1 iss1 lm2 iss2 lma1 is1 lma2 is2 L1 mk1 mr1 m1 isp1 L2 mk2 mr2 m2 isp2 iso1 iso2',12x,'rof_sd',23x,'Taull',22x, &
                   'Taull',18x,'State',7x,'State_i')
  130 format(/' Radial integral for the density of state:',/ &
              ' L1 mk1 mr1 m1 isp1 L2 mk2, mr2 m2 isp2',12x,'rof_sd',23x,'Taull',17x,'State',7x,'State_i')
  140 format(2(2i3,3i4,1x),2(i3,i4,i5,1x),1p,3(2x,2e13.5))
  150 format(i3,i5,i4,i5,i5,i4,i5,i4,1x,2(2i3,3i4,1x),2(i4,1x),1p,4(2x,2e13.5))
  160 format(2(2i3,3i4,1x),1p,3(2x,2e13.5))
  165 format(/' Multiple scattering amplitude:',/ &
              '  L  m isp ',100(8x,3i3,10x))
  166 format(/' Multiple scattering amplitude:',/ &
              '  L  m isp ',6x,'Taull_r',6x,'Taull_i')
  170 format(/' Radial integral for the density of state before singul:',/ &
              '  L  m isp ',100(8x,3i3,10x))
  175 format(/' Radial integral for the density of state before singul:',/ &
              ' L1 m1 isp1',7x,'State',7x,'State_i')
  180 format(3i3,1p,28(1x,2e13.5))
  190 format(/' drho before singul', /'   Radius   ',100(6x,i2,5x))
  200 format(/' drho before singul', /'   Radius   ',100(3x,i2,' up',8x,i2,' dn',5x))
  210 format(1p,300e13.5)
end

!******************************************************************************************

! Transformation of real (or tesseral) spherical harmonics to cubic spherical harmonics

subroutine Trans_TtoK(c_cubic,Harm_cubic,i_k,L,m_k,m_r)
    
  use declarations
  implicit none

  integer:: i_k, L, m_k, m_r

  logical:: Harm_cubic
  
  real(kind=db):: c_cubic 
  
  if( Harm_cubic .and. L >= 1 .and. L <= 3 ) then

    Select case(L)
      Case(1)

        select case(m_k)
          case(-1)
            m_r = 1; c_cubic = 1._db
          case(0)
            m_r = -1; c_cubic = 1._db
          case(1)
            m_r = 0; c_cubic = 1._db
        end select
  
      Case(2)

        select case(m_k)
          case(-2)
            m_r = 2; c_cubic = 1._db
          case(-1)
            m_r = 0; c_cubic = 1._db
          case(0)
            m_r = -1; c_cubic = 1._db
          case(1)
            m_r = 1; c_cubic = 1._db
          case(2)
            m_r = -2; c_cubic = 1._db
        end select
      
      Case(3)
      
        select case(m_k)
          case(-3)
            m_r = -2; c_cubic = 1._db
          case(-2)
            if( i_k == 1 ) then
              m_r = 1; c_cubic = - sqrt( 3._db / 8 )
            else
              m_r = 3; c_cubic = sqrt( 5._db / 8 )
            endif
          case(-1)
            if( i_k == 1 ) then
              m_r = -3; c_cubic = - sqrt( 5._db / 8 )
            else
              m_r = -1; c_cubic = - sqrt( 3._db / 8 )
            endif
          case(0)
            m_r = 0; c_cubic = 1._db
          case(1)
            if( i_k == 1 ) then
              m_r = 1; c_cubic = - sqrt( 5._db / 8 )
            else
              m_r = 3; c_cubic = - sqrt( 3._db / 8 )
            endif
          case(2)
            if( i_k == 1 ) then
              m_r = -3; c_cubic = - sqrt( 3._db / 8 )
            else
              m_r = -1; c_cubic = sqrt( 5._db / 8 )
            endif
          case(3)
            m_r = 2; c_cubic = 1._db
        end select
  
    end select 

  else
  
    c_cubic = 1._db
    m_r = m_k
  
  endif

  return
end

!***********************************************************************

! Transformation of complex spherical harmonics to real spherical (or tesseral) harmonics

subroutine Trans_CtoR(c_harm,isg,m,m_r,Comp_to_real)
    
  use declarations
  implicit none

  integer:: isg, m, m_r

  logical:: Comp_to_real
  
  complex(kind=db):: c_harm 
  
  if( Comp_to_real ) then
    if( m_r < 0 ) then
      if( isg == 1 ) then
        c_harm = cmplx( 0._db, sqrt_1o2, db )
      else
        c_harm = cmplx( 0._db, - (-1)**m_r * sqrt_1o2, db )
      endif
    elseif( m_r == 0 ) then
      c_harm = ( 1._db, 0._db )
    else
      if( isg == 1 ) then
        c_harm = cmplx( (-1)**m_r * sqrt_1o2, 0._db, db )
      else
        c_harm = cmplx( sqrt_1o2, 0._db, db )
      endif
    endif
    m = isg * m_r
  else
    m = m_r
    c_harm = ( 1._db, 0._db )
  endif

  return
end

!***********************************************************************

! Transformation of complex spherical harmonics to real spherical (or tesseral) harmonics

subroutine Trans_RtoC(c_harm,isg,m,m_r,Comp_to_real)
    
  use declarations
  implicit none

  integer:: isg, m, m_r

  logical:: Comp_to_real
  
  complex(kind=db):: c_harm 
  
  if( Comp_to_real ) then
    m = isg * m_r
    if( m < 0 ) then
      if( isg == 1 ) then
        c_harm = cmplx( 0._db, - sqrt_1o2, db )
      else
        c_harm = cmplx( sqrt_1o2, 0._db, db )
      endif
    elseif( m_r == 0 ) then
      c_harm = ( 1._db, 0._db )
    else
      if( isg == 1 ) then
        c_harm = cmplx( (-1)**m_r * sqrt_1o2, 0._db, db )
      else
        c_harm = cmplx( 0._db, (-1)**m_r * sqrt_1o2, db )
      endif
    endif
  else
    m = m_r
    c_harm = ( 1._db, 0._db )
  endif

  return
end

!***********************************************************************

! Calcul de la solution singuliere pour la densite d'etat

Subroutine Cal_Solsing_sd(Classic_irreg,drho,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Harm_cubic,Hubb_a,Hubb_d, &
          icheck,konde,LL,lmax,m_hubb,nlm,nlm_pot,nlm1,nlm2,nlma,nr,nrmtsd,nrs,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Self, &
          Spinorbite,State,Tau,ui,ur,V,V_hubb,Comp_to_real)

  use declarations
  implicit none

  integer:: i_k, icheck, ir, isg, iso, isp, iss, L, LL, lm, lmax, lp, m, m_hubb, m_k, m_r, mp, ms, mv, n_k, n_rk, nlm, &
    nlm_pot, nlm1, nlm2, n, n_r, n0, nlma, np, nr, nrmax, nrmtsd, nrs, nspin, nspino, nspinp, numat

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlm1,nspinp,nlm1,nspinp):: Tau
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  logical:: Classic_irreg, Ecomp, Full_potential, Harm_cubic, Hubb_a, Hubb_d, Radial_comp, Self, Spinorbite, Comp_to_real

  real(kind=db):: c_cubic, c_harm, Eimag, f_integr3, Rmtsd, Sing
  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nrmtsd):: fctt, ri2, rr
  real(kind=db), dimension(nrmtsd,nspinp):: fct
  real(kind=db), dimension(nr,nspin):: g0, gm, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur
  real(kind=db), dimension(nrs,nlm_pot,nspinp):: drho
  real(kind=db), dimension(nr,nlm,nlm,nspin):: V
  real(kind=db), dimension(nlma,nspinp,nlma,nspinp):: State
  
  real(kind=db), dimension(:,:,:,:,:), allocatable:: usi, usr

  rr(1:nrmtsd) = r(1:nrmtsd)

  nrmax = nrmtsd + 1
  allocate( usi(nrmax,nlm1,nlm2,nspinp,nspino) )
  allocate( usr(nrmax,nlm1,nlm2,nspinp,nspino) )

! gm and gp reverse in subroutine
  call Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde, &
       LL,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtsd,nrmax,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Spinorbite,Tau,usi,usr,V,V_hubb)

  if( icheck == 2 ) write(3,110)

  ri2(1:nrmtsd) = 1 / ( quatre_pi * r(1:nrmtsd)**2 )

! Here one takes r = r', because one calculates the trace of G(r,r')
! One calculates only the imaginary part
  n_k = 0
  do L = 0,lmax
    if( .not. Full_potential .and. L /= LL ) cycle
    n0 = L**2 + L + 1

    do m_k = -L,L
      if( nlm1 == 1 .and. m_k /= 0 ) cycle
      n_k = n_k + 1

! When nlm = 1, all radial wave function of same l are equal and the basis change for harm cubic is not useful
      if( nlm == 1 .or. .not. Harm_cubic .or. L /= 3 .or. ( m_k == -3 .or. m_k == 0 .or. m_k == 3 ) ) then
        n_rk = 1
      else
        n_rk = 2
      endif

      fct(:,:) = 0._db

      do isp = 1,nspinp

        do i_k = 1,n_rk

          if( nlm1 == 1 ) then
            m_r = m_k
            c_cubic = 1._db
          else
            call Trans_TtoK(c_cubic,Harm_cubic,i_k,L,m_k,m_r)
          endif
          
          n_r = n_k + m_r - m_k
          
          do isg = -1,1,2
            if( Comp_to_real .and. m_r /= 0 ) then
              c_harm = 0.5_db
              m = isg * m_r
              if( isg == 1 ) then
                n = n_r
              else
                n = n_r - 2 * m_r
              endif
            else
              if( isg == 1 ) cycle
              c_harm = 1._db
              m = m_r
              n = n_r
            endif

            c_harm = c_harm * c_cubic**2
            
            np = 0
            do lp = 0,lmax
              if( .not. Full_potential .and. lp /= L ) cycle
              do mp = -lp,lp
                if( nlm2 == 1 .and. mp /= m ) cycle
                np = np + 1
      
                do iso = 1,nspino
                  if( Spinorbite .and. nlm2 == 1 ) then
                    ms = m + iso - isp
                    if( ms > L .or. ms < -L ) cycle
                  endif
  
                  if( Spinorbite .and. nlm2 > 1 ) then
                    iss = iso
                  else
                    iss = isp
                  endif
                      
                  do ir = 1,nrmtsd
                    if( Radial_comp ) then
                      fct(ir,iss) = fct(ir,iss) + c_harm * r(ir) * ( ur(ir,n,np,isp,iso) * usi(ir,n,np,isp,iso) &
                                                           + ui(ir,n,np,isp,iso) * usr(ir,n,np,isp,iso) )
                    else
                      fct(ir,iss) = fct(ir,iss) + c_harm * r(ir) * ur(ir,n,np,isp,iso) * usi(ir,n,np,isp,iso)
                    endif
                  end do
      
                end do
              end do
            end do
          end do

        end do

      end do

      do iss = 1,nspinp

        fctt(:) = fct(:,iss)
               
        if( icheck > 2 .and. Self ) then
          write(3,115) L, m, iss
          do ir = 1,nrmtsd
            write(3,180) r(ir)*bohr, fctt(ir)
          end do
        endif
        
        Sing = f_integr3(rr,fctt,1,nrmtsd,Rmtsd)
        if( Self ) fctt(1:nrmtsd) = fctt(1:nrmtsd) * ri2(1:nrmtsd)

        if( nlm1 == 1 ) then
          do mv = -L,L
            lm = n0 + mv
            State(lm,iss,lm,iss) = State(lm,iss,lm,iss) + Sing
            if( Self ) drho(1:nrmtsd,1,iss) = drho(1:nrmtsd,1,iss) + fctt(1:nrmtsd)
            if( icheck > 2 ) write(3,110)
            if( icheck > 1 ) write(3,120) L, m, iss, Sing, State(lm,iss,lm,iss)
          end do
        else
          lm = n0 + m_k
          State(lm,iss,lm,iss) = State(lm,iss,lm,iss) + Sing
          if( Self ) drho(1:nrmtsd,1,iss) = drho(1:nrmtsd,1,iss) + fctt(1:nrmtsd)
          if( icheck > 2 ) write(3,110)
          if( icheck > 1 ) write(3,120) L, m, iss, Sing, State(lm,iss,lm,iss)
        endif

      end do
    end do
  end do

  deallocate( usi, usr )
  
  if( icheck > 2 .and. Self ) then
    write(3,170)
    do ir = 1,nrmtsd
      write(3,180) r(ir)*bohr, drho(ir,1,:)
    end do
  endif

  return
  110 format(/' Radial integral of singular solution for the density of state:',/'  L  m isp      Sing         State')
  115 format(/' Singular radial function for the density of state, L, m, isp =',3i3)
  120 format(3i3,1p,1x,e13.5,2x,2e13.5)
  170 format(/' drho after singul', /'   Radius         up           dn')
  180 format(1p,3e13.5)
end

!**********************************************************************

! Calculation of the multiple scttaering amplitude in the complex energy case with FDM
! Under construction

subroutine Cal_Tau_comp(Cal_xanes,Ecinetic,Eimag,Energ,Enervide,Full_atom, &
                Full_potential,Hubb,Hubb_diag,Hubb_diag_abs,iaabsi,iaprotoi,icheck,ie,itypei,itypepr,lmax_comp,lmax_pot, &
                lmaxat,m_hubb,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nenerg,nlm_comp,nlmagm, &
                nlm_pot,nrato,nrm,nspin,nspino,nspinp,ntype,numat,RadialIntegral,rato,Relativiste,Renorm,Rmtg,Rmtsd, &
                Spinorbite,State_all,Taull,Tau_comp,V_hubb,V_hubb_abs,V_intmax,V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer:: ia, iaabsi, iapr, icheck, ie, ipr, ir, isp1, isp2, it, je, L, lm1, lm2, lmax, lmax_comp, lmax_pot, &
    m_hubb, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, &
    n_atom_proto, natome, nenerg, nlm_comp, nlm_pot, nlma, nlma2, nlmagm, nr, nrm, nspin, nspino, nspinp, ntype, Z

  complex(kind=db), dimension(nenerg):: Int_step
  complex(kind=db), dimension(nlm_comp,nspinp,nlm_comp,nspinp,nenerg):: Tau_comp
  complex(kind=db), dimension(nlm_comp,1,nspinp,nlm_comp,1,nspinp,nenerg):: RadialIntegral
  complex(kind=db), dimension(nlm_comp,1,nspinp,nlm_comp,1,nspinp):: RadialIntegral_c
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,natome):: Taull
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb_abs
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb
  complex(kind=db), dimension(:,:,:,:), allocatable:: V_hubb_t
  complex(kind=db), dimension(:,:,:,:), allocatable:: Taulla

  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, lmaxat
  integer, dimension(natome):: iaprotoi, itypei

  logical:: Absorbeur, Cal_xanes, Full_atom, Full_potential, Hubb_a, Hubb_d, Hubb_diag_abs, Relativiste, &
    Renorm, Spinorbite, State_all, Ylm_comp
  logical, dimension(0:ntype):: Hubb
  logical, dimension(n_atom_0:n_atom_ind):: iapr_done
  logical, dimension(n_atom_0_self:n_atom_ind_self):: Hubb_diag

  real(kind=db):: a, b, dE, E1, E2, Eimag, Enervide, V_intmax
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(0:lmax_comp):: pop
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(0:n_atom_proto):: Rmtg, Rmtsd
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vrato
  real(kind=db), dimension(:), allocatable:: r
  real(kind=db), dimension(:,:,:), allocatable:: Vrato_t

  iapr_done(:) = .false.
  
  RadialIntegral_c(:,:,:,:,:,:) = (0._db,0._db)

  boucle_ia: do ia = 1,natome

    Absorbeur = ( ia == iaabsi .and. Cal_xanes )
    if( Cal_xanes .and. .not. ( State_all .or. Absorbeur ) ) cycle

    ipr = iaprotoi(ia)
    it = itypepr(ipr)
    Hubb_a = Hubb(it)

    Z = numat( it )
    if( Z < 19 ) then
      lmax = min(2,lmaxat(ipr))
    elseif( Z > 18 .and. Z < 55 ) then
      lmax = min(3,lmaxat(ipr))
    else
      lmax = min(4,lmaxat(ipr))
    endif
    lmax = min( lmax,lmax_comp)
    nlma = (lmax + 1)**2
    if( Full_potential ) then
      nlma2 = nlma
    else
      nlma2 = 1
    endif

    if( Full_atom ) then
      iapr = ia
    else
      iapr = ipr
      if( iapr_done(iapr) ) cycle
    endif
    it = itypei(ia)

    do ir = 1,nrato(it)
      if( rato(ir,it) > Rmtsd(ipr) + eps10 ) exit
    end do
    nr = ir + 1

    allocate( r(nr) )
    allocate( Vrato_t(nr,nlm_pot,nspin) )
    allocate( V_hubb_t(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp) )
    allocate( Taulla(nlma,nspinp,nlma,nspinp) )
    r(1:nr) = rato(1:nr,it)
    Vrato_t(1:nr,:,:) = Vrato(1:nr,:,:,iapr)
    if( Hubb(it) .and. iapr <= n_atom_ind_self ) then
      Hubb_a = .true.
      if( Absorbeur ) then
        V_Hubb_t(:,:,:,:) = V_Hubb_abs(:,:,:,:)
        Hubb_d = Hubb_diag_abs
      else
        V_Hubb_t(:,:,:,:) = V_Hubb(:,:,:,:,iapr)
        Hubb_d = Hubb_diag(iapr)
      endif
    else
      Hubb_a = .false.
      Hubb_d = .true.
    end if

    call radial_sd_comp(Ecinetic,Eimag,Enervide,Full_potential,Hubb_a,Hubb_d,icheck,ie,lmax, &
         lmax_pot,m_hubb,nenerg,nlm_comp,nlm_pot,nlma2,nr,nspin,nspino,nspinp,Z,r,RadialIntegral,RadialIntegral_c, &
         Relativiste,Renorm,Rmtg(ipr),Rmtsd(ipr),Spinorbite,V_hubb_t,V_intmax,V0bd,Vrato_t,Ylm_comp)

    deallocate( r, Taulla, V_hubb_t, Vrato_t )

!!!!!!!!!!!!!!!!!!!!!!
! Cas sans spinorbite pour L'instant, ni Full_poential, ni Hubbard
!!!!!!!!!!!!!!!!!!!!!!
   if( nenerg == 1 ) then
     dE = 1._db
   elseif( ie == 1 ) then
     dE = Energ(ie+1) - Energ(ie)
   elseif( ie == nenerg ) then
     dE = Energ(ie) - Energ(ie-1)
   else
     dE = ( Energ(ie+1) - Energ(ie-1) ) / 2
   endif

    do je = 1,nenerg

      if( ie == 1 ) then
        E1 = Energ(ie) - 3 * Eimag
      else
        E1 = ( Energ(ie) + Energ(ie-1) ) / 2
      endif
 
      if( ie == nenerg ) then
        E2 = Energ(ie) + 3 * Eimag
      else
        E2 = ( Energ(ie+1) + Energ(ie) ) / 2
      endif
 !     Int_step = ( img / pi ) * log( ( Energ(je) - E1 + img*Eimag/2 ) / ( Energ(je) - E2 + img*Eimag/2 ) )
      a = 0.5_db * log( ( (Energ(je) - E1)**2 + (Eimag/2)**2 ) / ( (Energ(je) - E2)**2 + (Eimag/2)**2 ) )
      b = atan( ( E1 - Energ(je) ) / (Eimag/2) ) - atan( ( E2 - Energ(je) ) / (Eimag/2) ) 
      Int_step(je) = ( img / pi ) * cmplx(a,b,db)
         
      do lm1 = 1,nlm_comp
        do isp1 = 1,nspinp
          do lm2 = 1,nlm_comp
            do isp2 = 1,nspinp
 
              Tau_comp(lm1,isp1,lm2,isp2,je) = Tau_comp(lm1,isp1,lm2,isp2,je) &
                                             + Int_step(je) * Taull(lm1,isp1,lm2,isp2,ia) * RadialIntegral_c(lm1,1,isp1,lm2,1,isp2)            
              
            end do
          end do      
        end do
      end do
      
    end do

    iapr_done(iapr) = .true.

  end do boucle_ia
  
  if( ie < nenerg .or. icheck < 1 ) return
 
  pop(:) = 0._db 
  write(3,110) ( L, L, L, L = 0,lmax_comp )  
  do je = 1,nenerg
      if( je == 1 ) then
        E1 = Energ(je) - 3 * Eimag
      else
        E1 = ( Energ(je) + Energ(je-1) ) / 2
      endif
 
      if( je == nenerg ) then
        E2 = Energ(je) + 3 * Eimag
      else
        E2 = ( Energ(je+1) + Energ(je) ) / 2
      endif

    do L = 0,lmax_comp 
      pop(L) = pop(L) - aimag( Tau_comp(L**2+L+1,1,L**2+L+1,1,je) ) * ( E2 - E1 )
    end do  
    write(3,120) Energ(je)*rydb, ( Tau_comp(L**2+L+1,1,L**2+L+1,1,je), pop(L),  L = 0,lmax_comp )  
  end do

  return
  110 format(/'      Energy ', 10(3x,'Tau_c_r(',i1,')',3x,'Tau_c_i(',i1,')',5x,'pop(',i1,')',4x) )
  120 format(f13.5,10(1x,2e13.5,1x,e13.5))
end

!**********************************************************************

! Calculation of radial integral in energy complex case for FDM
! Under construction

subroutine radial_sd_comp(Ecinetic,Eimag,Enervide,Full_potential,Hubb_a,Hubb_d,icheck,ie,lmax, &
         lmax_pot,m_hubb,nenerg,nlm_comp,nlm_pot,nlma2,nr,nspin,nspino,nspinp,numat,r,RadialIntegral,RadialIntegral_c, &
         Relativiste,Renorm,Rmtg,Rmtsd,Spinorbite,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, ie, ir, is1, is2, iso1, iso2, isp, isp1, isp2, L, l_hubbard, L2, lfin, lm, lm0, lm1, lm2, lmp, lmp1, lmp2, &
    lmax, lmax_pot, m, m_hubb, nenerg, nlm, nlm_comp, nlm_pot, nlm1, nlm2, nlma2, nr, &
    nrmtg, nrmtsd, nspin, nspino, nspinp, numat

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlm_comp,1,nspinp*nspino,nlm_comp,1,nspinp*nspino,nenerg):: RadialIntegral 
  complex(kind=db), dimension(nlm_comp,1,nspinp*nspino,nlm_comp,1,nspinp*nspino):: RadialIntegral_c 
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau

  logical:: Diagonale, Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m, Radial_comp, Relativiste, Renorm, &
    Spinorbite, Ylm_comp

  real(kind=db):: Eimag, Enervide, f_integr3, Hpsi_r, Hpsi_i, pHp_i, pHp_r, pp, pp_i, pp_r, radli, radlr, Rmtg, &
                  Rmtsd, td, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nr,nspin):: g0, gm, gp, gpi
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato
  real(kind=db), dimension(nr,nlm_comp,nlma2,nspinp,nspino):: uit, urt
  real(kind=db), dimension(:), allocatable:: dr, fci, fcr, r2, rr
  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur

  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

  uit(:,:,:,:,:) = 0._db
  urt(:,:,:,:,:) = 0._db

  if( Full_potential ) then
    nlm = ( lmax + 1 )**2
  else
    nlm = 1
  endif
  allocate( V(nr,nlm,nlm,nspin) )

  call mod_V(icheck,lmax,lmax_pot,nlm,nlm_pot,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  allocate( dr(nrmtsd) )
  allocate( rr(nrmtsd) )
  allocate( r2(nrmtsd) )
  allocate( fci(nrmtsd) )
  allocate( fcr(nrmtsd) )
  
  rr(1:nrmtsd) = r(1:nrmtsd)
  r2(1:nrmtsd) = r(1:nrmtsd)**2

  dr(1) = 0.5_db * ( r(2) + r(1) )
  do ir = 2,nrmtsd-1
    dr(ir) = 0.5_db * ( r(ir+1) - r(ir-1) )
  end do
  dr(nrmtsd) = Rmtsd - 0.5_db * ( r(nrmtsd-2) + r(nrmtsd-1) )

  call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)

  gpi(:,:) = 1 / gp(:,:)

  if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10 .or. Ecinetic(nspin) < eps10 ) then
    Ecomp = .true.
  else
    Ecomp = .false.
  endif
  Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

  if( Full_potential ) then
    lfin = 0
  else
    lfin = lmax
  endif

  do L = 0,lfin

    if( Hubb_a .and. L == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = nlm ! = ( lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*L + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*L + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif
    Diagonale = .not. ( Full_potential .or. Hubb_m )

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,L,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

! Energy calculation

    L2 = L*(L+1)

    php_r = 0._db
    php_i = 0._db
    pp_r = 0._db
    pp_i = 0._db
    isp = 1

    do ir = 1,nrmtsd
! one must add E because g0 contains V - E
      td = g0(ir,isp) + L2 * f2(ir) + Enervide
      if( ir == 1 ) then
        if( L == 0 ) then 
          Hpsi_r = gm(ir,isp) * ur(ir,1,1,1,1) +  td * ur(ir,1,1,1,1) + gp(ir,isp) * ur(ir+1,1,1,1,1)
          Hpsi_i = gm(ir,isp) * ui(ir,1,1,1,1) +  td * ui(ir,1,1,1,1) + gp(ir,isp) * ui(ir+1,1,1,1,1)
        else
          Hpsi_r = td * ur(ir,1,1,1,1) + gp(ir,isp) * ur(ir+1,1,1,1,1)
          Hpsi_i = td * ui(ir,1,1,1,1) + gp(ir,isp) * ui(ir+1,1,1,1,1)
        endif
      else
        Hpsi_r = gm(ir,isp) * ur(ir-1,1,1,1,1) +  td * ur(ir,1,1,1,1) + gp(ir,isp) * ur(ir+1,1,1,1,1)
        Hpsi_i = gm(ir,isp) * ui(ir-1,1,1,1,1) +  td * ui(ir,1,1,1,1) + gp(ir,isp) * ui(ir+1,1,1,1,1)
      endif
      php_r = php_r + ( ur(ir,1,1,1,1) * Hpsi_r + ui(ir,1,1,1,1) * Hpsi_i) * dr(ir) * r2(ir)
      php_i = php_i + ( ur(ir,1,1,1,1) * Hpsi_i - ui(ir,1,1,1,1) * Hpsi_r) * dr(ir) * r2(ir)
      fcr(ir) = ( ur(ir,1,1,1,1) * Hpsi_r + ui(ir,1,1,1,1) * Hpsi_i ) * r2(ir)
      fci(ir) = ( ur(ir,1,1,1,1) * Hpsi_i - ui(ir,1,1,1,1) * Hpsi_r ) * r2(ir)

      pp_r = pp_r + ur(ir,1,1,1,1)**2 * dr(ir) * r2(ir)
      pp_i = pp_i + ui(ir,1,1,1,1)**2 * dr(ir) * r2(ir)
    end do
    pp = pp_r + pp_i

    php_r = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
    php_i = f_integr3(rr,fci,1,nrmtsd,Rmtsd)

    do ir = 1,nrmtsd
      fcr(ir) = ( ur(ir,1,1,1,1)**2 + ui(ir,1,1,1,1)**2 ) * r2(ir)
    end do
    
    pp = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)

    php_r = ( php_r / pp ) * Rydb
    php_i = ( php_i / pp ) * Rydb
    
    do lm1 = 1,nlm1
      do lmp1 = 1,nlm2
        do isp1 = 1,nspinp
          do iso1 = 1,nspino

            if( Full_potential ) then
              urt(1:nrmtsd,lm1,lmp1,isp1,iso1) = ur(1:nrmtsd,lm1,lmp1,isp1,iso1)
              uit(1:nrmtsd,lm1,lmp1,isp1,iso1) = ui(1:nrmtsd,lm1,lmp1,isp1,iso1)
            elseif( Hubb_m .and. .not. Hubb_d ) then
              lm0 = L**2
              urt(1:nrmtsd,lm0+lm1,lm0+lmp1,isp1,iso1) = ur(1:nrmtsd,lm1,lmp1,isp1,iso1)
              uit(1:nrmtsd,lm0+lm1,lm0+lmp1,isp1,iso1) = ui(1:nrmtsd,lm1,lmp1,isp1,iso1)
            elseif( Spinorbite .or. Hubb_m ) then
              lm0 = L**2
              urt(1:nrmtsd,lm0+lm1,lmp1,isp1,iso1) = ur(1:nrmtsd,lm1,lmp1,isp1,iso1)
              uit(1:nrmtsd,lm0+lm1,lmp1,isp1,iso1) = ui(1:nrmtsd,lm1,lmp1,isp1,iso1)
            else
              lm0 = L**2 + L + 1
              do m = -L,L
                urt(1:nrmtsd,lm0+m,lmp1,isp1,iso1) = ur(1:nrmtsd,lm1,lmp1,isp1,iso1)
                uit(1:nrmtsd,lm0+m,lmp1,isp1,iso1) = ui(1:nrmtsd,lm1,lmp1,isp1,iso1)
              end do
            endif
                  
          end do
        end do
      end do
    end do

    deallocate( Tau, ui, ur )
    
  end do   ! end of loop over L

  deallocate( V )

! Radial integral for the DOS

  lm = (lmax + 1)**2

  if( .not. ( Full_potential .or. ( Hubb_a .and. .not. Hubb_d ) ) ) then
    lmp = 1
  else
    lmp = lm
  endif
  
  do lm1 = 1,lm
    do lmp1 = 1,lmp
      is1 = 0
      do isp1 = 1,nspinp
        do iso1 = 1,nspino
          is1 = is1 + 1 

          do lm2 = 1,lm
            do lmp2 = 1,lmp
              is2 = 0
              do isp2 = 1,nspinp
                do iso2 = 1,nspino
                  is2 = is2 + 1 
                  if( .not. Spinorbite .and. isp1 /= isp2 ) cycle
! u * up
                  fcr(1:nrmtsd) = r2(1:nrmtsd) * ( urt(1:nrmtsd,lm1,lmp1,isp1,iso1) * urt(1:nrmtsd,lm2,lmp2,isp2,iso2) &
                                                 - uit(1:nrmtsd,lm1,lmp1,isp1,iso1) * uit(1:nrmtsd,lm2,lmp2,isp2,iso2) )
                  radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
                  fci(1:nrmtsd) = r2(1:nrmtsd) * ( urt(1:nrmtsd,lm1,lmp1,isp1,iso1) * uit(1:nrmtsd,lm2,lmp2,isp2,iso2) &
                                                 + uit(1:nrmtsd,lm1,lmp1,isp1,iso1) * urt(1:nrmtsd,lm2,lmp2,isp2,iso2) )
                  radli = f_integr3(rr,fci,1,nrmtsd,Rmtsd)

                  RadialIntegral(lm1,lmp1,is1,lm2,lmp2,is2,ie) = cmplx( radlr, radli, db )
! u * conjg( up )                
                  fcr(1:nrmtsd) = r2(1:nrmtsd) * ( urt(1:nrmtsd,lm1,lmp1,isp1,iso1) * urt(1:nrmtsd,lm2,lmp2,isp2,iso2) &
                                                 + uit(1:nrmtsd,lm1,lmp1,isp1,iso1) * uit(1:nrmtsd,lm2,lmp2,isp2,iso2) )
                  radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
                  fci(1:nrmtsd) = r2(1:nrmtsd) * ( - urt(1:nrmtsd,lm1,lmp1,isp1,iso1) * uit(1:nrmtsd,lm2,lmp2,isp2,iso2) &
                                                 + uit(1:nrmtsd,lm1,lmp1,isp1,iso1) * urt(1:nrmtsd,lm2,lmp2,isp2,iso2) )
                  radli = f_integr3(rr,fci,1,nrmtsd,Rmtsd)

                  RadialIntegral_c(lm1,lmp1,is1,lm2,lmp2,is2) = cmplx( radlr, radli, db )
                
                end do
              end do
            end do
          end do

        end do
      end do
    end do
  end do

  deallocate( dr, fci, fcr, r2, rr )
  
  return
end

!***********************************************************************

! Calculation of the Crystal Orbital Overlap Population

subroutine Cal_COOP(Coop_z_along_bond,Coverlap,Density_comp_ext,Dist_coop,Ecinetic,Eimag,Energ,Enervide, &
            Full_atom,Full_potential,Harm_cubic_ext,Hubb,Hubb_diag,Hubb_diag_abs,ia_coop,iaabsi,iaprotoi,icheck,ie,index_coop, &
            iprabs,itypepr,lmax_pot,lmaxat,m_hubb,mpinodee,mpirank,n_atom_0,n_atom_0_self, &
            n_atom_coop,n_atom_ind,n_atom_ind_self,n_atom_proto,nab_coop,natome,nlm_pot,nlmagm,nrato,nrm, &
            nspin,nspino,nspinp,ntype,numat,Posi,rato,Relativiste,Renorm,Rmtg,Rmtg0,Rmtsd,Rot_atom,Spinorbite, &
            Tau_coop,V_hubb, &
            V_hubb_abs,V_intmax,V0bd,Vrato,Ylm_comp,Z_nospinorbite)

  use declarations
  implicit none

  integer:: i, i_rad, i_tr, ia, iaabsi, iab, iapr, iapra, iaprb, ib, icheck, ie, ipr, ipra, iprabs, iprb, isa, isb, &
    it, j, la, lb, lm, lma, lmax, lmax_coop, lmaxa, lmaxb, lmaxm, lmax_pot, &
    lmb, lmm, lmp, m_hubb, ma, mb, mpa, mpb, mpinodee, mpirank, &
    n, n_angle, n_atom_0, n_atom_0_self, n_atom_coop, n_atom_ind, n_atom_ind_self, n_atom_proto, n_pt_integr, &
    n_rad, n_rad_m, nab_coop, natome, nlm_pot, nlmagm, nlmc, nlmr, nr, nrm, nspin, nspino, nspinp, &
    ntype, Z, Z_nospinorbite, Za, Zb

  integer, dimension(n_atom_coop):: ia_coop
  integer, dimension(natome):: iaprotoi
  integer, dimension(n_atom_0:n_atom_ind):: n_radius
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, lmaxat
  integer, dimension(2,nab_coop):: index_coop
  integer, dimension(:), allocatable:: index_r_a, index_r_b

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb_abs
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,nspino):: Tau_coop_r
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,nab_coop,nspino):: Tau_coop
  complex(kind=db), dimension(:), allocatable:: Ylmc
  complex(kind=db), dimension(:,:,:), allocatable:: Dlmm_ia, Dlmm_ib
  complex(kind=db), dimension(:,:,:,:), allocatable:: Taull, V_hubb_t
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: um
  
  logical:: Absorbeur, Comp_to_real, Comp_to_real_out, Coop_z_along_bond, Density_comp_ext, Density_comp, &
            Density_real, &
            Full_atom, Full_potential, Harm_cubic, Harm_cubic_ext, Hubb_a, Hubb_d, Hubb_diag_abs, &
            Relativiste, Renorm, Spino, Rotation_bound, Spinorbite, Ylm_comp, Ylm_comp_out 
  logical, dimension(0:ntype):: Hubb
  logical, dimension(n_atom_0_self:n_atom_ind_self):: Hubb_diag

  real(kind=db):: Angle, Dist, Eimag, Energ, Enervide, R_circle, Rad_max, Radius_coop, Step, Step_rad, &
                  V_intmax, Vol
  real(kind=db), dimension(2):: Dist_coop
  real(kind=db), dimension(3):: V, Vx, Vy, Vz, Wx, Wy, Wz
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(0:n_atom_proto):: Rmtg, Rmtg0, Rmtsd
  real(kind=db), dimension(3,3):: Mat_rot, Rot_a, Rot_b
  real(kind=db), dimension(3,3,natome):: rot_atom
  real(kind=db), dimension(3,natome):: Posi
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vrato 
  real(kind=db), dimension(16,16,nspinp,nab_coop,0:mpinodee-1):: Coverlap
  real(kind=db), dimension(16,16,nspinp):: State

  real(kind=db), dimension(:), allocatable:: r, Ylmr
  real(kind=db), dimension(:,:), allocatable:: Radius_urm, Ylm_a, Ylm_b 
  real(kind=db), dimension(:,:,:), allocatable:: Vrato_e

  n_rad = 1
  Rad_max = 0.20_db
  Step_rad = Rad_max / max( n_rad, 1 )
  n_angle = 8
  Step = 2 * pi / max( n_angle, 1 )
 
! The following is because the real coop does not work when spinorbit...  
  if( Spinorbite ) then
    Density_comp = .true.
    Harm_cubic = .false.
  else
    Harm_cubic = Harm_cubic_ext
    Density_comp = Density_comp_ext
  endif

! Calculation of the radius where the radial functions must be calculated
  do i = 1,2
  
    n_radius(:) = 0
    
    do ia = 1,natome

     ipra = iaprotoi( ia )
      if( Full_atom ) then
        iapra = ia
      else
        iapra = ipra
      endif

      do ib = ia+1,natome

        if( n_atom_coop > 0 ) then
          do j = 1,n_atom_coop
            if( ia_coop(j) == 0 ) cycle
            if( ia == ia_coop(j) .or. ib == ia_coop(j) ) exit
          end do
          if( j > n_atom_coop ) cycle
        endif

        iprb = iaprotoi( ib )
        if( Full_atom ) then
          iaprb = ib
        else
          iaprb = iprb
        endif

        Dist = sqrt( sum( ( Posi(:,ia) - Posi(:,ib) )**2 ) )
        if( Dist < Dist_coop(1) - eps10 ) cycle
        if( Dist_coop(2) > eps10 ) then
          if( Dist > Dist_coop(2) + eps10 ) cycle
        else
          if( Dist > Rmtg( ipra ) + Rmtg( iprb ) ) cycle
        endif

        Radius_coop = ( Rmtg0( ipra ) + Rmtg0( iprb ) ) / 2 

        do i_rad = 0,n_rad
          n_radius(iapra) = n_radius(iapra) + 1
          n_radius(iaprb) = n_radius(iaprb) + 1
          if( i == 2 ) then
            if( i_rad == 0 ) then 
              Radius_urm(n_radius(iapra),iapra) = Rmtg0(ipra) 
              Radius_urm(n_radius(iaprb),iaprb) = Rmtg0(iprb)
            else 
              Radius_urm(n_radius(iapra),iapra) = sqrt( Rmtg0(ipra)**2 + ( i_rad * Step_rad * Radius_coop )**2 ) 
              Radius_urm(n_radius(iaprb),iaprb) = sqrt( Rmtg0(iprb)**2 + ( i_rad * Step_rad * Radius_coop )**2 )
            endif
          endif
        end do 
      end do
    end do
   
    if( i == 2 ) exit
    
    n_rad_m = 0
    do iapr = n_atom_0,n_atom_ind
      n_rad_m = max( n_rad_m, n_radius(iapr) )
    end do
    allocate( Radius_urm(n_rad_m,n_atom_0:n_atom_ind) )
    
  end do
  
  lmaxm = 0
  do ipr = 1,n_atom_proto
    lmm = min( lmaxat(ipr), lmax_coop( numat( itypepr(ipr) ) ) )
    lmaxm = max( lmaxm, lmm )
  end do
  lmm = (lmaxm + 1)**2

  allocate( um(n_rad_m,lmm,nspinp,nspino,n_atom_0:n_atom_ind) )
  um(:,:,:,:,:) = (0._db, 0._db)
    
! Loop over non equivalent absorbing atoms

  do iapr = n_atom_0,n_atom_ind

    if( Full_atom ) then
      ipr = iaprotoi( iapr )
    else
      ipr = iapr
    endif
    it = itypepr(ipr)
    Z = numat( it )
    if( Z == 0 ) cycle
    if( Z == Z_nospinorbite ) then
      Spino = .false.
    else
      Spino = Spinorbite
    endif

    allocate( V_hubb_t(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp))
    if( Hubb(it) .and. iapr <= n_atom_ind_self ) then
      Hubb_a = .true.
      if( Full_atom ) then
         Absorbeur = iapr == iaabsi
      else
         Absorbeur = iapr == iprabs
      endif
      if( Absorbeur ) then
        V_Hubb_t(:,:,:,:) = V_Hubb_abs(:,:,:,:)
        Hubb_d = Hubb_diag_abs
      else
        V_Hubb_t(:,:,:,:) = V_Hubb(:,:,:,:,iapr)
        Hubb_d = Hubb_diag(iapr)
      endif
    else
      Hubb_a = .false.
      Hubb_d = .true.
    end if

    nr = nrato(it)
    allocate( r(nr) )
    allocate( Vrato_e(nr,nlm_pot,nspin) )
    r(1:nr) = rato(1:nr,it)

    lmax = min( lmaxat(ipr), lmax_coop( Z ) )

    Vrato_e(1:nr,:,:) = Vrato(1:nr,:,:,iapr)

! In the following routine, when Hubbard, only the diagonal radial function is kept.
    call Cal_urm(Ecinetic,Eimag,Energ,Enervide, &
            Full_potential,Hubb_a,Hubb_d,iapr,icheck,lmax,lmax_pot,lmm,m_hubb,n_atom_0,n_atom_ind,n_rad_m,n_radius,nlm_pot, &
            nr,nspin,nspino,nspinp,r,Radius_urm,Relativiste,Renorm,Rmtg(ipr),Rmtsd(ipr),Spinorbite,um,V_hubb, &
            V_intmax,V0bd,Vrato_e,Ylm_comp,Z)

    deallocate( r )
    deallocate( Vrato_e, V_hubb_t )

  end do

  COverlap(:,:,:,:,:) = 0._db
  
  iab = 0
  n_radius(:) = 0

  Density_real = .not. Density_comp
  Comp_to_real = Ylm_comp .and. Density_real 

  do ia = 1,natome

    ipra = iaprotoi( ia )
    Za = numat( itypepr(ipra) )
    if( Full_atom ) then
      iapra = ia
    else
      iapra = ipra
    endif
    lmaxa = min( lmaxat( ipra ), lmax_coop( Za ) )
    allocate( Dlmm_ia(-lmaxa:lmaxa,-lmaxa:lmaxa,0:lmaxa) )

    do ib = ia+1,natome

      if( n_atom_coop > 0 ) then
        do j = 1,n_atom_coop
          if( ia_coop(j) == 0 ) cycle
          if( ia == ia_coop(j) .or. ib == ia_coop(j) ) exit
        end do
        if( j > n_atom_coop ) cycle
      endif

      iprb = iaprotoi( ib )
      if( Full_atom ) then
        iaprb = ib
      else
        iaprb = iprb
      endif

      Dist = sqrt( sum( ( Posi(:,ia) - Posi(:,ib) )**2 ) )
      if( Dist < Dist_coop(1) - eps10 ) cycle
      if( Dist_coop(2) > eps10 ) then
        if( Dist > Dist_coop(2) + eps10 ) cycle
      else
        if( Dist > Rmtg( ipra ) + Rmtg( iprb ) ) cycle
      endif

      iab = iab + 1

      Zb = numat( itypepr(iprb) )
      lmaxb = min( lmaxat( iprb ), lmax_coop( Zb ) )
      allocate( Dlmm_ib(-lmaxb:lmaxb,-lmaxb:lmaxb,0:lmaxb) )

      index_coop(1,iab) = ia
      index_coop(2,iab) = ib

      Radius_coop = ( Rmtg0( ipra ) + Rmtg0( iprb ) ) / 2 
      Vol = ( 4 * pi / 3._db ) * Radius_coop**3
     
      Rot_a(:,:) = rot_atom(:,:,ia)
      Rot_b(:,:) = rot_atom(:,:,ib)
      Rot_a = Transpose( Rot_a )
      Rot_b = Transpose( Rot_b )
      
! Local basis with z axis along the bond direction
      Vz(:) = Posi(:,ib) - Posi(:,ia)
      Vz(:) = Vz(:) / sqrt( sum( Vz(:)**2 ) )
            
      if( abs( Vz(3) - 1._db ) > eps10 ) then

        Vy(1) = - Vz(2) / sqrt( Vz(1)**2 + Vz(2)**2 ) 
        Vy(2) = Vz(1) / sqrt( Vz(1)**2 + Vz(2)**2 )
        Vy(3) = 0._db
              
        call prodvec( Vx, Vy, Vz )
        Vx(:) = Vx(:) / sqrt( sum( Vx(:)**2 ) )

      else
 
        Vx(1) = 1._db; Vx(2) = 0._db; Vx(3) = 0._db 
        Vy(1) = 0._db; Vy(2) = 1._db; Vy(3) = 0._db 
 
      endif

      if( Coop_z_along_bond .and. abs( Vz(3) - 1._db ) > eps10 ) then

        Mat_rot(:,1) = Vx(:)  
        Mat_rot(:,2) = Vy(:) 
        Mat_rot(:,3) = Vz(:)              

! It is the inverse rotation because it is the Ylm wich moves, not the space.
        Mat_rot = Transpose( Mat_rot )
        
        Rot_a = Matmul( Mat_rot, Rot_a )
        Rot_b = Matmul( Mat_rot, Rot_b )
        
      endif

      if( icheck > 1 .and. ie == 1 ) then
        write(3,'(/3(a6,i3),a7,f10.5,a4)') ' iab =',iab,', ia =',ia,', ib =', ib,', Vol =', Vol*bohr**3, ' A^3'
        write(3,'(/a44,3(/17x,3f10.5))') ' Local basis:     Vx        Vy        Vz', ( Vx(i), Vy(i), Vz(i), i = 1,3 )
        do i_tr = 1,nspino
          do isa = 1,nspinp
            do isb = 1,nspinp
              if( nspinp == 1 ) then
                write(3,110) iab, i_tr
              else
                write(3,120) iab, isa, isb, i_tr
              endif
              write(3,'(6x,30(10x,2i3,10x))') (( lb, mb, mb = -lb,lb), lb = 0,lmaxb )
              lma = 0 
              do la = 0,lmaxa
                do ma = -la,la
                  lma = lma + 1 
                  write(3,'(2i3,1p,100e13.5)') la, ma, (Tau_coop(lma,isa,lmb,isb,iab,i_tr), lmb = 1,(lmaxb + 1)**2)
                end do
              end do
            end do
          end do
        end do
      endif

      if( lmaxa == 10000 ) then   ! not used
!      if( Comp_to_real .and. .not. Spinorbite ) then
        lma = (lmaxa + 1)**2 
        lmb = (lmaxb + 1)**2 
        allocate( Taull(lma,nspinp,lmb,nspinp) )
        do i_tr = 1,nspino
          Taull(1:lma,:,1:lmb,:) = Tau_coop(1:lma,:,1:lmb,:,iab,i_tr)
          call Trans_Tau_ab(Comp_to_real,lmaxa,lmaxb,nspinp,Spinorbite,Taull)
          Tau_coop(1:lma,:,1:lmb,:,iab,i_tr) = Taull(1:lma,:,1:lmb,:)
          deallocate( Taull )
          if( icheck > 1 .and. ie == 1 ) then
            do isa = 1,nspinp
              do isb = 1,nspinp
                if( Comp_to_real ) then
                  write(3,'(/A)') ' Tau_coop in real Ylm basis'
                else
                  write(3,'(/A)') ' Tau_coop in complex Ylm basis'
                endif 
                lma = 0 
                do la = 0,lmaxa
                  do ma = -la,la
                    lma = lma + 1 
                    write(3,'(2i3,1p,100e13.5)') la, ma, (Tau_coop(lma,isa,lmb,isb,iab,i_tr), lmb = 1,(lmaxb + 1)**2)
                  end do
                end do
              end do
            end do  
          endif
        end do
        Comp_to_real_out = .false.
        Ylm_comp_out = .false.
      else
        Comp_to_real_out = Comp_to_real
        Ylm_comp_out = Ylm_comp
      endif
      
      Tau_coop_r(:,:,:,:,:) = (0._db, 0._db)

      Rotation_bound = &
            abs( Rot_a(1,1) - 1._db ) > eps10 .or. abs( Rot_a(2,2) - 1._db ) > eps10 .or. abs( Rot_a(3,3) - 1._db ) > eps10 &
       .or. abs( Rot_b(1,1) - 1._db ) > eps10 .or. abs( Rot_b(2,2) - 1._db ) > eps10 .or. abs( Rot_b(3,3) - 1._db ) > eps10
        
      if( Rotation_bound ) then

        call Dlmm_COOP(Dlmm_ia,icheck,lmaxa,Rot_a,Ylm_comp)
        call Dlmm_COOP(Dlmm_ib,icheck,lmaxb,Rot_b,Ylm_comp)

        do i_tr = 1,nspino
          do isa = 1,nspinp
            do isb = 1,nspinp
              if( .not. Spinorbite .and. isa /= isb ) cycle
              lma = 0
              do la = 0,lmaxa
                do ma = -la,la
                  lma = lma + 1
                  lmb = 0
                  do lb = 0,lmaxb
                    do mb = -lb,lb
                      lmb = lmb + 1
                      do mpa = -la,la
                        lm = la**2 + la + 1 + mpa
                        do mpb = -lb,lb
                          lmp = lb**2 + lb + 1 + mpb
                          if( i_tr == 1 ) then 
                            Tau_coop_r(lma,isa,lmb,isb,i_tr) = Tau_coop_r(lma,isa,lmb,isb,i_tr) + Dlmm_ia(ma,mpa,la) &
                                                           * Tau_coop(lm,isa,lmp,isb,iab,i_tr) * conjg( Dlmm_ib(mb,mpb,lb) )
                          else  ! rotation with spinorbit does not work
                            Tau_coop_r(lma,isa,lmb,isb,i_tr) = Tau_coop_r(lma,isa,lmb,isb,i_tr) + Dlmm_ia(ma,mpa,la) &
                                                           * Tau_coop(lm,isa,lmp,isb,iab,i_tr) * conjg( Dlmm_ib(mb,mpb,lb) )
                          endif  
                        end do
                      end do
                    end do  
                  end do 
                end do
              end do
            end do
          end do
        end do
        
      else
      
        Tau_coop_r(:,:,:,:,:) = Tau_coop(:,:,:,:,iab,:)
      
      endif
      
      deallocate( Dlmm_ib )        

      if( Rotation_bound .and. icheck > 1 .and. ie == 1 ) then
        do i_tr = 1,nspino
          do isa = 1,nspinp
            do isb = 1,nspinp
              if( nspinp == 1 ) then
                write(3,130) iab, i_tr
              else
                write(3,140) iab, isa, isb, i_tr
              endif 
              write(3,'(6x,30(10x,2i3,10x))') (( lb, mb, mb = -lb,lb), lb = 0,lmaxb )
              lma = 0 
              do la = 0,lmaxa
                do ma = -la,la
                  lma = lma + 1 
                  write(3,'(2i3,1p,100e13.5)') la, ma, (Tau_coop_r(lma,isa,lmb,isb,i_tr), lmb = 1,(lmaxb + 1)**2)
                end do
              end do
            end do
          end do
        end do
      endif
      
      if( Coop_z_along_bond ) then
! Points are in the rotated local basis
        Wx(1) = 1._db; Wx(2) = 0._db; Wx(3) = 0._db 
        Wy(1) = 0._db; Wy(2) = 1._db; Wy(3) = 0._db
        Wz(1) = 0._db; Wz(2) = 0._db; Wz(3) = 1._db
      else
! Points are in the cluster original basis
        Wx(:) = Vx(:)
        Wy(:) = Vy(:)
        Wz(:) = Vz(:)
      endif

      n = 1
      do i_rad = 1,n_rad
        n = n + i_rad * n_angle  
      end do
      Vol = Vol / n

      n_pt_integr = n

      nlmr = ( lmaxm + 1 )**2
      nlmc = ( ( lmaxm + 1 ) * ( lmaxm + 2 ) ) / 2      
      allocate( Ylm_a(nlmr,n_pt_integr) )
      allocate( Ylm_b(nlmr,n_pt_integr) )
      allocate( Ylmc(nlmc) )
      allocate( Ylmr(nlmr) )
      allocate( index_r_a(n_pt_integr) )
      allocate( index_r_b(n_pt_integr) )
      
      j = 0
      do i_rad = 0,n_rad
        n_radius(iapra) = n_radius(iapra) + 1
        n_radius(iaprb) = n_radius(iaprb) + 1

        R_circle = i_rad * Step_rad * Radius_coop
         
! One makes the average over i_rad * n_angle      
        do i = 1,max( i_rad * n_angle, 1 )
          j = j + 1
          index_r_a(j) = n_radius(iapra) 
          index_r_b(j) = n_radius(iaprb) 

          if( i_rad == 0 ) then
            V(:) = Rmtg0( ipra ) * Wz(:)          
          else
            Angle = i * Step / i_rad       
            V(:) = Rmtg0( ipra ) * Wz(:) + R_circle * ( cos( Angle ) * Wx(:) + sin( Angle ) * Wy(:) )
          endif
              
          Dist = sqrt( sum( V(:)**2 ) )
          call cYlm(lmaxm,V,Dist,Ylmc,nlmc)
          call Ylmcr(lmaxm,nlmc,nlmr,Ylmc,Ylmr)
          Ylm_a(1:nlmr,j) = Ylmr(1:nlmr)

          if( i_rad == 0 ) then
            V(:) = - Rmtg0( iprb ) * Wz(:)
          else  
            V(:) = - Rmtg0( iprb ) * Wz(:) + R_circle * ( cos( Angle ) * Wx(:) + sin( Angle ) * Wy(:) ) 
          endif
 
          Dist = sqrt( sum( V(:)**2 ) )
          call cYlm(lmaxm,V,Dist,Ylmc,nlmc)
          call Ylmcr(lmaxm,nlmc,nlmr,Ylmc,Ylmr)
          Ylm_b(1:nlmr,j) = Ylmr(1:nlmr)

        end do ! end of loop over n_angle
      end do ! end of loop over n_rad

      call Cal_COOP_state(Comp_to_real_out,Full_potential,Harm_cubic,iapra,iaprb,icheck,index_r_a,index_r_b,lmaxa,lmaxb,lmm, &
                          nlmagm,nlmr,n_atom_0,n_atom_ind,n_pt_integr,n_rad_m,nspino,nspinp,Spinorbite,State, &
                          Tau_coop_r,um,Vol,Ylm_a,Ylm_b,Ylm_comp_out)  

      Coverlap(:,:,:,iab,mpirank) = State(:,:,:)
  
      deallocate( index_r_a, index_r_b, Ylm_a, Ylm_b, Ylmc, Ylmr )
          
    end do ! end of loop over atoms ib

    deallocate( Dlmm_ia )

  end do ! end of loop over atoms ia 
  
  deallocate( Radius_urm, um )

  return
 110 format(/' Tau_coop, iab =',i3,', i_trans =',i2)
 120 format(/' Tau_coop, iab =',i3,', isa, isb =',2i2,', i_trans =',i2)
 130 format(/' Tau_coop rotated, iab =',i3,', i_trans =',i2)
 140 format(/' Tau_coop rotated, iab =',i3,', isa, isb =',2i2,', i_trans =',i2)

end

!***********************************************************************

! Spin rotation matrix
! R = cos(Gamma/2).I + i sin(Gamma/2) Sigma.Axe
! I: Identity matrix
! Sigma: Pauli matrices.

subroutine Cal_rot_spin(Gamma_rot,icheck,Mat_rot_spin,Axe)

  use declarations
  implicit none

  integer:: icheck, isp
  complex(kind=db), dimension(2,2):: Mat_rot_spin

  real(kind=db):: Gamma_rot, sin_g, cos_g
  real(kind=db), dimension(3):: Axe

  cos_g = cos( Gamma_rot / 2 )
  sin_g = sin( Gamma_rot / 2 )
  
  Mat_rot_spin(1,1) = cos_g + img * sin_g * Axe(3)
  Mat_rot_spin(1,2) = sin_g * ( img * Axe(1) + Axe(2) )   
  Mat_rot_spin(2,1) = sin_g * ( img * Axe(1) - Axe(2) )   
  Mat_rot_spin(2,2) = cos_g - img * sin_g * Axe(3)

  if( icheck > 1 ) then
    write(3,110) Gamma_rot*180/pi, Axe(:)
    Write(3,'(/A)') ' Spin rotation matrix ='  
    do isp = 1,2 
      write(3,120) Mat_rot_spin(isp,:)
    end do 
  endif

 110 format(/' Spin projection angle =',f9.3, / &
             ' Rotation axe =',3f9.5)
 120 format(3x,2(1x,2f9.5))  
  return
end

!***********************************************************************

subroutine Cal_COOP_state(Comp_to_real,Full_potential,Harm_cubic,iapra,iaprb,icheck,index_r_a,index_r_b,lmaxa,lmaxb,lmm, &
                          nlmagm,nlmr,n_atom_0,n_atom_ind,n_pt_integr,n_rad_m,nspino,nspinp,Spinorbite,State, &
                          Taull,um,Vol,Ylm_a,Ylm_b,Ylm_comp)

  use declarations
  implicit none

  integer:: i, icheck, iapra, iaprb, ira, irb, is, is1, is2, iso1, iso2, isp, &
    L1, L2, lm01, lm02, lm1, lm2, lmax, lmaxa, lmaxb, lmm, lmv1, lmv2, Lp1, Lp2, m, m1, m2, mp1, mp2, mv1, mv2, &
    n_atom_0, n_atom_ind, n_pt_integr, n_rad_m, nlma, nlmagm, nlmb, nlmr, nspino, nspinp

  integer, dimension(n_pt_integr):: index_r_a, index_r_b
  
  complex(kind=db):: dCoop, rof_sd, Tau_e, Tau_s, YaYb, Yc, Ycompa, Ycompb
  complex(kind=db), dimension(16,16,nspinp):: Coop
  complex(kind=db), dimension(n_rad_m,lmm,nspinp,nspino,n_atom_0:n_atom_ind):: um
  complex(kind=db), dimension(n_rad_m,lmm,nspinp,nspino,2):: ump
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,nspino):: Taull
  complex(kind=db), dimension(:,:), allocatable:: Mat, Trans, Trans_a, Trans_b, Trans_i
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau
  
  logical:: Comp_to_real, Ylm_comp, Full_potential, Harm_cubic, New_title, Radial_comp, Spinorbite

  real(kind=db):: Vol
  real(kind=db), dimension(16,16,nspinp):: State, State_i
  real(kind=db), dimension(nlmr,n_pt_integr):: Ylm_a, Ylm_b

  nlma = ( lmaxa + 1 )**2
  nlmb = ( lmaxb + 1 )**2

  ump(:,:,:,:,1) = um(:,:,:,:,iapra)
  ump(:,:,:,:,2) = um(:,:,:,:,iaprb)

  Radial_comp = abs( aimag( ump(1,1,1,1,1) ) ) > eps10 .or. abs( aimag( ump(1,1,1,1,2) ) ) > eps10
    
  Coop(:,:,:) = 0._db
  State(:,:,:) = 0._db
  State_i(:,:,:) = 0._db
  
  if( icheck > 1 ) write(3,110) iapra, iaprb

  if( Comp_to_real ) then
! Matrix real to complex  
    lmax = max( lmaxa, lmaxb )
  
    allocate( Trans(-lmax:lmax,-lmax:lmax) )
    allocate( Trans_i(-lmax:lmax,-lmax:lmax) )
    Trans(:,:) = (0._db, 0._db)
    Trans_i(:,:) = (0._db, 0._db)

    do m = -lmax,lmax
      if( m > 0 ) then
        is = (-1)**m
        Trans(m,m) = cmplx( is * sqrt_1o2, 0._db, db)
        Trans(m,-m) = cmplx( 0._db, is * sqrt_1o2, db)
      elseif( m == 0 ) then
        Trans(0,0) = (1._db, 0._db)
      else
        Trans(m,-m) = cmplx(sqrt_1o2, 0._db, db)
        Trans(m,m) = cmplx(0._db, -sqrt_1o2, db)
      endif
    end do
     
    Trans_i = Conjg( Transpose( Trans ) )

  endif
  
  do L1 = 0,lmaxa
    lm01 = L1**2 + L1 + 1

    if( Comp_to_real ) then
      allocate( Trans_a(-L1:L1,-L1:L1) )
      Trans_a(-L1:L1,-L1:L1) = Trans_i(-L1:L1,-L1:L1)
    endif

    do L2 = 0,lmaxb
      lm02 = L2**2 + L2 + 1

      if( Comp_to_real ) then
        allocate( Trans_b(-L2:L2,-L2:L2) )
        Trans_b(-L2:L2,-L2:L2) = Trans(-L2:L2,-L2:L2)
      endif

      allocate( Mat(-L1:L1,-L2:L2) )          

      do isp = 1,nspinp
        do iso1 = 1,nspino
          if( Spinorbite ) then
            is1 = iso1
          else
            is1 = isp
          endif
          do iso2 = 1,nspino
            if( Spinorbite ) then
              is2 = iso2
            else
              is2 = isp
            endif

            do i = 1,n_pt_integr
                                      
              ira = index_r_a(i)
              irb = index_r_b(i)

              New_title = .true.
              Mat(:,:) = (0._db, 0._db )
              
              do m1 = -L1,L1

                do Lp1 = 0,lmaxa     ! The Lp1, mp1, Lp2, mp2 loops are not active yet
                  if( .not. Full_potential .and. Lp1 /= L1 ) cycle
                  if( Lp1 /= L1 ) cycle  ! Full_potential not acive
                  do mp1 = -Lp1, Lp1
                    if( mp1 /= m1 ) cycle
            
                    if( Spinorbite ) then
                      mv1 = m1 + iso1 - isp
                      if( mv1 < -L1 .or. mv1 > L1 ) cycle
                    else
                      mv1 = m1
                    endif
                
                    lmv1 = lm01 + mv1
                    lm1 = lm01 + m1
                
                    do m2 = -L2,L2
                
                      do Lp2 = 0,lmaxa
                        if( .not. Full_potential .and. Lp2 /= L2 ) cycle
                        if( Lp2 /= L2 ) cycle  ! Full_potential not acive
                        do mp2 = -Lp2, Lp2
                          if( mp2 /= m2 ) cycle
                
                          if( Spinorbite ) then
                            mv2 = m2 + iso2 - isp
                            if( mv2 < -L2 .or. mv2 > L2 ) cycle
                          else
                            mv2 = m2
                          endif
                      
                          lmv2 = lm02 + mv2
                          lm2 = lm02 + m2
                      
                          rof_sd = ump(ira,lm1,isp,iso1,1) * ump(irb,lm2,isp,iso2,2)
                      
                          Tau_e = Taull(lmv1,is1,lmv2,is2,1)
                          if( Spinorbite ) then
                            Tau_s = Taull(lmv2,is2,lmv1,is1,2)
                          elseif( Ylm_comp ) then
                            Tau_s = (-1)**(mv1-mv2) * Taull(lmv1-2*mv1,is1,lmv2-2*mv2,is2,1)
                          else
                            Tau_s = Tau_e
                          endif
                      
                          Mat(m1,m2) = Mat(m1,m2) + 0.5_db * ( rof_sd * Tau_e - conjg( rof_sd * Tau_s ) ) 
 
                        end do
                      end do
               
                    end do
                  end do

                end do
              end do

              if( Comp_to_real ) then
              
                Mat = matmul( Trans_a, matmul( Mat, Trans_b ) )

! in order it works, I do not understand why
                do m1 = -L1,L1
                  do m2 = -L2,L2
                    if( ( m1 < 0 .and. m2 >= 0 ) .or. ( m2 < 0 .and. m1 >= 0 ) ) Mat(m1,m2) = - Mat(m1,m2)   
                  end do
                end do
                
              endif

              do m1 = -L1,L1
                if( Spinorbite ) then
                  mv1 = m1 + iso1 - isp
                  if( mv1 < -L1 .or. mv1 > L1 ) cycle
                else
                  mv1 = m1
                endif
                lmv1 = lm01 + mv1
                lm1 = lm01 + m1

                do m2 = -L2,L2
                  if( Spinorbite ) then
                    mv2 = m2 + iso2 - isp
                    if( mv2 < -L2 .or. mv2 > L2 ) cycle
                  else
                    mv2 = m2
                  endif
                  lmv2 = lm02 + mv2
                  lm2 = lm02 + m2

                  if( Ylm_comp .and. m1 /= 0 .and. .not. Comp_to_real ) then
                    Ycompa = Yc(m1,Ylm_a(lm1,i),Ylm_a(lm1-2*m1,i))
                  else
                    Ycompa = cmplx( Ylm_a(lm1,i), 0._db, db )
                  endif
                  if( Ylm_comp .and. m2 /= 0 .and. .not. Comp_to_real ) then
                    Ycompb = Yc(m2,Ylm_b(lm2,i),Ylm_b(lm2-2*m2,i))
                  else
                    Ycompb = cmplx( Ylm_b(lm2,i), 0._db, db )
                  endif
                  YaYb = conjg( Ycompa ) * Ycompb

! In order it works, I do not understand why:                 
                  if( Ylm_comp .and. .not. Comp_to_real ) YaYb = img**(m1-m2) * YaYb

                  dCoop = YaYb * Mat(m1,m2) * Vol

                  Coop(lm1,lm2,isp) = Coop(lm1,lm2,isp) + dCoop

                  if( icheck > 1 .and. abs( dCoop ) > 1.e-15_db ) then
                    if( New_title ) then
                      write(3,120)
                      New_title = .false.
                    endif

                    write(3,130) lm1, lm2, isp, L1, m1, mv1, iso1, L2, m2, mv2, iso2, i, &
                        img*Coop(lm1,lm2,isp), img*dCoop, Ycompa, Ycompb, YaYb, Mat(m1,m2) 
                  endif

                end do  ! end of loop over m2
              end do ! end of loop over m1
          
            end do ! end of loop over points i

          end do  ! iso2
        end do  ! iso1
      
      end do  ! end of loop over isp
      
      deallocate( Mat )
      if( Comp_to_real ) deallocate( Trans_b )
      
    end do   ! end of loop over L2
    
    if( Comp_to_real ) deallocate( Trans_a )

  end do  ! end of loop over L1

  if( Comp_to_real ) deallocate( Trans, Trans_i )
    
  if( Harm_cubic ) then
    allocate( Tau(nlma,nspinp,nlmb,nspinp) )
    do isp = 1,nspinp
      do lm2 = 1,nlmb
        Tau(1:nlma,isp,lm2,isp) = Coop(1:nlma,lm2,isp) 
      end do
    end do
    call Trans_Tau_ZtoK(lmaxa,lmaxb,nspinp,Spinorbite,Tau)
    do isp = 1,nspinp
      do lm2 = 1,nlmb
        Coop(1:nlma,lm2,isp) = Tau(1:nlma,isp,lm2,isp) 
      end do
    end do
    deallocate( Tau )
  endif

  if( icheck > 1 ) then
    write(3,'(/A)') ' Integral for the COOP:' 

    do isp = 1,nspinp
      if( nspinp == 2 ) write(3,200) isp
      write(3,215) (( L2, m2, m2 = -L2,L2), L2 = 0,lmaxb)
      lm1 = 0
      do L1 = 0,lmaxa
        do m1 = -L1,L1
          lm1 = lm1 + 1
          write(3,230) L1, m1, ( img * Coop(lm1,lm2,isp), lm2 = 1,(1+lmaxb)**2 )
        end do
      end do
    end do
  endif

  State(:,:,:) = - aimag( Coop(:,:,:) )
  State_i(:,:,:) = real( Coop(:,:,:), db )
  
  return
  110 format(/' ia =',i3,', ib =',i3,/)
  120 format(/' lm1 lm2 isp  L1  m1 mv1 is1  L2  m2 mv2 is2   i',5x,'Coop_r',7x,'Coop_i',7x,'dCoop_r',6x,'dCoop_i',15x,'Ylm_a', &
              22x,'Ylm_b',24x,'Yab',24x,'Mat')
  130 format(12i4,6(1x,2e13.5))
  200 format(/' isp =',i2)
  215 format(/' L1 m1  ',100(9x,2i3,12x))
  220 format(2i3,1p,28e13.5)
  230 format(2i3,1p,28(1x,2e13.5))

end

!***********************************************************************

! Calculation of the density of state for covelap, old version. Not used.

subroutine Cal_COOP_state_old(Comp_to_real,Full_potential,Harm_cubic,iapra,iaprb,icheck,index_r_a,index_r_b,lmaxa,lmaxb,lmm, &
                          nlmagm,nlmr,n_atom_0,n_atom_ind,n_pt_integr,n_rad_m,nspino,nspinp,Spinorbite,State, &
                          Taull,um,Vol,Ylm_a,Ylm_b,Ylm_comp)

  use declarations
  implicit none

  integer:: i, i_k1, i_k2, icheck, iapra, iaprb, ira, irb, is1, is2, isg1, isg2, iso1, iso2, isp, isp1, isp2, iss1, iss2, &
    L1, L2, lm01, lm02, &
    lm1, lm2, lma1, lma2, lmaxa, lmaxb, lmm, Lp1, Lp2, m_k1, m_k2, m_r1, m_r2, m1, m2, mp1, mp2, mv1, mv2, &
    n_atom_0, n_atom_ind, n_k1, n_k2, n_pt_integr, n_rad_m, n1, n2, nlmagm, nlmr, n_rk1, n_rk2, &
    nspino, nspinp

  integer, dimension(n_pt_integr):: index_r_a, index_r_b
  
  complex(kind=db):: c_harm, c_harm1, c_harm2, rof_sd, Sta, Tau_e, Tau_s, YaYb, Yc, Ycompa, Ycompb
  complex(kind=db), dimension(n_rad_m,lmm,nspinp,nspino,n_atom_0:n_atom_ind):: um
  complex(kind=db), dimension(n_rad_m,lmm,nspinp,nspino,2):: ump
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,nspino):: Taull
  
  logical:: Comp_to_real, Ylm_comp, Full_potential, Harm_cubic, New_title, Radial_comp, Spinorbite

  real(kind=db):: c_cubic1, c_cubic2, Sta_i, Sta_r, Vol
  real(kind=db), dimension(16,nspinp,16,nspinp):: State, State_i
  real(kind=db), dimension(nlmr,n_pt_integr):: Ylm_a, Ylm_b

  ump(:,:,:,:,1) = um(:,:,:,:,iapra)
  ump(:,:,:,:,2) = um(:,:,:,:,iaprb)

  Radial_comp = abs( aimag( ump(1,1,1,1,1) ) ) > eps10 .or. abs( aimag( ump(1,1,1,1,2) ) ) > eps10
    
  State(:,:,:,:) = 0._db
  State_i(:,:,:,:) = 0._db
  
  if( icheck > 1 ) write(3,110) iapra, iaprb

  n_k1 = 0
  lm1 = 0
  do L1 = 0,lmaxa
    lm01 = L1**2 + L1 + 1
    do m_k1 = -L1,L1
      n_k1 = n_k1 + 1
      lm1 = lm1 + 1

      do isp1 = 1,nspinp

        n_k2 = 0
        lm2 = 0
        do L2 = 0,lmaxb
          lm02 = L2**2 + L2 + 1
          do m_k2 = -L2,L2
            n_k2 = n_k2 + 1
            lm2 = lm2 + 1
 
            do isp2 = 1,nspinp
 
              New_title = .true.
                    
              if( .not. Harm_cubic .or. L1 /= 3 .or. ( m_k1 == -3 .or. m_k1 == 0 .or. m_k1 == 3 ) ) then
                n_rk1 = 1
              else
                n_rk1 = 2
              endif
          
              do i_k1 = 1,n_rk1
          
                call Trans_TtoK(c_cubic1,Harm_cubic,i_k1,L1,m_k1,m_r1)
          
                do isg1 = -1,1,2
          
                  if( ( .not. Comp_to_real .or. m_r1 == 0 ) .and. isg1 == -1 ) cycle
          
                  call Trans_RtoC(c_harm1,isg1,m1,m_r1,Comp_to_real)
          
                  do Lp1 = 0,lmaxa
                    if( .not. Full_potential .and. Lp1 /= L1 ) cycle
          
                    do mp1 = -Lp1, Lp1
                      if( mp1 /= m1 ) cycle
          
                      do iso1 = 1,nspino
          
                        if( Spinorbite ) then
                          is1 = iso1
                        else
                          is1 = isp1
                        endif
                        iss1 = isp1
          
                        if( Spinorbite ) then
                          mv1 = m1 + iso1 - isp1
                          if( mv1 > L1 .or. mv1 < - L1 ) cycle
                        else
                          mv1 = m1
                        endif                  
          
                        lma1 = lm01 + mv1
                        n1 = lm01 + m1
          
                        if( .not. Harm_cubic .or. L2 /= 3 .or. ( m_k2 == -3 .or. m_k2 == 0 .or. m_k2 == 3 ) ) then
                          n_rk2 = 1
                        else
                          n_rk2 = 2
                        endif
          
                        do i_k2 = 1,n_rk2
          
                          call Trans_TtoK(c_cubic2,Harm_cubic,i_k2,L2,m_k2,m_r2)
          
                          do isg2 = -1,1,2
          
                            if( ( .not. Comp_to_real .or. m_r2 == 0 ) .and. isg2 == -1 ) cycle
                     
                            call Trans_CtoR(c_harm2,isg2,m2,m_r2,Comp_to_real)

                            c_harm = conjg( c_harm1 ) * c_cubic1 * c_harm2 * c_cubic2
         
                            do Lp2 = 0,lmaxb
                              if( Lp2 /= L2 ) cycle
          
                              do mp2 = -Lp2,Lp2
                                if( mp2 /= m2 ) cycle
          
                                do iso2 = 1,nspino
                                  if( isp1 /= isp2 ) cycle ! in a scalar product there is no spin-crossing (but there is solution crossing)
                                  if( .not. Spinorbite .and. isp1 /= isp2 ) cycle

                                  if( Spinorbite ) then
                                    is2 = iso2
                                  else
                                    is2 = isp2
                                  endif
                                  iss2 = isp2

                                  if( Spinorbite ) then
                                    mv2 = m2 + iso2 - isp2  ! for tau
                                    if( mv2 > L2 .or. mv2 < - L2 ) cycle
                                  else
                                    mv2 = m2
                                  endif                  

                                  lma2 = lm02 + mv2
                                  n2 = lm02 + m2
    
                                  do i = 1, n_pt_integr
                                                                
                                    ira = index_r_a(i)
                                    irb = index_r_b(i)
                                      
                                    rof_sd = ump(ira,n1,isp1,iso1,1) * ump(irb,n2,isp2,iso2,2) * Vol
                                    
                                    if( Ylm_comp .and. m1 /= 0 ) then
                                      Ycompa = Yc(m1,Ylm_a(n1,i),Ylm_a(n1-2*m1,i))
                                    else
                                      Ycompa = cmplx( Ylm_a(lm01+m_r1,i), 0._db, db )
                                    endif
                                    if( Ylm_comp .and. m2 /= 0 ) then
                                      Ycompb = Yc(m2,Ylm_b(n2,i),Ylm_b(n2-2*m2,i))
                                    else
                                      Ycompb = cmplx( Ylm_b(lm02+m_r2,i), 0._db, db )
                                    endif
                                    YaYb = conjg( Ycompa ) * Ycompb
                                
                                    if( abs( YaYb ) < eps10 ) cycle

                                    Tau_e = Taull(lma1,is1,lma2,is2,1)
                                    if( Spinorbite ) then
                                      Tau_s = Taull(lma2,is2,lma1,is1,2)
                                    elseif( Ylm_comp ) then
                                      Tau_s = (-1)**(mv1-mv2) * Taull(lma1-2*mv1,is1,lma2-2*mv2,is2,1)
                                    else
                                      Tau_s = Tau_e
                                    endif

                                    Sta = c_harm * YaYb * 0.5_db * ( rof_sd * Tau_e - conjg( rof_sd * Tau_s ) )

! The following is in order it works, I do not understand why
!                                    if( Ylm_comp ) Sta = img**(mv1-mv2) * Sta
                                    if( Ylm_comp ) Sta = img**(m1-m2) * Sta

                                    Sta_i = - aimag( Sta )
                                    Sta_r = real( Sta, db )  
                                
                                    State(lm1,iss1,lm2,iss2) = State(lm1,iss1,lm2,iss2) + Sta_i
                                    State_i(lm1,iss1,lm2,iss2) = State_i(lm1,iss1,lm2,iss2) + Sta_r
                                                                    
                                    if( icheck > 1 .and. abs( Sta_i ) > 1.e-15_db ) then
                                      if( New_title ) then
                                        if( Radial_comp ) then
                                          write(3,120)
                                        else
                                          write(3,130)
                                        endif
                                        New_title = .false.
                                      endif

                                      if( Radial_comp ) then
                                        write(3,140) lm1, iss1, lm2, iss2, lma1, is1, lma2, is2, L1, m_k1, m_r1, m1, isp1, &
                                          L2, m_k2, m_r2, m2, isp2, iso1, iso2, i, State(lm1,iss1,lm2,iss2), Sta_i, Sta_r, &
                                          c_harm, Ycompa, Ycompb, & 
                                          ump(ira,n1,isp1,iso1,1), ump(irb,n2,isp2,iso2,2), &
                                          YaYb * c_harm * rof_sd, Tau_e, conjg( Tau_s )
                                      elseif( Ylm_comp ) then
                                        write(3,150) lm1, iss1, lm2, iss2, lma1, is1, lma2, is2, L1, m_k1, m_r1, m1, isp1, &
                                          L2, m_k2, m_r2, m2, isp2, iso1, iso2, i, State(lm1,iss1,lm2,iss2), Sta_i, Sta_r, &
                                          c_harm, Ycompa, Ycompb, & 
                                          real(ump(ira,n1,isp1,iso1,1),db), real(ump(irb,n2,isp2,iso2,2),db), &
                                          YaYb * c_harm * rof_sd, Tau_e, conjg( Tau_s ), Tau_e - conjg( Tau_s )
                                      else
                                        write(3,150) lm1, iss1, lm2, iss2, lma1, is1, lma2, is2, L1, m_k1, m_r1, m1, isp1, &
                                          L2, m_k2, m_r2, m2, isp2, iso1, iso2, i, State(lm1,iss1,lm2,iss2), Sta_i, Sta_r, &
                                          c_harm, Ycompa, Ycompb, & 
                                          real(ump(ira,n1,isp1,iso1,1),db), real(ump(irb,n2,isp2,iso2,2),db), &
                                          YaYb * c_harm * rof_sd, Tau_e, conjg( Tau_s )
                                      endif
                                    endif

                                  end do ! loop i (integral)
        
                                end do  ! boucle iso2
                              end do  ! boucle mp2
                            end do  ! boucle Lp2
                          end do  ! boucle isg2
                        end do  ! boucle ik_2

                      end do  ! boucle iso1
                    end do  ! boucle mp1
                  end do  ! boucle Lp1
                end do  ! boucle isg1
              end do  ! boucle i_k1

            end do  ! boucle isp2
          end do  ! boucle m_k2
        end do  ! boucle L2

      end do  ! boucle isp1
    end do  ! boucle m_k1
  end do  ! boucle L1

  if( icheck > 1 ) then
    write(3,'(/A)') ' Tau_coop:'

    if( .not. Full_potential ) then
      write(3,160) ((( L2, m2, isp, m2 = -L2,L2), L2 = 0,lmaxb), isp = 1,nspinp )
    else
      write(3,170)
    endif
    do isp1 = 1,nspinp
      lm1 = 0
      do L1 = 0,lmaxa
        do m1 = -L1,L1
          lm1 = lm1 + 1
          if( .not. Full_potential ) then
            write(3,180) L1, m1, isp1, ((Taull(lm1,isp1,lm2,isp2,1), lm2 = 1,(1+lmaxb)**2), isp2 = 1,nspinp )
          else
            write(3,180) L1, m1, isp1, (( Taull(lm1,isp1,lm2,isp2,1), lm2 = 1,(1+lmaxb)**2), isp2 = 1,nspinp )
          endif
        end do
      end do
    end do

    write(3,'(/A)') ' Integral for the COOP:' 

    write(3,190) ((( L2, m2, isp, m2 = -L2,L2), L2 = 0,lmaxb), isp = 1,nspinp )
    do isp1 = 1,nspinp
      lm1 = 0
      do L1 = 0,lmaxa
        do m1 = -L1,L1
          lm1 = lm1 + 1
          if( .not. Full_potential ) then
            write(3,200) L1, m1, isp1, (( State(lm1,isp1,lm2,isp2), lm2 = 1,(1+lmaxb)**2 ), isp2 = 1,nspinp )
          else
            write(3,180) L1, m1, isp1, (( State(lm1,isp1,lm2,isp2), State_i(lm1,isp1,lm2,isp2), lm2 = 1,(1+lmaxb)**2), &
                                                                                                isp2 = 1,nspinp )
          endif
        end do
      end do
    end do
  endif
  
  return
      
  110 format(/' ia =',i3,', ib =',i3,/)
  120 format(' lm1 iss1 lm2 iss2 lma1 is1 lma2 is2 L1 mk1 mr1 m1 isp1 L2 mk2 mr2 m2 isp2 iso1 iso2 i',8x, &
             'Coop',8x,'dcoop',9x,'Sta_i',15x,'c_harm',23x,'Ylma',23x,'Ylmb',24x,'ura',23x,'urb',11x, &
             'c_harm.Ylma^*.Ylmb.ura.urb.dV',10x,'Taull',17x,'conj(Taull_tr)',18x,'dif')
  130 format(' lm1 iss1 lm2 iss2 lma1 is1 lma2 is2 L1 mk1 mr1 m1 isp1 L2 mk2 mr2 m2 isp2 iso1 iso2 i',8x, &
             'Coop',8x,'dcoop',8x,'dcoop_i',14x,'c_harm',23x,'Ylma',23x,'Ylmb',16x,'ura',12x,'urb',4x, &
             'c_harm.Ylma^*.Ylmb.ura.urb.dV',10x,'Taull',17x,'conj(Taull_tr)',18x,'dif')
  140 format(i3,i5,i4,i5,i5,i4,i5,i4,i4,1x,2(2i3,3i4,1x),2i4,1x,1p,3(1x,e13.5),8(1x,2e13.5))
  150 format(i3,i5,i4,i5,i5,i4,i5,i4,i4,1x,2(2i3,3i4,1x),2i4,1x,1p,3(1x,e13.5),3(1x,2e13.5),2(1x,e13.5),4(1x,2e13.5))
  160 format('  L  m isp ',100(8x,3i3,10x))
  170 format('  L  m isp ',6x,'Taull_r',6x,'Taull_i')
  180 format(3i3,1p,28(1x,2e13.5))
  190 format(' L1 m1 isp1',100(3i3,4x))
  200 format(3i3,1p,28e13.5)
end

!***********************************************************************

! Writing of the COOP

subroutine Write_coop(Coverlap,Density_comp,Energ,Harm_cubic,iaprotoi,ie,ie_computer,index_coop,itypepr,mpinodee, &
                      n_atom_proto,nab_coop,natome,nomfich_s,nspin,nspinp,ntype,numat)

  use declarations
  implicit none
 
  integer:: ia, iab, ib, ie, ie_computer, index, isp, it, L, la, lam, lb, lbm, lma, lmax_coop, lmb, M, ma, mb, mpinodee, &
            n_atom_proto, nab_coop, natome, nspin, nspinp, ntype

  integer, dimension(natome):: iaprotoi
  integer, dimension(0:ntype):: numat
  integer, dimension(0:n_atom_proto):: itypepr
  integer, dimension(2,nab_coop):: index_coop

  character(len=9):: mot9
  character(len=15):: mot15
  character(len=Length_name):: nomfich_coop, nomfich_s
  character(len=2), dimension(2):: Sp
  character(len=1), dimension(0:4):: Orb_L
  character(len=9), dimension(16):: Orb_Lm, Orb_Lm_D, Orb_Lm_K, Orb_Lm_Z
  character(len=15), dimension((256+16)*nspinp):: Col_name
  
  logical:: Density_comp, Harm_cubic, Spin_out

  real(kind=db):: COverlap_T, Energ
  real(kind=db), dimension(16,16,nspinp,nab_coop,0:mpinodee-1):: Coverlap
  real(kind=db), dimension(:,:,:), allocatable:: COverlap_l

  data Orb_L/'s','p','d','f','g'/
  data Orb_Lm_K/'s','px','py','pz','dx2-y2','dz2','dyz','dxz','dxy','fxyz','fx3','fy3','fz3','fx(y2-z2)','fy(z2-x2)','fz(x2-y2)'/
  data Orb_Lm_Z/'s','py','pz','px','dxy','dyz','dz2','dxz','dx2-y2','(3,-3)','fxyz','(3,-1)','fz3','(3,1)','fz(x2-y2)','(3,3)'/
  data Orb_Lm_D/'(0,0)','(1,-1)','(1,0)','(1,1)','(2,-2)','(2,-1)','(2,0)','(2,1)','(2,2)','(3,-3)','(3,-2)','(3,-1)','(3,0)', &
                '(3,1)','(3,2)','(3,3)'/
  data Sp/'up','dn'/

  Spin_out = nspin == 2 .or. ( nspinp == 2 .and. Density_comp )
  
  if( Density_comp) then
    Orb_Lm(:) = Orb_Lm_D
  elseif( Harm_cubic ) then
    Orb_Lm(:) = Orb_Lm_K
  else
    Orb_Lm(:) = Orb_Lm_Z
  endif
                        
 ! Writing
  do iab = 1,nab_coop

    ia = index_coop(1,iab)
    ib = index_coop(2,iab)

    nomfich_coop = nomfich_s
    L = len_trim(nomfich_coop)
    nomfich_coop(L+1:L+5) = '_coop'
    L = len_trim(nomfich_coop)
    nomfich_coop(L+1:L+1) = '_'
    call ad_number(ia,nomfich_coop,Length_name)
    L = len_trim(nomfich_coop)
    nomfich_coop(L+1:L+1) = '_'
    call ad_number(ib,nomfich_coop,Length_name)
    L = len_trim(nomfich_coop)
    nomfich_coop(L+1:L+4) = '.txt'
 
    it = itypepr( iaprotoi(ia) ) 
    lam = lmax_coop( numat( it ) )
    it = itypepr( iaprotoi(ib) ) 
    lbm = lmax_coop( numat( it ) )

    allocate( COverlap_l(0:lam,0:lbm,nspinp) )
    COverlap_l(:,:,:) = 0._db
    COverlap_T = 0._db
     
    do isp = 1,nspinp
      lma = 0
      do la = 0,lam
        do ma = -la,la
          lma = lma + 1
          lmb = 0
          do lb = 0,lbm
            do mb = -lb,lb
              lmb = lmb + 1
              COverlap_l(la,lb,isp) = COverlap_l(la,lb,isp) + COverlap(lma,lmb,isp,iab,ie_computer) 
              COverlap_T = COverlap_T + COverlap(lma,lmb,isp,iab,ie_computer)
            end do
          end do
        end do
      end do
    end do
              
    if( ie == 1 ) then
      open(4, file = nomfich_coop )

      mot15 = 'Tot_COOP_'
      call ad_number(ia,mot15,15)
      L = len_trim(mot15) + 1
      if( L < 15 ) mot15(L:L) = '_' 
      call ad_number(ib,mot15,15)
      call Center_word(mot15,15)
      index = 1
      Col_name(index) = mot15 
           
      if( Spin_out ) then
        do isp = 1,nspinp
          do la = 0,lam
            do lb = 0,lbm
              mot15 = Orb_L(la) // ':' // Orb_L(lb) // '_' // Sp(isp)
              call Center_word(mot15,15)
              index = index + 1
              Col_name(index) = mot15 
            end do
          end do
        end do
      else
        do la = 0,lam
          do lb = 0,lbm
            mot15 = Orb_L(la) // ':' // Orb_L(lb)
            call Center_word(mot15,15)
            index = index + 1
            Col_name(index) = mot15 
          end do
        end do
      endif

      if( Spin_out ) then
        do isp = 1,nspinp
          lma = 0        
          do la = 0,lam
            do ma = -la,la
              lma = lma + 1
              lmb = 0
              do lb = 0,lbm
                do mb = -lb,lb
                  lmb = lmb + 1
                  mot15 = Orb_Lm(lma)
                  L = len_trim(mot15) + 1
                  mot15(L:L) = ':'
                  mot9 = Orb_Lm(lmb)
                  M = min( 14 - L, len_trim(mot9) ) 
                  mot15(L+1:L+M) = mot9(1:M)
                  L = L + M + 1 
                  if( L < 13 ) then
                    mot15(L:L) = '_'
                    L = L + 1
                  endif
                  if( L < 14 ) mot15(L:L+1) = Sp(isp)
                  call Center_word(mot15,15)
                  index = index + 1
                  Col_name(index) = mot15 
                end do
              end do
            end do
          end do
        end do                    
      else
        lma = 0
        do la = 0,lam
          do ma = -la,la
            lma = lma + 1
            lmb = 0
            do lb = 0,lbm
              do mb = -lb,lb
                lmb = lmb + 1
                mot15 = Orb_Lm(lma)
                L = len_trim(mot15) + 1
                mot15(L:L) = ':'
                mot9 = Orb_Lm(lmb)
                M = min( 14 - L, len_trim(mot9) ) 
                mot15(L+1:L+M) = mot9(1:M) 
                call Center_word(mot15,15)
                index = index + 1
                Col_name(index) = mot15 
              end do
            end do
          end do
        end do
      endif
          
      write(4,110) Col_name(1:index)
    else
      open(4, file = nomfich_coop, position='append')
    endif
    if( Spin_out ) then
      write(4,130) Energ*rydb, COverlap_T, & 
                     ((( COverlap_l(la,lb,isp), lb = 0,lbm), la = 0,lam), isp = 1,nspinp), & 
                     ((( COverlap(lma,lmb,isp,iab,ie_computer), &
                     lmb = 1,( 1 + lbm )**2), lma = 1,( 1 + lam )**2), isp = 1,nspinp)
    else
      write(4,130) Energ*rydb, 2*COverlap_T, & 
                     (( 2*COverlap_l(la,lb,1), lb = 0,lbm), la = 0,lam), & 
                     (( 2*COverlap(lma,lmb,1,iab,ie_computer), lmb = 1,( 1 + lbm )**2), lma = 1,( 1 + lam )**2)
    endif
     
    Close(4)
    
    deallocate( COverlap_l )
  end do
  
  return
  110 format('   Energy ',1000a15)
  130 format(f10.5,1p,1000e15.7)
end

!***********************************************************************

! Calculation of the radial wave function on the points  between the atoms
! One neglects the non spherical potential

subroutine Cal_urm(Ecinetic,Eimag,Energ,Enervide, &
            Full_potential,Hubb_a,Hubb_d,iapr,icheck,lmax,lmax_pot,lmm,m_hubb,n_atom_0,n_atom_ind,n_rad_m,n_radius,nlm_pot, &
            nr,nspin,nspino,nspinp,r,Radius_urm,Relativiste,Renorm,Rmtg,Rmtsd,Spinorbite,um,V_hubb,V_intmax,V0bd,&
            Vrato,Ylm_comp,Z)

  use declarations
  implicit none

  integer:: i, iapr, icheck, ir, iso, isp, L, l_hubbard, lfin, lm, lm0, lmax, lmax_pot, lmm, m, m_hubb, mp, n, n_atom_0, &
    n_atom_ind, n_rad_m, nlm, nlm_pot, nlm1, nlm2, np, nr, nrmtg, nrmtsd, nspin, nspino, nspinp, Z

  integer, dimension(n_atom_0:n_atom_ind):: n_radius
  integer, dimension(n_rad_m):: index_r

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(n_rad_m,lmm,nspinp,nspino,n_atom_0:n_atom_ind):: um
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  logical:: Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m, Radial_comp, Relativiste, Renorm, &
    Spinorbite, Ylm_comp

  real(kind=db):: Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(n_rad_m):: p
  real(kind=db), dimension(n_rad_m,n_atom_0:n_atom_ind):: Radius_urm
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato

  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nr,nspin):: g0, gm, gp
  real(kind=db), dimension(nr,nspino):: gso

  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur

  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag, db) )

  if( icheck > 1 ) write(3,110) iapr, Z, lmax
  if( icheck > 1 ) write(3,120) Energ*rydb, Enervide*rydb
  if( icheck > 2 ) then
    write(3,130) Ecinetic(:)*rydb
    write(3,140) V0bd(:)*rydb
    write(3,150) konde(:)
  endif

  if( Full_potential ) then
    nlm = ( lmax + 1 )**2
  else
    nlm = 1
  endif
  allocate( V(nr,nlm,nlm,nspin) )

  call mod_V(icheck,lmax,lmax_pot,nlm,nlm_pot,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,Z,r,Relativiste,Spinorbite,V)

  gp(1:nr-1,:) = 1 / gp(1:nr-1,:)

  if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10 .or. Ecinetic(nspin) < eps10 ) then
    Ecomp = .true.
  else
    Ecomp = .false.
  endif
  Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

  if( Full_potential ) then
    lfin = 0
  else
    lfin = lmax
  endif

  do i = 1,n_radius(iapr)
    do ir = 2,nr - 1
      if( r(ir) > Radius_urm(i,iapr) ) exit
    end do
    index_r(i) = ir
    p(i) = ( Radius_urm(i,iapr) - r(ir-1) ) / ( r(ir) - r(ir-1) )
  end do
    
  do L = 0,lfin

    if( Hubb_a .and. L == l_hubbard( Z ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = nlm   ! = (lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*L + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*L + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,L,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,Z,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

    do i = 1,n_radius(iapr)

      ir = index_r(i)
      if( Radial_comp ) then
        if( Full_potential ) then
          do lm = 1,nlm1
            um(i,lm,:,:,iapr) = cmplx( p(i) * ur(ir,lm,lm,:,:) + (1 - p(i)) * ur(ir-1,lm,lm,:,:), &
                                       p(i) * ui(ir,lm,lm,:,:) + (1 - p(i)) * ui(ir-1,lm,lm,:,:), db )
          end do
        else
          lm0 = L**2 + L + 1
          do m = -L,L
            if( nlm1 == 1 ) then
              n = 1
            else
              n = L + 1 + m
            endif
            do mp = -L,L
!              if( nlm2 == 1 .and. mp /= m ) cycle
              if( mp /= m ) cycle
              np = min(L + m + 1, nlm2 )
              um(i,lm0+m,:,:,iapr) = cmplx( p(i) * ur(ir,n,np,:,:) + (1 - p(i)) * ur(ir-1,n,np,:,:), &
                                            p(i) * ui(ir,n,np,:,:) + (1 - p(i)) * ui(ir-1,n,np,:,:), db )
            end do
          end do
        endif
      else
        if( Full_potential ) then
          do lm = 1,nlm1
            um(i,lm,:,:,iapr) = cmplx( p(i) * ur(ir,lm,lm,:,:) + (1 - p(i)) * ur(ir-1,lm,lm,:,:), 0._db, db )
          end do
        else
          lm0 = L**2 + L + 1
          do m = -L,L
            if( nlm1 == 1 ) then
              n = 1
            else
              n = L + 1 + m
            endif
            do mp = -L,L
!              if( nlm2 == 1 .and. mp /= m ) cycle
              if( mp /= m ) cycle
              np = min(L + m + 1, nlm2 )
              um(i,lm0+m,:,:,iapr) = cmplx( p(i) * ur(ir,n,np,:,:) + (1 - p(i)) * ur(ir-1,n,np,:,:), 0._db, db )
            end do
          end do
        endif
      endif
    
    end do
    
    deallocate( Tau, ui, ur )

  end do   ! end of loop over L

  if( icheck > 1 ) then
    write(3,155) iapr
    write(3,160) ( Radius_urm(i,iapr)*bohr, i = 1,n_radius(iapr) )
    Lm = 0
    do L = 0,lmax
      do m = -L,L
        Lm = Lm + 1
        do isp = 1, nspinp
          do iso = 1,nspino 
            write(3,170) L, m, isp, iso, ( um(i,Lm,isp,iso,iapr), i = 1,n_radius(iapr) )
          end do
        end do
      end do
    end do 
  endif
  
  deallocate( V )

  return
  110 format(/' iapr =',i3,', Z =',i3,', lmax =',i2)
  120 format(' Energ =',f10.3,' eV,  Enervide =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format(' V0bd  =',2f10.3)
  150 format(' konde =',2f12.5)
  155 format(/' Radial wave function for the coop versus radius, iapr =',i3)
  160 format('  L  m is io',10(10x,f9.5,8x))
  170 format(4i3,1p,20(1x,2e13.5))
end

!***********************************************************************

! Calculation of Dlmm for the COOP
! The other routine doing that (in mat.f90) seems not working 
! 6th January 2019

subroutine Dlmm_COOP(Dlmm,icheck,lmax,Mat_rot,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, is, L, lmax, m, mp, s
  
  integer, dimension(0:2*lmax):: Fact

  logical:: Ylm_comp
  
  complex(kind=db), dimension(-lmax:lmax,-lmax:lmax):: Mat, Trans, Transi
  complex(kind=db), dimension(-lmax:lmax,-lmax:lmax,0:lmax):: Dlmm

  real(kind=db):: a, alfa, b, beta, c, cos_b, f, gamma, sin_b
  real(kind=db), dimension(-lmax:lmax,-lmax:lmax,0:lmax):: R
  real(kind=db), dimension(3,3):: Mat_rot

  Dlmm(:,:,:) = (0._db, 0._db)
  
  call Euler_mat(icheck,Mat_rot,alfa,beta,gamma)

  Fact(0) = 1
  do L = 1,2*lmax
    Fact(L) = L * Fact(L-1)
  end do

  if( icheck > 2 ) write(3,100) alfa*180/pi, beta*180/pi, gamma*180/pi

  R(:,:,:) = 0._db

  cos_b = cos(beta/2)
  sin_b = sin(beta/2)
  
  do L = 0,lmax
    do m = -L,L
      do mp = -L,L
        f = sqrt( real( Fact(L+m) * Fact(L-m) * Fact(L+mp) * Fact(L-mp), db ) )
        
        do s = 0,2*L
          if( (L+mp-s) < 0 .or. (m-mp+s) < 0 .or. (L-m-s) < 0 ) cycle
          
          a = real( (-1)**(m-mp+s), db) / ( Fact(L+mp-s) * Fact(s) * Fact(m-mp+s) * Fact(L-m-s) )
          b = cos_b**(2*L+mp-m-2*s)
          c = Sin_b**(m-mp+2*s)
  
          R(m,mp,L) = R(m,mp,L) + f * a * b * c   

          if( icheck > 3 ) write(3,'(a14,4i3,f10.5)') ' L,m,mp,s,R =', L, m, mp, s, R(m,mp,L)     
        end do

      end do
    end do
  end do
  
  do L = 0,lmax
    do m = -L,L
      do mp = -L,L
        Dlmm(m,mp,L) = exp( -img * ( m*alfa + mp*gamma) ) * R(m,mp,L)
      end do
    end do
  end do

 ! Case of real harmonics:
  if( .not. Ylm_comp ) then
 
    Trans(:,:) = (0._db, 0._db)
    Transi(:,:) = (0._db, 0._db)
    do m = -lmax,lmax
      if( m > 0 ) then
        is = (-1)**m
        Trans(m,m) = cmplx( is * sqrt_1o2, 0._db, db)
        Trans(m,-m) = cmplx( 0._db, is * sqrt_1o2, db)
      elseif( m == 0 ) then
        Trans(0,0) = (1._db, 0._db)
      else
        Trans(m,-m) = cmplx( sqrt_1o2, 0._db, db)
        Trans(m,m) = cmplx(0._db, -sqrt_1o2, db)
      endif
    end do
    
    Transi = Conjg( Transpose( Trans ) )
      
    do L = 0,lmax

! I do not really understand the conjg
      do m = -L,L
        do mp = -L,L
          Mat(m,mp) = sum( conjg( Dlmm(m,-L:L,L) ) * Trans(-L:L,mp) )
        end do
      end do

      do m = -L,L
        do mp = -L,L
          Dlmm(m,mp,L) = sum( Transi(m,-L:L) * Mat(-L:L,mp) )
        end do
      end do

    end do
  
  endif
  
  if( icheck > 2 ) then
  
    do L = 0,lmax
      write(3,110) '    R, L =',L
      do m = -L,L
        Write(3,120) m, ( R(m,mp,L), mp = -L,L )
      end do

      write(3,110) ' Dlmm, L =',L
      do m = -L,L
        Write(3,130) m, ( Dlmm(m,mp,L), mp = -L,L )
      end do
    end do

  endif
  
  return
 100 format(/' alfa, beta, gamma =',3f10.3)
 110 format(/a10,i2)
 120 format(i3,5f9.5)
 130 format(i3,5(1x,2f8.5))
end

!***********************************************************************

! Called by main_tddft

subroutine Radial_wave(Ecinetic,Eimag,Energ,Enervide,Full_potential,Hubb_a,Hubb_d,icheck,initl,lmax,lmax_pot, &
            m_hubb,n_comp,nbseuil,ninit1,ninitt,nlm_pot,nlmam,nlmam2,nr,nrm,nspin,nspino,nspinp,numat,psii,r,Radial_comp, &
            Relativiste,Renorm,Rmtg,Rmtsd,rof_ph,Spinorbite,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp,zet)

  use declarations
  implicit none

  integer:: i_comp, icheck, initl, iseuil, iso, isp, L, l_hubbard, lfin, lm, lm0, lmax, lmax_pot, lmp, m, m_hubb, mp, &
    n, n_comp, nbseuil, ninit1, ninitt, nlm, nlm_pot, nlm1, nlm2, nlmam, nlmam2, np, nr, nrm, nrmtg, nrmtsd, nspin, &
    nspino, nspinp, numat

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(n_comp*nlmam,nlmam2,nspinp,nspino,0:0,1):: rof
  complex(kind=db), dimension(n_comp*nlmam,nspinp,nspino,ninitt):: rof_ph
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau

  logical:: Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m, NRIXS, Radial_comp, Relativiste, Renorm, Spinorbite, Ylm_comp

  real(kind=db):: Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nbseuil):: Vecond
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato

  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nr,nspin):: g0 , gm, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr,0:0):: r_or_bess
  real(kind=db), dimension(nr,nlmam,nlmam2,nspinp,nspino,ninitt):: zet

  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur
  
  rof(:,:,:,:,:,:) = ( 0._db, 0._db )
  
  NRIXS = .false.
  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

  if( nbseuil == 2 ) then
    if( nbseuil /= ninitt ) then
      ninit1 = ( ninitt - 2 ) / 2
      if( initl <= ninit1 ) then
        iseuil = 1
      else
        iseuil = 2
      endif
    else 
      iseuil = initl
    endif
  else
    iseuil = 1
  endif
 ! Fake value not used because only monopole calculated
  Vecond(:) = 0._db 

  r_or_bess(:,:) = 1._db
  
  if( icheck > 1 ) then
    write(3,110)
    write(3,120) Energ*rydb, Enervide*rydb
    write(3,130) Ecinetic(:)*rydb
    write(3,140) V0bd(:)*rydb
    write(3,150) konde(:)
  endif

  if( Full_potential ) then
    nlm = ( lmax + 1 )**2
  else
    nlm = 1
  endif
  allocate( V(nr,nlm,nlm,nspin) )

  call mod_V(icheck,lmax,lmax_pot,nlm,nlm_pot,nr,nrmtg,nrmtsd,nspin,r,Rmtg,Rmtsd,V,V_intmax,V0bd,Vrato,Ylm_comp)

  call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)

  gp(:,:) = 1 / gp(:,:)

  if( abs(Eimag) > eps10 .or. Ecinetic(1) < eps10 .or. Ecinetic(nspin) < eps10 ) then
    Ecomp = .true.
  else
    Ecomp = .false.
  endif

  if( Full_potential ) then
    lfin = 0
  else
    lfin = lmax
  endif

  do L = 0,lfin

    if( Hubb_a .and. L == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = ( lmax + 1 )**2
      nlm2 = nlm1
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*L + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*L + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,L,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

! Integrale radiale pour la regle d'or de Fermi
    call radial_matrix(.true.,1,0,0,iseuil,L,n_comp,nlm1,nlm2,nbseuil, &
           1,nlmam,nlmam2,nr,NRIXS,nrm,nspino,nspinp,psii,r,r_or_bess,Radial_comp,Rmtsd,rof,ui,ur,Vecond)

! Recopie
    if( Full_potential ) then
      zet(1:nr,1:nlm1,1:nlm2,:,:,initl) = ur(1:nr,1:nlm1,1:nlm2,:,:)
    else
      lm0 = L**2 + L + 1
      do m = -L,L
        if( nlm1 == 1 ) then
          n = 1
        else
          n = L + 1 + m
        endif
        do mp = -L,L
          if( nlm2 == 1 .and. mp /= m ) cycle
          np = min(L + m + 1, nlm2 )
          lmp = min(lm0 + mp, nlmam2 )
          zet(1:nr,lm0+m,lmp,:,:,initl) = ur(1:nr,n,np,:,:)
        end do
      end do
    endif

    deallocate( Tau )

    deallocate( ui, ur )

  end do   ! fin de la boucle sur L

  if( nlm2 == 1 ) then
    rof_ph(:,:,:,initl) = rof(:,1,:,:,0,1)
  else
    do iso = 1,nspino
      do isp = 1,nspinp
        do lm = 1,nlmam2
          do i_comp = 1,n_comp
            if( i_comp == 1 ) then
              rof_ph(lm,isp,iso,initl) =  sum( rof(1:nlmam,lm,isp,iso,0,1) )
            else
              rof_ph(lm+nlmam,isp,iso,initl) =  sum( rof(nlmam+1:2*nlmam,lm,isp,iso,0,1) )
            endif
          end do
        end do
      end do
    end do
  endif        

  deallocate( V )

  return
  110 format(/' ---- Radial_wave --',100('-'))
  120 format(/' Energ =',f10.3,' eV,  Enervide =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format(' V0bd  =',2f10.3)
  150 format(' konde =',2f12.5)
end

!***********************************************************************

! Calcul de Integrale_sur_r_et_r' de
! conj(phi(r) * f_reg( min(r,r') ) * f_irg( max(r,r') ) * phi(r') * dr * dr' )
! = Integrale_sur_r ( phi(r) * f_irg(r) * Integrale_sur_r'_de_0_a_r ( f_reg(r') * phi(r') * dr' )
!                   + phi(r) * f_reg(r) * Integrale_sur_r'_de_r_a_Rmax ( f_irg(r') * phi(r') * dr' ) * dr
!
! f_reg = r * solution reguliere
! f_irg = r * solution irreguliere
! phi = r * fonction initiale * r^p    ou   r * solution reguliere

function integr_sing(n,phi1,phi2,f_reg,f_irg,Rmtsd,r,icheck)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db):: integr_sing

  complex(kind=db), dimension(n):: f_irg, f_reg, fct, phi_irg, phi_reg, s_phi_irg, s_phi_reg
  integer icheck
  real(kind=db), dimension(n):: fct_r, phi1, phi2,  r


  phi_reg(:) = phi1(:) * f_reg(:)
  call ffintegr2(s_phi_reg,phi_reg,r,n,1,Rmtsd)

  phi_irg(:) = phi2(:) * f_irg(:)
  call ffintegr2(s_phi_irg,phi_irg,r,n,-1,Rmtsd)

  fct(:) = phi_irg(:) * s_phi_reg(:) + phi_reg(:) * s_phi_irg(:)

  fct_r(:) = real( fct(:),db )
  fr = f_integr3(r,fct_r,1,n,Rmtsd)

  fct_r(:) = aimag( fct(:) )
  fi = f_integr3(r,fct_r,1,n,Rmtsd)

  integr_sing = - cmplx(fr,fi,db)

  if( icheck > 3 ) then
    write(3,110)
    do ir = 1,n
      write(3,120) r(ir)*bohr, f_reg(ir), phi_reg(ir),s_phi_reg(ir), phi_irg(ir), s_phi_irg(ir), fct(ir)
      if( r(ir) > Rmtsd ) exit
    end do
    write(3,130) integr_sing
  endif

  return
  110 format(/5x,'Radius',15x,'f_reg',23x,'phi_reg',19x,'s_phi_reg',19x, 'phi_irg',21x,'s_phi_irg',22x,'fct')
  120 format(1p,13e14.6)
  130 format(/' Integr_sing =',1p,2e14.6)
end

!***********************************************************************

! Calcul L'integrale de 0 a r (is=1) ou r a Rmtsd (is=-1) de fct
! Cas complexe

subroutine ffintegr2(fint,fct,r,n,is,Rmtsd)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db):: a, a_tiers, b, b_demi, c, f0, fm, fp
  complex(kind=db), dimension(2):: dintegr
  complex(kind=db), dimension(n):: fct, fint

  real(kind=db), dimension(n):: r

  tiers = 1._db / 3._db

  if( is == 1 ) then
    i1 = 1
    i2 = n - 1
  else
    i1 = n - 1
    i2 = 1
    fint(n) = (0._db, 0._db)
  endif

  do i = i1,i2,is
    if( i == 1 ) then
      rm = r(i)
      r0 = r(i+1)
      rp = r(i+2)
      fm = fct(i)
      f0 = fct(i+1)
      fp = fct(i+2)
      xm = 0._db
      x0 = rm
      xp = 0.5_db * ( rm + r0 )
    else
      rm = r(i-1)
      r0 = r(i)
      rp = r(i+1)
      fm = fct(i-1)
      f0 = fct(i)
      fp = fct(i+1)
      xm = 0.5_db * ( rm + r0 )
      x0 = r0
      xp = 0.5_db * ( r0 + rp )
      if( r(i) > Rmtsd ) then
        rm = r(i-2); r0 = r(i-1); rp = r(i)
        fm = fct(i-2); f0 = fct(i-1); fp = fct(i)
      endif
   endif

    if( is == 1 .and. r0 > Rmtsd ) then
      if( r(i-1) > Rmtsd ) then
        fint(i) = fint(i-1)
        cycle
      elseif( rm > Rmtsd ) then
        fint(i) = fint(i-1) + dintegr(2)
        cycle
      else
        x0 = Rmtsd
      endif
    endif
    if( is == - 1 .and. r0 > Rmtsd ) then
      fint(i) = (0._db, 0._db)
      cycle
    endif
    if( xp > Rmtsd ) xp = Rmtsd

    a = ( fm * ( rp - r0 ) - f0 * ( rp - rm ) + fp * ( r0 - rm ) ) / ( ( r0 - rm ) * ( rp - r0 ) * ( rp - rm ) )
    b = ( f0 - fm ) / ( r0 - rm ) - a * ( r0 + rm )
    c = f0 - a * r0**2 - b * r0

    a_tiers = a * tiers
    b_demi = b * 0.5_db

    if( is == 1 ) then
      dintegr(1) = ( a_tiers * ( xm**2 + xm * x0  + x0**2 ) + b_demi * ( xm + x0 ) + c ) * ( x0 - xm )
    else
      dintegr(1) = ( a_tiers * ( x0**2 + x0 * xp  + xp**2 ) + b_demi * ( x0 + xp ) + c ) * ( xp - x0 )
    endif

    if( i == i1 ) then
      fint(i) = dintegr(1)
    else
      fint(i) = fint(i-is) + sum( dintegr(:) )
    endif

    if( is == 1 ) then
      dintegr(2) = ( a_tiers * ( x0**2 + x0 * xp  + xp**2 ) + b_demi * ( x0 + xp ) + c ) * ( xp - x0 )
    else
      dintegr(2) = ( a_tiers * ( xm**2 + xm * x0  + x0**2 ) + b_demi * ( xm + x0 ) + c ) * ( x0 - xm )
    endif

  end do

  if( is == 1 ) fint(n) = fint(n-1)

  return
end

!***********************************************************************

! Calculation of the integral from 0 to r (is=1) or from r to Rmtsd (is=-1) of fct

subroutine ffintegr2_r(fint,fct,r,n,is,Rmtsd)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db):: a, a_tiers, b, b_demi, c, f0, fm, fp
  real(kind=db), dimension(2):: dintegr
  real(kind=db), dimension(n):: fct, fint

  real(kind=db), dimension(n):: r

  tiers = 1._db / 3._db

  dintegr(:) = 0._db

  if( is == 1 ) then
    i1 = 1
    i2 = n - 1
  else
    i1 = n - 1
    i2 = 1
    fint(n) = 0._db
  endif

  do i = i1,i2,is
    if( i == 1 ) then
      rm = r(i)
      r0 = r(i+1)
      rp = r(i+2)
      fm = fct(i)
      f0 = fct(i+1)
      fp = fct(i+2)
      xm = 0._db
      x0 = rm
      xp = 0.5_db * ( rm + r0 )
    else
      rm = r(i-1)
      r0 = r(i)
      rp = r(i+1)
      fm = fct(i-1)
      f0 = fct(i)
      fp = fct(i+1)
      xm = 0.5_db * ( rm + r0 )
      x0 = r0
      xp = 0.5_db * ( r0 + rp )
    endif

    if( is == 1 .and. r0 > Rmtsd ) then
      if( r(i-1) > Rmtsd ) then
        fint(i) = fint(i-1)
        cycle
      elseif( rm > Rmtsd ) then
        fint(i) = fint(i-1) + dintegr(2)
        cycle
      else
        x0 = Rmtsd
      endif
    endif
    if( is == - 1 .and. r0 > Rmtsd ) then
      fint(i) = 0._db
      cycle
    endif
    if( xp > Rmtsd ) xp = Rmtsd

    a = ( fm * ( rp - r0 ) - f0 * ( rp - rm ) + fp * ( r0 - rm ) ) / ( ( r0 - rm ) * ( rp - r0 ) * ( rp - rm ) )
    b = ( f0 - fm ) / ( r0 - rm ) - a * ( r0 + rm )
    c = f0 - a * r0**2 - b * r0

    a_tiers = a * tiers
    b_demi = b * 0.5_db

    if( is == 1 ) then
      dintegr(1) = ( a_tiers * ( xm**2 + xm * x0  + x0**2 ) + b_demi * ( xm + x0 ) + c ) * ( x0 - xm )
    else
      dintegr(1) = ( a_tiers * ( x0**2 + x0 * xp  + xp**2 ) + b_demi * ( x0 + xp ) + c ) * ( xp - x0 )
    endif

    if( i == i1 ) then
      fint(i) = dintegr(1)
    else
      fint(i) = fint(i-is) + sum( dintegr(:) )
    endif

    if( is == 1 ) then
      dintegr(2) = ( a_tiers * ( x0**2 + x0 * xp  + xp**2 ) + b_demi * ( x0 + xp ) + c ) * ( xp - x0 )
    else
      dintegr(2) = ( a_tiers * ( xm**2 + xm * x0  + x0**2 ) + b_demi * ( xm + x0 ) + c ) * ( x0 - xm )
    endif

  end do

  if( is == 1 ) fint(n) = fint(n-1)

  return
end

!***********************************************************************

! Transformation harmo comp to Harmo real when Comp_to_real = .true.
!                                 and reverse when false.

subroutine Trans_Tau_ab(Comp_to_real,lmaxa,lmaxb,nspin,Spinorbite,Taull)

  use declarations
  implicit none

  integer:: is, isp1, isp2, la, lb, lm0a, lm0b, lmax, lmaxa, lmaxb, m, mp, nspin

  complex(kind=db), dimension(:,:), allocatable:: Tau_ab, Trans, Trans_i, Trans_a, Trans_b
  complex(kind=db), dimension((lmaxa+1)**2,nspin,(lmaxb+1)**2,nspin):: Taull

  logical:: Comp_to_real, Spinorbite

  lmax = max( lmaxa, lmaxb )
  allocate( Trans(-lmax:lmax,-lmax:lmax) )
  allocate( Trans_i(-lmax:lmax,-lmax:lmax) )
  Trans(:,:) = (0._db, 0._db)
  Trans_i(:,:) = (0._db, 0._db)

! Matrix real to complex  
  do m = -lmax,lmax
    if( m > 0 ) then
      is = (-1)**m
      Trans(m,m) = cmplx( is * sqrt_1o2, 0._db, db)
      Trans(m,-m) = cmplx( 0._db, is * sqrt_1o2, db)
    elseif( m == 0 ) then
      Trans(0,0) = (1._db, 0._db)
    else
      Trans(m,-m) = cmplx(sqrt_1o2, 0._db, db)
      Trans(m,m) = cmplx(0._db, -sqrt_1o2, db)
    endif
  end do
  
  if( .not. Comp_to_real ) Trans = Conjg( Transpose( Trans ) )
    
  Trans_i = Conjg( Transpose( Trans ) )
      
  do la = 0,lmaxa

    lm0a = la**2 + la + 1
    allocate( Trans_a(-la:la,-la:la) )
    Trans_a(-la:la,-la:la) = Trans_i(-la:la,-la:la)

    do lb = 0,lmaxb

      lm0b = lb**2 + lb + 1
      allocate( Trans_b(-lb:lb,-lb:lb) )
      Trans_b(-lb:lb,-lb:lb) = Trans(-lb:lb,-lb:lb)

      allocate( Tau_ab(-la:la,-lb:lb) )

      do isp1 = 1,nspin
        do isp2 = 1,nspin
          if( .not. Spinorbite .and. isp1 /= isp2 ) cycle

          Tau_ab(-la:la,-lb:lb) = Taull(lm0a-la:lm0a+la,isp1,lm0b-lb:lm0b+lb,isp2)

          Tau_ab = matmul( Trans_a, matmul( Tau_ab, Trans_b ) )

! in order it works, I do not understand why
          do m = -la,la
            do mp = -lb,lb
              if( ( m < 0 .and. mp >= 0 ) .or. ( mp < 0 .and. m >= 0 ) ) Tau_ab(m,mp) = - Tau_ab(m,mp)   
            end do
          end do

          Taull(lm0a-la:lm0a+la,isp1,lm0b-lb:lm0b+lb,isp2) = Tau_ab(-la:la,-lb:lb) 

        end do
      end do

      deallocate( Tau_ab, Trans_b )
    
    end do  

    deallocate( Trans_a )

  end do

  return
end

!***********************************************************************

! Transformation real (or tesseral) harmonics basis to cubic harmic basis

subroutine Trans_Tau_ZtoK(lmaxa,lmaxb,nspinp,Spinorbite,Taull)

  use declarations
  implicit none

  integer, parameter:: lmax_ZtoK = 4
   
  integer:: isp1, isp2, L, la, lb, lm0a, lm0b, lmax, lmaxa, lmaxb, nspinp

  complex(kind=db), dimension(:,:), allocatable:: Tau_ab
  complex(kind=db), dimension((lmaxa+1)**2,nspinp,(lmaxb+1)**2,nspinp):: Taull

  logical:: Spinorbite

  real(kind=db):: a, b, c, d, e, f
  real(kind=db), dimension(-lmax_ZtoK:lmax_ZtoK,-lmax_ZtoK:lmax_ZtoK,0:lmax_ZtoK):: ZtoK
  real(kind=db), dimension(:,:), allocatable:: Trans_a, Trans_b

  lmax = max( lmaxa, lmaxb )
  lmax = min( lmax, lmax_ZtoK )

  ZtoK(:,:,:) = 0.0_db

  do L = 0,lmax
    select case(L)
      case(0)
        ZtoK(0,0,0) = 1._db
      case(1)
        ZtoK(-1,1,L) = 1._db
        ZtoK(0,-1,L) = 1._db
        ZtoK(1,0,L)  = 1._db
      case(2)
        ZtoK(-2,2,L) = 1._db
        ZtoK(-1,0,L) = 1._db
        ZtoK(0,-1,L) = 1._db
        ZtoK(1,1,L) = 1._db
        ZtoK(2,-2,L) = 1._db
      case(3)
        a = - sqrt(3._db/8)
        b = sqrt(5._db/8)
        ZtoK(-3,-2,L) = 1._db
        ZtoK(-2,1,L) = a;       ZtoK(-2,3,L) = b 
        ZtoK(-2,1,L) = a;       ZtoK(-2,3,L) = b 
        ZtoK(-1,-3,L) = -b;     ZtoK(-1,-1,L) = a 
        ZtoK(0,0,L) = 1._db
        ZtoK(1,1,L) = -b;       ZtoK(1,3,L) = a 
        ZtoK(2,-3,L) = a;       ZtoK(2,-1,L) = b 
        ZtoK(3,2,L) = 1._db
      case(4)
        c = sqrt(7._db/12)
        d = sqrt(5._db/12)
        e = - sqrt(1._db/8)
        f = sqrt(7._db/8)
        ZtoK(-4,0,L) = c;     ZtoK(-4,4,L) = d
        ZtoK(-3,2,L) = 1._db
        ZtoK(-2,0,L) = -d;    ZtoK(-2,4,L) = c 
        ZtoK(-1,-3,L) = e;    ZtoK(-1,-1,L) = -f 
        ZtoK(0,1,L) = f;      ZtoK(0,3,L) = e 
        ZtoK(1,-4,L) = 1._db
        ZtoK(2,-3,L) = f;     ZtoK(2,-1,L) = e 
        ZtoK(3,1,L) = e;     ZtoK(3,3,L) = -f 
        ZtoK(4,-2,L) = 1._db
    end select
  end do  
  
  do la = 0,lmaxa

    lm0a = la**2 + la + 1
    allocate( Trans_a(-la:la,-la:la) )
    Trans_a(-la:la,-la:la) = ZtoK(-la:la,-la:la,la)

    do lb = 0,lmaxb

      lm0b = lb**2 + lb + 1
      allocate( Trans_b(-lb:lb,-lb:lb) )
      Trans_b(-lb:lb,-lb:lb) = ZtoK(-lb:lb,-lb:lb,lb)

      allocate( Tau_ab(-la:la,-lb:lb) )

      do isp1 = 1,nspinp
        do isp2 = 1,nspinp
          if( .not. Spinorbite .and. isp1 /= isp2 ) cycle

          Tau_ab(-la:la,-lb:lb) = Taull(lm0a-la:lm0a+la,isp1,lm0b-lb:lm0b+lb,isp2)

          Tau_ab = matmul( Trans_a, matmul( Tau_ab, Transpose( Trans_b ) ) )

          Taull(lm0a-la:lm0a+la,isp1,lm0b-lb:lm0b+lb,isp2) = Tau_ab(-la:la,-lb:lb) 

        end do
      end do

      deallocate( Tau_ab, Trans_b )
    
    end do  

    deallocate( Trans_a )

  end do

  return
end

!***********************************************************************

! Calculation and writing of the atomic electron density

! Similar to the first part of Cal_state

subroutine Write_Density(Energ,Density_comp,Full_atom,Harm_cubic,iaabsi,iaprotoi,icheck,ie,ie_computer,Int_statedens, &
                itypei,itypepr,lamstdens,lla_state,lla2_state,lmaxat,mpinodes,n_atom_0,n_atom_ind,n_atom_proto,natome,nenerg, &
                nomfich_s,nonexc_g,nrato,nrm,nspin,nspinp,ntype,numat,Rato,Rmtsd,State_all,Statedens)

  use declarations
  implicit none

  integer:: i_self, ia, iaabsi, iapr, icheck, ie, ie_computer, ipr, iprabs, ir, isp, it, L, &
    la, lamstdens, LL, lla_state, lla2_state, lm,lm1, lm2, lma, m, Long, mpinodes, &
    n_atom_0, n_atom_ind, n_atom_proto, natome, nenerg, nr, nrm, nspin, nspinp, ntype, Z

  character(len=Length_name) nomfich_s, nomficht
  character(len=28) mot28

  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, la_ipr, ll_ipr, lmaxat
  integer, dimension(natome):: iaprotoi, itypei

  logical:: Abs_exc, Density_comp, First, Full_atom, Harm_cubic, nonexc_g, Open_file, State_all
  logical, dimension(0:n_atom_proto):: proto_done

  real(kind=db):: de

  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(0:n_atom_proto):: Rmtsd
  real(kind=db), dimension(lla2_state,nspinp,n_atom_0:n_atom_ind):: Int_Statedens
  real(kind=db), dimension(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1):: Statedens
  real(kind=db), dimension(0:lla_state,nspinp,n_atom_0:n_atom_ind):: Int_Statedens_l
  real(kind=db), dimension(0:lla_state,nspinp):: Statedens_l
  real(kind=db), dimension(nspinp,n_atom_0:n_atom_ind):: Int_Statedens_t
  real(kind=db), dimension(0:nrm,0:ntype):: rato

  if( icheck > 2 )  then
    write(3,110)
    write(3,120) Energ(ie)*rydb
  endif

  proto_done(:) = .false.

  if( nenerg == 1 ) then
    de = 1._db
  elseif( ie == 1 ) then
    de = Energ(2) - Energ(1)
  elseif( ie == nenerg ) then
    de = Energ(ie) - Energ(ie-1)
  else
    de = 0.5_db * ( Energ(ie+1) -  Energ(ie-1) )
  endif

  do ipr = 0,n_atom_proto
    if( lamstdens > -1 ) then
      la = min(lamstdens,lmaxat(ipr))
    else
      Z = numat( itypepr(ipr) )
      if( Z < 19 ) then
        la = min(1,lmaxat(ipr))
        LL = min(2,lmaxat(ipr))
      elseif( Z > 18 .and. Z < 55 ) then
        la = min(2,lmaxat(ipr))
        LL = min(3,lmaxat(ipr))
      else
        la = min(3,lmaxat(ipr))
        LL = min(4,lmaxat(ipr))
      endif
    end if
    la_ipr(ipr) = la    ! max number of harmoniques for the writing of the DOS
    ll_ipr(ipr) = LL    
  end do

  First = .true.
  boucle_ia: do ia = 1,natome

    ipr = iaprotoi(ia)
    la = la_ipr(ipr)
    LL = ll_ipr(ipr)

    if( ia == iaabsi ) iprabs = ipr

    if( Full_atom ) then
      iapr = ia
    else
      iapr = ipr
    endif

    it = itypei(ia)
    Z = numat(it)
    if( proto_done(ipr) .and. .not. Full_atom ) cycle boucle_ia
    if( .not. ( ia == iaabsi .or. State_all ) ) cycle
     
! Calculation of the DOS integrated up to the Rmstd ratomic radius

    do ir = 1,nrato(it)
      if( rato(ir,it) > Rmtsd(ipr) + eps10 ) exit
    end do
    nr = ir

    lma = (LL + 1)**2

    Statedens_l(:,:) = 0._db
    do L = 0,la
      lm1 = L**2 + 1
      lm2 = ( L + 1 )**2
      do isp = 1,nspinp
        do lm = lm1,lm2
          Statedens_l(L,isp) = Statedens_l(L,isp) + Statedens(lm,isp,lm,isp,iapr,ie_computer)
        end do
      end do
    end do
    if( nspinp == 1 ) Statedens_l(0:la,1) = 2 * Statedens_l(0:la,1)

! Integral of the DOS
    do isp = 1,nspinp
      do lm = 1,lma
        Int_Statedens(lm,isp,iapr) = Int_Statedens(lm,isp,iapr) + de * Statedens(lm,isp,lm,isp,iapr,ie_computer)
      end do
    end do

    do L = 0,LL
      lm1 = L**2 + 1
      lm2 = ( L + 1 )**2
      do isp = 1,nspinp
        Int_Statedens_l(L,isp,iapr) = sum(Int_Statedens(lm1:lm2,isp,iapr))
      end do
      if( nspinp == 1 ) Int_Statedens_l(L,1,iapr) = 2 * Int_Statedens_l(L,1,iapr)
    end do

    Int_statedens_t(:,iapr) = 0._db
    do isp = 1,nspinp
      Int_Statedens_t(isp,iapr) = Int_Statedens_t(isp,iapr) + sum( Int_Statedens_l(0:la,isp,iapr) )
    end do

    if( icheck > 1 ) then
      write(3,180) iapr
      lm = 0
      do L = 0,LL
        do m = -L,L
          lm = lm + 1
          do isp = 1,nspinp
            write(3,190) L, m, isp, Statedens(lm,isp,lm,isp,iapr,ie_computer), Int_Statedens(lm,isp,iapr)
          end do
        end do
      end do
      write(3,200)
      do L = 0,la
        write(3,210) L, Int_Statedens_l(L,1:nspinp,iapr)
      end do
      write(3,220) Int_Statedens_t(1:nspinp,iapr)
    endif

    Open_file = ie == 1
    Abs_exc = Full_atom .and. ( ia == iaabsi )  .and. .not. nonexc_g

    if( Open_file .and. icheck > 0 ) then
      if( First ) then
        write(3,230)
        write(3,'(/A)') ' Output file name for the projected density of state:'
        write(3,240) 
      endif
      nomficht = ' '
      Long = len_trim(nomfich_s)
      do L = Long,1,-1
        if( nomfich_s(L:L) == '/' .or. nomfich_s(L:L) == '\' ) exit
      end do
      L = L + 1 
      nomficht(1:Long-L+1) = nomfich_s(L:Long+1)
      Long = len_trim(nomficht)
      nomficht(Long+1:Long+3) = '_sd'

      if( Abs_exc ) then
        call ad_number(0,nomficht,Length_name)
      else
        call ad_number(iapr,nomficht,Length_name)
      endif
      if( Abs_exc ) then
        mot28 = '     Excited absorbing atom:'
      elseif( ia == iaabsi ) then
        mot28 = ' Non-excited absorbing atom:'
      else
        mot28 = '                       Atom:'
      endif
      write(3,250) mot28, Z, ia, ipr, nomficht
    endif
    
  ! Writing file common with cal_state in SCF.f90
    call Write_stdens(Abs_exc,.true.,Density_comp,Energ(ie),Harm_cubic,i_self,iapr,ie_computer,la, &
           Int_Statedens_l,lla_state,lla2_state,mpinodes,n_atom_0,n_atom_ind,nomfich_s,nspin,nspinp, &
           Open_file,Statedens,Statedens_l,.false.)

    proto_done(ipr) = .true.
    First = .false.

  end do boucle_ia

  return
  110 format(/' ---- Write_density --------',100('-'))
  120 format(15x,' Energy =',f10.3,' eV')
  180 format(/'  L  m is  Density of states   Integral    ia =',i3)
  190 format(3i3,3f15.7)
  200 format(/'    L   sum_m(Integral)')
  210 format(i5,2f15.7)
  220 format(' Total =',f12.7,f15.7)
  230 format(/120('-'))
  240 format(/' See above, under "Atom_selec", index "ia" : the position of the atom in the symmetrized cluster',/ &
              '            under "Symsite",    index "ipr": the position of the equivalent atom site for the space group',/ )
  250 format(a28,' Z =',i3,', ia =',i3,', ipr =',i3,': ',A) 
end

!***********************************************************************

! Calculation of Bessel and neuman functions.

subroutine cbess(fnorm,z,lmax,lmaxm,bessel)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db):: bessel(0:lmaxm), fnorm, z

  if( abs( z ) < eps10 ) then
    bessel(0) = fnorm
  else
    bessel(0) = fnorm * sin(z) / z
  endif

  if( lmax == 0 ) return

  if( abs( z ) < eps10 ) then

    bessel(1:lmax) = (0._db, 0._db)

  else

    bessel(1) = bessel(0) / z - fnorm * cos(z) / z

    do L = 2,lmax
      L1 = 2*L - 1
      bessel(L) = L1 * bessel(L-1) / z - bessel(L-2)
    end do

  endif

  return
end

!***********************************************************************

! Calcul des fonctions de bessel et neuman (en fait, divisees par z )

subroutine cbessneu(fnorm,z,lmax,lmaxm,bessel,neuman)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db):: bessel(0:lmaxm), fnorm, neuman(0:lmaxm), z

  neuman(0) = - fnorm * cos(z) / z
  bessel(0) = fnorm * sin(z) / z

  if( lmax == 0 ) return

  neuman(1) = neuman(0) / z - bessel(0)
  bessel(1) = bessel(0) / z + neuman(0)

  do L = 2,lmax
    L1 = 2*L - 1
    neuman(L) = L1 * neuman(L-1) / z - neuman(L-2)
    bessel(L) = L1 * bessel(L-1) / z - bessel(L-2)
  end do

  return
end

!***********************************************************************

! Calculation of spherical Hankel function of first or second kind ( or Riccati-Hankel(+) or (-) function, divided by z )
! and its derivative versus r, not z

subroutine Cal_Hankel(d_Hankel,Hankel,konde,lmax,r,Second_kind)

  use declarations
  implicit none

  integer:: L, lmax

  complex(kind=db):: cfac, fnorm, konde, z
  complex(kind=db), dimension(0:lmax):: d_hankel, hankel
  complex(kind=db), dimension(0:lmax+1):: hankel_t

  logical:: Second_kind
  
  real(kind=db):: r

  z = konde * r

  if( Second_kind ) then
    hankel_t(0) = img
  else
    hankel_t(0) = - img
  endif

  if( Second_kind ) then
    hankel_t(1) = - 1._db + img / z
  else
    hankel_t(1) = - 1._db - img / z
  endif

  do L = 2,lmax+1
    hankel_t(L) = (2*L - 1) * hankel_t(L-1) / z - hankel_t(L-2)
  end do

  do L = 0,lmax
    d_hankel(L) = L * hankel_t(L) / z - hankel_t(L+1)
  end do

! For normalization by the squrare root of the density of state in vacuum
  fnorm = sqrt( konde / pi )
  if( Second_kind ) then
    cfac = fnorm * exp( - img * z ) / z
  else
    cfac = fnorm * exp( img * z ) / z
  endif

  hankel(0:lmax) = cfac * hankel_t(0:lmax)

! Multiplication by konde, because it is the derivative versus r, not z
  d_hankel(0:lmax) = cfac * konde * d_hankel(0:lmax)

  return
end

!***********************************************************************

! Calculation of spherical Bessel function ( or Riccati-Bessel function, divided by z )
! and its derivative versus r, not z

subroutine Cal_Bessel(d_Bessel,Bessel,konde,lmax,r)

  use declarations
  implicit none

  integer:: L, lmax

  complex(kind=db):: cfac, fnorm, konde, z
  complex(kind=db), dimension(0:lmax):: d_Bessel, Bessel
  complex(kind=db), dimension(0:lmax+1):: Bessel_t

  real(kind=db):: r

  z = konde * r

  Bessel_t(0) = sin( z )

  Bessel_t(1) = - cos( z ) + Bessel_t(0) / z

  do L = 2,lmax+1
    Bessel_t(L) = (2*L - 1) * Bessel_t(L-1) / z - Bessel_t(L-2)
  end do

  do L = 0,lmax
    d_Bessel(L) = L * Bessel_t(L) / z - Bessel_t(L+1)
  end do

! For normalization by the squrare root of the density of state in vacuum
  fnorm = sqrt( konde / pi )
  cfac = fnorm / z

  Bessel(0:lmax) = cfac * Bessel_t(0:lmax)

! Multiplication by konde, because it is the derivative versus r, not z
  d_Bessel(0:lmax) = cfac * konde * d_Bessel(0:lmax)

  return
end

!***********************************************************************

! Calcul des fonctions de bessel et neuman ( disisees par z )

subroutine cbessneur(fnorm,z,lmax,lmaxm,bessel,neuman)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db):: bessel(0:lmaxm), neuman(0:lmaxm), z

  neuman(0) = - fnorm * cos(z) / z
  bessel(0) = fnorm * sin(z) / z

  if( lmax == 0 ) return
  
  neuman(1) = neuman(0) / z - bessel(0)
  bessel(1) = bessel(0) / z + neuman(0)

  do L = 2,lmax
    L1 = 2*L - 1
    neuman(L) = L1 * neuman(L-1) / z - neuman(L-2)
    bessel(L) = L1 * bessel(L-1) / z - bessel(L-2)
  end do

  return
end

!***********************************************************************

! Calculation of Bessel function ( Riccati-Bessel / z )
! When Bessel < 10^(-7), one takes the origin limit: j_l = z^L / (2*L+1)!!
! because of unstability of the recursive formula.

subroutine cbessel(bess,ip0,lmax,nr,q,r)

  use declarations
  implicit none

  integer:: i, ip0, j, k, L, L1, lmax, nr

  real(kind=db):: fac, q, z_lim
  real(kind=db), dimension(nr):: r, z
  real(kind=db), dimension(nr,0:lmax):: bessel
  real(kind=db), dimension(nr,ip0:lmax):: bess

  z(:) = q * r(:)

  do L = 0,lmax

    select case(L)

      case(0)

        do i = 1,nr
          if( z(i) < 1.e-7_db ) then
            bessel(i,L) = 1._db - z(i)**2 / 6._db
          else
            bessel(i,L) = sin( z(i) ) / z(i)
          endif
        end do

      case(1)

        fac = 1._db / 3
        z_lim = 3.e-7_db  ! = ( 1.e-7_db / fac )**( 1._db / L )
        do i = 1,nr
          if( z(i) < z_lim ) then
            bessel(i,L) = fac * z(i)
          else
            bessel(i,L) = ( bessel(i,0) - cos( z(i) ) ) / z(i)
          endif
        end do

      case default

        k = 1
        do j = 3,2*L+1,2
          k = k * j
        end do
        fac = 1._db / k
        z_lim = ( 1.e-7_db / fac )**( 1._db / L )
        L1 = 2*L - 1
 ! L = 2, fac = 1/15, z_lim = 0.00122
 ! L = 3, fac = 1/105, z_lim = 0.0219
 ! L = 4, fac = 1/945, z_lim = 0.0986
 ! L = 5, fac = 1/10395, z_lim = 0.253

        do i = 1,nr
          if( z(i) < z_lim ) then
            bessel(i,L) = fac * z(i)**L
          else
            bessel(i,L) = L1 * bessel(i,L-1) / z(i) - bessel(i,L-2)
          endif
        end do

    end select

  end do

  do L = ip0,lmax
    bess(:,L) = bessel(:,L)
  end do

  return
end


