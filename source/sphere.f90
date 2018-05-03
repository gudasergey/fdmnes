! FDMNES subroutines
! File containing all concerning the solutions in the atomic sphere

!***********************************************************************

! Resolution of the spherical Schrodinger equation in the atoms

subroutine Sphere(Axe_Atom_grn,Ecinetic,Eimag,Energ,Enervide,Fac_tau,Full_atom, &
            Full_potential,Green,Hubb_a,Hubb_d,iaabsi,iapr,iaprotoi, &
            ibord,icheck,igreq,igroupi,iopsymr,lmax,lmax_pot,m_hubb, n_atom_0, &
            n_atom_ind,n_atom_proto,natome,nbord,nbtm,nbtm_fdm,neqm,ngroup_m,nlm_pot,nlmagm, &
            nlmmax,nphiato1,nphiato7,npsom,nr,nspin,nspino,nspinp,numat,phiato,posi,r,Relativiste,Renorm,Rmtg, &
            Rmtsd,Spinorbite,Tau_ato,V_hubb,V_intmax,V0bd,Vrato,xyz,Ylm_comp,Ylmato)

  use declarations
  implicit none

  integer:: ia, iaabsi, iang, iapr, ib, icheck, iga, igr, isp, kgr, l, &
    l_hubbard, lfin, lm, lm0, lmax, lmax_pot, lmp, m, m_hubb, mp, n, &
    n_atom_0, n_atom_ind, n_atom_proto, natome, nbtm, nbtm_fdm, neqm, ngroup_m, nlm, nlm_pot, &
    nlm1, nlm2, nlmagm, nlmmax, np, nphiato1, nphiato7, npsom, nr, nrmtg, nrmtsd, nspin, nspino, nspinp, numat
  integer, dimension(natome):: iaprotoi, igroupi, nbord
  integer, dimension(nopsm):: iopsymr
  integer, dimension(nbtm,natome):: ibord
  integer, dimension(0:n_atom_proto,neqm):: igreq

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(:,:,:,:), allocatable:: Fac_Wronsk, Tau
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,n_atom_0:n_atom_ind):: Fac_tau, Tau_ato
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  logical:: Ecomp, Full_atom, Full_potential, Green, Hubb_a, Hubb_d, Hubb_m, Radial_comp, Relativiste, Renorm, &
    Spinorbite, Ylm_comp

  real(kind=db):: cosang, Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_grn
  real(kind=db), dimension(nbtm_fdm,nlmmax,natome):: Ylmato

  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nr,nspin):: g0, gm, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(nphiato1,nlmagm,nspinp,nspino,natome,nphiato7):: phiato

  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur

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

  do l = 0,lfin

    if( Hubb_a .and. l == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = nlm   ! = (lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*l + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*l + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )
    allocate( Fac_Wronsk(nlm1,nlm2,nspinp,nspino) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,l,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

! Recopie
    if( Full_potential ) then
      do lm = 1,nlm1
        do lmp = 1,nlm2
          Tau_ato(lm,:,lmp,:,iapr) = Tau(lm,:,lmp,:)
        end do
      end do
    else
      lm0 = l**2 + l + 1
       do m = -l,l
        if( nlm1 == 1 ) then
          n = 1
        else
          n = l + 1 + m
        endif
        do mp = -l,l
          if( nlm1 == 1 ) then
            if( m /= mp ) cycle
            np = 1
          else
            np = l + 1 + mp
          endif
          Tau_ato(lm0+m,:,lm0+mp,:,iapr) = Tau(n,:,np,:)
        end do
      end do
    endif

    deallocate( Tau )

    if( .not. Green ) then

      call Wronsk_fac(Fac_Wronsk,Full_potential,icheck,konde,l,lmax,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,r, &
                  Radial_comp,Rmtg,ui,ur)

      if( Full_potential ) then
        do lm = 1,nlm1
          do lmp = 1,nlm2
            if( nspino == 1 ) then
              do isp = 1,nspinp
                Fac_tau(lm,isp,lmp,isp,iapr) = Fac_Wronsk(lm,lmp,isp,1)
              end do
            else ! not very sure
              Fac_tau(lm,:,lmp,:,iapr) = Fac_Wronsk(lm,lmp,:,:)
            endif
          end do
        end do
      else
        lm0 = l**2 + l + 1
        do m = -l,l
          if( nlm1 == 1 ) then
            n = 1
          else
            n = l + 1 + m
          endif
          do mp = -l,l
            if( nlm1 == 1 ) then
              if( m /= mp ) cycle
              np = 1
            else
              np = l + 1 + mp
            endif
            if( nspino == 1 ) then
              do isp = 1,nspinp
                Fac_tau(lm0+m,isp,lm0+mp,isp,iapr) = Fac_Wronsk(n,np,isp,1)
              end do
            else ! not very sure
              Fac_tau(lm0+m,:,lm0+mp,:,iapr) = Fac_Wronsk(n,1,:,:)
            endif
          end do
        end do
      endif
    endif
    
    deallocate( Fac_Wronsk )

 ! Calcul les fonctions radiales phiato.
    if( .not. Green ) then

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

        call cal_phiato(Full_potential,ia,iang,ibord,icheck,iopsymr,l,lmax,nlm1,nlm2,natome, &
             nbord(ia),nbtm,nbtm_fdm,nlmagm,nlmmax,nphiato1,nphiato7,npsom,nr,nspinp,nspino,phiato,posi,r, &
             Radial_comp,Spinorbite,ui,ur,xyz,Ylm_comp,Ylmato)
      end do
    endif

    deallocate( ui, ur )

  end do   ! end of loop over l

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

  integer:: icheck, ir, ispin, l1, l2, l3, lm, lm1, lm2, lm3, lmax, &
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

! n'a pas d'importance en FDM
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
  do l1 = 0,lmax
    do m1 = -l1,l1
      lm1 = lm1 + 1
      lm2 = 0
      do l2 = 0,lmax
        do m2 = -l2,l2
          lm2 = lm2 + 1
          if( lm1 == lm2 .or. nlm == 1 ) cycle

          lm3 = 1
          do l3 = 1,lmax_pot
            do m3 = -l3,l3
              lm3 = lm3 + 1

              if( Ylm_comp ) then
                g = gauntcp(l1,m1,l2,m2,l3,m3)
              else
                g = Gaunt_r(l1,m1,l2,m2,l3,m3)
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

! alfa_sf = constante de structure fine.
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

! Resolution de l'equation de Schrodinger radiale. Inclus le cas du potentiel non spherique.
! Applele par Sphere, radial, radial_sd
! Dans ur(ir,m,mp,isp,isol) : (mp, isol) definissent les etats de base
!                             (m, isp) sont le moment et le spin reels

subroutine Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,ll,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

  use declarations
  implicit none

  integer:: icheck, im, ip, ir, isol, isp, ispp, isq, isr, l, l_hubbard, l2, li, ll, lf, lmax, m, m_hubb, mp, ms, &
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
  real(kind=db), dimension(:,:), allocatable:: fac_m

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
      li = ll
      lf = ll
    endif

  do l = li,lf

    l2 = l * ( l + 1 )

    if( Hubb_a .and. l == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif
    Hubb_nd = Hubb_m .and. .not. Hubb_d

    ci = - Eimag / ( 4 * l + 6 )

    br = - numat / ( l + 1._db )
    cr(:) = - ( 2 * numat * br + Ecinetic(:) ) / ( 4 * l + 6 )

! Je fais l'hypothese que le terme de Hubbard ne croise pas les solutions 1 et 2 en cas de spinorbite.

! Developpement a l'origine
    do isol = 1,nspino
      do isp = 1,nspinp
        isr = min( isp, nspin )
        do m = -l,l

          if( Full_potential ) then
            n = l**2 + l + 1 + m
          elseif( Hubb_m .or. Spinorbite ) then
            n = l + 1 + m
          else
            n = 1
          endif

          if( (m == l .and. isp == 1) .or. (m == -l .and. isp == 2 ) ) then
            nsol = 1
          else
            nsol = 2
          endif
          if( Spinorbite .and. nsol == 1 .and. isp /= isol ) cycle

          if( Spinorbite .and. Relativiste ) then
            if( nsol == 1  .or. isol == 1 ) then
              p = sqrt( ( l + 1 )**2 - ( alfa_sf * numat )**2 )
            else
              p = sqrt( l**2 - ( alfa_sf * numat )**2 )
            endif
          elseif( Spinorbite ) then
            if( nsol == 1 ) then
              p = l + 1._db
            elseif( isol == 1 ) then
              p = 0.5_db + 0.5_db * sqrt( 1._db + 4*(l**2) + 8*l )
            else
              if( l == 1 ) then
! En fait, il n'y a pas de solution pour l=1 avec spin-orbite non nelativiste !
!                      p = 1._db * l
                p = l + 1._db
              else
                p = 0.5_db + 0.5_db * sqrt( -5._db + 4*(l**2) )
              endif
            endif
          elseif( Relativiste ) then
            p = sqrt( l**2 + l + 1 - ( alfa_sf * numat )**2 )
          else
            p = l + 1._db
          endif

 ! Attention, c'est du point de vue du m de l'orbitale
          if( Spinorbite .and. nsol /= 1 ) then
            if( isol == 1 .and. isp == 1 ) then
              fac = sqrt( ( l + m + 1 ) / ( 2*l + 1._db ) )
            elseif( isol == 1 .and. isp == 2 ) then
              fac = sqrt( ( l - m + 1 ) / ( 2*l + 1._db ) )
            elseif( isol == 2 .and. isp == 1 ) then
              fac = sqrt( ( l - m ) / ( 2*l + 1._db ) )
            else
              fac = - sqrt( ( l + m ) / ( 2*l + 1._db ) )
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
    l = ll
    allocate( fac_m(-l:l,2) )
    do m = -li,li
       fac_m(m,1) = sqrt( ( l - m ) * ( l + m + 1._db ) )
       fac_m(m,2) = sqrt( ( l + m ) * ( l - m + 1._db ) )
    end do
  endif

  do ir = 2,nr-1
    im = ir - 1
    ip = ir + 1

    do isol = 1,nspino
      do isp = 1,nspinp
        isr = min( isp, nspin )
        isq = 3 - isp

        do l = li,lf

          l2 = l * ( l + 1 )

          if( Hubb_a .and. l == l_hubbard( numat ) )  then
            Hubb_m = .true.
          else
            Hubb_m = .false.
          endif

          do m = -l,l

            if( Spino_simple .and. ( ( m == l .and. isp == 1 ) .or. ( m == -l .and. isp == 2 ) ) .and. isp /= isol ) cycle
  
            if( Full_potential ) then
              n = l**2 + l + 1 + m
            elseif( Hubb_m .or. Spinorbite ) then
              n = l + 1 + m
            else
              n = 1
              if( m /= 0 ) cycle
            endif

            td = g0(ir,isr) + l2 * f2(ir)

            if( Spinorbite ) then
              if( isp == 1 ) then
                ms = m
                mp = m + 1   ! m de l'autre spin
              else
                ms = - m
                mp = m - 1
              endif
              if( Full_potential ) then
                fac = sqrt( ( l - ms ) * ( l + ms + 1._db ) )
              else
                fac = fac_m(m,isp)
              endif
              td = td + ms * gso(ir,isp)
            else
              mp = l + 1 ! pour que le if en dessous soit faux
            endif
            if( Full_potential ) then
              np = l**2 + l + 1 + mp
            elseif( Hubb_nd .or. Spinorbite ) then
              np = l + 1 + mp
            else
              np = 1
            endif

! le deuxieme indice "m" est la composante non nulle a l'origine
            ur(ip,n,:,isp,isol) = td * ur(ir,n,:,isp,isol) + gm(ir,isr) * ur(im,n,:,isp,isol)
            if( Radial_comp ) then
              ui(ip,n,:,isp,isol) = td * ui(ir,n,:,isp,isol) + gm(ir,isr) * ui(im,n,:,isp,isol)
              if( Ecomp ) then 
                ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol) + Eimag * ui(ir,n,:,isp,isol)
                ui(ip,n,:,isp,isol) = ui(ip,n,:,isp,isol) - Eimag * ur(ir,n,:,isp,isol)
              endif
            endif

            if( abs(mp) <= l ) then
              ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol) + fac * gso(ir,isp) * ur(ir,np,:,isq,isol)
              if( Radial_comp ) ui(ip,n,:,isp,isol) = ui(ip,n,:,isp,isol) + fac * gso(ir,isp) * ui(ir,np,:,isq,isol)
            endif

            if( Hubb_m ) then
              do mp = -l,l
                if( Hubb_d .and. m /= mp ) cycle
                if( Full_potential ) then
                  np = l**2 + l + 1 + mp
                else
                  np = l + 1 + mp
                endif
                do ispp = 1,nspinp
                  V_h = V_hubb(m,mp,isp,ispp) 
                  if( Radial_comp ) then
                    ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol) + real( V_h, db) * ur(ir,np,:,ispp,isol) &
                                                              - aimag( V_h ) * ui(ir,np,:,ispp,isol)
                    ui(ip,n,:,isp,isol) = ui(ip,n,:,isp,isol) + real( V_h, db) * ui(ir,np,:,ispp,isol) &
                                                              + aimag( V_h ) * ur(ir,np,:,ispp,isol)
                  else
                    ur(ip,n,:,isp,isol) = ur(ip,n,:,isp,isol) + real( V_h, db) * ur(ir,np,:,ispp,isol)
                  endif
                end do
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

  if( Spinorbite .and. .not. Full_potential ) deallocate( fac_m )

  if( icheck > 2 ) call write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspinp,nspino,numat,r,Radial_comp,Rmtg,ui,ur,1)

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

  call Renormal(Radial_comp,Full_potential,icheck,konde,ll,lmax,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,r,Rmtg, &
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
    call write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspinp,nspino,numat,r,Radial_comp,Rmtg,ui,ur,2)
    do ir = 1,nr
      ur(ir,:,:,:,:) = ur(ir,:,:,:,:) / r(ir)
      if( Radial_comp ) ui(ir,:,:,:,:) = ui(ir,:,:,:,:) / r(ir)
    end do
  endif

  return
end

!***********************************************************************

! Ecriture de la fonction radiale

subroutine write_ur(Full_potential,li,lf,nlm1,nlm2,nr,nspinp,nspinpo,numat,r,Radial_comp,Rmtg,ui,ur,icom)

  use declarations
  implicit none

  integer:: i, icom, ir, isp, isol, k, kmax, l, li, lf, lp, m, nlm1, nlm2, mp, n, np, nr, nspinp, nspinpo, numat

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
  do l = li,lf

    if( nlm2 == 1 ) then
      do m = -l,l
        if( nlm1 == 1 .and. m /= 0 ) cycle
        if( nlm1 == 1 ) then
          if( icom == 1 ) then
            write(3,110) numat, l
          elseif( icom == 2 ) then
            write(3,120) numat, l
          else
            write(3,130) numat, l
          endif
        else
          if( icom == 1 ) then
            write(3,140) numat, l, m
          elseif( icom == 2 ) then
            write(3,150) numat, l, m
          else
            write(3,160) numat, l, m
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
      do mp = -l,l
        np = np + 1
        do isol = 1,nspinpo
          if( nspinpo == 1 ) then
            if( icom == 1 ) then
              write(3,140) numat, l, mp
            elseif( icom == 2 ) then
              write(3,150) numat, l, mp
            else
              write(3,160) numat, l, mp
            endif
          else
            if( icom == 1 ) then
              write(3,210) numat, l, mp, isol
            elseif( icom == 2 ) then
              write(3,220) numat, l, mp, isol
            else
              write(3,230) numat, l, mp, isol
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
  110 format(/' Radial wave function time r:  Z =',i3,', l =',i2,/)
  120 format(/' Radial wave function, time r, after normalization:', '  Z =',i3,', l =',i2,/)
  130 format(/' Radial singular function time r:  Z =',i3,', l =',i2,/)
  140 format(/' Radial wave function time r:  Z =',i3,', l =',i2, ', m =',i2,/)
  150 format(/' Radial wave function, time r, after normalization:', '  Z =',i3,', l =',i2,', m =',i2,/)
  160 format(/' Radial singular function time r:  Z =',i3,', l =',i2, ', m =',i2,/)
  170 format(37x,' Solution 1',46x,' Solution 2')
  180 format(23x,' Solution 1',18x,' Solution 2')
  190 format('     Radius    ',392a14)
  200 format(1p,491e14.6)
  210 format(/' Radial wave function time r:  Z =',i3,', l =',i2, ', m =',i2,', Solution',i2,/)
  220 format(/' Radial wave function, time r, after normalization:', '  Z =',i3,', l =',i2,', m =',i2,', Solution',i2,/)
  230 format(/' Radial singular function time r:  Z =',i3,', l =',i2, ', m =',i2,', Solution',i2,/)
end

!***********************************************************************

! Normalisation of radial functions using continuity at muffin-tin radius with solutions in vaccum
! Called by Sch_radial

subroutine Renormal(Radial_comp,Full_potential,icheck,konde,ll,lmax,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,r, &
                  Rmtg,Tau,ui,ur,Failed)

  use declarations
  implicit none

  integer icheck, i, inr, ir, isol, isp, isr, l, ll, lmax, lp, m, nlm1, nlm2, mp, n, np, nr, nrmtg, nspin, nspino, nspinp

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
  do l = 0,lmax
    if( .not. Full_potential .and. l /= ll ) cycle
    s1(:) = - img * hs(1,l,:)
    s2(:) = - img * hs(2,l,:)
    do mp = -l,l
      if( nlm2 == 1 .and. mp > -l ) exit
      np = np + 1
      do isp = 1,nspinp
        isr = min( isp, nspin )
        if( nlm2 == 1 ) then
          Wronsks(:,1,isp,:) = u1(:,1,isp,:) * s2(isr) - u2(:,1,isp,:) * s1(isr)
          Wronske(:,1,isp,:) = u1(:,1,isp,:) * bs(2,l,isr) - u2(:,1,isp,:) * bs(1,l,isr)
        else
! Dans ce cas le deuxieme indice de "u" est celui de l'attaque
! c'est aussi le premier indice de Wronsk
! et l est celui sortant
          Wronsks(:,np,isp,:) = u1(np,:,isp,:) * s2(isr) - u2(np,:,isp,:) * s1(isr)
          Wronske(:,np,isp,:) = u1(np,:,isp,:) * bs(2,l,isr) - u2(np,:,isp,:) * bs(1,l,isr)
        endif
      end do
    end do
  end do

  if( icheck > 2 ) then
    n = 0
    do l = 0,lmax
      if( .not. Full_potential .and. l /= ll ) cycle
      write(3,110) l
      write(3,120) - img * hs(1,l,:)
      write(3,130) - img * hs(2,l,:)
      write(3,140) bs(1,l,:)
      write(3,150) bs(2,l,:)
      write(3,160) Wronskout
      if( nspino == 2 ) then
        write(3,162)
      elseif( nspinp == 2 ) then
        write(3,164)
      else
        write(3,166)
      endif
      do m = -l,l
        if( nlm1 == 1 .and. m /= 0 ) cycle
        n = n + 1
        np = 0
        do lp = 0,lmax
          if( .not. Full_potential .and. lp /= ll ) cycle
          do mp = -lp,lp
            if( nlm2 == 1 .and. mp /= 0 ) cycle
            np = np + 1
            write(3,170) l, m, lp, mp, ' Wronsks   =', ( Wronsks(n,np,:,isol), isol = 1,nspino )
            write(3,170) l, m, lp, mp, ' u1        =', ( u1(n,np,:,isol), isol = 1,nspino )
            write(3,170) l, m, lp, mp, ' u2        =', ( u2(n,np,:,isol), isol = 1,nspino )
          end do
        end do
      end do
    end do
  endif

  call cal_ampl(Ampl,Full_potential,icheck,ll,lmax,nlm1,nlm2,nspino,nspinp,Tau,Wronske,Wronsks,Wronskout,Failed)

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
        do m = -ll,ll
          n = n + 1
          do isol = 1,nspino ! le spin d'attaque devient l'indice de solution

            mp = m + isol - isp
            if( mp > ll .or. mp < -ll ) cycle
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
          do isol = 1,nspino ! le spin d'attaque devient l'indice de solution
            do np = 1,nlm2

! np et n sont bien dans ce sens !
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
  110 format(/' Renormalisation parameters, l =',i2)
  120 format(/' -i.hankel =',1p,8e11.3)
  130 format('     deriv =',1p,8e11.3)
  140 format('    bessel =',1p,8e11.3)
  150 format('     deriv =',1p,8e11.3)
  160 format(' Wronskout =',1p,8e11.3)
  162 format('  l  m lp mp',21x,'up up',17x,'dn up',17x,'up dn',17x, 'dn dn')
  164 format('  l  m lp mp',  23x,'up',19x,'dn')
  166 format('  l  m lp mp')
  170 format(4i3,a12,1p,16e11.3)
end

!***********************************************************************

! Calcul des amplitudes des fonctions radiales

subroutine cal_ampl(Ampl,Full_potential,icheck,ll,lmax,nlm1,nlm2,nspino,nspinp,Tau,Wronske,Wronsks,Wronskout,Failed)

  use declarations
  implicit none

  integer:: i, icheck, isol, isp, ispp, j, l, ll, lmax, lp, m, nlm1, nlm2, mp, mq, n, nd, ndim, np, nspino, nspinp, nu

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
! Pour Ampl, isp est le spin d'attaque (donc l'indice de solution apres la renormalisation)
! Le m de Amp est m + 1/2 - is = m + isol - 1
      do m = -ll-1,ll
        if( m == -ll-1 ) then
          isp = 2
          isol = 2
          mq = m + isp - 1
          n = ll + 1 + mq
          Ampl(n,1,isp,isol) = Wronskout / Wronske(n,1,isp,isol)
          Tau(n,isp,n,isp) = - Wronske(n,1,isp,isol) / Wronsks(n,1,isp,isol)
        elseif( m == ll ) then
          isp = 1
          isol = 1
          mq = m + isp - 1
          n = ll + 1 + mq
          Ampl(n,1,isp,isol) = Wronskout / Wronske(n,1,isp,isol)
          Tau(n,isp,n,isp) = - Wronske(n,1,isp,isol) / Wronsks(n,1,isp,isol)
        else
          nu = ll + 1 + m
          nd = ll + 1 + m + 1
          Det = Wronsks(nu,1,1,1) * Wronsks(nd,1,2,2) - Wronsks(nd,1,2,1) * Wronsks(nu,1,1,2)
          Amp(nu,1,1,1) = Wronsks(nd,1,2,2) / Det
          Amp(nu,1,1,2) = - Wronsks(nd,1,2,1) / Det
          Amp(nd,1,2,1) = - Wronsks(nu,1,1,2) / Det
          Amp(nd,1,2,2) = Wronsks(nu,1,1,1) / Det

          do isp = 1,nspinp
            n = ll + 1 + m + isp - 1
            do ispp = 1,nspinp
              np = ll + 1 + m + ispp - 1
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
      do l = 0,lmax
        if( .not. Full_potential .and. l /= ll ) cycle
        do m = -l,l
          n = n + 1
          i = i + 1

          j = 0
          do ispp = 1,nspinp  ! solution
            if( nspino == 1 .and. ispp /= isp ) cycle
            isol = min(nspino,ispp)
            np = 0
            do lp = 0,lmax
              if( .not. Full_potential .and. lp /= ll ) cycle
              do mp = -lp,lp
                np = np + 1
                j = j + 1
! colonne harmonique reelle, ligne harmonique attaque
!                Mat_e(i,j) = Wronske(n,np,isp,isol)
!                Mat_s(i,j) = Wronsks(n,np,isp,isol)
                Mat_e(j,i) = Wronske(np,n,isp,isol)
                Mat_s(j,i) = Wronsks(np,n,isp,isol)
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

      Mat_A = Mat_s

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
        do l = 0,lmax
          if( .not. Full_potential .and. l /= ll ) cycle
          do m = -l,l
            n = n + 1
            i = i + 1

            j = 0
            do isol = 1,nspino
              np = 0
              do lp = 0,lmax
                if( .not. Full_potential .and. lp /= ll ) cycle
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
    if( ll == 0 .or. icheck > 1 ) then
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
      do l = 0,lmax
        if( .not. Full_potential .and. l /= ll ) cycle
        do m = -l,l
          if( nlm1 == 1 .and. m /= 0 ) cycle
          n = n + 1
          if( nlm1 == 1 .and. nspinp == 1 ) then
            write(3,220) l, Ampl(n,1,isp,min(isp,nspino))
          elseif( nlm1 == 1 ) then
            write(3,230) l, isp, Ampl(n,1,isp,min(isp,nspino))
          elseif( nlm2 == 1 .and. nspinp == 1 ) then
            write(3,230) l, m, Ampl(n,1,isp,:)
          elseif( nlm2 == 1 ) then
            write(3,240) l, m, isp, Ampl(n,1,isp,:)
          elseif( nspinp == 1 ) then
            write(3,230) l, m, ( Ampl(n,:,isp,ispp), ispp = 1,nspino )
          else
            write(3,240) l, m, isp, ( Ampl(n,:,isp,ispp), ispp = 1,nspino )
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
      do l = 0,lmax
        if( .not. Full_potential .and. l /= ll ) cycle
        do m = -l,l
          if( nlm1 == 1 .and. m /= 0 ) cycle
          n = n + 1
          if( nlm1 == 1 .and. nspinp == 1 ) then
            write(3,220) l, Tau(n,isp,n,isp)
          elseif( nlm1 == 1 ) then
            write(3,230) l, isp, Tau(n,isp,n,isp)
          elseif( nlm2 == 1 .and. nspinp == 1 ) then
            write(3,230) l, m, Tau(n,isp,n,isp)
          elseif( nlm2 == 1 .and. nspino == 1 ) then
            write(3,240) l, m, isp, Tau(n,isp,n,isp)
          else
            write(3,240) l, m, isp, ( Tau(n,isp,:,ispp), ispp = 1,nspinp )
          endif
        end do
        if( nlm2 == 1 .and. nspino == 1 ) exit
      end do
    end do
  endif

! Test sur le calcul
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
  140 format(/' Amplitude',/'  l',11x,'Ampl')
  150 format(/' Amplitude',/'  l isp',9x,'Ampl')
  155 format(/' Amplitude',/'  l  m isp',6x,'Ampl(Sol 1)',12x, 'Ampl(Sol 2)')
  160 format(/' Amplitude',/'  l  m',10x,'Ampl')
  170 format(/' Amplitude',/'  l  m isp',10x,'Ampl')
  180 format(/'  l',11x,'Tau')
  190 format(/'  l isp',10x,'Tau')
  200 format(/'  l  m',11x,'Tau')
  210 format(/'  l  m isp',10x,'Tau')
  220 format(i3,1p,16(1x,2e11.3))
  230 format(2i3,1p,16(1x,2e11.3))
  240 format(3i3,1p,16(1x,2e11.3))
  250 format(/'  Off-digonal Tau higher than diagonal for', ' lm =',i3,', lmp =',i3,'. Off-diagonal potential neglected !')
end

!***********************************************************************

! Calculation of Fac_tau = 0.5 * i * Wronskien(b,-i hankel^(-)) / Wronskien(b,bessel)
! Used by Mat when Eimag not zero (or FDM_comp_m true). Under construction. 

subroutine Wronsk_fac(Fac_Wronsk,Full_potential,icheck,konde,ll,lmax,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,r, &
                  Radial_comp,Rmtg,ui,ur)

  use declarations
  implicit none

  integer icheck, inr, isol, isp, isr, l, ll, lmax, lp, m, nlm1, nlm2, mp, n, np, nr, nrmtg, nspin, nspino, nspinp

  complex(kind=db), dimension(0:lmax):: bess, d_bess, d_hank, hank
  complex(kind=db), dimension(nspin):: e1, e2, konde, s1, s2
  complex(kind=db), dimension(nlm1,nlm2,nspinp,nspino):: Fac_Wronsk, u1, u2, Wronske, Wronsks
  complex(kind=db), dimension(2,0:lmax,nspin):: bs, hs

  logical:: Radial_comp, Full_potential

  real(kind=db):: a, ai, b, bi, c, ci, d0, d01, d02, d1, d2, d12, Rmtg
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
    call Cal_Hankel(d_hank,hank,konde(isp),lmax,Rmtg,.true.)
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

  np = 0
  do l = 0,lmax
    if( .not. Full_potential .and. l /= ll ) cycle
    e1(:) = bs(1,l,:)
    e2(:) = bs(2,l,:)
    s1(:) = - img * hs(1,l,:)
    s2(:) = - img * hs(2,l,:)
    do mp = -l,l
      if( nlm2 == 1 .and. mp > -l ) exit
      np = np + 1
      do isp = 1,nspinp
        isr = min( isp, nspin )
        if( nlm2 == 1 ) then
          Wronske(:,1,isp,:) = u1(:,1,isp,:) * e2(isr) - u2(:,1,isp,:) * e1(isr)
          Wronsks(:,1,isp,:) = u1(:,1,isp,:) * s2(isr) - u2(:,1,isp,:) * s1(isr)
        else
! Dans ce cas le deuxieme indice de "u" est celui de l'attaque
! c'est aussi le premier indice de Wronsk
! et l est celui sortant
          Wronske(:,np,isp,:) = u1(np,:,isp,:) * e2(isr) - u2(np,:,isp,:) * e1(isr)
          Wronsks(:,np,isp,:) = u1(np,:,isp,:) * s2(isr) - u2(np,:,isp,:) * s1(isr)
        endif
      end do
    end do
  end do

  do isp = 1,nspinp 
    Fac_Wronsk(:,:,isp,:) = 0.5_db * img * Wronsks(:,:,isp,:) / Wronske(:,:,isp,:)
  end do 
  
  if( icheck > 2 ) then
    n = 0
    do l = 0,lmax
      if( .not. Full_potential .and. l /= ll ) cycle
      write(3,110) l
      write(3,120) hs(1,l,:)
      write(3,130) hs(2,l,:)
      if( nspino == 2 ) then
        write(3,162)
      elseif( nspinp == 2 ) then
        write(3,164)
      else
        write(3,166)
      endif
      do m = -l,l
        if( nlm1 == 1 .and. m /= 0 ) cycle
        n = n + 1
        np = 0
        do lp = 0,lmax
          if( .not. Full_potential .and. lp /= ll ) cycle
          do mp = -lp,lp
            if( nlm2 == 1 .and. mp /= 0 ) cycle
            np = np + 1
            write(3,170) l, m, lp, mp, ' Fac_Wronsk =', ( Fac_Wronsk(n,np,:,isol), isol = 1,nspino )
            write(3,170) l, m, lp, mp, ' Wronske    =', ( Wronske(n,np,:,isol), isol = 1,nspino )
            write(3,170) l, m, lp, mp, ' Wronsks    =', ( Wronsks(n,np,:,isol), isol = 1,nspino )
            write(3,170) l, m, lp, mp, ' u1         =', ( u1(n,np,:,isol), isol = 1,nspino )
            write(3,170) l, m, lp, mp, ' u2         =', ( u2(n,np,:,isol), isol = 1,nspino )
          end do
        end do
      end do
    end do
  endif

  return
  110 format(/' Wronskian calculation, l =',i2)
  120 format(/' hankel(-) =',1p,8e15.7)
  130 format( '     deriv =',1p,8e15.7)
  162 format('  l  m lp mp',21x,'up up',17x,'dn up',17x,'up dn',17x, 'dn dn')
  164 format('  l  m lp mp',  23x,'up',19x,'dn')
  166 format('  l  m lp mp')
  170 format(4i3,a12,1p,16e15.7)
end

!**********************************************************************

! Calculation of radial integrals and irregular solution
! Called by Tenseur_car and S_nrixs_cal

subroutine radial(Classic_irreg,Ecinetic,Eimag,Energ,Enervide,Eseuil,Final_tddft,Full_potential,Hubb_a,Hubb_d,icheck, &
         initlv,ip_max,ip0,lmax,lmax_pot,m_hubb,nbseuil,ninit1,ninitlv,nlm_pot,nlma,nlma2,nr,NRIXS,nrm,nspin,nspino, &
         nspinp,numat,psii,r,r_or_bess,Relativiste,Renorm,Rmtg,Rmtsd,rof,Singul,Solsing,Spinorbite,V_hubb,V_intmax,V0bd,Vrato, &
         Ylm_comp)

  use declarations
  implicit none

  integer:: i, icheck, initl, initlv, ip, ip_max, ip0, ip1, ip2, iseuil, isp, l, l_hubbard, lfin, lm, lmax, lmax_pot, &
    lmp, lp, m, m_hubb, nlm1, nlm2, mp, nbseuil, ninit1, ninitlv, nlm, nlm_pot, nlma, nlma2, nr, nrm, &
    nrmtsd, nrmtg, nspin, nspino, nspinp, numat

  character(len=2):: mot1, mot2
  character(len=104):: mot

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlma,nspinp,ip0:ip_max,ip0:ip_max,ninitlv):: Singul
  complex(kind=db), dimension(nlma,nlma2,nspinp,nspino,ip0:ip_max,ninitlv):: rof
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
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur

  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

  if( nbseuil == 2 ) then
    if( Final_tddft .and. nbseuil /= ninitlv ) then
      ninit1 = ( ninitlv - 2 ) / 2
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
  if( Solsing ) gmi(1:nr-1,:) = 1 / gm(1:nr-1,:)

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

  do l = 0,lfin

    if( Hubb_a .and. l == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif
    Radial_comp = Ecomp .or. ( Hubb_m .and. Ylm_comp )

    if( Full_potential ) then
      nlm1 = nlm    ! = ( lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*l + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*l + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2, Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde, &
         l,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

! Radial integral for the golden rule
    call radial_matrix(Final_tddft,initlv,ip_max,ip0,iseuil,l,nlm1,nlm2,nbseuil, &
           ninitlv,nlma,nlma2,nr,NRIXS,nrm,nrmtsd,nspino,nspinp,psii,r,r_or_bess,Radial_comp,Rmtsd,rof,ui,ur,Vecond)

! Calculation of irregular solution
    if( Solsing ) call Cal_Solsing(Classic_irreg,Ecomp,Eimag,f2,Final_tddft,Full_potential,g0,gmi,gp,gso,Hubb_a,Hubb_d, &
         icheck,initlv,ip_max,ip0,iseuil,konde,l,lmax,m_hubb,nlm, &
         nlm1,nlm2,nbseuil,ninitlv,nlma,nr,NRIXS,nrm,nrmtsd,nspin,nspino,nspinp, &
         numat,psii,r,r_or_bess,Radial_comp,Rmtsd,Singul,Spinorbite,Tau,ui,ur,V,V_hubb,Vecond)

    deallocate( Tau, ui, ur )

  end do   ! end of loop over l

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
          write(3,170) '  l  m lp mp', mot
        else
          write(3,175) '  l  m', mot
        endif
        lm = 0
        do l = 0,lmax
          do m = -l,l
            lm = lm + 1
            if( .not. ( Hubb_a .or. Full_potential .or. Spinorbite ) .and. m /= 0 ) cycle
            lmp = 0
            do lp = 0,lmax
              if( .not. Full_potential .and. l /= lp ) cycle
              do mp = -lp,lp
                if( .not. ( Hubb_a .or. Full_potential .or. m == mp ) ) cycle
                lmp = lp**2 + lp + 1 + mp
                lmp = min(nlma2,lmp)
                if( Full_potential .or. Hubb_a ) then
                  if( Ecomp) then
                    write(3,180) l, m, lp, mp, ( rof(lm,lmp,isp,:,ip,initl), isp = 1,nspinp)
                  else
                    write(3,180) l, m, lp, mp, ( real(rof(lm,lmp,isp,:,ip,initl)), isp = 1,nspinp)
                  endif
                else
                  if( Ecomp) then
                    write(3,185) l, m, ( rof(lm,lmp,isp,:,ip,initl), isp = 1,nspinp )
                  else
                    write(3,185) l, m, ( real(rof(lm,lmp,isp,:,ip,initl)), isp = 1,nspinp )
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

          if( nspinp == 2 ) then
            write(3,210) '  l  m','   up  ','   dn  '
          else
            write(3,210) '  l  m','  Sing '
          endif

          lm = 0
          do l = 0,lmax
            do m = -l,l
              lm = lm + 1
              if( .not. ( Hubb_m .or. Full_potential .or. Spinorbite ) .and. m /= 0 ) cycle
              write(3,230) l, m, Singul(lm,:,ip1,ip2,initl)
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
  210 format(a6,9x,a7,3(19x,a7))
  230 format(2i3,1p,4e13.5)
end

!**********************************************************************

! Calculatim of the radial integral in the (Fermi) Golden rule.
! psii, u and ur are wave functions time r.

subroutine radial_matrix(Final_tddft,initlv,ip_max,ip0,iseuil,l,nlm1,nlm2,nbseuil, &
           ninitlv,nlma,nlma2,nr,NRIXS,nrm,nrmtsd,nspino,nspinp,psii,r,r_or_bess,Radial_comp,Rmtsd,rof,ui,ur,Vecond)

  use declarations
  implicit none

  integer initlv, ip, ip_max, ip0, iseuil, is, isol, isp, iss, l, lm, lm1, lm2, nlm1, nlm2, n1, n2, nbseuil, &
    ninitlv, nlma, nlma2, nr, nrm, nrmtsd, ns1, ns2, nspino, nspinp

  complex(kind=db), dimension(nlma,nlma2,nspinp,nspino,ip0:ip_max, ninitlv):: rof

  logical:: Final_tddft, NRIXS, Radial_comp

  real(kind=db):: f_integr3, fac, radlr, radli, Rmtsd
  real(kind=db), dimension(nbseuil):: Vecond
  real(kind=db), dimension(nr):: r, fct
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,ip0:ip_max):: r_or_bess
  real(kind=db), dimension(nrmtsd,nbseuil):: psiir
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur

  if( Final_tddft ) then
    ns1 = initlv
    ns2 = initlv
  else
    ns1 = 1
    ns2 = nbseuil
  endif

  do is = 1,nbseuil
    psiir(1:nrmtsd,is) =  psii(1:nrmtsd,is) * r(1:nrmtsd) 
  end do

  lm1 = l**2
  lm2 = (l + 1)**2

  do ip = ip0,ip_max
  
    do is = ns1,ns2

      if( Final_tddft ) then
        iss = nbseuil   ! c'est le 2eme seuil qui sert d'origine
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

              fct(1:nrmtsd) = psiir(1:nrmtsd,iseuil) * ur(1:nrmtsd,n1,n2,isp,isol) * r_or_bess(1:nrmtsd,ip)
              radlr = fac * f_integr3(r,fct,1,nr,Rmtsd)
              if( Radial_comp ) then
                fct(1:nrmtsd) = psiir(1:nrmtsd,iseuil) * ui(1:nrmtsd,n1,n2,isp,isol) * r_or_bess(1:nrmtsd,ip)
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

! Calcul de l'integrale radiale de la regle d'or de Fermi.
! ui et ur sont les fonctions d'onde fois r.

subroutine radial_matrix_optic(ip_max,ip0,ne,nlm1g,nlm_fp,nr,nrmtsd,nspinp,nspino,r,Radial_comp,Rmtsd,roff_rr,ui,ur,Vecond)

  use declarations
  implicit none

  integer:: ip, ip_max, ip0, ipp, ir, isol, isoli, isolf, isp, ispi, ispf, nlm1g, nlm_fp, n1i, n1f, n2i, n2f, ne, &
    nr, nrmtsd, nspinp, nspino

  logical:: Radial_comp

  real(kind=db):: f_integr3, fac, Rmtsd, Vecond
  real(kind=db), dimension(nr):: r, rp, fct
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

                    roff_rr(n1i,n1f,n2i,n2f,isp,isol,ip) = fac * f_integr3(r,fct,1,nr,Rmtsd)

                    if( Radial_comp ) then
                      do ir = 1,nrmtsd
                        fct(ir) = ur(ir,n1i,n2i,ispi,isoli,1) * ui(ir,n1f,n2f,ispf,isolf,2) * rp(ir)
                      end do

                      roff_ri(n1i,n1f,n2i,n2f,isp,isol,ip) = fac * f_integr3(r,fct,1,nr,Rmtsd)

                      do ir = 1,nrmtsd
                        fct(ir) = ui(ir,n1i,n2i,ispi,isoli,1) * ur(ir,n1f,n2f,ispf,isolf,2) * rp(ir)
                      end do

                      roff_ir(n1i,n1f,n2i,n2f,isp,isol,ip) = fac * f_integr3(r,fct,1,nr,Rmtsd)

                      do ir = 1,nrmtsd
                        fct(ir) = ui(ir,n1i,n2i,ispi,isoli,1) * ui(ir,n1f,n2f,ispf,isolf,2) * rp(ir)
                      end do

                      roff_ii(n1i,n1f,n2i,n2f,isp,isol,ip) = fac * f_integr3(r,fct,1,nr,Rmtsd)

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

! Calcul de la solution singuliere pour l'absorption
! Appele par Radial

Subroutine Cal_Solsing(Classic_irreg,Ecomp,Eimag,f2,Final_tddft,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d, &
         icheck,initlv,ip_max,ip0,iseuil,konde,ll,lmax,m_hubb,nlm, &
         nlm1,nlm2,nbseuil,ninitlv,nlma,nr,NRIXS,nrm,nrmtsd,nspin,nspino,nspinp, &
         numat,psii,r,r_or_bess,Radial_comp,Rmtsd,Singul,Spinorbite,Tau,ui,ur,V,V_hubb,Vecond)

  use declarations
  implicit none

  integer:: icheck, initlv, ip_max, ip0, ip1, ip2, is, ise, iseuil, isol, isp, iss,  l, ll, lmax, m_hubb, m, &
    ms, mv, n0, n, nbseuil, ninitlv, nlm, nlm1, nlm2, np, nlma, nr, nrm, nrmtsd, ns1, ns2, nspin, nspino, nspinp, numat, nv

  complex(kind=db):: integr_sing, Sing
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlma,nspinp,ip0:ip_max,ip0:ip_max, ninitlv):: Singul
  complex(kind=db), dimension(nrmtsd):: f_reg, f_irg
  complex(kind=db), dimension(nlm1,nspinp,nlm1,nspinp):: Tau
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  logical:: Classic_irreg, Ecomp, Final_tddft, Full_potential, Hubb_a, Hubb_d, NRIXS, Radial_comp, Spinorbite

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

!  real(kind=db), dimension(nrmtsd+1,nlm1,nlm2,nspinp,nspino):: usi, usr
  real(kind=db), dimension(:,:,:,:,:), allocatable:: usi, usr

  rr(1:nrmtsd) = r(1:nrmtsd)

  if( Final_tddft ) then
    ns1 = initlv
    ns2 = initlv
  else
    ns1 = 1
    ns2 = nbseuil
  endif

  allocate( usi(nrmtsd+1,nlm1,nlm2,nspinp,nspino) )
  allocate( usr(nrmtsd+1,nlm1,nlm2,nspinp,nspino) )

! gm et gp inverse dans le sousprogramme
  call Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde, &
       ll,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtsd,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Spinorbite,Tau,usi,usr,V,V_hubb)

  do is = ns1,ns2

    if( Final_tddft ) then
      iss = nbseuil   ! c'est le 2eme seuil qui sert d'origine
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

      do ip2 = ip0,ip_max  ! boucle sur dipole, quadrupole, Octupole

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

        n = 0
        do l = 0,lmax
          if( .not. Full_potential .and. l /= ll ) cycle
          n0 = l**2 + l + 1
          do m = -l,l
            if( nlm1 == 1 .and. m /= 0 ) cycle
            n = n + 1

            do isp = 1,nspinp

              do np = 1,nlm2

                do isol = 1,nspino
                  if( Spinorbite .and. nlm2 == 1 ) then
                    ms = m + isol - isp
                    if( ms > l .or. ms < -l ) cycle
                  endif

                  if( Radial_comp ) then
                    f_reg(1:nrmtsd) = r(1:nrmtsd) * cmplx( ur(1:nrmtsd,n,np,isp,isol), ui(1:nrmtsd,n,np,isp,isol), db )
                  else
                    f_reg(1:nrmtsd) = r(1:nrmtsd) * cmplx( ur(1:nrmtsd,n,np,isp,isol), 0._db, db )
                  endif

                  f_irg(1:nrmtsd) = cmplx( usr(1:nrmtsd,n,np,isp,isol), usi(1:nrmtsd,n,np,isp,isol), db)

                  Sing = fac1 * fac2 * integr_sing(nrmtsd,phi1,phi2,f_reg,f_irg,Rmtsd,rr,icheck)

                  if( nlm1 == 1 ) then
                    do mv = -l,l
                      nv = n0 + mv
                      Singul(nv,isp,ip1,ip2,is) = Singul(nv,isp,ip1,ip2,is) + Sing
                    end do
                  else
                    nv = n0 + m
                    Singul(nv,isp,ip1,ip2,is) = Singul(nv,isp,ip1,ip2,is) + Sing
                  endif

                end do
              end do

            end do    ! fin boucle isp
          end do   ! fin boucle m
        end do  ! fin boucle l

      end do ! fin boucle dipole, quadrupole
    end do ! fin boucle dipole, quadrupole

  end do ! fin boucle seuil

  deallocate( usi, usr )
  
  return
end

!***********************************************************************

! Resolution de l'equation de Schrodinger radiale pour la solution irreguliere

! Recopie de Sch_radial: gm et gp sont inverses dans l'appel

subroutine Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gp,gm,gso,Hubb_a,Hubb_d,icheck,konde, &
       ll,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtsd,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Spinorbite,Tau,usi,usr,V,V_hubb)

  use declarations
  implicit none

  integer:: icheck, im, ip, ir, isol, isp, ispp, isq, isr, l, l_hubbard, l2, li, ll, lf, lmax, lp, m, m_hubb, mp, ms, n, &
    nlm, nlm1, nlm2, np, nr, nrmtsd, ns, nspin, nspino, nspinp, numat

  logical:: Classic_irreg, Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m, Hubb_nd, Radial_comp, Spinorbite

  complex(kind=db):: fnormc, V_h, z
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(0:lmax):: bess, neum
  complex(kind=db), dimension(nlm1,nspinp,nlm1,nspinp):: Tau
  complex(kind=db), dimension(nrmtsd:nrmtsd+1,0:lmax,nspin):: Bessel, Hankel
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

  real(kind=db), dimension(nrmtsd+1,nlm1,nlm2,nspinp,nspino):: usi, usr
  
  allocate( us(nrmtsd+1,nlm1,nlm2,nspinp,nspino)  )

  us(:,:,:,:,:) = ( 0._db, 0._db )

! Hankel est en fait -i.hankel+
    do ir = nrmtsd,nrmtsd+1
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
    li = ll
    lf = ll
  endif

  n = 0
  do l = li,lf

    do m = -l,l
      if( nlm1 == 1 .and. m /= 0 ) cycle
      n = n + 1

      do isp = 1,nspinp  ! Spin reel
        isr = min( isp, nspin )
        
        np = 0
        do lp = li,lf
          if( .not. Full_potential .and. l /= lp ) cycle
          do mp = -lp,lp
            if( nlm2 == 1 .and. mp /= m ) cycle
            np = np  + 1

            do isq = 1,nspinp  ! spin d'attaque ( indice solution )
              isol = min( isq, nspino )
              if( nspino == 1 .and. isp /= isq ) cycle
              if( nlm2 == 1 ) then
                ms = m + isq - isp
                if( ms > l .or. ms < -l ) cycle
                ns = n + isq - isp
              else
                ns = np
              endif

              do ir = nrmtsd,nrmtsd+1
                if( Classic_irreg ) then
                  us(ir,n,np,isp,isol) = Tau(ns,isq,n,isp) * r(ir) * Hankel(ir,l,isr)
                else
                  us(ir,n,np,isp,isol) = r(ir) * Bessel(ir,l,isr)
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

      do l = li,lf

        l2 = l * ( l + 1 )

        if( Hubb_a .and. l == l_hubbard(numat) ) then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif
        Hubb_nd = Hubb_m .and. .not. Hubb_d

        do m = -l,l

          if( Full_potential ) then
            n = l**2 + l + 1 + m
          elseif( Hubb_m .or. Spinorbite ) then
            n = l + 1 + m
          else
            n = 1
            if( m /= 0 ) cycle
          endif

          td = g0(ir,isr) + l2 * f2(ir)

          if( Spinorbite ) then
            if( isp == 1 ) then
              ms = m
              mp = m + 1   ! m de l'autre spin
            else
              ms = - m
              mp = m - 1
            endif
            fac = sqrt( ( l - ms ) * ( l + ms + 1._db ) )
            td = td + ms * gso(ir,isp)
          else
            mp = l + 1 ! pour que le if en dessous soit faux
          endif
          if( Full_potential ) then
            np = l**2 + l + 1 + mp
          elseif( Hubb_nd .or. Spinorbite ) then
            np = l + 1 + mp
          else
            np = 1
          endif

! le deuxieme indice "m" est la composante non nule a l'origine
          us(ip,n,:,isp,:) = td * us(ir,n,:,isp,:) + gm(ir,isr) * us(im,n,:,isp,:)
          if( Ecomp ) us(ip,n,:,isp,:) = us(ip,n,:,isp,:) - img * Eimag * us(ir,n,:,isp,:)

          if( abs(mp) <= l ) us(ip,n,:,isp,:) = us(ip,n,:,isp,:) + fac * gso(ir,isr) * us(ir,np,:,isq,:)

          if( Hubb_m ) then
            do mp = -l,l
              if( Hubb_d .and. m /= mp ) cycle
              do ispp = 1,nspinp
                V_h = V_hubb(m,mp,isp,ispp) 
                if( Full_potential ) then
                  np = l**2 + l + 1 + mp
                else
                  np = l + 1 + mp
                endif
                us(ip,n,:,isp,:) = us(ip,n,:,isp,:) + V_h * us(ir,np,:,ispp,:)
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
    n = nrmtsd + 1
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

  integer:: i, icheck, ie, ie1, ie2, ief, ip, ip_max, ip0, iso1, iso2, isp, isp1, isp2, j, l, lf, l_hubbard, lfin, lm, lmf, &
    lm1, lm2, lmax, lmax_pot, lmp, lmpf, lp, lpf, m, m_hubb, mf, n_ph, nlm, nlm_pot, nlm1, nlm1g, nlm2, nlm_fp, mp, mpf, &
    n_Ec, n_V, nr, nr_zet, nrmtsd, nrmtg, nspin, nspino, nspinp, numat

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

  do l = 0,lfin

    if( Hubb_a .and. l == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = ( lmax + 1 )**2
      nlm2 = nlm1
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*l + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*l + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )
    allocate( us(nrmtsd+1,nlm1,nlm2,nspinp,nspino) )
    allocate( usi(nrmtsd+1,nlm1,nlm2,nspinp,nspino) )
    allocate( usr(nrmtsd+1,nlm1,nlm2,nspinp,nspino) )

    if( l == 0 ) allocate( us_t(nrmtsd+1,nlm1g,nlm_fp,nspinp,nspino,n_Ec) )

    do ie = 1,n_Ec

      if( icheck > 1 ) write(3,110) ie, Enervide(ie)*rydb, Ecinetic_t(:,ie)*rydb 
      if( icheck > 1 .and. n_V <= 2 ) write(3,120) V0bd_t(1,:)*rydb 
     
      Ecinetic(:) = Ecinetic_t(:,ie)
      Eimag = Eimag_t(ie)
      konde(:) = sqrt( cmplx(Ecinetic(:), Eimag, db) )

      call coef_sch_rad(Enervide(ie),f2,g0,gm,gp,gso,nlm,nr,nspin,nspino,numat,r,Relativiste,Spinorbite,V)

      gpi(:,:) = 1 / gp(:,:)

      call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,l,lmax,m_hubb, &
        nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
        ur,V,V_hubb)

! Calcul de la solution singuliere
! gm et gp inverse dans le sousprogramme
      if( Solsing ) then
        gmi(:,:) =  1 / gm(:,:)
        call Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gmi,gp,gso,Hubb_a,Hubb_d,icheck,konde, &
           l,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtsd,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Spinorbite,Tau,usi,usr,V,V_hubb)
      endif

      us(:,:,:,:,:) = cmplx( usr(:,:,:,:,:), usi(:,:,:,:,:), db )
      
      if( Radial_comp ) then
        if( Full_potential ) then
          ur_t(:,:,:,:,:,ie) = ur(:,:,:,:,:)
          ui_t(:,:,:,:,:,ie) = ui(:,:,:,:,:)
          us_t(:,:,:,:,:,ie) = us(:,:,:,:,:)
        elseif( Hubb_m .and. .not. Hubb_d ) then
          ur_t(:,l**2+1:l**2+nlm1,l**2+1:l**2+nlm1,:,:,ie) = ur(:,1:nlm1,1:nlm1,:,:)
          ui_t(:,l**2+1:l**2+nlm1,l**2+1:l**2+nlm1,:,:,ie) = ui(:,1:nlm1,1:nlm1,:,:)
          us_t(:,l**2+1:l**2+nlm1,l**2+1:l**2+nlm1,:,:,ie) = us(:,1:nlm1,1:nlm1,:,:)
        elseif( No_diag ) then
          do lm1 = l**2+1,(l+1)**2
            lm = min( lm1 - l**2, nlm1 )
            ur_t(:,lm1,lm1,:,:,ie) = ur(:,lm,1,:,:)
            ui_t(:,lm1,lm1,:,:,ie) = ui(:,lm,1,:,:)
            us_t(:,lm1,lm1,:,:,ie) = us(:,lm,1,:,:)
          end do
        elseif( Spinorbite .or. Hubb_m ) then
           do lm1 = l**2+1,(l+1)**2
             lm = lm1 - l**2
             ur_t(:,lm1,1,:,:,ie) = ur(:,lm,1,:,:)
             ui_t(:,lm1,1,:,:,ie) = ui(:,lm,1,:,:)
             us_t(:,lm1,1,:,:,ie) = us(:,lm,1,:,:)
           end do
        elseif( m_depend ) then
          do lm1 = l**2+1,(l+1)**2
            ur_t(:,lm1,1,:,:,ie) = ur(:,1,1,:,:)
            ui_t(:,lm1,1,:,:,ie) = ui(:,1,1,:,:)
            us_t(:,lm1,1,:,:,ie) = us(:,1,1,:,:)
          end do
        else
          ur_t(:,l,1,:,:,ie) = ur(:,1,1,:,:)
          ui_t(:,l,1,:,:,ie) = ui(:,1,1,:,:)
          us_t(:,l,1,:,:,ie) = us(:,1,1,:,:)
        endif
      else
        if( Full_potential ) then
          ur_t(:,:,:,:,:,ie) = ur(:,:,:,:,:)
        elseif( Hubb_m .and. .not. Hubb_d ) then
          ur_t(:,l**2+1:l**2+nlm1,l**2+1:l**2+nlm1,:,:,ie) = ur(:,1:nlm1,1:nlm1,:,:)
        elseif( No_diag ) then
          do lm1 = l**2+1,(l+1)**2
            lm = min( lm1 - l**2, nlm1 )
            ur_t(:,lm1,lm1,:,:,ie) = ur(:,lm,1,:,:)
          end do
        elseif( Spinorbite .or. Hubb_m ) then
           do lm1 = l**2+1,(l+1)**2
             lm = lm1 - l**2
             ur_t(:,lm1,1,:,:,ie) = ur(:,lm,1,:,:)
           end do
        elseif( m_depend ) then
          do lm1 = l**2+1,(l+1)**2
            ur_t(:,lm1,1,:,:,ie) = ur(:,1,1,:,:)
          end do
        else
          ur_t(:,l+1,1,:,:,ie) = ur(:,1,1,:,:)
        endif
      endif

    end do

    deallocate( Tau, ui, ur, us, usi, usr )

  end do   ! fin de la boucle sur l

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
        write(3,170) ' l1 m1 lp mp l2 m2 lp mp', mot
      else
        write(3,175) ' l1 m1 l2 m2', mot
      endif

      do l = 0,lmax
        do m = -l,l
          if( m_depend ) then
            lm = l**2 + l + 1 + m
          else
            if( m /= 0 ) cycle
            lm = l + 1
          endif
          do lp = 0,lmax
            if( .not. Full_potential .and. l /= lp ) cycle
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
                        write(3,180) l, m, lp, mp, lf, mf, lpf, mpf, (roff_rr(lm,lmf,lmp,lmpf,isp,:,ip), isp = 1,nspinp**2)
                      else
                        write(3,185) l, m, lf, mf, (roff_rr(lm,lmf,lmp,lmpf,isp,:,ip), isp = 1,nspinp**2)
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

! Calcul des fonctions de base sur les points du maillage

subroutine cal_phiato(Full_potential,ia,iang,ibord,icheck,iopsymr,ll,lmax,nlm1,nlm2,natome, &
             nbord,nbtm,nbtm_fdm,nlmagm,nlmmax,nphiato1,nphiato7,npsom,nr,nspinp,nspinpo,phiato,posi,r, &
             Radial_comp,Spinorbite,ui,ur,xyz,Ylm_comp,Ylmato)

  use declarations
  implicit none

  integer:: i, ia, iang, ib, icheck, isol, isp, isym, jsol, jsp, l, &
    ll, lm, lm0, lmax, lmv, lp, m, nlm1, nlm2, mp, mpp, mv, n, natome, nbord, nbtm, nbtm_fdm, nlmagm, nlmmax, np, &
    nphiato1, nphiato7, npsom, nr, nspinp, nspinpo
  integer, dimension(nopsm):: iopsymr
  integer, dimension(nbtm,natome):: ibord

  complex(kind=db):: phic, phid, yc

  logical:: Ylm_comp, Full_potential, Radial_comp, Spinorbite

  real(kind=db):: f_interp2, phi, phii, r0, rm, rp, rrel, u0, uc0, ucm, ucp, um, up
  real(kind=db), dimension(3):: ps, x, w
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspinpo):: ui, ur
  real(kind=db), dimension(4,npsom):: xyz
  real(kind=db), dimension(nbtm_fdm,nlmmax,natome):: Ylmato
  real(kind=db), dimension(nphiato1,nlmagm,nspinp,nspinpo,natome,nphiato7):: phiato

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
    do l = 0,lmax
      if( .not. Full_potential .and. l /= ll ) cycle
      lm0 = l**2 + l + 1
      do m = -l,l
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
            if( .not. Full_potential .and. lp /= l ) cycle
            do mp = -lp,lp
              if( nlm2 == 1 .and. mp /= m ) cycle
              np = np + 1
              do isol = 1,nspinpo
                if( iang == 1 ) then
                  jsol = isol
                elseif( iang == - 1 ) then
                  jsol = nspinpo + 1 - isol
                endif

                if( nlm2 > 1 ) then
                  lmv = lp**2 + lp + 1 + mp
                elseif( Spinorbite ) then
                  mv = m + isol - isp  ! m  faux
                  if( mv > l .or. mv < -l ) cycle
                  lmv = lm0 + m
                else
                  lmv = lm0 + m
                endif
                mv = m

                if( iang == 0 ) then
                  if( nspinpo == 2 .and. isp /= isol ) cycle
                  um = 0.5_db * ( ur(i-2,n,np,1,1) + ur(i-2,n,np,min(2,nspinp),min(2,nspinpo)) )
                  u0 = 0.5_db * ( ur(i-1,n,np,1,1) + ur(i-1,n,np,min(2,nspinp),min(2,nspinpo)) )
                  up = 0.5_db * ( ur(i,n,np,1,1) + ur(i,n,np,min(2,nspinp),min(2,nspinpo)) )
                else
                  um = ur(i-2,n,np,jsp,jsol)
                  u0 = ur(i-1,n,np,jsp,jsol)
                  up = ur(i,n,np,jsp,jsol)
                endif

                if( Radial_comp ) then
                  if( iang == 0 ) then
                    if( nspinpo == 2 .and. isp /= isol ) cycle
                    ucm = 0.5_db * ( ui(i-2,n,np,1,1) + ui(i-2,n,np,min(2,nspinp),min(2,nspinpo)) )
                    uc0 = 0.5_db * ( ui(i-1,n,np,1,1) + ui(i-1,n,np,min(2,nspinp),min(2,nspinpo)) )
                    ucp = 0.5_db * ( ui(i,n,np,1,1) + ui(i,n,np,min(2,nspinp),min(2,nspinpo)) )
                  else
                    ucm = ui(i-2,n,np,jsp,jsol)
                    uc0 = ui(i-1,n,np,jsp,jsol)
                    ucp = ui(i,n,np,jsp,jsol)
                  endif
                endif

                phi = f_interp2(rrel,rm,r0,rp,um,u0,up)
                if( Radial_comp ) then
                  phii = f_interp2(rrel,rm,r0,rp,ucm,uc0,ucp)
                else
                  phii = 0._db
                endif
                phic = cmplx( phi, phii, db )

                do mpp = -l,l
                  if( nlm1 /= 1 .and. mpp /= m ) cycle
                  lm = lm0 + mpp
                  if( nlm1 == 1 ) then
                    lmv = lm
                    mv = mpp
                  endif
                  if( mpp == 0 .or. .not. Ylm_comp ) then
                    phiato(ib,lmv,isp,isol,ia,1) = phiato(ib,lmv,isp,isol,ia,1) + phi * Ylmato(ib,lm,ia)
                    if( Radial_comp ) phiato(ib,lmv,isp,isol,ia,2) = phiato(ib,lmv,isp,isol,ia,2) + phii * Ylmato(ib,lm,ia)
                  else
                    phid = phic * yc(mv,Ylmato(ib,lm,ia),Ylmato(ib,lm-2*mv,ia))
                    phiato(ib,lmv,isp,isol,ia,1) = phiato(ib,lmv,isp,isol,ia,1) + Real( phid, db)
                    phiato(ib,lmv,isp,isol,ia,2) = phiato(ib,lmv,isp,isol,ia,2) +  aimag( phid )
                  endif

                end do
              end do
            end do
          end do

        end do
      end do
    end do
  end do

  if( icheck > 2 .and. ( ll == lmax .or. Full_potential ) ) then
    write(3,110) ia
    do l = 0,lmax
      lm0 = l**2 + l + 1
      do mp = -l,l
        if( nphiato7 == 1 ) then
          if( Spinorbite ) then
            write(3,120) l, mp
          else
            write(3,130) l, mp
          endif
        else
          if( Spinorbite ) then
            write(3,140) l, mp
          else
            write(3,150) l, mp
          endif
        endif
        lm = lm0 + mp
        do ib = 1,nbord
          write(3,160) ibord(ib,ia), ( ( phiato(ib,lm,isp,isol,ia,:), isp = 1,nspinp ), isol = 1,nspinpo )
        end do
      end do
    end do
  endif

  return
  110 format(/' ia =',i3)
  120 format(/' ibord phiato(u,1)  phiato(d,1)  phiato(u,2)  ', 'phiato(d,2)   l =',i3,', m =',i3)
  130 format(/' ibord  phiato(u)    phiato(d)   l =',i3,', m =',i3)
  140 format(/' ibord',6x,' phiato(u,1)',13x,'  phiato(d,1)',13x, '  phiato(u,2)  ',13x,'phiato(d,2)   l =',i3,', m =',i3)
  150 format(/' ibord',6x,'  phiato(u)',13x,'    phiato(d)   l =',i3, ', m =',i3)
  160 format(i5,1p,8e13.4)
end

!***********************************************************************

! Calculation of density of state

subroutine Cal_dens(Cal_xanes,Classic_irreg,Density_comp,drho_self,Ecinetic,Eimag,Energ,Enervide,Full_atom, &
            Full_potential,Hubb,Hubb_diag,iaabsi,iaprabs,iaprotoi,icheck,itypei,itypepr,lla2_state,lmax_pot,lmaxat,m_hubb, &
            mpinodes,mpirank,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nlmagm,nlm_pot, &
            nrato,nrm,nrm_self,nspin,nspino,nspinp,ntype,numat,rato,Relativiste,Renorm,Rmtg,Rmtsd,Solsing,Solsing_only, &
            Spinorbite,State_all,Statedens,Statedens_i,Taull,V_hubb,V_hubb_abs,V_intmax,V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer ia, iaabsi, iapr, iaprabs, icheck, ipr, ir, isp, isp1, isp2, it, lla2_state, lm, lmax, lmax_pot, &
    m_hubb, mpinodes, mpirank, n_atom_0, n_atom_0_self, n_atom_ind, n_atom_ind_self, &
    n_atom_proto, natome, nlm_pot, nlma, nlma2, nlmagm, nr, nrs, nrm, nrm_self, nspin, nspino, nspinp, ntype, Z

  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,natome):: Taull
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb_abs
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb
  complex(kind=db), dimension(:,:,:,:), allocatable:: V_hubb_t
  complex(kind=db), dimension(:,:,:,:), allocatable:: Taulla

  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, lmaxat
  integer, dimension(natome):: iaprotoi, itypei

  logical:: Absorbeur, Cal_xanes, Classic_irreg, Density_comp, Full_atom, Full_potential, Hubb_a, Hubb_d, Relativiste, &
    Renorm, Self, Solsing, Solsing_only, Spinorbite, State_all, Ylm_comp, Ylm_mod
  logical, dimension(0:ntype):: Hubb
  logical, dimension(n_atom_0:n_atom_ind):: iapr_done
  logical, dimension(n_atom_0_self:n_atom_ind_self):: Hubb_diag

  real(kind=db):: Eimag, Energ, Enervide, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(0:n_atom_proto):: Rmtg, Rmtsd
  real(kind=db), dimension(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1):: Statedens, Statedens_i
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vrato
  real(kind=db), dimension(nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self,0:mpinodes-1):: drho_self
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

    if( Self ) then
      nrs = nr
    else
      nrs = 0
    endif
    allocate( r(nr) )
    allocate( drho(nrs,nlm_pot,nspin) )
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
    if( Hubb(it) .and. iapr <= n_atom_ind_self ) then
      Hubb_a = .true.
      if( Absorbeur ) then
        V_Hubb_t(:,:,:,:) = V_Hubb_abs(:,:,:,:)
        Hubb_d = Hubb_diag(iaprabs)
      else
        V_Hubb_t(:,:,:,:) = V_Hubb(:,:,:,:,iapr)
        Hubb_d = Hubb_diag(iapr)
      endif
    else
      Hubb_a = .false.
      Hubb_d = .true.
    end if

! Ylm_mod = .true. correspond au calcul en harmo complexe qu'on projete sur des harmo reelles
!    if( ( Cal_xanes .and. Density_comp ) .or. ( Hubb_a .and. Spinorbite ) ) then
    if( ( Cal_xanes .and. Density_comp ) .or. .not. Cal_xanes ) then
      Ylm_mod = .false.
    else
      Ylm_mod = Ylm_comp
    endif

    do isp1 = 1,nspinp
      do lm = 1,nlma
        do isp2 = 1,nspinp
          Taulla(1:nlma,isp1,lm,isp2) = Taull(1:nlma,isp1,lm,isp2,ia)
        end do
      end do
    end do

    call radial_sd(Classic_irreg,drho,Ecinetic,Eimag,Energ,Enervide,Full_potential,Hubb_a,Hubb_d,ia,icheck,lmax, &
         lmax_pot,m_hubb,nlm_pot,nlma,nlma2,nr,nrs,nspin,nspino,nspinp,Z,r,Relativiste,Renorm,Rmtg(ipr), &
         Rmtsd(ipr),Self,Solsing,Solsing_only,Spinorbite,State,State_i,Taulla,V_hubb_t,V_intmax,V0bd,Vrato_t,Ylm_comp,Ylm_mod)

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
      do isp = 1,nspin
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

! Calcul des integrales radiales et de la solution singuliere pour la densite d'etat.
! Ylm_mod = .true. correspond au calcul en harmo complexe qu'on projete sur des harmo reelles
! Appele par cal_dens

subroutine radial_sd(Classic_irreg,drho,Ecinetic,Eimag,Energ,Enervide,Full_potential,Hubb_a,Hubb_d,ia,icheck,lmax, &
         lmax_pot,m_hubb,nlm_pot,nlma,nlma2,nr,nrs,nspin,nspino,nspinp,numat,r,Relativiste,Renorm,Rmtg, &
         Rmtsd,Self,Solsing,Solsing_only,Spinorbite,State,State_i,Taull,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp,Ylm_mod)

  use declarations
  implicit none

  integer:: ia, icheck, l, l_hubbard, lfin, lmax, lmax_pot, m_hubb, nlm, nlm_pot, nlm1, nlm2, nlma, nlma2, nr, &
    nrs, nrmtg, nrmtsd, nspin, nspino, nspinp, numat

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlma,nspinp,nlma,nspinp):: Taull
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(:,:,:,:), allocatable:: Tau

  logical:: Classic_irreg, Diagonale, Ecomp, Full_potential, Hubb_a, Hubb_d, Hubb_m, Radial_comp, Relativiste, Renorm, &
    Self, Solsing, Solsing_only, Spinorbite, Ylm_comp, Ylm_mod

  real(kind=db):: Eimag, Energ, Enervide, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nr,nspin):: g0, gm, gmi, gp, gpi
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato
  real(kind=db), dimension(nrs,nlm_pot,nspin):: drho
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
  Radial_comp = Ecomp .or. ( Hubb_a .and. Ylm_comp )

  if( Full_potential ) then
    lfin = 0
  else
    lfin = lmax
  endif

!      if( Ylm_comp .and. Hubb_a )
!     &    call Trans_Tau(icheck,0,lmax,nspin,Spinorbite,Taull,.true.)

  do l = 0,lfin

    if( Hubb_a .and. l == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = nlm ! = ( lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*l + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*l + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif
    Diagonale = .not. ( Full_potential .or. Hubb_m )

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,l,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

! Integrale radiale pour la densite d'etat
    if( .not. Solsing_only ) call Radial_matrix_sd(Diagonale,drho,Full_potential,icheck, &
           l,lmax,nlm_pot,nlm1,nlm2,nlma,nlma2,nr,nrs,nrmtsd,nspin,nspino,nspinp,r,Radial_comp,Rmtsd,State,State_i,Self, &
           Taull,ui,ur,Ylm_mod)

! Calcul de la solution singuliere
    if( Solsing ) call Cal_Solsing_sd(Classic_irreg,drho,Ecomp,Eimag,f2,Full_potential,g0,gmi,gp,gso,Hubb_a,Hubb_d,icheck,konde, &
          l,lmax,m_hubb,nlm,nlm_pot,nlm1,nlm2,nlma,nr,nrmtsd,nrs,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Self, &
          Spinorbite,State,Tau,ui,ur,V,V_hubb,Ylm_mod)

    deallocate( Tau )
    deallocate( ui, ur )

  end do   ! fin de la boucle sur l

  deallocate( V )

  return
  110 format(/' ---- Radial_sd ----',100('-'))
  120 format(/'     Atom =',i3,/'        Z =',i3,/ '     Rmtg =',f10.5,' A,  Rmtsd =',f10.5,' A',/ &
  '    Energ =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format('    V0bd  =',2f10.3)
  150 format('    konde =',2f12.5)
end

!**********************************************************************

! Calcul de l'integrale radiale pour la densite d'etat
! Ylm_mod = .true. correspond au calcul en harmo complexe qu'on projete sur des harmo reelles.

subroutine Radial_matrix_sd(Diagonale,drho,Full_potential,icheck,l,lmax,nlm_pot,nlm1,nlm2,nlma,nlma2,nr,nrs,nrmtsd, &
           nspin,nspino,nspinp,r,Radial_comp,Rmtsd,State,State_i,Self,Taull,ui,ur,Ylm_mod)

  use declarations
  implicit none

  integer:: icheck, ir, is1, is2, isg1, isg2, iso, iso1, iso2, isp, isp1, isp2, isr, l, l1, l2, l3, lm, lm0, lm1, lm2, lma1, &
    lma2, lmax, lp1, lp2, m1, m2, m3, mp1, mp2, mr1, mr2, mv, n1, n2, nlm_pot, nlm1, nlm2, nlma, nlma2, np1, np2, nr, &
    nrmtsd, nrs, nspin, nspino, nspinp

! The volatile statetement prevent for some compilation optimization.
! On linux, at times without it, one gets segmentation fault erros message during the running with magnetic system (YJ 2018/01/25)
  integer, volatile:: nr1, nr2
!  integer:: nr1, nr2
  
  complex(kind=db):: c_harm, c_harm1, c_harm2, rof_sd, Sta_e, Sta_s
  complex(kind=db), dimension(nspinp):: rof_sd_l
  complex(kind=db), dimension(nrmtsd):: fc
  complex(kind=db), dimension(nrmtsd,nspinp):: fc_l
  complex(kind=db), dimension(nlma,nspinp,nlma,nspinp):: Taull

  logical:: Diagonale, Full_potential, Radial_comp, Self, Ylm_mod

  real(kind=db):: f_integr3, g, Gaunt_r, gauntcp, rac2, radlr, radli, Rmtsd, Sta_i, Sta_r
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nrmtsd):: r2, rr, fcr, fci, ri2
  real(kind=db), dimension(nrs,nlm_pot,nspin):: drho
  real(kind=db), dimension(nlma,nspinp,nlma,nspinp):: State, State_i
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur

  rac2 = 1 / sqrt(2._db )

  if( Full_potential ) then
    lm0 = 0
  else
    lm0 = l**2
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
    elseif( nspino ==2 ) then
      write(3,120)
    else
      write(3,130)
    endif
  endif

  do isp1 = 1,nspinp
    isr = min(isp1,nspin)
    nr1 = 0
    do l1 = 0,lmax
      if( .not. Full_potential .and. l /= l1 ) cycle
      do mr1 = -l1,l1
        nr1 = nr1 + 1
        lm1 = lm0 + nr1

        do isp2 = 1,nspinp
          if( nspino == 1 .and. isp1 /= isp2 ) cycle
          nr2 = 0
          do l2 = 0,lmax
            if( .not. Full_potential .and. l /= l2 ) cycle
            do mr2 = -l2,l2
              nr2 = nr2 + 1

              if( Diagonale .and. mr1 /= mr2 ) cycle

              lm2 = lm0 + nr2

              do isg1 = -1,1,2
                if( Ylm_mod ) then
                  if( mr1 < 0 ) then
                    if( isg1 == 1 ) then
                      c_harm1 = cmplx( 0._db, rac2, db )
                    else
                      c_harm1 = cmplx( 0._db, - (-1)**mr1 * rac2, db )
                    endif
                  elseif( mr1 == 0 ) then
                    if( isg1 == -1 ) cycle
                    c_harm1 = ( 1._db, 0._db )
                  else
                    if( isg1 == 1 ) then
                      c_harm1 = cmplx( (-1)**mr1 * rac2, 0._db, db )
                    else
                      c_harm1 = cmplx( rac2, 0._db, db )
                    endif
                  endif
                  m1 = isg1*mr1
                else
                  if( isg1 == -1 ) cycle
                  m1 = mr1
                  c_harm1 = ( 1._db, 0._db )
                endif
                if( isg1 == 1 ) then
                  n1 = nr1
                else
                  n1 = nr1 - 2 * mr1
                endif

                do isg2 = -1,1,2
                  if( Ylm_mod ) then
                    if( mr2 < 0 ) then
                      if( isg2 == 1 ) then
                        c_harm2 = cmplx( 0._db, rac2,db)
                      else
                        c_harm2 = cmplx( 0._db, - (-1)**mr2 * rac2, db )
                      endif
                    elseif( mr2 == 0 ) then
                      if( isg2 == -1 ) cycle
                      c_harm2 = ( 1._db, 0._db )
                    else
                      if( isg2 == 1 ) then
                        c_harm2 = cmplx( (-1)**mr2 * rac2, 0._db, db )
                      else
                        c_harm2 = cmplx( rac2, 0._db, db )
                      endif
                    endif
                    m2 = isg2*mr2
                  else
                    if( isg2 == -1 ) cycle
                    m2 = mr2
                    c_harm2 = ( 1._db, 0._db )
                  endif
                  if( isg2 == 1 ) then
                    n2 = nr2
                  else
                    n2 = nr2 - 2 * mr2
                  endif

                  c_harm = c_harm1 * conjg( c_harm2 )

                  np1 = 0
                  do lp1 = 0,lmax
                    if( .not. Full_potential .and. lp1 /= l1 ) cycle

                    do mp1 = -lp1, lp1
                      if( nlm2 == 1 .and. mp1 /= m1 ) cycle
                      if( nlm2 == 1 ) then
                        np1 = 1
                      else
                        np1 = np1 + 1
                      endif

                      do iso1 = 1,nspino
!       if( iso1 /= isp1 ) cycle

                        if( nspino == 2 .and. nlm2 == 1 ) then
                          mv = mp1 + iso1 - isp1
                          if( mv > lp1 .or. mv < - lp1 ) cycle
                        else
                          mv = mp1
                        endif

                        if( nlma2 /= 1 ) then
                          is1 = iso1
                          lma1 = lp1**2 + lp1 + 1 + mp1
                        elseif( nspino == 2 ) then
                          is1 = iso1
                          lma1 = l1**2 + l1 + 1 + mv
                        else
                          is1 = isp1
                          lma1 = l1**2 + l1 + 1 + mv
                        endif

                        np2 = 0
                        do lp2 = 0,lmax
                          if( .not. Full_potential .and. lp2 /= l2 ) cycle

                          do mp2 = -lp2,lp2
                            if( nlm2 == 1 .and. mp2 /= m2 ) cycle
                            if( nlm2 == 1 ) then
                              np2 = 1
                            else
                              np2 = np2 + 1
                            endif

                            do iso2 = 1,nspino
!       if( iso2 /= isp2 ) cycle
                              if( nspino == 2 .and. nlm2 == 1 ) then
                                mv = mp2 + iso2 - isp2
                                if( mv > lp2 .or. mv < - lp2 ) cycle
                              else
                                mv = mp2
                              endif

                              if( nlma2 /= 1 ) then
                                is2 = iso2
                                lma2 = lp2**2 + lp2 + 1 + mp2
                              elseif( nspino == 2 ) then
                                is2 = iso2
                                lma2 = l2**2 + l2 + 1 + mv
                              else
                                is2 = isp2
                                lma2 = l2**2 + l2 + 1 + mv
                              endif

                              if( nlm1 == 1 ) then

                                rof_sd = rof_sd_l(isp1)
                                if( Self .and. isp1 == isp2 .and. lm1 == lm2 ) fc(:) = fc_l(:,isp1)

                              else

                                if( Radial_comp ) then
                                  fcr(1:nrmtsd) = r2(1:nrmtsd) * ( ur(1:nrmtsd,n1,np1,isp1,iso1) * ur(1:nrmtsd,n2,np2,isp2,iso2) &
                                                                 - ui(1:nrmtsd,n1,np1,isp1,iso1) * ui(1:nrmtsd,n2,np2,isp2,iso2) )
                                  radlr = f_integr3(rr,fcr,1,nrmtsd,Rmtsd)
                                  fci(1:nrmtsd) = r2(1:nrmtsd) * ( ur(1:nrmtsd,n1,np1,isp1,iso1) * ui(1:nrmtsd,n2,np2,isp2,iso2) &
                                                                 + ui(1:nrmtsd,n1,np1,isp1,iso1) * ur(1:nrmtsd,n2,np2,isp2,iso2) )
                                  radli = f_integr3(rr,fci,1,nrmtsd,Rmtsd)
                                else
                                  fcr(1:nrmtsd) = r2(1:nrmtsd) * ur(1:nrmtsd,n1,np1,isp1,iso1) &
                                                               * ur(1:nrmtsd,n2,np2,isp2,iso2)
                                  radlr = f_integr3(rr,fcr,1,nrmtsd, Rmtsd)
                                  fci(:) = 0._db
                                  radli = 0._db
                                endif

                                rof_sd = cmplx( radlr, radli, db)

!                                    if( Self .and. isp1 == isp2
!     &                                 .and. lm1 == lm2 ) fc(1:nrmtsd)
                                if( Self .and. isp1 == isp2 ) fc(1:nrmtsd) = cmplx( fcr(1:nrmtsd), &
                                   fci(1:nrmtsd),db) * ri2(1:nrmtsd)

                              endif

                              Sta_e = c_harm * rof_sd * Taull(lma1,is1,lma2,is2)
                              Sta_s = conjg( c_harm ) * rof_sd * Taull(lma2,is2,lma1,is1)

                              Sta_r = 0.5_db * real( Sta_e - Sta_s, db)
                              Sta_i = 0.5_db * aimag( Sta_e + Sta_s )

                              State(lm1,isp1,lm2,isp2) = State(lm1,isp1,lm2,isp2) - Sta_i
                              State_i(lm1,isp1,lm2,isp2) = State_i(lm1,isp1,lm2,isp2) + Sta_r

!                                  if( Self .and. isp1 == isp2
!     &                                .and. lm1 == lm2 ) then
                              if( Self .and. isp1 == isp2 ) then

                                lm = 0
                                boucle_l3: do l3 = 0,lmax
                                  do m3 = -l3,l3
                                    lm = lm + 1
                                    if( lm > nlm_pot) exit boucle_l3

                                    if( Ylm_mod ) then
                                      g = gauntcp(l1,m1,l2,m2,l3,m3)
                                    else
                                      g = Gaunt_r(l1,m1,l2,m2,l3,m3)
                                    endif

                                    if( abs( g ) < eps10 ) cycle

! plus tard, faire la multiplication par g
                                    if( lm == 1 ) g = 1._db
                                    do ir = 1,nrmtsd
!                                          drho(ir,lm,isr)
!     &                                      = drho(ir,lm,isr)
!     &                                      - g * aimag( c_harm * fc(ir)
!     &                                      * Taull(lma1,is1,lma2,is2) )

                                      Sta_e = c_harm * fc(ir) * Taull(lma1,is1,lma2,is2)
                                      Sta_s = conjg(c_harm) * fc(ir) * Taull(lma2,is2,lma1,is1)

                                      Sta_i = 0.5_db * aimag( Sta_e + Sta_s)

                                      drho(ir,lm,isr) = drho(ir,lm,isr) - g * Sta_i

                                    end do
                                  end do

                                end do boucle_l3

                              endif

                              if( icheck > 2 .and. abs( Taull(lma1,is1,lma2,is2) ) > 1.e-15_db ) then
                                if( nlm2 > 1 ) then
                                  write(3,140) l1, mr1, m1, isp1, l2, mr2, m2, isp2, lp1, mp1, &
                                    iso1, lp2, mp2, iso2, rof_sd, Taull(lma1,is1,lma2,is2), State(lm1,isp1,lm2,isp2), &
                                    State_i(lm1,isp1,lm2,isp2)
                                elseif( nspino == 2 ) then
                                  write(3,150) l1, mr1, m1, isp1, l2, mr2, m2, isp2, iso1, iso2, &
                                    rof_sd, Taull(lma1,is1,lma2,is2), State(lm1,isp1,lm2,isp2), State_i(lm1,isp1,lm2,isp2)
                                else
                                  write(3,160) l1, mr1, m1, isp1, l2, mr2, m2, isp2, rof_sd, Taull(lma1,is1,lma2,is2), &
                                    State(lm1,isp1,lm2,isp2), State_i(lm1,isp1,lm2,isp2)
                                endif
                              endif

                            end do  ! boucle iso2
                          end do  ! boucle mp2
                        end do  ! boucle lp2
                      end do  ! boucle iso1
                    end do  ! boucle mp1
                  end do  ! boucle lp1

                end do  ! boucle isg2
              end do  ! boucle isg1
            end do  ! boucle mr2
          end do  ! boucle l2
        end do  ! boucle isp2
      end do  ! boucle mr1
    end do  ! boucle l1
  end do  ! boucle isp1

  if( icheck > 1 ) then
    write(3,170)
    do l1 = 0,lmax
      if( .not. Full_potential .and. l /= l1 ) cycle
      do mr1 = -l1,l1
        lm1 = lm0 + 1
        do l2 = 0,lmax
          if( .not. Full_potential .and. l /= l2 ) cycle
          do mr2 = -l2,l2
            lm2 = lm0 + 1
            if( Diagonale .and. mr1 /= mr2 ) cycle
            do isp1 = 1,nspinp
              do isp2 = 1,nspinp
                if( nspino == 1 .and. isp1 /= isp2 ) cycle
                write(3,180) l1, mr1, isp1, l2, mr2, isp2, State(lm1,isp1,lm2,isp2), State_i(lm1,isp1,lm2,isp2)
              end do
            end do
          end do
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
    ' l1 mr1 m1 isp1 l2 mr2 m2 isp2 lp1 mp1 iso1 lp2 mp2 iso2',12x,'rof_sd',23x,'Taull',17x,'State',7x,'State_i')
  120 format(/' Radial integral for the density of state:',/ &
              ' l1 mr1 m1 isp1 l2 mr2 m2 isp2 iso1 iso2',12x,'rof_sd',23x,'Taull',17x,'State',7x,'State_i')
  130 format(/' Radial integral for the density of state:',/ &
              ' l1 mr1 m1 isp1 l2 mr2 m2 isp2',12x,'rof_sd',23x,'Taull',17x,'State',7x,'State_i')
  140 format(2(2i3,2i4,1x),2(i3,i4,i5,1x),1p,3(2x,2e13.5))
  150 format(2(2i3,2i4,1x),2(i4,1x),1p,3(2x,2e13.5))
  160 format(2(2i3,2i4,1x),1p,3(2x,2e13.5))
  170 format(/' Radial integral for the density of state before singul:',/ &
              ' l1 m1 isp1 l2 m2 isp2',7x,'State',7x,'State_i')
  180 format(2(2i3,i4,1x),1p,3(2x,2e13.5))
  190 format(/' drho before singul', /'   Radius   ',100(6x,i2,5x))
  200 format(/' drho before singul', /'   Radius   ',100(3x,i2,' up',8x,i2,' dn',5x))
  210 format(1p,300e13.5)
end

!***********************************************************************

! Calcul de la solution singuliere pour la densite d'etat

Subroutine Cal_Solsing_sd(Classic_irreg,drho,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,ll, &
          lmax,m_hubb,nlm,nlm_pot,nlm1,nlm2,nlma,nr,nrmtsd,nrs,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Self, &
          Spinorbite,State,Tau,ui,ur,V,V_hubb,Ylm_mod)

  use declarations
  implicit none

  integer:: icheck, ir, isg, isol, isp, isr, l, ll, lm, lmax, lp, m, m_hubb, mp, mr, ms, mv, nlm, nlm_pot, nlm1, nlm2, n, n0, &
    nlma, np, nr, nr1, nrmtsd, nrs, nspin, nspino, nspinp, numat

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlm1,nspinp,nlm1,nspinp):: Tau
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  logical:: Classic_irreg, Ecomp, Full_potential, Hubb_a, Hubb_d, Radial_comp, Self, Spinorbite, Ylm_mod

  real(kind=db):: c_harm, Eimag, f_integr3, Rmtsd, Sing
  real(kind=db), dimension(nr):: f2, r
  real(kind=db), dimension(nrmtsd):: fct, ri2, rr
  real(kind=db), dimension(nr,nspin):: g0, gm, gp
  real(kind=db), dimension(nr,nspino):: gso
  real(kind=db), dimension(nr,nlm1,nlm2,nspinp,nspino):: ui, ur
  real(kind=db), dimension(nrs,nlm_pot,nspin):: drho
  real(kind=db), dimension(nr,nlm,nlm,nspin):: V
  real(kind=db), dimension(nlma,nspinp,nlma,nspinp):: State
  
  real(kind=db), dimension(:,:,:,:,:), allocatable:: usi, usr

  rr(1:nrmtsd) = r(1:nrmtsd)

  allocate( usi(nrmtsd+1,nlm1,nlm2,nspinp,nspino) )
  allocate( usr(nrmtsd+1,nlm1,nlm2,nspinp,nspino) )

! gm et gp inverse dans le sousprogramme
  call Sch_radial_solsing(Classic_irreg,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde, &
       ll,lmax,m_hubb,nlm,nlm1,nlm2,nr,nrmtsd,nspin,nspino,nspinp,numat,r,Radial_comp,Rmtsd,Spinorbite,Tau,usi,usr,V,V_hubb)

  if( icheck == 2 ) write(3,110)

  ri2(1:nrmtsd) = 1 / ( quatre_pi * r(1:nrmtsd)**2 )

! Ici on prend r = r', car on calcule la trace de G(r,r')
! On ne calcule que la partie imaginaire.
  nr1 = 0
  do l = 0,lmax
    if( .not. Full_potential .and. l /= ll ) cycle
    n0 = l**2 + l + 1
    do mr = -l,l
      if( nlm1 == 1 .and. mr /= 0 ) cycle
      nr1 = nr1 + 1

      do isp = 1,nspinp
        isr = min(isp, nspin)

        fct(:) = 0._db

        do isg = -1,1,2
          if( Ylm_mod .and. mr /= 0 ) then
            c_harm = 0.5_db
            m = isg * mr
            if( isg == 1 ) then
              n = nr1
            else
              n = nr1 - 2 * mr
            endif
          else
            if( isg == 1 ) cycle
            c_harm = 1._db
            m = mr
            n = nr1
          endif

        np = 0
        do lp = 0,lmax
          if( .not. Full_potential .and. lp /= l ) cycle
          do mp = -lp,lp
            if( nlm2 == 1 .and. mp /= m ) cycle
            np = np + 1

            do isol = 1,nspino
              if( Spinorbite .and. nlm2 == 1 ) then
                ms = m + isol - isp
                if( ms > l .or. ms < -l ) cycle
              endif

              do ir = 1,nrmtsd
                if( Radial_comp ) then
                  fct(ir) = fct(ir) + c_harm * r(ir) * ( ur(ir,n,np,isp,isol) * usi(ir,n,np,isp,isol) &
                                                       + ui(ir,n,np,isp,isol) * usr(ir,n,np,isp,isol) )
                else
                  fct(ir) = fct(ir) + c_harm * r(ir) * ur(ir,n,np,isp,isol) * usi(ir,n,np,isp,isol)
                endif
              end do

            end do
          end do
        end do
      end do

        if( icheck > 2 .and. Self ) then
          write(3,115) l, m, isp
          do ir = 1,nrmtsd
            write(3,180) r(ir)*bohr, fct(ir)
          end do
        endif

        Sing = f_integr3(rr,fct,1,nrmtsd,Rmtsd)
        if( Self ) fct(1:nrmtsd) = fct(1:nrmtsd) * ri2(1:nrmtsd)

        if( nlm1 == 1 ) then
          do mv = -l,l
            lm = n0 + mv
            State(lm,isp,lm,isp) = State(lm,isp,lm,isp) + Sing
            if( Self ) drho(1:nrmtsd,1,isr) = drho(1:nrmtsd,1,isr) + fct(1:nrmtsd)
            if( icheck > 2 ) write(3,110)
            if( icheck > 1 ) write(3,120) l, m, isp, Sing, State(lm,isp,lm,isp)
          end do
        else
          lm = n0 + m
          State(lm,isp,lm,isp) = State(lm,isp,lm,isp) + Sing
          if( Self ) drho(1:nrmtsd,1,isr) = drho(1:nrmtsd,1,isr) + fct(1:nrmtsd)
          if( icheck > 2 ) write(3,110)
          if( icheck > 1 ) write(3,120) l, m, isp, Sing, State(lm,isp,lm,isp)
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
  110 format(/' Radial integral of singular solution for the density of state:',/'  l  m isp      Sing         State')
  115 format(/' Singular radial function for the density of state, l, m, isp =',3i3)
  120 format(3i3,1p,1x,e13.5,2x,2e13.5)
  170 format(/' drho after singul', /'   Radius         up           dn')
  180 format(1p,3e13.5)
end

!***********************************************************************

! Calculation of the number of couple of atoms where one calculate the Crystal Orbital Overlap 

function nab_coop_ev(COOP,Dist_coop,ia_coop,iaprotoi,n_atom_coop,n_atom_proto,natome,Posi,Rmtg)

  use declarations
  implicit none

  integer:: ia, iab, ib, j, n_atom_coop, nab_coop_ev, natome, ipra, iprb, n_atom_proto

  integer, dimension(n_atom_coop):: ia_coop
  integer, dimension(natome):: iaprotoi
  
  logical:: COOP
  
  real(kind=db):: Dist, Dist_coop
  real(kind=db), dimension(0:n_atom_proto):: Rmtg
  real(kind=db), dimension(3,natome):: Posi
  
  iab = 0
  
  if( COOP ) then
    do ia = 1,natome
      ipra = iaprotoi( ia )

      do ib = ia+1,natome
        do j = 1,n_atom_coop
          if( ia_coop(j) == 0 .or. ia == ia_coop(j) .or. ib == ia_coop(j) ) exit
        end do
        if( j > n_atom_coop ) cycle
 
        iprb = iaprotoi( ib )
        Dist = sqrt( sum( ( Posi(:,ia) - Posi(:,ib) )**2 ) )
        if( Dist_coop > eps10 ) then
          if( Dist > Dist_coop + eps10 ) cycle
        else
          if( Dist > Rmtg( ipra ) + Rmtg( iprb ) ) cycle
        endif  
        iab = iab + 1
      end do
    end do
  endif
  
  nab_coop_ev = iab

  return
end

!***********************************************************************

! Calculation of the Crystal Orbital Overlap Population

subroutine Cal_COOP(Density_comp,Dist_coop,Ecinetic,Eimag,Energ,Enervide,Full_atom,Full_potential,Hubb, &
            Hubb_diag,ia_coop,iabsorig,iaprotoi,icheck,ie,itypepr,lmax_pot,lmaxat,m_hubb,n_atom_0,n_atom_0_self, &
            n_atom_coop,n_atom_ind,n_atom_ind_self,n_atom_proto,n_multi_run,nab_coop,natome,nlm_pot,nlmagm,nomfich_s,nrato,nrm, &
            nspin,nspino,nspinp,ntype,numat,Posi,rato,Relativiste,Renorm,Rmtg,Rmtg0,Rmtsd,Rot_atom,Spinorbite,Tau_coop,V_hubb, &
            V_intmax,V0bd,Vrato,Ylm_comp,Z_nospinorbite)

  use declarations
  implicit none

  integer, parameter:: Zmax = 103
  
  integer:: i, i_rad, ia, iab, iabsorig, iapr, iapra, iaprb, ib, icheck, ie, ipr, ipra, iprb, ira, irb, isa, isb, &
    isoa, isob, ispa, ispb, it, j, l, la, lam, lb, lbm, lm, lm0a, lm0b, lma, lmax, lmaxa, lmaxb, lmaxm, lmax_pot, lmb, lmm, lmp, &
    m_hubb, ma, mb, mpa, mpb, mva, mvb, n, n_angle, n_atom_0, n_atom_0_self, n_atom_coop, n_atom_ind, n_atom_ind_self, &
    n_atom_proto, n_multi_run, n_rad, n_rad_m, nab_coop, natome, nlm_pot, nlmagm, nlmc, nlmr, nr, nrm, nspin, nspino, nspinp, &
    ntype, Z, Z_nospinorbite, Za, Zb

  integer, dimension(n_atom_coop):: ia_coop
  integer, dimension(natome):: iaprotoi
  integer, dimension(n_atom_0:n_atom_ind):: n_radius
  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, lmaxat
  integer, dimension(Zmax):: lmax_coop
  integer, dimension(:), allocatable:: indexa, indexb

  character(len=132):: nomfich_coop, nomfich_s
  character(len=2), dimension(-7:99):: ch

  complex(kind=db):: Yc, Ycompa, Ycompb, YlmaYlmb, YuaYub
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb
  complex(kind=db), dimension(nlmagm,nspinp,nlmagm,nspinp,nab_coop):: Tau_coop, Tau_coop_r
  complex(kind=db), dimension(:), allocatable:: Ylmc
  complex(kind=db), dimension(:,:,:), allocatable:: Dlmm_ia, Dlmm_ib
  complex(kind=db), dimension(:,:,:,:), allocatable:: Taull, V_hubb_t
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: um
  
  logical:: Density_comp, Full_atom, Full_potential, Hubb_a, Hubb_d, Relativiste, Renorm, Sens, Spino, Spinorbite, Ylm_comp 
  logical, dimension(0:ntype):: Hubb
  logical, dimension(n_atom_0_self:n_atom_ind_self):: Hubb_diag

  real(kind=db):: Angle, Dist, Dist_coop, dOverlap, Eimag, Energ, Enervide, Rad_max, Radius_coop, Step, Step_rad, V_intmax, Vol
  real(kind=db), dimension(3):: T, U, V, Va, Vb, W
  real(kind=db), dimension(nspin):: Ecinetic, V0bd
  real(kind=db), dimension(nab_coop):: COverlap_T
  real(kind=db), dimension(0:n_atom_proto):: Rmtg, Rmtg0, Rmtsd
  real(kind=db), dimension(3,3):: Rot_a, Rot_b
  real(kind=db), dimension(3,3,natome):: rot_atom
  real(kind=db), dimension(3,natome):: Posi
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vrato 

  real(kind=db), dimension(:), allocatable:: r, Ylmr
  real(kind=db), dimension(:,:), allocatable:: Radius_urm, Ylma, Ylmb 
  real(kind=db), dimension(:,:,:), allocatable:: Vrato_e
  real(kind=db), dimension(:,:,:,:,:), allocatable:: Coverlap, Coverlap_l

  data ch/'-7','-6','-5','-4','-3','-2','-1',',0', &
          ',1',',2',',3',',4',',5',',6',',7',',8',',9','10','11','12','13','14','15','16','17','18','19','20', & 
          '21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40', & 
          '41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60', & 
          '61','62','63','64','65','66','67','68','69','70','71','72','73','74','75','76','77','78','79','80', & 
          '81','82','83','84','85','86','87','88','89','90','91','92','93','94','95','96','97','98','99'/ 

  data lmax_coop/1,                                  1, &   ! He
                 1,1,                      2,2,2,2,2,2, &  ! Ne
                 2,2,                      2,2,2,2,2,2, &  ! Ar
                 2,2, 2,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2, &  ! Kr 36  
                 2,2, 2,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2, &  ! Xe 54   
                 3,3, 3,3,3,3,3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3, & ! Rn 86     
                 3,3, 3,3,3,3,3,3,3,3,3,3,3,3,3,3, 3/ ! Lw 103     
 
  COverlap_T(:) = 0._db
  n_rad = 1
  Rad_max = 0.20_db
  Step_rad = Rad_max / max( n_rad, 1 )
  n_angle = 8
  Step = 2 * pi / max( n_angle, 1 )

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

        do j = 1,n_atom_coop
          if( ia_coop(j) == 0 .or. ia == ia_coop(j) .or. ib == ia_coop(j) ) exit
        end do
        if( j > n_atom_coop ) cycle 

        iprb = iaprotoi( ib )
        if( Full_atom ) then
          iaprb = ib
        else
          iaprb = iprb
        endif

        Dist = sqrt( sum( ( Posi(:,ia) - Posi(:,ib) )**2 ) )
        if( Dist_coop > eps10 ) then
          if( Dist > Dist_coop + eps10 ) cycle
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
    lmaxm = max( lmaxm, lmaxat(ipr) )
    lmaxm = max( lmaxm, lmax_coop( numat( itypepr(ipr) ) ) )
  end do
  lmm = (lmaxm + 1)**2

  allocate( um(lmm,nspinp,nspino,n_rad_m,n_atom_0:n_atom_ind) )
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
      V_Hubb_t(:,:,:,:) = V_Hubb(:,:,:,:,iapr)
      Hubb_d = Hubb_diag(iapr)
    else
      Hubb_a = .false.
      Hubb_d = .true.
    end if

    nr = nrato(it)
    allocate( r(nr) )
    allocate( Vrato_e(nr,nlm_pot,nspin) )
    r(1:nr) = rato(1:nr,it)

    lmax = lmaxat(ipr)

    Vrato_e(1:nr,:,:) = Vrato(1:nr,:,:,iapr)

    call Cal_urm(Ecinetic,Eimag,Energ,Enervide, &
            Full_potential,Hubb_a,Hubb_d,iapr,icheck,lmax,lmax_pot,lmm,m_hubb,n_atom_0,n_atom_ind,n_rad_m,n_radius,nlm_pot, &
            nr,nspin,nspino,nspinp,r,Radius_urm,Relativiste,Renorm,Rmtg(ipr),Rmtsd(ipr),Spinorbite,um,V_hubb, &
            V_intmax,V0bd,Vrato_e,Ylm_comp,Z)

    deallocate( r )
    deallocate( Vrato_e, V_hubb_t )

  end do

  nlmr = ( lmaxm + 1 )**2
  nlmc = ( ( lmaxm + 1 ) * ( lmaxm + 2 ) ) / 2
  allocate( Ylmc(nlmc) )
  allocate( Ylmr(nlmr) )
  allocate( Ylma(nlmr,max( n_rad * n_angle, 1 )) )
  allocate( Ylmb(nlmr,max( n_rad * n_angle, 1 )) )

  allocate( COverlap(nlmr,nspinp,nlmr,nspinp,nab_coop) )
  allocate( COverlap_l(0:lmaxm,nspinp,0:lmaxm,nspinp,nab_coop) )
  COverlap(:,:,:,:,:) = 0._db
  COverlap_l(:,:,:,:,:) = 0._db
  allocate( indexa(nab_coop) )
  allocate( indexb(nab_coop) )
  
  iab = 0
  n_radius(:) = 0
  
  Tau_coop_r(:,:,:,:,:) = (0._db, 0._db)

  do ia = 1,natome
    ipra = iaprotoi( ia )
    Za = numat( itypepr(ipra) )
    lmaxa = min( lmaxat( ipra ), lmax_coop( Za ) )
    Rot_a(:,:) = rot_atom(:,:,ia)
    allocate( Dlmm_ia(-lmaxa:lmaxa,-lmaxa:lmaxa,0:lmaxa) )
    call Rotation_mat(Dlmm_ia,0,Rot_a,lmaxa,Ylm_comp)  ! 0 is for icheck

    do ib = ia+1,natome

      do j = 1,n_atom_coop
        if( ia_coop(j) == 0 .or. ia == ia_coop(j) .or. ib == ia_coop(j) ) exit
      end do
      if( j > n_atom_coop ) cycle 

      iprb = iaprotoi( ib )

      Dist = sqrt( sum( ( Posi(:,ia) - Posi(:,ib) )**2 ) )
      if( Dist_coop > eps10 ) then
        if( Dist > Dist_coop + eps10 ) cycle
      else
        if( Dist > Rmtg( ipra ) + Rmtg( iprb ) ) cycle
      endif  

      iab = iab + 1

      Zb = numat( itypepr(iprb) )
      lmaxb = min( lmaxat( iprb ), lmax_coop( Zb ) )
     
      Rot_b(:,:) = rot_atom(:,:,ib)
      allocate( Dlmm_ib(-lmaxb:lmaxb,-lmaxb:lmaxb,0:lmaxb) )
      call Rotation_mat(Dlmm_ib,0,Rot_b,lmaxb,Ylm_comp)

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
                      Tau_coop_r(lma,isa,lmb,isb,iab) = Tau_coop_r(lma,isa,lmb,isb,iab) &
                                                      + Dlmm_ia(mpa,ma,la) * Tau_coop(lm,isa,lmp,isb,iab) * Dlmm_ib(mpb,mb,lb)
                    end do
                  end do
                end do  
              end do 
            end do
          end do
        end do
      end do
      deallocate( Dlmm_ib )
      if( icheck > 1 .and. ie == 1 ) then
        do isa = 1,nspinp
          do isb = 1,nspinp
            write(3,'(/a26,3i3)') ' Tau_coop, iab, isa, isb =', iab, isa, isb
            lma = 0 
            do la = 0,lmaxa
              do ma = -la,la
                lma = lma + 1 
                write(3,'(2i3,1p,100e13.5)') la, ma, (Tau_coop(lma,isa,lmb,isb,iab), lmb = 1,(lmaxb + 1)**2)
              end do
            end do
            write(3,'(/A)') ' Tau_coop_r' 
            lma = 0 
            do la = 0,lmaxa
              do ma = -la,la
                lma = lma + 1 
                write(3,'(2i3,1p,100e13.5)') la, ma, (Tau_coop_r(lma,isa,lmb,isb,iab), lmb = 1,(lmaxb + 1)**2)
              end do
            end do
          end do
        end do
      endif

      if( ( Ylm_comp .and. .not. Density_comp ) .or. ( .not. Ylm_comp .and. Density_comp ) ) then
        Sens = Ylm_comp .and. .not. Density_comp 
        lma = (lmaxa + 1)**2 
        lmb = (lmaxb + 1)**2 
        allocate( Taull(lma,nspinp,lmb,nspinp) )
        Taull(1:lma,:,1:lmb,:) = Tau_coop_r(1:lma,:,1:lmb,:,iab)
        call Trans_Tau_ab(lmaxa,lmaxb,nspinp,Spinorbite,Taull,Sens)
        Tau_coop_r(1:lma,:,1:lmb,:,iab) = Taull(1:lma,:,1:lmb,:)
        deallocate( Taull )
        if( icheck > 1 .and. ie == 1 ) then
          do isa = 1,nspinp
            do isb = 1,nspinp
              if( Sens ) then
                write(3,'(/A)') ' Tau_coop_r in real Ylm basis'
              else
                write(3,'(/A)') ' Tau_coop_r in complex Ylm basis'
              endif 
              lma = 0 
              do la = 0,lmaxa
                do ma = -la,la
                  lma = lma + 1 
                  write(3,'(2i3,1p,100e13.5)') la, ma, (Tau_coop_r(lma,isa,lmb,isb,iab), lmb = 1,(lmaxb + 1)**2)
                end do
              end do
            end do
          end do
        endif
      endif 
 
    end do
    
    deallocate( Dlmm_ia )

  end do 
  
  iab = 0
  n_radius(:) = 0

  do ia = 1,natome
    ipra = iaprotoi( ia )
    Za = numat( itypepr(ipra) )
    if( Full_atom ) then
      iapra = ia
    else
      iapra = ipra
    endif

    do ib = ia+1,natome

      do j = 1,n_atom_coop
        if( ia_coop(j) == 0 .or. ia == ia_coop(j) .or. ib == ia_coop(j) ) exit
      end do
      if( j > n_atom_coop ) cycle 

      iprb = iaprotoi( ib )
      if( Full_atom ) then
        iaprb = ib
      else
        iaprb = iprb
      endif

      Dist = sqrt( sum( ( Posi(:,ia) - Posi(:,ib) )**2 ) )
      if( Dist_coop > eps10 ) then
        if( Dist > Dist_coop + eps10 ) cycle
      else
        if( Dist > Rmtg( ipra ) + Rmtg( iprb ) ) cycle
      endif  

      iab = iab + 1

      Zb = numat( itypepr(iprb) )

      indexa(iab) = ia
      indexb(iab) = ib

      Radius_coop = ( Rmtg0( ipra ) + Rmtg0( iprb ) ) / 2 
      Vol = ( 4 * pi / 3._db ) * Radius_coop**3

      if( icheck > 1 ) then 
        write(3,'(/3(a6,i3),a7,f10.5,a3)') ' iab =',iab,', ia =',ia,', ib =', ib,', Vol =', Vol, ' ua'
        write(3,'(/A,A)') '  ir  la ma  lb mb         Sum_Ya.Yb           um_r(ia)     um_i(ia)     um_r(ib)     um_i(ib)', &
                          '           Tau_coop            COverlap'
      endif

! Construction of the basis (Va (or Vb), U, T) to build n_angle points around a circle, + 1 at the center      
      V(:) = Posi(:,ib) - Posi(:,ia) 
      V(:) = V(:) / sqrt( sum( V(:)**2 ) )
      W(:) = 0._db
      if( abs( V(3) ) < 0.9_db ) then
        W(3) = 1._db
      else
        W(1) = 1._db
      endif
      call Prodvec(U,V,W)
      U(:) = Step_rad * Radius_coop * U(:) / sqrt( sum( U(:)**2 ) )
      call Prodvec(T,U,V)

      n = 1
      do i_rad = 1,n_rad
        n = n + i_rad * n_angle  
      end do
      Vol = Vol / n

      Va(:) = Rmtg0( ipra ) * V(:)
      Vb(:) = - Rmtg0( iprb ) * V(:)

      do i_rad = 0,n_rad
        n_radius(iapra) = n_radius(iapra) + 1
        n_radius(iaprb) = n_radius(iaprb) + 1

! One makes the average over i_rad * n_angle      
        do i = 1,max( i_rad * n_angle, 1 )
          Angle = i * Step / i_rad       

          if( i_rad == 0 ) then
            W(:) = Va(:)
          else  
            W(:) = Va(:) + i_rad * ( cos( Angle ) * U(:) + sin( Angle ) * T(:) ) 
          endif
          ira = n_radius(iapra) 
       
          Dist = sqrt( sum( W(:)**2 ) )
          call cYlm(lmaxm,W,Dist,Ylmc,nlmc)
! Conversion in real Ylm
          call Ylmcr(lmaxm,nlmc,nlmr,Ylmc,Ylmr)
          Ylma(1:nlmr,i) = Ylmr(1:nlmr)

          if( i_rad == 0 ) then
            W(:) = Vb(:)
          else  
            W(:) = Vb(:) + i_rad * ( cos( Angle ) * U(:) + sin( Angle ) * T(:) ) 
          endif
          irb = n_radius(iaprb) 
       
          Dist = sqrt( sum( W(:)**2 ) )
          call cYlm(lmaxm,W,Dist,Ylmc,nlmc)
! Conversion in real Ylm
          call Ylmcr(lmaxm,nlmc,nlmr,Ylmc,Ylmr)
          Ylmb(1:nlmr,i) = Ylmr(1:nlmr)
        end do

        do isa = 1,nspinp
          lma = 0
          do la = 0,min( lmaxat( ipra ), lmax_coop( Za ) )
            do ma = -la,la
              lma = lma + 1
              do ispa = 1,nspinp
                if( .not. Spinorbite .and. ispa /= isa ) cycle

                isoa = min( ispa, nspino )

                if( Spinorbite ) then
                  mva = ma - isa + isoa
                  if( mva > la .or. mva < -la ) cycle
                else
                  mva = ma
                endif
                if( Density_comp ) lm0a = la**2 + la + 1 + mva
              
                do isb = 1,nspinp
                  if( .not. Spinorbite .and. isa /= isb ) cycle
                  lmb = 0
                  do lb = 0,min( lmaxat( iprb ), lmax_coop( Zb ) )
                    do mb = -lb,lb
                      lmb = lmb + 1
         
                      do ispb = 1,nspinp
                        if( .not. Spinorbite .and. ispb /= isb ) cycle

                        isob = min( ispb, nspino )
  
                        if( Spinorbite ) then
                          mvb = mb - isb + isob
                          if( mvb > lb .or. mvb < -lb ) cycle
                        else
                          mvb = mb
                        endif

                        YlmaYlmb = (0._db, 0._db )
                        do i = 1,max( i_rad * n_angle, 1 )
                          if( Density_comp ) then
                            if( mva /= 0 ) then
                              Ycompa = Yc(mva,Ylma(lm0a,i),Ylma(lm0a-2*mva,i))
                            else
                              Ycompa = cmplx( Ylma(lm0a,i), 0._db, db )
                            endif
                            lm0b = lb**2 + lb + 1 + mvb
                            if( mvb /= 0 ) then
                              Ycompb = Yc(mvb,Ylmb(lm0b,i),Ylmb(lm0b-2*mvb,i))
                            else
                              Ycompb = cmplx( Ylmb(lm0b,i), 0._db, db )
                            endif
                            YlmaYlmb = YlmaYlmb + conjg( Ycompa ) * Ycompb
                          else
                            YlmaYlmb = YlmaYlmb + cmplx( Ylma(lma,i) * Ylmb(lmb,i), 0._db ) 
                          endif
                        end do
 
                        YuaYub = YlmaYlmb * um(lma,isa,isoa,ira,iapra) * um(lmb,isb,isob,irb,iaprb)
                        dOverlap = - aimag( YuaYub * Tau_coop_r(lma,ispa,lmb,ispb,iab) ) * Vol
                        COverlap(lma,isa,lmb,isb,iab) = COverlap(lma,isa,lmb,isb,iab) + dOverlap
                        COverlap_l(la,isa,lb,isb,iab) = COverlap_l(la,isa,lb,isb,iab) + dOverlap 
                        COverlap_T(iab) = COverlap_T(iab) + dOverlap
                    
                        if( icheck > 1 ) &
                          write(3,105) i_rad, la, ma, lb, mb, YlmaYlmb, um(lma,isa,isoa,ira,iapra), um(lmb,isb,isob,irb,iaprb), &
                                     Tau_coop_r(lma,ispa,lmb,ispb,iab), COverlap(lma,isa,lmb,isb,iab) 
                      end do
                    end do
                  end do
                end do
              end do
       
            end do
          end do

        end do ! end of loop over n_angle
      end do ! end of loop over n_rad
            
    end do
  end do
 
 ! Writing
  do iab = 1,nab_coop
    nomfich_coop = nomfich_s
    if( n_multi_run > 1 ) then
      l = len_trim(nomfich_coop)
      nomfich_coop(l+1:l+1) = '_'
      call ad_number(iabsorig,nomfich_coop,132)
    endif
    l = len_trim(nomfich_coop)
    nomfich_coop(l+1:l+1) = '_'
    call ad_number(indexa(iab),nomfich_coop,132)
    l = len_trim(nomfich_coop)
    nomfich_coop(l+1:l+1) = '_'
    call ad_number(indexb(iab),nomfich_coop,132)
    l = len_trim(nomfich_coop)
    nomfich_coop(l+1:l+9) = '_coop.txt'
 
    ia = indexa(iab)
    ib = indexb(iab)
    it = itypepr( iaprotoi(ia) ) 
    lam = lmax_coop( numat( it ) )
    it = itypepr( iaprotoi(ib) ) 
    lbm = lmax_coop( numat( it ) )
 
    if( ie == 1 ) then
      open(4, file = nomfich_coop )
      if( nspinp == 2 ) then
        write(4,110) ch(ia), ch(ib), (((( la, '..', isa, lb, '..', isb, &
                     la = 0,lam), isa = 1,nspinp), lb = 0,lbm), isb = 1,nspinp), &
                     (((((( la, ch(ma), isa, lb, ch(mb), isb, &
                      ma = -la,la), la = 0,lam), isa = 1,nspinp), mb = -lb,lb), lb = 0,lbm), isb = 1,nspinp)
      else
        write(4,120) ch(ia), ch(ib), (( la, '..', lb, '..', la = 0,lam), lb = 0,lbm), & 
                     (((( la, ch(ma), lb, ch(mb), ma = -la,la), la = 0,lam), mb = -lb,lb), lb = 0,lbm) 
      endif
    else
      open(4, file = nomfich_coop, position='append')
    endif
    write(4,130) Energ*rydb, COverlap_T(iab), & 
                     (((( COverlap_l(la,isa,lb,isb,iab), la = 0,lam), isa = 1,nspinp), lb = 0,lbm), isb = 1,nspinp), & 
                     (((( COverlap(lma,isa,lmb,isb,iab), &
                     lma = 1,( 1 + lam )**2), isa = 1,nspinp), lmb = 1,( 1 + lbm )**2), isb = 1,nspinp) 
    Close(4)
  end do
  
  deallocate( COverlap, Radius_urm, um, Ylma, Ylmb, Ylmc, Ylmr )
  
  return
  105 format(i4,2(1x,2i3),1p,10e13.5)
  110 format('   Energy  Total_COOP',2a2,' ',1000('  (',i1,a2,i1,',',i1,a2,i1,')  '))
  120 format('   Energy  Total_COOP',2a2,' ',1000('   (',i1,a2,',',i1,a2,')   '))
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

  integer:: i, iapr, icheck, ir, l, l_hubbard, lfin, lm, lm0, lmax, lmax_pot, lmm, m, m_hubb, mp, n, n_atom_0, &
    n_atom_ind, n_rad_m, nlm, nlm_pot, nlm1, nlm2, np, nr, nrmtg, nrmtsd, nspin, nspino, nspinp, Z

  integer, dimension(n_atom_0:n_atom_ind):: n_radius
  integer, dimension(n_rad_m):: index_r

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(lmm,nspinp,nspino,n_rad_m,n_atom_0:n_atom_ind):: um
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
    
  do l = 0,lfin

    if( Hubb_a .and. l == l_hubbard( Z ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = nlm   ! = (lmax + 1 )**2
      nlm2 = nlm
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*l + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*l + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,l,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,Z,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

    do i = 1,n_radius(iapr)

      ir = index_r(i)
      if( Radial_comp ) then
        if( Full_potential ) then
          do lm = 1,nlm1
            um(lm,:,:,i,iapr) = cmplx( p(i) * ur(ir,lm,lm,:,:) + (1 - p(i)) * ur(ir-1,lm,lm,:,:), &
                                       p(i) * ui(ir,lm,lm,:,:) + (1 - p(i)) * ui(ir-1,lm,lm,:,:), db )
          end do
        else
          lm0 = l**2 + l + 1
          do m = -l,l
            if( nlm1 == 1 ) then
              n = 1
            else
              n = l + 1 + m
            endif
            do mp = -l,l
!              if( nlm2 == 1 .and. mp /= m ) cycle
              if( mp /= m ) cycle
              np = min(l + m + 1, nlm2 )
              um(lm0+m,:,:,i,iapr) = cmplx( p(i) * ur(ir,n,np,:,:) + (1 - p(i)) * ur(ir-1,n,np,:,:), &
                                            p(i) * ui(ir,n,np,:,:) + (1 - p(i)) * ui(ir-1,n,np,:,:), db )
            end do
          end do
        endif
      else
        if( Full_potential ) then
          do lm = 1,nlm1
            um(lm,:,:,i,iapr) = cmplx( p(i) * ur(ir,lm,lm,:,:) + (1 - p(i)) * ur(ir-1,lm,lm,:,:), 0._db, db )
          end do
        else
          lm0 = l**2 + l + 1
          do m = -l,l
            if( nlm1 == 1 ) then
              n = 1
            else
              n = l + 1 + m
            endif
            do mp = -l,l
!              if( nlm2 == 1 .and. mp /= m ) cycle
              if( mp /= m ) cycle
              np = min(l + m + 1, nlm2 )
              um(lm0+m,:,:,i,iapr) = cmplx( p(i) * ur(ir,n,np,:,:) + (1 - p(i)) * ur(ir-1,n,np,:,:), 0._db, db )
            end do
          end do
        endif
      endif
    
    end do
    
    deallocate( Tau, ui, ur )

  end do   ! end of loop over l

  if( icheck > 1 ) then
    write(3,'(/a34,i3)') ' Radius, (um, l = 0,lmax),   atom', iapr
    do i = 1,n_radius(iapr) 
      write(3,160) Radius_urm(i,iapr)*bohr, ( um(l**2+l+1,1,1,i,iapr), l = 0,lmax )
    end do 
  endif
  
  deallocate( V )

  return
  110 format(/' iapr =',i3,', Z =',i3,', lmax =',i2)
  120 format(' Energ =',f10.3,' eV,  Enervide =',f10.3,' eV')
  130 format(' Ecinetic =',2f10.3)
  140 format(' V0bd  =',2f10.3)
  150 format(' konde =',2f12.5)
  160 format(' Radius =',f8.5,' A, um =',1p,10e13.5)
end

!***********************************************************************

subroutine Radial_wave(Ecinetic,Eimag,Energ,Enervide,Full_potential,Hubb_a,Hubb_d,icheck,initl,lmax,lmax_pot, &
            m_hubb,nbseuil,ninit1,ninitlt,nlm_pot,nlmam,nlmam2,nr,nrm,nspin,nspino,nspinp,numat,psii,r,Radial_comp, &
            Relativiste,Renorm,Rmtg,Rmtsd,rof_ph,Spinorbite,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp,zet)

  use declarations
  implicit none

  integer:: icheck, initl, iseuil, iso, isp, l, l_hubbard, lfin, lm, lm0, lmax, lmax_pot, lmp, m, m_hubb, mp, &
    n, nbseuil, ninit1, ninitlt, nlm, nlm_pot, nlm1, nlm2, nlmam, nlmam2, np, nr, nrm, nrmtg, nrmtsd, nspin, nspino, nspinp, numat

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nlmam,nlmam2,nspinp,nspino,0:0,1):: rof
  complex(kind=db), dimension(nlmam,nspinp,nspino,ninitlt):: rof_ph
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
  real(kind=db), dimension(nr,nlmam,nlmam2,nspinp,nspino,ninitlt):: zet

  real(kind=db), dimension(:,:,:,:), allocatable:: V
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur

  NRIXS = .false.
  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

  if( nbseuil == 2 ) then
    if( nbseuil /= ninitlt ) then
      ninit1 = ( ninitlt - 2 ) / 2
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

  do l = 0,lfin

    if( Hubb_a .and. l == l_hubbard( numat ) )  then
      Hubb_m = .true.
    else
      Hubb_m = .false.
    endif

    if( Full_potential ) then
      nlm1 = ( lmax + 1 )**2
      nlm2 = nlm1
    elseif( Hubb_m .and. .not. Hubb_d ) then
      nlm1 = 2*l + 1
      nlm2 = nlm1
    elseif( Spinorbite .or. Hubb_m ) then
      nlm1 = 2*l + 1
      nlm2 = 1
    else
      nlm1 = 1
      nlm2 = 1
    endif

    allocate( ui(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( ur(nr,nlm1,nlm2,nspinp,nspino) )
    allocate( Tau(nlm1,nspinp,nlm1,nspinp) )

    call Sch_radial(Ecinetic,Ecomp,Eimag,f2,Full_potential,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,l,lmax,m_hubb, &
         nlm,nlm1,nlm2,nr,nrmtg,nspin,nspino,nspinp,numat,r,Radial_comp,Relativiste,Renorm,Rmtg,Spinorbite,Tau,ui, &
         ur,V,V_hubb)

! Integrale radiale pour la regle d'or de Fermi
    call radial_matrix(.true.,1,0,0,iseuil,l,nlm1,nlm2,nbseuil, &
           1,nlmam,nlmam2,nr,NRIXS,nrm,nrmtsd,nspino,nspinp,psii,r,r_or_bess,Radial_comp,Rmtsd,rof,ui,ur,Vecond)

! Recopie
    if( Full_potential ) then
      zet(1:nr,1:nlm1,1:nlm2,:,:,initl) = ur(1:nr,1:nlm1,1:nlm2,:,:)
    else
      lm0 = l**2 + l + 1
      do m = -l,l
        if( nlm1 == 1 ) then
          n = 1
        else
          n = l + 1 + m
        endif
        do mp = -l,l
          if( nlm2 == 1 .and. mp /= m ) cycle
          np = min(l + m + 1, nlm2 )
          lmp = min(lm0 + mp, nlmam2 )
          zet(1:nr,lm0+m,lmp,:,:,initl) = ur(1:nr,n,np,:,:)
        end do
      end do
    endif

    deallocate( Tau )

    deallocate( ui, ur )

  end do   ! fin de la boucle sur l

  if( nlm2 == 1 ) then
    rof_ph(:,:,:,initl) = rof(:,1,:,:,0,1)
  else
    do iso = 1,nspino
      do isp = 1,nspinp
        do lm = 1,nlmam2
          rof_ph(lm,isp,iso,initl) =  sum( rof(:,lm,isp,iso,0,1) )
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

! Calcul l'integrale de 0 a r (is=1) ou r a Rmtsd (is=-1) de fct
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

! Calcul l'integrale de 0 a r (is=1) ou r a Rmtsd (is=-1) de fct

subroutine ffintegr2_r(fint,fct,r,n,is,Rmtsd)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  real(kind=db):: a, a_tiers, b, b_demi, c, f0, fm, fp
  real(kind=db), dimension(2):: dintegr
  real(kind=db), dimension(n):: fct, fint

  real(kind=db), dimension(n):: r

  tiers = 1._db / 3._db

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

! Transformation harmo comp vers Harmo reel si Sens = .true.
!                                 sinon inverse

subroutine Trans_Tau(icheck,lmin,lmax,nspin,Spinorbite,Taull,Sens)

  use declarations
  implicit none

  integer:: icheck, isp1, isp2, l, lm0, lmin, lmax, m, nspin

  complex(kind=db), dimension(:,:), allocatable:: Tau, Trans
  complex(kind=db), dimension((lmax+1)**2-lmin**2,nspin,(lmax+1)**2-lmin**2,nspin):: Taull

  logical:: Sens, Spinorbite

  do isp1 = 1,nspin
    do isp2 = 1,nspin
      if( .not. Spinorbite .and. isp1 /= isp2 ) cycle
      do l = lmin,lmax

        if( l == 0 ) cycle

        if( lmin == lmax ) then
          lm0 = l + 1
        else
          lm0 = l**2 + l + 1
        endif
        allocate( Tau(-l:l,-l:l) )
        allocate( Trans(-l:l,-l:l) )
        do m = -l,l
          Tau(-l:l,m) = Taull(lm0-l:lm0+l,isp1,lm0+m,isp2)
        end do

        call cal_trans_l(l,Trans)

        if( Sens ) then
          Tau = matmul( Trans, matmul( Tau, Transpose( conjg(Trans) ) ) )
        else
          Tau = matmul( Transpose( conjg(Trans) ), matmul( Tau, Trans ) )
        endif

        do m = -l,l
          Taull(lm0-l:lm0+l,isp1,lm0+m,isp2) = Tau(-l:l,m)
        end do

        if( icheck > 1 ) then
          if( nspin == 1 ) then
            if( Sens ) then
              write(3,110) l
            else
              write(3,120) l
            endif
          elseif( Spinorbite ) then
            if( Sens ) then
              write(3,130) l, isp1, isp2
            else
              write(3,140) l, isp1, isp2
            endif
          else
            if( Sens ) then
              write(3,150) l, isp1
            else
              write(3,160) l, isp1
            endif
          endif
          do m = -l,l
            write(3,170) m, Tau(m,:)
          end do
        endif

        deallocate( Tau, Trans )

      end do
    end do
  end do

  110 format(/' Tau in real basis, l =',i2)
  120 format(/' Tau in complex basis, l =',i2)
  130 format(/' Tau in real basis, l =',i2,',isp1, isp2 =',2i2)
  140 format(/' Tau in complex basis, l =',i2,',isp1, isp2 =',2i2)
  150 format(/' Tau in real basis, l =',i2,',isp =',i2)
  160 format(/' Tau in complex basis, l =',i2,',isp =',i2)
  170 format(i3,1p,7(1x,2e13.5))
  return
end

!***********************************************************************

! Transformation harmo comp vers Harmo reel si Sens = .true.
!                                 sinon inverse

subroutine Trans_Tau_ab(lmaxa,lmaxb,nspin,Spinorbite,Taull,Sens)

  use declarations
  implicit none

  integer:: isp1, isp2, la, lb, lm0a, lm0b, lmaxa, lmaxb, nspin

  complex(kind=db), dimension(:,:), allocatable:: Tau_ab, Trans_a, Trans_b
  complex(kind=db), dimension((lmaxa+1)**2,nspin,(lmaxb+1)**2,nspin):: Taull

  logical:: Sens, Spinorbite

  do isp1 = 1,nspin
    do isp2 = 1,nspin
      if( .not. Spinorbite .and. isp1 /= isp2 ) cycle

      do la = 0,lmaxa

        lm0a = la**2 + la + 1
        allocate( Trans_a(-la:la,-la:la) )
        call Cal_trans_l(la,Trans_a)

        do lb = 0,lmaxb

          lm0b = lb**2 + lb + 1
          allocate( Trans_b(-lb:lb,-lb:lb) )
          call Cal_trans_l(lb,Trans_b)

          allocate( Tau_ab(-la:la,-lb:lb) )

          Tau_ab(-la:la,-lb:lb) = Taull(lm0a-la:lm0a+la,isp1,lm0b-lb:lm0b+lb,isp2)

          if( Sens ) then
            Tau_ab = matmul( Trans_a, matmul( Tau_ab, Transpose( conjg(Trans_b) ) ) )
          else
            Tau_ab = matmul( Transpose( conjg(Trans_a) ), matmul( Tau_ab, Trans_b ) )
          endif

          Taull(lm0a-la:lm0a+la,isp1,lm0b-lb:lm0b+lb,isp2) = Tau_ab(-la:la,-lb:lb) 

          deallocate( Tau_ab, Trans_b )
        
        end do  

        deallocate( Trans_a )

      end do

    end do
  end do

  return
end

!***********************************************************************

! Transformation harmo reel vers Harmo comp pour un lh donne.
! La transformation inverse est le conjugue de la transpose

subroutine Cal_Trans_l(lh,Trans)

  use declarations
  implicit none

  integer:: is, lh, m1, m2

  complex(kind=db):: r2_r, r2_i
  complex(kind=db),dimension(-lh:lh,-lh:lh):: Trans

  real(kind=db):: r2

  Trans(:,:) = (0._db, 0._db)

  r2 = 1 / sqrt(2._db)
  r2_r = cmplx( r2,    0._db, db)
  r2_i = cmplx( 0._db, r2,    db)

  do m1 = -lh,lh
    is = (-1)**m1
    do m2 = -lh,lh

      if( m1 == m2 ) then

        if( m1 == 0 ) then
          Trans(m1,m2) = (1._db,0._db)
        elseif( m1 > 0 ) then
          Trans(m1,m2) = is * r2_r
        else
          Trans(m1,m2) = - r2_i
        endif

      elseif( m1 == - m2 ) then

        if( m1 > 0 ) then
          Trans(m1,m2) = r2_r
        else
          Trans(m1,m2) = is * r2_i
        endif

      endif

    end do
  end do

  return
end

!***********************************************************************

! Calculation and writing of the atomic electron density

! Similar to the first part of Cal_state

subroutine Cal_Density(Energ,Full_atom,iaabsi,iaprotoi,icheck,ie,ie_computer,Int_statedens,itypei,itypepr,lamstdens, &
                lla_state,lla2_state,lmaxat,mpinodes,n_atom_0,n_atom_ind,n_atom_proto,natome,nenerg,nomfich_s,nonexc_g,nrato, &
                nrm,nspinp,ntype,numat,Rato,Rmtsd,State_all,Statedens)

  use declarations
  implicit none

  integer:: i_self, ia, iaabsi, iapr, icheck, ie, ie_computer, ipr, iprabs, ir, isp, it, l, &
    la, lamstdens, ll, lla_state, lla2_state, lm,lm1, lm2, lma, m, mpinodes, &
    n_atom_0, n_atom_ind, n_atom_proto, natome, nenerg, nr, nrm, nspinp, ntype, Z

  character(len=132) nomfich_s

  integer, dimension(0:ntype):: nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr, la_ipr, ll_ipr, lmaxat
  integer, dimension(natome):: iaprotoi, itypei

  logical:: Abs_exc, Full_atom, nonexc_g, Open_file, State_all
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
        ll = min(2,lmaxat(ipr))
      elseif( Z > 18 .and. Z < 55 ) then
        la = min(2,lmaxat(ipr))
        ll = min(3,lmaxat(ipr))
      else
        la = min(3,lmaxat(ipr))
        ll = min(4,lmaxat(ipr))
      endif
    end if
    la_ipr(ipr) = la    ! nombre maximal d'harmoniques pour l'ecriture de la densite d'etat
    ll_ipr(ipr) = ll    ! nombre maximal d'harmoniques pour le calcul  de la densite d'etat
  end do

  boucle_ia: do ia = 1,natome

    ipr = iaprotoi(ia)
    la = la_ipr(ipr)
    ll = ll_ipr(ipr)

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
     
! Calcul de la densite d'etat integree jusqu'au rayon Rmstd de l'atome

    do ir = 1,nrato(it)
      if( rato(ir,it) > Rmtsd(ipr) + eps10 ) exit
    end do
    nr = ir

    lma = (ll + 1)**2

    Statedens_l(:,:) = 0._db
    do l = 0,la
      lm1 = l**2 + 1
      lm2 = ( l + 1 )**2
      do isp = 1,nspinp
        do lm = lm1,lm2
          Statedens_l(l,isp) = Statedens_l(l,isp) + Statedens(lm,isp,lm,isp,iapr,ie_computer)
        end do
      end do
    end do
    if( nspinp == 1 ) Statedens_l(0:la,1) = 2 * Statedens_l(0:la,1)

! Integrale de la densite d'etats
    do isp = 1,nspinp
      do lm = 1,lma
        Int_Statedens(lm,isp,iapr) = Int_Statedens(lm,isp,iapr) + de * Statedens(lm,isp,lm,isp,iapr,ie_computer)
      end do
    end do

    do l = 0,ll
      lm1 = l**2 + 1
      lm2 = ( l + 1 )**2
      do isp = 1,nspinp
        Int_Statedens_l(l,isp,iapr) = sum(Int_Statedens(lm1:lm2,isp,iapr))
      end do
      if( nspinp == 1 ) Int_Statedens_l(l,1,iapr) = 2 * Int_Statedens_l(l,1,iapr)
    end do

    Int_statedens_t(:,iapr) = 0._db
    do isp = 1,nspinp
      Int_Statedens_t(isp,iapr) = Int_Statedens_t(isp,iapr) + sum( Int_Statedens_l(0:la,isp,iapr) )
    end do

    if( icheck > 1 ) then
      write(3,180) iapr
      lm = 0
      do l = 0,ll
        do m = -l,l
          lm = lm + 1
          do isp = 1,nspinp
            write(3,190) l, m, isp, Statedens(lm,isp,lm,isp,iapr,ie_computer), Int_Statedens(lm,isp,iapr)
          end do
        end do
      end do
      write(3,200)
      do l = 0,la
        write(3,210) l, Int_Statedens_l(l,1:nspinp,iapr)
      end do
      write(3,220) Int_Statedens_t(1:nspinp,iapr)
    endif

    Open_file = ie == 1
    Abs_exc = Full_atom .and. ( ia == iaabsi )  .and. .not. nonexc_g

  ! Writting file common with cal_state in SCF.f90
    call write_stdens(Abs_exc,.true.,Energ(ie),i_self,iapr,ie_computer,la,Int_Statedens,Int_Statedens_t, &
           Int_Statedens_l,lla_state,lla2_state,mpinodes,n_atom_0,n_atom_ind,nomfich_s,nspinp,Open_file,Statedens, &
           Statedens_l)

    proto_done(ipr) = .true.

  end do boucle_ia

  return
  110 format(/' ---- Cal_density --------',100('-'))
  120 format(15x,' Energy =',f10.3,' eV')
  180 format(/'  l  m is  Density of states   Integral    ia =',i3)
  190 format(3i3,3f15.7)
  200 format(/'    l   sum_m(Integral)')
  210 format(i5,2f15.7)
  220 format(' Total =',f12.7,f15.7)
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

    do l = 2,lmax
      l1 = 2*l - 1
      bessel(l) = l1 * bessel(l-1) / z - bessel(l-2)
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

  do l = 2,lmax
    l1 = 2*l - 1
    neuman(l) = l1 * neuman(l-1) / z - neuman(l-2)
    bessel(l) = l1 * bessel(l-1) / z - bessel(l-2)
  end do

  return
end

!***********************************************************************

! Calculation of spherical Hankel function of first or second kind ( or Riccati-Hankel(+) or (-) function, divided by z )
! and its derivative versus r, not z

subroutine Cal_Hankel(d_Hankel,Hankel,konde,lmax,r,Second_kind)

  use declarations
  implicit none

  integer:: l, lmax

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

  do l = 2,lmax+1
    hankel_t(l) = (2*l - 1) * hankel_t(l-1) / z - hankel_t(l-2)
  end do

  do l = 0,lmax
    d_hankel(l) = l * hankel_t(l) / z - hankel_t(l+1)
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

  integer:: l, lmax

  complex(kind=db):: cfac, fnorm, konde, z
  complex(kind=db), dimension(0:lmax):: d_Bessel, Bessel
  complex(kind=db), dimension(0:lmax+1):: Bessel_t

  real(kind=db):: r

  z = konde * r

  Bessel_t(0) = sin( z )

  Bessel_t(1) = - cos( z ) + Bessel_t(0) / z

  do l = 2,lmax+1
    Bessel_t(l) = (2*l - 1) * Bessel_t(l-1) / z - Bessel_t(l-2)
  end do

  do l = 0,lmax
    d_Bessel(l) = l * Bessel_t(l) / z - Bessel_t(l+1)
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

  do l = 2,lmax
    l1 = 2*l - 1
    neuman(l) = l1 * neuman(l-1) / z - neuman(l-2)
    bessel(l) = l1 * bessel(l-1) / z - bessel(l-2)
  end do

  return
end

!***********************************************************************

! Calculation of Bessel function ( Riccati-Bessel / z )
! When Bessel < 10^(-7), one takes the origin limit: j_l = z^l / (2*l+1)!!
! because of unstability of the recursive formula.

subroutine cbessel(bess,ip0,lmax,nr,q,r)

  use declarations
  implicit none

  integer:: i, ip0, j, k, l, l1, lmax, nr

  real(kind=db):: fac, q, z_lim
  real(kind=db), dimension(nr):: r, z
  real(kind=db), dimension(nr,0:lmax):: bessel
  real(kind=db), dimension(nr,ip0:lmax):: bess

  z(:) = q * r(:)

  do l = 0,lmax

    select case(l)

      case(0)

        do i = 1,nr
          if( z(i) < 1.e-7_db ) then
            bessel(i,l) = 1._db - z(i)**2 / 6._db
          else
            bessel(i,l) = sin( z(i) ) / z(i)
          endif
        end do

      case(1)

        fac = 1._db / 3
        z_lim = 3.e-7_db  ! = ( 1.e-7_db / fac )**( 1._db / l )
        do i = 1,nr
          if( z(i) < z_lim ) then
            bessel(i,l) = fac * z(i)
          else
            bessel(i,l) = ( bessel(i,0) - cos( z(i) ) ) / z(i)
          endif
        end do

      case default

        k = 1
        do j = 3,2*l+1,2
          k = k * j
        end do
        fac = 1._db / k
        z_lim = ( 1.e-7_db / fac )**( 1._db / l )
        l1 = 2*l - 1
 ! l = 2, fac = 1/15, z_lim = 0.00122
 ! l = 3, fac = 1/105, z_lim = 0.0219
 ! l = 4, fac = 1/945, z_lim = 0.0986
 ! l = 5, fac = 1/10395, z_lim = 0.253

        do i = 1,nr
          if( z(i) < z_lim ) then
            bessel(i,l) = fac * z(i)**l
          else
            bessel(i,l) = l1 * bessel(i,l-1) / z(i) - bessel(i,l-2)
          endif
        end do

    end select

  end do

  do l = ip0,lmax
    bess(:,l) = bessel(:,l)
  end do

  return
end


