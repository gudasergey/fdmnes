!FDMNES subroutines

subroutine Lm_cluster_cal(i,icheck,iopsymr,Lm_cluster,Lmax,n,nlm_pot_r)

  use declarations
  implicit none
  
  integer:: i, icheck, im, L, Lm, Lmax, m, m2, m4, mlm, mm, n, nlm_pot_r
  integer, dimension(nopsm):: iopsym, iopsymr
  integer, dimension(nlm_pot_r,2):: Lm_cluster
 
  iopsym(:) = abs( iopsymr(:) )
    
  Lm = 0
  do L = 0,Lmax
! Centro symmetry
    if( iopsym(25) == 1 .and. mod(L,2) == 1 ) cycle
    do m = - L,L
    
      im = abs( m )

! C2x      
      mm = mod( L + im, 2 )
      if( iopsym(22) == 1 .and. .not. ( (mm == 0 .and. m >= 0) .or. (mm == 1 .and. m < 0) ) ) cycle     
! C2y      
      mm = mod( L + 2*im, 2 )
      if( iopsym(23) == 1 .and. .not. ( (mm == 0 .and. m >= 0) .or. (mm == 1 .and. m < 0) ) ) cycle     
! C2z      
      mm = mod( im, 2 )
      if( iopsym(24) == 1 .and. mm /= 0 ) cycle     

! C4 along (110) or (1-10)
      m4 = mod( im, 4 )
      mlm = mod( L + im, 2 )
      if( ( iopsym(10) == 1 .or. iopsym(11) == 1 ) .and. .not. &      
             ( ( m >= 0 .and. ( ( m4 == 0 .and. mlm == 0 ) .or. ( m4 == 2 .and. mlm == 1 ) ) ) &
          .or. ( m < 0  .and. ( ( m4 == 2 .and. mlm == 0 ) .or. ( m4 == 0 .and. mlm == 1 ) ) ) ) ) cycle

! C3z
      mm = m
      do
        if( mm >= 0 ) exit
        mm = mm + 3
      end do
      mm = mod( mm, 3 )
      if( ( iopsym(49) == 1 .or. iopsym(50) == 1 ) .and. mm /= 0 ) cycle

! S3z
      if( ( iopsym(53) == 1 .or. iopsym(54) == 1 ) .and. .not. ( iopsym(53) == iopsym(49) .and. mod(L+m,2) == 0 &
                                                                                                  .and. mm == 0 ) ) cycle

! C4z
      mm = m
      do
        if( mm >= 0 ) exit
        mm = mm + 4
      end do
      mm = mod(mm,4)
      if( ( iopsym(18) == 1 .or. iopsym(21) == 1 ) .and. mm /= 0 ) cycle
      
! S4   
      if( ( iopsym(28) == 1 .or. iopsym(31) == 1 ) .and. .not. ( iopsym(28) == iopsym(18) .and. mod(L+m,2) == 0 &
                                                                                                  .and. mm == 0 ) ) cycle

! C6
      mm = m
      do
        if( mm >= 0 ) exit
        mm = mm + 6
      end do
      mm = mod(mm,6)
      if( ( iopsym(51) == 1 .or. iopsym(52) == 1 ) .and. mm /= 0 ) cycle

! S6   
      if( ( iopsym(55) == 1 .or. iopsym(56) == 1 ) .and. .not. ( iopsym(55) == iopsym(51) .and. mod(L+m,2) == 0 &
                                                                                                  .and. mm == 0 ) ) cycle
                                                                                                  
! mx
      m2 = mod(im,2)
      if( iopsym(40) == 1 .and. ( (m > 0 .and. m2 == 1) .or. (m < 0 .and. m2 == 0) ) ) cycle

! my
      if( iopsym(41) == 1 .and. m < 0 ) cycle

! mz
      if( iopsym(42) == 1 .and. mod(L+m,2) == 1 ) cycle
      
! Diagonal planes parallel to Oz
      m4 = mod(abs(m),4)
      if( ( iopsym(45) == 1 .or. iopsym(48) == 1 ) .and. .not. ( (m4 == 0 .and. m >= 0) .or. (m4 == 2 .and. m < 0) ) ) cycle

      Lm = Lm + 1
      
      if( i == 2 ) then
        Lm_cluster(Lm,1) = L
        Lm_cluster(Lm,2) = m
      endif
      
    end do
  end do

  n = Lm

  if( i == 2 .and. icheck > 2 ) then
    write(3,'(/A)') ' Spherical hamonics for the potential expansion in the cluster'
    write(3,'(/A)') '  Lm   L   m'
    do Lm = 1,nlm_pot_r
      write(3,'(3i4)') Lm, Lm_cluster(Lm,:)
    end do
  endif
    
  return
end

!***********************************************************************

! Elaboration of the radial grid for the cluster for the Numerov'e method

subroutine Cluster_radial_grid_Numerov(i_cycle,icheck,natome,nr,nr_cluster,posi,r_cluster,Rsort) 

  use declarations
  implicit none

  integer, parameter:: nb_dbl_cluster_c = 12
  integer, parameter:: nb_dbl_cluster = 7

  integer:: i_cycle, ia, ib, icheck, id, idb, ir, jr, n, natome, nd, nr, nr_cluster  

  real(kind=db):: dr, dr_cluster, h, P, R, r_min, Ray_dbl_cluster, Rsort
  
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(nr_cluster):: r_cluster
  real(kind=db), dimension(0:nb_dbl_cluster_c):: dr_dbl, Ray_dbl
  real(kind=db), dimension(natome):: dist

! Evaluation of the radius where interpoint distances is multiplied by 2, for Numerov's method

  Ray_dbl_cluster = 2._db / bohr
  dr_cluster =  0.08_db / bohr   ! 0.04_db / bohr
  r_min = 0.0002_db / bohr 

  if( icheck > 1 .and. i_cycle == 2 ) then
    write(3,110)
    write(3 ,'(/A)') ' id      Radius         dr'
  endif

  P = 1 + dr_cluster / Ray_dbl_cluster

  do id = 0,nb_dbl_cluster_c

    R =  id * log(2._db) / log(P)
    Ray_dbl(id) = Ray_dbl_cluster / P**R
    dr_dbl(id) = dr_cluster / ( 2**id )   
    if( icheck > 1 .and. i_cycle == 2 ) write(3,'(i3,2f15.11)') id, Ray_dbl(id)*bohr, dr_dbl(id)*bohr
     
  end do

  do ia = 1,natome
    dist(ia) = sqrt( sum( posi(:,ia)**2 ) )
  end do
  
  do idb = 0,nb_dbl_cluster-1
    if( Rsort - dr_cluster - dist(natome) > Ray_dbl(idb) ) exit
  end do
  if( idb == 0 ) then
    h = dr_cluster
  else
    h = dr_cluster / 2**idb
  endif
  
  R = Rsort + 2 * h

  if( i_cycle == 1 ) then
    ir = 0
  else
    ir = nr_cluster + 1
  endif
  
  Loop_ia: do ia = natome,1,-1     

    if( ia < natome ) then
      if( dist(ia+1) - dist(ia) < dr_cluster / 2**nb_dbl_cluster ) cycle
    endif
    
    if( ia == 1 ) then
      nd = nb_dbl_cluster_c
    else
      nd = nb_dbl_cluster
    endif
    do id = idb,nd
      n = nint( ( R - dist(ia) - Ray_dbl(id) ) / h )
      do jr = 1,n   
        ir = ir - 1
        R = R - h
        if( i_cycle == 2 ) r_cluster(ir) = R
      end do
      h = h / 2  
    end do
    if( ia == 1 ) then
      n = nint( ( R - r_min ) / h )
    else
      n = nint( ( R - dist(ia) ) / h )
    endif
    do jr = 1,n   
      ir = ir - 1
      R = R - h
      if( i_cycle == 2 ) r_cluster(ir) = R
    end do
    if( ia == 1 ) exit

    do ib = ia-1,1,-1
      if( dist(ib) + h < dist(ia) ) exit
    end do    
    P = ( dist(ia) +  dist(ib) ) / 2
    
    do id = nd,0,-1
      do
        ir = ir - 1
        R = R - h
        if( i_cycle == 2 ) r_cluster(ir) = R
        if( R < P ) cycle Loop_ia
        if( R < dist(ia) - Ray_dbl(id) ) exit
      end do
      h = 2 * h
    end do
    idb = max( id, 0 )  
  
  end do Loop_ia
  
  if( i_cycle == 1 ) nr = abs( ir )
      
  if( icheck > 2 .and. i_cycle == 2 ) then
    write(3,'(/A)') ' index    radius       delta'
    do ir = 1,nr_cluster
      if( ir == 1 ) then
        dr = r_cluster(ir)
      else
        dr = r_cluster(ir) - r_cluster(ir-1)
      endif 
      write(3,'(i6,1p,2e13.5)') ir, r_cluster(ir) * bohr, dr * bohr 
    end do
  endif 
    
  return
  110 format(/' ----------- Cluster radial grid Numerov ',80('-'))
end

!***********************************************************************

! Elaboration of the radial grid for the cluster

subroutine Cluster_radial_grid(i_cycle,icheck,itypei,natome,nr,nr_cluster,nrato,nrm,ntype,posi,r_cluster,rato,Rsort) 

  use declarations
  implicit none

  integer:: i_cycle, j_cycle, ia, ib, icheck, ir, ir_cluster, ir_min, it, jr, kr, natome, nr, nr_cluster, nr_max, nrm, ntype  

  integer, dimension(natome):: itypei
  integer, dimension(0:ntype):: nrato

  real(kind=db):: dr, dr_min, P, R, Rsort
  
  real(kind=db), dimension(0:nrm,0:ntype):: Rato
  real(kind=db), dimension(3,natome):: posi
  real(kind=db), dimension(nr_cluster):: r_cluster
  real(kind=db), dimension(natome):: dist

  real(kind=db), dimension(:,:), allocatable:: Ra

! See also in Pot_SphereHarm
  dr_min = 0.00005_db / bohr
  
! Definition of atomic radial grid

  if( natome > 1 ) then
  
    do j_cycle = 1,2
    
      if( j_cycle == 2 ) then
        allocate( Ra(0:nr_max,0:ntype) )
        Ra(:,:) = 0._db
      endif
      
      Loop_it: do it = 0,ntype
      
        do ir_min = 2,nrato(it)
          if( rato(ir_min,it) - rato(ir_min-1,it) > dr_min ) exit
        end do

        do ir = 0,nrato(it)
          R = ir * dr_min
          if( j_cycle == 2 ) Ra(ir,it) = R 
          if( R > Rsort ) cycle Loop_it
          if( R > rato(ir_min,it) ) exit
        end do
        R = rato(ir_min,it)  
        if( j_cycle == 2 ) Ra(ir,it) = R 
        do jr = ir_min + 1, nrato(it)
          ir = ir + 1 
          R = rato(jr,it) 
          if( j_cycle == 2 ) then
            Ra(ir,it) = R
            if( jr == ir_min + 3 .and. ir > 7 ) then
              dr = ( R - Ra(ir-6,it) ) / 6
              do kr = ir-5,ir-1
                Ra(kr,it) = Ra(kr-1,it) + dr 
              end do
            endif
          endif 
          if( R > Rsort ) cycle Loop_it
        end do
        
        dr = rato(nrato(it),it) - rato(nrato(it)-1,it)
        do
          ir = ir + 1
          R = R + dr
          if( R > Rsort ) exit 
          if( j_cycle == 2 ) Ra(ir,it) = R 
        end do 
      end do Loop_it
      
      nr_max = max( nr_max, ir )
      
    end do
  
  else
  
    nr_max = nrato( itypei(1) )
     
  endif

  do ia = 1,natome
    dist(ia) = sqrt( sum( posi(:,ia)**2 ) )
  end do

  ir_cluster = 0
  R = 0._db
  do ia = 1,natome
    it = itypei(ia)

    if( natome == 1 .or. ia == natome ) then

      do ir = 1,nr_max
        ir_cluster = ir_cluster + 1
        if( ia == 1 ) then
          R = rato(ir,it)
        else
          R = R + Ra(ir,it) - Ra(ir-1,it)
        endif
        if( i_cycle == 2 ) r_cluster(ir_cluster) = R 
        if( R < Rsort - eps10 ) cycle
        if( abs( R - Rsort ) < eps10 ) exit
        R = Rsort
        if( i_cycle == 2 ) then
          r_cluster(ir_cluster) = R
          if( ir_cluster > 7 ) then
            dr = ( R - r_cluster(ir_cluster-6) ) / 6
            do kr = ir_cluster-5,ir_cluster-1
              r_cluster(kr) = r_cluster(kr-1) + dr 
            end do
          endif
        endif
        exit
      end do
      ir_cluster = ir_cluster + 1
      if( i_cycle == 2 ) r_cluster(ir_cluster) = 2 * r_cluster(ir_cluster-1) - r_cluster(ir_cluster-2)
      
    else

      if( ia > 1 ) then    
        if( dist(ia) == dist(ia-1) ) cycle
      endif 
 
      do ib = ia + 1, natome
        if( abs( dist(ib) - dist(ia) ) > eps10 ) exit
      end do
      if( ib > natome ) cycle
      P = dist(ia) + ( dist(ib) - dist(ia) ) / 2

      do ir = 1,max( nr_max, nrato(it) )
        ir_cluster = ir_cluster + 1
        if( ia == 1 ) then
          R = rato(ir,it)
        else
          R = R + Ra(ir,it) - Rato(ir-1,i_cycle)
        endif
        if( i_cycle == 2 ) r_cluster(ir_cluster) = R 
        if( R > P ) exit
      end do
      
      it = itypei(ib)
      do ir = 1,nr_max
        if( Ra(ir,it) > dist(ib) - R - eps10 ) exit
      end do
      jr = ir - 1
      do ir = jr,0,-1
        ir_cluster = ir_cluster + 1
        R = dist(ib) - Ra(ir,it)
        if( i_cycle == 2 ) then
          r_cluster(ir_cluster) = R
          if( ir == jr-2 .and. ( dist(ib) - dist(ia) > R - r_cluster(ir_cluster-6) )  ) then
            dr = ( R - r_cluster(ir_cluster-6) ) / 6
            do kr = ir_cluster-5,ir_cluster-1
              r_cluster(kr) = r_cluster(kr-1) + dr 
            end do
          endif
        endif 
      end do
     
    endif
  end do

  nr = ir_cluster
    
  if( icheck > 2 .and. i_cycle == 2 ) then
    write(3,110)
    write(3,'(/A)') ' index    radius       delta'
    do ir = 1,nr_cluster
      if( ir == 1 ) then
        dr = r_cluster(ir)
      else
        dr = r_cluster(ir) - r_cluster(ir-1)
      endif 
      write(3,'(i6,1p,2e13.5)') ir, r_cluster(ir) * bohr, dr * bohr 
    end do
  endif 
  
  if( natome > 1 ) deallocate( Ra )
    
  return
  110 format(/' ----------- Cluster radial grid ',80('-'))
end

!***********************************************************************

subroutine Core_Psi_radial_grid(i_cycle,icheck,nbseuil,nr,nr_cluster,nr_core,nrato_abs,nrm,Psi_core,Psii,r_cluster,r_core, &
                            Rato_abs)

  use declarations
  implicit none

  integer:: i_cycle, icheck, ir, is, jr, nbseuil, nr, nr_cluster, nr_core, nr_max, nrato_abs, nrm

  real(kind=db):: f, f1, f2, f3, f4, fp, r1, r2, r3, r4, Rmax 
  real(kind=db), dimension(nrato_abs):: Rato_abs 
  real(kind=db), dimension(nr_cluster):: r_cluster 
  real(kind=db), dimension(nr_core):: r_core 
  real(kind=db), dimension(nrm,nbseuil):: Psii
  real(kind=db), dimension(nr_core,nbseuil):: Psi_core

  do ir = 1,nrm
    if( abs( Psii(ir,1) ) < eps15 ) exit 
  end do
  nr_max = ir

  Rmax = Rato_abs(nr_max-1)

  nr = 1  
  do
    if( r_cluster(nr) > Rmax + eps10 ) exit
    if( nr == nr_cluster ) then
      Rmax = r_cluster(nr)
      exit 
    endif 
    nr = nr + 1
  end do
  
  if( i_cycle == 1 ) return
  
  r_core(1:nr_core) = r_cluster(1:nr_core)

  do ir = 1,nr_core
    do jr = 3,nr_max-1
      if( Rato_abs(jr) > r_core(ir) ) exit
    end do
    r1 = Rato_abs(jr-2)  
    r2 = Rato_abs(jr-1) 
    r3 = Rato_abs(jr) 
    r4 = Rato_abs(jr+1)
    do is = 1,nbseuil 
      f1 = Psii(jr-2,is)  
      f2 = Psii(jr-1,is) 
      f3 = Psii(jr,is) 
      f4 = Psii(jr+1,is)
      call interp3(f,fp,r_core(ir),r1,r2,r3,r4,f1,f2,f3,f4)
      Psi_core(ir,is) = f 
    end do
  end do  
  
  if( icheck > 2 ) then
    write(3,110)
    if( nbseuil == 1 ) then      
      write(3,'(/A)') '    Radius     Psi_core'
    else
      write(3,'(/A)') '    Radius     Psi_core(1)    Psi_core(2)'
    endif
    do ir = 1,nr_core
      write(3,'(f11.7,1p,2e13.5)') r_core(ir) * bohr, Psi_core(ir,:) 
    end do
  endif 

  return
  110 format(/' ----------- Core Psi radial grid ',80('-'))
end

!***********************************************************************

subroutine PotSup_SphereHarm(alfpot,E_max,iaproto,icheck,igrpt_nomag,iopsymr,itypep,itypepr,Lm_cluster, &
            Lmax_cluster,mpirank0,n_atom_proto,natomp,nlat,nlatm,nlm_pot_cp,nlm_pot_g,nlm_pot_r,nr_cluster,nrato,nrm, &
            nspin,ntype,numat,popatm,Popatv,Pos,Psival,r_cluster,Rho_cluster,Rho_cluster_cp,Rato,Rhoit,Vc_cluster,Vc_cluster_cp, &
            Vxc_cluster,Vxc_cluster_cp,Ylm_comp)

  use declarations
  implicit none

  integer:: igrpt_nomag, it, Lmax_cluster, mpirank0, n_atom_proto, natomp, nlatm, nlm_pot_cp, nlm_pot_g, nlm_pot_r, &
            nr_cluster, nrm, nspin, ntype

  integer, dimension(30):: icheck
  integer, dimension(0:ntype):: nlat, nrato, numat
  integer, dimension(0:n_atom_proto):: itypepr
  integer, dimension(natomp):: iaproto, itypep
  integer, dimension(nopsm):: iopsymr
  integer, dimension(nlm_pot_g,2):: Lm_cluster

  complex(kind=db), dimension(nr_cluster,nlm_pot_cp):: Vc_cluster_cp
  complex(kind=db), dimension(nr_cluster,nspin,nlm_pot_cp):: Rho_cluster_cp, Vxc_cluster_cp
  
  logical:: Ylm_comp

  real(kind=db):: alfpot, E_max
  real(kind=db), dimension(3,natomp):: Pos
  real(kind=db), dimension(nr_cluster):: r_cluster
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm
  real(kind=db), dimension(0:nrm,0:ntype):: Rato, Rhoit
  real(kind=db), dimension(0:nrm,0:n_atom_proto):: Vato
  real(kind=db), dimension(0:nrm,nspin,0:n_atom_proto):: Rhoigr
  real(kind=db), dimension(0:ntype,nlatm):: Popatv
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: Psival
  real(kind=db), dimension(nr_cluster,nlm_pot_r):: Vc_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_pot_r):: Rho_cluster, Vxc_cluster

  do it = 0,ntype
    call potato(icheck(10),it,itypepr,n_atom_proto,nlat,nlatm,nrato,nrm,nspin,ntype,numat, &
                popatm,popatv,psival,Rato,Rhoigr,Rhoit,Vato)
  end do

  call Pot_SphereHarm(alfpot,E_max,iaproto,icheck(13),igrpt_nomag,iopsymr,itypep,Lm_cluster,Lmax_cluster,mpirank0,n_atom_proto, &
               natomp,nlm_pot_cp,nlm_pot_g,nlm_pot_r,nr_cluster,nrato,nrm,nspin,ntype,numat,Pos,r_cluster,Rato,Rhoigr, &
               Rho_cluster,Rho_cluster_cp,Vato,Vc_cluster,Vc_cluster_cp,Vxc_cluster,Vxc_cluster_cp,Ylm_comp)


  return
end

!***********************************************************************

subroutine Pot_SphereHarm(alfpot,E_max,iaproto,icheck,igrpt_nomag,iopsymr,itypep,Lm_cluster,Lmax_cluster,mpirank0,n_atom_proto, &
               natomp,nlm_pot_cp,nlm_pot_g,nlm_pot_r,nr_cluster,nrato,nrm,nspin,ntype,numat,Pos,r_cluster,Rato,Rhoigr, &
               Rho_cluster,Rho_cluster_cp,Vato,Vc_cluster,Vc_cluster_cp,Vxc_cluster,Vxc_cluster_cp,Ylm_comp)

  use declarations
  implicit none

  integer:: i, ia, ib, icheck, icheck_Dlmm, igrpt_nomag, ipr, ir, ir_a, ir_a_p, ir_a_p_min, isp, isym, it, L, Length, Lm, Lmax, &
           Lmax_cluster, m, mpirank0, mult, multiplicity, n_atom_proto, natomp, nlm_pot_cp, nlm_pot_g, nlm_pot_r, &
           nr_cluster, nrm, nspin, ntype
  
  integer, dimension(32):: Symmetry
  integer, dimension(0:ntype):: nrato, numat

  integer, dimension(nlm_pot_g,2):: Lm_cluster
  integer, dimension(nopsm):: iopsymr
  integer, dimension(natomp):: iaproto, itypep
  
  character(len=10):: mot10
  character(len=10), dimension(:), allocatable:: Col_name, Col_name_sp

  complex(kind=db), dimension(nlm_pot_g):: Dlm0
  complex(kind=db), dimension(nr_cluster,nlm_pot_cp):: Vc_cluster_cp
  complex(kind=db), dimension(nr_cluster,nspin,nlm_pot_cp):: Rho_cluster_cp, Vxc_cluster_cp
  complex(kind=db), dimension(:,:,:), allocatable:: Dlmm
  
  logical:: Limitation, Outside_atom, Ylm_comp

  real(kind=db):: alfpot, dist, dr, dr_min, dTheta, Int_sin_theta_dtheta, E_max, Fac, p, Pot_lim, r_a, r_c, sin_theta_dtheta, &
                  Theta, Vc_a, Vc_out  
  real(kind=db), dimension(3):: Ang, dp, ps, V
  real(kind=db), dimension(nspin):: Rh, Rho_a, Vxc_a, Vxc_out
  real(kind=db), dimension(3,natomp):: Pos
  real(kind=db), dimension(3,3):: Matopsym, Rot
  real(kind=db), dimension(0:nrm,0:ntype):: Rato
  real(kind=db), dimension(nr_cluster,nlm_pot_r):: Vc_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_pot_r):: Rho_cluster, Vxc_cluster
  real(kind=db), dimension(nr_cluster):: r_cluster
  real(kind=db), dimension(0:nrm,0:n_atom_proto):: Vato
  real(kind=db), dimension(0:nrm,nspin,0:n_atom_proto):: Rhoigr

  real(kind=db), dimension(:), allocatable:: Vc, Ylm, YlmSintdomega
  real(kind=db), dimension(:,:), allocatable:: Rho

  if( icheck > 0 ) write(3,'(/A)') ' ------- Pot_SphereHarm ----------------------------------------------------------------------'
  
  icheck_Dlmm = icheck - 1
  
  dr_min = 0.00005_db / bohr
  
! Minimum potential for atoms just outside the outersphere
  Pot_lim = - 1._db
  
  if( Ylm_comp ) then
    Vc_cluster_cp(:,:) = ( 0._db, 0._db )
    Vxc_cluster_cp(:,:,:) = ( 0._db, 0._db )
    Rho_cluster_cp(:,:,:) = ( 0._db, 0._db )
  else
    Vc_cluster(:,:) = 0._db
    Vxc_cluster(:,:,:) = 0._db
    Rho_cluster(:,:,:) = 0._db
  endif

  allocate( Dlmm(-Lmax_cluster:Lmax_cluster,-Lmax_cluster:Lmax_cluster,0:Lmax_cluster) )
    
  do ia = 1,natomp

    ps(:) = pos(:,ia)
    call posequiv(mpirank0,ps,iopsymr,isym,igrpt_nomag)
    dp(:) = abs( ps(:) - pos(:,ia) )
    if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos ) cycle

    it = itypep(ia)
    ipr = iaproto( ia )

    dist = sqrt( sum( ps(:)**2 ) )

    Outside_atom = dist > r_cluster(nr_cluster)
    
! Looking for equivalent atoms
    mult = 1
    Symmetry(mult) = 1
    do ib = 1,natomp
      if( ib == ia .or. iaproto( ib ) /= ipr ) cycle
      V(:) = pos(:,ib)
      call posequiv(mpirank0,V,iopsymr,isym,igrpt_nomag)
      dp(:) = abs( V(:) - ps(:) )
      if( dp(1) > epspos .or. dp(2) > epspos .or. dp(3) > epspos ) cycle
      
      mult = mult + 1
      Symmetry(mult) = isym
    end do
    multiplicity = mult
    
    if( dist < eps10 ) then

      Fac = sqrt( quatre_pi )
      
! Only the isotropical term for the central atom
      Lm = 1       
      do ir = 1,nr_cluster-1
        r_c = r_cluster(ir)
        do ir_a = 2,nrato(it)
          if( rato(ir_a,it) > r_c ) exit
        end do
        if( ir_a > nrato(it) ) exit
        p = ( rato(ir_a,it) - r_c ) / ( rato(ir_a,it) - rato(ir_a-1,it) ) 
        Vc_a = p * Vato(ir_a-1,ipr) + ( 1 - p ) * Vato(ir_a,ipr)
        Rho_a(:) = p * Rhoigr(ir_a-1,:,ipr) + ( 1 - p ) * Rhoigr(ir_a,:,ipr)
        if( Ylm_comp ) then
          Vc_cluster_cp(ir,Lm) = Vc_cluster(ir,Lm) + Fac * cmplx( Vc_a, 0._db, db )
          Rho_cluster_cp(ir,:,Lm) = Rho_cluster_cp(ir,:,Lm) + Fac * cmplx( Rho_a(:), 0._db, db )
        else
          Vc_cluster(ir,Lm) = Vc_cluster(ir,Lm) + Fac * Vc_a
          Rho_cluster(ir,:,Lm) = Rho_cluster(ir,:,Lm) + Fac * Rho_a(:)
        endif 
      end do
      
      cycle
      
    endif

    do ir_a_p_min = 2,nrato(it)-1
      if( rato(ir_a_p_min,it) - rato(ir_a_p_min-1,it) > dr_min ) exit
    end do

    Dlm0(:) = ( 0._db, 0._db )
    do mult = 1,multiplicity
      V(:) = ps(:) / dist
      if( mult > 1 ) then    
        call opsym(Symmetry(mult),Matopsym)
        V = matmul( Matopsym, V )
      endif
      Ang(1) = atan2( V(2), V(1) )
      Ang(2) = - acos( V(3) )
      Ang(3) = 0._db 
      call Mat_Euler(Ang,Rot)
      Rot = Transpose( Rot )
      call Rotation_mat(Dlmm,icheck_Dlmm,Rot,Lmax_cluster,Ylm_comp)
      do Lm = 1,nlm_pot_g
        L = Lm_cluster(Lm,1)
        if( L > Lmax_cluster ) exit
        m = Lm_cluster(Lm,2)
        Dlm0(Lm) = Dlm0(Lm) + Dlmm(0,m,L)
      end do
    end do

    Loop_ir: do ir = 1,nr_cluster-1

      r_c = r_cluster(ir)

      if( abs( r_c - dist ) < eps10 ) then
        if( Ylm_comp ) then
          Vc_cluster_cp(ir,:) = 2 * Vc_cluster_cp(ir-1,:) - Vc_cluster_cp(ir-2,:) 
          Rho_cluster_cp(ir,:,:) = 2 * Rho_cluster_cp(ir-1,:,:) - Rho_cluster_cp(ir-2,:,:)
        else
          Vc_cluster(ir,:) = 2 * Vc_cluster(ir-1,:) - Vc_cluster(ir-2,:)
          Rho_cluster(ir,:,:) = 2 * Rho_cluster(ir-1,:,:) - Rho_cluster(ir-2,:,:)
        endif
        cycle
      endif 

      call clmax(E_max,r_c,-1,Lmax,-1,.true.)
      Lmax = max( Lmax, 8 )
      Lmax = min( Lmax, Lmax_cluster )
      
      allocate( YlmSintdomega(0:Lmax) )
      allocate( Vc(0:Lmax) )
      allocate( Rho(0:Lmax,nspin) )
      
      r_a = abs( dist - r_c )
      do ir_a = 1,nrato(it) 
        if( rato(ir_a,it) > r_a ) exit
      end do
      
      ir_a = min( ir_a, nrato(it) )
      ir_a_p = ir_a
      if( ir_a == 1 ) then
        dr = dr_min
      else  
        dr = max( rato(ir_a,it) - rato(ir_a-1,it), dr_min )
      endif 
      dtheta = min( dr / r_c, pi / 24 ) 
      
      Theta = 0._db
      Int_sin_theta_dtheta = 0._db
      Vc(:) = 0._db
      Rho(:,:) = 0._db

      do
      
        if( r_a < rato(nrato(it),it) ) then

          Limitation = .false.
          do ir_a = nrato(it),1,-1
            if( rato(ir_a,it) < r_a ) exit
            if( Outside_atom .and. Vato(ir_a,ipr) < Pot_lim ) then
              Limitation = .true.
              exit
            endif
          end do
          if( Limitation ) then
            Vc_a = Pot_lim
            p = ( Vato(ir_a+1,ipr) - Pot_lim ) / ( Vato(ir_a+1,ipr) - Vato(ir_a,ipr) )   
          elseif( ir_a == 0 ) then
            ir_a_p = 0
            p = ( rato(1,it) - r_a ) / rato(1,it)  
            if( r_a > eps10 ) then
         ! potential in - 2*Z/r + shift
              Vc_a = Vato(1,ipr) + 2 * Numat(it) / rato(1,it) - 2 * Numat(it) / r_a
            else
            ! in order the average in the disk is '- 2*Z'
            ! in fact not used because the value at this point is extrapolated 
              Vc_a = - 2 * Numat(it) / ( 1 - cos( dtheta / 2 ) )    
            endif
          else
            p = ( rato(ir_a+1,it) - r_a ) / ( rato(ir_a+1,it) - rato(ir_a,it) ) 
            Vc_a = p * Vato(ir_a,ipr) + ( 1 - p ) * Vato(ir_a+1,ipr)
          endif
          ir_a = max( ir_a, 1 )
          Rho_a(:) = p * Rhoigr(ir_a,:,ipr) + ( 1 - p ) * Rhoigr(ir_a+1,:,ipr) 
        
        else
        
          Vc_a = 0._db
          Rho_a(:) = 0._db  
        
        endif
! Only m = 0 calculated, thus spherical harmonics always real          
        call Cal_YlmSintdomega(dtheta,Lmax,sin_theta_dtheta,Theta,YlmSintdomega)
        
        Int_sin_theta_dtheta = Int_sin_theta_dtheta + sin_theta_dtheta
        Vc(:) = Vc(:) + Vc_a * YlmSintdomega(:)
        do isp = 1,nspin
          Rho(:,isp) = Rho(:,isp) + Rho_a(isp) * YlmSintdomega(:)
        end do
        
        if( Theta > pi - eps10 ) exit
        
        if( ir_a_p == 0 ) then
          dr = dr_min
          ir_a_p = ir_a_p_min
        else
          ir_a_p = min( ir_a_p + 1, nrato(it) )
          dr = max( rato(ir_a_p,it) - rato(ir_a_p-1,it), dr_min )
        endif 
        dtheta = dr / r_c 
        dtheta = min( dr / r_c, pi / 24 ) 

        Theta = Theta + dTheta
        if( Theta > pi + eps10 ) then
          Theta = Theta - dtheta
          dtheta = 2 * ( pi - Theta )
          Theta = pi
          r_a = dist + r_c
        else
          r_a = sqrt( dist**2 + r_c**2 - 2 * dist * r_c * cos( Theta ) )
        endif

      end do
      
      do Lm = 1,nlm_pot_g
        L = Lm_cluster(Lm,1)
        if( L > Lmax ) exit
        if( Ylm_comp ) then
          Vc_cluster_cp(ir,Lm) = Vc_cluster_cp(ir,Lm) + Dlm0(Lm) * Vc(L)
          Rho_cluster_cp(ir,:,Lm) = Rho_cluster_cp(ir,:,Lm) + Dlm0(Lm) * Rho(L,:)
        else
          Vc_cluster(ir,Lm) = Vc_cluster(ir,Lm) + real( Dlm0(Lm), db ) * Vc(L)
          Rho_cluster(ir,:,Lm) = Rho_cluster(ir,:,Lm) + real( Dlm0(Lm), db ) * Rho(L,:)
        endif 
      end do
      
      deallocate( Rho, Vc, YlmSintdomega )
      
    end do Loop_ir
  
  end do
  
  deallocate( Dlmm )

  if( Ylm_comp ) then
    Rho_cluster_cp(nr_cluster,:,1) = Rho_cluster_cp(nr_cluster-1,:,1)
    if( nlm_pot_g > 1 ) Vc_cluster_cp(nr_cluster-1,2:nlm_pot_g) = ( 0._db, 0._db )  
    Vc_cluster_cp(nr_cluster,1) = Vc_cluster_cp(nr_cluster-1,1)
  else
    Rho_cluster(nr_cluster,:,1) = Rho_cluster(nr_cluster-1,:,1)
    if( nlm_pot_g > 1 ) Vc_cluster(nr_cluster-1,2:nlm_pot_g) = 0._db  
    Vc_cluster(nr_cluster,1) = Vc_cluster(nr_cluster-1,1)
  endif
  
! Calculation of exchange-correlation potential
  call Exchange_correlation(alfpot,E_max,Lm_cluster,Lmax_cluster,nlm_pot_cp,nlm_pot_g,nlm_pot_r,nr_cluster,nspin, &
              r_cluster,Rho_cluster,Rho_cluster_cp,Vxc_cluster,Vxc_cluster_cp,Ylm_comp)

  if( icheck > 0 ) then
    if( Ylm_comp ) then
      Vc_out = real( Vc_cluster_cp(nr_cluster,1), db ) * rydb
      Vxc_out(:) = real( Vxc_cluster_cp(nr_cluster,:,1), db ) * rydb
    else
      Vc_out = real( Vc_cluster(nr_cluster,1), db ) * rydb
      Vxc_out(:) = real( Vxc_cluster(nr_cluster,:,1), db ) * rydb
    endif
    if( nspin == 1) then  
      write(3,100) Vc_out + Vxc_out(:)   
      write(3,102) Vc_out, Vxc_out(:)
    else
      write(3,105) Vc_out + Vxc_out(:)   
      write(3,107) Vc_out, Vxc_out(:)
    endif   
  endif
                  
  if( icheck > 1 ) then

    allocate( Col_name(nlm_pot_g) )
    do Lm = 1,nlm_pot_g
      mot10 = '(' 
      L = Lm_cluster(Lm,1)
      m = Lm_cluster(Lm,2)
      call ad_number(L,mot10,10)
      Length = len_trim(mot10) + 1
      mot10(Length:Length) = ','
      call ad_number(m,mot10,10)
      Length = len_trim(mot10) + 1
      if( Length <= 10 ) mot10(Length:Length) = ')'
      Col_name(Lm) = mot10
    end do

    if( nspin == 2 ) then
      allocate( Col_name_sp(nspin*nlm_pot_g) )
      i = 0
      do Lm = 1,nlm_pot_g
        do isp = 1,nspin
          i = i + 1 
          mot10 = '(' 
          L = Lm_cluster(Lm,1)
          m = Lm_cluster(Lm,2)
          call ad_number(L,mot10,10)
          Length = len_trim(mot10) + 1
          mot10(Length:Length) = ','
          call ad_number(m,mot10,10)
          Length = len_trim(mot10) + 1
          if( Length <= 9 ) mot10(Length:Length) = ')'
          Length = len_trim(mot10) + 1
          if( Length <= 10 ) then 
            if( isp == 1 ) then
              mot10(Length:Length) = 'u'
            else
              mot10(Length:Length) = 'u'
            endif
          endif
          Col_name_sp(i) = mot10
        end do
      end do
    endif
    
    do ir = 1,nr_cluster
      if( Ylm_comp ) then
        if( ir == 1 ) write(3,110) ( ' Vc', Col_name(Lm), Lm = 1,nlm_pot_g )
        write(3,120) r_cluster(ir)*bohr, ( Vc_cluster_cp(ir,Lm)*Rydb, Lm = 1,nlm_pot_g )
      else
        if( ir == 1 ) write(3,130) ( ' Vc', Col_name(Lm), Lm = 1,nlm_pot_g )
        write(3,140) r_cluster(ir)*bohr, ( Vc_cluster(ir,Lm)*Rydb, Lm = 1,nlm_pot_g )
      endif 
    end do

    do ir = 1,nr_cluster
      if( Ylm_comp ) then
        if( nspin == 1 ) then
          if( ir == 1 ) write(3,110) ( 'Rho', Col_name(Lm), Lm = 1,nlm_pot_g )
          write(3,120) r_cluster(ir)*bohr, ( Rho_cluster_cp(ir,:,Lm), Lm = 1,nlm_pot_g )
        else
          if( ir == 1 ) write(3,110) ( 'Rho', Col_name_sp(Lm), Lm = 1,nspin*nlm_pot_g )
          write(3,120) r_cluster(ir)*bohr, ( Rho_cluster_cp(ir,:,Lm), Lm = 1,nlm_pot_g )
        endif
      else
        if( nspin == 1 ) then
          if( ir == 1 ) write(3,130) ( 'Rho', Col_name(Lm), Lm = 1,nlm_pot_g )
          write(3,140) r_cluster(ir)*bohr, ( Rho_cluster(ir,:,Lm), Lm = 1,nlm_pot_g )
        else
          if( ir == 1 ) write(3,130) ( 'Rho', Col_name_sp(Lm), Lm = 1,nspin*nlm_pot_g )
          write(3,140) r_cluster(ir)*bohr, ( Rho_cluster(ir,:,Lm), Lm = 1,nlm_pot_g )
        endif 
      endif 
    end do
    
    do ir = 1,nr_cluster
      if( Ylm_comp ) then
        if( nspin == 1 ) then
          if( ir == 1 ) write(3,110) ( 'Vxc', Col_name(Lm), Lm = 1,nlm_pot_g )
          write(3,120) r_cluster(ir)*bohr, ( Vxc_cluster_cp(ir,:,Lm)*rydb, Lm = 1,nlm_pot_g )
        else
          if( ir == 1 ) write(3,110) ( 'Vxc', Col_name_sp(Lm), Lm = 1,nspin*nlm_pot_g )
          write(3,120) r_cluster(ir)*bohr, ( Vxc_cluster_cp(ir,:,Lm)*rydb, Lm = 1,nlm_pot_g )
        endif
      else
        if( nspin == 1 ) then
          if( ir == 1 ) write(3,130) ( 'Vxc', Col_name(Lm), Lm = 1,nlm_pot_g ) 
          write(3,140) r_cluster(ir)*bohr, ( Vxc_cluster(ir,:,Lm)*rydb, Lm = 1,nlm_pot_g )
        else
          if( ir == 1 ) write(3,130) ( 'Vxc', Col_name_sp(Lm), Lm = 1,nspin*nlm_pot_g )
          write(3,140) r_cluster(ir)*bohr, ( Vxc_cluster(ir,:,Lm)*rydb, Lm = 1,nlm_pot_g )
        endif 
      endif 
    end do
    
    deallocate( Col_name )
    if( nspin == 2 ) deallocate( Col_name_sp ) 
    
    Allocate( Ylm(0:Lmax_cluster) )
    Ylm(0) = sqrt( 1 / quatre_pi )
    if( Lmax_cluster > 0 ) Ylm(1) = sqrt( 3 / quatre_pi )
    do L = 2,Lmax_cluster
      Ylm(L) = ( sqrt( 4*L**2 - 1._db ) / L ) * ( Ylm(L-1) - ( ( L - 1._db ) / sqrt((2*L-1._db)*(2*L-3._db)) ) * Ylm(L-2) )
    end do

    write(3,'(/A)') ' Potential along Oz'
    
    if( nspin == 1 ) then
      write(3,'(/A)') '    Radius       V             Vc           Vxc         Rho      4.pi.r2.Rho'
    else
      write(3,'(/A,A)') '    Radius       V_up          V_dn          Vc         Vxc_up      Vxc_dn         Rho_up      Rho_dn', &
                      '   4.pi.r2.Rho_u 4.pi.r2.Rho_d'
    endif

! Writing along Oz
    do ir = 1,nr_cluster
      Vc_a = 0._db
      Rh(:) = 0._db
      Vxc_a(:) = 0._db
      do Lm = 1,nlm_pot_g
        L = Lm_cluster(Lm,1)
        if( Ylm_comp ) then 
          Vc_a = Vc_a + Ylm(L) * real( Vc_cluster_cp(ir,Lm), db )
          Vxc_a(:) = Vxc_a(:) + Ylm(L) * real( Vxc_cluster_cp(ir,:,Lm), db )  
          Rh(:) = Rh(:) + Ylm(L) * real( Rho_cluster_cp(ir,:,Lm), db )
        else
          Vc_a = Vc_a + Ylm(L) * Vc_cluster(ir,Lm) 
          Vxc_a(:) = Vxc_a(:) + Ylm(L) * Vxc_cluster(ir,:,Lm)  
          Rh(:) = Rh(:) + Ylm(L) * Rho_cluster(ir,:,Lm)
        endif 
      end do
      write(3,150) r_cluster(ir)*bohr, Vc_a*Rydb + Vxc_a(:)*Rydb, Vc_a*Rydb, Vxc_a(:)*Rydb, Rh(:), &
                                       quatre_pi * r_cluster(ir)**2 * Rh(:)
    end do
    deallocate( Ylm )
  endif

  return
  100 format(/' V0bdcF_out =',f8.3,' eV')
  102 format('     Vc_out =',f8.3,' eV',/ &
             '   VxcF_out =',f8.3,' eV')
  105 format(/' V0bdcF_out =',2f8.3,' eV')
  107 format('     Vc_out =',f8.3,/ &
             '   VxcF_out =',2f8.3,' eV')
  110 format(/'    Radius  ',2x,1000(7x,a3,a10,7x)) 
  120 format(f11.7,1p,1000(1x,2e13.5))
  130 format(/'    Radius  ',2x,1000(a3,a10)) 
  140 format(f11.7,1p,1000e13.5)
  150 format(f11.7,1p,9e13.5)
end

!***********************************************************************

subroutine Cal_YlmSintdomega(dtheta,Lmax,sin_theta_dtheta,Theta,YlmSintdomega)

  use declarations
  implicit none
  
  integer:: L, Lmax
  
  real(kind=db):: dtheta, cos_theta, ra1, ra2, ra3, sin_theta, sin_theta_dtheta, Theta 
  real(kind=db), dimension(0:Lmax):: YlmSintdomega

  sin_theta = sin( theta )
 
  if( abs( sin_theta ) < eps10 ) then
    sin_theta_dtheta = 1 - cos( dtheta / 2 ) 
  else
    sin_theta_dtheta = sin_theta * dtheta
  endif
 
  if( Lmax > 0 ) cos_theta = cos( theta )

  do L = 0,Lmax
    select case(L)
      case(0)
        YlmSintdomega(L) = sqrt( 1 / quatre_pi ) * sin_theta_dtheta
      case(1)
! Y(l,0)*Sin(theta)
        ra1 = sqrt( 3 / quatre_pi )
        YlmSintdomega(L) = ra1 * cos_theta * sin_theta_dtheta
      case(2)
        ra2 = 0.5_db * sqrt( 45 / quatre_pi )
        YlmSintdomega(L) = ra2 * ( cos_theta**2 - 1 / 3._db ) * sin_theta_dtheta
      case(3)
        ra3 = 0.5_db * sqrt( 175 / quatre_pi )
        YlmSintdomega(L) = ra3 * ( cos_theta**2 - 0.6_db ) * cos_theta * sin_theta_dtheta
      case default
        YlmSintdomega(L) = ( sqrt( 4*L**2 - 1._db ) / L ) * ( cos_theta * YlmSintdomega(L-1) &
                         - ( ( L - 1._db ) / sqrt((2*L-1._db)*(2*L-3._db)) ) * YlmSintdomega(L-2) )
      end select
    end do

    YlmSintdomega(:) = 2 * pi * YlmSintdomega(:)
     
  return
end

!***********************************************************************

subroutine Exchange_correlation(alfpot,E_max,Lm_cluster,Lmax_cluster,nlm_pot_cp,nlm_pot_g,nlm_pot_r,nr_cluster,nspin, &
              r_cluster,Rho_cluster,Rho_cluster_cp,Vxc_cluster,Vxc_cluster_cp,Ylm_comp)

  use declarations
  implicit none
  
  integer:: i, iphi, ir, isp, itheta, L, Lm, Lm_s, Lmax, Lmax_cluster, m, nlm, nlm_pot_cp, nlm_pot_g, nlm_pot_r, &
            nlm_cluster_s, nlm_cluster_s_cp, nlm_cluster_s_g, np, nphi, nr_cluster, nspin, ntheta
  
  integer, dimension(nlm_pot_g,2):: Lm_cluster
  integer, dimension(:,:), allocatable:: Lm_cluster_s

  complex(kind=db), dimension(nr_cluster,nspin,nlm_pot_cp):: Rho_cluster_cp, Vxc_cluster_cp
  complex(kind=db), dimension(:), allocatable:: Ylm_cp_i
  complex(kind=db), dimension(:,:), allocatable:: Ylm_cp
  
  logical:: Magnetic, Ylm_comp, Zero

  real(kind=db):: alfpot, Cos_p, Cos_t, domega_i, dphi, dTheta, E_max, Fac, Omega_t, Phi, r_c, Sin_p, Sin_t, Theta, Tiers
  real(kind=db), dimension(nr_cluster,nspin,nlm_pot_r):: Rho_cluster, Vxc_cluster
  real(kind=db), dimension(nr_cluster):: r_cluster

  real(kind=db), dimension(:), allocatable:: domega, Exc, Rs, Ylm_i
  real(kind=db), dimension(:,:), allocatable:: Rho, Vxc, Ylm

! Selection of the non zero Rho_L
  ir = nr_cluster - 2

  do i = 1,2
    Lm_s = 0
    do Lm = 1,nlm_pot_g
      if( Ylm_comp ) then
        Zero = abs( sum( Rho_cluster_cp(ir,:,Lm) ) ) < eps10
      else
        Zero = abs( sum( Rho_cluster(ir,:,Lm) ) ) < eps10 
      endif
      if( Zero ) cycle
      Lm_s = Lm_s + 1
      if( i == 1 ) cycle
      Lm_cluster_s(Lm_s,:) = Lm_cluster(Lm,:)       
    end do
    nlm_cluster_s_g = Lm_s
    if( i == 1 ) allocate( Lm_cluster_s(nlm_cluster_s_g,2) )
  end do
  
  if( Ylm_comp ) then
    nlm_cluster_s_cp = nlm_cluster_s_g
    nlm_cluster_s = 0
  else
    nlm_cluster_s_cp = 0
    nlm_cluster_s = nlm_cluster_s_g
  endif
    
  Magnetic = nspin == 2
  
  Fac = 0.75 / pi
  Tiers = 1 / 3._db
  
  ntheta = 720 ! 3600
  nphi = ntheta
  
  dTheta = pi / ntheta
  dphi = 2 * pi / nphi
  
  np = ( ntheta - 1 ) * nphi + 2

  allocate( Ylm(np,nlm_cluster_s) )
  allocate( Ylm_cp(np,nlm_cluster_s_cp) )
  
  allocate( Ylm_i(nlm_cluster_s) )
  allocate( Ylm_cp_i(nlm_cluster_s_cp) )
  
  allocate( Rs(np) )
  allocate( domega(np) )

  i = 0
  
  Omega_t = 0._db
    
  do itheta = 0,ntheta

    Theta = itheta * dtheta
    Cos_t = Cos( Theta )
    Sin_t = Sin( Theta )
    
    if( itheta == 0 .or. itheta == ntheta ) then
      domega_i = 2 * pi * ( 1 - cos( dtheta / 2 ) )
    else
      domega_i = sin( Theta ) * dtheta * dphi
    endif

    do iphi = 1,nphi

      i = i + 1
      
      domega(i) = domega_i
      
      Phi = ( iphi - 1 ) * dphi
      Cos_p = cos( Phi )
      Sin_p = sin( phi )
   
      call Cal_Ylm(Cos_p,Cos_t,Lm_cluster_s,Lmax_cluster,nlm_cluster_s,nlm_cluster_s_cp,nlm_cluster_s_g,Sin_p,Sin_t,Ylm_i, &
                   Ylm_comp,Ylm_cp_i)
      if( Ylm_comp ) then
        Ylm_cp(i,:) = Ylm_cp_i(:)
      else
        Ylm(i,:) = Ylm_i(:)
      endif

      Omega_t = Omega_t + domega(i)
      
      if( itheta == 0 .or. itheta == ntheta ) exit 

    end do
  end do

  Omega_t =  Omega_t / ( 4 * pi )
  domega(:) = domega(:) / Omega_t
  
  deallocate( Ylm_cp_i, Ylm_i )

  allocate( Exc(np) )
  allocate( Rho(np,nspin) )
  allocate( Vxc(np,nspin) )
   
  do ir = 1,nr_cluster-1
    r_c = r_cluster(ir)
      
    call clmax(E_max,r_c,-1,Lmax,-1,.true.)
    Lmax = max( Lmax, 8 )
    Lmax = min( Lmax, Lmax_cluster )

    Rho(:,:) = 0._db

    do Lm_s = 1,nlm_cluster_s_g
      L = Lm_cluster_s(Lm_s,1)
      if( L > Lmax ) exit
      m = Lm_cluster_s(Lm_s,2)
      do Lm = 1,nlm_pot_g
        if( L == Lm_cluster(Lm,1) .and. m == Lm_cluster(Lm,2) ) exit
      end do  
      
      do isp = 1,nspin
        if( Ylm_comp ) then
          Rho(:,isp) = Rho(:,isp) + real( Ylm_cp(:,Lm_s) * Rho_cluster_cp(ir,isp,Lm), db )
        else
          Rho(:,isp) = Rho(:,isp) + Ylm(:,Lm_s) * Rho_cluster(ir,isp,Lm)
        endif
      end do

    end do  
           
    Rho(:,:) = max( Rho(:,:), 0._db )

    if( Lm_s > nlm_cluster_s_g ) then
      nlm = nlm_cluster_s_g
    else
      nlm = Lm_s - 1
    endif
 
 ! Fermi radius
    if( nspin == 1 ) then
      Rs(:) = ( Fac / Rho(:,1) )**Tiers
    else
      Rs(:) = ( Fac / ( Rho(:,1) + Rho(:,2) ) )**Tiers
    endif

! Calculation of exchange-correlation potential
    call Potxc(Magnetic,np,np,nspin,alfpot,Rho,Rs,Vxc,Exc)

! We force Vxc to have the same expansion in spherical harmonics than Rho

    do isp = 1,nspin
      do Lm_s = 1,nlm          
        L = Lm_cluster_s(Lm_s,1)
        if( L > 0 .and. ir == nr_cluster - 1 ) cycle  
        m = Lm_cluster_s(Lm_s,2)
        do Lm = 1,nlm_pot_g
          if( L == Lm_cluster(Lm,1) .and. m == Lm_cluster(Lm,2) ) exit
        end do

        if( Ylm_comp ) then
          Vxc_cluster_cp(ir,isp,Lm) = Vxc_cluster_cp(ir,isp,Lm) + sum( domega(:) * conjg( Ylm_cp(:,Lm_s) ) * Vxc(:,isp) )
        else
          Vxc_cluster(ir,isp,Lm) = Vxc_cluster(ir,isp,Lm) + sum( domega(:) * Ylm(:,Lm_s) * Vxc(:,isp) )
        endif   
      end do
    end do    

  end do
  
  deallocate( domega, Exc, Lm_cluster_s, Rho, Rs, Vxc, Ylm, Ylm_cp )

  if( Ylm_comp ) then
    Vxc_cluster_cp(nr_cluster,:,1) = Vxc_cluster_cp(nr_cluster-1,:,1)
    if( nlm_pot_cp > 1 ) Vxc_cluster_cp(nr_cluster,:,2:nlm_pot_cp) = ( 0._db, 0._db )
  else
    Vxc_cluster(nr_cluster,:,1) = Vxc_cluster(nr_cluster-1,:,1)
    if( nlm_pot_r > 1 ) Vxc_cluster(nr_cluster,:,2:nlm_pot_r) = 0._db
  endif 
          
  return
end

!********************************************************************************************

! Calculation of Ylm

subroutine Cal_Ylm(Cos_p,Cos_t,Lm_cluster_s,Lmax,nlm_cluster_s,nlm_cluster_s_cp,nlm_cluster_s_g,Sin_p,Sin_t,Ylm, &
                   Ylm_comp,Ylm_cp)

  use declarations
  implicit none

  integer:: L, Lm, Lm_c, Lm0, Lm1, Lm2, Lmax, m, nlm_cluster_s, nlm_cluster_s_cp, nlm_cluster_s_g

  integer, dimension(nlm_cluster_s_g,2):: Lm_cluster_s

  complex(kind=db):: exp_p, Sin_t_exp_p
  complex(kind=db), dimension((Lmax+1)*(Lmax+2)/2):: Ylmc
  complex(kind=db), dimension(nlm_cluster_s_cp):: Ylm_cp
  
  logical:: Ylm_comp

  real(kind=db):: Cos_t, Cos_p, Cot_t, f, g, Sin_t, Sin_p
  real(kind=db), dimension(nlm_cluster_s):: Ylm

  if( abs( Sin_t ) < eps10 ) then
    exp_p = ( 1._db, 0._db )
    Sin_t_exp_p = ( 0._db, 0._db )
  else
    Cot_t = Cos_t / Sin_t
    exp_p = cmplx( Cos_p, Sin_p, db )
    Sin_t_exp_p = Sin_t * exp_p
  endif

! The index of Ylmc is lm = L(L+1)/2 + 1 + m
! Only m >= 0 are calculated
  
! Y(0,0)
  Ylmc(1) = 1._db / sqrt( quatre_pi )

! Y(L,L)
  do L = 1,Lmax
    Lm = ( ( L + 1 ) * ( L + 2 ) ) / 2
    Lm1 = ( L * ( L + 1 ) ) / 2
    f = - sqrt( 1 + 0.5_db / L )
    Ylmc(Lm) = f * Sin_t_exp_p * Ylmc(Lm1)
  end do

! Y(1,0)
  if( Lmax > 0 ) Ylmc(2) = sqrt( 3._db / quatre_pi ) * Cos_t

! Y(L,L-1)
  do L = 2,Lmax
    Lm = ( ( L + 1 ) * ( L + 2 ) ) / 2 - 1
    Lm1 = ( L**2 + L ) / 2 - 1
    f = - sqrt( (2*L + 1._db) / (2*L - 2) )
    Ylmc(Lm) = f * Sin_t_exp_p * Ylmc(Lm1)
  end do

! Y(L,m)
  exp_p = conjg( exp_p )

  do L = 2,Lmax
    Lm0 = ( L**2 + L ) / 2
    do m = L-2,1,-1
      Lm = Lm0 + m + 1
      if( abs( Sin_t ) < eps10 ) then
        Ylmc(Lm) = ( 0._db, 0._db )
      else
        f = - 2 * (m + 1._db) / sqrt( L*(L+1._db) - m*(m+1._db) )
        g = - sqrt( ( L*(L+1._db) - (m+1._db)*(m+2._db) ) / ( L*(L+1._db) - m*(m+1._db) ) )
        Ylmc(Lm) = exp_p * ( Cot_t * f * Ylmc(Lm+1) + exp_p * g * Ylmc(Lm+2) )
      endif
    end do
  end do

! Y(L,0)
  do L = 2,Lmax
    Lm = ( L * ( L + 1 ) ) / 2 + 1
    Lm1 = ( ( L - 1 ) * L ) / 2 + 1
    Lm2 = ( ( L - 2 ) * ( L - 1) ) / 2 + 1
    f = sqrt( (2*L + 1._db) / (2*L - 1._db) ) * (2*L - 1._db) / L
    g = - sqrt( (2*L + 1._db) / (2*L - 3._db) ) * (L - 1._db) / L
    Ylmc(Lm) = f * Cos_t * Ylmc(Lm1) + g * Ylmc(Lm2)
  end do

  do Lm = 1,nlm_cluster_s_g
    L = Lm_cluster_s(Lm,1)
    m = Lm_cluster_s(Lm,2)
    Lm_c = (L**2 + L ) / 2 + 1 + abs( m )
    
    if( Ylm_comp ) then

      if( m >= 0 ) then
        Ylm_cp(Lm) = Ylmc(Lm_c)
      else
        Ylm_cp(Lm) = (-1)**m * conjg( Ylmc(Lm_c) )
      endif

    else

      if( m < 0 ) then
        Ylm(Lm) = (-1)**m * sqrt_2 * aimag( Ylmc(Lm_c) )
      elseif( m == 0 ) then
        Ylm(Lm) = real( Ylmc(Lm_c), db )
      else
        Ylm(Lm) = (-1)**m * sqrt_2 * real( Ylmc(Lm_c), db )
      endif
      
    endif

  end do
  
  return
end

!********************************************************************************************************

! Calculation of potential above Fermi level, for cluster approach

subroutine Potex_cluster(Alfpot,Enervide,icheck,Lm_cluster,nlm_pot_cp,nlm_pot_g,nlm_pot_r,nr_cluster,nspin, &
                   Optic,r_cluster,Rho_cluster,Rho_cluster_cp,V_intmax,V0bdc_out,Vc_cluster,Vc_cluster_cp,Vr_cluster, &
                   Vr_cluster_cp,Vxc_cluster,Vxc_cluster_cp,Ylm_comp)

  use declarations
  implicit none

  integer:: i, icheck, ir, isp, L, Length, Lm, m, nlm_pot_cp, nlm_pot_g, nlm_pot_r, nr, nr_cluster, nspin

  integer, dimension(nlm_pot_g,2):: Lm_cluster

  character(len=10):: mot10
  character(len=10), dimension(:), allocatable:: Col_name, Col_name_sp

  complex(kind=db), dimension(nr_cluster,nlm_pot_cp):: Vc_cluster_cp
  complex(kind=db), dimension(nr_cluster,nspin,nlm_pot_cp):: Rho_cluster_cp, Vxc_cluster_cp
  complex(kind=db), dimension(nr_cluster,nspin,nlm_pot_cp):: Vr_cluster_cp
  
  logical:: Optic, Ylm_comp

  real(kind=db):: Alfpot, Enervide, F, Fac, Tiers, V_intmax, Vm  
  real(kind=db), dimension(nspin):: V0bdc_out
  real(kind=db), dimension(nr_cluster):: r_cluster, Rs, Vc, Vr, Vxc
  real(kind=db), dimension(nr_cluster,nlm_pot_r):: Vc_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_pot_r):: Rho_cluster, Vxc_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_pot_r):: Vr_cluster

  F = 1 / sqrt( 4 * pi )
  Fac = 0.75 / pi
  Fac = Fac / F  ! because Rho_cluster(:,1,1) would have to be multiplied by Y_00 to be the real rho 
  Tiers = 1 / 3._db

  if( Ylm_comp ) then
    if( nspin == 1 ) then
      Rs(:) = ( Fac / Rho_cluster_cp(:,1,1) )**Tiers
    else
      Rs(:) = ( Fac / ( Rho_cluster_cp(:,1,1) + Rho_cluster_cp(:,2,1) ) )**Tiers
    endif
  else
    if( nspin == 1 ) then
      Rs(:) = ( Fac / real( Rho_cluster(:,1,1), db ) )**Tiers
    else
      Rs(:) = ( Fac / ( real( Rho_cluster(:,1,1), db ) + real( Rho_cluster(:,2,1), db ) ) )**Tiers
    endif
  endif
    
  if( Alfpot < eps6 .and. .not. Optic )  then

    if( Ylm_comp ) then
      Vc(:) = real( Vc_cluster_cp(:,1), db )
    else
      Vc(:) = Vc_cluster(:,1)
    endif
    
    do isp = 1,nspin
      if( Ylm_comp ) then
        Vxc(:) = real( Vxc_cluster_cp(:,isp,1), db )
      else
        Vxc(:) = Vxc_cluster(:,isp,1)
      endif
      call subpotex(nr_cluster,Vr,Vc,Vxc,Rs,Enervide)
      if( Ylm_comp ) then
        Vr_cluster_cp(1:nr_cluster,isp,1) = cmplx( Vr(1:nr_cluster), 0._db, db )
      else
        Vr_cluster(1:nr_cluster,isp,1) = Vr(1:nr_cluster)
      endif
      
      nr = nr_cluster - 2
      if( Ylm_comp ) then
        do Lm = 2,nlm_pot_cp
          Vr_cluster_cp(1:nr,isp,Lm) = Vc_cluster_cp(1:nr,Lm) + ( ( Vr(1:nr) - real( Vc_cluster_cp(1:nr,1), db ) ) / Vxc(1:nr) ) &
                                                              * Vxc_cluster_cp(1:nr,isp,Lm)
          Vr_cluster_cp(nr+1:nr_cluster,isp,Lm) = ( 0._db, 0._db )
        end do
      else
        do Lm = 2,nlm_pot_r
          Vr_cluster(1:nr,isp,Lm) = Vc_cluster(1:nr,Lm) + ( ( Vr(1:nr) - Vc_cluster(1:nr,1) ) / Vxc(1:nr) ) &
                                                        * Vxc_cluster(1:nr,isp,Lm)
          Vr_cluster(nr+1:nr_cluster,isp,Lm) = 0._db
        end do
      endif
      
    end do

  else

    do isp = 1,nspin
      if( Ylm_comp ) then
        do Lm = 1,nlm_pot_cp
          Vr_cluster_cp(1:nr_cluster,isp,Lm) = Vc_cluster_cp(1:nr_cluster,Lm) + Vxc_cluster_cp(1:nr_cluster,isp,Lm)
        end do
      else
        do Lm = 1,nlm_pot_r
          Vr_cluster(1:nr_cluster,isp,Lm) = Vc_cluster(1:nr_cluster,Lm) + Vxc_cluster(1:nr_cluster,isp,Lm)
        end do
      endif
    end do

  endif

  if( V_intmax < 1000._db ) then
    Vm = V_intmax / F  
    if( Ylm_comp ) then
      Vr_cluster_cp(:,:,1) = cmplx( min( real( Vr_cluster_cp(:,:,1), db ), Vm ), 0._db, db )  
    else
      Vr_cluster(:,:,1) = min( Vr_cluster(:,:,1), Vm )
    endif
  endif  
  
  if( Ylm_comp ) then
    V0bdc_out(:) = F * real( Vr_cluster_cp(nr_cluster,:,1), db )
  else
    V0bdc_out(:) = F * Vr_cluster(nr_cluster,:,1)
  endif
  
  if( icheck > 2 )  then
  
    write(3,100)
    
    if( nspin == 1 ) then
      write(3,105) V0bdc_out(:)*rydb
    else
      write(3,107) V0bdc_out(:)*rydb
    endif
    
    allocate( Col_name(nlm_pot_g) )
    do Lm = 1,nlm_pot_g
      mot10 = '(' 
      L = Lm_cluster(Lm,1)
      m = Lm_cluster(Lm,2)
      call ad_number(L,mot10,10)
      Length = len_trim(mot10) + 1
      mot10(Length:Length) = ','
      call ad_number(m,mot10,10)
      Length = len_trim(mot10) + 1
      if( Length <= 10 ) mot10(Length:Length) = ')'
      Col_name(Lm) = mot10
    end do

    if( nspin == 2 ) then
      allocate( Col_name_sp(nspin*nlm_pot_g) )
      i = 0
      do Lm = 1,nlm_pot_g
        do isp = 1,nspin
          i = i + 1 
          mot10 = '(' 
          L = Lm_cluster(Lm,1)
          m = Lm_cluster(Lm,2)
          call ad_number(L,mot10,10)
          Length = len_trim(mot10) + 1
          mot10(Length:Length) = ','
          call ad_number(m,mot10,10)
          Length = len_trim(mot10) + 1
          if( Length <= 9 ) mot10(Length:Length) = ')'
          Length = len_trim(mot10) + 1
          if( Length <= 10 ) then 
            if( isp == 1 ) then
              mot10(Length:Length) = 'u'
            else
              mot10(Length:Length) = 'u'
            endif
          endif
          Col_name_sp(i) = mot10
        end do
      end do
    endif
      
    do ir = 1,nr_cluster
      if( Ylm_comp ) then
        if( nspin == 1 ) then
          if( ir == 1 ) write(3,110) ( ' Vr', Col_name(Lm), Lm = 1,nlm_pot_g )
          write(3,120) r_cluster(ir)*bohr, ( Vr_cluster_cp(ir,:,Lm)*rydb, Lm = 1,nlm_pot_g )
        else
          if( ir == 1 ) write(3,110) ( ' Vr', Col_name_sp(Lm), Lm = 1,nspin*nlm_pot_g )
          write(3,120) r_cluster(ir)*bohr, ( Vr_cluster_cp(ir,:,Lm)*rydb, Lm = 1,nlm_pot_g )
        endif
      else
        if( nspin == 1 ) then
          if( ir == 1 ) write(3,130) ( ' Vr', Col_name(Lm), Lm = 1,nlm_pot_g ) 
          write(3,140) r_cluster(ir)*bohr, ( Vr_cluster(ir,:,Lm)*rydb, Lm = 1,nlm_pot_g )
        else
          if( ir == 1 ) write(3,130) ( ' Vr', Col_name_sp(Lm), Lm = 1,nspin*nlm_pot_g )
          write(3,140) r_cluster(ir)*bohr, ( Vr_cluster(ir,:,Lm)*rydb, Lm = 1,nlm_pot_g )
        endif 
      endif 
    end do

  endif
  
  return
  100 format(/' ----- Potex_cluster ',80('-'))
  105 format(/' V0bdc_out =',f9.5,' eV')
  107 format(/' V0bdc_out =',2f9.5,' eV')
  110 format(/'    Radius  ',2x,1000(7x,a3,a10,7x)) 
  120 format(f11.7,1p,1000(1x,2e13.5))
  130 format(/'    Radius  ',2x,1000(a3,a10)) 
  140 format(f11.7,1p,1000e13.5)
end

!**********************************************************************

! Calculation of radial integrals for the FDM complex case
! Contain both regular and irregular solutions.
! Called by Tenseur_car and S_nrixs_cal

subroutine Radial_FDM_comp_cluster(E_comp,Ecinetic,Eimag,Energ,Enervide,Eseuil,Hubb_a, &
              Hubb_d,icheck,ip_max,ip0,iso,Lm_cluster,Lmaxso,Lso,m_hubb,mso,nbseuil,ngrph, &
              nlm_pot_cp,nlm_pot_g,nlm_pot_r,nlm_probe,nlmso0,nr_cluster,nr_core,NRIXS,nso1,nspin,nspino, &
              nspinp,numat_abs,Numerov,Psi_core,r_core,r_cluster,Rad_ioRIoi,Relativiste,Rmtg_abs,Spinorbite, &
              V_hubb,V0bdc_out,Vr_cluster,Vr_cluster_cp,Ylm_comp)

  use declarations
  implicit none

  integer:: i, icheck, ip, ip_max, ip0, Lmax_probe, Lmaxso, m_hubb, nbseuil, ngrph, nlm_E, nlm_pot_cp, nlm_pot_g, nlm_pot_r, &
    nlm_probe, nr_cluster, nr_core, nso1, nspin, nspino, nspinp, numat_abs

  integer, dimension(nlm_pot_g,2):: Lm_cluster
  integer, dimension(ngrph):: nlmso0
  integer, dimension(nso1,ngrph):: iso, Lso, mso
  integer, dimension(nso1,3,ngrph):: Lmso

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nlm_probe,nlm_probe,nspinp,ip_max-ip0+1,ip_max-ip0+1,nbseuil):: Rad_ioRIoi
  complex(kind=db), dimension(nr_cluster,nspin,nlm_pot_cp):: Vr_cluster_cp

  logical:: E_comp, Final_tddft, Hubb_a, Hubb_d, NRIXS, Numerov, Relativiste, Spinorbite, Ylm_comp

  real(kind=db):: Eimag, Energ, Enervide, Ephoton, Rmtg_abs
  real(kind=db), dimension(nspin):: Ecinetic, V0bdc_out
  real(kind=db), dimension(nr_cluster):: r_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_pot_r):: Vr_cluster
  real(kind=db), dimension(nr_core):: r_core
  real(kind=db), dimension(nbseuil):: Eseuil, Vecond
  real(kind=db), dimension(nr_core,nbseuil):: Psi_core

  real(kind=db), dimension(:,:), allocatable:: r_or_bess
  real(kind=db), dimension(:,:,:,:,:), allocatable:: u_irreg_i, u_irreg_r, u_reg_i, u_reg_r

  Rad_ioRIoi(:,:,:,:,:,:) = (0._db, 0._db )

  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag,db) )

  if( icheck > 1 ) then
    write(3,110) 
    write(3,120) Energ*rydb, Enervide*rydb
    if( nspin == 1 ) then
      write(3,130) Ecinetic(:)*rydb
      write(3,140) V0bdc_out(:)*rydb
      write(3,150) konde(:)
    else
      write(3,135) Ecinetic(:)*rydb
      write(3,145) V0bdc_out(:)*rydb
      write(3,155) konde(:)
    endif
  endif

  Lmax_probe = nint( sqrt( real( nlm_probe ) ) ) - 1
  nlm_E = ( Lmaxso + 1 )**2

  Lmso(:,1,:) = Lso(:,:)
  Lmso(:,2,:) = mso(:,:)
  Lmso(:,3,:) = iso(:,:)

  allocate( u_irreg_i(nr_core,nlm_E,nlm_probe,nspinp,nspino) )
  allocate( u_irreg_r(nr_core,nlm_E,nlm_probe,nspinp,nspino) )
  allocate( u_reg_i(nr_core,nlm_E,nlm_probe,nspinp,nspino) )
  allocate( u_reg_r(nr_core,nlm_E,nlm_probe,nspinp,nspino) )
  
  call Cluster_solutions(E_comp,Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d, &
            icheck,Lm_cluster,Lmaxso,Lmso,m_hubb,ngrph,nlm_E,nlm_pot_cp,nlm_pot_g,nlm_pot_r,nlm_probe,nlmso0, &
            nr_cluster,nr_core,nso1,nspin,nspino,nspinp,numat_abs,Numerov,r_cluster,Relativiste,Rmtg_abs, &
            Spinorbite,u_irreg_i,u_irreg_r,u_reg_i,u_reg_r,V_hubb,Vr_cluster,Vr_cluster_cp,Ylm_comp)
                 
  allocate( r_or_bess(nr_core,ip0:ip_max) )

  do ip = ip0,ip_max
    select case(ip)
      case(0)
        r_or_bess(:,ip) = 1._db
      case(1)
        r_or_bess(:,ip) = r_core(:)
      case default
        r_or_bess(:,ip) = r_core(:)**ip
    end select
  end do    

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

  call Integral_FDM_comp_cluster(icheck,ip_max,ip0,nbseuil,nlm_E,nlm_probe,nr_core,NRIXS, &
            nspino,nspinp,Psi_core,r_core,r_or_bess,Rad_ioRIoi,Spinorbite,u_irreg_i,u_irreg_r,u_reg_i,u_reg_r,Vecond)

  deallocate( r_or_bess, u_irreg_i, u_irreg_r, u_reg_i, u_reg_r )

  return
  110 format(/' ----- Radial_FDM_comp_cluster ',80('-'))
  120 format(/' Energ     =',f10.3,' eV,  Enervide =',f10.3,' eV')
  130 format(' Ecinetic  =',f10.3,' eV')
  135 format(' Ecinetic  =',2f10.3,' eV')
  140 format(' V0bdc_out =',f10.3,' eV')
  145 format(' V0bdc_out =',2f10.3,' eV')
  150 format(' konde     =',2f10.5,' /u.a.')
  155 format(' konde     =',4f10.5,' /u.a.')
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

!***********************************************************************

! Calculation of the cluster regular and irregular solutions

subroutine Cluster_solutions(E_comp,Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d, &
            icheck,Lm_cluster,Lmaxso,Lmso,m_hubb,ngrph,nlm_E,nlm_pot_cp,nlm_pot_g,nlm_pot_r,nlm_probe,nlmso0, &
            nr_cluster,nr_core,nso1,nspin,nspino,nspinp,numat_abs,Numerov,r_cluster,Relativiste,Rmtg_abs, &
            Spinorbite,u_irreg_i,u_irreg_r,u_reg_i,u_reg_r,V_hubb,Vr_cluster,Vr_cluster_cp,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, ir, Lm, Lmaxso, m_hubb, ngrph, nlm_E, nlm_cp_E, nlm_pot_cp, nlm_pot_g, nlm_pot_r, &
    nlm_r_E, nlm_probe, nr_cluster, nr_core, nso1, nspin, nspino, nspinp, numat_abs

  integer, dimension(nlm_pot_g,2):: Lm_cluster
  integer, dimension(ngrph):: nlmso0
  integer, dimension(nso1,3,ngrph):: Lmso
  
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nr_cluster+1,nspin,nlm_pot_cp):: Vr_cluster_cp
  complex(kind=db), dimension(:,:,:,:), allocatable:: V_cp
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: u, us

  logical:: E_comp, Gaussian, Hubb_a, Hubb_d, Numerov, Relativiste, Spinorbite, Ylm_comp

  real(kind=db):: Eimag, Enervide, Rmtsd, Rmtg_abs

  real(kind=db), dimension(nspin):: Ecinetic
  real(kind=db), dimension(nr_cluster):: f2, r_cluster
  real(kind=db), dimension(nr_cluster,nspin):: g0, gm, gmi, gp, gpi
  real(kind=db), dimension(nr_cluster,nspino):: gso
  real(kind=db), dimension(nr_cluster,nspin,nlm_pot_r):: Vr_cluster
  real(kind=db), dimension(nr_core,nlm_E,nlm_probe,nspinp,nspino):: u_irreg_i, u_irreg_r, u_reg_i, u_reg_r

  real(kind=db), dimension(:,:,:,:), allocatable:: V

  Rmtsd = R_cluster(nr_cluster)
  
  Gaussian = .true.
!  Gaussian = .false.

  konde(:) = sqrt( cmplx(Ecinetic(:), Eimag, db) )

  if( Ylm_comp ) then
    nlm_cp_E = nlm_E  
    nlm_r_E = 0  
  else
    nlm_r_E = nlm_E  
    nlm_cp_E = 0  
  endif
  
  allocate( V_cp(nr_cluster,nspin,nlm_cp_E,nlm_cp_E) )
  allocate( V(nr_cluster,nspin,nlm_r_E,nlm_r_E) )
 
  call Cal_V_LmLmp(icheck,Lm_cluster,nlm_E,nlm_cp_E,nlm_pot_cp,nlm_pot_r,nlm_pot_g,nlm_r_E,nr_cluster, &
                   nspin,r_cluster,V,V_cp,Vr_cluster,Vr_cluster_cp,Ylm_comp)

! Just real potential
  if( .not. Numerov .and. .not. Gaussian ) &
    call coef_sch_rad(Enervide,f2,g0,gm,gp,gso,nlm_r_E,nr_cluster,nspin,nspino,numat_abs,r_cluster,Relativiste,Spinorbite,V)

  allocate( u(nr_cluster,nlm_E,nlm_E,nspinp,nspino) )
  
  if( Gaussian ) then
!    call Sch_radial_reg_irreg_gaussian_old(Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d,icheck,.false.,konde,Lmso,Lmaxso,m_hubb, &
!         ngrph,nlm_cp_E,nlm_r_E,nlm_E,nlmso0,nr_cluster,nso1,nspin,nspino,nspinp,numat_abs,r_cluster, &
!         Relativiste,Rmtg_abs,Spinorbite,u,V,V_cp,V_hubb,Ylm_comp)

    call Sch_radial_reg_irreg_gaussian(Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d,icheck,.false.,konde,Lmso,Lmaxso,m_hubb, &
         ngrph,nlm_cp_E,nlm_r_E,nlm_E,nlmso0,nr_cluster,nso1,nspin,nspino,nspinp,numat_abs,Numerov,r_cluster, &
         Relativiste,Rmtg_abs,Spinorbite,u,V,V_cp,V_hubb,Ylm_comp)

  elseif( Numerov ) then  
    call Sch_radial_reg_Numerov(Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
         nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
         Relativiste,Rmtg_abs,Spinorbite,u,V,V_cp,V_hubb,Ylm_comp)
  
  else
    gpi(1:nr_cluster-1,:) = 1 / gp(1:nr_cluster-1,:)
    call Sch_radial_reg_cluster(Ecinetic,E_comp,Eimag,f2,g0,gm,gpi,gso,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
         nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
         Relativiste,Rmtg_abs,Spinorbite,u,V,V_cp,V_hubb,Ylm_comp)
  endif
  
  do ir = 1,nr_core
    do Lm = 1,min(nlm_probe, nlm_E)
      u_reg_r(ir,:,Lm,:,:) = real( u(ir,:,Lm,:,:), db )
      u_reg_i(ir,:,Lm,:,:) = aimag( u(ir,:,Lm,:,:) )
    end do
    do Lm = min(nlm_probe, nlm_E) + 1, nlm_probe 
      u_reg_r(ir,:,Lm,:,:) = 0._db
      u_reg_i(ir,:,Lm,:,:) = 0._db
    end do
  end do
  
  deallocate( u )

  allocate( us(nr_cluster,nlm_E,nlm_E,nspinp,nspino) )

  if( Gaussian ) then
!    call Sch_radial_reg_irreg_gaussian_old(Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d,icheck,.true.,konde,Lmso,Lmaxso,m_hubb, &
!         ngrph,nlm_cp_E,nlm_r_E,nlm_E,nlmso0,nr_cluster,nso1,nspin,nspino,nspinp,numat_abs,r_cluster, &
!         Relativiste,Rmtg_abs,Spinorbite,us,V,V_cp,V_hubb,Ylm_comp)

    call Sch_radial_reg_irreg_gaussian(Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d,icheck,.true.,konde,Lmso,Lmaxso,m_hubb, &
         ngrph,nlm_cp_E,nlm_r_E,nlm_E,nlmso0,nr_cluster,nso1,nspin,nspino,nspinp,numat_abs,Numerov,r_cluster, &
         Relativiste,Rmtg_abs,Spinorbite,us,V,V_cp,V_hubb,Ylm_comp)

  elseif( Numerov ) then  
    call Sch_radial_irreg_Numerov(Eimag,Enervide,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
       nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
       Relativiste,Rmtg_abs,Spinorbite,us,V,V_cp,V_hubb,Ylm_comp)
  
!    call Sch_radial_irreg_gaussian(Eimag,Enervide,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
!         nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
!         Relativiste,Rmtg_abs,Spinorbite,us,V,V_cp,V_hubb,Ylm_comp)

  else
    gmi(1:nr_cluster-1,:) = 1 / gm(1:nr_cluster-1,:)
    call Sch_radial_irreg_cluster(E_comp,Eimag,f2,g0,gmi,gp,gso,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
       nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
       Rmtg_abs,Spinorbite,us,V,V_cp,V_hubb,Ylm_comp)
  endif

  do ir = 1,nr_core
    do Lm = 1,min(nlm_probe, nlm_E)
      u_irreg_r(ir,:,Lm,:,:) = real( us(ir,:,Lm,:,:), db )
      u_irreg_i(ir,:,Lm,:,:) = aimag( us(ir,:,Lm,:,:) )
    end do
    do Lm = min(nlm_probe, nlm_E) + 1, nlm_probe 
      u_irreg_r(ir,:,Lm,:,:) = 0._db
      u_irreg_i(ir,:,Lm,:,:) = 0._db
    end do
  end do

  deallocate( us, V, V_cp )

  return
end

!***********************************************************************

Subroutine Cal_V_LmLmp(icheck,Lm_cluster,nlm_E,nlm_cp_E,nlm_pot_cp,nlm_pot_r,nlm_pot_g,nlm_r_E,nr_cluster, &
                   nspin,r_cluster,V,V_cp,Vr_cluster,Vr_cluster_cp,Ylm_comp)

  use declarations
  implicit none

  integer:: i_cycle, icheck, ir, L1, L2, L3, Lm, Lm1, Lm2, Lm3, Lmax, m1, m2, m3, &
    nlm_pot_r, nlm_pot_cp, nlm_cp_E, nlm_pot_g, nlm_E, nlm_r_E, nlma, nr_cluster, nspin

  integer, dimension(nlm_pot_g,2):: Lm_cluster
  integer, dimension(:,:), allocatable:: LmLmp

  logical:: Ylm_comp

  complex(kind=db), dimension(nr_cluster,nspin,nlm_pot_cp):: Vr_cluster_cp
  complex(kind=db), dimension(nr_cluster,nspin,nlm_cp_E,nlm_cp_E):: V_cp

  real(kind=db):: Fac, g, Gaunt_r, gauntcp
  real(kind=db), dimension(nr_cluster):: r_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_pot_r):: Vr_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_r_E,nlm_r_E):: V

  Fac = 1 / sqrt( quatre_pi )
  
  if( Ylm_comp ) then
    V_cp(:,:,:,:) = 0._db
  else
    V(:,:,:,:) = 0._db
  endif

  if( icheck > 2 ) write(3,'(/A)') ' L1 m1 L2 m2 L3 m3    Gaunt           V'          

  Lmax = nint( sqrt( real( nlm_E ) ) ) - 1
  do L1 = 0,Lmax
    do m1 = -L1, L1
      Lm1 = L1**2 + L1 + 1 + m1
      do L2 = 0,Lmax
        do m2 = -L2, L2
          Lm2 = L2**2 + L2 + 1 + m2

          do Lm3 = 1,nlm_pot_r
            L3 = Lm_cluster(Lm3,1)
            m3 = Lm_cluster(Lm3,2)
      
            if( Ylm_comp ) then
              g = gauntcp(L1,m1,L2,m2,L3,m3)
            else
              g = Gaunt_r(L1,m1,L2,m2,L3,m3)
            endif

            if( abs( g ) < eps10 ) cycle

            if( Ylm_comp ) then
              V_cp(:,:,Lm1,Lm2) = V_cp(:,:,Lm1,Lm2) + g * Vr_cluster_cp(:,:,Lm3)
            else
              V(:,:,Lm1,Lm2) = V(:,:,Lm1,Lm2) + g * Vr_cluster(:,:,Lm3)
            endif
            if( icheck > 2 .and. abs( Vr_cluster(nr_cluster-3,1,Lm3)) > eps10 ) &
                              write(3,'(6i3,f11.7,f13.5)') L1, m1, L2, m2, L3, m3, g, V(nr_cluster-3,1,Lm1,Lm2)*rydb        

          end do
        end do
      end do
    end do
  end do

  if( icheck > 2 ) then
    do i_cycle = 1,2
      if( i_cycle == 2 ) allocate( LmLmp(nlma,6) )
      Lm = 0
      Lm1 = 0
      do L1 = 0,Lmax
        do m1 = - L1,L1
          Lm1 = Lm1 + 1
          Lm2 = 0
          do L2 = 0,Lmax
            do m2 = - L2,L2
              Lm2 = Lm2 + 1
              if( Ylm_comp ) then
                if( abs( V_cp(nr_cluster-3,1,Lm1,Lm2) ) < eps10 ) cycle
              else
                if( abs( V(nr_cluster-3,1,Lm1,Lm2) ) < eps10 ) cycle
              endif
              Lm = Lm + 1
              if( i_cycle == 1 ) cycle
              LmLmp(Lm,1) = L1 
              LmLmp(Lm,2) = m1 
              LmLmp(Lm,3) = L2 
              LmLmp(Lm,4) = m2 
              LmLmp(Lm,5) = Lm1 
              LmLmp(Lm,6) = Lm2 
            end do
          end do
        end do
      end do
      nlma = Lm
    end do
    
    write(3,110)
    if( nspin == 2 .and. nlm_pot_g == 1 ) then
      write(3,'(/A)') '    Radius     V(up)      V(dn)'
    elseif( nspin == 2 .and. nlm_pot_g > 1 ) then
      write(3,'(/A,A)') '    Radius  V(up,Lm)   V(dn,Lm)    Lm = 1,nlm_pot_r'
    elseif( nspin == 1 .and. nlm_pot_g == 1 ) then
      write(3,'(/A)') '    Radius      V'
    else
      if( Ylm_comp ) then
        write(3,120) ( LmLmp(Lm,1:4), Lm = 1,nlma )
      else
        write(3,130) ( LmLmp(Lm,1:4), Lm = 1,nlma )
      endif
    endif
    do ir = 1,nr_cluster
      if( Ylm_comp ) then
        write(3,150) r_cluster(ir)*bohr, ( V_cp(ir,:,LmLmp(Lm,5),LmLmp(Lm,6))*rydb, Lm = 1,nlma )
      else
        write(3,150) r_cluster(ir)*bohr, ( V(ir,:,LmLmp(Lm,5),LmLmp(Lm,6))*rydb, Lm = 1,nlma )
      endif
    end do
    deallocate( LmLmp )
  endif

  return
  110 format(/' ------- Potential V(Lm1,Lm2) ', 80('-'))
  120 format(/'    Radius  ',10000(7x,'V(',i1,',',i2,';',i1,',',i2,')',8x))
  130 format(/'    Radius  ',10000('V(',i1,',',i2,';',i1,',',i2,')',1x))
  140 format(f11.7,1p,10000(1x,2e13.5))
  150 format(f11.7,1p,10000e13.5)
end

!***********************************************************************

! Resolution of the radial Schrodinger equation. Include the case of the non spherical potential, use Numerov method.
! In u(ir,Lmb,Lm,isp,isol) : (Lmb, isol) define the basis state
!                            (L, m, Lm = L**2 + L + 1 + m, isp) are the real moment and spin

subroutine Sch_radial_reg_Numerov(Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
         nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
         Relativiste,Rmtg_abs,Spinorbite,u,V,V_cp,V_hubb,Ylm_comp)

  use declarations
  implicit none

  integer:: i0, icheck, im, ip, ir, isol, isp, isq, isr, L, L_hubbard, L2, Lm, Lmax, Lmaxso, Lmp, m, m_hubb, mp, ms, &
    n, nlm_cp_E, nlm_r_E, nlm_E, np, nr_cluster, nsol, nspin, nspino, nspinp, numat_abs, Space

  logical:: Hubb_a, Hubb_d, Hubb_m, Relativiste, Spinorbite, Stop_job, Ylm_comp

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nr_cluster,nspin,nlm_cp_E,nlm_cp_E):: V_cp
  complex(kind=db), dimension(nspino*nlm_E,nspino*nlm_E):: Mat_0, Mat_m, Mat_p
  complex(kind=db), dimension(nr_cluster,nlm_E,nlm_E,nspinp,nspino):: u
  complex(kind=db), dimension(nr_cluster,nspino):: gso

  real(kind=db):: a2s4, br, ci, Eimag, Enervide, fac, h2, hm, hp, p, Rmtg_abs
  real(kind=db), dimension(nspin):: cr, Ecinetic
  real(kind=db), dimension(nr_cluster):: f2, r_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_r_E,nlm_r_E):: V

  if( icheck > 1 ) write(3,110)

! alfa_sf = fine structure constant
  a2s4 = 0.25_db * alfa_sf**2

  Lmax = Lmaxso
  
  u(:,:,:,:,:) = ( 0._db, 0._db )

  call Coef_rad_Numerov(Eimag,Enervide,f2,gso,nlm_cp_E,nlm_r_E,nr_cluster,nspin,nspino,numat_abs,r_cluster,V,V_cp, &
                        Ylm_comp)

  n = 0
  do L = 0,Lmax

    L2 = L * ( L + 1 )

    ci = - Eimag / ( 4 * L + 6 )

    br = - numat_abs / ( L + 1._db )
    cr(:) = - ( 2 * numat_abs * br + Ecinetic(:) ) / ( 4 * L + 6 )

! Je fais L'hypothese que le terme de Hubbard ne croise pas les solutions 1 et 2 en cas de spinorbite.

! Development at origin
    do m = -L,L
      n = n + 1
      
      do isol = 1,nspinp ! spin at origin
        do isp = 1,nspinp ! real spin
          if( .not. Spinorbite .and. isp /= isol ) cycle
          isr = min( isp, nspin )

          if( Spinorbite ) then
            if( (m == L .and. isp == 1) .or. (m == -L .and. isp == 2 ) ) then
              nsol = 1
            else
              nsol = 2
            endif
            if( nsol == 1 .and. isp /= isol ) cycle

            if( Relativiste ) then
              if( nsol == 1 .or. isol == 1 ) then
                p = sqrt( ( L + 1 )**2 - ( alfa_sf * numat_abs )**2 )
              else
                p = sqrt( L**2 - ( alfa_sf * numat_abs )**2 )
              endif
            else
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
            endif
            
          else
            if( Relativiste ) then
              p = sqrt( L**2 + L + 1 - ( alfa_sf * numat_abs )**2 )
            else
              p = L + 1._db
            endif
          endif

 ! Becareful, it is from the point of view of the m of the orbital
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
          if( .not. Spinorbite .or. isp == isol ) then
            np = n
          elseif( isp < isol ) then
            np = n + 1
          else
            np = n - 1
          endif
          u(1:2,n,np,isp,isol) = fac * ( 1 + br * r_cluster(1:2) + cr(isr) * r_cluster(1:2)**2 ) * r_cluster(1:2)**p &
                               + img * fac * ci * r_cluster(1:2)**(p+2)
        end do
      end do
    end do
  end do

  hm = r_cluster(2) - r_cluster(1)
    
  do ir = 2,nr_cluster - 1
    hp = r_cluster(ir+1) - r_cluster(ir)
    if( Space == - 1 ) then
      Space = 2
      h2 = hp**2 / 12
      im = ir - 1
      i0 = ir
      ip = ir + 1
    elseif( hp > hm + eps10 ) then
      Space = 1
      h2 = hp**2 / 12
      im = ir - 2
      i0 = ir
      ip = ir + 1
    elseif(  hp < hm - eps10 ) then
      Space = - 1
      h2 = hm**2 / 12
      im = ir - 1
      i0 = ir
      ip = ir + 2
    else
      Space = 0
      h2 = hp**2 / 12
      im = ir - 1
      i0 = ir
      ip = ir + 1
    endif

    hm = hp ! for the next iteration on ir

    do isp = 1,nspinp
      isr = min( isp, nspin )
      isq = 3 - isp

      if( .not. ( Spinorbite .and. isp == nspinp ) ) then
        Mat_m(:,:) = ( 0._db, 0._db )
        Mat_0(:,:) = ( 0._db, 0._db )
        Mat_p(:,:) = ( 0._db, 0._db )
      endif
        
      do L = 0,Lmax

        L2 = L * ( L + 1 )

        if( Hubb_a .and. L == L_hubbard( numat_abs ) .and. r_cluster(i0) < Rmtg_abs )  then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif

        do m = -L,L

          Lm = L2 + 1 + m
          
          if( Spinorbite .and. isp == 2 ) then
            n = Lm + nlm_E
          else
            n = Lm
          endif
          
          Mat_m(n,n) = 1 + h2 * ( Enervide + img * Eimag - V(im,isr,Lm,Lm) - L2 * f2(im) )
          Mat_0(n,n) = 1 - 5 * h2 * ( Enervide + img * Eimag - V(i0,isr,Lm,Lm) - L2 * f2(i0) )
          Mat_p(n,n) = 1 + h2 * ( Enervide + img * Eimag - V(ip,isr,Lm,Lm) - L2 * f2(ip) )
          
          do Lmp = 1,nlm_E
            if( Lm == Lmp ) cycle
            if( Spinorbite .and. isp == 2 ) then
              np = Lmp + nlm_E
            else
              np = Lmp
            endif
            if( Ylm_comp ) then
              Mat_m(n,np) = - h2 * V_cp(im,isr,Lm,Lmp)
              Mat_0(n,np) = 5 * h2 * V_cp(i0,isr,Lm,Lmp)
              Mat_p(n,np) = - h2 * V_cp(ip,isr,Lm,Lmp)
            else
              Mat_m(n,np) = - h2 * V(im,isr,Lm,Lmp)
              Mat_0(n,np) = 5 * h2 * V(i0,isr,Lm,Lmp)
              Mat_p(n,np) = - h2 * V(ip,isr,Lm,Lmp)
            endif
          end do
          
          if( Hubb_m ) then
            do mp = -L,L
              if( Hubb_d .and. m /= mp ) cycle
              Lmp = L2 + 1 + mp
              if( Spinorbite .and. isp == 2 ) then
                np = Lmp + nlm_E
              else
                np = Lmp
              endif
              Mat_m(n,np) = Mat_m(n,np) - h2 * V_hubb(m,mp,isp,isp)
              Mat_0(n,np) = Mat_0(n,np) + 5 * h2 * V_hubb(m,mp,isp,isp)
              Mat_p(n,np) = Mat_p(n,np) - h2 * V_hubb(m,mp,isp,isp)
            end do
          endif
          
          if( Relativiste ) then
            Mat_m(n,n) = Mat_m(n,n) - h2 * a2s4 * ( Enervide + img * Eimag - V(im,isr,Lm,Lm) )**2
            Mat_0(n,n) = Mat_0(n,n) + 5 * h2 * a2s4 * ( Enervide + img * Eimag - V(i0,isr,Lm,Lm) )**2
            Mat_p(n,n) = Mat_p(n,n) - h2 * a2s4 * ( Enervide + img * Eimag - V(ip,isr,Lm,Lm) )**2

            Mat_m(n,n) = Mat_m(n,n) - h2 * gso(im,isr)
            Mat_0(n,n) = Mat_0(n,n) + 5 * h2 * gso(i0,isr)
            Mat_p(n,n) = Mat_p(n,n) - h2 * gso(ip,isr)
          endif
          
          if( Spinorbite ) then
            if( isp == 1 ) then
              ms = m
              mp = m + 1   ! m of the other spin
            else
              ms = - m
              mp = m - 1
            endif
 
            Mat_m(n,n) = Mat_m(n,n) - h2 * ms * gso(im,isp)
            Mat_0(n,n) = Mat_0(n,n) + 5 * h2 * ms * gso(i0,isp)
            Mat_p(n,n) = Mat_p(n,n) - h2 * ms * gso(ip,isp)

            if( abs(mp) <= L ) then
              Lmp = L2 + 1 + mp
              if( Spinorbite .and. isq == 2 ) then
                np = Lmp + nlm_E
              else
                np = Lmp
              endif
              fac = sqrt( ( L - ms ) * ( L + ms + 1._db ) )
              Mat_m(n,np) = Mat_m(n,np) - h2 * fac * gso(im,isp)
              Mat_0(n,np) = Mat_0(n,np) + 5 * h2 * fac * gso(i0,isp)
              Mat_p(n,np) = Mat_p(n,np) - h2 * fac * gso(ip,isp)
            endif

          endif
          
        end do ! Loop on m
      end do ! Loop on L

      if( Spinorbite ) cycle

      if( icheck > 3 .and. ( ir < 4 .or. ir > nr_cluster - 3 ) ) then
        if( Space == 2 ) then
          write(3,'(/a30,i6)') '  Mat_0 before inversion, ir =', ir
        else
          write(3,'(/a13,i6)') '  Mat_0, ir =', ir
        endif   
        do n = 1,nlm_E
          write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_0(n,:)   
        end do
        write(3,'(/a13,i6)') '  Mat_m, ir =', ir   
        do n = 1,nlm_E
          write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_m(n,:)   
        end do
        if( Space == 2 ) then
          write(3,'(/a13,i6)') '  Mat_p, ir =', ir
        else
          write(3,'(/a30,i6)') '  Mat_p before inversion, ir =', ir
        endif   
        do n = 1,nlm_E
          write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_p(n,:)   
        end do
      endif

      if( Space == 2 ) then
        call invcomp(nlm_E,Mat_0,nlm_E,nlm_E,1,Stop_job)
      else
        call invcomp(nlm_E,Mat_p,nlm_E,nlm_E,1,Stop_job)
      endif

      if( icheck > 3 .and. ( ir < 4 .or. ir > nr_cluster - 3 ) ) then
        if( Space == 2 ) then
          write(3,'(/a29,i6)') '  Mat_0 after inversion, ir =', ir
          do n = 1,nlm_E
            write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_0(n,:)   
          end do
        else
          write(3,'(/a29,i6)') '  Mat_p after inversion, ir =', ir
          do n = 1,nlm_E
            write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_p(n,:)   
          end do
        endif   
      endif
      
      if( Space == 2 ) then
      
        Mat_p = Matmul( Mat_0, Mat_p ) 
        Mat_m = Matmul( Mat_0, Mat_m ) 

        do Lm = 1,nlm_E
          do Lmp = 1,nlm_E
            u(i0,:,Lm,isp,:) = u(i0,:,Lm,isp,:) + 0.5_db * ( Mat_p(Lm,Lmp) * u(ip,:,Lmp,isp,:) + Mat_m(Lm,Lmp) * u(im,:,Lmp,isp,:) ) 
          end do
        end do
        
      else
      
        Mat_0 = Matmul( Mat_p, Mat_0 ) 
        Mat_m = Matmul( Mat_p, Mat_m ) 

        do Lm = 1,nlm_E
          do Lmp = 1,nlm_E
            u(ip,:,Lm,isp,:) = u(ip,:,Lm,isp,:) + 2 * Mat_0(Lm,Lmp) * u(i0,:,Lmp,isp,:) - Mat_m(Lm,Lmp) * u(im,:,Lmp,isp,:) 
          end do
        end do
        
      endif
      
    end do ! Loop on isp

    if( .not. Spinorbite ) cycle
    
    call invcomp(2*nlm_E,Mat_p,2*nlm_E,2*nlm_E,0,Stop_job)
    
    Mat_0 = 2 * Matmul( Mat_p, Mat_0 ) 
    Mat_m = - Matmul( Mat_p, Mat_m ) 

    do isp = 1,nspino 
      do Lm = 1,nlm_E
        n = Lm + ( isp - 1 ) * nlm_E
        do isq = 1,nspino
          do Lmp = 1,nlm_E
            np = Lmp + ( isq - 1 ) * nlm_E
            u(ip,:,Lm,isp,:) = u(ip,:,Lm,isp,:) + Mat_0(n,np) * u(i0,:,lmp,isq,:) + Mat_m(n,np) * u(im,:,lmp,isq,:) 
          end do
        end do
      end do
    end do

  end do ! Loop on ir

  if( icheck > 1 ) call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,1)

  call Renormal_cluster(icheck,konde,Lmax,nlm_E,nr_cluster,nspin,nspino,nspinp,r_cluster,Spinorbite,u)

  if( icheck > 1 ) call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,2)

  return
  110 format(/'-------- Cluster regular solutions ',80('-'))
end

!***********************************************************************

subroutine Coef_rad_Numerov(Eimag,Enervide,f2,gso,nlm_cp_E,nlm_r_E,nr,nspin,nspino,numat,r,V,V_cp,Ylm_comp)

  use declarations
  implicit none

  integer:: ir, isp, nlm_r_E, nlm_cp_E, nr, nspin, nspino, numat

  complex(kind=db):: bder, dvr, VmE
  complex(kind=db), dimension(nr,nspino):: gso
  complex(kind=db), dimension(nr,nspin,nlm_cp_E,nlm_cp_E):: V_cp

  logical:: Ylm_comp

  real(kind=db):: a2s4, dr, Eimag, Enervide, r0, rm, rp

  real(kind=db), dimension(nr):: cgrad0, cgradm, cgradp, f2, r
  real(kind=db), dimension(nr,nspin,nlm_r_E,nlm_r_E):: V
  
  gso(:,:) = ( 0._db, 0._db )

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
    cgradm(ir) = ( rp - r0 ) / ( ( rm - r0 ) * ( rp - rm ) )
    cgradp(ir) = ( rm - r0 ) / ( ( rp - r0 ) * ( rm - rp ) )
    cgrad0(ir) = - cgradm(ir) - cgradp(ir)
    f2(ir) = 1 / r(ir)**2
  end do

! alfa_sf = fine structure constant
  a2s4 = 0.25_db * alfa_sf**2

  do ir = 1,nr-1
    do isp = 1,nspin

      if( Ylm_comp ) then
        Vme = V_cp(ir,isp,1,1) - Enervide - img * Eimag
      else
        Vme = V(ir,isp,1,1) - Enervide - img * Eimag
      endif

      bder = 1 / ( 1 - a2s4 * Vme )
      if( ir == 1 ) then
        dvr = cmplx( 2 * numat / r(1)**2, 0._db, db )
      else
        if( Ylm_comp ) then
          dvr = cgradm(ir) * V_cp(ir-1,isp,1,1) + cgrad0(ir) * V_cp(ir,isp,1,1) + cgradp(ir) * V_cp(ir+1,isp,1,1)
        else
          dvr = cgradm(ir) * V(ir-1,isp,1,1) + cgrad0(ir) * V(ir,isp,1,1) + cgradp(ir) * V(ir+1,isp,1,1)
        endif
      endif

      gso(ir,isp) = a2s4 * bder * dvr / r(ir)

    end do
  end do

  if( nspin < nspino ) gso(:,nspino) = gso(:,nspin)

  return
end

!***********************************************************************

! Resolution of the radial Schrodinger equation. Include the case of the non spherical potential
! Uses the matrix way with Gaussian elimination
! In ur(ir,m,mp,isp,isol) : (mp, isol) define the basis state
!                           (m, isp) are the real moment and spin

subroutine Sch_radial_reg_irreg_gaussian(Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d,icheck,Irregular,konde,Lmso,Lmaxso,m_hubb, &
         ngrph,nlm_cp_E,nlm_r_E,nlm_E,nlmso0,nr_cluster,nso1,nspin,nspino,nspinp,numat_abs,Numerov,r_cluster, &
         Relativiste,Rmtg_abs,Spinorbite,u,V,V_cp,V_hubb,Ylm_comp)

  use declarations
  implicit none

  integer:: i, i_sm, i_smp, ib, ib_min, icheck, igrph, im, ip, ir, iso, isop, isp, isp_no_so, isp_so, ispp, ispp_so, isr, j, &
    j_Lmp, j_max, jr, k, L, L_hubbard, LLp1, Lm, Lmax, Lmaxso, Lmp, Lms, Lmsp, Lmss, Lmssp, Lp, m, m_hubb, mp, n_band, n_line, &
    nb_m, nb_p, ngrph, nlm, nlm_cp_E, nlm_r_E, nlm_E, nlms, nr, nr_cluster, nso1, nspin, nspino, nspinp, numat_abs

  integer, dimension(ngrph):: nlmso0
  integer, dimension(nso1,3,ngrph):: Lmso
  integer, dimension(:), allocatable:: Band_p

  logical:: Hubb_a, Hubb_d, Hubb_m, Irregular, Numerov, Relativiste, Spinorbite, Ylm_comp

  complex(kind=db):: Coef_m, Coef_p, Diag, EmV_0, EmV_m, EmV_p, fnormc, Line_ib, z
  complex(kind=db), dimension(nspin):: cr, konde
  complex(kind=db), dimension(0:Lmaxso):: bess, neum
  complex(kind=db), dimension(0:Lmaxso,nspin):: Rap
  complex(kind=db), dimension(nr_cluster-1:nr_cluster,0:Lmaxso,nspin):: r_Bessel, r_Hankel
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nr_cluster,nspin,nlm_cp_E,nlm_cp_E):: V_cp
  complex(kind=db), dimension(nr_cluster,nlm_E,nlm_E,nspinp,nspino):: u

  complex(kind=db), dimension(:), allocatable:: Line
  complex(kind=db), dimension(:,:), allocatable:: Mat, Sm

  real(kind=db):: br, dr, h2, hm, hp, r0, rm, rp, Eimag, Enervide, Rmtg_abs
  real(kind=db), dimension(nspin):: Ecinetic
  real(kind=db), dimension(nr_cluster):: clapl0, claplm, claplp, f2, r_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_r_E,nlm_r_E):: V

  if( icheck > 1 ) then
    if( Irregular ) then
      Write(3,110)
    else
      Write(3,115)
    endif
  endif
  
  Lmax = Lmaxso
  
  if( Irregular ) then
    nr = nr_cluster - 2
  else
    nr = nr_cluster - 1
  endif
  
  u(:,:,:,:,:) = ( 0._db, 0._db )

  f2(:) = 1 / r_cluster(:)**2

    do ir = 1,nr
      rp = r_cluster(ir+1)
      if( ir == 1 ) then
    ! ir = 1 is not used
        rm = r_cluster(1)**2 / r_cluster(2)
      else
        rm = r_cluster(ir-1)
      endif
      r0 = r_cluster(ir)
      dr = 0.5_db * ( rp - rm )
      claplm(ir) = 1 / ( ( r0 - rm ) * dr )
      claplp(ir) = 1 / ( ( rp - r0 ) * dr )
      clapl0(ir) = - claplm(ir) - claplp(ir)
    end do
  
  do isp = 1,nspin
    fnormc = sqrt( konde(isp) / pi )
    do ir = nr_cluster-1,nr_cluster
      z = konde(isp) * r_cluster(ir)
      call cbessneu(fnormc,z,Lmax,Lmax,bess,neum)
      r_Bessel(ir,0:Lmax,isp) = r_cluster(ir) * bess(0:Lmax) 
      r_Hankel(ir,0:Lmax,isp) = - img * r_cluster(ir) * ( bess(0:Lmax) + img * neum(0:Lmax) ) 
    end do
  end do
  
  if( Irregular ) then
    call Isotrope_radial_irreg(clapl0,claplm,claplp,Eimag,Enervide,f2,icheck,Lmaxso,nlm_cp_E,nlm_r_E,nr_cluster, &
                                 nspin,r_cluster,r_Hankel,Rap,V,V_cp,Ylm_comp)
  else
    do L = 0,Lmax
      cr(:) = - ( 2 * numat_abs * br + Ecinetic(:) + img * Eimag  ) / ( 4 * L + 6 )
      Rap(L,:) = ( 1 + br * r_cluster(1) + cr(:) * r_cluster(1)**2 ) * r_cluster(1)**(L+1) / &
               ( ( 1 + br * r_cluster(2) + cr(:) * r_cluster(2)**2 ) * r_cluster(2)**(L+1) )
    end do
  endif

  Loop_igrph: do igrph = 1,ngrph

    do Lms = nlmso0(igrph),1,-1
      L = Lmso(Lms,1,igrph)
      if( L <= Lmax ) exit
    end do
    nlm = Lms
    if( nlm == 0 ) cycle

    nlms = nlm * nspino 
    n_line = ( nr - 1 ) * nlms
    
    if( Numerov ) then
   !   factor 3 and not 2 because of points where hp = 2*hm
      n_band = 3 * nlms - 1
    else
      n_band = nlms
    endif

    allocate( Band_p(2:nr) )
    allocate( Mat(n_line,n_band) )
    allocate( Sm(nlms,nlms) )
    allocate( Line(n_line) )

! Loop on spin used only when no spinorbite
    Loop_isp_no_so: do isp_no_so = 1,nspin   

      if( .not. Spinorbite ) isp = isp_no_so
      
      Mat(:,:) = ( 0._db, 0._db )
      Sm(:,:) = ( 0._db, 0._db )

      if( Numerov ) hm = r_cluster(2) - r_cluster(1)
           
      ir = 2
      i = 0

      do ir = 2,nr

        if( Numerov ) then
          hp = r_cluster(ir+1) - r_cluster(ir)
          if( hp > hm + eps10 ) then
            h2 = hp**2 / 12
            im = ir - 2
            ip = ir + 1
            nb_m = 2 * nlms
            nb_p = nlms
          elseif(  hp < hm - eps10 ) then
            h2 = hm**2 / 12
            im = ir - 1
            ip = ir + 2
            nb_m = nlms
            nb_p = 2 * nlms
          else
            h2 = hp**2 / 12
            im = ir - 1
            ip = ir + 1
            nb_m = nlms
            nb_p = nlms
          endif
 
          hm = hp ! for the next iteration on ir
          
        else
        
          nb_m = nlms
          nb_p = nlms
          
        endif

        Band_p(ir) = nb_p
          
        do isp_so = 1,nspino

          if( Spinorbite ) isp = isp_so
          isr = min( isp, nspin )

          do Lms = 1,nlm
            L = Lmso(Lms,1,igrph)
            m = Lmso(Lms,2,igrph)
            
            LLp1 = L * ( L + 1 )
  
            Lm = LLp1 + 1 + m
            i = i + 1
  
! Line filling 
            Line(:) = ( 0._db, 0._db )

            if( Ylm_comp ) then
              EmV_0 = Enervide + img * Eimag - V_cp(ir,isr,Lm,Lm) - LLp1 * f2(ir)
            else
              EmV_0 = cmplx( Enervide - V(ir,isr,Lm,Lm) - LLp1 * f2(ir), Eimag, db )
            endif

            if( Numerov ) then
              if( Ylm_comp ) then
                EmV_m = Enervide + img * Eimag - V_cp(im,isr,Lm,Lm) - LLp1 * f2(im) 
                EmV_p = Enervide + img * Eimag - V_cp(ip,isr,Lm,Lm) - LLp1 * f2(ip)
              else
                EmV_m = cmplx( Enervide - V(im,isr,Lm,Lm) - LLp1 * f2(im), Eimag, db )   
                EmV_p = cmplx( Enervide - V(ip,isr,Lm,Lm) - LLp1 * f2(ip), Eimag, db )
              endif
            endif

            if( Numerov ) then
              Diag = 2._db - 10 * h2 * EmV_0 
            else
              Diag = - clapl0(ir) - EmV_0 
            endif

            if( Numerov ) then
              Coef_m = - ( 1._db + h2 * EmV_m ) 
              Coef_p = - ( 1._db + h2 * EmV_p ) 
            else  
              Coef_m = cmplx( - claplm(ir), 0._db, db )
              Coef_p = cmplx( - claplp(ir), 0._db, db )
            endif

            if( ir > 2 )  Line( i - nb_m ) = Coef_m 
            if( ir < nr ) Line( i + nb_p ) = Coef_p

            if( ir == 2 ) then

              Diag = Diag + Coef_m * Rap(L,isr)

            elseif( ir == nr ) then

              if( Spinorbite ) then
                iso = min( Lmso(Lmsp,3,igrph), nspin )
              else
                iso = isp
              endif
              if( .not. Irregular ) Diag = Diag + Coef_p * r_Hankel(nr+1,L,iso) / r_Hankel(nr,L,iso)

! Second member filling
              if( iso == isp .or. .not. Spinorbite ) then
                if( isp_so == 2 ) then
                  i_sm = Lms + nlm
                else
                  i_sm = Lms
                endif
                if( Irregular ) then
                  Sm(i_sm,i_sm) = - Coef_p * r_Hankel(nr+1,L,iso)
                else 
                  Sm(i_sm,i_sm) = - Coef_p &
                                * ( r_Bessel(nr+1,L,iso) - r_Hankel(nr+1,L,iso) * r_Bessel(nr,L,iso) / r_Hankel(nr,L,iso) )
                endif
              endif
 
            endif

            Line(i) = Diag
                     
! Non diagonal potential components
! Potential taken diagonal at the first and last radius points
            if( ir > 2 .and. ir < nr_cluster-1 ) then         
  
              do Lmsp = 1,nlm
                Lp = Lmso(Lmsp,1,igrph)
                mp = Lmso(Lmsp,2,igrph)
                Lmp = Lp**2 + Lp + 1 + mp
            
                if( Lmsp == Lms ) cycle

                j_Lmp = i + Lmsp - Lms 
                
                if( Numerov ) then
  
                  do k = -1,1
                    select case(k)
                      case(-1)
                        if( Ylm_comp ) then
                          Line(j_Lmp - nb_m) = h2 * V_cp(im,isr,Lm,Lmp) 
                        else
                          Line(j_Lmp - nb_m) = h2 * V(im,isr,Lm,Lmp)
                        endif
                      case(0)
                        if( Ylm_comp ) then
                          Line(j_Lmp) = 10 * h2 * V_cp(ir,isr,Lm,Lmp) 
                        else
                          Line(j_Lmp) = 10 * h2 * V(ir,isr,Lm,Lmp)
                        endif
                      case(1)
                        if( Ylm_comp ) then
                          Line(j_Lmp + nb_p) = h2 * V_cp(ip,isr,Lm,Lmp) 
                        else
                          Line(j_Lmp + nb_p) = h2 * V(ip,isr,Lm,Lmp)
                        endif
                    end select
                 ! at nr_cluster - 1, potential is isotrope (and corresponding matrix line does not exist for irregular solution
                    if( ir == nr_cluster-2 .and. k == 0 ) exit 
                  end do
  
                else
  
                  if( Ylm_comp ) then
                    Line(j_Lmp) = V_cp(ir,isr,Lm,Lmp) 
                  else
                    Line(j_Lmp) = V(ir,isr,Lm,Lmp)
                  endif
                   
                endif
  
              end do
              
            endif
  
            if( icheck > 3 ) then
              if( ir == 2 ) then
                if( i == 1 )  write(3,'(/a16,i2,/A)') ' Matrix, igrph =', igrph, '   ir  Lm isp'
                write(3,120) ir, Lm, isp, Line(i), Line(i+nlms)
              elseif( ir == nr ) then
                if( Numerov ) then
                  write(3,130) ir, Lm, isp, Line(i-Lms+1+nb_m:n_line), Sm(i_sm,i_sm) 
                else
                  write(3,130) ir, Lm, isp, Line(i-nlms), Line(i+1-Lms:i+nlm-Lms), Sm(i_sm,i_sm) 
                endif
              else
                if( Numerov ) then
                  write(3,130) ir, Lm, isp, Line( i-Lms+1-nb_m:i+nlms-Lms+nb_p) 
                else
                  write(3,130) ir, Lm, isp, Line(i-nlms), Line(i+1-Lms:i+nlm-Lms), Line(i+nlms) 
                endif
              endif
            endif

! Triangularisation is immediatly done on the line
            if( ir == 2 ) then
              Mat(i,nlms) =  Line(i+nlms) / Line(i)   
              if( icheck > 3 ) write(3,'(27x,a9,1p,1000(1x,2e13.5))') '  Triang:', Mat(i,nlms)
              cycle
            endif

            if( Numerov ) then
              ib_min = i - Lms + 1 - nb_m
            else
              ib_min = i - nlms
            endif  
            do ib = ib_min,i-1 
              if( Numerov ) then
                if( ib <= i - Lms - nlms ) then
                  Lmp = ib - i + Lms + 2 * nlms
                  j_max = min( ib + nlms - Lmp + Band_p(ir-2), n_line )
                elseif( ib <= i - Lms ) then
                  Lmp = ib - i + Lms + nlms
                  j_max = min( ib + nlms - Lmp + Band_p(ir-1), n_line )
                else
                  Lmp = ib - i + Lms 
                  j_max = min( ib + nlms - Lmp + Band_p(ir), n_line )
                endif
              else
                j_max = min( ib + nlms, n_line )
              endif  
              Line_ib = Line(ib)  ! not divided by diagonal matrix term, because = 1 (and not not storaged)   
              do k = ib+1,j_max   ! starting at ib+1 and not ib because Line(ib) is not used anymore (and corresponding matrix term not storaged)
                j = k - ib
                Line(k) = Line(k) - Mat(ib,j) * Line_ib
              end do
     
              i_smp = ib - ( n_line - nlms )          
              if( i_smp > 0 ) Sm(i_sm,:) = Sm(i_sm,:) - Sm(i_smp,:) * Line_ib 
            end do
        
 ! normalisation by the diagonal matrix term, which is not storaged
            if( Numerov ) then
              j_max = min( nlms - Lms + Band_p(ir), n_line - i )
            else
              j_max = min( nlms, n_line - i )
            endif  
            do j = 1,j_max
              Mat(i,j) = Line(i+j) / Line(i)
            end do
            if( ir == nr ) Sm(Lms,:) = Sm(Lms,:) / Line(i)
          
            if( icheck > 3 ) then
              if( Numerov ) then
                write(3,'(27x,a9,1p,1000(1x,2e13.5))') '  Triang:', Mat(i,1:nlms-Lms+Band_p(ir))
              else
                write(3,'(27x,a9,1p,1000(1x,2e13.5))') '  Triang:', Mat(i,:)
              endif
     !!!!!!!
 !    if( igrph == 2 ) write(4,'(f11.7,2i4,1p,2e13.5)') r_cluster(ir)*bohr, nb_m, nb_p, Mat(i,1) 
              if( ir == nr ) write(3,'(27x,a9,1p,1000(1x,2e13.5))') '      Sm:', Sm(Lms,:)
            endif
            
          end do ! end of Loop on Lms
      
        end do ! end of Loop on isp_so
        
      end do ! end of Loop on ir
  
 ! Resolution
      i = n_line + 1 
      do ir = nr,2,-1

        Lmss = nlms + 1       
        do isp_so = nspino,1,-1     
          do Lms = nlm,1,-1
            Lmss = Lmss - 1
            if( Spinorbite ) isp = isp_so
              
            L = Lmso(Lms,1,igrph)
            m = Lmso(Lms,2,igrph)
            iso = Lmso(Lms,3,igrph)
            Lm = L**2 + L + 1 + m
        
            i = i - 1
            if( ir == nr ) then
              Lmssp = 0 
              do ispp_so = nspino,1,-1     
                do Lmsp = 1,nlm
                  Lmssp = Lmssp + 1
                  if( Spinorbite ) then
                    ispp = isp_so
                    isop = Lmso(Lmsp,3,igrph)
                  else
                    ispp = isp
                    isop = 1
                  endif
                  Lmp = Lmso(Lmsp,1,igrph)**2 + Lmso(Lmsp,1,igrph) + 1 + Lmso(Lmsp,2,igrph)
                  u(ir,Lmp,Lm,isp,iso) = Sm(Lmss,Lmssp)
                end do 
              end do 
            endif
        
            if( Numerov ) then
              j_max = i - Lms + nlms + Band_p(ir)
            else
              j_max = nlms
            endif
            
            do j = 1,j_max
              if( ir == nr .and. i + j > n_line ) exit
              Lmsp = Lms + j
              if( Lmsp > nlms ) then
                jr = ir + 1
                Lmsp = Lmsp - nlms
                if( Lmsp > nlms ) then
                  if( Band_p(ir) == nlms ) exit
                  jr = jr + 1
                  Lmsp = Lmsp - nlms
                  if( Lmsp > nlms ) exit 
                endif
                if( jr > nr ) exit
              else
                jr = ir
              endif
              if( Spinorbite ) then
                if( Lmsp > nlm ) then
                  ispp = 2
                  Lmsp = Lmsp - nlm
                else
                  ispp = 1
                endif
              else
                ispp = isp
              endif
                
              Lp = Lmso(Lmsp,1,igrph)
              mp = Lmso(Lmsp,2,igrph)
              isop = Lmso(Lmsp,3,igrph)
              Lmp = Lp**2 + Lp + 1 + mp
              u(ir,:,Lm,isp,iso) = u(ir,:,Lm,isp,iso) - Mat(i,j) * u(jr,:,Lmp,ispp,iso)  
            end do
          end do   
        end do ! Loop on Lms
      end do
      
    end do Loop_isp_no_so

    deallocate( Band_p, Line, Mat, Sm )
 
    if( Spinorbite ) exit
    
  end do Loop_igrph ! end of loop on representations

! Filling of the ir = 1 point using the ratio calculted previously

  do isp = 1,nspinp
    isr = min( isp, nspin )
    Lm = 0
    do L = 0,Lmax
      do m = -L,L
        Lm = Lm + 1
        u(1,:,Lm,isp,:) =  Rap(L,isr) * u(2,:,Lm,isp,:) 
      end do
    end do
  end do
  
  if( Irregular ) then
    do isp = 1,nspinp
      isr = min( isp, nspin )
      Lm = 0
      do L = 0,Lmax
        do m = -L,L
          Lm = Lm + 1
          u(nr_cluster-1,Lm,Lm,isp,:) =  r_Hankel(nr_cluster-1,L,isp) 
          u(nr_cluster,Lm,Lm,isp,:) =  r_Hankel(nr_cluster,L,isp) 
        end do
      end do
    end do
    
  else
  
    call Cal_tauLL(icheck,Lmaxso,r_Bessel,r_Hankel,nlm_E,nr_cluster,nspin,nspino,nspinp,Spinorbite,u)

  endif

  if( icheck > 1 ) then
    if( irregular ) then
      call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,3)
    else
      call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,4)
    endif
  endif

  return
  110 format(/'-------- Cluster irregular solutions ',80('-'))
  115 format(/'-------- Cluster regular solutions ',80('-'))
  120 format(i5,2i4,23x,1p,10000(1x,2e13.5))
  130 format(i5,2i4,1p,10000(1x,2e13.5))
end

!***********************************************************************

! Resolution of the radial Schrodinger equation. Include the case of the non spherical potential
! Uses the matrix way with Gaussian elimination
! In ur(ir,m,mp,isp,isol) : (mp, isol) define the basis state
!                           (m, isp) are the real moment and spin

subroutine Sch_radial_reg_irreg_gaussian_old(Ecinetic,Eimag,Enervide,Hubb_a,Hubb_d,icheck,Irregular,konde,Lmso,Lmaxso,m_hubb, &
         ngrph,nlm_cp_E,nlm_r_E,nlm_E,nlmso0,nr_cluster,nso1,nspin,nspino,nspinp,numat_abs,r_cluster, &
         Relativiste,Rmtg_abs,Spinorbite,u,V,V_cp,V_hubb,Ylm_comp)

  use declarations
  implicit none

  integer:: i, i_sm, i_smp, icheck, igrph, ir, iso, isop, isp, isp_no_so, isp_so, ispp, isr, j, j2, jr, k, L, L_hubbard, Lm, &
    Lmax, Lmaxso, Lmp, Lms, Lmsp, Lmss, Lmssp, Lp, m, m_hubb, mp, n_line, ngrph, nlm, nlm_cp_E, nlm_r_E, nlm_E, nlms, nr, &
    nr_cluster, nso1, nspin, nspino, nspinp, numat_abs

  integer, dimension(ngrph):: nlmso0
  integer, dimension(nso1,3,ngrph):: Lmso

  logical:: Hubb_a, Hubb_d, Hubb_m, Irregular, Relativiste, Spinorbite, Ylm_comp

  complex(kind=db):: Diag, fnormc, Line_j2, z
  complex(kind=db), dimension(nspin):: cr, konde
  complex(kind=db), dimension(0:Lmaxso):: bess, neum
  complex(kind=db), dimension(0:Lmaxso,nspin):: Rap
  complex(kind=db), dimension(nr_cluster-1:nr_cluster,0:Lmaxso,nspin):: r_Bessel, r_Hankel
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nr_cluster,nspin,nlm_cp_E,nlm_cp_E):: V_cp
  complex(kind=db), dimension(nr_cluster,nlm_E,nlm_E,nspinp,nspino):: u

  complex(kind=db), dimension(:), allocatable:: Line
  complex(kind=db), dimension(:,:), allocatable:: Mat, Sm

  real(kind=db):: br, dr, E_centrifuge, r0, rm, rp, Eimag, Enervide, Rmtg_abs
  real(kind=db), dimension(nspin):: Ecinetic
  real(kind=db), dimension(nr_cluster):: clapl0, claplm, claplp, f2, r_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_r_E,nlm_r_E):: V

  if( icheck > 1 ) then
    if( Irregular ) then
      Write(3,110)
    else
      Write(3,115)
    endif
  endif
  
  Lmax = Lmaxso
  
  if( Irregular ) then
    nr = nr_cluster - 2
  else
    nr = nr_cluster - 1
  endif
  
  u(:,:,:,:,:) = ( 0._db, 0._db )

  do ir = 1,nr
    rp = r_cluster(ir+1)
    if( ir == 1 ) then
  ! ir = 1 is not used
      rm = r_cluster(1)**2 / r_cluster(2)
    else
      rm = r_cluster(ir-1)
    endif
    r0 = r_cluster(ir)
    dr = 0.5_db * ( rp - rm )
    claplm(ir) = 1 / ( ( r0 - rm ) * dr )
    claplp(ir) = 1 / ( ( rp - r0 ) * dr )
    clapl0(ir) = - claplm(ir) - claplp(ir)
    f2(ir) = 1 / r_cluster(ir)**2
  end do
  
  do isp = 1,nspin
    fnormc = sqrt( konde(isp) / pi )
    do ir = nr_cluster-1,nr_cluster
      z = konde(isp) * r_cluster(ir)
      call cbessneu(fnormc,z,Lmax,Lmax,bess,neum)
      r_Bessel(ir,0:Lmax,isp) = r_cluster(ir) * bess(0:Lmax) 
      r_Hankel(ir,0:Lmax,isp) = - img * r_cluster(ir) * ( bess(0:Lmax) + img * neum(0:Lmax) ) 
    end do
  end do
  
  if( Irregular ) then
    call Isotrope_radial_irreg(clapl0,claplm,claplp,Eimag,Enervide,f2,icheck,Lmaxso,nlm_cp_E,nlm_r_E,nr_cluster, &
                                 nspin,r_cluster,r_Hankel,Rap,V,V_cp,Ylm_comp)
  else
    do L = 0,Lmax
      cr(:) = - ( 2 * numat_abs * br + Ecinetic(:) + img * Eimag  ) / ( 4 * L + 6 )
      Rap(L,:) = ( 1 + br * r_cluster(1) + cr(:) * r_cluster(1)**2 ) * r_cluster(1)**(L+1) / &
               ( ( 1 + br * r_cluster(2) + cr(:) * r_cluster(2)**2 ) * r_cluster(2)**(L+1) )
    end do
  endif

  Loop_igrph: do igrph = 1,ngrph

    do Lms = nlmso0(igrph),1,-1
      L = Lmso(Lms,1,igrph)
      if( L <= Lmax ) exit
    end do
    nlm = Lms
    if( nlm == 0 ) cycle

    nlms = nlm * nspino 
    n_line = ( nr - 1 ) * nlms

    allocate( Mat(n_line,nlms) )
    allocate( Sm(nlms,nlms) )
    allocate( Line(n_line) )

! Loop on spin used only when no spinorbite
    Loop_isp_no_so: do isp_no_so = 1,nspin   

      if( .not. Spinorbite ) isp = isp_no_so
      
      Mat(:,:) = ( 0._db, 0._db )
      Sm(:,:) = ( 0._db, 0._db )
    
      ir = 2
      i = 0

! Potential taken diagonal at the first radius point
      do isp_so = 1,nspino

        if( Spinorbite ) isp = isp_so
        isr = min( isp, nspin )
    
        do Lms = 1,nlm
          L = Lmso(Lms,1,igrph)
          m = Lmso(Lms,2,igrph)
          Lm = L**2 + L + 1 + m
          
          E_centrifuge = L * ( L + 1 ) * f2(ir)
  
          i = i + 1
  
          if( Ylm_comp ) then
            Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V_cp(ir,isr,Lm,Lm) 
          else
            Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V(ir,isr,Lm,Lm) 
          endif
          Diag = Diag - claplm(ir) * Rap(L,isr)
          Mat(i,nlms) = - claplp(ir) / Diag   
  
          if( icheck > 3 ) then
            if( i == 1 )  write(3,'(/A/A)') ' Matrix','   ir  Lm isp'       
            write(3,120) ir, Lm, isp, Diag, - claplp(ir)
          endif
          
        end do 
  
      end do
  
      do ir = 3,nr
  
        do isp_so = 1,nspino

          if( Spinorbite ) isp = isp_so
          isr = min( isp, nspin )

          do Lms = 1,nlm
            L = Lmso(Lms,1,igrph)
            m = Lmso(Lms,2,igrph)
  
            E_centrifuge = L * ( L + 1 ) * f2(ir)
  
            Lm = L**2 + L + 1 + m
            i = i + 1
  
! Line filling 
            Line(:) = ( 0._db, 0._db )
          
            if( Ylm_comp ) then
              Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V_cp(ir,isr,Lm,Lm) 
            else
              Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V(ir,isr,Lm,Lm) 
            endif
  
            Line(i-nlms) = - claplm(ir)
          
            do Lmsp = 1,nlm
              Lp = Lmso(Lmsp,1,igrph)
              mp = Lmso(Lmsp,2,igrph)
              Lmp = Lp**2 + Lp + 1 + mp
            
              if( Lmp == Lm ) cycle
              j = i + Lmsp - Lms
              if( Ylm_comp ) then
                Line(j) = V_cp(ir,isr,Lm,Lmp) 
              else
                Line(j) = V(ir,isr,Lm,Lmp) 
              endif
            end do
  
            if( ir == nr ) then
              if( Spinorbite ) then
                iso = min( Lmso(Lmsp,3,igrph), nspin )
              else
                iso = isp
              endif
              if( .not. Irregular ) Diag = Diag - claplp(ir) * r_Hankel(nr+1,L,iso) / r_Hankel(nr,L,iso)
            else
              Line(i+nlms) = - claplp(ir)
            endif
  
            Line(i) = Diag
                   
! Second member filling
            if( ir == nr .and. ( iso == isp .or. .not. Spinorbite ) ) then
              if( isp_so == 2 ) then
                i_sm = Lms + nlm
              else
                i_sm = Lms
              endif
              if( Irregular ) then
                Sm(i_sm,i_sm) = claplp(ir) * r_Hankel(nr+1,L,isp)
              else 
                Sm(i_sm,i_sm) = claplp(ir) &
                                       * ( r_Bessel(nr+1,L,isp) - r_Hankel(nr+1,L,isp) * r_Bessel(nr,L,isp) / r_Hankel(nr,L,isp) )
              endif
            endif
   
            if( icheck > 3 ) then
              if( ir == nr ) then
                write(3,130) ir, Lm, isp, Line(i-nlms), Line(i+1-Lms:i+nlm-Lms), Sm(i_sm,i_sm) 
              else
                write(3,130) ir, Lm, isp, Line(i-nlms), Line(i+1-Lms:i+nlm-Lms), Line(i+nlms) 
              endif
            endif

! Triangularisation is immediatly done on the line
            do j2 = i-nlms,i-1 
              Line_j2 = Line(j2)  ! not divided by diagonal matrix term, because = 1 (and not not storaged)   
              do k = j2+1,min(j2+nlms,n_line)   ! starting at j2+1 and not j2 because Line(j2) is not used anymore (and corresponding matrix term not storaged)
                j = k - j2
                Line(k) = Line(k) - Mat(j2,j) * Line_j2
              end do
     
              i_smp = j2 - ( n_line - nlms )          
              if( i_smp > 0 ) Sm(i_sm,:) = Sm(i_sm,:) - Sm(i_smp,:) * Line_j2 
            end do
        
 ! normalisation by the diagonal matrix term, which is not storaged
            do j = 1,nlms
              if( i + j > n_line ) exit 
              Mat(i,j) = Line(i+j) / Line(i)
            end do
            if( ir == nr ) Sm(Lms,:) = Sm(Lms,:) / Line(i)
          
            if( icheck > 3 ) write(3,'(23x,a9,1p,1000(1x,2e11.3))') '  Triang:', Mat(i,:)
            if( ir == nr .and. icheck > 3 ) write(3,'(23x,a9,1p,1000(1x,2e11.3))') '      Sm:', Sm(Lms,:)

          end do ! end of Loop on Lms
      
        end do ! end of Loop on isp_so
        
      end do ! end of Loop on ir
  
 ! Resolution
      i = n_line + 1 
      do ir = nr,2,-1
      
        do Lmss = nlms,1,-1
          if( Lmss > nlm ) then
            Lms = Lmss - nlm
            isp = 2
          else
            Lms = Lmss
            if( Spinorbite ) isp = 1
          endif
            
          L = Lmso(Lms,1,igrph)
          m = Lmso(Lms,2,igrph)
          iso = Lmso(Lms,3,igrph)
          Lm = L**2 + L + 1 + m
 
          i = i - 1
          if( ir == nr ) then
            do Lmssp = 1,nlms
              if( Lmssp > nlm ) then
                Lmsp = Lmssp - nlm
                iso = 2
              else
                Lmsp = Lmssp
                iso = 1
              endif
              Lmp = Lmso(Lmsp,1,igrph)**2 + Lmso(Lmsp,1,igrph) + 1 + Lmso(Lmsp,2,igrph)
              u(ir,Lmp,Lm,isp,iso) = Sm(Lmss,Lmssp)
            end do 
          endif
 
          do j = 1,nlms
            if( j > nlm ) then
              Lmsp = Lms + j - nlm
              ispp = 2
            else
              Lmsp = Lms + j
              if( Spinorbite ) then
                ispp = 1
              else
                ispp = isp
              endif
            endif

            if( Lmsp <= nlm ) then
              jr = ir
            else
              jr = ir + 1
              if( jr > nr ) exit
              Lmsp = Lmsp - nlm
            endif 
            Lp = Lmso(Lmsp,1,igrph)
            mp = Lmso(Lmsp,2,igrph)
            isop = Lmso(Lmsp,3,igrph)
            Lmp = Lp**2 + Lp + 1 + mp
            u(ir,:,Lm,isp,iso) = u(ir,:,Lm,isp,iso) - Mat(i,j) * u(jr,:,Lmp,ispp,isop)  
          end do
          
        end do
      end do
      
    end do Loop_isp_no_so

    deallocate( Line, Mat, Sm )
 
    if( Spinorbite ) exit
    
  end do Loop_igrph ! end of loop on representations

! Filling of the ir = 1 point using the ratio calculted previously

  do isp = 1,nspinp
    isr = min( isp, nspin )
    Lm = 0
    do L = 0,Lmax
      do m = -L,L
        Lm = Lm + 1
        u(1,:,Lm,isp,:) =  Rap(L,isr) * u(2,:,Lm,isp,:) 
      end do
    end do
  end do

  if( Irregular ) call Cal_tauLL(icheck,Lmaxso,r_Bessel,r_Hankel,nlm_E,nr_cluster,nspin,nspino,nspinp,Spinorbite,u)

  if( icheck > 1 ) then
    if( irregular ) then
      call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,3)
    else
      call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,4)
    endif
  endif

  return
  110 format(/'-------- Cluster irregular solutions ',80('-'))
  115 format(/'-------- Cluster regular solutions ',80('-'))
  120 format(i5,2i4,23x,1p,10000(1x,2e11.3))
  130 format(i5,2i4,1p,10000(1x,2e11.3))
end

!***********************************************************************

! Resolution of the radial Schrodinger equation. Include the case of the non spherical potential
! In ur(ir,m,mp,isp,isol) : (mp, isol) define the basis state
!                           (m, isp) are the real moment and spin

subroutine Cal_tauLL(icheck,Lmaxso,r_Bessel,r_Hankel,nlm_E,nr_cluster,nspin,nspino,nspinp,Spinorbite,u)

  use declarations
  implicit none

  integer, parameter:: Length_col = 27
  integer, parameter:: Length_n = 3

  integer:: icheck, iso, isp, isr, iss, L, Lm, Lm1, Lm2, Lmax, Lmaxso, m, nlm_E, nr, nr_cluster, nspin, nspino, nspinp
  
  character(len=3):: Name
  character(len=Length_col), dimension(nlm_E):: Col_name

  complex(kind=db), dimension(nlm_E,nspinp,nlm_E,nspinp):: TauLL
  complex(kind=db), dimension(nr_cluster-1:nr_cluster,0:Lmaxso,nspin):: r_Bessel, r_Hankel
  complex(kind=db), dimension(nr_cluster,nlm_E,nlm_E,nspinp,nspino):: u
  
  Logical:: L_index, Spinorbite
  
  Lmax = Lmaxso
  
  nr = nr_cluster - 1
  
  do Lm1 = 1,nlm_E
    do iso = 1,nspino
      do isp = 1,nspinp
        isr = min( isp,nspin)
        if( Spinorbite ) then
          iss = iso
        else
          iss = isp
        endif
        Lm2 = 0
        do L = 0,Lmax
          do m = -L,L
            Lm2 = Lm2 + 1
            TauLL(Lm1,iss,Lm2,isp) = u(nr,Lm1,Lm2,isp,iso) / r_Hankel(nr,L,isp)
            if( Lm1 == Lm2 ) TauLL(Lm1,iss,Lm2,isp) = TauLL(Lm1,iss,Lm2,isp) - r_Bessel(nr,L,isp) / r_Hankel(nr,L,isp)
            u(nr_cluster,Lm1,Lm2,isp,iso) = r_Hankel(nr_cluster,L,isp) * TauLL(Lm1,iss,Lm2,isp) 
            if( Lm1 == Lm2 ) u(nr_cluster,Lm1,Lm2,isp,iso) = u(nr_cluster,Lm1,Lm2,isp,iso) + r_Bessel(nr_cluster,L,isp)
          end do
        end do  
      end do  
    end do
  end do  

  if( icheck < 2 ) return
  
  Name = 'Tau'
  L_index = .true.
  
  call Create_col_name(Col_name,Length_col,Length_n,Lmax,L_index,Name,nlm_E,1)
  
  do iso = 1,nspinp
    do isp = 1,nspinp
      if( .not. Spinorbite .and. isp /= iso ) cycle
      if( Spinorbite ) then
        write(3,120) iso, isp
      elseif( nspinp > 1 ) then
        write(3,130) isp
      endif
      if( Spinorbite ) then
        iss = iso
      elseif( nspinp > 1 ) then
        iss = isp
      endif
      write(3,'(/7x,1000a27)') Col_name(:)
      Lm = 0
      do L = 0,Lmax
        do m = -L,L
          Lm = Lm + 1
          write(3,150) L, m, TauLL(Lm,iss,:,isp) 
        end do
      end do
    end do
  end do

  return
  120 format(/' isp_out =',i2,', isp =',i2)
  130 format(/' isp =',i2)
  150 format(2i3,1p,1000(2(1x,2e13.5)))
end

!***********************************************************************

! Resolution of the radial Schrodinger equation. Include the case of the non spherical potential
! In ur(ir,m,mp,isp,isol) : (mp, isol) define the basis state
!                           (m, isp) are the real moment and spin

subroutine Sch_radial_reg_cluster(Ecinetic,E_comp,Eimag,f2,g0,gm,gp,gso,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
         nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
         Relativiste,Rmtg_abs,Spinorbite,u,V,V_cp,V_hubb,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, im, ip, ir, isol, isp, isq, isr, L, L_hubbard, L2, Lmax, Lmaxso, m, m_hubb, mp, ms, &
    n, nlm_cp_E, nlm_r_E, nlm_E, np, nr_cluster, nsol, nspin, nspino, nspinp, numat_abs

  logical:: E_comp, Hubb_a, Hubb_d, Hubb_m, Relativiste, Spinorbite, Ylm_comp

  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nr_cluster,nspin,nlm_cp_E,nlm_cp_E):: V_cp
  complex(kind=db), dimension(nr_cluster,nlm_E,nlm_E,nspinp,nspino):: u

  real(kind=db):: br, ci, Eimag, fac, p, Rmtg_abs, td
  real(kind=db), dimension(nspin):: cr, Ecinetic
  real(kind=db), dimension(nr_cluster):: f2, r_cluster
  real(kind=db), dimension(nr_cluster,nspin):: g0, gm, gp
  real(kind=db), dimension(nr_cluster,nspino):: gso
  real(kind=db), dimension(nr_cluster,nspin,nlm_r_E,nlm_r_E):: V

  Lmax = Lmaxso
  
  u(:,:,:,:,:) = ( 0._db, 0._db )

  n = 0
  do L = 0,Lmax

    L2 = L * ( L + 1 )

    ci = - Eimag / ( 4 * L + 6 )

    br = - numat_abs / ( L + 1._db )
    cr(:) = - ( 2 * numat_abs * br + Ecinetic(:) ) / ( 4 * L + 6 )

! Je fais L'hypothese que le terme de Hubbard ne croise pas les solutions 1 et 2 en cas de spinorbite.

! Development at origin
    do m = -L,L
      n = n + 1
      
      do isol = 1,nspinp ! spin at origin
        do isp = 1,nspinp ! real spin
          if( .not. Spinorbite .and. isp /= isol ) cycle
          isr = min( isp, nspin )

          if( Spinorbite ) then
            if( (m == L .and. isp == 1) .or. (m == -L .and. isp == 2 ) ) then
              nsol = 1
            else
              nsol = 2
            endif
            if( nsol == 1 .and. isp /= isol ) cycle

            if( Relativiste ) then
              if( nsol == 1 .or. isol == 1 ) then
                p = sqrt( ( L + 1 )**2 - ( alfa_sf * numat_abs )**2 )
              else
                p = sqrt( L**2 - ( alfa_sf * numat_abs )**2 )
              endif
            else
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
            endif
            
          else
            if( Relativiste ) then
              p = sqrt( L**2 + L + 1 - ( alfa_sf * numat_abs )**2 )
            else
              p = L + 1._db
            endif
          endif

 ! Becareful, it is from the point of view of the m of the orbital
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
          if( .not. Spinorbite .or. isp == isol ) then
            np = n
          elseif( isp < isol ) then
            np = n + 1
          else
            np = n - 1
          endif
          u(1:2,n,np,isp,isol) = fac * ( 1 + br * r_cluster(1:2) + cr(isr) * r_cluster(1:2)**2 ) * r_cluster(1:2)**p &
                               + img * fac * ci * r_cluster(1:2)**(p+2)
        end do
      end do
    end do
  end do

  do ir = 2,nr_cluster-1
    im = ir - 1
    ip = ir + 1

    do isp = 1,nspinp
      isr = min( isp, nspin )
      isq = 3 - isp
        
      do L = 0,Lmax

        L2 = L * ( L + 1 )

        if( Hubb_a .and. L == L_hubbard( numat_abs ) .and. r_cluster(ir) < Rmtg_abs )  then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif

        do m = -L,L

          n = L2 + 1 + m

          td = g0(ir,isr) + L2 * f2(ir)

          if( Spinorbite ) then
            if( isp == 1 ) then
              ms = m
              mp = m + 1   ! m of the other spin
            else
              ms = - m
              mp = m - 1
            endif
            fac = sqrt( ( L - ms ) * ( L + ms + 1._db ) )
            td = td + ms * gso(ir,isp)
          else
            mp = L + 1 ! in order the "if" below be false
          endif
          np = L2 + 1 + mp

! The second index "m" is the non zero component at origin
          u(ip,:,n,isp,:) = td * u(ir,:,n,isp,:) + gm(ir,isr) * u(im,:,n,isp,:)
          if( E_comp ) u(ip,:,n,isp,:) = u(ip,:,n,isp,:) - img * Eimag * u(ir,:,n,isp,:)

          if( abs(mp) <= L ) u(ip,:,n,isp,:) = u(ip,:,n,isp,:) + fac * gso(ir,isp) * u(ir,:,np,isq,:)

          if( Hubb_m ) then
            do mp = -L,L
              if( Hubb_d .and. m /= mp ) cycle
              np = L2 + 1 + mp
              u(ip,:,n,isp,:) = u(ip,:,n,isp,:) + V_hubb(m,mp,isp,isp) * u(ir,:,np,isp,:)
            end do
          endif

          if( Ylm_comp ) then
            do np = 1,nlm_cp_E
              if( n == np ) cycle
              u(ip,:,n,isp,:) = u(ip,:,n,isp,:) + V_cp(ir,isr,n,np) * u(ir,:,np,isp,:)
            end do
          else
            do np = 1,nlm_r_E
              if( n == np ) cycle
              u(ip,:,n,isp,:) = u(ip,:,n,isp,:) + V(ir,isr,n,np) * u(ir,:,np,isp,:)
            end do
          endif

          u(ip,:,n,isp,:) = - u(ip,:,n,isp,:) * gp(ir,isr)

        end do ! Loop on m
      end do ! Loop on L

    end do ! Loop on isp

  end do ! Loop on ir

  if( icheck > 1 ) call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,1)

  call Renormal_cluster(icheck,konde,Lmax,nlm_E,nr_cluster,nspin,nspino,nspinp,r_cluster,Spinorbite,u)

  if( icheck > 1 ) call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,2)

  return
  110 format(/'-------- Cluster regular solutions ',80('-'))
end

!***********************************************************************

! Normalisation of radial functions using continuity at muffin-tin radius with solutions in vaccum

subroutine Renormal_cluster(icheck,konde,Lmax,nlm,nr_cluster,nspin,nspino,nspinp,r_cluster,Spinorbite,u)

  use declarations
  implicit none
  
  integer, parameter:: Length_col = 27
  integer, parameter:: Length_n = 3

  integer icheck, ir, iso, isp, isr, iss, L, Lm, Lm1, Lm2, Lmax, m, nlm, nr_cluster, nspin, nspino, nspinp

  character(len=Length_n):: Name
  character(len=Length_col), dimension(:), allocatable:: Col_name
  
  complex(kind=db):: a_cp, b_cp, c_cp
  complex(kind=db), dimension(0:Lmax):: d_bess, bess, d_hank, hank
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(nlm,nspinp):: Ampl, Wronsks, u1, u2, Wronske, Tau
  complex(kind=db), dimension(2,0:Lmax,nspin):: bs, hs
  complex(kind=db), dimension(nr_cluster,nlm,nlm,nspinp,nspino):: u

  complex(kind=db), dimension(:,:,:,:), allocatable:: TauLL, uu

  logical:: L_index, Spinorbite 

  real(kind=db):: d0, d01, d02, d1, d2, d12, Radius, Wronskout
  real(kind=db), dimension(nr_cluster):: r_cluster

  allocate( uu(nr_cluster-2:nr_cluster,nlm,nspinp,nspino) )
  
  Radius = r_cluster(nr_cluster-1)
   
  d12 = Radius - r_cluster(nr_cluster-2)
  d01 = r_cluster(nr_cluster) - Radius
  d02 = r_cluster(nr_cluster) - r_cluster(nr_cluster-2)
  d2 = 1 / ( d12 * d02 )
  d1 = 1 / ( d12 * d01 )
  d0 = 1 / ( d01 * d02 )

  Wronsks(:,:) = (0._db, 0._db)
  u1(:,:) = (0._db, 0._db)
  u2(:,:) = (0._db, 0._db)

! Hankel function and its derivative (versus r_cluster, not z) at Rmtg
  do isp = 1,nspin
    call Cal_Hankel(d_hank,hank,konde(isp),Lmax,Radius,.false.)
    hs(1,0:Lmax,isp) = - img * hank(0:Lmax)
    hs(2,0:Lmax,isp) = - img * d_hank(0:Lmax)

    call Cal_Bessel(d_bess,bess,konde(isp),Lmax,Radius)
    bs(1,0:Lmax,isp) = bess(0:Lmax)
    bs(2,0:Lmax,isp) = d_bess(0:Lmax)
  end do

  do ir = nr_cluster-2,nr_cluster
    do Lm = 1,nlm
      uu(ir,Lm,:,:) = u(ir,Lm,Lm,:,:) / r_cluster(ir)
    end do  
  end do

  do isp = 1,nspinp
    iso = min(nspino,isp)

    do Lm = 1,nlm
      a_cp = uu(nr_cluster-2,Lm,isp,iso) * d2 - uu(nr_cluster-1,Lm,isp,iso) * d1 + uu(nr_cluster,Lm,isp,iso) * d0
      b_cp = ( uu(nr_cluster-1,Lm,isp,iso) - uu(nr_cluster-2,Lm,isp,iso)) / d12 &
           - a_cp * ( r_cluster(nr_cluster-1) + r_cluster(nr_cluster-2) )
      c_cp = uu(nr_cluster,Lm,isp,iso) - a_cp * r_cluster(nr_cluster)**2 - b_cp * r_cluster(nr_cluster)

      u1(Lm,isp) = a_cp * Radius**2 + b_cp * Radius + c_cp
      u2(Lm,isp) = 2 * a_cp * Radius + b_cp
    end do
  end do

! The following formula for Wronskout comes from normalisation by sqrt(k/pi) of Bessel and Hankel functions
  Wronskout = 1 / ( pi * Radius**2 )

  Lm = 0
  do L = 0,Lmax
    do m = -L,L
      Lm = Lm + 1
      do isp = 1,nspinp
        isr = min( isp, nspin )
        Wronsks(Lm,isp) = u1(Lm,isp) * hs(2,L,isr) - u2(Lm,isp) * hs(1,L,isr)
        Wronske(Lm,isp) = u1(Lm,isp) * bs(2,L,isr) - u2(Lm,isp) * bs(1,L,isr)
      end do
    end do
  end do

! Normalisation by amplitude, and not amplitude / tau  
  Ampl(:,:) = Wronskout / Wronsks(:,:)
!  Ampl(:,:) = Wronskout / Wronske(:,:)
  
  Tau(:,:) = - Wronske(:,:) / Wronsks(:,:) 

  do Lm1 = 1,nlm  ! outer sphere state
    do Lm2 = 1,nlm
      if( nspino == 2 ) then
        do iso = 1,nspino
          u(:,Lm1,Lm2,:,iso) = u(:,Lm1,Lm2,:,iso) * Ampl(Lm1,iso)  
        end do  
      else
        do isp = 1,nspinp
          u(:,Lm1,Lm2,isp,1) = u(:,Lm1,Lm2,isp,1) * Ampl(Lm1,isp)  
        end do
      endif
    end do
  end do 

  if( icheck > 1 ) then
    write(3,100) konde(:)
    write(3,110)
    Lm = 0
    do L = 0,Lmax
      do m = -L,L
        Lm = Lm + 1
        write(3,'(2i3,1p,12(1x,2e13.5))') L, m, ( Ampl(Lm,isp), Tau(Lm,isp), Wronske(Lm,isp), Wronsks(Lm,isp), &
                                                                            u1(Lm,isp), u2(Lm,isp), isp = 1,nspinp ) 
      end do
    end do
    
    allocate( TauLL(nlm,nspinp,nlm,nspinp) )
    do Lm1 = 1,nlm
      do iso = 1,nspino
        do isp = 1,nspinp
          isr = min( isp,nspin)
          if( Spinorbite ) then
            iss = iso
          else
            iss = isp
          endif
          Lm2 = 0
          do L = 0,Lmax
            do m = -L,L
              Lm2 = Lm2 + 1
              if( Lm1 == Lm2 .and. iss == isp ) then
                TauLL(Lm1,iss,Lm2,isp) = Tau(Lm1,iss)
              else
                TauLL(Lm1,iss,Lm2,isp) = Ampl(Lm1,iss) * u(nr_cluster-1,Lm1,Lm2,isp,iso) / hs(1,L,isr)
              endif 
            end do
          end do
        end do
      end do
    end do

    allocate( Col_name(nlm) )
    Name = 'Tau'
    L_index = .true.
    
    call Create_col_name(Col_name,Length_col,Length_n,Lmax,L_index,Name,nlm,1)
    
    do iso = 1,nspinp
      do isp = 1,nspinp
        if( .not. Spinorbite .and. isp /= iso ) cycle
        if( Spinorbite ) then
          write(3,120) iso, isp
        elseif( nspinp > 1 ) then
          write(3,130) isp
        endif
        if( Spinorbite ) then
          iss = iso
        elseif( nspinp > 1 ) then
          iss = isp
        endif
        write(3,'(/7x,1000a27)') Col_name(:)
        Lm = 0
        do L = 0,Lmax
          do m = -L,L
            Lm = Lm + 1
            write(3,150) L, m, TauLL(Lm,iss,:,isp) 
          end do
        end do
      end do
    end do
      
    deallocate( Col_name, TauLL )
  endif
  
  deallocate( uu )
  
  return
  100 format(/' konde =',2e13.5)
  110 format(/' Amplitudes for normalisation',/'  L  m',13x,'Ampl',23x,'Tau',22x,'Wronske',20x,'Wronsks',23x,'u1',25x,'u2')
  120 format(/' isp_out =',i2,', isp =',i2)
  130 format(/' isp =',i2)
  150 format(2i3,1p,1000(2(1x,2e13.5)))
end

!****************************************************************************************************************************

subroutine Create_col_name(Col_name,Length_col,Length_n,Lmax,L_index,Name,nlm,nspin)

  use declarations
  implicit none
  
  integer:: i, isp, Length, Length_col, Length_n, L, Lm, Lmax, m, nlm, nspin
  
  character(len=Length_n):: Name, Name_b 
  character(len=Length_col):: Word 
  character(len=Length_col), dimension(nspin*nlm):: Col_name 

  Logical:: L_index
  
  i = 0
  
  Name_b = adjustl(Name)

  do isp = 1,nspin
    do Lm = 1,nlm
      do L = 0,Lmax
        if( .not. L_index .and. L > 0 ) exit 
          do m = -L, L
        
          i = i + 1
          
          Word = ' '
          Length = Len_trim(Name_b)
          Length = min( Length, Length_col-2 )
          Word(1:Length) = Name(1:Length)
          
          Length = len_trim(Word) + 1
          if( Length < Length_col - 2 ) Word(Length:Length) = '('
          
          if( L_index ) then
            call ad_number(L,Word,Length_col)
            Length = len_trim(Word) + 1
            if( Length < Length_col-1) Word(Length:Length) = ','
            call ad_number(m,Word,Length_col)
          else
            call ad_number(Lm,Word,Length_col)
          endif
          Length = len_trim(Word) + 1
          
          if( nspin == 2 ) then
            if( Length+2 < Length_col ) then 
              if( isp == 1 ) then
                Word(Length:Length+2) = ',u)'
              else
                Word(Length:Length+2) = ',d)'
              endif
            elseif( Length+1 < Length_col ) then
              if( isp == 1 ) then
                Word(Length:Length+1) = ',u'
              else
                Word(Length:Length+1) = ',d'
              endif
            elseif( Length < Length_col ) then
              if( isp == 1 ) then
                Word(Length:Length+1) = 'u'
              else
                Word(Length:Length+1) = 'd'
              endif
            endif
          elseif( Length < Length_col ) then
            Word(Length:Length) = ')'
          endif
          call center_word( Word, Length_col ) 
          Col_name(i) = Word
        end do
      end do
      if( L_index ) exit
    end do
  end do
    
  return
end
  
!***********************************************************************

! Resolution of the radial Schrodinger equation for the irregular solution by inward integration
! Use Numerov's method
! At R, continuity is with - i * Hankel
! when Classic_irreg = .false. irreg solution is multiplied by Tau_ato.

! Copy of Sch_radial: gm et gp are reversed in the call

! nr = nr_cluster

subroutine Sch_radial_irreg_Numerov(Eimag,Enervide,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
       nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
       Relativiste,Rmtg_abs,Spinorbite,us,V,V_cp,V_hubb,Ylm_comp)

  use declarations
  implicit none

  integer:: i0, icheck, im, ip, ir, isol, isp, isq, isr, L, l_hubbard, L2, Lm, Lmax, Lmaxso, Lmp, m, m_hubb, mp, ms, &
    n, nlm_cp_E, nlm_r_E, nlm_E, np, nr, nr_cluster, nspin, nspino, nspinp, numat_abs, Space

  logical:: Hubb_a, Hubb_d, Hubb_m, Relativiste, Spinorbite, Stop_job, Ylm_comp

  complex(kind=db):: fnormc, z
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(0:Lmaxso):: bess, neum
  complex(kind=db), dimension(nr_cluster-1:nr_cluster,0:Lmaxso,nspin):: Hankel
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nr_cluster,nspin,nlm_cp_E,nlm_cp_E):: V_cp
  complex(kind=db), dimension(nspino*nlm_E,nspino*nlm_E):: Mat_0, Mat_m, Mat_p
  complex(kind=db), dimension(nr_cluster,nspino):: gso
  complex(kind=db), dimension(nr_cluster,nlm_E,nlm_E,nspinp,nspino):: us

  real(kind=db):: a2s4, Eimag, Enervide, fac, h2, hm, hp, Rmtg_abs
  real(kind=db), dimension(nr_cluster):: f2, r_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_r_E,nlm_r_E):: V

  us(:,:,:,:,:) = ( 0._db, 0._db )

  call Coef_rad_Numerov(Eimag,Enervide,f2,gso,nlm_cp_E,nlm_r_E,nr_cluster,nspin,nspino,numat_abs,r_cluster,V,V_cp, &
                        Ylm_comp)

  nr = nr_cluster
  Lmax = Lmaxso
  
  do ir = nr_cluster-1,nr_cluster
    do isp = 1,nspin
      fnormc = sqrt( konde(isp) / pi )
      z = konde(isp) * r_cluster(ir)
      call cbessneu(fnormc,z,Lmax,Lmax,bess,neum)
      Hankel(ir,0:Lmax,isp) = bess(0:Lmax) + img * neum(0:Lmax) 
    end do
  end do

  n = 0
  do L = 0,Lmax

    do m = -L,L
      n = n + 1

      do isol = 1,nspinp  ! Outer sphere spin 
        isp = isol        ! isp is inside cluster spin
        isr = min( isp, nspin )

        do ir = nr_cluster-1,nr_cluster
          us(ir,n,n,isp,isol) = - img * r_cluster(ir) * Hankel(ir,L,isr)
        end do
      end do

    end do
  end do

  hm = r_cluster(nr_cluster) - r_cluster(nr_cluster-1)
 
  do ir = nr_cluster-1,2,-1
    hp = r_cluster(ir) - r_cluster(ir-1)
    if( Space == - 1 ) then
      Space = 2
      h2 = hp**2 / 12
      im = ir + 1
      i0 = ir
      ip = ir - 1
    elseif( hp > hm + eps10 ) then
      Space = 1
      h2 = hp**2 / 12
      im = ir + 2
      i0 = ir
      ip = ir - 1
    elseif(  hp < hm - eps10 ) then
      Space = - 1
      h2 = hm**2 / 12
      im = ir + 1
      i0 = ir
      ip = ir - 2
    else
      Space = 0
      h2 = hp**2 / 12
      im = ir + 1
      i0 = ir
      ip = ir - 1
    endif

    hm = hp ! for the next iteration on ir

    do isp = 1,nspinp   ! true spin inside cluster, isol is the "spin attack"
      isr = min( isp, nspin )
      isq = 3 - isp

      do L = 0,Lmax

        L2 = L * ( L + 1 )

        if( Hubb_a .and. L == L_hubbard(numat_abs) .and. r_cluster(ir) < Rmtg_abs ) then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif

        do m = -L,L

          Lm = L2 + 1 + m
          
          if( Spinorbite .and. isp == 2 ) then
            n = Lm + nlm_E
          else
            n = Lm
          endif
          
          Mat_m(n,n) = 1 + h2 * ( Enervide + img * Eimag - V(im,isr,Lm,Lm) - L2 * f2(im) )
          Mat_0(n,n) = 1 - 5 * h2 * ( Enervide + img * Eimag - V(i0,isr,Lm,Lm) - L2 * f2(i0) )
          Mat_p(n,n) = 1 + h2 * ( Enervide + img * Eimag - V(ip,isr,Lm,Lm) - L2 * f2(ip) )
          
          do Lmp = 1,nlm_E
            if( Lm == Lmp ) cycle
            if( Spinorbite .and. isp == 2 ) then
              np = Lmp + nlm_E
            else
              np = Lmp
            endif
            if( Ylm_comp ) then
              Mat_m(n,np) = - h2 * V_cp(im,isr,Lm,Lmp)
              Mat_0(n,np) = 5 * h2 * V_cp(i0,isr,Lm,Lmp)
              Mat_p(n,np) = - h2 * V_cp(ip,isr,Lm,Lmp)
            else
              Mat_m(n,np) = - h2 * V(im,isr,Lm,Lmp)
              Mat_0(n,np) = 5 * h2 * V(i0,isr,Lm,Lmp)
              Mat_p(n,np) = - h2 * V(ip,isr,Lm,Lmp)
            endif
          end do
          
          if( Hubb_m ) then
            do mp = -L,L
              if( Hubb_d .and. m /= mp ) cycle
              Lmp = L2 + 1 + mp
              if( Spinorbite .and. isp == 2 ) then
                np = Lmp + nlm_E
              else
                np = Lmp
              endif
              Mat_m(n,np) = Mat_m(n,np) - h2 * V_hubb(m,mp,isp,isp)
              Mat_0(n,np) = Mat_0(n,np) + 5 * h2 * V_hubb(m,mp,isp,isp)
              Mat_p(n,np) = Mat_p(n,np) - h2 * V_hubb(m,mp,isp,isp)
            end do
          endif
          
          if( Relativiste ) then
            Mat_m(n,n) = Mat_m(n,n) - h2 * a2s4 * ( Enervide + img * Eimag - V(im,isr,Lm,Lm) )**2
            Mat_0(n,n) = Mat_0(n,n) + 5 * h2 * a2s4 * ( Enervide + img * Eimag - V(i0,isr,Lm,Lm) )**2
            Mat_p(n,n) = Mat_p(n,n) - h2 * a2s4 * ( Enervide + img * Eimag - V(ip,isr,Lm,Lm) )**2

            Mat_m(n,n) = Mat_m(n,n) - h2 * gso(im,isr)
            Mat_0(n,n) = Mat_0(n,n) + 5 * h2 * gso(i0,isr)
            Mat_p(n,n) = Mat_p(n,n) - h2 * gso(ip,isr)
          endif
          
          if( Spinorbite ) then
            if( isp == 1 ) then
              ms = m
              mp = m + 1   ! m of the other spin
            else
              ms = - m
              mp = m - 1
            endif
 
            Mat_m(n,n) = Mat_m(n,n) - h2 * ms * gso(im,isp)
            Mat_0(n,n) = Mat_0(n,n) + 5 * h2 * ms * gso(i0,isp)
            Mat_p(n,n) = Mat_p(n,n) - h2 * ms * gso(ip,isp)

            if( abs(mp) <= L ) then
              Lmp = L2 + 1 + mp
              if( Spinorbite .and. isq == 2 ) then
                np = Lmp + nlm_E
              else
                np = Lmp
              endif
              fac = sqrt( ( L - ms ) * ( L + ms + 1._db ) )
              Mat_m(n,np) = Mat_m(n,np) - h2 * fac * gso(im,isp)
              Mat_0(n,np) = Mat_0(n,np) + 5 * h2 * fac * gso(i0,isp)
              Mat_p(n,np) = Mat_p(n,np) - h2 * fac * gso(ip,isp)
            endif

          endif

        end do ! Loop on m
      end do ! Loop on L

      if( Spinorbite ) cycle

      if( icheck > 3 .and. ( ir < 4 .or. ir > nr_cluster - 3 ) ) then
        if( Space == 2 ) then
          write(3,'(/a30,i6)') '  Mat_0 before inversion, ir =', ir
        else
          write(3,'(/a13,i6)') '  Mat_0, ir =', ir
        endif   
        do n = 1,nlm_E
          write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_0(n,:)   
        end do
        write(3,'(/a13,i6)') '  Mat_m, ir =', ir   
        do n = 1,nlm_E
          write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_m(n,:)   
        end do
        if( Space == 2 ) then
          write(3,'(/a13,i6)') '  Mat_p, ir =', ir
        else
          write(3,'(/a30,i6)') '  Mat_p before inversion, ir =', ir
        endif   
        do n = 1,nlm_E
          write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_p(n,:)   
        end do
      endif

      if( Space == 2 ) then
        call invcomp(nlm_E,Mat_0,nlm_E,nlm_E,1,Stop_job)
      else
        call invcomp(nlm_E,Mat_p,nlm_E,nlm_E,1,Stop_job)
      endif

      if( icheck > 3 .and. ( ir < 4 .or. ir > nr_cluster - 3 ) ) then
        if( Space == 2 ) then
          write(3,'(/a29,i6)') '  Mat_0 after inversion, ir =', ir
          do n = 1,nlm_E
            write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_0(n,:)   
          end do
        else
          write(3,'(/a29,i6)') '  Mat_p after inversion, ir =', ir
          do n = 1,nlm_E
            write(3,'(i3,1p,1000(1x,2e11.3))') n, Mat_p(n,:)   
          end do
        endif   
      endif
      
      if( Space == 2 ) then
      
        Mat_p = Matmul( Mat_0, Mat_p ) 
        Mat_m = Matmul( Mat_0, Mat_m ) 

        do Lm = 1,nlm_E
          do Lmp = 1,nlm_E
            us(i0,:,Lm,isp,:) = us(i0,:,Lm,isp,:) + 0.5_db * ( Mat_p(Lm,Lmp) * us(ip,:,Lmp,isp,:) &
                                                  + Mat_m(Lm,Lmp) * us(im,:,Lmp,isp,:) ) 
          end do
        end do
        
      else
      
        Mat_0 = Matmul( Mat_p, Mat_0 ) 
        Mat_m = Matmul( Mat_p, Mat_m ) 

        do Lm = 1,nlm_E
          do Lmp = 1,nlm_E
            us(ip,:,Lm,isp,:) = us(ip,:,Lm,isp,:) + 2 * Mat_0(Lm,Lmp) * us(i0,:,Lmp,isp,:) - Mat_m(Lm,Lmp) * us(im,:,Lmp,isp,:) 
          end do
        end do
        
      endif

    end do ! Loop on isp

    if( .not. Spinorbite ) cycle
    
    call invcomp(2*nlm_E,Mat_p,2*nlm_E,2*nlm_E,0,Stop_job)
    
    Mat_0 = 2 * Matmul( Mat_p, Mat_0 ) 
    Mat_m = - Matmul( Mat_p, Mat_m ) 

    do isp = 1,nspino 
      do Lm = 1,nlm_E
        n = Lm + ( isp - 1 ) * nlm_E
        do isq = 1,nspino
          do Lmp = 1,nlm_E
            np = Lmp + ( isq - 1 ) * nlm_E
            us(ip,:,Lm,isp,:) = us(ip,:,Lm,isp,:) + Mat_0(n,np) * us(i0,:,lmp,isq,:) + Mat_m(n,np) * us(im,:,lmp,isq,:) 
          end do
        end do
      end do
    end do
    
  end do

  if( icheck > 1 ) call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,us,3)

  return
  110 format(/'-------- Cluster irregular solutions ',80('-'))
end

!***********************************************************************

! Resolution of the radial Schrodinger equation to get the irregular solutions
! Include the case of the non spherical potential
! Uses the matrix way with Gaussian elimination
! In ur(ir,m,mp,isp,isol) : (mp, isol) define the basis state
!                           (m, isp) are the real moment and spin

subroutine Sch_radial_irreg_gaussian(Eimag,Enervide,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
         nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
         Relativiste,Rmtg_abs,Spinorbite,u,V,V_cp,V_hubb,Ylm_comp)

  use declarations
  implicit none

  integer:: i, icheck, ir, isol, isp, isr, j, j2, jr, k, L, L_hubbard, Lm, Lmax, Lmaxso, Lmp, m, m_hubb, &
    n_line, nlm_cp_E, nlm_r_E, nlm_E, nr, nr_cluster, nspin, nspino, nspinp, numat_abs

  logical:: Hubb_a, Hubb_d, Hubb_m, Relativiste, Spinorbite, Ylm_comp

  complex(kind=db):: Diag, fnormc, Line_j2, z
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(0:Lmaxso):: bess, neum
  complex(kind=db), dimension(0:Lmaxso,nspin):: Rap
  complex(kind=db), dimension(nr_cluster-1:nr_cluster,0:Lmaxso,nspin):: r_Hankel
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nr_cluster,nspin,nlm_cp_E,nlm_cp_E):: V_cp
  complex(kind=db), dimension(nr_cluster,nlm_E,nlm_E,nspinp,nspino):: u

  complex(kind=db), dimension(:), allocatable:: Line
  complex(kind=db), dimension(:,:), allocatable:: Mat, Sm

  real(kind=db):: dr, E_centrifuge, r0, rm, rp, Eimag, Enervide, Rmtg_abs

  real(kind=db), dimension(nr_cluster):: clapl0, claplm, claplp, f2, r_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_r_E,nlm_r_E):: V

  Lmax = Lmaxso
  nr = nr_cluster - 2
  
  u(:,:,:,:,:) = ( 0._db, 0._db )

  do ir = 1,nr_cluster - 1
    rp = r_cluster(ir+1)
    if( ir == 1 ) then
      rm = r_cluster(1)**2 / r_cluster(2)
    else
      rm = r_cluster(ir-1)
    endif
    r0 = r_cluster(ir)
    dr = 0.5_db * ( rp - rm )
    claplm(ir) = 1 / ( ( r0 - rm ) * dr )
    claplp(ir) = 1 / ( ( rp - r0 ) * dr )
    clapl0(ir) = - claplm(ir) - claplp(ir)
    f2(ir) = 1 / r_cluster(ir)**2
  end do
  
  do isp = 1,nspin
    fnormc = sqrt( konde(isp) / pi )
    do ir = nr_cluster-1,nr_cluster
      z = konde(isp) * r_cluster(ir)
      call cbessneu(fnormc,z,Lmax,Lmax,bess,neum)
      r_Hankel(ir,0:Lmax,isp) = - img * r_cluster(ir) * ( bess(0:Lmax) + img * neum(0:Lmax) )
    end do 
  end do

! Calculation of the ratio( Rap ) of the irregular solution between index 1 and 2 of ir.
! The guess is that at the origin the regular atomic solution is negligible in front of the irregular one.
  call Isotrope_radial_irreg(clapl0,claplm,claplp,Eimag,Enervide,f2,icheck,Lmaxso,nlm_cp_E,nlm_r_E,nr_cluster, &
                                 nspin,r_cluster,r_Hankel,Rap,V,V_cp,Ylm_comp)

  n_line = ( nr - 1 ) * nlm_E
  allocate( Mat(n_line,nlm_E) )
  allocate( Sm(nlm_E,nlm_E) )
  allocate( Line(n_line) )
  
  Mat(:,:) = ( 0._db, 0._db )
  Sm(:,:) = ( 0._db, 0._db )

  isp = 1
  
! potential taken diagonal
  
  ir = 2
  Lm = 0
  i = 0

  do L = 0,Lmax
  
    E_centrifuge = L * ( L + 1 ) * f2(ir)

    isr = min( isp, nspin )

    do m = -L,L

      Lm = Lm + 1
      i = i + 1
      if( Ylm_comp ) then
        Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V_cp(ir,isp,Lm,Lm) 
      else
        Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V(ir,isp,Lm,Lm) 
      endif
      Diag = Diag - claplm(ir) * Rap(L,isr)
      Mat(i,nlm_E) = - claplp(ir) / Diag   

      if( icheck > 3 ) then
        if( i == 1 )  write(3,'(/A/A)') ' Matrix','   ir  Lm'       
        write(3,120) ir, Lm, Diag, - claplp(ir)
      endif 

    end do
  end do

  do ir = 3,nr

    Lm = 0
    do L = 0,Lmax
      E_centrifuge = L * ( L + 1 ) * f2(ir)

      do m = -L,L
        Lm = Lm + 1
        i = i + 1

! Line filling 
        Line(:) = ( 0._db, 0._db )
        
        if( Ylm_comp ) then
          Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V_cp(ir,isp,Lm,Lm) 
        else
          Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V(ir,isp,Lm,Lm) 
        endif
 
        Line(i-nlm_E) = - claplm(ir)
        
        do Lmp = 1,nlm_E
          if( Lmp == Lm ) cycle
          j = i + Lmp - Lm
          if( Ylm_comp ) then
            Line(j) = V_cp(ir,isp,Lm,Lmp) 
          else
            Line(j) = V(ir,isp,Lm,Lmp) 
          endif
        end do

! Second member filling
        if( ir == nr ) then
          Sm(Lm,Lm) = claplp(ir) * r_Hankel(nr+1,L,isp)
        else
          Line(i+nlm_E) = - claplp(ir)
        endif   

        Line(i) = Diag
                 
        if( icheck > 3 ) then
          if( ir == nr ) then
            write(3,130) ir, Lm, Line(i-nlm_E), Line(i+1-Lm:i+nlm_E-Lm), Sm(Lm,Lm) 
          else
            write(3,130) ir, Lm, Line(i-nlm_E), Line(i+1-Lm:i+nlm_E-Lm), Line(i+nlm_E) 
          endif
        endif

! Triangularisation is immediatly done on the line
        do j2 = i-nlm_E,i-1 
          Line_j2 = Line(j2)  ! not divided by diagonal matrix term, because = 1 (and not not storaged)   
          do k = j2+1,min(j2+nlm_E,n_line)   ! starting at j2+1 and not j2 because Line(j2) is not used anymore (and corresponding matrix term not storaged)
            j = k - j2
            Line(k) = Line(k) - Mat(j2,j) * Line_j2
          end do

          Lmp = j2 - ( n_line - nlm_E )           
          if( Lmp > 0 ) Sm(Lm,:) = Sm(Lm,:) - Sm(Lmp,:) * Line_j2

        end do
        
  ! normalisation by the diagonal matrix term, which is not storaged
        do j = 1,nlm_E
          if( i + j > n_line ) exit 
          Mat(i,j) = Line(i+j) / Line(i)
        end do
        if( ir == nr ) Sm(Lm,:) = Sm(Lm,:) / Line(i)
        
        if( icheck > 3 ) write(3,'(23x,a9,1p,1000(1x,2e11.3))') '  Triang:', Mat(i,:)
        if( ir == nr .and. icheck > 3 ) write(3,'(23x,a9,1p,1000(1x,2e11.3))') '      Sm:', Sm(i,:)
      end do
    end do
    
  end do
  
  ! Resolution
  isp = 1
  isol = 1
  
  i = n_line + 1 
  do ir = nr,2,-1
    do Lm = nlm_E,1,-1
      i = i - 1
      if( ir == nr ) u(ir,:,Lm,isp,isol) = Sm(Lm,:) 
      do j = 1,nlm_E
        Lmp = Lm + j
        if( Lmp <= nlm_E ) then
          jr = ir
        else
          jr = ir + 1
          if( jr > nr ) exit
          Lmp = Lmp - nlm_E
        endif 
        u(ir,:,Lm,isp,isol) = u(ir,:,Lm,isp,isol) - Mat(i,j) * u(jr,:,Lmp,isp,isol)  
      end do
        
    end do
  end do
  
  
  if( icheck > 1 ) call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,3)

  return
  110 format(/'-------- Cluster irregular solutions ',80('-'))
  120 format(i5,i4,23x,1p,10000(1x,2e11.3))
  130 format(i5,i4,1p,10000(1x,2e11.3))
end

!***********************************************************************

! Calculation of the ratio between irregular solution at the ir = 1 and ir = 2 points

! Use the resolution of the radial Schrodinger equation to get the irregular solutions by inward integration in the isotropic case

subroutine Isotrope_radial_irreg(clapl0,claplm,claplp,Eimag,Enervide,f2,icheck,Lmaxso,nlm_cp_E,nlm_r_E,nr_cluster, &
                                 nspin,r_cluster,r_Hankel,Rap,V,V_cp,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, ir, isp, L, Lm, Lmax, Lmaxso, nlm_cp_E, nlm_r_E, nr, nr_cluster, nspin

  logical:: Ylm_comp

  complex(kind=db):: Diag
  complex(kind=db), dimension(nr_cluster-1:nr_cluster,0:Lmaxso,nspin):: r_Hankel
  complex(kind=db), dimension(nr_cluster,nspin,nlm_cp_E,nlm_cp_E):: V_cp
  complex(kind=db), dimension(nr_cluster,0:Lmaxso,nspin):: u

  complex(kind=db), dimension(0:Lmaxso,nspin):: Rap

  real(kind=db):: E_centrifuge, Eimag, Enervide
  real(kind=db), dimension(nr_cluster):: clapl0, claplm, claplp, f2, r_cluster
  real(kind=db), dimension(nr_cluster,nspin,nlm_r_E,nlm_r_E):: V

  Lmax = Lmaxso
  nr = nr_cluster - 1
  
  u(:,:,:) = ( 0._db, 0._db )

  do L = 0,Lmax
    Lm = L**2 + L + 1

    do isp = 1,nspin
    
      do ir = nr,nr_cluster
        u(ir,L,isp) = r_Hankel(ir,L,isp)
      end do

      do ir = nr,2,-1

        E_centrifuge = L * ( L + 1 ) * f2(ir)
        
        if( Ylm_comp ) then
          Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V_cp(ir,isp,Lm,Lm) 
        else
          Diag = - clapl0(ir) - Enervide - img * Eimag + E_centrifuge + V(ir,isp,Lm,Lm) 
        endif

        u(ir-1,L,isp) = ( Diag * u(ir,L,isp) - claplp(ir) * u(ir+1,L,isp) ) / claplm(ir)

      end do
      
      Rap(L,isp) = u(1,L,isp) / u(2,L,isp)
      
    end do
  end do

  if( icheck > 1 ) then
    write(3,110)
    if( nspin == 1 ) then
      write(3,120) ( L, L, L = 0,Lmax )
    else
      write(3,130) ( L, L, L, L, L = 0,Lmax )
    endif
    do ir = 1,nr_cluster
      write(3,140) r_cluster(ir)*bohr, ( u(ir,L,:), L = 0,Lmax )
    end do
  endif

  return
  110 format(/'-------- Isotrope irregular solutions ',80('-'))
  120 format(/'    Radius  ',100(5x,'ur(',i1,')',8x,'ui(',i1,')',4x))
  130 format(/'    Radius  ',100(4x,'ur(',i1,',u)',6x,'ui(',i1,',u)',7x,'ur(',i1,',d)',6x,'ui(',i1,',d)',3x))
  140 format(f11.7,1p,10000(1x,2e13.5))
end

!***********************************************************************

! Resolution of the radial Schrodinger equation for the irregular solution by inward integration
! At R, continuity is with - i * Hankel
! when Classic_irreg = .false. irreg solution is multiplied by Tau_ato.

! Copy of Sch_radial: gm et gp are reversed in the call

! nr = nr_cluster

subroutine Sch_radial_irreg_cluster(E_comp,Eimag,f2,g0,gp,gm,gso,Hubb_a,Hubb_d,icheck,konde,Lmaxso,m_hubb, &
       nlm_cp_E,nlm_r_E,nlm_E,nr_cluster,nspin,nspino,nspinp,numat_abs,r_cluster, &
       Rmtg_abs,Spinorbite,us,V,V_cp,V_hubb,Ylm_comp)

  use declarations
  implicit none

  integer:: icheck, im, ip, ir, isol, isp, isq, isr, L, l_hubbard, L2, Lmax, Lmaxso, m, m_hubb, mp, ms, n, &
    nlm_cp_E, nlm_r_E, nlm_E, np, nr, nr_cluster, nspin, nspino, nspinp, numat_abs

  logical:: E_comp, Hubb_a, Hubb_d, Hubb_m, Spinorbite, Ylm_comp

  complex(kind=db):: fnormc, z
  complex(kind=db), dimension(nspin):: konde
  complex(kind=db), dimension(0:Lmaxso):: bess, neum
  complex(kind=db), dimension(nr_cluster-1:nr_cluster,0:Lmaxso,nspin):: Hankel
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nr_cluster,nspin,nlm_cp_E,nlm_cp_E):: V_cp
  complex(kind=db), dimension(nr_cluster,nlm_E,nlm_E,nspinp,nspino):: us

  real(kind=db):: Eimag, fac, fnorm, konder, Rmtg_abs, td, zr
  real(kind=db), dimension(nr_cluster):: f2, r_cluster
  real(kind=db), dimension(0:Lmaxso):: bessr, neumr
  real(kind=db), dimension(nr_cluster,nspin):: g0, gm, gp
  real(kind=db), dimension(nr_cluster,nspino):: gso
  real(kind=db), dimension(nr_cluster,nspin,nlm_r_E,nlm_r_E):: V
  
  us(:,:,:,:,:) = ( 0._db, 0._db )

  nr = nr_cluster
  Lmax = Lmaxso
  
  do ir = nr_cluster-1,nr_cluster
    do isp = 1,nspin
      if( E_comp ) then
        fnormc = sqrt( konde(isp) / pi )
        z = konde(isp) * r_cluster(ir)
        call cbessneu(fnormc,z,Lmax,Lmax,bess,neum)
        Hankel(ir,0:Lmax,isp) = bess(0:Lmax) + img * neum(0:Lmax) 
      else
        konder = Real( konde(isp), db )
        fnorm = sqrt( konder / pi )
        zr = konder * r_cluster(ir)
        call cbessneur(fnorm,zr,Lmax,Lmax,bessr,neumr)
        Hankel(ir,0:Lmax,isp) =  bessr(0:Lmax) + img * neumr(0:Lmax)
      endif
    end do
  end do

  n = 0
  do L = 0,Lmax

    do m = -L,L
      n = n + 1

      do isol = 1,nspinp  ! Outer sphere spin 
        isp = isol        ! isp is inside cluster spin
        isr = min( isp, nspin )

        do ir = nr_cluster-1,nr_cluster
          us(ir,n,n,isp,isol) = - img * r_cluster(ir) * Hankel(ir,L,isr)
        end do
      end do

    end do
  end do

  do ir = nr_cluster-1,2,-1
    ip = ir - 1      ! inverse par rapport a Sch_radial
    im = ir + 1

    do isp = 1,nspinp   ! true spin inside cluster, isol is the "spin attack"
      isr = min( isp, nspin )
      isq = 3 - isp

      do L = 0,Lmax

        L2 = L * ( L + 1 )

        if( Hubb_a .and. L == L_hubbard(numat_abs) .and. r_cluster(ir) < Rmtg_abs ) then
          Hubb_m = .true.
        else
          Hubb_m = .false.
        endif

        do m = -L,L

          n = L2 + 1 + m

          td = g0(ir,isr) + L2 * f2(ir)

          if( Spinorbite ) then
            if( isp == 1 ) then
              ms = m
              mp = m + 1   ! m of the other spin
            else
              ms = - m
              mp = m - 1
            endif
            fac = sqrt( ( L - ms ) * ( L + ms + 1._db ) )
            td = td + ms * gso(ir,isp)
          else
            mp = L + 1 ! in order the "if" below be false
          endif
          np = L2 + 1 + mp

! the second index "m" is the non zero component at origin
          us(ip,:,n,isp,:) = td * us(ir,:,n,isp,:) + gm(ir,isr) * us(im,:,n,isp,:)
          if( E_comp ) us(ip,:,n,isp,:) = us(ip,:,n,isp,:) - img * Eimag * us(ir,:,n,isp,:)

          if( abs(mp) <= L ) us(ip,:,n,isp,:) = us(ip,:,n,isp,:) + fac * gso(ir,isr) * us(ir,:,np,isq,:)

          if( Hubb_m ) then
            do mp = -L,L
              if( Hubb_d .and. m /= mp ) cycle
              np = L2 + 1 + mp
              us(ip,:,n,isp,:) = us(ip,:,n,isp,:) + V_hubb(m,mp,isp,isp) * us(ir,:,np,isp,:)
            end do
          endif

          if( Ylm_comp ) then
            do np = 1,nlm_cp_E
              if( n == np ) cycle
              us(ip,:,n,isp,:) = us(ip,:,n,isp,:) + V_cp(ir,isr,n,np) * us(ir,:,np,isp,:)
            end do
          else
            do np = 1,nlm_r_E
              if( n == np ) cycle
              us(ip,:,n,isp,:) = us(ip,:,n,isp,:) + V(ir,isr,n,np) * us(ir,:,np,isp,:)
            end do
          endif

          us(ip,:,n,isp,:) = - us(ip,:,n,isp,:) * gp(ir,isr)

        end do
      end do
    end do
  end do

  if( icheck > 1 ) call Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,us,3)

  return
  110 format(/'-------- Cluster irregular solutions ',80('-'))
end

!***********************************************************************

! Calculation of the integral
! Called by Radial_integral_FDM_comp

Subroutine Integral_FDM_comp_cluster(icheck,ip_max,ip0,nbseuil,nlm_E,nlm_probe,nr_core,NRIXS, &
            nspino,nspinp,Psi_core,r_core,r_or_bess,Rad_ioRIoi,Spinorbite,u_irreg_i,u_irreg_r,u_reg_i,u_reg_r,Vecond)

  use declarations
  implicit none

  integer:: icheck, ip_max, ip0, ip1, ip2, is, iso, isp, isp_out, L, Lm, Lm1, Lm2, Lmax, m, nbseuil, &
    nlm_E ,nlm_probe, nr_core, nspino, nspinp

  complex(kind=db):: integr_sing, Sing
  complex(kind=db), dimension(nr_core):: f_reg, f_irg
  complex(kind=db), dimension(nlm_probe,nlm_probe,nspinp,ip_max-ip0+1,ip_max-ip0+1,nbseuil):: Rad_ioRIoi

  logical:: NRIXS, Reg_irreg, Spinorbite

  real(kind=db):: I_1r, I_1i, I_2r, I_2i
  real(kind=db):: f_integr3, Fac, Fac1, Fac2, Rmtg_abs
  real(kind=db), dimension(nr_core):: fct, r_core, phi1, phi2
  real(kind=db), dimension(nbseuil):: Vecond
  real(kind=db), dimension(nr_core,nbseuil):: Psi_core
  real(kind=db), dimension(nr_core,ip0:ip_max):: r_or_bess
  real(kind=db), dimension(nr_core,nlm_E,nlm_probe,nspinp,nspino):: u_irreg_i, u_irreg_r, u_reg_i,u_reg_r

  Reg_irreg = .true.
  
  if( icheck > 0 ) &
    write(3,'(/A)') ' ----- Integral_FDM_comp_cluster ----------------------------------------------------------------------------'

  Rmtg_abs = r_core(nr_core)

  do is = 1,nbseuil

    do ip1 = ip0,ip_max  ! Loop on dipole, quadrupole, Octupole

      if( NRIXS ) then
        fac1 = 1._db
      else
        select case(ip1)
          case(0)
            fac1 = - 0.5_db * alfa_sf
          case(1)
            fac1 = 1._db
          case(2)
            fac1 = 0.5_db * Vecond(is)
          case(3)
            fac1 = - ( 1._db / 6 ) * Vecond(is)**2
        end select
      endif

      phi1(:) = Psi_core(:,is) * r_or_bess(:,ip1)

      do ip2 = ip0,ip_max  ! loop on multipoles

        if( NRIXS ) then
          fac2 = 1._db
        else
          select case(ip2)
            case(0)
              fac2 = - 0.5_db * alfa_sf
            case(1)
              fac2 = 1._db
            case(2)
              fac2 = 0.5_db * Vecond(is)
            case(3)
              fac2 = - ( 1._db / 6 ) * Vecond(is)**2
          end select
        endif
        
        Fac = Fac1 * Fac2

        phi2(:) = Psi_core(:,is) * r_or_bess(:,ip2)

        do Lm1 = 1,nlm_probe

          do iso = 1,nspino ! solution index of regular solution
            
            do Lm2 = 1,nlm_probe
                      
              do isp = 1,nspinp  ! true spin which is also the index the irregular solution

                do Lm = 1,nlm_E  ! sum on the cluster solutions  
                      
                  do isp_out = 1,nspinp ! Attacking spin 
                      
                    if( .not. Spinorbite .and. isp /= isp_out ) cycle
 
                     if( Reg_irreg ) then
 !Reg and irreg are already multiplied by "r_core"
                      f_reg(:) = cmplx( u_reg_r(:,Lm,Lm1,isp,iso), u_reg_i(:,Lm,Lm1,isp,iso), db )
                      
                      f_irg(:) = cmplx( u_irreg_r(:,Lm,Lm2,isp,isp_out), u_irreg_i(:,Lm,Lm2,isp,isp_out), db )

    ! The minus is because cross section is related to - G(r_core,r_core')
                      Sing = - Fac * integr_sing(nr_core,phi1,phi2,f_reg,f_irg,Rmtg_abs,r_core,icheck)
                    else
                      fct(:) = phi1(:) * u_reg_r(:,Lm,Lm1,isp,iso)
                      I_1r = f_integr3(r_core,fct,1,nr_core,Rmtg_abs) 
                      fct(:) = phi1(:) * u_reg_i(:,Lm,Lm1,isp,iso)
                      I_1i = f_integr3(r_core,fct,1,nr_core,Rmtg_abs) 
                      fct(:) = phi1(:) * u_reg_r(:,Lm,Lm2,isp,iso)
                      I_2r = f_integr3(r_core,fct,1,nr_core,Rmtg_abs) 
                      fct(:) = phi1(:) * u_reg_i(:,Lm,Lm2,isp,iso)
                      I_2i = f_integr3(r_core,fct,1,nr_core,Rmtg_abs) 
                      Sing = - img * Fac * cmplx( I_1r * I_2r + I_1i * I_2i, I_1i * I_2r - I_1r * I_2i )
                   endif

! Sum on the outer sphere solutions
! The 2 solutions of regular solution are also summed up
                    Rad_ioRIoi(Lm1,Lm2,isp,ip1,ip2,is) = Rad_ioRIoi(Lm1,Lm2,isp,ip1,ip2,is) + Sing

                    if( icheck > 1 .and. abs( Sing ) > Fac * eps10 ) &
                        write(3,100) Lm1, Lm2, Lm, isp, isp_out, ip1, ip2, Fac, Sing, Rad_ioRIoi(Lm1,Lm2,isp,ip1,ip2,is) 

                  end do
            
                end do    ! end of loop Lm
              end do   ! end of loop isp 
            end do  ! end of loop Lm2
          end do   
        end do  ! end of loop Lm1

      end do ! end of loop dipole, quadrupole
    end do ! end of loop dipole, quadrupole

  end do ! end of loop seuil
  
  if( icheck > 0 ) then
    Lmax = nint( sqrt( real( nlm_probe ) ) ) - 1
    write(3,'(/A)') ' Radial integral'
    do is = 1,nbseuil
      if( nbseuil > 1 ) write(3,110) is
      do ip1 = ip0,ip_max
        do ip2 = ip0,ip_max
          if( ip1 /= ip2 .and. icheck == 1 ) cycle
          write(3,120) ip1, ip2
          write(3,130) (( L, m, m = -L,L), L = 0,Lmax) 
          do Lm1 = 1,nlm_probe
            do isp = 1,nspinp
              write(3,140) Lm1, isp, Rad_ioRIoi(Lm1,:,isp,ip1,ip2,is) 
            end do
          end do
        end do
      end do
    end do
    
  endif
  
  return
  100 format(/' Lm1 =',i2,', Lm2 =',i2,', Lm =',i3,', isp =',i2,', isp_out =',i2,', ip1 =',i2,', ip2 =',i2, &
              1p,', Fac =',e13.5,', Sing =',2e13.5,', Rad_ioRIoi =',2e13.5)
  110 format(/' Edge',i2)
  120 format(/' Multipole: ip1 =',i2,', ip2 =',i2)
  130 format(' Lm1 isp  (L2,m2) =',1000(1x,'(',i2,',',i2,')',19x))
  140 format(2i4,1x,1p,1000(1x,2e13.5))
end

!****************************************************************************************************************

! Calculation of the cartesian tensors

! < g_1 | o*( l_s m_s, k_s, j_s, l_s, irang ) | f_1 > < f_2 | o*( l_i, m_i, k_i, j_i, l_i, jrang  ) | g_2 >  

! n_Ec = ninit si TDDFT et XANES et Core_resolved
!      = nbseuil si TDDFT et XANES sans core_resolved
!      = 1 en DFT XANES
!      = 2 en Optic
! n_V = ninit si TDDFT et Core_resolved mais pas Optic
!     = nbseuil si TDDFT sans core_recolved mais pas Optic
!     = 1 en DFT ou optic

subroutine tenseur_car_cluster(coef_g,Core_resolved,Ecinetic, &
                Eimag,Energ,Enervide,Eseuil,FDM_comp,Final_optic,Final_tddft,Full_potential,Green,Green_int,Hubb_a,Hubb_d, &
                icheck,ie,ip_max,ip0,is_g,lmax,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi, &
                Multipole, &
                n_Ec,n_oo,n_rel,n_V,nbseuil,ns_dipmag,ndim2,nenerg_tddft,ninit1,ninit,ninitr,ninitv,nlm_pot,nlm_probe, &
                nlm_p_fp,nlmamax,nr,nrm,nspin,nspino,nspinp,numat,psii,r,Rad_ioRIoi,Relativiste,Renorm, &
                Rmtg,Rmtsd,rof0,rot_atom_abs,Rot_int, &
                secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m, &
                secqq,secqq_m,Solsing,Solsing_only,Spinorbite,Taull,Tddft,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp)

  use declarations
  implicit none

  integer:: h_s, hh_s, hh_i, h_i, i, icheck, ie, ief, initr, initlt, ip, ip_max, ip0, ipr, irang, irang1, is, isi, isol, &
    isp, isp1, isp2,isp3, isp4, ispfg, j, j_i, j_s, jh_i, jh_s, jj_i, jj_s, jjh_i, jjh_s, jrang, k, k_i, k_s, kk_i, kk_s, &
    l_i, l_s, lm, lm_i, lm_s, lmax, lmax_pot, lmomax, lomax, lseuil, m_hubb, m_i, m_s, mpinodes, mpirank, mpirank0, &
    n_Ec, n_oo, n_rel, n_V, nbseuil, ndim2, nenerg_tddft, ninit1, ninit, &
    ninitr, ninitv, nh_i, nh_s, nj_i, nj_s, nlm_p_fp, nlm_pot, nlm_probe, nlm1g, &
    nlm2g, nlmamax, nr, nr_zet, nrang, nrm, ns_dipmag, nspin, nspino, nspinp, numat

  parameter( lomax = 3, lmomax = ( lomax + 1 )**2 )

  complex(kind=db), dimension(ninitr):: Te, Te_m, Ten, Ten_m, Tens, Tens_m
  complex(kind=db), dimension(3,3) :: Mat2
  complex(kind=db), dimension(3,3,3) :: Mat3
  complex(kind=db), dimension(3,3,3,3) :: Mat4
  complex(kind=db), dimension(3,3,3,3,3,3) :: Mat6
  complex(kind=db), dimension(3,3,ninitr,0:mpinodes-1):: secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitr,0:mpinodes-1):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitr,0:mpinodes-1):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,3,n_rel,ninitr,0:mpinodes-1):: secdd, secdd_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:mpinodes-1):: secoo, secoo_m
  complex(kind=db), dimension(lmomax,lmomax,n_rel,ninitr):: Tens_lm, Tens_lm_m
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull
  complex(kind=db), dimension(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil):: rof0
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb
  complex(kind=db), dimension(nlm_probe,nlm_probe,nspinp,ip_max-ip0+1,ip_max-ip0+1,nbseuil):: Rad_ioRIoi
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: Singul
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: rof

  integer, dimension(ninit,2):: m_g
  integer, dimension(ninit):: is_g
  integer, dimension(3):: ldip
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3):: msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(3,n_oo,3,n_oo):: msymoo, msymooi

  logical:: Core_resolved, Dip_rel, E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, FDM_comp, Final_optic, &
    Final_tddft, Full_potential, Green, Green_int, Hubb_a, Hubb_d, lmoins1, lplus1, M_depend, M1M1, No_diag, NRIXS, Relativiste, &
    Renorm, Solsing, Solsing_only, Spinorbite, Tddft, Ylm_comp

  logical, dimension(10):: Multipole

  real(kind=db):: c, c0, c1, c12, c120, c3, c5, c8, clm_s, clm_i, Eimag, Energ, Enervide_t, Rmtg, Rmtsd, V_intmax
  real(kind=db), dimension(1):: Energ_t
  real(kind=db), dimension(nspin):: Ecinetic_e, V0bd_e
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitr):: Tensi, Tensr
  real(kind=db), dimension(0:lomax,3,3,3,lmomax):: clm
  real(kind=db), dimension(3,3):: rot_atom_abs, Rot_int, rot_tem
  real(kind=db), dimension(ninit,2):: coef_g
  real(kind=db), dimension(nspin,n_Ec):: Ecinetic
  real(kind=db), dimension(n_Ec):: Eimag_t, Enervide
  real(kind=db), dimension(nspin,n_V):: V0bd
  real(kind=db), dimension(nr):: r
  real(kind=db), dimension(nr,nlm_pot,nspin):: Vrato_e
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(nr,nlm_pot,nspin,n_V):: Vrato
  real(kind=db), dimension(nr,ip0:ip_max):: r_or_bess
  real(kind=db), dimension(:,:,:,:,:,:), allocatable:: zet
  real(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: roff_rr, roff0

  if( icheck > 1 ) write(3,100)

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  NRIXS = .false.
  Dip_rel = n_rel > 1

! Calculation of radial integrals for regular and irregular solutions
  if( Final_optic ) then

    No_diag = Full_potential .or. ( Hubb_a .and. .not. Hubb_d )
    M_depend = Spinorbite .or. Hubb_a .or. No_diag

    if( No_diag ) then
      nlm1g = ( lmax + 1 )**2
      nlm2g = nlm1g
    elseif( M_depend ) then
      nlm1g = ( lmax + 1 )**2
      nlm2g = 1
    else
      nlm1g = lmax + 1
      nlm2g = 1
    endif
    allocate( roff_rr(nlm1g,nlm1g,nlm2g,nlm2g,nspinp**2,nspino**2,ip0:ip_max) )
    nr_zet = 0
    allocate( zet(nr_zet,nlm1g,nlm2g,nspinp,nspino,n_Ec) )
    ief = 2
    allocate( roff0(ief-1,ief:n_Ec,nlm1g,nlm1g,nspinp,nspinp,nspino**2) )
    Eimag_t(:) = Eimag
    Energ_t(1) = Energ

    call radial_optic(Ecinetic,Eimag_t,Energ_t,Enervide,Full_potential,Hubb_a,Hubb_d,icheck,ief,ip_max,ip0, &
         lmax,lmax_pot,m_depend,m_hubb,1,n_Ec,n_V,nlm_pot,nlm1g,nlm2g,No_diag,nr,nr_zet,nspin,nspino,nspinp,numat,r,Relativiste, &
         Renorm,Rmtg,Rmtsd,roff_rr,roff0,Solsing,Spinorbite,V_hubb,V_intmax,V0bd,Vrato,Ylm_comp,zet)
    
    deallocate( roff0, zet )

  else

    do ip = ip0,ip_max
      select case(ip)
        case(0)
          r_or_bess(1:nr,ip) = 1._db
        case(1)
          r_or_bess(1:nr,ip) = r(1:nr)
        case default
          r_or_bess(1:nr,ip) = r(1:nr)**ip
      end select
    end do    
   
    allocate( Singul(nlm_probe,nspinp,nlm_probe,nspinp,ip0:ip_max,ip0:ip_max,ninitv) )
    allocate( rof(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,ninitv) )

    if( .not. FDM_comp ) then

      rof(:,:,:,:,:,:) = (0._db, 0._db)
      Singul(:,:,:,:,:,:,:) = (0._db, 0._db)

      do initlt = 1,n_Ec
        if( Final_tddft ) then
          Ecinetic_e(:) = Ecinetic(:,initlt)
        else
          Ecinetic_e(:) = Ecinetic(:,1)
        endif
        if( Final_tddft ) then
          do isp = 1,nspin
            V0bd_e(isp) = V0bd(isp,initlt)
            Vrato_e(1:nr,:,isp) = Vrato(1:nr,:,isp,initlt)
          end do
        else
          do isp = 1,nspin
            V0bd_e(isp) = V0bd(isp,1)
            Vrato_e(1:nr,:,isp) = Vrato(1:nr,:,isp,1)
          end do
        endif
        Enervide_t = Enervide(initlt)

        call Radial_integral(Ecinetic_e,Eimag,Energ,Enervide_t,Eseuil,Final_tddft,Full_potential,Hubb_a,Hubb_d,icheck, &
           initlt,ip_max,ip0,lmax,lmax_pot,m_hubb,nbseuil,ninit1,ninitv,nlm_pot,nlm_probe,nlm_p_fp,nr,NRIXS,nrm, &
           nspin,nspino,nspinp,numat,psii,r,r_or_bess,Relativiste,Renorm,Rmtg,Rmtsd,rof,Singul,Solsing,Spinorbite,V_hubb,V_intmax, &
           V0bd_e,Vrato_e,Ylm_comp)
      
      end do

    endif

    if( Tddft .and. .not. Final_tddft ) then
      do isol = 1,nspino
        do isp = 1,nspinp
          do lm = 1,nlm_probe
            do i = 1,nbseuil
              if( nlm_p_fp == 1 ) then
                rof0(ie,lm,isp,isol,i) = rof(lm,1,isp,isol,0,i)
              else
! ce qui sert de base en tddft, ce sont les (l_i,m_i,s) et pas les (l,m) vrais.
                rof0(ie,lm,isp,isol,:) = sum( rof(1:nlm_probe,lm,isp,isol,0,i) )
              endif
            end do
          end do
        end do
      end do
    endif

  endif

  if( icheck > 1 ) write(3,110)

  if( E1E3 .or. E3E3 ) then
    nrang = 3
  elseif( E1E2 .or. E2E2 ) then
    nrang = 2
  else
    nrang = 1
  endif
  if( E1M1 .or. M1M1 ) then
    irang1 = 0
  else
    irang1 = 1
  endif

  if( nrang > lomax .and. mpirank0 == 0 ) then
    call write_error
    do ipr = 6,9,3
      write(ipr,'(/A)') ' nrang > lomax in coabs.f !'
    end do
    stop
  endif

! Calcul des composantes de la transf base-cartesienne - base-spherique:

  do irang = irang1,nrang

    clm(irang,:,:,:,:) = 0._db

    select case(irang)

      case(0)
! Magnetic dipole
        c = sqrt( 4 * pi )
! Dans ce cas, correspond a l=0,m=0 mais sert a definir x, y ou z.
        clm(irang,1,:,:,4) = c
        clm(irang,2,:,:,2) = c
        clm(irang,3,:,:,3) = c

      case(1)
! Dipole
        c = sqrt( 4 * pi / 3 )
        clm(irang,1,:,:,4) = c
        clm(irang,2,:,:,2) = c
        clm(irang,3,:,:,3) = c

      case(2)
! Quadrupole
        c0 = sqrt( 4 * pi ) / 3
        c = sqrt( 4 * pi / 15 )
        c3 = c / sqrt( 3._db )

        clm(irang,1,1,:,1) = c0
        clm(irang,1,1,:,7) = - c3;   clm(irang,1,1,:,9) = c
        clm(irang,1,2,:,5) = c
        clm(irang,1,3,:,db) = c
        clm(irang,2,2,:,1) = c0
        clm(irang,2,2,:,7) = - c3;   clm(irang,2,2,:,9) = - c
        clm(irang,2,3,:,6) = c
        clm(irang,3,3,:,1) = c0
        clm(irang,3,3,:,7) = 2 * c3

        do i = 1,3
          do j = i+1,3
            clm(irang,j,i,:,:) = clm(irang,i,j,:,:)
          end do
        end do

      case(3)
! Octupole
        c1 = sqrt( 4 * pi / 75 )

        c = sqrt( 4 * pi / 35 )
        c5 = c / sqrt( 5._db )
        c8 = c / sqrt( 2._db )
        c12 = c / sqrt( 3._db )
        c120 = c / sqrt( 30._db )

        clm(irang,1,1,1,4) = 3 * c1
        clm(irang,1,1,1,14) = - 3 * c120; clm(irang,1,1,1,16) = c8
        clm(irang,2,2,2,2) = 3 * c1
        clm(irang,2,2,2,12) = - 3 * c120; clm(irang,2,2,2,10) = - c8
        clm(irang,3,3,3,3) = 3 * c1
        clm(irang,3,3,3,13) = 2 * c5
        clm(irang,1,1,2,2) = c1
        clm(irang,1,1,2,10) = c8;        clm(irang,1,1,2,12) = - c120
        clm(irang,1,2,2,4) = c1
        clm(irang,1,2,2,16) = - c8;      clm(irang,1,2,2,14) = - c120
        clm(irang,1,1,3,3) = c1
        clm(irang,1,1,3,13) = - c5;      clm(irang,1,1,3,15) = c12
        clm(irang,2,2,3,3) = c1
        clm(irang,2,2,3,13) = - c5;      clm(irang,2,2,3,15) = - c12
        clm(irang,2,3,3,2) = c1
        clm(irang,2,3,3,12) = 4 * c120
        clm(irang,1,3,3,4) = c1
        clm(irang,1,3,3,14) = 4 * c120
        clm(irang,1,2,3,11) = c12

        do i = 1,3
          do j = i,3
            do k = j,3
              if( i == j .and. i == k ) cycle
              clm(irang,i,k,j,:) = clm(irang,i,j,k,:)
              clm(irang,j,i,k,:) = clm(irang,i,j,k,:)
              clm(irang,k,i,j,:) = clm(irang,i,j,k,:)
              clm(irang,j,k,i,:) = clm(irang,i,j,k,:)
              clm(irang,k,j,i,:) = clm(irang,i,j,k,:)
            end do
          end do
        end do

    end select

  end do

  secmm(:,:,:,mpirank) = (0._db,0._db)
  secmd(:,:,:,mpirank) = (0._db,0._db)
  secdq(:,:,:,:,mpirank) = (0._db,0._db)
  secqq(:,:,:,:,:,mpirank) = (0._db,0._db)
  secdo(:,:,:,:,:,mpirank) = (0._db,0._db)
  secdd(:,:,:,:,mpirank) = (0._db,0._db)

  if( E3E3 ) secoo(:,:,:,:,:,mpirank) = (0._db,0._db)

  if( Green_int ) then
    secmm_m(:,:,:,mpirank) = (0._db,0._db)
    secmd_m(:,:,:,mpirank) = (0._db,0._db)
    secdq_m(:,:,:,:,mpirank) = (0._db,0._db)
    secqq_m(:,:,:,:,:,mpirank) = (0._db,0._db)
    secdo_m(:,:,:,:,:,mpirank) = (0._db,0._db)
    secdd_m(:,:,:,:,mpirank) = (0._db,0._db)
    if( E3E3 ) secoo_m(:,:,:,:,:,mpirank) = (0._db,0._db)
  endif

  if( icheck > 1 .and. .not. Final_optic ) write(3,120) lseuil

! Loops over tensor rank
  do irang = irang1,nrang
    do jrang = irang,nrang
      Tens_lm(:,:,:,:) = (0._db,0._db)
      Tens_lm_m(:,:,:,:) = (0._db,0._db)
      if( .not. E1E1 .and. irang == 1 .and. jrang == 1 ) cycle
      if( .not. E1E2 .and. irang == 1 .and. jrang == 2 ) cycle
      if( .not. E2E2 .and. irang == 2 .and. jrang == 2 ) cycle
      if( .not. E1E3 .and. irang == 1 .and. jrang == 3 ) cycle
      if( .not. M1M1 .and. irang == 0 .and. jrang == 0 ) cycle
      if( .not. E1M1 .and. irang == 0 .and. jrang == 1 ) cycle
      if( .not. E3E3 .and. irang == 3 .and. jrang == 3 ) cycle
      if( irang == 0 .and. jrang > 1 ) cycle
      if( jrang == 3 .and. irang == 2 ) cycle

      if( icheck > 1 ) write(3,130) irang, jrang

! Loops over tensor indices
      do k_s = 1,3

        if( irang == 1 .and. ldip(k_s) /= 1 ) cycle

        if( irang > 1 ) then
          nj_s = 3
        else
          nj_s = 1
        endif

        do j_s = 1,nj_s

          if( irang == 2 .and. lqua(k_s,j_s) /= 1 ) cycle

          if( irang == 3 ) then
            nh_s = 3
          else
            nh_s = 1
          endif

          do h_s = 1,nh_s

            if( irang == 3 .and. loct(k_s,j_s,h_s) /= 1 ) cycle
            jh_s = 3 * ( j_s - 1 ) + h_s

            do k_i = 1,3

              if( jrang == 1 .and. ldip(k_i) /= 1 ) cycle

              if( jrang > 1 ) then
                nj_i = 3
              else
                nj_i = 1
              endif

              do j_i = 1,nj_i

                if( jrang == 2 .and. lqua(k_i,j_i) /= 1 ) cycle

                if( jrang == 3 ) then
                  nh_i = 3
                else
                  nh_i = 1
                endif

                do h_i = 1,nh_i

                  if( jrang == 3 .and. loct(k_i,j_i,h_i) /= 1 ) cycle
                  jh_i = 3 * ( j_i - 1 ) + h_i

                  if( irang == 0 .and. jrang == 0 ) then
                    if( k_i < k_s ) cycle
                  elseif( irang == 1 .and. jrang == 1 ) then
                    if( k_i < k_s ) cycle
                    if( msymdd(k_s,k_i) == 0 .and. msymddi(k_s,k_i) == 0 ) cycle
                    if( sum( abs( secdd(k_s,k_i,1,:,mpirank) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 1 .and. jrang == 2 ) then
                    if( msymdq(k_s,k_i,j_i) == 0 .and. msymdqi(k_s,k_i,j_i) == 0 ) cycle
                    if( sum( abs( secdq(k_s,k_i,j_i,:,mpirank) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 2 .and. jrang == 2 ) then
                    if( msymqq(k_s,j_s,k_i,j_i) == 0 .and. msymqqi(k_s,j_s,k_i,j_i) == 0 ) cycle
                    if( sum( abs( secqq(k_s,j_s,k_i,j_i,:,mpirank) ) ) > 1.e-15_db ) cycle
                  elseif( irang == 1 .and. jrang == 3 ) then
                    if( msymdo(k_s,k_i,j_i,h_i) == 0 .and. msymdoi(k_s,k_i,j_i,h_i) == 0 ) cycle
                  elseif( irang == 3 .and. jrang == 3 ) then
                    if( msymoo(k_s,jh_s,k_i,jh_i) == 0 .and. msymooi(k_s,jh_s,k_i,jh_i) == 0 ) cycle
                    if( sum( abs( secoo(k_s,jh_s,k_i,jh_i,:,mpirank) ) ) > 1.e-15_db ) cycle
                  endif

! To take into account the relativistic transition channel (just for E1E1).
                  ispfg = 0
                  boucle_dip_pos: do isp1 = 1,2
                  do isp2 = 1,2
                  do isp3 = 1,2
                  do isp4 = 1,2
                    if( lseuil == 0 .and. isp1 /= isp4 ) cycle  ! at K edge, core states are mono-spin. 
                    ispfg = ispfg + 1
                    if( ispfg > 1 .and. ( ( irang /= 1 .or. jrang /= 1 ) .or. .not. dip_rel ) ) exit boucle_dip_pos 

                    Tens(:) = (0._db,0._db)
                    Tens_m(:) = (0._db,0._db)

                    if( icheck > 1 ) write(3,140) k_s, j_s, h_s, k_i, j_i, h_i, ispfg

! Loops over spherical components of the tensors
                  lm_s = 0
                  do l_s = 0,max(irang,1)
                    do m_s = -l_s,l_s
                      lm_s = lm_s + 1

                      clm_s = clm(irang,k_s,j_s,h_s,lm_s)
                      if( abs(clm_s) < eps10 ) cycle

                      lm_i = 0
                      do l_i = 0,max(jrang,1)
                        do m_i = -l_i,l_i
                          lm_i = lm_i + 1

                          clm_i = clm(jrang,k_i,j_i,h_i,lm_i)
                          if( abs(clm_i) < eps10 ) cycle

! Calcul de la composante du tenseur
                          if( sum( abs( Tens_lm(lm_s,lm_i,ispfg,:) ) ) < 1.e-15_db ) then

                            if( icheck > 1 ) write(3,150) l_s, m_s, l_i, m_i, clm_s, clm_i
                            
                            if( Final_optic ) then
                              call tens_op(Core_resolved,Final_tddft,icheck,ip_max,ip0,irang,jrang,l_s,m_s,l_i,m_i,lmax,lmoins1, &
                                lplus1,M_depend,ns_dipmag,ndim2,ninitr,nlm1g,nlm2g,nlm_probe,nlm_p_fp,nspinp,nspino, &
                                roff_rr,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)
                            else
                              call tens_ab_cluster(coef_g,Core_resolved,Dip_rel,FDM_comp,Final_tddft,Full_potential,Green, &
                                Green_int,icheck,ip_max, &
                                ip0,irang,is_g,isp1,isp2,isp3,isp4,jrang,l_s,m_s,l_i,m_g,m_i,lmax,lmoins1,lplus1,lseuil,nbseuil, &
                                ns_dipmag,ndim2,ninit1,ninit,ninitv,ninitr,nlm_probe,nlm_p_fp,NRIXS,nspinp,nspino,Rad_ioRIoi,rof, &
                                Singul,Solsing,Solsing_only,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)
                              if( numat == 1) then
                           ! For hydrogen there is only 1 core state
                                Ten(:) = Ten(:) / 2
                                Ten_m(:) = Ten_m(:) / 2
                              endif
                            endif

                            Tens_lm(lm_s,lm_i,ispfg,:) = Ten(:)
                            Tens_lm_m(lm_s,lm_i,ispfg,:) = Ten_m(:)

                          endif

                          Tens(:) = Tens(:) + clm_s * clm_i * Tens_lm(lm_s,lm_i,ispfg,:)
                          Tens_m(:) = Tens_m(:) + clm_s * clm_i * Tens_lm_m(lm_s,lm_i,ispfg,:)
                        end do
                      end do
                    end do
                  end do ! fin boucle l_s

! Remplissage de la valeur calculee dans tous les elements du tenseur
! equivalent par symetrie.

                  Tensr(:) = real( Tens(:),db )
                  Tensi(:) = aimag( Tens(:) )

! M1-M1 (dipole magnetique - dipole magnetique)
                  if( irang == 0 .and. jrang == 0 ) then
                    secmm(k_s,k_i,:,mpirank) = Tens(:)
                    if( Green_int ) then
                      secmm(k_i,k_s,:,mpirank) = Tens(:)
                      secmm_m(k_s,k_i,:,mpirank) = Tens_m(:)
                      secmm_m(k_i,k_s,:,mpirank) = - Tens_m(:)
                    else
                      secmm(k_i,k_s,:,mpirank) = conjg( Tens(:) )
                    endif
                  endif

! M1-E1 (dipole magnetique - dipole electrique)
                  if( irang == 0 .and. jrang == 1 ) then
                    secmd(k_s,k_i,:,mpirank) = Tens(:)
                    if( Green_int ) secmd_m(k_s,k_i,:,mpirank) = Tens_m(:)
                  endif

! E1-E1 (Dipole-dipole)
                  if( irang == 1 .and. jrang == 1 ) then

                    do kk_s = 1,3
                      do kk_i = kk_s,3
                        if( abs(msymdd(kk_s,kk_i)) /= abs(msymdd(k_s,k_i)) .or. abs(msymddi(kk_s,kk_i)) &
                           /= abs(msymddi(k_s,k_i)) ) cycle
                        if( msymdd(kk_s,kk_i) /= 0 ) then
                          is = msymdd(kk_s,kk_i) / msymdd(k_s,k_i)
                        else
                          is = 0
                        endif
                        if( msymddi(kk_s,kk_i) /= 0 ) then
                          isi = msymddi(kk_s,kk_i) / msymddi(k_s,k_i)
                        else
                          isi = 0
                        endif
                        if( Green_int ) then
                          secdd(kk_s,kk_i,ispfg,:,mpirank) = is * Tens(:)
                          secdd(kk_i,kk_s,ispfg,:,mpirank) = is * Tens(:)
                          secdd_m(kk_s,kk_i,ispfg,:,mpirank) = isi * Tens_m(:)
                          secdd_m(kk_i,kk_s,ispfg,:,mpirank) = -isi * Tens_m(:)
                        else
! tenseur hermitique
                          Te(:) = cmplx( is*Tensr(:), isi*Tensi(:), db)
                          secdd(kk_s,kk_i,ispfg,:,mpirank) = Te(:)
                          secdd(kk_i,kk_s,ispfg,:,mpirank) = conjg( Te(:) )
                        endif
                      end do
                    end do

! Dipole-quadrupole
                  elseif( irang == 1 .and. jrang == 2 ) then
                    do kk_s = 1,3
                      do kk_i = 1,3
                        do jj_i = kk_i,3
                          if( abs(msymdq(kk_s,kk_i,jj_i)) /= abs(msymdq(k_s,k_i,j_i)) .or. abs(msymdqi(kk_s,kk_i,jj_i)) &
                             /= abs(msymdqi(k_s,k_i,j_i)) ) cycle
                          if( msymdq(kk_s,kk_i,jj_i) /= 0 ) then
                            is = msymdq(kk_s,kk_i,jj_i) / msymdq(k_s,k_i,j_i)
                          else
                            is = 0
                          endif
                          if( msymdqi(kk_s,kk_i,jj_i) /= 0 ) then
                            isi = msymdqi(kk_s,kk_i,jj_i) / msymdqi(k_s,k_i,j_i)
                          else
                            isi = 0
                          endif
                          if( Green_int ) then
                            secdq(kk_s,kk_i,jj_i,:,mpirank) =is*Tens(:)
                            secdq(kk_s,jj_i,kk_i,:,mpirank) =is*Tens(:)
                            secdq_m(kk_s,kk_i,jj_i,:,mpirank) = isi*Tens_m(:)
                            secdq_m(kk_s,jj_i,kk_i,:,mpirank) = isi*Tens_m(:)
                          else
                            Te(:)=cmplx(is*Tensr(:),isi*Tensi(:),db)
                            secdq(kk_s,kk_i,jj_i,:,mpirank) = Te(:)
                            secdq(kk_s,jj_i,kk_i,:,mpirank) = Te(:)
                          endif
                        end do
                      end do
                    end do

! Dipole-octupole
                  elseif( irang == 1 .and. jrang == 3 ) then
                    do kk_s = 1,3
                      do kk_i = 1,3
                        do jj_i = 1,3
                          do hh_i = 1,3
                            if( abs( msymdo(kk_s,kk_i,jj_i,hh_i) ) /= abs( msymdo(k_s,k_i,j_i,h_i) )  &
                              .or. abs( msymdoi(kk_s,kk_i,jj_i,hh_i) ) /= abs( msymdoi(k_s,k_i,j_i,h_i) ) ) cycle
                            if( msymdo(kk_s,kk_i,jj_i,hh_i) /= 0 ) then
                              is = msymdo(kk_s,kk_i,jj_i,hh_i) / msymdo(k_s,k_i,j_i,h_i)
                            else
                              is = 0
                            endif
                            if( msymdoi(kk_s,kk_i,jj_i,hh_i)/=0) then
                              isi = msymdoi(kk_s,kk_i,jj_i,hh_i) / msymdoi(k_s,k_i,j_i,h_i)
                            else
                              isi = 0
                            endif
                            if( Green_int ) then
                              secdo(kk_s,kk_i,jj_i,hh_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,kk_i,hh_i,jj_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,jj_i,kk_i,hh_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,hh_i,kk_i,jj_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,jj_i,hh_i,kk_i,:,mpirank) = is * Tens(:)
                              secdo(kk_s,hh_i,jj_i,kk_i,:,mpirank) = is * Tens(:)
                              secdo_m(kk_s,kk_i,jj_i,hh_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,kk_i,hh_i,jj_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,jj_i,kk_i,hh_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,hh_i,kk_i,jj_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,jj_i,hh_i,kk_i,:,mpirank) = isi * Tens_m(:)
                              secdo_m(kk_s,hh_i,jj_i,kk_i,:,mpirank) = isi * Tens_m(:)
                            else
                              Te(:) = cmplx(is*Tensr(:),isi*Tensi(:),db)
                              secdo(kk_s,kk_i,jj_i,hh_i,:,mpirank) = Te(:)
                              secdo(kk_s,kk_i,hh_i,jj_i,:,mpirank) = Te(:)
                              secdo(kk_s,jj_i,kk_i,hh_i,:,mpirank) = Te(:)
                              secdo(kk_s,hh_i,kk_i,jj_i,:,mpirank) = Te(:)
                              secdo(kk_s,jj_i,hh_i,kk_i,:,mpirank) = Te(:)
                              secdo(kk_s,hh_i,jj_i,kk_i,:,mpirank) = Te(:)
                            endif
                          end do
                        end do
                      end do
                    end do

! Quadrupole-quadrupole
                  elseif( irang == 2 .and. jrang == 2 ) then
                    do kk_s = 1,3
                      do jj_s = 1,3
                        do kk_i = 1,3
                          do jj_i = 1,3
                            if( sum( abs( secqq(kk_s,jj_s,kk_i,jj_i,:,mpirank) ) ) > 1.e-15_db ) cycle
                            if( abs( msymqq(kk_s,jj_s,kk_i,jj_i) ) /= abs( msymqq(k_s,j_s,k_i,j_i) ) &
                                .or. abs( msymqqi(kk_s,jj_s,kk_i,jj_i) ) /= abs( msymqqi(k_s,j_s,k_i,j_i) ) ) cycle
                            if( msymqq(kk_s,jj_s,kk_i,jj_i) /= 0 ) then
                              is = msymqq(kk_s,jj_s,kk_i,jj_i) / msymqq(k_s,j_s,k_i,j_i)
                            else
                              is = 0
                            endif
                            if( msymqqi(kk_s,jj_s,kk_i,jj_i) /= 0 ) then
                              isi = msymqqi(kk_s,jj_s,kk_i,jj_i) / msymqqi(k_s,j_s,k_i,j_i)
                            else
                              isi = 0
                            endif
                            if( Green_int ) then
                              Te(:) = is * Tens(:)
                              secqq(kk_s,jj_s,kk_i,jj_i,:,mpirank)=Te(:)
                              secqq(jj_s,kk_s,kk_i,jj_i,:,mpirank)=Te(:)
                              secqq(kk_s,jj_s,jj_i,kk_i,:,mpirank)=Te(:)
                              secqq(jj_s,kk_s,jj_i,kk_i,:,mpirank)=Te(:)
                              secqq(kk_i,jj_i,kk_s,jj_s,:,mpirank)=Te(:)
                              secqq(jj_i,kk_i,kk_s,jj_s,:,mpirank)=Te(:)
                              secqq(kk_i,jj_i,jj_s,kk_s,:,mpirank)=Te(:)
                              secqq(jj_i,kk_i,jj_s,kk_s,:,mpirank)=Te(:)
                              Te(:) = isi * Tens_m(:)
                              secqq_m(kk_s,jj_s,kk_i,jj_i,:,mpirank) = Te(:)
                              secqq_m(jj_s,kk_s,kk_i,jj_i,:,mpirank) = Te(:)
                              secqq_m(kk_s,jj_s,jj_i,kk_i,:,mpirank) = Te(:)
                              secqq_m(jj_s,kk_s,jj_i,kk_i,:,mpirank) = Te(:)
                              secqq_m(kk_i,jj_i,kk_s,jj_s,:,mpirank) = - Te(:)
                              secqq_m(jj_i,kk_i,kk_s,jj_s,:,mpirank) = - Te(:)
                              secqq_m(kk_i,jj_i,jj_s,kk_s,:,mpirank) = - Te(:)
                              secqq_m(jj_i,kk_i,jj_s,kk_s,:,mpirank) = - Te(:)
                            else
                              Te(:)= cmplx( is*Tensr(:), isi*Tensi(:),db)
                              secqq(kk_s,jj_s,kk_i,jj_i,:,mpirank) = Te(:)
                              secqq(jj_s,kk_s,kk_i,jj_i,:,mpirank) = Te(:)
                              secqq(kk_s,jj_s,jj_i,kk_i,:,mpirank) = Te(:)
                              secqq(jj_s,kk_s,jj_i,kk_i,:,mpirank) = Te(:)
                              Te(:) = Conjg( Te(:) )
                              secqq(kk_i,jj_i,kk_s,jj_s,:,mpirank) = Te(:)
                              secqq(jj_i,kk_i,kk_s,jj_s,:,mpirank) = Te(:)
                              secqq(kk_i,jj_i,jj_s,kk_s,:,mpirank) = Te(:)
                              secqq(jj_i,kk_i,jj_s,kk_s,:,mpirank) = Te(:)
                            endif
                          end do
                        end do
                      end do
                    end do
! Octupole-octupole
                  elseif( irang == 3 .and. jrang == 3 ) then

                    do kk_s = 1,3
                      do jj_s = 1,3
                        do hh_s = 1,3
                          jjh_s = 3 * ( jj_s - 1 ) + hh_s
                          do kk_i = 1,3
                            do jj_i = 1,3
                              do hh_i = 1,3
                                jjh_i = 3 * ( jj_i - 1 ) + hh_i
                                i = msymoo(kk_s,jjh_s,kk_i,jjh_i)
                                j = msymoo(k_s,jh_s,k_i,jh_i)
                                if( abs(i) /= abs(j) ) cycle
                                if( j /= 0 ) then
                                  is = i / j
                                else
                                  is = 0
                                endif
                                i = msymooi(kk_s,jjh_s,kk_i,jjh_i)
                                j = msymooi(k_s,jh_s,k_i,jh_i)
                                if( abs(i) /= abs(j) ) cycle
                                if( j /= 0 ) then
                                  isi = i / j
                                else
                                  isi = 0
                                endif
                                if( Green_int ) then
                                  Te(:) = is * Tens(:)
                                  Te_m(:) = isi * Tens_m(:)
                                else
                                  Te(:)= cmplx( is * Tensr(:), isi * Tensi(:), db )
                                  Te_m(:) = Conjg( Te(:) )
                                endif

                                secoo(kk_s,jjh_s,kk_i,jjh_i,:,mpirank) = Te(:)
                                if( Green_int ) secoo_m(kk_s,jjh_s,kk_i,jjh_i,:,mpirank) = Te_m(:)

                              end do
                            end do
                          end do
                        end do
                      end do
                    end do

                  endif

                  end do
                  end do
                  end do 
                  end do boucle_dip_pos

                end do
              end do
            end do
          end do
        end do
      end do

      if( icheck < 2 ) cycle

      if( .not. E1E2 .and. irang == 1 .and. jrang == 2 ) cycle
      if( jrang == 3 .and. irang == 2 ) cycle
      if( irang == 0 .and. jrang > 1 ) cycle
      write(3,160)
      lm_s = 0
      do l_s = 0,max(irang,1)
        do m_s = -l_s,l_s
         lm_s = lm_s + 1
         lm_i = 0
         do l_i = 0,max(jrang,1)
            do m_i = -l_i,l_i
              lm_i = lm_i + 1
              do ispfg = 1,n_rel
                if( ( irang /= 1 .or. jrang /= 1 ) .and. ispfg > 1 ) exit
                if( sum( abs(Tens_lm(lm_s,lm_i,ispfg,:)) ) > eps15 ) write(3,170) irang, jrang, l_s, m_s, l_i, m_i, &
                    ispfg, ( Tens_lm(lm_s,lm_i,ispfg,initr), initr = 1,ninitr )
              end do
              if( .not. Green_int ) cycle
              do ispfg = 1,n_rel
                if( ( irang /= 1 .or. jrang /= 1 ) .and. ispfg > 1 ) exit
                if( sum( abs(Tens_lm_m(lm_s,lm_i,ispfg,:)) ) > eps15 ) write(3,180) irang, jrang, l_s, m_s, l_i, m_i, &
                   ispfg, ( Tens_lm_m(lm_s,lm_i,ispfg,initr), initr = 1,ninitr )
              end do
            end do
          end do
        end do
      end do

    end do
  end do

! Rotation to get the tensors in R1 basis

    rot_tem = matmul( rot_int, transpose(rot_atom_abs) )

    do initr = 1,ninitr

      if( E1E1 ) then
        do ispfg = 1,n_rel
          mat2(:,:) = secdd(:,:,ispfg,initr,mpirank)
          call rot_tensor_2( mat2, rot_tem )
          secdd(:,:,ispfg,initr,mpirank) = mat2(:,:)
        end do
      endif

      if( E1E2 ) then
        mat3(:,:,:) = secdq(:,:,:,initr,mpirank)
        call rot_tensor_3( mat3, rot_tem )
        secdq(:,:,:,initr,mpirank) = mat3(:,:,:)
      endif

      if( E2E2 ) then
        mat4(:,:,:,:) = secqq(:,:,:,:,initr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secqq(:,:,:,:,initr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1E3 ) then
        mat4(:,:,:,:) = secdo(:,:,:,:,initr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secdo(:,:,:,:,initr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1M1 ) then
        mat2(:,:) = secmd(:,:,initr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmd(:,:,initr,mpirank) = mat2(:,:)
      endif

      if( M1M1 ) then
        mat2(:,:) = secmm(:,:,initr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmm(:,:,initr,mpirank) = mat2(:,:)
      endif

      if( E3E3 ) then
        jh_s = 0
        do j_s = 1,3
          do h_s = 1,3
            jh_s = jh_s + 1
            jh_i = 0
            do j_i = 1,3
              do h_i = 1,3
                jh_i = jh_i + 1
                Mat6(:,j_s,h_s,:,j_i,h_i) = secoo(:,jh_s,:,jh_i,initr,mpirank)
              end do
            end do
          end do
        end do
        call rot_tensor_6( Mat6, Rot_tem )
        jh_s = 0
        do j_s = 1,3
          do h_s = 1,3
            jh_s = jh_s + 1
            jh_i = 0
            do j_i = 1,3
              do h_i = 1,3
                jh_i = jh_i + 1
                secoo(:,jh_s,:,jh_i,initr,mpirank) = Mat6(:,j_s,h_s,:,j_i,h_i)
              end do
            end do
          end do
        end do
      endif

      if( .not. Green_int ) Cycle

      if( E1E1 ) then
        do ispfg = 1,n_rel
          mat2(:,:) = secdd_m(:,:,ispfg,initr,mpirank)
          call rot_tensor_2( mat2, rot_tem )
          secdd_m(:,:,ispfg,initr,mpirank) = mat2(:,:)
        end do
      endif

      if( E1E2 ) then
        mat3(:,:,:) = secdq_m(:,:,:,initr,mpirank)
        call rot_tensor_3( mat3, rot_tem )
        secdq_m(:,:,:,initr,mpirank) = mat3(:,:,:)
      endif

      if( E2E2 ) then
        mat4(:,:,:,:) = secqq_m(:,:,:,:,initr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secqq_m(:,:,:,:,initr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1E3 ) then
        mat4(:,:,:,:) = secdo_m(:,:,:,:,initr,mpirank)
        call rot_tensor_4( mat4, rot_tem )
        secdo_m(:,:,:,:,initr,mpirank) = mat4(:,:,:,:)
      endif

      if( E1M1 ) then
        mat2(:,:) = secmd_m(:,:,initr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmd_m(:,:,initr,mpirank) = mat2(:,:)
      endif

      if( M1M1 ) then
        mat2(:,:) = secmm_m(:,:,initr,mpirank)
        call rot_tensor_2( mat2, rot_tem )
        secmm_m(:,:,initr,mpirank) = mat2(:,:)
      endif

      if( E3E3 ) then
        jh_s = 0
        do j_s = 1,3
          do h_s = 1,3
            jh_s = jh_s + 1
            jh_i = 0
            do j_i = 1,3
              do h_i = 1,3
                jh_i = jh_i + 1
                Mat6(:,j_s,h_s,:,j_i,h_i) = secoo_m(:,jh_s,:,jh_i,initr,mpirank)
              end do
            end do
          end do
        end do
        call rot_tensor_6( Mat6, Rot_tem )
        jh_s = 0
        do j_s = 1,3
          do h_s = 1,3
            jh_s = jh_s + 1
            jh_i = 0
            do j_i = 1,3
              do h_i = 1,3
                jh_i = jh_i + 1
                secoo_m(:,jh_s,:,jh_i,initr,mpirank) = Mat6(:,j_s,h_s,:,j_i,h_i)
              end do
            end do
          end do
        end do
      endif

    end do

  if( Final_optic ) then
    deallocate( roff_rr )
  else
    deallocate( rof, Singul )
  endif
  
  return
  100 format(/' ---- Tenseur_car Radial ---------',100('-'))
  110 format(/' ---- Tenseur_car Tens_ab --------',100('-'))
  120 format(/' lseuil =',i2)
  130 format(/' -- irang =',i2,', jrang =',i2,' --')
  140 format(/' k_s, j_s, h_s =',3i3,',  k_i, j_i, h_i =',3i3,',  ispfg =',i2)
  150 format(/' l_s, m_s =',2i3,',  l_i, m_i =',2i3,', clm_s, clm_i =',1p, 2e13.5)
  160 format(/' Tensor by harmonics (basis R4) :'/, ' ir jr  l_s m_s  l_i m_i ispfg     Tens(Ylm,i=1,ninitr)')
  170 format(2i3,i4,i3,2x,i4,i3,i5,1p,20e15.7)
  180 format('i',i2,i3,i4,i3,2x,i4,i3,1p,20e15.7)

end

!***********************************************************************

subroutine tens_ab_cluster(coef_g,Core_resolved,Dip_rel,FDM_comp,Final_tddft,Full_potential,Green, &
                              Green_int,icheck,ip_max, &
                              ip0,irang,is_g,isp1,isp2,isp3,isp4,jrang,l_s,m_s,l_i,m_g,m_i,lmax,lmoins1,lplus1,lseuil,nbseuil, &
                              ns_dipmag,ndim2,ninit1,ninit,ninitv,ninitr,nlm_probe,nlm_p_fp,NRIXS,nspinp,nspino,Rad_ioRIoi,rof, &
                              Singul,Solsing,Solsing_only,Spinorbite,Taull,Ten,Ten_m,Ylm_comp)

  use declarations
  implicit none

  integer, intent(in):: icheck, ip_max, ip0, irang, jrang, l_i, l_s, m_i, m_s, lseuil, nbseuil, ndim2, ninit1, &
      ninit, ninitv, ninitr, nlm_probe, nlm_p_fp, ns_dipmag, nspinp, nspino
  integer, dimension(ninit,2), intent(in):: m_g
  integer, dimension(ninit), intent(in):: is_g

  integer:: i_g_1, i_g_2, initl1, initl2, initr, is_dipmag, is_r1, is_r2, iseuil1, iseuil2, iso1, iso2, &
    isp1, isp2, isp3, isp4, ispf1, ispf2, ispinf1, ispinf2, isping1, isping2, &
    ispp_f1, ispp_f2, l_f1, l_f2, li, lm01, lm02, lm_f1, lm_f2, lmax, lmp01, &
    lmp02, lmp_f1, lmp_f2, lms_f1, lms_f2, lp_f1, lp_f2, m_f1, m_f2, mi1, mi2, mp_f1, mp_f2, mv1, mv2

  complex(kind=db):: cfe, cfs, Cg, rof_1, rof_2, Singul_e, Singul_s, Tau_rad, Tau_rad_i

  complex(kind=db):: Gaunt_i, Gaunt_nrixs, Gaunt_s, Gaunt_xrc, Gauntmag 
  complex(kind=db), dimension(ninitr):: Ten, Ten_m
  complex(kind=db), dimension(nlm_probe,nlm_p_fp,nspinp,nspino,ip0:ip_max,ninitv):: rof
  complex(kind=db), dimension(nlm_probe*nspino,nlm_probe*nspino,ndim2,ndim2,2,2,ns_dipmag):: Taull
  complex(kind=db), dimension(nlm_probe,nspinp,nlm_probe,nspinp,ip0:ip_max,ip0:ip_max,ninitv):: Singul
  complex(kind=db), dimension(nlm_probe,nlm_probe,nspinp,ip_max-ip0+1,ip_max-ip0+1,nbseuil):: Rad_ioRIoi

  logical:: Core_resolved, Dip_rel, FDM_comp, Final_tddft, Full_potential, Green, Green_int, lmoins1, lplus1, NRIXS, &
            Solsing, Solsing_only, Spinorbite, Titre, Ylm_comp

  real(kind=db):: Ci_1, Ci_2, Ci2, J_initl1, J_initl2, Jz1, Jz2
  real(kind=db), dimension(ninit,2):: coef_g
  
  li = lseuil

  Ten(:) = (0._db,0._db)
  Ten_m(:) = (0._db,0._db)

! Loop over core states < g_1 ...  > < ... >
  do i_g_1 = 1,ninit

    if( Final_tddft ) then
      initl1 = i_g_1
    else
      initl1 = 1
    endif

    J_initl1 = li + 0.5_db * is_g(i_g_1)
    if( i_g_1 <= ninit1 ) then
      Jz1 = - J_initl1 + i_g_1 - 1
    else
      Jz1 = - J_initl1 + i_g_1 - ninit1 - 1
    endif

    if( i_g_1 <= ninit1 ) then
      iseuil1 = 1
    else
      iseuil1 = 2
    endif

    if( Final_tddft ) then
      initr = 1
    elseif( Core_resolved ) then
      initr = i_g_1
    else
      initr = iseuil1
    endif

    if( ninitv == ninit .and. Final_tddft ) then
      is_r1 = i_g_1
    else
      is_r1 = iseuil1
    endif

! Loop over core state spin < g_1 ...  > < ... >
    do isping1 = 1,2 
    
      if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp1 /= isping1 ) cycle

      mi1 = m_g(i_g_1,isping1)
      Ci_1 = Coef_g(i_g_1,isping1)

      if( abs( Ci_1 ) < eps6 ) cycle

      do i_g_2 = 1,ninit

        if( .not. Final_tddft .and. i_g_1 /= i_g_2 ) cycle

        if( Final_tddft ) then
          initl2 = i_g_2
        else
          initl2 = 1
        endif

        J_initl2 = li + 0.5_db * is_g(i_g_2)
        if( i_g_2 <= ninit1 ) then
          Jz2 = - J_initl2 + i_g_2 - 1
        else
          Jz2 = - J_initl2 + i_g_2 - ninit1 - 1
        endif

        if( i_g_2 <= ninit1 ) then
          iseuil2 = 1
        else
          iseuil2 = 2
        endif

        if( ninitv == ninit .and. Final_tddft ) then
          is_r2 = i_g_2
        else
          is_r2 = iseuil2
        endif

! Loop over core state spin < ...  > < ... g_2 >
        do isping2 = 1,2  
        
!          if( .not. Final_tddft .and. isping1 /= isping2 ) cycle
          if( .not. ( Final_tddft .or. li /= 0 ) .and. isping1 /= isping2 ) cycle
          if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp4 /= isping2 ) cycle

          mi2 = m_g(i_g_2,isping2)
          Ci_2 = Coef_g(i_g_2,isping2)

          if( abs( Ci_2 ) < eps6 ) cycle

          Ci2 = Ci_1 * Ci_2

         if( icheck > 1 ) then
           write(3,110)
           write(3,120) i_g_1, isping1, Jz1, mi1, Ci_1, i_g_2, isping2, Jz2, mi2, Ci_2
         endif
         Titre = .true.

! Loop over spherical harmonics of final states < ... f_1 > < ... >
          do l_f1 = 0,lmax

            if( lplus1 .and. l_f1 < li + 1 ) cycle
            if( lmoins1 .and. l_f1 > li ) cycle
            if( l_f1 > li + irang .or. l_f1 < li - irang .or. mod(l_f1,2) /= mod(li+irang,2) ) cycle

            lm01 = l_f1**2 + l_f1 + 1
            do m_f1 = -l_f1,l_f1
              lm_f1 = lm01 + m_f1
              if( lm_f1 > nlm_probe ) cycle

! Loop over spin of final states < ... f_1 > < ... >
              do ispinf1 = 1,2

                if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp2 /= ispinf1 ) cycle

                ispf1 = min( ispinf1, nspinp )

! Il n'y a que pour l_s dipole magnetique qu'on peut avoir du spin-flip (ou pour l_s terme dipolaire relativiste)
                if( ispinf1 /= isping1 .and. ( NRIXS .or. irang > 1 .or. ( irang == 1 .and. &
                                                                           ( .not. Dip_rel .or. jrang /= 1 ) ) ) ) cycle

                if( NRIXS ) then
                  Gaunt_s = Gaunt_nrixs(l_f1,m_f1,l_s,m_s,li,mi1,Ylm_comp)
                elseif( irang == 0 ) then
                  Gaunt_s = Gauntmag(l_f1,m_f1,ispinf1,m_s,li,mi1,isping1,Ylm_comp,.true.)
                else
                  Gaunt_s = Gaunt_xrc(l_f1,m_f1,l_s,m_s,li,mi1,Ylm_comp)
                endif
                if( abs(Gaunt_s) < eps10 ) cycle

! Loop over spin of final states < ... > < f_2 ... >
                do l_f2 = 0,lmax
                  if( lplus1 .and. l_f2 < li + 1 ) cycle
                  if( lmoins1 .and. l_f2 > li ) cycle
                  if( l_f2 > li + jrang .or. l_f2 < li - jrang .or. mod(l_f2,2) /= mod(li+jrang,2) ) cycle

                  lm02 = l_f2**2 + l_f2 + 1
                  do m_f2 = -l_f2,l_f2
                    lm_f2 = lm02 + m_f2
                    if( lm_f2 > nlm_probe ) cycle

                    do ispinf2 = 1,2  ! spin etat final en sortie
                      ispf2 = min( ispinf2, nspinp )

                      if( irang == 1 .and. jrang == 1 .and. Dip_rel .and. isp3 /= ispinf2 ) cycle

                      if( ispinf2 /= isping2 .and. ( NRIXS .or. jrang > 1 .or. ( jrang == 1 &
                                                                             .and. ( .not. Dip_rel .or. irang /= 1 ) ) ) ) cycle 

                      if( .not. ( Final_tddft .or. Spinorbite ) .and. ispinf2 /= ispinf1 ) cycle

                      if( NRIXS ) then
                        Gaunt_i = Gaunt_nrixs(l_f2,m_f2,l_i,m_i,li,mi2,Ylm_comp)
                      elseif( jrang == 0 ) then
                        Gaunt_i = Gauntmag(l_f2,m_f2,ispinf2,m_i,li,mi2,isping2,Ylm_comp,.true.)
                      else
                        Gaunt_i = Gaunt_xrc(l_f2,m_f2,l_i,m_i,li,mi2,Ylm_comp)
                      endif
                      if( abs(Gaunt_i) < eps10 ) cycle

                      Cg = Ci2 * conjg( Gaunt_s ) * Gaunt_i

                      if( Final_tddft .and. ispinf1 /= isping1 ) then
                        is_dipmag = 2
                      else
                        is_dipmag = 1
                      endif

! Loop over harmonics of non spherical final states < ... f_1 > < ... >
          do lp_f1 = 0,lmax
            if( .not. Full_potential .and. lp_f1 /= l_f1 ) cycle
            lmp01 = lp_f1**2 + lp_f1 + 1

            do mp_f1 = -lp_f1,lp_f1
              if( ( nlm_p_fp == 1 .or. FDM_comp ) .and. mp_f1 /= m_f1 ) cycle
              lmp_f1 = min( lmp01 + mp_f1, nlm_p_fp ) 

! Loop over spin-orbit solution index < ... f_1 > < ... >
              do ispp_f1 = 1,nspinp
                if( .not. Spinorbite .and. ispp_f1 /= ispf1 ) cycle
                iso1 = min( ispp_f1, nspino )

                if( Spinorbite .and. nlm_p_fp == 1 ) then
                  mv1 = mp_f1 - ispinf1 + iso1
                  if( mv1 > l_f1 .or. mv1 < -l_f1 ) cycle
                else
                  mv1 = mp_f1
                endif
                lms_f1 = ( lmp01 + mv1 - 1 ) * nspino + iso1

! Loop over harmonics of non spherical final states < ... > < f_2 ... >
                do lp_f2 = 0,lmax
                  if( .not. Full_potential .and. lp_f2 /= l_f2 ) cycle
                  lmp02 = lp_f2**2 + lp_f2 + 1

                  do mp_f2 = -lp_f2,lp_f2
                    if( ( nlm_p_fp == 1 .or. FDM_comp ) .and. mp_f2 /= m_f2 ) cycle
                    lmp_f2 = min( lmp02 + mp_f2, nlm_p_fp ) 

! Loop over spin-orbit solution index < ... > < f_2 ... >
                    do ispp_f2 = 1,nspinp
                      if( .not. Spinorbite .and. ispp_f2 /= ispf2 ) cycle
                      iso2 = min( ispp_f2, nspino )

                      if( Spinorbite .and. nlm_p_fp == 1 ) then
                        mv2 = mp_f2 - ispinf2 + iso2
                        if( mv2 > l_f2 .or. mv2 < -l_f2 ) cycle
                      else
                        mv2 = mp_f2
                      endif
                      lms_f2 = ( lmp02 + mv2 - 1 ) * nspino + iso2

                      if( .not. Solsing_only ) then

                        cfe = ( 0._db, 0._db ) 
                          
                        if( FDM_comp ) then
                          cfe = Rad_ioRIoi(Lm_f1,Lm_f2,ispf1,irang,jrang,is_r1)
  
                          cfs = Rad_ioRIoi(Lm_f2,Lm_f1,ispf2,jrang,irang,is_r1)
                        
                        else

                          rof_1 = rof(lm_f1,lmp_f1,ispf1,iso1,irang,is_r1)
                          rof_2 = rof(lm_f2,lmp_f2,ispf2,iso2,jrang,is_r2)
 
                          if( Green ) then
                            cfe = rof_1 * rof_2 * Taull(lms_f1,lms_f2,initl1,initl2,ispinf1,ispinf2,is_dipmag)
                            cfs = rof_2 * rof_1 * Taull(lms_f2,lms_f1,initl2,initl1,ispinf2,ispinf1,is_dipmag)
                          else
                            cfe = conjg( rof_1 ) * rof_2 * Taull(lms_f1,lms_f2,initl1,initl2,ispinf1,ispinf2,is_dipmag)
                            cfs = conjg( rof_2 ) * rof_1 * Taull(lms_f2,lms_f1,initl2,initl1,ispinf2,ispinf1,is_dipmag)
                          endif
    
                        endif
                        
                      else

                        rof_1 = ( 0._db, 0._db )
                        rof_2 = ( 0._db, 0._db )
                        cfe = ( 0._db, 0._db )
                        cfs = ( 0._db, 0._db )

                      endif
                      
                      if( Solsing .and. i_g_1 == i_g_2 .and. lp_f1 == l_f1 .and. lp_f2 == l_f2  &
                                  .and. mp_f1 == m_f1 .and. mp_f2 == m_f2 .and. iso1 == iso2 .and. ispf1 == ispf2 &
                                  .and. ( .not. Spinorbite .or. iso1 == ispf1 ) ) then
                                                                                          ! index solution is diagonal
                        Singul_e = Singul(lm_f1+mv1-m_f1,ispf1,lm_f2+mv2-m_f2,ispf2,irang,jrang,is_r1)
                        Singul_s = Singul(lm_f2+mv2-m_f2,ispf2,lm_f1+mv1-m_f1,ispf1,jrang,irang,is_r1)
                        cfe = cfe + Singul_e
                        cfs = cfs + Singul_s
                      else
                        Singul_e = ( 0._db, 0._db ) 
                        Singul_s = ( 0._db, 0._db ) 
                      endif
                      
! Here one does not divide by pi, so convenient results seems multiplied by pi, when calculating atomic form factor.
! So, the normalization by pi must not be done in coabs to get the atomic form factor in DAFS.
                      if( Green_int ) then
                        Tau_rad = 0.5_db * ( cfe + cfs )
                        Tau_rad_i = 0.5_db * (cfe - cfs)
                        Ten(initr) = Ten(initr) + Real(Cg,db) * tau_rad + img * aimag(Cg)* tau_rad_i
                        Ten_m(initr) = Ten_m(initr) + Real(Cg,db) * tau_rad_i + img * aimag(Cg) * tau_rad
                      else
                        Tau_rad = 0.5_db * ( cfe - conjg( cfs ) )
                        Ten(initr) = Ten(initr) + Cg * Tau_rad
                      endif

                      if( icheck > 1 .and. abs( Tau_rad ) > eps15 ) then
                        if( Titre ) then
                          Titre = .false.
                          if( nlm_p_fp == 1 .and. FDM_comp ) then
                            write(3,125)
                          elseif( nlm_p_fp == 1 ) then
                            write(3,130)
                          elseif( FDM_comp ) then
                            write(3,132)
                          else
                            write(3,135)
                          endif
                        endif
                        if( nlm_p_fp == 1 .and. FDM_comp ) then
                          write(3,140) l_f1, m_f1, ispinf1, iso1, l_f2, m_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, cfe, cfs, Gaunt_i, Gaunt_s
                        elseif( nlm_p_fp == 1 ) then
                          write(3,140) l_f1, m_f1, ispinf1, iso1, l_f2, m_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, cfe, cfs, Gaunt_i, Gaunt_s, &
                            rof_1, rof_2, &
                            Taull(lms_f2,lms_f1,initl2,initl1,ispp_f2,ispp_f1,is_dipmag), &
                            Taull(lms_f1,lms_f2,initl1,initl2,ispp_f1,ispp_f2,is_dipmag), &
                            Singul_e, Singul_s
                        elseif( FDM_comp ) then
                          write(3,145) l_f1, m_f1, lp_f1, mp_f1, ispinf1, iso1, l_f2, m_f2, lp_f2, mp_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, cfe, cfs, Gaunt_i, Gaunt_s
                        else
                          write(3,145) l_f1, m_f1, lp_f1, mp_f1, ispinf1, iso1, l_f2, m_f2, lp_f2, mp_f2, ispinf2, iso2, &
                            Ten(initr), Cg*Tau_rad, cfe, cfs, Gaunt_i, Gaunt_s, &
                            rof_1, rof_2, &
                            Taull(lms_f2,lms_f1,initl2,initl1,ispp_f2,ispp_f1,is_dipmag), &
                            Taull(lms_f1,lms_f2,initl1,initl2,ispp_f1,ispp_f2,is_dipmag), &
                            Singul_e, Singul_s
                        endif 
                      endif

                    end do ! end of loop over spin-orbit solution f_i ( incoming )
                  end do ! end of loop over mp_f2 non-spherical
                end do ! end of loop over lp_f2 non-spherical

              end do ! end of loop over spin-orbit solution f_s ( outgoing )
            end do ! end of loop over mp_f1 non-spherical
          end do ! end of loop over lp_f1 non-spherical


                    end do ! end of loop over ispinf2, spin of final state ( incoming )
                  end do !  end of loop over  m_f2 final state
                end do !  end of loop over  l_f2 final state

              end do  !  end of loop over spin of final state ( outgoing )
            end do  !  end of loop over  m_f1 final state
          end do !  end of loop over  l_f1 final state

        end do    !  end of loop over spin of core states ( incoming )
      end do    !  end of loop over core states ( incoming )
    end do    !  end of loop over spin of core states ( outgoing )
  end do   !  end of loop over core states  ( outgoing )

! One multiplies by "i" in order the real part of the tensor is the absorption.
  Ten(:) = img * Ten(:)
  Ten_m(:) = img * Ten_m(:)

  return
  110 format(/' ini1 isg1 Jz1 mi1  Ci1  ini2 isg2 Jz2 mi2  Ci2')
  120 format(2(2i4,f6.1,i3,f7.3))
  125 format('  l1  m1 is1 io1  l2  m2 is2 io2',15x,'Ten',27x,'dTen',24x,'Rad_ioRIoi',19x,'Rad_ioRIoi_t',23x, &
                'Gaunt_i',24x,'Gaunt_s')
  130 format('  l1  m1 is1 io1  l2  m2 is2 io2',15x,'Ten',27x,'dTen',28x,'cfe',27x,'cfs',27x, &
                'Gaunt_i',24x,'Gaunt_s',25x,'rof_s',26x,'rof_i',26x,'Tau_s',26x,'Tau_i',25x,'Singul_e',24x,'Singul_s')
  132 format('  l1  m1 lp1 mp1 is1 io1  l2  m2 lp2 mp2 is2 io2',15x,'Ten',27x,'dTen',24x,'Rad_ioRIoi',19x,'Rad_ioRIoi_t',23x, &
                'Gaunt_i',24x,'Gaunt_s')
  135 format('  l1  m1 lp1 mp1 is1 io1  l2  m2 lp2 mp2 is2 io2',15x,'Ten',27x,'dTen',28x,'cfe',27x,'cfs',27x, &
                'Gaunt_i',24x,'Gaunt_s',25x,'rof_s',26x,'rof_i',26x,'Tau_s',26x,'Tau_i',25x,'Singul_e',24x,'Singul_s')
  140 format(8i4,1p,16(1x,2e15.7))
  145 format(12i4,1p,16(1x,2e15.7))
end

!***********************************************************************

! Preparation for writting or radial function

subroutine Write_u(icheck,Lmax,nlm_E,nr_cluster,nspinp,nspino,numat_abs,r_cluster,u,Writing_case)

  use declarations
  implicit none
  
  integer:: icheck, Lmax, n, nf, ni, nlm_E, nr_cluster, nspinp, nspino, numat_abs, Wr, Writing_case

  complex(kind=db), dimension(nr_cluster,nlm_E,nlm_E,nspinp,nspino):: u
  
  logical:: Full_potential, Radial_comp

  real(kind=db), dimension(nr_cluster):: r_cluster
  
  real(kind=db):: Rmax
  real(kind=db), dimension(:), allocatable:: r
  real(kind=db), dimension(:,:,:,:,:), allocatable:: ui, ur
  
  Full_potential = .true.
  Radial_comp = .true.

  if( icheck > 2 ) then
    n = nr_cluster
  else
    n = 6
  endif
  allocate( ui(n,nlm_E,nlm_E,nspinp,nspino) )
  allocate( ur(n,nlm_E,nlm_E,nspinp,nspino) )
  allocate( r(n) )
  if( icheck > 2 ) then
    ur(1:n,:,:,:,:) = real( u(1:n,:,:,:,:), db )
    ui(1:n,:,:,:,:) = aimag( u(1:n,:,:,:,:) )
    r(1:n) = r_cluster(1:n)
  else
    ni = min( 3, n )
    nf = n - ni
    ur(1:ni,:,:,:,:) = ur(1:ni,:,:,:,:)
    ui(1:ni,:,:,:,:) = ui(1:ni,:,:,:,:)
    r(1:ni) = r_cluster(1:ni)
    ur(ni+1:n,:,:,:,:) = real( u(nr_cluster-nf+1:nr_cluster,:,:,:,:), db )
    ui(ni+1:n,:,:,:,:) = aimag( u(nr_cluster-nf+1:nr_cluster,:,:,:,:) )
    r(ni+1:n) = r_cluster(nr_cluster-nf+1:nr_cluster)
  endif
  Radial_comp = .true.
  if( Writing_case == 4 ) then
    Rmax = r_cluster(nr_cluster-2)
    Wr = 2
  else  
    Rmax = r_cluster(nr_cluster)
    Wr = Writing_case
  endif
   
  call write_ur(Full_potential,0,Lmax,nlm_E,nlm_E,n,nspinp,nspino,numat_abs,r,Radial_comp,Rmax,ui,ur,Wr)

  deallocate( r, ui, ur )

  return
end

