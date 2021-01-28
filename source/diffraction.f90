! FDMNES subroutines
! Related to resonant x-ray diffraction and surface resonant x-ray diffraction

!----------------------------------------------------------------------------------------------------------------------------

! Preparation file for the film and surface slab position versus the bulk

Subroutine Prep_Film(angxyz,angxyz_bulk,angxyz_cap,angxyz_int,angxyz_sur,axyz,axyz_bulk,axyz_cap,axyz_int,axyz_sur, &
               Bulk_roughness,Cap_shift,Cap_roughness,Cap_thickness,Delta_bulk,Delta_cap,Delta_int,Delta_film, &
               Delta_roughness_film,Delta_sur,delta_z_bottom_cap,delta_z_top_bulk,delta_z_top_cap,delta_z_top_film,Film_shift, &
               Film_roughness,Film_thickness,hkl_film,ich,Interface_shift,itype,n_atom_bulk,n_atom_cap,n_atom_int,n_atom_per, &
               n_atom_sur,n_atom_uc,ngroup,ntype,numat,posn,posn_bulk,posn_cap,Surface_ref,Surface_shift)

  use declarations
  implicit none

  integer:: i, ich, icheck, igr, n_atom_bulk, n_atom_cap, n_atom_int, n_atom_per, n_atom_sur, n_atom_uc, ngroup, ntype, &
    Z_bottom_int, Z_bottom_film, Z_bottom_sur, Z_top_film, Z_top_int, Z_top_sur

  integer, dimension(ngroup):: itype
  integer, dimension(0:ntype):: numat
  
  logical:: hkl_film

  real(kind=db):: Atom_radius, Bulk_roughness, c_cos_z, c_cos_z_b, c_cos_z_c, c_cos_z_i, c_cos_z_s, Cap_roughness, Cap_shift, &
    Cap_thickness_used, Cap_thickness, cos_z, cos_z_b, cos_z_c, cos_z_i, cos_z_s, Delta_bulk, Delta_cap, Delta_Film, Delta_int, &
    Delta_roughness_film, Delta_sur, delta_z_bottom_cap, delta_z_top_bulk, delta_z_top_cap, delta_z_top_film, &
    Interface_thickness, Film_roughness, Film_thickness, Film_Thickness_used, Surface_ref, Surface_thickness, z_max, z_pos

  real(kind=db), dimension(4):: Film_shift, Interface_shift, Surface_shift
  real(kind=db), dimension(3):: angxyz, angxyz_bulk, angxyz_cap, angxyz_int, angxyz_sur, axyz, axyz_bulk, axyz_cap, axyz_int, &
                                axyz_sur
  real(kind=db), dimension(3,n_atom_bulk):: posn_bulk
  real(kind=db), dimension(3,n_atom_uc):: posn
  real(kind=db), dimension(3,n_atom_cap):: posn_cap

  icheck = ich

  if( hkl_film ) then
    Surface_ref = axyz(1) * axyz(2) * sin( angxyz(3) )  
  else
    Surface_ref = axyz_bulk(1) * axyz_bulk(2) * sin( angxyz_bulk(3) )  
  endif
  
! Bulk

  cos_z_b = sqrt( sin( angxyz_bulk(2) )**2 - ( ( cos( angxyz_bulk(1) ) - cos( angxyz_bulk(3) ) * cos( angxyz_bulk(2) ) ) &
        / sin( angxyz_bulk(3) ) )**2 )
  c_cos_z_b = axyz_bulk(3) * cos_z_b

  z_max = -100000._db
  do igr = 1,n_atom_bulk
    z_pos = posn_bulk(3,igr) * c_cos_z_b
    z_max = max( z_max, z_pos )
  end do

  Delta_bulk = - z_max

  delta_z_top_bulk = 0._db
  if( n_atom_per == 0 .and. n_atom_sur /= 0 ) then
    delta_z_top_bulk = Surface_shift(3) / 2
  elseif( n_atom_int /= 0 ) then
    delta_z_top_bulk = Interface_shift(3) / 2
  elseif( n_atom_int == 0 .and. n_atom_per /= 0 ) then
    delta_z_top_bulk = Film_shift(3) / 2
  endif  

! Second top most atom to define the half interlayer distance when there is no surface or film atoms
  if( delta_z_top_bulk < eps10 ) then 
    delta_z_top_bulk = c_cos_z_b / 2
    do igr = 1,n_atom_bulk
      z_pos = posn(3,igr) * c_cos_z_b
      if( z_pos > z_max - eps10 ) cycle
      delta_z_top_bulk = min( delta_z_top_bulk, ( z_max - z_pos ) / 2 )
    end do
  endif

! Film

  Z_bottom_film = 0
  Z_bottom_int = 0
  Z_bottom_sur = 0
  Z_top_film = 0
  Z_top_int = 0
  Z_top_sur = 0

  cos_z = sqrt( sin( angxyz(2) )**2 - ( ( cos( angxyz(1) ) - cos( angxyz(3) ) * cos( angxyz(2) ) ) &
        / sin( angxyz(3) ) )**2 )
  c_cos_z = axyz(3) * cos_z

  if( n_atom_int > 0 ) then
    cos_z_i = sqrt( sin( angxyz_int(2) )**2 - ( ( cos( angxyz_int(1) ) - cos( angxyz_int(3) ) * cos( angxyz_int(2) ) ) &
            / sin( angxyz_int(3) ) )**2 )
    c_cos_z_i = axyz_int(3) * cos_z_i
  endif
  if( n_atom_sur > 0 ) then
    cos_z_s = sqrt( sin( angxyz_sur(2) )**2 - ( ( cos( angxyz_sur(1) ) - cos( angxyz_sur(3) ) * cos( angxyz_sur(2) ) ) &
            / sin( angxyz_sur(3) ) )**2 )
    c_cos_z_s = axyz_sur(3) * cos_z_s
  endif

  Film_thickness_used = 0._db

  do igr = 1,n_atom_per

    do i = 0,1000000
      z_pos = ( posn(3,igr) + i ) * c_cos_z
      if( z_pos > Film_thickness - eps10 ) exit
      if( z_pos < Film_thickness_used ) cycle
      Film_thickness_used = max( Film_thickness_used, z_pos )
      Z_top_film = numat(itype(igr))
    end do

  end do

  Surface_thickness = 0._db
  do igr = n_atom_per + n_atom_int + 1, n_atom_uc
    z_pos = posn(3,igr) * c_cos_z_s
    if( z_pos < Surface_thickness ) cycle
    Surface_thickness = max( Surface_thickness, z_pos )
    Z_top_sur = numat(itype(igr))
  end do

  Interface_thickness = 0._db
  do igr = n_atom_per + 1, n_atom_per + n_atom_int
    z_pos = posn(3,igr) * c_cos_z_i
    Interface_thickness = max( Interface_thickness, z_pos )
    Z_top_int = numat(itype(igr))
  end do

  do igr = 1,n_atom_uc
    if( abs( posn(3,igr) ) < eps10 ) then
      if( igr <= n_atom_per ) then
        Z_bottom_film = numat( itype(igr) )
      elseif( igr <= n_atom_per + n_atom_int ) then
        Z_bottom_int = numat( itype(igr) )
      else
        Z_bottom_sur = numat( itype(igr) )
      endif
    endif
  end do

! Search of second top most atom
  z_max = 0._db
  do igr = 1,n_atom_per
    do i = 0,1000000
      z_pos = ( posn(3,igr) + i ) * c_cos_z
      if( z_pos > Film_thickness_used - eps10 ) exit
      z_max = max( z_max, z_pos )
    end do
  end do

  if( z_max > eps10 ) then
    delta_z_top_film = 0.5_db * ( Film_thickness_used - z_max )
  elseif( Z_top_film /= 0 .and. Z_top_film /= 0 ) then
    delta_z_top_film = Atom_radius( Z_top_film )
  else
    delta_z_top_film = 0._db
  endif

  if( Film_roughness > eps10 ) then
    Delta_roughness_film = Film_thickness_used
  else
    Delta_roughness_film = 0._db
  endif

! Cap layer

  cos_z_c = sqrt( sin( angxyz_cap(2) )**2 - ( ( cos( angxyz_cap(1) ) - cos( angxyz_cap(3) ) * cos( angxyz_cap(2) ) ) &
          / sin( angxyz_cap(3) ) )**2 )

  c_cos_z_c = cos_z_c * axyz_cap(3)

  Cap_thickness_used = 0._db
  do igr = 1,n_atom_cap

    do i = 0,1000000
      z_pos = ( posn_cap(3,igr) + i ) * c_cos_z_c
      if( z_pos > Cap_thickness - eps10 ) exit
      if( z_pos < Cap_thickness_used ) cycle
      Cap_thickness_used = max( Cap_thickness_used, z_pos )
    end do

  end do

! Search of second top most atom
  z_max = -100._db
  do igr = 1,n_atom_cap
    do i = 0,1000000
      z_pos = ( posn_cap(3,igr) + i ) * c_cos_z_c
      if( z_pos > Cap_thickness_used - eps10 ) exit
      z_max = max( z_max, z_pos )
    end do
  end do

  if( z_max > - eps10 ) then
    delta_z_top_cap = 0.5_db * ( Cap_thickness_used - z_max )
  else
    delta_z_top_cap = 0._db
  endif

  delta_z_bottom_cap = 100._db
  do igr = 1,n_atom_cap
    do i = 0,1
      z_pos = ( posn_cap(3,igr) + i ) * c_cos_z_c
      if( z_pos < eps10 ) cycle
      delta_z_bottom_cap = min( delta_z_bottom_cap, z_pos )
    end do
  end do
  delta_z_bottom_cap = 0.5_db * delta_z_bottom_cap

  Delta_int = Interface_shift(3)
  Delta_Film = Delta_int + Interface_thickness + Film_shift(3)
  Delta_sur = Delta_Film + Film_thickness_used + Surface_shift(3)
  Delta_cap = Delta_sur + Surface_thickness + Cap_shift

  if( icheck > 0 ) then
    write(3,100)
    if( n_atom_per > 0 ) then
      write(3,110) Film_thickness_used * bohr, Film_shift(3) * bohr, Film_thickness * bohr, Film_roughness * bohr
      write(3,120) delta_z_top_film * bohr, Z_bottom_film, Z_top_film
    endif
    if( n_atom_int > 0 ) write(3,130) Interface_thickness * bohr, Interface_shift(3) * bohr, Z_bottom_int, Z_top_int
    if( n_atom_sur > 0 ) write(3,140) Surface_thickness * bohr, Surface_shift(3) * bohr, Z_bottom_sur, Z_top_sur

    if( n_atom_cap > 0 ) then
      write(3,150) Cap_thickness_used * bohr, Cap_shift *bohr, Cap_thickness * bohr, Cap_roughness * bohr
      write(3,160) delta_z_top_cap * bohr, delta_z_bottom_cap * bohr
    endif
    if( n_atom_bulk > 0 ) write(3,170) Delta_bulk * bohr, delta_z_top_bulk * bohr, Bulk_roughness * bohr
    if( n_atom_int > 0 ) write(3,180) Delta_int * bohr
    if( n_atom_per > 0 ) write(3,190) Delta_Film * bohr
    if( n_atom_sur > 0 ) write(3,200) Delta_sur * bohr
    if( n_atom_cap > 0 ) write(3,210) Delta_cap * bohr
    write(3,220) Surface_ref * bohr**2 
  endif

  Film_Thickness = Film_thickness_used
  Cap_thickness = Cap_thickness_used

  return
  100 format(/'----- Prep_Film ', 100('-')/ )
  110 format(' Film thickness used =',f10.5,' A, Film shift       =',f10.5,' A,     Film thickness =',f10.5, &
             ' A, Film roughness =',f10.5,' A')
  120 format(                               36x,'Delta_z_top_film =',f10.5,' A,     Z bottom film  =',i3, &
             ',               Z top film =',i3)
  130 format(' Interface thickness =',f10.5,' A, Interface shift  =',f10.5,' A,     Z bottom int   =',i3, &
             ',               Z top int =',i3)
  140 format(' Surface thickness   =',f10.5,' A, Surface shift    =',f10.5,' A,     Z bottom sur   =',i3, &
             ',               Z top sur =',i3)
  150 format(' Cap thickness used  =',f10.5,' A, Cap shift        =',f10.5,' A,     Cap thickness  =',f10.5, &
             ',    Cap roughness =',f10.5,' A')
  160 format(                               36x,'Delta_z_top_cap  =',f10.5,' A, Delta_z_bottom_cap =',f10.5,' A')
  170 format(' Delta Bulk          =',f10.5,' A, Delta_z_top_bulk =',f10.5,' A,     Bulk roughness =',f10.5,' A')
  180 format(' Delta Interface     =',f10.5,' A')
  190 format(' Delta Film          =',f10.5,' A')
  200 format(' Delta Surface       =',f10.5,' A')
  210 format(' Delta Cap           =',f10.5,' A')
  220 format(' Reference unit cell surface =',f10.5,' A^2')
end

!*********************************************************************

! The bulk atoms equivalent by symmetry in 3D are different in 2D diffraction when not at the same z. 
! Indeed to avoid the singularities when L integer, one takes into account the absorbtion inside the bulk unit cell
! This one, exp(-mu*Length_rel), is included in the Bragg Term in coabs.
! n_bulk_z_abs is the number of different z positions for the prototypical bulk atom  
! n_bulk_zc_abs is the number of bulk atoms at a specific z position
! igr_bulk_z_abs is the index of the bulk atoms at a specific z position 

subroutine cal_n_bulk_z(First_bulk_step,i,igr_bulk_z,igr_bulk_z_abs,igreq,iprabs_nonexc,itype,n_atom_bulk,n_atom_proto, &
                        n_atom_proto_bulk,n_bulk_z,n_bulk_z_abs,n_bulk_zc,n_bulk_zc_abs,n_bulk_z_max,n_bulk_z_max_abs, &
                        neqm,ngreq,ngroup,ntype,numat,posn_bulk)

  use declarations
  implicit none

  integer:: i, igr, igr_bulk, ipr, iprabs_nonexc, jgr, jgr_bulk, jpr, m, m_abs, n, n_abs, n_atom_bulk, n_atom_proto, &
            n_atom_proto_bulk, n_bulk_z, n_bulk_z_abs, n_bulk_z_max, n_bulk_z_max_abs, neqm, ngroup, ntype, Z_abs
  
  integer, dimension(ngroup):: itype
  integer, dimension(0:ntype):: numat
  integer, dimension(0:n_atom_proto):: ngreq
  integer, dimension(n_bulk_z):: n_bulk_zc
  integer, dimension(n_bulk_z_abs):: n_bulk_zc_abs
  integer, dimension(n_bulk_z_max,n_bulk_z):: igr_bulk_z
  integer, dimension(n_bulk_z_max_abs,n_bulk_z_abs):: igr_bulk_z_abs
  integer, dimension(0:n_atom_proto,neqm):: igreq

  logical:: First_bulk_step  
  logical, dimension(n_atom_bulk):: ok
  
  real(kind=db):: p_z1, p_z2

  real(kind=db), dimension(3,n_atom_bulk):: posn_bulk

  Z_abs = numat( itype( igreq(iprabs_nonexc,1) ) )
  ok(:) = .false. 
  n = 0
  n_abs = 0

  do ipr = n_atom_proto-n_atom_proto_bulk+1,n_atom_proto

! The non absorbing atom are tacken into account only during the cycle of the first absorbing atom
    if( .not. First_bulk_step .and. ipr /= iprabs_nonexc ) cycle
    
    do igr = 1,ngreq(ipr)
      igr_bulk = igreq(ipr,igr) - ngroup + n_atom_bulk 
      if( ok(igr_bulk) ) cycle
      ok(igr_bulk) = .true.

      if( ipr == iprabs_nonexc ) then 
        n_abs = n_abs + 1
        m_abs = 1 
        if( i == 3 ) igr_bulk_z_abs(m_abs,n_abs) = igr
      elseif( numat( itype(igreq(ipr,igr)) ) == Z_abs ) then
! One must not consider here the absorbing atoms which corresponds to another cycle
        cycle
      else
        n = n + 1
        m = 1 
        if( i == 3 ) igr_bulk_z(m,n) = igreq(ipr,igr)
      endif
      p_z1 = posn_bulk( 3, igr_bulk )

      do jpr = n_atom_proto-n_atom_proto_bulk+1,n_atom_proto
        if( .not. First_bulk_step .and. jpr /= iprabs_nonexc ) cycle
        if( ( ipr == iprabs_nonexc .and. jpr /= iprabs_nonexc ) .or. ( ipr /= iprabs_nonexc .and. jpr == iprabs_nonexc ) ) cycle

        do jgr = 1,ngreq(jpr)
          jgr_bulk = igreq(jpr,jgr) - ngroup + n_atom_bulk 
          if( ok(jgr_bulk) ) cycle
          p_z2 = posn_bulk(3, jgr_bulk )
          if( abs( p_z1 - p_z2 ) > eps10 ) cycle
          ok(jgr_bulk) = .true.
          if( jpr == iprabs_nonexc ) then
            m_abs = m_abs + 1 
            if( i == 3 ) igr_bulk_z_abs(m_abs,n_abs) = jgr
          elseif( numat( itype( igreq(jpr,jgr) ) ) == Z_abs ) then
! One must not consider here the absorbing atoms which corresponds to another cycle
            cycle
          else 
            m = m + 1 
            if( i == 3 ) igr_bulk_z(m,n) = igreq(jpr,jgr)
          endif
        end do
        
      end do
 
      if( i == 2 ) then
        if( ipr == iprabs_nonexc ) then 
          n_bulk_zc_abs(n_abs) = m_abs
          n_bulk_z_max_abs = max( n_bulk_z_max_abs, m_abs )
        else
          n_bulk_zc(n) = m
          n_bulk_z_max = max( n_bulk_z_max, m )
        endif
      endif 
      
    end do
    if( i == 1 ) then
      n_bulk_z = n
      n_bulk_z_abs = n_abs
    endif
  end do

  return
end

!*********************************************************************

! Calculation of polarizations, Bragg factors and non resonant scattering amplitudes for RXD and SRXRD
! Outputs : Bragg_abs : Bragg term for the absorbing atoms x temp effect x occupancy
!           Bragg_rgh_bulk_abs : same thing but due to roughness when there is absorbing atoms in the bulk
!           Phdf0t : sum of the non resonant contributions
!           phdf0t_rgh_bulk : sum of the non resonant contributions from the bulk coming from the roughness
!                             when there is an absorbing atom in the bulk (for SRXRD)
!           phdt : sum of the bragg terms over all the absorbing atoms
!           Length_abs : path length in 1 nulk unit cell (for SRXRD)
!           Length_rel : fraction of length in bulk unit cell to avoid singularities ( for SRXRD )

subroutine Prepdafs(Abs_in_bulk,Angle_or,Angle_mode,Angpoldafs,Angxyz,Angxyz_bulk,Angxyz_cap,Angxyz_int,Angxyz_sur,Axe_atom_gr, &
          axyz,axyz_bulk,axyz_cap,axyz_int,axyz_sur,Bormann,Bragg_abs,Bragg_rgh_bulk_abs,Bulk,Bulk_roughness,Bulk_step,Cap_B_iso, &
          Cap_layer,Cap_roughness,Cap_thickness,Dafs_bio,Delta_bulk,Delta_cap,Delta_film,Delta_int, &
          Delta_roughness_film,Delta_sur,delta_z_bottom_cap,delta_z_top_bulk,delta_z_top_cap,delta_z_top_film,Eseuil,f_avantseuil, &
          f_no_res,Film,Film_roughness,Film_shift,Film_thickness,hkl_dafs,hkl_film,hkl_ref,ich,igr_bulk_z,igr_bulk_z_abs, &
          igreq,Interface_shift,iprabs_nonexc,isigpi,itabs,itypepr, &
          Length_abs,Length_rel,Length_rel_abs,lvval,Magnetic,Mat_or,Mat_UB,mpirank0,n_atom_bulk,n_atom_cap, &
          n_atom_int,n_atom_per,n_atom_proto,n_atom_proto_bulk,n_atom_proto_uc,n_atom_sur,n_atom_uc,n_bulk_z,n_bulk_z_abs, &
          n_bulk_z_max,n_bulk_z_max_abs, &
          n_bulk_zc,n_bulk_zc_abs,n_max,natomsym,nbseuil,neqm,ngreq,ngrm,ngroup,ngroup_m,ngroup_taux,ngroup_temp, &
          nlat,nlatm,nphi_dafs,nphim,npl_2d,npldafs,npldafs_2d,npldafs_e,npldafs_f,nrato,nrm,nspin,ntype,numat,Operation_mode, &
          Operation_mode_used,Orthmatt,phdf0t,phdf0t_rgh_bulk,phdt,phdt_rgh_Bulk,phi_0,Poldafse,Poldafsem,Poldafss,Poldafssm, &
          popatm,posn,posn_bulk,posn_cap,psival,rato,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_shift,Taux,Taux_cap,Taux_oc, &
          Temp,Temp_coef,Temp_B_iso,Vec_orig,Vecdafse,Vecdafsem,Vecdafss,Vecdafssm,Z_bulk,Z_cap)

  use declarations
  implicit none

  integer:: i_sens, i, ich, icheck, igr, ip, ipl, ipr, iprabs_nonexc, it, itabs, iwrite, j, jgr, mpirank0, n_atom_bulk, &
    n_atom_cap, n_atom_int, n_atom_per, n_atom_proto, n_atom_proto_bulk, n_atom_proto_uc, n_atom_sur, n_atom_uc, n_bulk_z, &
    n_bulk_z_abs, n_bulk_z_max, n_bulk_z_max_abs, n_max, n1_proto, n2_proto, &
    natomsym, nbseuil, neqm, ngrm, ngroup, ngroup_m, ngroup_taux, ngroup_temp, ni, nlatm, nphim, npldafs, npldafs_2d, npldafs_e, &
    npldafs_f, nrm, nspin, ntype, Z, Z_abs

  integer, dimension(2):: Mult_bulk, Mult_film
  integer, dimension(3):: hkl_ref
  integer, dimension(n_bulk_z):: n_bulk_zc
  integer, dimension(n_bulk_z_abs):: n_bulk_zc_abs
  integer, dimension(n_bulk_z_max,n_bulk_z):: igr_bulk_z
  integer, dimension(n_bulk_z_max_abs,n_bulk_z_abs):: igr_bulk_z_abs

  character(len=4):: mot4
  character(len=64):: mot64

  complex(kind=db):: f_avantseuil, cfac, ph_cjg, Sum_T_Bragg_f, Sum_T_Bragg_rgh_Bulk_f
  complex(kind=db), dimension(3):: vec_a, vec_b, pe, ps
  complex(kind=db), dimension(npldafs):: v_vec_a, v_vec_b
  complex(kind=db), dimension(npldafs):: Truncation
  complex(kind=db), dimension(npldafs,n_max):: Sum_Bragg_abs, Sum_Bragg_bulk_abs_f0
  complex(kind=db), dimension(npldafs,n_bulk_z):: Sum_Bragg_nonabs, Sum_Bragg_nonabs_f, Sum_Bragg_nonabs_f0, Sum_Bragg_nonabs_fa
  complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdf0t_rgh_bulk, phdt_rgh_Bulk 
  complex(kind=db), dimension(npldafs,nphim,n_max):: phdt
  complex(kind=db), dimension(natomsym,npldafs):: Bragg_abs, Bragg_rgh_bulk_abs
  complex(kind=db), dimension(n_atom_proto,npldafs):: Sum_Bragg, Sum_Bragg_f0, Sum_Bragg_f0_cj, Sum_Bragg_fa, Sum_Bragg_fa_cj, &
                                                      Sum_T_Bragg_fmag
  complex(kind=db), dimension(3,npldafs_e):: Poldafsem, Poldafssm
  complex(kind=db), dimension(3,npldafs,nphim):: Poldafse, Poldafss
  complex(kind=db), dimension(:), allocatable:: phd_f_bulk, phd_f_cap, Phase_bulk
  complex(kind=db), dimension(:,:), allocatable:: Bragg, Bragg_rgh_Bulk, Sum_Bragg_rgh_Bulk, Sum_Bragg_rgh_Bulk_f0, &
                                                  Sum_Bragg_rgh_Bulk_fa
  complex(kind=db), dimension(:,:,:), allocatable:: Bragg_f_mo, Bragg_f_ms, Bragg_rgh_Bulk_f_mo, Bragg_rgh_Bulk_f_ms

  integer, dimension(0:ntype):: nlat, nrato, numat
  integer, dimension(n_atom_bulk):: Z_bulk
  integer, dimension(n_atom_cap):: Z_cap
  integer, dimension(0:ntype,nlatm):: lvval
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(npldafs_2d):: npl_2d
  integer, dimension(5,npldafs_2d):: Operation_mode
  integer, dimension(2,npldafs):: isigpi
  integer, dimension(0:n_atom_proto):: itypepr, ngreq
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(neqm):: igr_tem

  logical:: Abs_in_bulk, Bulk, Bulk_step, Bormann, Cap_layer, Dafs_bio, Debye, Film, Film_periodical, hkl_film, Interlayer, &
            Magnetic, Operation_mode_used, Surface, Taux, Temp_B_iso

  real(kind=db):: abs_cap, arg, Bulk_roughness, c_cos_z, c_cos_z_b, c_cos_z_i, c_cos_z_s, Cal_volume_maille, Cap_B_iso, &
    Cap_roughness, Cap_thickness, cos_z, cos_z_b, cos_z_i, cos_z_s, Deb, Delta_2, Delta_bulk, Delta_cap, Delta_int, Delta_film, &
    Delta_roughness_film, Delta_sur, delta_z_bottom_cap, delta_z_top_bulk, delta_z_top_cap, delta_z_top_film, dpdeg, &
    DW, Energy_photon, Film_roughness, Film_thickness, fpp_bulk_tot, &
    fpp_cap_tot, konde, Orig_bottom_film, Orig_top_bulk, Orig_top_film, p_3, phi_0,  Q_mod_A, rap_lsur2s, Roughness_erfc, Taux_r, &
    Temp, Tempt, Volume_maille, x, z_pos, z_top

  real(kind=db), dimension(2):: f_no_res
  real(kind=db), dimension(3):: angxyz, angxyz_bulk, angxyz_cap, angxyz_int, angxyz_sur, axyz, axyz_bulk, axyz_cap, axyz_int, &
     axyz_sur, hkl, p, v, Vec_orig, we, ws
  real(kind=db), dimension(4):: Film_shift, Interface_shift, Surface_shift
  real(kind=db), dimension(npldafs):: Length_abs, Q_mod
  real(kind=db), dimension(npldafs_f):: Angle_or
  real(kind=db), dimension(3,npldafs_2d):: Angle_mode
  real(kind=db), dimension(n_atom_cap):: Taux_cap
  real(kind=db), dimension(ngroup_taux):: Taux_oc
  real(kind=db), dimension(ngroup_temp):: Temp_coef
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(n_bulk_z):: Length_rel
  real(kind=db), dimension(n_bulk_z_abs):: Length_rel_abs
  real(kind=db), dimension(3,3):: Mat_bulk, Mat_bulk_i, Mat_or, Mat_UB, Orthmatt, Orthmati
  real(kind=db), dimension(3,n_atom_uc):: posn
  real(kind=db), dimension(3,n_atom_bulk):: posn_bulk
  real(kind=db), dimension(3,n_atom_cap):: posn_cap
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_g, Axe_atom_gr
  real(kind=db), dimension(3,npldafs):: angpoldafs, hkl_dafs
  real(kind=db), dimension(3,npldafs_e):: Vecdafsem, Vecdafssm
  real(kind=db), dimension(3,npldafs,nphim):: Vecdafse, Vecdafss
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(n_atom_proto):: fp, fpp
  real(kind=db), dimension(n_atom_proto,npldafs):: f_ms, f_mo, f0
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm
  
  icheck = ich
  if( icheck > 0 ) write(3,110)

  phdf0t_rgh_Bulk(:,:) = (0._db, 0._db)
  Bragg_rgh_Bulk_abs(:,:) = (0._db, 0._db)
  Sum_Bragg_abs(:,:) = (0._db,0._db)
  Sum_Bragg_nonabs(:,:) = (0._db,0._db)
  Sum_Bragg_nonabs_f(:,:) = (0._db,0._db)
  Sum_Bragg_nonabs_f0(:,:) = (0._db,0._db)
  Sum_Bragg_nonabs_fa(:,:) = (0._db,0._db)

  call invermat(Orthmatt,Orthmati)

  Mult_film(:) = 0
  Mult_bulk(:) = 0

  Z_abs = numat(itabs)

  if( Bulk_step ) then
    n1_proto = n_atom_proto - n_atom_proto_bulk + 1
    n2_proto = n_atom_proto
  else
    n1_proto = 1
    n2_proto = n_atom_proto_uc
  endif

! lambda = 2 * pi / k = 2 * d * sintheta
! In S.I. vecond = k = E*alfa*4*pi*epsilon0 / (e*e)
! In ua and rydb : k = 0.5 * alfa * E
  Energy_photon = Eseuil(nbseuil) ! exact only for single edge, that is nbseuil = 1
  konde = 0.5_db * alfa_sf * max( Energy_photon, 1._db )

  if( Operation_mode_used ) then
    call SRXD(Angle_mode,angpoldafs,angxyz,angxyz_bulk,axyz,axyz_bulk,Bulk,Bulk_step,Energy_photon,Film,hkl_dafs,hkl_film, &
                hkl_ref,icheck,isigpi,Length_abs,Mat_UB,nphim,npl_2d,npldafs,npldafs_2d,Operation_mode,Orthmati,phi_0,Poldafse, &
                Poldafss,Q_mod,Vecdafse,Vecdafss)
  else
    call Pol_dafs(Angle_or,Angpoldafs,Angxyz,Angxyz_bulk,axyz,axyz_bulk,Bormann,Bulk,Bulk_step,Dafs_bio,Energy_photon, &
            Film,hkl_dafs,hkl_film,icheck,isigpi,Length_abs,Mat_or,mpirank0,nphi_dafs,nphim,npldafs,npldafs_e,npldafs_f,Orthmatt, &
            Orthmati,Poldafse,Poldafsem,Poldafss,Poldafssm,Q_mod,Vec_orig,Vecdafse,Vecdafsem,Vecdafss,Vecdafssm)
  endif

  call Atom_structure_factor(Energy_photon, f_ms, f_mo, f_no_res, f0, itypepr, lvval, Magnetic, n_atom_proto, nlat, nlatm, &
                                 npldafs, nrato, nrm, nspin, ntype, numat, Q_mod, rato, popatm, psival)

  if( Temp > 0.00001_db ) then
    Debye = .true.
    Tempt = max( 1._db, Temp )      ! default case: T = 1 K
  else
    Debye = .false.
  endif

  if( Bulk ) then

    if( n_atom_per == 0 ) then
      call bulk_base_tr(angxyz,angxyz_bulk,axyz,axyz_bulk,Surface_shift,icheck,Mat_bulk,Mult_bulk,Mult_film)
    else
      call bulk_base_tr(angxyz,angxyz_bulk,axyz,axyz_bulk,Film_shift,icheck,Mat_bulk,Mult_bulk,Mult_film)
    endif
    if( hkl_film .and. .not. Bulk_step ) call invermat( Mat_bulk, Mat_bulk_i )

    cos_z_b = sqrt( sin( angxyz_bulk(2) )**2 - ( ( cos( angxyz_bulk(1) ) - cos( angxyz_bulk(3) ) * cos( angxyz_bulk(2) ) ) &
            / sin( angxyz_bulk(3) ) )**2 )
    c_cos_z_b = axyz_bulk(3) * cos_z_b

  endif

  cos_z = sqrt( sin( angxyz(2) )**2 - ( ( cos( angxyz(1) ) - cos( angxyz(3) ) * cos( angxyz(2) ) ) &
          / sin( angxyz(3) ) )**2 )
  c_cos_z = axyz(3) * cos_z
  if( n_atom_sur > 0 ) then
    cos_z_s = sqrt( sin( angxyz_sur(2) )**2 - ( ( cos( angxyz_sur(1) ) - cos( angxyz_sur(3) ) * cos( angxyz_sur(2) ) ) &
            / sin( angxyz_sur(3) ) )**2 )
    c_cos_z_s = axyz_sur(3) * cos_z_s
  endif
  if( n_atom_int > 0 ) then
    cos_z_i = sqrt( sin( angxyz_int(2) )**2 - ( ( cos( angxyz_int(1) ) - cos( angxyz_int(3) ) * cos( angxyz_int(2) ) ) &
            / sin( angxyz_int(3) ) )**2 )
    c_cos_z_i = axyz_int(3) * cos_z_i
  endif

  if( Bulk_step ) then

! Conversion in linear absorption coefficient in micrometer^-1 (f_avanseuil is in Mbarn)
    Volume_maille = Cal_Volume_maille(axyz_bulk,angxyz_bulk)
    if( Abs_in_bulk ) then
      fpp_bulk_tot = 0._db    ! because it is added in the convolution step
    else
      fpp_bulk_tot = 100 * aimag( f_avantseuil ) / ( Volume_maille * bohr**3 )
    endif

    z_top = -100000._db
    do igr = 1,n_atom_bulk
      z_top = max( posn_bulk(3,igr), z_top )
    end do
    if( icheck > 1 ) write(3,'(/a17,1p,e13.5,a14)') ' fpp_bulk_tot   =',fpp_bulk_tot,' micrometer^-1'
  endif

  if( icheck > 1 .and. ( Temp_B_iso .or. Debye ) ) write(3,'(/A)') ' ( h, k, l)   Z   Debye-Waller Attenuation'

  if( Magnetic ) then
    allocate( Bragg_f_mo(n1_proto:n2_proto,ngrm,npldafs) )
    allocate( Bragg_f_ms(n1_proto:n2_proto,ngrm,npldafs) )
    Bragg_f_ms(:,:,:) = (0._db,0._db)
    Bragg_f_mo(:,:,:) = (0._db,0._db)
    if( Bulk_step .and. Bulk_roughness > eps10 ) then
      allocate( Bragg_rgh_Bulk_f_mo(n1_proto:n2_proto,ngrm,npldafs) )
      allocate( Bragg_rgh_Bulk_f_ms(n1_proto:n2_proto,ngrm,npldafs) )
      Bragg_rgh_Bulk_f_ms(:,:,:) = (0._db,0._db)
      Bragg_rgh_Bulk_f_mo(:,:,:) = (0._db,0._db)
    endif
  endif
  
  if( Bulk_step .and. Bulk_roughness > eps10 ) then
    allocate( Sum_Bragg_rgh_Bulk(n_atom_proto,npldafs) )
    allocate( Sum_Bragg_rgh_Bulk_f0(n_atom_proto,npldafs) )
    allocate( Sum_Bragg_rgh_Bulk_fa(n_atom_proto,npldafs) )
  endif

  Orig_top_film = Film_thickness + delta_z_top_film
  Orig_bottom_film = - delta_z_top_bulk
  Orig_top_bulk = delta_z_top_bulk - Delta_bulk  ! in the bulk origin unit cell 

! Loop over non equivalent atom sites in the unit cell
  do ipr = n1_proto,n2_proto

! Bragg term calculation
    it = itypepr( ipr )
    if( numat( it ) == 0 ) cycle
    igr_tem(:) = igreq(ipr,:)

    allocate( Bragg(ngreq(ipr),npldafs) )
    Bragg(:,:) = cmplx( 0._db, 0._db )
    if( Bulk_step .and. Bulk_roughness > eps10 ) then
      allocate( Bragg_rgh_Bulk(ngreq(ipr),npldafs) )
      Bragg_rgh_Bulk(:,:) = cmplx( 0._db, 0._db )
    endif

! Loop over polarizations (or Q direction)
    do ipl = 1,npldafs

      hkl(:) = hkl_dafs(:,ipl)
! When Bulk_step, Mat_bulk is transformation from film to bulk
      if( ( Bulk .and. .not. Bulk_step .and. .not. hkl_film ) .or. ( Bulk_step .and. hkl_film ) ) hkl = matmul( Mat_bulk, hkl )

      if( Debye ) then
        Q_mod_A = Q_mod(ipl) / bohr   ! one needs it in angstroem^-1
        Deb = DW(Q_mod_A,numat(it),tempt)
      elseif( Temp_B_iso ) then
! Delta_2 = ( Sin(Theta_Bragg)/Lambda )^2 (lambda in Angstroem)
! Temp_coef = 8*pi^2 * <u^2> is in Angtroem^2
        Delta_2 = ( Q_mod(ipl) / ( quatre_pi * bohr ) )**2
      endif

! Loop over equivalent atoms (related by symmetry operations)
      do igr = 1,ngreq(ipr)
        jgr = igr_tem(igr)
        if( Bulk_step ) then
          p(:) = posn_bulk(:,jgr-ngroup+n_atom_bulk)
        else
          p(:) = posn(:,jgr)
        endif
        Interlayer = Film .and. .not. Bulk_step .and. jgr > n_atom_per .and. jgr <= n_atom_per + n_atom_int
        Surface = Film .and. .not. Bulk_step .and. jgr > n_atom_per + n_atom_int
        Film_periodical = Film .and. .not. Bulk_step .and. jgr <= n_atom_per

        if( Interlayer ) then
          x = ( axyz_int(1) - axyz(1) ) * hkl(1) + ( axyz_int(2) - axyz(2) ) * hkl(2)
          if( abs( x ) > eps10 ) cycle
          p(3) = p(3) * c_cos_z_i / c_cos_z
        elseif( Surface ) then
          x = ( axyz_sur(1) - axyz(1) ) * hkl(1) + ( axyz_sur(2) - axyz(2) ) * hkl(2)
          if( abs( x ) > eps10 ) cycle
          p(3) = p(3) * c_cos_z_s / c_cos_z
        endif

 ! Shift between film and bulk
        if( Film .and. .not. Bulk_step ) then
          if( n_atom_int > 0 ) p(1:2) = p(1:2) + Interface_shift(1:2) / axyz_int(1:2)
          if( n_atom_per > 0 .and. .not. Interlayer ) p(1:2) = p(1:2) +  Film_shift(1:2) / axyz(1:2)
          if( n_atom_sur > 0 .and. .not. ( Interlayer .or. Film_periodical ) ) p(1:2) = p(1:2) + Surface_shift(1:2) / axyz_sur(1:2)
        endif

        if( Film_periodical .or. ( Bulk_step .and. Bulk_roughness > eps10 ) ) then
          ni = 1000000
        else
          ni = 0
        endif

        p_3 = p(3)
        do i_sens = 1,2  ! When i_sens = 2, it is to consider the bulk_roughness
          if( i_sens == 1 ) then
            p(3) = p_3 - 1._db
          else
            p(3) = p_3
          endif
          
          do i = 0,ni
 
            if( i_sens == 1 ) then
              p(3) = p(3) + 1._db
            elseif( .not. ( Bulk_roughness > eps10 .and. ( Film_periodical .or. Bulk_step ) ) ) then
              exit
            elseif( i == 0 ) then
              cycle
            else
              p(3) = p(3) - 1._db
            endif
            
            Taux_r = 1._db
            if( Film_periodical ) then
              z_pos = p(3) * c_cos_z
              if( Film_roughness > eps10 ) then
                if( z_pos - Orig_top_film > Delta_roughness_film + eps10 ) exit
                Taux_r = Roughness_erfc(z_pos,Orig_top_film,Film_roughness)
              else
                if( z_pos > Film_thickness + eps10 ) exit
              endif
              if( Bulk_roughness > eps10 ) then
                if( z_pos - Orig_bottom_film < - 4 * Bulk_roughness ) then
                  exit
                elseif( z_pos - Orig_bottom_film < 4 * Bulk_roughness ) then
                  Taux_r = Taux_r * ( 1 - Roughness_erfc(z_pos,Orig_bottom_film,Bulk_roughness) )
                endif
              endif
            endif

            arg = deux_pi * sum( p(:) * hkl(:) )
            if( Film_periodical ) then
              arg = arg + deux_pi * ( Delta_film / c_cos_z ) * hkl(3)
            elseif( Interlayer ) then
              arg = arg + deux_pi * ( Delta_int / c_cos_z ) * hkl(3)
            elseif( Surface ) then
              arg = arg + deux_pi * ( Delta_sur / c_cos_z ) * hkl(3)
            endif

            if( .not. ( Bulk_step .and. i /= 0 ) ) &
              Bragg(igr,ipl) = Bragg(igr,ipl) + Taux_r * cmplx( cos(arg), sin(arg), db )
            
            if( Bulk_step .and. Bulk_roughness > eps10 ) then
              z_pos = p(3) * c_cos_z
              if( abs( z_pos - Orig_top_bulk ) > 4 * Bulk_roughness ) then
                if( i == 0 ) then
                  cycle
                else
                  exit
                endif
              endif
              if( z_pos > Orig_top_bulk ) then
                Taux_r = Roughness_erfc(z_pos,Orig_top_bulk,Bulk_roughness)
              else
          ! Inside the bulk, one substracts the missing weight
                Taux_r = - 1 + Roughness_erfc(z_pos,Orig_top_bulk,Bulk_roughness) 
              endif
              Bragg_rgh_Bulk(igr,ipl) = Bragg_rgh_Bulk(igr,ipl) + Taux_r * cmplx( cos(arg), sin(arg), db ) 
            endif
                         
          end do

        end do
        
        if( Interlayer ) then
          Bragg(igr,ipl) = Bragg(igr,ipl) * axyz(1) * axyz(2) / ( axyz_int(1) * axyz_int(2) )
        elseif( Surface .and. n_atom_per == 0 ) then
          Bragg(igr,ipl) = Bragg(igr,ipl) * axyz(1) * axyz(2) / ( axyz_sur(1) * axyz_sur(2) )
        endif

! Taking into account of the ppcm between bulk and film surfaces
! For the bulk, it is done in the trucature calculation
        if( Bulk .and. .not. Bulk_step ) then
          cfac = ( 0._db, 0._db )
          do i = 0,Max( Mult_film(1)-1, 0 )
          ! When film and bulk are not rationnal, non specular relection contains just the bulk 
            if( .not. hkl_film .and. Mult_film(1) == 0 .and. hkl(1) /= 0 ) cycle
            do j = 0,Max( Mult_film(2)-1, 0 )
              if( .not. hkl_film .and. Mult_film(2) == 0 .and. hkl(2) /= 0 ) cycle
              arg = deux_pi * ( i * hkl(1) + j * hkl(2) )
              cfac = cfac + cmplx( cos(arg), sin(arg), db )
            end do
          end do
          if( abs( cfac ) < eps10 ) then
            Bragg(igr,ipl) = ( 0._db, 0._db )
          else
            Bragg(igr,ipl) = cfac * Bragg(igr,ipl)
          endif
        endif

        if( Bulk_step ) then
          boucle_i: do i = 1,n_bulk_z_abs
            do j = 1,n_bulk_zc_abs(i)
              if( igr_bulk_z_abs(j,i) /= igr ) cycle
              Length_rel_abs(i) = z_top - p(3)
              exit boucle_i
            end do
          end do boucle_i
          boucle_k: do i = 1,n_bulk_z
            do j = 1,n_bulk_zc(i)
              if( igr_bulk_z(j,i) /= igreq(ipr,igr) ) cycle
              Length_rel(i) = z_top - p(3)
              exit boucle_k
            end do
          end do boucle_k
        endif
        if( Taux ) Bragg(igr,ipl) = Bragg(igr,ipl) * Taux_oc(jgr)
        if( Temp_B_iso ) Deb = exp( - Temp_coef(jgr) * Delta_2 )
        if( icheck > 1 .and. ( Debye .or. Temp_B_iso ) ) write(3,165) nint( hkl_dafs(1:3,ipl) ), numat(it), Deb
        if( Debye .or. Temp_B_iso ) Bragg(igr,ipl) = Bragg(igr,ipl) * Deb
        if( ipr == iprabs_nonexc ) Bragg_abs(igr,ipl) = Bragg(igr,ipl)

        if( Bulk_step .and. Bulk_roughness > eps10 ) then
          Bragg_rgh_Bulk(igr,ipl) = Bragg_rgh_Bulk(igr,ipl) * exp( - ( z_top - p(3) ) * fpp_bulk_tot * Length_abs(ipl) )
          if( Taux ) Bragg_rgh_Bulk(igr,ipl) = Bragg_rgh_Bulk(igr,ipl) * Taux_oc(jgr)
          if( Temp_B_iso ) Deb = exp( - Temp_coef(igr) * Delta_2 )
          if( Debye .or. Temp_B_iso ) Bragg_rgh_Bulk(igr,ipl) = Bragg_rgh_Bulk(igr,ipl) * Deb
          if( ipr == iprabs_nonexc ) Bragg_rgh_Bulk_abs(igr,ipl) = Bragg_rgh_Bulk(igr,ipl)
        endif

      end do ! end of loop over igr

      Sum_Bragg(ipr,ipl) = sum( Bragg(1:ngreq(ipr),ipl) )
      if( ipr == iprabs_nonexc ) then
        if( Bulk_step ) then
          Sum_Bragg_abs(ipl,:) = (0._db,0._db)
          do i = 1,n_bulk_z_abs
            do j = 1,n_bulk_zc_abs(i)
              igr = igr_bulk_z_abs(j,i)
              Sum_Bragg_abs(ipl,i) = Sum_Bragg_abs(ipl,i) + Bragg(igr,ipl)
            end do
          end do
        else
          Sum_Bragg_abs(ipl,1) = Sum_Bragg(ipr,ipl)
        endif
      endif

      Sum_Bragg(ipr,ipl) = sum( Bragg(1:ngreq(ipr),ipl) )

      if( Bulk_step .and. Bulk_roughness > eps10 ) &
        Sum_Bragg_rgh_Bulk(ipr,ipl) = sum( Bragg_rgh_Bulk(1:ngreq(ipr),ipl) )

    end do ! end of loop over ipl

    Z = numat( itypepr( ipr ) )
    if( Z == 0 ) cycle
! We do not consider here the anomalous part of the absorbing atom because it is calculated further in the code !
    if( Z /= Z_abs ) then
      call fprime(Z,Energy_photon,fpp(ipr),fp(ipr))
    else
      fp(:) = 0._db; fpp(ipr) = 0._db
    endif

    if( Magnetic ) then

      do ipl = 1,npldafs
        if( abs(f_ms(ipr,ipl)) < eps10 ) cycle

        do igr = 1,ngreq(ipr)
          jgr = igreq(ipr,igr)
! Sign minus are to transform to crystallographic convention.
! This is conform to what is done in routine convolution.f90.
! Bragg is already with this convention.
          Bragg_f_ms(ipr,igr,ipl) = - img * Bragg(igr,ipl) * f_ms(ipr,ipl)
          Bragg_f_mo(ipr,igr,ipl) = - img * Bragg(igr,ipl) * f_mo(ipr,ipl)
          if( Bulk_step .and. Bulk_roughness > eps10 ) then
            Bragg_rgh_Bulk_f_ms(ipr,igr,ipl) = - img * Bragg_rgh_Bulk(igr,ipl) * f_ms(ipr,ipl)
            Bragg_rgh_Bulk_f_mo(ipr,igr,ipl) = - img * Bragg_rgh_Bulk(igr,ipl) * f_mo(ipr,ipl)
          endif
        end do

      end do

    endif

    do ipl = 1,npldafs
      Sum_Bragg_f0(ipr,ipl) = Sum_Bragg(ipr,ipl) * f0(ipr,ipl)
      Sum_Bragg_fa(ipr,ipl) = Sum_Bragg(ipr,ipl) * cmplx(fp(ipr), fpp(ipr), db)

      if( Bulk_step ) then
        if( ipr == iprabs_nonexc ) then
          do i = 1,n_bulk_z_abs
            Sum_Bragg_bulk_abs_f0(ipl,i) = Sum_Bragg_abs(ipl,i) * f0(ipr,ipl)
          end do 
        else
          do igr = 1,ngreq(ipr)
            do i = 1,n_bulk_z
              do j = 1,n_bulk_zc(i)
                if( igr_bulk_z(j,i) == igreq(ipr,igr) ) then
                  Sum_Bragg_nonabs_f0(ipl,i) = Sum_Bragg_nonabs_f0(ipl,i) + Bragg(igr,ipl) * f0(ipr,ipl)
                  Sum_Bragg_nonabs_fa(ipl,i) = Sum_Bragg_nonabs_fa(ipl,i) + Bragg(igr,ipl) * cmplx(fp(ipr), fpp(ipr), db)
                endif
              end do
            end do
          end do
        endif
        if( Bulk_roughness > eps10 ) then
          Sum_Bragg_rgh_Bulk_f0(ipr,ipl) = Sum_Bragg_rgh_Bulk(ipr,ipl) * f0(ipr,ipl)
          Sum_Bragg_rgh_Bulk_fa(ipr,ipl) = Sum_Bragg_rgh_Bulk(ipr,ipl) * cmplx(fp(ipr), fpp(ipr), db)
        endif
      endif
      
      if( Dafs_bio ) then
        Sum_Bragg_f0_cj(ipr,ipl) = conjg( Sum_Bragg(ipr,ipl) ) * f0(ipr,ipl)
        Sum_Bragg_fa_cj(ipr,ipl) = conjg( Sum_Bragg(ipr,ipl) ) * cmplx(fp(ipr), fpp(ipr), db)
      endif
    end do

    deallocate( Bragg )
    if( Bulk_step .and. Bulk_roughness > eps10 ) deallocate( Bragg_rgh_Bulk )

  end do  ! End of loop over prototypical atoms ipr

  if( icheck > 1 .and. Magnetic ) then
    write(3,'(/A)') ' ipr ipl      f_ms         f_mo      Bragg_f_ms(ipr,ipl,1..ngr)'
    do ipr = n1_proto,n2_proto
      if( ( Bulk_step .and. ipr <= n_atom_proto - n_atom_proto_bulk ) .or. ( .not. Bulk_step .and. ipr > n_atom_proto_uc ) ) cycle
      do ipl = 1,npldafs
        if( abs(f_ms(ipr,ipl)) < eps10 ) cycle
        write(3,'(2i4,1p,50e13.5)') ipr, ipl, f_ms(ipr,ipl), f_mo(ipr,ipl), Bragg_f_ms(ipr,1:ngreq(ipr),ipl)
      end do
    end do
  endif

  if( Film ) then
    allocate( phd_f_bulk(npldafs) )
    phd_f_bulk(:) = ( 0._db, 0._db )
    if( Bulk ) allocate( Phase_bulk(npldafs) )
  endif

  if( Bulk ) then

    call Cal_Phase_Bulk(angxyz,angxyz_bulk,axyz,axyz_bulk,c_cos_z,c_cos_z_b,Delta_bulk,Film_shift,hkl_dafs,hkl_film, &
                          Mult_bulk,npldafs,Phase_bulk)

    if( .not. Bulk_step .and. ( .not. Abs_in_bulk .or. Bulk_roughness > eps10 ) ) then
 ! Calculate the bulk scattering only when there is no absorbing atom in it
 ! In the other case, it only calculates the non resonant part coming from the bulk roughness
      call Bulk_scat(Abs_in_bulk,angxyz_bulk,axyz_bulk,Bulk_roughness,c_cos_z_b,Delta_bulk,delta_z_top_bulk,Energy_photon, &
        hkl_dafs,hkl_film,icheck,Length_abs,Mat_bulk_i,Mult_bulk,n_atom_bulk,ngroup_taux,ngroup_temp,npldafs,phd_f_bulk, &
        posn_bulk,Q_mod,Taux_oc,Temp_coef,Temp_B_iso,Truncation,Z_abs,Z_bulk)

      phd_f_bulk(:) = Phase_bulk(:) * phd_f_bulk(:)
    endif

  endif

  if( Film ) then
    allocate( phd_f_cap(npldafs) )
    phd_f_cap(:) = ( 0._db, 0._db )
  endif

  if( Cap_layer .and. .not. Bulk_step ) &
    call Cap_scat(angxyz,angxyz_bulk,angxyz_cap,axyz,axyz_bulk,axyz_cap,Cap_B_iso,Cap_roughness, &
              Cap_thickness,c_cos_z,c_cos_z_b,Delta_cap,Delta_roughness_film,delta_z_bottom_cap, &
              delta_z_top_cap,Energy_photon,Film_roughness,fpp_cap_tot,hkl_dafs,hkl_film, &
              Mult_bulk,Mult_film,n_atom_cap,npldafs,phd_f_cap,posn_cap,Q_mod,Taux_cap, &
              Z_cap)

  do ipl = npldafs,1,-1

    Sum_T_Bragg_f = sum( Sum_Bragg_f0(n1_proto:n2_proto,ipl) ) + sum( Sum_Bragg_fa(n1_proto:n2_proto,ipl) )
    if( Bulk_step .and. Bulk_roughness > eps10 ) &
      Sum_T_Bragg_rgh_Bulk_f = sum( Sum_Bragg_rgh_Bulk_f0(n1_proto:n2_proto,ipl) ) &
                             + sum( Sum_Bragg_rgh_Bulk_fa(n1_proto:n2_proto,ipl) )

    if( .not. Bulk_step ) then
      if( Bulk ) Sum_T_Bragg_f = Sum_T_Bragg_f + phd_f_bulk(ipl)

 ! Absorption through the cap. Well, it is extremely small...
      if( Cap_layer ) then
        x = Q_mod(ipl) / konde      ! = 2 * sin(theta)
        abs_cap = exp( - fpp_cap_tot * x )
        if( ipl == 1 .and. icheck == 1 )  write(3,'(a29,1p,e13.5)') ' Absorption through the cap =', 1._db - Abs_cap
        Sum_T_Bragg_f = abs_cap * Sum_T_Bragg_f
        Bragg_abs(1:natomsym,ipl) = abs_cap * Bragg_abs(1:natomsym,ipl)
        Sum_T_Bragg_f = Sum_T_Bragg_f + phd_f_cap(ipl)
      endif

      if( Dafs_bio ) ph_cjg = sum( Sum_Bragg_f0_cj(1:n_atom_proto_uc,ipl) ) + sum( Sum_Bragg_fa_cj(1:n_atom_proto_uc,ipl) )
    endif

! Multplication by epsilon_e * espilon_s
    do ip = nphi_dafs(ipl),1,-1
      pe(:) = Poldafse(:,ipl,ip)
      ps(:) = Poldafss(:,ipl,ip)
      cfac = sum( conjg( ps(:) ) * pe(:) )
      if( Dafs_bio .and. ip > 2 ) then
        phdf0t(ipl,ip) = cfac * ph_cjg
        phdt(ipl,ip,:) = cfac * conjg( Sum_Bragg_abs(ipl,:) )
      else
        phdf0t(ipl,ip) = cfac * Sum_T_Bragg_f
        phdt(ipl,ip,:) = cfac * Sum_Bragg_abs(ipl,:)
        if( Bulk_step .and. Bulk_roughness > eps10 ) then
          phdf0t_rgh_Bulk(ipl,ip) = cfac * Sum_T_Bragg_rgh_Bulk_f
          phdt_rgh_Bulk(ipl,ip) = cfac * Sum_Bragg_rgh_Bulk(iprabs_nonexc,ipl)
        endif
      endif
    end do

    Sum_Bragg_f0(n1_proto:n2_proto,ipl) = cfac * Sum_Bragg_f0(n1_proto:n2_proto,ipl)
    Sum_Bragg_fa(n1_proto:n2_proto,ipl) = cfac * Sum_Bragg_fa(n1_proto:n2_proto,ipl)
    if( Bulk_step ) then
      Sum_Bragg_bulk_abs_f0(ipl,:) = cfac * Sum_Bragg_bulk_abs_f0(ipl,:)
      Sum_Bragg_nonabs_f0(ipl,:) = cfac * Sum_Bragg_nonabs_f0(ipl,:)
      Sum_Bragg_nonabs_fa(ipl,:) = cfac * Sum_Bragg_nonabs_fa(ipl,:)
      if( Bulk_roughness > eps10 ) then
        Sum_Bragg_rgh_Bulk_f0(n1_proto:n2_proto,ipl) = cfac * Sum_Bragg_rgh_Bulk_f0(n1_proto:n2_proto,ipl)
        Sum_Bragg_rgh_Bulk_fa(n1_proto:n2_proto,ipl) = cfac * Sum_Bragg_rgh_Bulk_fa(n1_proto:n2_proto,ipl)
      endif
    endif

  end do ! end loop on ipl

! See M. Blume and Doon Gibbs, PRB, 37, 1779 (1988)
! Formula modified for vect A'' because error in the paper.
  Sum_T_Bragg_fmag(:,:) = (0._db,0._db)
  if( Magnetic ) then
    Axe_atom_g = matmul( orthmatt, Axe_atom_gr )
    do ipl = 1,npldafs
      do ip = nphi_dafs(ipl),1,-1
        pe(:) = Poldafse(:,ipl,ip)
        ps(:) = Poldafss(:,ipl,ip)
        we(:) = Vecdafse(:,ipl,ip)
        ws(:) = Vecdafss(:,ipl,ip)
! Approximation for orbital moment taken has parralel.
! See also appendix by G. T. Trammell, PR 92, 1387 (1953)
        call get_vec_b(pe,ps,vec_b,we,ws)
        call get_vec_a(pe,ps,vec_a,we,ws)
        do ipr = n1_proto,n2_proto

          do igr = 1,ngreq(ipr)
            jgr = igreq(ipr,igr)
            v(:) = Axe_atom_g(:,jgr)
            v_vec_b(ipl) = dot_product(v,vec_b)
            v_vec_a(ipl) = dot_product(v,vec_a)
            if( icheck > 1 ) then
              write(3,182) ipl
              do i = 1,3
                write(3,185) v(i), vec_b(i), vec_a(i)
              end do
              write(3,187) v_vec_b(ipl), v_vec_a(ipl)
            endif

            if( ip == 1 ) &
              Sum_T_Bragg_fmag(ipr,ipl) = Sum_T_Bragg_fmag(ipr,ipl) + v_vec_a(ipl) * Bragg_f_mo(ipr,igr,ipl) &
                                                                    + v_vec_b(ipl) * Bragg_f_ms(ipr,igr,ipl)

            phdf0t(ipl,ip) = phdf0t(ipl,ip) + v_vec_a(ipl) * Bragg_f_mo(ipr,igr,ipl) + v_vec_b(ipl) * Bragg_f_ms(ipr,igr,ipl)
            if( Bulk_step .and. Bulk_roughness > eps10 ) &
              phdf0t_rgh_Bulk(ipl,ip) = phdf0t_rgh_Bulk(ipl,ip) + v_vec_a(ipl) * Bragg_rgh_Bulk_f_mo(ipr,igr,ipl) &
                                                                + v_vec_b(ipl) * Bragg_rgh_Bulk_f_ms(ipr,igr,ipl)

            if( Bulk_step .and. ip == 1 ) then

              if( ipr == iprabs_nonexc ) then
                do i = 1,n_bulk_z_abs
                  do j = 1,n_bulk_zc_abs(i)
                    if( igr == igr_bulk_z_abs(j,i) ) Sum_Bragg_bulk_abs_f0(ipl,:) = Sum_Bragg_bulk_abs_f0(ipl,:) &
                                                                           + v_vec_a(ipl) * Bragg_f_mo(ipr,igr,ipl) &
                                                                           + v_vec_b(ipl) * Bragg_f_ms(ipr,igr,ipl)
                  end do
                end do

              else
                do i = 1,n_bulk_z_abs
                  do j = 1,n_bulk_zc_abs(i)
                    if( jgr == igr_bulk_z_abs(j,i) ) Sum_Bragg_nonabs_f0(ipl,:) = Sum_Bragg_nonabs_f0(ipl,:) &
                                                                           + v_vec_a(ipl) * Bragg_f_mo(ipr,igr,ipl) &
                                                                           + v_vec_b(ipl) * Bragg_f_ms(ipr,igr,ipl)
                  end do
                end do
              endif

            endif
 
          end do
        end do
      end do

    end do

  endif

  if( Bulk_step ) then
    do ipl = 1,npldafs

      do ip = 1,nphi_dafs(ipl)
        phdf0t(ipl,ip) = Phase_bulk(ipl) * phdf0t(ipl,ip)
        phdt(ipl,ip,:) = Phase_bulk(ipl) * phdt(ipl,ip,:)
      end do
      Bragg_abs(1:natomsym,ipl) = Phase_bulk(ipl) * Bragg_abs(1:natomsym,ipl)
      Sum_Bragg_bulk_abs_f0(ipl,:) = Phase_bulk(ipl) * Sum_Bragg_bulk_abs_f0(ipl,:)
      Sum_Bragg_nonabs_f0(ipl,:) = Phase_bulk(ipl) * Sum_Bragg_nonabs_f0(ipl,:)
      Sum_Bragg_nonabs_fa(ipl,:) = Phase_bulk(ipl) * Sum_Bragg_nonabs_fa(ipl,:)

      if( Bulk_roughness > eps10 ) then
        do ip = 1,nphi_dafs(ipl)
          phdf0t_rgh_Bulk(ipl,ip) = Phase_bulk(ipl) * phdf0t_rgh_Bulk(ipl,ip)
        end do
        Bragg_rgh_Bulk_abs(1:natomsym,ipl) = Phase_bulk(ipl) * Bragg_rgh_Bulk_abs(1:natomsym,ipl)
      endif

    end do
    
    Sum_Bragg_nonabs_f(:,:) = Sum_Bragg_nonabs_f0(:,:) + Sum_Bragg_nonabs_fa(:,:)
  endif

!----- Writing --------------------------------------------------------

  if( icheck > 1 ) then
    do ipl = 1,npldafs
      if( nphi_dafs(ipl) == 1 ) cycle
      dpdeg = 360._db / nphi_dafs(ipl)
      if( Dafs_bio ) then
        write(3,200)
      elseif( Film .or. Bulk_step ) then
        write(3,205) hkl_dafs(:,ipl), isigpi(:,ipl)
      else
        write(3,210) nint( hkl_dafs(:,ipl) ), isigpi(:,ipl)
      endif
      do ip = 1,nphi_dafs(ipl)
        pe(:) = Poldafse(:,ipl,ip)
        ps(:) = Poldafss(:,ipl,ip)
        we(:) = Vecdafse(:,ipl,ip)
        ws(:) = Vecdafss(:,ipl,ip)

        pe = matmul( orthmati, pe )
        ps = matmul( orthmati, ps )
        we = matmul( orthmati, we )
        ws = matmul( orthmati, ws )
        if( Dafs_bio ) then
          if( ip == 1 .or. ip == 3 ) then
            mot4 = ' sig'
          else
            mot4 = '  pi'
          endif
          if( ip == 1 .or. ip == 2 ) then
            write(3,215) nint(hkl_dafs(:,ipl)), mot4, real(pe(:)), real(ps(:)), we(:), ws(:)
          else
            write(3,215) - nint(hkl_dafs(:,ipl)), mot4, real(pe(:)), real(ps(:)), we(:), ws(:)
          endif
        else
          write(3,220)  ( ip - 1 ) * dpdeg, real(pe(:)), aimag(pe(:)), real(ps(:)), aimag(ps(:)), we(:), ws(:)
        endif
      end do
    end do
  endif

  if( icheck > 0 ) then

    if( Film .or. Bulk_step ) then
      write(3,'(/A)') '  (   h,      k,      l)  exp(i*Q.R_ia) (ia = 1,natomsym)'
    else
      write(3,'(/A)') '  (h,k,l)  exp(i*Q.R_ia) (ia = 1,natomsym)'
    endif

    do ipl = 1,npldafs
      if( Film .or. Bulk_step ) then
        if( mod(ipl,10) /= 1 .and. npldafs > 20 .and. icheck == 1 ) cycle
        write(3,235) hkl_dafs(:,ipl), Bragg_abs(1:natomsym,ipl)
      else
        write(3,240) nint(hkl_dafs(:,ipl)), Bragg_abs(1:natomsym,ipl)
      endif
    end do

    if( Bulk_step .and. Bulk_roughness > eps10 ) then
      write(3,'(/A)') '  (   h,      k,      l)  exp(i*Q.R_ia) (ia = 1,natomsym), roughness'
      do ipl = 1,npldafs
        write(3,235) hkl_dafs(:,ipl), Bragg_rgh_Bulk_abs(1:natomsym,ipl) 
      end do
    endif

  endif

  if( icheck > 0 .and. .not. Bulk_step ) then

    if( Magnetic ) then
      write(3,245) f_no_res(1)
      ipl = 1  ! The ratio does not depend on the reflexion.
      do ipr = n1_proto,n2_proto
        if( abs( f_ms(ipr,ipl) ) > eps10 ) then
          rap_lsur2s = 0.5_db * f_mo(ipr,ipl) / f_ms(ipr,ipl)
        else
          rap_lsur2s = 0._db
        endif
        write(3,246) ipr, rap_lsur2s
      end do
    endif

    write(3,250)
    do ipr = n1_proto,n2_proto
      if( nspin == 1 ) then
        if( Film .or. Bulk_step ) then
          write(3,260) ipr, numat( itypepr(ipr) )
        else
          write(3,265) ipr, numat( itypepr(ipr) )
        endif
        do ipl = 1,npldafs
          if( Film .or. Bulk_step ) then
            if( mod(ipl,10) /= 1 .and. npldafs > 20 .and. icheck == 1 ) cycle
            write(3,270) hkl_dafs(:,ipl), f0(ipr,ipl), fp(ipr), fpp(ipr), Sum_Bragg(ipr,ipl), Sum_Bragg_f0(ipr,ipl), &
                         Sum_Bragg_fa(ipr,ipl)
          else
            write(3,275) nint(hkl_dafs(:,ipl)), f0(ipr,ipl), fp(ipr), fpp(ipr), Sum_Bragg(ipr,ipl), Sum_Bragg_f0(ipr,ipl), &
                         Sum_Bragg_fa(ipr,ipl)
          endif
        end do
      else
        if( Film .or. Bulk_step ) then
          write(3,280) ipr, numat( itypepr(ipr) )
        else
          write(3,285) ipr, numat( itypepr(ipr) )
        endif
        do ipl = 1,npldafs
          if( Film .or. Bulk_step ) then
            if( mod(ipl,10) /= 1 .and. npldafs > 20 .and. icheck == 1 ) cycle
            write(3,290) hkl_dafs(:,ipl), f0(ipr,ipl), fp(ipr), fpp(ipr), f_ms(ipr,ipl), f_mo(ipr,ipl), Sum_Bragg(ipr,ipl), &
              Sum_Bragg_f0(ipr,ipl), Sum_Bragg_fa(ipr,ipl), Sum_T_Bragg_fmag(ipr,ipl), v_vec_b(ipl), v_vec_a(ipl)
          else
            write(3,295) nint(hkl_dafs(:,ipl)), f0(ipr,ipl), fp(ipr), fpp(ipr), f_ms(ipr,ipl), f_mo(ipr,ipl), Sum_Bragg(ipr,ipl), &
              Sum_Bragg_f0(ipr,ipl), Sum_Bragg_fa(ipr,ipl), Sum_T_Bragg_fmag(ipr,ipl), v_vec_b(ipl), v_vec_a(ipl)
          endif
        end do
      endif
    end do

    if( Film ) then
      mot64 =  '  (   h,      k,       l)     Total Structure factor          I'
      if( Bulk .and. .not. Abs_in_bulk ) then
        if( Cap_layer ) then
          write(3,'(/a64,18x,a4,3x,3(18x,a10))') mot64, 'Bulk', 'Truncation', 'Phase_Bulk', 'Cover cap '
        else
          write(3,'(/a64,18x,a4,3x,2(18x,a10))') mot64, 'Bulk', 'Truncation', 'Phase_Bulk'
        endif
      elseif( Cap_layer ) then
        write(3,'(/a64,18x,a10)') mot64, 'Cover cap '
      else
        write(3,'(/a64,18x,a10)') mot64
      endif
    else
      write(3,'(/A)') '  (h,k,l)   Total Structure factor       I'
    endif

    do ipl = 1,npldafs
      if( Film ) then
        if(  mod(ipl,10) /= 1 .and. npldafs > 20 .and. icheck == 1 ) cycle
        if( Bulk .and. .not. Abs_in_bulk ) then
        arg = deux_pi*hkl_dafs(3,ipl)*Film_shift(3)/c_cos_z
        Sum_T_Bragg_f = exp(arg*img)
          if( Cap_layer ) then
            write(3,300) hkl_dafs(:,ipl), phdf0t(ipl,1), abs( phdf0t(ipl,1) )**2, phd_f_bulk(ipl), Truncation(ipl), &
                         Phase_bulk(ipl), phd_f_cap(ipl)
          else
            write(3,300) hkl_dafs(:,ipl), phdf0t(ipl,1), abs( phdf0t(ipl,1) )**2, phd_f_bulk(ipl), Truncation(ipl), Phase_bulk(ipl)
          endif
        elseif( Cap_layer ) then
          write(3,300) hkl_dafs(:,ipl), phdf0t(ipl,1), abs( phdf0t(ipl,1) )**2, phd_f_cap(ipl)
        else
          write(3,300) hkl_dafs(:,ipl), phdf0t(ipl,1), abs( phdf0t(ipl,1) )**2
        endif
      else
        write(3,310) nint( hkl_dafs(:,ipl) ), phdf0t(ipl,1), abs( phdf0t(ipl,1) )**2
      endif
    end do

  endif

  if( icheck > 1 ) then
    iwrite = 1
    do ipl = 1,npldafs
      if( nphi_dafs(ipl) == 1 ) cycle
      if( iwrite == 1 ) then
        write(3,375)
        iwrite = 0
      endif
      write(3,376) ipl
      do ip = 1,nphi_dafs(ipl)
        write(3,*)
        do i = 1,3
         write(3,380) Poldafse(i,ipl,ip), Vecdafse(i,ipl,ip), Poldafss(i,ipl,ip), Vecdafss(i,ipl,ip)
        end do
      end do
    end do
  endif

  if( Film ) then
    deallocate( phd_f_bulk, phd_f_cap )
    if( Bulk ) deallocate( Phase_bulk )
  endif

  if( Magnetic ) then
    deallocate( Bragg_f_mo )
    deallocate( Bragg_f_ms )
  endif

  if( Bulk_step .and. Bulk_roughness > eps10 ) then
    deallocate( Sum_Bragg_rgh_Bulk_f0 )
    deallocate( Sum_Bragg_rgh_Bulk_fa )
    if( Magnetic ) then
      deallocate( Bragg_rgh_Bulk_f_mo )
      deallocate( Bragg_rgh_Bulk_f_ms )
    endif
  endif
  
  return
  110 format(/' ---- Prepdafs -----',100('-'))
  165 format(1x,3i3,i5,10x,f9.5)
  182 format(/' ipl =',i3,//'  Spin_axis          vec_b',17x,'vec_a')
  185 format(f10.5,2(2x,2f10.5))
  187 format(/' V.Vec_b =',2f10.5,',  V.Vec_a =',2f10.5)
  200 format(/'  (h, k, l)   pol',11x,'Poldafse',20x, 'Poldafss',20x,'Vecdafse',20x,'Vecdafss')
  205 format(/' (h,k,l) = (',3f8.3,'),  isigpi =',2i3,/' angle', &
        8x,'real(Poldafse)',13x,'imag(Poldafse)',11x,'real(Poldafss)',12x, 'imag(Poldafse)',15x,'Vecdafse',15x,'Vecdafss')
  210 format(/' (h,k,l) = (',3i3,'),  isigpi =',2i3,/' angle', &
        8x,'real(Poldafse)',13x,'imag(Poldafse)',11x,'real(Poldafss)',12x, 'imag(Poldafse)',15x,'Vecdafse',15x,'Vecdafss')
  215 format(3i4,1x,a4, 6(1x,3f9.5))
  220 format(f6.1, 6(1x,3f9.5))
  235 format(2f8.3,f8.4,1p,48(1x,2e13.5))
  240 format(3i3,1p,48(1x,2e13.5))
  245 format(/' Attenuation factor for non resonant magnetic structure', ' factor :'/,'      General : fma =',f7.3,/ &
        '      Orbital : Site    L/2S     (f_orb = (L/S)*f_spin)')
  246 format(10x,i9,f9.3)
  250 format(/' Value of the structure factors per atom site :')
  260 format(/' Site =',i3,', Z =',i3,//, &
        '  (     h,     k,      l)',6x,'f0',14x,'fp',10x,'fpp',14x, &
        'Sum_Bragg',12x,'Sum_Bragg*f0*(eps_e.eps_s)',5x,'Sum_Bragg*(fp+i*fpp)*(eps_e.eps_s)')
  265 format(/' Site =',i3,', Z =',i3,//, &
        '  (h,k,l)',6x,'f0',14x,'fp',10x,'fpp',14x, &
        'Sum_Bragg',12x,'Sum_Bragg*f0*(eps_e.eps_s)',5x,'Sum_Bragg*(fp+i*fpp)*(eps_e.eps_s)')
  270 format(2f8.3,f9.4,1p,e13.5,2x,2e13.5,5(2x,2e13.5))
  275 format(3i3,1p,e13.5,2x,2e13.5,5(2x,2e13.5))
  280 format(/' Site =',i3,', Z =',i3,//, &
        '  (     h,     k,      l)',6x,'f0',14x,'fp',10x,'fpp',11x,'f_spin',9x,'f_orb',14x, &
        'Sum_Bragg',12x,'Sum_Bragg*f0*(eps_e.eps_s)',5x,'Sum_Bragg*(fp+i*fpp)*(eps_e.eps_s)', &
        ' i*Bragg_m*fma*(f_spin*vs+f_orb*vo) vs = spin.vec_b',12x, 'vo = spin.vec_a')
  285 format(/' Site =',i3,', Z =',i3,//, &
        '  (h,k,l)',6x,'f0',14x,'fp',10x,'fpp',11x,'f_spin',9x,'f_orb',14x, &
        'Sum_Bragg',12x,'Sum_Bragg*f0*(eps_e.eps_s)',5x,'Sum_Bragg*(fp+i*fpp)*(eps_e.eps_s)', &
        ' i*Bragg_m*fma*(f_spin*vs+f_orb*vo) vs = spin.vec_b',12x, 'vo = spin.vec_a')
  290 format(2f8.3,f9.4,1p,e13.5,2x,2e13.5,2(2x,e13.5),6(2x,2e13.5))
  295 format(3i3,1p,e13.5,2x,2e13.5,2(2x,e13.5),6(2x,2e13.5))
  300 format(2f8.3,f9.4,1p,2x,2e13.5,2x,e13.5,8(2x,2e13.5))
  310 format(3i3,1p,2x,2e13.5,2x,e13.5,5(2x,2e13.5))
  375 format(/'  Polarization and wave vectors in the internal basis R1 ( orthogonal basis, z along c crystal )',// &
              '         pol            vec           pls',12x,'vos')
  376 format(/' ipl = ',i5)
  380 format(2(2f9.5,1x,f9.5,1x))

end

!*********************************************************************

! Calculation of polarization using the operation mode

subroutine SRXD(Angle_mode,angpoldafs,angxyz,angxyz_bulk,axyz,axyz_bulk,Bulk,Bulk_step,Energy_photon,Film,hkl_dafs,hkl_film, &
                hkl_ref,icheck,isigpi,Length_abs,Mat_UB,nphim,npl_2d,npldafs,npldafs_2d,Operation_mode,Orthmati,phi_0,Poldafse, &
                Poldafss,Q_mod,Vecdafse,Vecdafss)

  use declarations
  implicit none

  integer:: i, icheck, ipl, ipr, j, jpl, kpl, lpl, mpl, ni, npl, npldafs, npldafs_2d, nphim

  integer, dimension(3):: hkl_ref
  integer, dimension(5):: Operation_m
  integer, dimension(npldafs_2d):: npl_2d
  integer, dimension(2,npldafs):: isigpi
  integer, dimension(5,npldafs_2d):: Operation_mode

  complex(kind=db), dimension(3,npldafs,nphim):: Poldafse, Poldafss

  logical:: Alfa_Beta, Alfa_fixed, Below, Beta_fixed, Bulk, Bulk_step, Circular, Chi_fixed, Column_reference, Column_detector, &
            Delta_fixed, Detector_known, Eta_fixed, Eta_half_delta, Film, hkl_film, Mat_unit, Mu_fixed, Mu_half_nu, Naz_fixed, &
            Nu_fixed, Phi_fixed, Pol_pi, Psi_fixed, Q_par_n, Qaz_fixed, Ref_spec, Specular, UB_in

  real(kind=db):: a, b, alfa, Angfac, Angle, beta, c, c_cos_z, c_x, c_y, c_z, chi, cos_c, cos_d, cos_e, cos_k, &
        cos_m, cos_n, cos_p, cos_t, cos_tau, cos_tau_ref, cos_z, da, delta, Energy_photon, eta, Konde, ksi,  mu, Naz, nu, phi, &
        phi_0, psi, Qaz, sin_c, sin_e, sin_k, sin_m, sin_p, sin_t, sin_tau, sin_tau_ref, tau, tau_ref, Theta_Bragg, x
  real(kind=db), dimension(3):: ang, Angle_m, angxyz, angxyz_bulk, abc, axyz, axyz_bulk, cosdir, H_phi, H_phi_n, k_i, k_s, N_L, &
                                N_phi, Norm_phi, p, p_i, p_i_p, p_i_s, p_s_p, p_s_s, Q, Q_L, sindir, Vec, W
  real(kind=db), dimension(3,3):: Mat, Mat_B, Mat_K, Mat_K_i, Mat_NL, Mat_NN, Mat_Nphi, Mat_p, Mat_U, Mat_UB, Mat_UB_loc, &
                                  Mat_UB_t, Orthmati
  real(kind=db), dimension(npldafs):: Length_abs, Q_mod
  real(kind=db), dimension(3,npldafs):: angpoldafs, hkl_dafs
  real(kind=db), dimension(3,npldafs_2d):: Angle_mode
  real(kind=db), dimension(3,npldafs,nphim):: Vecdafse, Vecdafss

! lambda = 2 * pi / k = 2 * d * sintheta
! In S.I. vecond = k = E*alfa*4*pi*epsilon0 / (e*e)
! In ua and rydb : k = 0.5 * alfa * E
  konde = 0.5_db * alfa_sf * Energy_photon

  if( Bulk .and. .not. hkl_film ) then
 ! When bulk_step angxyz_bulk contains film data ! ( See fdm )
    abc(:) = axyz_bulk(:)
    ang(:) = angxyz_bulk(:)
  else
    abc(:) = axyz(:)
    ang(:) = angxyz(:)
  endif

  if( Film ) then
    if( Bulk_step ) then
      cos_z = sqrt( sin( angxyz(2) )**2 - ( ( cos( angxyz(1) ) - cos( angxyz(3) ) * cos( angxyz(2) ) ) / sin( angxyz(3) ) )**2 )
      c_cos_z = axyz(3) * cos_z
    else
      cos_z = sqrt( sin( angxyz_bulk(2) )**2 - ( ( cos( angxyz_bulk(1) ) - cos( angxyz_bulk(3) ) * cos( angxyz_bulk(2) ) ) &
                                           / sin( angxyz_bulk(3) ) )**2 )
      c_cos_z = axyz_bulk(3) * cos_z
    endif
! The factor is to convert Length_abs en micrometer
    c_cos_z = 0.0001_db * bohr * c_cos_z
  endif

  cosdir(:) = cos( ang(:) )
  sindir(:) = sin( ang(:) )

! Unit cell volume is a.b.c.A
! Angfac = 2 pi / A
  Angfac = deux_pi / sqrt( 1 - sum( cosdir(:)**2 ) + 2 * cosdir(1) * cosdir(2) * cosdir(3) )

! Transformation matrix of Q form (h,k,l) to cartesian coordinates (e_x, e_y, e_z)
! with e_y along b, e_z along c* and e_x = e_y x e_e_z.
  c_x = ( cosdir(2) - cosdir(1) * cosdir(3) ) / sindir(3)
  c_y = cosdir(1)
  c_z = sqrt( 1 - cosdir(1)**2 - c_x**2 )

  Mat_B(1,1) = c_z / abc(1);    Mat_B(1,2) = - cosdir(3) * c_z / abc(2);                      Mat_B(1,3) = 0._db
  Mat_B(2,1) = 0._db;           Mat_B(2,2) = sindir(3) * c_z / abc(2);                        Mat_B(2,3) = 0._db
  Mat_B(3,1) = - c_x / abc(1);  Mat_B(3,2) = ( cosdir(3) * c_x - sindir(3) * c_y ) / abc(2);  Mat_B(3,3) = sindir(3) / abc(3)

! Orientation matrix

  call Rotation_Matrix(-phi_0,3,Mat_U)
  Mat_UB_t = Matmul( Mat_U, Mat_B )

  if( sum( abs( Mat_UB(:,:) ) ) < eps10 ) then
    UB_in = .false.
    Mat_UB_loc(:,:) = Mat_UB_t(:,:)
  else
    UB_in = .true.
    Mat_UB_loc(:,:) = Mat_UB(:,:) / Angfac
  endif

  call invermat( Mat_UB_loc, Mat )
  Mat_K = matmul( Orthmati, Mat )
  Mat(:,:) = Mat_K(:,:)
  do i = 1,3
    Mat(i,i) = Mat(i,i) - 1
  end do
  x = sum( abs( Mat(:,:) ) )
  Mat_unit = x < eps10

  call invermat(Mat_K,Mat_K_i)

  Mat_UB_loc(:,:) = Angfac * Mat_UB_loc(:,:)

! Normal to the surface
  Norm_phi(1:2) = 0._db; Norm_phi(3) = 1._db
  Norm_phi = Matmul( Mat_UB_loc, Norm_phi )
  x = sqrt( sum( Norm_phi(:)**2 ) )
  Norm_phi(:) = Norm_phi(:) / x

! Reference vector which by default is the normal to the surface
  N_phi = Matmul( Mat_UB_loc, hkl_ref )
  x = sqrt( sum( N_phi(:)**2 ) )
  N_phi(:) = N_phi(:) / x

  Ref_spec = abs( 1 - dot_product( N_phi, Norm_phi ) ) < eps10

  if( icheck > 0 ) then
    if( .not. UB_in ) then
      write(3,'(/19x,a5,33x,a5,33x,a6)') 'Mat B', 'Mat U', 'Mat UB'
      do i = 1,3
        write(3,110) Mat_B(i,:), Mat_U(i,:), Mat_UB_loc(i,:)
      end do
    elseif( Mat_unit ) then
      write(3,'(/18x,a6,30x,a11)') 'Mat UB', 'Mat UB theo'
      do i = 1,3
        write(3,110) Mat_UB_loc(i,:),  Angfac * Mat_UB_t(i,:)
      end do
    else
      write(3,'(/18x,a6,30x,a11,29x,a7)') 'Mat UB', 'Mat UB theo', 'Mat Rot'
      do i = 1,3
        write(3,110) Mat_UB_loc(i,:), Angfac * Mat_UB_t(i,:), Mat_K(i,:)
      end do
    endif
    write(3,'(/2x,a21,3i3)') 'Reference axis, hkl =', hkl_ref(:)
  endif

  ipl = 0
  do jpl = 1,npldafs_2d

    Circular = Operation_mode(1,jpl) < 0 .or. Operation_mode(2,jpl) < 0
      
    Operation_m(:) = abs( Operation_mode(:,jpl) ) ! because negative values means circular polarization
    Angle_m(:) = Angle_mode(:,jpl)

    call Operation_mode_eval(alfa,Alfa_beta,Alfa_fixed,Angle_m,beta,Beta_fixed,chi,Chi_fixed,Column_detector, &
                Column_reference,delta,Delta_fixed,eta,Detector_known,Eta_fixed,Eta_half_delta,icheck,mu,Mu_fixed, &
                Mu_half_nu,Naz,Naz_fixed,nu,Nu_fixed,Operation_m,phi,Phi_fixed,psi,Psi_fixed,Qaz,Qaz_fixed)

    do kpl = 1,npl_2d(jpl)
      ipl = ipl + 1
      lpl = ipl + npl_2d(jpl)
      if( Circular ) then
        mpl = lpl + npl_2d(jpl)
        npl = mpl + npl_2d(jpl)
      else
        npl = lpl
      endif
      Q(:) = hkl_dafs(:,ipl)

      Specular = abs( Q(1) ) < eps10 .and. abs( Q(2) ) < eps10

      if( Specular .and. abs( Q(3) ) < eps10 ) Q(3) =  1.0001_db * eps10
      H_phi = Matmul( Mat_UB_loc, Q )

! Q modulus = mod(k_s - k_i) = 2 * k * sin( Theta_bragg )
      Q_mod(ipl) = sqrt( sum( H_phi(:)**2 ) )
      Q_mod(lpl) = Q_mod(ipl)
      if( Circular ) then
        Q_mod(mpl) = Q_mod(ipl)
        Q_mod(npl) = Q_mod(ipl)
      endif  
      H_phi_n(:) = H_phi(:) / Q_mod(ipl)

! Theta_Bragg determination
      x = Q_mod(ipl) / ( 2 * konde )
      if( x > 1._db ) then
        call write_error
        do ipr = 6,9,3
          if( ipr == 3 .and. icheck == 0 ) cycle
          write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
          write(ipr,'(/A//)') ' The reflexion does not exist at this energy !'
        end do
        stop
      endif

      Theta_Bragg = asin( x )
      cos_t = cos( Theta_Bragg )
      sin_t = sin( Theta_Bragg )
! Angle between reference vector (by default normal to the surface) and Q
      cos_tau_ref = dot_product( H_phi_n, N_phi )
      cos_tau = dot_product( H_phi_n, Norm_phi )
      if( abs( cos_tau - 1 ) < eps10 ) then
        tau = 0._db
      else
        tau = acos( cos_tau )
      endif
      if( abs( cos_tau_ref - 1 ) < eps10 ) then
        tau_ref = 0._db
      else
        tau_ref = acos( cos_tau_ref )
      endif
      sin_tau = sin( tau )
      sin_tau_ref = sin( tau_ref )
      Q_par_n = abs( sin_tau_ref ) < eps10
      if( Q_par_n .and. .not. Specular ) then
        call write_error
        do ipr = 6,9,3
          if( ipr == 3 .and. icheck == 0 ) cycle
          write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
          write(ipr,'(/a15,3i3)') ' Reference axis =', hkl_ref(:)
          write(ipr,'(/A//)') ' Q cannot be parallel to the reference axis for a non-specular reflexion !'
        end do
        stop
      endif

      if( Q_par_n .and. .not. ( &
        ( Alfa_beta .and. Phi_fixed .and. &
             ( ( Naz_fixed .or. Qaz_fixed ) .or. ( Chi_fixed .and. ( abs( sin(chi) ) < eps10 .or. abs( cos(chi) ) < eps10 ) ) ) ) &
        .or. ( Nu_fixed .and. abs( sin(nu) ) < eps10 .and.  Mu_fixed .and. abs( sin(mu) ) < eps10 .and.  Psi_fixed ) &
        .or. ( Delta_fixed .and. abs( sin(delta) ) < eps10 .and.  Eta_fixed .and. abs( sin(eta) ) < eps10 &
                                                                                                  .and.  Psi_fixed ) ) ) then
        call write_error
        do ipr = 6,9,3
          if( ipr == 3 .and. icheck == 0 ) cycle
          write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
          write(ipr,'(/a15,3i3)') ' Reference axis =', hkl_ref(:)
          write(ipr,'(/A)') ' When Q is parallel to the reference axis, the possible operation modes are : '
          write(ipr,'(10x,a15,/10x,a13)') '- "alfa = beta"','- "Phi fixed"'
          write(ipr,'(10x,A)')  '- "Naz fixed" or "Qaz fixed" or "Chi fixed = 0" or "Chi fixed = 90"'
          write(ipr,'(/10x,A)') 'or "Nu fixed = 0", "Psi fixed" and "Mu fixed = 0"  (vertical mode)'
          write(ipr,'(/10x,A//)') 'or "Delta fixed = 0", "Psi fixed" and "Eta fixed = 0"  (horizontal mode)'
        end do
        stop
      endif

! One term of the column reference is known
      if( Specular ) then
        alfa = Theta_Bragg
        beta = Theta_Bragg
      endif
      if( Column_reference ) then
        if( Alfa_fixed ) beta = asin( 2 * sin_t * cos_tau - sin( alfa ) )
        if( Beta_fixed ) alfa = asin( 2 * sin_t * cos_tau - sin( beta ) )
        if( Alfa_Beta ) then
          alfa = asin( sin_t * cos_tau )
          beta = alfa
        endif
        if( Psi_fixed .and. Ref_spec ) then
  ! When reference axis is not (0 0 1), calculated further on.
          alfa = asin( cos_tau * sin_t - cos_t * sin_tau * cos( psi ) )
          beta = asin( cos_tau * sin_t + cos_t * sin_tau * cos( psi ) )
        elseif( Ref_spec ) then
          if( abs( sin_tau ) < eps10 ) then
            psi = 0._db
          else
            x = ( cos_tau * sin_t - sin( alfa ) ) / ( sin_tau * cos_t )
            if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'cos(psi)')
            psi = acos( x )
          endif
        elseif( .not. Psi_fixed ) then
          psi = 0._db
        endif
      endif

! One of the detector angles is given
      if( Delta_fixed ) then

        x = cos( delta )
        if( abs( x ) > eps10 ) then
          cos_n = cos( 2 * Theta_Bragg ) / x
          if( abs( cos_n ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,cos_n,'cos(nu) ')
        else
          call write_error
          do ipr = 6,9,3
            if( ipr == 3 .and. icheck == 0 ) cycle
            write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
            write(ipr,'(/A//)') '   eta cannot be fixed at +/- 90 !'
          end do
        endif
        nu = acos( cos_n )

      elseif( Nu_fixed ) then

        x = cos( nu )
        if( abs( x ) > eps10 ) then
          cos_d = cos( 2 * Theta_Bragg ) / x
          if( abs( cos_d ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,cos_d,'cos(del)')
        else
          call write_error
          do ipr = 6,9,3
            if( ipr == 3 .and. icheck == 0 ) cycle
            write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
            write(ipr,'(/A//)') '    nu cannot be fixed at +/- 90 !'
          end do
        endif
        delta = acos( cos_d )

      elseif( Qaz_fixed ) then

        delta = asin( 2 * sin_t * cos_t * sin( Qaz ) )
        if( Specular ) eta = delta / 2
        if( abs( delta ) < eps10 ) then
          nu = 2 * Theta_Bragg
        elseif( abs( Qaz ) - pi < eps10 ) then
          nu = 0._db
        else 
          nu = asin( tan( delta ) / tan( Qaz ) )
        endif

      endif

      if( Delta_fixed .or. Nu_fixed ) then
        if( abs( nu ) < eps10 ) then
          Qaz = pi / 2
        else
          Qaz = atan2( sin( delta ), cos( delta ) * sin( nu ) )
        endif
      endif

      if( Column_reference .and. ( Delta_fixed .or. nu_fixed .or. Qaz_fixed ) ) then

        if( .not. ( Psi_fixed .and. .not. Ref_spec ) ) then
          x = ( cos_tau - sin( alfa ) * sin_t ) / ( cos( alfa ) * cos_t )
          if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'cos(Naz)')
          Naz = Qaz + acos( x )
        endif

      elseif( Column_reference .and. Naz_fixed ) then

        x = ( cos_tau - sin( alfa ) * sin_t ) / ( cos( alfa ) * cos_t )
        if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'cos(Qaz)')
        Qaz = Naz + acos( x )

        delta = asin( 2 * sin_t * cos_t * sin( Qaz ) )
        if( abs( delta ) > eps10 ) then
          nu = asin( tan( delta ) / tan( Qaz ) )
        else
          nu = 2 * Theta_Bragg
        endif

      endif

      if( ( Column_detector .or. Column_reference ) .and. .not. Q_par_n ) then
! Allways known
        Mat_Nphi(:,1) = H_phi_n(:)
        call Prodvec(Vec,H_phi_n,N_phi)
        x = sqrt( sum( Vec(:)**2 ) )
        Mat_Nphi(:,3) = Vec(:) / x
        call Prodvec(W,Vec,H_phi_n)
        x = sqrt( sum( W(:)**2 ) )
        Mat_Nphi(:,2) = W(:) / x
      endif

      if( Column_detector .and. Column_reference .and. .not. Q_par_n ) then

! Q_L is in the laboratory basis after rotation
        Q_L(1) = cos_t * sin( Qaz )
        Q_L(2) = - sin_t
        Q_L(3) = cos_t * cos( Qaz )

        if( Psi_fixed .or. .not. Ref_spec ) then
          N_L(1) = cos_tau_ref * Q_L(1) - sin_tau_ref * sin( psi ) * Q_L(3)
          N_L(2) = cos_tau_ref * Q_L(2) + sin_tau_ref * cos( psi )
          N_L(3) = cos_tau_ref * Q_L(3) + sin_tau_ref * sin( psi ) * Q_L(1)
          x = sqrt( sum( N_L(:)**2 ) )
          N_L(:) = N_L(:) / x
       else
          N_L(1) = cos( alfa ) * sin( Naz )
          N_L(2) = - sin( alfa )
          N_L(3) = cos( alfa ) * cos( Naz )
        endif

        Mat_NL(:,1) = Q_L(:)
        call Prodvec( Vec, Q_L, N_L )
        x = sqrt( sum( Vec(:)**2 ) )
        Mat_NL(:,3) = Vec(:) / x
        call Prodvec( W, Vec, Q_L )
        x = sqrt( sum( W(:)**2 ) )
        Mat_NL(:,2) = W(:) / x

        call invermat(Mat_Nphi,Mat )
        Mat_NN = Matmul( Mat_NL, Mat )

        if( Mu_fixed .or. Mu_half_nu ) then

          if( Mu_half_nu ) mu = nu / 2

          call Rotation_Matrix(-mu,1,Mat)
          Mat = Matmul( Mat, Mat_NN )

          phi = atan2( Mat(3,2), Mat(3,1) )
 ! Correction in order it works
          phi = phi + pi
          eta = atan2( - Mat(2,3), Mat(1,3) )
! When chi is between - pi/2 and pi/2
!          chi = atan2( sqrt( Mat(3,1)**2 + Mat(3,2)**2 ), Mat(3,3) )
! When chi is between 0 and pi
          chi = acos( Mat(3,3) )

        elseif( Phi_fixed ) then

          call Rotation_Matrix(phi,3,Mat)
          Mat = Matmul( Mat_NN, Mat )

          eta = atan2( Mat(1,2), sqrt( Mat(2,2)**2 + Mat(3,2)**2 ) )
          mu = atan2( Mat(3,2), Mat(2,2) )
! When chi is between - pi/2 and pi/2
!          chi = atan2( Mat(1,3), Mat(1,1) )
! When chi is between 0 and pi
          chi = acos( Mat(1,1) / cos( eta ) )

        elseif( Eta_fixed .or. Chi_fixed .or. Eta_half_delta ) then

           Mat(:,:) = Mat_NN(:,:)

           if( Eta_half_delta .or. Eta_fixed ) then

             if( Eta_half_delta ) eta = delta / 2
             cos_e = cos( eta )
             sin_e = sin( eta )
             if( abs( cos_e ) < eps10 ) then
 ! chi and mu play the same role
               chi = 0._db
             else
 ! When chi is between - pi/2 and pi/2
               sin_c = Mat(1,3) / cos_e
               if( abs( sin_c ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,sin_c,'sin(chi)')
               chi = asin( sin_c )
             endif
             cos_c = cos( chi )
             sin_c = sin( chi )

           else

             cos_c = cos( chi )
             sin_c = sin( chi )
             if( abs( sin_c ) < eps10 ) then
  ! eta and phi play the same role
               eta = 0._db
             else
               cos_e = Mat(1,3) / sin_c
               if( abs( cos_e ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,cos_e,'cos(eta)')
               eta = acos( cos_e )
             endif
             cos_e = cos( eta )
             sin_e = sin( eta )

           endif

           a = sin_e * sin_c
           mu  = atan2( a * Mat(3,3) + cos_c * Mat(2,3), - cos_c * Mat(3,3) + a * Mat(2,3) )
           mu = mu + pi

           a = cos_e * cos_c
           phi = atan2( a * Mat(1,2) - sin_e * Mat(1,1), sin_e * Mat(1,2) + a * Mat(1,1) )

        endif

      endif

! 1 reference and 2 sample angles ------------------------------------------------------------------------------------------

      if( Column_reference .and. .not. Column_detector .and. .not. Q_par_n ) then

        call Rotation_Matrix(Theta_Bragg,3,Mat)

        call Rotation_Matrix(-psi,1,Mat_p)

        Mat = matmul( Mat_p, Mat )

        Mat = matmul( Mat_Nphi, Mat )

        if( Phi_fixed .and. Chi_fixed ) then

          call Rotation_Matrix(-phi,3,Mat_p)

          Mat = matmul( Mat_p, Mat )

          call Rotation_Matrix(chi,2,Mat_p)

          Mat = matmul( Mat_p, Mat )

          eta = atan2( - Mat(1,2), Mat(2,2) )
          mu = asin( - Mat(3,2) )

          ksi = atan2( - Mat(3,1), Mat(3,3) )

        elseif( Phi_fixed .and. ( Mu_fixed .or. Eta_fixed ) ) then

          call Rotation_Matrix(-phi,3,Mat_p)

          Mat = matmul( Mat_p, Mat )

          if( Mu_fixed ) then
            cos_m = cos( mu );   sin_m = sin( mu )
            if( abs( cos_m ) < eps10 ) then
              chi = - atan2( Mat(1,2), Mat(3,2) )
              cos_c = cos( chi )
              a = acos( Mat(1,1) / cos_c )
              b = acos( Mat(1,3) / cos_c )
              eta = ( a - b ) / 2
              ksi = ( a + b ) / 2
            else
              cos_e = Mat(2,2) / cos_m
              eta = acos( cos_e )
              sin_e = sin( eta )

              a = - cos_m * sin_e
              b = sin_m
              chi = atan2( b * Mat(1,2) + a * Mat(3,2), a * Mat(1,2) - b * Mat(3,2) )
              sin_c = sin( chi )
              cos_c = cos( chi )

              a = sin_c * cos_e
              b = sin_c * sin_e * sin_m - cos_c * cos_m
              ksi = atan2( a * Mat(3,3) + b * Mat(3,1), a * Mat(3,1) - b * Mat(3,3) )

            endif
          else
            cos_e = cos( eta );  sin_e = sin( eta )
            if( abs( cos_e ) < eps10 ) then
              ksi = acos( Mat(2,1) / sin_e )
              a = - asin( Mat(3,2) )
 ! mu + chi = a when sin( eta ) = 1 and mu - chi = a when sin( eta ) = - 1
 ! Movement is horizontal, chi and mu have the same role
              mu = a ! ( = +/- alfa )
              chi = 0._db
            else
              cos_m = Mat(2,2) / cos_e
              if( abs( cos_m ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,cos_m,'cos(mu) ')
              mu = acos( cos_m )
              sin_m = sin( mu )

              a = - cos_m * sin_e
              b = sin_m
              chi = atan2( b * Mat(1,2) + a * Mat(3,2), a * Mat(1,2) - b * Mat(3,2) )
              sin_c = sin( chi )
              cos_c = cos( chi )

              a = sin_c * cos_e
              b = sin_c * sin_e * sin_m - cos_c * cos_m
              ksi = atan2( a * Mat(3,3) + b * Mat(3,1), a * Mat(3,1) - b * Mat(3,3) )

            endif
          endif

        else   ! phi is not fixed

          if( Eta_fixed .and. Mu_fixed ) then

            cos_e = cos( eta );  sin_e = sin( eta )
            cos_m = cos( mu );   sin_m = sin( mu )

            a = sin_m**2 + ( sin_e * cos_m )**2
            b = Mat(3,2) * sin_e * cos_m
            c = Mat(3,2)**2 - sin_m**2
            x = b**2 - a*c

            if( x < 0._db ) then
              call write_error
              do ipr = 6,9,3
                if( ipr == 3 .and. icheck == 0 ) cycle
                write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
                write(ipr,125) x
              end do
              stop
            endif

            sin_c = ( - b + sqrt( x ) ) / a
            if( abs(sin_c) > 1._db ) sin_c = ( - b - sqrt( x ) ) / a
            chi = asin( sin_c )
            cos_c = cos( chi )

            if( abs( Mat(3,2) + sin_c * sin_e * cos_m + cos_c * sin_m ) > eps10 ) then
              chi = pi - chi
              cos_c = - cos_c
            endif

          elseif( Eta_fixed .and. Chi_fixed ) then

            cos_c = cos( chi );  sin_c = sin( chi )
            cos_e = cos( eta );  sin_e = sin( eta )

            if( abs( cos_c ) < eps10 .and. abs( sin_e ) < eps10 ) then
              call write_error
              do ipr = 6,9,3
              if( ipr == 3 .and. icheck == 0 ) cycle
                write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
                write(ipr,'(/A//)') '   chi = +/- 90 with eta = 0, is not possible !'
              end do
              stop

            elseif( abs( sin_c ) < eps10 .or. abs( sin_e ) < eps10 ) then

              mu = - asin( Mat(3,2) / cos_c )
              cos_m = cos( mu )

              if( abs( cos_m ) < eps10 ) then
                call write_error
                do ipr = 6,9,3
                  if( ipr == 3 .and. icheck == 0 ) cycle
                  write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
                  write(ipr,'(/A//)') '   One gets, mu = +/- 90 with chi = 0 or eta = 0, what is not possible !'
                end do
                stop
              endif
              
              sin_m = sin( mu )

            elseif( abs( cos_c ) < eps10 ) then

              cos_m = - Mat(3,1) / ( sin_c * sin_e )
              if( abs( cos_m ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,cos_m,'cos(mu) ')
              mu = acos( cos_m )
              sin_m = sin( mu )

            else

              a = cos_c**2 + ( sin_c * sin_e )**2
              b = Mat(3,2) * cos_c
              c = Mat(3,2)**2 - ( sin_c * sin_e )**2
              x = b**2 - a*c

              if( x < 0._db ) then
                call write_error
                do ipr = 6,9,3
                  if( ipr == 3 .and. icheck == 0 ) cycle
                  write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
                  write(ipr,125) x
                end do
                stop
              endif

              sin_m = ( - b + sqrt( x ) ) / a
              if( abs(sin_m) > 1._db ) sin_m = ( - b - sqrt( x ) ) / a
              mu = asin( sin_m )
              cos_m = cos( mu )

              if( abs( Mat(3,2) + sin_c * sin_e * cos_m + cos_c * sin_m ) > eps10 ) then
                mu = pi - mu
                cos_m = - cos_m
              endif

            endif

          elseif( Mu_fixed .and. Chi_fixed ) then

            cos_c = cos( chi );  sin_c = sin( chi )
            cos_m = cos( mu );   sin_m = sin( mu )

            if( abs( sin_c ) < eps10 .or. abs( cos_m ) < eps10 ) then
!            if( abs( sin_c ) < eps10 .or. abs( sin_m ) < eps10 .or. abs( cos_c ) < eps10 .or. abs( cos_m ) < eps10 ) then
              call write_error
              do ipr = 6,9,3
              if( ipr == 3 .and. icheck == 0 ) cycle
                write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
                write(ipr,'(/A//)') &
                 '   When the reference is known, chi and mu cannot be both fixed with chi = 0 and mu = 90 degrees !'
              end do
              stop
            endif

            sin_e = - ( Mat(3,2) + cos_c * sin_m ) / ( sin_c * cos_m )
            if( abs( sin_e ) > 1._db) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,sin_e,'sin(eta)')
            eta = asin( sin_e )
            cos_e = cos( eta )

          endif

          a = sin_c * sin_m - cos_m * cos_c * sin_e
          b = cos_m * cos_e
 ! atan2 is better than atan to get value between - pi and pi and avoid problem with denominator zero.
          phi = atan2( a * Mat(2,2) - b * Mat(1,2), a * Mat(1,2) + b * Mat(2,2) )

          a = sin_c * cos_e
          b = sin_c * sin_e * sin_m - cos_c * cos_m
          ksi = atan2( Mat(3,3) * a + Mat(3,1) * b, Mat(3,1) * a - Mat(3,3) * b )

        endif  ! end of cases without phi fixed

        Qaz = ksi + pi / 2

        delta = asin( 2 * sin_t * cos_t * sin( Qaz ) )
        if( abs( delta ) < eps10 ) then
          nu = 2 * Theta_Bragg
        elseif( abs( abs( Qaz ) - pi ) < eps10 ) then
          nu = 0._db
        else
          x = tan( delta ) / tan( Qaz )
          if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(nu) ')
          nu = asin( x )
        endif

      endif

! 1 detector known and 2 sample angles ------------------------------------------------------------------------------------------

      if( Detector_known .and. .not. Column_reference ) then

! Qaz is necessarily known
        ksi = Qaz - pi / 2
        cos_k = cos( ksi )
        sin_k = sin( ksi )

        if( Mu_fixed .and. Eta_fixed ) then

          call Rotation_Matrix(-Theta_Bragg,3,Mat)
          call Rotation_Matrix(ksi,2,Mat_p)

          Mat = matmul( Mat_p, Mat )

          call Rotation_Matrix(-mu,1,Mat_p)

          Mat = matmul( Mat_p, Mat )

          call Rotation_Matrix(eta,3,Mat_p)

          Mat = matmul( Mat_p, Mat )

          x = - Mat(2,1) / sqrt( Mat_Nphi(1,1)**2 + Mat_Nphi(2,1)**2 )
          if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(phi)')
          phi = asin( x ) + atan2( Mat_Nphi(2,1), Mat_Nphi(1,1) )

          a = Mat_Nphi(3,1)
          b = Mat_Nphi(1,1) * cos( phi ) + Mat_Nphi(2,1) * sin( phi )
          chi = atan2( a * Mat(1,1) - b * Mat(3,1), a * Mat(3,1) + b * Mat(1,1) )

          a = tan( phi )
          psi = atan2( Mat_Nphi(1,2) * Mat(2,3) * a - Mat_Nphi(1,3) * Mat(2,2) * a &
                     + Mat_Nphi(2,3) * Mat(2,2) - Mat_Nphi(2,2) * Mat(2,3), &
                       Mat(2,3) * ( Mat_Nphi(1,3) * a - Mat_Nphi(2,3) ) + Mat(2,2) * ( Mat_Nphi(1,2) * a + Mat_Nphi(2,2) ) )

        elseif( Phi_fixed .and. Chi_fixed ) then

          call Rotation_Matrix(-phi,3,Mat_p)

          Mat = matmul( Mat_p, Mat_Nphi )

          call Rotation_Matrix(chi,2,Mat_p)

          Mat = matmul( Mat_p, Mat )

! Wrong formula 64 in You's paper
!          mu = asin( Mat(3,1) / sqrt( ( sin_k * cos_t )**2 + sin_t**2 ) ) + atan( sin_k * cos_t / sin_t )

          a = sin_t**2 + ( sin_k * cos_t )**2
          b = sin_k * cos_t * Mat(3,1)
          c = Mat(3,1)**2 - sin_t**2
          x = b**2 - a * c
          if( x < 0._db ) then
            call write_error
            do ipr = 6,9,3
              if( ipr == 3 .and. icheck == 0 ) cycle
              write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
              write(ipr,125) x
            end do
            stop
          endif
          cos_m = ( - b + sqrt( x ) ) / a
          if( abs( cos_m ) > 1._db  ) cos_m = ( - b - sqrt( x ) ) / a
          if( abs( cos_m ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,cos_m,'cos(mu) ')
          mu = acos( cos_m )

          sin_m = sin( mu )

          if( abs( Mat(3,1) + cos_m * sin_k * cos_t - sin_m * sin_t ) > eps10 ) then
            mu = - mu
            sin_m = - sin_m
          endif

          a = cos_t * cos_k
          b = cos_t * sin_m * sin_k  + cos_m * sin_t
          eta = atan2( a * Mat(2,1) + b * Mat(1,1), a * Mat(1,1) - b * Mat(2,1) )

          a = cos_m * sin_k * sin_t + sin_m * cos_t
          b = cos_m * cos_k
          psi = atan2( - a * Mat(3,3) - b * Mat(3,2), a * Mat(3,2) - b * Mat(3,3) )

        elseif( Mu_fixed .and. ( Phi_fixed .or. Chi_fixed ) ) then

          call Rotation_Matrix(-Theta_Bragg,3,Mat)
          call Rotation_Matrix(ksi,2,Mat_p)

          Mat = matmul( Mat_p, Mat )

          call Rotation_Matrix(-mu,1,Mat_p)

          Mat = matmul( Mat_p, Mat )

          if( Chi_fixed ) then
            sin_c = sin( chi )
            cos_c = cos( chi )
            a = - sin_c * Mat_Nphi(2,1)
            b = - sin_c * Mat_Nphi(1,1)
            c = Mat(3,1) - cos_c * Mat_Nphi(3,1)
          else
            sin_p = sin( phi )
            cos_p = cos( phi )
            a = - cos_p * Mat_Nphi(1,1) - sin_p * Mat_Nphi(2,1)
            b = Mat_Nphi(3,1)
            c = Mat(3,1)
          endif

          ni = 10
          da = 2 * pi / ni
          Angle = 0._db
          Below = b < c
          do j = 1,15
            do i = 1,ni
              Angle = Angle + da
 ! .EQV. instead of == is mandatory for logical variables in some compilers
              if( a * sin( Angle ) + b * cos( Angle ) < c .EQV. Below ) cycle
              Angle = Angle - da
              da = da / ni
              exit
            end do
            if( i == ni + 1 ) then
              Angle = Angle - ni * da
              da = da / 10
              ni = 10 * ni
            endif
          end do

          if( Chi_fixed ) then
            phi = Angle
            sin_p = sin( phi )
            cos_p = cos( phi )
          else
            chi = Angle
            sin_c = sin( chi )
            cos_c = cos( chi )
          endif

          a = cos_c * ( cos_p * Mat_Nphi(1,1) + sin_p * Mat_Nphi(2,1) ) + sin_c * Mat_Nphi(3,1)
          b = - sin_p * Mat_Nphi(1,1) + cos_p * Mat_Nphi(2,1)

          eta = atan2( b * Mat(1,1) - a * Mat(2,1), a * Mat(1,1) + b * Mat(2,1) )

! One calculates Psi^-1 =  N_phi^-1 Phi^-1 Chi^-1 Eta^-1 Mat     ( Mat = Mu^-1 F Theta )
          call Rotation_Matrix(eta,3,Mat_p)
          Mat = matmul( Mat_p, Mat )
          call Rotation_Matrix(-chi,2,Mat_p)
          Mat = matmul( Mat_p, Mat )
          call Rotation_Matrix(phi,3,Mat_p)
          Mat = matmul( Mat_p, Mat )
          call invermat(Mat_Nphi,Mat_p )
          Mat = matmul( Mat_p, Mat )
          psi = atan2( Mat(2,3), Mat(2,2) )

        endif

        x = cos_tau * sin_t - cos_t * sin_tau * cos( psi )
        if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(alf)')
        alfa = asin( x )
        x = cos_tau * sin_t + cos_t * sin_tau * cos( psi )
        if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(bet)')
        beta = asin( x )

      endif

! 3 sample angles known ------------------------------------------------------------------------------------

      if( .not. ( Column_detector .or. Column_reference ) ) then

        if( .not. Eta_fixed ) then

          cos_m = cos( mu );  sin_m = sin( mu )
          cos_c = cos( chi ); sin_c = sin( chi )
          cos_p = cos( phi ); sin_p = sin( phi )

          a = cos_m * ( - cos_c * ( cos_p * H_phi_n(1) + sin_p * H_phi_n(2) ) - sin_c * H_phi_n(3) )
          b = cos_m * ( - sin_p * H_phi_n(1) + cos_p * H_phi_n(2) )
          c = - sin_t - sin_m * sin_c * ( cos_p * H_phi_n(1) + sin_p * H_phi_n(2) ) + sin_m * cos_c * H_phi_n(3)

          x = c / sqrt( a**2 + b**2 )
          if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(eta)')
          eta = asin( x ) - atan2( b, a )
          cos_e = cos( eta ); sin_e = sin( eta )

        elseif( .not. Mu_fixed ) then

          cos_e = cos( eta ); sin_e = sin( eta )
          cos_c = cos( chi ); sin_c = sin( chi )
          cos_p = cos( phi ); sin_p = sin( phi )

          a = sin_c * ( cos_p * H_phi_n(1) + sin_p * H_phi_n(2) ) - cos_c * H_phi_n(3)
          b = ( - sin_e * cos_c * cos_p - cos_e * sin_p ) * H_phi_n(1) &
            + ( - sin_e * cos_c * sin_p + cos_e * cos_p ) * H_phi_n(2) &
            - sin_e * sin_c * H_phi_n(3)
          c = - sin_t

          x = c / sqrt( a**2 + b**2 )
          if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(mu) ')
          mu = asin( x ) - atan2( b, a )
          cos_m = cos( mu ); sin_m = sin( mu )

        elseif( .not. Chi_fixed ) then

          cos_m = cos( mu );  sin_m = sin( mu )
          cos_e = cos( eta ); sin_e = sin( eta )
          cos_p = cos( phi ); sin_p = sin( phi )

          a = sin_m * ( cos_p * H_phi_n(1) + sin_p * H_phi_n(2) ) - cos_m * sin_e * H_phi_n(3)
          b = - cos_m * sin_e * ( cos_p * H_phi_n(1) + sin_p * H_phi_n(2) ) - sin_m * H_phi_n(3)
          c = - sin_t + cos_m * cos_e * ( sin_p * H_phi_n(1) - cos_p * H_phi_n(2) )

          x = c / sqrt( a**2 + b**2 )
          if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(chi)')
          chi = asin( x ) - atan2( b, a )
          cos_c = cos( chi ); sin_c = sin( chi )

        else

          cos_m = cos( mu );  sin_m = sin( mu )
          cos_e = cos( eta ); sin_e = sin( eta )
          cos_c = cos( chi ); sin_c = sin( chi )

          a = - cos_m * cos_e * H_phi_n(1) + ( - cos_m * sin_e * cos_c + sin_m * sin_c ) * H_phi_n(2)
          b = cos_m * cos_e * H_phi_n(2) + ( - cos_m * sin_e * cos_c + sin_m * sin_c ) * H_phi_n(1)
          c = - sin_t + ( cos_m * sin_e * sin_c + sin_m * cos_c ) * H_phi_n(3)

          x = c / sqrt( a**2 + b**2 )
          if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(phi)')
          phi = asin( x ) - atan2( b, a )
          cos_p = cos( phi ); sin_p = sin( phi )

        endif

 ! Detector angles
        a = ( cos_e * cos_c * cos_p - sin_e * sin_p ) * H_phi_n(1) &
          + ( cos_e * cos_c * sin_p + sin_e * cos_p ) * H_phi_n(2) &
          + cos_e * sin_c * H_phi_n(3)

        b = ( ( - sin_m * sin_e * cos_c - cos_m * sin_c ) * cos_p - sin_m * cos_e * sin_p ) * H_phi_n(1) &
          + ( ( - sin_m * sin_e * cos_c - cos_m * sin_c ) * sin_p + sin_m * cos_e * cos_p ) * H_phi_n(2) &
          + ( - sin_m * sin_e * sin_c + cos_m * cos_c ) * H_phi_n(3)

        Qaz = atan2( a, b )

        delta = asin( 2 * sin_t * cos_t * sin( Qaz ) )
        if( abs( delta ) > eps10 ) then
          nu = asin( tan( delta ) / tan( Qaz ) )
        else
          nu = 2 * Theta_Bragg
        endif

      endif

      if( Q_par_n ) then

        p(:) = Mat_UB_loc(:,3)
        x = sqrt( sum( p(:)**2 ) )
        p(:) = p(:) / x
        b = asin( p(1) )
        a = - atan2( p(2), p(3) )

        if( Chi_fixed .and. Phi_fixed ) then
          if( abs( sin( chi ) ) < eps10 ) then
            mu = alfa
            eta = 0._db
          elseif( abs( cos( chi ) ) < eps10 ) then
            mu = 0.
            eta = alfa
          endif

        elseif( Nu_fixed .and. Mu_fixed .and. Psi_fixed ) then
          if( abs( sin( nu ) ) < eps10 .and. abs( sin( mu ) ) < eps10 ) then
            chi = pi / 2 - b
            eta = Theta_Bragg - a
            phi = 0._db
          endif

        elseif( Delta_fixed .and. Eta_fixed .and. Psi_fixed ) then
          if( abs( sin( delta ) ) < eps10 .and. abs( sin( eta ) ) < eps10 ) then
            chi = - b
            mu = Theta_Bragg - a
            phi = 0._db
          endif

        endif
      endif

! Final matrix = Z = M H X Sum_T_Bragg_f

      call Rotation_Matrix(-phi,3,Mat)
      call Rotation_Matrix(chi,2,Mat_p)
      Mat = matmul( Mat_p, Mat )
      call Rotation_Matrix(-eta,3,Mat_p)
      Mat = matmul( Mat_p, Mat )
      call Rotation_Matrix(mu,1,Mat_p)
      Mat = matmul( Mat_p, Mat )

      if( .not. ( Naz_fixed .or. ( Column_detector .and. Column_reference ) .or. Q_par_n  ) ) then
        if( Specular ) then
          Naz = Qaz
        else
! When reference axis is not the normal to the surface, it can be that this angle is related to the azimuth of the reference axis
          Naz = atan2( Mat(1,3), Mat(3,3) )
        endif
      endif

! Reference angles
      if( .not. ( Column_reference .or. Column_detector ) ) then
        x = - Mat(2,3)
        if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(alf)')
        alfa = asin( x )
        x = Mat(2,3) + 2 * sin_t * cos_tau
        if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(bet)')
        beta = asin( x )
        x = ( cos_tau * sin_t + Mat(2,3) ) / ( sin_tau * cos_t )
        if( abs( x ) > 1._db ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'sin(psi)')
        psi = acos( x )
      endif

      call invermat( Mat, Mat_p )

      if( Q_par_n .and. Phi_fixed .and. ( Naz_fixed .or. Qaz_fixed ) ) then

        k_i(1) = 0._db
        k_i(2) = cos( alfa )
        k_i(3) = - sin( alfa )

        p_i_s(1) = -1._db
        p_i_s(2:3) = 0._db

        p_i_p(1) = 0._db
        p_i_p(2) = sin( alfa )
        p_i_p(3) = cos( alfa )

        call Rotation_Matrix(phi,3,Mat)

        k_i = matmul( Mat, k_i )
        p_i_s = matmul( Mat, p_i_s )
        p_i_p = matmul( Mat, p_i_p )

        p_i(:) = cos( Naz ) * p_i_p(:) + sin( Naz ) * p_i_s(:)

      else

        k_i(:) = Mat_p(:,2)
        p_i(:) = Mat_p(:,3)

        p_i_s(:) = p_i(:)
        call prodvec(p_i_p,k_i,p_i_s)

      endif

      k_s(:) = H_phi(:) / konde + k_i(:)

      if( Column_detector .and. Column_reference .and. Psi_fixed ) then
        if( .not. Specular ) then
          alfa = - asin( k_i(3) )
          beta = asin( k_s(3) )

          alfa = - asin( dot_product(k_i,Norm_phi) )
          beta = asin( dot_product(k_s,Norm_phi) )

        endif

        if( .not. Naz_fixed .and. .not. Ref_spec ) Naz = Qaz
        if( .not. Naz_fixed .and. .not. Ref_spec ) then
          if( Specular ) then
            Naz = Qaz
          else
            x = ( cos_tau - sin( alfa ) * sin_t ) / ( cos( alfa ) * cos_t )
            if( abs( x ) > 1+eps10 ) call Error_angle(icheck,jpl,Q,Operation_m,Angle_m,x,'cos(Naz)')
! Numerical problem
            if( x > 1._db ) x = x / abs(x)
            Naz = Qaz + acos( x )
          endif
        endif

      endif

      if( Circular ) then
        call prodvec( p_i_p, p_i, k_i )
        Poldafse(:,ipl,1) = cmplx( p_i(:), p_i_p(:) ) / sqrt_2
        Poldafse(:,lpl,1) = cmplx( p_i(:), p_i_p(:) ) / sqrt_2
        Poldafse(:,mpl,1) = cmplx( p_i(:), - p_i_p(:) ) / sqrt_2
        Poldafse(:,npl,1) = cmplx( p_i(:), - p_i_p(:) ) / sqrt_2
        Vecdafse(:,ipl,1) = k_i(:)
        Vecdafss(:,ipl,1) = k_s(:)
        Vecdafse(:,lpl,1) = k_i(:)
        Vecdafss(:,lpl,1) = k_s(:)
        Vecdafse(:,mpl,1) = k_i(:)
        Vecdafss(:,mpl,1) = k_s(:)
        Vecdafse(:,npl,1) = k_i(:)
        Vecdafss(:,npl,1) = k_s(:)
      else
        Poldafse(:,ipl,1) = cmplx( p_i(:), 0._db )
        Poldafse(:,lpl,1) = cmplx( p_i(:), 0._db )
        Vecdafse(:,ipl,1) = k_i(:)
        Vecdafss(:,ipl,1) = k_s(:)
        Vecdafse(:,lpl,1) = k_i(:)
        Vecdafss(:,lpl,1) = k_s(:)
      endif
      
      if( Film .or. Bulk ) then
        Vec(1:2) = 0._db; Vec(3) = 1._db
        call prodvec(p_s_s,Vec,k_s)
        a = sqrt( sum( p_s_s(:)**2 ) )
        if( a < eps10 ) then
          p_s_s(1) = 1._db; p_s_s(2:3) = 0._db 
        else
          p_s_s(:) = p_s_s(:) / a
        endif
      else
        if( Q_mod(ipl) > eps10 ) then
          call prodvec(p_s_s,H_phi,k_s)
          a = sqrt( sum( p_s_s(:)**2 ) )
          p_s_s(:) = p_s_s(:) / a
        else
          p_s_s(:) = p_i(:)
        endif
      endif

      call prodvec(p_s_p,k_s,p_s_s)

! The following establishment of isigpi and angpoldafs is for the Col_dafs_name routine
      if( Film .or. Bulk ) then
        cos_p = abs( sum( p_i(:) * Vec(:) ) )
        call prodvec(W,k_i,p_i_s)
        a = sqrt( sum( ( W(:) - p_i_p(:) )**2 ) )
        Pol_pi = a < eps10
      else
        Pol_pi = .false.
        if( Q_mod(ipl) > eps10 ) then
          cos_p = abs( sum( p_i(:) * H_phi(:) ) ) / Q_mod(ipl)
        else
          cos_p = 0._db
        endif
      endif
      if( Circular ) then
        isigpi(1,ipl) = 3
        isigpi(1,lpl) = 3
        isigpi(1,mpl) = 4
        isigpi(1,npl) = 4
      elseif( cos_p < eps10 ) then
        isigpi(1,ipl) = 1
        isigpi(1,lpl) = 1
      elseif( ( ( Film .or. Bulk) .and. Pol_pi ) .or. ( .not. ( Film .or. Bulk) .and. abs( cos_p - cos_t ) < eps10 ) ) then
        isigpi(1,ipl) = 2
        isigpi(1,lpl) = 2
      else
        isigpi(1,ipl) = 5
        isigpi(1,lpl) = 5
      endif

      if( .not. ( Film .or. Bulk ) .or. ( ( Film .or. Bulk ) .and. Specular ) ) then
        isigpi(2,ipl) = 1
      else
        cos_p = abs( sum( p_s_s(:) * Vec(:) ) )
        call prodvec(W,p_s_s,k_s)
        a = sqrt( sum( ( W(:) - p_s_p(:) )**2 ) )
        Pol_pi = a < eps10
        if( cos_p < eps10 ) then
          isigpi(2,ipl) = 1
        elseif( Pol_pi ) then
          isigpi(2,ipl) = 2
        else
          isigpi(2,ipl) = 5
        endif
      endif

      if( .not. ( Film .or. Bulk ) .or. ( ( Film .or. Bulk ) .and. Specular ) ) then
        isigpi(2,lpl) = 2
      else
        cos_p = abs( sum( p_s_p(:) * Vec(:) ) )
        if( cos_p < eps10 ) then
          isigpi(2,lpl) = 1
!        elseif( cos_p > 0.99390826_db ) then
         ! cosp corresponds to 2 degres, it is the limit to write pi polarization... 
        elseif( cos_p > 0.999999_db ) then
          isigpi(2,lpl) = 2
        else
          isigpi(2,lpl) = 5
        endif
      endif
      if( Circular ) then
        isigpi(2,mpl) = isigpi(2,ipl)
        isigpi(2,npl) = isigpi(2,lpl)
      endif

      Poldafss(:,ipl,1) = cmplx( p_s_s(:), 0._db )
      Poldafss(:,lpl,1) = cmplx( p_s_p(:), 0._db )
      if( Circular ) then
        Poldafss(:,mpl,1) = cmplx( p_s_s(:), 0._db )
        Poldafss(:,npl,1) = cmplx( p_s_p(:), 0._db )
      endif
 
      if( Phi_fixed ) then
        angpoldafs(3,ipl) = phi
        angpoldafs(3,lpl) = phi
        if( Circular ) then
          angpoldafs(3,mpl) = phi
          angpoldafs(3,npl) = phi
        endif
      endif

      if( Film .and. alfa > eps10 .and. beta > eps10 ) then
        Length_abs(ipl) = c_cos_z * ( 1 / sin( alfa ) + 1 / sin( beta ) )
      else
        Length_abs(ipl) = 1.e+20_db
      endif
      Length_abs(lpl) = Length_abs(ipl)
      if( Circular ) then
        Length_abs(mpl) = Length_abs(ipl)
        Length_abs(npl) = Length_abs(ipl)
      endif

      if( phi > pi ) phi = phi - 2 * pi

      if( icheck > 0 ) then
        if( kpl == 1 .or. ( icheck > 1 .and. .not. Mat_unit ) ) then
          write(3,'(/A)') ' Orthonormalized Surface basis'
          if( Circular ) then
            write(3,150)
          else
            write(3,155)
          endif
        endif
        if( Q_par_n .and. ( Naz_fixed .or. Qaz_fixed ) ) then
          if( Circular ) then
            write(3,160) ipl, hkl_dafs(:,ipl), Theta_Bragg / radian, alfa / radian, beta / radian, delta / radian, &
                     nu / radian, Qaz / radian, Naz / radian, &
                     psi / radian, phi / radian, Vecdafse(:,ipl,1), Poldafse(:,ipl,1), &
                     Vecdafss(:,ipl,1), real( Poldafss(:,ipl,1), db ), real( Poldafss(:,lpl,1), db )
          else
            write(3,165) ipl, hkl_dafs(:,ipl), Theta_Bragg / radian, alfa / radian, beta / radian, delta / radian, &
                     nu / radian, Qaz / radian, Naz / radian, &
                     psi / radian, phi / radian, Vecdafse(:,ipl,1), real( Poldafse(:,ipl,1), db ), &
                     Vecdafss(:,ipl,1), real( Poldafss(:,ipl,1), db ), real( Poldafss(:,lpl,1), db )
          endif
        else
          if( Circular ) then
            write(3,170) ipl, hkl_dafs(:,ipl), Theta_Bragg / radian, alfa / radian, beta / radian, delta / radian, &
                     nu / radian, Qaz / radian, Naz / radian, &
                     psi / radian, eta / radian, mu / radian, chi / radian, phi / radian, &
                     Vecdafse(:,ipl,1), Poldafse(:,ipl,1), Vecdafss(:,ipl,1), &
                     real( Poldafss(:,ipl,1), db ), real( Poldafss(:,lpl,1), db )
          else
            write(3,175) ipl, hkl_dafs(:,ipl), Theta_Bragg / radian, alfa / radian, beta / radian, delta / radian, &
                     nu / radian, Qaz / radian, Naz / radian, &
                     psi / radian, eta / radian, mu / radian, chi / radian, phi / radian, &
                     Vecdafse(:,ipl,1), real( Poldafse(:,ipl,1), db ), Vecdafss(:,ipl,1), &
                     real( Poldafss(:,ipl,1), db ), real( Poldafss(:,lpl,1), db )
          endif
        endif
      endif

      if( icheck > 1 .and. .not. Mat_unit  ) then
        do i = ipl,npl,npl_2d(jpl)
          p(:) = real( Poldafse(:,i,1), db )
          p = matmul( Mat_K, p )
          Poldafse(:,i,1) = cmplx( p(:), 0._db )
          p(:) = real( Poldafss(:,i,1), db )
          p = matmul( Mat_K, p )
          Poldafss(:,i,1) = cmplx( p(:), 0._db )
          p(:) = Vecdafse(:,i,1)
          p = matmul( Mat_K, p )
          Vecdafse(:,i,1) = cmplx( p(:), 0._db )
          p(:) = Vecdafss(:,i,1)
          p = matmul( Mat_K, p )
          Vecdafss(:,i,1) = cmplx( p(:), 0._db )
        end do

        write(3,'(/A)') ' Internal basis R1, z along c'
        if( Circular ) then
          write(3,150)
        else
          write(3,155)
        endif
        
        if( Q_par_n .and. ( Naz_fixed .or. Qaz_fixed ) ) then
          if( Circular ) then
            write(3,160) ipl, hkl_dafs(:,ipl), Theta_Bragg / radian, alfa / radian, beta / radian, delta / radian, &
                     nu / radian, Qaz / radian, Naz / radian, &
                     psi / radian, phi / radian, Vecdafse(:,ipl,1), Poldafse(:,ipl,1), &
                     Vecdafss(:,ipl,1), real( Poldafss(:,ipl,1), db ), real( Poldafss(:,lpl,1), db )
          else
            write(3,165) ipl, hkl_dafs(:,ipl), Theta_Bragg / radian, alfa / radian, beta / radian, delta / radian, &
                     nu / radian, Qaz / radian, Naz / radian, &
                     psi / radian, phi / radian, Vecdafse(:,ipl,1), real( Poldafse(:,ipl,1), db ), &
                     Vecdafss(:,ipl,1), real( Poldafss(:,ipl,1), db ), real( Poldafss(:,lpl,1), db )
          endif
        else
          if( Circular ) then
            write(3,170) ipl, hkl_dafs(:,ipl), Theta_Bragg / radian, alfa / radian, beta / radian, delta / radian, &
                     nu / radian, Qaz / radian, Naz / radian, &
                     psi / radian, eta / radian, mu / radian, chi / radian, phi / radian, &
                     Vecdafse(:,ipl,1), Poldafse(:,ipl,1), Vecdafss(:,ipl,1), &
                     real( Poldafss(:,ipl,1), db ), real( Poldafss(:,lpl,1), db )
          else
            write(3,175) ipl, hkl_dafs(:,ipl), Theta_Bragg / radian, alfa / radian, beta / radian, delta / radian, &
                     nu / radian, Qaz / radian, Naz / radian, &
                     psi / radian, eta / radian, mu / radian, chi / radian, phi / radian, &
                     Vecdafse(:,ipl,1), real( Poldafse(:,ipl,1), db ), Vecdafss(:,ipl,1), &
                     real( Poldafss(:,ipl,1), db ), real( Poldafss(:,lpl,1), db )
          endif
        endif

        do i = ipl,npl,npl_2d(jpl)
          p(:) = real( Poldafse(:,i,1), db )
          p = matmul( Mat_K_i, p )
          Poldafse(:,i,1) = cmplx( p(:), 0._db )
          p(:) = real( Poldafss(:,i,1), db )
          p = matmul( Mat_K_i, p )
          Poldafss(:,i,1) = cmplx( p(:), 0._db )
          p(:) = Vecdafse(:,i,1)
          p = matmul( Mat_K_i, p )
          Vecdafse(:,i,1) = cmplx( p(:), 0._db )
          p(:) = Vecdafss(:,i,1)
          p = matmul( Mat_K_i, p )
          Vecdafss(:,i,1) = cmplx( p(:), 0._db )
        end do
      endif

! Test on results validity
      x = sqrt( sum( k_s(:)**2 ) )
      if( UB_in .and. Q_par_n .and. abs( x - 1 ) > eps10 .and. abs( x - 1 ) < 1.e-09_db) then
        k_s(:) = k_s(:) / x
        x = 1._db
      endif
      if( abs( x - 1 ) > eps10 .or. alfa < 0._db .or. beta < 0 ) then
        call write_error
        do ipr = 6,9,3
          if( ipr == 3 .and. icheck == 0 ) cycle
          write(ipr,120) jpl, hkl_dafs(:,ipl), Operation_m(:), Angle_m(:) / radian
          if( abs( x - 1 ) > eps10 ) write(ipr,180) x
          if( UB_in .and. Q_par_n ) then
            write(ipr,190)
          elseif( abs( x - 1 ) > eps10 ) then
            write(ipr,'(/A)') '   There is an error in the code, contact the author !'
          endif
          if( alfa < 0._db ) write(ipr,210) '   Incidence angle =', alfa / radian
          if( beta < 0._db ) write(ipr,210) '    Outgoing angle =', beta / radian
        end do
        stop
      endif

    end do

    ipl = npl

  end do
  
  return
 110 format(4(2x,3f12.7))
 120 format(//' For the reflection',i5,', (h,k,l) =',3f8.4,/ &
              '   with the operation mode =',5i2,/, &
              '   the angles',3f8.3,' or the current energy are not possible !')
 125 format(/'   One gets a negative discriminant =',1p,e11.3,','/ &
             '   from wich angle cannot be calculated !'//)
 150 format(/'   ipl    h      k      l     theta_B   alfa    beta   delta    nu     Qaz     Naz     psi', &
             '     eta      mu     chi      phi',16x,'k_i',36x,'Pol_i',36x,'k_s',23x,'Pol_s_1',22x,'Pol_s_2')
 155 format(/'   ipl    h      k      l     theta_B   alfa    beta   delta    nu     Qaz     Naz     psi', &
             '     eta      mu     chi      phi',16x,'k_i',24x,'Pol_i',23x,'k_s',23x,'Pol_s_1',21x,'Pol_s_2')
 160 format(i5,2f7.2,f9.4,8f8.3,4x,'x',7x,'x',7x,'x',3x,f8.3,1x,3f9.5,1x,3(f9.5,f8.4),3(1x,3f9.5))
 165 format(i5,2f7.2,f9.4,8f8.3,4x,'x',7x,'x',7x,'x',3x,f8.3,5(1x,3f9.5))
 170 format(i5,2f7.2,f9.4,12f8.3,1x,3f9.5,1x,3(f9.5,f8.5),3(1x,3f9.5))
 175 format(i5,2f7.2,f9.4,12f8.3,5(1x,3f9.5))
 180 format(/'   Normalized outgoing wave vector modulus =',f14.10,' different from one.')
 190 format(/'   When an experimental UB matrix is used,',/ &
             '   it is better to also use a reference axis not along 0 0 1, which is the default !', /'   ( keyword "setaz" )'//)
 210 format(/a20,f12.7,' < 0 !'//)
end

!*********************************************************************

! The operation modes come from You, J. Appl. Cryst. 32, 614 (1999),
! with the Spec inversion in sample angle order: eta is 1 and mu is 2
! Table is :     Detector  Reference    Sample   Sample   Sample
!              1  delta   alfa = beta     eta      eta      eta
!              2    nu       alfa          mu       mu       mu
!              3   Qaz       beta         chi      chi      chi
!              4   Naz       psi          phi      phi      phi
!              5    x         x      eta = delta/2  x        x
!              6    x         x       mu = nu/2     x        x

  subroutine Operation_mode_eval(alfa,Alfa_beta,Alfa_fixed,Angle_m,beta,Beta_fixed,chi,Chi_fixed,Column_detector, &
                Column_reference,delta,Delta_fixed,eta,Detector_known,Eta_fixed,Eta_half_delta,icheck,mu,Mu_fixed, &
                Mu_half_nu,Naz,Naz_fixed,nu,Nu_fixed,Operation_m,phi,Phi_fixed,psi,Psi_fixed,Qaz,Qaz_fixed)

  use declarations
  implicit none

  character(len=3):: mot3
  character(len=5), dimension(3):: Angle_name

  integer:: i, icheck, ipr, j, n_a, n_op
  integer, dimension(5):: Operation_m

  logical:: Alfa_Beta, Alfa_fixed, Beta_fixed, Chi_fixed, Column_detector, Column_reference, Delta_fixed, Detector_known, &
    Eta_fixed, Eta_half_delta, Mu_fixed, Mu_half_nu, Naz_fixed, Nu_fixed, Operation_m_error, Phi_fixed, Psi_fixed, Qaz_fixed

  real(kind=db):: alfa, beta, chi, delta, eta, mu, Naz, nu, phi, psi, Qaz
  real(kind=db), dimension(3):: Angle_m

  Alfa_Beta = .false.; Alfa_fixed = .false.;     Beta_fixed = .false.; Chi_fixed = .false.;  Delta_fixed = .false.
  Eta_fixed = .false.; Eta_half_delta = .false.; Mu_fixed = .false.;   Mu_half_nu = .false.; Naz_fixed = .false.
  Nu_fixed = .false.;  Phi_fixed = .false.;      Psi_fixed = .false.;  Qaz_fixed = .false.

  alfa = 0._db; beta = 0._db; chi = 0._db; delta = 0._db; eta = 0._db;
  mu = 0._db;    Naz = 0._db;  nu = 0._db;   phi = 0._db; psi = 0._db; Qaz = 0._db

  Operation_m_error = .false.

! j is the column in the You's table
  n_op = 0
  n_a = 0

  do j = 1,5

    if( Operation_m(j) == 0 ) cycle

    n_op = n_op + 1
    n_a = n_a + 1

    if( n_op > 3 ) cycle

    select case(j)

      case(1)

        select case( Operation_m(j) )
          case(1)
            Delta_fixed = .true.
            delta = Angle_m(n_a)
            Angle_name(n_a) = 'delta'
          case(2)
            Nu_fixed = .true.
            nu = Angle_m(n_a)
            Angle_name(n_a) = '   nu'
          case(3)
            Qaz_fixed = .true.
            Qaz = Angle_m(n_a)
            Angle_name(n_a) = '  Qaz'
          case(4)
            Naz_fixed = .true.
            Naz = Angle_m(n_a)
            Angle_name(n_a) = '  Naz'
          case default
            call write_error
            do ipr = 6,9,3
              if( ipr == 3 .and. icheck == 0 ) cycle
              write(ipr,100) Operation_m(:), Angle_m(:) / radian
              write(ipr,120) ' Detector', i, 4
            end do
            stop
        end select

      case(2)

        select case( Operation_m(j) )
          case(1)
            Alfa_beta = .true.
            n_a = n_a - 1
          case(2)
            Alfa_fixed = .true.
            alfa = Angle_m(n_a)
            Angle_name(n_a) = ' alfa'
          case(3)
            Beta_fixed = .true.
            beta = Angle_m(n_a)
            Angle_name(n_a) = ' beta'
          case(4)
            Psi_fixed = .true.
            psi = Angle_m(n_a)
            Angle_name(n_a) = '  psi'
          case default
            call write_error
            do ipr = 6,9,3
              if( ipr == 3 .and. icheck == 0 ) cycle
              write(ipr,100) Operation_m(:), Angle_m(:) / radian
              write(ipr,120) 'Reference', j, Operation_m(j), 4
            end do
            stop
        end select

      case(3,4,5)

        select case( Operation_m(j) )
          case(1)
            Operation_m_error = Eta_fixed
            mot3 = 'eta'
            Eta_fixed = .true.
            eta = Angle_m(n_a)
            Angle_name(n_a) = '  eta'
         case(2)
            Operation_m_error = Mu_fixed
            mot3 = ' mu'
            Mu_fixed = .true.
            mu = Angle_m(n_a)
            Angle_name(n_a) = '   mu'
          case(3)
            Operation_m_error = Chi_fixed
            mot3 = 'chi'
            Chi_fixed = .true.
            chi = Angle_m(n_a)
            Angle_name(n_a) = '  chi'
          case(4)
            Operation_m_error = Phi_fixed
            mot3 = 'phi'
            Phi_fixed = .true.
            phi = Angle_m(n_a)
            Angle_name(n_a) = '  phi'
          case(5)
            if( n_op /= 3 ) then
              call write_error
              do ipr = 6,9,3
                if( ipr == 3 .and. icheck == 0 ) cycle
                write(ipr,100) Operation_m(:), Angle_m(:) / radian
                write(ipr,110) ' eta = delta / 2'
              end do
              stop
            endif
            Operation_m_error = Eta_half_delta
            Eta_half_delta = .true.
          case(6)
            if( n_op /= 3 ) then
              call write_error
              do ipr = 6,9,3
                if( ipr == 3 .and. icheck == 0 ) cycle
                write(ipr,100) Operation_m(:), Angle_m(:) / radian
                write(ipr,110) '     mu = nu / 2'
              end do
              stop
            endif
            Operation_m_error = Mu_half_nu
            Mu_half_nu = .true.
          case default
            call write_error
            do ipr = 6,9,3
              if( ipr == 3 .and. icheck == 0 ) cycle
              write(ipr,100) Operation_m(:), Angle_m(:) / radian
              if( n_op == 3 ) then
                write(ipr,120) '   Sample', j, Operation_m(j), 6
              else
                write(ipr,120) '   Sample', j, Operation_m(j), 4
              endif
            end do
            stop
        end select

        if( Operation_m_error ) then
          call write_error
          do ipr = 6,9,3
            if( ipr == 3 .and. icheck == 0 ) cycle
            write(ipr,100) Operation_m(:), Angle_m(:) / radian
            write(ipr,130) mot3
          end do
          stop
        endif

    end select

  end do

  if( n_op /= 3 ) then
    call write_error
    do ipr = 6,9,3
      if( ipr == 3 .and. icheck == 0 ) cycle
      write(ipr,100) Operation_m(:), Angle_m(:) / radian
      write(ipr,140) n_op
    end do
    stop
  endif

  Column_reference = Alfa_fixed .or. Beta_fixed .or. Alfa_beta .or. Psi_fixed
  Detector_known = Delta_fixed .or. Nu_fixed .or. Qaz_fixed
  Column_detector = Detector_known .or. Naz_fixed

  if( icheck > 0 ) write(3,145) Operation_m(:), ( Angle_name(i), Angle_m(i) / radian, i = 1,n_a )
  if( icheck > 1 ) then
    write(3,150) Column_detector, Delta_fixed, Nu_fixed, Qaz_fixed, Naz_fixed
    write(3,160) Column_reference, Alfa_beta, Alfa_fixed, Beta_fixed, Psi_fixed
    write(3,170) Eta_fixed, Mu_fixed, Chi_fixed, Phi_fixed
  endif

  if( ( Naz_fixed .and. .not. Column_reference ) .or. &
      ( Detector_known .and. Eta_fixed .and. ( Chi_fixed .or. Phi_fixed ) ) ) then
    call write_error
    do ipr = 6,9,3
      if( ipr == 3 .and. icheck == 0 ) cycle
      write(ipr,100) Operation_m(:), ( Angle_name(i), Angle_m(i) / radian, i = 1,n_a )
      if( Naz_fixed ) then
        write(ipr,'(/A//)') '  Naz fixed with 2 sample angles fixed is not yet programmed, sorry !'
      elseif( Detector_known .and. Eta_fixed .and. Chi_fixed ) then
        write(ipr,'(/A//)') '  Detector position with eta and chi fixed is not yet programmed, sorry !'
      elseif( Detector_known .and. Eta_fixed .and. Phi_fixed ) then
        write(ipr,'(/A//)') '  Detector position with eta and phi fixed is not yet programmed, sorry !'
      endif
    end do
    stop
  endif

  return
100 format(//' The operation mode',5i2,' is not possible !',/ &
             ' Angles:',3(1x,a5,' =',f8.3,','))
110 format(/1x,a16,' is possible only without any other sample angle given !'//)
120 format(/1x,a9,' digit',i2,' is',i2,' > ',i1,/5x,' It is not possible !'//)
130 format(/'  Sample ',a3,' value at least 2 times given what is not possible !')
140 format(/'  Only 3 non-zero terms are possible, and there are',i2,' !'//)
145 format(/'  Operation mode =',5i2,',',3(1x,a5,' =',f8.3,','))
150 format(/'  Column detector  =',L2,2x,' delta_fixed, nu_fixed,  Qaz_fixed, Naz_fixed =',4L2)
160 format( '  Column reference =',L2,2x,' alfa=beta, alfa_fixed, beta_fixed, psi_fixed =',4L2)
170 format( '  Column sample',9x,        ' eta_fixed,   mu_fixed,  chi_fixed, phi_fixed =',4L2)
end

!***********************************************************************

subroutine Rotation_Matrix(angle,k,Mat)

  use declarations
  implicit none

  integer:: i, j, k

  real(kind=db):: angle
  real(kind=db), dimension(3,3):: Mat

  Mat(:,:) = 0._db
  Mat(k,k) = 1._db

  i = 1 + mod(k,3)
  j = 1 + mod(i,3)

  Mat(i,i) = cos( angle ); Mat(i,j) = - sin( angle )
  Mat(j,i) = - Mat(i,j);   Mat(j,j) = Mat(i,i)

  return
end

!***********************************************************************

subroutine Error_angle(icheck,jpl,Q,Operation_m,Angle_m,Argument,Angle_name)

  use declarations
  implicit none

  integer:: icheck, ipr, jpl
  integer, dimension(5):: Operation_m

  character(len=8):: Angle_name

  real(kind=db):: Argument
  real(kind=db), dimension(3):: Angle_m, Q

! Numerical imprecision cases
  if( abs( Argument - 1 ) < eps10 ) then
    Argument = 1._db
    return
  elseif( abs( Argument + 1 ) < eps10 ) then
    Argument = - 1._db
    return
  endif

  call write_error
  do ipr = 6,9,3
    if( ipr == 3 .and. icheck == 0 ) cycle
    write(ipr,110) jpl, Q(:), Operation_m(:), Angle_m(:) / radian, Angle_name, Argument, Angle_name(1:3)
  end do
  stop

  return
 110 format(//' For the reflection',i3,', (h,k,l) =',2f5.2,f7.4,/ &
              '   with the operation mode =',5i2,/, &
              '   the angles',3f8.3,' or the current energy are not possible !'// &
              '   One gets the absolute value of ',a8,' =',f10.3,' > 1,',/&
              '   from wich an arc',a3,' must be taken !'//)
end

!***********************************************************************

! Polarization calculation for DAFS

subroutine Pol_dafs(Angle_or,Angpoldafs,Angxyz,Angxyz_bulk,axyz,axyz_bulk,Bormann,Bulk,Bulk_step,Dafs_bio,Energy_photon, &
            Film,hkl_dafs,hkl_film,icheck,isigpi,Length_abs,Mat_or,mpirank0,nphi_dafs,nphim,npldafs,npldafs_e,npldafs_f,Orthmatt, &
            Orthmati,Poldafse,Poldafsem,Poldafss,Poldafssm,Q_mod,Vec_orig,Vecdafse,Vecdafsem,Vecdafss,Vecdafssm)

  use declarations
  implicit none

  integer:: i, icheck, ip, ipl, ipr, istop, j, mpirank0, nphim, npldafs, npldafs_e, npldafs_f

  complex(kind=db), dimension(3):: pe, ps
  complex(kind=db), dimension(3,npldafs_e):: Poldafsem, Poldafssm
  complex(kind=db), dimension(3,npldafs,nphim):: Poldafse, Poldafss

  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(2,npldafs):: isigpi

  logical:: Bormann, Bulk, Bulk_step, Film, hkl_film, Dafs_bio

  real(kind=db):: Angle, c_cos_z, cos_pe, cos_ps, cos_z, cosb, cosp, detmat, dhkl, dp, Energy_photon, &
    fac, konde, pp, psi, qkn, rad_i, sin_pe, sin_ps, sinb, sinp, Thetabragg, Vol, wn, wsn

  real(kind=db), dimension(3):: ang, angxyz, angxyz_bulk, axyz, axyz_bulk, cosdir, hklred, &
     qi, qj, qk, v, Vec_orig, vpie, vpe, vps, vpis, vsig, vx, vy, vz, w, we, ws, wx, wy, wz
  real(kind=db), dimension(npldafs):: Q_mod
  real(kind=db), dimension(npldafs_f):: Angle_or
  real(kind=db), dimension(0:3):: det
  real(kind=db), dimension(3,3):: Mat, Mat_or, Mat_ori, Orthmatt, Orthmati
  real(kind=db), dimension(3,npldafs):: angpoldafs, hkl_dafs
  real(kind=db), dimension(3,npldafs_e):: Vecdafsem, Vecdafssm
  real(kind=db), dimension(npldafs):: Length_abs
  real(kind=db), dimension(3,npldafs,nphim):: Vecdafse, Vecdafss

! lambda = 2 * pi / k = 2 * d * sintheta
! In S.I. vecond = k = E*alfa*4*pi*epsilon0 / (e*e)
! In ua and rydb : k = 0.5 * alfa * E
  konde = 0.5_db * alfa_sf * Energy_photon

  if( Bulk .and. .not. hkl_film ) then
 ! When bulk_step angxyz_bulk contains film data ! ( See fdm )
    cosdir(:) = cos( angxyz_bulk(:) )
  else
    cosdir(:) = cos( angxyz(:) )
  endif

  if( Film ) then
    if( Bulk_step ) then
      cos_z = sqrt( sin( angxyz(2) )**2 - ( ( cos( angxyz(1) ) - cos( angxyz(3) ) * cos( angxyz(2) ) ) / sin( angxyz(3) ) )**2 )
      c_cos_z = axyz(3) * cos_z
    else
      cos_z = sqrt( sin( angxyz_bulk(2) )**2 - ( ( cos( angxyz_bulk(1) ) - cos( angxyz_bulk(3) ) * cos( angxyz_bulk(2) ) ) &
                                             / sin( angxyz_bulk(3) ) )**2 )
      c_cos_z = axyz_bulk(3) * cos_z
    endif
! The factor is to convert Length_abs en micrometer
    c_cos_z = 0.0001_db * bohr * c_cos_z
  endif

  if( Bormann ) then
    Poldafse(:,:,:) = ( 0._db, 0._db )
    Poldafss(:,:,:) = ( 0._db, 0._db )
    Vecdafse(:,:,:) = 0._db
    Vecdafss(:,:,:) = 0._db
  endif

  istop = 0

  do ipl = 1,npldafs

    if( sum( abs( hkl_dafs(:,ipl) ) ) < eps10 ) then

      dhkl = 0._db
      thetabragg = 0._db

    else

      if( ( Bulk .and. .not. hkl_film ) .or. ( Bulk_step .and. hkl_film ) ) then
        hklred(:) = hkl_dafs(:,ipl) / axyz_bulk(:)
      else
        hklred(:) = hkl_dafs(:,ipl) / axyz(:)
      endif
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
! Inter-reticular distance
      dhkl = sqrt( det(0) / sum( hklred(1:3) * det(1:3) ) )
      fac = abs( pi / ( konde * dhkl ) )

      if( fac > 1._db .and. mpirank0 == 0 ) then
        if( istop == 0 ) call write_error
        do ipr = 6,9,3
          if( ipr == 3 .and. icheck == 0 ) cycle
          if( Film ) then
            write(ipr,120) ipl, hkl_dafs(:,ipl)
          else
            write(ipr,125) ipl, nint( hkl_dafs(:,ipl) )
          endif
        end do
        istop = 1
        cycle
      endif

      Thetabragg = asin( pi / ( konde * dhkl ) )

    endif

! When angpoldafs(3,ipl) = 10000, polarization and wave vector are imposed
    if( angpoldafs(3,ipl) > 9999._db .and. .not. Dafs_bio ) then

      Poldafse(:,ipl,1) = Poldafsem(:,ipl)
      Poldafss(:,ipl,1) = Poldafssm(:,ipl)
      Vecdafse(:,ipl,1) = Vecdafsem(:,ipl)
      Vecdafss(:,ipl,1) = Vecdafssm(:,ipl)

! One goes to the internal orthogonal basis
      v(:) = real( Poldafse(:,ipl,1), db )
      if( sum( v(:)**2 )  > eps6 ) then
        call trvec(mpirank0,Orthmatt,v,w)
        Poldafse(:,ipl,1) = cmplx( w(:), 0._db )
      endif

      v(:) = real( Poldafss(:,ipl,1), db )
      if( sum( v(:)**2 )  > eps6 ) then
        call trvec(mpirank0,Orthmatt,v,w)
        Poldafss(:,ipl,1) = cmplx( w(:), 0._db )
      endif

      v(:) = Vecdafse(:,ipl,1)
      if( sum( v(:)**2 )  > eps6 ) then
        call trvec(mpirank0,Orthmatt,v,w)
        Vecdafse(:,ipl,1) = w(:)
      endif

      v(:) = Vecdafss(:,ipl,1)
      if( sum( v(:)**2 )  > eps6 ) then
        call trvec(mpirank0,Orthmatt,v,w)
        Vecdafss(:,ipl,1) = w(:)
      endif

    elseif( Bormann ) then

      select case(ipl)
        case(1,10,19,28)
          Poldafse(1,ipl,1) = ( 1._db, 0._db )
          Poldafss(1,ipl,1) = ( 1._db, 0._db )
        case(2,11,20,29)
          Poldafse(1,ipl,1) = ( 1._db, 0._db )
          Poldafss(2,ipl,1) = ( 1._db, 0._db )
        case(3,12,21,30)
          Poldafse(1,ipl,1) = ( 1._db, 0._db )
          Poldafss(3,ipl,1) = ( 1._db, 0._db )
        case(4,13,22,31)
          Poldafse(2,ipl,1) = ( 1._db, 0._db )
          Poldafss(1,ipl,1) = ( 1._db, 0._db )
        case(5,14,23,32)
          Poldafse(2,ipl,1) = ( 1._db, 0._db )
          Poldafss(2,ipl,1) = ( 1._db, 0._db )
        case(6,15,24,33)
          Poldafse(2,ipl,1) = ( 1._db, 0._db )
          Poldafss(3,ipl,1) = ( 1._db, 0._db )
        case(7,16,25,34)
          Poldafse(3,ipl,1) = ( 1._db, 0._db )
          Poldafss(1,ipl,1) = ( 1._db, 0._db )
        case(8,17,26,35)
          Poldafse(3,ipl,1) = ( 1._db, 0._db )
          Poldafss(2,ipl,1) = ( 1._db, 0._db )
        case(9,18,27,36)
          Poldafse(3,ipl,1) = ( 1._db, 0._db )
          Poldafss(3,ipl,1) = ( 1._db, 0._db )
      end select

    else

      sinb = sin( thetabragg )
      cosb = cos( thetabragg )

! Diffraction vector in the internal orthogonal basis
      vx(:) = orthmatt(:,1)
      vy(:) = orthmatt(:,2)
      vz(:) = orthmatt(:,3)

! wx, wy, wz : reciprocal cell basis
      call prodvec(wx,vy,vz)

      vol = sum( wx(:) * vx(:) )
      wx(:) = wx(:) / vol
      call prodvec(wy,vz,vx)
      wy(:) = wy(:) / vol
      call prodvec(wz,vx,vy)
      wz(:) = wz(:) / vol

      qk(:) = hkl_dafs(1,ipl) * wx(:) + hkl_dafs(2,ipl) * wy(:) + hkl_dafs(3,ipl) * wz(:)
      qkn = sqrt( sum( qk(:)**2 ) )
      if( abs(qkn) > eps10 ) then
        qk(:) = qk(:) / qkn
      else
! Specular case
        qk(1:2) = 0._db; qk(3) = 1._db
      endif
! The origin of azimuth is along the plane (Q,Vec_orig)
! when not possible, one takes along (Q,Ox) or (Q,Oz).
      do i = 1,3
        select case(i)
          case(1)
            v(:) = Vec_orig(:)
          case(2)
            v(1) = 1._db; v(2:3) = 0._db
          case(3)
            v(1:2) = 0._db; v(3) = 1._db
        end select
        v = matmul( orthmatt, v )
        call prodvec(w,qk,v)
        wn = sqrt( dot_product(w,w) )
        if( wn > eps10 ) exit
      end do
      w(:) = w(:) / wn
      call prodvec(qi,w,qk)
      call prodvec(qj,qk,qi)

      if( icheck > 2 ) then
        if( Film ) then
          write(3,130) hkl_dafs(:,ipl)
        else
          write(3,135) nint( hkl_dafs(:,ipl) )
        endif
        write(3,140) ( wx(i)/bohr, wy(i)/bohr, wz(i)/bohr, i = 1,3 )
        write(3,150) ( qi(i), qj(i), qk(i), i = 1,3 )
      endif

      if( Dafs_bio ) then
! 4 polarizations are calculated : for Q and - Q and both with out sigma and pi.
! In Convolution routine, intensity is calculated for I(Q)_sigma + I(Q)_pi - ( I(-Q)_sigma + I(-Q)_pi )
        if( ipl == 1 ) then
          call invermat(Mat_or,Mat_ori)
          Mat_ori = bohr * Transpose( Mat_ori )
        endif

! Inverse of rotation matrix
        Angle = Angle_or(ipl) * pi / 180
        Mat(1,1) = 1._db; Mat(1,2:3) = 0._db
        Mat(2,1) = 0._db; Mat(2,2) = cos( Angle ); Mat(2,3) = sin( Angle )
        Mat(3,1) = 0._db; Mat(3,1) = - Mat(2,3);    Mat(3,3) = Mat(2,2)

        Mat = Matmul( Mat_ori, Mat )
        Mat = Matmul( Orthmatt, Mat )
! Polarization is along Ox
        v(:) = Mat(:,1)
! Wave vector is along Oz
        we(:) =  Mat(:,3)

        ws(:) = qk(:) - we(:)
        wsn = sqrt( sum( ws(:)**2 ) )
        ws(:) = ws(:) / wsn

        call prodvec(vsig,qk,we)
        wsn = sqrt( sum( vsig(:)**2 ) )
        vsig(:) = vsig(:) / wsn

        call prodvec(vpis,vsig,ws)

        Poldafse(:,ipl,1) = cmplx( v(:), 0._db )
        Vecdafse(:,ipl,1) = we(:)
        Poldafss(:,ipl,1) = cmplx( vsig(:), 0._db )
        Vecdafss(:,ipl,1) = ws(:)

        Poldafse(:,ipl,2) = cmplx( v(:), 0._db )
        Vecdafse(:,ipl,2) = we(:)
        Poldafss(:,ipl,2) = cmplx( vpis(:), 0._db )
        Vecdafss(:,ipl,2) = ws(:)

! Reflexion (-Q)
        Poldafse(:,ipl,3) = cmplx( v(:), 0._db )
        Vecdafse(:,ipl,3) = - we(:)
        Poldafss(:,ipl,3) = cmplx( vsig(:), 0._db )
        Vecdafss(:,ipl,3) = - ws(:)

        Poldafse(:,ipl,4) = cmplx( v(:), 0._db )
        Vecdafse(:,ipl,4) = - we(:)
        Poldafss(:,ipl,4) = cmplx( - vpis(:), 0._db )
        Vecdafss(:,ipl,4) = - ws(:)

      else

        dp = 2 * pi / nphi_dafs(ipl)

        if( angpoldafs(1,ipl) > -9999._db ) then
          cos_pe = cos( angpoldafs(1,ipl) )
          sin_pe = sin( angpoldafs(1,ipl) )
        endif
        if( angpoldafs(2,ipl) > -9999._db ) then
          cos_ps = cos( angpoldafs(2,ipl) )
          sin_ps = sin( angpoldafs(2,ipl) )
        endif

        do ip = 1,nphi_dafs(ipl)

          if( angpoldafs(1,ipl) < -9999._db ) then
            Angle = ( ip - 1 ) * dp
            cos_pe = cos( angle )
            sin_pe = sin( angle )
          elseif( angpoldafs(2,ipl) < -9999._db ) then
            Angle = ( ip - 1 ) * dp
            cos_ps = cos( angle )
            sin_ps = sin( angle )
          endif
          if( angpoldafs(3,ipl) < -9999._db ) then
            psi = ( ip - 1 ) * dp
          else
            psi = angpoldafs(3,ipl)
          endif

          if( ip == 1 .or. angpoldafs(3,ipl) < -9999._db ) then
            sinp = sin( psi )
            cosp = cos( psi )
            v(1) = cosb * cosp
            v(2) = cosb * sinp
            v(3) = - sinb
            we(:) = v(1) * qi(:) - v(2) * qj(:) + v(3) * qk(:)
            ws(:) = v(1) * qi(:) - v(2) * qj(:) - v(3) * qk(:)
            vsig(:) = sinp * qi(:) + cosp * qj(:)
            call prodvec(vpie,we,vsig)
            call prodvec(vpis,ws,vsig)
          endif

          vpe(:) = cos_pe * vsig(:) + sin_pe * vpie(:)
          vps(:) = cos_ps * vsig(:) + sin_ps * vpis(:)

! In diffraction all is defined with the complex conjugate, it is the same for circular polarizations
          select case( isigpi(1,ipl) )
            case(3)
              Poldafse(:,ipl,ip) = cmplx( vsig(:),-vpie(:), db) / sqrt(2._db)
            case(4)
              Poldafse(:,ipl,ip) = cmplx( vsig(:), vpie(:), db) / sqrt(2._db)
            case default
              Poldafse(:,ipl,ip) = cmplx( vpe(:), 0._db, db)
          end select
          Vecdafse(:,ipl,ip) = we(:)

          select case( isigpi(2,ipl) )
            case(3)
              Poldafss(:,ipl,ip) = cmplx( vsig(:),-vpis(:), db) / sqrt(2._db)
            case(4)
              Poldafss(:,ipl,ip) = cmplx( vsig(:), vpis(:), db) / sqrt(2._db)
            case default
              Poldafss(:,ipl,ip) = cmplx( vps(:), 0._db, db)
          end select
          Vecdafss(:,ipl,ip) = ws(:)

        end do

      endif
    endif

    if( sum( abs(hkl_dafs(:,ipl) ) ) < eps10 ) then

      if( Bormann ) then

        Q_mod(ipl) = 0._db

      else

        v(:) = Vecdafse(:,ipl,1)
        w(:) = Vecdafss(:,ipl,1)
        if( mpirank0 == 0 .and. ( sum(abs(v(:))) < eps10 .or. sum(abs(w(:))) < eps10 ) ) then
          call write_error
          do ipr = 6,9,3
            write(ipr,160)
          end do
          stop
        endif

! ThetaBragg is here the scattering angle
        pp = dot_product(v,w)
        if( abs( pp - 1._db ) < eps10 ) then
          Thetabragg = 0._db
        else
          Thetabragg = 0.5_db * acos( pp )
        endif
        Q_mod(ipl) = 2 * sin( thetabragg ) * konde
      endif

    else

      Q_mod(ipl) = 2 * pi / dhkl

    endif

    if( Film ) then
      if( ThetaBragg < eps10 ) then
        Length_abs(ipl) = 1.e+20_db
      else
        Length_abs(ipl) = c_cos_z * 2 / sin( thetabragg )
      endif
    endif

    if( icheck > 0 ) then
      if( ipl == 1 ) then
        if( Film ) then
          write(3,170)
        else
          write(3,175)
        endif
      endif
      if( Film .and. mod(ipl,10) /= 1 .and. npldafs > 20 .and. icheck == 1 ) cycle
      rad_i = 180._db / pi
      do i = 1,3
        if( abs( angpoldafs(i,ipl) ) > 9999._db ) then
          Ang(i) = 0._db
        else
          Ang(i) = angpoldafs(i,ipl)
        endif
      end do
      pe(:) = Poldafse(:,ipl,1)
      ps(:) = Poldafss(:,ipl,1)
      we(:) = Vecdafse(:,ipl,1)
      ws(:) = Vecdafss(:,ipl,1)
      pe = matmul( orthmati, pe )
      ps = matmul( orthmati, ps )
      we = matmul( orthmati, we )
      ws = matmul( orthmati, ws )
      if( Film ) then
        write(3,180) hkl_dafs(:,ipl), Thetabragg * rad_i, dhkl * bohr, Ang(1:3)*rad_i, ( Real(pe(i)), aimag(pe(i)), &
                 we(i), Real(ps(i)), aimag(ps(i)), ws(i), i = 1,3)
      else
        write(3,181) nint(hkl_dafs(:,ipl)), Thetabragg * rad_i, dhkl * bohr, Ang(1:3)*rad_i, ( Real(pe(i)), aimag(pe(i)), &
                 we(i), Real(ps(i)), aimag(ps(i)), ws(i), i = 1,3)
      endif
    endif

  end do  ! End of loop over ipl

  if( istop == 1 ) stop

  return
  120 format(//' The reflection number',i3,' : (h,k,l) = (',3f8.3,') does not exist at this photon energy !'//)
  125 format(//' The reflection number',i3,' : (h,k,l) = (',3i3,') does not exist at this photon energy !'//)
  130 format(/' (h,k,l) = (',3f8.3,')')
  135 format(/' (h,k,l) = (',3i3,')')
  140 format('  Reciprocal mesh base (A-1) :',/5x,'X',8x,'Y',8x,'Z', 3(/3f9.5))
  150 format('  Local base (A-1) :',/5x,'I',8x,'J',8x,'Q', 3(/3f9.5))
  160 format(' Calculations with (h,k,l) = (0,0,0)  need non zero values for the incoming and outgoing'/, &
        ' wave vector in order to calculate the scattering angle !')
  170 format(/'  (     h,     k,     l) ThetaBragg  d_hkl (A)    Polarization angles        Poldafse',5x, &
      'Vecdafse',7x,'Poldafss',5x,'Vecdafss'/)
  175 format(/'  (h,k,l) ThetaBragg  d_hkl (A)    Polarization angles        Poldafse',5x, &
      'Vecdafse',7x,'Poldafss',5x,'Vecdafss'/)
  180 format(3f8.3,2f10.3,2x,3f8.3,2(2x,2f8.5,2x,f8.5)/, 2(70x,2(2x,2f8.5,2x,f8.5)/))
  181 format(3i3,2f10.3,2x,3f8.3,2(2x,2f8.5,2x,f8.5)/, 2(55x,2(2x,2f8.5,2x,f8.5)/))
  end

!***********************************************************************

! Roughness model with erfc function

function Roughness_erfc(z,z0,Width)

  use declarations
  implicit none

  real(kind=db):: Roughness_erfc, Width, z, z0
  
  Roughness_erfc = 0.5_db * erfc( ( z - z0 ) / ( sqrt( 2._db ) * Width ) )  

  return
end

!***********************************************************************

! Atomic structure factor

subroutine Atom_structure_factor(Energy_photon, f_ms, f_mo, f_no_res, f0, itypepr, lvval, Magnetic, n_atom_proto, nlat, nlatm, &
                                 npldafs, nrato, nrm, nspin, ntype, numat, Q_mod, rato, popatm, psival)

  use declarations
  implicit none

! mc2_Rydb = mc2 in Rydberg = m c^2 / ( e Rydb ) = 2 / alfa^2
  real(kind=db), parameter:: mc2_Rydb = 2 / alfa_sf**2   ! = 37557.7301195846

  integer:: ipl, ipr, it, n_atom_proto, nlatm, npldafs, nrm, nspin, ntype
  integer, dimension(0:n_atom_proto):: itypepr
  integer, dimension(0:ntype):: nlat, nrato, numat
  integer, dimension(0:ntype,nlatm):: lvval

  character(len=4):: elemv
  character(len=2):: Chemical_symbol

  logical:: Magnetic

  real(kind=sg):: getf0, s
  real(kind=db):: Energy_photon, fmo, fms, x
  real(kind=db), dimension(2):: f_no_res
  real(kind=db), dimension(npldafs):: Q_mod
  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(n_atom_proto,npldafs):: f_ms, f_mo, f0
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm

  do ipl = 1,npldafs

    x = Q_mod(ipl) / ( 4 * pi )      ! = sin(theta) / lambda
    s = x / bohr                 ! one needs it in angstroem^(-1)

    do ipr = 1,n_atom_proto

      it = itypepr( ipr )
      elemv = ' '
      elemv(1:2) = Chemical_Symbol( numat(it) )
      f0(ipr,ipl) = getf0(elemv,s)

! Non-resonant scattering amplitude
! orbital moment is taken parralel to spin axis...
      if( Magnetic ) then
        call get_fmag(Q_mod(ipl),fmo,fms,ipr,it,lvval,n_atom_proto,nlat(it),nlatm,nrato(it),nrm,nspin,ntype,popatm,psival,rato)
        f_ms(ipr,ipl) = - ( Energy_photon / mc2_Rydb ) * f_no_res(1) * fms
        if( f_no_res(2) < -99._db ) then
          f_mo(ipr,ipl) = - ( Energy_photon / mc2_Rydb ) * f_no_res(1) * fmo
        else
! Factor 2 because f_ms corresponds to S et not to 2S.
          f_mo(ipr,ipl) = 2 * f_no_res(2) * f_ms(ipr,ipl)
        endif
      endif
    end do

  end do

  return
end

!***********************************************************************

! Basis change RR Bulk --> RR film.

subroutine Bulk_base_tr(angxyz,angxyz_bulk,axyz,axyz_bulk,Film_shift,icheck,Mat_bulk,Mult_bulk,Mult_film)

  use declarations
  implicit none

  integer:: i, icheck
  integer, dimension(2):: Mult_bulk, Mult_film

  real(kind=db):: Cal_Volume_maille, cos_g, sin_g, Vol
  real(kind=db), dimension(3):: angxyz, angxyz_bulk, axyz, axyz_bulk, cos_a, V, Vec_a, Vec_b, Vec_c
  real(kind=db), dimension(4):: Film_shift
  real(kind=db), dimension(3,3):: Mat, Mat_bulk, Mat_film, Mat_i

! Bulk

  cos_a(:) = cos( angxyz_bulk(:) )
  sin_g = sin( angxyz_bulk(3) )

  Vec_a(1) = axyz_bulk(1); Vec_a(2:3) = 0._db
  Vec_b(1) = axyz_bulk(2) * cos_a(3); Vec_b(2) = axyz_bulk(2) * sin_g; Vec_b(3) = 0._db
  Vec_c(1) = axyz_bulk(3) * cos_a(2); Vec_c(2) = axyz_bulk(3) * ( cos_a(1) - cos_a(3)*cos_a(2) ) / sin_g
  Vec_c(3) = sqrt( axyz_bulk(3)**2 - Vec_c(1)**2 - Vec_c(2)**2 )

  Vol = Cal_Volume_maille(axyz_bulk,angxyz_bulk)

! base RR bulk --> Matrice base ortho
  call prodvec(V,Vec_b,Vec_c)
  Mat_bulk(:,1) = V(:) / Vol
  call prodvec(V,Vec_c,Vec_a)
  Mat_bulk(:,2) = V(:) / Vol
  call prodvec(V,Vec_a,Vec_b)
  Mat_bulk(:,3) = V(:) / Vol

  if( abs( Film_shift(4) ) > eps10 ) then
    cos_g = cos( Film_shift(4) )
    sin_g = sin( Film_shift(4) )
    Mat(1,1) = cos_g; Mat(1,2) = - sin_g; Mat(1,3) = 0._db
    Mat(2,1) = sin_g; Mat(2,2) = cos_g;   Mat(2,3) = 0._db
    Mat(3,1) = 0._db; Mat(3,2) = 0._db;   Mat(3,3) = 1._db
    Mat_bulk = matmul( Mat, Mat_bulk )
  endif

! Film

  cos_a(:) = cos( angxyz(:) )
  sin_g = sin( angxyz(3) )

  Vec_a(1) = axyz(1); Vec_a(2:3) = 0._db
  Vec_b(1) = axyz(2) * cos_a(3); Vec_b(2) = axyz(2) * sin_g; Vec_b(3) = 0._db
  Vec_c(1) = axyz(3) * cos_a(2); Vec_c(2) = axyz(3) * ( cos_a(1) - cos_a(3)*cos_a(2) ) / sin_g
  Vec_c(3) = sqrt( axyz(3)**2 - Vec_c(1)**2 - Vec_c(2)**2 )

  Vol = Cal_Volume_maille(axyz,angxyz)

! Matrix RR --> Orthonormal basis
  call prodvec(V,Vec_b,Vec_c)
  Mat_film(:,1) = V(:) / Vol
  call prodvec(V,Vec_c,Vec_a)
  Mat_film(:,2) = V(:) / Vol
  call prodvec(V,Vec_a,Vec_b)
  Mat_film(:,3) = V(:) / Vol

  call invermat(Mat_film,Mat_i)

  Mat_bulk = matmul( Mat_i, Mat_bulk )

  if( icheck > 0 ) then
    write(3,'(/A)') ' Bulk - Film reciprocal space transformation Matrix'
    do i = 1,3
      write(3,'(3f10.5)') Mat_bulk(i,:)
    end do
  endif

! PPCM
  call mult_film_bulk(angxyz_bulk,axyz,axyz_bulk,Film_shift,Mult_bulk,Mult_film)

  if( icheck > 0 ) then
    write(3,*)
    write(3,110) 'film', Mult_film(:)
    write(3,110) 'bulk', Mult_bulk(:)
  endif

  return
  110 format(' Mult_',a4,' =',2i5)
end

!***********************************************************************

  subroutine mult_film_bulk(angxyz_bulk,axyz,axyz_bulk,Film_shift,Mult_bulk,Mult_film)

  use declarations
  implicit none

  integer:: i, j, k
  integer, dimension(2):: Mult_bulk, Mult_film

  logical:: Diagonal

  real(kind=db):: diag_bulk, Val_bulk, Val_film
  real(kind=db), dimension(3):: angxyz_bulk, axyz, axyz_bulk
  real(kind=db), dimension(4):: Film_shift

  Diagonal = abs( angxyz_bulk(3) - pi / 2 ) < eps10 .and. &
            ( abs( Film_shift(4) - pi / 4 ) <  eps10 .or. abs( Film_shift(4) + pi / 4 ) < eps10 )

  if( Diagonal ) diag_bulk = sqrt( axyz_bulk(1)**2 + axyz_bulk(2)**2 )

  do i = 1,2
    boucle_j: do j = 1,1000
      Val_film = j * axyz(i)
      do k = 1,1000
        if( Diagonal ) then
          Val_bulk = k * diag_bulk
        else
          Val_bulk = k * axyz_bulk(i)
        endif
        if( abs( Val_film - Val_bulk ) < eps6 ) exit boucle_j
        if( Val_bulk > Val_film ) exit
      end do
    end do boucle_j
    if( j > 1000 .or. k > 1000 ) then
      Mult_film(i) = 0
      Mult_bulk(i) = 0
    else
      Mult_film(i) = j
      Mult_bulk(i) = k
    endif
  end do

  return
end

!*********************************************************************************************************

! Calculation of bulk non resonant scattering amplitude
! When there is an absorbing atom in the bulk, the calculation is only for the non resonant scattering dure to the roughness 

subroutine Bulk_scat(Abs_in_bulk,angxyz_bulk,axyz_bulk,Bulk_roughness,c_cos_z_b,Delta_bulk,delta_z_top_bulk,Energy_photon, &
        hkl_dafs,hkl_film,icheck,Length_abs,Mat_bulk_i,Mult_bulk,n_atom_bulk,ngroup_taux,ngroup_temp,npldafs,phd_f_bulk, &
        posn_bulk,Q_mod,Taux_oc,Temp_coef,Temp_B_iso,Truncation,Z_abs,Z_bulk)

  use declarations
  implicit none

  integer:: i, i_sens, icheck, igr, ipl, j, jgr, n_atom_bulk, ngroup_taux, ngroup_temp, npldafs, Z, Z_abs

  integer, dimension(2):: Mult_bulk
  integer, dimension(n_atom_bulk):: Z_bulk

  character(len=4):: elemv
  character(len=2):: Chemical_symbol

  logical:: Abs_in_bulk, hkl_film, Temp_B_iso

  complex(kind=db):: Bragg_bulk, cfac
  complex(kind=db), dimension(npldafs):: phd_f_bulk, Truncation

  real(kind=sg):: getf0, s

  real(kind=db):: arg, Bulk_roughness, c_cos_z_b, Cal_Volume_maille, conv_mbarn_nelec, Delta_2, Delta_bulk, delta_z_top_bulk, &
    Energy_photon, fpp_bulk_tot, Orig_top_bulk, Roughness_erfc, Taux, Volume_maille, x, z_pos, z_top

  real(kind=db), dimension(3):: angxyz_bulk, axyz_bulk, hkl, p
  real(kind=db), dimension(3,3):: Mat_bulk_i
  real(kind=db), dimension(n_atom_bulk):: fp_bulk, fpp_bulk
  real(kind=db), dimension(n_atom_bulk,npldafs):: Abs_temp, f0_bulk
  real(kind=db), dimension(ngroup_taux):: Taux_oc
  real(kind=db), dimension(ngroup_temp):: Temp_coef
  real(kind=db), dimension(3,n_atom_bulk):: posn_bulk
  real(kind=db), dimension(npldafs):: abs_mesh, Length_abs, Q_mod
  real(kind=db), dimension(3,npldafs):: hkl_dafs

  Orig_top_bulk = delta_z_top_bulk - Delta_Bulk
  
  fp_bulk(:) = 0._db
  fpp_bulk(:) = 0._db
  fpp_bulk_tot = 0._db
  do igr = 1,n_atom_bulk
    Z = Z_bulk(igr)
    if( Z == Z_abs ) cycle
    call fprime(Z,Energy_photon,fpp_bulk(igr),fp_bulk(igr))
    fpp_bulk_tot = fpp_bulk_tot + fpp_bulk(igr)
  end do
! Conversion in Megabarn (= 10^-18 cm2 = 10^-22 m2 = 10^-2 A2)
  fpp_bulk_tot = fpp_bulk_tot / conv_mbarn_nelec(Energy_photon)
! Conversion in linear absorption coefficient in micrometer^-1
  Volume_maille = Cal_Volume_maille(axyz_bulk,angxyz_bulk)
  fpp_bulk_tot = 100 * fpp_bulk_tot / ( Volume_maille * bohr**3 )

  Truncation(:) = ( 0._db, 0._db )
  phd_f_bulk(:) = ( 0._db, 0._db )

  z_top = -100000._db
  do igr = 1,n_atom_bulk
    z_top = max( posn_bulk(3,igr), z_top )
  end do

! Calculation of Thomson scattering
  do ipl = 1,npldafs
    x = Q_mod(ipl) / ( 4 * pi )      ! = sin(theta) / lambda
    s = x / bohr                 ! on le veut en angstroem - 1
    do igr = 1,n_atom_bulk
      elemv = ' '
      if( igr > n_atom_bulk ) then
        jgr = igr - n_atom_bulk
        elemv(1:2) = Chemical_Symbol( Z_bulk(jgr) )
      else
        elemv(1:2) = Chemical_Symbol( Z_bulk(igr) )
      endif
      f0_bulk(igr,ipl) = getf0(elemv,s)
    end do
  end do

  if( Temp_B_iso ) then
    do ipl = 1,npldafs
      Delta_2 = ( Q_mod(ipl) / ( quatre_pi * bohr ) )**2
      do igr = 1,n_atom_bulk
        Abs_temp(igr,ipl) = exp( - Temp_coef(ngroup_temp-n_atom_bulk+igr) * Delta_2 )
      end do
    end do
  endif

  do ipl = 1,npldafs

    hkl(:) = hkl_dafs(:,ipl)
    if( hkl_film ) hkl = Matmul( Mat_bulk_i, hkl )
! Case of irrationnal unit cells
    if( ( abs( hkl(1) ) > eps10 .and. Mult_bulk(1) == 0 ) .or. ( abs( hkl(2) ) > eps10 .and. Mult_bulk(2) == 0 ) ) cycle

    if( .not. Abs_in_bulk ) then
    
      abs_mesh(ipl) = exp( - fpp_bulk_tot * Length_abs(ipl) )

      Truncation(ipl) = 1 / ( 1._db - abs_mesh(ipl) * exp( - img * deux_pi * hkl(3) ) )

      do igr = 1,n_atom_bulk

        arg = deux_pi * sum( posn_bulk(:,igr) * hkl(:) )
        Bragg_bulk = cmplx( cos(arg), sin(arg), db )

! Fraction of absorption (avoid the singularity when l close to 1)
        Bragg_bulk = Bragg_bulk * exp( - ( z_top - posn_bulk(3,igr) ) * fpp_bulk_tot * Length_abs(ipl) )

        if( Temp_B_iso ) Bragg_bulk = Bragg_bulk * Abs_temp(igr,ipl)
        if( ngroup_taux > 0 ) Bragg_bulk = Bragg_bulk * Taux_oc(ngroup_taux-n_atom_bulk+igr)

        phd_f_bulk(ipl) = phd_f_bulk(ipl) + Bragg_bulk * cmplx( f0_bulk(igr,ipl) + fp_bulk(igr), fpp_bulk(igr), db )

      end do

      phd_f_bulk(ipl) = phd_f_bulk(ipl) * Truncation(ipl)

    endif
    
! Taking into account of bulk roughness
! One substracts the part inside the semi-ininite bulk and add the part outwards
    if( Bulk_roughness > eps10 ) then

      do igr = 1,n_atom_bulk

        Bragg_bulk = (0._db, 0._db )

        do i_sens = -1,1,2

          p(:) = posn_bulk(:,igr) 
        
          do i = 0,10000
  
            if( i_sens == -1 .and. i == 0 ) cycle
            p(3) = posn_bulk(3,igr) + i_sens * i

            z_pos = p(3) * c_cos_z_b
            if( abs( z_pos - Orig_top_bulk ) > 4 * Bulk_roughness ) then
              if( i == 0 ) then
                cycle
              else
                exit
              endif
            endif

            arg = deux_pi * sum( p(:) * hkl(:) )

            if( z_pos > Orig_top_bulk ) then 
              Taux = Roughness_erfc(z_pos,Orig_top_bulk,Bulk_roughness)
            else
        ! One substracts the missing weight
              Taux = - 1 + Roughness_erfc(z_pos,Orig_top_bulk,Bulk_roughness)
            endif
                        
            Bragg_bulk = Bragg_bulk + Taux * cmplx( cos(arg), sin(arg), db ) 

          end do
        end do

        if( Temp_B_iso ) Bragg_bulk = Bragg_bulk * Abs_temp(igr,ipl)
        if( ngroup_taux > 0 ) Bragg_bulk = Bragg_bulk * Taux_oc(ngroup_taux-n_atom_bulk+igr)

        phd_f_bulk(ipl) = phd_f_bulk(ipl) + Bragg_bulk * cmplx( f0_bulk(igr,ipl) + fp_bulk(igr), fpp_bulk(igr), db ) 
    
      end do

    endif

! Taking into account of the superstructure
    cfac = ( 0._db, 0._db )
    do i = 0,Max( Mult_bulk(1)-1, 0 )
  ! When film and bulk are not rationnal, non specular relection contains just the film (when hkl_film) 
      if( hkl_film .and. Mult_bulk(1) == 0 .and. abs( hkl(1) ) > eps10 ) cycle
      do j = 0,Max( Mult_bulk(2)-1, 0 )
        if( hkl_film .and. Mult_bulk(2) == 0 .and. abs( hkl(2) ) > eps10 ) cycle
        arg = deux_pi * ( i * hkl(1) + j * hkl(2) )
        cfac = cfac + cmplx( cos(arg), sin(arg), db )
      end do
    end do

    phd_f_bulk(ipl) = cfac * phd_f_bulk(ipl)

  end do

  if( icheck > 0 .and. .not. Abs_in_bulk ) then
    write(3,100) ' mu_0_bulk      =', fpp_bulk_tot
    write(3,110) '   Length_abs(1,11...) =', ( Length_abs(ipl), ipl = 1,min(npldafs,201), 10 )
    write(3,110) ' 1 - abs_mesh(1,11...) =', ( 1._db - abs_mesh(ipl), ipl = 1,min(npldafs,201), 10 )
  endif

  return
  100 format(a17,1p,e13.5,' micrometer^-1')
  110 format(a24,1p,120e13.5)
end

!***********************************************************************************************************

! Phase between bulk and film.
! Contain also the ratio between the unit cell surfaces

  subroutine Cal_Phase_Bulk(angxyz,angxyz_bulk,axyz,axyz_bulk,c_cos_z,c_cos_z_b,Delta_bulk,Film_shift,hkl_dafs,hkl_film, &
                          Mult_bulk,npldafs,Phase_bulk)

  use declarations
  implicit none

  integer:: i, ipl, j, npldafs

  integer, dimension(2):: Mult_bulk

  complex(kind=db), dimension(npldafs):: Phase_Bulk

  logical:: hkl_film

  real(kind=db):: arg, c_cos_z, c_cos_z_b, Delta_bulk, Rap
  real(kind=db), dimension(3):: angxyz, angxyz_bulk, axyz, axyz_bulk, hkl
  real(kind=db), dimension(4):: Film_shift
  real(kind=db), dimension(3,npldafs):: hkl_dafs

  Phase_Bulk(:) = ( 0._db, 0._db )

  do ipl = 1,npldafs

! No need to transform hkl if hkl_film true, because below division by c_cos_z or c_ cos_z_b
    hkl(:) = hkl_dafs(:,ipl)

! Mult_bulk = 0 is for non rationnal ratio between bulk and film
    do i = 0,Max( Mult_bulk(1)-1, 0 )
      do j = 0,Max( Mult_bulk(2)-1, 0 )
        arg = deux_pi * ( i * hkl(1) + j * hkl(2) )
        Phase_Bulk(ipl) = Phase_Bulk(ipl) + cmplx( cos(arg), sin(arg), db )
      end do
    end do
    if( abs( Phase_Bulk(ipl) ) < eps10 ) Phase_Bulk(ipl) = ( 0._db, 0._db )

    if( hkl_film ) then
      arg = deux_pi * hkl(3) * Delta_bulk / c_cos_z
    else
      arg = deux_pi * hkl(3) * Delta_bulk / c_cos_z_b
    endif
    Phase_Bulk(ipl) = Phase_Bulk(ipl) * cmplx( cos(arg), sin(arg), db )

! Case of non rationnal ratio between bulk and film, or R45 (sqrt(2), sqrt(2))
    Rap = 1._db
    if( Mult_bulk(1) == 0 .and. Mult_bulk(2) == 0 ) then
      Rap = axyz(1) * axyz(2) * sin( angxyz(3) ) / ( axyz_bulk(1) * axyz_bulk(2) * sin( angxyz_bulk(3) ) )
    else
      do i = 1,2
        if( ( abs( hkl(i) ) < eps10 .and. Mult_bulk(i) == 0 ) .or. abs( Film_shift(4) - pi / 4 ) < eps10 ) &
          Rap = Rap * axyz(i) / axyz_bulk(i)
      end do
    endif

    Phase_Bulk(ipl) = Rap * Phase_Bulk(ipl)

  end do

  return
end

!***********************************************************************************************************

! Calculation of the scattering amplitude from the cap layer

subroutine Cap_scat(angxyz,angxyz_bulk,angxyz_cap,axyz,axyz_bulk,axyz_cap,Cap_B_iso,Cap_roughness, &
              Cap_thickness,c_cos_z,c_cos_z_b,Delta_cap,Delta_roughness_film,delta_z_bottom_cap, &
              delta_z_top_cap,Energy_photon,Film_roughness,fpp_cap_tot,hkl_dafs,hkl_film, &
              Mult_bulk,Mult_film,n_atom_cap,npldafs,phd_f_cap,posn_cap,Q_mod,Taux_cap, &
              Z_cap)

  use declarations
  implicit none

  integer:: i, i0, igr, ipl, n_atom_cap, n, npldafs, Z

  integer, dimension(n_atom_cap):: Z_cap

  integer, dimension(2):: Mult_bulk, Mult_film

  character(len=4):: elemv
  character(len=2):: Chemical_symbol

  complex(kind=db):: Bragg_cap
  complex(kind=db), dimension(npldafs):: phd_f_cap

  logical:: hkl_film

  real(kind=sg):: getf0, s

  real(kind=db):: arg, Conv_mbarn_nelec, c_cos_z, c_cos_z_b, c_cos_z_c, cos_z_c, Cap_amount, Cap_amount_ref, &
    Cap_B_iso, Cap_roughness, Cap_thickness, Cal_Volume_maille, Deb, Delta_cap, &
    Delta_roughness_film, delta_z_bottom_cap, delta_z_top_cap, Energy_photon, f0_cap,  &
    Film_roughness, fpp_cap_tot, Orig_bottom_cap, Orig_top_cap, Roughness_erfc, Taux_r, Taux_r2, Volume_maille, x, z_pos

  real(kind=db), dimension(3):: angxyz, angxyz_bulk, angxyz_cap, axyz, axyz_bulk, axyz_cap, hkl, p
  real(kind=db), dimension(n_atom_cap):: fp_cap, fpp_cap, Taux_cap
  real(kind=db), dimension(npldafs):: Q_mod
  real(kind=db), dimension(3,n_atom_cap):: posn_cap
  real(kind=db), dimension(3,npldafs):: hkl_dafs

  Orig_top_cap = Cap_thickness + delta_z_top_cap
  Orig_bottom_cap = - delta_z_bottom_cap

  cos_z_c = sqrt( sin( angxyz_cap(2) )**2 - ( ( cos( angxyz_cap(1) ) - cos( angxyz_cap(3) ) * cos( angxyz_cap(2) ) ) &
          / sin( angxyz_cap(3) ) )**2 )

  c_cos_z_c = cos_z_c * axyz_cap(3)

  if( Delta_roughness_film > eps10 ) then
    i0 = - nint( Delta_roughness_film / c_cos_z_c ) - 1
  else
    i0 = 0
  endif

  fpp_cap_tot = 0._db
  do igr = 1,n_atom_cap
    Z = Z_cap(igr)
    call fprime(Z,Energy_photon,fpp_cap(igr),fp_cap(igr))
    fpp_cap_tot = fpp_cap_tot + Taux_cap(igr) * fpp_cap(igr)
  end do

! Absorption through the layer
! Conversion in Megabarn (= 10^-18 cm2 = 10^-22 m2 = 10^-2 A2)
  fpp_cap_tot = fpp_cap_tot / conv_mbarn_nelec(Energy_photon)
! Conversion in linear absorption coefficient in micrometre^-1
  Volume_maille = Cal_Volume_maille(axyz_cap,angxyz_cap)
  fpp_cap_tot = 100 * fpp_cap_tot / ( Volume_maille * bohr**3 )
! Thickness to cross (go and return) (to multiply by 2/sin(theta_bragg) )
  fpp_cap_tot = 0.00001_db * bohr * fpp_cap_tot * Cap_thickness

  do ipl = 1,npldafs

    hkl(:) = hkl_dafs(:,ipl)

    if( abs( hkl(1) ) > eps10 .or. abs( hkl(2) ) > eps10 ) cycle

    x = Q_mod(ipl) / ( 4 * pi )      ! = sin(theta) / lambda
    s = x / bohr                      ! one want it in angstroem - 1

! ( Q_mod(ipl) / ( quatre_pi * bohr ) )**2 = ( Sin(Theta_Bragg)/Lambda )^2   (lambda in Angstroem)
! Cap_B_iso = 8*pi^2 * <u^2> is in Angtroem^2
    if( Cap_B_iso > eps10 ) then
      Deb = exp( - Cap_B_iso * ( Q_mod(ipl) / ( quatre_pi * bohr ) )**2 )
    else
      Deb = 1._db
    endif

    if( hkl_film ) then
      hkl(3) = hkl(3) * c_cos_z_c / c_cos_z
    else
      hkl(3) = hkl(3) * c_cos_z_c / c_cos_z_b
    endif

    Z = 0

    Cap_amount = 0._db
    Cap_amount_ref = 0._db

    do igr = 1,n_atom_cap

      elemv = ' '
      elemv(1:2) = Chemical_Symbol( Z_cap(igr) )
      f0_cap = getf0(elemv,s)

      Bragg_cap = ( 0._db, 0._db )

      do i = i0,1000000
        p(:) = posn_cap(:,igr)
        p(3) = p(3) + i
        z_pos = p(3) * c_cos_z_c
        if( z_pos < - delta_z_bottom_cap - Delta_roughness_film - eps10 ) cycle
        if( Cap_roughness > eps10 ) then
        ! no absolute value 
          if( z_pos - Orig_top_cap > 4 * Cap_roughness + eps10 ) exit
          Taux_r = Roughness_erfc(z_pos,Orig_top_cap,Cap_roughness)
        else
          if( z_pos > Cap_thickness + eps10 ) exit
          Taux_r = 1._db
        endif
        if( Film_roughness > eps10 .and. z_pos < Delta_roughness_film + eps10 ) then
          Taux_r2 = 1 - Roughness_erfc(z_pos,Orig_bottom_cap,Film_roughness)
          Taux_r = Taux_r * Taux_r2
        endif
        Cap_amount = Cap_amount + Taux_cap(igr) * Taux_r
        if( z_pos < Orig_top_cap .and. z_pos > Orig_bottom_cap ) Cap_amount_ref = Cap_amount_ref + Taux_cap(igr)
        arg = deux_pi * sum( p(:) * hkl(:) )
        Bragg_cap = Bragg_cap + Deb * Taux_r * cmplx( cos(arg), sin(arg), db )
      end do

      phd_f_cap(ipl) = phd_f_cap(ipl) + Taux_cap(igr) * Bragg_cap * cmplx( f0_cap + fp_cap(igr), fpp_cap(igr), db )

    end do

! This renormaization is to keep the same amount of cap atoms whatever is the roughness
    if( Cap_amount > eps10 ) phd_f_cap(ipl) = ( Cap_amount_ref / Cap_amount ) * phd_f_cap(ipl)

! hkl(3) is given in cap cell unit
    arg = deux_pi * hkl(3) * Delta_cap / c_cos_z_c
    phd_f_cap(ipl) = phd_f_cap(ipl) * cmplx( cos(arg), sin(arg), db )

! Case of non rationnal ratio between bulk and film, or R45 (sqrt(2), sqrt(2))

    if( hkl_film ) then
      n = max( 1, Mult_film(1) * Mult_film(2) )
      x = n * axyz(1) * axyz(2) * sin( angxyz(3) ) / ( axyz_cap(1) *  axyz_cap(2) * sin( angxyz_cap(3) ) )
    else
      n = max( 1, Mult_bulk(1) * Mult_bulk(2) )
      x = n * axyz_bulk(1) * axyz_bulk(2) * sin( angxyz_bulk(3) ) &
        / ( axyz_cap(1) *  axyz_cap(2) * sin( angxyz_cap(3) ) )
    endif
    phd_f_cap(ipl) = x * phd_f_cap(ipl)

  end do

  return
end

!***********************************************************************

! Establishment of the column names for dafs reflexions

subroutine Col_dafs_name(Angpoldafs,Bormann,Full_self_abs,hkl_dafs,isigpi,mpirank0,ncolm,ncolr,ncolt, &
         nomabs,npldafs,Operation_mode_used,Self_abs)

  use declarations
  implicit none

  integer:: i, icol, index, ipldafs, ipr, j, k, mpirank0, ncolm, ncolr, ncolt, npldafs
  integer, dimension(2,npldafs):: isigpi

  character(len=Length_word) nomab
  character(len=Length_word), dimension(ncolm):: nomabs

  logical:: Bormann, Full_self_abs, Operation_mode_used, Self_abs

  real(kind=db):: a
  real(kind=db), dimension(3):: hkl
  real(kind=db), dimension(3,npldafs):: angpoldafs, hkl_dafs

  icol = ncolr

  do ipldafs = 1,npldafs
    icol = icol + 1
    hkl(:) = hkl_dafs(1:3,ipldafs)
    nomab = ' r('
    j = 3
    call trnom_r(j,hkl,nomab)
    j = j + 1
    if( j < Length_word-1 ) nomab(j:j) = ')'

    if( .not. Operation_mode_used .or. abs( angpoldafs(3,ipldafs) ) < 9999._db ) then
      if( j < Length_word-1 ) then
        if( abs( angpoldafs(3,ipldafs) ) < 9999._db ) then
          a = angpoldafs(3,ipldafs) * 180 / pi
        else
          a = 0._db
        endif
        index = nint( a )
        call ad_number(index,nomab,Length_word)
      endif
      j = len_trim(nomab)
    endif

    if( Bormann ) then

      j = j + 1
      if( j > Length_word ) exit
      select case( ipldafs )
        case(1,2,3,10,11,12,19,20,21,28,29,30)
          nomab(j:j) = 'x'
        case(4,5,6,13,14,15,22,23,24,31,32,33)
          nomab(j:j) = 'y'
        case(7,8,9,16,17,18,25,26,27,34,35,36)
          nomab(j:j) = 'z'
      end select
      j = j + 1
      if( j > Length_word ) exit
      select case( ipldafs )
        case(1,4,7,10,13,16,19,22,25,28,31,34)
          nomab(j:j) = 'x'
        case(2,5,8,11,14,17,20,23,26,29,32,35)
          nomab(j:j) = 'y'
        case(3,6,9,12,15,18,21,24,27,30,33,36)
          nomab(j:j) = 'z'
      end select

    else

      do i = 1,2
        j = j + 1
        if( i == 1 ) k = j
        if( j > Length_word ) exit
        select case( isigpi(i,ipldafs) )
          case(1)
            nomab(j:j) = 's'
          case(2)
            nomab(j:j) = 'p'
          case(3)
            nomab(j:j) = 'r'
          case(4)
            nomab(j:j) = 'l'
          case default
            nomab(j:j) = 'a'
        end select
      end do

    endif

    nomabs(icol) = nomab

    icol = icol + 1
    nomab(2:2) = 'i'
    nomabs(icol) = nomab

    if( Full_self_abs ) then
      nomab = ' '
      icol = icol + 1
      select case( mod(ipldafs,4) )
        case(1)
          nomab(2:11) = 'mu_ss_in_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ss_in_i'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ss_ou_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ss_ou_i'
          nomabs(icol) = nomab
        case(2)
          nomab(2:11) = 'mu_sp_in_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_sp_in_i'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_sp_ou_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_sp_ou_i'
          nomabs(icol) = nomab
        case(3)
          nomab(2:11) = 'mu_ps_in_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ps_in_i'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ps_ou_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_ps_ou_i'
          nomabs(icol) = nomab
        case(0)
          nomab(2:11) = 'mu_pp_in_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_pp_in_i'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_pp_ou_r'
          nomabs(icol) = nomab
          icol = icol + 1
          nomab(2:11) = 'mu_pp_ou_i'
          nomabs(icol) = nomab
      end select
    elseif( Self_abs ) then
      icol = icol + 1
      nomab(2:2) = 'A'
      if( k < Length_word ) nomab(k:k+1) = 'in'
      if( k == Length_word ) nomab(k:k) = 'i'
      nomabs(icol) = nomab
      icol = icol + 1
      if( k < Length_word ) nomab(k:k+1) = 'ou'
      if( k == Length_word ) nomab(k:k) = 'o'
      nomabs(icol) = nomab
    endif
  end do

! ncolt = total number of columns for both XANES and DAFS
  ncolt = icol
  if( ncolt > ncolm .and. mpirank0 == 0 ) then
    call write_error
    do ipr = 6,9,3
      write(ipr,105) ncolt, ncolm
    end do
    stop
  endif

  105 format(/' ncolt =',i3,' > ncolm = ',i3, /' Bug in the program !')

  return
end

!***********************************************************************

subroutine get_fmag(Q_mod,fmo,fms,ipr,it,lvval,n_atom_proto,nl,nlatm,nr,nrm,nspin,ntype,popatm,psival,rato)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  integer, dimension(0:ntype,nlatm):: lvval

  real(kind=db), dimension(0:nrm,0:ntype):: rato
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival

  fms = 0._db
  fmo = 0._db

  do io = 1,nl

    popt = popatm(ipr,io,1) + popatm(ipr,io,nspin)
    spin = 0.5_db * ( popatm(ipr,io,1) - popatm(ipr,io,nspin) )
    pop = abs( popatm(ipr,io,1) - popatm(ipr,io,nspin) )

    if( abs(spin) < eps10 ) cycle

    l = lvval(it,io)
    hund_mo = ( pop * l - pop * ( pop - 1 ) / 2 ) / (2 * abs(spin))
    if( popt < 2 * l + 1 ) hund_mo = - hund_mo

    f = 0._db
! psival is in fact sqrt(4*pi)*r*psi
    do ir = 2,nr-1
      dr = rato(ir+1,it) - rato(ir-1,it)
      f = f + ( psival(ir,io,it)**2 / rato(ir,it) ) * sin( Q_mod * rato(ir,it) ) * dr
    end do
    f = 0.5_db * f / Q_mod

    fms = fms + spin * f
    fmo = fmo + hund_mo * 2 * spin * f

  end do

! Empirical factor
  fmo = 0.2_db * fmo

  return
end

!*********************************************************************

subroutine get_vec_b(pe,ps,vec_b,we,ws)

  use declarations
  implicit real(kind=db) (a-h,o-z)

  complex(kind=db), dimension(3):: vec_b, kve, kvs, p, pe, ps, q, ve, vs

  real(kind=db), dimension(3):: we, ws

  ve(:) = cmplx( we(:), 0._db, db )
  vs(:) = cmplx( ws(:), 0._db, db )

  call prodvec_cp( kve, ve, pe )
  call prodvec_cp( kvs, vs, ps )
  call prodvec_cp( p, ps, pe )
  call prodvec_cp( q, kvs, kve )

  vec_b(:) = p(:) + dot_product(vs,pe) * kvs(:) - dot_product(ve,ps) * kve(:) - q(:)

  return
end

!*********************************************************************

! Calculation of A" in Blume and Gibbs, PRB 37,1779 (1988).
! Mofied because of the error in this paper.

subroutine get_vec_a(pe,ps,vec_a,we,ws)

  use declarations
  implicit none

  complex(kind=db):: vv
  complex(kind=db), dimension(3):: vec_a, p, pe, ps, ve, ves, vs

  real(kind=db), dimension(3):: we, ws

  ve(:) = cmplx( we(:), 0._db, db )
  vs(:) = cmplx( ws(:), 0._db, db )

  vv = dot_product(ve,vs)

  call prodvec_cp( p, ps, pe )

! Formule Blume
!      vec_a(:) = 2 * ( 1 - vv ) * p(:) - dot_product(ve,ps) * kve(:)
!     &         + dot_product(vs,pe) * kvs(:)

! Formule modifiee
! La premiere partie de la formule en haut de la page 1781 est bonne
! A" = A' - (A'.K)K, mais le suite est fausse.
!      vec_a(:) = - 2 * ( 1 - vv )
!     &         * ( p(:) + ( dot_product(kve,ps) + dot_product(kvs,pe) )
!     &                  * ( ve(:) - vs(:) ) )

  ves(:) = ve(:) - vs(:)

  vec_a(:) = - 2 * ( 1 - vv ) * p(:) - dot_product( ves, p ) * ves(:)

  return
end

!*********************************************************************

! Oana Bunau 15 april 2007
! q = vecteur transfert d'impuls, en 1/A

function DW(q,Z,Temp)

  use declarations
  implicit none

  integer:: Z

! fonction de Debye, facteur de Debye Waller
  real(kind=db):: Debye_function, Debye_Temperature, DW, FD, Mass_atom, q, T_Debye, Temp, x

  T_Debye = Debye_Temperature(Z)
  
  FD = Debye_function(Temp,T_Debye)

  x = - 300 * (q*hbar)**2 * FD / ( Mass_atom(Z) * atom_mu * k_Boltzmann * T_Debye )

  DW = exp(x)

  return
end

!*********************************************************************

function Debye_function(Temp,T_Debye)

  use declarations
  implicit none
  
  integer:: i,imax

  real(kind=db):: a, Debye_function, Temp, T_Debye, dx, x

  a = Temp / T_Debye
  dx = 0.0001_db
  imax = nint(a/dx)
  dx = a/imax
  x = dx / 2

  Debye_function = 0._db

  do i = 1,imax
   x = x + dx
   Debye_function = Debye_function + ( dx*x /( exp(x)-1 ) )
  end do

  Debye_function = 0.25_db + ( T_Debye / Temp )**2 * Debye_function 

end

