! Subroutine of the FDMNES Package
! Contain the main of the FDM-MST part of the calculation

subroutine fdm(Ang_borm,Bormann,Check_extract,Comt,Convolution_cal,Delta_edge,E_cut_imp,E_cut_man,Ecent,Elarg, &
        Epsii_ref,Epsii_ref_man,Estart,Fit_cal, &
        Gamma_hole,Gamma_hole_man,Gamma_max,Gamma_tddft,hkl_borm,icheck,ifile_notskip,indice_par,iscratch, &
        itape1,itape4,Length_line,mpinodes0,mpirank0,n_atom_proto_p,ngamh,ngroup_par,nnotskip,nnotskipm, &
        nomfich,nomfichbav,npar,nparm,param,Scan_a,typepar,Use_FDMX,FDMX_only, &
        fdmnes_inp,cm2g,nobg,nohole,nodw,noimfp,imfp_inp,imfp_infile,elf_inp,elf_infile,dwfactor_inp,dwfactor,tdebye_inp, &
        tdebye,tmeas_inp,tmeas,expntl,expntlA,expntlB,victoreen,victA,victB,mermrank)

  use declarations
  implicit none
  include 'mpif.h'

  integer, parameter:: n_rout = 22
  integer, parameter:: nslapw_max = 48  ! Max number of symmetry matrix for the FLAPW data

  integer:: i, igr, igr_dop, iord, ipr, ipr_dop, ir, is, iscratch, istat, it, itape1, itape4, itph, itpm, itps, &
    itype_dop, j, jgr, jseuil, l, lamstdens, lecrantage, Length_line, Lin_gam, Lmax_nrixs, Lmax_pot, &
    Lmax_tddft_inp, Lmaxat0, Lmaxso_max, Lmaxso0, Lseuil, m_hubb_e, mermrank, mpierr, mpinodes, mpinodes0, mpirank, mpirank0, &
    multi_run, multi_run_e, multrmax, n, n_abs_rgh, n_atom_bulk, n_atom_cap, n_atom_coop, n_atom_int, n_atom_neq, n_atom_per, &
    n_atom_per_neq, &
    n_atom_proto, n_atom_proto_bulk, n_atom_proto_p, n_atom_proto_uc, n_atom_sur, n_atom_uc, n_devide, n_file_dafs_exp, &
    n_mat_polar, n_multi_run, n_multi_run_e, natomsym_max, n_adimp, n_bulk_sup, n_q_rixs, n_radius, n_range, n_theta_in, &
    n_two_theta, n_Z_abs, nb_atom_conf_m, nbseuil, nchemin, ncolm, necrantage, neimagent, nenerg_in_rixs, nenerg_s, neqm, &
    ngamh, ngamme, ngroup, ngroup_hubb, ngroup_lapw, ngroup_m, ngroup_nonsph, ngroup_par, ngroup_pdb, &
    ngroup_taux, ngroup_temp, nhybm, nklapw, nlatm, nlm_pot, nlmlapwm, nmatsym, nnotskip, nnotskipm, norbdil, normrmt, &
    nparm, nphim, npldafs, npldafs_2d, npldafs_e, npldafs_f, nple, nplrm, nq_nrixs, nrato_dirac, nrm, nself, nself_bulk, nseuil, &
    nslapwm, nspin, nspino, nspinp, ntype, ntype_bulk, ntype_conf, ntype_mod, Trace_k, Wien_save, Z_nospinorbite

  integer, dimension(3):: hkl_borm, hkl_ref
  integer, dimension(8):: n_atom
  integer, dimension(12):: Tensor_imp
  integer, dimension(30):: icheck
  integer, dimension(ngroup_par):: npar
  integer, dimension(nnotskipm):: ifile_notskip
  integer, dimension(ngroup_par,nparm):: indice_par

  integer, dimension(:), allocatable:: iabsm, iabsorig, icom, icom_temp, igr_coop, igr_i, igr_is, igr_proto, itdil, its_lapw, &
     itype, Kgroup, ldil, n_bulk_sup_c, natomeq_s, ngreq, ngreqm, nlat, nlat_temp, norbv, nphi_dafs, npl_2d, nposextract, nrato, &
     nrato_temp, nrato_lapw, nrato_lapw_temp, nsymextract, numat, numat_temp, numat_abs, Z_bulk, Z_cap

  integer, dimension(:,:), allocatable:: igreq, isigpi, isymqa, lvval, lvval_temp, nvval, nvval_temp, Operation_mode
  integer, dimension(:,:,:), allocatable:: Wien_matsym

  character(len=8):: PointGroup
  character(len=9):: keyword
  character(len=132):: comt, identmot, imfp_infile, elf_infile, fdmnes_inp132, mot, nomfich132
  character(len=Length_name):: nomfich, nom_fich_extract, nomfichbav, fdmnes_inp

  character(len=6), dimension(n_rout):: Name_rout
  character(len=9), dimension(ngroup_par,nparm):: typepar
  character(len=Length_name), dimension(9):: Wien_file

  character(len=35), dimension(:), allocatable:: Com, Com_temp
  character(len=Length_name), dimension(:), allocatable:: nomfich_cal_conv, nomfich_cal_tddft_conv

  complex(kind=db), dimension(:,:), allocatable:: poldafsem, poldafssm
  complex(kind=db), dimension(:,:,:), allocatable:: hybrid

  logical:: Absauto, Abs_in_bulk, All_nrixs, Allsite, ATA, Atom_nonsph, Atom_occ_hubb, Atomic_scr, Axe_loc, &
     Basereel, Base_ortho_int, Base_ortho_sur, Bormann, Bulk, Check_extract, Clementi, Cap_layer, Cartesian_tensor, Center_s, &
     Charge_free, Classic_irreg, Convolution_cal, COOP, Coop_z_along_bond, Core_energ_tot, Core_resolved, Coupelapw, Dafs, &
     Dafs_bio, Density, Density_comp, Dip_rel, Dipmag, Doping, Dyn_eg, Dyn_g, E_cut_man, Eneg_i, Eneg_n_i, Energphot, &
     Epsii_ref_man, Extract, Extract_ten, FDM_comp, FDMX_only, &
     Film, Fit_cal, Flapw, Flapw_new, Force_ecr, Full_atom_e, Full_potential, Full_self_abs, Gamma_hole_man, Gamma_tddft, &
     Green_bulk, Green_int, Green_s, Green_self, Harm_cubic, Helm_cos, hkl_film, Hubbard, Kern_fast, key_calc, korigimp, &
     lmaxfree, lmoins1, lplus1, Magnetic, Matper, Memory_save, Moment_conservation, Muffintin, No_DFT, No_solsing, Noncentre, &
     Nonexc, Normaltau, &
     NRIXS, Occupancy_first, Octupole, Old_zero, One_run, One_SCF, Operation_mode_used, Optic, Overad, Pdb, PointGroup_Auto, &
     Polarise, Powder, Quadmag, Quadrupole, Readfast, Relativiste, Renorm, RIXS, RIXS_core, RIXS_fast, RPALF, Rydberg, Scan_a, &
     Self_abs, Solsing_s, Solsing_only, Spherical_signal, Spherical_tensor, Spin_channel, Spinorbite, State_all, &
     State_all_out, Supermuf, Sym_2D, Symauto, Symmol, Taux, Tddft, &
     Temp_B_iso, Trace_format_wien, Use_FDMX, Xan_atom, Ylm_comp_inp, &
     cm2g, nobg, nohole, nodw, noimfp, imfp_inp, elf_inp, dwfactor_inp, tdebye_inp, tmeas_inp, expntl, victoreen

  logical, dimension(7):: SCF_log
  logical, dimension(10):: Multipole

  logical, dimension(:), allocatable:: Atom_with_axe, Atom_nsph_e, Hubb, Hubb_temp, ok, Run_done, Skip, Skip_run

  real(kind=db):: alfpot, Ang_borm, Bulk_roughness, Cap_B_iso, Cap_roughness, Cap_shift, Cap_thickness, &
     D_max_pot, Delta_edge, Delta_En_conv, Delta_Epsii, Delta_helm, E_cut_imp, Ecent, Eclie, Eclie_out, Elarg, &
     Ephot_min, Epsii_ref, Estart, Film_roughness, Film_thickness, Gamma_max, Kern_fac, Overlap, p_self_max, p_self0, &
     p_self0_bulk, p_z1, p_z2, Pas_SCF, phi_0, R_rydb, R_self, R_self_bulk, Ray_max_dirac, Roverad, Rpotmax, Rtph, Step_loss, &
     Temp, Test_dist_min, V_helm, V_intmax, dwfactor, tdebye, tmeas, expntlA, expntlB, victA, victB, Width_helm

  real(kind=db), dimension(2):: Dist_coop, f_no_res
  real(kind=db), dimension(3):: Ang_rotsup, angxyz, angxyz_bulk, angxyz_cap, angxyz_int, angxyz_sur, axyz, axyz_bulk, axyz_cap, &
                                axyz_int, axyz_sur, Centre, dpos, Vec_orig
  real(kind=db), dimension(4):: Film_shift, Interface_shift, Surface_shift
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(7):: Time_loc
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(3,3):: Mat_or, Mat_UB
  real(kind=db), dimension(n_rout):: Time_rout
  real(kind=db), dimension(ngroup_par,nparm):: param

  real(kind=db), dimension(:), allocatable:: Adimp_e, Angle_or, cdil, E_adimp, E_max_range, E_radius, Ecrantage, &
     Eeient, Egamme, Eimagent, Energ_rixs, r0_lapw, r0_lapw_temp, Rchimp, rchimp_temp, Rlapw, rlapw_temp, rmt, rmt_temp, &
     Rmtimp, rmtimp_temp, Rsorte_s, Taux_cap, Taux_oc, Temp_coef, V_hubbard, V_hubbard_temp, V0bdcFimp
  real(kind=db), dimension(:,:), allocatable :: Angle_mode, angpoldafs, Axe_atom_gr, &
     hkl_dafs, pdpolar, polar, pop_nonsph, posn, posn_bulk, posn_cap, q_nrixs, q_rixs, &
     Vecdafsem, Vecdafssm, veconde, Wien_taulap
  real(kind=db), dimension(:), allocatable:: Theta_in, Two_theta
  real(kind=db), dimension(:,:,:), allocatable:: popats, popval, popval_temp,rot_atom_gr, rotloc_lapw
  real(kind=db), dimension(:,:,:,:), allocatable:: occ_hubb_e

  real(kind=sg) time

  data Name_Rout/'Lectur','Symsit','Reseau','Potent','Ylm   ','Potex ','Sphere','Mat   ','Fillin','Triang', &
                 'Tensor','Densit','Coabs ','Radial','Chi_0 ','Kern  ','Chi   ','SCF   ','Xanes ','Optic ', &
                 'Tddft ','Total '/

  Time_rout(:) = 0._db
  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    Time_loc(1) = real(time,db)
  endif

  call lectdim(Absauto,Atom_occ_hubb,Atom_nonsph,Axe_loc,Bormann,Bulk,Cap_layer,Dafs_bio,Doping,Extract,Extract_ten,Film, &
    Flapw,Full_self_abs,Hubbard,itape4,Length_line,Magnetic,Memory_save,mpinodes0,mpirank0,n_adimp,n_atom,n_atom_coop, &
    n_file_dafs_exp,n_mat_polar,n_multi_run_e,n_q_rixs,n_radius,n_range,n_theta_in,n_two_theta,n_Z_abs,nb_atom_conf_m,ncolm, &
    neimagent,nenerg_in_rixs,nenerg_s,ngamme, &
    ngroup,nhybm,nklapw,nlatm,nlmlapwm,nmatsym,norbdil,npldafs,npldafs_2d,npldafs_e,nple,nplrm,nq_nrixs, &
    NRIXS,nspin,nspino,nspinp,ntype,ntype_bulk,ntype_conf,Occupancy_first,Pdb,Readfast,Self_abs,Taux,Temp_B_iso,Use_FDMX,Xan_atom)

  n_atom_bulk = n_atom(1); n_atom_cap = n_atom(2); n_atom_int = n_atom(3); n_atom_neq = n_atom(4); n_atom_per = n_atom(5)
  n_atom_per_neq = n_atom(6); n_atom_sur = n_atom(7); n_atom_uc = n_atom(8)

  if( Atom_nonsph ) then
    ngroup_nonsph = ngroup
  else
    ngroup_nonsph = 0
  endif
  if( Flapw ) then
    ngroup_lapw = n_atom_per
  else
    ngroup_lapw = 0
  endif
  if( Magnetic .or. Atom_nonsph ) then
    ngroup_m = ngroup
  else
    ngroup_m = 0
  endif
  if( Taux ) then
    ngroup_taux = ngroup
  else
    ngroup_taux = 0
  endif
  if( Atom_occ_hubb ) then
    m_hubb_e = 3
    ngroup_hubb = ngroup
  else
    m_hubb_e = 0
    ngroup_hubb = 0
  endif
  if( Pdb ) then
    ngroup_pdb = ngroup
  else
    ngroup_pdb = 0
  endif
  if( Temp_B_iso ) then
    ngroup_temp = ngroup
  else
    ngroup_temp = 0
  endif
  if( Flapw ) then
    nslapwm = nslapw_max
  else
    nslapwm = 0
  endif
  if( Dafs_bio ) then
    npldafs_f = npldafs
  else
    npldafs_f = 0
  endif

  allocate( Adimp_e(n_adimp) )
  allocate( Angle_mode(3,npldafs_2d) )
  allocate( Angle_or(npldafs_f) )
  allocate( Angpoldafs(3,npldafs) )
  allocate( Atom_nsph_e(ngroup) )
  allocate( Axe_atom_gr(3,ngroup_m) )
  allocate( Com(0:ntype) )
  allocate( cdil(norbdil) )
  allocate( E_adimp(n_adimp) )
  allocate( E_max_range(n_range) )
  allocate( E_radius(n_radius) )
  allocate( Ecrantage(nspin) )
  allocate( Egamme(ngamme) )
  allocate( Eeient(neimagent) )
  allocate( Eimagent(neimagent) )
  allocate( Energ_rixs(nenerg_in_rixs) )
  allocate( hkl_dafs(3,npldafs) )
  allocate( Hubb(0:ntype) )
  allocate( hybrid(nhybm,16,ngroup_nonsph) )
  allocate( iabsm(n_multi_run_e) )
  allocate( iabsorig(n_multi_run_e) )
  allocate( icom(0:ntype) )
  allocate( igr_coop(n_atom_coop) )
  allocate( isigpi(2,npldafs) )
  allocate( itdil(norbdil) )
  allocate( its_lapw(ngroup_lapw) )
  allocate( itype(ngroup) )
  allocate( Kgroup(ngroup_pdb) )
  allocate( ldil(norbdil) )
  allocate( lvval(0:ntype,nlatm) )
  allocate( natomeq_s(n_radius) )
  allocate( nlat(0:ntype) )
  allocate( norbv(0:ngroup_nonsph) )
  allocate( nphi_dafs(npldafs) )
  allocate( npl_2d(npldafs_2d) )
  allocate( nposextract(n_multi_run_e) )
  allocate( nrato(0:ntype) )
  allocate( nrato_lapw(0:ntype) )
  allocate( nsymextract(n_multi_run_e) )
  allocate( numat(0:ntype) )
  allocate( numat_abs(n_Z_abs) )
  allocate( nvval(0:ntype,nlatm) )
  allocate( occ_hubb_e(-m_hubb_e:m_hubb_e,-m_hubb_e:m_hubb_e,nspin,ngroup_hubb) )
  allocate( Operation_mode(5,npldafs_2d) )
  allocate( pdpolar(nple,2) )
  allocate( polar(3,nple) )
  allocate( poldafsem(3,npldafs_e) )
  allocate( poldafssm(3,npldafs_e) )
  allocate( pop_nonsph(nhybm,ngroup_nonsph) )
  allocate( popats(ngroup,nlatm,nspin) )
  allocate( popval(0:ntype,nlatm,nspin) )
  allocate( posn(3,n_atom_uc) )
  allocate( posn_bulk(3,n_atom_bulk) )
  allocate( posn_cap(3,n_atom_cap) )
  allocate( q_nrixs(4,nq_nrixs) )
  allocate( q_rixs(4,n_q_rixs) )
  allocate( r0_lapw(0:ntype) )
  allocate( Rchimp(0:ntype) )
  allocate( Rlapw(0:ntype) )
  allocate( rmt(0:ntype) )
  allocate( Rmtimp(0:ntype) )
  allocate( Rot_Atom_gr(3,3,ngroup_m) )
  allocate( Rotloc_lapw(3,3,ngroup_lapw) )
  allocate( Rsorte_s(n_radius) )
  allocate( Taux_cap(n_atom_cap) )
  allocate( Taux_oc(ngroup_taux) )
  allocate( Temp_coef(ngroup_temp) )
  allocate( Theta_in(n_theta_in) )
  allocate( Two_theta(n_two_theta) )
  allocate( V_hubbard(0:ntype) )
  allocate( V0bdcFimp(nspin) )
  allocate( Vecdafsem(3,npldafs_e) )
  allocate( Vecdafssm(3,npldafs_e) )
  allocate( Veconde(3,nple) )
  allocate( Wien_matsym(3,3,nslapwm) )
  allocate( Wien_taulap(3,nslapwm) )
  allocate( Z_bulk(n_atom_bulk) )
  allocate( Z_cap(n_atom_cap) )

  call lecture(Absauto,adimp_e,alfpot,All_nrixs,Allsite,Ang_borm,Ang_rotsup,Angle_mode,Angle_or,Angpoldafs,Angxyz,Angxyz_bulk, &
    Angxyz_cap,Angxyz_int,Angxyz_sur,ATA,Atom_occ_hubb,Atom_nonsph,Atom_nsph_e,Atomic_scr,Axe_atom_gr,Axe_loc,axyz,axyz_bulk, &
    axyz_cap,axyz_int,axyz_sur,Basereel,Bormann,Bulk,Bulk_roughness,Cap_B_iso,Cap_layer,Cap_roughness, &
    Cap_shift,Cap_thickness,Cartesian_tensor,cdil,Center_s,Centre,Charge_free,Check_extract,Classic_irreg,Clementi,Com,Comt, &
    COOP,Coop_z_along_bond,Core_energ_tot,Core_resolved,Coupelapw,D_max_pot,Dafs,Dafs_bio,Delta_En_conv,Delta_Epsii, &
    Delta_helm,Density, &
    Density_comp,Dip_rel,Dipmag,Dist_coop,Doping,dpos,Dyn_eg,Dyn_g,E_adimp,E_radius,E_max_range,Eclie,Eclie_out,Ecrantage,Eeient, &
    Egamme,Eimagent,Eneg_i,Eneg_n_i,Energ_rixs,Energphot,Ephot_min,Extract,Extract_ten,f_no_res,FDM_comp,FDMX_only,Film, &
    Film_roughness,Film_shift,Film_thickness,Fit_cal,Flapw,Flapw_new,Force_ecr,Full_atom_e,Full_potential,Full_self_abs, &
    Gamma_hole,Gamma_hole_man,Gamma_max,Gamma_tddft,Green_bulk,Green_int,Green_s,Green_self,Harm_cubic,Helm_cos, &
    hkl_borm,hkl_dafs,hkl_film,hkl_ref,Hubb, &
    Hubbard,hybrid,iabsm,iabsorig,icheck,icom,igr_coop,igr_dop,indice_par,Interface_shift,iscratch,isigpi,itdil,its_lapw,iord, &
    itape4,itype,itype_dop,jseuil,Kern_fac,Kern_fast,Kgroup,korigimp,Lmax_nrixs,lamstdens,ldil,lecrantage,Length_line,Lin_gam, &
    Lmax_pot,Lmax_tddft_inp,lmaxfree,Lmaxso_max,Lmaxso0,Lmaxat0,lmoins1,lplus1,Lseuil,lvval,m_hubb_e,Magnetic,Mat_or,Mat_UB, &
    Matper,Moment_conservation,mpinodes,mpinodes0,mpirank,mpirank0,Muffintin,Multipole, &
    multrmax,n_adimp,n_atom,n_atom_bulk,n_atom_cap,n_atom_coop,n_atom_uc,n_atom_proto,n_devide,n_file_dafs_exp,n_mat_polar, &
    n_multi_run_e,n_q_rixs,n_radius,n_range,n_theta_in,n_two_theta,n_Z_abs,nb_atom_conf_m,nbseuil,nchemin,necrantage, &
    neimagent,nenerg_in_rixs,nenerg_s,ngamh,ngamme,ngroup,ngroup_hubb,ngroup_lapw,ngroup_m,ngroup_nonsph,ngroup_par,ngroup_pdb, &
    ngroup_taux,ngroup_temp,nhybm,nlat,nlatm,No_DFT,No_solsing,nom_fich_extract, &
    nomfich,nomfichbav,Noncentre, &
    Nonexc,norbdil,norbv,Normaltau,normrmt,npar,nparm,nphi_dafs, &
    nphim,npl_2d,npldafs,npldafs_2d,npldafs_e,npldafs_f,nple,nposextract,nq_nrixs,nrato,nrato_dirac,nrato_lapw,nrm, &
    nself,nself_bulk,nseuil,nslapwm,nspin,nsymextract,ntype,ntype_bulk,ntype_conf,ntype_mod,numat,numat_abs, &
    nvval,occ_hubb_e,Occupancy_first,Octupole,Old_zero,One_run,One_SCF,Operation_mode,Operation_mode_used,Optic,Overad,Overlap, &
    p_self_max,p_self0,p_self0_bulk,param,Pas_SCF,pdpolar,phi_0,PointGroup,PointGroup_Auto,Polar,Polarise,poldafsem,poldafssm, &
    pop_nonsph,popats,popval,posn,posn_bulk,posn_cap,Powder,q_nrixs,q_rixs,Quadmag,Quadrupole,R_rydb,R_self,R_self_bulk, &
    r0_lapw,Ray_max_dirac,Rchimp,Readfast,Relativiste,Renorm,RIXS,RIXS_core,RIXS_fast,Rlapw,Rmt,Rmtimp,Rot_Atom_gr,rotloc_lapw, &
    roverad,RPALF,rpotmax,Rydberg,Rsorte_s,SCF_log,Self_abs, &
    Solsing_s,Solsing_only,Spherical_signal,Spherical_tensor,Spin_channel,Spinorbite,State_all,State_all_out, &
    Step_loss,Supermuf,Surface_shift,Symauto,Symmol,Taux,Taux_cap,Taux_oc,Tddft,Temp,Temp_coef,Temp_B_iso, &
    Tensor_imp,Test_dist_min,Theta_in,Trace_format_wien,Trace_k,Trace_p,Two_theta,Typepar,Use_FDMX,V_helm,V_hubbard,V_intmax, &
    Vec_orig,Vecdafsem,Vecdafssm,Veconde,V0bdcFimp,Width_helm,Wien_file,Wien_matsym,Wien_save,Wien_taulap,Ylm_comp_inp,Z_bulk, &
    Z_cap,Z_nospinorbite)

    allocate( Com_temp(0:ntype_mod) )
    allocate( Hubb_temp(0:ntype_mod) )
    allocate( icom_temp(0:ntype_mod) )
    allocate( lvval_temp(0:ntype_mod,nlatm) )
    allocate( nlat_temp(0:ntype_mod) )
    allocate( nrato_temp(0:ntype_mod) )
    allocate( nrato_lapw_temp(0:ntype_mod) )
    allocate( numat_temp(0:ntype_mod) )
    allocate( nvval_temp(0:ntype_mod,nlatm) )
    allocate( popval_temp(0:ntype_mod,nlatm,nspin) )
    allocate( r0_lapw_temp(0:ntype_mod) )
    allocate( rchimp_temp(0:ntype_mod) )
    allocate( rlapw_temp(0:ntype_mod) )
    allocate( rmt_temp(0:ntype_mod) )
    allocate( rmtimp_temp(0:ntype_mod) )
    allocate( V_hubbard_temp(0:ntype_mod) )

    do it = 0,ntype_mod
      com_temp(it) = Com(it)
      Hubb_temp(it) = Hubb(it)
      icom_temp(it) = icom(it)
      lvval_temp(it,:) = lvval(it,:)
      nlat_temp(it) = nlat(it)
      nrato_temp(it) = nrato(it)
      nrato_lapw_temp(it) = nrato_lapw(it)
      numat_temp(it) = numat(it)
      nvval_temp(it,:) = nvval(it,:)
      popval_temp(it,:,:) = popval(it,:,:)
      r0_lapw_temp(it) = r0_lapw(it)
      rchimp_temp(it) = Rchimp(it)
      rlapw_temp(it) = Rlapw(it)
      rmt_temp(it) = rmt(it)
      rmtimp_temp(it) = Rmtimp(it)
      V_hubbard_temp(it) = V_hubbard(it)
    end do

    deallocate( Com, Hubb, icom, lvval, nlat, nrato, nrato_lapw, numat, nvval, popval, r0_lapw, Rchimp, Rlapw, rmt, &
                Rmtimp,V_hubbard )

    ntype = ntype_mod

    allocate( Com(0:ntype) )
    allocate( Hubb(0:ntype) )
    allocate( icom(0:ntype_mod) )
    allocate( lvval(0:ntype,nlatm) )
    allocate( nlat(0:ntype) )
    allocate( nrato(0:ntype) )
    allocate( nrato_lapw(0:ntype) )
    allocate( numat(0:ntype) )
    allocate( nvval(0:ntype,nlatm) )
    allocate( popval(0:ntype,nlatm,nspin) )
    allocate( r0_lapw(0:ntype) )
    allocate( Rchimp(0:ntype) )
    allocate( Rlapw(0:ntype) )
    allocate( rmt(0:ntype) )
    allocate( Rmtimp(0:ntype) )
    allocate( V_hubbard(0:ntype) )

    com(:) = com_temp(:)
    Hubb(:) = Hubb_temp(:)
    icom(:) = icom_temp(:)
    lvval(:,:) = lvval_temp(:,:)
    nlat(:) = nlat_temp(:)
    nrato(:) = nrato_temp(:)
    nrato_lapw(:) = nrato_lapw_temp(:)
    numat(:) = numat_temp(:)
    nvval(:,:) = nvval_temp(:,:)
    popval(:,:,:) = popval_temp(:,:,:)
    r0_lapw(:) = r0_lapw_temp(:)
    Rchimp(:) = rchimp_temp(:)
    Rlapw(:) = rlapw_temp(:)
    rmt(:) = rmt_temp(:)
    Rmtimp(:) = rmtimp_temp(:)
    V_hubbard(:) = V_hubbard_temp(:)

    deallocate( Com_temp, Hubb_temp, icom_temp, lvval_temp, nlat_temp, nrato_temp, nrato_lapw_temp, numat_temp, nvval_temp, &
                popval_temp, r0_lapw_temp, rchimp_temp, rlapw_temp, rmt_temp, rmtimp_temp, V_hubbard_temp )

!*** JDB Sept. 2016
  if( FDMX_only .and. mpirank0 == 0 ) then
    if( nseuil /= 1 ) then
      write(*,*) 'ERROR: FDMX is currently only compatible with K-edge spectra'
    else
      write(*,*) 'CALLING FDMX'
      if( n_atom_uc /= ngroup ) then
        write(*,*) ' n_atom_uc /= ngroup,   contact Yves !'
        stop
      endif
      fdmnes_inp132 = fdmnes_inp
      nomfich132 = nomfich
      call fdmx(fdmnes_inp132,nomfich132,cm2g,nobg,nohole,nodw,noimfp,Gamma_hole,Gamma_hole_man,E_cut_imp*rydb,E_cut_man, &
        imfp_inp, imfp_infile, elf_inp, elf_infile, dwfactor_inp, dwfactor, tdebye_inp, tdebye, tmeas_inp, tmeas, Energphot, &
        expntl, expntlA, expntlB, victoreen, victA, victB, mermrank, ngroup, posn, itype, ntype, numat)
    endif
    return
  endif
!*** JDB

  Scan_a = nphim > 1

  if( mpirank0 /= 0 .and. Extract_ten ) then

    deallocate( Angle_mode, Angle_or, angpoldafs, Axe_atom_gr )
    deallocate( Com, cdil, Ecrantage, hkl_dafs, hybrid, iabsm )
    deallocate( Egamme,  Eeient, Eimagent )
    deallocate( iabsorig, icom, isigpi, itdil, itype, its_lapw )
    deallocate( Kgroup, ldil, lvval, nlat )
    deallocate( norbv, nphi_dafs, npl_2d, nposextract, nrato, nrato_lapw )
    deallocate( nsymextract, numat, nvval, occ_hubb_e, Operation_mode, pdpolar )
    deallocate( polar, poldafsem, poldafssm, pop_nonsph )
    deallocate( popats, popval, posn, posn_bulk, posn_cap, r0_lapw )
    deallocate( Rchimp, Rlapw, rmt, Rmtimp, Rot_Atom_gr )
    deallocate( rotloc_lapw, Taux_oc, Temp_coef, Hubb, V_hubbard )
    deallocate( V0bdcFimp, Vecdafsem, Vecdafssm, veconde )
    deallocate( Wien_matsym, Wien_taulap, Z_bulk, Z_cap )

    return

  endif

  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    Time_loc(2) = real(time,db)
    Time_rout(1) = Time_loc(2) - Time_loc(1)
  endif

  n_atom_proto = 0
  neqm = 0

  allocate( Atom_with_axe(0:ngroup_m) )
  allocate( igr_i(ngroup) )
  allocate( igr_is(ngroup) )
  allocate( igr_proto(ngroup) )

  Sym_2D = Film .and. n_atom_per == 0

! Evaluation of the symmetries of the (magnetic) space group
  call Symsite(Absauto,angxyz,angxyz_int,angxyz_sur,Atom_with_axe,Atom_nonsph,Atom_nsph_e,Axe_atom_gr, &
        axyz,axyz_int,axyz_sur,Base_ortho_int,Base_ortho_sur,.false.,Center_s,Centre, &
        Doping,Extract,Flapw,iabsm,icheck(3),igr_i,igr_is,igr_proto,itype,Magnetic,Matper, &
        Memory_save,n_atom_int,n_atom_per,n_atom_proto,n_atom_sur,n_atom_uc,n_multi_run_e,n_Z_abs,neqm, &
        ngroup,ngroup_m,ngroup_taux,nlat,nlatm,Noncentre,nspin,ntype,numat,numat_abs,popats,posn,Sym_2D,Taux_oc)

  n_atom_proto_uc = n_atom_proto
  if( Doping ) n_atom_proto = n_atom_proto + 1
  n = n_atom_proto

  if( Bulk ) then
    Sym_2D = .false.
    call Symsite(Absauto,angxyz_bulk,angxyz_int,angxyz_sur,Atom_with_axe,Atom_nonsph,Atom_nsph_e,Axe_atom_gr, &
        axyz_bulk,axyz_int,axyz_sur,Base_ortho_int,Base_ortho_sur,Bulk,.false.,Centre, &
        .false.,Extract,Flapw,iabsm,icheck(3),igr_i,igr_is,igr_proto,itype,Magnetic,.true., &
        Memory_save,n_atom_int,n_atom_bulk,n_atom_proto,n_atom_sur,n_atom_bulk,n_multi_run_e,n_Z_abs,neqm, &
        ngroup,ngroup_m,ngroup_taux,nlat,nlatm,.false.,nspin,ntype,numat,numat_abs,popats,posn_bulk,Sym_2D,Taux_oc)
  endif

  deallocate( Atom_nsph_e )

  n_atom_proto_bulk = n_atom_proto - n

  allocate( igreq(0:n_atom_proto,neqm) )
  allocate( ngreq(0:n_atom_proto) )

  ngreq(:) = 0
  igreq(:,:) = 0

  do igr = 1,ngroup
    if( Doping .and. igr == n_atom_uc + 1 ) cycle
    ipr = igr_proto(igr)
    ngreq(ipr) = ngreq(ipr) + 1
    i = igr_i(igr)
    igreq(ipr,i) = igr
  end do

  if( Doping ) then
    igreq(n,1) = n_atom_uc + 1
    ngreq(n) = 1
  endif

  if( Absauto ) then
    multi_run = 0
    do i = 1,n_Z_abs
      do ipr = 1,n_atom_proto
        it = abs( itype(igreq(ipr,1)) )
        if( numat( it ) /= numat_abs(i) ) cycle
        multi_run = multi_run + 1
      end do
    end do
    n_multi_run = multi_run
  else
    n_multi_run = n_multi_run_e
  endif

  if( n_multi_run > 1 .and. One_run .and. .not. nonexc .and. mpirank0 == 0 ) then
    call write_error
    do ipr = 3,9,3
      write(ipr,125)
    end do
    stop
  endif

  if( Absauto ) then
    deallocate( iabsm, iabsorig )
    deallocate( nsymextract, nposextract )
    allocate( iabsm(n_multi_run) )
    allocate( iabsorig(n_multi_run) )
    allocate( nsymextract(n_multi_run) )
    allocate( nposextract(n_multi_run) )
    nsymextract(:) = 1
    do multi_run = 1,n_multi_run
      nposextract(multi_run) = multi_run
    end do
    multi_run = 0
    do i = 1,n_Z_abs
      do ipr = 1,n_atom_proto
        it = abs( itype(igreq(ipr,1)) )
        if( numat( it ) /= numat_abs(i) ) cycle
        multi_run = multi_run + 1
        iabsm(multi_run) = igreq(ipr,1)
        iabsorig(multi_run) = ipr
      end do
    end do
  endif

  natomsym_max = 0
  do i = 1,n_Z_abs
    do igr = 1,ngroup
      if( Doping ) then
        if( numat( abs( itype(igr) ) ) /= numat( abs( itype(igr_dop) ) ) ) cycle
      else
        if( numat( abs( itype(igr) ) ) /= numat_abs(i) ) cycle
      endif
      natomsym_max = max( igr_i(igr), natomsym_max )
    end do
  end do

  allocate( isymqa(0:n_atom_proto,natomsym_max) )
  isymqa(:,:) = 0

  ipr_dop = 0

  do j = 1,n_Z_abs
    do igr = 1,ngroup

      if( Doping ) then
        if( numat( abs( itype(igr) ) ) /= numat( abs( itype(igr_dop) ) ) ) cycle
      else
        if( numat( abs( itype(igr) ) ) /= numat_abs(j) ) cycle
      endif

      if( igr == igr_dop ) ipr_dop = igr_proto(igr)

      ipr = igr_proto(igr)
      i = igr_i(igr)
      if( Doping .and. igr == n_atom_uc + 1 ) then
        isymqa(n_atom_proto_uc+1,:) = isymqa(ipr_dop,:)
      else
        isymqa(ipr,i) = igr_is(igr)
      endif

    end do
  end do

  allocate( n_bulk_sup_c(n_multi_run) )
  n_bulk_sup_c(:) = 0

  Abs_in_bulk = .false.
  do multi_run = 1,n_multi_run
    if( iabsm(multi_run) > ngroup - n_atom_bulk ) then
      Abs_in_bulk = .true.
      exit
    endif
  end do

  if( Coop .and. n_atom_coop == 0 ) then
    n_atom_coop = n_multi_run
    deallocate( igr_coop )
    allocate( igr_coop(n_atom_coop) )
    do i = 1,n_atom_coop
      igr_coop(i) = iabsm(i)
    end do
  endif

! Determination of supplementary output files when abs in bulk with different z
  n_abs_rgh = 0
  do multi_run_e = 1,n_multi_run
    if( Abs_in_bulk ) then
      multi_run = n_multi_run + 1 - multi_run_e
    else
      multi_run = multi_run_e
    endif
    if( iabsm(multi_run) <= ngroup - n_atom_bulk ) cycle
    if( Bulk_roughness > eps10 ) n_abs_rgh = n_abs_rgh + 1
    n_bulk_sup_c(multi_run) = - 1
    ipr = igr_proto( iabsm(multi_run) )
    allocate( ok( ngreq(ipr) ) )
    ok(:) = .false.
    do igr = 1,ngreq(ipr)
      if( ok(igr) ) cycle
      ok(igr) = .true.
      n_bulk_sup_c(multi_run) = n_bulk_sup_c(multi_run) + 1
      p_z1 = posn_bulk( 3, igreq(ipr,igr) - ngroup + n_atom_bulk )
      do jgr = 1,ngreq(ipr)
        if( ok(jgr) ) cycle
        p_z2 = posn_bulk( 3, igreq(ipr,jgr) - ngroup + n_atom_bulk )
        if( abs( p_z1 - p_z2 ) > eps10 ) cycle
        ok(jgr) = .true.
      end do
    end do
    deallocate( ok )
  end do
  n_bulk_sup = sum( n_bulk_sup_c(:) )

  deallocate( igr_i, igr_is, igr_proto )

  allocate( ngreqm(n_multi_run) )

  allocate( nomfich_cal_conv(n_multi_run+n_bulk_sup+n_abs_rgh) )
  allocate( nomfich_cal_tddft_conv(n_multi_run+n_bulk_sup+n_abs_rgh) )
  allocate( Run_done(n_multi_run) )
  allocate( Skip_run(n_multi_run) )
  allocate( Skip(n_multi_run) )
  Run_done(:) = .false.
  Skip_run(:) = .false.
  Skip(:) = .false.

  boucle_m: do multi_run = 1,n_multi_run

    boucle_ipr: do ipr = 1,n_atom_proto
      do i = 1,ngreq(ipr)
        if( igreq(ipr,i) == iabsm(multi_run) ) exit boucle_ipr
      end do
    end do boucle_ipr

    ngreqm(multi_run) = ngreq(ipr)

    do j = 1,multi_run - 1
      do i = 1,ngreq(ipr)
        if( iabsm(j) /= igreq(ipr,i) ) cycle
        run_done(multi_run) = .true.
        cycle boucle_m
      end do
    end do

  end do boucle_m

  if( nnotskip > 0 .and. Extract .and. n_atom_proto_p == n_atom_proto ) then
    do multi_run = 1,n_multi_run
      Skip_run(multi_run) = .true.
      do i = 1,nnotskip
        if( ifile_notskip(i) == iabsorig(multi_run) ) Skip_run(multi_run) = .false.
      end do
    end do
  endif

  n_atom_proto_p = n_atom_proto
  nlm_pot = ( Lmax_pot + 1 )**2

  if( mpirank0 == 0 ) write(6,130) n_multi_run

! ------- Calculation of all the non equivalent atoms signal ----------------------------------------

  skip(:) = Run_done(:) .or. Skip_run(:)

  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    Time_loc(3) = real(time,db)
    Time_rout(2) = Time_loc(3) - Time_loc(2)
  endif

  call Site_calculation(adimp_e,alfpot,All_nrixs,Allsite,Ang_rotsup,Angle_mode,Angle_or,Angpoldafs,Angxyz,Angxyz_bulk, &
      Angxyz_cap,Angxyz_int,Angxyz_sur,ATA,Atom_occ_hubb,Atom_nonsph,Atom_with_axe,Atomic_scr,Axe_atom_gr,axyz,axyz_bulk, &
      axyz_cap,axyz_int,axyz_sur,Base_ortho_int,Base_ortho_sur,Basereel,Bormann,Bulk,Bulk_roughness,Cap_B_iso, &
      Cap_layer,Cap_roughness,Cap_shift,Cap_thickness,Cartesian_tensor,cdil,Center_s,Centre,Charge_free,Classic_irreg, &
      Clementi,Com,COOP,Coop_z_along_bond,Core_energ_tot,Core_resolved,Coupelapw,D_max_pot,Dafs,Dafs_bio,Delta_edge,Delta_En_conv, &
      Delta_Epsii,Delta_helm,Density,Density_comp,Dist_coop,Dip_rel,Dipmag,Doping,dpos,Dyn_eg,Dyn_g,E_adimp,E_cut_imp,E_cut_man, &
      E_radius,E_max_range,Ecent,Eclie,Eclie_out,Ecrantage,Eeient,Egamme,Eimagent,Elarg,Eneg_i,Eneg_n_i,Energ_rixs,Energphot, &
      Ephot_min,Epsii_ref,Epsii_ref_man,Estart,Extract,Extract_ten,f_no_res, &
      FDM_comp,Film,Film_roughness,Film_shift,Film_thickness,Flapw,Flapw_new, &
      Force_ecr,Full_atom_e,Full_potential,Full_self_abs,Gamma_hole,Gamma_hole_man,Gamma_max,Gamma_tddft,Green_bulk,Green_int, &
      Green_s,Green_self,Harm_cubic,Helm_cos,hkl_dafs,hkl_film,hkl_ref,Hubb,Hubbard,hybrid,iabsm,iabsorig,icheck,icom,igr_coop, &
      igr_dop,igreq,Interface_shift,ipr_dop,iscratch,isigpi,itdil,its_lapw,iord,isymqa,itype,itype_dop,jseuil, &
      Kern_fac,Kern_fast,Kgroup,korigimp,Lmax_nrixs,lamstdens,ldil,lecrantage,Lin_gam, &
      Lmax_pot,Lmax_tddft_inp,lmaxfree,Lmaxso_max,Lmaxso0,Lmaxat0,lmoins1,lplus1,Lseuil,lvval,m_hubb_e, &
      Magnetic,Mat_or,Mat_UB,Matper,Moment_conservation,mpinodes,mpinodes0,mpirank,mpirank0,Muffintin,Multipole,multrmax, &
      n_abs_rgh,n_adimp,n_atom_bulk,n_atom_cap,n_atom_coop,n_atom_int,n_atom_per,n_atom_proto,n_atom_proto_bulk,n_atom_proto_uc, &
      n_atom_sur,n_atom_uc,n_bulk_sup,n_devide,n_mat_polar,n_multi_run,n_q_rixs,n_radius,n_range,n_rout,n_theta_in,n_two_theta, &
      natomeq_s,natomsym_max,nbseuil,nchemin,ncolm,necrantage, &
      neimagent,nenerg_in_rixs,nenerg_s,neqm,ngamh,ngamme,ngreq,ngroup,ngroup_hubb,ngroup_lapw, &
      ngroup_m,ngroup_nonsph,ngroup_pdb,ngroup_taux,ngroup_temp,nhybm,nklapw,nlm_pot,nlmlapwm,nlat,nlatm,nmatsym,No_DFT, &
      No_solsing,nom_fich_extract,nomfich,nomfich_cal_conv,nomfich_cal_tddft_conv,Noncentre, &
      Nonexc,norbdil,norbv,Normaltau,normrmt,nphi_dafs,nphim, &
      npl_2d,npldafs,npldafs_2d,npldafs_e,npldafs_f,nple,nplrm,nposextract,nq_nrixs,nrato,nrato_dirac,nrato_lapw,NRIXS,nrm, &
      nself,nself_bulk,nseuil,nslapwm,nspin,nspino,nspinp,nsymextract,ntype,numat,numat_abs(1), &
      nvval,occ_hubb_e,Octupole,Old_zero,One_run,One_SCF,Operation_mode,Operation_mode_used,Optic,Overad,Overlap,p_self_max, &
      p_self0,p_self0_bulk,Pas_SCF,pdpolar,phi_0,PointGroup,PointGroup_Auto,Polar,Polarise,poldafsem,poldafssm, &
      pop_nonsph,popats,popval,posn,posn_bulk,posn_cap,Powder,q_nrixs,q_rixs,Quadrupole,R_rydb,R_self,R_self_bulk, &
      r0_lapw,Ray_max_dirac,Rchimp,Relativiste,Renorm,RIXS,RIXS_core,RIXS_fast,Rlapw,rmt,Rmtimp,Rot_Atom_gr,rotloc_lapw, &
      roverad,RPALF,rpotmax,Rydberg,Rsorte_s,SCF_log,Self_abs,Skip,Solsing_s,Solsing_only, &
      Spherical_signal,Spherical_tensor,Spin_channel,Spinorbite,State_all,State_all_out,Step_loss,Supermuf,Surface_shift,Symmol, &
      Taux,Taux_cap,Taux_oc,Tddft,Temp,Temp_coef,Temp_B_iso,Tensor_imp,Test_dist_min,Theta_in,Time_rout,Trace_format_wien, &
      Trace_k,Trace_p,Two_theta,V_helm,V_hubbard,V_intmax,Vec_orig,Vecdafsem,Vecdafssm,Veconde,V0bdcFimp,Width_helm,Wien_file, &
      Wien_matsym,Wien_save,Wien_taulap,Xan_atom,Ylm_comp_inp,Z_bulk,Z_cap,Z_nospinorbite)

! -------------------------------------------------------------------------------------------------------------------------------------------

  if( Convolution_cal .and. mpirank0 == 0 ) then

    Rewind(itape1)
    key_calc = .false.
    do l = 1,10000
      read(itape1,'(A)',iostat=istat) mot
      if( istat /= 0 ) then
        backspace(itape1)
        exit
      endif
      keyword = identmot(mot,9)
      if( keyword == 'calculati' ) key_calc = .true.
      if( keyword == 'Run_done' ) then
        backspace(itape1)
        exit
      endif
    end do

    if( key_calc ) then
      write(itape1,'(A)') ' Run_done '
      do multi_run_e = 1,n_multi_run
        if( Abs_in_bulk ) then
          multi_run = n_multi_run + 1 - multi_run_e
        else
          multi_run = multi_run_e
        endif
        do i = 0,n_bulk_sup_c(multi_run)
          if( Skip_run(multi_run) ) then
            is = 0
          else
            is = 1
          endif
          if( Run_done(multi_run) ) then
            ir = 0
          else
            ir = 1
          endif
          write(itape1,'(2i2)') ir, is
        end do
      end do
    else
      write(itape1,'(A)') ' calculation'
      n = 0
      do multi_run_e = 1,n_multi_run
        if( Abs_in_bulk ) then
          multi_run = n_multi_run + 1 - multi_run_e
        else
          multi_run = multi_run_e
        endif
        do i = 0,n_bulk_sup_c(multi_run)
          n = n + 1
          if( Run_done(multi_run) ) cycle
          write(itape1,'(A)') nomfich_cal_conv(n)
        end do
        if( iabsm(multi_run) > ngroup - n_atom_bulk .and. Bulk_roughness > eps10 ) then
          n = n + 1
          if( Run_done(multi_run) ) cycle
          write(itape1,'(A)') nomfich_cal_conv(n)
        endif
      end do
      if( Tddft .and. .not. Extract .and. ( .not. Gamma_tddft .or. Dafs  .or. n_multi_run > 1 ) ) then
        write(itape1,'(A)') ' cal_tddft'
        n = 0
        do multi_run_e = 1,n_multi_run
          if( Abs_in_bulk ) then
            multi_run = n_multi_run + 1 - multi_run_e
          else
            multi_run = multi_run_e
          endif
          do i = 0,n_bulk_sup_c(multi_run)
            n = n + 1
            if( Run_done(multi_run) ) cycle
            write(itape1,'(A)') nomfich_cal_tddft_conv(n)
          end do
          if( iabsm(multi_run) > ngroup - n_atom_bulk .and. Bulk_roughness > eps10 ) then
            n = n + 1
            if( Run_done(multi_run) ) cycle
            write(itape1,'(A)') nomfich_cal_tddft_conv(n)
          endif
        end do
      endif
    endif
  endif

  deallocate( n_bulk_sup_c, nomfich_cal_conv, nomfich_cal_tddft_conv, Run_done, Skip_run )

  if( icheck(1) > 0 .and. mpirank0 == 0 ) then
    call CPU_TIME(time)
    Time_loc(4) = real(time,db)
    Time_rout(n_rout) = Time_loc(4) - Time_loc(1)
    itph = int( Time_rout(n_rout) / 3600 )
    rtph = Time_rout(n_rout) - itph * 3600
    itpm = int( rtph / 60 )
    itps = nint( rtph - itpm * 60 )
    write(3,210)
    if( Tddft) then
      write(3,220) (Name_rout(i), Time_rout(i), i = 1,17)
    else
      write(3,220) (Name_rout(i), Time_rout(i), i = 1,13)
    endif
    if( RIXS ) write(3,220) 'RIXS  ', Time_rout(20)
    if( Optic .and. Tddft ) then
      write(3,220) Name_rout(18),  Time_rout(18), 'Taull ', Time_rout(19), ( Name_rout(i),  Time_rout(i), i = 20,21 )
    elseif( Optic ) then
      write(3,220) 'Taull ', Time_rout(19), Name_rout(20),  Time_rout(20)
    elseif( Tddft ) then
      write(3,220) ( Name_rout(i),  Time_rout(i), i = 18,19 ), Name_rout(21),  Time_rout(21)
    else
      write(3,220) ( Name_rout(i),  Time_rout(i), i = 18,19 )
    endif
    if( itph > 0 .or. itpm > 0 ) then
      write(3,240) Name_rout(n_rout), Time_rout(n_rout), itph, itpm, itps
    else
      write(3,220) Name_rout(n_rout), Time_rout(n_rout)
    endif
  endif

!*** JDB Sept. 2016
  if( Use_FDMX .and. mpirank0 == 0 ) then
    if( nseuil /= 1 ) then
      write(*,*) 'ERROR: FDMX is currently only compatible with K-edge spectra'
    else
      write(*,*) 'CALLING FDMX'
      if( n_atom_uc /= ngroup ) then
        write(*,*) ' n_atom_uc /= ngroup, contact Yves !'
        stop
      endif
      fdmnes_inp132 = fdmnes_inp
      nomfich132 = nomfich
      call fdmx(fdmnes_inp132,nomfich132,cm2g,nobg,nohole,nodw,noimfp,Gamma_hole,Gamma_hole_man,E_cut_imp*rydb,E_cut_man, &
        imfp_inp, imfp_infile, elf_inp, elf_infile, dwfactor_inp, dwfactor, tdebye_inp, tdebye, tmeas_inp, tmeas, Energphot, &
        expntl, expntlA, expntlB, victoreen, victA, victB, mermrank, ngroup, posn, itype, ntype, numat)
    endif
  endif
!*** JDB

  deallocate( Atom_with_axe )

! Desallocation des tables allocated before subroutine lecture
  deallocate( Adimp_e, Angle_mode, Angle_or, Angpoldafs )
  deallocate( Axe_atom_gr )
  deallocate( Com, cdil )
  deallocate( E_adimp, E_max_range, E_radius )
  deallocate( Ecrantage )
  deallocate( Egamme,  Eeient, Eimagent, Energ_rixs )
  deallocate( hkl_dafs, Hubb, hybrid )
  deallocate( iabsm, iabsorig, icom, igr_coop, isigpi, itdil, its_lapw, itype )
  deallocate( Kgroup, ldil, lvval, natomeq_s )
  deallocate( nlat, norbv, nphi_dafs, npl_2d, nposextract, nrato, nrato_lapw, nsymextract, numat, numat_abs, nvval )
  deallocate( occ_hubb_e, Operation_mode )
  deallocate( pdpolar, polar, poldafsem, poldafssm, pop_nonsph, popats, popval, posn, posn_bulk, posn_cap, q_nrixs, q_rixs )
  deallocate( r0_lapw, Rchimp, Rlapw, rmt, Rmtimp, Rot_Atom_gr, Rotloc_lapw, Rsorte_s )
  deallocate( Taux_cap, Taux_oc, Temp_coef, Theta_in, Two_theta )
  deallocate( V_hubbard, V0bdcFimp, Vecdafsem, Vecdafssm, Veconde )
  deallocate( Wien_matsym, Wien_taulap, Z_bulk, Z_cap )

  return
  125 format(//' One_run process is not possible with absorbing atom in excited state !',// &
               ' Excited absorbing atom is the default for K, L1, M1, N1 edges.', // &
               ' Use non excited absorbing atom (keyword nonexc)',/ &
               ' or keep a multi run process.'//)
  130 format(/' Number of calculated non equivalent absorbing atom =', i3)
  210 format(/,1x,120('-')//' Subroutine times (s CPU)')
  220 format(4(4x,a6,' =',f12.2))
  240 format(4x,a6,' =',f12.2,' =',i4,' h',i3,' min',i3,' s CPU')
end

!*************************************************************************************************************************

! Calculation of all the non equivalent atoms signal

subroutine Site_calculation(adimp_e,alfpot,All_nrixs,Allsite,Ang_rotsup,Angle_mode,Angle_or,Angpoldafs,Angxyz,Angxyz_bulk, &
      Angxyz_cap,Angxyz_int,Angxyz_sur,ATA,Atom_occ_hubb,Atom_nonsph,Atom_with_axe,Atomic_scr,Axe_atom_gr,axyz,axyz_bulk, &
      axyz_cap,axyz_int,axyz_sur,Base_ortho_int, Base_ortho_sur,Basereel,Bormann,Bulk,Bulk_roughness,Cap_B_iso, &
      Cap_layer,Cap_roughness,Cap_shift,Cap_thickness,Cartesian_tensor,cdil,Center_s,Centre,Charge_free,Classic_irreg, &
      Clementi,Com,COOP,Coop_z_along_bond,Core_energ_tot,Core_resolved,Coupelapw,D_max_pot,Dafs,Dafs_bio,Delta_edge,Delta_En_conv, &
      Delta_Epsii,Delta_helm,Density,Density_comp,Dist_coop,Dip_rel,Dipmag,Doping,dpos,Dyn_eg,Dyn_g,E_adimp,E_cut_imp,E_cut_man, &
      E_radius,E_max_range,Ecent,Eclie,Eclie_out,Ecrantage,Eeient,Egamme,Eimagent,Elarg,Eneg_i,Eneg_n_i,Energ_rixs,Energphot, &
      Ephot_min,Epsii_ref,Epsii_ref_man,Estart,Extract,Extract_ten,f_no_res, &
      FDM_comp,Film,Film_roughness,Film_shift,Film_thickness,Flapw,Flapw_new, &
      Force_ecr,Full_atom_e,Full_potential,Full_self_abs,Gamma_hole,Gamma_hole_man,Gamma_max,Gamma_tddft,Green_bulk,Green_int, &
      Green_s,Green_self,Harm_cubic,Helm_cos,hkl_dafs,hkl_film,hkl_ref,Hubb,Hubbard,hybrid,iabsm,iabsorig,icheck,icom,igr_coop, &
      igr_dop,igreq,Interface_shift,ipr_dop,iscratch,isigpi,itdil,its_lapw,iord,isymqa,itype,itype_dop,jseuil, &
      Kern_fac,Kern_fast,Kgroup,korigimp,Lmax_nrixs,lamstdens,ldil,lecrantage,Lin_gam, &
      Lmax_pot,Lmax_tddft_inp,lmaxfree,Lmaxso_max,Lmaxso0,Lmaxat0,lmoins1,lplus1,Lseuil,lvval,m_hubb_e, &
      Magnetic,Mat_or,Mat_UB,Matper,Moment_conservation,mpinodes,mpinodes0,mpirank,mpirank0,Muffintin,Multipole,multrmax, &
      n_abs_rgh,n_adimp,n_atom_bulk,n_atom_cap,n_atom_coop,n_atom_int,n_atom_per,n_atom_proto,n_atom_proto_bulk,n_atom_proto_uc, &
      n_atom_sur,n_atom_uc,n_bulk_sup,n_devide,n_mat_polar,n_multi_run,n_q_rixs,n_radius,n_range,n_rout,n_theta_in,n_two_theta, &
      natomeq_s,natomsym_max,nbseuil,nchemin,ncolm,necrantage, &
      neimagent,nenerg_in_rixs,nenerg_s,neqm,ngamh,ngamme,ngreq,ngroup,ngroup_hubb,ngroup_lapw, &
      ngroup_m,ngroup_nonsph,ngroup_pdb,ngroup_taux,ngroup_temp,nhybm,nklapw,nlm_pot,nlmlapwm,nlat,nlatm,nmatsym,No_DFT, &
      No_solsing,nom_fich_extract,nomfich,nomfich_cal_conv,nomfich_cal_tddft_conv,Noncentre, &
      Nonexc,norbdil,norbv,Normaltau,normrmt,nphi_dafs,nphim, &
      npl_2d,npldafs,npldafs_2d,npldafs_e,npldafs_f,nple,nplrm,nposextract,nq_nrixs,nrato,nrato_dirac,nrato_lapw,NRIXS,nrm, &
      nself,nself_bulk,nseuil,nslapwm,nspin,nspino,nspinp,nsymextract,ntype,numat,numat_abs_1, &
      nvval,occ_hubb_e,Octupole,Old_zero,One_run,One_SCF,Operation_mode,Operation_mode_used,Optic,Overad,Overlap,p_self_max, &
      p_self0,p_self0_bulk,Pas_SCF,pdpolar,phi_0,PointGroup,PointGroup_Auto,Polar,Polarise,poldafsem,poldafssm, &
      pop_nonsph,popats,popval,posn,posn_bulk,posn_cap,Powder,q_nrixs,q_rixs,Quadrupole,R_rydb,R_self,R_self_bulk, &
      r0_lapw,Ray_max_dirac,Rchimp,Relativiste,Renorm,RIXS,RIXS_core,RIXS_fast,Rlapw,rmt,Rmtimp,Rot_Atom_gr,rotloc_lapw, &
      roverad,RPALF,rpotmax,Rydberg,Rsorte_s,SCF_log,Self_abs,Skip,Solsing_s,Solsing_only, &
      Spherical_signal,Spherical_tensor,Spin_channel,Spinorbite,State_all,State_all_out,Step_loss,Supermuf,Surface_shift,Symmol, &
      Taux,Taux_cap,Taux_oc,Tddft,Temp,Temp_coef,Temp_B_iso,Tensor_imp,Test_dist_min,Theta_in,Time_rout,Trace_format_wien, &
      Trace_k,Trace_p,Two_theta,V_helm,V_hubbard,V_intmax,Vec_orig,Vecdafsem,Vecdafssm,Veconde,V0bdcFimp,Width_helm,Wien_file, &
      Wien_matsym,Wien_save,Wien_taulap,Xan_atom,Ylm_comp_inp,Z_bulk,Z_cap,Z_nospinorbite)

  use declarations
  implicit none

  include 'mpif.h'

  integer, parameter:: nslapw_max = 48  ! Max number of symmetry matrix for the FLAPW data
  integer, parameter:: n_tens_dd = 9
  integer, parameter:: n_tens_dq = 15
  integer, parameter:: n_tens_qq = 25
  integer, parameter:: n_tens_dm = 9

  integer:: extract_nenerg, extract_nlmomax, i, i_adimp, i_radius, i_range, i_self, ia, iaabs, iaabsfirst, iabsfirst, &
    iaabsi, iabsorbeur, &
    iapr, ip0, ip00, ipr, iaprabs, iaprabs_nonexc, ich, ie, ie_computer, ie0, igr_dop, igrph, igrpt_nomag, igrpt0, &
    index_e, initl, iord, ip_max, ipr1, ipr_dop, iprabs, iprabs_nonexc, ir, iscratch, iso1, iso2, isp1, isp2, ispin, iss1, &
    iss2, it, itab, itabs, itabs_nonexc, itype_dop, j, je, jseuil, l, l_hubbard, l0_nrixs, &
    l1, l2, lamstdens, lecrantage, lh, Lin_gam, lla_state, lla2_state, Lm1, Lm2, Lmax, Lmax_abs, Lmax_abs_a, Lmax_abs_t, &
    Lmax_comp, Lmax_g, Lmax_nrixs, Lmax_pot, Lmax_probe, Lmax_tddft, Lmax_tddft_inp, Lmax_max, Lmaxat0, Lmaxso, &
    Lmaxso_m, Lmaxso_max, Lmaxso0, Lms1, Lms2, Lseuil, m, m_hubb, m_hubb_e, m1, m2, &
    mpierr, mpinodee, mpinodee0, mpinodes, mpinodes0, mpirank, mpirank0, mpirank_in_mumps_group, &
    multi_0, multi_imp, multi_run, multi_run_e, multrmax, n, n_abs_rgh, n_adimp, n_atom_0, n_atom_0_self, n_atom_bulk, &
    n_atom_cap, n_atom_coop, n_atom_ind, n_atom_ind_self, n_atom_int, n_atom_per, n_atom_proto, n_atom_proto_bulk, &
    n_atom_proto_uc, n_atom_sur, n_atom_uc, n_bulk_sup, n_bulk_z, n_bulk_z_abs, n_bulk_z_max, n_bulk_z_max_abs, n_devide, n_Ec, &
    n_mat_polar, n_max, n_multi_run, n_oo, n_orb_rel, n_orbexc, n_q_rixs, n_radius, n_range, n_rel, n_rout, n_tens_max, &
    n_theta_in, n_two_theta, n_V, &
    n_vr_0, n_vr_ind, n1, nab_coop, natome, natome_cal, natome_self, natomsym_max, natomeq, natomeq_coh, natomeq_self, &
    natomp, natomsym, nb_rep, nb_sym_op, nbbseuil, nbm, nbseuil, nbtm, nbtm_fdm, nchemin, ncolm, ncolr, ncolt, &
    ndim2, necrantage, neimagent, nenerg, nenerg_coh, nenerg_in_rixs, nenerg_rixs, nenerg_s, nenerg_tddft, nenerg0, neqm, ngamh, &
    ngamme, nge, ngrm, ngroup, ngroup_hubb, ngroup_lapw, ngroup_m, ngroup_nonsph, ngroup_pdb, &
    ngroup_taux, ngroup_temp, ngrph, nhybm, nicm, nim, ninit1, ninit, ninit_out, ninitr, &
    ninitlv, nklapw, nlatm, nlm, nlm_comp, nlm_pot, nlm_probe, nlm_p_fp, nlm_rixs, nlmagm, nlmamax, &
    nlmlapwm, nlmmax, nlmomax, nlms_pr, nlmsam, nlmsamax, nlmso, nmatsym, nnlm, norbdil, normrmt, &
    nphiato1, nphiato7, nphim, npldafs, npldafs_2d, npldafs_e, npldafs_f, nple, nplr, nplrm, npoint, npoint_ns, npr, &
    npso, npsom, nptmoy, nptmoy_out, nq_nrixs, nr, nrato_abs, nrato_dirac, nrm, nrm_self, ns_dipmag, nself, nself_bulk, nself_c, &
    nseuil, nslapwm, nso1, nsort, nsortf, nsm, nspin, nspino, nspinp, nspinorb, nstm, ntype, numat_abs, numat_abs_1, nvois, nx, &
    nxanout, Trace_k, Trace_l, Wien_save, Z, Z_nospinorbite

  integer, dimension(3):: hkl_ref, ldip
  integer, dimension(12):: Tensor_imp
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3):: msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(30):: icheck
  integer, dimension(n_atom_coop):: ia_coop, igr_coop
  integer, dimension(n_multi_run):: iabsm, iabsorig, nposextract, nsymextract
  integer, dimension(0:ntype):: icom, nlat, nrato, nrato_lapw, numat
  integer, dimension(0:n_atom_proto,neqm):: igreq
  integer, dimension(nopsm):: iopsymc, iopsymr
  integer, dimension(nrepm,2):: irep_util(nrepm,2)
  integer, dimension(2,npldafs):: isigpi
  integer, dimension(0:n_atom_proto,natomsym_max):: isymqa
  integer, dimension(norbdil):: itdil, ldil
  integer, dimension(ngroup_lapw):: its_lapw
  integer, dimension(ngroup):: itype
  integer, dimension(ngroup_pdb):: Kgroup
  integer, dimension(0:ntype,nlatm):: lvval, nvval
  integer, dimension(0:n_atom_proto):: iapot, itypepr, Lmaxat, ngreq
  integer, dimension(0:n_atom_proto,0:mpinodes-1):: lmax_dft
  integer, dimension(0:ngroup_nonsph):: norbv
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(npldafs_2d):: npl_2d
  integer, dimension(5,npldafs_2d):: Operation_mode
  integer, dimension(3,3,nslapwm):: Wien_matsym
  integer, dimension(n_atom_bulk):: Z_bulk
  integer, dimension(n_atom_cap):: Z_cap
  integer, dimension(n_radius):: natomeq_s
  integer, dimension(2,0:ntype):: lcoeur, ncoeur

  integer, dimension(:), allocatable:: ia_eq_inv, ia_eq_inv_self, iaproto, iaprotoi, &
     igroup, igroupi, imoy, imoy_out, is_g, ispin_maj, isrt, isymeq, &
     itypei, itypep, Lmaxa, Lval, n_bulk_zc, n_bulk_zc_abs, nb_eq, nb_eq_2D, nbord, nbordf, nlmsa, nlmso0, numia

  integer, dimension(:,:), allocatable:: ia_eq, ibord, igr_bulk_z, igr_bulk_z_abs, index_coop, indice, iopsym_atom, is_eq, &
     isbord, isvois, ivois, iso, lso, m_g, mso, nb_rpr, nlmsa0

  integer, dimension(:,:,:), allocatable:: ia_rep, iato, lato, mato, mpres, nb_rep_t

  integer, dimension(:,:,:,:), allocatable:: msymoo, msymooi

  character(len=5):: Struct
  character(len=8):: PointGroup
  character(len=Length_name):: nomfich, nom_fich_extract, nomfich_s

  character(len=35), dimension(0:ntype):: Com
  character(len=Length_name), dimension(9):: Wien_file
  character(len=Length_name), dimension(n_multi_run+n_bulk_sup+n_abs_rgh):: nomfich_cal_conv, nomfich_cal_tddft_conv
  character(len=Length_name), dimension(:), allocatable:: File_name_rixs

  character(len=5), dimension(:), allocatable:: ltypcal
  character(len=Length_word), dimension(:), allocatable:: nomabs

  complex(kind=db):: f_avantseuil, f_cal
  complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdf0t_rgh_bulk, phdt_rgh_Bulk
  complex(kind=db), dimension(3,npldafs_e):: poldafsem, poldafssm
  complex(kind=db), dimension(nhybm,16,ngroup_nonsph):: hybrid

  complex(kind=db), dimension(:,:), allocatable:: Bragg_abs, Bragg_rgh_bulk_abs, karact, pol, Sum_Bragg_bulk_abs_f0, &
                                                  Sum_Bragg_nonabs_f
  complex(kind=db), dimension(:,:,:), allocatable:: phdt, poldafse, poldafss
  complex(kind=db), dimension(:,:,:,:), allocatable:: secmd, secmd_m, secmm, secmm_m, V_hubb_abs, &
                                                      V_hubb_t
  complex(kind=db), dimension(:,:,:,:,:), allocatable:: Ampl_dft, Ampl_abs, S_nrixs, S_nrixs_m, secdd, secdd_m, Tau_comp, rof0, &
                                                        secdq, secdq_m, Tau_ato, Taull, Taull_dft, Taull_tdd, V_hubb, V_hubb_s
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: secdo, secdo_m, secoo, secoo_m, secqq, secqq_m, Tau_coop, Taull_stk
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: RadialIntegral, Taull_abs

  logical:: Abs_in_bulk, Abs_in_Bulk_roughness, Absorbeur, All_nrixs, Allsite, ATA, Bulk_atom_done, Atom_comp_cal, Atom_nonsph, &
     Atom_nonsph_loc, Atom_occ_hubb, Atomic_scr, Basereel, Base_hexa, Base_ortho_int, Base_ortho_sur, Bormann, Bulk, Bulk_step, &
     Clementi, Cal_xanes, Cap_layer, Cartesian_tensor, Center_s, Charge_free, Classic_irreg, &
     Convergence, COOP, Coop_z_along_bond, Core_energ_tot, Core_resolved, Coupelapw, Dafs, Dafs_bio, Density, &
     Density_comp, Devide_Ei, Dip_rel, Dipmag, Doping, Dyn_eg, Dyn_g, E_comp, E1E1, E1E2, E1E2e, E1E3, E1M1, E2E2, E3E3, &
     E_cut_man, Eneg, Eneg_i, Eneg_n_i, Energphot, Epsii_ref_man, Extract, &
     Extract_ten, FDM_comp, FDM_comp_m, Fermi, Fermi_first, Film, Final_optic, Final_tddft, First_E, First_Bulk_step, &
     Flapw, Flapw_new, Force_ecr, Full_atom, Full_atom_e, Full_atom_gen, Full_potential, Full_self_abs, &
     Gamma_hole_man, Gamma_tddft, Green, Green_bulk, Green_i, Green_int, Green_s, Green_self, &
     Harm_cubic, Helm_cos, hkl_film, Hubb_a, Hubb_abs, Hubb_d, Hubb_diag_abs, Hubbard, Kern_fast, &
     korigimp, Level_val_abs, Level_val_exc, lmaxfree, lmoins1, lplus1, M1M1, Magnetic, Matper, Moment_conservation, Moy_cluster, &
     Moy_loc, Moyenne, Muffintin, No_DFT, No_solsing, Noncentre, Nonexc_g, Nonexc, Nonsph, Normaltau, NRIXS, Octupole, &
     Old_zero, One_run, One_SCF, Operation_mode_used, Optic, Optic_xanes, Overad, PointGroup_Auto, Polarise, Powder, Proto_all, &
     Quadrupole, Recop, Relativiste, Renorm, Renorm_SCF, RIXS, RIXS_core, RIXS_fast, RPALF, Rydberg, SCF, SCF_bulk, SCF_c, &
     SCF_elecabs, &
     SCF_mag_fix, SCF_mag_free, Second_run, Second_SCF, Self_abs, Self_cons, Self_cons_bulk, Self_cons_c, Self_nonexc, &
     Solsing, Solsing_s, Solsing_o, Solsing_only, Spherical_signal, Spherical_tensor, Spin_channel, Spino, Spinorbite, State_all, &
     State_all_out, Supermuf, Sym_2D, Sym_4, Sym_cubic, Symmol, Tau_nondiag, Taux, Tddft, &
     Tddft_xanes, Temp_B_iso, Trace_format_wien, Xan_atom, Ylm_comp, Ylm_comp_inp

  logical, dimension(7):: SCF_log
  logical, dimension(10):: Multipole
  logical, dimension(n_multi_run):: Skip
  logical, dimension(0:ntype):: Hubb
  logical, dimension(0:ngroup_m):: Atom_with_axe
  logical, dimension(0:n_atom_proto):: Proto_calculated

  logical, dimension(:), allocatable:: Atom_axe, Hubb_diag, Repres_comp

  real(kind=db):: Abs_U_iso, adimp, alfpot, Bulk_roughness, Cal_Volume_maille, Cap_B_iso, Cap_roughness, &
    Cap_shift, Cap_thickness, Chargm, Chargm_SCF, d_ecrant, D_max_pot, Delta_bulk, Delta_cap, &
    Delta_E_conv, Delta_edge, Delta_En_conv, Delta_energ, Delta_energ_s, Delta_Epsii, Delta_Eseuil, &
    Delta_film, Delta_helm, Delta_int, Delta_roughness_film, Delta_sur, delta_z_bottom_cap, delta_z_top_bulk, delta_z_top_cap, &
    delta_z_top_film, E, E_cut, E_cut_imp, E_Fermi, E_Open_val, E_Open_val_exc, E_start, &
    E_zero, Ecent, Ecineticmax, Ecineticmax_out, Eclie, Eclie_out, Ecmax, Ecmax_out, Ei0, Eii, Elarg, Em, En_cluster, &
    En_cluster_s, Enervide, Ephot_min, Epsii_moy, Epsii_ref, Estart, Extract_E_cut, Extract_E_Fermi, Extract_V0bdcF, &
    Film_roughness, Film_thickness, Gamma_max, Kern_fac, Overlap, p_self, p_self_max, p_self0, p_self0_c, p_self0_bulk, Pas_SCF, &
    phi_0, R_SCF_c, R_rydb, R_self, R_self_bulk, Ray_max_dirac, Rmax, Rmtg_abs, Rmtsd_abs, Roverad, Rpotmax, rsbdc, rsbdc_out, &
    Rsort, Rsorte, Step_loss, &
    Surface_ref, Temp, Test_dist_min, Time_fill, Time_tria, tp_SCF_1, tp_SCF_2, tp_XANES_1, tp_XANES_2, V_helm, V_helm_t, &
    V_intmax, V0muf, Vhbdc, Vhbdc_init, Vhbdc_out, Volume_maille, Vsphere, Vsphere_cal, Width_helm, Workf, Workf_i, Workf_val

  real(kind=db), dimension(2):: chg_open_val, Dist_coop, f_no_res, pop_open_val
  real(kind=db), dimension(3):: Ang_rotsup, Angxyz, angxyz_bulk, angxyz_cap, angxyz_int, angxyz_sur, axyz, axyz_bulk, axyz_cap, &
                                axyz_int, axyz_sur, Centre, deccent, dpos, Vec_orig
  real(kind=db), dimension(4):: Film_shift, Interface_shift, Surface_shift
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(7):: Time_loc
  real(kind=db), dimension(10):: Gamma_hole
  real(kind=db), dimension(3,3):: Cubmat, Mat_or, Mat_UB, Orthmat, Orthmati, Orthmatt, Rot_atom_abs, Rot_int
  real(kind=db), dimension(n_rout):: Time_rout

  real(kind=db), dimension(n_adimp):: Adimp_e
  real(kind=db), dimension(npldafs_f):: Angle_or
  real(kind=db), dimension(3,npldafs):: Angpoldafs
  real(kind=db), dimension(3,npldafs_2d):: Angle_mode
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_gr
  real(kind=db), dimension(4,n_q_rixs):: q_rixs
  real(kind=db), dimension(0:n_atom_proto):: Chargat, Rmtg, Rmtg0, Rmtsd
  real(kind=db), dimension(norbdil):: cdil
  real(kind=db), dimension(n_adimp):: E_adimp
  real(kind=db), dimension(n_range):: E_max_range
  real(kind=db), dimension(n_radius):: E_radius
  real(kind=db), dimension(n_theta_in):: Theta_in
  real(kind=db), dimension(n_two_theta):: Two_theta
  real(kind=db), dimension(nspin):: dV0bdcF, Chg_reference, Ecinetic, Ecinetic_out, Ecrantage, V0bdc, V0bdcF, V0bdcF_out, &
                                    V0bdc_out, V0bdcFimp, VxcbdcF, VxcbdcF_out
  real(kind=db), dimension(neimagent):: Eeient, Eimagent
  real(kind=db), dimension(ngamme):: Egamme
  real(kind=db), dimension(nenerg_in_rixs):: Energ_rixs
  real(kind=db), dimension(npldafs):: Length_abs
  real(kind=db), dimension(3,npldafs):: hkl_dafs
  real(kind=db), dimension(3,npldafs_e):: Vecdafsem, Vecdafssm
  real(kind=db), dimension(-m_hubb_e:m_hubb_e,-m_hubb_e:m_hubb_e,nspin,ngroup_hubb):: occ_hubb_e
  real(kind=db), dimension(nple,2):: pdpolar
  real(kind=db), dimension(3,nple):: polar, Veconde
  real(kind=db), dimension(nhybm,ngroup_nonsph):: pop_nonsph
  real(kind=db), dimension(ngroup,nlatm,nspin):: popats
  real(kind=db), dimension(0:ntype,nlatm):: popatv
  real(kind=db), dimension(0:ntype,nlatm,nspin):: popval
  real(kind=db), dimension(3,n_atom_uc):: posn
  real(kind=db), dimension(3,n_atom_bulk):: posn_bulk
  real(kind=db), dimension(3,n_atom_cap):: posn_cap
  real(kind=db), dimension(4,nq_nrixs):: q_nrixs
  real(kind=db), dimension(0:ntype):: r0_lapw, Rchimp, Rlapw, rmt, Rmtimp, V_hubbard
  real(kind=db), dimension(3,3,ngroup_m):: Rot_Atom_gr
  real(kind=db), dimension(3,3,ngroup_lapw):: Rotloc_lapw
  real(kind=db), dimension(n_radius):: Rsorte_s
  real(kind=db), dimension(n_atom_cap):: Taux_cap
  real(kind=db), dimension(ngroup_taux):: Taux_oc
  real(kind=db), dimension(ngroup_temp):: Temp_coef
  real(kind=db), dimension(3,nslapwm):: Wien_taulap
  real(kind=db), dimension(3,ngroup_m):: Axe_atom_grn
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(0:nrm,2,0:ntype):: psi_coeur
  real(kind=db), dimension(nrm,nbseuil):: psii
  real(kind=db), dimension(0:nrm,nlatm,0:ntype):: psival
  real(kind=db), dimension(0:n_atom_proto,nlatm,nspin):: popatm
  real(kind=db), dimension(nrm,2):: psi_open_val
  real(kind=db), dimension(0:nrm,0:ntype):: Rato, rhoit, rho_coeur

  real(kind=db), dimension(:), allocatable:: Chg_coeur, cgrad, clapl, dExc_ex_nex, dista, distai, dvc_ex_nex, &
     E_KS_core, E_starta, Eimag, Eimag_coh, Eimag_s, Energ, Energ_coh, Energ_s, Energ_self, &
     Energ_self_s, Enervide_t, Epsii, Exc_abs_i, Length_rel, Length_rel_abs, poidso, poidsov, poidsov_out, r, &
     rato_abs, rsato_abs, rhons, rs, rvol, sec_atom, Taux_eq, Taux_ipr, Vc_abs_i, Vh, Vhns
  real(kind=db), dimension(:,:), allocatable:: coef_g, Axe_atom_clu, Chargat_init, chargat_init_bulk, Chargat_self, &
     chargat_self_bulk, chargat_self_s, drho_ex_nex, dv_ex_nex, Excato, Int_tens, pdp, poidsa, pop_orb_val, pop_orb_val_bulk, pos, &
     posi, posi_self, rho, rhoato_abs,rsato, V_abs_i, Vcato_abs, Vcato_init, Vcato_init_bulk, Vr, Vxc, vec, xyz, Ylmso
  real(kind=db), dimension(:,:,:), allocatable:: gradvr, Int_statedens, rho_chg, rho_chg_bulk, rho_self_t, rhoato, rhoato_init, &
     rhoato_init_bulk, rot_atom, Vcato, Vcato_bulk, Vecdafse, Vecdafss, Vra, Vrato_e, Vxcato_abs, Ylmato
  real(kind=db), dimension(:,:,:,:), allocatable:: rho_self, rho_self_bulk, rho_self_s, Vrato, Vxcato
  real(kind=db), dimension(:,:,:,:,:), allocatable:: Coverlap, drho_self, Occ_hubb, Occ_hubb_i
  real(kind=db), dimension(:,:,:,:,:,:), allocatable:: phiato, Statedens, Statedens_i

  real(kind=sg):: time

  call MPI_Comm_Rank(MPI_COMM_MUMPS,mpirank_in_mumps_group,mpierr)

  Delta_bulk = 0._db
  Delta_energ_s = 0._db
  First_Bulk_step = .true.
  Full_atom_gen = .true.
  ia_coop(:) = 0
! Moy_cluster is a test variable on the interatomic average potentiel calculation
! When false (and Green), the average is calculated in the close surrounding of the absorbing atom.
! Otherwise it is in all the area (cluster) of calculation
! True was found better
  Moy_cluster = .true.
  iaabsfirst = 0
  E1E2e = Multipole(2)
  SCF_elecabs = SCF_log(1); SCF_mag_free = SCF_log(2); Self_cons = SCF_log(3); Self_nonexc = SCF_log(4); SCF = SCF_log(5)
  Self_cons_bulk = SCF_log(6); SCF_bulk = SCF_log(7)
  Renorm_SCF = .true.
  Surface_ref = 0._db

  if( Dip_rel .and. Lseuil == 0 ) then ! For K edges, core state are pure spin states
    n_rel = 8
  elseif( Dip_rel ) then
    n_rel = 16
  else
    n_rel = 1
  endif

  l0_nrixs = 0
  if( Dipmag ) then
    ip0 = 0
  else
    ip0 = 1
  endif
  if( Octupole ) then
    ip_max = 3
  elseif( Quadrupole ) then
    ip_max = 2
  else
    ip_max = 1
  endif

! Lmax for density of state
  lla_state = max( 4, Lmax_pot )
  lla2_state = ( lla_state + 1 )**2

  ich = icheck(6)
! Film_thickness is reduced to be the distance between top most and botom most atoms
  if( Film ) call Prep_Film(angxyz,angxyz_bulk,angxyz_cap,angxyz_int,angxyz_sur,axyz,axyz_bulk,axyz_cap,axyz_int,axyz_sur, &
               Bulk_roughness,Cap_shift,Cap_roughness,Cap_thickness,Delta_bulk,Delta_cap,Delta_int,Delta_film, &
               Delta_roughness_film,Delta_sur,delta_z_bottom_cap,delta_z_top_bulk,delta_z_top_cap,delta_z_top_film,Film_shift, &
               Film_roughness,Film_thickness,hkl_film,ich,Interface_shift,itype,n_atom_bulk,n_atom_cap,n_atom_int,n_atom_per, &
               n_atom_sur,n_atom_uc,ngroup,ntype,numat,posn,posn_bulk,posn_cap,Surface_ref,Surface_shift)

! When Extract there is no parallelization in the energy loop.
! Parallelization remains for tddft
  if( Extract ) then
    mpinodee = 1
    mpinodee0 = 1
  else
    mpinodee = mpinodes
    mpinodee0 = mpinodes0
  endif

  if( Optic ) then
    ninit1 = 0; ninit = 0
  else
    call dim_init(jseuil,Lseuil,nbseuil,ninit1,ninit)
  endif
  allocate( coef_g(ninit,2) )
  allocate( m_g(ninit,2) )
  allocate( is_g(ninit) )
  if( .not. Optic ) call coef_init(coef_g,is_g,ninit,jseuil,Lseuil,m_g,nbseuil)

!-------------------------------------------------

  if( Extract ) then
    if( mpinodes > 1 ) then
    call MPI_Bcast(Green_s,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
    endif
  else
    allocate( Energ_s(nenerg_s) )
    allocate( Eimag_s(nenerg_s) )
    call grille_xanes(eeient,Eimag_s,eimagent,Egamme,Energ_s,icheck(4),Lin_gam,ngamme,neimagent,nenerg_s)
  endif

  if( Multipole(7) ) then
    n_oo = 9
  else
    n_oo = 0
  endif
  allocate( msymoo(3,n_oo,3,n_oo) )
  allocate( msymooi(3,n_oo,3,n_oo) )

! Variable to know if an absorbing atom is in the bulk
! In that case, there is an energy dependant truncature factor calculated in coabs/convolution and not in prepdafs.
  Abs_in_bulk = .false.
  do multi_run = 1,n_multi_run
    iabsorbeur = iabsm(multi_run)
    if( iabsorbeur > ngroup - n_atom_bulk ) then
      Abs_in_bulk = .true.
      exit
    endif
  end do

  multi_0 = 0

!- Loop over non equivalent atoms ---------------------------------------------------------------------------

  boucle_multi: do multi_run_e = 1,n_multi_run

    if( Abs_in_bulk ) then
      multi_run = n_multi_run + 1 - multi_run_e
    else
      multi_run = multi_run_e
    endif

    if( Skip(multi_run) ) then
      multi_0 = multi_0 + 1
      cycle
    endif

    if( icheck(4) > 0 ) write(3,110) iabsorig(multi_run)

    Second_run = multi_run_e > 1 .and. One_run
    Second_SCF = multi_run_e > 1 .and. One_SCF

    multi_imp = nposextract(multi_run_e)

    if( Spinorbite ) then
      Tau_nondiag = .true.
    else
      Tau_nondiag = .false.
    endif
    iabsorbeur = iabsm(multi_run)
    if( Doping ) then
      itabs_nonexc = abs( itype_dop )
    else
      itabs_nonexc = abs(itype(iabsorbeur))
    endif
    numat_abs = numat( itabs_nonexc )

    call Esdata(Eseuil,icheck(4),jseuil,nbseuil,nseuil,numat_abs,mpirank0)

! ninitr become 1 for the Tddft cycle.
    if( Core_resolved ) then
      if( Optic ) then
       ! In optics, there is one output per "l" plus the total.
        if( numat_abs_1 <= 4 ) then
          ninitr = 1
        elseif( numat_abs_1 > 4 .and. numat_abs_1 <= 20 ) then
          ninitr = 2
        elseif( numat_abs_1 > 21 .and. numat_abs_1 <= 56 ) then
          ninitr = 3
        else
          ninitr = 4
        endif
        ninitr = ninitr + 1
      else
        ninitr = ninit
      endif
    else
      ninitr = nbseuil
    endif

    allocate( Epsii(ninitr) )
    Epsii(:) = 0._db
    Epsii_moy = 0._db

    if( Old_zero .and. .not. Flapw ) then
      Workf_i = Workf_val(icheck(4),numat_abs)
    else
      Workf_i = 0._db
    endif

    nomfich_s = nomfich
    Final_tddft = .false.
    Final_optic = .false.

! Maximum number of atomic orbitals
    nnlm = 0
    do it = 1,ntype
      n = n_orb_rel( numat(it) )
      if( numat(it) == numat_abs )then
        n = n + 2 * max( 1, nlat(it) ) + 2
      elseif( nlat(it) > 0 ) then
        n = n + 2 * nlat(it)
      endif
      nnlm = max( nnlm, n )
    end do

    Bulk_step = Bulk .and.  iabsorbeur > ngroup - n_atom_bulk
    Sym_2D = iabsorbeur > n_atom_per .and. .not. Bulk_step .and. .not. Doping
    Abs_in_Bulk_roughness = Bulk_step .and. Bulk_roughness > eps10
    Bulk_atom_done = Sym_2D .and. Abs_in_bulk .and. .not. Full_atom_gen .and. SCF .and. Self_nonexc

    if( Bulk_step ) then
      SCF_c = SCF_bulk
      Self_cons_c = Self_cons_bulk
      p_self0_c = p_self0_bulk
      nself_c = nself_bulk
    else
      SCF_c = SCF
      Self_cons_c = Self_cons
      p_self0_c = p_self0
      nself_c = nself
    endif

    if( Sym_2D .or. .not. Matper ) then
      V_helm_t = V_helm
    else
      V_helm_t = 0._db
    endif
    call init_run(cdil,Chargat,Charge_free,Chargm,Chargm_SCF,Clementi,Com,Doping,Ecrantage,Flapw,Force_ecr,Hubb,iabsorbeur, &
      iabsorig(multi_run),icheck(4),icom,igreq,iprabs,iprabs_nonexc,itabs,itdil,itype,itype_dop,itypepr,jseuil,lcoeur,ldil, &
      lecrantage,Lseuil,lvval,mpinodee0,mpirank0,n_atom_proto,n_multi_run,n_orbexc,nbseuil,ncoeur,necrantage,neqm,ngreq,ngroup, &
      nlat,nlatm,nnlm,nomfich_s,Nonexc,norbdil,nrato,nrato_dirac,nrato_lapw,nrm,nseuil,nspin,ntype,numat,numat_abs,nvval, &
      pop_open_val,popatm,popats,popatv,popval,psi_coeur,psii,psi_open_val,psival,r0_lapw,Rato,Ray_max_dirac,Rchimp,Relativiste, &
      rho_coeur,rhoit,Rlapw,Rmt,Rmtimp,Self_nonexc,V_hubbard,Wien_file(8))

    nrato_abs = nrato(itabs)

    if( Doping ) then
      natomsym = ngreq(ipr_dop)
    else
      natomsym = ngreq(iprabs_nonexc)
    endif

    allocate( isymeq(natomsym) )
    allocate( Taux_eq(natomsym) )

    if( Doping ) then
      isymeq(1:natomsym) = isymqa(ipr_dop,1:natomsym)
    else
      isymeq(1:natomsym) = isymqa(iprabs_nonexc,1:natomsym)
    endif

    if( Taux ) then
      Taux_eq(1:natomsym) = Taux_oc( abs( igreq(iprabs_nonexc,1:natomsym) ) )
    else
      Taux_eq(1:natomsym) = 1._db
    endif

    if( Temp_B_iso ) then
      Abs_U_iso = Temp_coef( abs( igreq(iprabs_nonexc,1) ) ) / ( 8 * pi**2 )
    else
      Abs_U_iso = 0._db
    endif

    if( .not. Extract_ten ) then
      if( Hubbard ) then
        m_hubb = 0
        do it = 0,ntype
          lh = l_hubbard(numat(it))
          if( Hubb(it) ) m_hubb = max( m_hubb, lh)
        end do
      else
        m_hubb = 0
      endif
      allocate( V_hubb_abs(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp) )
    end if

    if( Extract_ten ) then

      if( .not. Old_zero ) E_Fermi = extract_E_Fermi(multi_imp,nom_fich_extract)
      E_cut = extract_E_cut(multi_imp,nom_fich_extract)
      if( Optic ) then
        V0muf = 0._db
        Close(1)
      else
        V0bdcF(1) = extract_V0bdcF()
        if( Old_zero ) then
          V0muf = WorkF_i + V0bdcF(1)
        else
          V0muf = WorkF_i - E_Fermi + V0bdcF(1)
        endif
        call extract_Epsii(Core_resolved,Delta_Epsii,Delta_Eseuil,Epsii,Epsii_moy,icheck(14),nbseuil, &
                 ninit1,ninit,ninitr)
      endif
      nenerg_s = extract_nenerg(multi_imp,nom_fich_extract,.false.,Tddft)
      allocate( Energ_s(nenerg_s) )
      allocate( Eimag_s(nenerg_s) )
      if( Tddft ) then
        nbbseuil = 1
      else
        nbbseuil = nbseuil
      endif
      call extract_energ(Energ_s,Eseuil,multi_imp,nbbseuil,nbseuil,nenerg_s,nom_fich_extract,.false.,Tddft)

    elseif( Extract ) then

      allocate( rato_abs(nrato_abs) )
      allocate( rsato_abs(nrato_abs) )
      allocate( rhoato_abs(nrato_abs,nspin) )
      allocate( Vcato_abs(nrato_abs,nlm_pot) )
      allocate( Vxcato_abs(nrato_abs,nlm_pot,nspin) )
      call Potential_reading(Delta_Eseuil,dV0bdcf,E_cut,E_Fermi,Ecineticmax,Epsii,Epsii_moy,Hubb_abs,Hubb_diag_abs, &
                             m_hubb,mpinodes,mpirank0,multi_run_e,ninitr,nlm_pot,nom_fich_extract,nrato_abs,nspin,nspinp, &
                             rato_abs,rhoato_abs,Rmtg_abs,Rmtsd_abs,rsato_abs,rsbdc,V_hubb_abs,V0muf,Vcato_abs,Vhbdc, &
                             Vxcato_abs,VxcbdcF)

      if( multi_run_e == 1 ) then
        if( mpirank0 == 0 ) nenerg_s = extract_nenerg(multi_imp,nom_fich_extract,Optic,.false.)
        if( mpinodes > 1 ) then
          call MPI_Bcast(nenerg_s,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
        endif
        allocate( Energ_s(nenerg_s) )
        allocate( Eimag_s(nenerg_s) )
        if( mpirank0 == 0 ) call extract_energ(Energ_s,Eseuil,multi_imp,nbseuil,nbseuil,nenerg_s,nom_fich_extract,Optic,.false.)
        if( mpinodes > 1 ) then
          call MPI_Bcast(Energ_s,nenerg_s,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(Eimag_s,nenerg_s,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
        endif
      endif

! Just useful for iaabsi
      Rmtsd(:) = Rmtsd_abs
      Rmtg(:) = Rmtg_abs

    endif

    if( Doping ) then
! Used just for the atom position
      iabsfirst = igr_dop
    elseif( One_run ) then
      iabsfirst = iabsm(1)
    else
      iabsfirst = iabsorbeur
    endif

    if( Bulk_step ) then
      call Cal_Cubmat(angxyz_bulk,Cubmat,Struct)
    else
      call Cal_Cubmat(angxyz,Cubmat,Struct)
    endif

    if( ( .not. One_run ) .or. multi_run_e == 1 ) then
      d_ecrant = 1 - sum( Ecrantage(:) )
      if( Bulk_step ) then
        R_SCF_c = R_self_bulk
      else
        R_SCF_c = R_self
      endif
      call natomp_cal(angxyz,angxyz_bulk,angxyz_int,angxyz_sur,ATA,axyz,axyz_bulk,axyz_int,axyz_sur,Base_ortho_int, &
        Base_ortho_sur,Bulk,Bulk_step,Center_s,Centre,Chargat,d_ecrant,deccent,Delta_bulk,Delta_film,Delta_int,Delta_sur,Doping, &
        dpos,Film_shift,Film_thickness,Flapw,iabsorbeur,iabsfirst,icheck(5),igr_dop,igreq,Interface_shift, &
        itabs,itype,Kgroup,Matper,mpirank0,multi_run_e,multrmax,n_atom_bulk,n_atom_int,n_atom_per,n_atom_proto,n_atom_proto_uc, &
        n_atom_sur,n_atom_uc,n_radius,natomeq_s,natomeq_coh,natomp,neqm,ngreq,ngroup,ngroup_pdb,ngroup_taux, &
        Noncentre,posn,posn_bulk,One_run,Proto_all,R_SCF_c,rsorte_s,rmax,rpotmax,Self_cons_c,Surface_shift,Sym_2D, &
        Taux_oc)
    endif

    allocate( Axe_atom_clu(3,natomp) )
    allocate( iaproto(natomp) )
    allocate( igroup(natomp) )
    allocate( itypep(natomp) )
    allocate( dista(natomp) )
    allocate( pos(3,natomp) )
    allocate( karact(nopsm,nrepm) )

    ich = icheck(5)

! Le groupe ponctuel de l'agregat est calcule avec l'absorbeur eventuellement excite meme pour le calcul SCF non excite.
! L'abaissement de symmetrie n'est effectif que si l'absorbeur est non centre
! De meme "Atom_nonsph" n'est utile que pour le Xanes
    call agregat(angxyz,angxyz_bulk,angxyz_int,angxyz_sur,ATA,Atom_with_axe,Atom_nonsph,Axe_atom_clu,Axe_atom_gr,Axe_atom_grn, &
          axyz,axyz_bulk,axyz_int,axyz_sur,Base_hexa,Base_ortho_int,Base_ortho_sur,Bulk,Bulk_step,Center_s, &
          Chargat,Cubmat,deccent,Delta_bulk,Delta_film,Delta_int,Delta_sur,dista, &
          Doping,dpos,Film_shift,Film_thickness,iaabs,iaabsfirst,iabsorbeur,iabsfirst, &
          iaproto,iapot,ich,igr_dop,igreq,igroup,igrpt_nomag,igrpt0,Interface_shift,iopsymc,iopsymr,itabs,itype,itypep, &
          karact,Kgroup,Magnetic,Matper,mpirank0,multi_run_e,n_atom_bulk,n_atom_int,n_atom_per,n_atom_proto,n_atom_sur,n_atom_uc, &
          natomp,nb_rep,nb_sym_op,neqm,ngreq,ngroup,ngroup_m,ngroup_pdb,ngroup_taux,nlat,nlatm,Noncentre,nspin,ntype,numat, &
          One_run,Orthmat,Orthmati,Orthmatt,PointGroup,PointGroup_Auto,popats,pos,posn,posn_bulk,Rmax,Rot_int,Self_nonexc, &
          Spinorbite,Rot_Atom_gr,Struct,Surface_shift,Sym_2D,Sym_4,Sym_cubic,Symmol,Taux,Taux_oc,Test_dist_min)

    if( Ylm_comp_inp .and. .not. ( Atom_comp_cal(igrpt0) .or. igrpt_nomag <= 5 ) ) then
      if( mpirank0 == 0 ) then
        call write_error
        do ipr = 6,9,3
          write(ipr,'(//A//)') ' Complex spherical harmonics not compatible with the point group !'
        end do
      endif
      stop
    endif

    Ylm_comp = Ylm_comp_inp .or. Atom_comp_cal(igrpt0)

!    if( .not. One_run .or. ( One_run .and. multi_run_e == 1 ) ) iaabsfirst = iaabs

    FDM_comp_m = FDM_comp
! FDM_comp_m = .false. ! annule l'effet FDM_comp
!                         sur traitement cas complexe sur calcul Tau_LL dans Mat
!                         sur complexe conjugue dans Tenseur_car sur integrales radiales
!         --> calcul classique FDM avec complexe conjuge
!             On garde Eneg et Singul

! nspinorb = 2 when magnetic without spin-orbit
! nspinorb = 1 when non-magnetic or spin-orbit
! nspino = 2 when spin-orbit
    if( nspin == 1 .or. Spinorbite ) then
      nspinorb = 1
    else
      nspinorb = 2
    endif

    i_radius = 1
    i_adimp = 1
    i_range = 1

!-------------------------------------------------------------------------------------------------

    natomeq = natomeq_coh
    if( Bulk_step ) then
      Rsorte = R_self_bulk
    else
      Rsorte = R_self
    endif
    adimp = adimp_e(1)

! Calculation of the number of atom in the small symmetrized cluster for the XANES calculation
    natome = natome_cal(igrpt_nomag,iopsymr,mpirank0,natomeq,natomp,pos)

! The atoms used to describe the potential can be the non-equivalent ones of the unit cell, or the non-equivalent ones
! of the cluster.
    if( multi_run_e == 1 ) then
 !   Full_atom = ( Full_atom_e  .or. ( Self_cons_c .and. .not. ( Self_nonexc .and. Proto_all ) ) &
 !                     .or. ( natome <= n_atom_proto + 1 .and. .not. Bulk_step ) &
 !                    .or. ( natome <= n_atom_proto_bulk + 1 .and. Bulk_step ) ) .and. .not. Flapw
      Full_atom = ( Full_atom_e  .or. ( Self_cons_c .and. .not. Self_nonexc ) &
                     .or. ( natome <= n_atom_proto - n_atom_proto_bulk .and. .not. Bulk_step ) &
                     .or. ( natome <= n_atom_proto_bulk .and. Bulk_step ) ) .and. .not. Flapw
      Full_atom_gen = Full_atom
    else
      Full_atom = Full_atom_gen
    endif

    if( Full_atom ) then
      n_atom_0 = 1
      n_atom_ind = natome
    else
      if( Nonexc ) then
        n_atom_0 = 1
      else
        n_atom_0 = 0
      endif
      n_atom_ind = n_atom_proto
    endif

    if( Self_cons_c ) then
      if( Self_nonexc .or. Full_atom ) then
         n_atom_0_self = 1
      else
         n_atom_0_self = 0
      endif
      n_atom_ind_self = n_atom_ind
      natome_self = natome
      natomeq_self = natomeq
      nrm_self = nrm
    else
      n_atom_0_self = 0
      n_atom_ind_self = 0
      natome_self = 0
      natomeq_self = 0
      nrm_self = 0
    endif

    if( .not. Extract ) then
      if( multi_run_e == 1 .and. Self_nonexc ) then
        allocate( chargat_init_bulk(n_atom_proto-n_atom_proto_bulk+1:n_atom_proto,nspin))
        allocate( chargat_self_bulk(n_atom_proto-n_atom_proto_bulk+1:n_atom_proto,nspin))
        allocate( pop_orb_val_bulk(n_atom_proto-n_atom_proto_bulk+1:n_atom_proto,nspin) )
        allocate( rhoato_init_bulk(0:nrm_self,nspin,n_atom_proto-n_atom_proto_bulk+1:n_atom_proto) )
        allocate( rho_chg_bulk(0:nrm_self,nspin,n_atom_proto-n_atom_proto_bulk+1:n_atom_proto))
        allocate( rho_self_bulk(0:nrm_self,nlm_pot,nspin,n_atom_proto-n_atom_proto_bulk+1:n_atom_proto))
        allocate( Vcato_bulk(0:nrm_self,nlm_pot,n_atom_proto-n_atom_proto_bulk+1:n_atom_proto) )
        allocate( Vcato_init_bulk(0:nrm_self,n_atom_proto-n_atom_proto_bulk+1:n_atom_proto) )
        chargat_init_bulk(:,:) = 0._db
        chargat_self_bulk(:,:) = 0._db
        pop_orb_val_bulk(:,:) = 0._db
        rhoato_init_bulk(:,:,:) = 0._db
        rho_chg_bulk(:,:,:) = 0._db
        rho_self_bulk(:,:,:,:) = 0._db
        Vcato_bulk(:,:,:) = 0._db
        Vcato_init_bulk(:,:) = 0._db
      endif
      if( .not. Second_SCF ) then
        allocate( Chargat_init(n_atom_0_self:n_atom_ind_self,nspin))
        allocate( Chargat_self(n_atom_0_self:n_atom_ind_self,nspin))
        allocate( rho_chg(0:nrm_self,nspin,n_atom_0_self:n_atom_ind_self) )
        allocate( rho_self(0:nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self) )
        allocate( rhoato_init(0:nrm_self,nspin,n_atom_0_self:n_atom_ind_self) )
        allocate( Vcato_init(0:nrm_self,n_atom_0_self:n_atom_ind_self) )
        allocate( drho_ex_nex(nrm,nspin) )
        allocate( dvc_ex_nex(nrm) )
        allocate( dv_ex_nex(nrm,nspin) )
        allocate( dExc_ex_nex(nrm) )
        Chargat_init(:,:) = 0._db
        Chargat_self(:,:) = 0._db
        rho_chg(:,:,:) = 0._db
        rho_self(:,:,:,:) = 0._db
        rhoato_init(:,:,:) = 0._db
        Vcato_init(:,:) = 0._db
        drho_ex_nex(:,:) = 0._db
        dExc_ex_nex(:) = 0._db
        dv_ex_nex(:,:) = 0._db
        dvc_ex_nex(:) = 0._db
      endif
      allocate( Exc_abs_i(nrm) )
      allocate( V_abs_i(nrm,nspin) )
      allocate( Vc_abs_i(nrm) )
      if( .not. Bulk_step .and. .not. Full_atom .and. Self_nonexc ) then
        do ipr = n_atom_proto-n_atom_proto_bulk+1, n_atom_proto
          Chargat_init(ipr,:) = chargat_init_bulk(ipr,:)
          Chargat_self(ipr,:) = chargat_self_bulk(ipr,:)
          rho_chg(:,:,ipr) = rho_chg_bulk(:,:,ipr)
          rho_self(:,:,:,ipr) = rho_self_bulk(:,:,:,ipr)
          rhoato_init(:,:,ipr) = rhoato_init_bulk(:,:,ipr)
          Vcato_init(:,ipr) = Vcato_init_bulk(:,ipr)
        end do
      endif
    endif

    allocate( Hubb_diag(n_atom_0_self:n_atom_ind_self) )
    Hubb_diag(:) = .true.

    if( .not. Extract_ten ) then
      allocate( V_hubb(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self) )
      V_hubb(:,:,:,:,:) = (0._db, 0._db)
    endif

!------ SCF ----------------------------------------------------------------------------------------------------------------------------------

! Self_cons_c = nself_c > 0
    if( Self_cons_c .and. .not. Second_run .and. .not. Second_SCF ) then

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        tp_SCF_1 = real(time,db)
      endif

! Preambule coherence:
      Convergence = .false.

      allocate( chargat_self_s(n_atom_0_self:n_atom_ind_self,nspin) )
      allocate( E_KS_core(n_atom_0_self:n_atom_ind_self) )
      E_KS_core(:) = 0._db
      allocate( Energ_self(n_atom_0_self:n_atom_ind_self) )
      allocate( Energ_self_s(n_atom_0_self:n_atom_ind_self) )
      allocate( Occ_hubb(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self) )
      allocate( Occ_hubb_i(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self) )
      allocate( pop_orb_val(n_atom_0_self:n_atom_ind_self,nspin) )
      allocate( rho_self_s(0:nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self) )
      allocate( rho_self_t(0:nrm_self,nlm_pot,n_atom_0_self:n_atom_ind_self) )
      allocate( V_hubb_s(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self) )
      Energ_self_s(:) = 0._db
      V_hubb_s(:,:,:,:,:) = (0._db, 0._db)

      if( .not. Bulk_step .and. .not. Full_atom ) then
        do ipr = n_atom_proto-n_atom_proto_bulk+1, n_atom_proto
          pop_orb_val(ipr,:) = pop_orb_val_bulk(ipr,:)
        end do
      endif

      if( Optic .or. tddft ) then
        Fermi_first = .true.
        Devide_Ei = .true.
      else
        Fermi_first = .false.
        Devide_Ei = .false.
      endif

      Cal_xanes = .false.
      Nonsph = .false.
      Green = Green_self
      Eneg = .true.
      Moy_loc = Green .and. .not. Moy_cluster
      Nonexc_g = Self_nonexc
      Recop = .false.

      if( Nonexc_g ) then
        ipr1 = 1
        itab = itabs_nonexc
      else
        ipr1 = 0
        itab = itabs
      endif
      WorkF = WorkF_i

      if( .not. Nonexc_g ) iaproto(iaabs) = 0
      itypep(iaabs) = itab
      Vsphere = Vsphere_cal(Chargat,Chargm_SCF,dista(natomp),Flapw,iaproto,icheck(5),Matper,n_atom_proto,natomp)

! Pour le calcul auto-coherent, il faut calculer toutes les representations
      if( icheck(27) > 1 ) then
        ich = max( icheck(8), icheck(27) )
      else
        ich = 0
      endif
      call irep_util_all(ich,iopsymr,irep_util,karact,nb_rep,ngrph,nspino,Spinorbite)
      allocate( Repres_comp(ngrph) )
      call irep_comp(Ylm_comp,Green,iopsymr,irep_util,karact,ngrph,nspino,Repres_comp)

! L'indice 0, correspond a l'agregat en entier
      allocate( Atom_axe(natome) )
      allocate( distai(natome) )
      allocate( ia_eq(nb_sym_op,natome) )
      allocate( ia_eq_inv(natomeq) )
      allocate( ia_rep(nb_sym_op,natome,natome) )
      allocate( iaprotoi(natome) )
      allocate( igroupi(natome) )
      allocate( iopsym_atom(nopsm,natome) )
      allocate( is_eq(nb_sym_op,natome) )
      allocate( itypei(natome) )
      allocate( nb_eq(natome) )
      allocate( nb_eq_2D(0:n_atom_proto) )
      nb_eq_2D(:) = 0
      allocate( nb_rpr(natome,natome) )
      allocate( nb_rep_t(nb_sym_op,natome,natome) )
      allocate( posi(3,natome) )
      allocate( posi_self(3,natome) )
      allocate( rot_atom(3,3,natome) )

      allocate( ia_eq_inv_self(natomeq_self) )

      if( icheck(27) > 1 ) then
        ich = max( icheck(7), icheck(27) )
      else
        ich = 0
      endif

      i_self = 1

      if( Sym_2D ) call Cal_nb_eq_2D(iaproto,n_atom_proto,nb_eq_2D,natomeq,natomp)

      call Atom_selec(adimp,Atom_axe,Atom_with_axe,Nonsph,Atom_occ_hubb,Axe_atom_clu, &
        dista,distai,Full_atom,Green,Hubbard,i_self,ia_eq,ia_eq_inv,ia_rep,iaabs,iaabsi,iabsorbeur,iaproto,iaprotoi,ich,igreq, &
        igroup,igroupi,igrpt_nomag,igrpt0,iopsym_atom,iopsymr,iord,ipr1,is_eq,isymeq,itype,itypei,itypep,itypepr,Magnetic, &
        m_hubb,m_hubb_e,mpirank0,natome,n_atom_0_self,n_atom_ind_self,n_atom_proto,natomeq,natomp,natomsym,nb_eq,nb_rpr, &
        nb_rep_t,nb_sym_op,neqm,ngreq,ngroup,ngroup_hubb,ngroup_m,nlat,nlatm,nspin,nspinp,ntype,numat,nx,occ_hubb_e,Overad, &
        popats,pos,posi,Proto_calculated,rmt,rot_atom,roverad,rsort,rsorte,Spinorbite,Symmol,V_hubb,V_hubbard,Ylm_comp)

      ia_eq_inv_self(:) = ia_eq_inv(:)
      posi_self(:,:) = posi(:,:)

! in the SCF part even if SCFexc, one calls it iaprabs_nonexc anyway
      if( Full_atom ) then
        iaprabs_nonexc = iaabsi
      else
        iaprabs_nonexc = iprabs_nonexc
      endif

      nbm = 0;     nbtm = 0
      npoint = 0;  npr = 0
      nsort = 0
      nsm = 0;     nstm = 0

! Calculation of the dimensions for the FDM grid of points
      call nbpoint(Adimp,Base_hexa,D_max_pot,Green,iaabsfirst,igrpt_nomag,iopsymr,iord,Moy_loc, &
                 mpirank0,natomp,npoint,npso,nvois,nx,pos,rsort)

      if( Nonsph ) then
        npoint_ns = npoint
      else
        npoint_ns = 0
      endif
      nim = npoint
      npsom = npso

      allocate( clapl(0:nvois) )
      allocate( cgrad(nvois) )
      allocate( ivois(npsom,nvois) )
      allocate( isvois(npsom,nvois) )
      allocate( numia(npsom) )
      allocate( rvol(nim) )
      allocate( xyz(4,npsom) )

      allocate( indice(npsom,3) )
      allocate( mpres(-nx:nx,-nx:nx,-nx:nx) )

! Making of the FDM grid of points. Even with Green, it is used for the calculation of the interstitial potential
      call reseau(Adimp,Base_hexa,D_max_pot,Green,iaabsfirst,icheck(9), &
           igrpt_nomag,indice,iopsymr,iord,itypei,Moy_loc,mpirank0,mpres,natome,natomp,nim, &
           npoint,npr,npso,npsom,ntype,numia,nx,pos,posi,rmt,rsort,rvol,xyz)

      if( .not. Green ) then
        call laplac(Adimp,Base_hexa,cgrad,clapl,icheck(9),igrpt_nomag,indice,iopsymr,iord,ivois,isvois, &
            mpirank0,mpres,npso,npsom,nvois,nx)
        call bordure(Green,icheck(9),iopsymr,iord,iscratch,ivois,mpirank0,natome,nbm,nbtm,nim, &
            npoint,npso,npsom,nsm,nstm,numia,nvois,posi,rvol,xyz)
      endif
      deallocate( indice, mpres )

      allocate( Chg_coeur(n_atom_0_self:n_atom_ind_self) )
      allocate( E_starta(n_atom_0_self:n_atom_ind_self) )
      allocate( Excato(nrm,n_atom_0:n_atom_ind) )
      allocate( ibord(nbtm,natome) )
      allocate( imoy(npoint) );  allocate( imoy_out(npoint) )
      allocate( isbord(nbtm,natome) )
      allocate( ispin_maj(n_atom_0_self:n_atom_ind_self) )
      allocate( isrt(nstm) )
      allocate( nbord(natome) )
      allocate( nbordf(natome) )
      allocate( poidsa(nbm,natome) )
      allocate( poidso(nsm) )
      allocate( poidsov(npoint) ); allocate( poidsov_out(npoint) )
      allocate( rho(npoint,nspin) )
      allocate( rhoato(0:nrm,nspin,n_atom_0:n_atom_ind) )
      allocate( rhons(npoint_ns) )
      allocate( rs(npoint) )
      allocate( rsato(0:nrm,n_atom_0:n_atom_ind) )
      allocate( Vcato(0:nrm,nlm_pot,n_atom_0:n_atom_ind) )
      allocate( Vh(npoint) )
      allocate( Vhns(npoint_ns) )
      allocate( Vr(npoint,nspin) )
      allocate( Vrato(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) )
      allocate( Vxc(npoint,nspin) )
      allocate( Vxcato(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) )
      allocate( Int_statedens(lla2_state,nspinp,n_atom_0:n_atom_ind) )

      if( .not. Green ) call recup_bordure(ibord,isbord,iscratch,isrt,mpinodes0,mpirank0,natome,nbord,nbordf,nbm,nbtm,nsm, &
                                           nsort,nsortf,nstm,poidsa,poidso)

! Boucle coherence
      boucle_coh: do i_self = 1,nself_c

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(1) = real(time,db)
        endif

        if( icheck(4) > 1 ) write(3,150) i_self

        E_Fermi = -5._db / Rydb
        E_cut = E_Fermi
        E_Open_val = E_Fermi
        E_Open_val_exc = E_Fermi
        Fermi = .false.
        Level_val_abs = .false.
        Level_val_exc = .false.

        Int_statedens(:,:,:) = 0._db
        rhoato(:,:,:) = 0._db
        rsato(:,:) = 0._db
        Vcato(:,:,:) = 0._db
        Vxcato(:,:,:,:) = 0._db

        if( Bulk_atom_done ) then
          do ipr = n_atom_proto-n_atom_proto_bulk+1, n_atom_proto
            Vcato(:,:,ipr) = Vcato_bulk(:,:,ipr)
          end do
        endif

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(2) = real(time,db)
          Time_rout(3) = Time_rout(3) + Time_loc(2) - Time_loc(1)
        endif

        call potsup(alfpot,Axe_atom_gr,Bulk_atom_done,Cal_xanes,Chargat,Chargat_init, &
          Chargat_self,Delta_helm,dExc_ex_nex,drho_ex_nex,dv_ex_nex,dvc_ex_nex,Exc_abs_i,Excato,Full_atom,Helm_cos,hybrid, &
          i_self,ia_eq_inv,ia_eq_inv_self,iaabs,iaproto,iaprotoi,iapot,icheck,igreq,igroup,iprabs_nonexc,ipr1,itab, &
          itypei,itypep,itypepr,Lmax_pot,lvval,Magnetic,mpirank0,n_atom_0,n_atom_0_self,n_atom_ind, &
          n_atom_ind_self,n_atom_proto,n_atom_proto_bulk,natome,natome_self,natomeq,natomeq_self,natomp,neqm, &
          ngroup_m,ngroup_nonsph, &
          nhybm,nlat,nlatm,nlm_pot,Nonexc,Nonsph,norbv,normrmt,npoint,npoint_ns,npsom,nrato,nrm,nrm_self,nspin, &
          ntype,numat,Overlap,pop_nonsph,popatm,popatv,pos,posi,posi_self,Proto_calculated,psival,Rato,rho,rho_chg, &
          rho_self,rhoato,rhoato_init,rhoit,rhons,Rmtg,Rmtimp,Rmtg0,Rmtsd,Rot_Atom_gr,Rot_int,rs,rsato,rsort, &
          SCF_c,Self_nonexc,Sym_2D,V_abs_i,V_helm_t,V_intmax,Vc_abs_i,Vcato,Vcato_init,Vh,Vhns,Vsphere,Vxc,Vxcato, &
          V0bdcFimp(1),Width_helm,xyz)

! Operations complementaires sur le potentiel
        call Potential_comp(Bulk_atom_done,Cal_xanes,distai(natome),dV0bdcF,Full_atom, &
          Green,i_range,i_self,iaabs,iapot,iaproto,iaprotoi,icheck(13),imoy,imoy_out,iopsymr,ipr1,isrt,itypepr,korigimp,Magnetic, &
          Moy_loc,mpirank0,n_atom_0,n_atom_ind,n_atom_proto,n_atom_proto_bulk,natome,natomp,nim,nlm_pot,npoint,npsom,nptmoy, &
          nptmoy_out,nrato,nrm,nsortf,nspin,nstm,ntype,numat,poidsov,poidsov_out,pos, &
          Rato,Rchimp,rhoato,Rmtg,Rmtg0,Rmtsd,rs,rsbdc,rsbdc_out,rsort,rvol,SCF_c,Sym_2D,V_intmax,V0bdcF,V0bdcF_out,V0bdcFimp, &
          V0muf,Vcato,Vh,Vhbdc,Vhbdc_init,Vhbdc_out,Vr,Vxc,Vxcato,VxcbdcF,VxcbdcF_out,xyz,Workf)


        if( i_self == 1 ) then

          if( ipr1 == 1 ) Rmtg(0) = Rmtg(iprabs_nonexc)

! Initialisation of rho_self, when new potential in potsup is calculated
! Beyon Rmtsd initialisation is at zero.
          rho_self_s(:,1,:,:) = rhoato_init(:,:,:)

! Calculation of the starting energy for the SCF
! It is calculated only at the first cycle eventhough states can be shifted along the cycles
          call Starting_energ(Bulk_step,E_start,E_starta,Full_atom,iaprotoi,icheck(27),itypepr,lcoeur,n_atom_0, &
              n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,n_atom_proto_uc,natome,ncoeur,nenerg_coh,nlm_pot, &
              nrato,nrm,numat,nspin,ntype,Pas_SCF,Proto_calculated,psi_coeur,Rato,Vcato,Vxcato, &
              Workf)

! Calculation of the reference charge for self-consistency
          call Chg_agr(Bulk_step,Chargat,Chargat_init,Chg_coeur,Chg_reference,chg_open_val,Doping,Full_atom, &
              iaprotoi,ipr_dop,iprabs_nonexc,ispin_maj,itabs,icheck(27),itypepr,mpirank0, &
              n_atom_0_self,n_atom_ind_self,n_atom_proto,n_atom_proto_bulk,natome, &
              nb_eq,nb_eq_2D,ngreq,nrato,nrm,nrm_self,nspin,ntype,numat,pop_open_val, &
              Proto_calculated,psi_open_val,Rato,rho_chg,rho_coeur,rhoato_init,Rmtsd,SCF_mag_fix,SCF_mag_free,Sym_2D)

          chargat_self_s(:,:) = Chargat_init(:,:)
! Establishment of the energy grid for self-consistency
          allocate( Energ_coh(nenerg_coh) )
          allocate( Eimag_coh(nenerg_coh) )
          call grille_coh(Eimag_coh,Energ_coh,E_start,Green,icheck(27),nenerg_coh,Pas_SCF)

        endif

        if( Fermi_first .or. .not. Devide_Ei ) then
          nenerg = nenerg_coh
          allocate( Energ(nenerg) )
          allocate( Eimag(nenerg) )
          Energ(:) = Energ_coh(:)
          Eimag(:) = Eimag_coh(:)
        else
          nenerg = n_devide * nenerg_coh
          allocate( Energ(nenerg) )
          allocate( Eimag(nenerg) )
          Energ(1) = energ_coh(1)
          do ie = 2,nenerg
            Energ(ie) = Energ(ie-1) + Pas_SCF / n_devide
          end do
          Eimag(:) = Eimag_coh(1) / n_devide
        endif

        call Check_Ecmax(Ecineticmax,Ecineticmax_out,Eclie,Eclie_out,Eneg,Energ(nenerg),mpirank,nspin,V0bdcF,V0bdcF_out,Workf)

        do iapr = n_atom_0_self,n_atom_ind_self
          if( Full_atom ) then
            ipr = iaprotoi(iapr)
          else
            ipr = iapr
            if( Bulk_atom_done .and. ipr > n_atom_proto - n_atom_proto_bulk ) cycle
          endif
          it = itypepr(ipr)
          nr = nrato(it)
          do n = 1,nr
            if( Rato(n,it) > Rmtsd(ipr) ) exit
          end do
          do ispin = 1,nspin
            do ir = 0,n
              rho_self(ir,1,ispin,iapr) = rho_chg(ir,ispin,iapr) + ( rho_coeur(ir,it) / nspin )
            end do
            rho_self(n+1:nr,1,:,iapr) = rhoato_init(n+1:nr,:,iapr)
          end do
        end do

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(3) = real(time,db)
          Time_rout(4) = Time_rout(4) + Time_loc(3) - Time_loc(2)
        endif

        if( i_self == 1 ) then

! Lmax calculation
          if( .not. Extract .and. .not. Green ) then
            call clmax(Ecineticmax_out,Rsort,Lmaxso0,Lmaxso_m,-1,.true.)
            Lmaxso_m = min( Lmaxso_m, Lmaxso_max )
          else
            Lmaxso_m = 0
          endif

          Lmax_max = 0
          Lmax_abs = 0
          do ipr = 0,n_atom_proto
            Z = numat( itypepr(ipr) )
            call clmax(Ecineticmax,Rmtg(ipr),Lmaxat0,Lmaxat(ipr),Z,lmaxfree)
            if( Z == numat_abs ) Lmax_abs = max( Lmaxat(ipr), Lmax_abs )
            Lmax_max = max(Lmax_max,Lmaxat(ipr))
          end do
          do ipr = 0,n_atom_proto
            Z = numat( itypepr(ipr) )
            if( Z == numat_abs ) Lmaxat(ipr) = Lmax_abs
          end do

          nlmmax = (Lmax_max + 1 )**2
          nlmsam = nspinp * nlmmax
          nlmomax = ( Lmaxso_m + 1 )**2
          nso1 = nspino * nlmomax

          allocate( iato(nlmsam,natome,ngrph) )
          allocate( lato(nlmsam,natome,ngrph) )
          allocate( iso(nso1,ngrph) )
          allocate( lso(nso1,ngrph) )
          allocate( mato(nlmsam,natome,ngrph) )
          allocate( mso(nso1,ngrph) )
          allocate( nlmso0(ngrph) )
          allocate( nlmsa0(natome,ngrph) )

! Determination des (l,m) de chaque representation
          call lmrep(Green,iaprotoi,iato,icheck(15),iopsym_atom,iopsymr,irep_util,iso,itypei,karact,lato, &
               Lmaxat,Lmaxso_m,lso,mato,mpirank0,mso,n_atom_proto,natome,nbordf,ngrph,nlmsa0,nlmsam,nlmso0, &
               nso1,nsortf,nspino,ntype,numat,Orthmati,posi,rot_atom)

! Calcul of the Ylm on the FDM grid points
          if( Green ) then
            nbtm_fdm = 0
          else
            nbtm_fdm = nbtm
          endif
          allocate( Ylmato(nbtm_fdm,nlmmax,natome) )
          if( .not. Green )  then
            allocate( Ylmso(nsort,nlmomax) )
            call ylmpt(iaprotoi,ibord,icheck(15),iopsymr,isrt,Lmaxat,Lmaxso_m,n_atom_proto,natome, &
               nbord,nbtm,nlmmax,nlmomax,nsort,nstm,npsom,posi,rot_atom,xyz,Ylmato,Ylmso)
          endif

          if( .not. Green ) then
            if( Spinorbite .or. Relativiste ) then
              nicm = nim
            else
              nicm = 1
            endif
            allocate( gradvr(nicm,3,nspin))
          endif
          allocate( Lmaxa(natome) )
          allocate( nlmsa(natome) )
          allocate( drho_self(nrm_self,nlm_pot,nspinp,n_atom_0_self:n_atom_ind_self,0:mpinodes-1) )
          allocate( Statedens(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1) )
          allocate( Statedens_i(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1) )

          nlm_RIXS = 0
          allocate( Ampl_dft(nlm_RIXS,nspinp,nlmomax,nspinp,0:mpinodes-1) )

        endif

! Initialisation of cluster energy
        Energ_self(:) = 0._db

        ninitlv = nbseuil
        n_vr_0 = n_atom_0
        n_vr_ind = n_atom_ind

        Green_i = Green_int

        Tau_nondiag = .false.
        if( Spinorbite .or. Full_potential ) then
          Tau_nondiag = .true.
        elseif( Hubbard ) then
          do iapr = n_atom_0_self,n_atom_ind_self
            if( Hubb_diag(iapr) ) cycle
            Tau_nondiag = .true.
            exit
          end do
        endif

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(4) = real(time,db)
          Time_rout(5) = Time_rout(5) + Time_loc(4) - Time_loc(3)
        endif

        ie0 = 1
        nenerg0 = nenerg

! Loop over energy
        nge = ( nenerg0 - ie0 ) / mpinodes + 1

        index_e = 0

        boucle_energ: do je = 1,nge

          ie = ie0 + ( je - 1 ) * mpinodes + mpirank

          if( ie <= nenerg0 ) then

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              Time_loc(1) = real(time,db)
            endif

            index_e = index_e + 1
            if( Old_zero ) then
              Enervide = Energ(ie) - Workf
            else
              Enervide = Energ(ie)
            endif

            ich = min( icheck(16), icheck(27) )
            if( ich < 2 .and. ie > 1 ) then
              Trace_l = 0
            else
              Trace_l = Trace_k
            endif

! Potentiel calculation in the excited state
            call potex(Nonsph,axyz,alfpot,dV0bdcF,Energ(ie),Enervide,Full_atom, &
              iaabsi,iapot,iaprotoi,ich,iprabs_nonexc,iprabs,itab,itypepr,Magnetic, &
              n_atom_0,n_atom_ind,n_atom_proto,n_vr_0,n_vr_ind,natome,nlm_pot,Nonexc_g,npoint,npoint_ns,npsom,nptmoy, &
              nptmoy_out,nrato,nrm,nspin,ntype,Optic,Rato,rho,rhons,Rmtg,rs,rsato,rsbdc,rsbdc_out,Trace_l,Trace_p, &
              V_intmax,Vcato,Vh,Vhbdc,Vhbdc_out,Vhns,Vr,Vxc,VxcbdcF,VxcbdcF_out,Vrato,Vxcato,V0bdc,V0bdc_out,xyz)

            do ispin = 1,nspin

              V0bdc_out(ispin) = V0bdc(ispin)

              if( Muffintin ) call mdfmuf(Nonsph,Axe_Atom_gr,axyz, &
                Full_atom,iaabs,iaproto,iaprotoi,ich,igreq,igroup,ispin,itypep,n_atom_0,n_atom_ind, &
                n_atom_proto,natome,natomp,neqm,ngroup_m,nlm_pot,npoint,npoint_ns,npsom,nrato,nrm,nspin,ntype,pos,Rato,rho, &
                rhons,Rmtg,Trace_l,Trace_p,Vhns,V0bdc(ispin),Vr,Vrato,xyz)

              if( Supermuf ) call modmuf(Full_atom,iaprotoi,itypepr,ich, &
                ispin,n_atom_0,n_atom_ind,n_atom_proto,natome,nlm_pot,nrato,nrm,nspin,ntype,Rato,Rmtg,V0bdc(ispin),Vrato)

              if( .not. Green .and. ( Spinorbite .or. Relativiste ) ) &
                call gradpot(Base_hexa,cgrad,gradvr,ich,iord,ispin,ivois,nicm,npoint,npsom,nspin,nvois,Vr)

              Ecinetic_out(ispin) = Enervide - V0bdc_out(ispin)
              Ecinetic(ispin) = Enervide - V0bdc(ispin)

              if( Eneg ) then
! Numerical problem around 0.
                em = 0.01_db / rydb
                ei0 = 0.01_db / rydb
                if( Ecinetic_out(ispin) < 0._db ) then
                  Eimag(ie) = max( Eimag(ie), ei0 )
                elseif( Ecinetic_out(ispin) < em ) then
                  eii = ei0 * ( 1 - Ecinetic_out(ispin) / em )
                  Eimag(ie) = max( Eimag(ie), eii )
                endif
              else
                Ecinetic_out(ispin) = max( Ecinetic_out(ispin), Eclie_out)
                Ecinetic(ispin) = max( Ecinetic(ispin), Eclie )
              endif

              if( ( Ecinetic(ispin) < eps10 .or. Ecinetic_out(ispin) < eps10 ) .and. .not. Eneg ) then
                if( mpirank0 == 0 ) then
                  call write_error
                  do ipr = 6,9,3
                    write(ipr,170) Ecinetic(ispin) * rydb
                    if( nptmoy_out > 0 ) write(ipr,180) Ecinetic_out(ispin) * rydb
                  end do
                endif
                stop
              endif

            end do  ! End of loop over spin

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              Time_loc(2) = real(time,db)
              Time_rout(6) = Time_rout(6) + Time_loc(2) - Time_loc(1)
            endif

            if( abs( Eimag(ie) ) > eps10 .or. Ecinetic_out(1) < eps10 .or. Ecinetic_out(nspin) < eps10 ) then
              E_comp = .true.
            else
              E_comp = .false.
            endif

            if( No_solsing ) then
              Solsing = .false.
            elseif( E_comp .and. ( Green .or. FDM_comp ) ) then
              Solsing = .true.
            else
              Solsing = Solsing_s
            endif

            Ecmax = 0._db
            Ecmax_out = 0._db
            do ispin = 1,nspin
              Ecmax_out = max( Ecmax_out, Ecinetic_out(ispin) )
              Ecmax = max( Ecmax, Ecinetic(ispin) )
            end do

            Lmax_g = 0
            Lmax_abs_a = 0
            do ipr = 0,n_atom_proto
              Z = numat( itypepr(ipr) )
              call clmax(Ecmax,Rmtg(ipr),Lmaxat0,Lmaxat(ipr),Z,lmaxfree)
              Lmaxat(ipr) = min(Lmaxat(ipr), Lmax_max)
  ! For SCF, limited to 1 for H and He
              if( Z < 3 ) Lmaxat(ipr) = min( Lmaxat(ipr), 1 )
              if( Z == numat_abs ) Lmax_abs_a = max( Lmax_abs_a, Lmaxat(ipr) )
              Lmax_g = max( Lmaxat(ipr), Lmax_g)
            end do
            do ipr = 0,n_atom_proto
              Z = numat( itypepr(ipr) )
              if( Z == numat_abs ) Lmaxat(ipr) = Lmax_abs_a
            end do
            nlmagm = ( Lmax_g + 1 )**2

            allocate( Tau_ato(nlmagm,nspinp,nlmagm,nspinp,n_atom_0:n_atom_ind) )
            if( Green ) then
              nphiato1 = 0; nphiato7 = 0
            else
              nphiato1 = nbtm
              nphiato7 = 1
              if( abs(Eimag(ie)) > eps10 .or. Spinorbite .or. Ylm_comp ) then
                nphiato7 = 2
              else
                do igrph = 1,ngrph
                  if( .not. Repres_comp(igrph) ) cycle
                  nphiato7 = 2
                  exit
                end do
              endif
            endif
            nab_coop = 0
            allocate( phiato(nphiato1,nlmagm,nspinp,nspino,natome,nphiato7) )
            allocate( Taull(nlmagm,nspinp,nlmagm,nspinp,natome) )
            allocate( Tau_coop(nlmagm,nspinp,nlmagm,nspinp,nab_coop,nspino) )
            phiato(:,:,:,:,:,:) = 0._db
            Tau_ato(:,:,:,:,:) = (0._db,0._db)
            Taull(:,:,:,:,:) = (0._db,0._db)
            Tau_coop(:,:,:,:,:,:) = (0._db,0._db)

            if( .not. Green ) then
              call clmax(Ecmax_out,Rsort,Lmaxso0,Lmaxso,-1,.true.)
              Lmaxso = min( Lmaxso, Lmaxso_max )
            endif

            ich = min(icheck(18), icheck(27))
            if( ich > 1 ) write(3,190)

! Loop over non equivalent atoms
            n1 = n_atom_0_self

            do iapr = n1,n_atom_ind

              if( Full_atom ) then
                ipr = iaprotoi( iapr )
              else
                ipr = iapr
              endif
              Z = numat( itypepr(ipr) )
              if( Z == 0 .and. .not. Green ) cycle
              if( Z == Z_nospinorbite ) then
                Spino = .false.
              else
                Spino = Spinorbite
              endif
              it = itypepr(ipr)
              Absorbeur = .false.

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
              r(1:nr) = Rato(1:nr,it)

              Lmax = Lmaxat(ipr)

              Vrato_e(1:nr,:,:) = Vrato(1:nr,:,:,iapr)

              call Sphere(Axe_atom_grn,Ecinetic,Eimag(ie),Energ(ie),Enervide,Full_atom, &
              Full_potential,Green,Hubb_a,Hubb_d,iaabsi,iapr,iaprotoi,ibord,ich,igreq,igroupi,iopsymr,Lmax,Lmax_pot,m_hubb, &
              n_atom_0,n_atom_ind,n_atom_proto,natome,nbord,nbtm,nbtm_fdm,neqm,ngroup_m,nlm_pot, &
              nlmagm,nlmmax,nphiato1,nphiato7,npsom,nr,nspin,nspino,nspinp,Z,phiato,posi,r,Relativiste,Renorm_SCF,Rmtg(ipr), &
              Rmtsd(ipr),Spino,Tau_ato,V_hubb_t,V_intmax,V0bdc,Vrato_e,xyz,Ylm_comp,Ylmato)

              deallocate( r )
              deallocate( Vrato_e, V_hubb_t )

            end do

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              Time_loc(3) = real(time,db)
              Time_rout(7) = Time_rout(7) + Time_loc(3) - Time_loc(2)
            endif

            Lmax_g = 0
            do ia = 1,natome
              Lmaxa(ia) = Lmaxat( iaprotoi(ia) )
              Lmax_g = max(Lmax_g,Lmaxa(ia))
            end do

            Time_fill = 0._db; Time_tria = 0._db

            if( .not. ( Second_run .or. Second_SCF ) ) then

              do ispin = 1,nspinorb
! igrph: index of representations
                do igrph = 1,ngrph

                  nlmsamax = 0
                  nlmsa(:) = 0
                  do ia = 1,natome
                    Z = numat(itypei(ia))
                    if( .not. Green .and. Z == 0 ) cycle
                    n = nlmsa0(ia,igrph)
                    allocate( Lval(n) )
                    Lval(1:n) = lato(1:n,ia,igrph)
                    call cnlmmax(Lmaxa(ia),Lval,n,nlmsa(ia))
                    nlmsamax = max(nlmsamax,nlmsa(ia))
                    deallocate( Lval )
                  end do
                  if( nlmsamax == 0 .and. igrph /= ngrph ) cycle
                  ich = min(icheck(19), icheck(27))

                  if( Green ) then
                    call msm(Axe_Atom_grn,Cal_xanes,Classic_irreg,Dist_coop,Ecinetic,Eimag(ie),Full_atom, &
                      ia_coop,ia_eq,ia_rep, &
                      iaabsi,iaprotoi,iato,ich,igroupi,igrph,iopsymr,irep_util,is_eq,ispin,karact,lato,Lmaxa,Lmax_g,mato, &
                      n_atom_0,n_atom_coop,n_atom_ind,n_atom_proto,nab_coop,natome,natomp,nb_eq,nb_rpr,nb_rep_t,nb_sym_op,nchemin, &
                      ngroup_m,ngrph,nlmagm,nlmsa,nlmsam,nlmsamax,Normaltau,nspin,nspino,nspinp,pos,posi,Recop,Repres_comp(igrph), &
                      Rmtg,rot_atom,Solsing,Spinorbite,State_all,Sym_cubic,Tau_ato,Tau_nondiag,Tau_coop,Taull,Time_fill, &
                      Time_tria,Ylm_comp)
                  else
                    n = nlmso0(igrph)
                    allocate( Lval(n) )
                    Lval(1:n) = lso(1:n,igrph)
                    call cnlmmax(Lmaxso,Lval,n,nlmso)
                    deallocate( Lval )
                    call mat(Adimp,Ampl_dft,Atom_axe,Axe_Atom_grn,Base_hexa,Basereel,Cal_xanes,cgrad, &
                      clapl,Classic_irreg,Dist_coop,distai,E_comp,Ecinetic_out,Eclie_out,Eimag(ie),Eneg,Enervide, &
                      FDM_comp_m,Full_atom,gradvr,ia_coop,iaabsi,iaprotoi,iato,ibord,ich,ie,igroupi,igrph,irep_util, &
                      isbord,iso,ispin,isrt,ivois,isvois,karact,lato,Lmaxa,Lmaxso,lso,mato,mpinodes,mpirank0,mso,nab_coop, &
                      natome,n_atom_0,n_atom_coop,n_atom_ind,n_atom_proto,nbm,nbord,nbordf,nbtm,ngroup_m,ngrph,nim, &
                      nicm,nlm_rixs,nlmagm,nlmmax,nlmomax,nlmsa,nlmsam,nlmso,nphiato1,nphiato7,npoint,npr,npsom,nsm, &
                      nsort,nsortf,nso1,nspin,nspino,nspinp,nstm,numia,nvois,phiato,poidsa,poidso,Posi,R_rydb,Recop,Relativiste, &
                      Repres_comp(igrph),RIXS,Rmtg,Rsort,rvol,Rydberg,Solsing,Spinorbite,State_all,Sym_cubic,Tau_ato,Tau_coop, &
                      Taull,Time_fill,Time_tria,V0bdc_out,Vr,xyz,Ylm_comp,Ylmato,Ylmso)
                  endif

                end do  ! end of loop over representations
              end do  ! end of loop over spin

            endif

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              Time_loc(4) = real(time,db)
              Time_rout(8) = Time_rout(8) + Time_loc(4) - Time_loc(3)
              Time_rout(9) = Time_rout(9) + Time_fill ! Filling
              Time_rout(10) = Time_rout(10) + Time_tria ! Triang
            endif

            deallocate( Tau_ato )
            deallocate( phiato )

            ich = icheck(27)
            Solsing_o = .false.

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              Time_loc(5) = real(time,db)
              Time_rout(11) = Time_rout(11) + Time_loc(5) - Time_loc(4)
            endif

            call Cal_dens(Cal_xanes,Classic_irreg,Density_comp,drho_self,Ecinetic,Eimag(ie),Energ(ie),Enervide,Full_atom, &
              Full_potential,.false.,Hubb,Hubb_diag,Hubb_diag_abs,iaabsi,iaprotoi,ich,itypei,itypepr,lla2_state, &
              Lmax_pot,Lmaxat,m_hubb,mpinodes,mpirank,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,natome, &
              nlmagm,nlm_pot, &
              nrato,nrm,nrm_self,nspin,nspino,nspinp,ntype,numat,Rato,Relativiste,Renorm_SCF,Rmtg,Rmtsd,Solsing,Solsing_o, &
              Spinorbite,State_all_out,Statedens,Statedens_i,Taull,V_hubb,V_hubb_abs,V_intmax,V0bdc,Vrato,Ylm_comp)

            deallocate( Taull )
            deallocate( Tau_coop )

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              Time_loc(6) = real(time,db)
              Time_rout(12) = Time_rout(12) + Time_loc(6) - Time_loc(5)
            endif

          endif ! end of condition on ie <= nenerg, for parallel processing
!if( ie == 1 ) stop
          if( mpinodes > 1 ) then

            call MPI_RECV_statedens(lla2_state,Lmaxat,mpinodes,mpirank,mpirank0,n_atom_0, &
                                      n_atom_ind,n_atom_proto,nspinp,Statedens,Statedens_i)

            call MPI_RECV_self(drho_self,nlm_pot,mpinodes,mpirank,mpirank0,n_atom_0_self, &
                                 n_atom_ind_self,nrm_self,nspinp)
          endif

          do ie_computer = 0,mpinodes-1

            ie = ( je - 1 ) * mpinodes + ie_computer + ie0

            if( ie > nenerg0 ) exit

            if( mpirank0 == 0 ) then

              call CPU_TIME(time)
              Time_loc(1) = real(time,db)

              call CPU_TIME(time)
              Time_loc(2) = real(time,db)
              Time_rout(13) = Time_rout(13) + Time_loc(2) - Time_loc(1)

              ich = icheck(27)

              call Search_Fermi(Bulk_atom_done,Bulk_step,Chg_reference,chg_open_val,Chargat_self,Density,Doping,drho_self,E_cut, &
                E_Open_val,E_Open_val_exc,E_starta,Energ,E_Fermi,Energ_self,Fermi,Full_atom,Hubb,Hubb_diag,iaabsi,iaprotoi, &
                i_self,ich,ie,ie_computer,Int_statedens,ipr_dop,ispin_maj,itypei,itypepr,lamstdens, &
                Level_val_abs,Level_val_exc,lla_state,lla2_state,Lmaxat,m_hubb,mpinodes,n_atom_0,n_atom_0_self, &
                n_atom_ind,n_atom_ind_self,n_atom_proto,n_atom_proto_bulk,n_atom_proto_uc,natome,nb_eq,nb_eq_2D,nenerg,ngreq, &
                nlm_pot,nomfich_s,nrato, &
                nrm,nrm_self,nspin,nspinp,ntype,numat,Occ_hubb,Occ_hubb_i,pop_orb_val,Proto_calculated,Rato,rho_self,rho_self_t, &
                Rmtsd,SCF_elecabs,SCF_mag_fix,Self_nonexc,Statedens,Statedens_i,Sym_2D,V_hubb,V_hubbard,Ylm_comp)

              call CPU_TIME(time)
              Time_loc(3) = real(time,db)
              Time_rout(14) = Time_rout(14) + Time_loc(3) - Time_loc(2)

            endif

              call MPI_Bcast(Fermi,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
              call MPI_Bcast(E_cut,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
              call MPI_Bcast(E_Fermi,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)

! When reaching Fermi one exits the loop
            if( Fermi ) exit boucle_energ

          end do

        end do boucle_energ   ! End of loop over energy


        if( Hubbard .and. mpinodes0 > 1 ) call MPI_Bcast_Hubb(Hubb_diag,m_hubb,mpirank0, &
                                                              n_atom_0_self,n_atom_ind_self,nspinp,V_hubb)

! Cluster energy calculation

        if( mpirank0 == 0 ) then
          call Energ_KS_core(Bulk_step,E_KS_core,Full_atom,iaprotoi,icheck(27),itypepr,Jseuil,lcoeur,Lseuil,n_atom_0, &
              n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,n_atom_proto_uc,natome,ncoeur,nlm_pot,nrato, &
              nrm,nspin,nseuil,ntype,numat,Proto_calculated,psi_coeur,Rato,Vcato,Vxcato)

          call Energ_DFT(Bulk_step,Doping,En_cluster,Energ_self,E_KS_core,Excato,Full_atom,Hubb,iaprotoi,icheck(27),ipr_dop, &
            itypepr,m_hubb,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,n_atom_proto_uc,natome,nb_eq,nb_eq_2D, &
            ngreq,nlm_pot,nrm,nrm_self,nrato,nspin,nspinp,ntype,numat,occ_hubb,Proto_calculated,Rato,rho_self,rmtsd,Sym_2D, &
            V_hubb,Vcato,Vxcato)

        end if

! Exchange of tables for SCF
        if( mpinodes0 > 1 ) then
          m = n_atom_ind_self - n_atom_0_self + 1
          m2 = nspin * m
          n = ( 1 + nrm_self ) * nspin * m * nlm_pot

            call MPI_Bcast(rho_self,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
            call MPI_Bcast(Chargat_self,m2,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
            call MPI_Bcast(Energ_self,m,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
            call MPI_Bcast(En_cluster,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
        endif

! Preparation of next iteration for self-consistence
        if( Fermi_first .and. (TDDFT .or. Optic) ) then
          Delta_E_conv = 2 * Delta_En_conv
        else
          Delta_E_conv = Delta_En_conv
        endif

        call prep_next_iter(Bulk_step,Chargat_self,chargat_self_s,Convergence,Delta_E_conv,Delta_energ,Delta_energ_s, &
             Doping,En_cluster,En_cluster_s,Energ_self,Energ_self_s,Fermi,Fermi_first,Full_atom,Hubb, &
             Hubbard,i_self,iaprotoi,icheck(27),ipr_dop,itypepr,Lmaxat,m_hubb,mpirank0,n_atom_0_self,n_atom_ind_self, &
             n_atom_proto,n_atom_proto_uc,n_devide,natome,nb_eq,nb_eq_2D,ngreq,nlm_pot,nrm_self,nself_c,nspin,nspinp, &
             ntype,numat,Occ_hubb,Occ_hubb_i,p_self,p_self_max,p_self0_c,pop_orb_val,Proto_calculated,rho_self,rho_self_s, &
             Rmtsd,Sym_2D,V_hubb,V_hubb_s,V_hubbard)

        deallocate( Energ, Eimag )

        if( Convergence ) exit

      end do boucle_coh  ! End of loop over iteration to converge in SCF

      if( Bulk_step .and. .not. Full_atom .and. Self_nonexc ) then
   ! One keeps the bulk atom occupancies
        pop_orb_val_bulk(iprabs_nonexc,:) = pop_orb_val(iprabs_nonexc,:)
        do ipr = n_atom_proto - n_atom_proto_bulk + 1,n_atom_proto
          if( numat( itypepr(ipr) ) == numat_abs ) cycle
          pop_orb_val_bulk(ipr,:) = pop_orb_val(ipr,:)
        end do
      endif

      deallocate( Ampl_dft, Chg_coeur, chargat_self_s, clapl, cgrad, drho_self )
      deallocate( iato, ibord, imoy, imoy_out, isbord, ispin_maj, isrt, Int_statedens, iso, ivois, isvois )
      deallocate( lato, Lmaxa, lso )
      deallocate( E_KS_core, E_starta, Energ_self, Energ_self_s, Excato )
      deallocate( mato, mso )
      deallocate( nbord, nbordf, nlmsa, nlmsa0, nlmso0, numia )
      deallocate( Occ_hubb, Occ_hubb_i, poidsa, poidso, poidsov, poidsov_out, pop_orb_val )
      deallocate( Repres_comp, rho, rhoato, rhons, rho_self_s, rho_self_t, rs, rsato, rvol )
      deallocate( Statedens, Statedens_i, V_hubb_s, Vcato, Vh, Vhns, Vr, Vrato, Vxc, Vxcato, xyz, Ylmato )
      if( .not. Green ) deallocate( gradvr, Ylmso )

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        tp_SCF_2 = real(time,db)
        Time_rout(18) = Time_rout(18) + tp_SCF_2 - tp_SCF_1
      endif

    elseif( .not. ( Extract .or. Second_run .or. Second_SCF ) ) then

      if( Flapw ) then
        E_Fermi = 0._db
      else
        E_Fermi = -5._db / Rydb
      endif
      E_cut = E_Fermi

    endif

!------ XANES----------------------------------------------------------------------------------------------------------------------------------

    if( mpirank0 == 0 ) then
      call CPU_TIME(time)
      tp_XANES_1 = real(time,db)
      Time_loc(1) = tp_XANES_1
    endif

    Cal_xanes = .true.
    if( Bulk_step ) then
      Green = Green_bulk
    else
      Green = Green_s
    endif
    Moy_loc = Green .and. .not. Moy_cluster

    Lmax_tddft = Lmax_tddft_inp
    if( Optic ) then
      if( Lmax_tddft_inp >= 0 ) then
        Lmax_probe = Lmax_tddft_inp
      else
        if( numat_abs < 21 ) then
          Lmax_probe = 2
        else
          Lmax_probe = 3
        endif
        if( Lmaxat0 > - 1 ) Lmax_probe = min( Lmax_probe, Lmaxat0 )
      endif
    else
      if( Octupole ) then
        Lmax_probe = Lseuil + 3
      elseif( Quadrupole ) then
        Lmax_probe = Lseuil + 2
      elseif( Multipole(1) ) then
        Lmax_probe = Lseuil + 1
      else
        Lmax_probe = Lseuil
      endif
    endif
    Lmax_tddft = max( Lmax_tddft, Lmax_probe )
    nlm_probe = (Lmax_probe + 1)**2
    nlms_pr = nlm_probe * nspino

    if( .not. Self_cons_c ) i_self = 1
    Nonsph = .not. Green .and. Atom_nonsph

    if( Old_zero .or. Extract_ten ) then
      WorkF = WorkF_i
    else
 ! For XANES cycle, zero energy is set at Fermi energy
 ! E_cut is different from E_Fermi when cutting is set on the absorbing atom occupancy l orbital.
      if( .not. ( Second_run .or. Second_SCF .or. Extract ) ) E_cut = E_cut - E_Fermi
      WorkF = WorkF_i - E_Fermi
    endif

    if( Tddft .or. Optic ) then
      Eneg = .false.
    elseif( Green .or. FDM_comp ) then
      Eneg = .not. Eneg_n_i
    else
      Eneg = Eneg_i
    endif

    if( icheck(4) > 0 ) write(3,140)

    if( icheck(27) > 0 ) write(3,160) E_cut * Rydb, E_Fermi * Rydb

    Nonexc_g = Nonexc
    if( Nonexc_g ) then
      ipr1 = 1
      itab = itabs_nonexc
    else
      ipr1 = 0
      itab = itabs
    endif
    if( .not. nonexc_g ) iaproto(iaabs) = 0
    itypep(iaabs) = itab

    Vsphere = Vsphere_cal(Chargat,Chargm,dista(natomp),Flapw,iaproto,icheck(5),Matper,n_atom_proto,natomp)

! Evaluation of the tensor shape of the cartesian tensors
    call Tensor_shape(Atom_with_axe,Nonsph,Axe_atom_clu,Dipmag, &
           E1E2e,Green,iaabs,icheck(6),igroup,igrpt0,iopsymc,iopsymr,itype,itypep,ldip,loct,lqua,Lseuil,Magnetic,mpirank0, &
           msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole,n_oo,natomp,ngroup, &
           ngroup_m,nlat,nlatm,nspin,ntype,numat,Octupole,popats,pos,Quadrupole,rot_atom_abs,Rot_int, &
           Spinorbite,State_all,Symmol,Tensor_imp)
    E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
    E1M1 = Multipole(4); E2E2 = Multipole(6);
    E3E3 = Multipole(7); M1M1 = Multipole(8)

    allocate( Bragg_abs(natomsym,npldafs) )
    allocate( Bragg_rgh_bulk_abs(natomsym,npldafs) )
    allocate( poldafse(3,npldafs,nphim) )
    allocate( poldafss(3,npldafs,nphim) )
    allocate( Vecdafse(3,npldafs,nphim) )
    allocate( Vecdafss(3,npldafs,nphim) )

    ngrm = 0
    do ipr = 1,n_atom_proto
      ngrm = max( ngrm, ngreq(ipr) )
    end do

    allocate( ltypcal(nplrm) )
    allocate( nomabs(ncolm) )
    allocate( pol(3,nplrm) )
    allocate( vec(3,nplrm) )
    allocate( pdp(nplrm,2) )

! Evaluation of xanes polarisazions and wave vectors
    if( multi_run == 1 .or. Bulk_step ) then
      ich = icheck(6)
    else
      ich = 0
    endif
    call polond(axyz,Dipmag,ich,ltypcal,Moyenne,mpirank0,msymdd,msymqq,n_mat_polar,ncolm,ncolr,nomabs,nple,nplr, &
        nplrm,nxanout,Octupole,Orthmatt,pdp,pdpolar,pol,polar,Polarise,Quadrupole,Veconde,Vec,Xan_atom)

    ncolt = ncolr

    if( .not. Matper ) then
      Volume_maille = 0._db
    elseif( Bulk_step ) then
      Volume_maille = Cal_Volume_maille(axyz_bulk,angxyz_bulk)
    else
      Volume_maille = Cal_Volume_maille(axyz,angxyz)
    endif

    allocate( Taux_ipr(n_atom_proto) )
    do ipr = 1,n_atom_proto
      if( Taux ) then
        Taux_ipr(ipr) = sum( Taux_oc(abs( igreq(ipr,1:ngreq(ipr)) )) ) / ngreq(ipr)
      else
        Taux_ipr(ipr) = 1._db
      endif
    end do
    f_avantseuil = f_cal(Bulk_step,Doping,Eseuil,icheck(6),itypepr,n_atom_proto,n_atom_proto_uc,nbseuil, &
                         ngreq,ntype,numat,numat_abs,Taux_ipr,Volume_maille)
    deallocate( Taux_ipr )

    n_bulk_z = 0; n_bulk_z_max = 0
    n_bulk_z_abs = 0; n_bulk_z_max_abs = 0
    if( Dafs .and. Bulk_step ) then
      do i = 1,3
        if( i == 2 ) allocate( n_bulk_zc(n_bulk_z) )
        if( i == 2 ) allocate( n_bulk_zc_abs(n_bulk_z_abs) )
        if( i == 3 ) allocate( igr_bulk_z(n_bulk_z_max,n_bulk_z) )
        if( i == 3 ) allocate( igr_bulk_z_abs(n_bulk_z_max_abs,n_bulk_z_abs) )
        call cal_n_bulk_z(First_bulk_step,i,igr_bulk_z,igr_bulk_z_abs,igreq,iprabs_nonexc,itype,n_atom_bulk,n_atom_proto, &
                        n_atom_proto_bulk,n_bulk_z,n_bulk_z_abs,n_bulk_zc,n_bulk_zc_abs,n_bulk_z_max,n_bulk_z_max_abs, &
                        neqm,ngreq,ngroup,ntype,numat,posn_bulk)
      end do
    else
      allocate( n_bulk_zc(n_bulk_z) )
      allocate( n_bulk_zc_abs(n_bulk_z_abs) )
      allocate( igr_bulk_z(n_bulk_z_max,n_bulk_z) )
      allocate( igr_bulk_z_abs(n_bulk_z_max_abs,n_bulk_z_abs) )
    endif
    n_max = max(1,n_bulk_z_abs)
    allocate( phdt(npldafs,nphim,n_max) )
    allocate( Sum_Bragg_bulk_abs_f0(npldafs,n_max) )
    allocate( Sum_Bragg_nonabs_f(npldafs,n_bulk_z) )
    allocate( Length_rel(n_bulk_z) )
    allocate( Length_rel_abs(n_bulk_z_abs) )

    if( Dafs ) then

      call Prepdafs(Abs_in_bulk,Angle_or,Angle_mode,Angpoldafs,Angxyz,Angxyz_bulk,Angxyz_cap,Angxyz_int,Angxyz_sur,Axe_atom_gr, &
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
          popatm,posn,posn_bulk,posn_cap,psival,Rato,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_shift,Taux,Taux_cap,Taux_oc, &
          Temp,Temp_coef,Temp_B_iso,Vec_orig,Vecdafse,Vecdafsem,Vecdafss,Vecdafssm,Z_bulk,Z_cap)

      call Col_dafs_name(angpoldafs,Bormann,Full_self_abs,hkl_dafs,isigpi,mpirank0,ncolm,ncolr,ncolt, &
                           nomabs,npldafs,Operation_mode_used,Self_abs)
    endif

    n_tens_max = 0
    if( mpirank0 == 0 .and. Spherical_tensor ) then
      if( E1E1 ) n_tens_max = n_tens_max + n_tens_dd
      if( E1E2 ) n_tens_max = n_tens_max + n_tens_dq
      if( E1E2 .and. ( nspin == 2 .or. Spinorbite ) ) n_tens_max = n_tens_max + n_tens_dq
      if( E2E2 ) n_tens_max = n_tens_max + n_tens_qq
      if( E1M1 ) n_tens_max = n_tens_max + n_tens_dm
      if( E1M1 .and. ( nspin == 2 .or. Spinorbite ) ) n_tens_max = n_tens_max + n_tens_dm
    endif
    if( mpirank0 == 0 ) allocate( Int_tens(n_tens_max*ninitr,0:natomsym) )

! When true, only the cartesian tensors are extracted
    if( Extract_ten ) then

      deallocate( karact )

      if( Tddft ) then
        Final_tddft = .true.
        deallocate( Epsii )
        ninitr = 1
        allocate( Epsii(ninitr) )
        Epsii(1) = Epsii_moy
      elseif( Core_resolved ) then
        if( Optic ) then
       ! In optic, the is one output per "l" plus the total
          if( numat_abs <= 4 ) then
            ninitr = 1
          elseif( numat_abs > 4 .and. numat_abs <= 20 ) then
            ninitr = 2
          elseif( numat_abs > 21 .and. numat_abs <= 56 ) then
            ninitr = 3
          else
            ninitr = 4
          endif
          ninitr = ninitr + 1
        else
          ninitr = ninit
        endif
      else
        ninitr = nbseuil
      endif

      call Extract_write_coabs(Abs_in_Bulk_roughness,Abs_U_iso,Allsite,Ang_rotsup,angxyz,axyz,Bragg_abs,Bragg_rgh_bulk_abs, &
         Bulk_step,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio,E_cut,Energ_s,Energphot,Epsii,Eseuil,Final_tddft,f_avantseuil, &
         Full_self_abs,hkl_dafs,i_range,iabsorig(multi_run),icheck(21),igr_bulk_z_abs,Int_tens,isigpi,isymeq, &
         nsymextract(multi_run),jseuil,Length_abs,Length_rel,Length_rel_abs,ltypcal,Matper,Moyenne,mpinodee0,multi_0,multi_imp, &
         Multipole, &
         n_abs_rgh,n_bulk_sup,n_bulk_z,n_bulk_z_abs,n_bulk_z_max_abs,n_bulk_zc_abs,n_max,n_multi_run,n_oo, &
         n_rel,n_tens_max,natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ninit1,ninitr,nomabs,nomfich,nomfich_cal_conv, &
         nom_fich_extract,nomfich_s,nphi_dafs,nphim,npldafs,nplr,nplrm,nseuil,nspinp,numat_abs,nxanout,pdp, &
         phdf0t,phdf0t_rgh_bulk,phdt,phdt_rgh_Bulk,pol,poldafse,poldafss,Self_abs,Spherical_signal,Spherical_tensor,Spinorbite, &
         Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_ref,Taux_eq,Tddft,V0muf,Vecdafse,Vecdafss,Vec,Volume_maille,Xan_atom)

      deallocate( Energ_s, Eimag_s )
      deallocate( Axe_atom_clu, dista, iaproto )
      deallocate( igroup, isymeq, itypep, pos )
      deallocate( Taux_eq, Int_tens, poldafse, poldafss )
      deallocate( Bragg_abs, Bragg_rgh_bulk_abs, phdt, Sum_Bragg_bulk_abs_f0, Sum_Bragg_nonabs_f, Vecdafse, Vecdafss )
      deallocate( ltypcal, nomabs, pdp, pol, vec )

      multi_0 = multi_0 + max( 1, n_Bulk_z_abs )
      if( Bulk_roughness > 0 .and. Bulk_step ) multi_0 = multi_0 + 1

      cycle boucle_multi

    endif

    Tddft_xanes = Tddft .and. .not. Optic
    Optic_xanes = Optic

! Calculation of useful representations, that is probed by the photoelectron
    if( .not. ( Density .or. Optic .or. NRIXS .or. COOP .or. RIXS ) ) then
      call Etafin(E1M1,icheck(8),iopsymr,irep_util,jseuil,karact,ldip,lmoins1,loct,lplus1,lqua,Lseuil, &
                  M1M1,mpirank0,nb_rep,nbseuil,ngrph,nspino,Spinorbite,State_all,Sym_cubic,Ylm_comp)
      Recop = ( Sym_cubic .or. Sym_4 ) .and. .not. State_all
    else
      if( icheck(27) > 1 ) then
        ich = max( icheck(8), icheck(27) )
      else
        ich = 0
      endif
      call irep_util_all(ich,iopsymr,irep_util,karact,nb_rep,ngrph,nspino,Spinorbite)
      Recop = .false.
    endif
    allocate( Repres_comp(ngrph) )
    call irep_comp(Ylm_comp,Green,iopsymr,irep_util,karact,ngrph,nspino,Repres_comp)

    if( Self_cons_c .and. .not. ( Second_run .or. Second_SCF .or. Extract ) ) then
      deallocate( Atom_axe )
      deallocate( distai )
      deallocate( ia_eq, ia_eq_inv, ia_rep, iaprotoi, igroupi, iopsym_atom, is_eq, itypei )
      deallocate( nb_eq, nb_eq_2D, nb_rpr, nb_rep_t )
      posi_self(:,:) = posi(:,:)
      deallocate( posi )
      deallocate( rot_atom )
    endif

! Loop over the cluster radius - interpoint distances values
    boucle_i_range: do i_range = 1,n_range

      natomeq = natomeq_s(i_radius)
      Rsorte = Rsorte_s(i_radius)
      adimp = adimp_e(i_adimp)

! Calculation of the non-equivalent atoms inside the cluster of calculation
      natome = natome_cal(igrpt_nomag,iopsymr,mpirank0,natomeq,natomp,pos)

! The atoms used to describe the potential can be only the non-equivalent atoms of the unit cell, or the non-equivalent atoms
! inside the cluster
      if( i_range == 1 .and. multi_run_e == 1 .and. .not. self_cons_c ) then
        Full_atom = ( Full_atom_e  &
                     .or. ( natome <= n_atom_proto -  n_atom_proto_bulk .and. .not. Bulk_step ) &
                     .or. ( natome <= n_atom_proto_bulk .and. Bulk_step ) ) .and. .not. Flapw
        Full_atom_gen = Full_atom
      else
        Full_atom = Full_atom_gen
      endif

      if( Full_atom ) then
        n_atom_0 = 1
        n_atom_ind = natome
      else
        if( Nonexc ) then
          n_atom_0 = 1
        else
          n_atom_0 = 0
        endif
        n_atom_ind = n_atom_proto
      endif

      if( mpirank0 == 0 .and. i_range == 1 ) then
        allocate( Int_statedens(lla2_state,nspinp,n_atom_0:n_atom_ind) )
        Int_statedens(:,:,:) = 0._db
      endif

      allocate( Atom_axe(natome) )
      allocate( distai(natome) )
      allocate( ia_eq(nb_sym_op,natome) )
      allocate( ia_eq_inv(natomeq) )
      allocate( ia_rep(nb_sym_op,natome,natome) )
      allocate( iaprotoi(natome) )
      allocate( igroupi(natome) )
      allocate( iopsym_atom(nopsm,natome) )
      allocate( is_eq(nb_sym_op,natome) )
      allocate( itypei(natome) )
      allocate( nb_eq(natome) )
      allocate( nb_rpr(natome,natome) )
      allocate( nb_rep_t(nb_sym_op,natome,natome) )
      allocate( posi(3,natome) )
      allocate( rot_atom(3,3,natome) )

      Rsorte = Rsorte_s(i_radius)

      ich = icheck(7)

      call Atom_selec(adimp,Atom_axe,Atom_with_axe,Nonsph,Atom_occ_hubb,Axe_atom_clu, &
        dista,distai,Full_atom,Green,Hubbard,i_self,ia_eq,ia_eq_inv,ia_rep,iaabs,iaabsi,iabsorbeur,iaproto,iaprotoi,ich,igreq, &
        igroup,igroupi,igrpt_nomag,igrpt0,iopsym_atom,iopsymr,iord,ipr1,is_eq,isymeq,itype,itypei,itypep,itypepr,Magnetic, &
        m_hubb,m_hubb_e,mpirank0,natome,n_atom_0_self,n_atom_ind_self,n_atom_proto,natomeq,natomp,natomsym,nb_eq,nb_rpr, &
        nb_rep_t,nb_sym_op,neqm,ngreq,ngroup,ngroup_hubb,ngroup_m,nlat,nlatm,nspin,nspinp,ntype,numat,nx,occ_hubb_e,Overad, &
        popats,pos,posi,Proto_calculated,rmt,rot_atom,roverad,Rsort,Rsorte,Spinorbite,Symmol,V_hubb,V_hubbard,Ylm_comp)

      Hubb_abs = Hubb(itabs)
      if( Hubb_abs .and. i_range == 1 .and. .not. extract ) then
        if( Full_atom ) then
          iapr = iaabsi
        elseif( Self_nonexc ) then
          iapr = iprabs_nonexc
        else
          iapr = iprabs
        endif
        V_hubb_abs(:,:,:,:) = V_hubb(:,:,:,:,iapr)
        Hubb_diag_abs = Hubb_diag(iapr)
      elseif( .not. i_range > 1 .and. .not. extract ) then  ! read in the bav file
        Hubb_diag_abs = .true.
      endif

      if( .not. ( Self_cons_c .or. Second_run .or. Extract ) .and. i_range == 1 ) then
        allocate( ia_eq_inv_self(natomeq_self) )
        ia_eq_inv_self(:) = 0
      endif

      nbm = 0;    nbtm = 0
      npoint = 0; npr = 0
      npoint_ns = 0
      nsort = 0
      nsm = 0;    nstm = 0
      nim = 0
      npsom = 0

      if( .not. ( second_run .or. Extract ) ) then

! Calculation of the table dimensions for the FDM grid of points
        call nbpoint(Adimp,Base_hexa,D_max_pot,Green,iaabsfirst,igrpt_nomag,iopsymr,iord,Moy_loc, &
                   mpirank0,natomp,npoint,npso,nvois,nx,pos,rsort)

        if( Nonsph ) npoint_ns = npoint
        nim = npoint
        npsom = npso

        if( i_range > 1 ) deallocate( clapl, cgrad, ivois, isvois, numia, rvol, xyz )

        allocate( clapl(0:nvois) )
        allocate( cgrad(nvois) )
        allocate( ivois(npsom,nvois) )
        allocate( isvois(npsom,nvois) )
        allocate( numia(npsom) )
        allocate( rvol(nim) )
        allocate( xyz(4,npsom) )
        allocate( indice(npsom,3) )
        allocate( mpres(-nx:nx,-nx:nx,-nx:nx) )

! Elaboration of the FDM grid of point. Even in Green, one defines a grid of point to calculate the average interstitial potential
        call reseau(Adimp,Base_hexa,D_max_pot,Green,iaabsfirst,icheck(9), &
             igrpt_nomag,indice,iopsymr,iord,itypei,Moy_loc,mpirank0,mpres,natome,natomp,nim, &
             npoint,npr,npso,npsom,ntype,numia,nx,pos,posi,rmt,rsort,rvol,xyz)

        if( .not. Green ) then
          call laplac(Adimp,Base_hexa,cgrad,clapl,icheck(9),igrpt_nomag,indice,iopsymr,iord,ivois,isvois, &
              mpirank0,mpres,npso,npsom,nvois,nx)
          call bordure(Green,icheck(9),iopsymr,iord,iscratch,ivois,mpirank0,natome,nbm,nbtm,nim, &
              npoint,npso,npsom,nsm,nstm,numia,nvois,posi,rvol,xyz)
        endif
        deallocate( indice, mpres )

        allocate( Excato(nrm,n_atom_0:n_atom_ind) )
        allocate( rs(npoint) )
        allocate( rhoato(0:nrm,nspin,n_atom_0:n_atom_ind) )
        allocate( rsato(0:nrm,n_atom_0:n_atom_ind) )
        allocate( Vcato(0:nrm,nlm_pot,n_atom_0:n_atom_ind) )
        allocate( Vh(npoint) )
        allocate( Vhns(npoint_ns) )
        allocate( Vxc(npoint,nspin) )
        allocate( Vxcato(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) )
        Excato(:,:) = 0._db
        rsato(:,:) = 0._db
        rhoato(:,:,:) = 0._db
        Vcato(:,:,:) = 0._db
        Vxcato(:,:,:,:) = 0._db
        if( Bulk_atom_done ) then
          do ipr = n_atom_proto-n_atom_proto_bulk+1, n_atom_proto
            Vcato(:,:,ipr) = Vcato_bulk(:,:,ipr)
          end do
        endif

      endif

      allocate( ibord(nbtm,natome) )
      allocate( isbord(nbtm,natome) )
      allocate( imoy(npoint) );  allocate( imoy_out(npoint) )
      allocate( isrt(nstm) )
      allocate( nbord(natome) )
      allocate( nbordf(natome) )
      allocate( poidsa(nbm,natome) )
      allocate( poidso(nsm) )
      allocate( poidsov(npoint) ); allocate( poidsov_out(npoint) )
      allocate( rho(npoint,nspin) )
      allocate( rhons(npoint_ns) )
      allocate( Vr(npoint,nspin) )
      allocate( Vrato(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) )

      if( .not. ( Green .or. Second_run .or. Extract ) ) &
        call recup_bordure(ibord,isbord,iscratch,isrt,mpinodes0,mpirank0,natome,nbord, &
                           nbordf,nbm,nbtm,nsm,nsort,nsortf,nstm,poidsa,poidso)

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        Time_loc(2) = real(time,db)
        Time_rout(3) = Time_rout(3) + Time_loc(2) - Time_loc(1)
      endif

      if( .not. ( Second_run .or. Extract ) ) then

! Potential calculation
        if( Flapw ) then
          call potlapw(axyz,Chargat,Coupelapw,deccent,Flapw_new,Full_atom,iapot, &
            iaproto,iaprotoi,icheck(13),igroup,iprabs_nonexc,ipr1,itabs,its_lapw,itypei,itypep,itypepr,Magnetic,mpinodes0, &
            mpirank0,n_atom_0,n_atom_ind,n_atom_proto,natome,natomeq,natomp,ngroup,ngroup_lapw,nklapw,nlm_pot, &
            nlmlapwm,nmatsym,normrmt,npoint,npsom,nrato,nrato_lapw,nrm,nslapwm,nspin,ntype, &
            numat,Orthmat,Overlap,pos,Rato,rho,rhoato,Rlapw,Rmtg,Rmtg0,Rmtimp,rmtsd,Rot_int,rotloc_lapw,rs,rsato,Rsort, &
            Trace_format_wien,Trace_k,Trace_p,V0bdcFimp(1),Vcato,Vh,Vxc,Vxcato,Wien_file,Wien_matsym, &
            Wien_save,Wien_taulap,xyz)
        else

          if( .not. Self_cons_c ) allocate( posi_self(3,natome) )

          call potsup(alfpot,Axe_atom_gr,Bulk_atom_done,Cal_xanes,Chargat,Chargat_init, &
            Chargat_self,Delta_helm,dExc_ex_nex,drho_ex_nex,dv_ex_nex,dvc_ex_nex,Exc_abs_i,Excato,Full_atom,Helm_cos,hybrid, &
            i_self,ia_eq_inv,ia_eq_inv_self,iaabs,iaproto,iaprotoi,iapot,icheck,igreq,igroup,iprabs_nonexc,ipr1,itab, &
            itypei,itypep,itypepr,Lmax_pot,lvval,Magnetic,mpirank0,n_atom_0,n_atom_0_self,n_atom_ind, &
            n_atom_ind_self,n_atom_proto,n_atom_proto_bulk,natome,natome_self,natomeq,natomeq_self,natomp,neqm, &
            ngroup_m,ngroup_nonsph, &
            nhybm,nlat,nlatm,nlm_pot,Nonexc,Nonsph,norbv,normrmt,npoint,npoint_ns,npsom,nrato,nrm,nrm_self,nspin, &
            ntype,numat,Overlap,pop_nonsph,popatm,popatv,pos,posi,posi_self,Proto_calculated,psival,Rato,rho,rho_chg, &
            rho_self,rhoato,rhoato_init,rhoit,rhons,Rmtg,Rmtimp,Rmtg0,Rmtsd,Rot_Atom_gr,Rot_int,rs,rsato,rsort, &
            SCF_c,Self_nonexc,Sym_2D,V_abs_i,V_helm_t,V_intmax,Vc_abs_i,Vcato,Vcato_init,Vh,Vhns,Vsphere,Vxc,Vxcato, &
            V0bdcFimp(1),Width_helm,xyz)

          if( .not. Self_cons_c ) deallocate( posi_self )

        endif

        call Potential_comp(Bulk_atom_done,Cal_xanes,distai(natome),dV0bdcF,Full_atom, &
          Green,i_range,i_self,iaabs,iapot,iaproto,iaprotoi,icheck(13),imoy,imoy_out,iopsymr,ipr1,isrt,itypepr,korigimp,Magnetic, &
          Moy_loc,mpirank0,n_atom_0,n_atom_ind,n_atom_proto,n_atom_proto_bulk,natome,natomp,nim,nlm_pot,npoint,npsom,nptmoy, &
          nptmoy_out,nrato,nrm,nsortf,nspin,nstm,ntype,numat,poidsov,poidsov_out,pos, &
          Rato,Rchimp,rhoato,Rmtg,Rmtg0,Rmtsd,rs,rsbdc,rsbdc_out,rsort,rvol,SCF_c,Sym_2D,V_intmax,V0bdcF,V0bdcF_out,V0bdcFimp, &
          V0muf,Vcato,Vh,Vhbdc,Vhbdc_init,Vhbdc_out,Vr,Vxc,Vxcato,VxcbdcF,VxcbdcF_out,xyz,Workf)

      endif

      call Atom_coop_selec(COOP,Dist_coop,ia_coop,iaprotoi,icheck(28),igr_coop,igreq,igroupi,iprabs_nonexc,itypei,n_atom_coop, &
                         n_atom_proto,nab_coop,natome,neqm,ngreq,nomfich_s,ntype,numat,Posi,Rmtg)

      if( i_range == 1 ) then
        nenerg = nenerg_s
        allocate( Energ(nenerg) )
        allocate( Eimag(nenerg) )
        if( Optic ) then
          if( E_cut_man ) then
            Energ(:) = Energ_s(:) + E_cut_imp
          elseif( Old_zero ) then
            Energ(:) = Energ_s(:) + E_cut
          else
            Energ(:) = Energ_s(:)
          endif
        else
          Energ(:) = Energ_s(:)
        endif
        Eimag(:) = Eimag_s(:)
      end if

      if( Full_atom ) then
        if( nonexc ) then
          iaprabs_nonexc = iaabsi
        else
          iaprabs_nonexc = 0
        endif
      else
        iaprabs_nonexc = iprabs_nonexc
      endif

      if( ipr1 == 1 ) Rmtg(0) = Rmtg(iprabs_nonexc)

      if( .not. Extract ) then

! Complementary work on potential
        if( .not. Second_run ) &
          call Check_Ecmax(Ecineticmax,Ecineticmax_out,Eclie,Eclie_out,Eneg,Energ(nenerg),mpirank,nspin,V0bdcF,V0bdcf_out,Workf)

! Non excited absorbing atom potential (used in Energseuil)
        if( Second_run .or. Second_SCF ) then
          do ispin = 1,nspin
            V_abs_i(1:nrm,ispin) = Vcato(1:nrm,1,iaprabs_nonexc) + Vxcato(1:nrm,1,ispin,iaprabs_nonexc)
          end do
          Vc_abs_i(1:nrm) = Vcato(1:nrm,1,iaprabs_nonexc)
          Exc_abs_i(1:nrm) = Excato(1:nrm,iaprabs_nonexc)
        endif

! Calculation of the core level (Kohn-Sham) energy (Epsii)
        if( i_range == 1 ) then
          if( Optic ) then
            Epsii(:) = 0._db; Epsii_moy = 0._db; Delta_Eseuil = 0._db
          else
            if( Old_zero ) then
              E_zero = E_cut
            else
              E_zero = 0._db
            endif
            call Energseuil(Core_energ_tot,Core_resolved,Delta_Epsii,Delta_Eseuil,E_zero,Epsii,Epsii_moy,Eseuil,Exc_abs_i, &
              icheck(14),is_g,itabs_nonexc,Jseuil,Lseuil,m_g,mpirank0,nbseuil,ninit1,ninit,ninitr,nrato(itabs_nonexc),nrm,nseuil, &
              nspin,ntype,numat_abs,psii,Rato,Rmtg(iprabs_nonexc),Rmtsd(iprabs_nonexc),V_abs_i,V_intmax,Vc_abs_i,V0bdcf,Workf)
          endif
        endif

      endif

! Determination of Lmax
      if( .not. Green ) then
        call clmax(Ecineticmax_out,Rsort,Lmaxso0,Lmaxso_m,-1,.true.)
        Lmaxso_m = min( Lmaxso_m, Lmaxso_max )
      else
        Lmaxso_m = 0
      endif
      Lmax_max = 0
      if( Full_atom ) then
        iaprabs = iaabsi
      else
        iaprabs = iprabs
      endif

      nlm_p_fp = 1
      if( Full_Potential ) then
        nlm_p_fp = nlm_probe
      elseif( Hubbard ) then
        if( Hubb(itabs) .and. .not. Hubb_diag_abs ) nlm_p_fp = nlm_probe
      endif

      Lmax_abs = 0
      do ipr = 0,n_atom_proto
        Z = numat( itypepr(ipr) )
        call clmax(Ecineticmax,Rmtg(ipr),Lmaxat0,Lmaxat(ipr),Z,lmaxfree)
        if( Z == numat_abs ) then
          Lmaxat(ipr) = max( Lmax_probe, Lmaxat(ipr) )
          Lmax_abs = max( Lmaxat(ipr), Lmax_abs )
        endif
        Lmax_max = max(Lmax_max,Lmaxat(ipr))
      end do

      do ipr = 0,n_atom_proto
        Z = numat( itypepr(ipr) )
        if( Z == numat_abs ) Lmaxat(ipr) = Lmax_abs
      end do

      nlmmax = (Lmax_max + 1 )**2
      if( Tddft .or. Optic ) then
! nlmamax is used only for Tddft, it is limited to save memory space
        Lmax_abs_t = min( Lmax_abs, Lmax_tddft)
        nlmamax = (Lmax_abs_t + 1 )**2
      else
        Lmax_abs_t = 0
        nlmamax = 0
      endif

      nlmsam = nspinp * nlmmax
      nlmomax = ( Lmaxso_m + 1 )**2
      nso1 = nspino * nlmomax

      if( .not. Green .and. FDM_comp ) then
        Lmax_comp = Lmax_probe
        nlm_comp = ( Lmax_probe + 1 )**2
      else
        nlm_comp = 0
        Lmax_comp = 0
      endif
      allocate( Tau_comp(nlm_comp,nspinp,nlm_comp,nspinp,nenerg) )
      Tau_comp(:,:,:,:,:) = (0._db,0._db)
      allocate( RadialIntegral(nlm_comp,1,nspinp,nlm_comp,1,nspinp,nenerg) )
      RadialIntegral(:,:,:,:,:,:,:) = (0._db,0._db)

! Absorbing atom potential and density storage
      if( i_range == 1 .and. .not. Extract ) then
        if( Full_atom ) then
          iaprabs = iaabsi
        else
          iaprabs = iprabs
        endif
        Rmtg_abs = Rmtg(iprabs)
        Rmtsd_abs = Rmtsd(iprabs)
        allocate( rato_abs(nrato_abs) )
        allocate( rhoato_abs(nrato_abs,nspin) )
        allocate( rsato_abs(nrato_abs) )
        allocate( Vcato_abs(nrato_abs,nlm_pot) )
        allocate( Vxcato_abs(nrato_abs,nlm_pot,nspin) )
        rato_abs(1:nrato_abs) = Rato(1:nrato_abs,itabs)
        rhoato_abs(1:nrato_abs,:) = rhoato(1:nrato_abs,:,iaprabs)
        rsato_abs(1:nrato_abs) = rsato(1:nrato_abs,iaprabs)
        Vcato_abs(1:nrato_abs,:) = Vcato(1:nrato_abs,:,iaprabs)
        Vxcato_abs(1:nrato_abs,:,:) = Vxcato(1:nrato_abs,:,:,iaprabs)
        if( icheck(13) > 0 .and. mpirank0 == 0 ) call Potential_writing(Delta_Eseuil,dV0bdcf,E_cut,E_Fermi,Ecineticmax,Epsii, &
                           Epsii_moy,Hubb_abs,Hubb_diag_abs,Lmax_abs,m_hubb,ninitr,nlm_pot,nrato_abs,nspin,nspinp, &
                           rato_abs,rhoato_abs,Rmtg_abs,Rmtsd_abs,rsato_abs,rsbdc,V_Hubb_abs,V0muf,Vcato_abs,Vhbdc,Vxcato_abs, &
                           VxcbdcF)
      endif

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        Time_loc(3) = real(time,db)
        Time_rout(4) = Time_rout(4) + Time_loc(3) - Time_loc(2)
      endif

      allocate( iato(nlmsam,natome,ngrph) )
      allocate( lato(nlmsam,natome,ngrph) )
      allocate( iso(nso1,ngrph) )
      allocate( lso(nso1,ngrph) )
      allocate( mato(nlmsam,natome,ngrph) )
      allocate( mso(nso1,ngrph) )
      allocate( nlmso0(ngrph) )
      allocate( nlmsa0(natome,ngrph) )

! Determination of the (l,m) of each representation
      if( .not. ( Second_run .or. Extract ) ) &
        call lmrep(Green,iaprotoi,iato,icheck(15),iopsym_atom,iopsymr,irep_util,iso,itypei,karact,lato, &
               Lmaxat,Lmaxso_m,lso,mato,mpirank0,mso,n_atom_proto,natome,nbordf,ngrph,nlmsa0,nlmsam,nlmso0, &
               nso1,nsortf,nspino,ntype,numat,Orthmati,posi,rot_atom)

! Ylm calculation on the FDM grid points
      if( Green ) then
        nbtm_fdm = 0
      else
        nbtm_fdm = nbtm
      endif
      allocate( Ylmato(nbtm_fdm,nlmmax,natome) )
      if( .not. ( Green .or. Second_run .or. Extract ) )  then
        allocate( Ylmso(nsort,nlmomax) )
        call ylmpt(iaprotoi,ibord,icheck(15),iopsymr,isrt,Lmaxat,Lmaxso_m,n_atom_proto,natome, &
           nbord,nbtm,nlmmax,nlmomax,nsort,nstm,npsom,posi,rot_atom,xyz,Ylmato,Ylmso)
      endif

      if( .not. ( Green .or. Second_run .or. Extract ) ) then
        if( Spinorbite .or. Relativiste ) then
          nicm = nim
        else
          nicm = 1
        endif
        allocate( gradvr(nicm,3,nspin))
      endif
      allocate( Lmaxa(natome) )
      allocate( nlmsa(natome) )
      allocate( drho_self(nrm_self,nlm_pot,nspinp,n_atom_0_self:n_atom_ind_self,0:mpinodee-1) )
      allocate( Statedens(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodee-1) )
      allocate( Statedens_i(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodee-1) )
      allocate( Coverlap(16,16,nspinp,nab_coop,0:mpinodee-1) )
      allocate( index_coop(2,nab_coop) )

      ninitlv = nbseuil
      n_vr_0 = n_atom_0
      n_vr_ind = n_atom_ind

      if( i_range == 1 ) then
        if( Tddft_xanes ) then
          nenerg_tddft = nenerg
        else
          nenerg_tddft = 0
        endif
!  nlmmax = (Lmax_max + 1 )**2   max dimension
        allocate( rof0(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil) )
        allocate( Taull_dft(nlmmax,nspinp,nlmmax,nspinp,0:mpinodee-1) )
        Taull_dft(:,:,:,:,:) = ( 0._db, 0._db )

        if( Tddft_xanes .or. Optic_xanes ) then
          rof0(:,:,:,:,:) = (0._db,0._db)
          allocate( Taull_tdd(nenerg,nlmamax,nspinp,nlmamax,nspinp) )
          Taull_tdd(:,:,:,:,:) = (0._db, 0._db)
        endif

        if( RIXS ) then
          nenerg_RIXS = nenerg
          nlm_rixs = min( nlmmax, 16 )
          if( Extract ) then
            nlmomax = extract_nlmomax(multi_run,nom_fich_extract)
            Lmaxso_m = nint( sqrt( real( nlmomax, db ) ) ) - 1
          endif
        else
          nenerg_RIXS = 0
          nlm_rixs = 0
        endif

        allocate( Ampl_dft(nlm_rixs,nspinp,nlmomax,nspinp,0:mpinodee-1) )
        Ampl_dft(:,:,:,:,:) = ( 0._db, 0._db )

        allocate( Ampl_abs(nenerg_RIXS,nlm_probe,nspinp,nlmomax,nspinp) )
        Ampl_abs(:,:,:,:,:) = (0._db, 0._db)

        allocate( sec_atom(ninitr) )
        allocate( secdd(3,3,n_rel,ninitr,0:mpinodee-1) )
        allocate( secmd(3,3,ninitr,0:mpinodee-1) )
        allocate( secmm(3,3,ninitr,0:mpinodee-1) )
        allocate( secdq(3,3,3,ninitr,0:mpinodee-1) )
        allocate( secdo(3,3,3,3,ninitr,0:mpinodee-1) )
        allocate( secoo(3,n_oo,3,n_oo,ninitr,0:mpinodee-1) )
        allocate( secqq(3,3,3,3,ninitr,0:mpinodee-1) )
        allocate( secdd_m(3,3,n_rel,ninitr,0:mpinodee-1) )
        allocate( secmd_m(3,3,ninitr,0:mpinodee-1) )
        allocate( secmm_m(3,3,ninitr,0:mpinodee-1) )
        allocate( secdq_m(3,3,3,ninitr,0:mpinodee-1) )
        allocate( secdo_m(3,3,3,3,ninitr,0:mpinodee-1) )
        allocate( secoo_m(3,n_oo,3,n_oo,ninitr,0:mpinodee-1) )
        allocate( secqq_m(3,3,3,3,ninitr,0:mpinodee-1) )
        allocate( S_nrixs(nq_nrixs,(Lmax_nrixs+1)**2,(Lmax_nrixs+1)**2,ninitr,0:mpinodee-1) )
        allocate( S_nrixs_m(nq_nrixs,(Lmax_nrixs+1)**2,(Lmax_nrixs+1)**2,ninitr,0:mpinodee-1) )

      endif

      Green_i = Green_int

      Tau_nondiag = .false.
      if( Spinorbite .or. Full_potential ) then
        Tau_nondiag = .true.
      elseif( Hubbard ) then
        do iapr = n_atom_0_self,n_atom_ind_self
          if( Hubb_diag(iapr) ) cycle
          Tau_nondiag = .true.
          exit
        end do
      endif

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        Time_loc(4) = real(time,db)
        Time_rout(5) = Time_rout(5) + Time_loc(4) - Time_loc(3)
      endif

      if( i_range == 1 ) then
        ie0 = 1
      else
        ie0 = nenerg0 + 1
      endif
      if( i_range == n_range ) then
        nenerg0 = nenerg
      else
        do ie = ie0,nenerg
          if( Energ(ie) > E_max_range(i_range) + eps10 ) exit
        end do
        nenerg0 = ie - 1
      endif

      if( n_range /= 1 ) then
        if( abs( E_radius(i_radius) - E_max_range(i_range) ) < eps10 ) i_radius = i_radius + 1
        if( abs( E_adimp(i_adimp) - E_max_range(i_range) ) < eps10 ) i_adimp = i_adimp + 1
      endif

! Loop over energy
      nge = ( nenerg0 - ie0 ) / mpinodee + 1

! One makes the storage only when mpinodee > 1, otherwise the writing is a temporary file
      if( One_run .and. multi_run_e == 1 .and. i_range == 1 ) &
        allocate( Taull_stk(nlm_probe,nspinp,nlm_probe,nspinp,2:n_multi_run,nge) )
      index_e = 0

      First_E = .true.
      boucle_energ_xanes: do je = 1,nge

        if( mpinodee /= mpinodes .and. mpirank0 /= 0 ) cycle

        ie = ie0 + ( je - 1 ) * mpinodee + mpirank

        if( ie <= nenerg0 ) then

          if( mpirank0 == 0 ) then
            call CPU_TIME(time)
            Time_loc(1) = real(time,db)
          endif

          index_e = index_e + 1
          if( Old_zero ) then
            Enervide = Energ(ie) - Workf
          else
            Enervide = Energ(ie) + E_Fermi
          endif

          if( icheck(16) < 2 .and. ie > 1 ) then
            Trace_l = 0
          else
            Trace_l = Trace_k
          endif

! Calculation of potential in the excited state
          if( Extract ) then
            call E_reading(icheck(19),Eimag(ie),Energ(ie),ie,multi_run_e,nom_fich_extract)
            allocate( Vra(nrato_abs,nlm_pot,nspin) )
            call Potex_abs(alfpot,dV0bdcF,Enervide,icheck(16),nlm_pot,nrato_abs,nspin,Optic,rato_abs,rsato_abs,rsbdc,V_intmax, &
                           V0bdc,Vcato_abs,Vhbdc,Vra,Vxcato_abs,VxcbdcF)
          else
            call potex(Nonsph,axyz,alfpot,dV0bdcF,Energ(ie),Enervide,Full_atom, &
              iaabsi,iapot,iaprotoi,icheck(16),iprabs_nonexc,iprabs,itab,itypepr,Magnetic, &
              n_atom_0,n_atom_ind,n_atom_proto,n_vr_0,n_vr_ind,natome,nlm_pot,Nonexc_g,npoint,npoint_ns,npsom,nptmoy, &
              nptmoy_out,nrato,nrm,nspin,ntype,Optic,Rato,rho,rhons,Rmtg,rs,rsato,rsbdc,rsbdc_out,Trace_l,Trace_p, &
              V_intmax,Vcato,Vh,Vhbdc,Vhbdc_out,Vhns,Vr,Vxc,VxcbdcF,VxcbdcF_out,Vrato,Vxcato,V0bdc,V0bdc_out,xyz)
          endif

          do ispin = 1,nspin

            V0bdc_out(ispin) = V0bdc(ispin)

            if( .not. Extract ) then

              if( Muffintin ) call mdfmuf(Nonsph,Axe_Atom_gr,axyz, &
                Full_atom,iaabs,iaproto,iaprotoi,icheck(16),igreq,igroup,ispin,itypep,n_atom_0,n_atom_ind, &
                n_atom_proto,natome,natomp,neqm,ngroup_m,nlm_pot,npoint,npoint_ns,npsom,nrato,nrm,nspin,ntype,pos,Rato,rho, &
                rhons,Rmtg,Trace_l,Trace_p,Vhns,V0bdc(ispin),Vr,Vrato,xyz)

              if( Supermuf ) call modmuf(Full_atom,iaprotoi,itypepr,icheck(16), &
                ispin,n_atom_0,n_atom_ind,n_atom_proto,natome,nlm_pot,nrato,nrm,nspin,ntype,Rato,Rmtg,V0bdc(ispin),Vrato)

              if( .not. ( Green .or. Second_run ) .and. ( Spinorbite .or. Relativiste ) ) &
                call gradpot(Base_hexa,cgrad,gradvr,icheck(16),iord,ispin,ivois,nicm,npoint,npsom,nspin,nvois,Vr)

            endif

            Ecinetic_out(ispin) = Enervide - V0bdc_out(ispin)
            Ecinetic(ispin) = Enervide - V0bdc(ispin)

            if( .not. Eneg ) then
              Ecinetic_out(ispin) = max( Ecinetic_out(ispin), Eclie_out )
              Ecinetic(ispin) = max( Ecinetic(ispin), Eclie )

              if( Ecinetic(ispin) < eps10 .or. Ecinetic_out(ispin) < eps10 ) then
                if( mpirank0 == 0 ) then
                  call write_error
                  do ipr = 6,9,3
                    write(ipr,170) Ecinetic(ispin) * rydb
                    if( nptmoy_out > 0 ) write(ipr,180) Ecinetic_out(ispin) * rydb
                  end do
                endif
                stop
              endif
            endif

          end do  ! end of loop over spin

          if( Eneg ) then
! Numerical problem around zero
            em = 0.01_db / rydb
            ei0 = 0.01_db / rydb
            E =  min( Ecinetic_out(1), Ecinetic_out(nspin) )
            if( E < 0._db ) then
              Eimag(ie) = max( Eimag(ie), ei0 )
            elseif( E < em ) then
              eii = ei0 * ( 1 - E / em )
              Eimag(ie) = max( Eimag(ie), eii )
            endif
          endif

          if( mpirank0 == 0 ) then
            call CPU_TIME(time)
            Time_loc(2) = real(time,db)
            Time_rout(6) = Time_rout(6) + Time_loc(2) - Time_loc(1)
          endif

          if( abs( Eimag(ie) ) > eps10 .or. Ecinetic_out(1) < eps10 .or. Ecinetic_out(nspin) < eps10 ) then
            E_comp = .true.
          else
            E_comp = .false.
          endif

          if( No_solsing ) then
            Solsing = .false.
          elseif( E_comp .and. ( Green .or. FDM_comp ) ) then
!          elseif( E_comp .and. ( Green .or. FDM_comp_m ) ) then
            Solsing = .true.
          else
            Solsing = Solsing_s
          endif

          Ecmax = 0._db
          Ecmax_out = 0._db
          do ispin = 1,nspin
            Ecmax_out = max( Ecmax_out, Ecinetic_out(ispin) )
            Ecmax = max( Ecmax, Ecinetic(ispin) )
          end do

          Lmax_g = 0
          Lmax_abs_a = 0
          do ipr = 0,n_atom_proto
            Z = numat( itypepr(ipr) )
            call clmax(Ecmax,Rmtg(ipr),Lmaxat0,Lmaxat(ipr),Z,lmaxfree)
            Lmaxat(ipr) = min(Lmaxat(ipr), Lmax_max)
            if( Z == numat_abs ) then
               Lmaxat(ipr) = max( Lmax_probe, Lmaxat(ipr) )
               Lmax_abs_a = max( Lmax_abs_a, Lmaxat(ipr) )
            end if
            Lmax_g = max( Lmaxat(ipr), Lmax_g)
          end do
          do ipr = 0,n_atom_proto
            Z = numat( itypepr(ipr) )
            if( Z == numat_abs ) Lmaxat(ipr) = Lmax_abs_a
          end do
          nlmagm = ( Lmax_g + 1 )**2

          allocate( Tau_ato(nlmagm,nspinp,nlmagm,nspinp,n_atom_0:n_atom_ind) )
          if( Green ) then
            nphiato1 = 0; nphiato7 = 0
          else
            nphiato1 = nbtm
            nphiato7 = 1
            if( abs(Eimag(ie)) > eps10 .or. Spinorbite .or. Ylm_comp ) then
              nphiato7 = 2
            else
              do igrph = 1,ngrph
                if( .not. Repres_comp(igrph) ) cycle
                nphiato7 = 2
                exit
              end do
            endif
          endif

          allocate( phiato(nphiato1,nlmagm,nspinp,nspino,natome,nphiato7) )
          allocate( Taull(nlmagm,nspinp,nlmagm,nspinp,natome) )
          allocate( Tau_coop(nlmagm,nspinp,nlmagm,nspinp,nab_coop,nspino) )
          phiato(:,:,:,:,:,:) = 0._db
          Tau_ato(:,:,:,:,:) = (0._db,0._db)
          Taull(:,:,:,:,:) = (0._db,0._db)
          Tau_coop(:,:,:,:,:,:) = (0._db,0._db)

          if( .not. Green ) then
            call clmax(Ecmax_out,Rsort,Lmaxso0,Lmaxso,-1,.true.)
            Lmaxso = min( Lmaxso, Lmaxso_m )
            if( icheck(18) > 0 ) write(3,'(/a9,i4)') ' Lmaxso =', Lmaxso
          endif

          if( icheck(18) > 1 ) write(3,190)

! Loop over non-equivalent atoms
          n1 = n_atom_0

          do iapr = n1,n_atom_ind

            if( Second_run .or. Extract ) exit

            if( Full_atom ) then
              ipr = iaprotoi( iapr )
            else
              ipr = iapr
            endif
            Z = numat( itypepr(ipr) )
            if( Z == 0 .and. .not. Green ) cycle
            if( Z == Z_nospinorbite ) then
              Spino = .false.
            else
              Spino = Spinorbite
            endif
            it = itypepr(ipr)
            if( Full_atom ) then
              Absorbeur = iapr == iaabsi
            else
              Absorbeur = iapr == iprabs
            endif

            allocate( V_hubb_t(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp))
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

            nr = nrato(it)
            allocate( r(nr) )
            allocate( Vrato_e(nr,nlm_pot,nspin) )
            r(1:nr) = Rato(1:nr,it)

            Lmax = Lmaxat(ipr)

            Vrato_e(1:nr,:,:) = Vrato(1:nr,:,:,iapr)

            call Sphere(Axe_atom_grn,Ecinetic,Eimag(ie),Energ(ie),Enervide,Full_atom, &
              Full_potential,Green,Hubb_a,Hubb_d,iaabsi,iapr,iaprotoi,ibord,icheck(18),igreq,igroupi,iopsymr,Lmax,Lmax_pot, &
              m_hubb,n_atom_0,n_atom_ind,n_atom_proto,natome,nbord,nbtm,nbtm_fdm,neqm,ngroup_m, &
              nlm_pot,nlmagm,nlmmax,nphiato1,nphiato7,npsom,nr,nspin,nspino,nspinp,Z,phiato,posi,r,Relativiste,Renorm, &
              Rmtg(ipr),Rmtsd(ipr),Spino,Tau_ato,V_hubb_t,V_intmax,V0bdc,Vrato_e,xyz,Ylm_comp,Ylmato)

            deallocate( r )
            deallocate( Vrato_e, V_hubb_t )

          end do

          if( mpirank0 == 0 ) then
            call CPU_TIME(time)
            Time_loc(3) = real(time,db)
            Time_rout(7) = Time_rout(7) + Time_loc(3) - Time_loc(2)
          endif

! nspinorb = 2 if magnetic without Spin-orbit
! nspinorb = 1 if non-magnetic or Spin-orbit
! nspino = 2 if Spin-orbit

! igrph: loop over representations

          Lmax_g = 0
          do ia = 1,natome
            Lmaxa(ia) = Lmaxat( iaprotoi(ia) )
            Lmax_g = max(Lmax_g,Lmaxa(ia))
          end do
          Time_fill = 0._db; Time_tria = 0._db

          if( .not. ( Second_run .or. Extract ) ) then

            do ispin = 1,nspinorb
              do igrph = 1,ngrph

                nlmsamax = 0
                nlmsa(:) = 0
                do ia = 1,natome
                  Z = numat(itypei(ia))
                  if( .not. Green .and. Z == 0 ) cycle
                  n = nlmsa0(ia,igrph)
                  allocate( Lval(n) )
                  Lval(1:n) = lato(1:n,ia,igrph)
                  call cnlmmax(Lmaxa(ia),Lval,n,nlmsa(ia))
                  nlmsamax = max(nlmsamax,nlmsa(ia))
                  deallocate( Lval )
                end do
                if( nlmsamax == 0 .and. igrph /= ngrph ) cycle
                ich = icheck(19)

                if( Green ) then
                  call msm(Axe_Atom_grn,Cal_xanes,Classic_irreg,Dist_coop,Ecinetic,Eimag(ie),Full_atom, &
                    ia_coop,ia_eq,ia_rep, &
                    iaabsi,iaprotoi,iato,ich,igroupi,igrph,iopsymr,irep_util,is_eq,ispin,karact,lato,Lmaxa,Lmax_g,mato, &
                    n_atom_0,n_atom_coop,n_atom_ind,n_atom_proto,nab_coop,natome,natomp,nb_eq,nb_rpr,nb_rep_t,nb_sym_op,nchemin, &
                    ngroup_m,ngrph,nlmagm,nlmsa,nlmsam,nlmsamax,Normaltau,nspin,nspino,nspinp,pos,posi,Recop,Repres_comp(igrph), &
                    Rmtg,rot_atom,Solsing,Spinorbite,State_all,Sym_cubic,Tau_ato,Tau_nondiag,Tau_coop,Taull,Time_fill, &
                    Time_tria,Ylm_comp)
                else
                  n = nlmso0(igrph)
                  allocate( Lval(n) )
                  Lval(1:n) = lso(1:n,igrph)
                  call cnlmmax(Lmaxso,Lval,n,nlmso)
                  deallocate( Lval )
                  call mat(Adimp,Ampl_dft,Atom_axe,Axe_Atom_grn,Base_hexa,Basereel,Cal_xanes,cgrad, &
                      clapl,Classic_irreg,Dist_coop,distai,E_comp,Ecinetic_out,Eclie_out,Eimag(ie),Eneg,Enervide, &
                      FDM_comp_m,Full_atom,gradvr,ia_coop,iaabsi,iaprotoi,iato,ibord,ich,ie,igroupi,igrph,irep_util, &
                      isbord,iso,ispin,isrt,ivois,isvois,karact,lato,Lmaxa,Lmaxso,lso,mato,mpinodes,mpirank0,mso,nab_coop, &
                      natome,n_atom_0,n_atom_coop,n_atom_ind,n_atom_proto,nbm,nbord,nbordf,nbtm,ngroup_m,ngrph,nim, &
                      nicm,nlm_rixs,nlmagm,nlmmax,nlmomax,nlmsa,nlmsam,nlmso,nphiato1,nphiato7,npoint,npr,npsom,nsm, &
                      nsort,nsortf,nso1,nspin,nspino,nspinp,nstm,numia,nvois,phiato,poidsa,poidso,Posi,R_rydb,Recop,Relativiste, &
                      Repres_comp(igrph),RIXS,Rmtg,Rsort,rvol,Rydberg,Solsing,Spinorbite,State_all,Sym_cubic,Tau_ato,Tau_coop, &
                      Taull,Time_fill,Time_tria,V0bdc_out,Vr,xyz,Ylm_comp,Ylmato,Ylmso)
                endif

              end do  ! end of loop over representations
            end do  ! end of loop over spin

          endif

          if( One_run .and. .not. Extract ) &
            call Data_one_run(iabsm,icheck(19),igreq,index_e,igroupi,ipr1,Lmax_probe,Lmaxa, &
              multi_run_e,n_atom_proto,n_multi_run,natome,neqm,nge,ngreq,nlm_probe,nlmagm,nspinp, &
              Rot_atom,Rot_int,Spinorbite,Taull,Taull_stk,Ylm_comp)

          if( Extract ) call Tau_reading(icheck(19),ie,iaabsi,Lmaxa(iaabsi),natome,nenerg,nlmagm,nspinp,RIXS,Taull)

          if( Extract .and. RIXS ) call Ampl_reading(Ampl_dft,icheck(19),ie,Lmaxso,mpinodee,mpirank0,nenerg,nlm_probe,nlm_rixs, &
                                                     nlmomax,nspinp)

          nlm = ( Lmaxa(iaabsi) + 1 )**2
          Taull_dft(1:nlm,:,1:nlm,:,mpirank) = Taull(1:nlm,:,1:nlm,:,iaabsi)
          lmax_dft(:,mpirank) = Lmaxat(:)

          if( RIXS ) Ampl_abs(ie,1:nlm_probe,:,1:(Lmaxso+1)**2,:) = Ampl_dft(1:nlm_probe,:,1:(Lmaxso+1)**2,:,mpirank0)

          if( Optic_xanes .or. Tddft_xanes ) then
! Can be done like this because Eimag = 0 in tddft and optic
            Lm1 = 0
            do l1 = 0,min(Lmaxa(iaabsi),Lmax_tddft)
              do m1 = -l1,l1
                Lm1 = Lm1 + 1
                Lm2 = 0
                do l2 = 0,min(Lmaxa(iaabsi),Lmax_tddft)
                  do m2 = -l2,l2
                    Lm2 = Lm2 + 1
                    do isp1 = 1,nspinp
                      do isp2 = 1,nspinp
! In FDM, next line does not change nothing
                        Taull_tdd(ie,Lm1,isp1,Lm2,isp2) = 0.5_db * ( Taull(Lm1,isp1,Lm2,isp2,iaabsi) &
                                                            - conjg( Taull(Lm2,isp2,Lm1,isp1,iaabsi) ) )
                      end do
                    end do
                  end do
                end do
              end do
            end do
            if( Optic_Xanes ) then
              if( ie == 1 ) write(6,'(/A)') '   Energy   -Tau_i(l,m=0,l,m=0) l = 0,Lmax_probe'
              write(6,'(f10.3,1p,6e13.5)') Energ(ie)*rydb, ( - aimag( Taull_tdd(ie,l**2+l+1,1,l**2+l+1,1) ), &
                   l = 0,Lmax_probe )
            endif
          endif

          if( mpirank0 == 0 ) then
            call CPU_TIME(time)
            Time_loc(4) = real(time,db)
            Time_rout(8) = Time_rout(8) + Time_loc(4) - Time_loc(3)
            Time_rout(9) = Time_rout(9) + Time_fill ! Filling
            Time_rout(10) = Time_rout(10) + Time_tria ! Triang
          endif

! Calculation of cartesian tensors
          if( .not. Optic ) then

            ns_dipmag = 1
            ndim2 = 1
            allocate( taull_abs(nlms_pr,nlms_pr,ndim2,ndim2,2,2,ns_dipmag) )

            do i = 1,2

              if( i == 1 .and. .not. Xan_atom ) cycle

              taull_abs(:,:,:,:,:,:,:) = (0._db,0._db)
              if( i == 1 ) then
                if( Full_atom ) then
                  iapr = iaabsi
                elseif( Nonexc ) then
                  iapr = iprabs_nonexc
                else
                  iapr = 0
                endif
              endif

              Lms1 = 0
              do Lm1 = 1,nlm_probe
! For the atomic part, all signal is in the singular solution
                if( i == 1 .and. Solsing ) exit
                do iso1 = 1,nspino
                  Lms1 = Lms1 + 1
                  Lms2 = 0
                  do Lm2 = 1,nlm_probe
                    do iso2 = 1,nspino
                      Lms2 = Lms2 + 1
                      do isp1 = 1,2
                        do isp2 = 1,2
                          if( Spinorbite ) then
                            iss1 = iso1; iss2 = iso2
                          else
                            if( isp1 /= isp2 ) cycle
                            iss1 = min(isp1,nspin)
                            iss2 = min(isp2,nspin)
                          endif

                          if( i == 1 ) then
                            taull_abs(Lms1,Lms2,1,1,isp1,isp2,1) = 0.5_db * &
                                                ( Tau_ato(Lm1,iss1,Lm2,iss2,iapr) -  conjg( Tau_ato(Lm2,iss2,Lm1,iss1,iapr) ) )
                          else
                            taull_abs(Lms1,Lms2,1,1,isp1,isp2,1) = Taull(Lm1,iss1,Lm2,iss2,iaabsi)
                          endif
                        end do
                      end do
                    end do
                  end do
                end do
              end do


              if( .not. Extract .and. ( i == 1 .or. .not. Xan_atom ) ) then
                allocate( Vra(nrato_abs,nlm_pot,nspin) )
                if( Full_atom ) then
                  Vra(1:nrato_abs,:,:) = Vrato(1:nrato_abs,:,:,iaabsi)
                else
                  Vra(1:nrato_abs,:,:) = Vrato(1:nrato_abs,:,:,iprabs)
                endif
              endif
              if( Tddft ) then
                ip00 = 0
              else
                ip00 = ip0
              endif
              n_V = 1
              n_Ec = 1
              allocate( Enervide_t(n_Ec) )
              Enervide_t(1) = Enervide
              if( mpirank_in_mumps_group == 0 ) then
                if( NRIXS ) &
                  call S_nrixs_cal(Classic_irreg,coef_g,Core_resolved,Ecinetic, &
                    Eimag(ie),Energ(ie),Enervide_t,Eseuil,FDM_comp_m,Final_tddft,Full_potential,Green,Green_i,Hubb_abs, &
                    Hubb_diag_abs,icheck(20),l0_nrixs,Lmax_nrixs,is_g,Lmax_probe,Lmax_pot,lmoins1,lplus1,Lseuil,m_g,m_hubb, &
                    mpinodee,mpirank,n_Ec,n_V,nbseuil,ns_dipmag,ndim2, &
                    ninit1,ninit,ninitr,ninitlv,nlm_pot,nlm_probe,nlm_p_fp,nq_nrixs,nrato_abs,nrm,nspin,nspino,nspinp, &
                    numat(itabs),psii,q_nrixs,rato_abs,Relativiste,Renorm,Rmtg_abs,Rmtsd_abs, &
                    S_nrixs,S_nrixs_m,Solsing,Solsing_only,Spinorbite,Taull_abs, &
                    V_hubb_abs,V_intmax,V0bdc,Vra,Ylm_comp)

                call tenseur_car(Classic_irreg,coef_g,Core_resolved,Ecinetic,Eimag(ie), &
                    Energ(ie),Enervide_t,Eseuil,FDM_comp_m,Final_optic,Final_tddft,Full_potential,Green,Green_i,Hubb_abs, &
                    Hubb_diag_abs,icheck(20),ie,ip_max,ip00,is_g,Lmax_probe,Lmax_pot,ldip,lmoins1,loct,lplus1,lqua,Lseuil,m_g, &
                    m_hubb,mpinodee,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi, &
                    Multipole,n_Ec,n_oo,n_rel,n_V,nbseuil,ns_dipmag,ndim2,nenerg_tddft,ninit1,ninit,ninitr,ninitlv,nlm_pot, &
                    nlm_probe,nlm_p_fp,nlmamax,nrato_abs,nrm,nspin,nspino,nspinp,numat(itabs),psii,rato_abs,Relativiste,Renorm, &
                    .false., &
                    Rmtg_abs,Rmtsd_abs,rof0,rot_atom_abs,Rot_int,secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm, &
                    secmm_m,secoo,secoo_m,secqq,secqq_m,Solsing,Solsing_only,Spinorbite,Taull_abs,Tddft,V_hubb_abs,V_intmax, &
                    V0bdc,Vra,Ylm_comp)
              endif
              deallocate( Enervide_t )

              if( i == 1 ) then

                sec_atom(:) = 0._db
                do initl = 1,ninitr
                  do j = 1,3
                    sec_atom(initl) = sec_atom(initl) + Real(secdd(j,j,1,initl,mpirank), db)
                  end do
                end do
                sec_atom(:) = sec_atom(:) / 3

              endif

            end do

            deallocate( taull_abs )

          endif

          deallocate( Tau_ato )
          deallocate( phiato )

          if( mpirank0 == 0 ) then
            call CPU_TIME(time)
            Time_loc(5) = real(time,db)
            Time_rout(11) = Time_rout(11) + Time_loc(5) - Time_loc(4)
          endif

          ich = icheck(28)
          Solsing_o = Solsing_only

          if( Density .or. FDM_comp ) then
            if( Extract ) then
              if ( Full_atom ) then
                Vrato(1:nrato_abs,:,:,iaabsi) = Vra(1:nrato_abs,:,:)
              else
                Vrato(1:nrato_abs,:,:,iprabs) = Vra(1:nrato_abs,:,:)
              endif
            endif
            if( Full_atom ) then
              iaprabs = iaabsi
            else
              iaprabs = iprabs
            endif
            if( Density ) &
              call Cal_dens(Cal_xanes,Classic_irreg,Density_comp,drho_self,Ecinetic,Eimag(ie),Energ(ie),Enervide,Full_atom, &
                Full_potential,Harm_cubic,Hubb,Hubb_diag,Hubb_diag_abs,iaabsi,iaprotoi,ich,itypei,itypepr,lla2_state, &
                Lmax_pot,Lmaxat, &
                m_hubb,mpinodee,mpirank,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nlmagm,nlm_pot, &
                nrato,nrm,nrm_self,nspin,nspino,nspinp,ntype,numat,Rato,Relativiste,Renorm,Rmtg,Rmtsd,Solsing,Solsing_o, &
                Spinorbite,State_all_out,Statedens,Statedens_i,Taull,V_hubb,V_hubb_abs,V_intmax,V0bdc,Vrato,Ylm_comp)
            if( .not. Green .and. FDM_comp ) &    ! under construction
              call Cal_Tau_comp(Cal_xanes,Ecinetic,Eimag(ie),Energ,Enervide,Full_atom, &
                Full_potential,Hubb,Hubb_diag,Hubb_diag_abs,iaabsi,iaprotoi,ich,ie,itypei,itypepr,Lmax_comp,Lmax_pot, &
                Lmaxat,m_hubb,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nenerg,nlm_comp,nlmagm, &
                nlm_pot,nrato,nrm,nspin,nspino,nspinp,ntype,numat,RadialIntegral,Rato,Relativiste,Renorm,Rmtg,Rmtsd, &
                Spinorbite,State_all_out,Taull,Tau_comp,V_hubb,V_hubb_abs,V_intmax,V0bdc,Vrato,Ylm_comp)
          endif

          if( COOP ) call Cal_COOP(Coop_z_along_bond,Coverlap,Density_comp,Dist_coop,Ecinetic,Eimag(ie),Energ(ie),Enervide, &
            Full_atom,Full_potential,Harm_cubic,Hubb,Hubb_diag,Hubb_diag_abs,ia_coop,iaabsi,iaprotoi,icheck(28),ie,index_coop, &
            iprabs,itypepr,Lmax_pot,Lmaxat,m_hubb,mpinodee,mpirank,n_atom_0,n_atom_0_self, &
            n_atom_coop,n_atom_ind,n_atom_ind_self,n_atom_proto,nab_coop,natome,nlm_pot,nlmagm,nrato,nrm, &
            nspin,nspino,nspinp,ntype,numat,Posi,Rato,Relativiste,Renorm,Rmtg,Rmtg0,Rmtsd,Rot_atom,Spinorbite,Tau_coop,V_hubb, &
            V_hubb_abs,V_intmax,V0bdc,Vrato,Ylm_comp,Z_nospinorbite)

          deallocate( Tau_coop, Taull )
          if( Extract .or. ( .not. Optic .and. .not. Extract ) ) deallocate( Vra )

          if( mpirank0 == 0 ) then
            call CPU_TIME(time)
            Time_loc(6) = real(time,db)
            Time_rout(12) = Time_rout(12) + Time_loc(6) - Time_loc(5)
          endif

        endif ! end of condition on ie <= nenerg

        if( mpinodee > 1 ) then

          call MPI_RECV_Tau(lmax_dft,mpinodee,mpirank,mpirank0,n_atom_proto,nlmmax,nspinp,Taull_dft)

          if( Density ) &
            call MPI_RECV_statedens(lla2_state,Lmaxat,mpinodes,mpirank,mpirank0,n_atom_0, &
                                      n_atom_ind,n_atom_proto,nspinp,Statedens,Statedens_i)

          if( COOP ) &
            call MPI_RECV_COOP(mpinodes,mpirank,mpirank0,nab_coop,nspinp,Coverlap)

          if( .not. Optic ) &
            call MPI_RECV_all(l0_nrixs,Lmax_nrixs,mpinodee,mpirank,mpirank0,Multipole,n_oo,n_rel, &
                              ninitr,nq_nrixs,S_nrixs,secdd,secdo,secdq,secmd,secmm,secoo,secqq)
          if( Green_i ) &
            call MPI_RECV_all(l0_nrixs,Lmax_nrixs,mpinodee,mpirank,mpirank0,Multipole,n_oo,n_rel, &
                              ninitr,nq_nrixs,S_nrixs_m,secdd_m,secdo_m,secdq_m,secmd_m,secmm_m,secoo_m,secqq_m)

          if( Tddft_xanes .or. Optic_xanes ) &
            call MPI_Bcast_tddft(je,mpinodee,mpirank,nbseuil,nenerg,nenerg_tddft,nlmamax,nspinp,nspino,rof0,Taull_tdd)

          if( RIXS ) call MPI_Bcast_RIXS(Ampl_abs,je,mpinodee,mpirank,nenerg_RIXS,nlm_probe,nlmomax,nspinp)

       endif

        do ie_computer = 0,mpinodee-1

          ie = ( je - 1 ) * mpinodee + ie_computer + ie0

          if( ie > nenerg0 ) exit

          if( mpinodee > 1 ) call MPI_Bcast(Eimag(ie),1,MPI_REAL8,ie_computer,MPI_COMM_GATHER,mpierr)

          if( mpirank0 == 0 ) then

            call CPU_TIME(time)
            Time_loc(1) = real(time,db)

            if( icheck(19) > 1 ) &
              call Tau_wout(iabsorig(multi_run),Energ(ie),ie,ie_computer,lmax_dft(iprabs,ie_computer),mpinodee, &
                           n_multi_run,nlmmax,nomfich,nspinp,Taull_dft)

            if( icheck(19) > 0 ) &
              call Tau_writing(Eimag(ie),Energ(ie),iprabs,ie_computer,itypepr,lmax_dft,Lmax_probe,mpinodee,n_atom_proto,nlmmax, &
                               nspinp,ntype,numat,Taull_dft)

            if( RIXS .and. icheck(19) > 0 ) &
              call Ampl_writing(Ampl_dft,ie_computer,Lmaxso,mpinodes,nlm_rixs,nlmomax,nspinp)

            if( .not. Optic ) then

              if( Abs_in_Bulk_roughness ) then

                call Write_coabs(Abs_U_iso,Allsite,angxyz,axyz,Bragg_rgh_bulk_abs,.false.,Cartesian_tensor, &
                Core_resolved,Dafs,Dafs_bio,E_cut,Energ,Energphot,.false.,Epsii,Eseuil,Final_tddft,First_E,f_avantseuil, &
                Full_self_abs,Green_i,hkl_dafs,i_range,iabsorig(multi_run),icheck(21),ie,ie_computer,igr_bulk_z_abs, &
                Int_tens,isigpi,isymeq,jseuil,Length_abs,Length_rel,Length_rel_abs,ltypcal,Matper,Moyenne,mpinodee,multi_0, &
                Multipole,n_abs_rgh, &
                n_bulk_sup,n_multi_run,n_bulk_z,n_bulk_z_abs,n_bulk_z_max_abs,n_bulk_zc_abs,n_max,n_oo, &
                n_rel,n_tens_max,natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg,ninit1,ninitr,nomabs,nomfich,nomfich_cal_conv, &
                nomfich_s,nphi_dafs,npldafs,nphim,nplr,nplrm,nseuil,nspinp,numat_abs,nxanout,pdp,phdf0t_rgh_bulk,phdt_rgh_Bulk, &
                pol,poldafse,poldafss,sec_atom,secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, &
                secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Self_abs,Spherical_signal, &
                Spherical_tensor,Spinorbite,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_ref,Taux_eq,V0muf,Vecdafse,Vecdafss, &
                Vec,Volume_maille,Xan_atom)

                multi_0 = multi_0 + 1
              endif

              call Write_coabs(Abs_U_iso,Allsite,angxyz,axyz,Bragg_abs,Bulk_step,Cartesian_tensor, &
                Core_resolved,Dafs,Dafs_bio,E_cut,Energ,Energphot,.false.,Epsii,Eseuil,Final_tddft,First_E,f_avantseuil, &
                Full_self_abs,Green_i,hkl_dafs,i_range,iabsorig(multi_run),icheck(21),ie,ie_computer,igr_bulk_z_abs, &
                Int_tens,isigpi,isymeq,jseuil,Length_abs,Length_rel,Length_rel_abs,ltypcal,Matper,Moyenne,mpinodee,multi_0, &
                Multipole,n_abs_rgh, &
                n_bulk_sup,n_multi_run,n_bulk_z,n_bulk_z_abs,n_bulk_z_max_abs,n_bulk_zc_abs,n_max,n_oo, &
                n_rel,n_tens_max,natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg,ninit1,ninitr,nomabs,nomfich,nomfich_cal_conv, &
                nomfich_s,nphi_dafs,npldafs,nphim,nplr,nplrm,nseuil,nspinp,numat_abs,nxanout,pdp,phdf0t,phdt, &
                pol,poldafse,poldafss,sec_atom,secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, &
                secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Self_abs,Spherical_signal, &
                Spherical_tensor,Spinorbite,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_ref,Taux_eq,V0muf,Vecdafse,Vecdafss, &
                Vec,Volume_maille,Xan_atom)

              if( Abs_in_Bulk_roughness ) multi_0 = multi_0 - 1

              if( NRIXS ) call Write_nrixs(Abs_U_iso,All_nrixs,Allsite,axyz,Core_resolved,E_cut,Energ,Energphot,.false., &
                Epsii,Eseuil,f_avantseuil,Final_tddft,First_E,Green_i,i_range, &
                iabsorig(multi_run),icheck(21),ie,ie_computer,isymeq,jseuil,l0_nrixs,Lmax_nrixs,mpinodee,n_multi_run, &
                natomsym,nbseuil,nenerg,ninit1,ninitr,nomfich,nomfich_cal_conv(multi_run),nq_nrixs,nseuil,nspinp, &
                numat_abs,Orthmatt,q_nrixs,Rot_atom_abs,Rot_int,S_nrixs,S_nrixs_m,Spinorbite,Taux_eq,V0muf,Volume_maille)

              First_E = .false.
            endif

            call CPU_TIME(time)
            Time_loc(2) = real(time,db)
            Time_rout(13) = Time_rout(13) + Time_loc(2) - Time_loc(1)

            if( Density .and. i_range == 1 ) &
              call Write_Density(Energ,Density_comp,Full_atom,Harm_cubic,iaabsi,iaprotoi,icheck(28),ie,ie_computer,Int_statedens, &
                itypei,itypepr,lamstdens,lla_state,lla2_state,Lmaxat,mpinodee,n_atom_0,n_atom_ind,n_atom_proto,natome,nenerg, &
                nomfich_s,nonexc_g,nrato,nrm,nspin,nspinp,ntype,numat,Rato,Rmtsd,State_all_out,Statedens)

            if( COOP ) then
              if( Spinorbite ) then  ! because coop with real harmonics and spinorbit does not work
                call Write_coop(Coverlap,.true.,Energ(ie),.false.,iaprotoi,ie,ie_computer,index_coop,itypepr, &
                                       mpinodee,n_atom_proto,nab_coop,natome,nomfich_s,nspin,nspinp,ntype,numat)
              else
                call Write_coop(Coverlap,Density_comp,Energ(ie),Harm_cubic,iaprotoi,ie,ie_computer,index_coop,itypepr, &
                                       mpinodee,n_atom_proto,nab_coop,natome,nomfich_s,nspin,nspinp,ntype,numat)
              endif
            endif

            call CPU_TIME(time)
            Time_loc(3) = real(time,db)
            Time_rout(14) = Time_rout(14) + Time_loc(3) - Time_loc(2)

          endif

        end do

      end do boucle_energ_xanes   ! End of loop over energy

      if( Bulk_step .and. .not. Full_atom .and. Self_nonexc ) then
   ! One keeps the bulk atom potentials
        Vcato_bulk(:,:,iprabs_nonexc) = Vcato(:,:,iprabs_nonexc)
        do ipr = n_atom_proto - n_atom_proto_bulk + 1,n_atom_proto
          Vcato_bulk(:,:,ipr) = Vcato(:,:,ipr)
        end do
      endif

      deallocate( Atom_axe, distai, ia_eq, ia_eq_inv, ia_rep, iaprotoi, igroupi, iopsym_atom, is_eq, itypei, nb_eq )
      deallocate( nb_rpr, nb_rep_t, posi, rot_atom )
      deallocate( RadialIntegral, Tau_comp )

      if( .not. ( Green .or. Second_run .or. Extract ) ) then
        deallocate( gradvr )
        deallocate( Ylmso )
      endif
      deallocate( Ylmato )

      deallocate( Coverlap, drho_self, index_coop, Lmaxa, nlmsa, Statedens, Statedens_i )
      deallocate( iato, lato, iso, lso, mato, mso, nlmso0, nlmsa0 )
      deallocate( ibord, isbord, imoy, imoy_out, isrt, nbord, nbordf )
      deallocate( poidsa, poidso, poidsov, poidsov_out, rho, rhons, Vr, Vrato )

      if( .not. ( ( One_run .and. multi_run_e /= n_multi_run ) .or. Extract ) ) &
        deallocate( Excato, rs, rhoato, rsato, Vcato, Vh, Vhns, Vxc, Vxcato )

    end do boucle_i_range

    if( One_run .and. multi_run_e == n_multi_run ) deallocate( Taull_stk )

    deallocate( S_nrixs, S_nrixs_m )
    deallocate( secdd, secmd, secmm, secdq, secdo, secoo, secqq )
    deallocate( secdd_m, secmd_m, secmm_m, secdq_m, secdo_m, secoo_m, secqq_m )
    deallocate( sec_atom )
    deallocate( Energ, Eimag )
    deallocate( Ampl_dft, Taull_dft )

    if( mpirank0 == 0 ) deallocate( Int_statedens )

    deallocate( Repres_comp )

    if( mpirank0 == 0 ) then
      call CPU_TIME(time)
      tp_XANES_2 = real(time,db)
      Time_rout(19) = Time_rout(19) + tp_XANES_2 - tp_XANES_1 ! XANES
    endif

!----------------------------------------------------------------------------------------------------------------------------------------------

    if( RIXS ) then

      if( multi_run == 1 ) allocate( File_name_rixs(n_multi_run) )

      call main_RIXS(Ampl_abs,Coef_g,Core_resolved,dV0bdcF,E_cut,E_cut_imp,E_Fermi,E_cut_man,Eclie,Eneg,Energ_rixs, &
        Energ_s,Epsii,Epsii_ref,Epsii_ref_man,Eseuil,Estart,File_name_rixs(multi_run),FDM_comp_m,Full_potential,Gamma_hole, &
        Gamma_hole_man,Hubb_abs,Hubb_diag_abs,iabsorig(multi_run), &
        icheck(29),igreq,iprabs_nonexc,isymeq,ip0,ip_max,is_g,jseuil,lmoins1,lplus1,Lseuil,m_g,m_hubb, &
        Moment_conservation,multi_run, &
        Multipole,n_atom_proto,n_atom_uc,n_multi_run,n_oo,n_q_rixs,n_rel,n_theta_in,n_two_theta,natomsym,nbseuil, &
        nenerg_in_rixs,nenerg_s,neqm,ngamh,ninit1,ninit,ninitr,nlm_pot,nlm_probe,nlmomax,nomfich, &
        nrato_abs,nrm,ns_dipmag,nseuil,nspin,nspino,nspinp,numat_abs,Orthmati,Posn,Powder, &
        psii,q_rixs,rato_abs,Relativiste,Renorm,RIXS_fast,Rmtg_abs,Rmtsd_abs,Rot_atom_abs,Rot_int,Spin_channel,Spinorbite, &
        Step_loss,Theta_in,Two_theta,V_intmax,V_hubb,Vcato_abs,Vhbdc,VxcbdcF,Vxcato_abs,Ylm_comp)

      if( multi_run == n_multi_run ) then
        if( mpirank0 == 0 ) call Write_Int_rixs(File_name_rixs,n_multi_run,RIXS_core)
        deallocate( File_name_rixs )
      endif

      if( mpirank0 == 0 ) then
        call CPU_TIME(time)
        Time_loc(2) = real(time,db)
        Time_rout(20) = Time_rout(20) + Time_loc(2) - tp_XANES_2
      endif

    endif

    deallocate( Ampl_abs )

    if( Tddft .or. Optic ) then

      Atom_nonsph_loc = .false.

      if( Optic ) then

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(1) = real(time,db)
        endif

        if( .not. No_DFT ) call main_optic(Abs_U_iso,angxyz,Allsite,axyz,Bragg_abs,Cartesian_tensor,Classic_irreg, &
          Core_resolved,Dafs,Dafs_bio,dV0bdcF,E_cut,E_cut_imp,E_cut_man,Eclie,Eneg,Energ_s,Ephot_min, &
          Extract_ten,Eseuil,Full_potential,Full_self_abs,Green,hkl_dafs,Hubb_abs,Hubb_diag_abs,icheck, &
          iabsorig(multi_run),ip_max,ip0,isigpi,isymeq, &
          jseuil,ldip,Length_abs,Lmax_pot,Lmax_probe,Lmax_abs_t,Lmaxat0,lmaxfree,lmoins1,loct,lplus1,lqua, &
          Lseuil,ltypcal,m_hubb,Matper,Moyenne,mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq, &
          msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,multi_0,Multipole,n_bulk_sup,n_multi_run, &
          n_oo,n_rel,n_rout,n_tens_max,natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ninit1,ninit,ninitr,nlm_pot, &
          nlm_probe,nlmamax,nomabs,nomfich,nomfich_cal_conv,nomfich_s,nphi_dafs,nphim,npldafs,nplr,nplrm, &
          nrato_abs,nrm,nseuil,nspin,nspino,nspinp, &
          numat_abs,nxanout,pdp,phdf0t,phdt,pol,poldafse,poldafss,psii, &
          rato_abs,Relativiste,Renorm,Rmtg_abs,Rmtsd_abs,rot_atom_abs,Rot_int,Self_abs,Solsing_only,Spherical_signal, &
          Spherical_tensor,Spinorbite,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Taull_tdd,Taux_eq,Tddft,Time_rout,V_intmax, &
          V_hubb_abs,V0muf,Vcato_abs,vec,Vecdafse,Vecdafss,Vhbdc,Volume_maille,VxcbdcF,Vxcato_abs,Workf,Ylm_comp)

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(2) = real(time,db)
          Time_rout(20) = Time_rout(20) + Time_loc(2) - Time_loc(1)
        endif

        if( Core_resolved ) then
          ninit_out = ninitr
        else
          ninit_out = 1
        endif

        if( Tddft ) call main_tddft_optic(Abs_U_iso,alfpot,angxyz,Allsite,Atomic_scr,axyz,Bragg_abs, &
          Classic_irreg,coef_g,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio,Dipmag, &
          dV0bdcF,E_cut,E_cut_imp,E_cut_man,Eclie,Eneg, &
          Energ_s,Energphot,Ephot_min,Eseuil,f_avantseuil,Full_potential,Full_self_abs, &
          Gamma_tddft,hkl_dafs,Hubb_abs,Hubb_diag_abs,icheck,iabsorig(multi_run),is_g,isigpi,isymeq, &
          jseuil,Kern_fac,Kern_fast,ldip,Length_abs,Lmax_pot,Lmax_abs_t,lmoins1,loct,lplus1,lqua,Lseuil, &
          ltypcal,m_g,m_hubb,Magnetic,Matper,Moyenne,mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq, &
          msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,multi_0,Multipole, &
          n_bulk_sup,n_multi_run,n_oo,n_rel,n_rout,n_tens_max, &
          natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,nenerg_tddft,ninit1,ninit,ninit_out,ninitr,nlm_pot,nlmamax, &
          nomabs,nomfich,nomfich_cal_tddft_conv,nomfich_s, &
          nphi_dafs,nphim,npldafs,nplr,nplrm,nrato_abs,nrm,nseuil,nspin,nspino,nspinp, &
          numat_abs,nxanout,Octupole,pdp,phdf0t,phdt,pol,poldafse,poldafss, &
          psii,rato_abs,Relativiste,Renorm,rhoato_abs,Rmtg_abs,Rmtsd_abs, &
          rot_atom_abs,Rot_int,RPALF,rsato_abs,Self_abs,Solsing_only, &
          Spherical_signal,Spherical_tensor,Spinorbite,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Taull_tdd,Taux_eq,Time_rout, &
          V_intmax,V_hubb_abs,V0muf,Vcato_abs,Vec,Vecdafse,Vecdafss,Vhbdc,Volume_maille,VxcbdcF, &
          Vxcato_abs,Workf,Ylm_comp)

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(3) = real(time,db)
          Time_rout(21) = Time_rout(21) + Time_loc(3) - Time_loc(2)
        endif

      elseif( Tddft ) then

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(2) = real(time,db)
        endif

        ninit_out = 1

        call main_tddft(Abs_in_Bulk_roughness,Abs_U_iso,alfpot,All_nrixs,angxyz,Allsite,Atomic_scr,axyz,Bragg_abs, &
          Bragg_rgh_bulk_abs,Bulk_step,Classic_irreg,coef_g,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio,Delta_edge,Delta_Eseuil, &
          Dipmag,dV0bdcF,Dyn_eg,Dyn_g,E_cut,E_cut_imp,E_Fermi,E_cut_man,Ecent,Eclie,Elarg,Eneg, &
          Energ_s,Energphot,Epsii,Epsii_moy,Eseuil,Estart,f_avantseuil,Full_potential,Full_self_abs, &
          Gamma_hole,Gamma_hole_man,Gamma_max,Gamma_tddft,hkl_dafs,Hubb_abs,Hubb_diag_abs,icheck, &
          iabsorig(multi_run),igr_bulk_z_abs,iopsymc(25),is_g,isigpi,isymeq,jseuil,Kern_fac, &
          l0_nrixs,ldip,Length_abs,Length_rel,Length_rel_abs,Lmax_pot,Lmax_nrixs,Lmax_abs_t,Lmaxat0,lmaxfree,lmoins1,loct,lplus1, &
          lqua,Lseuil,ltypcal,m_g,m_hubb,Magnetic,Matper,Moyenne,mpinodes,mpirank,mpirank0,msymdd,msymddi, &
          msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,multi_0,multi_run,Multipole, &
          n_abs_rgh,n_bulk_sup,n_bulk_z,n_bulk_z_abs,n_bulk_z_max_abs,n_bulk_zc_abs,n_max,n_multi_run, &
          n_oo,n_rel,n_rout,n_tens_max, &
          natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,nenerg_tddft,ngamh,ninit1,ninit,ninit_out,ninitr,nlm_pot,nlmamax, &
          nomabs,nomfich,nomfich_cal_tddft_conv,nomfich_s, &
          nphi_dafs,nphim,npldafs,nplr,nplrm,nq_nrixs,nrato_abs,NRIXS,nrm,nseuil,nspin,nspino,nspinp, &
          numat_abs,nxanout,Octupole,Old_zero,Orthmatt,pdp,phdf0t,phdf0t_rgh_bulk,phdt,phdt_rgh_Bulk,pol,poldafse,poldafss, &
          psii,q_nrixs,Quadrupole,rato_abs,Relativiste,Renorm,rhoato_abs,Rmtg_abs,Rmtsd_abs, &
          rof0,Rot_atom_abs,Rot_int,RPALF,rsato_abs,rsbdc,Self_abs,Solsing_only,Spherical_signal,Spherical_tensor, &
          Spinorbite,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_ref,Taull_tdd,Taux_eq,Time_rout,V_intmax,V_hubb_abs,V0muf, &
          Vcato_abs,Vec,Vecdafse,Vecdafss,Vhbdc,Volume_maille,VxcbdcF, &
          Vxcato_abs,Workf,Ylm_comp)

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          Time_loc(3) = real(time,db)
          Time_rout(21) = Time_rout(21) + Time_loc(3) - Time_loc(2)
        endif

      endif

      deallocate( Taull_tdd )

    endif

    deallocate( rof0 )

    if( Self_cons_c .and. .not. Second_run .and. .not. Second_SCF ) deallocate( Energ_coh, Eimag_coh, posi_self  )

    if( .not. Extract ) then
      if( Bulk_step .and. .not. Full_atom .and. Self_nonexc ) then
   ! One keeps the bulk atom potentials
        chargat_init_bulk(iprabs_nonexc,:) = Chargat_init(iprabs_nonexc,:)
        chargat_self_bulk(iprabs_nonexc,:) = Chargat_self(iprabs_nonexc,:)
        rho_chg_bulk(:,:,iprabs_nonexc) = rho_chg(:,:,iprabs_nonexc)
        rho_self_bulk(:,:,:,iprabs_nonexc) = rho_self(:,:,:,iprabs_nonexc)
        rhoato_init_bulk(:,:,iprabs_nonexc) = rhoato_init(:,:,iprabs_nonexc)
        Vcato_init_bulk(:,iprabs_nonexc) = Vcato_init(:,iprabs_nonexc)
        do ipr = n_atom_proto - n_atom_proto_bulk + 1,n_atom_proto
          if( numat( itypepr(ipr) ) == numat_abs ) cycle
          chargat_init_bulk(ipr,:) = Chargat_init(ipr,:)
          chargat_self_bulk(ipr,:) = Chargat_self(ipr,:)
          rho_chg_bulk(:,:,ipr) = rho_chg(:,:,ipr)
          rho_self_bulk(:,:,:,ipr) = rho_self(:,:,:,ipr)
          rhoato_init_bulk(:,:,ipr) = rhoato_init(:,:,ipr)
          Vcato_init_bulk(:,ipr) = Vcato_init(:,ipr)
        end do
      endif
      if( .not. ( One_SCF .and. multi_run_e /= n_multi_run ) ) then
        deallocate( Chargat_init )
        deallocate( Chargat_self )
        deallocate( rho_chg, rho_self )
        deallocate( rhoato_init, Vcato_init )
        deallocate( dExc_ex_nex, drho_ex_nex, dvc_ex_nex, dv_ex_nex )
      endif
      deallocate( Exc_abs_i, Vc_abs_i, V_abs_i )
      if( multi_run_e == n_multi_run .and. Self_nonexc ) &
        deallocate( chargat_init_bulk, chargat_self_bulk, pop_orb_val_bulk, rho_chg_bulk, rho_self_bulk, rhoato_init_bulk, &
                    Vcato_bulk, Vcato_init_bulk )

    endif
    if( .not. ( Second_run .or. One_SCF .or. Extract ) ) deallocate( ia_eq_inv_self )

    deallocate( Hubb_diag, V_hubb_abs )
    deallocate( V_hubb )

    if( ( .not. One_run .or. multi_run_e == n_multi_run ) .and. .not. Extract ) then
      deallocate( xyz )
      deallocate( cgrad, clapl, ivois, isvois )
      deallocate( numia, rvol )
    endif

    deallocate( Bragg_abs, Bragg_rgh_bulk_abs, karact, Axe_atom_clu, dista )
    deallocate( Epsii )
    deallocate( iaproto, igroup, isymeq, itypep )
    deallocate( pos, Taux_eq )
    deallocate( poldafse, poldafss, phdt, Sum_Bragg_bulk_abs_f0, Sum_Bragg_nonabs_f )
    deallocate( Vecdafse, Vecdafss, ltypcal, nomabs, pdp, pol, vec )
    deallocate( igr_bulk_z, igr_bulk_z_abs, Length_rel, Length_rel_abs, n_bulk_zc, n_bulk_zc_abs )

    if( mpirank0 == 0 ) deallocate( Int_tens )

    deallocate( rato_abs, rhoato_abs, rsato_abs, Vcato_abs, Vxcato_abs )

    multi_0 = multi_0 + max( 1, n_Bulk_z_abs )
    if( Bulk_roughness > 0 .and. Bulk_step ) multi_0 = multi_0 + 1
    if( Bulk_step ) First_bulk_step = .false.

  end do boucle_multi ! End of loop over non-equivalent atoms

  if( One_SCF ) deallocate( ia_eq_inv_self )
  deallocate( coef_g, is_g, m_g )
  if( .not. Extract ) deallocate( Eimag_s, Energ_s )
  deallocate( msymoo, msymooi )

  return

  110 format(/' ----------',110('-'),//' Absorbing atom',i4,' begining')
  140 format(/1x,120('-')//,' Last cycle, XANES calculation')
  150 format(/1x,120('-')//,' Begining of cycle',i3)
  155 format(/1x,120('-')//,' Begining of cycle, calculation of excited absorbing atom Fermi energy')
  160 format(/' E_cut =',f11.5,' eV, E_Fermi =',f11.5,' eV')
  170 format(///' E_kinetic =',f7.3,' eV < 0.',/ ' Start the calculation at higher energy !'///)
  180 format(///' or E_kinetic_ext =',f7.3,' eV < 0.',/ ' Start the calculation at higher energy !'///)
  190 format(/' ---- Sphere -------',100('-'))
end

!*************************************************************************************************************************

subroutine Extract_write_coabs(Abs_in_Bulk_roughness,Abs_U_iso,Allsite,Ang_rotsup,angxyz,axyz,Bragg_abs,Bragg_rgh_bulk_abs, &
         Bulk_step,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio,E_cut,Energ_s,Energphot,Epsii,Eseuil,Final_tddft,f_avantseuil, &
         Full_self_abs,hkl_dafs,i_range,iabsorig,icheck,igr_bulk_z_abs,Int_tens,isigpi,isymeq, &
         isymext,jseuil,Length_abs,Length_rel,Length_rel_abs,ltypcal,Matper,Moyenne,mpinodee,multi_0,multi_run_e,Multipole, &
         n_abs_rgh,n_bulk_sup,n_bulk_z,n_bulk_z_abs,n_bulk_z_max_abs,n_bulk_zc_abs,n_max,n_multi_run,n_oo, &
         n_rel,n_tens_max,natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ninit1,ninitr,nomabs,nomfich,nomfich_cal_conv, &
         nom_fich_extract,nomfich_s,nphi_dafs,nphim,npldafs,nplr,nplrm,nseuil,nspinp,numat_abs,nxanout,pdp, &
         phdf0t,phdf0t_rgh_bulk,phdt,phdt_rgh_Bulk,pol,poldafse,poldafss,Self_abs,Spherical_signal,Spherical_tensor,Spinorbite, &
         Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_ref,Taux_eq,Tddft,V0muf,Vecdafse,Vecdafss,Vec,Volume_maille,Xan_atom)

  use declarations
  implicit none

  integer:: i_range, iabsorig, icheck, ie, ie_computer, isymext, jseuil, mpinodee, multi_0, multi_run_e, n_abs_rgh, n_bulk_sup, &
    n_bulk_z, n_bulk_z_abs, n_bulk_z_max_abs, n_max, n_multi_run, n_oo, n_rel, n_tens_max, natomsym, nbseuil, &
    ncolm, ncolr, ncolt, ninit1, ninitr, nenerg_s, nseuil, nspinp, nphim, npldafs, nplr, nplrm, numat_abs, nxanout

  integer, dimension(2,npldafs):: isigpi
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(natomsym):: isymeq
  integer, dimension(n_bulk_z_abs):: n_bulk_zc_abs
  integer, dimension(n_bulk_z_max_abs,n_bulk_z_abs):: igr_bulk_z_abs

  character(len=Length_name):: nomfich, nom_fich_extract, nomfich_s
  character(len=5), dimension(nplrm):: ltypcal
  character(len=Length_word), dimension(ncolm):: nomabs
  character(len=Length_name), dimension(n_multi_run+n_bulk_sup+n_abs_rgh):: nomfich_cal_conv

  complex(kind=db):: f_avantseuil
  complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdf0t_rgh_bulk, phdt_rgh_Bulk
  complex(kind=db), dimension(npldafs,nphim,n_max):: phdt
  complex(kind=db), dimension(npldafs,n_max):: Sum_Bragg_bulk_abs_f0
  complex(kind=db), dimension(npldafs,n_bulk_z):: Sum_Bragg_nonabs_f
  complex(kind=db), dimension(natomsym,npldafs):: Bragg_abs, Bragg_rgh_bulk_abs
  complex(kind=db), dimension(3,nplrm):: pol
  complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss
  complex(kind=db), dimension(3,3,n_rel,ninitr,0:0):: secdd, secdd_m
  complex(kind=db), dimension(3,3,ninitr,0:0):: secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitr,0:0):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitr,0:0):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:0):: secoo, secoo_m

  logical:: Abs_in_Bulk_roughness, Allsite, Bulk_step, Cartesian_tensor, Core_resolved, Dafs, Dafs_bio, &
     Energphot, Final_tddft, First_E, Full_self_abs, Green_int, Matper, Moyenne, Self_abs, Spherical_signal, &
     Spherical_tensor, Spinorbite, Tddft, Tensor_rot, Xan_atom

  logical, dimension(10):: Multipole

  real(kind=db):: Abs_U_iso, E_cut, Surface_ref, V0muf, Volume_maille
  real(kind=db), dimension(3):: Ang_rotsup, angxyz, axyz
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitr):: Epsii, Sec_atom
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(n_tens_max*ninitr,0:natomsym):: Int_tens
  real(kind=db), dimension(3,3):: Rotsup
  real(kind=db), dimension(nplrm,2):: pdp
  real(kind=db), dimension(natomsym):: Taux_eq
  real(kind=db), dimension(3,nplrm):: vec
  real(kind=db), dimension(npldafs):: Length_abs
  real(kind=db), dimension(n_bulk_z):: Length_rel
  real(kind=db), dimension(n_bulk_z_abs):: Length_rel_abs
  real(kind=db), dimension(3,npldafs):: hkl_dafs
  real(kind=db), dimension(3,npldafs,nphim):: Vecdafse, Vecdafss

! Sec_atom is not extracted from the bav file
  Sec_atom(:) = 0._db

  ie_computer = 0
  do ie = 1,nenerg_s
    First_E = ie == 1
    call extract_coabs(Ang_rotsup,Core_resolved,Green_int,icheck,ie,isymext,multi_run_e,Multipole, &
            n_oo,n_rel,nenerg_s,ninit1,ninitr,nom_fich_extract,Rotsup,secdd,secdd_m,secdo,secdo_m,secdq, &
            secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Tensor_rot,Tddft)

    if( Abs_in_Bulk_roughness ) then

      call Write_coabs(Abs_U_iso,Allsite,angxyz,axyz,Bragg_rgh_bulk_abs,.false.,Cartesian_tensor, &
            Core_resolved,Dafs,Dafs_bio,E_cut,Energ_s,Energphot,.false.,Epsii,Eseuil,Final_tddft,First_E, &
            f_avantseuil,Full_self_abs,Green_int,hkl_dafs,i_range,iabsorig,icheck,ie,ie_computer,igr_bulk_z_abs, &
            Int_tens,isigpi,isymeq,jseuil,Length_abs,Length_rel,Length_rel_abs,ltypcal,Matper,Moyenne,mpinodee, &
            multi_0,Multipole,n_abs_rgh,n_bulk_sup,n_multi_run, &
            n_bulk_z,n_bulk_z_abs,n_bulk_z_max_abs,n_bulk_zc_abs,n_max,n_oo,n_rel,n_tens_max,natomsym,nbseuil, &
            ncolm,ncolr,ncolt,nenerg_s,ninit1,ninitr,nomabs,nomfich,nomfich_cal_conv,nomfich_s,nphi_dafs,npldafs, &
            nphim,nplr,nplrm,nseuil,nspinp,numat_abs,nxanout,pdp,phdf0t_rgh_bulk,phdt_rgh_Bulk,pol,poldafse,poldafss, &
            sec_atom,secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, &
            secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Self_abs,Spherical_signal, &
            Spherical_tensor,Spinorbite,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_ref,Taux_eq,V0muf,Vecdafse,Vecdafss, &
            Vec,Volume_maille,Xan_atom)

      multi_0 = multi_0 + 1
    endif

    call Write_coabs(Abs_U_iso,Allsite,angxyz,axyz,Bragg_abs,Bulk_step,Cartesian_tensor, &
            Core_resolved,Dafs,Dafs_bio,E_cut,Energ_s,Energphot,.true.,Epsii,Eseuil,Final_tddft,First_E, &
            f_avantseuil,Full_self_abs,Green_int,hkl_dafs,i_range,iabsorig,icheck,ie,ie_computer,igr_bulk_z_abs, &
            Int_tens,isigpi,isymeq,jseuil,Length_abs,Length_rel,Length_rel_abs,ltypcal,Matper,Moyenne,mpinodee, &
            multi_0,Multipole,n_abs_rgh,n_bulk_sup,n_multi_run, &
            n_bulk_z,n_bulk_z_abs,n_bulk_z_max_abs,n_bulk_zc_abs,n_max,n_oo,n_rel,n_tens_max,natomsym,nbseuil, &
            ncolm,ncolr,ncolt,nenerg_s,ninit1,ninitr,nomabs,nomfich,nomfich_cal_conv,nomfich_s,nphi_dafs,npldafs, &
            nphim,nplr,nplrm,nseuil,nspinp,numat_abs,nxanout,pdp,phdf0t,phdt,pol,poldafse,poldafss, &
            sec_atom,secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, &
            secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Self_abs,Spherical_signal, &
            Spherical_tensor,Spinorbite,Sum_Bragg_bulk_abs_f0,Sum_Bragg_nonabs_f,Surface_ref,Taux_eq,V0muf,Vecdafse,Vecdafss, &
            Vec,Volume_maille,Xan_atom)

    if( Abs_in_Bulk_roughness ) multi_0 = multi_0 - 1

  end do

  return
end

!***********************************************************************

subroutine MPI_RECV_Tau(lmax_dft,mpinodes,mpirank,mpirank0,n_atom_proto,nlmmax,nspinp,Taull_dft)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: idim, ie, mpirank, mpirank0, n_atom_proto, nlmmax, nspinp, mpirank_in_mumps_group
  integer, dimension(0:n_atom_proto,0:mpinodes-1):: lmax_dft

  complex(kind=db), dimension(nlmmax,nspinp,nlmmax,nspinp,0:mpinodes-1):: Taull_dft

  real(kind=db), dimension(nlmmax,nspinp,nlmmax,nspinp,0:mpinodes-1):: Taull_dft_i, Taull_dft_r

  call MPI_Comm_Rank(MPI_COMM_MUMPS,mpirank_in_mumps_group,mpierr)
  idim = n_atom_proto + 1
  if( mpirank0 == 0 ) then
   call MPI_GATHER(MPI_IN_PLACE,idim,MPI_INTEGER,lmax_dft,idim,MPI_INTEGER,0,MPI_COMM_GATHER,mpierr)
  elseif( mpirank_in_mumps_group == 0 ) then
	  call MPI_GATHER(lmax_dft(0,mpirank),idim,MPI_INTEGER,lmax_dft,idim,MPI_INTEGER,0,MPI_COMM_GATHER,mpierr)
  endif

  Taull_dft_r(:,:,:,:,mpirank) = real( Taull_dft(:,:,:,:,mpirank),db )
  Taull_dft_i(:,:,:,:,mpirank) = aimag( Taull_dft(:,:,:,:,mpirank) )

  idim = ( nlmmax * nspinp )**2

  if( mpirank0 == 0 ) then
   call MPI_GATHER(MPI_IN_PLACE,idim,MPI_REAL8,Taull_dft_r,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
   call MPI_GATHER(MPI_IN_PLACE,idim,MPI_REAL8,Taull_dft_i,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  elseif( mpirank_in_mumps_group == 0 ) then
	  call MPI_GATHER(Taull_dft_r(1,1,1,1,mpirank),idim,MPI_REAL8,Taull_dft_r,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	  call MPI_GATHER(Taull_dft_i(1,1,1,1,mpirank),idim,MPI_REAL8,Taull_dft_i,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  endif



  if( mpirank0 == 0 ) then
    do ie = 1,mpinodes-1
      Taull_dft(:,:,:,:,ie) = cmplx( Taull_dft_r(:,:,:,:,ie), Taull_dft_i(:,:,:,:,ie),db )
    end do
  end if

  return
end

!***********************************************************************

subroutine MPI_RECV_all(l0_nrixs,Lmax_nrixs,mpinodes,mpirank,mpirank0,Multipole,n_oo,n_rel,ninitr, &
                  nq_nrixs,S_nrixs,secdd,secdo,secdq,secmd,secmm,secoo,secqq)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: idim0, idim1, idim2, idim3, idim4, idim5, mpirank_in_mumps_group, mpierr, mpinodes, mpirank, mpirank0, n_oo, n_rel, &
            ninitr, rang

  complex(kind=db), dimension(3,3,n_rel,ninitr,0:mpinodes-1):: secdd
  complex(kind=db), dimension(3,3,ninitr,0:mpinodes-1):: secmd, secmm
  complex(kind=db), dimension(3,3,3,ninitr,0:mpinodes-1):: secdq
  complex(kind=db), dimension(3,3,3,3,ninitr,0:mpinodes-1):: secdo, secqq
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:mpinodes-1):: secoo
  complex(kind=db), dimension(nq_nrixs,(Lmax_nrixs+1)**2,(Lmax_nrixs+1)**2,ninitr,0:mpinodes-1):: S_nrixs

  logical:: E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, M1M1

  logical, dimension(10):: Multipole

  real(kind=db), dimension(nq_nrixs,(Lmax_nrixs+1)**2,(Lmax_nrixs+1)**2,ninitr,0:mpinodes-1):: S_nrixs_ei, S_nrixs_er
  real(kind=db), dimension(3,3,n_rel,ninitr,0:mpinodes-1):: secdd_er, secdd_ei
  real(kind=db), dimension(3,3,ninitr,0:mpinodes-1):: secmd_er, secmm_er, secmd_ei, secmm_ei
  real(kind=db), dimension(3,3,3,ninitr,0:mpinodes-1):: secdq_er, secdq_ei
  real(kind=db), dimension(3,3,3,3,ninitr,0:mpinodes-1):: secdo_ei, secdo_er, secqq_ei, secqq_er
  real(kind=db), dimension(3,n_oo,3,n_oo,ninitr,0:mpinodes-1):: secoo_ei, secoo_er

  call MPI_Comm_Rank(MPI_COMM_MUMPS,mpirank_in_mumps_group,mpierr)
  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  idim0 = 9 * n_rel * ninitr
  idim1 = 9 * ninitr
  idim2 = 3 * idim1
  idim3 = 3 * idim2
  idim4 = ( n_oo**2 ) * idim2
  idim5 = nq_nrixs * ninitr * (Lmax_nrixs+1)**4

! Copies to store the real and imaginary part
! The processor index is the last index

  if( E1E1 )  then
    secdd_er(:,:,:,:,mpirank) = real( secdd(:,:,:,:,mpirank),db )
    secdd_ei(:,:,:,:,mpirank) = aimag( secdd(:,:,:,:,mpirank) )
  endif

  if( E1M1 )  then
    secmd_er(:,:,:,mpirank) = real( secmd(:,:,:,mpirank),db )
    secmd_ei(:,:,:,mpirank) = aimag( secmd(:,:,:,mpirank) )
  endif

  if( M1M1 )  then
    secmm_er(:,:,:,mpirank) = real( secmm(:,:,:,mpirank),db )
    secmm_ei(:,:,:,mpirank) = aimag( secmm(:,:,:,mpirank) )
  endif

  if( E1E2 )  then
    secdq_er(:,:,:,:,mpirank) = real( secdq(:,:,:,:,mpirank),db)
    secdq_ei(:,:,:,:,mpirank) = aimag( secdq(:,:,:,:,mpirank) )
  endif

  if( E1E3 )  then
    secdo_er(:,:,:,:,:,mpirank) = real( secdo(:,:,:,:,:,mpirank),db )
    secdo_ei(:,:,:,:,:,mpirank) = aimag( secdo(:,:,:,:,:,mpirank) )
  endif

  if( E3E3 )  then
    secoo_er(:,:,:,:,:,mpirank) = real( secoo(:,:,:,:,:,mpirank),db )
    secoo_ei(:,:,:,:,:,mpirank) = aimag( secoo(:,:,:,:,:,mpirank) )
  endif

  if( E2E2 )  then
    secqq_er(:,:,:,:,:,mpirank) = real( secqq(:,:,:,:,:,mpirank),db )
    secqq_ei(:,:,:,:,:,mpirank) = aimag( secqq(:,:,:,:,:,mpirank) )
  endif

  S_nrixs_er(:,:,:,:,mpirank) = real( S_nrixs(:,:,:,:,mpirank),db)
  S_nrixs_ei(:,:,:,:,mpirank) = aimag( S_nrixs(:,:,:,:,mpirank) )

  if( mpirank0 == 0 ) then
    if( E1E1 )  then
      call MPI_GATHER(MPI_IN_PLACE,idim0,MPI_REAL8,secdd_er,idim0,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      call MPI_GATHER(MPI_IN_PLACE,idim0,MPI_REAL8,secdd_ei,idim0,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
    endif
    if( E1M1 )  then
      call MPI_GATHER(MPI_IN_PLACE,idim1,MPI_REAL8,secmd_er,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      call MPI_GATHER(MPI_IN_PLACE,idim1,MPI_REAL8,secmd_ei,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
    end if
    if( M1M1 )  then
      call MPI_GATHER(MPI_IN_PLACE,idim1,MPI_REAL8,secmm_er,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      call MPI_GATHER(MPI_IN_PLACE,idim1,MPI_REAL8,secmm_ei,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
    end if
    if( E1E2 )  then
      call MPI_GATHER(MPI_IN_PLACE,idim2,MPI_REAL8,secdq_er,idim2,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      call MPI_GATHER(MPI_IN_PLACE,idim2,MPI_REAL8,secdq_ei,idim2,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
    end if
    if( E2E2 )  then
      call MPI_GATHER(MPI_IN_PLACE,idim3,MPI_REAL8,secqq_er,idim3,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      call MPI_GATHER(MPI_IN_PLACE,idim3,MPI_REAL8,secqq_ei,idim3,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
    end if
    if( E1E3 )  then
      call MPI_GATHER(MPI_IN_PLACE,idim3,MPI_REAL8,secdo_er,idim3,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      call MPI_GATHER(MPI_IN_PLACE,idim3,MPI_REAL8,secdo_ei,idim3,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
    end if
    if( E3E3 )  then
      call MPI_GATHER(MPI_IN_PLACE,idim3,MPI_REAL8,secoo_er,idim4,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      call MPI_GATHER(MPI_IN_PLACE,idim3,MPI_REAL8,secoo_ei,idim4,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
    end if

    if( idim5 > 0 ) then
      call MPI_GATHER(MPI_IN_PLACE,idim5,MPI_REAL8,S_nrixs_er,idim5,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      call MPI_GATHER(MPI_IN_PLACE,idim5,MPI_REAL8,S_nrixs_ei,idim5,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
    endif

  elseif( mpirank_in_mumps_group == 0 ) then

	  if( E1E1 ) then
	    call MPI_GATHER(secdd_er(1,1,1,1,mpirank),idim0,MPI_REAL8,secdd_er,idim0,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	    call MPI_GATHER(secdd_ei(1,1,1,1,mpirank),idim0,MPI_REAL8,secdd_ei,idim0,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	  endif
	  if( E1M1 ) then
	    call MPI_GATHER(secmd_er(1,1,1,mpirank),idim1,MPI_REAL8,secmd_er,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	    call MPI_GATHER(secmd_ei(1,1,1,mpirank),idim1,MPI_REAL8,secmd_ei,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	  end if

	  if( M1M1 ) then
	    call MPI_GATHER(secmm_er(1,1,1,mpirank),idim1,MPI_REAL8,secmm_er,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	    call MPI_GATHER(secmm_ei(1,1,1,mpirank),idim1,MPI_REAL8,secmm_ei,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	  end if

	  if( E1E2 ) then
	    call MPI_GATHER(secdq_er(1,1,1,1,mpirank),idim2, MPI_REAL8,secdq_er,idim2,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	    call MPI_GATHER(secdq_ei(1,1,1,1,mpirank),idim2, MPI_REAL8,secdq_ei,idim2,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	  end if

	  if( E2E2 ) then
	    call MPI_GATHER(secqq_er(1,1,1,1,1,mpirank),idim3, MPI_REAL8,secqq_er,idim3,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	    call MPI_GATHER(secqq_ei(1,1,1,1,1,mpirank),idim3, MPI_REAL8,secqq_ei,idim3,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	  end if

	  if( E1E3 ) then
	    call MPI_GATHER(secdo_er(1,1,1,1,1,mpirank),idim3, MPI_REAL8,secdo_er,idim3,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	    call MPI_GATHER(secdo_ei(1,1,1,1,1,mpirank),idim3, MPI_REAL8,secdo_ei,idim3,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	  end if

	  if( E3E3 ) then
	    call MPI_GATHER(secoo_er(1,1,1,1,1,mpirank),idim4, MPI_REAL8,secoo_er,idim4,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	    call MPI_GATHER(secoo_ei(1,1,1,1,1,mpirank),idim4, MPI_REAL8,secoo_ei,idim4,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	  end if

      if( idim5 > 0 ) then
	    call MPI_GATHER(S_nrixs_er(1,1,1,1,mpirank),idim5,MPI_REAL8,S_nrixs_er,idim5,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	    call MPI_GATHER(S_nrixs_ei(1,1,1,1,mpirank),idim5,MPI_REAL8,S_nrixs_ei,idim5,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      endif

  end if

! One rebuilds the secXX; Now the central processor has the results from all the others
  if( mpirank0 == 0 ) then
    do rang = 1,mpinodes-1
      if( E1E1 ) secdd(:,:,:,:,rang) = cmplx( secdd_er(:,:,:,:,rang), secdd_ei(:,:,:,:,rang),db )
      if( E1M1 ) secmd(:,:,:,rang) = cmplx( secmd_er(:,:,:,rang), secmd_ei(:,:,:,rang),db )
      if( M1M1 ) secmm(:,:,:,rang) = cmplx( secmm_er(:,:,:,rang), secmm_ei(:,:,:,rang),db )
      if( E1E2 ) secdq(:,:,:,:,rang) = cmplx( secdq_er(:,:,:,:,rang), secdq_ei(:,:,:,:,rang),db )
      if( E2E2 ) secqq(:,:,:,:,:,rang) = cmplx( secqq_er(:,:,:,:,:,rang), secqq_ei(:,:,:,:,:,rang),db )
      if( E1E3 ) secdo(:,:,:,:,:,rang) = cmplx( secdo_er(:,:,:,:,:,rang), secdo_ei(:,:,:,:,:,rang),db )
      if( E3E3 ) secdo(:,:,:,:,:,rang) = cmplx( secoo_er(:,:,:,:,:,rang), secoo_ei(:,:,:,:,:,rang),db )
      S_nrixs(:,:,:,:,rang) = cmplx( S_nrixs_er(:,:,:,:,rang), S_nrixs_ei(:,:,:,:,rang),db )
    end do
  end if

  return
end

!***********************************************************************

subroutine MPI_RECV_statedens(lla2_state,Lmaxat,mpinodes,mpirank,mpirank0,n_atom_0, &
                              n_atom_ind,n_atom_proto,nspinp,Statedens,Statedens_i)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: mpirank_in_mumps_group
  integer, dimension(0:n_atom_proto):: Lmaxat
  integer, dimension(0:n_atom_proto,0:mpinodes-1):: lmaxat_e

  real(kind=db), dimension(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1):: Statedens, Statedens_i

  call MPI_Comm_Rank(MPI_COMM_MUMPS,mpirank_in_mumps_group,mpierr)
  lmaxat_e(:,mpirank) = Lmaxat(:)
  idim = n_atom_proto + 1

! MPI_GATHER: le choix lorsque tous les ordinateurs envoyent le meme nombre d'elements a l'ordinateur central

  if( mpirank0 == 0 ) then
    call MPI_GATHER(MPI_IN_PLACE,idim,MPI_INTEGER,lmaxat_e,idim,MPI_INTEGER,0,MPI_COMM_GATHER,mpierr)
  elseif ( mpirank_in_mumps_group == 0 ) then
    call MPI_GATHER(lmaxat_e(0,mpirank),idim,MPI_INTEGER,lmaxat_e,idim,MPI_INTEGER,0,MPI_COMM_GATHER,mpierr)
  end if

! It is the Lmax corresponding to the highest energy, which are the biggest

  if( mpirank0 == 0 ) Lmaxat(:) = lmaxat_e(:,mpinodes-1)

  idim = ( n_atom_ind - n_atom_0 + 1 ) * (lla2_state * nspinp)**2

  if( mpirank0 == 0 ) then
    call MPI_GATHER(MPI_IN_PLACE,idim,MPI_REAL8,Statedens,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  elseif( mpirank_in_mumps_group == 0 ) then
    call MPI_GATHER(Statedens(1,1,1,1,n_atom_0,mpirank),idim,MPI_REAL8,Statedens,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  end if

  if( mpirank0 == 0 ) then
    call MPI_GATHER(MPI_IN_PLACE,idim,MPI_REAL8,Statedens_i,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  elseif( mpirank_in_mumps_group == 0 ) then
    call MPI_GATHER(Statedens_i(1,1,1,1,n_atom_0,mpirank),idim, MPI_REAL8,Statedens_i,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  end if

  return
end

!***********************************************************************

subroutine MPI_RECV_COOP(mpinodes,mpirank,mpirank0,nab_coop,nspinp,Coverlap)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: mpirank_in_mumps_group
  integer:: idim, mpinodes, mpirank, mpirank0, nab_coop, nspinp

  real(kind=db), dimension(16,16,nspinp,nab_coop,0:mpinodes-1):: Coverlap

  call MPI_Comm_Rank(MPI_COMM_MUMPS,mpirank_in_mumps_group,mpierr)

! MPI_GATHER: le choix lorsque tous les ordinateurs envoyent le meme nombre d'elements a l'ordinateur central

  idim = nab_coop * nspinp * (16)**2

  if( mpirank0 == 0 ) then
    call MPI_GATHER(MPI_IN_PLACE,idim,MPI_REAL8,Coverlap,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  elseif( mpirank_in_mumps_group == 0 ) then
    call MPI_GATHER(Coverlap(1,1,1,1,mpirank),idim,MPI_REAL8,Coverlap,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  end if

  return
end

!***********************************************************************
! IMPORTANT: depending on the operating system, the executable may run out
!     of virtual memory as the dummy arrays (such as drho_self) and the temporary ones are created in
!     the stack (i.e. static memory), which is subject to a limited available space.
!     Whereas a sequential application would normally return an error message,
!     a parallel MPI one would crash without any indication.
!     In order to avoid such problems it is advisable to force your compiler use
!     dynamic allocation even for temporary arrays, at least for those whose size
!     exceeds a certain limit that you may indicate at compilation time.
!     Should you use this trick, make sure that you compile and run the program
!     on the very same machine. This operation might have a price to pay in terms of performance.

!     Linux Intel compiler for Itanium based applications: -heap-array[:size]

subroutine MPI_RECV_self(drho_self,nlm_pot,mpinodes,mpirank,mpirank0,n_atom_0_self,n_atom_ind_self, &
                         nrm_self,nspinp)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: idim1, lm, mpirank_in_mumps_group, mpinodes, mpirank, mpirank0, n_atom_0_self, n_atom_ind_self, nlm_pot, nrm_self, &
            nspinp, rang

  real(kind=db), dimension(nrm_self,nlm_pot,nspinp,n_atom_0_self:n_atom_ind_self,0:mpinodes-1):: drho_self
  real(kind=db), dimension(nrm_self,nspinp,0:mpinodes-1):: drho_self_e

  call MPI_Comm_Rank(MPI_COMM_MUMPS,mpirank_in_mumps_group,mpierr)
  idim1 = nrm_self * nspinp

  do lm = 1,nlm_pot
    do ia = n_atom_0_self,n_atom_ind_self

      drho_self_e(:,:,mpirank) = drho_self(:,lm,:,ia,mpirank)

! Barrier is here very important: without it, it can be that we use the same
!    buffer rhov_self_eX by 2 different process with different iapr

      if( mpirank0 == 0 ) then
        call MPI_GATHER(MPI_IN_PLACE,idim1,MPI_REAL8,drho_self_e,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      elseif( mpirank_in_mumps_group == 0 ) then
        call MPI_GATHER(drho_self_e(1,1,mpirank),idim1,MPI_REAL8,drho_self_e,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      end if

      if( mpirank0 == 0 ) then
        do rang = 1,mpinodes-1
          drho_self(:,lm,:,ia,rang) = drho_self_e(:,:,rang)
        end do
      end if

! This barrier is important, because drho_self_e eis used as
! buffer and must not be filled by 2 different atom ia

    end do
  end do

  return
end

!***********************************************************************

subroutine MPI_Bcast_Hubb(Hubb_diag,m_hubb,mpirank,n_atom_0_self,n_atom_ind_self,nspinp,V_hubb)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: m_hubb, mpirank, n_atom_0_self, n_atom_ind_self, ndim, nspinp

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb

  logical, dimension(n_atom_0_self:n_atom_ind_self):: Hubb_diag

  real(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self):: V_hubb_i, V_hubb_r

  ndim = ( n_atom_ind_self - n_atom_0_self + 1 ) * nspinp**2 * ( 2 * m_hubb + 1 )**2

  if( mpirank == 0 ) then
    V_hubb_r(:,:,:,:,:) = real( V_hubb(:,:,:,:,:), db )
    V_hubb_i(:,:,:,:,:) = aimag( V_hubb(:,:,:,:,:) )
  endif

  call MPI_Bcast(V_hubb_r,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(V_hubb_i,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)



  if( mpirank /= 0 ) V_hubb(:,:,:,:,:) = cmplx( V_hubb_r(:,:,:,:,:), V_hubb_i(:,:,:,:,:), db )

  ndim = n_atom_ind_self - n_atom_0_self + 1

  call MPI_Bcast(Hubb_diag,ndim,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)

  return
end

!***********************************************************************

subroutine MPI_Bcast_tddft(je,mpinodes,mpirank,nbseuil,nenerg,nenerg_tddft,nlmamax,nspinp,nspino,rof0,Taull_tdd)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: ie, ie_computer, je, mpinodes, mpirank, ndim, nbseuil, nenerg, nenerg_tddft, nlmamax, nspinp, nspino

  complex(kind=db), dimension(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil):: rof0
  complex(kind=db), dimension(nenerg,nlmamax,nspinp,nlmamax,nspinp):: Taull_tdd

  real(kind=db), dimension(nlmamax,nspinp,nlmamax,nspinp):: Taull_tdd_i, Taull_tdd_r
  real(kind=db), dimension(nlmamax,nspinp,nspino,nbseuil):: rof0_i, rof0_r

  ndim = ( nlmamax * nspinp )**2

  do ie_computer = 0,mpinodes-1

    ie = ( je - 1 ) * mpinodes + ie_computer + 1
    if( ie > nenerg ) exit

    if( ie_computer == mpirank ) then
      Taull_tdd_r(:,:,:,:) = real( Taull_tdd(ie,:,:,:,:), db )
      Taull_tdd_i(:,:,:,:) = aimag( Taull_tdd(ie,:,:,:,:) )
    endif

    call MPI_Bcast(Taull_tdd_r,ndim,MPI_REAL8,ie_computer,MPI_COMM_GATHER,mpierr)
    call MPI_Bcast(Taull_tdd_i,ndim,MPI_REAL8,ie_computer,MPI_COMM_GATHER,mpierr)



    if( ie_computer /= mpirank ) Taull_tdd(ie,:,:,:,:) = cmplx( Taull_tdd_r(:,:,:,:), Taull_tdd_i(:,:,:,:), db )

  end do

  ndim = nlmamax * nspinp * nspino * nbseuil

  if( nenerg_tddft > 0 ) then

    do ie_computer = 0,mpinodes-1

      ie = ( je - 1 ) * mpinodes + ie_computer + 1
      if( ie > nenerg_tddft ) exit

      if( ie_computer == mpirank ) then
        rof0_r(:,:,:,:) = real( rof0(ie,:,:,:,:),db )
        rof0_i(:,:,:,:) = aimag( rof0(ie,:,:,:,:) )
      endif

      call MPI_Bcast(rof0_r,ndim,MPI_REAL8,ie_computer,MPI_COMM_GATHER,mpierr)

      call MPI_Bcast(rof0_i,ndim,MPI_REAL8,ie_computer,MPI_COMM_GATHER,mpierr)



      if( ie_computer /= mpirank ) rof0(ie,:,:,:,:) = cmplx( rof0_r(:,:,:,:), rof0_i(:,:,:,:), db )

    end do

  endif

  return
end

!***********************************************************************

subroutine MPI_Bcast_RIXS(Ampl_Abs,je,mpinodee,mpirank,nenerg_RIXS,nlm_probe,nlmomax,nspinp)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: ie, ie_computer, je, mpinodes, mpirank, ndim, nenerg_RIXS, nlm_probe, nlmomax, nspinp

  complex(kind=db), dimension(nenerg_RIXS,nlm_probe,nspinp,nlmomax,nspinp):: Ampl_abs

  real(kind=db), dimension(nlm_probe,nspinp,nlmomax,nspinp):: Ampl_abs_i, Ampl_abs_r

  ndim = nlm_probe * nspinp**2 * nlmomax

  do ie_computer = 0,mpinodes-1

    ie = ( je - 1 ) * mpinodes + ie_computer + 1
    if( ie > nenerg ) exit

    if( ie_computer == mpirank ) then
      Ampl_abs_r(:,:,:,:) = real( Ampl_abs(ie,:,:,:,:), db )
      Ampl_abs_i(:,:,:,:) = aimag( Ampl_abs(ie,:,:,:,:) )
    endif

    call MPI_Bcast(Ampl_abs_r,ndim,MPI_REAL8,ie_computer,MPI_COMM_GATHER,mpierr)
    call MPI_Bcast(Ampl_abs_i,ndim,MPI_REAL8,ie_computer,MPI_COMM_GATHER,mpierr)

    if( ie_computer /= mpirank ) Ampl_abs(ie,:,:,:,:) = cmplx( Ampl_abs_r(:,:,:,:), Ampl_abs_i(:,:,:,:), db )

  end do

  return
end

!************************************************************************************************************

! When extract is used, reading of the potential and other data in the bav file

subroutine Potential_reading(Delta_Eseuil,dV0bdcf,E_cut,E_Fermi,Ecineticmax,Epsii,Epsii_moy,Hubb_abs,Hubb_diag_abs, &
                             m_hubb,mpinodes,mpirank0,multi_run,ninitr,nlm_pot,nom_fich_extract,nrato_abs,nspin,nspinp, &
                             rato_abs,rhoato_abs,Rmtg_abs,Rmtsd_abs,rsato_abs,rsbdc,V_hubb_abs,V0muf,Vcato_abs,Vhbdc, &
                             Vxcato_abs,VxcbdcF)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: i, ir, isp, jsp, Lmax_abs, m, m_hubb, mp, mpierr, mpinodes, mpirank0, multi_run, n, ndim, ninitr, nlm_pot, nnombre, &
            nrato_abs, nspin, nspinp

  character(len=15):: mot
  character(len=Length_name):: nom_fich_extract

  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb_abs

  logical:: Hubb_abs, Hubb_diag_abs

  real(kind=db):: Delta_Eseuil, E_cut, E_Fermi, Ecineticmax, Rmtg_abs, Epsii_moy, Rmtsd_abs, rsbdc, V0muf, Vhbdc
  real(kind=db), dimension(nspin):: dV0bdcf, VxcbdcF
  real(kind=db), dimension(ninitr):: Epsii
  real(kind=db), dimension(nrato_abs):: rato_abs, rsato_abs
  real(kind=db), dimension(nrato_abs,nspin):: rhoato_abs
  real(kind=db), dimension(nrato_abs,nlm_pot):: Vcato_abs
  real(kind=db), dimension(nrato_abs,nlm_pot,nspin):: Vxcato_abs
  real(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb_abs_i
  real(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb_abs_r

  if( mpirank0 == 0 ) then

    Open(1,file = nom_fich_extract )

    i = 0
    do
      read(1,'(A)' ) mot
      if( mot /= ' ---- Potential' ) cycle
      i = i + 1
      if( i == multi_run ) exit
    end do

    read(1,*); read(1,*)
    read(1,*) nspin, nspinp, nrato_abs, nlm_pot, m_hubb, Lmax_abs
    read(1,*); read(1,*)
    read(1,*) E_cut, E_Fermi, V0muf, Vhbdc, VxcbdcF(:), dV0bdcF(:), rsbdc, Ecineticmax
    read(1,*); read(1,*)
    read(1,*) Hubb_abs, Hubb_diag_abs
    if( Hubb_abs ) then
      read(1,*); read(1,*)
      do isp = 1,nspinp
        do jsp = 1,nspinp
          do m = -m_hubb,m_hubb  ! it is the transpose but also in the writing
            read(1,*) ( V_hubb_abs_r(mp,m,isp,jsp), V_hubb_abs_i(mp,m,isp,jsp), mp = -m_hubb,m_hubb )
          end do
        end do
      end do
      V_hubb_abs(:,:,:,:) = cmplx( V_hubb_abs_r(:,:,:,:), V_hubb_abs_i(:,:,:,:), db )
    endif

    read(1,*); read(1,*)
    read(1,*) Rmtg_abs, Rmtsd_abs

    read(1,*); read(1,*)
    n = nnombre(1,300) - 2
    read(1,*) Epsii_moy, Delta_Eseuil, Epsii(1:n)
    if( n < ninitr ) then  ! when core_resolved is used in the extract and not in bav
      if( n == 1 ) then
        Epsii(2:ninitr) = Epsii(1)
      elseif( ninitr == 6 ) then ! L23, M23...
        Epsii(3:ninitr) = Epsii(2)
        Epsii(2) = Epsii(1)
      elseif( ninitr == 10 ) then ! M45
        Epsii(5:ninitr) = Epsii(2)
        Epsii(2:4) = Epsii(1)
      endif
    endif

    read(1,*); read(1,*)

    do ir = 1,nrato_abs
     read(1,*) rato_abs(ir), Vcato_abs(ir,:),  Vxcato_abs(ir,:,:), rhoato_abs(ir,:), rsato_abs(ir)
    end do

    Delta_Eseuil = Delta_Eseuil / Rydb
    E_cut = E_cut / Rydb
    E_Fermi = E_Fermi / Rydb
    Ecineticmax = Ecineticmax / Rydb
    Epsii(:) = Epsii(:) / Rydb
    Epsii_moy = Epsii_moy / Rydb
    V0muf = V0muf / Rydb
    Vhbdc = Vhbdc / Rydb
    VxcbdcF(:) = VxcbdcF(:) / Rydb
    dV0bdcF(:) = dV0bdcF(:) / Rydb

    if( Hubb_abs ) V_hubb_abs(:,:,:,:) = V_hubb_abs(:,:,:,:) / Rydb

    Rmtg_abs = Rmtg_abs / bohr
    Rmtsd_abs =  Rmtsd_abs / bohr
    rato_abs(:) = rato_abs(:) / bohr
    Vcato_abs(:,:) = Vcato_abs(:,:) / Rydb
    Vxcato_abs(:,:,:) = Vxcato_abs(:,:,:) / Rydb

  endif

  if( mpinodes > 1 ) then
      call MPI_Bcast(Delta_Eseuil,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(dV0bdcf,nspin,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(E_cut,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(E_Fermi,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Ecineticmax,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Epsii,ninitr,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Epsii_moy,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Hubb_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Hubb_diag_abs,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(rato_abs,nrato_abs,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(rsato_abs,nrato_abs,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(rhoato_abs,nrato_abs*nspin,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Vcato_abs,nrato_abs*nlm_pot,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Vxcato_abs,nrato_abs*nlm_pot*nspin,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Rmtg_abs,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Rmtsd_abs,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(rsbdc,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(V0muf,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(Vhbdc,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(VxcbdcF,nspin,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    ndim = ( ( 2 * m_hubb + 1 ) * nspinp )**2

    if( mpirank0 == 0 ) then
      V_hubb_abs_r(:,:,:,:) = real( V_hubb_abs(:,:,:,:), db )
      V_hubb_abs_i(:,:,:,:) = aimag( V_hubb_abs(:,:,:,:) )
    endif

      call MPI_Bcast(V_hubb_abs_i,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(V_hubb_abs_r,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    if( mpirank0 /= 0 ) V_hubb_abs(:,:,:,:) = cmplx( V_hubb_abs_r(:,:,:,:), V_hubb_abs_i(:,:,:,:), db )

  endif

  Close(1)

  return
end

