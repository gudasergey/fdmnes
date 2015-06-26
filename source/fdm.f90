! Subroutine of the FDMNES Package
! Contain the main of the FDM-MST part of the calculation

subroutine fdm(Ang_borm,Bormann,comt,Convolution_cal,Delta_edge,E_cut_imp,E_Fermi_man,Ecent,Elarg,Estart,Fit_cal, &
        Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,Green_plus,hkl_borm,icheck,ifile_notskip,indice_par,iscratch, &
        itape1,itape4,MPI_host_num_for_mumps,mpinodes0,mpirank,mpirank0,n_atom_proto_p,ngamh,ngroup_par,nnotskip,nnotskipm, &
        nomfich,nomfichbav,npar,nparm,param,Scan_a,Solver,Space_file,typepar,xsect_file)

  use declarations
  implicit none
  include 'mpif.h'

  integer, parameter:: nslapw_max = 48  ! Max number of symmetry matrix for the FLAPW data

  integer:: extract_nenerg, i, i_range, i_self, ia, iaabs, iaabsfirst, iabsfirst, iaabsi, iabsorbeur, iapr, iaprabs, ich, &
    ie, ie_computer, igr, igr_dop, igrph, igrpt_nomag, igrpt0, imparite, index_e, initl, iord, ip_max, ip0, ip00, ipr, ipr1, &
    ipr_dop, iprabs, iprabs_nonexc, ir, is, iscratch, iso1, iso2, isp1, isp2, ispin, iss1, iss2, istat, &
    it, itab, itabs, itabs_nonexc, itape1, itape4, itph, itpm, itps, itype_dop, j, je, jseuil, l, l_hubbard, l_selec_max, l1, &
    l2, lamstdens, lecrantage, lh, lin_gam, lla_state, lla2_state, lm1, lm1g, lm2, lm2g, &
    lmax, lmax_pot, lmax_t, lmaxabs, lmaxabs_t, lmaxabsa, lmaxat0, lmaxg, lmaxmax, lmaxso, lmaxso_max, lmax_probe, lmaxso0, &
    lmv1, lmv2, lseuil, m, m_hubb, m_hubb_e, m1, m2, MPI_host_num_for_mumps, mpierr, mpinodee, mpinodes, mpinodes0, mpirank, &
    mpirank0, mpirank_in_mumps_group, multi_imp, multi_run, multrmax, n, n_oo, n_atom_0, n_atom_0_t, n_atom_0_self, &
    n_atom_ind, n_atom_ind_t, n_atom_ind_self, n_atom_proto, n_atom_proto_p, n_devide, n_energy_range, &
    n_file_dafs_exp, n_multi_run, n_multi_run_e, n_orbexc, n_vr_0, n_vr_ind, n_tens_dd, n_tens_dq, n_tens_qq, &
    n_tens_t, n_tens_max, n1, natome, natome_cal, natome_self, natomeq, natomeq_coh, natomeq_s, natomeq_self, natomp, &
    natomsym, natomsym_max, nb_atom_conf_m, nb_rep, nb_sym_op, nbbseuil, nbm, nbseuil, nbtm, nbtm_fdm, nchemin, ncolm, ncolr, &
    ncolt, ndim1, ndim2, necrantage, neimagent, nenerg, nenerg_coh, nenerg_s, nenerg_tddft, neqm, ngamh, ngamme, nge, &
    ngrm, ngroup, ngroup_hubb, ngroup_lapw, ngroup_m, ngroup_neq, &
    ngroup_nonsph, ngroup_par, ngroup_pdb, ngroup_taux, ngroup_temp, ngrph, nhybm, nicm, nim, ninit1, ninitl, ninitlr, &
    ninitlt, ninitlv, nklapw, nlatm, nlm_pot, nlm_probe, nlma, nlma2, nlmagm, nlmamax, nlmam_u, nlmamax_u, &
    nlmlapwm, nlmmax, nlmomax, nlmsam, nlmsamax, nlmso, nmatsym, nnlm, nnotskip, nnotskipm, norbdil, normrmt, &
    nparm, nphiato1, nphiato7, nphim, npldafs, nple, nplr, nplrm, npoint, npoint_ns, npr, &
    npso, npsom, nptmoy, nptmoy_out, nr, nrato_dirac, nrm, nrm_self, nself, nseuil, nslapwm, &
    nso1, nsort, nsortf, nsm, nspin, nspino, nspinp, nspinorb, nstm, ntype, ntype_conf, numat_abs, nvois, nx, &
    nxanout, Trace_k, Wien_save, Z, Z_nospinorbite

  character(len=5):: Solver, Struct
  character(len=8):: PointGroup
  character(len=9):: keyword
  character(len=132):: comt, identmot, mot, nomfich, nom_fich_extract, nomfichbav, nomfich_cal_convt, &
    nomfich_s, nomfich_optic_data, nomfich_tddft_data, Space_file, xsect_file

  character(len=9), dimension(ngroup_par,nparm):: typepar
  character(len=6), dimension(10):: nomspr
  character(len=132), dimension(9):: Wien_file

  character(len=13), dimension(:), allocatable:: ltypcal
  character(len=35), dimension(:), allocatable:: com
  character(len=Length_word), dimension(:), allocatable:: nomabs
  character(len=132), dimension(:), allocatable:: nomfile_atom, nomfich_cal_conv, nomfich_cal_tddft_conv

  complex(kind=db):: f_avantseuil, f_cal
  complex(kind=db), dimension(:,:), allocatable:: karact, phdafs, phdf0t, phdt, pol, poldafsem, poldafssm
  complex(kind=db), dimension(:,:,:), allocatable:: hybrid, poldafse, poldafss
  complex(kind=db), dimension(:,:,:,:), allocatable:: secdd, secdd_m, secmd, secmd_m, secmm, secmm_m, V_hubb_abs, V_hubb_t
  complex(kind=db), dimension(:,:,:,:,:), allocatable::  rof0, secdq, secdq_m, Tau_ato, Taull, Taull_opt, V_hubb, V_hubb_s
  complex(kind=db), dimension(:,:,:,:,:,:), allocatable:: secdo, secdo_m, secoo, secoo_m, secqq, secqq_m, taull_stk
  complex(kind=db), dimension(:,:,:,:,:,:,:), allocatable:: Taull_abs

  integer, dimension(3):: hkl_borm, ldip
  integer, dimension(12):: Tensor_imp
  integer, dimension(3,3):: lqua, msymdd, msymddi
  integer, dimension(3,3,3):: loct, msymdq, msymdqi
  integer, dimension(3,3,3,3):: msymdo, msymdoi, msymqq, msymqqi
  integer, dimension(30):: icheck
  integer, dimension(nopsm):: iopsymc, iopsymr
  integer, dimension(ngroup_par) :: npar
  integer, dimension(nnotskipm) :: ifile_notskip
  integer, dimension(ngroup_par,nparm) :: indice_par
  integer, dimension(nrepm,2):: irep_util(nrepm,2)

  integer, dimension(:), allocatable:: ia_eq_inv, ia_eq_inv_self, iabsm,  iabsorig, iapot, iaproto, iaprotoi, icom, igr_i, &
     igr_is, igr_proto, igroup, igroupi, imoy, imoy_out, is_g, ispin_maj, isrt, isymeq, itdil, &
     its_lapw, itype, itypei, itypep, itypepr, Kgroup, ldil, lmaxa, lmaxat, &
     lqnexc, lval, nb_eq, nbord, nbordf, ngreq, ngreqm, nlat, nlmsa, nlmso0, norbv, nqnexc, nphi_dafs, nposextract, nrato, &
     nrato_lapw, nsymextract, numat, numia

  integer, dimension(:,:), allocatable:: hkl_dafs, ia_eq, ibord, igreq, indice, iopsym_atom, is_eq, isbord, &
     isigpi, isymqa, isvois, ivois, iso, lcoeur, lso, lvval, m_g, mso, nb_rpr, ncoeur, nlmsa0, nvval

  integer, dimension(:,:,:), allocatable:: ia_rep, iato, lato, mato, mpres, nb_rep_t, Wien_matsym

  integer, dimension(:,:,:,:), allocatable:: msymoo, msymooi

  logical:: Absauto, Absorbeur, Allsite, ATA, Atom_nonsph, Atom_nonsph_loc, Atom_occ_hubb, Atomic_scr, Axe_loc, Basereel, &
     Base_hexa, Base_ortho, Base_spin, Bormann, BSE, Clementi, Cal_xanes, Cartesian_tensor, Charge_free, Convergence, &
     Core_resolved, Convolution_cal, Coupelapw, Dafs, Dafs_bio, Density, Density_comp, Devide_Ei, Dipmag, &
     Doping, Dyn_eg, Dyn_g, E_comp, E1E2e, E_Fermi_man, Eneg, Eneg_i, Eneg_n_i, Energphot, Extract, &
     Extract_Green, Fermi, Fermi_first, Final_optic, Final_tddft, Fit_cal, Flapw, Flapw_new, Force_ecr, &
     Full_atom, Full_atom_e, Full_potential, Full_self_abs, Gamma_hole_imp, Gamma_tddft, Green, Green_i, Green_int, &
     Green_plus, Green_s, Green_self, Hubb_a, Hubb_d, Hubbard, key_calc, korigimp, &
     Level_val_abs, Level_val_exc, lmaxfree, lmoins1, lplus1, Magnetic, Matper, Memory_save, Moy_cluster, &
     Moy_loc, Moyenne, Muffintin, No_solsing, Noncentre, Nonexc_g, Nonexc, Normaltau, Octupole, Old_reference, One_run, &
     Optic, Optic_xanes, Overad, Pdb, PointGroup_Auto, Polarise, Proto_all, Quadmag, &
     Quadrupole, Readfast, Recop, Recup_optic, Recup_tddft_data, Relativiste, RPALF, &
     Rydberg, Save_optic, Save_tddft_data, Scan_a, Scf_elecabs, SCF_mag_fix, &
     SCF_mag_free, Second_run, Self_abs, Self_cons, Self_nonexc, &
     Solsing, Solsing_s, Solsing_o, Solsing_only, Spherical_signal, Spherical_tensor, Spino, Spinorbite, State_all, &
     State_all_out, Supermuf, Sym_4, Sym_cubic, Symauto, Symmol, Tau_nondiag, Taux, Tddft, &
     Tddft_so, Tddft_xanes, Temperature, Trace_format_wien, Xan_atom, Ylm_comp, Ylm_comp_inp

  logical, dimension(4):: SCF_log
  logical, dimension(10):: Multipole

  logical, dimension(:), allocatable:: Atom_axe, Atom_with_axe, Atom_nsph_e, Hubb, Hubb_diag, Repres_comp, &
     Run_done, Skip_run

  real(kind=db):: adimp, alfpot, Ang_borm, Cal_Volume_maille, chargm, d_ecrant, D_max_pot, Delta_E_conv, Delta_edge, &
     Delta_En_conv, Delta_energ, Delta_energ_s, Delta_energ_t, Delta_Epsii, Delta_Eseuil, Densite_atom, E_cut, &
     E_cut_imp, E_Fermi, E_Open_val, E_Open_val_exc, E_start, Ecent, Ecineticmax, Ecineticmax_out, &
     Eclie, Ecmax, Ecmax_out, Ei0, Eii, Elarg, Em, En_cluster, En_cluster_s, En_cluster_t, &
     Energ_max, Enervide, Enragr, Epsii_moy, Estart, Extract_E_cut, Extract_V0bdcF, &
     Gamma_max, Kern_fac, overlap, p_self, p_self_s, p_self_t, p_self0, Pas_SCF, R_rydb, R_self, Rmax, &
     Rmtsd_abs, Roverad, Rpotmax, rsbdc, rsbdc_out, Rsort, Rsorte, Rsorte_s, Rtph, &
     Temp, Test_dist_min, tp_init, tp_SCF_1, tp_SCF_2, tp1, tp2, tpt_tot, tpt1, tpt2, tptt, V_intmax, V0muf, &
     Vhbdc, Vhbdc_out, Volume_maille, Vsphere, Workf

  real(kind=db), dimension(2):: chg_open_val, f_no_res, pop_open_val
  real(kind=db), dimension(3):: Ang_rotsup, angxyz, axyz, dcosxyz, deccent, dpos, Vec_orig
  real(kind=db), dimension(6):: Trace_p
  real(kind=db), dimension(7):: tp
  real(kind=db), dimension(10) :: Gamma_hole
  real(kind=db), dimension(12) :: tpt
  real(kind=db), dimension(3,3):: Cubmat, Orthmat, Mat_or, Orthmati, Rot_atom_abs, Rot_int
  real(kind=db), dimension(ngroup_par,nparm):: param

  real(kind=db), dimension(:), allocatable:: Angle_or, cdil, ch_coeur, chargat, chg_cluster, cgrad, clapl, dista, distai, &
     dV0bdcF, dvc_ex_nex, dv_ex_nex, E_starta, Ecinetic, Ecinetic_out, Ecrantage, Eeient, Energ, Egamme, Eimag, &
     Eimag_coh, Eimag_s, Eimagent, Energ_coh, Energ_s, Energ_self, Energ_self_s, En_coeur, En_coeur_s, Epsii, Eseuil, &
     poidso, poidsov, poidsov_out,  popatc, r0_lapw, r, rchimp, &
     rhons, rlapw, rmt, rmtimp, rmtg, rmtg0, rmtsd, rs, rsato_abs, rvol, sec_atom, Taux_eq, Taux_ipr, Taux_oc, Temp_coef, &
     V_hubbard, Vh, Vhns, V0bdc, V0bdcF, V0bdc_out, V0bdcFimp, VxcbdcF, VxcbdcF_out
  real(kind=db), dimension(:,:), allocatable :: coef_g, angpoldafs, Axe_atom_clu, Axe_atom_gr, Axe_atom_grn, &
     chargat_init, chargat_self, chargat_self_s, drho_ex_nex, Excato, Int_tens, pdp, pdpolar, poidsa, &
     polar, pop_nonsph, pop_orb_val, popatv, popexc, pos, poseq, posi, posi_self, posn, psi_open_val, psii, &
     rato, rho, rhoa, rho_coeur, rhoato_abs, rhoit, rsato, rsato_t, V_abs_i, Vcato_init, vecdafsem, &
     vecdafssm, veconde, Vr, Vxc, vec, Wien_taulap, xyz, Ylmso
  real(kind=db), dimension(:,:,:), allocatable:: gradvr, Int_statedens, popatm, &
     popats, popval, posq, psi_coeur, psival, rho_chg, rho_self_t, rhoato_init, rot_atom, rot_atom_gr, &
     rotloc_lapw, Vcato, Vcato_t, vecdafse, vecdafss, Vra, Vrato_e, Ylmato
  real(kind=db), dimension(:,:,:,:), allocatable:: occ_hubb_e, rho_self, rho_self_s, Vrato, Vxcato, Vxcato_t
  real(kind=db), dimension(:,:,:,:,:), allocatable:: drho_self, imag_taull, occ_hubb,  occ_hubb_i
  real(kind=db), dimension(:,:,:,:,:,:), allocatable:: phiato, Statedens, Statedens_i

  real(kind=sg) time

  parameter( n_tens_dd=9, n_tens_dq=15, n_tens_qq=25, n_tens_t = n_tens_dd + n_tens_dq + n_tens_qq, &
           n_tens_max = 8 + 2 * n_tens_t + 2 * n_tens_dq )

  data nomspr/'Lectur','Reseau','Potent','Ylm   ','Potex ', 'Sphere','Mat   ','Tensor','Densit','Coabs '/

  mpirank_in_mumps_group = mod( mpirank0, MPI_host_num_for_mumps )
  mpinodes = mpinodes0 / MPI_host_num_for_mumps

  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    tp1 = real(time,db)
  endif

  E_Fermi = -5._db / Rydb
  E_cut = E_Fermi
  E_Open_val = E_Fermi
  E_Open_val_exc = E_Fermi
  Moy_cluster = .true.
  iaabsfirst = 0
  ipr_dop = 0

  call lectdim(Absauto,Atom_occ_hubb,Atom_nonsph,Axe_loc,Bormann,Doping,Extract,Flapw,Full_self_abs,Hubbard,itape4, &
    Magnetic,Memory_save,mpinodes0,mpirank0,n_file_dafs_exp,n_multi_run_e,nb_atom_conf_m,ncolm,neimagent, &
    nenerg_s,ngamme,ngroup,ngroup_neq,nhybm,nklapw,nlatm,nlmlapwm,nmatsym,norbdil,npldafs,nple,nplrm,nspin,nspino,nspinp,ntype, &
    ntype_conf,Pdb,Readfast,Self_abs,Space_file,Taux,Temperature,Xan_atom)

  if( Atom_nonsph ) then
    ngroup_nonsph = ngroup
  else
    ngroup_nonsph = 1
  endif
  if( Flapw ) then
    ngroup_lapw = ngroup_neq
  else
    ngroup_lapw = 1
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
  if( Temperature ) then
    ngroup_temp = ngroup
  else
    ngroup_temp = 0
  endif
  if( Flapw ) then
    nslapwm = nslapw_max
  else
    nslapwm = 0
  endif

  allocate( Angle_or(npldafs) )
  allocate( egamme(ngamme) )
  allocate( eeient(neimagent) )
  allocate( eimagent(neimagent) )
  allocate( angpoldafs(3,npldafs) )
  allocate( Atom_nsph_e(ngroup) )
  allocate( Axe_atom_gr(3,ngroup_m) )
  allocate( Axe_atom_grn(3,ngroup_m) )
  allocate( com(0:ntype) )
  allocate( cdil(norbdil) )
  allocate( ecrantage(nspin) )
  allocate( hkl_dafs(3,npldafs) )
  allocate( hybrid(nhybm,16,ngroup_nonsph) )
  allocate( iabsm(n_multi_run_e) )
  allocate( iabsorig(n_multi_run_e) )
  allocate( icom(0:ntype) )
  allocate( isigpi(npldafs,2) )
  allocate( itdil(norbdil) )
  allocate( its_lapw(ngroup_lapw) )
  allocate( itype(ngroup) )
  allocate( Kgroup(ngroup_pdb) )
  allocate( ldil(norbdil) )
  allocate( lvval(0:ntype,nlatm) )
  allocate( nlat(0:ntype) )
  allocate( nomfile_atom(0:ntype) )
  allocate( norbv(0:ngroup_nonsph) )
  allocate( nphi_dafs(npldafs) )
  allocate( nrato(0:ntype) )
  allocate( nrato_lapw(0:ntype) )
  allocate( nposextract(n_multi_run_e) )
  allocate( nsymextract(n_multi_run_e) )
  allocate( numat(0:ntype) )
  allocate( nvval(0:ntype,nlatm) )
  allocate( occ_hubb_e(-m_hubb_e:m_hubb_e,-m_hubb_e:m_hubb_e,nspin,ngroup_hubb) )
  allocate( pdpolar(nple,2) )
  allocate( polar(3,nple) )
  allocate( poldafsem(3,npldafs) )
  allocate( poldafssm(3,npldafs) )
  allocate( pop_nonsph(nhybm,ngroup_nonsph) )
  allocate( popatc(0:ntype) )
  allocate( popats(ngroup,nlatm,nspin) )
  allocate( popatv(0:ntype,nlatm) )
  allocate( popval(0:ntype,nlatm,nspin) )
  allocate( posn(3,ngroup) )
  allocate( r0_lapw(0:ntype) )
  allocate( rchimp(0:ntype) )
  allocate( rlapw(0:ntype) )
  allocate( rmt(0:ntype) )
  allocate( rmtimp(0:ntype) )
  allocate( rotloc_lapw(3,3,ngroup_lapw) )
  allocate( Rot_Atom_gr(3,3,ngroup_m) )
  allocate( Taux_oc(ngroup_taux) )
  allocate( Temp_coef(ngroup_temp) )
  allocate( V_hubbard(0:ntype) )
  allocate( Hubb(0:ntype) )
  allocate( V0bdcFimp(nspin) )
  allocate( Vecdafsem(3,npldafs) )
  allocate( Vecdafssm(3,npldafs) )
  allocate( Veconde(3,nple) )
  allocate( Wien_matsym(3,3,nslapwm) )
  allocate( Wien_taulap(3,nslapwm) )

  if( mpinodes0 > 1 ) call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  call lecture(Absauto,adimp,alfpot,Allsite,Ang_borm,Ang_rotsup,Angle_or,Angpoldafs,Angxyz,ATA,Atom_occ_hubb,Atom_nonsph, &
    Atom_nsph_e,Atomic_scr,Axe_atom_gr,Axe_loc,axyz,Base_spin,Basereel,Bormann,BSE,Cartesian_tensor,Charge_free,Clementi, &
    com,comt,Core_resolved,Coupelapw,Cubmat,D_max_pot,Dafs,Dafs_bio,Delta_En_conv, &
    Delta_Epsii,Density,Density_comp,Dipmag,Doping,dpos,Dyn_eg,Dyn_g,Eclie,Ecrantage,Eeient,Egamme,Eimagent, &
    Eneg_i,Eneg_n_i,Energphot,Extract,f_no_res,Fit_cal,Flapw,Flapw_new, &
    Force_ecr,Full_atom_e,Full_potential,Full_self_abs,Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,Green_int,Green_s, &
    Green_self,hkl_borm,hkl_dafs,Hubb,Hubbard,hybrid,iabsm,iabsorig,icheck,icom,igr_dop,indice_par,iscratch, &
    isigpi,itdil,its_lapw,iord,itape4,itype,itype_dop,jseuil,Kern_fac,Kgroup,korigimp,l_selec_max, &
    lamstdens,ldil,lecrantage,lin_gam,lmax_pot,lmaxfree,lmaxso0,lmaxat0,lmoins1,lplus1,lseuil,lvval,m_hubb_e, &
    Magnetic,Mat_or,Matper,MPI_host_num_for_mumps,mpinodes,mpinodes0,mpirank0,Muffintin, &
    Multipole,multrmax,n_atom_proto,n_devide,n_file_dafs_exp,n_multi_run_e, &
    nb_atom_conf_m,nbseuil,nchemin,necrantage,neimagent,nenerg_s,ngamh,ngamme,ngroup,ngroup_hubb,ngroup_lapw,ngroup_m, &
    ngroup_neq,ngroup_nonsph,ngroup_par,ngroup_pdb,ngroup_taux, &
    ngroup_temp,nhybm,nlat,nlatm,nnlm,No_solsing,nom_fich_extract, &
    nomfich,nomfich_optic_data,nomfich_tddft_data,nomfichbav,nomfile_atom,Noncentre, &
    Nonexc,norbdil,norbv,Normaltau,normrmt,npar,nparm,nphi_dafs, &
    nphim,npldafs,nple,nposextract,nrato,nrato_dirac,nrato_lapw,nrm, &
    nself,nseuil,nslapwm,nspin,nsymextract,ntype,ntype_conf,numat,numat_abs, &
    nvval,occ_hubb_e,Octupole,Old_reference,One_run,Optic,Overad,Overlap,p_self0, &
    param,Pas_SCF,pdpolar,PointGroup,PointGroup_Auto,Polar,Polarise,poldafsem,poldafssm, &
    pop_nonsph,popatc,popats,popatv,popval,posn,Quadmag,Quadrupole,R_rydb, &
    r0_lapw,rchimp,Readfast,Recup_optic,Recup_tddft_data,Relativiste,r_self,rlapw,rmt,rmtimp,Rot_Atom_gr,rotloc_lapw, &
    roverad,RPALF,rpotmax,Rydberg,rsorte_s,Save_optic,Save_tddft_data,SCF_log,Self_abs, &
    Solsing_s,Solsing_only,Solver,Space_file,Spherical_signal,Spherical_tensor, &
    Spinorbite,State_all,State_all_out,Struct,Supermuf,Symauto,Symmol,Taux,Taux_oc,Tddft,Tddft_so, &
    Temp,Temp_coef,Temperature,Tensor_imp,Test_dist_min,Trace_format_wien,Trace_k,Trace_p,Typepar,V_hubbard,V_intmax,Vec_orig, &
    Vecdafsem,Vecdafssm,Veconde,V0bdcFimp,Wien_file,Wien_matsym,Wien_save,Wien_taulap,Ylm_comp_inp,Z_nospinorbite)

  E1E2e = Multipole(2)

  SCF_elecabs = SCF_log(1)
  SCF_mag_free = SCF_log(2)
  Self_cons = SCF_log(3)
  Self_nonexc = SCF_log(4)

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

  if( mpirank0 /= 0 .and. Extract ) then

    deallocate( Angle_or, angpoldafs, Axe_atom_gr, Axe_atom_grn )
    deallocate( com, cdil, Ecrantage, hkl_dafs, hybrid, iabsm )
    deallocate( egamme,  eeient, eimagent )
    deallocate( iabsorig, icom, isigpi, itdil, itype, its_lapw )
    deallocate( Kgroup, ldil, lvval, nlat, nomfile_atom )
    deallocate( norbv, nphi_dafs, nposextract, nrato, nrato_lapw )
    deallocate( nsymextract, numat, nvval, occ_hubb_e, pdpolar )
    deallocate( polar, poldafsem, poldafssm, pop_nonsph )
    deallocate( popatc, popats, popatv, popval, posn, r0_lapw )
    deallocate( rchimp, rlapw, rmt, rmtimp, Rot_Atom_gr )
    deallocate( rotloc_lapw, Taux_oc, Temp_coef, Hubb, V_hubbard )
    deallocate( V0bdcFimp, vecdafsem, vecdafssm, veconde )
    deallocate( Wien_matsym, Wien_taulap )

    return

  endif

! lmax pour la densite d'etat
  lla_state = max( 4, lmax_pot )
  lla2_state = ( lla_state + 1 )**2
  nlm_pot = ( lmax_pot + 1 )**2

  if( Extract ) then
    mpinodee = 1
  else
    mpinodee = mpinodes0
  endif

  if( nphim > 1 ) Scan_a = .true.

  allocate( Eseuil(nbseuil) )
  allocate( lcoeur(2,0:ntype) )
  allocate( ncoeur(2,0:ntype) )
  allocate( popexc(nnlm,nspin) )
  allocate( psi_coeur(0:nrm,2,0:ntype) )
  allocate( psii(nrm,nbseuil) )
  allocate( psival(0:nrm,nlatm,0:ntype) )
  allocate( psi_open_val(nrm,2) )
  allocate( rato(0:nrm,0:ntype) )
  allocate( rhoit(0:nrm,0:ntype) )
  allocate( rho_coeur(0:nrm,0:ntype) )

  call atom(Clementi,com,icheck(2),icom,itype,jseuil,lcoeur,lseuil,lvval,mpinodee,mpirank0,nbseuil,ncoeur,ngroup,nlat, &
        nlatm,nnlm,nomfile_atom,Nonexc,nrato,nrato_dirac,nrato_lapw,nrm,nseuil,nspin,ntype,numat,nvval,popatc,popats, &
        popatv,popexc,popval,psi_coeur,psii,psival,r0_lapw,rato,Relativiste,rho_coeur,rhoit,rlapw,rmt)

  allocate( Atom_with_axe(0:ngroup_m) )
  allocate( igr_i(ngroup) )
  allocate( igr_is(ngroup) )
  allocate( igr_proto(ngroup) )

  call symsite(absauto,angxyz,Atom_with_axe,Atom_nonsph,Atom_nsph_e,Axe_atom_gr,axyz,Base_ortho,Cubmat,dcosxyz, &
        Doping,Extract,Flapw,iabsm,icheck(3),igr_i,igr_is,igr_proto,itype, &
        Magnetic,Matper,Memory_save,n_atom_proto,n_multi_run_e,neqm, &
        ngroup,ngroup_m,nlat,nlatm,nspin,ntype,numat,numat_abs,popats,posn,Struct)

  deallocate( Atom_nsph_e )

  if( Doping ) n_atom_proto = n_atom_proto + 1

  allocate( igreq(0:n_atom_proto,neqm) )
  allocate( ngreq(0:n_atom_proto) )

  ngreq(:) = 0
  igreq(:,:) = 0

  do igr = 1,ngroup
    if( Doping .and. igr == ngroup ) exit
    ipr = igr_proto(igr)
    ngreq(ipr) = ngreq(ipr) + 1
    i = igr_i(igr)
    igreq(ipr,i) = igr
  end do

  if( Doping ) then
    igreq(n_atom_proto,1) = ngroup
    ngreq(n_atom_proto) = 1
  endif

  if( absauto ) then
    multi_run = 0
    do ipr = 1,n_atom_proto
      it = abs( itype(igreq(ipr,1)) )
      if( numat( it ) /= numat_abs ) cycle
      multi_run = multi_run + 1
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

  if( absauto ) then
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
    do ipr = 1,n_atom_proto
      it = abs( itype(igreq(ipr,1)) )
      if( numat( it ) /= numat_abs ) cycle
      multi_run = multi_run + 1
      iabsm(multi_run) = igreq(ipr,1)
      iabsorig(multi_run) = ipr
    end do
  endif

  natomsym_max = 0
  do igr = 1,ngroup
    if( Doping ) then
      if( numat( abs( itype(igr) ) ) /= numat( abs( itype(igr_dop) ) ) ) cycle
    else
      if( numat( abs( itype(igr) ) ) /= numat_abs ) cycle
    endif
    natomsym_max = max( igr_i(igr), natomsym_max )
  end do

  allocate( isymqa(0:n_atom_proto,natomsym_max) )
  allocate( posq(3,n_atom_proto,natomsym_max) )
  isymqa(:,:) = 0
  posq(:,:,:) = 0._db

  do igr = 1,ngroup

    if( Doping .and. igr == ngroup ) then
      posq(:,n_atom_proto,:) = posq(:,ipr_dop,:)
      isymqa(n_atom_proto,:) = isymqa(ipr_dop,:)
      exit
    endif
    if( Doping ) then
      if( numat( abs( itype(igr) ) ) /= numat( abs( itype(igr_dop) ) ) ) cycle
    else
      if( numat( abs( itype(igr) ) ) /= numat_abs ) cycle
    endif

    if( igr == igr_dop ) ipr_dop = igr_proto(igr)

    ipr = igr_proto(igr)
    i = igr_i(igr)
    posq(:,ipr,i) = posn(:,igr)
    isymqa(ipr,i) = igr_is(igr)

  end do

  deallocate( igr_i, igr_is, igr_proto )

  Volume_maille = Cal_Volume_maille(axyz,angxyz)
  Densite_atom = ngroup / Volume_maille

  call esdata(Eseuil,icheck(4),jseuil,nbseuil,nseuil,numat_abs,Old_reference,Workf,mpirank0)

  if( Flapw ) Workf = 0._db

  if( Optic ) then
    ninit1 = 0; ninitl = 0
  else
    call dim_init(jseuil,lseuil,nbseuil,ninit1,ninitl)
  endif
  allocate( coef_g(ninitl,2) )
  allocate( m_g(ninitl,2) )
  allocate( is_g(ninitl) )
  if( .not. Optic ) call coef_init(coef_g,is_g,ninitl,jseuil,lseuil,m_g,nbseuil)

  allocate( chargat(0:n_atom_proto) )
  allocate( ngreqm(n_multi_run) )
  allocate( popatm(0:n_atom_proto,nlatm,nspin) )
  allocate( Run_done(n_multi_run) )
  Run_done(:) = .false.
  allocate( Skip_run(n_multi_run) )
  Skip_run(:) = .false.

  call pop_proto(chargat,Charge_free,chargm,Flapw,iabsm,icheck(4),igreq,itype,mpirank0,n_atom_proto,n_multi_run,neqm, &
        ngreq,ngreqm,ngroup,nlat,nlatm,nspin,ntype,numat,popatc,popatm,popats,Run_done)

  if( nnotskip > 0 .and. Extract .and. n_atom_proto_p == n_atom_proto ) then
    do multi_run = 1,n_multi_run
      Skip_run(multi_run) = .true.
      do i = 1,nnotskip
        if( ifile_notskip(i) == iabsorig(multi_run) ) Skip_run(multi_run) = .false.
      end do
    end do
  endif

  n_atom_proto_p = n_atom_proto

  if( Extract ) then
    Green_s = Extract_Green(nom_fich_extract)
  else
    if( Recup_tddft_data ) call Recup_nenerg(Energ_max,mpirank0,nenerg_s,nomfich_tddft_data)
    if( Recup_optic ) call Recup_optic_dim(E_cut,Full_atom,iaabsi,iprabs,mpirank0,natome,nenerg_s,nlm_pot,nlm_probe, &
                                           nomfich_optic_data)
    allocate( Energ_s(nenerg_s) )
    allocate( Eimag_s(nenerg_s) )
    if( Recup_tddft_data .or. Recup_optic ) then
      if( Recup_tddft_data ) Energ_s(nenerg_s) = Energ_max
    else
      call grille_xanes(eeient,eimag_s,eimagent,egamme,Energ_s,icheck(4),lin_gam,ngamme,neimagent,nenerg_s)
    endif
  endif

  deallocate( egamme, eeient, eimagent )
  if( Multipole(7) ) then
    n_oo = 9
  else
    n_oo = 0
  endif
  allocate( msymoo(3,n_oo,3,n_oo) )
  allocate( msymooi(3,n_oo,3,n_oo) )

  if( mpirank0 == 0 ) then
    call CPU_TIME(time)
    tp2 = real(time,db)
    tpt(1) = tp2 - tp1
  endif

  allocate( nomfich_cal_conv(n_multi_run) )
  if( Tddft ) allocate( nomfich_cal_tddft_conv(n_multi_run) )

! Boucle sur tous les absorbeurs nonequivalents selon le groupe d'espace
  if( mpirank0 == 0 ) write(6,130) n_multi_run

  boucle_multi: do multi_run = 1,n_multi_run

    if( mpirank0 == 0 ) then
      tpt(2:12) = 0._db
      tpt1 = 0._db
      tpt2 = 0._db
      call CPU_TIME(time)
      tp(1) = real(time,db)
      tp_init = tp(1)
    endif

    Second_run = multi_run > 1 .and. One_run

    multi_imp = nposextract(multi_run)

! ninitlr devient 1 pour le final Tddft
    if( Core_resolved ) then
      if( Optic ) then
     ! En optique il y a une sortie par "l" plus le total
        if( numat_abs <= 4 ) then
          ninitlr = 1
        elseif( numat_abs > 4 .and. numat_abs <= 20 ) then
          ninitlr = 2
        elseif( numat_abs > 21 .and. numat_abs <= 56 ) then
          ninitlr = 3
        else
          ninitlr = 4
        endif
        ninitlr = ninitlr + 1
      else
        ninitlr = ninitl
      endif
    else
      ninitlr = nbseuil
    endif

    if( multi_run > 1 )  deallocate( Epsii )
    allocate( Epsii(ninitlr) )

    if( Run_done(multi_run) .or. Skip_run(multi_run) ) cycle

    if( Extract ) then
      if( Tddft ) then
        nbbseuil = 1
      else
        nbbseuil = nbseuil
      endif
    endif

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

    nomfich_s = nomfich
    Final_tddft = .false.
    Final_optic = .false.
    m = min( multi_run, n_multi_run_e )

    allocate( lqnexc(nnlm) ); allocate( nqnexc(nnlm) )

    call init_run(com,Doping,Ecrantage,Force_ecr,Hubb,iabsm(multi_run),iabsorig(multi_run),icheck(4),icom, &
      itabs,itype,itype_dop,jseuil,lcoeur,lecrantage,lqnexc,lseuil,lvval,mpinodee,mpirank0,n_multi_run,n_orbexc,nbseuil,ncoeur, &
      necrantage,ngroup,nlat,nlatm,nnlm, nomfich_s,Nonexc,nqnexc, &
      nrato,nrato_dirac,nrm,nseuil,nspin,ntype,numat,nvval,pop_open_val,popatc,popats,popatv,popexc,popval, &
      psi_coeur,psii,psi_open_val,psival,rato,rchimp,Relativiste,rho_coeur,rhoit,rmt,rmtimp,V_hubbard,Wien_file(8))

    boucle_1: do iprabs_nonexc = 1,n_atom_proto
      do i = 1,ngreq(iprabs_nonexc)
        if( abs(igreq(iprabs_nonexc,i)) == iabsorbeur) exit boucle_1
      end do
    end do boucle_1

    if( Nonexc ) then
      iprabs = iprabs_nonexc
    else
      iprabs = 0
    endif

    if( Doping ) then
      natomsym = ngreq(ipr_dop)
    else
      natomsym = ngreq(iprabs_nonexc)
    endif
    allocate( poseq(3,natomsym) )
    allocate( isymeq(natomsym) )
    allocate( Taux_eq(natomsym) )

    poseq(:,1:natomsym) = posq(:,iprabs_nonexc,1:natomsym)
    isymeq(1:natomsym) = isymqa(iprabs_nonexc,1:natomsym)
    if( Taux ) then
      Taux_eq(1:natomsym) = Taux_oc( abs( igreq(iprabs_nonexc,1:natomsym) ) )
    else
      Taux_eq(1:natomsym) = 1._db
    endif

    if( Extract ) then
      E_cut = extract_E_cut(multi_imp,nom_fich_extract)
      if( Optic ) then
        V0muf = 0._db
        Epsii(:) = 0._db
        Epsii_moy = 0._db
        Close(1)
      else
        allocate( V0bdcF(1) )
        V0bdcF(1) = extract_V0bdcF()
        V0muf = WorkF + V0bdcF(1)
        deallocate( V0bdcf )
        call extract_Epsii(Core_resolved,Delta_Epsii,Delta_Eseuil,Epsii,Epsii_moy,icheck(14),nbseuil, &
                 ninit1,ninitl,ninitlr)
      endif
      nenerg_s = extract_nenerg(multi_imp,nom_fich_extract,Tddft)
      allocate( Energ_s(nenerg_s) )
      allocate( eimag_s(nenerg_s) )
      call extract_energ(Energ_s,Eseuil,multi_imp,nbbseuil,nbseuil,nenerg_s,nom_fich_extract,Tddft)
    endif

    allocate( iapot( 0:n_atom_proto ) )
    allocate( itypepr( 0:n_atom_proto ) )

    call Pop_mod(chargat,chargm,Doping,Flapw,iabsorbeur,icheck(4),igreq,itabs,iprabs_nonexc,itype,itype_dop,itypepr,lqnexc,lvval, &
        n_atom_proto,n_orbexc,neqm,ngreq,ngroup,nlat,nlatm,nnlm,nqnexc,nspin,ntype,numat_abs,nvval,popatm,popexc)

    deallocate( lqnexc, nqnexc )

    d_ecrant = 1 - sum( ecrantage(:) )
    if( One_run ) then
      iabsfirst = iabsm(1)
    else
      iabsfirst = iabsorbeur
    endif

    if( ( .not. One_run ) .or. multi_run == 1 ) &
      call natomp_cal(angxyz,ATA,axyz,Base_ortho,chargat,d_ecrant,dcosxyz,deccent,Doping,dpos,Flapw,iabsorbeur,iabsfirst, &
        icheck(5),igr_dop,igreq,itabs,itype,Kgroup,Matper,mpirank0,multrmax,n_atom_proto,natomeq_s,natomeq_coh, &
        natomp,neqm,ngreq,ngroup,ngroup_pdb,ngroup_taux,Noncentre,posn,Proto_all,r_self,rsorte_s,rmax,rpotmax, &
        Self_cons,Taux,Taux_oc,Test_dist_min)

    allocate( Axe_atom_clu(3,natomp) )
    allocate( iaproto(natomp) )
    allocate( igroup(natomp) )
    allocate( itypep(natomp) )
    allocate( dista(natomp) )
    allocate( pos(3,natomp) )

    if( .not. Extract ) then
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
      allocate( VxcbdcF(nspin) )
      allocate( VxcbdcF_out(nspin) )
    end if

    if( Optic .and. .not. Extract ) then
! En optique, on traite differemment les états non excites sous Fermi
      n_energy_range = 1
    else
      n_energy_range = 1
    endif

! Boucle etats occupes - etats inoccupes
    boucle_energy_range: do i_range = 1,n_energy_range

      if( n_energy_range == 2 ) then
        do ipr = 3,6,3
          if( icheck(1) == 0 .and. ipr == 3 ) cycle
          if( i_range == 1 ) then
            write(ipr,'(/a21)') '  Occupied state part'
          else
            write(ipr,'(/a20)') '  Excited state part'
          endif
        end do
      endif

! Preambule coherence:
      if( Self_cons .and. .not. Second_run .and. .not. Recup_optic ) then
        Convergence = .false.
      else
        Convergence = .true.
      end if
      if( Self_cons .and. .not. Second_run .and. (Optic .or. tddft) ) then
        Fermi_first = .true.
        Devide_Ei = .true.
      else
        Fermi_first = .false.
        Devide_Ei = .false.
      endif

! Boucle coherence
      boucle_coh: do i_self = 1,nself+1

        if( Convergence .and. i_self < nself + 1 ) cycle

        if( i_self == nself + 1 ) then
          Cal_xanes = .true.
        else
          Cal_xanes = Convergence
          Fermi = .false.
          E_Fermi = -5._db / Rydb
          E_cut = -5._db / Rydb
        end if

        if( i_self == 1 ) then
          call CPU_TIME(time)
          tp_SCF_1 = real(time,db)
        endif
        if( Cal_xanes ) then
          if( Recup_optic ) then
            tpt(11) = 0._db
          else
            call CPU_TIME(time)
            tp_SCF_2 = real(time,db)
            tpt(11) = tp_SCF_2 - tp_SCF_1
          endif
        endif

        Level_val_abs = .false.
        Level_val_exc = .false.

! Par defaut, le calcul auto-coherent se fait en diffusion multiple
! Un calcul non auto-coherent est mene en mode donne en entree
        if( Cal_xanes ) then
          natomeq = natomeq_s
          Green = Green_s
          rsorte = rsorte_s
         else
          natomeq = natomeq_coh
          Green = Green_self
          rsorte = r_self
        endif
        if( Cal_xanes .and. ( Tddft .or. Optic ) ) then
          Eneg = .false.
        elseif( Green ) then
          Eneg = .not. Eneg_n_i
        else
          Eneg = Eneg_i
        endif

        if( icheck(4) > 0 ) then
          if( Cal_xanes ) then
            write(3,140)
          else
            write(3,150) i_self
          endif
        endif

        if( Cal_xanes .and. icheck(27) > 0 ) write(3,160) E_cut * Rydb

        if( i_self > 1 .and. Cal_xanes .and. .not. ( Second_run .or. Recup_optic ) ) deallocate( karact )
        if( i_self == 1 .or. Cal_xanes ) then
          allocate( karact(nopsm,nrepm) )
          if( Cal_xanes ) then
            Nonexc_g = Nonexc
          else
            Nonexc_g = Self_nonexc
          endif
          if( Nonexc_g ) then
            ipr1 = 1
            itab = itabs_nonexc
          else
            ipr1 = 0
            itab = itabs
          endif
          if( Cal_xanes .and. i_range == 1 ) then
            ich = icheck(5)
          elseif( icheck(27) > 1 ) then
            ich = max( icheck(5), icheck(27) )
          else
            ich = 0
          endif

          call agregat(angxyz,ATA,Atom_with_axe,Atom_nonsph,Axe_atom_clu,Axe_atom_gr,Axe_atom_grn,axyz,Base_hexa,Base_ortho, &
            chargat,chargm,Cubmat,dcosxyz,deccent,dista,Doping,dpos,Flapw,iaabs,iaabsfirst,iabsorbeur,iaproto,iapot, &
            ich,igr_dop,igreq,igroup,igrpt_nomag,igrpt0,iopsymc,iopsymr,itab,itype,itypep,karact,Kgroup,Magnetic,Matper, &
            mpirank0,multi_run,n_atom_proto,natomp,nb_rep,nb_sym_op,neqm,ngreq,ngroup,ngroup_m,ngroup_pdb,ngroup_taux,nlat, &
            nlatm,Noncentre,Nonexc_g,nspin,ntype,numat,One_run,Orthmat,Orthmati,PointGroup,PointGroup_Auto, &
            popats,pos,posn,rmax,Rot_int,Self_nonexc,Spinorbite,Rot_Atom_gr,Sym_4,Struct,Sym_cubic,Symmol,Taux,Taux_oc,Vsphere)

          if( .not. One_run .or. ( One_run .and. multi_run == 1 ) ) iaabsfirst = iaabs

        endif

        if( Cal_xanes ) then

! Evaluation de la forme des tenseurs cartesiens
          call Tensor_shape(Atom_with_axe,Atom_nonsph,Axe_atom_clu,Base_ortho,dcosxyz,Dipmag, &
               E1E2e,Green,iaabs,icheck(6),igroup,igrpt0,iopsymc,iopsymr,itype,itypep,ldip,loct,lqua,lseuil,Magnetic, &
               msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole,n_oo,natomp,ngroup, &
               ngroup_m,nlat,nlatm,nspin,ntype,numat,Octupole,popats,pos,Quadrupole,rot_atom_abs,Rot_int, &
               Spinorbite,State_all,Symmol,Tensor_imp)

          allocate( poldafse(3,npldafs,nphim) )
          allocate( poldafss(3,npldafs,nphim) )
          allocate( phdafs(natomsym,npldafs) )
          allocate( phdf0t(npldafs,nphim) )
          allocate( phdt(npldafs,nphim) )
          allocate( vecdafse(3,npldafs,nphim) )
          allocate( vecdafss(3,npldafs,nphim) )

          ngrm = 0
          do ipr = 1,n_atom_proto
            ngrm = max( ngrm, ngreq(ipr) )
          end do

          allocate( ltypcal(nplrm) )
          allocate( nomabs(ncolm) )
          allocate( pol(3,nplrm) )
          allocate( vec(3,nplrm) )
          allocate( pdp(nplrm,2) )

! Evaluation des polarisations et vecteurs d'onde.
          call polond(Base_spin,Dipmag,Green_plus,icheck(6),ltypcal,Moyenne,mpirank0,msymdd, &
            msymqq,ncolm,ncolr,nomabs,nple,nplr,nplrm,nxanout,Octupole,Orthmat,Orthmati,pdp, &
            pdpolar,pol,polar,Polarise,Quadrupole,Rot_int,veconde,vec,Xan_atom)

          ncolt = ncolr

          allocate( Taux_ipr(n_atom_proto) )
          do ipr = 1,n_atom_proto
            if( Taux ) then
              Taux_ipr(ipr) = sum( Taux_oc(abs( igreq(ipr,1:ngreq(ipr)) )) ) / ngreq(ipr)
            else
              Taux_ipr(ipr) = 1._db
            endif
          end do
          f_avantseuil = f_cal(Doping,Eseuil,Full_self_abs,icheck(6),itypepr,n_atom_proto,nbseuil,ngreq,ntype,numat,numat_abs, &
             Self_abs,Taux_ipr,volume_maille,xsect_file)
          deallocate( Taux_ipr )

          if( Dafs ) then
            call prepdafs(Angle_or,angpoldafs,angxyz,Axe_atom_gr,axyz,Base_spin,Bormann,Dafs_bio,Eseuil,f_no_res,hkl_dafs, &
              icheck(6),igreq,iprabs_nonexc,isigpi,itabs,itypepr,lvval,Magnetic,Mat_or,mpirank0,n_atom_proto,natomsym,nbseuil, &
              neqm,ngreq,ngrm,ngroup,ngroup_m,ngroup_taux,ngroup_temp,nlat,nlatm,nphi_dafs,nphim,npldafs, &
              nrato,nrm,nspin,ntype,numat,Orthmat,Orthmati,phdafs,phdf0t,phdt,poldafse,poldafsem,poldafss,poldafssm, &
              popatm,posn,psival,rato,Rot_int,Taux,Taux_oc,temp,Temp_coef,Temperature,Vec_orig,vecdafse,vecdafsem, &
              vecdafss,vecdafssm,xsect_file)

            call col_dafs_name(angpoldafs,Bormann,Full_self_abs,hkl_dafs,isigpi,mpirank0,ncolm,ncolr,ncolt, &
                               nomabs,npldafs,Self_abs)
          endif

          if( mpirank0 == 0 ) allocate(Int_tens(n_tens_max*ninitlr,0:natomsym))

          if( Recup_optic ) exit boucle_coh

        endif

        if( Extract ) then

          deallocate( karact )

          if( Tddft ) then
            Final_tddft = .true.
            deallocate( Epsii )
            ninitlr = 1
            allocate( Epsii(ninitlr) )
            Epsii(1) = Epsii_moy
          elseif( Core_resolved ) then
            if( Optic ) then
           ! En optique il y a une sortie par "l" plus le total
              if( numat_abs <= 4 ) then
                ninitlr = 1
              elseif( numat_abs > 4 .and. numat_abs <= 20 ) then
                ninitlr = 2
              elseif( numat_abs > 21 .and. numat_abs <= 56 ) then
                ninitlr = 3
              else
                ninitlr = 4
              endif
              ninitlr = ninitlr + 1
            else
               ninitlr = ninitl
            endif
          else
            ninitlr = nbseuil
          endif

          call extract_write_coabs(Allsite,Ang_rotsup,angxyz,axyz,Base_spin,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
            Densite_atom,E_cut,Energ_s,Energphot,Epsii,Eseuil,Final_tddft, &
            f_avantseuil,Full_self_abs,Green_i,Green_plus,hkl_dafs,iabsorig(multi_run),icheck(21),Int_tens,isigpi,isymeq, &
            nsymextract(multi_run),jseuil,ltypcal,Moyenne,mpinodee,multi_imp,Multipole,n_multi_run,n_oo,n_tens_max, &
            natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ninit1,ninitlr,nomabs,nomfich,nomfich_cal_conv(multi_run), &
            nomfich_cal_convt,nom_fich_extract,nomfich_s,nphi_dafs,nphim,npldafs,nplr,nplrm,nseuil,nspin,numat_abs,nxanout,pdp, &
            phdafs,phdf0t,phdt,pol,Polarise,poldafse,poldafss,Rot_int,Self_abs,Spherical_signal,Spherical_tensor,Spinorbite, &
            Taux_eq,Tddft,V0muf,vecdafse,vecdafss,Vec,Volume_maille,Xan_atom)

          deallocate( Energ_s, Eimag_s )
          deallocate( Axe_atom_clu, dista, iapot, iaproto )
          deallocate( igroup, isymeq, itypep, itypepr, pos, poseq )
          deallocate( Taux_eq, Int_tens, poldafse, poldafss )
          deallocate( phdafs, phdf0t, phdt, vecdafse, vecdafss )
          deallocate( ltypcal, nomabs, pdp, pol, vec )

          cycle boucle_multi

        endif

! Calcul du nombre d'atomes du petit agregat symmetrise
        natome = natome_cal(igrpt_nomag,iopsymr,mpirank0,natomeq,natomp,pos)

! Le nombre d'atomes decrivant le potentiel est defini soit par les  atomes prototypiques,
! soit par les atomes a l'interieur du petit agregat
        if( i_self == 1 .or. Second_run ) Full_atom = ( Full_atom_e &
         .or. ( Self_cons .and. .not. (Self_nonexc .and. Proto_all) ) .or. ( natome <= n_atom_proto + 1 ) .or. Full_potential &
         .or. Hubbard ) .and. .not. Flapw

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

        if( mpirank0 == 0 ) then
          allocate( Int_statedens(lla2_state,nspinp,n_atom_0:n_atom_ind) )
          Int_statedens(:,:,:) = 0._db
        endif

        if( Self_cons .and. .not. Second_run ) then
          nrm_self = nrm
        else
          nrm_self = 0
        endif

        if( i_self == 1 .or. Second_run ) then
          if( Self_nonexc .or. Full_atom ) then
            n_atom_0_self = 1
          else
            n_atom_0_self = 0
          endif
          n_atom_ind_self = n_atom_ind
          natome_self = natome
          natomeq_self = natomeq
          allocate( chargat_init(n_atom_0_self:n_atom_ind_self,nspin))
          allocate( chargat_self(n_atom_0_self:n_atom_ind_self,nspin))
          allocate( chargat_self_s(n_atom_0_self:n_atom_ind_self,nspin) )
          allocate( drho_ex_nex(nrm,nspin) )
          allocate( dvc_ex_nex(nrm) )
          allocate( dv_ex_nex(nrm) )
          allocate( Energ_self(n_atom_0_self:n_atom_ind_self) )
          allocate( Energ_self_s(n_atom_0_self:n_atom_ind_self) )
          allocate( En_coeur(n_atom_0_self:n_atom_ind_self) )
          allocate( En_coeur_s(n_atom_0_self:n_atom_ind_self) )
          allocate( pop_orb_val(n_atom_0_self:n_atom_ind_self,nspin) )
          allocate( rho_chg(0:nrm_self,nspin,n_atom_0_self:n_atom_ind_self) )
          allocate( rho_self(0:nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self) )
          allocate( rho_self_t(0:nrm_self,nlm_pot,n_atom_0_self:n_atom_ind_self) )
          allocate( rho_self_s(0:nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self) )
          allocate( rhoato_init(0:nrm_self,nspin,n_atom_0_self:n_atom_ind_self) )
          allocate( vcato_init(0:nrm_self,n_atom_0_self:n_atom_ind_self) )

          allocate( occ_hubb(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self) )
          allocate( occ_hubb_i(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self) )
          allocate( V_hubb(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self) )
          allocate( V_hubb_s(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp,n_atom_0_self:n_atom_ind_self) )
          allocate( Hubb_diag(n_atom_0_self:n_atom_ind_self) )
          if( Hubbard .and. .not. Cal_xanes ) then
            V_hubb(:,:,:,:,:) = (0._db, 0._db)
            V_hubb_s(:,:,:,:,:) = (0._db, 0._db)
          endif
          Hubb_diag(:) = .true.

        endif

        allocate( dV0bdcF(nspin) )

        if( i_self > 1 .and. Cal_xanes .and. .not. Second_run ) then
          deallocate( Atom_axe )
          deallocate( distai )
          deallocate( ia_eq, ia_eq_inv, ia_rep, iaprotoi, igroupi, iopsym_atom, is_eq, itypei )
          deallocate( nb_eq, nb_rpr, nb_rep_t )
          posi_self(:,:) = posi(:,:)
          deallocate( posi )
          deallocate( rot_atom )
        endif

        if( i_self == 1 .or. Second_run ) allocate( posi_self(3,natome) )

        if( i_self == 1 .or. Cal_xanes ) then

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
          allocate( nb_rpr(natome,natome) )
          allocate( nb_rep_t(nb_sym_op,natome,natome) )
          allocate( posi(3,natome) )
          allocate( rot_atom(3,3,natome) )
          Energ_self_s(:) = 0._db

          if( Cal_xanes ) then
            ich = icheck(7)
          elseif( icheck(27) > 1 ) then
            ich = max( icheck(7), icheck(27) )
          else
            ich = 0
          endif

          call Atom_selec(adimp,Atom_axe,Atom_with_axe,Atom_nonsph,Atom_occ_hubb,Axe_atom_clu,Base_ortho,dcosxyz, &
            dista,distai,Full_atom,Green,Hubbard,i_self,ia_eq,ia_eq_inv,ia_rep,iaabs,iaabsi,iaproto,iaprotoi,ich,igreq, &
            igroup,igroupi,igrpt_nomag,igrpt0,iopsym_atom,iopsymr,iord,is_eq,itype,itypei,itypep,itypepr,Magnetic,m_hubb, &
            m_hubb_e,mpirank0,natome,n_atom_0_self,n_atom_ind_self,n_atom_proto,natomeq,natomp,nb_eq,nb_rpr, &
            nb_rep_t,nb_sym_op,neqm,ngroup,ngroup_hubb,ngroup_m,nlat,nlatm,nspin,nspinp,ntype,numat,nx,occ_hubb_e,Overad,popats, &
            pos,posi,rmt,rot_atom,roverad,rsort,rsorte,Spinorbite,Symmol,V_hubb,V_hubbard,Ylm_comp,Ylm_comp_inp)

        end if

        if( Cal_xanes .and. Hubb(itabs) ) then
          if( Full_atom ) then
            iapr = iaabsi
          elseif( Self_nonexc ) then
            iapr = iprabs_nonexc
          else
            iapr = iprabs
          endif
          V_hubb_abs(:,:,:,:) = V_hubb(:,:,:,:,iapr)
        endif

        if( i_self == 1 .or. Second_run ) then
          allocate( ia_eq_inv_self(natomeq_self) )
          ia_eq_inv_self(:) = ia_eq_inv(:)
        endif

        Tddft_xanes = Tddft .and. Cal_xanes
        Optic_xanes = Optic .and. Cal_xanes

! Calcul des representations utiles
        Recop = .false.
        if( Cal_xanes .and. .not. Optic ) then
          call etafin(Ylm_comp,icheck(8),iopsymr,irep_util,jseuil,karact,ldip,lmoins1,loct,lplus1,lqua,lseuil, &
                mpirank0,nb_rep,nbseuil,ngrph,nspino,Spinorbite,State_all,Sym_cubic)
          if( ( Sym_cubic .or. Sym_4 ) .and. .not. State_all ) Recop = .true.
        elseif( i_self == 1 ) then
! Pour le calcul auto-coherent, il faut calculer toutes les representations
          if( icheck(27) > 1 ) then
            ich = max( icheck(8), icheck(27) )
          else
            ich = 0
          endif
          call irep_util_all(ich,iopsymr,irep_util,karact,nb_rep,ngrph,nspino,Spinorbite)
        endif
        if( i_self > 1 .and. Cal_xanes .and. .not. Second_run ) deallocate( Repres_comp )
        if( i_self == 1 .or. Cal_xanes ) then
           allocate( Repres_comp(ngrph) )
           call irep_comp(Ylm_comp,Green,iopsymr,irep_util,karact,ngrph,nspino,Repres_comp)
        end if

        Moy_loc = Green .and. .not. Moy_cluster

        if( ( i_self == 1 .or. Cal_xanes ) .and. .not. second_run ) then

          nbm = 0;     nbtm = 0
          npoint = 0;  npr = 0
          nsort = 0
          nsm = 0;     nstm = 0

! Calcul des dimensions de tableaux pour le maillage

          call nbpoint(Adimp,Base_hexa,Base_ortho,D_max_pot,dcosxyz,Green,iaabsfirst,igrpt_nomag,iopsymr,iord,Moy_loc, &
                       mpirank0,natomp,npoint,npso,nvois,nx,pos,rsort)

          if( Atom_nonsph ) then
            npoint_ns = npoint
          else
            npoint_ns = 0
          endif
          nim = npoint
          npsom = npso

          if( Cal_xanes .and. i_self > 1 ) deallocate( clapl, cgrad, ivois, isvois, numia, rvol, xyz )

          allocate( clapl(0:nvois) )
          allocate( cgrad(nvois) )
          allocate( ivois(npsom,nvois) )
          allocate( isvois(npsom,nvois) )
          allocate( numia(npsom) )
          allocate( rvol(nim) )
          allocate( xyz(4,npsom) )

! Elaboration du maillage

          allocate( indice(npsom,3) )
          allocate( mpres(-nx:nx,-nx:nx,-nx:nx) )

! Meme en Green, on definit des points afin de calculer le potentiel moyen
          call reseau(Adimp,Base_hexa,Base_ortho,D_max_pot,dcosxyz,Green,iaabsfirst,icheck(9), &
                 igrpt_nomag,indice,iopsymr,iord,itypei,Moy_loc,mpirank0,mpres,natome,natomp,nim, &
                 npoint,npr,npso,npsom,ntype,numia,nx,pos,posi,rmt,rsort,rvol,xyz)

          if( .not. Green ) then
            call laplac(Adimp,Base_hexa,cgrad,clapl,icheck(9),igrpt_nomag,indice,iopsymr,iord,ivois,isvois, &
                  mpirank0,mpres,npso,npsom,nvois,nx)
            call bordure(Base_ortho,dcosxyz,Green,icheck(9),iopsymr,iord,iscratch,ivois,mpirank0,natome,nbm,nbtm,nim, &
                  npoint,npso,npsom,nsm,nstm,numia,nvois,posi,rvol,xyz)
          endif
          deallocate( indice, mpres )

        endif

! On calcule les rayons seulement aux premiere et derniere iterations

        if( i_self == 1 .or. Second_run ) then
          allocate( V_abs_i(nrm,nspin) )
          allocate( rhoato_abs(nrm,nspin) )
        endif

        if( .not. Second_run ) then
          if( i_self == 1 ) then
            allocate( lmaxat(0:n_atom_proto) )
            allocate( rmtg(0:n_atom_proto) )
            allocate( rmtg0(0:n_atom_proto) )
            allocate( rmtsd(0:n_atom_proto) )
          end if
          allocate( rs(npoint) )
          allocate( rsato(0:nrm,n_atom_0:n_atom_ind) )
          allocate( Vcato(0:nrm,nlm_pot,n_atom_0:n_atom_ind) )
          allocate( Vh(npoint) )
          allocate( Vhns(npoint_ns) )
          allocate( Vxc(npoint,nspin) )
          allocate( Vxcato(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) )
          rsato(:,:) = 0._db
          Vcato(:,:,:) = 0._db
          Vxcato(:,:,:,:) = 0._db
        endif

        allocate( Excato(nrm,n_atom_0:n_atom_ind) )
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

        if( .not. Green ) call recup_bordure(ibord,isbord,iscratch,isrt,mpinodes0,mpirank0,natome,nbord,nbordf,nbm,nbtm,nsm, &
            nsort,nsortf,nstm,poidsa,poidso)

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          tp(2) = real(time,db)
        endif

        if( .not. Second_run ) then

            if( i_self == 1 ) then
              allocate( chg_cluster(nspin) )
              allocate( E_starta(n_atom_0_self:n_atom_ind_self) )
              allocate( ispin_maj(n_atom_0_self:n_atom_ind_self) )
            endif

! Calcul du potentiel
          if( Flapw ) then

            call potlapw(axyz,Base_ortho,chargat,Coupelapw,dcosxyz,deccent,Flapw_new,Full_atom,iapot, &
              iaproto,iaprotoi,icheck(13),igroup,iprabs_nonexc,ipr1,itabs,its_lapw,itypei,itypep,itypepr,Magnetic,mpinodes0, &
              mpirank0,n_atom_0,n_atom_ind,n_atom_proto,natome,natomeq,natomp,ngreq,ngroup,ngroup_lapw,nklapw,nlm_pot, &
              nlmlapwm,nmatsym,normrmt,npoint,npsom,nrato,nrato_lapw,nrm,nslapwm,nspin,ntype, &
              numat,Orthmat,overlap,pos,rato,rchimp,rho,rlapw,rmtg,rmtg0,rmtimp,rmtsd,Rot_int,rotloc_lapw,rs,rsato,rsort, &
              Trace_format_wien,Trace_k,Trace_p,V_abs_i,V_intmax,V0bdcFimp(1),Vcato,Vh,Vxc,Vxcato,Wien_file,Wien_matsym, &
              Wien_save,Wien_taulap,xyz)

          else

            call potsup(alfpot,Atom_nonsph,Axe_atom_gr,Base_ortho,Cal_xanes,cdil,chargat,chargat_init, &
              chargat_self,dcosxyz,drho_ex_nex,dv_ex_nex,dvc_ex_nex,Excato,Full_atom, hybrid,i_self, &
              ia_eq_inv,ia_eq_inv_self,iaabs,iaproto,iaprotoi,iapot,icheck,igreq,igroup,iprabs_nonexc,iprabs,ipr1,itab,itdil, &
              itypei,itypep,itypepr,ldil,lmax_pot,lvval,Magnetic,mpirank0,n_atom_0,n_atom_0_self,n_atom_ind, &
              n_atom_ind_self,n_atom_proto,natome,natome_self,natomeq,natomeq_self,natomp,neqm,ngreq,ngroup_m,ngroup_nonsph, &
              nhybm,nlat,nlatm,nlm_pot,Nonexc,norbdil,norbv,normrmt,npoint,npoint_ns,npsom,nrato,nrm,nrm_self,nspin,ntype, &
              numat,overlap,pop_nonsph,popatm,popatv,pos,posi,posi_self,psival,r_self,rato,rchimp,rho,rho_chg, &
              rho_self,rhoato_abs,rhoato_init,rhoit,rhons,Rmtg,Rmtimp,Rmtg0,Rmtsd,Rot_Atom_gr,Rot_int,rs, &
              rsato,rsort,Self_nonexc,Tddft_xanes,V_abs_i,V_intmax,Vcato,Vcato_init,Vh,Vhns,Vsphere,Vxc,Vxcato,V0bdcFimp(1),xyz)

! Initialisation de rho_self, une fois avoir calcule le nouveau potentiel en potsup.
! Au dela du rayon rmtsd on initialise a zero.
            if( .not. Cal_xanes ) then

              if( i_self == 1 ) then
                allocate( ch_coeur(n_atom_0_self:n_atom_ind_self) )

                rho_self_s(:,1,:,:) = rhoato_init(:,:,:)

! Calcul de l'energie de depart du comptage d'electrons; on la calcule  une seule fois,
! bien que les Vcato et Vxcato changent a chaque iteration
                call En_dep(E_start,E_starta,Full_atom,iaprotoi,icheck(27),itypepr,lcoeur,n_atom_0,n_atom_0_self, &
                   n_atom_ind,n_atom_ind_self,n_atom_proto,natome,ncoeur,nenerg_coh,nlm_pot,nrato,nrm,numat,nspin, &
                   ntype,Pas_SCF,psi_coeur,Relativiste,rato,Rmtg,Rmtsd,V_intmax,Vcato,Vxcato,workf)

! Calcul de la charge de reference en cas de calcul auto_coherent
                call chg_agr(chargat,chargat_init,ch_coeur,chg_cluster,chg_open_val,Doping,Full_atom,iaprotoi,ipr_dop, &
                   iprabs_nonexc,ispin_maj,itabs,icheck(27),itypepr,mpirank0,natome,n_atom_0_self,n_atom_ind_self, &
                   n_atom_proto,nb_eq,ngreq,nrato,nrm,nrm_self,nspin,ntype,numat,pop_open_val, &
                   psi_open_val,rato,rho_chg,rho_coeur,rhoato_init,rmtsd,SCF_mag_fix,SCF_mag_free)

                chargat_self_s(:,:) = chargat_init(:,:)
! Mise en place de la grille en energie pour l'autocoherence
                allocate( Energ_coh(nenerg_coh) )
                allocate( Eimag_coh(nenerg_coh) )
                call grille_coh(Eimag_coh,Energ_coh,E_start,Green,icheck(27),nenerg_coh,Pas_SCF)

              endif

              do iapr = n_atom_0_self,n_atom_ind_self
                if( Full_atom ) then
                  ipr = iaprotoi(iapr)
                else
                  ipr = iapr
                endif
                it = itypepr(ipr)
                nr = nrato(it)
                do n = 1,nr
                  if( rato(n,it) > Rmtsd(ipr) ) exit
                end do
                do ispin = 1,nspin
                  do ir = 0,n
                    rho_self(ir,1,ispin,iapr) = rho_chg(ir,ispin,iapr) + ( rho_coeur(ir,it) / nspin )
                  end do
                  rho_self(n+1:nr,1,:,iapr) = rhoato_init(n+1:nr,:,iapr)
                end do
              end do

            endif

          endif
        endif

        if( Full_atom ) then
          iaprabs = iaabsi
        else
          iaprabs = iprabs_nonexc
        endif
! On stocke le potentiel non excite de l'absorbeur (sert au calcul de l'energie du niveau initial )
        if( Second_run ) then
          do ispin = 1,nspin
            V_abs_i(1:nrm,ispin) = Vcato(1:nrm,1,iaprabs) + Vxcato(1:nrm,1,ispin,iaprabs)
          end do
        endif

        if( ipr1 == 1 ) rmtg(0) = rmtg(iprabs_nonexc)

        if( Cal_xanes ) then
          nenerg = nenerg_s
          allocate( Energ(nenerg) )
          allocate( Eimag(nenerg) )
          Energ(:) = Energ_s(:)
          if( Optic ) then
            if( E_Fermi_man ) then
              Energ(:) = Energ(:) + E_cut_imp
            else
              Energ(:) = Energ(:) + E_cut
            endif
          endif
          Eimag(:) = eimag_s(:)
        else
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
        end if

        allocate( V0bdcF(nspin) )
! Operations complementaires sur le potentiel
        call potential_comp(Base_ortho,Cal_xanes,dcosxyz,distai(natome),dV0bdcF,ecineticmax,ecineticmax_out, &
            Eclie,Eneg,Energ(nenerg),Green,iaabs,iaproto,icheck(13),imoy,imoy_out,iopsymr,isrt,korigimp,Magnetic, &
            Moy_loc,mpirank0,n_atom_proto,natomp,nim,npoint,npsom,nptmoy,nptmoy_out,nsortf,nspin,nstm,poidsov,poidsov_out,pos, &
            rmtg0,rs,rsbdc,rsbdc_out,rsort,rvol,V0bdcF,V0bdcFimp,V0muf,Vh,Vhbdc,Vhbdc_out,Vr,Vxc,VxcbdcF,VxcbdcF_out,xyz,Workf)

! Calcul de l'energie du niveau de coeur initial.
        if( Cal_xanes ) then
          if( Optic ) then
            Epsii(:) = 0._db; Epsii_moy = 0._db
          else
            call Energseuil(Core_resolved,Delta_Epsii,Delta_Eseuil,Epsii,Epsii_moy,Eseuil,icheck(14),is_g, &
              itabs_nonexc,lseuil,m_g,mpirank0,nbseuil,ninit1,ninitl,ninitlr,nrato(itabs_nonexc),nrm,nseuil,nspin,ntype, &
              numat_abs,psii,rato,Rmtg(iprabs_nonexc),Rmtsd(iprabs_nonexc),V_abs_i,V_intmax,V0bdcf)
          endif
        endif

        deallocate( V0bdcF )

! Determination du lmax
        if( .not. Green ) then
          call clmax(Ecineticmax_out,rsort,lmaxso0,lmaxso_max,2,.true.)
        else
          lmaxso_max = 0
        endif
        lmaxmax = 0
        if( Optic ) then
          if( numat_abs < 21 ) then
            lmax_probe = 2
          else
            lmax_probe = 3
          endif
          if( lmaxat0 > - 1 ) lmax_probe = min( lmax_probe, lmaxat0 )
        else
          if( Octupole ) then
            lmax_probe = lseuil + 3
          elseif( Quadrupole ) then
            lmax_probe = lseuil + 2
          else
            lmax_probe = lseuil + 1
          endif
        endif
        nlm_probe = (lmax_probe + 1)**2

        lmaxabs = 0
        do ipr = 0,n_atom_proto
          Z = numat( itypepr(ipr) )
          call clmax(Ecineticmax,rmtg(ipr),lmaxat0,lmaxat(ipr),Z,lmaxfree)
          if( Z == numat_abs ) then
            if( Cal_xanes ) lmaxat(ipr) = max( lmax_probe, lmaxat(ipr) )
            lmaxabs = max( lmaxat(ipr), lmaxabs )
          endif
          lmaxmax = max(lmaxmax,lmaxat(ipr))
        end do
        do ipr = 0,n_atom_proto
          Z = numat( itypepr(ipr) )
          if( Z == numat_abs ) lmaxat(ipr) = lmaxabs
        end do

        nlmmax = (lmaxmax + 1 )**2
        nlmamax = (lmaxabs + 1 )**2
        if( Tddft ) then
! nlmamax ne sert que pour la Tddft, on la limite pour espace memoire:
          lmaxabs_t = min( lmaxabs, 4)
          nlmamax = (lmaxabs_t + 1 )**2
          if( iopsymc(25) == 0 .or. Quadrupole .or. Dipmag ) then
            imparite = 2
          elseif( mod(lseuil,2) == 0 ) then
            imparite = 1
          else
            imparite = 0
          endif
          nlmamax_u = 0
          do l = 0,lmaxabs_t
            if( (imparite /= mod(l,2)) .and. (imparite /= 2) ) cycle
            nlmamax_u = nlmamax_u + 2 * l + 1
          end do
        else
          lmaxabs_t = 0
          nlmamax = 0
          nlmamax_u = 0
          imparite = 2
        endif

        nlmsam = nspinp * nlmmax
        nlmomax = ( lmaxso_max + 1 )**2
        nso1 = nspino * nlmomax

        if( mpirank0 == 0 ) then
          call CPU_TIME(time)
          tp(3) = real(time,db)
        endif

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
               lmaxat,lmaxso_max,lso,mato,mpirank0,mso,n_atom_proto,natome,nbordf,ngrph,nlmsa0,nlmsam,nlmso0, &
               nso1,nsortf,nspino,ntype,numat,Orthmati,posi,rot_atom)

! Calcul des Ylm sur les points du maillage
        if( Green ) then
          nbtm_fdm = 0
        else
          nbtm_fdm = nbtm
        endif
        allocate( Ylmato(nbtm_fdm,nlmmax,natome) )
        if( .not. Green )  then
          allocate( Ylmso(nsort,nlmomax) )
          call ylmpt(Base_ortho,dcosxyz,iaprotoi,ibord,icheck(15),iopsymr,isrt,lmaxat,lmaxso_max,n_atom_proto,natome, &
               nbord,nbtm,nlmmax,nlmomax,nsort,nstm,npsom,posi,rot_atom,rot_atom_abs,xyz,Ylmato,Ylmso)
        endif

        if( nspin == 1 .or. Spinorbite ) then
          nspinorb = 1
        else
          nspinorb = 2
        endif

        allocate( Ecinetic_out(nspin) )
        allocate( lmaxa(natome) )
        allocate( nlmsa(natome) )
        allocate( V0bdc_out(nspin) )
        if( .not. Green ) then
          if( Spinorbite .or. Relativiste ) then
            nicm = nim
          else
            nicm = 1
          endif
          allocate( gradvr(nicm,3,nspin))
        endif
        allocate( drho_self(nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self,0:mpinodes-1) )
        allocate( Statedens(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1) )
        allocate( Statedens_i(lla2_state,nspinp,lla2_state,nspinp,n_atom_0:n_atom_ind,0:mpinodes-1) )
        allocate( sec_atom(ninitlr) )

! Initialisation de l'energie de l'agregat
        Enragr = 0._db
        Energ_self(:) = 0._db

        if( Cal_xanes ) then
          if( Tddft_xanes ) then
            nenerg_tddft = nenerg
          else
            nenerg_tddft = 0
          endif
          allocate( rof0(nenerg_tddft,nlmamax,nspinp,nspino,nbseuil) )
        endif
        if( Tddft_xanes ) then
          rof0(:,:,:,:,:) = (0._db,0._db)
          allocate( imag_taull(nenerg,nlmamax_u,nspinp,nlmamax_u,nspinp) )
          imag_taull(:,:,:,:,:) = 0._db
        end if
        if( Optic_xanes ) then
          allocate( Taull_opt(nenerg,nlm_probe,nlm_probe,nspinp,nspinp) )
          Taull_opt(:,:,:,:,:) = (0._db,0._db)
        end if

        nbbseuil = nbseuil
        ninitlv = nbseuil
        n_vr_0 = n_atom_0
        n_vr_ind = n_atom_ind

        allocate( secdd(3,3,ninitlr,0:mpinodes-1) )
        allocate( secmd(3,3,ninitlr,0:mpinodes-1) )
        allocate( secmm(3,3,ninitlr,0:mpinodes-1) )
        allocate( secdq(3,3,3,ninitlr,0:mpinodes-1) )
        allocate( secdo(3,3,3,3,ninitlr,0:mpinodes-1) )
        allocate( secoo(3,n_oo,3,n_oo,ninitlr,0:mpinodes-1) )
        allocate( secqq(3,3,3,3,ninitlr,0:mpinodes-1) )
        allocate( secdd_m(3,3,ninitlr,0:mpinodes-1) )
        allocate( secmd_m(3,3,ninitlr,0:mpinodes-1) )
        allocate( secmm_m(3,3,ninitlr,0:mpinodes-1) )
        allocate( secdq_m(3,3,3,ninitlr,0:mpinodes-1) )
        allocate( secdo_m(3,3,3,3,ninitlr,0:mpinodes-1) )
        allocate( secoo_m(3,n_oo,3,n_oo,ninitlr,0:mpinodes-1) )
        allocate( secqq_m(3,3,3,3,ninitlr,0:mpinodes-1) )
        allocate( Ecinetic(nspin) )
        allocate( V0bdc(nspin) )

        Green_i = Green_int

! Boucle sur l'energie
        nge = ( nenerg - 1 ) / mpinodes + 1

        if( One_run .and. mpinodes0 > 1 .and. multi_run == 1 .and. Cal_xanes ) &
          allocate( taull_stk(nlm_probe,nspinp,nlm_probe,nspinp,2:n_multi_run,nge) )
        index_e = 0

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
          tp(4) = real(time,db)
          do i = 2,4
            tpt(i) = tp(i) - tp(i-1)
          end do
        endif

        boucle_energ: do je = 1,nge

          if( Recup_tddft_data .and. Cal_xanes ) exit

          ie = ( je - 1 ) * mpinodes + mpirank + 1

          if( ie <= nenerg ) then

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              tp(1) = real(time,db)
            endif

            index_e = index_e + 1
            Enervide = Energ(ie) - Workf

! Calcul du potentiel dans l'etat excite
            call potex(Atom_nonsph,axyz,alfpot,Cal_xanes,dV0bdcF,Energ(ie),Enervide,Final_tddft,Full_atom, &
              iaabsi,iapot,iaprotoi,icheck(16),1,iprabs_nonexc,iprabs,itab,itypepr,Magnetic, &
              n_atom_0,n_atom_ind,n_atom_proto,n_vr_0,n_vr_ind,natome,nlm_pot,Nonexc_g,npoint,npoint_ns,npsom,nptmoy, &
              nptmoy_out,nrato,nrm,nspin,ntype,rato,rho,rhons,rmtg,rs,rsato,rsbdc,rsbdc_out,Trace_k,Trace_p, &
              V_intmax,Vcato,Vh,Vhbdc,Vhbdc_out,Vhns,Vr,Vxc,VxcbdcF,VxcbdcF_out,Vrato,Vxcato,V0bdc,V0bdc_out,xyz)

            do ispin = 1,nspin

              V0bdc_out(ispin) = V0bdc(ispin)

              if( Muffintin ) call mdfmuf(Atom_nonsph,Axe_Atom_gr,axyz,Base_ortho, &
                dcosxyz,Full_atom,iaabs,iaproto,iaprotoi,icheck(16),igreq,igroup,ispin,itypep,n_atom_0,n_atom_ind, &
                n_atom_proto,natome,natomp,neqm,ngroup_m,nlm_pot,npoint,npoint_ns,npsom,nrato,nrm,nspin,ntype,pos,rato,rho,rhons, &
                Rmtg,Trace_k,Trace_p,Vhns,V0bdc(ispin),Vr,Vrato,xyz)

              if( Supermuf ) call modmuf(Full_atom,iaprotoi,itypepr,icheck(16), &
                ispin,n_atom_0,n_atom_ind,n_atom_proto,natome,nlm_pot,nrato,nrm,nspin,ntype,rato,rmtg,V0bdc(ispin),Vrato)

              if( .not. ( Green .or. Second_run ) .and. ( Spinorbite .or. Relativiste ) ) &
                call gradpot(Base_hexa,cgrad,gradvr,icheck(16),iord,ispin,ivois,nicm,npoint,npsom,nspin,nvois,Vr)

              Ecinetic_out(ispin) = Enervide - V0bdc_out(ispin)
              Ecinetic(ispin) = Enervide - V0bdc(ispin)

              if( Eneg ) then
! Problemes numeriques autour de 0
                em = 0.01_db / rydb
                ei0 = 0.01_db / rydb
                if( Ecinetic_out(ispin) < 0._db ) then
                  Eimag(ie) = max( Eimag(ie), ei0 )
                elseif( ecinetic_out(ispin) < em ) then
                  eii = ei0 * ( 1 - ecinetic_out(ispin) / em )
                  Eimag(ie) = max( Eimag(ie), eii )
                endif
              else
                Ecinetic_out(ispin) = max( Ecinetic_out(ispin), Eclie)
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

            end do  ! fin boucle sur spin

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              tp(2) = real(time,db)
              tpt(5) = tpt(5) + tp(2) - tp(1)
            endif

            if( abs( Eimag(ie) ) > eps10 .or. Ecinetic_out(1) < eps10 .or. Ecinetic_out(nspin) < eps10 ) then
              E_comp = .true.
            else
              E_comp = .false.
            endif

            if( No_solsing ) then
              Solsing = .false.
            elseif( E_comp ) then
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

            lmaxg = 0
            lmaxabsa = 0
            do ipr = 0,n_atom_proto
              Z = numat( itypepr(ipr) )
              call clmax(Ecmax,Rmtg(ipr),lmaxat0,lmaxat(ipr),Z,lmaxfree)
              lmaxat(ipr) = min(lmaxat(ipr), lmaxmax)
              if( Z == numat_abs ) then
                 if( Cal_xanes ) lmaxat(ipr) = max( lmax_probe, lmaxat(ipr) )
                 lmaxabsa = max( lmaxabsa, lmaxat(ipr) )
              end if
              lmaxg = max( lmaxat(ipr), lmaxg)
            end do
            do ipr = 0,n_atom_proto
              Z = numat( itypepr(ipr) )
              if( Z == numat_abs ) lmaxat(ipr) = lmaxabsa
            end do
            nlmagm = ( lmaxg + 1 )**2
            if( Tddft ) then
              nlmam_u = 0
              lmax_t = min(lmaxabsa,lmaxabs_t)
              do l = 0,lmax_t
                if( (imparite /= mod(l,2)) .and. (imparite/=2) ) cycle
                nlmam_u = nlmam_u + 2 * l + 1
              end do
            endif

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
            phiato(:,:,:,:,:,:) = 0._db
            Tau_ato(:,:,:,:,:) = (0._db,0._db)
            Taull(:,:,:,:,:) = (0._db,0._db)

            if( .not. Green ) then
              Z = 2
              call clmax(Ecmax_out,rsort,lmaxso0,lmaxso,Z,.true.)
              if( Cal_xanes .and. icheck(18) > 0 ) write(3,'(/a9,i4)') ' lmaxso =', lmaxso
              lmaxso = min( lmaxso, lmaxso_max )
            endif

            if( Cal_xanes ) then
              ich = icheck(18)
            else
              ich = min(icheck(18), icheck(27))
            endif
            if( ich > 1 ) write(3,190)

! Boucle sur les atomes nonequivalents dans la maille
            if( Cal_xanes ) then
              n1 = n_atom_0
            else
              n1 = n_atom_0_self
            endif

            do iapr = n1,n_atom_ind

              if( Full_atom ) then
                ipr = iaprotoi( iapr )
              else
                ipr = iapr
! pour sauter l'absorbeur non excite et les atomes non calcule par SCF
! Commente car bug dans compilateur linux !
!              if( Cal_xanes ) then
!                do ia = 1,natome
!                  if( iaprotoi(ia) == iapr ) exit
!                end do
!                if( ia == natome+1 ) cycle
!              endif
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
                Absorbeur = iapr == iaabsi .and. Cal_xanes
              else
                Absorbeur = iapr == iprabs .and. Cal_xanes
              endif

              allocate( V_hubb_t(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp))
              if( Hubb(it) .and. iapr <= n_atom_ind_self ) then
                Hubb_a = .true.
                if( Absorbeur ) then
                  V_Hubb_t(:,:,:,:) = V_Hubb_abs(:,:,:,:)
                else
                  V_Hubb_t(:,:,:,:) = V_Hubb(:,:,:,:,iapr)
                endif
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

              call Sphere(Axe_atom_grn,Base_ortho,dcosxyz,Ecinetic,Eimag(ie),Energ(ie),Enervide,Full_atom, &
              Full_potential,Green,Hubb_a,Hubb_d,iaabsi,iapr,iaprotoi,ibord,ich,igreq,igroupi,iopsymr,lmax,lmax_pot,m_hubb, &
              n_atom_0,n_atom_ind,n_atom_proto,natome,nbord,nbtm,nbtm_fdm,neqm,ngroup_m, &
              nlm_pot,nlmagm,nlmmax,nphiato1,nphiato7,npsom,nr,nspin,nspino,nspinp,Z,phiato,posi,r,Relativiste,Rmtg(ipr), &
              Rmtsd(ipr),Spino,Tau_ato,V_hubb_t,V_intmax,V0bdc,Vrato_e,xyz,Ylm_comp,Ylmato)

              deallocate( r )
              deallocate( Vrato_e, V_hubb_t )

            end do

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              tp(3) = real(time,db)
              tpt(6) = tpt(6) + tp(3) - tp(2)
            endif

! nspinorb = 2 si magnetique sans Spinorbite
! nspinorb = 1 si nonmagnetique ou Spinorbite
! nspino = 2 si Spinorbite

! igrph: boucle sur les representations

            lmaxg = 0
            do ia = 1,natome
              lmaxa(ia) = lmaxat( iaprotoi(ia) )
              lmaxg = max(lmaxg,lmaxa(ia))
            end do

            if( multi_run == 1 .or. .not. One_run ) then

              do ispin = 1,nspinorb
                do igrph = 1,ngrph

                  nlmsamax = 0
                  nlmsa(:) = 0
                  do ia = 1,natome
                    Z = numat(itypei(ia))
                    if( .not. Green .and. Z == 0 ) cycle
                    lmaxa(ia) = lmaxat( iaprotoi(ia) )
                    n = nlmsa0(ia,igrph)
                    allocate( lval(n) )
                    lval(1:n) = lato(1:n,ia,igrph)
                    call cnlmmax(lmaxa(ia),lval,n,nlmsa(ia))
                    nlmsamax = max(nlmsamax,nlmsa(ia))
                    deallocate( lval )
                  end do
                  if( nlmsamax == 0 .and. igrph /= ngrph ) cycle
                  if( Cal_xanes ) then
                    ich = icheck(19)
                  else
                    ich = min(icheck(19), icheck(27))
                  endif
                  if( Green ) then
                    call msm(Axe_Atom_grn,Base_ortho,Cal_xanes,dcosxyz,Ecinetic,Eimag(ie),Full_atom,ia_eq,ia_rep,iaabsi,iaprotoi, &
                      iato,ich,igreq,igroupi,igrph,iopsymr,irep_util,is_eq,ispin,karact,lato,lmaxa,lmaxg,mato,n_atom_0, &
                      n_atom_ind,n_atom_proto,natome,natomp,nb_eq,nb_rpr,nb_rep_t,nb_sym_op,nchemin,neqm,ngroup, &
                      ngrph,nlmagm,nlmsa,nlmsam,nlmsamax,Normaltau,nspin,nspino,nspinp,pos,posi,Recop,Repres_comp(igrph), &
                      rot_atom,Solsing,Spinorbite,State_all,Sym_cubic,Tau_ato,Tau_nondiag,Taull,tpt1,tpt2,Ylm_comp)
                  else
                    n = nlmso0(igrph)
                    allocate( lval(n) )
                    lval(1:n) = lso(1:n,igrph)
                    call cnlmmax(lmaxso,lval,n,nlmso)
                    deallocate( lval )
                    if( Cal_xanes .and. ie > 1 ) ich = ich - 1
                    call mat(Adimp,Atom_axe,Axe_Atom_grn,Base_hexa,Base_ortho,Basereel,Cal_xanes,cgrad, &
                    clapl,dcosxyz,distai,E_comp,Ecinetic_out,Eclie,Eimag(ie),Eneg,Enervide,Full_atom, &
                    gradvr,iaabsi,iaprotoi,iato,ibord,ich,igreq,igroupi,igrph,irep_util,isbord,iso,ispin,isrt,ivois,isvois, &
                    karact,lato,lmaxa,lmaxso,lso,mato,MPI_host_num_for_mumps,mpirank0,mso, &
                    natome,n_atom_0,n_atom_ind,n_atom_proto,nbm,nbord,nbordf,nbtm,neqm,ngroup,ngrph,nim,nicm, &
                    nlmagm,nlmmax,nlmomax,nlmsa,nlmsam,nlmso,nphiato1,nphiato7,npoint,npr,npsom,nsm, &
                    nsort,nsortf,nso1,nspin,nspino,nspinp,nstm,numia,nvois,phiato,poidsa,poidso,R_rydb,Recop,Relativiste, &
                    Repres_comp(igrph),Rsort,rvol,Rydberg,Solsing,Spinorbite,State_all,Sym_cubic,Tau_ato,Taull, &
                    tpt1,tpt2,V0bdc_out,Vr,xyz,Ylm_comp,Ylmato,Ylmso)
                  endif

                end do  ! fin de la boucle sur les representations
              end do  ! fin de la boucle sur le spin

            endif

! C'est bien mpinodes et non mpinodes0 car correspond a la boucle sur l'energie.
            if( One_run .and. Cal_xanes ) call Data_one_run(iabsm,iaprotoi,icheck(19), &
                igreq,index_e,igroupi,lmax_probe,lmaxa,mpinodes0,multi_run,n_atom_proto,n_multi_run,natome, &
                neqm,nge,ngreq,nlm_probe,nlmagm,nspin,Rot_atom,Rot_int,Spinorbite,taull,taull_stk,Ylm_comp)

            if( Tddft_xanes ) then
              lm1 = 0
              do l1 = 0,lmax_t
                if( mod(l1,2) /= imparite .and. imparite /= 2 ) cycle
                do m1 = -l1,l1
                  lm1 = lm1 + 1
                  lmv1 = l1**2 + l1 + 1 + m1
                  lm2 = 0
                  do l2 = 0,lmax_t
                    if( mod(l2,2)/= imparite .and. imparite /= 2 ) cycle
                    do m2 = -l2,l2
                      lm2 = lm2 + 1
                      lmv2 = l2**2 + l2 + 1 + m2
                      imag_taull(ie,lm1,1:nspin,lm2,1:nspin) = - aimag( taull(lmv1,1:nspin,lmv2,1:nspin,iaabsi) )
                    end do
                  end do
                end do
              end do
            endif
            if( Optic_xanes ) then
              do lm1 = 1,nlm_probe
                do lm2 = 1,nlm_probe
                  do isp1 = 1,nspinp
                    Taull_opt(ie,lm1,lm2,isp1,:) = Taull(lm1,isp1,lm2,:,iaabsi)
                  end do
                end do
              end do
              do ipr = 3,6,3
                if( icheck(19) == 0 .and. ipr == 3 ) cycle
                if( ie == 1 .or. ipr == 3 ) write(ipr,'(/A)') '   Energy   -Tau_i(l,m=0,l,m=0) l = 0,lmax_probe'
                write(ipr,'(f10.3,1p,6e13.5)') Energ(ie)*rydb, ( - aimag( Taull_opt(ie,l**2+l+1,l**2+l+1,1,1) ), &
                     l = 0,lmax_probe )

              end do
            endif

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              tp(4) = real(time,db)
              tpt(7) = tpt(7) + tp(4) - tp(3)
            endif

! Calcul des tenseurs cartesiens
            if( Cal_xanes .and. .not. Optic ) then

              nlma = nlm_probe
              ndim1 = 1
              ndim2 = 1
              allocate( taull_abs(nlma*nspino,nlma*nspino,2,2,ndim1,ndim2,ndim2) )

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

                do lm1 = 1,nlma
! Pour la partie atomique, tous le signal est dans la solution singuliere
                  if( i == 1 .and. Solsing ) exit
                  do iso1 = 1,nspino
                    lm1g = nlma*(iso1-1) + lm1
                    do lm2 = 1,nlma
                      do iso2 = 1,nspino
                        lm2g = nlma*(iso2-1) + lm2
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
                              taull_abs(lm1g,lm2g,isp1,isp2,1,1,1) = Tau_ato(lm1,iss1,lm2,iss2,iapr)
                            else
                              taull_abs(lm1g,lm2g,isp1,isp2,1,1,1) = taull(lm1,iss1,lm2,iss2,iaabsi)
                            endif
                          end do
                        end do
                      end do
                    end do
                  end do
                end do

                Hubb_a = Hubb(itabs)
                Hubb_d = Hubb_diag(iaprabs)
                nr = nrato(itabs)
                allocate( Vra(nr,nlm_pot,nspin) )
                allocate( r(nr) )
                if( Full_atom ) then
                  Vra(1:nr,:,:) = Vrato(1:nr,:,:,iaabsi)
                else
                  Vra(1:nr,:,:) = Vrato(1:nr,:,:,iprabs)
                endif
                r(1:nr) = rato(1:nr,itabs)
                if( Tddft ) then
                  ip00 = 0
                else
                  ip00 = ip0
                endif
                ninitlt = 1
                if(Full_Potential.or.(Hubb_a .and. .not. Hubb_d)) then
                  nlma2 = nlma
                else
                  nlma2 = 1
                endif
                if( mpirank_in_mumps_group == 0 ) then
                  call tenseur_car(Base_spin,coef_g,Core_resolved,Ecinetic, &
                  Eimag(ie),Energ(ie),Enervide,Eseuil,Final_optic,Final_tddft,Full_potential,Green_i,Green_plus,Hubb_a,Hubb_d, &
                  icheck(20),ie,ip_max,ip00,is_g,lmax_probe,lmax_pot,ldip,lmoins1,loct,lplus1,lqua,lseuil,m_g,m_hubb, &
                  mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq,msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole, &
                  n_oo,nbseuil,ndim1,ndim2,nenerg_tddft,ninit1,ninitl,ninitlr,ninitlt,ninitlv,nlm_pot,nlma,nlma2,nlmamax, &
                  nr,nrm,nspin,nspino,nspinp,numat(itabs),psii,r,Relativiste,Rmtg(iprabs),Rmtsd(iprabs),rof0,rot_atom_abs, &
                  Rot_int,secdd,secdd_m,secdo,secdo_m,secdq,secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m, &
                  secqq,secqq_m,Solsing,Solsing_only,Spinorbite,Taull_abs,Tddft,V_hubb_abs,V_intmax,V0bdc,Vra,Ylm_comp)
                endif
                deallocate( r, Vra )

                if( i == 1 ) then

                  sec_atom(:) = 0._db
                  do initl = 1,ninitlr
                    do j = 1,3
                      sec_atom(initl) = sec_atom(initl) + Real(secdd(j,j,initl,mpirank), db)
                    end do
                  end do
                  sec_atom(:) = sec_atom(:) / 3

                endif

              end do

              deallocate( taull_abs )

            endif

            deallocate( Tau_ato )
            deallocate( phiato )

            if( .not. Convergence ) then
              ich = icheck(27)
            else
              ich = icheck(28)
            endif
            if( Cal_xanes ) then
              Solsing_o = Solsing_only
            else
              Solsing_o = .false.
            endif

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              tp(5) = real(time,db)
              tpt(8) = tpt(8) + tp(5) - tp(4)
            endif

            if( Density .or. .not. Cal_xanes ) &
              call cal_dens(Cal_xanes,Density_comp,drho_self,Ecinetic,Eimag(ie),Energ(ie),Enervide,Full_atom,Full_potential, &
              Hubb,Hubb_diag,iaabsi,iaprabs,iaprotoi,ich,itypei,itypepr,lla2_state,lmax_pot,lmaxat,m_hubb, &
              mpinodes,mpirank,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nlmagm,nlm_pot, &
              nrato,nrm,nrm_self,nspin,nspino,nspinp,ntype,numat,rato,Relativiste,Rmtg,Rmtsd,Solsing,Solsing_o,Spinorbite, &
              State_all_out,Statedens,Statedens_i,Taull,V_hubb,V_hubb_abs,V_intmax,V0bdc,Vrato,Ylm_comp)

            deallocate( taull )

            if( mpirank0 == 0 ) then
              call CPU_TIME(time)
              tp(6) = real(time,db)
              tpt(9) = tpt(9) + tp(6) - tp(5)
            endif

          endif ! fin de la condition sur ie <= nenerg

          if( mpinodes > 1 ) then

            if( Cal_xanes .and. .not. Optic ) call MPI_RECV_all(MPI_host_num_for_mumps,mpinodes,mpirank,mpirank0,Multipole,n_oo, &
                                                              ninitlr,secdd,secdo,secdq,secmd,secmm,secoo,secqq)
            if( Cal_xanes .and. Green_i ) call MPI_RECV_all(MPI_host_num_for_mumps,mpinodes,mpirank,mpirank0,Multipole,n_oo, &
                                                          ninitlr,secdd_m,secdo_m,secdq_m,secmd_m,secmm_m,secoo_m,secqq_m)

            if( .not. Cal_xanes .or. Density ) then
              call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
              call MPI_RECV_statedens(MPI_host_num_for_mumps,lla2_state,lmaxat,mpinodes,mpirank,mpirank0,n_atom_0, &
                                      n_atom_ind,n_atom_proto,nspin,Statedens,Statedens_i)
            endif

            if( .not. Cal_xanes ) then
              call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
              call MPI_RECV_self(drho_self,MPI_host_num_for_mumps,nlm_pot,mpinodes,mpirank,mpirank0,n_atom_0_self, &
                                 n_atom_ind_self,nrm_self,nspin)
            endif

            if( Tddft_xanes ) then
              call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
              call MPI_Bcast_tddft(imag_taull,je,mpinodes,mpirank,nbseuil,nenerg,nlmamax,nlmamax_u,nspinp,nspino,rof0)
            endif
            if( Optic_xanes ) then
              call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
              call MPI_Bcast_Optic(Taull_opt,je,mpinodes,mpirank,nenerg,nlm_probe,nspinp)
            endif

          endif

          do ie_computer = 0,mpinodes-1

            ie = ( je - 1 ) * mpinodes + ie_computer + 1

            if( ie > nenerg ) exit

            if( mpirank0 == 0 ) then

              call CPU_TIME(time)
              tp(4) = real(time,db)

              if( Cal_xanes .and. .not. Optic ) then

                call write_coabs(Allsite,angxyz,axyz,Base_spin, Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
                  Densite_atom,E_cut,Energ,Energphot,Extract,Epsii,Eseuil,Final_tddft, &
                  f_avantseuil, Full_self_abs,Green_i,Green_plus,hkl_dafs,iabsorig(multi_run),icheck(21),ie,ie_computer, &
                  Int_tens,isigpi,isymeq,jseuil,ltypcal,Moyenne,mpinodee,Multipole,n_multi_run,n_oo,natomsym,nbseuil, &
                  ncolm,ncolr,ncolt,nenerg,ninit1,ninitlr,nomabs,nomfich,nomfich_cal_convt,nomfich_s,nphi_dafs,npldafs, &
                  nphim,nplr,nplrm,nseuil,nspinp,numat_abs,nxanout,pdp,phdafs,phdf0t, &
                  phdt,pol,Polarise,poldafse,poldafss,Rot_int,sec_atom,secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, &
                  secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Self_abs,Spherical_signal, &
                  Spherical_tensor,Spinorbite,Taux_eq,v0muf,vecdafse,vecdafss,vec,Volume_maille,Xan_atom)

                if( ie == 1 ) nomfich_cal_conv(multi_run) = nomfich_cal_convt
              endif

              call CPU_TIME(time)
              tp(5) = real(time,db)
              tpt(10) = tpt(10) + tp(5) - tp(4)

              if( ( Density .and. Cal_xanes  ) .or. .not. Convergence ) then

                if( .not. Convergence ) then
                  ich = icheck(27)
                else
                  ich = icheck(28)
                endif

                call Cal_State(chg_cluster,chg_open_val,Cal_xanes,chargat_self,Density,Doping,drho_self,E_cut,E_Open_val, &
                E_Open_val_exc,E_starta,Energ,E_Fermi,Enragr,Energ_self,Fermi,Full_atom,Hubb,Hubb_diag,iaabsi,iaprotoi, &
                i_self,ich,ie,ie_computer,Int_statedens,ipr_dop,ispin_maj,itypei,itypepr,lamstdens, &
                Level_val_abs,Level_val_exc,lla_state,lla2_state,lmaxat,m_hubb,mpinodes,n_atom_0,n_atom_0_self, &
                n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nb_eq,nenerg,ngreq,nlm_pot,nomfich_s,Nonexc_g,nrato, &
                nrm,nrm_self,nspin,nspinp,ntype,numat,occ_hubb,occ_hubb_i,pop_orb_val,rato,rho_self,rho_self_t,Rmtsd,SCF_elecabs, &
                SCF_mag_fix,Self_nonexc,State_all_out,Statedens,Statedens_i,V_hubb,V_hubbard)

              endif

              call CPU_TIME(time)
              tp(6) = real(time,db)
              tpt(9) = tpt(9) + tp(6) - tp(5)

            endif

            if( .not. Cal_xanes ) then

              call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
              call MPI_Bcast(Fermi,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
              call MPI_Bcast(E_cut,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
! Si on a atteint le niveau de Fermi on sort de la boucle
              if( Fermi )  exit boucle_energ

            endif

          end do

        end do boucle_energ   ! Fin de la boucle sur l'energie.

        if( Cal_xanes .and. .not. TDDFT_xanes ) deallocate( rof0 )

        if( One_run .and. mpinodes > 1 .and. multi_run == n_multi_run .and. Cal_xanes ) deallocate( taull_stk )

        if( Save_tddft_data .and. mpirank0 == 0 .and. Cal_xanes ) call Write_tddft_data(E_cut,imag_taull,Energ_s,multi_run, &
                    n_multi_run,nbseuil,nenerg,nlmamax,nlmamax_u,nomfich,nspin,nspino,rof0)

        if( .not. Cal_xanes .and. Hubbard .and. mpinodes0 > 1 ) call MPI_Bcast_Hubb(Hubb_diag,m_hubb,mpirank0, &
                   n_atom_0_self,n_atom_ind_self,nspinp,V_hubb)

! Calcul de l'energie de l agregat

        if( .not. Cal_xanes .and. mpirank0 == 0 ) then

          call eps_coeur(ch_coeur,En_coeur,Full_atom,iaprotoi,icheck(27),itypepr,lcoeur,n_atom_0,n_atom_0_self, &
              n_atom_ind,n_atom_ind_self,n_atom_proto,natome,ncoeur,nlm_pot,nrato,nrm,nspin,ntype,numat,psi_coeur, &
              rato,Relativiste,Rmtg,Rmtsd,V_intmax,Vcato,Vxcato)

          if( i_self == 1 ) then
            En_coeur_s(:) = En_coeur(:)
            En_coeur(:) = 0._db
          else
            En_coeur(:) = En_coeur(:) - En_coeur_s(:)
          end if

          call Energ_DFT(Doping,En_cluster,Energ_self,En_coeur,excato,Full_atom,Hubb,iaprotoi,icheck(27),ipr_dop,itypepr, &
            m_hubb,n_atom_0,n_atom_0_self,n_atom_ind,n_atom_ind_self,n_atom_proto,natome,nb_eq,ngreq,nlm_pot,nrm,nrm_self, &
            nrato,nspin,nspinp,ntype,numat,rato,rho_self,rmtsd,V_hubb,Vcato,Vxcato)

          if( i_self == 1 ) then
            En_coeur(:) = 0._db
          else
            En_coeur(:) = En_coeur(:) - En_coeur_s(:)
          end if

        end if

! Echange des grandeurs liees a l'auto-coherence
        if( .not. Cal_xanes .and. mpinodes0 > 1 ) then
          m = n_atom_ind_self - n_atom_0_self + 1
          m2 = nspin * m
          n = ( 1 + nrm_self ) * nspin * m * nlm_pot
          call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(rho_self,n,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(chargat_self,m2,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(Energ_self,m,MPI_REAL8,0, MPI_COMM_WORLD,mpierr)
          call MPI_Bcast(En_cluster,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
        endif

! Preparation de l'iteration suivante pour la self-consistence
        if( .not. Cal_xanes ) then
          if( Fermi_first .and. (TDDFT .or. Optic) ) then
            Delta_E_conv = 2 * Delta_En_conv
          else
            Delta_E_conv = Delta_En_conv
          endif
          call prep_next_iter(chargat_self,chargat_self_s,Convergence,Delta_E_conv,Delta_energ,Delta_energ_s, &
               Delta_energ_t,Doping,En_cluster,En_cluster_s,En_cluster_t,Energ_self,Energ_self_s,Fermi,Fermi_first,Full_atom, &
               Hubbard,i_self,icheck(27),ipr_dop,m_hubb,mpirank0,n_atom_0_self, &
               n_atom_ind_self,n_atom_proto,n_devide,natome,natomeq,nb_eq,ngreq,nlm_pot,nrm_self,nself,nspin,nspinp, &
               p_self,p_self_s,p_self_t,p_self0,rho_self,rho_self_s,V_hubb,V_hubb_s)
        endif

        deallocate( drho_self, Statedens, Statedens_i )
        deallocate( secdd, secmd, secmm, secdq, secdo, secoo, secqq )
        deallocate( secdd_m, secmd_m, secmm_m, secdq_m, secdo_m, secoo_m, secqq_m )
        deallocate( sec_atom )
        deallocate( Energ, Eimag, Ecinetic, Ecinetic_out, Excato )
        deallocate( ibord, isbord, isrt, iato, lato )
        deallocate( lmaxa )
        deallocate( iso, lso, mato, mso )
        deallocate( nbord, nbordf )
        deallocate( nlmsa, nlmsa0, nlmso0 )
        deallocate( poidsa, poidso )
        deallocate( rho, rhons )
        deallocate( V0bdc, V0bdc_out, Vr, Vrato )
        deallocate( imoy, imoy_out, poidsov, poidsov_out )

        if( .not. ( Cal_xanes .and. ( Tddft .or. Optic ) ) ) deallocate( dV0bdcF )

        if( .not. ( ( One_run .and. Cal_xanes .and. multi_run /= n_multi_run ) .or. ( Cal_xanes .and. ( Tddft .or. Optic ) ) ) ) &
               deallocate( rsato, Vcato, Vxcato  )

        if( .not. ( ( One_run .and. Cal_xanes .and. multi_run /= n_multi_run )  ) ) &
               deallocate( rs, Vh, Vhns, Vxc)

        if( .not. Green ) then
          deallocate( gradvr )
          deallocate( Ylmso )
        endif
        deallocate( Ylmato )

        if( mpirank0 == 0 ) deallocate( Int_statedens )

      end do boucle_coh  ! fin de la boucle coherente

    end do boucle_energy_range  ! fin de la boucle gamme d energie pour l'optic

    if( .not. Recup_optic ) deallocate( Repres_comp )

    if( mpirank0 == 0 ) then
      call CPU_TIME(time)
      tp1 = real(time,db)
    endif

    if( Tddft .or. Optic ) then

      if( .not. recup_optic ) then
        Hubb_a = Hubb(itabs)
        Hubb_d = Hubb_diag(iaprabs)
        Rmtsd_abs = Rmtsd(iprabs)
      endif

      nr = nrato(itabs)
      if( Full_atom ) then
        iapr = iaabsi
      else
        iapr = iprabs
      end if
      allocate( r(nr) )
      allocate( rsato_abs(nr) )
      r(1:nr) = rato(1:nr,itabs)

      n_atom_0_t = iapr
      n_atom_ind_t = iapr
      allocate( rsato_t(0:nrm,n_atom_0_t:n_atom_ind_t) )
      allocate( Vcato_t(0:nrm,nlm_pot,n_atom_0_t:n_atom_ind_t) )
      allocate( Vxcato_t(0:nrm,nlm_pot,nspin,n_atom_0_t:n_atom_ind_t) )
      if( .not. Recup_optic ) then
        rsato_abs(1:nr) = rsato(1:nr,iapr)
        rsato_t(0:nrm,n_atom_0_t:n_atom_ind_t) = rsato(0:nrm,n_atom_0_t:n_atom_ind_t)
        Vcato_t(0:nrm,nlm_pot,n_atom_0_t:n_atom_ind_t) = Vcato(0:nrm,nlm_pot,n_atom_0_t:n_atom_ind_t)
        Vxcato_t(0:nrm,nlm_pot,nspin,n_atom_0_t:n_atom_ind_t) = Vxcato(0:nrm,nlm_pot,nspin,n_atom_0_t:n_atom_ind_t)
        deallocate( rsato )
        deallocate( Vcato )
        deallocate( Vxcato )
      endif

      Atom_nonsph_loc = .false.

      if( Tddft ) then

        allocate( rhoa(nr,nspin) )
        rhoa(1:nr,:) = rhoato_abs(1:nr,:)

        call main_tddft(alfpot,angxyz,Allsite,Atom_nonsph_loc,Atomic_scr,axyz,Base_spin,BSE,coef_g, &
        Cartesian_tensor,Core_resolved,Dafs,Dafs_bio,Delta_edge,Delta_Eseuil,Densite_atom,Dipmag, &
        dv_ex_nex,dV0bdcF,Dyn_eg,Dyn_g,E_cut,E_cut_imp,E_Fermi_man,Ecent,Elarg, &
        Energ_s,Energphot,Epsii,Extract,Epsii_moy,Eseuil,Estart,f_avantseuil,Full_atom,Full_potential,Full_self_abs, &
        Gamma_hole,Gamma_hole_imp,Gamma_max,Gamma_tddft,Green_int,Green_plus,hkl_dafs,Hubb_a,Hubb_d,icheck,iaabsi, &
        iabsorig(multi_run),iapot,iaprotoi,imag_taull,imparite,iprabs,iprabs_nonexc,is_g,isigpi,isymeq, &
        itabs,itypepr,jseuil,Kern_fac,ldip,lmax_pot,lmaxabs_t,lmaxat0,lmaxfree,lmoins1,loct,lplus1, &
        lqua,lseuil,ltypcal,m_g,m_hubb,Magnetic,Moyenne,MPI_host_num_for_mumps,mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq, &
        msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,multi_run,Multipole,n_atom_0_t, &
        n_atom_ind_t,n_atom_proto,n_multi_run,n_oo,n_tens_max,natome, &
        natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ngamh,ninit1,ninitl,ninitlr,nlm_pot,nlmamax,nlmamax_u, &
        nomabs,nomfich,nomfich_cal_tddft_conv(multi_run),nomfich_s,nomfich_tddft_data,Nonexc_g, &
        nphi_dafs,nphim,npldafs,nplr,nplrm,npsom,nptmoy,nptmoy_out,nr,nrato,nrm,nseuil,nspin,nspino,nspinp, &
        ntype,numat_abs,nxanout,Octupole,pdp,phdafs,phdf0t,phdt,pol,Polarise,poldafse,poldafss, &
        psii,Quadrupole,r,rato,Recup_tddft_data,Relativiste,rhoa,Rmtg,Rmtsd_abs, &
        rof0,rot_atom_abs,Rot_int,RPALF,rsato_t,rsato_abs,rsbdc,rsbdc_out,Self_abs,Solsing_only, &
        Spherical_signal,Spherical_tensor,Spinorbite,Taux_eq,Tddft_so,tpt,V_intmax,V_hubb_abs,V0muf, &
        Vcato_t,Vec,vecdafse,vecdafss,Vhbdc,Vhbdc_out,Volume_maille,VxcbdcF,VxcbdcF_out, &
        Vxcato_t,Workf,xsect_file,xyz,Ylm_comp)

        deallocate( imag_taull, rhoa, rof0 )

      else

        if( Recup_optic ) then
          allocate( dv0bdcF(nspin) )
          allocate( Rmtg(0:n_atom_proto) )
          allocate( Taull_opt(nenerg_s,nlm_probe,nlm_probe,nspinp,nspinp) )
          allocate( iaprotoi(natome) )
          npsom = 0
          allocate( xyz(4,npsom) )
          call Read_optic_data(dv0bdcF, Energ_s, Hubb_a, Hubb_d, iaprotoi, lmax_probe, lmaxabs_t, m_hubb, mpinodes, mpirank0, &
               multi_run, n_atom_0_t, n_atom_ind_t, n_atom_proto, n_multi_run, natome, nenerg_s, nlm_pot, nlm_probe, &
               nomfich_optic_data, nptmoy, nptmoy_out, nrm, nspin, nspinp, Rmtg, Rmtsd_abs, rsato_t, rsbdc, rsbdc_out, Taull_opt, &
               V_hubb_abs, V0muf, Vcato_t, Vhbdc, Vhbdc_out, Vxcato_t, VxcbdcF, VxcbdcF_out)
! Energ_s lu au dessus est la vraie energie (ce qui n'est pas le cas sans recup_optic).
! On enleve E_cut en dessous qui est remis en entree de main optic.
! On conserve ainsi le E_cut_imp seulement pour le point de coupure et pour aucun décalage.
          if( E_Fermi_man ) then
            Energ_s(:) = Energ_s(:) - E_cut_imp
          else
            Energ_s(:) = Energ_s(:) - E_cut
          endif

        elseif( Save_optic .and. mpirank0 == 0 ) then
          call Write_optic_data(dv0bdcF, E_cut, Energ_s, Full_atom, Hubb_a, Hubb_d, iaabsi, iaprotoi, iprabs, lmax_probe, &
               lmaxabs_t, m_hubb, multi_run, n_atom_0_t, n_atom_ind_t, n_atom_proto, n_multi_run, natome, nenerg_s, nlm_pot, &
               nlm_probe, nomfich, nptmoy, nptmoy_out, nrm, nspin, nspinp, Rmtg, Rmtsd_abs, rsato_t, rsbdc, rsbdc_out, Taull_opt, &
               V_hubb_abs, V0muf, Vcato_t, Vhbdc, Vhbdc_out, Vxcato_t, VxcbdcF, VxcbdcF_out)
        endif

        call main_optic(alfpot,angxyz,Allsite,Atom_nonsph_loc,axyz,Base_spin,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
          Densite_atom,dV0bdcF,E_cut,E_cut_imp,E_Fermi_man,Eclie,Eneg,Energ_s, &
          Extract,Eseuil,Full_atom,Full_potential,Full_self_abs,Green_int,Green_plus,hkl_dafs,Hubb_a,Hubb_d,icheck,iaabsi, &
          iabsorig(multi_run),iapot,iaprotoi,ip_max,ip0,iprabs,iprabs_nonexc,isigpi,isymeq, &
          itabs,itypepr,jseuil,ldip,lmax_pot,lmax_probe,lmaxabs_t,lmaxat0,lmaxfree,lmoins1,loct,lplus1,lqua, &
          lseuil,ltypcal,m_hubb,Magnetic,Moyenne,MPI_host_num_for_mumps,mpinodes,mpirank,mpirank0,msymdd,msymddi,msymdq, &
          msymdqi,msymdo,msymdoi,msymoo,msymooi,msymqq,msymqqi,Multipole,n_atom_0_t,n_atom_ind_t,n_atom_proto,n_multi_run, &
          n_oo,n_tens_max,natome,natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ninit1,ninitl,ninitlr,nlm_pot, &
          nlm_probe,nomabs,nomfich,nomfich_s,Nonexc_g,nphi_dafs,nphim,npldafs,nplr,nplrm, &
          npsom,nptmoy,nptmoy_out,nr,nrato,nrm,nseuil,nspin,nspino,nspinp,ntype, &
          numat_abs,nxanout,pdp,phdafs,phdf0t,phdt,pol,Polarise,poldafse,poldafss,psii, &
          r,rato,Relativiste,Rmtg,Rmtsd_abs,rot_atom_abs,Rot_int,rsato_t,rsbdc,rsbdc_out,Self_abs,Solsing_only, &
          Spherical_signal,Spherical_tensor,Spinorbite,Taull_opt,Taux_eq,tpt,V_intmax,V_hubb_abs,V0muf, &
          Vcato_t,vec,vecdafse,vecdafss,Vhbdc,Vhbdc_out,Volume_maille,VxcbdcF,VxcbdcF_out,Vxcato_t,Workf,xyz,Ylm_comp)

        deallocate( Taull_opt )
      endif

      deallocate( r, rsato_abs )
      deallocate( dV0bdcF )
      deallocate( rsato_t, Vcato_t, Vxcato_t );

    endif

    if( mpirank0 == 0 ) then
      call CPU_TIME(time)
      tp2 = real(time,db)
      tpt(12) = tp2 - tp1
    endif

    if( .not. ( Second_run .or. Recup_optic ) ) deallocate( Chg_cluster, E_starta, ispin_maj )

    if( Self_cons .and. nself > 0 .and. .not. ( Second_run .or. Recup_optic ) ) then
      deallocate( ch_coeur )
      deallocate( Energ_coh, Eimag_coh )
    end if

    if( .not. Recup_optic ) then
      deallocate( Atom_axe, chargat_init )
      deallocate( chargat_self, chargat_self_s )
      deallocate( distai, drho_ex_nex, dvc_ex_nex, dv_ex_nex )
      deallocate( Energ_self, Energ_self_s )
      deallocate( En_coeur, En_coeur_s )
      deallocate( ia_eq, ia_eq_inv, ia_eq_inv_self, ia_rep )
      deallocate( igroupi, iopsym_atom, is_eq, itypei )
      deallocate( nb_eq, nb_rpr, nb_rep_t )
      deallocate( pop_orb_val, posi, posi_self )
      deallocate( rho_chg, rho_self, rho_self_s, rho_self_t )
      deallocate( rhoato_init, rot_atom, V_abs_i, Vcato_init )
      deallocate( Hubb_diag, occ_hubb, occ_hubb_i, V_hubb, V_hubb_s )
    endif

    deallocate( iaprotoi )
    deallocate( V_hubb_abs, VxcbdcF, VxcbdcF_out )
    deallocate( rhoato_abs )

    if( .not. One_run .or. multi_run == n_multi_run ) then
      deallocate( xyz )
      if( .not. Recup_optic ) then
        deallocate( cgrad, clapl, ivois, isvois )
        deallocate( lmaxat, numia, rmtg, rmtg0, rmtsd, rvol )
      endif
    endif

    deallocate( karact, Axe_atom_clu, dista, iapot )
    deallocate( iaproto, igroup, isymeq, itypep, itypepr )
    deallocate( pos, poseq, Taux_eq )
    deallocate( poldafse, poldafss, phdafs, phdf0t, phdt )
    deallocate( vecdafse, vecdafss, ltypcal, nomabs, pdp, pol, vec )

    if( mpirank0 == 0 ) deallocate( Int_tens )

    if( icheck(1) > 0 .and. mpirank0 == 0 ) then
      call CPU_TIME(time)
      tp1 = real(time,db)
      tpt_tot = tp1 - tp_init + tpt(1)
      tptt = sum( tpt(1:10) )
      do i = 1,10000
        if( tpt_tot >= tptt - 0.1_db ) exit
        tpt_tot = tpt_tot + 86400.
      end do
      itph = int( tpt_tot / 3600 )
      rtph = tpt_tot - itph * 3600
      itpm = int( rtph / 60 )
      itps = nint( rtph - itpm * 60 )
      write(3,210) multi_run
      write(3,220) (nomspr(i), tpt(i), i = 1,10)
      if( Optic ) then
        write(3,220) 'Rempli', tpt1, 'Triang', tpt2, 'SCF   ', tpt(11), 'Optic', tpt(12)
      elseif( TDDFT ) then
        write(3,220) 'Rempli', tpt1, 'Triang', tpt2, 'SCF   ', tpt(11), 'Tddft', tpt(12)
      else
        write(3,220) 'Rempli', tpt1, 'Triang', tpt2, 'SCF   ', tpt(11)
      endif
      if( itph > 0 .or. itpm > 0 ) then
        write(3,240) 'Total ', tpt_tot, itph, itpm, itps
      else
        write(3,220) 'Total ', tpt_tot
      endif
    endif

  end do boucle_multi ! Fin boucle sur sites non equivalents

  if( Recup_tddft_data ) Convolution_cal = .false.

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
      do i = 1,n_multi_run
        if( Skip_run(i) ) then
          is = 0
        else
          is = 1
        endif
        if( Run_done(i) ) then
          ir = 0
        else
          ir = 1
        endif
        write(itape1,'(2i2)') ir, is
      end do
    else
      write(itape1,'(A)') ' calculation'
      do multi_run = 1,n_multi_run
        if( Run_done(multi_run) ) cycle
        write(itape1,'(A)') nomfich_cal_conv(multi_run)
      end do
      if( Tddft .and. .not. Extract .and. ( .not. Gamma_tddft .or. Dafs  .or. n_multi_run > 1 ) ) then
      write(itape1,'(A)') ' cal_tddft'
      do multi_run = 1,n_multi_run
        if( Run_done(multi_run) ) cycle
        write(itape1,'(A)') nomfich_cal_tddft_conv(multi_run)
      end do
      endif
    endif
  endif

! desallocation des tableaux attribues avant la boucle multi_run

  deallocate( Atom_with_axe)
  deallocate( Chargat, Epsii, Eseuil )
  deallocate( coef_g, m_g, is_g )
  deallocate( igreq, isymqa )
  deallocate( lcoeur )
  deallocate( msymoo, msymooi )
  deallocate( ncoeur, ngreq, ngreqm, nomfich_cal_conv )
  deallocate( popatm, popexc, posq, psi_coeur, psii, psival, psi_open_val )
  deallocate( rato, rhoit, rho_coeur )
  deallocate( Run_done, Skip_run )
  if( .not. Extract ) deallocate( Energ_s, Eimag_s )
  if( Tddft ) deallocate( nomfich_cal_tddft_conv )

! Desallocation des tableaux alloues avant le sous-programme lecture
  deallocate( Angle_or, Angpoldafs, Axe_atom_gr, Axe_atom_grn )
  deallocate( com, cdil )
  deallocate( Ecrantage )
  deallocate( Hubb, hkl_dafs, hybrid )
  deallocate( iabsm, iabsorig, icom, isigpi, itdil, itype, its_lapw )
  deallocate( Kgroup )
  deallocate( ldil, lvval )
  deallocate( nlat, nomfile_atom, norbv, nphi_dafs, nposextract, nrato, nrato_lapw, nsymextract, numat, nvval )
  deallocate( occ_hubb_e )
  deallocate( pdpolar, polar, poldafsem, poldafssm, pop_nonsph, popatc, popats, popatv, popval, posn )
  deallocate( r0_lapw, rchimp, rlapw, rmt, rmtimp, Rot_Atom_gr, rotloc_lapw )
  deallocate( Taux_oc, Temp_coef )
  deallocate( V_hubbard, V0bdcFimp, Vecdafsem, Vecdafssm, Veconde )
  deallocate( Wien_matsym, Wien_taulap )

  return
  125 format(//' One_run process is not possible with absorbing atom in excited state !',// &
               ' Excited absorbing atom is the default for K, L1, M1, N1 edges.', // &
               ' Use non excited absorbing atom (keyword nonexc)',/ &
               ' or keep a multi run process.'//)
  130 format(/' Number of calculated non equivalent absorbing atom =', i3)
  140 format(/1x,120('-')//,' Last cycle, XANES calculation')
  150 format(/1x,120('-')//,' Begining of cycle',i3)
  160 format(/' E_cut =',f11.5,' eV')
  170 format(///' E_kinetic =',f7.3,' eV < 0.',/ ' Start the calculation at higher energy !'///)
  180 format(///' or E_kinetic_ext =',f7.3,' eV < 0.',/ ' Start the calculation at higher energy !'///)
  190 format(/' ---- Sphere -------',100('-'))
  210 format(/,1x,120('-')//' Subroutine times for absorbing atom',i3,' (sCPU)')
  220 format(4(4x,a6,' =',f12.2))
  240 format(4x,a6,' =',f12.2,' =',i4,' h,',i3,' min,',i3,' sCPU')
end

!***********************************************************************

subroutine extract_write_coabs(Allsite,Ang_rotsup,angxyz,axyz,Base_spin,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
         Densite_atom,E_cut,Energ_s,Energphot,Epsii,Eseuil,Final_tddft, &
         f_avantseuil,Full_self_abs,Green_int,Green_plus,hkl_dafs,iabsorig,icheck,Int_tens,isigpi,isymeq, &
         isymext,jseuil,ltypcal,Moyenne,mpinodee,multi_run,Multipole,n_multi_run,n_oo,n_tens_max, &
         natomsym,nbseuil,ncolm,ncolr,ncolt,nenerg_s,ninit1,ninitlr,nomabs,nomfich,nomfich_cal_conv, &
         nomfich_cal_convt,nom_fich_extract,nomfich_s,nphi_dafs,nphim,npldafs,nplr,nplrm,nseuil,nspin,numat_abs,nxanout,pdp, &
         phdafs,phdf0t,phdt,pol,Polarise,poldafse,poldafss,Rot_int,Self_abs,Spherical_signal,Spherical_tensor,Spinorbite, &
         Taux_eq,Tddft,V0muf,vecdafse,vecdafss,Vec,Volume_maille,Xan_atom)

  use declarations
  implicit none

  integer:: iabsorig, icheck, ie, ie_computer, isymext, jseuil, &
    mpinodee, multi_run, n_multi_run, n_oo, n_tens_max, natomsym, nbseuil, ncolm, &
    ncolr, ncolt, ninit1, ninitlr, nenerg_s, nseuil, nspin, nphim, npldafs, nplr, nplrm, numat_abs, nxanout
  integer, dimension(3,npldafs):: hkl_dafs
  integer, dimension(npldafs,2):: isigpi
  integer, dimension(npldafs):: nphi_dafs
  integer, dimension(natomsym):: isymeq

  character(len=132):: nomfich, nomfich_cal_conv, nomfich_cal_convt, nom_fich_extract, nomfich_s
  character(len=13), dimension(nplrm):: ltypcal
  character(len=Length_word), dimension(ncolm):: nomabs

  complex(kind=db):: f_avantseuil
  complex(kind=db), dimension(npldafs,nphim):: phdf0t, phdt
  complex(kind=db), dimension(natomsym,npldafs):: phdafs
  complex(kind=db), dimension(3,nplrm):: pol
  complex(kind=db), dimension(3,npldafs,nphim):: poldafse, poldafss
  complex(kind=db), dimension(3,3,ninitlr,0:0):: secdd, secdd_m, secmd, secmd_m, secmm, secmm_m
  complex(kind=db), dimension(3,3,3,ninitlr,0:0):: secdq, secdq_m
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:0):: secdo, secdo_m, secqq, secqq_m
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:0):: secoo, secoo_m

  logical:: Allsite, Base_spin, Cartesian_tensor, Core_resolved, Dafs, Dafs_bio, &
     Energphot, Final_tddft, Full_self_abs, Green_int, Green_plus, Moyenne, Polarise, Self_abs, Spherical_signal, &
     Spherical_tensor, Spinorbite, Tddft, Tensor_rot, Xan_atom

  logical, dimension(10):: Multipole

  real(kind=db):: Densite_atom, E_cut, V0muf, Volume_maille
  real(kind=db), dimension(3):: Ang_rotsup, angxyz, axyz
  real(kind=db), dimension(nbseuil):: Eseuil
  real(kind=db), dimension(ninitlr):: Epsii, Sec_atom
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(n_tens_max*ninitlr,0:natomsym):: Int_tens
  real(kind=db), dimension(3,3):: Rot_int, Rotsup
  real(kind=db), dimension(nplrm,2):: pdp
  real(kind=db), dimension(natomsym):: Taux_eq
  real(kind=db), dimension(3,nplrm):: vec
  real(kind=db), dimension(3,npldafs,nphim):: vecdafse, vecdafss

! Sec_atom n'est pas extrait du fichier bav.
  Sec_atom(:) = 0._db

  ie_computer = 0
  do ie = 1,nenerg_s
    call extract_coabs(Ang_rotsup,Core_resolved,Green_int,icheck,ie,isymext,multi_run,Multipole, &
            n_oo,nenerg_s,ninit1,ninitlr,nom_fich_extract,Rotsup,secdd,secdd_m,secdo,secdo_m,secdq, &
            secdq_m,secmd,secmd_m,secmm,secmm_m,secoo,secoo_m,secqq,secqq_m,Tensor_rot,Tddft)
    call write_coabs(Allsite,angxyz,axyz,Base_spin,Cartesian_tensor,Core_resolved,Dafs,Dafs_bio, &
            Densite_atom,E_cut,Energ_s, &
            Energphot,.true.,Epsii,Eseuil,Final_tddft,f_avantseuil,Full_self_abs,Green_int,Green_plus,hkl_dafs, &
            iabsorig,icheck,ie,ie_computer,Int_tens,isigpi,isymeq,jseuil,ltypcal, &
            Moyenne,mpinodee,Multipole,n_multi_run,n_oo,natomsym,nbseuil, ncolm,ncolr,ncolt,nenerg_s,ninit1,ninitlr,nomabs, &
            nomfich,nomfich_cal_convt,nomfich_s,nphi_dafs,npldafs, nphim,nplr, &
            nplrm,nseuil,nspin,numat_abs,nxanout,pdp,phdafs,phdf0t, phdt,pol,Polarise,poldafse,poldafss,Rot_int,sec_atom, &
            secdd,secdd_m,secdq,secdq_m,secdo,secdo_m, secmd,secmd_m,secmm,secmm_m,secoo,secoo_m, &
            secqq,secqq_m,Self_abs,Spherical_signal, Spherical_tensor,Spinorbite,Taux_eq,V0muf, &
            Vecdafse,Vecdafss,Vec,Volume_maille,Xan_atom)
    if( ie == 1 ) nomfich_cal_conv = nomfich_cal_convt
  end do

  return
end

!***********************************************************************

subroutine MPI_RECV_all(MPI_host_num_for_mumps,mpinodes,mpirank,mpirank0,Multipole,n_oo,ninitlr,secdd,secdo, &
                  secdq,secmd,secmm,secoo,secqq)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: idim1, idim2, idim3, idim4, MPI_host_num_for_mumps, mpierr, mpinodes, mpirank, mpirank0, n_oo, ninitlr, rang

  complex(kind=db), dimension(3,3,ninitlr,0:mpinodes-1):: secdd, secmd, secmm
  complex(kind=db), dimension(3,3,3,ninitlr,0:mpinodes-1):: secdq
  complex(kind=db), dimension(3,3,3,3,ninitlr,0:mpinodes-1):: secdo, secqq
  complex(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:mpinodes-1):: secoo

  logical:: E1E1, E1E2, E1E3, E1M1, E2E2, E3E3, M1M1

  logical, dimension(10):: Multipole

  real(kind=db), dimension(3,3,ninitlr,0:mpinodes-1):: secdd_er, secmd_er, secmm_er, secdd_ei, secmd_ei, secmm_ei
  real(kind=db), dimension(3,3,3,ninitlr,0:mpinodes-1):: secdq_er, secdq_ei
  real(kind=db), dimension(3,3,3,3,ninitlr,0:mpinodes-1):: secdo_ei, secdo_er, secqq_ei, secqq_er
  real(kind=db), dimension(3,n_oo,3,n_oo,ninitlr,0:mpinodes-1):: secoo_ei, secoo_er

  E1E1 = Multipole(1); E1E2 = Multipole(2); E1E3 = Multipole(3);
  E1M1 = Multipole(4); E2E2 = Multipole(6);
  E3E3 = Multipole(7); M1M1 = Multipole(8)

  idim1 = 9 * ninitlr
  idim2 = 3 * idim1
  idim3 = 3 * idim2
  idim4 = ( n_oo**2 ) * idim2

! Recopies afin de stoquer les parties reele et imaginaire; l'indice sur
! les processeurs est en dernier

  if( E1E1 )  then
    secdd_er(:,:,:,mpirank) = real( secdd(:,:,:,mpirank),db )
    secdd_ei(:,:,:,mpirank) = aimag( secdd(:,:,:,mpirank) )
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

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  if( mpirank0 == 0 ) then
    if( E1E1 )  then
      call MPI_GATHER(MPI_IN_PLACE,idim1,MPI_REAL8,secdd_er,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      call MPI_GATHER(MPI_IN_PLACE,idim1,MPI_REAL8,secdd_ei,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
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

  elseif( mod(mpirank0,MPI_host_num_for_mumps) == 0 ) then

	  if( E1E1 ) then
	    call MPI_GATHER(secdd_er(1,1,1,mpirank),idim1,MPI_REAL8,secdd_er,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
	    call MPI_GATHER(secdd_ei(1,1,1,mpirank),idim1,MPI_REAL8,secdd_ei,idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
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

  end if

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

! On reconstruit les secXX; maintenant le processus central voit les resultats de tous les autres
  if( mpirank0 == 0 ) then
    do rang = 1,mpinodes-1
      if( E1E1 ) secdd(:,:,:,rang) = cmplx( secdd_er(:,:,:,rang), secdd_ei(:,:,:,rang),db )
      if( E1M1 ) secmd(:,:,:,rang) = cmplx( secmd_er(:,:,:,rang), secmd_ei(:,:,:,rang),db )
      if( M1M1 ) secmm(:,:,:,rang) = cmplx( secmm_er(:,:,:,rang), secmm_ei(:,:,:,rang),db )
      if( E1E2 ) secdq(:,:,:,:,rang) = cmplx( secdq_er(:,:,:,:,rang), secdq_ei(:,:,:,:,rang),db )
      if( E2E2 ) secqq(:,:,:,:,:,rang) = cmplx( secqq_er(:,:,:,:,:,rang), secqq_ei(:,:,:,:,:,rang),db )
      if( E1E3 ) secdo(:,:,:,:,:,rang) = cmplx( secdo_er(:,:,:,:,:,rang), secdo_ei(:,:,:,:,:,rang),db )
      if( E3E3 ) secdo(:,:,:,:,:,rang) = cmplx( secoo_er(:,:,:,:,:,rang), secoo_ei(:,:,:,:,:,rang),db )
    end do
  end if

  return
end

!***********************************************************************

subroutine MPI_RECV_statedens(MPI_host_num_for_mumps,lla2_state,lmaxat,mpinodes,mpirank,mpirank0,n_atom_0, &
                              n_atom_ind,n_atom_proto,nspin,Statedens,Statedens_i)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: MPI_host_num_for_mump
  integer, dimension(0:n_atom_proto):: lmaxat
  integer, dimension(0:n_atom_proto,0:mpinodes-1):: lmaxat_e

  real(kind=db), dimension(lla2_state,nspin,lla2_state,nspin,n_atom_0:n_atom_ind,0:mpinodes-1):: Statedens, Statedens_i

  lmaxat_e(:,mpirank) = lmaxat(:)
  idim = n_atom_proto + 1

  if( mpirank0 == 0 ) then
    call MPI_GATHER(MPI_IN_PLACE,idim,MPI_INTEGER,lmaxat_e,idim,MPI_INTEGER,0,MPI_COMM_GATHER,mpierr)
  elseif ( mod(mpirank0,MPI_host_num_for_mumps) == 0 ) then
    call MPI_GATHER(lmaxat_e(0,mpirank),idim,MPI_INTEGER,lmaxat_e,idim,MPI_INTEGER,0,MPI_COMM_GATHER,mpierr)
  end if

! Ceux sont les lmax correspondant a l'energie la plus grande qui sont
! les plus grands.

  if( mpirank0 == 0 ) lmaxat(:) = lmaxat_e(:,mpinodes-1)

! MPI_GATHER: le choix lorsque tous les ordinateurs envoyent le meme nombre
! d'elements a l'ordinateur central

  idim = ( n_atom_ind - n_atom_0 + 1 ) * (lla2_state * nspin)**2

  if( mpirank0 == 0 ) then
    call MPI_GATHER(MPI_IN_PLACE,idim,MPI_REAL8,Statedens,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  elseif( mod(mpirank0,MPI_host_num_for_mumps) == 0 ) then
    call MPI_GATHER(Statedens(1,1,1,1,n_atom_0,mpirank),idim,MPI_REAL8,Statedens,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  end if

  if( mpirank0 == 0 ) then
    call MPI_GATHER(MPI_IN_PLACE,idim,MPI_REAL8,Statedens_i,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
  elseif( mod(mpirank0,MPI_host_num_for_mumps) == 0 ) then
    call MPI_GATHER(Statedens_i(1,1,1,1,n_atom_0,mpirank),idim, MPI_REAL8,Statedens_i,idim,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
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

subroutine MPI_RECV_self(drho_self,MPI_host_num_for_mumps,nlm_pot,mpinodes,mpirank,mpirank0,n_atom_0_self,n_atom_ind_self, &
                         nrm_self,nspin)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: idim1, lm, MPI_host_num_for_mumps, mpinodes, mpirank, mpirank0, n_atom_0_self, n_atom_ind_self, nlm_pot, nrm_self, &
            nspin, rang

  real(kind=db), dimension(nrm_self,nlm_pot,nspin,n_atom_0_self:n_atom_ind_self,0:mpinodes-1):: drho_self
  real(kind=db), dimension(nrm_self,nspin,0:mpinodes-1):: drho_self_e

  idim1 = nrm_self * nspin

  do lm = 1,nlm_pot
    do ia = n_atom_0_self,n_atom_ind_self

      drho_self_e(:,:,mpirank) = drho_self(:,lm,:,ia,mpirank)

! Ici la barriere est tres importante: si on l'ommet on risque d'utiliser le
!    meme buffer rhov_self_eX par deux processus differents aux iapr differents

      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

      if( mpirank0 == 0 ) then
        call MPI_GATHER(MPI_IN_PLACE,idim1,MPI_REAL8, drho_self_e,idim1,MPI_REAL8, 0,MPI_COMM_GATHER,mpierr)
      elseif( mod(mpirank0,MPI_host_num_for_mumps) == 0 ) then
        call MPI_GATHER(drho_self_e(1,1,mpirank), idim1,MPI_REAL8,drho_self_e, idim1,MPI_REAL8,0,MPI_COMM_GATHER,mpierr)
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

      if( mpirank0 == 0 ) then
        do rang = 1,mpinodes-1
          drho_self(:,lm,:,ia,rang) = drho_self_e(:,:,rang)
        end do
      end if

! Cette barriere est importante, car drho_self_e est utilise en tant que
! buffer et ne devrait pas etre remplie pour deux ia differents
      call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

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

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  if( mpirank /= 0 ) V_hubb(:,:,:,:,:) = cmplx( V_hubb_r(:,:,:,:,:), V_hubb_i(:,:,:,:,:), db )

  ndim = n_atom_ind_self - n_atom_0_self + 1

  call MPI_Bcast(Hubb_diag,ndim,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)

  return
end

!***********************************************************************

subroutine MPI_Bcast_tddft(imag_taull,je,mpinodes,mpirank,nbseuil,nenerg,nlmamax,nlmamax_u,nspin,nspino,rof0)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: ie, ie_computer, je, mpinodes, mpirank, ndim, nbseuil, nenerg, nlmamax, nlmamax_u, nspin, nspino

  complex(kind=db), dimension(nenerg,nlmamax,nspin,nspino,nbseuil):: rof0

  real(kind=db), dimension(nenerg,nlmamax_u,nspin,nlmamax_u,nspin):: imag_taull
  real(kind=db), dimension(nlmamax_u,nspin,nlmamax_u,nspin):: imag_taull_e
  real(kind=db), dimension(nlmamax,nspin,nspino,nbseuil):: rof0_i, rof0_r

  ndim = ( nlmamax_u * nspin )**2

  do ie_computer = 0,mpinodes-1

    ie = ( je - 1 ) * mpinodes + ie_computer + 1
    if( ie > nenerg ) exit

    if( ie_computer == mpirank ) imag_taull_e(:,:,:,:) = imag_taull(ie,:,:,:,:)

    call MPI_Bcast(imag_taull_e,ndim,MPI_REAL8,ie_computer,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if( ie_computer /= mpirank ) imag_taull(ie,:,:,:,:) = imag_taull_e(:,:,:,:)

  end do

  ndim = nlmamax * nspin * nspino * nbseuil

  do ie_computer = 0,mpinodes-1

    ie = ( je - 1 ) * mpinodes + ie_computer + 1
    if( ie > nenerg ) exit

    if( ie_computer == mpirank ) then
      rof0_r(:,:,:,:) = real( rof0(ie,:,:,:,:),db )
      rof0_i(:,:,:,:) = aimag( rof0(ie,:,:,:,:) )
    endif

    call MPI_Bcast(rof0_r,ndim,MPI_REAL8,ie_computer,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(rof0_i,ndim,MPI_REAL8,ie_computer,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if( ie_computer /= mpirank ) rof0(ie,:,:,:,:) = cmplx( rof0_r(:,:,:,:), rof0_i(:,:,:,:),db )

  end do

  return
end

!***********************************************************************

subroutine MPI_Bcast_Optic(Taull_opt,je,mpinodes,mpirank, nenerg,nlm_probe,nspin)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  integer:: ie, ie_computer, je, mpinodes, mpirank, ndim, nlm_probe, nenerg, nspin

  complex(kind=db), dimension(nenerg,nlm_probe,nlm_probe,nspin, nspin):: Taull_opt

  real(kind=db), dimension(nlm_probe,nlm_probe,nspin,nspin):: Taull_i, Taull_r

  ndim = ( nlm_probe * nspin )**2

  do ie_computer = 0,mpinodes-1

    ie = ( je - 1 ) * mpinodes + ie_computer + 1
    if( ie > nenerg ) exit

    if( ie_computer == mpirank ) then
      Taull_r(:,:,:,:) = real( Taull_opt(ie,:,:,:,:),db )
      Taull_i(:,:,:,:) = aimag( Taull_opt(ie,:,:,:,:) )
    endif

    call MPI_Bcast(Taull_r,ndim,MPI_REAL8,ie_computer,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Taull_i,ndim,MPI_REAL8,ie_computer,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if( ie_computer /= mpirank ) Taull_opt(ie,:,:,:,:) = cmplx( Taull_r(:,:,:,:), Taull_i(:,:,:,:),db )

  end do

  return
end

!***********************************************************************

subroutine Write_tddft_data(E_cut,imag_taull,Energ,multi_run,n_multi_run,nbseuil,nenerg,nlmamax, &
                    nlmamax_u,nomfich,nspin,nspino,rof0)

  use declarations
  implicit none
  include 'mpif.h'

  character(len=132):: nomfich, nomfich_tddft_data

  integer:: ie, long, nbseuil,multi_run, n_multi_run, nenerg, nlmamax, nlmamax_u, nspin, nspino

  complex(kind=db), dimension(nenerg,nlmamax,nspin,nspino, nbseuil):: rof0

  real(kind=db):: E_cut
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(nenerg,nlmamax_u,nspin,nlmamax_u,nspin):: imag_taull

  nomfich_tddft_data = nomfich

  long = len_trim(nomfich_tddft_data)
  nomfich_tddft_data(long+1:long+15) = '_tddft_data.txt'

  if( multi_run == 1 ) open(20, file = nomfich_tddft_data )

  write(20,*) E_cut*Rydb, nenerg, Energ(nenerg)*Rydb

  do ie = 1,nenerg
    write(20,*) Energ(ie)*Rydb
    write(20,*) imag_taull(ie,:,:,:,:)
  end do
  do ie = 1,nenerg
    write(20,*) Energ(ie)*Rydb
    write(20,*) rof0(ie,:,:,:,:)
  end do

  if( multi_run == n_multi_run ) close(20)

  return
end

!***********************************************************************

subroutine Write_optic_data( dv0bdcF, E_cut, Energ_s, Full_atom, Hubb_a, Hubb_d, iaabsi, iaprotoi, iprabs, lmax_probe, &
               lmaxabs_t, m_hubb, multi_run, n_atom_0, n_atom_ind, n_atom_proto, n_multi_run, natome, nenerg_s, nlm_pot, &
               nlm_probe, nomfich, nptmoy, nptmoy_out, nrm, nspin, nspinp, Rmtg, Rmtsd, rsato, rsbdc, rsbdc_out, Taull_opt, &
               V_hubb, V0muf, Vcato, Vhbdc, Vhbdc_out, Vxcato, VxcbdcF, VxcbdcF_out)

  use declarations
  implicit none
  include 'mpif.h'

  character(len=132):: nomfich, nomfich_optic_data

  integer:: iaabsi, ie, iprabs, lmax_probe, lmaxabs_t, long, m_hubb, multi_run, n_atom_0, n_atom_ind, n_atom_proto, n_multi_run, &
            natome, nenerg_s, nlm_pot, nlm_probe, nptmoy, nptmoy_out, nrm, nspin, nspinp

  integer, dimension(natome):: iaprotoi

  complex(kind=db), dimension(nenerg_s,nlm_probe,nlm_probe,nspinp,nspinp):: Taull_opt
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  logical:: Full_atom, Hubb_a, Hubb_d

  real(kind=db):: E_cut, Rmtsd, rsbdc, rsbdc_out, Vhbdc, Vhbdc_out, V0muf
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(0:n_atom_proto):: Rmtg
  real(kind=db), dimension(nspin):: dv0bdcF, VxcbdcF, VxcbdcF_out
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind) :: Vxcato
  real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: rsato

  nomfich_optic_data = nomfich

  long = len_trim(nomfich_optic_data)
  nomfich_optic_data(long+1:long+15) = '_optic_data.txt'

  if( multi_run == 1 ) open(20, file = nomfich_optic_data )

  write(20,*) E_cut, Full_atom
  write(20,*) iaabsi, iprabs, natome, nenerg_s, nlm_pot, nlm_probe

  write(20,*) lmax_probe, lmaxabs_t, nptmoy, nptmoy_out
  write(20,'(10i5)') iaprotoi(:)

  do ie = 1,nenerg_s
    write(20,*) ( Energ_s(ie) + E_cut )*Rydb
    write(20,*) Taull_opt(ie,:,:,:,:)
  end do

  write(20,*) Hubb_a, Hubb_d
  write(20,*) Rmtg(:), Rmtsd
  write(20,*) dv0bdcF(:), rsbdc, rsbdc_out, Vhbdc, Vhbdc_out, VxcbdcF(:), VxcbdcF_out(:), V0muf
  write(20,*) rsato(0:nrm,:)
  write(20,*) Vcato(0:nrm,:,:)
  write(20,*) Vxcato(0:nrm,:,:,:)
  if( Hubb_a ) write(20,*) V_hubb(:,:,:,:)

  if( multi_run == n_multi_run ) close(20)

  return
end

!***********************************************************************

subroutine Read_optic_data(dv0bdcF, Energ_s, Hubb_a, Hubb_d, iaprotoi, lmax_probe, lmaxabs_t,  m_hubb, mpinodes, mpirank0, &
               multi_run, n_atom_0, n_atom_ind, n_atom_proto, n_multi_run, natome, nenerg_s, nlm_pot, nlm_probe, &
               nomfich_optic_data, nptmoy, nptmoy_out, nrm, nspin, nspinp, Rmtg, Rmtsd, rsato, rsbdc, rsbdc_out, Taull_opt, &
               V_hubb, V0muf, Vcato, Vhbdc, Vhbdc_out, Vxcato, VxcbdcF, VxcbdcF_out)

  use declarations
  implicit none
  include 'mpif.h'

  character(len=132):: nomfich_optic_data

  integer:: ie, lmax_probe, lmaxabs_t,  m_hubb, mpierr, mpinodes, mpirank0, multi_run, n_atom_0, n_atom_ind, n_atom_proto, &
            n_multi_run, natome, ndim, nenerg_s, nlm_pot, nlm_probe, nptmoy, nptmoy_out, nrm, nspin, nspinp

  integer, dimension(natome):: iaprotoi

  complex(kind=db), dimension(nenerg_s,nlm_probe,nlm_probe,nspinp,nspinp):: Taull_opt
  complex(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb

  logical:: Hubb_a, Hubb_d

  real(kind=db):: Rmtsd, rsbdc, rsbdc_out, V0muf, Vhbdc, Vhbdc_out
  real(kind=db), dimension(nenerg_s):: Energ_s
  real(kind=db), dimension(0:n_atom_proto):: Rmtg
  real(kind=db), dimension(nspin):: dv0bdcF, VxcbdcF, VxcbdcF_out
  real(kind=db), dimension(0:nrm,nlm_pot,n_atom_0:n_atom_ind):: Vcato
  real(kind=db), dimension(0:nrm,nlm_pot,nspin,n_atom_0:n_atom_ind):: Vxcato
  real(kind=db), dimension(0:nrm,n_atom_0:n_atom_ind):: rsato
  real(kind=db), dimension(-m_hubb:m_hubb,-m_hubb:m_hubb,nspinp,nspinp):: V_hubb_i, V_hubb_r
  real(kind=db), dimension(nlm_probe,nlm_probe,nspinp,nspinp):: Taull_opt_i, Taull_opt_r

  if( multi_run == 1 ) open(20, file = nomfich_optic_data )

  if( mpirank0 == 0 ) then

    read(20,*)
    read(20,*)
    read(20,*) lmax_probe, lmaxabs_t, nptmoy, nptmoy_out
    read(20,'(10i5)') iaprotoi(:)

    do ie = 1,nenerg_s
      read(20,*) Energ_s(ie)
      read(20,*) Taull_opt(ie,:,:,:,:)
    end do
    Energ_s(:) =  Energ_s(:) / Rydb

    read(20,*) Hubb_a, Hubb_d
    read(20,*) Rmtg(:), Rmtsd
    read(20,*) dv0bdcF(:), rsbdc, rsbdc_out, Vhbdc, Vhbdc_out, VxcbdcF(:), VxcbdcF_out(:), V0muf
    read(20,*) rsato(0:nrm,:)
    read(20,*) Vcato(0:nrm,:,:)
    read(20,*) Vxcato(0:nrm,:,:,:)
    if( Hubb_a ) read(20,*) V_hubb(:,:,:,:)

    if( multi_run == n_multi_run ) close(20)

 endif

  if( mpinodes == 1 ) return

  ndim = ( nlm_probe * nspinp )**2

  do ie = 1,nenerg_s

    if( mpirank0 == 0 ) then
      Taull_opt_r(:,:,:,:) = real( Taull_opt(ie,:,:,:,:), db )
      Taull_opt_i(:,:,:,:) = aimag( Taull_opt(ie,:,:,:,:) )
    endif

    call MPI_Bcast(Taull_opt_r,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(Taull_opt_i,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if( mpirank0 /= 0) Taull_opt(ie,:,:,:,:) = cmplx( Taull_opt_r(:,:,:,:), Taull_opt_i(:,:,:,:), db )

  end do

    call MPI_Bcast(nptmoy,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(nptmoy_out,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(lmax_probe,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(lmaxabs_t,1,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(iaprotoi,natome,MPI_INTEGER,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Hubb_a,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Hubb_d,1,MPI_LOGICAL,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Energ_s,nenerg_s,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
  ndim = n_atom_proto + 1
    call MPI_Bcast(Rmtg,ndim,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Rmtsd,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(dv0bdcF,nspin,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(rsbdc,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(rsbdc_out,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(V0muf,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Vhbdc,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(Vhbdc_out,1,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(VxcbdcF,nspin,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(VxcbdcF_out,nspin,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)

  ndim = ( nrm + 1 ) * ( n_atom_ind - n_atom_ind + 1 )
    call MPI_Bcast(rsato,ndim,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
  ndim = ( nrm + 1 ) * nlm_pot * ( n_atom_ind - n_atom_ind + 1 )
    call MPI_Bcast(Vcato,ndim,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
  ndim = ( nrm + 1 ) * nlm_pot * nspin * ( n_atom_ind - n_atom_ind + 1 )
    call MPI_Bcast(Vxcato,ndim,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)

  ndim = ( (2 * m_hubb + 1) * nspinp )**2
  if( mpirank0 == 0 ) then
    V_hubb_r(:,:,:,:) = real( V_hubb(:,:,:,:), db )
    V_hubb_i(:,:,:,:) = aimag( V_hubb(:,:,:,:) )
  endif
    call MPI_Bcast(V_hubb_i,ndim,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
    call MPI_Bcast(V_hubb_r,ndim,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)
  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  if( mpirank0 /= 0) V_hubb(:,:,:,:) = cmplx( V_hubb_r(:,:,:,:), V_hubb_i(:,:,:,:), db )

  return
end

!***********************************************************************

subroutine Recup_nenerg(Energ_max,mpirank0,nenerg_s,nomfich_tddft_data)

  use declarations
  implicit none
  include 'mpif.h'

  integer:: istat, mpierr, mpirank0, nenerg_s

  character(len=132):: nomfich_tddft_data

  real(kind=db):: E_cut, Energ_max

  if( mpirank0 == 0 ) then

    open(20, file = nomfich_tddft_data, status='old', iostat=istat)
    if( istat /= 0 ) call write_open_error(nomfich_tddft_data,istat,1)

    read(20,*) E_cut, nenerg_s, Energ_max
    Energ_max = Energ_max / Rydb

    Close(20)

  endif

  call MPI_Bcast(nenerg_s,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(Energ_max,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  return
end

!***********************************************************************

subroutine Recup_optic_dim(E_cut,Full_atom,iaabsi,iprabs,mpirank0,natome,nenerg_s,nlm_pot,nlm_probe,nomfich_optic_data)

  use declarations
  implicit none
  include 'mpif.h'

  character(len=132):: nomfich_optic_data

  integer:: iaabsi, iprabs, istat, mpierr, mpirank0, natome, nenerg_s, nlm_pot, nlm_probe

  logical:: Full_atom

  real(kind=db):: E_cut

  if( mpirank0 == 0 ) then

    open(20, file = nomfich_optic_data, status='old', iostat=istat)
    if( istat /= 0 ) call write_open_error(nomfich_optic_data,istat,1)

    read(20,*) E_cut, Full_atom
    read(20,*) iaabsi, iprabs, natome, nenerg_s, nlm_pot, nlm_probe

    Close(20)

  endif

  call MPI_Bcast(E_cut,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(Full_atom,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(iaabsi,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(iprabs,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(natome,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(nenerg_s,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(nlm_pot,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(nlm_probe,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  return
end

!***********************************************************************

subroutine Read_tddft_data(E_cut,Energ,imag_taull,mpirank,multi_run,n_multi_run,nbseuil,nenerg,nlmamax, &
                   nlmamax_u,nomfich_tddft_data,nspin,nspino,rof0)

  use declarations
  implicit real(kind=db) (a-h,o-z)
  include 'mpif.h'

  character(len=132):: nomfich_tddft_data

  integer:: ie, istat, mpirank, multi_run, n_multi_run, nbseuil, ndim, nenerg, nlmamax, nlmamax_u, nspin, nspino

  complex(kind=db), dimension(nenerg,nlmamax,nspin,nspino,nbseuil):: rof0

  real(kind=db):: E_cut
  real(kind=db), dimension(nenerg):: Energ
  real(kind=db), dimension(nlmamax,nspin,nspino,nbseuil):: rof0_i, rof0_r
  real(kind=db), dimension(nenerg,nlmamax_u,nspin,nlmamax_u,nspin):: imag_taull
  real(kind=db), dimension(nlmamax_u,nspin,nlmamax_u,nspin):: imag_taull_e

  if( mpirank == 0 ) then

    if( multi_run == 1 ) then
      open(20, file = nomfich_tddft_data, status='old', iostat=istat)
      if(istat/=0) call write_open_error(nomfich_tddft_data,istat,1)
    endif

    read(20,*) E_cut

    do ie = 1,nenerg
      read(20,*) Energ(ie)
      read(20,*) imag_taull(ie,:,:,:,:)
    end do

    do ie = 1,nenerg
      read(20,*)
      read(20,*) rof0(ie,:,:,:,:)
    end do

    E_cut = E_cut / Rydb
    Energ(:) = Energ(:) / Rydb

    if( multi_run == n_multi_run ) close(20)

  endif

  call MPI_Bcast(E_cut,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  call MPI_Bcast(Energ,nenerg,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

  call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

  ndim = ( nlmamax_u * nspin )**2

  do ie = 1,nenerg

    if( mpirank == 0 ) imag_taull_e(:,:,:,:)= imag_taull(ie,:,:,:,:)

    call MPI_Bcast(imag_taull_e,ndim,MPI_REAL8,0,MPI_COMM_WORLD, mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if( mpirank /= 0) imag_taull(ie,:,:,:,:) = imag_taull_e(:,:,:,:)

  end do

  ndim = nlmamax * nspin * nspino * nbseuil

  do ie = 1,nenerg

    if( mpirank == 0 ) then
      rof0_r(:,:,:,:) = real( rof0(ie,:,:,:,:),db )
      rof0_i(:,:,:,:) = aimag( rof0(ie,:,:,:,:) )
    endif

    call MPI_Bcast(rof0_r,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(rof0_i,ndim,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)

    call MPI_BARRIER(MPI_COMM_WORLD,mpierr)

    if( mpirank /= 0 ) rof0(ie,:,:,:,:) = cmplx( rof0_r(:,:,:,:), rof0_i(:,:,:,:),db )

  end do

  return
end

