module lmd_kpp_mod

#include "cppdefs.opt"

  use open_boundary_mod, only: ew_periodic, ns_periodic
  ! The modules of lmd_kpp_tile.h: here, so that the build's dependency scan sees them
  use alphabeta_mod, only: alphabeta_tile
  use surf_flux, only: stflx, srflx, svstr, sustr, svstr_r, sustr_r, bulk_frc_on
!     bulk_frc-> only: sustr_r, svstr_r, stflx, srflx, only:
  use grid, only: f, rmask, umask, vmask
!     flux_frc-> only: sustr, svstr, stflx, srflx
  use roms_mpi, only: exchange_xxx
  use dimensions, only: inode, jnode
  use param, only: ieast, itemp, iwest, jnorth, jsouth, np_xi, np_eta&
  &, isalt, salinity
  use mixing, only: akt, akv
  use mixing, only: bvf
  use mixing, only: hbls, swr_frac   ! lmd_kpp
  use mixing, only: ghat             ! lmd_nonlocal
  use mixing, only: hbbl             ! lmd_bkpp
  use ocean_vars, only: z_w, hz, v, u, z_r
  use scalars, only: nz, forw_start, g, iic, nstp, vonkar, zob
  implicit none
  private

  ! Called only with lmd_kpp or lmd_bkpp (VERTICAL_MIXING_SETTINGS; lmd_vmix).
  public :: lmd_kpp_tile
  public :: check_kpp_switches

contains

  subroutine lmd_kpp_tile (istr,iend,jstr,jend,       Kv,Kt,Ks,&
  &ustar, Bo,Bosol, hbl,bbl, FX,FE,FE1,&
  &Cr,FC, wrk1,wrk2,&
  &Gm1,dGm1dS,  Gt1,dGt1dS,  Gs1,dGs1dS,&
  &kbls, kmo,&
  &tind)
    ! The KPP boundary layers of the settings lmd_kpp, lmd_bkpp and lmd_nonlocal: one
    ! copy of lmd_kpp_tile.h each (below), chosen here.
    use scalars, only: lmd_kpp, lmd_bkpp, lmd_nonlocal
    implicit none
    integer(kind=4) istr,iend,jstr,jend
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY,0:nz) :: Kv, Kt, Ks
    real(kind=8), dimension(PRIVATE_2D_SCRATCH_ARRAY) :: ustar, Bo, Bosol&
    &, hbl, bbl&
    &, FX, FE, FE1
    real(kind=8), dimension(PRIVATE_1D_SCRATCH_ARRAY,0:nz) :: Cr,FC, wrk1,wrk2
    real(kind=8), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS
    integer(kind=4), dimension(PRIVATE_1D_SCRATCH_ARRAY) :: kbls, kmo
    integer(kind=4) :: tind

    if (lmd_kpp .and. lmd_bkpp .and. lmd_nonlocal) then
      call lmd_kpp_tile_sbn(istr,iend,jstr,jend, Kv,Kt,Ks, ustar, Bo,Bosol, hbl,bbl, FX,FE,FE1,&
      &Cr,FC, wrk1,wrk2, Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS, kbls, kmo, tind)
    elseif (lmd_kpp .and. lmd_bkpp) then
      call lmd_kpp_tile_sb(istr,iend,jstr,jend, Kv,Kt,Ks, ustar, Bo,Bosol, hbl,bbl, FX,FE,FE1,&
      &Cr,FC, wrk1,wrk2, Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS, kbls, kmo, tind)
    elseif (lmd_kpp .and. lmd_nonlocal) then
      call lmd_kpp_tile_sn(istr,iend,jstr,jend, Kv,Kt,Ks, ustar, Bo,Bosol, hbl,bbl, FX,FE,FE1,&
      &Cr,FC, wrk1,wrk2, Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS, kbls, kmo, tind)
    elseif (lmd_kpp) then
      call lmd_kpp_tile_s(istr,iend,jstr,jend, Kv,Kt,Ks, ustar, Bo,Bosol, hbl,bbl, FX,FE,FE1,&
      &Cr,FC, wrk1,wrk2, Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS, kbls, kmo, tind)
    else
      call lmd_kpp_tile_b(istr,iend,jstr,jend, Kv,Kt,Ks, ustar, Bo,Bosol, hbl,bbl, FX,FE,FE1,&
      &Cr,FC, wrk1,wrk2, Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS, kbls, kmo, tind)
    endif
  end subroutine lmd_kpp_tile

! lmd_kpp, lmd_bkpp and lmd_nonlocal
#define KPP_SURF
#define KPP_BOT
#define KPP_NONLOC
#define LMD_KPP_TILE lmd_kpp_tile_sbn
#include "lmd_kpp_tile.h"
#undef LMD_KPP_TILE
#undef KPP_NONLOC
! lmd_kpp and lmd_bkpp
#define LMD_KPP_TILE lmd_kpp_tile_sb
#include "lmd_kpp_tile.h"
#undef LMD_KPP_TILE
#undef KPP_BOT
! lmd_kpp and lmd_nonlocal
#define KPP_NONLOC
#define LMD_KPP_TILE lmd_kpp_tile_sn
#include "lmd_kpp_tile.h"
#undef LMD_KPP_TILE
#undef KPP_NONLOC
! lmd_kpp
#define LMD_KPP_TILE lmd_kpp_tile_s
#include "lmd_kpp_tile.h"
#undef LMD_KPP_TILE
#undef KPP_SURF
! lmd_bkpp
#define KPP_BOT
#define LMD_KPP_TILE lmd_kpp_tile_b
#include "lmd_kpp_tile.h"
#undef LMD_KPP_TILE
#undef KPP_BOT


subroutine check_kpp_switches

! This routine is to keep track of the status of CPP-switches defined
! locally in this file above. It must be placed here rather than stored
! in a separate file in order to be exposed to relevant CPP-settings.
! This routine does not affect model results, other than signature in
! global attribute "CPPS" in output netCDF files.
  use strings, only: cpps, max_opt_size
  use param, only: mynode
  use error_handling_mod, only: error_log
  use utils_mod, only: lenstr

  implicit none
  integer(kind=4) is,ie
  character(len=19) :: sr_name = "check_kpp_switches"
  character(len=1024) :: error_info = ""

  ie=lenstr(cpps)
  is=ie+2 ; ie=is+10

  write(error_info,*)&
  &'Insufficient length of string "cpps" in file "strings".',&
  &'Increase parameter "max_opt_size" it and recompile.'

  if (ie > max_opt_size) then
    call error_log%raise_global(&
    &info=error_info,&
    &context=sr_name)
  end if
  cpps(is:ie)='<lmd_kpp.F>'
  is=ie+2 ; ie=is+16
  if (ie > max_opt_size) then
    call error_log%raise_global(&
    &info=error_info,&
    &context=sr_name)
  end if
  cpps(is:ie)='INT_AT_RHO_POINTS'
  is=ie+2 ; ie=is+9
  if (ie > max_opt_size) then
    call error_log%raise_global(&
    &info=error_info,&
    &context=sr_name)
  end if
  cpps(is:ie)='SMOOTH_HBL'
  call error_log%abort_check()
end subroutine check_kpp_switches
end module lmd_kpp_mod
