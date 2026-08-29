module mpi_masking_mod
#include "cppdefs.opt"
! Runtime MPI-masking tiling, ported from New-tools/optimal_part.F
! and New-tools/roms_part.F (check_mask, partition_mask, neighbors).
!
! When MPI_MASKING is enabled, ROMS chooses NP_XI/NP_ETA at startup from
! the land mask, then MPI_Comm_split so ocean_grid_comm has nparts ranks.
! Extra launched ranks exit.  MPI_COMM_WORLD itself cannot be resized.

  implicit none
  private

#ifdef MPI_MASKING
  public :: find_optimal_tiling
  public :: build_mpi_masking_map
  public :: get_mpi_masking_topology
#endif

  character(len=15) :: module_name = "mpi_masking_mod"

#ifdef MPI_MASKING
  integer(kind=4), allocatable :: node_map(:,:)
  integer(kind=4), allocatable :: npi_c(:), npj_c(:)
  integer(kind=4) :: nparts_saved = 0
  integer(kind=4), parameter :: invalid_coarse = 100000
#endif

contains

#ifdef MPI_MASKING
!----------------------------------------------------------------------
  subroutine find_optimal_tiling(coarse_frc)
    ! Search candidate (NP_XI, NP_ETA) pairs as in optimal_part, pick the
    ! tiling with the lowest horizontal points per tile that does not
    ! trigger the coarse-grid / Parallel-IO warnings, then split MPI so
    ! ocean_grid_comm has that many ranks.
    ! coarse_frc: the run reads coarse-resolution forcing (interp_*_frc),
    ! so tilings with coarse-grid caveats are unusable rather than merely
    ! warned about.

    use param, only: mynode, nsize, nnodes, NP_XI, NP_ETA, ocean_grid_comm,&
    &LLm, MMm
    use dimensions, only: grdname, analytical_grid
    use error_handling_mod, only: error_log
    use mpi_f08, only: MPI_Comm_split, MPI_Comm_size, MPI_Comm_rank,&
    &MPI_UNDEFINED, MPI_Finalize, MPI_COMM_WORLD
    use netcdf, only: nf90_nowrite, nf90_noerr, nf90_open, nf90_close,&
    &nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var
    implicit none

    logical, intent(in) :: coarse_frc
    character(len=19) :: sr_name = "find_optimal_tiling"
    integer(kind=4) :: ncid, dimid, varid, ierr, color
    integer(kind=4) :: gnx, gny, gnx_rho, gny_rho
    integer(kind=4) :: npx0, npx1, npy0, npy1, npx, npy
    integer(kind=4) :: nx, ny, nxc, nyc, i, k, best
    integer(kind=4) :: launched
    real(kind=8) :: cores, ocean, prct
    real(kind=8), allocatable :: msk(:,:), active(:,:)
    integer(kind=4) :: nparts_loc
    logical :: fallback

#ifndef PARALLEL_IO
    call error_log%raise_global(&
    &context=module_name//'/'//sr_name,&
    &info='MPI_MASKING requires PARALLEL_IO (joined grid file).')
    call error_log%abort_check()
    return
#endif
    if (analytical_grid) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name,&
      &info='MPI_MASKING cannot be used with ANA_GRID.')
      call error_log%abort_check()
      return
    endif

    launched = nsize
    cores = real(launched, kind=8)
    prct = 0.0_8

    ierr = nf90_open(trim(grdname), nf90_nowrite, ncid)
    if (ierr /= nf90_noerr) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name,&
      &info='error opening grid file '//trim(grdname))
      call error_log%abort_check()
      return
    endif

    ierr = nf90_inq_dimid(ncid, 'xi_rho', dimid)
    ierr = nf90_inquire_dimension(ncid, dimid, len=gnx_rho)
    ierr = nf90_inq_dimid(ncid, 'eta_rho', dimid)
    ierr = nf90_inquire_dimension(ncid, dimid, len=gny_rho)
    gnx = gnx_rho - 2
    gny = gny_rho - 2

    allocate(msk(0:gnx+1, 0:gny+1))
    msk = 1.0_8
    ierr = nf90_inq_varid(ncid, 'mask_rho', varid)
    ierr = nf90_get_var(ncid, varid, msk)
    ierr = nf90_close(ncid)
    if (ierr /= nf90_noerr) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name,&
      &info='error reading mask_rho from '//trim(grdname))
      call error_log%abort_check()
      return
    endif

    ocean = sum(msk) / real(gnx*gny, kind=8)
    if (ocean <= 0.0_8) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name,&
      &info='mask_rho has no ocean points.')
      call error_log%abort_check()
      return
    endif

    npx1 = ceiling(cores/ocean)
    npy1 = ceiling(cores/ocean)
    npx0 = 1
    npy0 = 1

    allocate(active(1:(npx1*npy1), 1:9))
    active(:,:) = 0.0_8
    i = 0
    do npy = npy0, npy1
      do npx = npx0, npx1
        if ((npx*npy >= cores) .and. (npx*npy <= cores/ocean) .and.&
        &(npx < gnx) .and. (npy < gny)) then
          call check_mask(npx, npy, gnx, gny, msk, nx, ny, nxc, nyc,&
          &nparts_loc)
          ! nparts_loc == 0 flags a decomposition check_mask could not
          ! partition at all; it must never enter the candidate list.
          if ((nparts_loc >= 1) .and. (nparts_loc <= cores) .and.&
          &(nparts_loc >= (cores*prct))) then
            i = i + 1
            active(i,1) = real(nparts_loc, kind=8)
            active(i,2) = real(npx, kind=8)
            active(i,3) = real(npy, kind=8)
            active(i,4) = real(npx*npy, kind=8) / real(nparts_loc, kind=8)
            active(i,5) = real(nx*ny, kind=8)
            active(i,6) = real(nx, kind=8)
            active(i,7) = real(ny, kind=8)
            active(i,8) = real(nxc, kind=8)
            active(i,9) = real(nyc, kind=8)
          endif
        endif
      enddo
    enddo

    if (i > 0) call sort_active(active(1:i,:))

    ! After sort, largest tiles are first; smallest (most efficient) are last.
    ! Prefer a tiling free of coarse-grid caveats. Those caveats only matter
    ! for runs reading coarse-resolution forcing (interp_*_frc), so when the
    ! run does not (coarse_frc false) and every candidate carries one
    ! (typical for small grids), fall back to the most efficient candidate
    ! rather than refusing to run.
    best = 0
    do k = i, 1, -1
      if (tiling_has_warning(active(k,:))) cycle
      best = k
      exit
    enddo

    fallback = (best == 0) .and. (i > 0) .and. (.not. coarse_frc)
    if (fallback) then
      do k = i, 1, -1
        if ((nint(active(k,8)) == invalid_coarse) .or.&
        &(nint(active(k,9)) == invalid_coarse)) cycle
        best = k
        exit
      enddo
      if (best == 0) best = i
    endif

    if (best == 0) then
      if (mynode == 0) call report_tiling_scan(active, i, 0, launched)
      if (i > 0) then
        call error_log%raise_global(&
        &context=module_name//'/'//sr_name,&
        &info='every candidate MPI_MASKING tiling is incompatible with'//&
        &' coarse-resolution forcing (interp_*_frc enabled); disable'//&
        &' forcing interpolation or relaunch with a different number'//&
        &' of ranks.')
      else
        call error_log%raise_global(&
        &context=module_name//'/'//sr_name,&
        &info='no candidate MPI_MASKING tiling for the launched core'//&
        &' count; relaunch with a different number of ranks.')
      endif
      call error_log%abort_check()
      return
    endif

    nparts_saved = int(active(best,1))
    NP_XI  = int(active(best,2))
    NP_ETA = int(active(best,3))
    nnodes = nparts_saved
    LLm = gnx
    MMm = gny

    if (mynode == 0) then
      call report_tiling_scan(active, i, best, launched)
      if (fallback) then
        write(*,'(a)') 'MPI_MASKING WARNING: every candidate tiling has'//&
        &' coarse-grid caveats; selected the most efficient one anyway.'
        write(*,'(a)') 'MPI_MASKING WARNING: this run cannot read'//&
        &' coarse-grid inputs with this tiling.'
      endif
      if (nparts_saved < launched) then
        write(*,'(a,i0,a,i0,a)')&
        &'MPI_MASKING :: splitting communicator to ', nparts_saved,&
        &' ranks (launched ', launched, '). Extra ranks will exit.'
      endif
    endif

    deallocate(msk, active)

    if (mynode < nparts_saved) then
      color = 0
    else
      color = MPI_UNDEFINED
    endif
    call MPI_Comm_split(MPI_COMM_WORLD, color, mynode, ocean_grid_comm, ierr)

    if (mynode >= nparts_saved) then
      call MPI_Finalize(ierr)
      stop
    endif

    call MPI_Comm_size(ocean_grid_comm, nsize, ierr)
    call MPI_Comm_rank(ocean_grid_comm, mynode, ierr)
    nnodes = nsize
  end subroutine find_optimal_tiling

!----------------------------------------------------------------------
  logical function tiling_has_warning(row)
    real(kind=8), intent(in) :: row(:)
    tiling_has_warning = .false.
    if ((nint(row(8)) == invalid_coarse) .or.&
    &(nint(row(9)) == invalid_coarse)) then
      tiling_has_warning = .true.
    else if ((row(8) <= 2.0_8) .or. (row(9) <= 2.0_8)) then
      tiling_has_warning = .true.
    endif
  end function tiling_has_warning

!----------------------------------------------------------------------
  subroutine sort_active(active)
    ! Insertion sort by horizontal points per tile, descending
    ! (same as New-tools/optimal_part.F).
    real(kind=8), intent(inout) :: active(:,:)
    integer(kind=4) :: i, j
    real(kind=8) :: key(size(active, 2))

    do i = 2, size(active, 1)
      key = active(i,:)
      j = i - 1
      do while (j >= 1)
        if (active(j,5) < key(5)) then
          active(j+1,:) = active(j,:)
          j = j - 1
        else
          exit
        endif
      enddo
      active(j+1,:) = key
    enddo
  end subroutine sort_active

!----------------------------------------------------------------------
  subroutine report_tiling_scan(active, nfound, best, launched)
    ! Print the optimal_part scan (most efficient candidates first)
    ! and the selected NP_XI, NP_ETA, nparts.
    real(kind=8), intent(in) :: active(:,:)
    integer(kind=4), intent(in) :: nfound, best, launched
    integer(kind=4) :: k, nt, nshow

    write(*,*)
    write(*,'(a)') 'MPI_MASKING :: optimal tiling scan'
    write(*,'(a,i0)') '  launched MPI ranks = ', launched
    write(*,'(a,i0)') '  candidate tilings  = ', nfound
    write(*,*)

    nshow = min(10, nfound)
    nt = 0
    do k = nfound, nfound-nshow+1, -1
      nt = nt + 1
      write(*,'(a8,i0,a1)') 'Tiling #', nt, ':'
      write(*,'(a9,i0)') 'NP_XI  = ', int(active(k,2))
      write(*,'(a9,i0)') 'NP_ETA = ', int(active(k,3))
      write(*,'(i0,a36,i0,a15)') int(active(k,1)),&
      &' cores used for ocean points out of ',&
      &int(active(k,2)*active(k,3)), ' tiles overall.'
      write(*,'(a28,i0)') 'Horizontal points per tile: ', int(active(k,5))
      write(*,'(a27,i0,a1)') "Please run with 'mpirun -n ",&
      &int(active(k,1)), "'"
      if ((nint(active(k,8)) == invalid_coarse) .or.&
      &(nint(active(k,9)) == invalid_coarse)) then
        write(*,'(a)') 'WARNING: This decomposition will '//&
        &'not work with coarse grid variables.'
      else if ((active(k,8) <= 2.0_8) .or. (active(k,9) <= 2.0_8)) then
        write(*,'(a)') 'WARNING: Coarse grid inputs will '//&
        &'not work with Parallel IO.'
      endif
      if (k == best) write(*,'(a)') '  *** selected ***'
      write(*,*)
    enddo

    if (best > 0) then
      write(*,'(a)') 'MPI_MASKING :: selected tiling written to log:'
      write(*,'(a,i0)') 'NP_XI  = ', int(active(best,2))
      write(*,'(a,i0)') 'NP_ETA = ', int(active(best,3))
      write(*,'(a,i0)') 'nparts = ', int(active(best,1))
      write(*,*)
    endif
  end subroutine report_tiling_scan

!----------------------------------------------------------------------
  subroutine check_mask(npx, npy, gnx, gny, mask, nx, ny, nxc, nyc,&
  &nparts_out)
    ! From New-tools/roms_part.F :: check_mask
    integer(kind=4), intent(in) :: npx, npy, gnx, gny
    real(kind=8), intent(in) :: mask(0:,0:)
    integer(kind=4), intent(out) :: nx, ny, nxc, nyc, nparts_out

    real(kind=8) :: msk_mx
    integer(kind=4) :: surplus_x, surplus_y, loc_x, loc_y
    integer(kind=4) :: gnxc, gnyc, surplus_xc, surplus_yc, loc_xc, loc_yc
    integer(kind=4) :: npi, npj, count
    logical :: edge_w, edge_e, edge_s, edge_n
    integer(kind=4), allocatable :: iloc(:,:), jloc(:,:), ilcu(:,:), jlcv(:,:)
    integer(kind=4), allocatable :: ilocc(:,:), jlocc(:,:), ilcuc(:,:), jlcvc(:,:)

    nx = ceiling(1.0_8*gnx/npx)
    ny = ceiling(1.0_8*gny/npy)
    surplus_x = nx*npx - gnx
    surplus_y = ny*npy - gny
    gnxc = gnx/2
    gnyc = gny/2
    nxc = (gnxc+npx-1)/npx
    nyc = (gnyc+npy-1)/npy
    surplus_xc = nxc*npx - gnxc
    surplus_yc = nyc*npy - gnyc

    allocate(iloc(npx,3), jloc(npy,3))
    loc_x = 1
    do npi = 1, npx
      if (npi == 1) then
        iloc(npi,1) = nx - surplus_x/2
      else if (npi == npx) then
        iloc(npi,1) = nx - (surplus_x+1)/2
      else
        iloc(npi,1) = nx
      endif
      iloc(npi,2) = loc_x
      iloc(npi,3) = loc_x + iloc(npi,1) - 1
      loc_x = loc_x + iloc(npi,1)
    enddo
    loc_y = 1
    do npj = 1, npy
      if (npj == 1) then
        jloc(npj,1) = ny - surplus_y/2
      else if (npj == npy) then
        jloc(npj,1) = ny - (surplus_y+1)/2
      else
        jloc(npj,1) = ny
      endif
      jloc(npj,2) = loc_y
      jloc(npj,3) = loc_y + jloc(npj,1) - 1
      loc_y = loc_y + jloc(npj,1)
    enddo

    allocate(ilcu(npx,3), jlcv(npy,3))
    ilcu = iloc
    jlcv = jloc
    ilcu(1,1) = iloc(1,1) - 1
    jlcv(1,1) = jloc(1,1) - 1
    ilcu(2:npx,2) = ilcu(2:npx,2) - 1
    ilcu(:,3) = ilcu(:,3) - 1
    jlcv(2:npy,2) = jlcv(2:npy,2) - 1
    jlcv(:,3) = jlcv(:,3) - 1

    if ((any(iloc(:,2) < 1)) .or. (any(jloc(:,2) < 1)) .or.&
    &(any(ilcu(:,2) < 1)) .or. (any(jlcv(:,2) < 1))) then
      deallocate(iloc, jloc, ilcu, jlcv)
      nx = invalid_coarse
      ny = invalid_coarse
      nxc = invalid_coarse
      nyc = invalid_coarse
      nparts_out = 0
      return
    endif

    allocate(ilocc(npx,3), jlocc(npy,3))
    loc_xc = 1
    do npi = 1, npx
      if (npi == 1) then
        ilocc(npi,1) = nxc - surplus_xc/2
      else if (npi == npx) then
        ilocc(npi,1) = nxc - (surplus_xc+1)/2
      else
        ilocc(npi,1) = nxc
      endif
      ilocc(npi,2) = loc_xc
      ilocc(npi,3) = loc_xc + ilocc(npi,1) - 1
      loc_xc = loc_xc + ilocc(npi,1)
    enddo
    loc_yc = 1
    do npj = 1, npy
      if (npj == 1) then
        jlocc(npj,1) = nyc - surplus_yc/2
      else if (npj == npy) then
        jlocc(npj,1) = nyc - (surplus_yc+1)/2
      else
        jlocc(npj,1) = nyc
      endif
      jlocc(npj,2) = loc_yc
      jlocc(npj,3) = loc_yc + jlocc(npj,1) - 1
      loc_yc = loc_yc + jlocc(npj,1)
    enddo

    allocate(ilcuc(npx,3), jlcvc(npy,3))
    ilcuc = ilocc
    jlcvc = jlocc
    ilcuc(1,1) = ilocc(1,1) - 1
    jlcvc(1,1) = jlocc(1,1) - 1
    ilcuc(2:npx,2) = ilcuc(2:npx,2) - 1
    ilcuc(:,3) = ilcuc(:,3) - 1
    jlcvc(2:npy,2) = jlcvc(2:npy,2) - 1
    jlcvc(:,3) = jlcvc(:,3) - 1

    if ((any(ilocc(:,2) < 1)) .or. (any(jlocc(:,2) < 1)) .or.&
    &(any(ilcuc(:,2) < 1)) .or. (any(jlcvc(:,2) < 1))) then
      nxc = invalid_coarse
      nyc = invalid_coarse
    endif

    count = 0
    edge_w = .false.
    edge_e = .false.
    edge_s = .false.
    edge_n = .false.
    do npj = 1, npy
      do npi = 1, npx
        msk_mx = maxval(mask(iloc(npi,2):iloc(npi,3), jloc(npj,2):jloc(npj,3)))
        if (msk_mx > 0.0_8) then
          count = count + 1
          if (npi == 1)   edge_w = .true.
          if (npi == npx) edge_e = .true.
          if (npj == 1)   edge_s = .true.
          if (npj == npy) edge_n = .true.
        endif
      enddo
    enddo
    ! The PIO boundary decompositions require at least one active (ocean)
    ! tile on each physical edge of the domain; a tiling whose edge tiles
    ! are all land is unusable (PIO aborts with totiosize <= 0).
    if (edge_w .and. edge_e .and. edge_s .and. edge_n) then
      nparts_out = count
    else
      nparts_out = 0
    endif

    deallocate(iloc, jloc, ilcu, jlcv, ilocc, jlocc, ilcuc, jlcvc)
  end subroutine check_mask

!----------------------------------------------------------------------
  subroutine build_mpi_masking_map(npx, npy)
    ! From New-tools/roms_part.F :: partition_mask (mask / node_map half)
    use dimensions, only: grdname
    use error_handling_mod, only: error_log
    use netcdf, only: nf90_nowrite, nf90_noerr, nf90_open, nf90_close,&
    &nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var
    implicit none
    integer(kind=4), intent(in) :: npx, npy
    character(len=22) :: sr_name = "build_mpi_masking_map"
    integer(kind=4) :: ncid, dimid, varid, ierr
    integer(kind=4) :: gnx, gny, gnx_rho, gny_rho
    integer(kind=4) :: nx, ny, surplus_x, surplus_y, loc_x, loc_y
    integer(kind=4) :: npi, npj, count
    integer(kind=4), allocatable :: iloc(:,:), jloc(:,:)
    real(kind=8), allocatable :: mask(:,:)
    real(kind=8) :: msk_mx

    ierr = nf90_open(trim(grdname), nf90_nowrite, ncid)
    if (ierr /= nf90_noerr) then
      call error_log%raise_global(&
      &context=module_name//'/'//sr_name,&
      &info='error opening grid file '//trim(grdname))
      return
    endif
    ierr = nf90_inq_dimid(ncid, 'xi_rho', dimid)
    ierr = nf90_inquire_dimension(ncid, dimid, len=gnx_rho)
    ierr = nf90_inq_dimid(ncid, 'eta_rho', dimid)
    ierr = nf90_inquire_dimension(ncid, dimid, len=gny_rho)
    gnx = gnx_rho - 2
    gny = gny_rho - 2

    nx = ceiling(1.0_8*gnx/npx)
    ny = ceiling(1.0_8*gny/npy)
    surplus_x = nx*npx - gnx
    surplus_y = ny*npy - gny

    allocate(iloc(npx,3), jloc(npy,3))
    loc_x = 1
    do npi = 1, npx
      if (npi == 1) then
        iloc(npi,1) = nx - surplus_x/2
      else if (npi == npx) then
        iloc(npi,1) = nx - (surplus_x+1)/2
      else
        iloc(npi,1) = nx
      endif
      iloc(npi,2) = loc_x
      iloc(npi,3) = loc_x + iloc(npi,1) - 1
      loc_x = loc_x + iloc(npi,1)
    enddo
    loc_y = 1
    do npj = 1, npy
      if (npj == 1) then
        jloc(npj,1) = ny - surplus_y/2
      else if (npj == npy) then
        jloc(npj,1) = ny - (surplus_y+1)/2
      else
        jloc(npj,1) = ny
      endif
      jloc(npj,2) = loc_y
      jloc(npj,3) = loc_y + jloc(npj,1) - 1
      loc_y = loc_y + jloc(npj,1)
    enddo

    allocate(mask(0:gnx+1, 0:gny+1))
    ierr = nf90_inq_varid(ncid, 'mask_rho', varid)
    ierr = nf90_get_var(ncid, varid, mask)
    ierr = nf90_close(ncid)

    if (allocated(node_map)) deallocate(node_map)
    if (allocated(npi_c)) deallocate(npi_c)
    if (allocated(npj_c)) deallocate(npj_c)
    allocate(node_map(0:npy+1, 0:npx+1))
    node_map = -98
    allocate(npi_c(npx*npy), npj_c(npx*npy))
    npi_c = 0
    npj_c = 0

    count = 0
    do npj = 1, npy
      do npi = 1, npx
        msk_mx = maxval(mask(iloc(npi,2):iloc(npi,3), jloc(npj,2):jloc(npj,3)))
        if (msk_mx > 0.0_8) then
          count = count + 1
          node_map(npj, npi) = count
          npi_c(count) = npi
          npj_c(count) = npj
        endif
      enddo
    enddo
    nparts_saved = count
    deallocate(iloc, jloc, mask)
  end subroutine build_mpi_masking_map

!----------------------------------------------------------------------
  subroutine get_mpi_masking_topology(rank, neighbors, subdompos, nnodes_out)
    ! neighbors: N,NE,E,SE,S,SW,W,NW as in make_partial_files (0-based ranks)
    integer(kind=4), intent(in) :: rank
    integer(kind=4), intent(out) :: neighbors(8)
    integer(kind=4), intent(out) :: subdompos(2)
    integer(kind=4), intent(out) :: nnodes_out
    integer(kind=4) :: part, i, j
    integer(kind=4) :: loc_neigh(8)

    part = rank + 1
    i = npi_c(part)
    j = npj_c(part)
    subdompos = (/ i-1, j-1 /)
    loc_neigh = (/ node_map(j+1,i), node_map(j+1,i+1),&
    &node_map(j,i+1), node_map(j-1,i+1),&
    &node_map(j-1,i), node_map(j-1,i-1),&
    &node_map(j,i-1), node_map(j+1,i-1) /)
    neighbors = loc_neigh - 1
    nnodes_out = nparts_saved
  end subroutine get_mpi_masking_topology

#endif /* MPI_MASKING */

end module mpi_masking_mod
