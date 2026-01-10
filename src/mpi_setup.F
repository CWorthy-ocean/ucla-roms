#include "cppdefs.opt"

      subroutine mpi_setup(ierr)
      use param, only:
     &     ieast, iwest, jnorth, jsouth, llm,
     &     lm, mm, mmm, ocean_grid_comm, np_xi, np_eta, nsize
      use hidden_mpi_vars, only:
     &     east_msg_exch, exc_call_count, inode,
     &     isw_corn, jnode, jsw_corn, mynode, nnodes, north_msg_exch,
     &     p_e, p_n, p_ne, p_nw, p_s, p_se, p_sw, p_w,
     &     south_msg_exch, west_msg_exch
      use mpi_f08, only: mpi_comm_world
      use dimensions, only:
     &     eta_rho, eta_v, xi_rho, xi_u, gnx, gny,
     &     npx, npy, grdname, i0, i1, j0, j1, nx, ny
      use netcdf, only:
     &     nf90_open, nf90_global, nf90_noerr, nf90_get_att, nf90_close,
     &     nf90_nowrite
      implicit none
      integer :: ierr,newpart,ncid
      integer,dimension(8) :: neighbors
      integer,dimension(2) :: subdompos
      integer,dimension(4) :: partition
      character(len=128) :: fname
      logical :: new_part
      integer :: offx,offy
      integer :: ijnode2rank


      ! read in gridfile to check if mpi masking is present
# ifdef ANA_GRID
      nnodes=np_xi*np_eta
      new_part = .false.
# else
      ierr=nf90_open(grdname,nf90_nowrite,ncid)
c$$$       if (ierr/=nf90_noerr) call handle_ierr(ierr)
        ierr=nf90_get_att(ncid,nf90_global,"neighbors",neighbors)
        ! new style partitioning has the neighbors global attribute
        if (ierr.eq.nf90_noerr) then
           if (mynode.eq.0) print * ,'neighbors att found',grdname
           new_part = .true.
           ierr=nf90_get_att(ncid,nf90_global,"partition",partition)
           nnodes=partition(2)
           ierr=nf90_get_att(ncid,nf90_global,"neighbors",neighbors)
           ierr=nf90_get_att(ncid,nf90_global,"subdompos",subdompos)
        else
           if (mynode.eq.0) print * ,'neighbors att not found',grdname
           new_part = .false.
           nnodes = np_xi*np_eta
        endif
        ierr=nf90_close(ncid)
#endif

      ! Check whether the number of expected nodes
      ! matches the mpirun or mpiexec argument

      ! Also check npx,npy somewhere

      if (nsize/=nnodes) then
        if (mynode.eq.0) write(*,'(/1x,A,I4,1x,A,I3,A/)')
     &   '### ERROR: mpi_setup: number of MPI-nodes should be',
     &               nnodes, 'instead of ', nsize, '.'
        ierr=99
        return
      endif

!     Determine position of subdomain in logical grid. inodes, jnodes
!     will be needed to determine whether the subdomain is at a physical
!     boundary.
!     Afterwards, determine the mpi ranks of the neighboring subdomains.
!     These are needed to designate MPI ranks for sources of
!     incoming and targets for outgoing messages. Negative numbers are invalid
!     and indicate that there is not a valid subdomain in that direction

      if (new_part) then
        inode=subdompos(1)
        jnode=subdompos(2)
        p_N  = neighbors(1)
        p_NE = neighbors(2)
        p_E  = neighbors(3)
        p_SE = neighbors(4)
        p_S  = neighbors(5)
        p_SW = neighbors(6)
        p_W  = neighbors(7)
        p_NW = neighbors(8)
      else
        jnode=mynode/npx
        inode=mynode-jnode*npx
        p_N = ijnode2rank(inode  ,jnode+1,npx,npy)
        p_NE= ijnode2rank(inode+1,jnode+1,npx,npy)
        p_E = ijnode2rank(inode+1,jnode  ,npx,npy)
        p_SE= ijnode2rank(inode+1,jnode-1,npx,npy)
        p_S = ijnode2rank(inode  ,jnode-1,npx,npy)
        p_SW= ijnode2rank(inode-1,jnode-1,npx,npy)
        p_W = ijnode2rank(inode-1,jnode  ,npx,npy)
        p_NW= ijnode2rank(inode-1,jnode+1,npx,npy)
      endif

!     Determine the sizes of subdomains.
!     The first and last subdomain are possibly smaller to add
!     to the exact size of the domain. Half of the surplus is
!     removed from those. If the surplus (offx,offy) is odd, the
!     last subdomain will be the smallest.

      nx = ceiling(1.0*gnx/npx)
      ny = ceiling(1.0*gny/npy)

        ! Check against original weird integer devision
        if (nx.ne.(gnx+npx-1)/npx )
     &       stop 'integer devision error in mpi_setup'
        if (ny.ne.(gny+npy-1)/npy )
     &       stop 'integer devision error in mpi_setup'


        offx = npx*nx-gnx
        if (inode.eq.0) then
          iSW_corn = 0
        else
          iSW_corn = inode*nx - offx/2
        endif

        if (inode.eq.0) then
          nx = nx - offx/2
        endif
        if (inode.eq.npx-1) then
          nx = nx -(offx+1)/2
        endif

        offy= npy*ny-gny
        if (jnode == 0) then
          jSW_corn = 0
        else
          jSW_corn = jnode*ny - offy/2
        endif
        if (jnode == 0) then
          ny = ny - offy/2
        endif
        if (jnode == NP_ETA-1) then
          ny = ny -(offy+1)/2
        endif

        iwest = 1;  ieast = nx
        jsouth = 1; jnorth = ny

        ! Set sizes for output files. The first and last partial files
        ! include the buffers (for the rho-dimension), the others are
        ! all sized: [1:nx],[1:ny]

        i0 = 1; i1 = nx
        j0 = 1; j1 = ny

        xi_rho = nx
        xi_u   = nx
        if (inode.eq.0) then
          i0 = 0
          xi_rho = nx+1
        endif
        if (inode.eq.npx-1) then
          i1 = nx+1
          xi_rho = nx+1
          xi_u   = nx+1
        endif

        eta_rho = ny
        eta_v   = ny
        if (jnode.eq.0) then
          j0=0
          eta_rho = ny+1
        endif
        if (jnode.eq.npy-1) then
          j1=ny+1
          eta_rho = ny+1
          eta_v   = ny+1
        endif

      end subroutine mpi_setup !]
!--------------------------------------------------------
      function ijnode2rank(inode,jnode,npx,npy) result(rank) ![

! Convert processor-grid indices inode,jnode into MPI rank. Check whether
! the indices are within the bounds of processor grid. Return them back
! into periodicity domain, if there is periodicity in either direction,
! but only if the number of subdomains is greater than one (hence
! periodicity must go through message passing).
      use param, only: np_xi, np_eta
      implicit none

      ! import/export
      integer,intent(in) :: inode,jnode
      integer,intent(in) :: npx,npy

      integer :: rank

      ! local
      integer :: i,j

      i=inode ; j=jnode

#ifdef EW_PERIODIC
        if (i < 0) then
          i=i+npx
        elseif (i > npx-1) then
          i=i-np_xi
        endif
#endif

#ifdef NS_PERIODIC
        if (j < 0) then
          j=j+np_eta
        elseif (j > np_eta-1) then
          j=j-np_eta
        endif
#endif

      if (0<=i .and. i<npx .and. 0<=j .and. j<npy) then
        rank = i + j*npx
      else
        rank = -1 !<- meaning that the rank does not exist
      endif
      end function ijnode2rank
