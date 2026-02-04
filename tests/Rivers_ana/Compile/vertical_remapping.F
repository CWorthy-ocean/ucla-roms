      module vertical_remapping

      ! vertical remapping module
      ! --------------------

#include "cppdefs.opt"


      use dimensions, only: nz

      implicit none

      private

      public remap_src_to_grid

      contains

! ----------------------------------------------------------------------
      subroutine remap_src_to_grid(N_src, H_src, t_src, N_tgt, H_tgt, t_tgt)

      ! Perform conservative remapping of a depth profile from the source grid to the
      ! target (current) grid using piecewise parabolic reconstruction.

      implicit none

      integer, intent(in)                 :: N_src
      real, dimension(N_src), intent(in)  :: H_src
      real, dimension(N_src), intent(in)  :: t_src
      integer, intent(in)                 :: N_tgt
      real, dimension(N_tgt), intent(in)  :: H_tgt
      real, dimension(N_tgt), intent(out) :: t_tgt

      ! local
      integer :: k,k_new,trc,idx
      real, dimension(N_src+1) :: interface_values
      real :: total_depth_src, total_depth_tgt
      real, dimension(N_src) :: a0, a1, a2
      real, dimension(N_src+1) :: z_orig_interfaces
      real, dimension(N_src) :: H_orig
      real :: H_src_tot

      integer :: curr_src_index
      integer, dimension(N_tgt) :: tgt_index_start
      integer, dimension(N_tgt) :: tgt_index_end

      ! Keep track of how far up we are in the water column
      real :: curr_tgt_height, curr_src_height

      ! Keep track of where in each H_orig cell the interfaces of H_tgt are
      real, dimension(N_tgt) :: tgt_frac_start
      real, dimension(N_tgt) :: tgt_frac_end

      ! This holds the integrated_values
      real, dimension(N_tgt) :: definite_integral

      ! Temporary array to hold the remapped fluxes
      real, dimension(N_tgt) :: t_tmp

      ! Scalars to hold the tracer mass in each cell
      real :: total_t_src
      real :: total_t_tgt
      real :: total_t_diff
      real :: cellwise_mass_frac


      tgt_frac_start(:) = 0.0
      tgt_frac_end(:)   = 0.0

      definite_integral(:) = 0.0
      t_tmp(:) = 0.0
      t_tgt(:) = 0.0

      call calc_interface_values(N_src, H_src, t_src, interface_values)

      ! Compare total depths and expand/contract current grid to match original depth
      total_depth_src = 0
      total_depth_tgt = 0
      do k=1,N_tgt
        total_depth_tgt = total_depth_tgt + H_tgt(k)
      enddo

      do k=1,N_src
        total_depth_src = total_depth_src + H_src(k)
      enddo

      total_t_src = 0.0
      do k=1,N_src
          a0(k) = interface_values(k)
          a1(k) = 6*t_src(k) - 4*interface_values(k) - 2*interface_values(k+1)
          a2(k) = 3*(interface_values(k) + interface_values(k+1) - 2*t_src(k))
          total_t_src = total_t_src + t_src(k)*H_src(k)
      enddo

      ! These are vector copies of the source and target column.
      do k = 1,N_src
        H_orig(k) = H_src(k)
      enddo

      ! Match the thicknesses of the source and target columns
      do k = 1,N_src
        H_orig(k) = H_orig(k) * (total_depth_tgt / total_depth_src)
      enddo
      !! The division and multiplication in the above step may not make the total column thicknesses
      !! match exactly, so one more really small correction will be added to the top layer.
      H_src_tot = 0
      do k = 1,N_src
        H_src_tot = H_src_tot + H_orig(k)
      enddo
      H_orig(N_src) = H_orig(N_src) + total_depth_tgt - H_src_tot

      !! Also keep track of the original interface depths.
      z_orig_interfaces(1) = 0
      do k = 1,N_src
        z_orig_interfaces(k+1) = z_orig_interfaces(k) + H_orig(k)
      enddo

      ! Keep track of how far up we are in the water column
      tgt_index_start = 1
      tgt_index_end = 1
      curr_tgt_height = H_tgt(1)
      curr_src_height = H_orig(1)
      curr_src_index = 1

      ! Sort the source and target cell heights so we know how to perform the integration
      do k_new = 1,N_tgt-1

        do while (curr_tgt_height > curr_src_height)
          curr_src_index = curr_src_index + 1
          curr_src_height = curr_src_height + H_orig(curr_src_index)
        enddo

        tgt_index_end(k_new) = curr_src_index
        tgt_index_start(k_new+1) = curr_src_index

        tgt_frac_end(k_new) = (curr_tgt_height - z_orig_interfaces(curr_src_index)) / H_orig(curr_src_index)
        tgt_frac_start(k_new+1) = tgt_frac_end(k_new)

        curr_tgt_height = curr_tgt_height + H_tgt(k_new+1)

      enddo

      ! Remap by integrating the piecewise parabolas from the source grid
      tgt_index_end(N_tgt) = N_src
      tgt_frac_end(N_tgt)  = 1

      total_t_tgt = 0
      do k = 1,N_tgt

        do idx = tgt_index_start(k),tgt_index_end(k)

          if (tgt_index_start(k) == tgt_index_end(k)) then
            definite_integral(k) = integrate(a0(idx), a1(idx), a2(idx), tgt_frac_start(k), tgt_frac_end(k))*H_src(idx)
          elseif ((idx == tgt_index_start(k)) .and. (idx < tgt_index_end(k))) then
            definite_integral(k) = integrate(a0(idx), a1(idx), a2(idx), tgt_frac_start(k), 1.0)*H_src(idx)
          elseif ((idx < tgt_index_end(k)) .and. (idx > tgt_index_start(k))) then
            definite_integral(k) = definite_integral(k) + integrate(a0(idx), a1(idx), a2(idx), 0.0, 1.0)*H_src(idx)
          else
            definite_integral(k) = definite_integral(k) +
     &        integrate(a0(idx), a1(idx), a2(idx), 0.0, tgt_frac_end(k))*H_src(idx)
          endif

        enddo
        t_tmp(k) = definite_integral(k) / H_tgt(k)
        total_t_tgt = total_t_tgt + t_tmp(k)*H_tgt(k)
      enddo

      ! Correction to ensure exact tracer conservation in the remapping
      ! We distribute the correction over the entire column, where the cellwise amount
      ! depends on how much of the column mass the cell contains
      total_t_diff = total_t_tgt - total_t_src

      if (total_t_tgt /= 0.0) then
        do k = 1,N_tgt
          cellwise_mass_frac = t_tmp(k) / total_t_tgt
          t_tgt(k) = t_tmp(k) - total_t_diff*cellwise_mass_frac
        enddo
      endif

      end subroutine remap_src_to_grid
! ----------------------------------------------------------------------
      real function integrate(a0, a1, a2, z0, z1) ![
      ! Definite integration of a parabola

      implicit none

      real :: a0, a1, a2, z0, z1
      real :: one_third = 0.3333333333

      integrate = a0 * (z1 - z0) + 0.5*a1*(z1**2 - z0**2) +
     &  one_third*a2*(z1**3 - z0**3)

      end function integrate !]
! ----------------------------------------------------------------------
      subroutine calc_interface_values(N_src, H, arr, int_vals)  ![

      ! Interpolate/extrapolate values from cell centers to the cell interfaces

      implicit none

      integer, intent(in) :: N_src
      real, dimension(N_src), intent(in)    :: H
      real, dimension(N_src), intent(in)    :: arr
      real, dimension(N_src+1), intent(out) :: int_vals

      real, dimension(0:3) :: B
      real, dimension(0:3,0:3) :: M
      real :: h_b
      real :: h_t
      real :: iH
      real :: Ts_bot
      real :: Ts_top

      integer :: k, kk, kkp, kkp1, trc
      integer :: ord = 3

      !!                  Get tracer values at layer interfaces                !!
      !!  Calculated with Implicit Fourth-order scheme using Thomas algorithm  !!
      !!                   (White and Adcroft, 2008, Eq. 46)                   !!

      !! Bottom boundary extrapolation
      !! Calculated using a polynomial of order "ord" over the bottom (ord+1) cells

      h_b = 0
      h_t = H(1)
      iH  = 1.0 / (h_t - h_b)

      do k = 0,ord
        do kk = 0,ord
          kkp1 = kk+1
          M(k,kk) = (1.0 / kkp1) * iH * (h_t**(kkp1) - h_b**(kkp1))
        enddo

        h_b = h_b + H(k+1)
        h_t = h_t + H(k+2)
        iH = 1.0/(h_t - h_b)
        B(k) = arr(k+1)
      enddo

      Ts_bot = gauss(M,B,ord)

      !! Top boundary extrapolation
      h_b = 0
      h_t = H(N_src)
      iH = 1.0 / (h_t - h_b)
      do k = 0,ord
        do kk = 0,ord
          kkp1 = kk+1
          M(k,kk) = (1.0 / kkp1) * iH * (h_t**kkp1 - h_b**kkp1)
        enddo

        h_b = h_b + H(N_src-k)
        h_t = h_t + H(N_src-1-k)
        iH = 1.0 / (h_t - h_b)
        B(k) = arr(N_src-k)
      enddo

      Ts_top = gauss(M,B,ord)
      !!!!!!!!!!!!!!!!!!!!!!

      call thomas_PPM(N_src, Ts_bot, Ts_top, arr, H, int_vals)

      end subroutine calc_interface_values  !]
! ----------------------------------------------------------------------
      real function gauss(M, b, ord) ![

      !! This routine uses Gaussian Elimination to reduce the matrix M to row-echelon form (lower triangular).
      !! It only returns the first element of the solution vector.

      implicit none

      real, dimension(0:3,0:3), intent(inout) :: M
      real, dimension(0:3), intent(inout) :: b
      integer, intent(in) :: ord
      real :: ratio
      integer :: i,j,k

      !! Starting from the last row, we're going to work upwards to
      !! make M a lower triangular matrix
      do i = ord,1,-1
      !! Iterate over all rows above i
        do j = 0,(i-1)
          ratio = M(j,i) / M(i,i)

          !! Iterate over all columns up to i.
          !! This operation basically is multiplying the ith row by "ratio" and adding it to the jth row
          do k = 0,i
            M(j,k) = M(j,k) - ratio*M(i,k)
          enddo
          b(j) = b(j) - ratio*b(i)
        enddo
      enddo

      gauss = (b(0) / M(0,0))

      end function gauss !]
! ----------------------------------------------------------------------
      subroutine thomas_PPM(N_src, Ts_bot, Ts_top, arr, H, int_vals)  ![

      ! Extremely efficient tridiagonal matrix solver for cases where we
      ! have a diagonally dominant input matrix

      implicit none

      integer, intent(in) :: N_src
      real, intent(in) :: Ts_bot
      real, intent(in) :: Ts_top
      real, dimension(N_src),   intent(in)  :: arr
      real, dimension(N_src),   intent(in)  :: H
      real, dimension(N_src+1), intent(out) :: int_vals

      integer :: n,t,k
      real, dimension(N_src+1) :: a, b, c, d, cp, dp, x
      real :: h0, h1, alpha, beta, d1, d2

      n = N_src+1

      ! Solution for top and bottom boundaries
      b(1) = 1.0
      b(n) = 1.0
      d(1) = Ts_bot
      d(n) = Ts_top
      a(n) = 0
      c(1) = 0

      ! Fill interior of tridiagonal matrix
      do k=2,(n-1)
        h0 = H(k-1)
        h1 = H(k)

        alpha = (h1**2) / ((h0 + h1)**2)
        beta  = (h0**2) / ((h0 + h1)**2)
        d1 = 2*(h1**2)*(h1**2 + 2*h0**2 + 3*h0*h1) / ((h0+h1)**4)
        d2 = 2*(h0**2)*(h0**2 + 2*h1**2 + 3*h0*h1) / ((h0+h1)**4)

        a(k) = alpha
        b(k) = 1.0
        c(k) = beta

        d(k) = d1*arr(k-1) + d2*arr(k)
      enddo

      ! Solve using Thomas algorithm
      cp(1) = c(1) / b(1)
      dp(1) = d(1) / b(1)

      do k = 2,(n-1)
        cp(k) = c(k) / (b(k) - a(k)*cp(k-1))
        dp(k) = (d(k) - a(k)*dp(k-1)) / (b(k) - a(k)*cp(k-1))
      enddo
      dp(n) = (d(n) - a(n)*dp(n-1)) / (b(n) - a(n)*cp(n-1))

      int_vals(n) = d(n)
      do k = (n-1),1,-1
        int_vals(k) = dp(k) - cp(k)*int_vals(k+1)
      enddo

      end subroutine thomas_PPM !]
!-------------------------------------------------------------------

      end module vertical_remapping
