#include "cppdefs.h"
#ifdef ANA_INITIAL
                                  ! Set initial conditions for momentum
      subroutine ana_init (tile)  ! free surface, and tracer variables
      implicit none               ! using analytical expressions.
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call ana_init_tile( istr,iend,jstr,jend)

      end

      subroutine ana_init_tile (istr,iend,jstr,jend)
      implicit none
# include "param.h"
# include "scalars.h"

      integer istr,iend,jstr,jend, i,j,k,itrc

      real x0,y0

# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "eos_vars.h"

      real  amp,Lt
# if defined TANKINT
      real drho, thickn
# endif

# include "compute_auxiliary_bounds.h"

      do j=jstrR,jendR         ! Set everything (except temperature
        do i=istrR,iendR       ! and salinity) to all-zero state, then
          zeta(i,j,1)=0.       ! modify some of the variables, if a
          ubar(i,j,1)=0.       ! non-trivial initialization required.
          vbar(i,j,1)=0.       ! Note: A code to initialize T [and S]
        enddo                  ! must always be supplied for 3D
      enddo                    ! applications.
# ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            u(i,j,k,1)=0.
            u(i,j,k,2)=0.
            v(i,j,k,1)=0.
            v(i,j,k,2)=0.
#  ifdef NHMG
            w(i,j,k,1)=0.
            w(i,j,k,2)=0.
#  endif
            t(i,j,k,1,:)= 10.
            t(i,j,k,2,:)= 10.
          enddo
        enddo
      enddo
# endif

# if defined TANK
      amp = 0.001
      Lt  = 10.0
      do j=jstrR,jendR
        do i=istrR,iendR
          zeta(i,j,1)=amp*cos(pi*xr(i,j)/Lt)
        enddo
      enddo

      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=20.
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo

# elif defined TANKINT
      drho=10
      thickn=2
      cff=2.*(log((0.99+1.)/(-0.99+1.))/2.)/thickn
      amp=0.5
      Lt=10.
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            t(i,j,k,1,itemp)=rho0-0.5*drho*
     &      tanh(cff*(z_r(i,j,k)+h(i,j)/2-amp*cos(pi*xr(i,j)/Lt)))
     &      - 1000.
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
          enddo
        enddo
      enddo
# else

      'No analytical init defined'

# endif  /* end of configuration selection switch */

! The rest of the code is independent of the specific configuration.

! Analytical problems with open boundaries may require external forcing
! data to be applied at inflow side.  For the problems of with time
! independent inflow conditions it is sufficient to save initial
! condition values at the boundary point into boundary forcing arrays,
! which will remain constant thereafter and be used as the source of
! forcing data.  Note that this automatically guarantees consistency
! between the initial state and side boundary conditions.


# if defined OBC_WEST  || defined OBC_EAST ||\
     defined OBC_SOUTH || defined OBC_NORTH
#  ifdef ANA_BRY
#   ifdef OBC_WEST
      if (WESTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR                     ! Save perimeter values
          zeta_west(j)=zeta(istr-1,j,1)      ! of the initial fields
          ubar_west(j)=ubar(istr  ,j,1)      ! to be used as external
          vbar_west(j)=vbar(istr-1,j,1)      ! forcing data for side
        enddo                                ! boundaries.
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_west(j,k)=u(istr  ,j,k,1)
            v_west(j,k)=v(istr-1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_west(j,k,itemp)=t(istr-1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_EAST
      if (EASTERN_EDGE) then
#    ifdef M2_FRC_BRY
        do j=jstrR,jendR
          zeta_east(j)=zeta(iend+1,j,1)
          ubar_east(j)=ubar(iend+1,j,1)
          vbar_east(j)=vbar(iend+1,j,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do j=jstrR,jendR
#     ifdef M3_FRC_BRY
            u_east(j,k)=u(iend+1,j,k,1)
            v_east(j,k)=v(iend+1,j,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_east(j,k,itemp)=t(iend+1,j,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif


#   ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_south(i)=zeta(i,jstr-1,1)
          ubar_south(i)=ubar(i,jstr-1,1)
          vbar_south(i)=vbar(i,jstr  ,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_south(i,k)=u(i,jstr-1,k,1)
            v_south(i,k)=v(i,jstr  ,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_south(i,k,itemp)=t(i,jstr-1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif

#   ifdef OBC_NORTH
      if (NORTHERN_EDGE) then
#    ifdef M2_FRC_BRY
        do i=istrR,iendR
          zeta_north(i)=zeta(i,jend+1,1)
          ubar_north(i)=ubar(i,jend+1,1)
          vbar_north(i)=vbar(i,jend+1,1)
        enddo
#    endif
#    ifdef SOLVE3D
        do k=1,N
          do i=istrR,iendR
#     ifdef M3_FRC_BRY
            u_north(i,k)=u(i,jend+1,k,1)
            v_north(i,k)=v(i,jend+1,k,1)
#     endif
#     ifdef T_FRC_BRY
            t_north(i,k,itemp)=t(i,jend+1,k,1,itemp)
#     endif
          enddo
        enddo
#    endif
      endif
#   endif
#  else /* alternative to ANA_BRY */
#   ifdef UCLIMATOLOGY
      do j=jstrR,jendR                   ! Save initial data into
        do i=istrR,iendR                 ! climatology for the
          ubclm(i,j)=ubar(i,j,1)         ! subsequent use as inflow
          vbclm(i,j)=vbar(i,j,1)         ! boundary conditions
        enddo
      enddo
#   endif
#   ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
#    ifdef UCLIMATOLOGY
            uclm(i,j,k)=u(i,j,k,1)
            vclm(i,j,k)=v(i,j,k,1)
#    endif
#    ifdef TCLIMATOLOGY
            tclm(i,j,k,itemp)=t(i,j,k,1,itemp)
#    endif
          enddo
        enddo
      enddo
#   endif /* SOLVE3D */
#  endif /* ANA_BRY */
# endif /* at least one of OBC_XXXXX */





# ifdef EXCHANGE
      call exchange2d_3_tile(istr,iend,jstr,jend,
     &                       zeta(START_2D_ARRAY,1),
     &                       ubar(START_2D_ARRAY,1),
     &                       vbar(START_2D_ARRAY,1))
#  ifdef SOLVE3D
      call exchange_2_tile(istr,iend,jstr,jend,
     &                     u(START_2D_ARRAY,1,1), N,
     &                     v(START_2D_ARRAY,1,1), N)
      do itrc=1,NT
        call exchange_tile(istr,iend,jstr,jend,
     &                t(START_2D_ARRAY,1,1,itrc), N)
      enddo
#  endif
#  ifdef ISWAKE
#   ifdef UCLIMATOLOGY
      call exchange2d_2_tile(istr,iend,jstr,jend, ubclm,vbclm)
#   endif
#   ifdef SOLVE3D
#    ifdef UCLIMATOLOGY
      call exchange_2_tile(istr,iend,jstr,jend, uclm,N, vclm,N)
#    endif
#    ifdef TCLIMATOLOGY
      call exchange_tile(istr,iend,jstr,jend,
     &                    tclm(START_2D_ARRAY,1,itemp), N)
#    endif
#   endif
#  elif defined COLD_FILAMENT
      call exchange_tile(istr,iend,jstr,jend, Akv,N+1)

#  endif
# endif  /* EXCHANGE */


# if defined OBC_WEST && (defined NS_PERIODIC || defined MPI)
      if (WESTERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, zeta_west,1)
        call exch_bry_EW_tile (jstr,jend, ubar_west,1)
        call exch_bry_EW_tile (jstr,jend, vbar_west,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, u_west,N)
        call exch_bry_EW_tile (jstr,jend, v_west,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_EW_tile (jstr,jend, t_west(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_EAST && (defined NS_PERIODIC || defined MPI)
      if (EASTERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, zeta_east,1)
        call exch_bry_EW_tile (jstr,jend, ubar_east,1)
        call exch_bry_EW_tile (jstr,jend, vbar_east,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_EW_tile (jstr,jend, u_east,N)
        call exch_bry_EW_tile (jstr,jend, v_east,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_EW_tile (jstr,jend, t_east(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_SOUTH && (defined EW_PERIODIC || defined MPI)
      if (SOUTHERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_NS_tile (istr,iend, zeta_south,1)
        call exch_bry_NS_tile (istr,iend, ubar_south,1)
        call exch_bry_NS_tile (istr,iend, vbar_south,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_NS_tile (istr,iend, u_south,N)
        call exch_bry_NS_tile (istr,iend, v_south,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_NS_tile (istr,iend, t_south(0,1,itrc),N)
        enddo
#  endif
      endif
# endif

# if defined OBC_NORTH && (defined EW_PERIODIC || defined MPI)
      if (NORTHERN_EDGE) then
#  ifdef M2_FRC_BRY
        call exch_bry_NS_tile (istr,iend, zeta_north,1)
        call exch_bry_NS_tile (istr,iend, ubar_north,1)
        call exch_bry_NS_tile (istr,iend, vbar_north,1)
#  endif
#  ifdef M3_FRC_BRY
        call exch_bry_NS_tile (istr,iend, u_north,N)
        call exch_bry_NS_tile (istr,iend, v_north,N)
#  endif
#  ifdef T_FRC_BRY
        do itrc=1,NT
          call exch_bry_NS_tile (istr,iend, t_north(0,1,itrc),N)
        enddo
#  endif
      endif
# endif
      end
#else
      subroutine ana_init_empty
      end
#endif /* ANA_INITIAL */
