      module step3d_uv_mod
#include "cppdefs.opt"
      implicit none
      private
#ifdef SOLVE3D
# define UPSTREAM_UV
      public :: step3d_uv1, check_step_uv1_switches
      public :: step3d_uv2, check_step_uv2_switches

      contains

      subroutine step3d_uv1(tile)
      use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
      use private_scratch, only: a3d, a2d
      implicit none
      integer tile

# include "compute_tile_bounds.h"
      call step3d_uv1_tile( istr,iend,jstr,jend,    A3d(1,1),A3d(1,2),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                            A2d(1,1),A2d(1,2),A2d(1,3),A2d(1,4),
     &                                              A2d(1,5),A2d(1,6)
# ifdef NHMG
     &                                             ,A3d(1,5)
# endif
     &                                                              )
      end

      subroutine step3d_uv1_tile( istr,iend,jstr,jend, ru,rv,
     &                                            WC,FC,CF,DC,
     &                             UFx,UFe,VFx,VFe, wrk1,wrk2
# ifdef NHMG
     &                                            ,rw
# endif
     &     )
      use dimensions, only: inode, jnode
# ifdef NHMG
      use ocean_vars, only:
     &     nh_ubar, nh_vbar, nh_wcor, nhdu,
     &     nhdv, nhdw, w
      use nhmg, only:
     &     nhmg_solve,  iprec1, iprec2, halo
      use mg_namelist, only: surface_neumann
      use mg_grids, only: mggrid => grid
# endif
      use diagnostics, only: u_prev, v_prev, calc_diag, udiag,
     &  vdiag, icori, dxdyi_u, dxdyi_v, idiss, ipgr, iadv, ivmix,
     &     set_diags_u_4th_adv, set_diags_v_4th_adv, diag_uv
      use surf_flux, only: sustr, svstr
      use coupling, only: r_d, rufrc, rvfrc
      use grid, only:
     &     nx,ny,nz, fomn, umask, vmask, pn, pm, dn_u,
     &     dm_u, dn_v, dm_v, dn_r, dm_r, pn_u, pm_v
#if (defined CURVGRID && defined UV_ADV)
     &     , dmde,dndx
#endif
      use mixing, only: akv
      use param, only:
     &     ieast, iwest, jnorth, jsouth,
     &     lm, mm, padd_e, padd_x, np_xi, np_eta
      use ocean_vars, only: u, v, hz, flxu, flxv, we, wi, z_w
      use scalars, only: n, dt, nnew, nrhs, ntstart
#ifdef EXCHANGE
      use roms_mpi, only: exchange_xxx
#endif
# if defined ADV_ISONEUTRAL && defined EXCHANGE
      use eos_vars, only: drde, drdx
#endif
      implicit none
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,N) :: ru,rv
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: WC,FC,CF,DC
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: UFx,UFe,VFx,VFe,
     &                                                   wrk1,wrk2
      real cff
      real, parameter ::  delta=0.1666666666666667 ! delta=0.125
# ifdef UPSTREAM_UV
     &                  , gamma=0.3333333333333333 ! gamma=0.25
# endif
# ifdef NHMG
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY,0:N) :: rw
      real  Flxw,Uflxw,Vflxw
      real  ub,ut
      real  div,dmax

#endif

# include "compute_auxiliary_bounds.h"

# ifdef DIAGNOSTICS
#  ifdef NHMGDIAG
      if (diag_uv.and.calc_diag) then
        Wdiag(istr:iend,jstr:jend,:,iwprsgr) =rw(istr:iend,jstr:jend,:)
      endif
#  endif
# endif

# ifdef NHMG

#  define FX UFx
#  define FE VFe

      do k=1,N
#  include "compute_horiz_rhs_w_terms.h"
      enddo

#  undef FE
#  undef FX

      do j=jstr,jend

#  include "compute_vert_rhs_w_terms.h"
        do i=istr,iend
          DC(i,0)=dt*pm(i,j)*pn(i,j)
        enddo
        do k=1,N           !! w is dz*w after this
          do i=istr,iend
            w(i,j,k,nnew)=w(i,j,k,nnew) +DC(i,0)*rw(i,j,k)
          enddo
        enddo

      enddo
#  ifdef DIAGNOSTICS_NHMG
      if (calc_diag) call set_diags_w_at_uv1( istr, iend, jstr, jend, DC, rw )
#  endif
# endif  /* NHMG */

      do k=1,N
# include "compute_horiz_rhs_uv_terms.h"
      enddo

# ifdef DIAGNOSTICS
        if (diag_uv .and. calc_diag) then
          call set_diags_u_4th_adv
          call set_diags_v_4th_adv
        endif
# endif

      do j=jstr,jend

# include "compute_vert_rhs_uv_terms.h"

#define IMPLICIT_BOTTOM_DRAG
# ifdef IMPLICIT_BOTTOM_DRAG
#ifdef DIAGNOSTICS
        if (diag_uv) then
          u_prev = u(1:nx,j,:,nnew)
          v_prev = v(1:nx,j,:,nnew)
        endif
#endif
        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i-1,j,N-1))
     &                      /(  Hz(i,j,N  )+Hz(i-1,j,N  )
     &                         +Hz(i,j,N-1)+Hz(i-1,j,N-1))

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i-1,j,N-1))

          cff=1./( 0.5*(Hz(i,j,N)+Hz(i-1,j,N))
     &                   +FC(i,N-1)-min(WC(i,N-1),0.) )

          CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

          DC(i,N)=cff*( u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
     &                                     +dt*sustr(i,j)
     &                                                  )
        enddo
        do k=N-1,2,-1      !--> forward elimination
          do i=istrU,iend
            FC(i,k-1)= 2.*dt*(Akv(i,j,k-1)+Akv(i-1,j,k-1))
     &                        /(  Hz(i,j,k  )+Hz(i-1,j,k  )
     &                           +Hz(i,j,k-1)+Hz(i-1,j,k-1))

            WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i-1,j,k-1))

            cff=1./( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))
     &                           +FC(i,k-1)-min(WC(i,k-1),0.)
     &                             +FC(i,k)+max(WC(i,k),0.)
     &                    -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                      )
            CF(i,k-1)=cff*(   FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=cff*( u(i,j,k,nnew) +DC(i,0)*ru(i,j,k)
     &                 +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
          enddo
        enddo
        do i=istrU,iend
          DC(i,1)=( u(i,j,1,nnew) +DC(i,0)*ru(i,j,1)
     &                      +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                          )/( 0.5*(Hz(i,j,1)+Hz(i-1,j,1))
     &                            +0.5*dt*(r_D(i,j)+r_D(i-1,j))
     &                                +FC(i,1)+max(WC(i,1),0.)
     &                       -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                       )

          u(i,j,1,nnew)=DC(i,1) * 0.5*(Hz(i,j,1)+Hz(i-1,j,1))

          rufrc(i,j)=ru(i,j,1) +dm_u(i,j)*dn_u(i,j)*(
     &                                      sustr(i,j)
     &            -0.5*(r_D(i-1,j)+r_D(i,j))*DC(i,1) )
        enddo
        do k=2,N,+1          !--> backsubstitution
          do i=istrU,iend
            DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)

            u(i,j,k,nnew)=DC(i,k) * 0.5*(Hz(i,j,k)+Hz(i-1,j,k))

            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
          enddo
        enddo
# ifdef DIAGNOSTICS
        if (diag_uv.and.calc_diag) then
          do k=1,nz
           do i=1,nx
            Udiag(i,j,k,idiss) = ru(i,j,k)*dxdyi_u(i,j)*umask(i,j)
     &          - Udiag(i,j,k,ipgr)
     &          - Udiag(i,j,k,icori)
     &          - Udiag(i,j,k,iadv)
            Udiag(i,j,k,ivmix) = (u(i,j,k,nnew)-u_prev(i,k))/dt -
     &              ru(i,j,k)*dxdyi_u(i,j)
           enddo
          enddo
        endif
# endif

        if (j >= jstrV) then
          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i,j-1,N-1))
     &                        /(  Hz(i,j,N  )+Hz(i,j-1,N  )
     &                           +Hz(i,j,N-1)+Hz(i,j-1,N-1))

            WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i,j-1,N-1))

            cff=1./( 0.5*(Hz(i,j,N)+Hz(i,j-1,N))
     &                     +FC(i,N-1)-min(WC(i,N-1),0.) )

            CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

            DC(i,N)=cff*( v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
     &                                       +dt*svstr(i,j)
     &                                                    )
          enddo
          do k=N-1,2,-1      !--> forward elimination
            do i=istr,iend
              FC(i,k-1)= 2.*dt*(Akv(i,j,k-1)+Akv(i,j-1,k-1))
     &                          /(  Hz(i,j,k  )+Hz(i,j-1,k  )
     &                             +Hz(i,j,k-1)+Hz(i,j-1,k-1))

              WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i,j-1,k-1))

              cff=1./( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                              +FC(i,k-1)-min(WC(i,k-1),0.)
     &                                +FC(i,k)+max(WC(i,k),0.)
     &                       -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                        )
              CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

              DC(i,k)=cff*( v(i,j,k,nnew) +DC(i,0)*rv(i,j,k)
     &                   +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
            enddo
          enddo
          do i=istr,iend
            DC(i,1)=( v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)
     &                       +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                            )/( 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
     &                              +0.5*dt*(r_D(i,j)+r_D(i,j-1))
     &                                  +FC(i,1)+max(WC(i,1),0.)
     &                         -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                          )
            v(i,j,1,nnew)=DC(i,1) * 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
            rvfrc(i,j)=rv(i,j,1) +dm_v(i,j)*dn_v(i,j)*(
     &                                        svstr(i,j)
     &              -0.5*(r_D(i,j-1)+r_D(i,j))*DC(i,1) )
          enddo
          do k=2,N,+1          !--> backsubstitution
            do i=istr,iend
              DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)

              v(i,j,k,nnew)=DC(i,k) * 0.5*(Hz(i,j,k)+Hz(i,j-1,k))
              rvfrc(i,j)=rvfrc(i,j) + rv(i,j,k)
            enddo
          enddo

# ifdef DIAGNOSTICS
         if (diag_uv.and.calc_diag) then
          do k=1,nz
           do i=1,nx
            Vdiag(i,j,k,idiss) = rv(i,j,k)*dxdyi_v(i,j)*vmask(i,j)
     &          - Vdiag(i,j,k,ipgr)
     &          - Vdiag(i,j,k,icori)
     &          - Vdiag(i,j,k,iadv)
            Vdiag(i,j,k,ivmix) = (v(i,j,k,nnew)-v_prev(i,k))/dt -
     &              rv(i,j,k)*dxdyi_v(i,j)
           enddo
          enddo
         endif
# endif

        endif   !<-- j >= jstrV

# else  /* --> not IMPLICIT_BOTTOM_DRAG */

        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          u(i,j,N,nnew)=u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
     &                                     +dt*sustr(i,j) !! add brk2dx
          u(i,j,1,nnew)=u(i,j,1,nnew)+DC(i,0)*ru(i,j,1)

          rufrc(i,j)=ru(i,j,1) +ru(i,j,N) +dn_u(i,j)*dm_u(i,j)*(
     &                                                 sustr(i,j)
     &                 -0.5*(r_D(i-1,j)+r_D(i,j))*u(i,j,1,nstp) )!! add brk2dx
        enddo
        do k=2,N-1
          do i=istrU,iend
            u(i,j,k,nnew)=u(i,j,k,nnew)+DC(i,0)*ru(i,j,k)  ! u*dz
            rufrc(i,j)=rufrc(i,j) +ru(i,j,k)
          enddo
        enddo
# ifdef DIAGNOSTICS
        if (diag_uv.and.calc_diag) then
          do k=1,nz
            Udiag(:,j,k,idiss) = ru(1:nx,j,k)*dxdyi_u(:,j)*umask(i,j)
     &          - Udiag(:,j,k,ipgr)
     &          - Udiag(:,j,k,icori)
     &          - Udiag(:,j,k,iadv)
          enddo
        endif
# endif

        if (j >= jstrV) then
          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            v(i,j,N,nnew)=v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
     &                                       +dt*svstr(i,j)
            v(i,j,1,nnew)=v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)

            rvfrc(i,j)=rv(i,j,1) +rv(i,j,N) +dm_v(i,j)*dn_v(i,j)*(
     &                                                   svstr(i,j) ! add brk2de
     &                   -0.5*(r_D(i,j-1)+r_D(i,j))*v(i,j,1,nstp) )
          enddo
          do k=2,N-1
            do i=istr,iend
              v(i,j,k,nnew)=v(i,j,k,nnew)+DC(i,0)*rv(i,j,k)
              rvfrc(i,j)=rvfrc(i,j) +rv(i,j,k)
            enddo
          enddo
# ifdef DIAGNOSTICS
        if (diag_uv.and.calc_diag) then
          do k=1,nz
            Vdiag(:,j,k,idiss) = rv(1:nx,j,k)*dxdyi_v(:,j)*vmask(i,j)
     &          - Vdiag(:,j,k,ipgr)
     &          - Vdiag(:,j,k,icori)
     &          - Vdiag(:,j,k,iadv)
          enddo
        endif
# endif
        endif  !<-- if j>=jstrV

# endif /* IMPLICIT_BOTTOM_DRAG */
      enddo  !<-- j
      !! velocities are in the form:  u*dz, v*dz, w*dz

# if defined NHMG

#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            u(istr,j,k,nnew) = u(istr,j,k,nnew)*
     &         0.5*(Hz(istr-1,j,k)+Hz(istr,j,k))
          enddo
        enddo
      endif
      if (EASTERN_EDGE) then
        do k=1,N
          do j=jstr,jend
            u(iend+1,j,k,nnew) = u(iend+1,j,k,nnew)*
     &         0.5*(Hz(iend+1,j,k)+Hz(iend,j,k))
          enddo
        enddo
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            v(i,jstr,k,nnew) = v(i,jstr,k,nnew)*
     &         0.5*(Hz(i,jstr-1,k)+Hz(i,jstr,k))
          enddo
        enddo
      endif
      if (NORTHERN_EDGE) then
        do k=1,N
          do i=istr,iend
            v(i,jend+1,k,nnew) = v(i,jend+1,k,nnew)*
     &         0.5*(Hz(i,jend+1,k)+Hz(i,jend,k))
          enddo
        enddo
      endif
#  endif


      ! Translate to fluxes, different than in pre_step, we have dz*u
      ! instead of u

      do k=1,N
        do j=jstr,jend
          do i=istr,iend+1
            u(i,j,k,nnew) = u(i,j,k,nnew)*dn_u(i,j)
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
            v(i,j,k,nnew) = v(i,j,k,nnew)*dm_v(i,j)
          enddo
        enddo
        if (k==N) then
          do j=Jstr,Jend
            do i=istr,iend
              w(i,j,k,nnew) = w(i,j,k,nnew)*dm_r(i,j)*dn_r(i,j)
     &            /(0.5*(            Hz(i,j,k)))
            enddo
          enddo
        else
          do j=Jstr,Jend
            do i=istr,iend
              w(i,j,k,nnew) = w(i,j,k,nnew)*dm_r(i,j)*dn_r(i,j)
     &          /(0.5*(Hz(i,j,k+1)+Hz(i,j,k)))
            enddo
          enddo
        endif
      enddo

#  ifdef EXCHANGE
      !! Todo: we should be able to avoid this. We only need u and v at the boundaries
      !! Todo: computational ranges for u and v need to be extended to Iend+1, Jend+1
      call exchange_xxx(u(:,:,:,nnew),v(:,:,:,nnew) )
#  endif

      !! only do this for neumann nhmg conditions
      ! Compute a 'barotropic' correction to w such that it matches
      ! the current ubar,vbar divergence (see Molemaker et al., 2018)
      if (surface_neumann) then

        nh_ubar = 0.
        nh_vbar = 0.
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend+1
              nh_ubar(i,j) = nh_ubar(i,j) + u(i,j,k,nnew)
            enddo
          enddo
          do j=Jstr,Jend+1
            do i=Istr,Iend
              nh_vbar(i,j) = nh_vbar(i,j) + v(i,j,k,nnew)
            enddo
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            nh_wcor(i,j) = w(i,j,N,nnew) +
     &         (nh_ubar(i+1,j)-nh_ubar(i,j)+nh_vbar(i,j+1)-nh_vbar(i,j))
          enddo
        enddo
        do k=1,N
          do j=Jstr,Jend
            do i=Istr,Iend
              w(i,j,k,nnew) = w(i,j,k,nnew) - nh_wcor(i,j)
     &                       * (z_w(i,j,k)-z_w(i,j,0))
     &                       / (z_w(i,j,N)-z_w(i,j,0))
            enddo
          enddo
        enddo
      endif

      call nhmg_solve(Lm,Mm,N,halo,padd_X,padd_E,
     &           u(:,:,:,nnew),v(:,:,:,nnew),w(:,:,:,nnew) )


      ! At this point u,v,w are fluxes with the dimension of m^3/s
      do k=1,N
        do j=jstr,jend
          do i=Istr,Iend
            u(i,j,k,nnew) = (u(i,j,k,nnew) + mggrid(1)%du(k,j,i) )*pn_u(i,j)

            rufrc(i,j) = rufrc(i,j) + dn_u(i,j)*mggrid(1)%du(k,j,i)/dt
            nhdu(i,j,k,iprec1) = mggrid(1)%du(k,j,i)/dt !stored as a trend
          enddo
        enddo

        do j=jstr,jend
          do i=Istr,Iend
            v(i,j,k,nnew) = (v(i,j,k,nnew) + mggrid(1)%dv(k,j,i) )*pm_v(i,j)

            rvfrc(i,j) = rvfrc(i,j) + dm_v(i,j)*mggrid(1)%dv(k,j,i)/dt
            nhdv(i,j,k,iprec1) = mggrid(1)%dv(k,j,i)/dt !stored as a trend
          enddo
        enddo

        if (k==N) then
          do j=jstr,jend
            do i=Istr,Iend
              w(i,j,k,nnew) = (w(i,j,k,nnew) + mggrid(1)%dw(k+1,j,i) )
     &              * (pm(i,j)*pn(i,j))*(0.5*(            Hz(i,j,k)))

              nhdw(i,j,k,iprec1) = mggrid(1)%dw(k+1,j,i)/dt !stored as a trend
            enddo
          enddo
        else
          do j=jstr,jend
            do i=Istr,Iend
              w(i,j,k,nnew) = (w(i,j,k,nnew) + mggrid(1)%dw(k+1,j,i) )
     &              * (pm(i,j)*pn(i,j))*(0.5*(Hz(i,j,k+1)+Hz(i,j,k)))

              nhdw(i,j,k,iprec1) = mggrid(1)%dw(k+1,j,i)/dt !stored as a trend
            enddo
          enddo
        endif
      enddo

      iprec1 = 1+mod(iprec1,2)
      iprec2 = 1+mod(iprec2,2)

# endif /* NHMG */


       !!!! JM not clear what this does either
# undef r_D
# if defined ADV_ISONEUTRAL && defined EXCHANGE
       !!!! JM not clear what this does !!!
      call exchange_xxx(dRdx,dRde)
# endif

      end

      subroutine check_step_uv1_switches()

! This code must be placed here rather than in a separate file in order
! to be exposed to the relevant CPP-macros defined locally in this file
! above.  It does not affect any model results, other than signature in
! global attribute "CPPS" in output netCDF files.

      use param, only: mynode
      use strings, only: cpps, max_opt_size
      use error_handling_mod, only: error_log
      use utils_mod, only: lenstr

      implicit none
      integer  is,ie
      character(len=24) :: sr_name = "check_step_uv1_switches"
      character(len=1024) :: error_info =""

      write(error_info,*)
     &  'Insufficient length of string "cpps" in file "strings".',
     &        'Increase parameter "max_opt_size" and recompile.'


      ie=lenstr(cpps)
      is=ie+2 ; ie=is+13
      if (ie > max_opt_size) then
         call error_log%raise_global(
     &        info=error_info,
     &        context=sr_name)
      end if
      cpps(is:ie)='<step3d_uv1.F>'

# include "track_advec_switches.h"

      return
      end subroutine check_step_uv1_switches

      subroutine step3d_uv2(tile)
      use param, only: ieast, iwest, jnorth, jsouth, nsub_e, nsub_x
      use private_scratch, only: a2d
      implicit none
      integer tile

# include "compute_tile_bounds.h"
      call step3d_uv2_tile( istr,iend,jstr,jend, A2d(1,1), A2d(1,2),
     &                                           A2d(1,3), A2d(1,4))
      end

      subroutine step3d_uv2_tile(istr,iend,jstr,jend, WC,FC,CF,DC)
      use river_frc, only:
     &     iriver, riv_depth, riv_uvel, riv_vvel,
     &     riv_uflx, riv_vol, riv_vflx, river_source
      use param, only: ieast, iwest, jnorth, jsouth, np_xi, np_eta
      use dimensions, only: inode, jnode

#ifdef ADV_ISONEUTRAL
      use eos_vars, only: drdx, drde, idrz
#endif

#ifdef SPLIT_EOS
      use eos_vars, only: qp2, qp1, rho1
#else
      use eos_vars, only: rho
#endif

      use coupling, only: du_avg1, dv_avg1, du_avg2, dv_avg2
      use grid, only:
     &     nx,ny,nz, dn_u, umask, dm_v, vmask, dm_u, dn_v, f,
     &     dm_r, pn, pm
      use ocean_vars, only: z_r, wi
#ifdef NHMG
     &     , w
#endif
      use mixing, only: akv
#ifdef ADV_ISONEUTRAL
     &     , diff3u, diff3v
#endif
#ifdef LMD_KPP
     &     , hbls
#endif
#ifdef LMD_BKPP
     &     , hbbl
#endif
#ifdef DIAGNOSTICS
      use diagnostics, only: diag_uv, calc_diag, udiag, vdiag, icoup
#endif
      use ocean_vars, only: u, ubar, v, vbar, hz, flxu, flxv, z_w
      use scalars, only: n, knew, nnew, nstp, g, rho0, dt
      use roms_mpi, only: exchange_xxx
#ifdef SOLVE3D
      use u3dbc_mod, only: u3dbc_tile
      use v3dbc_mod, only: v3dbc_tile
#ifdef NHMG
      use w3dbc_mod, only: w3dbc_tile
#endif
#endif

      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real, dimension(PRIVATE_1D_SCRATCH_ARRAY,0:N) :: WC,FC,CF,DC
      real cff !,dmax,div
      real ub,ut
      real :: uint,vint
      real :: huint,hvint

# ifdef ADV_ISONEUTRAL
      real dRz,dRx_max, dpth, r0g, cfs,cfb
      real, parameter :: epsil=1.E-33
      !  Alphamax = 4 leads to mixing that is worse than
      ! up3 over steep slopes
#  define  ALPHA_MAX 2.
#  define GAMMA 0.0833333333333
# endif

# include "compute_auxiliary_bounds.h"

# define IMPLICIT_BOTTOM_DRAG

      do j=jstr,jend

# ifdef DIAGNOSTICS
#  ifndef IMPLICIT_BOTTOM_DRAG
        if (diag_uv) then
          u_prev = u(1:nx,j,:,nnew)
          v_prev = v(1:nx,j,:,nnew)
        endif
#  endif
# endif


# ifndef IMPLICIT_BOTTOM_DRAG

! Apply implicit time step for vertical viscosity term to horizontal
! velocity components,
!
!   du(k)     1     [         u(k+1)-u(k)             u(k)-u(k-1) ]
!  ------ = ----- * [ Akv(k)* ----------- - Akv(k-1)* ----------- ]
!    dt     Hz(k)   [           dz(k)                   dz(k-1)   ]
!
! where u(k) represents u,v(:,:,:,nnew) velocity components; Hz(k) and
! dz(k) are grid spacing between vertical W and RHO-points respectively
! interpolated to horizontal U,V locations. Backward Euler implicit
! time step requires solution of a tri-diagonal system,
!
!     -FC(k-1)*u_new(k-1) +[Hz(k)+FC(k-1)+FC(k)]*u_new(k)
!                       -FC(k)*u_new(k+1) = Hz(k)*u_old(k)
!
!                dt*Akv(k)
! where FC(k) = ----------- is normalized vertical viscosity defined
!                  dz(k)
!
! at W-points; u_new(k) is the new-time-step (unknown) velocity;
! u_old(k) is old-time-step tracer (known).  As long as vertical
! viscosity Ak(k) is nonnegative, the tri-diagonal matrix is
! diagonally dominant which guarantees stability of Gaussian
! elimination procedure, (e.g., Richtmeyer annd  Morton,  1967).
! Since top (wind) and bottom (drag) boundary stresses are applied
! explicitly, for the purpose of solving the tri-diagonal system
! boundary conditions are effectively no-stress Akv(N)=Akt(0)=0, so
! that the corresponding FC(N)=FC(1)=0. This leads to equations for
! top and bottom grid boxes,
!
!   -FC(N-1)*u_new(N-1) +[Hz(N)+FC(N-1)]*u_new(N) = Hz(N)*u_old(N)
!
!          [Hz(1)+FC(1)]*u_new(1) -FC(1)*u_new(2) = Hz(1)*u_old(1)
!
! Basically FC(N)=FC(1)=0 conditions should be interpreted as the
! implicit step merely redistributes the horizontal momentum
! components u,v throughout the vertical column, while conserving
! its vertical integral. At this moment the content of arrays
! u,v(:,:,:,nnew) has meaning of Hz*velocity. After the implicit
! step they becomes just velocity, which is temporarily loaded into
! array DC(:,1:N)
!>
! Also perform coupling between the barotropic and baroclinic modes.
! For the purpose of the computational efficiency (cache utilization)
! the implicit time step and coupling are partially overlapped. To
! distinguish the two operations in the code below we place separator
! (!>): lines above the separator (!>) belong to the tri-diagonal
! solver, while operations below it belong to the coupling.
! The coupling procedure makes vertical integration of the newly
! computed velocities [which ends up being stored in the scratch
! array DC(:,1:N) after the tri-diagonal problem is resolved] to
! compute its vertical mean [stored as DC(:,0)], and replace it with
! that from the 2D barotropic sub-model ubar=DU_avg1/(D*dn_u), where
! D=CF(:,0) is the total depth of water column.
!
!
! Time step momentum equation in the XI-direction
!----- ---- -------- -------- -- --- ------------

        do i=istrU,iend
          DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

          FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i-1,j,N-1))
     &                      /(  Hz(i,j,N  )+Hz(i-1,j,N  )
     &                         +Hz(i,j,N-1)+Hz(i-1,j,N-1))

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i-1,j,N-1))

          cff=1./( 0.5*(Hz(i,j,N)+Hz(i-1,j,N))
     &                   +FC(i,N-1)-min(WC(i,N-1),0.) )

          CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

          DC(i,N)=cff*u(i,j,N,nnew)
        enddo
        do k=N-1,2,-1      !--> forward elimination
          do i=istrU,iend
            FC(i,k-1)= 2.*dt*(Akv(i,j,k-1)+Akv(i-1,j,k-1))
     &                        /(  Hz(i,j,k  )+Hz(i-1,j,k  )
     &                           +Hz(i,j,k-1)+Hz(i-1,j,k-1))

            WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i-1,j,k-1))

            cff=1./( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))
     &                           +FC(i,k-1)-min(WC(i,k-1),0.)
     &                             +FC(i,k)+max(WC(i,k),0.)
     &                    -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                      )
            CF(i,k-1)=cff*(   FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=cff*( u(i,j,k,nnew)
     &                 +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
          enddo
        enddo       !--> discard WC,DC(i,0)
        do i=istrU,iend
          DC(i,1)=( u(i,j,1,nnew) +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                                )/( 0.5*(Hz(i,j,1)+Hz(i-1,j,1))
#  ifdef IMPLCT_NO_SLIP_BTTM_BC
     &                                  +0.5*dt*(r_D(i,j)+r_D(i-1,j))
#  endif
     &                                      +FC(i,1)+max(WC(i,1),0.)
     &                             -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                       )
!>
          CF(i,0)=0.5*(Hz(i,j,1)+Hz(i-1,j,1))    !--> coupling
          DC(i,0)=CF(i,0)*DC(i,1)
        enddo
        do k=2,N,+1          !--> backsubstitution
          do i=istrU,iend
            DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1) !<-- backsubstitution
!>
            cff=0.5*(Hz(i,j,k)+Hz(i-1,j,k))      !--> coupling
            CF(i,0)=CF(i,0)+cff                  !    [ vertical    ]
            DC(i,0)=DC(i,0)+cff*DC(i,k)          !    [ integration ]
          enddo
        enddo
        do i=istrU,iend
          DC(i,0)=( DC(i,0)*dn_u(i,j)-DU_avg1(i,j)  !<-- mismatch
     &                       )/(CF(i,0)*dn_u(i,j))
        enddo
        do k=N,1,-1
          do i=istrU,iend
            u(i,j,k,nnew)=(DC(i,k)-DC(i,0))      !<-- correction
#  ifdef MASKING
     &                          *umask(i,j)
#  endif
          enddo
        enddo        !--> discard DC,CF

#  ifdef DIAGNOSTICS
      if (diag_uv .and. calc_diag) then
       do k = 1,nz
        do i = 1,nx
          Udiag(i,j,1:nz,ivmix) = ( DC(i,k)*
     &          0.5*( Hz(i,j-1,k)+Hz(i,j,k) )-v_prev(i,k) ) / dt
          Udiag(i,j,k,icoup) = -DC(i,0)*umask(i,j)
     &           *0.5*(Hz(i,j,k)+Hz(i-1,j,k))/dt
        enddo
       enddo
      endif
#  endif

# else /* yes IMPLICIT_BOTTOM_DRAG */
        ! The only thing we do is 'coupling'.
        ! Correcting the 3D flow to be in line with the 2D prediction

        do i=istrU,iend                          !--> coupling only
          CF(i,0)=0.5*(Hz(i,j,N)+Hz(i-1,j,N))
          DC(i,0)=u(i,j,N,nnew)
          u(i,j,N,nnew)=u(i,j,N,nnew)/CF(i,0)
        enddo
        do k=N-1,1,-1
          do i=istrU,iend
            cff=0.5*(Hz(i,j,k)+Hz(i-1,j,k))
            CF(i,0)=CF(i,0) +cff
            DC(i,0)=DC(i,0) +u(i,j,k,nnew)
            u(i,j,k,nnew)=u(i,j,k,nnew)/cff
          enddo
        enddo
        do i=istrU,iend
          DC(i,0)=( DC(i,0)*dn_u(i,j)-DU_avg1(i,j)  !<-- mismatch
     &                       )/(CF(i,0)*dn_u(i,j))
        enddo
        do k=1,nz
          do i=istrU,iend
            u(i,j,k,nnew)=(u(i,j,k,nnew)-DC(i,0))
#  ifdef MASKING
     &                                *umask(i,j)
#  endif
          enddo
        enddo
#  ifdef DIAGNOSTICS
        if (diag_uv .and. calc_diag) then
         do k=1,nz
          do i=1,nx
           Udiag(i,j,k,icoup) = -DC(i,0)*umask(i,j)
     &           *0.5*(Hz(i,j,k)+Hz(i-1,j,k))/dt
          enddo
         enddo
        endif
#  endif

# endif

        if (j>=jstrV) then

# ifndef IMPLICIT_BOTTOM_DRAG

! Time step momentum equation in the ETA-direction.
!----- ---- -------- -------- -- --- --------------
! (All steps are similar to XI-component above.)

          do i=istr,iend
            DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

            FC(i,N-1)= 2.*dt*(Akv(i,j,N-1)+Akv(i,j-1,N-1))
     &                        /(  Hz(i,j,N  )+Hz(i,j-1,N  )
     &                           +Hz(i,j,N-1)+Hz(i,j-1,N-1))

            WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N-1)+Wi(i,j-1,N-1))

            cff=1./( 0.5*(Hz(i,j,N)+Hz(i,j-1,N))
     &                     +FC(i,N-1)-min(WC(i,N-1),0.) )

            CF(i,N-1)=cff*( FC(i,N-1)+max(WC(i,N-1),0.) )

            DC(i,N)=cff*v(i,j,N,nnew)
          enddo
          do k=N-1,2,-1      !--> forward elimination
            do i=istr,iend
              FC(i,k-1)= 2.*dt*(Akv(i,j,k-1)+Akv(i,j-1,k-1))
     &                          /(  Hz(i,j,k  )+Hz(i,j-1,k  )
     &                             +Hz(i,j,k-1)+Hz(i,j-1,k-1))

              WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k-1)+Wi(i,j-1,k-1))

              cff=1./( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))
     &                              +FC(i,k-1)-min(WC(i,k-1),0.)
     &                                +FC(i,k)+max(WC(i,k),0.)
     &                       -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                        )
              CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

              DC(i,k)=cff*( v(i,j,k,nnew)
     &                   +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)) )
            enddo
          enddo
          do i=istr,iend
            DC(i,1)=( v(i,j,1,nnew) +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                                  )/( 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
     &                                    +0.5*dt*(r_D(i,j)+r_D(i,j-1))
     &                                        +FC(i,1)+max(WC(i,1),0.)
     &                               -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                              )
!>
            CF(i,0)=0.5*(Hz(i,j,1)+Hz(i,j-1,1))   !--> coupling
            DC(i,0)=CF(i,0)*DC(i,1)
          enddo
          do k=2,N,+1          !--> backsubstitution
            do i=istr,iend
              DC(i,k)=DC(i,k)+CF(i,k-1)*DC(i,k-1) !<-- backsubstitution
!>
              cff=0.5*(Hz(i,j,k)+Hz(i,j-1,k))     !--> coupling
              CF(i,0)=CF(i,0)+cff                 !    [ vertical    ]
              DC(i,0)=DC(i,0)+cff*DC(i,k)         !    [ integration ]
            enddo
          enddo
          do i=istr,iend
            DC(i,0)=( DC(i,0)*dm_v(i,j)-DV_avg1(i,j)  !<-- mismatch
     &                         )/(CF(i,0)*dm_v(i,j))
          enddo
          do k=N,1,-1
            do i=istr,iend
              v(i,j,k,nnew)=(DC(i,k)-DC(i,0))     !<-- correction
#  ifdef MASKING
     &                            *vmask(i,j)
#  endif
            enddo
          enddo       !--> discard DC,CF(:,0)

#  ifdef DIAGNOSTICS
      if (diag_uv .and. calc_diag) then
       do k=1,nz
        do i=1,nx
          Vdiag(i,j,k,ivmix) = ( DC(i,k)*
     &          0.5*( Hz(i,j-1,k)+Hz(i,j,k) )-v_prev(i,k) ) / dt
          Vdiag(i,j,k,icoup) = -DC(i,0)*vmask(i,j)
     &           *0.5*(Hz(i,j,k)+Hz(i,j-1,k))/dt
        enddo
       enddo
      endif
#  endif


# else  /* yes IMPLICIT_BOTTOM_DRAG */

          do i=istr,iend
            CF(i,0)=0.5*(Hz(i,j,N)+Hz(i,j-1,N))
            DC(i,0)=v(i,j,N,nnew)
            v(i,j,N,nnew)=v(i,j,N,nnew)/CF(i,0)
          enddo
          do k=N-1,1,-1
            do i=istr,iend
              cff=0.5*(Hz(i,j,k)+Hz(i,j-1,k))
              CF(i,0)=CF(i,0) +cff
              DC(i,0)=DC(i,0) +v(i,j,k,nnew)
              v(i,j,k,nnew)=v(i,j,k,nnew)/cff
            enddo
          enddo      !--> discard FC, keep DC,CF(:,0)
          do i=istr,iend
            DC(i,0)=( DC(i,0)*dm_v(i,j)-DV_avg1(i,j)
     &                         )/(CF(i,0)*dm_v(i,j))
          enddo
          do k=1,N
            do i=istr,iend
              v(i,j,k,nnew)=(v(i,j,k,nnew)-DC(i,0))
#  ifdef MASKING
     &                                  *vmask(i,j)
#  endif
            enddo
          enddo       !--> discard DC,CF(:,0)
# endif
#  ifdef DIAGNOSTICS
      if (diag_uv .and. calc_diag) then
       do k=1,nz
        do i=1,nx
          Vdiag(i,j,k,icoup) = -DC(i,0)*vmask(i,j)
     &           *0.5*(Hz(i,j,k)+Hz(i,j-1,k))/dt
        enddo
       enddo
      endif
#  endif

        endif  !<-- j>=jstrV

# ifdef NHMG
      !!! JM fix this wrt to IMPLICIT BOTTOM STRESS
        do i=istr,iend
          DC(i,0) =dt*pm(i,j)*pn(i,j)
          FC(i,N-1)= 0.5*dt*(Akv(i,j,N)+Akv(i,j,N-1))
     &                      /Hz(i,j,N)

          WC(i,N-1)= DC(i,0)*0.5*(Wi(i,j,N)+Wi(i,j,N-1))

          cff=1./( 0.5*(Hz(i,j,N)                )
     &                    +FC(i,N-1)-min(WC(i,N-1),0.) )

          CF(i,N-1)=cff*(  FC(i,N-1)+max(WC(i,N-1),0.) )

          DC(i,N)=cff*w(i,j,N,nnew)
        enddo

        do k=N-1,2,-1      !--> forward elimination for w
          do i=istr,iend
            FC(i,k-1)= 0.5*dt*(Akv(i,j,k)+Akv(i,j,k-1))!! flux coef between w(k) and w(k-1)
     &                   /Hz(i,j,k)

            WC(i,k-1)= DC(i,0)*0.5*(Wi(i,j,k)+Wi(i,j,k-1))

            cff=1./( 0.5*(Hz(i,j,k)+Hz(i,j,k-1))
     &                            +FC(i,k-1)-min(WC(i,k-1),0.)
     &                            +FC(i,k  )+max(WC(i,k),0.)
     &                     -CF(i,k)*(FC(i,k)-min(WC(i,k),0.))
     &                                                       )
            CF(i,k-1)=cff*( FC(i,k-1)+max(WC(i,k-1),0.) )

            DC(i,k)=cff*( w(i,j,k,nnew)
     &                 +DC(i,k+1)*(FC(i,k)-min(WC(i,k),0.)))
          enddo
        enddo
        !! change this you feel energetic. Instead of no-flux at the bottom we
        !! have w(0) = 0 (different BC for tri-diag system)
        do i=istr,iend
          w(i,j,1,nnew)=( w(i,j,k,nnew)  +DC(i,2)*(FC(i,1)-min(WC(i,1),0.))
     &                         )/( 0.5*(Hz(i,j,2)+Hz(i,j,1))
     &                                       +FC(i,1)+max(WC(i,1),0.)
     &                              -CF(i,1)*(FC(i,1)-min(WC(i,1),0.))
     &                                                               )
        enddo
        do k=2,N,+1          !--> backsubstitution for w
          do i=istr,iend
            w(i,j,k,nnew)=DC(i,k) +CF(i,k-1)*w(i,j,k-1,nnew)
          enddo
        enddo
        !------- end computing of w(:,:,:,nnew) -----

# endif /* NHMG */

      enddo      ! <-- j

! Set PHYSICAL lateral boundary conditions.

      call u3dbc_tile (istr,iend,jstr,jend, DC)
      call v3dbc_tile (istr,iend,jstr,jend, DC)


! Coupling 2D and 3D momentum equations:
!--------- -- --- -- -------- ----------
! Compute inaccurate vertical integral of 3D velocity field, and
! replace it with the one computed by barotropic mode. After that
! compute mass fluxes through grid box faces.
!
! Meaning of scratch variables below:
!
! DC(i,k) [where k=1,N]  height of grid box for U- or V-cell.
! DC(i,0) total depth of water column at horizontal U- or V-points.
! FC(i,0) vertically integrated mass flux, error, and correction for
!                       the new time step velocity u,v(:,:,:,nnew)
! CF(i,0) vertically integrated mass flux/error/correction for the
!             intermediate time step [n+1/2] mass fluxes FlxU,FlxV.
!
! This procedure also replaces instantaneous (with respect to the
! fast time step) vertically integrated (barotropic) velocities with
! their values based on fast-time-averaged mass fluxes. These are to
! be used as initial conditions for the barotropic mode at the new
! time step.
!
! Explanation of horizontal loop indices: in the case of periodic
! boundaries (in either direction) the coupling and computation of
! mass fluxes FlxU,FlxV is performed within the internal range of
! indices (excluding ghost zones), after that the ghost points for
! the newly computed arrays are exchanged; in the case of nonperiodic
! boundaries the coupling is done over the extended range of indices
! (that is including boundary points).

# ifdef EW_PERIODIC
#  define IU_RANGE istr,iend
#  define IV_RANGE istr,iend
# else
#  define IU_RANGE istr,iendR
#  define IV_RANGE istrR,iendR
# endif

# ifdef NS_PERIODIC
#  define J_RANGE jstr,jend
# else
#  define J_RANGE jstrR,jendR
# endif
      ! DevinD & JM: seems un-necessary to 2D-3D couple the interior when only
      ! the boundary has been update since the first 2D-3D correction above.
      ! We only add the boundary changes to the diagnostics and our budget still closes,
      ! as the interior correction is only in the order of e-17.
      do j=J_RANGE                  !--> XI-component
        do i=IU_RANGE
          DC(i,N)=0.5*(Hz(i,j,N)+Hz(i-1,j,N))*dn_u(i,j)
          DC(i,0)=DC(i,N)
          FC(i,0)=DC(i,N)*u(i,j,N,nnew)
        enddo
        do k=N-1,1,-1
          do i=IU_RANGE
            DC(i,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            FC(i,0)=FC(i,0)+DC(i,k)*u(i,j,k,nnew)
          enddo
        enddo
        ! FC(0) is int(dz*U)*dy, DC(0) is int(dz)*dy
        do i=IU_RANGE
          DC(i,0)=1./DC(i,0)  !!1/H
          ubar(i,j,knew)=DC(i,0)*DU_avg1(i,j)  !! avg1 is at n+1?

          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg1(i,j)) !! FC(0) is now correction velocity
!>
          CF(i,0)=0.   !<-- initialize vertical integral
        enddo          !--> keep DC(:,0) for later

        do k=N,1,-1
          do i=IU_RANGE
            u(i,j,k,nnew)=(u(i,j,k,nnew)-FC(i,0))
# ifdef MASKING
     &                                  *umask(i,j)
# endif
c--#define EPSIL 0.125
c--#define DELTA 0.75

#  define EPSIL 0.36
#  define DELTA 0.28

c>># define EPSIL 0.2105
c>># define DELTA 0.5790

            CF(i,k)=DELTA*FlxU(i,j,k) + EPSIL*DC(i,k)
     &                         *(u(i,j,k,nstp)+u(i,j,k,nnew))
            CF(i,0)=CF(i,0)+CF(i,k)    !<-- vertical integration
          enddo
        enddo
        do i=IU_RANGE
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg2(i,j))  !<-- mismatch
        enddo

        do k=1,N,+1
          do i=IU_RANGE
            FlxU(i,j,k)=CF(i,k)-DC(i,k)*CF(i,0)   !<-- correction

# ifdef ADV_ISONEUTRAL
            diff3u(i,j,k)=sqrt(GAMMA*abs(u(i,j,k,nnew))*dm_u(i,j))
     &                                                  *dm_u(i,j)
# endif
          enddo
        enddo

        if (j >= jstr) then         !--> ETA-component
          do i=IV_RANGE
            DC(i,N)=0.5*(Hz(i,j,N)+Hz(i,j-1,N))*dm_v(i,j)
            DC(i,0)=DC(i,N)
            FC(i,0)=DC(i,N)*v(i,j,N,nnew)
          enddo
          do k=N-1,1,-1
            do i=IV_RANGE
              DC(i,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              FC(i,0)=FC(i,0)+DC(i,k)*v(i,j,k,nnew)
            enddo
          enddo
          do i=IV_RANGE
            DC(i,0)=1./DC(i,0)
            vbar(i,j,knew)=DC(i,0)*DV_avg1(i,j)

            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg1(i,j))
            CF(i,0)=0.
          enddo        !--> keep DC(:,0) for later
          do k=N,1,-1
            do i=IV_RANGE
              v(i,j,k,nnew)=(v(i,j,k,nnew)-FC(i,0))
# ifdef MASKING
     &                                    *vmask(i,j)

              CF(i,k)=DELTA*FlxV(i,j,k) + EPSIL*DC(i,k)
     &                           *(v(i,j,k,nstp)+v(i,j,k,nnew))
              CF(i,0)=CF(i,0)+CF(i,k)
            enddo
          enddo
          do i=IV_RANGE
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg2(i,j))
          enddo
          do k=1,N,+1
            do i=IV_RANGE
              FlxV(i,j,k)=CF(i,k)-DC(i,k)*CF(i,0)
# ifdef ADV_ISONEUTRAL
              diff3v(i,j,k)=sqrt(GAMMA*abs(v(i,j,k,nnew))*dn_v(i,j))
     &                                                   *dn_v(i,j)
# endif
            enddo
          enddo
# ifdef ADV_ISONEUTRAL
          if (j <= jend) then
            r0g=rho0/g
            do k=1,N-1
              do i=istr,iend
#  define ISO
#  ifdef ISO
#    ifdef SPLIT_EOS
                dpth=-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
                dRz=rho1(i,j,k)-rho1(i,j,k+1)+(qp1(i,j,k)-qp1(i,j,k+1))
     &                                      *dpth*(1.-2.*qp2*dpth)
#    else
                dRz=rho(i,j,k)-rho(i,j,k+1)
#    endif
                dRz=max(dRz,0.) +r0g*f(i,j)**2
     &                        *(z_r(i,j,k+1)-z_r(i,j,k))
#  else
                dRz=z_r(i,j,k+1)-z_r(i,j,k)
#  endif


c*              dRx_max=(z_r(i,j,k+1)-z_r(i,j,k))*max(
c*   &                        abs(dRdx(i,j,k  )), abs(dRdx(i+1,j,k  )),
c*   &                        abs(dRdx(i,j,k+1)), abs(dRdx(i+1,j,k+1)),
c*
c*   &                        abs(dRde(i,j,k  )), abs(dRde(i,j+1,k  )),
c*   &                        abs(dRde(i,j,k+1)), abs(dRde(i,j+1,k+1))
c*   &                                                               )

               dRx_max=max(

     &       dm_u(i  ,j)*max(abs(dRdx(i  ,j,k)), abs(dRdx(i  ,j,k+1))),
     &       dm_u(i+1,j)*max(abs(dRdx(i+1,j,k)), abs(dRdx(i+1,j,k+1))),

     &       dn_v(i,j  )*max(abs(dRde(i,j  ,k)), abs(dRde(i,j  ,k+1))),
     &       dn_v(i,j+1)*max(abs(dRde(i,j+1,k)), abs(dRde(i,j+1,k+1)))
     &                                                               )

# ifdef LMD_KPP
                cfs=min(1., (z_w(i,j,N)-z_w(i,j,k))/max( 50.,
     &                                    hbls(i,j) ))
# else
                cfs=min(1., (z_w(i,j,N)-z_w(i,j,k))/50.)
# endif
# ifdef LMD_BKPP
                cfb=min(1., (z_w(i,j,k)-z_w(i,j,0))/max( 50.,
     &                                    hbbl(i,j) ))
# else
                cfb=min(1., (z_w(i,j,k)-z_w(i,j,0))/50.)
# endif

                cff=ALPHA_MAX*cfs*(2.-cfs) * cfb*(2.-cfb)


c               idRz(i,j,k)=1./dRz

               idRz(i,j,k)=cff/max(cff*dRz, dRx_max, epsil)

c             idRz(i,j,k)=ALPHA_MAX/max(ALPHA_MAX*dRz,dRx_max,epsil)

              enddo
            enddo
          endif      !<-- j <= jend
# endif
        endif    !<-- j >= jstr
      enddo    ! <-- j

      if (river_source) then
        do j=jstr,jend
          do i=istr,iend

            if (abs(riv_uflx(i,j)).gt.1e-3) then
              riv_depth = 0.5*( z_w(i-1,j,N)-z_w(i-1,j,0)
     &                     +z_w(i  ,j,N)-z_w(i  ,j,0) )
              iriver = nint(riv_uflx(i,j)/10)
              riv_uvel = riv_vol(iriver)*(riv_uflx(i,j)-10*iriver)/
     &                   ( dn_u(i,j)*riv_depth)
              do k= 1,N
                u(i,j,k,nnew) = riv_uvel
              enddo
            endif

            if (abs(riv_vflx(i,j)).gt.1e-3) then
              riv_depth = 0.5*( z_w(i,j-1,N)-z_w(i,j-1,0)
     &                     +z_w(i,j  ,N)-z_w(i,j  ,0) )
              iriver = nint(riv_vflx(i,j)/10)
              riv_vvel = riv_vol(iriver)*(riv_vflx(i,j)-10*iriver)/
     &                   ( dm_v(i,j)*riv_depth)
              do k= 1,N
                v(i,j,k,nnew) = riv_vvel
              enddo
            endif

          enddo
        enddo
      endif ! <-- river_source

# undef IU_RANGE
# undef IV_RANGE
# undef J_RANGE

! Exchange periodic boundaries and computational margins.

# ifdef EXCHANGE
      call exchange_xxx(FlxU,u(:,:,:,nnew),ubar(:,:,knew) )

      call exchange_xxx(FlxV,v(:,:,:,nnew),vbar(:,:,knew) )

#  ifdef ADV_ISONEUTRAL
      call exchange_xxx(idRz,diff3u,diff3v)
#  endif
# endif


# ifdef NHMG
!--------------------------------------------------------------------
! Coupling w
!--------------------------------------------------------------------
!
      ! We've modified u and v, w needs to be made consistent
!      do j=istr,iend  ! BUG HERE
      do j=jstr,jend  ! BUG HERE
        do k=1,N
          do i=Istr,Iend
            DC(i,0)=pm(i,j)*pn(i,j)
            ! at this time, Huon and Hvom CANNOT used : they are
            ! defined at n+1/2 instead of n+1
            w(i,j,k,nnew)=w(i,j,k-1,nnew)
     &      - DC(i,0)*(
     &       +0.5*(Hz(i,j,k)+Hz(i+1,j,k))*dn_u(i+1,j)*u(i+1,j,k,nnew)
     &       -0.5*(Hz(i-1,j,k)+Hz(i,j,k))*dn_u(i,j)*u(i,j,k,nnew)
     &       +0.5*(Hz(i,j,k)+Hz(i,j+1,k))*dm_v(i,j+1)*v(i,j+1,k,nnew)
     &       -0.5*(Hz(i,j-1,k)+Hz(i,j,k))*dm_v(i,j)*v(i,j,k,nnew) )
          enddo
        enddo
      enddo
# ifdef DIAGNOSTICS_NHMG
      if (calc_diag) call set_diags_w_at_uv2( istr, iend, jstr, jend )
# endif
!

!--------------------------------------------------------------------
! Update physical boundary values
!--------------------------------------------------------------------
#ifdef SOLVE3D
      call w3dbc_tile (Istr,Iend,Jstr,Jend, DC)
#endif

#  ifdef DIAGNOSTICS_NHMG
      if (calc_diag) call set_diags_w_at_bc( istr, iend, jstr, jend )
#  endif
!--------------------------------------------------------------------
! Exchange periodic boundaries and computational margins.
!--------------------------------------------------------------------
#  ifdef EXCHANGE
      call exchange_xxx(w(:,:,:,nnew))
#  endif

# endif /* NHMG */

# ifdef DIAGNOSTICS
#  ifdef NHMGDIAG
      if (calc_diag.and.diag_uv)  call set_diags_w_at_uv2_end
#  endif
# endif

      end


      subroutine check_step_uv2_switches()

! This code must be placed here rather than in a separate file in order
! to be exposed to the relevant CPP-macros defined locally in this file
! above.  It does not affect any model results, other than signature in
! global attribute "CPPS" in output netCDF files.

      use param, only: mynode
      use strings, only: cpps, max_opt_size
      use error_handling_mod, only: error_log
      use utils_mod, only: lenstr

      implicit none
      integer is,ie

      integer lstr
      character(len=32) tmpstr
      character(len=24) :: sr_name = "check_step_uv2_switches"
      character(len=1024) :: error_info =""

      write(error_info,*)
     &  'Insufficient length of string "cpps" in file "strings".',
     &        'Increase parameter "max_opt_size" and recompile.'

      ie=lenstr(cpps)
      is=ie+2 ; ie=is+13
      if (ie > max_opt_size) then
         call error_log%raise_global(
     &        info=error_info,
     &        context=sr_name)
      end if

      cpps(is:ie)='<step3d_uv2.F>'
# ifdef DELTA
      write(tmpstr,*) DELTA ;  lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr)=='0' .and. lstr>3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+5+lstr
      if (ie > max_opt_size) then
         call error_log%raise_global(
     &        info=error_info,
     &        context=sr_name)
      end if

      cpps(is:ie)='DELTA='/ /tmpstr(1:lstr)
# endif
# ifdef EPSIL
      write(tmpstr,*) EPSIL  ;  lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr)=='0' .and. lstr>3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+5+lstr
      if (ie > max_opt_size) then
         call error_log%raise_global(
     &        info=error_info,
     &        context=sr_name)
      end if

      cpps(is:ie)='EPSIL='/ /tmpstr(1:lstr)
# endif
# ifdef GAMMA
      write(tmpstr,*) GAMMA ; lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr)=='0' .and. lstr>3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+5+lstr
      if (ie > max_opt_size) then
         call error_log%raise_global(
     &        info=error_info,
     &        context=sr_name)
      end if

      cpps(is:ie)='GAMMA='/ /tmpstr(1:lstr)
# endif
# ifdef ALPHA_MAX
      write(tmpstr,'(F16.12)') ALPHA_MAX
      lstr=lenstr(tmpstr)
      do while(tmpstr(lstr:lstr) == '0' .and. lstr > 3)
        lstr=lstr-1
      enddo
      is=ie+2 ; ie=is+9+lstr
      if (ie > max_opt_size) then
         call error_log%raise_global(
     &        info=error_info,
     &        context=sr_name)
      end if

      cpps(is:ie)='ALPHA_MAX='/ /tmpstr(1:lstr)
# endif
      return
      end subroutine check_step_uv2_switches

#else
      contains
      subroutine step3d_uv_empty
      end
#endif
      end module step3d_uv_mod
