! Flux components of the horizontal divergence of the stress tensor at
! level k (visc3d_mod.F90, visc3d_s_tile), included twice there: UU, VV
! are the velocities the Laplacian acts on, the stock u, v(nstp), or (imex_igw
! with imex_mix_split 1, 2) the split-step fields Uk, Vk.
      do j=jstrV-1,jend
        do i=istrU-1,iend
          cff=0.5_8*Hz(i,j,k)*visc2_r(i,j)*(&

          &dn_r(i,j)*pm(i,j)*( (pn(i  ,j)+pn(i+1,j))*UU(i+1,j)&
          &-(pn(i-1,j)+pn(i  ,j))*UU(i  ,j)&
          &)&

          &-dm_r(i,j)*pn(i,j)*( (pm(i,j  )+pm(i,j+1))*VV(i,j+1)&
          &-(pm(i,j-1)+pm(i,j  ))*VV(i,j  )&
          &))
          UFx(i,j)=  cff*dn_r(i,j)*dn_r(i,j)
          VFe(i,j)= -cff*dm_r(i,j)*dm_r(i,j)
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend+1
          cff=0.125_8*(Hz(i-1,j,k)+Hz(i,j,k)+Hz(i-1,j-1,k)&
          &+Hz(i,j-1,k))*visc2_p(i,j)*(&

          &0.25_8*(pm(i-1,j)+pm(i,j)+pm(i-1,j-1)+pm(i,j-1))*dn_p(i,j)&
          &*( (pn(i  ,j-1)+pn(i  ,j))*VV(i  ,j)&
          &-(pn(i-1,j-1)+pn(i-1,j))*VV(i-1,j)&
          &)&

          &+0.25_8*(pn(i-1,j)+pn(i,j)+pn(i-1,j-1)+pn(i,j-1))*dm_p(i,j)&
          &*( (pm(i-1,j  )+pm(i,j  ))*UU(i,j  )&
          &-(pm(i-1,j-1)+pm(i,j-1))*UU(i,j-1)&
          &))&
# ifdef MASKING
          &*pmask(i,j)
# endif
          UFe(i,j)= cff*dm_p(i,j)*dm_p(i,j)
          VFx(i,j)= cff*dn_p(i,j)*dn_p(i,j)
        enddo
      enddo
