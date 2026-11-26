! Forward sweep of the implicit vertical solve of tracer itrc on row j
! (step3d_t_mod.F90), included twice there: with ISO_AKZ defined the
! diffusivity adds the isoneutral Akz (adv_isoneutral), without it the
! loops are those of a build without isoneutral mixing.

    do i=istr,iend
      DC(i,0)=dt*pm(i,j)*pn(i,j) !<--horizontal metric

# ifdef ISO_AKZ
      FC(i,1)=2._8*dt*(Akt(i,j,1,iAkt)+Akz(i,j,1))&
# else
      FC(i,1)=2._8*dt* Akt(i,j,1,iAkt)&
# endif
      &/(Hz(i,j,1)+Hz(i,j,2))
      WC(i,1)=DC(i,0)*Wi(i,j,1)

      cff=1._8/(Hz(i,j,1) +FC(i,1)+max(WC(i,1),0._8))
      CF(i,1)=cff*(      FC(i,1)-min(WC(i,1),0._8))
      DC(i,1)=cff*t(i,j,1,nnew,itrc)
      if (itrc == iTandS) then
        td_piv(i,1)=cff
        td_cf (i,1)=CF(i,1)
        td_a  (i,1)=FC(i,1)+max(WC(i,1),0._8)
      endif
    enddo
    do k=2,nz-1,+1
      do i=istr,iend
# ifdef ISO_AKZ
        FC(i,k)=2._8*dt*(Akt(i,j,k,iAkt)+Akz(i,j,k))&
# else
        FC(i,k)=2._8*dt* Akt(i,j,k,iAkt)&
# endif
        &/(Hz(i,j,k)+Hz(i,j,k+1))

        WC(i,k)=Wi(i,j,k)*DC(i,0)

        cff=1._8/( Hz(i,j,k) +FC(i,k)+max(WC(i,k),0._8)&
        &+FC(i,k-1)-min(WC(i,k-1),0._8)&
        &-CF(i,k-1)*(FC(i,k-1)+max(WC(i,k-1),0._8))&
        &)

        CF(i,k)=cff*(       FC(i,k)-min(WC(i,k),0._8) )

        DC(i,k)=cff*( t(i,j,k,nnew,itrc) +DC(i,k-1)*(&
        &FC(i,k-1)+max(WC(i,k-1),0._8) ))
        if (itrc == iTandS) then
          td_piv(i,k)=cff
          td_cf (i,k)=CF(i,k)
          td_a  (i,k)=FC(i,k)+max(WC(i,k),0._8)
        endif
      enddo
    enddo          !--> discard DC(:,0)
