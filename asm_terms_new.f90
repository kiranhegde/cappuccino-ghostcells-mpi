!***********************************************************************
!
      subroutine additional_algebraic_stress_terms
!
!***********************************************************************
!     Calculates the additional agebraic stress terms for momentum eq.
!     and ads them to su, sv, sw rhs. vectors.
!
!***********************************************************************
      use types
      use parameters
      use indexes
      use geometry
      use coef
      use coefb
      use variables
      use gradients

      implicit none
!
!***********************************************************************
!
!
!     Local variables
!
      integer :: i, j, k, inp,ine,inw,inn,ins,int,inb
      real(prec) :: two_thirds, &
                    term1e, term1w, term1n, term1s, term1t, term1b, dterm1dx, &   !
                    term2e, term2w, term2n, term2s, term2t, term2b, dterm2dy, &   !
                    term3e, term3w, term3n, term3s, term3t, term3b, dterm3dz      ! end vars. for afm


!     the additional algebraic stress terms
!
      two_thirds=2./3.0_dp
      ine = inp+nj
      inw = inp-nj
      inn = inp+1
      ins = inp-1
      int = inp+nij
      inb = inp-nij
!-----------------------------------------------------
!     [u-velocity component: ]
!-----------------------------------------------------
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      term1e=(den(ine)*uu(ine)+(gradu(1,ine)+gradu(1,ine))*(vis(ine)-viscos))*fx(inp)+ &
             (den(inp)*uu(inp)+(gradu(1,inp)+gradu(1,inp))*(vis(inp)-viscos))*(1.-fx(inp))

      term1w=(den(inp)*uu(inp)+(gradu(1,inp)+gradu(1,inp))*(vis(inp)-viscos))*fx(inw)+ &
             (den(inw)*uu(inw)+(gradu(1,inw)+gradu(1,inw))*(vis(inw)-viscos))*(1.-fx(inw))

      term1n=(den(inn)*uu(inn)+(gradu(1,inn)+gradu(1,inn))*(vis(inn)-viscos))*fy(inp)+ &
             (den(inp)*uu(inp)+(gradu(1,inp)+gradu(1,inp))*(vis(inp)-viscos))*(1.-fy(inp))

      term1s=(den(inp)*uu(inp)+(gradu(1,inp)+gradu(1,inp))*(vis(inp)-viscos))*fy(ins)+ &
             (den(ins)*uu(ins)+(gradu(1,ins)+gradu(1,ins))*(vis(ins)-viscos))*(1.-fy(ins))

      term1t=(den(int)*uu(int)+(gradu(1,int)+gradu(1,int))*(vis(int)-viscos))*fz(inp)+ &
             (den(inp)*uu(inp)+(gradu(1,inp)+gradu(1,inp))*(vis(inp)-viscos))*(1.-fz(inp))

      term1b=(den(inp)*uu(inp)+(gradu(1,inp)+gradu(1,inp))*(vis(inp)-viscos))*fz(inb)+ &
             (den(inb)*uu(inb)+(gradu(1,inb)+gradu(1,inw))*(vis(inb)-viscos))*(1.-fz(inb))

      dterm1dx=((term1e-term1w)*ar1x(inp)+ &
                (term1n-term1s)*ar2x(inp)+ &
                (term1t-term1b)*ar3x(inp))
!-----------------------------------------------------------------
!
      term2e=(den(ine)*uv(ine)+(gradu(2,ine)+gradv(1,ine))*(vis(ine)-viscos))*fx(inp)+ &
             (den(inp)*uv(inp)+(gradu(2,inp)+gradv(1,inp))*(vis(inp)-viscos))*(1.-fx(inp))

      term2w=(den(inp)*uv(inp)+(gradu(2,inp)+gradv(1,inp))*(vis(inp)-viscos))*fx(inw)+ &
             (den(inw)*uv(inw)+(gradu(2,inw)+gradv(1,inw))*(vis(inw)-viscos))*(1.-fx(inw))

      term2n=(den(inn)*uv(inn)+(gradu(2,inn)+gradv(1,inn))*(vis(inn)-viscos))*fy(inp)+ &
             (den(inp)*uv(inp)+(gradu(2,inp)+gradv(1,inp))*(vis(inp)-viscos))*(1.-fy(inp))

      term2s=(den(inp)*uv(inp)+(gradu(2,inp)+gradv(1,inp))*(vis(inp)-viscos))*fy(ins)+ &
             (den(ins)*uv(ins)+(gradu(2,ins)+gradv(1,ins))*(vis(ins)-viscos))*(1.-fy(ins))

      term2t=(den(int)*uv(int)+(gradu(2,int)+gradv(1,int))*(vis(int)-viscos))*fz(inp)+ &
             (den(inp)*uv(inp)+(gradu(2,inp)+gradv(1,inp))*(vis(inp)-viscos))*(1.-fz(inp))

      term2b=(den(inp)*uv(inp)+(gradu(2,inp)+gradv(1,inp))*(vis(inp)-viscos))*fz(inb)+ &
             (den(inb)*uv(inb)+(gradu(2,inb)+gradv(1,inw))*(vis(inb)-viscos))*(1.-fz(inb))

      dterm2dy=((term2e-term2w)*ar1y(inp)+ &
                (term2n-term2s)*ar2y(inp)+ &
                (term2t-term2b)*ar3y(inp))
!------------------------------------------------------------------
!
      term3e=(den(ine)*uw(ine)+(gradu(3,ine)+gradw(1,ine))*(vis(ine)-viscos))*fx(inp)+ &
             (den(inp)*uw(inp)+(gradu(3,inp)+gradw(1,inp))*(vis(inp)-viscos))*(1.-fx(inp))

      term3w=(den(inp)*uw(inp)+(gradu(3,inp)+gradw(1,inp))*(vis(inp)-viscos))*fx(inw)+ &
             (den(inw)*uw(inw)+(gradu(3,inw)+gradw(1,inw))*(vis(inw)-viscos))*(1.-fx(inw))

      term3n=(den(inn)*uw(inn)+(gradu(3,inn)+gradw(1,inn))*(vis(inn)-viscos))*fy(inp)+ &
             (den(inp)*uw(inp)+(gradu(3,inp)+gradw(1,inp))*(vis(inp)-viscos))*(1.-fy(inp))

      term3s=(den(inp)*uw(inp)+(gradu(3,inp)+gradw(1,inp))*(vis(inp)-viscos))*fy(ins)+ &
             (den(ins)*uw(ins)+(gradu(3,ins)+gradw(1,ins))*(vis(ins)-viscos))*(1.-fy(ins))

      term3t=(den(int)*uw(int)+(gradu(3,int)+gradw(1,int))*(vis(int)-viscos))*fz(inp)+ &
             (den(inp)*uw(inp)+(gradu(3,inp)+gradw(1,inp))*(vis(inp)-viscos))*(1.-fz(inp))

      term3b=(den(inp)*uw(inp)+(gradu(3,inp)+gradw(1,inp))*(vis(inp)-viscos))*fz(inb)+ &
             (den(inb)*uw(inb)+(gradu(3,inb)+gradw(1,inw))*(vis(inb)-viscos))*(1.-fz(inb))

      dterm3dz=((term3e-term3w)*ar1z(inp)+ &
                (term3n-term3s)*ar2z(inp)+ &
                (term3t-term3b)*ar3z(inp))
!
!------------------------------------------------------------------

      su(inp)=su(inp)-dterm1dx-dterm2dy-dterm3dz+two_thirds*gradte(1,inp)

      end do !!i-loop
      end do !!j-loop
      end do !!k-loop

!
!-----------------------------------------------------
!     [v-velocity component: ]
!-----------------------------------------------------
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      term1e=(den(ine)*uv(ine)+(gradv(1,ine)+gradu(2,ine))*(vis(ine)-viscos))*fx(inp)+ &
             (den(inp)*uv(inp)+(gradv(1,inp)+gradu(2,inp))*(vis(inp)-viscos))*(1.-fx(inp))

      term1w=(den(inp)*uv(inp)+(gradv(1,inp)+gradu(2,inp))*(vis(inp)-viscos))*fx(inw)+ &
             (den(inw)*uv(inw)+(gradv(1,inw)+gradu(2,inw))*(vis(inw)-viscos))*(1.-fx(inw))

      term1n=(den(inn)*uv(inn)+(gradv(1,inn)+gradu(2,inn))*(vis(inn)-viscos))*fy(inp)+ &
             (den(inp)*uv(inp)+(gradv(1,inp)+gradu(2,inp))*(vis(inp)-viscos))*(1.-fy(inp))

      term1s=(den(inp)*uv(inp)+(gradv(1,inp)+gradu(2,inp))*(vis(inp)-viscos))*fy(ins)+ &
             (den(ins)*uv(ins)+(gradv(1,ins)+gradu(2,ins))*(vis(ins)-viscos))*(1.-fy(ins))

      term1t=(den(int)*uv(int)+(gradv(1,int)+gradu(2,int))*(vis(int)-viscos))*fz(inp)+ &
             (den(inp)*uv(inp)+(gradv(1,inp)+gradu(2,inp))*(vis(inp)-viscos))*(1.-fz(inp))

      term1b=(den(inp)*uv(inp)+(gradv(1,inp)+gradu(2,inp))*(vis(inp)-viscos))*fz(inb)+ &
             (den(inb)*uv(inb)+(gradv(1,inb)+gradu(2,inw))*(vis(inb)-viscos))*(1.-fz(inb))

      dterm1dx=(term1e-term1w)*ar1x(inp)+ &
               (term1n-term1s)*ar2x(inp)+ &
               (term1t-term1b)*ar3x(inp)
!-----------------------------------------------------------------
!
      term2e=(den(ine)*vv(ine)+(gradv(2,ine)+gradv(2,ine))*(vis(ine)-viscos))*fx(inp)+ &
             (den(inp)*vv(inp)+(gradv(2,inp)+gradv(2,inp))*(vis(inp)-viscos))*(1.-fx(inp))

      term2w=(den(inp)*vv(inp)+(gradv(2,inp)+gradv(2,inp))*(vis(inp)-viscos))*fx(inw)+ &
             (den(inw)*vv(inw)+(gradv(2,inw)+gradv(2,inw))*(vis(inw)-viscos))*(1.-fx(inw))

      term2n=(den(inn)*vv(inn)+(gradv(2,inn)+gradv(2,inn))*(vis(inn)-viscos))*fy(inp)+ &
             (den(inp)*vv(inp)+(gradv(2,inp)+gradv(2,inp))*(vis(inp)-viscos))*(1.-fy(inp))

      term2s=(den(inp)*vv(inp)+(gradv(2,inp)+gradv(2,inp))*(vis(inp)-viscos))*fy(ins)+ &
             (den(ins)*vv(ins)+(gradv(2,ins)+gradv(2,ins))*(vis(ins)-viscos))*(1.-fy(ins))

      term2t=(den(int)*vv(int)+(gradv(2,int)+gradv(2,int))*(vis(int)-viscos))*fz(inp)+ &
             (den(inp)*vv(inp)+(gradv(2,inp)+gradv(2,inp))*(vis(inp)-viscos))*(1.-fz(inp))

      term2b=(den(inp)*vv(inp)+(gradv(2,inp)+gradv(2,inp))*(vis(inp)-viscos))*fz(inb)+ &
             (den(inb)*vv(inb)+(gradv(2,inb)+gradv(2,inw))*(vis(inb)-viscos))*(1.-fz(inb))

      dterm2dy=(term2e-term2w)*ar1y(inp)+ &
               (term2n-term2s)*ar2y(inp)+ &
               (term2t-term2b)*ar3y(inp)

!------------------------------------------------------------------
!
      term3e=(den(ine)*vw(ine)+(gradv(3,ine)+gradw(2,ine))*(vis(ine)-viscos))*fx(inp)+ &
             (den(inp)*vw(inp)+(gradu(3,inp)+gradw(1,inp))*(vis(inp)-viscos))*(1.-fx(inp))

      term3w=(den(inp)*vw(inp)+(gradv(3,inp)+gradw(2,inp))*(vis(inp)-viscos))*fx(inw)+ &
             (den(inw)*vw(inw)+(gradv(3,inw)+gradw(2,inw))*(vis(inw)-viscos))*(1.-fx(inw))

      term3n=(den(inn)*vw(inn)+(gradv(3,inn)+gradw(2,inn))*(vis(inn)-viscos))*fy(inp)+ &
             (den(inp)*vw(inp)+(gradv(3,inp)+gradw(2,inp))*(vis(inp)-viscos))*(1.-fy(inp))

      term3s=(den(inp)*vw(inp)+(gradv(3,inp)+gradw(2,inp))*(vis(inp)-viscos))*fy(ins)+ &
             (den(ins)*vw(ins)+(gradv(3,ins)+gradw(2,ins))*(vis(ins)-viscos))*(1.-fy(ins))

      term3t=(den(int)*vw(int)+(gradv(3,int)+gradw(2,int))*(vis(int)-viscos))*fz(inp)+ &
             (den(inp)*vw(inp)+(gradv(3,inp)+gradw(2,inp))*(vis(inp)-viscos))*(1.-fz(inp))

      term3b=(den(inp)*vw(inp)+(gradv(3,inp)+gradw(2,inp))*(vis(inp)-viscos))*fz(inb)+ &
             (den(inb)*vw(inb)+(gradv(3,inb)+gradw(2,inw))*(vis(inb)-viscos))*(1.-fz(inb))

      dterm3dz=(term3e-term3w)*ar1z(inp)+ &
               (term3n-term3s)*ar2z(inp)+ &
               (term3t-term3b)*ar3z(inp)
!
!------------------------------------------------------------------

      sv(inp)=sv(inp)-dterm1dx-dterm2dy-dterm3dz+two_thirds*gradte(2,inp)

      end do !!i-loop
      end do !!j-loop
      end do !!k-loop
!
!-----------------------------------------------------
!     [w-velocity component: ]
!-----------------------------------------------------
      do k=3,nkmm
      do i=3,nimm
      do j=3,njmm

      inp=lk(k)+li(i)+j

      term1e=(den(ine)*uw(ine)+(gradw(1,ine)+gradu(3,ine))*(vis(ine)-viscos))*fx(inp)+ &
             (den(inp)*uw(inp)+(gradw(1,inp)+gradu(3,inp))*(vis(inp)-viscos))*(1.-fx(inp))

      term1w=(den(inp)*uw(inp)+(gradw(1,inp)+gradu(3,inp))*(vis(inp)-viscos))*fx(inw)+ &
             (den(inw)*uw(inw)+(gradw(1,inw)+gradu(3,inw))*(vis(inw)-viscos))*(1.-fx(inw))

      term1n=(den(inn)*uw(inn)+(gradw(1,inn)+gradu(3,inn))*(vis(inn)-viscos))*fy(inp)+ &
             (den(inp)*uw(inp)+(gradw(1,inp)+gradu(3,inp))*(vis(inp)-viscos))*(1.-fy(inp))

      term1s=(den(inp)*uw(inp)+(gradw(1,inp)+gradu(3,inp))*(vis(inp)-viscos))*fy(ins)+ &
             (den(ins)*uw(ins)+(gradw(1,ins)+gradu(3,ins))*(vis(ins)-viscos))*(1.-fy(ins))

      term1t=(den(int)*uw(int)+(gradw(1,int)+gradu(3,int))*(vis(int)-viscos))*fz(inp)+ &
             (den(inp)*uw(inp)+(gradw(1,inp)+gradu(3,inp))*(vis(inp)-viscos))*(1.-fz(inp))

      term1b=(den(inp)*uw(inp)+(gradw(1,inp)+gradu(3,inp))*(vis(inp)-viscos))*fz(inb)+ &
             (den(inb)*uw(inb)+(gradw(1,inb)+gradu(3,inw))*(vis(inb)-viscos))*(1.-fz(inb))

      dterm1dx=((term1e-term1w)*ar1x(inp)+ &
                (term1n-term1s)*ar2x(inp)+ &
                (term1t-term1b)*ar3x(inp))
!-----------------------------------------------------------------
!
      term2e=(den(ine)*vw(ine)+(gradw(2,ine)+gradv(3,ine))*(vis(ine)-viscos))*fx(inp)+ &
             (den(inp)*vw(inp)+(gradw(2,inp)+gradv(3,inp))*(vis(inp)-viscos))*(1.-fx(inp))

      term2w=(den(inp)*vw(inp)+(gradw(2,inp)+gradv(3,inp))*(vis(inp)-viscos))*fx(inw)+ &
             (den(inw)*vw(inw)+(gradw(2,inw)+gradv(3,inw))*(vis(inw)-viscos))*(1.-fx(inw))

      term2n=(den(inn)*vw(inn)+(gradw(2,inn)+gradv(3,inn))*(vis(inn)-viscos))*fy(inp)+ &
             (den(inp)*vw(inp)+(gradw(2,inp)+gradv(3,inp))*(vis(inp)-viscos))*(1.-fy(inp))

      term2s=(den(inp)*vw(inp)+(gradw(2,inp)+gradv(3,inp))*(vis(inp)-viscos))*fy(ins)+ &
             (den(ins)*vw(ins)+(gradw(2,ins)+gradv(3,ins))*(vis(ins)-viscos))*(1.-fy(ins))

      term2t=(den(int)*vw(int)+(gradw(2,int)+gradv(3,int))*(vis(int)-viscos))*fz(inp)+ &
             (den(inp)*vw(inp)+(gradw(2,inp)+gradv(3,inp))*(vis(inp)-viscos))*(1.-fz(inp))

      term2b=(den(inp)*vw(inp)+(gradw(2,inp)+gradv(3,inp))*(vis(inp)-viscos))*fz(inb)+ &
             (den(inb)*vw(inb)+(gradw(2,inb)+gradv(3,inw))*(vis(inb)-viscos))*(1.-fz(inb))


      dterm2dy=(term2e-term2w)*ar1y(inp)+ &
               (term2n-term2s)*ar2y(inp)+ &
               (term2t-term2b)*ar3y(inp)
!------------------------------------------------------------------
!
      term3e=(den(ine)*ww(ine)+(gradw(3,ine)+gradw(3,ine))*(vis(ine)-viscos))*fx(inp)+ &
             (den(inp)*ww(inp)+(gradw(3,inp)+gradw(3,inp))*(vis(inp)-viscos))*(1.-fx(inp))

      term3w=(den(inp)*ww(inp)+(gradw(3,inp)+gradw(3,inp))*(vis(inp)-viscos))*fx(inw)+ &
             (den(inw)*ww(inw)+(gradw(3,inw)+gradw(3,inw))*(vis(inw)-viscos))*(1.-fx(inw))

      term3n=(den(inn)*ww(inn)+(gradw(3,inn)+gradw(3,inn))*(vis(inn)-viscos))*fy(inp)+ &
             (den(inp)*ww(inp)+(gradw(3,inp)+gradw(3,inp))*(vis(inp)-viscos))*(1.-fy(inp))

      term3s=(den(inp)*ww(inp)+(gradw(3,inp)+gradw(3,inp))*(vis(inp)-viscos))*fy(ins)+ &
             (den(ins)*ww(ins)+(gradw(3,ins)+gradw(3,ins))*(vis(ins)-viscos))*(1.-fy(ins))

      term3t=(den(int)*ww(int)+(gradw(3,int)+gradw(3,int))*(vis(int)-viscos))*fz(inp)+ &
             (den(inp)*ww(inp)+(gradw(3,inp)+gradw(3,inp))*(vis(inp)-viscos))*(1.-fz(inp))

      term3b=(den(inp)*ww(inp)+(gradw(3,inp)+gradw(3,inp))*(vis(inp)-viscos))*fz(inb)+ &
             (den(inb)*ww(inb)+(gradw(3,inb)+gradw(3,inw))*(vis(inb)-viscos))*(1.-fz(inb))

      dterm3dz=(term3e-term3w)*ar1z(inp)+ &
               (term3n-term3s)*ar2z(inp)+ &
               (term3t-term3b)*ar3z(inp)
!
!------------------------------------------------------------------

      sw(inp)=sw(inp)-dterm1dx-dterm2dy-dterm3dz+two_thirds*gradte(3,inp)

      end do !!j-loop
      end do !!i-loop
      end do !!k-loop


      return
      end
