      subroutine cgstab_sip(fi,ifi)
!
!***********************************************************************
!
!    This routine incorporates the CGSTAB solver for seven-diagonal,
!    non-symmetric coefficient matrices (suitable for convection/
!    diffusion problems). See Sect. 5.3.7 for details. Array index
!    IJK calculated from indices I, J, and K according to Table 3.1.
!
!    Writen by Samir Muzaferija, Institut fuer Schiffbau, Hamburg, 1995.
!
!    Original routine modified to incorporate SIP as a preconditioner
!    by Nikola Mirkov, 28.01.2014. nmirkov@vinca.rs
!
!***********************************************************************
!
      use types
      use parameters
      use indexes
      use coef
      use coefb
      use title_mod

      implicit none
!
!***********************************************************************
!
      integer, intent(in) :: ifi
      real(prec), dimension(nxyza) :: fi 

!
!     Local variables
!
      integer :: i, j, k, ijk, ns, l
      real(prec), dimension(nxyza) :: reso,pk,uk,zk,vk
      real(prec) :: res0, rsm, resmax, resl, p1, p2, p3
      real(prec) :: alf, beto, gam, bet, om, vres, vv, ukreso

!.....max no. of inner iters
      resmax = sor(ifi) 
!
!.....calculate initial residual vector
!
      res0=0.0d0
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
             ijk=lk(k)+li(i)+j
             res(ijk)=ae(ijk)*fi(ijk+nj)+aw(ijk)*fi(ijk-nj)+an(ijk)* &
             fi(ijk+1)+as(ijk)*fi(ijk-1)+at(ijk)*fi(ijk+nij)+ &
             ab(ijk)*fi(ijk-nij)+su(ijk)-ap(ijk)*fi(ijk)
             res0=res0+abs(res(ijk))
          end do
        end do
      end do
!
      if(ltest) write(66,*) '                   0  SWEEP, RES0 = ',RES0

!.....calculate coefficients of  l  and  u  matrices using sip
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
          ijk=lk(k)+li(i)+j
          bb(ijk)=-ab(ijk)/(1.+alfa*(bn(ijk-nij)+be(ijk-nij)))
          bw(ijk)=-aw(ijk)/(1.+alfa*(bn(ijk-nj)+bt(ijk-nj)))
          bs(ijk)=-as(ijk)/(1.+alfa*(be(ijk-1)+bt(ijk-1)))
          p1=alfa*(bb(ijk)*bn(ijk-nij)+bw(ijk)*bn(ijk-nj))
          p2=alfa*(bb(ijk)*be(ijk-nij)+bs(ijk)*be(ijk-1))
          p3=alfa*(bw(ijk)*bt(ijk-nj)+bs(ijk)*bt(ijk-1))
          bp(ijk)=1./(ap(ijk)+p1+p2+p3 &
            -bb(ijk)*bt(ijk-nij) &
            -bw(ijk)*be(ijk-nj) &
            -bs(ijk)*bn(ijk-1)+small)
           bn(ijk)=(-an(ijk)-p1)*bp(ijk)
           be(ijk)=(-ae(ijk)-p2)*bp(ijk)
           bt(ijk)=(-at(ijk)-p3)*bp(ijk)
          end do
        end do
      end do
!
!.....initialize working arrays and constants
!
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            reso(ijk)=res(ijk)
            pk(ijk)=0.0d0
            uk(ijk)=0.0d0
            zk(ijk)=0.0d0
            vk(ijk)=0.0d0
          end do
        end do
      end do
      alf=1.0d0
      beto=1.0d0
      gam=1.0d0
!
!.....start inner iterations
!
      ns=nsw(ifi)
      do l=1,ns
!
!..... calculate beta and omega
!
      bet=0.0d0
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            bet=bet+res(ijk)*reso(ijk)
          end do
        end do
      end do
      om=bet*gam/(alf*beto+small)
      beto=bet
!
!..... calculate pk
!
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            pk(ijk)=res(ijk)+om*(pk(ijk)-alf*uk(ijk))
          end do
        end do
      end do
!
!.....solve (m zk = pk) - forward substitution
!
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            zk(ijk)=(pk(ijk)-bb(ijk)*zk(ijk-nij)-bw(ijk)*zk(ijk-nj)- &
            bs(ijk)*zk(ijk-1))*bp(ijk)
          end do
        end do
      end do

!
!..... backward elimination
!
      do k=nkm,2,-1
        do i=nimm,3,-1
          do j=njmm,3,-1
            ijk=lk(k)+li(i)+j
            zk(ijk)=zk(ijk)-bn(ijk)*zk(ijk+1)-be(ijk)*zk(ijk+nj)- &
                    bt(ijk)*zk(ijk+nij)
          end do
        end do
      end do
!
!.....calculate uk = a.pk
!
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            uk(ijk)=ap(ijk)*zk(ijk)-ae(ijk)*zk(ijk+nj) -    &
                     aw(ijk)*zk(ijk-nj)-an(ijk)*zk(ijk+1)-  &
                     as(ijk)*zk(ijk-1)-at(ijk)*zk(ijk+nij)- &
                     ab(ijk)*zk(ijk-nij)
          end do
        end do
      end do
!
!..... calculate scalar product uk.reso and gamma
!
      ukreso=0.0d0
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            ukreso=ukreso+uk(ijk)*reso(ijk)
          end do
        end do
      end do
      gam=bet/ukreso
!
!.....update (fi) and calculate w (overwrite res - it is res-update)
!
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            fi(ijk)=fi(ijk)+gam*zk(ijk)   
            res(ijk)=res(ijk)-gam*uk(ijk) !w
          end do
        end do
      end do
!
!.....solve (m y = w); y overwrites zk; forward substitution
!
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            zk(ijk)=(res(ijk)-bb(ijk)*zk(ijk-nij)-bw(ijk)*zk(ijk-nj)- &
            bs(ijk)*zk(ijk-1))*bp(ijk)
          end do
        end do
      end do
!
!.....backward substitution
!
      do k=nkm,2,-1
        do i=nimm,3,-1
          do j=njmm,3,-1
            ijk=lk(k)+li(i)+j
            zk(ijk)=zk(ijk)-bn(ijk)*zk(ijk+1)-be(ijk)*zk(ijk+nj)- &
                    bt(ijk)*zk(ijk+nij)
          end do
        end do
      end do
!
!.....calculate v = a.y (vk = a.zk)
!
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            vk(ijk)=ap(ijk)*zk(ijk)   -ae(ijk)*zk(ijk+nj)-  &
                    aw(ijk)*zk(ijk-nj)-an(ijk)*zk(ijk+1)-   &
                    as(ijk)*zk(ijk-1) -at(ijk)*zk(ijk+nij)- &
                    ab(ijk)*zk(ijk-nij)
          end do
        end do
      end do
!
!..... calculate alpha (alf)
!
      vres=0.0d0
      vv=0.0d0
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            vres=vres+vk(ijk)*res(ijk)
            vv=vv+vk(ijk)*vk(ijk)
          end do
        end do
      end do

      alf=vres/(vv+small)
!
!.....update variable (fi) and residual (res) vectors
!
      resl=0.0d0
      do k=2,nkm
        do i=3,nimm
          do j=3,njmm 
            ijk=lk(k)+li(i)+j
            fi(ijk)=fi(ijk)+alf*zk(ijk)
            res(ijk)=res(ijk)-alf*vk(ijk)
            resl=resl+abs(res(ijk))
          end do
        end do
      end do
!
!.....check convergence
!
      if(l.eq.1) resor(ifi)=res0
      rsm=resl/(resor(ifi)+small)
      if(ltest) write(66,'(19x,3a,I4,a,1PE10.3,a,1PE10.3)') ' FI=',CHVAR(IFI),' SWEEP = ',L,' RESL = ',RESL,' RSM = ',RSM
      if(rsm.lt.resmax) exit
!
!.....END OF ITERATION LOOP
!
      END DO

!.....Write linear solver report:
      write(66,'(3a,i1,a,i1,2a,1PE10.3,a,1PE10.3,a,I0)') &
      'BiCGStab(SIP):  Solving for ',trim(chvarSolver(IFI)),'(icorr=',icorr,',npcor=',ipcorr,')', &
      ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L

      return
      end

