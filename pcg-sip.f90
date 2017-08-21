!***********************************************************************
!
  subroutine pcgsip(fi,ifi)
!
!***********************************************************************
!
!    This routine incorporates the sip preconditioned 
!    conjugate gradient solver for symmetric matrices in 3d problems
!    with seven-diagonal matrix structure.
!
!    Written by Nikola Mirkov, 28.01.2014. nmirkov@vinca.rs
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
! Local variables
!
  integer :: i, j, k, ijk, ns, l
  real(prec), dimension(nxyza) :: pk,zk
  real(prec) :: rsm, resmax, res0, resl, p1, p2, p3
  real(prec) :: s0, sk, alf, bet, pkapk

! max no. of inner iters
  resmax = sor(ifi)
!
! initalize working arrays
!
  pk = 0.0_dp
  zk = 0.0_dp
  res = 0.0_dp
!
! Calculate initial residual vector and the norm
!
  res0=0.0d0
  do k=3,nkmm
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

  call global_sum(res0)
!
! if ltest=true, print the norm 
!
  if(ltest) write(66,'(a,1pe10.3)') '                    res0 = ',res0

! alculate coefficients of  l  and  u  matrices using sip
  do k=3,nkmm
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
  s0=1.e20
!
! Start inner iterations
!
  ns=nsw(ifi)
  do l=1,ns
!
! solve for zk(ijk) -- forward elimination
!
  do k=3,nkmm
    do i=3,nimm
      do j=3,njmm 
        ijk=lk(k)+li(i)+j
        zk(ijk)=(res(ijk)-bb(ijk)*zk(ijk-nij)-bw(ijk)*zk(ijk-nj)- &
        bs(ijk)*zk(ijk-1))*bp(ijk)
      end do
    end do
  end do

!
!  backward substitution; calculate scalar product sk
!
  sk=0.0d0
  do k=nkmm,3,-1
    do i=nimm,3,-1
      do j=njmm,3,-1
        ijk=lk(k)+li(i)+j
        zk(ijk)=zk(ijk)-bn(ijk)*zk(ijk+1)-be(ijk)*zk(ijk+nj)- &
                bt(ijk)*zk(ijk+nij)
        sk=sk+res(ijk)*zk(ijk)
      end do
    end do
  end do

  call global_sum(sk)
!
! calculate beta
!
  bet=sk/s0
!
! calculate new search vector pk
!
  do k=3,nkmm
    do i=3,nimm
      do j=3,njmm 
        ijk=lk(k)+li(i)+j
        pk(ijk)=zk(ijk)+bet*pk(ijk)
      end do
    end do
  end do
!
!  calculate scalar product (pk.a pk) and alpha (overwrite zk)
!
  pkapk=0.0d0
  do k=3,nkmm
    do i=3,nimm
      do j=3,njmm 
        ijk=lk(k)+li(i)+j
        zk(ijk)=ap(ijk)*pk(ijk)-ae(ijk)*pk(ijk+nj)-                &
          aw(ijk)*pk(ijk-nj)-an(ijk)*pk(ijk+1)-as(ijk)*pk(ijk-1)-  &
          at(ijk)*pk(ijk+nij)-ab(ijk)*pk(ijk-nij)
        pkapk=pkapk+pk(ijk)*zk(ijk)
      end do
    end do
  end do

  call global_sum(pkapk)

  alf=sk/pkapk
!
! calculate variable correction, new residual vector, and norm
!
  resl=0.0d0
  do k=3,nkmm
    do i=3,nimm
      do j=3,njmm 
        ijk=lk(k)+li(i)+j

        fi(ijk)=fi(ijk)+alf*pk(ijk)

        res(ijk)=res(ijk)-alf*zk(ijk)

        resl=resl+abs(res(ijk))

      end do
    end do
  end do

  call global_sum(resl)
  call exchange(fi)
  call exchange(res)

  s0=sk
!
! Check convergence
!
  if(l.eq.1) resor(ifi)=res0
  rsm=resl/(resor(ifi)+small)
  if(ltest) write(66,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',chvar(ifi),' sweep = ',l,' resl = ',resl,' rsm = ',rsm
  if(rsm.lt.resmax) exit
!
! End of iteration loop
!
  end do

! Write linear solver report:
  write(66,'(3a,1pe10.3,a,1pe10.3,a,i0)') &
  'PCG(SIP):  Solving for ',trim(chvarsolver(ifi)), &
  ', Initial residual = ',res0,', fFinal residual = ',resl,', No Iterations ',l
!
  return
  end subroutine
