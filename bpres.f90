!***********************************************************************
! 
  subroutine bpres(p)
!
!***********************************************************************
! 
!     Set pressures at ghost cells by extrapolation from inner cells.
!
!***********************************************************************
! 
  use types
  use parameters
  use indexes
  use geometry

  implicit none
!
!***********************************************************************
! 
  real(prec), dimension(nxyza) :: p

! Locals:
  integer :: i, j, k, ijk


!---------------------------------------------------------
! bottom boundary
!---------------------------------------------------------
  do i=3,nimm
  do j=3,njmm
    ijk=lk(3)+li(i)+j
    p(ijk-nij)=p(ijk)-(p(ijk+nij)-p(ijk))*fz(ijk)
  enddo
  enddo
!---------------------------------------------------------
! top boundary
!---------------------------------------------------------
  do i=3,nimm
  do j=3,njmm
    ijk=lk(nkmm)+li(i)+j
    p(ijk+nij)=p(ijk)+(p(ijk)-p(ijk-nij))*(1.0_dp-fz(ijk-nij))
  enddo
  enddo
!---------------------------------------------------------
! west boundary
!---------------------------------------------------------
  do j=3,njmm
  do k=3,nkmm
    ijk=lk(k)+li(3)+j
    p(ijk-nj)=p(ijk)-(p(ijk+nj)-p(ijk))*fx(ijk)
  enddo
  enddo
!---------------------------------------------------------
! east boundary
!---------------------------------------------------------
  do j=3,njmm
  do k=3,nkmm
    ijk=lk(k)+li(nimm)+j
    p(ijk+nj)=p(ijk)+(p(ijk)-p(ijk-nj))*(1.0_dp-fx(ijk-nj))
  enddo
  enddo
!---------------------------------------------------------
! south boundary
!---------------------------------------------------------
  do i=3,nimm
  do k=3,nkmm
    ijk=li(i)+lk(k)+3
    p(ijk-1)=p(ijk)-(p(ijk+1)-p(ijk))*fy(ijk)
  enddo
  enddo
!---------------------------------------------------------
! north boundary
!---------------------------------------------------------
  do i=3,nimm
  do k=3,nkmm
    ijk=lk(k)+li(i)+njmm
    p(ijk+1)=p(ijk)+(p(ijk)-p(ijk-1))*(1.0_dp-fy(ijk-1))
  enddo
  enddo

  end subroutine
