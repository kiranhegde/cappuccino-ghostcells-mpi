!***********************************************************************
!
   subroutine writehistory
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use variables
  use title_mod
  use time_mod
  use omega_turb_models

  implicit none
!
!***********************************************************************
!
 integer :: i, j, k, ijk

!---------------------------------------------------------
!  Results at each time-step for transient simulation
!---------------------------------------------------------

  if(ltransient) then
    do imon=1,mpoints
      read(89,*) i,j,k
      ijk=lk(k)+li(i)+j
      write(91+imon,'(2x,1p7e14.5,2x)') time,u(ijk),v(ijk),w(ijk),te(ijk),ed(ijk)
    end do
    rewind 89
  end if

  return
  end
