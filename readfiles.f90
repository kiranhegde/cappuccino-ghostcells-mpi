!***********************************************************************
!
subroutine readfiles
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use variables
  use title_mod
  use buoy
  use time_mod
  use statistics
  use utils

  implicit none  
!
!***********************************************************************
! 
  integer :: restart_unit
  character( len = 5) :: nproc_char 

  ! NOTE: nproc_char <- this (=myid + 1) zapisan levo u vidu stringa.
  call i4_to_s_left ( this, nproc_char )
  print*, 'Opening restart file: ',adjustl(trim(restart_file))//'-'//trim(nproc_char)

  call get_unit ( restart_unit )

  open(unit=restart_unit,file=adjustl(trim(restart_file))//'-'//trim(nproc_char),form='unformatted')
  rewind restart_unit


  read(restart_unit) itime,time
  if(const_mflux) read(restart_unit) gradpcmf
  read(restart_unit) f1
  read(restart_unit) f2
  read(restart_unit) f3
  read(restart_unit) u
  read(restart_unit) v
  read(restart_unit) w
  read(restart_unit) p
  read(restart_unit) te
  read(restart_unit) ed
  read(restart_unit) t
  read(restart_unit) vis
  !read(restart_unit) vart
  !read(restart_unit) edd
  !read(restart_unit) ret
  !read(restart_unit) den
  !read(restart_unit) utt
  !read(restart_unit) vtt
  !read(restart_unit) wtt
  read(restart_unit) uu
  read(restart_unit) vv
  read(restart_unit) ww
  read(restart_unit) uv
  read(restart_unit) uw
  read(restart_unit) vw
  read(restart_unit) uo
  read(restart_unit) vo
  read(restart_unit) wo
  !read(restart_unit) to
  read(restart_unit) teo
  read(restart_unit) edo
  !read(restart_unit) varto
  !read(restart_unit) con
  !read(restart_unit) cono
  !read(restart_unit) alph

  rewind restart_unit
  close (restart_unit)

  if (ltransient) then
  !------------------------------------------------
  !     [read statistics after first collection: ]
  !------------------------------------------------
  open(unit=85,file=trim(out_folder_path)//'/statistics1')   ! <- n_sample is here, statistics restart file 1
  open(unit=86,file=trim(out_folder_path)//'/statistics2')   ! <- u_aver, v_aver,... are here, statistics restart file 2
  rewind 85
  rewind 86

  read(85,*) n_sample
  read(86,*) u_aver,v_aver,w_aver, &
            uu_aver,vv_aver,ww_aver, &
            uv_aver,uw_aver,vw_aver, &
            te_aver
  close (85)
  close (86)
  endif

end subroutine
