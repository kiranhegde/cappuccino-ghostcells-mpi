!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
program cappuccino_ghostcells_mpi
!
!***********************************************************************
!
!  Three dimensional finite volume solver for Navier-Stokes equations
!  for structured grids, with collocated variable arrangement.
!  Boundary conditions are treated using ghost cells.
!  
! ./channel <input_file> <inlet_file> <grid_file> <monitor_file> <restart_file> <out_folder_path>      
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  use types
  use parameters
  use indexes
  use geometry
  use coef
  use variables
  use title_mod
  use buoy
  use time_mod
  use gradients
  use fieldManipulation
  use utils

  implicit none

  include 'mpif.h'

  integer :: i, j, k, inp
  integer :: itimes, itimee
  integer :: narg
  real(prec) :: source
  real(prec) :: magUbarStar, rUAw, gragPplus, flowDirection
  real(prec) :: suma,dt
  real :: start, finish
  character(len=2) :: trpn
  character(len=6) :: timechar
!                                                                       
!***********************************************************************
!

  !----------------------------------------------------------------------
  ! MPI start up

    call MPI_INIT(ierr)                   

    call MPI_COMM_SIZE(MPI_COMM_WORLD, & 
                                 nproc,&  
                                 ierr  )  

    call MPI_COMM_RANK(MPI_COMM_WORLD, & 
                                 myid, &  
                                 ierr  )  
    
    this = myid + 1

    if(nproc == 1) then
      nproc = 0
      this = 0
    endif

    write(*,'(2(a,i2))') ' np = ', nproc, ' myid = ', myid
  !----------------------------------------------------------------------



  ! Check command line arguments
  narg=command_argument_count()
  if (narg==0.or.narg<5) then
    write(*,*) 'Not enough arguments - exiting!'
    stop
  endif
  call get_command_argument(1,input_file)
  call get_command_argument(2,inlet_file)
  call get_command_argument(3,grid_file)
  call get_command_argument(4,monitor_file)
  call get_command_argument(5,restart_file)
  call get_command_argument(6,out_folder_path)

  if (myid .eq. 0) then

    ! Open files
    call openfiles

    ! Print cappuccino logo to log file.
    call show_logo

  endif

  ! Read input file
  call read_input

  ! Read grid file
  call read_grid

  ! Initialisation
  call init

  ! Open files for data at monitoring points 
  if(ltransient) then
  open(unit=89,file=trim(out_folder_path)//'/transient_monitoring_points')
  open(unit=91,file=trim(out_folder_path)//'/transient_monitoring_points_names')
  rewind 89
  rewind 91
    do imon=1,mpoints
      read(91, *) trpn
      open(91+imon,file=trim(out_folder_path)//"/transient_monitor_point_"//trpn,access='append')
      ! rewind(91+imon)
    end do
  end if

  ! Initial output
  if (myid .eq. 0) then
    call print_header
  endif

!
!===============================================
!     T i m e   l o o p : 
!===============================================

  if (myid .eq. 0) then
    write(66,'(a)') ' '
    write(66,'(a)') '  Start iteration!'
    write(66,'(a)') ' '
  endif

  itimes=itime+1 
  itimee=itime+numstep

  time_loop: do itime=itimes,itimee 
!
!===============================================
!   Update variables : 
!===============================================

    if(bdf) then
      uoo = uo 
      voo = vo 
      woo = wo 
      teoo = teo 
      edoo = edo
      if (lcal(ien)) too = to 
      if (lcal(ivart)) vartoo = varto 
      if (lcal(icon)) conoo = cono 
    endif
    uo = u 
    vo = v 
    wo = w 
    teo = te 
    edo = ed 
    if (lcal(ien)) to = t         
    if (lcal(ivart)) varto = vart 
    if (lcal(icon)) cono = con 
!
!===============================================
!.....Set inlet boundary conditions at every timestep
!===============================================
    if(itime.eq.itimes) call bcin
!
!===============================================
!.....ITERATION CONTROL MONITOR
!===============================================
!
    ! Courant number report:
    include 'CourantNo.h'

    call abort_mission

! 
!===============================================
!.....ITERATION loop
!===============================================
!
    iteration_loop: do iter=1,maxit

      if (myid .eq. 0 ) then
        call cpu_time(start)
      endif

      ! Update values at ghost cells
      call update_values_at_ghostcells

      ! Calculate velocities. Momentum predictor for PISO.
      if(lcal(iu))    call calcuvw
      ! Update OUTLET BC.
      if(lcal(iu).and..not.const_mflux)    call outbc  
      ! Pressure-velocity coupling. Two options: SIMPLE and PISO
      if(lcal(ip).and.simple)   CALL CALCP
      if(lcal(ip).and.piso)     CALL PISO_multiple_correction
      if(lcal(ip).and.pimple)   CALL PIMPLE_multiple_correction
      ! Turbulence equations
      if(lcal(ite))   call calcscm(te,gradte  ,ite) 
      if(lcal(ied))   call calcscm(ed,graded  ,ied)
      ! Temperature , temperature variance, and concentration eqs.
      if(lcal(ien))   call calcsct(t,   gradt     ,ien)
      if(lcal(ivart)) call calcsct(vart,gradvart  ,ivart)
      if(lcal(icon))  call calcsct(con, gradcon   ,icon)
      ! Update eddy viscosity
      if( lcal(ivis) .and. .not. lles)  call modvis
      ! The SGS viscosity of les model
      if(lcal(ivis) .and. lles)         call calc_vis_les

      call synchronize_processes 

      if (myid .eq. 0) then
        call cpu_time(finish)
        write(66,'(a,g0.3,a)') 'ExecutionTime = ',finish-start,' s'
        write(66,*)
      endif

      !---------------------------------------------------------------------------------------------
      ! Residual normalization, convergence check  
      !---------------------------------------------------------------------------------------------
      do i=1,nphi
      resor(i)=resor(i)*snorin(i)
      end do 

      source=max(resor(iu),resor(iv),resor(iw),resor(ip)) 

      ! Now find global maximum
      call global_max( source )

      if(source.gt.slarge) then
        if ( myid .eq. 0 ) write(66,"(//,10x,a)") "*** Program terminated -  iterations diverge ***" 
        call abort_mission
      endif

      !---------------------------------------------------------------------------------------------
      ! False time steping: jump to next line after the end of time loop
      !---------------------------------------------------------------------------------------------
      if(.not.ltransient) then
          if(source.lt.sormax) exit time_loop
      end if
!      
!--------------------------------------------------------------------------------------------------
!     nestacionar: 
!     imas dve mogucnosti ovde:
! 1) nije konvergirao a potrosio je predvidje ni broj iteracija za ovaj vremenski korak
! 2) konvergirao je pre nego sto je zavrsio sa svim predvidjenim iteracijama 
!--------------------------------------------------------------------------------------------------

      if(ltransient) then 

          ! konverigao u okviru timestep-a ili potrosio sve iteracije
          if(source.lt.sormax.or.iter.eq.maxit) then 

            if(const_mflux) then
             include 'const_massflux_correction.f90'
            endif

            if(mod(itime,nzapis).eq.0) call writefiles
            call writehistory !<- write monitoring points
            call calc_statistics 

            if(mod(itime,50).eq.0) then 

               Open(Unit=87,File=Trim(Out_Folder_Path)//'/tecplot-vel.plt') 
               Rewind 87
               Write(87,*) 'Title     = " Velocity field - snapshot"'
               Write(87,*) 'Variables = "X"'
               Write(87,*) '"Y"'
               Write(87,*) '"Z"'
               Write(87,*) '"U"'
               Write(87,*) '"V"'
               Write(87,*) '"W"'
               Write(87,*) 'Zone T=" "'
               Write(87,*) 'I=',Nimm, ' ,J=',Njmm, ' ,K=',Nkmm,', F=Point'
               do k=2,nkm
               do j=2,njm
               do i=2,nim
               Inp=Lk(K)+Li(I)+J
               Write(87,*) Xc(Inp),Yc(Inp),Zc(Inp),U(Inp),V(Inp),W(Inp)
               Enddo
               Enddo
               Enddo 
               Close(87)

               write(timechar,'(i6)') itime
               call execute_command_line("tec360 -b -p makro-film.mcr")
               call execute_command_line("mv untitled100.png "//timechar//".png")
            endif
           
            cycle time_loop 
          endif

      end if 

    end do iteration_loop

    ! Zapis poslije svakih [nzapis] iteracija:
    if(.not.ltransient) then
      if(mod(itime,nzapis).eq.0.and.itime.ne.numstep) call writefiles
    endif

    if( ltransient .and. myid.eq.0 ) call flush(66)
 
  end do time_loop

  ! False time stepping comes here after time loop with exit command

  if(loute) call report_wall_stresses !<- vrati kad ti bude trebalo y+, Tau na zidovima

  ! Write field values for the next run 
  if(lwrite) call writefiles

 ! call deallocate_arrays

  !----------------------------------------------------------------------
  ! MPI final call
    call MPI_Finalize(ierr)

end program

