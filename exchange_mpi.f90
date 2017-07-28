!
!***********************************************************************
!
  subroutine exchange(phi) 
!
!***********************************************************************
!
!   Exchanges field values between connected processes from their
!   respective buffers. 
!   Executes a blocking send and receive. The same buffer is used both 
!   for the send and for the receive, so that the message sent is 
!   replaced by the message received.                   
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use mpi_exchange

  implicit none

  include 'mpif.h'

  integer :: i,iDomain,iDFriend,iStart,iEnd
  integer :: rectag,sendtag,iarbitr
  integer :: length
  integer :: status(mpi_status_size)


  real(dp), dimension(nxyza) :: phi

  iarbitr = 20 ! Helps create message tag (an integer btw)
 
  ! >> Fill the buffers with new values

  ! ! Idi po svim domenima sa kojima je ovaj konektovan
  ! do iDomain = 1,size(neighbour)

  !   len = ioffset_buf(iDomain+1)-ioffset_buf(iDomain)

  !   do i=1,len
  !     buffer(i) = phi( bufind(i) )
  !   enddo

  ! enddo


  ! Moze i da se pojendostavi punjenje buffera jer je sve vec podeseno
  do i=1,lenbuf
    buffer(i) = phi( bufind(i) )
  enddo 


! >> Exchange the values

  ! Idi po svim domenima sa kojima je ovaj konektovan
  do iDomain = 1,size(neighbour)

    iDFriend = neighbour(iDomain)

    sendtag = 123 + this + iDFriend !iarbitr*this + iDFriend    ! tag for sending
    rectag  = sendtag !iarbitr*iDFriend + this    ! tag for receiving

    iStart = ioffset_buf(iDomain)
    iEnd = ioffset_buf(iDomain+1)-1

    length = ioffset_buf(iDomain+1)-ioffset_buf(iDomain) ! ...also iEnd-iStart+1

    call MPI_SENDRECV_REPLACE & 
     (buffer(iStart),   &   ! buffer  salje donju vrednost jer je ona najjniza negativna odatle ide u plus za jedan po jedan, ali vrednosti mora da su contiguous u phi
      length,           &     ! length   
      MPI_DOUBLE_PRECISION, & ! datatype  
      iDFriend,          &    ! dest,      
      sendtag,          &     ! sendtag,    
      iDFriend,         &    ! source,      
      rectag,           &     ! recvtag,      
      MPI_COMM_WORLD,   &     ! communicator
      status,           &     ! status
      ierr)                   ! error

  end do

  ! Prebaci iz buffera u phi polje 
  do i=1,lenbuf
    phi( bufind(i) ) = buffer(i)
  enddo
 
  end subroutine exchange