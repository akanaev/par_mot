      subroutine compare(ar_p,ar,im,jm,km,var)
       
      include 'mpif.h'             
!        include 'par.inc'
!        include 'com.inc'
      character(len=*) var
      COMMON / MPI311 / IPROC,NUMPROC, KB1,KE1,KB2,KE2
      COMMON / MPI313 / KBA1(0:16384),KEA1(0:16384),
     *                  KBA2(0:16384),KEA2(0:16384),
     *                  KBA3(0:16384),KEA3(0:16384)
      real
     # ar_p(KB1:KE1,KB2:KE2,km),ar(im,jm,km)
     #, buf(im,jm,km)
      integer,allocatable:: rcounts(:),displs(:)     
      real,allocatable:: rbuf2(:),sbuf2(:)
       integer :: alloc_status
      allocate(rcounts(numproc),displs(numproc))
      allocate(sbuf2((ke1-kb1+1)*(ke2-kb2+1)*km))
      allocate(rbuf2(im*jm*km))
!          allocate(rbuf2(im*jm*km), stat=alloc_status)
!      if ( alloc_status /= 0 ) then
!	print *, 'Allocate failed.  status = ', alloc_status
!	stop
!      endif
      
      displs(1)=0
      do i=1,numproc
        i0=kea1(i-1)-kba1(i-1)+1
        j0=kea2(i-1)-kba2(i-1)+1
        k0=km
        itmp=i0*j0*k0
        rcounts(i)=itmp
        if(i<numproc)then
          displs(i+1)=displs(i)+itmp
!           if(iproc.eq.0) print *,'rank=',IPROC,'i',i
!     *     ,"displs(",i+1,")=",displs(i+1)
!     *     ,"displs(",i,")=",displs(i)
!     *     ,"itmp=",itmp  
        endif
      enddo
      m=0
      do k=1,km
        do j=kb2,ke2
          do i=kb1,ke1
          sbuf2(m+1)=ar_p(i,j,k)
          m=m+1
          enddo
        enddo
      enddo
      itmp=(ke1-kb1+1)*(ke2-kb2+1)*km
      
!      if(iproc.eq.0) print *,'rank=',IPROC
!     *,'kba1(0)=',kba1(0),'kba1(1)=',kba1(1),
!     *'kea1(0)=',kea1(0),'kea1(1)=',kea1(1)
!     *,'kba2(0)=',kba2(0),'kba2(1)=',kba2(1),
!     *'kea2(0)=',kea2(0),'kea2(1)=',kea2(1)
     
      CALL MPI_GATHERV(sbuf2(1:itmp),itmp,MPI_REAL,rbuf2
     #,RCOUNTS, DISPLS, MPI_REAL ,0, MPI_COMM_WORLD, ierr)
 

      
      if(iproc.eq.0)then
      
!!      open(file=var, unit=888, access="stream",form="unformatted")
!      open(file=var, unit=888)
!      write(888,2012)rbuf2      
!2012  format(e12.6)      
!!      write(888) rbuf2
!      close(888)
      
      m2=0
      do ii=0,numproc-1
        do k=1,km
            do j=kba2(ii),kea2(ii)
                do i=kba1(ii),kea1(ii)
                  buf(i,j,k)=rbuf2(m2+1)
                   if(buf(i,j,k).ne.ar(i,j,k)) then
                    write(*,*) 'Error in ',var,' compare: i=',i,' j=',j
     *               ,' k=',k
                    write(*,*) 'parallel(i,j,k) = ', buf(i,j,k)
                    write(*,*) 'single(i,j,k) =  ', ar(i,j,k)
                    stop
                   endif                  
                  m2=m2+1
                enddo
            enddo
        enddo
      enddo
      endif
      
!!      open(file=var, unit=888, access="stream",form="unformatted")
!      open(file=var, unit=888)
!      write(888,2012)buf    
!2012  format(e12.6)      
!!      write(888) rbuf2
!      close(888)      

      return
      end
      
      subroutine gather2d(sbuf,rbuf,im,jm,km,var)
      include 'mpif.h'             
!        include 'par.inc'
!        include 'com.inc'
      character(len=*) var
      COMMON / MPI311 / IPROC,NUMPROC, KB1,KE1,KB2,KE2
      COMMON / MPI313 / KBA1(0:16384),KEA1(0:16384),
     *                  KBA2(0:16384),KEA2(0:16384),
     *                  KBA3(0:16384),KEA3(0:16384)
      real
     * sbuf(KB1:KE1,KB2:KE2),rbuf(im,jm)

      integer,allocatable:: rcounts(:),displs(:)     
      real,allocatable:: sbuf2(:),rbuf2(:)
       integer :: alloc_status
      allocate(rcounts(numproc),displs(numproc))
      allocate(sbuf2((ke1-kb1+1)*(ke2-kb2+1)))
      allocate(rbuf2(im*jm))     
      displs(1)=0
      do i=1,numproc
        i0=kea1(i-1)-kba1(i-1)+1
        j0=kea2(i-1)-kba2(i-1)+1
        itmp=i0*j0
        rcounts(i)=itmp
        if(i<numproc)then
          displs(i+1)=displs(i)+itmp
!           if(iproc.eq.0) print *,'rank=',IPROC,'i',i
!     *     ,"displs(",i+1,")=",displs(i+1)
!     *     ,"displs(",i,")=",displs(i)
!     *     ,"itmp=",itmp  
        endif
      enddo
      m=0
        do j=kb2,ke2
          do i=kb1,ke1
          sbuf2(m+1)=sbuf(i,j)
          m=m+1
          enddo
        enddo
      itmp=(ke1-kb1+1)*(ke2-kb2+1)

      CALL MPI_ALLGATHERV(sbuf2(1:itmp),itmp,MPI_REAL,rbuf2
     #,RCOUNTS, DISPLS, MPI_REAL , MPI_COMM_WORLD, ierr)
     
      m2=0
      do ii=0,numproc-1
            do j=kba2(ii),kea2(ii)
                do i=kba1(ii),kea1(ii)
                  rbuf(i,j)=rbuf2(m2+1)         
                  m2=m2+1
                enddo
            enddo
        enddo

      return
      end
