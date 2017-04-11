c      program parallel motovskiy bay with asm

      include 'mpif.h'
      include 'par.inc'
      include 'com.inc'
      include 'asm.inc'

      common ns,i1,i2,j1,j2,hours
      common /const/ c2,c4,c6,c7,bet1,bet2,rfcr

      COMMON / MPI311 / IPROC,NUMPROC, KB1,KE1,KB2,KE2
      COMMON / MPI312 / NPX,NPY,IPX,IPY
      COMMON / MPI313 / KBA1(0:16384),KEA1(0:16384),
     *                  KBA2(0:16384),KEA2(0:16384)
      COMMON / MPI314 / IPXA(0:16384),IPYA(0:16384)
      COMMON / MPI315 / ileft,iright,itop,ibottom,inear,ifar
      COMMON / MPI316 / itopleft,itopright
     *                 ,ibottomleft,ibottomright 

      real buf
      integer num_lq
      character(len=1024) :: fname
      integer(kind=MPI_OFFSET_KIND)::OFFSET      
            
!    real*8 progtime
      real*8 timestep(1200000)
      real a,calendar(12),dlta(19),startD(4),curentD(4)
     *,t1(10,19),s1(10,19),p1(10),energ(10000)
     &,emid(700),emax(700),emin(700),tmid(700)
     &,tam(6),txm(6),tym(6),sg(6),tg(6),vd(6),vo(6),tr(6)
     &,dzg(12),taa(24),saa(24),avp(24)
     &,levelsolo(2000),levsolo(504),pare(5)
 
      real,allocatable::  qss(:,:),sa(:,:),modu(:,:)
     *,thx(:,:),thy(:,:),work(:,:)
     *,tax(:,:),tay(:,:),work1(:,:)
     *,u(:,:,:),v(:,:,:),w(:,:,:)
     *,u_p(:,:,:),v_p(:,:,:),w_p(:,:,:)  !velocity components     
     *,ut(:,:,:),vt(:,:,:)
     *,ut_p(:,:,:),vt_p(:,:,:)     
     *,s(:,:,:),t(:,:,:),at(:,:,:),st(:,:,:)
     *,t_p(:,:,:),s_p(:,:,:),at_p(:,:,:),st_p(:,:,:)
     *,ph(:,:),php(:,:),pht(:,:),phtt(:,:)
     *,ta(:,:),qt(:,:),qb(:,:),eph(:,:)
     &,uss(:,:,:),wt(:,:,:),pz(:,:,:)
     &,ssp(:,:,:),ssp_p(:,:,:)     
     *,ae(:,:),aw(:,:),as(:,:),an(:,:),ap(:,:)
     *,af(:,:),aff(:,:),x1(:),uh(:,:),vh(:,:)
     *,wr1(:,:),wr(:,:)
     *,fu(:,:,:), fv(:,:,:),us(:,:,:)
     *,hp(:,:),rot(:,:,:),romt(:,:,:)
     &,zxt(:,:),zxs(:,:),zxro(:,:),zxahz(:,:)
     &,dzita(:,:),ttemp(:,:),stemp(:,:)
     &,rotemp(:,:),bottom(:,:),sect(:,:)
     &,dzit1(:,:),aft(:,:)                    
     &,ahz(:,:,:),pr(:,:,:),frf(:,:,:)
     &,tb(:,:,:),te(:,:,:)
     &,tbt(:,:,:),tet(:,:,:)
     &,a3(:),b3(:),c3(:),fuu(:),fvv(:)
     &,e2(:),d2(:),r2(:),g2(:),ef(:),azz(:)
     &,ua(:),va(:)
     &,ala(:),fi(:),pph(:,:)
     &,pmax(:,:),pc(:,:),rmax(:,:),rc(:,:)
     &,ro(:,:,:)
     &,aa(:)
     &,hhbot(:,:),hbot(:,:),bot(:,:),ptbt(:,:)
     &,ptb(:,:)
     &,pi(:,:), ppa(:,:,:),tta(:,:,:)
     &,tam1(:,:),tam2(:,:),scl(:,:,:,:)
     &,pam1(:,:),pam2(:,:),sintrp(:,:,:)
     &,soct(:,:,:),snov(:,:,:)
     &,saug(:,:,:),ssep(:,:,:),res(:,:),smay(:,:,:)
     &,taug(:,:,:),tsep(:,:,:),tmay(:,:,:)
     &,toct(:,:,:),tnov(:,:,:)
     &,tdec(:,:,:),sdec(:,:,:)
     &,tjan(:,:,:),sjan(:,:,:)
     &,Tst(:,:),Sst(:,:),rmseT(:),rmseS(:)
     &,ppT(:),ppS(:),Rlong(:,:),Rlatt(:,:)
     &,h(:,:)
     &,Uf(:,:,:),Vf(:,:,:),Wf(:,:,:)
!      &,att(:,:,:),tf(:,:,:)
!      &,tft(:,:,:),stt(:,:,:)
!      &,sf(:,:,:),sft(:,:,:)

       integer,allocatable:: kh(:,:),iiu(:,:)

       integer inx(12),iny(12)


      real Tsta(70,500),Ssta(70,500),dmin(80),xyz(3,500),Timsta(70)
      real zTS(70,500)
      real Slong(70),Slatt(70)

      integer iAi(500),iAj(500),iAns(1000),kTS(70)


!!!!!!!!!!!!!!!!!!!!!!!!!!! for  POLY_WS_2007  !!!!!!!!!!!!!!!!!
      Character*8   file11(8)
      Character*8   file12(8)
      Character*8   file13(8)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Blagodat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      Character*8   file_bl
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Blagodat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      integer bufsize
      real sbuf,rbuf
      integer, dimension (MPI_STATUS_SIZE) :: istatus
      integer im,jm,km
       
      data file11/'U1','U2','U3','U4','U5','U6','U7','U8'/
      data file12/'V1','V2','V3','V4','V5','V6','V7','V8'/
      data file13/'W1','W2','W3','W4','W5','W6','W7','W8'/
!!!!!!!!!!!!!!!!!!!!!!!!!!! for POLY_WS_2007  !!!!!!!!!!!!!!!!!!
      parameter (fpi=3.141592)
      character file1*11
      character file2*11
      character file3*11
      character*2 fnum(25)
      character*1 name


!c--------------------------------------------------

      call MPI_INIT( ierr )
      if (ierr.ne.MPI_SUCCESS) call mpi_abort (MPI_COMM_WORLD, -1001)
      call MPI_COMM_RANK( MPI_COMM_WORLD, IPROC,  ierr )
      if (ierr.ne.MPI_SUCCESS) call mpi_abort (MPI_COMM_WORLD, -1002)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROC,ierr ) 
      if (ierr.ne.MPI_SUCCESS) call mpi_abort (MPI_COMM_WORLD, -1003)
      
      n1=in
      n2=jn
      
ctext decomposition
      if(IPROC.eq.0) call InputDecompNum (NPX,NPY,NUMPROC)
      call MPI_BCAST(NPX,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (ierr.ne.MPI_SUCCESS) call mpi_abort (MPI_COMM_WORLD, -1007)
      call MPI_BCAST(NPY,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (ierr.ne.MPI_SUCCESS) call mpi_abort (MPI_COMM_WORLD, -1008)
      call MPI_BCAST(NPZ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (ierr.ne.MPI_SUCCESS) call mpi_abort (MPI_COMM_WORLD, -1009)
      
      call Decomp(n1,n2,n3)      

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        print *,'RANK=',IPROC,'IPX=',IPX,'IPY=',IPY,'KB1=',KB1,
     *'KE1=',KE1,'KB2=',KB2,'KE2=',KE2     
        print *,'----------------------'
c--------------------------------------------------  

c     neighbours
      ileft = IPROC-1
      if (IPX.eq.0) ileft = MPI_PROC_NULL
      iright = IPROC+1
      if (IPX+1.eq.NPX) iright = MPI_PROC_NULL

      ibottom = IPROC-NPX
      if (IPY.eq.0) ibottom = MPI_PROC_NULL
      itop = IPROC+NPX
      if (IPY+1.eq.NPY) itop = MPI_PROC_NULL

      itopleft = MPI_PROC_NULL
      if(itop.ne.MPI_PROC_NULL.and.ileft.ne.MPI_PROC_NULL)
     * itopleft = itop - 1

      itopright = MPI_PROC_NULL
      if(itop.ne.MPI_PROC_NULL.and.iright.ne.MPI_PROC_NULL)
     * itopright = itop + 1
     
      ibottomleft = MPI_PROC_NULL
      if(ibottom.ne.MPI_PROC_NULL.and.ileft.ne.MPI_PROC_NULL)
     * ibottomleft = ibottom - 1      

      ibottomright = MPI_PROC_NULL
      if(ibottom.ne.MPI_PROC_NULL.and.iright.ne.MPI_PROC_NULL)
     * ibottomright = ibottom + 1  
           
!!c     neighbours
!        call MPI_Barrier(MPI_COMM_WORLD,ierr)
!      print *,'RANK=',IPROC,' ileft=',ileft,' iright=',iright
!     *,' ibottom=',ibottom,' itop=',itop,' itopleft=',itopleft 
!     *,' itopright=',itopright,' ibottomleft=',ibottomleft
!     *,' ibottomright=',ibottomright

   
c----------------------------------------------------------

      allocate(  qss(in,jn),sa(in,jn),modu(in,jn)
     *,thx(in,jn),thy(in,jn),work(in,jn)
     *,tax(in,jn),tay(in,jn),work1(in,jn)
     *,u(in,jn,kn),v(in,jn,kn),w(in,jn,kn)
     *     ,u_p(KB1-2:KE1+2,KB2-2:KE2+2,kn)
     *     ,v_p(KB1-2:KE1+2,KB2-2:KE2+2,kn)
     *     ,w_p(KB1-1:KE1+1,KB2-1:KE2+1,kn) 
     *,ut(in,jn,kn),vt(in,jn,kn)
     *     ,ut_p(KB1-2:KE1+2,KB2-2:KE2+2,km)
     *     ,vt_p(KB1-2:KE1+2,KB2-2:KE2+2,km)      
     *,s(in,jn,kn),t(in,jn,kn),at(in,jn,kn),st(in,jn,kn)
     *,t_p(KB1-1:KE1+1,KB2-1:KE2+1,kn)
     *,s_p(KB1-1:KE1+1,KB2-1:KE2+1,kn)
     *,at_p(KB1-1:KE1+1,KB2-1:KE2+1,kn)
     *,st_p(KB1-1:KE1+1,KB2-1:KE2+1,kn)     
     *,ph(in,jn),php(in,jn),pht(in,jn),phtt(in,jn)
     *,ta(in,jn),qt(in,jn),qb(in,jn),eph(in,jn)
     &,uss(in,jn,kn),wt(in,jn,kn),pz(in,jn,kn)
     &,ssp(in,jn,kn),ssp_p(KB1:KE1,KB2:KE2,kn)    
     *,ae(in,jn),aw(in,jn),as(in,jn),an(in,jn),ap(in,jn)
     *,af(in,jn),aff(in,jn),x1(in),uh(in,jn),vh(in,jn)
     *,wr1(in,jn),wr(in,jn)
     *,fu(in,jn,kn), fv(in,jn,kn),us(in,jn,kn)
     *,hp(in,jn),rot(in,jn,kn),romt(in,jn,kn)
     &,zxt(in,kn),zxs(in,kn),zxro(in,kn),zxahz(in,kn)
     &,dzita(in,jn),ttemp(in,jn),stemp(in,jn)
     &,rotemp(in,jn),bottom(in,jn),sect(in,kn)
     &,dzit1(in,jn),aft(in,jn)                    
     &,ahz(in,jn,kn),pr(in,jn,kn),frf(in,jn,kn)
     &,tb(in,jn,kn),te(in,jn,kn)
     &,tbt(in,jn,kn),tet(in,jn,kn)
     &,a3(kn),b3(kn),c3(kn),fuu(kn),fvv(kn)
     &,e2(kn),d2(kn),r2(kn),g2(kn),ef(kn),azz(kn)
     &,ua(kn),va(kn)
     &,ala(in),fi(jn),pph(in,jn)

     &,pmax(in,jn),pc(in,jn),rmax(in,jn),rc(in,jn)
     &,ro(in,jn,kn)
     &,aa(in)
     &,hhbot(in,jn),hbot(in,109),bot(in,55),ptbt(in,jn)
     &,ptb(in,jn)
     &,pi(in,jn), ppa(in,jn,800),tta(in,jn,800)
     &,tam1(in,jn),tam2(in,jn),scl(12,in,jn,kn)
     &,pam1(in,jn),pam2(in,jn),sintrp(in,jn,kn)
     &,soct(in,jn,kn),snov(in,jn,kn)
     &,saug(in,jn,kn),ssep(in,jn,kn),res(in,jn),smay(in,jn,kn)
     &,taug(in,jn,kn),tsep(in,jn,kn),tmay(in,jn,kn)
     &,toct(in,jn,kn),tnov(in,jn,kn)
     &,tdec(in,jn,kn),sdec(in,jn,kn)
     &,tjan(in,jn,kn),sjan(in,jn,kn)
!      &,att(in,jn,kn),tf(in,jn,kn)
!      &,tft(in,jn,kn),stt(in,jn,kn)
!      &,sf(in,jn,kn),sft(in,jn,kn)
     &,kh(inn,jnn),iiu(in,jn)
c      MOTO aug
     &,h(in,jn)
c      MOTO aug
     &,Tst(70,kn),Sst(70,kn),rmseT(kn),rmseS(kn)
     &,ppT(kn),ppS(kn),Rlong(in,jn),Rlatt(in,jn)
     &,Uf(8,7,kn),Vf(8,7,kn),Wf(8,7,kn)
     &,)




     
      ntype=1
      ln=1           
      im=in
      jm=jn
      km=kn
      kmm=km

      data calendar/31,28,31,30,31,30,31,31,30,31,30,31/

      i1=1
      i2=in
      j1=1
      j2=jn

        zw(1)=0.
        do k=2,41
        zw(k)=zw(k-1)+7.
      enddo

!!      goto 65
!      do k=1,41
!       print*,zw(k),k
!      enddo
!!      pause 567
!65    continue

      zw=zw*100.

       
      z(1)=zw(1)       !!!! was deleted
      do 1 k=2,kn
      z(k)=(zw(k)+zw(k-1))/2.
1     continue


!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW WHAT do you need WWWWWWWWWWWWW

!         ntype=1   BARKAR

          ntype=2   ! poly 2012 (test for IC dim4bar)

!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW WHAT WWWWWWWWWWWWWWWWWWWWWWWWW





        ln=1 
c--------------------------------------------------
       ptbt=0.
       n2=10



! WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW HAVE TO READ  WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!   MAIN assimilations comments:

!   iASM=1, that meens assimilations only if dt<0, every time

!   itPAR = 1 means direction LOOPs forward (dt>0 and dt<0), and =0 => otherwise

!   nk - & station in '9 circle' for iRMSE=1

!   ICpar=1  ASM whithout depending of time => 3DVAR for IC (iRMSE=0)
    
!   calculations RSME only for direct task whithout assimilations (iASM=0)

!   all modifications in one place, begining of main (9) circle

!   iDT=1 means direct task  


! STEPs:
! 0.      calculation  init.  boundary  cond.                         iRMSE =0 (iASM=1, dt<0,dt>0)  maybe, do not need, if climate IC !!

! 1.      calculation   i,j,n points are nearest  to stations            iRMSE =1 (iASM=0, dt<0)

! 2.      calculation  with ASM,  any varients                           iRMSE =2 (iASM=1)

! 3.      direct task with    RMSE                                       iRMSE =2 (iASM=0,  iDT=1 - direct task) 
! WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW HAVE TO READ  WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW 



         open(19,file='it')
         read(19,110)it
         close (19)

         open(18,file='kp')
         read(18,117)kp
         close (18)
         
110     format(310i1)
117     format(310i3)         

        open(49,file='Ai')
        open(50,file='Aj')
        open(51,file='Ans')

      icl=0
      rmseT=0.
      rmseS=0.
      ppT=0.
      ppS=0.
     
      ns=1
           
           
      itpar=1  !   parameter of direction the global iterations
      iter=1   !  modify by dump  from COND

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      if(iproc.eq.0) then
      call cond(nstep,nprin,nk,al,al1,ctah,ct,cadv,eps
     & ,qwa,cu1,x1,energ,timestep,ns,taa,saa,itpar,iter,idd
     & ,kh,s,t,ta,u,v,w,tx,ty,thx,thy,ph,pht,ssp
     & ,hp,ahz,aday,p0d,p0md,cst,str,h,relax,af
     & ,mon2,fday,hours,fminutes)
      endif


      call cond_p(nstep,nprin,nk,al,al1,ctah,ct,cadv,eps
     & ,qwa,cu1,x1,energ,timestep,ns,taa,saa,itpar,iter,idd
     & ,kh,s_p,t_p,ta,u_p,v_p,w_p,tx,ty,thx,thy,ph,pht,ssp_p
     & ,hp,ahz,aday,p0d,p0md,cst,str,h,relax,af
     & ,mon2,fday,hours,fminutes)



      call compare(u_p(KB1:KE1,KB2:KE2,1:km),u,im,jm,km,'u')
      call compare(v_p(KB1:KE1,KB2:KE2,1:km),v,im,jm,km,'v')
      call compare(w_p(KB1:KE1,KB2:KE2,1:km),w,im,jm,km,'w')
      call compare(ssp_p(KB1:KE1,KB2:KE2,1:km),ssp,im,jm,km,'ssp')
      call compare(t_p(KB1:KE1,KB2:KE2,1:km),t,im,jm,km,'t')
      call compare(s_p(KB1:KE1,KB2:KE2,1:km),s,im,jm,km,'s')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
!!!!MAKE IC TCr_p FILE FOR PARALLEL READ

!!!      if(iproc.eq.0) then
!!!       open(441,file='ph')
!!!       do j=1,jm
!!!       do i=1,im
!!!        write(441,*)ph(i,j)
!!!       enddo
!!!       enddo
!!!       close(441)    
!!!      endif
!!!
!!!      do k=1,kn
!!!      do j=KB2,KE2
!!!      do i=KB1,KE1
!!!       u_p(i,j,k)=u(i,j,k)
!!!       v_p(i,j,k)=v(i,j,k)
!!!       w_p(i,j,k)=w(i,j,k)
!!!       ssp_p(i,j,k)=ssp(i,j,k)
!!!       t_p(i,j,k)=t(i,j,k)
!!!       s_p(i,j,k)=s(i,j,k)                   
!!!      enddo
!!!      enddo
!!!      enddo

!!!      OFFSET=0
!!!      call pwrite(u_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
!!!      
!!!      OFFSET=OFFSET+in*jn*kn*4
!!!      call pwrite(v_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
!!!      
!!!      OFFSET=OFFSET+in*jn*kn*4
!!!      call pwrite(w_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
!!!      
!!!      OFFSET=OFFSET+in*jn*kn*4
!!!      call pwrite(ssp_p,in,jn,kn,'IC TCr_p',OFFSET)
!!!      
!!!      OFFSET=OFFSET+in*jn*kn*4
!!!      call pwrite(t_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
!!!      
!!!      OFFSET=OFFSET+in*jn*kn*4
!!!      call pwrite(s_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        nswr=1

      if(iproc.eq.0) then
      ut=u
      vt=v
      pht=ph
      at=t
      st=s
      endif
      
      ut_p=u_p
      vt_p=v_p
      pht=ph
      at_p=t_p
      st_p=s_p
      
      !       print*,' ns = ',ns

     
!        ns=102240  !  1 aug. CONTROL => MODIFY

!        ns=ns*180./dt  ! if dt ne 3 min

        nsstart=ns                !  2012
              
!        print*,ns

!         + LOOPS 
         istep=nstep/nprin
         nnstep=3*nprin*(istep-2)+2*nprin
!         - LOOPS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! + ASM  !!!!!!!!!!!!!

!WWWWWWWWWWWWWWWWWWWWWW INIT=0  WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW






      west=31.97  ! MOTO here
      east=34.75

      south=69.2                    
      pnorth=70.0


!    print*,'west,east,south,pnorth',west,east,south,pnorth
!        pause 1    
     
       dyi=(pnorth-south)/real(jn-1)
       dxi=(east-west)/real(in-1)
        
!        print*,'dx(km)=',dxi*grad*rad*sin((90.-(69.+26./60.))*grad)
!     & /100000.,'dy=(km)',dyi*grad*rad/100000.

      do j=1,jn
      do i=1,in
         Rlong(i,j)=west +real(i-1)*dxi
         Rlatt(i,j)=pnorth-real(j-1)*dyi    !  RRR  18.08.2013
      enddo
      enddo

!       goto 416
!      print*,west
!      do i=1,in
!      print*,Rlong(i,10),i   !   RRR 20.08.2013
!      enddo
!      print*,east
!      !pause 1

!      print*,pnorth
!      do j=1,jn
!      print*,Rlatt(10,j),j   !   RRR 20.08.2013
!      enddo
!      print*,south
!      !pause 2
!416   continue

!        if(init.eq.0)then  ! need 0 for start   !   if  1

!c                      2011     MOTO
!       open(1,file='prepstat MOTO R.txt')     !  RRR  18.08.2013



!        Tst=100.  ! IMP  need to determine bottom values, see later
!        Sst=100.

!        zmax=0.

!         nst=70
!         

!         do k=1,28  !  1     &  RRR  18.08.2013

!      print*,'  '
!      print*,' stat & ',k

!      kk=0
!      kkk=1
!      read(1,*)Alt,fLat,num,kz,mon,day,year,hours,fmin,sek,dirW,velW
!!      print*,Alt,fLat,num,kz,mon,day,year,hours,fmin,sek,dirW,velW

!       if(kz.gt.zmax)stHmax=k
!       if(kz.gt.zmax)zmax=kz


!           TimeTotal=(day)*24*3600.+hours*3600.
!     &    +fmin*60. +sek
!           Timsta(k)=TimeTotal-10*24*3600-17*3600-11*60 ! start poly -> time=0
!     

!           Slong(k)=Alt
!           Slatt(k)=fLat

!2      kk=kk+1    
!       read(1,*)iz,tt,ss
!!       print*,kkk,tt,ss
!!       pause 2

!        if(kk.eq.1)Tst(k,1)=tt
!        if(kk.eq.1)Sst(k,1)=ss
!        if(kk.eq.1)print*,' 0 m  ',Tst(k,1),Sst(k,1)
!        if(kk.eq.1)goto 2

!!        print*,mod(kk,iz),kkk,'  = mod, kkk'
!!        pause 3

!        if(mod(iz,7).eq.0)kkk=kkk+1
!        if(mod(iz,7).eq.0)Tst(k,kkk)=tt         
!        if(mod(iz,7).eq.0)Sst(k,kkk)=ss 
!      
!       if(mod(iz,7).eq.0)print*,iz,'<-zw  ',Timsta(k)/3600.,
!     & Tst(k,kkk),Sst(k,kkk)

!       if(kk.ne.kz)goto 2

!!         pause 1

!         enddo  ! 1



!         do k=29,nst  !  2

!      print*,'  '
!      print*,' stat & ',k

!      kkk=1
!      kk=0
!      read(1,*)Alt,fLat,num,kz,mon,day,year,hours,fmin,dirW,velW
!!      print*,Alt,fLat,num,kz,mon,day,year,hours,fmin,dirW,velW


!       if(kz.gt.zmax)stHmax=k
!       if(kz.gt.zmax)zmax=kz

!           TimeTotal=(day)*24*3600.+hours*3600.+fmin*60.+sek
!           Timsta(k)=TimeTotal-10*24*3600-17*3600-11*60 ! start poly -> time=0

!!           if(Timsta(k).lt.100.)then
!           print*,k,Timsta(k)
!!           pause 1
!!                                endif
! 



!           Slong(k)=Alt
!           Slatt(k)=fLat
! 


!3      kk=kk+1    
!       read(1,*)iz,tt,ss
!!       print*,kkk,tt,ss

!        if(kk.eq.1)Tst(k,1)=tt
!        if(kk.eq.1)Sst(k,1)=ss
!        if(kk.eq.1)print*,' 0 m  ',Tst(k,1),Sst(k,1)
!        if(kk.eq.1)goto 3



!        if(mod(iz,7).eq.0)kkk=kkk+1
!        if(mod(iz,7).eq.0)Tst(k,kkk)=tt         
!        if(mod(iz,7).eq.0)Sst(k,kkk)=ss

!       if(mod(iz,7).eq.0)print*,iz,'<-zw  ',Timsta(k)/3600.,
!     & Tst(k,kkk),Sst(k,kkk)


!       if(kk.ne.kz)goto 3

!!        pause 2

!         enddo  ! 2      &  RRR  18.08.2013

!       



!!     CONTROL   IMP
!        Timsta(43)=127140.  !  ? WHY it's need, I don't no !
!        do k=1,-nst
!        print*,Timsta(k),k
!        enddo
!!        pause 1



! 
!!         print*,stHmax  !  Hmax here  st & 45
!!         pause 1




!!        goto 32   !  USIAL  'C', only for IC calc.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  create IC'0' TS !!!!!!!!!!!!!!!!!!!!!!! START
!        t=100.
!        s=100.

!        do i=1,in
!        do j=1,jn
!        ke=kp(i,j)
!        do k=1,ke
!        t(i,j,k)=Tst(45,k)
!        s(i,j,k)=Sst(45,k)
!        enddo
!        enddo
!        enddo


!        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  create IC TS !!!!!!!!!!!!!!!!!!!!!!! END






!!         goto 32  ! need to be out
!        do i=1,in
!        do j=1,jn
!        kk=0
!        do k=2,kn
!        if(kk.eq.1)goto 665
!        if(s(i,j,k).gt.90.)kk=1 
!        if(t(i,j,k).gt.90.)t(i,j,k)=t(i,j,k-1)
!        if(s(i,j,k).gt.90.)s(i,j,k)=s(i,j,k-1)
!665     continue       
!        enddo
!        if(kp(i,j).eq.kn)kp(i,j)=kn-1
!        enddo
!        enddo
!32      continue


!        open(1,file='t0')
!        open(2,file='s0')
!        write(1,*)t
!        write(2,*)s
!        close (1)
!        close (2)
!!        pause 1


!        open(211,file='Rlong')
!        open(222,file='Rlatt')
!        write(211,*)Rlong
!        write(222,*)Rlatt
!        close (211) 
!        close (222)

!        open(90,file='Slong')
!        open(91,file='Slatt')
!        write(90,*)Slong
!        write(91,*)Slatt
!        close (91) 
!        close (90)



!        open(21,file='Tst.dat')
!        write(21,*)Tst
!        close (21) 
!        open(22,file='Sst.dat')
!        write(22,*)Sst
!        close (22)


!        open(25,file='nstat.dat')
!        write(25,*)nst
!        close (25)



!        open(26,file='timsta.dat')
!        write(26,*)Timsta
!        close (26)
!        
!        

!        open(1,file='t0')
!        open(2,file='s0')
!        read(1,*)t
!        read(2,*)s
!        close (1)
!        close (2)
!        goto 69  !  control   !   RRR 20.08.2013
!        print*,kp(10,10),kp(100,100) 
!        do k=1,kn
!        print*,k
!        print*,t(10,10,k),t(100,100,k),Tst(45,k),zw(k)/100.
!        enddo
!        !pause 1
!69      continue

!        print*,' now   init=0 !!!'



!         endif !  if 1                                   &  end of init=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  end of init=0    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



!    next  ->  IC

!   then   ->  SOL 





!WWWWWWWWWWWWWWWWWWWWWWWWW init=2 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

!  now  ->    READING  all arrays


      read(49,66)iAi
      read(50,66)iAj
      read(51,66)iAns
66        format(80i5)
!          goto 44
!      do n=1,71
!      print*,iAi(n),iAj(n),n   !   RRR 20.08.2013
!      enddo
!      !pause 1
!44        continue


        open(211,file='Rlong')
        open(222,file='Rlatt')
        read(211,*)Rlong
        read(222,*)Rlatt
        close (211) 
        close (222)


        open(90,file='Slong')
        open(91,file='Slatt')
        read(90,*)Slong
        read(91,*)Slatt
        close (91) 
        close (90)
!        goto 58
!        do k=1,nst    !  RRR  18.08.2013
!        print*,Slong(k),Slatt(k),k
!        enddo
!        !pause 1
!58      continue






!        at=t
!        st=s


        open(21,file='Tst.dat')
        read(21,*)Tst
        close (21)
        

        open(22,file='Sst.dat')
        read(22,*)Sst
        close (22)

        
!        goto 67   ! control    !   RRR 20.08.2013
!        do n=1,nst
!        print*,n
!        do k=1,kn
!        print*,Tst(n,k),Sst(n,k),zw(k)/100.,k
!        enddo
!        !pause 1
!        enddo
!        !pause  3
!67      continue


        open(25,file='nstat.dat')   !   RRR 20.08.2013
        read(25,*)nst
        close (25)
!        print*,nst

        open(26,file='timsta.dat')
        read(26,*)timsta
        close (26)
!          goto 442
!      do n=1,71
!      print*,Timsta(n)/3600.,n   !   RRR 20.08.2013
!      enddo
!      !pause 21
!442        continue
!        pause 1

!WWWWWWWWWWWWWWWWWWWWWWWWW init=2 WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! + ASM  !!!!!!!!!!!!!

!      ku=kp
!      kv=kp
!      do i=2,in-1
!      do j=2,jn-1

!      if(it(i,j).eq.1) goto 438

!      ipa=kp(i,j)
!      ipw=kp(i-1,j)
!      ipn=kp(i,j-1)
!      ipe=kp(i+1,j)
!      ips=kp(i,j+1)

!      if(ipa.lT.ipw)then
!       ku(i,j)=ipa
!      else
!       ku(i,j)=ipw
!      endif

!      if(ipa.lT.ipe)then
!       ku(i+1,j)=ipa
!      else
!       ku(i+1,j)=ipe
!      endif

!      if(ipa.lT.ipn)then
!       kv(i,j)=ipa
!      else
!       kv(i,j)=ipn
!      endif

!      if(ipa.lT.ips)then
!       kv(i,j+1)=ipa
!      else
!       kv(i,j+1)=ips
!      endif

!438   continue
!      enddo
!      enddo

!!         open(80,file='ku')
!!         open(81,file='kv')

!!         write(80,117)ku
!!         write(81,117)kv

!313    continue   ! only out BK 

!!      ititititititititititititititititititititititititititititit
!!      KPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKPKP


!!!!!!!!!!!!!!!!!!!!!!!!!!  control t s  relief it  kp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!  control t s  relief it  kp  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!         goto 560   ! becouse IC homogenious from single  stat, later from row   measurements

!      if(iproc.eq.0) then
!       do i=1,in 
!       do j=1,jn
!       do k=1,kp(i,j) 
!c       if(it(i,j).eq.1)goto 220
!       if(kp(i,j).le.2)goto 229
!        if(k.le.kp(i,j))then
!       if(t(i,j,k).gt.20.or.t(i,j,k).lt.-2.5)print*,t(i,j,k),i,j,k,
!     &     kp(i,j),it(i,j)
!       if(t(i,j,k).gt.20.or.t(i,j,k).lt.-2.5)pause 871

!       if(s(i,j,k).gt.40.or.s(i,j,k).lt.8.)print*,s(i,j,k),t(i,j,k)
!     &,i,j,k,kp(i,j),it(i,j)
!       if(s(i,j,k).gt.40.or.s(i,j,k).lt.8.)pause 872
!                     endif
!       enddo

!       k1=kp(i,j)+1
!       do k=k1,kn
!       t(i,j,k)=t(i,j,k-1)
!       s(i,j,k)=s(i,j,k-1)
!      enddo

!229   continue

!      enddo
!      enddo
!      endif !iproc=0
!      
!560   continue

!       do j=kb2,ke2
!       do i=kb1,ke1 
!       k1=kp(i,j)+1
!       do k=k1,kn
!       t_p(i,j,k)=t_p(i,j,k-1)
!       s_p(i,j,k)=s_p(i,j,k-1)
!      enddo
!      enddo
!      enddo
!      
!!       pause 1
!!!!!!!!!!!!!!!!!!!!!!!!!!  control t s   and relief  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       ahzt=0.01
!       alb=5.e5          !!! eq BK
!       alb=1.e5          !!! ne BK
!       
!       
!      call incbe(tb,te,frf,pr,prt,ahz,ef,
!     &           tb0,te0,
!     &           wt,drdt,pare,alz,z0,rfcr,qb, !initial
!     &           l0t,lbt,l0b,lbb,l0e,lbe,ql,prb,ns,alb)



      call MPI_FINALIZE(ierr)
      
      return
      end  

ctext  *****************************************************************
       subroutine InputDecompNum (NPX,NPY,NUMPROC)
ctext         called by : conv3d
ctext         calls     : ijke_skip_str
ctext         prescription: 
ctext            read dimensions from 'decomp_mg.inp'
ctext         input:
ctext         output:     
ctext            NPX,NPY,NPZ
ctext         note:

         idfile=17
         open(idfile,file='decomp_mg.inp')

          call ijke_skip_str ( idfile )
         read(idfile,*) NPX,NPY
         if(NPX*NPY.ne.NUMPROC) then
          write(*,*) 'Error Number of Decomposition !!!'
          write(*,*) 'NPX*NPY = ', NPX*NPY
          write(*,*) 'NUMPROC =  ', NUMPROC
          stop
         endif

         close(idfile)

         return
         end
ctext  *****************************************************************
      subroutine ijke_skip_str ( 
     i                             idfile )
ctext         called by : ijke_inp_geo; ijke_inp_key; ijke_inp_prop;
ctext                     ijke_inp_mapf
ctext         calls     : 
ctext         prescription: 
ctext           pass commentary('$'or'#' begining) string in input files
ctext         input:
ctext            idfile  - identificator of file
ctext         output:     
ctext           pass string
ctext         note: 
      character*256 s

10       continue
          read(idfile,101) s
c           if ( .not.EOF(idfile) ) then 
           if ( (s(1:1).ne.'$').and.(s(1:1).ne.'#') ) then 
            backspace(17)
           else
            goto 10
           endif
c           endif

101      format(a)

      return
      end
ctext  *****************************************************************
      subroutine Decomp(n1,n2,n3)

      include 'mpif.h'
      COMMON / MPI311 / IPROC,NUMPROC, KB1,KE1,KB2,KE2
      COMMON / MPI312 / NPX,NPY,IPX,IPY
      COMMON / MPI313 / KBA1(0:16384),KEA1(0:16384),
     *                  KBA2(0:16384),KEA2(0:16384)
      COMMON / MPI314 / IPXA(0:16384),IPYA(0:16384)
      integer, dimension (MPI_STATUS_SIZE) :: istatus

c
c     get processor index in 2D grid (IPX, IPY)
c
      IPX = mod(IPROC, NPX)
      ll = (IPROC - IPX)/NPX
      IPY = mod(ll, NPY)
c
c     get local grid size (nlx, nly, nlz)
c
      nlx = (n1-1+NPX-1)/NPX
      nly = (n2-1+NPY-1)/NPY
c
c     get local grid indices (KBP1,KEP1,KBP2,KEP2,KBP3,KEP3)
c
      KB1 = nlx * IPX + 1
      KB2 = nly * IPY + 1

      if (IPX+1 .eq. NPX) then
         KE1 = n1
      else
!         KE1 = nlx*(IPX+1)+1
         KE1 = nlx*(IPX+1)
      endif
      if (IPY+1 .eq. NPY) then
         KE2 = n2
      else
!         KE2 = nly*(IPY+1)+1
         KE2 = nly*(IPY+1)
      endif
      
      if(IPROC.ne.0) then
        call MPI_SEND(KB1,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierr)
        call MPI_SEND(KE1,1,MPI_INTEGER,0,2,MPI_COMM_WORLD,ierr)
        call MPI_SEND(KB2,1,MPI_INTEGER,0,3,MPI_COMM_WORLD,ierr)
        call MPI_SEND(KE2,1,MPI_INTEGER,0,4,MPI_COMM_WORLD,ierr)
        call MPI_SEND(IPX,1,MPI_INTEGER,0,7,MPI_COMM_WORLD,ierr)
        call MPI_SEND(IPY,1,MPI_INTEGER,0,8,MPI_COMM_WORLD,ierr)
      endif                         

      if (IPROC.eq.0) then
         KBA1(0) = KB1
         KEA1(0) = KE1
         KBA2(0) = KB2
         KEA2(0) = KE2
         IPXA(0) = IPX
         IPYA(0) = IPY
         do 30 i = 1, NUMPROC-1
            call MPI_RECV(KBA1(i),1,MPI_INTEGER,i,1,MPI_COMM_WORLD,
     &                                                 istatus,ierr)
            call MPI_RECV(KEA1(i),1,MPI_INTEGER,i,2,MPI_COMM_WORLD,
     &                                                 istatus,ierr)
            call MPI_RECV(KBA2(i),1,MPI_INTEGER,i,3,MPI_COMM_WORLD,
     &                                                 istatus,ierr)
            call MPI_RECV(KEA2(i),1,MPI_INTEGER,i,4,MPI_COMM_WORLD,
     &                                                 istatus,ierr)
            call MPI_RECV(IPXA(i),1,MPI_INTEGER,i,7,MPI_COMM_WORLD,
     &                                                 istatus,ierr)
            call MPI_RECV(IPYA(i),1,MPI_INTEGER,i,8,MPI_COMM_WORLD,
     &                                                 istatus,ierr)
  30     continue
      endif

      call MPI_BCAST(KBA1(0),NUMPROC,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(KEA1(0),NUMPROC,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(KBA2(0),NUMPROC,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(KEA2(0),NUMPROC,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPXA(0),NUMPROC,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPYA(0),NUMPROC,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      
      return
      end
 
ctext  *****************************************************************


      subroutine cond(nstep,nprin,nk,al,al1,ctah,ct,cadv
     & ,eps,qwa,cu1,x1,energ,timestep,ns,taa,saa,itpar,iter,idd
     & ,kh,s,t,ta,u,v,w,tx,ty,thx,thy,ph,pht,ssp
     & ,hp,ahz,aday,p0d,p0md,cstep,str,h,relax,af
     & ,mon2,fday,hours,fminutes)

        include 'par.inc'
        include 'com.inc'
     
       integer  kh(in-1,jn-1)
      real s(in,jn,kn),t(in,jn,kn),fi(jn),ta(in,jn)
     *     ,thx(in,jn),thy(in,jn),x1(in),taa(kn),saa(kn)
     *     ,u(in,jn,kn),v(in,jn,kn),tax(in,jn),tay(in,jn)
     *     ,w(in,jn,kn),ph(in,jn),pht(in,jn),phtt(in,jn)
     *     ,energ(120000),ssp(in,jn,kn),ahz(in,jn,kn)
     *     ,hp(in,jn),af(in,jn),aa(kn),cc(kn)


      real*8 timestep(120000)

cc      parameter (fpi=3.1415926)

!!!!!!!!!!!!!!!!!!!!!!!!!  MOTO cond  special !!!!!!!!!!
      DOUBLE PRECISION c,d,e,f,ff,fff      ! IMP  for reading files.grd
      integer*2 a,b                        ! IMP  for reading files.grd
      real*4 g(4101,4051)                  ! IMP  for reading files.grd
      real h(in,jn)
!!!!!!!!!!!!!!!!!!!!!!!!!  MOTO cond  special !!!!!!!!!!

      parameter (fpi=3.1415926)

      open(8,file='par0.dat')

      read(8,*)init
      read(8,*)dt,vis,p0d,p0md
      read(8,*)rad,ql,al,al1,grad,omega,ddd,ddw
      read(8,*)nstep,nprin,nk,eps,niter          ! niter new
      read(8,*)u0,v0,t0,tx,ty,ctah,ct,cadv
      read(8,*)tan,tas,qwa,cu1
      read(8,*)cstep



!    print*,'lat0,latl,fn,fs',almin,almax,(90.-teta0),(90.-tetmax)





      dt=dt*60.
      aga=ddd/dt
      ag=ddw/dt
      et=vis/dt


      west=31.97  ! MOTO here
      east=34.75

      south=69.2                    
      pnorth=70.0



      almin=31.97  ! MOTO here
      almax=34.75

      tetmax=90-69.2                    
      tetmin=90-70.0


      dy=(tetmax-tetmin)/real(jn-1)*grad
      dx=(almax-almin)/real(in-1)*grad
    
!      print*,dx*rad/1.e5*sin(teta0*grad),dy*rad/1.e5


      do 2 i=1,in
         x(i)=dx*real(i-1)
 2    continue


      do 4 j=1,jn
         y(j)=90.-tetmax+(tetmax-tetmin)/(jnn)*(1.*j-0.5)

         TETA=tetmin*grad+real(j-1)*dy
         sn(J)=sin(TETA)
         cs(J)=cos(TETA)
  4   continue

!*******************************************************
!      print*,'Lx (km)',real(in-1)*dx*rad*cos(69.*grad)/100000.,
!     &       'Ly (km)',real(jn-1)*dy*rad/100000.
!      print*,'dx (km)',dx*rad*sn(1)/100000.,'dy (km)'
!     &         ,dy*rad/100000.
!      print*,'dx (km min)',dx*rad*sin(tetmin*grad)/100000.,'dy (km)'
!     &         ,dy*rad/100000.
!      print*,'dx (km max)',dx*rad*sin(tetmax*grad)/100000.,'dy (km)'
!     &         ,dy*rad/100000.

!!    pause 1




      do 5 j=1,jn
         acs(j)=cs(j)
         asn(j)=sn(j)
         yy=real(j-0.5)/real(jn-1.)*fpi
         fi(j)=cos(yy)
  5   continue

!      et=(ql*omega*acs(6))**0.5*ctah/z(kn)!!!!!! 04.08.02 delete c  

!      print*,'init = ', init

      if(init.eq.0)then

       goto 3
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       nm=11
       rewind 11

      open(nm,file='a0')

      read(nm,*)ns,mon2,fday,hours,fminutes,itpar,iter,idd 
!      read(nm,*)ns,iter
      read(nm,*)u
      read(nm,*)v
      read(nm,*)w
      read(nm,*)ph
      read(nm,*)ssp
      read(nm,*)t
      read(nm,*)s
3      continue

            ns=1
            u=0.
            v=0.
            w=0.
            ph=0.
            hour=0.      
            hours=0.


         else     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  init = 2    -> reading

       nm=11
       rewind 11

c        print *,'filename to read'
c        read(*,1100)file1
c        open(11,file=file1)
c1100  format(a3)

      open(nm,file='IC TSr')
      
!      open(nm,file='a21')

!       read(nm,*)ns,iter 
      read(nm,*)ns,mon2,fday,hours,fminutes
      
!      print*,ns,mon2,fday,hours,fminutes   ! standart
!      print*,itpar,iter,idd
!      pause 1
      
      read(nm,*)u
      read(nm,*)v
      read(nm,*)w
      read(nm,*)ph
      read(nm,*)ssp
      read(nm,*)t
      read(nm,*)s
      
!      print*,' n= ',t(50,10,9),t(50,10,15)
      
           goto 7      
            hour=0.      
            hours=0.
            u=0.
            v=0.
            w=0.
            ph=0.
7        continue


!      read(nm,*)ahz
!      read(nm,*)af
!       print*,' ns start =',ns
!       print*,'in cond'       
!       call prARRAYy(t,1)
!     pause 222

cc      write(nm,*)aff
cc      write(nm,*)aff
cc      write(nm,*)aff
cc      write(nm,*)aff


      endif


      close(8)
      close(2)
      close(3)
      close(4)
!      close(6)
      close(11)
      close(41)
      close(42)

2012  format(301f12.6)

!124    continue


      do j=1,jnn
      do i=1,inn
      iit=(it(i,jn-j)-2)*(it(i,jn-j)-1)*(it(i,jn-j)-3)
      if(iit.eq.0) kh(i,j)=99
      enddo
      enddo


ccc        ns=1  !!!!  TEMP !


      return
      end
c---------------------------------------------

ctext  *****************************************************************


      subroutine cond_p(nstep,nprin,nk,al,al1,ctah,ct,cadv
     & ,eps,qwa,cu1,x1,energ,timestep,ns,taa,saa,itpar,iter,idd
     & ,kh,s_p,t_p,ta,u_p,v_p,w_p,tx,ty,thx,thy,ph,pht,ssp_p
     & ,hp,aday,p0d,p0md,cstep,str,h,relax,af
     & ,mon2,fday,hours,fminutes)

      include 'par.inc'
      include 'com.inc'
      include 'mpif.h'        
      COMMON / MPI311 / IPROC,NUMPROC, KB1,KE1,KB2,KE2
     
       integer  kh(in-1,jn-1)
       integer(kind=MPI_OFFSET_KIND)::OFFSET      

      real s_p(KB1-1:KE1+1,KB2-1:KE2+1,kn)
     *,t_p(KB1-1:KE1+1,KB2-1:KE2+1,kn)
     *,fi(jn),ta(in,jn)
     *     ,thx(in,jn),thy(in,jn),x1(in),taa(kn),saa(kn)
     *     ,u_p(KB1-2:KE1+2,KB2-2:KE2+2,kn)
     *     ,v_p(KB1-2:KE1+2,KB2-2:KE2+2,kn)
     *     ,w_p(KB1-1:KE1+1,KB2-1:KE2+1,kn)
     *     ,tax(in,jn),tay(in,jn)
     *     ,ph(in,jn),pht(in,jn),phtt(in,jn)
     * ,energ(120000),ssp_p(KB1:KE1,KB2:KE2,kn)
     * ,hp(in,jn),af(in,jn),aa(kn),cc(kn)


      real*8 timestep(120000)

cc      parameter (fpi=3.1415926)

!!!!!!!!!!!!!!!!!!!!!!!!!  MOTO cond  special !!!!!!!!!!
      DOUBLE PRECISION c,d,e,f,ff,fff      ! IMP  for reading files.grd
      integer*2 a,b                        ! IMP  for reading files.grd
      real*4 g(4101,4051)                  ! IMP  for reading files.grd
      real h(in,jn)
!!!!!!!!!!!!!!!!!!!!!!!!!  MOTO cond  special !!!!!!!!!!

      parameter (fpi=3.1415926)

      open(8,file='par0.dat')

      read(8,*)init
      read(8,*)dt,vis,p0d,p0md
      read(8,*)rad,ql,al,al1,grad,omega,ddd,ddw
      read(8,*)nstep,nprin,nk,eps,niter          ! niter new
      read(8,*)u0,v0,t0,tx,ty,ctah,ct,cadv
      read(8,*)tan,tas,qwa,cu1
      read(8,*)cstep
      
      close(8)

      dt=dt*60.
      aga=ddd/dt
      ag=ddw/dt
      et=vis/dt


      west=31.97  ! MOTO here
      east=34.75

      south=69.2                    
      pnorth=70.0



      almin=31.97  ! MOTO here
      almax=34.75

      tetmax=90-69.2                    
      tetmin=90-70.0


      dy=(tetmax-tetmin)/real(jn-1)*grad
      dx=(almax-almin)/real(in-1)*grad

      do 2 i=1,in
         x(i)=dx*real(i-1)
 2    continue


      do 4 j=1,jn
         y(j)=90.-tetmax+(tetmax-tetmin)/(jnn)*(1.*j-0.5)
         TETA=tetmin*grad+real(j-1)*dy
         sn(J)=sin(TETA)
         cs(J)=cos(TETA)
  4   continue

*******************************************************
!      print*,'Lx (km)',real(in-1)*dx*rad*cos(69.*grad)/100000.,
!     &       'Ly (km)',real(jn-1)*dy*rad/100000.
!      print*,'dx (km)',dx*rad*sn(1)/100000.,'dy (km)'
!     &         ,dy*rad/100000.
!      print*,'dx (km min)',dx*rad*sin(tetmin*grad)/100000.,'dy (km)'
!     &         ,dy*rad/100000.
!      print*,'dx (km max)',dx*rad*sin(tetmax*grad)/100000.,'dy (km)'
!     &         ,dy*rad/100000.



      do 5 j=1,jn
         acs(j)=cs(j)
         asn(j)=sn(j)
         yy=real(j-0.5)/real(jn-1.)*fpi
         fi(j)=cos(yy)
  5   continue

!      et=(ql*omega*acs(6))**0.5*ctah/z(kn)!!!!!! 04.08.02 delete c  


      if(init.eq.0)then

!       goto 3
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       nm=11
!       rewind 11

!      open(nm,file='a0')

!      read(nm,*)ns,mon2,fday,hours,fminutes,itpar,iter,idd 
!!      read(nm,*)ns,iter
!      read(nm,*)u
!      read(nm,*)v
!      read(nm,*)w
!      read(nm,*)ph
!      read(nm,*)ssp
!      read(nm,*)t
!      read(nm,*)s
!3      continue

            ns=1
            u_p=0.
            v_p=0.
            w_p=0.
            ph=0.
            hour=0.      
            hours=0.


         else     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  init = 2    -> reading

       nm=11
       rewind 11

c        print *,'filename to read'
c        read(*,1100)file1
c        open(11,file=file1)
c1100  format(a3)

      open(nm,file='IC TSr')
      
!      open(nm,file='a21')

!       read(nm,*)ns,iter 
      read(nm,*)ns,mon2,fday,hours,fminutes
      
      close(nm)
      
!      print*,ns,mon2,fday,hours,fminutes   ! standart
!      print*,itpar,iter,idd
!      pause 1
      
      OFFSET=0
      call pread(u_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
      call compare(u_p(KB1:KE1,KB2:KE2,1:km),u,im,jm,km,'u pread')
      
      OFFSET=OFFSET+in*jn*kn*4
      call pread(v_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
      call compare(v_p(KB1:KE1,KB2:KE2,1:km),v,im,jm,km,'v pread')

      OFFSET=OFFSET+in*jn*kn*4
      call pread(w_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
      call compare(w_p(KB1:KE1,KB2:KE2,1:km),w,im,jm,km,'w pread')

      OFFSET=OFFSET+in*jn*kn*4
      call pread(ssp_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
      call compare(ssp_p(KB1:KE1,KB2:KE2,1:km),ssp,im,jm,km,'ssp pread')
      
      OFFSET=OFFSET+in*jn*kn*4
      call pread(t_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
      call compare(t_p(KB1:KE1,KB2:KE2,1:km),t,im,jm,km,'t pread')
      
      OFFSET=OFFSET+in*jn*kn*4
      call pread(s_p(KB1:KE1,KB2:KE2,1:kn),in,jn,kn,'IC TCr_p',OFFSET)
      call compare(s_p(KB1:KE1,KB2:KE2,1:km),s,im,jm,km,'s pread')
      
      open(41,file='ph')
      read(41,*)ph
      close(41)
      
!      print*,' n= ',t(50,10,9),t(50,10,15)
      
!           goto 7      
!            hour=0.      
!            hours=0.
!            u=0.
!            v=0.
!            w=0.
!            ph=0.
!7        continue

      endif

!      close(2)
!      close(3)
!      close(4)
!      close(11)
!      close(42)

!2012  format(301f12.6)

!124    continue


      do j=1,jnn
      do i=1,inn
      iit=(it(i,jn-j)-2)*(it(i,jn-j)-1)*(it(i,jn-j)-3)
      if(iit.eq.0) kh(i,j)=99
      enddo
      enddo


ccc        ns=1  !!!!  TEMP !


      return
      end
c---------------------------------------------

c---------------------------------------------
      subroutine pwrite(ar1,im,jm,km,fname,OFFSET)
      include 'par.inc'
      include 'mpif.h'       
c        include 'com.inc'
      COMMON / MPI311 / IPROC,NUMPROC, KB1,KE1,KB2,KE2
      
      real,dimension(KB1:KE1,KB2:KE2,1:km)::ar1
      CHARACTER*(*) fname
      integer(kind=MPI_OFFSET_KIND)::OFFSET
      
!      parameter (NBUFSIZE=3) 
!      real*8 buf(NBUFSIZE) 
!      integer(kind=MPI_OFFSET_KIND) disp 
    
       integer NDIMS
       integer ARRAY_OF_SIZES(3)
       integer ARRAY_OF_SUBSIZES(3)
       integer ARRAY_OF_STARTS(3)
       integer count
 
       
      call MPI_FILE_OPEN(MPI_COMM_WORLD, fname
     *,MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr) 
      
       NDIMS=3
       ARRAY_OF_SIZES(1)=im
       ARRAY_OF_SIZES(2)=jm
       ARRAY_OF_SIZES(3)=km

       ARRAY_OF_SUBSIZES(1)=KE1-KB1+1
       ARRAY_OF_SUBSIZES(2)=KE2-KB2+1
       ARRAY_OF_SUBSIZES(3)=km

       ARRAY_OF_STARTS(1)=KB1-1
       ARRAY_OF_STARTS(2)=KB2-1
       ARRAY_OF_STARTS(3)=0

       call MPI_TYPE_CREATE_SUBARRAY(NDIMS,ARRAY_OF_SIZES,
     *ARRAY_OF_SUBSIZES,ARRAY_OF_STARTS,MPI_ORDER_FORTRAN,
     *MPI_REAL, filetype_scalar, ierr)
      call MPI_TYPE_COMMIT(filetype_scalar, ierr)

      call MPI_FILE_SET_VIEW(fh,OFFSET,MPI_REAL,filetype_scalar,
     *'native',MPI_INFO_NULL,ierr)
     
      count=(KE1-KB1+1)*(KE2-KB2+1)*km

      call MPI_FILE_WRITE_ALL(fh,ar1(KB1:KE1,KB2:KE2,1:km), 
     *count,MPI_REAL,MPI_STATUS_IGNORE,ierr)



ccccc

!      OFFSET=IPROC*4
!      call MPI_FILE_SET_VIEW(fh,OFFSET,MPI_INTEGER4,MPI_INTEGER4,
!     *'native',MPI_INFO_NULL,ierr)
!      count=1
!      call MPI_FILE_WRITE_ALL(fh,IPROC, 
!     *1,MPI_INTEGER4,MPI_STATUS_IGNORE,ierr)      
            
            
     
!      do i = 0, NBUFSIZE 
!        buf(i) = IPROC * NBUFSIZE + i 
!      enddo 
!      
!      disp = IPROC * NBUFSIZE * 8 
!      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER,
!     &                      MPI_DOUBLE, 'native',
!     &                      MPI_INFO_NULL, ierr) 
!      call MPI_FILE_WRITE_ALL(fh, buf, NBUFSIZE, MPI_DOUBLE,
!     &                  MPI_STATUS_IGNORE, ierr)      
     
      call MPI_FILE_CLOSE(fh, ierr)
      
      return
      end

c---------------------------------------------
      subroutine pread(ar1,im,jm,km,fname,OFFSET)
      include 'par.inc'
      include 'mpif.h'       
c        include 'com.inc'
      COMMON / MPI311 / IPROC,NUMPROC, KB1,KE1,KB2,KE2
      
      real,dimension(KB1:KE1,KB2:KE2,1:km)::ar1
      CHARACTER*(*) fname
      integer(kind=MPI_OFFSET_KIND)::OFFSET
      
!      parameter (NBUFSIZE=3) 
!      real*8 buf(NBUFSIZE) 
!      integer(kind=MPI_OFFSET_KIND) disp 
    
       integer NDIMS
       integer ARRAY_OF_SIZES(3)
       integer ARRAY_OF_SUBSIZES(3)
       integer ARRAY_OF_STARTS(3)
       integer count
 
       
      call MPI_FILE_OPEN(MPI_COMM_WORLD, fname
     *,MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr) 
      
       NDIMS=3
       ARRAY_OF_SIZES(1)=im
       ARRAY_OF_SIZES(2)=jm
       ARRAY_OF_SIZES(3)=km

       ARRAY_OF_SUBSIZES(1)=KE1-KB1+1
       ARRAY_OF_SUBSIZES(2)=KE2-KB2+1
       ARRAY_OF_SUBSIZES(3)=km

       ARRAY_OF_STARTS(1)=KB1-1
       ARRAY_OF_STARTS(2)=KB2-1
       ARRAY_OF_STARTS(3)=0

       call MPI_TYPE_CREATE_SUBARRAY(NDIMS,ARRAY_OF_SIZES,
     *ARRAY_OF_SUBSIZES,ARRAY_OF_STARTS,MPI_ORDER_FORTRAN,
     *MPI_REAL, filetype_scalar, ierr)
      call MPI_TYPE_COMMIT(filetype_scalar, ierr)

      call MPI_FILE_SET_VIEW(fh,OFFSET,MPI_REAL,filetype_scalar,
     *'native',MPI_INFO_NULL,ierr)
     
      count=(KE1-KB1+1)*(KE2-KB2+1)*km

      call MPI_FILE_READ_ALL(fh,ar1(KB1:KE1,KB2:KE2,1:km), 
     *count,MPI_REAL,MPI_STATUS_IGNORE,ierr)

     
      call MPI_FILE_CLOSE(fh, ierr)
      
      return
      end
                            
