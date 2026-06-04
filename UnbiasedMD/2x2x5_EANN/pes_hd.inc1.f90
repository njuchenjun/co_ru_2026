     integer(kind=intype) :: yy,mm,dd,hh,mi,ss,ms     ! for print the system time
     integer(kind=intype) :: nloop,num_der,table_coor
     integer(kind=intype) :: start_force,start_wb,start_init,table_grid
     integer(kind=intype) :: numpoint,nsurf,nt,nv,ntraw,nvraw,nbatch,maxbatch
     integer(kind=intype) :: maxnumatom,maxneff,maxnumforce,maxtpoint,maxnpoint,dimjac
     integer(kind=intype) :: mnl,mnhid,mnw,nkpoint,outputneuron,wave_allnw,nn_allnw,allnw,ncycle
!----------------------------------real-------------------------------------------------
     real(kind=retype) :: mu,biase_rmse,force_perindex
!---------------------------------integer allocatable array---------------------------------------------------
     integer(kind=intype),allocatable :: addpoint_t(:),addpoint_v(:),addfpoint_t(:),addfpoint_v(:),numth(:)
     integer(kind=intype),allocatable :: ntpoint(:),nvpoint(:)
     integer(kind=intype),allocatable :: numatom(:),npoint(:)
     integer(kind=intype),allocatable :: nl(:,:),nhid(:),nw(:),elmnw(:),hidneu(:),last_neu(:)
     integer(kind=intype),allocatable :: index_table(:,:),index_ele(:,:)
     integer(kind=intype),allocatable :: neff(:),numforce(:)
     integer(kind=intype),allocatable :: index_surf_t(:),index_surf_v(:),index_num(:)
     integer(kind=intype),allocatable :: batch(:),dimbatch(:),index_batch(:),sumbatch(:),countbatch(:)
     integer(kind=intype),allocatable :: ipiv(:)
!----------------------------------real allocatable array-------------------------------------------------------
     real(kind=retype),allocatable :: w(:,:,:,:,:),w_back(:,:,:,:,:),w_save(:,:,:,:,:),perindex(:)
     real(kind=retype),allocatable :: dw(:),dw_back(:)
     real(kind=retype),allocatable :: y(:,:,:),abforce(:,:,:,:,:),error_weight(:)
     real(kind=retype),allocatable :: yt(:,:,:),yv(:,:,:),abforcet(:,:,:,:,:),abforcev(:,:,:,:,:),ewt(:),ewv(:)
     real(kind=retype),allocatable :: yft(:),yfv(:),nnyft(:),nnyfv(:)
     real(kind=retype),allocatable :: z(:,:),totet(:,:,:),totev(:,:,:)
     real(kind=retype),allocatable :: et(:),ev(:),evr(:),etr(:)
     real(kind=retype),allocatable :: forcet(:,:,:,:),forcev(:,:,:,:)
     real(kind=retype),allocatable :: JAC(:,:),JTJ(:,:),JTJM(:,:)
     real(kind=retype),allocatable :: dfdw(:,:,:),dedw(:,:,:,:),dode(:,:,:),dgdw(:,:,:)
     real(kind=retype),allocatable :: errt(:),errv(:),errtr(:),errvr(:),errall(:),tot_err(:)
     real(kind=retype),allocatable :: th(:,:,:)
     real(kind=retype),allocatable :: work(:)
     real(kind=retype),allocatable :: storeyt(:,:),jpvt(:),storeth(:,:),elm_yt(:,:)
     real(kind=retype),allocatable :: allcoor(:,:,:),allcoort(:,:,:),allcoorv(:,:,:)
!----------------------------------character allocatable array----------------------------------------------
     character(len=5),allocatable :: atom(:,:),atomtype(:)
!---------------------------------logical------------------------------------------------
     logical,allocatable :: m00(:,:)
!
!
!------------------------------------------wave parameter----------------------------------
     real(kind=retype),allocatable :: wf(:,:),dwfdcoor(:,:,:,:)
     integer(kind=intype) :: norbit
     integer(kind=intype) :: maxnumtype,maxnpara,atomwave,maxnwave
     integer(kind=intype) :: ipsin
!-----------------------------integer array-------------------------------------------
!-----------------------------integer allocatable array--------------------------------------------------
     integer(kind=intype),allocatable :: nimage(:,:),length(:)
     integer(kind=intype),allocatable :: nipsin(:),npara(:),index_power(:,:,:),nwave(:)
     integer(kind=intype),allocatable :: index_type(:),index_orbit(:)
     integer(kind=intype),allocatable :: index_para(:)
!------------------------------real array--------------------------------------------------
     real(kind=retype),allocatable :: rc(:)
     real(kind=retype),allocatable :: inta(:,:),rs(:,:)
     real(kind=retype),allocatable :: factorial(:),factor_wave(:,:,:),maxwf(:,:)
     real(kind=retype) :: scal(3,3)
!    real(kind=retype),allocatable :: scalmatrix(:,:,:)
     real(kind=retype),allocatable :: weight_wave(:,:),weight_wave_div_maxwf(:,:)

contains

subroutine init_pes1 !!@@@@ rename
    implicit none

    call readinput
    call allocate_all
    call get_index

    ! calculate one point
    nsurf=1
    numatom(1)=maxnumatom
    neff(1)=maxnumatom
    npoint(1)=1

    call nt_nv
    call allocate_nt_nv
    call readnet

    return
end subroutine

subroutine evaluate_pes1(ifrac,iforce,natoms,symbols,coord,e,f) !!@@@@ rename
      implicit none
      integer,intent(in) :: ifrac,iforce,natoms
      character(kind=1,len=2) :: symbols(natoms)
      real(kind=8),intent(in) :: coord(3,natoms)
      real(kind=8),intent(out) :: e,f(3*natoms)
      real(kind=8) :: cart(3,natoms),frac(3,natoms)
      integer :: i,k,p

      e=0.d0
      f=0.d0

      if (ifrac.eq.0) then
         cart=coord
        !frac=matmul(inv_matrix,coord)
      else
        !frac=coord
         cart=matmul(scal,coord)
      endif

      ! transfer inputs to ea-nn
      numatom(1)=natoms
      neff(1)=natoms
      npoint(1)=1
      do i=1,natoms
      atom(i,1)=symbols(i)
      allcoor(1:3,i,1)=cart(1:3,i)
      enddo

      call readcoor

      do p=1,maxnumtype
         do i=1,norbit
            k=(i-1)*maxnumtype+1
            weight_wave(k:k+maxnumtype,p)=weight_wave_div_maxwf(k:k+maxnumtype,p)*maxwf(i,p)
         end do
      end do

      call random_divide

      ! evaluate ea-nn
      z(0,:)=1.d0 !! this is important
      if (iforce.eq.0) then
         call feedforward
      else
         call get_force
      end if

      ! get results from ea-nn
      e=totet(1,1,1)
      if (start_force.ne.1) then
         f=0.d0
      else
         do i=1,natoms
            f(i*3-2:i*3)=forcet(i*3-2:i*3,1,1,1)
         enddo
      endif

      return
end subroutine

subroutine deallocate_pes1 !!@@@@ rename
     implicit none
       deallocate(nwave)
       deallocate(inta)
       deallocate(rs)
       deallocate(factorial)
       deallocate(factor_wave)
       deallocate(nipsin)
       deallocate(npara)
       deallocate(atom)
       deallocate(atomtype)
       deallocate(npoint)
       deallocate(nl)
       deallocate(nhid)
       deallocate(nw)
       deallocate(elmnw)
       deallocate(w)
       deallocate(w_back)
       deallocate(w_save)
       deallocate(neff)
       deallocate(addpoint_t)
!      deallocate(addpoint_v)
       deallocate(addfpoint_t)
!      deallocate(addfpoint_v)
       deallocate(numth)
       deallocate(ntpoint)
!      deallocate(nvpoint)
       deallocate(hidneu)
       deallocate(last_neu)
       deallocate(numatom) 
       deallocate(numforce) 
       deallocate(batch)
       deallocate(dimbatch)
       deallocate(index_batch)
       deallocate(sumbatch)
       deallocate(countbatch)
       deallocate(m00)
       deallocate(index_table)
       deallocate(dode)
       deallocate(dedw)
       deallocate(dfdw)
       deallocate(dgdw)
       deallocate(z)
       deallocate(wf)
       deallocate(error_weight)
       deallocate(y)
       deallocate(index_ele)
       deallocate(dw)
       deallocate(dw_back)
       deallocate(allcoor)
       if(start_force.eq.1) then
         deallocate(dwfdcoor)
         deallocate(abforce)
         deallocate(abforcet)
!        deallocate(abforcev)
         deallocate(forcet)
!        deallocate(forcev)
       end if
!      deallocate(JTJ)
!      deallocate(JTJM)
       deallocate(errtr)
!      deallocate(errvr)
       deallocate(errt)
!      deallocate(errv)
       deallocate(errall)
       deallocate(tot_err)
       deallocate(index_para)
       deallocate(index_type)
       deallocate(index_orbit)
       deallocate(weight_wave)
       deallocate(weight_wave_div_maxwf)
       deallocate(index_power)
       deallocate(maxwf)
       deallocate(ipiv)
       deallocate(index_surf_t)
!      deallocate(index_surf_v)
       deallocate(index_num)
       deallocate(yt)
!      deallocate(yv)
       deallocate(allcoort)
!      deallocate(allcoorv)
       deallocate(totet)
!      deallocate(totev)
       deallocate(yft)
!      deallocate(yfv)
       deallocate(nnyft)
!      deallocate(nnyfv)
       deallocate(et)
!      deallocate(ev)
       deallocate(ewt)     
!      deallocate(ewv)
       deallocate(etr)
!      deallocate(evr)
       deallocate(nimage)
       deallocate(length)
!      deallocate(scalmatrix)
       deallocate(perindex)
       deallocate(rc)
       if(start_wb.eq.1) then
         deallocate(th)
         deallocate(storeth)
         deallocate(elm_yt)
       end if
     return
end subroutine

subroutine allocate_all
     implicit none
       allocate(m00(maxnpoint,nsurf))
       allocate(index_table(5,allnw))
       allocate(dode(mnl,outputneuron,mnhid+1))
       allocate(dedw(mnl,mnl,mnhid+1,mnhid+1))
       allocate(dfdw(mnw,norbit,outputneuron))
       allocate(dgdw(norbit,norbit,outputneuron))
       allocate(wf(maxnpara,atomwave))
       allocate(z(0:mnl,0:mnhid+1))
       allocate(error_weight(numpoint))
       allocate(y(outputneuron,nkpoint,numpoint))
       allocate(index_ele(maxnumatom,nsurf))
       allocate(dw(allnw))
       allocate(dw_back(allnw))
       allocate(allcoor(atomdim,maxnumatom,numpoint))
       if(start_force==1) then
         allocate(abforce(atomdim,maxneff,outputneuron,nkpoint,numpoint))
         allocate(dwfdcoor(atomdim,maxnumatom,maxnpara,atomwave))
       end if
!      allocate(JTJ(allnw,allnw))
!      allocate(JTJM(allnw,allnw))
       allocate(ipiv(allnw))
       allocate(errtr(0:nloop))
!      allocate(errvr(0:nloop))
       allocate(errt(0:nloop))
!      allocate(errv(0:nloop))
       allocate(errall(0:nloop))
       allocate(tot_err(0:nloop))
       allocate(index_para(norbit))
       allocate(index_type(atomwave))
       allocate(index_orbit(atomwave))
       allocate(weight_wave(atomwave+1,maxnumtype)) !! @@ 2021/8/11
       allocate(weight_wave_div_maxwf(atomwave+1,maxnumtype)) !! @@ 2021/8/11
       allocate(index_power(3,maxnpara,norbit))
       allocate(maxwf(norbit,maxnumtype))
!      allocate(scalmatrix(atomdim,atomdim,nsurf))
       allocate(nimage(atomdim,nsurf))
       allocate(length(nsurf))
     return
end subroutine

subroutine allocate_nn
     implicit none
       allocate(numatom(nsurf)) 
       allocate(numforce(nsurf)) 
       allocate(atom(maxnumatom,nsurf))
       allocate(atomtype(maxnumtype))
       allocate(nl(0:mnhid+1,maxnumtype))
       allocate(nhid(maxnumtype))
       allocate(nw(maxnumtype))
       allocate(elmnw(maxnumtype))
       allocate(w(0:mnl,1:mnl,mnhid+1,maxnumtype,nkpoint))
       allocate(w_back(0:mnl,1:mnl,mnhid+1,maxnumtype,nkpoint))
       allocate(w_save(0:mnl,1:mnl,mnhid+1,maxnumtype,nkpoint))
       allocate(neff(nsurf))
       allocate(addpoint_t(nsurf))
!      allocate(addpoint_v(nsurf))
       allocate(addfpoint_t(nsurf))
!      allocate(addfpoint_v(nsurf))
       allocate(numth(0:nsurf))
       allocate(ntpoint(nsurf))
!      allocate(nvpoint(nsurf))
       allocate(hidneu(0:maxnumtype))
       allocate(last_neu(0:maxnumtype))
       allocate(batch(nbatch))
       allocate(dimbatch(nbatch))
       allocate(index_batch(maxtpoint))
       allocate(sumbatch(0:nbatch))
       allocate(countbatch(0:nbatch))
     return
end subroutine
 
subroutine allocate_nt_nv
     implicit none
       allocate(index_surf_t(nt))
!      allocate(index_surf_v(nv))
       allocate(index_num(0:nt))
       allocate(yt(outputneuron,nkpoint,nt))
!      allocate(yv(outputneuron,nkpoint,nv))
       allocate(allcoort(atomdim,maxnumatom,nt))
!      allocate(allcoorv(atomdim,maxnumatom,nv))
       if(start_force==1) then
         allocate(abforcet(atomdim,maxneff,outputneuron,nkpoint,nt))
!        allocate(abforcev(atomdim,maxneff,outputneuron,nkpoint,nv))
         allocate(forcet(maxnumforce,outputneuron,nkpoint,nt))
!        allocate(forcev(maxnumforce,outputneuron,nkpoint,nv))
       end if
       allocate(totet(outputneuron,nkpoint,nt))
!      allocate(totev(outputneuron,nkpoint,nv))
       allocate(yft(ntraw))
!      allocate(yfv(nvraw))
       allocate(nnyft(ntraw))
!      allocate(nnyfv(nvraw))
       allocate(et(ntraw))
!      allocate(ev(nvraw))
       allocate(ewt(nt))
!      allocate(ewv(nv)) 
       allocate(etr(nt))
!      allocate(evr(nv))
       if(start_wb==1) then
         allocate(th(numth(nsurf),last_neu(maxnumtype),nkpoint))
         allocate(storeth(numth(nsurf),last_neu(maxnumtype)))
         allocate(elm_yt(numth(nsurf),outputneuron))
       end if
     return
end subroutine

subroutine allocate_wave
     implicit none
       allocate(nwave(maxnumtype))
       allocate(inta(maxnumtype,0:maxnwave))
       allocate(rs(maxnumtype,0:maxnwave))
       allocate(factorial(0:ipsin))
       allocate(factor_wave(0:ipsin,0:ipsin,0:ipsin))
       allocate(nipsin(0:ipsin))
       allocate(npara(0:ipsin))
       allocate(npoint(nsurf))
       allocate(perindex(0:nsurf)) 
       allocate(rc(maxnumtype))
     return
end subroutine

subroutine readinput
      implicit none
      integer(kind=intype) :: i,j,k,l,ntype
      real(kind=retype) :: dier_rs,hyperpara

      open(870711,file=trim(pesdir)//'/scalmatrix',status='old',action='read')
      read(870711,*) scal(1:3,1)
      read(870711,*) scal(1:3,2)
      read(870711,*) scal(1:3,3)
      close(870711)

      open(880725,file=trim(pesdir)//'/eann.inp',status='old',action='read')
      read(880725,*) start_force
      read(880725,*) start_wb
      read(880725,*) start_init
      read(880725,*) table_coor
      read(880725,*) table_grid
      read(880725,*) biase_rmse
      read(880725,*) maxnumtype
      read(880725,*) maxnumatom,maxneff
      read(880725,*) nsurf; nsurf=1 !! nsurf
      read(880725,*) ipsin
      read(880725,*) maxnwave
      call allocate_wave
      read(880725,*) nwave(1:maxnumtype)
      read(880725,*) ! cutnum,ncell
      read(880725,*) ncycle
      read(880725,*) nloop; nloop=1
      read(880725,*) nbatch
      read(880725,*) numpoint,maxnpoint
      read(880725,*) maxtpoint ; maxtpoint=numpoint
      read(880725,*) perindex(0),force_perindex
      read(880725,*) rc(1:maxnumtype)
      read(880725,*) mnl
      read(880725,*) mnhid
      read(880725,*) nkpoint
      read(880725,*) outputneuron
!     read(880725,*) num_der
      call allocate_nn
!-----------read the wave_nn structure----------------
       do i=0,ipsin
         nipsin(i)=i
       end do
       factorial(0)=1.d0
       do i=1,ipsin
         factorial(i)=factorial(i-1)*dble(i)
       end do
       npara=0 
       do i=0,ipsin
         do j=0,nipsin(i)
           do k=0,nipsin(i)-j
             npara(i)=npara(i)+1
             l=nipsin(i)-j-k
             factor_wave(k,j,i)=factorial(nipsin(i))/(factorial(j)*factorial(k)*factorial(l)) 
           end do
         end do
       end do
       maxnpara=maxval(npara)
       nw=0
       norbit=(1+maxnwave)*(1+ipsin)
       nl(0,:)=norbit
       hidneu=0
       last_neu=0
       mnw=0
       allnw=0
       if(table_grid.eq.0) then
         do i=1,maxnumtype
           read(880725,*) atomtype(i),nhid(i),nl(1:nhid(i),i)
           nl(nhid(i)+1,i)=outputneuron
           do k=1,nhid(i)+1
             nw(i)=nw(i)+(nl(k-1,i)+1)*nl(k,i)
           end do
           if(mnw<nw(i)) mnw=nw(i)
           elmnw(i)=nw(i)
           if(start_wb.eq.1) elmnw(i)=nw(i)-(nl(nhid(i),i)+1)*nl(nhid(i)+1,i)
           nw(i)=elmnw(i)
           hidneu(i)=hidneu(i-1)+nw(i)
           last_neu(i)=last_neu(i-1)+nl(nhid(i),i)
           do j=0,nwave(i)
             read(880725,*) l,inta(i,j),rs(i,j)
           end do
         end do
       else
         rs=0d0
         do i=1,maxnumtype
           read(880725,*) atomtype(i),nhid(i),nl(1:nhid(i),i),hyperpara
           dier_rs=rc(i)/(nwave(i)+1.d0/3.d0)
           nl(nhid(i)+1,i)=outputneuron
           do k=1,nhid(i)+1
             nw(i)=nw(i)+(nl(k-1,i)+1)*nl(k,i)
           end do
           if(mnw<nw(i)) mnw=nw(i)
           elmnw(i)=nw(i)
           if(start_wb.eq.1) elmnw(i)=nw(i)-(nl(nhid(i),i)+1)*nl(nhid(i)+1,i)
           nw(i)=elmnw(i)
           hidneu(i)=hidneu(i-1)+nw(i)
           last_neu(i)=last_neu(i-1)+nl(nhid(i),i)
           inta(i,:)=hyperpara/dier_rs**2
           do j=1,nwave(i)
             rs(i,j)=rs(i,j-1)+dier_rs
           end do
         end do
       end if
       wave_allnw=0
       do i=0,ipsin
         do k=0,maxnwave
           do ntype=1,maxnumtype
             wave_allnw=wave_allnw+1
           end do
         end do
       end do
       atomwave=wave_allnw
       wave_allnw=atomwave*maxnumtype
       allnw=hidneu(maxnumtype)
       nn_allnw=allnw*nkpoint
       allnw=nn_allnw+wave_allnw
       if(start_wb.eq.1) allnw=allnw+maxnumtype*outputneuron*nkpoint
       if(start_force.eq.1) maxnumforce=maxneff*atomdim
       i=0
       l=floor(dble(maxtpoint/nbatch))
       countbatch=0
       do j=1,nbatch-1
         batch(j)=l
         index_batch(i+1:i+l)=j
         i=i+l
         countbatch(j)=i
       end do 
       countbatch=maxtpoint
       index_batch(i+1:maxtpoint)=nbatch
       batch(nbatch)=maxtpoint-i
       close(880725)
     return
end subroutine

subroutine readnet
      implicit none
      character(kind=1,len=3) :: cc,ckpoint
      integer(kind=intype) :: ilayer,jneuron,jp,i,p,k,ikpoint

      write(cc,'(i3.3)') pesid

      do ikpoint=1,nkpoint
         write(ckpoint,'(i3)') ikpoint
         do p=1,maxnumtype
            open(123,file=trim(pesdir)//'/pes'//cc//'-weight-'//trim(atomtype(p))//trim(adjustl(ckpoint))//'.txt',status='old',action='read')
            read(123,*) !! NN structure
            read(123,*) !! RMSE

            do ilayer=1,nhid(p)+1
               ! weights
               do jneuron=1,nl(ilayer,p)
               do jp=1,nl(ilayer-1,p)
                  read(123,*) w(jp,jneuron,ilayer,p,ikpoint)
               enddo
               enddo
               ! blank line
               ! biases
               do jneuron=1,nl(ilayer,p)
                  read(123,*) w(0,jneuron,ilayer,p,ikpoint)
               enddo
               ! blank line
            end do
            ! following are transfer functions
            close(123)
         end do
      end do

      do p=1,maxnumtype
         open(123,file=trim(pesdir)//'/pes'//cc//'-wave-'//trim(adjustl(atomtype(p)))//'.txt',status='old',action='read')
         do i=1,norbit
            k=(i-1)*maxnumtype+1
            read(123,'(e25.18)') weight_wave_div_maxwf(k:k+maxnumtype,p)
!           weight_wave(k:k+maxnumtype,p)=weight_wave_div_maxwf(k:k+maxnumtype,p)*maxwf(i,p)
         end do
      end do
      close(123)

!     write(*,*) "a: ",weight_wave_div_maxwf(1:3,2)

      return
end subroutine 

subroutine readcoor
     implicit none
     integer(kind=intype) :: isurf,i,j,n,natom,k,num,ntype
     real(kind=retype) :: matrix(atomdim,atomdim),tmp(norbit,maxnumatom),maxrc

     error_weight=1d0
     num=0
    !npoint=0
    !numatom=0
    !neff=0
     maxwf=1.d0
     tmp=0.d0
     maxrc=maxval(rc)

     do isurf=1,nsurf

         call inverse_matrix(scal(:,:),matrix)
         nimage(1,isurf)=nint(maxrc/scal(1,1))
         nimage(2,isurf)=nint(maxrc/scal(2,2))
         nimage(3,isurf)=nint(maxrc/scal(3,3))
         length(isurf)=(2*nimage(1,isurf)+1)*(2*nimage(2,isurf)+1)*(2*nimage(3,isurf)+1)
      !---------------------------save the index_ele for element---------------------------------
         do i=1,numatom(isurf)
           do j=1,maxnumtype
             if(atom(i,isurf)==atomtype(j)) index_ele(i,isurf)=j
           end do
         end do
      !--------------------------------------------------------------------------------------------

           num=num+1

           if(table_coor==0) allcoor(:,:,num)=matmul(matrix,allcoor(:,:,num))
!$OMP parallel default(shared) 
!$OMP do private (natom,k,i,j,wf) schedule(dynamic,1)
           do natom=1,numatom(isurf)
             call get_wave(natom,isurf,numatom(isurf),allcoor(:,1:numatom(isurf),num),index_ele(1:numatom(isurf),isurf),wf)
             do i=1,atomwave
               k=index_orbit(i)
               do j=1,index_para(k)
                 if(dabs(wf(j,i)) .gt. tmp(k,natom)) tmp(k,natom)=dabs(wf(j,i))
               end do
             end do
           end do
!$OMP end do
!$OMP end parallel
         enddo
         101 continue


       maxwf=0.d0
       do isurf=1,nsurf
         do natom=1,numatom(isurf)
           ntype=index_ele(natom,isurf)
           do k=1,norbit
             if(maxwf(k,ntype) .lt. tmp(k,natom)) maxwf(k,ntype)=tmp(k,natom)
           end do
         end do
       end do
       do ntype=1,maxnumtype
         do k=1,norbit
           if(maxwf(k,ntype).lt.1d-20) maxwf(k,ntype)=1.d0
         end do
       end do
       numforce=0
       if(start_force==1) then
         numforce=neff*atomdim
       end if
       perindex(1:nsurf)=dsqrt(force_perindex/(neff*3d0))
     return
end subroutine 

subroutine random_divide
      implicit none
      integer(kind=intype) :: i,j,k,m,p,n,isurf,num

      index_surf_t=0
      m=0
      dimbatch=0

      isurf=1; i=1
      num=1

      index_surf_t(1)=1
      yt(1:outputneuron,1:nkpoint,1)=y(1:outputneuron,1:nkpoint,1)

      dimbatch(index_batch(1))=dimbatch(index_batch(1))+1
      allcoort(:,:,1)=allcoor(:,:,1)

        !do i=1+addpoint_t(isurf),ntpoint(isurf)+addpoint_t(isurf)
        !  index_num(i-1)=m
        !  do j=1,nkpoint
        !    do k=1,outputneuron 
        !      m=m+1
        !      yft(m)=yt(k,j,i)*ewt(i)*perindex(0)
        !      if(start_force==1) then
        !        do p=1,neff(isurf)
        !          do n=1,atomdim
        !            m=m+1
        !            yft(m)=abforcet(n,p,k,j,i)*perindex(isurf)
        !          end do
        !        end do 
        !      end if
        !    end do
        !  end do
        !end do

      index_num(nt)=ntraw
      dimjac=maxval(dimbatch)
      sumbatch=0
      do i=1,nbatch
         sumbatch(i)=sumbatch(i-1)+dimbatch(i)
      end do

      return
end subroutine

subroutine cart_to_frac(matrix,cart,fcoor)
    implicit none
    real(kind=retype) :: matrix(3,3),fcoor(3),cart(3)
    fcoor=matmul(matrix,cart)
    return
end subroutine

subroutine inverse_matrix(matrix,inv_matrix)
    implicit none
    real(kind=retype) :: matrix(3,3),inv_matrix(3,3)
    real(kind=retype) :: tmp(3,3)
    real(kind=retype),allocatable :: work(:)
    integer :: lwork,ipiv(3),info
    lwork=-1
222 tmp=matrix
    ipiv=0
    allocate(work(max(1,lwork)))
    call dgetrf(3,3,tmp,3,ipiv,info)
    call dgetri(3,tmp,3,ipiv,work,lwork,info)
    if(lwork==-1) then
      lwork=work(1)
      deallocate(work)
      goto 222
    end if
    deallocate(work)
    inv_matrix=tmp
    return
end subroutine

subroutine get_index
     implicit none
     integer(kind=intype) :: i,j,m,ntype,k,num,l,n
       m=0
       num=0
       do i=0,ipsin
         do j=0,maxnwave
           num=num+1
           index_orbit((num-1)*maxnumtype+1:num*maxnumtype)=num
           do ntype=1,maxnumtype
             index_type((num-1)*maxnumtype+ntype)=ntype
           end do
           l=0
           do k=0,nipsin(i)
             do m=0,nipsin(i)-k
               l=l+1
               index_power(1:3,l,num)=[m,k,i]
             end do
           end do
           index_para(num)=npara(i)
         end do
       end do
     return
end subroutine

subroutine get_z(ntype,wf0,z0)
     implicit none
     real(kind=retype) :: tmp1(maxnpara,norbit),z0(norbit),wf0(maxnpara,atomwave)
     integer(kind=intype) :: ntype,i,j,k,h,p,n
       tmp1=0.d0
       z0(1:norbit)=0.d0
       do i=1,atomwave
         k=index_orbit(i)
         do j=1,index_para(k)
           tmp1(j,k)=tmp1(j,k)+weight_wave(i,ntype)*wf0(j,i)
         end do
       end do
       do k=1,norbit
         do j=1,index_para(k)
           h=index_power(1,j,k)
           p=index_power(2,j,k)
           n=index_power(3,j,k)
           z0(k)=z0(k)+tmp1(j,k)*tmp1(j,k)*factor_wave(h,p,n)
         end do
       end do
     return
end subroutine

subroutine get_force
    implicit none
    integer(kind=intype) :: ksample,isurf,natom,ntype,flayer,neuron,i,j,fneuron,num,ikpoint,h,p,n,k,num1
    real(kind=retype) :: doutdg(mnl,mnl),c(mnl,mnl),tmp(maxnumforce,norbit),tmp1(maxnpara,norbit),tmp2

      forcet=0.d0
      totet=0.d0
!$omp parallel default(shared)
!$omp do private(ksample,isurf,natom,ntype,flayer,doutdg,neuron,k,i,j,h,p,n,fneuron,num,num1,c,ikpoint,tmp,tmp1,tmp2,wf,dwfdcoor)   &
!$omp firstprivate(z) schedule(dynamic,1)
      do ksample=1,nt 
        isurf=index_surf_t(ksample)
        do ikpoint=1,nkpoint
          do natom=1,numatom(isurf)
            ntype=index_ele(natom,isurf)
            call get_wave_force(natom,isurf,numatom(isurf),allcoort(:,1:numatom(isurf),ksample),  &
                  index_ele(1:numatom(isurf),isurf),wf,dwfdcoor(:,1:numatom(isurf),:,:))
            call get_z(ntype,wf,z(1:norbit,0))
            call NN(mnl,nhid(ntype),nl(0:nhid(ntype)+1,ntype),w(0:mnl,1:mnl,1:nhid(ntype)+1,ntype,ikpoint),z(0:mnl,0:nhid(ntype)+1))
            totet(1:outputneuron,ikpoint,ksample)=totet(1:outputneuron,ikpoint,ksample)+z(1:outputneuron,nhid(ntype)+1)
            doutdg=0d0
            do neuron=1,nl(nhid(ntype)+1,ntype)
              doutdg(neuron,neuron)=1d0
            end do
            do flayer=nhid(ntype)+1,1,-1
              call dgemm('N','N',nl(flayer-1,ntype),nl(nhid(ntype)+1,ntype),nl(flayer,ntype),1d0,     &
                  w(1:nl(flayer-1,ntype),1:nl(flayer,ntype),flayer,ntype,ikpoint),nl(flayer-1,ntype), &
                  doutdg(1:nl(flayer,ntype),1:nl(nhid(ntype)+1,ntype)),nl(flayer,ntype),0d0,          &
                  c(1:nl(flayer-1,ntype),1:nl(nhid(ntype)+1,ntype)),nl(flayer-1,ntype))
              do fneuron=1,nl(nhid(ntype)+1,ntype)
                do neuron=1,nl(flayer-1,ntype)
                  if(flayer>1) then
                    doutdg(neuron,fneuron)=c(neuron,fneuron)*(1d0-z(neuron,flayer-1)*z(neuron,flayer-1))
                  else 
                    doutdg(neuron,fneuron)=c(neuron,fneuron)
                  end if
                end do
              end do
            end do
            tmp1=0d0
            do i=1,atomwave
              k=index_orbit(i)
              do j=1,index_para(k)
                tmp1(j,k)=tmp1(j,k)+weight_wave(i,ntype)*wf(j,i)
              end do
            end do
            tmp=0d0
            do i=1,atomwave
              k=index_orbit(i)
              do j=1,index_para(k)
                h=index_power(1,j,k)
                p=index_power(2,j,k)
                n=index_power(3,j,k)
                tmp2=2d0*tmp1(j,k)*factor_wave(h,p,n)*weight_wave(i,ntype)
                num=0 
                do fneuron=1,neff(isurf)
                  do neuron=1,atomdim
                    num=num+1
                    tmp(num,k)=tmp(num,k)+tmp2*dwfdcoor(neuron,fneuron,j,i)
                  end do
                end do
              end do
            end do
            call dgemm('N','N',numforce(isurf),outputneuron,norbit,-1d0,tmp(1:numforce(isurf),1:norbit),&
                  numforce(isurf),doutdg(1:norbit,1:outputneuron),norbit,1d0,&
                  forcet(1:numforce(isurf),1:outputneuron,ikpoint,ksample),numforce(isurf))
          end do
        end do
      end do
!$omp end do
!$omp end parallel

    return
end subroutine

subroutine get_wave(i,isurf,natom,coor,eindex,wf0)
     implicit none
     integer(kind=intype) :: i,isurf,natom,eindex(natom)
     real(kind=retype) :: coor(atomdim,natom)
     integer(kind=intype) :: l,num,j,k,m,n,p,ntype,countnum,num1,h,ntype1,num2
     integer(kind=intype) :: index_num_dis(length(isurf)),index_dis(natom,length(isurf))
     real(kind=retype) :: vec(atomdim,natom,length(isurf)),power_vec(0:ipsin,0:ipsin,0:ipsin,natom,length(isurf))
     real(kind=retype) :: tmp1,tmp2
     real(kind=retype) :: wf0(maxnpara,atomwave),gauss_wave(natom,length(isurf),0:maxnwave)
       call gauss(i,isurf,natom,coor,eindex,index_num_dis,index_dis,vec,gauss_wave)
       call wave_factor(isurf,natom,index_num_dis,index_dis,vec,power_vec)
       wf0=0d0
       countnum=0
       ntype1=eindex(i) 
       num2=0
      !====================================simulate the s orbital which is the radial function==========================
       do m=0,ipsin
         do k=0,maxnwave
           num2=num2+1
           do num=1,length(isurf)
             do l=1,index_num_dis(num)
               j=index_dis(l,num)
               ntype=eindex(j)
               if(k<=nwave(ntype)) then
                 tmp2=gauss_wave(j,num,k)
                 num1=0
                 do n=0,nipsin(m)
                   do p=0,nipsin(m)-n
                     num1=num1+1
                     h=nipsin(m)-n-p
                     wf0(num1,countnum+ntype)=wf0(num1,countnum+ntype)+tmp2*power_vec(h,p,n,j,num)
                   end do
                 end do
               end if
             end do
           end do
           wf0(:,countnum+1:countnum+maxnumtype)=wf0(:,countnum+1:countnum+maxnumtype)/maxwf(num2,ntype1)
           countnum=countnum+maxnumtype
         end do
       end do
     return
end subroutine

subroutine get_wave_force(i,isurf,natom,coor,eindex,wf0,dwfdcoor0)
     implicit none
     integer(kind=intype) :: i,isurf,natom,eindex(natom)
     real(kind=retype) :: coor(atomdim,natom)
     integer(kind=intype) :: l,num,j,k,m,n,p,ntype,countnum,num1,h,ntype1,num2
     integer(kind=intype) :: index_num_dis(length(isurf)),index_dis(natom,length(isurf))
     real(kind=retype) :: vec(atomdim,natom,length(isurf))
     real(kind=retype) :: tmp(atomdim),tmp1,tmp2,tmp3,tmp4(atomdim),tmp5
     real(kind=retype) :: dwfdcoor0(atomdim,natom,maxnpara,atomwave)
     real(kind=retype) :: gauss_wave(natom,length(isurf),0:maxnwave),wf0(maxnpara,atomwave)
     real(kind=retype) :: der_gauss_wave(atomdim,natom,length(isurf),0:maxnwave)
     real(kind=retype) :: power_vec(-1:ipsin,-1:ipsin,-1:ipsin,natom,length(isurf))
     real(kind=retype) :: der_power_vec(atomdim,0:ipsin,0:ipsin,0:ipsin,natom,length(isurf))
     call der_gauss(i,isurf,natom,coor,eindex,index_num_dis,index_dis,vec,gauss_wave,der_gauss_wave)
     call der_wave_factor(isurf,natom,index_num_dis,index_dis,vec,power_vec,der_power_vec)
     dwfdcoor0=0.d0
     countnum=0
     ntype1=eindex(i)
     num2=0
     wf0=0.d0
     do m=0,ipsin
       do k=0,maxnwave
         num2=num2+1
         tmp5=maxwf(num2,ntype1)
         do num=1,length(isurf)
           do l=1,index_num_dis(num)
             j=index_dis(l,num)
             ntype=eindex(j)
             if(k<=nwave(ntype)) then 
               tmp=der_gauss_wave(:,j,num,k)
               tmp1=gauss_wave(j,num,k)
               num1=0
               do n=0,nipsin(m)
                 do p=0,nipsin(m)-n
                   num1=num1+1
                   h=nipsin(m)-n-p
                   wf0(num1,countnum+ntype)=wf0(num1,countnum+ntype)+tmp1*power_vec(h,p,n,j,num)
                   tmp4=(tmp*power_vec(h,p,n,j,num)+tmp1*der_power_vec(:,h,p,n,j,num))/tmp5
                   dwfdcoor0(:,j,num1,countnum+ntype)=dwfdcoor0(:,j,num1,countnum+ntype)+tmp4
                   dwfdcoor0(:,i,num1,countnum+ntype)=dwfdcoor0(:,i,num1,countnum+ntype)-tmp4
                 end do
               end do 
             end if
           end do
         end do
         wf0(:,countnum+1:countnum+maxnumtype)=wf0(:,countnum+1:countnum+maxnumtype)/tmp5
         countnum=countnum+maxnumtype
       end do
     end do
     return
end subroutine

subroutine NN(mnl,nhid,nl,w,z)
     implicit none
     integer(kind=intype) :: mnl,nhid,nl(0:nhid+1)
     real(kind=retype) :: w(0:mnl,mnl,nhid+1),z(0:mnl,0:nhid+1)
     real(kind=retype),external :: ddot
     integer(kind=intype) :: ilayer,neu,num
     do ilayer=1,nhid+1
        do neu=1,nl(ilayer) 
           num=nl(ilayer-1)+1
           z(neu,ilayer)=ddot(num,z(0:nl(ilayer-1),ilayer-1),1,w(0:nl(ilayer-1),neu,ilayer),1)
           if(ilayer<nhid+1) z(neu,ilayer)=dtanh(z(neu,ilayer))
        end do
     end do
     return
end subroutine

subroutine feedforward
     implicit none
     real(kind=retype) :: tmp
     integer(kind=intype) :: natom,ntype,ksample,isurf,ikpoint

     totet=0.d0
!$omp parallel default(shared)
!$omp do private(ksample,natom,ntype,isurf,ikpoint,wf) firstprivate(z)
     do ksample=1,nt
        isurf=index_surf_t(ksample)
        do ikpoint=1,nkpoint
           do natom=1,numatom(isurf) 
              ntype=index_ele(natom,isurf)
              call get_wave(natom,isurf,numatom(isurf),allcoort(:,1:numatom(isurf),ksample),index_ele(1:numatom(isurf),isurf),wf)
              call get_z(ntype,wf,z(1:norbit,0))
              call NN(mnl,nhid(ntype),nl(0:nhid(ntype)+1,ntype),w(0:mnl,1:mnl,1:nhid(ntype)+1,ntype,ikpoint),z(0:mnl,0:nhid(ntype)+1))
              totet(1:outputneuron,ikpoint,ksample)=totet(1:outputneuron,ikpoint,ksample)+z(1:outputneuron,nhid(ntype)+1)
           end do
        end do
     end do
!$omp end do
!$omp end parallel

     return
end subroutine

subroutine nt_nv
     implicit none
     integer(kind=intype) :: isurf,num,i,j
     logical :: m11(numpoint)
       nt=0
       ntraw=0
       addpoint_t=0
       addfpoint_t=0
       numth=0

       m11=.true.

       num=0
       ntpoint=0

       do isurf=1,nsurf
         do i=1,npoint(isurf)
           m00(i,isurf)=m11(i+num)
           if(m11(i+num).eq..true.) ntpoint(isurf)=ntpoint(isurf)+1
         end do
         num=num+npoint(isurf)
       end do

       do isurf=1,nsurf
         addpoint_t(isurf)=nt
         addfpoint_t(isurf)=ntraw
         nt=nt+ntpoint(isurf) ! number of training data
         numth(isurf)=numth(isurf-1)+ntpoint(isurf)*(numforce(isurf)+1)
         if(start_force==1) then
           ntraw=ntraw+ntpoint(isurf)*((1+numforce(isurf))*outputneuron*nkpoint)
         else 
           ntraw=ntraw+ntpoint(isurf)
         end if
       end do
     return
end subroutine

subroutine period(k,isurf,natom,coor,vec)
     implicit none
     integer(kind=intype) :: k,isurf,natom
     real(kind=retype) :: coor(atomdim,natom),vec(atomdim,natom,length(isurf))
     real(kind=retype) :: ycoor(atomdim,natom),fcoor(atomdim,natom),vector(atomdim)
     integer(kind=intype) :: i,j,m,n,num,ninit
     real(kind=retype) :: bond,sca
       ycoor=coor
       ninit=(length(isurf)+1)/2
       do i=1,natom
         if(i/=k) then
           do j=1,atomdim
             bond=ycoor(j,i)-ycoor(j,k)
             sca=anint(bond)*1d0
             fcoor(j,i)=ycoor(j,i)-sca
           end do
         else
           fcoor(:,i)=ycoor(:,i)
         end if
       end do
       vector=0d0
       num=0
       do n=-nimage(3,isurf),nimage(3,isurf)
         vector(3)=1d0*n
         do j=-nimage(2,isurf),nimage(2,isurf)
           vector(2)=1d0*j
           do m=-nimage(1,isurf),nimage(1,isurf)
             vector(1)=1d0*m
             num=num+1
             do i=1,natom
               if(num/=ninit.or.i/=k) then
                 vec(:,i,num)=fcoor(:,i)-fcoor(:,k)+vector
                 vec(:,i,num)=matmul(scal(:,:),vec(:,i,num)) ! chenjun
               end if
             end do
           end do
         end do
       end do
     return
end subroutine

subroutine distance(i,isurf,ntype,natom,vec,index_num_dis,index_dis,dis)
     implicit none
     integer(kind=intype) :: i,isurf,ntype,natom,index_num_dis(length(isurf)),index_dis(natom,length(isurf))
     real(kind=retype) :: vec(atomdim,natom,length(isurf))
     real(kind=retype) :: dis(natom,length(isurf))
     real(kind=retype),external :: dnrm2
     integer(kind=intype) :: num,j,ninit,k
       index_num_dis=0
       ninit=(length(isurf)+1)/2
       index_dis=0
       do num=1,length(isurf)
         do j=1,natom
           do k=1,atomdim
             if(rc(ntype)<vec(k,j,num)) goto 124
           end do
           if(num/=ninit.or.i/=j) then
             dis(j,num)=dnrm2(atomdim,vec(:,j,num),1)
             if(dis(j,num)<=rc(ntype)) then
               index_num_dis(num)=index_num_dis(num)+1
               index_dis(index_num_dis(num),num)=j
             end if
           end if
124        continue
         end do
       end do
     return
end subroutine

subroutine effectcos(isurf,ntype,natom,index_num_dis,index_dis,dis,effect)
     implicit none
     integer(kind=intype) :: isurf,ntype,natom,index_num_dis(length(isurf)),index_dis(natom,length(isurf))
     real(kind=retype) :: dis(natom,length(isurf)),tmp
     real(kind=retype) :: effect(natom,length(isurf))
     integer(kind=intype) :: j,num,l!i is the centre atom
       do num=1,length(isurf)
         do l=1,index_num_dis(num)
           j=index_dis(l,num)
           tmp=0.5d0*dcos(pi*dis(j,num)/rc(ntype))+0.5d0
           effect(j,num)=tmp*tmp
         end do
       end do
     return
end subroutine

subroutine deffect(isurf,ntype,natom,dis,index_num_dis,index_dis,effect,dedr)
     implicit none
     integer(kind=intype) :: isurf,ntype,natom,index_num_dis(length(isurf)),index_dis(natom,length(isurf))
     real(kind=retype) :: dis(natom,length(isurf)),tmp
     real(kind=retype) :: dedr(natom,length(isurf)),effect(natom,length(isurf))
     integer(kind=intype) :: j,l,num !k is the centre atom
       do num=1,length(isurf)
         do l=1,index_num_dis(num)
           j=index_dis(l,num)
           tmp=0.5d0*dcos(pi*dis(j,num)/rc(ntype))+0.5d0
           effect(j,num)=tmp*tmp
           dedr(j,num)=-tmp*dsin(pi*dis(j,num)/rc(ntype))*pi/rc(ntype)
         end do
       end do
     return
end subroutine

subroutine dbondlendx(isurf,natom,vec,dis,index_num_dis,index_dis,drdx)
     implicit none
     integer(kind=intype) :: isurf,natom,index_num_dis(length(isurf)),index_dis(natom,length(isurf))
     real(kind=retype) :: dis(natom,length(isurf))
     real(kind=retype) :: drdx(atomdim,natom,length(isurf)),vec(atomdim,natom,length(isurf))
     integer(kind=intype) :: j,l,num! dr/dxi
       do num=1,length(isurf)
         do l=1,index_num_dis(num)
           j=index_dis(l,num)
           drdx(:,j,num)=vec(:,j,num)/dis(j,num)
         end do
       end do
     return
end subroutine

subroutine gauss(i,isurf,natom,coor,eindex,index_num_dis,index_dis,vec,gauss_wave)
     implicit none
     integer(kind=intype) :: i,isurf,natom,index_num_dis(length(isurf)),index_dis(natom,length(isurf)),eindex(natom)
     real(kind=retype) :: coor(atomdim,natom)
     real(kind=retype) :: dis(natom,length(isurf)),effect(natom,length(isurf))
     real(kind=retype) :: gauss_wave(natom,length(isurf),0:maxnwave),vec(atomdim,natom,length(isurf))
     real(kind=retype) :: tmp
     integer(kind=intype) :: j,l,num,m,k,ntype
       ntype=eindex(i)
       call period(i,isurf,natom,coor,vec)
       call distance(i,isurf,ntype,natom,vec,index_num_dis,index_dis,dis)
       call effectcos(isurf,ntype,natom,index_num_dis,index_dis,dis,effect)
       do k=0,maxnwave
         do num=1,length(isurf)
           do l=1,index_num_dis(num)
             j=index_dis(l,num)
             ntype=eindex(j)
             if(k<=nwave(ntype)) then
               tmp=dis(j,num)-rs(ntype,k)
               gauss_wave(j,num,k)=dexp(-inta(ntype,k)*tmp*tmp)*effect(j,num)
             end if
           end do 
         end do
       end do
     return
end subroutine

subroutine der_gauss(i,isurf,natom,coor,eindex,index_num_dis,index_dis,vec,gauss_wave,der_gauss_wave)
     implicit none
     integer(kind=intype) :: i,isurf,natom,index_num_dis(length(isurf)),index_dis(natom,length(isurf)),eindex(natom)
     real(kind=retype) :: coor(atomdim,natom)
     real(kind=retype) :: dis(natom,length(isurf)),effect(natom,length(isurf))
     real(kind=retype) :: dedr(natom,length(isurf)),drdx(atomdim,natom,length(isurf))
     real(kind=retype) :: gauss_wave(natom,length(isurf),0:maxnwave),vec(atomdim,natom,length(isurf))
     real(kind=retype) :: der_gauss_wave(atomdim,natom,length(isurf),0:maxnwave)
     real(kind=retype) :: tmp,tmp1,tmp2
     integer(kind=intype) :: j,l,num,m,k,ntype
       ntype=eindex(i)
       call period(i,isurf,natom,coor,vec)
       call distance(i,isurf,ntype,natom,vec,index_num_dis,index_dis,dis)
       call deffect(isurf,ntype,natom,dis,index_num_dis,index_dis,effect,dedr)
       call dbondlendx(isurf,natom,vec,dis,index_num_dis,index_dis,drdx)
       do k=0,maxnwave
         do num=1,length(isurf)
           do l=1,index_num_dis(num)
             j=index_dis(l,num)
             ntype=eindex(j)
             if(k<=nwave(ntype)) then
               tmp=dis(j,num)-rs(ntype,k)
               tmp1=-inta(ntype,k)*tmp
               tmp2=dexp(tmp1*tmp)
               gauss_wave(j,num,k)=tmp2*effect(j,num)
               der_gauss_wave(:,j,num,k)=(2d0*tmp1*effect(j,num)+dedr(j,num))*tmp2*drdx(:,j,num)
             end if
           end do 
         end do
       end do
     return
end subroutine

subroutine wave_factor(isurf,natom,index_num_dis,index_dis,vec,power_vec)
     implicit none
     integer(kind=intype) :: isurf,natom,index_num_dis(length(isurf)),index_dis(natom,length(isurf))
     real(kind=retype) :: vec(atomdim,natom,length(isurf)),tmp_vec(atomdim,0:ipsin,natom,length(isurf))
     real(kind=retype) :: power_vec(0:ipsin,0:ipsin,0:ipsin,natom,length(isurf))
     integer(kind=intype) :: i,j,k,l,m,n,p,num
       tmp_vec=1d0
       do num=1,length(isurf)
         do l=1,index_num_dis(num)
           j=index_dis(l,num)
           do m=1,ipsin
             do k=1,atomdim
               tmp_vec(k,m,j,num)=tmp_vec(k,m-1,j,num)*vec(k,j,num)
             end do
           end do
         end do
       end do  
       do i=0,ipsin 
         do num=1,length(isurf)
           do l=1,index_num_dis(num)
             j=index_dis(l,num)
             do m=0,nipsin(i)
               do n=0,nipsin(i)-m
                 p=nipsin(i)-m-n
                 power_vec(p,n,m,j,num)=tmp_vec(1,m,j,num)*tmp_vec(2,n,j,num)*tmp_vec(3,p,j,num)
               end do
             end do
           end do
         end do
       end do
     return
end subroutine

subroutine der_wave_factor(isurf,natom,index_num_dis,index_dis,vec,power_vec,der_power_vec)
     implicit none
     integer(kind=intype) :: isurf,natom,index_num_dis(length(isurf)),index_dis(natom,length(isurf))
     real(kind=retype) :: vec(atomdim,natom,length(isurf)),tmp_vec(atomdim,0:ipsin,natom,length(isurf))
     real(kind=retype) :: power_vec(-1:ipsin,-1:ipsin,-1:ipsin,natom,length(isurf))
     real(kind=retype) :: der_power_vec(atomdim,0:ipsin,0:ipsin,0:ipsin,natom,length(isurf))
     integer(kind=intype) :: i,j,k,l,m,n,p,num
       tmp_vec=1d0
       do num=1,length(isurf)
         do l=1,index_num_dis(num)
           j=index_dis(l,num)
           do m=1,ipsin
             do k=1,atomdim
               tmp_vec(k,m,j,num)=tmp_vec(k,m-1,j,num)*vec(k,j,num)
             end do
           end do
         end do
       end do  
       power_vec=1d0
       do i=0,ipsin 
         do num=1,length(isurf)
           do l=1,index_num_dis(num)
             j=index_dis(l,num)
             do m=0,nipsin(i)
               do n=0,nipsin(i)-m
                 p=nipsin(i)-m-n
                 power_vec(p,n,m,j,num)=tmp_vec(1,m,j,num)*tmp_vec(2,n,j,num)*tmp_vec(3,p,j,num)
                 der_power_vec(1,p,n,m,j,num)=dble(m)*power_vec(p,n,m-1,j,num)
                 der_power_vec(2,p,n,m,j,num)=dble(n)*power_vec(p,n-1,m,j,num)
                 der_power_vec(3,p,n,m,j,num)=dble(p)*power_vec(p-1,n,m,j,num)
               end do
             end do
           end do
         end do
       end do
     return
end subroutine

