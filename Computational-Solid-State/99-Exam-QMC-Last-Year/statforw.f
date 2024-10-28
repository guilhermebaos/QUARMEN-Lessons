c read a stream of data a(i) from stdin; 
c perform statistical analysis, see statfor.pdf
      implicit none
      integer i,m
      parameter(m=1000000)
      real*8 a(m),w(m)
      do i=1,m
       read(*,*,end=1)a(i),w(i)
      enddo
    1 i=i-1
      call corr(a,w,i)
      call blocking(a,i)
      call histo(a,i)
      stop
      end

      subroutine corr(a,w,n)
c estimated mean error and autocorrelation
      implicit none
      integer mc,i,k,n,l
      parameter(mc=200)
      real*8 a(n),w(n),ave,peso,neff,v,rkappa,c(mc),f
c mean
      ave=0.d0
      peso=0.d0
      do i=1,n
       ave=ave+a(i)*w(i)
       peso=peso+w(i)
      enddo
      ave=ave/peso
      write(*,'(''average   '',f40.10)')ave
c neff
      neff=0.d0
      do i=1,n
       neff=neff+w(i)**2
      enddo
      neff=peso**2/neff
      write(*,'(''n,neff    '',2f20.6)')float(n),neff
c variance
      v=0.d0
      do i=1,n
       v=v+w(i)*(a(i)-ave)**2
      enddo
      v=v/peso
      v=v*neff/(neff-1.d0)
      write(*,'(''variance '',f40.10)')v
c autocorrelation (up to a mc steps) and correlation time
      open(2,file='corr.out')
      l=min(n-1,mc)
      do i=1,l
       c(i)=0.d0
       peso=0.d0
       do k=1,n-i
        c(i)=c(i)+sqrt(w(k)*w(k+i))*(a(k)-ave)*(a(k+i)-ave)
        peso=peso+sqrt(w(k)*w(k+i))
       enddo
       c(i)=c(i)/peso/v
      enddo
      rkappa=1.d0
      f=1
      do i=1,l
       write(2,*)i,c(i)
       if(c(i).lt.0.d0)f=0
       rkappa=rkappa+2*c(i)*f
      enddo  
      close(2)
      rkappa=max(1.d0,rkappa)
      write(*,'(''t corr   '',f40.10)')rkappa
c effective number of data
      write(*,'(''n eff    '',f40.10)')neff/rkappa
c error of mean
      write(*,'(''sigma    '',f40.10)')sqrt(v*rkappa/neff)
      return
      end

      subroutine blocking(a,n)
c blocking analysis
      implicit none
      integer i,k,l,n,nblk,large,isize,isize_step,minleft,nsizes
      parameter(minleft=20,nsizes=100)
      real*8 a(n),ab,average,average2,error
      open(2,file='blocking.out')
      large=n/minleft ! want at least minleft blocks left
      isize_step=max(1,large/nsizes) ! want at most ~nsizes block sizes
c loop on block size
      do isize=1,large,isize_step
c # blocks
       nblk=n/isize
       k=0
       average=0.d0
       average2=0.d0
       do i=1,nblk
c block averages
        ab=0.d0
        do l=1,isize
         k=k+1
         ab=ab+a(k)
        enddo
        ab=ab/isize
        average=average+ab
        average2=average2+ab*ab
       enddo
       average=average/nblk
       average2=average2/nblk
c estimated error of mean at this block size
       error=sqrt((average2-average**2)/(nblk-1))
       write(2,*)isize,error
      enddo
      close(2)
      return
      end

      subroutine histo(a,n)
      implicit none
      integer i,j,n,m,nbin
      parameter(nbin=21)
      real*8 a(n),h(0:nbin+1),a_min,a_max,delta
      do j=0,nbin+1
       h(j)=0.d0
      enddo
c min and max
      a_min=a(1)
      a_max=a(1)
      do i=2,n
       a_min=min(a_min,a(i))
       a_max=max(a_max,a(i))
      enddo
c bin size
      delta=(a_max-a_min)/nbin
c histogram
      do i=1,n
       j=nint(0.5d0+(a(i)-a_min)/delta)
       h(j)=h(j)+1.d0
      enddo
      h(1)=h(1)+h(0)
      h(nbin)=h(nbin)+h(nbin+1)
c write
      open(2,file='histo.out')
      do j=1,nbin
       write(2,*)a_min+(j-0.5d0)*delta,h(j)/(n*delta)
      enddo
      close(2)
      return
      end
