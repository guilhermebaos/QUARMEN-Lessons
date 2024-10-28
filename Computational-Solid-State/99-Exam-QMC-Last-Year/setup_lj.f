c setup program for 2d bosons with lj potential
      implicit none
      integer mnk,mdim,msites,mbasis,mnts
      parameter(mnk=9000,mdim=2,msites=100,mbasis=2,mnts=501  )
      integer nt,np,len,i,j,ndim,ng,ngfound,ifv0
      integer ifex,iexp,iunit,i_cub(mdim),nbasis
      integer nkappa
      integer jp,m_parm
      parameter(m_parm=20)
      real*8 p(m_parm)
      common /c_wrkparm/p,jp
      common /c_v/vshift,vtail
      common /c_lj/epsilon,sigma
      real*8 drt,rs,hbs2m,gvect(mdim,mnk),a(mdim,mdim),cut
      real*8 rho,el(2),l2,tail,d,pi,v0,sites(mdim,msites)
      real*8 aux(mdim),delta,basis(mdim,mbasis)
     &      ,factor,abrav(mdim),shift(mdim),gg,gg0,gstart,gdelta
      real*8 vshift,vtail,epsilon,sigma,alat
      integer phase_flag
      character*30 runid,filename,routinename,te_name
      character word*200
      logical l_ex
      common /c_rs/rs
      external lj2d,mcmillan,csi,gauss
      pi=acos(-1.d0)
      iunit=10
c
      write(*,*)'Setup program for VMC and DMC simulation of 2d bosons'
      write(*,*)'with mass 4 a.m.u. interacting with a LJ potential.'
      write(*,*)'Units: Kelvin for energies, Angstrom for lengths.'
      write(*,*)' '
      write(*,*)'enter run id'
      read(*,'(a)')runid
c     write(*,*)'enter LJ parameters epsilon, sigma'
c     read(*,*)epsilon,sigma
      epsilon=10.22d0
      sigma=2.556d0
      write(*,*)'enter density'
      read(*,*)rho
c     write(*,*)'number of particles (choose either 16 or 30)'
c     read(*,*)np
      np=16
      write(*,*)'choose wave function: Jastrow or Nosanow (0/1)'
      read(*,*)phase_flag

c runid file
      write(*,*)'creating file "runid" with your chosen run id...'
      open(2,file='runid')
      write(2,*)trim(runid)
      close(2)

c .sy file
      write(*,*)'creating file .sy with system description...'
      ndim=2
      len=index(runid,' ')-1
      open(iunit,file=runid(1:len)//'.sy',status='unknown')
      write(word,*)'ndim 2'
      call writestring(word,iunit)
      hbs2m=6.0596
      alat=sqrt(2.d0/(sqrt(3.d0)*rho))
      if(np.eq.16)then
       i_cub(1)=4
       i_cub(2)=2
      elseif(np.eq.30)then
       i_cub(1)=5
       i_cub(2)=3
      else
       stop'stop: choose either 16 or 30 particles'
      endif
      el(1)=i_cub(1)*alat
      el(2)=i_cub(2)*alat*sqrt(3.d0)

c up particle positions
      write(*,*)'creating file .lj.x with particle positions...'
      write(word,*)'type lj',  np,  hbs2m,' ',runid(1:len)//'.lj.x'
      call writestring(word,iunit)
      ifex=1
      filename=runid(1:len)//'.lj.x'
      nbasis=2
      basis(1,1)=0.d0
      basis(2,1)=0.d0
      basis(1,2)=0.5d0*alat
      basis(2,2)=0.5d0*alat*sqrt(3.d0)
      abrav(1)=alat
      abrav(2)=alat*sqrt(3.d0)
      delta=alat*0.2
      call xtal_sites(mdim,msites,abrav,ndim,basis,nbasis,rho
     &               ,np,i_cub,shift,delta,aux,sites,factor,1)
      open(2,file=filename,status='unknown')
      do i=1,np
       write(2,*)(sites(j,i),j=1,ndim)
      enddo
      close(2)

c lattice sites
      if(phase_flag.eq.1)then
       write(*,*)'creating file .s with lattice sites...'
       write(word,*)'sites s',np,' ',runid(1:len)//'.s'
       call writestring(word,iunit)
       filename=runid(1:len)//'.s'
       nbasis=2
       basis(1,1)=0.d0
       basis(2,1)=0.d0
       basis(1,2)=0.5d0*alat
       basis(2,2)=0.5d0*alat*sqrt(3.d0)
       abrav(1)=alat
       abrav(2)=alat*sqrt(3.d0)
       delta=0.d0
       call xtal_sites(mdim,msites,abrav,ndim,basis,nbasis,rho
     &                ,np,i_cub,shift,delta,aux,sites,factor,1)
       open(2,file=filename,status='unknown')
       do i=1,np
        write(2,*)(sites(j,i),j=1,ndim)
       enddo
       close(2)
      endif

c numero di punti e passo delle tabelle
      nt=1000
      drt=0.5d0*min(el(1),el(2))/nt

c lj potential
      write(*,*)'creating file .v with tabulated potential'
      filename=runid(1:len)//'.v'
      write(routinename,'(2a,28x)')'lj2d'
      call tgen(lj2d,filename,nt,drt,1,routinename)
      write(word,*)'v2 lj lj ',filename,' 0'
      call writestring(word,iunit)
      write(word,*)'vshift ',vshift
      call writestring(word,iunit)
      write(word,*)'vtail ',vtail
      call writestring(word,iunit)

c pair pseudopotential
      write(*,*)'creating file .u with tabulated pair correlation...'
      filename=runid(1:len)//'.u'
      write(routinename,'(8a,22x)')'mcmillan'
      jp=2
c     p(1)=sigma**5
c     p(2)=5.d0
      if(phase_flag.eq.0)then
       call mcmillan_parm_liq(rho,p)
      else
       call mcmillan_parm_sol(rho,p)
      endif
      call tgen(mcmillan,filename,nt,drt,0,routinename)
      write(word,*)'u2 lj lj ',filename
      call writestring(word,iunit)

      write(word,*)'pbc',el(1),el(2)
      call writestring(word,iunit)

      if(phase_flag.eq.0)then
c 3body
       write(*,*)'creating file .u3 with tabulated 3body correlation...'
       filename=runid(1:len)//'.u3'
       write(routinename,'(a3,27x)')'csi'
       jp=3
c      p(1)=1.d0
c      p(2)=sigma*0.1d0
c      p(3)=1.d0/sigma
       call csi_parm_liq(rho,p)
       call tgen(csi,filename,nt,drt,0,routinename)
       write(word,*)'u3 lj lj ',filename
       call writestring(word,iunit)
      else
c nosanow
       write(*,*)'creating file .n with tabulated nosanow function...'
       filename=runid(1:len)//'.n'
       write(routinename,'(a5,25x)')'gauss'
       jp=1
c      p(1)=1.d0/sigma
       call gauss_parm_sol(rho,p)
       call tgen(gauss,filename,nt,drt,0,routinename)
       write(word,*)'nosanow lj s ',filename
       call writestring(word,iunit)
      endif
c sofk
      write(*,*)'turn on calculation of S(k)...'
      write(word,*)'sofk'
      call writestring(word,iunit)
c gofr
      write(*,*)'turn on calculation of g(r)...'
      write(word,*)'gofr'
      call writestring(word,iunit)

      stop
      end
 
      subroutine tgen(a,file,nt,drt,stdin_flag,routinename)
      implicit none
      integer mgrid,i,j,nt,m_parm,jp,stdin_flag,spline_flag
      parameter(mgrid=10001,m_parm=20)
      real*8 t(0:mgrid,4),drt,p(m_parm)
      character*20 file
      character*30 routinename
      common /c_wrkparm/p,jp
      external a
c     jp=0 ! OKKIO
c make table
      call a(0.d0,drt,0,nt,t(0,1),t(0,2),stdin_flag,spline_flag)
      call spline(mgrid,nt,drt,t,spline_flag)
c write table
      open(2,file=file,status='unknown')
      write(2,*)nt,drt
      do i=0,nt
       write(2,'(4e20.12)')(t(i,j),j=1,4)
      enddo
c this stuff is needed for optimization
      write(2,*)routinename
      write(2,*)jp
      j=1
      do i=1,jp
       write(2,*)p(i),j
      enddo
      close(2)
      return
      end

      subroutine spline(mgrid,nt,drt,t,flag)
c scale and spline
      implicit none
      integer mgrid,nt,i,flag
      real*8 drt,t(0:mgrid,4)
      do i=0,nt
       t(i,2)=t(i,2)*drt
       t(i,3)=t(i,3)*drt*drt*0.5d0
      enddo
c this cubic spline has the correct values of f and df/dr at the grid points.
      do i=0,nt-1
       t(i,3)= 3.d0*(t(i+1,1)-t(i,1))  -(t(i+1,2)+2.d0*t(i,2))
       t(i,4)=-2.d0*(t(i+1,1)-t(i,1))  +(t(i+1,2)+     t(i,2))
      enddo
c zero last point
      if(flag.ne.0)then
       t(nt,1)=0.d0
       t(nt,2)=0.d0
       t(nt,3)=0.d0
       t(nt,4)=0.d0
      endif
      return
      end

      subroutine xtal_sites(mdim,msites,abrav,ndim,basis,nbasis,rho,np
     &                     ,i_cub,shift,delta,el,sites,factor,verbose)
c abrav(ndim)      sono i vettori primitivi del reticolo di Bravais;
c                  indice = indice componente cartesiana = indice vettore
c basis(mdim,mbasis) sono i vettori di base, ibasis=2,mbasis
c i_cubetto(ndim) sono le dimensioni (in cubetti lungo idim) della cella;
      implicit none
      integer mdim,msites,ndim,nbasis,i_cub(ndim),np,n_sites
     &       ,i,j,npts,nx(3),ix(3),idim,jdim,ipts,verbose,n_sites_left
     &       ,n_p_left
      real*8 abrav(ndim),basis(mdim,nbasis),el(ndim)
     &      ,vol0,rho,shift(ndim),delta,factor
     &      ,sites(mdim,msites),aux,rannyu
      if(verbose.ne.0)write(6,*)'=========>> xtal_sites <<=========='
c no more than 3d
      if(ndim.gt.3)stop'ndim.gt.3: stop...'
c consistency between i_cub and np
      n_sites=i_cub(1)
      do i=2,ndim
       n_sites=n_sites*i_cub(i)
      enddo
      n_sites=n_sites*nbasis
      if(verbose.ne.0)write(6,*)'n_sites = ',n_sites,'  np = ',np
      if(n_sites.lt.np)stop'n_sites.lt.np: stop...'
c volume della cella primitiva
      vol0=abrav(1)
      do i=2,ndim
       vol0=vol0*abrav(i)
      enddo
c aggiusta la lunghezza degli abrav in modo che vol0 --> nbasis/rho
      factor=(nbasis/(rho*vol0))**(1.d0/ndim)
      do i=1,ndim
       abrav(i)=abrav(i)*factor
      enddo
      vol0=nbasis/rho
      do i=1,nbasis
       do idim=1,ndim
        basis(idim,i)=basis(idim,i)*factor
       enddo
      enddo
      do i=1,ndim
       shift(i)=shift(i)*factor
      enddo
c lati della cella di simulazione
      do i=1,ndim
       el(i)=abrav(i)*i_cub(i)
      enddo
c numero di siti del reticolo di bravais
      npts=n_sites/nbasis
c quanti siti ci sono tra tutti gli idim piu' piccoli
      do idim=1,ndim
       nx(idim)=1
       do jdim=2,idim
        nx(idim)=nx(idim)*i_cub(jdim-1)
       enddo
      enddo
c loop sui siti
      n_sites_left=n_sites
      n_p_left=np
      do ipts=1,npts
c ricostruisce gli indici
       i=ipts-1
       do idim=ndim,2,-1
        ix(idim)=(i)/nx(idim)
        i=mod(i,nx(idim))
       enddo
       ix(1)=i
c loop sui vettori di base
       do j=1,nbasis
        if(float(n_p_left)/float(n_sites_left).gt.rannyu())then
         n_p_left=n_p_left-1
c particle position (+shift +random diplacement + pbc)
         do idim=1,ndim
          aux=ix(idim)*abrav(idim)+basis(idim,j)
     &       +shift(idim)+2.d0*delta*(0.5d0-rannyu())
          sites(idim,np-n_p_left)=aux-el(idim)*nint(aux/el(idim))
         enddo
        endif
        n_sites_left=n_sites_left-1
       enddo
      enddo
      if(verbose.ne.0)then
       write(6,*)'nsites np nvac ',n_sites,np,n_sites-np
       write(6,*)'-------> end xtal_sites <--------'
      endif
      return
      end

      subroutine writestring(string,iunit)
      implicit none
      integer iunit,i,l,j
      character*79 record
      character*200 string
      l=200
      record(1:1)=string(1:1)
      do i=2,79
       record(i:i)=' '
      enddo
      i=2
      do j=2,l
       if(i.gt.80)stop'max record length reached'
       if(string(j:j).ne.' ')then
        record(i:i)=string(j:j)
        i=i+1
       elseif(string(j-1:j-1).ne.' ')then
        record(i:i)=string(j:j)
        i=i+1
       endif
      enddo
      write(iunit,*)record
      return
      end

      subroutine r_set(n,a,b)
      implicit none
      integer i,n
      real*8 a(n),b
      do i=1,n
       a(i)=b
      enddo
      return
      end

      double precision function rannyu()
      real*8 twom12
      parameter(twom12=0.000244140625d0)
      common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4,nbit,irnyuc
c
c generic statement functions
c
      ishft12(ii)=ii/4096
      mask12(ii)=mod(ii,4096)
c
c unix f77 statement functions
c
c     ishft12(ii)=rshift(ii,12)
c     mask12(ii)=and(ii,4095)
c
c fps statement functions
c
c     ishft12(ii)=shift(ii,-12)
c     mask12(ii)=and(ii,4095)
c
c cray cft statement functions
c
c     ishft12(ii)=shiftr(ii,12)
c     mask12(ii)=and(ii,4095)
c
c vms fortran and convex fc statement functions
c
c     ishft12(ii)=ishft(ii,-12)
c     mask12(ii)=iand(ii,4095)
      i1=l1*m4+l2*m3+l3*m2+l4*m1
      i2=l2*m4+l3*m3+l4*m2
      i3=l3*m4+l4*m3
      i4=l4*m4
      l4=mask12(i4)
      i3=i3+ishft12(i4)
      l3=mask12(i3)
      i2=i2+ishft12(i3)
      l2=mask12(i2)
      l1=mask12(i1+ishft12(i2))
      rannyu=twom12*(l1+
     +       twom12*(l2+
     +       twom12*(l3+
     +       twom12*(l4))))
      return
      end
      subroutine setrn(iseed)
      common /rnyucm/ m(4),l(4),nbit,irnyuc
      integer iseed(4)
      ishft12(ii)=ii/4096
      mask12(ii)=mod(ii,4096)
      do 10 i=1,4
   10 l(i)=mod(iseed(i),4096)
      l(4)=2*(l(4)/2)+1
c
c shift everything to the left if not 48 bit
c
      if (nbit.lt.48) then
         do 20 i=1,48-nbit
         i1=l(1)*2
         i2=l(2)*2
         i3=l(3)*2
         i4=l(4)*2
         l(4)=mask12(i4)
         i3=i3+ishft12(i4)
         l(3)=mask12(i3)
         i2=i2+ishft12(i3)
         l(2)=mask12(i2)
         l(1)=mask12(i1+ishft12(i2))
   20    continue
         endif
      return
      end
      subroutine savern(iseed)
      common /rnyucm/ m(4),l(4),nbit,irnyuc
      integer iseed(4)
      do 10 i=1,4
      iseed(i)=l(i)
   10 continue
c
c shift everything to the right if not 48 bit
c
      if (nbit.lt.48) then
         do 20 i=1,48-nbit
         i1=iseed(1)/2
         ir=iseed(1)-i1*2
         iseed(2)=iseed(2)+4096*ir
         i2=iseed(2)/2
         ir=iseed(2)-i2*2
         iseed(3)=iseed(3)+4096*ir
         i3=iseed(3)/2
         ir=iseed(3)-i3*2
         iseed(4)=iseed(4)+4096*ir
         i4=iseed(4)/2
         iseed(1)=i1
         iseed(2)=i2
   20    continue
         endif
      return
      end

      block data randat
      common /rnyucm/ m(4),l(4),nbit,irnyuc
      data m / 502,1521,4071,2107/
      data l /   0,   0,   0,   1/
      data nbit /48/
      end

      subroutine mcmillan(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      implicit none
      integer i,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm),x,a,b,alpha,beta,aux,c
      common /c_wrkparm/p,jp
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'mcmillan pseudopotential -a/r**b'
       write(6,*)'enter a, b'
       read(*,*)p(1),p(2)
       jp=2
      endif
      a=p(1)
      b=p(2)
      x=drt*nt
      c=2*a/x**b
      x=drt*10
      alpha=a*(1+0.5*b)/x**b
      beta=0.5*a*b/x**(b+2)
      do i=i0,nt
       r=r0+i*drt
       if(r.lt.x)then
        f(i)=alpha-beta*r**2
        df(i)=-2*beta*r
       else
        aux=a/r**b
        f(i)=aux
        df(i)=-b/r*aux
       endif
       r=(2*nt-i)*drt
       aux=a/r**b
       f(i)=f(i)+aux-c
       df(i)=df(i)+b/r*aux
      enddo
      return
      end

      subroutine csi(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
c eta function for backflow correlation
      implicit none
      integer i,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm),x,aux,p0,p1,q0,q1,aux1
      common /c_wrkparm/p,jp
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'a*exp(-((r-b**2)*c)**2) '
       write(6,*)'enter a,b,c'
       read(*,*)p(1),p(2),p(3)
       jp=3
      endif
      do i=i0,nt
       r=r0+drt*i
       x=(r-p(2)**2)*p(3)
       aux=p(1)*exp(-x*x)
       q0=aux
       aux=aux*p(3)
       q1=-aux*2.d0*x
       aux=aux*p(3)
       r=max(r,1.d-10)
       x=1.d0/r
       f(i)=q0
       df(i)=q1
      enddo
      aux1=2*f(nt)
      do i=i0,nt
       r=r0+drt*(2*nt-i)
       x=(r-p(2)**2)*p(3)
       aux=p(1)*exp(-x*x)
       q0=aux
       aux=aux*p(3)
       q1=-aux*2.d0*x
       aux=aux*p(3)
       r=max(r,1.d-10)
       x=1.d0/r
       f(i)=f(i)+q0-aux1
       df(i)=df(i)-q1
      enddo
      return
      end

      subroutine gauss(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      implicit none
      integer i,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm)
      common /c_wrkparm/p,jp
      spline_flag=0
      if(stdin_flag.ne.0)then
       write(6,*)'gaussian pseudopotential -c*r**2'
       write(6,*)'enter c'
       read(*,*)p(1)
       jp=1
      endif
      do i=i0,nt
       r=r0+i*drt
       f(i)=p(1)*r**2
       df(i)=2.d0*p(1)*r
      enddo
      return
      end

      subroutine lj2d(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
c lj potential
      implicit none
      integer i,i0,nt,mt,stdin_flag,spline_flag
      real*8 r0,drt,f(0:nt),df(0:nt),r,v,dv
      real*8 vshift,vtail,epsilon,sigma
      common /c_lj/epsilon,sigma
      common /c_v/vshift,vtail
      spline_flag=1
      call v_lj(nt*drt,vshift,dv)
      call r_set(4*(mt+1),f(0),0.d0)
      do i=i0,nt-1
       r=r0+drt*i
       call v_lj(r,f(i),df(i))
       f(i)=f(i)-vshift
      enddo
c     call spline(mt,nt,drt,f)
      i=2
      call tail_correction(r,drt,i,v_lj,vtail)
      return
      end

      subroutine v_lj(r,f,df)
c LJ potential
      implicit none
      real*8 r,f,df,epsilon,sigma,rm1,f6,f12,a,b
      common /c_lj/epsilon,sigma
      if(r.lt.sigma*0.5d0)then
       a=4*epsilon*(2.d0**12-2.d0**6+12*2.d0**11-6*2.d0**5)
       b=16*epsilon/sigma**2*(12*2.d0**11-6*2.d0**5)
       f=a-b*r**2
       df=-2*b*r
      else
       rm1=1.d0/r
       f6=(sigma*rm1)**6
       f12=f6**2
       f=4*epsilon*(f12-f6)
       df=4*epsilon*(-12*f12+6*f6)*rm1
      endif
      return
      end

      subroutine tail_correction(r,drt,ndim,f,tail)
      implicit none
      integer i,ndim
      real*8 r,drt,tail,v,dv,s
      external f
      r=r+0.5d0*drt
      s=0.d0
      do i=1,10000
       call f(r,v,dv)
       s=s+v*r**(ndim-1)
       r=r+drt
      enddo
      if(ndim.eq.1)tail=0.5d0                 *drt*s
      if(ndim.eq.2)tail=0.5d0*2*acos(-1.d0)   *drt*s
      if(ndim.eq.3)tail=0.5d0*4*acos(-1.d0)**2*drt*s
      return
      end

      subroutine mcmillan_parm_liq(rho,p)
      implicit none
      real*8 p(2),rho,a1,b1,c1,a2,b2,c2
	a1 = 269.425071428199
	b1 = -6488.15357141639
	c1 = 108433.928571336
	a2 = 4.88518411709897
	b2 = 15.1652198542003
	c2 = 114.821385767455
      p(1)=a1+b1*rho+c1*rho**2
      p(2)=a2+b2*rho+c2*rho**2
      return
      end
      
      subroutine csi_parm_liq(rho,p)
      implicit none
      real*8 p(3),rho,a3,b3,c3,a4,b4,c4,a5,b5,c5
	a3 = 0.53218869281834
	b3 = 24.7688531318064
	c3 = -76.5713977242002
	a4 = 0.546405871451597
	b4 = -22.5512664293173
	c4 = 242.710357148502
	a5 = 0.525837451085429
	b5 = -2.16323053092823
	c5 = 33.0160590625688
      p(1)=a3+b3*rho+c3*rho**2
      p(2)=a4+b4*rho+c4*rho**2
      p(3)=a5+b5*rho+c5*rho**2
      return
      end

      subroutine mcmillan_parm_sol(rho,p)
      implicit none
      real*8 p(2),rho,a6,b6,c6,a7,b7,c7
	a6 = -834.479142846381
	b6 = 28064.2714282939
	c6 = -143778.571426838
	a7 = 2.72693428574214
	b7 = 93.6078571421391
	c7 = -426.642857138374
      p(1)=a6+b6*rho+c6*rho**2
      p(2)=a7+b7*rho+c7*rho**2
      return
      end
      
      subroutine gauss_parm_sol(rho,p)
      implicit none
      real*8 p(1),rho,a8,b8,c8 
	a8 = 0.499657714272663
	b8 = -16.7376571425207
	c8 = 213.942857140757
      p(1)=a8+b8*rho+c8*rho**2
      return
      end
