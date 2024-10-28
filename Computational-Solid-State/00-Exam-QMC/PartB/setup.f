c setup program for 2d
      implicit none
      integer mnk,mdim,msites,mbasis,mnts
      parameter(mnk=9000,mdim=2,msites=100,mbasis=1,mnts=501  )
      integer nt,np,nup,ndown,len,i,j,ndim,ng,ngfound,ifv0
      integer ifex,iexp,iunit,i_cub(mdim),nbasis
      integer nkappa
      integer jp,m_parm
      parameter(m_parm=20)
      real*8 p(m_parm)
      common /c_wrkparm/p,jp
      real*8 drt,rs,hbs2m,gvect(mdim,mnk),a(mdim,mdim),cut
      real*8 rho,el,l2,tail,d,pi,v0,sites(mdim,msites)
      real*8 aux(mdim),delta,basis(mdim,mbasis)
     &      ,factor,abrav(mdim),shift(mdim),gg,gg0,gstart,gdelta
      real*8 const,konst,self,zelf,usr(mnts,2),vsr(mnts,2)
     &      ,ulr(mnk),vlr(mnk)
      character*30 runid,filename,routinename,te_name
      character word*200
      logical l_ex
      common /c_rs/rs
      external coulomb,u2_ob_cos,csi,polycusp
      pi=acos(-1.d0)
      iunit=10
c
      write(*,*)'run id'
      read(*,'(a)')runid
      open(2,file='runid',status='unknown')
      write(2,'(a)')runid
      close(2)
      write(*,*)'rs'
      read(*,*)rs
      write(*,*)'nup'
      read(*,*)nup
      write(*,*)'ndown'
      read(*,*)ndown
      if(nup.lt.ndown)stop'nup.lt.ndown'
      np=nup+ndown
c
      ndim=2
      len=index(runid,' ')-1
      open(iunit,file=runid(1:len)//'.sy',status='unknown')
      write(word,*)'ndim 2'
      call writestring(word,iunit)
      hbs2m=1.d0/rs**2
      rho=1.d0/pi
      el=((np)/rho)**(1.d0/ndim)

c up particle positions
      write(word,*)'type up',  nup,  hbs2m,' ',runid(1:len)//'.up.x'
      call writestring(word,iunit)
      ifex=1
      filename=runid(1:len)//'.up.x'
c     inquire(file=filename,exist=l_ex)
c     if(l_ex)then
c      write(*,*)'replace existing file ',filename,'? (1/0)'
c      read(*,*)ifex
c     endif
c     if(ifex.ne.0)then
       nbasis=1
       delta=1.d-2
       do i=1,ndim
        basis(i,1)=0.d0
        shift(i)=0.d0
        abrav(i)=1.d0
        i_cub(i)=(float(nup))**(1.d0/ndim)+1
       enddo
       call xtal_sites(mdim,msites,abrav,ndim,basis,nbasis,rho
     &                ,nup,i_cub,shift,delta,aux,sites,factor,1)
       open(2,file=filename,status='unknown')
       do i=1,nup
        write(2,*)(sites(j,i),j=1,ndim)
       enddo
       close(2)
c     endif

c down particle positions
      if(ndown.ne.0)then
       write(word,*)
     & 'type down',  ndown,  hbs2m,' ',runid(1:len)//'.down.x'
       call writestring(word,iunit)
       ifex=1
       filename=runid(1:len)//'.down.x'
       inquire(file=filename,exist=l_ex)
c      if(l_ex)then
c       write(*,*)'replace existing file ',filename,'? (1/0)'
c       read(*,*)ifex
c      endif
c      if(ifex.ne.0)then
        nbasis=1
        delta=1.d-2
        do i=1,ndim
         basis(i,1)=0.d0
         shift(i)=0.5d0
         abrav(i)=1.d0
         i_cub(i)=(float(nup))**(1.d0/ndim)+1
        enddo
        call xtal_sites(mdim,msites,abrav,ndim,basis,nbasis,rho
     &                 ,ndown,i_cub,shift,delta,aux,sites,factor,1)
        open(2,file=filename,status='unknown')
        do i=1,ndown
         write(2,*)(sites(j,i),j=1,ndim)
        enddo
        close(2)
c      endif
      endif

c numero di punti e passo delle tabelle
      nt=500
      drt=0.5d0*el/nt

c coulomb potential
      ifv0=1
      filename=runid(1:len)//'.v'
      call ocpset(np,rs,nt
     &           ,vsr,vlr,self,const,usr,ulr,zelf,konst,nkappa)
      te_name='te.vrk'
      call traduci(te_name,np,nt,drt,v0,ifv0,filename,ng
     &           ,vsr,vlr,self,const,nkappa)
      tail=0 ! somme di ewald
      iexp=1 ! moltiplica la tabella in spazio r per 1/r**iexp
      factor=1 ! fattore moltiplicativo del potenziale
      write(word,*)'v2 up up ',filename,iexp
      call writestring(word,iunit)
      if(ndown.ne.0)then
       write(word,*)'v2 down down ',filename,iexp
       call writestring(word,iunit)
       write(word,*)'v2 up down ',filename,iexp
       call writestring(word,iunit)
      endif
      if(ifv0.ne.0)then
       write(word,*)'v0 ',v0*np
       call writestring(word,iunit)
      endif

c pair pseudopotential
      write(*,*)'2body: yes or no? (1/0)'
      read(*,*)i
      if(i.ne.0)then
       filename=runid(1:len)//'.u'
c     write(*,*)'pseudo short range o rpa  (0/1)'
c     read(*,*)i
       write(*,*)'pseudo rpa'
       i=1
       if(i.eq.0)then
        write(routinename,'(8a,22x)')'polycusp'
        write(*,*)'cusp unpol ',-rs,'  pol ',-rs/3
        call tgen(polycusp,filename,nt,drt,1,routinename)
       else
        te_name='te.urk'
        call traduci(te_name,np,nt,drt,v0,ifv0,filename,ng
     &              ,usr,ulr,zelf,konst,nkappa)
       endif
       write(word,*)'u2 up up ',filename
       call writestring(word,iunit)
       if(ndown.ne.0)then
        write(word,*)'u2 down down ',filename
        call writestring(word,iunit)
        write(word,*)'u2 up down ',filename
        call writestring(word,iunit)
       endif
      endif
c rlv of the simulation box (need k vectors for plane waves)
      cut=2.d0
c     ng=256 ! questo ora e' in output da ocpset
      do i=1,ndim
       do j=1,ndim
        a(j,i)=0.d0
       enddo
       a(i,i)=el
      enddo
c loop with increasing cut until ngfound.ge.ng
      do i=1,10
       call krlv(cut,a,ndim,mdim,gvect,mnk,ngfound,0)
       if(ngfound.ge.ng)go to 1
       cut=cut*2.d0**(1.d0/ndim)
      enddo
c write ng rlv's
    1 filename=runid(1:len)//'.k'
      open(2,file=filename,status='unknown')
      gg0=0.d0
      do i=1,ng
       gg=sqrt(gvect(1,i)**2+gvect(2,i)**2)
       write(2,*)(gvect(j,i),j=1,ndim),gg
      enddo
      close(2)
      write(word,*)'pbc',el,el
      call writestring(word,iunit)
      write(word,*)'kspace ',filename
      call writestring(word,iunit)

c 3body
       write(*,*)'3body: yes or no? (1/0)'
       read(*,*)i
       if(i.ne.0)then
        jp=3
        if(rs.le.1.5)then
         p(1)=0.02
        elseif(rs.le.5.5)then
         p(1)=0.1
        else
         p(1)=0.2
        endif
        p(2)=0.01 ! 0.1
        p(3)=1.0
        filename=runid(1:len)//'.u3'
        write(routinename,'(a3,27x)')'csi'
        call tgen(csi,filename,nt,drt,0,routinename)
        write(word,*)'u3 up up ',filename
        call writestring(word,iunit)
        if(ndown.ne.0)then
         write(word,*)'u3 up down ',filename
         call writestring(word,iunit)
         write(word,*)'u3 down down ',filename
         call writestring(word,iunit)
        endif
       endif
c rhok
      write(word,*)'rhok up'
      call writestring(word,iunit)
      if(ndown.ne.0)then
       write(word,*)'rhok down'
       call writestring(word,iunit)
      endif
c backflow
       write(*,*)'backflow: yes or no? (1/0)'
       read(*,*)i
       if(i.ne.0)then
        jp=3
        p(1)=0.1 ! 0.02 ! 0.2
        p(2)=0.01 ! 0.1
        p(3)=1.0
        filename=runid(1:len)//'.b'
        write(routinename,'(a3,27x)')'csi'
        call tgen(csi,filename,nt,drt,0,routinename)
        write(word,*)'backflow up up ',filename
        call writestring(word,iunit)
        if(ndown.ne.0)then
         write(word,*)'backflow up down ',filename
         call writestring(word,iunit)
         write(word,*)'backflow down down ',filename
         call writestring(word,iunit)
        endif
       else
c onde piane
        write(word,*)'plane-wave up '
        call writestring(word,iunit)
        if(ndown.ne.0)then
         write(word,*)'plane-wave down '
         call writestring(word,iunit)
        endif
       endif

      close(iunit)

      stop
      end

      subroutine ocpset(nparts,rs,nt
     &                 ,vsr,vlr,self,const,usr,ulr,zelf,konst,nkappa)
c program to generate qucu input for the 2-d electron gas
c version of march 7,1990 for qucu2 input
c added extra cosines-4/12/90
      implicit real*8(a-h,o-z)
      parameter (mnp=400,mnkv=9000,mnsh=9000,ndim=2,mnts=501)
      dimension tpiell(ndim),rkcomp(ndim,mnkv),rknorm(mnsh),kmult(mnsh)
      dimension work1(400),work2(5500)
     &         ,vkbare(mnsh),ukbare(mnsh),wtk(mnsh)
      dimension x(ndim,mnp),ncell(3),rv(mnts),vsr(mnts,2),vlr(mnkv)
      dimension usr(mnts,2),ulr(mnkv)
      real*8 konst
      dimension nppss(2),wsites(ndim,mnp)
c%
      dimension gvctr(ndim),wfoval(2)
      character qid*14,filen*14,filenn*14
      common/copenfn/filen
 

c periodic boundary conditions fucntions. ell and el2 must be defined elsewhere.
c returns the displacement vector between -L/2 and +L/2
      dimension tpbc(3) ! (also, tpbc -> t in last line)
      common/cbc/ell(3),el2(3)
      fpbc(d,l)=cvmgp(d-sign(ell(l),d),d,abs(d)-el2(l))
c returns the square of the displacement vector
      fabc(d,l)=(el2(l)-abs(el2(l)-abs(d)))**2
c this is when you know which side one particle is on
      fpbct(d,l)=cvmgp(d-tpbc(l),d,abs(d)-el2(l))

      i=12312
      call ranset(i)
      pi=3.1415 92653 58979d0
      sangle=sang(ndim)
c     write (*,*) 'input name of file'
      qid='te'
50    format(a8)
           ln=index(qid,' ')-1
c     open(1,file=qid(1:ln)//'.sy',status='unknown')
c     rewind(1)
c     open(66,file=qid(1:ln)//'.os',status='unknown')
c     rewind(66)
 
c        write (1,*) 'UNITS ENERGY Rydb'
c        write (1,*) 'UNITS LENGTH a'
 
c     write (*,*) ' input: nparts,  rs '
c     read (*,*) nparts,rs
c      write (66,*)  'rs = ',rs
c      write(66,'(''rs = '',e20.13)')rs
      ro=ndim/sangle
      vol=nparts/ro
c     write (*,*) 'input crystal, number of unit cells in the direction'
c     write (*,*) ' crystal:1=sc,2=bcc,3=hcp,4=fcc(3d),3=triangular(2d)'
c     read (*,*) nxtalp,(ncell(l),l=1,ndim)
      nxtalp=0
      ncell(1)=0
      ncell(2)=0
c%
c     write (*,*) 'input ng'
c     read  (*,*) ng
c     write (*,*) 'input vext'
c     read  (*,*) vext
c     write (*,*) 'input ifmat'
c     read  (*,*) ifmat
c     veff=0.d0
c     if(ifmat.ne.0)then
c        write(*,*)'input veff '
c        read(*,*)veff
c     endif
      ifmat=0
 
      call sites(1,wsites,nparts,ndim,nxtalp,vol,ncell,ell,rnn)
      hbar=1.d0/rs
      enorm=nparts
c     write (1,*) 'PARAMETER ENORM ',enorm,' HBAR ',hbar
c     write 
c    +(1,'(''PARAMETER ENORM '',e20.13,'' HBAR '',e20.13)')enorm,hbar
c     write (1,88) (ell(l),l=1,ndim)
88    format(' PARAMETER BOXSIZE ',3e18.10)
 
      rmass=.5d0
c     write (*,*) ' input: nspins nppss '
c     read (*,*) nspins,(nppss(j),j=1,nspins-1)
      nspins=2
      nppss(1)=nparts/2
      nppss(2)=nparts-nppss(1)
 
      if(nspins.gt.0) then
      nbands=2
      nppss(2)=nparts-nppss(1)
      npmax=max(nppss(1),nppss(2))
 
      else
      nbands=1
      endif
      if(nspins.gt.0) then
c        write (1,'(''TYPE e '',2i4,e20.13,2i4,'' '',a20)')
c    +   nspins,nparts,rmass,nppss(1),nppss(2),qid(1:ln)//'.e_ic'
c        write (1,*) 'TYPE e ',nspins,nparts,rmass
c    +   ,nppss(1),' ',nppss(2),' ',qid(1:ln)//'.e_ic'
c%
         if(ifmat.ne.1)then
 
c  plane wave bands for the electrons
c first a k-point at the origin
c           write (1,*)'DEFINE CONSTANT e'
c then fill up fermi sea
c           write (1,*)'DEFINE  PWAVE e ',npmax-1,' 1 '
c start filling with first wave vector
c%
         endif
 
      else
c        write (1,'(''TYPE e '',2i4,e20.13,'' '',a20)')
c    +   nspins,nparts,rmass,qid(1:ln)//'e_ic'
c     write (1,*) 'TYPE e ',nspins,nparts,rmass,' '
c    +,qid(1:ln)//'.e_ic'
c for bosons localize them
         write (*,*) 'input c for localization'
         read (*,*) cp
      filen=qid(1:ln)//'sites'
      call writecon(wsites,ndim,nparts,ell)
c     write (1,*)'WF GAUSS e ',filen,' cp ',cp
       endif
 
c     lptable=500
c     write(*,*)'lptable'
c      read*,lptable
      lptable=nt+1
       if(lptable.gt.mnts) stop
c cutoff radial tables at 1/2 the smallest box dimension
      cutr=.5d0*ell(1)
      do 910 l=2,ndim
910      cutr=min(cutr,.5d0*ell(l))
c     write (*,*) ' input cutk mpoly argek ,cutr'
      cutk = 10.d0
c this is the order of the polynomial
      mpoly=8
      argek=14.d0
c     write (*,*) ' current values', cutk, mpoly ,argek ,cutr
c     read (*,*) cutk,mpoly,argek,cutr
 
      eps=1.d-6
      csiv=cutr/(float(lptable-1)+eps)
      call ggrid(rv,'LINEAR',lptable,eps,cutr)
c     write(*,*) ' input  nlamb '
c     read (*,*) nlamb
      nlamb=20
c     nlamb=20
      do 66 l=1,ndim
      el2(l)=.5d0*ell(l)
66    tpiell(l)=2*pi/ell(l)
 
c generate the set of kvectors
      call shells(ndim,tpiell,cutk,nshlls,rkcomp,rknorm(2),kmult
     +,nvects,mnkv,mnsh)
      hbs2m=hbar**2/(2.d0*rmass)
      call fpke(nppss,rknorm(2),kmult,vol,hbs2m,enorm,nspins,ndim)
      cutko=max(1.5d0,rknorm(nlamb+2))
 
c     write (1,*) 'PARAMETER LPTABLE ',lptable,' CUTR ',cutr
c     write (1,*) 'PARAMETER NLAMB ',nlamb,' CUTK ',cutko,' IEXPONV ',1
c     write (1,'(''PARAMETER LPTABLE '',i5,'' CUTR '',e20.13)')
c    +lptable,cutr
c     write (1,'(''PARAMETER NLAMB '',i5,'' CUTK '',e20.13,
c    +'' IEXPONV 1'')')nlamb,cutko
 
      call fillk(rknorm,wtk,kmult,nshlls,nshex1,argek,ndim,vol,mnsh) 
 
      e2=2.d0/rs
      if(nspins.eq.0) then
       fermiwv=0.d0
      else
       fermiwv=2*pi*((ndim*nparts)/(sangle*nspins*vol))**(1.d0/ndim)
       endif
c set up potential and rpa pseudopotential
c     write (66,*) ' fermi wavevector = ',fermiwv
 
c put in uniform background by setting k=0 term to zero
      vkbare(1)=0.d0
       ukbare(1)=0.d0
      do 111 k=2,nshex1
      vkbare(k)=(e2*sangle)/(vol*rknorm(k)**(ndim-1))
      if(nspins.gt.0) then
c compute structure factor for ideal fermi gas
      if(rknorm(k).lt.2.d0*fermiwv) then
        y=.5d0*rknorm(k)/fermiwv
        if(ndim.eq.2)ski=1.5708d0/(asin(y)+y*sqrt(1-y*y))
        if(ndim.eq.3)ski=2.d0/(y*(3-y*y))
      else
         ski=1.d0
      endif
      s1=-ski
      s2=ski**2
      else
c for the crystal
      s1=-1-4*cp/rknorm(k)**2
      s2=1.d0+8*cp/rknorm(k)**2
      endif
       ax=4.d0*nparts*rmass*vkbare(k)/(hbar*rknorm(k))**2
111   ukbare(k)=(s1+sqrt(s2+ax))/(2.d0*nparts)
 
      do 312 l=1,lptable
      vsr(l,1)=0.d0
      vsr(l,2)=0.d0
312   continue
 
      if(ndim.eq.3)vmad=-2.837297479d0*e2/vol**(1.d0/3.d0)
      if(ndim.eq.2)vmad=-3.90026492d0*e2/vol**(1.d0/2.d0)
      call mdlng(ndim,1.d0,ell,e2,vmad2)
c     write (66,*) ' ideal square madelung constant',vmad,vmad2
c this prints out the potential fits
      filen=qid(1:ln)//'.vrk'
      call fitpn(mpoly,vkbare,rknorm,wtk,nshex1,nlamb,e2,work1
     +,work2,0,vol,rv,lptable,vsr,vlr,self,const,nkappa,mnts,ndim,1)
 
c     write (1,*) 'POT PAIR RKSPACE e e ',filen
c write out coordinates
c add a random displacement so determinants won't be zero
      delta=.20d0*ell(1)
      do 90 i=1,nparts
      do 90 l=1,ndim
90    x(l,i)=fpbc(wsites(l,i)+delta*(ranf()-.5d0),l)
c     filen= qid(1:ln)//'.e_ic'
c     call writecon(x,ndim,nparts,ell)
 
c these are the symbolic name of the correlation functions
c     write (*,*) ' input uee  '
c     read (*,*) uee
      uee=1.d0
 
      do 1112 l=1,lptable
      usr(l,2)=0.d0
1112  usr(l,1)=0
       cusp=-rs/float(ndim-1)
       filen=qid(1:ln)//'.urk'
       call fitpn(mpoly,ukbare,rknorm,wtk,nshex1,nlamb,cusp,work1
     +,work2,1,vol,rv,lptable,usr,ulr,zelf,konst,nkappa,mnts,ndim,1)
c      write (1,*) 'WF PAIR RKSPACE e e ',filen,' uee ',uee
c      write (1,'(''WF PAIR RKSPACE e e '',a15,'' uee '',e20.13)')
c    + filen,uee
c%
c     write (1,*) 'PARAMETER VEXT ',vext
c     write (1,'(''PARAMETER VEXT '',e20.13)')vext
 
       nextra=0
       nlambp=0
c       write (*,*)' type number of extra functions  and lamb'
c       read (*,*) nextra,nlambp
       nlambp=min(nlambp,nlamb)
       do 110 ni=1,nextra
      ln=index(filen,' ')-1
      filenn=filen(1:ln)//char(48+ni)
c      write (1,*)'WF PAIR RKSPACE e e ',filenn,' uee',ni,'  1.'
      open(21,file=filenn(1:ln+1))
      rewind(21)
c     call wcos(ni,lptable,rv,nlamb,nlambp,vol,rknorm,wtk)
110   continue
       return
       end
      subroutine ranset(idum)
c initialize variables in /cran/ (from ran3(idum) of numerical recipes)
      implicit none
      integer idum,mbig,mseed,mz,ma,mj,mk,i,ii,k,inext,inextp
      real*8 fac
      dimension ma(55)
      common /cran/fac,ma,mbig,mz,inext,inextp
      mbig=1000000000
      mseed=161803398
      mz=0
      fac=1.d0/mbig
c
      mj=mseed-iabs(idum)
      mj=mod(mj,mbig)
      ma(55)=mj
      mk=1
      do i=1,54
         ii=mod(21*i,55)
         ma(ii)=mk
         mk=mj-mk
         if(mk.lt.mz)mk=mk+mbig
         mj=ma(ii)
      enddo
      do k=1,4
         do i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.mz)ma(i)=ma(i)+mbig
         enddo
      enddo
      inext=0
      inextp=31
      return
      end
      function sang(ndim)
      implicit real*8(a-h,o-z)
      save
c solid angle in ndim dimensions
      pi=3.1415 92653 58979d0
      sang=2*pi*(ndim-1)
      return
      end
      subroutine sites(ir,x,natoms,ndim,nxtalp,vol,ncell,ell,rnn)
      implicit real*8(a-h,o-z)
      save
      dimension x(ndim,natoms),ncell(3),ell(3),q(3,16),icount(3)
      dimension b(3)
c******************************************************
c  sites computes natoms crystal sites and puts them in x
c  ndim is the spatial dimensionality
c  nxtal is the crystal type see 1,2,3,4 below
c      if zero will chose crystal type to minimize number of vacancies
c  ro is the natoms per unit volume-used to compute ell
c  ncell(3) are the number of unit cells in each direction
c    if product(ncell)*npuc.lt.natoms ncell is increased so that
c    there at least as many lattice sites as particles and for hcp
c    the box is roughly cubic
c ell(3) --computed--is the size of the simulation box =ncell*cell size
c*********************************************************************
      ro=natoms/vol
c     if(ir.ne.0)write (66,15) natoms,ndim,ro
15    format(/' computing',i4,' lattice sites  dimensionality',i2
     +,' number density ',e12.5)
      do 16 ii=1,8
      do 16 l=1,3
16    q(l,ii)=0.0d0
      do 17 l=1,ndim
17    b(l)=1.0d0
      nxtal=iabs(nxtalp)
      if(nxtal.gt.0) go to 45
c determine lattice type by minimizing 2**l*ncell**ndim-natoms
      nvac=natoms
      do 460 l=1,ndim
      npuc=2**(l-1)
      nc=int((dble(natoms)/dble(npuc))**(1.d0/ndim)+1.0d0-small)
      nvact=npuc*nc**ndim-natoms
      if(nvact.ge.nvac) go to 460
      nvac=nvact
      nxtal=npuc
460   continue
c     if(ir.ne.0)write (66,461)
461   format('  crystal type chosen by default to minimize vacancies ')
45    npuc=nxtal
 
       if(nxtal.eq.1) then
c     if(ir.ne.0)write (66,*)'  simple cubic lattice'
 
      elseif(nxtal.eq.2) then
c     if(ir.ne.0)write (66,*)'  body centered cubic lattice'
      d=.5d0
      do 112 l=1,ndim
112   q(l,2)=d
 
      elseif(nxtal.eq.3.or.nxtal.eq.9) then
c     if(ir.ne.0)write (66,*)'   hexagonal close packed lattice'
      npuc=4
      if(ndim.eq.2) npuc=2
      if(ndim.eq.2) b(1)=3.d0**(-0.25d0)
      if(ndim.eq.3) b(1)=1.d0/sqrt(2.d0)
      b(2)=b(1)*sqrt(3.d0)
      b(3)=b(2)*sqrt(2.d0)/1.5d0
      q(1,2)=0.5d0
      q(2,2)=0.5d0
      q(1,3)=.5d0
      q(2,3)=5.d0/6.d0
      q(3,3)=0.5d0
      q(1,4)=1.0d0
      q(2,4)=1.d0/3.d0
      q(3,4)=.5d0
 
      if(nxtal.eq.3) go to 10
c     if(ir.ne.0)write(66,*)'mhcp lattice'
      npuc=8
c scale by the c/a ratio
      factor=(x(2,1)*sqrt(.375d0))**(-1.d0/3.d0)
      b(1)=b(1)*factor
      b(2)=b(2)*factor
      b(3)=b(3)/factor**2
      ds=x(1,1)*ro**(1.d0/3.d0)/(4.d0*b(3))
      do 512 is=1,4
      do 512 l=1,3
      qq=q(l,is)
      if(l.le.2) then
        d=0.d0
      else
        d=ds
      endif
      q(l,is)=qq-d
512   q(l,is+4)=qq+d
c b is the size of the unit cell--volume of unit cell is one
c q(npuc,ndim)*b(ndim) are the vector dispacements of sites withincell
 
      elseif(nxtal.ge.4.and.nxtal.le.7) then
c     if(ir.ne.0)write (66,*)'  face-centered cubic lattice'
      do 114 ii=1,3
      do 114 l=1,3
114   if(ii.ne.l) q(l,ii+1)=.5d0
      if(nxtal.eq.4) go to 10
      npuc=8
      if(nxtal.eq.6) go to 810
      if(nxtal.eq.5) go to 801
 
c now add alpha nitrogen displacements
c     if(ir.ne.0)write (66,*)' alpha nitrogen fcc lattice'
      ds=x(1,1)*ro**(1.d0/3.d0)/sqrt(48.0d0)
      do 500 l=1,3
      do 500 is=1,4
      d=-ds
      if(q(l,is).gt..2d0) d=ds
      qq=q(l,is)
      q(l,is)=qq+d
500   q(l,is+4)=qq-d
      go to 10
 
c     if(ir.ne.0) write(66,*)' diamond lattice'
801   continue
      ds=.125d0
      do 803 l=1,3
      do 803 is=1,4
      qq=q(l,is)
      q(l,is)=qq-ds
803   q(l,is+4)=qq+ds
      go to 10
 
c     if(ir.ne.0) write (66,*)' diamond bond lattice'
810   continue
      npuc=16
      k=4
      ds=.25d0
      do 813 is=1,4
      do 813 it=1,3
      k=k+1
      do 813 l=1,3
      q(l,k)=q(l,is)
813   if(l.ne.it)q(l,k)=q(l,k)+ds
 
      elseif(nxtal.eq.8) then
c     if(ir.ne.0)write(66,*)'simple hexagonal lattice'
      npuc=  2
      alpha=1.d0/(x(2,1)*sqrt(3.d0))**(1.d0/3.d0)
      b(1)=alpha
      b(2)=sqrt(3.d0)*alpha
      b(3)=x(2,1)*alpha
      do 910 l=1,2
910   q(l,2)=.5d0
 
      endif
10    continue
c     if(ir.ne.0)write (66,505) ((q(l,is),l=1,ndim),is=1,npuc)
505   format(' q displs ' 3f10.5)
      npts=1
      do 20 l=1,ndim
20    npts=npts*ncell(l)
      if(npts*npuc.ge.natoms) go to 30
c recalculate ncell since there are too few lattice points
      xcell=(dble(natoms)/dble(npuc))**(1.d0/dble(ndim))
      npts=1
      do 56 l=1,ndim
      ncell(l)=int((xcell/b(l))+1.0d0-small)
56    npts=npts*ncell(l)
30    a=(vol/npts)**(1.d0/ndim)
      do 31 l=1,ndim
      icount(l)=0
31    ell(l)=ncell(l)*a*b(l)
      nvac=npts*npuc-natoms
c     if(ir.ne.0)write (66,255) npuc,nvac,(ell(l),l=1,ndim)
255   format(' npuc ',i2,' nvacancies ',i5,' box size',3e12.5)
c     if(ir.ne.0)write (66,256) (ncell(l),l=1,ndim)
256   format(' number of cells in each direction',3i5)
c put in nvac vacancies by flags in x array
      data small/1.0d-3/
      do 5 i=1,natoms
5     x(1,i)=0.d0
      if (nvac.eq.0) go to 9
c compute index  for skipping
      do 7 i=1,nvac
      j=int((natoms-1)*ranf()+2.d0)
c note that the first site will always be filled
7     x(1,j)=x(1,j)+1.d0
9     continue
      i=1
c  the particles are confined abs(x(l)).le.el2(l)
      ellsq=ell(1)**2
      rnn=ellsq
      rn2=ellsq
c   loop over the different points inthe unit cell
      do 451 is=1,npuc
c loop over all the unit cells
      do 451 j=1,npts
c skip over lattice site if vacant
      if(x(1,i).lt.small) go to 8
      x(1,i)=x(1,i)-1.d0
      go to 300
8     rsq=0.0d0
      do 88 l=1,ndim
      xt=a*b(l)*(icount(l)+q(l,is))
      if(xt.ge.0.5d0*ell(l)) xt=xt-ell(l)
      if(xt.lt.-.5d0*ell(l)) xt=xt+ell(l)
      x(l,i)=xt
88    rsq=rsq+(xt-x(l,1))**2
      i=i+1
c find the nearest and next nearest distance from the origin
c exclude the origin
      if(rsq.lt.small*rnn) go to 300
c is this point greater than previously found 2nd minima
      if(rsq-rn2.gt.-small*rnn) go to 300
      if(abs(rsq-rnn).lt.rnn*small) go to 300
c rnn and rsq must be first and second smallest distances.
      rn2=rnn
      if(rsq.gt.rn2) rn2=rsq
      if(rsq.lt.rnn) rnn=rsq
300    continue
      do 101 l=1,ndim
      icount(l)=icount(l)+1
      if(icount(l).lt.ncell(l)) go to 450
101   icount(l)=0
450   continue
451   continue
      rnn=sqrt(rnn)
      rn2=sqrt(rn2)
      if(ir.ne.0)then
c      write (66,100) rnn,rn2
100   format(' nearest and next nearest neighbor distance  ',2e12.5/)
c      do 110 i=1,natoms
c110    write (66,111) i,(x(l,i),l=1,ndim)
c111    format(' lattice sites ',i5,3e12.5)
       endif
      return
      end
      subroutine writecon(x,ndim,ncomps,ell)
      implicit real*8(a-h,o-z)
      save
c due to peculiarities of UNICOS filename must be passed thru copenfn
      dimension x(ndim,ncomps),ell(ndim)
      character file*14
      common/copenfn/file
      ln=index(file,' ') -1
      open(21,file=file(1:ln),status='unknown',form='formatted')
      rewind(21)
      write (21,*) 'RANK ',2,' ',ndim,' ',ncomps
      write (21,22) (ell(l),l=1,ndim)
22    format(' SIZE ',3e18.10)
      write (21,*) 'BEGIN  coordinates'
      do 10 i=1,ncomps
10    write (21,*) (x(l,i),l=1,ndim)
      close(21)
      return
      end
      subroutine ggrid(r,style,n,r0,r1)
      implicit real*8(a-h,o-z)
      save
c generates a mesh of "style"; n=number of points
c r0 is first point r1 last point
      dimension r(n)
      character style*(*)
 
      if(style.eq.'LINEAR') then
       dr=(r1-r0)/(n-1)
       do 10 i=1,n
10     r(i)=r0+(i-1)*dr
 
      elseif(style.eq.'LOG') then
       dr=(r1/r0)**(1.d0/(n-1))
       rr=r0
       do 20 i=1,n
       r(i)=rr
20     rr=rr*dr
      endif
      return
      end
      subroutine shells(ndim,a,cut,nshlls,rkcomp,rknorm,kmult
     +,nvects,mnkv,mnsh)
      implicit real*8(a-h,o-z)
      save
      dimension a(ndim),rkcomp(ndim,mnkv),rknorm(mnsh),kmult(mnsh)
     +, icount(3),nkspan(3),x(3)
c computes the vectors x(ndim)=(a(1)*n(1),..,a(ndim)*n(ndim))
c where n(i) are integers and x(1)**2+..+x(ndim)**2.le.cut**2
c the vectors x(i) are stored in rkcomp in
c the order given by the values of their norms.
c  also nshlls gives the number of
c different values of the norms ( the relative square norms
c differ by less than 1.e-5) and cc(lknorm+i) gives
c these nshlls norms and kmult(i) is the last vector
c whose norm is given by knorm(i). hence the total
c number of vectors of magnitude less than cut is kmult(nshlls).
      c2=cut**2
      npts=1
      do 2 l=1,ndim
      nkspan(l)=int(0.00001d0+abs(cut/a(l)))
c range of search is (-nkspan(l),+nkspan(l))
      icount(l)=-nkspan(l)
2     npts=(2*nkspan(l)+1)*npts
      nvects=0
      do 3 i=1,npts
      rsq=0.0d0
c nzero will be the first nonzero entry in icount
      nzero=0
      do 4 l=1,ndim
      if(nzero.eq.0)nzero=icount(l)
      x(l)=icount(l)*a(l)
      rsq=rsq+x(l)**2
      if(rsq.gt.c2) go to 30
4      continue
c we only take half of the vectors and exclude the one at the origin
      if(nzero.le.0) go to 30
c     if(nvects.gt.mnkv)
c    +call mcheck(nvects,mnkv,'nvects','mnkv','shells')
c we have found a vector
      nvects=nvects+1
c go thru previous vectors. if they have a greater norm move
c them up one slot
      do 6 j=1,nvects-1
      kj=j
       rks=0.d0
      do 66 l=1,ndim
66    rks=rks+rkcomp(l,j)**2
6     if(rks.ge.rsq) go to 7
      kj=nvects
7     do 8 jp=nvects,kj+1,-1
      do 8 l=1,ndim
8     rkcomp(l,jp)=rkcomp(l,jp-1)
      do 9 l=1,ndim
9     rkcomp(l,kj)=x(l)
c increase counters with carries for the next vector
30     do 31 l=1,ndim
      icount(l)=icount(l)+1
      if(icount(l).le.nkspan(l)) go to 3
31     icount(l)=-nkspan(l)
3     continue
      nshlls=0
c count number of different norms and find pointers
      rnow=0.d0
      do 51 i=1,nvects
       rsq=0.d0
       do 514 l=1,ndim
514    rsq=rsq+rkcomp(l,i)**2
      if(rsq-rnow.gt..001d0*rnow)nshlls=nshlls+1
c     if(nshlls.gt.mnsh)
c    +call mcheck(nshlls,mnsh,'nshlls','mnsh','shells')
      rnow=rsq
      rknorm(nshlls)=sqrt(rnow)
51    kmult(nshlls)=i
      return
      end
      subroutine fpke(nppss,rknorm,kmult,vol,hbs2m,enorm,nspins,ndim)
      implicit real*8(a-h,o-z)
      save
      dimension rknorm(2),kmult(2),nppss(nspins)
c determine the infinite system and finite system free particle kinetic energies
      pi=3.14159265d0
c nppss=number of states occupied. rknorm=list of k magnitudes
c kmult are the multiplicites of shells (as order by shell)
c vol is the volume of the box. hbs2m=hbar**2/2*mass, enorm is the
c energy conversion. nspins=number of spin states, ndim= dimensionality
      sangle=2*pi*(ndim-1)
      tktot=0.d0
      tkitot=0.d0
      do 1 i=1,nspins
      tkf=0.d0
c fill up one at the origin
      needed=nppss(i)-1
      do 2 ks=1,300
      if(ks.eq.1) mult=2*kmult(1)
      if(ks.gt.1) mult=2*(kmult(ks)-kmult(ks-1))
      mult=min0(needed,mult)
      tkf=tkf+hbs2m*rknorm(ks)**2*mult
      needed=needed-mult
      kl=ks
c     write (66,*) ' shell',ks,'mult ',mult,' k ',rknorm(ks)
      if(needed.le.0) go to 3
2     continue
      write (*,*)' too few states in fpke '
      stop
3     continue
      fermiwv=2*pi*(ndim*nppss(i)/(vol*sangle))**(1.d0/ndim)
      tkinf=hbs2m*sangle*vol*
     +       fermiwv**(ndim+2)/((2.d0*pi)**ndim*(ndim+2))
c      write (66,*)'  spin ',i,' ke finite ',tkf,' infinite ',tkinf
c      write (66,*) ' fermiwv ',kl,rknorm(kl),fermiwv
      tktot=tktot+tkf
      tkitot=tkitot+tkinf
1     continue
      tktot=tktot/enorm
      tkitot=tkitot/enorm
c     write (66,*)' free particle ke finite ',tktot,' infinite ',tkitot
c    +,' difference ',tktot-tkitot
      end
      subroutine fillk(rknorm,wtk,kmult,nshlls,nshex1,argek,ndim,vol
     +,mnsh)
      implicit real*8(a-h,o-z)
      save
c sets up rknorm wtk for fitp routines by extending the grid to larger k's
      dimension rknorm(mnsh),wtk(mnsh),kmult(mnsh)
      rknorm(1)=0.d0
c this is the weight of a given k vector
      wtk(1)=1.d0
      wtk(2)=2*kmult(1)
      do 110 k=2,nshlls
110   wtk(k+1)=2*(kmult(k)-kmult(k-1))
c fill in region between cut and argek with a uniform grid
      dk=rknorm(2)*.25d0
      cutk=rknorm(nshlls)
      nshex=nshlls+(argek-cutk)/dk
c     call mcheck(nshex+1,mnsh,'nshex','mnsh','fillk')
      pi=3.1415 92653 58979d0
      con=vol*sang(ndim)/(ndim*(2*pi)**ndim)
      vdown=1+2*kmult(nshlls)
      do 150 k=nshlls+1,nshex
      vup=con*(cutk+dk*(k-nshlls))**ndim
      rknorm(k+1)=cutk+dk*(k-nshlls-.5d0)
      wtk(k+1)=vup-vdown
150   vdown=vup
      nshex1=nshex+1
      return
      end
      subroutine mdlng(ndim,p,ell,e2,v)
      implicit real*8(a-h,o-z)
      save
      dimension ell(ndim),icount(3),nkspan(3)
      pi=3.1415926535d0
      eps=1.d-9
      limit=1+sqrt(-log(eps)/pi)
c omputes the madelung constant for r**-p interaction in ndim dim
      vol=1.d0
      do 1 l=1,ndim
      vol=vol*ell(l)
1     icount(l)=-limit
      alpha=pi*vol**(-2.d0/ndim)
      v=-2*alpha**(.5d0*p)/(p*gam(.5d0*p))
     +-2*pi**(.5d0*ndim)/((ndim-p)*vol*gam(.5d0*p)*
     +alpha**(.5d0*(ndim-p)))
      do 2 i=1,(2*limit+1)**ndim
      r2=0.d0
      rk2=0.d0
      do 3 l=1,ndim
      r2=r2+(icount(l)*ell(l))**2
3     rk2=rk2+(2*pi*icount(l)/ell(l))**2
      if(r2.gt.0)  then
      p2=.5d0*p
      call gammi(gr,p2,alpha*r2,gr0)
      call gammi(gk,.5d0*(ndim-p),rk2/(4.d0*alpha),gk0)
      v=v+pi**(.5d0*ndim)*(4.d0/rk2)**(.5d0*(ndim-p))*gk/(vol*gr0)
     +    +gr/(gr0*r2**(p2)) 
      endif
      do 5 l=1,ndim
      icount(l)=icount(l)+1
      if(icount(l).le.limit) go to 2
5     icount(l)=-limit
2     continue
      v=v*e2
      return
      end
      subroutine fitpn(n,v,rk,wt,nk,nf,b0,b,hs,i0,vol,r,lptable
     +,vsr,vlr,self,const,nkappa,mnts,ndim,ifcon)
      implicit real*8(a-h,o-z)
      save
c version of 3/7/90 DMC
      dimension v(9),rk(9),wt(9),b(9),hs(n,9),vsr(mnts,2),r(lptable)
      dimension vlr(*)
c fits best spherical polynomial b to v(k); hs must have dim n*(n+2)+np
c first term is r**(i0-1). b0 is the value of r**(-1+2*i0) term.
c vsr is a preexisting short range potential and its derivative
c ifcon=1 use the constraint b0 otherwise not
      a=r(lptable)
      pi=3.1415 92653 58979d0
      pn=2*pi*(ndim-1)*a**ndim/vol
c np is the number of factors of (r-a) we multiply by
       np=3 
c ifalt is a flag to decide which way to compute the lhs matrix elements
       ifalt=0
c m is the total number of free parameters
      m=n-np
      if(ifcon.ne.1) then
       mp=m
       i2=1
       beta1=0.d0
       beta2=0.d0
      else
c one more equation because of cusp condition
       mp=m-1
       i2=2
c We are setting up cusp condition constraints
       anp=1-2*mod(np,2)
       if(i0.ne.0) then
        beta1=-anp*b0*a/dble(np)
        beta2=1.d0/dble(np)
       else  
         beta1=anp*b0/a
         beta2=0.d0
       endif 
      endif
c     write (66,*) ' beginning fitpn n nk nf i0= ',n,nk,nf,i0
c     write (66,*) ' ifcon ndim  ',ifcon,ndim
c     write (66,*) ' mnts lptable vol b0  ',mnts,lptable,vol,b0
c     write (66,*) ' a pn np m mp i2 beta ',a,pn,np,m,mp,i2,beta1,beta2
c zero right and left hand side of fitting equation
      do 1 i=1,n+2
      do 1 j=1,n
1     hs(j,i)=0.d0
      chisq=0.d0
c go over k values larger than those explicitly used
      do 20 k=nf+2,nk
c the matrix elements are different in 2 and 3 dimensions
      if(ndim.eq.3)call plint3(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      if(ndim.eq.2)call plint2(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      vv=v(k)-beta1*hs(1,m+2)
      chisq=chisq+wt(k)*vv**2
c for the derivative constraint
      hs(2,m+2)=hs(2,m+2)+beta2*hs(1,m+2)
c add to right hand side
      do 20 i=i2,m
      hs(i,m+1)=hs(i,m+1)+wt(k)*vv*hs(i,m+2)
c add to left hand side
      if(ifalt.eq.0) then
      do 24 j=i2,m
24    hs(j,i)=hs(j,i)+wt(k)*hs(i,m+2)*hs(j,m+2)
      endif
20    continue
      if(ifalt.eq.1) then
c this is an experiment to compute lhs more directly
      do 90 k=1,nf+1
      if(ndim.eq.3)call plint3(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      if(ndim.eq.2)call plint2(rk(k)*a,m,hs(1,m+2),i0,pn,np)
      do 90 i=i2,m
      do 90 j=i2,m
90    hs(j,i)=hs(j,i)-wt(k)*hs(i,m+2)*hs(j,m+2)
      do 95 i=i2,m
      call prodi(hs(1,m+2),pn,ndim,m,i0,np,i)
      do 95 j=i2,m
95    hs(j,i)=hs(j,i)+hs(j,m+2)
      endif
 
c invert right hand side
      call spoco(hs(i2,i2),n,mp,rcond,hs(1,m+2),info)
      if(info.ne.0) then
       write (*,*) ' trouble in fitpn info=',info
       stop
      endif
c make a spare copy of right hand side
      do 40 i=i2,m
40    b(i)=hs(i,m+1)
c solve linear equations
      call sposl(hs(i2,i2),n,mp,hs(i2,m+1))
      do 50 i=i2,m
50    chisq=chisq-b(i)*hs(i,m+1)
c     if(chisq.gt.0) then
c         write (66,*) ' rms error ',sqrt(chisq)
c      else
c         write (66,*) ' chisq negative ? ',chisq
c      endif
c this is cusp constraint
      if(ifcon.eq.1)hs(1,m+1)=beta1+beta2*hs(2,m+1)
c subtract effect of short range potential on fourier components
      do 60 k=1,nk
      if(ndim.eq.3)call plint3(rk(k)*a,m,hs,i0,pn,np)
      if(ndim.eq.2)call plint2(rk(k)*a,m,hs,i0,pn,np)
      do 60 i=1,m
60    v(k)=v(k)-hs(i,m+1)*hs(i,1)
c     write (66,6) n,i0,rcond
6     format(/' fitp m=',2i5,' rcond ',e12.4)
      do 69 i=1,m
69    b(i)=hs(i,m+1)
c     write(66,7) (b(i),i=1,m)
c     if(nf.gt.0)write (66,*) ' k=   ',(rk(k),k=1,nf+1)
c     if(nf.gt.0)write (66,*)' f=   ',(v(k),k=1,nf+1)
7     format(' poly= ',5e14.6)
cvax71    format(' k=    ',10f14.6)
cvax72    format(' f=    ',10e14.6)
c take the derivative of a polyonomial
       call polyd(m,b,hs)
      iexp=1-i0
      ai=a**iexp
c write out table of potential and first derivative
      do 120 i=1,lptable
c avoid the origin for singular potentials
      x=r(i)/a
      v1=stfun(x,b,m)
      v2=stfun(x,hs,m)
      if(iexp.ne.0) then
        v1=v1*ai
        v2=v2*ai
       endif
       if(np.gt.0) then
       v2=v2*(x-1.d0)**np+np*(x-1.d0)**(np-1)*v1
       v1=v1*(x-1.d0)**np
       endif
       vsr(i,1)=vsr(i,1)+v1
       vsr(i,2)=vsr(i,2)+v2/a
120    continue
c     call writetb(vsr,lptable,mnts,2,1,0,a,-1)
c     write (21,*) 'RANK 3 1 1 1 '
c     write (21,*) 'BEGIN self-energy'
c this madelung constant is just for checking in some cases
c (cubic lattice 1/r potential)
      vmad=(-1.d0)**np*(b(2)-np*b(1))
c       write (21,*) vmad
        self=vmad
      do 140 k=1,nk
140   vmad=vmad+wt(k)*v(k)
c     write(66,*) ' computed madelung constant',vmad
c     write (21,*) 'RANK 3 1 1 1 '
c     write (21,*) 'BEGIN constant energy'
c     write (21,*) v(1)
      const=v(1)
      nkpair=0
      do 123 k=2,nf+1
123   nkpair=nkpair+.5d0*wt(k)+1.d-4
c     write (21,*) 'RANK 3 ',nkpair,' 1 1 '
c     write (21,*) 'BEGIN k-space energy'
       nkappa=0
       do 121 k=2,nf+1
c expand out with multiplicity of each k-point
       mult=.5d0*wt(k)+1.d-4
       do 121 iult=1,mult
        nkappa=nkappa+1
        vlr(nkappa)=v(k)
c      write(21,*)v(k)
121    continue
c121    write(21,*) mult,'*',v(k)
      close(21)
      return
      end
      function ranf()
c generates a uniform deviate in (0,1) (from ran3 on numerical recipes)
      implicit none
      integer mbig,mz,ma,mj,inext,inextp
      real*8 ranf,fac
      dimension ma(55)
      common /cran/fac,ma,mbig,mz,inext,inextp
c
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz)mj=mj+mbig
      ma(inext)=mj
      ranf=1.d0-mj*fac
c
      return
      end

      subroutine mathieu(nppss,v,gvctr,tpiell,hbs2m,iwrite)
      implicit real*8(a-h,o-z)
c lowest-energy configuration and the coefficients of mathieu functions
c [2d; eff.pot. v(x)=v*cos(ng*tpiell(1)*x)]
 
      parameter (np=201)
      dimension tpiell(2)
      dimension cmf(np,np),icmf(np,np),imf(0:np,2)
      dimension a(np,np),d(np),e(np)
      dimension l(np),k(np),zero(np)
      save
c     character filen*14
c     common/copenfn/filen
 
      do i=1,np
         do j=1,np
            a(j,i)=0.d0
            cmf(j,i)=0.d0
            icmf(j,i)=0
         end do
         imf(i,1)=0
         imf(i,2)=0
         l(i)=0
         k(i)=0
         d(i)=0.d0
         e(i)=0.d0
         zero(i)=0.d0
      end do
      imf(0,1)=0
      imf(0,2)=0
 
c constants
      data small/.0001d0/
      sq2i=1.d0/sqrt(2.d0)
      ng=nint(gvctr/tpiell(1))
      e0=tpiell(2)**2*hbs2m
      v2=v/2
 
c matrix elements and diagonalization
c     m0=max(2,nint(1.d0+nppss/2.d0/ng),nint(sqrt(v/e0/small)/ng))
      m0=max(10,nint(1.d0+nppss/2.d0/ng),nint(sqrt(v/e0/small)/ng))
      do m=m0,0,-1
         m2p1=m*2+1
         n2p1=ng*m2p1
         if(n2p1.le.np)go to 3
      end do
3     n=(n2p1-1)/2
      np1=n+1
      kplus=np1
      kmins=np1
      imf(1,1)=kplus
      do i=1,n
         kplus=kplus+1
         kmins=kmins-1
         imf(2*i,1)=kplus
         imf(2*i+1,1)=kmins
      end do
      kplus=m+1
      kmins=m+1
      k(1)=kplus
      do i=1,m
         kplus=kplus+1
         kmins=kmins-1
         k(2*i)=kplus
         k(2*i+1)=kmins
      end do
      do ig=1,ng
         istart=isign(ig/2,-mod(ig,2))
         kplus=istart
         kmins=istart
         l(1)=istart
         d(k(1))=e0*istart**2
         do i=1,m
            i2=i*2
            i2p1=i2+1
            kplus=kplus+ng
            kmins=kmins-ng
            l(i2)=kplus
            l(i2p1)=kmins
            d(k(i2))=e0*kplus**2
            d(k(i2p1))=e0*kmins**2
         end do
         do i=2,m2p1
            e(i)=v2
         end do
         do i=1,m2p1
            do j=1,m2p1
               a(j,i)=0.d0
            end do
            a(i,i)=1.d0
         end do
         call tqli(d,e,m2p1,np,a)
         ishift=m2p1*(ig-1)
         do i=1,m2p1
            zero(i+ishift)=d(i)
            do j=1,m2p1
               cmf(l(j)+np1,i+ishift)=a(k(j),i)
            end do
         end do
      end do
      call indexx(n2p1,zero,l)
      do i=1,n2p1
         d(i)=zero(l(i))
         do j=1,n2p1
            a(j,i)=cmf(imf(j,1),l(i))
         end do
      end do
 
c lowest nppss one-particle levels -- start with zero pw 
      do i=1,nppss
         e(i)=d(i)
         imf(i,1)=i
         imf(i,2)=1
      end do
c                                  -- add pairs of pw 
      imf(0,2)=nppss/2
      imf(0,1)=nppss
      do jpw=1,imf(0,2)
         epw=e0*jpw**2 
         do jmf=1,imf(0,1)
            etot=d(jmf)+epw 
            do i=1,nppss+1
               iffind=0
               if(etot.le.e(i))then
                  iffind=1
                  do j=nppss+1,i+2,-1
                     e(j)=e(j-2)
                     imf(j,1)=imf(j-2,1)
                     imf(j,2)=imf(j-2,2)
                  end do
                  e(i)=etot
                  e(i+1)=etot
                  imf(i,1)=jmf
                  imf(i+1,1)=jmf
                  imf(i,2)=jpw*2
                  imf(i+1,2)=imf(i,2)+1
                  go to 1
               endif
            end do
1           continue 
            if(iffind.eq.0)go to 2
         end do
2        continue 
      end do
 
c number of needed mathieu functions and plane waves
      imf(0,1)=0
      imf(0,2)=0
      do i=1,nppss
         imf(0,1)=max0(imf(0,1),imf(i,1))
         imf(0,2)=max0(imf(0,2),imf(i,2))
      end do
 
c basis transf.
      do i=1,imf(0,1)
         do j=2,n2p1,2
            jp1=j+1
            x=a(j,i)
            a(j,i)=(a(j,i)-a(jp1,i))*sq2i
            a(jp1,i)=(x+a(jp1,i))*sq2i
         end do
      end do
c     if(ng.ne.1)then
       do i=3,imf(0,1)
          x=abs(d(i)-d(i-1))
          if(x.lt.1.d-10)then
             im1=i-1
             do j=2,n2p1
                x=a(j,i)
                a(j,i)=(a(j,i)+a(j,im1))*sq2i
                a(j,im1)=(x-a(j,im1))*sq2i
             end do
          endif
       end do
c     endif
 
c needed coefficients and pointers
      ncmf=0
      mcmf=0
      do i=1,imf(0,1)
         jj=0
         do j=1,n2p1
            x=abs(a(j,i))
            if(x.gt..00001d0)then
               jj=jj+1
               icmf(jj,i)=j
               mcmf=max0(mcmf,j)
            endif
         end do
         ncmf=max0(ncmf,jj)
      end do
      do i=1,imf(0,1)
         do j=1,mcmf
            cmf(j,i)=a(icmf(j,i),i)
         end do
      end do
 
      if(iwrite.eq.1)then
         q=2.d0*v/hbs2m/gvctr**2
         write (66,'(''m0      '',i10,/
     +              ''m       '',i10,/
     +              ''q       '',f14.3,/
     +              ''veff    '',f14.3,/
     +              ''ng      '',i10)')m0,m,q,v,ng
         write (66,'(/''configuration'')')
         do i=1,nppss
            write (66,'(i10,i5)')imf(i,1),imf(i,2)
         end do
         x=abs(e(nppss+1)-e(nppss))
         write (66,'(/''to next one particle level: '',e15.3,'' Ry'')')x
         write (66,'(/''eigenvalues, pw energies'')')
         do i=1,imf(0,1)
            x=e0*(i/2)**2
            write (66,'(i4,2f15.5)')i,d(i),x
         end do 
         write (66,'(/''eigenvectors'')')
         do i=1,imf(0,1)
            write (66,'(i4,''    ('',f10.5,'')'')')i,d(i)
            write (66,'(f22.5)')a(1,i)
            write (66,'(2f14.5)')(a(j,i),j=2,mcmf)
         end do
 
         x=1.d0*sqrt(2.d0)
         gx=tpiell(1)*x
         vx=v*cos(gx*ng)
         e(1)=1.d0/sqrt(2.d0)
         do j=1,mcmf/2
            j2=j*2
            j2p1=j2+1
            gxj=gx*j
            e(j2)=sin(gxj)
            e(j2+1)=cos(gxj)
         end do
         do i=1,imf(0,1)
            f=0.d0
            ddf=0.d0
            do j=1,mcmf
               gj2=(tpiell(1)*(j/2))**2
               f=f+a(j,i)*e(j)
               ddf=ddf-a(j,i)*gj2*e(j)
            end do
            d(i)=-hbs2m*ddf+(vx-d(i))*f
         end do
         write (66,'(/''check'')')
         do i=1,imf(0,1)
            write (66,'(i4,f15.5)')i,d(i)
         end do
 
      endif
 
      return
      end
      subroutine mcheck(iv,mv,civ,cmv,sub)
      implicit real*8(a-h,o-z)
      save
      character civ*(*),cmv*(*),sub*(*)
      write (77,1) iv,mv,civ,cmv,sub
1     format(2i8,3a10)
      if(iv.le.mv) return
      write(*,2) sub,civ,cmv,iv,mv
2     format(' memory overflow in ',a,' variable ',2a,' values ',2i7)
      stop
      end
      function gam(a)
      implicit real*8(a-h,o-z)
      save
      dimension p(7),q(7),p1(5,2)
      save p,q,p1
      data(p(i),i=1,7)/.2514886580600251d+05,.5298348466186016d+04,
     1.6177655268060726d+04,.2509926126029017d+03,.5785107455981657d+03,
     2-.2316113056472773d+02,.2021018352970918d+02/
      data(q(i),i=1,7)/.1000000000000000d+01,.2513925177055871d+05,
     1.1985721823627555d+05,-.7373357770095750d+04,
     x -.5047154372852621d+03,
     2.3632214014257158d+03,-.3119624946361091d+02/
      z=a
      fctr=1.0d0
303   if(z.ge.1.0d0) go to 304
      fctr=fctr/z
      z=z+1.0d0
      go to 303
304   if(z.le.2.0d0) go to 301
      z=z-1.0d0
      fctr=fctr*z
      go to 304
301   tnum=(((((p(7)*z+p(6))*z+p(5))*z+p(4))*z+p(3))*z+p(2))*z+p(1)
      tden=(((((q(7)*z+q(6))*z+q(5))*z+q(4))*z+q(3))*z+q(2))*z+q(1)
      gam=tnum/tden*fctr
      return
      end
      subroutine gammi(g,a,x,gm)
      implicit real*8(a-h,o-z)
      save
c incomplete gamma function=g(a,x)  gm=g(a,0)
      save xcut,ns,nf
      data xcut,ns,nf/1.2d0,15,20/
      y=abs(x)
      xa=y**a
      gm=gam(a)
      if(y.gt.xcut) go to 2
      gami=gm-xa/a
      fn=1.d0 
      term=-xa
      do 3 n=1,ns
      term=-term*y/n
3     gami=gami+term/(a+n)
       g=gami
      return
2     term=2*y
      m=nf  
      do 4 n=1,nf
      term=y+(m-a)*term/(term+m)
4     m=m-1
      gami=exp(-y)*xa/term
       g=gami
      return
      end
      subroutine plint3(r,n,c,i0,pn,np)
      implicit real*8(a-h,o-z)
      save
c integrates sin(r*x)*x**i for i=i0 to n-1+i0 and x from0 to 1
c pn=solid angle /volume, multiply by (r-1)**np
c result goes into c
      dimension c(9)
      complex*16 ti,et,em
 
      if(abs(r).gt.1.d-10) then
      ri=1.d0/r
      ti=cmplx(0.d0,-ri)
      et=cmplx(sin(r)*ri,-cos(r)*ri)
      em=ti*(et-ti)
      do 1 i=1,n+i0+np
      if(i.gt.i0)c(i-i0)=dreal(em)
1     em=ti*(et-i*em)
 
      else
       do 11 i=1,n+np
11     c(i)=1.d0/(i+1+i0)
      endif
      call pmult(c,n,np,pn)
      return
      end
      subroutine plint2(r,n,c,i0,pn,np)
      implicit real*8(a-h,o-z)
      save
      parameter (mbf=999)
      dimension c(9),bess(mbf+2)
c integrates besj0(r*x)*x**i for i from i0-1 to n+i0-2
c using formula 11.1.1 from Abramowitz and Stegun
 
      if(r.ge.1.d-5) then
c determine mbf bessel functions
       nbf=max0(100,6*int(r))
       if(nbf.gt.mbf) then 
         write (*,*) ' danger not enough space in plint2'
         nbf=mbf
       endif
       call mmbsjn(r,nbf,bess,ier)
c      if(ier.ne.0) write (*,*) 'problem in plint2 ',ier
       fact=1.d-10
       con=2.d0/(exp(1.d0)*r)
      endif
 
      do 1 i=1,n+np
      nn=i+i0-2
 
      if(r.lt.1.d-5) then
      c(i)=1.d0
      else
 
      c(i)=bess(2)
      rat=1.d0
      do 2 k=4,nbf,2
      rat=rat*dble(-nn+k-4)/dble(nn+k)
      term=rat*bess(k)*(k-1)
      c(i)=c(i)+term
c stopping criterea: next term will be small relative to total.
2     if(abs(rat).lt.(con*k)**k*fact) go to 5
c     write (*,*) ' loop does not terminate in plint ',r,i,nbf
5     c(i)=2.d0*c(i)/r
      endif
 
1      c(i)=c(i)/(i+i0)
      call pmult(c,n,np,pn)
      return
      end
      subroutine prodi(c,pn,ndim,n,i0,np,i)
      implicit real*8(a-h,o-z)
      save
      dimension c(9)
      do 1 j=1,n+2*np
1     c(j)=1.d0/dble(j+2*(i0-2)+i+ndim)
      call pmult(c,n,2*np,pn)
      return
      end 
      SUBROUTINE SPOCO(A,LDA,N,RCOND,Z,INFO)
      INTEGER LDA,N,INFO
      REAL*8 A(LDA,1),Z(1)
      REAL*8 RCOND
      REAL*8 dDOT,EK,T,WK,WKM
      REAL*8 ANORM,S,dASUM,SM,YNORM
      INTEGER I,J,JM1,K,KB,KP1
      DO 30 J = 1, N
         Z(J) = dASUM(J,A(1,J),1)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            Z(I) = Z(I) + ABS(A(I,J))
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
      ANORM = 0.0d0
      DO 40 J = 1, N
         ANORM = MAX(ANORM,Z(J))
   40 CONTINUE
      CALL SPOFA(A,LDA,N,INFO)
      IF (INFO .NE. 0) GO TO 180
         EK = 1.0d0
         DO 50 J = 1, N
            Z(J) = 0.0d0
   50    CONTINUE
         DO 110 K = 1, N
            IF (Z(K) .NE. 0.0d0) EK = SIGN(EK,-Z(K))
            IF (ABS(EK-Z(K)) .LE. A(K,K)) GO TO 60
               S = A(K,K)/ABS(EK-Z(K))
               CALL dSCAL(N,S,Z,1)
               EK = S*EK
   60       CONTINUE
            WK = EK - Z(K)
            WKM = -EK - Z(K)
            S = ABS(WK)
            SM = ABS(WKM)
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
            KP1 = K + 1
            IF (KP1 .GT. N) GO TO 100
               DO 70 J = KP1, N
                  SM = SM + ABS(Z(J)+WKM*A(K,J))
                  Z(J) = Z(J) + WK*A(K,J)
                  S = S + ABS(Z(J))
   70          CONTINUE
               IF (S .GE. SM) GO TO 90
                  T = WKM - WK
                  WK = WKM
                  DO 80 J = KP1, N
                     Z(J) = Z(J) + T*A(K,J)
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
            Z(K) = WK
  110    CONTINUE
         S = 1.0d0/dASUM(N,Z,1)
         CALL dSCAL(N,S,Z,1)
         DO 130 KB = 1, N
            K = N + 1 - KB
            IF (ABS(Z(K)) .LE. A(K,K)) GO TO 120
               S = A(K,K)/ABS(Z(K))
 
  120       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL dAXPY(K-1,T,A(1,K),1,Z(1),1)
  130    CONTINUE
         S = 1.0d0/dASUM(N,Z,1)
         CALL dSCAL(N,S,Z,1)
         YNORM = 1.0d0
         DO 150 K = 1, N
            Z(K) = Z(K) - dDOT(K-1,A(1,K),1,Z(1),1)
            IF (ABS(Z(K)) .LE. A(K,K)) GO TO 140
               S = A(K,K)/ABS(Z(K))
               CALL dSCAL(N,S,Z,1)
               YNORM = S*YNORM
  140       CONTINUE
            Z(K) = Z(K)/A(K,K)
  150    CONTINUE
         S = 1.0d0/dASUM(N,Z,1)
         CALL dSCAL(N,S,Z,1)
         YNORM = S*YNORM
         DO 170 KB = 1, N
            K = N + 1 - KB
            IF (ABS(Z(K)) .LE. A(K,K)) GO TO 160
               S = A(K,K)/ABS(Z(K))
               CALL dSCAL(N,S,Z,1)
               YNORM = S*YNORM
  160       CONTINUE
            Z(K) = Z(K)/A(K,K)
            T = -Z(K)
            CALL dAXPY(K-1,T,A(1,K),1,Z(1),1)
  170    CONTINUE
         S = 1.0d0/dASUM(N,Z,1)
         CALL dSCAL(N,S,Z,1)
         YNORM = S*YNORM
         IF (ANORM .NE. 0.0d0) RCOND = YNORM/ANORM
         IF (ANORM .EQ. 0.0d0) RCOND = 0.0d0
  180 CONTINUE
      RETURN
      END
      SUBROUTINE SPOSL(A,LDA,N,B)
      INTEGER LDA,N
      REAL*8 A(LDA,1),B(1)
      REAL*8 dDOT,T
      INTEGER K,KB
      DO 10 K = 1, N
         T = dDOT(K-1,A(1,K),1,B(1),1)
         B(K) = (B(K) - T)/A(K,K)
   10 CONTINUE
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL dAXPY(K-1,T,A(1,K),1,B(1),1)
   20 CONTINUE
      RETURN
      END
      subroutine polyd(m,a,ap)
      implicit real*8(a-h,o-z)
      save
c takes the derivative of a polynomial
      dimension a(m),ap(m)
      do 1 i=2,m
1     ap(i-1)=(i-1.d0)*a(i)
      ap(m)=0.d0
      return
      end
      function stfun(r,a,mpoly)
      implicit real*8(a-h,o-z)
      save
c evaluates a polynomial
      dimension a(mpoly)
      stfun=0.d0
      do 20 i=1,mpoly
20    stfun=a(mpoly-i+1)+r*stfun
      return
      end
        subroutine writetb(x,n1,m1,n2,n3,locat,up,ifopcl)
      implicit real*8(a-h,o-z)
      save
c due to peculiarities of UNICOS filename must be passed thru copenfn
c writes table file, the addresses of 3rd index are given in locat
c if abs(ifopcl)=1 open the file
c if ifopcl>0 close file
c 
         dimension x(m1,n2),locat(n3)
      character file*14
      common/copenfn/file
      ln=index(file,' ')-1
      if(iabs(ifopcl).eq.1)then 
      open(21,file=file(1:ln),status='unknown',form='formatted')
      rewind(21)
      endif
      write (21,*) 'RANK 3 ',n1,' ',n2,' ',n3
c     write (21,*) 'GRID ',1,' LINEAR ',0.d0,up
      write (21,'(''GRID 1 LINEAR 0.d0 '',e18.11)')up
      write (21,*) 'BEGIN  rspace table'
      do 10 i3=1,n3
      if(n3.gt.0)then
         l0=locat(i3)
      else
          l0=0
      endif
      do 11 i2=1,n2
      do 12 i1=1,n1
12    write (21,*) x(l0+i1,i2)
11    continue
10    continue
      if(ifopcl.gt.0)close(21)
      return
      end
      SUBROUTINE pTQLI(D,E,N,NP,Z)
      implicit real*8(a-h,o-z)
      DIMENSION D(NP),E(NP),Z(NP,NP)
      DO 11 I=2,N
        E(I-1)=E(I)
11    CONTINUE
      E(N)=0.d0
      DO 15 L=1,N
        ITER=0
1       DO 12 M=L,N-1
          DD=ABS(D(M))+ABS(D(M+1))
          IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12      CONTINUE
        M=N
2       IF(M.NE.L)THEN
          IF(ITER.EQ.30)PAUSE 'too many iterations'
          ITER=ITER+1
          G=(D(L+1)-D(L))/(2.d0*E(L))
          R=SQRT(G**2+1.d0)
          G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
          S=1.d0
          C=1.d0
          P=0.d0
          DO 14 I=M-1,L,-1
            F=S*E(I)
            B=C*E(I)
            IF(ABS(F).GE.ABS(G))THEN
              C=G/F
              R=SQRT(C**2+1.d0)
              E(I+1)=F*R
              S=1.d0/R
              C=C*S
            ELSE
              S=F/G
              R=SQRT(S**2+1.d0)
              E(I+1)=G*R
              C=1.d0/R
              S=S*C
            ENDIF
            G=D(I+1)-P
            R=(D(I)-G)*S+2.d0*C*B
            P=S*R
            D(I+1)=G+P
            G=C*R-B
C     Omit lines from here ...
            DO 13 K=1,N
              F=Z(K,I+1)
              Z(K,I+1)=S*Z(K,I)+C*F
              Z(K,I)=C*Z(K,I)-S*F
13          CONTINUE
C     ... to here when finding only eigenvalues.
14        CONTINUE
          D(L)=D(L)-P
          E(L)=G
          E(M)=0.d0
          GO TO 1
        ENDIF
15    CONTINUE
      RETURN
      END
      SUBROUTINE pINDEXX(N,ARRIN,INDX)
      implicit real*8(a-h,o-z)
c Indexes an array ARRIN of length N, i.e. outputs the array INDX 
c such that ARRIN(INDX(J)) is in ascending order for J=1,2,...,N. 
c The input quantities N and ARRIN are not changed.
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
         INDX(J)=J
11    CONTINUE
      IF(N.EQ.1)RETURN
      L=N/2+1
      IR=N
10    CONTINUE
         IF(L.GT.1)THEN
            L=L-1
            INDXT=INDX(L)
            Q=ARRIN(INDXT)
         ELSE
            INDXT=INDX(IR)
            Q=ARRIN(INDXT)
            INDX(IR)=INDX(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
               INDX(1)=INDXT
               RETURN
            ENDIF
         ENDIF
         I=L
         J=L+L
20       IF(J.LE.IR)THEN
            IF(J.LT.IR)THEN
               IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
            ENDIF
            IF(Q.LT.ARRIN(INDX(J)))THEN
               INDX(I)=INDX(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
            GO TO 20
         ENDIF
         INDX(I)=INDXT
      GO TO 10
      END
      subroutine pmult(c,n,np,pn)
      implicit real*8(a-h,o-z)
      save
      dimension c(n+np)
c now multiply by (r-1)**np
      do 22 k=1,np
      do 22 i=1,n+np-k
22    c(i)=c(i+1)-c(i)
c multiply by pn
      do 24 i=1,n
24    c(i)=pn*c(i)
      return
      end
      subroutine mmbsjn (arg,n,b,ier)
      implicit real*8(a-h,o-z)
      save
c   purpose             - bessel function of the first kind of
c                           nonnegative integer order for
c                           real arguments
c
c   arguments    arg    - input argument. the absolute value of arg must
c                           be less than or equal to 100000. arg must be
c                n      - input parameter specifying the number of
c                           function values to be computed.
c                b      - output vector of length n containing the
c                           computed function values. b must be typed
c                           appropriately in the calling program.
c                           b(1) will contain the computed value for
c                           order zero, b(2) will contain the computed
c                           value for order 1, b(3) for order 2, etc.
c                ier    - error parameter. (output)
c                           ier = 129 + j indicates that b(i), (i=1,j)
c                             are computed to machine precision, but
c                             precision is lost for b(i), (i=j+1,n.)
c                             see the programming notes.
c
c     dimension               b(1)
      dimension b(n)
c                                  first executable statement
      ier = 0
      tempa = abs(arg)
      magx = int(tempa)
      if(n.gt.0 .and. magx.le.100000) go to 10
c                                  error return -- arg,n is out of range
      ier = 129
   10 rsign = 1.d0
      ncalc = n
c                                  use 2-term ascending series for
c                                    small arg
      tmpa4 = tempa**4.d0
      smallx = 1.d-14
      if(tmpa4.ge.smallx) go to 20
c                                  two-term ascending series for
c                                    small arg
      tempa = 1.d0
      tempb = -.25d0*arg*arg*rsign
      b(1) = 1.d0+tempb
      if(n.eq.1) go to 9005
      do 15 nn=2,n
         tempa = tempa*arg/(dble(2*nn-2))
         b(nn) = tempa*(1.d0+tempb/(dble(nn)))
   15 continue
      go to 9005
c                                  initialize the calculation of p*s
   20 nbmx = n-magx
      nn = magx+1
      plast = 1.d0
      p = (dble(2*nn))/tempa
c                                  calculate general significance test
      test = 2.d14
      m = 0
      if(nbmx.lt.3) go to 30
c                                  calculate p*s until nn=n-1.
c                                    check for possible overflow.
      tover = 1.d35
      nstart = magx+2
      nend = n-1
      do 25 nn=nstart,nend
         pold = plast
         plast = p
         p = (dble(2*nn))*plast/tempa-rsign*pold
         if(p-tover) 25, 25, 35
   25 continue
      nn = nend
c                                  calculate special significance test
c                                    for nbmx.gt.2.
c
      test = max(test,sqrt(plast*1.d14)*sqrt(2.d0*p))
c
c                                  calculate p*s until significance
c                                    test passes
   30 nn = nn+1
      pold = plast
      plast = p
      p = (dble(2*nn))*plast/tempa-rsign*pold
      if(p.lt.test) go to 30
      if(m.eq.1) go to 55
c                                  for j*s, a strong variant of the test
c                                    is necessary. calculate it, and
c                                    calculate p*s until this test is
c                                    passed.
      m = 1
      tempb = p/plast
      tempc = (dble(nn+1))/tempa
      if(tempb+1.d0/tempb.gt.2.d0*tempc)tempb=tempc+sqrt(tempc**2-1.d0)
      test = test/sqrt(tempb-1.d0/tempb)
      if(p-test) 30, 55, 55
c                                  to avoid overflow, divide p*s by
c                                    tover.  calculate p*s until
c                                    abs(p).gt.1.
   35 tover = 1.d35
      p = p/tover
      plast = plast/tover
      psave = p
      psavel = plast
      nstart = nn+1
   40 nn = nn+1
      pold = plast
      plast = p
      p = (dble(2*nn))*plast/tempa-rsign*pold
      if(p.le.1.d0) go to 40
      tempb = (dble(2*nn))/tempa
      tempc = .5d0*tempb
      tempb = plast/pold
      if(tempb+1.d0/tempb.gt.2.d0*tempc)tempb=tempc+sqrt(tempc**2-1.d0)
c
c                                  calculate backward test, and find
c                                    ncalc, the highest nn such that the
c                                    test is passed.
      test = .5d0*pold*plast*(1.d0-1.d0/tempb**2)*1.d-14
      p = plast*tover
      nn = nn-1
      nend = min0(n,nn)
      do 45 ncalc=nstart,nend
         pold = psavel
         psavel = psave
         psave = (dble(2*nn))*psavel/tempa-rsign*pold
         if(psave*psavel-test) 45, 45, 50
   45 continue
      ncalc = nend+1
   50 ncalc = ncalc-1
c                                  the sum b(1)+2b(3)+2b(5)... is used
c                                    to normalize. m, the coefficient of
c                                    b(nn), is initialized to 2 or 0.
   55 nn = nn+1
      m = 2*nn-4*(nn/2)
c                                  initialize the backward recursion and
c                                    the normalization sum
      tempb = 0.d0
      tempa = 1.d0/p
      sum = (dble(m))*tempa
      nend = nn-n
      if(nend) 80, 70, 60
c                                  recur backward via difference
c                                    equation, calculating (but not
c                                    storing) b(nn), until nn=n.
   60 do 65 l=1,nend
         nn = nn-1
         tempc = tempb
         tempb = tempa
         tempa = ((dble(2*nn))*tempb)/arg-rsign*tempc
         m = 2-m
         sum = sum+(dble(m))*tempa
   65 continue
c                                  store b(nn)
   70 b(nn) = tempa
      if(n.gt.1) go to 75
c                                  n=1.  since 2*tempa is added to the
c                                    sum, tempa must be subtracted
      sum = sum-tempa
      go to 110
c                                  calculate and store b(nn-1)
   75 nn = nn-1
      b(nn) = ((dble(2*nn))*tempa)/arg-rsign*tempb
      if(nn.eq.1) go to 105
      m = 2-m
      sum = sum+(dble(m))*b(nn)
      go to 90
c                                  nn.lt.n, so store b(nn) and set
c                                  higher orders to zero
   80 b(nn) = tempa
      nend = -nend
      do 85 l=1,nend
         itemp = nn+l
         b(itemp) = 0.0d0
   85 continue
   90 nend = nn-2
      if(nend.eq.0) go to 100
c                                  calculate via difference equation and
c                                    store b(nn), until nn=2
      do 95 l=1,nend
         nn = nn-1
         b(nn) = ((dble(2*nn))*b(nn+1))/arg-rsign*b(nn+2)
         m = 2-m
         sum = sum+(dble(m))*b(nn)
   95 continue
c                                  calculate b(1)
  100 b(1) = 2.d0*b(2)/arg-rsign*b(3)
  105 sum = sum+b(1)
c                                  normalize--if ize=1, divide sum by
c                                    cosh(arg). divide all b(nn) by sum.
  110 continue
      do 115 nn=1,n
  115 b(nn) = b(nn)/sum
      if(ncalc.eq.n) go to 9005
      ier = 129+ncalc
 9005 return
      end
 
      FUNCTION dASUM(N,SX,INCX)
      REAL*8 SX(1),dASUM
      dASUM = 0.0d0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
      NS = N*INCX
          DO 10 I=1,NS,INCX
          dASUM = dASUM + ABS(SX(I))
   10     CONTINUE
      RETURN
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        dASUM = dASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        dASUM = dASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))
     1  + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END
      SUBROUTINE SPOFA(A,LDA,N,INFO)
      INTEGER LDA,N,INFO
      REAL*8 A(LDA,1)
      REAL*8 dDOT,T
      REAL*8 S
      INTEGER J,JM1,K
         DO 30 J = 1, N
            INFO = J
            S = 0.0d0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - dDOT(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            S = A(J,J) - S
            IF (S .LE. 0.0d0) GO TO 40
            A(J,J) = SQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
 
      SUBROUTINE pAXPY(N,SA,SX,INCX,SY,INCY)
      REAL*8 SX(1),SY(1),SA
      IF(N.LE.0.OR.SA.EQ.0.d0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END
 
      FUNCTION pDOT(N,SX,INCX,SY,INCY)
      REAL*8 SX(1),SY(1),pDOT
      pDOT = 0.0d0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        pDOT = pDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        pDOT = pDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        pDOT = pDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     1   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
      RETURN
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        pDOT = pDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END
      function cvmgp(x1,x2,x3)
      real*8 x1,x2,x3,cvmgp
      if(x3.ge.0.d0)cvmgp=x1
      if(x3.lt.0.d0)cvmgp=x2
      return
      end

      subroutine traduci(fname,np,nt,drt,v0,ifv0,filename,nk
     &                  ,vsr,vlr,self,const,nkappa)
      implicit none

      integer mnt,mnk,mnts
      parameter (mnt=900,mnk=9000,mnts=501)

      integer i,j,k,nt,ivexp,np,ifv0,nk,nkappa
      real*8 drt,vt(0:mnt,4),vkt(mnk),v0,aux
      real*8 vsr(mnts,2),vlr(mnk),self,const
      character*30 fname
      character*80 record
      character*30 filename
      write(*,*)'TRADUCI: np nt drt ',np,nt,drt
      write(*,*)'fname = ',fname

c parte in spazio r
      do i=0,nt
       j=i+1
       vt(i,1)=vsr(j,1)
      enddo
      do i=0,nt
       j=i+1
       vt(i,2)=vsr(j,2)
      enddo
      call pspline(mnt,nt,drt,vt)
      do i=1,4
       vt(nt,i)=0.d0
      enddo

      i=index(filename,' ')-1
      open(11,file=filename(1:i))
      write(11,*)nt,drt
      do i=0,nt
       write(11,*)(vt(i,j),j=1,4)
      enddo

      write(11,*)'fitpn'
      write(11,*)'0'

      write(11,'(a)')'kspace'
c parte in spazio k
      aux=self
      v0=const
      v0=(v0*np+aux)*0.5d0

      do i=1,nkappa
       vkt(i)=vlr(i)
      enddo

      write(11,*)i-1
      nk=i-1
      do j=1,i-1
       write(11,*)vkt(j)
      enddo
      if(ifv0.eq.0)then
       write(11,*)v0
       v0=0.d0
      else
       write(11,*)'0.d0'
      endif
      close(11)

      return
      end

      subroutine pspline(mnt,nt,drt,t)
c scale and spline
      implicit none
      integer mnt,nt,i
      real*8 drt,t(0:mnt,4)
      call dscal(nt+1,drt,t(0,2),1)
      call dscal(nt+1,drt*drt*0.5,t(0,3),1)
      do i=0,nt-1
       t(i,3)= 3.d0*(t(i+1,1)-t(i,1))  -(t(i+1,2)+2.d0*t(i,2))
       t(i,4)=-2.d0*(t(i+1,1)-t(i,1))  +(t(i+1,2)+     t(i,2))
      enddo
      return
      end
 
      subroutine krlv(cut,a,ndim,mdim,gvect,mnk,ng,verbose)
c vettori del reticolo reciproco
      implicit none
      integer idim,jdim,ndim,mdim,ix(3),nx(3),nkspan(3),i,j
     &       ,ng,mnk,npts,ipts,nzero,verbose
      real*8 a(mdim,ndim),arlv(3,3),vol,pi,tpiba
     &      ,g2,gvect(mdim,mnk),g(3),cut,cut1,cut2,k2,small,phi
     &      ,norm(3)

      if(verbose.ne.0)write(6,*)'=========>> krlv <<=========='
      if(ndim.lt.1.or.ndim.gt.3)then
       write(6,*)'    ndim = ',ndim,' ...stop.'
       stop
      endif
      pi=acos(-1.d0)
      small=1.d-7
c genera i vettori primitivi del rlv
      if(ndim.eq.3)then
       vol=a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3))
     &    -a(2,1)*(a(1,2)*a(3,3)-a(3,2)*a(1,3))
     &    +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
       arlv(1,1)= (a(2,2)*a(3,3)-a(3,2)*a(2,3))*2*pi/vol
       arlv(2,1)=-(a(1,2)*a(3,3)-a(3,2)*a(1,3))*2*pi/vol
       arlv(3,1)= (a(1,2)*a(2,3)-a(2,2)*a(1,3))*2*pi/vol
       arlv(1,2)= (a(2,3)*a(3,1)-a(3,3)*a(2,1))*2*pi/vol
       arlv(2,2)=-(a(1,3)*a(3,1)-a(3,3)*a(1,1))*2*pi/vol
       arlv(3,2)= (a(1,3)*a(2,1)-a(2,3)*a(1,1))*2*pi/vol
       arlv(1,3)= (a(2,1)*a(3,2)-a(3,1)*a(2,2))*2*pi/vol
       arlv(2,3)=-(a(1,1)*a(3,2)-a(3,1)*a(1,2))*2*pi/vol
       arlv(3,3)= (a(1,1)*a(2,2)-a(2,1)*a(1,2))*2*pi/vol
      elseif(ndim.eq.2)then
       vol=a(1,1)*a(2,2)-a(2,1)*a(1,2)
       arlv(1,1)= a(2,2)*2*pi/vol
       arlv(2,1)=-a(1,2)*2*pi/vol
       arlv(1,2)=-a(2,1)*2*pi/vol
       arlv(2,2)= a(1,1)*2*pi/vol
      elseif(ndim.eq.1)then
       vol=a(1,1)
       arlv(1,1)= 2*pi/vol
      endif
      if(verbose.ne.0)write(6,*)'    primitive vectors:'
      do idim=1,ndim
       if(verbose.ne.0)write(6,'(4x,3e20.5)')(arlv(i,idim),i=1,ndim)
      enddo
c smallest primitive vector
      norm(1)=0.d0
      do j=1,ndim
       norm(1)=norm(1)+arlv(j,1)**2
      enddo
      norm(1)=sqrt(norm(1))
      tpiba=norm(1)
      do i=2,ndim
       norm(i)=0.d0
       do j=1,ndim
        norm(i)=norm(i)+arlv(j,i)**2
       enddo
       norm(i)=sqrt(norm(i))
       tpiba=min(tpiba,norm(i))
      enddo
      if(verbose.ne.0)write(6,*)'    norm of shortest rlv: ',tpiba
c genera i vettori del rlv di modulo meno di cut
      cut1=cut*tpiba
      cut2=cut1**2
      if(verbose.ne.0)write(6,*)'    cut: ',cut1
c zona dove cercare
      if(ndim.eq.1)then
       nkspan(1)=int(abs((cut1+small)/arlv(1,1)))
      elseif(ndim.eq.2)then
       norm(1)=sqrt(arlv(1,1)**2+arlv(2,1)**2)
       norm(2)=sqrt(arlv(1,2)**2+arlv(2,2)**2)
       phi=acos((arlv(1,1)*arlv(1,2)+arlv(2,1)*arlv(2,2))
     &                              /(norm(1)*norm(2)))
       nkspan(1)=int(abs((cut1+small)/(norm(1)*sin(phi))))
       nkspan(2)=int(abs((cut1+small)/(norm(2)*sin(phi))))
      elseif(ndim.eq.3)then
       norm(1)=sqrt(arlv(1,1)**2+arlv(2,1)**2+arlv(3,1)**2)
       norm(2)=sqrt(arlv(1,2)**2+arlv(2,2)**2+arlv(3,2)**2)
       norm(3)=sqrt(arlv(1,3)**2+arlv(2,3)**2+arlv(3,3)**2)
       nkspan(1)=  int(abs((cut1+small)/
     &       ((arlv(1,1)*( arlv(2,2)*arlv(3,3)-arlv(3,2)*arlv(2,3))
     &        +arlv(2,1)*(-arlv(1,2)*arlv(3,3)+arlv(3,2)*arlv(1,3))
     &        +arlv(3,1)*( arlv(1,2)*arlv(2,3)-arlv(2,2)*arlv(1,3)))
     &                              /(norm(2)*norm(3)))))
       nkspan(2)=  int(abs((cut1+small)/
     &       ((arlv(1,2)*( arlv(2,3)*arlv(3,1)-arlv(3,3)*arlv(2,1))
     &        +arlv(2,2)*(-arlv(1,3)*arlv(3,1)+arlv(3,3)*arlv(1,1))
     &        +arlv(3,2)*( arlv(1,3)*arlv(2,1)-arlv(2,3)*arlv(1,1)))
     &                              /(norm(3)*norm(1)))))
       nkspan(3)=  int(abs((cut1+small)/
     &       ((arlv(1,3)*( arlv(2,1)*arlv(3,2)-arlv(3,1)*arlv(2,2))
     &        +arlv(2,3)*(-arlv(1,1)*arlv(3,2)+arlv(3,1)*arlv(1,2))
     &        +arlv(3,3)*( arlv(1,1)*arlv(2,2)-arlv(2,1)*arlv(1,2)))
     &                              /(norm(1)*norm(2)))))
      endif
c quanti punti ci sono
      npts=1
      do i=1,ndim
       npts=npts*(2*nkspan(i)+1)
      enddo
      if(verbose.ne.0)write(6,*)'    nkspan '
     &               ,(nkspan(i),i=1,ndim),' npts ',npts
c quanti punti ci sono tra tutti gli idim piu' piccoli
      do idim=1,ndim
       nx(idim)=1
       do jdim=2,idim
        nx(idim)=nx(idim)*(2*nkspan(jdim-1)+1)
       enddo
      enddo
      ng=0
c loop sui punti
      do ipts=1,npts
c ricostruisce gli indici
       i=ipts-1
       do idim=ndim,2,-1
        ix(idim)=i/nx(idim)-nkspan(idim)
        i=mod(i,nx(idim))
       enddo
       ix(1)=i-nkspan(1)
c nzero serve a togliere g=0 e tenere solo la meta' degli alrti vettori
       nzero=0
       g2=0.d0
       call r_set(ndim,g,0.d0)
       do idim=1,ndim
        if(nzero.eq.0)nzero=ix(idim)
        do jdim=1,ndim
         g(jdim)=g(jdim)+ix(idim)*arlv(jdim,idim)
        enddo
       enddo
       do jdim=1,ndim
        g2=g2+g(jdim)*g(jdim)
       enddo
       if(g2.gt.cut2)go to 1 
       if(nzero.le.0)go to 1
c trovato un vettore
       ng=ng+1
       if(ng.gt.mnk)stop'krlv: nrlv.gt.mnk'
c ordinamento
       do idim=1,ndim
        gvect(idim,ng)=g(idim)
       enddo    
       do i=1,ng-1
        k2=0.d0
        do idim=1,ndim
         k2=k2+gvect(idim,i)**2
        enddo
        if(g2.lt.k2)then
         do j=ng,i+1,-1
          do idim=1,ndim
           gvect(idim,j)=gvect(idim,j-1)
          enddo
         enddo
         do idim=1,ndim
          gvect(idim,i)=g(idim)
         enddo
         go to 1
        endif
       enddo
    1 enddo
      if(verbose.ne.0)write(6,*)'    nrlv = ',ng
      if(verbose.ne.0)write(6,*)'-------> end krlv <--------'
      return
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

      subroutine eta_yk(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
c backflow pseudopotential from kwon et al. 1998
      implicit none
      integer i,i0,nt,j,jp,m_parm,stdin_flag,spline_flag
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),c,xi,g,dg,ddg,beta0,p(m_parm)
      common /c_obsmooth/c,xi,g,dg,ddg,beta0
      common /c_wrkparm/p,jp
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'yk backflow lambda_B (1+s_Br)/(r_B+w_Br+r^(7/2))'
       write(6,*)'from kwon et al 98'
       write(6,*)'enter lambda_B, s_B, r_B, w_B'
       read(*,*)p(1),p(2),p(3),p(4)
       jp=4
      endif
      p(3)=sqrt(p(3)**2)
c     c=drt*nt*(4.d0/5.d0)
c     xi=c/4.d0
      c=drt*(nt*5/6)
      xi=drt*nt-c
      g=p(1)*(1+p(2)*c)/(p(3)+p(4)*c+c**3.5)
      dg=(p(2)*g-(p(4)+3.5*c**2.5)*g**2/p(1))/(1+p(2)*c)
      ddg=-(12*c**2*g**2+2*g*dg*(p(4)+3.5*c**2.5))/(1+p(2)*c)/p(1)
      do i=i0,nt
       r=r0+i*drt
       if(r.lt.c)then
        f(i)=p(1)*(1+p(2)*r)/(p(3)+p(4)*r+r**(3.5))
        df(i)=(p(2)*f(i)-(p(4)+3.5*r**2.5)*f(i)**2/p(1))/(1+p(2)*r)
       else
        go to 1
       endif
      enddo
    1 call obsmooth(r0,drt,i,nt,f,df)
      do j=0,i-1
       f(j)=f(j)-g+beta0
      enddo
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
      subroutine u2_ob_cos(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
      implicit none
      integer i,j,k,n,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm),x,xi,g,dg,ddg,beta0
      real*8 pl,c,dc,ddc,kk,p3,p4
      common /c_wrkparm/p,jp
      common /c_obsmooth/x,xi,g,dg,ddg,beta0
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'u2_ob_cos pseudo: [(ar+br2)/(1+c2r+d2r2) + e]'
       write(6,*)'                 *[1+sum_i e_i cos(k_i*r)]'
       write(6,*)'enter cusp b c d e'
       read(*,*)p(1),p(2),p(3),p(4),p(5)
       write(6,*)'enter # cos, e_i'
       read(*,*)n,(p(5+i),i=1,n)
       jp=5+n
      endif
      n=jp-5
      p3=p(3)**2
      p4=p(4)**2
c     x=drt*nt*(5.d0/6.d0)
c     xi=x/5.d0
      x=drt*(nt*5/6)
      xi=drt*nt-x
      pl=acos(-1.d0)/(2.d0*drt*nt)
      c=1
      dc=0
      ddc=0
      do k=1,n
       kk=2*k*pl
       c=c+p(5+k)*(1-cos(kk*x))
       dc=dc+p(5+k)*sin(kk*x)*kk
       ddc=ddc+p(5+k)*cos(kk*x)*kk**2
      enddo
      g=(p(1)*x+p(2)*x**2)                   /(1+p3*x+p4*x**2) + p(5)
      dg=(p(1)+2*p(2)*x)                     /(1+p3*x+p4*x**2)
     &  -(p(1)*x+p(2)*x**2)*(p3+2*p4*x)      /(1+p3*x+p4*x**2)**2
      ddg=2*p(2)                             /(1+p3*x+p4*x**2)
     &   -2*(p(1)+2*p(2)*x)*(p3+2*p4*x)      /(1+p3*x+p4*x**2)**2
     &   -(p(1)*x+p(2)*x**2)*2*p4            /(1+p3*x+p4*x**2)**2
     &   +2*(p(1)*x+p(2)*x**2)*(p3+2*p4*x)**2
     &                                       /(1+p3*x+p4*x**2)**3
      ddg=ddg+ddc+2*dg*dc
      dg=g*dc+c*dg
      g=g*c
      do i=i0,nt
       r=r0+i*drt
       r=max(r,1.d-10)
       if(r.lt.x)then
        c=1.d0
        dc=0.d0
        do k=1,n
        kk=2*k*pl
       c=c+p(5+k)*(1-cos(kk*r))
       dc=dc+p(5+k)*sin(kk*r)*kk
c      ddc=ddc+p(5+k)*cos(kk*r))*kk**2

        enddo

       f(i)=(p(1)*r+p(2)*r**2)              /(1+p3*r+p4*r**2) + p(5)
        df(i)=(p(1)+2*p(2)*r)               /(1+p3*r+p4*r**2)
     &       -(p(1)*r+p(2)*r**2)*(p3+2*p4*r)/(1+p3*r+p4*r**2)**2
        df(i)=df(i)*c+f(i)*dc
        f(i)=f(i)*c
       else
        goto 1
       endif
      enddo
    1 call obsmooth(r0,drt,i,nt,f,df)
      do j=0,i-1
       f(j)=f(j)-g+beta0
      enddo
      return
      end

      subroutine obsmooth(r0,drt,i0,nt,f,df)
      implicit none
      integer i,i0,nt
      real*8 r0,drt,f(0:nt),df(0:nt),c,xi,g,dg,ddg
     &      ,r,beta0,beta1,beta2,beta3,beta4
      common /c_obsmooth/c,xi,g,dg,ddg,beta0
      beta0=-0.5d0*dg*xi-ddg*xi**2/12.d0
      beta1=dg
      beta2=0.5d0*ddg
      beta3=-dg/xi**2-2.d0/3.d0*ddg/xi
      beta4=0.5d0*dg/xi**3+0.25d0*ddg/xi**2
      do i=i0,nt
       r=r0+i*drt
       f(i)=beta0
     &     +beta1*(r-c)
     &     +beta2*(r-c)**2
     &     +beta3*(r-c)**3
     &     +beta4*(r-c)**4
       df(i)=     beta1
     &      +2.d0*beta2*(r-c)
     &      +3.d0*beta3*(r-c)**2
     &      +4.d0*beta4*(r-c)**3
      enddo
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

      SUBROUTINE TQLI(D,E,N,NP,Z)
      implicit real*8(a-h,o-z)
      DIMENSION D(NP),E(NP),Z(NP,NP)
      DO 11 I=2,N
        E(I-1)=E(I)
11    CONTINUE
      E(N)=0.d0
      DO 15 L=1,N
        ITER=0
1       DO 12 M=L,N-1
          DD=ABS(D(M))+ABS(D(M+1))
          IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12      CONTINUE
        M=N
2       IF(M.NE.L)THEN
          IF(ITER.EQ.30)PAUSE 'too many iterations'
          ITER=ITER+1
          G=(D(L+1)-D(L))/(2.d0*E(L))
          R=SQRT(G**2+1.d0)
          G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
          S=1.d0
          C=1.d0
          P=0.d0
          DO 14 I=M-1,L,-1
            F=S*E(I)
            B=C*E(I)
            IF(ABS(F).GE.ABS(G))THEN
              C=G/F
              R=SQRT(C**2+1.d0)
              E(I+1)=F*R
              S=1.d0/R
              C=C*S
            ELSE
              S=F/G
              R=SQRT(S**2+1.d0)
              E(I+1)=G*R
              C=1.d0/R
              S=S*C
            ENDIF
            G=D(I+1)-P
            R=(D(I)-G)*S+2.d0*C*B
            P=S*R
            D(I+1)=G+P
            G=C*R-B
C     Omit lines from here ...
            DO 13 K=1,N
              F=Z(K,I+1)
              Z(K,I+1)=S*Z(K,I)+C*F
              Z(K,I)=C*Z(K,I)-S*F
13          CONTINUE
C     ... to here when finding only eigenvalues.
14        CONTINUE
          D(L)=D(L)-P
          E(L)=G
          E(M)=0.d0
          GO TO 1
        ENDIF
15    CONTINUE
      RETURN
      END

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      implicit real*8(a-h,o-z)
c Indexes an array ARRIN of length N, i.e. outputs the array INDX 
c such that ARRIN(INDX(J)) is in ascending order for J=1,2,...,N. 
c The input quantities N and ARRIN are not changed.
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
         INDX(J)=J
11    CONTINUE
      IF(N.EQ.1)RETURN
      L=N/2+1
      IR=N
10    CONTINUE
         IF(L.GT.1)THEN
            L=L-1
            INDXT=INDX(L)
            Q=ARRIN(INDXT)
         ELSE
            INDXT=INDX(IR)
            Q=ARRIN(INDXT)
            INDX(IR)=INDX(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
               INDX(1)=INDXT
               RETURN
            ENDIF
         ENDIF
         I=L
         J=L+L
20       IF(J.LE.IR)THEN
            IF(J.LT.IR)THEN
               IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
            ENDIF
            IF(Q.LT.ARRIN(INDX(J)))THEN
               INDX(I)=INDX(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
            GO TO 20
         ENDIF
         INDX(I)=INDXT
      GO TO 10
      END

      FUNCTION DDOT(N,SX,INCX,SY,INCY)
      REAL*8 SX(n*incx+1),SY(n*incy+1),DDOT
      DDOT = 0.0d0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DDOT = DDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DDOT = DDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     1   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
      RETURN
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        DDOT = DDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END

      SUBROUTINE DSCAL(N,SA,SX,INCX)
      REAL*8 SA,SX(n)
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          SX(I) = SA*SX(I)
   10     CONTINUE
      RETURN
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END

      SUBROUTINE DAXPY(N,SA,SX,INCX,SY,INCY)
      REAL*8 SX(n),SY(n),SA
      IF(N.LE.0.OR.SA.EQ.0.d0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END

      subroutine coulomb(r0,drt,i0,nt,f,df,stdin_flag)
      implicit none
      integer i,i0,nt,stdin_flag
      real*8 r,r0,drt,f(0:nt),df(0:nt),rs
      common /c_rs/rs
      do i=i0,nt
       f(i)=2.d0/rs
       df(i)=0.d0
      enddo
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
         iseed(3)=i3
         iseed(4)=i4
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

      subroutine polycusp(r0,drt,i0,nt,f,df,stdin_flag,spline_flag)
c eta function for backflow correlation
      implicit none
      integer i,i0,nt,stdin_flag,spline_flag,m_parm,jp
      parameter(m_parm=20)
      real*8 r,r0,drt,f(0:nt),df(0:nt),p(m_parm),p0,p1,q0,q1
      real*8 a0,a1,a2,a3,a4,a5,a6,a7,a8,rcut,cusp
      common /c_wrkparm/p,jp
c lopez-rios
      spline_flag=1
      if(stdin_flag.ne.0)then
       write(6,*)'(cut-r)**3 * (sum_i a_i r**i)'
       write(6,*)'enter a_0,a_2,a_3,a_4,a_5,a_6,a_7,a_8,cut,cusp'
       write(6,*)'a_1 from cusp condition'
       read(*,*)p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10)
       jp=10
      endif
      rcut=p(9)
      cusp=p(10)
      a0=p(1)
      a1=cusp+3*a0/rcut
      a2=p(2)
      a3=p(3)
      a4=p(4)
      a5=p(5)
      a6=p(6)
      a7=p(7)
      a8=p(8)
      do i=i0,nt
       r=r0+drt*i
       if(r.lt.rcut)then
        p0=((rcut-r)/rcut)**3
        q0=a0 +a1*r   +a2*r**2+a3*r**3+a4*r**4
     &        +a5*r**5+a6*r**6+a7*r**7+a8*r**8
        p1=-3*((rcut-r)/rcut)**2/rcut
        q1=    a1     +2*a2*r   +3*a3*r**2+4*a4*r**3
     &      +5*a5*r**4+6*a6*r**5+7*a7*r**6+8*a8*r**7
        f(i)=p0*q0
        df(i)=p1*q0+p0*q1
       else
        f(i)=0.d0
        df(i)=0.d0
       endif
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
