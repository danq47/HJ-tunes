c A trial analysis to test out the HJ generator

c Initialise histograms
      subroutine init_hist
      implicit none
      include 'LesHouches.h'
      include 'pwhg_math.h'
      integer 			 ixx,jxx,kxx,len_pt,len_y,t1,t2,ls1
C     integer maxbins1,maxbins2
C     parameter (maxbins1=30,maxbins2=8)
      real*8           pt_bins1(100),y_bins1(100)
      character * 20   s1(100)
      character * 20   tmp,tmp2
      common/bins/pt_bins1,y_bins1,len_pt,len_y
      integer          lenocc
      external         lenocc
      
      call inihists
      
c     Initialise everything correctly  
      do ixx=1,100
         pt_bins1(ixx) = -1d6
         y_bins1(ixx)  = -1d6
         s1(ixx)       = ''
      enddo
      
cccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                c
c           INPUT CUSTOM BINS HERE               c
c                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c 						 c
c     Enter custom bin edges (from the lowest	 c
c     to the highest. Enter the last bin edge 	 c
c     twice so we don't get any funny business	 c
c     from its original value of -1m.		 c
c     Enter bins in the form pt_bins1(x)	 c
c 						 c
cccccccccccccccccccccccccccccccccccccccccccccccccc

      pt_bins1(1)=0d0
      pt_bins1(2)=10d0
      pt_bins1(3)=20d0
      pt_bins1(4)=30d0
      pt_bins1(5)=40d0
      pt_bins1(6)=50d0
      pt_bins1(7)=60d0
      pt_bins1(8)=80d0
      pt_bins1(9)=100d0
      pt_bins1(10)=100d0

      y_bins1(1)=-4d0
      y_bins1(2)=-3d0
      y_bins1(3)=-2d0
      y_bins1(4)=-1.5d0
      y_bins1(5)=-0.75d0
      y_bins1(6)=-0.5d0
      y_bins1(7)=-0.25d0
      y_bins1(8)=0d0
      y_bins1(9)=0.25d0
      y_bins1(10)=0.5d0
      y_bins1(11)=0.75d0
      y_bins1(12)=1d0
      y_bins1(13)=1.5d0
      y_bins1(14)=2d0
      y_bins1(15)=3d0
      y_bins1(16)=4d0
      y_bins1(17)=4d0

      len_pt=0                  ! length of pt and y bins array 
      len_y=0                   ! (ignoring any terms initialised to -1m)
      do jxx=1,100
         if(pt_bins1(jxx).gt.-1d5) len_pt=len_pt + 1
         if(y_bins1(jxx).gt.-1d5) len_y=len_y + 1
      enddo
      
c     total plots
      call bookupeqbins('sigmatot',1d0,-1d0,2d0)
      call bookup('H-pt',len_pt - 1,pt_bins1)    ! nbins=len_x - 1
      call bookup('H-y',len_y - 1,y_bins1)       ! because we entered
      call bookup('ptj1',len_pt - 1,pt_bins1)    ! both the upper and
      call bookup('Yj1',len_y - 1,y_bins1)       ! lower edges
      call bookup('ptj2',len_pt - 1,pt_bins1)    ! both the upper and
      call bookup('Yj2',len_y - 1,y_bins1)       ! lower edges
      
c     pt-j1 and pt-H cut into y_H bins
      
      do kxx=1,len_y
         
         write(tmp,"(F6.2)") y_bins1(kxx) ! Write bin edges as strings
         
         t1=lenocc(tmp)         ! trim the strings down
         t2=lenocc(tmp2)
         
         s1(kxx) = tmp2(1:t1)//'-yH-'//tmp(1:t2) ! plot suffices
         ls1=lenocc(s1(kxx))    ! trim down plot suffix
         
         if(kxx.eq.1) then      ! initial cut - take from Y-higgs=-inft up to the lowest bin edge
            call bookup('H-pt'//'--inf-yH-'//tmp,len_pt - 1,pt_bins1)
            call bookup('ptj1'//'--inf-yH-'//tmp,len_pt - 1,pt_bins1)
!     call bookup('ptj2'//'--inf-yH-'//tmp,len_pt - 1,pt_bins1)
            tmp2=tmp
         elseif(kxx.gt.1.and.kxx.lt.len_y) then ! This includes all the different cuts in the middle
            call bookup('H-pt'//s1(kxx)(1:ls1),len_pt - 1,pt_bins1)
            call bookup('ptj1'//s1(kxx)(1:ls1),len_pt - 1,pt_bins1)
!     call bookup('ptj2'//s1(kxx)(1:ls1),len_pt - 1,pt_bins1)
            tmp2=tmp
         elseif(kxx.eq.len_y) then ! and this cuts Y-higgs > highest bin edge
            call bookup('H-pt'//tmp2(1:t2)//'-yH-inf',len_pt - 1,pt_bins1)
            call bookup('ptj1'//tmp2(1:t2)//'-yH-inf',len_pt - 1,pt_bins1)
!     call bookup('ptj2'//tmp2(1:t2)//'-yH-inf',len_pt - 1,pt_bins1)
         endif
         
      enddo      	
      end
      
c     Analysis subroutine      
      subroutine analysis(dsig0)
      implicit none
      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_math.h'
      include 'pwhg_weights.h'
      include 'pwhg_lhrwgt.h'
      real * 8     dsig0,dsig(7)
      integer      nweights
      logical      iniwgts
      data         iniwgts/.true./
      save         iniwgts
      character*6  WHCPRG
      common/cWHCPRG/WHCPRG
      data         WHCPRG/'NLO   '/
      integer      ihep,ixx,i_higgs,i_j1,i_j2,id,kxx
      integer      mjets,njets,maxjet
      parameter   (maxjet=2048)
      logical      IsForClustering(maxjet)
      real * 8     ktj(maxjet),etaj(maxjet),rapj(maxjet),
     1     phij(maxjet),pj(4,maxjet),jetRadius,ptrel(4),ptmin
      real * 8     ph(4)
      real * 8     ptj1,ptj2,pt_test,pt_higgs,y_higgs
      real * 8     y,eta,pt,m
      real * 8     powheginput
      real * 8 	 pt_bins1(100),y_bins1(100)
      integer      len_pt,len_y,lenocc,t1,t2,ls1
      character*20 s1(100),tmp,tmp2
      common/bins/ pt_bins1,y_bins1,len_pt,len_y
      external     powheginput,lenocc
      
      if (iniwgts) then
         write(*,*) '*********************'
         if(whcprg.eq.'NLO    ') then
            write(*,*) ' NLO ANALYSIS      '
            weights_num=0
         elseif(WHCPRG.eq.'LHE    ') then
            write(*,*) ' LHE ANALYSIS      '
         elseif(WHCPRG.eq.'HERWIG ') then
            write(*,*) ' HERWIG ANALYSIS   '
         elseif(WHCPRG.eq.'PYTHIA ') then
            write(*,*) ' PYTHIA ANALYSIS   '
         elseif(WHCPRG.eq.'PYTHIA8') then
            write(*,*) ' PYTHIA8 ANALYSIS   '
         endif
         write(*,*) '*********************'
         if(weights_num.eq.0) then
            call setupmulti(1)
         else
            call setupmulti(weights_num)
         endif
         iniwgts=.false.
      endif
      
      dsig=0
      if(weights_num.eq.0) then
         dsig(1)=dsig0
      else
         dsig(1:weights_num)=weights_val(1:weights_num)
      endif
      if(sum(abs(dsig)).eq.0) return
      
c     Initialise everything to zero
      i_higgs=0 		! index of the higgs
      i_j1=0			! index of the hardest jet
      i_j2=0                    ! index of 2nd hardest jet
      pj(:,1)=0
      pj(:,2)=0
      do ixx=1,100
         s1(ixx) = ''
      enddo
      do ixx=1,4
         pj(:,1)=0.0
         pj(:,2)=0.0
         ph(:)  =0.0
      enddo

c     Find the Higgs - we will choose the final higgs
      do ihep=1,nhep
         id=abs(idhep(ihep))
         if(idhep(ihep).eq.25) i_higgs = ihep
         if(whcprg.eq.'NLO') then
            if(ihep.eq.4) then  ! The two jets are at 4 and 5, but
               i_j1 = ihep      ! we don't know yet which one is harder
            elseif(ihep.eq.5) then
               i_j2 = ihep
            endif
         else
            i_j1 = 0
            i_j2 = 0
         endif
      enddo
      
c     We still don't know (in the NLO case) which jet is
c     harder, the 'Born' one or the POWHEG one
      if(whcprg.eq.'NLO') then
         if(i_j1.gt.0) pj(:,1)=phep(1:4,i_j1)
         if(i_j2.gt.0) pj(:,2)=phep(1:4,i_j2)
         
         call getyetaptmass(pj(:,1),y,eta,pt,m)
         ptj1=pt
         call getyetaptmass(pj(:,2),y,eta,pt,m)
         ptj2=pt
         if(ptj2.gt.ptj1) then
            pj(:,1)=phep(1:4,i_j2)
            pj(:,2)=phep(1:4,i_j1)
         endif
      endif
      
      ph=phep(1:4,i_higgs)
      
c     Call Fastjet to build jets for LHEF and PYTHIA case      
      if(whcprg.ne.'NLO') then
         jetRadius= 0.4d0       
         ptmin = 1d0
         call buildjets(1,jetRadius,ptmin,mjets,ktj,
     1        etaj,rapj,phij,ptrel,pj) 	
      endif
      
cccccccccccccccccccccc
c                    c
c     Make plots     c
c                    c
cccccccccccccccccccccc

c     Total plots
      call getyetaptmass(ph,y,eta,pt,m)
      pt_higgs=pt
      y_higgs=y
      call filld('sigmatot',0.5d0,dsig)
      call filld('H-pt',pt_higgs,dsig)
      call filld('H-y',y_higgs,dsig)
      call filld('ptj1',ktj(1),dsig)
      call filld('Yj1',rapj(1),dsig)
      call filld('ptj2',ktj(2),dsig)
      call filld('Yj2',rapj(2),dsig)
      
c     Transverse momentum plots split into yHiggs bins
      do kxx=1,len_y
         
         write(tmp,"(F6.2)") y_bins1(kxx) ! Write bin edges as strings
         
         t1=lenocc(tmp)         ! trim the strings down
         t2=lenocc(tmp2)
         
         s1(kxx) = tmp2(1:t1)//'-yH-'//tmp(1:t2) ! plot suffices
         ls1=lenocc(s1(kxx))    ! trim down plot suffix
         
         if(kxx.eq.1) then      ! initial cut - take from Y-higgs=-inft up to the lowest bin edge
            if(y_higgs.lt.y_bins1(kxx)) then
               call filld('H-pt'//'--inf-yH-'//tmp,pt_higgs,dsig)
               call filld('ptj1'//'--inf-yH-'//tmp,ktj(1),dsig)
!     call filld('ptj2'//'--inf-yH-'//tmp,ktj(2),dsig)
            endif
            tmp2=tmp
         elseif(kxx.gt.1.and.kxx.lt.len_y) then ! This includes all the different cuts in the middle
            if(y_higgs.gt.y_bins1(kxx - 1).and.y_higgs.lt.y_bins1(kxx)) then	
               call filld('H-pt'//s1(kxx)(1:ls1),pt_higgs,dsig)
               call filld('ptj1'//s1(kxx)(1:ls1),ktj(1),dsig)
!     call filld('ptj2'//s1(kxx)(1:ls1),ktj(2),dsig)
            endif
            tmp2=tmp
         elseif(kxx.eq.len_y) then ! and this cuts Y-higgs > highest bin edge
            if(y_higgs.gt.y_bins1(kxx)) then
               call filld('H-pt'//tmp2(1:t2)//'-yH-inf',pt_higgs,dsig)
               call filld('ptj1'//tmp2(1:t2)//'-yH-inf',ktj(1),dsig)
!     call filld('ptj2'//tmp2(1:t2)//'-yH-inf',ktj(2),dsig)
            endif
         endif
      enddo      	
      end
      
      subroutine getyetaptmass(p,y,eta,pt,mass)
      implicit none
      real * 8 p(4),y,eta,pt,mass,pv
      real *8 tiny
      parameter (tiny=1.d-5)
      y=0.5d0*log((p(4)+p(3))/(p(4)-p(3)))
      pt=sqrt(p(1)**2+p(2)**2)
      pv=sqrt(pt**2+p(3)**2)
      if(pt.lt.tiny)then
         eta=sign(1.d0,p(3))*1.d8
      else
         eta=0.5d0*log((pv+p(3))/(pv-p(3)))
      endif
      mass=sqrt(abs(p(4)**2-pv**2))
      end
      
      
      subroutine buildjets(iflag,rr,ptmin,mjets,kt,eta,rap,phi,
     $     ptrel,pjet)
c     arrays to reconstruct jets, radius parameter rr
      implicit none
      integer iflag,mjets
      real * 8  rr,ptmin,kt(*),eta(*),rap(*),
     1     phi(*),ptrel(3),pjet(4,*)
      include   'hepevt.h'
      include  'LesHouches.h'
      integer   maxtrack,maxjet
      parameter (maxtrack=2048,maxjet=2048)
      real * 8  ptrack(4,maxtrack),pj(4,maxjet)
      integer   jetvec(maxtrack),itrackhep(maxtrack)
      integer   ntracks,njets
      integer   j,k,mu,jb,i
      real * 8 r,palg,pp,tmp
      logical islept
      external islept
      real * 8 vec(3),pjetin(0:3),pjetout(0:3),beta,
     $     ptrackin(0:3),ptrackout(0:3)
      real * 8 get_ptrel
      external get_ptrel
C     - Initialize arrays and counters for output jets
      do j=1,maxtrack
         do mu=1,4
            ptrack(mu,j)=0d0
         enddo
         jetvec(j)=0
         ptrel(j) = 0d0
      enddo      
      ntracks=0
      do j=1,maxjet
         do mu=1,4
            pjet(mu,j)=0d0
            pj(mu,j)=0d0
         enddo
         kt(j)  = 0d0
         eta(j) = 0d0
         rap(j) = 0d0
         phi(j) = 0d0
      enddo
      if(iflag.eq.1) then
C     - Extract final state particles to feed to jet finder
         do j=1,nhep
c     all but the Higgs
            if (isthep(j).eq.1.and..not.idhep(j).eq.25) then
               if(ntracks.eq.maxtrack) then
                  write(*,*) 'analyze: need to increase maxtrack!'
                  write(*,*) 'ntracks: ',ntracks
                  stop
               endif
               ntracks=ntracks+1
               do mu=1,4
                  ptrack(mu,ntracks)=phep(mu,j)
               enddo
               itrackhep(ntracks)=j
            endif
         enddo
      else
         do j=1,nup
            if (istup(j).eq.1.and..not.islept(idup(j))) then
               if(ntracks.eq.maxtrack) then
                  write(*,*) 'analyze: need to increase maxtrack!'
                  write(*,*) 'ntracks: ',ntracks
                  stop
               endif
               ntracks=ntracks+1
               do mu=1,4
                  ptrack(mu,ntracks)=pup(mu,j)
               enddo
               itrackhep(ntracks)=j
            endif
         enddo
      endif
      if (ntracks.eq.0) then
         mjets=0
         return
      endif
C     --------------------------------------------------------------------- C
C     R = 0.7   radius parameter
c     palg=1 is standard kt, -1 is antikt
      palg=-1
      r=rr
c     ptmin=20d0 
      call fastjetppgenkt(ptrack,ntracks,r,palg,ptmin,pjet,njets,
     $     jetvec)
      mjets=njets
      if(njets.eq.0) return
c     check consistency
      do k=1,ntracks
         if(jetvec(k).gt.0) then
            do mu=1,4
               pj(mu,jetvec(k))=pj(mu,jetvec(k))+ptrack(mu,k)
            enddo
         endif
      enddo
      tmp=0
      do j=1,mjets
         do mu=1,4
            tmp=tmp+abs(pj(mu,j)-pjet(mu,j))
         enddo
      enddo
      if(tmp.gt.1d-4) then
         write(*,*) ' bug!'
      endif
C     --------------------------------------------------------------------- C
C     - Computing arrays of useful kinematics quantities for hardest jets - C
C     --------------------------------------------------------------------- C
      do j=1,mjets
         call getyetaptmass(pjet(:,j),rap(j),eta(j),kt(j),tmp)
         phi(j)=atan2(pjet(2,j),pjet(1,j))
      enddo
      
c     loop over the hardest 3 jets
      do j=1,min(njets,3)
         do mu=1,3
            pjetin(mu) = pjet(mu,j)
         enddo
         pjetin(0) = pjet(4,j)         
         vec(1)=0d0
         vec(2)=0d0
         vec(3)=1d0
         beta = -pjet(3,j)/pjet(4,j)
         call mboost(1,vec,beta,pjetin,pjetout)         
c     write(*,*) pjetout
         ptrel(j) = 0
         do i=1,ntracks
            if (jetvec(i).eq.j) then
               do mu=1,3
                  ptrackin(mu) = ptrack(mu,i)
               enddo
               ptrackin(0) = ptrack(4,i)
               call mboost(1,vec,beta,ptrackin,ptrackout) 
               ptrel(j) = ptrel(j) + get_ptrel(ptrackout,pjetout)
            endif
         enddo
      enddo
      end
      
      function islept(j)
      implicit none
      logical islept
      integer j
      if(abs(j).ge.11.and.abs(j).le.15) then
         islept = .true.
      else
         islept = .false.
      endif
      end
      
      function get_ptrel(pin,pjet)
      implicit none
      real * 8 get_ptrel,pin(0:3),pjet(0:3)
      real * 8 pin2,pjet2,cth2,scalprod
      pin2  = pin(1)**2 + pin(2)**2 + pin(3)**2
      pjet2 = pjet(1)**2 + pjet(2)**2 + pjet(3)**2
      scalprod = pin(1)*pjet(1) + pin(2)*pjet(2) + pin(3)*pjet(3)
      cth2 = scalprod**2/pin2/pjet2
      get_ptrel = sqrt(pin2*abs(1d0 - cth2))
      end

      
