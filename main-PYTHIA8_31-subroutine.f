      subroutine main_pythia8
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer j,k,l,m,iret
      integer maxev
      common/mcmaxev/maxev
      real * 8 powheginput,scalupfac
      external powheginput

      character * 20 pwgprefix
      integer lprefix,ihep,id
      real *8 y,eta,pt,mass
      common/cpwgprefix/pwgprefix,lprefix

      logical analysis_jetveto
      common/canalysis_jetveto/analysis_jetveto

      logical ini
      data    ini/.true./
      save    ini

      character *20 eventfile
      integer iun

      WHCPRG='PYTHIA'

c Only want to 
      if(ini) then

         call init_hist

         call opencountunit(maxev,iun)
         write(*,*) 'iun:',iun

         call lhefreadhdr(iun)

         call pythia_init
         
         ini=.false.

      endif 

      call lhefreadev(iun)
      write(*,*) 'iret:',iret
      call pythia_next(iret)
      call pythia_to_hepevt(nmxhep,nhep,isthep,idhep,jmohep,
     1           jdahep,phep,vhep)
      


      do ihep=1,nhep
         id=abs(idhep(ihep))
         if(idhep(ihep).eq.25) then 
         	ihiggs = ihep    				 		
!         	write(*,*) 'higgs at ',ihiggs - 1 
         endif										 
      enddo											 

!       write(*,*) 'final higgs at',ihiggs - 1

! for some reason, idhep is shifted forward one, so the higgs is actually
! at i_higgs - 1. However, the higgs momentum is phep(1:4,i_higgs)

      write(*,*) 'higgs 4 momentum:',phep(1:4,ihiggs)

      call getyetaptmass(phep(1:4,ihiggs),y,eta,pt,mass)
      write(*,*) 'hpt after showering:',pt
      write(*,*) 'hy after showering:',y
      write(*,*)

C       write(*,*) 'maxev:',maxev

C       nevhep=0

C       do l=1,maxev

C          call lhefreadev(iun)

C          nevhep=nevhep+1

C          do m=1,1 !ER:
C ! Insist to shower this event;
C             call pythia_next(iret)

C             if(iret.eq.-1) then
C                nevhep=nevhep-1
C                print*, 'EOF'
C                print*, nevhep
C                print*, l
C                goto 123
C             endif
           
C             if(iret.ne.1) then
C                write(*,*) ' return code ',iret
C                if(m.eq.1) then
C                   write(*,*) ' Pythia could not shower this event'
C                   call printleshouches
C                endif
C c               write(*,*) ' retry: ',l
C             else
C c               write(*,*) ' done: ',l
C                exit
C             endif
C          enddo

C          if(iret.eq.1) then
C             call pythia_to_hepevt(nmxhep,nhep,isthep,idhep,jmohep,
C      1           jdahep,phep,vhep)
C C             if(nevhep.lt.6) then
C C                do j=1,nhep
C C                   write(*,100)j,isthep(j),idhep(j),jmohep(1,j),
C C      1           jmohep(2,j),jdahep(1,j),jdahep(2,j), (phep(k,j),k=1,5)
C C                enddo
C C             endif

C             call pyanal

C             if(nevhep.gt.0.and.mod(nevhep,20000).eq.0) then
C                write(*,*)'# of events processed=',nevhep
C                write(*,*)'# of events generated=',nevhep
C                call pyaend
C             endif 
C             if(nevhep.gt.0.and.mod(nevhep,20000).eq.0) then
C                write(*,*)'# of events processed=',nevhep
C                write(*,*)'# of events generated=',nevhep
C             endif 
C          endif


C       enddo

C  123  continue

C       write(*,*) 'At the end NEVHEP is ',nevhep

C       call pythia_stat
C C !:      write(*,*) 'At the end: #warnings= ',mstu(27),' #errors= ',mstu(23)
C c---user's terminal calculations


C       call pyaend

 100  format(i4,2x,i5,2x,i5,2x,i4,1x,i4,2x,i4,1x,i4,2x,5(d10.4,1x))
      end

      subroutine pyanal
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
c     check parameters
      logical verbose
      parameter (verbose=.false.)
      real * 8 powheginput
      external powheginput
c     !ER:
c      nevhep=nevhep+1
c     print*, xwgtup,xsecup(1),nup
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end

      subroutine getmaxev(maxev)
      integer maxev
C--- Opens input file and counts number of events, setting MAXEV;
      call opencount(maxev)
      end

      subroutine pyaend
      character * 6 vetoname
      character * 20 pwgprefix
      character * 100 filename
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      include 'pwhg_rnd.h'
      logical analysis_jetveto
      common/canalysis_jetveto/analysis_jetveto
      
      if(analysis_jetveto) then
         vetoname=trim(adjustl('_veto'))
      else
         vetoname=trim(adjustl(''))
      endif
      if(rnd_cwhichseed.ne.'none') then
         filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA8-output-'
     1        //rnd_cwhichseed//vetoname
      else
         filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA8-output'//vetoname
      endif

      call pwhgsetout
      call pwhgtopout(filename)
      end


