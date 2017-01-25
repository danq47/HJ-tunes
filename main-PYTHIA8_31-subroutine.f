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
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix

      logical analysis_jetveto
      common/canalysis_jetveto/analysis_jetveto

      logical ini
      data    ini/.true./
      save    ini

      character *20 eventfile
      integer iun

      WHCPRG='PYTHIA'

      if(ini) then

         call init_hist

         call opencountunit(maxev,iun)
         write(*,*) 'iun:',iun

         call lhefreadhdr(iun)

         call pythia_init
         
         ini=.false.

      endif

      nevhep=0

      do l=1,maxev

         call lhefreadev(iun)

         nevhep=nevhep+1

         do m=1,1 !ER:
c Insist to shower this event;
            call pythia_next(iret)

            if(iret.eq.-1) then
               nevhep=nevhep-1
               print*, 'EOF'
               print*, nevhep
               print*, l
               goto 123
            endif
           
            if(iret.ne.1) then
               write(*,*) ' return code ',iret
               if(m.eq.1) then
                  write(*,*) ' Pythia could not shower this event'
                  call printleshouches
               endif
c               write(*,*) ' retry: ',l
            else
c               write(*,*) ' done: ',l
               exit
            endif
         enddo

         if(iret.eq.1) then
            call pythia_to_hepevt(nmxhep,nhep,isthep,idhep,jmohep,
     1           jdahep,phep,vhep)
C             if(nevhep.lt.6) then
C                do j=1,nhep
C                   write(*,100)j,isthep(j),idhep(j),jmohep(1,j),
C      1           jmohep(2,j),jdahep(1,j),jdahep(2,j), (phep(k,j),k=1,5)
C                enddo
C             endif

            call pyanal

            if(nevhep.gt.0.and.mod(nevhep,20000).eq.0) then
               write(*,*)'# of events processed=',nevhep
               write(*,*)'# of events generated=',nevhep
               call pyaend
            endif 
            if(nevhep.gt.0.and.mod(nevhep,20000).eq.0) then
               write(*,*)'# of events processed=',nevhep
               write(*,*)'# of events generated=',nevhep
            endif 
         endif


      enddo

 123  continue

      write(*,*) 'At the end NEVHEP is ',nevhep

      call pythia_stat
!:      write(*,*) 'At the end: #warnings= ',mstu(27),' #errors= ',mstu(23)
c---user's terminal calculations


      call pyaend

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


