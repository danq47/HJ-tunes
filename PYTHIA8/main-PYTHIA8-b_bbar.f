      program main_pythia8
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
      include 'nlegborn.h'
      include 'pwhg_rad.h'
      include 'pwhg_rwl.h'
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer j,k,l,m,iret,kloop,kpy8
      integer maxev,maxevin
      common/mcmaxev/maxev
      integer canveto
      integer allrad,nlowhich
      real * 8 vetoscaletp,vetoscaletm,
     1         vetoscalewp,vetoscalewm
      common/resonancevetos/vetoscaletp,vetoscaletm,
     1                      vetoscalewp,vetoscalewm,canveto
      integer ntpdec,ntmdec,tpiddec(8),tmiddec(8)
      real * 8 tpdecsc,tmdecsc,tppdec(4,8),tmpdec(4,8)
      real * 8 wpdecsc,wmdecsc
      common/ctptmdec/tpdecsc,tmdecsc,tppdec,tmpdec,
     1     ntpdec,ntmdec,tpiddec,tmiddec
      integer py8tune,nohad
      common/cpy8tune/py8tune,nohad
      real * 8 powheginput,scalupfac,phepdot
      logical weveto,guessres,stripres,nores,hvq_gen
      real * 8 ub_btilde_corr, ub_remn_corr, ub_corr
      external powheginput
      common/guessresStats/n1,n2,n3,n4,n5,n6
      integer n1,n2,n3,n4,n5,n6
      integer iun

      WHCPRG='PYTHIA'

      py8tune = powheginput("#py8tune")
      nohad = powheginput("#nohad")
c     canveto 1 means Pythia is vetoig radiation in resonance decay. We
c     decided not to use this.
      canveto = 0
c     read allrad (default is 1)
      allrad = powheginput("#allrad")
c     whether we are analyzing a run with no resonance info
      nores = powheginput("#nores").eq.1
c     whether we are introducing resonance information guessed from kinematics
      guessres = powheginput("#guessres").ne.0
c     whether we cancel or not resonance information from the LH event
      stripres =  powheginput("#stripres").eq.1
c     this flag must be present if we are running the hvq generator.
c     the hvq program generates charm, bottom or top depending upon the qmass value
      hvq_gen = powheginput("#qmass") > 0
c     read nlowhich (default is 0). This flag ment:
c     0 all that can radiate is active
c     1 only radiation in production
c     2 only radiation in ??? (check in ttb_NLO_dec
c     It is relevant for the ttb_NLO_dec generator. Here
c     we keep it zero.
      nlowhich=0

      if(powheginput("#LOevents") == 1) then
         allrad = 0
         weveto = .false.
         canveto = 0
         guessres = .false.
         stripres = .false.
      elseif(nores) then
c     this covers the nores case. allrad must be 0 now (see previous if)
         allrad = powheginput("#allrad")
         if(allrad == 1) then
            write(*,*) " you can't use allrad=1 with nores=1, exiting ..."
         endif
         allrad = 0
         if(stripres) then
            write(*,*) " you can't use stripres with nores=1, exiting ..."
         endif
         if(hvq_gen) then
            write(*,*)
     1           " you can't use nores=1 with the hvq generator, exiting ..."
            call exit(-1)
         endif
         if(guessres) then
            weveto = .true.
         else
            weveto = .false.
         endif
      elseif(hvq_gen) then
c guessres by default should be off
         guessres = powheginput("#guessres").eq.1
         if(guessres) then
            write(*,*) " hvq_gen = 1, guessres = 1: can't be, exiting ..."
            call exit(-1)
         endif
         if(stripres) then
            write(*,*) " hvq_gen = 1, stripres = 1: can't be, exiting ..."
            call exit(-1)
         endif
         weveto = .false.
c     The hvq program does not generate radiation in decay.
c     In all cases Pythia should not veto radiation in decays. The
c     following should not matter.
         nlowhich = 1
      elseif(allrad.ne.0) then
c     allrad 1 is for running bbar4l in default mode, and also for ttb_NLO_dec
         allrad=1
c     weveto means the we examine radiation in resonance decays to set Pythia
c     showering scales for resonances. If false, Pythia does not veto radiation
c     in resonance decays. Only production radiation is vetoed according to scalup
         weveto = .true.
         if(stripres) then
            write(*,*) " you can't use stripres with allrad=1, exiting ..."
            call exit(-1)
         endif
         guessres = .false.
      else
c     this covers the allrad 0 case
c     Now guessres is by default false
         guessres = powheginput("#guessres").eq.1
         if(guessres .and. .not. stripres) then
            write(*,*) " allrad = 0, guessres = 1 but stripres = 0: can't be, exiting ..."
            call exit(-1)
         endif
         if(stripres .and. .not. guessres) then
            weveto = .false.
         else
            weveto = .true.
         endif
      endif


      if(powheginput('#pythiaveto') == 1) then
         weveto=.false.
         canveto=1
      endif
      

      scalupfac=powheginput('#scalupfac')
      if(scalupfac.lt.0) scalupfac=1

c read in btilde and remn corrections factors (used together with ubexcess_correct at the generation stage)

      ub_btilde_corr = powheginput('#ub_btilde_corr')
      if (ub_btilde_corr < 0d0) then
        ub_btilde_corr = 1d0
      endif
      ub_remn_corr = powheginput('#ub_remn_corr')
      if (ub_remn_corr < 0d0) then
        ub_remn_corr = 1d0
      endif

      call init_hist

      call opencountunit(maxev,iun)
      write(*,*) 'iun:',iun

      call lhefreadhdr(iun)

      if(powheginput("#pyMEC").eq.0) then
         call pythia_option("TimeShower:MEcorrections = off");
      endif

      if(powheginput("#pyMEaf").eq.0) then
         call pythia_option("TimeShower:MEafterFirst = off");
      endif

      call pythia_init



      nevhep=0
      kpy8 = 0
      maxevin = powheginput('#maxev')
      if (maxevin>0.and.maxevin<=maxev) then
        maxev = maxevin
      endif
      do l=1,maxev

         call lhefreadev(iun)
c rescale the weight of the event depending on the rad_type (1..btilde, 2..remn)
c   using the ub_..._corr factors
         if (rad_type == 1) then
            ub_corr = ub_btilde_corr
         else if (rad_type == 2) then
            ub_corr = ub_remn_corr
         else 
            ub_corr = 1d0
         endif
         rwl_weights(1:rwl_num_weights)=
     c         ub_corr * rwl_weights(1:rwl_num_weights)

c strip the resonance assignment        
         if(stripres) then
           call strip_resonances
         endif

c guess resonance assignment
         if(guessres) then
            call guess_resonances
         endif
       
         scalup=scalupfac*scalup
         if(nlowhich.eq.0) then
            if(allrad.eq.1) then
               call findresscale( 6,vetoscaletp)
               call findresscale(-6,vetoscaletm)
c 7 and 9 are the fermions from W+ or W-; if they are hadrons,
c find the veto scales.
               if(abs(idup(7)).lt.5) then
                  call findresscale(24,vetoscalewp)
               else
                  vetoscalewp = 1d30
               endif
               if(abs(idup(9)).lt.5) then
                  call findresscale(-24,vetoscalewm)
               else
                  vetoscalewm = 1d30
               endif
            else
               vetoscaletp = scalup
               vetoscaletm = scalup
            endif
         elseif(nlowhich.eq.1) then
c            canveto = 0
         elseif(nlowhich.eq.2) then
c            canveto = 1
            vetoscaletp = scalup
            vetoscaletm = scalup
         elseif(nlowhich.eq.3) then
            call findresscale(-24,vetoscalewm)
            if(dabs(scalup-vetoscalewm).gt.0.01d0) then
               if(l.lt.10) then
                  write(*,*) 'nlowhich = 3: problem with veto scales'
                  write(*,*) scalup,vetoscalewm
               endif
            endif
            vetoscalewm=scalup
            vetoscalewp=1d30
            scalup = sqrt(2*phepdot(pup(:,1),pup(:,2)))
         elseif(nlowhich.eq.4) then
c            canveto = 1
            vetoscaletp = scalup
c            call findresscale( 6,scalup)
c            write(*,*) 'vetoscaletp/scalup', vetoscaletp/scalup
            vetoscaletm = pup(5,4)
            scalup = sqrt(2*phepdot(pup(:,1),pup(:,2)))
            !ER: try to set this to a (much) bigger value and see what happens
            ! on tbar observables
         else
            write(*,*) 'nlowhich = ',nlowhich,' not yet implemented'
            call pwhg_exit(-1)
         endif

c useful to check that py8 is doing what we are asking
         if(l.lt.10) then
            write(*,*) 'veto scales: ',l,scalup,vetoscaletp,vetoscalewp,
     1           vetoscaletm,vetoscalewm,canveto
         endif
 
         m=1
         call copylh
         do kloop=1,1000000
c Insist to shower this event;
            call pythia_next(iret)
            call resetlh
            call checklh
            kpy8 = kpy8+1
            
            if(iret.ne.1) then
               write(*,*) ' return code ',iret
               if(m.eq.1) then
                  write(*,*) ' Pythia could not shower this event'
                  call printleshouches
               endif
               write(*,*) ' retry'
               if(m.gt.4) then
                  write(*,*) ' after 5 attempts'
                  write(*,*) ' abandoning the event'
                  exit
               endif
               m=m+1
            else
               call pythia_to_hepevt(nmxhep,nhep,isthep,
     1              idhep,jmohep,jdahep,phep,vhep)
               if(weveto) then
                  call getdechardness(1,tpdecsc,ntpdec,tpiddec,tppdec)
                  call getdechardness(-1,tmdecsc,ntmdec,tmiddec,tmpdec)
c the following is probably not needed (was used at some point to
c produce some plots)
                  if(tpdecsc>0) then
                     call boost2reson4(tppdec,ntpdec,tppdec,tppdec)
                  endif
                  if(tmdecsc>0) then
                     call boost2reson4(tmpdec,ntmdec,tmpdec,tmpdec)
                  endif
c 7 and 9 are the fermions from W+ or W-; if they are hadrons,
c find the veto scales.
                  if(abs(idup(7)).lt.5) then
                     call getdechardnessw(1,wpdecsc)
                  else
                     wpdecsc = 0
                  endif
                  if(abs(idup(9)).lt.5) then
                     call getdechardnessw(1,wmdecsc)
                  else
                     wmdecsc = 0
                  endif

                  if(nlowhich.eq.0) then
                     if(allrad.eq.1) then
                        if(tpdecsc.gt.vetoscaletp .or.
     1                     tmdecsc.gt.vetoscaletm .or.
     2                     wpdecsc.gt.vetoscalewp .or.
     3                     wmdecsc.gt.vetoscalewm ) then
                           continue
                        else
                           exit
                        endif
                     else
                        if(tpdecsc.gt.scalup .or.
     1                     tmdecsc.gt.scalup .or.
     1                     wpdecsc.gt.scalup .or.
     1                     wmdecsc.gt.scalup ) then
                           continue
                        else
                           exit
                        endif
                     endif                     
                  elseif(nlowhich.eq.1) then
                     exit
                  elseif(nlowhich.eq.4) then
                     if(tpdecsc.gt.vetoscaletp) then
                        continue
                     else
                        exit
                     endif
                  endif
               else
                  exit
               endif
            endif
         enddo

         if(iret.eq.1) then
            if(nevhep.lt.10) then
               do j=1,nhep
                  write(*,100)j,isthep(j),idhep(j),jmohep(1,j),
     1           jmohep(2,j),jdahep(1,j),jdahep(2,j), (phep(k,j),k=1,5)
               enddo
            endif
            call pyanal
            if(nevhep.gt.0.and.mod(nevhep,1000).eq.0) then
               write(*,*)'# of events processed=',kpy8
               write(*,*)'# of events generated=',nevhep
               call pyaend
            endif 
         endif
      enddo

      write(*,*) 'At the end NEVHEP is ',nevhep
      if (guessres) then
         write(*,*) 'guessres stats: ', n1, n2, n3, n4, n5, n6
      endif
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
      nevhep=nevhep+1
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup 
      end

      subroutine pyaend
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      include 'pwhg_rnd.h'

      call pwhgsetout

      if(rnd_cwhichseed.eq.'none') then
         call pwhgtopout(pwgprefix(1:lprefix)//
     1                   'POWHEG+PYTHIA8-output')
      else
         call pwhgtopout(pwgprefix(1:lprefix)//
     1                   '-'//rnd_cwhichseed //'-'//
     2                   'POWHEG+PYTHIA8-output')
      endif

      end

      subroutine printleshouches
c useful for debugging
      call lhefwritev(6)
      end

c...lhefeader(nlf)
c...writes event information to a les houches events file on unit nlf. 
      subroutine lhefwritev(nlf)
      implicit none
      integer nlf
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      integer i,j
      write(nlf,'(a)')'<event>'
      write(nlf,210) nup,idprup,xwgtup,scalup,aqedup,aqcdup
      do 200 i=1,nup
         write(nlf,220) idup(i),istup(i),mothup(1,i),
     & mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     & vtimup(i),spinup(i)
 200  continue
      write(nlf,'(a)')'</event>'      
 210  format(1p,2(1x,i8),4(1x,e12.5))
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)
      end


      subroutine findresscale(idres,scale)
      implicit none
      include "LesHouches.h"
      integer idres
      real * 8 scale
      real * 8 p(0:3,3),p0(0:3),dotp
      integer nres,ids(3),idb,idg,idw,idq,ida
      real * 8 yq,ya,csi,q2
      integer j,k
      nres=0
      do j=3,nup
         if(idup(j).eq.idres) then
            p0(1:3) = pup(1:3,j)
            p0(0) = pup(4,j)
            do k=3,nup
               if(mothup(1,k).eq.j) then
                  if(nres.ge.3) then
                     write(*,*)
     1                    ' findresscale: error: more than 3 '
     1                    //'particles in resonance decay'
                     call exit(-1)
                  endif
                  nres = nres+1
                  p(1:3,nres) = pup(1:3,k)
                  p(0,nres) = pup(4,k)
                  ids(nres) = idup(k)
               endif
            enddo
            goto 10
         endif
      enddo
 10   continue
      if(nres == 0) then
c     No resonance found, set scale to high value
c     Pythia will shower any MC generated resonance unrestricted
         scale = 1d30
         return
      elseif(nres.lt.3) then
c No radiating resonance found
         scale = 0.8d0
         return
      endif
      call boost2reson(p0,nres,p,p)
      if(abs(idres).eq.6) then
         do j=1,3
            if(abs(ids(j)).eq.24) then
               idw=j
            elseif(abs(ids(j)).eq.5) then
               idb=j
            elseif(abs(ids(j)).eq.21) then
               idg=j
            endif
         enddo
c         scale = sqrt( (
c     1        (p(1,idg)*p(2,idw)-p(2,idg)*p(1,idw))**2+
c     2        (p(2,idg)*p(3,idw)-p(3,idg)*p(2,idw))**2+
c     3        (p(3,idg)*p(1,idw)-p(1,idg)*p(3,idw))**2 )
c     4         /
c     4        (p(1,idw)**2+p(2,idw)**2+p(3,idw)**2)  )
c
c
c         if(scale.lt.5d0) then
c            scale = min(5d0,p(0,idg))
c         endif
         scale = sqrt( 2 * dotp(p(:,idg),p(:,idb))*
     1        p(0,idg)/p(0,idb) )
      elseif(abs(idres).eq.24) then
         do j=1,3
            if(ids(j).eq.21) then
               idg=j
            elseif(ids(j).gt.0) then
               idq=j
            elseif(ids(j).lt.0) then
               ida=j
            endif
         enddo
         p0 = p(:,idg)+p(:,idq)+p(:,ida)
         q2 = dotp(p0,p0)
         csi = 2*p(0,idg)/sqrt(q2)
         yq = 1 - dotp(p(:,idg),p(:,idq))/(p(0,idg)*p(0,idq))
         ya = 1 - dotp(p(:,idg),p(:,ida))/(p(0,idg)*p(0,ida))
         scale = sqrt(min(1-yq,1-ya)*csi**2*q2/2)
c         print*, 'scaleq,scalea = ',sqrt((1-ya)*csi**2*q2/2),sqrt((1-yq)*csi**2*q2/2)
c         print*, 'scale,scalup = ',scale,scalup
      endif
      end


      subroutine getdechardness(ichtop,hardness,nmoms,iddec,pmoms)
c ichtop = +- 1, is the charge of the top.
c it returns the hardness of b radiation
      implicit none
      integer ichtop
      real * 8 hardness
      integer nmoms
      integer iddec(8)
      real * 8 pmoms(4,8)
      include 'hepevt.h'
      integer jhep
      integer i_top,i_b,i_w,i_g,j,k,wid,bid,tid
      real * 8 pchain(4,3)
      real * 8 phepdot
      nmoms = 0
      pmoms = 0
      iddec = 0
      tid = 6*ichtop
      wid = 24*ichtop
      bid = 5*ichtop
c     find last top in record
      i_top = -1
      do jhep=1,nhep
         if(idhep(jhep).eq.tid) then
            i_top = jhep
         endif
      enddo
      if(i_top == -1) then
         hardness = 0
         return
      endif
      pchain(:,1)=phep(1:4,i_top)
      pmoms(:,1) = pchain(:,1)
      iddec(1)=tid
      nmoms = 1
c look for top direct sons
      if(jdahep(2,i_top)-jdahep(1,i_top).eq.1) then
         i_w = jdahep(1,i_top)
         i_b = jdahep(2,i_top)
         if(idhep(i_w).ne.wid) then
            write(*,*) ' top did not go in W!'
            goto 998
         endif
         if(idhep(i_b).ne.bid) then
            write(*,*) ' top did not go in b!'
            goto 998
         endif
         nmoms = nmoms+1
         pmoms(:,nmoms) = phep(1:4,i_w)
         iddec(nmoms) = wid
         nmoms = nmoms+1
         pmoms(:,nmoms) = phep(1:4,i_b)
         iddec(nmoms) = bid
         if(jdahep(2,i_b)-jdahep(1,i_b).gt.1) then
            write(*,*) ' found b-> more than 2 particles'
            goto 998
         elseif(idhep(jdahep(1,i_b)).eq.bid
     1        .and.idhep(jdahep(2,i_b)).eq.21) then
c     the b has radiated a gluon
            pchain(:,2) = phep(1:4,jdahep(1,i_b))
            pchain(:,3) = phep(1:4,jdahep(2,i_b))
            nmoms = nmoms+1
            pmoms(:,nmoms) = pchain(:,2)
            iddec(nmoms) = bid
            nmoms = nmoms+1
            pmoms(:,nmoms) = pchain(:,3)
            iddec(nmoms) = 21
         else
            hardness = -1
            return
         endif
c now pchain contains the 4-momenta of the top, and the b-g pair
         call boost2reson4(pchain,3,pchain,pchain)
         
         hardness = sqrt( 2 * phepdot(pchain(:,2),pchain(:,3))
     1        * pchain(4,3)/pchain(4,2) )
         return
      elseif(jdahep(2,i_top)-jdahep(1,i_top).eq.2) then         
c here we have W b g 
         if(.not.(idhep(jdahep(1,i_top)).eq.wid
     1        .and.idhep(jdahep(1,i_top)+1).eq.bid
     2        .and.idhep(jdahep(2,i_top)).eq.21)) then
            write(*,*) ' was not expecting this!'
            goto 998
         endif
         i_w = jdahep(1,i_top)
         i_b = i_w+1
         i_g = i_b+1
         nmoms = nmoms+1
         pmoms(:,nmoms) = phep(1:4,i_w)
         iddec(nmoms) = wid
         nmoms = nmoms+1
         pmoms(:,nmoms) = phep(1:4,i_b)
         iddec(nmoms) = bid
         nmoms = nmoms+1
         pmoms(:,nmoms) = phep(1:4,i_g)
         iddec(nmoms) = 21
c see if b goes into b g
         if(jdahep(2,i_b)-jdahep(1,i_b).gt.1) then
            write(*,*) ' found b-> more than 2 particles'
            goto 998            
         endif
         if(idhep(jdahep(1,i_b)).eq.bid
     1        .and.idhep(jdahep(2,i_b)).eq.21) then
c     the b has radiated a gluon
            pchain(:,2) = phep(1:4,jdahep(1,i_b))
            pchain(:,3) = phep(1:4,jdahep(2,i_b))
            nmoms = nmoms+1
            pmoms(:,nmoms) = pchain(1:4,2)
            iddec(nmoms) = bid
            nmoms = nmoms+1
            pmoms(:,nmoms) = pchain(1:4,3)
            iddec(nmoms) = 21
            call boost2reson4(pchain,2,pchain(1,2),pchain(1,2))
            hardness = sqrt( 2 * phepdot(pchain(:,2),pchain(:,3))
     1        * pchain(4,3)/pchain(4,2) )
         else
            hardness = -1
         endif
         if(jdahep(2,i_g)-jdahep(1,i_g).eq.1) then
            pchain(:,2) = phep(1:4,jdahep(1,i_g))
            pchain(:,3) = phep(1:4,jdahep(2,i_g))
            nmoms = nmoms+1
            pmoms(:,nmoms) = pchain(1:4,2)
            iddec(nmoms) = idhep(jdahep(1,i_g))
            nmoms = nmoms+1
            pmoms(:,nmoms) = pchain(1:4,3)
            iddec(nmoms) = idhep(jdahep(2,i_g))
            call boost2reson4(pchain,2,pchain(1,2),pchain(1,2))
            hardness = max(hardness,
     1           sqrt( 2 * phepdot(pchain(:,2),pchain(:,3))
     2        * (pchain(4,3)*pchain(4,2))
     3           /(pchain(4,3)**2+pchain(4,2)**2)))
         elseif(jdahep(2,i_g)-jdahep(1,i_g).gt.1) then
            write(*,*) ' found g-> more than 2 particles'
            goto 998            
         endif
         return
      else
         write(*,*) ' was not expecting this!'
         goto 998
      endif
      goto 999
 998  continue
      write(*,*) 'top=',i_top
      do j=1,nhep
         write(*,100)j,isthep(j),idhep(j),jmohep(1,j),
     1        jmohep(2,j),jdahep(1,j),jdahep(2,j), (phep(k,j),k=1,5)
      enddo
      call exit(-1)
 100  format(i4,2x,i5,2x,i5,2x,i4,1x,i4,2x,i4,1x,i4,2x,5(d10.4,1x))
 999  end


      subroutine getdechardnessw(ichw,hardness)
c ichtop = +- 1, is the charge of the top.
c it returns the hardness of b radiation
      implicit none
      integer ichw
      real * 8 hardness
      integer nmoms
      integer iddec(8)
      real * 8 pw(4),h1,h2,h3
      include 'hepevt.h'
      integer jhep
      integer i_w,j,k,wid
      wid = 24*ichw
      i_w=0
c find last W in record
      do jhep=1,nhep
         if(idhep(jhep).eq.wid) then
            i_w = jhep
         endif
      enddo
      if(i_w.eq.0) then
         write(*,*) 'getdechardnessw: could not find the W! exiting ...'
         call exit(-1)
      endif
      pw(:)=phep(1:4,i_w)
c look for W direct sons
      if(jdahep(2,i_w)-jdahep(1,i_w).eq.1) then
         call findpy8dec(jdahep(2,i_w),pw,h1)
         call findpy8dec(jdahep(1,i_w),pw,h2)
         hardness=max(h1,h2)
         if(hardness.gt.1d0) hardness = 1d25
      elseif(jdahep(2,i_w)-jdahep(1,i_w).eq.2) then
         call findpy8dec(jdahep(1,i_w),pw,h1)
         call findpy8dec(jdahep(1,i_w)+1,pw,h2)
         call findpy8dec(jdahep(2,i_w),pw,h3)
         hardness = max(h1,h2,h3)
      else
         write(*,*) 'getdechardnessw: was not expecting this!'
         goto 998
      endif
      goto 999
 998  continue
      write(*,*) 'getdechardnessw: wid=',wid
      do j=1,nhep
         write(*,100)j,isthep(j),idhep(j),jmohep(1,j),
     1        jmohep(2,j),jdahep(1,j),jdahep(2,j), (phep(k,j),k=1,5)
      enddo
      call exit(-1)
 100  format(i4,2x,i5,2x,i5,2x,i4,1x,i4,2x,i4,1x,i4,2x,5(d10.4,1x))
 999  end


      subroutine findpy8dec(j,p0,h)
      implicit none
      include 'hepevt.h'
      integer j
      real * 8 p0(4),h
      real * 8 pchain(4,3)
      real * 8 phepdot
      if(jdahep(2,j).eq.jdahep(1,j)) then
         h = 0
         return
      endif

      pchain(:,1) = p0
      pchain(:,2) = phep(1:4,jdahep(1,j))
      pchain(:,3) = phep(1:4,jdahep(2,j))

      if(jdahep(2,j)-jdahep(1,j).eq.1) then
         call boost2reson4(p0,3,pchain(:,:),pchain(:,:))
         h = sqrt( 2 * phepdot(pchain(:,2),pchain(:,3))
     2        * (pchain(4,3)*pchain(4,2))
     3        /(pchain(4,3)**2+pchain(4,2)**2))
      else
         write(*,*) 'findpy8dec: was not expecting this!'
         call exit(-1)
      endif

      end



      subroutine copylh
      implicit none
      include 'LesHouches.h'
      integer idbmupz,pdfgupz,pdfsupz,idwtupz,nprupz,lprupz
      double precision ebmupz,xsecupz,xerrupz,xmaxupz
      common /heprups/ idbmupz(2),ebmupz(2),pdfgupz(2),pdfsupz(2),
     &                idwtupz,nprupz,xsecupz(maxpup),xerrupz(maxpup),
     &                xmaxupz(maxpup),lprupz(maxpup)
      integer nupz,idprupz,idupz,istupz,mothupz,icolupz
      double precision xwgtupz,scalupz,aqedupz,aqcdupz,pupz,vtimupz,spinupz
      common/hepeups/nupz,idprupz,xwgtupz,scalupz,aqedupz,aqcdupz,
     &              idupz(maxnup),istupz(maxnup),mothupz(2,maxnup),
     &              icolupz(2,maxnup),pupz(5,maxnup),vtimupz(maxnup),
     &              spinupz(maxnup)
      idbmupz = idbmup
      pdfgupz = pdfgup
      pdfsupz = pdfsup
      idwtupz = idwtup
      nprupz  = nprup
      lprupz  = lprup



      ebmupz  =       ebmup   
      xsecupz =       xsecup  
      xerrupz =       xerrup  
      xmaxupz =       xmaxup   



      nupz       =      nup   
      idprupz    =      idprup
      idupz      =      idup  
      istupz     =      istup 
      mothupz    =      mothup
      icolupz    =      icolup
                              
      xwgtupz    =      xwgtup
      scalupz    =      scalup
      aqedupz    =      aqedup
      aqcdupz    =      aqcdup
      pupz       =      pup   
      vtimupz    =      vtimup
      spinupz    =      spinup

      end


      
      
      subroutine checklh
      implicit none
      include 'LesHouches.h'
      integer idbmupz,pdfgupz,pdfsupz,idwtupz,nprupz,lprupz
      double precision ebmupz,xsecupz,xerrupz,xmaxupz
      common /heprups/ idbmupz(2),ebmupz(2),pdfgupz(2),pdfsupz(2),
     &                idwtupz,nprupz,xsecupz(maxpup),xerrupz(maxpup),
     &                xmaxupz(maxpup),lprupz(maxpup)
      integer nupz,idprupz,idupz,istupz,mothupz,icolupz
      double precision xwgtupz,scalupz,aqedupz,aqcdupz,pupz,vtimupz,spinupz
      common/hepeups/nupz,idprupz,xwgtupz,scalupz,aqedupz,aqcdupz,
     &              idupz(maxnup),istupz(maxnup),mothupz(2,maxnup),
     &              icolupz(2,maxnup),pupz(5,maxnup),vtimupz(maxnup),
     &              spinupz(maxnup)
      if(sum(abs(idbmupz - idbmup)).ne.0
     1 .or. sum(abs(pdfgupz - pdfgup)).ne.0
     2 .or. sum(abs(pdfsupz - pdfsup)).ne.0
     3 .or. idwtupz - idwtup .ne.0
     4 .or. nprupz  - nprup .ne.0
     5 .or. sum(abs(lprupz  - lprup)).ne.0 ) goto 998



      if(    sum(abs(ebmupz  -       ebmup  ))  .ne.0
     1 .or.  sum(abs(xsecupz -       xsecup ))  .ne.0
     1 .or.  sum(abs(xerrupz -       xerrup ))  .ne.0
     1 .or.  sum(abs(xmaxupz -       xmaxup ))  .ne.0) goto 998


      if( nupz   -    nup     .ne.0
     1 .or. idprupz    -      idprup   .ne.0
     1 .or. sum(abs( idupz      -      idup  )).ne.0
     1 .or. sum(abs( istupz     -      istup )).ne.0
     1 .or. sum(abs( mothupz    -      mothup)).ne.0
     1 .or. sum(abs( icolupz    -      icolup)).ne.0 ) goto 998

      if(    xwgtupz    -      xwgtup   .ne.0
     1 .or.  scalupz    -      scalup   .ne.0
     1 .or.  aqedupz    -      aqedup   .ne.0
     1 .or.  aqcdupz    -      aqcdup   .ne.0
     1 .or.  sum(abs( pupz       -      pup   )).ne.0
     1 .or.  sum(abs( vtimupz    -      vtimup)).ne.0
     1 .or.  sum(abs( spinupz    -      spinup)).ne.0) goto 998
      return
 998  write(*,*) ' checklh: fails ...'
      call exit(-1)

      end


      
      
      
      subroutine resetlh
      implicit none
      include 'LesHouches.h'
      integer idbmupz,pdfgupz,pdfsupz,idwtupz,nprupz,lprupz
      double precision ebmupz,xsecupz,xerrupz,xmaxupz
      common /heprups/ idbmupz(2),ebmupz(2),pdfgupz(2),pdfsupz(2),
     &                idwtupz,nprupz,xsecupz(maxpup),xerrupz(maxpup),
     &                xmaxupz(maxpup),lprupz(maxpup)
      integer nupz,idprupz,idupz,istupz,mothupz,icolupz
      double precision xwgtupz,scalupz,aqedupz,aqcdupz,pupz,vtimupz,spinupz
      common/hepeups/nupz,idprupz,xwgtupz,scalupz,aqedupz,aqcdupz,
     &              idupz(maxnup),istupz(maxnup),mothupz(2,maxnup),
     &              icolupz(2,maxnup),pupz(5,maxnup),vtimupz(maxnup),
     &              spinupz(maxnup)
      idprup = idprupz
      end


      
      
      subroutine pythia_option(string)
      character * (*) string
      character * 1 null
      null=char(0)
      call  pythia_option0(trim(string)//null)
      end

      subroutine strip_resonances
      implicit none
      include 'LesHouches.h'
      integer i
      integer nuptmp,iduptmp(maxnup),istuptmp(maxnup),mothuptmp(2,maxnup),
     #     icoluptmp(2,maxnup)
      real * 8 puptmp(5,maxnup),vtimuptmp(maxnup),spinuptmp(maxnup)
      integer idupnew
      nuptmp=0
      do i=1,nup
        if (abs(idup(i)).ne.6.and.abs(idup(i)).ne.24.and.abs(idup(i)).ne.23) then
          nuptmp = nuptmp+1
          iduptmp(nuptmp) = idup(i)
          istuptmp(nuptmp) = istup(i)
          if (i > 2) then
            mothuptmp(1,nuptmp) = 1
            mothuptmp(2,nuptmp) = 2
          else 
            mothuptmp(:,nuptmp) = mothup(:,i)
          endif
          icoluptmp(:,nuptmp) = icolup(:,i)
          puptmp(:,nuptmp) = pup(:,i)
          vtimuptmp(nuptmp) = vtimup(i)
          spinuptmp(nuptmp) = spinup(i)
        endif
      enddo
      nup=nuptmp
      do i=1,nuptmp
        idup(i) = iduptmp(i)
        istup(i) = istuptmp(i)
        mothup(:,i) = mothuptmp(:,i)
        icolup(:,i) = icoluptmp(:,i)
        pup(:,i) = puptmp(:,i)
        vtimup(i) = vtimuptmp(i)
        spinup(i) = spinuptmp(i)
      enddo
      end

      subroutine guess_resonances
      implicit none
      include 'LesHouches.h'
      integer j,iwp,iwm,ib,iab,ir,it,iat,iret
      real * 8 ptop(1:5),patop(1:5),pz(1:5),ptemp(1:5),pcm(5),
     1  ptopr(1:5),patopr(1:5),ecm,
     1  resfacZ,resfacTTB,resfacTTB_T,resfacTTB_AT, 
     1  ptisr2,ptfsr2b,ptfsr2ab,pttr2b,pttr2ab,r,w1,w2,w3,w4,norm,erad,eem
      real * 8 phepdot,random
      integer guess_resonances_std_res_factors
      real * 8, save :: ph_twidth,ph_tmass,ph_wwidth,ph_wmass,ph_zmass,ph_zwidth
      real * 8 alpha,gfermi,two,four,rt2,pi
      parameter ( alpha= 1/132.50698d0, gfermi = 0.1166390d-4)
      parameter( Two = 2.0d0, Four = 4.0d0 )
      parameter( Rt2   = 1.414213562d0 )
      parameter( Pi = 3.14159265358979323846d0 )
      logical, parameter :: debug = .false.
      real * 8 powheginput
      external powheginput
      logical, save :: ini = .true.
      common/guessresStats/n1,n2,n3,n4,n5,n6
      integer n1,n2,n3,n4,n5,n6
c This attempts to reconstruct the resonance assignment based on the
c  initial, final states and the kinematics of the event
      if (ini) then
         ph_zmass = 91.188d0
         ph_zwidth=2.441d0
 
         ph_tmass=powheginput('#tmass')
         if(ph_tmass<0d0) ph_tmass=172.5d0 ! default value 
         ph_twidth=powheginput('#twidth')
         if(ph_twidth<0d0) ph_twidth=1.3285137647990974 ! as calculated by MCFM (with mb=4.75)
 
         ph_wmass=sqrt(ph_zmass**2/Two+sqrt(ph_zmass**4/Four-Pi/Rt2*alpha/gfermi*ph_zmass**2))
         ph_wwidth=2.04807d0
 
         ini = .false.
      endif
c nup == 8 for events with no radiation and no resonance assignment
c nup == 9 for events with radiation and no resonance assignment
c add -24>-14,13
      call add_resonance(-24,3,5,6) ! at this stage e-, nu, mu+ and nu~ are at positions 3, 4, 5, 6
c add 24>-11,12
      call add_resonance(24,3,4,5) ! at this stage e-, nu 4, 5       
c now -24 and 24 on the position 3 and 4; 5 and -5 on position 9 and 10;
c   the eventual radiation on position 11
      iwp = 3
      iwm = 4
      ib = 9
      iab = 10
      ir = 11
c     event without radiation
      pcm = pup(:,1)+pup(:,2)
      ecm = sqrt( pcm(4)**2 - pcm(1)**2 - pcm(2)**2 - pcm(3)**2 )
      if (nup.eq.10) then ! no radiation
         ir = guess_resonances_std_res_factors()
         n1 = n1 + 1
c quark-gluon channel in which the f.s. quark always comes from the production
      elseif (nup.eq.11.and.idup(ir).ne.21) then
         ir = guess_resonances_std_res_factors()
         n2 = n2 + 1
      elseif (nup.eq.11.and.idup(ir).eq.21) then
c calculate the kT,isr and kT,fsr (for both b's, since the gluon could
c   potentially be radiated from both)
         ptisr2 = pup(1,ir)**2+pup(2,ir)**2
         if(debug) print*,"t: ptisr = ", sqrt(ptisr2)
         erad = phepdot(pcm,pup(:,ir))/ecm
         eem = phepdot(pcm,pup(:,ib))/ecm
         ptfsr2b = 2 * phepdot(pup(:,ir),pup(:,ib)) *
     1     eem*erad/(eem**2+erad**2)
         if(debug) print*,"t: ptfsr(b) = ", sqrt(ptfsr2b)
         eem = phepdot(pcm,pup(:,iab))/ecm
         ptfsr2ab = 2 * phepdot(pup(:,ir),pup(:,iab)) *
     1     eem*erad/(eem**2+erad**2)
         if(debug) print*,"t: ptfsr(b~) = ", sqrt(ptfsr2ab)
         ptop = pup(:,iwp) + pup(:,ib)
         patop = pup(:,iwm) + pup(:,iab)
         ptopr = ptop + pup(:,ir)
         patopr = patop + pup(:,ir)
c find radiation energy in top frame (assuming it is coming from the b)
         erad = phepdot(ptopr,pup(:,ir))/sqrt(phepdot(ptopr,ptopr))
         eem  = phepdot(ptopr,pup(:,ib))/sqrt(phepdot(ptopr,ptopr))
         pttr2b = 2 * phepdot(pup(:,ib),pup(:,ir))* erad*eem/(eem**2+erad**2)
         if(debug) print*,"t: pttr2(b) = ", sqrt(pttr2b)
c find radiation energy in anti-top frame (assuming it is coming from the b)
         erad = phepdot(patopr,pup(:,ir))/sqrt(phepdot(patopr,patopr))
         eem  = phepdot(patopr,pup(:,iab))/sqrt(phepdot(patopr,patopr))
         pttr2ab = 2 * phepdot(pup(:,iab),pup(:,ir))* erad*eem/(eem**2+erad**2)
         if(debug) print*,"t: pttr2(b~) = ", sqrt(pttr2ab)
c calculate the resonance factor assuming it is Z(>WW),b,b resonance history
         pz = pup(:,iwm) + pup(:,iwp)
         resfacZ = (ph_zmass**2)**2/((phepdot(pz,pz)-ph_zmass**2)**2+(ph_zmass*ph_zwidth)**2)
         if(debug) print*,"t: virt(z) = ", sqrt(phepdot(pz,pz))
c calculate the resonance factor assuming it is tt(>WWbb) resonance history
         resfacTTB = (ph_tmass**2)**2/((phepdot(ptop,ptop)-ph_tmass**2)**2+(ph_tmass*ph_twidth)**2)*
     1     (ph_tmass**2)**2/((phepdot(patop,patop)-ph_tmass**2)**2+(ph_tmass*ph_twidth)**2)  
         if(debug) print*,"t: virt(top) = ", sqrt(phepdot(ptop,ptop))
         if(debug) print*,"t: virt(top+rad) = ", sqrt(phepdot(ptopr,ptopr))
c calculate the resonance factor assuming it is tt(>WWbb) resonance
c     history and that the gluon is radiated from the top system
         resfacTTB_T = (ph_tmass**2)**2/((phepdot(ptopr,ptopr)-ph_tmass**2)**2+(ph_tmass*ph_twidth)**2)*
     1     (ph_tmass**2)**2/((phepdot(patop,patop)-ph_tmass**2)**2+(ph_tmass*ph_twidth)**2)
c calculate the resonance factor assuming it is tt(>WWbb) resonance
c     history and that the gluon is radiated from the anti-top system
         if(debug) print*,"t: virt(atop) = ", sqrt(phepdot(patop,patop))
         if(debug) print*,"t: virt(atop+rad) = ", sqrt(phepdot(patopr,patopr))
         resfacTTB_AT = (ph_tmass**2)**2/((phepdot(ptop,ptop)-ph_tmass**2)**2+(ph_tmass*ph_twidth)**2)*
     1     (ph_tmass**2)**2/((phepdot(patopr,patopr)-ph_tmass**2)**2+(ph_tmass*ph_twidth)**2) 
c construct the weights 
c        tt(>WWbb); gluon from ISR
         w1 = (1/ptisr2) * resfacTTB
c        tt(>WWbb); gluon from top
         w2 = (1/pttr2b) * resfacTTB_T
c        tt(>WWbb); gluon from anti top
         w3 = (1/pttr2ab) * resfacTTB_AT
c        Z(>WW),b,b; gluon from production (ISR or FSR from b(b~))
         w4 = (1/ptisr2+1/ptfsr2b+1/ptfsr2ab) * resfacZ 
c colour may forbid one of the options above, so if that is the case
c    let us set the corresponding weight to zero
c gluon in-consistent with the colour of b
         if (.not.(icolup(1,ib).eq.icolup(2,ir))) then
            w2 = 0
         endif
         if (.not.(icolup(2,iab).eq.icolup(1,ir))) then
            w3 = 0
         endif
c normalize the weights         
         norm = w1+w2+w3+w4
         w1 = w1/norm
         w2 = w2/norm
         w3 = w3/norm
         w4 = w4/norm         
c let the random number generator decide using the corresponding probabilities
         r = random()
         if(debug) print*,'t: w1,w2,w3,w4,r=',w1,w2,w3,w4,r
         if(r.lt.w1+w2+w3) then
c add both tops
            call add_resonance(-6,3,iwm,iab)
            iwp = iwp+1
            ib = ib+1
            call add_resonance( 6,3,iwp,ib)
            ir = ir + 2
            it = 3
            iat = 4
c and decide on where the radiation comes from
            if(r.lt.w1) then
              if(debug) print*,'t: ttbar ISR'
              n3 = n3 + 1
c gluon from ISR, do nothing               
            elseif(r.lt.w1+w2) then
              if(debug) print*,'t: ttbar g(t)'
c gluon from top
              call add_to_resonance(it,ir)
              n4 = n4 + 1
            else
c gluon from anti-top
              if(debug) print*,'t: ttbar g(t~)'
              call add_to_resonance(iat,ir)
              n5 = n5 + 1
            endif
         else ! w4
c add 23>-24,24 and do nothing else
            if(debug) print*,'t: z ISR'
            call add_resonance(23,3,iwp,iwm)
            n6 = n6 + 1
         endif
      else 
         write(*,*) " guess_resonances: INSANE! Case that is not covered."
         call print_colors
         write(*,*) " Exiting ...!"
         call exit(-1)
      endif

c complete the resonance_colours
      call complete_resonance_colours(iret)
      if(debug) call print_colors
      if(iret.lt.0) then
         write(*,*) " guess_resonances: Impossible to complete resonance colours!"
         call print_colors
         write(*,*) " Exiting ...!"
         call exit(-1)
      endif

      end

      function guess_resonances_std_res_factors()
      integer guess_resonances_std_res_factors
      include 'LesHouches.h'
      integer iwp,iwm,ib,iab,ir
      real * 8 ptop(1:5),patop(1:5),pz(1:5),resfacZ,resfacTTB,r
      real * 8 phepdot,random
      external phepdot,random
c these should be the positions of corresponding particles at this stage      
      iwp = 3
      iwm = 4
      ib = 9
      iab = 10
      ir = 11
c calculate the resonance factor assuming it is Z(>WW),b,b resonance history
      pz = pup(:,iwm) + pup(:,iwp)
      resfacZ = (ph_zmass**2)**2/((phepdot(pz,pz)-ph_zmass**2)**2+(ph_zmass*ph_zwidth)**2)
c calculate the resonance factor assuming it is tt(>WWbb) resonance history
      ptop = pup(:,iwp) + pup(:,ib)
      patop = pup(:,iwm) + pup(:,iab)
      resfacTTB = (ph_tmass**2)**2/((phepdot(ptop,ptop)-ph_tmass**2)**2+(ph_tmass*ph_twidth)**2)*
     1  (ph_tmass**2)**2/((phepdot(patop,patop)-ph_tmass**2)**2+(ph_tmass*ph_twidth)**2)
c let the random number generator decide using the corresponding probabilities
      r = random()
      if(r.lt.resfacZ/(resfacZ+resfacTTB)) then
c Z resonance history
c add 23>-24,24
         call add_resonance(23,3,iwm,iwp)
         guess_resonances_std_res_factors = ir + 1
      else
c top-pair resonance history
         call add_resonance(-6,3,iwm,iab)
c don't forget to shift iwp and ib            
         iwp = iwp+1
         ib = ib+1
         call add_resonance( 6,3,iwp,ib)
         guess_resonances_std_res_factors = ir + 2
      endif
      end

      subroutine add_resonance(idpdg,i0,i1,i2)
      implicit none
      include 'LesHouches.h'
      include 'PhysPars.h'
      real * 8 iduptmp(maxnup),istuptmp(maxnup),mothuptmp(2,maxnup),
     #     icoluptmp(2,maxnup),puptmp(5,maxnup),vtimuptmp(maxnup),
     #     spinuptmp(maxnup)
      integer idpdg,i0,i1,i2,i,j
      if (i0.gt.i1.or.i1.ge.i2) then
         write(*,*) 'wrong sequence in add_resonance'
         write(*,*) 'i0=', i0, 'i1=', i1, 'i2=', i2
         call exit(-1)
         stop
      endif
      do i=i0,nup
         iduptmp(i) = idup(i)
         istuptmp(i) = istup(i)
         mothuptmp(1,i) = mothup(1,i)
         mothuptmp(2,i) = mothup(2,i)
         icoluptmp(1,i) = icolup(1,i)
         icoluptmp(2,i) = icolup(2,i)
         do j=1,5
            puptmp(j,i)=pup(j,i)
         enddo
         vtimuptmp(i)=vtimup(i)
         spinuptmp(i)=spinup(i)
      enddo
      idup(i0)=idpdg
      istup(i0)=+2
      mothup(1,i0)=1
      mothup(2,i0)=2
      icolup(1,i0)=0
      icolup(2,i0)=0
      do j=1,4
         pup(j,i0)=puptmp(j,i1)+puptmp(j,i2)
      enddo
      pup(5,i0)=sqrt(pup(4,i0)**2-
     #     (pup(1,i0)**2+pup(2,i0)**2+pup(3,i0)**2))
      vtimup(i0)=0
      spinup(i0)=9
c     renumber all the other particle mother if their mother is shifting up by one
      do i=i0+1, nup
        if (mothuptmp(1,i).ge.i0) then
          mothuptmp(1,i)=mothuptmp(1,i)+1
          mothuptmp(2,i)=mothuptmp(2,i)+1
        endif
      enddo
c     change mothers of decaying particles
      mothuptmp(1,i1)=i0
      mothuptmp(2,i1)=i0
      mothuptmp(1,i2)=i0
      mothuptmp(2,i2)=i0
      nup=nup+1      
      do i=i0+1, nup
         idup(i) = iduptmp(i-1)
         istup(i) = istuptmp(i-1)
         mothup(1,i) = mothuptmp(1,i-1)
         mothup(2,i) = mothuptmp(2,i-1)
         icolup(1,i) = icoluptmp(1,i-1)
         icolup(2,i) = icoluptmp(2,i-1)
         do j=1,5
            pup(j,i)=puptmp(j,i-1)
         enddo
         vtimup(i)=vtimuptmp(i-1)
         spinup(i)=spinuptmp(i-1)         
      enddo
      end

      subroutine add_to_resonance(i0,ir)
      implicit none
      integer i0, ir
      include 'LesHouches.h'
      real * 8 phepdot
      pup(:,i0) = pup(:,i0) + pup(:,ir)
      pup(5,i0) = sqrt(phepdot(pup(:,i0),pup(:,i0)))
      mothup(:,ir) = i0
      end

      subroutine complete_resonance_colours(iret)
      implicit none
      integer iret
      include 'LesHouches.h'
      integer iup,kup,moth,col(nup),acol(nup),ncol,nacol,jcol,jacol
      logical is_coloured
      external is_coloured
c Look for coloured resonances      
      do iup = 3,nup
         if(istup(iup).gt.1) then
            ncol = 0
            nacol = 0
            do kup = 3,nup
               if(istup(kup).eq.1) then
                  moth = mothup(1,kup)
                  do while (moth.ne.0)
                     if(moth.eq.iup) exit
                     moth = mothup(1,moth)
                  enddo
                  if(moth.ne.0) then
                     if(icolup(1,kup).ne.0) then
                        ncol = ncol + 1
                        col(ncol) = icolup(1,kup)
                     endif
                     if(icolup(2,kup).ne.0) then
                        nacol = nacol + 1
                        acol(nacol) = icolup(2,kup)
                     endif
                  endif
               endif
            enddo
            do jacol = 1,nacol
               do jcol = 1,ncol
                  if(col(jcol).eq.acol(jacol)) then
                     col(jcol) = 0
                     acol(jacol) = 0
                  endif
               enddo
            enddo
c at most one colour and one anticolour should be left
            icolup(:,iup) = 0
            do jcol = 1,ncol
               if(col(jcol).ne.0) then
                  icolup(1,iup) = col(jcol)
                  col(jcol) = 0
                  if(.not.all(col(1:ncol).eq.0)) then
                     iret = -1
                     return
                  else
                     exit
                  endif
               endif
            enddo
            do jacol = 1,nacol
               if(acol(jacol).ne.0) then
                  icolup(2,iup) = acol(jacol)
                  acol(jacol) = 0
                  if(.not.all(acol(1:nacol).eq.0)) then
                     iret = -1
                     return
                  else
                     exit
                  endif
               endif
            enddo
c TODO : check that the colour assignment is consistent with flavour
            if(( idup(iup).eq. 6 .and. .not. (icolup(1,iup).ne.0.and.icolup(2,iup).eq.0) ) .or.
     1         ( idup(iup).eq.-6 .and. .not. (icolup(2,iup).ne.0.and.icolup(1,iup).eq.0) ) .or.
     2         ( abs(idup(iup)) .ne. 6 .and. .not. (icolup(1,iup).eq.0.and.icolup(2,iup).eq.0) ) ) then
               iret = -1
               return
            endif
         endif
      enddo
      iret = 0
      end

      subroutine print_colors
      include 'LesHouches.h'
      integer j
      write (*,'(xA)',advance='no') "t: "
      do j=1,nup-1
        write(*,'(I3AI2AI2I2AI3AI3A)',advance='no') idup(j),'|',j,'(',mothup(1,j),mothup(2,j),'): (',icolup(1,j), ',',icolup(2,j),'), '
      enddo
        write(*,'(I3AI2AI2I2AI3AI3A)')idup(j),'|',j,'(',mothup(1,j),mothup(2,j),'): (',icolup(1,j), ',',icolup(2,j),')'
      end 
