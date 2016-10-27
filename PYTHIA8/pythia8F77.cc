// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2012 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk. 
// It studies the charged multiplicity distribution at the LHC.

//#include "Pythia8/Pythia.h"
//#include "Pythia8/LHAFortran.h"

#include "Pythia.h"
#include "LHAFortran.h"

using namespace Pythia8; 

// extern "C" {
//   extern struct {
//     double vetoscale;
//     int canveto; // 1: can veto, #1 cannot
//   } resonancevetos_;
// }

// class MyUserHooks : public UserHooks {
  
// public:
  
//   MyUserHooks() {cout << "Setting up Hook\n";}
//   // Destructor
//   ~MyUserHooks() {}

//   // only used to print out and check that scalup is effectively changed
//   // on an event-by-event basis
//   // bool canVetoMPIEmission() { return true ; }
//   // bool doVetoMPIEmission(int, const Event &e) {
//   //   cout<<"here"<<endl;
//   //   cout<<infoPtr->scalup()<<endl;
//   //   return false;
//   // }

//   //   // Allow process cross section to be modified..
//   //   virtual bool canSetResonanceScale() {
//   //     if(resonancevetos_.canveto == 1) return true;
//   //     else return false; }
  
//   //   virtual double scaleResonance( const int iRes, const Event& event) {
//   //     //    cout << "now in scaleResonance"; cout << resonancevetos_.vetoscale;
//   //     //    cout << "\n";
//   //     return  resonancevetos_.vetoscale;
//   //   }    

// private:

// };

// uncommet to check that scalup is changed on an event-by-event basis.
//MyUserHooks* MyHook=new MyUserHooks();

Pythia pythia;
LHAupFortran LHAinstance;

extern "C" {
  // F77 interface to pythia8
  void pythia_init_() {
    // Generator. Process selection. LHC initialization. Histogram.

    // uncommet to check that scalup is changed on an event-by-event basis, if
    // the changescalup flag is activated in the input file.
    //pythia.setUserHooksPtr(MyHook); 
    
    pythia.readString("Tune:pp = 14"); // Monash2013 tune
    //    pythia.readString("Tune:pp = 5"); // Default

    pythia.readString("SpaceShower:QEDshowerByQ = off"); // From quarks.        
    pythia.readString("SpaceShower:QEDshowerByL = off"); // From Leptons.       
    pythia.readString("TimeShower:QEDshowerByQ = off"); // From quarks.         
    pythia.readString("TimeShower:QEDshowerByL = off"); // From Leptons.

    pythia.readString("PartonLevel:MPI = off"); 
    //pythia.readString("PartonLevel:ISR = off"); 
    //pythia.readString("PartonLevel:FSR = off"); 
    pythia.readString("HadronLevel:all = off"); 
    pythia.readString("BeamRemnants:primordialKT = off");


    //make pi0 stable
    pythia.readString("111:mayDecay = off");

    // ER: leave the following commented
    //pythia.readString("SpaceShower:pTmaxMatch = 0");
    //pythia.readString("TimeShower:pTmaxMatch = 0");

    pythia.init(&LHAinstance);
  }

  void pythia_next_(int & iret){
    // Begin event loop. Generate event. Skip if error. List first one.
    iret = pythia.next();
    //    cout<<" EVENT SIZE IS "<<pythia.event.size()<<" "<<iret<<endl;
  }

  void pythia_to_hepevt_(const int &nmxhep, int & nhep, int * isthep,
			 int * idhep,
			 int  (*jmohep)[2], int (*jdahep)[2],
			 double (*phep)[5], double (*vhep)[4]) {
    nhep = pythia.event.size();
    if(nhep>nmxhep) {cout << "too many particles!" ; exit(-1); }
    for (int i = 0; i < pythia.event.size(); ++i) {
      *(isthep+i) = pythia.event.statusHepMC(i);
      *(idhep+i) = pythia.event[i].id();
      // All pointers should be increased by 1, since here we follow
      // the c/c++ convention of indeces from 0 to n-1, but i fortran
      // they are from 1 to n.
      (*(jmohep+i))[0] = pythia.event[i].mother1() + 1;
      (*(jmohep+i))[1] = pythia.event[i].mother2() + 1;
      (*(jdahep+i))[0] = pythia.event[i].daughter1() + 1;
      (*(jdahep+i))[1] = pythia.event[i].daughter2() + 1;
      (*(phep+i))[0] = pythia.event[i].px();
      (*(phep+i))[1] = pythia.event[i].py();
      (*(phep+i))[2] = pythia.event[i].pz();
      (*(phep+i))[3] = pythia.event[i].e();
      (*(phep+i))[4] = pythia.event[i].m();
    }
    // override mother of very first event, set to 0
    *(jmohep)[0] = 0 ;
    *(jmohep)[1] = 0 ;
  }

  void pythia_stat_() {
    pythia.stat(); 
  }


}





// User-written routine that does the event generation and fills hepeup.
// We do nothing here; we fill the common block on the fortran side
// before we start pythia and before we call each event.

bool LHAupFortran::fillHepRup() {
  return true;
}

bool  LHAupFortran::fillHepEup(){
  return true;
}

