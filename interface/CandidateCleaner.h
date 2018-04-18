#ifndef ANALYSIS_VLQANA_CANDIDATECLEANER_HH
#define ANALYSIS_VLQANA_CANDIDATECLEANER_HH
#include <vector>
#include <algorithm>
#include "AnalysisDataFormats/BoostedObjects/interface/Candidate.h"

///
/// Class to clean one set of candidates T1 from another T2:
/// Configurables:
///   dr_ = Separation in DR(eta, phi) between the two cadidates.
///   ptrel_ = relative pt between the two candidates (used for lepton 2D isolation). 
///   If ptrel_ < 0 (negative), 2D isolation is not applied, and an ordinary DR cleaning is used.
///   If ptrel_ > 0, T1 collection are cleaned from T2 if dr(t1, t2) > dr_ *or* if pt(t1)/pt(t2) > ptrel_
///

class CandidateCleaner {
  public:
    CandidateCleaner(double dr, double ptrel) : dr_(dr), ptrel_(ptrel) {}
    ~CandidateCleaner() {}
    template <class T1, class T2>
      /*    
  void operator () (T1& cleanedcands, T2 othercands) {
      
      // cout<< " pt rel = " << ptrel_ << endl;
      //cout << " cleanedcands size (Before ^^^^^^^^^^^^^) = " << cleanedcands.size()<<endl;
        for (typename  T1::iterator icand = cleanedcands.begin(); icand < cleanedcands.end(); ++icand) {
          bool isclean(true); 
          for (auto othercand : othercands) {
            double dr((icand->getP4()).DeltaR(othercand.getP4()));

            if (ptrel_ < 0) {
              if ( dr < dr_ ) {
		cout<<" %%%% inside Candidate Cleaner %%% " <<endl;
		cout <<"dR = " << dr <<", isclean =" << isclean << endl; 
		
		isclean = false; break ;}
	    }
	      else { //// Apply 2D isolation
		double ptrel ( ( ((icand->getP4()).Vect()).Cross((othercand.getP4()).Vect()) ).Mag()/(othercand.getP4()).Mag() ) ;
		if ( dr < dr_ && ptrel < ptrel_ ) { isclean = false; break ; }  
	      }
	    }
	  if ( !isclean ) cleanedcands.erase(icand) ; 
	  }
	cout<< "%%%% Exiting Candidate Cleaner %%% " <<endl;
	//cout << " cleanedcands size (After ******************) = " << cleanedcands.size()<<endl;
	}
      */
      
                  
      //bug fixed
      
      void operator () (std::vector<T1>& cleanedcands, std::vector<T2> othercands) {
      cleanedcands.erase( std::remove_if(cleanedcands.begin(), cleanedcands.end(),
					 [othercands, this](T1 t) { 
					   // bool isclean(false);
					   bool toremove(false);
					   for (auto cand : othercands) {
					     double dr((t.getP4()).DeltaR(cand.getP4())) ;
					     double ptrel( ( ((t.getP4()).Vect()).Cross((cand.getP4()).Vect()) ).Mag()/(cand.getP4()).Mag() ) ; 
					     // if ( dr > dr_ && ( ptrel_ > 0 ? ptrel > ptrel_ : true)  ) { isclean = true; }
					     // else { isclean = false ; } 
					     // if ( dr < dr_ && ( ptrel_ > 0 ? ptrel < ptrel_ : true)  ) { toremove = true; }
					     // else { toremove = false ; } 

					     if ( dr < dr_ && ( ptrel_ > 0 ? ptrel < ptrel_ : true)  ) { toremove = true; return toremove; }
					     //  cout<<" %%%% inside Candidate Cleaner %%% " <<endl;                                                                         
					     //			     cout <<"dR = " << dr <<", toremove =" << toremove << endl;  
					   }
					   //  return isclean; 
					   return toremove; 

					 }
					 )
			  ,cleanedcands.end()  ) ; 
    }
        
      
    /*  
      //bug reintroduced 

      void operator () (T1& cleanedcands, T2 othercands) {
      for (typename  T1::iterator icand = cleanedcands.begin(); icand < cleanedcands.end(); ++icand) {
	bool isclean(true);
	for (auto othercand : othercands) {
	  double dr((icand->getP4()).DeltaR(othercand.getP4()));
	  if (ptrel_ < 0) {
	    if ( dr < dr_ ) { isclean = false; break ; }
	  }
	  else { //// Apply 2D isolation                                                                                                                   
	    double ptrel ( ( ((icand->getP4()).Vect()).Cross((othercand.getP4()).Vect()) ).Mag()/(othercand.getP4()).Mag() ) ;
	    if ( dr < dr_ && ptrel < ptrel_ ) { isclean = false; break ; }
	  }
	}
	if ( !isclean ) cleanedcands.erase(icand) ;
      }
    }
    */
    
  private:
    double dr_;
    double ptrel_;
}; 
#endif
