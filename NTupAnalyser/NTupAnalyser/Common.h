#pragma once

// EDM include(s):
#ifndef __CINT__
  #include "PathResolver/PathResolver.h"
  #include "xAODBase/IParticle.h"
  #include "xAODRootAccess/tools/Message.h"
#endif

// ROOT include(s):
#include "TH1.h"
#include "TLorentzVector.h"
#include "TObject.h"
#include "AthContainers/ConstDataVector.h"
#include "AthContainers/DataVector.h"
#include "xAODJet/JetContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticleAuxContainer.h"


//! \name   A few general helper methods and definitions
//@{


//! \brief Hcgamma+gamma namespace
namespace HC {


  typedef const xAOD::TruthParticleContainer            TruthParticles;
  typedef ConstDataVector<xAOD::TruthParticleContainer> TruthContainer;
  // constants and typedefs

  /// Helper macro to check xAOD::TReturnCode return values
  /// See Attila's slide 8-9 at: https://indico.cern.ch/event/362819/
  /// TODO: replace with ANA_CHECK (JR)
#define EL_CHECK( COMMENT, EXP )                        \
  do {                                                  \
    if ( ! EXP.isSuccess() ) {                          \
      Error( COMMENT,                                   \
             XAOD_MESSAGE("\n  Failed to execute %s"),   \
             #EXP );                                     \
      return EL::StatusCode::FAILURE;                   \
    }                                                   \
  } while( false );

#define CP_CHECK( COMMENT, EXP )                        \
  do {                                                  \
    if ( EXP != CP::SystematicCode::Ok ) {              \
      Error( COMMENT,                                   \
             XAOD_MESSAGE("\n  Failed to execute %s"),   \
             #EXP );                                     \
      return CP::SystematicCode::Unsupported;           \
    }                                                   \
  } while( false );

#define CC_CHECK( COMMENT, EXP )                        \
  do {                                                  \
    if ( EXP != CP::CorrectionCode::Ok ) {              \
      Error( COMMENT,                                   \
             XAOD_MESSAGE("\n  Failed to execute %s"),   \
             #EXP );                                     \
    }                                                   \
  } while( false );

#define HG_CHECK( COMMENT, EXP )                        \
  do {                                                  \
    if ( ! EXP.isSuccess() ) {                          \
      Error( COMMENT,                                   \
             XAOD_MESSAGE("\n  Failed to execute %s"),   \
             #EXP );                                     \
    }                                                   \
  } while( false );

  //! \brief method to abort program with error message
  void fatal(TString msg);

  //! \brief typedef for a vector of doubles (to save some typing)
  typedef std::vector<double>  NumV;

  //! \brief typedef for a vector of ints (to save some typing)
  typedef std::vector<int>     IntV;

  //! \brief typedef for a vector of TStrings (to save some typing)
  typedef std::vector<TString> StrV;

  //! \brief typedef for a vector of std::strings (to save some typing)
  typedef std::vector<std::string> StdStrV;

  //! \brief Converts a text line to a vector of words
  //  \param str input string with words
  //  \param sep separator to define where a word ends or starts
  StrV vectorize(TString str, TString sep = " ");

  //! \brief print 4-vector as string
  TString fourVecAsText(const TLorentzVector &p4);
  TString fourVecAsText(const xAOD::IParticle *p);





  template <typename T> // what is this used for? JR
  struct Identity {
    typedef T type;
  };

  //1*GeV = 1000*MeV
  static const double GeV(1000), invGeV(1.0 / GeV);



  //! \brief calculates DeltaR in (y,phi)-space instead of (eta,phi) given by p4().DeltaR()
  inline double DRrap(const TLorentzVector &p1, const TLorentzVector &p2)
  {
    double dy(p1.Rapidity() - p2.Rapidity()), dphi(p1.DeltaPhi(p2));
    return sqrt(dy * dy + dphi * dphi);
  }

  //! \brief returns true if a given file or directory exist
  bool fileExist(TString fn);

  TH1 *getHistogramFromFile(TString fname, TString hname);

  template<typename T>
  std::unique_ptr<T> getHistogramPtrFromFile(TString fname, TString hname)
  {
    std::unique_ptr<T> ptr(dynamic_cast<T *>(getHistogramFromFile(fname, hname)));
    return ptr;
  }




#ifndef __MAKECINT__


  //! \brief calculates DeltaR in (y,phi)-space instead of (eta,phi) given by p4().DeltaR()
  inline double DRrap(const xAOD::IParticle *p1, const xAOD::IParticle *p2)
  {
    return DRrap(p1->p4(), p2->p4());
  }

  //! \brief calculates DeltaR in (eta,phi)-space
  inline double DR(const xAOD::IParticle *p1, const xAOD::IParticle *p2)
  {
    return p1->p4().DeltaR(p2->p4());
  }

  //! \brief returns smallest DR between ptcl and any of the objects in ptcls
  //!        if ptcl occurs in the list of particles, it is ignored
  template <class T> double minDR(const xAOD::IParticle *ptcl, T ptcls)
  {
    double mindr = 99;

    for (auto p : ptcls) if (p != ptcl && DR(ptcl, p) < mindr) { mindr = DR(ptcl, p); }

    return mindr;
  }

  //! \brief returns smallest DR between ptcl and any of the objects in ptcls in (y,phi)-space
  //!        if ptcl occurs in the list of particles, it is ignored
  template <class T> double minDRrap(const xAOD::IParticle *ptcl, T ptcls)
  {
    double mindr = 99;

    for (auto p : ptcls) if (p != ptcl && DRrap(ptcl, p) < mindr) { mindr = DRrap(ptcl, p); }

    return mindr;
  }

  //! \brief check if the given object is inside a given container. Returns a bool.
  template <class T1, class T2> bool isObjInCont(T1 &ptclSearch, T2 &ptcls)
  {
    bool found = false;

    for (auto ptcl : ptcls)
      if (ptcl == ptclSearch)
      { found = true; }

    return found;
  }

#endif

  //@}
}
