// STL include(s):
#include <vector>

// EDM include(s):
#include "MCUtils/PIDUtils.h"
#include "xAODBase/IParticleHelpers.h"

// Local include(s):
#include "NTupAnalyser/Common.h"
#include "NTupAnalyser/TruthHelper.h"

// for documentaiton, see the header file

namespace HC {

  /// print details about the truth particle to the screen
  void printTruthPtcl(const xAOD::TruthParticle *ptcl, TString comment,
                      int childDepth, int parentDepth, int currentDepth)
  {
    // indentation in case we print decay chain. Three spaces per level
    TString indent(Form(Form("%%%ds", 3 * currentDepth), ""));

    if (ptcl == NULL) { printf("%sNULL\n", indent.Data()); return; }

    printf("%sTruth part. ID:%5d, status: %2d, %s  %s\n",
           indent.Data(), ptcl->pdgId(), ptcl->status(), fourVecAsText(ptcl).Data(), comment.Data());

    if (childDepth > 0 || parentDepth > 0) {
      int npar = ptcl->nParents(), nchild = ptcl->nChildren();
      printf("%s-> %d parent and %d children\n", indent.Data(), npar, nchild);

      if (parentDepth > 0)
        for (int ip = 0; ip < npar; ++ip) printTruthPtcl(ptcl->parent(ip), Form("parent %d of ", ip + 1) + comment,
                                                           childDepth - 1, parentDepth - 1, currentDepth + 1);

      if (childDepth > 0)
        for (int ic = 0; ic < nchild; ++ic) printTruthPtcl(ptcl->child(ic), Form("child %d of ", ic + 1) + comment,
                                                             childDepth - 1, parentDepth - 1, currentDepth + 1);
    }
  }

  bool isStable(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("isStable FATAL: particle is NULL"); }

    return ptcl->status() == 1 && ptcl->barcode() < 200000;
  }

  bool isDalitz(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("isDalitz FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    for (auto ptcl : *truthPtcls) // if H -> y*
      if (abs(ptcl->pdgId()) == 25 && (ptcl->status() == 62 || ptcl->status() == 52) && ptcl->nChildren() >= 2 &&
          ((ptcl->child(0) && ptcl->child(0)->pdgId() == 22 && ptcl->child(0)->status() != 1) ||
           (ptcl->child(1) && ptcl->child(1)->pdgId() == 22 && ptcl->child(1)->status() != 1)))
      { return true; }

    return false;
  }

  // Return true if not from hadron
  bool notFromHadron(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("notFromHadron FATAL: particle is NULL"); }

    int ID = ptcl->pdgId();

    // if the particle is a hadron, return false
    if (MCUtils::PID::isHadron(ID)) { return false; }

    // if there are no parents, not from hadron
    if (ptcl->nParents() == 0) { return true; }

    const xAOD::TruthParticle *parent = ptcl->parent(0);
    int parentID = parent->pdgId();

    if (MCUtils::PID::isHadron(parentID)) { return false; } // from hadron!

    if (MCUtils::PID::isTau(parentID) || parentID == ID) { return notFromHadron(parent); }

    // if we get here, all is good
    return true;
  }

  // adds up 4-vectors of all stable particles:
  //  should always give E=m=sqrt(s) \vec{p}=\vec{0} !
  // unless partciels are missing from the file
  TLorentzVector getStableParticle4VectorSum(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getStableParticle4VectorSum FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TLorentzVector sum;

    for (const xAOD::TruthParticle *ptcl : *truthPtcls)
      if (isStable(ptcl)) { sum += ptcl->p4(); }

    return sum;
  }

  double getTruthIsolation(const xAOD::TruthParticle *ptcl,
                           const xAOD::TruthParticleContainer *truthPtcls,
                           double dr,
                           bool chargeOnly,
                           std::vector<int> ignorePdgIds,
                           double ptcut)
  {
    if (truthPtcls == nullptr || ptcl == nullptr) { HC::fatal("getTruthIsolation FATAL: truthPtcls or ptcl is NULL (probably not saved in input sample)"); }

    // Pointer to be used in this function
    const xAOD::TruthParticle *_ptcl = ptcl;

    // Check if this points back to an original particle
    static SG::AuxElement::ConstAccessor<ElementLink<xAOD::IParticleContainer> > acc("originalObjectLink");

    if (acc.isAvailable(*ptcl)) {
      _ptcl = dynamic_cast<const xAOD::TruthParticle *>(xAOD::getOriginalObject(*ptcl));
    }

    // Calculate isolation
    TLorentzVector iso(0, 0, 0, 0);

    for (auto p : *truthPtcls) {
      // Don't count the particles own energy
      if (p->barcode() == _ptcl->barcode())
      { continue; }

      // Only consider stable particles
      if (not isStable(p))
      { continue; }

      // Must be withing the dR cone
      if (HC::DR(p, _ptcl) >= dr)
      { continue; }

      // Only include charge particles?
      if (chargeOnly && p->threeCharge() == 0)
      { continue; }

      // Don't consider muons or neutrinos
      if (std::find(ignorePdgIds.begin(), ignorePdgIds.end(), abs(p->pdgId())) != ignorePdgIds.end())
      { continue; }

      // only consider particles above 1 GeV
      if (ptcut > 0 && p->pt() < ptcut)
      { continue; }

      iso += p->p4();
    }

    if (iso.Px() == 0 && iso.Py() == 0)
    { return 0.0; }

    return iso.Et();
  }

  bool isFinalHiggs(const xAOD::TruthParticle *part)
  {
    if (part == nullptr) { HC::fatal("isFinalHiggs FATAL: particle is NULL"); }

    if (!MCUtils::PID::isHiggs(part->pdgId())) { return false; }

    if (part->child() == nullptr) { return false; }

    if (MCUtils::PID::isHiggs(part->child()->pdgId())) { return false; }

    return true;
  }

  TruthPtcls getFinalHiggsBosons(const xAOD::TruthParticleContainer *truthParticles)
  {
    if (truthParticles == nullptr) { HC::fatal("getFinalHiggsBosons FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls higgs(SG::VIEW_ELEMENTS);

    for (auto part : *truthParticles) {
      if (isFinalHiggs(part))
      { higgs.push_back(part); }
    }

    return higgs;
  }

  bool isFromHiggs(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("isFromHiggs FATAL: particle is NULL"); }

    if (MCUtils::PID::isHiggs(ptcl->pdgId())) { return true; }

    if (ptcl->parent() == nullptr) { return false; }

    return isFromHiggs(ptcl->parent());
  }

  bool isFromZ(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("isFromZ FATAL: particle is NULL"); }

    if (MCUtils::PID::isZ(ptcl->pdgId())) { return true; }

    if (ptcl->parent() == nullptr) { return false; }

    return isFromZ(ptcl->parent());
  }

  bool isFromBhadron(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("isFromBhadron FATAL: particle is NULL"); }

    if (MCUtils::PID::isBottomHadron(ptcl->pdgId())) { return true; }

    if (ptcl->parent() == nullptr) { return false; }

    return isFromBhadron(ptcl->parent());
  }

  bool isGoodTruthPhoton(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("isGoodTruthPhoton FATAL: particle is NULL"); }

    return isStable(ptcl) && MCUtils::PID::isPhoton(ptcl->pdgId()) && notFromHadron(ptcl);
  }

  bool isGoodTruthElectron(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("isGoodTruthElectron FATAL: particle is NULL"); }

    return isStable(ptcl) && MCUtils::PID::isElectron(ptcl->pdgId()) && notFromHadron(ptcl);
  }

  bool isGoodTruthMuon(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("isGoodTruthMuon FATAL: particle is NULL"); }

    return isStable(ptcl) && MCUtils::PID::isMuon(ptcl->pdgId()) && notFromHadron(ptcl);
  }

  bool isZdecayLepton(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("isZdecayLepton FATAL: particle is NULL"); }

    return isStable(ptcl) && (MCUtils::PID::isElectron(ptcl->pdgId()) || MCUtils::PID::isMuon(ptcl->pdgId())) && isFromZ(ptcl);
  }

  std::vector<const xAOD::TruthParticle *> getGoodTruthPhotonsOld(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getGootTruthPhotonsOld FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    ::std::vector<const xAOD::TruthParticle *> truthPhotons;

    for (const xAOD::TruthParticle *ptcl : *truthPtcls)
      if (isGoodTruthPhoton(ptcl)) { truthPhotons.push_back(ptcl); }

    return truthPhotons;
  }

  //! /brief returns all stable photons that do not originate from hadrons
  TruthPtcls getGoodTruthPhotons(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getGootTruthPhotons FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls ys(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls) if (isGoodTruthPhoton(ptcl)) { ys.push_back(ptcl); }

    return ys;
  }

  //! /brief returns all stable electrons that do not originate from hadrons
  TruthPtcls getGoodTruthElectrons(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getGootTruthElectrons FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls es(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls) if (isGoodTruthElectron(ptcl)) { es.push_back(ptcl); }

    return es;
  }

  //! /brief returns stable DM particles (PDGID==1000022)
  TruthPtcls getGoodTruthDM(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getGootTruthDM FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls es(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls) {
      if (ptcl->status() != 1) { continue; }

      if (ptcl->absPdgId() != 1000022) { continue; }

      es.push_back(ptcl);
    }

    return es;
  }

  //! /brief returns all stable muons that do not originate from hadrons
  TruthPtcls getGoodTruthMuons(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getGootTruthMuons FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls mus(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls) if (isGoodTruthMuon(ptcl)) { mus.push_back(ptcl); }

    return mus;
  }

  //! /brief returns all stable electrons or muons originating from a Z boson
  TruthPtcls getZdecayLeptons(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getZdecayLeptons FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls ls(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls) if (isZdecayLepton(ptcl)) { ls.push_back(ptcl); }

    return ls;
  }

  TruthPtcls getHadronsAndTheirDecay(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getHadronsAndTheirDecay FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls hadrons(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls) {
      if (!isStable(ptcl)) { continue; }

      if (isGoodTruthPhoton(ptcl)) { continue; }

      if (isGoodTruthElectron(ptcl)) { continue; }

      if (isGoodTruthMuon(ptcl)) { continue; }

      hadrons.push_back(ptcl);
    }

    return hadrons;
  }

  TruthPtcls getStableDecayProducts(const xAOD::TruthParticle *ptcl)
  {
    if (ptcl == nullptr) { HC::fatal("getStableDecayProducts FATAL: particle is NULL"); }

    TruthPtcls decay(SG::VIEW_ELEMENTS);

    if (HC::isStable(ptcl)) { decay.push_back(ptcl); return decay; }

    for (size_t ichild = 0; ichild < ptcl->nChildren(); ++ichild)
      if (ptcl->child(ichild)) for (auto p : getStableDecayProducts(ptcl->child(ichild))) { decay.push_back(p); }

    return decay;
  }


  //! /brief returns all stable electrons that do not originate from hadrons
  TruthPtcls getBHadrons(const xAOD::TruthParticleContainer *truthPtcls, double pTcut)
  {
    if (truthPtcls == nullptr) { HC::fatal("getBHadrons FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls Bs(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls)
      if (MCUtils::PID::isBottomHadron(ptcl->pdgId()) && (pTcut < 0 || ptcl->pt() > pTcut))
      { Bs.push_back(ptcl); }

    return Bs;
  }

  TruthPtcls getDHadrons(const xAOD::TruthParticleContainer *truthPtcls, double pTcut)
  {
    if (truthPtcls == nullptr) { HC::fatal("getDHadrons FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls Ds(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls)
      if (MCUtils::PID::isCharmHadron(ptcl->pdgId()) && (pTcut < 0 || ptcl->pt() > pTcut))
      { Ds.push_back(ptcl); }

    return Ds;
  }

  TruthPtcls getPhotonsFromHiggs(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getPhotonsFromHiggs FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls ys(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls)
      if (isGoodTruthPhoton(ptcl) && isFromHiggs(ptcl)) { ys.push_back(ptcl); }

    return ys;
  }

  TruthPtcls getHiggsDecayProducts(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getHiggsDecayProducts FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls decay(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls)
      if (isStable(ptcl) && isFromHiggs(ptcl)) { decay.push_back(ptcl); }

    return decay;
  }


  TruthPtcls getMuonsFromBs(const xAOD::TruthParticleContainer *truthPtcls)
  {
    if (truthPtcls == nullptr) { HC::fatal("getMuonsFromBs FATAL: truthPtcls is NULL (probably not saved in input sample)"); }

    TruthPtcls mus(SG::VIEW_ELEMENTS);

    for (auto ptcl : *truthPtcls)
      if (isStable(ptcl) && MCUtils::PID::isMuon(ptcl->pdgId()) && isFromBhadron(ptcl))
      { mus.push_back(ptcl); }

    return mus;
  }

  TruthParticleStruct identifyTruthParticles(xAOD::TEvent *event,
                                             double /*jet_pTcut*/)
  {

    const xAOD::TruthParticleContainer *truthParticles = nullptr;

    if (event->retrieve(truthParticles, "TruthParticle").isFailure())
    { HC::fatal("Cannot access TruthParticle"); }

    const xAOD::JetContainer *truthJets = nullptr;

    if (event->retrieve(truthJets, "AntiKt4TruthJets").isFailure())
    { HC::fatal("Cannot access AntiKt4TruthJets"); }

    return HC::identifyTruthParticles(truthParticles, truthJets);
  }

  void removeTruthOverlap(DataVector<xAOD::IParticle> &photons,
                          DataVector<xAOD::IParticle> &electrons,
                          DataVector<xAOD::IParticle> &/*muons*/,
                          DataVector<xAOD::IParticle> &jets,
                          double jet_pTcut)
  {
    for (auto tj = jets.begin(); tj != jets.end();) {
      // apply a pT cut, if requested
      if (jet_pTcut > 0 && (*tj)->pt() < jet_pTcut) {
        tj = jets.erase(tj);
        continue;
      }

      // ignore jets overlapping with good photons, electrons or muons
      if (HC::minDRrap(*tj, photons) < 0.4) {
        tj = jets.erase(tj);
        continue;
      }

      if (HC::minDRrap(*tj, electrons) < 0.4) {
        tj = jets.erase(tj);
        continue; // <<== WZ jets should not do this
      }

      // if (HC::minDRrap(tj,    muonss)<0.4) continue; // ??
    }
  }

  TruthParticleStruct identifyTruthParticles(const xAOD::TruthParticleContainer *truthPtcls,
                                             const xAOD::JetContainer *truthJets,
                                             double jet_pTcut)
  {
    TruthParticleStruct tp;
    TruthPtcls ys  = getGoodTruthPhotons(truthPtcls);
    TruthPtcls es  = getGoodTruthElectrons(truthPtcls);
    TruthPtcls mus = getGoodTruthMuons(truthPtcls);
    TruthPtcls hads =  getHadronsAndTheirDecay(truthPtcls);

    // TO-DO
    // Dressing should happen here !
    tp.electrons = es;
    tp.muons = mus;
    // some ys should probably go to hads here
    tp.photons = ys;
    tp.hadrons = hads;

    tp.photonsFromHiggs = getPhotonsFromHiggs(truthPtcls);
    // this one might be slow ... ?
    tp.HiggsDecay = getHiggsDecayProducts(truthPtcls);

    tp.Bhadrons    = getBHadrons(truthPtcls);
    tp.Dhadrons    = getDHadrons(truthPtcls);
    tp.muonsFromBs = getMuonsFromBs(truthPtcls);

    // TruthJets jets(truthJets->begin(), truthJets->end(), SG::VIEW_ELEMENTS);
    TruthJets jets(SG::VIEW_ELEMENTS);
    TruthJets bjets(SG::VIEW_ELEMENTS);
    TruthJets cjets(SG::VIEW_ELEMENTS);
    TruthJets lightJets(SG::VIEW_ELEMENTS);

    // Here applying a 5 GeV cut for the jet labelling
    TruthPtcls Bs    = getBHadrons(truthPtcls, 5.0 * HC::GeV);
    TruthPtcls Ds    = getDHadrons(truthPtcls, 5.0 * HC::GeV);

    // removeTruthOverlap(tp.photons, tp.electrons, tp.muons, jets, jet_pTcut);

    for (const xAOD::Jet *tjet : *truthJets) {
      // apply a pT cut, if requested
      if (jet_pTcut > 0 && tjet->pt() < jet_pTcut) { continue; }

      // ignore jets overlapping with good photons, electrons or muons
      if (HC::minDRrap(tjet, ys) < 0.4) { continue; }

      if (HC::minDRrap(tjet, es) < 0.4) { continue; } // <<== WZ jets should not do this

      // if (HC::minDRrap(tjet,mus)<0.4) continue; ??
      jets.push_back(tjet);

      // classify all jets into b, c or light
      if (HC::minDRrap(tjet, Bs) < 0.4) { bjets.push_back(tjet); }
      else if (HC::minDRrap(tjet, Ds) < 0.4) { cjets.push_back(tjet); }
      else { lightJets.push_back(tjet); }
    }

    // later: further split light jets into: LQ, gluon, unmatched
    tp.jets  = jets;
    tp.bJets = bjets;
    tp.cJets = cjets;
    tp.lightJets = lightJets;
    return tp;
  }

  void printTruthParticles(const TruthParticleStruct &tp)
  {
    printf("Identified truth particles:\n");
    printf("  %lu photons\n", tp.photons.size());
    printf("  %lu electrons, %lu muons\n", tp.electrons.size(), tp.muons.size());
    printf("  %lu photons from Higgs\n", tp.photonsFromHiggs.size());
    printf("  %lu B- and %lu D-hadrons\n", tp.Bhadrons.size(), tp.Dhadrons.size());
    printf("  %lu muons from B-hadrons\n", tp.muonsFromBs.size());
    printf("  %lu jets, of which\n", tp.jets.size());
    printf("  %lu b-, %lu c- and %lu light jets\n",
           tp.bJets.size(), tp.cJets.size(), tp.lightJets.size());
  }

} // namespace HC
