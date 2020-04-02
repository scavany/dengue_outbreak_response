#include "Location.h"
#include <sstream>
#include <iostream>
#include <cmath>

string Location::getLocID() const {
  return locID;
}

string Location::getLocType() const {
  return locType;
}

string Location::getMoHID() const{
  return MoHID;
}

double Location::getLocX() const {
  return xCor;
}

double Location::getLocY() const {
  return yCor;
}

double Location::getEmergenceRate() const {
  return emergenceRate;
}

double Location::getEggs() const {
  return eggs;
}

double Location::getLarvae() const {
  return larvae;
}

double Location::getPupae() const {
  return pupae;
}

double Location::getLarvalCapacity() const {
  return larvalCapacity;
}

double Location::getInitialAdults() const {
  return initialAdults;
}

std::deque<unsigned> Location::getRecentAdults() const {
  return recentAdults;
}

std::deque<double> Location::getRecentDensityDependenceTerms() const {
  return recentDensityDependenceTerms;
}


Location::Location(string lID, string lType, string mID, double x, double y, double e) {
  locID = lID;
  locType = lType;
  xCor = x;
  yCor = y;
  MoHID = mID;
  emergenceRate = e;
  closeLocs.reset(new vector<string>());
  radiusLocs.reset(new vector<string>());
  infectedVisitor = false;
  visitorBites.clear();
  bitesCounterEnabled = false;
  insecticideDecayRate = 1.0;
  insecticideEfficacy = 0.0;
  insecticideStartDay = 0;
  insecticideSprayed = false;
}


// For delay differential equations model, with global density dependence:
Location::Location(string lID, string lType, double x, double y, std::deque<unsigned> adults, double a) {
  locID = lID;
  locType = lType;
  xCor = x;
  yCor = y;
  initialAdults = a;
  recentAdults = adults;
  closeLocs.reset(new vector<string>());
  radiusLocs.reset(new vector<string>());
  infectedVisitor = false;
  visitorBites.clear();
  bitesCounterEnabled = false;
  insecticideDecayRate = 1.0;
  insecticideEfficacy = 0.0;
  insecticideStartDay = 0;
  insecticideSprayed = false;
}

//extra mortality, extra locations
Location::Location(string lID, string lType, string mID, double x, double y, double e, double l, double p, double kl) {
  locID = lID;
  locType = lType;
  MoHID = mID;
  xCor = x;
  yCor = y;
  pupae = p;
  larvae = l;
  eggs = e;
  larvalCapacity = kl;
  closeLocs.reset(new vector<string>());
  radiusLocs.reset(new vector<string>());
  infectedVisitor = false;
  visitorBites.clear();
  bitesCounterEnabled = false;
  insecticideDecayRate = 1.0;
  insecticideEfficacy = 0.0;
  insecticideStartDay = 0;
  insecticideSprayed = false;
}


string Location::getRandomCloseLoc(RandomNumGenerator& rGen) {
  int i = closeLocs->size();
  if (i > 0)
  return (*closeLocs)[rGen.getMozNextLoc(i)];
  else return "TOO_FAR_FROM_ANYWHERE";
}

double Location::getDistanceFromLoc(Location& loc) const {
  return sqrt((xCor - loc.getLocX()) * (xCor - loc.getLocX()) + (yCor - loc.getLocY()) * (yCor - loc.getLocY()));
}

vector<string> Location::getRadiusLocations(){
  return(*(radiusLocs.get()));
}

void Location::addCloseLoc(string loc) {
  closeLocs->push_back(loc);
}

void Location::addRadiusLoc(string loc) {
  radiusLocs->push_back(loc);
}

void Location::sprayAdultInsecticide(unsigned currDay,double efficacy_,double residuality_, unsigned efficacyLength, unsigned residualityLag){
  insecticideStartDay = currDay;
  insecticideEndDay = currDay + efficacyLength;
  insecticideResidualityStartDay = currDay + residualityLag;
  insecticideSprayed = true;
  insecticideEfficacy = efficacy_;
  insecticideDecayRate = residuality_;
}

double Location::getIncreasedMortalityInsecticide(unsigned currDay, double mozMortality_){
  if(insecticideSprayed){
    double mortality_rate = -log(1 - mozMortality_);
    double decayProportion = currDay > insecticideResidualityStartDay ? exp(-insecticideDecayRate * (currDay - insecticideResidualityStartDay)) : 1;
    mortality_rate += insecticideEfficacy * decayProportion;
    if (currDay > insecticideEndDay) {
      insecticideSprayed = false;
    }
    return(1-exp(-mortality_rate));
  }else{
    return(mozMortality_);
  }
}

void Location::addHuman(sp_human_t h) {
  humans.insert(h);
}

void Location::removeHuman(sp_human_t h){
  // set.erase(val) requires no check
  humans.erase(h);
}

//Location::Location() {
//}

//Location::Location(const Location& orig) {
//}

//Location::~Location() {
//}

void Location::updatePupae(double flux){
  pupae = flux;
}

void Location::updateLarvae(double flux){
  larvae = flux;
}

void Location::updateEggs(double flux){
  eggs = flux;
}

void Location::updateRecentAdults(unsigned newAdults){
  recentAdults.push_back(newAdults);
  recentAdults.pop_front();
}

void Location::updateRecentDensityDependenceTerms(double newTerm){
  recentDensityDependenceTerms.push_back(newTerm);
  recentDensityDependenceTerms.pop_front();
}


void Location::updateInfectedVisitor(){
  infectedVisitor = false;
  for(auto itHum = humans.begin(); itHum != humans.end(); itHum++){
    if((*itHum)->infection != nullptr){
      infectedVisitor = true;
      return;
    }
  }
}

void Location::printHumans(){
  for(auto itHum = humans.begin(); itHum != humans.end(); itHum++){
    printf("Human %s in Location %s\n", (*itHum)->getPersonID().c_str(), locID.c_str());
  }
}



void Location::increaseBites(string personID){
  if(bitesCounterEnabled == true){
    visitorBites[personID]++;
  }
}

double Location::calculateGiniIndex(){
  double num = 0;
  double den = 0;
  for(auto it = visitorBites.begin(); it != visitorBites.end();){
    den += (*it).second;
    for(auto jt = visitorBites.begin(); jt != visitorBites.end();){
      num += std::abs((*it).second - (*jt).second);
      ++jt;
    }
    ++it;
  }
  if(den > 0){
    return(num / (2.0 * visitorBites.size() * den));
  }else{
    return(-1.0);
  }
}
