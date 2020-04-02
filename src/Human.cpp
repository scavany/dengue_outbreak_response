#include "Human.h"
#include <sstream>

Human::Human(string hID,
	     int hMemID,
	     string zID,
	     char gen,
	     int birthYear,
	     int deathYear,
	     RandomNumGenerator& rGen)
{
  houseID = hID;
  infLocID = "";
  zoneID = zID;
  houseMemNum = hMemID;
  personID = hID + std::to_string(hMemID);
  visitorID = "";
  gender = gen;
  bday = birthYear * 365 + rGen.getRandomNum(365);
  if(deathYear >= 0){
    dday = deathYear * 365 + rGen.getRandomNum(365);
    if(dday <= bday){
      dday = bday + 1; // give them a chance to live at least one day
    }
  }else{
    dday = -1;
  }
  infection.reset(nullptr);
  vaccinated = false;
  dead = false;
  resetRecent();
  cohort = 0;
  tAge = 0;
  vday = -1;
  trialDay = 0;
  vaxWaning_pos = 0;
  vaxWaning_neg = 0;
  vaccineComplete = false;
  enrolledInTrial = false;
  seroStatusAtVaccination = false;
  dateOfExposures.clear();
  bites = 0;
  for(int i = 0;i < 4; i++){
    preExposureAtVaccination[i] = 0;
    exposedCount[i] = 0;
    dateOfExposures.push_back("");
    R0[i+1] = -1;
  }

  infected = false;
  infectedImport = false;
  symptomatic = false;
  hospitalized = false;
  vaccineImmunity = false;
  reportSymptoms = false;
  vaccineDosesReceived = 0;
  lastDayContactedByTrial = 0;
  firstContactWithTrial = 0;
  selfReportProb = 0.0;
  immunity_temp = false;
  trajectories.reset(nullptr);
  for(int i = 1; i < N_SERO + 1; i++){
    updateImmunityPerm(i,false);
  }
  trajDay = 0;
  immEndDay = 0;
}

void Human::initializeHuman(unsigned currDay, vector<double> FOI, RandomNumGenerator& rGen){
  initiateBodySize(currDay,rGen);

  attractiveness = rGen.getAttractiveness();
  if(attractiveness == -1){
    updateAttractiveness(currDay);
  }
  immunity_temp = false;
  if(currDay == 0){
    // Set the initial conditions for the immune profile by serotype
    for(int i = 0; i < N_SERO; i++){
      updateImmunityPerm(i+1, rGen.getHumanSeropositivity(FOI[i], double(currDay - bday)));
    }
  }
}

void Human::setTrajectories(unique_ptr<traject_t> & paths){
  if(!trajectories){
    trajectories = std::move(paths);
    trajDay = 0;
  }else{
    printf("Trajectories for %s already set\n", personID.c_str());
  }
}

void Human::checkRecovered(int currDay){
  if(infection->getEndDay() <= currDay){
    infection.reset(nullptr);
    infected = false;
    hospitalized = false;
    symptomatic = false;
    reportSymptoms = false;
    if(infectedImport){
      infectedImport = false;
      visitorID = "";
    }
  }
}

double Human::getAttractiveness() const {
  return attractiveness;
}

int Human::getAgeDays(unsigned currDay) const {
  return currDay - bday;
}

double Human::getBodySize() const {
  return bodySize;
}

const string & Human::getCurrentLoc(double time){
  const auto & today = (*trajectories)[trajDay];
  auto itrLoc = today.begin();

  for(double runSum = 0; runSum < time && itrLoc != today.end(); itrLoc++){
    runSum += itrLoc->second;
  }
  if(today.size() == 1){
    itrLoc = today.begin();
  } else {
    itrLoc--;
  }

  return itrLoc->first;
}

char Human::getGender() const {
  return gender;
}

unsigned Human::getImmStartDay() const {
  return immStartDay;
}

unsigned Human::getVaxImmStartDay() const {
  return vaxImmStartDay;
}

unsigned Human::getImmEndDay() const {
  return immEndDay;
}

unsigned Human::getVaxImmEndDay() const {
  return vaxImmEndDay;
}

const string & Human::getZoneID() const{
  return zoneID;
}

const string & Human::getHouseID() const {
  return houseID;
}

const string & Human::getInfLocID() const {
  return infLocID;
}

const string & Human::getPersonID() const{
  return personID;
}

const string & Human::getVisitorID() const{
  return visitorID;
}

int Human::getHouseMemNum() const {
  return houseMemNum;
}

const std::set<string> & Human::getLocsVisited(){
  locsVisited.clear();
  for(auto & the_traj : (*trajectories)){
    for(auto & place : the_traj) {
      // set.insert includes check
      locsVisited.insert(place.first);
    }
  }
  return locsVisited;
}



int Human::getPreviousInfections(){
  int previnf = 0;

  if(immunity_perm[1])
    previnf++;
  if(immunity_perm[2])
    previnf++;
  if(immunity_perm[3])
    previnf++;
  if(immunity_perm[4])
    previnf++;

  return previnf;
}

//traject_t const& Human::getTrajectory(unsigned i) const {
//   return (*trajectories.get())[i];
//}

void Human::infectImport(
			 int currentDay,
			 unsigned infectionType,
			 RandomNumGenerator * rGen,
			 int visit_)
{
  // Imported cases have an infection but do not report any symptoms... They should be invisible for our surveillance system
  infectedImport = true;
  visitorR0[infectionType] = 0;
  infection.reset(new Infection(
				currentDay + 1, currentDay + 15, 0.0, infectionType, getPreviousInfections() == 0, recent_dis > 0, exp(rGen->getRandomNormal() * 0.2701716 + 1.750673)));
  visitorID = "IQVIS" + std::to_string(visit_);
}

void Human::infect(
		   int currentDay,
		   unsigned infectionType,
		   RandomNumGenerator * rGen,
		   map<unsigned,double> * disRates,
		   map<unsigned,double> * hospRates,
		   sp_human_t humInf,
		   string infLocID_)
{
  double RR = 1.0;
  double RRInf = 1.0;
  double RRDis = 1.0;
  double RRHosp = 1.0;

  int vaxAdvancement = 0;

  if(vaccinated){
    //printf("Infecting human vaccinated with profile %s...age: %d\n",vaccineProfile.getMode().c_str(),getAgeDays(currentDay));
    if(vaccineProfile.getVaccineID() == -1){
      printf("Human::infect vaccineProfile is -1 and person is vaccinated\n");
      exit(1);
    }

    // There are multiple vaccines supported: GSK, advance (Sanofi), or age
    // We have to specify the effects for each of these vaccine modes
    if(vaccineProfile.getMode() == "advance"){
      vaxAdvancement = 1;
    }else if(vaccineProfile.getMode() == "age"){
      RR =  vaccineProfile.getRR(getPreviousInfections(), double(getAgeDays(currentDay)));
      RRInf = pow(RR, vaccineProfile.getPropInf());
      RRDis = pow(RR, 1.0 - vaccineProfile.getPropInf());
      //	    printf("Human ID %s vaccinated age %d. RR: %.2f RRInf: %.2f RRDis: %.2f\n",getHouseID().c_str(), getAgeDays(currentDay),RR,RRInf,RRDis);
    }else if(vaccineProfile.getMode() == "GSK"){
      // After the waning period there's no effect of the vaccine in the reduction of the relative risk of infection
      // The waning time is approximately tau * 4, being waning = exp(-t/tau), the RR should go up from RR(0) to 1
      double vaxWaning = getPreviousInfections() > 0 ? vaxWaning_pos : vaxWaning_neg;
      double wan_ = exp(-double(currentDay - vday) / (vaxWaning));
      RRInf = 1 - (1 - vaccineProfile.getRRInf(getPreviousInfections() > 0,infectionType - 1)) * wan_ ;
      RRDis = 1 - (1 - vaccineProfile.getRRDis(getPreviousInfections() > 0, infectionType - 1)) * wan_;
      RRHosp = 1 - (1 - vaccineProfile.getRRHosp(getPreviousInfections() > 0, infectionType - 1)) * wan_;
      //printf("SEROTYPE %u RRInf %.4f RRDis %.4f\n", infectionType - 1, RRInf, RRDis);
    }
  }

  double vax_protection = 1.0;
  if(isImmuneVax() == true && vaccineProfile.getMode() == "advance"){
    vax_protection = 1.0 - vaccineProfile.getVaccineProtection();
  }

  exposedCount[infectionType - 1]++;
  // Records the date of each exposure
  string ss = dateOfExposures[infectionType - 1].empty() ? std::to_string(currentDay) : ";" + std::to_string(currentDay);
  dateOfExposures[infectionType - 1] += ss;


  if(rGen->getEventProbability() < RRInf * vax_protection){
    infected = true;
    infLocID = infLocID_;
    recent_inf = infectionType;
    recent_dis = 0;
    recent_hosp = 0;
    last_serotype = infectionType;
    R0[infectionType] = 0;
    if(humInf != nullptr){
      if(humInf->isImported()){
	humInf->increaseVisitorR0(infectionType);
      }else{
	humInf->increaseR0(infectionType);
      }
    }

    if(getPreviousInfections() + vaxAdvancement == 0){
      if(rGen->getEventProbability() < (*disRates)[0] * RRDis){
	recent_dis = infectionType;
	symptomatic = true;
	if(rGen->getEventProbability() < (*hospRates)[0] * RRHosp){
	  recent_hosp = infectionType;
	  hospitalized = true;
	}
      }
    } else if(getPreviousInfections() + vaxAdvancement == 1) {
      if(rGen->getEventProbability() < (*disRates)[1] * RRDis){
	recent_dis = infectionType;
	symptomatic = true;
	if(rGen->getEventProbability() < (*hospRates)[1] * RRHosp){
	  recent_hosp = infectionType;
	  hospitalized = true;
	}
      }
    } else {
      if(rGen->getEventProbability() < (*disRates)[2] * RRDis){
	recent_dis = infectionType;
	symptomatic = true;
	if(rGen->getEventProbability() < (*hospRates)[2] * RRHosp){
	  recent_hosp = infectionType;
	  hospitalized = true;
	}
      }
    }
    infection.reset(new Infection(
				  currentDay + 1, currentDay + 15, 0.0, infectionType, getPreviousInfections() == 0, recent_dis > 0, exp(rGen->getRandomNormal() * 0.2701716 + 1.750673)));

    if(symptomatic == true && enrolledInTrial == true){
      if(hospitalized == true){
	reportSymptoms = true;
      }else if(rGen->getEventProbability() < selfReportProb){
	reportSymptoms = true;
      }
    }
    preExposureAtVaccination[infectionType - 1] = getPreviousInfections();

    for(unsigned i = 0;i < N_SERO; i++){
      if(i != (infectionType - 1) && immunity_perm[i + 1] == false){
	preExposureAtVaccination[i]++;
      }
    }

    updateImmunityPerm(infectionType, true);
    setImmunityTemp(true);
    setImmStartDay(currentDay);
    setImmEndDay(currentDay + 15 + rGen->getHumanImmunity());
  }
}


void Human::initiateBodySize(unsigned currDay, RandomNumGenerator& rGen){
  bodySizeBirth = -1.0;
  bodySizeAdult = -1.0;
  while(bodySizeBirth <= 0.0 || bodySizeAdult <= 0.0){
    if(gender == 'F'){
      bodySizeBirth = 0.3085318 + 0.09302602 * rGen.getRandomNormal();
      bodySizeAdult = 1.505055 + 0.12436 * rGen.getRandomNormal();
      bodySizeSlope = (bodySizeAdult - bodySizeBirth) / 6030.0;
    } else {
      bodySizeBirth = 0.3114736 + 0.1532624 * rGen.getRandomNormal();
      bodySizeAdult = 1.712391 + 0.1523652 * rGen.getRandomNormal();
      bodySizeSlope = (bodySizeAdult - bodySizeBirth) / 6809.0;
    }
  }
  updateBodySize(currDay);
}

bool Human::isImmune(unsigned serotype) const {
  bool immunity = false;

  if(immunity_temp){
    immunity = true;
  } else if(immunity_perm.at(serotype)) {
    immunity = true;
  }
  return immunity;
}

bool Human::isPermImmune(unsigned serotype) const{
  if(immunity_perm.at(serotype)){
    return true;
  }else{
    return false;
  }
}

void Human::resetRecent(){
  recent_inf = 0;
  recent_dis = 0;
  recent_hosp = 0;
  last_serotype = -1;
}

void Human::setImmEndDay(unsigned d) {
  immEndDay = d;
}

void Human::setVaxImmEndDay(unsigned d){
  vaxImmEndDay = d;
}


void Human::setImmStartDay(unsigned d) {
  immStartDay = d;
}

void Human::setVaxImmStartDay(unsigned d){
  vaxImmStartDay = d;
}

void Human::setImmunityTemp(bool status) {
  immunity_temp = status;
}

void Human::setVaxImmunity(bool status) {
  vaccineImmunity = status;
}

void Human::setSeroStatusAtVaccination(){
  if(getPreviousInfections() > 0){
    seroStatusAtVaccination = true;
  }
}

int Human::getPreExposureAtVaccination(unsigned sero){
  return preExposureAtVaccination[sero];
}

void Human::updateAttractiveness(unsigned currDay){
  //printf("Updating attractiveness\n");
  updateBodySize(currDay);
  attractiveness = pow(bodySize, 1.541);
}

void Human::updateBodySize(unsigned currDay){
  if(gender == 'F'){
    if(currDay - bday >= 6030){
      bodySize = bodySizeAdult;
    } else{
      bodySize = bodySizeBirth + bodySizeSlope * (currDay - bday);
    }
  } else {
    if(currDay - bday >= 6809){
      bodySize = bodySizeAdult;
    } else{
      bodySize = bodySizeBirth + bodySizeSlope * (currDay - bday);
    }
  }
}

void Human::updateImmunityPerm(unsigned serotype, bool status) {
  immunity_perm.erase(serotype);
  immunity_perm.insert(std::make_pair(serotype,status));
}

void Human::updateRecent(int infIn, int disIn, int hospIn){
  if(infIn > 0){
    recent_inf = infIn;
  }
  if(disIn > 0){
    recent_dis = disIn;
  }
  if(hospIn > 0){
    recent_hosp = hospIn;
  }
}

void Human::vaccinate(int currDay)
{
  vaccinated = true;
  vday = currDay;
}

void Human::vaccinateAdvanceMode(int currDay, RandomNumGenerator& rGen)
{
  vaccinated = true;
  vday = currDay;
  setVaxImmunity(true);
  setVaxImmStartDay(currDay);
  setVaxImmEndDay(currDay + 365.0 + rGen.getVaxHumanImmunity(rGen.getWaningTime(vaccineProfile.getWaning())));
}

void Human::vaccinateGSKMode(int currDay, RandomNumGenerator& rGen)
{
  vaccinated = true;
  vday = currDay;
  vaxWaning_neg = rGen.getVaxHumanImmunity(vaccineProfile.getWaning(false));
  vaxWaning_pos = rGen.getVaxHumanImmunity(vaccineProfile.getWaning(true));
  //    printf("wanpos, %d, wanneg, %d\n", vaxWaning_pos, vaxWaning_neg);
}

void Human::vaccinateWithProfile(int currDay, RandomNumGenerator * rGen, Vaccine  vax){
  vaccineProfile = vax;
  //    vaccineProfile.printVaccine();
  if(vaccineProfile.getVaccineID() != -1){
    for(int i = 0;i < 4; i++){
      exposedCount[i] = 0;
      dateOfExposures[i].clear();
    }
    vaccineDosesReceived = 1;
    vday = currDay;
    if(vaccineProfile.getMode() == "advance"){
      this->vaccinateAdvanceMode(currDay, (*rGen));
    }else if(vaccineProfile.getMode() == "age"){
      this->vaccinate(currDay);
    }else if(vaccineProfile.getMode() == "GSK"){
      this->vaccinateGSKMode(currDay, (*rGen));
    }
    if(vaccineProfile.getDoses() == vaccineDosesReceived){
      vaccineComplete = true;
    }
  }else{
    printf("VaccineProfile is -1 in Human::vaccinateWithProfile\n");
    exit(1);
  }
}

void Human::boostVaccine(int currDay, RandomNumGenerator * rGen){
  if(vaccineProfile.getVaccineID() != -1){
    vaccineDosesReceived++;
    if(vaccineProfile.getMode() == "advance"){
      this->vaccinateAdvanceMode(currDay, (*rGen));
    }else if(vaccineProfile.getMode() == "age"){
      this->vaccinate(currDay);
    }else if(vaccineProfile.getMode() == "GSK"){
      this->vaccinateGSKMode(currDay, (*rGen));
    }
    if(vaccineProfile.getDoses() == vaccineDosesReceived){
      vaccineComplete = true;
    }
  }else{
    printf("VaccineProfile is -1 in Human::boostVaccine\n");
    exit(1);
  }
}

void Human::enrollInTrial(int currDay, string arm_){
  tAge = getAgeDays(currDay);
  trialDay = currDay;
  enrolledInTrial = true;
  trialArm = arm_;
}

int Human::getNextDoseDay(){
  fflush(stdout);
  if(vaccineComplete == true){
    return -1;
  }else{
    if(vaccineProfile.getVaccineID() != -1){
      return vaccineProfile.getNextDoseTime(vday,vaccineDosesReceived);
    }else{
      printf("VaccineProfile is -1 in Human::getNextDoseDay\n");
      exit(1);
    }
  }
}

void Human::updateVaccineEfficacy(int currDay){
  // Sanofi-like vaccine includes a temporary complete immunity that wanes with time
  if(vaccineProfile.getMode() == "advance" && currDay == (int)this->getVaxImmEndDay()){
    this->setVaxImmunity(false);
  }
}

bool Human::testSeropositivity(double sensitivity_in, double specificity_in, RandomNumGenerator& rGen){
  if(getPreviousInfections() > 0){
    if(rGen.getEventProbability() < sensitivity_in){
      return true;
    }else{
      return false;
    }
  }else{
    if(rGen.getEventProbability() < specificity_in){
      return false;
    }else{
      return true;
    }
  }
}
