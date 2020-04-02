#include <fstream>
#include <string>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include <list>
#include <algorithm>
#include <utility>
#include <sstream>
#include <cmath>
#include "Simulation.h"

using std::runtime_error;
using std::stringstream;

void Simulation::simulate() {
  simEngine();
  out.close();
  outpop.close();
  outprevac.close();
}


string Simulation::readInputs() {
  readSimControlFile(configLine);
  outbreakIntervention.setup(outbreakFile);
  if (earlyStages){
    std::cout << "Reading early stage parameters" << std::endl;
    readFixedParameters(fixedParameterFile);
    readAegyptiImmatureFile(aegyptiRatesFile);
  }
  else if (delayedEmergence){
    std::cout << "Reading Aegypti file for delayed emergence" << std::endl;
    readFixedParameters(fixedParameterFile);
    readAegyptiDelayedFile(aegyptiRatesFile);
    readDensityFile(densityFile);
    readBurnPopnFile(burnPopnFile);
    //Make the last term correspond to current population, as is the case for location numbers:
    densityDependenceTerms.push_back(std::exp(-alpha*std::pow(mosquitoes.size()/totalInitialPopulation, beta)));
    densityDependenceTerms.pop_front();
  }
  else {
    std::cout << "Reading Aegypti file - not including early stages" << std::endl;
    readAegyptiFile(aegyptiRatesFile);
  }

  outputPath.erase(remove(outputPath.begin(), outputPath.end(), '\"'), outputPath.end());

  RandomNumGenerator rgen(rSeed, huImm, emergeFactor, 1 / mozDailyDeathRate.back(), firstBiteRate.back(), halflife, attractShape);
  rGen = rgen;

  RandomNumGenerator rgen2(rSeedInf, huImm, emergeFactor, 1 / mozDailyDeathRate.back(), firstBiteRate.back(), halflife, attractShape);
  rGenInf = rgen2;

  RandomNumGenerator rgen3(rSeedInf, huImm, emergeFactor, 1 / mozDailyDeathRate.back(), firstBiteRate.back(), halflife, attractShape);
  rGenControl = rgen3;


  readDiseaseRatesFile();

  readVaccineSettingsFile();

  readVaccineProfilesFile();

  if(vaccinationStrategy == "random_trial"){
    recruitmentTrial.setupRecruitment(trialSettingsFile, vaccines, outputPath, simName, &rGen);
    //	printf("Vax Sample: %d, Plac Sample: %d\n", recruitmentTrial.getVaccineSampleSize(), recruitmentTrial.getPlaceboSampleSize());
    //	printf("Recruitment day: %d\n",recruitmentTrial.getRecruitmentStartDay());
  }
  if(vaccinationStrategy == "sanofi_trial"){
    readVaccinationGroupsFile();
  }
  readLocationFile(locationFile);

  if (closeLocDist > 1e-5) {
    setLocNeighborhoodGrid();
  }
  readBirthsFile(birthsFile);
  readTrajectoryFile(trajectoryFile);

  outputReport.setupReport(reportsFile,outputPath,simName);

  if(outbreakIntervention.isOutbreakResponseEnabled()){
    if (outbreakIntervention.outbreakResponseStrategy() == "ring"){
      //This function gets the response radius from the input file, then pairs locations within that radius of each other.
      setRadiusLocations(outbreakIntervention.getResponseRadius());
    }
    //This function passes the locations to outbreakIntervention
    outbreakIntervention.setLocations(&locations);
  }  
  return simName;
}



void Simulation::simEngine() {
  deathMoz = 0;
  lifeMoz = 0;
  visitorsCounter = 0;
  bool finishTrial = false;
  while(currentDay < numDays){
    //rGen.showState(48, "# Daily RNG state excerpt: ");
    humanDeaths = 0;
    if(ceil(double(currentDay + 1) / 365.0) != ceil(double(currentDay) / 365.0)){
      year++; //increase year if in last day of year
      if(attractShape < 0){
        updatePop();//demographics
      }
    }
    printf("day %d year %d Humans %lu Mosquitoes %lu visitors %d\n",currentDay, year,humans.size(), mosquitoes.size(), visitorsCounter);
    for(auto itLoc = locations.begin(); itLoc != locations.end(); itLoc++){
      itLoc->second->updateInfectedVisitor();
      if(currentDay == recruitmentTrial.getRecruitmentStartDay() && vaccinationStrategy == "random_trial"){
        itLoc->second->enableBitesCounter();
      }
      if(finishTrial == true){
        itLoc->second->disableBitesCounter();
      }
    }
    if(vaccinationStrategy == "random_trial"){
      //	    printf("Day: %d\n", currentDay);
      if(currentDay == recruitmentTrial.getRecruitmentStartDay()){
        printf("Current Day : %d is recruitment Start Day\n",currentDay);
        selectEligibleTrialParticipants();
      }
      int trialEndDay = recruitmentTrial.update(currentDay);
      if(trialEndDay == (int)currentDay && vaccinationStrategy == "random_trial"){
        finishTrial = true;
      }
    }
    humanDynamics();//movement routine
    outbreakIntervention.update(currentDay,&rGenControl,&outputReport);
    outputReport.printReport(currentDay);
    mosquitoDynamics();	//bitings and population dynamics
    if (delayedEmergence) {
      //update density dependence terms:
      densityDependenceTerms.push_back(std::exp(-alpha*std::pow(mosquitoes.size()/totalInitialPopulation, beta)));
      densityDependenceTerms.pop_front();
      for (auto& iloc : locations){
        iloc.second->updateRecentAdults(mosquitoes.count(iloc.first));
        //only need the below if doing location density dependence:
        //double newTermTemp = std::exp(-alpha*std::pow(mosquitoes.count(iloc.first)/iloc.second->getInitialAdults(), beta));
        //iloc.second->updateRecentDensityDependenceTerms(newTermTemp);

      }
    }
    if (outputReport.printMozHeterogeneityReport(currentDay)) {
      string mozHetFileStr = outputPath+"/heterogeneity_"+simName
	+"_day_"+std::to_string(currentDay)+".csv";
      std::ofstream mozHetFile;
      mozHetFile.open(mozHetFileStr);
      mozHetFile << "locID,xCor,yCor,mozzes\n";
      if (!mozHetFile.good()) {
	printf("cannot create mozHetFile %s\n", mozHetFileStr.c_str());
	exit(1);
      }
      for (auto& x : locations) { //locations is not in this file!!!! Nor is mosquitoes
	mozHetFile << x.first;
	mozHetFile << ",";
	mozHetFile << x.second->getLocX();
	mozHetFile << ",";
	mozHetFile << x.second->getLocY();
	mozHetFile << ",";
	mozHetFile << mosquitoes.count(x.first);
	mozHetFile << "\n";
      }
      mozHetFile.close();
    }
    currentDay++;
    //       	printf("tempStats, %d, %d, %lu\n", currentDay, humanDeaths, mosquitoes.size());
    
  }
  outputReport.finalizeReport(currentDay);


  if(vaccinationStrategy == "random_trial"){
    recruitmentTrial.finalizeTrial(currentDay);
    if(outputReport.isReportHouseBitesEnabled() == true){
      string outputBitesFile = outputPath + "/" + simName + "_houseBites.csv";
      std::ofstream outBites;
      outBites.open(outputBitesFile);

      if (!outBites.good()) {
        printf("cannot create household output bites %s\n", outputBitesFile.c_str());
        exit(1);
      }
      outBites << "houseID,gini\n";
      for(auto itLoc = locations.begin(); itLoc != locations.end(); itLoc++){
        double giniIndex = itLoc->second->calculateGiniIndex();
        outBites << itLoc->second->getLocID().c_str() << "," << giniIndex <<"\n";
      }
      outBites.close();
    }
  }

}

//This just iterates through the humans and updates their attractiveness using a function in Human.cpp
void Simulation::updatePop(){
  for(auto itHum = humans.begin(); itHum != humans.end(); itHum++){
    itHum->second->updateAttractiveness(currentDay);
  }
}


void Simulation::humanDynamics() {
  int age;
  int cohort = 0;
  if(currentDay >= vaccineDay){
    cohort = floor(double(currentDay - vaccineDay) / 365.0) + 1;
  }
  map<unsigned,double>ForceOfImportation;
  ForceOfImportation.clear();
  ForceOfImportation = currentDay < dailyForceOfImportation.size() ? dailyForceOfImportation[currentDay] : dailyForceOfImportation.back();

  map<unsigned,double>ForceOfLocalImportation;
  ForceOfLocalImportation.clear();
  ForceOfLocalImportation = currentDay < dailyForceOfLocalImportation.size() ? dailyForceOfLocalImportation[currentDay] : dailyForceOfLocalImportation.back();
  //printf("size %lu current day %u LocalFOI(1): %.2f (2):%.2f\n", dailyForceOfLocalImportation.size(),currentDay,ForceOfLocalImportation.at(1), ForceOfLocalImportation.at(2));
  for(auto locIt = zones.begin(); locIt != zones.end();){
    std::string tmpstr = (*locIt);
    zones_counts[tmpstr] = 0;
    ++locIt;
  }

  // Newborns !!!
  auto newbornsRange = future_humans.equal_range(currentDay);
  // walk through today's births
  // add to humans, locations
  for (auto it = newbornsRange.first; it != newbornsRange.second; it++){
    //printf("Human %s is born in day %d\n", it->second->getPersonID().c_str(), currentDay);
    // copy-construct shared pointer to human
    sp_human_t h = it->second;
    if(locations.find(h->getHouseID()) != locations.end()){
      for(const auto & loc_str : h->getLocsVisited()){
        // set.insert already performs check for presence
        auto the_loc = locations.find(loc_str);
        if (the_loc != locations.end() ) {
          the_loc->second->addHuman(h);
        }
      }
      h->initializeHuman(currentDay, InitialConditionsFOI,rGen);
      // deep copy of string ref
      // unnecessary?
      string loc_copy(h->getHouseID());
      humans.insert(make_pair(loc_copy, h));
    }
  }

  // finally, delete range
  future_humans.erase(newbornsRange.first, newbornsRange.second);


  //update alive humans
  for(auto & phum : humans) {
    // check logic!! increment moved to loop
    if(phum.second->isDead()){
      // skip this person
      continue;
    }
    //printf("Human %s is born in day %d\n", phum.second->getPersonID().c_str(), phum.second->getBirthday());
    // daily mortality for humans by age
    if(phum.second->getDeathday() == int(currentDay)){
      if(vaccinationStrategy == "random_trial" && phum.second->isEnrolledInTrial() == true){
        printf("Human %s removed from trial\n", phum.second->getPersonID().c_str());
        recruitmentTrial.removeParticipant(phum.second.get(),currentDay,true);
      }
      // name of each location
      for(const auto & loc_str : phum.second->getLocsVisited() ){
        // please check logic
        // grab iterator to map element
        auto loc_it = locations.find(loc_str);
        // error??
        if (loc_it != locations.end() ) {
          // reemove human from the location pointed to
          loc_it->second->removeHuman(phum.second);
        }
      }
      //	    printf("Human is dead %s\n", phum.second->getPersonID().c_str());
      phum.second->kill();
      humanDeaths++;
      continue;
    }
    // update temporary cross-immunity status if necessary
    //this can be called for humans who immEndDay is not set - is this a problem?
    if(currentDay == phum.second->getImmEndDay() && currentDay != 0)
    phum.second->setImmunityTemp(false);

    // update infection status if necessary
    if(phum.second->infection != nullptr){
      if(outbreakIntervention.isOutbreakResponseEnabled()){
        //add to outbreak symptomatic list, if symptomatic and outbreak response enabled.
        if(phum.second->isSymptomatic() && floor(phum.second->infection->getSymptomOnset()) == currentDay){

          //printf("Symptomatic added to outbreak response Day: %d onset: %.2f\n", currentDay, phum.second->infection->getSymptomOnset());
          outbreakIntervention.addSymptomatic(currentDay, (phum.second).get(), &rGenControl, &outputReport);

        }
      }
      phum.second->checkRecovered(currentDay);
    }

    // select movement trajectory for the day
    phum.second->setTrajDay(rGen.getRandomNum(N_TRAJECTORY));

    // simulate possible imported infection, ignore vaccinated individuals
    for(unsigned serotype = 1; serotype <= (N_SERO); serotype++){
      if(rGen.getEventProbability() < ForceOfImportation.at(serotype)){
        if(phum.second->isImported() == false && phum.second->infection == nullptr){
          phum.second->infectImport(currentDay, serotype,&rGenInf,visitorsCounter);
          visitorsCounter++;
          if(phum.second->infection != nullptr){
            outputReport.addImportation(currentDay, serotype,(phum.second).get());
          }
        }
      }
      if(rGen.getEventProbability() < ForceOfLocalImportation.at(serotype)){
        if(phum.second->infection == nullptr){
          phum.second->infect(currentDay, serotype, &rGenInf, &disRates, &hospRates, nullptr, phum.second->getHouseID());
          if(phum.second->infection != nullptr){
            outputReport.addImportation(currentDay, serotype,(phum.second).get());
          }
        }
      }
    }

    //update vaccine immunity if necessary
    if(phum.second->isVaccinated()){
      phum.second->updateVaccineEfficacy(currentDay);
    }


    // Routine vaccination: first record the cohorts depending on the vaccination strategy, then vaccinate.
    age = phum.second->getAgeDays(currentDay);
    if(vaccinationStrategy == "sanofi_trial"){
      if(currentDay == vaccineDay){
        phum.second->setAgeTrialEnrollment(age);
        phum.second->setSeroStatusAtVaccination();
        phum.second->setCohort(cohort);
        //one-time vaccination by age groups in the trial
        if(checkAgeToVaccinate(age)){
          if(rGenInf.getEventProbability() < vaccineCoverage){
            phum.second->vaccinateWithProfile(currentDay, &rGenInf, vaccines.at(vaccineID));
          }
        }
      }
    }else if(vaccinationStrategy == "routine" && currentDay >= vaccineDay){
      // routine vaccination by age if testBeforeVaccine is enabled, then only vaccinate people that test seropositive
      if(age == (int)vaccineAge * 365){
        phum.second->setAgeTrialEnrollment(age);
        phum.second->setSeroStatusAtVaccination();
        phum.second->setCohort(cohort);
        if(rGenInf.getEventProbability() < vaccineCoverage){
          if(testBeforeVaccine == true){
            bool testPositive = phum.second->testSeropositivity(routineTestSensitivity, routineTestSpecificity, rGenInf);
            if(testPositive == true){
              //printf("Vaccinating due to positive test\n");
              phum.second->vaccinateWithProfile(currentDay, &rGenInf, vaccines.at(vaccineID));
            }
            outputReport.updatePossibleVaccineeStatus(testPositive,phum.second->getPreviousInfections() > 0);
          }else{
            phum.second->vaccinateWithProfile(currentDay, &rGenInf, vaccines.at(vaccineID));
          }
        }
      }
    }
    // catchup vaccination
    if(catchupFlag == true && vaccineDay == currentDay && vaccinationStrategy == "routine"){
      if(age > (int)vaccineAge * 365 && age < 18 * 365){
        if(rGenInf.getEventProbability() < vaccineCoverage){
          if(testBeforeVaccine == true){
            bool testPositive = phum.second->testSeropositivity(routineTestSensitivity, routineTestSpecificity, rGenInf);
            if(testPositive == true){
              phum.second->vaccinateWithProfile(currentDay, &rGenInf, vaccines.at(vaccineID));
            }
          }else{
            phum.second->vaccinateWithProfile(currentDay, &rGenInf, vaccines.at(vaccineID));
          }
        }
      }
    }
    outputReport.updateReport(currentDay,(phum.second).get());
    Location * housetmp = locations[phum.second->getHouseID()].get();
    zones_counts[housetmp->getMoHID()]++;
  }

}

void Simulation::mosquitoDynamics(){
  double biteTime, dieTime;
  //bool biteTaken;

  generateMosquitoes();
  // Read entomological parameters that depend on temperature
  // If there are not enough values, take the last one

  double mozEIP = currentDay < meanDailyEIP.size() ? meanDailyEIP[currentDay] : meanDailyEIP.back();
  //double mozDeathBaseline = currentDay < mozDailyDeathRate.size() ? mozDailyDeathRate[currentDay] : mozDailyDeathRate.back();
  //If moz death is a rate, convert to a risk instead:
  double mozDeathRateBaseline = currentDay < mozDailyDeathRate.size() ? mozDailyDeathRate[currentDay] : mozDailyDeathRate.back();
  double mozDeathBaseline = 1 - std::exp(-mozDeathRateBaseline);
  double mozHalfDeathBaseline = 0;//1 - mozDeathBaseline/mozDeathRateBaseline; //See methods for this formulation - it describes the probability of death for a mosquito on the day it was born, assuming births evenly distributed through the day.
  double mozSBiteRate = currentDay < secondBiteRate.size() ? secondBiteRate[currentDay] : secondBiteRate.back();
  for(auto it = mosquitoes.begin(); it != mosquitoes.end();){
    double mozDeath = locations[it->second->getLocationID()]->getIncreasedMortalityInsecticide(currentDay, mozDeathBaseline);
    double mozHalfDeath = locations[it->second->getLocationID()]->getIncreasedMortalityInsecticide(currentDay, mozHalfDeathBaseline);
    if(it->second->infection != nullptr){
      if(it->second->infection->getStartDay() < 0 && rGenInf.getEventProbability() < rGenInf.getMozLatencyRate(mozEIP)){
        //            if(currentDay == it->second->infection->getStartDay())
        //		printf("Setting infection for mosquito day %d\n",currentDay);
        it->second->infection->setInfectiousnessMosquito(mozInfectiousness, currentDay);
      }
    }
    // determine if the mosquito will bite and/or die today, and if so at what time
    biteTime = double(numDays + 1);
    dieTime = double(numDays + 1);

    if(it->second->getBiteStartDay() <= double(currentDay + 1)){
      biteTime = it->second->getBiteStartDay() - double(currentDay);
      if(biteTime < 0.0){
        biteTime = rGen.getEventProbability();
      }
    }
    // If the mosquito dies today, then set a time to day during the day
    // Right now that time is being set with an uniform distribution -- Double check!!!!
    if (it->second->getBirthDay() == (int)currentDay && currentDay > 0){
      if (rGen.getEventProbability() < mozHalfDeath) {
        dieTime =rGen.getEventProbability();
      }
    }
    else if(rGen.getEventProbability() < mozDeath){
      dieTime =rGen.getEventProbability();
      // it->second->getDDay() - double(currentDay);
    }

    // if the mosquito dies first, then kill it
    if(dieTime <= biteTime && dieTime <= 1.0){
      auto it_temp = it;
      lifeMoz += currentDay - it->second->getBirthDay();
      deathMoz++;
      ++it;
      mosquitoes.erase(it_temp);
      continue;
    }

    // if the mosquito bites first, then let it bite and then see about dying
    while(biteTime < dieTime && biteTime <= 1.0){
      //biteTaken = it->second->takeBite(biteTime,locations[it->second->getLocationID()].get(),&rGen,&rGenInf,&disRates,&hospRates,currentDay,numDays,&out, mozEIP);
      it->second->takeBite(biteTime,locations[it->second->getLocationID()].get(),&rGen,&rGenInf,&disRates,&hospRates,currentDay,numDays,&out, mozEIP);
      //Always assign the secondary biting rate, don't keep trying
      it->second->setBiteStartDay(currentDay + rGen.getMozRestDays(mozSBiteRate));
      biteTime = it->second->getBiteStartDay() - double(currentDay);
    }
    if(dieTime < 1.0){
      auto it_temp = it;
      lifeMoz += currentDay - it->second->getBirthDay();
      deathMoz++;
      ++it;
      mosquitoes.erase(it_temp);
      continue;
    }
    //update report for mosquitoes
    outputReport.updateMosquitoReport(currentDay,(it->second).get());
    // let the mosquito move if that happens today
    double moveProb = rGen.getEventProbability();
    if(moveProb < mozMoveProbability){
      string newLoc = locations.find(it->first)->second->getRandomCloseLoc(rGen);
      if(newLoc != "TOO_FAR_FROM_ANYWHERE"){
	it->second->setLocation(newLoc);
        mosquitoes.insert(make_pair(newLoc, move(it->second)));
        auto it_temp = it;
        ++it;
        mosquitoes.erase(it_temp);
      } else {
        ++it;
      }
    } else {
      ++it;
    }
  }
  if (earlyStages){
    implicitEarlyStages();
  }
}

void Simulation::generateMosquitoes(){
  int mozCount = 0;
  double mozFBiteRate = currentDay < firstBiteRate.size() ? firstBiteRate[currentDay] : firstBiteRate.back();
  //have an extra parameter for early stages
  if (earlyStages) {
    //double eggDeathRate = (currentDay) < eggDailyDeathRate.size() ? eggDailyDeathRate[currentDay+1] : eggDailyDeathRate.back();
    //double larvaeDeathRate = (currentDay) < larvaeDailyDeathRate.size() ? larvaeDailyDeathRate[currentDay+1] : larvaeDailyDeathRate.back();
    //double pupaeDeathRate = (currentDay) < pupaeDailyDeathRate.size() ? pupaeDailyDeathRate[currentDay+1] : pupaeDailyDeathRate.back();
    //double extraDeathRate = correctionTermA*std::pow(correctionTermB, currentDay)*((currentDay) < extraDailyDeathRate.size() ? extraDailyDeathRate[currentDay+1] : extraDailyDeathRate.back());

    //double eggDevelopmentRate = (currentDay) < dailyEggDevelopmentRate.size() ? dailyEggDevelopmentRate[currentDay+1] : dailyEggDevelopmentRate.back();
    //double larvalDevelopmentRate = (currentDay) < dailyLarvalDevelopmentRate.size() ? dailyLarvalDevelopmentRate[currentDay+1] : dailyLarvalDevelopmentRate.back();
    double pupalDevelopmentRate = (currentDay) < dailyPupalDevelopmentRate.size() ? dailyPupalDevelopmentRate[currentDay+1] : dailyPupalDevelopmentRate.back();
    //double gonotrophicCycleRate = (currentDay) < dailyGonotrophicCycleRate.size() ? dailyGonotrophicCycleRate[currentDay+1] : dailyGonotrophicCycleRate.back();
    //double newEggs,newLarvae, newPupae, larvalCapacity;
    for (auto& x : locations){
      //larvalCapacity = x.second->getLarvalCapacity();
      //newEggs = x.second->getEggs()*(1 - eggDeathRate - eggDevelopmentRate) + gonotrophicCycleRate*mosquitoes.count(x.first);
      //newLarvae = x.second->getLarvae()*(1 - larvaeDeathRate - larvalDevelopmentRate - extraDeathRate - x.second->getLarvae()/larvalCapacity) + eggDevelopmentRate*x.second->getEggs();
      //newPupae = x.second->getPupae()*(1 - pupaeDeathRate - pupalDevelopmentRate - extraDeathRate) + larvalDevelopmentRate*x.second->getLarvae();
      
      mozCount = rGen.getMozDevelopment(pupalDevelopmentRate*x.second->getPupae());
      for(int i = 0; i < mozCount; i++){
        unique_ptr<Mosquito> moz(new Mosquito(double(currentDay) + rGen.getMozLifeSpan(), double(currentDay) + rGen.getMozRestDays(mozFBiteRate), x.first));
        moz->setBirthDay(currentDay);
        mosquitoes.insert(make_pair(x.first, move(moz)));
      }
      
      //x.second->updateEggs(newEggs); x.second->updateLarvae(newLarvae); x.second->updatePupae(newPupae);
    }
  } else if (delayedEmergence){
    double meanMoz;
    double emergenceFactor = currentDay < dailyEmergenceFactor.size() ? dailyEmergenceFactor[currentDay] : dailyEmergenceFactor.back();
    vector<double> densityVector = currentDay < densityMatrix.size() ? densityMatrix[currentDay] : densityMatrix.back();
    for(auto& x : locations){
      //This term needs to depend on the integral of past mosquitoes
      //The commented out bits are for including density dependence at the location level
      //Note that densityVector refers to probability density (see equation), and unrelated to mortality
      std::deque<unsigned> recentAdultsTemp = x.second->getRecentAdults();
      //std::deque<double> recentDensDepTermsTemp = x.second->getRecentDensityDependenceTerms();
      auto adults = recentAdultsTemp.begin();
      auto densDepTerm = densityDependenceTerms.begin();
      auto density = densityVector.begin();
      double integral = 0.0;
      while (adults != recentAdultsTemp.end()){
        integral = integral + (*adults) * (*densDepTerm) * (*density);
        ++adults; ++densDepTerm; ++density;
      }
      integral = integral - 0.5*(recentAdultsTemp.front()*densityDependenceTerms.front()*densityVector.front() + recentAdultsTemp.back()*densityDependenceTerms.back()*densityVector.back());
      meanMoz = emergenceFactor*integral;
      mozCount = rGen.getMozEmerge(meanMoz);
      for(int i = 0; i < mozCount; i++){
        unique_ptr<Mosquito> moz(new Mosquito(double(currentDay) + rGen.getMozLifeSpan(), double(currentDay) + rGen.getMozRestDays(mozFBiteRate), x.first));
        moz->setBirthDay(currentDay);
        mosquitoes.insert(make_pair(x.first, move(moz)));

      }
    }


  } else {
    double seasFactor = currentDay < dailyEmergenceFactor.size() ? dailyEmergenceFactor[currentDay] : dailyEmergenceFactor.back();

    for(auto& x : locations){
      mozCount = rGen.getMozEmerge(x.second->getEmergenceRate(), seasFactor * emergeFactor);

      for(int i = 0; i < mozCount; i++){
        unique_ptr<Mosquito> moz(new Mosquito(double(currentDay) + rGen.getMozLifeSpan(), double(currentDay) + rGen.getMozRestDays(mozFBiteRate), x.first));
        moz->setBirthDay(currentDay);
        mosquitoes.insert(make_pair(x.first, move(moz)));
      }
    }
  }
}

void Simulation::implicitEarlyStages(){
  double eggDeathRate = (currentDay +1) < eggDailyDeathRate.size() ? eggDailyDeathRate[currentDay+1] : eggDailyDeathRate.back();
  double larvaeDeathRate = (currentDay +1) < larvaeDailyDeathRate.size() ? larvaeDailyDeathRate[currentDay+1] : larvaeDailyDeathRate.back();
  double pupaeDeathRate = (currentDay +1) < pupaeDailyDeathRate.size() ? pupaeDailyDeathRate[currentDay+1] : pupaeDailyDeathRate.back();
  double extraDeathRate = correctionTermA*std::pow(correctionTermB, currentDay)*((currentDay +1) < extraDailyDeathRate.size() ? extraDailyDeathRate[currentDay+1] : extraDailyDeathRate.back());

  double eggDevelopmentRate = (currentDay +1) < dailyEggDevelopmentRate.size() ? dailyEggDevelopmentRate[currentDay+1] : dailyEggDevelopmentRate.back();
  double larvalDevelopmentRate = (currentDay +1) < dailyLarvalDevelopmentRate.size() ? dailyLarvalDevelopmentRate[currentDay+1] : dailyLarvalDevelopmentRate.back();
  double pupalDevelopmentRate = (currentDay +1) < dailyPupalDevelopmentRate.size() ? dailyPupalDevelopmentRate[currentDay+1] : dailyPupalDevelopmentRate.back();
  double gonotrophicCycleRate = (currentDay +1) < dailyGonotrophicCycleRate.size() ? dailyGonotrophicCycleRate[currentDay+1] : dailyGonotrophicCycleRate.back();
  double emergenceFactor = (currentDay +1) < dailyEmergenceFactor.size() ? dailyEmergenceFactor[currentDay+1] : dailyEmergenceFactor.back();
  double newEggs, newLarvae, newPupae, larvalCapacity, q_param, r_param, s_param, t_param;
  for (auto& x : locations){
    larvalCapacity = emergenceFactor*emergeFactor*x.second->getLarvalCapacity();
    newEggs = (gonotrophicCycleRate*mosquitoes.count(x.first) + x.second->getEggs())/(1 + eggDeathRate + eggDevelopmentRate);
    q_param = larvalCapacity*(1+larvalDevelopmentRate+larvaeDeathRate+extraDeathRate)/3;
    r_param = larvalCapacity*(eggDevelopmentRate*newEggs + x.second->getLarvae())/2;
    s_param = std::cbrt(r_param + std::sqrt(std::pow(q_param, 3.0) + std::pow(r_param, 2.0)));
    t_param = std::cbrt(r_param - std::sqrt(std::pow(q_param, 3.0) + std::pow(r_param, 2.0)));
    newLarvae = s_param + t_param;
    newPupae = (larvalDevelopmentRate*newLarvae + x.second->getPupae())/(1 + pupaeDeathRate + pupalDevelopmentRate + extraDeathRate);
    x.second->updateEggs(newEggs); x.second->updateLarvae(newLarvae); x.second->updatePupae(newPupae);
    outputReport.updateImmatureReport(currentDay, (x.second).get());
    if (mosquitoes.count(x.first) == 0) {
      outputReport.addZeroMozLocation(currentDay);
    }
    // if ((x.first.find("_0001") != std::string::npos) && (currentDay % 10 == 0)){
    //   if (x.first.find("01") != std::string::npos) {
    // 	  std::cout << "Location\t Adults\t Pupae\t Larvae\t Eggs" << std::endl;
    //   }
    //   std::cout << x.first << '\t' << mosquitoes.count(x.first) << '\t' << newPupae << '\t';
    //   std::cout << newLarvae << '\t' << newEggs << std::endl;
    // }
  }
}

void Simulation::setLocNeighborhood(double dist){
  int locCounter = 0;
  std::cout << "Setting up closeLocs" << std::endl;
  for(auto it1 = locations.begin(); it1 != locations.end(); ++it1){
    auto it2 = it1;
    it2++;
    for(; it2 != locations.end(); ++it2){
      if(it1->second->getDistanceFromLoc(*it2->second) <= dist){
        it1->second->addCloseLoc(it2->first);
        it2->second->addCloseLoc(it1->first);
      }
    }
    if (locCounter % 10000 == 0) {
      std::cout << locCounter << " locations setup." <<std::endl;
    }
    locCounter++;
  }
}

void Simulation::setLocNeighborhoodGrid(){
  if(neighborhoodGrid.empty()){
    printf("Grid is empty, please make sure that neighborhoodgrid is properly initialized\n");
    exit(1);
  }    
  for(auto it1 = locations.begin(); it1 != locations.end(); ++it1){
    int x_d = floor(it1->second->getLocX() / closeLocDist);
    int y_d = floor(it1->second->getLocY() / closeLocDist);
    // Find neighbor cells
    unsigned emptyGridCount = 0;
    for(int ix = x_d - 1; ix <= x_d + 1; ix++){
      for(int iy = y_d - 1; iy <= y_d + 1; iy++){
	// Is the cell empty?
	if ((neighborhoodGrid.find(ix) != neighborhoodGrid.end())
	    && (neighborhoodGrid[ix].find(iy) != neighborhoodGrid[ix].end())) {
	  for(unsigned iloc = 0; iloc < neighborhoodGrid[ix][iy].size(); iloc++){
	    if (it1->first != neighborhoodGrid[ix][iy][iloc]) {
	      it1->second->addCloseLoc(neighborhoodGrid[ix][iy][iloc]);
	    }
	  }
	} else if ((ix == x_d) && (iy == y_d)) {
	    throw std::runtime_error("The neighborhood grid is set up wrong");
	} else {
	  emptyGridCount++;
	}
      }
    }
    if ((emptyGridCount == 8) && (neighborhoodGrid[x_d][y_d].size() == 1)) {
      std::cout << "House " << it1->first << " has " << emptyGridCount;
      std::cout << " adjacent cells that are empty, and ";
      std::cout << neighborhoodGrid[x_d][y_d].size() - 1 << " other houses in its grid";
      std::cout << std::endl;
    }
  }
}

void Simulation::setRadiusLocations(double dist){
  printf("setting up locations within %.2f\n", dist);
  for(auto it1 = locations.begin(); it1 != locations.end(); ++it1){
    auto it2 = it1;
    it2++;
    for(; it2 != locations.end(); ++it2){
      if(it1->second->getDistanceFromLoc(*it2->second) <= dist){
        it1->second->addRadiusLoc(it2->first);
        it2->second->addRadiusLoc(it1->first);
      }
    }
  }
}

void Simulation::selectEligibleTrialParticipants(){
  for(auto & hum : humans) {
    recruitmentTrial.addPossibleParticipant(hum.second.get(),currentDay);
  }
  //    printf("In total %lu participants are eligible out of %lu\n",recruitmentTrial.getEligibleParticipantsSize(),humans.size());
  recruitmentTrial.shuffleEligibleParticipants();
  //    recruitmentTrial.printEligibleGroups();
}

void Simulation::readAegyptiFile(string file){
  ifstream infile(file);
  if (!infile.good()) {
    //exit(1);
    throw runtime_error("In Simulation.cpp, Aegypti file not good");
  }
  string line;
  getline(infile,line);
  firstBiteRate.clear();
  secondBiteRate.clear();
  mozDailyDeathRate.clear();
  meanDailyEIP.clear();
  dailyEmergenceFactor.clear();

  while (getline(infile, line, ',')) {
    double eip_temp = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double fb = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double sb = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double dr = strtod(line.c_str(), NULL);
    getline(infile, line, '\n');
    double ef = strtod(line.c_str(), NULL);
    if(eip_temp + fb + sb + dr + ef > 0){
      firstBiteRate.push_back(fb);
      secondBiteRate.push_back(sb);
      mozDailyDeathRate.push_back(dr);
      meanDailyEIP.push_back(eip_temp);
      dailyEmergenceFactor.push_back(ef);
    }
  }
  if(firstBiteRate.empty() || secondBiteRate.empty() || mozDailyDeathRate.empty() || meanDailyEIP.empty() || dailyEmergenceFactor.empty()){
    throw runtime_error("In Simulation.cpp, Aegypti file has empty values");
    //exit(1);
  }
  infile.close();
}

void Simulation::readAegyptiImmatureFile(string file){
  ifstream infile(file);
  if (!infile.good()) {
    //exit(1);
    throw runtime_error("In Simulation.cpp, Aegypti file not good");
  }
  string line;
  getline(infile,line);
  firstBiteRate.clear();
  secondBiteRate.clear();
  meanDailyEIP.clear();
  eggDailyDeathRate.clear();
  larvaeDailyDeathRate.clear();
  pupaeDailyDeathRate.clear();
  mozDailyDeathRate.clear();
  dailyEggDevelopmentRate.clear();
  dailyLarvalDevelopmentRate.clear();
  dailyPupalDevelopmentRate.clear();
  dailyGonotrophicCycleRate.clear();
  dailyEmergenceFactor.clear();

  while (getline(infile, line, ',')) {
    double eip_temp = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double fb = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double sb = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double edr = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double ldr = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double pdr = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double mdr = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double xdr = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double edev = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double ldev = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double pdev = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double gonrate = strtod(line.c_str(), NULL);
    getline(infile, line, '\n');
    double efactor = strtod(line.c_str(), NULL);
    if(eip_temp + fb + sb + edr +ldr + pdr + mdr + xdr + edev + ldev + pdev + gonrate + efactor > 0){
      firstBiteRate.push_back(fb);
      secondBiteRate.push_back(sb);
      meanDailyEIP.push_back(eip_temp);
      eggDailyDeathRate.push_back(edr);
      larvaeDailyDeathRate.push_back(ldr);
      pupaeDailyDeathRate.push_back(pdr);
      mozDailyDeathRate.push_back(mdr);
      extraDailyDeathRate.push_back(xdr);
      dailyEggDevelopmentRate.push_back(edev);
      dailyLarvalDevelopmentRate.push_back(ldev);
      dailyPupalDevelopmentRate.push_back(pdev);
      dailyGonotrophicCycleRate.push_back(gonrate);
      dailyEmergenceFactor.push_back(efactor);
    }
  }
  if(firstBiteRate.empty() || secondBiteRate.empty() || mozDailyDeathRate.empty() || meanDailyEIP.empty() || eggDailyDeathRate.empty()|| larvaeDailyDeathRate.empty()|| pupaeDailyDeathRate.empty()|| extraDailyDeathRate.empty()|| dailyEggDevelopmentRate.empty()|| dailyLarvalDevelopmentRate.empty()|| dailyPupalDevelopmentRate.empty()||dailyGonotrophicCycleRate.empty()||dailyEmergenceFactor.empty()){
    throw runtime_error("In Simulation.cpp, Aegypti file has empty values");
    //exit(1);
  }
  infile.close();
}

void Simulation::readAegyptiDelayedFile(string file){
  ifstream infile(file);
  if (!infile.good()) {
    //exit(1);
    throw runtime_error("In Simulation.cpp, Aegypti file not good");
  }
  string line;
  getline(infile,line);
  meanDailyEIP.clear();
  firstBiteRate.clear();
  secondBiteRate.clear();
  mozDailyDeathRate.clear();
  dailyEmergenceFactor.clear();


  while (getline(infile, line, ',')) {
    double eip_temp = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double fb = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double sb = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double mdr = strtod(line.c_str(), NULL);
    getline(infile, line, '\n');
    double emerge = strtod(line.c_str(), NULL);
    if(eip_temp + fb + sb + mdr + emerge > 0){
      firstBiteRate.push_back(fb);
      secondBiteRate.push_back(sb);
      meanDailyEIP.push_back(eip_temp);
      mozDailyDeathRate.push_back(mdr);
      dailyEmergenceFactor.push_back(emerge);
    }
  }
  if(firstBiteRate.empty() || secondBiteRate.empty() || mozDailyDeathRate.empty() ||dailyEmergenceFactor.empty()){
    throw runtime_error("In Simulation.cpp, Aegypti file has empty values");
    //exit(1);
  }
  infile.close();
}

void Simulation::readDensityFile(string file){
  ifstream infile(file);
  if (!infile.good()) {
    //exit(1);
    throw runtime_error("In Simulation.cpp, Density file not good");
  }
  string line;
  while(getline(infile, line)){
    std::vector<double> tempvector;
    tempvector.clear();
    std::stringstream cols(line);
    string value;
    while(getline(cols, value, ',')){
      tempvector.push_back(strtod(value.c_str(), NULL));
    }
    densityMatrix.push_back(tempvector);
  }
}

void Simulation::readBurnPopnFile(string file){
  ifstream infile(file);
  if (!infile.good()) {
    //exit(1);
    throw runtime_error("In Simulation.cpp, Density file not good");
  }
  string line;
  getline(infile, line);
  normalizedEarlyAdults.clear();
  densityDependenceTerms.clear();

  while (getline(infile, line, ',')) {
    normalizedEarlyAdults.push_back(strtod(line.c_str(), NULL));
    getline(infile, line, '\n');
    densityDependenceTerms.push_back(strtod(line.c_str(), NULL));
  }

}


void Simulation::readFixedParameters(string file) {
  if (file.length() == 0)
  {
    exit(1);
    std::cout << "Fixed parameters file is not good" << std::endl;
  }
  string line;
  std::ifstream infile(file);
  if(!infile.good())
  {
    exit(1);
    std::cout << "Fixed parameters file is not good" << std::endl;
  }

  while(getline(infile,line,'\n')){
    this->addParameter(line);
  }
  infile.close();

  normalizedInitialPupae = this->readParameter("initial_pupae_normalized",0.0);
  normalizedInitialLarvae = this->readParameter("initial_larvae_normalized",0.0);
  normalizedInitialEggs = this->readParameter("initial_eggs_normalized",0.0);
  normalizedLarvalCapacity = this->checkParameterInf("larval_capacity_normalized") ? std::numeric_limits<double>::infinity() : this->readParameter("larval_capacity_normalized",0.0);
  larvalPower = this->readParameter("larval_power",0.0);
  eggsPerGonCycle = this->readParameter("eggs_per_gon_cycle",0.0);
  maxn = this->readParameter("maxn",0);
  this->readParameter("burnin_population", &burnPopnFile);
  this->readParameter("density_values",&densityFile);
  alpha = this->readParameter("alpha",0.0);
  beta = this->readParameter("beta", 0.0);
  totalInitialPopulation = this->readParameter("total_initial_population",0.0);
  extraMortalityRate = this->readParameter("extra_mortality", 0.0);
}

void Simulation::readSimControlFile(string line) {
  std::stringstream infile;
  infile << line;

  getline(infile, line, ',');
  simName = line;

  getline(infile, line, ',');
  rSeed = strtol(line.c_str(), NULL, 10);

  getline(infile, line, ',');
  rSeedInf = strtol(line.c_str(), NULL, 10);

  getline(infile, line, ',');
  numDays = strtol(line.c_str(), NULL, 10);

  getline(infile, line, ',');
  vaccineSettingsFile = line;

  getline(infile, line, ',');
  outputPath = line;

  getline(infile, line, ',');
  reportsFile = line;

  getline(infile, line, ',');
  diseaseRatesFile = line;

  getline(infile, line, ',');
  locationFile = line;

  getline(infile, line, ',');
  trajectoryFile = line;

  getline(infile, line, ',');
  birthsFile = line;

  getline(infile, line, ',');
  string dailyFoiFile = line;

  getline(infile, line, ',');
  string dailyResidentsFoiFile = line;

  getline(infile, line, ',');
  string foiFile = line;

  getline(infile, line, ',');
  huImm = strtol(line.c_str(), NULL, 10);

  getline(infile, line, ',');
  emergeFactor = strtod(line.c_str(), NULL);

  getline(infile, line, ',');
  mozInfectiousness = strtod(line.c_str(), NULL);

  getline(infile, line, ',');
  mozMoveProbability = strtod(line.c_str(), NULL);

  getline(infile, line, ',');
  closeLocDist = strtod(line.c_str(), NULL);
  
  getline(infile, line, ',');
  aegyptiRatesFile = line;

  getline(infile, line, ',');
  fixedParameterFile = line;

  getline(infile, line, ',');
  attractShape = strtod(line.c_str(), NULL);

  getline(infile, line, ',');
  outbreakFile = line;

  getline(infile, line, ',');
  std::istringstream(line) >> earlyStages;

  getline(infile, line, ',');
  std::istringstream(line) >> delayedEmergence;

  getline(infile, line, ',');
  correctionTermA = strtod(line.c_str(), NULL);
 
  getline(infile, line, ',');
  correctionTermB = strtod(line.c_str(), NULL);



  readInitialFOI(foiFile);
  readDailyFOI(dailyFoiFile);
  readDailyLocalFOI(dailyResidentsFoiFile);
}

void Simulation::readInitialFOI(string fileIn){
  if(fileIn.length() == 0){
    throw runtime_error("In Simulation.cpp, InitialFoI file has length 0");
    //exit(1);
  }
  ifstream infile(fileIn);
  string line;
  if(!infile.good()){
    throw runtime_error("In Simulation.cpp, Initial FoI file is not good");
    //exit(1);
  }
  getline(infile, line);
  for(int i = 0; i < N_SERO; i++){
    getline(infile, line, ',');
    double foiTemp = strtod(line.c_str(),NULL);
    InitialConditionsFOI.push_back(foiTemp);
  }
  if(InitialConditionsFOI.size() != N_SERO){
    throw runtime_error("InitialConditionsFOI not set properly\n");
    //exit(1);
  }
}

void Simulation::readDailyFOI(string fileIn){
  if(fileIn.length() == 0){
    printf("In Simulation.cpp, something missing in readDailyFoI %s\n", fileIn.c_str());
    throw runtime_error("In Simulation.cpp, something missing in readDailyFoI");
    //exit(1);
  }
  ifstream infile(fileIn);
  string line;
  if(!infile.good()){
    printf("In Simulation.cpp, something missing, file %s is not good in readDailyFoI\n", fileIn.c_str());
    throw runtime_error("In Simulation.cpp, something missing, file is not good in readDailyFoI");
    //exit(1);
  }
  dailyForceOfImportation.clear();
  getline(infile, line);
  while (getline(infile, line, ',')) {
    double d1 = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double d2 = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double d3 = strtod(line.c_str(), NULL);
    getline(infile, line, '\n');
    double d4 = strtod(line.c_str(), NULL);
    map<unsigned,double> map_temp;
    map_temp.clear();
    map_temp.insert(make_pair(1, d1));
    map_temp.insert(make_pair(2, d2));
    map_temp.insert(make_pair(3, d3));
    map_temp.insert(make_pair(4, d4));
    dailyForceOfImportation.push_back(map_temp);
  }
  if(dailyForceOfImportation.size() == 0){
    throw runtime_error("ForceOfImportation not set\n");
    //exit(1);
  }
}

void Simulation::readDailyLocalFOI(string fileIn){
  if(fileIn.length() == 0){
    printf("In Simulation.cpp, something missing in readDailyLocalFoI %s\n", fileIn.c_str());
    throw runtime_error("In Simulation.cpp, something missing in readDailyLocalFoI");
    //exit(1);
  }
  ifstream infile(fileIn);
  string line;
  if(!infile.good()){
    printf("In Simulation.cpp, something missing, file %s is not good in readDailyLocalFoI\n", fileIn.c_str());
    throw runtime_error("In Simulation.cpp, something missing, file is not good in readDailyLocalFoI");
    //exit(1);
  }
  dailyForceOfLocalImportation.clear();
  getline(infile, line);
  while (getline(infile, line, ',')) {
    double d1 = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double d2 = strtod(line.c_str(), NULL);
    getline(infile, line, ',');
    double d3 = strtod(line.c_str(), NULL);
    getline(infile, line, '\n');
    double d4 = strtod(line.c_str(), NULL);
    map<unsigned,double> map_temp;
    map_temp.clear();
    map_temp.insert(make_pair(1, d1));
    map_temp.insert(make_pair(2, d2));
    map_temp.insert(make_pair(3, d3));
    map_temp.insert(make_pair(4, d4));
    dailyForceOfLocalImportation.push_back(map_temp);
  }
  if(dailyForceOfLocalImportation.size() == 0){
    throw runtime_error("ForceOfLocalImportation not set\n");
    //exit(1);
  }
}


void Simulation::readLocationFile(string locFile) {
  if (locFile.length() == 0) {
    throw runtime_error("In Simulation.cpp, Locations file has length 0");
    //exit(1);
  }
  string line, locID, locType, mohID;
  double x, y, mozzes;

  ifstream infile(locFile);
  if (!infile.good()) {
    throw runtime_error("In Simulation.cpp, locations file is not good to read");
    //exit(1);
  }
  getline(infile, line);
  while (getline(infile, line, ',')) {
    x = strtod(line.c_str(), NULL);

    getline(infile, line, ',');
    y = strtod(line.c_str(), NULL);

    getline(infile, line, ',');
    locType = line;

    getline(infile, line, ',');
    mozzes = strtod(line.c_str(), NULL);

    getline(infile, line, ',');
    mohID = line;

    getline(infile, line);
    locID = line;

    while (infile.peek() == '\n'){
      infile.ignore(1, '\n');
    }
    if (earlyStages){
      double emergenceFactor = (currentDay +1) < dailyEmergenceFactor.size() ? dailyEmergenceFactor[currentDay+1] : dailyEmergenceFactor.back();
      //initializationFactor is to incorporate emergeFactor and dailyEmergenceFactor
      //(which control total average abundance at overal and daily level, respectively,
      //via the larval capacity) into the initial conditions. It is not included in the larval capacity
      //here, as that is done each day
      double initializationFactor = std::pow(emergenceFactor*emergeFactor, 1/larvalPower)*mozzes;
      unique_ptr<Location> location(new Location(locID, locType, mohID, x, y, 
						 normalizedInitialEggs*initializationFactor, 
						 normalizedInitialLarvae*initializationFactor, 
						 normalizedInitialPupae*initializationFactor, 
						 normalizedLarvalCapacity*std::pow(mozzes, 
										   larvalPower)));
      locations.insert(make_pair(locID, move(location)));
      //initialize mosquitoes:
      int initialMozCount = initializationFactor+0.5; //add 0.5 so that rounds to nearest int
      double mozFBiteRate = currentDay < firstBiteRate.size() ? firstBiteRate[currentDay] : firstBiteRate.back();
      for(int i = 0; i < initialMozCount; i++){
        unique_ptr<Mosquito> moz(new Mosquito(double(currentDay) + rGen.getMozLifeSpan(),
					      double(currentDay) + rGen.getMozRestDays(mozFBiteRate), locID));
        moz->setBirthDay(currentDay);
        mosquitoes.insert(make_pair(locID, move(moz)));
      }
    }
    else if (delayedEmergence){
      unsigned initialMozCount =(unsigned)(mozzes+0.5);
      //denormalize vector
      vector<double> earlyAdultsTemp = normalizedEarlyAdults;
      for (auto term = earlyAdultsTemp.begin(); term != earlyAdultsTemp.end(); ++term) {
        *term = mozzes * (*term) + 0.5;
      }
      std::deque<unsigned> earlyAdults(earlyAdultsTemp.begin(), earlyAdultsTemp.end());

      //The density dependent term for each location for each time step is needed if doing location based density dependence rather than population level
      //std::deque<double> earlyDensDepTerms (earlyAdults.begin(), earlyAdults.end());
      //for (auto term = earlyDensDepTerms.begin(); term != earlyDensDepTerms.end(); ++term) {
      //  *term = std::exp(-alpha * std::pow(*term / mozzes, beta));
      //}

      double mozFBiteRate = currentDay < firstBiteRate.size() ? firstBiteRate[currentDay] : firstBiteRate.back();
      for(int i = 0; i < (int)initialMozCount; i++){
        unique_ptr<Mosquito> moz(new Mosquito(double(currentDay) + rGen.getMozLifeSpan(), double(currentDay) + rGen.getMozRestDays(mozFBiteRate), locID));
        moz->setBirthDay(currentDay);
        mosquitoes.insert(make_pair(locID, move(moz)));
      }
      unique_ptr<Location> location(new Location(locID, locType, x, y, earlyAdults, mozzes));
      location->updateRecentAdults(initialMozCount);
      locations.insert(make_pair(locID, move(location)));

    }
    else {
      unique_ptr<Location> location(new Location(locID, locType, mohID, x, y, mozzes));
      locations.insert(make_pair(locID, move(location)));
    }
    zones.insert(mohID);
    if (closeLocDist > 1e-5) {
      // Discretize x and y and add to the grid
      int x_d = floor(x / closeLocDist);
      int y_d = floor(y / closeLocDist);
      neighborhoodGrid[x_d][y_d].push_back(locID);
    }
  }
  infile.close();
  outputReport.setupZones(zones);
  if (closeLocDist > 1e-5) {
    printf("There are %lu columns in the grid\n", neighborhoodGrid.size());
  } else {
    std::cout << "There will be no mosquito movement" << std::endl;
  }
  for(auto locIt = zones.begin(); locIt != zones.end();){
    std::string tmpstr = (*locIt);
    zones_counts.insert(make_pair(tmpstr,0));
    ++locIt;
  }
}



void Simulation::readDiseaseRatesFile(){
  if(diseaseRatesFile.length() == 0){
    throw runtime_error("In Simulation.cpp, disease rates file has length 0");
    //exit(1);
  }
  string line;
  unsigned par = 0;
  double parDis;
  double parHosp;

  ifstream infile(diseaseRatesFile);
  if(!infile.good()){
    throw runtime_error("In Simulation.cpp, disease rates file is not good");
    //exit(1);
  }
  while(getline(infile, line, ',')){
    parDis = strtod(line.c_str(), NULL);
    getline(infile, line, '\n');
    parHosp = strtod(line.c_str(), NULL);
    if(par <= 2){
      disRates.insert(make_pair(par, parDis));
      hospRates.insert(make_pair(par, parHosp));
    }
    par++;
  }
  infile.close();
}

void Simulation::readVaccineSettingsFile(){
  vaccinationStrategy = "none";
  catchupFlag = false;
  testBeforeVaccine = false;
  vaccineCoverage = 0.0;
  routineTestSpecificity = 0.0;
  routineTestSensitivity = 0.0;
  if(vaccineSettingsFile.length() == 0){
    printf("In Simulation.cpp, readvaccinesettings, something missing length is 0 for %s\n", vaccineSettingsFile.c_str());
    throw runtime_error("In Simulation.cpp, readvaccinesettings, something missing length is 0 for vaccine file");
    //exit(1);
  }
  string line;
  ifstream infile(vaccineSettingsFile);
  if(!infile.good()){
    //exit(1);
    printf("In Simulation.cpp (readVaccineSettingsFile) something missing. file %s not good\n", vaccineSettingsFile.c_str());
    throw runtime_error("In Simulation.cpp (readVaccineSettingsFile) something missing. file of vaccine settings is not good\n");
  }
  printf("Reading vaccine settings file %s\n", vaccineSettingsFile.c_str());
  while(getline(infile,line,'\n')){
    string line2,line3;
    vector<string>param_line = getParamsLine(line);
    if(param_line.empty()){
      continue;
    }
    line2 = param_line[0];
    line3 = param_line[1];
    //	printf("%s:%s|\n",line2.c_str(),line3.c_str());
    if(line2 == "vaccination_strategy"){
      vaccinationStrategy = this->parseString(line3);
      printf("VACCINATION STRATEGY:|%s|\n", vaccinationStrategy.c_str());
    }
    if(line2 == "vaccine_day"){
      vaccineDay = this->parseInteger(line3);
    }
    if(line2 == "vaccine_age"){
      vaccineAge = this->parseInteger(line3);
    }
    if(line2 == "vaccine_coverage"){
      vaccineCoverage = this->parseDouble(line3);
    }
    if(line2 == "vaccine_routine_test_specificity"){
      routineTestSpecificity = this->parseDouble(line3);
    }
    if(line2 == "vaccine_routine_test_sensitivity"){
      routineTestSensitivity = this->parseDouble(line3);
    }
    if(line2 == "vaccine_catchup"){
      printf("%s\n", line2.c_str());
      catchupFlag = this->parseInteger(line3) == 1 ? true : false;
    }
    if(line2 == "vaccine_routine_test"){
      printf("%s\n", line2.c_str());
      testBeforeVaccine = this->parseInteger(line3) == 1 ? true : false;
    }
    if(line2 == "vaccine_groups_file"){
      vaccinationGroupsFile = this->parseString(line3);
    }
    if(line2 == "vaccine_profiles_file"){
      vaccineProfilesFile = this->parseString(line3);
    }
    if(line2 == "trial_settings_file"){
      trialSettingsFile = this->parseString(line3);
    }
    if(line2 == "vaccine_ID"){
      vaccineID = this->parseInteger(line3);
    }
  }
  infile.close();
  printf("file %s read\n", vaccineSettingsFile.c_str());
  printf("vStrategy %s day %d age %d coverage: %.2f: groups_file: %s: profilesFile: %s: ID: %d\n",vaccinationStrategy.c_str(), vaccineDay, vaccineAge,vaccineCoverage, vaccinationGroupsFile.c_str(), vaccineProfilesFile.c_str(), vaccineID);
  if(testBeforeVaccine){
    printf("Test before vaccine enabled. Sensitivity %.2f, Specificity %.2f\n", routineTestSensitivity, routineTestSpecificity);
  }else{
    printf("Test before vaccine disabled\n");
  }
}


void Simulation::readVaccineProfilesFile(){
  printf("reading vaccine profiles file\n");
  printf("file = %s\n", vaccineProfilesFile.c_str());

  if(vaccineProfilesFile.length() == 0){
    //exit(1);
    throw runtime_error("In Simulation.cpp, vaccine profile file has length 0");
  }
  vaccines.clear();
  string line;
  ifstream infile(vaccineProfilesFile);
  if(!infile.good()){
    //exit(1);
    throw runtime_error("In Simulation.cpp, vaccine profile file is not good");
  }
  int numVaccines = 0;
  // First find the number of vaccines
  printf("reading vaccine profiles %s\n", vaccineProfilesFile.c_str());
  while(getline(infile,line,'\n')){
    string line2,line3;
    vector<string>param_line = getParamsLine(line);
    if(param_line.empty()){
      continue;
    }
    line2 = param_line[0];
    line3 = param_line[1];
    if(line2 == "vaccine_ids"){
      numVaccines = this->parseInteger(line3);
    }
  }

  if(numVaccines > 0){
    //	printf("There are %d vaccines to read in the file\n",numVaccines);
    // Now we should read the rest of the parameters and store them in the appropiate vaccine structure
    for(int i = 0;i < numVaccines; i++){
      infile.clear();
      infile.seekg(0);
      Vaccine vaxTemp;
      vaxTemp.setID(i);
      while(getline(infile,line,'\n')){
        string line2,line3;
        vector<string>param_line = getParamsLine(line);
        if(param_line.empty()){
          continue;
        }
        line2 = param_line[0];
        line3 = param_line[1];
        string s_id = std::to_string(i);
        if(line2 == "vaccine_mode_" + s_id){
          vaxTemp.setMode(this->parseString(line3));
        }
        if(line2 == "vaccine_name_" + s_id){
          vaxTemp.setName(this->parseString(line3));
        }
        if(line2 == "vaccine_waning_" + s_id){
          vaxTemp.setWaning(this->parseDouble(line3));
        }
        if(line2 == "vaccine_protection_" + s_id){
          vaxTemp.setProtection(this->parseDouble(line3));
        }
        if(line2 == "vaccine_ve_" + s_id){
          vaxTemp.setTotalVE(this->parseDouble(line3));
        }
        if(line2 == "vaccine_vepos_a_" + s_id){
          vaxTemp.addVE_pos(this->parseDouble(line3),0);
        }
        if(line2 == "vaccine_veneg_a_" + s_id){
          vaxTemp.addVE_neg(this->parseDouble(line3),0);
        }
        if(line2 == "vaccine_vepos_b_" + s_id){
          vaxTemp.addVE_pos(this->parseDouble(line3),1);
        }
        if(line2 == "vaccine_veneg_b_" + s_id){
          vaxTemp.addVE_neg(this->parseDouble(line3),1);
        }
        if(line2 == "vaccine_vepos_c_" + s_id){
          vaxTemp.addVE_pos(this->parseDouble(line3),2);
        }
        if(line2 == "vaccine_veneg_c_" + s_id){
          vaxTemp.addVE_neg(this->parseDouble(line3),2);
        }
        if(line2 == "vaccine_veneg_" + s_id){
          vaxTemp.setVaccineEfficacy(false,this->parseDouble(line3));
        }
        if(line2 == "vaccine_vepos_" + s_id){
          vaxTemp.setVaccineEfficacy(true,this->parseDouble(line3));
        }
        if(line2 == "vaccine_RRInfneg_" + s_id){
          vector<double> rrVector;
          this->parseVector(line3, &(rrVector));
          if(rrVector.empty()){
            throw runtime_error("vaccine_RRInfneg has been specified but no value was found\n");
          }else if(rrVector.size() < 4){
            vaxTemp.setRRInf(false,rrVector[0]);
          }else{
            for(unsigned ii = 0; ii < rrVector.size(); ii++){
              vaxTemp.setRRInf(false,rrVector[ii],ii);
            }
          }
        }
        if(line2 == "vaccine_RRInfpos_" + s_id){
          vector<double> rrVector;
          this->parseVector(line3, &(rrVector));
          if(rrVector.empty()){
            throw runtime_error("vaccine_RRInfpos has been specified but no value was found\n");
          }else if(rrVector.size() < 4){
            vaxTemp.setRRInf(true,rrVector[0]);
          }else{
            for(unsigned ii = 0; ii < rrVector.size(); ii++){
              vaxTemp.setRRInf(true,rrVector[ii],ii);
            }
          }
        }
        if(line2 == "vaccine_RRDisneg_" + s_id){
          vector<double> rrVector;
          this->parseVector(line3, &(rrVector));
          if(rrVector.empty()){
            throw runtime_error("vaccine_RRDisneg has been specified but no value was found\n");
          }else if(rrVector.size() < 4){
            vaxTemp.setRRDis(false,rrVector[0]);
          }else{
            for(unsigned ii = 0; ii < rrVector.size(); ii++){
              vaxTemp.setRRDis(false,rrVector[ii],ii);
            }
          }
        }
        if(line2 == "vaccine_RRDispos_" + s_id){
          vector<double> rrVector;
          this->parseVector(line3, &(rrVector));
          if(rrVector.empty()){
            throw runtime_error("vaccine_RRDispos has been specified but no value was found\n");
          }else if(rrVector.size() < 4){
            vaxTemp.setRRDis(true,rrVector[0]);
          }else{
            for(unsigned ii = 0; ii < rrVector.size(); ii++){
              vaxTemp.setRRDis(true,rrVector[ii],ii);
            }
          }
        }
        if(line2 == "vaccine_RRHospneg_" + s_id){
          vector<double> rrVector;
          this->parseVector(line3, &(rrVector));
          if(rrVector.empty()){
            throw runtime_error("vaccine_RRHospneg has been specified but no value was found\n");
          }else if(rrVector.size() < 4){
            vaxTemp.setRRHosp(false,rrVector[0]);
          }else{
            for(unsigned ii = 0; ii < rrVector.size(); ii++){
              vaxTemp.setRRHosp(false,rrVector[ii],ii);
            }
          }
        }
        if(line2 == "vaccine_RRHosppos_" + s_id){
          vector<double> rrVector;
          this->parseVector(line3, &(rrVector));
          if(rrVector.empty()){
            throw runtime_error("vaccine_RRHosppos has been specified but no value was found\n");
          }else if(rrVector.size() < 4){
            vaxTemp.setRRHosp(true,rrVector[0]);
          }else{
            for(unsigned ii = 0; ii < rrVector.size(); ii++){
              vaxTemp.setRRHosp(true,rrVector[ii],ii);
            }
          }
        }
        if(line2 == "vaccine_waning_pos_" + s_id){
          vaxTemp.setWaning(true,this->parseDouble(line3));
        }
        if(line2 == "vaccine_waning_neg_" + s_id){
          vaxTemp.setWaning(false,this->parseDouble(line3));
        }
        if(line2 == "vaccine_prop_inf_" + s_id){
          vaxTemp.setPropInf(this->parseDouble(line3));
        }
        if(line2 == "vaccine_normdev_pos_" + s_id){
          vaxTemp.setNormdevPos(this->parseDouble(line3));
        }
        if(line2 == "vaccine_normdev_neg_" + s_id){
          vaxTemp.setNormdevNeg(this->parseDouble(line3));
        }
        if(line2 == "vaccine_schedule_" + s_id){
          vector<int> rSchedule;
          this->parseVector(line3, &(rSchedule));
          vaxTemp.setRelativeSchedule(rSchedule);
          rSchedule.clear();
        }
      }
      vaccines.insert(make_pair(i,vaxTemp));
      vaxTemp.printVaccine();
    }
    for(unsigned j = 0;j < vaccines.size(); j++){
      vaccines.at(j).printVaccine();
    }
  }else{
    string msg("There are no vaccines to read in ");
    msg += vaccineProfilesFile;
    throw runtime_error(msg);
    //exit(1);
  }
  infile.close();
}
/*
vector<string> Simulation::getParamsLine(string line_){
std::stringstream linetemp;
string line2_,line3_;
linetemp.clear();
linetemp << line_;
getline(linetemp,line2_,'=');
getline(linetemp,line3_,'=');
linetemp.clear();
linetemp << line2_;
getline(linetemp,line2_,' ');
vector<string> params;
params.push_back(line2_);
params.push_back(line3_);
return params;
}*/
vector<string> Simulation::getParamsLine(string line){
  vector<string> params;
  if(line.size() > 0 && line[0] != '#' && line[0] != ' '){
    string param_name, param_value;
    size_t pos_equal = line.find_first_of('=');
    if(pos_equal != string::npos){
      param_name = line.substr(0,pos_equal);
      param_value = line.substr(pos_equal + 1);
      // trim trailing spaces and weird stuff for param_name
      pos_equal = param_name.find_first_of(" \t");
      if(pos_equal != string::npos){
        param_name = param_name.substr(0,pos_equal);
      }
      // trim trailing and leading spaces and weird stuff from param_value
      pos_equal = param_value.find_first_not_of(" \t");
      if(pos_equal != string::npos){
        param_value = param_value.substr(pos_equal);
      }
      pos_equal = param_value.find_first_of("#");
      if(pos_equal != string::npos){
        param_value = param_value.substr(0,pos_equal);
      }
      params.push_back(param_name);
      params.push_back(param_value);
    }
  }
  return params;
}

int Simulation::parseInteger(string line){
  return strtol(line.c_str(), NULL, 10);
}

double Simulation::parseDouble(string line){
  return strtod(line.c_str(), NULL);
}

void Simulation::parseVector(string line, vector<int> * vector_temp){
  std::stringstream linetemp;
  string line2;
  linetemp.clear();
  linetemp << line;
  vector_temp->clear();

  while(getline(linetemp,line2,',')){
    int temp = strtol(line2.c_str(), NULL, 10);
    if(temp >= 0){
      vector_temp->push_back(temp);
    }
  }

  if(vector_temp->empty()){
    throw runtime_error("Parsevector Vector_temp is empty\n");
    //exit(1);
  }
}
void Simulation::parseVector(string line, vector<double> * vector_temp){
  std::stringstream linetemp;
  string line2;
  linetemp.clear();
  linetemp << line;
  vector_temp->clear();

  while(getline(linetemp,line2,',')){
    double temp = strtod(line2.c_str(), NULL);
    if(temp >= 0){
      vector_temp->push_back(temp);
    }
  }

  if(vector_temp->empty()){
    throw runtime_error("Parsevector Vector_temp is empty\n");
    //exit(1);
  }
}
string Simulation::parseString(string line){
  size_t first_ = line.find_first_not_of(" \t#");
  size_t last_ = line.find_last_not_of(" \t#");
  return line.substr(first_,(last_ - first_ + 1));
}

void Simulation::readBirthsFile(string bFile){
  if(bFile.length() == 0){
    throw runtime_error("Incorrect births file\n");
    //exit(1);
  }
  ifstream infile(bFile);
  if(!infile.good()){
    std::cout << "births file is empty: " << bFile.c_str() << "\n";
    throw runtime_error("In Simulation.cpp, something missing");
    //exit(1);
  }
  string line, houseID;
  int bday,dday;
  unsigned hMemID;
  char gen;

  while (getline(infile, line, ',')) {
    houseID = line;
    getline(infile, line, ',');
    hMemID = strtol(line.c_str(), NULL, 10);
    getline(infile, line, ',');
    gen = line[0];
    getline(infile, line, ',');
    bday = strtol(line.c_str(), NULL, 10);
    getline(infile, line, '\n');
    dday = strtol(line.c_str(), NULL, 10);
    // should this be an error condition?
    if(locations.find(houseID) == locations.end()){
      printf("HouseID: %s not found\n", houseID.c_str());
    } else {
      Location * housetmp = locations[houseID].get();
      sp_human_t h(new Human(houseID, hMemID, housetmp->getMoHID(), gen, bday, dday, rGen));
      //if(total_humans_by_id.find(h->getPersonID()) == total_humans_by_id.end()){
      // map.insert checks, only inserts if none present
      total_humans_by_id.insert(make_pair(h->getPersonID(), h));
      //}else{
      //		printf("Human is repeated %s\n", h->getPersonID().c_str());
      //}
    }
  }
  infile.close();
  /*    for(auto it = total_humans_by_id.begin(); it != total_humans_by_id.end(); ++it){
  printf("humans by id %s -> house %s\n", it->first.c_str(),it->second->getHouseID().c_str() );
}*/
}

void Simulation::readTrajectoryFile(string trajFile){
  if(trajFile.length() == 0){
    //exit(1);
    throw runtime_error("In Simulation.cpp, trajectories file has length 0");
  }
  printf("reading %s file with trajectories\n", trajFile.c_str());
  string line, houseID, personID;
  unsigned hMemID;
  //
  ifstream infile(trajFile);
  if(!infile.good()){
    throw runtime_error("In Simulation.cpp, trajectories file is not good");
    //exit(1);
  }
  while(getline(infile, line, ',')){
    unique_ptr<traject_t> ptrajectories(new traject_t);
    for(int itraj = 0; itraj < N_TRAJECTORY; itraj++){
      houseID = line;
      getline(infile, line, ',');
      hMemID = strtol(line.c_str(), NULL, 10);
      personID =houseID + std::to_string(hMemID);
      vpath_t path;
      getline(infile, line);
      std::stringstream ss;
      ss << line;
      while(getline(ss, line, ',')){
        string hID = line;
        getline(ss, line, ',');
        double timeSpent = strtod(line.c_str(), NULL);
        path.push_back(make_pair(hID, timeSpent));
      }
      // directly insert path
      (*ptrajectories)[itraj] = move(path);
      if(itraj < N_TRAJECTORY-1){
        getline(infile, line, ',');
      }
    }

    //error condition? -> just ignore if not valid trajectory
    if(locations.find(houseID) != locations.end()){
      auto ithum = total_humans_by_id.find(personID);
      if( ithum != total_humans_by_id.end()){
        //auto tmpIt = ithum;
        //++ithum;
        sp_human_t h = ithum->second;
        // set
        h->setTrajectories(ptrajectories);
        if(h->getBirthday() < 0){
          for(const auto & loc_str : h->getLocsVisited()) {
            auto itloc = locations.find(loc_str);
            if( itloc != locations.end()){
              itloc->second->addHuman(h);
            }
          }
          h->initializeHuman(currentDay, InitialConditionsFOI,rGen);
          humans.insert(make_pair(houseID, h));
        }else{
          future_humans.insert(make_pair(h->getBirthday(),h));
        }
        // cleanup at end
        // necessary? shouldn't each person be unique? -> I don't understand this question
        total_humans_by_id.erase(ithum);
      }
    }

    while (infile.peek() == '\n'){
      infile.ignore(1, '\n');
    }
  }
  infile.close();
  /*    for(auto itHum = humans.begin(); itHum != humans.end(); itHum++){
  printf("Human %s attractiveness %f\n", itHum->second->getPersonID().c_str(), itHum->second->getAttractiveness());
}*/
total_humans_by_id.clear();
printf("Trajectories have finished successfully\n");
}

void Simulation::readVaccinationGroupsFile(){
  if (vaccinationGroupsFile.length() == 0) {
    //exit(1);
    throw runtime_error("In Simulation.cpp, vaccination groups file has length 0");
  }
  string line;
  int maxAge;
  int minAge;
  ifstream infile(vaccinationGroupsFile);
  if(!infile.good()){
    //exit(1);
    throw runtime_error("In Simulation.cpp,  vaccination groups file is not good");
  }
  while(getline(infile, line, ',')){
    minAge = strtol(line.c_str(), NULL, 10);
    getline(infile, line, '\n');
    maxAge = strtol(line.c_str(), NULL,10);
    ageGroups.insert(make_pair(minAge,maxAge));
  }
  infile.close();
}


bool Simulation::checkAgeToVaccinate(int age_){
  //map<int,int> & the_age
  for(const auto & the_age : ageGroups) {
    //      for(int k = (*itAge).first; k <= (*itAge).second; k++){
    if(age_ >= the_age.first * 365 && age_ <= the_age.second * 365){
      return true;
    }
  }
  return false;
}

int Simulation::readParameter(string param_name, int vtemp){
  map<string, string>::iterator it;
  int values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    values_ = this->parseInteger(it->second);
  }
  return values_;
}

double Simulation::readParameter(string param_name, double vtemp){
  map<string, string>::iterator it;
  double values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    values_ = this->parseDouble(it->second);
  }
  return values_;
}

bool Simulation::checkParameterInf(string param_name){
  map<string, string>::iterator it;
  it = parameters.find(param_name);
  bool infValue = false;
  if(it != parameters.end()){
    if(this->parseString(it->second)=="Inf"){
      infValue = true;
    }
  }
  return infValue;
}

bool Simulation::readParameter(string param_name, bool vtemp){
  map<string, string>::iterator it;
  bool values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    values_ = this->parseBoolean(it->second);
  }
  return values_;
}

string Simulation::readParameter(string param_name, string vtemp){
  map<string, string>::iterator it;
  string values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    values_ = this->parseString(it->second);
  }
  return values_;
}

void Simulation::readParameter(string param_name, string * param_var){
  map<string, string>::iterator it;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    *param_var = this->parseString(it->second);
  }
}


bool Simulation::parseBoolean(string line){
  bool flag_temp = (std::stoi(line.c_str(), NULL, 10) == 0 ? false : true);
  return flag_temp;
}

void Simulation::addParameter(string line){
  if(line.size() > 0 && line[0] != '#' && line[0] != ' '){
    string param_name, param_value;
    size_t pos_equal = line.find_first_of('=');
    if(pos_equal != string::npos){
      param_name = line.substr(0,pos_equal);
      param_value = line.substr(pos_equal + 1);
      // trim trailing spaces and weird stuff for param_name
      pos_equal = param_name.find_first_of(" \t");
      if(pos_equal != string::npos){
        param_name = param_name.substr(0,pos_equal);
      }
      // trim trailing and leading spaces and weird stuff from param_value
      pos_equal = param_value.find_first_not_of(" \t");
      if(pos_equal != string::npos){
        param_value = param_value.substr(pos_equal);
      }
      pos_equal = param_value.find_first_of("#");
      if(pos_equal != string::npos){
        param_value = param_value.substr(0,pos_equal);
      }
      // Add the parameter name and value to the map
      parameters.insert(make_pair(param_name,param_value));
    }
  }
}


Simulation::Simulation(string line) {
  currentDay = 0;
  year = 0;
  configLine = line;
}

//Simulation::Simulation() {
//}



//Simulation::Simulation(const Simulation & orig) {
//}

//Simulation::~Simulation() {
//}
