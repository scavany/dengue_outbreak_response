#include <fstream>
#include <string>
#include <limits>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>
#include "OutbreakResponse.h"
#include "Report.h"

OutbreakResponse::OutbreakResponse(){
  surveillanceEffort = 0.0;
  responseThreshold = 0.0;
  aggressiveness = 0.0;
  thoroughness = 0.0;
  spatialRadius = 0.0;
  residuality = 0.0;
  compliance = 0.0;
  //Add secondary compliance and training efficacy and training length, remove radius

  symptomaticCases = 0;
  totalHousesSprayed = 0;
  surveillanceDelay = 0;
  maxHouses = 0;
  maxHousesPerDay = 0;
  startDay = 0;
  endDay = 0;
  cycleCounter = 0;

  responseStrategy = "none";
  activeOutbreakResponse = false;
}

//setup reads the outbreak response input file and sets the parameters
void OutbreakResponse::setup(string file){
  if (file.length() == 0)
  {
    exit(1);
    std::cout << "Outbreak response file is not good" << std::endl;
  }
  string line;
  std::ifstream infile(file);
  if(!infile.good())
  {
    exit(1);
    std::cout << "Outbreak response file is not good" << std::endl;
  }

  while(getline(infile,line,'\n')){
    this->addParameter(line);
  }
  infile.close();

  surveillanceEffort = this->readParameter("outbreak_surveillance_effort",0.0);
  responseThreshold = this->readParameter("outbreak_response_threshold", 0.0);
  aggressiveness = this->readParameter("outbreak_aggressiveness", 0.0);
  thoroughness = this->readParameter("outbreak_thoroughness",0.0);
  spatialRadius = this->readParameter("outbreak_spatial_radius",0.0);
  residuality = this->readParameter("outbreak_residuality",0.0);
  compliance = this->readParameter("outbreak_compliance",0.0);

  surveillanceDelay = this->readParameter("outbreak_surveillance_delay",0);
  maxHouses = this->readParameter("outbreak_max_houses",std::numeric_limits<int>::max());
  maxHousesPerDay = this->readParameter("outbreak_max_houses_per_day",0);
  startDay = this->readParameter("outbreak_start_day",0);
  endDay = this->readParameter("outbreak_end_day",  std::numeric_limits<int>::max());
  cycleEndDay = endDay;
  simultaneousZones = this->readParameter("outbreak_simultaneous_zones", 0);
  numberSprayCycles = this->readParameter("outbreak_number_spray_cycles", 0);
  totalZones = this->readParameter("outbreak_total_number_zones", 0);
  trainingLength = this->readParameter("outbreak_training_length", 0);
  sprayCycleDelay = this->readParameter("outbreak_spray_cycle_delay", 0);
  returnVisits = this ->readParameter("outbreak_return_visits", 0);
  weekLength = this->readParameter("outbreak_week_length", 0);
  numberYearsForMean = this->readParameter("outbreak_number_mean_years", 3);
  numberSigmas = this->readParameter("outbreak_number_sigmas", 1);
  spraysPerYear = this->readParameter("outbreak_sprays_per_year", 2);
  repeatSprayGap = this->readParameter("outbreak_repeat_spray_gap", 0);
  seasonStartDay = this->readParameter("outbreak_season_start_day", 1);
  efficacyLength = this->readParameter("outbreak_efficacy_length", std::numeric_limits<int>::max());
  residualityLag = this->readParameter("outbreak_residuality_lag", 0);
  
  this->readParameter("outbreak_response_strategy",&responseStrategy);
  this->readParameter("outbreak_zonal_ordering", &zonalOrderFile);
  this->readParameter("outbreak_starting_method", &startingMethod);
  this->readParameter("outbreak_past_incidence", &incidenceFile);
  
  repeatYearly = (responseStrategy == "zonal" &&
		  startingMethod.find("yearly") != std::string::npos);
  if (repeatYearly) {
    if (startingMethod.find("once") != std::string::npos ||
	startingMethod == "yearly") {
      repeatSprayGap = 365;
    } else if (startingMethod.find("twice") != std::string::npos) {
      if (repeatSprayGap < 30 || repeatSprayGap > 335) {
	throw std::runtime_error("In OutbreakResponse.cpp, period between sprays is less than a month, or not specified");
      }
    } else {
      throw std::runtime_error("In OutbreakResponse.cpp, unknown starting method specified");
    }
  }
  
  if (startingMethod.find("threshold") != std::string::npos && responseStrategy == "zonal") {
    //in the threshold method the surveillanceDelay is simply added to the training length
    //so the start of the response is delayed by the length of the training and the surveillance gap
    trainingLength += surveillanceDelay;
    if (startingMethod.find("week") != std::string::npos) {
      subdivisionsOfYear = 52;
      lengthOfSubdivisions = 7;
    } else if (startingMethod.find("month") != std::string::npos) {
      subdivisionsOfYear = 12;
      lengthOfSubdivisions = 30;
    } else {
      throw std::runtime_error("In OutbreakResponse.cpp, threshold response specified, but no frequency given");
    }
    if (startingMethod.find("fixed") != std::string::npos) {
      startingMethod = "fixed threshold";
    } else {
      startingMethod = "threshold";
      this->readIncidenceFile(incidenceFile);
      for (int i = 0; i<(int)subdivisionsOfYear; ++i){
	double mean = calculateMean(i);
	meanIncidence.push_back(mean);
	sigmaIncidence.push_back(calculateSigma(i, mean)); 
      }
    }
    startDay = std::numeric_limits<int>::max();   
  }
  
  printf("surveillance effort: %.2f | response threshold: %.2f | aggressiveness: %.2f | thoroughness: %.2f | spatial radius: %.2f\n", surveillanceEffort, responseThreshold, aggressiveness, thoroughness, spatialRadius);
  printf("residuality: %.2f | compliance: %.2f | surveillance delay: %u | max houses: %u | max houses per day: %u\n", residuality, compliance, surveillanceDelay, maxHouses,maxHousesPerDay);
  printf("training length: %u | simultaneous zones: %u | spray cycle delay: %u | number of cycles: %u\n", trainingLength, simultaneousZones, sprayCycleDelay, numberSprayCycles);
  printf("start day: %u \n", startDay);
  std::cout << "Reading zonal file" << std::endl;
  this->readZonalFile(zonalOrderFile); 
}

void OutbreakResponse::setLocations(map < string,std::unique_ptr<Location> > * locations_in){
  locations_ptr = locations_in;
}



bool OutbreakResponse::isOutbreakResponseActive(unsigned currDay){
  if(this->isOutbreakResponseEnabled() && currDay >= startDay && currDay <= endDay && totalHousesSprayed < maxHouses){
    return(true);
  }else{
    return(false);
  }
}


bool OutbreakResponse::isOutbreakResponseEnabled(){
  if(responseStrategy == "ring"|| responseStrategy == "zonal"){
    return(true);
  }else{
    return(false);
  }
}

string OutbreakResponse::outbreakResponseStrategy()
{
  if(responseStrategy == "ring"|| responseStrategy == "zonal"){
    return(responseStrategy);
  }else{
    return("none");
  }
}

void OutbreakResponse::addSymptomatic(unsigned currDay, Human * h, RandomNumGenerator * rGen, Report * r){
  //Add a case if surveilled correctly (only called if cases occurs):
  if(rGen->getEventProbability() < surveillanceEffort){
    symptomaticCases++; //surveilled cases
    todaySympCounts++; //surveilled cases today
    //Does symptomaticCases account for the delay??? No, delay is from surveillance to response
    if (this->outbreakResponseStrategy() == "ring" && symptomaticCases >= responseThreshold)
    {
      if(rGen->getEventProbability() < aggressiveness)
      {
        r->addOutbreakSymptomatic(currDay,h,true);
        if(this->isOutbreakResponseActive(currDay))
        {
          todaySymptomatics.push_back(h); //pointer to the symptomatics from today, past the threshold
        }
      }
      else
      {
        r->addOutbreakSymptomatic(currDay,h,false);
      }
    }
    if (this->outbreakResponseStrategy() == "zonal")
    {
      r->addOutbreakSymptomatic(currDay,h,false);
    }
  }
}

//EDIT THIS ONE
void OutbreakResponse::update(unsigned currDay, RandomNumGenerator * rGen, Report * r){
  if (this->isOutbreakResponseEnabled()){
    if (currDay % 365 == seasonStartDay) {
      sprayCounter = 0;
    }
    if (!this->isOutbreakResponseActive(currDay)
	&& startingMethod.find("threshold") != std::string::npos
	&& thresholdPassed) {
      startDay = currDay;
    }
    if (this->isOutbreakResponseActive(currDay)){
      if (this->outbreakResponseStrategy() == "ring"){
        printf("Updating Outbreakresponse is active day %d records: %lu symptomatics today: %u beyond threshold: %lu\n",
	       currDay, symptomaticRecords.size(), todaySympCounts, todaySymptomatics.size());
        // Check symptomatic records and get the delayed records, erasing the ones from surveillanceDelay days ago
        if(symptomaticRecords.size() == surveillanceDelay){
          delayedSymptomatics = symptomaticRecords.front();//delayedSymptomatics is those past the responseThreshold surveillanceDelay ago
          symptomaticRecords.erase(symptomaticRecords.begin()); //symptomaticRecords just contains those symptomaticRecords from the last surveillanceDelay days, past the responseThreshold
        }
        symptomaticRecords.push_back(todaySymptomatics);
        todaySymptomatics.clear();

        //Mark houses for spraying

        if(!delayedSymptomatics.empty()){
          printf("Day: %d There are %lu symptomatics delayed\n",currDay, delayedSymptomatics.size());
          for(auto it_h:delayedSymptomatics){
            string house_id_tmp = (*it_h).getHouseID();//get house of the delayed symptomatic
            Location * house_tmp = (*locations_ptr)[house_id_tmp].get();
            size_t n_nearby_locs = house_tmp->getNumberRadiusLocations();
            vector<string> nearby_locs = house_tmp->getRadiusLocations();
            printf("Day: %d house: %s No Locations: %lu\n", currDay,house_id_tmp.c_str(), n_nearby_locs);
            if(rGen->getEventProbability() < 1-std::pow((1-compliance),(returnVisits+1))){//add the house of the case, if compliant
              if(housesEnlisted.empty()){//if no houses yet enlisted
                housesEnlisted.insert(house_id_tmp);//insert to unordered_set of enlisted houses
                housesToSpray.push_back(house_id_tmp);//add to queue of houses to spray
              }else if(housesEnlisted.find(house_id_tmp) == housesEnlisted.end()){ //if house is not already in the unordered set
                housesEnlisted.insert(house_id_tmp);
                housesToSpray.push_back(house_id_tmp);
              }
            }
            for(auto it_loc: nearby_locs){//add houses within radius, if compliant
              if(rGen->getEventProbability() < 1-std::pow((1-compliance),(returnVisits+1))){
                if(housesEnlisted.find(it_loc) == housesEnlisted.end()){
                  housesEnlisted.insert(it_loc);
                  housesToSpray.push_back(it_loc);
                }
              }
            }
          }
        }
        delayedSymptomatics.clear();
      } else if (this->outbreakResponseStrategy() == "zonal") {
        //make sure end the intervention correctly
        //make sure increase counter and record end of cycle correctly
        if (zonesToSpray.empty() && housesToSpray.empty()){
          if (cycleCounter == 0){
            //The next expression ensures training does not happen on weekends
            if (currDay >= startDay + 7*((int) (trainingLength/weekLength)) + trainingLength%weekLength + ((int) ((currDay+firstDay)%7)/weekLength)) {
              printf("Day %u training is now complete, enlisting the first set of zones and beginning first cycle\n", currDay);
              //enqueue zones to spray
              zonesToSpray = orderedZones;
              cycleStartDay = currDay;
	      sprayCounter++;
            }
          }
          else if (cycleCounter < numberSprayCycles) {
            //account for the gap between cycles (this gap includes weekends, unlike the training delay):
            if (currDay > cycleEndDay + sprayCycleDelay){
              printf("Day %u: beginning cycle number %u\n", currDay, cycleCounter+1);
              //enqueue later cycles to spray
              zonesToSpray = orderedZones;
              cycleStartDay = currDay;
            }
          }
        }
        if (!zonesToSpray.empty() && housesToSpray.empty()){
          zonesToEnlist.clear();
          while (zonesToEnlist.size() < simultaneousZones && !zonesToSpray.empty()) {
            zonesToEnlist.insert(zonesToSpray.front());
            zonesToSpray.pop();
          }
          for (auto it_loc= locations_ptr->begin(); it_loc != locations_ptr->end(); ++it_loc){
            //Add houses if compliant, assuming returnVisits are made and each visit has independent compliance probability
            if(zonesToEnlist.find(it_loc->second->getMoHID()) != zonesToEnlist.end() && rGen->getEventProbability() < 1-std::pow((1-compliance),(returnVisits+1))){
	      housesEnlisted.insert(it_loc->first);
              housesToSpray.push_back(it_loc->first);
            }
          }
          std::random_shuffle(housesToSpray.begin(), housesToSpray.end());
        }
      }
      if (responseStrategy != "zonal" || (currDay >= startDay + 7*((int) (trainingLength/weekLength)) + trainingLength%weekLength + ((int) ((currDay+firstDay)%7)/weekLength))){
        printf("Day %u in total %lu unique houses enlisted, with %lu houses to spray in this batch and %lu zones to spray in this cycle\n", currDay, housesEnlisted.size(), housesToSpray.size(), zonesToSpray.size());
        this->sprayTodayLocations(currDay, r);
      }
      else if ((int) ((currDay+firstDay)%7)/weekLength){
        printf("Day %u no training (weekend)\n", currDay);
      }
      else {
        printf("Day %u training ongoing\n", currDay);
      }
      if (responseStrategy=="zonal" && housesToSpray.empty() && zonesToSpray.empty()){
        if ((cycleCounter == 0 && (currDay >= startDay + 7*((int) (trainingLength/weekLength)) + trainingLength%weekLength + ((int) ((currDay+firstDay)%7)/weekLength)))||(cycleCounter > 0 && currDay > cycleEndDay + sprayCycleDelay)){
          printf("Day %u completed cycle number %u, which took %u days\n", currDay, cycleCounter+1, 1+currDay-cycleStartDay);
          cycleEndDay = currDay;
          cycleCounter++;
          if (cycleCounter == numberSprayCycles){
            if(repeatYearly){
	      startDay += repeatSprayGap;
	      //the next line changes the repeatSprayGap, if spraying is twice yearly:
	      repeatSprayGap = 365 - repeatSprayGap % 365;
	      cycleCounter = 0;
	    }
	    else if (startingMethod.find("threshold") != std::string::npos){
	      startDay = std::numeric_limits<int>::max();
	      thresholdPassed = false;
	      cycleCounter = 0;
	    }else{
	      endDay = currDay;
	    }
          }
        }
      }
    }
    if(totalHousesSprayed >= maxHouses || currDay == endDay){
      if (currDay < endDay){
        printf("Day %u we have sprayed %u houses, which is the maximum. Spraying will now cease.\n", currDay, totalHousesSprayed);
        endDay = currDay;
      }
      else if (currDay==endDay){
        printf("Day %u we reached the end of the intervention. Spraying will now cease. \nIn total we sprayed %u houses and %lu unique houses\nThe intervention lasted %u days\n", currDay, totalHousesSprayed, housesEnlisted.size(), currDay-startDay+1);
      }
      symptomaticRecords.clear();
      housesEnlisted.clear();
      while(!housesToSpray.empty()){
        housesToSpray.pop_front();
      }
    } 
    if (startingMethod.find("threshold") != std::string::npos){
      thresholdPassed = this->updateCheckIncidence(currDay, todaySympCounts);
    }
    todaySympCounts = 0;
  }
}

void OutbreakResponse::sprayTodayLocations(unsigned currDay, Report * r){
  if ((int) ((currDay+firstDay)%7)/weekLength){
    printf("Day %u no houses sprayed (weekend), total sprayed %u houses to spray %lu\n", currDay, totalHousesSprayed,housesToSpray.size());
  }
  else {
    unsigned houses_sprayed_today = 0;
    while(houses_sprayed_today < maxHousesPerDay && totalHousesSprayed < maxHouses && !housesToSpray.empty()){
      string loc_id_tmp = housesToSpray.front();
      housesToSpray.pop_front();
      Location * loc_tmp = (*locations_ptr)[loc_id_tmp].get();
      loc_tmp->sprayAdultInsecticide(currDay, thoroughness, residuality,
				     efficacyLength, residualityLag);
      houses_sprayed_today++;
      totalHousesSprayed++;
      r->addOutbreakResponseLocations(currDay, loc_tmp);
      r->updateSprayReport(currDay, loc_id_tmp);
    }
    printf("Day %u houses sprayed %u total sprayed %u houses to spray %lu\n", currDay, houses_sprayed_today, totalHousesSprayed,housesToSpray.size());
  }
}

// Parameters parsing
void OutbreakResponse::addParameter(string line){
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

int OutbreakResponse::readParameter(string param_name, int vtemp){
  map<string, string>::iterator it;
  int values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    values_ = this->parseInteger(it->second);
  }
  return values_;
}
double OutbreakResponse::readParameter(string param_name, double vtemp){
  map<string, string>::iterator it;
  double values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    values_ = this->parseDouble(it->second);
  }
  return values_;
}

bool OutbreakResponse::readParameter(string param_name, bool vtemp){
  map<string, string>::iterator it;
  bool values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    values_ = this->parseBoolean(it->second);
  }
  return values_;
}

bool OutbreakResponse::parseBoolean(string line){
  bool flag_temp = (std::stoi(line.c_str(), NULL, 10) == 0 ? false : true);
  return flag_temp;
}

void OutbreakResponse::readParameter(string param_name, string * param_var){
  map<string, string>::iterator it;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    *param_var = this->parseString(it->second);
  }
}

string OutbreakResponse::parseString(string line){
  size_t first_ = line.find_first_not_of(' ');
  size_t last_ = line.find_last_not_of(' ');
  return line.substr(first_,(last_ - first_ + 1));
}


int OutbreakResponse::parseInteger(string line){
  return strtol(line.c_str(), NULL, 10);
}

double OutbreakResponse::parseDouble(string line){
  return strtod(line.c_str(), NULL);
}

void OutbreakResponse::readZonalFile(string file){
  if (file.length() == 0 && responseStrategy == "zonal") {
    throw std::runtime_error("In OutbreakResponse.cpp, zonalOrdering file has length 0");
  }

  string line, zone;

  std::ifstream infile(file);
  if (!infile.good() && responseStrategy == "zonal") {
    throw std::runtime_error("In OutbreakResponse.cpp, zonalOrdering file is not good to read");
  }

  getline(infile, line);
  while (getline(infile, line)) {
    zone = line;

    while (infile.peek() == '\n'){
      infile.ignore(1, '\n');
    }
    orderedZones.push(zone);
  }
  infile.close();
}

bool OutbreakResponse::updateCheckIncidence(unsigned currDay, unsigned todaySympCounts){
  unsigned dayOfYear = currDay % 365;
  unsigned dayOfSubdivision = dayOfYear % lengthOfSubdivisions;
  currentIncidence += todaySympCounts;
  bool passed = false;
  unsigned subdivision = dayOfYear/lengthOfSubdivisions;
  if (dayOfSubdivision == (lengthOfSubdivisions - 1) && subdivision != (subdivisionsOfYear - 1)){
    std::cout << "The incidence in subdivision " << subdivision << " was " << currentIncidence;
    if (startingMethod.find("fixed") != std::string::npos) {
      std::cout << ". The threshold is " << responseThreshold;
      std::cout << ". We have sprayed " << sprayCounter << " times this season." << std::endl;
      if (currentIncidence > responseThreshold && sprayCounter < spraysPerYear){
	passed = true;
      }
    } else {
      std::cout << ". The threshold is " << meanIncidence[subdivision] << " + ";
      std::cout << numberSigmas << " * " << sigmaIncidence[subdivision] << " = ";
      std::cout << meanIncidence[subdivision] + numberSigmas*sigmaIncidence[subdivision];
      std::cout << ". We have sprayed " << sprayCounter << " times this season." << std::endl;
      if (currentIncidence > meanIncidence[subdivision] + numberSigmas*sigmaIncidence[subdivision]
	  && sprayCounter < spraysPerYear){
	passed = true;
      }
      subdividedIncidence[subdivision].push_back(currentIncidence);
      meanIncidence[subdivision] = calculateMean(subdivision);
      sigmaIncidence[subdivision] = calculateSigma(subdivision, meanIncidence[subdivision]);
    }
    currentIncidence = 0;
  } else if (dayOfYear == 364){
    std::cout << "The incidence in subdivision " << subdivision - 1  << " was " << currentIncidence;
    if (startingMethod.find("fixed") != std::string::npos) {
      std::cout << ". The threshold is " << responseThreshold;
      std::cout << ". We have sprayed " << sprayCounter << " times this season." << std::endl;
      if (currentIncidence > responseThreshold && sprayCounter < spraysPerYear){
	passed = true;
      }
    } else {
      std::cout << ". The threshold is " << meanIncidence[subdivision - 1] << " + ";
      std::cout << numberSigmas << " * " << sigmaIncidence[subdivision - 1] << " = ";
      std::cout << meanIncidence[subdivision - 1] + numberSigmas*sigmaIncidence[subdivision - 1];
      std::cout << ". We have sprayed " << sprayCounter << " times this season." << std::endl;
      if (currentIncidence > meanIncidence[subdivision - 1] + numberSigmas*sigmaIncidence[subdivision - 1]
	  && sprayCounter < spraysPerYear){
	passed = true;
      }
      subdividedIncidence[subdivision - 1].push_back(currentIncidence);
      meanIncidence[subdivision - 1] = calculateMean(subdivision - 1);
      sigmaIncidence[subdivision - 1] = calculateSigma(subdivision - 1, meanIncidence[subdivision - 1]);
    }
    currentIncidence = 0;
  }
  return passed;
}

double OutbreakResponse::calculateMean(unsigned subdivision){
  int startIndex = std::max(0, (int) (subdividedIncidence[subdivision].size() - numberYearsForMean));
  int sum = std::accumulate(subdividedIncidence[subdivision].begin()+startIndex, subdividedIncidence[subdivision].end(),0);
  return (double) sum/std::min(numberYearsForMean, (unsigned) subdividedIncidence[subdivision].size());
}

double OutbreakResponse::calculateSigma(unsigned subdivision, double mean){
  int startIndex = std::max(0, (int)(subdividedIncidence[subdivision].size() - numberYearsForMean));
  double sum = 0;
  for (std::vector<unsigned>::iterator it = subdividedIncidence[subdivision].begin()+startIndex;
       it != subdividedIncidence[subdivision].end(); ++it){
    sum += std::pow(*it - mean, 2);
  }
  if (sum == 0) return 0;
  else {
    return std::sqrt(sum/(std::min(numberYearsForMean,
				   (unsigned) subdividedIncidence[subdivision].size()) - 1));
  }
}
 
void OutbreakResponse::readIncidenceFile(string file){
  std::ifstream infile(file);
  if (!infile.good()) {
    //exit(1);
    throw std::runtime_error("In OutbreakResponse.cpp, incidence file not good");
  }
  string line;
  while(getline(infile, line)){
    std::vector<unsigned> tempvector;
    tempvector.clear();
    std::stringstream cols(line);
    string value;
    while(getline(cols, value, ',')){
      tempvector.push_back(strtoul(value.c_str(), NULL, 0));
    }
    subdividedIncidence.push_back(tempvector);
  }
  if (subdividedIncidence.size() < subdivisionsOfYear ){
    throw std::runtime_error("In OutbreakResponse.cpp, incidence file too short");
  } else if (subdividedIncidence.size() > subdivisionsOfYear) {
    throw std::runtime_error("In OutbreakResponse.cpp, incidence file too long");
  }
}
