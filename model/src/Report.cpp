#include <fstream>
#include <string>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <sstream>
#include "Report.h"

using std::stringstream;

Report::Report(){
  reportHouseBites = false;
  reportCohort = false;
  reportAges = false;
  reportGroups = false;
  reportFOI = false;
  reportFOIOutbreakSymptomatics = false;
  reportImmature = false;

  printR0 = false;
  reportSpatial = false;
  printCohortPop = false;
  printAgesPop = false;
  printGroupsPop = false;
  printZonesFOI = false;
  printGroupsAgeFirst = false;
  printGroupsTotalAges = false;
  reportSpray = false;
  reportSprayLocs = false;
  spatialMosquitoes = false;
  spatialOutbreakReport = false;
  spatialSymptomatics = false;
  groupsMaxIndex = 0;
  cohortMaxIndex = 0;
  ageMaxIndex = 0;
  groupsAges.clear();
  groupsStats.clear();
  groupsAvgAgeFirst = 0;
  cohortAges.clear();
  cohortStats.clear();
  secondaryCases.clear();
  dailyFOI.clear();
  dailyImmature.clear();
  dailySprayLocs.clear();
  discreteAges.min = 0;
  discreteAges.max = 0;
  truePositives = 0;
  falsePositives = 0;
  falseNegatives = 0;
  trueNegatives = 0;
  ageStats.clear();
  parameters.clear();
  spatialData.clear();

  for(int i = 0;i < 5;i++){
    cohortEvents[i] = 0;
    groupsEvents[i] = 0;
    ageEvents[i] = 0;
  }
  for(int i = 0; i < 2; i++){
    cohortStatus[i] = 0;
    groupsStatus[i] = 0;
    ageStatus[i] = 0;
  }
  for(int i = 0; i < 3; i++){
    groupsReportPeriod[i] = 0;
    cohortReportPeriod[i] = 0;
    ageReportPeriod[i] = 0;
    foiReportPeriod[i] = 0;
    spatialReportPeriod[i] = 0;
  }
  for(int i = 0; i < 4; i++){
    foiTypes[i] = 0;
    newInfections[i] = 0;
    susceptibles[i] = 0;
    susceptibles_temp[i] = 0;
    mozSusceptibles[i] = 0;
    mozExposed[i] = 0;
    mozInfectious[i] = 0;
    importations[i] = 0;
  }
  for(int i = 0; i < 5; i++){
    // There are four status VacSero+ VacSero- PlacSero+ PlacSero-
    for(int j = 0; j < 4; j++){
      groupsTotalAgeStats.status[j].events[i] = 0;
      groupsTotalAgeStats.status[j].nonevents[i] = 0;
    }
    groupsTotalAgeStats.total.events[i] = 0;
    groupsTotalAgeStats.total.nonevents[i] = 0;
  }

  outCohort.close();
  outAges.close();
  outGroups.close();
  events = {"inf", "dis", "hosp", "serop", "vac"};
  status = {"vac", "plac", "serop", "seron"};
}

void Report::setupZones(set<string> zonesIn){
  zonesToPrint = zonesIn;
  for(auto locIt = zonesToPrint.begin(); locIt != zonesToPrint.end();){
    std::string tmpstr = (*locIt);
    printf("Zone %s in Report\n", tmpstr.c_str());
    ++locIt;
  }
}

void Report::setupReport(string file, string outputPath_, string simName_) {
  //       Read the reporting file and assign variable values
  if (file.length() == 0) {
    exit(1);
  }
  string line;
  std::ifstream infile(file);
  if(!infile.good()){
    exit(1);
  }
  while(getline(infile,line,'\n')){
    this->addParameter(line);
  }
  infile.close();

  this->readParameter("groups_events", "events",groupsEvents);
  this->readParameter("cohort_events","events",cohortEvents);
  this->readParameter("age_events","events", ageEvents);
  this->readParameter("groups_status","status",groupsStatus);
  this->readParameter("cohort_status","status",cohortStatus);
  this->readParameter("age_status","status",ageStatus);
  this->readParameter("groups_report_period","period",groupsReportPeriod);
  this->readParameter("cohort_report_period","period",cohortReportPeriod);
  this->readParameter("age_report_period","period",ageReportPeriod);
  this->readParameter("foi_report_period","period",foiReportPeriod);
  this->readParameter("foi_serotypes","serotypes",foiTypes);
  this->readParameter("spatial_report_period", "period", spatialReportPeriod);

  this->readParameter("groups_ages","ages",&groupsAges);
  this->readParameter("cohort_ages","ages",&cohortAges);

  discreteAges = this->readParameter("age_ages","ages", discreteAges);
  printGroupsPop = this->readParameter("groups_complement", printGroupsPop);
  printCohortPop = this->readParameter("cohort_complement", printCohortPop);
  printAgesPop = this->readParameter("age_complement", printAgesPop);
  printGroupsAgeFirst = this->readParameter("groups_avg_first", printGroupsAgeFirst);
  printGroupsTotalAges = this->readParameter("groups_print_total_ages", printGroupsTotalAges);
  printZonesFOI = this->readParameter("foi_print_zones", printZonesFOI);
  reportGroups = this->readParameter("groups_print", reportGroups);
  reportCohort = this->readParameter("cohort_print", reportCohort);
  reportHouseBites = this->readParameter("housebites_print", reportHouseBites);
  reportAges = this->readParameter("age_print", reportAges);
  reportFOI = this->readParameter("foi_print", reportFOI);
  reportFOIOutbreakSymptomatics = this->readParameter("foi_outbreak_symptomatics", reportFOIOutbreakSymptomatics);
  reportZeroMozLocations = this->readParameter("zero_moz_locations", reportZeroMozLocations);
  printR0 = this->readParameter("foi_print_R0", printR0);
  reportSpatial = this->readParameter("spatial_print", reportSpatial);
  reportImmature = this->readParameter("immature_print", reportImmature);
  reportMozAges = this->readParameter("moz_ages_print", reportMozAges);
  reportSpray = this->readParameter("spray_print", reportSpray);
  reportSprayLocs = this->readParameter("spray_locs_print", reportSprayLocs);
  spatialMosquitoes = this->readParameter("spatial_mosquitoes", spatialMosquitoes);
  spatialSymptomatics = this->readParameter("spatial_symptomatics", spatialSymptomatics);
  spatialOutbreakReport = this->readParameter("spatial_outbreak_report", spatialOutbreakReport);
  reportMozHeterogeneity = this->readParameter("moz_heterogeneity_print", reportMozHeterogeneity);
  mozHeterogeneityPeriod = this->readParameter("moz_heterogeneity_period", 100);

  if(reportFOI == true){
    outputFOIFile = outputPath_ + "/" + simName_ + "_foi.csv";
  }
  if(reportSpatial == true){
    string outputSpatialFile = outputPath_ + "/" + simName_ + "_spatial.csv";
    outSpatial.open(outputSpatialFile);
    if (!outSpatial.good()) {
      exit(1);
    }
  }
  if(reportGroups == true){
    outputGroupsFile = outputPath_ + "/" + simName_ + "_pop.csv";
    outGroups.open(outputGroupsFile);
    if (!outGroups.good()) {
      exit(1);
    }
  }
  if(reportCohort == true){
    outputCohortFile = outputPath_ + "/" + simName_ + "_cohort.csv";
    outCohort.open(outputCohortFile);
    if (!outCohort.good()) {
      exit(1);
    }
  }
  if(reportAges == true){
    outputAgesFile = outputPath_ + "/" + simName_ + "_ages.csv";
    outAges.open(outputAgesFile);
    if (!outAges.good()) {
      exit(1);
    }
  }
  if (reportImmature){
    outputImmatureFile = outputPath_ + "/" + simName_ + "_immature.csv";
  }

  if (reportMozAges){
    outputMozAgesFile = outputPath_ + "/" + simName_ + "_moz_ages.csv";
  }

  if (reportSpray || reportSprayLocs){
    outputSprayFile = outputPath_ + "/" + simName_ + "_spray.csv";
  }
  resetReports();
  printHeaders();
}

void Report::readParameter(string param_name,string param_type, vector<rangeStruct> * values_){
  unordered_map<string, string>::iterator it;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    if(param_type == "ages"){
      if(param_name != "age_ages"){
        parseGroupsAges(it->second,values_);
      }
    }
  }
}

Report::rangeStruct Report::readParameter(string param_name,string param_type, rangeStruct vtemp){
  unordered_map<string, string>::iterator it;
  rangeStruct values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    if(param_type == "ages" && param_name == "age_ages"){
      values_ = parseDiscreteAges(it->second);
    }
  }
  return values_;
}

bool Report::readParameter(string param_name, bool vtemp){
  unordered_map<string, string>::iterator it;
  bool values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    values_ = parseBoolean(it->second);
  }
  return values_;
}

void Report::readParameter(string param_name, string param_type, int * values_){
  unordered_map<string, string>::iterator it;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    if(param_type == "events"){
      parseEvents(it->second,values_,5);
    }
    if(param_type == "status"){
      parseEvents(it->second,values_,2);
    }
    if(param_type == "period"){
      parsePeriod(it->second,values_);
    }
    if(param_type == "serotypes"){
      parseEvents(it->second,values_,4);
    }
  }
}

void Report::readParameter(string param_name, string param_type, unsigned * values_){
  unordered_map<string, string>::iterator it;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    if(param_type == "events"){
      parseEvents(it->second,values_,5);
    }
    if(param_type == "status"){
      parseEvents(it->second,values_,2);
    }
    if(param_type == "period"){
      parsePeriod(it->second,values_);
    }
    if(param_type == "serotypes"){
      parseEvents(it->second,values_,4);
    }
  }
}

int Report::readParameter(string param_name, int vtemp){
  unordered_map<string, string>::iterator it;
  int values_ = vtemp;
  it = parameters.find(param_name);
  if(it != parameters.end()){
    values_ = strtol((it->second).c_str(), NULL, 10);
  }
  return values_;
}

void Report::addParameter(string line){
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

bool Report::parseBoolean(string line){
  bool flag_temp = (std::stoi(line.c_str(), NULL, 10) == 0 ? false : true);
  return flag_temp;
}

void Report::parsePeriod(string line, int * period_temp){
  stringstream linetemp;
  string line2;
  linetemp.clear();
  int count =0;
  linetemp << line;
  while(getline(linetemp,line2,',')){
    if(count > 2){
      break;
    }
    period_temp[count] = strtol(line2.c_str(), NULL, 10);
    if(period_temp[count] < 0){
      period_temp[count] = 0;
    }
    count++;
  }
  // If there are less than 3 values in the period, or the the start:increase:end don't make sense
  if(count < 3 || !(period_temp[0] < period_temp[1] + period_temp[2] && period_temp[1] < period_temp[2])){
    printf("Reporting issue: The period has not been setup properly\n");
    exit(1);
  }
}

void Report::parsePeriod(string line, unsigned * period_temp){
  stringstream linetemp;
  string line2;
  linetemp.clear();
  int count =0;
  linetemp << line;
  while(getline(linetemp,line2,',')){
    if(count > 2){
      break;
    }
    period_temp[count] = strtol(line2.c_str(), NULL, 10);
    if(period_temp[count] < 0){
      period_temp[count] = 0;
    }
    count++;
  }
  // If there are less than 3 values in the period, or the the start:increase:end don't make sense
  if(count < 3 || !(period_temp[0] < period_temp[1] + period_temp[2] && period_temp[1] < period_temp[2])){
    printf("Reporting issue: The period has not been setup properly\n");
    exit(1);
  }
}

void Report::parseGroupsAges(string line, vector<rangeStruct> * ages_temp){
  stringstream linetemp;
  string line2;
  linetemp.clear();
  linetemp << line;
  ages_temp->clear();
  while(getline(linetemp,line2,';')){
    stringstream lTemp; lTemp << line2;
    string line3;
    rangeStruct rangeTemp;
    getline(lTemp,line3,',');
    rangeTemp.min = strtol(line3.c_str(), NULL, 10);
    getline(lTemp,line3,',');
    rangeTemp.max = strtol(line3.c_str(), NULL, 10);
    if(rangeTemp.min + rangeTemp.max > 0){
      ages_temp->push_back(rangeTemp);
    }
  }
  if(ages_temp->empty()){
    printf("Report::ParseGroupsAges ages_temp is empty\n");
    exit(1);
  }
}

Report::rangeStruct Report::parseDiscreteAges(string line){
  stringstream linetemp;
  string line2;
  linetemp.clear();
  linetemp << line;
  rangeStruct rangeTemp;
  getline(linetemp,line2,',');
  rangeTemp.min = strtol(line2.c_str(), NULL, 10);
  getline(linetemp,line2,',');
  rangeTemp.max = strtol(line2.c_str(), NULL, 10);
  if(rangeTemp.min + rangeTemp.max > 0 && rangeTemp.min < rangeTemp.max && rangeTemp.min >= 0){
    return rangeTemp;
  }else{
    printf("Report::parseDiscreteAges minimum and maximum ages do not make sense\n");
    exit(1);
  }
}

void Report::parseEvents(string line, int * Events_, int len){
  stringstream linetemp;
  string line2;
  linetemp.clear();
  int count =0;
  linetemp << line;
  while(getline(linetemp,line2,',')){
    if(count >= len){
      break;
    }
    Events_[count] = strtol(line2.c_str(), NULL, 10);
    if(Events_[count] > 1){
      Events_[count] = 1;
    }
    if(Events_[count] < 0){
      Events_[count] = 0;
    }
    count++;
  }
}

void Report::parseEvents(string line, unsigned * Events_, int len){
  stringstream linetemp;
  string line2;
  linetemp.clear();
  int count =0;
  linetemp << line;
  while(getline(linetemp,line2,',')){
    if(count >= len){
      break;
    }
    Events_[count] = strtol(line2.c_str(), NULL, 10);
    if(Events_[count] > 1){
      Events_[count] = 1;
    }
    if(Events_[count] < 0){
      Events_[count] = 0;
    }
    count++;
  }
}


void Report::updateMosquitoReport(int currDay, Mosquito * m){
  if (reportImmature){
    mosquito_counts[currDay]++;
  }
  if (reportMozAges) {
    int mozAge = currDay - m->getBirthDay();
    mosquito_ages[mozAge][currDay]++;
  }
  if(reportFOI){
    if(m->infection == nullptr){
      for(unsigned i = 1; i < 5; i++){
        dailyFOI[currDay][i]["mozsus"]++;
      }
    }else{
      unsigned sero = m->infection->getInfectionType();
      if(m->infection->getInfectiousness() > 0.0){
        dailyFOI[currDay][sero]["mozinf"]++;
      }else{
        dailyFOI[currDay][sero]["mozexp"]++;
      }
    }
  }
  if(reportSpatial == true && spatialMosquitoes == true){
    if(currDay >= spatialReportPeriod[0] && currDay <= spatialReportPeriod[2] && (currDay - spatialReportPeriod[0]) % spatialReportPeriod[1] == 0){
      if(m->infection != nullptr){
        unsigned sero = m->infection->getInfectionType();
        if(m->infection->getInfectiousness() > 0.0){
          if(currDay == m->infection->getStartDay()){
            string tmp_str = "Mosquitoes," + m->getLocationID()  + "," + std::to_string(currDay) + "," + std::to_string(sero);
            spatialData.push_back(tmp_str);
          }
        }
      }
    }
  }
}

void Report::updateImmatureReport(int currDay, Location * loc){
  if (reportImmature == true){
    dailyImmature[currDay]["eggs"] += loc->getEggs();
    dailyImmature[currDay]["larvae"] += loc->getLarvae();
    dailyImmature[currDay]["pupae"] += loc->getPupae();
  }
}

bool Report::printMozHeterogeneityReport (int currDay) {
  if (reportMozHeterogeneity && (currDay % mozHeterogeneityPeriod == 0)) {
    return true;
  } else {
    return false;
  }
}

void Report::updateSprayReport(int currDay, string loc){
  if (reportSpray){
    spray_counts[currDay]++;
  }
  if (reportSprayLocs){
    dailySprayLocs[currDay].push_back(loc);
  }
}

void Report::updateReport(int currDay, Human * h){
  int reportNum = 0;
  if(reportGroups == true){
    if(currDay >= groupsReportPeriod[0] && currDay <= groupsReportPeriod[2] && (currDay - groupsReportPeriod[0]) % groupsReportPeriod[1] == 0){
      updateGroupsReport(currDay, h);
      reportNum++;
    }
  }
  if(reportCohort == true){
    if(currDay >= cohortReportPeriod[0] && currDay <= cohortReportPeriod[2] && (currDay - cohortReportPeriod[0]) % cohortReportPeriod[1] == 0){
      updateCohortReport(currDay, h);
      reportNum++;
    }
  }
  if(reportAges == true){
    if(currDay >= ageReportPeriod[0] && currDay <= ageReportPeriod[2] && (currDay - ageReportPeriod[0]) % ageReportPeriod[1] == 0){
      updateAgesReport(currDay, h);
      reportNum++;
    }
  }
  if(reportFOI == true){
    if(printR0 == true){
      updateSecondaryCases(currDay,h);
    }
    updateFOI(currDay, h);
  }
  if(reportSpatial == true && spatialMosquitoes == false){
    if(currDay >= spatialReportPeriod[0] && currDay <= spatialReportPeriod[2] && (currDay - spatialReportPeriod[0]) % spatialReportPeriod[1] == 0){
      updateSpatialReport(currDay, h);
    }
  }
  if(reportNum > 0){
    h->resetRecent();
  }
}

void Report::printReport(int currDay){
  if(reportGroups == true){
    if(currDay >= groupsReportPeriod[0] && currDay <= groupsReportPeriod[2] && (currDay - groupsReportPeriod[0]) % groupsReportPeriod[1] == 0){
      printGroupsReport(currDay);
      resetGroupStats();
    }
  }
  if(reportCohort == true){
    if(currDay >= cohortReportPeriod[0] && currDay <= cohortReportPeriod[2] && (currDay - cohortReportPeriod[0]) % cohortReportPeriod[1] == 0){
      printCohortReport(currDay);
      resetCohortStats();
    }
  }
  if(reportAges == true){
    if(currDay >= ageReportPeriod[0] && currDay <= ageReportPeriod[2] && (currDay - ageReportPeriod[0]) % ageReportPeriod[1] == 0){
      printAgesReport(currDay);
      resetAgeStats();
    }
  }
  if(reportSpatial == true){
    if(currDay >= spatialReportPeriod[0] && currDay <= spatialReportPeriod[2] && (currDay - spatialReportPeriod[0]) % spatialReportPeriod[1] == 0){
      printSpatialReport(currDay);
      resetSpatialStats();
    }
  }
}

void Report::updateSpatialReport(int currDay, Human * h){
  if(h->infection != nullptr){
    unsigned sero = h->infection->getInfectionType();
    if(currDay == h->infection->getStartDay()){
      string tmp_str = "InfectiousHum," + h->getHouseID() + "," + std::to_string(currDay) + "," + std::to_string(sero);
      spatialData.push_back(tmp_str);
    }
    if(h->isSymptomatic() && currDay == floor(h->infection->getSymptomOnset())){
      string tmp_str = "SymptomaticHum," + h->getHouseID() + "," + std::to_string(currDay) + "," + std::to_string(sero);
      spatialData.push_back(tmp_str);
    }
  }
}
void Report::addOutbreakResponseLocations(int currDay, Location * locNow){
  dailyOutbreakResponseLocations[currDay]++;
  if(reportSpatial == true && spatialOutbreakReport){
    string tmp_str = "OutbreakResponseLocation," + locNow->getLocID() + "," + std::to_string(currDay) + ",NaN";
    spatialData.push_back(tmp_str);
  }
}
void Report::addOutbreakSymptomatic(int currDay, Human * h, bool respond_){
  dailyOutbreakResponseSymp[currDay]++;
  if(reportSpatial == true && spatialOutbreakReport){
    string tmp_str = "OutbreakReportedHum," + h->getHouseID() + "," + std::to_string(currDay) + ",NaN";
    spatialData.push_back(tmp_str);
    tmp_str = "OutbreakReportedHumInfLoc," + h->getInfLocID() + "," + std::to_string(currDay) + ",NaN";
    spatialData.push_back(tmp_str);
    if(respond_){
      tmp_str = "OutbreakRespondToHum," +  h->getHouseID() + "," + std::to_string(currDay) + ",NaN";
      spatialData.push_back(tmp_str);
      tmp_str = "OutbreakRespondToHumInfLoc," + h->getInfLocID() + "," + std::to_string(currDay) + ",NaN";
      spatialData.push_back(tmp_str);
    }
  }
}

void Report::addZeroMozLocation(int currDay) {
  dailyZeroMozLocations[currDay]++;
}
  
void Report::updateAgesReport(int currDay, Human * h){

  int age_temp = floor((double) h->getAgeDays(currDay) / 365.0);
  if(age_temp > discreteAges.max || age_temp < discreteAges.min){
    return;
  }
  int group = age_temp - discreteAges.min;

  if(h->getRecentInf() > 0){
    ageStats[group].total.events[0]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[0].events[0]++;
      }else{
	ageStats[group].status[1].events[0]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[2].events[0]++;
      }else{
	ageStats[group].status[3].events[0]++;
      }
    }
    if(h->getRecentDis() > 0){
      ageStats[group].total.events[1]++;
      if(h->isVaccinated() == true){
	if(h->getSeroStatusAtVaccination() == true){
	  ageStats[group].status[0].events[1]++;
	}else{
	  ageStats[group].status[1].events[1]++;
	}
      }else{
	if(h->getSeroStatusAtVaccination() == true){
	  ageStats[group].status[2].events[1]++;
	}else{
	  ageStats[group].status[3].events[1]++;
	}
      }
      if(h->getRecentHosp() > 0){
	ageStats[group].total.events[2]++;
	if(h->isVaccinated() == true){
	  if(h->getSeroStatusAtVaccination() == true){
	    ageStats[group].status[0].events[2]++;
	  }else{
	    ageStats[group].status[1].events[2]++;
	  }
	}else{
	  if(h->getSeroStatusAtVaccination() == true){
	    ageStats[group].status[2].events[2]++;
	  }else{
	    ageStats[group].status[3].events[2]++;
	  }
	}
      } else {
	ageStats[group].total.nonevents[2]++;
	if(h->isVaccinated() == true){
	  if(h->getSeroStatusAtVaccination() == true){
	    ageStats[group].status[0].nonevents[2]++;
	  }else{
	    ageStats[group].status[1].nonevents[2]++;
	  }
	}else{
	  if(h->getSeroStatusAtVaccination() == true){
	    ageStats[group].status[2].nonevents[2]++;
	  }else{
	    ageStats[group].status[3].nonevents[2]++;
	  }
	}
      }
    } else {
      ageStats[group].total.nonevents[1]++;
      if(h->isVaccinated() == true){
	if(h->getSeroStatusAtVaccination() == true){
	  ageStats[group].status[0].nonevents[1]++;
	}else{
	  ageStats[group].status[1].nonevents[1]++;
	}
      }else{
	if(h->getSeroStatusAtVaccination() == true){
	  ageStats[group].status[2].nonevents[1]++;
	}else{
	  ageStats[group].status[3].nonevents[1]++;
	}
      }
    }
  } else {
    ageStats[group].total.nonevents[0]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[0].nonevents[0]++;
      }else{
	ageStats[group].status[1].nonevents[0]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[2].nonevents[0]++;
      }else{
	ageStats[group].status[3].nonevents[0]++;
      }
    }
  }

  if(h->getPreviousInfections() > 0){
    ageStats[group].total.events[3]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[0].events[3]++;
      }else{
	ageStats[group].status[1].events[3]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[2].events[3]++;
      }else{
	ageStats[group].status[3].events[3]++;
      }
    }
  } else {
    ageStats[group].total.nonevents[3]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[0].nonevents[3]++;
      }else{
	ageStats[group].status[1].nonevents[3]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[2].nonevents[3]++;
      }else{
	ageStats[group].status[3].nonevents[3]++;
      }
    }
  }
  if(h->isVaccinated() == true){
    ageStats[group].total.events[4]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[0].events[4]++;
      }else{
	ageStats[group].status[1].events[4]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[2].events[4]++;
      }else{
	ageStats[group].status[3].events[4]++;
      }
    }
  } else {
    ageStats[group].total.nonevents[4]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[0].nonevents[4]++;
      }else{
	ageStats[group].status[1].nonevents[4]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	ageStats[group].status[2].nonevents[4]++;
      }else{
	ageStats[group].status[3].nonevents[4]++;
      }
    }
  }
}

void Report::updatePossibleVaccineeStatus(bool test_, bool serostatus_){
  if(test_ == true){
    if(serostatus_ == true){
      truePositives++;
    }else{
      falsePositives++;
    }
  }else{
    if(serostatus_ == false){
      trueNegatives++;
    }else{
      falseNegatives++;
    }
  }
}

void Report::updateGroupsReport(int currDay, Human * h){

  if(h->getRecentInf() && h->getPreviousInfections() == 1){
    groupsAvgAgeFirst += h->getAgeDays(currDay);
    groupsTotalFirstInf++;
  }

  int group = getGroup(h->getAgeDays(currDay), groupsAges);

  if(h->getRecentInf() > 0){
    groupsTotalAgeStats.total.events[0]++;
    if(group >= 0){
      groupsStats[group].total.events[0]++;
    }
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[0].events[0]++;
	if(group >= 0){
	  groupsStats[group].status[0].events[0]++;
	}
      }else{
	groupsTotalAgeStats.status[1].events[0]++;
	if(group >= 0){
	  groupsStats[group].status[1].events[0]++;
	}
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[2].events[0]++;
	if(group >= 0){
	  groupsStats[group].status[2].events[0]++;
	}
      }else{
	groupsTotalAgeStats.status[3].events[0]++;
	if(group >= 0){
	  groupsStats[group].status[3].events[0]++;
	}
      }
    }
    if(h->getRecentDis() > 0){
      groupsTotalAgeStats.total.events[1]++;
      if(group >= 0){
	groupsStats[group].total.events[1]++;
      }
      if(h->isVaccinated() == true){
	if(h->getSeroStatusAtVaccination() == true){
	  groupsTotalAgeStats.status[0].events[1]++;
	  if(group >= 0){
	    groupsStats[group].status[0].events[1]++;
	  }
	}else{
	  groupsTotalAgeStats.status[1].events[1]++;
	  if(group >= 0){
	    groupsStats[group].status[1].events[1]++;
	  }
	}
      }else{
	if(h->getSeroStatusAtVaccination() == true){
	  groupsTotalAgeStats.status[2].events[1]++;
	  if(group >= 0){
	    groupsStats[group].status[2].events[1]++;
	  }
	}else{
	  groupsTotalAgeStats.status[3].events[1]++;
	  if(group >= 0){
	    groupsStats[group].status[3].events[1]++;
	  }
	}
      }
      if(h->getRecentHosp() > 0){
	groupsTotalAgeStats.total.events[2]++;
	if(group >= 0){
	  groupsStats[group].total.events[2]++;
	}
	if(h->isVaccinated() == true){
	  if(h->getSeroStatusAtVaccination() == true){
	    groupsTotalAgeStats.status[0].events[2]++;
	    if(group >= 0){
	      groupsStats[group].status[0].events[2]++;
	    }
	  }else{
	    groupsTotalAgeStats.status[1].events[2]++;
	    if(group >= 0){
	      groupsStats[group].status[1].events[2]++;
	    }
	  }
	}else{
	  if(h->getSeroStatusAtVaccination() == true){
	    groupsTotalAgeStats.status[2].events[2]++;
	    if(group >= 0){
	      groupsStats[group].status[2].events[2]++;
	    }
	  }else{
	    groupsTotalAgeStats.status[3].events[2]++;
	    if(group >= 0){
	      groupsStats[group].status[3].events[2]++;
	    }
	  }
	}
      } else {
	groupsTotalAgeStats.total.nonevents[2]++;
	if(group >= 0){
	  groupsStats[group].total.nonevents[2]++;
	}
	if(h->isVaccinated() == true){
	  if(h->getSeroStatusAtVaccination() == true){
	    groupsTotalAgeStats.status[0].nonevents[2]++;
	    if(group >= 0){
	      groupsStats[group].status[0].nonevents[2]++;
	    }
	  }else{
	    groupsTotalAgeStats.status[1].nonevents[2]++;
	    if(group >= 0){
	      groupsStats[group].status[1].nonevents[2]++;
	    }
	  }
	}else{
	  if(h->getSeroStatusAtVaccination() == true){
	    groupsTotalAgeStats.status[2].nonevents[2]++;
	    if(group >= 0){
	      groupsStats[group].status[2].nonevents[2]++;
	    }
	  }else{
	    groupsTotalAgeStats.status[3].nonevents[2]++;
	    if(group >= 0){
	      groupsStats[group].status[3].nonevents[2]++;
	    }
	  }
	}
      }
    } else {
      groupsTotalAgeStats.total.nonevents[1]++;
      if(group >= 0){
	groupsStats[group].total.nonevents[1]++;
      }
      if(h->isVaccinated() == true){
	if(h->getSeroStatusAtVaccination() == true){
	  groupsTotalAgeStats.status[0].nonevents[1]++;
	  if(group >= 0){
	    groupsStats[group].status[0].nonevents[1]++;
	  }
	}else{
	  groupsTotalAgeStats.status[1].nonevents[1]++;
	  if(group >= 0){
	    groupsStats[group].status[1].nonevents[1]++;
	  }
	}
      }else{
	if(h->getSeroStatusAtVaccination() == true){
	  groupsTotalAgeStats.status[2].nonevents[1]++;
	  if(group >= 0){
	    groupsStats[group].status[2].nonevents[1]++;
	  }
	}else{
	  groupsTotalAgeStats.status[3].nonevents[1]++;
	  if(group >= 0){
	    groupsStats[group].status[3].nonevents[1]++;
	  }
	}
      }
    }
  } else {
    groupsTotalAgeStats.total.nonevents[0]++;
    if(group >= 0){
      groupsStats[group].total.nonevents[0]++;
    }
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[0].nonevents[0]++;
	if(group >= 0){
	  groupsStats[group].status[0].nonevents[0]++;
	}
      }else{
	groupsTotalAgeStats.status[1].nonevents[0]++;
	if(group >= 0){
	  groupsStats[group].status[1].nonevents[0]++;
	}
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[2].nonevents[0]++;
	if(group >= 0){
	  groupsStats[group].status[2].nonevents[0]++;
	}
      }else{
	groupsTotalAgeStats.status[3].nonevents[0]++;
	if(group >= 0){
	  groupsStats[group].status[3].nonevents[0]++;
	}
      }
    }
  }


  if(h->getPreviousInfections() > 0){
    groupsTotalAgeStats.total.events[3]++;
    if(group >= 0){
      groupsStats[group].total.events[3]++;
    }
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[0].events[3]++;
	if(group >= 0){
	  groupsStats[group].status[0].events[3]++;
	}
      }else{
	groupsTotalAgeStats.status[1].events[3]++;
	if(group >= 0){
	  groupsStats[group].status[1].events[3]++;
	}
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[2].events[3]++;
	if(group >= 0){
	  groupsStats[group].status[2].events[3]++;
	}
      }else{
	groupsTotalAgeStats.status[3].events[3]++;
	if(group >= 0){
	  groupsStats[group].status[3].events[3]++;
	}
      }
    }
  } else {
    groupsTotalAgeStats.total.nonevents[3]++;
    if(group >= 0){
      groupsStats[group].total.nonevents[3]++;
    }
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[0].nonevents[3]++;
	if(group >= 0){
	  groupsStats[group].status[0].nonevents[3]++;
	}
      }else{
	groupsTotalAgeStats.status[1].nonevents[3]++;
	if(group >= 0){
	  groupsStats[group].status[1].nonevents[3]++;
	}
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[2].nonevents[3]++;
	if(group >= 0){
	  groupsStats[group].status[2].nonevents[3]++;
	}
      }else{
	groupsTotalAgeStats.status[3].nonevents[3]++;
	if(group >= 0){
	  groupsStats[group].status[3].nonevents[3]++;
	}
      }
    }
  }
  if(h->isVaccinated() == true){
    groupsTotalAgeStats.total.events[4]++;
    if(group >= 0){
      groupsStats[group].total.events[4]++;
    }
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[0].events[4]++;
	if(group >= 0){
	  groupsStats[group].status[0].events[4]++;
	}
      }else{
	groupsTotalAgeStats.status[1].events[4]++;
	if(group >= 0){
	  groupsStats[group].status[1].events[4]++;
	}
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[2].events[4]++;
	if(group >= 0){
	  groupsStats[group].status[2].events[4]++;
	}
      }else{
	groupsTotalAgeStats.status[3].events[4]++;
	if(group >= 0){
	  groupsStats[group].status[3].events[4]++;
	}
      }
    }
  } else {
    groupsTotalAgeStats.total.nonevents[4]++;
    if(group >= 0){
      groupsStats[group].total.nonevents[4]++;
    }
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[0].nonevents[4]++;
	if(group >= 0){
	  groupsStats[group].status[0].nonevents[4]++;
	}
      }else{
	groupsTotalAgeStats.status[1].nonevents[4]++;
	if(group >= 0){
	  groupsStats[group].status[1].nonevents[4]++;
	}
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	groupsTotalAgeStats.status[2].nonevents[4]++;
	if(group >= 0){
	  groupsStats[group].status[2].nonevents[4]++;
	}
      }else{
	groupsTotalAgeStats.status[3].nonevents[4]++;
	if(group >= 0){
	  groupsStats[group].status[3].nonevents[4]++;
	}
      }
    }
  }
}



void Report::updateCohortReport(int currDay, Human * h){
  int cohortNum = h->getCohort();
  if(cohortNum != 1){
    return;
  }
  // get group based on age at trial enrollment
  int cohortAgeGroup = getGroup(h->getAgeTrialEnrollment(),cohortAges);
  if(cohortAgeGroup < 0){
    return;
  }
  if(h->getRecentInf() > 0){
    cohortStats[cohortAgeGroup].total.events[0]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[0].events[0]++;
      }else{
	cohortStats[cohortAgeGroup].status[1].events[0]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[2].events[0]++;
      }else{
	cohortStats[cohortAgeGroup].status[3].events[0]++;
      }
    }
    if(h->getRecentDis() > 0){
      cohortStats[cohortAgeGroup].total.events[1]++;
      if(h->isVaccinated() == true){
	if(h->getSeroStatusAtVaccination() == true){
	  cohortStats[cohortAgeGroup].status[0].events[1]++;
	}else{
	  cohortStats[cohortAgeGroup].status[1].events[1]++;
	}
      }else{
	if(h->getSeroStatusAtVaccination() == true){
	  cohortStats[cohortAgeGroup].status[2].events[1]++;
	}else{
	  cohortStats[cohortAgeGroup].status[3].events[1]++;
	}
      }
      if(h->getRecentHosp() > 0){
	cohortStats[cohortAgeGroup].total.events[2]++;
	if(h->isVaccinated() == true){
	  if(h->getSeroStatusAtVaccination() == true){
	    cohortStats[cohortAgeGroup].status[0].events[2]++;
	  }else{
	    cohortStats[cohortAgeGroup].status[1].events[2]++;
	  }
	}else{
	  if(h->getSeroStatusAtVaccination() == true){
	    cohortStats[cohortAgeGroup].status[2].events[2]++;
	  }else{
	    cohortStats[cohortAgeGroup].status[3].events[2]++;
	  }
	}
      } else {
	cohortStats[cohortAgeGroup].total.nonevents[2]++;
	if(h->isVaccinated() == true){
	  if(h->getSeroStatusAtVaccination() == true){
	    cohortStats[cohortAgeGroup].status[0].nonevents[2]++;
	  }else{
	    cohortStats[cohortAgeGroup].status[1].nonevents[2]++;
	  }
	}else{
	  if(h->getSeroStatusAtVaccination() == true){
	    cohortStats[cohortAgeGroup].status[2].nonevents[2]++;
	  }else{
	    cohortStats[cohortAgeGroup].status[3].nonevents[2]++;
	  }
	}
      }
    } else {
      cohortStats[cohortAgeGroup].total.nonevents[1]++;
      if(h->isVaccinated() == true){
	if(h->getSeroStatusAtVaccination() == true){
	  cohortStats[cohortAgeGroup].status[0].nonevents[1]++;
	}else{
	  cohortStats[cohortAgeGroup].status[1].nonevents[1]++;
	}
      }else{
	if(h->getSeroStatusAtVaccination() == true){
	  cohortStats[cohortAgeGroup].status[2].nonevents[1]++;
	}else{
	  cohortStats[cohortAgeGroup].status[3].nonevents[1]++;
	}
      }
    }
  } else {
    cohortStats[cohortAgeGroup].total.nonevents[0]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[0].nonevents[0]++;
      }else{
	cohortStats[cohortAgeGroup].status[1].nonevents[0]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[2].nonevents[0]++;
      }else{
	cohortStats[cohortAgeGroup].status[3].nonevents[0]++;
      }
    }
  }


  if(h->getPreviousInfections() > 0){
    cohortStats[cohortAgeGroup].total.events[3]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[0].events[3]++;
      }else{
	cohortStats[cohortAgeGroup].status[1].events[3]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[2].events[3]++;
      }else{
	cohortStats[cohortAgeGroup].status[3].events[3]++;
      }
    }
  } else {
    cohortStats[cohortAgeGroup].total.nonevents[3]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[0].nonevents[3]++;
      }else{
	cohortStats[cohortAgeGroup].status[1].nonevents[3]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[2].nonevents[3]++;
      }else{
	cohortStats[cohortAgeGroup].status[3].nonevents[3]++;
      }
    }
  }
  if(h->isVaccinated() == true){
    cohortStats[cohortAgeGroup].total.events[4]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[0].events[4]++;
      }else{
	cohortStats[cohortAgeGroup].status[1].events[4]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[2].events[4]++;
      }else{
	cohortStats[cohortAgeGroup].status[3].events[4]++;
      }
    }
  } else {
    cohortStats[cohortAgeGroup].total.nonevents[4]++;
    if(h->isVaccinated() == true){
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[0].nonevents[4]++;
      }else{
	cohortStats[cohortAgeGroup].status[1].nonevents[4]++;
      }
    }else{
      if(h->getSeroStatusAtVaccination() == true){
	cohortStats[cohortAgeGroup].status[2].nonevents[4]++;
      }else{
	cohortStats[cohortAgeGroup].status[3].nonevents[4]++;
      }
    }
  }
}



int Report::getGroup(int age_, vector<rangeStruct> groups_temp){
  vector<rangeStruct>::iterator itAge = groups_temp.begin();
  int count = 0;
  for(; itAge != groups_temp.end(); itAge++){
    if((double )age_ / 365.0 >= (*itAge).min && (double) age_ / 365.0 < (*itAge).max){
      return count;
    }
    count++;
  }
  return -1;
}

void Report::printSpatialReport(int currDay){
  if(!spatialData.empty()){
    for(auto it = spatialData.begin(); it != spatialData.end(); ++it){
      outSpatial << (*it) << "\n";
    }
  }
}

void Report::printAgesReport(int currDay){
  outAges << currDay << ",";
  for(int i = 0; i < 5 ; i++){
    if(ageEvents[i] == 1){
      if(ageStatus[0] > 0 || ageStatus[1] > 0){
	if(ageStatus[0] > 0 && ageStatus[1] > 0){
	  for(int j = 0; j < 4; j++){
	    for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
	      if(printAgesPop == true){
		outAges << ageStats[k].status[j].nonevents[i]<<",";
	      }
	      outAges << ageStats[k].status[j].events[i];
	      if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min && j == 3){
		outAges << "\n";
	      }else{
		outAges << ",";
	      }
	    }
	  }
	}else{
	  int inc_ = 0;
	  int sum_ = 0;
	  if(ageStatus[0] > 0){
	    inc_ = 2;
	    sum_ = 1;
	  }else if(ageStatus[1] > 0){
	    inc_ = 1;
	    sum_ = 2;
	  }
	  for(int j = 0;j < 1 +  inc_;j += inc_){
	    for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
	      if(printAgesPop == true){
		outAges << ageStats[k].status[j].nonevents[i] + ageStats[k].status[j + sum_].nonevents[i]<<",";
	      }
	      outAges << ageStats[k].status[j].events[i] + ageStats[k].status[j + sum_].events[i];
	      if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min && j == inc_){
		outAges << "\n";
	      }else{
		outAges << ",";
	      }
	    }
	  }
	}

      }else{
	for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
	  if(printAgesPop == true){
	    outAges << ageStats[k].total.nonevents[i]<<",";
	  }
	  outAges << ageStats[k].total.events[i];
	  if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min){
	    outAges << "\n";
	  }else{
	    outAges << ",";
	  }
	}
      }
    }
  }
}



void Report::printGroupsReport(int currDay){
  outGroups << currDay << ",";
  if(printGroupsAgeFirst == true){
    outGroups << (double) groupsAvgAgeFirst / (double) groupsTotalFirstInf / 365.0 << ",";
  }
  vector<string> outdata; string outstring;
  outdata.clear();
  if(groupsStatus[0] > 0 || groupsStatus[1] > 0){
    if(groupsStatus[0] > 0 && groupsStatus[1] > 0){
      for(int i = 0; i < 5 ; i++){
	if(groupsEvents[i] == 1){
	  for(int j = 0;j < 4 ;j++){
	    for(unsigned k = 0; k < groupsAges.size(); k++){
	      if(printGroupsPop == true){
		outdata.push_back(std::to_string(groupsStats[k].status[j].nonevents[i]));
	      }
	      outdata.push_back(std::to_string(groupsStats[k].status[j].events[i]));
	      if(printGroupsTotalAges == true && k == groupsAges.size() - 1){
		if(printGroupsPop == true){
		  outdata.push_back(std::to_string(groupsTotalAgeStats.status[j].nonevents[i]));
		}
		outdata.push_back(std::to_string(groupsTotalAgeStats.status[j].events[i]));
	      }
	    }
	  }
	}
      }
    }else{
      int inc_ = 0;
      int sum_ = 0;
      if(groupsStatus[0] > 0){
	inc_ = 2;
	sum_ = 1;
      }else if(groupsStatus[1] > 0){
	inc_ = 1;
	sum_ = 2;
      }
      for(int i = 0; i < 5 ; i++){
	if(groupsEvents[i] == 1){
	  for(int j = 0;j < 1 + inc_;j+=inc_){
	    for(unsigned k = 0; k < groupsAges.size(); k++){
	      if(printGroupsPop == true){
		outdata.push_back(std::to_string(groupsStats[k].status[j].nonevents[i] + groupsStats[k].status[j + sum_].nonevents[i]));
	      }
	      outdata.push_back(std::to_string(groupsStats[k].status[j].events[i] + groupsStats[k].status[j + sum_].events[i]));

	      if(printGroupsTotalAges == true && k == groupsAges.size() - 1){
		if(printGroupsPop == true){
		  outdata.push_back(std::to_string(groupsTotalAgeStats.status[j].nonevents[i] + groupsTotalAgeStats.status[j + sum_].nonevents[i]));
		}
		outdata.push_back(std::to_string(groupsTotalAgeStats.status[j].events[i] + groupsTotalAgeStats.status[j + sum_].events[i]));
	      }
	    }
	  }
	}
      }
    }
  }else{
    outdata.clear();
    for(int i = 0; i < 5 ; i++){
      if(groupsEvents[i] == 1){
	for(unsigned k = 0; k < groupsStats.size(); k++){
	  if(printGroupsPop == true){
	    outdata.push_back(std::to_string(groupsStats[k].total.nonevents[i]));
	  }
	  outdata.push_back(std::to_string(groupsStats[k].total.events[i]));
	  if(printGroupsTotalAges == true && k == groupsAges.size() - 1){
	    if(printGroupsPop == true){
	      outdata.push_back(std::to_string(groupsTotalAgeStats.total.nonevents[i]));
	    }
	    outdata.push_back(std::to_string(groupsTotalAgeStats.total.events[i]));
	  }
	}
      }
    }
  }
  outdata.push_back(std::to_string(truePositives));
  outdata.push_back(std::to_string(trueNegatives));
  outdata.push_back(std::to_string(falsePositives));
  outdata.push_back(std::to_string(falseNegatives));
  Report::join(outdata,',',outstring);
  outGroups << outstring;
}



void Report::printCohortReport(int currDay){
  outCohort << currDay << ",";

  for(int i = 0; i < 5 ; i++){
    if(cohortEvents[i] == 1){
      if(cohortStatus[0] > 0 || cohortStatus[1] > 0){
	if(cohortStatus[0] > 0 && cohortStatus[1] > 0){
	  for(int j = 0;j < 4 ;j++){
	    for(unsigned k = 0; k < cohortAges.size(); k++){
	      if(printCohortPop == true){
		outCohort << cohortStats[k].status[j].nonevents[i]<<",";
	      }
	      outCohort << cohortStats[k].status[j].events[i];
	      if(i == cohortMaxIndex && k == cohortAges.size() - 1 && j == 3){
		outCohort << "\n";
	      }else{
		outCohort << ",";
	      }
	    }
	  }
	}else{
	  int inc_ = 0;
	  int sum_ = 0;
	  if(cohortStatus[0] > 0){
	    inc_ = 2;
	    sum_ = 1;
	  }else if(cohortStatus[1] > 0){
	    inc_ = 1;
	    sum_ = 2;
	  }
	  // ?? CHECK j = j = j
	  for(int j = 0;j < 2 * inc_;j = sum_ + inc_){
	    for(unsigned k = 0; k < cohortAges.size(); k++){
	      if(printCohortPop == true){
		outCohort << cohortStats[k].status[j].nonevents[i] + cohortStats[k].status[j + sum_].nonevents[i]<<",";
	      }
	      outCohort << cohortStats[k].status[j].events[i] + cohortStats[k].status[j + sum_].events[i];
	      if(i == cohortMaxIndex && k == cohortAges.size() - 1 && j == inc_){
		outCohort << "\n";
	      }else{
		outCohort << ",";
	      }
	    }
	  }
	}

      }else{
	for(unsigned k = 0; k < cohortStats.size(); k++){
	  if(printCohortPop == true){
	    outCohort << cohortStats[k].total.nonevents[i]<<",";
	  }
	  outCohort << cohortStats[k].total.events[i];
	  if(i == cohortMaxIndex && k == cohortStats.size() - 1){
	    outCohort << "\n";
	  }else{
	    outCohort << ",";
	  }
	}
      }
    }
  }
}



void Report::printHeaders(){
  if(reportGroups == true){
    printGroupsHeader();
  }
  if(reportCohort == true){
    printCohortHeader();
  }
  if(reportAges == true){
    printAgesHeader();
  }
  if(reportSpatial == true){
    printSpatialHeader();
  }
}

void Report::join(const vector<string>& v, char c, string& s) {

  s.clear();

  for (vector<string>::const_iterator p = v.begin(); p != v.end(); ++p) {
    s += *p;
    if (p != v.end() - 1){
      s += c;
    }else{
      s +="\n";
    }
  }

}

void Report::printSpatialHeader(){
  outSpatial << "Variable,code,Day,Serotype\n";
}


void Report::printAgesHeader(){
  outAges << "day,";
  for(int i = 0; i < 5 ; i++){
    if(ageEvents[i] == 1){
      if(ageStatus[0] > 0 || ageStatus[1] > 0){
	if(ageStatus[0] > 0 && ageStatus[1] > 0){
	  for(int j = 0;j < 2 ;j++){
	    for(int jj = 2; jj < 4; jj++){
	      for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
		if(printAgesPop == true){
		  outAges << status[j].c_str() <<"_"<< status[jj].c_str() << "_age_" << k + discreteAges.min << "_";
		  outAges << "no"<<events[i].c_str()<<",";
		}
		outAges << status[j].c_str() << "_" << status[jj].c_str() << "_age_" << k + discreteAges.min << "_";
		outAges << events[i].c_str();
		if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min && j == 1 && jj == 3){
		  outAges << "\n";
		}else{
		  outAges << ",";
		}
	      }
	    }
	  }
	}else{
	  int ind_ = 0;
	  if(ageStatus[0] > 0){
	    ind_ = 0;
	  }else if(ageStatus[1] > 0){
	    ind_ = 2;
	  }
	  for(int j = ind_;j < ind_ + 2;j++){
	    for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
	      if(printAgesPop == true){
		outAges << status[j].c_str() << "_age_" << k + discreteAges.min << "_";
		outAges << "no"<<events[i].c_str()<<",";
	      }
	      outAges << status[j].c_str() << "_age_" << k + discreteAges.min << "_";
	      outAges << events[i].c_str();
	      if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min && j == ind_ + 1){
		outAges << "\n";
	      }else{
		outAges << ",";
	      }
	    }
	  }
	}

      }else{
	for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
	  if(printAgesPop == true){
	    outAges << "age_" << k + discreteAges.min << "_";
	    outAges << "no"<<events[i].c_str()<<",";
	  }
	  outAges << "age_" << k + discreteAges.min << "_";
	  outAges << events[i].c_str();
	  if(i == ageMaxIndex && k == discreteAges.max - discreteAges.min){
	    outAges << "\n";
	  }else{
	    outAges << ",";
	  }
	}
      }
    }
  }
}

void Report::printGroupsHeader(){
  outGroups << "day,";
  if(printGroupsAgeFirst == true){
    outGroups << "avg_age_first,";
  }
  for(int i = 0; i < 5 ; i++){
    if(groupsEvents[i] == 1){
      if(groupsStatus[0] > 0 || groupsStatus[1] > 0){
	if(groupsStatus[0] > 0 && groupsStatus[1] > 0){
	  for(int j = 0;j < 2 ;j++){
	    for(int jj = 2; jj < 4; jj++){
	      for(unsigned k = 0; k < groupsAges.size(); k++){
		if(printGroupsPop == true){
		  outGroups << status[j].c_str() <<"_"<< status[jj].c_str() << "_age_" << groupsAges[k].min << "_" << groupsAges[k].max << "_";
		  outGroups << "no"<<events[i].c_str()<<",";
		}
		outGroups << status[j].c_str() << "_" << status[jj].c_str() << "_age_" << groupsAges[k].min << "_" << groupsAges[k].max << "_";
		outGroups << events[i].c_str();

		if(printGroupsTotalAges == true && k == groupsAges.size() - 1){
		  if(printGroupsPop == true){
		    outGroups << "," << status[j].c_str() <<"_"<< status[jj].c_str() << "_all_";
		    outGroups << "no"<<events[i].c_str();
		  }
		  outGroups << "," << status[j].c_str() << "_" << status[jj].c_str() << "_all_";
		  outGroups << events[i].c_str();
		}
		if(i == groupsMaxIndex && k == groupsAges.size() - 1 && j == 1 && jj == 3){
		  outGroups << ",truePositives,trueNegatives,falsePositives,falseNegatives\n";
		}else{
		  outGroups << ",";
		}
	      }
	    }
	  }
	}else{
	  int ind_ = 0;
	  if(groupsStatus[0] > 0){
	    ind_ = 0;
	  }else if(groupsStatus[1] > 0){
	    ind_ = 2;
	  }
	  for(int j = ind_;j < ind_ + 2;j++){
	    for(unsigned k = 0; k < groupsAges.size(); k++){
	      if(printGroupsPop == true){
		outGroups << status[j].c_str() << "_age_" << groupsAges[k].min << "_" << groupsAges[k].max << "_";
		outGroups << "no"<<events[i].c_str()<<",";
	      }
	      outGroups << status[j].c_str() << "_age_" << groupsAges[k].min << "_" << groupsAges[k].max << "_";
	      outGroups << events[i].c_str();

	      if(printGroupsTotalAges == true && k == groupsAges.size() - 1){
		if(printGroupsPop == true){
		  outGroups << "," << status[j].c_str() << "_all_";
		  outGroups << "no"<<events[i].c_str();
		}
		outGroups << "," << status[j].c_str() << "_all_";
		outGroups << events[i].c_str();
	      }
	      if(i == groupsMaxIndex && k == groupsAges.size() - 1 && j == ind_ + 1){
		outGroups << ",truePositives,trueNegatives,falsePositives,falseNegatives\n";
	      }else{
		outGroups << ",";
	      }
	    }
	  }
	}

      }else{
	for(unsigned k = 0; k < groupsAges.size(); k++){
	  if(printGroupsPop == true){
	    outGroups << "age_" << groupsAges[k].min << "_" << groupsAges[k].max << "_";
	    outGroups << "no"<<events[i].c_str()<<",";
	  }
	  outGroups << "age_" << groupsAges[k].min << "_" << groupsAges[k].max << "_";
	  outGroups << events[i].c_str();

	  if(printGroupsTotalAges == true && k == groupsAges.size() - 1){
	    if(printGroupsPop == true){
	      outGroups << "," <<  "all_" ;
	      outGroups << "no"<<events[i].c_str();
	    }
	    outGroups << "," << "all_";
	    outGroups << events[i].c_str();
	  }
	  if(i == groupsMaxIndex && k == groupsAges.size() - 1){
	    outGroups << ",truePositives,trueNegatives,falsePositives,falseNegatives\n";
	  }else{
	    outGroups << ",";
	  }
	}
      }
    }
  }
}



void Report::printCohortHeader(){
  outCohort << "day,";
  for(int i = 0; i < 5 ; i++){
    if(cohortEvents[i] == 1){
      if(cohortStatus[0] > 0 || cohortStatus[1] > 0){
	if(cohortStatus[0] > 0 && cohortStatus[1] > 0){
	  for(int j = 0;j < 2 ;j++){
	    for(int jj = 2; jj < 4; jj++){
	      for(unsigned k = 0; k < cohortAges.size(); k++){
		if(printCohortPop == true){
		  outCohort << status[j].c_str() <<"_"<< status[jj].c_str() << "_age_" << cohortAges[k].min << "_" << cohortAges[k].max << "_";
		  outCohort << "no"<<events[i].c_str()<<",";
		}
		outCohort << status[j].c_str() << "_" << status[jj].c_str() << "_age_" << cohortAges[k].min << "_" << cohortAges[k].max << "_";
		outCohort << events[i].c_str();
		if(i == cohortMaxIndex && k == cohortAges.size() - 1 && j == 1 && jj == 3){
		  outCohort << "\n";
		}else{
		  outCohort << ",";
		}
	      }
	    }
	  }
	}else{
	  int ind_ = 0;
	  if(cohortStatus[0] > 0){
	    ind_ = 0;
	  }else if(cohortStatus[1] > 0){
	    ind_ = 2;
	  }
	  for(int j = ind_;j < ind_ + 2;j++){
	    for(unsigned k = 0; k < cohortAges.size(); k++){
	      if(printCohortPop == true){
		outCohort << status[j].c_str() << "_age_" << cohortAges[k].min << "_" << cohortAges[k].max << "_";
		outCohort << "no"<<events[i].c_str()<<",";
	      }
	      outCohort << status[j].c_str() << "_age_" << cohortAges[k].min << "_" << cohortAges[k].max << "_";
	      outCohort << events[i].c_str();
	      if(i == cohortMaxIndex && k == cohortAges.size() - 1 && j == ind_ + 1){
		outCohort << "\n";
	      }else{
		outCohort << ",";
	      }
	    }
	  }
	}

      }else{
	for(unsigned k = 0; k < cohortAges.size(); k++){
	  if(printCohortPop == true){
	    outCohort << "age_" << cohortAges[k].min << "_" << cohortAges[k].max << "_";
	    outCohort << "no"<<events[i].c_str()<<",";
	  }
	  outCohort << "age_" << cohortAges[k].min << "_" << cohortAges[k].max << "_";
	  outCohort << events[i].c_str();
	  if(i == cohortMaxIndex && k == cohortAges.size() - 1){
	    outCohort << "\n";
	  }else{
	    outCohort << ",";
	  }
	}
      }
    }
  }
}

void Report::resetReports(){
  if(reportGroups == true){
    resetGroupStats();
  }
  if(reportCohort == true){
    resetCohortStats();
  }
  if(reportAges == true){
    resetAgeStats();
  }
  if(reportSpatial == true){
    resetSpatialStats();
  }
}

void Report::resetSpatialStats(){
  spatialData.clear();
}


void Report::resetAgeStats(){
  ageStats.clear();
  if(discreteAges.min + discreteAges.max <= 0 || discreteAges.min >= discreteAges.max || discreteAges.min < 0){
    exit(1);
  }

  for(int k = 0; k <= discreteAges.max - discreteAges.min; k++){
    reportStats tempStats;
    for(int i = 0; i < 5; i++){
      // There are four status VacSero+ VacSero- PlacSero+ PlacSero-
      for(int j = 0; j < 4; j++){
	tempStats.status[j].events[i] = 0;
	tempStats.status[j].nonevents[i] = 0;
      }
      tempStats.total.events[i] = 0;
      tempStats.total.nonevents[i] = 0;
    }
    ageStats.push_back(tempStats);
  }

  for(int i = 0; i < 5; i++){
    if(ageEvents[i] == 1){
      ageMaxIndex = i;
    }
  }
}

void Report::resetCohortStats(){
  cohortStats.clear();
  if(cohortAges.empty()){
    exit(1);
  }

  for(unsigned k = 0; k < cohortAges.size(); k++){
    reportStats tempStats;
    for(int i = 0; i < 5; i++){
      // There are four status VacSero+ VacSero- PlacSero+ PlacSero-
      for(int j = 0; j < 4; j++){
	tempStats.status[j].events[i] = 0;
	tempStats.status[j].nonevents[i] = 0;
      }
      tempStats.total.events[i] = 0;
      tempStats.total.nonevents[i] = 0;
    }
    cohortStats.push_back(tempStats);
  }

  for(int i = 0; i < 5; i++){
    if(cohortEvents[i] == 1){
      cohortMaxIndex = i;
    }
  }
}


void Report::resetGroupStats(){
  groupsAvgAgeFirst = 0;
  groupsTotalFirstInf = 0;
  groupsStats.clear();
  falsePositives = 0;
  falseNegatives = 0;
  truePositives = 0;
  trueNegatives = 0;
  if(groupsAges.empty()){
    exit(1);
  }
  for(unsigned k = 0; k < groupsAges.size(); k++){
    reportStats tempStats;
    for(int i = 0; i < 5; i++){
      // There are four status VacSero+ VacSero- PlacSero+ PlacSero-
      for(int j = 0; j < 4; j++){
	tempStats.status[j].events[i] = 0;
	tempStats.status[j].nonevents[i] = 0;
      }
      tempStats.total.events[i] = 0;
      tempStats.total.nonevents[i] = 0;
    }
    groupsStats.push_back(tempStats);
  }

  for(int i = 0; i < 5; i++){
    if(groupsEvents[i] == 1){
      groupsMaxIndex = i;
    }
  }

  for(int i = 0; i < 5; i++){
    // There are four status VacSero+ VacSero- PlacSero+ PlacSero-
    for(int j = 0; j < 4; j++){
      groupsTotalAgeStats.status[j].events[i] = 0;
      groupsTotalAgeStats.status[j].nonevents[i] = 0;
    }
    groupsTotalAgeStats.total.events[i] = 0;
    groupsTotalAgeStats.total.nonevents[i] = 0;
  }
}


void Report::updateFOI(int currDay, Human * h){
  human_counts[currDay]++;
  string tmpzone = h->getZoneID();
  if(h->isInfected()){
    if(h->infection != nullptr){
      if(h->infection->getStartDay() == currDay){
	int sero = h->infection->getInfectionType();
	if(sero > 0){
	  dailyFOI[currDay][sero]["newinf"]++;
	  dailyFOI[currDay][sero][tmpzone + "newinf"]++;
	}
      }
      if(h->isSymptomatic() && floor(h->infection->getSymptomOnset()) == currDay){
	int sero = h->infection->getInfectionType();
	if(sero > 0){
	  dailyFOI[currDay][sero]["newsymp"]++;
	}
      }
    }
  }
  // Check for immunity against all the serotypes
  for(unsigned i = 1; i < 5; i++){
    if(!h->isPermImmune(i)){
      dailyFOI[currDay][i]["sus"]++;
      dailyFOI[currDay][i][tmpzone + "sus"]++;
    }
    if(!h->isImmune(i)){
      dailyFOI[currDay][i]["sustmp"]++;
    }
  }
}

void Report::updateSecondaryCases(int currDay, Human * h){
  if(h != nullptr){
    if(h->isInfected()){
      if(h->infection != nullptr){
	unsigned sero = h->infection->getInfectionType();
	int infDay = h->infection->getStartDay();
	string id = h->getPersonID();
	secondaryCases[infDay][sero][id] = h->getR0(sero);
      }
    }
    if(h->isImported()){
      if(h->infection != nullptr){
	unsigned sero = h->infection->getInfectionType();
	int infDay = h->infection->getStartDay();
	string id = h->getVisitorID();
	secondaryCases[infDay][sero][id] = h->getVisitorR0(sero);
      }
    }
  }
}

void Report::addImportation(int currDay, int sero_in, Human * h){
  if(sero_in > 0 && sero_in < 5){
    dailyFOI[currDay][sero_in]["imports"]++;
  }
}

void Report::printFOIReport(int lastDay){
  outFOI.open(outputFOIFile);
  if (!outFOI.good()) {
    exit(1);
  }
  std::vector<string> headerstr; string headout;
  headerstr.clear();
  for(int s = 1; s < 5; s++){
    if(foiTypes[s - 1]){
      if(printR0 == true){
	headerstr.push_back("R0_Denv"+std::to_string(s));
      }
      headerstr.push_back("Denv"+std::to_string(s));
      headerstr.push_back("Susceptible_" + std::to_string(s));
      headerstr.push_back("Susceptible_temp_" + std::to_string(s));
      headerstr.push_back("Infectious_" + std::to_string(s));
      headerstr.push_back("Symptomatic_" + std::to_string(s));
      headerstr.push_back("MozSusceptible_" + std::to_string(s));
      headerstr.push_back("MozExposed_" + std::to_string(s));
      headerstr.push_back("MozInfectious_" + std::to_string(s));
      headerstr.push_back("Importations_" + std::to_string(s));
    }
    if(printZonesFOI == true){
      for(auto locIt = zonesToPrint.begin(); locIt != zonesToPrint.end();){
	headerstr.push_back("FOI_" + *locIt + std::to_string(s));
	++locIt;
      }
    }
  }
  if(reportFOIOutbreakSymptomatics){
    headerstr.push_back("ReportedSymptomatics");
    headerstr.push_back("OutbreakResponseLocations");
  }
  headerstr.push_back("Humans");

  if(!headerstr.empty()){
    outFOI << "day,";
    Report::join(headerstr,',',headout);
    outFOI << headout;
  }

  for(int d = 0; d < lastDay; d++){
    string outstring;
    vector<string> foi_values; foi_values.clear();
    for(unsigned s = 1; s < 5; s++){
      if(foiTypes[s - 1] != 1){
	continue;
      }
      if(printR0 == true){
	string R0str = "0.000";
	if(secondaryCases.find(d) != secondaryCases.end()){
	  if(secondaryCases[d].find(s) != secondaryCases[d].end()){
	    int dailyR0 = 0;
	    int sumInfectors = 0;
	    for(auto idt = secondaryCases[d][s].begin(); idt != secondaryCases[d][s].end();idt++){
	      if(idt->second >= 0){
		dailyR0 += idt->second;
		sumInfectors++;
	      }
	    }
	    R0str = sumInfectors > 0 ? std::to_string((double) dailyR0 / (double) sumInfectors) : 0;
	  }
	}
	foi_values.push_back(R0str);
      }

      string foistr = "0.000";
      if(dailyFOI.find(d) != dailyFOI.end()){
	if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	  double foi_temp = dailyFOI[d][s]["sustmp"] > 0 ? (double) dailyFOI[d][s]["newinf"] / (double) dailyFOI[d][s]["sustmp"] : 0;
	  foistr = std::to_string(foi_temp);
	}
      }
      foi_values.push_back(foistr);

      foistr = "0.000";
      if(dailyFOI.find(d) != dailyFOI.end()){
	if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	  foistr = std::to_string(dailyFOI[d][s]["sus"]);
	}
      }
      foi_values.push_back(foistr);

      foistr = "0.000";
      if(dailyFOI.find(d) != dailyFOI.end()){
	if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	  foistr = std::to_string(dailyFOI[d][s]["sustmp"]);
	}
      }
      foi_values.push_back(foistr);

      foistr = "0.000";
      if(dailyFOI.find(d) != dailyFOI.end()){
	if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	  foistr = std::to_string(dailyFOI[d][s]["newinf"]);
	}
      }
      foi_values.push_back(foistr);

      foistr = "0.000";
      if(dailyFOI.find(d) != dailyFOI.end()){
	if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	  foistr = std::to_string(dailyFOI[d][s]["newsymp"]);
	}
      }
      foi_values.push_back(foistr);

      foistr = "0.000";
      if(dailyFOI.find(d) != dailyFOI.end()){
	if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	  foistr = std::to_string(dailyFOI[d][s]["mozsus"]);
	}
      }
      foi_values.push_back(foistr);

      foistr = "0.000";
      if(dailyFOI.find(d) != dailyFOI.end()){
	if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	  foistr = std::to_string(dailyFOI[d][s]["mozexp"]);
	}
      }
      foi_values.push_back(foistr);

      foistr = "0.000";
      if(dailyFOI.find(d) != dailyFOI.end()){
	if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	  foistr = std::to_string(dailyFOI[d][s]["mozinf"]);
	}
      }
      foi_values.push_back(foistr);

      foistr = "0.000";
      if(dailyFOI.find(d) != dailyFOI.end()){
	if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	  foistr = std::to_string(dailyFOI[d][s]["imports"]);
	}
      }
      foi_values.push_back(foistr);

      if(printZonesFOI == true){
	for(auto locIt = zonesToPrint.begin(); locIt != zonesToPrint.end();){
	  foistr = "0.000";
	  if(dailyFOI[d].find(s) != dailyFOI[d].end()){
	    if(dailyFOI[d][s].find(*locIt + "newinf") != dailyFOI[d][s].end()){
	      double foi_temp = dailyFOI[d][s][*locIt + "sus"] > 0 ? (double) dailyFOI[d][s][*locIt + "newinf"] / (double) dailyFOI[d][s][*locIt + "sus"] : 0;
	      foistr = std::to_string(foi_temp);
	    }
	  }
	  ++locIt;
	  foi_values.push_back(foistr);
	}
      }
    }
    if(reportFOIOutbreakSymptomatics){
      string outbreakstr= "0";
      if(dailyOutbreakResponseSymp.find(d) != dailyOutbreakResponseSymp.end()){
	outbreakstr = std::to_string(dailyOutbreakResponseSymp[d]);
      }
      foi_values.push_back(outbreakstr);
      outbreakstr= "0";

      if(dailyOutbreakResponseLocations.find(d) != dailyOutbreakResponseLocations.end()){
	outbreakstr = std::to_string(dailyOutbreakResponseLocations[d]);
      }
      foi_values.push_back(outbreakstr);
    }

    string humstr = "0";
    if(human_counts.find(d) != human_counts.end()){
      humstr = std::to_string(human_counts[d]);
    }
    foi_values.push_back(humstr);

    if(!foi_values.empty()){
      Report::join(foi_values,',',outstring);
      outFOI << d << ",";
      outFOI << outstring;
    }
  }
  outFOI.close();
}

void Report::printImmatureReport(int lastDay){
  outImmature.open(outputImmatureFile);
  if (!outImmature.good()) {
    exit(1);
  }
  std::vector<string> headerstr; string headout;
  headerstr.clear();
  headerstr.push_back("Eggs");
  headerstr.push_back("Larvae");
  headerstr.push_back("Pupae");
  headerstr.push_back("Adults");
  if (reportZeroMozLocations) {
    headerstr.push_back("Zero Locations");
  }
  if(!headerstr.empty()){
    outImmature << "day,";
    Report::join(headerstr,',',headout);
    outImmature << headout;
  }

  for(int d = 0; d < lastDay; d++){
    string outstring;
    vector<string> immature_values; immature_values.clear();
    
    string immstr = "0.000";
    if(dailyImmature.find(d) != dailyImmature.end()){
      immstr = std::to_string(dailyImmature[d]["eggs"]);
    }
    immature_values.push_back(immstr);

    immstr = "0.000";
    if(dailyImmature.find(d) != dailyImmature.end()){
      immstr = std::to_string(dailyImmature[d]["larvae"]);
    }
    immature_values.push_back(immstr);

    immstr = "0.000";
    if(dailyImmature.find(d) != dailyImmature.end()){
      immstr = std::to_string(dailyImmature[d]["pupae"]);
    }
    immature_values.push_back(immstr);

    string mozstr = "0";
    if(mosquito_counts.find(d) != mosquito_counts.end()){
      mozstr = std::to_string(mosquito_counts[d]);
    }
    immature_values.push_back(mozstr);

    if (reportZeroMozLocations) {
      string locstr = "0";
      if(dailyZeroMozLocations.find(d) != dailyZeroMozLocations.end()){
	locstr = std::to_string(dailyZeroMozLocations[d]);
      }
      immature_values.push_back(locstr);
    }
    
    if(!immature_values.empty()){
      Report::join(immature_values,',',outstring);
      outImmature << d << ",";
      outImmature << outstring;
    }
  }
  outImmature.close();
}

void Report::printMozAgesReport(int lastDay){
  outMozAges.open(outputMozAgesFile);
  if (!outMozAges.good()) {
    exit(1);
  }
  std::vector<string> headerstr;
  string headout;
  headerstr.clear();
  headerstr.push_back("Age");
  for (int d=0; d<lastDay; d++) {
    headerstr.push_back("Day " + std::to_string(d));
  }
  if(!headerstr.empty()){
    Report::join(headerstr,',',headout);
    outMozAges << headout;
  }
  for (auto it = mosquito_ages.begin(); it != mosquito_ages.end(); ++it) {
    string outstring;
    vector<string> moz_age_values; moz_age_values.clear();
    moz_age_values.push_back(std::to_string(it->first));
    for(int d = 0; d < lastDay; d++){
      string mozagestr = "0.000";
      if(it->second.find(d) != it->second.end()){
	mozagestr = std::to_string(it->second[d]);
      }
      moz_age_values.push_back(mozagestr);
    }
    if(!moz_age_values.empty()){
      Report::join(moz_age_values,',',outstring);
      outMozAges << outstring;
    }
  }
  outMozAges.close();
}

void Report::printSprayReport(int lastDay){
  outSpray.open(outputSprayFile);
  if (!outSpray.good()) {
    exit(1);
  }
  std::vector<string> headerstr; string headout;
  headerstr.clear();
  if (reportSpray) headerstr.push_back("# sprayed");
  if (reportSprayLocs) headerstr.push_back("locIDs");
  if(!headerstr.empty()){
    outSpray << "day,";
    Report::join(headerstr,',',headout);
    outSpray << headout;
  }

  for(int d = 0; d < lastDay; d++){
    string outstring;
    vector<string> spray_values; spray_values.clear();
    
    string sprstr = "0";
    if(spray_counts.find(d) != spray_counts.end()){
      sprstr = std::to_string(spray_counts[d]);
    }
    spray_values.push_back(sprstr);

    if(dailySprayLocs.find(d) != dailySprayLocs.end()){
      for (auto & loc: dailySprayLocs[d]){
	spray_values.push_back(loc);
      }
    }
    
    if(!spray_values.empty()){
      Report::join(spray_values,',',outstring);
      outSpray << d << ",";
      outSpray << outstring;
    }
  }
  outSpray.close();
}


void Report::finalizeReport(int currDay){
  outCohort.close();
  outAges.close();
  outGroups.close();
  outSpatial.close();
  if(reportFOI == true){
    this->printFOIReport(currDay);
  }
  if (reportImmature) this->printImmatureReport(currDay);
  if (reportMozAges) this->printMozAgesReport(currDay);
  if (reportSpray || reportSprayLocs) this->printSprayReport(currDay);
}



//Report::~Report() {
//}

