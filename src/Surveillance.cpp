#include <fstream>
#include <string>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <sstream>
#include "Surveillance.h"


Surveillance::Surveillance(){
    contactFrequency = 0;
    firstContactDelay = 0;
    selfReportProb = 0.0;
    reportTodayProb = 0.0;
    printExposure = false;
    recordsDatabase.clear();
    parameters.clear();
}

void Surveillance::setup(string file){
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

    contactFrequency = this->readParameter("surveillance_contact_frequency",0);
    firstContactDelay = this->readParameter("surveillance_first_contact_delay",0);
    selfReportProb = this->readParameter("surveillance_self_report_probability",0.0);
    printExposure = this->readParameter("surveillance_print_exposure",printExposure);
    double avgDelay = this->readParameter("surveillance_avg_report_delay",0.0);

    reportTodayProb = avgDelay > 0.0 ? (double) 1.0 / avgDelay : 1.0;
    reportTodayProb = reportTodayProb > 1.0 ? 1.0 : reportTodayProb;
    printf("Surveillance setup contact frecuency %d first delay %d selfReport probability %.2f prob report %.6f\n", contactFrequency, firstContactDelay, selfReportProb, reportTodayProb);
    if(printExposure){
	printf("Printing specific date of exposure\n");
    }
}

void Surveillance::initialize_human_surveillance(Human * h, int currDay){
    h->setFirstContactWithTrial(currDay + firstContactDelay);
    h->setContactByTrial(currDay + firstContactDelay);
    h->setSeroStatusAtVaccination();
    h->setSelfReportProb(selfReportProb);
    string id(h->getPersonID());
    hRecord tempRecord;
    tempRecord.ageDaysAtVaccination = h->getAgeTrialEnrollment();
    tempRecord.seroStatusAtVaccination = h->getSeroStatusAtVaccination();
    tempRecord.trialArm = h->getTrialArm();
    tempRecord.houseID = h->getHouseID();
    tempRecord.houseMemNum = h->getHouseMemNum();
    tempRecord.dropoutDay = -1;
    tempRecord.enrollmentDay = currDay;
    tempRecord.firstTTL = currDay;
    tempRecord.firstTTR = -1;
    tempRecord.firstExp = 0;
    tempRecord.firstPCR = "NA";
    tempRecord.pcr.clear();
    tempRecord.primary.clear();
    tempRecord.dateExposure.clear();
    tempRecord.dropout = false;
    tempRecord.bites = 0;
    tempRecord.attractiveness = 0.0;
    tempRecord.first_real_infected = -1;
    tempRecord.first_real_symptomatic = -1;
    for(int i = 0;i < 4; i++){
	tempRecord.TTL[i] = currDay;
	tempRecord.TTR[i] = -1;
	tempRecord.onset[i] = -1.0;
	tempRecord.symptoms[i] = -1;
	tempRecord.hosp[i] = -1;
	tempRecord.real_infected_day[i] = -1;
	tempRecord.real_symptomatic_day[i] = -1;
	tempRecord.numExp[i] = 0;
	tempRecord.previousExposure[i] = 0;       
	tempRecord.pcr.push_back("NA");
	tempRecord.pcrDay[i] = -1;
	tempRecord.primary.push_back("NA");
	tempRecord.dateExposure.push_back("NA");
    }
    recordsDatabase.insert(make_pair(id,tempRecord));
}
void Surveillance::track_infected(Human * h, int currDay){
    string id(h->getPersonID());
    //    printf("track_infected %s day %d\n", id.c_str(),currDay);
    if(h->isInfected() && h->infection != nullptr){	
	if(recordsDatabase.find(id)->second.first_real_infected == -1){
	    //	    printf("FIRST INFECTION IN THE TRIAL: %s day %d\n", id.c_str(),currDay);
	    recordsDatabase.find(id)->second.first_real_infected = h->infection->getStartDay();
	    recordsDatabase.find(id)->second.first_real_symptomatic = (h->isSymptomatic() == true) ? h->infection->getSymptomOnset() : -1;
	    //	    if(recordsDatabase.find(id)->second.first_real_symptomatic >= 0){
		//		printf("FIRST SYMPTOMATIC IN THE TRIAL: %s day %d onset %d\n", id.c_str(),currDay,recordsDatabase.find(id)->second.first_real_symptomatic);
	    //	    }
	}
	unsigned serotype = h->infection->getInfectionType() - 1;
	if(recordsDatabase.find(id)->second.real_infected_day[serotype] == -1){
	    //	    printf("FIRST INFECTION IN THE TRIAL BY SEROTYPE %d %s day %d\n", serotype + 1, id.c_str(),currDay);
	    recordsDatabase.find(id)->second.real_infected_day[serotype] = h->infection->getStartDay();
	    recordsDatabase.find(id)->second.real_symptomatic_day[serotype] = (h->isSymptomatic() == true) ? h->infection->getSymptomOnset() : -1;
	}
    }
}

int Surveillance::update_human_surveillance(Human * h, int currDay, RandomNumGenerator * rGen){
    int pcr_result = -1;
    string id(h->getPersonID());
    if(h->getFirstContactWithTrial() > currDay){
	return pcr_result;
    }
    //    printf("update surveillance %s day %d\n", id.c_str(), currDay);
    if(recordsDatabase.find(id) != recordsDatabase.end() && h->isEnrolledInTrial()){
	if( ( currDay - recordsDatabase.find(id)->second.enrollmentDay ) >= 30){
	    recordsDatabase.find(id)->second.firstExp = 0;
	    recordsDatabase.find(id)->second.bites = h->getNumBites();
	    recordsDatabase.find(id)->second.attractiveness = h->getAttractiveness();	    
	    for(unsigned i = 0; i < 4; i++){
		recordsDatabase.find(id)->second.numExp[i] = h->getExposedCount(i);
		recordsDatabase.find(id)->second.firstExp += h->getExposedCount(i);
		recordsDatabase.find(id)->second.previousExposure[i] = h->getPreExposureAtVaccination(i);
		recordsDatabase.find(id)->second.dateExposure[i] = h->getExposureDate(i);
	    }
	    if(h->isInfected() && h->infection != nullptr){
		// set IIP and serotype and symptoms
		unsigned serotype = h->infection->getInfectionType() - 1;
		recordsDatabase.find(id)->second.onset[serotype] = h->infection->getSymptomOnset();
		recordsDatabase.find(id)->second.symptoms[serotype] = (h->isSymptomatic() == true) ? 1 : 0;
		recordsDatabase.find(id)->second.hosp[serotype] = (h->isHospitalized() == true) ? 1 : 0;
		recordsDatabase.find(id)->second.primary[serotype] = (h->infection->isPrimary()) ? "primary" : "secondary";
		
		if(h->getReportSymptoms()){
		    if(currDay > recordsDatabase.find(id)->second.onset[serotype] + 2 && currDay <= recordsDatabase.find(id)->second.onset[serotype] + 7){
			if(recordsDatabase.find(id)->second.TTR[h->infection->getInfectionType() - 1] == -1){
			    if(rGen->getEventProbability() < reportTodayProb){
				pcr_result = this->PCR_test(h,currDay,rGen);
				recordsDatabase.find(id)->second.pcrDay[serotype] = currDay;
				if(pcr_result >= 0){
				    recordsDatabase.find(id)->second.TTR[pcr_result] = currDay;
				    recordsDatabase.find(id)->second.pcr[pcr_result] = "POSITIVE"; 
				    if(recordsDatabase.find(id)->second.firstTTR == -1){
					recordsDatabase.find(id)->second.firstTTR = currDay;
					recordsDatabase.find(id)->second.firstPCR = "POSITIVE";
				    }
				}
			    }
			    //Bookkeeping when reporting
			    for(int i = 0; i < 4; i++){
				if(recordsDatabase.find(id)->second.TTR[i] == -1){
				    recordsDatabase.find(id)->second.TTL[i] = currDay;
				}
			    }
			    if(recordsDatabase.find(id)->second.firstTTR == -1){
				recordsDatabase.find(id)->second.firstTTL = currDay;
			    }
			    h->setContactByTrial(currDay);
			}
		    }else if(currDay > recordsDatabase.find(id)->second.onset[serotype] + 7){
			h->setReportSymptoms(false);
		    }
		}
	    }
	}

	if(h->getLastContactByTrial() + contactFrequency == currDay){
	    //	    printf("Contacting person\n");
	    this->contactPerson(h, currDay, rGen);
	    //Bookkeeping for each contact 
	    for(int i = 0; i < 4; i++){
		if(recordsDatabase.find(id)->second.TTR[i] == -1){
		    recordsDatabase.find(id)->second.TTL[i] = currDay;
		}
	    }
	    if(recordsDatabase.find(id)->second.firstTTR == -1){
		recordsDatabase.find(id)->second.firstTTL = currDay;
	    }
	    h->setContactByTrial(currDay);
	}

	//annual visits -> detect one or more infections in a small sample -- TO BE IMPLEMENTED !!!!
	
    }
    return pcr_result;
}

void Surveillance::finalize_human_surveillance(Human *h, int currDay, bool drop_in){
    string id(h->getPersonID());
    recordsDatabase.find(id)->second.dropout = drop_in;
    recordsDatabase.find(id)->second.dropoutDay = currDay;
    /*      
    recordsDatabase.find(id)->second.firstExp = 0;
    for(unsigned i = 0; i < 4; i++){
	recordsDatabase.find(id)->second.numExp[i] = h->getExposedCount(i);
	recordsDatabase.find(id)->second.firstExp += h->getExposedCount(i);
	recordsDatabase.find(id)->second.previousExposure[i] = h->getPreExposureAtVaccination(i);
	recordsDatabase.find(id)->second.dateExposure[i] = h->getExposureDate(i);
    }
    */
}

void Surveillance::contactPerson(Human * h, int currDay, RandomNumGenerator * rGen){
    string id(h->getPersonID());
    if(h->isInfected() && h->infection != nullptr){
	int serotype =  h->infection->getInfectionType() -1;
	if(recordsDatabase.find(id)->second.TTR[serotype] == -1){
	    // If the symptoms happened between last contact and today
	    if(recordsDatabase.find(id)->second.symptoms[serotype] > 0 && recordsDatabase.find(id)->second.onset[serotype] <= currDay && recordsDatabase.find(id)->second.onset[serotype] > h->getLastContactByTrial()){
		h->setReportSymptoms(true);
	    }
	}
    }
}

int Surveillance::PCR_test(Human * h, int currDay, RandomNumGenerator * rGen){
    string id(h->getPersonID());
    int pcr_result = -1;
    double sensitivity = 0.0;
    if(h->isInfected() && h->infection != NULL){
	// Is this a primary or secondary infection? 
	// Vaccinees will have a secondary-like viral curve
	double b1 = 0.0;
	double b2 = 0.0;
	if(h->infection->isPrimary() == false || h->isVaccinated() == true){
	//if(h->infection->isPrimary() == false){
	    // secondary infection
	    b1 = 6.834631;
	    b2 = -1.166282;
	}else{
	    b1 = 13.185066;
	    b2 = -1.665468;
	}
	sensitivity = (1.0 / (1.0 + exp(-1 * (b1 + b2 * (currDay - h->infection->getSymptomOnset())))));
	if(rGen->getEventProbability() < sensitivity){
	    pcr_result = h->infection->getInfectionType() - 1;
	}
    }
    return pcr_result;
}


void Surveillance::printRecords(string file, int currDay){
    if (file.length() == 0) {
	exit(1);
    }
    std::ofstream outSurveillance;
    outSurveillance.open(file);
    if (!outSurveillance.good()) {
	exit(1);
    }
    
    vector<string> headers; string outstring;

    headers.push_back("ID,Age,Agegroup,Arm,Serostatus,Enrollment_day,Last_day,Dropout");
    for(unsigned i = 0; i < 4; i++){
	string nn = std::to_string(i + 1);
	if(printExposure){
	    headers.push_back("previous_exposure_"+ nn + ",onset_" + nn + ",symptoms_" + nn + ",severity_" + nn + ",PCRday_" + nn + ",PCR_" + nn + ",TYPE_" + nn + ",TTEL_" + nn + ",TTER_" + nn + ",numexp_" + nn + ",dateexp_" + nn + 
			      ",first_real_inf_" + nn + ",first_real_symp_" + nn);
	}else{
	    headers.push_back("previous_exposure_"+ nn + ",onset_" + nn + ",symptoms_" + nn + ",severity_" + nn + ",PCRday_" + nn + ",PCR_" + nn + ",TYPE_" + nn + ",TTEL_" + nn + ",TTER_" + nn + ",numexp_" + nn +
			      ",first_real_inf_" + nn + ",first_real_symp_" + nn);
	}
    }
    headers.push_back("firstTTL, firstTTR, firstPCR, firstNumExp,firstRealInf,firstRealSymp,bites,attractiveness");
    Surveillance::join(headers,',',outstring);

    outSurveillance << outstring;

    map<string, hRecord>::iterator it;
    for(it = recordsDatabase.begin(); it != recordsDatabase.end(); ++it){
	string serostatus = ((*it).second.seroStatusAtVaccination == true) ? "POSITIVE" : "NEGATIVE";
	string drop_str = ((*it).second.dropout == true) ? "TRUE" : "FALSE";
	string age_g = (*it).second.ageDaysAtVaccination / 365 < 18 ? "Children" : "Adults";
	int lastDay = (*it).second.dropoutDay > -1 ? (*it).second.dropoutDay : currDay;
	outSurveillance << (*it).second.houseID.c_str() << "_" << (*it).second.houseMemNum << "," << (*it).second.ageDaysAtVaccination / 365 << "," << age_g <<","<< (*it).second.trialArm.c_str() << "," << serostatus.c_str() << ",";
	outSurveillance << (*it).second.enrollmentDay << "," << lastDay << "," << drop_str << ",";
	for(int i = 0; i < 4; i++){
	    string sympt = "NA";
	    string hosp = "NA";
	    string real_inf = "NA";
	    string real_symp = "NA";
	    string onset = ((*it).second.onset[i] < 0) ? "NA" : std::to_string((*it).second.onset[i]);
	    if((*it).second.real_infected_day[i] >= 0){
		real_inf = std::to_string((*it).second.real_infected_day[i]);
		if((*it).second.real_symptomatic_day[i] >= 0){
		    real_symp = std::to_string((*it).second.real_symptomatic_day[i]);
		}
	    }	    
	    if((*it).second.symptoms[i] == 1){
		sympt =  "symptomatic";
		hosp = ((*it).second.hosp[i] == 1) ? "hospitalized" : "mild";
	    }
	    string ttr = (*it).second.TTR[i] >= 0 ? std::to_string((*it).second.TTR[i]) : "NA";
	    string pcr_day = (*it).second.pcrDay[i] >= 0 ? std::to_string((*it).second.pcrDay[i]) : "NA";
	    outSurveillance << (*it).second.previousExposure[i] << "," << onset << "," << sympt << "," << hosp << "," << pcr_day << ","<< (*it).second.pcr[i].c_str() << "," << (*it).second.primary[i] << "," << (*it).second.TTL[i] << "," << ttr;
	    if(printExposure){
		outSurveillance << "," << (*it).second.numExp[i] << "," << (*it).second.dateExposure[i] << "," << real_inf << "," << real_symp << ",";
	    }else{
		outSurveillance << "," << (*it).second.numExp[i] << "," << real_inf << "," << real_symp << ",";
	    }
	}
	string ttr = (*it).second.firstTTR >= 0 ? std::to_string((*it).second.firstTTR) : "NA";
	string first_real_inf = (*it).second.first_real_infected >= 0 ? std::to_string((*it).second.first_real_infected) : "NA";
	string first_real_symp = (*it).second.first_real_symptomatic >= 0 ? std::to_string((*it).second.first_real_symptomatic) : "NA";
	outSurveillance << (*it).second.firstTTL << "," << ttr << "," << (*it).second.firstPCR << "," << (*it).second.firstExp << "," << first_real_inf << "," << first_real_symp << "," << (*it).second.bites << "," << (*it).second.attractiveness <<"\n"; 
    }
    outSurveillance.close();

}

void Surveillance::addParameter(string line){
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

int Surveillance::readParameter(string param_name, int vtemp){
    map<string, string>::iterator it;
    int values_ = vtemp;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	values_ = this->parseInteger(it->second);
    }
    return values_;
}
double Surveillance::readParameter(string param_name, double vtemp){
    map<string, string>::iterator it;
    double values_ = vtemp;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	values_ = this->parseDouble(it->second);
    }
    return values_;
}

bool Surveillance::readParameter(string param_name, bool vtemp){
    map<string, string>::iterator it;
    bool values_ = vtemp;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	values_ = this-parseBoolean(it->second);
    }
    return values_;
}

bool Surveillance::parseBoolean(string line){
    bool flag_temp = (std::stoi(line.c_str(), NULL, 10) == 0 ? false : true);
    return flag_temp;
}


int Surveillance::parseInteger(string line){
    return strtol(line.c_str(), NULL, 10);
}

double Surveillance::parseDouble(string line){
    return strtod(line.c_str(), NULL);
}

void Surveillance::join(const vector<string>& v, char c, string& s) {
    s.clear();
    for (vector<string>::const_iterator p = v.begin(); p != v.end(); ++p) {
	s += *p;
	if (p != v.end() - 1){
	    s += c;
	}else{
	    s += "\n";
	}
    }
}

//Surveillance::~Surveillance() {
//}
