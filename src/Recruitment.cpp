#include <fstream>
#include <string>
#include <iostream>
#include <utility>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <sstream>
#include "Surveillance.h"
#include "Recruitment.h"
#include "Human.h"

Recruitment::Recruitment(){
    vaccineSampleSize = 0;
    placeboSampleSize = 0;
    recruitmentTimeFrame = 0;
    dailyVaccineRecruitmentRate = 0;
    dailyPlaceboRecruitmentRate = 0;
    recruitmentStartDay = 0;
    trialDurationDays = 0;
    trialMaximumDays = 0;
    trialMinimumCases = 0;
    pcr_cases = 0;
    dropoutRate = 0.0;

    recruitmentStrategy = "none";
    recruitmentZone = "";
    ageGroups.clear();

    vaccineProfile = 0;
    placeboProfile = 0;
    trialLastDay = -1;
}

int Recruitment::update(unsigned currDay){
    printf("Updating trial day %d\n", currDay);
    if(currDay >= recruitmentStartDay){
	if(currDay < recruitmentStartDay + recruitmentTimeFrame){
	    this->enrollTodayParticipants(currDay);
	}
	if(currDay > recruitmentStartDay && currDay <= ( recruitmentStartDay + recruitmentTimeFrame + trialMaximumDays)){
	    bool trialEnd = this->updateParticipants(currDay);
	    if(trialEnd == true && trialLastDay == -1){
		trialLastDay = currDay;
	    }
	}else if(currDay == recruitmentStartDay + recruitmentTimeFrame){
	    printf("Recruitment finished\n");
	    for(unsigned i = 0; i < ageGroups.size(); i ++){
		if(ageGroups[i].vaccine.size() != vaccineSampleSize || ageGroups[i].placebo.size() != placeboSampleSize ){
		    printf("Recruitment finished: take a look at the size  age group %d of vaccine %lu and placebo %lu\n",
			   i, ageGroups[i].vaccine.size(),ageGroups[i].placebo.size());
		    exit(1);
		}else{
		    printf("age group: %d - %d, size vaccine: %lu, size placebo: %lu. Eligibles: %lu\n",
			   ageGroups[i].min, ageGroups[i].max, ageGroups[i].vaccine.size(), ageGroups[i].placebo.size(), ageGroups[i].eligible.size());
		    ageGroups[i].eligible.clear();
		}
	    }
	}
    }
    printf("Updating finished day %d\n", currDay);
    return(trialLastDay);
}

void Recruitment::removeParticipant(Human * h, int currDay, bool drop_in){
    trialSurveillance.finalize_human_surveillance(h, currDay,drop_in);
    h->unenrollTrial();
}

void Recruitment::finalizeTrial(unsigned currDay){
    printf("Finalizing trial. Day: %d, PCR cases: %d\n", currDay, pcr_cases);
    trialSurveillance.printRecords(outSurveillance, currDay);
}

bool Recruitment::updateParticipants(unsigned currDay){
    // update doses if needed, remove dead people  and test for denv: self-reported and calls
    long int active_participants = 0;
    for(unsigned i = 0;i < ageGroups.size(); i++){
	updateArm(vaccineProfile, ageGroups[i].vaccine, currDay);
	updateArm(placeboProfile, ageGroups[i].placebo, currDay);
	active_participants += ageGroups[i].vaccine.size();
	active_participants += ageGroups[i].placebo.size();
    }
    if(active_participants == 0){
	printf("trial finished in day %u\n", currDay);
	return(true);
    }else{
	return(false);
    }
}

void Recruitment::updateArm(unsigned vaxID, recruit_t & arm, int currDay){
    // boost vaccine, decide if dropout, remove death people from the list
    //    printf("Update Arm\n");
    if(vaccinesPtr.at(vaxID).getDoses() > 1){
        // added redundant continue statements
        // to ensure that "it" only changes (delete or increment) once per for()
	for(auto it = arm.begin(); it != arm.end(); ){
            // readability, reset after each continue
            Human * phum = (*it);
	    if(phum == nullptr){
		printf("We found a null pointer\n");
		it = arm.erase(it);
                continue;
	    } else {
		if(phum->isEnrolledInTrial() == true){
		    if(phum->isInfected() == true){
			// count symptomatic infections even if they don't get caught by pcr
			trialSurveillance.track_infected(phum,currDay);
		    }
		    if((currDay < phum->getTrialEnrollmentDay() + (int)trialMaximumDays ) && (currDay < phum->getTrialEnrollmentDay() + trialDurationDays ||  pcr_cases < (int)trialMinimumCases)){
			if(currDay > (int)recruitmentStartDay + (int)recruitmentTimeFrame && rGen->getEventProbability() < dropoutRate){
			    removeParticipant(phum,currDay,true);
			    it = arm.erase(it);
                            continue;
			}else{
			    int pcr = trialSurveillance.update_human_surveillance(phum, currDay, rGen);
			    if(pcr >= 0){
                                //printf("PERSON %s removed from trial for real\n", phum->getPersonID().c_str());
				pcr_cases++;
				removeParticipant(phum,currDay,false);
				it = arm.erase(it);
                                continue;
			    }else{
				if(phum->isFullyVaccinated() == false && phum->getNextDoseDay() == currDay){
				    //				    printf("human boosting %s\n", phum->getPersonID().c_str());
				    phum->boostVaccine(currDay, rGen);
				}
				it++;
                                continue;
			    }
			}
		    }else{
			//			printf("Participant removed because time is over\n");
			removeParticipant(phum,currDay,false);
			it = arm.erase(it);
                        continue;
		    }
		}else{
		    //		    printf("Participant is not really enrolled in trial\n");
		    it = arm.erase(it);
                    continue;
		}
            }
	}
    }
    //    printf("Arm updated\n");
}
void Recruitment::printEligibleGroups(){
    for(unsigned i = 0; i < ageGroups.size(); i ++){
	printf(" Eligibles age-group: %d - %d, (eligible size %zu).\n", ageGroups[i].min, ageGroups[i].max, ageGroups[i].eligible.size() );
    }
}

void Recruitment::enrollTodayParticipants(int currDay){
    if(dailyVaccineRecruitmentRate <= 0 || dailyPlaceboRecruitmentRate <= 0){
	printf("Daily recruitment rate <= 0\n");
	exit(1);
    }
    for(unsigned i = 0; i < ageGroups.size(); i ++){
	//
	//Vaccine enrollment
	enrollArmParticipants(ageGroups[i].vaccine, ageGroups[i].eligible, "vaccine", currDay,vaccineSampleSize,dailyVaccineRecruitmentRate,ageGroups[i].min, ageGroups[i].max,vaccineProfile);
	enrollArmParticipants(ageGroups[i].placebo, ageGroups[i].eligible, "placebo", currDay,placeboSampleSize,dailyPlaceboRecruitmentRate,ageGroups[i].min, ageGroups[i].max,placeboProfile);
	//	printf("%zu vaccine and %zu placebo participants successfully enrolled at day %d\n", ageGroups[i].vaccine.size(), ageGroups[i].placebo.size(), currDay);
    }
}

void Recruitment::enrollArmParticipants(
					recruit_t & arm,
					eligible_t & eligible,
					string arm_str, int currDay,
					unsigned sample_size,
					int rec_rate,
					int min_,
					int max_,
					unsigned vProfile)
{
    //    printf("Enroll participants today: %d, min: %d, max: %d, arm: %s size: %zu, (eligible %zu)\n", currDay, min_, max_, arm_str.c_str(), arm.size(), eligible.size());
    int nrecruit = 0;
    // process from end, possibly deleting as we go
    auto it = eligible.end();
    while(it >= eligible.begin() && nrecruit < rec_rate && arm.size() < sample_size){
        // back up from end
        it--;
        // convenience / clarity var
        // reset after each continue
        Human * phum = (*it);
        //
	if(phum == nullptr){
	    printf("Found a dead participant in enrollment day %d\n", currDay);
        // remove and restart loop
	    it = eligible.erase(it);
	    continue;
	}
	double temp_age = (double) phum->getAgeDays(currDay) / 365.0;
	if(temp_age >= max_){
        // remove and restart loop
	    it = eligible.erase(it);
	    continue;
	} else if(temp_age < min_){
        // skip, reshuffle at end? sounds good
	    continue;
            /*
	    Human * temp_h = eligible_vector->back();
	    eligible_vector->pop_back();
	    vector<Human *>::iterator it = eligible_vector->begin();
	    long unsigned pos = rGen->getRandomNum(eligible_vector->size());
	    eligible_vector->insert(it+pos,temp_h);
            */
	}else{
	    if(phum->isEnrolledInTrial() == false){
		phum->enrollInTrial(currDay, arm_str);
		arm.insert(phum);
		trialSurveillance.initialize_human_surveillance(phum, currDay);
		phum->vaccinateWithProfile(currDay, rGen, vaccinesPtr.at(vProfile));
		nrecruit++;
		it = eligible.erase(it);
	    }
        continue; // added for clarity
	}
    }
    // reshuffle at the end of each enrollment
    shuffleEligibleParticipants();
}

void Recruitment::setupRecruitment(string file, map<unsigned,Vaccine> vaccines_, string outputPath, string simName_, RandomNumGenerator * _rGen){
    // store reference to simulator-wide rng
    rGen = _rGen;
    // Initialize vaccines profiles
    vaccinesPtr = vaccines_;
    if(vaccinesPtr.size() == 0){
	printf("Recruitment::setupRecruitment there are no vaccines to copy\n");
	exit(1);
    }
    //Setup the surveillance system for the trial
    trialSurveillance.setup(file);

    if (file.length() == 0) {
		exit(1);
    }
    outSurveillance = outputPath + "/" + simName_ + "_trial.csv";

    string line;
    std::ifstream infile(file);

    if(!infile.good()){
		exit(1);
    }

    printf("Reading recruitment setup file %s\n", file.c_str());

    // Read the trial recruitment parameters
    while(getline(infile,line,'\n')){
	this->addParameter(line);
    }

    double tmp;
    this->readParameter("trial_age_groups",&ageGroups);
    this->readParameter("trial_recruitment_strategy", &recruitmentStrategy);
    this->readParameter("trial_recruitment_zone", &recruitmentZone);
    this->readParameter("trial_vaccine_profile", &vaccineProfile);
    this->readParameter("trial_placebo_profile", &placeboProfile);
    this->readParameter("trial_recruitment_start_day", &recruitmentStartDay);
    this->readParameter("trial_recruitment_timeframe", &recruitmentTimeFrame);
    this->readParameter("trial_vaccine_sample_size", &vaccineSampleSize);
    this->readParameter("trial_placebo_sample_size", &placeboSampleSize);
    this->readParameter("trial_length_days", &trialDurationDays);
    this->readParameter("trial_maximum_days", &trialMaximumDays);
    this->readParameter("trial_minimum_cases", &trialMinimumCases);
    this->readParameter("trial_avg_enrollment_days",&tmp);
    dropoutRate = (double) 1.0 / tmp;

    printf("TRIAL RECRUITMENT SETTINGS----------------\n");
    if(recruitmentStrategy == "zones"){
	printf("starting development of zones... be careful, this is experimental\n");
	printf("Recruitment zone selected %s\n", recruitmentZone.c_str());
    }
    printf("Recruitment strategy: |%s| duration: %d dropoutRate: %.6f recruitment_time_frame: %d Max. Days %d, min cases: %d\n",
	   recruitmentStrategy.c_str(), trialDurationDays, dropoutRate,recruitmentTimeFrame, trialMaximumDays, trialMinimumCases);
    printf("Vaccine Profile ID: %d, placebo profile ID: %d\n",vaccineProfile, placeboProfile);
    if(recruitmentStrategy == "none"){
	printf("Please specify a recruitment Strategy\n");
	exit(1);
    }
    for(unsigned i = 0; i < ageGroups.size(); i ++){
	printf("Age group: %d - %d\n",ageGroups[i].min, ageGroups[i].max - 1);
    }

    infile.close();
}

void Recruitment::addParameter(string line){
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
	    // trim trailing and leading spaces and weird stuff from param_value
	    pos_equal = param_value.find_last_not_of(" \t");
	    if(pos_equal != string::npos){
		param_value = param_value.substr(0,pos_equal+1);
	    }
	    // Add the parameter name and value to the map
	    parameters.insert(make_pair(param_name,param_value));
	}
    }
}

void Recruitment::readParameter(string param_name, vector<Recruitment::groupStruct> * param_var){
   map<string, string>::iterator it;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	this->parseAges(it->second,param_var);
    }
}

void Recruitment::readParameter(string param_name, unsigned * param_var){
   map<string, string>::iterator it;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	*param_var = parseInteger(it->second);
    }
}

void Recruitment::readParameter(string param_name, int * param_var){
   map<string, string>::iterator it;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	*param_var = parseInteger(it->second);
    }
}

void Recruitment::readParameter(string param_name, double * param_var){
   map<string, string>::iterator it;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	*param_var = parseDouble(it->second);
    }
}

void Recruitment::readParameter(string param_name, string * param_var){
   map<string, string>::iterator it;
    it = parameters.find(param_name);
    if(it != parameters.end()){
	*param_var = parseString(it->second);
    }
}

void Recruitment::addPossibleParticipant(Human * h, int currDay){
    // First verify that the age is in some group
    int group_ = getPossibleAgeGroup(h->getAgeDays(currDay),ageGroups,recruitmentTimeFrame);
    if(group_ < 0){
	return;
    }
    // If it's random recruitment just push it into the list without any other requirement
    if(recruitmentStrategy == "random"){
	ageGroups[group_].eligible.push_back(h);
    }else if(recruitmentStrategy == "zones"){
	if(h->getZoneID() == recruitmentZone){
	    ageGroups[group_].eligible.push_back(h);
	}
    }else{
	printf("Please modify recruitment strategy %s is not a valid one: only random and zones are supported at this momment\n", recruitmentStrategy.c_str());
	exit(1);
    }
}

void Recruitment::shuffleEligibleParticipants(){
    // counter
    int ii{};
    // process eligibles for each group
    for(auto & grp : ageGroups ){
        ii++;
	if(grp.eligible.size() > 0){
	    // First make sure eligible participants have not died
	    for(auto it = grp.eligible.begin(); it !=  grp.eligible.end();){
		if( (*it) == nullptr){
		    printf("NULL pointer in shuffle eligible participants\n");
                    // erase and advance (vector safe)
		    it = grp.eligible.erase(it);
                    continue; // redund
		}else{
                    // skip
		    it++;
                    continue; // redund
		}
	    }
	    if(grp.eligible.size() > 0){
		rGen->shuffle(grp.eligible);
	    }
	}else{
	    printf("Eligibles are 0, there are no more to recruit. vaccine arm size: %lu, placebo arm size: %lu\n", grp.vaccine.size(), grp.placebo.size());

	}
    }
    dailyVaccineRecruitmentRate = ceil( (double) vaccineSampleSize / (double) recruitmentTimeFrame);
    dailyPlaceboRecruitmentRate = ceil( (double) placeboSampleSize / (double) recruitmentTimeFrame);
    //    printf("Shuffle participants finished\n");
}

int Recruitment::getAgeGroup(int age_, vector<groupStruct> groups_temp){
    vector<groupStruct>::iterator itAge = groups_temp.begin();
    int count = 0;
    for(; itAge != groups_temp.end(); itAge++){
	if((double )age_ / 365.0 >= (*itAge).min && (double) age_ / 365.0 < (*itAge).max){
	    return count;
	}
	count++;
    }
    return -1;
}

int Recruitment::getPossibleAgeGroup(int age_, vector<groupStruct> groups_temp, int time_temp){
    if(time_temp < 0){
	time_temp = 0;
    }
    vector<groupStruct>::iterator itAge = groups_temp.begin();
    int count = 0;
    for(; itAge != groups_temp.end(); itAge++){
	if((double ) (age_ + time_temp) / 365.0 >= (*itAge).min && (double) age_ / 365.0 < (*itAge).max){
	    return count;
	}
	count++;
    }
    return -1;
}

long int Recruitment::getEligibleParticipantsSize(){
    long int s = 0;
    for(unsigned i = 0; i < ageGroups.size(); i++ ){
	s += ageGroups[i].eligible.size();
    }
    return s;
}

int Recruitment::parseInteger(string line){
    return strtol(line.c_str(), NULL, 10);
}

double Recruitment::parseDouble(string line){
    return strtod(line.c_str(), NULL);
}

string Recruitment::parseString(string line){
    size_t first_ = line.find_first_not_of(' ');
    size_t last_ = line.find_last_not_of(' ');
    return line.substr(first_,(last_ - first_ + 1));
}

void Recruitment::parseAges(string line, vector<groupStruct> * ages_temp){
    std::stringstream linetemp;
    string line2;
    linetemp.clear();
    linetemp << line;
    ages_temp->clear();
    while(getline(linetemp,line2,';')){
		std::stringstream lTemp; lTemp << line2;
		string line3;
		groupStruct rangeTemp;
		rangeTemp.vaccine.clear();
		rangeTemp.placebo.clear();
		rangeTemp.eligible.clear();
		getline(lTemp,line3,',');
		rangeTemp.min = strtol(line3.c_str(), NULL, 10);
		getline(lTemp,line3,',');
		rangeTemp.max = strtol(line3.c_str(), NULL, 10);
		if(rangeTemp.min + rangeTemp.max > 0){
		    ages_temp->push_back(rangeTemp);
		}
    }

    if(ages_temp->empty()){
		exit(1);
    }
}
void Recruitment::parseVector(string line, vector<double> * vector_temp){
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
		exit(1);
    }
}

//Recruitment::~Recruitment() {
//}
