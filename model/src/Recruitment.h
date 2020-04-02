#ifndef RECRUITMENT_H
#define	RECRUITMENT_H
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <utility>
#include <algorithm>
#include <memory>
#include <cmath>
#include "RandomNumGenerator.h"
#include "Vaccine.h"
#include "Surveillance.h"
#include "Human.h"

using std::string;
using std::vector;
using std::map;

// allow for shuffle without random access
using eligible_t = vector<Human *>;
// filled from above
using recruit_t = std::set<Human *>;

class Recruitment {
    struct groupStruct{
        int min;
        int max;
        eligible_t eligible;
        recruit_t placebo;
        recruit_t vaccine;
    };

 private:
    RandomNumGenerator * rGen;
    Surveillance trialSurveillance;

    map<string, string> parameters;

    unsigned vaccineSampleSize;
    unsigned placeboSampleSize;
    unsigned recruitmentTimeFrame;
    unsigned recruitmentStartDay;
    int dailyVaccineRecruitmentRate;
    int dailyPlaceboRecruitmentRate;
    int trialDurationDays;
    int trialLastDay;
    unsigned trialMaximumDays;
    unsigned trialMinimumCases;
    int pcr_cases;

    string recruitmentStrategy;
    string outSurveillance;
    string recruitmentZone;
    vector<groupStruct> ageGroups;

    unsigned vaccineProfile;
    unsigned placeboProfile;

    double dropoutRate;

    int getAgeGroup(int, vector<groupStruct>);
    int getPossibleAgeGroup(int, vector<groupStruct>, int);
    int parseInteger(string);
    double parseDouble(string);

    void addParameter(string line);
    void readParameter(string param_name, int * param_var);
    void readParameter(string param_name, unsigned * param_var);
    void readParameter(string param_name, double * param_var);
    void readParameter(string param_name, string * param_var);
    void readParameter(string param_name, vector<Recruitment::groupStruct> * param_var);

    void parseAges(string line, vector<groupStruct> *);
    void parseVector(string line, vector<double> *);
    void enrollTodayParticipants(int);
    void enrollArmParticipants(recruit_t &, eligible_t &, string, int, unsigned, int, int, int,unsigned);
    string parseString(string);
    map<unsigned, Vaccine> vaccinesPtr;

 public:

    Recruitment(string);
    Recruitment();
    Recruitment(const Recruitment& orig);
    //virtual ~Recruitment();

    int update(unsigned);
    void updateArm(unsigned, recruit_t &, int);
    bool updateParticipants(unsigned);
    // pseudo-ctor
    void setupRecruitment(string, map<unsigned,Vaccine>, string, string, RandomNumGenerator * _rGen);
    void addPossibleParticipant(Human *, int);
    void shuffleEligibleParticipants();
    void finalizeTrial(unsigned);
    void removeParticipant(Human *, int, bool);
    void printEligibleGroups();
    int getVaccineSampleSize(){return vaccineSampleSize;}
    int getPlaceboSampleSize(){return placeboSampleSize;}
    unsigned getRecruitmentStartDay(){return recruitmentStartDay;}
    long int getEligibleParticipantsSize();
};

#endif	/* RECRUITMENT_H */
