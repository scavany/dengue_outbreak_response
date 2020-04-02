#ifndef SURVEILLANCE_H
#define	SURVEILLANCE_H

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
#include "Human.h"

using std::vector;
using std::string;
using std::map;

class Surveillance {
    struct hRecord{
	int ageDaysAtVaccination;
	int TTL[4];
	int TTR[4];
	int numExp[4];
	int firstTTL;
	int firstTTR;
	int firstExp;
	int first_real_infected;
	int first_real_symptomatic;
	int symptoms[4];
	int hosp[4];
	int pcrDay[4];
	int previousExposure[4];
	int real_infected_day[4];
	int real_symptomatic_day[4];
	int dropoutDay;
	int enrollmentDay;
	int houseMemNum;
	int bites;
	double attractiveness;
	double onset[4];
	bool seroStatusAtVaccination;
	bool dropout;
	string trialArm;
	vector<string> dateExposure;      
	vector<string> primary;
	vector<string> pcr;
	string firstPCR;
	string houseID;
    };
 public:
    Surveillance();
    Surveillance(const Surveillance& orig);
    //virtual ~Surveillance();

    int update_human_surveillance(Human *, int, RandomNumGenerator *);
    void track_infected(Human *, int);
    void initialize_human_surveillance(Human *, int);
    void finalize_human_surveillance(Human *, int, bool);
    void setup(string);
    void printRecords(string, int);

 private:
    int contactFrequency;
    int firstContactDelay;
    double selfReportProb;
    double reportTodayProb;
    bool printExposure;
    
    map<string, hRecord> recordsDatabase;
    map<string, string> parameters;

    void contactPerson(Human *, int, RandomNumGenerator *);
    void addParameter(string);
    void join(const vector<string>& , char , string& );
    int PCR_test(Human *, int, RandomNumGenerator *);
    int parseInteger(string);
    int readParameter(string, int );
    bool parseBoolean(string);
    double parseDouble(string);
    double readParameter(string, double );
    bool readParameter(string , bool );
};

#endif	/* SURVEILLANCE_H */
