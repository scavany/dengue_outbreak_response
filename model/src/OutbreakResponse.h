#ifndef OUTBREAKRESPONSE_H
#define	OUTBREAKRESPONSE_H

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <unordered_set>
#include <iostream>
#include <utility>
#include <algorithm>
#include <memory>
#include <cmath>
#include "RandomNumGenerator.h"
#include "Human.h"
#include "Location.h"
#include "Report.h"

using std::string;
using std::queue;
using std::map;
using std::unordered_set;
using std::vector;
using std::deque;
using std::min;

using symptomatic_t = vector<Human *>;

class OutbreakResponse {
public:
  OutbreakResponse();
  OutbreakResponse(const OutbreakResponse& orig);

  void setup(string file);
  void setLocations(map <string,std::unique_ptr<Location>> *);
  void update(unsigned currDay, RandomNumGenerator * rGen, Report * r);
  void addSymptomatic(unsigned currDay, Human * h, RandomNumGenerator * rGen, Report * r);

  bool isOutbreakResponseActive(unsigned currDay);
  bool isOutbreakResponseEnabled();
  string outbreakResponseStrategy();

  double getResponseRadius(){return spatialRadius;}
  unsigned getTodaySymptomatics(){return todaySympCounts;}

  void finalizeOutbreakResponse();

  void readZonalFile(string);

  bool updateCheckIncidence(unsigned, unsigned);
  double calculateMean(unsigned);
  double calculateSigma(unsigned, double);

private:

  double surveillanceEffort;
  double responseThreshold;
  double aggressiveness;
  double thoroughness;
  double spatialRadius;
  double residuality;
  double compliance;

  unsigned surveillanceDelay;
  unsigned symptomaticCases;
  unsigned maxHouses;
  unsigned maxHousesPerDay;
  unsigned startDay;
  unsigned endDay;
  unsigned totalHousesSprayed;
  unsigned todaySympCounts = 0;
  unsigned simultaneousZones;
  unsigned numberSprayCycles;
  unsigned totalZones;
  unsigned trainingLength;
  unsigned sprayCycleDelay;
  unsigned cycleCounter;
  unsigned sprayCounter = 0;
  unsigned cycleEndDay;
  unsigned cycleStartDay;
  unsigned returnVisits;
  unsigned weekLength;
  unsigned numberSigmas;
  unsigned numberYearsForMean;
  unsigned currentIncidence = 0;
  unsigned subdivisionsOfYear;
  unsigned lengthOfSubdivisions;
  unsigned spraysPerYear;
  unsigned repeatSprayGap;
  unsigned efficacyLength;
  unsigned residualityLag;
  //firstDay governs the day of the week that day 0 is on; 0 = Monday etc. 1 Jan 2000 was a Saturday
  unsigned firstDay = 5;
  unsigned seasonStartDay = 1;
  string zonalOrderFile;
  string startingMethod;

  string responseStrategy;
  string incidenceFile;
  
  bool activeOutbreakResponse;
  bool repeatYearly;
  bool thresholdPassed = false;
  symptomatic_t todaySymptomatics;
  symptomatic_t delayedSymptomatics;

  vector<symptomatic_t> symptomaticRecords;
  deque<string> housesToSpray;
  queue<string> zonesToSpray;
  queue<string> orderedZones;
  unordered_set<string> housesEnlisted;
  unordered_set<string> zonesToEnlist;
  map < string,std::unique_ptr<Location> > * locations_ptr;

  map<string, string> parameters;

  vector<vector<unsigned>> subdividedIncidence;
  vector<double> meanIncidence;
  vector<double> sigmaIncidence;
  
  void sprayTodayLocations(unsigned currDay, Report * r);

  // Routines to read parameters
  void addParameter(string);
  int parseInteger(string);
  int readParameter(string, int );
  bool parseBoolean(string);
  double parseDouble(string);
  double readParameter(string, double);
  bool readParameter(string , bool);
  void readParameter(string , string *);
  string parseString(string);
  void readIncidenceFile(string);
};

#endif	/* OUTBREAKRESPONSE_H */
