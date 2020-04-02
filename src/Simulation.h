#ifndef SIMULATION_H
#define SIMULATION_H

#include <map>
#include <set>
#include <memory>
#include <fstream>
#include <queue>
#include <stdexcept>
#include <array>

#include "defines.h"
#include "Location.h"
#include "Mosquito.h"
#include "RandomNumGenerator.h"
#include "Report.h"
#include "Vaccine.h"
#include "Recruitment.h"
#include "Human.h"
#include "OutbreakResponse.h"

using std::string;
using std::vector;
using std::ifstream;
using std::make_pair;
using traject_t = std::array<vpath_t, N_TRAJECTORY >;



class Simulation {
public:
  RandomNumGenerator rGen;
  RandomNumGenerator rGenInf;
  RandomNumGenerator rGenControl;

  Simulation(string);
  Simulation() = delete;
  Simulation(const Simulation& orig) = delete;
  string readInputs();
  int readParameter(string, int);
  double readParameter(string, double);
  bool readParameter(string, bool);
  string readParameter(string, string);
  void readParameter(string , string *);
  bool checkParameterInf(string);
  void addParameter(string);
  void readTrajectoryFile(string);
  void readBirthsFile(string);
  void readInitialInfectionsFile(string);
  void readSimControlFile(string);
  void readLocationFile(string);
  void readDiseaseRatesFile();
  void readVaccineProfilesFile();
  void readVaccinationGroupsFile();
  void readVaccineSettingsFile();
  void readAegyptiFile(string);
  void readAegyptiImmatureFile(string);
  void readAegyptiDelayedFile(string);
  void readFixedParameters(string);
  void readInitialFOI(string);
  void readDailyFOI(string);
  void readDailyLocalFOI(string);
  void setLocNeighborhood(double);
  void setLocNeighborhoodGrid();
  void setRadiusLocations(double);
  void simEngine();
  void humanDynamics();
  void tests();
  void mosquitoDynamics();
  void generateMosquitoes();
  void implicitEarlyStages();
  unsigned setInitialInfection(double, unsigned);
  void simulate();
  void updatePop();
  void selectEligibleTrialParticipants();
  bool checkAgeToVaccinate(int age_);
  void readBurnPopnFile(string);
  void readDensityFile(string);
  //virtual ~Simulation();
private:
  vector<string>getParamsLine(string);
  int parseInteger(string);
  string parseString(string);
  double parseDouble(string);
  void parseVector(string line, vector<int> *);
  void parseVector(string line, vector<double> *);
  bool parseBoolean(string);

  map<string, string> parameters;
  map <string,std::unique_ptr<Location>> locations;
  unordered_map<int, unordered_map<int,vector<string>>> neighborhoodGrid;
  map<unsigned,double> halflife;
  map<unsigned,double> disRates;
  map<unsigned,double> hospRates;
  map<int,int> ageGroups;
  map<unsigned, Vaccine> vaccines;
  std::multimap<string,std::unique_ptr<Mosquito>> mosquitoes;
  std::multimap<string,sp_human_t> humans;
  std::multimap<int,sp_human_t> future_humans;
  map<string, sp_human_t> total_humans_by_id;

  std::set<std::string> zones;
  std::map<std::string,unsigned> zones_counts;
  unsigned currentDay;
  unsigned numDays;
  string trajectoryFile;
  string birthsFile;
  string configLine;
  string locationFile;
  string vaccineProfilesFile;
  string vaccineSettingsFile;
  string trialSettingsFile;
  string vaccinationStrategy;
  string reportsFile;
  string diseaseRatesFile;
  string vaccinationGroupsFile;
  string aegyptiRatesFile;
  string fixedParameterFile;
  double alpha;
  double beta;

  bool routineVaccination;
  bool catchupFlag;
  bool trialVaccination;
  bool randomTrial;
  bool testBeforeVaccine;
  double routineTestSpecificity;
  double routineTestSensitivity;
  int vaccineID;
  unsigned vaccineDay;
  unsigned vaccineAge;
  double vaccineCoverage;
  Report outputReport;
  Recruitment recruitmentTrial;
  string outputFile;
  string outputPopFile;
  string outputPrevacFile;
  string simName;
  unsigned rSeed;
  unsigned rSeedInf;
  string outputPath;
  unsigned humanInfectionDays;
  unsigned huImm;
  double mlife;
  double mbite;
  double emergeFactor;
  double correctionTermA;
  double correctionTermB;
  double biteProbablity;
  double mozInfectiousness;
  double mozMoveProbability;
  double closeLocDist;
  double attractShape;
  double selfReportProb;
  unsigned year;

  bool earlyStages;
  bool delayedEmergence;

  double normalizedInitialPupae;
  double normalizedInitialLarvae;
  double normalizedInitialEggs;
  double normalizedLarvalCapacity;
  double larvalPower;
  double eggsPerGonCycle;
  unsigned maxn;
  string burnPopnFile;
  string densityFile;
  double totalInitialPopulation;
  double extraMortalityRate;

  double gonotrophicCycleRate;

  OutbreakResponse outbreakIntervention;
  string outbreakFile;
  
  int deathMoz;
  int lifeMoz;
  int newMoz;
  int humanDeaths;
  int visitorsCounter;
  vector<double> meanDailyEIP;
  vector<double> firstBiteRate;
  vector<double> secondBiteRate;
  vector<double> mozDailyDeathRate;
  vector<double> dailyEmergenceFactor;
  vector<double> dailyEggDevelopmentRate;
  vector<double> dailyLarvalDevelopmentRate;
  vector<double> dailyPupalDevelopmentRate;
  vector<double> dailyGonotrophicCycleRate;
  vector<double> dailyNormalizedLarvalCapacity;
  vector<double> eggDailyDeathRate;
  vector<double> larvaeDailyDeathRate;
  vector<double> pupaeDailyDeathRate;
  vector<double> extraDailyDeathRate;
  vector<double> dailyDensitySum;
  vector<double> InitialConditionsFOI;
  vector<map<unsigned, double>> dailyForceOfImportation;
  vector<map<unsigned, double>> dailyForceOfLocalImportation;
  vector< vector<double> > densityMatrix;
  vector<double> normalizedEarlyAdults;
  std::deque<double> densityDependenceTerms;

  std::ofstream out;
  std::ofstream outpop;
  std::ofstream outprevac;

};

#endif  /* SIMULATION_H */
