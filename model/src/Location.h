#ifndef LOCATION_H
#define	LOCATION_H

#include <string>
#include <vector>
#include <deque>
#include <set>
#include <forward_list>
#include <memory>
#include <unordered_map>
#include <fstream>
#include "RandomNumGenerator.h"
#include "Human.h"

using std::string;
using std::vector;

class Location {
private:
    string locID;
    double xCor;
    double yCor;
    string MoHID;

    double initialAdults;
    double pupae;
    double larvae;
    double eggs;
    double larvalCapacity;

    std::deque<unsigned> recentAdults;
    std::deque<double> recentDensityDependenceTerms;


    double emergenceRate;
    bool infectedVisitor;
    string locType;
    bool bitesCounterEnabled;
    double insecticideDecayRate;
    double insecticideEfficacy;
    unsigned insecticideStartDay;
  unsigned insecticideEndDay;
  bool insecticideSprayed;
  unsigned insecticideResidualityStartDay;

    std::unique_ptr<vector<string>> closeLocs;
    std::unique_ptr<vector<string>> radiusLocs;
    std::set<sp_human_t,Human::sortid> humans;
    std::unordered_map<string,int>visitorBites;
public:
    string getRandomCloseLoc(RandomNumGenerator&);
    size_t getNumberRadiusLocations(){return radiusLocs->size();}
    void addHuman(sp_human_t);
    void increaseBites(string);
    void enableBitesCounter(){bitesCounterEnabled = true;}
    void disableBitesCounter(){bitesCounterEnabled = false;}
  void sprayAdultInsecticide(unsigned,double,double, unsigned, unsigned);

    double getIncreasedMortalityInsecticide(unsigned, double);
    double calculateGiniIndex();
    void removeHuman(sp_human_t);
    std::set<sp_human_t,Human::sortid> & getHumans(){return humans;}
    void addCloseLoc(string);
    void addRadiusLoc(string);
    void printHumans();
    double getDistanceFromLoc(Location &) const;
    double getLocX() const;
    double getLocY() const;
    double getEmergenceRate() const;
    double getEggs() const;
    double getLarvae() const;
    double getPupae() const;
    double getInitialAdults() const;
    std::deque<unsigned> getRecentAdults() const;
    std::deque<double> getRecentDensityDependenceTerms() const;
    void updatePupae(double);
    void updateLarvae(double);
    void updateEggs(double);
    void updateRecentAdults(unsigned);
    void updateRecentDensityDependenceTerms(double);
    double getLarvalCapacity() const;
    bool getInfectedVisitor(){return infectedVisitor;}
    string getLocID() const;
    string getLocType() const;
    string getMoHID() const;
    vector<string> getRadiusLocations();

  Location(string, string, string, double, double, double);
    Location(string, string, double, double, std::deque<unsigned>, double);
    Location(string, string, string, double, double, double, double, double, double);
    Location(const Location& orig);
    //virtual ~Location();
    void updateInfectedVisitor();

private:

};

#endif	/* LOCATION_H */
