#ifndef RANDOMNUMGENERATOR_H
#define	RANDOMNUMGENERATOR_H

#include <random>
#include <map>
#include <sstream>
#include <iostream>
#include <algorithm>

using std::string;
using std::map;

class RandomNumGenerator {
public:
    template <class RandIt>
    void shuffle(RandIt & obj) {
        std::shuffle(obj.begin(), obj.end(), gen);
    }
    unsigned getMozEmerge(double);
    unsigned getMozEmerge(double, double);
    unsigned getHumanTrajectory();
    unsigned getHumanImmunity();
    unsigned getVaxHumanImmunity(unsigned);
    unsigned getRandomNum(unsigned);
    unsigned getMozNextLoc(unsigned);
    unsigned getImmatureMozDeath(double);
    unsigned getMozDevelopment(double);

    bool getHumanSeropositivity(double, double);

    double getWaningTime(double);
    double getMozLifeSpan();
    double getMozLifeSpan(double);
    double getMozDeathRate(double);
    double getRandomNormal();
    double getEventProbability();
    double getMozLatencyDays(double);
    double getMozLatencyRate(double);
    double getMozRestDays();
    double getMozRestDays(double);
    double getAttractiveness();
    double getHypoExpDensity(double, double, double, double);

    int getSelfReportDay(double);

    void setSeed(unsigned);
    string toString() const;
    RandomNumGenerator(unsigned, unsigned, double, double, double, map<unsigned,double>,double);
    RandomNumGenerator();
    RandomNumGenerator(const RandomNumGenerator& orig);
    //virtual ~RandomNumGenerator();

    // print complete rng state
    // absurd amounts of output
    void showAllState() {
        std::cout << gen << std::endl;
    }
    // print RNG state summary (first nchar of full state)
    void showState(unsigned nchar, const string prefix="## ") {
        char buff[nchar];
        std::stringstream the_state;
        the_state << gen;
        the_state.read(buff,nchar);
        std::cout << prefix << buff << std::endl;
    };

private:
    unsigned seed;
    std::mt19937 gen;
    unsigned huImmunity;
    double emergeFactor;
    double mozLife;
    double mozRest;
    double attractShape;
    map<unsigned,double> halflife;
};

#endif	/* RANDOMNUMGENERATOR_H */
