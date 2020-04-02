#include "RandomNumGenerator.h"
#include <iostream>
#include <cmath>
#include <sstream>

double RandomNumGenerator::getWaningTime(double wan) {
    std::exponential_distribution<> d(1./wan);
    return d(gen);
}

double RandomNumGenerator::getAttractiveness(){
    if(attractShape > 0){
	std::gamma_distribution <> d(attractShape,1.0);
	return d(gen);
    }else if (attractShape == 0){
	return 1;
    }else{
	return -1;
    }
}

unsigned RandomNumGenerator::getMozDevelopment(double devMean) {
    std::poisson_distribution<> dis(devMean);
    return dis(gen);
}

double RandomNumGenerator::getHypoExpDensity(double x, double l1, double l2, double l3){
  return l1*l2*l3*((l3-l1)*std::exp(-x*l2) + (l1-l2)*std::exp(-x*l3) + (l2-l3)*std::exp(-x*l1))/((l1-l3)*(l2-l1)*(l3-l2));
}

unsigned RandomNumGenerator::getImmatureMozDeath(double deathMean) {
    std::poisson_distribution<> dis(deathMean);
    return dis(gen);
}

unsigned RandomNumGenerator::getMozEmerge(double mozMean) {
    std::poisson_distribution<> dis(emergeFactor * mozMean);
    return dis(gen);
}

unsigned RandomNumGenerator::getMozEmerge(double mozMean, double seasonalFactor) {
    std::poisson_distribution<> dis( mozMean * seasonalFactor);
    return dis(gen);
}

double RandomNumGenerator::getMozLifeSpan() {
    std::exponential_distribution<> d(1./mozLife);
    return d(gen);
}

double RandomNumGenerator::getMozDeathRate(double deathIn) {
    std::exponential_distribution<> d(deathIn);
    return (1 / d(gen));
}

double RandomNumGenerator::getMozLifeSpan(double deathIn) {
    std::exponential_distribution<> d(deathIn);
    return d(gen);
}

double RandomNumGenerator::getMozLatencyDays(double muIn) {
    //    lognormal_distribution<> d(1.648721, 0.451754);
    std::lognormal_distribution<> d(muIn,0.451754);
    return d(gen);
}

double RandomNumGenerator::getMozLatencyRate(double muIn) {
    //    lognormal_distribution<> d(1.648721, 0.451754);
    std::lognormal_distribution<> d(muIn,0.451754);
    double EIPrate = 1 / d(gen);
    return EIPrate;
}

double RandomNumGenerator::getMozRestDays() {
    std::exponential_distribution<> d(mozRest);
    return d(gen);
}

double RandomNumGenerator::getMozRestDays(double mozRestIn) {
    std::exponential_distribution<> d(mozRestIn);
    return d(gen);
}

unsigned RandomNumGenerator::getMozNextLoc(unsigned num) {
    std::uniform_int_distribution<> dis(0, num-1);
    return dis(gen);
}

unsigned RandomNumGenerator::getHumanTrajectory() {
    std::uniform_int_distribution<> dis(0, 4);
    return dis(gen);
}

unsigned RandomNumGenerator::getHumanImmunity() {
    std::exponential_distribution<> d(1./huImmunity);
    return ceil(d(gen));
}

unsigned RandomNumGenerator::getVaxHumanImmunity(unsigned immdays) {
  std::exponential_distribution<> d(1./immdays);
  return ceil(d(gen));
}


bool RandomNumGenerator::getHumanSeropositivity(double FOI, double age) {
    if(getEventProbability() < 1 - exp(-FOI * age)){
        return true;
    } else {
        return false;
    }
}

unsigned RandomNumGenerator::getRandomNum(unsigned num) {
    std::uniform_int_distribution<> dis(0, num-1);
    return dis(gen);
}

int RandomNumGenerator::getSelfReportDay(double IIP){
    return floor(IIP);
}

double RandomNumGenerator::getRandomNormal(){
    std::normal_distribution<> d(0.0, 1.0);
    return d(gen);
}

double RandomNumGenerator::getEventProbability() {
    std::uniform_real_distribution<> dis(0, 1);
    return dis(gen);
}

void RandomNumGenerator::setSeed(unsigned s) {
    seed = s;
    gen.seed(s);
}

string RandomNumGenerator::toString() const {
    std::stringstream ss;
    ss <<" huImmunity:" << huImmunity;
    ss <<" emergeFactor:" << emergeFactor;
    ss <<" mozLife:" << mozLife;
    ss <<" mozRest:" << mozRest;
    return ss.str();
}

RandomNumGenerator::RandomNumGenerator(
    unsigned s, unsigned huImm, double efactor, double mlife,
    double mbite, map<unsigned,double> hlife, double atShape)
{
    seed = s;
    gen.seed(s);
    huImmunity = huImm;
    emergeFactor = efactor;
    mozLife = mlife;
    mozRest = mbite;
    halflife = hlife;
    attractShape = atShape;
}

RandomNumGenerator::RandomNumGenerator() {
}

RandomNumGenerator::RandomNumGenerator(const RandomNumGenerator& orig) {
}

//RandomNumGenerator::~RandomNumGenerator() {
//}
