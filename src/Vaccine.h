#ifndef VACCINE_H
#define	VACCINE_H

#include <string>
#include <iostream>
#include <utility>
#include <memory>
#include <cmath>
#include <map>
#include <vector>
#include <cstdlib>

#include "RandomNumGenerator.h"

using std::map;
using std::vector;
using std::string;

class Vaccine {

 private:
    int vaccineID;
    int doses;
    double waning;
    double protection;
    double total_VE;
    double propInf;
    double normdev_pos;
    double normdev_neg;
    double seroposVE;
    double seronegVE;
    double seroposWaning;
    double seronegWaning;
    double RRInf_seropos[4];
    double RRInf_seroneg[4];
    double RRDis_seropos[4];
    double RRDis_seroneg[4];
    double RRHosp_seropos[4];
    double RRHosp_seroneg[4];

    string mode;
    string name;

    map<unsigned,double> VE_pos;
    map<unsigned,double> VE_neg;
    vector<int> relative_schedule;

public:
    Vaccine();
    //virtual ~Vaccine();

    void init();
    void setID(unsigned id){vaccineID = id;}
    void setMode(string m){mode = m;}
    void setName(string n){name = n;}
    void setWaning(double w){waning = w;}
    void setWaning(bool seroposIn, double w);
    void setVaccineEfficacy(bool seroposIn, double w);
    void setRRInf(bool seroposIn, double rr, unsigned sero_);
    void setRRDis(bool seroposIn, double rr, unsigned sero_);
    void setRRHosp(bool seroposIn, double rr, unsigned sero_);
    void setRRInf(bool seroposIn, double rr);
    void setRRDis(bool seroposIn, double rr);
    void setRRHosp(bool seroposIn, double rr);
    void setProtection(double p){protection = p;}
    void setTotalVE(double ve){total_VE = ve;}
    void addVE_pos(double v, unsigned a){VE_pos.at(a) = v;}
    void addVE_neg(double v, unsigned a){VE_neg.at(a) = v;}
    void setPropInf(double p){this->propInf = p;}
    void setNormdevPos(double n){this->normdev_pos = n;}
    void setNormdevNeg(double n){this->normdev_neg = n;}
    void setRelativeSchedule(vector<int>);
    void printVaccine();

    int getVaccineID(){return vaccineID;}
    string getMode(){return mode;}

    int getDoses(){return doses;}
    int getNextDoseTime(int, int);

    double getPropInf(){return propInf;}
    double getVaccineProtection(){return protection;}
    double getWaning(){return waning;}
    double getRR(double, double);
    double getVE(){return total_VE;}
    double getVaccineEfficacy(bool seroposIn);
    double getRRInf(bool seroposIn, unsigned sero_);
    double getRRDis(bool seroposIn, unsigned sero_);
    double getRRHosp(bool seroposIn, unsigned sero_);    
    double getWaning(bool seroposIn);
};

#endif	/* VACCINE_H */
