#ifndef INFECTION_H
#define	INFECTION_H

#include <cmath>
#include <string>
#include <iostream>

class Infection {
public:
    int getStartDay() const;
    int getEndDay() const;
    double getSymptomOnset() const;
    double getInfectiousness() const;
    unsigned getInfectionType() const;
    bool isPrimary() const;
    void setInfectiousnessHuman(int);
    void setInfectiousnessMosquito(double);
    void setInfectiousnessMosquito(double, int);
    std::string toString() const;
    Infection(int, unsigned, double, unsigned, bool, bool, double);
    Infection();
    Infection(const Infection& orig);
    //virtual ~Infection();
private:
    int startDay;
    int endDay;
    double infectiousness;
    unsigned infType;
    bool primary;
    bool symptomatic;
    double IIP;
};

#endif	/* INFECTION_H */
