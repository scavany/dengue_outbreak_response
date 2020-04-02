#include "Infection.h"
#include <sstream>

int Infection::getStartDay() const {
    return startDay;
}

int Infection::getEndDay() const {
    return endDay;
}

double Infection::getSymptomOnset() const{
    return IIP + startDay;
}

bool Infection::isPrimary() const{
    return primary;
}

double Infection::getInfectiousness() const {
    return infectiousness;
}

unsigned Infection::getInfectionType() const {
    return infType;
}

void Infection::setInfectiousnessMosquito(double in){
    infectiousness = in;
}
void Infection::setInfectiousnessMosquito(double in, int sday){
    infectiousness = in;
    startDay = sday;
}



void Infection::setInfectiousnessHuman(int currentDay){
    if(symptomatic){
        if(primary){
            infectiousness = exp(-0.2784750 * pow(double(currentDay - startDay) - (IIP - 0.2880570), 2)) / 0.9107737;
        }
        else{
            infectiousness = exp(-0.2784750 * pow(double(currentDay - startDay) - (IIP - 0.2880570), 2)) / 0.9107737;
        }
    }
    else{
        if(primary){
            infectiousness = exp(-0.2784750 * pow(double(currentDay - startDay) - (IIP - 0.2880570), 2)) / 0.9107737;
        }
        else{
            infectiousness = exp(-0.2784750 * pow(double(currentDay - startDay) - (IIP - 0.2880570), 2)) / 0.9107737;
        }
    }
}

std::string Infection::toString() const {
    std::stringstream ss;
    ss << infType << " " << startDay <<" " << endDay << " " << infectiousness;
    return ss.str();
}

Infection::Infection(int sd, unsigned ed, double infn, unsigned ityp, bool prim, bool symp, double IIPin) {
    startDay = sd;
    endDay = ed;
    infectiousness = infn;
    infType = ityp;
    primary = prim;
    symptomatic = symp;
    IIP = IIPin;
}

Infection::Infection() {
    /*    startDay = -1;
    endDay = -1;
    infectiousness = -1;
    infType = -1;
    primary = false;
    symptomatic = false;
    IIP = -1;*/
}

Infection::Infection(const Infection& orig) {
    startDay = orig.startDay;
    endDay = orig.endDay;
    infectiousness = orig.infectiousness;
    infType = orig.infType;
    primary = orig.primary;
    symptomatic = orig.symptomatic;
    IIP = orig.IIP;
}

//Infection::~Infection() {
//}
