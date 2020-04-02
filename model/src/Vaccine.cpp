#include "Vaccine.h"
#include <map>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <string>


Vaccine::Vaccine(){
    vaccineID = -1;
    mode = "none";
    name = "none";
    waning = 0.0;
    protection = 0.0;
    propInf = 0.0;
    normdev_pos = 0.0;
    normdev_neg = 0.0;
    seroposVE = 0.0;
    seronegVE = 0.0;
    seroposWaning = 0.0;
    seronegWaning = 0.0;
    for(unsigned s = 0; s < 4; s++){
	RRInf_seropos[s] = 1.0;
	RRInf_seroneg[s] = 1.0;
	RRDis_seropos[s] = 1.0;
	RRDis_seroneg[s] = 1.0;
	RRHosp_seropos[s] = 1.0;
	RRHosp_seroneg[s] = 1.0;
    }
    for(unsigned k = 0; k < 3; k++){
	VE_pos.insert(std::make_pair(k,0.0));
	VE_neg.insert(std::make_pair(k,0.0));
    }
    total_VE = 0.0;
    doses = 1;
    relative_schedule.clear();
}

void Vaccine::init(){
    vaccineID = -1;
    mode = "none";
    name = "none";
    waning = 0.0;
    protection = 0.0;
    propInf = 0.0;
    normdev_pos = 0.0;
    normdev_neg = 0.0;
    seroposVE = 0.0;
    seronegVE = 0.0;
    seroposWaning = 0.0;
    seronegWaning = 0.0;
    for(unsigned s = 0; s < 4; s++){
	RRInf_seropos[s] = 1.0;
	RRInf_seroneg[s] = 1.0;
	RRDis_seropos[s] = 1.0;
	RRDis_seroneg[s] = 1.0;
	RRHosp_seropos[s] = 1.0;
	RRHosp_seroneg[s] = 1.0;
    }
    for(unsigned k = 0; k < 3; k++){
	VE_pos.insert(std::make_pair(k,0.0));
	VE_neg.insert(std::make_pair(k,0.0));
    }
    total_VE = 0.0;
    doses = 1;
    relative_schedule.clear();
}

int Vaccine::getNextDoseTime(int vday, int received){
    if(relative_schedule.size() == 0){
	printf("Relative schedule is empty id: %d\n",vaccineID);
	return vday;
    }
    if(received < doses){
	return (relative_schedule[received] + vday);
    }else{
	return (relative_schedule.back() + vday);
    }
}
double Vaccine::getRR(double prevInf, double ageDays){
    double normdev = 0;
    if(prevInf > 0){
	normdev = normdev_pos;
	total_VE = 1.0 - VE_pos.at(0) / (1.0 + exp(VE_pos.at(1) * (ageDays / 365.0 - VE_pos.at(2))));
    }else{
	normdev = normdev_neg;
	total_VE = 1.0 - VE_neg.at(0) / (1.0 + exp(VE_neg.at(1) * (ageDays / 365.0 - VE_neg.at(2))));
    }    
    double RR = exp(log(1.0 - total_VE) + normdev * pow(1.0 / 100.5 + 1.0 / (100.0 * (1.0 - total_VE) + 0.5), 0.5));
    //    printf("getRR: total_VE: %.2f RR: %.2f previnf: %.0f\n", total_VE, RR, prevInf);
    return RR;
}

void Vaccine::printVaccine(){
    printf("Vaccine ID: %d Name: %s Mode: %s Waning %.2f Protection %.2f Efficacy %.2f PropInf %.2f NormDevPos %.2f NormDevNeg %.2f Doses %d\n",vaccineID, name.c_str(), mode.c_str(), waning,protection,total_VE,
	   propInf, normdev_pos,normdev_neg, doses);
    printf("VE seropos: %.2f VE seroneg %.2f Waning seropos %.2f Waning seroneg %.2f\n",seroposVE,seronegVE,seroposWaning,seronegWaning);
    if(mode == "GSK"){
	for(unsigned s = 0; s < 4; s++){
	    printf("serotype: %u. RRInf: pos %.2f neg %.2f RRDis: pos %.2f neg %.2f RRHosp: pos %.2f neg %.2f\n",s,RRInf_seropos[s],RRInf_seroneg[s],RRDis_seropos[s],RRDis_seroneg[s],RRHosp_seropos[s],RRHosp_seroneg[s]);
	}
    }
    if(mode == "age"){
	for(unsigned j = 0; j < 3; j++){
	    printf("PAR %d VE pos %.2f VE neg %.2f\n",j,VE_pos.at(j), VE_neg.at(j));
	}
    }
    for(unsigned j = 0; j < relative_schedule.size(); j++){
	printf("Dose %u Day %d\t", j + 1, relative_schedule[j]);
    }
    printf("\n");
}

void Vaccine::setRelativeSchedule(vector<int> rs){
    relative_schedule = rs;
    if(relative_schedule.size() > 0){
	doses = relative_schedule.size();
    }
}

void Vaccine::setWaning(bool seroposIn, double w){
    if(seroposIn){
	seroposWaning = w;
    }else{
	seronegWaning = w;
    }
}

double Vaccine::getWaning(bool seroposIn){
    if(seroposIn){
	return seroposWaning;
    }else{
	return seronegWaning;
    }
}

void Vaccine::setVaccineEfficacy(bool seroposIn, double veIn){
    if(seroposIn){
	seroposVE = veIn;
    }else{
	seronegVE = veIn;
    }
}
void Vaccine::setRRInf(bool seroposIn, double rr, unsigned sero_){
    if(seroposIn){
	RRInf_seropos[sero_] = rr;
    }else{
	RRInf_seroneg[sero_] = rr;
    }
}

void Vaccine::setRRDis(bool seroposIn, double rr, unsigned sero_){
    if(seroposIn){
	RRDis_seropos[sero_] = rr;
    }else{
	RRDis_seroneg[sero_] = rr;
    }
}

void Vaccine::setRRHosp(bool seroposIn, double rr, unsigned sero_){
    if(seroposIn){
	RRHosp_seropos[sero_] = rr;
    }else{
	RRHosp_seroneg[sero_] = rr;
    }
}

void Vaccine::setRRInf(bool seroposIn, double rr){
    for(unsigned sero_ = 0; sero_ < 4; sero_++){
	if(seroposIn){
	    RRInf_seropos[sero_] = rr;
	}else{
	    RRInf_seroneg[sero_] = rr;
	}
    }
}

void Vaccine::setRRDis(bool seroposIn, double rr){
    for(unsigned sero_ = 0; sero_ < 4; sero_++){
	if(seroposIn){
	    RRDis_seropos[sero_] = rr;
	}else{
	    RRDis_seroneg[sero_] = rr;
	}
    }
}

void Vaccine::setRRHosp(bool seroposIn, double rr){
    for(unsigned sero_ = 0; sero_ < 4; sero_++){
	if(seroposIn){
	    RRHosp_seropos[sero_] = rr;
	}else{
	    RRHosp_seroneg[sero_] = rr;
	}
    }
}


double Vaccine::getVaccineEfficacy(bool seroposIn){
    if(seroposIn){
	return seroposVE;
    }else{
	return seronegVE;
    }
}

double Vaccine::getRRInf(bool seroposIn, unsigned sero_){
    if(seroposIn){
	return RRInf_seropos[sero_];
    }else{
	return RRInf_seroneg[sero_];
    }    
}

double Vaccine::getRRDis(bool seroposIn, unsigned sero_){
    if(seroposIn){
	return RRDis_seropos[sero_];
    }else{
	return RRDis_seroneg[sero_];
    }    
}

double Vaccine::getRRHosp(bool seroposIn, unsigned sero_){
    if(seroposIn){
	return RRHosp_seropos[sero_];
    }else{
	return RRHosp_seroneg[sero_];
    }    
}

//Vaccine::~Vaccine() {
//}

