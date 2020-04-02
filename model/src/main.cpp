#include <cstdlib>
//#include <random>
#include <iostream>
#include <string>
#include <ctime>
#include <fstream>
#include <chrono>
#include "ThreadPool.h"
#include "Simulation.h"

//using namespace std;
using std::chrono::nanoseconds;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

int main(int argc, char** argv) {
  //srand(time(0));
  char *configFileName;
  if (argc != 2 && argc != 3 && argc != 5) { // Check for number of input arguments
    exit(1); // Exit if not enough arguments
  } else if (argc == 2) {
    std::chrono::time_point<high_resolution_clock> begin, end, newSimBegin;
    begin = high_resolution_clock::now();
    configFileName = argv [1];

    string conFile(configFileName);
    //double min, sec;

    std::ifstream infile(conFile);
    if (!infile.good()) {
      exit(1);
    }
    string line;
    getline(infile, line);
    while (getline(infile, line)) {
      newSimBegin = high_resolution_clock::now();
      Simulation sim(line);
      string simName = sim.readInputs();
      end = high_resolution_clock::now();
      //min = duration_cast<nanoseconds> (end-newSimBegin).count()/ 60000000000.;
      //sec = (min - (int)min)*60.0;

      sim.simulate();
      end = high_resolution_clock::now();
      //min = duration_cast<nanoseconds> (end-newSimBegin).count()/ 60000000000.;
      //sec = (min - (int)min)*60.0;
      while (infile.peek() == '\n') {
        infile.ignore(1, '\n');
      }

    }
    infile.close();

    end = high_resolution_clock::now();
    //min = duration_cast<nanoseconds> (end-begin).count()/ 60000000000.;
    //sec = (min - (int)min)*60.0;


  } else if (argc == 3) { // to parallelize
    printf("arguments = 3\n");
    std::chrono::time_point<high_resolution_clock> begin, end, newSimBegin;
    begin = high_resolution_clock::now();
    unsigned maxThreads = strtol(argv[2], NULL, 10);
    configFileName = argv [1];
    Simulation sim(configFileName);
    string conFile(configFileName);

    std::ifstream infile(conFile);
    if (!infile.good()) {
      exit(1);
    }
    string line;
    getline(infile, line);

    ThreadPool pool(maxThreads);
    std::vector< std::future<int> > results;

    while (getline(infile, line)) {
      results.push_back(
        pool.enqueue([line,&infile,&begin] {
          std::chrono::time_point<high_resolution_clock> newSimBegin, newSimEnd;
          newSimBegin = high_resolution_clock::now();
          Simulation sim(line);
          string simName = sim.readInputs();
          newSimEnd = high_resolution_clock::now();
          double min = duration_cast<nanoseconds> (newSimEnd-newSimBegin).count()/ 60000000000.;
          //double sec = (min - (int)min)*60.0;
          sim.simulate();
          newSimEnd = high_resolution_clock::now();
          min = duration_cast<nanoseconds> (newSimEnd-newSimBegin).count()/ 60000000000.;
          //sec = (min - (int)min)*60.0;
          while (infile.peek() == '\n') {
            infile.ignore(1, '\n');
          }
          std::cout<< min << std::endl;
          return 0;
        })
      );
    }
    infile.close();

    for(size_t i = 0;i<results.size();++i) {
      results[i].get();
    }
    end = high_resolution_clock::now();
    double min = duration_cast<nanoseconds> (end-begin).count()/ 60000000000.;
    //double sec = (min - (int)min)*60.0;
    std::cout<< min << std::endl;

  } else if (argc == 5) {
    if (string(argv[1]) != "-d") {
      exit(1);
    }
    //Simulation sim;
  }
  return 0;
}
