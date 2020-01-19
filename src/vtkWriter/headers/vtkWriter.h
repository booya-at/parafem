#ifndef vtkWriter_H
#define vtkWriter_H

#include "case.h"

#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <string>


using namespace std;

namespace parafem{

class VtkWriter{
private:
    string file_name;
    int count = 0;
public:
    VtkWriter(const char* file_name);
    int writeCase(parafem::FemCasePtr c, double coordSysSize=0);
};


};

#endif