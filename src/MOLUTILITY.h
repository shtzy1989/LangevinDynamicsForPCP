#ifndef MD_MOLUTILITY_H
#define MD_MOLUTILITY_H

#include "cppheader.h"
#include "text.h"

namespace MOLUTILITY{
    void findNextBackUpFileName(std::string filename);
    bool file_exist(std::string filename);

    void linearRegression(std::vector<double> *x, std::vector<double> *y, double& slope, double& intercept, double& R);

    bool InBound(double value, double lowBound, double highBound);
    bool InBound(int value, int lowBound, int highBound);
    std::vector<size_t> GenerateLog10Scale(double timestep, size_t numberOfFrame, size_t minimalSampleSize, size_t samplePerOrder, bool bIncludeZero);

    std::tuple<bool, double, double, double, int, int> AutoLinearFitting(std::vector<double>& msdTime, std::vector<double>& msdValue, 
        bool bConvertToLog = false, double eps = 0.01, int num = 3, text label = "", int ncore = 1, bool bPrint = true);


};

#endif