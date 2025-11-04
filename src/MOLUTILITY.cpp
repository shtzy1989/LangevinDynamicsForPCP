#include "MOLUTILITY.h"
#include "DataLocator.h"
#include "STATISTICS.h"

namespace MOLUTILITY{
    void findNextBackUpFileName(std::string filename){
        if( file_exist(filename) ){
            int count = 1;
            while( true ){
                text path = text(filename).parsePath();
                text name = text(filename).parseFileName();
                char buffer[256];
                sprintf(buffer, "%d", count);
                std::string backUpFileName = path.string() + "#" + name.string() + "_bk." + std::string(buffer) + "#";
                if( !file_exist(backUpFileName) ){
                    rename(filename.c_str(), backUpFileName.c_str());
                    break;
                }
                count++;
            }
        }
    };

    bool file_exist(std::string filename){
        return access(filename.c_str(), F_OK) != -1;
    };

    void linearRegression(std::vector<double> *x, std::vector<double> *y, double& slope, double& intercept, double& R){
        int n = x->size();
        double Sx = STATISTICS::SUM(x);
        double Sy = STATISTICS::SUM(y);
        double xavg = Sx / n;
        double yavg = Sy / n;
        double Sxy = 0.0;
        double Sxx = 0.0;
        double Syy = 0.0;

        for (int i = 0; i < n; i++) {
            double xi = (*x)[i];
            double yi = (*y)[i];
            Sxy += (xi - xavg) * (yi - yavg); // Covariance
            Sxx += (xi - xavg) * (xi - xavg); // Variance of x
            Syy += (yi - yavg) * (yi - yavg); // Variance of y
        }

        slope = Sxy / Sxx;
        intercept = yavg - slope * xavg;
        R = Sxy / (sqrt(Sxx) * sqrt(Syy)); // Pearson correlation coefficient

    };

    bool InBound(double value, double lowBound, double highBound){
        return lowBound <= value && value < highBound;
    }
    bool InBound(int value, int lowBound, int highBound){
        return lowBound <= value && value < highBound;
    }

    std::vector<size_t> GenerateLog10Scale(double timestep, size_t numberOfFrame, size_t minimalSampleSize, size_t samplePerOrder, bool bIncludeZero){
        const double step = 1.0 / samplePerOrder;
        double maxTime = timestep * (numberOfFrame - minimalSampleSize);
        double timeOrder0 = log10(timestep);
        double timeOrder1 = log10(maxTime);
        std::vector<size_t> result;
        if( bIncludeZero ) result.push_back(0);
        for(double t=timeOrder0;t<=timeOrder1;t+=step){
            size_t interval = round(pow(10.0, t) / timestep);
            if( !result.size() || result.back() != interval ){
                result.push_back(interval);
            }
        }
        return result;
    }

    std::tuple<bool, double, double, double, int, int> AutoLinearFitting(std::vector<double>& msdTime, std::vector<double>& msdValue, bool bConvertToLog, double eps, int num, text label, int ncore, bool bPrint){
        // fitted, slope, intercept, correlation, leftBound, rightBound

        std::vector<double> msdSlope;
        std::vector<double> msdToFit;

        std::vector<double> msdLogTime;
        std::vector<double> msdLogValue;
        
        if( bConvertToLog ){
            msdLogTime.resize(msdTime.size());
            msdLogValue.resize(msdTime.size());

            #pragma omp parallel for num_threads(ncore)
            for(int i=0;i<msdTime.size();i++){
                msdLogTime[i] = std::log10(msdTime[i]);
                msdLogValue[i] = std::log10(msdValue[i]);
            }
        }

        msdSlope.resize(msdTime.size());
        msdToFit.resize(msdTime.size());

        int half = (num - 1) / 2;
        #pragma omp parallel for num_threads(ncore)
        for(int i=0;i<msdTime.size();i++){
            std::vector<double> x;
            std::vector<double> y;
            for(int j=i-num;j<=i+num;j++){
                if( InBound(j, 0, msdTime.size()) ){
                    if( bConvertToLog ){
                        x.push_back(msdLogTime[j]);
                        y.push_back(msdLogValue[j]);
                    }else{
                        x.push_back(msdTime[j]);
                        y.push_back(msdValue[j]);
                    }
                }
            }
            double slope, intercept, correl;
            linearRegression(&x, &y, slope, intercept, correl);
            msdSlope[i] = slope;

            if( std::abs(slope - 1.0) < eps ){
                msdToFit[i] = true;
            }else{
                msdToFit[i] = false;
            }
        }

        // find groups (segments of continuous points)
        std::vector<dataLocator::Element2Directional<int> > group;
        int s = -1;
        for(int i=0;i<msdTime.size();i++){
            if( msdToFit[i] && s == -1 ){
                s = i;
            }else if( s != -1 && !msdToFit[i] ){
                group.push_back(dataLocator::Element2Directional<int>(s, i-1));
                s = -1;
            }else if( s != -1 && i == msdTime.size() - 1 ){
                group.push_back(dataLocator::Element2Directional<int>(s, i));
                s = -1;
            }
        }

        bool bFitted = true;
        double resultSlope = 0.0;
        double resultIntercept = 0.0;
        double resultCorrel = 0.0;
        int leftBound = -1;
        int rightBound = -1;
        // for(int i=0;i<msdTime.size();i++){
        //     fprintf(stderr, "%3d %23.16E %23.16E %23.16E %23.16E %23.16E %s\n", i, msdTime[i], msdValue[i], msdLogTime[i], msdLogValue[i], msdSlope[i], msdToFit[i] ? "True" : "");
        // }

        // for(int i=0;i<group.size();i++){
        //     fprintf(stderr, "Group %d: %d - %d\n", i+1, group[i].Value(0), group[i].Value(1));
        // }

        if( group.size() == 0 ){
            if( bPrint ){
                fprintf(stderr, "%23s %s\n", "Not ready", label.c_str());
            }
            bFitted = false;
        }else{
            // find the largest group
            int groupIndex = -1;
            double groupSize = 0;
            for(int i=0;i<group.size();i++){
                double groupSize0 = msdLogTime[group[i].Value(1)] - msdLogTime[group[i].Value(0)];
                if( groupSize0 > groupSize ){
                    groupSize = groupSize0;
                    groupIndex = i;
                }
            }

            if( groupIndex == -1 ){
                bFitted = false;
            }else{
                bFitted = true;
                
                std::vector<double> x;
                std::vector<double> y;
                for(int i=group[groupIndex].Value(0);i<=group[groupIndex].Value(1);i++){
                    x.push_back(msdTime[i]);
                    y.push_back(msdValue[i]);
                }

                linearRegression(&x, &y, resultSlope, resultIntercept, resultCorrel);
                
                leftBound = group[groupIndex].Value(0);
                rightBound = group[groupIndex].Value(1);
            }
        }

        return std::make_tuple(bFitted, resultSlope, resultIntercept, resultCorrel, leftBound, rightBound);
    }
};