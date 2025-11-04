
#include "STATISTICS.h"

namespace STATISTICS{
    double MIN(std::vector<double>* in){
        if( in->size() ){
            double _min = DBL_MAX;
            for(int i=0;i<in->size();i++){
                _min = std::min<double>((*in)[i], _min);
            }
            return _min;
        }else{
            return 0.0;
        }
    };
    double MAX(std::vector<double>* in){
        if( in->size() ){
            double _max = -DBL_MAX;
            for(int i=0;i<in->size();i++){
                _max = std::max<double>((*in)[i], _max);
            }
            return _max;
        }else{
            return 0.0;
        }
    };
    int MINOFX(std::vector<double>* in){
        if( in->size() ){
            double _min = DBL_MAX;
            int index = -1;
            for(int i=0;i<in->size();i++){
                if( (*in)[i] < _min ){
                    _min = (*in)[i];
                    index = i;
                }
            }
            return index;
        }else{
            return 0.0;
        }
    };
    int MAXOFX(std::vector<double>* in){
        if( in->size() ){
            double _max = -DBL_MAX;
            int index = -1;
            for(int i=0;i<in->size();i++){
                if( (*in)[i] > _max ){
                    _max = (*in)[i];
                    index = i;
                }
            }
            return index;
        }else{
            return 0.0;
        }
    };
    int MIN(std::vector<int>* in){
        if( in->size() ){
            int _min = INT_MAX;
            for(int i=0;i<in->size();i++){
                _min = std::min<int>((*in)[i], _min);
            }
            return _min;
        }else{
            return 0.0;
        }
    };
    int MAX(std::vector<int>* in){
        if( in->size() ){
            int _max = -INT_MAX;
            for(int i=0;i<in->size();i++){
                _max = std::max<int>((*in)[i], _max);
            }
            return _max;
        }else{
            return 0.0;
        }
    };
    double SUM(std::vector<double>* in){
        double sum = 0;
        for(int i=0;i<in->size();i++){
            sum += (*in)[i];
        }
        return sum;
    };
    double SUM2(std::vector<double>* in){
        double sum = 0;
        for(int i=0;i<in->size();i++){
            sum += (*in)[i] * (*in)[i];
        }
        return sum;
    };
    int SUM(std::vector<int>* in){
        int sum = 0;
        for(int i=0;i<in->size();i++){
            sum += (*in)[i];
        }
        return sum;
    };
    double AVERAGE(std::vector<double>* in){
        if( in->size() ){
            return SUM(in) / double(in->size());
        }else{
            return 0;
        }
    };
    double AVERAGE2(std::vector<double>* in){
        if( in->size() ){
            return SUM2(in) / double(in->size());
        }else{
            return 0;
        }
    };
    void CUMMULATIVE(std::vector<double>* out, std::vector<double>* in){
        out->resize(in->size());
        double sum =0;
        for(int i=0;i<in->size();i++){
            sum += (*in)[i];
            (*out)[i] = sum / double(i+1);
        }
    };
    double RMSD(std::vector<double>* in){
        if( in->size() ){
            double rmsd = 0;
            double average = AVERAGE(in);
            for(int i=0;i<in->size();i++){
                rmsd += ((*in)[i] - average) * ((*in)[i] - average);
            }
            rmsd = sqrt(rmsd / in->size());
            return rmsd;
        }else{
            return 0;
        }
    };
    double RMSD2(std::vector<double>* in){
        if( in->size() ){
            double sum = 0.0;
            double sum2 = 0.0;
            for(int i=0;i<in->size();i++){
                sum += (*in)[i];
                sum2 += (*in)[i] * (*in)[i];
            }
            double avg = sum / double(in->size());
            double avg2 = sum2 / double(in->size());
            return sqrt(avg2 - avg * avg);
        }else{
            return 0;
        }
    };
    std::tuple<double, double> AVGERAGEANDRMSD(std::vector<double>* in){
        if( in->size() ){
            double sum = 0.0;
            double sum2 = 0.0;
            for(int i=0;i<in->size();i++){
                sum += (*in)[i];
                sum2 += (*in)[i] * (*in)[i];
            }
            double avg = sum / double(in->size());
            double avg2 = sum2 / double(in->size());
            return std::tuple(avg, sqrt(avg2 - avg * avg));
        }else{
            return std::tuple(0, 0);
        }
    }
    std::tuple<double, double> MINANDMAX(std::vector<double>* in){
        if( in->size() ){
            double minValue = DBL_MAX;
            double maxValue = 0.0;
            for(int i=0;i<in->size();i++){
                minValue = std::min<double>(minValue, (*in)[i]);
                maxValue = std::max<double>(maxValue, (*in)[i]);
            }
            return std::tuple(minValue, maxValue);
        }else{
            return std::tuple(0, 0);
        }
    }
    double SD(std::vector<double>* in){
        if( in->size() - 1 > 0 ){
            double rmsd = 0;
            double average = AVERAGE(in);
            for(int i=0;i<in->size();i++){
                rmsd += ((*in)[i] - average) * ((*in)[i] - average);
            }
            rmsd = sqrt(rmsd / (in->size() - 1));
            return rmsd;
        }else{
            return 0;
        }
    };
    double SDM(std::vector<double>* in){
        double rmsd = 0;
        double average = AVERAGE(in);
        for(int i=0;i<in->size();i++){
            rmsd += ((*in)[i] - average) * ((*in)[i] - average);
        }
        rmsd = sqrt(rmsd) / double(in->size());
        return rmsd;
    };
    int SUBARRAY(std::vector<double>* out, std::vector<double>* in, int begin, int length){
        if( begin + length > in->size() ){
            length = in->size() - begin;
        }
        out->resize(length);
        for(int i=0;i<length;i++){
            (*out)[i] = (*in)[i+begin];
        }
        return length;
    };
    void PRINT(std::vector<double>* in){
        for(int i=0;i<in->size();i++){
            fprintf(stderr, "v[%-4d] = %13.7f\n", i, (*in)[i]);
        }
    }
    void COMBINEARRAY(std::vector<double>* out, std::vector<double>* in){
        for(int i=0;i<in->size();i++){
            out->push_back((*in)[i]);
        }
    }
    void COMBINEARRAY(std::vector<double>* out, std::vector<double>* in1, std::vector<double>* in2){
        out->clear();
        for(int i=0;i<in1->size();i++){
            out->push_back((*in1)[i]);
        }
        for(int i=0;i<in2->size();i++){
            out->push_back((*in2)[i]);
        }
    }
    double VARIANCE(std::vector<double>* in){
        if( in->size() > 0 ){
            double variance = 0;
            double average = AVERAGE(in);
            for(int i=0;i<in->size();i++){
                double deviation = ((*in)[i] - average);
                variance += deviation * deviation;
            }
            return variance / (in->size() - 1);
        }else{
            return 0;
        }
    };
    double MOMENT4(std::vector<double>* in){
        if( in->size() > 0 ){
            double moment4 = 0;
            double average = AVERAGE(in);
            for(int i=0;i<in->size();i++){
                double deviation = ((*in)[i] - average);
                double deviation2 = deviation * deviation;
                moment4 += deviation2 * deviation2;
            }
            return moment4 / (in->size());
        }else{
            return 0;
        }
    };
    std::tuple<double, double, double, double>  STATISTICSET(std::vector<double>* in){
        if( in->size() ){
            double _min = DBL_MAX;
            double _max = -DBL_MAX;
            double sum = 0;
            double sum2 = 0;

            for(int i=0;i<in->size();i++){
                _min = std::min<double>((*in)[i], _min);
                _max = std::max<double>((*in)[i], _max);
                sum += (*in)[i];
                sum2 += (*in)[i] * (*in)[i];
            }
            double avg = sum / in->size();
            double avg2 = sum2 / in->size();
            double rmsd = sqrt(avg2 - avg * avg);
            return std::tuple<double, double, double, double>(avg, _min, _max, rmsd);
        }else{
            return std::tuple<double, double, double, double>(0.0, 0.0, 0.0, 0.0);
        }
    }
    void SUBTRACTAVERAGE(std::vector<double>* in){
        double avg = AVERAGE(in);
        for(int i=0;i<in->size();i++){
            (*in)[i] -= avg;
        }
    }

    double Average(std::vector<double>* in, int first, int last){
        if( last < 0 || last > in->size() ) last = in->size();
        if( first < 0 || first >= in->size() ) first = 0;
        int size = last - first;
        if( size != 0 ){
            double value = 0.0;
            for(int i=first;i<last;i++){
                value += (*in)[i];
            }
            value /= double(size);
            return value;
        }else{
            return 0;
        }
    }

    std::tuple<double, double> AverageAndSD(std::vector<double>* in){
        if( in->size() ){
            double rmsd = 0;
            double average = AVERAGE(in);
            for(int i=0;i<in->size();i++){
                rmsd += ((*in)[i] - average) * ((*in)[i] - average);
            }
            rmsd = sqrt(rmsd / in->size());
            return { average, rmsd };
        }else{
            return { 0.0, 0.0 };
        }
    }

    double CorrelationCoefficient(std::vector<double>* x, std::vector<double>* y){
        double avgX = Average(x);
        double avgY = Average(y);
        double cov = 0.0;
        double sdX = 0.0;
        double sdY = 0.0;
        for(int i=0;i<x->size();i++){
            cov += ((*x)[i] - avgX) * ((*y)[i] - avgY);
            sdX += ((*x)[i] - avgX) * ((*x)[i] - avgX);
            sdY += ((*y)[i] - avgY) * ((*y)[i] - avgY);
        }
        return cov / std::sqrt(sdX * sdY);
        
    }

    double mean(const std::vector<double>& data) {
        return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    }

    // Function to compute variance
    double variance(const std::vector<double>& data, double mean) {
        double sum = 0.0;
        for (double x : data) {
            sum += (x - mean) * (x - mean);
        }
        return sum / data.size();
    }

    double skewness(const std::vector<double>& data, double mean, double std_dev) {
        double sum = 0.0;
        for (double x : data) {
            sum += std::pow((x - mean) / std_dev, 3);
        }
        return sum / data.size();
    }

    // Function to compute kurtosis
    double kurtosis(const std::vector<double>& data, double mean, double std_dev) {
        double sum = 0.0;
        for (double x : data) {
            sum += std::pow((x - mean) / std_dev, 4);
        }
        return sum / data.size() - 3.0;  // Excess kurtosis
    }

    double chi_square_p_value(double chi_square_stat, int df) {
        // Approximate p-value using the Gamma function
        double k = df / 2.0;
        double gamma_k = std::tgamma(k);
        double x = chi_square_stat / 2.0;
        
        // Series expansion approximation of the incomplete gamma function
        double sum = 1.0, term = 1.0;
        for (int i = 1; i < 100; ++i) {
            term *= x / (k + i);
            sum += term;
            if (term < 1e-10) break;
        }
        
        double p_value = std::exp(-x) * std::pow(x, k) / gamma_k * sum;
        return 1.0 - p_value;  // Survival function
    }

    double DagostinoPearsonTest(const std::vector<double>& data){
        double mu = mean(data);
        double var = variance(data, mu);
        double sigma = std::sqrt(var);

        // Compute skewness and kurtosis
        double skew = skewness(data, mu, sigma);
        double kurt = kurtosis(data, mu, sigma);

        // Compute D'Agostino & Pearson's test statistic
        double test_stat = skew * skew + kurt * kurt / 4.0;
        
        // Get p-value from chi-square distribution with df=2
        double p_value = chi_square_p_value(test_stat, 2);

        return p_value;  // Return p-value
    }

    double empirical_cdf(const std::vector<double>& data, double x) {
        return std::count_if(data.begin(), data.end(), [&](double v) { return v <= x; }) / (double)data.size();
    }

    // Gaussian CDF approximation (using error function)
    double gaussian_cdf(double x, double mean, double stddev) {
        return 0.5 * (1 + std::erf((x - mean) / (stddev * std::sqrt(2))));
    }

    // Kolmogorov-Smirnov Test Statistic
    double ks_test_statistic(const std::vector<double>& data, double mean, double stddev) {
        double max_diff = 0.0;
        for (double x : data) {
            double ecdf = empirical_cdf(data, x);
            double gaussian_cdf_value = gaussian_cdf(x, mean, stddev);
            max_diff = std::max(max_diff, std::abs(ecdf - gaussian_cdf_value));
        }
        return max_diff;
    }

    // Approximate K-S p-value using Smirnov's formula
    double ks_p_value(double ks_stat, int n) {
        double sqrt_n = std::sqrt(n);
        double lambda = (sqrt_n + 0.12 + 0.11 / sqrt_n) * ks_stat;
        double p = 2 * std::exp(-2 * lambda * lambda);
        return (p > 1.0) ? 1.0 : p; // Ensure p-value is at most 1
    }

    double KSTest(const std::vector<double>& data){
        double mu = mean(data);
        double var = variance(data, mu);
        double stddev = std::sqrt(var);

        // Compute K-S test statistic and p-value
        double ks_stat = ks_test_statistic(data, mu, stddev);
        double p_value = ks_p_value(ks_stat, data.size());

        return p_value;  // Return p-value
    }
};
