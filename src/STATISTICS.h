#ifndef MD_STATISTICS_H
#define MD_STATISTICS_H

#include "cppheader.h"

namespace STATISTICS{
    double MIN(std::vector<double>* in);
    double MAX(std::vector<double>* in);
    int MINOFX(std::vector<double>* in);
    int MAXOFX(std::vector<double>* in);
    double SUM(std::vector<double>* in);
    double SUM2(std::vector<double>* in);
    int MIN(std::vector<int>* in);
    int MAX(std::vector<int>* in);
    int SUM(std::vector<int>* in);
    double AVERAGE(std::vector<double>* in);
    double AVERAGE2(std::vector<double>* in);
    void CUMMULATIVE(std::vector<double>* out, std::vector<double>* in);
    double RMSD(std::vector<double>* in);
    double RMSD2(std::vector<double>* in);
    //         avg     rmsd     
    std::tuple<double, double> AVGERAGEANDRMSD(std::vector<double>* in);
    //         min     max     
    std::tuple<double, double> MINANDMAX(std::vector<double>* in);
    double SD(std::vector<double>* in);
    double SDM(std::vector<double>* in);
    int SUBARRAY(std::vector<double>* out, std::vector<double>* in, int begin, int length);
    void PRINT(std::vector<double>* in);
    void COMBINEARRAY(std::vector<double>* out, std::vector<double>* in);
    void COMBINEARRAY(std::vector<double>* out, std::vector<double>* in1, std::vector<double>* in2);
    double VARIANCE(std::vector<double>* in);
    double MOMENT4(std::vector<double>* in);    
    //         avg     min     max     std
    std::tuple<double, double, double, double> STATISTICSET(std::vector<double>* in);
    void SUBTRACTAVERAGE(std::vector<double>* in);

    double Average(std::vector<double>* in, int first = 0, int last = -1);

    std::tuple<double, double> AverageAndSD(std::vector<double>* in);

    double CorrelationCoefficient(std::vector<double>* x, std::vector<double>* y);

    double mean(const std::vector<double>& data);
    double variance(const std::vector<double>& data, double mean);
    double skewness(const std::vector<double>& data, double mean, double std_dev);
    double kurtosis(const std::vector<double>& data, double mean, double std_dev);
    double chi_square_p_value(double chi_square_stat, int df);
    double DagostinoPearsonTest(const std::vector<double>& data);

    double empirical_cdf(const std::vector<double>& data, double x);
    double gaussian_cdf(double x, double mean, double stddev);
    double ks_test_statistic(const std::vector<double>& data, double mean, double stddev);
    double ks_p_value(double ks_stat, int n);
    double KSTest(const std::vector<double>& data);

    template<typename T>
    class EntryCounter{
    protected:
        class Entry{
        public:
            T m_Key;
            int m_Index;
        public:
            Entry(T key, int index){
                m_Key = key;
                m_Index = index;
            };
        public:
            bool operator<(const Entry& rhs) const{
                return m_Key < rhs.m_Key;
            };
        };
        std::set<Entry> m_Set;
        std::vector<int> m_Count;
        // 
        std::vector<std::pair<T, int> > m_SortedData;
    public:
        void Add(T value, int count = 1){
            Entry newEntry(value, m_Count.size());
            auto result = m_Set.insert(newEntry);
            if( result.second ){
                // did not exist 
                m_Count.push_back(count);
            }else{
                // already exists
                m_Count[result.first->m_Index] += count;
            }
        }
        void Sort(){
            m_SortedData.clear();
            for(auto it : m_Set){
                int count = m_Count[it.m_Index];
                m_SortedData.push_back(std::make_pair(it.m_Key, count));
            }
        }
        int GetNumberOfEntry(){
            if( m_SortedData.size() == 0 ){
                fprintf(stderr, "Error: Sort() must be called first. \n");
            }
            return m_SortedData.size();
        }
        T GetEntry(int index){
            if( m_SortedData.size() == 0 ){
                fprintf(stderr, "Error: Sort() must be called first. \n");
            }
            return m_SortedData[index].first;
        }
        int GetCount(int index){
            if( m_SortedData.size() == 0 ){
                fprintf(stderr, "Error: Sort() must be called first. \n");
            }
            return m_SortedData[index].second;
        }
        void AddFrom(EntryCounter& rhs){
            if( rhs.m_SortedData.size() == 0 ){
                fprintf(stderr, "Error: Sort() must be called first for RHS. \n");
            }
            m_SortedData.clear();
            for(auto it : rhs.m_SortedData){
                Add(it.first, it.second);
            }
        }
    };
};

#endif
