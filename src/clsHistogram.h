#ifndef MD_HISTOGRAM_H
#define MD_HISTOGRAM_H

#include "cppheader.h"
#include "streamer.h"

namespace HISTOGRAM{
    std::tuple<double, double> ExtendRange(double bound0, double bound1, double scale = 1.05, double minimalRange = -1.0);

    class HistogramSetting{
    public:
        union{
            double m_HistHSetting[4];
            struct{
                double m_HistHMin, m_HistHMax, m_HistHBin, m_HistHBinInverse;
            };
        };
        int m_HistNBin;
    public:
        HistogramSetting();
        HistogramSetting(double hmin, double hmax, double hbin);
        HistogramSetting(double hmin, double hmax, int nbin);
    public:
        void Set(double hmin, double hmax, double hbin);
        void Set(double hmin, double hmax, int nbin);
        int GetIndex(double value);
        bool InBound(double value);
        void MeasureFromVector(std::vector<double>* data, int nbin, double extendRange = 1.0);
        void MeasureFromVector(std::vector<double>* data, double hbin, int minNBin = 10, double extendRange = 1.0);
    public:
        void Print(FILE* fout = stderr);
    };
    bool operator== (const HistogramSetting& lhs, const HistogramSetting& rhs);
    bool operator!= (const HistogramSetting& lhs, const HistogramSetting& rhs);

    struct HISTOGRAM{
    public:
        double min;
        double max;
        double step;
        int size;
        double stepInverse;
        std::vector<double> count;
        std::vector<double> percent;
    public:
        void NormalizePercent();
        void NormalizePercent(double sum);
        double GetTotalCount();
        void Print(FILE* fout = stderr);
        int GetNBin();
        double GetHBin();
        void SetNBin(int value);
        void SetHBin(double value);
        void Allocate(double histMin, double histMax, double histHbin);
        void Allocate(double histMin, double histMax, int histNbin);
        void Allocate(HistogramSetting histogramSetting);
        void CopyFrom(HISTOGRAM& rhs);
        bool ResizeToFitNewBound(double value, size_t* zeroIndexInNewHistogram = 0);
        bool ReBoundBySetting(HistogramSetting *setting, int ncore = 1, double eps = 0.000001, bool bPrint = true);
        void HistogramVector(std::vector<double>* data, int first = 0, int last = -1, int freq = 1);
        bool AddValue(double value);
    public:
        HISTOGRAM& operator= (HISTOGRAM& rhs);
    };

    class HistogramMultiple{
    public:
        HistogramSetting m_HistogramSetting;
        std::vector<std::vector<double> > m_Data;
    public:
        void NormalizePercent(int index, double sum = 1.0);
        void NormalizePercent(double sum = 1.0);
        double GetSum(int index);
        void Print(FILE* fout = stderr);
        void Allocate(int numberOfData, HistogramSetting histogramSetting);
        void Allocate(int numberOfData, double histMin, double histMax, double histHbin);
        void Allocate(int numberOfData, double histMin, double histMax, int histNbin);
        void InitializeValue(int index, double value);
        void CopyFrom(HistogramMultiple& rhs);
        // 0 = +, 1 = min, 2 = max
        bool SumOMP(std::vector<HistogramMultiple>& rhs, std::vector<int> operation);
        double HistHBin();
        double HistHBinInverse();
        double HistHMin();
        double HistHMax();
        int HistNBin();
        bool AddValue(double valueBin, double value, int dataIndex);
        bool AddValue(int binIndex, double value, int dataIndex);
        void AddValueNoCheck(int binIndex, double value, int dataIndex);
        bool MaxValue(double valueBin, double value, int dataIndex);
        bool MaxValue(int binIndex, double value, int dataIndex);
        bool MinValue(double valueBin, double value, int dataIndex);
        bool MinValue(int binIndex, double value, int dataIndex);
        int GetValueIndex(double value);
        bool InBound(double value);
        bool InBound(int value);
        void HistogramVector(std::vector<double>* data, int dataIndex, int first = 0, int last = -1, int freq = 1);
    public:
        HistogramMultiple& operator= (HistogramMultiple& rhs);
    };

};

#endif
