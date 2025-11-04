#include "clsHistogram.h"
#include "MOLUTILITY.h"
#include "STATISTICS.h"

std::tuple<double, double> HISTOGRAM::ExtendRange(double bound0, double bound1, double scale, double minimalRange){
    double minValue = std::min<double>(bound0, bound1);
    double maxValue = std::max<double>(bound0, bound1);

    double range = maxValue - minValue;
    range = std::max<double>(range, minimalRange);

    double mean = (minValue + maxValue) * 0.5;
    double half = (range) * 0.5 * scale;

    return std::tuple<double, double>(mean - half, mean + half);
}   

HISTOGRAM::HistogramSetting::HistogramSetting(){

};

HISTOGRAM::HistogramSetting::HistogramSetting(double hmin, double hmax, double hbin){
    Set(hmin, hmax, hbin);
};

HISTOGRAM::HistogramSetting::HistogramSetting(double hmin, double hmax, int nbin){
    Set(hmin, hmax, nbin);
};

void HISTOGRAM::HistogramSetting::Set(double hmin, double hmax, double hbin){
    m_HistHMin = hmin;
    m_HistHMax = hmax;
    m_HistNBin = round((hmax - hmin) / hbin);
    m_HistHBin = (hmax - hmin) / m_HistNBin;
    m_HistHBinInverse = 1.0 / m_HistHBin;
};

void HISTOGRAM::HistogramSetting::Set(double hmin, double hmax, int nbin){
    m_HistHMin = hmin;
    m_HistHMax = hmax;
    m_HistNBin = nbin;
    m_HistHBin = (hmax - hmin) / m_HistNBin;
    m_HistHBinInverse = 1.0 / m_HistHBin;
};

int HISTOGRAM::HistogramSetting::GetIndex(double value){
    int index = (value - m_HistHMin) * m_HistHBinInverse;
    return index;
};

bool HISTOGRAM::HistogramSetting::InBound(double value){
    return MOLUTILITY::InBound(value, m_HistHMin, m_HistHMax);
};

void HISTOGRAM::HistogramSetting::MeasureFromVector(std::vector<double>* data, int nbin, double extendRange){
    double minValue = STATISTICS::MIN(data);
    double maxValue = STATISTICS::MAX(data);
    auto [minBound, maxBound] = ExtendRange(minValue, maxValue, extendRange);
    Set(minBound, maxBound, nbin);
};

void HISTOGRAM::HistogramSetting::MeasureFromVector(std::vector<double>* data, double hbin, int minNBin, double extendRange){
    double minValue = STATISTICS::MIN(data);
    double maxValue = STATISTICS::MAX(data);
    auto [minBound, maxBound] = ExtendRange(minValue, maxValue, extendRange, minNBin * hbin);

    int nbin = ceil((maxBound - minBound) / hbin);

    Set(minBound, maxBound, nbin);
};

void HISTOGRAM::HistogramSetting::Print(FILE* fout){
    fprintf(fout, "# HMin: %13.7f\n", m_HistHMin);
    fprintf(fout, "# HMax: %13.7f\n", m_HistHMax);
    fprintf(fout, "# HBin: %13.7f\n", m_HistHBin);
    fprintf(fout, "# NBin: %13d\n", m_HistNBin);
};
bool HISTOGRAM::operator== (const HistogramSetting& lhs, const HistogramSetting& rhs){
    return lhs.m_HistHMin == rhs.m_HistHMin &&
           lhs.m_HistHMax == rhs.m_HistHMax &&
           lhs.m_HistHBin == rhs.m_HistHBin &&
           lhs.m_HistNBin == rhs.m_HistNBin;
};
bool HISTOGRAM::operator!= (const HistogramSetting& lhs, const HistogramSetting& rhs){
    return !(lhs == rhs);
};

double HISTOGRAM::HISTOGRAM::GetTotalCount(){
    double value = 0;
    for(int i=0;i<count.size();i++){
        value += count[i];
    }
    return value;
};

void HISTOGRAM::HISTOGRAM::NormalizePercent(){
    double sum = 0.0;
    for(int i=0;i<percent.size();i++){
        sum += count[i];
    }
    if( sum != 0.0 ){
        for(int i=0;i<percent.size();i++){
            percent[i] = count[i] / sum;
        }
    }
};

void HISTOGRAM::HISTOGRAM::NormalizePercent(double sum){
    if( sum != 0.0 ){
        for(int i=0;i<percent.size();i++){
            percent[i] = count[i] / sum;
        }
    }
};

void HISTOGRAM::HISTOGRAM::Print(FILE* fout){
    fprintf(fout, "# HMin %23.16E\n", min);
    fprintf(fout, "# HMax %23.16E\n", max);
    fprintf(fout, "# HBin %23.16E\n", step);
    fprintf(fout, "# NBin %23d\n", size);

    for(int i=0;i<size;i++){
        fprintf(fout, "%23.16E %23.16E %23.16E %23.0f %23.16E\n",
            min + step * i,
            min + step * (i + 0.5),
            min + step * (i + 1.0),
            count[i],
            percent[i]);
    }
};

int HISTOGRAM::HISTOGRAM::GetNBin(){
    return size;
};

double HISTOGRAM::HISTOGRAM::GetHBin(){
    return step;
};

void HISTOGRAM::HISTOGRAM::SetNBin(int value){
    size = value;
};

void HISTOGRAM::HISTOGRAM::SetHBin(double value){
    step = value;
};

void HISTOGRAM::HISTOGRAM::Allocate(double histMin, double histMax, double histHbin){
    min = histMin;
    max = histMax;
    step = histHbin;
    size = round((max - min) / step);
    stepInverse = 1.0 / step;
    percent.resize(size);
    count.resize(size);
};

void HISTOGRAM::HISTOGRAM::Allocate(double histMin, double histMax, int histNbin){
    min = histMin;
    max = histMax;
    size = histNbin;
    step = (histMax - histMin) / size;
    stepInverse = 1.0 / step;
    percent.resize(size);
    count.resize(size);
};

void HISTOGRAM::HISTOGRAM::Allocate(HistogramSetting histogramSetting){
    Allocate(histogramSetting.m_HistHMin, histogramSetting.m_HistHMax, histogramSetting.m_HistNBin);
};

void HISTOGRAM::HISTOGRAM::CopyFrom(HISTOGRAM& rhs){
    this->min = rhs.min;
    this->max = rhs.max;
    this->step = rhs.step;
    this->size = rhs.size;
    this->stepInverse = rhs.stepInverse;
    this->count = rhs.count;
    this->percent = rhs.percent;
};

bool HISTOGRAM::HISTOGRAM::ResizeToFitNewBound(double value, size_t* zeroIndexInNewHistogram){
    // fix step and stepInverse
    if( value < min){
        int index = floor((value - min) * stepInverse);
        min = min + index * step;
        auto countTemp = count;
        auto percentTemp = percent;
        count.clear();
        percent.clear();
        count.resize(size - index);
        percent.resize(size - index);
        for(int i=0;i<size;i++){
            count[i - index] = countTemp[i];
            percent[i - index] = percentTemp[i];
        }
        size = size - index;
        if( zeroIndexInNewHistogram ) *zeroIndexInNewHistogram = -index;
        return true;
    }else if( value >= max ){
        int index = (value - min) * stepInverse;
        count.resize(index + 1);
        percent.resize(index + 1);
        size = index + 1;
        max = min + step * size;
        if( zeroIndexInNewHistogram ) *zeroIndexInNewHistogram = 0;
        return true;
    }else{
        if( zeroIndexInNewHistogram ) *zeroIndexInNewHistogram = 0;
        return false;
    }
};

bool HISTOGRAM::HISTOGRAM::ReBoundBySetting(HistogramSetting *setting, int ncore, double eps, bool bPrint){
    bool bReturn = fabs(setting->m_HistHBin - step) < eps;

    std::vector<double> countOld = count;
    int index0 = (min - setting->m_HistHMin) / step;

    Allocate(*setting);

    memset(&count[0], 0, sizeof(double) * setting->m_HistNBin);

    #pragma omp parallel for num_threads(ncore)
    for(int i=0;i<countOld.size();i++){
        count[i + index0] = countOld[i];
    }

    if( bPrint && !bReturn ){
        fprintf(stderr, "Warning: hbin changes more than eps new %23.16E old %23.16E\n", setting->m_HistHBin, step);
    }

    return bReturn;
};

    void VerifySelectedNumberOfFrame(int* selection, int nframe){
        if( selection[0] < 0 ) selection[0] = 0;
        if( selection[1] < 0 || selection[1] > nframe ) selection[1] = nframe;
        if( selection[2] <= 0 ) selection[2] = 1;
    }

void HISTOGRAM::HISTOGRAM::HistogramVector(std::vector<double>* data, int first, int last, int freq){
    int frame[3] = { first, last, freq };

    if( first < 0 ) first = 0;
    if( last < 0 || last > data->size() ) last = data->size();
    if( freq <= 0 ) freq = 1;

    for(int i=frame[0];i<frame[1];i+=frame[2]){
        double value = (*data)[i];
        int index = floor(double(value - min) * stepInverse);
        if( index >= 0 && index < size ){
            count[index]++;
        }
    }
}
bool HISTOGRAM::HISTOGRAM::AddValue(double value){
    int index = (value - min) * stepInverse;
    if( MOLUTILITY::InBound(index, 0, size) ){
        count[index] ++;
        return true;
    }else{
        return false;
    }
};

HISTOGRAM::HISTOGRAM& HISTOGRAM::HISTOGRAM::operator= (HISTOGRAM& rhs){
    this->min = rhs.min;
    this->max = rhs.max;
    this->step = rhs.step;
    this->size = rhs.size;
    this->stepInverse = rhs.stepInverse;
    this->count = rhs.count;
    this->percent = rhs.percent;
    return *this;
};

void HISTOGRAM::HistogramMultiple::NormalizePercent(int index, double sum){
    double sum0 = 0.0;
    for(int i=0;i<m_Data[index].size();i++){
        sum0 += m_Data[index][i];
    }
    if( sum0 != 0.0 ){
        double scaler = sum / sum0;
        for(int i=0;i<m_Data[index].size();i++){
            m_Data[index][i] *= scaler;
        }
    }
};
void HISTOGRAM::HistogramMultiple::NormalizePercent(double sum){
    for(int i=0;i<m_Data.size();i++){
        NormalizePercent(i, sum);
    }
};
double HISTOGRAM::HistogramMultiple::GetSum(int index){
    double sum0 = 0.0;
    for(int i=0;i<m_Data[index].size();i++){
        sum0 += m_Data[index][i];
    }
    return sum0;
};
void HISTOGRAM::HistogramMultiple::Print(FILE* fout){
    fprintf(fout, "# HMin %23.16E\n", m_HistogramSetting.m_HistHMin);
    fprintf(fout, "# HMax %23.16E\n", m_HistogramSetting.m_HistHMax);
    fprintf(fout, "# HBin %23.16E\n", m_HistogramSetting.m_HistHBin);
    fprintf(fout, "# NBin %23d\n", m_HistogramSetting.m_HistNBin);

    for(int i=0;i<m_HistogramSetting.m_HistNBin;i++){
        fprintf(fout, "%23.16E %23.16E %23.16E ",
            m_HistogramSetting.m_HistHMin + m_HistogramSetting.m_HistHBin * i,
            m_HistogramSetting.m_HistHMin + m_HistogramSetting.m_HistHBin * (i + 0.5),
            m_HistogramSetting.m_HistHMin + m_HistogramSetting.m_HistHBin * (i + 1.0));
        for(int j=0;j<m_Data.size();j++){
            fprintf(fout, "%23.16E ", m_Data[j][i]);
        }
        fprintf(fout, "\n");
    }
};
void HISTOGRAM::HistogramMultiple::Allocate(int numberOfData, HistogramSetting histogramSetting){
    m_HistogramSetting = histogramSetting;
    m_Data.resize(numberOfData);
    for(int i=0;i<numberOfData;i++){
        m_Data[i].resize(m_HistogramSetting.m_HistNBin);
    }
};
void HISTOGRAM::HistogramMultiple::Allocate(int numberOfData, double histMin, double histMax, double histHbin){
    m_HistogramSetting.Set(histMin, histMax, histHbin);
    m_Data.resize(numberOfData);
    for(int i=0;i<numberOfData;i++){
        m_Data[i].resize(m_HistogramSetting.m_HistNBin);
    }
};
void HISTOGRAM::HistogramMultiple::Allocate(int numberOfData, double histMin, double histMax, int histNbin){
    m_HistogramSetting.Set(histMin, histMax, histNbin);
    m_Data.resize(numberOfData);
    for(int i=0;i<numberOfData;i++){
        m_Data[i].resize(m_HistogramSetting.m_HistNBin);
    }
};
void HISTOGRAM::HistogramMultiple::InitializeValue(int index, double value){
    for(int i=0;i<m_Data[index].size();i++){
        m_Data[index][i] = value;
    }
};
void HISTOGRAM::HistogramMultiple::CopyFrom(HistogramMultiple& rhs){
    m_HistogramSetting = rhs.m_HistogramSetting;
    m_Data = rhs.m_Data;
};
bool HISTOGRAM::HistogramMultiple::SumOMP(std::vector<HistogramMultiple>& rhs, std::vector<int> operation){
    if( rhs.size() ){
        bool bSame = true;
        for(int i=0;i<rhs.size();i++){
            if( m_HistogramSetting != rhs[i].m_HistogramSetting ||
                m_Data.size() != rhs[i].m_Data.size() ){
                bSame = false;
                break;
            }
        }
        if( !bSame ){
            fprintf(stderr, "Error: OMP histograms are different\n");
            return false;
        }else{
            if( m_Data.size() != operation.size() ){
                fprintf(stderr, "Error: operation is not properly specified\n");
                return false;
            }else{
                for(int i=0;i<rhs.size();i++){
                    for(int j=0;j<m_Data.size();j++){
                        if( operation[j] == 0 ){
                            for(int k=0;k<m_Data[j].size();k++){
                                m_Data[j][k] += rhs[i].m_Data[j][k];
                            }
                        }else if( operation[j] == 1 ){
                            for(int k=0;k<m_Data[j].size();k++){
                                m_Data[j][k] = std::min<double>(m_Data[j][k], rhs[i].m_Data[j][k]);
                            }
                        }else if( operation[j] == 2 ){
                            for(int k=0;k<m_Data[j].size();k++){
                                m_Data[j][k] = std::max<double>(m_Data[j][k], rhs[i].m_Data[j][k]);
                            }
                        }else{
                            fprintf(stderr, "Error: undefined operation %d\n", operation[j]);
                            return false;
                        }
                    }
                }
                return true;
            }
        }
    }else{
        fprintf(stderr, "Error: empty OMP histogram parts\n");
        return false;
    }
};
HISTOGRAM::HistogramMultiple& HISTOGRAM::HistogramMultiple::operator= (HistogramMultiple& rhs){
    m_HistogramSetting = rhs.m_HistogramSetting;
    m_Data = rhs.m_Data;
    return *this;
};
double HISTOGRAM::HistogramMultiple::HistHBin(){
    return m_HistogramSetting.m_HistHBin;
};
double HISTOGRAM::HistogramMultiple::HistHBinInverse(){
    return m_HistogramSetting.m_HistHBinInverse;
};
double HISTOGRAM::HistogramMultiple::HistHMin(){
    return m_HistogramSetting.m_HistHMin;
};
double HISTOGRAM::HistogramMultiple::HistHMax(){
    return m_HistogramSetting.m_HistHMax;
};
int HISTOGRAM::HistogramMultiple::HistNBin(){
    return m_HistogramSetting.m_HistNBin;
};
bool HISTOGRAM::HistogramMultiple::AddValue(double valueBin, double value, int dataIndex){
    int index = (valueBin - m_HistogramSetting.m_HistHMin) * m_HistogramSetting.m_HistHBinInverse;
    if( MOLUTILITY::InBound(index, 0, m_HistogramSetting.m_HistNBin) ){
        m_Data[dataIndex][index] += value;
        return true;
    }else{
        return false;
    }
};
bool HISTOGRAM::HistogramMultiple::AddValue(int binIndex, double value, int dataIndex){
    if( MOLUTILITY::InBound(binIndex, 0, m_HistogramSetting.m_HistNBin) ){
        m_Data[dataIndex][binIndex] += value;
        return true;
    }else{
        return false;
    }
};
bool HISTOGRAM::HistogramMultiple::MaxValue(double valueBin, double value, int dataIndex){
    int index = (valueBin - m_HistogramSetting.m_HistHMin) * m_HistogramSetting.m_HistHBinInverse;
    if( MOLUTILITY::InBound(index, 0, m_HistogramSetting.m_HistNBin) ){
        m_Data[dataIndex][index] = std::max<double>(m_Data[dataIndex][index], value);
        return true;
    }else{
        return false;
    }
};
bool HISTOGRAM::HistogramMultiple::MaxValue(int binIndex, double value, int dataIndex){
    if( MOLUTILITY::InBound(binIndex, 0, m_HistogramSetting.m_HistNBin) ){
        m_Data[dataIndex][binIndex] = std::max<double>(m_Data[dataIndex][binIndex], value);
        return true;
    }else{
        return false;
    }
};
bool HISTOGRAM::HistogramMultiple::MinValue(double valueBin, double value, int dataIndex){
    int index = (valueBin - m_HistogramSetting.m_HistHMin) * m_HistogramSetting.m_HistHBinInverse;
    if( MOLUTILITY::InBound(index, 0, m_HistogramSetting.m_HistNBin) ){
        m_Data[dataIndex][index] = std::min<double>(m_Data[dataIndex][index], value);
        return true;
    }else{
        return false;
    }
};
bool HISTOGRAM::HistogramMultiple::MinValue(int binIndex, double value, int dataIndex){
    if( MOLUTILITY::InBound(binIndex, 0, m_HistogramSetting.m_HistNBin) ){
        m_Data[dataIndex][binIndex] = std::min<double>(m_Data[dataIndex][binIndex], value);
        return true;
    }else{
        return false;
    }
};
void HISTOGRAM::HistogramMultiple::AddValueNoCheck(int binIndex, double value, int dataIndex){
    m_Data[dataIndex][binIndex] += value;
};
int HISTOGRAM::HistogramMultiple::GetValueIndex(double value){
    int index = (value - m_HistogramSetting.m_HistHMin) * m_HistogramSetting.m_HistHBinInverse;
    return index;
}
bool HISTOGRAM::HistogramMultiple::InBound(double value){
    if( MOLUTILITY::InBound(value, m_HistogramSetting.m_HistHMin, m_HistogramSetting.m_HistHMax) ){
        return true;
    }else{
        return false;
    }
}
bool HISTOGRAM::HistogramMultiple::InBound(int value){
    if( MOLUTILITY::InBound(value, 0, m_HistogramSetting.m_HistNBin) ){
        return true;
    }else{
        return false;
    }
}
void HISTOGRAM::HistogramMultiple::HistogramVector(std::vector<double>* data, int dataIndex, int first, int last, int freq){
    int frame[3] = { first, last, freq };
    if( first < 0 ) first = 0;
    if( last < 0 || last > data->size() ) last = data->size();
    if( freq <= 0 ) freq = 1;
    for(int i=frame[0];i<frame[1];i+=frame[2]){
        double value = (*data)[i];
        int index = floor(double(value - HistHMin()) * HistHBinInverse());
        if( InBound(index) ){
            m_Data[dataIndex][index]++;
        }
    }
}
