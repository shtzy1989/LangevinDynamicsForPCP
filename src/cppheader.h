#ifndef MD_UTILITY_CPPHEADER_H
#define MD_UTILITY_CPPHEADER_H

#include <limits>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <cfloat>
#include <unordered_map>
#include <iomanip>
#include <functional>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <sys/time.h>
#include "unistd.h"
#include <time.h>
#include <omp.h>
#include <cstdarg>
#include <utility>
#include <complex>
#include <thread>
#include <execution>
#ifdef FILESYSTEM
#include <filesystem>
#endif
#include <cstring>
#include "array.h"
#include <bitset>
#include <queue>
#include <stack>
#include <initializer_list>
#include <random>
#include <regex>
#include <type_traits>
#include <chrono>
#include <unordered_set>
#if __cplusplus >= 202300L
#include <any>
#include <format>
#include <string_view>
// #include <print>
#endif

#ifndef __SIZEOF_INT128__
    typedef long long __int128;              // fallback, 64-bit only
    typedef unsigned long long unsigned___int128;
#endif

typedef std::string String;
#ifndef MD_ARRAY_H
template <typename T> using Array = std::vector<T>;
#endif
template <typename T> using Set = std::set<T>;
template <typename T> using LinkedArray = std::list<T>;
template <typename T1, typename T2> using Map = std::map<T1, T2>;
template <typename T> using Complex = std::complex<T>;

#define SAFE_DELETE(a) if( (a) != NULL ) delete (a); (a) = NULL;
#define SAFE_DELETE_ARRAY(a) if( (a) != NULL ) delete [] (a); (a) = NULL;
#define SAFE_RELEASE(p) { if ( (p) ) { (p)->Release(); (p) = 0; } }

template<typename T>
const char* GetPrintfFormat() {
    if constexpr (std::is_same_v<T, int>) {
        return "%d";
    } else if constexpr (std::is_same_v<T, unsigned int>) {
        return "%u";
    } else if constexpr (std::is_same_v<T, double>) {
        return "%f";
    } else if constexpr (std::is_same_v<T, float>) {
        return "%f";
    } else if constexpr (std::is_same_v<T, char>) {
        return "%c";
    } else if constexpr (std::is_same_v<T, char*>) {
        return "%s";
    } else if constexpr (std::is_same_v<T, const char*>) {
        return "%s";
    }
    // Add more types as needed

    // Fallback or error
    return "Unsupported type";
}

namespace VectorAuxiliary{
    template <typename T>
    std::vector<T> LinearSpace(T start, T end, int numberOfPoints, bool endPoint = true, T* step = nullptr) {
        if (numberOfPoints <= 0) {
            return {};
        }

        T step0;
        if( endPoint && numberOfPoints == 1 ){
            step0 = T(0);
        }else if( endPoint ){
            step0 = (end - start) / static_cast<T>(numberOfPoints - 1);
        }else{
            step0 = (end - start) / static_cast<T>(numberOfPoints);
        }

        std::vector<T> result;
        result.reserve(numberOfPoints);

        for(int i=0;i<numberOfPoints;i++){
            result.push_back(start + step0 * i);
        }

        if( step ) *step = step0;
        return result;
    }


    template <typename T>
    std::vector<T> LogSpace(T start, T end, int numberOfPoints, bool endPoint = true, T* step = nullptr) {
        if( start <= 0 || end <= 0 ){
            throw std::domain_error("LogSpace: start and end must be positive.");
        }
        if( numberOfPoints <= 0 ){
            return {};
        }

        T logStart = std::log(start);
        T logEnd   = std::log(end);
        T step0 = endPoint && numberOfPoints > 1
                ? (logEnd - logStart) / (numberOfPoints - 1)
                : (logEnd - logStart) / numberOfPoints;

        std::vector<T> result;
        result.reserve(numberOfPoints);

        for(int i=0;i<numberOfPoints;i++){
            result.push_back(std::exp(logStart + step0 * i));
        }

        if( step ) *step = step0;
        return result;
    }

    // math functions
    template <typename T, typename UnaryFunc>
    std::vector<T> ApplyFunction(const std::vector<T>& data, UnaryFunc func) {
        std::vector<T> result(data.size());
        std::transform(data.begin(), data.end(), result.begin(), func);
        return result;
    }

    template <typename T, typename UnaryFunc>
    void ApplyFunctionInPlace(std::vector<T>& data, UnaryFunc func) {
        std::transform(data.begin(), data.end(), data.begin(), func);
    }

    // log 
    template <typename T>
    std::vector<T> Log(const std::vector<T>& data) {
        std::vector<T> result(data.size());
        auto log_func = static_cast<T(*)(T)>(std::log);  // Disambiguate
        std::transform(data.begin(), data.end(), result.begin(), log_func);
        return result;
    }

    template <typename T>
    void LogInPlace(std::vector<T>& data) {
        auto log_func = static_cast<T(*)(T)>(std::log);  // Disambiguate
        std::transform(data.begin(), data.end(), data.begin(), log_func);
    }

    // log10 
    template <typename T>
    std::vector<T> Log10(const std::vector<T>& data) {
        std::vector<T> result(data.size());
        auto log_func = static_cast<T(*)(T)>(std::log10);  // Disambiguate
        std::transform(data.begin(), data.end(), result.begin(), log_func);
        return result;
    }

    template <typename T>
    void Log10InPlace(std::vector<T>& data) {
        auto log_func = static_cast<T(*)(T)>(std::log10);  // Disambiguate
        std::transform(data.begin(), data.end(), data.begin(), log_func);
    }

    // pow
    template <typename T>
    std::vector<T> Pow(const std::vector<T>& data, T exponent) {
        std::vector<T> result(data.size());
        std::transform(data.begin(), data.end(), result.begin(),
                    [exponent](T x) { return std::pow(x, exponent); });
        return result;
    }

    template <typename T>
    void PowInPlace(std::vector<T>& data, T exponent) {
        std::transform(data.begin(), data.end(), data.begin(),
                    [exponent](T x) { return std::pow(x, exponent); });
    }

    // other
    template <typename T>
    void VectorReindex(std::vector<T>& data, const std::vector<int>& newIndices) {
        std::vector<T> temp(data.size());

        for (size_t i = 0; i < newIndices.size(); ++i) {
            temp[i] = std::move(data[newIndices[i]]);
        }

        data = std::move(temp);
    }

    template<typename T>
    T CumulativeProduct(std::vector<T> * lhs, int first = 0, int last = -1, int freq = 1){
        if( last == -1 ){
            last = lhs->size();
        }
        T product = 1;
        for(int i=first;i<last;i+=freq){
            product *= (*lhs)[i];
        }
        return product;
    };

    template<typename T>
    std::vector<T> Slice(std::vector<T>& lhs, int start, int stop, int stride){
        int len = lhs.size();
        
        // Normalize indices beyond the boundaries
        if( stride > 0 ){
            if( start < 0 ) start = len + start;
            if( stop < 0 ) stop = len + stop;
            start = std::max(0, start);
            stop = std::min(len, stop);
        }else if( stride == 0 ){
            fprintf(stderr, "Error: stride cannot be 0\n");
            return {};
        }else{
            if( start >= len ) start = len - 1;
            stop = stop < -1 ? -1 : stop; // Allowing stop to be -1 for negative stride
        }

        std::vector<T> result;
        if (stride > 0) {
            for (int i = start; i < stop; i += stride) {
                result.push_back(lhs[i]);
            }
        } else { // Handle negative stride
            for (int i = start; i > stop; i += stride) {
                result.push_back(lhs[i]);
            }
        }
        return result;
    }

    template<typename T>
    std::vector<T> AdjointVector(std::vector<T>* lhs, std::vector<T>* rhs){
        std::vector<T> result;
        result.reserve(lhs->size() + rhs->size());
        for(int i=0;i<lhs->size();i++){
            result.push_back((*lhs)[i]);
        }
        for(int i=0;i<rhs->size();i++){
            result.push_back((*rhs)[i]);
        }
        return result;
    };

    template<typename T>
    void AdjointVectorToLHS(std::vector<T>* lhs, std::vector<T>* rhs){
        lhs->reserve(lhs->size() + rhs->size());
        for(int i=0;i<rhs->size();i++){
            lhs->push_back((*rhs)[i]);
        }
    };

    template<typename T>
    size_t SumVectorGroupSize(std::vector<std::vector<T> >* vectorGroup){
        size_t value = 0;
        for(int i=0;i<vectorGroup->size();i++){
            value += (*vectorGroup)[i].size();
        }
        return value;
    }

    template<typename T>
    size_t MaxVectorGroupSize(std::vector<std::vector<T> >* vectorGroup){
        size_t value = -ULONG_MAX;
        for(int i=0;i<vectorGroup->size();i++){
            value = std::max<size_t>(value, (*vectorGroup)[i].size());
        }
        return value;
    }

    template<typename T>
    size_t MinVectorGroupSize(std::vector<std::vector<T> >* vectorGroup){
        size_t value = ULONG_MAX;
        for(int i=0;i<vectorGroup->size();i++){
            value = std::min<size_t>(value, (*vectorGroup)[i].size());
        }
        return value;
    }

    template<typename T>
    void NormalizedVector(std::vector<T> * rhs, T value = 1.0){
        if( rhs->size() ){
            T sum = 0.0;
            for(int i=0;i<rhs->size();i++){
                sum += (*rhs)[i] * (*rhs)[i];
            }
            T scaler = sqrt(1.0 / sum) * value;
            for(int i=0;i<rhs->size();i++){
                (*rhs)[i] *= scaler;
            }
        }
    }
    template<typename T>
    void ScaleVectorSum(std::vector<T> * rhs, T value = 1.0){
        if( rhs->size() ){
            T sum = 0.0;
            for(int i=0;i<rhs->size();i++){
                sum += (*rhs)[i];
            }
            T scaler = value / sum;
            for(int i=0;i<rhs->size();i++){
                (*rhs)[i] *= scaler;
            }
        }
    }
    template<typename T>
    size_t MinValueInVector(std::vector<T> * rhs, size_t lowBound = 0, size_t highBound = -1){
        if( rhs->size() ){
            if( highBound < 0 || highBound > rhs->size() ) highBound = rhs->size();
            if( highBound <= lowBound ){
                return -1;
            }
            size_t index = lowBound;
            T minvalue = (*rhs)[lowBound];
            for(size_t i=lowBound+1;i<highBound;i++){
                if( (*rhs)[i] < minvalue ){
                    index = i;
                    minvalue = (*rhs)[i];
                }
            }
            return index;
        }else{
            return -1;
        }
    }
    template<typename T>
    size_t MaxValueInVector(std::vector<T> * rhs, size_t lowBound = 0, size_t highBound = -1){
        if( rhs->size() ){
            if( highBound < 0 || highBound > rhs->size() ) highBound = rhs->size();
            if( highBound <= lowBound ){
                return -1;
            }
            size_t index = lowBound;
            T minvalue = (*rhs)[lowBound];
            for(size_t i=lowBound+1;i<highBound;i++){
                if( (*rhs)[i] > minvalue ){
                    index = i;
                    minvalue = (*rhs)[i];
                }
            }
            return index;
        }else{
            return -1;
        }
    }

    template<typename T>
    std::vector<T> GenerateRange(T start, T end, T freq){
        std::vector<T> result;
        for(T i=start;i<end;i+=freq){
            result.push_back(i);
        }
        return result;
    }

    template<typename T>
    std::vector<T> GenerateRange(int end){
        std::vector<T> result;
        for(int i=0;i<end;i++){
            result.push_back(i);
        }
        return result;
    }

    template<typename T>
    std::vector<T> MergeVectorOfVector(std::vector<std::vector<T> >*rhs){
        std::vector<T> result;
        for(int i=0;i<rhs->size();i++){
            result.insert(result.end(), (*rhs)[i].begin(), (*rhs)[i].end());
        }
        return result;
    }

    template<typename T>
    void SumOMPVectorOfVector(std::vector<T>* result, std::vector<std::vector<T> >*rhs){
        for(int i=0;i<rhs->size();i++){
            for(int j=0;j<result->size();j++){
                (*result)[j] += (*rhs)[i][j];
            }
        }
    }

    template<typename T>
    void ZeroVectorOfVector(std::vector<std::vector<T> > *rhs){
        for(int i=0;i<rhs->size();i++){
            for(int j=0;j<(*rhs)[i].size();j++){
                (*rhs)[i][j] = 0;
            }
        }
    }

    template<typename T>
    void PrintVector(std::vector<T>* rhs, const char* format, FILE* fout = stderr){
        for(int i=0;i<rhs->size();i++){
            fprintf(fout, format, (*rhs)[i]);
        }
    }

    template<typename T>
    std::string PrintVectorToString(std::vector<T>* rhs, const char* format){
        char buffer[256];
        std::string out = "";
        for(int i=0;i<rhs->size();i++){
            sprintf(buffer, format, (*rhs)[i]);
            out = out + buffer;
        }
        return out;
    }

    template<typename T>
    void PrintVectorFormatted(std::vector<T>* rhs, const char* format, int nperline, FILE* fout = stderr){
        for(int i=0;i<rhs->size();i++){
            fprintf(fout, format, (*rhs)[i]);
            if( (i + 1) % nperline == 0 || i == rhs->size() -1 ){
                fprintf(fout, "\n");
            }
        }
    }

    template<typename T>
    std::string Vector2String(std::vector<T> rhs, std::string delimiter = " ", std::string to = " to "){
        std::string result = "";
        if( rhs.size() ){
            std::sort(rhs.begin(), rhs.end());
            auto last = std::unique(rhs.begin(), rhs.end());
            rhs.erase(last, rhs.end());
            T initialValue = rhs[0];
            T endValue = rhs[0];
            for(int i=1;i<rhs.size();i++){
                if( rhs[i] != rhs[i-1] + 1 ){
                    if( initialValue == endValue ){
                        result += std::to_string(initialValue) + delimiter;
                    }else{
                        result += std::to_string(initialValue) + to + std::to_string(endValue) + delimiter;
                    }
                    initialValue = rhs[i];
                    endValue = rhs[i];
                }else{
                    endValue = rhs[i];
                }
                if( i == rhs.size() - 1 ){
                    if( initialValue == endValue ){
                        result += std::to_string(initialValue) + delimiter;
                    }else{
                        result += std::to_string(initialValue) + to + std::to_string(endValue) + delimiter;
                    }
                }
            }
        }
        if( result.size() >= delimiter.size() && result.substr(result.size() - delimiter.size(), delimiter.size()) == delimiter ){
            return result.substr(0, result.size() - delimiter.size());
        }else{
            return result;
        }
    }
    template<typename T>
    void AllocateVector2D(std::vector<std::vector<T> >* rhs, size_t d1, size_t d2){
        rhs->resize(d1);
        for(int i=0;i<d1;i++){
            (*rhs)[i].resize(d2);
        }
    }
    template<typename T>
    void AllocateVector3D(std::vector<std::vector<std::vector<T> > >* rhs, size_t d1, size_t d2, size_t d3){
        rhs->resize(d1);
        for(int i=0;i<d1;i++){
            (*rhs)[i].resize(d2);
            for(int j=0;j<d2;j++){
                (*rhs)[i][j].resize(d3);
            }
        }
    }
    template<typename T>
    void AllocateVector4D(std::vector<std::vector<std::vector<std::vector<T> > > >* rhs, size_t d1, size_t d2, size_t d3, size_t d4){
        rhs->resize(d1);
        for(int i=0;i<d1;i++){
            (*rhs)[i].resize(d2);
            for(int j=0;j<d2;j++){
                (*rhs)[i][j].resize(d3);
                for(int k=0;k<d3;k++){
                    (*rhs)[i][j][k].resize(d4);
                }
            }
        }
    }
    template<typename T>
    void AllocateVector5D(std::vector<std::vector<std::vector<std::vector<std::vector<T> > > > >* rhs, size_t d1, size_t d2, size_t d3, size_t d4, size_t d5){
        rhs->resize(d1);
        for(int i=0;i<d1;i++){
            (*rhs)[i].resize(d2);
            for(int j=0;j<d2;j++){
                (*rhs)[i][j].resize(d3);
                for(int k=0;k<d3;k++){
                    (*rhs)[i][j][k].resize(d4);
                    for(int l=0;l<d4;l++){
                        (*rhs)[i][j][k][l].resize(d5);
                    }
                }
            }
        }
    }

    template<typename T>
    void AllocateVector6D(std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<T> > > > > >* rhs, size_t d1, size_t d2, size_t d3, size_t d4, size_t d5, size_t d6){
        rhs->resize(d1);
        for(int i=0;i<d1;i++){
            (*rhs)[i].resize(d2);
            for(int j=0;j<d2;j++){
                (*rhs)[i][j].resize(d3);
                for(int k=0;k<d3;k++){
                    (*rhs)[i][j][k].resize(d4);
                    for(int l=0;l<d4;l++){
                        (*rhs)[i][j][k][l].resize(d5);
                        for(int m=0;m<d5;m++){
                            (*rhs)[i][j][k][l][m].resize(d6);
                        }
                    }
                }
            }
        }
    }

    template<typename T>
    void AllocatePoint2D(T**& pointer, size_t d1, size_t d2){
        pointer = new T*[d1];
        for(int i=0;i<d1;i++){
            pointer[i] = new T[d2];
        }
    }
    template<typename T>
    void DeallocatePoint2D(T**& pointer, size_t d1, size_t d2){
        for(int i=0;i<d1;i++){
            SAFE_DELETE_ARRAY(pointer[i]);
        }
        SAFE_DELETE_ARRAY(pointer);
    }
    template<typename T>
    void AllocatePoint3D(T***& pointer, size_t d1, size_t d2, size_t d3){
        pointer = new T**[d1];
        for(int i=0;i<d1;i++){
            pointer[i] = new T*[d2];
            for(int j=0;j<d2;j++){
                pointer[i][j] = new T[d3];
            }
        }
    }
    template<typename T>
    void DeallocatePoint3D(T**& pointer, size_t d1, size_t d2, size_t d3){
        for(int i=0;i<d1;i++){
            for(int j=0;j<d2;j++){
                SAFE_DELETE_ARRAY(pointer[i][j]);
            }
            SAFE_DELETE_ARRAY(pointer[i]);
        }
        SAFE_DELETE_ARRAY(pointer);
    }
    template<typename T>
    int ValuePositionInSortedVector(T value, std::vector<T>* array){
        if( value >= array->back() ){
            return array->size() - 1;
        }else if( value < array->front() ){
            return -1;
        }else{
            int left = 0;
            int right = array->size();
            int mid = (left + right) / 2;
            while( right - left != 1 ){
                if( value >= (*array)[mid] ){
                    left = mid;
                }else{
                    right = mid;
                }
                mid = (left + right) / 2;
            }
            return mid;
        }
    };
    template<typename T>
    int ValuePositionInSortedVectorExist(T value, std::vector<T>* array){
        int index = ValuePositionInSortedVector(value, array);
        if( index >= 0 && index < array->size() && (*array)[index] == value ){
            return index;
        }else{
            return -1;
        }
    };
    template<typename T>
    int ValuePositionInSortedVectorDouble(T value, std::vector<T>* array, double tolerance = 0.00001){
        if( value >= array->back() ){
            return array->size() - 1;
        }else if( value < array->front() ){
            return -1;
        }else{
            int left = 0;
            int right = array->size();
            int mid = (left + right) / 2;
            while( right - left != 1 ){
                if( value > (*array)[mid] || fabs(value - (*array)[mid]) < tolerance ){
                    left = mid;
                }else{
                    right = mid;
                }
                mid = (left + right) / 2;
            }
            return mid;
        }
    };
    template<typename T>
    int ValuePositionInSortedVectorDoubleExist(T value, std::vector<T>* array, double tolerance = 0.00001){
        int index = ValuePositionInSortedVectorDouble(value, array);
        if( index >= 0 && index < array->size() && fabs((*array)[index] - value) < tolerance ){
            return index;
        }else{
            return -1;
        }
    };
    template<typename T>
    int ValuePositionInSortedVectorDoubleRigorous(T value, std::vector<T>* array){
        if( value >= array->back() ){
            return array->size() - 1;
        }else if( value < array->front() ){
            return -1;
        }else{
            int left = 0;
            int right = array->size();
            int mid = (left + right) / 2;
            while( right - left != 1 ){
                if( value > (*array)[mid] ){
                    left = mid;
                }else{
                    right = mid;
                }
                mid = (left + right) / 2;
            }
            return mid;
        }
    };
    template<typename T>
    std::vector<std::vector<T>* > GenerateVectorOfVectorOfPointer(std::vector<std::vector<T> >* input){
        std::vector<std::vector<T>* > result(input->size());
        for(int i=0;i<input->size();i++){
            result[i] = &(*input)[i];
        }
        return result;
    }
};

static inline void memset_zero(char *buf, size_t size)
{
    memset(buf, 0, size);
}

static inline void memset_zero_omp(char *buf, size_t size)
{
    size_t my_start, my_size;
 
    if (omp_in_parallel())
    {
        int id = omp_get_thread_num();
        int num = omp_get_num_threads();
 
        my_start = (id*size)/num;
        my_size = ((id+1)*size)/num - my_start;
    }
    else
    {
        my_start = 0;
        my_size = size;
    }
 
    memset(buf + my_start, 0, my_size);
}

static inline void memset_one_omp(char *buf, size_t size)
{
    size_t my_start, my_size;

    if (omp_in_parallel())
    {
        int id = omp_get_thread_num();
        int num = omp_get_num_threads();

        my_start = (id*size)/num;
        my_size = ((id+1)*size)/num - my_start;
    }
    else
    {
        my_start = 0;
        my_size = size;
    }

    memset(buf + my_start, 255, my_size);
}

static inline void memcpy_omp(char *dest, char *source, size_t size)
{
    size_t my_start, my_size;

    if (omp_in_parallel())
    {
        int id = omp_get_thread_num();
        int num = omp_get_num_threads();

        my_start = (id*size)/num;
        my_size = ((id+1)*size)/num - my_start;
    }
    else
    {
        my_start = 0;
        my_size = size;
    }

    memcpy(dest + my_start, source + my_start, my_size);
}

template <typename T>
class Array2D {
protected:
    std::vector<T> m_Data;
    size_t m_Size;
    size_t m_D1;
    size_t m_D2;
public:
    Array2D operator=(const Array2D& rhs){
        m_Data = rhs.m_Data;
        m_Size = rhs.m_Size;
        m_D1 = rhs.m_D1;
        m_D2 = rhs.m_D2;
        return *this;
    }
    bool Allocate(int d1, int d2) {
        m_D1 = d1;
        m_D2 = d2;
        m_Size = d1 * d2;

        try{
            m_Data.resize(m_Size);
        }catch(...){
            return false;
        }

        return true;
    }
    int D1() {
        return m_D1;
    }
    int D2() {
        return m_D2;
    }
    int Size() {
        return m_Size;
    }
    void Set(int i, int j, T value) {
        m_Data[GetIndex(i, j)] = value;
    }
    void Set(int i, T value) {
        m_Data[i] = value;
    }
    T& Get(int i, int j) {
        return m_Data[GetIndex(i, j)];
    }
    T& Get(int i) {
        return m_Data[i];
    }                                                                                        
    size_t GetIndex(int i, int j) {
        return i * m_D2 + j;
    }
    void ZeroAll(){
        memset(&m_Data[0], 0, sizeof(T) * m_Size);
    }
    void SetAll(T value){
        for(int i=0;i<m_Size;i++){
            m_Data[i] = value;
        }
    }
    T& operator [] (int index){
        return m_Data[index];
    }
    bool IndexInRange(int i, int j){
        return 0 <= i && i < m_D1 && 0 <= j && j < m_D2;
    };
    std::vector<T>* GetAddress(){
        return &m_Data;
    };
};

template <typename T>
class Array3D {
protected:
    std::vector<T> m_Data;
    size_t m_Size;
    int m_D1;
    int m_D2;
    int m_D3;
    size_t m_D2D3;
public:
    Array3D operator=(const Array3D& rhs){
        m_Data = rhs.m_Data;
        m_Size = rhs.m_Size;
        m_D1 = rhs.m_D1;
        m_D2 = rhs.m_D2;
        m_D3 = rhs.m_D3;
        m_D2D3 = rhs.m_D2D3;
        return *this;
    }
    bool Allocate(int d1, int d2, int d3) {
        m_D1 = d1;
        m_D2 = d2;
        m_D3 = d3;
        m_D2D3 = (size_t)m_D2 * (size_t)m_D3;
        m_Size = (size_t)d1 * (size_t)d2 * (size_t)d3;
        try{
            m_Data.resize(m_Size);
        }catch(...){
            return false;
        }
        return true;
    }
    int D1() {
        return m_D1;
    }
    int D2() {
        return m_D2;
    }
    int D3() {
        return m_D3;
    }
    int Dimension(int index){
        if( index == 0 ){
            return m_D1;
        }else if( index == 1 ){
            return m_D2;
        }else{
            return m_D3;
        }
    }
    size_t Size() {
        return m_Size;
    }
    void Set(int i, int j, int k, T value) {
        m_Data[GetIndex(i, j, k)] = value;
    }
    void Set(size_t i, T value) {
        m_Data[i] = value;
    }
    void Add(int i, int j, int k, T value) {
        m_Data[GetIndex(i, j, k)] += value;
    }
    void Add(size_t i, T value) {
        m_Data[i] += value;
    }
    T& Get(int i, int j, int k) {
        return m_Data[GetIndex(i, j, k)];
    }
    T& Get(size_t i) {
        return m_Data[i];
    }
    size_t GetIndex(int i, int j, int k) {
        return (size_t)i * (size_t)m_D2D3 + (size_t)j * (size_t)m_D3 + (size_t)k;
    }
    size_t GetIndexWrap(int i, int j, int k) {
        i = (i % m_D1 + m_D1) % m_D1;
        j = (j % m_D1 + m_D1) % m_D1;
        k = (k % m_D1 + m_D1) % m_D1;
        return (size_t)i * m_D2D3 + j * m_D3 + k;
    }
    bool InBound(int i, int j, int k){
        return 0 <= i && i < m_D1 &&
               0 <= j && j < m_D2 &&
               0 <= k && k < m_D3;
    }
    std::tuple<int, int, int> GetIndex(size_t index){
        int x = int(index / m_D2D3);
        size_t residue = index % m_D2D3;
        int y = int(residue / (size_t)m_D3);
        int z = int(residue % (size_t)m_D3);
        return std::tuple<int, int, int>(x, y, z);
    }
    void ZeroAll(){
        memset(&m_Data[0], 0, sizeof(T) * m_Size);
    }
    void SetAll(T value){
        for(size_t i=0;i<m_Size;i++){
            m_Data[i] = value;
        }
    }
    T& operator [] (size_t index){
        return m_Data[index];
    }
    std::vector<T>* GetAddress(){
        return &m_Data;
    };
};

/*
template<typename T>
class Pointer{
protected:
    T* m_Data;
    size_t m_Size;
public:
    Pointer() : m_Data(static_cast<T*>(0)){};
    Pointer(T* value) : m_Data(value){};
    operator T*&(){ return m_Data; };
    operator T*() const{ return m_Data; };
public:
    void NewArray(size_t size){
        m_Data = new T[size];
        m_Size = size;
    };
    void New(){
        m_Data = new T;
        m_Size = 1;
    };
    void Release(){
        SAFE_DELETE(m_Data);
    };
    void ReleaseArray(){
        SAFE_DELETE_ARRAY(m_Data); 
    }
    T& operator [](int index){
        return m_Data[index];
    }
    const T operator [](int index) const{
        return m_Data[index];
    }
    Pointer operator=(const Pointer& rhs){
        m_Data = rhs.m_Data;
        m_Size = rhs.m_Size;
    }
};
*/

#endif
