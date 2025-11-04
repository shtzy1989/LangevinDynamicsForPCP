#ifndef MD_RANDOMNUMBER_H
#define MD_RANDOMNUMBER_H

#include "cppheader.h"
#include "simd.h"
#include "mkl.h"

namespace RANDOMNUMBER{
    int SEED(int seed = -1);
    std::vector<int> RANDOM_NON_REPEAT_INT(int size);
    void RANDOM_VECTOR_NON_REPEAT(std::vector<unsigned int> *vector, int numberToChoose, int availableNumber);

    class UniformRandomNumberSequence{
    protected:
        size_t m_Process;
        long long m_TimeOfRenew;
        VSLStreamStatePtr m_Stream;
        std::vector<double> m_Bank;
        double m_LowBound;
        double m_HighBound;
        int m_BankSize;
    public:
        ~UniformRandomNumberSequence();
    public:
        int Initialize(size_t bankSize, int seed, double lowBound = 0.0, double highBound = 1.0);
        double GetRandomNumber();
        void Renew();
        void ResetRenewTime();
        long long GetRenewTime();
        size_t GetProcess();
        void SetRenewTime(long long value);
        void SetProcess(size_t value);
        void Advance(long long count, size_t progress);
        void AdvanceFromZero(long long count, size_t progress);
        int GetCoreBankSize();
        void PrintStatus();
        double GetLowBound();
        double GetHighBound();
    };

    class GaussianRandomNumberSequenceSIMD{
    protected:
        VSLStreamStatePtr m_Stream;
        alignas(SIMDBYTESIZE) std::vector<double> m_Bank;
        int m_BankSize;
        size_t m_Process;
    protected:
        long long m_TimeOfRenew;
        double m_Average;
        double m_Sigma;
    public:
        GaussianRandomNumberSequenceSIMD();
        ~GaussianRandomNumberSequenceSIMD();
    public:
        void Initialize(size_t bankSize, int seed, double average = 0.0, double sigma = 1.0);
        void Renew();
        void ResetRenewTime();
        long long GetRenewTime();
        void SetRenewTime(long long value);
        void SetProcess(size_t value);
        void Advance(long long count, size_t progress);
        void AdvanceFromZero(long long count, size_t progress);
        __mSIMDd GetRandomNumber();
        __mSIMDd GetRandomNumberSafe();
        size_t GetProcess();
        int GetCoreBankSize();
        void PrintStatus();
    };

    class UniformRandomNumberSequenceSIMD{
    protected:
        VSLStreamStatePtr m_Stream;
        alignas(SIMDBYTESIZE) std::vector<double> m_Bank;
        int m_BankSize;
        size_t m_Process;
    protected:
        long long m_TimeOfRenew;
        double m_LowBound;
        double m_HighBound;
    public:
        UniformRandomNumberSequenceSIMD();
        ~UniformRandomNumberSequenceSIMD();
    public:
        void Initialize(size_t bankSize, int seed, double lowBound = 0.0, double highBound = 1.0);
        void Renew();
        void ResetRenewTime();
        long long GetRenewTime();
        void SetRenewTime(long long value);
        void SetProcess(size_t value);
        void Advance(long long count, size_t progress);
        void AdvanceFromZero(long long count, size_t progress);
        __mSIMDd GetRandomNumber();
        __mSIMDd GetRandomNumberSafe();
        size_t GetProcess();
        int GetCoreBankSize();
        void PrintStatus();
    };

};

#endif
