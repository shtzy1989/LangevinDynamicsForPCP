#include "RANDOMNUMBER.h"
#include "simd.h"

namespace RANDOMNUMBER{
    int SEED(int seed){
        int value;
        if( seed == - 1 ){
            value = time(0);
        }else{
            value = seed;
        }
        srand(value);
        return value;
    };

    std::vector<int> RANDOM_NON_REPEAT_INT(int size){
        std::vector<int> result;
        result.reserve(size);
        for(int i=0;i<size;i++){
            bool bUnique = false;
            int newValue;
            while( !bUnique ){
                newValue = rand() % RAND_MAX;
                bUnique = true;
                for(int j=0;j<result.size();j++){
                    if( result[j] == newValue ){
                        bUnique = false;
                    } 
                }
            }
            result.push_back(newValue);
        }
        return result;
    }
    void RANDOM_VECTOR_NON_REPEAT(std::vector<unsigned int> *vector, int n, int m){
        // n = size of vector
        // m = max value in vector
        if( n != 0 ){
            std::vector<unsigned int> value(m);
            unsigned int r, t;

            if( n > m ){
                for(int i=0;i<m;i++){
                    vector->push_back(i);
                }
            }else{
                for(int i=0;i<m;i++){
                    value[i] = i;
                }
                vector->clear();
                r = (rand() % RAND_MAX) % m;
                vector->push_back(value[r]);
                t = value[r];
                value[r] = value[m-1];
                value[m-1] = t;
                m--;
                for(int i=1;i<n;i++){
                    r = (r + (rand() % RAND_MAX)) % m;
                    vector->push_back(value[r]);
                    t = value[r];
                    value[r] = value[m-1];
                    value[m-1] = t;
                    m--;
                }
            }
        }
    };

    // UniformRandomNumberSequence
    UniformRandomNumberSequence::~UniformRandomNumberSequence(){
        vslDeleteStream(&m_Stream);
        m_Bank.clear();
        m_Bank.shrink_to_fit();
    }
    int UniformRandomNumberSequence::Initialize(size_t bankSize, int seed, double lowBound, double highBound){
        m_LowBound = lowBound;
        m_HighBound = highBound;
        m_TimeOfRenew = 0;
        m_BankSize = bankSize;
        m_Bank.resize(m_BankSize);
        if( seed == -1 ){
            seed = std::random_device{}();
        }
        vslNewStream(&m_Stream, VSL_BRNG_MT19937, seed);
        m_Process = m_Bank.size();

        Renew();
        return seed;
    }
    void UniformRandomNumberSequence::Renew(){
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, m_Stream, m_BankSize, &m_Bank[0], m_LowBound, m_HighBound);
        m_Process = 0;
        m_TimeOfRenew++;
    };
    void UniformRandomNumberSequence::ResetRenewTime(){
        m_TimeOfRenew = 0;
    };
    long long UniformRandomNumberSequence::GetRenewTime(){
        return m_TimeOfRenew;
    };
    void UniformRandomNumberSequence::SetRenewTime(long long value){
        m_TimeOfRenew = value;
    };
    void UniformRandomNumberSequence::SetProcess(size_t value){
        m_Process = value;
    }
    void UniformRandomNumberSequence::Advance(long long count, size_t progress){
        for(long long i=0;i<count-1;i++){
            vslSkipAheadStream(m_Stream, (long long)(m_BankSize));
        }
        vslSkipAheadStream(m_Stream, (long long)(progress));
    };
    void UniformRandomNumberSequence::AdvanceFromZero(long long count, size_t progress){
        if( count == 0 ) count = 1;
        if( count == 1 ){
            m_Process = progress;
        }else{
            // when doing nothing, count is 1 (renew'ed in Initialization)
            // so subtract 2, one is the initialization renew, the other is the following renew
            for(long long i=0;i<count-2;i++){
                vslSkipAheadStream(m_Stream, (long long)(m_BankSize));
            }
            // subtract 1, which is the following renew
            m_TimeOfRenew = count-1;
            Renew();
            m_Process = progress;
        }
    };
    double UniformRandomNumberSequence::GetRandomNumber(){
        if( m_Process == m_Bank.size() ){
            Renew();
        }
        double value = m_Bank[m_Process];
        m_Process++;
        return value;
    }
    size_t UniformRandomNumberSequence::GetProcess(){
        return m_Process;
    };
    int UniformRandomNumberSequence::GetCoreBankSize(){
        return m_BankSize;
    };
    void UniformRandomNumberSequence::PrintStatus(){
        fprintf(stderr, "Time:    %10lld\n", m_TimeOfRenew);
        fprintf(stderr, "Process: %10zd\n", m_Process);
        fprintf(stderr, "First    %13.7f\n", m_Bank[0]);
    }
    double UniformRandomNumberSequence::GetLowBound(){
        return m_LowBound;
    };
    double UniformRandomNumberSequence::GetHighBound(){
        return m_HighBound;
    };

    // =======================================================================

    // =======================================================================
    GaussianRandomNumberSequenceSIMD::GaussianRandomNumberSequenceSIMD(){
    };
    GaussianRandomNumberSequenceSIMD::~GaussianRandomNumberSequenceSIMD(){
        vslDeleteStream(&m_Stream);
        m_Bank.clear();
        m_Bank.shrink_to_fit();
    };
    void GaussianRandomNumberSequenceSIMD::Initialize(size_t bankSize, int seed, double average, double sigma){
        m_Average = average;
        m_Sigma = sigma;
        m_TimeOfRenew = 0;
        m_BankSize = bankSize;
        m_Bank.resize(m_BankSize);
        // this time random number is taken care of in the runnner
        if( seed == -1 ){
            seed = std::random_device{}();
        }
        vslNewStream(&m_Stream, VSL_BRNG_MT19937, seed);
        m_Process = 0;

        Renew();
        // return seed;
    }
    void GaussianRandomNumberSequenceSIMD::Renew(){
        vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, m_Stream, m_BankSize, &m_Bank[0], m_Average, m_Sigma);
        m_Process = 0;
        m_TimeOfRenew++;
    }
    void GaussianRandomNumberSequenceSIMD::ResetRenewTime(){
        m_TimeOfRenew = 0;
    }
    long long GaussianRandomNumberSequenceSIMD::GetRenewTime(){
        return m_TimeOfRenew;
    }
    void GaussianRandomNumberSequenceSIMD::SetRenewTime(long long value){
        m_TimeOfRenew = value;
    }
    void GaussianRandomNumberSequenceSIMD::SetProcess(size_t value){
        m_Process = value;
    }
    void GaussianRandomNumberSequenceSIMD::Advance(long long count, size_t progress){
        for(long long i=0;i<count;i++){
            vslSkipAheadStream(m_Stream, (long long)(m_BankSize));
        }
        vslSkipAheadStream(m_Stream, (long long)(progress));
    }
    void GaussianRandomNumberSequenceSIMD::AdvanceFromZero(long long count, size_t progress){
        if( count == 0 ) count = 1;
        if( count == 1 ){
            m_Process = progress;
        }else{
            // when doing nothing, count is 1 (renew'ed in Initialization)
            // so subtract 2, one is the initialization renew, the other is the following renew
            for(long long i=0;i<count-2;i++){
                vslSkipAheadStream(m_Stream, (long long)(m_BankSize));
            }
            // subtract 1, which is the following renew
            m_TimeOfRenew = count - 1;
            Renew();
            m_Process = progress;
        }
    }
    __mSIMDd GaussianRandomNumberSequenceSIMD::GetRandomNumber(){
        if( m_Process + SIMDWIDTH > m_BankSize ){
            fprintf(stderr, "Error: random number bank is exhausted, increase bank size %10zd %10d\n",
                m_Process, m_BankSize);
            exit(0);
        }
        __mSIMDd _value = _mmSIMD_load_pd(&m_Bank[m_Process]);
        m_Process += SIMDWIDTH;
        return _value;
    }
    __mSIMDd GaussianRandomNumberSequenceSIMD::GetRandomNumberSafe(){
        if( m_Process + SIMDWIDTH > m_BankSize ){
            Renew();
        }
        __mSIMDd _value = _mmSIMD_load_pd(&m_Bank[m_Process]);
        m_Process += SIMDWIDTH;
        return _value;
    }
    size_t GaussianRandomNumberSequenceSIMD::GetProcess(){
        return m_Process;
    }
    int GaussianRandomNumberSequenceSIMD::GetCoreBankSize(){
        return m_BankSize;
    }
    void GaussianRandomNumberSequenceSIMD::PrintStatus(){
        fprintf(stderr, "Time:    %10lld\n", m_TimeOfRenew);
        fprintf(stderr, "Process: %10zd\n", m_Process);
        fprintf(stderr, "First    %13.7f\n", m_Bank[0]);
    }

    // =======================================================================
    UniformRandomNumberSequenceSIMD::UniformRandomNumberSequenceSIMD(){
    };
    UniformRandomNumberSequenceSIMD::~UniformRandomNumberSequenceSIMD(){
        vslDeleteStream(&m_Stream);
        m_Bank.clear();
        m_Bank.shrink_to_fit();
    };
    void UniformRandomNumberSequenceSIMD::Initialize(size_t bankSize, int seed, double lowBound, double highBound){
        m_LowBound = lowBound;
        m_HighBound = highBound;
        m_TimeOfRenew = 0;
        m_BankSize = bankSize;
        m_Bank.resize(m_BankSize);
        // this time random number is taken care of in the runnner
        if( seed == -1 ){
            seed = std::random_device{}();
        }
        vslNewStream(&m_Stream, VSL_BRNG_MT19937, seed);
        m_Process = m_Bank.size();
        Renew();
        // return seed;
    }
    void UniformRandomNumberSequenceSIMD::Renew(){
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, m_Stream, m_BankSize, &m_Bank[0], m_LowBound, m_HighBound);
        m_Process = 0;
        m_TimeOfRenew++;
    }
    void UniformRandomNumberSequenceSIMD::ResetRenewTime(){
        m_TimeOfRenew = 0;
    }
    long long UniformRandomNumberSequenceSIMD::GetRenewTime(){
        return m_TimeOfRenew;
    }
    void UniformRandomNumberSequenceSIMD::SetRenewTime(long long value){
        m_TimeOfRenew = value;
    }
    void UniformRandomNumberSequenceSIMD::SetProcess(size_t value){
        m_Process = value;
    }
    void UniformRandomNumberSequenceSIMD::Advance(long long count, size_t progress){
        for(long long i=0;i<count-1;i++){
            vslSkipAheadStream(m_Stream, (long long)(m_BankSize));
        }
        vslSkipAheadStream(m_Stream, (long long)(progress));
    }
    void UniformRandomNumberSequenceSIMD::AdvanceFromZero(long long count, size_t progress){
        if( count == 0 ) count = 1;
        if( count == 1 ){
            m_Process = progress;
        }else{
            // when doing nothing, count is 1 (renew'ed in Initialization)
            // so subtract 2, one is the initialization renew, the other is the following renew
            for(long long i=0;i<count-2;i++){
                vslSkipAheadStream(m_Stream, (long long)(m_BankSize));
            }
            // subtract 1, which is the following renew
            m_TimeOfRenew = count - 1;
            Renew();
            m_Process = progress;
        }
    }
    __mSIMDd UniformRandomNumberSequenceSIMD::GetRandomNumber(){
        if( m_Process + SIMDWIDTH > m_BankSize ){
            fprintf(stderr, "Error: random number bank is exhausted, increase bank size %10zd %10d\n",
                m_Process, m_BankSize);
            exit(0);
        }
        __mSIMDd _value = _mmSIMD_load_pd(&m_Bank[m_Process]);
        m_Process += SIMDWIDTH;
        return _value;
    }
    __mSIMDd UniformRandomNumberSequenceSIMD::GetRandomNumberSafe(){
        if( m_Process + SIMDWIDTH > m_BankSize ){
            Renew();
        }
        __mSIMDd _value = _mmSIMD_load_pd(&m_Bank[m_Process]);
        m_Process += SIMDWIDTH;
        return _value;
    }
    size_t UniformRandomNumberSequenceSIMD::GetProcess(){
        return m_Process;
    }
    int UniformRandomNumberSequenceSIMD::GetCoreBankSize(){
        return m_BankSize;
    }
    void UniformRandomNumberSequenceSIMD::PrintStatus(){
        fprintf(stderr, "Time:    %10lld\n", m_TimeOfRenew);
        fprintf(stderr, "Process: %10zd\n", m_Process);
        fprintf(stderr, "First    %13.7f\n", m_Bank[0]);
    }

};

